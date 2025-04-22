library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
library(reshape2)

all.days.df <- read.csv('data/all_days_df.csv') %>%
  subset(select = -X) 

source('code/functions.R')
source('code/prior_sampling.R')
nList <- lme4:::namedList

provs <- rev(sort(unique(all.days.df$province)))
stages <- unique(all.days.df$stage)

nstage <- length(stages)
ncolony <- length(provs)

##### Run model #####
basename <- "regniere_model"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

parms <- list(
  'rho' = 0.045*(-1.05*(0:4)^2 + 4.22*(0:4) + 4.08),
  'HL' = 8.1,
  'HH' = 24,
  'HA' = 0.72,
  'TL' = 284,
  'TH' = 305,
  'rho_mult' = rep(1, nstage*ncolony),
  'HL_mult' = rep(1, ncolony),
  'HH_mult' = rep(1, ncolony),
  'HA_mult' = rep(1, ncolony), 
  'TL_mult' = rep(1, ncolony),
  'TH_mult' = rep(1, ncolony),
  's_eps' = rep(0.3, nstage*ncolony),
  's_upsilon' = rep(0.08, 35))
parms1 <- parms
diagnostics <- list()
tries <- 30
chains <- 8

set.seed(125)
parm.lst.orig <- prior.samp.reduced(tries, nprov = 6, ab = TRUE)

parm.test <- sapply(parm.lst.orig, function(parm) {
  df <- as.data.frame(parm[1:6])
  gc35 <- get_curves(df, temp = 35)
  finite <- ifelse(any(log(gc35) < -5), FALSE, TRUE)
  return(finite)
})

parm.lst <- parm.lst.orig[parm.test][1:chains]
if (any(sapply(parm.lst, is.null))) {"Incomplete set of starting values!!!"}

all.days.df$province <- factor(all.days.df$province, levels = provs)

all.data <- all.days.df
all.data$colony <- as.numeric(factor(all.data$province)) - 1
all.data$block <- paste(all.data$stage, all.data$temp1, sep = '_')
lev.ord <- paste(rep(stages, each = 7), rep(seq(5, 35, by = 5), 5), sep = '_')
all.data$block <- factor(all.data$block)
all.data$block <- factor(all.data$block, levels = lev.ord)
w1 <- which(all.data$temp1 %in% c(15, 20, 25) & all.data$time2 == 0)
all.data$time2[w1] <- 1
w2 <- which(all.data$temp1 %in% c(5, 10, 30, 35) & all.data$time1 == 0)
all.data$time1[w2] <- 1

all.data$time1.orig <- all.data$time1
all.data$time2.orig <- all.data$time2

all.data$l2 <- sapply(all.data$stage, function(x) {ifelse(x == 'L2', 1, 0)})
all.data$cur0 <- sapply(all.data$time2.orig, function(x) {ifelse(x == 0, 1, 0)})
all.data$prev0 <- c(0, all.data$cur0[1:(nrow(all.data) - 1)])

all.data$time1 <- all.data$time1.orig + (1-all.data$l2)*all.data$prev0
all.data$time2 <- all.data$time2.orig + (1-all.data$l2)*(1-all.data$prev0)
all.data$time1d <- all.data$time1.orig - all.data$cur0
all.data$time2d <- all.data$time2.orig - (1 - all.data$cur0)
all.data <- subset(all.data, select = c(nobs, temp1, temp2, stage, time1, 
                                        time2, time1d, time2d, block, colony))

dd <- all.data
dd$stagename <- dd$stage
dd$stage <- as.numeric(factor(dd$stage)) - 1

ntreat <- length(unique(dd$temp1))

dd$t_block1 <- (dd$colony*ntreat + (40 - dd$temp1)/5 - 1)*nstage + dd$stage
dd$t_block2 <- (dd$colony*ntreat + (40 - dd$temp2)/5 - 1)*nstage + dd$stage

dd$eps_block <- dd$colony*nstage + dd$stage
dd$ups_block <- ((40 - dd$temp1)/5 - 1)*nstage + dd$stage

dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,colony,
                      t_block1, t_block2, eps_block, ups_block,
                      nobs=as.integer(nobs)))



# minval <- min(c(dd2$t_block1, dd2$t_block2)) 
# dd2$t_block1 <- dd2$t_block1 - minval
# dd2$t_block2 <- dd2$t_block2 - minval

dd2$use_prior <- 1
sm1 <- 0:4
dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08

#nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
nblock <- length(unique(dd2$t_block1))

parms$upsilon <- rep(0, nblock)

pnames <- c("rho",
            "HA","TL","HL","TH","HH",
            'rho_mult',
            "HA_mult","TL_mult","HL_mult","TH_mult","HH_mult",
            "s_eps", 's_upsilon', 
            'upsilon')

parm.lst <- lapply(parm.lst, function(x) {
  x$HL <- abs(x$HL)
  #x$upsilon <- rep(0, nstage*ntreat)
  return(x)
})

test.vals <- sapply(parm.lst, function(x) {
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(x[pnames]),
                  DLL=basename,
                  random=c('upsilon'),
                  silent=TRUE)
  ff$fn(ff$par)
})

if (all(is.finite(test.vals))) {print('starting values all good')}

ff <- MakeADFun(data=dd2,
                parameters=as.list(parms[pnames]),
                DLL=basename,
                random=c('upsilon'),
                silent=TRUE)

stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE,
                 chains = length(parm.lst), iter = 1375, warmup = 700,
                 control = list(max_treedepth = 15, adapt_delta = 0.99))
saveRDS(stan1, 'code/output/stanfit_sups_st.rds')
# stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE,
#                  chains = 7, iter = 750, warmup = 500)
post.df <- as.data.frame(stan1)
# write.csv(post.df, 'code/output/all_provs_post.csv', row.names = FALSE)

stan.array <- as.array(stan1)
Rhat <- max(apply(stan.array, 3, Rhat))
ESS.bulk <- min(apply(stan.array, 3, ess_bulk))
ESS.tail <- min(apply(stan.array, 3, ess_tail))

g <- get_sampler_params(stan1, inc_warmup = FALSE)
g.df <- do.call('rbind', g)
div <- sum(g.df[,'divergent__'])
diag.df <- data.frame('Rhat' = Rhat, 
                      'ESS.bulk' = ESS.bulk,
                      'ESS.tail' = ESS.tail,
                      'divergences' = div)
write.csv(diag.df, 'code/output/all_provs_diag.csv', row.names = FALSE)

