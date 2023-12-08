library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)

all.days.df <- read.csv('data/all_days_df.csv')
source('code/functions.R')
source('code/noss_data_simulation.R')
source('code/likelihood.R')
nList <- lme4:::namedList

##### Run model #####
basename <- "noss_regniere_structured_reduced"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())

parms.init <- list(
  'rho' = 0.045*(-1.05*(0:4)^2 + 4.22*(0:4) + 4.08),
  'HL' = 8.4,
  'HH' = 9.3,
  'HA' = 1.4, 
  'TL' = 284,
  'TH' = 304,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5))
parms1 <- parms.init
diagnostics <- list()
chains <- 4
pl <- prior.samp(chains)

set.seed(123)
provs <- unique(all.days.df$province)
stages <- unique(all.days.df$stage)

model.run <- function(provs, parm.lst = pl, parms = parms.init) {
  all.data <- subset(all.days.df, province %in% provs)
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
                                          time2, time1d, time2d, block))
  all.data <- aggregate(data = all.data, nobs ~ ., sum)
  
  dd <- all.data
  dd$stagename <- dd$stage
  dd$stage <- as.numeric(factor(dd$stage)) - 1
  
  dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                        nobs=as.integer(nobs)))
  ntreat <- length(unique(dd2$temp1))
  dd2$t_block1 <- dd2$stage*ntreat + dd2$temp1/5 - 1
  dd2$t_block2 <- dd2$stage*ntreat + dd2$temp2/5 - 1
  
  minval <- min(c(dd2$t_block1, dd2$t_block2)) 
  dd2$t_block1 <- dd2$t_block1 - minval
  dd2$t_block2 <- dd2$t_block2 - minval
  
  dd2$use_prior <- 1
  sm1 <- 0:4
  dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08
  
  nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
  nstage <- length(unique(dd2$stage))
  
  parms$upsilon <- rep(0, nblock)
  
  pnames <- c("rho",
              "HA","TL","HL","TH","HH",
              "s_eps", 's_upsilon', 
              'upsilon')
  
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(parms[pnames]),
                  DLL=basename,
                  random=c('upsilon'),
                  silent=TRUE)
  
  parm.lst <- lapply(parm.lst, function(x) {
    x$HL <- abs(x$HL)
    x$upsilon <- rep(0, nstage*ntreat)
    return(x)
  })
  stan1 <- tmbstan(ff, init = parm.lst, silent = TRUE, 
                   chains = chains, iter = 1500, warmup = 750,
                   control = list(max_treedepth = 15, adapt_delta = 0.99))
  p.write <- paste(provs, collapse = '_')
  write.csv(as.data.frame(stan1), paste0('code/output/', p.write, '_post.csv'), row.names = FALSE)
  
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
  
  write.csv(diag.df, paste0('code/output/', p.write, '_diagnostics.csv'), row.names = FALSE)
  
  return(stan1)
}

combo.vec <- vector()
stan.lst <- list()
p.remove <- vector()
for (p in provs) {
  p.remove <- c(p.remove, p)
  pr.vec <- provs[!(provs %in% p.remove)]
  if (length(pr.vec) == 0) {break}
  for (p2 in pr.vec) {
    ### UNCOMMENT TO RUN ALL PAIRED MODELS ###
    # if (p %in% c('ON', 'IPU') & p2 %in% c('ON', 'IPU')) {next}
    # mod.run <- model.run(provs = c(p, p2))
    # stan.lst <- append(stan.lst, list(mod.run))
    # names(stan.lst)[length(stan.lst)] <- paste0(p, '_', p2)
    combo.vec <- c(combo.vec, paste0(p, '_', p2))
  }
}

##### Analyze results #####
post.lst <- lapply(provs, function(p) {
  read.csv(paste0('code/output/', p, '_post.csv'))
})
names(post.lst) <- provs

lk.lst <- lapply(combo.vec, function(pvs) {
  print(pvs)
  prov.vec <- strsplit(pvs, split = '_')[[1]]
  p1 <- prov.vec[1]
  p2 <- prov.vec[2]
  m.data <- data.adjust(subset(all.days.df, province %in% prov.vec))
  m.post <- param.extract(read.csv(paste0('code/output/', pvs, '_post.csv')))
  
  data1 <- data.adjust(subset(all.days.df, province == p1))
  data2 <- data.adjust(subset(all.days.df, province == p2))
  post1 <- param.extract(post.lst[[p1]])
  post2 <- param.extract(post.lst[[p2]])
  post1 <- post1[1:pmin(nrow(post1), nrow(post2)),]
  post2 <- post2[1:pmin(nrow(post1), nrow(post2)),]
  
  m1.1.lk <- likelihood.fn(data1, post1)
  m1.2.lk <- likelihood.fn(data2, post2)
  m1.1.prior <- sapply(1:nrow(post1), function(x) {prior.eval(post1[x,])})
  m1.2.prior <- sapply(1:nrow(post2), function(x) {prior.eval(post2[x,])})
  m1.post <- m1.1.lk + m1.2.lk + m1.1.prior + m1.2.prior
  
  m2.lk <- likelihood.fn(m.data, m.post)
  m2.prior <- sapply(1:nrow(m.post), function(x) {prior.eval(m.post[x,])})
  m2.post <- m2.lk + m2.prior
  
  df <- data.frame('pvs' = pvs,
                   'l.bayesFactor' = mean(m1.post) - mean(m2.post))
 
})
lk.df <- bind_rows(lk.lst)
