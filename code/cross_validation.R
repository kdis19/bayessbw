library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)

all.days.df <- read.csv('data/all_days_df.csv')
all.days.df.ind <- read.csv('data/all_days_df_ind.csv')
all.days.df.ind <- subset(all.days.df.ind, select = -X)
source('code/functions.R')
source('code/noss_data_simulation.R')
nList <- lme4:::namedList

##### Separate datasets #####
K <- 10
all.days.df.ind$index <- paste0(all.days.df.ind$province, '_',
                                all.days.df.ind$temp1, '_',
                                all.days.df.ind$cup)
set.seed(543)
group.lst <- lapply(unique(all.days.df.ind$province), function(prov) {
  sub1 <- subset(all.days.df.ind, province == prov)
  n.ind <- 1
  p.lst <- list()
  for (tmp in unique(sub1$temp1)) {
    sub2 <- subset(sub1, temp1 == tmp)
    u.cups <- sort(unique(sub2$cup))
    n.cups <- length(u.cups)
    m.nper <- floor(n.cups/K)
    rmd <- n.cups %% K
    gp.vec <- rep(1:K, m.nper)
    if (rmd != 0) {
      for (i in 1:rmd) {
        gp.vec <- c(gp.vec, n.ind)
        n.ind <- ifelse(n.ind == 10, 1, n.ind + 1)
      }
    }
    
    gp <- sample(gp.vec, size = n.cups, replace = FALSE)
    gp.key <- data.frame('cup' = u.cups, 'group' = gp)
    
    sub2 <- merge(sub2, gp.key)
    p.lst <- append(p.lst, list(sub2))
  }
  sub.new <- bind_rows(p.lst)
  return(sub.new)
})
group.df <- bind_rows(group.lst)
group.df$run.index <- paste0(group.df$province, '_', group.df$group)

group.df$nobs <- 1
group.df <- subset(group.df, select = -cup)
ag.group.df <- aggregate(data = group.df, nobs ~ ., sum)
write.csv(ag.group.df, 'data/all_data_grouped.csv', row.names = FALSE)

basename <- "noss_regniere_structured_reduced"
basename_loc <- paste0('code/', basename)
cc <- compile(paste0(basename_loc, '.cpp'))
try(dyn.unload(dynlib(basename_loc)),silent=TRUE)
dyn.load(dynlib(basename_loc))
options(mc.cores = parallel::detectCores())


parms <- list(
  'rho' = 0.045*(-1.05*(0:4)^2 + 4.22*(0:4) + 4.08),
  'HL' = 8.4,
  'HH' = 9.3,
  'HA' = 1.4, 
  'TL' = 284,
  'TH' = 304,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5))
parms1 <- parms
diagnostics <- list()
chains <- 4
parm.lst <- prior.samp(chains)

set.seed(123)
gis <- unique(ag.group.df$run.index)
stages <- unique(all.days.df$stage)

for (gi in gis) {
  prov <- strsplit(gi, split = '_')[[1]][1]
  prov.df <- subset(ag.group.df, province == prov)
  all.data <- subset(prov.df, run.index != gi)
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
  write.csv(as.data.frame(stan1), paste0('code/output/cv/', gi, '_post.csv'), row.names = FALSE)
  
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
  
  diag.df$run.index <- gi
  write.csv(diag.df, paste0('code/output/cv/', gi, '_diagnostics.csv'), row.names = FALSE)
}

diag.lst <- lapply(gis, function(x) {
  read.csv(paste0('code/output/cv/', x, '_diagnostics.csv'))
})
write.csv(bind_rows(diag.lst), 'code/output/cv/diagnostics.csv', row.names = FALSE)
