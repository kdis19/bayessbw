library(tidyverse)
library(DPQ)
library(TMB)
library(reshape2)
source('code/functions.R')
source('code/noss_data_simulation.R')
all.days.df <- read.csv('data/all_days_df.csv')
nList <- lme4:::namedList
stages <- paste0('L', 2:6)
param.nms <- c('HA', 'TL', 'HL', 'TH', 'HH', 
               'rho', 's_eps', 's_upsilon', 'upsilon')

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
  'TL' = 285.90432,
  'TH' = 306.47530,
  's_eps' = rep(0.3, 5),
  's_upsilon' = rep(0.08, 5))


## run on subset of all.days.df
data.adjust <- function(data) {
  w1 <- which(data$temp1 %in% c(15, 20, 25) & data$time2 == 0)
  data$time2[w1] <- 1
  w2 <- which(data$temp1 %in% c(5, 10, 30, 35) & data$time1 == 0)
  data$time1[w2] <- 1
  
  data$time1.orig <- data$time1
  data$time2.orig <- data$time2
  
  data$l2 <- sapply(data$stage, function(x) {ifelse(x == 'L2', 1, 0)})
  data$cur0 <- sapply(data$time2.orig, function(x) {ifelse(x == 0, 1, 0)})
  data$prev0 <- c(0, data$cur0[1:(nrow(data) - 1)])
  
  data$time1 <- data$time1.orig + (1-data$l2)*data$prev0
  data$time2 <- data$time2.orig + (1-data$l2)*(1-data$prev0)
  data$time1d <- data$time1.orig - data$cur0
  data$time2d <- data$time2.orig - (1 - data$cur0)
  data <- subset(data, select = c(nobs, temp1, temp2, stage, time1, 
                                          time2, time1d, time2d))
  
  dd <- data
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
  
  return(as.data.frame(dd2))
}

## run on posteriors in code/output folder
param.extract <- function(post, list = FALSE) {
  samp.lst <- lapply(1:nrow(post), function(samp) {
    sub <- post[samp,]
    par.lst <- lapply(param.nms, function(x) {
      df <- sub[,grep(paste0('^', x), names(sub))]
      vec <- unlist(df)
      names(vec) <- rep(x, length(vec))
      return(vec)
    })
    if (list) {
      names(par.lst) <- param.nms
      return(par.lst)
    }
    else {
      upl <- unlist(par.lst)
      return(upl)
    }
  })
  if (!list) {
    samp.df <- do.call('rbind', samp.lst)
    return(samp.df)
  }
  else {
    return(samp.lst)
  }
}

ind.lk <- function(data, param) {
  llk <- with(param, {
    stage <- data$stage + 1
    t_block1 <- data$t_block1 + 1
    t_block2 <- data$t_block2 + 1
    
    TA <- ta.fun(-HL, HH, TL, TH)
    tpred1 <- calc_pred3(data$temp1, rho[stage], HA, TL, -HL, TH, HH, TA)
    tpred2 <- calc_pred3(data$temp2, rho[stage], HA, TL, -HL, TH, HH, TA)
    
    c_ups1 <- exp(upsilon[t_block1]*s_upsilon[stage])
    c_ups2 <- exp(upsilon[t_block2]*s_upsilon[stage])
    
    tpred1 <- tpred1*c_ups1
    tpred2 <- tpred2*c_ups2
    
    epsm1 <- log(data$time1d/tpred1 + data$time2d/tpred2)
    epsij <- log(data$time1/tpred1 + data$time2/tpred2)
    
    epsm1_std <- epsm1/s_eps[stage]
    epsij_std <- epsij/s_eps[stage]
    
    pnorm_m1 <- pnorm(epsm1_std, log.p = TRUE)
    pnorm_ij <- pnorm(epsij_std, log.p = TRUE)
    
    timesum <- data$time1d + data$time2d
    
    pnormdiff <- sapply(1:length(timesum), function(x) {
      pij <- pnorm_ij[x]
      pm1 <- pnorm_m1[x]
      ifelse(x == 0, pij,
             logspace.sub(pij, pm1))})
    
    lp.lst <- lapply(1:nrow(data), function(r) {
      rep(pnormdiff[r], data$nobs[r])
    })
    nlogp <- unlist(lp.lst)
    #nlogp <- data$nobs*pnormdiff
    return(nlogp)
  })
  return(llk)
}

## run on outputs of data.adjust and param.extract
likelihood.fn <- function(data, param, prior = FALSE, reduce = TRUE, expand = TRUE) {
  dd2 <- as.list(data)
  dd2$use_prior <- ifelse(prior, 1, 0)

  sm1 <- 0:4
  dd2$k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08
  
  nblock <- length(unique(dd2$stage))*length(unique(dd2$temp1))
  parms$upsilon <- rep(0, nblock)
  
  if (nblock == 30 & reduce) {
    g <- grep('^upsilon', colnames(param))
    w <- g[(0:4)*7 + 1]
    param <- param[,-w]
  }
  
  ff <- MakeADFun(data=dd2,
                  parameters=as.list(parms[param.nms]),
                  DLL=basename,
                  silent=TRUE)
  
  lk <- sapply(1:nrow(param), function(x) {
    px <- param[x,]
    g <- grep('^upsilon', names(px))
    if (expand & length(g) == 30) {
      ups <- px[g]
      ups.new <- vector()
      for (i in 1:length(ups)) {
        if ((i %% 6) == 1) {ups.new <- c(ups.new, 0)}
        ups.new <- c(ups.new, ups[i])
      }
      names(ups.new) <- rep('upsilon', length(ups.new))
      px <- px[-g]
      px <- c(px, ups.new)
    }
    ff$fn(px)
  })
  #lk <- lk/(sum(dd2$nobs)/5)
  return(-lk)
}

prior.eval <- function(param) {
  u.names <- unique(names(param))
  p.lst <- lapply(u.names, function(nm) {
    param[names(param) == nm]
  })
  names(p.lst) <- u.names
  sm1 <- 0:4
  k_vec <- -1.05*sm1^2 + 4.22*sm1 + 4.08
  pr <- with(p.lst, {
    jnll <- 0
    TA <- ta.fun(-HL, HH, TL, TH)
    r40 <- 1/calc_pred3(40, rho[2], HA, TL, -HL, TH, HH, TA)
    tdiff <- TH - TL
    jnll <- jnll - sum(dnorm(upsilon), log = TRUE)
    jnll <- jnll - dexp(r40, 40, log = TRUE)
    jnll <- jnll - dgamma(HL, 3.6, scale = 2.253, log = TRUE)
    jnll <- jnll - dgamma(HA, 5.4, scale = 0.134, log = TRUE)
    jnll <- jnll - dgamma(HH, 7.6, scale = 3.12, log = TRUE)
    jnll <- jnll - dnorm(log(TL), 5.64, 0.0067, log = TRUE)
    jnll <- jnll - dgamma(tdiff, 112, scale = 0.226, log = TRUE)
    jnll <- jnll - sum(dnorm(log(s_eps), -1.5, 0.1, log = TRUE))
    jnll <- jnll - sum(dnorm(log(s_upsilon), -2.5, 0.05, log = TRUE))
    for (i in 1:5) {
      jnll <- jnll - dgamma(rho[i], k_vec[i], scale = 0.045, log = TRUE)
    }
    return(jnll)
  })
  names(pr) <- NULL
  return(pr)
}

### Calculate Jensen-Shannon divergence 

## input two lists; eg. prov1 = list('data' = data1, 'post' = post1)
### post1 is in the form of code/output/[prov]_post.csv
### data1 is in the form of data.adjust(dat) where
###  dat is in the form of subset(all.days.df, province == [prov])
js.divergence <- function(prov1, prov2, N = 6000) {
  data <- list(prov1$data, prov2$data)
  
  post <- bind_rows(prov1$post, prov2$post)
  post.ps <- param.extract(post)
  
  p.lst <- lapply(1:2, function(p) {
    dat <- data[[p]]
    lk <- likelihood.fn(dat, post.ps)
    return(lk)
  })
  df <- as.data.frame(p.lst)
  names(df) <- c('prov1', 'prov2')
  
  df$mixture <- apply(df, 1, mean)
  df1 <- df[1:N,]
  df2 <- df[(N/2 + 1):N,]
  
  js <- 0.5*mean(df1[,'prov1'] - df1$mixture) + 0.5*mean(df2[,'prov2'] - df2$mixture)
  return(js)
}
