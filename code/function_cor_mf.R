prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    #tl <- exp(rnorm(1, 5.64, 0.0067))
    tl <- exp(rnorm(1, 5.67, 0.01))
    tdiff <- rgamma(1, 112, scale = 0.226) 
    th <- tl + tdiff
    s_upsilon <- exp(rnorm(5, -2.5, 0.05))
    sm1 <- 0:4
    k <- -1.05*sm1^2 + 4.22*sm1 + 4.08
    lst <- list(
      'rho' = rgamma(5, k, scale = 0.045),
      'HA' = rgamma(1, 5.4, scale = 0.134),
      'TL' = rnorm(1, 284, 2),
      #'TL' = tl,
      'HL' = -rgamma(1, 3.6, scale = 2.253),
      'TH' = rnorm(1, 304, 2),
      #'TH' = th,
      'HH' = rgamma(1, 7.6, scale = 3.12),
      's_eps' = exp(rnorm(5, -1.5, 0.1)),
      's_upsilon' = s_upsilon,
      #'upsilon' = rlnorm(35, 0, rep(s_upsilon, each = 7)),
      'upsilon' = exp(rnorm(35, 0, rep(s_upsilon, 7))),
      'male_l6' = rgamma(1, 3.25, scale = 0.05)
    )
    return(lst)
  })
  # wvec <- sapply(l, function(x) {
  #   row <- as.data.frame(x)[3,]
  #   row$rho25 <- row$rho
  #   gc <- get_curves(row, temp = c(0, 40))
  #   probs <- dexp(gc$rate, rate = 20)
  # })
  # wt2 <- apply(wvec, 1, function(x) {x/sum(x)})
  # wt <- apply(wt2, 1, prod)
  # wt.st <- wt/sum(wt)
  # samp <- sample(1:length(wt.st), size = chains, replace = TRUE, prob = wt.st)
  return(l)
}

data.adjust <- function(data, cup = FALSE) {
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
  if (!cup) {
    data <- subset(data, select = c(nobs, temp1, temp2, stage, time1, 
                                    time2, time1d, time2d))
    dd <- data
    dd$stagename <- dd$stage
    dd$stage <- as.numeric(factor(dd$stage)) - 1
    dd2 <- with(dd, nList(temp1,time1,temp2,time2,time1d,time2d,stage,
                          nobs=as.integer(nobs)))
  }
  else {
    data <- subset(data, select = c(nobs, temp1, temp2, stage, time1, 
                                    time2, time1d, time2d, index))
    dd <- data
    dd$stagename <- dd$stage
    dd$stage <- as.numeric(factor(dd$stage)) - 1
    dd2 <- with(dd, nList(temp1,time1,temp2,time2,
                          time1d,time2d,stage,index,
                          nobs=as.integer(nobs)))
  }
  ntreat <- length(unique(dd2$temp1))
  dd2$t_block1 <- dd2$stage*ntreat + dd2$temp1/5 - 1
  dd2$t_block2 <- dd2$stage*ntreat + dd2$temp2/5 - 1
  
  minval <- min(c(dd2$t_block1, dd2$t_block2)) 
  dd2$t_block1 <- dd2$t_block1 - minval
  dd2$t_block2 <- dd2$t_block2 - minval
  
  return(as.data.frame(dd2))
}