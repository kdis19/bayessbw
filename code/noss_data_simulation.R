source("code/functions.R")

all.days.df <- read.csv('data/all_days_df.csv')
nList <- lme4:::namedList
stages <- paste0('L', 2:6)

##### Prior Sampling #####
prior.samp <- function(chains) {
  l <- lapply(1:chains, function(x) {
    tl <- exp(rnorm(1, 5.64, 0.0067))
    tdiff <- rgamma(1, 112, scale = 0.226) 
    th <- tl + tdiff
    s_upsilon <- exp(rnorm(5, -2.5, 0.05))
    sm1 <- 0:4
    k <- -1.05*sm1^2 + 4.22*sm1 + 4.08
    lst <- list(
      'rho' = rgamma(5, k, scale = 0.045),
      'HA' = rgamma(1, 5.4, scale = 0.134),
      'TL' = tl,
      'HL' = -rgamma(1, 3.6, scale = 2.253),
      'TH' = th,
      'HH' = rgamma(1, 7.6, scale = 3.12),
      's_eps' = exp(rnorm(5, -1.5, 0.1)),
      's_upsilon' = s_upsilon,
      'upsilon' = exp(rnorm(35, 0, rep(s_upsilon, 7)))
    )
    return(lst)
  })
  return(l)
}

##### Data Generation #####
on.add <- subset(all.days.df, province == 'ON' & transfer)
on.add$time1 <- round(on.add$time1)
treat.times <- aggregate(data = on.add, time1 ~ temp1 + stage, 
                         function(x) {which.max(tabulate(x))})

## Individual progress for a given set of parameters
ind.progress <- function(param.df, dev.df) {
  sustainable <- ifelse(dev.df$treatment[1] == dev.df$sustainable[1], TRUE, FALSE)
  days <- vector()
  eps <- rnorm(5, 0, param.df$s_eps)
  deltas <- exp(eps)
  names(deltas) <- stages
  dev <- 0
  i <- 1
  dev.lst <- list()
  extra.dev <- 0
  dev.df$r.treatment <- dev.df$r.treatment/deltas
  dev.df$r.sustainable <- dev.df$r.sustainable/deltas
  
  for (s in 1:length(stages)) {
    st <- stages[s]
    rdf <- subset(dev.df, stage == st)
    trt.time <- (1/rdf$r.treatment)*(1-extra.dev)
    trt.time <- pmax(trt.time, 0)
    
    if (sustainable) {
      obs.time <- ceiling(trt.time)
      rdf$time.treatment <- 0
      rdf$time.sustainable <- obs.time
      if (st != 'L6') {
        extra.time <- obs.time - trt.time
        rdf.new <- subset(dev.df, stage == stages[s+1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
      dev.lst <- append(dev.lst, list(rdf))
      next
    }
    
    if (trt.time > rdf$time.treatment) {
      rem.prop <- 1 - extra.dev - rdf$r.treatment*rdf$time.treatment
      sus.time <- rem.prop/rdf$r.sustainable
      obs.time <- ceiling(sus.time)
      rdf$time.sustainable <- obs.time
      
      if (st != 'L6') {
        extra.time <- obs.time - sus.time
        rdf.new <- subset(dev.df, stage == stages[s + 1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
    }
    
    else {
      rdf$time.treatment <- ceiling(trt.time)
      rdf$time.sustainable <- 0
      if (st != 'L6') {
        if (extra.dev > 1) {extra.time <- extra.time - 1/rdf$r.sustainable}
        else {extra.time <- ceiling(trt.time) - trt.time}
        rdf.new <- subset(dev.df, stage == stages[s + 1])
        extra.dev <- extra.time*rdf.new$r.sustainable
      }
    }
    dev.lst <- append(dev.lst, list(rdf))
  }
  data <- bind_rows(dev.lst)
  data$nobs <- 1
  return(data)
  
}

## Calculate development rates from given set of parameters
pop.dev <- function(t.treat, param.df, size = 250) {
  if (t.treat %in% c(5, 10, 30, 35)) {all.temps <- c('treatment' = t.treat, 'sustainable' = 20)}
  else {all.temps <- c('treatment' = t.treat, 'sustainable' = t.treat)}
  dev.lst <- lapply(stages, function(st) {
    s <- which(stages == st)
    params <- subset(param.df, stage == st & temp == t.treat)
    params$ups_sus <- subset(param.df, stage == st & temp == all.temps['sustainable'])$upsilon
    times <- with(params, {
      TA <- ta.fun(HL, HH, TL, TH)
      cp1 <- calc_pred3(all.temps[1], rho, HA, TL, HL, TH, HH, TA)*upsilon
      cp2 <- calc_pred3(all.temps[2], rho, HA, TL, HL, TH, HH, TA)*ups_sus
      cvec <- c(unique(cp1), unique(cp2))
      names(cvec) <- names(all.temps)
      return(cvec)})
    rates <- 1/times
    names(rates) <- paste('r', names(all.temps))
    s.df <- as.data.frame(append(as.list(all.temps), as.list(rates)))
    s.df$stage <- st
    if (t.treat %in% c(5, 10, 30, 35))
    {s.df$time.treatment <- subset(treat.times, temp1 == all.temps['treatment'] & 
                                     stage == st)$time1}
    return(s.df)
  })
  dev.df <- bind_rows(dev.lst)
  l <- replicate(size, {
    ind.progress(param.df, dev.df)
  }, simplify = FALSE)
  df <- bind_rows(l)
  nms <- c('treatment', 'sustainable', 'stage', 'time.treatment', 'time.sustainable', 'nobs')
  df <- df[,nms]
  df$cup <- rep(1:size, each = 5)
  return(df)
}

## Wrapper function to generate populations
gen.pops <- function(N, ps = NULL) {
  if (is.null(ps)) {
    ps <- prior.samp(N)
  }
  else {N <- length(ps)}
  lst <- lapply(1:N, function(i) {
    p <- ps[[i]]
    pdf <- as.data.frame(p)
    pdf$stage <- rep(stages, 7)
    pdf$temp <- rep(seq(5, 35, by = 5), each = 5)
    t.lst <- lapply(seq(5, 35, by = 5), function(tmp) {
      pdat <- pop.dev(tmp, pdf)
    })
    pd <- bind_rows(t.lst)
    pd$prior.samp <- i
    return(list('priors' = p, 'data' = pd))
  })
  return(lst)
}
