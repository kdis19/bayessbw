library(tidyverse)
source('code/functions.R')
all.days.df <- read.csv('data/all_days_df.csv')
provs <- unique(all.days.df$province)
stages <- unique(all.days.df$stage)
nstage <- length(stages)
pnames <- c("rho",
            "HA","TL","HL","TH","HH",
            "s_eps", 's_upsilon', 
            'upsilon')

##### Organize results #####
process.standf <- function(province) {
  data <- read.csv(paste0('code/output/', province, '_post.csv'))
  data$HL <- -abs(data$HL)
  
  d.post.lst <- lapply(pnames, function(nm) {
    sub <- data[,grep(paste0('^', nm), names(data))]
    if (length(sub) == nstage) {
      new <- as.data.frame(unlist(lapply(1:nrow(sub), function(x) {unlist(sub[x,])})))
      names(new) <- nm
    }
    else if (length(sub) <= 35) {
      ntreat <- length(sub)/nstage
      lst <- lapply(1:ntreat, function(i) {
        inds <- ntreat*(0:4) + i
        ups <- sub[,inds]
        df <- as.data.frame(unlist(lapply(1:nrow(ups), function(x) {unlist(ups[x,])})))
        names(df) <- paste0('upsilon.', i + 7 - ntreat, '.')
        return(df)
      })
      new <- bind_cols(lst)
      if (ntreat < 7) {
        new <- bind_cols(data.frame('upsilon.1.' = rep(0, nrow(new))), new)
      }
    }
    else {
      new <- as.data.frame(rep(unlist(sub), each = nstage))
      names(new) <- nm
    }
    return(as.data.frame(new))
  })
  d.post.df <- bind_cols(d.post.lst)
  
  d.post.df$stage <- rep(stages, nrow(data))
  d.post.df$index <- rep(1:nrow(data), each = nstage)
  d.post.df$province <- province
  return(d.post.df)
}

pprovs.lst <- lapply(provs, function(p) {process.standf(p)})
pprovs.df <- bind_rows(pprovs.lst)
write.csv(pprovs.df, 'code/output/post_cleaned.csv', row.names = FALSE)

##### Calculate Dev Curves #####
pprovs.df$index <- paste0(pprovs.df$province, '_', pprovs.df$stage)

gc.lst <- lapply(unique(pprovs.df$index), function(ind) {
  sub <- subset(pprovs.df, index == ind)
  gc <- get_curves(sub, spline = TRUE)
  ind.vec <- strsplit(ind, split = '_')[[1]]
  gc$stage <- ind.vec[2]
  gc$province <- ind.vec[1]
  return(gc)
})
gc.df <- bind_rows(gc.lst)
gc.df$dev <- gc.df$rate/24

ag.gc.orig <- aggregate(data = gc.df, rate ~ temp + province + stage,
                        function(x) {
                          q <- quantile(x, probs = c(0.05, 0.5, 0.95))
                          return(q)})
ag.gc1 <- subset(ag.gc.orig, select = c(temp, province, stage))
ag.gc2 <- ag.gc.orig$rate
ag.gc2 <- as.data.frame(ag.gc2)
names(ag.gc2) <- c('low', 'mid', 'high')
ag.gc <- bind_cols(ag.gc1, ag.gc2)

ag.gc$province[ag.gc$province == 'NB2'] <- 'NB'

write.csv(ag.gc, 'code/output/ag_devcurves.csv', row.names = FALSE)

total.samples <- nrow(pprovs.df)/(length(stages)*length(provs))
samp <- sample(1:total.samples, 500)

gc.df.sub <- subset(gc.df, index %in% samp)
write.csv(gc.df.sub, 'code/output/all_posterior_curves.csv', row.names = FALSE)
