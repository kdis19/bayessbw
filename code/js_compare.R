source('code/noss_data_simulation.R')
source('code/likelihood.R')
post.clean <- read.csv('code/output/post_cleaned.csv')
all.days.df <- read.csv('data/all_days_df.csv')
provs <- unique(all.days.df$province)

### Test with simulated data ###
js.sims.lst <- list()
for (p in provs) {
  print(p)
  data <- subset(all.days.df, province == p)
  adj <- data.adjust(data)
  post.orig <- read.csv(paste0('code/output/forcerate/', p, '_post.csv'))
  post.c <- subset(post.clean, province == p)
  samp <- sample(unique(post.c$index), 25)
  post.sub <- subset(post.c, index %in% samp)
  
  gp.lst <- gen.pops(ps = post.sub)
  sims <- lapply(gp.lst, function(x) {
    dat <- x$data
    if (p == 'AB') {
      dat <- subset(dat, treatment != 5)
    }
    ind <- unique(dat$prior.samp)
    newdat <- adjust.gp(dat)
    return(newdat)
  })
  
  js.sims <- sapply(1:length(sims), function(i) {
    prov1 <- list('data' = adj, 'post' = post.orig)
    prov2 <- list('data' = sims[[i]], 'post' = post.orig)
    js <- js.divergence(prov1, prov2)
    return(js)
  })
  js.sims.lst <- append(js.sims.lst, list(js.sims))
}

prov.lst <- lapply(provs, function(p) {
  print(p)
  data <- subset(all.days.df, province == p)
  adj <- data.adjust(data)
  post.orig <- read.csv(paste0('code/output/forcerate/', p, '_post.csv'))
  post.c <- subset(post.clean, province == p)
  samp <- sample(unique(post.c$index), 25)
  post.sub <- subset(post.c, index %in% samp)
  
  gp.lst <- gen.pops(ps = post.sub)
  sims <- lapply(gp.lst, function(x) {
    dat <- x$data
    if (p == 'AB') {
      dat <- subset(dat, treatment != 5)
    }
    ind <- unique(dat$prior.samp)
    newdat <- adjust.gp(dat)
    return(newdat)
  })
  
  js.sims <- sapply(1:length(sims), function(i) {
    prov1 <- list('data' = adj, 'post' = post.orig)
    prov2 <- list('data' = sims[[i]], 'post' = post.orig)
    js <- js.divergence(prov1, prov2)
    return(js)
  })
  return(js.sims)
})

js.other <- sapply(provs[!(provs %in% c('ON', 'AB'))], function(p) {
  print(p)
  prov1 <- list('data' = on.adj, 'post' = on.post.orig)
  dat2.orig <- subset(all.days.df, province == p)
  dat2 <- data.adjust(dat2.orig)
  post2 <- read.csv(paste0('code/output/forcerate/', p, '_post.csv'))
  prov2 <- list('data' = dat2, 'post' = post2)
  js <- js.divergence(prov1, prov2)
  return(js)
})
names(js.other) <- provs[!(provs %in% c('ON', 'AB'))]
