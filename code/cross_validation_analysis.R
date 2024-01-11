source('code/likelihood.R')

group.df <- read.csv('data/all_data_grouped.csv')
provs <- unique(group.df$province)

post.files <- list.files('code/output/cv', pattern = 'post')
nms <- sapply(post.files, function(pf) {
  ss <- strsplit(pf, split = '_')[[1]][1:2]
  nm <- paste0(ss[1], '_', ss[2])
  return(nm)
})
names(post.files) <- nms

diag.files <- list.files('code/output/cv', pattern = 'diagnostics')
diag.lst <- lapply(diag.files, function(x) {
  vec <- strsplit(x, split = '_')[[1]][1:2]
  dat <- read.csv(paste0('code/output/cv/', x))
  dat$province <- vec[1]
  dat$group <- vec[2]
  return(dat)
})
diag.df <- bind_rows(diag.lst)
write.csv(diag.df, 'code/output/cv/diagnostics.csv', row.names = FALSE)

##### Evaluate elpd within colonies #####
prov.lst <- list()
for (p in provs) {
  print(p)
  gps.lst <- lapply(1:10, function(x) {
    print(x)
    data <- subset(group.df, run.index == paste0(p, '_', x))
    post <- read.csv(paste0('code/output/cv/', p, '_', x, '_post.csv'))
    
    adj.data <- data.adjust(data)
    pe.post <- param.extract(post, list = TRUE)
    
    llk.lst <- lapply(pe.post, function(pe) {ind.lk(adj.data, pe)})
    llk.mat <- do.call('rbind', llk.lst)
    lk.mat <- exp(llk.mat)
    
    ind.vec <- apply(lk.mat, 2, function(x) {log(mean(x))})
    
    return(ind.vec)
  })
  gps <- unlist(gps.lst)
  prov.lst <- append(prov.lst, list(gps))
}
names(prov.lst) <- provs
saveRDS(prov.lst, file = 'code/output/cv/prov_individual.rds')

col.post.lst <- lapply(provs, function(x) {
  read.csv(paste0('code/output/', x, '_post.csv'))
})
names(col.post.lst) <- provs

prov.comp.lst <- list()
pcl.names <- vector()
comp.lst <- list()
for (p in provs) {
  dat <- subset(group.df, province == p)
  provs2 <- provs[provs != p]
  for (p2 in provs2) {
    nm <- paste(c(p, p2), collapse = '_')
    pcl.names <- c(pcl.names, nm)
    print(nm)
    post <- col.post.lst[[p2]]
    pe.post <- param.extract(post, list = TRUE)
    
    adj.data <- data.adjust(dat)
    llk.lst <- lapply(pe.post, function(pe) {
      if (p2 == 'AB') {
        ups <- pe$upsilon
        ups.new <- vector()
        for (i in 1:length(ups)) {
          if ((i %% 6) == 1) {ups.new <- c(ups.new, 0)}
          ups.new <- c(ups.new, ups[i])
        }
        names(ups.new) <- rep('upsilon', length(ups.new))
        pe$upsilon <- ups.new
      }
      lk <- ind.lk(adj.data, pe)
      return(lk)
    })
    llk.mat <- do.call('rbind', llk.lst)
    lk.mat <- exp(llk.mat)
    ind.vec <- apply(lk.mat, 2, function(x) {log(mean(x))})
    
    prov.comp.lst <- append(prov.comp.lst, list(ind.vec))
    
    iv.orig <- post.lst[[p]]
    m.comp <- sum(iv.orig) - sum(ind.vec)
    se.comp <- sqrt(length(ind.vec)*var(iv.orig - ind.vec))
    df <- data.frame('prov.orig' = p,
                     'prov.comp' = p2,
                     'diff' = m.comp,
                     'se.diff' = se.comp)
    comp.lst <- append(comp.lst, list(df))
  }
}
names(prov.comp.lst) <- pcl.names
comp.df <- bind_rows(comp.lst)
saveRDS(prov.comp.lst, file = 'code/output/cv/prov_compare.rds')
#write.csv(comp.df, 'code/output/cv/prov_compare.csv', row.names = FALSE)

##### Summarize #####
s.prov.lst <- lapply(provs, function(p) {
  x <- prov.lst[[p]]
  mu <- mean(x)
  sigma <- sd(x)
  se <- sigma/sqrt(length(x))
  return(data.frame('prov.orig' = p,
                    'prov.comp' = p,
                    'mean' = mu,
                    'se' = se))
})
s.prov.df <- bind_rows(s.prov.lst)

s.prov.comp.lst <- lapply(names(prov.comp.lst), function(nm) {
  x <- prov.comp.lst[[nm]]
  pvec <- strsplit(nm, split = '_')[[1]]
  mu <- mean(x)
  sigma <- sd(x)
  se <- sigma/sqrt(length(x))
  return(data.frame('prov.orig' = pvec[1],
                    'prov.comp' = pvec[2],
                    'mean' = mu,
                    'se' = se))
})
s.prov.comp.df <- bind_rows(s.prov.comp.lst)

s.all.df <- bind_rows(s.prov.df, s.prov.comp.df)
write.csv(s.all.df, 'code/output/cv/s_col_compare.csv', row.names = FALSE)

ad.lst <- lapply(provs, function(p) {
  dat <- subset(group.df, province == p)
  adj.data <- data.adjust(dat, cup = TRUE)
  adj.data$province <- p
  return(adj.data)
})
ad.df <- bind_rows(ad.lst)

pl <- lapply(provs, function(p) {
  #print(p)
  x <- prov.lst[[p]]
  g.lst <- lapply(1:10, function(g) {
    #print(g)
    sub <- subset(group.df, province == p & group == g)
    da <- data.adjust(sub, cup = TRUE)
    #print(nrow(sub) == nrow(da))
    df <- data.frame('prov.orig' = p, 
                     'prov.comp' = p,
                     'stage' = stages[da$stage + 1], 
                     'temp' = da$temp1,
                     'index' = da$index)
    return(df)
  })
  g.df <- bind_rows(g.lst)
  g.df$elpd <- x
  return(g.df)
})
pl.df <- bind_rows(pl)

pcl <- lapply(names(prov.comp.lst), function(nm) {
  #print(nm)
  x <- prov.comp.lst[[nm]]
  pvec <- strsplit(nm, split = '_')[[1]]
  p1 <- pvec[1]
  p2 <- pvec[2]
  sub <- subset(ad.df, province == p1)
  df <- data.frame('prov.orig' = p1, 
                   'prov.comp' = p2,
                   'elpd' = x,
                   'stage' = stages[sub$stage + 1], 
                   'temp' = sub$temp1,
                   'index' = sub$index)
  return(df)
})
pcl.df <- bind_rows(pcl)
pcl.df <- pcl.df[,names(pl.df)]

p.df <- bind_rows(pl.df, pcl.df)
p.df$cup <- sapply(p.df$index, function(x) {
  strsplit(x, split = '_')[[1]][3]})

ag.pdf <- aggregate(data = p.df, elpd ~ prov.orig + 
                      prov.comp + temp + cup, sum)
write.csv(ag.pdf, 'code/output/cv/full_elpd.csv', row.names = FALSE)

ag2.pdf <- aggregate(data = ag.pdf, elpd ~ prov.orig + 
                       prov.comp + temp, 
                     function(x) {
                         mu <- mean(x)
                         se <- sd(x)/sqrt(length(x))
                         return(c(mu, se))
                       })
elpd.df <- as.data.frame(ag2.pdf$elpd)
names(elpd.df) <- c('mean', 'se')
elpd.temp.df <- bind_cols(ag2.pdf[,1:3], elpd.df)
write.csv(elpd.temp.df, 'code/output/cv/elpd_temp.csv', row.names = FALSE)

##### Paired Comparison #####
paired.p.lst <- lapply(provs, function(p) {
  sub <- subset(ag.pdf, prov.orig == p)
  elpd.orig.df <- subset(sub, prov.comp == p)
  elpd.orig.df <- elpd.orig.df[order(elpd.orig.df$index),]
  elpd.orig <- elpd.orig.df$elpd
  provs2 <- provs[which(provs != p)]
  lst2 <- lapply(provs2, function(p2) {
    sub2 <- subset(sub, prov.comp == p2)
    sub2 <- sub2[order(sub2$index),]
    elpd.comp <- sub2$elpd
    paired.mean <- mean(elpd.orig - elpd.comp)
    paired.se <- sqrt(var(elpd.orig - elpd.comp)/length(elpd.orig))
    df <- data.frame('prov.data' = p,
                     'prov.post' = p2,
                     'mean' = paired.mean,
                     'se' = paired.se)
    return(df)
  })
  df2 <- bind_rows(lst2)
  return(df2)
})
paired.p <- bind_rows(paired.p.lst)