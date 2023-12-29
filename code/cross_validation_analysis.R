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
prov.lst <- lapply(provs, function(p) {
  print(p)
  gps <- sapply(1:10, function(x) {
    print(x)
    data <- subset(group.df, run.index == paste0(p, '_', x))
    post <- read.csv(paste0('code/output/cv/', p, '_', x, '_post.csv'))
    
    adj.data <- data.adjust(data)
    pe.post <- param.extract(post, list = TRUE)
    
    lk.vec <- sapply(pe.post, function(pe) {
      lk <- exp(ind.lk(adj.data, pe))
      m <- mean(lk)
      return(log(m))
    })
    elpd.i <- log(mean(exp(lk.vec)))
    
    # lk <- likelihood.fn(data = adj.data, param = pe.post, 
    #                     reduce = FALSE, expand = FALSE)
    # elpd.i <- log(mean(exp(lk)))
    return(elpd.i)
  })
  return(gps)
})
prov.df <- data.frame('post.prov' = rep(provs, each = 10),
                      'group' = rep(1:10, length(provs)),
                      'elpd' = unlist(prov.lst),
                      'prov' = rep(provs, each = 10))
prov.df.ag <- aggregate(data = prov.df, elpd ~ prov + post.prov, mean)
write.csv(prov.df.ag, 'code/output/cv/prov_individual.csv', row.names = FALSE)

prov.comp.lst <- lapply(provs, function(p) {
  print(p)
  dat <- subset(all.days.df, province == p)
  wp <- grep(p, post.files)
  pf.vec <- post.files[-wp]
  post.lst <- lapply(pf.vec, function(y) {
    st.vec <- strsplit(y, split = '_')[[1]][1:2]
    print(paste(st.vec, collapse = '_'))
    post <- read.csv(paste0('code/output/cv/', y))
    adj.data <- data.adjust(dat)
    pe.post <- param.extract(post, list = TRUE)
    lk.vec <- sapply(pe.post, function(pe) {
      if (st.vec[1] == 'AB') {
        ups <- pe$upsilon
        ups.new <- vector()
        for (i in 1:length(ups)) {
          if ((i %% 6) == 1) {ups.new <- c(ups.new, 0)}
          ups.new <- c(ups.new, ups[i])
        }
        names(ups.new) <- rep('upsilon', length(ups.new))
        pe$upsilon <- ups.new
      }
      lk <- exp(ind.lk(adj.data, pe))
      m <- mean(lk)
      return(log(m))
    })
    elpd.i <- log(mean(exp(lk.vec)))
    df <- data.frame('post.prov' = st.vec[1],
                     'group' = st.vec[2],
                     'elpd' = elpd.i)
    return(df)
  })
  post.df <- bind_rows(post.lst)
  post.df$prov <- p
  return(post.df)
})
prov.comp.df <- bind_rows(prov.comp.lst)
prov.comp.df.all <- rbind(prov.comp.df, prov.df)

write.csv(prov.comp.df.all, 'code/output/cv/prov_compare.csv', row.names = FALSE)

pcdf.ag <- aggregate(data = prov.comp.df, elpd ~ prov + post.prov, mean)
pcdf.all.ag <- rbind(pcdf.ag, prov.df.ag)

lks <- lapply(1:length(post.files), function(x) {
  ind <- names(post.files)[x]
  data <- subset(group.df, run.index == ind)
  files <- post.files[-x]
  post.lst <- lapply(files, function(y) {
    post <- read.csv(paste0('code/output/cv/', y))
    
  })
})