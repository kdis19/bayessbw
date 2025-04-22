library(tidyverse)
library(parallel)
library(dunn.test)
library(ggpubr)
source('code/functions.R')

col.ord <- c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2')
col.ord.lab <- c('Lab-Reared', 'Northwest Territories', 'Alberta', 'Ontario',
                 'Québec', 'New Brunswick')


all.data <- read.csv('data/all_days_df.csv')
all.data <- all.data %>% 
  rename(colony = province) %>%
  subset(select = -X) %>%
  # mutate(colony = factor(colony, levels = col.ord, labels = col.ord.lab),
  #        dd = time1*temp1 + time2*temp2)
  mutate(colony = factor(colony, levels = col.ord),
         dd = time1*temp1 + time2*temp2,
         colony.full = factor(colony, levels = col.ord, labels = col.ord.lab))


posterior.orig <- read.csv('code/output/all_provs_clean.csv') %>%
  rename('temp1' = treatment)
# posterior <- read.csv('code/output/all_provs_clean.csv')

t.time.df <- all.data %>% 
  group_by(stage, temp1) %>%
  summarise(t.time = which.max(tabulate(time1))) %>%
  mutate(opt.temp = ifelse(temp1 %in% 15:25, 0, 1),
         t.time = t.time*opt.temp) %>%
  #subset(!(temp1 %in% 15:25)) %>%
  subset(select = -opt.temp) %>%
  arrange(stage, temp1)

s.all.data <- all.data %>%
  subset(stage == 'L2') %>%
  group_by(colony, temp1) %>%
  summarise(count = sum(nobs))

posterior <- posterior.orig %>% 
  merge(s.all.data) %>%
  merge(t.time.df)


age.fun <- function(df) {
  
  age <- 0
  time1 <- vector()
  time2 <- vector()
  extra <- 0
  
  # u.temp <- unique(df$temp)
  # o.temp <- u.temp %in% 15:25
  
  for (i in 1:nrow(df)) {
    row <- df[i,]
    
    u.temp <- row$temp
    o.temp <- u.temp %in% 15:25
    
    extra.dev <- extra*ifelse(o.temp, row$new.rate, row$new.r20)
    #print(extra.dev)
    treat.dev <- cumsum(c(extra.dev, rep(row$new.rate, row$t.time)))
    #print(treat.dev)
    if (any(treat.dev > 1)) {
      ind <- min(which(treat.dev > 1))
      extra <- treat.dev[ind] - 1
      time1 <- c(time1, ind - 1)
      time2 <- c(time2, 0)
      age <- 0
    }
    else {
      time1 <- c(time1, row$t.time)
      age <- row$t.time*row$new.rate
      rate <- ifelse(o.temp, row$new.rate, row$new.r20)
      time <- (1 - age)/rate
      c.time <- ceiling(time)
      extra <- c.time - time
      time2 <- c(time2, c.time)
    }
  }
  new.df <- subset(df, select = c(colony, stage, cup, index))
  new.df <- df %>%
    mutate(temp1 = temp,
           temp2 = ifelse(o.temp, u.temp, 20),
           transfer = !o.temp) %>%
    subset(select = c(colony, stage, cup, temp1, temp2)) %>%
    mutate(time1 = time1, time2 = time2)
  return(new.df)
    
}

single.draw <- function(ind, post.in = posterior) {
  post <- post.in %>%
    subset(index == ind) %>%
    mutate(#curve.index = 1:nrow(.),
           gp.index = paste0(colony, '_', stage),
           curve.index = as.numeric(factor(gp.index))) %>%
    rename(treatment = temp1) %>%
    arrange(curve.index, treatment)
  
  print('get curves')
  
  upi <- unique(post$gp.index)
  post.lst <- lapply(seq(upi), function(ind) {
    x <- upi[ind]
    sub <- subset(post, gp.index == x)
    #sub$curve.index <- ind
    return(sub)
  })
  # curves <- get_curves_pp(post.lst, spline = TRUE,
  #                         temp = seq(5, 35, by = 5)) %>%
  curves <- get_curves_pp(post.lst, spline = TRUE) %>%
    subset(temp %in% seq(5, 35, by = 5)) %>%
    rename(curve.index = index) %>%
    group_by(curve.index) %>%
    mutate(r20 = rate[temp == 20]) %>%
    rename(treatment = temp)
  
  print('cycle through individuals')
  eps.lst <- lapply(1:nrow(post), function(r) {
    row <- post[r,]
    epsilon <- rnorm(row$count, 0, row$s_eps)
    delta <- exp(epsilon)
    new.row <- row %>%
      select(curve.index, colony, treatment, stage, t.time)
    df <- data.frame(delta = delta,
                     cup = 1:(row$count)) %>%
      merge(new.row)
    # df <- data.frame(curve.index = row$curve.index,
    #                  colony = row$colony,
    #                  temp = row$treatment,
    #                  stage = row$stage,
    #                  t.time = row$t.time,
    #                  delta = delta,
    #                  cup = 1:(row$count))
    return(df)
  })
  eps.df <- bind_rows(eps.lst)
  
  print('merge data frames')
  ind.curves <- curves %>%
    merge(eps.df) %>%
    mutate(new.rate = delta*rate,
           new.r20 = delta*r20,
           index = paste(colony, treatment, cup, sep = '_')) %>%
    arrange(colony, treatment, cup, stage) %>%
    rename(temp = treatment)
  
  print('generate data')
  new.data <- ind.curves %>%
    group_by(index) %>%
    age.fun() %>%
    # count(colony, stage, temp1, temp2, time1, time2) %>%
    # rename(nobs = n) %>%
    mutate(temp2 = ifelse(time1 == 0, temp1, temp2),
           dd = temp1*time1 + temp2*time2)
  
  return(new.data)
}

post.sub <- subset(posterior, index <= 50)

start.time <- proc.time()

set.seed(857)
s100 <- sample(1:length(unique(posterior$index)), 150)

cl <- makeCluster(50)
clusterEvalQ(cl, {
  library(tidyverse)
  source('code/functions.R')
})
clusterExport(cl, c('posterior', 'age.fun', 'single.draw', 's100'))

data.gen.lst <- parLapply(cl, 1:150, function(i) {
  sd <- single.draw(ind = s100[i], post.in = posterior)
  sd$iteration <- i
  return(sd)
})
data.gen.df <- bind_rows(data.gen.lst)

stopCluster(cl)

end.time <- proc.time()
full.time <- end.time - start.time

col.ord <- c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2')
col.ord.lab <- c('Lab-Reared', 'Northwest Territories', 'Alberta', 
                 'Ontario', 'Québec', 'New Brunswick')

dg.prov.lst <- lapply(col.ord, function(col) {
  sub <- subset(data.gen.df, colony == col)
  write.csv(sub, paste0('code/output/datagen/dg_', col, '.csv'), 
            row.names = FALSE)
})
write.csv(data.gen.df, 'code/output/datagen/dg_all.csv', row.names = FALSE)


data.gen.gp <- data.gen.df %>% 
  count(colony, stage, temp1, dd, name = 'nobs') %>%
  mutate(colony.full = factor(colony, levels = col.ord, labels = col.ord.lab))

 
dgd.cs <- data.gen.df %>%
  group_by(colony, cup, temp1) %>%
  arrange(stage) %>%
  reframe(dd = cumsum(dd))

cbP <- c("#E69F00", "#56B4E9",  "#F0E442", "#009E73",
               "#0072B2", "#D55E00")

ggplot(data = data.gen.gp) +
  geom_violin(aes(x = factor(temp1), y = dd, fill = colony.full, 
                  col = colony.full, weight = nobs), alpha = 0.8) +
  geom_point(data = all.data, aes(x = factor(temp1), y = dd,
                                  size = nobs), alpha = 0.3) +
  facet_grid(rows = vars(stage), cols = vars(colony), scales = 'free') +
  scale_fill_manual(values = cbP) +
  scale_color_manual(values = cbP) +
  guides(fill = 'none', col = 'none') +
  labs(size = '# of Obs', y = expression('Degree Days ' (degree~C)),
       x = expression('Rearing Temperature ' (degree~C))) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, 
                                    margin = margin(t = 0.5, unit = 'cm')),
        axis.title.y = element_text(size = 14, 
                                    margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 13))

pars <- c('HA', 'TL', 'HL', 'TH', 'HH')

post.df <- readRDS('code/output/stanfit_sups_st.rds') %>%
  as.data.frame()

plots <- list()
for (p in pars) {
  sub <- post.df[,grep(p, names(post.df))]
  sub2 <- sub[,2:length(sub)]
  names(sub2) <- provs
  
  m.sub <- melt(sub2)
  levels(m.sub$variable) <- provs
  m.sub$variable <- factor(m.sub$variable, levels = c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2'))
  gg <- ggplot(data = m.sub) +
    geom_violin(aes(x = variable, y = value, fill = variable, col = variable),
                draw_quantiles = 0.5, alpha = 0.5) +
    scale_fill_manual(values = cbP) +
    scale_color_manual(values = cbP) +
    theme_minimal() +
    theme(legend.position = 'none') +
    geom_hline(yintercept = 1) +
    labs(title = p)
  plots <- append(plots, list(gg))
  print(dunn.test(as.list(sub2)))
  print(round(apply(sub2, 2, mean), 3))
  print(apply(sub2, 2, 
              function(x) {
                wt <- wilcox.test(x, mu = 1, alternative = 'two.sided')
                round(wt$p.value, 2)}))
}
ggarrange(plotlist = plots, nrow = 2, ncol = 3)
