library(tidyverse)
library(GGally)
library(lubridate)
source('code/likelihood.R')

prov.ord <- c('IN', 'AB', 'ON', 'QC', 'NB2', 'IPU')
cbPalette <- c("#E69F00", "#56B4E9",  "#F0E442", "#009E73",
               "#0072B2", "#D55E00")

c.post <- read.csv('code/output/post_cleaned.csv')
all.data <- read.csv('data/all_days_df.csv')
provs <- unique(all.data$province)
dg.lst <- lapply(provs, function(p) {
  dat <- read.csv(paste0('code/output/datagen/dg_', p, '.csv'))
  dat$province <- p
  return(dat)
})
dg.df <- bind_rows(dg.lst)

##### Figure 3 #####
## Illustration of how parameters change curve shape ##
tdiff <- qgamma(0.5, 112, scale = 0.226)
tl <- exp(5.67)
sm1 <- 0:4
k <- -1.05*sm1^2 + 4.22*sm1 + 4.08
qs <- seq(0.1, 0.9, by = 0.1)
params.orig <- list('HA' = qgamma(qs, 5.4, scale = 0.134),
                    'HL' = -qgamma(qs, 3.6, scale = 2.253),
                    'HH' = qgamma(qs, 7.6, scale = 3.12),
                    'TL' = qnorm(qs, 284, 2),
                    'TH' = qnorm(qs, 304, 2),
                    'rho' = qgamma(qs, median(k), scale = 0.045))
gc <- get_curves(as.data.frame(params.orig))

gc.lst <- lapply(names(params.orig), function(par) {
  params <- lapply(names(params.orig), function(p2) {
    if (p2 == par) {
      return(params.orig[[p2]])
    }
    else {
      return(params.orig[[p2]][5])
    }
  })
  names(params) <- names(params.orig)
  df <- as.data.frame(params)
  gc <- get_curves(df)
  gc$parameter <- par
  gc$quantile <- qs[gc.ha$index]
  return(gc)
})
gc.df <- bind_rows(gc.lst)
gc.df$parameter <- factor(gc.df$parameter,
                          labels = c('H[A]', 'H[H]', 'H[L]', 
                                     'rho[L4]', 'T[H]', 'T[L]'))

ggplot(data = gc.df) +
  geom_line(aes(x = temp, y = rate, group = index, 
                col = quantile), linewidth = 1.5) +
  facet_wrap(vars(parameter), labeller = 'label_parsed') +
  theme_minimal() +
  labs(x = expression('Input Temperature ' (degree~C)),
       y = 'Development Rate', col = 'Quantile of Prior') +
  theme(strip.text = element_text(face = 'bold', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13, hjust = 0.5))

##### Figure 4 #####
## Violin plots of model predictions vs observations ##
a.data.lst <- lapply(provs, function(p) {
  sub <- subset(all.data, province == p)
  a.sub <- data.adjust(sub)
  a.sub$dd <- a.sub$time1*a.sub$temp1 + a.sub$time2*a.sub$temp2
  a.sub$province <- p
  return(a.sub)
})
a.data <- bind_rows(a.data.lst)

dg.df$province <- factor(dg.df$province, levels = prov.ord)
a.data$province <- factor(a.data$province, levels = prov.ord)

dg.df$stage <- stages[dg.df$stage + 1]
a.data$stage <- stages[a.data$stage + 1]
ggplot(data = dg.df) +
  geom_violin(aes(x = factor(temp1), y = dd, weight = nobs,
                  fill = province, col = province)) +
  geom_point(data = a.data, 
             aes(x = factor(temp1), y = dd, size = nobs), 
             alpha = 0.3) +
  facet_grid(rows = vars(province), cols = vars(stage), scales = 'free') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  guides(fill = 'none', col = 'none') +
  labs(size = '# of Obs', y = expression('Degree Days ' (degree~C)),
       x = expression('Rearing Temperature ' (degree~C))) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13))

##### Figure 5 #####
prov.comp.df.all <- read.csv('code/output/cv/prov_compare.csv')
prov.df.ag <- read.csv('code/output/cv/prov_individual.csv')

prov.comp.df.all$prov <- factor(prov.comp.df.all$prov, levels = prov.ord)
prov.comp.df.all$post.prov <- factor(prov.comp.df.all$post.prov, 
                                     levels = prov.ord)

prov.df.ag$prov <- factor(prov.df.ag$prov, levels = prov.ord)

ggplot(data = prov.comp.df.all) +
  geom_hline(data = prov.df.ag, aes(yintercept = elpd, col = prov)) +
  geom_boxplot(aes(x = post.prov, y = elpd, col = post.prov)) +
  scale_color_manual(values = cbPalette) +
  guides(col = 'none') +
  facet_wrap(vars(prov)) +
  theme_minimal()

##### Figure 6 #####
cps <- c('HA', 'TL', 'HL', 'TH', 'HH')
ps.lst <- lapply(prior.samp(3000), function(ps) {
  as.data.frame(ps[cps])
})
ps.df <- bind_rows(ps.lst)
ps.df$province <- 'Prior'

c.post.l2 <- subset(c.post, stage == 'L2',
                    select = c('HA', 'TL', 'HL', 'TH', 'HH', 'province'))

c.post.l2 <- bind_rows(ps.df, c.post.l2)
c.post.l2$province <- factor(c.post.l2$province, levels = c('Prior', prov.ord))
names(c.post.l2) <- c('H[A]', 'T[L]', 'H[L]', 'T[H]', 'H[H]', 'province')

ggpairs(data = c.post.l2, columns = 1:5, 
        # upper = list(continuous = 'blank'),
        upper = list(continuous = wrap('points')),
        diag = list(continuous = wrap('densityDiag', alpha = 0.5)),
        lower = list(continuous = wrap('points')),
        labeller = 'label_parsed',
        legend = c(1, 1),
        aes(fill = province, colour = province)) +
  scale_fill_manual(values = c('black', cbPalette)) +
  scale_colour_manual(values = c('black',cbPalette)) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 14),
        legend.position = 'right') +
  labs(fill = 'Colony', col = 'Colony')


##### Figure 7 #####
vars <- c('rho', 's_eps')
vmap <- c('rho', 'sigma[epsilon]')
c.post.sw <- subset(c.post, select = c(vars, 'stage', 'province'))
cp.sw.ag.mat <- aggregate(data = c.post.sw, cbind(rho, s_eps) ~ ., 
                          function(x) {quantile(x, probs = c(0.05, 0.5, 0.95))})

inds <- cp.sw.ag.mat[,c('stage', 'province')]

var.lst <- lapply(vars, function(v) {
  mat <- as.data.frame(cp.sw.ag.mat[,v])
  names(mat) <- c('low', 'mid', 'high')
  mat$variable <- v
  mat <- cbind(inds, mat)
  return(mat)
})
cp.sw.ag <- bind_rows(var.lst)

cp.sw.ag$province <- factor(cp.sw.ag$province, levels = prov.ord)
cp.sw.ag$var2 <- factor(factor(cp.sw.ag$variable), 
                        labels = c('`Multiplicative Intercept` ~~ (rho)', 
                                   '`Scale of Individual Variation` ~~ (sigma[epsilon])'))


ggplot(data = cp.sw.ag) +
  geom_ribbon(aes(x = factor(stage),
                  ymin = low, ymax = high,
                  fill = province, group = province),
              alpha = 0.5) + 
  geom_line(aes(x = factor(stage), y = mid,
                col = province, group = province)) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  facet_wrap(vars(var2), scales = 'free', labeller = 'label_parsed') +
  theme_minimal() +
  labs(x = 'Stage', y = 'Value', col = 'Colony', fill = 'Colony') +
  theme(strip.text = element_text(face = 'bold', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13))

##### Figure 8 #####
on.weather <- read.csv('data/on_real_weather.csv')
on.weather$DATE <- as.Date(on.weather$DATE)
on.weather$datetime <- on.weather$DATE + hours(on.weather$Hour)
daily <- aggregate(data = on.weather, Temp.orig ~ DATE, mean)
daily$DATE <- as.POSIXct(daily$DATE)

pdates <- read.csv('code/output/pupal_dates.csv')


