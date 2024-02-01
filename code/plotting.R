library(tidyverse)
library(GGally)
library(lubridate)
library(ggpubr)
source('code/likelihood.R')

prov.ord <- c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2')
prov.ord.lab <- c('Lab-Reared', 'Northwest Territories', 'Alberta', 
                  'Ontario', 'Québec', 'New Brunswick')
prov.ord.lab.short <- c('IPU', 'NWT', 'AB', 'ON', 'QC', 'NB')
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
dg.df$province <- factor(dg.df$province, levels = prov.ord,
                         labels = prov.ord.lab)
dg.df$stage <- stages[dg.df$stage + 1]

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
  gc$quantile <- qs[gc$index]
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

a.data$province <- factor(a.data$province, levels = prov.ord,
                          labels = prov.ord.lab)
a.data$stage <- stages[a.data$stage + 1]

ggplot(data = dg.df) +
  geom_violin(aes(x = factor(temp1), y = dd, weight = nobs,
                  fill = province, col = province)) +
  geom_point(data = a.data, 
             aes(x = factor(temp1), y = dd, size = nobs), 
             alpha = 0.3) +
  facet_grid(cols = vars(province), rows = vars(stage), scales = 'free') +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  guides(fill = 'none', col = 'none') +
  labs(size = '# of Obs', y = expression('Degree Days ' (degree~C)),
       x = expression('Rearing Temperature ' (degree~C))) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13))

##### Figure 5 #####
## Cross-Validation results ##
s.all.df <- read.csv('code/output/cv/s_col_compare.csv')
elpd.temp.df <- read.csv('code/output/cv/elpd_temp.csv')

self.df <- subset(elpd.temp.df, prov.orig == prov.comp)
s.prov.df <- subset(s.all.df, prov.comp == prov.orig)

s.all.df$prov.comp <- factor(s.all.df$prov.comp, levels = prov.ord,
                             labels = prov.ord.lab.short)
s.all.df$prov.orig <- factor(s.all.df$prov.orig, levels = prov.ord,
                             labels = prov.ord.lab)

s.prov.df$prov.comp <- factor(s.prov.df$prov.comp, levels = prov.ord,
                              labels = prov.ord.lab.short)
s.prov.df$prov.orig <- factor(s.prov.df$prov.orig, levels = prov.ord,
                              labels = prov.ord.lab)

elpd.temp.df$prov.comp <- factor(elpd.temp.df$prov.comp, levels = prov.ord,
                                 labels = prov.ord.lab)
elpd.temp.df$prov.orig <- factor(elpd.temp.df$prov.orig, levels = prov.ord,
                             labels = prov.ord.lab)
self.df$prov.comp <- factor(self.df$prov.comp, levels = prov.ord,
                                 labels = prov.ord.lab)
self.df$prov.orig <- factor(self.df$prov.orig, levels = prov.ord,
                              labels = prov.ord.lab)

ggplot(data = s.all.df) +
  geom_hline(data = s.prov.df, aes(yintercept = mean)) +
  geom_point(aes(x = prov.comp, y = mean, col = prov.comp), size = 3) +
  geom_segment(aes(x = prov.comp, xend = prov.comp, 
                   y = mean - se, yend = mean + se,
                   col = prov.comp)) +
  facet_wrap(vars(prov.orig)) +
  scale_y_continuous(trans = 'exp') +
  scale_color_manual(values = cbPalette) +
  theme_minimal() +
  labs(x = 'Posterior Colony', y = 'Mean Individual ELPD') +
  theme(legend.position = 'none', 
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 14))

ggplot(data = elpd.temp.df) +
  geom_line(data = self.df, aes(x = temp, y = mean, col = prov.orig), 
            show.legend = FALSE) +
  geom_point(aes(x = temp, y = mean, col = prov.comp), size = 2) +
  geom_segment(aes(x = temp, xend = temp, y = mean - se, yend = mean + se, 
                   col = prov.comp)) +
  facet_wrap(vars(prov.orig)) +
  scale_color_manual(values = cbPalette) +
  theme_minimal() +
  labs(x = 'Posterior Colony', y = 'Mean Individual ELPD', col = 'Data Colony') +
  theme(strip.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13))

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
c.post.l2$province <- factor(c.post.l2$province, levels = c('Prior', prov.ord),
                             labels = c('Prior', prov.ord.lab.short))
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

cp.sw.ag$province <- factor(cp.sw.ag$province, levels = prov.ord,
                            labels = prov.ord.lab.short)
v2labs <- c('`Multiplicative Intercept` ~~ (rho)', 
            '`Scale of Individual Variation` ~~ (sigma[epsilon])')
cp.sw.ag$var2 <- factor(factor(cp.sw.ag$variable), 
                        labels = v2labs)


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

par.lst <- lapply(provs, function(p) {
  post <- read.csv(paste0('code/output/', p, '_post.csv'))
  w <- which.max(post$`lp__`)
  sub <- subset(c.post, province == p & index == w)
  return(sub)
})
par.df <- bind_rows(par.lst)
par.df2 <- subset(par.df,
                  select = c(rho, HA, TL, HL, TH, HH, 
                             province, stage))
gc2 <- get_curves(par.df2)
gc2$province <- sapply(gc2$index, function(x) {par.df2$province[x]})
gc2$stage <- sapply(gc2$index, function(x) {par.df2$stage[x]})
gc2$province <- factor(gc2$province, levels = prov.ord,
                       labels = prov.ord.lab)

ggplot(data = gc2) +
  geom_line(aes(x = temp, y = rate, col = province), linewidth = 1.25) +
  theme_minimal() +
  scale_color_manual(values = cbPalette) +
  facet_wrap(vars(stage)) +
  labs(x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate', col = 'Colony') +
  theme(strip.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.position = c(0.8, 0.28))

##### Figure 8 #####
pdates <- read.csv('code/output/pupal_dates.csv')
pdates$province <- factor(pdates$province, levels = prov.ord,
                          labels = prov.ord.lab)

ggplot(data = pdates) +
  geom_boxplot(aes(x = province, y = days(pup.jd) + as.Date('2018-12-31'), 
                   fill = province, col = province), 
               alpha = 0.5, linewidth = 1.5, size = 2) +
  labs(x = 'Colony', y = 'Estimated Date of Pupation') +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 12)) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)
  
## Currently only have ON 2019 and IN 1997 ##
wp.fun <- function(Province, year, tags = c('A', 'B', 'C')) {
  province <- tolower(Province)
  yr <- substr(as.character(year), 3, 4)
  
  
  on.weather <- read.csv(paste0('data/', province, yr, '_real_weather.csv'))
  on.weather$DATE <- as.Date(on.weather$DATE)
  on.weather$datetime <- on.weather$DATE + hours(on.weather$Hour)
  on.weather <- subset(on.weather, between(Month, 4, 7))
  
  daily <- aggregate(data = on.weather, Temp ~ DATE, mean)
  daily$DATE <- as.POSIXct(daily$DATE)
  
  on.ages.orig <- read.csv(paste0('code/output/weathersims/', Province, 
                                  yr, '_sim_age.csv'))
  on.ages.ag <- aggregate(data = on.ages.orig, age ~ DATE + province, 
                          function(x) {quantile(x, probs = c(0.05, 0.25, 0.5, 
                                                             0.75, 0.95))})
  age.mat <- as.data.frame(on.ages.ag$age)
  names(age.mat) <- c('lowest', 'low', 'mid', 'high', 'highest')
  on.ages <- bind_cols(on.ages.ag[,c('DATE', 'province')], age.mat)
  on.ages$DATE <- as.Date(on.ages$DATE)
  on.ages <- subset(on.ages, between(month(DATE), 4, 7))
  on.ages$province <- factor(on.ages$province, levels = prov.ord, 
                             labels = prov.ord.lab.short)
  
  on19 <- ggplot(data = on.weather) +
    geom_line(aes(x = datetime, y = Temp), alpha = 0.3) +
    geom_line(data = daily, aes(x = DATE, y = Temp), linewidth = 1.5) +
    theme_minimal() +
    labs(x = 'Date', y = expression('Temperature ' (degree~C)),
         tag = tags[1], title = paste0(Province, ' ', year)) +
    theme(axis.title = element_text(size = 14),
          plot.tag = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold')) +
    ylim(c(0, 35))
  
  ages <- ggplot(data = on.ages) +
    # geom_ribbon(aes(x = DATE, ymin = lowest, ymax = highest,
    #                 fill = province), alpha = 0.3) +
    geom_ribbon(aes(x = DATE, ymin = low, ymax = high,
                    fill = province), alpha = 0.7) +
    geom_line(aes(x = DATE, y = mid, col = province)) +
    scale_fill_manual(values = cbPalette) +
    scale_color_manual(values = cbPalette) +
    scale_y_continuous(labels = c(paste0('L', 2:6), 'Pupa')) +
    labs(col = 'Colony', fill = 'Colony', x = 'Date', 
         y = 'Larval Stage', tag = tags[2]) +
    theme_minimal() +
    theme(axis.title = element_text(size = 14),
          legend.title = element_text(size = 13),
          #legend.position = 'bottom',
          legend.position = 'none',
          plot.tag = element_text(size = 14),
          axis.text = element_text(size = 12)) #+
    # guides(col = guide_legend(nrow = 1),
    #        fill = guide_legend(nrow = 1))
  
  pdates <- read.csv(paste0('code/output/weathersims/', Province, yr,
                            '_pupal_dates.csv'))
  pdates$province <- factor(pdates$province, levels = prov.ord,
                            labels = prov.ord.lab)
  
  pdplot <- ggplot(data = pdates) +
    geom_boxplot(aes(x = province, 
                     y = days(pup.jd) + as.Date(paste0(year - 1, '-12-31')), 
                     fill = province, col = province), 
                 alpha = 0.5, linewidth = 1.5, size = 2) +
    labs(x = 'Colony', y = 'Pupation Date', tag = tags[3]) +
    theme_minimal() +
    theme(legend.position = 'none',
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 12),
          plot.tag = element_text(size = 14)) +
    scale_color_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette)
  
  return(list(on19, ages, pdplot))
}

on.lst <- wp.fun('ON', 2019, tags = c('A', 'C'))
in.lst <- wp.fun('IN', 1997, tags = c('B', 'D'))

ggarrange(on.lst[[1]], in.lst[[1]], on.lst[[2]], in.lst[[2]],
          nrow = 2, ncol = 2, heights = c(0.3, 0.7), 
          common.legend = TRUE, legend = 'bottom')
