library(tidyverse)
library(GGally)
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
params <- list('HA' = qgamma(0.5, 5.4, scale = 0.134),
               'HL' = -qgamma(0.5, 3.6, scale = 2.253),
               'HH' = qgamma(0.5, 7.6, scale = 3.12),
               'TL' = tl,
               'TH' = tl + tdiff,
               'rho' = qgamma(0.5, median(k), scale = 0.045))
gc <- get_curves(as.data.frame(params))

tdiff <- qgamma(0.5, 112, scale = 0.226)
tl <- exp(5.67)
sm1 <- 0:4
k <- -1.05*sm1^2 + 4.22*sm1 + 4.08
params <- list('HA' = qgamma(0.5, 5.4, scale = 0.134),
               'HL' = -qgamma(0.5, 4, scale = 2),
               'HH' = qgamma(0.5, 5, scale = 0.5),
               'TL' = tl,
               'TH' = tl + tdiff,
               'rho' = qgamma(0.5, median(k), scale = 0.045))
gc <- get_curves(as.data.frame(params))


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
  labs(size = '# of Obs', y = 'Degree Days',
       x = expression('Rearing Temperature ' (degree~C))) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 13))

##### Figure 5 #####
##### Figure 6 #####
c.post.l2 <- subset(c.post, stage == 'L2',
                    select = c('HA', 'TL', 'HL', 'TH', 'HH', 'province'))
c.post.l2$province <- factor(c.post.l2$province, levels = prov.ord)
ggpairs(data = c.post.l2, columns = 1:5, 
        # upper = list(continuous = 'blank'),
        upper = list(continuous = 'points'),
        aes(fill = province, colour = province, alpha = 0.4)) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 14))


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



