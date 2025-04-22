library(tidyverse)
library(reshape2)
library(factoextra)
library(ggpubr)
library(GGally)
library(dunn.test)

source('code/functions.R')
source('code/prior_sampling.R')

prov.ord <- c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2')
prov.ord.lab <- c('Lab-Reared', 'Northwest Territories', 'Alberta', 
                  'Ontario', 'Québec', 'New Brunswick')
prov.ord.lab.short <- c('IPQL', 'NWT', 'AB', 'ON', 'QC', 'NB')

c.post <- read.csv('code/output/all_provs_clean.csv') %>%
  #rename(colony = 'province') %>%
  mutate(colony = factor(colony, levels = prov.ord, 
                         labels = prov.ord.lab.short))
all.data <- read.csv('data/all_days_df.csv') %>%
  select(-X) %>%
  rename(colony = province) %>%
  mutate(colony = factor(colony, levels = prov.ord, labels = prov.ord.lab),
         dd = time1*temp1 + time2*temp2)

pars.col <- c('HA', 'TL', 'HL', 'TH', 'HH')
ps.lst <- lapply(prior.samp.reduced(3000), function(x) {
  as.data.frame(x[pars.col])
})
ps.df <- bind_rows(ps.lst)
prior.df <- ps.df %>%
  mutate(TA = ta.fun(HL, HH, TL, TH) - 273.15,
         TL = TL - 273.15,
         TH = TH - 273.15) %>%
  mutate(colony = 'Prior') %>%
  pivot_longer(!colony, names_to = 'parameter')

cbP <- c("#E69F00", "#56B4E9",  "#F0E442", 
         "#009E73", "#0072B2", "#D55E00")
cbPalette <- c('black', 'grey', cbP)
cbP2 <- c('black', cbP)

##### Population Counts #####
pop.counts <- all.data %>%
  mutate(colony = factor(colony, 
                         levels = c('Northwest Territories', 'Ontario',
                                    'Québec', 'Alberta', 'Lab-Reared', 
                                    'New Brunswick'))) %>%
  subset(stage == 'L6') %>%
  group_by(temp1, colony) %>%
  summarize(nobs = sum(nobs))

ggplot(data = pop.counts) +
  geom_bar(aes(x = factor(temp1), y = nobs, fill = colony), 
           stat = 'identity') +
  geom_text(aes(x = factor(temp1), y = nobs + 15, label = nobs), size = 3.75) +
  scale_fill_manual(values = cbP[c(2, 4, 5, 3, 1, 6)]) +
  facet_wrap(vars(colony)) +
  theme_minimal() +
  labs(x = expression('Treatment Temperature ' (degree~C)),
       y = 'Number of Individuals') +
  theme(legend.position = 'none',
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 13, 
                                    margin = margin(r = 0.5, unit = 'cm')),
        axis.title.x = element_text(size = 13, 
                                    margin = margin(t = 0.5, unit = 'cm')),
        axis.text = element_text(size = 11),
        panel.spacing = unit(15, units = 'pt')) +
  ylim(0, 250)

##### Schoolfield marginal boxplots #####
final.bp <- c.post %>%
  subset(stage == 'L2', select = c(HA, TL, HL, TH, HH, colony)) %>%
  mutate('TA' = ta.fun(HL, HH, TL, TH) - 273.15,
         'HL' = -abs(HL),
         'TL' = TL - 273.15,
         'TH' = TH - 273.15) %>%
  pivot_longer(!colony, names_to = 'parameter', cols_vary = 'slowest') %>%
  bind_rows(prior.df) %>%
  mutate(parameter = factor(parameter, 
                           levels = c('HL', 'HA', 'HH',
                                      'TL', 'TA' ,'TH'),
                           labels = c('H[L]', 'H[A]', 'H[H]', 
                                      'T[L]', 'T[A]', 'T[H]')),
         colony = factor(colony, levels = c('Prior', prov.ord.lab.short)))


ggplot(data = final.bp) +
  geom_boxplot(aes(x = colony, y = value, fill = colony)) +
  scale_fill_manual(values = cbP2) +
  facet_wrap(vars(parameter), scales = 'free', labeller = 'label_parsed') +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 15, 
                                  margin = margin(b = 10)),
        axis.title.x = element_text(size = 14,
                                    margin = margin(t = 10)),
        axis.text = element_text(size = 11),
        panel.spacing = unit(10, units = 'pt')) +
  guides(fill = 'none') +
  labs(fill = 'Colony', col = 'Colony', y = '', x = 'Colony')


##### Colony x stage-wise boxplots of sigma epsilon and rho #####
sw.post <- c.post %>%
  select(rho, s_eps, stage, index, colony) %>%
  melt(id.vars = c('stage', 'index', 'colony'),
       variable.name = 'parameter', value.name = 'Value') %>%
  mutate(parameter = factor(parameter, levels = c('rho', 's_eps'), 
                            labels = c('rho', 'sigma[epsilon]')))

ggplot(data = sw.post) +
  geom_boxplot(aes(x = colony, y = Value, fill = colony)) +
  scale_fill_manual(values = cbP) +
  facet_grid(rows = vars(parameter), cols = vars(stage), scales = 'free',
             labeller = 'label_parsed') +
  labs(y = 'Parameter Value', x = '', fill = 'Colony') +
  theme_minimal() +
  theme(strip.text.y = element_text(size = 18, face = 'bold', 
                                    margin = margin(l = 20, r = 15), angle = 0),
        strip.text.x = element_text(size = 12, margin = margin(b = 10)),
        axis.title.y = element_text(size = 12, margin = margin(r = 10)),
        axis.text.x = element_text(size = 8.5, angle = 45),
        legend.title = element_text(size = 12, margin = margin(b = 10)),
        panel.spacing = unit(15, units = 'pt'))

##### PCA Plots #####
pc.post <- c.post %>%
  select(all_of(c('colony', pars.col))) %>%
  distinct()
pc.post.rho <- c.post %>%
  select(c(all_of(c('colony', 'stage', pars.col)),
           contains('rho'))) %>%
  distinct() %>%
  mutate(stage.rho = paste0('rho', substr(stage, 2, 2))) %>%
  dcast(colony + HA + TL + HL + TH + HH ~ stage.rho,
        value.var = 'rho')

pca.1 <- prcomp(select(pc.post.rho, c(all_of(pars.col),
                                  contains('rho'))), scale. = TRUE)
p1 <- fviz_pca_biplot(pca.1,
                      axes=c(1,2),
                      col.ind = pc.post$colony,
                      #select.ind = list(cos2 = 0.90),
                      select.ind = list(cos2 = 0.8),
                      geom = c("point"),
                      col.var = 'black',
                      addEllipses = TRUE,
                      ellipse.level = 0.95) +
  scale_color_manual(values = cbP) +
  scale_fill_manual(values = cbP) +
  guides(fill = 'none', col = 'none', shape = 'none') +
  scale_shape_manual(values =rep(19, 6)) +
  labs(title = 'Dimension 1 vs. Dimension 2')
p2 <- fviz_pca_biplot(pca.1,
                      axes=c(1,3),
                      col.ind = pc.post$colony,
                      select.ind = list(cos2 = 0.80),
                      geom = c("point"),
                      col.var = 'black',
                      addEllipses = TRUE,
                      ellipse.level = 0.95) +
  scale_color_manual(values = cbP) +
  scale_fill_manual(values = cbP) +
  scale_shape_manual(values = rep(19, 6)) +
  labs(col = 'Colony', fill = 'Colony', 
       title = 'Dimension 1 vs. Dimension 3') +
  guides(shape = 'none')
ggarrange(p1,p2, widths = c(0.43, 0.57))

##### Posterior pairs #####
pp.post <- ps.df %>%
  mutate(colony = 'Prior') %>%
  relocate(colony) %>%
  bind_rows(pc.post) %>%
  mutate(colony = factor(colony, levels = c('Prior', prov.ord.lab.short)),
         TH = TH - 273.15,
         TL = TL - 273.15,
         TA = ta.fun(HL, HH, TL, TH))

pars.col.ps <- sapply(c(pars.col, 'TA'), function(x) {
  c1 <- substr(x, 1, 1)
  c2 <- substr(x, 2, 2)
  paste0(c1, '[', c2, ']')
})
names(pp.post)[2:length(pp.post)] <- pars.col.ps
ggpairs(data = pp.post, columns = 2:7, 
        # upper = list(continuous = 'blank'),
        upper = list(continuous = wrap('points', size = 1, alpha = 0.5)),
        diag = list(continuous = wrap('densityDiag', alpha = 0.5)),
        #diag = list(continuous = wrap('barDiag')),
        lower = list(continuous = wrap('points', size = 1, alpha = 0.5)),
        labeller = 'label_parsed',
        legend = c(1, 1),
        aes(fill = colony, colour = colony)) +
  scale_fill_manual(values = cbP2) +
  scale_colour_manual(values = cbP2) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold', size = 14),
        legend.position = 'right') +
  labs(fill = 'Colony', col = 'Colony')

##### Median dev curves NO SPLINE #####
gc.df <- read.csv('code/output/all_parametric_curves.csv')

gc.ag.df <- gc.df %>%
  group_by(temp, stage, colony) %>%
  summarize(rate = median(rate)) %>%
  mutate(colony = factor(colony, levels = prov.ord.lab))

gc.ag2.df <- gc.ag.df %>%
  subset(colony != 'Lab-Reared') %>%
  group_by(temp, stage) %>%
  summarize(rate.range = max(rate) - min(rate)) %>%
  mutate(stage.ind = as.numeric(factor(stage)) + 1)

ggplot(data = gc.ag2.df) +
  geom_line(aes(x = temp, y = rate.range, group = stage,
                col = stage.ind), linewidth = 1.25) +
  theme_minimal() +
  labs(x = expression('Rearing Temperature ' (degree~C)),
       y = 'Range of Development Rates (% of stage per day)',
       col = 'Larval Stage') +
  scale_color_gradient(breaks = 2:6, labels = paste0('L', 2:6)) +
  theme(axis.title.x = element_text(size = 14, 
                                     margin = margin(t = 0.5, unit = 'cm')),
         axis.title.y = element_text(size = 14, 
                                     margin = margin(r = 0.5, unit = 'cm')),
         legend.title = element_text(size = 13,
                                     margin = margin(b = 0.5, unit = 'cm')),
         axis.text = element_text(size = 12),
         legend.text = element_text(size = 12))

ggplot(data = gc.ag.df) +
  geom_line(aes(x = temp, y = rate, col = colony), linewidth = 1.25) +
  scale_color_manual(values = cbP) +
  facet_wrap(vars(stage)) +
  theme_minimal() +
  labs(x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate (% of stage per day)', col = 'Colony') +
  theme(strip.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 14, 
                                    margin = margin(t = 0.5, unit = 'cm')),
        axis.title.y = element_text(size = 14, 
                                    margin = margin(r = 0.5, unit = 'cm')),
        legend.position = 'inside',
        legend.position.inside = c(0.8, 0.28))

##### Curve Features #####
lowdev.t <- 1/30
final.s <- gc.df %>%
  mutate(dev.ind = ifelse(rate > lowdev.t, 1, 0)) %>%
  group_by(index, colony, stage) %>%
  summarise(area = sum(0.1*rate),
            optimal = temp[which.max(rate)],
            max.rate = max(rate),
            temp.range = sum(dev.ind)/10) %>%
  melt(id.vars = c('index', 'colony', 'stage'),
       variable.name = 'feature', value.name = 'value') %>%
  mutate(#colony = factor(colony, 
    #levels = c.ab.levels,
    #labels = c.ab.levels2),
    feature = factor(feature, levels = c('area', 'temp.range',
                                         'optimal', 'max.rate'),
                     labels = c('Area', 'Temperature Range', 
                                'Optimal Temperature',
                                'Max Development Rate')),
    colony = factor(colony, levels = prov.ord.lab, labels = prov.ord.lab.short)) 


## Facet wrap (L4 only) ##

ggplot(data = subset(final.s, stage == 'L4')) +
  geom_boxplot(aes(x = colony, y = value, fill = colony)) +
  scale_fill_manual(values = cbP) +
  facet_wrap(vars(feature), scales = 'free') +
  #facet_grid(rows = vars(feature), scales = 'free_y', cols = vars(stage)) +
  theme_minimal() +
  labs(x = '', y = '', fill = 'Colony') +
  guides(fill = 'none') +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12))

## Facet grid (all stages) ##
ggplot(data = final.s) +
  geom_boxplot(aes(x = colony, y = value, fill = colony)) +
  scale_fill_manual(values = cbP) +
  facet_grid(rows = vars(feature), cols = vars(stage), scales = 'free') +
  #facet_grid(rows = vars(feature), scales = 'free_y', cols = vars(stage)) +
  theme_minimal() +
  labs(x = '', y = '', fill = 'Colony') +
  guides(fill = 'none') +
  theme(strip.text.y = element_text(size = 12, 
                                    margin = margin(l = 0.5, unit = 'cm')),
        strip.text.x = element_text(size = 14, 
                                    margin = margin(b = 0.5, unit = 'cm')),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12))


##### Median dev curves SPLINE #####
gcs.df <- read.csv('code/output/all_splined_curves.csv')

gcs.ag.df <- gcs.df %>%
  group_by(temp, stage, colony) %>%
  summarize(rate = median(rate)) %>%
  mutate(colony = factor(colony, levels = prov.ord.lab))

ggplot(data = gcs.ag.df) +
  geom_line(aes(x = temp, y = rate, col = colony), linewidth = 1.25) +
  scale_color_manual(values = cbP) +
  facet_wrap(vars(stage)) +
  theme_minimal() +
  labs(x = expression('Rearing Temperature ' (degree~C)),
       y = 'Development Rate (% of stage per day)', col = 'Colony') +
  theme(strip.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.title.x = element_text(size = 14, 
                                    margin = margin(t = 0.5, unit = 'cm')),
        axis.title.y = element_text(size = 14, 
                                    margin = margin(r = 0.5, unit = 'cm')),
        legend.position = 'inside',
        legend.position.inside = c(0.8, 0.28))


##### Goodness of fit violin plots #####
data.gen.df <- read.csv('code/output/datagen/dg_all.csv')
data.gen.gp <- data.gen.df %>% 
  count(colony, stage, temp1, dd, name = 'nobs') %>%
  mutate(colony = factor(colony, levels = prov.ord, labels = prov.ord.lab))
  
ggplot(data = data.gen.gp) +
  geom_violin(aes(x = factor(temp1), y = dd, fill = colony, 
                  col = colony, weight = nobs), alpha = 0.8) +
  geom_point(data = all.data, aes(x = factor(temp1), y = dd,
                                  size = nobs), alpha = 0.3) +
  facet_grid(rows = vars(stage), cols = vars(colony), scales = 'free') +
  scale_fill_manual(values = cbP) +
  scale_color_manual(values = cbP) +
  guides(fill = 'none', col = 'none') +
  labs(size = '# of Obs', y = expression('Degree Days ' (degree~C)),
       x = expression('Rearing Temperature ' (degree~C))) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, 
                                    margin = margin(t = 0.5, unit = 'cm')),
        axis.title.y = element_text(size = 14, 
                                    margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 13))

##### Dunn test/comparison of multipliers #####
pars <- c('HA', 'TL', 'HL', 'TH', 'HH')
prov.run.ord <- rev(sort(prov.ord))
post.df <- read.csv('code/output/all_provs_post.csv')

plots <- list()
for (p in pars) {
  sub <- post.df[,grep(p, names(post.df))]
  sub2 <- sub[,2:length(sub)]
  names(sub2) <- prov.run.ord
  sub2 <- select(sub2, -IPU)
  sub3 <- sub2*sub[,1]

  m.sub <- melt(sub2, variable.name = 'colony') %>%
    mutate(colony = factor(colony, levels = prov.ord, 
           labels = prov.ord.lab.short),
           variable = p)
  gg <- ggplot(data = m.sub) +
    geom_violin(aes(x = colony, y = value, fill = colony, col = colony),
                draw_quantiles = 0.5, alpha = 0.5) +
    scale_fill_manual(values = cbP) +
    scale_color_manual(values = cbP) +
    theme_minimal() +
    theme(legend.position = 'none') +
    geom_hline(yintercept = 1) +
    labs(title = p)
  plots <- append(plots, list(gg))
  #print(dunn.test(as.list(sub2)))
  print(dunn.test(as.list(sub3)))
  print(round(apply(sub2, 2, mean), 3))
  print(apply(sub2, 2,
              function(x) {
                wt <- wilcox.test(x, mu = 1, alternative = 'two.sided')
                round(wt$p.value, 2)}))
}
ggarrange(plotlist = plots, nrow = 2, ncol = 3)

value.df <- bind_rows(value.lst)
dc.vals <- dcast(value.df, colony + index ~ variable) %>%
  arrange(index, colony) %>%
  mutate(TA = ta.fun(HL, HH, TL, TH))

par.df <- select(c.post, all_of(c(pars, 'TA', 'colony'))) %>%
  distinct() %>%
  subset(colony != 'IPQL')
meds <- par.df %>%
  group_by(colony) %>%
  summarize_all(median) %>%
  as.data.frame() %>%
  mutate(TL = TL - 273.15,
         TA = TA - 273.15,
         TH = TH - 273.15)
dunn.cols <- unique(par.df$colony)
dunn.lst <- lapply(c(pars, 'TA'), function(p) {
  print(p)
  c.lst <- lapply(dunn.cols, function(cl) {
    subset(par.df, colony == cl)[,p]
  })
  dunn.test(c.lst)
})

##### Weather Simulation Plots #####
wp.fun <- function(Province, year, tags = c('A', 'B', 'C')) {
  province <- tolower(Province)
  yr <- substr(as.character(year), 3, 4)
  
  
  on.weather <- read.csv(paste0('data/', province, yr, '_real_weather.csv'))
  on.weather$DATE <- as.Date(on.weather$DATE)
  on.weather$datetime <- on.weather$DATE + hours(on.weather$Hour)
  on.weather <- subset(on.weather, between(Month, 4, 7))
  
  full.weather <- read.csv(paste0('data/', Province,
                                  'real_weather_full.csv'))
  full.weather.ag <- aggregate(data = full.weather, Temp ~ date.orig,
                               function(x) {
                                 q <- quantile(x, probs = c(0.25, 0.75))
                                 r <- range(x)
                                 return(c(r, q))})
  full.weather.ag$Temp <- as.data.frame(full.weather.ag$Temp)
  names(full.weather.ag$Temp) <- c('min', 'max', 'lq', 'uq')
  full.weather.ag$date.orig <- as.Date(full.weather.ag$date.orig)
  plot.year <- unique(year(full.weather.ag$date.orig))
  
  daily <- aggregate(data = on.weather, Temp ~ DATE, mean)
  #daily$DATE <- as.POSIXct(daily$DATE)
  daily$DATE.back <- daily$DATE - years(unique(year(daily$DATE)) - plot.year)
  date.range <- range(daily$DATE.back)
  
  full.weather.ag <- subset(full.weather.ag, between(date.orig, date.range[1], 
                                                     date.range[2]))
  
  on.ages.orig <- read.csv(paste0('code/output/weathersims/', Province, 
                                  yr, '_sim_age.csv'))
  on.ages.ag <- aggregate(data = on.ages.orig, age ~ DATE + colony, 
                          function(x) {quantile(x, probs = c(0.05, 0.25, 0.5, 
                                                             0.75, 0.95))})
  age.mat <- as.data.frame(on.ages.ag$age)
  names(age.mat) <- c('lowest', 'low', 'mid', 'high', 'highest')
  on.ages <- bind_cols(on.ages.ag[,c('DATE', 'colony')], age.mat)
  on.ages$DATE <- as.Date(on.ages$DATE)
  on.ages <- subset(on.ages, between(month(DATE), 4, 7))
  on.ages$colony <- factor(on.ages$colony, levels = prov.ord.lab, 
                             labels = prov.ord.lab.short)
  
  on19 <- ggplot() +
    geom_ribbon(data = full.weather.ag, aes(x = date.orig, ymin = Temp$min, 
                                            ymax = Temp$max), fill = 'grey90') +
    geom_ribbon(data = full.weather.ag, aes(x = date.orig, ymin = Temp$lq, 
                                            ymax = Temp$uq), fill = 'grey65') +
    #geom_line(data = on.weather, aes(x = datetime, y = Temp), alpha = 0.3) +
    geom_line(data = daily, aes(x = DATE.back, y = Temp), linewidth = 1.5) +
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
                    fill = colony), alpha = 0.7) +
    geom_line(aes(x = DATE, y = mid, col = colony)) +
    scale_fill_manual(values = cbP) +
    scale_color_manual(values = cbP) +
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
                            '_pupal_dates.csv')) %>%
    mutate(colony = factor(colony, levels = prov.ord.lab),
           pup.jd = ceiling(pup.jd)) 
  # pdates$colony <- factor(pdates$colony, levels = prov.ord.lab,
  #                         labels = prov.ord.lab.short)
  
  pdplot <- ggplot(data = pdates) +
    geom_boxplot(aes(x = colony, 
                     y = days(pup.jd) + as.Date(paste0(year - 1, '-12-31')), 
                     fill = colony, col = colony), 
                 alpha = 0.5, linewidth = 1.5, size = 2) +
    labs(x = 'Colony', y = 'Pupation Date', tag = tags[3]) +
    theme_minimal() +
    theme(legend.position = 'none',
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 12),
          plot.tag = element_text(size = 14)) +
    scale_color_manual(values = cbP) +
    scale_fill_manual(values = cbP)
  ggarrange(on19, ages, pdplot, nrow = 3, heights = c(0.25, 0.5, 0.25))
  
  # return(list(on19, ages, pdplot))
}

on.lst <- wp.fun('ON', 2019, tags = c('A', 'B', 'C'))
on.lst
# in.lst <- wp.fun('IN', 1997, tags = c('B', 'D'))
# 
# ggarrange(on.lst[[1]], in.lst[[1]], on.lst[[2]], in.lst[[2]],
#           nrow = 2, ncol = 2, heights = c(0.3, 0.7), 
#           common.legend = TRUE, legend = 'bottom')

##### Pupal Date PCA #####
pup.cast <- read.csv('code/output/all_rel_pupdates.csv')
pca.pup <- prcomp(select(pup.cast, -colony), scale. = FALSE)
p1 <- fviz_pca_biplot(pca.pup,
                      axes=c(1,2),
                      col.ind = pup.cast$colony,
                      select.ind = list(cos2 = 0.8),
                      geom = c("point"),
                      geom.var = c('arrow'),
                      col.var = 'black',
                      alpha.var = 0.2,
                      addEllipses = TRUE,
                      ellipse.level = 0.95) +
  scale_color_manual(values = cbP) +
  scale_fill_manual(values = cbP) +
  guides(shape = 'none') +
  scale_shape_manual(values =rep(19, 6)) +
  labs(title = 'Dimension 1 vs. Dimension 2',
       fill = 'Colony', col = 'Colony')