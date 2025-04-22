library(tmbstan)
library(TMB)
library(tidyverse)
library(lme4)
library(reshape2)
library(factoextra)
library(ggpubr)

source('code/functions.R')
source('code/prior_sampling.R')

prov.plot.ord <- c('IPU', 'IN', 'AB', 'ON', 'QC', 'NB2')
provs <- rev(sort(unique(prov.plot.ord)))

stages <- paste0('L', 2:6)

post.df <- read.csv('code/output/all_provs_post.csv')

##### Clean stan output #####
stage.df <- data.frame('stage' = rep(stages, each = length(provs)),
                       'colony' = rep(provs, length(stages)))
inds <- nrow(post.df)
### Main curve parameters ###
pars.col <- c('HA', 'TL', 'HL', 'TH', 'HH')
pars.col.lst <- lapply(pars.col, function(p) {
  p.col <- post.df[,p]
  p.mult <- post.df[,grep(paste0(p, '_mult'), names(post.df))]
  
  pcol.mult <- p.mult*p.col
  pcol.mult$index <- 1:nrow(pcol.mult)
  m.par <- melt(pcol.mult, id.vars = 'index',
                variable.name = 'colony', value.name = p)
  
  m.par$colony <- provs[as.numeric(factor(m.par$colony))]
  return(m.par)
})
pars.col.df <- pars.col.lst[[1]]
for (i in 2:length(pars.col.lst)) {
  print(i)
  pars.col.df <- merge(pars.col.df, pars.col.lst[[i]])
}
m.pc.df <- merge(stage.df, pars.col.df)
m.pc.df$colony <- factor(m.pc.df$colony, levels = provs)

all.rows <- inds*length(stages)*length(provs)

temps <- seq(35, 5, by = -5)
add.spcols <- function(data, ups = FALSE) {
  data <- data %>%
    arrange(index) %>%
    mutate(stage = rep(stages, 
                       length(provs)*inds),
           colony = rep(rep(provs,
                            each = length(stages)),
                        inds)) %>%
    subset(select = -block)
}
### rho parameter ###
w.rho <- grep('rho', names(post.df))
w.rho.mult <- grep('rho_mult', names(post.df))

rho.df <- post.df[,w.rho[!(w.rho %in% w.rho.mult)]]
rho.df$index <- 1:nrow(rho.df)

rho.mult <- post.df[,w.rho.mult]
rho.mult$index <- 1:nrow(rho.mult)

m.rho <- melt(rho.df, id.vars = 'index',
              variable.name = 'stage', value.name = 'rho.orig')
m.rho$stage <- stages[as.numeric(factor(m.rho$stage))]

m.rho.mult <- melt(rho.mult, id.vars = 'index',
                   variable.name = 'block', value.name = 'rho.mult')
m.rho.mult <- m.rho.mult %>% add.spcols

### s_eps parameter ###
seps.df <- post.df[,grep('s_eps', names(post.df))]
seps.df$index <- 1:inds

m.seps.df <- melt(seps.df, id.vars = 'index',
                  variable.name = 'block', value.name = 's_eps')
m.seps.df <- m.seps.df %>% add.spcols

### upsilon parameter ###
sups.names <- names(post.df)[grep('s_upsilon', names(post.df))]
sups.block.df <- data.frame('block' = sups.names,
                            'stage' = rep(stages, length(temps)),
                            'treatment' = rep(temps, each = length(stages)))


m.sups <- post.df[,sups.names] %>%
  mutate(index = row_number()) %>%
  melt(id.vars = 'index', value.name = 's_upsilon',
       variable.name = 'block') %>%
  merge(sups.block.df)


ups.df.orig <- post.df[,grep('^upsilon.', names(post.df))]

ups.df <- ups.df.orig

melt.rows <- prod(dim(ups.df))
ups.df$index <- 1:inds

m.ups.df <- ups.df %>% 
  melt(id.vars = 'index', 
       variable.name = 'block', 
       value.name = 'upsilon') %>%
  mutate(stage = rep(rep(rep(stages, each = inds), 
                         length(temps)), length(provs))[1:melt.rows],
         treatment = rep(rep(temps, each = inds*nstage), 
                         length(provs))[1:melt.rows],
         colony = rep(provs,
                      each = inds*nstage*length(temps))[1:melt.rows]) %>%
  subset(select = -block) 


final <- m.pc.df %>%
  merge(m.rho.mult) %>%
  merge(m.rho) %>%
  mutate(rho = rho.mult*rho.orig) %>%
  merge(m.seps.df) %>%
  merge(m.ups.df) %>% 
  merge(m.sups) %>%
  arrange(index, colony, stage) %>%
  mutate(HL = -abs(HL),
         #upsilon.5 = ifelse(is.na(upsilon.5), 0, upsilon.5),
         TA = ta.fun(HL, HH, TL, TH)) 

write.csv(final, 'code/output/all_provs_clean.csv', row.names = FALSE)

gc.post <- final %>%
  select(all_of(c(pars.col, 'rho', 'stage', 'colony'))) %>%
  distinct() %>%
  mutate(index = 1:nrow(.),
         colony = factor(colony, levels = prov.ord.lab.short, 
                         labels = prov.ord.lab)) 

gc.df <- get_curves(gc.post) %>%
  merge(gc.post)

write.csv(gc.df, 'code/output/all_parametric_curves.csv', row.names = FALSE)

gcs.post <- c.post %>%
  mutate(run.index = index,
         colony = factor(colony, levels = prov.ord.lab.short, 
                         labels = prov.ord.lab),
         gc.index = paste0(colony, '_', stage, '_', run.index),
         index = as.numeric(factor(gc.index))) %>%
  arrange(gc.index, treatment)

gcs.post.lst <- gcs.post %>%
  split(., .$gc.index)

gcs.df <- get_curves_pp(gcs.post.lst, spline = TRUE) %>%
  merge(gcs.post)

set.seed(100)
ind.samp <- sample(unique(gcs.df$run.index), size = 100)
gcs.sub <- subset(gcs.df, run.index %in% ind.samp)

write.csv(gcs.sub, 'code/output/all_splined_curves.csv', row.names = FALSE)
