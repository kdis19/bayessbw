library(tidyverse)
library(lubridate)
library(reshape2)
source('code/functions.R')
stages <- paste0('L', 2:6)

weather <- read.csv('data/real_weather.csv')
weather$DATE <- as.Date(weather$DATE)
weather <- subset(weather, Month >= 3 & Month <= 9)
weather$Temp.orig <- weather$Temp
weather$Temp <- round(weather$Temp, 1)
weather$Temp <- pmin(pmax(weather$Temp, 0), 40)

on.weather <- subset(weather, Location == 'ON' & Year == 2019)
write.csv(on.weather, 'data/on_real_weather.csv', row.names = FALSE)

in.weather <- subset(weather, Location == 'IN' & Year == 1997)
write.csv(in.weather, 'data/in_real_weather.csv', row.names = FALSE)

# ag.gc <- read.csv('code/output/ag_devcurves.csv')
# ag.gc <- subset(ag.gc, select = -c(low, high))
# names(ag.gc)[1] <- c('Temp')
# ag.gc$dev <- ag.gc$mid/24

gc.df <- read.csv('code/output/all_posterior_curves.csv')

gc.cast <- dcast(gc.df, temp + province + index ~ stage, value.var = 'dev')
names(gc.cast)[names(gc.cast) == 'temp'] <- 'Temp'
mg <- merge(on.weather, gc.cast)

dev.traverse <- function(df) {
  #dates <- df$DATE
  df <- df[,stages]
  ind <- 0
  st.ind <- 1
  dev.vec <- vector()
  for (s in 1:length(stages)) {
    stage <- stages[s]
    vec <- df[,s]
    cs.vec <- cumsum(vec)
    
    ind.prev <- ind
    ind <- min(which(cs.vec > 1))
    
    cs.vec <- cs.vec + st.ind - 1
    st.ind <- st.ind + 1
    
    dev.vec <- c(dev.vec, cs.vec[(ind.prev + 1):ind])
    df[1:ind, (1:length(stages)) > s] <- 0
  }
  dev.vec <- c(dev.vec, rep(5, pmax(0, nrow(df) - length(dev.vec))))
  #plot(dates, dev.vec, type = 'l')
  return(dev.vec)
}

mg$full.index <- paste0(mg$province, '_', mg$index)
dt.lst <- lapply(unique(mg$full.index), function(ind) {
  sub <- subset(mg, full.index == ind)
  sub <- sub[order(sub$Hour),]
  sub <- sub[order(sub$JDay),]
  dt <- dev.traverse(sub)
  sub$age <- dt
  return(sub)
})
age.df <- bind_rows(dt.lst)
write.csv(age.df, 'code/output/sim_age.csv')
ag.age.df <- aggregate(data = age.df, age ~ DATE + province + index + full.index, max)
write.csv(ag.age.df, 'code/output/sim_age_daily.csv', row.names = FALSE)

##### Compare pupal dates #####
p.lst <- lapply(unique(age.df$full.index), function(ind) {
  sub <- subset(age.df, full.index == ind)
  w <- min(which(sub$age == 5))
  pupation <- sub$JDay[w]
  ind.split <- strsplit(ind, split = '_')[[1]]
  pup.df <- data.frame('province' = ind.split[1], 
                       'index' = ind.split[2],
                       'pup.jd' = pupation)
  return(pup.df)
})
p.df <- bind_rows(p.lst)

write.csv(p.df, 'code/output/pupal_dates.csv', row.names = FALSE)

yr.lst <- lapply(unique(p.df$index), function(ind) {
  sub <- subset(p.df, index == ind)
  mean.jd <- mean(sub$pup.jd)
  sub$rel.pup.jd <- sub$pup.jd/mean.jd
  return(sub)
})
yr.df <- bind_rows(yr.lst)

yr.df$province <- factor(yr.df$province, levels = c('IN', 'AB', 'ON', 'QC', 'NB', 'IPU'))
