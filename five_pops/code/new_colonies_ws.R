library(tidyverse)
#library(lubridate)
library(reshape2)
library(parallel)

source('code/functions.R')
stages <- paste0('L', 2:6)

weather <- read.csv('data/real_weather.csv')
weather$DATE <- as.Date(weather$DATE)
weather <- subset(weather, Month >= 3 & Month <= 9)
weather$Temp.orig <- weather$Temp
weather$Temp <- round(weather$Temp, 1)
weather$Temp <- pmin(pmax(weather$Temp, 0), 40)

# provs <- unique(weather$Location)
# for (p in provs) {
#   p.sub <- subset(weather, Location == p)
#   p.sub$date.orig <- p.sub$DATE - years(as.numeric(factor(p.sub$Year)))
#   ag.p.sub <- aggregate(data = p.sub, Temp ~ date.orig + Year, mean)
#   write.csv(ag.p.sub, 
#             paste0('data/', p, 'real_weather_full.csv'), row.names = FALSE)
# }

# on.weather <- subset(weather, Location == 'ON' & Year == 2019)
# write.csv(on.weather, 'data/on19_real_weather.csv', row.names = FALSE)

gcs.df <- read.csv('code/output/all_splined_curves.csv') %>%
  select(temp, rate, stage, colony, run.index) %>%
  distinct() %>%
  mutate(dev = rate/24) %>%
  arrange(run.index)

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

age.fun <- function(year, curves = gcs.df, weather.df = on.weather) {
  yr <- substr(as.character(year), 3, 4)
  wthr <- subset(weather.df, Year == year)
  Prov <- unique(wthr$Location)
  
  gc.cast <- dcast(curves, temp + colony + run.index ~ stage,
                   value.var = 'dev')
  names(gc.cast)[names(gc.cast) == 'temp'] <- 'Temp'
  mg <- merge(wthr, gc.cast)
  
  mg$full.index <- paste0(mg$colony, '_', mg$run.index)
  dt.lst <- lapply(unique(mg$full.index), function(ind) {
    sub <- subset(mg, full.index == ind) %>%
      arrange(JDay, Hour)
    # sub <- sub[order(sub$Hour),]
    # sub <- sub[order(sub$JDay),]
    dt <- dev.traverse(sub)
    sub$age <- dt
    return(sub)
  })
  age.df <- bind_rows(dt.lst)
  ag.age.df <- aggregate(data = age.df, age ~ DATE + colony + 
                           run.index + full.index, max)
  write.csv(ag.age.df, paste0('code/output/weathersims/', 
            Prov, yr, '_sim_age.csv'), row.names = FALSE)
  
  p.lst <- lapply(unique(age.df$full.index), function(ind) {
    sub <- subset(age.df, full.index == ind)
    w <- min(which(sub$age == 5))
    pupation <- sub$JDay[w] + sub$Hour[w]/24
    ind.split <- strsplit(ind, split = '_')[[1]]
    pup.df <- data.frame('colony' = ind.split[1], 
                         'index' = ind.split[2],
                         'pup.jd' = pupation)
    return(pup.df)
  })
  p.df <- bind_rows(p.lst)
  
  write.csv(p.df, paste0('code/output/weathersims/', 
                         Prov, yr, '_pupal_dates.csv'), row.names = FALSE)
  return(age.df)
}

#on19 <- age.fun('ON', 2019)
on.weather <- subset(weather, Location == 'ON')
in.weather <- subset(weather, Location == 'IN')
nb.weather <- subset(weather, Location == 'NB')
ab.weather <- subset(weather, Location == 'AB')
qc.weather <- subset(weather, Location == 'QC')

cl <- makeCluster(30)
clusterExport(cl, c('gcs.df', 'qc.weather', 'dev.traverse', 
                    'age.fun', 'stages'))
clusterEvalQ(cl, {
  library(tidyverse)
  library(reshape2)
})
qc.dev <- parLapply(cl, 1991:2020, function(year) {
  af <- age.fun(year, weather.df = qc.weather)
  return(af)
})
stopCluster(cl)

qc.pup.lst <- lapply(1991:2020, function(year) {
  yr.short <- substr(as.character(year), 3, 4)
  file.name <- paste0('QC', yr.short, '_pupal_dates.csv')
  data <- read.csv(file.path('code/output/weathersims', file.name))
  data$year <- year
  middle <- mean(data$pup.jd)
  data$pup.jd.adj <- data$pup.jd - mean(data$pup.jd)
  return(data)
})
qc.pup.df <- bind_rows(qc.pup.lst)
qc.pup.ag <- qc.pup.df %>%
  group_by(colony, year) %>%
  summarize(pup.jd.adj = median(pup.jd.adj)) %>%
  mutate(colony = factor(colony, levels = prov.ord.lab))

ggplot(data = qc.pup.ag) +
  geom_violin(aes(x = colony, y = pup.jd.adj, fill = colony, col = colony),
              scale = 'width', alpha = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = cbP) +
  scale_color_manual(values = cbP) +
  labs(x = '', y = 'Days from Annual Mean Pupation Date', fill = 'Colony', 
       col = 'Colony', title = 'Alberta Weather Simulations') 


on.pup.df$Location <- 'ON'
nb.pup.df$Location <- 'NB'
in.pup.df$Location <- 'IN'
ab.pup.df$Location <- 'AB'
qc.pup.df$Location <- 'QC'
pup.df <- bind_rows(on.pup.df, nb.pup.df, in.pup.df, 
                    ab.pup.df, qc.pup.df) %>%
  mutate(siteyear = paste0(Location, '_', year))
pup.cast <- dcast(data = pup.df, colony + index ~ Location + year, 
                  value.var = 'pup.jd.adj') %>%
  select(-index) %>%
  mutate(colony = factor(colony, levels = prov.ord.lab))
write.csv(pup.cast, 'code/output/all_rel_pupdates.csv', row.names = FALSE)


