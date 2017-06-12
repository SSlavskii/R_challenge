load("spec_dt.RData")

library(dplyr)
start <- Sys.time()
result <- species.dt %>% group_by(species) %>% mutate(rg = rank(value))  #%>% arrange(species, rg)
Sys.time() - start

samp <- sample(1:1000000, 50)
samp <- species.dt[samp, ]
result <- samp %>% group_by(species) %>% mutate(rg = rank(value))

start <- Sys.time()
l <- lapply(unique(species.dt$species), function(x) {
  temp = subset(species.dt, species == x)
  temp$rg = rank(temp$value)
  return(temp[order(temp$rg), ] )
})
result <- do.call(rbind, l)
Sys.time() - start