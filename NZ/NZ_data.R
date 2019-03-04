rm(list=ls())

################
## Directories
################

nz_dir <- "C:\\merrill\\stream_networks\\NZ"

data_dir <- file.path(nz_dir, "data")
data_dir2 <- file.path("C:\\merrill\\StreamUtils\\data")

fig_dir <- file.path(nz_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

# devtools::install_github("james-thorson/VAST", ref="development")
# library(VAST)
# # devtools::install_github("merrillrudd/StreamUtils")
# library(StreamUtils)

# library(TMB)
library(tidyverse)
# library(RColorBrewer)
library(proj4)
library(RuddR)

# ################
## Load data
################

## all network
load(file.path(data_dir, "REC2.4fromGDB.Rdata"))
network_raw <- REC2.4fromGDB

# names1 <- colnames(network_raw)
# write.csv(names1, file.path(data_dir, "Network_REC2.csv"))

covariates <- read.csv(file.path(data_dir, "longfin_covariates_network.csv"), header=TRUE, stringsAsFactors=FALSE)
covar_toUse <- covariates[which(covariates$toUse==1),"x"]

network <- network_raw %>%
	select(c('CatName','nzsegment','fnode','tnode','Shape_Leng', 'upcoordX', 'upcoordY', 'downcoordX', 'downcoordY','NextDownID','Headwater','REC2_TerminalSegment', covar_toUse)) %>%
	rename('parent_s' = tnode, 'child_s' = fnode, 'dist_s'=Shape_Leng, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY,'NextDownSeg'=NextDownID, 'Headwater'=Headwater) %>%
	mutate('dist_s' = dist_s / 1000)
network <- network[-which(is.na(network$child_s)),]

hab <- network %>% 
		dplyr::select('nzsegment', 'parent_s','child_s', covar_toUse) %>%
		tidyr::gather(key = covariate, value = value, covar_toUse[1]:covar_toUse[length(covar_toUse)])


## all observations
# load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
# obs_raw <- NZFFD.REC2.Diad.EF
obs_raw <- read.csv(file.path(data_dir, "NZFFD_Joined_REC2_fulldataset.csv"))


# p <- ggplot(network %>% filter(grepl("Waitaki",CatName))) +
# 	geom_point(aes(x = northing_child, y = easting_child, color = Dist2Coast_FromBottom)) +
# 	mytheme()

# p <- ggplot(network %>% filter(grepl("Waitaki",CatName))) +
# 	geom_point(aes(x = northing_child, y = easting_child, color = Q50Cumecs)) +
# 	mytheme()

# p <- ggplot(network %>% filter(grepl("Waitaki",CatName))) +
# 	geom_point(aes(x = northing_child, y = easting_child, color = Q5_normCumecs)) +
# 	mytheme()


# cont <- hab %>% filter(covariate == "Contingency") %>% group_by(child_s) %>% summarise("med_value"=median(value))
# cont2 <- inner_join(network, cont)

# p <- ggplot(cont2 %>% filter(grepl("Waitaki",CatName))) +
# 	geom_point(aes(x = northing_child, y = easting_child, color = med_value)) +
# 	mytheme()

## additional data types
dens_raw <- read.csv(file.path(data_dir, "longfin_density_data.csv"))
length_raw <- read.csv(file.path(data_dir, "longfin_length_data.csv"))
age_raw <- read.csv(file.path(data_dir, "Waitaki_aging_data_DONOTPUBLISH.csv"))

# names2 <- colnames(obs_raw)
# write.csv(names2, file.path(data_dir, "Observations_REC2.csv"))

obs <- obs_raw %>% 
	select(c('catchname', 'catch','nzsegment', 'fishmeth','angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y')) %>%
	rename('catchment'=catchname, 'catch_number'=catch, 'fishmethod'=fishmeth, 'present'=angdie, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY, 'year'=y) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit() %>%
	# mutate('fishmethod' = 'ef') %>%
	mutate('data_type'='encounter') %>%
	rename('data_value'='present') %>%
	mutate('pass'=0) %>%
	mutate('source'='NZFFD')

# length_data <- length_raw %>% 
# 	select('year','fishmeth','pass','nzsegment','locality','length1','length2','length3','length4','length5','length6','length7','length8') %>%
# 	tidyr::gather(key = 'lnum', value = 'length', length1:length8) %>%
# 	rename('fishmethod'='fishmeth','catchment'='locality') %>%
# 	select(-lnum)
# length_obs <- length_data[-which(is.na(length_data$length)),] %>%
# 	mutate('source'='length_data_only')

# age_data <- age_raw %>%
# 	rename('catchment'=colnames(age_raw)[1], 'year'='Year', "length"="Length", "age"="Age") %>%
# 	select('catchment','year','length','age','nzsegment','source')

# al <- full_join(length_obs, age_data)
# agelength <- al[-which(is.na(al$nzsegment)),]

# lengths_all <- agelength[which(is.na(agelength$length)==FALSE),] %>%
# 	select(-age) %>%
# 	mutate('data_type'='length') %>%
# 	rename('data_value'=length)
# ages_all <- agelength[which(is.na(agelength$age)==FALSE),] %>%
# 	select(-length) %>%
# 	mutate('data_type'='age') %>%
# 	rename('data_value'=age)

# al2 <- rbind.data.frame(lengths_all, ages_all)
# al3 <- inner_join(network, al2) %>%
# 	select(-catchment) %>%
# 	rename("catchment"=CatName) %>%
# 	select(names(obs))

# obs_all <- full_join(obs, al3)


# ## doesn't have year
# obs_dens <- dens_raw %>%
# 	select('nzsegment','fish.m2',"AREA",'east','north') %>%
# 	rename('density'='fish.m2', 'Area_m'=AREA, 'easting_child'=east, 'northing_child'=north)

#############################
## latitude and longitude
#############################
calc_NZ_latlon <- function(northing, easting){
	proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
	p <- project(matrix(c(northing, easting),nrow=1), proj=proj4string, inv=T)
	colnames(p) <- c('long', 'lat')
	return(p)
}

## network
network_ll_parent <- lapply(1:nrow(network), function(x){
	p <- calc_NZ_latlon(northing = network$northing_parent[x], easting = network$easting_parent[x])
	return(p)
})
network_ll_parent <- do.call(rbind, network_ll_parent)
network_ll_parent <- data.frame(network_ll_parent) %>% dplyr::rename('long_parent'=long, 'lat_parent'=lat)
network_ll_child <- lapply(1:nrow(network), function(x){
	p <- calc_NZ_latlon(northing = network$northing_child[x], easting = network$easting_child[x])
	return(p)
})
network_ll_child <- do.call(rbind, network_ll_child)
network_ll_child <- data.frame(network_ll_child) %>% dplyr::rename('long_child'=long, 'lat_child'=lat)

network <- cbind.data.frame(network, network_ll_parent, network_ll_child)
network <- network %>% filter(lat_child > -50) 

nzmap <- ggplot(network) +
		geom_point(aes(x = long_child, y = lat_child)) +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
# ggsave(file.path(fig_dir, "NZmap.png"), nzmap)

## observations
obs_ll_parent <- lapply(1:nrow(obs), function(x){
	p <- calc_NZ_latlon(northing = obs$northing_parent[x], easting = obs$easting_parent[x])
	return(p)
})
obs_ll_parent <- do.call(rbind, obs_ll_parent)
obs_ll_parent <- data.frame(obs_ll_parent) %>% dplyr::rename('long_parent'=long, 'lat_parent'=lat)
obs_ll_child <- lapply(1:nrow(obs), function(x){
	p <- calc_NZ_latlon(northing = obs$northing_child[x], easting = obs$easting_child[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
obs_ll_child <- data.frame(obs_ll_child) %>% dplyr::rename('long_child'=long, 'lat_child'=lat)

obs <- cbind.data.frame(obs, obs_ll_parent, obs_ll_child)
obs <- obs %>% filter(nzsegment %in% network$nzsegment == TRUE)

obsmap <- ggplot() +
		geom_point(data=network, aes(x = long_child, y = lat_child), col = "black") +
		geom_point(data=obs %>% filter(data_type=="encounter"), aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
# ggsave(file.path(fig_dir, "NZmap_obs.png"), obsmap)

#############################
## format all data
#############################

## network
network_reformat <- network %>% select('CatName', 'nzsegment','parent_s', 'child_s', 'dist_s', 'lat_child', 'long_child', "NextDownSeg")

root_nodes <- which(network_reformat$parent_s %in% network_reformat$child_s == FALSE)
true_root_node <- which(network_reformat$NextDownSeg==-1)

root_toUse <- lapply(1:length(root_nodes), function(x){
	sub <- network_reformat[root_nodes[x],]
	df <- sub %>% mutate('child_s'=parent_s) %>% mutate('parent_s'=0) %>% mutate('dist_s'=Inf)
	return(df)
})
root_toUse <- do.call(rbind, root_toUse)

network_all <- rbind.data.frame(network_reformat, unique(root_toUse)) %>% 
	rename('lat'=lat_child, 'long'=long_child)

# ## rename nodes
# nodes <- unique(c(network2$child_s, network2$parent_s))
# inodes <- seq_along(nodes)

# net_parents <- sapply(1:nrow(network2), function(x){
#   if(network2$parent_s[x] != 0) new_node <- inodes[which(nodes == network2$parent_s[x])]
#   if(network2$parent_s[x] == 0) new_node <- 0
#   return(new_node)
# })
# net_children <- sapply(1:nrow(network2), function(x) inodes[which(nodes == network2$child_s[x])])

# network_toUse <- network2
# network_toUse$parent_s <- net_parents
# network_toUse$child_s <- net_children

# obs_parents <- sapply(1:nrow(obs), function(x){
#   if(obs$parent_i[x] != 0) new_node <- inodes[which(nodes == obs$parent_i[x])]
#   if(obs$parent_i[x] == 0) new_node <- 0
#   return(new_node)  
# })
# obs_children <- sapply(1:nrow(obs), function(x) inodes[which(nodes == obs$child_i[x])])

# obs_toUse <- obs
# obs_toUse$parent_i <- obs_parents
# obs_toUse$child_i <- obs_children


## observations
obs_reformat <- inner_join(network_all, obs) %>% filter(parent_s !=0) %>% select(-c('catchment','northing_child', 'northing_parent', 'easting_child', 'easting_parent'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})

obs_all <- obs_reformat %>% select(-c('lat_parent', 'long_parent', 'lat_child','long_child', 'NextDownSeg')) %>%
			rename('parent_i' = parent_s, 'child_i' = child_s, 'dist_i' = dist_s)

saveRDS(obs_all, file.path(data_dir, "NZ_observations.rds"))
saveRDS(network_all, file.path(data_dir, "NZ_network.rds"))

obs_all <- readRDS(file.path(data_dir, "NZ_observations.rds"))
network_all <- readRDS(file.path(data_dir, "NZ_network.rds"))
#############################
## subset Waitaki catchment
#############################

network_sub <- network_all %>% filter(grepl("aitaki", CatName))
obs_sub <- obs_all %>% filter(grepl("711",catch_number))
obs_sub <- obs_sub %>% filter(nzsegment %in% network_sub$nzsegment == TRUE)
hab_sub <- hab %>% filter(child_s %in% network_sub$child_s == TRUE)

catchmap <- ggplot() +
		geom_point(data=network_sub, aes(x = long, y = lat), col="black") +
		geom_point(data=obs_sub, aes(x = long, y = lat), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
# ggsave(file.path(fig_dir, "Waitaki_map.png"), catchmap)

#############################
## format
#############################

## rename nodes
nodes <- unique(c(network_sub$child_s, network_sub$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network_sub), function(x){
  if(network_sub$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub$parent_s[x])]
  if(network_sub$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network_sub), function(x) inodes[which(nodes == network_sub$child_s[x])])

network_sub$parent_s <- net_parents
network_sub$child_s <- net_children

obs_parents <- sapply(1:nrow(obs_sub), function(x){
  if(obs_sub$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub$parent_i[x])]
  if(obs_sub$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs_sub), function(x) inodes[which(nodes == obs_sub$child_i[x])])

obs_sub$parent_i <- obs_parents
obs_sub$child_i <- obs_children

saveRDS(obs_sub, file.path(data_dir, "Waitaki_observations.rds"))
saveRDS(network_sub, file.path(data_dir, "Waitaki_network.rds"))

obs_sub <- readRDS(file.path(data_dir, "Waitaki_observations.rds"))
network_sub <- readRDS(file.path(data_dir, "Waitaki_network.rds"))

#########################################
## Waitaki catchment downstream segments
#########################################

obs_child <- unique(obs_sub$child_i)

net_obs <- network_sub %>% filter(child_s %in% obs_child)
nextdown <- network_sub %>% filter(child_s %in% net_obs$parent_s)
save <- rbind.data.frame(net_obs,nextdown)
for(i in 1:100){
  nextdown <- network_sub %>% filter(child_s %in% nextdown$parent_s)
  save <- unique(rbind.data.frame(save, nextdown))
  print(nrow(save))
}
network_sub2 <- save

hab_sub2 <- hab_sub %>% filter(child_s %in% network_sub2$child_s)

## rename nodes
nodes <- unique(c(network_sub2$child_s, network_sub2$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network_sub2), function(x){
  if(network_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub2$parent_s[x])]
  if(network_sub2$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network_sub2), function(x) inodes[which(nodes == network_sub2$child_s[x])])

network_sub2$parent_s <- net_parents
network_sub2$child_s <- net_children

obs_parents <- sapply(1:nrow(obs_sub), function(x){
  if(obs_sub$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub$parent_i[x])]
  if(obs_sub$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs_sub), function(x) inodes[which(nodes == obs_sub$child_i[x])])

obs_sub2 <- obs_sub
obs_sub2$parent_i <- obs_parents
obs_sub2$child_i <- obs_children

saveRDS(obs_sub2, file.path(data_dir, "Waitaki_observations_downstreamOnly.rds"))
saveRDS(network_sub2, file.path(data_dir, "Waitaki_network_downstreamOnly.rds"))


obsfull <- readRDS(file.path(data_dir, "NZ_observations.rds"))
obssub <- readRDS(file.path(data_dir, "Waitaki_observations.rds"))
obssub2 <- readRDS(file.path(data_dir, "Waitaki_observations_downstreamOnly.rds"))

netfull <- readRDS(file.path(data_dir, "NZ_network.rds"))
netsub <- readRDS(file.path(data_dir, "Waitaki_network.rds"))
netsub2 <- readRDS(file.path(data_dir, "Waitaki_network_downstreamOnly.rds"))



## save rda
nz_longfin_eel <- list()
nz_longfin_eel$observations <- obsfull
nz_longfin_eel$network <- netfull
save(nz_longfin_eel, file=file.path(data_dir2, "nz_longfin_eel.rda"))

nz_waitaki_longfin_eel <- list()
nz_waitaki_longfin_eel$observations <- obssub
nz_waitaki_longfin_eel$network <- netsub
save(nz_waitaki_longfin_eel, file=file.path(data_dir2, "nz_waitaki_longfin_eel.rda"))

nz_waitaki_longfin_eel_downstream <- list()
nz_waitaki_longfin_eel_downstream$observations <- obssub2
nz_waitaki_longfin_eel_downstream$network <- netsub2
save(nz_waitaki_longfin_eel_downstream, file=file.path(data_dir2, "nz_waitaki_longfin_eel_downstream.rda"))
