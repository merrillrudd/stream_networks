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
# library(proj4)
# library(RuddR)

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
		dplyr::select('nzsegment', covar_toUse) %>%
		tidyr::gather(key = covariate, value = value, covar_toUse[1]:covar_toUse[length(covar_toUse)]) %>%


## all observations
load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
obs_raw <- NZFFD.REC2.Diad.EF

## covariates for longfins
cov_list <- read.csv(file.path(data_dir, "longfin_covariates.csv"), stringsAsFactors=FALSE)
covs <- cov_list[,2]


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
	select(c('catchname', 'nzsegment', 'angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y')) %>%
	rename('catchment'=catchname, 'present'=angdie, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY, 'year'=y) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit() %>%
	mutate('fishmethod' = 'ef') %>%
	mutate('data_type'='encounter') %>%
	rename('data_value'='present') %>%
	mutate('pass'=0) %>%
	mutate('source'='encounter')

length_data <- length_raw %>% 
	select('year','fishmeth','pass','nzsegment','locality','length1','length2','length3','length4','length5','length6','length7','length8') %>%
	tidyr::gather(key = 'lnum', value = 'length', length1:length8) %>%
	rename('fishmethod'='fishmeth','catchment'='locality') %>%
	select(-lnum)
length_obs <- length_data[-which(is.na(length_data$length)),] %>%
	mutate('source'='length_data_only')

age_data <- age_raw %>%
	rename('catchment'=colnames(age_raw)[1], 'year'='Year', "length"="Length", "age"="Age") %>%
	select('catchment','year','length','age','nzsegment','source')

al <- full_join(length_obs, age_data)
agelength <- al[-which(is.na(al$nzsegment)),]

lengths_all <- agelength[which(is.na(agelength$length)==FALSE),] %>%
	select(-age) %>%
	mutate('data_type'='length') %>%
	rename('data_value'=length)
ages_all <- agelength[which(is.na(agelength$age)==FALSE),] %>%
	select(-length) %>%
	mutate('data_type'='age') %>%
	rename('data_value'=age)

al2 <- rbind.data.frame(lengths_all, ages_all)
al3 <- inner_join(network, al2) %>%
	select(-catchment) %>%
	rename("catchment"=CatName) %>%
	select(names(obs))

obs_all <- full_join(obs, al3)


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
ggsave(file.path(fig_dir, "NZmap.png"), nzmap)

## observations
obs_ll_parent <- lapply(1:nrow(obs_all), function(x){
	p <- calc_NZ_latlon(northing = obs_all$northing_parent[x], easting = obs_all$easting_parent[x])
	return(p)
})
obs_ll_parent <- do.call(rbind, obs_ll_parent)
obs_ll_parent <- data.frame(obs_ll_parent) %>% dplyr::rename('long_parent'=long, 'lat_parent'=lat)
obs_ll_child <- lapply(1:nrow(obs_all), function(x){
	p <- calc_NZ_latlon(northing = obs_all$northing_child[x], easting = obs_all$easting_child[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
obs_ll_child <- data.frame(obs_ll_child) %>% dplyr::rename('long_child'=long, 'lat_child'=lat)

obs_all <- cbind.data.frame(obs_all, obs_ll_parent, obs_ll_child)
obs_all <- obs_all %>% filter(nzsegment %in% network$nzsegment == TRUE)

obsmap <- ggplot() +
		geom_point(data=network, aes(x = long_child, y = lat_child), col = "black") +
		geom_point(data=obs_all %>% filter(data_type=="encounter"), aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
# ggsave(file.path(fig_dir, "NZmap_obs.png"), obsmap)

#############################
## format all data
#############################

## network
network_reformat <- network %>% select('parent_s', 'child_s', 'dist_s', 'lat_child', 'long_child', 'lat_parent', 'long_parent', "NextDownSeg")

root_nodes <- which(network_reformat$parent_s %in% network_reformat$child_s == FALSE)
true_root_node <- which(network_reformat$NextDownSeg==-1)

root_toUse <- lapply(1:length(root_nodes), function(x){
	sub <- network_reformat[root_nodes[x],]
	df <- data.frame('parent_s' = 0, 'child_s' = sub$parent_s, 'dist_s' = Inf, "lat" = sub$lat_parent, 'long' = sub$long_parent)
	return(df)
})
root_toUse <- do.call(rbind, root_toUse)

network_toUse <- network_reformat %>% select(-c('lat_parent','long_parent','NextDownSeg')) %>%
	rename('lat'=lat_child, 'long'=long_child)
network_toUse <- rbind.data.frame(network_toUse, unique(root_toUse))

## observations
obs_reformat <- inner_join(obs_all, network) %>% select(-c('catchment','nzsegment', 'northing_child', 'northing_parent', 'easting_child', 'easting_parent', 'Headwater', 'CatName', 'REC2_TerminalSegment'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})

obs_toUse <- obs_reformat %>% select(-c('lat_parent', 'long_parent', 'NextDownSeg'))%>%
			rename('lat'=lat_child, 'long'=long_child) %>%
			rename('parent_i' = parent_s, 'child_i' = child_s, 'dist_i' = dist_s)

obs_list <- NULL
obs_list$observations <- obs_toUse
obs_list$habitat <- hab

saveRDS(obs_list, file.path(data_dir, "NZ_observations.rds"))
saveRDS(network_toUse, file.path(data_dir, "NZ_network.rds"))

#############################
## subset Waitaki catchment
#############################

network_sub <- network %>% filter(grepl("Waitaki", CatName))
obs_sub <- obs_all %>% filter(grepl("Waitaki", catchment))
obs_sub <- obs_sub %>% filter(nzsegment %in% network_sub$nzsegment == TRUE)
hab_sub <- hab %>% filter(child_s %in% network_sub$child_s == TRUE)

submap <- obsmap +
		geom_point(data=network_sub, aes(x = long_child, y = lat_child), col="blue") +
		geom_point(data=obs_sub, aes(x = long_child, y = lat_child), col = "yellow")

catchmap <- ggplot() +
		geom_point(data=network_sub, aes(x = long_child, y = lat_child), col="black") +
		geom_point(data=obs_sub, aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
# ggsave(file.path(fig_dir, "Waitaki_map.png"), catchmap)

# network_sub_cut <- network_sub %>% 
# 			filter(long_parent >= min(c(obs_sub$long_child,obs_sub$long_parent))) %>%
# 			filter(lat_parent <= max(c(obs_sub$lat_child,obs_sub$lat_parent)))
# all(obs_sub$nzsegment %in% network_sub_cut$nzsegment)
network_sub_cut <- network_sub

#############################
## format
#############################

## rename nodes
nodes <- unique(c(network_sub_cut$parent_s, network_sub_cut$child_s))
# nodes <- unique(c(network_sub$parent_s, network_sub$child_s))
inodes <- seq_along(nodes)

## network
network_reformat <- network_sub_cut %>% select('parent_s', 'child_s', 'dist_s', 'lat_child', 'long_child', 'lat_parent', 'long_parent', "NextDownSeg")
network_parents <- sapply(1:nrow(network_reformat), function(x) inodes[which(nodes == network_reformat$parent_s[x])])
network_children <- sapply(1:nrow(network_reformat), function(x) inodes[which(nodes == network_reformat$child_s[x])])
network_reformat$parent_s <- network_parents
network_reformat$child_s <- network_children

root_nodes <- which(network_reformat$parent_s %in% network_reformat$child_s == FALSE)
true_root_node <- which(network_reformat$NextDownSeg==-1)

root_toUse <- lapply(1:length(root_nodes), function(x){
	sub <- network_reformat[root_nodes[x],]
	df <- data.frame('parent_s' = 0, 'child_s' = sub$parent_s, 'dist_s' = Inf, "lat" = sub$lat_parent, 'long' = sub$long_parent)
	return(df)
})
root_toUse <- do.call(rbind, root_toUse)

network_toUse <- network_reformat %>% select(-c('lat_parent','long_parent','NextDownSeg')) %>%
	rename('lat'=lat_child, 'long'=long_child)
network_toUse <- rbind.data.frame(network_toUse, unique(root_toUse))

## observations
obs_reformat <- inner_join(obs_sub, network_sub_cut) %>% select(-c('catchment','nzsegment', 'northing_child', 'northing_parent', 'easting_child', 'easting_parent', 'Headwater', 'CatName', 'REC2_TerminalSegment'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})

obs_parents <- sapply(1:nrow(obs_reformat), function(x) inodes[which(nodes == obs_reformat$parent_s[x])])
obs_children <- sapply(1:nrow(obs_reformat), function(x) inodes[which(nodes == obs_reformat$child_s[x])])
obs_reformat$parent_s <- obs_parents
obs_reformat$child_s <- obs_children

obs_toUse <- obs_reformat %>% select(-c('lat_parent', 'long_parent', 'NextDownSeg'))%>%
			rename('lat'=lat_child, 'long'=long_child) %>%
			rename('parent_i' = parent_s, 'child_i' = child_s, 'dist_i' = dist_s)

hab_children <- sapply(1:nrow(hab_sub), function(x) inodes[which(nodes == hab_sub$child_s[x])])
hab_sub$child_s <- hab_children

obs_list <- NULL
obs_list$observations <- obs_toUse
obs_list$habitat <- hab_sub

saveRDS(obs_list, file.path(data_dir, "Waitaki_observations.rds"))
saveRDS(network_toUse, file.path(data_dir, "Waitaki_network.rds"))

nz1 <- readRDS(file.path(data_dir, "NZ_observations.rds"))
nz2 <- readRDS(file.path(data_dir, "Waitaki_observations.rds"))

nz3 <- readRDS(file.path(data_dir, "NZ_network.rds"))
nz4 <- readRDS(file.path(data_dir, "Waitaki_network.rds"))

## save rda
nz_longfin_eel <- list()
nz_longfin_eel$observations <- nz1$observations %>% filter(data_type!="age")
nz_longfin_eel$habitat <- nz1$habitat
nz_longfin_eel$network <- nz3
save(nz_longfin_eel, file=file.path(data_dir2, "nz_longfin_eel.rda"))

nz_waitaki_longfin_eel <- list()
nz_waitaki_longfin_eel$observations <- nz2$observations %>% filter(data_type!="age")
nz_waitaki_longfin_eel$habitat <- nz2$habitat
nz_waitaki_longfin_eel$network <- nz4
save(nz_waitaki_longfin_eel, file=file.path(data_dir2, "nz_waitaki_longfin_eel.rda"))

### figure - observations by year
obs <- nz2 %>% select('lat','long','present','year') %>% mutate(observation = ifelse(present == 0, 'absent','present'))
net <- nz4 %>% select('lat','long')

years <- unique(obs$year)[order(unique(obs$year))]
net_byYear <- lapply(1:length(years), function(x){
	netyr <- data.frame(net, "year"=years[x])
	return(netyr)
})
net_byYear <- do.call(rbind, net_byYear)

p <- ggplot() +
	geom_point(data = net_byYear, aes(x = long, y = lat)) +
	geom_point(data = obs, aes(x = long, y = lat, color = observation), cex=2, alpha=0.8) +
	scale_color_brewer(palette = "Set1") +
	facet_wrap(~year) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear.png"), width=15, height=10)