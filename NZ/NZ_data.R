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

## list of covariates to consider initially at least
covariates <- read.csv(file.path(data_dir, "longfin_covariates_network.csv"), header=TRUE, stringsAsFactors=FALSE)
covar_toUse <- covariates[which(covariates$toUse==1),"x"]

## filter important things from network and adjust distance to km
network <- network_raw %>%
	select(c('CatName','nzsegment','fnode','tnode','Shape_Leng', 'upcoordX', 'upcoordY', 'downcoordX', 'downcoordY','NextDownID','Headwater','REC2_TerminalSegment', covar_toUse)) %>%
	rename('parent_s' = tnode, 'child_s' = fnode, 'dist_s'=Shape_Leng, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX,'NextDownSeg'=NextDownID, 'Headwater'=Headwater) %>%
	mutate('dist_s' = dist_s / 1000)

## remove NAs from network
network <- network[-which(is.na(network$child_s)),]

## network
network_reformat <- network %>% 
	select('CatName', 'nzsegment','parent_s', 'child_s', 'dist_s', 'easting_child','northing_child', "NextDownSeg", covar_toUse)  %>%
	rename('easting'=easting_child, 'northing'=northing_child) %>%
	filter(easting != 0)

root_nodes <- which(network_reformat$parent_s %in% network_reformat$child_s == FALSE)
true_root_node <- which(network_reformat$NextDownSeg==-1)

root_list <- lapply(1:length(root_nodes), function(x){
	sub <- network_reformat[root_nodes[x],]
	df <- sub %>% mutate('child_s'=parent_s) %>% 
		mutate('parent_s'=0) %>% 
		mutate('dist_s'=Inf) 
	return(df)
})
roots <- do.call(rbind, root_list)

child_roots <- unique(roots$child_s)
root_byChild <- lapply(1:length(child_roots), function(x){
	sub <- roots %>% filter(child_s == child_roots[x])
	return(sub)
})
ii <- sapply(1:length(root_byChild), function(x) nrow(root_byChild[[x]]))
root_single <- root_byChild[which(ii==1)]
root_multi <- root_byChild[which(ii > 1)]
multi_to_single <- lapply(1:length(root_multi), function(x){
	sub <- root_multi[[x]]
	if(all(sub$parent_s == sub$parent_s[1]) & all(sub$child_s == sub$child_s[1])){
		out <- sub[1,]
	} else {
		out <- NULL
	}
	return(out)
})
any(is.null(multi_to_single))
root_single2 <- do.call(rbind, root_single)
root_single3 <- do.call(rbind, multi_to_single)
root_toUse <- unique(rbind.data.frame(root_single2, root_single3))

network_all <- rbind.data.frame(network_reformat, unique(root_toUse))
nrow(network_all)
length(unique(network_all$child_s))
nrow(unique(network_all %>% select(easting,northing)))


easting_unique <- unique(network_all$easting)
xx <- table(network_all$easting)
easting_single <- names(xx)[which(xx == 1)]
easting_2 <- names(xx)[which(xx == 2)]
easting_more <- names(xx)[which(xx > 2)] 
net_byEastSingle <- network_all %>% filter(easting %in% easting_single)

net_byEast2 <- network_all %>% filter(easting %in% easting_2) %>%
				mutate(easting_new = ifelse(parent_s == 0, easting + 1, easting)) %>%
				# mutate(northing_new = ifelse(parent_s == 0, northing + 1, northing)) %>%
				select(-c(easting)) %>%
				rename('easting'=easting_new)

network_adj <- rbind.data.frame(net_byEastSingle, net_byEast2)
nrow(network_adj)
length(unique(network_adj$child_s))
nrow(unique(network_adj %>% select(easting,northing)))

xx2 <- table(network_adj$easting)
yy2 <- table(network_adj$northing)
easting_22 <- names(xx2)[which(xx2 > 1)]
northing_22 <- names(yy2)[which(yy2 > 1)]

net2 <- network_adj %>% filter(easting %in% easting_22) %>% filter(northing %in% northing_22)
net3 <- network_adj %>% filter(easting != unique(net2$easting)) %>% filter(northing != unique(net2$northing))
net2$easting <- c(net2$easting[1], net2$easting[2]+1)
net2$northing <- c(net2$northing[1], net2$northing[2]+1)

network_adj2 <- rbind.data.frame(net3, net2)

nrow(network_adj2)
length(unique(network_adj2$child_s))
nrow(unique(network_adj2 %>% select(easting,northing)))

# child_count <- table(network_adj2$child_s)
# child_multi <- child_count[which(child_count > 1)]
# find <- network_adj2 %>% filter(child_s == names(child_multi))
# network_adj3 <- rbind.data.frame(network_adj2 %>% filter(child_s != names(child_multi)), find[1,])

# nrow(network_adj3)
# length(unique(network_adj3$child_s))
# nrow(unique(network_adj3 %>% select(easting,northing)))


sub1 <- network_all %>% filter(easting == easting_more[1])
net4 <- sub1 %>% 
	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1, easting[4])) %>%
	mutate(northing = c(northing[1:3],northing[4]+1))

sub2 <- network_all %>% filter(easting == easting_more[2])
net5 <- sub2 %>% 
	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1)) %>%
	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1))

sub3 <- network_all %>% filter(easting == easting_more[3])
net6 <- sub3 %>% 
	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1, easting[4])) %>%
	mutate(northing = c(northing[1:3], northing[4]-1))

sub4 <- network_all %>% filter(easting == easting_more[4])
net7 <- sub4 %>% 
	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1)) %>%
	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1))

network_adj3 <- rbind.data.frame(network_adj2, net4, net5, net6, net7)

nrow(network_adj3)
length(unique(network_adj3$child_s))
nrow(unique(network_adj3 %>% select(easting,northing)))


# network_adj4 <- rbind.data.frame(network_adj3 %>% filter(child_s != names(child_multi)), find[1,])

# nrow(network_adj4)
# length(unique(network_adj4$child_s))
# nrow(unique(network_adj4 %>% select(easting,northing)))

#############################
## latitude and longitude
#############################

## function to calculate latitude and longitude from eastings and northings
calc_NZ_latlon <- function(northing, easting){
	proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
	p <- project(matrix(c(easting, northing),nrow=1), proj=proj4string, inv=T)
	colnames(p) <- c('long', 'lat')
	return(p)
}

## latitude and longitude for parent nodes in network
# network_ll_parent <- lapply(1:nrow(network_adj3), function(x){
# 	p <- calc_NZ_latlon(northing = network_adj3$northing_parent[x], easting = network_adj3$easting_parent[x])
# 	return(p)
# })
# network_ll_parent <- do.call(rbind, network_ll_parent)
# network_ll_parent <- data.frame(network_ll_parent) %>% dplyr::rename('long_parent'=long, 'lat_parent'=lat)

## latitude and longitude for child nodes in network
network_ll_child <- lapply(1:nrow(network_adj3), function(x){
	p <- calc_NZ_latlon(northing = network_adj3$northing[x], easting = network_adj3$easting[x])
	return(p)
})
network_ll_child <- do.call(rbind, network_ll_child)
# network_ll_child <- data.frame(network_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

## attach latitude and longtiude to network
network_full <- cbind.data.frame(network_adj3, network_ll_child)

## map of full NZ network
nzmap <- ggplot(network_full) +
		geom_point(aes(x = easting, y = northing), cex=0.2) +
		xlab("Easting") + ylab("Northing") +
		mytheme()
# ggsave(file.path(fig_dir, "NZmap.png"), nzmap)

## all observations
# load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
# obs_raw <- NZFFD.REC2.Diad.EF
obs_raw <- read.csv(file.path(data_dir, "NZFFD_Joined_REC2_fulldataset.csv"))

## additional data types
dens_raw <- read.csv(file.path(data_dir, "longfin_density_data.csv"))
length_raw <- read.csv(file.path(data_dir, "longfin_length_data.csv"))
age_raw <- read.csv(file.path(data_dir, "Waitaki_aging_data_DONOTPUBLISH.csv"))

## filter information we need from observations, do some renaming, and label encounter data
obs <- obs_raw %>% 
	select(c('catchname', 'catch','nzsegment', 'fishmeth','angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y')) %>%
	rename('catchment'=catchname, 'catch_number'=catch, 'fishmethod'=fishmeth, 'present'=angdie, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX, 'year'=y) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit() %>%
	# mutate('fishmethod' = 'ef') %>%
	mutate('data_type'='encounter') %>%
	rename('data_value'='present') %>%
	mutate('pass'=0) %>%
	mutate('source'='NZFFD') %>%
	select(-c(easting_parent, northing_parent)) %>%
	rename('easting'=easting_child, 'northing'=northing_child)

obs_ll_child <- lapply(1:nrow(obs), function(x){
	p <- calc_NZ_latlon(northing = obs$northing[x], easting = obs$easting[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
# obs_ll_child <- data.frame(obs_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

obs <- cbind.data.frame(obs, obs_ll_child)
all(obs$nzsegment %in% network_full$nzsegment)
# obs <- obs %>% filter(nzsegment %in% network_full$nzsegment == TRUE)

## observations
network_sz <- network_full %>% select('nzsegment','parent_s','child_s','dist_s')
obs_reformat <- inner_join(network_sz, obs, by='nzsegment') %>% filter(parent_s !=0) # %>% select(-c('catchment','northing', 'easting'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})

obs_full <- obs_reformat %>% 
			rename('parent_i' = parent_s, 'child_i' = child_s, 'dist_i' = dist_s)

obsmap <- ggplot() +
		geom_point(data=network_full, aes(x = easting, y = northing), col = "black", cex=0.2) +
		geom_point(data=obs %>% filter(data_type=="encounter"), aes(x = easting, y = northing), col = "red") +
		xlab("Easting") + ylab("Northing") +
		mytheme()
ggsave(file.path(fig_dir, "NZmap_obs.png"), obsmap)

## select habitat data from network separately
hab_full <- network_full %>% 
		dplyr::select('nzsegment', 'parent_s','child_s', covar_toUse) %>%
		tidyr::gather(key = covariate, value = value, covar_toUse[1]:covar_toUse[length(covar_toUse)])


#### code to potentially add length, age, and later density data to all observations
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
## format all data
#############################
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



saveRDS(obs_full, file.path(data_dir, "NZ_observations.rds"))
saveRDS(network_full, file.path(data_dir, "NZ_network.rds"))
saveRDS(hab_full, file.path(data_dir, "NZ_habitat.rds"))

obs_full <- readRDS(file.path(data_dir, "NZ_observations.rds"))
network_full <- readRDS(file.path(data_dir, "NZ_network.rds"))
hab_full <- readRDS(file.path(data_dir, "NZ_habitat.rds"))
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
