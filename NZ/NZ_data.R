rm(list=ls())

################
## Directories
################

nz_dir <- "C:\\merrill\\stream_networks\\NZ"

data_dir <- file.path(nz_dir, "data")

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
# library(sf)

# ################
## Load data
################

## all network
load(file.path(data_dir, "REC2.4fromGDB.Rdata"))
network_raw <- REC2.4fromGDB

load(file.path(data_dir, "FINAL_REC2_FOR_PREDICTIONS.Rdata"))
network_old <- REC2

choose <- network_old %>% select("nzsegment", "FWENZ_DSDamAffected", "REC1_WidthHUC.MALF_cumecs") %>%
	rename('DamAffected'=FWENZ_DSDamAffected, "width"=REC1_WidthHUC.MALF_cumecs) %>%
	mutate('width' = width / 1000)

## list of covariates to consider initially at least
covariates <- read.csv(file.path(data_dir, "longfin_covariates_network.csv"), header=TRUE, stringsAsFactors=FALSE)
covar_toUse <- covariates[which(covariates$toUse==1),"x"]

## filter important things from network and adjust distance to km
network <- network_raw %>%
	select(c('CatName','nzsegment','fnode','tnode','Shape_Leng', 'upcoordX', 'upcoordY', 'downcoordX', 'downcoordY','NextDownID','Headwater','REC2_TerminalSegment', covar_toUse)) %>%
	rename('parent_s' = tnode, 'child_s' = fnode, 'length'=Shape_Leng, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX,'NextDownSeg'=NextDownID, 'Headwater'=Headwater) %>%
	mutate('length' = length / 1000)

## remove NAs from network
network <- network[-which(is.na(network$child_s)),]

network_combine <- left_join(network, choose)
sapply(1:ncol(network_combine), function(x) length(which(is.na(network_combine[,x])))/nrow(network_combine))


## wait to adjust widths for smaller areas
# ## order by northing
# width <- akima::interp(x = network_combine$easting_child[-which(is.na(network_combine$width))], y = network_combine$northing_child[-which(is.na(network_combine$width))], z = network_combine$width[-which(is.na(network_combine$width))], xo = network_combine$easting_child[which(is.na(network_combine$width))], yo = network_combine$northing_child[which(is.na(network_combine$width))], duplicate = "median")

# # sub <- network_combine %>% filter(grepl('aitaki',CatName))

covar_toUse <- c(covar_toUse, "DamAffected")

## network
network_reformat <- network_combine %>% 
	select('CatName', 'nzsegment','parent_s', 'child_s', 'length', 'width','easting_child','northing_child', "NextDownSeg", covar_toUse)  %>%
	rename('easting'=easting_child, 'northing'=northing_child) %>%
	filter(easting != 0) %>%
	mutate('dist_s' = length * width)

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

e_mult <- names(table(network_all$easting))[which(table(network_all$easting)>1)]
e_uni <- network_all %>% filter(easting %in% e_mult == FALSE)
set.seed(123)
e_rep <- network_all %>% filter(easting %in% e_mult) %>% mutate(easting = easting + runif(length(easting),-1,1)) 

network_all2 <- rbind.data.frame(e_uni, e_rep)
max(table(network_all2$easting))

n_mult <- names(table(network_all2$northing))[which(table(network_all2$northing)>1)]
n_uni <- network_all2 %>% filter(northing %in% n_mult == FALSE)
set.seed(456)
n_rep <- network_all2 %>% filter(northing %in% n_mult) %>% mutate(northing = northing + runif(length(northing),-1,1))

network_all3 <- rbind.data.frame(n_uni, n_rep)
max(table(network_all3$northing))

# easting_unique <- unique(network_all$easting)
# xx <- table(network_all$easting)
# easting_single <- names(xx)[which(xx == 1)]
# easting_2 <- names(xx)[which(xx == 2)]
# easting_more <- names(xx)[which(xx > 2)] 
# net_byEastSingle <- network_all %>% filter(easting %in% easting_single)

# net_byEast2 <- network_all %>% filter(easting %in% easting_2) %>%
# 				mutate(easting_new = ifelse(parent_s == 0, easting + 1, easting)) %>%
# 				# mutate(northing_new = ifelse(parent_s == 0, northing + 1, northing)) %>%
# 				select(-c(easting)) %>%
# 				rename('easting'=easting_new)

# network_adj <- rbind.data.frame(net_byEastSingle, net_byEast2)
# nrow(network_adj)
# length(unique(network_adj$child_s))
# nrow(unique(network_adj %>% select(easting,northing)))

# xx2 <- table(network_adj$easting)
# yy2 <- table(network_adj$northing)
# easting_22 <- names(xx2)[which(xx2 > 1)]
# northing_22 <- names(yy2)[which(yy2 > 1)]

# net2 <- network_adj %>% filter(easting %in% easting_22) %>% filter(northing %in% northing_22)
# net3 <- network_adj %>% filter(easting != unique(net2$easting)) %>% filter(northing != unique(net2$northing))
# net2$easting <- c(net2$easting[1], net2$easting[2]+1)
# net2$northing <- c(net2$northing[1], net2$northing[2]+1)

# network_adj2 <- rbind.data.frame(net3, net2)

# nrow(network_adj2)
# length(unique(network_adj2$child_s))
# nrow(unique(network_adj2 %>% select(easting,northing)))

# # child_count <- table(network_adj2$child_s)
# # child_multi <- child_count[which(child_count > 1)]
# # find <- network_adj2 %>% filter(child_s == names(child_multi))
# # network_all3 <- rbind.data.frame(network_adj2 %>% filter(child_s != names(child_multi)), find[1,])

# # nrow(network_all3)
# # length(unique(network_all3$child_s))
# # nrow(unique(network_all3 %>% select(easting,northing)))


# sub1 <- network_all %>% filter(easting == easting_more[1])
# net4 <- sub1 %>% 
# 	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1, easting[4]-1)) %>%
# 	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1, northing[4]+1))

# sub2 <- network_all %>% filter(easting == easting_more[2])
# net5 <- sub2 %>% 
# 	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1)) %>%
# 	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1))

# sub3 <- network_all %>% filter(easting == easting_more[3])
# net6 <- sub3 %>% 
# 	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1, easting[4]+1)) %>%
# 	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1, northing[4]-1))

# sub4 <- network_all %>% filter(easting == easting_more[4])
# net7 <- sub4 %>% 
# 	mutate(easting = c(easting[1], easting[2]-1, easting[3]+1)) %>%
# 	mutate(northing = c(northing[1], northing[2]-1, northing[3]+1))

# network_all3 <- rbind.data.frame(network_adj2, net4, net5, net6, net7)

# nrow(network_all3)
# length(unique(network_all3$child_s))
# nrow(unique(network_all3 %>% select(easting,northing)))


# dam_gaps <- network_all3 %>% filter(is.na(DamAffected))
# dam_gap_child <- unique(dam_gaps$child_s)

# dam_gap_down <- network_all3 %>% filter(child_s %in% dam_gaps$parent_s) %>% filter(is.na(DamAffected)==FALSE)
# dam_gap_up <- network_all3 %>% filter(parent_s %in% dam_gaps$child_s) %>% filter(is.na(DamAffected)==FALSE)

# dam_new <- lapply(1:nrow(dam_gaps), function(x){
# 	sub <- dam_gaps[x,]
# 	nextdown <- dam_gap_down




# 	up_info <- dam_gap_up %>% filter(parent_s == sub$child_s[x])
# 	down_info <- dam_gap_down %>% filter(child_s == sub$parent_s[x])
# 	if(nrow(up_info)==0 & nrow(down_info)==0) return(data.frame("child_s"=dam_gap_child[x], "DamAffected"=NA))
# 	if(nrow(up_info)!=0) return(data.frame('child_s'=dam_gap_child[x], 'DamAffected'=up_info$DamAffected))
# 	if(nrow(down_info)!=0) return(data.frame('child_s'=dam_gap_child[x], "DamAffected"=down_info$DamAffected))
# })
# dam_check <- do.call(rbind, dam_new)


# save <- rbind.data.frame(net_obs,nextdown)
# for(i in 1:100){
#   nextdown <- network_sub %>% filter(child_s %in% nextdown$parent_s)
#   save <- unique(rbind.data.frame(save, nextdown))
#   print(nrow(save))
# }
# network_sub2 <- save



# lapply(1:length(dam_gap_child), function(x){
# 	sub <- dam_gaps %>% filter(child_s == dam_gap_child[x])
# 	dam_up <- network_all3 %>% filter(parent_s == sub$child_s)
# 	dam_down <- network_all3 %>% filter(child_s == sub$parent_s)
# 	up_info <- ifelse(nrow(dam_up)==0, NA, dam_up$DamAffected)
# 	down_info <- ifelse(nrow(dam_down)==0, NA, dam_down$DamAffected)
# 	if(is.na(up_info) & is.na(down_info)){
# 		if(nrow(dam_down)>0){
# 			dam_down <- network_all3 %>% filter(child_s == dam_down$parent_s)
# 			down_info <- ifelse(nrow(dam_down)==0, NA, dam_down$DamAffected)
# 		}
# 		if(nrow(dam_up)>0){

# 		}
# 	}
# })

# dam_next <- network_all3 %>% filter(nzsegment %in% dam_gaps$NextDownSeg)
# dam_up <- network_all3 %>% filter(NextDownSeg %in% dam_gaps$nzsegment)




# network_adj4 <- rbind.data.frame(network_all3 %>% filter(child_s != names(child_multi)), find[1,])

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
# network_ll_parent <- lapply(1:nrow(network_all3), function(x){
# 	p <- calc_NZ_latlon(northing = network_all3$northing_parent[x], easting = network_all3$easting_parent[x])
# 	return(p)
# })
# network_ll_parent <- do.call(rbind, network_ll_parent)
# network_ll_parent <- data.frame(network_ll_parent) %>% dplyr::rename('long_parent'=long, 'lat_parent'=lat)

## latitude and longitude for child nodes in network
network_ll_child <- lapply(1:nrow(network_all3), function(x){
	p <- calc_NZ_latlon(northing = network_all3$northing[x], easting = network_all3$easting[x])
	return(p)
})
network_ll_child <- do.call(rbind, network_ll_child)
# network_ll_child <- data.frame(network_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

## attach latitude and longtiude to network
network_full <- cbind.data.frame(network_all3, network_ll_child)
nrow(network_full)
nrow(unique(network_full))
nrow(network_full %>% select('lat','long'))
nrow(network_full %>% select('easting','northing'))

## map of full NZ network

## all observations
# load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
# obs_raw <- NZFFD.REC2.Diad.EF
obs_raw <- read.csv(file.path(data_dir, "NZFFD_Joined_REC2_fulldataset.csv"))

## additional data types
dens_raw <- read.csv(file.path(data_dir, "longfin_density_data.csv"))
length_raw <- read.csv(file.path(data_dir, "longfin_length_data.csv"))
age_raw <- read.csv(file.path(data_dir, "Waitaki_aging_data_DONOTPUBLISH.csv"))

waikato_dens_raw <- read.csv(file.path(data_dir, "Waikato Abundance.REC.csv"), stringsAsFactors=FALSE)

## filter information we need from observations, do some renaming, and label encounter data
obs_enc <- obs_raw %>% 
	select(c('catchname','nzsegment', 'fishmeth','angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y', 'org.groups')) %>%
	rename('catchment'=catchname,'fishmethod'=fishmeth, 'present'=angdie, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX, 'year'=y, 'agency'=org.groups) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit() %>%
	# mutate('fishmethod' = 'ef') %>%
	mutate('data_type'='encounter') %>%
	rename('data_value'='present') %>%
	mutate('pass'=0) %>%
	mutate('source'='NZFFD') %>%
	select(-c(easting_parent, northing_parent)) %>%
	rename('easting'=easting_child, 'northing'=northing_child)

obs_ll_child <- lapply(1:nrow(obs_enc), function(x){
	p <- calc_NZ_latlon(northing = obs_enc$northing[x], easting = obs_enc$easting[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
# obs_ll_child <- data.frame(obs_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

obs_enc <- cbind.data.frame(obs_enc, obs_ll_child)
all(obs_enc$nzsegment %in% network_full$nzsegment)
# obs <- obs %>% filter(nzsegment %in% network_full$nzsegment == TRUE)

## observations
network_sz <- network_full %>% select('CatName','nzsegment','parent_s','child_s','width')
obs_reformat <- inner_join(network_sz, obs_enc, by='nzsegment') %>% filter(parent_s !=0) # %>% select(-c('catchment','northing', 'easting'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})
obs_reformat <- obs_reformat %>%
	mutate('length' = 125 / 1000) %>%
	mutate('dist_i' = length * width) %>%
	select(-catchment)

## waikato densities 
waikato_dens <- waikato_dens_raw %>%
	select("Sample.Date_fish",'nzsegment', "E.NZTM", "N.NZTM", "average_measured_stream_channel_width_m", "Angdie.ALL") 
waikato_dens$Sample.Date_fish <- as.character(waikato_dens$Sample.Date_fish)
waikato_dens$year <- sapply(1:nrow(waikato_dens), function(x) strsplit(waikato_dens$Sample.Date_fish[x], "/")[[1]][3])
waikato_dens <- waikato_dens %>%
	select(-Sample.Date_fish) %>%
	rename('easting'=E.NZTM, "northing"=N.NZTM, 'width' = average_measured_stream_channel_width_m, 'count'=Angdie.ALL) %>%
	mutate(width = width / 1000) %>%
	mutate(length = 150 / 1000) %>%
	mutate(dist_i = width * length) %>%
	mutate('data_type' = "count") %>%
	rename('data_value' = count)

obs_ll_child <- lapply(1:nrow(waikato_dens), function(x){
	p <- calc_NZ_latlon(northing = waikato_dens$northing[x], easting = waikato_dens$easting[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
# obs_ll_child <- data.frame(obs_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

waikato_dens_ll <- cbind.data.frame(waikato_dens, obs_ll_child)
all(waikato_dens_ll$nzsegment %in% network_full$nzsegment)

waikato_reformat <- inner_join(network_sz, waikato_dens_ll, by="nzsegment") %>% filter(parent_s != 0)
waikato_reformat <- waikato_reformat %>% 
	select(-width.x) %>%
	rename('width' = width.y) %>%
	mutate('fishmethod'=unique(obs_reformat$fishmethod)[grepl("Electric",unique(obs_reformat$fishmethod))]) %>%
	mutate('agency'=unique(obs_reformat$agency)[grepl('council',unique(obs_reformat$agency))]) %>%
	mutate(pass = 0) %>%
	mutate(source = "Waikato_region")

obs_all <- rbind.data.frame(obs_reformat, waikato_reformat)

obs_full <- obs_all %>% 
			rename('parent_i' = parent_s, 'child_i' = child_s)


## select habitat data from network separately
hab_full <- network_full %>% 
		dplyr::select('nzsegment', 'parent_s','child_s', covar_toUse, easting, northing) %>%
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

nzmap <- ggplot(network_full) +
		geom_point(aes(x = long, y = lat), cex=0.2) +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
ggsave(file.path(fig_dir, "NZmap.png"), nzmap)

obsmap <- ggplot() +
		geom_point(data=network_full, aes(x = long, y = lat), col = "black", cex=0.2) +
		geom_point(data=obs_full %>% filter(data_type=="encounter"), aes(x = long, y = lat, color = data_type)) +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()
ggsave(file.path(fig_dir, "NZmap_obs_encounter.png"), obsmap)

#############################
## subset Waitaki catchment
#############################

network_sub <- network_full %>% filter(grepl("aitaki", CatName))

width_new <- network_sub$width
width_new[which(is.na(network_sub$width))] <- median(network_sub$width, na.rm=TRUE)
network_sub$width <- width_new

# net_sub2 <- network_full %>% filter(grepl("hitney", CatName))
# obs_sub <- obs_full %>% filter(grepl("711",catch_number))
obs_sub <- obs_full %>% filter(grepl("aitaki", CatName))

all(obs_sub$nzsegment %in% network_sub$nzsegment)
# obs_sub <- obs_sub %>% filter(nzsegment %in% network_sub$nzsegment == TRUE)
hab_sub <- hab_full %>% filter(child_s %in% network_sub$child_s == TRUE)


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

hab_parents <- sapply(1:nrow(hab_sub), function(x){
  if(hab_sub$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub$parent_s[x])]
  if(hab_sub$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
hab_children <- sapply(1:nrow(hab_sub), function(x) inodes[which(nodes == hab_sub$child_s[x])])

hab_sub$parent_s <- hab_parents
hab_sub$child_s <- hab_children

sapply(1:length(covar_toUse), function(x){
	sub <- hab_sub %>% filter(covariate == covar_toUse[x])
	any(is.na(sub$value))
})

hab_sub2 <- lapply(1:length(covar_toUse), function(x){
	sub <- hab_sub %>% filter(covariate == covar_toUse[x])
	# any(is.na(sub$value))
	if(any(is.na(sub$value))){

		if(covar_toUse[x]!="DamAffected"){
			interp_east <- sub$easting[which(is.na(sub$value)==FALSE)]
			interp_north <- sub$northing[which(is.na(sub$value)==FALSE)]
			interp_z <- sub$value[which(is.na(sub$value)==FALSE)]	

			find_df <- data.frame('east' = sub$easting[which(is.na(sub$value))], 'north' = sub$northing[which(is.na(sub$value))])	

			east <- sub$easting[order(sub$easting)]
			north <- sub$northing[order(sub$northing)]
			# mat2 <- zoo::na.approx(object = mat)
			compute <- akima::interp(x = interp_east, y = interp_north, z = interp_z, xo=east, yo=north, extrap=TRUE)
			mat2 <- compute$z	

			vals <- sapply(1:nrow(find_df), function(y){
				mat2[which(compute$x == find_df$east[y]), which(compute$y == find_df$north[y])]
			})	

			inp_vals <- sub$value
			inp_vals[which(is.na(inp_vals))] <- vals	

			sub$value <- inp_vals	

			if(length(which(is.na(sub$value)))==1){
				xx <- sub[(which(is.na(sub$value))-5):(which(is.na(sub$value))+5),]
				xx2 <- xx[order(xx$easting),]
				val_inp <- median(xx$value, na.rm=TRUE)
				sub$value[which(is.na(sub$value))] <- val_inp
			}
		}
		if(covar_toUse[x]=="DamAffected"){


			inp <- sub$value
			inp[which(is.na(inp))] <- 2
			ggplot(sub) + geom_point(aes(x = easting, y = northing, color = factor(inp)))

			input_val <- sub$value
			input_val[which(is.na(input_val) & sub$northing > 5025000)] <- 1
			input_val[which(is.na(input_val) & sub$northing < 5025000)] <- 0

			sub$value <- input_val
		}
	}
	return(sub)
})
check <- sapply(1:length(hab_sub2), function(x) any(is.na(hab_sub2[[x]]$value)))
all(check == FALSE)

hab_sub2 <- do.call(rbind, hab_sub2)

saveRDS(obs_sub, file.path(data_dir, "Waitaki_observations.rds"))
saveRDS(network_sub, file.path(data_dir, "Waitaki_network.rds"))
saveRDS(hab_sub2, file.path(data_dir, "Waitaki_habitat.rds"))

obs_sub <- readRDS(file.path(data_dir, "Waitaki_observations.rds"))
network_sub <- readRDS(file.path(data_dir, "Waitaki_network.rds"))
hab_sub <- readRDS(file.path(data_dir, "Waitaki_habitat.rds"))

catchmap <- ggplot() +
		geom_point(data=network_sub, aes(x = long, y = lat), col="gray") +
		geom_point(data=obs_sub, aes(x = long, y = lat, fill = data_type), pch=22, alpha=0.6) +
		scale_fill_brewer(palette = "Set1") +
		xlab("Longitude") + ylab("Latitude") +
		guides(fill = FALSE) +
		mytheme()
ggsave(file.path(fig_dir, "Waitaki_map.png"), catchmap)

catchmap2 <- ggplot() +
		geom_point(data=network_full, aes(x = long, y = lat), col = "black", cex=0.2) +
		geom_point(data = network_sub, aes(x = long, y = lat), col = "gray") +
		# geom_point(data=obs_full %>% filter(data_type=="encounter"), aes(x = long, y = lat, fill=data_type), pch=22, alpha=0.6) +
		xlab("Longitude") + ylab("Latitude") +
		# scale_fill_brewer(palette = "Set1") +
		mytheme()
ggsave(file.path(fig_dir, "Waitaki_on_NZ.png"), catchmap2)

#########################################
## Waitaki catchment downstream segments
#########################################

# loc_df <- hab_sub %>% select('easting', 'northing')
# loc_mat <- as.matrix(loc_df)
# loc <- st_linestring(loc_mat)
# loc_simp <- st_simplify(loc, dTolerance = 5)

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

hab_parents <- sapply(1:nrow(hab_sub2), function(x){
  if(hab_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub2$parent_s[x])]
  if(hab_sub2$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
hab_children <- sapply(1:nrow(hab_sub2), function(x) inodes[which(nodes == hab_sub2$child_s[x])])

hab_sub2$parent_s <- hab_parents
hab_sub2$child_s <- hab_children


saveRDS(obs_sub2, file.path(data_dir, "Waitaki_observations_downstreamOnly.rds"))
saveRDS(network_sub2, file.path(data_dir, "Waitaki_network_downstreamOnly.rds"))
saveRDS(hab_sub2, file.path(data_dir, "Waitaki_habitat_downstreamOnly.rds"))

obs_sub2 <- readRDS(file.path(data_dir, "Waitaki_observations_downstreamOnly.rds"))
network_sub2 <- readRDS(file.path(data_dir, "Waitaki_network_downstreamOnly.rds"))

		l2 <- lapply(1:nrow(network_sub2), function(x){
			parent <- network_sub2$parent_s[x]
			find <- network_sub2 %>% filter(child_s == parent)
			if(nrow(find)>0) out <- cbind.data.frame(network_sub2[x,], 'long2'=find$long, 'lat2'=find$lat)
			if(nrow(find)==0) out <- cbind.data.frame(network_sub2[x,], 'long2'=NA, 'lat2'=NA)
			# if(nrow(find)>0) out <- cbind.data.frame(network_sub2[x,], 'long2'=find$long, 'lat2'=find$lat)
			# if(nrow(find)==0) out <- cbind.data.frame(network_sub2[x,], 'long2'=NA, 'lat2'=NA)
			return(out)
		})
		l2 <- do.call(rbind, l2)

catchmap <- ggplot() +
		geom_point(data=network_sub2, aes(x = long, y = lat), col="gray") +
		geom_segment(data=l2, aes(x = long2,y = lat2, xend = long, yend = lat), col="gray") +
		geom_point(data=obs_sub2, aes(x = long, y = lat, fill = data_type), pch=22, alpha=0.6) +
		scale_fill_brewer(palette = "Set1") +
		xlab("Longitude") + ylab("Latitude") +
		guides(fill = FALSE) +
		mytheme()
ggsave(file.path(fig_dir, "Waitaki_map_downstream.png"), catchmap)


#############################
## subset Waikato catchment
#############################

network_sub_wk <- network_full %>% filter(grepl("aikato", CatName))

# net_sub2 <- network_full %>% filter(grepl("hitney", CatName))
# obs_sub <- obs_full %>% filter(grepl("711",catch_number))
obs_sub_wk <- obs_full %>% filter(grepl("aikato", CatName))

all(obs_sub_wk$nzsegment %in% network_sub_wk$nzsegment)
hab_sub_wk <- hab_full %>% filter(child_s %in% network_sub_wk$child_s == TRUE)

catchmap2 <- ggplot() +
		geom_point(data=network_sub_wk, aes(x = easting, y = northing), col="black") +
		geom_point(data=obs_sub_wk, aes(x = easting, y = northing, color = data_type)) +
		xlab("Easting") + ylab("Northing") +
		mytheme()
ggsave(file.path(fig_dir, "Waikato_map.png"), catchmap)

#############################
## format
#############################

## rename nodes
nodes <- unique(c(network_sub_wk$child_s, network_sub_wk$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network_sub_wk), function(x){
  if(network_sub_wk$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub_wk$parent_s[x])]
  if(network_sub_wk$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network_sub_wk), function(x) inodes[which(nodes == network_sub_wk$child_s[x])])

network_sub_wk$parent_s <- net_parents
network_sub_wk$child_s <- net_children

obs_parents <- sapply(1:nrow(obs_sub_wk), function(x){
  if(obs_sub_wk$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub_wk$parent_i[x])]
  if(obs_sub_wk$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs_sub_wk), function(x) inodes[which(nodes == obs_sub_wk$child_i[x])])

obs_sub_wk$parent_i <- obs_parents
obs_sub_wk$child_i <- obs_children

hab_parents <- sapply(1:nrow(hab_sub_wk), function(x){
  if(hab_sub_wk$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub_wk$parent_s[x])]
  if(hab_sub_wk$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
hab_children <- sapply(1:nrow(hab_sub_wk), function(x) inodes[which(nodes == hab_sub_wk$child_s[x])])

hab_sub_wk$parent_s <- hab_parents
hab_sub_wk$child_s <- hab_children


saveRDS(obs_sub_wk, file.path(data_dir, "Waikato_observations.rds"))
saveRDS(network_sub_wk, file.path(data_dir, "Waikato_network.rds"))
saveRDS(hab_sub_wk, file.path(data_dir, "Waikato_habitat.rds"))

obs_sub_wk <- readRDS(file.path(data_dir, "Waikato_observations.rds"))
network_sub_wk <- readRDS(file.path(data_dir, "Waikato_network.rds"))
hab_sub_wk <- readRDS(file.path(data_dir, "Waikato_habitat.rds"))

#########################################
## Waikato catchment downstream segments
#########################################

obs_sub_wk2 <- obs_sub_wk %>% filter(is.na(dist_i)==FALSE)
obs_child <- unique(obs_sub_wk2$child_i)

net_obs <- network_sub_wk %>% filter(child_s %in% obs_child)
nextdown <- network_sub_wk %>% filter(child_s %in% net_obs$parent_s)
save <- rbind.data.frame(net_obs,nextdown)
for(i in 1:100){
  nextdown <- network_sub_wk %>% filter(child_s %in% nextdown$parent_s)
  save <- unique(rbind.data.frame(save, nextdown))
  print(nrow(save))
}
network_sub_wk2 <- save

network_sub_wk2 <- network_sub_wk2 %>% filter(is.na(dist_s) == FALSE)

hab_sub_wk2 <- hab_sub_wk %>% filter(child_s %in% network_sub_wk2$child_s)

## rename nodes
nodes <- unique(c(network_sub_wk2$child_s, network_sub_wk2$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network_sub_wk2), function(x){
  if(network_sub_wk2$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub_wk2$parent_s[x])]
  if(network_sub_wk2$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network_sub_wk2), function(x) inodes[which(nodes == network_sub_wk2$child_s[x])])

network_sub_wk2$parent_s <- net_parents
network_sub_wk2$child_s <- net_children

obs_parents <- sapply(1:nrow(obs_sub_wk), function(x){
  if(obs_sub_wk$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub_wk$parent_i[x])]
  if(obs_sub_wk$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs_sub_wk), function(x) inodes[which(nodes == obs_sub_wk$child_i[x])])

# obs_sub_wk2 <- obs_sub_wk
obs_sub_wk2$parent_i <- obs_parents
obs_sub_wk2$child_i <- obs_children

all(obs_sub_wk2$child_i %in% network_sub_wk2$child_s)

hab_parents <- sapply(1:nrow(hab_sub_wk2), function(x){
  if(hab_sub_wk2$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub_wk2$parent_s[x])]
  if(hab_sub_wk2$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
hab_children <- sapply(1:nrow(hab_sub_wk2), function(x) inodes[which(nodes == hab_sub_wk2$child_s[x])])

hab_sub_wk2$parent_s <- hab_parents
hab_sub_wk2$child_s <- hab_children



catchmap <- ggplot() +
		geom_point(data=network_sub_wk2, aes(x = easting, y = northing), col="black") +
		geom_point(data=obs_sub_wk2, aes(x = easting, y = northing, color = data_type)) +
		geom_point(data=network_sub_wk2 %>% filter(parent_s == 0), aes(x = easting, y = northing), fill = "goldenrod", pch=21, cex=3) +
		xlab("Easting") + ylab("Northing") +
		mytheme()
ggsave(file.path(fig_dir, "Waikato_map_downstream.png"), catchmap)


saveRDS(obs_sub_wk2, file.path(data_dir, "Waikato_observations_downstreamOnly.rds"))
saveRDS(network_sub_wk2, file.path(data_dir, "Waikato_network_downstreamOnly.rds"))
saveRDS(hab_sub_wk2, file.path(data_dir, "Waikato_habitat_downstreamOnly.rds"))







obsfull <- readRDS(file.path(data_dir, "NZ_observations.rds"))
obssub <- readRDS(file.path(data_dir, "Waitaki_observations.rds"))
obssub2 <- readRDS(file.path(data_dir, "Waitaki_observations_downstreamOnly.rds"))
# obssubwk <- readRDS(file.path(data_dir, "Waikato_observations.rds"))
# obssubwk2 <- readRDS(file.path(data_dir, "Waikato_observations_downstreamOnly.rds"))

netfull <- readRDS(file.path(data_dir, "NZ_network.rds"))
netsub <- readRDS(file.path(data_dir, "Waitaki_network.rds"))
netsub2 <- readRDS(file.path(data_dir, "Waitaki_network_downstreamOnly.rds"))
# netsubwk <- readRDS(file.path(data_dir, "Waikato_network.rds"))
# netsubwk2 <- readRDS(file.path(data_dir, "Waikato_network_downstreamOnly.rds"))


habfull <- readRDS(file.path(data_dir, "NZ_habitat.rds"))
habsub <- readRDS(file.path(data_dir, "Waitaki_habitat.rds"))
habsub2 <- readRDS(file.path(data_dir, "Waitaki_habitat_downstreamOnly.rds"))
# habsubwk <- readRDS(file.path(data_dir, "Waikato_habitat.rds"))
# habsubwk2 <- readRDS(file.path(data_dir, "Waikato_habitat_downstreamOnly.rds"))

data_dir2 <- file.path("C:\\merrill\\stream_networks\\NZ\\data_save")
## save rda
nz_longfin_eel <- list()
nz_longfin_eel$observations <- obsfull
nz_longfin_eel$network <- netfull
save(nz_longfin_eel, file=file.path(data_dir2, "nz_longfin_eel.rda"))

nz_longfin_eel_habitat <- habfull
save(nz_longfin_eel_habitat, file=file.path(data_dir2, 'nz_longfin_eel_habitat.rda'))

nz_waitaki_longfin_eel <- list()
nz_waitaki_longfin_eel$observations <- obssub
nz_waitaki_longfin_eel$network <- netsub
nz_waitaki_longfin_eel$habitat <- habsub
save(nz_waitaki_longfin_eel, file=file.path(data_dir2, "nz_waitaki_longfin_eel.rda"))

nz_waitaki_longfin_eel_downstream <- list()
nz_waitaki_longfin_eel_downstream$observations <- obssub2
nz_waitaki_longfin_eel_downstream$network <- netsub2
nz_waitaki_longfin_eel_downstream$habitat <- habsub2
save(nz_waitaki_longfin_eel_downstream, file=file.path(data_dir2, "nz_waitaki_longfin_eel_downstream.rda"))

nz_waikato_longfin_eel <- list()
nz_waikato_longfin_eel$observations <- obssubwk
nz_waikato_longfin_eel$network <- netsubwk
nz_waikato_longfin_eel$habitat <- habsubwk
save(nz_waikato_longfin_eel, file=file.path(data_dir2, "nz_waikato_longfin_eel.rda"))

nz_waikato_longfin_eel_downstream <- list()
nz_waikato_longfin_eel_downstream$observations <- obssubwk2
nz_waikato_longfin_eel_downstream$network <- netsubwk2
nz_waikato_longfin_eel_downstream$habitat <- habsubwk2
save(nz_waikato_longfin_eel_downstream, file=file.path(data_dir2, "nz_waikato_longfin_eel_downstream.rda"))
