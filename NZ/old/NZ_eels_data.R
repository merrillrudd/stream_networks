rm(list=ls())

#######################
## Directories
#######################

main_dir <- "C:\\merrill\\stream_networks"
data_dir <- file.path(main_dir, "data")
fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#######################
## Packages
#######################

devtools::install_github("james-thorson/VAST", ref="development")
# devtools::install_github("merrillrudd/VAST")
library(VAST)

library(TMB)
library(dplyr)
library(RColorBrewer)
library(proj4)
# library(maps)
# library(mapproj)



#######################
## Load data
#######################

# ## Network
# load(file.path(data_dir, "FINAL_REC2_FOR_PREDICTIONS.Rdata"))
# rec2 <- REC2
# rm(REC2)

# ## Observations
# load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
# obs_raw <- NZFFD.REC2.Diad.EF


#######################
## Some explorations
#######################

# unique nodes
children <- unique(rec2$nz_fnode) ## upstream
parents <- unique(rec2$nz_tnode)  ## downstream
nodes <- unique(c(parents, children)) ## all nodes

## check proportion of parent nodes that are not child nodes (downstream nodes with no nodes further downstream)
length(which(parents %in% children == TRUE))/length(parents)
## proportion of child nodes that are not parent nodes (upstream nodes that do not have nodes further upstream)
length(which(children %in% parents == TRUE))/length(children)

## subset headings of interest from full REC2 dataset
rec2sub <- rec2 %>% 
	select(nz_tnode, nz_fnode, Shape_length, upcoordX, upcoordY, downcoordX, downcoordY, nzsegment, NextDownSeg, Headwater) %>%
	rename('parent_s'=nz_tnode, 'child_s'=nz_fnode, 'dist_s'=Shape_length, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY, 'segment'=nzsegment, "NextDownSeg"=NextDownSeg, "Headwater"=Headwater)

## find child nodes and rename columns
csub <- rec2sub %>% 
	select(-c('northing_parent', 'easting_parent')) %>%
	rename('northing'=northing_child, 'easting'=easting_child)

## find parent nodes that are not child nodes and set up data frame
psub <- rec2sub %>% filter(parent_s %in% child_s == FALSE)
psub2 <- data.frame('parent_s'=0, 'child_s'=psub$parent_s, 'dist_s'=Inf, 'northing'=psub$northing_parent, 'easting'=psub$easting_parent, 'segment'=psub$segment, "NextDownSeg"=psub$NextDownSeg, "Headwater"=psub$Headwater)

## with segments
Network_sz_ne_seg <- rbind.data.frame(csub, psub2)

### convert UTM to lat/long using PBSmapping package
calc_NZ_latlon <- function(northing, easting){
	proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
	p <- project(matrix(c(northing, easting),nrow=1), proj=proj4string, inv=T)
	colnames(p) <- c('long', 'lat')
	return(p)
}
ll <- lapply(1:nrow(Network_sz_ne_seg), function(x){
	p <- calc_NZ_latlon(northing = Network_sz_ne_seg$northing[x], easting = Network_sz_ne_seg$easting[x])
	return(p)
})
ll <- do.call(rbind, ll)

Network_all <- cbind.data.frame(Network_sz_ne_seg, ll)

saveRDS(Network_all, file.path(data_dir, "Network_segments_all.rds"))
Network_sz_all <- readRDS(file.path(data_dir, "Network_segments_all.rds"))

###########################
## Attach to observations
###########################

obs <- obs_raw %>% 
	select(c('locality','catchname', 'nzsegment', 'angdie', 'y',"Headwater")) %>%
	rename('catchment'=catchname, 'segment'=nzsegment, 'present'=angdie, 'year'=y,"Headwater"=Headwater) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit()


Network_iz_ll <- inner_join(obs, Network_sz_all) %>%
	rename('parent_i'=parent_s, 'child_i'=child_s, 'dist_i'=dist_s)
Network_iz_ll$present <- sapply(1:nrow(Network_iz_ll), function(x) ifelse(Network_iz_ll$present[x]==FALSE, 0, 1))


saveRDS(Network_iz_ll, file.path(data_dir, "Network_observations_all.rds"))
Network_iz_ll <- readRDS(file.path(data_dir, "Network_observations_all.rds"))

###############################
## plot observations
##############################
plot(x = Network_iz_ll$long, y = Network_iz_ll$lat)

###############################
## Subset Waitaki catchment
##############################
##filter observations from Waitaki catchment
Network_iz_sub <- Network_iz_ll %>% filter(grepl("Waitaki",catchment))

## unique observed nodes
Network_iz_sub2 <- unique(Network_iz_sub %>% select(c('parent_i','child_i','dist_i','lat','long','NextDownSeg','Headwater','segment')))

## nodes in specified catchment
wnodes <- unique(c(Network_iz_sub2$parent_i, Network_iz_sub2$child_i))

## subset total network to Waitaki catchment
## find observed nodes in total network
Network_sz_sub <- Network_sz_all[c(which(Network_sz_all$parent_s %in% wnodes),which(Network_sz_all$child_s %in% wnodes)),]
Network_sz_sub2 <- unique(Network_sz_sub %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long', 'NextDownSeg', 'Headwater', 'segment')))

###find where observations are within 5 segments of observations
keep <-  lapply(1:nrow(Network_iz_sub2), function(x){
	sub <- Network_iz_sub2[x,] %>% rename('parent_s'=parent_i, 'child_s'=child_i, 'dist_s'=dist_i)

	for(i in 1:10){
		up <- Network_sz_all %>% filter(NextDownSeg %in% sub$segment) %>% select(-c('northing','easting'))
		sub <- rbind.data.frame(sub, up)
	}

	for(i in 1:10){
		down <- Network_sz_all %>% filter(segment %in% sub$NextDownSeg) %>% select(-c('northing','easting'))
		sub <- rbind.data.frame(sub, down)
	}

	if(any(sub$segment %in% Network_iz_ll$segment)){
		return(sub)
	} else{ return(NULL) }
})
keep2 <- do.call(rbind,keep)
keep3 <- unique(keep2) %>% filter(lat >= -45.1)
lower_waitaki <- keep3

png(file.path(fig_dir, "LowerWaitaki_10connections.png"))
plot(x = lower_waitaki$long , y = lower_waitaki$lat, col="#AAAAAA50", pch=19)
points(x = Network_iz_sub$long, y = Network_iz_sub$lat, col="red", pch=19)
legend("topright", legend=c("Connected segments", 'Observed segments'), pch=19, col=c("#AAAAAA50", "red"))
dev.off()

upstream <- NULL
for(i in 1:400){ ## same number of nodes added with 400, 500, etc.
	if(i==1) add <- Network_sz_all %>% filter(NextDownSeg %in% Network_sz_sub2$segment)
	if(i>1) add <- Network_sz_all %>% filter(NextDownSeg %in% add$segment)
	# points(x = add$long, y = add$lat, pch=19, col="yellow", cex=0.5)
	upstream <- rbind.data.frame(upstream, add)
}
upstream <- unique(upstream %>% select(-c("northing", "easting")))

# find all downstream segments of nodes related to observed nodes
downstream <- NULL
for(i in 1:200){ ## same number of nodes added with 200, 300, etc.
	if(i==1) add <- Network_sz_all %>% filter(segment %in% Network_sz_sub2$NextDownSeg)
	if(i>1) add <- Network_sz_all %>% filter(segment %in% add$NextDownSeg)
	# points(x = add$long,  y = add$lat, pch=19, col="blue", cex=0.5)
	downstream <- rbind.data.frame(downstream,add)
}
downstream <- unique(downstream %>% select(-c("northing","easting")))

png(file.path(fig_dir, "NZ_Waitaki.png"))
plot( x = Network_sz_all$long, y = Network_sz_all$lat, pch=".", xlab="Longitude", ylab="Latitude")#, xlim=c(170.25,171.25), ylim=c(-45.4,-44.3))
points(x = upstream$long, y = upstream$lat, pch=19, col="yellow", cex=0.5)
points(x = downstream$long, y = downstream$lat, pch=19, col="yellow", cex=0.5)
points( x = Network_iz_sub2$long, y = Network_iz_sub2$lat, pch=20, col="red", cex=0.5)
legend("bottomright", legend=c("REC2 segments", "Waitaki catchment", "Observed segments"), col = c("black", "yellow", "red"), pch=19)
dev.off()

waitaki <- unique(rbind.data.frame(Network_sz_sub2, downstream, upstream)) 
waitaki <- waitaki %>% filter(lat > -45.10) %>% filter(parent_s != 0)
Network_iz_sub2 <- Network_iz_sub2 %>% filter(lat > -45.10)
Network_iz_sub <- Network_iz_sub %>% filter(lat > -45.10)

png(file.path(fig_dir, "Waitaki.png"))
plot(x = waitaki$long, y = waitaki$lat, pch=19, col="#AAAAAA50", xlab="Longitude", ylab="Latitude")#, xlim=c(170.2,171.25), ylim=c(-45.2,-44.2))
points(x = Network_iz_sub2$long, y = Network_iz_sub2$lat, col="red", pch=19)
legend("topright", legend=c("Waitaki catchment segments", 'Observed segments'), pch=19, col=c("#AAAAAA50", "red"))
dev.off()

### observations
waitaki_obs <- Network_iz_sub

## all observations are in network
all(waitaki_obs$segment %in% lower_waitaki$segment)

## all observed nodes are in network
all(unique(c(waitaki_obs$parent_i,waitaki_obs$child_i)) %in% unique(c(lower_waitaki$parent_s, lower_waitaki$child_s)))
length(which(unique(c(waitaki_obs$parent_i,waitaki_obs$child_i)) %in% unique(c(lower_waitaki$parent_s, lower_waitaki$child_s)) == FALSE))

saveRDS(waitaki, file.path(data_dir, "waitaki_network.rds"))
saveRDS(lower_waitaki, file.path(data_dir, "lower_waitaki_network.rds"))
saveRDS(waitaki_obs, file.path(data_dir, "waitaki_observations.rds"))

lower_waitaki <- readRDS(file.path(data_dir, "lower_waitaki_network.rds"))
waitaki_obs <- readRDS(file.path(data_dir, "waitaki_observations.rds"))
all(waitaki_obs$segment %in% lower_waitaki$segment)

## renumber network nodes
net_nodes <- unique(c(lower_waitaki$parent_s, lower_waitaki$child_s))
i_net_nodes <- seq_along(net_nodes)


waitaki2 <- lower_waitaki
waitaki2$parent_s <- sapply(1:nrow(lower_waitaki), function(x) i_net_nodes[which(net_nodes ==lower_waitaki$parent_s[x])])
waitaki2$child_s <- sapply(1:nrow(lower_waitaki), function(x) i_net_nodes[which(net_nodes ==lower_waitaki$child_s[x])])
waitaki2 <- waitaki2 %>% select('parent_s', 'child_s', 'dist_s', 'lat', 'long')

find_heads <- which(waitaki2$parent_s %in% waitaki2$child_s ==  FALSE)
add_heads <- data.frame('parent_s' = 0, 'child_s'=waitaki2$parent_s[find_heads], 'dist_s'=Inf, 'lat'=waitaki2$lat[find_heads]+1e-5, 'long'=waitaki2$long[find_heads]+1e-5)
waitaki2 <- rbind.data.frame(waitaki2, add_heads)

waitaki_obs2 <- waitaki_obs
waitaki_obs2$parent_i <- sapply(1:nrow(waitaki_obs), function(x) i_net_nodes[which(net_nodes ==waitaki_obs$parent_i[x])])
waitaki_obs2$child_i <- sapply(1:nrow(waitaki_obs), function(x) i_net_nodes[which(net_nodes ==waitaki_obs$child_i[x])])
waitaki_obs2 <- unique(waitaki_obs2 %>% select('present','year','parent_i', 'child_i', 'dist_i', 'lat', 'long'))

saveRDS(waitaki2, file.path(data_dir, "lower_waitaki_network_renumbered.rds"))
saveRDS(waitaki_obs2, file.path(data_dir, "waitaki_observations_renumbered.rds"))

waitaki_network <- readRDS(file.path(data_dir, 'lower_waitaki_network_renumbered.rds'))
waitaki_observations <- readRDS(file.path(data_dir, 'waitaki_observations_renumbered.rds'))

plot(x = waitaki_network$long, y = waitaki_network$lat, col = "gray")
points(x = waitaki_observations$long, y = waitaki_observations$lat, col="red", pch=19)





###############################
## Subset greater Otago catchment
##############################
## Otago
Network_iz1 <- Network_iz_ll %>% filter(grepl("Waitaki",catchment))
Network_iz2 <- Network_iz_ll %>% filter(grepl("Taieri",catchment))
Network_iz3 <- Network_iz_ll %>% filter(grepl("Clutha",catchment))

Network_iz <- rbind.data.frame(Network_iz1, Network_iz2, Network_iz3) %>% select(-c('northing','easting','locality','catchment'))

wnodes <- unique(c(Network_iz$parent_i, Network_iz$child_i))

Network_sz_sub <- Network_sz_all[c(which(Network_sz_all$parent_s %in% wnodes), which(Network_sz_all$child_s %in% wnodes)),] %>% select(-c('northing','easting'))

## find observed nodes in total network
# upstream <- NULL
# for(i in 1:1000){
# 	if(i==1) add <- Network_sz_all %>% filter(NextDownSeg %in% Network_iz$segment)
# 	if(i>1) add <- Network_sz_all %>% filter(NextDownSeg %in% add$segment)
# 	# points(x = add$long, y = add$lat, pch=19, col="yellow", cex=0.5)
# 	upstream <- rbind.data.frame(upstream, add)
# }
# upstream <- unique(upstream %>% select(-c('northing','easting')))

# downstream <- NULL
# for(i in 1:400){ ## same number of nodes added with 200, 300, etc.
# 	if(i==1) add <- Network_sz_all %>% filter(segment %in% Network_iz$NextDownSeg)
# 	if(i>1) add <- Network_sz_all %>% filter(segment %in% add$NextDownSeg)
# 	# points(x = add$long,  y = add$lat, pch=19, col="blue", cex=0.5)
# 	downstream <- rbind.data.frame(downstream,add)
# }
# downstream <- unique(downstream %>% select(-c("northing","easting")))
# # points(x = Network_iz$long, y = Network_iz$lat, pch=19, col="red")

keep <-  lapply(1:nrow(Network_iz), function(x){
	sub <- Network_iz[x,] %>% rename('parent_s'=parent_i, 'child_s'=child_i, 'dist_s'=dist_i) %>% select(-c('present', 'year'))

	for(i in 1:10){
		up <- Network_sz_all %>% filter(NextDownSeg %in% sub$segment) %>% select(-c('northing','easting'))
		sub <- rbind.data.frame(sub, up)
	}

	for(i in 1:10){
		down <- Network_sz_all %>% filter(segment %in% sub$NextDownSeg) %>% select(-c('northing','easting'))
		sub <- rbind.data.frame(sub, down)
	}

	if(any(sub$segment %in% Network_iz_ll$segment)){
		return(sub)
	} else{ return(NULL) }
})
keep2 <- do.call(rbind,keep)
lower_otago <- unique(keep2)


# Network_sz_sub2 <- unique(rbind.data.frame(Network_sz_sub, upstream, downstream))
all(Network_iz$segment %in% lower_otago$segment)


png(file.path(fig_dir, "NZ_Otago.png"))
plot( x = Network_sz_all$long, y = Network_sz_all$lat, pch=".", xlab="Longitude", ylab="Latitude")#, xlim=c(170.25,171.25), ylim=c(-45.4,-44.3))
points(x = lower_otago$long, y = lower_otago$lat, pch=19, col="yellow", cex=0.5)
points( x = Network_iz$long, y = Network_iz$lat, pch=20, col="red", cex=0.5)
legend("bottomright", legend=c("REC2 segments", "Otago catchment", "Observed segments"), col = c("black", "yellow", "red"), pch=19)
dev.off()


otago <- Network_sz_sub2
otago_obs <- Network_iz

years <- unique(otago_obs$year)
sapply(1:length(years), function(x){
	sub <- otago_obs %>% filter(year == years[x])
	return(sum(sub$present))
})

# ## renumber network nodes
net_nodes <- unique(c(lower_otago$parent_s, lower_otago$child_s))
i_net_nodes <- seq_along(net_nodes)


otago2 <- lower_otago
otago2$parent_s <- sapply(1:nrow(lower_otago), function(x) i_net_nodes[which(net_nodes ==lower_otago$parent_s[x])])
otago2$child_s <- sapply(1:nrow(lower_otago), function(x) i_net_nodes[which(net_nodes ==lower_otago$child_s[x])])
otago2 <- otago2 %>% select('parent_s', 'child_s', 'dist_s', 'lat', 'long')

find_heads <- which(otago2$parent_s %in% otago2$child_s ==  FALSE)
add_heads <- data.frame('parent_s' = 0, 'child_s'=otago2$parent_s[find_heads], 'dist_s'=Inf, 'lat'=otago2$lat[find_heads]+1e-5, 'long'=otago2$long[find_heads]+1e-5)
otago2 <- rbind.data.frame(otago2, add_heads)

otago_obs2 <- otago_obs
otago_obs2$parent_i <- sapply(1:nrow(otago_obs), function(x) i_net_nodes[which(net_nodes ==otago_obs$parent_i[x])])
otago_obs2$child_i <- sapply(1:nrow(otago_obs), function(x) i_net_nodes[which(net_nodes ==otago_obs$child_i[x])])
otago_obs2 <- unique(otago_obs2 %>% select('present','year','parent_i', 'child_i', 'dist_i', 'lat', 'long'))

saveRDS(otago2, file.path(data_dir, "lower_otago_network_renumbered.rds"))
saveRDS(otago_obs2, file.path(data_dir, "otago_observations_renumbered.rds"))



