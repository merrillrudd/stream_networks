rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks\\OR"
data_dir <- file.path(main_dir, "data")

fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

devtools::install_github("james-thorson/VAST", ref="development")
library(VAST)
devtools::install_github("merrillrudd/FishStatsUtils", ref = "stream")
library(FishStatsUtils)

library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

################
## Load data
################

network <- read.csv(file.path(data_dir, "Siletz", "Network_WGS.csv"))
nodes <- read.csv(file.path(data_dir, "Siletz", "Nodes_WGS.csv"))
hab_raw <- read.csv(file.path(data_dir, "Siletz", "Habitat_Siletz_WGS_line_final.csv"))
juv <- read.csv(file.path(data_dir, "Siletz", "Juvenile_Siletz_WGS_line_final.csv"))
spawn <- read.csv(file.path(data_dir, "Siletz", "Spawning_Siletz_WGS_line_final.csv"))



##############
##### NETWORK
##############
net <- network %>% 
	select("Shape_Leng", "ChildNode", "ParentNode", "lat", "long") %>%
	rename("dist_s"="Shape_Leng", "child_s"="ChildNode", "parent_s"="ParentNode") %>%
	mutate("dist_s" = dist_s / 1000)
nod <- nodes %>% select("NodeID", "lat", "long")

if(length(net$parent_s %in% net$child_s) > 0){
	iroot <- which(net$parent_s %in% net$child_s == FALSE)
	new <- lapply(1:length(iroot), function(x){
		sub <- net[iroot[x],]
		nodsub <- nod %>% filter(NodeID == sub$parent_s)
		df <- data.frame("dist_s"=Inf, "child_s"=sub$parent_s, "parent_s"=0, "lat"=nodsub$lat, "long"=nodsub$long)
		return(df)
	})
	add <- do.call(rbind, new)
	net_toUse <- rbind.data.frame(net, add)
} else {
	net_toUse <- net
}

### check for unique locations
nrow(net_toUse)
length(unique(net_toUse$lat))
length(unique(net_toUse$long))

## change one of the longitudes to be a unique location by adding a tiny amount
net_toUse %>% filter(long == names(table(net_toUse$long))[which(table(net_toUse$long)==max(table(net_toUse$long)))])
index <- which(net_toUse$long == names(table(net_toUse$long))[which(table(net_toUse$long)==max(table(net_toUse$long)))])
net_toUse$long[index[2]] <- net_toUse$long[index[2]] + 1e-4
net_toUse$long[index]

## recheck
nrow(net_toUse)
length(unique(net_toUse$lat))
length(unique(net_toUse$long))


##############
##### SURVEYS
##############

obs_spawn <- spawn %>% 
		select("SpawnYear", "adults_per_km", "lat", "long", "ChildNode") %>%
		rename("year"="SpawnYear", "density"="adults_per_km","child_i"="ChildNode") %>%
		mutate("survey"="spawners") %>%
		mutate('surveynum' = 1)

obs_juv <- juv %>%
		select("JuvYear", "parr_per_k", "lat", "long","ChildNode") %>%
		rename("year"="JuvYear", "density"="parr_per_k","child_i"="ChildNode") %>%
		na.omit() %>%
		mutate("survey"="juveniles") %>% 
		mutate("surveynum" = 2)

obs_dens <- rbind.data.frame(obs_spawn, obs_juv)


##############
##### HABITAT
##############

hab <- hab_raw %>% 
	select("lat","long","YEAR_", "ChildNode","PCTSCCHNLA","ACW","PCTSWPOOL","POOL1P_KM","POOLS100","RIFFLEDEP","LWDVOL1","WGTED_SLOPE_GRAVEL","WGTED_ALLUNITS_BEDROCK") %>%
	rename("child_i"="ChildNode", "year"="YEAR_")
dcovar <- c("PCTSCCHNLA","ACW","PCTSWPOOL","POOL1P_KM","POOLS100","RIFFLEDEP","LWDVOL1","WGTED_SLOPE_GRAVEL","WGTED_ALLUNITS_BEDROCK")

# lapply(1:length(dcovar), function(x){

# 	## subset covariate
# 	sub <- hab %>% select(c('lat','long','year','child_i',dcovar[x]))

# 	## identify latitude and longitude of child node to match network, not habitat observation
# 	lat_new <- unlist(sapply(1:nrow(sub), function(y){
# 		lat_old <- sub[y,] %>% select('lat')
# 		lat_new <- net_toUse %>% filter(child_s == sub[y,'child_i']) %>% select('lat')
# 		return(lat_new)
# 	}))
# 	long_new <- unlist(sapply(1:nrow(sub), function(y){
# 		long_old <- sub[y,] %>% select('long')
# 		long_new <- net_toUse %>% filter(child_s == sub[y,'child_i']) %>% select('long')
# 		return(long_new)
# 	}))
# 	sub$lat2 <- lat_new
# 	sub$long2 <- long_new

# 	lapply(1:length(years), function(y){
# 		sub2 <- sub %>% filter(year==years[y])
# 		sub2 <- sub2[order(sub2$lat2),]

# 		## vectors with information to interpolate between
# 		interp_lat <- sub2$lat2
# 		interp_long <- sub2$long2
# 		interp_z <- sub2[,dcovar[x]]

# 		find_lat <- net_toUse$lat[which(net_toUse$lat %in% interp_lat == FALSE)]
# 		find_lat <- find_lat[order(find_lat)]
# 		find_long <- net_toUse$long[which(net_toUse$long %in% interp_long == FALSE)]
# 		find_long <- find_long[order(find_long)]

# 		compute <- akima::interp(x = interp_long, y = interp_lat, z = interp_z, xo = find_lat, yo = find_long, extrap = TRUE, duplicate = "median", linear=FALSE)
# 		mat <- compute$z

# 		vals <- sapply(1:length(find_lat), function(z){
# 			mat[which(compute$x == find_long[z]), which(compute$y == find_lat[z])]
# 		})

# 		df <- data.frame('lat'=find_lat, 'long'=find_long, 'value'=vals)

	

# 	})






	## network matrix: lat x long x year
	netmat <- array(NA, dim=c(nrow(net_toUse),ncol(net_toUse),length(years)))

	## names of rows and columns = index in network of that lat/long
	rownames(netmat) <- order(net_toUse$lat)
	colnames(netmat) <- order(net_toUse$long)

	## lat and long in order
	lat_order <- net_toUse$lat[order(net_toUse$lat)]
	long_order <- net_toUse$long[order(net_toUse$long)]




	p <- ggplot(sub) +
		geom_point(aes(x = long, y = lat, color = PCTSCCHNLA))

})



obs_hab <- hab %>%
		tidyr::gather(key = "covariate", value = "value", PCTSCCHNLA:WGTED_ALLUNITS_BEDROCK)

####################
##### PLOT NETWORK
####################
years <- unique(obs_dens$year)[order(unique(obs_dens$year))]
net_byYear <- lapply(1:length(years), function(x){
	netyr <- data.frame(net_toUse, "year"=years[x])
	return(netyr)
})
net_byYear <- do.call(rbind, net_byYear)








saveRDS(obs_hab, file.path(data_dir, "Siletz_habitat.rds"))
saveRDS(obs_dens, file.path(data_dir, "Siletz_observations_density.rds"))
saveRDS(net_toUse, file.path(data_dir, "Siletz_network.rds"))
write.csv(net_toUse, file.path(data_dir, "Siletz_network.csv"), row.names=FALSE)


# or_siletz_coho <- list()
# or_siletz_coho$observations <- obs_dens
# or_siletz_coho$network <- net_toUse
# save(or_siletz_coho, file=file.path(data_dir2, "or_siletz_coho.rda"))


#############################
## Figures - data
##############################
### look for latitude and longitude of parent nodes
l2 <- lapply(1:nrow(net_toUse), function(x){
	parent <- net_toUse$parent_s[x]
	find <- net_toUse %>% filter(child_s == parent)
	if(nrow(find)>0) out <- cbind.data.frame(net_toUse[x,], 'lat2'=find$lat, 'long2'=find$lon)
	if(nrow(find)==0) out <- cbind.data.frame(net_toUse[x,], 'lat2'=NA, 'long2'=NA)
	return(out)
})
l2 <- do.call(rbind, l2) %>% mutate('region'="Siletz")

aa <- ggplot(l2) +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
    geom_point(aes(x = long, y = lat), color = "black") +
    geom_segment(aes(x = long2,y = lat2, xend = long, yend = lat), arrow=arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), col="gray", alpha=0.6) +
    xlab("Longitude") + ylab("Latitude") +  
    mytheme()
ggsave(file.path(fig_dir, "Network.png"), aa, width=8, height=8)

bb <- ggplot(l2) +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
    geom_point(aes(x = long, y = lat), color = "black") +
    geom_segment(aes(x = long2,y = lat2, xend = long, yend = lat), arrow=arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), col="gray", alpha=0.6) +
    geom_point(data = obs_dens, aes(x = long, y = lat, fill=survey), pch=22, cex=3) +
    scale_fill_brewer(palette = "Set1")+
    xlab("Longitude") + ylab("Latitude") +  
    mytheme()
ggsave(file.path(fig_dir, "Network_observations.png"), bb, width=8, height=6)

cc <- aa +
	geom_point(data = obs_hab, aes(x = long, y = lat, fill = value), cex=3, pch=22) +
	scale_fill_viridis_c() +
	facet_wrap(.~covariate)

years <- unique(obs_dens$year)[order(unique(obs_dens$year))]
net_byYear <- lapply(1:length(years), function(x){
	netyr <- data.frame(net_toUse, "year"=years[x])
	return(netyr)
})
net_byYear <- do.call(rbind, net_byYear)

mapbyyr <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat)) +
	geom_point(data=obs_dens, aes(x = long, y = lat, fill=survey), cex=3, pch=22) +
	facet_wrap(~year) +
	scale_fill_brewer(palette = "Set1") +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear.png"), mapbyyr, width=10, height=8)

mapbyyrjuv <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat)) +
	geom_point(data=obs_dens %>% filter(survey=="juveniles"), aes(x = long, y = lat, size=density), fill=brewer.pal(3,"Set1")[1], pch=22) +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	facet_wrap(~year) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear_juveniles.png"), mapbyyrjuv, width=10, height=8)

mapbyyrspawn <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat), cex=2) +
	geom_point(data=obs_dens %>% filter(survey=="spawners"), aes(x = long, y = lat, size=density), fill=brewer.pal(3,"Set1")[2], pch=22) +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	facet_wrap(~year) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear_spawners.png"), mapbyyrspawn, width=10, height=8)