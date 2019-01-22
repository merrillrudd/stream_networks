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

# devtools::install_github("james-thorson/VAST", ref="development")
library(VAST)
devtools::install_github("merrillrudd/StreamUtils")
library(StreamUtils)

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

net <- network %>% 
	select("Shape_Leng", "ChildNode", "ParentNode", "lat", "long") %>%
	rename("dist_s"="Shape_Leng", "child_s"="ChildNode", "parent_s"="ParentNode")
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

obs_spawn <- read.csv(file.path(data_dir, "Siletz", "Spawn_Density.csv"))
obs_juv <- read.csv(file.path(data_dir, "Siletz", "Juv_Density.csv"))

obs_spawn <- obs_spawn %>% 
		select("SpawningYear", "Adlt.Mi", "lat", "long") %>%
		rename("year"="SpawningYear", "density"="Adlt.Mi") %>%
		mutate("density"=density / 1.6) %>%
		mutate("survey"="spawners") %>%
		mutate('surveynum' = 1)

obs_juv <- obs_juv %>%
		select("Year", "CohoPerKilometer", "lat", "long") %>%
		rename("year"="Year", "density"="CohoPerKilometer") %>%
		na.omit() %>%
		mutate("survey"="juveniles") %>% 
		mutate("surveynum" = 2)

obs_toUse <- rbind.data.frame(obs_spawn, obs_juv)

years <- unique(obs_toUse$year)[order(unique(obs_toUse$year))]
net_byYear <- lapply(1:length(years), function(x){
	netyr <- data.frame(net_toUse, "year"=years[x])
	return(netyr)
})
net_byYear <- do.call(rbind, net_byYear)

saveRDS(obs_toUse, file.path(data_dir, "Siletz_observations_density.rds"))
saveRDS(net_toUse, file.path(data_dir, "Siletz_network.rds"))


# or_siletz_coho <- list()
# or_siletz_coho$observations <- obs_toUse
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
    geom_point(data = obs_toUse, aes(x = long, y = lat, fill=survey), pch=22, cex=3) +
    scale_fill_brewer(palette = "Set1")+
    xlab("Longitude") + ylab("Latitude") +  
    mytheme()
ggsave(file.path(fig_dir, "Network_observations.png"), bb, width=8, height=6)

years <- unique(obs_toUse$year)[order(unique(obs_toUse$year))]
net_byYear <- lapply(1:length(years), function(x){
	netyr <- data.frame(net_toUse, "year"=years[x])
	return(netyr)
})
net_byYear <- do.call(rbind, net_byYear)

mapbyyr <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat)) +
	geom_point(data=obs_toUse, aes(x = long, y = lat, fill=survey), cex=3, pch=22) +
	facet_wrap(~year) +
	scale_fill_brewer(palette = "Set1") +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear.png"), mapbyyr, width=10, height=8)

mapbyyrjuv <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat)) +
	geom_point(data=obs_toUse %>% filter(survey=="juveniles"), aes(x = long, y = lat, size=density), fill=brewer.pal(3,"Set1")[1], pch=22) +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	facet_wrap(~year) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear_juveniles.png"), mapbyyrjuv, width=10, height=8)

mapbyyrspawn <- ggplot() +
    geom_point(data = l2 %>% filter(parent_s==0), aes(x = long, y = lat), color="goldenrod", cex=5) +
	geom_point(data=net_byYear, aes(x = long, y = lat), cex=2) +
	geom_point(data=obs_toUse %>% filter(survey=="spawners"), aes(x = long, y = lat, size=density), fill=brewer.pal(3,"Set1")[2], pch=22) +
	scale_x_continuous(breaks = quantile(net_byYear$long,prob=c(0.1,0.5,0.9)), labels=round(quantile(net_byYear$long,prob=c(0.1,0.5,0.9))[1:3],1)) +
	facet_wrap(~year) +
	xlab("Longitude") + ylab("Latitude") +
	mytheme()
ggsave(file.path(fig_dir, "Observations_byYear_spawners.png"), mapbyyrspawn, width=10, height=8)