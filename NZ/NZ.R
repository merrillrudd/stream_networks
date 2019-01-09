rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks\\NZ"
data_dir <- file.path(main_dir, "data")

Date <- Sys.Date()
DateFile <- paste0(main_dir, "/", Date)
dir.create(DateFile, showWarnings=FALSE)
setwd(DateFile)


#################
## Packages
#################

devtools::install_github("james-thorson/VAST", ref="development")
library(VAST)
devtools::install_github("merrillrudd/FishStatsUtils")
library(FishStatsUtils)

library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

################
## Load data
################

## all network
load(file.path(data_dir, "REC2.4fromGDB.Rdata"))
network_raw <- REC2.4fromGDB

network <- network_raw %>%
	select(c('CatName','nzsegment','fnode','tnode','Shape_Leng', 'upcoordX', 'upcoordY', 'downcoordX', 'downcoordY','NextDownID','Headwater','REC2_TerminalSegment')) %>%
	rename('parent_s' = tnode, 'child_s' = fnode, 'dist_s'=Shape_Leng, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY, 'segment'=nzsegment, 'NextDownSeg'=NextDownID, 'Headwater'=Headwater)

## all observations
load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
obs_raw <- NZFFD.REC2.Diad.EF

obs <- obs_raw %>% 
	select(c('catchname', 'nzsegment', 'angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y',"Headwater")) %>%
	rename('catchment'=catchname, 'segment'=nzsegment, 'present'=angdie, 'northing_child'=upcoordX, 'easting_child'=upcoordY, 'northing_parent'=downcoordX, 'easting_parent'=downcoordY, 'year'=y,"Headwater"=Headwater) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit()

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
obs <- obs %>% filter(segment %in% network$segment == TRUE)

obsmap <- ggplot() +
		geom_point(data=network, aes(x = long_child, y = lat_child), col = "black") +
		geom_point(data=obs, aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()

#############################
## subset Waitaki catchment
#############################

network_sub <- network %>% filter(grepl("Waitaki", CatName))
obs_sub <- obs %>% filter(grepl("Waitaki", catchment))
obs_sub <- obs_sub %>% filter(segment %in% network_sub$segment == TRUE)

submap <- obsmap +
		geom_point(data=network_sub, aes(x = long_child, y = lat_child), col="blue") +
		geom_point(data=obs_sub, aes(x = long_child, y = lat_child), col = "yellow")

catchmap <- ggplot() +
		geom_point(data=network_sub, aes(x = long_child, y = lat_child), col="black") +
		geom_point(data=obs_sub, aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()

network_sub_cut <- network_sub %>% 
			filter(long_parent >= min(c(obs_sub$long_child,obs_sub$long_parent))) %>%
			filter(lat_parent <= max(c(obs_sub$lat_child,obs_sub$lat_parent)))
all(obs_sub$segment %in% network_sub_cut$segment)

catchmap2 <- ggplot() +
		geom_point(data=network_sub_cut, aes(x = long_child, y = lat_child), col="black") +
		geom_point(data=obs_sub, aes(x = long_child, y = lat_child), col = "red") +
		xlab("Longitude") + ylab("Latitude") +
		mytheme()

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
obs_reformat <- inner_join(obs_sub, network_sub_cut) %>% select('parent_s', 'child_s', 'dist_s', 'lat_child', 'long_child', 'lat_parent', 'long_parent', "NextDownSeg", 'present', 'year')
obs_reformat$present <- sapply(1:nrow(obs_reformat), function(x) ifelse(obs_reformat$present[x]==FALSE, 0, 1))

obs_parents <- sapply(1:nrow(obs_reformat), function(x) inodes[which(nodes == obs_reformat$parent_s[x])])
obs_children <- sapply(1:nrow(obs_reformat), function(x) inodes[which(nodes == obs_reformat$child_s[x])])
obs_reformat$parent_s <- obs_parents
obs_reformat$child_s <- obs_children

obs_toUse <- obs_reformat %>% select(-c('lat_parent', 'long_parent', 'NextDownSeg'))%>%
			rename('lat'=lat_child, 'long'=long_child) %>%
			rename('parent_i' = parent_s, 'child_i' = child_s, 'dist_i' = dist_s)

#######################
## VAST model settings
#######################
setwd(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network_toUse)  # Specify number of stations (a.k.a. "knots")

## spatial intercept, keep "Omega2"=0 and "Epsilon2"=0 for presence/absence data
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## lognormally-distributed
## link = 3 because some years have no observations (recommended in user manual)
ObsModel = c("PosDist"=2, "Link"=1)

## region and dataset
Data_Set = "NZ_longfin_eel"
Region = switch(Data_Set, "NZ_longfin_eel"="New_Zealand_streams")
strata.limits <- data.frame('STRATA'="All_areas")

## presence-absence observations
b_i <- obs_toUse$present

## add a small amount to presence observations
set.seed(123)
b_i_new <- b_i + rnorm(length(b_i), 0, 0.0001)
b_i_new[which(b_i==0)] <- 0

# save.image(file="NZeels.Rdata")

## setup geostatistical data
## observations
Data_Geostat <- data.frame( "Catch_KG" = b_i_new, 
							"Year" = obs_toUse$year,
							 "Vessel" = "missing", 
							 "AreaSwept_km2" = obs_toUse$dist_i, 
							 "Lat" = obs_toUse$lat, 
							 "Lon" = obs_toUse$long, 
							 "Pass" = 0)

## network - input grid latitude and longitude by node
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
									input_grid=cbind("Lat"=obs_toUse$lat, 
													"Lon"=obs_toUse$long, 
													"Area_km2"=obs_toUse$dist_i), 
									strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
												  Method=Method, 
												  Lon_i=Data_Geostat[,'Lon'], 
												  Lat_i=Data_Geostat[,'Lat'], 
												  "LAT_intensity"=network_toUse$lat, 
												  "LON_intensity"=network_toUse$long, 
												  Extrapolation_List=Extrapolation_List, 
												  DirPath=DateFile, 
												  Save_Results=TRUE )


## number of knots now equal to number of nodes
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save(Data_Geostat, file=file.path(data_dir, "Data_Geostat.Rdata"))

## plot data
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile, pch=19, cex=2 )

# Network_sz <- unique(network %>% select(c('parent_s','child_s','dist_s')))
## data file for TMB
Data = Data_Fn( "Version"=Version, 
				"FieldConfig"=FieldConfig, 
				"OverdispersionConfig"=OverdispersionConfig, 
				"RhoConfig"=RhoConfig, 
				"ObsModel"=ObsModel, 
				"c_iz"=rep(0,nrow(Data_Geostat)), 
				"b_i"=Data_Geostat[,'Catch_KG'], 
				"a_i"=Data_Geostat[,'AreaSwept_km2'], 
				"v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
				"s_i"=Data_Geostat[,'knot_i']-1, 
				"t_iz"=Data_Geostat[,'Year'], 
				"a_xl"=Spatial_List$a_xl, 
				"MeshList"=Spatial_List$MeshList, 
				"GridList"=Spatial_List$GridList, 
				"Method"=Spatial_List$Method, 
				"Options"=Options, 
				"Network_sz"=network_toUse %>% select(c('parent_s','child_s','dist_s')) )

## First model run
TmbList1 = Build_TMB_Fn("TmbData"=Data, 
						"RunDir"=DateFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj1 <- TmbList1[["Obj"]]

check <- TMBhelper::Check_Identifiable(Obj1)

## add parameters to be fixed
Map_custom <- TmbList1$Map
	# Map_custom[["logkappa1"]] <- factor(NA)
	# Map_custom[["beta1_ct"]] <- factor(rep(NA, length(TmbList1$Parameters$beta1_ct)))

TmbList = Build_TMB_Fn("TmbData"=Data, 
						"Map" = Map_custom,
						"RunDir"=DateFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj = TmbList[["Obj"]]
Opt = tryCatch(TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") ), error=function(e) NA)

i <- 0
while(any(grepl("opt", names(Opt))) | all(is.na(Opt))){
  i <- i + 1
  TmbList = Build_TMB_Fn("TmbData"=Data, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Opt = tryCatch(TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") ), error=function(e) NA)
  # names <- names(Opt_orig)
}

Report = Obj$report()
SaveResults = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=Data)
save(SaveResults, file=paste0(DateFile,"SaveResults.RData"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
	FileName_Phist="Posterior_Predictive-Histogram", 
	FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile )

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

## Model selection
Opt$AIC

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.2, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( DirName=DateFile, TmbData=Data, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

