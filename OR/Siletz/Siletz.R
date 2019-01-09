rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks\\OR"
data_dir <- file.path(main_dir, "data")

proj_dir <- file.path(main_dir, "Siletz")

Date <- Sys.Date()
DateFile <- paste0(proj_dir, "/", Date)
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

netmap <- ggplot(net_toUse) +
		geom_point(aes(x = long, y = lat)) +
		geom_point(data=net_toUse %>% filter(parent_s == 0), aes(x = long, y = lat), col = "yellow") +
		mytheme()

obs_spawn <- read.csv(file.path(data_dir, "Siletz", "Spawn_Density.csv"))
obs_juv <- read.csv(file.path(data_dir, "Siletz", "Juv_Density.csv"))

obsmap <- netmap + 
		geom_point(data = obs_spawn, aes(x = long, y = lat), col="red", cex=2) +
		geom_point(data = obs_juv, aes(x = long, y = lat), col="blue", cex=2)

obs_toUse <- obs_spawn %>% 
		select("SpawningYear", "Adlt.Mi", "lat", "long") %>%
		rename("year"="SpawningYear", "adults_per_mile"="Adlt.Mi")

#######################
## VAST model settings
#######################

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Grid", "Mesh", "Stream_network")[3]
grid_size_km = 1
n_x = nrow(net_toUse)  # Specify number of stations (a.k.a. "knots")


#######################
## TEMPORAL ONLY
#######################
ModFile <- file.path(DateFile, "Temporal")
dir.create(ModFile)
setwd(ModFile)
## first linear predictor = encounter probability
## second linear predictor = catch rates
## omega = spatial variation
## epsilon = spatiotemporal variation 
## keep "Omega2"=0 and "Epsilon2"=0 for presence/absence data
## on or off with 1 or 0
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)

## beta = temporal structure of intercepts 
## epsilon = temporal structure on spatiotemporal variation
## see user manual for options
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## overdispersion of encounter probability and catch rates
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

## observation error distribution [1] see ?Data_Fn for options
## link function for first linear predictor = 0 logit-link, 1 poisson
ObsModel = c("PosDist"=2, "Link"=1)

Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)


## region and dataset
Data_Set = "OR_salmonids"
Region = switch(Data_Set, "OR_salmonids"="OR_streams")
strata.limits <- data.frame('STRATA'="All_areas")


# save.image(file="NZeels.Rdata")

## setup geostatistical data
## observations
Data_Geostat <- data.frame( "Catch_KG" = obs_toUse$adults_per_mile, 
							"Year" = obs_toUse$year,
							 "Vessel" = "missing", 
							 "AreaSwept_km2" = 1, 
							 "Lat" = obs_toUse$lat, 
							 "Lon" = obs_toUse$long, 
							 "Pass" = 0)

## network - input grid latitude and longitude by node
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
									input_grid=cbind("Lat"=obs_toUse$lat, 
													"Lon"=obs_toUse$long, 
													"Area_km2"=1), 
									strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
												  Method=Method, 
												  Lon_i=Data_Geostat[,'Lon'], 
												  Lat_i=Data_Geostat[,'Lat'], 
												  "LAT_intensity"=net_toUse$lat, 
												  "LON_intensity"=net_toUse$long, 
												  Extrapolation_List=Extrapolation_List, 
												  DirPath=ModFile, 
												  Save_Results=TRUE )


## number of knots now equal to number of nodes
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save(Data_Geostat, file=file.path(ModFile, "Data_Geostat.Rdata"))

## plot data
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=ModFile, pch=19, cex = 1 )

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
				"Network_sz"=net_toUse %>% select(c('parent_s','child_s','dist_s')) )

## First model run
TmbList1 = Build_TMB_Fn("TmbData"=Data, 
						"RunDir"=ModFile, 
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
						"RunDir"=ModFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj = TmbList[["Obj"]]
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=ModFile, bias.correct=TRUE, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
SaveResults = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=Data)
save(SaveResults, file=paste0(ModFile,"SaveResults.RData"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=ModFile)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
	FileName_Phist="Posterior_Predictive-Histogram", 
	FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=ModFile )

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q, savedir=ModFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

## Model selection
Opt$AIC

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( DirName=ModFile, TmbData=Data, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )


#######################
## SPATIAL EFFECTS
#######################
ModFile <- file.path(DateFile, "Spatial")
dir.create(ModFile)
setwd(ModFile)
## first linear predictor = encounter probability
## second linear predictor = catch rates
## omega = spatial variation
## epsilon = spatiotemporal variation 
## keep "Omega2"=0 and "Epsilon2"=0 for presence/absence data
## on or off with 1 or 0
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)

## beta = temporal structure of intercepts 
## epsilon = temporal structure on spatiotemporal variation
## see user manual for options
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## overdispersion of encounter probability and catch rates
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

## observation error distribution [1] see ?Data_Fn for options
## link function for first linear predictor = 0 logit-link, 1 poisson
ObsModel = c("PosDist"=1, "Link"=1)

Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)


## region and dataset
Data_Set = "OR_salmonids"
Region = switch(Data_Set, "OR_salmonids"="OR_streams")
strata.limits <- data.frame('STRATA'="All_areas")


# save.image(file="NZeels.Rdata")

## setup geostatistical data
## observations
Data_Geostat <- data.frame( "Catch_KG" = obs_toUse$adults_per_mile, 
							"Year" = obs_toUse$year,
							 "Vessel" = "missing", 
							 "AreaSwept_km2" = 1, 
							 "Lat" = obs_toUse$lat, 
							 "Lon" = obs_toUse$long, 
							 "Pass" = 0)

## network - input grid latitude and longitude by node
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
									input_grid=cbind("Lat"=obs_toUse$lat, 
													"Lon"=obs_toUse$long, 
													"Area_km2"=1), 
									strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
												  Method=Method, 
												  Lon_i=Data_Geostat[,'Lon'], 
												  Lat_i=Data_Geostat[,'Lat'], 
												  "LAT_intensity"=net_toUse$lat, 
												  "LON_intensity"=net_toUse$long, 
												  Extrapolation_List=Extrapolation_List, 
												  DirPath=ModFile, 
												  Save_Results=TRUE )


## number of knots now equal to number of nodes
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save(Data_Geostat, file=file.path(ModFile, "Data_Geostat.Rdata"))

## plot data
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=ModFile, pch=19, cex = 1 )

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
				"Network_sz"=net_toUse %>% select(c('parent_s','child_s','dist_s')) )

## First model run
TmbList1 = Build_TMB_Fn("TmbData"=Data, 
						"RunDir"=ModFile, 
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
						"RunDir"=ModFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj = TmbList[["Obj"]]
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=ModFile, bias.correct=FALSE, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
SaveResults = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=Data)
save(SaveResults, file=paste0(ModFile,"SaveResults.RData"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=ModFile)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
	FileName_Phist="Posterior_Predictive-Histogram", 
	FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=ModFile )

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q, savedir=ModFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

## Model selection
Opt$AIC

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( DirName=ModFile, TmbData=Data, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

