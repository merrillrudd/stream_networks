#################################################
## Many years with no observations in Waitaki catchment
## Turn off spatiotemporal variation
## Use RhoConfig to fix the intercept at the same value for every year
## (Required to do smething with intercepts when there are many years with no data)
## Stripping out time as a place to start
## also needed to fix logkappa1 in order for model to converge/get standard errors
## 		check details around fixing logkappa1
##################################################



rm(list=ls())

#######################
## Directories
#######################

main_dir <- "C:\\merrill\\stream_networks"
data_dir <- file.path(main_dir, "data")
fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

# Save
Date = Sys.Date()
DateFile = paste0(main_dir,'/',Date,'_EEL_Waitaki/')
dir.create(DateFile)
setwd(DateFile)

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
library(ggplot2)
library(RuddR)


#######################
## Load data
#######################

network <- readRDS(file.path(data_dir, 'lower_waitaki_network_renumbered.rds'))
observations <- readRDS(file.path(data_dir, 'waitaki_observations_renumbered.rds'))

years <- unique(observations$year)
sapply(1:length(years), function(x){
	sub <- observations %>% filter(year==years[x])
	sum(sub$present)
})


plot(x = network$long, y = network$lat, col = "gray")
points(x = observations$long, y = observations$lat, col="red", pch=19)

## number of sites with observations
obs_sub <- observations %>% select("parent_i",'child_i','dist_i')
nrow(unique(obs_sub))

## total number of observations
nrow(observations)

## number of segments
nrow(network %>% filter(parent_s!=0))

## number of nodes
nrow(network)

network$long2 <- sapply(1:nrow(network), function(x) ifelse(network$parent_s[x]==0, NA, network$long[which(network$child_s == network$parent_s[x])]))
network$lat2 <- sapply(1:nrow(network), function(x) ifelse(network$parent_s[x]==0, NA, network$lat[which(network$child_s == network$parent_s[x])]))

aa <- ggplot(network %>% filter(lat < -44.9) %>% filter(long > 171)) +
    geom_segment(aes(x = long,y = lat,xend = long2,yend = lat2),arrow=arrow()) +
    # geom_point(aes(x = long, y = lat, color = c_i), cex=5) +
    geom_label(aes(x = long, y = lat, label = child_s), fill = "black", color = "white") + #, nudge_x = -0.02) + 
    # guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()

#######################
## VAST model settings
#######################

#######################################
### Spatial and temporal effects
#######################################
setwd(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Grid", "Mesh", "Stream_network")[3]
grid_size_km = 1
n_x = nrow(network)  # Specify number of stations (a.k.a. "knots")

## spatial intercept, keep "Omega2"=0 and "Epsilon2"=0 for presence/absence data
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)

## lognormally-distributed
## link = 3 because some years have no observations (recommended in user manual)
ObsModel = c("PosDist"=1, "Link"=1)

## region and dataset
Data_Set = "NZ_longfin_eel"
Region = switch(Data_Set, "NZ_longfin_eel"="New_Zealand_streams")
strata.limits <- data.frame('STRATA'="All_areas")

## presence-absence observations
b_i <- observations$present

## add a small amount to presence observations
set.seed(123)
b_i_new <- b_i + rnorm(length(b_i), 0, 0.00001)
b_i_new[which(b_i==0)] <- 0

## saveimage here
## send Jim the following lines to reproduce the problem
# save.image(file="NZeels.Rdata")

## setup geostatistical data
## observations
Data_Geostat <- data.frame( "Catch_KG" = b_i_new, 
							"Year" = observations$year,
							 "Vessel" = "missing", 
							 "AreaSwept_km2" = observations$dist_i, 
							 "Lat" = observations$lat, 
							 "Lon" = observations$long, 
							 "Pass" = 0)

## network - input grid latitude and longitude by node
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
									input_grid=cbind("Lat"=observations$lat, 
													"Lon"=observations$long, 
													"Area_km2"=observations$dist_i), 
									strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
												  Method=Method, 
												  Lon_i=Data_Geostat[,'Lon'], 
												  Lat_i=Data_Geostat[,'Lat'], 
												  "LAT_intensity"=network$lat, 
												  "LON_intensity"=network$long, 
												  Extrapolation_List=Extrapolation_List, 
												  DirPath=DateFile, 
												  Save_Results=TRUE )

## number of knots now equal to number of nodes
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save(Data_Geostat, file=file.path(data_dir, "Data_Geostat.Rdata"))

## plot data
plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile, pch=19 )


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
				"Network_sz"=network %>% select(c('parent_s','child_s','dist_s')) )

# ### First model run
TmbList1 = Build_TMB_Fn("TmbData"=Data, 
						"RunDir"=DateFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj1 <- TmbList1[["Obj"]]

check <- TMBhelper::Check_Identifiable(Obj1)

	# # ## add parameters to be fixed
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
	Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=FALSE, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
SaveResults = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=Data)
save(SaveResults, file=paste0(DateFile,"SaveResults.RData"))

## plot data
# plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )

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
