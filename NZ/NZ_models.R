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

devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("merrillrudd/FishStatsUtils")
devtools::install_github("merrillrudd/RuddR")
devtools::install_github("merrillrudd/StreamUtils")

library(VAST)
library(FishStatsUtils)
library(StreamUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

#########################
## read in data
##########################

# network <- readRDS(file.path(data_dir, "Waitaki_network.rds"))
# obs <- readRDS(file.path(data_dir, "Waitaki_observations_encounters.rds"))

data <- data("nz_waitaki_longfin_eel", package="StreamUtils")

network <- nz_waitaki_longfin_eel[["network"]]
obs <- nz_waitaki_longfin_eel[["observations"]]

Network_sz = network %>% select(-c("long","lat"))

##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters")
dir.create(nz_enc_dir)
setwd(nz_enc_dir)


##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.00001)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "eel_encounters" = present_new, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0)


## include latitude and longitude for user-supplied area
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs$lat, 
                    "Lon"=obs$long,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$lat, 
                          "LON_intensity"=network$long, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# save.image("NZ_test.Rdata")

########################
## No spatial variation
########################
temp_dir <- file.path(nz_enc_dir, "temporal_only")
dir.create(temp_dir)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)


present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]>0.9, present[x]+devs[x], present[x]))

Data_Geostat[,"eel_encounters"] <- present_new

Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'eel_encounters'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="User", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(1,3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

# ## index of abundance
Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

##################
# I'll try turning on spatial variation
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_effects")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)

OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0,
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'eel_encounters'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)
#############
# It has L_omega1_z going to 0 and gradient for kappa is going to zero, so the model is hitting a bound and wants to turn on spatial variation.  So obviously there's too little informatino here for a spatial model.
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]

Report <- Obj$report()

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="User", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

# ## index of abundance
Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

#####################
## spatiotemporal

space_dir <- file.path(nz_enc_dir, "spatiotemporal_effects")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)

OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0,
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'eel_encounters'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=Spatial_List$a_xl,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)
#############
# It has L_omega1_z going to 0 and gradient for kappa is going to zero, so the model is hitting a bound and wants to turn on spatial variation.  So obviously there's too little informatino here for a spatial model.
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]

Report <- Obj$report()

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="User", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

# ## index of abundance
Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )
