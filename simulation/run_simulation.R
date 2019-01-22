rm(list=ls())

# devtools::install_github("james-thorson/VAST", ref = "development")
# devtools::install_github("merrillrudd/FishStatsUtils")
library(VAST)
library(FishStatsUtils)
library(TMB)
library(tidyverse)
library(RuddR)

wd <- "C:\\merrill\\stream_networks\\examples"
setwd(wd)

data_dir <- file.path(wd, "data")
dir.create(data_dir)

fig_dir <- file.path(wd, "figures")
dir.create(fig_dir)

res_dir <- file.path(wd, "results")
dir.create(res_dir)


#########################
## read in simulated data
##########################

network <- readRDS(file.path(data_dir, "Simulated_network.rds"))
obs_catch <- readRDS(file.path(data_dir, "Simulated_observations_catch.rds"))
obs_enc <- readRDS(file.path(data_dir, "Simulated_observations_encounters.rds"))
obs_catch_samp <- readRDS(file.path(data_dir, "Simulated_observations_catch_sampled.rds"))
obs_enc_samp <- readRDS(file.path(data_dir, "Simulated_observations_encounters_sampled.rds"))

###############################
## model - catch observations
###############################
sim_catch_dir <- file.path(res_dir, "sim_catch")
dir.create(sim_catch_dir)
setwd(sim_catch_dir)


Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")


FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=1, "Epsilon2"=1)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,1)


Data_Geostat <- data.frame( "Catch_KG" = obs_catch$catch, 
              "Year" = obs_catch$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_catch$Lat, 
               "Lon" = obs_catch$Lon, 
               "Pass" = 0)

## typically use observations latitude, longitude, and distance sampled
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs_catch$Lat, 
                    "Lon"=obs_catch$Lon,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$Lat, 
                          "LON_intensity"=network$Lon, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

Network_sz = network %>% rename("parent_s"=parent_x, "child_s"=child_x, "dist_s"=dist_x) %>% select(-c("Lon","Lat"))
Data = Data_Fn("Version"=Version, 
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
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=2, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

###############################
## model - catch observations - fewer sites sampled
###############################
sim_catch_samp_dir <- file.path(res_dir, "sim_catch_sampled")
dir.create(sim_catch_samp_dir)
setwd(sim_catch_samp_dir)


Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")


FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=1, "Epsilon2"=1)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,1)


Data_Geostat <- data.frame( "Catch_KG" = obs_catch_samp$catch, 
              "Year" = obs_catch_samp$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_catch_samp$Lat, 
               "Lon" = obs_catch_samp$Lon, 
               "Pass" = 0)

## typically use observations latitude, longitude, and distance sampled
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs_catch_samp$Lat, 
                    "Lon"=obs_catch_samp$Lon,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$Lat, 
                          "LON_intensity"=network$Lon, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

Network_sz = network %>% rename("parent_s"=parent_x, "child_s"=child_x, "dist_s"=dist_x) %>% select(-c("Lon","Lat"))
Data = Data_Fn("Version"=Version, 
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
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=2, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

 
##################################
## model - encounter observations
##################################
sim_enc_dir <- file.path(res_dir, "sim_encounters")
dir.create(sim_enc_dir)
setwd(sim_enc_dir)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")


FieldConfig = c("Omega1"=0, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=1, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,1)


Data_Geostat <- data.frame( "Catch_KG" = obs_enc$catch, 
              "Year" = obs_enc$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_enc$Lat, 
               "Lon" = obs_enc$Lon, 
               "Pass" = 0)

## typically use observations latitude, longitude, and distance sampled
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs_enc$Lat, 
                    "Lon"=obs_enc$Lon,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$Lat, 
                          "LON_intensity"=network$Lon, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

Network_sz = network %>% rename("parent_s"=parent_x, "child_s"=child_x, "dist_s"=dist_x) %>% select(-c("Lon","Lat"))
Data = Data_Fn("Version"=Version, 
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
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=FALSE, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
# Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

##################################
## model - encounter observations -- sampled
##################################
sim_enc_samp_dir <- file.path(res_dir, "sim_encounters_sampled")
dir.create(sim_enc_samp_dir)
setwd(sim_enc_samp_dir)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")


FieldConfig = c("Omega1"=0, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,1)


Data_Geostat <- data.frame( "Catch_KG" = obs_enc_samp$catch, 
              "Year" = obs_enc_samp$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_enc_samp$Lat, 
               "Lon" = obs_enc_samp$Lon, 
               "Pass" = 0)

## typically use observations latitude, longitude, and distance sampled
Extrapolation_List = FishStatsUtils::make_extrapolation_info( Region="User", 
  input_grid=cbind("Lat"=obs_enc_samp$Lat, 
                    "Lon"=obs_enc_samp$Lon,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$Lat, 
                          "LON_intensity"=network$Lon, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

Network_sz = network %>% rename("parent_s"=parent_x, "child_s"=child_x, "dist_s"=dist_x) %>% select(-c("Lon","Lat"))
Data = Data_Fn("Version"=Version, 
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
                  "Network_sz"=Network_sz )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE)#, bias.correct=FALSE, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, pch=19, cex=3 )

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist")

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data, Report=Report, Q=Q,  MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
# Index = plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )
