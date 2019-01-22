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

network <- readRDS(file.path(data_dir, "Siletz_network.rds"))
obs <- readRDS(file.path(data_dir, "Siletz_observations_density.rds"))

Network_sz = network %>% select(-c("long","lat"))

############################################
## model - density observations, 2 surveys
############################################
sil_dir <- file.path(main_dir, "siletz_density")
dir.create(sil_dir)
setwd(sil_dir)

##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

test_dir <- file.path(main_dir, "test")
dir.create(test_dir, showWarnings=FALSE)

# save.image(file.path(test_dir,"OR_test.Rdata"))

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = obs$density, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Category" = obs$surveynum)

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

##################
## 2 surveys
##################
surv2_dir <- file.path(sil_dir, "spatial_effect")
dir.create(surv2_dir, showWarnings = FALSE)
setwd(surv2_dir)

#####################################
## try turning on spatial variation
#####################################
FieldConfig = c("Omega1"="IID", "Epsilon1"=0, "Omega2"="IID", "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,1)

Data = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
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

category_names <- unique(obs$survey)
p1 <- plot_network(Extrapolation_List = Extrapolation_List, Spatial_List = Spatial_List, TmbData=Data, Data_Geostat = Data_Geostat, category_names=category_names, observations=TRUE, arrows=TRUE, root=TRUE, plot_type=1)
p2 <- plot_network(Extrapolation_List = Extrapolation_List, Spatial_List = Spatial_List, TmbData=Data, Data_Geostat = Data_Geostat, category_names=category_names, observations=TRUE, root=TRUE, plot_type=2)


TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

# Change starting values for kappa - was leading to decorrelation at much smaller distances than the minimum distance between sites
Obj$par[grep("logkappa",names(Obj$par))]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

# check that all parameters have a reasonable gradient
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Initial run and hessian check
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=paste0(getwd(),"/"), bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
H = optimHess( par=Opt1$par, fn=Obj$fn, gr=Obj$gr )

# Re-run from last MLE
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

Report = Obj$report()

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data, DirName=NULL)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, save_dir=NULL, plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

plot_maps(plot_set=3, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, alpha=0.7, cex=0.8)

plot_maps(plot_set=12, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, alpha=0.7, cex=0.8)
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, DirName=NULL, strata_names="Longfin eel", "category_names"=category_names )

#####################################
## try turning on spatiotemporal variation
#####################################
mod2 <- file.path(sil_dir, "spatiotemporal_effect")
dir.create(mod2)
setwd(mod2)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

Data = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
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

# Change starting values for kappa - was leading to decorrelation at much smaller distances than the minimum distance between sites
Obj$par[grep("logkappa",names(Obj$par))]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

# check that all parameters have a reasonable gradient
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Initial run and hessian check
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=paste0(getwd(),"/"), bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
H = optimHess( par=Opt1$par, fn=Obj$fn, gr=Obj$gr )

# Re-run from last MLE
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]

Report = Obj$report()

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report, Data=Data, DirName=NULL)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, save_dir=NULL, plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=1 )

## variance for positive catch rates NA
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, savedir=NULL, plot_type=2 )

plot_maps(plot_set=3, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, alpha=0.7, cex=0.8)

plot_maps(plot_set=12, Report=Report, TmbData=Data, Spatial_List=Spatial_List, savedir=NULL, category_names = category_names, alpha=0.7, cex=0.8)
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, DirName=NULL, strata_names="Longfin eel", "category_names"=category_names )
