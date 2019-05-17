rm(list=ls())

###########################
## set working directory
###########################
dir <- "C:\\merrill\\stream_networks\\OR\\Siletz_demo"
setwd(dir)

## read data prepared for analysis
data <- load("or_siletz_coho.rda")

network <- or_siletz_coho[["network"]]
obs <- or_siletz_coho[["observations"]]

####################
## load packages
####################

devtools::install_github("merrillrudd/FishStatsUtils", ref="stream")
library(FishStatsUtils)

devtools::install_github("james-thorson/VAST")
library(VAST)

library(TMB)
library(tidyverse)

###########################
## Observations for VAST
###########################

## counts per kilometer, keep AreaSwept = 1
Data_Geostat <- data.frame( "Catch_KG" = obs$value, 
              "Year" = as.numeric(obs$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = obs$survey,
               "CategoryNum"=obs$surveynum)
category_names <- unique(Data_Geostat$Category)

Network_sz_LL <- network %>% rename("Lon"=long, "Lat"=lat)
Network_sz <- Network_sz_LL %>% select('parent_s','child_s','dist_s')

#######################################
## Spawners only
#######################################
path <- file.path(dir, "Spawners")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.o"), to = path)


Data <- Data_Geostat %>% filter(CategoryNum == 1)

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=2,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0.nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

## make sure none of the parameters are hitting their bounds
fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                "Lat_i"=Data[,"Lat"], 
                "Lon_i"=Data[,"Lon"], 
                "t_iz"=Data[,'Year'], 
                "c_i"=rep(0.nrow(Data)), 
                "b_i"=Data[,'Catch_KG'], 
                "a_i"=Data[,'AreaSwept_km2'], 
                "v_i"=Data[,'Vessel'], 
                working_dir = path,
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

 saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds"))    

 
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=path)

Index = plot_biomass_index( DirName=path, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1] )

plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=path, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=path )

Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=path, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=path, FileName=path, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)

#######################################
## Juveniles only
#######################################
path <- file.path(dir, "Juveniles")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings=FALSE)

ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(dir, "VAST_v8_0_0.o"), to = path)

# Density of spawners in juvenile year - 1
spawn <- fit$Report$D_gcy

Data <- Data_Geostat %>% filter(CategoryNum == 2)

n_x <- nrow(Network_sz_LL)
n_t <- length(unique(Data$Year))
n_i <- nrow(Data)

X_gtp_input <- array(0, dim=c(n_x,n_t,1))
X_gtp_input[,2:n_t,1] <- spawn[,1,]

X_itp_input <- array(X_gtp_input[which(Data$CategoryNum==2),,1], dim=c(n_i,n_t,1))

## turn on spatial and spatiotemporal effects
FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")

## IID structure on temporal intercepts
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## gamma distribution, conventional delta link model
ObsModel = c("PosDist"=2,"Link"=0)

## other options
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_input, X_itp = X_itp_input,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,'Vessel'], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_input, X_itp = X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

## make sure none of the parameters are hitting their bounds
fit1$parameter_estimates$diagnostics

## run the model
fit = fit_model( "settings"=settings, 
                "Lat_i"=Data[,"Lat"], 
                "Lon_i"=Data[,"Lon"], 
                "t_iz"=Data[,'Year'], 
                "c_i"=rep(0,nrow(Data)), 
                "b_i"=Data[,'Catch_KG'], 
                "a_i"=Data[,'AreaSwept_km2'], 
                "v_i"=Data[,'Vessel'], 
                working_dir = path,
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_input, X_itp = X_itp_input,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

 saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds"))    

 
map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=path)

Index = plot_biomass_index( DirName=path, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names[2] )

plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=path, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[2])

Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=path )

Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=path, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[2], Cex=0.5)

plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=path, FileName=path, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)
