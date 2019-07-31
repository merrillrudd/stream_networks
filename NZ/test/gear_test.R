rm(list=ls())
setwd("C:\\merrill\\stream_networks\\NZ\\test")

load("NZ_eel_gear.Rdata")

library(VAST)
library(FishStatsUtils)
library(TMB)

  # check model setup and parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run --- system is exactly singular
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))

  ## doesn't get this far bc system is exactly singular for first model run
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
