rm(list=ls())

###########################
## set working directory
###########################
setwd("C:\\merrill\\stream_networks\\OR\\Siletz_demo")

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

Network_sz_LL <- network %>% rename("Lon"=long, "Lat"=lat)
Network_sz <- Network_sz_LL %>% select('parent_s','child_s','dist_s')
#######################################
## Spawners only
#######################################
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
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)



