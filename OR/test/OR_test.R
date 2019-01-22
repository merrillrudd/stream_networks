rm(list=ls())

devtools::install_github("james-thorson/VAST", ref = "development")

library(VAST)
library(FishStatsUtils)
library(TMB)

setwd("C:\\merrill\\stream_networks\\OR\\test")

load("OR_test.Rdata")

##### setup data frame
## density in adults or juveniles per kilometer
Data_Geostat <- data.frame( "Catch_KG" = obs$density, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "spp" = obs$surveynum)

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
                  "c_i"=as.numeric(Data_Geostat[,"spp"])-1, 
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
check <- Check_Identifiable2(Obj)

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Opt1$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

##################
# Problem
#  kappa1 and kapp2 parameters are not identifiable and final gradients are very tiny (1e-61)
#
