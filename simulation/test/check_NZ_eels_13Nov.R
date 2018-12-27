rm(list=ls())
devtools::install_github("james-thorson/VAST", ref="development")

library(VAST)
library(TMB)
library(dplyr)

setwd( "C:\\merrill\\stream_networks\\test" )

load("NZeels.Rdata")


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
# TmbList1 = Build_TMB_Fn("TmbData"=Data, 
# 						"RunDir"=DateFile, 
# 						"Version"=Version, 
# 						"RhoConfig"=RhoConfig, 
# 						"loc_x"=Spatial_List$loc_x, 
# 						"Method"=Method)

# Obj1 <- TmbList1[["Obj"]]

# check <- TMBhelper::Check_Identifiable(Obj1)

# # ## add parameters to be fixed
# Map_custom <- TmbList1$Map

TmbList = Build_TMB_Fn("TmbData"=Data, 
						"RunDir"=DateFile, 
						"Version"=Version, 
						"RhoConfig"=RhoConfig, 
						"loc_x"=Spatial_List$loc_x, 
						"Method"=Method)

Obj = TmbList[["Obj"]]
Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
SaveResults = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "Data"=Data)
save(SaveResults, file=paste0(DateFile,"SaveResults.RData"))
