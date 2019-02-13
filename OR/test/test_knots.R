rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks\\OR"
test_dir <- file.path(main_dir, "test")
setwd(test_dir)


data_dir <- file.path(main_dir, "data")
data_dir2 <- "C:\\merrill\\StreamUtils\\data"

fig_dir <- file.path(main_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("merrillrudd/RuddR")
devtools::install_github("merrillrudd/StreamUtils")

library(VAST)
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
hab <- readRDS(file.path(data_dir, "Siletz_habitat.rds"))

Network_sz = network %>% select(-c("long","lat")) 
Network_sz <- Network_sz[order(Network_sz$child_s),]
Network_sz <- Network_sz %>% mutate(dist_s = dist_s / 1000)

category_names <- unique(obs$survey)

############################################
## model - density observations, 2 surveys
############################################
##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

### for testing purposes
# test_dir <- file.path(main_dir, "test")
# dir.create(test_dir, showWarnings=FALSE)
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

# include latitude and longitude for user-supplied area
Extrapolation_List1 = make_extrapolation_info( Region="Stream", 
  stream_info=cbind("Lat"=obs$lat, 
                    "Lon"=obs$long,
                    "child_i"=obs$child_i,
                    "Area_km2"=1), 
                  strata.limits=strata.limits )

Extrapolation_List = FishStatsUtils::make_extrapolation_info(Region = "User",
    input_grid = cbind("Lat"=obs$lat,
                       "Lon"=obs$long,
                       "Area_km2"=1),
    strata.limits = strata.limits)

Spatial_List1 = make_spatial_info( n_x=n_x,
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          Lat_x=network$lat, 
                          Lon_x=network$long, 
                          Extrapolation_List=Extrapolation_List1, 
                          Save_Results=TRUE)

# ## change latitude and longitude by node, not using Kmeans
Spatial_List = FishStatsUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          "LAT_intensity"=network$lat, 
                          "LON_intensity"=network$long, 
                          Extrapolation_List=Extrapolation_List, 
                          Save_Results=TRUE )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i, "knot_i1"=Spatial_List1$knot_i )

## check that assigned node is either parent or child node for stream segment associated with the observation
check <- correct <- rep(NA, nrow(Data_Geostat))
for(i in 1:nrow(Data_Geostat)){
  obs1 <- obs[i,]
  dat1 <- Data_Geostat[i,]
  net1 <- network %>% filter(child_s == obs1$child_i)
  check[i] <- (dat1$knot_i == net1$child_s | dat1$knot_i == net1$parent_s)
  correct[i] <- dat1$knot_i == obs1$child_i
}
all(check)
length(which(check==TRUE))/length(check)
length(which(correct==TRUE))/length(correct)

####### temporal
temp_dir <- file.path(test_dir, "temporal")
dir.create(temp_dir, showWarnings = FALSE)
setwd(temp_dir)

#####################################
## try turning on spatial variation
#####################################
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)


nodes <- unique(network$child_s)[order(unique(network$child_s))]
n_t <- length(unique(Data_Geostat$Year))
n_x <- nrow(network)
X_xtp <- array(0,dim=c(n_x,n_t,1))
acw_xtp <- sapply(1:length(nodes), function(x){
	cov <- hab %>% filter(covariate == "ACW")
	find <- cov %>% filter(child_i == nodes[x])
	if(nrow(find)==0)

})

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

Data1 = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
                  "s_i"=Data_Geostat[,'knot_i1']-1, 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "a_xl"=Spatial_List1$a_xl, 
                  "MeshList"=Spatial_List1$MeshList, 
                  "GridList"=Spatial_List1$GridList, 
                  "Method"=Spatial_List1$Method, 
                  "Options"=Options, 
                  "Network_sz"=Network_sz )

## set1
TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

# check that all parameters have a reasonable gradient
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Initial run and hessian check
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, savedir=paste0(getwd(),"/"), bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt$diagnostics

## set1
TmbList1 = Build_TMB_Fn("TmbData"=Data1, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List1$loc_x, "Method"=Method)
Obj1 = TmbList1[["Obj"]]

# check that all parameters have a reasonable gradient
Obj1$fn( Obj1$par )
Obj1$gr( Obj1$par )

# Initial run and hessian check
Opt2 = TMBhelper::Optimize( obj=Obj1, lower=TmbList1[["Lower"]], upper=TmbList1[["Upper"]], getsd=FALSE, savedir=paste0(getwd(),"/"), bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt3 = TMBhelper::Optimize( obj=Obj1, startpar=Opt2$par, lower=TmbList1[["Lower"]], upper=TmbList1[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt3$diagnostics
Opt$diagnostics

#### spatial 
space_dir <- file.path(test_dir, "spatial")
dir.create(space_dir, showWarnings = FALSE)
setwd(space_dir)

#####################################
## try turning on spatial variation
#####################################
FieldConfig = c("Omega1"="IID", "Epsilon1"=0, "Omega2"="IID", "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
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

Data1 = Data_Fn("Version"=Version, 
                  "FieldConfig"=FieldConfig, 
                  "OverdispersionConfig"=OverdispersionConfig, 
                  "RhoConfig"=RhoConfig, 
                  "ObsModel"=ObsModel, 
                  "c_i"=as.numeric(Data_Geostat[,"Category"])-1, 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, 
                  "s_i"=Data_Geostat[,'knot_i1']-1, 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "a_xl"=Spatial_List1$a_xl, 
                  "MeshList"=Spatial_List1$MeshList, 
                  "GridList"=Spatial_List1$GridList, 
                  "Method"=Spatial_List1$Method, 
                  "Options"=Options, 
                  "Network_sz"=Network_sz )

## set1
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
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt$diagnostics

## set1
TmbList1 = Build_TMB_Fn("TmbData"=Data1, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List1$loc_x, "Method"=Method)
Obj1 = TmbList1[["Obj"]]

# Change starting values for kappa - was leading to decorrelation at much smaller distances than the minimum distance between sites
Obj1$par[grep("logkappa",names(Obj1$par))]
Obj1$par[grep("logkappa",names(Obj1$par))] = log(1/median(Network_sz[,'dist_s']))

# check that all parameters have a reasonable gradient
Obj1$fn( Obj1$par )
Obj1$gr( Obj1$par )

# Initial run and hessian check
Opt2 = TMBhelper::Optimize( obj=Obj1, lower=TmbList1[["Lower"]], upper=TmbList1[["Upper"]], getsd=FALSE, savedir=paste0(getwd(),"/"), bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt3 = TMBhelper::Optimize( obj=Obj1, startpar=Opt2$par, lower=TmbList1[["Lower"]], upper=TmbList1[["Upper"]], getsd=TRUE, savedir=paste0(getwd(),"/"), bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Opt3$diagnostics

check <- TMBhelper::Check_Identifiable(Obj1)

Report <- Obj1$report()
## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data1, savedir=NULL)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data1, Report=Report, FileName_QQ="Q-Q_plot", DateFile=NULL, savedir=NULL, plot=2) 

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

StreamUtils::plot_maps(plot_set=3, Report=Report, TmbData=Data1, Spatial_List=Spatial_List1, savedir=NULL, category_names = category_names, cex=0.5, alpha=0.8)
