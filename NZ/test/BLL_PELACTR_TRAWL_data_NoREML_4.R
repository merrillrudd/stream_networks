###############################################################################################################################
##
## This script fits a VAST model to BLL + PELACTR + TRAWL data and does not use REML - Version 4
##
###############################################################################################################################
rm(list=ls()) 

######## Install VAST and TMB
#devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("james-thorson/VAST")
devtools::install_github("kaskr/adcomp/TMB")

######## Define working directories
RootDir = "E:/Combined_dataset_study/BLL_PELACTR_TRAWL_data_NoREML_4/"

######## Load required libraries
library(TMB)
library(VAST)

######## Define the dataset to use
Data_Set = c("Snapper")[1]

######## Define the CPP version to use
Version = get_latest_version( package="VAST" )
#Version = "VAST_v5_3_0"

######## Define the spatial resolution for the model and whether to use a grid or mesh approximation
Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 25
#### Specify the number of stations (a.k.a. "knots")
n_x = 200
#### Specify k-means settings
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )    

######## Define whether to include spatial and spatio-temporal variation, the rank of the covariance among species, 
######## whether its autocorrelated, and whether there is overdispersion
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel_ez = cbind( c(13,14,2), 1 )

######## Decide on which post-hoc calculations to include in the output
Options =  c("SD_site_density"=0, "SD_site_logdensity"=1, "Calculate_Range"=1, "Calculate_evenness"=0, 
	"Calculate_effective_area"=1, "Calculate_Cov_SE"=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0, 
	'normalize_GMRF_in_CPP'=TRUE)

######## Define any potential stratification of results
strata.limits <- data.frame('STRATA'="All_areas")

######## Set the "Region" parameter to "Gulf_of_Mexico"
Region = switch( Data_Set, "Snapper"="Gulf_of_Mexico" )

######## Set the location for saving files
Date = Sys.Date()
DateFile = paste0(RootDir,Date,'/')
dir.create(DateFile)

######## Save all settings for later reference
Record = ThorsonUtilities::bundlelist( c("Version","Method","grid_size_km","n_x","FieldConfig",
	"RhoConfig","OverdispersionConfig","ObsModel_ez","Kmeans_Config","Options") )
save( Record, file=file.path(DateFile,"Record.RData"))
capture.output( Record, file=paste0(DateFile,"Record.txt"))

######## Read in and process data
DF = read.csv( paste0(RootDir,"Red_snapper_dataset.csv") )
Data_Geostat = data.frame( "Catch_KG"=DF[,'Value'], "Year"=DF[,'Year'], "Vessel"=1, 
	"AreaSwept_km2"=DF[,"AreaSwept_km2"], "Lat"=DF[,'Latitude'], "Lon"=DF[,'Longitude'], 
	Type=sapply(as.character(DF[,'Datatype']),FUN=switch,"Encounter"=1,"Count"=2,"Biomass"=3) )
Data_Geostat = na.omit( Data_Geostat )

######## Generate an extrapolation grid
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, 
	observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample=15 )

######## Generate the information used for conducting spatio-temporal parameter estimation, bundled in list `Spatial_List`
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, 
	Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, 
	randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], 
	iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
save( Data_Geostat, file=file.path(DateFile,"Data_Geostat.RData"))
MyPlotDF = cbind( Extrapolation_List[["Data_Extrap"]][,c('Lat','Lon')],
                  'x2i'=Spatial_List$PolygonList$NN_Extrap$nn.idx[,1] )
save( MyPlotDF, file=file.path(DateFile,"Extrapolation_grid_information.RData"))

######## Explore data
tapply( Data_Geostat[,'Catch_KG'], INDEX=list(Data_Geostat[,'Year'],Data_Geostat[,'Type']), FUN=mean )

######## Modify data
TypeSet = c(1,2,3)
Data_Geostat = Data_Geostat[ which(Data_Geostat[,'Type']%in%TypeSet), ]

#### Relabel Type in same order as TypeSet
Data_Geostat[,'Type'] = match(Data_Geostat[,'Type'], TypeSet)

######## Set up the "catchability matrix"
#if( length(TypeSet)==1 ){
#	Q_ik = NULL
#}else{
#	Q_ik = ThorsonUtilities::vector_to_design_matrix( Data_Geostat[,'Type'] )[,-1,drop=FALSE]
#}
Q_ik = ThorsonUtilities::vector_to_design_matrix( Data_Geostat[,'Type'] )[,-3,drop=FALSE]

######## Modify settings
if( 3 %in% TypeSet ){
	TurnOff2 = FALSE
}else{
  	TurnOff2 = TRUE
  	FieldConfig[3:4] = 0
}
ObsModel_ez = ObsModel_ez[TypeSet,]
Aniso = TRUE
Use_REML = FALSE

######## To estimate parameters, build a list of data inputs used for parameter estimation
Data_orig = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, 
	"RhoConfig"=RhoConfig, "ObsModel_ez"=ObsModel_ez, "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], 
	"a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "e_i"=Data_Geostat[,'Type']-1, 
	"s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, 
	"MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method, 
	"Options"=Options, "Q_ik"=Q_ik, "CheckForErrors"=FALSE, "Aniso"=Aniso )

######## Build the TMB object 
TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, 
	"loc_x"=Spatial_List$loc_x, "Method"=Method, "Use_REML"=Use_REML)

######## Modify the Map list 
Map = TmbList_orig$Map
ParHat = TmbList_orig$Obj$env$parList()
Map[["lambda2_k"]] = factor(rep(NA,length(ParHat$lambda2_k)))
if( TurnOff2==TRUE ){
	Map[["beta2_ct"]] = factor(rep(NA,length(ParHat$beta2_ct)))
}

######## Modify the "Random" list (because REML is having trouble for some reason)
Random = TmbList_orig$Random
if( Use_REML==TRUE ) Random = union(Random,"beta1_ct")
if( Use_REML==TRUE & TurnOff2==FALSE ) Random = union(Random,"beta2_ct")

######## Rebuild the TMB object
TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig,
	"loc_x"=Spatial_List$loc_x, "Method"=Method, "Map"=Map, "Random"=Random)
Obj_orig = TmbList_orig[["Obj"]]

######## Conduct simple bug checks
#### Check that no data is missing a likelihood component
Report = Obj_orig$report()
Which = which( Report$LogProb1_i+Report$LogProb2_i == 0 )
if(length(Which)>0) stop("Something is wrong")

#### Check that no fixed effect has a zero gradient
Gr = Obj_orig$gr(Obj_orig$par)
if(any(Gr==0)) stop("Something is wrong")

######## Use a gradient-based nonlinear minimizer to identify maximum likelihood estimates for fixed effects
Opt_orig = TMBhelper::Optimize( obj=Obj_orig, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], 
	getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=1, 
	bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#### Save results
Report_orig = Obj_orig$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj_orig$env$parList(Opt_orig$par), 
	"Data"=Data_orig)
save(Save_orig, file=paste0(DateFile,"Save_orig.RData"))

######## Visualize the spatial distribution of data
SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, 
	Data_Geostat=Data_Geostat, PlotDir=DateFile )

##### Check whether observed encounter frequencies for either low or high probability samples are within the 95% 
##### predictive interval for predicted encounter probability
Enc_prob = SpatialDeltaGLMM::Check_encounter_prob( Report=Report_orig, Data_Geostat=Data_Geostat, DirName=DateFile)

######## Get region-specific settings for plots
MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, 
	"Extrapolation_List"=Extrapolation_List )

######## Define some settings
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

######## Visualize fit to residuals of catch-rates given encounters using a Q-Q plot. A good Q-Q plot will have residuals 
######## along the one-to-one line
setwd(DateFile)
Q = SpatialDeltaGLMM::QQ_Fn( TmbData=Data_orig, Report=Report_orig, FileName_PP="Posterior_Predictive", 
	FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist" )

######## Plot Pearson residuals. If there are visible patterns (areas with consistently positive or negative residuals accross 
######## or within years), then this is an indication of the model "overshrinking" results towards the intercept, and model 
######## results should then be treated with caution
#### 2018-09-04: Failed: "Error in FishStatsUtils::plot_residuals(...) : 
####  Something is wrong with `Q` input"
SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data_orig, 
	Report=Report_orig, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], 
	PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], 
	Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, 
	Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], 
	Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

######## Visualize which direction has faster or slower decorrelation (termed "geometric anisotropy")
SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report_orig, TmbData=Data_orig )

######## Produce an index of abundance
Index = SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=Data_orig, Sdreport=Opt_orig[["SD"]], Year_Set=Year_Set, 
	Years2Include=Years2Include, use_biascorr=TRUE )

######## Visualize outputs from the model
Dens_xt = SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], 
	Report=Report_orig, Sdreport=Opt_orig$SD, PlotDF=MapDetails_List[["PlotDF"]], 
	MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], 
	FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], 
	Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], 
	mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)
CV_xt = SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=c(10), MappingDetails=MapDetails_List[["MappingDetails"]], 
	Report=Report_orig, Sdreport=Opt_orig$SD, PlotDF=MapDetails_List[["PlotDF"]], 
	MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], 
	FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], 
	Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], 
	mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE)

#### Convert to values at each "extrapolation_grid" cell
Dens_et = Dens_xt[MapDetails_List[["PlotDF"]][,'x2i'],]
save(Dens_et, file=paste0(DateFile,"Dens_et_orig.RData"))
CV_et = CV_xt[MapDetails_List[["PlotDF"]][,'x2i'],]
save(CV_et, file=paste0(DateFile,"CV_et_orig.RData"))

######## See if you can detect shifts in distribution or range expansion/contraction
SpatialDeltaGLMM::Plot_range_shifts(Report=Report_orig, TmbData=Data_orig, Sdreport=Opt_orig[["SD"]], 
	Znames=colnames(Data_orig$Z_xm), PlotDir=DateFile, Year_Set=Year_Set)


