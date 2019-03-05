rm(list=ls())

################
## Directories
################

### set to your directory where you will save results
nz_dir <- "C:\\merrill\\stream_networks\\NZ"

### create a general figure directory to save data figures, model comparisons, etc.
fig_dir <- file.path(nz_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

## download 2 required packages for running stream network models
## will first need to install devtools
# install.packages('devtools')

## 1) VAST
## 2) StreamUtils -- supplementary functions to run stream network models in VAST
devtools::install_github("james-thorson/VAST", ref="master")
devtools::install_github("merrillrudd/StreamUtils")

## additional package for aesthetic purposes only
## e.g. plotting function mytheme()
devtools::install_github("merrillrudd/RuddR")

## read in required libraries
library(VAST)
library(StreamUtils)

## should automatically be installed with VAST, but if not can install via CRAN
# install.packages('TMB')
## additional installation tips for Windows available at: https://github.com/kaskr/adcomp/wiki/Windows-installation
library(TMB)	

## libraries for manipulating datasets, plotting, etc.
library(tidyverse)
library(RuddR)

#########################
## Read in data
##########################

## several pre-formatted datasets available in StreamUtils package
## for this example we'll use the Waitaki data with network nodes downstream of observations only
ignore <- data("nz_waitaki_longfin_eel_downstream", package="StreamUtils")

## network information
network <- nz_waitaki_longfin_eel_downstream[["network"]]
network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## Network data in format to input into VAST
## where n_s is the number of network segments including the root node
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))

## observations
obs <- nz_waitaki_longfin_eel_downstream[["observations"]] %>%
   dplyr::filter(data_type=="encounter") %>% 	## use encounter data only
   dplyr::select(-data_type) %>%				## remove data label --> using encounter only
   dplyr::rename('present' = data_value)		## rename data values to 'present'

## 
p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=1, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  xlab("Longitude") + ylab("Latitude") + 
  mytheme()
ggsave(file.path(fig_dir, "Waitaki.png"), p1)

p2 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.01, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat, color = factor(present)), pch=19, cex=1.5,alpha=0.8) +
  scale_color_brewer(palette = "Set1") +
  xlab("Longitude") + ylab("Latitude") + 
  guides(color=guide_legend(title="Encounter")) +
  facet_wrap(~year) +
  mytheme()
ggsave(file.path(fig_dir, "Waitaki_byYear.png"), p2)


##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters_example")
dir.create(nz_enc_dir, showWarnings=FALSE)
setwd(nz_enc_dir)

## save dataset used in this round of model runs within `nz_enc_dir`
## so that all model runs within this folder use the same dataset
saveRDS(obs, file.path(nz_enc_dir, "observations.rds"))
saveRDS(network, file.path(nz_enc_dir, "network.rds"))

## setting up list to look at AIC as we go
AIC_list <- NULL

##### General settings
Version = "VAST_v7_0_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
n_x = nrow(network)   # Specify number of network nodes
## vestigial settings for 2D spatial model --- not actually used here but need to write something
strata.limits <- data.frame('STRATA'="All_areas")
grid_size_km = 1

## add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### vessel effect formatting
vname_i <- unique(obs$fishmethod)
v_i <- sapply(1:nrow(obs), function(x) which(vname_i == obs$fishmethod[x]))
table(obs$fishmethod)

## setup input data frame
## all about observations
Data_Geostat <- data.frame( "Value" = present_new, 
              "Year" = obs$year,
               "Vessel" = obs$fishmethod, ## if no vessel information, write "missing"
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Category" = "Longfin_eels")


## include latitude and longitude for user-supplied area
Extrapolation_List = StreamUtils::make_extrapolation_info( Region="Stream", 
  stream_info=cbind("Lat"=obs$lat, 
                    "Lon"=obs$long,
                    "child_i"=obs$child_i,
                    "Area_km2"=obs$dist_i), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = StreamUtils::make_spatial_info( n_x=n_x, 
                          Method=Method,
                          Network_sz=Network_sz, 
                          # stream_width = xx, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          Lat_x=network$lat, 
                          Lon_x=network$long, 
                          Extrapolation_List=Extrapolation_List )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

########################
## TEMPORAL LOGNORMAL DELTA LINK
########################
temp_dir <- file.path(nz_enc_dir, "temporal_lognormal_deltalink_vesseleffect")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

## turns on or off the spatial (omega) or spatiotemporal (epsilon) variation
## in component 1 (encounter probability for deltalink or numbers density for poisson-link)
## and component 2 (positive catch rates for delta link or average weight for poisson-link)
## encounter/non-encounter data will leave Omega2 and Epsilon2 = 0 because there is no information on catch rates
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)

## governs how temporal (beta) and spatiotemporal (epsilon) intercepts are structured
## for component 1 and component 2
## 0 = fixed effect
## 1 = random effect
## 2 = random walk
## 3 = fixed across years
## with encounter data, Beta2 = 3 to fix intercept for positive catch rates over time
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)


## specifies the observation model distribution function (PosDist)
## and link function - relationship between linear predictor and mean of the distribution function
## PosDist = 0 (normal), 1 (lognormal), or 2 (gamma), other discrete functions in ?VAST::make_data
## Link = 0 (conventional delta-link function) or 1 (Poisson link function) or other options at ?VAST::make_data
ObsModel = c("PosDist"=1, "Link"=0)

## governing any correlated overdispersion among categories for each level of v_i, where eta1 is for encounter probability, and eta2 is for positive catch rates, where 0 is off, "AR1" is an AR1 process, and >0 is the number of elements in a factor-analysis covariance
OverdispersionConfig = c("Eta1"=1, "Eta2"=0)

## derived values
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1,
            "Calculate_proportion"=1)

### outputs data input for VAST
Data = VAST::make_data("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel"=ObsModel,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Value'],
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

plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=TRUE, savedir=fig_dir, cex=0.2)

## and prints structure for each model parameter
## information on model parameters found at ?VAST::Param
TmbList = VAST::make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)

Obj = TmbList[["Obj"]]

## first run of optimization without calculating standard errors or bias correction -- just running the optimizer once
Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

## convergence -- check final gradients
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## convergence -- check standard errors
Opt[["SD"]]
Opt$Convergence_check

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

AIC_list$temporal_lognormal_deltalink <- as.numeric(Opt$AIC)

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
## plot_type = 1 plots residuals for encounter probability
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 ) 

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- StreamUtils::plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, cex_network=0.0001)
# Dens_xt <- StreamUtils::plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.00001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year")

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))


