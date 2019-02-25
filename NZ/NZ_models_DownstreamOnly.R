rm(list=ls())

################
## Directories
################

nz_dir <- "C:\\merrill\\stream_networks\\NZ"

fig_dir <- file.path(nz_dir, "figures")
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

data <- data("nz_waitaki_longfin_eel", package="StreamUtils")

network <- nz_waitaki_longfin_eel[["network"]]
obs <- nz_waitaki_longfin_eel[["observations"]] %>%
    dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value)

Network_sz = network %>% select(-c("long","lat"))

p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.2, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  mytheme()

### prune upstream segments
## do not automatically include observed segment and then make sure it still shows up on network when tracing the upstream and downstream segments from the observed segment
## start from observation
obs_child <- unique(obs$child_i)

finddown <- lapply(1:length(obs_child), function(x){
  net_obs <- network %>% filter(child_s %in% obs_child[x])
  nextdown <- network %>% filter(child_s == net_obs$parent_s)
  save <- rbind.data.frame(net_obs,nextdown)
  for(i in 1:200){
    nextdown <- network %>% filter(child_s == nextdown$parent_s)
    save <- rbind.data.frame(save, nextdown)
    if(any(save$parent_s == 0)) break
  }
  return(save) 
})
check1 <- sapply(1:length(finddown), function(x) nrow(finddown[[x]]))
check2 <- sapply(1:length(finddown), function(x) any(finddown[[x]]$parent_s == 0))
prunedown <- do.call(rbind, finddown)
prunedown2 <- unique(prunedown)

# findup <- lapply(1:length(obs_child), function(x){
#   net_obs <- network %>% filter(child_s %in% obs_child[x])
#   nextup <- network %>% filter(parent_s == net_obs$child_s)
#   save <- rbind.data.frame(net_obs, nextup)
#   for(i in 1:100){
#     nextup <- network %>% filter(parent_s %in% nextup$child_s)
#     save <- rbind.data.frame(save, nextup)
#   }
#   return(save)
# })
# pruneup <- do.call(rbind, findup)
# pruneup2 <- unique(pruneup)
# head_children <- network$child_s[which(network$child_s %in% network$parent_s ==FALSE)]
# all(head_children %in% pruneup2$child_s)
# length(which(head_children %in% pruneup2$child_s))/length(head_children)

p1 <- ggplot() +
  geom_point(data = prunedown2, aes(x = long, y = lat), pch=19, cex=0.5, alpha=0.8) +
  # geom_point(data = pruneup2, aes(x = long, y = lat), pch=19, cex=0.5, alpha=0.8) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  mytheme()

network2 <- rbind.data.frame(prunedown2)#, pruneup2)

## rename nodes
nodes <- unique(c(network2$child_s, network2$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network2), function(x){
  if(network2$parent_s[x] != 0) new_node <- inodes[which(nodes == network2$parent_s[x])]
  if(network2$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network2), function(x) inodes[which(nodes == network2$child_s[x])])

network_toUse <- network2
network_toUse$parent_s <- net_parents
network_toUse$child_s <- net_children

obs_parents <- sapply(1:nrow(obs), function(x){
  if(obs$parent_i[x] != 0) new_node <- inodes[which(nodes == obs$parent_i[x])]
  if(obs$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs), function(x) inodes[which(nodes == obs$child_i[x])])

obs_toUse <- obs
obs_toUse$parent_i <- obs_parents
obs_toUse$child_i <- obs_children


##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters_DownstreamOnly")
dir.create(nz_enc_dir, showWarnings=FALSE)
setwd(nz_enc_dir)

saveRDS(obs_toUse, file.path(nz_enc_dir, "observations.rds"))
saveRDS(network_toUse, file.path(nz_enc_dir, "network.rds"))

Network_sz <- network_toUse %>% select(-c("long", "lat"))


AIC_list <- NULL

##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs_toUse$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs_toUse$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = obs_toUse$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = obs_toUse$dist_i, 
               "Lat" = obs_toUse$lat, 
               "Lon" = obs_toUse$long, 
               "Pass" = 0,
               "Category" = "Longfin_eels")


## include latitude and longitude for user-supplied area
Extrapolation_List = StreamUtils::make_extrapolation_info( Region="Stream", 
  stream_info=cbind("Lat"=obs_toUse$lat, 
                    "Lon"=obs_toUse$long,
                    "child_i"=obs_toUse$child_i,
                    "Area_km2"=obs_toUse$dist_i), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = StreamUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          Lat_x=network_toUse$lat, 
                          Lon_x=network_toUse$long, 
                          Extrapolation_List=Extrapolation_List )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )


########################
## TEMPORAL GAMMA
########################
temp_dir <- file.path(nz_enc_dir, "temporal_gamma")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

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
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=FALSE, savedir=nz_enc_dir, cex=0.2)


TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$temporal_gamma <- as.numeric(Opt$AIC)

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 )

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, cex_network=0.0001)
# Dens_xt <- StreamUtils::plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.00001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year")

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )

###########################
### TEMPORAL LOGNORMAL
###########################
temp_dir <- file.path(nz_enc_dir, "temporal_lognormal")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(1,0)

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
                  # "X_xtp"=X_xtp_inp,
                  # "Xconfig_zcp"=Xconfig_zcp,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$temporal_lognormal <- as.numeric(Opt$AIC)

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 )

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, cex_network=0.0001)
# Dens_xt <- StreamUtils::plot_maps(plot_set=c(3), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.00001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year")

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )

##################
## SPATIAL GAMMA
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_gamma")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

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
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatial_temporal_gamma <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 )

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- StreamUtils::plot_maps(plot_set=3, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.000001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001)

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )

##################
## SPATIAL LOGNORMAL
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_lognormal")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(1,0)

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
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatial_temporal_lognormal <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 )

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- StreamUtils::plot_maps(plot_set=3, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.000001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001)

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )


#####################
## spatiotemporal
st_dir <- file.path(nz_enc_dir, "spatiotemporal_spatial_lognormal")
dir.create(st_dir)
setwd(st_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(1,0)

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
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatiotemporal_spatial_lognormal <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## diagnostics for encounter probability component
Enc_prob = StreamUtils::plot_encounter_diagnostic( Report=Report, Data=Data)

## diagnostics for positive catch rate component
Q = StreamUtils::plot_quantile_diagnostic( TmbData=Data, Report=Report, FileName_QQ="Q-Q_plot", plot=2) #StreamUtils::

## Plot Pearson residuals
StreamUtils::plot_residuals(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, TmbData=Data, Data_Geostat=Data_Geostat, Report=Report, Q=Q, plot_type=1 )

# ## index of abundance
# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

##### Model output
## density surface for each year
Dens_xt <- StreamUtils::plot_maps(plot_set=3, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.000001)
# Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001)

# ## index of abundance
Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )


models <- c("temporal_gamma", "temporal_lognormal", "spatial_gamma", "spatial_lognormal", "spatiotemporal_gamma", "spatiotemporal_spatial_gamma")
AIC <- lapply(1:length(models), function(x){
  res <- readRDS(file.path(nz_enc_dir, models[x], "Save.rds"))
  aic <- as.numeric(res$Opt$AIC)
  if(length(aic)==0) aic <- as.numeric(res$Opt$opt$AIC)
  if(length(aic)>0) df <- data.frame("model"=models[x], "AIC"=aic)
  return(df)
})
AIC <- do.call(rbind, AIC)
AIC <- AIC %>% mutate(dAIC = AIC - min(AIC))