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

devtools::install_github("james-thorson/VAST", ref="master")
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

data <- data("nz_waitaki_longfin_eel_downstream", package="StreamUtils")

network <- nz_waitaki_longfin_eel_downstream[["network"]]
obs <- nz_waitaki_longfin_eel_downstream[["observations"]] %>%
    dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value)
# hab <- nz_waitaki_longfin_eel_downstream[["habitat"]]

Network_sz = network %>% select(c('parent_s','child_s','dist_s'))

p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.2, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  xlab("Longitude") + ylab("Latitude") + 
  mytheme()
# ggsave(file.path(fig_dir, "Waitaki.png"), p1)

p2 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.01, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat, color = factor(present)), pch=19, cex=1.5,alpha=0.8) +
  scale_color_brewer(palette = "Set1") +
  xlab("Longitude") + ylab("Latitude") + 
  guides(color=guide_legend(title="Encounter")) +
  facet_wrap(~year) +
  mytheme()
# ggsave(file.path(fig_dir, "Waitaki_byYear.png"), p2)


####################################
## habitat information
## treated as density covariates
###################################

# nodes <- network$child_s[order(network$child_s)]
# years <- min(obs$year):max(obs$year)
# covar <- unique(hab$covariate)
# n_x <- length(nodes)
# n_t <- length(years)
# n_p <- length(covar)

# df <- lapply(1:length(covar), function(x){
#   sub <- hab %>% filter(covariate == covar[x])
#   df <- data.frame(sub$value)
#   return(df)
# })
# df <- do.call(cbind, df)
# colnames(df) <- covar

# df2 <- cor(df)
# panel.cor <- function(x, y, ...){
#   par(usr = c(0, 1, 0, 1))
#   txt <- as.character(format(cor(x, y), digits=2))
#   text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
# }
# png(file.path(fig_dir, "Habitat_pairs.png"), width=15, height=10, units='in', res=200)
# pairs(df, upper.panel=panel.cor)
# dev.off()



# hab_toUse <- hab_all

# df2 = cor(df)
# fc <- findCorrelation(df2, cutoff=0.5) # 
# fc <- sort(fc)
# reduced_Data = df[,-c(fc)]
# pairs(reduced_Data, upper.panel=panel.cor)

# covar_toUse <- colnames(df)
# # hab_toUse <- hab %>% filter(covariate %in% covar_toUse)

# X_xtp <- array(0, dim=c(n_x, n_t, n_p))
# for(p in 1:n_p){
#   psub <- hab_toUse %>% filter(covariate == covar_toUse[p])
#   mat <- matrix(0, nrow=n_x, ncol = 1)
#   mat[psub$child_s,1] <- psub$value
#   mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
#   X_xtp[,,p] <- mat_sd
# }


# hab_toPlot <- inner_join(network_toUse, hab_toUse %>% select(-'parent_s'))

# for(i in 1:length(covar_toUse)){
#   p <- ggplot(hab_toPlot %>% filter(covariate == covar_toUse[i])) +
#   geom_point(aes(x = long, y = lat, color = value)) +
#   guides(color=guide_legend(title=covar_toUse[i])) +
#   scale_color_viridis_c() +
#   mytheme()
#   ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar_toUse[i],".png")),p)
# }
##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters_DownstreamOnly")
dir.create(nz_enc_dir, showWarnings=FALSE)
setwd(nz_enc_dir)

## single vessel 
vname_i <- unique(obs$fishmethod)
v_i <- sapply(1:nrow(obs), function(x) which(vname_i == obs$fishmethod[x]))

network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001


saveRDS(obs, file.path(nz_enc_dir, "observations.rds"))
saveRDS(network, file.path(nz_enc_dir, "network.rds"))
# saveRDS(hab_toUse, file.path(nz_enc_dir, "habitat.rds"))

AIC_list <- NULL

##### General settings
Version = "VAST_v7_0_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
vname_i <- unique(obs$fishmethod)
v_i <- sapply(1:nrow(obs), function(x) which(vname_i == obs$fishmethod[x]))
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = obs$year,
               "Vessel" = v_i, 
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
temp_dir <- file.path(nz_enc_dir, "temporal_lognormal_deltalink")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

### outputs data input for VAST
## and prints structure for each model parameter
## information on model parameters found at ?VAST::Param
Data = make_data("Version"=Version,
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
                  "a_xl"=matrix(Spatial_List$a_xl[,1], nrow=n_x),
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

# plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=FALSE, savedir=fig_dir, cex=0.2)


TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=FALSE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$temporal_lognormal_deltalink <- as.numeric(Opt$AIC)

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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))

########################
## SPATIAL + TEMPORAL LOGNORMAL DELTA LINK
########################
temp_dir <- file.path(nz_enc_dir, "spatial_temporal_lognormal_deltalink")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=1, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

Data = make_data("Version"=Version,
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
                  "a_xl"=matrix(Spatial_List$a_xl[,1], nrow=n_x),
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

# plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=FALSE, savedir=fig_dir, cex=0.2)


TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$spatial_temporal_lognormal_deltalink <- as.numeric(Opt$AIC)

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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))




#########################
## SPATIOTEMPORAL LOGNORMAL DELTA LINK
########################
st_dir <- file.path(nz_enc_dir, "spatiotemporal_spatial_lognormal_deltalink")
dir.create(st_dir)
setwd(st_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=2, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

Data = make_data("Version"=Version,
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

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatiotemporal_spatial_lognormal_deltalink <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))



##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters_DownstreamOnly_EF")
dir.create(nz_enc_dir, showWarnings=FALSE)
setwd(nz_enc_dir)

## single vessel 
vname_i <- unique(obs$fishmethod)
v_i <- sapply(1:nrow(obs), function(x) which(vname_i == obs$fishmethod[x]))
obs_ef <- obs %>% filter(fishmethod == "Electric fishing")

network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001


saveRDS(obs_ef, file.path(nz_enc_dir, "observations.rds"))
saveRDS(network, file.path(nz_enc_dir, "network.rds"))
# saveRDS(hab_toUse, file.path(nz_enc_dir, "habitat.rds"))

##### General settings
Version = "VAST_v7_0_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs_ef$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs_ef$present <- present_new

##### setup data frame
vname_i <- unique(obs$fishmethod)
v_i <- sapply(1:nrow(obs), function(x) which(vname_i == obs$fishmethod[x]))
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = obs_ef$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = obs_ef$dist_i, 
               "Lat" = obs_ef$lat, 
               "Lon" = obs_ef$long, 
               "Pass" = 0,
               "Category" = "Longfin_eels")


## include latitude and longitude for user-supplied area
Extrapolation_List = StreamUtils::make_extrapolation_info( Region="Stream", 
  stream_info=cbind("Lat"=obs_ef$lat, 
                    "Lon"=obs_ef$long,
                    "child_i"=obs_ef$child_i,
                    "Area_km2"=obs_ef$dist_i), 
                  strata.limits=strata.limits )

## change latitude and longitude by node, not using Kmeans
Spatial_List = StreamUtils::make_spatial_info( n_x=n_x, 
                          Method=Method, 
                          Network_sz=Network_sz,
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          Lat_x=network$lat, 
                          Lon_x=network$long, 
                          Extrapolation_List=Extrapolation_List )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )


########################
## TEMPORAL GAMMA DELTA LINK
########################
temp_dir <- file.path(nz_enc_dir, "temporal_gamma_deltalink")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

Data = make_data("Version"=Version,
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
                  "a_xl"=matrix(Spatial_List$a_xl[,1],nrow=n_x),
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=FALSE, savedir=fig_dir, cex=0.2)

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$temporal_gamma_deltalink_ef <- as.numeric(Opt$AIC)

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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))


########################
## TEMPORAL LOGNORMAL DELTA LINK
########################
temp_dir <- file.path(nz_enc_dir, "temporal_lognormal_deltalink")
dir.create(temp_dir, showWarnings=FALSE)
setwd(temp_dir)

FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

Data = make_data("Version"=Version,
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

# plot_network(Spatial_List=Spatial_List, Extrapolation_List=Extrapolation_List, TmbData=Data, Data_Geostat=Data_Geostat, observations=TRUE, arrows=TRUE, root=FALSE, savedir=fig_dir, cex=0.2)


TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=0, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

## converges without estimating standard error - turn on and estimate
Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report = Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
saveRDS(Save, file.path(getwd(),"Save.rds"))

## convergence
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt[["SD"]]
AIC_list$temporal_lognormal_deltalink_ef <- as.numeric(Opt$AIC)

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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))

##################
## SPATIAL GAMMA DELTA LINK
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_gamma_deltalink")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

Data = make_data("Version"=Version,
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

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatial_temporal_gamma_deltalink_ef <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))

##################
## SPATIAL LOGNORMAL DELTA LINK
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_lognormal_deltalink")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

Data = make_data("Version"=Version,
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

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatial_temporal_lognormal_deltalink_ef <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))

#########################
## SPATIOTEMPORAL GAMMA DELTA LINK
########################
st_dir <- file.path(nz_enc_dir, "spatiotemporal_spatial_gamma_deltalink")
dir.create(st_dir)
setwd(st_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=2, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

Data = make_data("Version"=Version,
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

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatiotemporal_spatial_gamma_deltalink_ef <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))



#########################
## SPATIOTEMPORAL GAMMA DELTA LINK
########################
st_dir <- file.path(nz_enc_dir, "spatiotemporal_spatial_lognormal_deltalink")
dir.create(st_dir)
setwd(st_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=2, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(1,0)

Data = make_data("Version"=Version,
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

TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList[["Obj"]]
Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))

Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
check <- TMBhelper::Check_Identifiable(Obj)

Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

#############
Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

Opt[["SD"]]
AIC_list$spatiotemporal_spatial_lognormal_deltalink_ef <- as.numeric(Opt$AIC)

Report <- Obj$report()
Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
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

Range <- FishStatsUtils::plot_range_index(Sdreport=Opt$SD, Report=Report, TmbData=Data, Year_Set=Year_Set, PlotDir=paste0(getwd(),"/"))




















##################
## SPATIAL LOGNORMAL HABITAT
##################
space_dir <- file.path(nz_enc_dir, "spatial_temporal_lognormal_habitat")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(1,0)

for(p in 1:n_p){
  new_dir <- file.path(space_dir, covar[p])
  dir.create(new_dir, showWarnings=FALSE)
  setwd(new_dir)

  files <- list.files(space_dir)[grep("VAST",list.files(space_dir))]
  ignore <- sapply(1:length(files), function(x) file.copy(from=file.path(space_dir, files[x]), to=new_dir))

  # Xconfig_zcp_inp <- array(0, dim=c(2,1,dim(X_xtp)[3]))
  # Xconfig_zcp_inp[,,p] <- 1
  Data = make_data("Version"=Version,
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
                    "X_xtp" = array(X_xtp[,,p], dim=c(n_x,n_t,1)), 
                    # "Xconfig_zcp"=Xconfig_zcp_inp,
                    "MeshList"=Spatial_List$MeshList,
                    "GridList"=Spatial_List$GridList,
                    "Method"=Spatial_List$Method,
                    "Options"=Options,
                    "Network_sz"=Network_sz ) 

  TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))  

  Opt1 = TMBhelper::Optimize( obj=Obj, starpar=Obj$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  check <- TMBhelper::Check_Identifiable(Obj) 

  Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )  

  #############
  Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] 

  Opt[["SD"]]

  name <- paste0("spatial_temporal_lognormal_", covar[p])
  AIC_list[[name]] <- as.numeric(Opt$AIC) 

  Report <- Obj$report()
  Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
  saveRDS(Save, file.path(getwd(),"Save.rds"))  

  # Save <- readRDS(file.path(getwd(), "Save.rds"))
  # TmbList <- Save$TmbList
  # Data <- Save$Data
  # Obj <- Save$Obj
  # Opt <- Save$Opt
  # Report <- Save$Report 

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
  Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001) 

  # ## index of abundance
  Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" ) 

}




#####################
## spatiotemporal
st_dir <- file.path(nz_enc_dir, "spatiotemporal_spatial_lognormal_habitat")
dir.create(st_dir)
setwd(st_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=2, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(1,0)

for(p in 1:n_p){
  new_dir <- file.path(st_dir, covar[p])
  dir.create(new_dir, showWarnings=FALSE)
  setwd(new_dir)

  files <- list.files(st_dir)[grep("VAST",list.files(st_dir))]
  ignore <- sapply(1:length(files), function(x) file.copy(from=file.path(st_dir, files[x]), to=new_dir))

  # Xconfig_zcp_inp <- array(0, dim=c(2,1,dim(X_xtp)[3]))
  # Xconfig_zcp_inp[,,p] <- 1
  Data = make_data("Version"=Version,
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
                    "X_xtp" = array(X_xtp[,,p], dim=c(n_x,n_t,1)), 
                    # "Xconfig_zcp"=Xconfig_zcp_inp,
                    "MeshList"=Spatial_List$MeshList,
                    "GridList"=Spatial_List$GridList,
                    "Method"=Spatial_List$Method,
                    "Options"=Options,
                    "Network_sz"=Network_sz ) 

  TmbList = make_model("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s']))  

  Opt1 = TMBhelper::Optimize( obj=Obj, starpar=Obj$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  check <- TMBhelper::Check_Identifiable(Obj) 

  Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )  

  #############
  Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')] 

  Opt[["SD"]]

  name <- paste0("spatiotemporal_spatial_lognormal_", covar[p])
  AIC_list[[name]] <- as.numeric(Opt$AIC) 

  Report <- Obj$report()
  Save <- list("TmbList"=TmbList, "Data"=Data, "Obj"=Obj, "Opt"=Opt, "Report"=Report)
  saveRDS(Save, file.path(getwd(),"Save.rds"))  

  # Save <- readRDS(file.path(getwd(), "Save.rds"))
  # TmbList <- Save$TmbList
  # Data <- Save$Data
  # Obj <- Save$Obj
  # Opt <- Save$Opt
  # Report <- Save$Report 

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
  Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001) 

  # ## index of abundance
  Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" ) 

}



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