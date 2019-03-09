rm(list=ls())

################
## Directories
################

nz_dir <- "C:\\merrill\\stream_networks\\NZ\\waikato"
dir.create(nz_dir, showWarnings=FALSE)

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

########################
## read in data
##########################

data <- data("nz_waikato_longfin_eel_downstream", package="StreamUtils")

network <- nz_waikato_longfin_eel_downstream[["network"]] 

Network_sz = network %>% select(c('parent_s','child_s','dist_s'))

obs <- nz_waikato_longfin_eel_downstream[["observations"]] 

# obs_sub <- obs %>% filter(is.na(dist_i))
# up <- obs %>% filter(parent_i == obs_sub$child_i[2])
# down <- obs %>% filter(child_i == obs_sub$parent_i[2])


    # dplyr::filter(data_type=="encounter") %>%
    # select(-data_type) %>%
    # rename('present' = data_value)

hab <- nz_waikato_longfin_eel_downstream[["habitat"]]
hab_toUse <- hab %>% filter(covariate %in% c("upElev", "us_tmin") == FALSE)
# dam <- hab %>% filter(covariate == "DamAffected")

# ## NAs
# dam_gap <- dam %>% filter(is.na(value))

# ## values
# dam_good <- dam %>% filter(is.na(value)==FALSE)

# all(dam_gap$parent_s %in% dam$child_s)

# ## check if segment upstream was affected
# dam_fill <- dam %>% filter(parent_s %in% dam_gap$child_s)
# dam_check <- lapply(1:nrow(dam_gap), function(x){
# 	sub <- dam_gap[x,]
# 	repl <- dam_fill %>% filter(parent_s == sub$child_s)
# 	if(nrow(repl)>0) sub$value <- repl$value
# 	if(nrow(repl)==0) sub$value <- NA
# 	return(sub)
# })
# dam_check <- do.call(rbind, dam_check)

# dam_gap2 <- dam_check %>% filter(is.na(value))
# dam_good2 <- dam_check %>% filter(is.na(value)==FALSE)

# dam_good <- rbind.data.frame(dam_good, dam_good2)

# check1 <- dam_gap2[1,]
# prev <- check1
# save <- prev
# for(i in 1:10){
# 	up <- dam %>% filter(parent_s == prev$child_s)
# 	down <- dam %>% filter(child_s == prev$parent_s)

# 	save <- rbind.data.frame(save, up, down)
# 	if(nrow(up)>1) prev <- up
# 	if(nrow(up)==0) prev <- down
# 	if(any(is.na(save))==FALSE) break
# }
# check1$value <- 1

# check2 <- dam_gap2[2,]
# prev <- check2
# save <- prev
# for(i in 1:10){
# 	up <- dam %>% filter(parent_s == prev$child_s)
# 	down <- dam %>% filter(child_s == prev$parent_s)

# 	save <- rbind.data.frame(save, up, down)
# 	if(nrow(up)>1) prev <- up
# 	if(nrow(up)==0) prev <- down
# 	if(any(is.na(save))==FALSE) break
# }
# check2$value <- 1

# check3 <- dam_gap2[3,]
# prev <- check3
# save <- prev
# for(i in 1:10){
# 	up <- dam %>% filter(parent_s == prev$child_s)
# 	down <- dam %>% filter(child_s == prev$parent_s)

# 	save <- rbind.data.frame(save, up, down)
# 	if(nrow(up)>1) prev <- up
# 	if(nrow(up)==0) prev <- down
# 	if(any(is.na(save))==FALSE) break
# }
# check3$value <- 0

# dam_good <- rbind.data.frame(dam_good, check1, check2, check3)

# hab_rmdam <- hab %>% filter(covariate != "DamAffected")
# hab_wdam <- dam_good

# hab_toUse <- rbind.data.frame(hab_rmdam, hab_wdam)
####################################
## habitat information
## treated as density covariates
###################################

nodes <- network$child_s[order(network$child_s)]
years <- min(obs$year):max(obs$year)
covar <- unique(hab_toUse$covariate)
n_x <- length(nodes)
n_t <- length(years)
n_p <- length(covar)

df <- lapply(1:length(covar), function(x){
  sub <- hab_toUse %>% filter(covariate == covar[x])
  df <- data.frame(sub$value)
  return(df)
})
df <- do.call(cbind, df)
colnames(df) <- covar

df2 <- cor(df)
panel.cor <- function(x, y, ...){
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
png(file.path(fig_dir, "Waikato_habitat_pairs.png"), width=15, height=10, units='in', res=200)
pairs(df, upper.panel=panel.cor)
dev.off()

covar_toUse <- colnames(df)
hab_toUse <- hab_toUse %>% filter(covariate %in% covar_toUse)

hab_toPlot <- inner_join(network, hab_toUse %>% select(-'parent_s'))

for(i in 1:length(covar_toUse)){
  p <- ggplot(hab_toPlot %>% filter(covariate == covar_toUse[i])) +
  geom_point(aes(x = long, y = lat, color = value)) +
  guides(color=guide_legend(title=covar_toUse[i])) +
  scale_color_viridis_c() +
  mytheme()
  ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar_toUse[i],".png")),p)
}

X_xtp <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab_toUse %>% filter(covariate == covar_toUse[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_s,1] <- psub$value
  if(covar_toUse[p]=="DamAffected"){
    X_xtp[,,p] <- mat
  } else {
      mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
      X_xtp[,,p] <- mat_sd
  }
}

## years since dam impact
X_choose <- X_xtp[,,which(covar_toUse == "DamAffected")]
X_xtp1 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- years[x] - 1948
  return(sub)
})
X_xtp1_sd <- (X_xtp1 - mean(X_xtp1))/sd(X_xtp1)

## years since impact squared
X_xtp2 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- (years[x] - 1948)^2
  return(sub)
})
X_xtp2_sd <- (X_xtp2 - mean(X_xtp2))/sd(X_xtp2)

covar_toUse <- covar_toUse[-which(covar_toUse=="DamAffected")]
covar_toUse <- c(covar_toUse, "YearsSinceDam","YearsSinceDam2")
X_xtp_new <- array(0, dim=c(n_x,n_t,n_p+1))
for(p in 1:(n_p+1)){
  if(p < n_p) X_xtp_new[,,p] <- X_xtp[,,p]
  if(p == n_p) X_xtp_new[,,p] <- X_xtp1_sd
  if(p > n_p) X_xtp_new[,,p] <- X_xtp2_sd
}

##################################
## model - count & encounter observations
##################################
res_dir <- file.path(nz_dir, "Counts_Encounters_DownstreamOnly_VesselEffect")
dir.create(res_dir, showWarnings=FALSE)
setwd(res_dir)

saveRDS(obs, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab_toUse, file.path(res_dir, "habitat.rds"))

##### General settings
Version = "VAST_v7_0_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = "Stream_network"
grid_size_km = 1
n_x = nrow(network)   # Specify number of stations (a.k.a. "knots")
strata.limits <- data.frame('STRATA'="All_areas")

##### add small value to encounter observations
present <- obs %>% filter(data_type=="encounter") %>% select(data_value)
devs <- rnorm(nrow(present), 0, 0.01)
present_new <- sapply(1:nrow(present), function(x) ifelse(present[x,1]==1, present[x,1]+devs[x], present[x,1]))
obs[which(obs$data_type == "encounter"),"data_value"] <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = obs$data_value, 
              "Year" = as.numeric(obs$year),
               "Vessel" = obs$fishmethod, 
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Type"=sapply(as.character(obs[,'data_type']),FUN=switch,"encounter"=1,"count"=2),
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

# ######## Modify data
# TypeSet = c(1,2,3)
# Data_Geostat = Data_Geostat[ which(Data_Geostat[,'Type']%in%TypeSet), ]

# #### Relabel Type in same order as TypeSet
# Data_Geostat[,'Type'] = match(Data_Geostat[,'Type'], TypeSet)

AIC_list <- NULL

######## Define whether to include spatial and spatio-temporal variation, the rank of the covariance among species, 
######## whether its autocorrelated, and whether there is overdispersion
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) 
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel_ez = cbind( c(2,14), 1 )

######## Decide on which post-hoc calculations to include in the output
Options =  c("SD_site_density"=0, "SD_site_logdensity"=1, "Calculate_Range"=1, "Calculate_evenness"=0, 
	"Calculate_effective_area"=1, "Calculate_Cov_SE"=0, 'Calculate_Synchrony'=0, 'Calculate_Coherence'=0, 
	'normalize_GMRF_in_CPP'=TRUE)

Data = make_data("Version"=Version,
                  "FieldConfig"=FieldConfig,
                  "OverdispersionConfig"=OverdispersionConfig,
                  "RhoConfig"=RhoConfig,
                  "ObsModel_ez"=ObsModel_ez,
                  "c_iz"=rep(0,nrow(Data_Geostat)),
                  "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'],
                  "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1,
                  "s_i"=Data_Geostat[,'knot_i']-1,
                  "e_i"=Data_Geostat[,'Type']-1, 
                  "t_iz"=Data_Geostat[,'Year'],
                  "a_xl"=matrix(Spatial_List$a_xl[,1],nrow=n_x),
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
AIC_list$temporal_gamma_deltalink <- as.numeric(Opt$AIC)

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
