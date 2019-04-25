rm(list=ls())

################
## Directories
################

nz_dir <- "C:\\merrill\\stream_networks\\NZ"

res_dir <- file.path(nz_dir, "Waitaki")
dir.create(res_dir, showWarnings=FALSE)

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
##remember to pull upstream development branches

devtools::install_github("james-thorson/VAST", ref="development")
devtools::install_github("merrillrudd/FishStatsUtils", ref='stream')

library(VAST)
library(FishStatsUtils)
# library(StreamUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

#########################
## read in data
##########################

data <- data("nz_waitaki_longfin_eel_downstream", package="FishStatsUtils")

network <- nz_waitaki_longfin_eel_downstream[["network"]]

## add a small value to the latitude and longitude where we don't have updated location for the nodes emptying into the ocean
network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- nz_waitaki_longfin_eel_downstream[["observations"]] %>%
    dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value)

## habitat data
hab <- nz_waitaki_longfin_eel_downstream[["habitat"]]
hab <- hab %>% filter(covariate %in% c("upElev", "us_tmin") == FALSE)


###################################
## plot network and observations
###################################


## plot network and observations
p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.5, alpha=0.5) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  xlab("Longitude") + ylab("Latitude") +
  mytheme()
# ggsave(file.path(fig_dir, "Waitaki.png"), p1)

## plot encounters and non-encounters on network by year
p2 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.001, alpha=0.2) +
  geom_point(data = obs, aes(x = long, y = lat, fill = factor(present)), pch=21, cex=1.6,alpha=0.5) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = round(seq(min(network$long),max(network$long), length.out=3),0)) +
  xlab("Longitude") + ylab("Latitude") + 
  guides(fill=guide_legend(title="Encounter")) +
  facet_wrap(~year) +
  mytheme()
# ggsave(file.path(fig_dir, "Waitaki_byYear.png"), p2)


####################################
## habitat information
## treated as density covariates
###################################

nodes <- network$child_s[order(network$child_s)]
years <- min(obs$year):max(obs$year)
covar <- unique(hab$covariate)
n_x <- length(nodes)
n_t <- length(years)
n_p <- length(covar)

df <- lapply(1:length(covar), function(x){
  sub <- hab %>% filter(covariate == covar[x])
  df <- data.frame(sub$value)
  return(df)
})
df <- do.call(cbind, df)
colnames(df) <- covar


# png(file.path(fig_dir, "Habitat_pairs.png"), width=15, height=10, units='in', res=200)
# pairs(df, upper.panel=panel.cor)
# dev.off()


hab_toUse <- hab %>% filter(covariate %in% covar)

# hab_toPlot <- inner_join(network, hab_toUse %>% select(-'parent_s'))

# for(i in 1:length(covar_toUse)){
#   p <- ggplot(hab_toPlot %>% filter(covariate == covar_toUse[i])) +
#   geom_point(aes(x = long, y = lat, color = value)) +
#   guides(color=guide_legend(title=covar_toUse[i])) +
#   scale_color_viridis_c() +
#   mytheme()
#   ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar_toUse[i],".png")),p)
# }

covar1 <- covar
X_gtp_1 <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab_toUse %>% filter(covariate == covar[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_s,1] <- psub$value
  if(covar[p]=="DamAffected"){
    X_gtp_1[,,p] <- mat
  } else {
      mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
      X_gtp_1[,,p] <- mat_sd
  }
}

## years since dam impact
X_choose <- X_gtp_1[,,which(covar == "DamAffected")]
X_gtp1 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- years[x] - 1935
  return(sub)
})
X_gtp1_sd <- (X_gtp1 - mean(X_gtp1))/sd(X_gtp1)

## years since impact squared
X_gtp2 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- (years[x] - 1935)^2
  return(sub)
})
X_gtp2_sd <- (X_gtp2 - mean(X_gtp2))/sd(X_gtp2)

covar2 <- covar[-which(covar=="DamAffected")]
covar2 <- c(covar2, "YearsSinceDam")#,"YearsSinceDam2")
X_gtp_2 <- array(0, dim=c(n_x,n_t,n_p))#+1))
for(p in 1:(n_p)){#+1)){
  if(p < n_p) X_gtp_2[,,p] <- X_gtp_1[,,p]
  if(p == n_p) X_gtp_2[,,p] <- X_gtp1_sd
  # if(p > n_p) X_gtp_2[,,p] <- X_gtp2_sd
}

X_gtp_2_2000 <- X_gtp_2[,which(years >= 2000),]
X_gtp_2_1980 <- X_gtp_2[,which(years >= 1980),]

Xconfig_zcp_inp <- array(0, dim=c(2,1,n_p))
Xconfig_zcp_inp[1,,] <- 1 

##################################
## save data used for model runs
##################################

saveRDS(obs, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab_toUse, file.path(res_dir, "habitat.rds"))


##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = as.numeric(obs$year),
               "Vessel" = obs$fishmethod, 
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = "Longfin_eels")


####### temporal model settings
## turn off spatial and spatio-temporal variation
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
## temporal intercept on encounters independent among years, constant intercept for positive catch rate
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
## one random effect estimated for each unique level of Data_Geostat$Vessel for component 1
OverdispersionConfig = c("Eta1"=1, "Eta2"=0)
## more options
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
## lognormal probability density function, conventional delta-link function
ObsModel = c("PosDist"=1, "Link"=0)

settings_temp <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings_temp$Method <- "Stream_network"
settings_temp$grid_size_km <- 1

######## spatial model settings
FieldConfig = c("Omega1"="IID", "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=1, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c("PosDist"=1, "Link"=0)

settings_st <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings_st$Method <- "Stream_network"
settings_st$grid_size_km <- 1

####### spatiotemporal model settings
FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=3, "Epsilon1"=2, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=1, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c("PosDist"=1, "Link"=0)

settings_sst <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings_sst$Method <- "Stream_network"
settings_sst$grid_size_km <- 1

#############################
## directories
#############################

temp_dir <- file.path(res_dir, "Temporal")
dir.create(temp_dir, showWarnings=FALSE)
tmfig_dir <- file.path(temp_dir, "figures")
dir.create(tmfig_dir, showWarnings=FALSE)

st_dir <- file.path(res_dir, "Spatial+Temporal")
dir.create(st_dir, showWarnings=FALSE)
stfig_dir <- file.path(st_dir, "figures")
dir.create(stfig_dir, showWarnings=FALSE)

sst_dir <- file.path(res_dir, "Spatiotemporal+Spatial+Temporal")
dir.create(sst_dir, showWarnings=FALSE)
sstfig_dir <- file.path(sst_dir, "figures")
dir.create(sstfig_dir, showWarnings=FALSE)

hst_dir <- file.path(res_dir, "Habitat+Spatial+Temporal")
dir.create(hst_dir, showWarnings=FALSE)
hstfig_dir <- file.path(hst_dir, "figures")
dir.create(hstfig_dir, showWarnings=FALSE)

hsst_dir <- file.path(res_dir, "Habitat+Spatiotemporal+Spatial+Temporal")
dir.create(hsst_dir, showWarnings=FALSE)
hsstfig_dir <- file.path(hsst_dir, "figures")
dir.create(hsstfig_dir, showWarnings=FALSE)

hsst_dir1 <- file.path(res_dir, "Habitat+Spatiotemporal+Spatial+Temporal_1p")
dir.create(hsst_dir1, showWarnings=FALSE)
hsstfig_dir1 <- file.path(hsst_dir1, "figures")
dir.create(hsstfig_dir1, showWarnings=FALSE)

sst2000_dir <- file.path(res_dir, "Spatiotemporal+Spatial+Temporal_2000")
dir.create(sst2000_dir, showWarnings=FALSE)
sst2000fig_dir <- file.path(sst2000_dir, "figures")
dir.create(sst2000fig_dir, showWarnings=FALSE)

hsst2000_dir <- file.path(res_dir, "Habitat+Spatiotemporal+Spatial+Temporal_2000")
dir.create(hsst2000_dir, showWarnings=FALSE)
hsst2000fig_dir <- file.path(hsst2000_dir, "figures")
dir.create(hsst2000fig_dir, showWarnings=FALSE)

sst1980_dir <- file.path(res_dir, "Spatiotemporal+Spatial+Temporal_1980")
dir.create(sst1980_dir, showWarnings=FALSE)
sst1980fig_dir <- file.path(sst1980_dir, "figures")
dir.create(sst1980fig_dir, showWarnings=FALSE)

hsst1980_dir <- file.path(res_dir, "Habitat+Spatiotemporal+Spatial+Temporal_1980")
dir.create(hsst1980_dir, showWarnings=FALSE)
hsst1980fig_dir <- file.path(hsst1980_dir, "figures")
dir.create(hsst1980fig_dir, showWarnings=FALSE)

##################
## Temporal only
##################
setwd(temp_dir)

# ## wrapper function
fit1 = fit_model( "settings"=settings_temp, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=temp_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE) )
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit = fit_model( "settings"=settings_temp, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=temp_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3)
saveRDS(fit, file.path(temp_dir, "Fit.rds"))

fit <- readRDS(file.path(temp_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_temp, fit=fit, working_dir=tmfig_dir, category_names="Longfin eels", strata_names="Longfin eels")

##################
## Spatial + temporal
##################
setwd(st_dir)


## wrapper function
fit1 = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=st_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE) )
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

## wrapper function
fit = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=st_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3)
saveRDS(fit, file.path(st_dir, "Fit.rds"))

fit <- readRDS(file.path(st_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_st, fit=fit, working_dir=stfig_dir, category_names="Longfin eels", strata_names = "Longfin eels" )


##################
## Spatiotemporal + spatial + temporal
##################
setwd(sst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=sst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=sst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3)#,
               # optimize_args = list(obj=Obj))
saveRDS(fit, file.path(sst_dir, "Fit.rds"))

fit <- readRDS(file.path(sst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_sst, fit=fit, working_dir=sstfig_dir,category_names="Longfin eels", strata_names = "Longfin eels" )


##################
## Spatial + temporal + habitat
##################
setwd(hst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE),
                Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_2, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2:3])))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                 Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_2, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2:3])),
                newtonsteps=3)#,
               # optimize_args = list(obj=Obj))
saveRDS(fit, file.path(hst_dir, "Fit.rds"))

fit <- readRDS(file.path(hst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_st, fit=fit, working_dir=hstfig_dir, category_names="Longfin eels", strata_names = "Longfin eels", plot_set=c(3,15), covar_names=covar2 )



##################
## Spatiotemporal + spatial + temporal + habitat
##################
setwd(hsst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE),
                X_gtp = X_gtp_1, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_1)[2:3])))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_1, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_1)[2:3])),
                newtonsteps=3)#,
               # optimize_args = list(obj=Obj))
saveRDS(fit, file.path(hsst_dir, "Fit.rds"))

fit <- readRDS(file.path(hsst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_sst, fit=fit, working_dir=hsstfig_dir, category_names="Longfin eels", strata_names = "Longfin eels", plot_set=c(3,15), covar_names = covar1 )

##################
## Spatiotemporal + spatial + temporal + habitat
##################
setwd(hsst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE),
                X_gtp = X_gtp_2, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2:3])))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_2, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2:3])),
                newtonsteps=3)#,
               # optimize_args = list(obj=Obj))
saveRDS(fit, file.path(hsst_dir, "Fit.rds"))

fit <- readRDS(file.path(hsst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_sst, fit=fit, working_dir=hsstfig_dir, category_names="Longfin eels", strata_names = "Longfin eels", plot_set=c(3,15), covar_names = covar2 )


##################
## Spatiotemporal + spatial + temporal + 1 habitat covariate at a time
##################
setwd(hsst_dir1)

for(p in 1:(n_p-1)){

  subdir <- file.path(hsst_dir1, covar2[p])
  dir.create(subdir, showWarnings=FALSE)
  setwd(subdir)

  # Xconfig_zcp_inp = array(0, dim=c(2,1,n_p))
  # Xconfig_zcp_inp[,,p] <- 1

  ## wrapper function
  fit1 = fit_model( "settings"=settings_sst, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=subdir, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE),
                  X_gtp = array(X_gtp_2[,,p], dim=c(dim(X_gtp_2)[1:2],1)), X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2],1)))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
  # Obj <- fit1$tmb_list$Obj
  # Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))  

  fit = fit_model( "settings"=settings_sst, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=subdir, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=obs$lat, "Lon"=obs$long,"child_i"=obs$child_i,"Area_km2"=obs$dist_i), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                 X_gtp = array(X_gtp_2[,,p], dim=c(dim(X_gtp_2)[1:2],1)), X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_2)[2],1)),
                  newtonsteps=3)#,
                 # optimize_args = list(obj=Obj))
  saveRDS(fit, file.path(subdir, "Fit.rds")) 

  fit <- readRDS(file.path(subdir, "Fit.rds")) 

  # Plot results
  plot_results( settings=settings_sst, fit=fit, working_dir=subdir, category_names="Longfin eels", strata_names = "Longfin eels", plot_set=c(3,15), covar_names = covar2 )

}




############################
## from 2000
############################

# ##################
# ## Spatiotemporal + spatial + temporal
# ##################
# setwd(sst2000_dir)

# Data_Geostat_2000 <- Data_Geostat %>% filter(Year >= 2000)

# ## wrapper function
# fit1 = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_2000[,"Lat"], 
#                 "Lon_i"=Data_Geostat_2000[,"Lon"], 
#                 "t_iz"=Data_Geostat_2000[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_2000)), 
#                 "b_i"=Data_Geostat_2000[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_2000[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_2000[,'Vessel'], 
#                 working_dir=sst2000_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_2000[,"Lat"], "Lon"=Data_Geostat_2000[,"Lon"],"child_i"=Data_Geostat_2000[,"Knot"],"Area_km2"=Data_Geostat_2000[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 optimize_args = list(getsd=FALSE))
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# # Obj <- fit1$tmb_list$Obj
# # Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

# fit = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_2000[,"Lat"], 
#                 "Lon_i"=Data_Geostat_2000[,"Lon"], 
#                 "t_iz"=Data_Geostat_2000[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_2000)), 
#                 "b_i"=Data_Geostat_2000[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_2000[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_2000[,'Vessel'], 
#                 working_dir=sst2000_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_2000[,"Lat"], "Lon"=Data_Geostat_2000[,"Lon"],"child_i"=Data_Geostat_2000[,"Knot"],"Area_km2"=Data_Geostat_2000[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 newtonsteps=3)#,
#                # optimize_args = list(obj=Obj))
# saveRDS(fit, file.path(sst2000_dir, "Fit.rds"))

# fit <- readRDS(file.path(sst2000_dir, "Fit.rds"))

# # Plot results
# plot_results( settings=settings_sst, fit=fit, working_dir=sst2000fig_dir,category_names="Longfin eels", strata_names = "Longfin eels" )


# ##################
# ## Spatiotemporal + spatial + temporal + habitat
# ##################
# setwd(hsst2000_dir)

# Data_Geostat_2000 <- Data_Geostat %>% filter(Year >= 2000)

# ## wrapper function
# fit1 = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_2000[,"Lat"], 
#                 "Lon_i"=Data_Geostat_2000[,"Lon"], 
#                 "t_iz"=Data_Geostat_2000[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_2000)), 
#                 "b_i"=Data_Geostat_2000[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_2000[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_2000[,'Vessel'], 
#                 working_dir=hsst2000_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_2000[,"Lat"], "Lon"=Data_Geostat_2000[,"Lon"],"child_i"=Data_Geostat_2000[,"Knot"],"Area_km2"=Data_Geostat_2000[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 X_gtp = X_gtp_2_2000, X_itp = array(0, dim=c(nrow(Data_Geostat_2000),dim(X_gtp_2_2000)[2:3])),
#                 optimize_args = list(getsd=FALSE))
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# # Obj <- fit1$tmb_list$Obj
# # Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

# fit = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_2000[,"Lat"], 
#                 "Lon_i"=Data_Geostat_2000[,"Lon"], 
#                 "t_iz"=Data_Geostat_2000[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_2000)), 
#                 "b_i"=Data_Geostat_2000[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_2000[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_2000[,'Vessel'], 
#                 working_dir=hsst2000_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_2000[,"Lat"], "Lon"=Data_Geostat_2000[,"Lon"],"child_i"=Data_Geostat_2000[,"Knot"],"Area_km2"=Data_Geostat_2000[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 X_gtp = X_gtp_2_2000, X_itp = array(0, dim=c(nrow(Data_Geostat_2000),dim(X_gtp_2_2000)[2:3])),
#                 newtonsteps=3)
# saveRDS(fit, file.path(hsst2000_dir, "Fit.rds"))

# fit <- readRDS(file.path(hsst2000_dir, "Fit.rds"))

# # Plot results
# plot_results( settings=settings_sst, fit=fit, working_dir=hsst2000fig_dir,category_names="Longfin eels", strata_names = "Longfin eels" )


# ############################
# ## from 1980
# ############################

# ##################
# ## Spatiotemporal + spatial + temporal
# ##################
# setwd(sst1980_dir)

# Data_Geostat_1980 <- Data_Geostat %>% filter(Year >= 1980)

# ## wrapper function
# fit1 = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_1980[,"Lat"], 
#                 "Lon_i"=Data_Geostat_1980[,"Lon"], 
#                 "t_iz"=Data_Geostat_1980[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_1980)), 
#                 "b_i"=Data_Geostat_1980[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_1980[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_1980[,'Vessel'], 
#                 working_dir=sst1980_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_1980[,"Lat"], "Lon"=Data_Geostat_1980[,"Lon"],"child_i"=Data_Geostat_1980[,"Knot"],"Area_km2"=Data_Geostat_1980[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 optimize_args = list(getsd=FALSE))
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# # Obj <- fit1$tmb_list$Obj
# # Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

# fit = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_1980[,"Lat"], 
#                 "Lon_i"=Data_Geostat_1980[,"Lon"], 
#                 "t_iz"=Data_Geostat_1980[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_1980)), 
#                 "b_i"=Data_Geostat_1980[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_1980[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_1980[,'Vessel'], 
#                 working_dir=sst1980_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_1980[,"Lat"], "Lon"=Data_Geostat_1980[,"Lon"],"child_i"=Data_Geostat_1980[,"Knot"],"Area_km2"=Data_Geostat_1980[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 newtonsteps=3)#,
#                # optimize_args = list(obj=Obj))
# saveRDS(fit, file.path(sst1980_dir, "Fit.rds"))

# fit <- readRDS(file.path(sst1980_dir, "Fit.rds"))

# # Plot results
# plot_results( settings=settings_sst, fit=fit, working_dir=sst1980fig_dir, category_names="Longfin eels", strata_names = "Longfin eels" )


# ##################
# ## Spatiotemporal + spatial + temporal + habitat
# ##################
# setwd(hsst1980_dir)

# Data_Geostat_1980 <- Data_Geostat %>% filter(Year >= 1980)

# ## wrapper function
# fit1 = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_1980[,"Lat"], 
#                 "Lon_i"=Data_Geostat_1980[,"Lon"], 
#                 "t_iz"=Data_Geostat_1980[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_1980)), 
#                 "b_i"=Data_Geostat_1980[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_1980[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_1980[,'Vessel'], 
#                 working_dir=hsst1980_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_1980[,"Lat"], "Lon"=Data_Geostat_1980[,"Lon"],"child_i"=Data_Geostat_1980[,"Knot"],"Area_km2"=Data_Geostat_1980[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 X_gtp = X_gtp_2_1980, X_itp = array(0, dim=c(nrow(Data_Geostat_1980),dim(X_gtp_2_1980)[2:3])),
#                 optimize_args = list(getsd=FALSE))
# check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# # Obj <- fit1$tmb_list$Obj
# # Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

# fit = fit_model( "settings"=settings_sst, 
#                 "Lat_i"=Data_Geostat_1980[,"Lat"], 
#                 "Lon_i"=Data_Geostat_1980[,"Lon"], 
#                 "t_iz"=Data_Geostat_1980[,'Year'], 
#                 "c_i"=rep(0,nrow(Data_Geostat_1980)), 
#                 "b_i"=Data_Geostat_1980[,'Catch_KG'], 
#                 "a_i"=Data_Geostat_1980[,'AreaSwept_km2'], 
#                 "v_i"=Data_Geostat_1980[,'Vessel'], 
#                 working_dir=hsst1980_dir, 
#                 extrapolation_args=list(
#                   input_grid=cbind("Lat"=Data_Geostat_1980[,"Lat"], "Lon"=Data_Geostat_1980[,"Lon"],"child_i"=Data_Geostat_1980[,"Knot"],"Area_km2"=Data_Geostat_1980[,"AreaSwept_km2"]), 
#                   Network_sz_LL=Network_sz_LL),
#                 Network_sz = Network_sz,
#                 X_gtp = X_gtp_2_1980, X_itp = array(0, dim=c(nrow(Data_Geostat_1980),dim(X_gtp_2_1980)[2:3])),
#                 newtonsteps=3)
# saveRDS(fit, file.path(hsst1980_dir, "Fit.rds"))

# fit <- readRDS(file.path(hsst1980_dir, "Fit.rds"))

# # Plot results
# plot_results( settings=settings_sst, fit=fit, working_dir=hsst1980fig_dir,category_names="Longfin eels", strata_names = "Longfin eels" )








###############
## AIC
###############
# models <- c("Temporal","Spatial+Temporal","Habitat+Spatial+Temporal","Spatiotemporal+Spatial+Temporal","Habitat+Spatiotemporal+Spatial+Temporal")
# df <- data.frame("Model"=models, "Directory"= c(temp_dir, st_dir, hst_dir2, sst_dir, hsst_dir2))

models <- c("Temporal","Spatial+Temporal","Spatiotemporal+Spatial+Temporal","Habitat+Spatial+Temporal","Habitat+Spatiotemporal+Spatial+Temporal")
df <- data.frame("Model"=models, "Directory"= c(temp_dir, st_dir, sst_dir, hst_dir, hsst_dir))


AIC <- NULL
for(i in 1:length(models)){
  fit <- readRDS(file.path(df[i,"Directory"], "Fit.rds"))
  AIC[[i]] <- fit$parameter_estimates$AIC
}
df$AIC <- unlist(AIC)
df <- df %>% mutate(dAIC = AIC - min(AIC))
