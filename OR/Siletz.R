rm(list=ls())

################
## Directories
################

or_dir <- "C:\\merrill\\stream_networks\\OR"

res_dir <- file.path(or_dir, "Siletz")
dir.create(res_dir, showWarnings=FALSE)

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
##remember to pull upstream development branches

devtools::install_github("merrillrudd/VAST", ref="stream")
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

data <- data("oregon_siletz_coho", package="FishStatsUtils")

network <- or_siletz_coho[["network"]]

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- or_siletz_coho[["observations"]] 
category_names <- unique(obs$survey)

## habitat data
hab <- or_siletz_coho[["habitat"]]


###################################
## plot network and observations
###################################


## plot network and observations
p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=1, alpha=0.5) +
  geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  geom_point(data = network %>% filter(parent_s == 0), aes(x = long, y = lat), pch=19, cex=3, color="goldenrod") +
  xlab("Longitude") + ylab("Latitude") +
  mytheme()
ggsave(file.path(fig_dir, "Siletz.png"), p1)

## plot encounters and non-encounters on network by year
p2 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.001, alpha=0.2) +
  geom_point(data = obs %>% filter(survey == "spawners"), aes(x = long, y = lat, fill = density), pch=21, cex=2, alpha=0.8) +
  geom_point(data = obs %>% filter(survey == "juveniles"), aes(x = long, y = lat, color = density), pch=22, cex=2) +
  # scale_fill_brewer(palette = "Set1") +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = round(seq(min(network$long),max(network$long), length.out=3),0)) +
  xlab("Longitude") + ylab("Latitude") + 
  guides(fill=guide_legend(title="spawners"), color=guide_legend(title="juveniles")) +
  facet_wrap(~year) +
  mytheme()
ggsave(file.path(fig_dir, "Siletz_byYear.png"), p2)


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

# hab_toPlot <- inner_join(obs, hab_toUse )

for(i in 1:length(covar)){
  p <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), alpha = 0.3, cex=0.8) +
  geom_point(data = hab_toUse %>% filter(covariate == covar[i]), aes(x = long, y = lat, fill = value), cex=4, pch=21) +
  guides(fill=guide_legend(title=covar[i])) +
  scale_fill_viridis_c() +
  mytheme()
  ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p)
}

X_gtp_inp <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab_toUse %>% filter(covariate == covar[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_i,1] <- psub$value
      mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
      X_gtp_inp[,,p] <- mat_sd
}

Xconfig_zcp_inp <- array(1, dim=c(2,length(category_names),n_p))

##################################
## save data used for model runs
##################################

saveRDS(obs, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab_toUse, file.path(res_dir, "habitat.rds"))


##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = obs$density, 
              "Year" = as.numeric(obs$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = obs$surveynum)



######## spatial model settings
FieldConfig = c("Omega1"="IID", "Epsilon1"=0, "Omega2"="IID", "Epsilon2"=0)
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c("PosDist"=2, "Link"=1)

settings_st <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, n_categories = length(unique(Data_Geostat$Category)) )
settings_st$Method <- "Stream_network"
settings_st$grid_size_km <- 1

####### spatiotemporal model settings
FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c("PosDist"=2, "Link"=1)

settings_sst <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE, n_categories = length(unique(Data_Geostat$Category)) )
settings_sst$Method <- "Stream_network"
settings_sst$grid_size_km <- 1

#############################
## directories
#############################

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

##################
## Spatial + temporal
##################
setwd(st_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=st_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE) )
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)


## wrapper function
fit = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=st_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3,
                optimize_args = list(startpar = fit1$tmb_list$Obj$par))
saveRDS(fit, file.path(st_dir, "Fit.rds"))

fit <- readRDS(file.path(st_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_st, fit=fit, working_dir=stfig_dir, category_names=category_names, strata_names = "Coho salmon" )



##################
## Spatiotemporal + spatial + temporal
##################
setwd(sst_dir)


## wrapper function
fit1 = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=sst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE) )
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)


## wrapper function
fit = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=sst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3)
saveRDS(fit, file.path(sst_dir, "Fit.rds"))

fit <- readRDS(file.path(sst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_sst, fit=fit, working_dir=sstfig_dir, category_names=category_names, strata_names = "Coho salmon" )


##################
## Spatial + temporal + habitat
##################
setwd(hst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE),
                Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_inp, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_inp)[2:3])))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_st, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_inp, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_inp)[2:3])),
                newtonsteps=3,
               optimize_args = list(startpar = fit1$tmb_list$Obj$par))
saveRDS(fit, file.path(hst_dir, "Fit.rds"))

fit <- readRDS(file.path(hst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_st, fit=fit, working_dir=hstfig_dir, category_names = category_names, strata_names="Coho salmon", plot_set=c(3,15), covar_names=covar )


##################
## Spatiotemporal + spatial + temporal + habitat
##################
setwd(hsst_dir)

## wrapper function
fit1 = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(getsd=FALSE),
                X_gtp = X_gtp_inp, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_inp)[2:3])))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)
# Obj <- fit1$tmb_list$Obj
# Obj$par[grep("logkappa",names(Obj$par))] = 0.3 #log(1/median(Network_sz[,'dist_s']))

fit = fit_model( "settings"=settings_sst, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=Data_Geostat[,"Category"]-1, 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=hsst_dir, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"], "child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = array(0, dim=c(nrow(Data_Geostat),dim(X_gtp_inp)[2:3])),
                newtonsteps=3)#,
               # optimize_args = list(obj=Obj))
saveRDS(fit, file.path(hsst_dir, "Fit.rds"))

fit <- readRDS(file.path(hsst_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_sst, fit=fit, working_dir=hsstfig_dir, category_names = category_names, strata_names="Coho salmon", plot_set=c(3,15), covar_names=covar )




models <- c("Spatial+Temporal","Spatiotemporal+Spatial+Temporal", "Habitat+Spatial+Temporal", "Habitat+Spatiotemporal+Spatial+Temporal")
df <- data.frame("Model"=models, "Directory"= c( st_dir, sst_dir, hst_dir, hsst_dir))


AIC <- NULL
for(i in 1:length(models)){
  fit <- readRDS(file.path(df[i,"Directory"], "Fit.rds"))
  AIC[[i]] <- fit$parameter_estimates$AIC
}
df$AIC <- unlist(AIC)
df <- df %>% mutate(dAIC = AIC - min(AIC))
