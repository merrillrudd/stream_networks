rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks"
R_dir <- file.path(main_dir, "R")
R_files <- list.files(R_dir)
readr <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

nz_dir <- file.path(main_dir, "NZ")

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
library(foreach)
library(doParallel)

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
    rename('present' = data_value) %>%
    mutate('vessel_agency' = paste0(fishmethod, "_", agency)) %>%
    mutate('fishmethod2' = ifelse(fishmethod != "Electric fishing", "Other", "Electric fishing")) %>%
    rename("Year"=year)

  # ### by year
  # years <- unique(obs$Year)[order(unique(obs$Year))]

  # Network_sz_LL_wYear <- lapply(1:length(years), function(x){
  #   out <- cbind.data.frame(Network_sz_LL, "Year"=years[x])
  #   return(out)
  # })
  # Network_sz_LL_wYear <- do.call(rbind, Network_sz_LL_wYear)

  # bb <- ggplot(Network_sz_LL_wYear) +
  #     geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", cex=0.5, alpha=0.6) +
  #     geom_point(data = obs, aes(x = long, y = lat, fill = factor(round(present,0))), cex=1.8, pch=22, alpha=0.6) +
  #     scale_fill_brewer(palette = "Set1") +
  #     scale_x_continuous(breaks = round(quantile(obs$long,prob=c(0.2,0.5,0.8)),0), labels=round(quantile(obs$long,prob=c(0.2,0.5,0.8)),0)) +
  #     xlab("Longtiude") + ylab("Latitude") +
  #     facet_wrap(~Year) +
  #     guides(fill = guide_legend(title = "Encounter")) +
  #     mytheme()
  # ggsave(file.path(fig_dir, "Network_byYear_byEncounter.png"), bb, width = 10, height = 8)

##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = as.numeric(obs$Year),
               "Vessel" = obs$fishmethod,
               "Vessel2" = obs$vessel_agency, 
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = "Longfin_eels")

##### setup data frame
Data_Geostat_ef <- Data_Geostat %>% filter(Vessel == "Electric fishing")
obs_ef <- obs %>% filter(fishmethod == "Electric fishing")

###### remove angling
# Data_Geostat <- Data_Geostat %>% filter(Vessel != "Angling")
# obs <- obs %>% filter(fishmethod != "Angling")

## habitat data
hab <- nz_waitaki_longfin_eel_downstream[['habitat']]

covar_toUse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm",'DamAffected')

hab <- hab %>% filter(covariate %in% covar_toUse)

####################################
## habitat information
## treated as density covariates
###################################

nodes <- network$child_s[order(network$child_s)]
years <- min(obs$Year):max(obs$Year)
covar <- unique(hab$covariate)
n_x <- length(nodes)
n_t <- length(years)
n_p <- length(covar)
n_i <- nrow(obs)
n_i_ef <- nrow(obs_ef)

for(i in 1:length(covar)){
  p <- ggplot(hab %>% filter(covariate == covar[i])) +
  geom_point(aes(x = easting, y = northing, color = value)) +
  guides(color=guide_legend(title=covar[i])) +
  scale_color_viridis_c() +
  mytheme()
  ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p, width=10, height=8)
}

X_gtp_input1 <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab %>% filter(covariate == covar[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_s,1] <- psub$value
  if(covar[p]=="DamAffected"){
    X_gtp_input1[,,p] <- mat
  } else {
      mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
      X_gtp_input1[,,p] <- mat_sd
  }
}

## years since dam impact
X_choose <- X_gtp_input1[,,which(covar == "DamAffected")]
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

covar2 <- c(covar, "YearsSinceDam","YearsSinceDam2")[-which(covar=="DamAffected")]
n_p <- length(covar2)
X_gtp_input <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:(n_p)){
  ## skip dam affected
  if(p < length(covar)) X_gtp_input[,,p] <- X_gtp_input1[,,p]

  ## in place of dam affected, years since dam
  if(p == length(covar)) X_gtp_input[,,p] <- X_gtp1_sd

  ## additional covariate, years since dam squared
  if(p ==length(covar)+1) X_gtp_input[,,p] <- X_gtp2_sd
}

## match habitat covariates to observations
## double check the indices will match up properly
X_itp_input <- array(0, dim=c(n_i,n_t,n_p))
for(i in 1:n_i){
  for(p in 1:n_p){
    child_i <- obs$child_i[i]
    index <- which(nodes == child_i)
    X_itp_input[i,,p] <- X_gtp_input[index,,p]
  }
}

X_itp_ef <- array(0, dim=c(n_i_ef,n_t,n_p))
for(i in 1:n_i_ef){
  for(p in 1:n_p){
    child_i <- obs_ef$child_i[i]
    index <- which(nodes == child_i)
    X_itp_ef[i,,p] <- X_gtp_input[index,,p]
  }
}

# X_itp <- array(0, dim=c(n_i,n_t,n_p))
# for(i in 1:n_i){
#   for(p in 1:n_p){
#     child_i <- obs$child_i[i]
#     index <- which(nodes == child_i)
#     X_itp[i,,p] <- X_gtp_input[index,,p]
#   }
# }

####################################
## sampling information
## treated as catchability covariates
###################################

method <- unique(obs$fishmethod)
cmethod <- sapply(1:length(method), function(x) length(which(obs$fishmethod ==method[x])) / nrow(obs))
names(cmethod) <- method
cmethod[order(cmethod)]

Q_ik_method <- ThorsonUtilities::vector_to_design_matrix(obs $fishmethod )[,-2,drop=FALSE]
Q_ik_method2 <- ThorsonUtilities::vector_to_design_matrix( obs$fishmethod2 )[,-1,drop=FALSE]


agency <- unique(obs$agency)
cagency <- sapply(1:length(agency), function(x) length(which(obs$agency == agency[x])) / nrow(obs))
names(cagency) <- agency
cagency[order(cagency)]
## agency2
obs$agency <- as.character(obs$agency)
obs$agency2 <- sapply(1:nrow(obs), function(x) ifelse(obs$agency[x] %in% c("doc",'university'), obs$agency[x], "other"))
Q_ik_agency2 <- ThorsonUtilities::vector_to_design_matrix( obs$agency2 )[,-1,drop=FALSE]

metag <- unique(obs$vessel_agency)
cmetag <- sapply(1:length(metag), function(x) length(which(obs$vessel_agency == metag[x])) / nrow(obs))
names(cmetag) <- metag
cmetag[order(cmetag)]
obs$vessel_agency2 <- sapply(1:nrow(obs), function(x) ifelse(obs$vessel_agency[x] %in% c("Electric fishing_doc","Electric fishing_university"), obs$vessel_agency[x], "Other"))
Q_ik_vesselagency2 <- ThorsonUtilities::vector_to_design_matrix( obs$vessel_agency2 )[,-1,drop=FALSE]

# vessels <- unique(obs$fishmethod)
# count_vessels <- sapply(1:length(vessels), function(x) length(which(obs$fishmethod == vessels[x])))
# Q_ik_method <- ThorsonUtilities::vector_to_design_matrix( obs$fishmethod )[,-1,drop=FALSE]

# agency <- unique(obs_ef$agency)
# count_agency <- sapply(1:length(agency), function(x) length(which(obs_ef$agency == agency[x])))
# names(count_agency) <- agency

# obs_ef$agency2 <- sapply(1:nrow(obs_ef), function(x) ifelse(obs_ef$agency[x] %in% c('fish&game','council','consultants'), 'other', as.character(obs_ef$agency[x])))
# agency2 <- unique(obs_ef$agency2)
# count_agency2 <- sapply(1:length(agency2), function(x) length(which(obs_ef$agency2 == agency2[x])))
# names(count_agency2) <- agency2

# Q_ik_agency_ef <- ThorsonUtilities::vector_to_design_matrix( obs_ef$agency2 )[,-1,drop=FALSE]

##################################
## save data used for model runs
##################################

saveRDS(obs, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab, file.path(res_dir, "habitat.rds"))


##################################
## Models
##################################
#########################
## SST_IID_ALL_HAB
#########################
path <- file.path(res_dir, "SST_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # save.image(file.path(nz_dir, "test", "NZ_eel_gear.Rdata"))

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map))#,
                # optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    



#########################
## ST_IID_ALL_HAB
#########################
path <- file.path(res_dir, "ST_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.3

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## ST_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "ST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # save.image(file.path(nz_dir, "test", "NZ_eel_gear.Rdata"))

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  ## adjust starting value
  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.3
  # newpar[["L_epsilon1_z"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## SST_RW_ALL_HAB
#########################
path <- file.path(res_dir, "SST_RW_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- -1

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_RW_GEAR_HAB
#########################
path <- file.path(res_dir, "SST_RW_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## SST_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## ST_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "ST_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


#########################
## SST_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## ST_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "ST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## T_IID_ALL_HAB
#########################
path <- file.path(res_dir, "T_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## T_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "T_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## T_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "T_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## T_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "T_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_AGENCY_HAB
#########################
path <- file.path(res_dir, "SST_IID_AGENCY_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(grepl("EF", msubx[3])==FALSE){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_GEARAGENCY_HAB
#########################
path <- file.path(res_dir, "SST_IID_GEARAGENCY_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(grepl("EF", msubx[3])==FALSE){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


df <- data.frame("Model"=c("SST_IID_AGENCY_HAB",
                          "SST_IID_ALL_HAB",
                          "SST_IID_ALL_NOHAB",
                          "SST_IID_GEAR_HAB",
                          "SST_IID_GEAR_NOHAB",
                          "SST_IID_GEARAGENCY_HAB",
                          "SST_RW_ALL_HAB",
                          "SST_RW_GEAR_HAB",
                          "ST_IID_ALL_HAB",
                          "ST_IID_ALL_NOHAB",
                          "ST_IID_GEAR_HAB",
                          "ST_IID_GEAR_NOHAB",
                          "T_IID_ALL_HAB",
                          "T_IID_ALL_NOHAB",
                          "T_IID_GEAR_HAB",
                          "T_IID_GEAR_NOHAB"))
df$AIC <- NA

for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df <- df %>% mutate(dAIC = AIC - min(AIC, na.rm=TRUE))
df[order(df$dAIC),]


gamma <- list()
for(i in 1:nrow(df)){
    res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
    sum <- summary(res$parameter_estimates$SD)
    sum2 <- sum[grepl('gamma1',rownames(sum)),]
    if(nrow(sum2)>0){
      gamma[[i]] <- data.frame("Model"=df[i,"Model"], "gamma"=sum2[,"Estimate"])
    } else{ gamma[[i]] <- NULL}
}


par <- list()
for(i in 1:nrow(df)){
    res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
    dfx <- data.frame(Model = df[i,"Model"], 'par' = res$parameter_estimates$par, "names"=names(res$parameter_estimates$par))
    par[[i]] <- dfx
}








#########################
## SST_IID_EF_HAB
#########################
path <- file.path(res_dir, "SST_IID_EF_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  ## adjust starting value
  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- -0.5
  newpar[["L_epsilon1_z"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## SST_IID_EFAG_HAB
#########################
path <- file.path(res_dir, "SST_IID_EFAG_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_agency_ef
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  ## adjust starting value
  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 1
  newpar[["L_epsilon1_z"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## ST_IID_EF_HAB
#########################
path <- file.path(res_dir, "ST_IID_EF_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # Map <- fit0$tmb_list$Map
  # Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 10
  # newpar[["L_epsilon1_z"]] <- 0.01

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 # model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3)#,
                # model_args = list(Map = Map),
                # optimize_args = list(startpar= fit1$parameter_estimates))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## ST_IID_EFAG_HAB
#########################
path <- file.path(res_dir, "ST_IID_EFAG_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_agency_ef
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # ## adjust starting value
  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.1
  # newpar[["L_epsilon1_z"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3)#,
                # model_args = list(Map = Map),
                # optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    




#########################
## SST_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 0.3

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## SST_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 0.3

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## SST_RW_ALL_NOHAB
#########################
path <- file.path(res_dir, "SST_RW_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## SST_RW_GEAR_NOHAB
#########################
path <- file.path(res_dir, "SST_RW_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    



#########################
## ST_IID_ALL_HAB
#########################
path <- file.path(res_dir, "ST_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # Map <- fit0$tmb_list$Map
  # Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## ST_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "ST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 0.3

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## ST_RW_ALL_HAB
#########################
path <- file.path(res_dir, "ST_RW_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # Map <- fit0$tmb_list$Map
  # Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## ST_RW_GEAR_HAB
#########################
path <- file.path(res_dir, "ST_RW_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "EF"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- NULL 
  }
  if(msubx[3] == "EFAG"){ 
    Data_Geostat_inp <- Data_Geostat_ef
    Q_ik_inp <- Q_ik_efagency 
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 0.3

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    



df <- data.frame("Model"=c("SST_IID_ALL_HAB", "SST_IID_GEAR_HAB", "SST_IID_ALL_NOHAB", "SST_IID_GEAR_NOHAB", "SST_RW_ALL_HAB", "SST_RW_GEAR_HAB", "ST_IID_ALL_HAB", "ST_IID_GEAR_HAB","ST_RW_ALL_HAB","ST_RW_GEAR_HAB", "SST_IID_EF_HAB"))
df$AIC <- NA

for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df <- df %>% mutate(dAIC = AIC - min(AIC, na.rm=TRUE))
df[order(df$dAIC),]















#########################
## BASE, but RW
#########################
path <- file.path(res_dir, "SST_RW_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  # ## adjust starting value
  # newpar <- fit1$parameter_estimates$par
  # newpar[["logkappa1"]] <- 0.3

  # fit2 = fit_model( "settings"=settings, 
  #                 "Lat_i"=Data_Geostat_inp[,"Lat"], 
  #                 "Lon_i"=Data_Geostat_inp[,"Lon"], 
  #                 "t_iz"=Data_Geostat_inp[,'Year'], 
  #                 "c_i"=rep(0,nrow(Data_Geostat_inp)), 
  #                 "b_i"=Data_Geostat_inp[,'Catch_KG'], 
  #                 "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
  #                 "v_i"=Data_Geostat_inp[,'Vessel'], 
  #                 working_dir=path, 
  #                 extrapolation_args=list(
  #                   input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
  #                   Network_sz_LL=Network_sz_LL),
  #                 Network_sz = Network_sz,
  #                 Xconfig_zcp = Xconfig_zcp_inp,
  #                 X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
  #                 Q_ik = Q_ik_inp, 
  #                 model_args = list(Map = Map),
  #                 optimize_args = list(getsd=FALSE, startpar = newpar))
  # check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## Base but NOHAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


#########################
## Base but RW NOHAB
#########################
path <- file.path(res_dir, "SST_RW_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


#########################
## SPATIAL ONLY
#########################
path <- file.path(res_dir, "ST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    



#########################
## ELECTROFISHING ONLY
#########################
path <- file.path(res_dir, "SST_IID_EF_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # Map <- fit0$tmb_list$Map
  # Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3)
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15))    




#########################
## GEAR * AGENCY
#########################
path <- file.path(res_dir, "SST_IID_GEARAG_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(grepl("EF", msubx[3])){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }
  if(msubx[3] == "GEAR"){ Q_ik_inp <- Q_ik_method }
  if(msubx[3] == "ALL"){ Q_ik_inp <- NULL}
  if(msubx[3] == "GEARAG"){ Q_ik_inp <- Q_ik_methodagency }
  if(msubx[3] == "EF"){ Q_ik_inp <- NULL }
  if(msubx[3] == "EFAG"){ Q_ik_inp <- Q_ik_efagency }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))
  Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15))    






















#########################
## model run directories
#########################

config_vec <- c("Temporal", "Spatial+Temporal", "Spatiotemporal+Spatial+Temporal")
smooth_vec <- c("IID", "RandomWalk")
vessel_vec <- c("AllGears", "Gear", "Electrofishing", "ElectrofishingAgency")#, "VesselAgency", )
hab_vec <- c("Habitat","NoHabitat")

models <- expand.grid("SpatialConfig"=config_vec, "TemporalSmoother"=smooth_vec, "Vessel"=vessel_vec, "Habitat"=hab_vec)

model_sub <- models %>% filter(SpatialConfig=="Spatial+Temporal")

ignore <- sapply(1:length(config_vec), function(aa){
  config_mod <- file.path(res_dir, config_vec[aa])
  dir.create(config_mod, showWarnings=FALSE)

  ignore2 <- sapply(1:length(smooth_vec), function(bb){
    smooth_mod <- file.path(config_mod, smooth_vec[bb])
    dir.create(smooth_mod, showWarnings=FALSE)

    ignore3 <- sapply(1:length(vessel_vec), function(cc){
      vessel_mod <- file.path(smooth_mod, vessel_vec[cc])
      dir.create(vessel_mod, showWarnings=FALSE)

      ignore4 <- sapply(1:length(hab_vec), function(dd){
        hab_mod <- file.path(vessel_mod, hab_vec[dd])
        hab_fig <- file.path(hab_mod, "figures")
        dir.create(hab_mod, showWarnings=FALSE)
        dir.create(hab_fig, showWarnings=FALSE)

        if(file.exists(file.path(res_dir, "VAST_v8_0_0.cpp"))){
          file.copy(from=file.path(res_dir, "VAST_v8_0_0.cpp"), to=hab_mod)
        }
        if(file.exists(file.path(res_dir, "VAST_v8_0_0.dll"))){
          file.copy(from=file.path(res_dir, "VAST_v8_0_0.dll"), to=hab_mod)
        }
        if(file.exists(file.path(res_dir, "VAST_v8_0_0.o"))){
          file.copy(from=file.path(res_dir, "VAST_v8_0_0.o"), to=hab_mod)
        }
      })
    })
  })
})



# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat, FileName = fig_dir)


###################################
## plot network and observations
###################################


# ## plot network and observations
# p1 <- ggplot() +
#   geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.5, alpha=0.5) +
#   geom_point(data = obs, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
#   xlab("Longitude") + ylab("Latitude") +
#   mytheme()
# # ggsave(file.path(fig_dir, "Waitaki.png"), p1)

# ## plot encounters and non-encounters on network by year
# p2 <- ggplot() +
#   geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.001, alpha=0.2) +
#   geom_point(data = obs, aes(x = long, y = lat, fill = factor(present)), pch=21, cex=1.6,alpha=0.5) +
#   scale_fill_brewer(palette = "Set1") +
#   scale_x_continuous(breaks = round(seq(min(network$long),max(network$long), length.out=3),0)) +
#   xlab("Longitude") + ylab("Latitude") + 
#   guides(fill=guide_legend(title="Encounter")) +
#   facet_wrap(~year) +
#   mytheme()
# # ggsave(file.path(fig_dir, "Waitaki_byYear.png"), p2)

#######################################
### VESSEL: All Gears
#######################################

#################################
## TEMPORAL SMOOTHER: Random walk
#################################

########################
## HABITAT: None
########################
### ALLGEARS/ST/RW/NOHAB
### spatial + temporal
## rw temporal smoother
## no habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

### ALLGEARS/SST/RW/NOHAB
### spatiotemporal + spatial + temporal
## rw temporal smoother
## no habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

########################
## HABITAT: Included
########################
### ALLGEARS/ST/RW/HAB
### spatial + temporal
## rw temporal smoother
## habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)
  Map <- fit0$tmb_list$Map
  Map[["gamma2_ctp"]] <- factor(rep(NA, length(Map[["gamma2_ctp"]])))
  Map[["beta2_ft"]] <- factor(rep(NA, length(Map[["beta2_ft"]])))

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[['logkappa1']] <- 0.1

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

### ALLGEARS/SST/RW/HAB
### spatiotemporal + spatial + temporal
## rw temporal smoother
## habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  new <- fit1$parameter_estimates$par
  new[["logkappa1"]] <- 0.01

  # ## wrapper function
  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, startpar = new))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=new))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

##################################
## TEMPORAL SMOOTHER: IID
##################################

###################
## HABITAT: None
###################

### ALLGEARS/ST/IID/NOHAB
### spatial + temporal
## iid temporal smoother
## no habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

### ALLGEARS/SST/ID/NOHAB
### spatiotemporal + spatial + temporal
## iid temporal smoother
## no habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

###################
## HABITAT: Included
###################
### ALLGEARS/ST/IID/HAB
### spatial + temporal
## iid temporal smoother
## habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

### ALLGEARS/SST/IID/HAB
### spatiotemporal + spatial + temporal
## iid temporal smoother
## habitat info
## all gears pooled
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "AllGears") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    



#######################################
### VESSEL: Gear effect
#######################################

#################################
### TEMPORAL SMOOTHER: Random walk
#################################

######################
## HABITAT: None
######################

### GEAR/ST/RW/NOHAB
### spatial + temporal
## rw temporal smoother
## no habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }


settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model=FALSE)
  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

### GEAR/SST/RW/NOHAB
### spatiotemporal + spatial + temporal
## rw temporal smoother
## no habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

##########################
### HABITAT: Included
##########################

### spatial + temporal
## rw temporal smoother
## habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }


  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

  new <- fit1$parameter_estimates$par
  new[["logkappa1"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, startpar=new))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj)

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=new))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

### spatiotemporal + spatial + temporal
## rw temporal smoother
## habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  new <- fit1$parameter_estimates$par
  new[["logkappa1"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, startpar=new))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj)

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=new))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

####################################
## TEMPORAL SMOOTHER: IID
####################################

#######################
### HABITAT: None
#######################

### spatial + temporal
## iid temporal smoother
## no habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


### spatiotemporal + spatial + temporal
## iid temporal smoother
## no habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

##########################
## HABITAT: Included
##########################

### spatial + temporal
## iid temporal smoother
## habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)


  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,15,16), covar_names=covar2)    


### spatiotemporal + spatial + temporal
## iid temporal smoother
## habitat info
## gear effect
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Gear") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }


  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  new <- fit1$parameter_estimates$par
  new[["logkappa1"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, startpar=new))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj)

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit2$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,15,16), covar_names=covar2)    


#######################################
### VESSEL: Electrofishing only
#######################################

#################################
## TEMPORAL SMOOTHER: Random walk
#################################

########################
## HABITAT: None
########################

### spatial + temporal
## rw temporal smoother
## no habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


### spatiotemporal + spatial + temporal
## rw temporal smoother
## no habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

########################
## HABITAT: Included
########################

### spatial + temporal
## rw temporal smoother
## habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    


### spatiotemporal + spatial + temporal
## rw temporal smoother
## habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "RandomWalk") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    

##################################
## TEMPORAL SMOOTHER: IID
##################################

###################
## HABITAT: None
###################

### spatial + temporal
## iid temporal smoother
## no habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


### spatiotemporal + spatial + temporal
## iid temporal smoother
## no habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "NoHabitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

###################
## HABITAT: Included
###################

### spatial + temporal
## iid temporal smoother
## habitat info
## electrofishing only
msub <- models %>% 
        filter(SpatialConfig == "Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    


### spatiotemporal + spatial + temporal
## iid temporal smoother
## habitat info
## electrofishing gear only
msub <- models %>% 
        filter(SpatialConfig == "Spatiotemporal+Spatial+Temporal") %>%
        filter(TemporalSmoother == "IID") %>%
        filter(Vessel == "Electrofishing") %>%
        filter(Habitat == "Habitat")
path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
setwd(path)
fig <- file.path(path, "figures")

  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }

  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      X_itp_inp <- X_itp_ef
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }

  settings <- model_setup(model_info = msub, Method = "Stream_network")

  # ## wrapper function
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat[,"Lat"], 
                  "Lon_i"=Data_Geostat[,"Lon"], 
                  "t_iz"=Data_Geostat[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat)), 
                  "b_i"=Data_Geostat[,'Catch_KG'], 
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat[,"Lat"], 
                "Lon_i"=Data_Geostat[,"Lon"], 
                "t_iz"=Data_Geostat[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat)), 
                "b_i"=Data_Geostat[,'Catch_KG'], 
                "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names=covar2)    





































noc <- NULL
for(loop in 1:nrow(models)){

  msub <- models[loop,]
  path <- file.path(res_dir, msub$SpatialConfig, msub$TemporalSmoother, msub$Vessel, msub$Habitat)
  setwd(path)
  fig <- file.path(path, "figures")

  if(file.exists(file.path(path, "Fit.rds"))){
      fit <- readRDS(file.path(path, "Fit.rds"))

      if(msub$Habitat == "NoHabitat"){
        plot_set <- c(1,6,8)
      } else { plot_set <- c(1,6,8,11,13,15)}
      
      plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=plot_set)
      next
  }


  if(grepl("Electrofishing", msub$Vessel)){
    Data_Geostat_inp <- Data_Geostat %>% filter(Vessel == "Electric fishing")
  } else {
    Data_Geostat_inp <- Data_Geostat
  }

  if(msub$Vessel == "Gear"){
    Q_ik_inp <- Q_ik_method
  }
  if(msub$Vessel == "AllGears"){
    Q_ik_inp <- NULL
  }
  if(msub$Vessel == "ElectrofishingAgency"){
    Q_ik_inp <- Q_ik_efagency
  }
  if(msub$Vessel == "Electrofishing"){
    Q_ik_inp <- NULL
  }


  if(msub$Habitat == "Habitat"){
    X_gtp_inp <- X_gtp_input

    if(grepl("Electrofishing", msub$Vessel)){
      n_i <- nrow(Data_Geostat_inp)
      X_itp_inp <- array(0, dim=c(n_i,n_t,n_p))
      for(i in 1:n_i){
        for(p in 1:n_p){
          child_i <- Data_Geostat_inp$Knot[i]
          index <- which(nodes == child_i)
          X_itp_inp[i,,p] <- X_gtp_input[index,,p]
        }
      }
    } else {
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
  }


  settings <- model_setup(model_info = msub, Method = "Stream_network")

  ## first try
  fit0 = tryCatch(fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE)), error = function(e) NA)

  ## if first try has issues, skip for now
  if(all(is.na(fit1))==FALSE){
    check1 <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
  } else { 
    noc <- c(noc, paste0(loop,"_fit1"))
    next 
  }

  ## if not all parameters are identifiable, fit again with a new starting value for logkappa1
  if(length(check1$WhichBad)>0){
    new <- fit1$parameter_estimates$par
    new[["logkappa1"]] <- 0.01  

    fit2 = tryCatch(fit_model( "settings"=settings, 
                    "Lat_i"=Data_Geostat_inp[,"Lat"], 
                    "Lon_i"=Data_Geostat_inp[,"Lon"], 
                    "t_iz"=Data_Geostat_inp[,'Year'], 
                    "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                    "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                    "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                    "v_i"=Data_Geostat_inp[,'Vessel'], 
                    working_dir=path, 
                    extrapolation_args=list(
                      input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                      Network_sz_LL=Network_sz_LL),
                    Network_sz = Network_sz,
                    X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                    Q_ik = Q_ik_inp, 
                    optimize_args = list(getsd=FALSE, startpar=new)), error = function(e) NA)

    ## if second try has issues still, skip for now
    if(all(is.na(fit2))==FALSE){
      check2 <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj)
    } else {
      noc <- c(noc, paste0(loop, "_fit2"))
      next
    }

    ## if second try has non-identifiable parameters, skip for now
    if(length(check2$WhichBad)>0){
      noc <- c(noc, paste0(loop,"_check2"))
      next
    }
  } else {
    ## if they are identifiable, starting values for model run are the first run's estimates
    new <- fit1$parameter_estimates$par
  }

  fit = tryCatch(fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                optimize_args = list(startpar=new)), error = function(e) NA)
  if(all(is.na(fit))==FALSE){
      saveRDS(fit, file.path(path, "Fit.rds"))    
      fit <- readRDS(file.path(path, "Fit.rds"))

      if(msub$Habitat == "NoHabitat"){
        plot_set <- c(1,6,8)
      } else { plot_set <- c(1,6,8,11,13,15)}
      
      plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=plot_set)    
  } else {
    noc <- c(noc, paste0(loop, "_fit"))
  }
}















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
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
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
                  input_grid=cbind("Lat"=Data_Geostat[,"Lat"], "Lon"=Data_Geostat[,"Lon"],"child_i"=Data_Geostat[,"Knot"],"Area_km2"=Data_Geostat[,"AreaSwept_km2"]), 
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                newtonsteps=3)
saveRDS(fit, file.path(st_dir, "Fit.rds"))

fit <- readRDS(file.path(st_dir, "Fit.rds"))

# Plot results
plot_results( settings=settings_st, fit=fit, working_dir=stfig_dir, category_names="Longfin eels", strata_names = "Longfin eels", plot_set=1 )


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
plot_results( settings=settings_sst, fit=fit, working_dir=sstfig_dir,category_names="Longfin eels", strata_names = "Longfin eels", plot_set=1 )


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
