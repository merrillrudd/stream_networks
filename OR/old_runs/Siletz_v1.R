rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks"
R_dir <- file.path(main_dir, "R")
R_files <- list.files(R_dir)
readr <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

or_dir <- file.path(main_dir, "OR")

res_dir <- file.path(or_dir, "Siletz")
dir.create(res_dir, showWarnings=FALSE)

data_dir <- file.path(or_dir, "data")

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
##remember to pull upstream development branches

devtools::install_github("merrillrudd/FishStatsUtils", ref='stream')
library(FishStatsUtils)


# devtools::install_github("james-thorson/VAST", ref="master")
library(VAST)

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

network_raw <- read.csv(file.path(data_dir, "Siletz", "Network_WGS.csv"))
nodes <- read.csv(file.path(data_dir, "Siletz", "Nodes_WGS.csv"))
hab_raw <- read.csv(file.path(data_dir, "Siletz", "Habitat_Siletz_WGS_line_final.csv"))
juv <- read.csv(file.path(data_dir, "Siletz", "Juvenile_Siletz_WGS_line_final.csv"))
spawn <- read.csv(file.path(data_dir, "Siletz", "Spawning_Siletz_WGS_line_final.csv"))

##############
##### NETWORK
##############
net <- network_raw %>% 
  select("Shape_Leng", "ChildNode", "ParentNode", "lat", "long") %>%
  rename("dist_s"="Shape_Leng", "child_s"="ChildNode", "parent_s"="ParentNode") %>%
  mutate("dist_s" = dist_s / 1000)
nod <- nodes %>% select("NodeID", "lat", "long")

if(length(net$parent_s %in% net$child_s) > 0){
  iroot <- which(net$parent_s %in% net$child_s == FALSE)
  new <- lapply(1:length(iroot), function(x){
    sub <- net[iroot[x],]
    nodsub <- nod %>% filter(NodeID == sub$parent_s)
    df <- data.frame("dist_s"=Inf, "child_s"=sub$parent_s, "parent_s"=0, "lat"=nodsub$lat, "long"=nodsub$long)
    return(df)
  })
  add <- do.call(rbind, new)
  network <- rbind.data.frame(net, add)
} else {
  network <- net
}

### check for unique locations
nrow(network)
length(unique(network$lat))
length(unique(network$long))

## change one of the longitudes to be a unique location by adding a tiny amount
network %>% filter(long == names(table(network$long))[which(table(network$long)==max(table(network$long)))])
index <- which(network$long == names(table(network$long))[which(table(network$long)==max(table(network$long)))])
network$long[index[2]] <- network$long[index[2]] + 1e-4
network$long[index]

## recheck
nrow(network)
length(unique(network$lat))
length(unique(network$long))

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)



##############
##### SURVEYS
##############

obs_spawn <- spawn %>% 
    select("SpawnYear", "adults_per_km", "lat", "long", "ChildNode") %>%
    rename("year"="SpawnYear", "observation"="adults_per_km","child_i"="ChildNode") %>%
    mutate("obs_type"="spawners") %>%
    mutate('surveynum' = 1) %>%
    mutate('survey'="spawners")

obs_juv <- juv %>%
    select("JuvYear", "parr_per_k", "lat", "long","ChildNode") %>%
    rename("year"="JuvYear", "observation"="parr_per_k","child_i"="ChildNode") %>%
    na.omit() %>%
    mutate("obs_type"="juveniles") %>% 
    mutate("surveynum" = 2) %>%
    mutate("survey" = "juveniles")

obs_dens <- rbind.data.frame(obs_spawn, obs_juv)
category_names <- unique(obs_dens$survey)
years <- unique(obs_dens$year)
n_t <- length(years)
n_x <- nrow(Network_sz)

hab <- hab_raw %>% 
  select("lat","long","YEAR_", "ChildNode","PCTSCCHNLA","ACW","PCTSWPOOL","POOL1P_KM","POOLS100","RIFFLEDEP","LWDVOL1","WGTED_SLOPE_GRAVEL","WGTED_ALLUNITS_BEDROCK") %>%
  rename("child_i"="ChildNode", "year"="YEAR_")
dcovar <- c("PCTSCCHNLA","ACW","PCTSWPOOL","POOL1P_KM","POOLS100","RIFFLEDEP","LWDVOL1","WGTED_SLOPE_GRAVEL","WGTED_ALLUNITS_BEDROCK")

hab <- hab %>% tidyr::gather(key="obs_type", value="observation", PCTSCCHNLA:WGTED_ALLUNITS_BEDROCK) %>% mutate('survey'="habitat")
hab_addYr <- lapply(1:length(years), function(x){
  out <- hab
  out$year <- years[x]
  return(out)
})
hab_addYr <- do.call(rbind, hab_addYr)
hab <- hab_addYr %>% mutate('surveynum'=3)

  obs_all <- rbind.data.frame(obs_dens, hab)



###################################
## plot network and observations
###################################

## plot network and observations
p1 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=1, alpha=0.5) +
  geom_point(data = obs_dens, aes(x = long, y = lat), pch=19, cex=2, color="red", alpha=0.8) +
  geom_point(data = network %>% filter(parent_s == 0), aes(x = long, y = lat), pch=19, cex=3, color="goldenrod") +
  xlab("Longitude") + ylab("Latitude") +
  mytheme()
# ggsave(file.path(fig_dir, "Siletz.png"), p1)

## plot encounters and non-encounters on network by year
p2 <- ggplot() +
  geom_point(data = network, aes(x = long, y = lat), pch=19, cex=0.001, alpha=0.2) +
  geom_point(data = obs_dens %>% filter(survey == "spawners"), aes(x = long, y = lat, fill = density), pch=21, cex=2, alpha=0.8) +
  geom_point(data = obs_dens %>% filter(survey == "juveniles"), aes(x = long, y = lat, color = density), pch=22, cex=2) +
  # scale_fill_brewer(palette = "Set1") +
  scale_fill_viridis_c() +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = round(seq(min(network$long),max(network$long), length.out=3),0)) +
  xlab("Longitude") + ylab("Latitude") + 
  guides(fill=guide_legend(title="spawners"), color=guide_legend(title="juveniles")) +
  facet_wrap(~year) +
  mytheme()
# ggsave(file.path(fig_dir, "Siletz_byYear.png"), p2)


##################################
## save data used for model runs
##################################

saveRDS(obs_dens, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))


##### setup data frame
## density = individuals per kilometer
Data_Geostat <- data.frame( "Catch_KG" = obs_dens$observation, 
              "Year" = as.numeric(obs_dens$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_dens$lat, 
               "Lon" = obs_dens$long, 
               "Pass" = 0,
               "Knot" = obs_dens$child_i,
               "Category" = obs_dens$survey,
               "CategoryNum"=obs_dens$surveynum)

Data_Geostat_all <- data.frame( "Catch_KG" = obs_all$observation, 
              "Year" = as.numeric(obs_all$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_all$lat, 
               "Lon" = obs_all$long, 
               "Pass" = 0,
               "Knot" = obs_all$child_i,
               "Category" = obs_all$survey)

Data_Geostat_hab <- data.frame( "Catch_KG" = hab$observation, 
              "Year" = as.numeric(hab$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = hab$lat, 
               "Lon" = hab$long, 
               "Pass" = 0,
               "Knot" = hab$child_i,
               "Category" = hab$survey)

plot_network(Network_sz_LL, Data_Geostat_all, FileName = fig_dir, root=TRUE)
plot_network(Network_sz_LL, Data_Geostat, FileName = fig_dir, root=TRUE)
plot_network(Network_sz_LL, Data_Geostat_hab, FileName = fig_dir, root=TRUE)


years <- unique(Data_Geostat$Year)[order(unique(Data_Geostat$Year))]
## check data
find <- sapply(1:length(years), function(x){
  sub <- Data_Geostat %>% filter(Year==years[x])
  sub1 <- sub %>% filter(Category==1)
  sub2 <- sub %>% filter(Category==2)
  df <- data.frame("Spawner"=length(which(sub1$Catch_KG == 0))/nrow(sub1), "Juvenile"=length(which(sub2$Catch_KG == 0))/nrow(sub2))
  return(df)
})

##################################
## Models
##################################

#########################
## SST_MULTIVARIATE_LOGNORMAL_TEMPIID
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_LOGNORMAL_TEMPIID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=1,"Omega2"=0)
RhoConfig_inp <- c("Beta1"=1,"Beta2"=1,"Omega1"=1,"Omega2"=1)
ObsModel_inp <- c(1,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
  ### GRADIENTS ARE HIGH

#########################
## SST_MULTIVARIATE_GAMMA_TEMPIID
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPIID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=1,"Omega2"=0)
RhoConfig_inp <- c("Beta1"=1,"Beta2"=1,"Omega1"=1,"Omega2"=1)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
  ## HIGH GRADIENTS

#########################
## SST_MULTIVARIATE_LOGNORMAL_TEMPRW
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_LOGNORMAL_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
# RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=2,"Omega2"=0)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=2,"Omega1"=2,"Omega2"=2)
ObsModel_inp <- c(1,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
  ### GRADIENTS ARE HIGH

  Map <- fit0$tmb_list$Map
  Map[["logSigmaM"]] <- factor(rep(NA, length(Map[["logSigmaM"]])))

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  newpar <- fit2$parameter_estimates$par
  newpar[["logkappa2"]] <- -3

  fit3 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar=newpar))
  check <- TMBhelper::Check_Identifiable(fit3$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=year_labels, Years2Include=years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))


#########################
## SST_MULTIVARIATE_GAMMA_TEMPRW
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=1,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=1,"Omega1"=1,"Omega2"=1)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=year_labels, Years2Include=years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

#########################
## SST_SPAWNERS_GAMMA_TEMPRW
#########################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=1,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=1,"Omega1"=1,"Omega2"=1)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1])

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)

#########################
## SST_JUVENILES_GAMMA_TEMPRW
#########################
path <- file.path(res_dir, "SST_JUVENILES_GAMMA_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum == 2)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=1,"Omega2"=0)
# RhoConfig_inp <- c("Beta1"=1,"Beta2"=1,"Omega1"=1,"Omega2"=1)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-2, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-2, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-2, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[2] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[2])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[2])

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#########################
## SST_MULTIVARIATE_GAMMA_TEMPAR1
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPAR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=4,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("SD_site_logdensity"=1, "Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.7

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar=newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj)   

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    
 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,10,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)


  Dens_xt = plot_maps(plot_set=c(10, 3,12), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Panel="Year")

  # MappingDetails=map_list[["MappingDetails"]]
  # TmbData=fit$data_list
  # spatial_list=fit$spatial_list
  # Report=fit$Report
  # Sdreport=fit$parameter_estimates$SD
  # PlotDF=map_list[["PlotDF"]]
  # MapSizeRatio=map_list[["MapSizeRatio"]]
  # Xlim=map_list[["Xlim"]]
  # Ylim=map_list[["Ylim"]]
  # FileName=fig
  # Year_Set=fit$year_labels
  # Years2Include=fit$years_to_plot
  # Rotate=map_list[["Rotate"]]
  # Cex=map_list[["Cex"]]
  # Legend=map_list[["Legend"]]
  # zone=map_list[["Zone"]]
  # mar=c(0,0,2,0)
  # oma=c(3.5,3.5,0,0)
  # cex=1.8
  # plot_legend_fig=FALSE
  # category_names=category_names
  # Panel="Year"

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))



#########################
## SST_MULTIVARIATE_GAMMA_TEMPRW_FACTOR
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPRW_FACTOR")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"=2,"Omega2"=2,"Epsilon1"=2,"Epsilon2"=2)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=2,"Omega1"=2,"Omega2"=2)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    
 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))


#########################
## SST_MULTIVARIATE_GAMMA_TEMPAR1_FACTOR
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPAR1_FACTOR")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"=2,"Omega2"=2,"Epsilon1"=2,"Epsilon2"=2)
RhoConfig_inp <- c("Beta1"=4,"Beta2"=4,"Omega1"=4,"Omega2"=4)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.7

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar=newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    
 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

  plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))


#########################
## ST_MULTIVARIATE_GAMMA_TEMPAR1
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_TEMPAR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=4,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    
 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

#########################
## ST_MULTIVARIATE_GAMMA_TEMPRW
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    
 
    map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,12), MappingDetails=map_list[["MappingDetails"]], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, plot_legend_fig=FALSE,category_names=category_names)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, MappingDetails=map_list[["MappingDetails"]], PlotDF=map_list[["PlotDF"]], MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Rotate=map_list[["Rotate"]], Cex=map_list[["Cex"]], Legend=map_list[["Legend"]], zone=map_list[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))



df <- data.frame('Model'=c('SST_MULTIVARIATE_LOGNORMAL_TEMPRW',
                           'SST_MULTIVARIATE_GAMMA_TEMPRW',
                           "SST_MULTIVARIATE_GAMMA_TEMPAR1",
                           "SST_MULTIVARIATE_GAMMA_TEMPRW_FACTOR",
                           "SST_MULTIVARIATE_GAMMA_TEMPAR1_FACTOR",
                           "ST_MULTIVARIATE_GAMMA_TEMPAR1",
                           "ST_MULTIVARIATE_GAMMA_TEMPRW"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- sapply(1:nrow(df), function(x) df$AIC[x] - min(df$AIC))
df[order(df$deltaAIC),]













#########################
## SST_SPAWNERS_LOGNORMAL_Beta1IID
#########################
path <- file.path(res_dir, "SST_SPAWNERS_LOGNORMAL_Beta1IID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(1,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names[1], plot_set=c(1,2,3,4,5,6,7,8,9,12))  

#########################
## SST_SPAWNERS_GAMMA_Beta1IID
#########################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_Beta1IID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=1,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names[1], plot_set=c(1,2,3,4,5,6,7,8,9,12))  


#########################
## SST_SPAWNERS_GAMMA_Beta1IID
#########################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_Beta1RW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"="IID","Epsilon2"="IID")
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names[1], plot_set=c(1,2,3,4,5,6,7,8,9,12))  




df <- data.frame("Model" = c("SST_SPAWNERS_LOGNORMAL_Beta1IID",
                            "SST_SPAWNERS_GAMMA_Beta1IID",
                            "SST_SPAWNERS_GAMMA_Beta1RW"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df$Model[i], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- sapply(1:nrow(df), function(x) df$AIC[x]-min(df$AIC))
df[order(df$deltaAIC),]





#########################
## ST_MULTIVARIATE_GAMMA_Beta1RW
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_Beta1RW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names, plot_set=c(1,2,3,4,5,6,7,8,9,12))  


#########################
## ST_MULTIVARIATE_GAMMA_Beta1AR1
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_Beta1AR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=4,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names, plot_set=c(1,2,3,4,5,6,7,8,9,12))  



#########################
## ST_MULTIVARIATE_GAMMA_Beta1AR1_SurveyType
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_Beta1AR1_SurveyType")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=4,"Beta2"=4,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Q_ik = Q_ik_inp,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names, plot_set=c(1,2,3,4,5,6,7,8,9,12))  



#########################
## ST_MULTIVARIATE_GAMMA_Beta1RW_SurveyType
#########################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_Beta1RW_SurveyType")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"="IID","Omega2"="IID","Epsilon1"=0,"Epsilon2"=0)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=2,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_inp, 
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Q_ik = Q_ik_inp,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names, plot_set=c(1,2,3,4,5,6,7,8,9,12))  



#########################
## SST_MULTIVARIATE_GAMMA_Beta1RW_Factor
#########################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_Beta1RW_Factor")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(Category == 1)

FieldConfig_inp <- c("Omega1"=2,"Omega2"=2,"Epsilon1"=2,"Epsilon2"=2)
RhoConfig_inp <- c("Beta1"=2,"Beta2"=0,"Omega1"=0,"Omega2"=0)
ObsModel_inp <- c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=c(0,0), Options=c("Calculate_Range"=1,"Calculate_effective_area"=1), ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=Data_Geostat_inp[,"CategoryNum"]-1, 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names=category_names, plot_set=c(1,2,3,4,5,6,7,8,9,12))  


df <- data.frame("Model" = c("SST_MULTIVARIATE_GAMMA_BETA1RW_FACTOR",
                            "SST_MULTIVARIATE_GAMMA_Beta1RW",
                            "SST_MULTIVARIATE_GAMMA_Beta1AR1",
                            "ST_MULTIVARIATE_GAMMA_Beta1RW",
                            "ST_MULTIVARIATE_GAMMA_Beta1AR1"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df$Model[i], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- sapply(1:nrow(df), function(x) df$AIC[x]-min(df$AIC))
df[order(df$deltaAIC),]













































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
