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

### create previously saved data into Rdata file
obs <- readRDS(file.path(res_dir, "observations.rds"))
net <- readRDS(file.path(res_dir, "network.rds"))
hab <- readRDS(file.path(res_dir, "habitat.rds"))

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
    rename("year"="SpawnYear", "value"="adults_per_km","child_i"="ChildNode") %>%
    mutate('surveynum' = 1) %>%
    mutate('survey'="spawners")

obs_juv <- juv %>%
    select("JuvYear", "parr_per_k", "lat", "long","ChildNode") %>%
    rename("year"="JuvYear", "value"="parr_per_k","child_i"="ChildNode") %>%
    na.omit() %>%
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
dcovar_names <- c("Percent secondary channel area", "Active channel width", "Percent slack water pools","Deep pools per km","Pools per 100m reach","Riffle depth", "Volume large wood 100m","Weighted percent gravel in riffles","Weighted percent bedrock")

hab_df <- hab %>% tidyr::gather(key="variable", value="value", PCTSCCHNLA:WGTED_ALLUNITS_BEDROCK) %>% mutate('survey'="habitat")
saveRDS(hab_df, file.path(res_dir, "habitat.rds"))

hab_addYr <- lapply(1:length(years), function(x){
  out <- hab_df
  out$year <- years[x]
  return(out)
})
hab_addYr <- do.call(rbind, hab_addYr)
hab_addYr <- hab_addYr %>% mutate('surveynum'=3)


##################################
## save data used for model runs
##################################

saveRDS(obs_dens, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))

##### setup data frame
## density = individuals per kilometer
Data_Geostat <- data.frame( "Catch_KG" = obs_dens$value, 
              "Year" = as.numeric(obs_dens$year),
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
               "Lat" = obs_dens$lat, 
               "Lon" = obs_dens$long, 
               "Pass" = 0,
               "Knot" = obs_dens$child_i,
               "Category" = obs_dens$survey,
               "CategoryNum"=obs_dens$surveynum)

plot_network(Network_sz_LL, Data_Geostat, FileName=fig_dir, root=TRUE)

#######################################
## Habitat by network node
#######################################

hab_interp <- lapply(1:length(dcovar), function(p){
	## smooth across space
	habsub <- hab_df %>% filter(variable == dcovar[p])	

	interp_lat <- habsub$lat
	interp_lon <- habsub$long
	interp_z <- habsub$value	

	p1 <- ggplot(habsub) +
	geom_point(data=Network_sz_LL, aes(x=Lon,y=Lat), color="gray", alpha=0.5)+
	geom_point(aes(x = long, y = lat, color=value), cex=3) +
	ggtitle(paste0(dcovar_names[p], " (", dcovar[p], ")")) +
	xlab("Longitude") + ylab("Latitude") +
	scale_color_viridis_c() +
	mytheme()
	ggsave(file.path(fig_dir, paste0(dcovar[p], ".png")),p1)	

	find_lat <- Network_sz_LL$Lat[which(Network_sz_LL$child_s %in% habsub$child_i == FALSE)]
	find_lon <- Network_sz_LL$Lon[which(Network_sz_LL$child_s %in% habsub$child_i == FALSE)]
	find_child <- Network_sz_LL$child_s[which(Network_sz_LL$child_s %in% habsub$child_i == FALSE)]	

	# compute <- akima::interp(x = interp_lon, y = interp_lat, z = interp_z, xo = find_lon, yo=find_lat, linear=FALSE, extrap=TRUE, duplicate = "mean")
	compute <- akima::interpp(x = interp_lat, y = interp_lon, z = interp_z, xo=find_lat, yo=find_lon, duplicate = "mean", extrap=TRUE)	

	interp_df <- data.frame('lat'=find_lat, 'lon'=find_lon, 'value'=compute$z, 'child_i'=find_child)		

	p2 <- ggplot(interp_df)+
	geom_point(aes(x=lon,y=lat,color=value)) +
	geom_point(data=habsub, aes(x=long, y=lat, fill=value), pch=22, cex=3) +
	guides(fill = FALSE) +
	ggtitle(paste0(dcovar_names[p], " (", dcovar[p], ")", " interpolated (v1)")) +
	xlab("Longitude") + ylab("Longitude") +
	scale_color_viridis_c() +
	scale_fill_viridis_c() +
	mytheme()	
	ggsave(file.path(fig_dir, paste0(dcovar[p], "_interpolated_v1.png")),p2)

	habsub2 <- habsub %>% select('lat','long','child_i','value') %>% rename('lon'=long)
	hab_info <- rbind.data.frame(interp_df, habsub2)
	obs_new <- sapply(1:nrow(hab_info), function(x) ifelse(is.na(hab_info$value[x]),mean(hab_info$value,na.rm=TRUE),hab_info$value[x]))
	hab_info$value <- obs_new	

	p3 <- ggplot(hab_info)+
	geom_point(aes(x=lon,y=lat,color=value)) +
	geom_point(data=habsub, aes(x=long, y=lat, fill=value), pch=22, cex=3) +
	guides(fill = FALSE) +
	ggtitle(paste0(dcovar_names[p], " (", dcovar[p], ")", " interpolated")) +
	xlab("Longitude") + ylab("Longitude") +
	scale_color_viridis_c() +
	scale_fill_viridis_c() +
	mytheme()	
	ggsave(file.path(fig_dir, paste0(dcovar[p], "_interpolated.png")),p3)

	hab_new <- lapply(1:nrow(Network_sz_LL), function(x){
		child <- Network_sz_LL$child_s[x]
		find_hab <- hab_info %>% filter(child_i==child)
		if(nrow(find_hab)==1) return(find_hab)
		if(nrow(find_hab)>1){
			find_hab$value <- mean(find_hab$value)
			return(find_hab[1,])
		}
	})
	hab_new <- do.call(rbind, hab_new)

	return(hab_new)
})

# saveRDS(hab_interp, file.path(res_dir, 'habitat_interp.rds'))
hab_interp <- readRDS(file.path(res_dir, "habitat_interp.rds"))

nodes <- Network_sz_LL$child_s
n_x <- nrow(Network_sz_LL)
n_t <- length(years)
n_p <- length(dcovar)
n_i <- nrow(obs_dens)

X_gtp_input <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:n_p){
	psub <- hab_interp[[p]]
	mat <- matrix(0, nrow=n_x, ncol=1)
	mat[psub$child_i,1] <- psub$observation
	mat_sd <- (mat - mean(mat))/sd(mat)
	X_gtp_input[,,p] <- mat_sd
}

X_itp_input <- array(0, dim=c(n_i,n_t,n_p))
for(i in 1:n_i){
	for(p in 1:n_p){
		child_i <- obs_dens$child_i[i]
		index <- which(nodes == child_i)
		X_itp_input[i,,p] <- X_gtp_input[index,,p]
	}
}


##################################
## Models
##################################
#######################################
## SST_SPAWNERS_GAMMA_TEMPIID
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPIID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)

#######################################
## SST_SPAWNERS_GAMMA_TEMPRW
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_SPAWNERS_GAMMA_TEMPRW_STRW
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPRW_STRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)



#######################################
## SST_SPAWNERS_GAMMA_TEMPIID_STRW
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPIID_STRW")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)



#######################################
## SST_SPAWNERS_GAMMA_TEMPAR1
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPAR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_SPAWNERS_GAMMA_TEMPAR1_STAR1
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPAR1_STAR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_SPAWNERS_GAMMA_TEMPIID_STAR1
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPIID_STAR1")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)



#######################################
## SST_SPAWNERS_GAMMA_TEMPIID_STRW_HAB
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPIID_STRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)
X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)

#######################################
## SST_SPAWNERS_GAMMA_TEMPRW_STRW_HAB
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPRW_STRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)
X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)



#######################################
## SST_SPAWNERS_GAMMA_TEMPIID_HAB
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPIID_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)
X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_SPAWNERS_GAMMA_TEMPRW_HAB
#######################################
path <- file.path(res_dir, "SST_SPAWNERS_GAMMA_TEMPRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)
X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                  X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
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
                X_gtp = X_gtp_inputs, X_itp=X_itp_inputs,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


 df <- data.frame("Model"=c("SST_SPAWNERS_GAMMA_TEMPIID",
 							"SST_SPAWNERS_GAMMA_TEMPRW",
 							"SST_SPAWNERS_GAMMA_TEMPRW_STRW",
 							"SST_SPAWNERS_GAMMA_TEMPIID_STRW",
 							"SST_SPAWNERS_GAMMA_TEMPAR1",
 							"SST_SPAWNERS_GAMMA_TEMPIID_HAB",
 							"SST_SPAWNERS_GAMMA_TEMPIID_STRW_HAB",
 							"SST_SPAWNERS_GAMMA_TEMPRW_HAB",
 							"SST_SPAWNERS_GAMMA_TEMPRW_STRW_HAB",
 							"ST_SPAWNERS_GAMMA_TEMPIID",
 							"ST_SPAWNERS_GAMMA_TEMPIID_HAB"))
 df$AIC <- NA
 for(i in 1:nrow(df)){
 	res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
 	df$AIC[i] <- res$parameter_estimates$AIC
 }
 df$deltaAIC <- df$AIC - min(df$AIC)
 df[order(df$deltaAIC),]



 df <- data.frame("Model"=c("SST_MULTIVARIATE_GAMMA_TEMPIID_HAB",
 							"SST_MULTIVARIATE_GAMMA_TEMPIID_STRW_HAB",
 							"SST_MULTIVARIATE_GAMMA_TEMPRW_STRW",
 							"SST_MULTIVARIATE_GAMMA_TEMPRW_STRW_HAB",
 							"ST_MULTIVARIATE_GAMMA_TEMPRW",
 							"SST_MULTIVARIATE_GAMMA_TEMPIID",
 							"SST_MULTIVARIATE_GAMMA_TEMPRW_HAB"))
 df$AIC <- NA
 for(i in 1:nrow(df)){
 	res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
 	df$AIC[i] <- res$parameter_estimates$AIC
 }
 df$deltaAIC <- df$AIC - min(df$AIC)
 df[order(df$deltaAIC),]


#######################################
## ST_SPAWNERS_GAMMA_TEMPIID
#######################################
path <- file.path(res_dir, "ST_SPAWNERS_GAMMA_TEMPIID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)



#######################################
## ST_SPAWNERS_GAMMA_TEMPIID_HAB
#######################################
path <- file.path(res_dir, "ST_SPAWNERS_GAMMA_TEMPIID_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat %>% filter(CategoryNum==1)
X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp=X_gtp_inputs, X_itp=X_itp_inputs,
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
                  X_gtp=X_gtp_inputs, X_itp=X_itp_inputs,
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
                X_gtp=X_gtp_inputs, X_itp=X_itp_inputs,
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1] )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names[1])

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,13,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5,covar_names=dcovar)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names[1], Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_MULTIVARIATE_GAMMA_TEMPIID_HAB
#######################################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPIID_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  Map <- fit0$tmb_list$Map
  Map$logSigmaM <- factor(c(1,NA,NA,NA,NA,NA))

  Par <- fit0$tmb_list$Parameters
  Par$logSigmaM[2,1] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  Map$L_epsilon1_z <- factor(c(1,NA))
  Map$L_epsilon2_z <- factor(c(1,NA))
  Map$L_beta2_z <- factor(c(1,NA))

  Par$L_epsilon1_z[2] <- 0
  Par$L_epsilon2_z[2] <- 0
  Par$L_beta2_z[2] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
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
                X_gtp = X_gtp_input, X_itp=X_itp_input,
                optimize_args = list(startpar= fit3$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_MULTIVARIATE_GAMMA_TEMPIID_STRW_HAB
#######################################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPIID_STRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat #%>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  Map <- fit0$tmb_list$Map
  Map$logSigmaM <- factor(c(1,NA,NA,NA,NA,NA))

  Par <- fit0$tmb_list$Parameters
  Par$logSigmaM[2,1] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  Map$L_epsilon2_z <- factor(c(1,NA))
  Map$L_beta2_z <- factor(c(1,NA))

  Par$L_epsilon2_z[2] <- 0
  Par$L_beta2_z[2] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
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
                X_gtp = X_gtp_input, X_itp=X_itp_input,
                model_args = list(Map = Map, Par=Par),
                optimize_args = list(startpar= fit3$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)




#######################################
## SST_MULTIVARIATE_GAMMA_TEMPRW_HAB
#######################################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  Map <- fit0$tmb_list$Map
  Map$logSigmaM <- factor(c(1,NA,NA,NA,NA,NA))

  Par <- fit0$tmb_list$Parameters
  Par$logSigmaM[2,1] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  Map$L_epsilon2_z <- factor(c(1,NA))

  Par$L_epsilon2_z[2] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
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
                X_gtp = X_gtp_input, X_itp=X_itp_input,
                optimize_args = list(startpar= fit3$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)


#######################################
## SST_MULTIVARIATE_GAMMA_TEMPRW_STRW_HAB
#######################################
path <- file.path(res_dir, "SST_MULTIVARIATE_GAMMA_TEMPRW_STRW_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat #%>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"="IID", "Epsilon2"="IID")
RhoConfig = c("Beta1"=2, "Beta2"=2, "Epsilon1"=2, "Epsilon2"=2)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                X_gtp = X_gtp_input, X_itp=X_itp_input,
                model_args = list(Map = Map, Par=Par),
                optimize_args = list(startpar= fit1$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)

#######################################
## ST_MULTIVARIATE_GAMMA_TEMPIID
#######################################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_TEMPIID")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  # X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  # X_gtp = X_gtp_input, X_itp=X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  Map <- fit0$tmb_list$Map
  Map$logSigmaM <- factor(c(1,NA,NA,NA,NA,NA))
  Map$L_beta2_z <- factor(c(1,NA))

  Par <- fit0$tmb_list$Parameters
  Par$logSigmaM[2,1] <- 0
  Par$L_beta2_z[2] <- 0

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
                  # X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
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
                # X_gtp = X_gtp_input, X_itp=X_itp_input,
                optimize_args = list(startpar= fit2$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)





#######################################
## ST_MULTIVARIATE_GAMMA_TEMPIID_HAB
#######################################
path <- file.path(res_dir, "ST_MULTIVARIATE_GAMMA_TEMPIID_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Siletz/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

Data_Geostat_inp <- Data_Geostat# %>% filter(CategoryNum==1)
# X_itp_inputs <- X_itp_input[which(Data_Geostat$CategoryNum==1),-length(years),]
# X_gtp_inputs <- X_gtp_input[,-length(years),]

FieldConfig = c("Omega1"="IID", "Epsilon1"="IID", "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)
ObsModel = c(2,0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  Map <- fit0$tmb_list$Map
  Map$logSigmaM <- factor(c(1,NA,NA,NA,NA,NA))

  Par <- fit0$tmb_list$Parameters
  Par$logSigmaM[2,1] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  Map$L_epsilon1_z <- factor(c(1,NA))
  Map$L_epsilon2_z <- factor(c(1,NA))
  Map$L_beta2_z <- factor(c(1,NA))

  Par$L_epsilon1_z[2] <- 0
  Par$L_epsilon2_z[2] <- 0
  Par$L_beta2_z[2] <- 0

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
                  X_gtp = X_gtp_input, X_itp=X_itp_input,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
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
                X_gtp = X_gtp_input, X_itp=X_itp_input,
                optimize_args = list(startpar= fit3$parameter_estimates$par))

   saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

 
  map_list = make_map_info( "Region"=settings$Region, "spatial_list"=fit$spatial_list, "Extrapolation_List"=fit$extrapolation_list )

  Enc_prob = plot_encounter_diagnostic( Report=fit$Report, Data_Geostat=cbind("Catch_KG"=fit$data_frame[,'b_i']), DirName=fig)

  Index = plot_biomass_index( DirName=fig, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, use_biascorr=FALSE, category_names=category_names )

  plot_range_index(Report=fit$Report, TmbData=fit$data_list, Sdreport=fit$parameter_estimates$SD, Znames=colnames(fit$data_list$Z_xm), PlotDir=fig, Year_Set=fit$years_to_plot, use_biascorr=TRUE, category_names=category_names)

  Q = plot_quantile_diagnostic( TmbData=fit$data_list, Report=fit$Report, FileName_PP="Posterior_Predictive", FileName_Phist="Posterior_Predictive-Histogram", FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", save_dir=fig )

  Dens_xt = plot_maps(plot_set=c(1,2,3,4,5,6,7,8,9,11,15), TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5, covar_names=dcovar_names)

  Dens_xt = plot_maps(plot_set=12, TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Sdreport=fit$parameter_estimates$SD, MapSizeRatio=map_list[["MapSizeRatio"]], Xlim=map_list[["Xlim"]], Ylim=map_list[["Ylim"]], FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, category_names=category_names, Cex=0.5)

    plot_residuals(Lat_i=fit$data_frame[,'Lat_i'], Lon_i=fit$data_frame[,'Lon_i'], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, FileName=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot)






















