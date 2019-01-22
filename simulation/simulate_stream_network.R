rm(list=ls())

# devtools::install_github("james-thorson/VAST", ref = "development")
# devtools::install_github("merrillrudd/FishStatsUtils")
library(VAST)
library(FishStatsUtils)
library(TMB)
library(tidyverse)
library(RuddR)

wd <- "C:\\merrill\\stream_networks\\examples"
setwd(wd)

data_dir <- file.path(wd, "data")
dir.create(data_dir)

fig_dir <- file.path(wd, "figures")
dir.create(fig_dir)


##################################
## simulated stream network data
##################################
# Network parentage
set.seed(1)
parent_x = c( 0, 1, 2, 3, 3, 5, 5, 5, 8, 9, 10 )
dist_x = c( Inf, rlnorm(length(parent_x)-1, meanlog=0, sdlog=1) )
Network_xn = cbind( parent_x=parent_x, child_x=1:length(parent_x), dist_x=dist_x )
Obs_i = 1:nrow(Network_xn)
loc_xy = matrix(c( 0,0, 1,0, 2,0, 2,-1, 2,1, 1,1, 2,2, 3,1, 3,0, 3,-1, 3,-2), byrow=TRUE, ncol=2, dimnames=list(NULL,c("x","y")) )
loc_LL = matrix(loc_xy*0.1 + outer(rep(1,nrow(loc_xy)),c(20,40)), ncol=2, dimnames=list(NULL,c("Lon","Lat")))
xlab <- c(NA, sapply(2:nrow(loc_xy), function(x) mean(c(loc_xy[parent_x[x],'x'],loc_xy[x,'x']))))
ylab <- c(NA, sapply(2:nrow(loc_xy), function(x) mean(c(loc_xy[parent_x[x],'y'],loc_xy[x,'y']))))

# Distance matrix
D = array(NA, dim=rep(length(dist_x),2) )
D[ Network_xn[-1,c("parent_x","child_x")] ] = Network_xn[-1,'dist_x']
D[ Network_xn[-1,c("child_x","parent_x")] ] = Network_xn[-1,'dist_x']
Dsparse = as(ifelse(is.na(D),0,D),"dsCMatrix")

sparse_exp = function(A){
  B = exp(as.matrix(A))
  B = ifelse( B==1, 0, B )
  return( as(B,"dsCMatrix") )
}
exp_neg_D = sparse_exp(-D)

# Parameter values
Rho=0.3           # Correlation at distance = 1
SDmarg=0.3        # Marginal SD of epsilon_b
theta = -log(Rho)  # Decorrelation rate per distance
log_mean = 4

# Derived
SDinput = SDmarg * sqrt(2*theta)
SDcond = SDmarg * sqrt(1-Rho^2)

# Simulate random effect
years <- 1:10
pow = function(a,b) a^b
Type = function(a) a
# epsilon_x = rep(NA, length(Obs_i))
epsilon_yx <- matrix(NA, nrow=length(years), ncol=length(Obs_i))
# rho_b = SDinput_b = rep(NA, length(Obs_i))
rho_yb <- SDinput_yb <- matrix(NA, nrow=length(years), ncol=length(Obs_i))
for(y in 1:length(years)){
  for( b in 1:length(Obs_i) ){
    if( is.na(dist_x[b]) || dist_x[b]==Inf ){
      # Correlation between i and parent(i) as distance -> INF
      rho_yb[y,b] = 0;
      # SD of Ornstein-Uhlenbeck process as distance -> INF
      SDinput_yb[y,b] = SDinput / pow(2*theta, 0.5);
      # conditional probability
      epsilon_yx[y,b] = rnorm(n=1, Type(0.0), SDinput_yb[y,b] );
    }else{
      # Correlation between i and parent(i)
      rho_yb[y,b] = exp(-theta * dist_x[b]);
      # SD of O-U process
      SDinput_yb[y,b] = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_x[b])), 0.5 );
      # conditional probability
      epsilon_yx[y,b] = rnorm(n=1, rho_yb[y,b]*epsilon_yx[y,parent_x[b]], SDinput_yb[y,b] );
    }
  }
}

# Simulate data
# lambda_i = exp( epsilon_x + log_mean )
lambda_yi <- exp(epsilon_yx + log_mean)

## catch
# catch_i = rpois( n=length(lambda_i), lambda=lambda_i )
# catch_i[4] = 0
catch_yi <- t(sapply(1:length(years), function(x) rpois( n = length(lambda_yi[x,]), lambda=lambda_yi[x,])))
setzero <- sample(1:length(Obs_i), length(years))
for(y in 1:length(years)){
  catch_yi[y,setzero[y]] <- 0
}
dfcatch <- lapply(1:length(years), function(y){
  df <- data.frame("year" = years[y], "child_x"=Obs_i, "catch"=catch_yi[y,])
  return(df)
})
dfcatch <- do.call(rbind, dfcatch)

## presence/absence
# pres_i <- sapply(1:length(catch_i), function(x) ifelse(catch_i[x] > 0, 1, 0))
# dev <- rnorm(length(catch_i), 0, 0.0001)
# pres_i_dev <- pres_i + dev
# pres_i_dev[which(pres_i == 0)] <- 0
pres_yi <- t(sapply(1:length(years), function(y){
    sapply(1:length(Obs_i), function(x){
      out <- ifelse(catch_yi[y,x] > 0, 1, 0)
      return(out)
    })
}))
dev <- matrix(rnorm(length(Obs_i)*length(years), 0, 0.0001), nrow=length(years), ncol=length(Obs_i))
pres_yi_dev <- pres_yi + dev
pres_yi_dev[which(pres_yi == 0)] <- 0
dfpres <- lapply(1:length(years), function(y){
  df <- data.frame("year" = years[y], "child_x"=Obs_i, "catch"=pres_yi_dev[y,])
  return(df)
})
dfpres <- do.call(rbind, dfpres)


dfnet <- cbind.data.frame(Network_xn, loc_LL)

# plotdf <- cbind.data.frame(Network_xn, loc_LL, catch_i, pres_i, pres_i_dev)
dfnet$Lon2 <- sapply(1:nrow(dfnet), function(x) ifelse(parent_x[x]==0, NA, dfnet$Lon[which(dfnet$child_x == dfnet$parent_x[x])]))
dfnet$Lat2 <- sapply(1:nrow(dfnet), function(x) ifelse(parent_x[x]==0, NA, dfnet$Lat[which(dfnet$child_x == dfnet$parent_x[x])]))

dfobs <- inner_join(dfcatch, dfnet)
dfpres <- inner_join(dfpres, dfnet)


## sampling locations
samp <- sapply(1:length(years), function(x){
  out <- sample(1:length(Obs_i), round(length(Obs_i)*0.75))
  return(out)
})

dfcatch_samp <- lapply(1:length(years), function(x){
  sub <- dfcatch %>% filter(year == years[x])
  sample <- sapply(1:nrow(sub), function(y) ifelse(sub$child_x[y] %in% samp[,x], TRUE, FALSE))
  sub$Sampled <- sample
  catch <- sapply(1:nrow(sub), function(y) ifelse(sub$Sampled[y] == TRUE, sub$catch[y], NA))
  sub$catch <- catch
  return(sub)
})
dfcatch_samp <- do.call(rbind, dfcatch_samp)

dfpres_samp <- lapply(1:length(years), function(x){
  sub <- dfpres %>% filter(year == years[x])
  sample <- sapply(1:nrow(sub), function(y) ifelse(sub$child_x[y] %in% samp[,x], TRUE, FALSE))
  sub$Sampled <- sample
  catch <- sapply(1:nrow(sub), function(y) ifelse(sub$Sampled[y] == TRUE, sub$catch[y], NA))
  sub$catch <- catch
  return(sub)
})
dfpres_samp <- do.call(rbind, dfpres_samp)

dfobs_samp <- inner_join(dfcatch_samp, dfnet)
dfpres_samp <- inner_join(dfpres_samp, dfnet)



dfobs_samp_toUse <- dfobs_samp %>% select("year","catch","Lon","Lat") %>% na.omit()
dfpres_samp_toUse <- dfpres_samp %>% select("year","catch", "Lon","Lat") %>% na.omit()
dfobs_toUse <- dfobs %>% select("year","catch","Lon","Lat")
dfpres_toUse <- dfpres %>% select("year","catch","Lon","Lat")
dfnet_toUse <- dfnet %>% select("parent_x","child_x","dist_x","Lon","Lat")

saveRDS(dfnet_toUse, file.path(data_dir, "Simulated_network.rds"))
saveRDS(dfobs_toUse, file.path(data_dir, "Simulated_observations_catch.rds"))
saveRDS(dfpres_toUse, file.path(data_dir, "Simulated_observations_encounters.rds"))
saveRDS(dfobs_samp_toUse, file.path(data_dir, "Simulated_observations_catch_sampled.rds"))
saveRDS(dfpres_samp_toUse, file.path(data_dir, "Simulated_observations_encounters_sampled.rds"))


aa <- ggplot(dfnet) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_x), fill = "black", color = "white") + #, nudge_x = -0.02) + 
    # guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()
ggsave(file.path(fig_dir, "Simulated_network.png"), aa)

a <- ggplot(dfobs) +
    facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = catch), color = "white") + #, nudge_x = -0.02) + 
    guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()
ggsave(file.path(fig_dir, "Simulated_catch.png"), a)

a1 <- ggplot(dfobs %>% filter(year==1)) +
    # facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = catch), color = "white") + #, nudge_x = -0.02) + 
    guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()
ggsave(file.path(fig_dir, "Simulated_catch_oneYear.png"), a1)


b <- ggplot(dfpres) +
    facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = factor(round(catch))), color = "white") + #, nudge_x = -0.02) + 
    # geom_point(aes(x = Lon, y = Lat, color = factor(pres_i)), cex=5) +
    guides(fill = guide_legend(title = "Presence")) +
    mytheme()
ggsave(file.path(fig_dir, "Simulated_encounters.png"), b)

b1 <- ggplot(dfpres %>% filter(year==1)) +
    # facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = factor(round(catch))), color = "white") + #, nudge_x = -0.02) + 
    # geom_point(aes(x = Lon, y = Lat, color = factor(pres_i)), cex=5) +
    guides(fill = guide_legend(title = "Presence")) +
    mytheme()
ggsave(file.path(fig_dir, "Simulated_encounters_oneYear.png"), b1)


a2 <- ggplot(dfobs_samp) +
    facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = catch), color = "white") + #, nudge_x = -0.02) + 
    guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()
ggsave(file.path(fig_dir, "Simulated_catch_sampled.png"), a2)

a21 <- ggplot(dfobs_samp %>% filter(year==1)) +
    # facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = catch), color = "white") + #, nudge_x = -0.02) + 
    guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()
ggsave(file.path(fig_dir, "Simulated_catch_sampled_oneYear.png"), a21)

b2 <- ggplot(dfpres_samp) +
    facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = factor(round(catch))), color = "white") + #, nudge_x = -0.02) + 
    # geom_point(aes(x = Lon, y = Lat, color = factor(pres_i)), cex=5) +
    guides(fill = guide_legend(title = "Presence")) +
    mytheme()
ggsave(file.path(fig_dir, "Simulated_encounters_sampled.png"), b2)

b21 <- ggplot(dfpres_samp %>% filter(year==1)) +
    # facet_wrap(~year) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    geom_label(aes(x = Lon, y = Lat, label = child_x, fill = factor(round(catch))), color = "white") + #, nudge_x = -0.02) + 
    # geom_point(aes(x = Lon, y = Lat, color = factor(pres_i)), cex=5) +
    guides(fill = guide_legend(title = "Presence")) +
    mytheme()
ggsave(file.path(fig_dir, "Simulated_encounters_sampled_oneYear.png"), b21)



########################
# Fit in VAST -- Prsence/absence
########################
# Save
DateFile = paste0(wd,'/simulation/Sim_Stream_Presence/')
dir.create(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Mesh", "Stream_network")[2]
grid_size_km = 1
n_x = nrow(Network_xn)   # Specify number of stations (a.k.a. "knots")
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)

RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = c(1,1)

# Default
strata.limits <- data.frame('STRATA'="All_areas")

# Load data
Data_Geostat = data.frame( "Catch_KG"=dfpres_toUse$catch, "Year"=dfpres_toUse$year, "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=dfpres_toUse$Lat, "Lon"=dfpres_toUse$Lon, "Pass"=0)
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="User", input_grid=cbind(dfpres_toUse %>% select("Lat","Lon"),"Area_km2"=1), strata.limits=strata.limits )
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Run first time
Data_orig = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "c_iz"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method, "Options"=Options, "Network_sz"=Network_sz )

TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList_orig[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt_orig = tryCatch(TMBhelper::Optimize( obj=Obj, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=FALSE, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") ), error=function(e) NA)

i <- 0
while(any(grepl("opt", names(Opt_orig))) | all(is.na(Opt_orig))){
  i <- i + 1
  TmbList = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Opt_orig = tryCatch(TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") ), error=function(e) NA)
  # names <- names(Opt_orig)
}

Opt_orig = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=5, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )


Report_orig = Obj$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj$env$parList(Opt_orig$par), "Data"=Data_orig)
save(Save_orig, file=paste0(DateFile,"Save_orig.RData"))

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt_orig$SD

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile, pch=19 )

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report_orig, Data_Geostat=Data_Geostat, DirName=DateFile)

## Model selection
Opt_orig$AIC






################## 
## model 2
##################
DateFile = paste0(wd,'/simulation/Sim_Stream_Catch_Spatial/')
dir.create(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Grid", "Mesh", "Stream_network")[3]
grid_size_km = 1
n_x = nrow(Network_xn)   # Specify number of stations (a.k.a. "knots")


## Omega1 = Spatial variation in encounter probability
## Epsilon1 = Spatiotemporal variation in encounter probability
## Omega2 = Spatial variation in positive catch rates
## Epislon2 = Spatiotemporal variation in positive catch rates
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)

## Optional vector
## Eta1 = Overdispersion in encounter probability
## Eta2 = Overdispersion in positive catch rates
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

## Optional vector
## Calculate_Range = 1 calculates center of gravity
## Calculate_effective_area = 1 calculates effective area occupied
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)

## Optional vector
## Specifies how temporal intercepts (Beta) or
## spatiotemporal intercepts (Epsilon) are structured
## 0 = Each year is fixed effect
## 1 = Each year is random effect
## 2 = Each year follows random walk
## 3 = Constant among years as fixed effect
## 4 = Each year is random following AR1 process
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## [1] Specifies distribution of positive catch rates
## [2] Specifies functional form for encounter probabilities
## example: [1] 7 = Poisson, [2] 0 = Delta model with logit-link
ObsModel = c(2,1)

# Default
strata.limits <- data.frame('STRATA'="All_areas")

# Load data
Data_Geostat = data.frame( "Catch_KG"=dfobs$catch, "Year"=dfobs$year, "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=dfobs$Lat, "Lon"=dfobs$Lon, "Pass"=0)
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="User", input_grid=cbind(loc_LL,"Area_km2"=1), strata.limits=strata.limits )
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Run first time
Network_sz = cbind( "parent_s"=Network_xn[,'parent_x'], "child_s"=Network_xn[,'child_x'], "dist_s"=Network_xn[,'dist_x'] )

Data_orig = Data_Fn("Version"=Version, 
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

TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList_orig[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt_orig = TMBhelper::Optimize( obj=Obj, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report_orig = Obj$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj$env$parList(Opt_orig$par), "Data"=Data_orig)
save(Save_orig, file=paste0(DateFile,"Save_orig.RData"))

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt_orig$SD

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile, pch=19 )

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report_orig, Data_Geostat=Data_Geostat, DirName=DateFile)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data_orig, Report=Report_orig, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile )

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data_orig, Report=Report_orig, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

## Model selection
Opt_orig$AIC

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report_orig, Sdreport=Opt_orig$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( DirName=DateFile, TmbData=Data_orig, Sdreport=Opt_orig$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )

################## 
## model 3
##################
DateFile = paste0(wd,'/simulation/Sim_Stream_Catch_Temporal/')
dir.create(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Grid", "Mesh", "Stream_network")[3]
grid_size_km = 1
n_x = nrow(Network_xn)   # Specify number of stations (a.k.a. "knots")


## Omega1 = Spatial variation in encounter probability
## Epsilon1 = Spatiotemporal variation in encounter probability
## Omega2 = Spatial variation in positive catch rates
## Epislon2 = Spatiotemporal variation in positive catch rates
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)

## Optional vector
## Eta1 = Overdispersion in encounter probability
## Eta2 = Overdispersion in positive catch rates
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

## Optional vector
## Calculate_Range = 1 calculates center of gravity
## Calculate_effective_area = 1 calculates effective area occupied
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)

## Optional vector
## Specifies how temporal intercepts (Beta) or
## spatiotemporal intercepts (Epsilon) are structured
## 0 = Each year is fixed effect
## 1 = Each year is random effect
## 2 = Each year follows random walk
## 3 = Constant among years as fixed effect
## 4 = Each year is random following AR1 process
RhoConfig = c("Beta1"=1, "Beta2"=1, "Epsilon1"=0, "Epsilon2"=0)

## [1] Specifies distribution of positive catch rates
## [2] Specifies functional form for encounter probabilities
## example: [1] 7 = Poisson, [2] 0 = Delta model with logit-link
ObsModel = c(2,1)

# Default
strata.limits <- data.frame('STRATA'="All_areas")

# Load data
Data_Geostat = data.frame( "Catch_KG"=dfobs$catch, "Year"=dfobs$year, "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=dfobs$Lat, "Lon"=dfobs$Lon, "Pass"=0)
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="User", input_grid=cbind(loc_LL,"Area_km2"=1), strata.limits=strata.limits )
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Run first time
Network_sz = cbind( "parent_s"=Network_xn[,'parent_x'], "child_s"=Network_xn[,'child_x'], "dist_s"=Network_xn[,'dist_x'] )

Data_orig = Data_Fn("Version"=Version, 
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

TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj = TmbList_orig[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj)

Opt_orig = TMBhelper::Optimize( obj=Obj, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )

Report_orig = Obj$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj$env$parList(Opt_orig$par), "Data"=Data_orig)
save(Save_orig, file=paste0(DateFile,"Save_orig.RData"))

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]
Opt_orig$SD

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile, pch=19 )

## convergence
Opt_orig$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]

## diagnostics for encounter probability component
Enc_prob = plot_encounter_diagnostic( Report=Report_orig, Data_Geostat=Data_Geostat, DirName=DateFile)

## diagnostics for positive catch rate component
Q = plot_quantile_diagnostic( TmbData=Data_orig, Report=Report_orig, FileName_PP="Posterior_Predictive",
  FileName_Phist="Posterior_Predictive-Histogram", 
  FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile )

## Diagnostics for plotting residuals on a map
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"="Sim", "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

# Decide which years to plot                                                   
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

## Plot Pearson residuals
plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Data_orig, Report=Report_orig, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=ModFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=3, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

## Model selection
Opt_orig$AIC

##### Model output
## density surface for each year
Dens_xt = plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report_orig, Sdreport=Opt_orig$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=1.5, pch=19, Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), plot_legend_fig=TRUE)

## density predictions at different locations - output in UTM
Dens_DF = cbind( "Density"=as.vector(Dens_xt), "Year"=Year_Set[col(Dens_xt)], "E_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'E_km'], "N_km"=Spatial_List$MeshList$loc_x[row(Dens_xt),'N_km'] )

## index of abundance
Index = plot_biomass_index( DirName=DateFile, TmbData=Data_orig, Sdreport=Opt_orig$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE )



