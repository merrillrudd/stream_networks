rm(list=ls())

# devtools::install_github("james-thorson/VAST", ref = "development")
library(VAST)
library(TMB)
library(tidyverse)
library(RuddR)


wd <- "C:\\merrill\\stream_networks"
setwd(wd)

##################################
## simulated stream network data
##################################
# Network parentage
set.seed(1)
parent_x = c( 0, 1, 2, 3, 3, 5, 5, 5, 8, 9, 10 )
dist_x = c( Inf, rlnorm(length(parent_x)-1, meanlog=0, sdlog=1) )
Network_xn = cbind( parent_x=parent_x, child_b=1:length(parent_x), dist_x=dist_x )
Obs_i = 1:nrow(Network_xn)
loc_xy = matrix(c( 0,0, 1,0, 2,0, 2,-1, 2,1, 1,1, 2,2, 3,1, 3,0, 3,-1, 3,-2), byrow=TRUE, ncol=2, dimnames=list(NULL,c("x","y")) )
loc_LL = matrix(loc_xy*0.1 + outer(rep(1,nrow(loc_xy)),c(20,40)), ncol=2, dimnames=list(NULL,c("Lon","Lat")))
xlab <- c(NA, sapply(2:nrow(loc_xy), function(x) mean(c(loc_xy[parent_x[x],'x'],loc_xy[x,'x']))))
ylab <- c(NA, sapply(2:nrow(loc_xy), function(x) mean(c(loc_xy[parent_x[x],'y'],loc_xy[x,'y']))))

# Distance matrix
D = array(NA, dim=rep(length(dist_x),2) )
D[ Network_xn[-1,c("parent_x","child_b")] ] = Network_xn[-1,'dist_x']
D[ Network_xn[-1,c("child_b","parent_x")] ] = Network_xn[-1,'dist_x']
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
epsilon_x = rep(NA, length(Obs_i))
pow = function(a,b) a^b
Type = function(a) a
rho_b = SDinput_b = rep(NA, length(Obs_i))
for( b in 1:length(epsilon_x) ){
  if( is.na(dist_x[b]) || dist_x[b]==Inf ){
    # Correlation between i and parent(i) as distance -> INF
    rho_b[b] = 0;
    # SD of Ornstein-Uhlenbeck process as distance -> INF
    SDinput_b[b] = SDinput / pow(2*theta, 0.5);
    # conditional probability
    epsilon_x[b] = rnorm(n=1, Type(0.0), SDinput_b[b] );
  }else{
    # Correlation between i and parent(i)
    rho_b[b] = exp(-theta * dist_x[b]);
    # SD of O-U process
    SDinput_b[b] = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_x[b])), 0.5 );
    # conditional probability
    epsilon_x[b] = rnorm(n=1, rho_b[b]*epsilon_x[parent_x[b]], SDinput_b[b] );
  }
}

# Simulate data
lambda_i = exp( epsilon_x + log_mean )

## catch
catch_i = rpois( n=length(lambda_i), lambda=lambda_i )
catch_i[4] = 0

## presence/absence
pres_i <- sapply(1:length(catch_i), function(x) ifelse(catch_i[x] > 0, 1, 0))
dev <- rnorm(length(catch_i), 0, 0.0001)
pres_i_dev <- pres_i + dev
pres_i_dev[which(pres_i == 0)] <- 0

plotdf <- cbind.data.frame(Network_xn, loc_LL, catch_i, pres_i, pres_i_dev)
plotdf$Lon2 <- sapply(1:nrow(plotdf), function(x) ifelse(parent_x[x]==0, NA, plotdf$Lon[which(plotdf$child_b == plotdf$parent_x[x])]))
plotdf$Lat2 <- sapply(1:nrow(plotdf), function(x) ifelse(parent_x[x]==0, NA, plotdf$Lat[which(plotdf$child_b == plotdf$parent_x[x])]))

aa <- ggplot(plotdf) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_b), fill = "black", color = "white") + #, nudge_x = -0.02) + 
    # guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()

a <- ggplot(plotdf) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    # geom_point(aes(x = Lon, y = Lat, color = catch_i), cex=5) +
    geom_label(aes(x = Lon, y = Lat, label = child_b, fill = catch_i), color = "white") + #, nudge_x = -0.02) + 
    guides(fill = guide_legend(title = "Catch")) +
    xlab("Longitude") + ylab("Latitude") + 
    mytheme()

b <- ggplot(plotdf) +
    geom_segment(aes(x = Lon,y = Lat,xend = Lon2,yend = Lat2),arrow=arrow()) +
    geom_label(aes(x = Lon, y = Lat, label = child_b, fill = factor(pres_i)), color = "white") + #, nudge_x = -0.02) + 
    # geom_point(aes(x = Lon, y = Lat, color = factor(pres_i)), cex=5) +
    guides(fill = guide_legend(title = "Presence")) +
    mytheme()


########################
# Fit in VAST -- Catch data
########################

# Save
DateFile = paste0(wd,'/examples/Sim_Stream_Catch/')
dir.create(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Grid", "Mesh", "Stream_network")[3]
grid_size_km = 1
n_x = nrow(Network_xn)   # Specify number of stations (a.k.a. "knots")


## Omega1 = Spatial variation in encounter probability
## Epsilon1 = Spatiotemporal variation in encounter probability
## Omega2 = Spatial variation in positive catch rates
## Epislon2 = Spatiotemporal variation in positive catch rates
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=0)

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
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)

## [1] Specifies distribution of positive catch rates
## [2] Specifies functional form for encounter probabilities
## example: [1] 7 = Poisson, [2] 0 = Delta model with logit-link
ObsModel = c(7,0)

# Default
strata.limits <- data.frame('STRATA'="All_areas")

# Load data
Data_Geostat = data.frame( "Catch_KG"=catch_i, "Year"=2011, "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=loc_LL[,'Lat'], "Lon"=loc_LL[,'Lon'], "Pass"=0)
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="User", input_grid=cbind(loc_LL,"Area_km2"=1), strata.limits=strata.limits )
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Run first time
Network_sz = cbind( "parent_s"=Network_xn[,'parent_x'], "child_s"=Network_xn[,'child_b'], "dist_s"=Network_xn[,'dist_x'] )

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
Obj_orig = TmbList_orig[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj_orig)

Opt_orig = TMBhelper::Optimize( obj=Obj_orig, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Report_orig = Obj_orig$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj_orig$env$parList(Opt_orig$par), "Data"=Data_orig)
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

########################
# Fit in VAST -- Prsence/absence
########################
# Save
DateFile = paste0(wd,'/examples/Sim_Stream_Presence/')
dir.create(DateFile)

Version = "VAST_v5_2_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
Method = c("Mesh", "Stream_network")[2]
grid_size_km = 1
n_x = nrow(Network_xn)   # Specify number of stations (a.k.a. "knots")
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, "Calculate_effective_area"=0)

RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = c(2,3)

# Default
strata.limits <- data.frame('STRATA'="All_areas")

# Load data
Data_Geostat = data.frame( "Catch_KG"=pres_i_dev, "Year"=2011, "Vessel"="missing", "AreaSwept_km2"=1, "Lat"=loc_LL[,'Lat'], "Lon"=loc_LL[,'Lon'], "Pass"=0)
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region="User", input_grid=cbind(loc_LL,"Area_km2"=1), strata.limits=strata.limits )
Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, DirPath=DateFile, Save_Results=TRUE )
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# Run first time
Network_sz = cbind( "parent_s"=Network_xn[,'parent_x'], "child_s"=Network_xn[,'child_b'], "dist_s"=Network_xn[,'dist_x'] )
Data_orig = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "c_iz"=rep(0,nrow(Data_Geostat)), "pres_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method, "Options"=Options, "Network_sz"=Network_sz )
TmbList_orig = Build_TMB_Fn("TmbData"=Data_orig, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
Obj_orig = TmbList_orig[["Obj"]]
check <- TMBhelper::Check_Identifiable(Obj_orig)

Opt_orig = TMBhelper::Optimize( obj=Obj_orig, lower=TmbList_orig[["Lower"]], upper=TmbList_orig[["Upper"]], getsd=TRUE, savedir=DateFile, newtonsteps=5, bias.correct=FALSE, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
Report_orig = Obj_orig$report()
Save_orig = list("Opt"=Opt_orig, "Report"=Report_orig, "ParHat"=Obj_orig$env$parList(Opt_orig$par), "Data"=Data_orig)
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
