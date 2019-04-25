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
library(RuddR)
library(caret)


#########################
## read in data
##########################

data <- data("nz_waitaki_longfin_eel", package="StreamUtils")

network <- nz_waitaki_longfin_eel[["network"]]
obs <- nz_waitaki_longfin_eel[["observations"]]
hab <- nz_waitaki_longfin_eel[["habitat"]]

Network_sz = network %>% select(-c("long","lat"))

Network_sz$dist_s <- Network_sz$dist_s/1000
obs$dist_i <- obs$dist_i/1000

##################################
## model - encounter observations
##################################
nz_enc_dir <- file.path(nz_dir, "waitaki_encounters")
dir.create(nz_enc_dir, showWarnings=FALSE)
setwd(nz_enc_dir)

AIC_list <- NULL

##### General settings
Version = "VAST_v5_3_0" # SpatialDeltaGLMM::get_latest_version( package="VAST" )
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
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = obs$year,
               "Vessel" = "missing", 
               "AreaSwept_km2" = 1, 
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
                          Lon_i=Data_Geostat[,'Lon'], 
                          Lat_i=Data_Geostat[,'Lat'], 
                          Lat_x=network$lat, 
                          Lon_x=network$long, 
                          Extrapolation_List=Extrapolation_List )

## add locations to dataset
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

####################################
## habitat information
## treated as density covariates
###################################

nodes <- Network_sz$child_s[order(Network_sz$child_s)]
years <- min(Data_Geostat$Year):max(Data_Geostat$Year)
covar <- unique(hab$covariate)

df <- lapply(1:length(covar), function(x){
  sub <- hab %>% filter(covariate == covar[x])
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
pairs(df, upper.panel=panel.cor)

df2 = cor(df)
fc <- findCorrelation(df2, cutoff=0.5) # 
fc <- sort(fc)
reduced_Data = df[,-c(fc)]
pairs(reduced_Data, upper.panel=panel.cor)

covar_toUse <- colnames(reduced_Data)
hab_toUse <- hab %>% filter(covariate %in% covar_toUse)

n_t <- length(years)
n_x <- length(nodes)
n_p <- length(covar_toUse)
X_xtp <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:n_p){
    psub <- hab_toUse %>% filter(covariate == covar_toUse[p])
    mat <- matrix(NA, nrow=n_x, ncol=n_t)
    for(t in 1:n_t){
      tsub <- psub %>% filter(year==years[t])
      mat[tsub$child_s,t] <- tsub$value
    }
    mat_std <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
    for(t in 1:n_t){
      sub <- mat_std[,t]
      sub[which(is.na(sub))] <- 0
      mat_std[,t] <- sub
    }
    X_xtp[,,p] <- mat_std
}

# Xconfig_zcp <- array(0, dim=c(2,1,n_p))

##################
## SPATIAL GAMMA HABITAT
##################
space_dir <- file.path(nz_enc_dir, "spatial_gamma_habitat")
dir.create(space_dir)
setwd(space_dir)

FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=0, 
            "Calculate_effective_area"=0)
ObsModel = c(2,0)

new_dir <- file.path(space_dir, "all")
dir.create(new_dir, showWarnings=FALSE)
setwd(new_dir)

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
                  "X_xtp"=X_xtp,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

  TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s'])) 

  Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  check <- TMBhelper::Check_Identifiable(Obj) 

  Opt = TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  
  if(all(is.na(Opt))==FALSE){
    AIC_list[[p]] <- as.numeric(Opt$AIC)  

    #############
    Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]   

    Opt[["SD"]]   

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
    Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001)   

    # ## index of abundance
    Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )


for(p in 1:n_p){
  new_dir <- file.path(space_dir, covar_toUse[p])
  dir.create(new_dir, showWarnings=FALSE)
  setwd(new_dir)

  if(file.exists("Save.rds")) next

  X_inp <- X_xtp[,,p]
  X_inp <- array(X_inp, dim=c(dim(X_inp),1))

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
                  "X_xtp"=X_inp,
                  "MeshList"=Spatial_List$MeshList,
                  "GridList"=Spatial_List$GridList,
                  "Method"=Spatial_List$Method,
                  "Options"=Options,
                  "Network_sz"=Network_sz,
                  "CheckForErrors"=FALSE )

  TmbList = Build_TMB_Fn("TmbData"=Data, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj = TmbList[["Obj"]]
  Obj$par[grep("logkappa",names(Obj$par))] = log(1/median(Network_sz[,'dist_s'])) 

  Opt1 = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=FALSE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  check <- TMBhelper::Check_Identifiable(Obj) 

  Opt = tryCatch(TMBhelper::Optimize( obj=Obj, startpar=Opt1$par, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, bias.correct=TRUE, newtonsteps=3, bias.correct.control=list(sd=TRUE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") ), error=function(e) NA)
  
  if(all(is.na(Opt))==FALSE){
    AIC_list[[p]] <- as.numeric(Opt$AIC)  

    #############
    Opt$diagnostics[,c('Param','Lower','MLE','Upper','final_gradient')]   

    Opt[["SD"]]   

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
    Cov_xt <- StreamUtils::plot_maps(plot_set=c(11), TmbData=Data, Report=Report, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, Panel="year", cex=0.0001)   

    # ## index of abundance
    Index = StreamUtils::plot_biomass_index( TmbData=Data, Sdreport=Opt$SD, Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=FALSE, strata_names = "Longfin eels" )
  } else {
    Save <- list("TmbList"=TmbList, "Obj"=Obj, "Opt"=Opt, "Report"=NA)
    saveRDS(Save, file.path(getwd(), "Save.rds"))
  }
}

aic <- sapply(1:n_p, function(p){
  read <- readRDS(file.path(space_dir, covar_toUse[p]), "Save.rds")
}
