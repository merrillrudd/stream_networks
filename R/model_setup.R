
model_setup <- function(model_info, Method, grid_size_km = 1, type, ncategories = 0, Options = c("Calculate_Range"=1, "Calculate_effective_area" = 1), link=0){

	FieldConfig_temp <- c(0,0,0,0)
	if(type == "encounter"){
    if(any(grepl("FieldIID",model_info))){
      FieldConfig_st <- c("IID",0,0,0)
	    FieldConfig_sst <- c("IID","IID",0,0)	
    }
    if(any(grepl("FieldFACTOR",model_info))){
      FieldConfig_st <- c(ncategories,0,0,0)
      FieldConfig_sst <- c(ncategories,ncategories,0,0)
    }
  }
  if(type != "encounter"){
    if(any(grepl("FieldIID",model_info))){
      FieldConfig_st <- c("IID", 0, "IID", 0)
      FieldConfig_sst <- rep("IID", 4)
    }
    if(any(grepl("FieldFACTOR",model_info))){
      FieldConfig_st <- rep(c(ncategories,0),2)
      FieldConfig_sst <- rep(ncategories,4)
    }
  }

  if(type == "encounter"){
    RhoConfig0 <- c(0,3,0,0)
    RhoConfig1_iid <- c(1,3,0,0)
    RhoConfig1_rw <- c(2,3,0,0)
    RhoConfig2_iid <- c(1,3,1,0)
    RhoConfig2_rw <- c(2,3,2,0)
  }
  if(type != "density"){
    RhoConfig0 <- c(0,0,0,0)
    RhoConfig1_iid <- c(1,1,0,0)
    RhoConfig1_rw <- c(2,2,0,0)
    RhoConfig2_iid <- c(1,1,1,1)
    RhoConfig2_rw <- c(2,2,2,2)
  }

  if(model_info[1] == "T"){ FieldConfig_inp <- FieldConfig_temp }
  if(model_info[1] == "ST"){ FieldConfig_inp <- FieldConfig_st }
  if(model_info[1] == "SST"){ FieldConfig_inp <- FieldConfig_sst }

  if(any(grepl("SmoothIID",model_info))){
    if(model_info[1] == "SST"){
      RhoConfig_inp <- RhoConfig2_iid
    } else { RhoConfig_inp <- RhoConfig1_iid }
  }
  if(any(grepl("SmoothRW",model_info))){
    if(model_info[1] == "SST"){
      RhoConfig_inp <- RhoConfig2_rw
    } else { RhoConfig_inp <- RhoConfig1_rw }
  }
  if(link==3) RhoConfig_inp <- RhoConfig0

  OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

if(any(grepl("LOGNORMAL",path))) ObsModel1 <- 1
if(any(grepl("GAMMA",path))) ObsModel1 <- 2
if(any(grepl("NORMAL",path))) ObsModel1 <- 0
ObsModel_inp <- c(ObsModel1, link)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel_inp, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  return(settings)

}