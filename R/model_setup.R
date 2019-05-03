
model_setup <- function(model_info, Method, grid_size_km = 1){

	FieldConfig_temp <- c(0,0,0,0)
	FieldConfig_st <- c("IID",0,0,0)
	FieldConfig_sst <- c("IID","IID",0,0)	

	RhoConfig1_iid <- c(1,3,0,0)
	RhoConfig1_rw <- c(2,3,0,0)
	RhoConfig2_iid <- c(1,3,1,0)
	RhoConfig2_rw <- c(2,3,2,0)

  if(model_info[1] == "T"){ FieldConfig_inp <- FieldConfig_temp }
  if(model_info[1] == "ST"){ FieldConfig_inp <- FieldConfig_st }
  if(model_info[1] == "SST"){ FieldConfig_inp <- FieldConfig_sst }

  if(model_info[2] == "IID"){
    if(model_info[1] == "SST"){
      RhoConfig_inp <- RhoConfig2_iid
    } else { RhoConfig_inp <- RhoConfig1_iid }
  }
  if(model_info[2] == "RW"){
    if(model_info[1] == "SST"){
      RhoConfig_inp <- RhoConfig2_rw
    } else { RhoConfig_inp <- RhoConfig1_rw }
  }

  OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
  Options =  c("Calculate_Range"=1, 
              "Calculate_effective_area"=1)
  ObsModel = c("PosDist"=1, "Link"=0)

  settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig_inp, RhoConfig=RhoConfig_inp, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
  settings$Method <- "Stream_network"
  settings$grid_size_km <- 1

  return(settings)

}