#' @title
#' Plot stream network and optional observations
#' 
#' @description
#' \code{plot_network} generates figure showing network nodes, directional arrows, and observations
#' 
#' @param Network_sz_LL
#' @param Data_Geostat default TRUE to include observations on the figure, FALSE to show only network
#' @param root default TRUE to show root nodes, FALSE in case there are other root nodes that are not meaningful.
#' @param FileName Directory (absolute path) and base for filenames of plots
#' @return Figure plotting stream network and observations
#' @export
plot_network <- function(Network_sz_LL, Data_Geostat,  FileName, root = FALSE){

	### across years
	aa <- ggplot(Network_sz_LL) +
			mytheme()

	## add roots underneath points
	if(root == TRUE) aa <- aa + geom_point(data = Network_sz_LL %>% filter(parent_s==0), aes(x = Lon, y = Lat), color="goldenrod", cex=5)

	## option to add arrows
	# if(arrows == TRUE){
		l2 <- lapply(1:nrow(Network_sz_LL), function(x){
			parent <- Network_sz_LL$parent_s[x]
			find <- Network_sz_LL %>% filter(child_s == parent)
			if(nrow(find)>0) out <- cbind.data.frame(Network_sz_LL[x,], 'Lon2'=find$Lon, 'Lat2'=find$Lat)
			if(nrow(find)==0) out <- cbind.data.frame(Network_sz_LL[x,], 'Lon2'=NA, 'Lat2'=NA)
			# if(nrow(find)>0) out <- cbind.data.frame(Network_sz_LL[x,], 'long2'=find$long, 'lat2'=find$lat)
			# if(nrow(find)==0) out <- cbind.data.frame(Network_sz_LL[x,], 'long2'=NA, 'lat2'=NA)
			return(out)
		})
		l2 <- do.call(rbind, l2)
		aa <- aa + geom_segment(data=l2, aes(x = Lon2,y = Lat2, xend = Lon, yend = Lat), col="gray")
	# }

	## add points and complete figure
	aa <- aa +
	    geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", alpha=0.6) +
	    xlab("Longitude") + ylab("Latitude")

	## option to add observations
	if(all(is.null(Data_Geostat))==FALSE){
		 aa <- aa + 
				geom_point(data = Data_Geostat, aes(x = Lon, y = Lat, fill=Category), cex=2, pch=22) + 
				scale_fill_brewer(palette = "Set1")
	}
	if(length(unique(Data_Geostat$Category))==1){
		aa <- aa + guides(fill = FALSE)
		width <- 9
		height <- 8
	} else {
		width = 10
		height = 8
	}
	ggsave(file.path(FileName, "Network.png"), aa, width = width, height = height)


	### by year
	years <- unique(Data_Geostat$Year)[order(unique(Data_Geostat$Year))]

	Network_sz_LL_wYear <- lapply(1:length(years), function(x){
		out <- cbind.data.frame(Network_sz_LL, "Year"=years[x])
		return(out)
	})
	Network_sz_LL_wYear <- do.call(rbind, Network_sz_LL_wYear)

	bb <- ggplot(Network_sz_LL_wYear) +
			facet_wrap(~Year) +
			mytheme() 

	## add roots underneath points
	if(root == TRUE) bb <- bb + geom_point(data = Network_sz_LL %>% filter(parent_s==0), aes(x = Lon, y = Lat), color="goldenrod", cex=5)

	## add points and complete figure
	bb <- bb +
	    geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", cex=0.5, alpha=0.6) +
	    xlab("Longitude") + ylab("Latitude")

	## option to add observations
	## option to add observations
		 bb <- bb + 
				geom_point(data = Data_Geostat, aes(x = Lon, y = Lat, fill=Category), cex=1.8, pch=22, alpha=0.6) + 
				scale_fill_brewer(palette = "Set1") +
				scale_x_continuous(breaks=round(quantile(Data_Geostat$Lon,prob=c(0.2,0.5,0.8)),0), labels=round(quantile(Data_Geostat$Lon,prob=c(0.2,0.5,0.8)),0)) 
	if(length(unique(Data_Geostat$Category))==1){
		bb <- bb + guides(fill = FALSE)
		width <- 8
		height = 8
	} else {
		width = 10
		height = 8
	}

	ggsave(file.path(FileName, "Network_byYear.png"), bb, width = width, height = height)

	
}