# Some parts of code adapted from https://philippgaertner.github.io/2019/12/earth-engine-rstudio-reticulate/
library(reticulate)
library(sf)
library(raster)
library(stars)
library(ggplot2)
library(viridis)
library(brms)
AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Set up GEE session #####
# use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
# ee <- import("ee")          # Import the Earth Engine library
# ee$Initialize()             # Trigger the authentication
# ##### Get elevation data #####
# # Load raster files
# usDEM <- ee$ImageCollection("JAXA/ALOS/AW3D30/V3_2")$select("DSM")
# usDEM2 <- ee$ImageCollection$mean(usDEM)
# exportRegion <- ee$Geometry$Rectangle(-125, 24, -66, 50)
# usDEM_image <- ee$Image$pixelLonLat()$addBands(usDEM2)
# task <- ee$batch$Export$image$toDrive(image = usDEM_image,
#                                       scale = 50000, # higher resolution available if desired
#                                       region = exportRegion)
# task$start()
# task$status() # Run as many times as you'd like to check on the task. Status should progress from "READY" to "RUNNING" to "COMPLETED"

##### Plot raster #####
elev <- raster::brick("/Users/jacobsocolar/Downloads/US_DEM.tiff")$DSM
elev_f <- raster::focal(elev, w = matrix(data = 1, nrow = 3, ncol = 3), na.rm = TRUE)
elev[is.na(elev)] <- elev_f[is.na(elev)]
raster::plot(elev)

##### Extract raster values #####
new_data_0 <- raster::as.data.frame(elev, xy = T)


##### Get distances #####
setwd("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/")
n_am_map <- readRDS("n_am_map_buffin.Rdata")
warbler_array <- readRDS('warbler_2018_array.RDS')
species <- warbler_array$species

sites_prelim <-  st_as_sf(new_data_0, 
                          coords = c('x', 'y'), crs = 4326)
sites <- st_transform(sites_prelim, AEAstring)

# Load distances
distances_allpoints <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/distance_allpoints.RDS")

distances_2 <- as.data.frame(boot::inv.logit(3*distances_allpoints/400000))
names(distances_2) <- species$SCINAME

new_data_0 <- cbind(new_data_0, distances_2)

new_data <- new_data_0

flattened_data <- readRDS('flattened_data_2018.RDS')
elev_mean <- mean(flattened_data$elev)
elev_sd <- sd(flattened_data$elev)

new_data$elev_scaled <- (new_data$DSM - elev_mean)/elev_sd
new_data$elev_scaled2 <- new_data$elev_scaled^2

##### Get predictions #####
setwd("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM")

bbs_naive <- readRDS("stan_outputs/brms_bbs/naive.RDS")
bbs_dist <- readRDS("stan_outputs/brms_bbs/dist.RDS")
bbs_clip <- readRDS("stan_outputs/brms_bbs/dist_clip_400K.RDS")


get_predictions <- function (sp_code) {
  scientific <- species$SCINAME[species$code == sp_code]
  bbs_naive_predictOcc <- boot::inv.logit(-1 * posterior_linpred(bbs_naive, dpar = "zi", 
                                                                 newdata = 
                                                                   data.frame(trials = 5, species = sp_code, 
                                                                              elev_scaled = new_data$elev_scaled)))
  bbs_dist_predictOcc <- boot::inv.logit(-1 * posterior_linpred(bbs_dist, dpar = "zi", 
                                                                newdata = 
                                                                  data.frame(trials = 5, species = sp_code, 
                                                                             elev_scaled = new_data$elev_scaled,
                                                                             distance_transformed = new_data[[scientific]])))
  bbs_clip_predictOcc <- boot::inv.logit(-1 * posterior_linpred(bbs_clip, dpar = "zi", 
                                                                newdata = 
                                                                  data.frame(trials = 5, species = sp_code, 
                                                                             elev_scaled = new_data$elev_scaled,
                                                                             distance_transformed = new_data[[scientific]])))
  naive_med <- apply(bbs_naive_predictOcc, 2, median)
  dist_med <- apply(bbs_dist_predictOcc, 2, median)
  clip_med <- apply(bbs_clip_predictOcc, 2, median)
  naive_raster <- raster::disaggregate(raster::rasterFromXYZ(data.frame(x=new_data$x, y = new_data$y, z = naive_med)), fact = 10)
  dist_raster <- raster::disaggregate(raster::rasterFromXYZ(data.frame(x=new_data$x, y = new_data$y, z = dist_med)), fact = 10)
  clip_raster <- raster::rasterFromXYZ(data.frame(x=new_data$x, y = new_data$y, z = clip_med))
  clip_raster_clipped <- clip_raster
  clip_raster_clipped[new_data[[scientific]] > boot::inv.logit(3)] <- NA
  clip_raster <- raster::disaggregate(clip_raster, fact = 10)
  clip_raster_clipped <- raster::disaggregate(clip_raster_clipped, fact = 10)
  
  naive <- stars::st_as_stars(naive_raster)
  dist <- stars::st_as_stars(dist_raster)
  clip_raw <- stars::st_as_stars(clip_raster)
  clip_clipped <- stars::st_as_stars(clip_raster_clipped)
  st_crs(naive) <- st_crs(dist) <- st_crs(clip_raw) <- st_crs(clip_clipped) <- 4326
  
  out <- list(naive = naive, dist = dist, 
              clip_raw = clip_raw, clip_clipped = clip_clipped)
}

conus <- spData::us_states
conus <- st_transform(conus, 4326)

maxprob <- 1
out <- get_predictions("PROW")
PROW1 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$naive[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white") + theme(legend.position = "none")
PROW2 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$dist[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white") + theme(legend.position = "none")
PROW3 <- ggplot() + geom_sf(data = conus, fill = "grey80",
                            colour = "grey80", size = 0.2) + 
  geom_stars(data = out$clip_clipped[conus])+ theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "transparent") + theme(legend.position = "none")

out <- get_predictions("RFWA")
RFWA1 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$naive[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white") + theme(legend.position = "none")
RFWA2 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$dist[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white") + theme(legend.position = "none")
RFWA3 <- ggplot() + geom_sf(data = conus, fill = "grey80",
                            colour = "grey80", size = 0.2) + 
  geom_stars(data = out$clip_clipped[conus])+ theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "transparent") + theme(legend.position = "none")

out <- get_predictions("MOWA")
maxprob = 1
MOWA1 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$naive[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white")+ theme(legend.position = "none")
MOWA2 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$dist[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white")+ theme(legend.position = "none")
MOWA3 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$clip_raw[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "transparent") + theme(legend.position = "none")
MOWA4 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$clip_clipped[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "transparent") + theme(legend.position = "none")

dev.off()
pdf(file = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/bbs_maps.pdf",width=10, height=7, family="sans")
ggpubr::ggarrange(OVEN1, BTYW1, MGWA1, 
                  OVEN2, BTYW2, MGWA2,
                  OVEN3, BTYW3, MGWA3, 
                  ncol = 3, nrow = 3)
dev.off()

pdf(file = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/bbs_legend_1.pdf",width=10, height=7, family="sans")
ggplot() + geom_stars(data = out$naive[conus]) + coord_equal() + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white")
dev.off()


pdf(file = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/figures/bbs_PROW_RFWA.pdf",width=7, height=4.95, family="sans")
ggpubr::ggarrange(PROW1,  RFWA1,
                  PROW2, RFWA2,
                  PROW3, RFWA3,
                  ncol = 2, nrow = 3)
dev.off()


pdf(file = "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/figures/bbs_MOWA.pdf",width=7, height=3.3, family="sans")
ggpubr::ggarrange(MOWA1,  
                  MOWA2, 
                  MOWA4, 
                  MOWA3, 
                  ncol = 2, nrow = 2)
dev.off()

maxprob = 1
for(i in 1:nrow(species)){
  print(i)
  out <- get_predictions(species$code[i])
  sp_name <- species$English[i]
  if(grepl("Myrtle", sp_name)){sp_name <- "Myrtle Warbler"}
  if(grepl("Audubon", sp_name)){sp_name <- "Audubon's Warbler"}
  plot1 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$naive[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white")+ theme(legend.position = "none") + ggtitle(paste0(sp_name, ": traditional MSOM")) + theme(plot.title = element_text(size=26))
  plot2 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$dist[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "white")+ theme(legend.position = "none") + ggtitle(paste0(sp_name, ": bMSOM")) + theme(plot.title = element_text(size=26))
  plot3 <- ggplot() + geom_sf(data = conus, fill = "grey80", colour = "grey80", size = 0.2) + geom_stars(data = out$clip_clipped[conus]) + theme_void() + scale_fill_viridis(limits = c(0,maxprob), na.value = "transparent") + theme(legend.position = "none") + ggtitle(paste0(sp_name, ": clipped bMSOM")) + theme(plot.title = element_text(size=26))
  png(file = paste0("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/figures/bbs_all_spp/", species$code[i], ".png"),width=700, height=1000, family="sans")
  print(ggpubr::ggarrange(plot1, plot2, plot3, ncol = 1, nrow = 3))
  dev.off()
}
