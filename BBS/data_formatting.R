##### BBS download and import #####
setwd('/Users/JacobSocolar/Desktop/useful_datasets/BBS/BBS_50_stop')
filenames <- paste0('Fifty', 1:10, '.zip')
for(i in 1:10){
  download.file(
    paste0('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/50-StopData/1997ToPresent_SurveyWide/', filenames[i]),
    destfile = filenames[i])
  unzip(filenames[i])
  file.remove(filenames[i])
}
download.file('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/SpeciesList.txt', 'SpeciesList.txt')
download.file('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/routes.zip', 'routes.zip')
unzip('routes.zip')
file.remove('routes.zip')
`%ni%` <- Negate(`%in%`)
filenames <- paste0('Fifty', 1:10, '.csv')
year <- 2018
bbs_list <- list()
for(i in 1:10){
  bbsImport <- read.csv(filenames[i])
  bbs_list[[i]] <- bbsImport[bbsImport$Year == year, ]
  rm(bbsImport)
}
colnames(bbs_list[[8]])[3] <- 'StateNum'  # typo in column names from BBS, with third column given as 'statenum' instead of 'StateNum'
bbs <- do.call(rbind, bbs_list)
# restrict to continental US (this gives us more data-poor species by excluding Canadian detections)
bbs <- bbs[(bbs$CountryNum == 840) & (bbs$StateNum != 3), ] 
# Aggregate data to five ten-stop blocks per route
bbs$b1 <- as.numeric(rowSums(bbs[,8:17]) > 0)
bbs$b2 <- as.numeric(rowSums(bbs[,18:27]) > 0)
bbs$b3 <- as.numeric(rowSums(bbs[,28:37]) > 0)
bbs$b4 <- as.numeric(rowSums(bbs[,38:47]) > 0)
bbs$b5 <- as.numeric(rowSums(bbs[,48:57]) > 0)
bbs <- bbs[,c(1:7, 58:62)]
# Route data
bbs$routeID <- paste(bbs$CountryNum, bbs$StateNum, bbs$Route, 
                     sep = '_') # create unique route IDs
rte_list <- unique(bbs$routeID) # unique routes run in 2018
routes <- read.csv('routes.csv')
routes$routeID <- paste(routes$CountryNum, routes$StateNum, routes$Route, 
                        sep = '_') # unique route IDs
routes_year <- routes[routes$routeID %in% rte_list, ]
routes_year <- routes_year[routes_year$RouteTypeDetailID == 1, ] # just use the random 50-stop routes
bbs <- bbs[bbs$routeID %in% routes_year$routeID,]
# Species data
species <- read.fwf('SpeciesList.txt', skip = 10, strip.white = TRUE, 
                    header = FALSE,
                    colClasses = c("integer", "integer", rep("character", 7)),
                    widths = c(6, -1, 5, -1, 50, -1, 50, -1, 50, -1, 50, -1, 
                               50, -1, 50, -1, 50),
                    fileEncoding = "iso-8859-1")
species <- species[,c(2,3,6,7,8,9)]
colnames(species) <- c('AOU', 'English', 'Order', 'Family', 'Genus', 'Species')

##### Get species-site-visit array for Parulids #####
warblers <- species[species$Family == "Parulidae", ]
sum(bbs$AOU == 06556) # All YRWA from 2018 are classified as MYWA or AUWA, which is convenient for HBW matching
warblers <- warblers[warblers$AOU %ni% c(6412, 6413, 6556, 6686, 6685), ] # Remove slashes and hybrids
warblers$Species[warblers$Species == 'coronata coronata'] <- 'coronata'
warblers$Species[warblers$Species == 'coronata audoboni'] <- 'auduboni'
# Add two species that regularly breed in lower 48 but never sampled in BBS
extra_warblers <- data.frame(AOU = NA, 
                             English = c("Colima Warbler", 
                                         "Rufous-capped Warbler"), 
                             Order = "Passeriformes",
                             Family = "Parulidae", 
                             Genus = c("Leiothlypis", "Basileuterus"), 
                             Species = c("crissalis", "rufifrons"))
warblers <- rbind(warblers, extra_warblers)
warblers$Genus <- gsub('Oreothlypis', 'Leiothlypis', warblers$Genus) # Update harmonize to range map taxonomy (BirdLife)
warblers$SCINAME <- paste(warblers$Genus, warblers$Species, sep = ' ')
warblers$code <- c("OVEN", "WEWA", "LOWA", "NOWA", "GWWA", "BWWA", "BAWW", 
                   "PROW", "SWWA", "TEWA", "OCWA", "LUWA", "NAWA", "VIWA", 
                   "CONW", "MGWA", "MOWA", "KEWA", "COYE", "HOWA", "AMRE", 
                   "KIWA", "CMWA", "CERW", "NOPA", "TRPA", "MAWA", "BBWA", 
                   "BLBW", "YEWA", "CSWA", "BLPW", "BTBW", "PAWA", "PIWA", 
                   "MYWA", "AUWA", "YTWA", "PRAW", "GRWA", "BTYW", "TOWA", 
                   "HEWA", "GCWA", "BTNW", "CAWA", "WIWA", "RFWA", "PARE", 
                   "COWA", "RCWA")
bbs_warblers <- bbs[bbs$AOU %in% warblers$AOU, ] # COWA and RCWA were never detected in BBS, so we don't need to worry about their codes.
detection_array <- array(data = 0, dim = c(nrow(routes_year), 5, nrow(warblers)))
for(k in 1:nrow(warblers)){
  spdata <- bbs_warblers[bbs_warblers$AOU == warblers$AOU[k], ]
  spdata <- spdata[!duplicated(spdata$routeID),]
  for(i in 1:nrow(routes_year)){
    if(routes_year$routeID[i] %in% spdata$routeID){
      detection_array[i, ,k] <- as.numeric(spdata[spdata$routeID == routes_year$routeID[i], 8:12])
    }
  }
}
warbler_array <- list(detection_array = detection_array, sites = routes_year, species = warblers)
saveRDS(warbler_array, file = paste0('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/warbler_', year, '_array.RDS'))

##### Range maps and distance to range #####
library('sf')
library('ggplot2')
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')
AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
states <- raster::getData(country="USA", level=1)
provinces <- raster::getData(country="Canada", level=1)
estados <- raster::getData(country="Mexico", level=1)
states2 <- states[states$NAME_1 != "Hawaii", ]
states_sf <- st_as_sf(states2)
us_sf <- st_union(states_sf)
provinces_sf <- st_as_sf(provinces)
ca_sf <- st_union(provinces_sf)
estados_sf <- st_as_sf(estados)
mx_sf <- st_union(estados_sf)
n_am_sf <- st_union(us_sf, mx_sf)
n_am_sf <- st_union(n_am_sf, ca_sf)
n_am_AEA <- st_transform(n_am_sf, AEAstring)
n_am_buffin <- st_buffer(n_am_AEA, -30000)
# Read in HydroLAKES, from https://www.hydrosheds.org/products/hydrolakes
lakes <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
big_lakes <- lakes[lakes$Continent == "North America" & lakes$Lake_area > 1000, ]
big_lakes <- st_union(st_make_valid(big_lakes))
big_lakes <- st_transform(big_lakes, AEAstring)
big_lakes_buffout <- st_buffer(big_lakes, 30000)
bc <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/coastal_bc.kml")
la <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/coastal_la.kml")
sl <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/st_lawrence_snippet.kml")
extra_clip <- st_union(bc, la)
extra_clip <- st_union(extra_clip, sl)
extra_clip <- st_transform(extra_clip, AEAstring)
n_am_map <- st_difference(n_am_buffin, big_lakes_buffout)
n_am_map <- st_difference(n_am_map, extra_clip)

##### Birdlife maps #####
# BirdLife range map files from BirdLife via credentials in email from Mark Balman 
# sent 19 February 2020.  See that email for details of proper citation and terms 
# of use.
botw <- "Data/GIS/birdlife_maps/BOTW/BOTW.gdb"
st_layers(botw)
orig_range_maps <- st_read(dsn=botw,layer="All_Species")
recast_range_maps <- st_cast(orig_range_maps, to = "MULTIPOLYGON")

###########
# Load warbler BBS data, update taxonomy, and create sf object for sites
warbler_array <- readRDS(paste0('warbler_', year, '_array.RDS'))
detection_array <- warbler_array$detection_array
species <- warbler_array$species

sites_prelim <-  st_as_sf(warbler_array$sites, 
                          coords = c('Longitude', 'Latitude'), crs = 4326)
sites <- st_transform(sites_prelim, AEAstring)
warbler_ranges <- recast_range_maps[recast_range_maps$SCINAME %in% species$SCINAME, ]
warbler_breeding_prelim <- warbler_ranges[warbler_ranges$SEASONAL %in% c(1,2), ]
warbler_breeding <- st_transform(warbler_breeding_prelim, AEAstring)

# Get distances from each point to the range of each species (based on BirdLife maps)
distances_allpoints <- matrix(data = NA, nrow = nrow(sites), ncol = nrow(species))
for(k in 1:nrow(species)){
  print(k)
  species_range <- st_union(st_make_valid(warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]))
  species_border <- st_cast(species_range, to = "MULTILINESTRING")
  species_border_crop <- st_intersection(species_border, n_am_map)
  distances <- (as.numeric(st_distance(sites, species_range)) > 0) * as.numeric(st_distance(sites, species_border)) -
    (as.numeric(st_distance(sites, species_range)) <= 0) * as.numeric(st_distance(sites, species_border_crop))
  
  # Positive distances for outside-of-range and negative distances for in-range.  Turns 
  # out that as.numeric(st_distance(sites, species_range))>0) is 
  # much faster than st_within(sites, species_range)
  distances_allpoints[, k] <- distances
}

rangemap_distances <- distances_allpoints
saveRDS(rangemap_distances, file = paste0('rangemap_distances_2way_coastclip', year, '.RDS'))

##### Format for occupancy modeling #####
library(sf)
library(ggplot2)
library(reticulate)
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')

year <- 2018

AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


warbler_array <- readRDS(paste0('warbler_', year, '_array.RDS'))
detection_array <- warbler_array$detection_array
species <- warbler_array$species
sites_prelim <-  st_as_sf(warbler_array$sites, coords = c('Longitude', 'Latitude'), crs = 4326)
sites_prelim$lon <- warbler_array$sites$Longitude
sites_prelim$lat <- warbler_array$sites$Latitude

# get elevations for all sites
# get ALOS elevations
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
# Featurecollection of point geometries
geompts <- sapply(1:nrow(sites_prelim),function(x)ee$Geometry$Point(c(sites_prelim$lon[x],sites_prelim$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
# Extract ALOS elevations for all points - combine into dataframe
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_elev$features)),function(x)pts_elev$features[[x]]$properties$AVE_DSM)

sites_prelim$elev_ALOS <- ALOSelev
# Check that this looks right
ggplot(sites_prelim) + geom_sf(aes(col=elev_ALOS))

sites <- st_transform(sites_prelim, AEAstring)

rangemap_distances <- readRDS(paste0('rangemap_distances_2way_coastclip', year, '.RDS'))
distances <- rangemap_distances

detection_slice <- list()
for(k in 1:nrow(species)){
  detection_slice[[k]] <- detection_array[,,k]
}
flattened_data <- as.data.frame(do.call(rbind, detection_slice))
names(flattened_data) <- c('v1', 'v2', 'v3', 'v4', 'v5')
flattened_data$species <- rep(species$code, each = nrow(sites))
flattened_data$sp_id <- rep(1:nrow(species), each = nrow(sites))
flattened_data$site <- rep(sites$routeID, nrow(species))
flattened_data$distance <- as.vector(distances)
# flattened_data$distance_updated <- as.vector(distances_updated)
flattened_data$elev <- rep(sites$elev_ALOS, nrow(species))
flattened_data$Q <- rowSums(flattened_data[,c(1:5)]) > 0
flattened_data$N <- rowSums(flattened_data[,c(1:5)])
flattened_data$distance_scaled <- flattened_data$distance/400000
p <- dd <- vector()
dev.off()
for(sp in unique(flattened_data$species)){
  fd_sp <- flattened_data#[flattened_data$species == sp, ]
  for(i in 1:50){
    j <- -2 + 7*i/50
    j0 <- j - 7/50
    p[i] <- mean(fd_sp$Q[fd_sp$distance_scaled>j0 & fd_sp$distance_scaled < j])
    dd[i] <- j
  }
  dd2 <- boot::inv.logit(dd*3)
  plot(boot::logit(p) ~ dd2, main = sp)
}
flattened_data$distance_transformed <- boot::inv.logit(flattened_data$distance_scaled*3)


saveRDS(flattened_data, file = paste0('flattened_data_', year, '.RDS'))
