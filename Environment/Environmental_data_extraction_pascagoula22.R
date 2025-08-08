#### Code to download/extract and process environmental data for Pyron et al 2023 test dataset

## Load libraries
library(terra)
library(geodata)
library(ClimateNAr)
library(whitebox)
library(dplyr)


## Set working directory and name of final file
work_dir_main <- "Test data"
final_csv_file <- file.path(work_dir_main, "Pascagoula_environmental.csv")


## Input dataset with Latitude and Longitude columns and ID as rownames
environmental_dataset <- Pascagoula_spatial[, c("Latitude", "Longitude")]


## Define extend of study area
study_area <- terra::ext(min(environmental_dataset$Longitude) - 0.1 * diff(range(environmental_dataset$Longitude)), 
                         max(environmental_dataset$Longitude) + 0.1 * diff(range(environmental_dataset$Longitude)),
                         min(environmental_dataset$Latitude) - 0.1 * diff(range(environmental_dataset$Latitude)),
                         max(environmental_dataset$Latitude) + 0.1 * diff(range(environmental_dataset$Latitude)))


## Extract land mask
countries <- geodata::world(path = work_dir_main)
North_America <- countries[countries$GID_0 %in% c("CAN", "USA")]
North_America <- disagg(North_America)
North_America_areas <- expanse(North_America)
north_america_land_mask <- North_America[order(North_America_areas, decreasing = TRUE)[1:6]]


## Create function to crop and mask raster
crop.and.mask.land <- function(raster_layer, land_mask) {
  cropped <- terra::crop(raster_layer, terra::ext(land_mask))
  masked  <- terra::mask(cropped, land_mask)
  return(masked)
}


## Download and process elevation raster
elevation_raster_file <- file.path(work_dir_main, "elevation_raster.tif")
elevation_raster <- geodata::elevation_global(res = 0.5, path = work_dir_main, download = T) #download and load global elevation raster
if (terra::crs(elevation_raster) != "EPSG:4326") terra::crs(elevation_raster) <- "EPSG:4326" #harmonize CRS
elevation_raster <- crop(elevation_raster, study_area) #crop elevation raster to study area
elevation_raster <- crop.and.mask.land(raster_layer = elevation_raster, land_mask = north_america_land_mask) #mask raster

extract_elevation <- function(df, raster, crs_obj) { #create function to extract raster values
  coords <- terra::vect(df, geom = c("Longitude", "Latitude"), crs = crs_obj) #convert presence data to spatial object (SpatVector) using longitude and latitude
  vals <- terra::extract(raster, coords) #extract values from raster using nearest-cell method
  vals_bilin <- terra::extract(raster, coords, method = "bilinear") #extract elevation values using bilinear interpolation
  final_vals <- vals[, 2]
  na_idx <- which(is.na(final_vals)) #identify which values are NA
  final_vals[na_idx] <- vals_bilin[na_idx, 2] #add to main data frame
  final_vals[is.nan(final_vals)] <- NA #replace NaN values with NA
  return(final_vals)}

environmental_dataset$Elevation <- extract_elevation(environmental_dataset, elevation_raster, "EPSG:4326") #add raster values
names(elevation_raster) <- "elevation" #rename layer
terra::writeRaster(elevation_raster, filename = elevation_raster_file, overwrite = T) #export raster



## Download and process ClimateNAr data (https://climatena.ca/; scale-free)
climateNA_dir <- file.path(work_dir_main, "ClimateNA") #set directory
if (!dir.exists(climateNA_dir)) {dir.create(climateNA_dir, recursive = T)} #create output directory if it does not exist

Lat_Long_Elev <- data.frame(ID1 = rownames(environmental_dataset), #prepare dataframes with selected columns for required format 
                            ID2 = seq_len(nrow(environmental_dataset)),
                            lat = as.numeric(environmental_dataset$Latitude),
                            long = as.numeric(environmental_dataset$Longitude),
                            el = as.numeric(environmental_dataset$Elevation))
write.csv(Lat_Long_Elev, file = file.path(climateNA_dir, "Lat_Long_Elev.csv"), row.names = F) #save dataframes as CSV files
ClimateNAr(inputFile = file.path(climateNA_dir, "Lat_Long_Elev.csv"), #run ClimateNAr
           varList = "YM",
           periodList = "Decade_2011_2020.dcd",
           outDir = climateNA_dir)
file.rename(file.path(work_dir_main, "ClimateNALat_Long_Elev_Decade_2011_2020.csv"), 
            file.path(climateNA_dir, "Lat_Long_Elev_Decade_2011_2020.csv")) #move and rename outfiles
ClimateNAr_data <- read.table(file.path(climateNA_dir, "Lat_Long_Elev_Decade_2011_2020.csv"), sep = ",", header = T) %>% #read results
  dplyr::select(-ID2, -lat, -long, -el) %>% #remove columns
  rename(ID = ID1) #rename column
environmental_dataset$ID <- rownames(environmental_dataset)
environmental_dataset <- merge(environmental_dataset, ClimateNAr_data, by = "ID", all.x = T) #merge into original data
rownames(environmental_dataset) <- environmental_dataset$ID
environmental_dataset$ID <- NULL


## Download and process Enhanced Vegetation Index (EVI) data (https://stac.openlandmap.org/evi_mod13q1.tmwm.inpaint/collection.json; 250m resolution)
output_dir <- file.path(work_dir_main, "EVI_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist
EVI_raster_file <- file.path(output_dir, "EVI_raster.tif") #define EVI raster path
EVI_urls <- c( #specify EVI 2020 URLs
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200101_20200228_go_epsg.4326_v20230608.tif", #Jan-Feb 2020
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200301_20200430_go_epsg.4326_v20230608.tif", #Mar-Apr 2020
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200501_20200630_go_epsg.4326_v20230608.tif", #May-Jun 2020
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200701_20200831_go_epsg.4326_v20230608.tif", #Jul-Aug 2020
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20200901_20201031_go_epsg.4326_v20230608.tif", #Sep-Oct 2020
  "https://s3.openlandmap.org/arco/evi_mod13q1.tmwm.inpaint_p.90_250m_s_20201101_20201231_go_epsg.4326_v20230608.tif") #Nov-Dec 2020

save_paths <- file.path(output_dir, basename(EVI_urls)) #combine directory and filenames
for (i in seq_along(EVI_urls)) { #download files if not already present
  if (!file.exists(save_paths[i])) {download.file(EVI_urls[i], destfile = save_paths[i], mode = "wb") #download raster files
  } else {message(paste("File already exists, skipping:", save_paths[i]))}}

EVI_files <- list.files(output_dir, pattern = "evi_mod13q1.*\\.tif$", full.names = T) #extract file names
EVI_raster <- rast(EVI_files) #stack raster files
names(EVI_raster) <- paste0("EVI_", seq_along(names(EVI_raster))) #rename layers
if (crs(EVI_raster) != "EPSG:4326") crs(EVI_raster) <- "EPSG:4326" #harmonize CRS

EVI_raster <- crop(EVI_raster, study_area) #crop to extent of study area
EVI_raster <- crop.and.mask.land(raster_layer = EVI_raster, land_mask = north_america_land_mask) #mask using land mask

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_EVI_values <- terra::extract(EVI_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_EVI_values <- terra::extract(EVI_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_EVI_layers <- matrix(NA, nrow = nrow(full_EVI_values), ncol = length(names(EVI_raster))) #initialize empty vector for final EVI values
for (i in seq_along(names(EVI_raster))) { #for each of EVI layers, handle NA replacement with bilinear interpolation
  final_EVI_layers[, i] <- full_EVI_values[, i + 1] #get EVI values from nearest-cell
  na_positions <- which(is.na(final_EVI_layers[, i])) #identify positions with NA values
  final_EVI_layers[na_positions, i] <- bilinear_EVI_values[na_positions, i + 1] #replace NA with bilinear interpolated values
  remaining_NAs <- sum(is.na(final_EVI_layers[, i]))  #count remaining NAs
  print(paste("Remaining NAs for EVI", i, ":", remaining_NAs))} #print number of NAs remaining for current EVI variable
for (i in seq_along(names(EVI_raster))) {environmental_dataset[[ paste0("EVI_", i) ]] <- final_EVI_layers[, i]} #combine all variables with original dataframe
EVI_median <- apply(final_EVI_layers, 1, function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = T))) #calculate median EVI for each row across layers
EVI_min <- apply(final_EVI_layers, 1, function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = T))) #calculate min EVI for each row across layers
EVI_max <- apply(final_EVI_layers, 1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = T))) #calculate max EVI for each row across layers
environmental_dataset$EVI_median <- EVI_median #add median EVI values to dataframe
environmental_dataset$EVI_min <- EVI_min #add min EVI values to dataframe
environmental_dataset$EVI_max <- EVI_max #add max EVI values to dataframe 


## Calculate and process terrain metrics (Wilson et al. 2007) (using raster package)
output_dir <- file.path(work_dir_main, "Terrain_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist

TRI_raster <- terrain(elevation_raster, v = "TRI", unit = "degrees", neighbors = 8) #Terrain ruggedness index (TRI) – mean of absolute differences
TPI_raster <- terrain(elevation_raster, v = "TPI", unit = "degrees", neighbors = 8) #Topographic position index (TPI) – difference from mean of 8 neighbors
roughness_raster <- terrain(elevation_raster, v = "roughness", unit = "degrees", neighbors = 8) ##Terrain roughness index – difference between max and min of 8 neighbors
slope_raster  <- terrain(elevation_raster, v = "slope", unit = "radians") #extract slope from elevation
aspect_raster <- terrain(elevation_raster, v = "aspect", unit = "radians") #extract aspect from elevation
HLI_raster <- (0.339 + 0.806 * cos(slope_raster) - 0.196 * sin(slope_raster) * cos(aspect_raster) - 0.482 * sin(slope_raster) * sin(aspect_raster)) #calculate heat load index (following McCune & Keon 2002 formula)

TRI_raster[is.na(TRI_raster)] <- focal(TRI_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(TRI_raster)] #fill NA values in rasters using focal median of neighboring cells
TPI_raster[is.na(TPI_raster)] <- focal(TPI_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(TPI_raster)] #fill NA values in rasters using focal median of neighboring cells
roughness_raster[is.na(roughness_raster)] <- focal(roughness_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(roughness_raster)] #fill NA values in rasters using focal median of neighboring cells
slope_raster[is.na(slope_raster)] <- focal(slope_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(slope_raster)] #fill NA values in rasters using focal median of neighboring cells
aspect_raster[is.na(aspect_raster)] <- focal(aspect_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(aspect_raster)] #fill NA values in rasters using focal median of neighboring cells
HLI_raster[is.na(HLI_raster)] <- focal(HLI_raster, w = matrix(1, 3, 3), fun = median, na.rm = T)[is.na(HLI_raster)] #fill NA values in rasters using focal median of neighboring cells

zfactor <- 1 / (111320 * cos(mean(ext(elevation_raster)[c(3, 4)]) * pi / 180)) #compute zfactor for slope (meters to degrees conversion)
wbt_init() #initialize WhiteboxTools
slope_file <- file.path(output_dir, "slope_raster_degrees.tif")
wbt_slope(dem = elevation_raster_file, output = slope_file, #compute slope in degrees using zfactor
          zfactor = zfactor, units = "degrees")
SCA_file <- file.path(output_dir, "SCA_raster.tif")
wbt_d8_flow_accumulation(input = elevation_raster_file, output = SCA_file, #compute Specific Contributing Area (SCA) — D8 flow accumulation
                         out_type = "specific contributing area") #important for TWI
TWI_file <- file.path(output_dir, "TWI_raster.tif")
wbt_wetness_index(sca = SCA_file, slope = slope_file, output = TWI_file) #compute Topographic Wetness Index (TWI)
TWI_raster <- rast(TWI_file) #load TWI raster
TWI_raster[is.na(TWI_raster)] <- focal(TWI_raster, w = matrix(1, 3, 3), fun = median, na.rm = TRUE)[is.na(TWI_raster)] #fill NA values using focal median

Terrain_raster <- c(TRI_raster, TPI_raster, roughness_raster, HLI_raster, TWI_raster) #combine rasters to stack
names(Terrain_raster) <- c("TRI", "TPI", "roughness", "HLI", "TWI")
Terrain_raster <- lapply(Terrain_raster, function(x) {if (crs(x) != "EPSG:4326") crs(x) <- "EPSG:4326"; x}) #harmonize CRS across all rasters
Terrain_raster <- terra::rast(Terrain_raster) #rebuild SpatRaster from list
names(Terrain_raster) <- c("TRI", "TPI", "roughness", "HLI", "TWI")

Terrain_raster <- crop(Terrain_raster, study_area) #crop to extent of study area
Terrain_raster <- crop.and.mask.land(raster_layer = Terrain_raster, land_mask = north_america_land_mask)

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_Terrain_values <- terra::extract(Terrain_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_Terrain_values <- terra::extract(Terrain_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_Terrain_layers <- matrix(NA, nrow = nrow(full_Terrain_values), ncol = length(names(Terrain_raster))) #initialize empty matrix for final Terrain values
for (i in seq_len(nlyr(Terrain_raster))) { #for each of Terrain layers, handle NA replacement with bilinear interpolation
  final_Terrain_layers[, i] <- full_Terrain_values[, i + 1] #get Terrain values from nearest-cell
  na_positions <- which(is.na(final_Terrain_layers[, i])) #identify positions with NA values
  final_Terrain_layers[na_positions, i] <- bilinear_Terrain_values[na_positions, i + 1] #replace NA with bilinear interpolated values
  remaining_NAs <- sum(is.na(final_Terrain_layers[, i]))  #count remaining NAs
  print(paste("Remaining NAs for Terrain", i, ":", remaining_NAs))} #print number of NAs remaining for current Terrain variable
for (i in seq_len(nlyr(Terrain_raster))) {environmental_dataset[[names(Terrain_raster)[i]]] <- final_Terrain_layers[, i]} #combine all variables with original dataframe


## Download and process ENVIREM variables (Title & Bemmels 2018, Ecography: https://envirem.github.io/; resolution: 30 arcsec)
output_dir <- file.path(work_dir_main, "ENVIREM_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist
ENVIREM_raster_file <- file.path(output_dir, "ENVIREM_raster.tif") #define raster path

envirem_urls <- c("https://deepblue.lib.umich.edu/data/downloads/jh343s427",
                  "https://deepblue.lib.umich.edu/data/downloads/2f75r815f",
                  "https://deepblue.lib.umich.edu/data/downloads/x633f1111",
                  "https://deepblue.lib.umich.edu/data/downloads/9c67wm919")
zip_paths <- file.path(output_dir, paste0("envirem_", seq_along(envirem_urls), ".zip")) #specify local zip file paths
for (i in seq_along(envirem_urls)) {
  if (!file.exists(zip_paths[i])) {
    download.file(envirem_urls[i], destfile = zip_paths[i], mode = "wb") #download and unzip rasters
    unzip(zip_paths[i], exdir = output_dir)}}

envirem_tifs <- list.files(path = file.path(work_dir_main, "ENVIREM_data"), 
                           pattern = "^current_.*\\.tif$", 
                           full.names = TRUE)
layer_names <- sub("\\.tif$", "", basename(envirem_tifs))
layer_names <- sub("^current_30arcsec_", "", layer_names)
envirem_layers <- setNames(lapply(envirem_tifs, rast), layer_names)
envirem_layers <- lapply(envirem_layers, function(x) {if (crs(x) != "EPSG:4326") crs(x) <- "EPSG:4326"; x}) #harmonize CRS across all rasters
ENVIREM_raster <- rast(envirem_layers) #combine into one SpatRaster stack
ENVIREM_raster <- crop(ENVIREM_raster, study_area) #crop to extent of study area
ENVIREM_raster <- crop.and.mask.land(raster_layer = ENVIREM_raster, land_mask = north_america_land_mask) #mask raster

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_ENVIREM_values <- terra::extract(ENVIREM_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_ENVIREM_values <- terra::extract(ENVIREM_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_ENVIREM_layers <- matrix(NA, nrow = nrow(full_ENVIREM_values), ncol = length(names(ENVIREM_raster))) #initialize empty vector for final ENVIREM values
for (i in 1:length(names(ENVIREM_raster))) { #for each of ENVIREM layers, handle NA replacement with bilinear interpolation
  final_ENVIREM_layers[, i] <- full_ENVIREM_values[, i + 1] #get ENVIREM values from nearest-cell
  na_positions <- which(is.na(final_ENVIREM_layers[, i])) #identify positions with NA values
  final_ENVIREM_layers[na_positions, i] <- bilinear_ENVIREM_values[na_positions, i + 1] #replace NA with bilinear interpolated values
  remaining_NAs <- sum(is.na(final_ENVIREM_layers[, i]))  #count remaining NAs
  print(paste("Remaining NAs for ENVIREM", i, ":", remaining_NAs))} #print number of NAs remaining for current ENVIREM variable
for (i in seq_len(nlyr(ENVIREM_raster))) {environmental_dataset[[ names(ENVIREM_raster)[i] ]] <- final_ENVIREM_layers[, i]} #combine all variables with original dataframe


#### Download and process human footprint data for 2009 (using geodata package; 30-seconds resolution)
output_dir <- file.path(work_dir_main, "footprint_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist
footprint_raster_file <- file.path(output_dir, "footprint_raster.tif") #define raster path


footprint_raster <- geodata::footprint(year = 2009, path = output_dir) #download raster
if (crs(footprint_raster) != "EPSG:4326") crs(footprint_raster) <- "EPSG:4326"
names(footprint_raster) <- "Human_footprint"
footprint_raster <- crop(footprint_raster, study_area) #crop to extent of study area
footprint_raster <- crop.and.mask.land(raster_layer = footprint_raster, land_mask = north_america_land_mask) #mask raster

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_footprint_values <- terra::extract(footprint_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_footprint_values <- terra::extract(footprint_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_footprint_layers <- matrix(NA, nrow = nrow(full_footprint_values), ncol = 1) #initialize final values matrix
final_footprint_layers[, 1] <- full_footprint_values[, 2] #get values from nearest-cell
na_positions <- which(is.na(final_footprint_layers[, 1])) #identify NA positions
final_footprint_layers[na_positions, 1] <- bilinear_footprint_values[na_positions, 2] #replace NAs with bilinear interpolated values
print(paste("Remaining NAs for footprint:", sum(is.na(final_footprint_layers[, 1])))) #print number of NAs remaining
environmental_dataset$footprint_value <- final_footprint_layers[, 1] #combine all variables with original dataframe


## Download and process landcover variables (using geodata package)
output_dir <- file.path(work_dir_main, "landcover_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist
landcover_raster_file <- file.path(output_dir, "landcover_raster.tif") #define raster path

landcover_vars <- c("trees", "grassland", "shrubs", "cropland", "built", "bare", #specify landcover variables
                    "snow", "water", "wetland", "mangroves", "moss")
landcover_layers <- list()
for (var in landcover_vars) {
  landcover_path <- file.path(output_dir, paste0(var, ".tif"))
  if (!file.exists(landcover_path)) {
    message(paste("Downloading landcover layer:", var))
    lc <- geodata::landcover(var = var, path = output_dir, download = T) #download and save .tif
  } else {
    message(paste("Loading existing layer:", var))
    lc <- rast(landcover_path)}
  landcover_layers[[var]] <- lc}

landcover_layers <- lapply(landcover_layers, function(x) {if (crs(x) != "EPSG:4326") crs(x) <- "EPSG:4326"; x}) #harmonize CRS across all rasters
landcover_raster <- rast(landcover_layers) #combine into one SpatRaster stack
names(landcover_raster) <- landcover_vars
landcover_raster <- terra::crop(landcover_raster, study_area) #crop raster
landcover_raster <- terra::mask(landcover_raster, north_america_land_mask) #mask raster

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_landcover_values <- terra::extract(landcover_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_landcover_values <- terra::extract(landcover_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_landcover_layers <- matrix(NA, nrow = nrow(full_landcover_values), ncol = length(names(landcover_raster))) #initialize empty vector for final landcover values
for (i in 1:length(names(landcover_raster))) { #for each of landcover layers, handle NA replacement with bilinear interpolation
  final_landcover_layers[, i] <- full_landcover_values[, i + 1] #get landcover values from nearest-cell
  na_positions <- which(is.na(final_landcover_layers[, i])) #identify positions with NA values
  final_landcover_layers[na_positions, i] <- bilinear_landcover_values[na_positions, i + 1] #replace NA with bilinear interpolated values
  remaining_NAs <- sum(is.na(final_landcover_layers[, i]))  #count remaining NAs
  print(paste("Remaining NAs for landcover", i, ":", remaining_NAs))} #print number of NAs remaining for current landcover variable
for (i in seq_len(nlyr(landcover_raster))) {environmental_dataset[[names(landcover_raster)[i]]] <- final_landcover_layers[, i]} #combine all variables with original dataframe


## Download and process soil variables (using geodata package: derived from the SoilGRIDS database; 30second resolution)
output_dir <- file.path(work_dir_main, "soil_data") #specify directory
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = T)} #create output directory if it does not exist
soil_raster_file <- file.path(output_dir, "soil_raster.tif") #define raster path

Soil_variables <- c("bdod", "cfvo", "clay", "nitrogen", "ocd", "phh2o", "sand", "silt", "soc") #specify soil variables to download
depth <- 5 #specify soil depths for downloading
depth_range <- "0-5" #specify depth ranges for filenames

soil_data_list <- list() #initialize empty list to store individual rasters

for (var in Soil_variables) {
  soil_file <- file.path(output_dir, "soil_world", paste0(var, "_", depth_range, "cm_mean_30s.tif"))
  if (file.exists(soil_file)) {
    message("File for ", var, " at ", depth_range, "cm already exists — skipping download")
    soil_data <- terra::rast(soil_file)
  } else {
    message("Downloading file for ", var, " at ", depth_range, "cm")
    soil_data <- geodata::soil_world(var = var, depth = depth, stat = "mean", path = output_dir, download = T)} #download rasters
  soil_data_list[[paste(var, depth_range, "cm", sep = "_")]] <- soil_data}

soil_data_list <- lapply(soil_data_list, function(x) {if (crs(x) != "EPSG:4326") crs(x) <- "EPSG:4326"; x})
soil_raster <- rast(soil_data_list) #combine into one SpatRaster stack
names(soil_raster) <- c("Bulk_density", "Coarse_fragments_volume", "Clay_fraction", 
                        "Nitrogen_content", "Organic_carb_density", "pH_H2O", 
                        "Sand_fraction", "Silt_fraction", "Soil_organic_carb")

na_only_layers <- sapply(1:nlyr(soil_raster), function(i) {all(is.na(terra::minmax(soil_raster[[i]], compute = T)[, 1]))}) #check if min is NA
if (any(na_only_layers)) {removed_layers <- names(soil_raster)[na_only_layers]
message("Removing raster layers with only NA values: ", paste(removed_layers, collapse = ", "))}
soil_raster <- soil_raster[[which(!na_only_layers)]] #keep only layers that are not all NA

soil_raster <- crop(soil_raster, study_area) #crop to extent of study area
soil_raster <- crop.and.mask.land(raster_layer = soil_raster, land_mask = north_america_land_mask) #mask raster

coordinate_vector <- terra::vect(environmental_dataset[, c("Longitude", "Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326") #convert presence data to spatial object (SpatVector) using longitude and latitude
full_soil_values <- terra::extract(soil_raster, coordinate_vector) #extract values from raster using nearest-cell method
bilinear_soil_values <- terra::extract(soil_raster, coordinate_vector, method = "bilinear") #extract values using bilinear interpolation
final_soil_layers <- matrix(NA, nrow = nrow(full_soil_values), ncol = length(names(soil_raster))) #initialize empty vector for final soil values
for (i in seq_along(names(soil_raster))) { #for each of soil layers, handle NA replacement with bilinear interpolation
  final_soil_layers[, i] <- full_soil_values[, i + 1] #get soil values from nearest-cell
  na_positions <- which(is.na(final_soil_layers[, i])) #identify positions with NA values
  final_soil_layers[na_positions, i] <- bilinear_soil_values[na_positions, i + 1] #replace NA with bilinear interpolated values
  remaining_NAs <- sum(is.na(final_soil_layers[, i])) #count remaining NAs
  print(paste("Remaining NAs for soil", i, ":", remaining_NAs))} #print number of NAs remaining for current soil variable
for (i in seq_len(nlyr(soil_raster))) {environmental_dataset[[names(soil_raster)[i]]] <- final_soil_layers[, i]} #combine all variables with original dataframe


## Delete intermediate files and raster files
file.remove(elevation_raster_file)
unlink(file.path(work_dir_main, "elevation"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "ENVIREM_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "footprint_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "landcover_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "soil_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "Terrain_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "EVI_data"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "ClimateNA"), recursive = TRUE, force = TRUE)
unlink(file.path(work_dir_main, "gadm"), recursive = TRUE, force = TRUE)


## Save final dataset as csv file
write.csv(environmental_dataset, file = final_csv_file, row.names = TRUE)

