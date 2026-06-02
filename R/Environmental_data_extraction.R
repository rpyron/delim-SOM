#### Install NicheDiv ##########################################################
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv", force = T, upgrade = "never")




#### Pyron 2023 ################################################################

## Clear environment
rm(list = ls())
invisible(gc())


## Set base directory
base_dir <- "C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Test data/Pyron_2023"
setwd(base_dir)


## Set results directories
results_dir <- base_dir
intermediate_files_dir_name <- "Intermediate_files"


## Set main input parameters
occurrence_data_file <- "monticola71.csv"
csv_occurrence_out_file <- "Monticola71_environmental.csv"

Latitude_col <- "Lat"
Longitude_col <- "Long"
ID_column <- "Sample"
buffer_km <- 5



## Import occurrences
occurrence_data <- read.csv(occurrence_data_file)
occurrence_data <- occurrence_data[, c(ID_column, Latitude_col, Longitude_col)]



## Extract environmental data and background
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data, #input data.frame with coords (rownames need to be unique IDs)
                                     longitude.col = Longitude_col, #column name for longitude
                                     latitude.col = Latitude_col,  #column name for latitude
                                     generate.background.data = FALSE, #whether to sample background points
                                     buffer.km = buffer_km, #buffer in km around convex hull for accessible area
                                     env.datasets = c("elevation", "ClimateNA", "EVI", "terrain", 
                                                      "ENVIREM", "footprint", "landcover", "soil", "forest_height", 
                                                      "atmosphere", "nightlight", "burned_area", "snow_water_equivalent", 
                                                      "daylength", "soil_moisture"),
                                     csv.occurrence.out.file = csv_occurrence_out_file, #output CSV name for occurrence+env values
                                     output.dir = results_dir, #main output directory
                                     intermediate.files.dir = intermediate_files_dir_name, #name for subfolder for cached intermediate CSVs
                                     overwrite = TRUE) #overwrite existing occurrence/background CSVs


## Evaluate extracted data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
head(Env_data_occurrences)





#### Pyron et al 2022 ##########################################################

## Clear environment
rm(list = ls())
invisible(gc())


## Set base directory
base_dir <- "C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Test data/Pyron_et_al_2022"
setwd(base_dir)


## Set results directories
results_dir <- base_dir
intermediate_files_dir_name <- "Intermediate_files"


## Set main input parameters
occurrence_data_file <- "pascagoula22.csv"
csv_occurrence_out_file <- "Pascagoula22_environmental.csv"

Latitude_col <- "Lat"
Longitude_col <- "Long"
ID_column <- "Sample"
buffer_km <- 5



## Import occurrences
occurrence_data <- read.csv(occurrence_data_file)
occurrence_data <- occurrence_data[, c(ID_column, Latitude_col, Longitude_col)]



## Extract environmental data and background
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data, #input data.frame with coords (rownames need to be unique IDs)
                                     longitude.col = Longitude_col, #column name for longitude
                                     latitude.col = Latitude_col,  #column name for latitude
                                     generate.background.data = FALSE, #whether to sample background points
                                     buffer.km = buffer_km, #buffer in km around convex hull for accessible area
                                     env.datasets = c("elevation", "ClimateNA", "EVI", "terrain", 
                                                      "ENVIREM", "footprint", "landcover", "soil", "forest_height", 
                                                      "atmosphere", "nightlight", "burned_area", "snow_water_equivalent", 
                                                      "daylength", "soil_moisture"),
                                     csv.occurrence.out.file = csv_occurrence_out_file, #output CSV name for occurrence+env values
                                     output.dir = results_dir, #main output directory
                                     intermediate.files.dir = intermediate_files_dir_name, #name for subfolder for cached intermediate CSVs
                                     overwrite = TRUE) #overwrite existing occurrence/background CSVs


## Evaluate extracted data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
head(Env_data_occurrences)




#### Pyron et al 2024 ##########################################################

## Clear environment
rm(list = ls())
invisible(gc())


## Set base directory
base_dir <- "C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Test data/Pyron_et_al_2024"
setwd(base_dir)


## Set results directories
results_dir <- base_dir
intermediate_files_dir_name <- "Intermediate_files"


## Set main input parameters
occurrence_data_file <- "aeneus56.csv"
csv_occurrence_out_file <- "Aeneus56_environmental.csv"

Latitude_col <- "Lat"
Longitude_col <- "Long"
buffer_km <- 5
ID_column <- "Sample"


## Import occurrences
occurrence_data <- read.csv(occurrence_data_file)
occurrence_data <- occurrence_data[, c(ID_column, Latitude_col, Longitude_col)]



## Extract environmental data and background
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data, #input data.frame with coords (rownames need to be unique IDs)
                                     longitude.col = Longitude_col, #column name for longitude
                                     latitude.col = Latitude_col,  #column name for latitude
                                     generate.background.data = FALSE, #whether to sample background points
                                     buffer.km = buffer_km, #buffer in km around convex hull for accessible area
                                     env.datasets = c("elevation", "ClimateNA", "EVI", "terrain", 
                                                      "ENVIREM", "footprint", "landcover", "soil", "forest_height", 
                                                      "atmosphere", "nightlight", "burned_area", "snow_water_equivalent", 
                                                      "daylength", "soil_moisture"),
                                     csv.occurrence.out.file = csv_occurrence_out_file, #output CSV name for occurrence+env values
                                     output.dir = results_dir, #main output directory
                                     intermediate.files.dir = intermediate_files_dir_name, #name for subfolder for cached intermediate CSVs
                                     overwrite = TRUE) #overwrite existing occurrence/background CSVs


## Evaluate extracted data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
head(Env_data_occurrences)




#### Dupuis et al 2018 ##########################################################

## Clear environment
rm(list = ls())
invisible(gc())


## Set base directory
base_dir <- "C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Test data/Dupuis_et_al_2018"
setwd(base_dir)


## Set results directories
results_dir <- base_dir
intermediate_files_dir_name <- "Intermediate_files"


## Set main input parameters
occurrence_data_file <- "Polygonia_metadata.csv"
csv_occurrence_out_file <- "Polygonia_environmental.csv"

Latitude_col <- "Latitude"
Longitude_col <- "Longitude"
ID_column <- "ID"
buffer_km <- 5



## Import occurrences
occurrence_data <- read.csv(occurrence_data_file, sep = ";")
occurrence_data <- occurrence_data[, c(ID_column, Latitude_col, Longitude_col)]



## Extract environmental data and background
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data, #input data.frame with coords (rownames need to be unique IDs)
                                     longitude.col = Longitude_col, #column name for longitude
                                     latitude.col = Latitude_col,  #column name for latitude
                                     generate.background.data = FALSE, #whether to sample background points
                                     buffer.km = buffer_km, #buffer in km around convex hull for accessible area
                                     env.datasets = c("elevation", "ClimateNA", "EVI", "terrain", 
                                                      "ENVIREM", "footprint", "landcover", "soil", "forest_height", 
                                                      "atmosphere", "nightlight", "burned_area", "snow_water_equivalent", 
                                                      "daylength", "soil_moisture"),
                                     csv.occurrence.out.file = csv_occurrence_out_file, #output CSV name for occurrence+env values
                                     output.dir = results_dir, #main output directory
                                     intermediate.files.dir = intermediate_files_dir_name, #name for subfolder for cached intermediate CSVs
                                     overwrite = TRUE) #overwrite existing occurrence/background CSVs


## Evaluate extracted data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
head(Env_data_occurrences)




#### van Elst et al 2024 ##########################################################
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Daniel-1232/NicheDiv", force = T, upgrade = "never")

## Clear environment
rm(list = ls())
invisible(gc())


## Set base directory
base_dir <- "C:/Users/danie/Desktop/PhD research/Manuscripts/SOM package/Test data/van_Elst_et_al_2024"
setwd(base_dir)


## Set results directories
results_dir <- base_dir
intermediate_files_dir_name <- "Intermediate_files"


## Set main input parameters
occurrence_data_file <- "Microcebus_multiple_data_combined.csv"
csv_occurrence_out_file <- "Microcebus_environmental.csv"

Latitude_col <- "latitude"
Longitude_col <- "longitude"
ID_column <- "Individual.ID"
buffer_km <- 5



## Import occurrences
occurrence_data <- read.csv(occurrence_data_file)
occurrence_data <- occurrence_data[, c(ID_column, Latitude_col, Longitude_col)]



## Extract environmental data and background
NicheDiv::extract.env.and.background(occurrence.data = occurrence_data, #input data.frame with coords (rownames need to be unique IDs)
                                     longitude.col = Longitude_col, #column name for longitude
                                     latitude.col = Latitude_col,  #column name for latitude
                                     generate.background.data = FALSE, #whether to sample background points
                                     buffer.km = buffer_km, #buffer in km around convex hull for accessible area
                                     env.datasets = c("elevation", "EVI", "terrain", 
                                                      "ENVIREM", "footprint", "landcover", "soil", "forest_height", 
                                                      "atmosphere", "nightlight", "burned_area", "soil_moisture"),
                                     csv.occurrence.out.file = csv_occurrence_out_file, #output CSV name for occurrence+env values
                                     output.dir = results_dir, #main output directory
                                     intermediate.files.dir = intermediate_files_dir_name, #name for subfolder for cached intermediate CSVs
                                     overwrite = TRUE) #overwrite existing occurrence/background CSVs


## Evaluate extracted data
Env_data_occurrences <- read.csv(file.path(results_dir, csv_occurrence_out_file))
dim(Env_data_occurrences)
head(Env_data_occurrences)

