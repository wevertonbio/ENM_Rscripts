                          ###################################
                          #### PREPARE ENVIRONMENTAL DATA ###
                          ###################################
                          
#Download the bioclimatic variables from WorldClim
#Here, we are using the resolution of 2.5 arcmin
#We are also using topographic variables (Slope) from ErthEnv: https://www.earthenv.org/topography
#And soil variabls (Clay and SoilType) from soilgrids: https://soilgrids.org/
                          
#Others variables that can be included:

#Biome stability in the last 30k years
# Reference: https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12694
# Link to download: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgeb.12694&file=geb12694-sup-0001-suppinfo1.zip
# Original resolution: 2.5arcmin (~4.5km x 4.5km)

# CHELSA: Climatologies at high resolution for the earth's land surface areas
# Reference: https://www.nature.com/articles/sdata2017122
# Link to download: https://chelsa-climate.org
# Original resolution: 30arc-sec (~1km x 1km)

# Cloud cover:
#   Reference: http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002415
# Link to download: https://www.earthenv.org/cloud
# Original resolution: 30arc-sec (~1km x 1km)

# Topography (elevation, slope, ruggedness, roughness, and others)
# Reference: https://www.nature.com/articles/sdata201840
# Link to download: https://www.earthenv.org/topography
# Original resolution: 30arc-sec (~1km x 1km)

# ENVIRONMENTAL RASTERS FOR ECOLOGICAL MODELING (annual potential evapotranspiration, aridity index, continentality and others)
# Reference: http://onlinelibrary.wiley.com/doi/10.1111/ecog.02880/full
# Link to download: https://envirem.github.io
# Original resolution: 30 arc-seconds, 2.5arcmin, 5arcmin and 10arcmin

# SOILGRIDS - Chemical (pH, Nitrogen, Carbon, and others) and physical (clay, silt, sand, and others) soil variables
# Reference: https://soil.copernicus.org/articles/7/217/2021/
#   Link to download: https://files.isric.org

#We need to mask and resample some variables, making sure all of them have the exact same extent and resolution

#Load packages
library(terra)
library(dplyr)

#Import vector from Neotropic
nt <- vect("https://github.com/wevertonbio/Get_and_Filter_Points/raw/main/Vectors/Neotropic.gpkg")
plot(nt)

####World clim variables####
#Set folder with worldclim variables:
wc_folder <- "C:/Users/wever/Documents/EnvironmentalLayers/WorldClim2.5/"
                          
#Import WORLDCLIM variables
wc <- list.files(path = wc_folder,
                 pattern = ".tif$", full.names = TRUE, recursive = FALSE) %>%
  rast()
wc

#Mask variables to neotropic
wc.neo <- mask(crop(wc, nt), nt)
plot(wc.neo[[1]])

####Topographic variables####

#Set folder with topographic variables:
topo_folder <- "C:/Users/wever/Documents/EnvironmentalLayers/EarthEnvTopographic/"
topo <- list.files(path = topo_folder,
                   pattern = ".tif", full.names = TRUE) %>% rast()
#Cut variables using Neotropic mask
topo.neo <- mask(crop(topo, nt), nt)
plot(topo.neo[[1]])
#See resolution
res(topo.neo)
#It's different from resolution of worldclim variables
res(wc.neo)

#Resampling
#Define factor to aggregate (how many pixels need to be aggregated)
fact <- round(res(wc.neo)/res(topo.neo)[1],0)
topo.r <- terra::aggregate(topo.neo, fact = fact, method="bilinear")
#Resampling to WorldClim, to make sure they have the same extent
topo.r <- resample(topo.r, wc.neo, method = "bilinear")
plot(topo.r[[2]])

#### Soil variables - Continuous ####
#Set folder with soil variables
soil_folder <- "C:/Users/wever/Documents/EnvironmentalLayers/SoilGrids/"

soil <- list.files(path = soil_folder,
                   pattern = ".tif$", full.names = TRUE) %>% rast()
res(soil[[1]])
names(soil)
#Cut variables using Neotropic mask
soil.neo <- mask(crop(soil, nt), nt)
plot(soil.neo[[2]])

#Resample
fact <- round(res(wc.neo)/res(soil.neo)[1],0)
soil.r <- terra::aggregate(soil.neo, fact = fact, method="bilinear")

#Resampling to WorldClim, to make sure they have the same extent
soil.r <- resample(soil.r, wc.neo, method = "bilinear")
plot(soil.r[[1]])

#The data from soilgrids has a lot of missing values. Let's filled this NA values considering the values of the neighboor cells
soil_complete <- focal(soil.r, w=9, fun=mean, na.policy="only", na.rm=T)
plot(soil_complete$Clay.0)


#### Soil variables - Categorical SoilType####
#Variables MostProbable. needed to be downloaded from: https://files.isric.org/soilgrids/latest/data/wrb/
#Open the variable in qgis and resample there (virtual files do not work in R)
soiltype <- rast("C:/Users/wever/Documents/EnvironmentalLayers/SoilGridsResampled/SoilType.tiff")
res(soiltype) #It's not necessary aggregate
plot(soiltype)
names(soiltype)
#Cut variables using Neotropic mask
soiltype.neo <- crop(soiltype, nt)
plot(soiltype.neo)

#Resampling to WorldClim, to make sure they have the same extent
soiltype.r <- resample(soiltype.neo, wc.neo,
                       method = "near") #Use near because of the categorical values
plot(soiltype.r)

#Join all variables in a unique raster
var <- c(wc.neo, topo.r, soil.r, soiltype.r)
#Rename variables
  #Worldclim variables
names(var) <- gsub("wc2.1_2.5m_bio", "Bio", names(var))
names(var) <- gsub("Bio_(\\d)$", "Bio_0\\1", names(var))
names(var) <- gsub("_", "", names(var))
  #Other variables
names(var) %>% dput #Get variables as character
names(var) <- c("Bio01", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", 
                "Bio15", "Bio16", "Bio17", "Bio18", "Bio19", "Bio02", "Bio03", 
                "Bio04", "Bio05", "Bio06", "Bio07", "Bio08", "Bio09", "Altitude", 
                "Slope", "Ruggedness", "TopographicPositionIndex", 
                "Clay", "Nitrogen", "OrgCarbon", "pH", 
                "SoilType")


#Remove Nitrogen, Orgcarbon, pH, Ruggedness, and TopographicPositionIndex
var <- var[[c("Bio01", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", 
                "Bio15", "Bio16", "Bio17", "Bio18", "Bio19", "Bio02", "Bio03", 
                "Bio04", "Bio05", "Bio06", "Bio07", "Bio08", "Bio09", "Altitude",
                "Slope", "Clay", "SoilType")]]
names(var)
#Salvar variaveis
dir.create("Data_example//Current_Neotropic")

writeRaster(var, "Data_example/Current_Neotropic/Variables.tiff", overwrite=TRUE)

