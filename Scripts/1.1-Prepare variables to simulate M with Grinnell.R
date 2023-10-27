#############################################
#### PREPARE LGM VARIABLES TO SIMULATE M ####
#############################################

#To create the M with the Grinnell package, we need the variables from Last Glacial Maximum
#Download the variables from here: http://www.worldclim.com/past
#We are using the CCSM GCM.

#Current variables at 2.5arcmin already prepared with script 1

#Load packages
library(terra)
library(dplyr)

#Import vector from Neotropic
nt <- vect("https://github.com/wevertonbio/Get_and_Filter_Points/raw/main/Vectors/Neotropic.gpkg")
plot(nt)

#Import LGM variables (2.5 arcmin)
lgm2.5 <- list.files("C:/Users/wever/Documents/EnvironmentalLayers/wc_2_5m_CCSM_21k_bio/2_5m/",
                     pattern = ".bil$", full.names = T) %>%
  rast()
#Rename variables
names(lgm2.5) <- gsub("wc_2_5m_CCSM_21k_bio", "Bio", names(lgm2.5))
names(lgm2.5) <- gsub("Bio_(\\d)$", "Bio0\\1", names(lgm2.5))
names(lgm2.5) <- gsub("_", "", names(lgm2.5))

#Cut variables to neotropic
lgm2.5_neo <- crop(lgm2.5, nt, mask = TRUE)


#Temperature variables are in ºC x 10. Fix it:
plot(lgm2.5_neo$Bio01)
lgm10.2.5 <- lgm2.5_neo[[c("Bio01", "Bio02", "Bio04",
                       "Bio05", "Bio06","Bio07",
                       "Bio08", "Bio09", "Bio10",
                       "Bio11")]] / 10
lgm2.5final <- c(lgm10.2.5, lgm2.5_neo[[c("Bio03", "Bio12", "Bio13", "Bio14",
                                      "Bio15", "Bio16","Bio17", "Bio18",
                                      "Bio19")]])
plot(lgm2.5final$Bio01)

#Write rasters
dir.create("Data_example/LGMcc_Neotropic/")
writeRaster(lgm2.5final, "Data_example/LGMcc_Neotropic/LGMcc_2.5arcmin.tif",
            overwrite=TRUE)