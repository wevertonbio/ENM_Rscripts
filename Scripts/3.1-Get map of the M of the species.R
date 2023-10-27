########################################
#### Get map of the M of the specie ####
########################################

#Load packages
library(terra)
library(dplyr)
library(pbapply)
library(dplyr)
library(ggplot2)
library(marmap)
library(raster)
library(geobr)
library(metR)
library(future.apply)
library(progressr)
library(data.table)
library(ggnewscale)
library(tidyterra)


#Create directory to save images of M simulations
dir.create("Data_example/M_simulations/Maps", recursive = T)

#Import all occurrences
occ <- read.csv("Data_example/Occurrences.csv")
unique(occ$species) %>% length()
#Get species
all_spp <- unique(occ$species)

#Get all Ms
path_to_m <- "Data_example/M_simulations/Shapefiles/"
all_m <- list.files(path = path_to_m,
                    pattern = ".gpkg$", full.names = TRUE, recursive = T)
all_m <- all_m[grepl("M_OK", all_m)]
#Get species
spp_with_m <- gsub(".gpkg", "", basename(all_m))
#Species with m
spp <- intersect(all_spp, spp_with_m)
#Species withou m - Check what happened if != character ()0!
spp_without_m <- setdiff(all_spp, spp_with_m)

#Get map of altitude
alt <- rast("Data_example/Current_Neotropic/Variables.tiff")
alt <- alt$Altitude

#Import Brasil map
br <- geobr::read_state(code_state = "all")

#Import south america map
sa <- vect("Vectors/South_America.gpkg")

#Start looping
pblapply(seq_along(spp), function(i){
  #Get specie
  sp <- spp[1]
  
  #Get M of the specie
  m <- paste0("Data_example/M_simulations/Shapefiles/",
              sp, "/M_OK/", sp, ".gpkg") %>% vect()
  
  
  
  #Cut altraster using limits
  crs(m) <- crs(alt)
  alt_m <- crop(alt, buffer(m, width = 1000*1000), mask = TRUE)
  
  # Convert alt to data frame
  alt_df <- as.data.frame(alt_m, xy = TRUE, na.rm = T)
  
  #Get limits to plot
  bb <- ext(m)
  
  # Create base map
  map <- ggplot() +
    geom_tile(data = alt_df, aes(x = x, y = y, fill = Altitude)) +
    scale_fill_etopo() +
    geom_sf(data = sa, linewidth = 0.2, colour = "grey40", alpha = 0) +
    geom_sf(data = br, linewidth = 0.2, colour = "grey40", alpha = 0) +
    geom_sf(data = m, colour = "darkred", alpha = 0, linewidth = 1.2) +
    new_scale_fill() +
    geom_point(data = occ, aes(x = x, y = y),
               size = 1.4,shape = 21, colour = "black", fill = "yellow1") +
    coord_sf(xlim = c(bb[1] - 1, bb[2] + 1),
             ylim = c(bb[3] - 1, bb[4] + 1),
             expand = F) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle(sp) +
    metR::scale_x_longitude(ticks = 2) + metR::scale_y_latitude(ticks = 2) +
    theme_classic() +
    theme(legend.position = "none",
          panel.border = element_rect(colour = "black", fill=NA, size=2),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  #map
  ggsave(paste0("Data_example/M_simulations/Maps/", sp, ".png"), map, dpi = 300,
         width = 6, height = 7)
  
  
})



#plot(dem.m)

#Get m in sf to plot
m_sf <- st_as_sf(m_final)

# Create base map
map <- ggplot() +
  geom_tile(data = dem_df, aes(x = x, y = y, fill = elevation_1KMmn_SRTM)) +
  scale_fill_etopo() +
  geom_sf(data = sa, linewidth = 0.2, colour = "grey40", alpha = 0) +
  geom_sf(data = br, linewidth = 0.2, colour = "grey40", alpha = 0) +
  geom_sf(data = m_sf, colour = "darkred", alpha = 0, linewidth = 1.2) +
  new_scale_fill() +
  geom_point(data = occ_m, aes(x = longitude, y = latitude,
                               fill = Colored), size = 1.4,
             shape = 21, colour = "black") +
  scale_fill_manual(values = occ_m$Colored, guide = "none") +
  coord_sf(xlim = c(bb[1], bb[2]),
           ylim = c(bb[3], bb[4]),
           expand = F) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle(sp_name) +
  metR::scale_x_longitude(ticks = 2) + metR::scale_y_latitude(ticks = 2) +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        panel.background = element_rect(fill = 'aliceblue', colour = NA),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
#map
ggsave(paste0("M/PNG_M/", save_dir, "/", sp_name, ".png"), map, dpi = 300,
       width = 6, height = 7)
gc()}, #Print progress
error = function(e) NULL)
})


