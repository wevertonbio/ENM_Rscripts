#### 2-Filter by distance and spatial autocorrelation ####
#With this function, we are going to test different distances to filter/thinning the occurrences.
#We are going to select the distance that generate a final data set with low average spatial autocorrelation with the highest number occurrences.

#Load packages
library(dplyr)
library(flexsdm)
library(terra)
library(pbapply)
library(data.table)
library(parallel)
library(moranfast)


#Import function to crop variable and do PCA from: https://github.com/wevertonbio/pre_kuenm2
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/prepare_var.R")

#Import occurrences: here, we are going to work with the occurrences of Araucaria angustifolia from:
df <- read.csv("https://raw.githubusercontent.com/wevertonbio/Get_and_Filter_Points/main/Examples/Check_Points/8-Rescued_from_Sea.csv")
#Get and rename columns
df <- df %>% dplyr::select(species,
                           x = decimalLongitude.new,
                           y = decimalLatitude.new,
                           year = year.new)

#Remove duplicates in space and time
df <- df %>% distinct(species, x, y, year)

#Create ID
df$ID <- row.names(df)

#Import variables (only worldclim)
var <- rast("Data_example/Current_Neotropic/Variables.tiff")
var <- var[[grepl("Bio", names(var))]]
names(var)
plot(var$Bio01)

#Wrap variables to run in parallel
var_wrap <- wrap(var)

#### Geographical thinning and calculate spatial autocorrelation####
#See more in: https://github.com/sjevelazco/Ligustrum_lucidum_modeling/blob/main/R/02_Clean_coordinate_Env_filtering.R

#Set distances to be tested (in km)
d <- c(0, 5, 6, 7, 8, 9, 10, 15, 20, 25)

#Get species
spp <- unique(df$species)

#Make cluster
parallel::detectCores()
cl <- makeCluster(1)
clusterExport(cl, varlist= c("spp", "d", "occ", "var_wrap", "df",
                             "prepare_var"), #Get all objects created untill now
              envir=environment())
clusterEvalQ(cl, {
  library(dplyr)
  library(terra)
  library(flexsdm)
  library(pbapply)
  library(moranfast)
})

#Looping with pbapply
moran_df <- pblapply(seq_along(spp), function(i) {
  tryCatch({ #Avoid errors
  #Get species i
    sp <- spp[i]
    print(sp)
    occ <- df %>% filter(species == sp)
    
    #Unwrap var
    var <- unwrap(var_wrap)
    
    #Create pseudo-m: a minimum convex polygon with a buffer of 1500km
    #We set this initial buffer to decrease the time of filtering
    ca <- flexsdm::calib_area(occ, x = "x", y = "y",
                              method = c('bmcp', width=1500*1000),
                              crs = "+init=epsg:4326")
    
    #Run PCA of variables in the pseudo_M and select axis that explain 90%
    var_sp <- prepare_var(variables = var, calib_area = ca, do_PCA = TRUE,
                          var_portion = 90)
    
    
    #Filter using distances
    filtered <- suppressMessages(lapply(d, function(x){
      occfilt_geo(
        data = occ,
        x = "x",
        y = "y",
        env_layer = var_sp,
        method = c("defined", d = x),
        prj = crs(var_sp)
      )
    }))
    
    
    #Rename list with distances
    names(filtered) <- d
    
    # Calculate spatial autoccorelation (Moran I)
    imoran <- lapply(seq_along(filtered), function(x){
      tryCatch({ #Avoid errors
      coord <- filtered[[x]] %>% dplyr::select(x, y)
      data <- data.frame(terra::extract(var_sp, coord, ID = FALSE))
      coord <- coord[!is.na(data[,1]),]
      data <- data[!is.na(data[,1]),]
      imoran_x <- apply(data, 2, function(x)
        moranfast(x, coord$x, coord$y)$observed)
      },
      error=function(e) NULL) #Avoid errors
    })
    names(imoran) <- d
    
    #Convert list to dataframe
    imorandf <- do.call("rbind", imoran) %>% as.data.frame() %>% 
      mutate(Distance = d, .before = 1)
    
    
    #Get mean of imoran across PCA variables
    imorandf <- imorandf %>% dplyr::mutate(
      median_moran=apply(imorandf[, names(var_sp)], 1, median))
    
    #Get number of records remained in each distance
    imorandf$n_filtered <- sapply(filtered, nrow)

    #Put name of the specie in the dataframe and total number of records
    imorandf <- imorandf %>% mutate(species = sp, .before = 1) %>% 
      mutate(all_records = nrow(occ))
    #Propotion of lost records
    imorandf$prop_lost <- (imorandf$all_records - imorandf$n_filtered)/imorandf$all_records
    
    #Filtering distances: select lower autocorrelation (first quantile) which keeps the maximum number of occurrences
    finalfilter <- imorandf %>%
      filter(Distance > 0) %>% 
      mutate(median_moran = round(median_moran, 2)) %>% #Keep only 2 decimal places
      dplyr::filter(median_moran<=quantile(median_moran)[2]) %>% # Select 25th lower spatial autocorrelation
      dplyr::filter(n_filtered==max(n_filtered)) %>% # Select distance with higher number of records
      dplyr::sample_n(1) # Select a random value if more than one is selected
    
    #Get final points
    final_points <- filtered[[finalfilter$Distance %>% as.character()]]
    
    #Return final points
    return(final_points)
    
     },
    error=function(e) NULL) #Avoid errors
}, cl = cl)
#stop cluster
stopCluster(cl)

#Convert list to dataframe
occ_filt <- bind_rows(moran_df)

#Save
write.csv(occ_filt, "Data_example/Occurrences.csv", row.names = F)

# ####Plot maps of species####
# library(ggplot2)
# library(geobr)
# library(tidyterra)
# 
# #Import altitude raster
# dem <- rast("DEM/DEM.tiff")
# 
# #Import Brasil map
# br <- read_state(code_state = "all") %>% vect()
# #Import south america map
# sa <- vect("Vetores/South_America.shp")
# 
# 
# ####Looping####
# dir.create("Selected_species")
# dir.create("Selected_species/Maps")
# #Get species
# selected_spp <- sel_spp$species
# 
# #Unpack var
# var2 <- unwrap(var_wrap)
# 
# #Test
# i <- selected_spp[1]
# 
# pblapply(selected_spp, function(i){
# #Get species
# sp <- i
# d_sp <- sel_spp %>% filter(species == sp) %>% pull(Distance)
# occ_sp <- df200 %>% filter(species == sp)
# 
# 
# #Filter
# filtered_sp <-  occfilt_geo(
#   data = occ_sp,
#   x = "x",
#   y = "y",
#   env_layer = var2,
#   method = c("defined", d = d_sp),
#   prj = crs(var2)
# )
# 
# #Save filtered records
# write.csv(filtered_sp,
#           paste0("Selected_species/", sp, ".csv"), row.names = F)
# 
# #Cut map
# ca_sp <- calib_area(filtered_sp, x = "x", y = "y", 
#                     method = c('bmcp', width= 250*1000), crs = crs(dem))
# dem_sp <- crop(dem, ca_sp, mask = TRUE)
# 
# # Convert DEM to data frame
# dem_df <- as.data.frame(dem_sp, xy = TRUE, na.rm = T)
# 
# #Get limits to plot
# bb <- ext(dem_sp)
# 
# # Create base map
# map <- ggplot() +
#   geom_tile(data = dem_df, aes(x = x, y = y, fill = elevation_1KMmn_SRTM)) +
#   scale_fill_etopo() +
#   geom_sf(data = sa, linewidth = 0.2, colour = "grey40", alpha = 0) +
#   geom_sf(data = br, linewidth = 0.2, colour = "grey40", alpha = 0) +
#   new_scale_fill() +
#   geom_point(data = occ_sp, aes(x = x, y = y), size = 1.75,
#              shape = 21, colour = "black", fill = "black") +
#   geom_point(data = filtered_sp, aes(x = x, y = y), size = 1.75,
#              shape = 21, colour = "black", fill = "yellow", alpha = 0.75) +
#   coord_sf(xlim = c(bb[1] + 0.5, bb[2] - 0.5),
#            ylim = c(bb[3] + 0.5, bb[4] - 0.5),
#            expand = F) +
#   xlab("Longitude") + ylab("Latitude") +
#   ggtitle(sp,
#           subtitle = paste0(nrow(occ_sp), " records in total\n",
#                             nrow(filtered_sp), " filtered records")) +
#   metR::scale_x_longitude(ticks = 2) + metR::scale_y_latitude(ticks = 2) +
#   theme_classic() +
#   theme(legend.position = "none",
#         panel.border = element_rect(colour = "black", fill=NA, size=2),
#         panel.background = element_rect(fill = 'aliceblue', colour = NA),
#         axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8),
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5))
# #map
# ggsave(paste0("Selected_species/Maps/", sp, ".png"), map, dpi = 300, 
#        width = 7, height = 7)
#     
# 
# })
