#### M simulation with grinnell package####
#See more in: https://github.com/fmachados/grinnell

#LOad packages
library(grinnell)
library(dplyr)
library(raster)
library(data.table)
library(sf)
library(terra)
library(flexsdm)

#Create directory to save M simulations
dir.create("Data_example/M_simulations/Shapefiles", recursive = T)

#Import all occurrences
occ <- read.csv("Data_example/Occurrences.csv")
unique(occ$species)

#Import variables from current and lgm and select only worldclim variables
current <- rast("Data_example/Current_Neotropic/Variables.tiff")
current <- current[[grepl("Bio", names(current))]]
#LGM
lgm <- rast("Data_example/LGMcc_Neotropic/LGMcc_2.5arcmin.tif")
lgm <- current[[grepl("Bio", names(lgm))]]
res(lgm)
#Remove variables 08, 09, 18 and 19
current <- current[[c(setdiff(names(current), c("Bio08", "Bio09", "Bio18", "Bio19")))]]
lgm <- lgm[[c(setdiff(names(lgm), c("Bio08", "Bio09", "Bio18", "Bio19")))]]

#To run faster, let's change to a coarser resolution
  #Convert to 10 arc min
c5 <- terra::aggregate(current, fact = 2, method = "bilinear")
res(c5)
lgm5 <- terra::aggregate(lgm, fact = 2, method = "bilinear")
res(lgm5)
  #Convert to 10 arc min
c10 <- terra::aggregate(current, fact = 4, method = "bilinear")
res(c10)
lgm10 <- terra::aggregate(lgm, fact = 4, method = "bilinear")
res(lgm10)

#Import neotropic islands - Necessary to check umber of polygons only in mainland
neotropic_islands <- vect("Vectors/Americas_Islands.gpkg")


#Create grid of combinations of parameters to increase KS
ks <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15) #Kernel spread
suit.t <- c(5, 1) #Suitability threshold
was_projected <- c(TRUE, FALSE) #Project to past?
#Combination grid
df_comb <- expand.grid(ks = ks,
                       suit.t = suit.t,
                       was_projected = was_projected)

# In this grid, the dispersal capability of the species is increasing. We are going to try these combinations until we meet the following requirements:
# All records must be inside the 'M'
# In the mainland, there must be only one continuous polygon.

#Remove points inside same pixel
#If has less than 200 records, filter using 5arc-min
if(nrow(occ) <= 200) {
  c.var <- c5
  lgm.var <- lgm5
}

#If has more than 200 records, filter using 10arc-min
if(nrow(occ) > 200) {
  c.var <- c10
  lgm.var <- lgm10
}

#Remove records with NA values in the worlclim variables
df.sp <- extract(c.var$Bio01, occ[,c("x","y")], cells=TRUE, ID = FALSE)
df.sp <- cbind(occ, df.sp)
df.sp <- df.sp %>% distinct(cell, .keep_all = TRUE) %>% na.omit()
df.sp <- df.sp %>% dplyr::select(species, longitude = x, latitude = y)

#Plot map, if you want:
# pts <- df.sp %>% st_as_sf(., coords = c(x = "x", y = "y"), crs = 4326)
# mapview(pts)

#Get spname - Use this only with a dataframe with a unique specie
sp.name <- occ[1,"species"]

#Create directory to save results
sp.path <- file.path("Data_example/M_simulations/Shapefiles", sp.name)
dir.create(sp.path)
dir.create(file.path(sp.path, "Simulations_details"))
dir.create(file.path(sp.path, "M_OK"))
dir.create(file.path(sp.path, "M_NOT_OK"))


#Set initial M: a minimum convex polygon with a buffer of 1500km
#We set this initial buffer to decrease the time of simulation
ca <- flexsdm::calib_area(df.sp, x = "longitude", y = "latitude",
                 method = c('bmcp', width=1500*1000),
                 crs = "+init=epsg:4326") %>% st_as_sf()

#Cut variables and transform in stack raster
c.var <- crop(c.var, ca, mask = TRUE) %>% raster::stack() #Current
lgm.var <- crop(lgm.var, ca) %>% raster::stack() #LGM
#Make sure they have the same extent
extent(c.var) <- extent(lgm.var)


#Simulate M
for (z in 1:nrow(df_comb)) {
  message("\nTesting combination ", z)
  #Get the combination of parameters i
  df_comb_i <- df_comb[z,]
  ks <- df_comb_i$ks
  suit_t <- df_comb_i$suit.t
  to_project <- df_comb_i$was_projected
  
  #M_simulationR
  m <- try(M_simulationR(data = df.sp,
                         current_variables = c.var,
                         starting_proportion = 0.5,
                         sampling_rule = "random",
                         barriers = NULL,
                         scale = TRUE,
                         center = TRUE,
                         project = to_project,
                         projection_variables = lgm.var,
                         dispersal_kernel = "normal",
                         kernel_spread = ks, #kernel spread
                         max_dispersers = 5,
                         suitability_threshold = suit_t, #Ellipsoid threshold
                         replicates = 10,
                         dispersal_events = 25, #One dispersed event each 40 years
                         access_threshold = 5,
                         simulation_period = 50,
                         stable_lgm = 7,
                         transition_to_lgm = 100,
                         lgm_to_current = 7,
                         stable_current = 13,
                         scenario_span = 1,
                         out_format = "GTiff",
                         set_seed = 42,
                         write_all_scenarios = F,
                         output_directory = sp.path)
  )
  #Delete all results that cause errors
  all_files <- list.files(sp.path, recursive = F,
                          pattern = ".tif$|accessible_area_M|png", full.names = T)
  unlink(all_files, recursive = T, force = T)
  unlink(file.path(sp.path, "PCA_results"), recursive = T, force = T)
  unlink(file.path(sp.path, "Suitability_results"), recursive = T, force = T)

  #Check if all points fall inside M
  if(class(m) == "list") {
    #Disaggregate M
    new_m <- disagg(vect(m$A_polygon))
    new_m$binary <- 1:length(new_m$binary)

    #See where the occurrences falling
    occ_m <- terra::extract(new_m, df.sp[,2:3])
    occ_m <- cbind(df.sp, "feature" = occ_m$binary)

    #Remove features without occurrences 
    m_final <- subset(new_m, new_m$binary %in% unique(occ_m$feature))
    crs(m_final) <- "+init=epsg:4326"
    
    #To check number of polygons, does not include islands and put buffer of 5 km
    m_no_islands <- terra::erase(m_final, neotropic_islands) %>% 
      buffer(., width = 5*1000) %>% terra::aggregate() %>% disagg()

    #Calculate number of polygons now
    n_pol <- length(m_no_islands)
    #Are there points outside the M?
    n_occ_outside <- sum(is.na(occ_m$feature))
    outside_occ <- ifelse(any(is.na(occ_m$feature)), "Yes", "No")
  }

  #If the result is a error, inform that there is multipolygons and there is occurrences outside the M (to run analysis with the next combination of parameters)
  if(class(m) == "try-error") {
    #Delete all results
    try(dev.off())
    unlink(sp.path, recursive = T, force = T)
    n_pol <- 2
    outside_occ <- "Yes"
  }
  
  #Write vector
  if(class(m) == "list") {
    ####Save results to check####
    #Identify if there is a problem
    if(n_pol == 1 & n_occ_outside == 0) {p <- sp.name} #Sp name if there is no problem
    if(n_pol > 1 & n_occ_outside == 0) {p <- "Multiple_Polygons"}
    if(n_pol == 1 & n_occ_outside > 0) {p <- "Outside_Points"}
    if(n_pol > 1 & n_occ_outside > 0) {p <- "Multiple_Polygons_and_Outside_Points"}
    
  if(p == sp.name){
    writeVector(m_final, filename = paste0(sp.path,
                                           "/M_OK/",
                                           p, ".gpkg"),
                overwrite = TRUE)
  } else {
    writeVector(m_final, filename = paste0(sp.path,
                                           "/M_NOT_OK/",
                                           z, "_", p, ".gpkg"),
                overwrite = TRUE)
  }
    
    #Identify when specie is not projected
    if (to_project == TRUE) {
      simul_period = 50
    }
    if (to_project == FALSE) {
      simul_period = "Only_Current"
    }
    
    #Save results in a dataframe
    sp.info <- data.frame("Combination" = z,
                          "species" = sp.name,
                          "N_records" = nrow(df.sp),
                          "Outside_M" =  n_occ_outside,
                          "N-polygons" = n_pol,
                          "kernel_spread" = ks,
                          "suitability_threshold" = suit_t,
                          "replicates" = 10,
                          "dispersal_events" = 25, #One dispersed event each 40 years
                          "access_threshold"= 5,
                          "simulation_period"= simul_period,
                          "stable_lgm" = 7,
                          "transition_to_lgm" = 100,
                          "lgm_to_current" = 7,
                          "stable_current" = 13,
                          "scenario_span" = 1,
                          "starting_proportion" = 0.5)
    #Write dataframe
    write.csv(sp.info, paste0(sp.path,
                              "/Simulations_details/Model_info_",
                              z, ".csv"),
              row.names = F)
      }
  
  #If the function generate a M with a unique polygon including all occurrences, stop function :)
  if(n_pol == 1 & outside_occ == "No") {
    break
  } else {
    message("Combination ", z, " failed because of: ", p)
  }

} #End of for in