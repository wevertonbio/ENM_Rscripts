###############################
#### PREPARE DATA TO MODEL ####
###############################

#Load packages
library(terra)

#Import functions from github
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/prepare_var.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/prepare_data.R")

#Create a folder to save the models
dir.create("Data_example/Models")

#Load the variables
var <- rast(list.files("Data_example/Current_Neotropic/", full.names = TRUE))
names(var)
#Delete Altitude variables that mix temperature and precipitation:
var <- var[[setdiff(names(var), c("Bio08", "Bio09", "Bio18", "Bio18",
                                  "Altitude"))]]
names(var)
#Lets run the pca only with the worldclima variables
#To do that, we need to identify the raw version of the others variables
others_vars <- names(var)[!grepl("Bio", names(var))]

#For each specie, we are going to cut the variables using the M, run a PCA and get the x axis that explain 95% of the variance
#Import occurrences
occ <- read.csv("Data_example/Occurrences.csv")
sp <- occ[1,1]

#Create directory to model specie
sp_dir <- file.path("Data_example/Models", sp)
dir.create(sp_dir)


#Get M of the specie
m_sp <- vect(paste0("Data_example/M_simulations/Shapefiles/",
            sp, "/M_OK/", sp, ".gpkg"))

#Function to prepare var
pca_var <- prepare_var(variables = var, #Variables to fit the model
            calib_area = m_sp, #M or another mask to crop the variables
            do_PCA = T, #Do PCA?
            exclude_from_PCA = others_vars, #Variables to exclude from PCA
            var_portion = 95, #Get the x axis that explain 95% of the variance
            writeFiles = T, #Write rasters?
            write_PCA_files = T, #Write PCA files? (Model, rotation and importance)
            out_dir = sp_dir, #Outuput directory
            overwrite = T, #Overwrite when rasters exists?
            project_PCA = F, #Project PCA?
            proj_folder = NULL, #Folder with variables to project PCA
            return_proj = FALSE, #Return projections as a object in the environment?
            verbose = TRUE) #Show messages?

#Now, let's prepare the table to models, which is a file with the values of the variables extracted to each occurrence record of the species, and to n (10.000) background points
data <- prepare_data(occ = occ,
        species = "species",
        x = "x",
        y = "y",
        spat_variables = pca_var,
        categorical_variables = "SoilType",
        nbg = 10000,
        kfolds = 4,
        include_xy = TRUE,
        write_files = T,
        file_name = file.path(sp_dir, "Data_to_model"),
        seed = 42,
        verbose = TRUE) #Output directory