##########################
#### Fit final models ####
##########################

#Load packages
library(terra)
library(foreach)

#Import functions from Github
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/calibration_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/helpers_calibration_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/calibration_grid_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/glmnet_mx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/helpers_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/evaluation_functions.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/fit_selected_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/part_data.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/predict_selected_glmnetmx.R")

#Get species
spp <- list.dirs("Data_example/Models/", recursive = F, full.names = F)

#Get species directory
sp_dir <- file.path("Data_example/Models", spp)

#Get candidate results
candidate_results <- readRDS(file.path(sp_dir, "candidate_results.RDS"))
#See best models
#candidate_results$selected_models

#Fit best model with replicates (10 subsamples)
res_subsample <- fit_selected_glmnetmx(calibration_results = candidate_results,
                                       n_replicates = 10,
                                       rep_type = "subsample",
                                       train_portion = 0.7,
                                       write_models = TRUE, #Write files?
                                       file_name = file.path(sp_dir, "best_models"), #Name of the folder to write final models
                                       parallel = TRUE,
                                       ncores = 5,
                                       parallelType = "doSNOW",
                                       progress_bar = TRUE,
                                       verbose = TRUE) #Show messages?

####Predict####
#Get variables of the specie
var <- rast(file.path(sp_dir, "PCA_var.tiff"))

p <- predict_selected_glmnetmx(models = res_subsample,
                               spat_var = var,
                               write_files = TRUE,
                               write_replicates = FALSE,
                               out_dir = file.path(sp_dir, "Predictions"),
                               consensus_per_model = TRUE,
                               consensus_general = TRUE,
                               consensus = c("median", "range", "mean", "stdev"), #weighted mean
                               type = "cloglog",
                               overwrite = TRUE)
# #Plot consensus general
# plot(p$General_consensus$median)



