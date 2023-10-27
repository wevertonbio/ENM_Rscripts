##########################
#### Fit final models ####
####    IN PARALELL   ####
##########################

#Load packages
library(terra)
library(foreach)
library(pbapply)
library(parallel)

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


#Make cluster
parallel::detectCores()
cl <- makeCluster(4)
#Get necessary objects and send to nodes
clusterExport(cl, varlist= c("spp", "calibration_glmnetmx",
                             "helpers_calibration_glmnetmx",
                             "calibration_grid_glmnetmx",
                             "glmnet_mx", "helpers_glmnetmx",
                             "evaluation_functions", 
                             "fit_selected_glmnetmx", "part_data",
                             "predict_selected_glmnetmx"),
              envir=environment())
#Send necessary package to nodes
clusterEvalQ(cl, {
  library(terra)
  library(enmpa)
})

#Start looping
pblapply(seq_along(spp), function(i) {
  tryCatch({ #Avoid errors
    #Get species i
    sp <- spp[i]
    
    #Get species directory
    sp_dir <- file.path("Data_example/Models", sp)
    
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
  },
  error=function(e) NULL) #Avoid errors
}, cl = cl)
#stop cluster
stopCluster(cl)




