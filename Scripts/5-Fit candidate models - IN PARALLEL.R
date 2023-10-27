##############################
#### FIT CANDIDATE MODELS ####
####      IN PARALELL     ####
#############################

#Load packages
library(terra)
library(foreach)
library(pbapply)
library(parallel)
library(enmpa)

#Import functions from Github
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/calibration_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/helpers_calibration_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/calibration_grid_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/glmnet_mx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/helpers_glmnetmx.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/evaluation_functions.R")
source("https://raw.githubusercontent.com/wevertonbio/pre_kuenm2/main/Functions/select_best_models.R")

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
                             "evaluation_functions", "select_best_models"),
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
    
    #Get swd file
    data <- readRDS(file.path(sp_dir, "Data_to_model.RDS"))
    
    #Let's create a grid of formulas, combining differents sets of variables, features and regularization multipliers
    f_grid <- calibration_grid_glmnetmx(swd = data$calibration_data, #Swd dataframe
                                        x = "x", y = "y", #If swd dataframe cointains x and y, set here; if not, use NULL
                                        min.number = 2, #Minimum number of variables in each set of variavles
                                        categorical_var = "SoilType", #Identify categorical variables, if exists
                                        features = c("l", "lq", "lqp"), #Features
                                        regm = c(0.1, 1, 5)) #Regularization multipliers (only 3 because it's an example)
    #See how many candidate models will be testes
    nrow(f_grid)
    
    #Fit candidate models
    #USe 70% of the available cores
    ncores <- round(parallel::detectCores()* 0.7, 0)
    
    m_sp <- calibration_glmnetmx(data = data, #Data prepared with the script 4
                                 formula_grid = f_grid, #Grid with formulas
                                 test_convex = TRUE, #Test concave curves in quadratic models?
                                 parallel = TRUE,
                                 ncores = ncores,
                                 progress_bar = TRUE, #Show progress bar? Only works if parallelType = "doSNOW"
                                 write_summary = FALSE, #Write candidate evaluations?
                                 out_dir = NULL, #Name of the folder to write candidate evaluations
                                 parallel_type = "doSNOW",
                                 return_replicate = TRUE,
                                 omrat_thr = c(5, 10),
                                 omrat_threshold = 10,
                                 skip_existing_models = FALSE, #Only works if writeFiles = TRUE
                                 verbose = TRUE)
    
    #Save candidate models
    saveRDS(m_sp,
            file.path(sp_dir, "candidate_results.RDS"))
  },
  error=function(e) NULL) #Avoid errors
}, cl = cl)
#stop cluster
stopCluster(cl)
