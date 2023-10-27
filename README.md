# ENM_Rscripts

This repository serves to store R scripts to prepare data and develop Ecological Niche Models using glmnet (like maxnet)

Authors: Weverton Carlos Ferreira Trindade\*, Márcia Marques

Contact: [wevertonf1993\@gmail.com](wevertonf1993@gmail.com)

## Description

The 'Scripts' folder contains the code to run the analysis, including both example scripts and scripts for run the analysis in parallel (indicated by "- IN PARALLEL""). The 'Data_example' folder contains the data necessary to run the examples with *Araucaria angustifolia*.

## Workflow

### 1 - Prepare environmental variables

We obtained bioclimatic variables from [WorldClim](https://www.worldclim.org/data/worldclim21.html), topographic variables from [EarthEnv](https://www.earthenv.org/topography), and edaphic (soil) variables from [SoilGrids](https://files.isric.org). Additional environmental variables from various sources are referenced in the script.

To ensure consistency in spatial extent and resolution, we used this script to crop all variables to the Neotropic region and resampled them to a 2.5 arc-minute resolution. The script also handled the filling of missing values in the soil variables.

With the script **1.1-Prepare variables to simulate M with Grinnell**, we crop and fix the [variables from the Last Glacial Maximum obtained from WorldClim](http://www.worldclim.com/past) (GCM: CCSM), which are necessary to simulate the M of the species.

### 2 - Filter occurrences by distance and spatial autocorrelation

To mitigate the effects of spatial sampling biases, we adapt the `occfilt_geo` function from the [flexsdm package](https://sjevelazco.github.io/flexsdm/index.html). We examined various distance thresholds between pairs of points, and for each distance, we computed spatial autocorrelation using a Moran's semivariogram. Instead of opting for the distance that reduces Moran's I to values below 0.1, we chose the distance that yielded the lowest average spatial autocorrelation (25th percentile) while maximizing the number of occurrences.

### 3 - M simulation with grinnell package

M represents the areas explored by the species of interest throughout its recent history. To delineate M, we employed a simulation approach that accounts for processes of dispersal, colonization, and extinction in both constant current and glacial-interglacial climate change scenarios. This approach is implemented within the [grinnell package](https://escholarship.org/uc/item/8hq04438).

The parameters of the simulations depend on the biology of the species but also on the available occurrences and the resolution of environmental variables. To address this, we created a grid of parameter combinations for the simulations. In summary, each simulation increases the species' dispersal capability. The simulation process halts when the simulated M meets the following requirements:
* All records must fall within the 'M' boundary.
* On the mainland, there should be only one continuous polygon, ensuring there are no discontinuous M areas.

### 4 - Prepare data for modelling
To reduce multicollinearity among the bioclimatic variables, we conduct a Principal Component Analysis (PCA) and spatialize the principal components that collectively explain 95% of the total variance. The topographic and soil variables were incorporated into these principal components. For the occurrences and background points, we generated 10,000 background points and partitioned the data into 4 k-folds.

### 5 - Fit candidate models
Each candidate model represents a unique combination of variables, features (linear, quadratic, and product), and regularization multipliers (0.1, 1, 3, 5). For each candidate model, we assessed AIC, pROC, omission rate at 10%, and the presence of concave curves. Our selection of the best models is based on the following criteria:

* Significant pROC.
* Omission rate for test points below 10% (with training points also at 10%).
* Absence of concave curves.
* Lowest AIC (delta AIC less than or equal to 2)

### 6 - Fit and predict best models
The final results were derived from 10 replicates of the best model, each using a subsample method with 70% of the records assigned as training and 30% as test.