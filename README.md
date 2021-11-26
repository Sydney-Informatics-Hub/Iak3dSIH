![Sydney Informatics Hub Logo](U_Sydney_Informatics_Hub_mono_300_152.png)
```
      /  /\        ___          /__/\   
     /  /:/_      /  /\         \  \:\  
    /  /:/ /\    /  /:/          \__\:\ 
   /  /:/ /::\  /__/::\      ___ /  /::\
  /__/:/ /:/\:\ \__\/\:\__  /__/\  /:/\:\
  \  \:\/:/~/:/    \  \:\/\ \  \:\/:/__\/
   \  \::/ /:/      \__\::/  \  \::/    
    \__\/ /:/       /__/:/    \  \:\    
      /__/:/ please \__\/      \  \:\   
      \__\/ acknowledge your use\__\/   
```

# Increment-averaged kriging for 3D prediction of soil properties.
Rewrite of Tom Orton's [iak3d package](https://github.com/ortont/iak3d) by Kristian Maras of the [Sydney Informatics Hub](https://www.sydney.edu.au/informatics-hub)

IAK is a framework that accounts for vertical and spatial correlation between observations. Regression parameters are obtained from either Cubist or Spline models. Model outputs are returned and optional plots can also be generated from the outputs.

Example datasets from Edgeroi package dataset are given to illustrate using
this model. However, any datasets can be used provided the required column
names are present (mentioned below). 

## Installation


1. In R, install the package 

```r
#install.packages("devtools") #install devtools if required

devtools::install_github("Sydney-Informatics-Hub/Iak3dSIH")

```

3. If this doesn't work you can instead download the package as a zip file from here, then install it locally from source, with

```r
install.packages("~/Downloads/Iak3dSIH-master.zip", 
                 repos=NULL, 
                 type="source")
```
as long as you change the path to match where the downloaded zip file is on your system.


## Usage: 

To use the package

```r
library(Iak3dSIH)
```

```r
Cubistdata <- CubistIAK(fit_data = EdgeroiFitData,
                validate_data = EdgeroiValidationData, 
                spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
```

OR:

```r
Splinedata <- SplineIAK(fit_data = EdgeroiFitData,
                validate_data = EdgeroiValidationData, 
                spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
```

 These return fitted model results, and paramaters that have been assigned during
model run. At the moment, dIMidPts is a required parameter with the spatialCovs argument along with column in the fit and/or validate data.

A special mention on modelX within the lmm.fit.selected variable,
which holds the fitted model which can be used for prediction purposes (below). Information on negative log-likelihood and semivariogram parameters are also given lmm.fit.selected    variable. zkVal and vkVal variables are also returned as outputs, where zkVal is the horizontal predictions made using the validation set, and vkVal are the vertical predictions made using the validation set.

Along with the above necessary parameters for running spline and cubist, optional parameters include : 
proportion_crossvalidation - Proportion of data used for cross validation. Between 0 and 1, and defaults to 0.7.
layerMidPts - User defined mid points for the validation statistics. 

An example usage with option parameters:
```r
Cubistdata <- CubistIAK(fit_data = EdgeroiFitData,
                        validate_data = EdgeroiValidationData, 
                        spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'), 
                        proportion_crossvalidation = 0.7, 
                        layerMidPts = c(0.025, 0.1, 0.225, 0.45, 0.8, 1.5))
                        
```

To run further plots saved in working directory using the IAK model outputs:

```r
RunPlots(fit = Cubistdata$lmm.fit.selected, 
            InputParamatersList = Cubistdata$InputParamatersList,
            chooseToPlot = c(1,2,3,4,5,6),
            distance_semivariogram = 20,
            bins = 500,
            depth_interval_plots = data.frame('dL' = c(0,15 , 30 , 60)/100 , 'dU' = c(15,30 , 60 , 90)/100))
```

## Example workflow on Prediction using IAK results with Edgeroi data.
This also reflects currently required variables in the data the package expects.
```r
FitData <- Iak3dSIH::EdgeroiFitData
ValidateData <- Iak3dSIH::EdgeroiValidationData

predictors <- c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
required <- c('x','y','lowerDI','upperDI','z','profIDFit')
exclude <- c('z')
```

Extra work to set up prediction dataset for demo purposes.
```r
in_train_set <- sample(1:nrow(FitData), floor(.8*nrow(FitData)))
train_pred <- FitData[ in_train_set, c(predictors,required)]
train_resp <- train_pred$z # repose variable for training
test_pred  <- FitData[-in_train_set, c(predictors,required)]
test_resp  <- test_pred$z #response variable
test_pred <- test_pred[-which(names(test_pred) %in% exclude)] 

model_IAK <- Iak3dSIH::CubistIAK(fit_data = train_pred,
                                 validate_data = ValidateData,
                                 spatialCovs = predictors)

model_tree <- model_IAK$lmm.fit.selected$modelX
summary(model_tree)
predict_using_IAK_model <- test_pred[predictors]
```

Using model result as prediction on new data
```r
model_tree_pred <- predict(model_tree, predict_using_IAK_model) 
```

Test set RMSE from using prediction
```r
sqrt(mean((model_tree_pred - test_resp)^2))
```

## Data Inputs
Data inputs (either Fitting or Validation data) at the minimum expect some column names being:

* Dataframe columns
    + x : x dimension (in Km)
    + y : y dimension (in Km)
    + lowerDI : lower Depth Interval of soil range in metres
    + upperDI : upper Depth Interval of soil range in metres
    + z:  response variable for valid set
    + dIMidPts :  Mid point of Depth Interval

Other columns can be added that relate to covariates. These column names are expected to be referenced with the spatialCovs argument.


