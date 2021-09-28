```
#      /  /\        ___          /__/\   
#     /  /:/_      /  /\         \  \:\  
#    /  /:/ /\    /  /:/          \__\:\ 
#   /  /:/ /::\  /__/::\      ___ /  /::\
#  /__/:/ /:/\:\ \__\/\:\__  /__/\  /:/\:\
#  \  \:\/:/~/:/    \  \:\/\ \  \:\/:/__\/
#   \  \::/ /:/      \__\::/  \  \::/    
#    \__\/ /:/       /__/:/    \  \:\    
#      /__/:/ please \__\/      \  \:\   
#      \__\/ acknowledge your use\__\/   
#
```


Increment-averaged kriging for 3D prediction of soil properties.

Example datasets from Edgeroi package dataset are given to illustrate using
this model. However, any datasets can be used provided the required column
names are present. Please refer to the EdgeroiValidationData and EdgeroiFitData
documentation for more information on this.

Usage: 

```r
Cubistdata <- CubistIAK(fit_data = EdgeroiFitData,
                validate_data = EdgeroiValidationData, 
                spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
```

OR:

```r
Splinedata <- SplineIAK(fit_data = EdgeroiFitData,
                validate_data = EdgeroiValidationData, 
                spatialCovs = c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
```

These return fitted model results, and paramaters that have been assigned during
model run. A special mention on modelX within the lmm.fit.selected variable,
which holds the fitted model which can be used for prediction purposes (below).

To run further plots using model output

```r
RunPlots(fit = Cubistdata$lmm.fit.selected, 
            InputParamatersList = Cubistdata$InputParamatersList,
            chooseToPlot = c(1,2,3,4,5,6))
```

# Example workflow on Prediction using IAK results with Edgeroi data.
This also reflects currently required variables in the data the package expects.
```r
FitData <- Iak3dSIH::EdgeroiFitData
ValidateData <- Iak3dSIH::EdgeroiValidationData

predictors <- c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
required <- c('x','y','lowerDI','upperDI','z','profIDFit')
exclude <- c('z')

#Extra work to set up prediction dataset for demo purposes.
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

#using model result as prediction on new data
model_tree_pred <- predict(model_tree, predict_using_IAK_model) 
## Test set RMSE from using prediction
sqrt(mean((model_tree_pred - test_resp)^2))
```


