##############################################################
# Main package improvements
#1. Turn Edgeroi into a portable package
#2. Run Edgeroi with input parameters in more modular way, including model selection and plots using outputs from model run.
#3. Expand to support different data inputs.

#Notes on usages: 
#Splinedata <- SplineIAK(fit_data = EdgeroiFitData,validate_data = EdgeroiValidationData, spatialCovs = c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
#Cubistdata <- CubistIAK(fit_data = EdgeroiFitData,validate_data = EdgeroiValidationData, spatialCovs = c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
#RunPlots(fit = Cubistdata$lmm.fit.selected, InputParamatersList = Cubistdata$InputParamatersList,c(1,2,3,4,5,6))

##############################################################


#' Represents the underlying data and required format needed to run IAK.
#' Example data has been included in the package, sourced from the Edgeroi Package (i.e data(edgeroi) and data(edgeroiCovariates), 
#' named Iak3dSIH::Uniform_Data_Edgeroi. Any data can be inputted provided the structure is as follows:
#'
#' Model Fit set includes:
#' cFit : Coordinates (in Km) of calib set, in two dimensions x y. 
#' dIFit : Depth intervals of fit set (metres) in two dimensions representing upper and lower soil sample range.
#' covsFit : covariates of fit set, with column names that represent covariates aligned with param spatialCovs. Within this, dIMidPts is Required (represents the mid depth i.e. average between the upper and lower depth of an observation.)
#' zFit : response variable of fit set
#' profIDFit : profile ID of fit set

#' Validation Set - Required if validation flag is TRUE, includes:
#' cVal :  Coordinates (in Km) of valid set
#' dIVal : depth intervals of valid set (in metres)
#' covsVal :  Covariates for valid set, with column names that represent covariates aligned with param spatialCovs, Within this, dIMidPts Required (represents the mid depth i.e. average between the upper and lower depth of an observation.)
#' zVal :  response variable for valid set
#' profIDVal :  profile ID representing unique strings or numbers
#' rList :  raster list of covariate - OPTIONAL
#' @name Uniform_Data_Edgeroi
#' @docType data
#' @keywords data
NULL

#' Validation data example for Edgeroi Area. Expected to be DataFrame with column names:
#' x : x dimension (in Km)
#' y : y dimension (in Km)
#' lowerDI : lower Depth Interval of soil range in metres
#' upperDI : upper Depth Interval of soil range in metres
#' z:  response variable for valid set
#' dIMidPts :  Mid point of Depth Interval
#' Other columns relate to covariates and column names should reflect spatialCovs argument
#' @name EdgeroiValidationData
#' @docType data
#' @keywords data
NULL

#' Fit data example for Edgeroi Area. Expected to be DataFrame with column names:
#' x : x dimension (in Km)
#' y : y dimension (in Km)
#' lowerDI : lower Depth Interval of soil range in metres
#' upperDI : upper Depth Interval of soil range in metres
#' z:  response variable for valid set
#' dIMidPts :  Mid point of Depth Interval
#' Other columns relate to covariates and column names should reflect spatialCovs argument

#' @name EdgeroiFitData
#' @docType data
#' @keywords data
NULL


FitSplineModel <- function(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList,spatialCovs) {
  #################################################################################################
  ### set knots for sdfd spline fn (if used)
  #################################################################################################
  ### alternatively, set up for fitting a spline model.
  ### include interactions between depth and spatial covariates (but here::here not between different spatial covariates)
  ###################################################################################
  print("FitSplineModel is engaged")
  modelX <- list('type' = 'gam2')
  scaleCovs <- TRUE
  nIntKnotsd <- 4 # number of internal knots for the spline function (nat spline, clamped to have grad=0 at upper bdry) of depth; if this is complex enough, probably no need for the depth component in prod-sum covariance model
  nIntKnotss <- 4 # number of internal knots for the spline functions (nat spline, clamped to have grad=0 at upper and lower bdries) of covariates
  
  ### don't include depth here::here.   
  
  if(scaleCovs){
    print("Scaled covariates created")
    ### to work with scaled covariates
    spatialCovs <- paste0(spatialCovs , '_SCALED')
  }else{
    ### to work with unscaled (raw) covariates
    spatialCovs <- spatialCovs # no change here::here
  } 
  
  ### add any scaled variables (_SCALED) to covs dfs...
  tmp <- addScaledCovs2df(dfFit = covsFit , dfPred = covsVal , covNames = spatialCovs) #new but deleted
  covsFit <- tmp$dfFit #new to pass
  covsVal <- tmp$dfPred #new to pass
  q4BdryKnots <- c(0.05 , rep(0.01 , length(spatialCovs))) 
  q4BdryKnots <- cbind(q4BdryKnots , 1 - q4BdryKnots)
  nIntKnots <- c(nIntKnotsd , rep(nIntKnotss , length(spatialCovs)))
  sType <- c('nscug' , rep('nsclug' , length(spatialCovs)))
  modelX$listfefdKnots <- makelistfefdKnots(dfFit = covsFit , covNames = c('dIMidPts' , spatialCovs) , nIntKnots = nIntKnots , q4BdryKnots = q4BdryKnots , sType = sType)
  ### to include interactions between depth and the other spatial covariates...  
  modelX$incInts <- list('dIMidPts' , spatialCovs)
  modelX$intMthd <- 0
  XcnsTmp <- makeXcns(dfCovs = covsFit , dIData = dIFit , listfefdKnots = modelX$listfefdKnots , incInts = modelX$incInts , colnamesX = NULL , intMthd = modelX$intMthd) # intMthd = 1 for now. 0 = simpler
  modelX$colnamesX <- colnames(XcnsTmp)
  
  rm(tmp , XcnsTmp , q4BdryKnots , nIntKnots , sType)
  return(list(modelX=modelX,cFit=cFit, dIFit=dIFit, covsFit=covsFit, zFit=zFit, profIDFit=profIDFit, cVal=cVal, dIVal=dIVal, covsVal=covsVal, zVal=zVal, profIDVal=profIDVal, rList=rList))
  
}

FitCubistModel <- function(paramaters, tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList,spatialCovs) {
  ##############################################################
  ### an algorithm to select number of rules for cubist model
  
  print("entering cubist model")
  #print(paramaters)
  if(is.na(paramaters$otherparamaters$nRules) | is.na(paramaters$otherparamaters$refineCubistModel)){
    if(is.na(paramaters$otherparamaters$nRules)){ nRulesVec <- seq(20) }else{ nRulesVec <- nRules }
    if(is.na(paramaters$otherparamaters$refineCubistModel)){ refineCubistModelVec <- c(FALSE , TRUE) }else{ refineCubistModelVec <- refineCubistModel }
    tmp <- selectCubistOptsXV(cFit = cFit , zFit = zFit , covsFit = covsFit , allKnotsd = paramaters$otherparamaters$allKnotsd , nRulesVec = nRulesVec , refineCubistModelVec = refineCubistModelVec)
    nRules <- tmp$nRules                            # new variable
    refineCubistModel <- tmp$refineCubistModel      # new variable even though initially set as boolean
    rmseMatList <- tmp$rmseMatList                  # new variable
    warningFlagFitList <- tmp$warningFlagFitList    # new variable
  }else{}
  
  ###############################################################
  ### fit or load the Cubist model... related to above
  
  
  if(paramaters$FitCubits){
    ### fit cubist model...
    print("fitCubistModelNow is run. ......................")
    
    cmFit <- Cubist::cubist(x = covsFit , y = zFit , committees = 1 , Cubist::cubistControl(rules = paramaters$otherparamaters$nRules))
    
    ### convert to des mtx

    tmp <- cubist2X(cubistModel = cmFit, dataFit = covsFit , zFit = zFit , profIDFit = profIDFit , allKnotsd = paramaters$otherparamaters$allKnotsd , refineCubistModel = paramaters$otherparamaters$refineCubistModel)
    #browser()
    cmFit <- tmp$cubistModel
    XFit <- tmp$X
    matRulesFit <- tmp$matRuleData
   # save(cmFit , file = paste0(dataDir , '/cmFit.RData'))
    
  }
  
  modelX <- cmFit
  modelX$type = "NOTgam2"
  # cfit something different now.......
  return(list(modelX=modelX,cFit=cFit, dIFit=dIFit, covsFit=covsFit, zFit=zFit, profIDFit=profIDFit, cVal=cVal, dIVal=dIVal, covsVal=covsVal, zVal=zVal, profIDVal=profIDVal, rList=rList))
}


LoadModelDirectly <- function(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) {
  modelX <- load(file = paste0(dataDir , '/cmFit.RData'))
  return(list(modelX=modelX,cFit=cFit, dIFit=dIFit, covsFit=covsFit, zFit=zFit, profIDFit=profIDFit, cVal=cVal, dIVal=dIVal, covsVal=covsVal, zVal=zVal, profIDVal=profIDVal, rList=rList))
  
}
#spatialCovs <- c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
LoadModel <- function(paramaters,spatialCovs) {
  #  input : list(FitCubits=paramaters$fitCubistModelNow, 
  #       useCubistForTrend = paramaters$useCubistForTrend, 
  #       LoadModel = paramaters$fitModelNow, #NEGATE THIS
  #       otherparamaters=paramaters$otherparamaters,
  #       data=tmp)
  #   output: ModelX plus parameters
  
  print("in LoadModel......................")
  
  tmp <- paramaters$data
  cFit <- tmp$cFit
  dIFit <- tmp$dIFit
  covsFit <- tmp$covsFit
  zFit <- tmp$zFit
  profIDFit <- tmp$profIDFit
  
  cVal <- tmp$cVal
  dIVal <- tmp$dIVal
  covsVal <- tmp$covsVal
  zVal <- tmp$zVal
  profIDVal <- tmp$profIDVal
  rList <- tmp$rList
  
  if(paramaters$FitCubits){
    print("doing cubist")
    ModelOutput <- FitCubistModel(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList,spatialCovs) #................
  }
  else if (paramaters$LoadModel) {
    print("loading directly from a file")
    #load directly from a file that is expected 
    ModelOutput <- LoadModelDirectly(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList)
  }
  else {
    print("doing spline")
    #spline model
    ModelOutput <- FitSplineModel(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList,spatialCovs) # first seperation into Spline and return stuff
  }
  return(ModelOutput)
}

LoadData <- function(paramaters,datafeedin,spatialCovs){
  
  ##############################################
  ### load the edgeroi dataset (from GSIF package) and put into format for iak3d...
  ##############################################
  print("now in loadData................................")
  
  tmp <- datafeedin
  print("Print constructed ModelOptions")
  ModelOptions <- list(FitCubits=paramaters$fitCubistModelNow, 
                       useCubistForTrend = paramaters$useCubistForTrend, 
                       LoadModel = !paramaters$fitModelNow, 
                       otherparamaters=paramaters$otherparamaters,
                       data=tmp)
  output <- LoadModel(ModelOptions,spatialCovs)
  
  return(output)
  
}

#' Run Iak3d project with building spline model. Use after running either SplineIAK OR CubistIAK functions.
#' Example: RunPlots(fit = Cubistdata$lmm.fit.selected, InputParamatersList = Cubistdata$InputParamatersList,c(1,2,3,4,5,6))
#' @param lmm.fit.selected Included in the result of running either Cubist or Spline Mode.
#' @param InputParamatersList Included in the result of running either Cubist or Spline Models. Attaches as input run to an output result and enables plots on the outputs for analysis.
#' @param chooseToPlot Select a maximum of six co-ordinate points to plot for inspection, example: c(1,2,3,4,5,6). Plot is calebrated for 6 plots. These are rrediction plots that are saved in plotVal4Plot.pdf. Upper limit Limited by unique number of cVal. 
#' @return saves plots in the working directory
#' @export
RunPlots <- function(fit = lmm.fit.selected, InputParamatersList = InputParamatersList,chooseToPlot = chooseToPlot) {
  lmm.fit.selected <- fit
  dIPlot <- data.frame('dU' = c(0 , 20 , 50 , 90 , 150 , 190)/100 , 'dL' = c(10 , 30 , 60 , 100 , 160 , 200)/100)
  hx <- seq(0 , 20 , 1)
  grDevices::pdf(file = paste0(getwd() , '/varioFitgam22.pdf'))
  tmp <- plotCovx(lmm.fit = lmm.fit.selected , hx = hx , dIPlot = dIPlot , addExpmntlV = TRUE , hzntlUnits = 'km')
  grDevices::dev.off()
  
  hdPlot <- seq(0 , 2 , 0.01)
  grDevices::pdf(file = paste0(getwd() , '/cordFit.pdf'))
  qwe <- plotCord(lmm.fit = lmm.fit.selected , hdPlot = hdPlot, vrtclUnits = 'm')
  grDevices::dev.off()
  
  dTmp <- seq(0 , 2 , 0.1)
  dIPlot <- data.frame('dU' = dTmp[-length(dTmp)] , 'dL' = dTmp[-1])
  grDevices::pdf(file = paste0(getwd() , '/covardFit.pdf'))
  qwe <- plotCovd(lmm.fit = lmm.fit.selected , dIPlot = dIPlot , vrtclUnits = 'm')
  grDevices::dev.off()
  
  ### plot of the variances...
  grDevices::pdf(file = paste0(getwd() , '/varComps.pdf'))
  dPlot <- seq(0 , 2 , 0.01)
  plotVarComps(lmm.fit = lmm.fit.selected , dPlot = dPlot)
  grDevices::dev.off()
  #Prediction Plots
  ModelOutput <- InputParamatersList
  LastSeperation(ModelOutput,lmm.fit.selected ,chooseToPlot)

}
RunValidation <- function(ModelOutput,dataDir,namePlot,lmm.fit.selected,rqrBTfmdPreds,constrainX4Pred,fnamezkVal,fnamevkVal) {
  
  nVal <- nrow(ModelOutput$cVal) #YES
  iU <- Matrix::which(!duplicated(ModelOutput$cVal)) #YES
  cValU <- ModelOutput$cVal[iU,,drop=FALSE] #NO
  covsValU <- ModelOutput$covsVal[iU,,drop=FALSE]
  zkVal <- vkVal <- NA * numeric(nVal) 
  for(i in 1:nrow(cValU)){
    iTmp <- Matrix::which(ModelOutput$cVal[,1] == cValU[i,1] & ModelOutput$cVal[,2] == cValU[i,2])
    tmp <- profilePredictIAK3D(xMap = cValU[i,,drop=FALSE] , covsMap = ModelOutput$covsVal[iTmp,,drop=FALSE] , dIMap = ModelOutput$dIVal[iTmp,,drop=FALSE] , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
    zkVal[iTmp] <- tmp$zMap
    vkVal[iTmp] <- tmp$vMap
  }
  
  
  dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
  tmp <- predictIAK3D(xMap = cValU , dIMap = dIPred , covsMap = covsValU , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
  
  zkProfPred <- tmp$zMap
  vkProfPred <- tmp$vMap
  pi90LkProfPred <- tmp$pi90LMap
  pi90UkProfPred <- tmp$pi90UMap
  
  ### calc and print val stats...
  tmp <- calcValStats(zVal = ModelOutput$zVal , dIVal = ModelOutput$dIVal , zkVal = zkVal , vkVal = vkVal , layerMidPts = c(0.025 , 0.1 , 0.225 , 0.45 , 0.8 , 1.5) , printValStats = TRUE)
  valStatsAllLayers <- tmp$valStatsAllLayers
  valStatsTot <- tmp$valStatsTot
  
  #Turned off for now - seperate plottting elsewhere
  #tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = ModelOutput$cVal , dIData = ModelOutput$dIVal , zData = ModelOutput$zVal ,
  #                         xPred = cValU , dIPred = dIPred , zPred = zkProfPred , pi90LPred = pi90LkProfPred , pi90UPred = pi90UkProfPred ,
  #                         zhatxv = zkVal , pi90Lxv = zkVal - 1.64 * sqrt(vkVal) , pi90Uxv = zkVal + 1.64 * sqrt(vkVal))
  
  

  
  return(list(zkVal=zkVal,vkVal=vkVal))
}

LastSeperation <- function(ModelOutput,lmm.fit.selected , rand6ForPlot,rqrBTfmdPreds = FALSE, constrainX4Pred = FALSE){
  
  #browser()
  #rand6ForPlot <- c(6 , 19 , 49 , 41 , 3 , 24) # For edg
  #rand6ForPlot <- c(6 , 8, 15 , 5 , 3 , 4) # For OOd
  i4PlotU <- Matrix::which(!duplicated(ModelOutput$cVal))[rand6ForPlot]
  cVal4PlotU <- ModelOutput$cVal[i4PlotU,,drop=FALSE]
  covsVal4PlotU <- ModelOutput$covsVal[i4PlotU,,drop=FALSE]
  
  iVal4Plot <- c()
  for(i in 1:nrow(cVal4PlotU)){
    iTmp <- Matrix::which(ModelOutput$cVal[,1] == cVal4PlotU[i,1] & ModelOutput$cVal[,2] == cVal4PlotU[i,2])
    iVal4Plot <- c(iVal4Plot , iTmp)
  }
  
  cVal4Plot <- ModelOutput$cVal[iVal4Plot,,drop=FALSE]
  dIVal4Plot <- ModelOutput$dIVal[iVal4Plot,,drop=FALSE]
  covsVal4Plot <- ModelOutput$covsVal[iVal4Plot,,drop=FALSE]
  zVal4Plot <- ModelOutput$zVal[iVal4Plot]
  
  zkVal4Plot <- vkVal4Plot <- NA * numeric(length(zVal4Plot))
  #browser() #for Ood at this point some nulls in dataframe cVal4PlotU
  for(i in 1:nrow(cVal4PlotU)){
    iTmp <- Matrix::which(cVal4Plot[,1] == cVal4PlotU[i,1] & cVal4Plot[,2] == cVal4PlotU[i,2])
    tmp <- profilePredictIAK3D(xMap = cVal4PlotU[i,,drop=FALSE] , covsMap = covsVal4Plot[iTmp,,drop=FALSE] , dIMap = dIVal4Plot[iTmp,,drop=FALSE] , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
    zkVal4Plot[iTmp] <- tmp$zMap 
    vkVal4Plot[iTmp] <- tmp$vMap 
  }
  
  dIPred <- cbind(seq(0 , 1.98 , 0.02) , seq(0.02 , 2 , 0.02))
  tmp <- predictIAK3D(xMap = cVal4PlotU , dIMap = dIPred , covsMap = covsVal4PlotU , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
  
  zkProfPred <- tmp$zMap
  vkProfPred <- tmp$vMap
  pi90LkProfPred <- tmp$pi90LMap
  pi90UkProfPred <- tmp$pi90UMap
  
  zVal4Plot_PLOT <- zVal4Plot
  
  zkVal4Plot_PLOT <- zkVal4Plot
  pi90LkVal4Plot_PLOT <- zkVal4Plot - 1.64 * sqrt(vkVal4Plot)
  pi90UkVal4Plot_PLOT <- zkVal4Plot + 1.64 * sqrt(vkVal4Plot)
  
  zkProfPred_PLOT <- zkProfPred
  pi90LkProfPred_PLOT <- pi90LkProfPred
  pi90UkProfPred_PLOT <- pi90UkProfPred
  
  xlab <- "Response Variable"
  vecTmp <- c(zVal4Plot_PLOT , zkVal4Plot_PLOT , as.numeric(zkProfPred_PLOT))
  xlim <- c(min(vecTmp) , max(vecTmp))
  
  namePlot = paste0(getwd() , '/plotVal4Plot.pdf')
  
  tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal4Plot , dIData = dIVal4Plot , zData = zVal4Plot_PLOT , 
                           xPred = cVal4PlotU , dIPred = dIPred , zPred = zkProfPred_PLOT , pi90LPred = pi90LkProfPred_PLOT , pi90UPred = pi90UkProfPred_PLOT , 
                           zhatxv = zkVal4Plot_PLOT , pi90Lxv = pi90LkVal4Plot_PLOT , pi90Uxv = pi90UkVal4Plot_PLOT , 
                           profNames = paste0('Profile ' , rand6ForPlot) , xlim = xlim , xlab = xlab) 
  
}

# to integrate
recombine_data <- function(fit_data = fit_data, validate = TRUE, validate_data = NULL) {
  #check if validate is true that validate_data is present
  is_validate_data_missing <- is.data.frame(validate_data) && nrow(validate_data)==0
  
  if (is_validate_data_missing == TRUE & validate == TRUE) {
    stop("expect validation data passed if validate set to TRUE")
  }
  #combining fitted data
  labels_non_covariates <- c("x","y","lowerDI","upperDI","z", "profIDFit")
  core_data_fit <- fit_data[, which(names(fit_data) %in% labels_non_covariates)]
  covariates <- fit_data[ , -which(names(fit_data) %in% labels_non_covariates)]
  new_fit <- list()
  new_fit$cFit <- array(c(core_data_fit$x, core_data_fit$y),dim = c(nrow(core_data_fit),2))
  new_fit$dIFit <- array(c(core_data_fit$lowerDI, core_data_fit$upperDI),dim = c(nrow(core_data_fit),2))
  new_fit$covsFit <- covariates
  new_fit$zFit <- core_data_fit$z
  new_fit$profIDFit <- core_data_fit$profIDFit
  
  total <- c()
  if (validate == TRUE){
    # combine validation data in a like way and join data with fit
    core_data_val <- validate_data[, which(names(validate_data) %in% labels_non_covariates)]
    covariates <- validate_data[ , -which(names(validate_data) %in% labels_non_covariates)]
    new_val = list()
    new_val$cVal <- array(c(core_data_val$x, core_data_val$y),dim = c(nrow(core_data_val),2))
    new_val$dIVal <- array(c(core_data_val$lowerDI, core_data_val$upperDI),dim = c(nrow(core_data_val),2))
    new_val$covsVal <- covariates
    new_val$zVal <- core_data_val$z
    new_val$profIDVal <- core_data_val$profIDFit
    total <- list()
    total <- c(new_fit,new_val)
    
  } else {
    
    total <- new_fit
  }
  
  return(total)
}

reduce_data_based_on_covariate_selection <- function(data,covariate_subset) {
  required <- c('x','y','lowerDI','upperDI','z','profIDFit')
  select <- c(required,covariate_subset)
  for (col in select) {
    if(!(col %in% colnames(data)))
      {
        stop(paste0("missing column expected ",col))
      } 
  }
  data <- data[select]
  return(data)
}
#' This function builds spline model from a uniform data structure available in package (Uniform_Data_Edgeroi).
#' @param fit_data Fit data to feed in such as EdgeroiFitData. Expected to be a dataframe.
#' @param validate_data OPTIONAL Validation data to feed in such as EdgeroiValidationData. Expected to be a dataframe.
#' @param spatialCovs list of spatial covariates with 'dIMidPts' required. example c(''dIMidPts',elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4').
#' @return an object with lmm.fit.selected,xkVal and vkVal
#' @export
#' @examples
#' Splinedata <- SplineIAK(fit_data = EdgeroiFitData,validate_data = EdgeroiValidationData, spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
SplineIAK <- function(fit_data = fit_data,validate_data = validate_data, spatialCovs = spatialCovs) {
  
  fit_data <- reduce_data_based_on_covariate_selection(fit_data,spatialCovs)
  if (exists('validate_data') && is.data.frame(get('validate_data'))) {
    validation <- TRUE
    validate_data <- reduce_data_based_on_covariate_selection(validate_data,spatialCovs)
    } else {
      validation <- FALSE
      validate_data <- NULL
    }
  
  data <- recombine_data(fit_data = fit_data, validate = validation, validate_data = validate_data)
  return(RunEdgeroi(fitCubistModelNow = FALSE,LoadModel = FALSE,validation,data,spatialCovs))
}

#' Run Iak3d project with building Cubist model
#'
#' This function builds cubist model
#' @param fit_data Fit data to feed in such as EdgeroiFitData. Expected to be a dataframe.
#' @param validate_data OPTIONAL Validation data to feed in such as EdgeroiValidationData. Expected to be a dataframe.
#' @param spatialCovs list of spatial covariates with 'dIMidPts' required. example c(''dIMidPts',elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
#' @return an object with lmm.fit.selected,xkVal and vkVal
#' @export
#' @examples
#' Cubistdata <- CubistIAK(fit_data = EdgeroiFitData,validate_data = EdgeroiValidationData, spatialCovs = c('dIMidPts','elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
CubistIAK <- function(fit_data = fit_data,validate_data = validate_data, spatialCovs = spatialCovs) {
  #original usage
  #Cubistdata <- CubistIAK(Uniform_Data_Edgeroi,TRUE,c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4'))
  fit_data <- reduce_data_based_on_covariate_selection(fit_data,spatialCovs)
  if (exists('validate_data') && is.data.frame(get('validate_data'))) {
    validation <- TRUE
    validate_data <- reduce_data_based_on_covariate_selection(validate_data,spatialCovs)
    } else {
      validation <- FALSE
      validate_data <- NULL
    }
  
  data <- recombine_data(fit_data = fit_data, validate = validation, validate_data = validate_data)
  return(RunEdgeroi(fitCubistModelNow = TRUE,LoadModel = FALSE,validation, data,spatialCovs))
}


ModelFromFile <- function(datafeedin){
  #expect a cmFit.RData file to load 
  return(RunEdgeroi(fitCubistModelNow = TRUE,LoadModel = TRUE,validation,data))
}
RunEdgeroi <- function(fitCubistModelNow,LoadModel,validation, datafeedin,spatialCovs){
  assign("last.warning", NULL, envir = baseenv())
  ##############################################################
  ### Model paramaters 
  ##############################################################
  #fitCubistModelNow <- FALSE # fit cubist model if TRUE, spline if no LoadModel given
  #LoadModel <- FALSE # expect a cmFit.RData file to load 
  useCubistForTrend <- fitCubistModelNow # an algorithm to select number of rules for cubist model
  fitModelNow <- TRUE #  runs fitIAK3D assume is generally true. 
  #creates object lmm.fit.selected.RData, otherwise lmmFitFile (lmm.fit.selected.RData) is expected and loaded
  #browser()
  #other paramaters
  plotVargiogramFit <- TRUE
  valNow <- validation
  val4PlotNow <- validation
  mapNow <- FALSE 
  printnllTime <<- FALSE
  CRSAusAlbers <- sp::CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  CRSAusAlbersNEW <- sp::CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  CRSLongLat <- sp::CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
  #options for iak
  nRules <- 5 # number of rules for Cubist model. If NA, then a x val routine to select nRules will be used. 
  refineCubistModel <- TRUE # use a stepwise algorithm (refineXIAK3D in cubist2XIAK3D.R) to remove predictors from the fitted Cubist model? # If NA, a x val routine will be used to select T/F
  constrainX4Pred <- FALSE # use the limits of each column of X (non-zeros only) from the fitting data to constrain the columns of X for prediction?
  prodSum <- TRUE # product sum model? If FALSE, product model is used.  
  nud <- 0.5 #NULL if running the modelSelectIAK3D function.
  lnTfmdData <- FALSE # something about log whether data is log transformed or not
  rqrBTfmdPreds <- FALSE # only relevant if data were log-transformed - back transforms via exp(z + 0.5 * v), ie to minimize expected sqd err
  useReml <- TRUE
  testCL <- FALSE # Keep False for now until associated compLikMats.RData sourced and fed as 
  allKnotsd <- c()
  #incdSpline <- FALSE # wasnt passed through but set afterwards depended on the allKnotsd condition of length 
  # load all in R package development used instead of sourcing various files
  
  otherparamaters <- list(plotVargiogramFit=plotVargiogramFit,valNow=valNow,val4PlotNow=val4PlotNow,mapNow=mapNow,printnllTime=printnllTime,
                          CRSAusAlbers=CRSAusAlbers,CRSAusAlbersNEW=CRSAusAlbersNEW,CRSLongLat=CRSLongLat, nRules= nRules,refineCubistModel=refineCubistModel,
                          constrainX4Pred=constrainX4Pred,prodSum=prodSum,nud=nud,lnTfmdData=lnTfmdData,rqrBTfmdPreds=rqrBTfmdPreds,useReml=useReml,
                          testCL=testCL,allKnotsd=allKnotsd)
  
  #Create main parameter list
  paramaters <- list(fitCubistModelNow=fitCubistModelNow,useCubistForTrend=useCubistForTrend,fitModelNow=fitModelNow, otherparamaters=otherparamaters)
  ModelOutput <- LoadData(paramaters,datafeedin,spatialCovs)
  #iftestCL logic was here::here
  wDir <- here::here()
  lmm2Dir <- here::here('R/fLMM2')
  dataDir <- here::here('tests/run_results') # change this to /src or /data #https://r-pkgs.org/package-structure-state.html
  #dataDir <- here::here('data')
  # when incorporating R package structure
  setwd(wDir)
  
  compLikMats <- list()
  compLikMats$compLikOptn <- 0
  
  lmmFitFile <- paste0(dataDir , '/lmm.fit.selected.RData') # for the fitted model
  nmplt <- paste0(dataDir , '/plot.selected.gam2.pdf') # for a plot with the internal 'predictions' = predictions through profiles of sampled profiles (not validation, can be a check of what's going on)
  
  if(fitModelNow){
    print("fitModelNow now run .........")
    ### refit cubist model as lmm...if selectCovIAK3D was run don't need to do this bit
    start_time <- Sys.time()
    
    if(paramaters$fitCubistModelNow) { # cubist model
      #sdfdTypeANDcmeInit <- c(-9 , -1 , -1 , 1)
      sdfdTypeANDcmeInit <- c(0 , -1 , -1 , 1)
      
      sdfdKnots <- setKnots4sdfd(ModelOutput$dIFit , sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3])
    } else { #spline model
      sdfdTypeANDcmeInit <- c(-9 , -1 , -1 , 1)  # Not sure on logic when this is supposed to be used.
      sdfdKnots <- setKnots4sdfd(ModelOutput$dIFit , sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3])
    }
    # True XData same here
    #Feed model output to main function...
    tmp <- fitIAK3D(xData = ModelOutput$cFit , dIData = ModelOutput$dIFit , zData = ModelOutput$zFit , covsData = ModelOutput$covsFit , 
                    modelX = ModelOutput$modelX , modelx = 'matern' , nud = paramaters$otherparamaters$nud , 
                    allKnotsd = paramaters$otherparamaters$allKnotsd , 
                    sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3] , 
                    cmeOpt = sdfdTypeANDcmeInit[4] , sdfdKnots = sdfdKnots , prodSum = paramaters$otherparamaters$prodSum , lnTfmdData = lnTfmdData , useReml = useReml , compLikMats = compLikMats ,
                    namePlot = nmplt , rqrBTfmdPreds = rqrBTfmdPreds) 
    
    end_time <- Sys.time()
    print('Time to fit was:')
    print(end_time - start_time)
    
    lmm.fit.selected <- tmp$lmmFit
    #save(lmm.fit.selected , file = lmmFitFile)
    
  }else{
    load(file = lmmFitFile)
  }
  
  ###########################################################################
  ### some plots of the fitted covariance model...
  # saves pdfs to getwd()
  ###########################################################################
  # Run plots are exported function requiring result of this run
  #if(plotVargiogramFit){
  #  RunPlots(lmm.fit.selected)
  #}else{}
  if (valNow == TRUE) {
    #VAIDATIONS BIT
    fnamezkVal <- paste0(getwd() , '/zkVal.RData')
    fnamevkVal <- paste0(getwd() , '/vkVal.RData')
    namePlot = paste0(getwd(), '/plotVal.pdf')
    xkvkVal <- RunValidation(ModelOutput,dataDir,namePlot,lmm.fit.selected,rqrBTfmdPreds,constrainX4Pred,fnamezkVal,fnamevkVal)
   
    #Move this to plotting
    #LastSeperation(ModelOutput,lmm.fit.selected ,rand6ForPlot)
  } else {
    xkvkVal <- NULL
  }

  return(list(lmm.fit.selected=lmm.fit.selected,xkvkVal=xkvkVal,InputParamatersList=ModelOutput))
}


PostFitValidate <- function() {
  # this is what was after stop - when the fitting and validation completed
  ########################################################################
  ### set up covaraites for mapping, based on the grid that all covariates are on in rList...
  ########################################################################
  if(mapNow){
    
    #  dIMap <- data.frame('dU' = c(0 , 0.05 , 0.15 , 0.3 , 0.6 , 1.0) , 'dL' = c(0.05 , 0.15 , 0.3 , 0.6 , 1.0 , 2.0))
    dIMap <- data.frame('dU' = c(0 , 0.15 , 0.6) , 'dL' = c(0.05 , 0.3 , 1.0))
    
    firstRow <- 1
    lastRow <- 10
    
    ### get which rows are needed in this batch...
    rowsToDo <- firstRow:lastRow
    
    for(irow in rowsToDo){
      
      ### cell centres for this row...
      xVecMap <- seq(xFromCol(rList[[1]] , 1) , xFromCol(rList[[1]] , ncol(rList[[1]])) , res(rList[[1]])[1])
      yVecMap <- yFromRow(rList[[1]] , irow)
      
      ### get sp::coordinates and covariates for this row...
      cMap <- data.frame('Eastings' = xVecMap , 'Northings' = yVecMap)
      
      ### define cMap and raster::extract covsMap for this row...
      covsMap <- data.frame(matrix(NA , ncol(rList[[1]]) , length(rList)))
      for (icov in 1:length(rList)){
        covsMap[,icov] <- raster::extract(rList[[icov]] , cMap)
      }
      names(covsMap) <- names(rList)
      iIn <- Matrix::which(!is.na(raster::rowSums((covsMap))))
      
      covsMap[['dIMidPts']] <- NA
      
      ################################################
      ### calculate the predictions for the mapping depth...
      ################################################
      lmm.map.tmp <- predictIAK3D(xMap = cMap[iIn,,drop=FALSE] , covsMap = covsMap[iIn,,drop=FALSE] , dIMap = dIMap , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
      lmm.map <- list()
      for (j in 1:length(lmm.map.tmp)){
        lmm.map[[names(lmm.map.tmp)[j]]] <- matrix(NA , nrow(dIMap) , ncol(rList[[1]]))
        lmm.map[[names(lmm.map.tmp)[j]]][,iIn] <- lmm.map.tmp[[names(lmm.map.tmp)[j]]]
      }
      
      #save(lmm.map , file = paste0(dataDir , '/map.row' , irow , '.RData'))
      
    }
    
    print(paste0('Mapped for row ' , irow))
    
  } ### done
}

