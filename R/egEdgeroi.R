##############################################################
#
# Underlying data taken from data(edgeroi) and data(edgeroiCovariates) packages
# in GSIF
# datafeed <- list(edgeroi=edgeroi,elevation=elevation,landsat_b3=landsat_b3,landsat_b4=landsat_b4,radK=radK,twi=twi)
#Usage: 
# result <- Iak3dSIH::CubistIAK(datafeedin = Uniform_Data_Edgeroi) OR
# result <- Iak3dSIH::SplineIAK(datafeedin = Uniform_Data_Edgeroi)

##############################################################


#' Represents the underlying data needed to run Spline or Cubist Models within this package.
#' Underlying data taken from data(edgeroi) and data(edgeroiCovariates) packages
#' in GSIF. Functions that build on this data (i.e. getEdgeroiData() or getOodnadattaData()
#' have been applied to give this data i.e. Uniform_Data_Edgeroi. 
#' You can supply other data objects within the model (i.e.Uniform_Data_Ood ) taken 
#' from other regions, provided it fits similar structures.
#' @name Uniform_Data_Edgeroi
#' @docType data
#' @keywords data
NULL



FitSplineModel <- function(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) {
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
  spatialCovs <- c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
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

FitCubistModel <- function(paramaters, tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) {
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
LoadModel <- function(paramaters) {
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
    ModelOutput <- FitCubistModel(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) #................
  }
  else if (paramaters$LoadModel) {
    print("loading directly from a file")
    #load directly from a file that is expected 
    ModelOutput <- LoadModelDirectly(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList)
  }
  else {
    print("doing spline")
    #spline model
    ModelOutput <- FitSplineModel(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) # first seperation into Spline and return stuff
  }
  return(ModelOutput)
}

LoadData <- function(paramaters,datafeedin){
  
  ##############################################
  ### load the edgeroi dataset (from GSIF package) and put into format for iak3d...
  ##############################################
  print("now in loadData................................")
  
  tmp <- datafeedin
  #tmp <- getEdgeroiData(datafeed$edgeroi, datafeed$elevation , datafeed$twi , datafeed$radK , datafeed$landsat_b3 , datafeed$landsat_b4)
  #saveRDS(tmp,here::here("Original_Edgeroi.rds"))
  #determin model options from flags and be able to pass all other paramaters needed
  print("Print constructed ModelOptions")
  ModelOptions <- list(FitCubits=paramaters$fitCubistModelNow, 
                       useCubistForTrend = paramaters$useCubistForTrend, 
                       LoadModel = !paramaters$fitModelNow, 
                       otherparamaters=paramaters$otherparamaters,
                       data=tmp)
  output <- LoadModel(ModelOptions)
  
  return(output)
  
}

#' Run Iak3d project with building spline model
#'
#' RunPlots(iakdata$lmm.fit.selected)
#' @param lmm.fit.selected The result of running either Cubist or Spline Models
#' @return saves plots in the working directory
#' @export
RunPlots <- function(lmm.fit.selected) {
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

LastSeperation <- function(ModelOutput,dataDir,lmm.fit.selected , rqrBTfmdPreds , constrainX4Pred){
  
  
  rand6ForPlot <- c(6 , 19 , 49 , 41 , 3 , 24) # For edg
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
  
  xlab <- "Clay percent"
  vecTmp <- c(zVal4Plot_PLOT , zkVal4Plot_PLOT , as.numeric(zkProfPred_PLOT))
  xlim <- c(min(vecTmp) , max(vecTmp))
  
  namePlot = paste0(getwd() , '/plotVal4Plot.pdf')
  
  # turned off for now.  returns null
  
  tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal4Plot , dIData = dIVal4Plot , zData = zVal4Plot_PLOT , 
                           xPred = cVal4PlotU , dIPred = dIPred , zPred = zkProfPred_PLOT , pi90LPred = pi90LkProfPred_PLOT , pi90UPred = pi90UkProfPred_PLOT , 
                           zhatxv = zkVal4Plot_PLOT , pi90Lxv = pi90LkVal4Plot_PLOT , pi90Uxv = pi90UkVal4Plot_PLOT , 
                           profNames = paste0('Profile ' , rand6ForPlot) , xlim = xlim , xlab = xlab) 
  
}


#' This function builds spline model from a uniform data structure available in package (Uniform_Data_Edgeroi).
#' @param data data to feed in such as Uniform_Data_Edgeroi
#' @param vaidation Boolean TRUE or FALSE
#' @return an object with lmm.fit.selected,xkVal and vkVal
#' @export
#' @examples
#' Splinedata <- SplineIAK(Uniform_Data_Edgeroi,TRUE)
SplineIAK <- function(data,validation) {
  return(RunEdgeroi(fitCubistModelNow = FALSE,LoadModel = FALSE,validation,data))
}

#' Run Iak3d project with building Cubist model
#'
#' This function builds cubist model
#' @param data data to feed in such as Uniform_Data_Edgeroi
#' @param vaidation Boolean TRUE or FALSE
#' @return an object with lmm.fit.selected,xkVal and vkVal
#' @export
#' @examples
#' Cubistdata <- CubistIAK(Uniform_Data_Edgeroi,TRUE)
CubistIAK <- function(data,validation) {
  return(RunEdgeroi(fitCubistModelNow = TRUE,LoadModel = FALSE,validation, data))
}


ModelFromFile <- function(datafeedin){
  #expect a cmFit.RData file to load 
  return(RunEdgeroi(fitCubistModelNow = TRUE,LoadModel = TRUE,validation,data))
}
RunEdgeroi <- function(fitCubistModelNow,LoadModel,validation, datafeedin){
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
  ModelOutput <- LoadData(paramaters,datafeedin)
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
    browser()
    LastSeperation(ModelOutput,dataDir,lmm.fit.selected , rqrBTfmdPreds , constrainX4Pred)
  } else {
    xkvkVal <- NULL
  }

  return(list(lmm.fit.selected=lmm.fit.selected,xkvkVal=xkvkVal))
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

