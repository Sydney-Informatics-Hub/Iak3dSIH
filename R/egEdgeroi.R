##############################################################
# libraries required
library(sp)
library(raster)
library(rgdal)
library(Cubist)
library(mgcv)
library(Matrix)
library(MASS)
library(splines)
library(deldir)
library(lme4)
library(aqp)
library(GSIF)
library(ithir)
library(parallel)
library(here)
##############################################################
FitSplineModel <- function(paramaters,tmp,cFit, dIFit, covsFit, zFit, profIDFit, cVal, dIVal, covsVal, zVal, profIDVal, rList) {
  #################################################################################################
  ### set knots for sdfd spline fn (if used)
  #################################################################################################
  ### alternatively, set up for fitting a spline model.
  ### include interactions between depth and spatial covariates (but here not between different spatial covariates)
  ###################################################################################
  print("FitSplineModel is engaged")
  modelX <- list('type' = 'gam2')
  scaleCovs <- TRUE
  nIntKnotsd <- 4 # number of internal knots for the spline function (nat spline, clamped to have grad=0 at upper bdry) of depth; if this is complex enough, probably no need for the depth component in prod-sum covariance model
  nIntKnotss <- 4 # number of internal knots for the spline functions (nat spline, clamped to have grad=0 at upper and lower bdries) of covariates

    ### don't include depth here.   
  spatialCovs <- c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
  if(scaleCovs){
    print("Scaled covariates created")
    ### to work with scaled covariates
    spatialCovs <- paste0(spatialCovs , '_SCALED')
  }else{
    ### to work with unscaled (raw) covariates
    spatialCovs <- spatialCovs # no change here
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
    
    saveRDS(covsFit, file = here("tests/run_results/check_x_cubits.rds"))
    saveRDS(zFit, file = here("tests/run_results/check_y_cubits.rds"))
    saveRDS(paramaters$otherparamaters$nRules, file = here("tests/run_results/check_nRules.rds"))
    
    cmFit <- cubist(x = covsFit , y = zFit , committees = 1 , cubistControl(rules = paramaters$otherparamaters$nRules))
   
    
    saveRDS(cmFit, file = here("tests/run_results/check_cubistModel.rds")) # just make sure
    saveRDS(covsFit, file = here("tests/run_results/check_dataFit.rds")) # just make sure
    saveRDS(zFit, file = here("tests/run_results/check_zFit.rds")) # just make sure
    saveRDS(profIDFit, file = here("tests/run_results/check_profIDFit.rds"))
    saveRDS(paramaters$otherparamaters$allKnotsd, file = here("tests/run_results/check_allKnotsd.rds"))
    saveRDS(paramaters$otherparamaters$refineCubistModel, file = here("tests/run_results/check_refineCubistModel.rds"))
    
    
    ### convert to des mtx
    tmp <- cubist2X(cubistModel = cmFit, dataFit = covsFit , zFit = zFit , profIDFit = profIDFit , allKnotsd = paramaters$otherparamaters$allKnotsd , refineCubistModel = paramaters$otherparamaters$refineCubistModel)
    cmFit <- tmp$cubistModel
    XFit <- tmp$X
    matRulesFit <- tmp$matRuleData
    save(cmFit , file = paste0(dataDir , '/cmFit.RData'))

    saveRDS(cmFit, file = here("tests/run_results/check_cmFit.rds")) #problem
    
  }
    
  modelX <- cmFit
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

LoadData <- function(paramaters){

  ##############################################
  ### load the edgeroi dataset (from GSIF package) and put into format for iak3d...
  ##############################################
  print("now in loadData................................")
  
  tmp <- getEdgeroiData()
  
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

# ModelOutput
#list(modelX=modelX,cFit=cFit, dIFit=dIFit, covsFit=covsFit, zFit=zFit, profIDFit=profIDFit, cVal=cVal, dIVal=dIVal, covsVal=covsVal, zVal=zVal, profIDVal=profIDVal, rList=rList))

RunPlots <- function(ModelOutput,dataDir,lmm.fit.selected) {
  dIPlot <- data.frame('dU' = c(0 , 20 , 50 , 90 , 150 , 190)/100 , 'dL' = c(10 , 30 , 60 , 100 , 160 , 200)/100)
  hx <- seq(0 , 20 , 1)
  pdf(file = paste0(dataDir , '/varioFitgam22.pdf'))
  tmp <- plotCovx(lmm.fit = lmm.fit.selected , hx = hx , dIPlot = dIPlot , addExpmntlV = TRUE , hzntlUnits = 'km')
  dev.off()
  
  hdPlot <- seq(0 , 2 , 0.01)
  pdf(file = paste0(dataDir , '/cordFit.pdf'))
  qwe <- plotCord(lmm.fit = lmm.fit.selected , hdPlot = hdPlot, vrtclUnits = 'm')
  dev.off()
  
  dTmp <- seq(0 , 2 , 0.1)
  dIPlot <- data.frame('dU' = dTmp[-length(dTmp)] , 'dL' = dTmp[-1])
  pdf(file = paste0(dataDir , '/covardFit.pdf'))
  qwe <- plotCovd(lmm.fit = lmm.fit.selected , dIPlot = dIPlot , vrtclUnits = 'm')
  dev.off()
  
  ### plot of the variances...
  pdf(file = paste0(dataDir , '/varComps.pdf'))
  dPlot <- seq(0 , 2 , 0.01)
  plotVarComps(lmm.fit = lmm.fit.selected , dPlot = dPlot)
  dev.off()
  
}
RunEdgeroi <- function(){
  assign("last.warning", NULL, envir = baseenv())
  ##############################################################
  ### Model paramaters 
  ##############################################################
  fitCubistModelNow <- TRUE # fit cubist model, spline if no LoadModel given
  LoadModel <- FALSE # expect a cmFit.RData file to load 
  useCubistForTrend <- TRUE # an algorithm to select number of rules for cubist model
  fitModelNow <- TRUE #  runs fitIAK3D assume is generally true. 
                      #creates object lmm.fit.selected.RData, otherwise lmmFitFile (lmm.fit.selected.RData) is expected and loaded

  #other paramaters
  plotVargiogramFit <- TRUE
  valNow <- TRUE
  val4PlotNow <- TRUE
  mapNow <- FALSE 
  printnllTime <<- FALSE
  crsAusAlbers <- CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  crsAusAlbersNEW <- CRS("+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  crsLongLat <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
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
                    crsAusAlbers=crsAusAlbers,crsAusAlbersNEW=crsAusAlbersNEW,crsLongLat=crsLongLat, nRules= nRules,refineCubistModel=refineCubistModel,
                    constrainX4Pred=constrainX4Pred,prodSum=prodSum,nud=nud,lnTfmdData=lnTfmdData,rqrBTfmdPreds=rqrBTfmdPreds,useReml=useReml,
                    testCL=testCL,allKnotsd=allKnotsd)
  
  #Create main parameter list
  paramaters <<- list(fitCubistModelNow=fitCubistModelNow,useCubistForTrend=useCubistForTrend,fitModelNow=fitModelNow, otherparamaters=otherparamaters)
  ModelOutput <- LoadData(paramaters)
  #iftestCL logic was here
  wDir <- here()
  lmm2Dir <- here('R/fLMM2')
  dataDir <- here('tests/run_results') # change this to /src or /data #https://r-pkgs.org/package-structure-state.html
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
      
      sdfdKnots <- setKnots4sdfd(dIFit , sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3])
    } else { #spline model
      sdfdTypeANDcmeInit <- c(-9 , -1 , -1 , 1)  # Not sure on logic when this is supposed to be used.
      sdfdKnots <- setKnots4sdfd(dIFit , sdfdType_cd1 = sdfdTypeANDcmeInit[1] , sdfdType_cxd0 = sdfdTypeANDcmeInit[2] , sdfdType_cxd1 = sdfdTypeANDcmeInit[3])
    }
    print("check final params before fitIAK3D........")
 
    print("For Testing only - all parameters of fitIAK in order are ..........")
    saveRDS(ModelOutput$cFit, file = here("tests/run_results/check_cFit.rds"))
    saveRDS(ModelOutput$dIFit, file = here("tests/run_results/check_dIFit.rds"))
    saveRDS(ModelOutput$zFit, file = here("tests/run_results/check_zFit.rds"))
    saveRDS(ModelOutput$covsFit, file = here("tests/run_results/check_covsFit.rds"))
    saveRDS(ModelOutput$modelX, file = here("tests/run_results/check_modelX.rds"))
    saveRDS(paramaters$otherparamaters$nud, file = here("tests/run_results/check_nud.rds"))
    saveRDS(paramaters$otherparamaters$allKnotsd, file = here("tests/run_results/check_allKnotsd.rds"))
    saveRDS(sdfdTypeANDcmeInit, file = here("tests/run_results/check_sdfdTypeANDcmeInit.rds"))
    saveRDS(sdfdKnots, file = here("tests/run_results/check_sdfdKnots.rds"))
    saveRDS(paramaters$otherparamaters$prodSum, file = here("tests/run_results/check_prodSum.rds"))
    saveRDS(lnTfmdData, file = here("tests/run_results/check_lnTfmdData.rds"))
    saveRDS(useReml, file = here("tests/run_results/check_useReml.rds"))
    saveRDS(compLikMats, file = here("tests/run_results/check_compLikMats.rds"))
    saveRDS(rqrBTfmdPreds, file = here("tests/run_results/check_rqrBTfmdPreds.rds"))
    saveRDS(nmplt, file = here("tests/run_results/check_nmplt.rds"))
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
    save(lmm.fit.selected , file = lmmFitFile)

  }else{
    load(file = lmmFitFile)
  }

  ###########################################################################
  ### some plots of the fitted covariance model...
  ###########################################################################
  if(plotVargiogramFit){
    RunPlots(ModelOutput,dataDir,lmm.fit.selected)
  }else{}

  #########################################################
  ### validation bit...
  #########################################################
  fnamezkVal <- paste0(dataDir , '/zkVal.RData')
  fnamevkVal <- paste0(dataDir , '/vkVal.RData')
  namePlot = paste0(dataDir , '/plotVal.pdf')

  if(valNow){

    nVal <- nrow(cVal)
    iU <- which(!duplicated(cVal))
    cValU <- cVal[iU,,drop=FALSE]
    covsValU <- covsVal[iU,,drop=FALSE]
    zkVal <- vkVal <- NA * numeric(nVal)
    for(i in 1:nrow(cValU)){
        iTmp <- which(cVal[,1] == cValU[i,1] & cVal[,2] == cValU[i,2])
        tmp <- profilePredictIAK3D(xMap = cValU[i,,drop=FALSE] , covsMap = covsVal[iTmp,,drop=FALSE] , dIMap = dIVal[iTmp,,drop=FALSE] , lmmFit = lmm.fit.selected , rqrBTfmdPreds = rqrBTfmdPreds , constrainX4Pred = constrainX4Pred)
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
    tmp <- calcValStats(zVal = zVal , dIVal = dIVal , zkVal = zkVal , vkVal = vkVal , layerMidPts = c(0.025 , 0.1 , 0.225 , 0.45 , 0.8 , 1.5) , printValStats = TRUE)
    valStatsAllLayers <- tmp$valStatsAllLayers 
    valStatsTot <- tmp$valStatsTot
    
    tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal , dIData = dIVal , zData = zVal , 
                    xPred = cValU , dIPred = dIPred , zPred = zkProfPred , pi90LPred = pi90LkProfPred , pi90UPred = pi90UkProfPred , 
                    zhatxv = zkVal , pi90Lxv = zkVal - 1.64 * sqrt(vkVal) , pi90Uxv = zkVal + 1.64 * sqrt(vkVal)) 

    save(zkVal , file = fnamezkVal)
    save(vkVal , file = fnamevkVal)

  }else{

    load(file = fnamezkVal)
    load(file = fnamevkVal)

  }
  
  ### keeping this bit separate. 
  if(val4PlotNow){

    rand6ForPlot <- c(6 , 19 , 49 , 41 , 3 , 24)

    i4PlotU <- which(!duplicated(cVal))[rand6ForPlot]
    cVal4PlotU <- cVal[i4PlotU,,drop=FALSE]
    covsVal4PlotU <- covsVal[i4PlotU,,drop=FALSE]

    iVal4Plot <- c()
    for(i in 1:nrow(cVal4PlotU)){
        iTmp <- which(cVal[,1] == cVal4PlotU[i,1] & cVal[,2] == cVal4PlotU[i,2])
        iVal4Plot <- c(iVal4Plot , iTmp)
    }

    cVal4Plot <- cVal[iVal4Plot,,drop=FALSE]
    dIVal4Plot <- dIVal[iVal4Plot,,drop=FALSE]
    covsVal4Plot <- covsVal[iVal4Plot,,drop=FALSE]
    zVal4Plot <- zVal[iVal4Plot]

    zkVal4Plot <- vkVal4Plot <- NA * numeric(length(zVal4Plot))
    for(i in 1:nrow(cVal4PlotU)){
        iTmp <- which(cVal4Plot[,1] == cVal4PlotU[i,1] & cVal4Plot[,2] == cVal4PlotU[i,2])
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

    namePlot = paste0(dataDir , '/plotVal4Plot.pdf')

    tmp <- plotProfilesIAK3D(namePlot = namePlot , xData = cVal4Plot , dIData = dIVal4Plot , zData = zVal4Plot_PLOT , 
                    xPred = cVal4PlotU , dIPred = dIPred , zPred = zkProfPred_PLOT , pi90LPred = pi90LkProfPred_PLOT , pi90UPred = pi90UkProfPred_PLOT , 
                    zhatxv = zkVal4Plot_PLOT , pi90Lxv = pi90LkVal4Plot_PLOT , pi90Uxv = pi90UkVal4Plot_PLOT , 
                    profNames = paste0('Profile ' , rand6ForPlot) , xlim = xlim , xlab = xlab) 

  }else{}

  setUptests(lmm.fit.selected,vkVal,zkVal)
  
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
      
      ### get coordinates and covariates for this row...
      cMap <- data.frame('Eastings' = xVecMap , 'Northings' = yVecMap)
      
      ### define cMap and extract covsMap for this row...
      covsMap <- data.frame(matrix(NA , ncol(rList[[1]]) , length(rList)))
      for (icov in 1:length(rList)){
        covsMap[,icov] <- extract(rList[[icov]] , cMap)
      }
      names(covsMap) <- names(rList)
      iIn <- which(!is.na(rowSums((covsMap))))
      
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
      
      save(lmm.map , file = paste0(dataDir , '/map.row' , irow , '.RData'))
      
    }
    
    print(paste0('Mapped for row ' , irow))
    
  } ### done
}

