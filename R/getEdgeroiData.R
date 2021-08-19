getEdgeroiData <- function(edgeroi,elevation , twi , radK , landsat_b3 , landsat_b4){

  #data(edgeroi)
  edgeroi$sites[edgeroi$sites$SOURCEID=="399_EDGEROI_ed095_1",]
  edgeroi$horizons[edgeroi$horizons$SOURCEID=="399_EDGEROI_ed095_1",]
## spPoints:
  sites <- edgeroi$sites
  sp::coordinates(sites) <- ~ LONGDA94 + LATGDA94
  sp::proj4string(sites) <- sp::CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
  sites <- sp::spTransform(sites, sp::CRS("+init=epsg:28355"))

## plot points and grids:
  pnts <- list("sp.points", sites, pch="+", col="black")

  edgeroi$horizons <- edgeroi$horizons[Matrix::which(!is.na(edgeroi$horizons$CLYPPT)),]
  
  idTmp <- NA * numeric(nrow(edgeroi$horizons))
  dITmp <- cbind(edgeroi$horizons$UHDICM , edgeroi$horizons$LHDICM) / 100
  cTmp <- matrix(NA , nrow(edgeroi$horizons) , 2)
  zTmp <- edgeroi$horizons$CLYPPT

  for(i in 1:nrow(edgeroi$horizons)){
    iThis <- Matrix::which(sites$SOURCEID == edgeroi$horizons$SOURCEID[i])  
    if(length(iThis) == 1){
      idTmp[i] <- iThis

      cTmp[i,] <- sp::coordinates(sites)[iThis,]

      # 899_Forrest_48_1 and 899_Forrest_49_1 have exatly the same coords - shift 899_Forrest_49_1 east by 1m...
      if(as.character(edgeroi$horizons$SOURCEID[i]) == '899_Forrest_49_1'){
        cTmp[i,1] <- cTmp[i,1] + 1
      }else{}

    }else if(length(iThis) == 0){
      print('Site not found for:')
      print(edgeroi$horizons[i,])
    }else{
      print('Multiple sites found for')
      stop(edgeroi$horizons[i,])
    }
  }

# ### find any duplicates...
#   iDup <- Matrix::which(duplicated(cbind(cTmp , dITmp[,1])))
# 
#   {
#     for(i in iDup){
#       print('Duplicated data locations found:')
#       iop = Matrix::which(cTmp[,1] == cTmp[i,1] & cTmp[,2] == cTmp[i,2])
#       print(cbind(cTmp , dITmp , zTmp)[iop,])
#     }
#   }

##############################################################  
### get the ithir::edgeroiCovariates...
##############################################################  
  #data(edgeroiCovariates)
  rList <- list(elevation , twi , radK , landsat_b3 , landsat_b4)
  names(rList) <- c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
  covsTmp <- data.frame('elevation' = NA * numeric(nrow(edgeroi$horizons)) , 'twi' = NA , 'radK' = NA , 'landsat_b3' = NA , 'landsat_b4' = NA)
  
  # elevation : numeric; topographic variable of bare earth ground elevation. Derived from digital elevation model
  # twi : numeric; topographic wetness index. Secondary derivative of the digital elevation model
  # radK : numeric; gamma radiometric data
  # landsat_b3 : numeric; band 3 reflectance of the Landsat 7 satelite 
  # landsat_b4 : numeric; band 4 reflectance of the Landsat 7 satelite
  
  for (j in 1:ncol(covsTmp)){
    covsTmp[,j] <- raster::extract(rList[[j]] , cTmp)
  }

### sort by id then dU...
  iOrder <- order(idTmp , dITmp[,1])
  idTmp <- idTmp[iOrder]
  cTmp <- cTmp[iOrder,]
  dITmp <- dITmp[iOrder,]
  zTmp <- zTmp[iOrder]
  covsTmp <- covsTmp[iOrder,]

### put sp::coordinates into km
  cTmp <- cTmp / 1000

########################################################
### split into Fit and Val (to be predicted) data...
########################################################
  set.seed(123)
  nProfVal <- 60
  uidTmp <- unique(idTmp)
  uidVal <- sample(uidTmp , nProfVal)
  uidFit <- setdiff(uidTmp , uidVal)
  
  iVal <- Matrix::which(is.element(idTmp , uidVal))
  iVal <- iVal[order(iVal)]
  iFit <- setdiff(seq(nrow(covsTmp)) , iVal)
  iFit <- iFit[order(iFit)]
  
### and split c/dI/covs/z into Fit/Val...  
  cFit <- cTmp[iFit,,drop = FALSE] 
  dIFit <- dITmp[iFit,,drop = FALSE] 
  covsFit <- covsTmp[iFit,,drop = FALSE] 
  zFit <- zTmp[iFit]

  cVal <- cTmp[iVal,,drop = FALSE] 
  dIVal <- dITmp[iVal,,drop = FALSE] 
  covsVal <- covsTmp[iVal,,drop = FALSE] 
  zVal <- zTmp[iVal]

##################################################################  
### give Fit data and Val data profIDs, each starting from 1 ...
##################################################################  
  cFitU <- cFit[Matrix::which(!duplicated(cFit)),,drop=FALSE]
  profIDFit <- NA * numeric(nrow(cFit))
  for (i in 1:nrow(cFitU)){
    iThis <- Matrix::which(cFit[,1] == cFitU[i,1] & cFit[,2] == cFitU[i,2])
    profIDFit[iThis] <- i
  }
  
  cValU <- cVal[Matrix::which(!duplicated(cVal)),,drop=FALSE]
  profIDVal <- NA * numeric(nrow(cVal))
  for (i in 1:nrow(cValU)){
    iThis <- Matrix::which(cVal[,1] == cValU[i,1] & cVal[,2] == cValU[i,2])
    profIDVal[iThis] <- i
  }
  
##################################################################  
### for fitting the initial Cubist model, add column dIMidPts to both Fit and Val data
### note that proper support will be used when Cubist regression parameters refitted with iak
### and when predicting (obtained from dIFit and dIVal)
##################################################################  
  covsFit$dIMidPts <- rowMeans(dIFit)
  covsVal$dIMidPts <- rowMeans(dIVal)
  
  return(list('cFit' = cFit , 'dIFit' = dIFit , 'covsFit' = covsFit , 'zFit' = zFit , 'profIDFit' = profIDFit ,
              'cVal' = cVal , 'dIVal' = dIVal , 'covsVal' = covsVal , 'zVal' = zVal , 'profIDVal' = profIDVal , 'rList' = rList))
  
}


