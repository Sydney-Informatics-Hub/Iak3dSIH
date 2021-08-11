getOodnadattaData <- function(path = "./data/Oodnadatta_soil_covariates.csv"){

 # want: 
  # tmp <- getOodnadattaData() # list of data, function organises, gets rid of duplicates, appends covariates
  # cFit <- tmp$cFit    # sp::coordinates (in Km) of calib set
  # dIFit <- tmp$dIFit      # Depth intervals of fit set (metres)
  # covsFit <- tmp$covsFit    # covariates of fit set
  # zFit <- tmp$zFit          # response variable of fit set
  # profIDFit <- tmp$profIDFit # profile ID of fit set
  # 
  # cVal <- tmp$cVal       # sp::coordinates (in Km) of valid set
  # dIVal <- tmp$dIVal      # depth intervals of valid set (in metres)
  # covsVal <- tmp$covsVal  # Covariates for valid set
  # zVal <- tmp$zVal         # response variable for valid set
  # profIDVal <- tmp$profIDVal  # profile ID
  # 
  # rList <- tmp$rList  # raster list of covariates
  #
  
  
  ########################################################
  ### Read in the data...
  ########################################################
  
  soil = read.csv(path)
  
  soil$Site_ID = paste(soil$Field, soil$Sample.ID, sep="_")
  
  ##############################################################  
  ### Convert sp::coordinates to km...
  ##############################################################  
  
  # So far only need to transform into km
  soil$x_km <- soil$x / 1000
  soil$y_km <- soil$y / 1000  
  
  ##############################################################  
  ### Convert depth intervals to m...
  ##############################################################  
  
  soil$depth_upper_m = soil$DepthMin.cm./100
  soil$depth_lower_m = soil$DepthMax.cm./100
  # soil$dIMidPts = soil$mid_depth/100   
  depths = cbind(soil$depth_upper_m, soil$depth_lower_m)
  soil$dIMidPts <- rowMeans(depths)
  
  
  # sort by ID, then depth upper
  iOrder = order(soil$Site_ID , soil$depth_lower_m)
  soil = soil[iOrder,]
  
  ########################################################
  ### split into Fit and Val (to be predicted) data...
  ########################################################
  
  set.seed(123)
  nProfVal <- round(nrow(soil[unique(soil$Site_ID),])*0.25, 0)    #25/75 split for now 
  ids <- unique(soil$Site_ID)
  uidVal <- sample(ids, nProfVal)
  pos = Matrix::which(soil$Site_ID%in%uidVal)
  calib = soil[-pos,]
  valid = soil[pos,]
  
  # Prep for export, make similar to Edgeroi output
  
  cFit = cbind(calib$x_km,calib$y_km)  
  dIFit = cbind(calib$depth_upper_m,calib$depth_lower_m) 
  covsFit = as.data.frame(cbind(calib$DEM_30, calib$NDVI_5,
                                calib$NDVI_50, calib$NDVI_95,
                                 calib$gamma_k,calib$EM_100, calib$MrRTF, 
                                calib$rad_dose, calib$dIMidPts))
  colnames(covsFit) = c("DEM_30", "NDVI_5", "NDVI_50", "NDVI_95", 
                        "gamma_k", "EM_100", "MrRTF", 
                        "rad_dose", "dIMidPts")
  zFit <- calib$ESP         
  profIDFit <- calib$Site_ID 
  
  cVal <- cbind(valid$x_km,valid$y_km)       
  dIVal <- cbind(valid$depth_upper_m,valid$depth_lower_m)      
  covsVal <- as.data.frame(cbind(valid$DEM_30, valid$NDVI_5,
                                 valid$NDVI_50,valid$NDVI_95,
                                 valid$gamma_k, valid$EM_100, valid$MrRTF, 
                                 valid$rad_dose, valid$dIMidPts))  
  colnames(covsVal) = c("DEM_30","NDVI_5", "NDVI_50", "NDVI_95", 
                        "gamma_k", "EM_100", "MrRTF", 
                        "rad_dose", "dIMidPts")
  zVal <- valid$ESP         
  profIDVal <- valid$Site_ID  
  
  
##############################################################  
### get the Covariates...
##############################################################  
#   data(ithir::edgeroiCovariates)
#   rList <- list(elevation , twi , radK , landsat_b3 , landsat_b4)
#   names(rList) <- c('elevation' , 'twi' , 'radK' , 'landsat_b3' , 'landsat_b4')
#   covsTmp <- data.frame('elevation' = NA * numeric(nrow(edgeroi$horizons)) , 'twi' = NA , 'radK' = NA , 'landsat_b3' = NA , 'landsat_b4' = NA)
#   
#   # elevation : numeric; topographic variable of bare earth ground elevation. Derived from digital elevation model
#   # twi : numeric; topographic wetness index. Secondary derivative of the digital elevation model
#   # radK : numeric; gamma radiometric data
#   # landsat_b3 : numeric; band 3 reflectance of the Landsat 7 satelite 
#   # landsat_b4 : numeric; band 4 reflectance of the Landsat 7 satelite
#   
#   for (j in 1:ncol(covsTmp)){
#     covsTmp[,j] <- raster::extract(rList[[j]] , cTmp)
#   }
# 
# ### sort by id then dU...
#   iOrder <- order(idTmp , dITmp[,1])
#   idTmp <- idTmp[iOrder]
#   cTmp <- cTmp[iOrder,]
#   dITmp <- dITmp[iOrder,]
#   zTmp <- zTmp[iOrder]
#   covsTmp <- covsTmp[iOrder,]
# 
# ### put sp::coordinates into km --> moved this to top as did not need to raster::extract variables
#   cTmp <- cTmp / 1000
# 

 ##################################################################  
 ### give Fit data and Val data profIDs, each starting from 1 ...
 ##################################################################  
   # cFitU <- cFit[Matrix::which(!duplicated(cFit)),,drop=FALSE]
   # profIDFit <- NA * numeric(nrow(cFit))
   # for (i in 1:nrow(cFitU)){
   #   iThis <- Matrix::which(cFit[,1] == cFitU[i,1] & cFit[,2] == cFitU[i,2])
   #   profIDFit[iThis] <- i
   # }
   # 
   # cValU <- cVal[Matrix::which(!duplicated(cVal)),,drop=FALSE]
   # profIDVal <- NA * numeric(nrow(cVal))
   # for (i in 1:nrow(cValU)){
   #   iThis <- Matrix::which(cVal[,1] == cValU[i,1] & cVal[,2] == cValU[i,2])
   #   profIDVal[iThis] <- i
   # }
  
  # Liana comment: Have tried doing this vs provide profile IDs, doesn't change output
  
##################################################################  
### for fitting the initial Cubist model, add column dIMidPts to both Fit and Val data
### note that proper support will be used when Cubist regression parameters refitted with iak
### and when predicting (obtained from dIFit and dIVal)
##################################################################  
  # covsFit$dIMidPts <- rowMeans(dIFit)
  # covsVal$dIMidPts <- rowMeans(dIVal)
  # 

  # raster list= may be needed later for mapping, but null for now
   rList=c()
   
return(list('cFit' = cFit , 'dIFit' = dIFit , 'covsFit' = covsFit , 'zFit' = zFit , 'profIDFit' = profIDFit ,
            'cVal' = cVal , 'dIVal' = dIVal , 'covsVal' = covsVal , 'zVal' = zVal , 'profIDVal' = profIDVal , 'rList' = rList))
  
  #return(list('Fit' = calib, 'Val' = valid))
  
}

