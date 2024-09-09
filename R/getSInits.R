#' Generation of activity centers for SCR and OPSCR initial values.
#'
#' R utility function to generate initial AC locations using compiled Nimble functions for efficiency.
#'
#' The \code{getSInits} function is used in advance of model building.
#'
#' @param AllDetections a 2- or 3-dimensional array of individual detections. 
#' This argument can be provided in two formats: 
#' (i) with the \emph{y} object as returned by \code{\link{getSparseY}} 
#' (ii) with the \emph{yCombined} object as returned by \code{\link{getSparseY}}.
#' The \emph{yCombined} object combines \emph{detNums}, \emph{y}, and \emph{detIndices} (in that order). When such consolidated 
#' representation of the detection data x is used, \emph{detIndices} and \emph{detNums} arguments should not be specified.
#' @param Id.vector A matrix giving the x- and y-coordinates of each trap (scaled to the habitat; see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param idAugmented A matrix giving the x- and y-coordinates of each trap (scaled to the habitat; see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' One row for each window. Each window should be of size 1x1.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in \code{lowerCoords} and \code{upperCoords}. 
#' When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.
#' @param sd Standard deviation of the isotropic bivariate normal distribution.
#' @param movementMethod a character denoting the nimbleSCR to be used.
#' 
#' @return This function returns a 2- or 3-dimensional array of initial AC locations.
#' 
#' @author Pierre Dupont
#'
#' @import nimble 
#' @importFrom nimbleSCR rbernppACmovement_normal rbernppAC
#' 
#' @rdname getSInits
#' @export
getSInits <- function( AllDetections = AllDetections,
                       Id.vector = Id.vector, # vector with the id of the individuals in the same order than in Y
                       idAugmented = idAugmented,#Id in the Id.vector that are augemented individuals 
                       lowerCoords = lowerCoords,
                       upperCoords = upperCoords,
                       habitatGrid = habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")
  {
  #compile R nimble function for faster runs 
  if(movementMethod == "dbernppACmovement_normal"){
  crbernppACmovement_normal <- compileNimble(rbernppACmovement_normal)
  }
  crbernppAC <- compileNimble(rbernppAC)
  
  #GET NUMBER OF YEARS 
  AllDetections$Id <- as.character(AllDetections$Id) 
  years <- sort(unique(AllDetections[,"Year"]))
  n.years <- length(years)
 
  
  detected <- table(AllDetections[,"Id"],AllDetections[,"Year"])
  #reorder the table to follow the vector
  Id.vectorOriginal <- Id.vector
  idAugmentedOriginal <- idAugmented
  if(sum(!Id.vector %in% row.names(detected))>0){# if some ids are present in the y but not in the detections, then assume they are augmented ids
    whichNoDetections <- Id.vector %in% row.names(detected)
    idAugmented <- c(idAugmented, which(!whichNoDetections))
    Id.vector <- Id.vector[whichNoDetections]
    
  }else{
  detected <- detected[Id.vector,]
  }
  
  n.detected <- length(Id.vector)
  n.augmentedId <- length(idAugmented) 
  n.individuals <- n.detected + n.augmentedId
  
  #INITIALIZE ACS LOCATION 
  s <- array(NA, c(n.individuals, 2, n.years))
  dimnames(s)[[2]] <- c("x","y")
  dimnames(s)[[1]] <- c(Id.vectorOriginal, paste("Augmented",1:length(idAugmentedOriginal)))
  aug.vect <- dimnames(s)[[1]][idAugmented]
  ##Get object ready for pp functions 
  numGridRows <- dim(habitatGrid)[1]
  numGridCols <- dim(habitatGrid)[2]
  numHabWindows <- nrow(lowerCoords)
  
  #Intensity vector/matrix
  if(is.null(intensity)){
    intensity <- matrix(2,nrow=numHabWindows,ncol=n.years)
      for(t in 1:n.years){
      intensity[,t] <- rep(2,numHabWindows)
      }
  }
  if(is.null(dim(intensity))){
    intensity1 <- matrix(2,nrow=numHabWindows,ncol=n.years)
      for(t in 1:n.years){
        intensity1[,t] <- intensity
      }
        intensity <- intensity1
  }
  
  logIntensities <- log(intensity)
  logSumIntensity <- apply(logIntensities, 2, sum)
  
  
  ##ID DETECTED, AVERAGE LOCATION OF DETECTIONS. 
  for(i in 1:n.detected){
    whereDetections <- which(detected[i,]>0)
    detNums <- detected[Id.vector[i],whereDetections]
    for(t in whereDetections){
        tmp <- AllDetections[AllDetections$Id %in% Id.vector[i] &
                         AllDetections$Year %in% years[t],]
        s[Id.vector[i],,t] <- colMeans(tmp[,c("x","y")])
        
        #if out of habitat, find the closest habitat cell
        sxyID <- habitatGrid[trunc(s[Id.vector[i],2,t])+1, trunc(s[Id.vector[i],1,t])+1]
        if(sxyID==0){
          which.minhab <- which.min(abs(s[Id.vector[i],c("x"),t]-lowerCoords[,c("x")]) + 
                                      abs(s[Id.vector[i],c("y"),t]-lowerCoords[,c("y")]))
          s[Id.vector[i],,t] <- lowerCoords[which.minhab,]
        }
    }
    yrs <- 1:n.years
    whichYearNotDet <- yrs[which(!(detected[i,]>0))]
    whichYearDet <- yrs[which((detected[i,]>0))]
    minYearDet <- min(whichYearDet)
    # if the first year is detected
    
    if(sum(whichYearDet%in%1)==1){
      for(t in whichYearNotDet){
      s[Id.vector[i],,t] <- crbernppACmovement_normal(
        n = 1,
        lowerCoords = lowerCoords,
        upperCoords = upperCoords, 
        s = s[Id.vector[i],1:2,t-1],
        sd = 4,
        baseIntensities = intensity[,t],
        habitatGrid =  habitatGrid,
        numGridRows = numGridRows,
        numGridCols = numGridCols,
        numWindows= numHabWindows)
      }#t
    } else {
      #backward
      for(t in minYearDet:2){
      s[Id.vector[i],,t-1] <- crbernppACmovement_normal(
        n = 1,
        lowerCoords = lowerCoords,
        upperCoords = upperCoords,
        s = s[Id.vector[i], 1:2, t],
        sd = 4,
        baseIntensities = intensity[,t],
        habitatGrid =  habitatGrid,
        numGridRows = numGridRows,
        numGridCols = numGridCols,
        numWindows = numHabWindows)
      }#t
     # Forward
      #remove the first year
      whichYearNotDet <- whichYearNotDet[-1]
      
      for(t in whichYearNotDet){
        s[Id.vector[i],,t] <- crbernppACmovement_normal(
          n = 1,
          lowerCoords = lowerCoords,
          upperCoords = upperCoords,
          s = s[Id.vector[i], 1:2, t-1],
          sd = sd,
          baseIntensities = intensity[,t],
          habitatGrid =  habitatGrid,
          numGridRows = numGridRows,
          numGridCols = numGridCols,
          numWindows = numHabWindows)
      }#t
    }
  }
    
#Augmented ids 
  for(i in 1:n.augmentedId){
    s[aug.vect[i],,1] <- crbernppAC(
      n = 1,
      lowerCoords = lowerCoords,
      upperCoords = upperCoords,
      logIntensities = logIntensities[,1],
      logSumIntensity = logSumIntensity[1],
      habitatGrid = habitatGrid,
      numGridRows = numGridRows,
      numGridCols = numGridCols)
  #t>1
    for(t in 2:n.years){
      s[aug.vect[i],,t] <- crbernppACmovement_normal(
        n = 1,
        lowerCoords = lowerCoords,
        upperCoords = upperCoords,
        s = s[aug.vect[i], 1:2, t - 1],
        sd = sd,
        baseIntensities = intensity[,t],
        habitatGrid =  habitatGrid,
        numGridRows = numGridRows,
        numGridCols = numGridCols,
        numWindows= numHabWindows)
    }#t
}#i
  
  return(s)
  
}
