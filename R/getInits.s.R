#' Activity center generator for SCR initial values.
#'
#' R utility function to generate initial AC locations that are compatible with 
#' a given SCR or OPSCR model.
#'
#' The \code{getInits.s} function is used in advance of model building.
#'
#' @param y a 2- or 3-dimensional array of individual detections. 
#' This argument can be provided in two formats: 
#' (i) with the \emph{y} object as returned by \code{\link{getSparseY}} 
#' (ii) with the \emph{yCombined} object as returned by \code{\link{getSparseY}}.
#' The \emph{yCombined} object combines \emph{detNums}, \emph{y}, and \emph{detIndices} (in that order). When such consolidated 
#' representation of the detection data x is used, \emph{detIndices} and \emph{detNums} arguments should not be specified.
#' @param trapCoords A matrix giving the x- and y-coordinates of each trap (scaled to the habitat; see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' One row for each window. Each window should be of size 1x1.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in \code{lowerCoords} and \code{upperCoords}. 
#' When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.
#' @param known.s coordinates of known AC locations.
#' @param baseIntensities Vector of baseline habitat intensities for all habitat windows.
#' @param sd Standard deviation of the isotropic bivariate normal distribution.
#' @param detNums an \code{integer} denoting the number of detections of the focal individual.
#' @param detIndices a vector of length \emph{detNums} denoting the IDs of the detectors where the focal individual was detected.
#' 
#' @return This function returns a 2- or 3-dimensional array of initial AC locations.
#' 
#' @author Pierre Dupont
#'
#' @importFrom nimbleSCR rbernppACmovement_normal rbernppAC
#' 
#' @rdname getInits.s
#' @export
getInits.s <- function( y,
                        trapCoords,
                        lowerCoords,
                        upperCoords,
                        habitatGrid,
                        known.s = NULL,
                        baseIntensities = NULL,
                        sd = NULL,
                        detNums = NULL,
                        detIndices = NULL)
  {
  if(is.null(detNums)){
    detNums <- y[ ,1, ]
    nMaxDetectors <- (dim(y)[2]-1)/2
    detIndices <- y[ ,(nMaxDetectors + 2):dim(y)[2], ]
  } 
  
  ##-- Set-up objects  
  n.years <- dim(y)[3]
  n.individuals <- dim(y)[1]
  numGridRows <- dim(habitatGrid)[1]
  numGridCols <- dim(habitatGrid)[2]
  numHabWindows <- nrow(lowerCoords)
  if(is.null(baseIntensities)){baseIntensities <- rep(2,numHabWindows)}
  logIntensities <- log(baseIntensities)
  logSumIntensity <- log(sum(baseIntensities))
  s <- array(NA, c(n.individuals, 2, n.years))
  
  ##-- Check if trap coords vary among years
  if(length(dim(trapCoords)) == 2){
    trapCoords.xy <- array(NA, c(dim(trapCoords), n.years))
    for(t in 1:n.years){
      trapCoords.xy[ ,1,t] <- trapCoords[ ,1]
      trapCoords.xy[ ,2,t] <- trapCoords[ ,2]
    }#t
  } else {
    trapCoords.xy <- trapCoords
  }
  
  ##-- Find out who's detected at least once
  detected.any <- which(apply(y, 1, function(x)any(x[1, ] > 0)))
  if(!is.null(known.s)){ 
    detected.any <- unique(c(detected.any,
                             which(apply(known.s, 1, function(x)any(!is.na(x))))))
  }
  
  ##-- Find out who's never detected
  detected.never <- which(!(1:dim(y)[1]) %in% detected.any)
  
  ##-- Find out who's detected each year
  detected.t <- lapply(1:n.years, function(x) which(y[ ,1,x] > 0))
  
  
  ##-- Detected IDs : 
  ##-- Use the centroid of detections each year the ID is detected
  ##-- Sample a random location centered on the centroids of all detections 
  ##-- when the individual is not detected.
  for(i in detected.any){
    ##-- Fill in known locations
    if(!is.null(known.s)){s[i, , ] <- known.s[i, , ]}
    
    ##-- Calculate locations of detection centroids each year
    centroids.t <- lapply(1:n.years, function(t){
      if(detNums[i,t] > 0){
        detInd <- detIndices[i,detIndices[i, ,t] > 0,t]
        if(detNums[i,t] > 1){
          colMeans(trapCoords.xy[detInd, ,t])
        } else {
          trapCoords.xy[detInd, ,t]
        }
      }
    })
    
    ##-- For years detected:
    years.det <- which(sapply(detected.t, function(x)i %in% x))
    for(t in years.det){ s[i, ,t] <- centroids.t[[t]] }#t
      
    ##-- For years not detected:
    ##-- Identify the closest year 
    years.not.det <- which(is.na(s[i,1, ]))
    years.det <- which(!is.na(s[i,1, ]))
    for(t in years.not.det){
        closestInTime <- years.det[which.min(abs(years.det - t))]
        s[i, ,t] <- rbernppACmovement_normal(
          n = 1,
          lowerCoords = as.matrix(lowerCoords),
          upperCoords = as.matrix(upperCoords),
          s = s[i, ,closestInTime],
          sd = sd,
          baseIntensities = baseIntensities,
          habitatGrid = habitatGrid,
          numGridRows = numGridRows,
          numGridCols = numGridCols,
          numWindows = numHabWindows)
    }#t
  }#i
  
  
  ##-- Augmented IDs : 
  ##-- Sample a random walk 
  for(i in detected.never){
    s[i, ,1] <- rbernppAC(
      n = 1,
      lowerCoords = as.matrix(lowerCoords),
      upperCoords = as.matrix(upperCoords),
      logIntensities = logIntensities,
      logSumIntensity = logSumIntensity,
      habitatGrid = habitatGrid,
      numGridRows = numGridRows,
      numGridCols = numGridCols)
    ##-- t > 1 
    if(n.years>1){
      for(t in 2:n.years){
        s[i, ,t] <- rbernppACmovement_normal(
          n = 1,
          lowerCoords = as.matrix(lowerCoords),
          upperCoords = as.matrix(upperCoords),
          s = s[i, ,t-1],
          sd = sd,
          baseIntensities = baseIntensities,
          habitatGrid = habitatGrid,
          numGridRows = numGridRows,
          numGridCols = numGridCols,
          numWindows= numHabWindows)
      }#t
    }
  }#i
  
  return(s)
}