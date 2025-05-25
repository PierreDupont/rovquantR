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
#' @importFrom RANN nn2
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
                        sd,
                        detNums = NULL,
                        detIndices = NULL,
                        radius = NULL){
  
  ## ----- 1. Function set-up -----
  
  ##-- Inner function to create equally spaced points 
  splitNequal <- function(start,end,n) {
    numSegments <- n + 1
    deltax <- (end[1]-start[1])/numSegments
    deltay <- (end[2]-start[2])/numSegments
    xcoords <- start[1] + deltax + (1:n-1)*deltax
    ycoords <- start[2] + deltay + (1:n-1)*deltay
    out <- rbind(xcoords,ycoords)
    dimnames(out) <- list("coords" = c("x","y"),
                          "years" = 1:n)
    return(out)
  }
  
  ##-- Check detection matrix format
  if (is.null(detNums)) {
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
  if (is.null(baseIntensities)) { baseIntensities <- rep(2,numHabWindows) }
  logIntensities <- log(baseIntensities)
  logSumIntensity <- log(sum(baseIntensities))
  
  ##-- Check if trap coords vary among years
  if (length(dim(trapCoords)) == 2) {
    trapCoords.xy <- array(NA, c(dim(trapCoords), n.years))
    for(t in 1:n.years){
      trapCoords.xy[ ,1,t] <- trapCoords[ ,1]
      trapCoords.xy[ ,2,t] <- trapCoords[ ,2]
    }#t
  } else {
    trapCoords.xy <- trapCoords
  }
  
  ##-- Fill in known locations (if any)
  if(!is.null(known.s)) { s <- known.s } else {s <- array(NA, c(n.individuals, 2, n.years))}
  
  
  
  ## ----- 2. For detected IDs -----
  
  ##-- Find out who's detected at least once
  detected.any <- which(apply(y, 1, function(x)any(x[1, ] > 0)))
  if (!is.null(known.s)) { 
    detected.any <- unique(c(detected.any,
                             which(apply(known.s, 1, function(x)any(!is.na(x))))))
  }
  
  ##-- Find out who's detected each year
  detected.t <- lapply(1:n.years, function(x) which(y[ ,1,x] > 0))
  
  for(i in detected.any){
    
    ## -----    2.1. For years detected ... -----
    
    ##-- We use the location of detection centroids each year
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
    for(t in which(sapply(detected.t, function(x)i %in% x))){ s[i, ,t] <- centroids.t[[t]] }#t
    
    
    
    ## -----    2.2. Before first detection ... -----
    years.det <- which(!is.na(s[i,1,]))
    ##-- We sample backward from first detection
    years.before.det <- which(1:n.years < min(years.det)) 
    if(length(years.before.det) > 0){
      for(t in rev(years.before.det)){
        s[i, ,t] <- rbernppACmovement_normal(
          n = 1,
          lowerCoords = as.matrix(lowerCoords),
          upperCoords = as.matrix(upperCoords),
          s = s[i, ,t+1],
          sd = sd,
          baseIntensities = baseIntensities,
          habitatGrid = habitatGrid,
          numGridRows = numGridRows,
          numGridCols = numGridCols,
          numWindows = numHabWindows)
      }#t
    }#if
    
    
    
    ## -----    2.3. After last detection ... -----
    
    ##-- We sample forward from last detection
    years.after.det <- which(1:n.years > max(years.det))
    if(length(years.after.det) > 0){
      for(t in years.after.det){
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
          numWindows = numHabWindows)
      }#t
    }#if
    
    
    
    ## -----    2.4. In between ... -----
    ##-- We sample activity center locations on a straight line
    if(length(years.det)>1) {
      ##-- First, we need to identify gaps >= 2 years
      start.years <- years.det[1:(length(years.det)-1)] 
      end.years <- years.det[2:length(years.det)] 
      gaps <- end.years - start.years - 1
      
      start.years <- start.years[gaps > 0] 
      end.years <- end.years[gaps > 0] 
      gaps <- gaps[gaps > 0]
      
      ##-- Identify the closest year ...
      if(length(gaps) > 0){
        for(g in 1:length(gaps)){
          
          ##-- Split the trajectory between the two extremities into equal segments
          tmp <- splitNequal( 
            start = s[i, ,start.years[g]],
            end = s[i, ,end.years[g]],
            n = gaps[g])
          
          ##-- Check if new points are in valid habitat
          outsideHab <- apply(tmp, 2, function(p){
            habitatGrid[trunc(p[2])+1, trunc(p[1])+1] <= 0
          })
          
          ##-- If not ...
          if(any(outsideHab)){
            ##-- Find out the closest habitat grid cell
            whichOut <- which(outsideHab)
            closest <- RANN::nn2( data = lowerCoords+0.5,
                                  query = data.frame( tmp[1,whichOut],
                                                      tmp[2,whichOut]),
                                  k = 1, 
                                  searchtype = "radius",
                                  radius = radius)
            
            ##-- Sample uniform AC location in closest habitat grid cell
            tmp[1,which(outsideHab)] <- runif( n = length(whichOut),
                                               min = lowerCoords[closest$nn.idx,1],
                                               max = upperCoords[closest$nn.idx,1])
            tmp[2,which(outsideHab)] <- runif( n = length(whichOut),
                                               min = lowerCoords[closest$nn.idx,2],
                                               max = upperCoords[closest$nn.idx,2])
          }#if
          
          s[i, ,(start.years[g]+1):(end.years[g]-1)] <- tmp
          
        }#g
      }#if
    }#if

  }#i
  
  
  ## ----- 3. For augmented IDs ------
  
  ##-- First, we find out who's never detected
  detected.never <- which(!(1:dim(y)[1]) %in% detected.any)
  
  ##-- ...then, we sample a random walk trajectory
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
  
  
  
  ## ----- 4. Output ------
  
  return(s)
  
}