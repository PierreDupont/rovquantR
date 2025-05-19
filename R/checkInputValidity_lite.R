#' @title Detectors Assignment Function
#'
#' @description
#' \code{checkInputValidity_lite} tests the validity of the input data.
#' 
#' @param data A \code{sf} object containing the detection data
#' @param detectors A \code{list} of sf objects containing the detector locations
#' @param subDetectors Optional; a \code{list} of sf objects containing the subdetectors locations (for PAB models only)
#' @param radius a \code{numeric} value denoting the maximum radius to consider for detector assignment. 
#' Any detection further than \code{radius}m from any detector will not be assigned to any detector. Instead it will get a 'NA'.
#'
#' @return A \code{sf} object with the x and y locations of the
#'  different detections with associated detector's Id 
#'
#' @rdname checkInputValidity_lite
#' @export
checkInputValidity_lite <- function(
    id = NULL,
    s,
    y,
    z,
    lowerHabCoords,
    upperHabCoords,
    localHabWindowIndices = NULL,
    numLocalHabWindows = NULL,
    trapCoords,
    localTrapsIndices = NULL,
    localTrapsNum = NULL,
    resizeFactor = 1,
    habitatGrid,
    habitatGridLocal = NULL,
    silent = FALSE,
    printReport = FALSE,
    pathReport = NULL){
  
  ##-- Initial set-up & checks
  # finalReport1 <- finalReport2 <- NULL
  n.years <- dim(s)[3]
  lengthYCombined <- dim(y)[2]
  nMaxDetectors <- (lengthYCombined-1)/2
  #myCol <- grDevices:hcl.colors(n.years)
  
  # if (printReport) {
  #   silent <- TRUE 
  #   if(is.null(pathReport))pathReport <- getwd()
  #   dir.create(pathReport,showWarnings = FALSE, recursive = TRUE)
  #   grDevices::pdf(file.path(pathReport, "CheckDetections.pdf"), width = 10, height = 10)
  # }
  if (is.null(id)) { id <- 1:dim(s)[1] }
  
  # ##-- Make habitat grid
  # habGrid <- sf::st_multipolygon(
  #   lapply(1:nrow(lowerHabCoords), function(p) {
  #     sf::st_polygon(list(
  #       matrix(c(lowerHabCoords[p,1],lowerHabCoords[p,2],
  #                lowerHabCoords[p,1],upperHabCoords[p,2],
  #                upperHabCoords[p,1],upperHabCoords[p,2],
  #                upperHabCoords[p,1],lowerHabCoords[p,2],
  #                lowerHabCoords[p,1],lowerHabCoords[p,2]),
  #              byrow = T, ncol = 2))
  #     )}))
  # 
  
  ##-- Loop over individuals
  for (i in id) {
    ##----- 1. Check the validity of the AC trajectory -----
    
    ##-- Initialize checks
    # report <- NULL
    isNotOk_sID <- isNotOk_sID_resize <- isNotOk_local <- rep(FALSE, n.years)
    localWindows <- list()
    sID <- NULL
    
    ##-- First occasion
    sID[1] <- habitatGrid[trunc(s[i,2,1])+1, trunc(s[i,1,1])+1]
    isNotOk_sID[1] <- sID[1] <= 0
    if(isNotOk_sID[1]){
      cat(paste0("\nHabitat grid cell not valid for Id: ", i, " ; t = 1"))
    } 
    
    ##-- Following occasions
    for(t in 2:n.years){
      ##-- Check if proposed AC location is valid, i.e. proposed AC falls in a 
      ##-- valid habitat grid cell
      sID[t] <- habitatGrid[trunc(s[i,2,t])+1, trunc(s[i,1,t])+1]
      isNotOk_sID[t] <- sID[t] <= 0
      if(isNotOk_sID[t]){
        cat(paste0("\nHabitat grid cell not valid for Id: ", i, " ; t = ", t))
      }    
      
      ##-- Check if previous AC location is valid, i.e. previous AC falls in a 
      ##-- valid resized habitat grid cell 
      sID_resize <- habitatGridLocal[trunc(s[i,2,t-1]/resizeFactor)+1,
                                     trunc(s[i,1,t-1]/resizeFactor)+1]
      isNotOk_sID_resize[t] <- sID_resize <= 0
      if (isNotOk_sID_resize[t]) {
        cat(paste0(" \nRezised habitat grid cell not valid for Id: ", i, " ; t = ", t))
      } else {
        if(!is.null(localHabWindowIndices)){
          ##-- Check if local evaluation is valid, i.e. proposed AC falls in one 
          ##-- of the local habitat grid cells around previous AC
          localWindows[[t]] <- localHabWindowIndices[sID_resize, 1:numLocalHabWindows[sID_resize]]
          isNotOk_local[t] <- !sID[t] %in% localWindows[[t]] 
          if(isNotOk_local[t]){
            cat(paste0("\nLocal evaluation of movement not valid for Id: ", i, " ; t = ", t))
          }
        }
      }
    }#t
    

    ##----- 2. Check individual detections -----

    ##-- List detections per year
    numDets <- y[i,1, ]
    detIndices <- y[i,(nMaxDetectors+2):lengthYCombined, ]
    
    ##-- Initialize checks
    isNotOk_sID_resize <- isNotOk_localTraps <- isNotOk_detIndices <- rep(FALSE,n.years)
    localTraps <- theseDets <- badDets <- list()
    
    ##-- Loop over years with detections only
    theseYears <- which(numDets > 0)
    for(t in theseYears){
      ##-- Identify detector indices
      theseDets[[t]] <- detIndices[1:numDets[t],t]
      
      ##-- Identify resized habitat grid cell
      sID_resize <- habitatGridLocal[trunc(s[i,2,t]/resizeFactor)+1, trunc(s[i,1,t]/resizeFactor)+1]
      isNotOk_sID_resize[t] <- sID_resize <= 0
      if(isNotOk_sID_resize[t]){
        cat(paste0("\nRezised habitat grid cell not valid for detections of Id = ", i, " ; t = ", t))
      } else {
        ##-- Identify corresponding number of local detectors
        isNotOk_localTraps[t] <- localTrapsNum[sID_resize,t] <= 0
        if(isNotOk_localTraps[t]){
          cat(paste0("\nNo local detectors available for detections of Id = ", i, " ; t = ", t))
        } else {
          ##-- Identify corresponding local detectors
          localTraps[[t]] <- localTrapsIndices[sID_resize, 1:localTrapsNum[sID_resize,t],t]
          ##-- Check if all detections happen at local detectors 
          isNotOk_detIndices[t] <- !all(theseDets[[t]] %in% localTraps[[t]])
          if(isNotOk_detIndices[t]){
            ##-- Identify "bad" detections
            badDets[[t]] <- theseDets[[t]][!(theseDets[[t]] %in% localTraps[[t]])]
            cat(paste0("\nDetections of Id = ", i, " ; t = ", t,
                             " at detectors = ", badDets[[t]], " are not valid!"))
          }
        }
      }
    }#t
  }#i
  

}