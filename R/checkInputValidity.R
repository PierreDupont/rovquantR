#' @title Detectors Assignment Function
#'
#' @description
#' \code{checkInputValidity} tests the validity of the input data.
#'
#' @name checkInputValidity
#'
#' @param data A \code{sf} object containing the detection data
#' @param detectors A \code{list} of sf objects containing the detector locations
#' @param subDetectors Optional; a \code{list} of sf objects containing the subdetectors locations (for PAB models only)
#' @param radius a \code{numeric} value denoting the maximum radius to consider for detector assignment. 
#' Any detection further than \code{radius}m from any detector will not be assigned to any detector. Instead it will get a 'NA'.
#'    id = NULL,
#' @param s array of activity center coordinates.
#' @param y detection histories.
#' @param z matrix of individual states.
#' @param lowerHabCoords lower coordinates of habitat grid cells.
#' @param upperHabCoords upper coordinates of habitat grid cells.
#' @param localHabWindowIndices matrix of IDs of the local habitat grid cells.
#' @param numLocalHabWindows vector of numbers of local habitat grid cells.
#' @param trapCoords matrix of detector coordinates
#' @param localTrapsIndices matrix of IDs of the local detectors.
#' @param localTrapsNum vector of numbers of local detectors.
#' @param resizeFactor default value= 1.
#' @param habitatGrid matrix of habitat grid cell IDs.
#' @param habitatGridLocal matrix of resized habitat grid cell IDs.
#' @param silent Logical; whether to print out outputs or not
#' @param printReport Logical; whether to print out a .pdf report or not
#' @param pathReport (optional) path to print the report
#'
#' @return A \code{.pdf} report with potentially problematic individual detections and/or movements. 
#'
#' @importFrom sf st_geometry st_multipolygon st_polygon
#' @importFrom grDevices pdf hcl.colors
#' @importFrom graphics points arrows
#' 
NULL
#' @rdname checkInputValidity
#' @export
checkInputValidity <- function(
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
  finalReport1 <- finalReport2 <- NULL
  n.years <- dim(s)[3]
  lengthYCombined <- dim(y)[2]
  nMaxDetectors <- (lengthYCombined-1)/2
  myCol <- grDevices:hcl.colors(n.years)
  
  if (printReport) {
    silent <- TRUE 
    if(is.null(pathReport))pathReport <- getwd()
    dir.create(pathReport,showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file.path(pathReport, "CheckDetections.pdf"), width = 10, height = 10)
  }
  if(is.null(id)){id <- 1:dim(s)[1]}
  
  ##-- Make habitat grid
  habGrid <- sf::st_multipolygon(
    lapply(1:nrow(lowerHabCoords),function(p){
      sf::st_polygon(list(
        matrix(c(lowerHabCoords[p,1],lowerHabCoords[p,2],
                 lowerHabCoords[p,1],upperHabCoords[p,2],
                 upperHabCoords[p,1],upperHabCoords[p,2],
                 upperHabCoords[p,1],lowerHabCoords[p,2],
                 lowerHabCoords[p,1],lowerHabCoords[p,2]),
               byrow = T, ncol = 2))
      )}))
  
  
  ##-- Loop over individuals
  for (i in id) {
    ##----- 1. Check the validity of the AC trajectory -----
    ##-- Initialize checks
    report <- NULL
    isNotOk_sID <- isNotOk_sID_resize <- isNotOk_local <- rep(FALSE, n.years)
    localWindows <- list()
    sID <- NULL
    
    ##-- First occasion
    sID[1] <- habitatGrid[trunc(s[i,2,1])+1, trunc(s[i,1,1])+1]
    isNotOk_sID[1] <- sID[1] <= 0
    if(isNotOk_sID[1]){
      report <- paste0(report, "\nHabitat grid cell not valid for Id: ", i, " ; t = 1")
    } 
    
    ##-- Following occasions
    for(t in 2:n.years){
      ##-- Check if proposed AC location is valid, i.e. proposed AC falls in a 
      ##-- valid habitat grid cell
      sID[t] <- habitatGrid[trunc(s[i,2,t])+1, trunc(s[i,1,t])+1]
      isNotOk_sID[t] <- sID[t] <= 0
      if(isNotOk_sID[t]){
        report <- paste0(report, "\nHabitat grid cell not valid for Id: ", i, " ; t = ", t)
      }    
      
      ##-- Check if previous AC location is valid, i.e. previous AC falls in a 
      ##-- valid resized habitat grid cell 
      sID_resize <- habitatGridLocal[trunc(s[i,2,t-1]/resizeFactor)+1,
                                     trunc(s[i,1,t-1]/resizeFactor)+1]
      isNotOk_sID_resize[t] <- sID_resize <= 0
      if(isNotOk_sID_resize[t]){
        report <- paste0(report, " \nRezised habitat grid cell not valid for Id: ", i, " ; t = ", t)
      } else {
        if(!is.null(localHabWindowIndices)){
          ##-- Check if local evaluation is valid, i.e. proposed AC falls in one 
          ##-- of the local habitat grid cells around previous AC
          localWindows[[t]] <- localHabWindowIndices[sID_resize, 1:numLocalHabWindows[sID_resize]]
          isNotOk_local[t] <- !sID[t] %in% localWindows[[t]] 
          if(isNotOk_local[t]){
            report <- paste0(report, "\nLocal evaluation of movement not valid for Id: ", i, " ; t = ", t)
          }
        }
      }
    }#t
    
    
    ##-- Print report and plot only if something is wrong 
    if(!is.null(report)){ 
      if(!silent){cat(paste0("\n      ===> Invalid AC movement found for Id ", i, ":"))}
      
      ##-- Plot habitat grid
      sRange <- apply(s[i,,], 1, range)
      maxDiff <- max(apply(sRange, 2, diff))
      sideSize <- maxDiff*2
      xMargin <- (sideSize - diff(sRange[ ,1]))/2
      yMargin <- (sideSize - diff(sRange[ ,2]))/2
      plot(-100,
           main = paste0("Invalid AC movement found for Id: ", i, " ; t = 1 to ", n.years),
           xlab ="x", ylab = "y", 
           xlim = c(sRange[1,1]-xMargin, sRange[2,1]+xMargin),
           ylim = c(sRange[1,2]-yMargin, sRange[2,2]+yMargin))
      plot(habGrid, col = "gray80", border = "gray40", add= T)
      
      ##-- Plot initial AC location
      if(isNotOk_sID[1]) {
        graphics::points(s[i,1,1], s[i,2,1], pch = 19, cex = 2, col = "red")
      } else {
        graphics::points(s[i,1,1], s[i,2,1], pch = 19, cex = 2, col = "green")
      }
      
      for(t in 2:n.years){ 
        ##-- Plot local habitat grid cells and proposed habitat habitat grid cell
        if(!isNotOk_sID[t]){
          if(isNotOk_local[t]){
            lapply(localWindows[[t]], function(p){
              plot(habGrid[[p]], add = T,
                   col = adjustcolor(col = "blue", alpha.f = 0.1))
            })
            plot(habGrid[[sID[t]]], add = T, col = adjustcolor(col = "red",alpha.f = 0.2))
          } else {
            plot(habGrid[[sID[t]]], add = T, col = adjustcolor(col = "forestgreen",alpha.f = 0.1))
          }
        }
        
        ##-- Plot AC movement
        graphics::arrows( x0 = s[i,1,t-1],
                x1 = s[i,1,t],
                y0 = s[i,2,t-1],
                y1 = s[i,2,t], col = myCol[t], length = 0.08, lwd = 1.5)
        
        ##-- Plot proposed AC location
        if(isNotOk_sID[t] | isNotOk_local[t]){
          graphics::points(s[i,1,t], s[i,2,t], pch = 19, col = "red")
        } else { 
          graphics::points(s[i,1,t], s[i,2,t], pch = 19, col = "green")
        }
        
        if(!silent){invisible(readline(prompt = "Press [enter] to continue"))}
      }#t
      
      ##-- Print report
      if(!silent){cat(report)}
      if(is.null(finalReport1)){
        finalReport1 <- paste0("\n      ===> Invalid AC movements found: ", report) 
      } else {
        finalReport1 <- paste0(finalReport1, report)
      }
    }#if
    
    
    
    
    ##----- 2. Check individual detections -----
    report <- NULL
    
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
        report <- paste0(report, "\nRezised habitat grid cell not valid for detections of Id = ", i, " ; t = ", t)
      } else {
        ##-- Identify corresponding number of local detectors
        isNotOk_localTraps[t] <- localTrapsNum[sID_resize,t] <= 0
        if(isNotOk_localTraps[t]){
          report <- paste0(report, "\nNo local detectors available for detections of Id = ", i, " ; t = ", t)
        } else {
          ##-- Identify corresponding local detectors
          localTraps[[t]] <- localTrapsIndices[sID_resize, 1:localTrapsNum[sID_resize,t],t]
          ##-- Check if all detections happen at local detectors 
          isNotOk_detIndices[t] <- !all(theseDets[[t]] %in% localTraps[[t]])
          if(isNotOk_detIndices[t]){
            ##-- Identify "bad" detections
            badDets[[t]] <- theseDets[[t]][!(theseDets[[t]] %in% localTraps[[t]])]
            report <- paste0(report, "\nDetections of Id = ", i, " ; t = ", t,
                             " at detectors = ", badDets[[t]], " are not valid!")
          }
        }
      }
    }#t
    
    
    ##-- Print report and plot only if something is wrong 
    if(!is.null(report)){
      
      if(!silent){cat(paste0("\n      ===> Invalid detections for Id ", i, ":"))}
      
      for(t in theseYears){
        
        ##-- Only if there is an issue
        if(isNotOk_sID_resize[t] | isNotOk_localTraps[t] | isNotOk_detIndices[t]){
          
          ##-- 1. Plot habitat grid
          if(isNotOk_sID_resize[t]){
            sRange <- apply(s[i, , ], 1, range)
          } else { 
            if(isNotOk_localTraps[t]){
              sRange <- apply(rbind(s[i, ,t],
                                    trapCoords[theseDets[[t]], ,t]), 2, range)
            } else {
              sRange <- apply(rbind(s[i, ,t],
                                    trapCoords[theseDets[[t]], ,t],
                                    trapCoords[localTraps[[t]], ,t]), 2, range) 
            }
          }
          maxDiff <- max(apply(sRange,2,diff))
          sideSize <- maxDiff*2
          xMargin <- (sideSize - diff(sRange[ ,1]))/2
          yMargin <- (sideSize - diff(sRange[ ,2]))/2
          plot( habGrid, col = "gray80", border = "gray40",
                main = paste0("Individual detections for Id = ", i, " ; t = ", t),
                xlab ="x", ylab = "y",
                xlim = c(sRange[1,1]-xMargin, sRange[2,1]+xMargin),
                ylim = c(sRange[1,2]-yMargin, sRange[2,2]+yMargin))
          
          ##-- 2. Plot detectors
          graphics::points(trapCoords[ , ,t], pch = 3, cex = 0.5)
          
          ##-- 3. Plot local detectors
          if (!isNotOk_localTraps[t] & !isNotOk_sID_resize[t]) {
            graphics::points( x = trapCoords[localTraps[[t]],1,t],
                              y = trapCoords[localTraps[[t]],2,t],
                              pch = 3, cex = 0.5, col = "blue")
          }
          
          ##-- 4. Plot detections
          if (isNotOk_localTraps[t]) { thisCol = "red" } else { thisCol = "green" }
          graphics::points( x = trapCoords[theseDets[[t]],1,t],
                            y = trapCoords[theseDets[[t]],2,t],
                            pch = 16, col = thisCol)
          ##-- If some detections outside local detectors 
          if (isNotOk_detIndices[t]) {
            graphics::points( x = trapCoords[badDets[[t]],1,t],
                              y = trapCoords[badDets[[t]],2,t],
                              pch = 16, col = "red")
          }
          
          ##-- 5. Plot AC location
          if (isNotOk_sID_resize[t]) { thisCol = "red" } else { thisCol = "green" }
          graphics::points( x = s[i,1,t], y = s[i,2,t],
                            pch = 19, cex = 1.5, col = thisCol)
          
          if (!silent) {invisible(readline(prompt = "Press [enter] to continue"))}
        }#if
      }#t
      
      ##-- Print report
      if (!silent) {cat(report);cat("\n")}
      if (is.null(finalReport2)) {
        finalReport2 <- paste0("\n\n      ===> Invalid detections found:", report)
      } else {
        finalReport2 <- paste0(finalReport2, report)
      }
    }
  }#i
  
  ##-- Close pdf report
  if (printReport) {
    if (!is.null(finalReport1) | !is.null(finalReport2)) {
      cat( paste0(finalReport1,finalReport2),
           file = file.path(pathReport, "CheckDetections.txt"))
    } else {
      print("All good! nothing to print out... ")
    }
    dev.off()
  }
}

NULL
#' @rdname checkInputValidity
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
    habitatGridLocal = NULL){
  
  ##-- Initial set-up & checks
  n.years <- dim(s)[3]
  lengthYCombined <- dim(y)[2]
  nMaxDetectors <- (lengthYCombined-1)/2

  if (is.null(id)) { id <- 1:dim(s)[1] }
  
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