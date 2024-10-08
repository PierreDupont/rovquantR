#' @title Detectors Assignment Function
#'
#' @description
#' \code{AssignDetectors} identifies and assigns the closest detector to each detection in \code{myData}.
#'
#' @param myData A \code{sf} object containing the detection data
#' @param myDetectors A \code{list} of sf objects containing the detector locations
#' @param mysubDetectors Optional; a \code{list} of sf objects containing the subdetectors locations (for PAB models only)
#' @param radius a \code{numeric} value denoting the maximum radius to consider for detector assignment. 
#' Any detection further than \code{radius}m from any detector will not be assigned to any detector. Instead it will get a 'NA'.
#'
#' @return A \code{sf} object with the x and y locations of the
#'  different detections with associated detector's Id 
#'  
#' @importFrom RANN nn2 
#' @importFrom sf st_geometry 
#'
#' @rdname AssignDetectors
#' @export
AssignDetectors <- function( 
    myData,                              
    myDetectors,
    mysubDetectors = NULL,
    radius = 20000)
{
  #myData$record.id <- 1:dim(myData)[1]
  
  ##-- RETRIEVE NUMBER OF YEARS 
  if(is.null(myData$Year)){myData$Year <- 1}
  years <- min(unique(myData$Year)):max(unique(myData$Year))
  
  ##-- SET UP THE DETECTORS LIST
  if(!inherits(myDetectors,"list")){myDetectors <- list(myDetectors)}
  
  ##-- CHECK THAT THE NUMBER OF YEARS MATCH THE LENGTH OF THE DETECTORS LIST
  if(length(myDetectors) != length(years)){
    myDetectors <- lapply(1:length(years), function(x) myDetectors[[1]])
  }
  
  ##-- IF PAB PROCESS
  if(!is.null(mysubDetectors)){
    ##-- SET UP THE SUB-DETECTORS LIST
    if(!inherits(mysubDetectors,"list")){
      mysubDetectors <- list(mysubDetectors)
    }
    
    if(length(mysubDetectors) != length(years)){
      mysubDetectors <- lapply(1:length(years), function(x) mysubDetectors[[1]])
    }
    ##-- INITIALIZE DETECTORS AND SUB-DETECTORS ID COLUMNS
    myData$sub.detector <- NA 
    myData$Detector <- NA
    for(t in 1:length(years)){
      
      ##-- Get coordinate matrices
      detector_coords <- do.call(rbind, sf::st_geometry(mysubDetectors[[t]]))
      detection_coords <- do.call(rbind, sf::st_geometry(myData[myData$Year == years[t], ]))
      
      if(!is.null(detection_coords)){
        
        ##-- Fast nearest neighbour search
        closest <- RANN::nn2(data = detector_coords,
                             query = data.frame( detection_coords[ ,1],
                                                 detection_coords[ ,2]),
                             k = 1, searchtype = "radius", radius = radius)
        closest$nn.idx[closest$nn.idx==0] <- NA
        myData$sub.detector[myData$Year == years[t]] <- closest$nn.idx
        
        main.cell.id <- mysubDetectors[[t]]$main.cell.id[closest$nn.idx]     
        myData$Detector[myData$Year == years[t]] <- unlist(
          lapply(1:length(main.cell.id), function(x){
            out <- NA
            if(!is.na(main.cell.id[x])){
              out <- which(myDetectors[[t]]$main.cell.id %in% main.cell.id[x])
            }
            return(out)
          }))
      }#if
    }#t
  } else {
    myData$Detector <- NA
    for(t in 1:length(years)){
      ##-- Get coordinate matrices
      detector_coords <- do.call(rbind, sf::st_geometry(myDetectors[[t]]))
      detection_coords <- do.call(rbind, sf::st_geometry(myData[myData$Year == years[t],]))
      
      if(!is.null(detection_coords)){
        ##-- Fast nearest neighbour search
        closest <- RANN::nn2( data = detector_coords,
                              query = data.frame(detection_coords[,1],
                                                 detection_coords[,2]),
                              k = 1, searchtype = "radius", radius = radius)
        
        closest$nn.idx[closest$nn.idx==0] <- NA
        myData$Detector[myData$Year == years[t]] <- closest$nn.idx
      }#if 
    }
  }
  
  myData <- myData[!is.na(myData$Detector),]
  myData$Id <- as.character(myData$Id)
  
  ##-- Output return
  if(!is.null(mysubDetectors)){                                  
    n.trials <- lapply(mysubDetectors,function(x)table(x$main.cell.id))
    return(list( myData.sp = myData,
                 n.trials = n.trials))
  } else {
    return(myData)
  }
}
