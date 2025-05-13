#' @title Detectors Assignment Function
#'
#' @description
#' \code{assignDetectors} identifies and assigns the closest detector to each detection in \code{data}.
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
#' @importFrom RANN nn2 
#' @importFrom sf st_geometry 
#'
#' @rdname assignDetectors
#' @export
assignDetectors <- function( 
    data,                              
    detectors,
    subDetectors = NULL,
    radius = 20000)
{

  ##-- RETRIEVE NUMBER OF YEARS 
  if(is.null(data$Year)){data$Year <- 1}
  years <- min(unique(data$Year)):max(unique(data$Year))
  
  ##-- SET UP THE DETECTORS LIST
  if(!inherits(detectors,"list")){detectors <- list(detectors)}
  
  ##-- CHECK THAT THE NUMBER OF YEARS MATCH THE LENGTH OF THE DETECTORS LIST
  if(length(detectors) != length(years)){
    detectors <- lapply(1:length(years), function(x) detectors[[1]])
  }
  
  ##-- IF PAB PROCESS
  if(!is.null(subDetectors)){
    ##-- SET UP THE SUB-DETECTORS LIST
    if(!inherits(subDetectors,"list")){
      subDetectors <- list(subDetectors)
    }
    
    if(length(subDetectors) != length(years)){
      subDetectors <- lapply(1:length(years), function(x) subDetectors[[1]])
    }
    ##-- INITIALIZE DETECTORS AND SUB-DETECTORS ID COLUMNS
    data$sub.detector <- NA 
    data$Detector <- NA
    for(t in 1:length(years)){
      
      ##-- Get coordinate matrices
      detector_coords <- do.call(rbind, sf::st_geometry(subDetectors[[t]]))
      detection_coords <- do.call(rbind, sf::st_geometry(data[data$Year == years[t], ]))
      
      if(!is.null(detection_coords)){
        
        ##-- Fast nearest neighbour search
        closest <- RANN::nn2(data = detector_coords,
                             query = data.frame( detection_coords[ ,1],
                                                 detection_coords[ ,2]),
                             k = 1, searchtype = "radius", radius = radius)
        closest$nn.idx[closest$nn.idx==0] <- NA
        data$sub.detector[data$Year == years[t]] <- closest$nn.idx
        
        main.cell.id <- subDetectors[[t]]$main.cell.id[closest$nn.idx]     
        data$Detector[data$Year == years[t]] <- unlist(
          lapply(1:length(main.cell.id), function(x){
            out <- NA
            if(!is.na(main.cell.id[x])){
              out <- which(detectors[[t]]$main.cell.id %in% main.cell.id[x])
            }
            return(out)
          }))
      }#if
    }#t
  } else {
    data$Detector <- NA
    for(t in 1:length(years)){
      ##-- Get coordinate matrices
      detector_coords <- do.call(rbind, sf::st_geometry(detectors[[t]]))
      detection_coords <- do.call(rbind, sf::st_geometry(data[data$Year == years[t],]))
      
      if(!is.null(detection_coords)){
        ##-- Fast nearest neighbour search
        closest <- RANN::nn2( data = detector_coords,
                              query = data.frame(detection_coords[,1],
                                                 detection_coords[,2]),
                              k = 1, searchtype = "radius", radius = radius)
        
        closest$nn.idx[closest$nn.idx==0] <- NA
        data$Detector[data$Year == years[t]] <- closest$nn.idx
      }#if 
    }
  }
  
  data <- data[!is.na(data$Detector),]
  data$Id <- as.character(data$Id)
  
  ##-- Output return
  if(!is.null(subDetectors)){                                  
    n.trials <- lapply(subDetectors,function(x)table(x$main.cell.id))
    return(list( data.sp = data,
                 n.trials = n.trials))
  } else {
    return(data)
  }
}
