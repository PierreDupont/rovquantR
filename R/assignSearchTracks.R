#' Assign DNA samples to search tracks 
#' 
#' The \code{assignSearchtTracks} function is a R function used to match RovBase
#' DNA samples to the corresponding search track. Matching is based on date 
#' (exact match) and distance (within a given distance; see the \code{dist} argument below).
#' 
#' @name assignSearchtTracks
#' 
#' @param data A \code{sf} object containing the detection data.
#' @param tracks A \code{sf} object containing the search tracks.
#' @param dist a \code{numeric} value denoting the maximum distance (in meters) to consider for tracks assignment. 
#' Any detection further than \code{dist} m from any track will not be assigned to any track. Instead it will get a 'NA'.
#' @param progress.bar a \code{logical}, whether to display a progress bar or not. 
#' 
#' @return A \code{sf} object with the x and y locations of the different detections with associated tracks ID. 
#'
#' @author Pierre Dupont
#' 
#' @importFrom sf st_intersection st_as_sf st_buffer st_distance
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
NULL
#' @rdname assignSearchtTracks
#' @export
assignSearchtTracks <- function( data,
                                 tracks,
                                 dist = 500,
                                 progress.bar = T)
{
  ##-- Initialise track ID and distance
  data$track_RovbaseID <- NA
  data$track_dist <- NA
  
  ##-- Buffer DNA detections
  buffData <- st_buffer(st_as_sf(data), dist = dist)
  
  ##-- Set-up progress bar 
  if(progress.bar){
    pb = utils::txtProgressBar( min = 1,
                                max = nrow(data),
                                initial = 0,
                                style = 3) 
  }
  
  ##-- ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
  for(i in 1:nrow(data)){
    
    ##-- INTERSECT BUFFERED POINT WITH TRACKS
    tmpTRACKS <- st_intersection(filter(tracks, Dato == data$Date[i]), buffData[i, ])
    
    ##-- If no match, move on...
    if(nrow(tmpTRACKS) == 0){
      next
    } else {
      print(i)
      ##-- ... find the closest track
      dist <- st_distance(data[i, ], tmpTRACKS, by_element = F)
      data$track_RovbaseID[i] <- tmpTRACKS$RovbaseID[which.min(dist)]
      data$track_dist[i] <- min(dist)
    }
    ## Print progress bar
    if(progress.bar){ utils::setTxtProgressBar(pb,i) }
  }#i
  if(progress.bar)close(pb)  
return(data)
}
  
  
  
  
  