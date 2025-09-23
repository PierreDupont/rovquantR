#' Assign DNA samples to search tracks 
#' 
#' the \code{assignSearchtTracks} function is a R function used to match RovBase
#' DNA samples to the corresponding GPS search track. Matching is based on date 
#' (exact match) and distance (within a given threshold).
#' 
#' @name assignSearchTracks
#' 
#' @param data A \code{sf} object containing the coordinates of the detections to be designated to tracks.
#' @param tracks A \code{sf} object containing the search tracks 
#' @param progress.bar a \code{logical}, whether to display a progress bar or not. 
#' 
#' @return A \code{sf} object with the x and y locations of the different detections with associated tracks Id 
#'
#' @author Pierre Dupont
#' 
#' @importFrom sf st_distance 
#' @importFrom dplyr filter 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
NULL
#' @rdname assignSearchTracks
#' @export
assignSearchTracks <- function( 
    data,
    tracks,
    #dist = 500,
    progress.bar = T)
{
  ##-- Initialise track ID and distance
  data$trackID <- NA
  data$trackDist <- NA
  
  # ##-- Buffer DNA detections
  # buffData <- st_buffer(data, dist = dist)
  
  ##-- Set-up progress bar 
  if(progress.bar){
    pb = utils::txtProgressBar( min = 1,
                                max = nrow(data),
                                initial = 0,
                                style = 3) 
  }
  
  ##-- ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
  for(i in 1:nrow(data)){
    
    ##-- Print progress bar
    if(progress.bar){ utils::setTxtProgressBar(pb,i) }
    
    ##-- Subset search tracks to the day the sample was collected
    theseTracks <- tracks %>% filter(Dato == as.character(data$Date[i]))
    
    ##-- If no tracks on the same day, move to next sample...
    if(nrow(theseTracks) == 0){next}
    
    ##-- INTERSECT BUFFERED POINT WITH TRACKS RAN THE SAME DAY
    # theseTracks <- st_intersection(
    #   x = tracks[whichSameDate, ],
    #   Y = buffData[i, ])
    # ##-- If no match, move on...
    # if(nrow(theseTracks) == 0){next}
    
    ##-- ... Else, find the closest track, and retrieve track ID and distance
    dist <- st_distance(data[i, ], theseTracks, by_element = F)
    data$trackID[i] <- theseTracks$RovbaseID[which.min(dist)]
    data$trackDist[i] <- min(dist)
  }#i
  
  if(progress.bar){ close(pb) }
  return(data)
}




