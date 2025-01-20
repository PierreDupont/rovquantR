#' Assign DNA samples to search tracks bite-size NIMBLE MCMC algorithms
#' 
#' the \code{assignSearchtTracks} function is a R function used to match RovBase
#' DNA samples to the corresponding search track. Matching is based on date 
#' (exact match) and distance (within a given threshold).
#' samples on the fly.
#' 
#' @name assignSearchtTracks
#' 
#' @param mcmc a \code{NIMBLE MCMC algorithm} object. 
#' @param model a \code{NIMBLE model} object.
#' @param conf a \code{NIMBLE MCMC configuration} object.
#' @param bite.size an \code{integer} denoting the number of MCMC iterations in 
#' each MCMC bite. 
#' @param bite.number an \code{integer} denoting the number of MCMC bites to run.
#' @param path a \code{character string} denoting the path where MCMC bite outputs are
#'  to be saved as .RData files (when using \code{runMCMCbites} or \code{restartMCMCbites}) 
#'  or looked for (when using \code{collectMCMCbites}). 
#' @param save.rds (optional) a \code{logical value}. if TRUE, the state of the
#' nimble model is saved at the end of every MCMC bite. This is later used to restart 
#' the MCMC sampling process if necessary using \code{restartMCMCbites}.
#' @param path.rds (optional) a \code{path} to the folder containing the last MCMC model state file.
#' @param burnin an \code{integer} denoting the number of MCMC bites to be removed 
#' as burn-in. 
#' @param pattern a \code{character string} denoting the name of the object containing
#' MCMC samples in each .R Data bite file.
#' @param param.omit a \code{character vector} denoting the names of parameters 
#' that are to be ignored when combining MCMC bites.
#' @param progress.bar a \code{logical value}. if TRUE, a separate progress bar is 
#' printed for each MCMC chain.
#' @param samplerDef a  \code{NIMBLE sampler definition} object.
#' @param modelState \code{NIMBLE model state} object.
#' @param mcmcState \code{NIMBLE MCMC state} object.
#' 
#' @return \code{runMCMCbites} and \code{restartMCMCbites} run and save locally 
#' nimble MCMC posterior samples. \code{collectMCMCbites} combines multiple MCMC
#'  bites (and MCMC chains) and returns them in the commnonly used \code{mcmc.list}
#'  format from the \code{coda} package.
#'
#' @author Pierre Dupont
#' 
#'@importFrom sf st_intersection st_as_sf st_buffer
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
NULL
#' @rdname assignSearchtTracks
#' @export
assignSearchtTracks <- function( data = myFilteredData.sp$alive,
                                 tracks = TRACKS_YEAR,
                                 dist = 750,
                                 progress.bar = T)
{
  ##-- Initialise track ID and distance
  data$track_RovbaseID <- NA
  data$track_dist <- NA
  
  ##-- Buffer DNA detections
  buffData <- st_as_sf(data) %>% st_buffer(., dist = dist)
  
  ##-- Set-up progress bar 
  if(progress.bar){
    pb = utils::txtProgressBar( min = 1,
                                max = nrow(data),
                                initial = 0,
                                style = 3) 
  }
  
  ##-- ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
  for(i in 1:nrow(data)){
    
    # ##-- INTERSECT BUFFERED POINT WITH TRACKS
    # t <- which(years %in% data$Year[i])
    # whichSameDate <- which(as.character(tracks[[t]]$Dato) == as.character(data$Date[i]))
    # tmpTRACKS <- st_intersection(tracks[[t]][whichSameDate, ],
    #                              buffData[i, ])
    
    ##-- INTERSECT BUFFERED POINT WITH TRACKS
    tmpTRACKS <- tracks %>%
      filter(., Dato == data$Date[i]) %>%
      st_intersection(., buffData[i, ])
    
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
    if(progress.bar){ utils::setTxtProgressBar(pb,b) }
  }#i
  if(progress.bar)close(pb)  
return(data)
}
  
  
  
  
  