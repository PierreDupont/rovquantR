#' @title Create a trap response covariate
#'
#' @description
#' \code{makeTrapResponseCov} ide and assigns the closest detector to each detection in \code{data}.
#'
#' @param data A \code{sf} object containing the detection data
#' @param data.dead OPTIONAL: A second \code{sf} object containing the dead recovery data
#' @param IDs OPTIONAL: A \code{vector} of IDs to return the trap response covariate for.
#' Any detection further than \code{radius}m from any detector will not be assigned to any detector. Instead it will get a 'NA'.
#'
#' @return A \code{matrix} object with the x and y locations of the
#'  different detections with associated detector's Id 
#'  
#' @importFrom RANN nn2 
#' @importFrom sf st_geometry 
#'
#' @rdname makeTrapResponseCov
#' @export

makeTrapResponseCov <- function(
    data,
    data.dead = NULL,
    IDs = NULL
    ){
  
  ##-- LIST INDIVIDUALS
  if(is.null(IDs)){
    IDs <- unique(c(as.character(data$Id), as.character(data.dead$Id)))
  }
  
  ##-- LIST YEARS
  #Years <- sort(unique(c(data$Year, data.dead$Year)))#[RB] ADDED SORTING 2019-07-30 - reason: ran into problems with bear data
  Years <- min(data$Year):max(data$Year)
  
  ##-- CREATE A DUMMY DATASET
  dummy <- data[ ,c("Id","Year")]
  dummy$Id <- factor(dummy$Id, IDs)
  dummy$Year <- factor(dummy$Year, Years)
  temp <- table(dummy$Id, dummy$Year)
  
  ##-- CREATE THE TRAP RESPONSE COVARIATE MATRIX
  ## [PD] CHANGE TO   tab.ar <- matrix(NA, dim(temp)[1], dim(temp)[2]) ???
  tab.ar <- matrix(0, dim(temp)[1], dim(temp)[2])
  tab.ar[ ,2:dim(temp)[2]] <- temp[ ,1:(dim(temp)[2]-1)]
  tab.ar[tab.ar > 0] <- 1
  dimnames(tab.ar) <- dimnames(temp)
  
  return(tab.ar)
}

