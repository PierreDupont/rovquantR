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
    IDs = NULL)
{
  
  # ##-- LIST ALL INDIVIDUALS
  # if(is.null(IDs)){
  #   dummyId <- unique(as.character(data$Id))
  #   if(!is.null(data.dead)){
  #     dummyId <- unique(c(as.character(data$Id), as.character(data.dead$Id)))
  #   }
  # } else {
  #   dummyId <- IDs
  #   data <- data[data$Id %in% dummyId, ]
  #   data$Id <- factor(as.character(data$Id), levels = unique(as.character(data$Id)))
  # }
  # N_dummyId <- length(dummyId)
  # 
  # ##-- LIST ALL YEARS
  # Years <- range(unique(data$Year))
  # dummyYears <- Years[1]:Years[2]
  # N_dummyYears <- length(dummyYears)
  # 
  # ##-- CREATE A DUMMY DATASET
  # D_Id <- rep(dummyId, N_dummyYears)
  # D_Year <- rep.int(dummyYears, rep(N_dummyId, N_dummyYears))
  # dummy <- cbind.data.frame(Id = D_Id, Year = D_Year)
  # 
  # ##-- CREATE THE COMBINED DETECTION ARRAY OF ALIVE & DUMMY DETECTIONSdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
  # dummy <- rbind(dummy, data.frame(data)[ ,c("Id", "Year")])
  # tab <- table(dummy$Id, dummy$Year)
  # tab.ar <- matrix(0, dim(tab)[1], dim(tab)[2])
  # tab.ar[ ,2:dim(tab)[2]] <- tab[ ,1:(dim(tab)[2]-1)] - 1
  # tab.ar[tab.ar > 0] <- 1
  # dimnames(tab.ar) <- dimnames(tab)
  
  
  
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

