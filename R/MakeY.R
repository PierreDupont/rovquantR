#' @title Generate Detection Histories Arrays
#'
#' @description
#' \code{makeY} returns a \code{array} object with the detections (presence/absence or counts) of individuals over detectors and time.
#'
#' @param data A \code{SpatialPointsDataFrame} of detections with ID, Year and Detector columns.
#' @param detectors A \code{list} of SpatialPointsDataFrames containing the detector locations.
#' @param method A \code{Character} object defining the array to be returned ("Binomial" , "Poisson" , "Bernoulli").
#' @param data2 OPTIONAL: A second \code{SpatialPointsDataFrame} of detections with ID, Year and Detector columns.
#' @param detectors2 OPTIONAL: A \code{SpatialPointsDataFrame}of SpatialPointsDataFrames containing the detector locations.
#' @param method2 OPTIONAL: A \code{Character} object defining the array to be returned ("Binomial" , "Poisson" , "Bernoulli").
#' @param IDs OPTIONAL: A \code{vector} of IDs to return the y array for.
#' @param returnIdvector OPTIONAL: A \code{logical} to return the vector of individual ID.
#' 
#' @return A \code{array} object containing the number of detections per individual (rows), detector (columns), and years (third dimensions).
#' 
#' @author Pierre Dupont
#' 
#' @export
makeY <- function( 
    data,                 
    detectors,           
    method = "Bernoulli",
    data2 = NULL, 
    detectors2 = NULL,
    method2 = "Bernoulli",
    IDs = NULL,
    returnIdvector = FALSE)
{
  
  ##-- PREPARE THE DATA 
  data <- as.data.frame(data)
  data2 <- as.data.frame(data2)
  
   # if(sum(duplicated(data2$Id))>0){stop("INDIVIDUALS CANNOT BE DEAD TWICE!!!!")}
  if(!sum(class(detectors) %in% "list")==1){detectors <- list(detectors)}
  if(!sum(class(detectors2) %in% "list")==1){detectors2 <- list(detectors2)}
  
  ##-- AGGREGATE MULTIPLE DETECTIONS AT THE SAME SUBDETECTOR TO ONE DETECTION ONLY
  if(method == "Binomial"){
    data <- data[!duplicated(data[ ,c("sub.detector","Year","Id")]), ] 
    data <- droplevels(data)
  }
  
  ##-- NUMBER OF INDIVIDUALS
  if(is.null(IDs)){
    IDs <- unique(c(as.character(data$Id), as.character(data2$Id)))
  }
  #N_Id <- length(IDs)
  
  ##-- NUMBER OF YEARS
  Years <- sort(unique(c(data$Year, data2$Year)))#[RB] ADDED SORTING 2019-07-30 - reason: ran into problems with bear data
  #N_Years <- length(Years)
  
  ##-- NUMBER OF DETECTORS 
  ## [PD] : need to check if this does not introduce a problem (e.g. if detectors have names and are not simply numbered)
  ## Here we assume detector 'names' are always 1:n.detectors[t] each year
  ## This is the most likely scenario but we mau have situations where this is not the case.
  ## Suggested fix: use the actual detector names, whether they are reset every year or not.
  ## Detectors <- unique(data$detectors) 
  N_Detec <- max(unlist(lapply(detectors, function(x) dim(x)[1])))
  Detectors <- 1:N_Detec
  
  ##-- CREATE A DUMMY DATASET
  dummy <- data[ ,c("Id","Year","Detector")]
  dummy$Id <- factor(dummy$Id, IDs)
  dummy$Year <- factor(dummy$Year, Years)
  dummy$Detector <- factor(dummy$Detector, Detectors)
  temp <- table(dummy$Id, dummy$Detector, dummy$Year)
  
  ##-- CREATE THE DETECTION ARRAY
  id_names <- row.names(temp)
  Count.ar <- array(NA, c(dim(temp)))
  Count.ar[] <- temp#-1
  dimnames(Count.ar) <- dimnames(temp)
  if(method =="Bernoulli"){Count.ar[Count.ar>=1] <- 1}
  
  ##-- REORDER BY INDIVIDUALS TO MATCH THE PREVIOUS VERSION OF MAKE Y 
  Count.ar <- Count.ar[order(dimnames(Count.ar)[[1]]), , ]
  id_names <- id_names[order(id_names)]
  
  ##-- OUTPUT SETTINGS
  if(returnIdvector){
    output <- list( y.ar = Count.ar,
                    Id.vector = id_names)
  } else {
    output <- list(y.ar = Count.ar)
  }
  
  ##-- IF SECOND DATASET (e.g. DEAD RECOVERIES)
  if(!is.null(data2)){
    ##-- AGGREGATE MULTIPLE DETECTIONS AT THE SAME SUBDETECTOR TO ONE DETECTIONS ONLY
    if(method2 == "Binomial"){
      data2 <- data2[!duplicated(data2[ ,c("sub.detector","Year", "Id" )]), ] 
      data2 <- droplevels(data2)
    }
    
    ##-- NUMBER OF DETECTORS 
    N_Detec2 <- max(unlist(lapply(detectors2, function(x) length(x))))
    Detectors2 <- 1:N_Detec2
    
    ##-- CREATE A DUMMY DATASET
    dummy2 <- data2[ ,c("Id","Year","Detector")]
    dummy2$Id <- factor(dummy2$Id,IDs)
    dummy2$Year <- factor(dummy2$Year,Years)
    dummy2$Detector <- factor(dummy2$Detector,Detectors)
    temp2 <- table(dummy2$Id, dummy2$Detector, dummy2$Year)
    
    Count.ar2 <- array(NA,c(dim(temp2)))
    Count.ar2[] <- temp2
    dimnames(Count.ar2) <- dimnames(temp2)
    if(method=="Bernoulli"){Count.ar2[Count.ar2 >= 1] <- 1}
    
    
    ##-- OUTPUT SETTINGS
    Count.ar <- as.array(Count.ar)
    Count.ar2 <- as.array(Count.ar2)
    
    ##-- reorder to match previous MAKE Y 
    Count.ar <- Count.ar[order(dimnames(Count.ar)[[1]]), , ]
    Count.ar2 <- Count.ar2[order(dimnames(Count.ar2)[[1]]), , ]
    id_names <- id_names[order(id_names)]
    
    
    ##-- RETURN OUTPUT
    if(returnIdvector){
      output <- list( y.ar = as.array(Count.ar),
                      y.ar2 = as.array(Count.ar2),
                      Id.vector = id_names)
    } else {
      output <- list( y.ar = Count.ar,
                      y.ar2 = Count.ar2)
    }
  }#if
  return(output)
}

