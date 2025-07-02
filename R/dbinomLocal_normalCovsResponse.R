#' Local evaluation of a binomial SCR detection process 
#'
#' The \code{dbinomLocal_normalCovsResponse} distribution is a NIMBLE custom 
#' distribution which can be used to model and simulate binomial observations 
#' (\emph{x}) of a single individual over a set of detectors defined by their 
#' coordinates (\emph{trapCoords}). The distribution assumes that an individual’s 
#' detection probability at any detector follows a half-normal function of the 
#' distance between the individual's activity center (\emph{s}) and the detector location. 
#' All coordinates (\emph{s} and \emph{trapCoords}) should be scaled to the habitat 
#' (see \code{\link[nimbleSCR]{scaleCoordsToHabitatGrid}})
#'
#' The \code{dbinomLocal_normalCovsResponse} distribution incorporates three 
#' features to increase computation efficiency (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item A local evaluation of the detection probability calculation (see Milleret 
#' et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item A sparse matrix representation (\emph{x}, \emph{detIndices} and \emph{detNums})
#' of the observation data to reduce the size of objects to be processed.
#' \item An indicator (\emph{indicator}) to shortcut calculations for individuals
#' unavailable for detection.
#' }
#' 
#' The \code{dbinomLocal_normalCovsResponse} distribution requires x- and y- 
#' detector coordinates (\emph{trapCoords}) and activity centers coordinates 
#' (\emph{s}) to be scaled to the habitat grid (\emph{habitatGrid}) using the 
#' (\code{\link[nimbleSCR]{scaleCoordsToHabitatGrid}} function.)
#'
#' When the aim is to simulate detection data: 
#' \enumerate{
#' \item \emph{x} should be provided using the \emph{yCombined} object as returned 
#' by \code{\link[nimbleSCR]{getSparseY}}, 
#' \item arguments \emph{detIndices} and \emph{detNums} should not be provided, 
#' \item argument \emph{lengthYCombined} should be provided using the \emph{lengthYCombined} 
#' object as returned by \code{\link[nimbleSCR]{getSparseY}}.
#' }
#' 
#' 
#' 
#' @name dbinomLocal_normalCovsResponse
#'
#' @param x Vector of individual detection frequencies. This argument can be provided in two formats: (i) with the \emph{y} object as returned by the \code{\link[nimbleSCR]{getSparseY}} function; (ii) with the \emph{yCombined} object as returned by \code{\link[nimbleSCR]{getSparseY}}. 
#' Note that when the random generation functionality is used (\code{rbinomLocal_normal}), only the \emph{yCombined} format can be used. 
#' The \emph{yCombined} object combines \emph{detNums}, \emph{x}, and \emph{detIndices} (in that order).  When such consolidated representation of the detection data \emph{x} is used, \emph{detIndices} and \emph{detNums} arguments shouldn’t be specified.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param detNums Number of detections recorded in \emph{x}, as returned by the \emph{detNums} object from the \code{\link[nimbleSCR]{getSparseY}} function. This argument should not be specified when the \emph{yCombined} object (returned by \code{\link[nimbleSCR]{getSparseY}}) is provided as \emph{x}, and when detection data are simulated.
#' @param detIndices Vector of indices of traps where the detections in x were recorded, as returned by the \emph{detIndices} object from the \code{\link[nimbleSCR]{getSparseY}} function. This argument should not be specified when \emph{x} is provided as the \emph{yCombined} object (returned by \code{\link[nimbleSCR]{getSparseY}}) and when detection data are simulated.
#' @param size Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).
#' @param p0State Vector of baseline detection probabilities for each trap used in the half-normal detection function. When \emph{p0Traps} is used, \emph{p0} should not be provided. 
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Individual activity center x- and y-coordinates scaled to the habitat (see \code{\link[nimbleSCR]{scaleCoordsToHabitatGrid}}).
#' @param trapCoords Matrix of x- and y-coordinates of all traps scaled to the habitat (see \code{\link[nimbleSCR]{scaleCoordsToHabitatGrid}}).
#' @param localTrapsIndices Matrix of indices of local traps around each habitat grid cell, as returned by the \code{\link[nimbleSCR]{getLocalObjects}} function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the \code{\link[nimbleSCR]{getLocalObjects}} function.
#' @param resizeFactor Aggregation factor used in the \code{\link[nimbleSCR]{getLocalObjects}} function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid Matrix of habitat grid cells indices, as returned by the \code{\link[nimbleSCR]{getLocalObjects}} function.
#' @param indicator Logical argument specifying whether the individual is available for detection.
#' @param lengthYCombined The length of the  x argument when the (\emph{yCombined}) format of the detection data is provided (as returned by the \emph{lengthYCombined} object from \code{\link[nimbleSCR]{getSparseY}}). 
#' @param trapCountries To allow the possibility that some habitat grid cells have no local traps in the surroundings (default to FALSE). 
#' @param trapCovs To allow the possibility that some habitat grid cells have no local traps in the surroundings (default to FALSE). 
#' @param trapBetas To allow the possibility that some habitat grid cells have no local traps in the surroundings (default to FALSE). 
#' @param responseCovs To allow the possibility that some habitat grid cells have no local traps in the surroundings (default to FALSE). 
#' @param responseBetas To allow the possibility that some habitat grid cells have no local traps in the surroundings (default to FALSE). 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / \sigma^2)}.
#'
#' @author Pierre Dupont, Cyril Milleret 
#'
#' @import nimble 
#' @importFrom stats dbinom rbinom
#' 
#' @export
NULL

#' @rdname dbinomLocal_normalCovsResponse
#' @export
dbinomLocal_normalCovsResponse <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0State = double(1),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),                    
                  indicator = double(0, default = 1.0),
                  lengthYCombined = double(0, default = 0),
                  trapCountries = double(1),
                  trapCovs = double(2),
                  trapBetas = double(1),
                  responseCovs = double(0),
                  responseBetas = double(0),
                  log = integer(0, default = 0)
  ){
    
    ##-- Return type declaration
    returnType(double(0))
    
    ##-- Deal with cases where detection info is combined in one vector 
    if(detNums == -999){
      detNums <- x[1]
      nMaxDetectors <- (lengthYCombined-1)/2
      detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
      x1 <- x[2:(nMaxDetectors+1)]
    } else {
      x1 <- x
      detIndices1 <- detIndices
    }
    
    ##-- Shortcut if the current individual is not available for detection
    if(indicator == 0){
      if(detNums == 0){
        if(log == 0) return(1.0)
        else return(0.0)
      } else {
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
    
    ##-- Retrieve the index of the habitat cell where the current AC is
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ##-- Retrieve the indices of the local traps surrounding the selected habitat grid cell
    theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
    
    ##-- Check if detections are within the list of local traps
    if(detNums > 0){
      for(r in 1:detNums){
        if(sum(detIndices1[r] == theseLocalTraps) == 0){
          if(log == 0) return(0.0)
          else return(-Inf)
        }
      }#r
    }
    
    ##-- Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logitp0 <- logit(p0State)
    logProb <- 0.0 
    detIndices1 <- c(detIndices1, 0)
    indCov <- responseCovs * responseBetas
    count <- 1 
    
    for(r in 1:localTrapsNum[sID]){
      ##-- Calculate distance to local traps
      d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) +
        pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
      ##-- Calculate p0 
      pZero <- ilogit(logitp0[trapCountries[theseLocalTraps[r]]] + 
                        inprod(trapBetas,trapCovs[theseLocalTraps[r], ]) +
                        indCov)
      p <- pZero * exp(alpha * d2)
      
      if(theseLocalTraps[r] == detIndices1[count]){ 
        logProb <- logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
        count <- count + 1
      } else {
        logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
      }
    }#r
    
    if(log)return(logProb)
    return(exp(logProb))
  })


NULL
#' @rdname dbinomLocal_normalCovsResponse
#' @export
rbinomLocal_normalCovsResponse <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0State = double(1),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),                    
                  indicator = double(0, default = 1.0),
                  lengthYCombined = double(0, default = 0),
                  trapCountries = double(1),
                  trapCovs = double(2),
                  trapBetas = double(1),
                  responseCovs = double(0),
                  responseBetas = double(0)
  ) {
    ##-- Specify return type
    returnType(double(1))
    
    ##-- Initial checks
    if(detNums >= 0) stop("Random generation for the dbinomLocal_normalCovsResponse distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    if(n!=1){print("dbinomLocal_normalCovsResponse only allows n = 1; using n = 1")}
    
    ##-- Get necessary info
    alpha <- -1.0 / (2.0 * sigma * sigma)
    nMAxDetections <- (lengthYCombined-1)/2
    
    ##-- SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
    if(indicator == 0){return(rep(0.0, lengthYCombined))}
    
    ##-- RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ##-- RETRIEVE THE IDs OF THE RELEVANT DETECTORS
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
    
    ##-- INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
    detectOut <- rep(0, localTrapsNum[sID])
    ys <- rep(-1, nMAxDetections)
    dets <- rep(-1, nMAxDetections)
    logitp0 <- logit(p0State)
    indCov <- responseCovs * responseBetas
    count <- 1
    
    ##-- SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
    for(r in 1:localTrapsNum[sID]){
      ##-- Calculate distance to local traps
      d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) +
        pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
      
      ## Calculate p0 
      pZero <- ilogit(logitp0[trapCountries[theseLocalTraps[r]]] + 
                        inprod(trapBetas,trapCovs[theseLocalTraps[r], ]) +
                        indCov)
      p <- pZero * exp(alpha * d2)
      
      ##-- Draw the observation at detector j from a binomial distribution with probability p
      detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], p)
      if(detectOut[r] > 0){
        if(nMAxDetections < count){
          stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                          You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")
        }
        ys[count] <- detectOut[r]
        dets[count] <- theseLocalTraps[r]
        count <- count + 1
      }#if
    }#r 
    
    ##-- Format output vector
    count <- count - 1
    out <- rep(-1, lengthYCombined)
    out[1] <- count
    if(count >= 1){
      out[2:(count+1)] <- ys[1:count]
      out[(nMAxDetections+2):(nMAxDetections+count+1)] <- dets[1:count]
    }
    
    ##-- OUTPUT
    return(out)
  })
