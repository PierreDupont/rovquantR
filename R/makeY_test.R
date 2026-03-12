#' Detection history formatting
#'
#' R utility function to format SCR and OPSCR detection data into detection history arrays of different formats.
#' The format of the detection history is determined by the name of the \code{nimble} distribution used in the model (provided by the 'detection.fun' argument)
#'
#' The \code{makeDetectionHistory} function is used in the data2nimbleSCR function to help
#' formatting SCR data for use with nimbleSCR.
#'
#' @param data a data.frame (or list of data.frames if the model integrates multiple datasets) containing individual detections (in rows).
#' Mandatory column names are "id" and "detector",
#' Optional column names are "x", "y", "session" and "occasion".
#' @param detectors a data.frame (or list of data.frames if the model integrates multiple datasets) containing the coordinates of all possible detection locations.
#' @param detection.fun a character string indicating the nimbleSCR (or nimble) function used to model individual detections.
#' Must be one of :
#' \itemize{
#' \item "dbern", "dbinom", "dpois": Bernoulli, Binomial or Poisson distributions to model binary detections or detection frequencies on detection location at a time.
#' \item "dbinomLocal_normal", "dpoisLocal_normal": Local-evaluation functions; used to model a vector of individual detections at multiple discrete detectors. 
#' \item "dbernppDetection_normal", "dpoisppDetection_normal": Bernoulli and Poisson point-process functions; used to model SCR detections in continuous space.
#' \item "dbernppLocalDetection_normal", "dpoisppLocalDetection_normal": local evaluation versions of the two funcctions above for faster implementation.
#' }
#' @param all.ids vector of individual identifiers to be included in the detection history array. Can be used to either 
#' @param diff.occ An aggregation factor to reduce the number of habitat cells to 
#' retrieve local objects for. Defaults to 1; no aggregation.
#' @param M Augmented population size. 
#'
#' @return This function returns a detection history array in one of the following 
#' format:
#' \itemize{
#' \item y[]: ...
#' \item y[]: ...
#' \item y[]: ...
#' \item y[]: ...
#'}
#' 
#' @author Pierre Dupont
#' 
#' @import nimble
#' @importFrom stats dnorm 
#' 
#' @references
#' nimbleSCR manual: ... 
#'
#' @examples 
#' n.samples <- 150
#' n.ids <- 40
#' n.detectors <- 10
#' n.sessions <- 3
#' n.occasions <- 5
#' 
#' detectors <- expand.grid(list("id" = 1:n.detectors,
#'                               #"occasion" = 1:n.occasions,
#'                               "session"=1:n.sessions))
#' detectors <- cbind.data.frame("id" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                               "session" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                               "occasion" = sample(n.occasions,size = n.samples,replace = TRUE))
#' head(detectors)
#' 
#' data <- cbind.data.frame("id" = sample(n.ids,size = n.samples,replace = TRUE),
#'                          "detector" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                          "session" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                          "occasion" = sample(n.occasions,size = n.samples,replace = TRUE))
#' head(data)
#' 
#' all.ids <- 1:n.ids
#'
#' makeDetectionHistory(data, detectors, "dbern", all.ids)
#' makeDetectionHistory(data, detectors, "dbinom", all.ids)
#' makeDetectionHistory(data, detectors, "dbinom_vector", all.ids)
#' makeDetectionHistory(data, detectors, "dpois", all.ids)
#' makeDetectionHistory(data, detectors, "dcat", all.ids)
#' makeDetectionHistory(data, detectors, "dbinomLocal_normal", all.ids)
#' makeDetectionHistory(data, detectors, "dpoisLocal_normal", all.ids)
#' makeDetectionHistory(data, detectors, "dbernppLocalDetection_normal", all.ids)
#' makeDetectionHistory(data, detectors, "dpoisppLocalDetection_normal", all.ids)


NULL
#' @rdname makeDetectionHistory
#' @export
makeDetectionHistory <- function( data,   
                                  detectors,
                                  detection.fun,
                                  all.ids = NULL,
                                  diff.occ = FALSE,
                                  M = NULL){
  ##---- 1. Initial checks ----
  
  ##-- 1. Make sure requested detection functions are available in nimbleSCR
  availFunctions <- c("dcat", "dbern", "dbinom", "dpois",
                      "dbinom_vector", 
                      "dbinomLocal_normal", "dbinomLocal_exp", "dbinomLocal_normalPlateau",
                      "dpoisLocal_normal",
                      "dbernppDetection_normal", "dpoisppDetection_normal",
                      "dbernppLocalDetection_normal", "dpoisppLocalDetection_normal")
  if(!all(detection.fun %in% availFunctions)) {
    stop(paste0(detection.fun, " is not a detection function available in the nimbleSCR package!"))
  }
  
  
  ##-- 2. Sanitize `data` and `detectors` 
  sanitizeInput(data, detectors)
  
  
  ##-- 3.Check if the number of `detection.fun` matches the number of `data` 
  if(length(detection.fun) != length(data)){
    if(length(detection.fun) == 1){
      warning("There is only one `detection.fun` provided for multiple `data`.\n Recycling the same `detection.fun` for all `data`.")
      detection.fun <- rep(detection.fun, length(data))
    } else {
      stop("The number of `detection.fun` is different from the number of `data` provided.")
    }
  }
  
  
  ##-- 4. Get the total number of ids 
  if(is.null(all.ids)){
    all.ids <- sort(unique(unlist(lapply(data, function(x)x[ ,"id"]))))
  }
  if(!is.null(M)){
    if(M <= length(all.ids)){ 
      warnings("M is lower than the number of individuals detected!")
    }
    n.aug <- M - length(all.ids)
    all.ids <- as.factor(c(as.character(all.ids),
                           paste0("augmented.",1:n.aug)))
  }
  n.ids <- length(all.ids)
  
  
  
  ##---- 2. Create detection histories ----
  
  for(l in 1:length(data)){
    
    thisData <- data[[l]]
    thisDetectors <- detectors[[l]]
    
    ##-- Tabulate # of secondary occasions per detector per primary session
    tab <- as.array(table(thisDetectors$session,
                          thisDetectors$detector,
                          thisDetectors$occasion))
    
    ##-- Get number of primary sessions
    all.sess <- dimnames(tab)[[1]]
    n.sess <- length(all.sess)
    
    ##-- Get number of detectors per primary session
    all.dets <- dimnames(tab)[[2]]
    n.dets <- rowSums(apply(tab, c(1,2), function(x)any(x > 0)))
    n.dets.max <- max(n.dets)
    
    ##-- Get number of secondary occasions per detector and primary session
    all.occ <- dimnames(tab)[[3]]
    n.occ <- apply(tab, c(1,2), function(x)sum(x > 0))
    n.occ.max <- max(n.occ)

    
    
    ##----    2.1. "Classic" formats ----
    
    if(detection.fun %in% c("dbern","dpois")){
      if(detection.fun %in% c("dbern")){
        ##-- If Bernoulli-type detection function, keep only one detection per secondary
        ##-- occasion (e.g. weekly sampling) 
        thisData <- thisData[!duplicated(thisData[ ,c("detector","session","occasion","id")]), ] 
        thisData <- droplevels(thisData)
      }
      
      ##-- Tabulate individual frequencies per detector, occasion and session
      thisData$id <- factor(thisData$id, all.ids)
      thisData$detector <- factor(thisData$detector, all.dets)
      thisData$occasion <- factor(thisData$occasion, all.occ)
      thisData$session <- factor(thisData$session, all.sess)
      y <- as.array(table(thisData$id, thisData$detector, thisData$occasion, thisData$session))
      dimnames(y) <- list("id" = all.ids,
                          "detector" = all.dets,
                          "occasion" = all.occ,
                          "session" = all.sess)
    }
    
    if(detection.fun %in% c("dbinom","dbinom_vector")){
      if("subdetector" %in% names(thisData)){
        ##-- If PAB-type detection function, keep only one detection per subdetector
        thisData <- thisData[!duplicated(thisData[ ,c("detector","session","occasion","subdetector","id")]), ] 
        thisData <- droplevels(thisData) 
        
        ##-- Tabulate individual frequencies per detector, occasion and session
        thisData$id <- factor(thisData$id, all.ids)
        thisData$detector <- factor(thisData$detector, all.dets)
        thisData$occasion <- factor(thisData$occasion, all.occ)
        thisData$session <- factor(thisData$session, all.sess)
        y <- as.array(table(thisData$id, thisData$detector, thisData$occasion, thisData$session))
        dimnames(y) <- list("id" = all.ids,
                            "detector" = all.dets,
                            "occasion" = all.occ,
                            "session" = all.sess)
      } else {
        ##-- If binomial-type detection function, keep only one detection per secondary
        ##-- occasion (e.g. weekly sampling) 
        thisData <- thisData[!duplicated(thisData[ ,c("detector","session","occasion","id")]), ] 
        thisData <- droplevels(thisData)
        
        ##-- Tabulate individual frequencies per detector and session
        thisData$id <- factor(thisData$id, all.ids)
        thisData$detector <- factor(thisData$detector, all.dets)
        thisData$occasion <- factor(thisData$occasion, all.occ)
        thisData$session <- factor(thisData$session, all.sess)
        y <- as.array(table(thisData$id, thisData$detector, thisData$session))
        dimnames(y) <- list("id" = all.ids,
                            "detector" = all.dets,
                            "session" = all.sess)
      }
    }
    
    if(detection.fun %in% c("dcat")){
      ##-- Tabulate # of detections per individual, occasion and session
      detNums <- as.array(table(thisData$id, thisData$occasion, thisData$session))
      if(any(detNums > 1)){
        stop("Individuals can be detected at only one detector per sampling occasion when using `dcat`")
      }
      detI <- dimnames(detNums)[[1]]
      detK <- dimnames(detNums)[[2]]
      detT <- dimnames(detNums)[[3]]
      
      ##-- Set-up array to store individual detector indices
      ##-- ´max(n.dets) + 1´ is used for individuals not detected
      y <- array(max(n.dets) + 1, c(n.ids, n.occ.max, n.sess))
      dimnames(y) <- list("id" = all.ids,
                          "session" = all.sess,
                          "occasion" = all.occ)
      
      ##-- Fill in the array
      detIndex <- which(detNums == 1, arr.ind = T)
      for(d in 1:dim(detIndex)[1]){
        thisI <- which(all.ids == detI[d,1])
        thisK <- which(all.occ == detK[d,2])
        thisT <- which(all.sess == detT[d,3])
        y[thisI,thisK,thisT] <- thisData[thisData$id == detI[d,1] &
                                       thisData$occasion == detK[d,2] & 
                                       thisData$session == detT[d,3], "detector"]
        
      }#d
    }
    
    
    
    ##----    2.2. "Sparse" formats ----
    
    if(detection.fun %in% c("dpoisLocal_normal", "dpoisLocal_exp")){
      ##-- Tabulate # of individual detections per detector, occasion and session
      detMat <- as.array(table(thisData$id, thisData$detector, thisData$occasion, thisData$session))
      
      ##-- Get # of individuals detected
      n.detected <- dim(detMat)[1]
      
      ##-- Get # of detectors at which each individual is detected
      detNums <- apply(detMat, c(1,3,4), function(x) sum(x>0))
      
      ##-- Set-up array to store individual detection frequencies and indices
      detIndices <- detFreq <- array(-1, c(n.ids, max(detNums), n.occ.max, n.sess))
      y <- array(0, c(n.ids, 2*max(detNums)+1, n.occ.max, n.sess)) 
      dimnames(y) <- list("id" = all.ids,
                          "detections" = c("detNums",
                                           rep("detFreq", max(detNums)),
                                           rep("detIndex", max(detNums))),
                          "occasion" = all.occ,
                          "session" = all.sess)    
      
      ##-- Fill in the array
      for(i in 1:n.detected){
        thisID <- which(all.ids == dimnames(detMat)[[1]][i])
        for(k in 1:n.occ.max){
          for(t in 1:n.sess){
            if(detNums[i,k,t] > 0){
              detIndices[thisID,1:detNums[i,k,t],k,t] <- as.numeric(names(which(detMat[i, ,k,t] > 0)))
              detFreq[thisID,1:detNums[i,k,t],k,t] <- detMat[i,which(detMat[i, ,k,t] > 0),k,t]
              y[thisID, ,k,t] <- c(detNums[i,k,t],
                                   detFreq[thisID, ,k,t],
                                   detIndices[thisID, ,k,t])
            }
          }#t
        }#k
      }#i
    }
    
    if(detection.fun %in% c("dbinomLocal_normal", "dbinomLocal_exp",
                            "dbinomLocal_normalPlateau", "dmultiLocal_normal")){
      if("subdetector" %in% names(thisData)){ 
        ##-- If PAB-type detection function, keep only one detection per subdetector
        thisData <- thisData[!duplicated(thisData[ ,c("detector","session","occasion","subdetector","id")]), ] 
        thisData <- droplevels(thisData) 
        
        ##-- Tabulate # of individual detections per detector, occasion and session
        detMat <- as.array(table(thisData$id, thisData$detector, thisData$occasion, thisData$session))
        
        ##-- Get # of individuals detected
        n.detected <- dim(detMat)[1]
        
        ##-- Get # of detectors at which each individual is detected
        detNums <- apply(detMat, c(1,3,4), function(x) sum(x>0))
        
        ##-- Set-up array to store individual detection frequencies and indices
        detIndices <- detFreq <- array(-1, c(n.ids, max(detNums), n.occ.max, n.sess))
        y <- array(0, c(n.ids, 2*max(detNums)+1, n.occ.max, n.sess)) 
        dimnames(y) <- list("id" = all.ids,
                            "detections" = c("detNums",
                                             rep("detFreq", max(detNums)),
                                             rep("detIndex", max(detNums))),
                            "occasion" = all.occ,
                            "session" = all.sess)    
        
        ##-- Fill in the array
        for(i in 1:n.detected){
          thisID <- which(all.ids == dimnames(detMat)[[1]][i])
          for(k in 1:n.occ.max){
            for(t in 1:n.sess){
              if(detNums[i,k,t] > 0){
                detIndices[thisID,1:detNums[i,k,t],k,t] <- as.numeric(names(which(detMat[i, ,k,t] > 0)))
                detFreq[thisID,1:detNums[i,k,t],k,t] <- detMat[i,which(detMat[i, ,k,t] > 0),k,t]
                y[thisID, ,k,t] <- c(detNums[i,k,t],
                                     detFreq[thisID, ,k,t],
                                     detIndices[thisID, ,k,t])
              }
            }#t
          }#k
        }#i
      } else { 
        ##-- If binomial-type detection function, keep only one detection per secondary
        ##-- occasion (e.g. weekly sampling) 
        thisData <- thisData[!duplicated(thisData[ ,c("detector","session","occasion","id")]), ] 
        thisData <- droplevels(thisData) 
        
        if(diff.occ){
          ##-- Tabulate # of individual detections per detector and session
          detMat <- as.array(table(thisData$id,thisData$detector,thisData$occasion,thisData$session))
          
          ##-- Get # of individuals detected
          n.detected <- dim(detMat)[1]
          
          ##-- Get # of detectors at which each individual is detected
          detNums <- apply(detMat, c(1,3,4), function(x) sum(x>0))
          
          ##-- Set-up array to store individual detection frequencies and indices
          detIndices <- detFreq <- array(-1, c(n.ids, max(detNums), n.occ.max, n.sess))
          y <- array(0, c(n.ids, 2*max(detNums) + 1, n.occ.max, n.sess)) 
          dimnames(y) <- list("id" = all.ids,
                              "detections" = c("detNums",
                                               rep("detFreq", max(detNums)),
                                               rep("detIndex", max(detNums))),
                              "occasion" = all.occ,
                              "session" = all.sess)    
          
          ##-- Fill in the array
          for(i in 1:n.detected){
            thisID <- which(all.ids == dimnames(detMat)[[1]][i])
            for(t in 1:n.sess){
              for(o in 1:n.occ.max){
                if(detNums[i,o,t] > 0){
                  detIndices[thisID,1:detNums[i,o,t],o,t] <- as.numeric(names(which(detMat[i, ,o,t] > 0)))
                  detFreq[thisID,1:detNums[i,o,t],o,t] <- detMat[i,which(detMat[i, ,o,t] > 0),o,t]
                  y[thisID, ,o,t] <- c(detNums[i,o,t],
                                       detFreq[thisID, ,o,t],
                                       detIndices[thisID, ,o,t])
                } 
              }#o
            }#t
          }#i
        } else { 
          ##-- Tabulate # of individual detections per detector and session
          detMat <- as.array(table(thisData$id, thisData$detector, thisData$session))
          
          ##-- Get # of individuals detected
          n.detected <- dim(detMat)[1]
          
          ##-- Get # of detectors at which each individual is detected
          detNums <- apply(detMat, c(1,3), function(x) sum(x>0))
          
          ##-- Set-up array to store individual detection frequencies and indices
          detIndices <- detFreq <- array(-1, c(n.ids, max(detNums), n.sess))
          y <- array(0, c(n.ids, 2*max(detNums) + 1, n.sess)) 
          dimnames(y) <- list("id" = all.ids,
                              "detections" = c("detNums",
                                               rep("detFreq", max(detNums)),
                                               rep("detIndex", max(detNums))),
                              "session" = all.sess)    
          
          ##-- Fill in the array
          for(i in 1:n.detected){
            thisID <- which(all.ids == dimnames(detMat)[[1]][i])
            for(t in 1:n.sess){
              if(detNums[i,t] > 0){
                detIndices[thisID,1:detNums[i,t],t] <- as.numeric(names(which(detMat[i, ,t] > 0)))
                detFreq[thisID,1:detNums[i,t],t] <- detMat[i,which(detMat[i, ,t] > 0),t]
                y[thisID, ,t] <- c(detNums[i,t],
                                   detFreq[thisID, ,t],
                                   detIndices[thisID, ,t])
              }
            }#t
          }#i 
        }
      }
    }
    
    
    
    ##----    2.3. "Point-process" formats ----
    
    if(detection.fun %in% c("dbernppLocalDetection_normal", "dbernppDetection_normal",
                            "dpoisppLocalDetection_normal", "dpoisppDetection_normal")){
      ##-- Check that detection "x" and "y" coordinates are provided
      if(!all(c("x","y") %in% names(thisData))){
        stop("`x` and `y` coordinates of detections must be provided in `data` when 
      using `dbernppLocalDetection_normal`, `dbernppDetection_normal`, 
           `dpoisppLocalDetection_normal` or `dpoisppDetection_normal`.")
      }
      
      ##-- Tabulate # of detections per individual, occasion and session
      detNums <- as.array(table(thisData$id, thisData$occasion, thisData$session))
      if(detection.fun %in% c("dbernppLocalDetection_normal", "dbernppDetection_normal")){
        if(any(detNums>1)){
          stop("Individuals can be detected only once per sampling occasion when using `dbernppLocalDetection_normal` or `dbernppDetection_normal`")
        }
      } 
      
      ##-- Get # of individuals detected
      detIDs <- dimnames(detNums)[[1]]
      
      y <- array(0, c(n.ids, 3, max(detNums), n.occ.max, n.sess))
      dimnames(y) <- list("id" = all.ids,
                          "coords" = c("x","y","index"),
                          "detections" = 1:max(detNums),
                          "occasion" = all.occ,
                          "session" = all.sess)
      
      ##-- Fill in the array
      for(i in 1:length(detIDs)){
        thisID <- which(all.ids == detIDs[i])
        for(k in 1:n.occ.max){
          for(t in 1:n.sess){
            if(detNums[i,k,t] > 0){
              y[thisID,1:3,1:detNums[i,k,t],k,t] <- unlist(thisData[thisData$id == detIDs[i],
                                                                c("x","y","detector")])
            }
          }#t
        }#k
      }#i
    }
    
    
    
    ##----    2.4. Drop redundant dimensions & return y ----
    
    out[[l]] <- drop(y)
    
  }#l
  
  
  ##-- Output 
  return(out)
  
  # stop("You've reached an area under construction...\n (please come back soon!)")
}
