#' @title Check individual distances between all detections of an individual and the centroid of its detection
#'
#' @description
#' \code{CheckDistanceDetections} returns a \code{list} object with the distances between detections and centroid of detections and which detection from which id is outside the max.distance if specified 
#'
#' @param y A \code{matrix} or \code{array} of detections 
#'	@param detector.xy A \code{matrix} object with xy coordinates of detectors. 
#' @param max.distance A \code{numeric} distance (radius) used to identify  which detections are too far
#' @param plot.check A \code{logical} plot distance histogram if TRUE
#' @return  a \code{list} object with the distances between detections and centroid of detections and which detection from which id is outside the max.distance if specified 
#' 

CheckDistanceDetections <- function(  
    y = y,
    detector.xy = detector.xy,
    max.distance = NULL,
    method = "pairwise",
    plot.check = TRUE){
  
  ##-- MAKE CONTAINER FOR FLAGS
  y.flagged <- y.distance <- y
  y.flagged[] <- 0
  y.distance[] <- NA
  
  ##-- SET THE DECTECTOR COORDINATES
  if(is.matrix(detector.xy)){detector.xy <- as.data.frame(detector.xy)}
  
  ##-- Convert to lower capital
  colnames(detector.xy) <- tolower(colnames(detector.xy))
  detector.sp <-  st_as_sf(detector.xy, coords = c("x", "y"))
  
  
  detector.index <- apply(y, 1, function(x) detector.sp[which(x>0),])
  detections.sp.ls <- apply(y, 1 ,function(x) detector.sp[which(x>0),])
  
  
  if(method == "pairwise"){
    for(i in 1:dim(y)[1]){
      detector.index <- which(y[i,]>0)
      this.det <- detector.sp[detector.index,]
      dist.det <- st_distance(this.det, this.det, byid = T)
      n.too.far <- apply(dist.det,2,function(x){
        sum(x>= max.distance)
      })
      
      to.remove <- NULL  
      
      if(any(n.too.far > 0)){
        if(dim(dist.det)[1] == 2){
          to.remove <- 2 ##-- remove detections at the second detector if only 2 detectors and they are too far from each other
        }
        
        if(dim(dist.det)[1] >2 ){
          combos <- lapply(1:(nrow(this.det)-1),function(x){
            combn(1:nrow(this.det),x)
          })
          
          temp3 <- lapply(1:length(combos), function(x) {
            xx <- combos[[x]][,1]
            out <- apply(combos[[x]], 2, function(xx){
              this.det2 <- this.det[-xx,]
              
              if(!any(st_distance(this.det2, this.det2, byid = T)>=max.distance)){
                return(xx)
              }
            })
            try(out<-out[,1],silent=TRUE)
            return(out)
          })
          
          to.remove <- temp3[!unlist(lapply(temp3, is.null))][[1]]
          if(is.list(to.remove)){
            to.remove <- to.remove[!unlist(lapply(to.remove, is.null))][[1]]
          }
        }
        
      }
      if(plot.check){
        if(!is.null(to.remove)){
          ext <- data.frame(t(apply(st_coordinates(this.det), 2, min)))
          ext[2,2] <- ext[1,2] + max.distance
          ext[2,1] <- ext[1,1]
          
          ext <-  st_as_sf(ext, coords = c("X", "Y"))
          
          plot(st_geometry(this.det), main=i, pch=19, col="grey", cex=1)
          lines(st_coordinates(ext), main=i, col="green")
          plot(this.det[to.remove,], add=TRUE, pch=19, col="red", cex=0.8)
          box()
        }
      }
      
      y.flagged[i,detector.index[to.remove]] <- 1
      y.distance[i,detector.index[to.remove]] <- unlist(lapply(to.remove, function(x) max(dist.det[x,])))
    }
  }
  
  
  print(paste("Detections removed: ", sum(y * y.flagged), " of ",sum(y),sep=""))
  print(paste("Individuals affected: ", sum(apply(y.flagged, 1, function(x) sum(x)>0)), " of ", sum(apply(y, 1, function(x)sum(x)>0)), sep=""))
  
  out <- list(y.flagged=y.flagged, y.distance=y.distance) 
  return(out)
  
}
