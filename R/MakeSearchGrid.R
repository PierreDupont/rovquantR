#' @title Function to create a detector grid
#' 
#' @description
#' \code{makeSearchGrid} creates equally spaced detectors within a defined area (data) using 
#' a given resolution (i.e. distance between two detectors).
#' If specified, each detector grid can be divived into multiple equally spaced sub-detectors (e.g. for use in the PAB model).  
#' 
#' @param data A \code{sf} object. can be a \code{SpatialPolygons} or \code{SpatialPoints} or \code{SpatialLines} or \code{raster}. Spatial Lines objects takes time to compute
#' @param resolution Numeric variable denoting the size of grid cells in units of (\code{polygon}).
#' @param div Numeric variable denoting the number of equally-sized subplots grid cell is to be divided into.
#' @param center Logical variable denoting whether nodes of the resulting search grid are to be centered within their respective subplots.
#' @param plot Logical for whether (\code{TRUE}) or not (\code{FALSE}) check plots are to be generated.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no} and Pierre Dupont
#' 
#' @import sf 
#' @import raster
#' @importFrom dplyr group_by summarise
#' @importFrom graphics plot points
#' @importFrom stars st_as_stars
#' @importFrom methods as
#'    
#' @rdname makeSearchGrid 
#' @export
makeSearchGrid <- function(
    data,
    resolution,
    div = 1,
    center = T,
    plot = TRUE){
  
  ##-- if the data is already a raster   
  subdetector.r <- data

  ##-- Check that 'data' is a raster  
  if(!inherits(data,"Raster")){
    stop("Data must be a raster")
  }
  
  ##-- AGGREGATE TO MAIN DETECTORS
  if(sqrt(div) > 1){
    maindetector.r <- raster::aggregate( subdetector.r,
                                         fact = sqrt(div), 
                                         fun = sum)
  } else {
    maindetector.r <- subdetector.r
  }
  
  ##-- CREATE POLYGONS FROM RASTER
  ##-- Main detectors 
  temp.r <- raster::raster(maindetector.r)
  temp.r[] <- 1:length(temp.r)
  maindetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                    as_points = FALSE, merge = TRUE)
  
  
  ##-- Sub-detectors
  temp.r <- subdetector.r
  temp.r[] <- 1:length(temp.r)
  subdetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                   as_points = FALSE, merge = TRUE)
  
  
  ##-- OBTAIN SPATIALPOINTS FROM DETECTORS 
  ##-- Main detectors 
  main.detector.xy <- as.data.frame(raster::xyFromCell(maindetector.r,
                                                       1:raster::ncell(maindetector.r)))
  
  main.detector.sp <- sf::st_as_sf(main.detector.xy, coords = c("x", "y"))
  st_crs(main.detector.sp) <- sf::st_crs(data)
  main.detector.sp$main.cell.x <- main.detector.xy[,"x"]
  main.detector.sp$main.cell.y <- main.detector.xy[,"y"]
  main.detector.sp$main.cell.id <- 1:nrow(main.detector.sp)
  
  ##-- Sub-detectors 
  if(center){
    ##-- ALTERNATIVE 1: center points of subdetector detector cells
    detector.xy <- as.data.frame(raster::xyFromCell(subdetector.r,
                                                    1:raster::ncell(subdetector.r)))
    
    sub.detector.sp <-  sf::st_as_sf(detector.xy, coords = c("x", "y"))
    sf::st_crs(sub.detector.sp) <- sf::st_crs(data)
  } else {
    stop("Only center is available")
  }
  sub.detector.sp$Id <- 1:length(sub.detector.sp)
  
  ##-- SELECT ACTIVE DECTECTORS (SEARCHED) 
  sub.detector.sp$main.cell.id <- as.numeric(sf::st_intersects(sub.detector.sp, maindetector.poly))#over(sub.detector.sp, maindetector.poly)
  merge.df <- sf::st_join(sub.detector.sp, main.detector.sp, by = "main.cell.id")
  colnames(merge.df)[which(colnames(merge.df) %in% "main.cell.id.x" )] <- "main.cell.id"
  
  merge.df <- merge.df[order(merge.df$Id), ]
  sub.detector.sp <- merge.df
  
  ##-- Subset main detectors 
  values.main.r <- raster::extract(maindetector.r, main.detector.sp)
  values.main.r[is.na(values.main.r)] <- 0
  main.detector.sp <- main.detector.sp[values.main.r>0, ]
  main.detector.sp$count <- values.main.r[values.main.r>0]
  
  ##-- Subset sub detectors 
  values.sub.r <- raster::extract(subdetector.r, sub.detector.sp)
  values.sub.r[is.na(values.sub.r)] <- 0
  sub.detector.sp <- sub.detector.sp[values.sub.r>0, ]
  sub.detector.sp$count <- values.sub.r[values.sub.r>0]
  
  ##-- MAKE A PLOTTING FUNCTION 
  if(plot){
    plot(data, col= "gray20")
    plot(st_geometry(maindetector.poly), add=TRUE, lwd=3)
    plot(st_geometry(subdetector.poly), add=TRUE)
    plot(st_geometry(sub.detector.sp), pch=19, cex=0.6, col=as.numeric(sub.detector.sp$main.cell.id), add=TRUE)
    points(main.cell.y ~ main.cell.x, data=sub.detector.sp, pch=19, cex=1, col=as.numeric(sub.detector.sp$main.cell.id))
  }
  
  ##-- OUTPUT
  out <- list( detector.sp = sub.detector.sp,
               main.detector.sp = main.detector.sp,
               sub.grid.poly = subdetector.poly,
               grid.poly = maindetector.poly)
}


