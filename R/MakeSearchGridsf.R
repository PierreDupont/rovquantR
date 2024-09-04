#' @title Function to create detectors
#' #'
#' @description
#' \code{MakeSearchGrid} creates  detectors within a defined area (data) and resolution between detectors. 
#' Division precises if subspatial division should be performed (Following PAB model).  
#' 
#' @param data A \code{sp} object. can be a \code{SpatialPolygons} or \code{SpatialPoints} or \code{SpatialLines} or \code{raster}. Spatial Lines objects takes time to compute
#' @param resolution Numeric variable denoting the size of grid cells in units of (\code{polygon}).
#' @param div Numeric variable denoting the number of equally-sized subplots grid cell is to be divided into.
#' @param center Logical variable denoting whether nodes of the resulting search grid are to be centered within their respective subplots.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) check plots are to be generated.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeSearchGrid.R
#' @keywords simul

MakeSearchGridsf <- function(data = data,
                           resolution = resolution,
                           div = 1,
                           center = T,
                           plot = TRUE,
                           fasterize = FALSE){
   
   ## if the data is already a raster   
   subdetector.r <- data
   
   ## Obtain the resolution of the subdetectors 
   res1 <- resolution/sqrt(div)  
   
   ## Convert other data types to easier simpler format (no dataframe). 
   if(inherits(data,"SpatialPolygonsDataFrame")){data <- as(data,"SpatialPolygons")}
   if(inherits(data,"SpatialPointsDataFrame")){data <- as(data,"SpatialPoints")}
   if(inherits(data,"SpatialLinesDataFrame")){data <- as(data,"SpatialLines")}
   
   ### ==== CREATE THE SUBDETECTORS  ====
   ## if SPATIALLINES
   if(inherits(data,"SpatialLines")){
     stop("SpatialLines not available; data should be a raster")
   }
   
   ## if SPATIALPOINTS 
   if(inherits(data,"SpatialPoints")){
     stop("SpatialPoints not available; data should be a raster")
   }
   
   ## if SPATIALPOLYGONS
   if(inherits(data,"SpatialPolygons")){
     stop("SpatialPolygons not available; data should be a raster")
      detector.r1 <- raster(extent(data), resolution = res1, crs = proj4string(data))

   }
   
   ### ==== AGGREGATE TO MAIN DETECTORS  ====
   if(sqrt(div)>1){
      maindetector.r <- aggregate(subdetector.r, fact = sqrt(div), fun = sum)
   }else{
      maindetector.r <- subdetector.r
   }
   
   ### ==== CREATE POLYGONS FROM RASTER  ====
   ## Main detectors 
   temp.r <- raster(maindetector.r)
   temp.r[] <- 1:length(temp.r)
   maindetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                                   as_points = FALSE, merge = TRUE)
    
   
   ## Sub-detectors
   temp.r <- subdetector.r
   temp.r[] <- 1:length(temp.r)

   subdetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                            as_points = FALSE, merge = TRUE)
  
   
   ### ==== OBTAIN SPATIALPOINTS FROM DETECTORS ====
   ## Main detectors 
   main.detector.xy <- as.data.frame(xyFromCell(maindetector.r, 1:ncell(maindetector.r)))
   
   main.detector.sp <-  st_as_sf(main.detector.xy, coords = c("x", "y"))
   st_crs(main.detector.sp) <- st_crs(data)
   main.detector.sp$main.cell.x <- main.detector.xy[,"x"]
   main.detector.sp$main.cell.y <- main.detector.xy[,"y"]
   main.detector.sp$main.cell.id <- 1:nrow(main.detector.sp)
   
   ## Sub-detectors 
   if(center){
      ## ALTERNATIVE 1: center points of subdetector detector cells
      detector.xy <- as.data.frame(xyFromCell(subdetector.r, 1:ncell(subdetector.r)))
      
      sub.detector.sp <-  st_as_sf(detector.xy, coords = c("x", "y"))
      st_crs(sub.detector.sp) <- st_crs(data)
   }else{
     stop("Only center is available")
      
   }
   
   # col(sub.detector.sp) <- c("x","y")
   sub.detector.sp$Id <- 1:length(sub.detector.sp)
   
   ### ==== SELECT ACTIVE DECTECTORS (SEARCHED) ====
   sub.detector.sp$main.cell.id <- as.numeric(st_intersects(sub.detector.sp, maindetector.poly))#over(sub.detector.sp, maindetector.poly)
   merge.df <- st_join(sub.detector.sp, main.detector.sp, by="main.cell.id")
   colnames(merge.df)[which(colnames(merge.df) %in%"main.cell.id.x" )] <- "main.cell.id"
   
   merge.df <- merge.df[order(merge.df$Id), ]
   sub.detector.sp <- merge.df
   
   ## Subset main detectors 
   values.main.r <- raster::extract(maindetector.r, main.detector.sp)
   values.main.r[is.na(values.main.r)] <- 0
   main.detector.sp <- main.detector.sp[values.main.r>0, ]
   main.detector.sp$count <- values.main.r[values.main.r>0]
   
   ## Subset sub detectors 
   values.sub.r <- raster::extract(subdetector.r, sub.detector.sp)
   values.sub.r[is.na(values.sub.r)] <- 0
   sub.detector.sp <- sub.detector.sp[values.sub.r>0, ]
   sub.detector.sp$count <- values.sub.r[values.sub.r>0]
   
   ### ==== MAKE A PLOTTING FUNCTION ====
   if(plot){
      if(sum(class(data) %in% c("SpatialPolygons","SpatialPoints","SpatialLines"))>0){
         plot(data, col=grey(0.3))
         ## else: raster layer
      }else{plot(data)}
      plot(st_geometry(maindetector.poly), add=TRUE, lwd=3)
      plot(st_geometry(subdetector.poly), add=TRUE)
      plot(st_geometry(sub.detector.sp), pch=19, cex=0.6, col=as.numeric(sub.detector.sp$main.cell.id), add=TRUE)
      points(main.cell.y ~ main.cell.x, data=sub.detector.sp, pch=19, cex=1, col=as.numeric(sub.detector.sp$main.cell.id))
   }
   
   ### ==== OUTPUT ====
   out <- list( detector.sp = sub.detector.sp,
                main.detector.sp = main.detector.sp,
                sub.grid.poly = subdetector.poly,
                grid.poly = maindetector.poly)
}


