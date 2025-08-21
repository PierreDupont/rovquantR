#' @title Create a list of objects to define the habitat from a predefined habitat raster 
#'
#' @description
#' \code{makeHabitatFromRaster} returns a \code{List} of objects necessary to define the habitat in SCR (can be slow for large habitat...)
#'
#' @param poly A \code{sf} object with the study area
#' @param buffer A \code{Numeric} with the size of the buffer area around the study area (in meters in poly is UTM).
#' @param habitat.r A \code{raster} defining suitable habitat.
#' @param plot.check A \code{Logical} to display checking plots (if TRUE).
#' @param tolerance A \code{Numeric}. Default is (res(habitat.r)[1]/5). 
#' This tolerance threshold is used to identify cells that should be considered in the buffer and not searched. 
#' Sometimes raster cells will be identified as not-searched but should actually be searched.
#' This is because the center of the cell does not fall in the polygon used to create the searched area.
#' Instead we identify cells that are in the habitat but "far" from the poly.
#' Those cells are then considered as buffer.
#' 
#' @return A \code{List} object with the coordinates and attributes of the habitat
#'  @return habitat.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat
#'  @return habitat.clip.sp: A \code{SpatialPointsDataFrame} object with the coordinates of the center of cells of the habitat clipped to the suitable habitat.
#'  @return habitat.xy: A \code{DataFrame} object with the x and y coordinates of the habitat cells.
#'  @return IDCells.mx: A \code{matrix} with the ID of each cell of the raster
#'  @return habitat.r: A \code{raster} with the suitable habitat (1) and non suitable (0)
#'  @return habitat.mx: A \code{matrix} with the suitable habitat (1) and non suitable (0)
#'  @return resolution: A \code{numeric} with the resolution used for the raster 
#'  @return buffered.habitat.poly: A \code{SpatialPolygons} that includes the buffer. 
#'  @return habitat.poly: A \code{SpatialPolygons} with the original polygon used. 
#'  @return habitat.index: A \code{vector} with ID cell of the raster that fall within the suitable habitat. 
#'  
#' @author Cyril Milleret
#' 
#' @import sf 
#' @importFrom raster mask crop xyFromCell res as.matrix ncell plot
#' @importFrom dplyr group_by summarise
#' @importFrom graphics plot
#'    
#' @rdname makeHabitatFromRaster
#' @export
makeHabitatFromRaster <- function( 
    poly,
    habitat.r, 			  
    buffer = NULL, 		 
    plot.check = TRUE, 
    tolerance = NULL){ 
  
  if(is.null(tolerance)){
    tolerance <- (raster::res(habitat.r)[1]/5)
    }
  
  ##-- Create buffer around study area polygon
  polyBuffered <- sf::st_buffer(poly, dist = buffer)

  ##-- Mask and crop habitat with the buffer 
  habitat.r <- raster::mask(habitat.r, polyBuffered)
  buffered.habitat.poly <- sf::st_as_sf(stars::st_as_stars(habitat.r), 
                                    as_points = FALSE, merge = TRUE)
  buffered.habitat.poly <- buffered.habitat.poly[buffered.habitat.poly$Habitat>0,]
  habitat.r <- raster::crop(habitat.r, buffered.habitat.poly)
  
  ##-- Create Habitat matrix (habitat : 1 and non-habitat: 0) 
  habitat.r[is.na(habitat.r)] <- 0                                     ## Give 0 values to raster cells outside the study area
  habitat.mx <- as.matrix(habitat.r)                                   ## Convert to matrix
  
  ##-- Give unique IDs to cells 
  IDCells.r <- habitat.r
  IDCells.r[] <- 1:length(IDCells.r)				  ## Cell ID starts from the top left corner 
  IDCells.mx <- raster::as.matrix(IDCells.r)	## Convert to matrix
  
  ##-- Obtain xy coordinates of cells 
  habitat.xy <- raster::xyFromCell(habitat.r, 1:raster::ncell(habitat.r))
  dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
  habitat.xy <- as.data.frame(habitat.xy)
  habitat.sp <- sf::st_as_sf(habitat.xy, coords = c("x", "y"))
  sf::st_crs(habitat.sp) <- sf::st_crs(poly)

  poly$id <- 1
  polyAggregated <- poly %>%
    dplyr::group_by(id) %>%
    dplyr::summarise()
  
  habitat.index <- which(!as.numeric(unlist(sf::st_intersects(habitat.sp, polyAggregated))))
  habitat.clip.sp <- habitat.sp[habitat.index, ]
  
  ##-- Obtain lower and upper cell coordinates
  resolution <- raster::res(habitat.r)[1]
  lower.hab.sp <- data.frame(sf::st_coordinates(habitat.sp) - resolution/2)
  upper.hab.sp <- data.frame(sf::st_coordinates(habitat.sp) + resolution/2)
  colnames(lower.hab.sp) <- colnames(upper.hab.sp) <- c("x", "y")
  upper.hab.sp <- sf::st_as_sf(upper.hab.sp, coords = c("x", "y"))
  lower.hab.sp <- sf::st_as_sf(lower.hab.sp, coords = c("x", "y"))
  sf::st_crs(lower.hab.sp) <- sf::st_crs(poly)
  sf::st_crs(upper.hab.sp) <- sf::st_crs(poly)

  ##-- Create an habitat raster without buffer
  polyBuffered$id <- 1
  polyBuffAggregated <- polyBuffered %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise()
  CellsInBuffer  <- as.numeric(sf::st_intersects(habitat.sp, polyBuffAggregated, sparse = F))
  CellsInBuffer[CellsInBuffer > 0] <- 1
  whichCellsInBuffer <- which(CellsInBuffer %in% 1)
  
  ##-- Get cells in buffer that are "close" from the searched area
  dist <- as.numeric(sf::st_distance(habitat.sp[whichCellsInBuffer, ], polyAggregated))
  whichFarFromBuffer <- whichCellsInBuffer[dist > tolerance]
  habitat.rSearchedAndBuffer <- habitat.rWthBuffer <- habitat.r
  whichFarFromBufferInHabitat <- whichFarFromBuffer[habitat.rWthBuffer[whichFarFromBuffer]%in%1]
  
  ##-- Give cells in the buffer a NA
  habitat.rWthBuffer[whichFarFromBufferInHabitat] <- NA
  habitat.rSearchedAndBuffer[whichFarFromBufferInHabitat] <- 2

  ##-- Visual plotting to check if everything is right 
  if(plot.check){
    raster::plot(habitat.r, col = c("white","green"), legend = F)							    
    raster::plot(habitat.rWthBuffer, add = T, col = c("white","brown"), legend = F)
    graphics::plot(sf::st_geometry(poly), add = TRUE)						
  }
  
  ##-- List of output objects 
  output <- list( habitat.sp = habitat.sp,
                  habitat.clip.sp = habitat.clip.sp,
                  habitat.xy = habitat.xy,
                  IDCells.mx = IDCells.mx,
                  habitat.r = habitat.r,
                  habitat.mx = habitat.mx,
                  resolution = resolution,
                  buffered.habitat.poly = buffered.habitat.poly,
                  habitat.poly = poly,
                  habitat.index = habitat.index,
                  upper.hab.sp = upper.hab.sp,
                  lower.hab.sp = lower.hab.sp,
                  habitat.rWthBuffer = habitat.rWthBuffer,
                  habitat.rSearchedAndBuffer = habitat.rSearchedAndBuffer)
  return(output)
}

