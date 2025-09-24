#' Add population size estimates to plot
#'
#' Function to add population size estimates and country flags to an existing plot.
#'
#' The \code{addPopSize} function uses the outputs from the \code{getDensity} function as input.
#'
#' @param x A named \code{list} containing the population size estimates. 
#' If names are provided, they will be used to match country flags.
#'
#' @return This function adds population size estimates as text and country flags.
#' 
#' @author Pierre Dupont
#'
#' @importFrom raster res image extent writeRaster
#' @importFrom graphics plot layout par segments mtext text rasterImage
#' @importFrom sf st_geometry 
#' @importFrom grDevices png colorRampPalette
#' @importFrom png readPNG
#'
#' @rdname addPopSize
#' @export
addPopSize <- function(x, y, size){ 
  
  ##-- Initial checks
  if(!inherits(x,"list")){ x <- list(x) }
  
  ##-- Check if flag names are given
  if(sum(grep("nor", species, ignore.case = T))>0) { 
    norFlag <- png::readPNG( system.file("images", "nor.png", package = "rovquantR"))
  }
  if(sum(grep("swe", species, ignore.case = T))>0) { 
    sweFlag <- png::readPNG( system.file("images", "swe.png", package = "rovquantR"))
  }
  
  ##-- Set x- and y-positions
  plotSize <- dev.size()
  par("usr")
  
  xLims <- raster::extent(background)[1:2]
  xRange <- diff(xLims)
  yLims <- raster::extent(background)[3:4]
  yRange <- diff(yLims)
  legend.x <- xLims[1] + 0.65 * xRange
  legend.y <- yLims[1] + 0.2 * yRange
  
  ##-- Add Norwegian flag 
  norFlag <- png::readPNG( system.file("images", "nor.png", package = "rovquantR"))
  norSize <- dim(norFlag)
  xPos <- xLims[1] + 0.05 * xRange
  xSize <- 0.1 * xRange
  yPos <- yLims[1] + 0.8 * yRange
  ySize <- xSize*picSize[1]/picSize[2]
  rasterImage( norFlag,
               xleft = xPos,
               xright = xPos + xSize,
               ybottom = yPos,
               ytop = yPos + ySize)
  
  ##-- Add abundance estimate
  text(x = xPos + xSize + 0.1 * xRange,
       y = yPos + ySize/2, 
       labels = paste0(round(q95[1]), "-", round(q95[2])),
       cex = 1.2, font = 2)
  
  ##-- Add caption
  mtext(text = paste0("Density map and ranges of abundance \nestimated for wolves in ",
                      names(estimates)[length(density)]),
        side = 1,line = 2, adj = 0.5, cex = 1.2, font = 2)
} 

