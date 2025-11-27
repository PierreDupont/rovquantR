#' Add population size estimates to plot
#'
#' Function to add population size estimates and country flags to an existing plot.
#'
#' The \code{addPopSize} function uses the outputs from the \code{getDensity} function as input.
#'
#' @param text A named \code{list} containing the population size estimates. 
#' If names are provided, they will be used to match country flags.
#' @param x,y \code{Numeric} vectors of values between 0 and 1, one value for each item in the \code{text} list, denoting the relative position(s) along the x- and y-axes.
#' 
#' @return This function adds population size estimates as text and country flags.
#' 
#' @author Pierre Dupont
#'
#' @importFrom graphics par text rasterImage
#' @importFrom png readPNG
#' 
#' @examples
#' plot(sf::st_geometry(COUNTRIES), border = NA, col = "gray80")
#' 
#' labels <- list(NOR = c(90,123), SWEden = c(20,40))
#' x <- c(0.1,0.6)
#' y <- c(0.7,0.3)
#' addPopSize(x,y,labels)
#'
#' @rdname addPopSize
#' @export
addPopSize <- function(x,y,labels){
  
  ##-- Initial checks
  if(!inherits(labels,"list")){ labels <- list(labels) }
  if(length(x) != length(labels) | length(y) != length(labels)){
    stop("The length of the 'x' and 'y' coordinates must match the length of the list of labels to be printed.")
  }
  
  ##-- Get plot dimensions
  xLims <- par("usr")[1:2]
  xRange <- diff(xLims)
  yLims <- par("usr")[3:4]
  yRange <- diff(yLims)
  base_cex <- par("din")[1] / 6 
  
  
  ##-- Loop over the different labels in the list
  for(l in 1:length(labels)){
    
    ##-- Check if a flag is needed
    flag <- flag2 <- NULL
    if(sum(grep("nor", names(labels)[l], ignore.case = T))>0) {
      flag <- png::readPNG( system.file("images", "nor.png", package = "rovquantR"))
    }
    if(sum(grep("swe", names(labels)[l], ignore.case = T))>0) {
      flag <- png::readPNG( system.file("images", "swe.png", package = "rovquantR"))
    } 
    if(sum(grep("both", names(labels)[l], ignore.case = T))>0) {
      flag <- png::readPNG( system.file("images", "swe.png", package = "rovquantR"))
      flag2 <- png::readPNG( system.file("images", "nor.png", package = "rovquantR"))
    } 
    
    ##-- If not flag
    if(is.null(flag)){
      ##-- Get position of the label
      xPos <- xLims[1] + x[l] * xRange
      yPos <- yLims[1] + y[l] * yRange
      ##-- Add label
      graphics::text( x = xPos,
                      y = yPos, 
                      labels = paste0(labels[[l]], collapse = "-"),
                      cex = base_cex,
                      font = 2)
    } else {
      
      ##-- Get position of the flag
      xPos <- xLims[1] + x[l] * xRange
      yPos <- yLims[1] + y[l] * yRange
      
      ##-- Get size of the flag
      flagSize <- dim(flag)
      xSize <- 0.08 * xRange
      ySize <- xSize*flagSize[1]/flagSize[2]
      
      ##-- add flag(s) to plot
      graphics::rasterImage( flag,
                             xleft = xPos - 0.6*xSize,
                             xright = xPos + 0.4*xSize,
                             ybottom = yPos - 0.5*ySize,
                             ytop = yPos + 0.5*ySize)
      if(!is.null(flag2)){
        graphics::rasterImage( flag2,
                               xleft = xPos - 1.6*xSize,
                               xright = xPos - 0.6*xSize,
                               ybottom = yPos - 0.5*ySize,
                               ytop = yPos + 0.5*ySize)
      }
      
      ##-- Add label
      graphics::text( x = xPos + 0.5*xSize + 0.06*xRange,
                      y = yPos, 
                      labels = paste0(labels[[l]], collapse = "-"),
                      cex = base_cex,
                      font = 2)
    }
  }#l
}

