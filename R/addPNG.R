#' Add a .png image to a plot
#'
#' The \code{addPNG} function adds one or several .png images to an already existing plot.
#' The position and size of the image is given by its coordinates relative to the size of the plot.
#'
#' @param x,y \code{Numeric} vectors of values between 0 and 1, one value for each item in the \code{name} list, denoting the relative position(s) along the x- and y-axes.
#' @param name A \code{character} vector to be matched with the name of the image. 
#' @param size A \code{numeric} of values between 0 and 1, denoting the size of the image relative to the width of the plot.
#' @param path (optional) A \code{path} to the directory containing the image file(s).
#' 
#' @return This function adds .png images to an existing plot.
#' 
#' @author Pierre Dupont
#'
#' @importFrom graphics par rasterImage
#' @importFrom png readPNG
#' 
#' @examples
#' graphics::par(mar = c(0,0,0,0))
#' plot(sf::st_geometry(COUNTRIES), border = NA, col = "gray80")
#' addPNG(x = c(0.5,0.8), y = c(0.5,0.15), name = c("prrt","bear"), size = c(0.8,0.2))
#'
#' @rdname addPNG
#' @export
addPNG <- function(x, y, name, size, path = NULL){
  
  ##-- Identify directory containing the images
  if(is.null(path)){
    path <- system.file("images", package = "rovquantR")
  }
  ##-- Initial checks
  if(length(x) != length(name) | length(y) != length(name) | length(size) != length(name)){
    stop("The length of the 'x', 'y' and 'size' arguments must match the length of 'name'.")
  }
  
  ##-- Get dimensions of existing plot
  xLims <- par("usr")[1:2]
  xRange <- diff(xLims)
  yLims <- par("usr")[3:4]
  yRange <- diff(yLims)
  
  ##-- Loop over image names
  for(l in 1:length(name)){
    
    ##-- Identify the image 
    imageIndex <- grep(name[l],list.files(path), ignore.case = TRUE)
    
    ##-- Print warning if image not found
    if(sum(imageIndex)==0){
      warning(paste0("Image '", name[l], "' not found in '", path, "'."))
      next
    }
    
    ##-- Read the image
    pic <- png::readPNG(list.files(path,full.names = TRUE)[imageIndex])
    
    ##-- Get the image size
    picSize <- dim(pic)
    xSize <- size[l] * xRange
    ySize <- xSize*picSize[1]/picSize[2]
    
    ##-- Get the image position
    xPos <- xLims[1] + x[l] * xRange - (xSize/2)
    yPos <- yLims[1] + y[l] * yRange - (ySize/2)
    
    ##-- Print the image
    graphics::rasterImage( pic,
                           xleft = xPos,
                           xright = xPos + xSize,
                           ybottom = yPos,
                           ytop = yPos + ySize)
  }#l
}


