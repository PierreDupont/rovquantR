#' Add a scale to a map plot
#'
#' The \code{addScale} function adds a scale with text to a map.
#' The position, size and orientation of the scale are given.
#'
#' @param x,y \code{Numeric} vectors of values between 0 and 1, one value for each item in the \code{name} list, denoting the relative position of the center of the scale along the x- and y-axes.
#' @param size A \code{numeric} of values between 0 and 1, denoting the size of the image relative to the width of the plot.
#' @param dir A \code{character} denoting the orientation of the scale; \code{"v"} (the default value).
#' 
#' @return This function adds a scale (in km) to an existing map.
#' 
#' @author Pierre Dupont
#'
#' @importFrom graphics par segments text
#' 
#' @examples
#' graphics::par(mar = c(0,0,0,0))
#' plot(sf::st_geometry(COUNTRIES), border = NA, col = "gray80")
#' addScale(x = 0.8, y = 0.2, size = 350000, dir = "h")
#'
#' @rdname addScale
#' @export
addScale <- function(x, y, size, dir = "v"){
  
  ##-- Get dimensions of existing plot
  xLims <- graphics::par("usr")[1:2]
  xRange <- diff(xLims)
  yLims <- graphics::par("usr")[3:4]
  yRange <- diff(yLims)
  base_cex <- graphics::par("din")[1]/6
  
  ##-- Coordinates of the segments
  legend.x <- xLims[1] + x * xRange
  legend.y <- yLims[1] + y * yRange
  
  if(dir == "v"){
    graphics::segments(
      x0 = legend.x, x1 = legend.x,
      y0 = legend.y-size/2, y1 = legend.y + size/2,
      col = "gray30", lwd = base_cex * 3, lend = 2)
    graphics::text(
      x = legend.x - 0.03 * xRange,
      y = legend.y,
      labels = paste0(round(size/1000)," km"),
      srt = 90,
      cex = base_cex)
  } else {
    graphics::segments(
      x0 = legend.x-size/2, x1 = legend.x + size/2,
      y0 = legend.y, y1 = legend.y,
      col = "gray30", lwd = base_cex * 3, lend = 2)
    graphics::text(
      x = legend.x,
      y = legend.y + 0.03 * yRange,
      labels = paste0(round(size/1000)," km"),
      srt = 0,
      cex = base_cex)
  }
}


