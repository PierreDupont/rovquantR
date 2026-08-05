#' Add a color bar to a plot
#'
#' The \code{addColorBar} function adds a color bar to an already existing plot.
#' It is used to e.g. plot a continuous legend on a map.
#' The position and size of the bar are given by its position relative to the size of the plot (starting in the bottom left corner).
#'
#' @param zlim \code{Numeric} values giving the min and max values of the color bar.
#' @param cols A \code{character} vector of colors to display on the bar. 
#' @param bar_width, bar_height \code{Numeric} values between 0 and 1, denoting the size of the bar relative to the size of the plot.
#' @param x_pad, y_pad \code{Numeric} values between 0 and 1, denoting the position of the bar relative to the size of the plot.
#' @param lab_name (optional) A \code{character} string to display on the side of the bar.
#' @param lab_cex (optional) A \code{numeric} value denoting the relative size of the label text.
#' @param lab_namee (optional) A \code{character} string to display on the side of the bar.
#' @param tick_cex, tick_num, tick digits (optional) \code{numeric values denoting the relative size, number and digits to display for the values alongside the bar.

#' @return This function adds a color bar to an existing plot.
#' 
#' @author Pierre Dupont
#'
#' @importFrom graphics par rect text
#' 
#' @examples
#' graphics::par(mar = c(0,0,0,0))
#' plot(sf::st_geometry(COUNTRIES), border = NA, col = "gray80")
#' zlim <-  c(0,1)
#' cols <- hcl.colors(100, "Viridis")
#' addColorBar( zlim = zlim, cols = cols, y_pad = 0.5, x_pad = 0.15, lab_name = "something", lab_cex = 1.5)
#'
#' @rdname addColorBar
#' @export
addColorBar <- function(
    zlim,
    cols,
    bar_width = 0.035,
    bar_height = 0.30,
    x_pad = 0.04,
    y_pad = 0.05,
    lab_name = NULL,
    lab_cex = 0.9,
    tick_cex = 0.9,
    tick_num = 5,
    tick_digits = 1) {
  
  usr <- par("usr")
  xlim <- usr[1:2]
  ylim <- usr[3:4]
  dx <- diff(xlim)
  dy <- diff(ylim)
  
  # Color bar position in lower-left
  x_left <- xlim[1] + x_pad * dx
  x_right <- x_left + bar_width * dx
  y_bot <- ylim[1] + y_pad * dy
  y_top <- y_bot + bar_height * dy
  
  # Draw the color bar
  yy <- seq(y_bot, y_top, length.out = length(cols) + 1)
  for (i in seq_along(cols)) {
    rect( xleft = x_left,
          ybottom = yy[i],
          xright = x_right,
          ytop = yy[i + 1],
          col = cols[i],
          border = cols[i],
          xpd = NA)
  }#i
  
  # Outline around color bar
  rect( xleft = x_left,
        ybottom = y_bot,
        xright = x_right,
        ytop = y_top,
        border = "black",
        lwd = 0.8,
        xpd = NA)
  
  # Optional variable name vertically on the left side of the color bar
  if (!is.null(lab_name)) {
    text( x = x_left - 0.02 * dx,
          y = (y_bot + y_top) / 2,
          labels = lab_name,
          font = 2,
          srt = 90,
          adj = c(0.5, 0.5),
          cex = lab_cex,
          xpd = NA)
  }
  
  # Tick labels:
  vals <- seq(zlim[1], zlim[2], length.out = tick_num)
  ypos <- seq(y_bot, y_top, length.out = tick_num)
  text( x = x_right + 0.01 * dx,
        y = ypos,
        labels = format(round(vals, tick_digits), trim = TRUE),
        adj = c(0, 0.5),
        cex = tick_cex,
        xpd = NA)
}
