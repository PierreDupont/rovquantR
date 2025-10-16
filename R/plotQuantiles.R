#' Plot quantiles
#'
#' R utility function to plot the quantiles of a distribution of samples as a bar.
#'
#' @param x Samples for which to plot the quantiles.
#' @param at x-axis coordinates of the bar. 
#' @param width width of the bar. 
#' @param quantiles a vector of two numeric values denoting the quantiles to be represented.
#' @param quantiles2 a second vector of numeric values denoting the second set of quantiles to be represented.
#' @param col color to be used for plotting.
#'
#' @return This function plots quantile bars with the specified characteristics.
#'
#' @author Pierre Dupont
#'
#' @importFrom graphics polygon 
#' @importFrom stats quantile
#' @importFrom grDevices adjustColor
#' 
#' @examples   
#' plot(1, xlim = c(0, 2), ylim = c(0,100), type = "n")
#' plotQuantiles(x = rnorm(1000, 60, sd = 12), at = 1, width = 0.18, quantiles2 = c(0.25,0.75))
#' 
#' @rdname plotQuantiles
#' @export
plotQuantiles <- function(
    x, 
    at,
    width = 0.15, 
    quantiles = c(0.0275,0.975),
    quantiles2 = c(0.25,0.75),
    col = adjustcolor("firebrick3",0.5))
{
  q1 <- quantile(x, prob = quantiles)
  polygon(x = c(at - width, at + width, at + width, at - width),
          y = c(q1[1], q1[1],  q1[2], q1[2]), 
          col = col,
          border = NA)
  
  if(!is.null(quantiles2)){
    q2 <- quantile(x,prob = quantiles2)
    polygon(x = c(at - width, at + width, at + width, at - width),
            y = c(q2[1], q2[1],  q2[2], q2[2]), 
            col = col,
            border = NA)
  }
}