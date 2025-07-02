#' Print Bayesian model estimates in a clean way.
#'
#' R utility function to print summary metrics from a vector of posterior samples.
#'
#' The \code{getCleanEstimates} function is used to print the posterior mean (or median) of a vector of posterior MCMC samples, together with the credible interval of the specified width (in parentheses).
#'
#' @param x a vector of posterior MCMC samples. 
#' @param moment an \code{character string} denoting the type of metric to print; can be one of "mean" or "median".
#' @param quantiles a numeric vector of length 2 denoting the boundaries of the credible interval (e.g. c(0.025,0.975), the default value for the 95 percent credible interval).
#' 
#' @return A formated character string with the mean (or median) and the associated credible interval in parentheses.
#' 
#' @author Pierre Dupont
#'
#' @importFrom stats quantile median
#' 
#' @rdname getCleanEstimates
#' @export
getCleanEstimates <- function(
    x,
    moment = "mean",
    quantiles = c(0.025,0.975))
{
  if(moment == "mean"){
    paste0(format(round(mean(x),digits = 2),nsmall = 2)," (",
           format(round(quantile(x, probs = quantiles[1]), digits = 2),nsmall = 2), "-",
           format(round(quantile(x, probs = quantiles[2]), digits = 2),nsmall = 2), ")")
  } else {
    paste0(format(round(median(x), digits = 2),
                  nsmall = 2)," (",
           format(round(quantile(x, probs = quantiles[1]), digits = 2),nsmall = 2), "-",
           format(round(quantile(x, probs = quantiles[2]), digits = 2),nsmall = 2), ")")
  }# else
}