#' Character subsetting
#' 
#' Identifies and returns the \code{n} last characters from the character string \code{x}. 
#'  
#' @param x a \code{character} string to substract.
#' @param n a \code{integer} indicating the number of characters to return
#' 
#' @return A vector of window sizes.
#' 
#' @author Pierre Dupont
#' 
#' @examples
#' x <- "this.character.string"
#' n <- 6
#' substrRight(x, n)
#'
#' @rdname substrRight
#' @export
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
