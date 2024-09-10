#' Scandinavian character translation
#' 
#' Replaces Scandinavian characters for easier processing. 
#'  
#' @param data a character vector where matches are sought, or an object which can be coerced by as.character to a character vector.
#' @param dir.translation (optional) the path to the file containing special characters and their "non-special" translation.
#' 
#' @return the same object than \code{data} with all Scandinavian characters translated.
#' 
#' @author Richard Bischof and Pierre Dupont
#' 
#' @examples
#' translateForeignCharacters("bjørn")
#'
#' @rdname translateForeignCharacters
#' @export
translateForeignCharacters <- function( data,
                                        dir.translation = NULL) 
{
  ##-- Load .RData file with special characters
  if(!is.null(dir.translation)){
   # data(fromto, envir = environment()) 
    #load(system.file("extdata", "CharacterTranslation.RData", package = "rovquantR"))
  #} else {
    load(dir.translation)
  }
  
  ##-- Replace special characters in data
  for(i in 1:nrow(fromto)) {
    data <- gsub( pattern = fromto$from[i],
                  replacement = fromto$to[i], 
                  x = data)
  }#i
  return(data)
}

