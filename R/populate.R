#' @title Populate coda output 
#'
#' @description
#' \code{populate} and \code{get.dim} are internal functions used in \code{ProcessCodaOutput}.
#' 
#' @name populate
#'
#' @param params A \code{character} object.
#' @param input A \code{logical} if DIC statistics should be computed or not. 
#' @param dim A \code{vector} of parameters for which statistics should not be computed.
#' @param simslist A \code{logical}. If TRUE, messages are displayed.
#' @param samples A \code{logical}. If TRUE, messages are displayed.
#' 
#' @return A \code{list} with sims lists, parameter summary values and rhat values.
#' 
#' 
NULL
#' @rdname populate
#' @export
get.dim <- function(params)
{
  
  ##-- Get all unique parameters (i.e., collapse indexed non-scalars)
  ps <- unique(sapply(strsplit(params, "\\["), "[", 1)) 
  ##-- Slice indexes from non-scalar parameter entries
  test <- sapply(strsplit(params, "\\["), "[", 1)
  
  ##-- Calculate dimension for each parameter i
  dim <- lapply(ps, function(i){
    
    ##-- Extract indices from each element j of parameter i
    w <- params[test==i]
    getinds <- lapply(w,FUN=function(j){
      
      w2 <- strsplit(j,'\\[')[[1]][2]
      w3 <- strsplit(w2,"\\]")[[1]] 
      w4 <- as.numeric(unlist(strsplit(w3,",")))
      return(w4)
      
    })
    
    ##-- Get max value from each dimension of i
    collapsedinds <- do.call(rbind,getinds)
    apply(collapsedinds,2,max)  
  })
  
  names(dim) = ps
  dim
}

NULL
#' @rdname populate
#' @export
populate <- function( input,
                      dim,
                      simslist = FALSE,
                      samples = NULL)
  {
  
  if(!simslist){
    charinds <- sub(".*\\[(.*)\\].*", "\\1", names(input), perl=TRUE) 
    fill <- array(NA, dim = dim)
    
    for(i in 1:length(input)){
      ind <- lapply(strsplit(charinds[i], ','), as.integer)[[1]]
      fill[matrix(ind,1)] <- input[i]
    }
  } else {
    charinds <- sub(".*\\[(.*)\\].*", "\\1", colnames(input), perl=TRUE) 
    fill <- array(NA,dim=c(samples,dim))
    for (i in 1:length(charinds)){
      eval(parse(text=paste('fill[','1:',samples,',',charinds[i],']','<- input[,i]',sep="")))
    }
  }
  
  return(fill)
}

