#' @title Sanitize nimbleSCR input data
#'
#' @description
#' R utility function to sanitize one or multiple individual detection datasets together with the corresponding detectors information.
#' This function ensures data format and names are compatible with rovquantR and nimbleSCR functions.
#' 
#' The \code{sanitizeData} function is used among other in the data2nimbleSCR function to help formatting SCR and OPSCR data for use with nimbleSCR.
#'
#' @name sanitize
#'
#' @param data a data.frame (or list of data.frames if the model integrates multiple datasets) containing individual detections (in rows).
#' @param name (optional) a character string with the name of the data frame tested (used to report potential errors).
#' 
#' @author Pierre Dupont
#' 
#' @references
#' nimbleSCR manual: ... 
#'
#' @examples 
#' n.samples <- 150
#' n.ids <- 40
#' n.detectors <- 10
#' n.sessions <- 3
#' n.occasions <- 5
#' 
#' ## Example 1: single dataset - OPSCR
#' dataset1 <- cbind.data.frame( "iD" = sample(n.ids,size = n.samples,replace = TRUE),
#'                              "detECtor" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                              "sub.detector" = sample(n.occasions, size = n.samples,replace = TRUE),
#'                              "ids.covs" = sample(n.occasions, size = n.samples,replace = TRUE),
#'                              "Year" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                              "occasion" = sample(n.occasions,size = n.samples,replace = TRUE))
#' head(dataset1)
#' 
#' detectors <- expand.grid(list( "id" = 1:n.detectors,
#'                                #"occasion" = 1:n.occasions,
#'                                "session" = 1:n.sessions))
#' DETS <- cbind.data.frame( "detectors" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                           "session" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                           "occasion" = sample(n.occasions,size = n.samples,replace = TRUE))
#' head(DETS)
#' 
#' sanitizeInput(dataset1, DETS)
#' lapply(dataset1, head)
#' lapply(DETS, head)
#'
#' 
#' ## Example 2: multiple datasets - OPSCR
#' myData <- list( data_1 =  cbind.data.frame( "IDs" = sample(n.ids,size = n.samples,replace = TRUE),
#'                                             "detECtor" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                                             "Years" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                                             "Occasion" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                                             "sub.detector" = sample(n.occasions,size = n.samples,replace = TRUE)),
#'                 SCR = cbind.data.frame( "id" = sample(n.ids,size = n.samples,replace = TRUE),
#'                                         "DeteCtors" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                                         "Occasion" = sample(n.sessions,size = n.samples,replace = TRUE)))
#' lapply(myData, head)
#' 
#' DETS <- list( cbind.data.frame( "detectors" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                                 "session" = sample(n.sessions,size = n.samples,replace = TRUE),
#'                                 "occasion" = sample(n.occasions,size = n.samples,replace = TRUE),
#'                                 "cov1" = rnorm(n.samples)),
#'               cbind.data.frame( "detectors" = sample(n.detectors,size = n.samples,replace = TRUE),
#'                                 "occasion" = sample(n.occasions,size = n.samples,replace = TRUE)))
#' lapply(DETS,head)
#' 
#' sanitizeInput(myData,DETS)
#' lapply(myData, head)
#' lapply(DETS, head)

NULL
#' @rdname sanitize
#' @export
sanitizeInput <- function(
    data,
    detectors)
{
  ##-- 1. Get input names (needed as this function modifies the `data` and `detectors` objects "in place)
  data_name <- deparse(substitute(data))
  detectors_name <- deparse(substitute(detectors))
  
  
  ##-- 2.Turn inputs into lists (even if only one dataset)
  if(!inherits(data, "list")){ 
    data <- list(data)
    names(data) <- "data"
  }
  if(!inherits(detectors, "list")){ 
    detectors <- list(detectors)
    names(detectors) <- "detectors"
  }
  
  
  ##-- 3.Check if the number of 'detectors' matches the number of 'data' 
  if(length(detectors) != length(data)){
    if(length(detectors) == 1){
      warning("There is only one `detectors` dataframe provided for multiple `data`.\n Recycling the same `detectors` for all `data`.")
      detectors <- lapply(1:length(data), function(x)detectors[[1]])
    } else {
      stop("The number of `data` datasets is different from the number of `detectors` provided.")
    }
  }
  
  
  ##-- 4. Check list names
  if(is.null(names(data))){ dataNames <- paste0("data", 1:length(data)) } else { dataNames <- names(data) }
  if(is.null(names(detectors))){ detNames <- paste0("detectors", 1:length(detectors)) } else { detNames <- names(detectors)} 
  
  
  ##-- 5. Check column names 
  data <- lapply( 1:length(data), function(x)fixNames(data = data[[x]], name = paste0(data_name,"[[",x,"]]")))
  names(data) <- dataNames
  detectors <- lapply( 1:length(detectors), function(x)fixNames(data = detectors[[x]], name = paste0(detectors_name,"[[",x,"]]"), is.data = FALSE))
  names(detectors) <- detNames
  
  
  ##-- 6. Check compatibility of `data` and `detectors`
  for(l in 1:length(data)){
    
    ##-- Make sure `data` is a dataframe
    if(!inherits(data[[l]], "dataframe")){data[[l]] <- as.data.frame(data[[l]])}
    ##-- Make sure `detector` is a dataframe
    if(!inherits(detectors[[l]], "dataframe")){detectors[[l]] <- as.data.frame(detectors[[l]])}
    
    ##-- Make sure the numbers of sessions & occasions are consistent between `data` & `detectors`
    ##-- Check if multiple primary sessions
    if("session" %in% names(data[[l]])) {
      dat.sess <- unique(data[[l]]$session)
      
      ##-- Check that 'detectors' also contains multiple sessions
      if("session" %in% names(detectors[[l]])){
        det.sess <- unique(detectors[[l]]$session)
        if(!all(dat.sess %in% det.sess)){
          stop(paste0("Some primary sessions in `data` are not present in `detectors`."))
        }
      } else {
        warning("`data` contains multiple primary sessions but `detectors` does not.\nUsing the same `detectors` for all primary sessions in `data`.")
        detectors[[l]] <- do.call(rbind, lapply(dat.sess, function(x)cbind.data.frame(detectors[[l]], "session" = x)))
      }
    } else {
      if("session" %in% names(detectors[[l]])){
        stop(paste0("`detectors` contains multiple primary sessions but `data` does not.")) 
      } else { 
        detectors[[l]]$session <- 1
        data[[l]]$session <- 1
        dat.sess <- unique(data[[l]]$session)
      }
    }
    
    ##-- Check if multiple secondary occasions
    if("occasion" %in% names(data[[l]])) {
      dat.occ <- unique(data[[l]]$occasion)
      
      ##-- Check that 'detectors' also contains multiple occasions
      if("occasion" %in% names(detectors[[l]])){
        det.occ <- unique(detectors[[l]]$occasion)
        if(!all(dat.occ %in% det.occ)){
          stop(paste0("Some occasions in `data` are not present in `detectors`."))
        }
      } else {
        stop(paste0("`data` contains multiple occasions but `detectors` does not."))
      }
    } else {
      if("occasion" %in% names(detectors[[l]])){
        stop(paste0("`detectors` contains multiple secondary occasions but `data` does not.")) 
      } else { 
        detectors[[l]]$occasion <- 1
        data[[l]]$occasion <- 1
      }
    }  
    
  }#l
  
  
  ##-- 7.Output
  ##-- Modify `data` and `detectors` directly in the global environment
  assign(data_name, data, envir = .GlobalEnv)
  assign(detectors_name, detectors, envir = .GlobalEnv)
  
  invisible(TRUE)
}


NULL
#' @rdname sanitize
#' @export
##-- Function to fix names 
fixNames <- function( 
    data,
    name = NULL,
    is.data = TRUE)
{
  
  ##-- 1 - "id"
  if(is.data){
    idTest <- grep("^ids?$", names(data), ignore.case = T)
    if(length(idTest) == 0)stop(paste0("Missing column 'id' in '",name,"'."))
    if(length(idTest) > 1)stop(paste0("More than one column matching 'id' in :'",name,"'."))
    names(data)[idTest] <- "id"
  }
  
  ##-- 2 - "detector"
  detTest <- grep("^detectors?$", names(data), ignore.case = T)
  if(length(detTest) == 0)stop(paste0("Missing column 'detector' in '",name,"'."))
  if(length(detTest) > 1)stop(paste0("More than one column matching 'detector' in '",name,"'."))
  names(data)[detTest] <- "detector"
  
  ##-- 3 - "session"
  sessTest <- grep("^sessions?$", names(data), ignore.case = T)
  if(length(sessTest) > 1)stop(paste0("More than one column matching 'session' in '",name,"'."))
  if(length(sessTest) == 0){
    sessTest <- grep("^years?$", names(data), ignore.case = T)
    if(length(sessTest) > 1)stop(paste0("More than one column matching 'year' in '",name,"'."))
  }
  names(data)[sessTest] <- "session"
  
  ##-- 4 - "occasion"
  occTest <- grep("^occasions?$", names(data), ignore.case = T)
  if(length(occTest) > 1)stop(paste0("More than one column matching 'occasion' in '",name,"'."))
  names(data)[occTest] <- "occasion"
  
  ##-- 5 - "sub.detector"
  subTest <- c( grep("^sub.detectors?$", names(data), ignore.case = T),
                grep("^subdetectors?$", names(data), ignore.case = T))
  if(length(subTest) > 1)stop(paste0("More than one column matching 'subdetector' in '",name,"'."))
  names(data)[subTest] <- "subdetector"
  
  ##-- Output
  return(data)
}


# sanitizeData <- function(data)
# {
#   data_name <- deparse(substitute(data))
#   detectors_name <- deparse(substitute(detectors))
#   
#   
#   ##-- Turn into a list (even if only one dataset)
#   if(!inherits(data, "list")){ 
#     data <- list(data)
#     names(data) <- "data"
#   }
#   
#   ##-- Check list names
#   if(is.null(names(data))){
#     dataNames <- paste0("data", 1:length(data))
#   } else {
#     dataNames <- names(data)
#   }
#   
#   ##-- Check names 
#   data <- lapply( 1:length(data), function(x)checkNames(data = data[[x]], name = dataNames[x]))
#   names(data) <- dataNames
#   
#   ##-- Output
#   assign(input_name, data, envir = .GlobalEnv)
#   
#   invisible(TRUE)
# }
# 
# sanitizeDetectors <- function(
    #     detectors,
#     data)
# {
#   ##-- Turn into a list (even if only one dataset)
#   if(!inherits(detectors, "list")){ 
#     detectors <- list(detectors)
#     names(detectors) <- "detectors"
#   }
#   
#   ##-- Check if the number of 'detectors' equals the number of 'data' 
#   if(length(detectors) != length(data)){
#     if(length(detectors) > length(data)){
#       stop("The number of 'data' datasets is less than the number of 'detectors' provided.")
#     } else {
#       detectors <- lapply(1:length(data), function(x)detectors)
#     }
#   }
#   
#   ##-- Check list names
#   if(is.null(names(detectors))){
#     detNames <- paste0("detectors", 1:length(detectors))
#   }
#   
#   ##-- Check names 
#   detectors <- lapply( 1:length(detectors), function(x)checkNames(data = detectors[[x]], name = names(detectors)[x]))
#   names(detectors) <- detNames
#   
#   ##-- Output
#   return(detectors)
# }
