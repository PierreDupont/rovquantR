#' Read and combine multiple files
#'
#' The \code{readMultiples} function identifies and loads multiple files matching a set of characteristics, before combining them using rbind.
#' 
#' @name readMultiples
#'
#' @param path A character string denoting the data directory path.
#' @param extension A character string denoting the type of file to be loaded (can be one of ".csv", ".RData", "xls", "xlsx")
#' @param pattern (Optional) An additional character string to be matched with the file name.
#' @param sep The field separator character used when loading \code{".csv"} files. Values on each line of the file are separated by this character (default is ","),.
#' @param dec The character used in the file for decimal points.
#' @param ... additional optional parameters.
#' 
#' @return The data loaded
#'
#' @author Pierre Dupont
#'
#' @importFrom readxl read_excel
#' @importFrom utils read.csv 
#' 
NULL
#' @rdname readMultiples
#' @export
##-- Generic function to read the most recent file with a given extension (and optionally pattern)
readMultiples <- function( 
    path,
    extension = ".csv",
    pattern = NULL,
    sep = ",",
    dec = ".",
    ...)
{
  
  ##-- List all files with the requested extension (including in sub-directories)
  infiles <- list.files( path = path,
                         pattern = extension,
                         recursive = TRUE,
                         ignore.case = TRUE)
  if(length(infiles) < 1){stop(paste0("No file was found with a '", extension,"' extension in the requested directory."))}
  
  ##-- Further subset to files that match the requested pattern 
  if(!is.null(pattern)){
    whichFiles <- grep(pattern, infiles)
    if(length(infiles) < 1){stop(paste0("No file was found that matched the pattern: '", pattern,"'."))}
    infiles <- infiles[whichFiles] 
  }

  ##-- Identify files that match the requested extension and pattern
  message(paste0('Loading files ', infiles, '...\n'))

  ##-- read the different .csv files
  if(length(grep("csv", extension, ignore.case = T)) > 0){
    data_list <- lapply( infiles,
                         function(x){
                           utils::read.csv( file = file.path(path,x),
                                            header = TRUE,
                                            sep = sep,
                                            dec = dec,
                                            ...)})
  }
  
  ##-- function to read the most recent .xls or .xlsx file
  if(length(grep("xls", extension, ignore.case = T)) > 0){
    data_list <- lapply( infiles,
                         function(x){
                           readxl::read_excel(path = file.path(path, x), ...)
                         })
  }
  
  ##-- function to load and return the most recent .RData file
  if(length(grep("RData", extension, ignore.case = T)) > 0){
    fileName <- file.path(path,infiles)
    readRData <- function(fileName, ...){
      load(fileName, ...)
      mget(ls()[ls() != "fileName"])
    }
    data_list <- lapply( fileName,
                         function(x){
                           readRData(x)
                           })
  }

  ##-- combine data files
  data <- do.call(rbind, data_list)
  
  ##-- Output
  return(data)
}
