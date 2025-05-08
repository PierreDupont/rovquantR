#' Read most recent version of a file
#'
#' The \code{readMostRecent} functions provide a suite of utility functions to 
#' identify and load the most recent file matching a set of characteristics.
#' 
#' @name readMostRecent
#'
#' @param path Vector of individual detection frequencies.
#' @param extension Integer specifying the number of realizations to generate. 
#' Only n = 1 is supported.
#' @param pattern Vector of indices of traps where the detections in \emph{x}
#'  were recorded; from the \emph{detIndices} object returned by the \code{\link{getSparseY}} function. 
#' @param returnDate Number of traps with at least one detection recorded in 
#' \emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function. 
#' @param sep Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).
#' @param dec Baseline detection probability (scalar) used in the half-normal 
#' detection function. For trap-specific baseline detection probabilities use 
#' argument \emph{p0Traps} (vector) instead.
#' @param fileEncoding Vector of baseline detection probabilities for each trap 
#' used in the half-normal detection function. When \emph{p0Traps} is used, 
#' \emph{p0} should not be provided. 
#' @param x Scale parameter of the half-normal detection function.
#' @param file Individual activity center x- and y-coordinates scaled to the 
#' habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param ... additional optional parameters.
#' 
#' @return The data loaded
#'
#' @author Pierre Dupont
#'
#' @importFrom readxl read_excel
#' @importFrom readr guess_encoding 
#' @importFrom utils read.csv write.csv
#' 
NULL
#' @rdname readMostRecent
#' @export
##-- Generic function to read the most recent file with a given extension (and optionally pattern)
readMostRecent <- function( 
    path,
    extension = ".csv",
    pattern = NULL,
    returnDate = FALSE,
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
  
  ##-- Extract time last modified for all retained files
  modTime <- lapply(file.path(path, infiles), file.mtime)
  
  ##-- Identify most recent file that matches the requested extension and pattern
  lastFile <- which.max(unlist(modTime))
  message(paste0('Loading file ', infiles[lastFile], ' last modified on ', as.character(modTime[[lastFile]]), '...\n'))
  date <- as.Date(modTime[[lastFile]])
  
  ##-- read the most recent .csv file
  if(length(grep("csv", extension, ignore.case = T)) > 0){
    data <- utils::read.csv( file = file.path(path,infiles[lastFile]),
                      header = TRUE,
                      sep = sep,
                      dec = dec,
                      ...)
  }
  
  ##-- function to read the most recent .xls or .xlsx file
  if(length(grep("xls", extension, ignore.case = T)) > 0){
    data <- readxl::read_excel(path = file.path(path, infiles[lastFile]), ...)
  }
  
  ##-- function to load and return the most recent .RData file
  if(length(grep("RData", extension, ignore.case = T)) > 0){
    fileName <- file.path(path,infiles[lastFile])
    readRData <- function(fileName, ...){
      load(fileName, ...)
      mget(ls()[ls() != "fileName"])
    }
    data <- readRData(fileName)
  }
  
  ##-- Output
  if(returnDate) {
    return(list("data" = data,
                "date" = date))
  } else {
    return(data)
  }
  
}

NULL
#' @rdname readMostRecent
#' @export
##-- Function to read the most recent .csv file
##-- (More restrictive than 'readMostRecent()')
readMostRecent.csv <- function( 
    path,
    returnDate = F,
    fileEncoding = NULL,
    ...)
  {
  infiles <- list.files(path = path, pattern = ".csv", recursive = TRUE)
  modTime <- lapply(file.path(path, infiles), file.mtime)
  lastFile <- which.max(unlist(modTime))
  message(paste0('Loading file ', infiles[lastFile], ' last modified on ', as.character(modTime[[lastFile]]),'...\n'))
  date <- as.Date(modTime[[lastFile]])
  
  if(!is.null(fileEncoding)){
    data <- read.csv(file = file.path(path,infiles[lastFile]), fileEncoding = fileEncoding)
  } else{
    bestEncoding <- guess_encoding(file.path(path,infiles[lastFile]), n_max = 1000)$encoding[1]
    data <- read.csv(file = file.path(path,infiles[lastFile]), fileEncoding = bestEncoding) 
  }
  
  if(returnDate) {
    return(list("data" = data,
                "date" = date))
  } else {
    return(data)
  }
}

NULL
#' @rdname readMostRecent
#' @export
##-- function to read the most recent .xls or .xlsx file
##-- (More restrictive than 'readMostRecent()')
readMostRecent.excel <- function( 
    path,
    returnDate = F,
    ...)
  {
  infiles <- list.files(path = path, pattern = ".xls", recursive = TRUE)
  modTime <- lapply(file.path(path, infiles), file.mtime)
  lastFile <- which.max(unlist(modTime))
  message(paste0('Loading file ', infiles[lastFile], ' last modified on ', as.character(modTime[[lastFile]]),'...\n'))
  date <- as.Date(modTime[[lastFile]])
  data <- readxl::read_excel(path = file.path(path,infiles[lastFile]), ...)
  if(returnDate) {
    return(list("data" = data,
                "date" = date))
  } else {
    return(data)
  }
}

NULL
#' @rdname readMostRecent
#' @export
##-- function to load the most recent .RData file
##-- (More restrictive than 'readMostRecent()')
readMostRecent.RData <- function( 
    path,
    pattern = ".RData",
    returnDate = F,
    ...)
  {
  infiles <- list.files(path = path, pattern = pattern, recursive = TRUE)
  modTime <- lapply(file.path(path, infiles), file.mtime)
  lastFile <- which.max(unlist(modTime))
  message(paste0('Loading file ', infiles[lastFile], ' last modified on ', as.character(modTime[[lastFile]]),'...\n'))
  date <- as.Date(modTime[[lastFile]])
  fileName <- file.path(path,infiles[lastFile])
  loadRData <- function(fileName,...){
    load(fileName,...)
    mget(ls()[ls() != "fileName"])
  }
  data <- loadRData(fileName)
  if(returnDate) {
    return(list("data" = data,
                "date" = date))
  } else {
    return(data)
  }
}

NULL
#' @rdname readMostRecent
#' @export
##-- function to get the most recent modification date
getMostRecent <- function(
    path,
    pattern)
  {
  
  dirTest <- dir( path = path,
                  pattern = pattern,
                  recursive = TRUE,
                  ignore.case = TRUE,
                  full.names = T)
  
  if (length(dirTest) < 1) {
    stop(paste0("No file matching the pattern '", 
                pattern, "' was found in the specified directory: \n",  path))
  }
  
  dirTest |>
    file.mtime() |>
    max() |>
    as.Date()
}

NULL
#' @rdname readMostRecent
#' @export
##-- function to write ".csv" (but first checks if the file exists) 
writeMostRecent.csv <- function(
    x,
    file)
  {
  if(file.exists(file)){
    message(paste0("A file named '", file, "' already exists in the specified directory."))
    message("Are you sure you want to proceed and overwrite existing data? (y/n) ")
    question1 <- readLines(n = 1)
    if(regexpr(question1, 'y', ignore.case = TRUE) != 1){
      message("Not overwriting existing file...\n")
    } else {
      message(paste0("Now overwriting '", file, "'.\n"))
      utils::write.csv( x = x, file = file, row.names = F)
    }
  } else {
    utils::write.csv( x = x, file = file, row.names = F)
  }
}