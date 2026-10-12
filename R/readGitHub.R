#' Read files from GitHub 
#'
#' The \code{readGitHub} function identifies and loads specific files hosted on GitHub.
#' 
#' @name readGitHub
#'
#' @param repo A character string with the name of the target repository.
#' @param subdirectory (optional) A character string denoting the subfolders where the file(s) are located.
#' @param extension A character string denoting the type of file extensions targeted. Currently, this function works with .tif and .csv files only.
#' @param base_url A character string denoting the raw url of the GitHub repository.
#' 
#' @return The data loaded (in a list if multiple matched files).
#'
#' @author Pierre Dupont
#'
#' @importFrom gh gh
#' @importFrom terra rast 
#' @importFrom utils download.file 
#' 
#' @examples
#' example1 <- readGitHub(
#'   owner = "richbi",
#'   repo = "RovQuantPublic",
#'   target = "fagrapport112_Wolverine2016-2025/output",
#'   extension = "Table1.csv")
#'   
#' example2 <- readGitHub(
#'   owner = "PierreDupont",
#'   repo = "rovquantR",
#'   target = "COUNTIES.rda")
#'   
#' example3 <- readGitHub(
#'   owner = "PierreDupont",
#'   repo = "Microcebe_Public",
#'   target = "Microcebe_Mandena_data.RData")
#'   
#' example4 <- readGitHub(
#'   owner = "richbi",
#'   repo = "RovQuantPublic",
#'   target = "fagrapport112_Wolverine2016-2025/output/rasters/RasterForRovbase",
#'   extension = ".tif")
#' 
NULL
#' @rdname readGitHub
#' @export
readGitHub <- function(
    owner = "richbi",
    repo = "RovQuantPublic",
    target = NULL,
    extension = NULL
){
  ##-- Create raw_url base for that repository
  base_url <- paste0("https://raw.github.com/",owner,"/",repo,"/master")
  
  ##-- Fetch the recursive tree of the repository
  repo_tree <- gh::gh("GET /repos/:owner/:repo/git/trees/:tree_sha",
                      owner = owner,
                      repo = repo,
                      tree_sha = "main", # Can be a branch name or commit SHA
                      recursive = 1)     # 1 tells GitHub to return all nested subfolders
  
  ##-- Extract the file paths into a vector
  gh_paths <- sapply(repo_tree$tree, function(x) x$path) 
  
  ##-- Identify the paths to the correct directory
  if(!is.null(target)){
    gh_paths <- gh_paths[grep(target, gh_paths)]
  }
  
  ##-- Identify the file(s)
  if(!is.null(extension)){
    gh_paths <- gh_paths[grep(extension, gh_paths)]
  }
  
  ##-- List raw GitHub URLs
  gh_url <- file.path(base_url, gh_paths)
  
  
  ##-- Loop over files 
  out <- list()
  for( f in 1:length(gh_url) ){
    
    ##-- Load .csv files from GitHub
    if( grepl(pattern =  ".csv", x = gh_paths[f]) ){
      out[[f]] <- read.csv(
        gh_url[f],
        header = TRUE,
        sep = ",",
        quote = "\"",
        dec = ".",
        fill = TRUE, 
        comment.char = "")
    } 
    
    ##-- Load .tif files from GitHub
    if( grepl(pattern =  ".tif", x = gh_paths[f]) ){
      out[[f]] <- readGitHub_tif(gh_url[f])
    }
    
    ##-- Load .rda files from GitHub
    if( grepl(pattern = ".rda", x = gh_paths[f]) ){
      out[[f]] <- readGitHub_rda(raw_url = gh_url[f])
    }
    
    ##-- Load .RData files from GitHub
    if( grepl(pattern = ".RData", x = gh_paths[f]) ){
      out[[f]] <- readGitHub_RData(raw_url = gh_url[f])
    }
  }#f
  
  ##-- Return output
  if (length(out) == 1){ out <- out[[1]] } 
  return(out)
}

NULL
#' @rdname readGitHub
#' @export
readGitHub_tif <- function(raw_url){
  
  ## Create a temporary file path
  temp_tif <- tempfile(fileext = ".tif")
  
  ## Download the binary file from GitHub
  utils::download.file(raw_url,
                       destfile = temp_tif,
                       mode = "wb",
                       quiet = TRUE)
  
  ## Read and return the GeoTIFF data
  spatial_raster <- terra::rast(temp_tif)
  return(spatial_raster)
}

NULL
#' @rdname readGitHub
#' @export
readGitHub_rda <- function(raw_url){
  
  ## Create a temporary file path
  temp_rda <- tempfile(fileext = ".rda")
  
  ## Download the binary file from GitHub
  utils::download.file(raw_url,
                       destfile = temp_rda,
                       mode = "wb",
                       quiet = TRUE)
  
  ## Load .Rdata file into a R object
  loadRData <- function(fileName,...){
    load(fileName,...)
    mget(ls()[ls() != "fileName"])
  }
  data <- loadRData(temp_rda)

  ## Read and return the .rda data
  return(data)
}

NULL
#' @rdname readGitHub
#' @export
readGitHub_RData <- function(raw_url){
  
  ## Create a temporary file path
  temp_RData <- tempfile(fileext = ".RData")
  
  ## Download the binary file from GitHub
  utils::download.file(raw_url,
                       destfile = temp_RData,,
                       mode = "wb",
                       quiet = TRUE)
  
  ## Load .Rdata file into a R object
  loadRData <- function(fileName,...){
    load(fileName,...)
    mget(ls()[ls() != "fileName"])
  }
  data <- loadRData(temp_RData)

  ## Read and return the .rda data
  return(data)
}