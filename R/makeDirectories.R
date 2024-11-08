#' @title Prepare folder structure
#'
#' @description \code{makeDirectories} creates a set of folders to store all the necessary
#' data for hte OPSCR ROvQuant analysis. 
#' 
#' @param path A path indicating where to create the different folders.
#' @param twoSex A \code{logical}. If TRUE, sex-specific folders are created for the female and male inputs and outputs files, respectively.
#'  
#' @return boolean invisible(FALSE) if nothing was created, invisible(TRUE) if folders were created in \emph{path}.
#'
#' @importFrom fs dir_tree
#' 
#' @examples \dontrun{makeDirectories()}
#' 
#' @rdname makeDirectories
#' @export
makeDirectories <- function( path = NULL,
                             two.sex = T,
                             show.dir){
  ##-- Check if the directory already exists 
  if(is.null(path)){path <- getwd()}
  
  if(dir.exists(path)){
    split <- unlist(strsplit(path, split = "/"))
    message(paste0("A folder named '", split[length(split)], "' already exists in the specified directory:"))
    message(paste(split[-length(split)], collapse = "/"))
    } 

  ##-- Set-up directory structure
  if(two.sex){
  dir.create( file.path(path, "nimbleInFiles/Hann"), showWarnings = F, recursive = T)
  dir.create( file.path(path, "nimbleInFiles/Hunn"), showWarnings = F, recursive = T)
  dir.create( file.path(path, "nimbleOutFiles/Hann"), showWarnings = F, recursive = T)
  dir.create( file.path(path, "nimbleOutFiles/Hunn"), showWarnings = F, recursive = T)
  } else {
    dir.create( file.path(path, "nimbleInFiles"), showWarnings = F)
    dir.create( file.path(path, "nimbleOutFiles"), showWarnings = F)  
  }
  dir.create( file.path(path, "data"), showWarnings = F)
  dir.create( file.path(path, "figures"), showWarnings = F)
  dir.create( file.path(path, "tables"), showWarnings = F)
  dir.create( file.path(path, "rasters"), showWarnings = F)
  dir.create( file.path(path, "reports"), showWarnings = F)
  
  ##-- Display contents of directories in a tree-like format
  if(show.dir){
    message("\n The following directory  structure was created for this analysis:\n")
    fs::dir_tree(path)
  }
}
