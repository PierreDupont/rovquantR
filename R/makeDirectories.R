#' @title Prepare folder structure
#'
#' @description \code{makeDirectories} creates a set of folders to store all the necessary
#' data for the OPSCR ROvQuant analysis. 
#' 
#' @param path A path indicating where to create the different folders.
#' @param subFolders A \code{character vector} with the names of the subfolders to be created for the nimble input and output files.
#' @param show.dir A \code{logical}. (Optional; requires package  \emph{fs}) If TRUE, the directory structure created is displayed in a tree-like format.
#'
#' @return boolean invisible(FALSE) if nothing was created, invisible(TRUE) if folders were created in \emph{path}.
#' 
#' @examples \dontrun{makeDirectories()}
#' 
#' @author Pierre Dupont
#' 
#' @rdname makeDirectories
#' @export
makeDirectories <- function( path = NULL,
                             subFolders = NULL,
                             show.dir = TRUE){
  ##-- Check if the directory already exists 
  if (is.null(path)) {path <- getwd()}
  if (dir.exists(path)) {
    split <- unlist(strsplit(path, split = "/"))
    message(paste0("A folder named '", split[length(split)], "' already exists in the specified directory:"))
    message(paste(split[-length(split)], collapse = "/"))
  } else {
    dir.create(path, showWarnings = F, recursive = T)
    }

  ##-- Set-up directory structure
  if (!is.null(subFolders)) {
   for (f in subFolders) {
     dir.create( file.path(path, "nimbleInFiles", f), showWarnings = F, recursive = T)
     dir.create( file.path(path, "nimbleInFiles", f), showWarnings = F, recursive = T)
     dir.create( file.path(path, "nimbleOutFiles", f), showWarnings = F, recursive = T)
     dir.create( file.path(path, "nimbleOutFiles", f), showWarnings = F, recursive = T)
   }#f
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
  if (show.dir) {
    if (requireNamespace("fs", quietly = TRUE)) {
    message("\n The following directory structure was created for this analysis:\n")
    fs::dir_tree(path)
    } else {
      warning("The 'fs' package must be installed to use this functionality.")
    }
  }
}
