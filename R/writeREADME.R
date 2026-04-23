#' @title Write README file for the analysis
#'
#' @description \code{writeREADME} creates and opens a README.rmd file to document the analysis.
#' 
#' @param path A \code{path} indicating where to create the README file.
#' @param template (Optional) A \code{path} indicating where to find the template of the README file. Default is to use the template provided in the \code{rovquantR} package.
#' @param A \code{logical}. (Optional; requires package  \emph{fs}) If TRUE, the directory structure created is displayed in a tree-like format.
#'
#' @return boolean invisible(FALSE) if nothing was created, invisible(TRUE) if the README file was created created in \emph{path}.
#' 
#' @examples \dontrun{writeREADME()}
#' 
#' @author Pierre Dupont
#' 
#' @importFrom whisker whisker.render 
#' 
#' @rdname writeREADME
#' @export
writeREADME <- function(
    path = NULL,
    template = NULL,
    name = "rovquantR analysis",
    save_as = "README.Rmd",
    overwrite = FALSE){
  if (is.null(path)) {path <- getwd()}
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  if(!overwrite) {
    existTest <- file.exists(file.path(path,save_as))
    if (any(existTest)) {
      message(paste0("A file named '", save_as, "' already exists in the specified directory:"))
      message(path)
      message("Are you sure you want to proceed and overwrite existing file? (y/n) ")
      question1 <- readLines(n = 1)
      if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
        message("Not overwriting existing file...")
        return(invisible(FALSE))
      } else {
        message(paste0("Now overwriting '", save_as,"'.\n"))
      }
    }
  }
  
  ##-- Find the template for the README.rmd file
  if(is.null(template)) {
    template_path <- system.file("rmd", "template-README", package = "rovquantR")
    if(!file.exists(template_path)) {
      stop("Can not find the template for the README file.\n You must provide the path to the README template through the \"Rmd_template\" argument.")
    }
  }
  
  ##-- Create list of data used to render the template
  data <- list( Analysis = name,
                Rmd = TRUE)
  
  ##-- Render the template with custom data list
  template_contents <- strsplit(whisker::whisker.render(readLines(template_path, encoding = 'UTF-8'),
                                                        data),
                                "\n")[[1]]

  ##-- convert embedded newlines
  lines <- gsub("\r?\n", "\r\n", template_contents)
  
  ##-- write file
  writeLines( text = enc2utf8(lines),
              con = file.path(path, save_as),
              sep = "\r\n",
              useBytes = TRUE)
  
  ##-- Open file for editing if possible
  if(Sys.getenv("RSTUDIO")  == "1"){
    rstudioapi::navigateToFile(file.path(path, save_as))
  }
  
  invisible(TRUE)
}









