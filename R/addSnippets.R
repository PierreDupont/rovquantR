#' @title Add R snippets
#'
#' @description \code{addSnippets} copies all (missing) snippet definitions in 
#' \code{snippetsFile} to the RStudio's user snippet location.
#' (adapted for RovQuant from :
#' https://stackoverflow.com/questions/62220596/can-r-packages-add-code-snippets-to-users-snippet-files)
#'  
#' @param path A path indicating where to search for the snippet file.
#' @param snippetsFile A \code{character} object with the name of thefile containing the snippets definitions.
#'  
#' @return boolean invisible(FALSE) if nothing was added, invisible(TRUE) if snippet definitions were added
#' 
#' @importFrom rstudioapi versionInfo
#'
#' @examples \dontrun{addSnippets()}
#' 
#' @rdname addSnippets
#' @export
addSnippets <- function( path = getwd(),
                         snippetsFile = "snippets.R"){
  
  added <- FALSE
  
  ##-- If not on RStudio or RStudioServer exit
  if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
    return(NULL)
  }
  
  ##-- Path to directory for RStudios user files depends on OS
  if (rstudioapi::versionInfo()$version < "1.3") {
    rstudioSnippetsPath <- file.path(path.expand('~'),".R", "snippets")
  } else {
    if (.Platform$OS.type == "windows") {
      rstudioSnippetsPath <- file.path(Sys.getenv("APPDATA"), "RStudio", "snippets")
    } else {
      rstudioSnippetsPath <- file.path(path.expand('~'), ".config/rstudio", "snippets")
    }
  }
  
  ##-- Try to get template, if template is not found skip it
  snippetsFilesPath <- file.path(path, snippetsFile)#  system.file("rstudio", snippetsFile, lib.loc = path)
  if (snippetsFilesPath == "") {
    next()
  }
  
  ##-- load package snippets definitions
  snippetsFileContent <- readLines(snippetsFilesPath, warn = FALSE)
  
  ##-- Extract names of package snippets
  snippetsFileDefinitions <- snippetsFileContent[grepl("^snippet (.*)", snippetsFileContent)]
  
  ##-- Construct path for destination file
  rstudioSnippetsFilePath <- file.path(rstudioSnippetsPath, "r.snippets")
  
  ##-- If targeted RStudios user file does not exist, raise error (otherwise we would 'remove'
  ##-- the default snippets from the 'user file')
  if (!file.exists(rstudioSnippetsFilePath)) {
    stop(paste0( "'", rstudioSnippetsFilePath, "' does not exist yet\n.",
                 "Use RStudio -> Tools -> Global Options -> Code -> Edit Snippets\n",
                 "To initalize user defined snippets file by adding dummy snippet\n"))
  }
  
  ##-- Extract 'names' of already existing snitppets
  rstudioSnippetsFileContent <- readLines(rstudioSnippetsFilePath)
  rstudioSnippetDefinitions <- rstudioSnippetsFileContent[grepl("^snippet (.*)", rstudioSnippetsFileContent)]
  
  ##-- replace two spaces with tab, ONLY at beginning of string
  snippetsFileContentSanitized <- gsub("(?:^ {2})|\\G {2}|\\G\t", "\t", snippetsFileContent, perl = TRUE)
  
  ##-- find defintions appearing in custom snippets but not in rstudioSnippets
  ##-- if no snippets are missing go to next file
  snippetsToCopy <- setdiff(trimws(snippetsFileDefinitions), trimws(rstudioSnippetDefinitions))
  snippetsNotToCopy <- base::intersect(trimws(snippetsFileDefinitions), trimws(rstudioSnippetDefinitions))

  ##-- Inform user about changes, ask to confirm action
  if (interactive()) {
    cat(paste0("You are about to add the following ", length(snippetsToCopy),
               " snippets to '", rstudioSnippetsFilePath, "':\n",
               paste0(paste0("-", snippetsToCopy), collapse="\n")))
    if (length(snippetsNotToCopy) > 0) {
      cat(paste0("\n(The following snippets will NOT be added because there is already a snippet with that name:\n",
                 paste0(snippetsNotToCopy, collapse=", ") ,")"))
    }
    answer <- readline(prompt="Do you want to procedd (y/n): ")
    if (substr(answer, 1, 1) == "n") { next() }
  }
  
  ##-- Create list of line numbers where snippet definitons start
  ##-- This list is used to determine the end of each definition block
  allSnippetDefinitonStarts <- grep("^snippet .*", snippetsFileContentSanitized)
  
  for (s in snippetsToCopy) {
    startLine <- grep(paste0("^", s, ".*"), snippetsFileContentSanitized)
    
    ##-- Find last line of snippet definition:
    ##-- First find start of next defintion and return
    ##-- previous line number or lastline if already in last definiton
    endLine <- allSnippetDefinitonStarts[allSnippetDefinitonStarts > startLine][1] -1
    if(is.na(endLine)) {
      endLine <- length(snippetsFileContentSanitized)
    }
    
    snippetText <- paste0(snippetsFileContentSanitized[startLine:endLine], collapse = "\n")
    
    ##-- Make sure there is at least one empty line between entries
    if(tail(readLines(rstudioSnippetsFilePath), n=1) != "") {
      snippetText <- paste0("\n", snippetText)
    }
    
    ##-- Append snippet block, print message
    cat(paste0(snippetText, "\n"), file = rstudioSnippetsFilePath, append = TRUE)
    cat(paste0("* Added '", s, "' to '", rstudioSnippetsFilePath, "'\n"))
    added <- TRUE
  }


if (added) {
  cat("Restart RStudio to use new snippets")
}

return(invisible(added))

}  