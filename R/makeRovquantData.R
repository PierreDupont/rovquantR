#' @title RovQuant OPSCR data preparation.
#'
#' @description
#' The \code{makeRovquantData} function calls a custom Rmarkdown template that identifies 
#' and loads the most recent Rovbase data available for the specified species 
#' and performs the OPSCR data preparation for model fitting. 
#' 
#' The data preparation process is composed of three main steps:
#' \enumerate{
#' \item Defining and formatting habitat characteristics
#' \item Defining and formatting detectors characteristics
#' \item Defining and formatting individual detection histories
#' }
#' 
#' The data preparation process is specific to each species and incorporates the
#' different developments from project RovQuant to allow fitting large-scale
#' SCR and OPSCR models (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item Rescaling spatial coordinates to a grid for fast look-up assignment.
#' \item A local evaluation approach (see Milleret et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item A sparse matrix representation of the observation data to reduce the size of objects to be processed.
#' \item An indicator to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' @param species A \code{character} denoting the species; can be any of "bear", "wolf" or "wolverine.
#' @param data.dir A \code{path} to the directory containing the clean Rovbase data,
#'  as prepared by \code{cleanRovBaseData}.
#' @param working.dir A \code{path} to the directory for this analysis containing
#'  the \code{nimbleInFiles} folder to store the prepared data.  
#' @param years A \code{numeric} vector denoting which years of data will be used for this analysis.
#' @param sex A \code{character} denoting for which sex the data should be prepared ("female", "male" or both).
#' @param aug.factor A \code{numeric} denoting the augmentation factor.
#' @param sampling.months A \code{list} denoting the months of the monitoring period. Can accomodate sampling periods over two calendar years, e.g. 'list(c(10:12),c(1:3))' for a monitoring period extending from October to March (included).
#' @param habitat.res A \code{numeric} denoting the resolution (in meters) to be used for the habitat definition.
#' @param buffer.size A \code{numeric} denoting the size (in meters) of the buffer around the searched area.
#' @param max.move.dist A \code{numeric} denoting the maximum allowed distance for inter-annual movement.
#' @param detector.res A \code{numeric} denoting the resolution (in meters) to be used for the detectors definition.
#' @param subdetector.res A \code{numeric} denoting the resolution (in meters) to be used for the sub-detectors.
#' @param max.det.dist A \code{numeric} denoting the maximum allowed distance for intra-annual movement.
#' @param resize.factor A \code{numeric} denoting the aggregation factor to be used for the local evaluation.
#' @param print.report A \code{logical} denoting whether to print out a \code{.html} report summarizing the data preparation process or not.
#' @param Rmd.template (optional) A \code{path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
#' @param output.dir (optional) A \code{path} denoting where to print the .html report. Default is in a folder named 'reports' in the working directory.
#' @param overwrite A \code{logical} (default = FALSE) to force overwriting of previously existing data.
#'  If FALSE, the function checks for any pre-existing data files and asks whether to overwrite it or not.
#'  
#' @return This function returns:
#' \enumerate{
#' \item Multiple \code{.RData} nimble input files containing the input data, model code and initial values necessary for fitting the model in nimble (one input file per MCMC chain).
#' \item A \code{.html} report summarizing the data preparation process. 
#' \item Additional \code{.png} images and \code{.csv} tables of the figures and tables present in the report.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom rmarkdown render
#' 
#' @rdname makeRovquantData
#' @export
makeRovquantData <- function(
  ##-- paths
  data.dir = "./Data",
  working.dir = NULL,
  
  ##-- data
  species = c("bear","wolf","wolverine"),
  years = NULL,
  sex = c("female","male"),
  aug.factor = NULL,
  sampling.months = NULL,
  
  ##-- habitat
  habitat.res = NULL, 
  buffer.size = NULL,
  max.move.dist = NULL,
  
  ##-- detectors
  detector.res = NULL,
  subdetector.res = NULL,
  max.det.dist = NULL,
  resize.factor = NULL,
  
  ##-- miscellanious
  rename.list = NULL,
  print.report = TRUE,
  Rmd.template = NULL,
  output.dir = NULL,
  overwrite = FALSE
) {
  
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  if (!overwrite) {
    if (length(list.files( file.path(working.dir, "nimbleInFiles"),
                              recursive = TRUE,
                              pattern = "nimbleInput.+RData"))>0) {
      message(paste0("Nimble input files already exist in: \n",
                     file.path(working.dir, "nimbleInFiles")))
      message("Are you sure you want to proceed and overwrite existing input files? (y/n) ")
      question1 <- readLines(n = 1)
      if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
        message("Not overwriting existing files...")
        return(invisible(NULL))
      } else {
        message("Now overwriting existing nimble input files.\n")
      }
    }
  }

  
  ##---- 1. BROWN BEAR DATA PREPARATION -----
  
  ##-- Check species and use corresponding function
  if (sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0) {
  
    ##-- Prepare the data
    out <- makeRovquantData_bear(
      ##-- paths
      data.dir,
      working.dir,
      ##-- data
      years,
      sex,
      aug.factor,
      sampling.months,
      ##-- habitat
      habitat.res, 
      buffer.size,
      max.move.dist,
      ##-- detectors
      detector.res,
      subdetector.res,
      max.det.dist,
      resize.factor)
  }
  
  
  
  ##---- 2. WOLF DATA PREPARATION -----
  
  # ##-- Check species and use corresponding function
  # if(sum(grep("wolf", species, ignore.case = T))>0|
  #    sum(grep("wolves", species, ignore.case = T))>0|
  #    sum(grep("ulv", species, ignore.case = T))>0){
  #   
  #   ##-- Get clean name for the report
  #   SPECIES <- "Gray wolf"
  #   
  #   ##-- Extract date from the last cleaned data file
  #   DATE <- getMostRecent( 
  #     path = file.path(working.dir,"data"),
  #     pattern = "CleanData_wolf")
  #   
  #   ##-- Prepare the data
  #   out <- makeRovquantData_wolf(
  #     ##-- paths
  #     data.dir,
  #     working.dir,
  #     ##-- data
  #     years,
  #     sex,
  #     aug.factor,
  #     sampling.months,
  #     ##-- habitat
  #     habitat.res, 
  #     buffer.size,
  #     max.move.dist,
  #     ##-- detectors
  #     detector.res,
  #     subdetector.res,
  #     max.det.dist,
  #     resize.factor)
  # }

  
  
  ##---- 3. WOLVERINE DATA PREPARATION -----

  ##-- Check species and use corresponding function
  if (sum(grep("wolverine", species, ignore.case = T))>0|
     sum(grep("jerv", species, ignore.case = T))>0|
     sum(grep("järv", species, ignore.case = T))>0) {

    ##-- Get clean name for the report
    SPECIES <- "Wolverine"

    ##-- Extract date from the last cleaned data file
    DATE <- getMostRecent(
      path = file.path(working.dir,"data"),
      pattern = "CleanData_wolverine")

    ##-- Prepare the data
    out <- makeRovquantData_wolverine(
      ##-- paths
      data.dir,
      working.dir,
      ##-- data
      years,
      sex,
      aug.factor,
      sampling.months,
      ##-- habitat
      habitat.res,
      buffer.size,
      max.move.dist,
      ##-- detectors
      detector.res,
      subdetector.res,
      max.det.dist,
      resize.factor,
      renam.list)
  }
  
  
  
  ##---- 4. PRINT REPORT -----
  
  if (print.report) {
    
    ##-- Find the .rmd template for the report.
    if (is.null(Rmd.template)) {
      Rmd.template <- system.file("rmd", "RovQuant_DataReport.Rmd", package = "rovquantR")
      if (!file.exists(Rmd.template)) {
        stop('Can not find a .rmd template called "RovQuant_DataReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
      } 
    }
    
    ##-- Find the directory to print the report.
    if (is.null(output.dir)) { output.dir <- file.path(working.dir, "reports") }
    
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = out$SPECIES,
                     years = out$YEARS,
                     sex = out$SEX,
                     date = out$DATE,
                     working.dir = working.dir),
      output_dir = output.dir,
      output_file = paste0("ProcessedData_", out$engSpecies, "_", out$DATE,".html"))
  }
  
  }