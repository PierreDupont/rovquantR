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
#'  the \code{nimbleInputFiles} folder to store the prepared data.  
#' @param years A \code{list}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param sex A \code{character} denoting the species; can be any of "bear", "wolf" or "wolverine.
#' @param aug.factor A \code{numeric} to the directory containing the clean Rovbase data,
#' @param sampling.months A \code{list}.
#' @param habitat.res A \code{numeric}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param buffer.size A \code{numeric} denoting the species; can be any of "bear", "wolf" or "wolverine.
#' @param max.move.dist A \code{numeric}.
#' @param detector.res A \code{numeric} to the directory for this analysis containing
#' @param subdetector.res A \code{numeric}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param max.det.dist A \code{numeric} denoting the species; can be any of "bear", "wolf" or "wolverine.
#' @param resize.factor A \code{numeric}.
#' @param print.report A \code{logical}.
#' @param Rmd.template A \code{path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
#' @param output.dir A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' 
#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects
#'  for the species and period specified.
#' \item A \code{.html} report summarizing the data cleaning process. 
#' \item Additional \code{.png} images that can be reused somewhere else.
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
  species = c("bear","wolf","wolverine"),
  data.dir = getwd(),
  working.dir = getwd(),
  
  ##-- data
  years = NULL,
  sex = c("Hann","Hunn"),
  aug.factor = 2,
  sampling.months = NULL,
  
  ##-- habitat
  habitat.res = NULL, 
  buffer.size = NULL,
  max.move.dist = NULL,
  
  ##-- detectors
  detector.res = NULL,
  subdetector.res = NULL,
  max.det.dist = NULL,
  resize.factor = 1,
  
  ##-- miscellanious
  print.report = TRUE, 
  Rmd.template = NULL,
  output.dir = NULL
) {
  
  ##---- 1. BROWN BEAR DATA PREPARATION -----
  
  ##-- Check species and use corresponding function
  if(sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0){
  
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
  #
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
  #
  # ##-- Check species and use corresponding function
  # if(sum(grep("wolverine", species, ignore.case = T))>0|
  #    sum(grep("jerv", species, ignore.case = T))>0|
  #    sum(grep("järv", species, ignore.case = T))>0){
  #   
  #   ##-- Get clean name for the report
  #   SPECIES <- "Wolverine"
  #   
  #   ##-- Extract date from the last cleaned data file
  #   DATE <- getMostRecent( 
  #     path = file.path(working.dir,"data"),
  #     pattern = "CleanData_wolverine")
  # 
  #   ##-- Prepare the data
  #   out <- makeRovquantData_wolverine(
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
  
  
  
  ##---- 4. PRINT REPORT -----
  
  if(print.report){
    
    ##-- Find the .rmd template for the report.
    if(is.null(Rmd.template)){
      Rmd.template <- system.file("rmd", "RovQuant_DataReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop('Can not find a .rmd template called "RovQuant_DataReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
      } 
    }
    
    
    ##-- Find the directory to print the report.
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    
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