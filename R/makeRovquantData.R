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
#' @param data_dir A \code{path} to the directory containing the clean Rovbase data,
#'  as prepared by \code{cleanRovBaseData}.
#' @param working_dir A \code{path} to the directory for this analysis containing
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
#' @param Rmd_template A \code{path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
#' @param output_dir A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
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
#' @rdname makeRovquantData
#' @export
makeRovquantData <- function(
  ##-- paths
  species = c("bear","wolf","wolverine"),
  data_dir = "./Data",
  working_dir,
  
  ##-- data
  years = NULL,
  sex = c("Hann","Hunn"),
  aug.factor = 2,
  sampling.months = list(4,5,6,7,8,9,10,11),
  
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 50000,
  max.move.dist = 250000,
  
  ##-- detectors
  detector.res = 5000,
  subdetector.res = 1000,
  max.det.dist = 70000,
  resize.factor = 1,
  
  ##-- miscellanious
  print.report = TRUE, 
  Rmd_template = NULL,
  output_dir = NULL
) {
  
  ##-- Check species and set corresponding sampling period
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
    
    SPECIES <- "Brown bear"
    
    ##-- Extract the date from the last cleaned data file
    DATE <- getMostRecent( 
      path = file.path(working_dir,"data"),
      pattern = "Data_bear")
    
    ##-- Prepare the data
    out <- makeRovquantData_bear(
      ##-- paths
      data_dir,
      working_dir,
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
      detector.res = 5000,
      subdetector.res = 1000,
      max.det.dist = 70000,
      resize.factor = 1)
  }
  
  # if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolf(...)
  # }
  
  # if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolverine(...)
  # }
  
  if(print.report){
    
    ##-- Find the .rmd template for the report.
    if(is.null(Rmd_template)){
      Rmd_template <- system.file("rmd", "RovQuant_DataReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd_template)) {
        stop('Can not find a .rmd template called "RovQuant_DataReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd_template" argument.')
      } 
    }
    
    
    ##-- Find the directory to print the report.
    if(is.null(output_dir)){
      output_dir <- working_dir
    }
    
    
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd_template,
      params = list( species = SPECIES,
                     years = years,
                     sex = sex,
                     working_dir = working_dir),
      output_dir = output_dir,
      output_file = paste0("Data_", SPECIES, "_", DATE,".html"))
  }
  
  return(out)
}