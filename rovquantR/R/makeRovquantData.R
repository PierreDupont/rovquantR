#' @title Data preparation.
#'
#' @description
#' \code{makeRovquantData} calls a custom Rmarkdown template that identifies 
#' and loads the most recent Rovbase data available for the specified species 
#' and performs the OPSCR data preparation for model fitting. 
#' The data preparation is specific to each species and incorporates the
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
#'  as prepared by \code{\link{cleanRovBaseData}}.
#' @param working_dir A \code{path} to the directory for this analysis containing
#'  the \code{nimbleInputFiles} folder to store the prepared data.  
#' @param keep_dead A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
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
#' @examples
#' makeRovquantData_bear(2012:2021)
#' 
#' 
#' @rdname makeRovquantData_bear
#' @export
#' 
makeRovquantData <- function(
    ##-- paths
  species = c("bear","wolf","wolverine"),
  data_dir = "./Data",
  working_dir = NULL,
  
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
  plot.check = FALSE,
  print.report = TRUE
) {
  
  ##-- Check species and set corresponding sampling period
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
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
      resize.factor = 1,
      
      ##-- miscellanious
      plot.check = FALSE,
      print.report = FALSE)
  }
  # if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolf(...)
  # }
  # if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
  #   out <- makeRovquantData_wolverine(...)
  # }
  
  return(out)
}