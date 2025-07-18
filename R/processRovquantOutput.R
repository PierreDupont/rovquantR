#' @title RovQuant OPSCR output processing
#'
#' @description
#' The \code{processRovquantOutput} function calls a custom Rmarkdown template that identifies 
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
#' @param data.dir A \code{path} to the directory containing the clean Rovbase data,
#'  as prepared by \code{cleanRovBaseData}.
#' @param working.dir A \code{path} to the directory for this analysis containing
#'  the \code{nimbleInputFiles} folder to store the prepared data.  
#' @param nburnin An \code{integer} denoting the number of MCMC  dead recovery should be included (TRUE) or not (FALSE)
#' @param niter A \code{logical} Whether dead recovery should be included (TRUE) or not (FALSE)
#' @param extraction.res A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param print.report A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param output_dir A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' @param Rmd_template A \code{path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
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
#' @rdname processRovquantOutput
#' @export
processRovquantOutput <- function(
  ##-- paths
  data.dir = "./Data",
  working.dir = NULL,
  
  ##-- MCMC processing
  species = c("bear","wolf","wolverine"),
  nburnin = 5,
  
  ##-- Density extraction
  niter = 100,
  extraction.res = 5000,
  
  ##-- miscellanious
  print.report = TRUE,
  Rmd.template = NULL,
  output.dir = NULL,
  overwrite = FALSE
) {
  
  ##---- 1. BROWN BEAR RESULTS PROCESSING -----
  
  if(sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0){
    
    ##-- Process the model output
    out <- processRovquantOutput_bear(
      data.dir,
      working.dir,
      nburnin,
      niter,
      extraction.res,
      overwrite)
  }
  
  

  ##---- 2. WOLF RESULTS PROCESSING -----
  
  if(sum(grep("wolf", species, ignore.case = T))>0|
     sum(grep("ulv", species, ignore.case = T))>0){

    ##-- Process the model output
    out <- processRovquantOutput_wolf(
      data.dir,
      working.dir,
      nburnin,
      niter,
      extraction.res,
      overwrite)
  }

  
  
  ##---- 3. WOLVERINE RESULTS PROCESSING -----
  
  if(sum(grep("wolverine", species, ignore.case = T))>0|
     sum(grep("jerv", species, ignore.case = T))>0){
    
    ##-- Process the model output
    out <- processRovquantOutput_wolverine(
      data.dir,
      working.dir,
      nburnin,
      niter,
      extraction.res,
      overwrite)
  }
  
  

  ##---- 4. PRINT REPORT -----
  
  if(print.report){
    
    ##-- Find the .rmd template for the report.
    if(is.null(Rmd.template)){
      Rmd.template <- system.file("rmd", "RovQuant_OutputReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop('Can not find a .rmd template called "RovQuant_OutputReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
      } 
    }
    
    ##-- Find the directory to print the report.
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = out$SPECIES,
                     years = out$YEARS,
                     date = out$DATE,
                     working.dir = working.dir),
      output_dir = output.dir,
      output_file = paste0("Results_", out$engSpecies, "_", out$DATE,".html"))
  }
}