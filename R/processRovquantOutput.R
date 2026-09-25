#' @title RovQuant model output processing
#'
#' @description
#' \code{processRovquantOutput_bear} calls a custom Rmarkdown template that combines 
#' and processes MCMC outputs from NIMBLE models and produces figures,
#' tables and rasters of interest (e.g. population density maps)
#' 
#' @param working.dir \code{Path} to the directory containing the \code{nimbleOutFiles} folder where MCMC outputs are stored. 
#' @param species \code{Character} string denoting the species; can be one of "bear", "wolf" or "wolverine.
#' @param nburnin Number of initial, pre-thinning, MCMC bites to discard. Default value is 0.
#' @param thin Thinning interval for collecting MCMC samples, corresponding to monitors. Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' @param thin2 Thinning interval for collecting MCMC samples, corresponding to the second, optional set of monitors2. Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' @param niter Number of MCMC iterations to be used for density extraction. Default value is 100.
#' @param extraction.res Resolution (in meters) to use for density extraction. Default value is 5000m.
#' @param years Which years to use for density extraction.
#' @param print.report \code{Logical}. Whether an .html report summarizing the results should be printed (TRUE) or not (FALSE).
#' @param Rmd.template (optional) \code{Path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
#' @param output.dir (optional) \code{Path} to the location of the final .html report. Default is to print in the 'reports' folder of the working directory.
#' @param overwrite \code{logical} Whether to silently overwrite (TRUE) or ask before overwriting potentially existing .html reports (FALSE).
#'
#' @return This function returns:
#' \enumerate{
#' \item Multiple \code{.RData} files with the processed MCMC outputs and density outputs.
#' \item A \code{.html} report summarizing the data cleaning process. 
#' \item Additional \code{.png} images and \code{.csv} files that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom  rmarkdown render
#' 
#' @rdname processRovquantOutput
#' @export
processRovquantOutput <- function(
  ##-- paths
  working.dir = NULL,
  
  ##-- MCMC processing
  species = c("bear","wolf","wolverine"),
  nburnin = 0,
  thin = 1,
  thin2 = 1,
  
  ##-- Density extraction
  niter = 100,
  extraction.res = 5000,
  years = NULL,
  
  ##-- miscellanious
  print.report = TRUE,
  Rmd.template = NULL,
  output.dir = NULL,
  overwrite = FALSE, 
  full.report = FALSE
) {
  
  ##---- 1. BROWN BEAR RESULTS PROCESSING -----
  
  if(sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0){
    
    ##-- Process the model output
    out <- processRovquantOutput_bear(
      working.dir = working.dir,
      nburnin = nburnin,
      niter = niter,
      thin = thin, 
      thin2 = thin2,
      extraction.res = extraction.res,
      years = years,
      overwrite = overwrite)
  }
  
  

  ##---- 2. WOLF RESULTS PROCESSING -----
  
  if(sum(grep("wolf", species, ignore.case = T))>0|
     sum(grep("ulv", species, ignore.case = T))>0){

    ##-- Process the model output
    out <- processRovquantOutput_wolf(
      working.dir = working.dir,
      nburnin = nburnin,
      niter = niter,
      thin = thin,
      thin2 = thin2,
      extraction.res = extraction.res,
      years = years,
      overwrite = overwrite)
  }

  
  
  ##---- 3. WOLVERINE RESULTS PROCESSING -----
  
  if(sum(grep("wolverine", species, ignore.case = T))>0|
     sum(grep("jerv", species, ignore.case = T))>0){
    
    ##-- Process the model output
    out <- processRovquantOutput_wolverine_SCR(
      working.dir = working.dir,
      nburnin = nburnin,
      niter = niter,
      thin = thin,
      thin2 = thin2,
      extraction.res = extraction.res,
      years = years,
      overwrite = overwrite)
  }
  
  

  ##---- 4. PRINT REPORT -----
  
  if(print.report){
    
    ##-- Find the correct .rmd template for the report.
    if(is.null(Rmd.template)){
      if(full.report){
        Rmd.template <- system.file("rmd", "RovQuant_FullReport.Rmd", package = "rovquantR")
        if(!file.exists(Rmd.template)) {
          stop('Can not find the default .rmd template called "RovQuant_FullReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
        } 
      } else {
        if(sum(grep("wolverine", species, ignore.case = T))>0|
           sum(grep("jerv", species, ignore.case = T))>0){
          Rmd.template <- system.file("rmd", "RovQuant_OutputReport_SCR.Rmd", package = "rovquantR")
        } else {
          Rmd.template <- system.file("rmd", "RovQuant_OutputReport.Rmd", package = "rovquantR")
        }
        if(!file.exists(Rmd.template)) {
          stop('Can not find the .rmd template called "RovQuant_OutputReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
        } 
      }
    }
    
    ##-- Find the directory to print the report.
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = out$SPECIES,
                     years = out$YEARS,
                     seasons = out$SEASONS,
                     date = out$DATE,
                     working.dir = working.dir),
      output_dir = output.dir,
      output_file = paste0("Results_", out$engSpecies, "_", out$DATE,".html"))
  }
}