#' @title Data set clean-up.
#'
#' @description
#' \code{cleanRovbaseData} calls a custom Rmarkdown template that identifies and loads
#'  the most recent RovBase data available on dropbox for the specified species
#'  and conducts a set of data cleaning steps that include:
#'  - removing un-identified samples
#'  - checking sex-assignment
#'  - removing flagged samples from RovData/Mike
#'
#' @name cleanRovbaseData
#' 
#' @param species A \code{character} string with the name of the focal species
#'  ("bear", "wolf", or "wolverine").
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param samplingMonths A \code{list} containing the sampling period months.
#' If the sampling period overlaps two calendar years, the list should contain one element per year
#' (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param data_dir the \code{path} pointing to the directory containing the raw 
#' data from Rovbase.
#' @param output_dir the \code{path} pointing to the directory where the cleaned 
#' data will be stored. The clean data, as well as a \code{.html} report describing 
#' the content of the clean data, will be placed in a folder with the species name,
#' inside a subfolder with the date of extraction of the raw Rovbase data.
#' @param Rmd_template the \code{path} to the \code{.rmd} template to be used for
#'  cleaning the data. By default, the \code{.rmd} template provided with the 
#'  \code{rovquantR} package is used.  
#' @param overwrite A \code{logical} Whether previously existing clean data should
#' be overwritten or not (default = FALSE).
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
NULL
#' @rdname cleanRovbaseData
#' @export
cleanRovbaseData <- function( species,
                              years = NULL, 
                              samplingMonths = NULL,
                              data_dir = "./Data",
                              output_dir = "./Data",
                              Rmd_template = NULL,
                              overwrite = FALSE)
  {
  ##-- Check if years are provided (if not, uses the period 1990 until now)
  if(is.null(years)){ years <- 1990:as.numeric(format(Sys.Date(), "%Y")) }
  
  
  ##-- BROWN BEARS 
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
    
    ##-- Call bear-specific data cleaning function 
    cleanRovbaseData_bear( years,
                           samplingMonths,
                           data_dir,
                           output_dir,
                           Rmd_template,
                           overwrite = FALSE)
  }

  
  # ##-- WOLVES 
  # if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
  #   
  #   ##-- Set monitoring season for the bear
  #   if(is.null(samplingMonths)){ samplingMonths <- list(10:12,1:4) }
  #   
  #   ##-- Call bear-specific data cleaning function 
  #   cleanRovbaseData_wolf( years,
  #                          samplingMonths,
  #                          data_dir,
  #                          output_dir,
  #                          Rmd_template,
  #                          overwrite = FALSE)  
  # }

  
  # ##-- WOLVERINES 
  # if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
  #   
  #   ##-- Set monitoring season for the bear
  #   if(is.null(samplingMonths)){ samplingMonths <- list(12,1:5) }
  #   
  #   ##-- Call bear-specific data cleaning function 
  #   cleanRovbaseData_wolverine( years,
  #                               samplingMonths,
  #                               data_dir,
  #                               output_dir,
  #                               Rmd_template,
  #                               overwrite = FALSE) 
  # }
  
 
  ##-- Find the .rmd template used to clean the data and print out the report.
  if(is.null(Rmd_template)){
    Rmd_template <- system.file("rmd", "RovQuant_CleaningReport.Rmd", package = "rovquantR")
    if(!file.exists(Rmd_template)) {
      stop('Can not find the Rmarkdown document to use for cleaning Rovbase.3.0 data.\n You must provide the path to the Rmarkdown template through the "Rmd_template" argument.')
    } 
  }
  
  
  ##-- Check output directory 
  output_folder <- file.path(output_dir, SPECIES, DATE)

  
  ##-- Check that the "report" does not already exist to avoid overwriting
  if(file.exists(file.path(output_folder, paste0("Data_", SPECIES, "_", DATE,".html")))){
    message(paste0("A file named 'Data_", SPECIES, "_", DATE,
                   ".html' already exists in the specified directory."))
    message("Are you sure you want to proceed and overwrite existing clean data? (y/n) ")
    question1 <- readLines(n = 1)
    if(regexpr(question1, 'y', ignore.case = TRUE) != 1){
      message("Not overwriting existing files...")
    } else {
      message(paste0("Now overwriting 'Data_", SPECIES, "_", DATE,".html'."))
      
      ##-- Clean the data and print report
      rmarkdown::render(
        input = Rmd_template,
        params = list( species = SPECIES,
                       years = years,
                       samplingMonths = SP,
                       dir.in = data_dir,
                       dir.out = output_folder,
                       modDate = DATE),
        output_dir = output_folder,
        output_file = paste0("Data_", SPECIES, "_", DATE,".html"))
    }
  } else {
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd_template,
      params = list( species = SPECIES,
                     years = years,
                     samplingMonths = SP,
                     dir.in = data_dir,
                     dir.out = output_folder,
                     modDate = DATE),
      output_dir = output_folder,
      output_file = paste0("Data_", SPECIES, "_", DATE,".html"))
  }
}