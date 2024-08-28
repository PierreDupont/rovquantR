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
#' @param species A \code{character} string with the name of the focal species
#' @param dead_recoveries A \code{DataFrame} object containing the raw dead recoveries data.
#' @param species_id A \code{Numeric} object containing the id of the focal species (1=Bears;2=Wolverines;3:Wolves).
#' @param country_polygon A \code{SpatialPointsDataFrame} with Country polygon for correct assignment of samples to countries
#' @param threshold_month A \code{Numeric} with initial month of the biological year: 1:January...12=December. all samples with months<threshold get year-1. so they get similar year.  
#' @param keep_dead A \code{logical}  Whether dead recovery should be included (TRUE) or not(FALSE)
#' 
#' @return 
#' A \code{.RData} file with the clean NGS and dead recovery data objects
#' for the species and period specified.
#' A \code{html} report summarizing the data cleaning process
#' Additional \code{.png} images that can be reused somewhere else.
#'
#' @author Pierre Dupont
#' 
#' @examples
#' cleanRovbaseData("bear", 2012:2021)
#' cleanRovbaseData("Jerv", 2020)
#' cleanRovbaseData("Gray wolf", 2018:2022)
#' 
#' @importFrom rmarkdown::render
#' 
#' @rdname cleanRovbaseData
#' @export
#' 
cleanRovbaseData <- function( species,
                              years = NULL, 
                              data_dir = "./Data",
                              Rmd_template = NULL,
                              output_dir = "./Data",
                              overwrite = FALSE){
  ##-- Check if years are provided -----
  ##-- (if not, uses the period 1990 until now)
  if(is.null(years)){
    years <- 1990:as.numeric(format(Sys.Date(), "%Y"))
  }
  
  
  ##-- Check species and set corresponding sampling period
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
    SPECIES <- "bear"
    SP <- list(4:11)
  }
  if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
    SPECIES <- "wolf"
    SP <- list(10:12,1:4)
  }
  if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
    SPECIES <- "wolverine"
    SP <- list(12,1:5)
  }
  
  
  ##-- Extract the date from the last .csv data file
  DATE <- getMostRecent( path = data_dir,
                         pattern = paste0("_",SPECIES,".csv"))

  
  ##-- Find the .rmd template used to clean the data and print out the report.
  if(is.null(Rmd_template)){
    Rmd_template <- system.file("rmd", "RovBase_DataCleaning.Rmd", package = "rovquantR")
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