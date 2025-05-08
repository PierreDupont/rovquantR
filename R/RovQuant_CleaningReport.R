#' @title Data set clean-up.
#'
#' @description
#' \code{cleanRovbaseData} calls a custom Rmarkdown template that identifies and loads
#' the most recent RovBase data available for the specified species in the specified
#'  \code{data.dir} location, and conducts a set of data cleaning steps that include:
#'  \itemize{
#'  \item{"parameter 1"}{removing un-identified samples}
#'  \item{"parameter 2"}{checking sex-assignment}
#'  \item{"parameter 3"}{removing samples flagged as unusable by RovData/Mike}
#' }
#'  
#'  Additionally, it can produce a \code{html} report describing the content of the data in terms of number of samples, individuals, etc... 
#'
#' @name cleanRovbaseData
#' 
#' @param data.dir the \code{path} pointing to the directory containing the raw 
#' data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. By default,
#'  the cleaned data will be stored in a subfolder of this working directory called 'data'.
#' @param species A \code{character} string with the name of the focal species
#'  ("bear", "wolf", or "wolverine").
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param sex A \code{character} vector containing the sex of interest. 
#' Can be "Hunn" for females or "Hann" for males. Default is both sexes (c("Hunn","Hann")).
#' @param sampling.months (Optional) A \code{list} containing the sampling period months.
#' If the sampling period overlaps two calendar years, the list should contain one element per year.
#' (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param rename.list (Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.
#' @param Rmd.template (Optional) The \code{path} to the \code{.rmd} template to be used for
#'  cleaning the data. By default, the \code{.rmd} template provided with the 
#'  \code{rovquantR} package is used.  
#' @param overwrite A \code{logical} (default = FALSE) to force ovewriting of previously existing clean data.
#'  If FALSE, the function checks for any pre-existing clean data files and ask whether to overwrite it or not.
#' @param output.dir (Optional) the \code{path} pointing to the directory where the \code{.html} report will be printed.
#' By default, the \code{.html} report describing the content of the clean data will 
#' be placed in a subfolder of the working directory (\code{working.dir}) called 'reports'.

#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects 
#' for the species and period specified. The clean data file is saved as an \code{.RData} 
#' file named using the species name and the date of extraction of the raw Rovbase data 
#' to facilitate replicability (e.g. 'CleanData_bear_2024-08-10.RData').
#' \item A \code{.html} report summarizing the data cleaning process. 
#' The \code{.RData} report is using the same naming convention as the clean \code{.RData} (e.g. 'CleanData_bear_2024-08-10.html').
#' \item Additional \code{.png} images that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom rmarkdown render
#' @import ggplot2 
#' @import kableExtra 

NULL
#' @rdname cleanRovbaseData
#' @export
cleanRovbaseData <- function(
    ##-- paths
  data.dir,
  working.dir,
  
  ##-- data
  species = c("bear","wolf","wolverine"),
  years = NULL, 
  sex = c("female","male"),
  sampling.months = NULL,
  rename.list = NULL,
  
  ##-- miscellanious
  Rmd.template = NULL,
  overwrite = FALSE,
  output.dir = NULL
) {
  
  ##----- 1. INITIAL CHECKS -----
  
  ##-- Make sure directory structure exists
  makeDirectories( path = working.dir,
                   subFolders = sex,
                   show.dir = TRUE)
  
  ##-- Years
  if(is.null(years)){ years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }
  
  ##-- Renaming list
  if(is.null(rename.list)){
    rename.list = c(
      Age_estimated = "Alder, vurdert",
      Age = "Alder, verifisert",
      Age_verif_by = "Alder, verifisert av",
      Age_class = "Alder på dødt individ",
      Age_class_verif = "Aldersklasse verifisert SVA",
      Analyzed_by = "AnalysertAv",
      Analysis_priority = "Analyseprioritet",
      Approved_by = "Godkjent av",
      Approved_date = "Godkjentdato",
      Assessment = "Vurdering",
      Barcode_sample = "Strekkode (Prøve)",
      Barcode = "Strekkode (Analyse)",
      Birth_territory = "Født revir",
      CITES = "CITES-nummer",
      Collected_by = "Hvem samlet inn",
      Collector_name = "Samlet selv - Navn",
      Collector_phone = "Samlet selv - Telefon",
      Collector_email = "Samlet selv - E-post",
      Collector_role = "Samlet selv - Rolle",
      Collector_other_name = "Annen innsamler - Navn" ,
      Collector_other_phone = "Annen innsamler - Telefon",
      Collector_other_email = "Annen innsamler - E-post",
      Collector_other_role = "Annen innsamler - Rolle",
      Comments_sample = "Merknad (Prøve)",
      Comments = "Merknad (Analyse)",
      Control_status = "Kontrollstatus",
      Coordinate_system = "Koordinatsystem",
      Counted_off_against_decision = "Regnes av mot vedtak",
      County_number = "Fylkenummer",
      County = "Fylke",
      Date = "Funnetdato",
      Date = "Dødsdato",
      Death_cause = "Bakgrunn/årsak",
      Death_method = "Bakgrunn/årsak metode",
      Death_purpose = "Bakgrunn/årsak formål",
      DNAID_sample = "DNAID (Prøve)",
      DNAID = "DNAID (Analyse)",
      EventID = "HendelseID",
      East_Original = "Øst (opprinnelig)",
      East_RT90 = "Øst (RT90)",
      East_UTM33 = "Øst (UTM33/SWEREF99 TM)",
      Felling_site_verif = "Kontroll av fellingsted",
      Field_personnel ="Feltpersonell",
      Hunting_date = "Observasjons/Jaktdato",
      Id = "Individ",
      Juvenile = "Yngling",
      Mountain_area = "Fjellområde",
      Method = "Metode",
      Municipality_number = "Kommunenummer",
      Municipality = "Kommune",
      North_original = "Nord (opprinnelig)",
      North_RT90 = "Nord (RT90)",
      North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
      Origin = "Opprinnelse",
      Outcome = "Utfall",
      Last_saved_by_sample = "Sist lagret av (Prøve)",
      Last_saved_sample = "Sist lagret dato (Prøve)",
      Last_saved_by = "Sist lagret av (Analyse)",
      Last_saved = "Sist lagret dato (Analyse)",
      Last_saved_by = "Sist lagret av",
      Last_saved =  "Sist lagret dato",
      Locality = "Lokalitet",
      Location = "Funnsted",
      Lansstyrelsen_number = "Länsstyrelsens nr",
      Quality_checked = "Kvalitetssikret av feltpersonell",
      Quality_check_name = "Kvalitetssikrer - navn",
      Quality_check_orga = "Kvalitetssikrer - Organisasjon",
      Release_Date = "Frigivelsesdato",
      Sample_type = "Prøvetype",
      Sensitivity = "Følsomhet",
      Species_sample = "Art (Prøve)",
      Site_quality = "Stedkvalitet",
      Time_of_death = "Dødstidspunkt",
      Tips_name = "Tipser - Navn",
      Tips_phone = "Tipser - Telefon",
      Tips_email = "Tipser - E-post",
      Tips_role = "Tipser - Rolle",
      Tissue_sample = "Vevsprøve tatt",
      Release_Date = "Frigivelsesdato",
      RovbaseID = "RovbaseID (Analyse)",
      RovbaseID_sample = "RovbaseID (Prøve)",
      Species = "Art (Analyse)",
      Species = "Art",
      Sample_status = "Prøvestatus",
      Sensitivity = "Følsomhet",
      Sex_analysis = "Kjønn (Analyse)",
      Sex = "Kjønn (Individ)",
      Sex = "Kjønn",
      Site_quality = "Stedkvalitet",
      SVAID = "SVAID",
      Uncertain_date = "Usikker dødsdato",
      Weight_slaughter = "Slaktevekt",
      Weight_total =  "Helvekt")
  }
  
  ##-- Extract the date from the last .xlsx data file
  DATE <- getMostRecent( path = data.dir, pattern = "DNA")
  
  ##-- Species
  if (length(species)>1){
    stop('This function can only deal with one species at a time... \nPlease, use one of "bear", "wolf", or "wolverine" for the n target species.')
  }
  
  if(sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0){
    species <- "bear"
  } else {
    if(sum(grep("wolf", species, ignore.case = T))>0|
       sum(grep("ulv", species, ignore.case = T))>0){
      species <- "wolf"
    } else {
      if(sum(grep("wolverine", species, ignore.case = T))>0|
         sum(grep("jerv", species, ignore.case = T))>0){
        species <- "Wolverine"
      }
    }
  }
  
  
  ##-- Set file name for clean data
  fileName <- paste0("CleanData_", species, "_", DATE, ".RData")
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  if (!overwrite) {
    existTest <- file.exists(file.path(working.dir, "data", fileName))
    if (any(existTest)) {
      message(paste0("A file named '", fileName[existTest], "' already exists in: \n",
                     file.path(working.dir, "data")))
      message("Are you sure you want to proceed and overwrite existing clean data? (y/n) ")
      question1 <- readLines(n = 1)
      if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
        message("Not overwriting existing files...")
        return(invisible(NULL))
      } else {
        message(paste0("Now overwriting '", fileName[existTest],"'.\n"))
      }
    }
  }
  
  
  
  ##----- 2. CLEAN THE DATA -----
  
  ##-----   2.1. CLEAN BEAR DATA -----
  
  ##-- Check species and use corresponding function
  if("bear" %in% species){
    
    ##-- Prepare the data
    out <- cleanRovbaseData_bear(
      ##-- paths
      data.dir,
      working.dir,
      ##-- data
      years,
      sex,
      sampling.months,
      rename.list)
  }
  
  
  
  ##-----   2.2. CLEAN WOLF DATA -----
  
  ##-- Check species and use corresponding function
  if(species %in% "wolf"){
    
    ##-- Prepare the data
    out <- cleanRovbaseData_wolf(
      ##-- paths
      data.dir,
      working.dir,
      ##-- data
      years,
      sex,
      sampling.months, 
      rename.list)
  }
  
  
  
  ##-----   2.3. CLEAN WOLVERINE DATA -----
  
  ##-- Check species and use corresponding function
  if(species %in% "wolverine"){
    
    ##-- Prepare the data
    out <- cleanRovbaseData_wolverine(
      ##-- paths
      data.dir,
      working.dir,
      ##-- data
      years,
      sex,
      sampling.months, 
      rename.list)
  }
  
  
  
  ##---- 3. PRINT THE REPORT -----
  
  if(print.report){
    
    ##-- Find the .rmd template for the report.
    if(is.null(Rmd.template)){
      Rmd.template <- system.file("rmd", "RovQuant_CleaningReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop('Can not find a .rmd template called "RovQuant_CleaningReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
      } 
    }
    
    ##-- Find the directory to print the report.
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    
    ##-- Render .Rmd report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = out$SPECIES,
                     years = out$YEARS,
                     sex = out$SEX,
                     date = out$DATE,
                     info = out$INFO,
                     working.dir = working.dir),
      output_dir = output.dir,
      output_file = paste0("CleanData_", out$SPECIES, "_", out$DATE,".html"))
  }
}
