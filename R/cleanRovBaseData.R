#' @title Data set clean-up.
#'
#' @description
#' \code{cleanRovbaseData} identifies and loads the most recent RovBase data available 
#' for the specified species in the specified \code{data.dir} location, and conducts a set of data cleaning steps that include:
#' \itemize{
#'  \item{removing un-identified samples}
#'  \item{checking sex-assignment}
#'  \item{removing samples flagged as unusable by RovData/Mike}
#' }
#'  
#' Additionally, it can produce a \code{html} report describing the content of the data in terms of number of samples, individuals, etc... 
#'
#' @name cleanRovbaseData
#' 
#' @param data.dir the \code{path} pointing to the directory containing the raw data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. By default, the cleaned data will be stored in a subfolder of this working directory called 'data'.
#' @param species A \code{character} string with the name of the focal species
#'  ("bear", "wolf", or "wolverine").
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param two.sex A \code{logical} determining whether the analysis will be done by sex (two.sex = T) or both together (two.sex = F).
#' @param sampling.months (Optional) A \code{list} containing the sampling period months. If the sampling period overlaps two calendar years, the list should contain one element per year (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param rename.list (Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.
#' @param print.report A \code{logical} denoting whether to print out a \code{.html} report summarizing the cleaning process or not.
#' @param Rmd.template (Optional) The \code{path} to the \code{.rmd} template to be used for cleaning the data. By default, the \code{.rmd} template provided with the \code{rovquantR} package is used.  
#' @param overwrite A \code{logical} (default = FALSE) to force overwriting of previously existing clean data.
#'  If FALSE, the function checks for any pre-existing clean data files and ask whether to overwrite it or not.
#' @param output.dir (Optional) the \code{path} pointing to the directory where the \code{.html} report will be printed.
#' By default, the \code{.html} report describing the content of the clean data will 
#' be placed in a subfolder of the working directory (\code{working.dir}) called 'reports'.
#' 
#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects 
#' for the species and period specified. The clean data file is saved as an \code{.RData} 
#' file named using the species name and the date of extraction of the raw Rovbase data 
#' to facilitate replicability (e.g. 'CleanData_bear_2024-08-10.RData').
#' \item A \code{.html} report summarizing the data cleaning process. 
#' The \code{.RData} report is using the same naming convention as the clean \code{.RData} (e.g. 'CleanData_bear_2024-08-10.html').
#' \item Additional \code{.png} images and summary \code{.csv} tables that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom rmarkdown render
#' @importFrom grDevices png
#' @importFrom graphics mtext 
#' @import sf 
#' @import ggplot2
#' @import patchwork
#' @import dplyr 
#' 
#' @rdname cleanRovbaseData
#' @export
cleanRovbaseData <- function(
  ##-- paths
  data.dir = "./Data",
  working.dir = NULL,
  
  ##-- data
  species = c("bear","wolf","wolverine"),
  years = NULL, 
  two.sex = TRUE,
  sampling.months = NULL,
  rename.list = NULL,
  
  ##-- miscellanious
  print.report = TRUE,
  Rmd.template = NULL,
  output.dir = NULL,
  overwrite = FALSE
) {
  
  ##----- 1. INITIAL CHECKS -----
  
  ##-- Make sure directory structure exists
  if(two.sex){
  makeDirectories( path = working.dir,
                   subFolders = c("female","male"),
                   show.dir = TRUE)
    } else {
    makeDirectories( path = working.dir,
                     show.dir = TRUE)
  }
  
  ##-- Species
  if (length(species)>1) {
    stop('This function can only deal with one species at a time... \nPlease, use one of "bear", "wolf", or "wolverine" for the n target species.')
  }
  if (sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0) {
    SPECIES <- "Brown bear"
    engSpecies <- "bear"
    norSpecies <- c("Bjørn", "BjÃ¸rn")
  } else {
    if (sum(grep("wolf", species, ignore.case = T))>0|
       sum(grep("ulv", species, ignore.case = T))>0) {
      SPECIES <- "Gray wolf"
      engSpecies <- "wolf"
      norSpecies <- "Ulv"
    } else {
      if (sum(grep("wolverine", species, ignore.case = T))>0|
          sum(grep("järv", species, ignore.case = T))>0|
         sum(grep("jerv", species, ignore.case = T))>0) {
        SPECIES <- "wolverine"
        engSpecies <- "wolverine"
        norSpecies <- "Jerv"
      } else {
        SPECIES <- engSpecies <- norSpecies <- species 
      }
    }
  }
  
  ##-- Years
  if (is.null(years)) { years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }
  
  ##-- Sampling months
  if (is.null(sampling.months)) {
    if (engSpecies == "bear") {
      sampling.months <- list(4:11)
    } else {
      if (engSpecies == "wolf") {
        sampling.months <- list(c(10:12),c(1:3))
      } else {
        if (engSpecies == "wolverine") {
          sampling.months <- list(c(10:12),c(1:4))
        } else {
          stop("No default setting available for the monitoring period of this species. \n You must specify the monitoring season months through the 'sampling.months' argument.")
        }
      }
    }
  }
  
  ##-- Renaming list
  if (is.null(rename.list)) {
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
  
  ##-- Load pre-processed habitat shapefiles
  data(COUNTRIES, envir = environment()) 
  
  ##-- data info
  DATE <- getMostRecent(path = data.dir, pattern = "DNA")
  
  ##-- Set file name for clean data
  fileName <- paste0("CleanData_", engSpecies, "_", DATE, ".RData")
  
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
  
  ##-----   2.1. RAW NGS DATA -----
  
  ##-- NGS data
  DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                          extension = ".xls",
                                          pattern = "DNA")) %>%
    ##-- Rename columns to facilitate manipulation
    dplyr::rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    dplyr::filter(., Species %in% norSpecies) %>%
    ##-- Remove any duplicates
    dplyr::distinct(., .keep_all = TRUE) %>%
    ##-- Add some columns
    dplyr::mutate( 
      ##-- Add "Country" column
      Country_sample = substrRight(County, 3),
      ##-- Change date format
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      ##-- Extract year
      Year = as.numeric(format(Date,"%Y")),
      ##-- Extract month
      Month = as.numeric(format(Date,"%m")),
      ##-- Extract sampling season
      ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
      ##-- Set all months in given sampling period to the same year)
      Season = ifelse( length(sampling.months) > 1, 
                       ifelse( Month < unlist(sampling.months)[1],
                               Year,
                               Year-1), 
                       Year),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", "unknown", Sex),
      Sex = ifelse(is.na(Sex), "unknown", Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    dplyr::filter(., Year %in% years) 
  
  
  ##-- Number of NGS samples
  NGS_samples <- table(DNA$Sex, DNA$Year, useNA = "ifany")
  NGS_samples <- rbind(NGS_samples, "Total" = colSums(NGS_samples))
  NGS_samples <- cbind(NGS_samples, "Total" = rowSums(NGS_samples))
  write.csv( NGS_samples,
             file = file.path( working.dir, "tables",
                               paste0(engSpecies, "_Raw NGS Samples_",
                                      years[1]," to ", years[length(years)],
                                      ".csv")))
  
  
  ##-- Number of individuals detected alive
  NGS_ids <- apply(table(DNA$Sex, DNA$Year, DNA$Id, useNA = "ifany"), c(1,2), function(x)sum(x>0))
  NGS_ids <- rbind(NGS_ids, "Total" = apply(table(DNA$Year,DNA$Id, useNA = "ifany"), 1, function(x)sum(x>0)))
  NGS_ids <- cbind(NGS_ids, "Total" = c(apply(table(DNA$Sex,DNA$Id, useNA = "ifany"), 1, function(x)sum(x>0)),length(unique(DNA$Id))))
  write.csv( NGS_ids,
             file = file.path( working.dir, "tables",
                               paste0(engSpecies, "_Raw NGS Ids_",
                                      years[1]," to ", years[length(years)],
                                      ".csv")))
  
  
  
  ##-----   2.2. RAW DEAD RECOVERY DATA -----
  
  ##-- Load raw excel file imported from rovbase 
  DR <- suppressWarnings(readMostRecent( path = data.dir,
                                         extension = ".xls",
                                         pattern = "dead")) %>%
    ##-- Rename columns to facilitate manipulation
    dplyr::rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    dplyr::filter(., Species %in% norSpecies) %>%
    ##-- Remove any duplicates
    dplyr::distinct(., .keep_all = TRUE) %>%
    ##-- Add some columns
    dplyr::mutate( 
      ##-- Add "Country" column
      Country_sample = substrRight(County, 3),
      ##-- Change date format
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      ##-- Extract year
      Year = as.numeric(format(Date,"%Y")),
      ##-- Extract month
      Month = as.numeric(format(Date,"%m")),
      ##-- Extract sampling season
      ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
      ##-- Set all months in given sampling period to the same year)
      Season = ifelse( length(sampling.months) > 1, 
                       ifelse( Month < unlist(sampling.months)[1],
                               Year,
                               Year-1), 
                       Year),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", "unknown", Sex),
      Sex = ifelse(is.na(Sex), "unknown", Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    dplyr::filter(., Year %in% years) 
  
  ##-- Number of DR samples
  DR_samples <- table(DR$Sex, DR$Year, useNA = "ifany")
  DR_samples <- rbind(DR_samples, "Total" = colSums(DR_samples))
  DR_samples <- cbind(DR_samples, "Total" = rowSums(DR_samples))
  write.csv( DR_samples,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Raw DR Samples_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  ##-- Number of individuals recovered dead
  DR_ids <- apply(table(DR$Sex, DR$Year, DR$Id, useNA = "ifany"), c(1,2), function(x)sum(x>0))
  DR_ids <- rbind(DR_ids, "Total" = apply(table(DR$Year,DR$Id, useNA = "ifany"), 1, function(x)sum(x>0)))
  DR_ids <- cbind(DR_ids, "Total" = c(apply(table(DR$Sex,DR$Id, useNA = "ifany"), 1, function(x)sum(x>0)),length(unique(DR$Id))))
  write.csv( DR_ids,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Raw DR Ids_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  
  
  ##-----   2.3. CHECKS -----
  
  ##-- Make sure all dead recoveries in DNA are in DR
  check1 <- all(DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$DNAID) 
  probs_DR_in_DNA <- NULL
  if(!check1){
    tmp <- DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"]
    probs_DR_in_DNA <- tmp[!tmp %in% DR$DNAID]
  } 
  
  tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
  test <- dplyr::anti_join(tmp,DR, by = names(tmp)[names(tmp) %in% names(DR)])
  
  ##-- Make sure that only "Dead recoveries" are in DR 
  check2 <- all(substr(DR$RovbaseID,1,1) %in% "M")
  probs_DNA_in_DR <- NULL
  if(!check2){
    probs_DNA_in_DR <- DR$DNAID[!substr(DR$RovbaseID,1,1) %in% "M"]
  }
  
  
  
  ##-----   2.4. MERGE -----
  
  ##-- Merge DNA and dead recoveries files using all shared names columns
  #DATA <- full_join(DNA, DR, by = names(DNA)[names(DNA) %in% names(DR)]) 
  DATA <- merge( DR, DNA,
                 by = c("Id","RovbaseID","DNAID","Species","Sex",
                        "Date","Year","Month","Season",
                        "East_UTM33","North_UTM33",
                        "County","Country_sample"),
                 all = TRUE)
  
  ##-- Determine Death and Birth Years
  DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
  DATA$RovbaseID <- as.character(DATA$RovbaseID)
  DATA$Death <- NA
  DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
  DATA$Birth <- DATA$Death - DATA$Age
  
  ##-- Extract useful numbers
  numNoID <- sum(is.na(DATA$Id))              ## number of samples without ID
  numNoDate <- sum(is.na(DATA$Year))          ## number of samples without Date
  numNoCoords <- sum(is.na(DATA$East_UTM33))  ## number of samples without Coords  
  
  ##-- Filter out unusable samples
  DATA <- DATA %>%
    dplyr::filter(., 
                  ##-- Filter out samples with no ID
                  !is.na(Id),
                  ##-- Filter out samples with no Coordinates
                  !is.na(East_UTM33),
                  ##-- Filter out samples with no dates  
                  !is.na(Year)) %>%
    droplevels(.)
  
  
  
  ##-----   2.5. SEX ASSIGNMENT -----
  
  ID <- unique(as.character(DATA$Id))
  DATA$Sex <- as.character(DATA$Sex)
  doubleSexID <- IdDoubleSex <- NULL  
   
  counter <- 1
  for(i in 1:length(ID)){
    ##-- Subset data to individual i
    tmp <- DATA$Sex[DATA$Id == ID[i]]
    
    ##-- Number of times individual i was assigned to each sex
    tab <- table(tmp[tmp %in% c("female","male")])
    
    ##-- If conflicting sexes (ID identified as both "female" and "male")
    if(length(tab) == 2){
      ##-- If ID assigned the same number of times to the 2 sexes, assign to unknown
      if(tab[1] == tab[2]){
        DATA$Sex[DATA$Id == ID[i]] <- "unknown"
      } else {
        ##-- Otherwise pick the most common sex
        DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
      }
      # print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))])) 
      IdDoubleSex[counter] <- ID[i]
      counter <- counter + 1
    }
    # ##-- If only one of "female" or "male" registered
    # if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
    # ##-- If anything else registered : "unknown"
    # if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "unknown"}
    doubleSexID[i] <- length(tab)
  }#i
  
  
  
  ##-----   2.6. SPLIT DATA -----
  
  ##-- Split DATA into alive and dead.recovery datasets
  alive <- DATA[is.na(DATA$Death), ]
  dead.recovery <- DATA[!is.na(DATA$Death), ]
  
  ##-- Add earlier detection index
  alive$detected.earlier <-
    unlist(lapply(1:nrow(alive),
                  function(i){
                    this.id <- alive[i,"Id"]
                    this.date <- alive[i,"Date"]
                    any(alive$Id %in% this.id & alive$Date < this.date)
                  }))
  
  dead.recovery$detected.earlier <-
    unlist(lapply(1:nrow(dead.recovery),
                  function(i){
                    this.id <- dead.recovery[i,"Id"]
                    this.date <- dead.recovery[i,"Date"]
                    any(alive$Id %in% this.id & alive$Date < this.date)
                  }))
  
  
  
  ##----- 3. SPECIES-SPECIFIC CLEANING STEPS ------
  
  ##-----   3.1. WOLVERINE -----
  
  if(engSpecies == "wolverine"){
    ##-- Remove un-verified dead recoveries [HB] 
    ##-- ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
    dead.recovery <- dead.recovery[!grepl(pattern = "Påskutt",
                                          x = as.character(dead.recovery$Outcome)), ]
    
    ##-- Remove suspect NGS samples according to Henrik
    SUSPECT_NGS_SAMPLES <- readMostRecent(
      path = data.dir,
      extension = ".xls",
      pattern = "Remove ngs samples list wolverine")
    alive$DNAID <- as.character(alive$DNAID)
    alive <- alive[!(alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
    
    ##-- Remove suspect dead recoveries according to Henrik
    SUSPECT_DeadRecoSAMPLES <- readMostRecent(
      path = data.dir,
      extension = ".xls",
      pattern = "Remove dead recoveries list wolverine")
    dead.recovery$DNAID <- as.character(dead.recovery$DNAID)
    dead.recovery <- dead.recovery[!(dead.recovery$RovbaseID %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
    
    
    ##-- Remove pups killed before recruitment based on weight (cf. Henrik)
    ##-- 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
    youngDeads <- which(dead.recovery$Age_class %in% "Unge" &
                              dead.recovery$Month > 2 &
                              dead.recovery$Month < 12)
    if(length(youngDeads) > 0){
      dead.recovery <- dead.recovery[-youngDeads, ]
    }
    
    
    ##-- 2) remove individuals with 0 <= weight < 4kg between March and November 
    ##-- Format the weight correctly 
    dead.recovery$Weight_total <- as.character(dead.recovery$Weight_total)
    dead.recovery$Weight_slaughter <- as.character(dead.recovery$Weight_slaughter)
    ##-- Convert to decimals
    dead.recovery$Weight_total <- as.numeric(gsub(",", ".", dead.recovery$Weight_total))
    dead.recovery$Weight_slaughter <- as.numeric(gsub(",", ".", dead.recovery$Weight_slaughter))
    ##-- Get the two weight columns together.
    dead.recovery$weight <- ifelse(!is.na(dead.recovery$Weight_total),
                                   dead.recovery$Weight_total,
                                   dead.recovery$Weight_slaughter)
    ##-- Assign negative values to nas to avoid issues
    dead.recovery$weight[is.na(dead.recovery$weight)] <- -999
    
    
    ##-- Check with Henrik (this step does not remove dead recoveries on id with weight==0 should it?)
    ##-- Check how many dead reco we remove and remove if more than 0
    lowWeightDeads <- which(dead.recovery$weight > 0 & dead.recovery$weight < 4 &
                                  dead.recovery$Month > 2 & dead.recovery$Month < 12)
    if(length(lowWeightDeads) > 0){
      dead.recovery <- dead.recovery[-lowWeightDeads, ]
    }
    
    ##-- Check how many dead reco with a weight of 0 kg and recovered between March and November
    zeroWeightDeads <- which(dead.recovery$Age %in% 0 &
                                   dead.recovery$Month > 2 &
                                   dead.recovery$Month < 12)
  }
  
  
  
  ##-----   3.2. WOLF -----
  
  if(engSpecies == "wolf"){
    ##-- Load most recent Micke's file
    INDIVIDUAL_ID <- readMostRecent.csv(
      path = data.dir,
      pattern = "_ID Grouping ",
      fileEncoding = "latin1")  
    
    ##-- Translate Scandinavian characters
    colnames(INDIVIDUAL_ID) <- translateForeignCharacters( data = colnames(INDIVIDUAL_ID))
    
    ##-- Overwrite gender from Micke's data when available
    micke.sex <- as.character(unlist(lapply(DATA$Id,
                                            function(i){ 
                                              INDIVIDUAL_ID[as.character(INDIVIDUAL_ID$Individ..Rovbase.) %in% i,"Sex"][1]
                                            })))
    micke.sex[micke.sex %in% "0"] <- NA
    micke.sex[micke.sex %in% names(table(micke.sex))[3]] <- NA
    micke.sex[micke.sex %in% "Hona"] <- "female"
    micke.sex[micke.sex %in% "Hane"] <- "male"
    new.sex <- ifelse(!is.na(micke.sex), as.character(micke.sex), as.character(DATA$Sex))
    DATA$Sex <- new.sex
    
    numOverwiteSex <- sum(unique(as.character(INDIVIDUAL_ID$Individ..Rovbase.)) %in% DATA$Id)
  }
  
  
  
  ##-----   3.3. BEAR -----
  
  if(engSpecies == "bear"){
    ##-- Load most recent "flagged" file from HB
    flagged <- readMostRecent( 
      path = data.dir,
      extension = ".csv",
      pattern = "dna_bear_to_remove", 
      fileEncoding = "Latin1") 
    
    ##-- Remove flagged samples 
    remove.alive <- !alive$Barcode_sample %in% flagged$Strekkode
    alive <- alive[remove.alive, ]
    remove.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
    dead.recovery <- dead.recovery[remove.dead, ]
    dead.recovery$Missing <- NA
    dead.recovery$Individ <- NA
  }
  
  
  
  ##----- 4. TURN INTO .sf OBJECTS -----
  
  ##-- Turn into sf points dataframe
  alive <- sf::st_as_sf( x = alive,
                         coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(., sf::st_crs(32633)) 
  
  ##-- Intersect and extract country name
  #alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
  alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  ##-- Turn into sf points dataframe
  dead.recovery <- sf::st_as_sf( x = dead.recovery,
                                 coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633))
  
  ##-- Intersect and extract country name
  #dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]
  dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  
  ##----- 5. DATA SUMMARY -----
  
  ##-----   5.1. DATA SUMMARY - TABLES -----
  
  ##-- Number of NGS samples per year and country (date,rovbase)
  samples <- table(alive$Country_sample, alive$Year)
  samples <- rbind(samples, "Total" = colSums(samples))
  samples <- cbind(samples, "Total" = rowSums(samples))
  write.csv( samples,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Clean NGS Samples_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  ##-- Number of individuals detected alive
  ids <- apply(table(alive$Country_sample, alive$Year, alive$Id), c(1,2), function(x)sum(x>0))
  ids <- rbind(ids, "Total" = apply(table(alive$Year, alive$Id), 1, function(x)sum(x>0)))
  ids <- cbind(ids, "Total" = c(apply(table(alive$Country_sample,alive$Id), 1, function(x)sum(x>0)), length(unique(alive$Id))))
  write.csv( ids,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Clean NGS Ids_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  ##-- Number of DR samples
  deadSamples <- table(dead.recovery$Country_sample, dead.recovery$Year)
  deadSamples <- rbind(deadSamples, "Total" = colSums(deadSamples))
  deadSamples <- cbind(deadSamples, "Total" = rowSums(deadSamples))
  write.csv( deadSamples,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Clean DR Samples_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  ##-- Number of individuals recovered
  deadIds <- apply(table(dead.recovery$Country_sample,dead.recovery$Year,dead.recovery$Id),c(1,2),function(x)sum(x>0))
  deadIds <- rbind(deadIds, "Total" = apply(table(dead.recovery$Year,dead.recovery$Id), 1, function(x)sum(x>0)))
  deadIds <- cbind(deadIds, "Total" = c(apply(table(dead.recovery$Country_sample,dead.recovery$Id), 1, function(x)sum(x>0)), length(unique(dead.recovery$Id))))
  write.csv( deadIds,
             file = file.path( working.dir, "tables",
                               paste0( engSpecies, "_Clean DR Ids_",
                                       years[1]," to ", years[length(years)],
                                       ".csv")))
  
  
  
  ##-----   5.2. NUMBER OF SAMPLES - FIGURE -----
  
  ##-- Number of NGS per month
  dat.alive <- alive %>%
    dplyr::mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = dplyr::n())
  dat.alive$type = "NGS"
  
  ##-- Number of dead recoveries per month
  dat.dead <- dead.recovery %>%
    dplyr::mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = -dplyr::n())
  dat.dead$type = "dead.recovery"
  
  ##-- Combine NGS and dead recoveries
  dat <- rbind(dat.alive, dat.dead)
  dat$Date <- as.Date(dat$Date)
  
  ##-- Plot time series of number of samples per month
  NGS_plot <- ggplot(dat) +
    geom_col(aes(x = Date, y = n, fill = type)) +
    ylab("Number of samples") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position.inside = c(0.1,0.9),
          axis.text.x = element_text(angle = 60,
                                     hjust = 1)) +
    scale_x_date( date_breaks = "years",
                  date_labels = "%Y") 

  ##-- Export as .png
  ggsave(filename = file.path(working.dir, "figures",
                              paste0( engSpecies, "_Clean Rovbase Samples_",
                                      years[1]," to ", years[length(years)],
                                      ".png")),
         plot = NGS_plot,
         dpi = 300, height = 6, width = 12, device = "png")
  # grDevices::png( filename = file.path(working.dir, "figures",
  #                                      paste0( engSpecies, "_Clean Rovbase Samples_",
  #                                              years[1]," to ", years[length(years)],
  #                                              ".png")),
  #                 width = 8, height = 6,
  #                 units = "in", pointsize = 12,
  #                 res = 300, bg = NA)  
  # NGS_plot
  # dev.off()

  
  
  ##-----   5.3. NUMBER OF INDIVIDUALS - FIGURE -----
  
  ##-- Number of IDs
  dat.alive <- alive %>% 
    dplyr::group_by(Year) %>% 
    dplyr::summarise(n = length(unique(Id)))
  dat.alive$type = "NGS"
  
  ##-- Number of Dead Recoveries
  dat.dead <- dead.recovery %>% 
    dplyr::group_by(Year) %>% 
    dplyr::summarise(n = -length(unique(Id)))
  dat.dead$type = "dead.recovery"
  
  ##-- Combine NGS and dead recoveries
  dat <- rbind(dat.alive, dat.dead)
  
  ##-- Plot time series of number of IDs per year
  IDS_plot <- ggplot2::ggplot(dat) +
    geom_col(aes(x = Year, y = n, fill = type)) +
    ylab("Number of individuals") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position.inside = c(0.1,0.9),
          axis.text.x = element_text(angle = 60,
                                     hjust = 1)) +
    scale_x_continuous( breaks = years,
                        labels = years)
  ##-- Export as .png
  ggsave(filename = file.path(working.dir, "figures",
                              paste0( engSpecies, "_Clean Rovbase Ids_",
                                      years[1]," to ", years[length(years)],
                                      ".png")),
         plot = IDS_plot,
         dpi = 300, height = 6, width = 12, device = "png")
  
  # grDevices::png( filename = file.path(working.dir, "figures",
  #                                      paste0(engSpecies, "_Clean Rovbase Ids_",
  #                                             years[1]," to ", years[length(years)],
  #                                             ".png")),
  #                 width = 8, height = 6,
  #                 units = "in", pointsize = 12,
  #                 res = 300, bg = NA)
  # IDS_plot
  # dev.off()
  
  
  
  ##-----   5.4. SAMPLING MAPS - FIGURE ------

  ##-- Maps layout
  L <- length(years)
  if(L < 6){ nrows <- 1 } else{
    if(L < 13){ nrows <- 2 } else {
      if(L < 22){ nrows <- 3 } else {
        if(L < 33){ nrows <- 4 } else {
          nrows <- 5
        }}}}
  ncols <- ceiling(L/nrows)
  
  ##-- Save maps as .png
  grDevices::png(filename = file.path( working.dir, "figures",
                                       paste0( engSpecies, "_Clean Rovbase Maps_",
                                               years[1]," to ", years[length(years)],
                                               ".png")),
                 width = ncols*2, height = nrows*4,
                 units = "in", pointsize = 12,
                 res = 300, bg = NA)

  ##-- Maps layout
  mx <- matrix(NA, nrow = nrows*2, ncol =  (ncols*2)+1)
  for(r in 1:nrows){
    mx[r*2-1, ] <- c(1,rep(1:ncols, each = 2)) + (r-1)*ncols
    mx[r*2, ] <- c(rep(1:ncols, each = 2),ncols) + (r-1)*ncols
  }#r
  nf <- graphics::layout(mx,
                         widths = c(rep(1,ncol(mx))),
                         heights = rep(1,2))
  par(mar = c(0,0,0,0))
  
  for(t in 1:length(years)){
    ##-- Plot maps
    plot( sf::st_geometry(COUNTRIES), border = NA, col = c("gray80","gray60"))
    try(
      plot( sf::st_geometry(alive[alive$Year == years[t], ]), add = TRUE, col = "orange", pch = 3),
      silent = TRUE)
    try(
      plot( sf::st_geometry(dead.recovery[dead.recovery$Year == years[t], ]), add = TRUE, col = "slateblue", pch = 3),
      silent = TRUE)
    plot( sf::st_geometry(COUNTRIES), border = "gray40", col = NA, add = TRUE)
    
    ##-- Add year
    graphics::mtext(text = years[t],
          side = 1, line = -18,
          adj = 0.18, cex = 1.2)
  }#t
  dev.off()
  
  
  
  ##-----   5.5. PREVIOUSLY DETECTED - FIGURE -----
  
  ##-- Plot number of individuals with previous NGS detections
  plot1 <- alive %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id))) %>%
    ggplot2::ggplot() +
    geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
    labs( tag = "A",
          x = "years",
          y = "Number of individuals detected through NGS") +
    theme( legend.position = "none",
           legend.title = element_blank(),
           axis.text.x = element_text(angle = 60,
                                      hjust = 1)) +
    scale_x_continuous(breaks = years, labels = years) +
    scale_fill_manual(values = c("gray20", "gray60"))
  
  ##-- Plot number of dead recoveries with previous detections
  plot2 <- dead.recovery %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id))) %>%
    ggplot2::ggplot() +
    ggplot2::geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
    labs( tag = "B",
          x = "years",
          y = "Number of individuals recovered dead") +
    theme( legend.title = element_blank(),
           legend.position.inside = 2,
           axis.text.x = element_text(angle = 60,
                                      hjust = 1)) +
    scale_x_continuous(breaks = years, labels = years) +
    scale_fill_manual(values = c("gray20", "gray60"))
  plot_total <- plot1 + plot2
  
  #-- Save figure
  ggsave(filename = file.path(working.dir, "figures",
                              paste0( engSpecies, "_Previous Detection_",
                                      years[1]," to ", years[length(years)],
                                      ".png")),
         plot = plot_total,
         dpi = 300, height = 6, width = 12, device = "png")
  # grDevices::png( filename = file.path( working.dir, "figures",
  #                                       paste0( engSpecies, "_Previous Detection_",
  #                                               years[1]," to ", years[length(years)],
  #                                               ".png")),
  #                 width = 12, height = 6,
  #                 units = "in", pointsize = 12,
  #                 res = 300, bg = NA)
  # dev.off()
  
  
  
  ##----- 6. DATA ISSUES -----
  
  ##-----   6.1. MULTIPLE DEATHS ------
  
  ##-- Identify and count individuals dead "more than once"
  ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
  multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
  for (i in 1:length(ID)) {
    tmp <- dead.recovery[dead.recovery$Id == ID[i], ] 
    ##-- Multiple death dates
    if(length(unique(tmp$Date)) > 1){
      multiDeathDate <- c(multiDeathDate, ID[i])
    }
    ##-- Multiple death years
    if(length(unique(tmp$Year)) > 1){
      multiDeathYear <- c(multiDeathYear, ID[i])
    }
    ##-- Multiple death locations
    if(length(unique(tmp$East)) > 1 | length(unique(tmp$North)) > 1){
      multiDeathLocs <- c(multiDeathLocs, ID[i])
    }
  }#i
  
  ##-- Remove individuals that died more than once
  dead.recovery$Id <- as.character(dead.recovery$Id)
  IdDoubleDead <- names(table(dead.recovery$Id))[table(dead.recovery$Id) > 1]
  if (length(IdDoubleDead) > 0) {
    for(i in IdDoubleDead){
      ##-- Identify repeated deaths
      tmp <- which(dead.recovery$Id %in% i) 
      ##-- Try to keep death with known death cause 
      tmp2 <- which(!is.na(dead.recovery$DeathCause_2[tmp]))[1]
      if(length(tmp2) == 0){tmp <- tmp[-1]} else {tmp <- tmp[!tmp %in% tmp2]}
      ##-- Remove repeated deaths
      dead.recovery <- dead.recovery[-tmp, ]
    }#i
  }#if
  
  
  
  ##-----   6.2. GHOST INDIVIDUALS ------
  
  id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
  ghosts <- unlist(lapply(id.list, function(id) {
    out <- NULL
    try({
      if(id %in% dead.recovery$Id){
        mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
        this.alive <- alive[alive$Id == id, ]
        ##-- Was it detected alive in any season after death?
        temp <- this.alive[this.alive$Year > mort.year, ]
        if(length(temp) > 0){
          out <- rownames(temp)
          names(out) <- id
        }
      }
    }, silent = TRUE)
    return(out)
  }))
  samples.to.remove <- unlist(ghosts)
  
  ##-- Remove flagged NGS detections after dead recovery
  alive <- alive[!rownames(alive) %in% samples.to.remove, ]
  
  
  
  ##----- 7. SAVE DATA ------
  
  save( alive, 
        dead.recovery,
        file = file.path( working.dir, "data", fileName))
                         
  
  
  ##----- 8. PRINT REPORT -----
  
  if (print.report) {
    
    ##-- List of rendering parameters for the .Rmd report
    info.ls <- list()
    info.ls$IdDoubleSex <- IdDoubleSex
    info.ls$doubleSexID <- doubleSexID
    info.ls$numNoID <- numNoID
    info.ls$numNoDate <- numNoDate
    info.ls$numNoCoords <- numNoCoords
    info.ls$multiDeathDate <- multiDeathDate
    info.ls$multiDeathDate <- multiDeathDate
    info.ls$multiDeathYear <- multiDeathYear
    info.ls$multiDeathLocs <- multiDeathLocs
    info.ls$samples.to.remove <- samples.to.remove
    if (engSpecies == "bear") {
      info.ls$remove.alive <- remove.alive
      info.ls$remove.dead <- remove.dead
    }
    if (engSpecies == "wolf") {}
    if (engSpecies == "wolverine") {}
    
    
    ##-- Find the .rmd template for the report
    if(is.null(Rmd.template)) {
      Rmd.template <- system.file("rmd", "RovQuant_CleaningReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop("Can not find the Rmarkdown document to use for cleaning Rovbase.3.0 data.\n You must provide the path to the Rmarkdown template through the \"Rmd_template\" argument.")
      }
    }
    
    ##-- Find the directory to print the report
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    ##-- Render .Rmd report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = SPECIES, 
                     years = years,
                     sampling.months = sampling.months,
                     data.dir = data.dir,
                     working.dir = working.dir,
                     date = DATE,
                     info.ls = info.ls),
      output_dir = output.dir,
      output_file = paste0("CleanData_", engSpecies, "_", DATE,".html"))
  }
}
