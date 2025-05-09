#' @title Brown bear Data set clean-up.
#'
#' @description
#' \code{cleanRovbaseData_bear} identifies and loads the most recent Brown bear data
#'  extracted from Rovbase.no in the specified \code{data.dir} location, and 
#'  conducts a set of data cleaning steps.  @seealso [cleanRovbaseData()] for more details on the data cleaning process.
#'
#' @name cleanRovbaseData_bear
#' 
#' @param data.dir the \code{path} pointing to the directory containing the raw 
#' data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. By default,
#'  the cleaned data will be stored in a subfolder of this working directory called 'data'.
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param sex A \code{character} vector containing the sex of interest. 
#' Can be "Hunn" for females or "Hann" for males. Default is both sexes (c("Hunn","Hann")).
#' @param sampling.months (Optional) A \code{list} containing the sampling period months.
#' If the sampling period overlaps two calendar years, the list should contain one element per year.
#' (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param rename.list (Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.


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
#' @importFrom sf st_intersects
#' @import ggplot2 
#' @import dplyr 


NULL
#' @rdname cleanRovbaseData_bear
#' @export
cleanRovbaseData_bear <- function(
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
  
  ##-- Species
  if (length(species)>1) {
    stop('This function can only deal with one species at a time... \nPlease, use one of "bear", "wolf", or "wolverine" for the n target species.')
  }
  if(sum(grep("bear", species, ignore.case = T))>0|
     sum(grep("bjørn", species, ignore.case = T))>0|
     sum(grep("bjorn", species, ignore.case = T))>0){
    species <- "bear"
    engSpecies <- "Brown bear"
    norSpecies <- c("Bjørn", "BjÃ¸rn")
  } else {
    if(sum(grep("wolf", species, ignore.case = T))>0|
       sum(grep("ulv", species, ignore.case = T))>0){
      species <- "wolf"
      engSpecies <- "Gray wolf"
      norSpecies <- "Ulv"
    } else {
      if(sum(grep("wolverine", species, ignore.case = T))>0|
         sum(grep("jerv", species, ignore.case = T))>0){
        species <- "wolverine"
        engSpecies <- "wolverine"
        norSpecies <- "Jerv"
      } else {
        engSpecies <- norSpecies <- species 
      }
    }
  }
  
  ##-- Years
  if (is.null(years)) { years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }
  
  ##-- Sampling months
  if (is.null(sampling.months)) {
    if (species == "bear") {
      sampling.months <- list(4:11)
    } else {
      if (species == "wolf") {
        sampling.months <- list(c(10:12),c(1:3))
      } else {
        if (species == "wolverine") {
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
  
  ##-----   2.1. RAW NGS DATA -----
  
  ##-- NGS data
  DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                          extension = ".xls",
                                          pattern = "DNA")) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    filter(., Species %in% norSpecies) %>%
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
      Season = ifelse( Month < unlist(sampling.months)[1],
                       Year,
                       Year-1),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", NA, Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    filter(., Year %in% years) 
  
  
  ##-- Number of NGS samples
  NGS_samples <- table(DNA$Sex, DNA$Year, useNA = "ifany")
  NGS_samples <- rbind(NGS_samples, "Total" = colSums(NGS_samples))
  NGS_samples <- cbind(NGS_samples, "Total" = rowSums(NGS_samples))
  write.csv(NGS_samples,
            file = file.path(working.dir, "tables",
                             paste0(species, "_NGS samples raw_",
                                    years[1]," to ", years[length(years)], ".csv")))
  

  ##-- Number of individuals detected alive
  NGS_ids <- apply(table(DNA$Sex, DNA$Year, DNA$Id, useNA = "ifany"),
                   c(1,2), function(x)sum(x>0))
  NGS_ids <- rbind(NGS_ids,
                   "Total" = apply(table(DNA$Year,DNA$Id, useNA = "ifany"),
                                   1, function(x)sum(x>0)))
  NGS_ids <- cbind(NGS_ids,
                   "Total" = c(apply(table(DNA$Sex,DNA$Id, useNA = "ifany"),
                                     1, function(x)sum(x>0)),length(unique(DNA$Id))))
  write.csv(NGS_ids,
            file = file.path(working.dir, "tables",
                             paste0(species, "_NGS ids raw_",
                                    years[1]," to ", years[length(years)], ".csv")))

  
  
  ##-----   2.2. RAW DEAD RECOVERY DATA -----
  
  ##-- Load raw excel file imported from rovbase 
  DR <- suppressWarnings(readMostRecent( path = data.dir,
                                         extension = ".xls",
                                         pattern = "dead")) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    filter(., Species %in% norSpecies) %>%
    ##-- Remove any duplicates
    distinct(., .keep_all = TRUE) %>%
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
      Season = ifelse( Month < unlist(sampling.months)[1],
                       Year,
                       Year-1),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", NA, Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    filter(., Year %in% years) 
  
  ##-- Number of DR samples
  DR_samples <- table(DR$Sex, DR$Year, useNA = "ifany")
  DR_samples <- rbind(DR_samples, "Total" = colSums(DR_samples))
  DR_samples <- cbind(DR_samples, "Total" = rowSums(DR_samples))
  write.csv(DR_samples,
            file = file.path( working.dir, "tables",
                              paste0( species, "_DR samples raw_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  ##-- Number of individuals detected alive
  DR_ids <- apply(table(DR$Sex, DR$Year, DR$Id, useNA = "ifany"),
                  c(1,2),
                  function(x)sum(x>0))
  DR_ids <- rbind(DR_ids,
                  "Total" = apply(table(DR$Year,DR$Id, useNA = "ifany"),
                                  1,
                                  function(x)sum(x>0)))
  DR_ids <- cbind(DR_ids,
                  "Total" = c(apply(table(DR$Sex,DR$Id, useNA = "ifany"),
                                    1,
                                    function(x)sum(x>0)),length(unique(DR$Id))))
  write.csv(DR_ids,
            file = file.path( working.dir, "tables",
                              paste0( species, "_DR ids raw_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  
  
  ##-----   2.3. CHECKS -----
  
  ##-- Make sure all dead recoveries in DNA are in DR
  check1 <- all(DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$DNAID) 
  probs_DR_in_DNA <- NULL
  if(!check1){
    tmp <- DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"]
    probs_DR_in_DNA <- tmp[!tmp %in% DR$DNAID]
  } 
  
  tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
  test <- anti_join(tmp,DR, by = names(tmp)[names(tmp) %in% names(DR)])
  
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
                  !is.na(Year),
                  ##-- Filter out samples with 
                  Year %in% years) %>%
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
    
    ##-- If only one of "female" or "male" registered
    if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
    
    ##-- If anything else registered : "Ukjent"
    if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "unknown"}
    
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
  
  
  
  ##-----   2.7. WOLVERINE -----
  
  if(engSpecies == "wolverine"){
    ##-- Remove un-verified dead recoveries [HB] 
    ##-- ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
    dead.recovery <- dead.recovery[!grepl(pattern = "Påskutt",
                                          x = as.character(dead.recovery$Outcome)), ]
    
    ##-- Remove suspect NGS samples according to Henrik
    SUSPECT_NGS_SAMPLES <- readMostRecent(
      path = dir.in,
      extension = ".xls",
      pattern = "Remove ngs samples list wolverine")
    alive$DNAID <- as.character(alive$DNAID)
    alive <- alive[!(alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
    
    ##-- Remove suspect dead recoveries according to Henrik
    SUSPECT_DeadRecoSAMPLES <- readMostRecent(
      path = dir.in,
      extension = ".xls",
      pattern = "Remove dead recoveries list wolverine")
    dead.recovery$DNAID <- as.character(dead.recovery$DNAID)
    dead.recovery <- dead.recovery[!(dead.recovery$RovbaseID %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
    
    
    ##-- Remove pups killed before recruitment based on weight (cf. Henrik)
    ##-- 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
    out$youngDeads <- which(dead.recovery$Age_class %in% "Unge" &
                          dead.recovery$Month > 2 &
                          dead.recovery$Month < 12)
    if(length(youngDeads) > 0){
      dead.recovery <- dead.recovery[-out$youngDeads, ]
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
    out$lowWeightDeads <- which(dead.recovery$weight > 0 & dead.recovery$weight < 4 &
                              dead.recovery$Month > 2 & dead.recovery$Month < 12)
    if(length(lowWeightDeads) > 0){
      dead.recovery <- dead.recovery[-out$lowWeightDeads, ]
    }
    
    ##-- Check how many dead reco with a weight of 0 kg and recovered between March and November
    out$zeroWeightDeads <- which(dead.recovery$Age %in% 0 &
                               dead.recovery$Month > 2 &
                               dead.recovery$Month < 12)
  }
  
  

  ##-----   2.8. WOLF -----
  
  if(engSpecies == "wolf"){
    ##-- Load most recent Micke's file
    INDIVIDUAL_ID <- readMostRecent.csv(
      path = dir.in,
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
    
    out$numOverwiteSex <- sum(unique(as.character(INDIVIDUAL_ID$Individ..Rovbase.)) %in% DATA$Id)
  }
  
  
  
  ##-----   2.9. BEAR -----
  
  if(engSpecies == "bear"){
    ##-- Load most recent "flagged" file from HB
    flagged <- readMostRecent( 
      path = dir.in,
      extension = ".csv",
      pattern = "dna_bear_to_remove", 
      fileEncoding = "Latin1") 
    
    ##-- Remove flagged samples 
    out$remove.alive <- !alive$Barcode_sample %in% flagged$Strekkode
    alive <- alive[out$remove.alive, ]
    out$remove.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
    dead.recovery <- dead.recovery[out$remove.dead, ]
    dead.recovery$Missing <- NA
    dead.recovery$Individ <- NA
  }
  
  
  
  ##-----   2.10. TURN INTO .sf OBJECTS -----
  
  ##-- Turn into sf points dataframe
  alive <- sf::st_as_sf( x = alive,
                         coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633)) 
  
  ##-- Intersect and extract country name
  #alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
  alive$Country_sf[!is.na(as.numeric(st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  alive$Country_sf[!is.na(as.numeric(st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  ##-- Turn into sf points dataframe
  dead.recovery <- sf::st_as_sf( x = dead.recovery,
                                 coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633))
  
  ##-- Intersect and extract country name
  #dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]
  dead.recovery$Country_sf[!is.na(as.numeric(st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  dead.recovery$Country_sf[!is.na(as.numeric(st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  
  ##-----   2.11. DATA SUMMARY -----
  
  ##-----     2.11.1. DATA SUMMARY - TABLES -----
  
  ##-- Number of NGS samples
  samples <- table(alive$Country_sample, alive$Year)
  samples2 <- table(alive$Country_sf, alive$Year)
  samples <- rbind(samples, "Total" = colSums(samples))
  samples <- cbind(samples, "Total" = rowSums(samples))
  write.csv(samples,
            file = file.path( working.dir, "tables",
                              paste0( species, "_NGS samples clean_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  
  ##-- Number of individuals detected alive
  ids <- apply(table(alive$Country_sample,
                     alive$Year,
                     alive$Id),
               c(1,2),
               function(x)sum(x>0))
  
  ids <- rbind(ids,
               "Total" = apply(table(alive$Year,
                                     alive$Id),
                               1,
                               function(x)sum(x>0)))
  
  ids <- cbind(ids,
               "Total" = c(apply(table(alive$Country_sample,alive$Id),
                                 1,
                                 function(x)sum(x>0)),
                           length(unique(alive$Id))))
  write.csv(ids,
            file = file.path( working.dir, "tables",
                              paste0( species, "_NGS ids clean_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  
  ##-- Number of DR samples
  deadSamples <- table(dead.recovery$Country_sample, dead.recovery$Year)
  deadSamples <- rbind(deadSamples, "Total" = colSums(deadSamples))
  deadSamples <- cbind(deadSamples, "Total" = rowSums(deadSamples))
  write.csv(ids,
            file = file.path( working.dir, "tables",
                              paste0( species, "_DR samples clean_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  
  ##-- Number of individuals recovered
  deadIds <- apply(table(dead.recovery$Country_sample,
                         dead.recovery$Year,
                         dead.recovery$Id),
                   c(1,2),
                   function(x)sum(x>0))
  
  deadIds <- rbind(deadIds,
                   "Total" = apply(table(dead.recovery$Year,
                                         dead.recovery$Id),
                                   1,
                                   function(x)sum(x>0)))
  
  deadIds <- cbind(deadIds,
                   "Total" = c(apply(table(dead.recovery$Country_sample,
                                           dead.recovery$Id),
                                     1,
                                     function(x)sum(x>0)),
                               length(unique(dead.recovery$Id))))
  write.csv(ids,
            file = file.path( working.dir, "tables",
                              paste0( species, "_DR ids clean_",
                                      years[1]," to ", years[length(years)], ".csv")))
  
  
  
  ##-----     2.11.2. NUMBER OF SAMPLES - FIGURE -----
  
  ##-- Number of NGS per month
  dat.alive <- alive %>%
    mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = n())
  dat.alive$type = "NGS"
  
  ##-- Number of dead recoveries per month
  dat.dead <- dead.recovery %>%
    mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = -n())
  dat.dead$type = "dead.recovery"
  
  ##-- Combine NGS and dead recoveries
  dat <- rbind(dat.alive, dat.dead)
  dat$Date <- as.Date(dat$Date)
  
  ##-- Plot time series of number of samples per month
  grDevices::png( filename = file.path(working.dir, "figures",
                                       paste0( species, "_clean data samples_",
                                               years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  ggplot(dat) +
    geom_col(aes(x = Date, y = n, fill = type)) +
    ylab("Number of samples") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position.inside = c(0.1,0.9),
          axis.text.x=element_text(angle = 60,
                                   hjust = 1)) +
    scale_x_date( date_breaks = "years",
                  date_labels = "%Y") 
  graphics.off()
  

  
  ##-----     2.11.3. NUMBER OF INDIVIDUALS - FIGURE -----
  
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
  grDevices::png( filename = file.path(working.dir, "figures",
                                       paste0(species, "_clean data ids_",
                                              years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  ggplot(dat) +
    geom_col(aes(x = Year, y = n, fill = type)) +
    ylab("Number of individuals") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position.inside = c(0.1,0.9)) +
    scale_x_continuous( breaks = years,
                        labels = years)
  graphics.off()
  
  

  ##-----     2.11.4. PREVIOUSLY DETECTED - FIGURES -----

  ##-- Plot number of individuals with previous NGS detections
  grDevices::png( filename = file.path(working.dir, "figures",
                                       paste0(species, "_NGS previously detected_",
                                              years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  alive %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id)))%>%
    ggplot() +
    geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
    ylab("Number of individuals detected through NGS") +
    theme(legend.position = c(0.1,0.9)) +
    scale_x_continuous(breaks = years, labels = years) +
    scale_fill_manual(values = c("gray20", "gray60"))
  graphics.off()
  

  ##-- Plot number of dead recoveries with previous detections
  grDevices::png( filename = file.path(working.dir, "figures",
                                       paste0(species, "_DR previously detected_",
                                              years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  dead.recovery %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id))) %>%
    ggplot() +
    geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
    ylab("Number of individuals recovered dead") +
    theme(legend.position = c(0.1,0.9)) +
    scale_x_continuous(breaks = years, labels = years) +
    scale_fill_manual(values = c("gray20", "gray60"))
  graphics.off()
  
  
  
  ##-----   2.12. DATA ISSUES -----
  
  sexTab <- cbind.data.frame(
    "problems" = c("Unknown sex", "both 'female' and 'male'"),
    "number of individuals" = as.numeric(table(doubleSexID)[c(1,3)]))
  
  kable(sexTab, align = "lc") %>%
    kable_styling(full_width = F)
  
  
  
  ## ---- multiple deaths --------------------------------------------------------
  ##-- Identify and count individuals dead "more than once"
  ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
  multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
  for(i in 1:length(ID)){
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
  if(length(IdDoubleDead) > 0){
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
  
  
  
  ## ---- ghost individuals ------------------------------------------------------
  id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
  ghosts <- unlist(lapply(id.list, function(id){
    out <- NULL
    try({
      if(id %in% dead.recovery$Id){
        mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
        this.alive <- alive[alive$Id == id, ]
        ## Was it detected alive in any season after death?
        temp <- this.alive[this.alive$Year > mort.year, ]
        if(length(temp) > 0){
          out <- rownames(temp)
          names(out) <- id
        }## FLAG THOSE FOR HENDRIK
      }
    }, silent = TRUE)
    return(out)
  }))
  samples.to.remove <- unlist(ghosts)
  
  ##-- Remove flagged NGS detections after dead recovery
  alive <- alive[!rownames(alive) %in% samples.to.remove, ]
  
  
  
  ## ---- NGS maps ----------------------------------
  ##-- NGS map
  numRows <- ceiling(length(years)/5)
  numCols <- 5
  NGS_map <- ggplot(data = alive) +
    geom_sf(data = COUNTRIES,
            aes(fill = ISO),
            alpha = 0.3,
            color = NA) +
    geom_sf(color = "black",
            alpha = 0.3, size = 0.8, pch = 3) +
    facet_wrap(~Year, nrow = numRows, ncol = numCols) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  
  ##-- Save maps as .png
  grDevices::png(filename = file.path(working.dir, "figures",
                                      paste0(species, "_NGS_", years[1]," to ", years[length(years)], ".png")),
                 width = 8, height = 6, units = "in", res = 300)
  NGS_map
  graphics.off()
  
  ## ---- Dead recovery maps ------------------------
  dead_map <- ggplot(data = dead.recovery) +
    geom_sf(data = COUNTRIES, 
            aes(fill = ISO),
            alpha = 0.3,
            color = NA) + 
    geom_sf(color = "black", alpha = 0.5, size = 0.8, pch = 3) +
    facet_wrap(~Year, nrow = numRows, ncol = numCols) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  
  ##-- Save maps as .png
  grDevices::png(filename = file.path(working.dir, "figures",
                                      paste0(engSpecies, "_DEAD_", years[1]," to ", years[length(years)], ".png")),
                 width = 8, height = 6, units = "in", res = 300)
  dead_map
  graphics.off()
  
  
  ## ---- save data ---------------------------------------------------------------------------------------
  fileName <-  paste0("CleanData_", engSpecies, "_",DATE,".RData")
  
  save( alive, 
        dead.recovery,
        IdDoubleSex,
        file = file.path(working.dir, "data", fileName))
  
  
  
  ## ----- OUTPUT ------
  ##-- List of outputs for the .Rmd report
  out <- list()
  out$SPECIES <- "bear"
  out$YEARS <- years
  out$SEX <- sex
  out$DATE <- DATE
  out$info.ls <- list()
  out$info.ls$IdDoubleSex <- IdDoubleSex
  out$info.ls$doubleSexID <- doubleSexID
  out$info.ls$numNoID <- numNoID
  out$info.ls$numNoDate <- numNoDate
  out$info.ls$numNoCoords <- numNoCoords
  if(species == "bear"){}
  if(species == "wolf"){}
  if(species == "wolverine"){}
  
  
  return(out)
  
}