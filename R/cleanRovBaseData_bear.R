#' @title Brown bear dataset clean-up 
#'
#' @description
#' \code{cleanRovbaseData_bear} is an internal function called by the \code{cleanRovbaseData} function.
#' It identifies and loads the most recent Brown bear data downloaded from Rovbase
#' and conducts a set of cleaning steps that include:
#'  - removing un-identified samples
#'  - checking sex-assignment
#'  - removing problematic samples flagged by RovData.
#'
#' @name cleanRovbaseData_bear
#' 
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
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
#' @importFrom dplyr distinct rename filter
#' @importFrom sf st_as_sf st_crs st_set_crs st_intersects

NULL
#' @rdname cleanRovbaseData_bear
#' @export
cleanRovbaseData_bear <- function(
    years = NULL, 
    samplingMonths = NULL,
    data_dir = "./Data",
    output_dir = "./Data",
    Rmd_template = NULL,
    overwrite = FALSE){
  
  ##----- 0. INITIAL CHECKS -----
  ##-- Extract the date from the last .xlsx data file
  DATE <- getMostRecent( path = data_dir, pattern = "DNA.xlsx")
  
  ##-- Get file name for clean data
  fileName <- paste0("Data_Bear_", DATE,".RData")
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  if(file.exists(file.path(output_dir, fileName))){
    message(paste0("A file named '",fileName,"' already exists in the specified directory."))
    message("Are you sure you want to proceed and overwrite existing clean data? (y/n) ")
    question1 <- readLines(n = 1)
    if(regexpr(question1, 'y', ignore.case = TRUE) != 1){
      message("Not overwriting existing files...")
      overwrite <- FALSE
      return(overwrite)
    } else {
      message(paste0("Now overwriting 'Data_", SPECIES, "_", DATE,".html'."))
    }
  }
  
  ##----- 1. LOAD RAW NGS DATA -----
  ##-- Load raw excel file imported from rovbase 
  DNA <- suppressWarnings(readMostRecent( path = data_dir,
                                          extension = ".xls",
                                          pattern = "DNA")) %>%
    ##-- Remove any duplicates
    distinct(., .keep_all = TRUE) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(c( DNAID_sample = "DNAID (Prøve)",
                        Barcode_sample = "Strekkode (Prøve)",
                        RovbaseID_sample = "RovbaseID (Prøve)",
                        EventID = "HendelseID",
                        Analysis_priority = "Analyseprioritet",
                        Species_sample = "Art (Prøve)",
                        Sample_type = "Prøvetype", 
                        Date = "Funnetdato",
                        Mountain_area = "Fjellområde",
                        Sensitivity = "Følsomhet",
                        Release_Date = "Frigivelsesdato", 
                        Locality = "Lokalitet",
                        Location = "Funnsted",
                        North_original = "Nord (opprinnelig)",
                        East_Original = "Øst (opprinnelig)",
                        North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
                        East_UTM33 = "Øst (UTM33/SWEREF99 TM)" ,
                        North_RT90 = "Nord (RT90)",
                        East_RT90 = "Øst (RT90)" ,
                        Coordinate_system = "Koordinatsystem" ,
                        Site_quality = "Stedkvalitet",
                        Collected_by = "Hvem samlet inn",
                        Collector_name = "Samlet selv - Navn",
                        Collector_phone = "Samlet selv - Telefon", 
                        Collector_email = "Samlet selv - E-post",
                        Collector_role = "Samlet selv - Rolle",
                        Collector_other_name = "Annen innsamler - Navn" ,
                        Collector_other_phone = "Annen innsamler - Telefon", 
                        Collector_other_email = "Annen innsamler - E-post",
                        Collector_other_role = "Annen innsamler - Rolle",
                        Tips_name = "Tipser - Navn",
                        Tips_phone = "Tipser - Telefon",
                        Tips_email = "Tipser - E-post",
                        Tips_role = "Tipser - Rolle",
                        Quality_checked = "Kvalitetssikret av feltpersonell",
                        Quality_check_name = "Kvalitetssikrer - navn",
                        Quality_check_orga = "Kvalitetssikrer - Organisasjon",
                        Comments_sample = "Merknad (Prøve)",
                        Last_saved_by_sample = "Sist lagret av (Prøve)",
                        Last_saved_sample = "Sist lagret dato (Prøve)",
                        DNAID = "DNAID (Analyse)",
                        Barcode = "Strekkode (Analyse)",
                        RovbaseID = "RovbaseID (Analyse)",
                        Species = "Art (Analyse)",
                        Sample_status = "Prøvestatus",
                        Id = "Individ",
                        Sex_analysis = "Kjønn (Analyse)",
                        Sex = "Kjønn (Individ)",
                        Method = "Metode",
                        Analyzed_by = "AnalysertAv",
                        Comments = "Merknad (Analyse)",
                        Last_saved_by = "Sist lagret av (Analyse)" ,
                        Last_saved = "Sist lagret dato (Analyse)" ,
                        Municipality_number = "Kommunenummer",
                        Municipality = "Kommune", 
                        County_number = "Fylkenummer",
                        County = "Fylke"))) %>%
    ##-- Add some columns
    dplyr::mutate( 
      ##-- Add "Country" column
      Country_sample = substrRight(County, 3),
      ##-- Change date format
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      ##-- Extract year
      Year = as.numeric(format(Date,"%Y")),
      ##-- Extract month
      Month = as.numeric(format(Date,"%m"))) %>%
  ##-- Filter data
  filter(.,
         ##-- Filter to the focal species
         Species %in% "Bjørn",
         ##-- Filter out samples with 
         Year %in% years) 
  
  ##############################################################################
  numDetDNA <- nrow(DNA)              ## Number of NGS samples in Rovbase
  numIdDNA <- length(unique(DATA$Id)) ## Number of NGS individuals in Rovbase 
  tabDetDNA_sexyear <- table(DNA$Year, DNA$Sex)
  ##############################################################################
  
  
  
  ##----- 2. LOAD RAW DEAD RECOVERY DATA -----
  
  ##-- Load raw excel file imported from rovbase 
  DR <- suppressWarnings(readMostRecent( path = data_dir,
                                         extension = ".xls",
                                         pattern = "dead")) %>%
    ##-- Remove any duplicates
    distinct(., .keep_all = TRUE) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(c( Species = "Art",
                        Death_cause = "Bakgrunn/årsak",
                        Death_method = "Bakgrunn/årsak metode",
                        Death_purpose = "Bakgrunn/årsak formål", 
                        Date = "Dødsdato",
                        Uncertain_date = "Usikker dødsdato",
                        Time_of_death = "Dødstidspunkt",
                        Age_class = "Alder på dødt individ",
                        Age_class_verif = "Aldersklasse verifisert SVA",
                        Juvenile = "Yngling",
                        Sex = "Kjønn", 
                        Age_estimated = "Alder, vurdert",
                        Age = "Alder, verifisert",
                        Age_verif_by = "Alder, verifisert av",
                        Full_weight =  "Helvekt",
                        Slaughter_weight =  "Slaktevekt",
                        CITES = "CITES-nummer",
                        Felling_site_verif = "Kontroll av fellingsted",
                        Tissue_sample = "Vevsprøve tatt",
                        Hunting_date = "Observasjons/Jaktdato",
                        Field_personnel ="Feltpersonell",
                        Control_status = "Kontrollstatus",
                        Assessment = "Vurdering",
                        Location = "Funnsted",
                        North_original = "Nord (opprinnelig)",
                        East_Original = "Øst (opprinnelig)",
                        North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
                        East_UTM33 = "Øst (UTM33/SWEREF99 TM)",
                        North_RT90 = "Nord (RT90)",
                        East_RT90 = "Øst (RT90)",
                        Coordinate_system = "Koordinatsystem",
                        Site_quality = "Stedkvalitet",
                        Approved_date = "Godkjentdato",
                        Outcome = "Utfall", 
                        Counted_off_against_decision = "Regnes av mot vedtak", 
                        Approved_by = "Godkjent av",
                        Sensitivity = "Følsomhet",
                        Release_Date = "Frigivelsesdato", 
                        Id = "Individ", 
                        SVAID = "SVAID",
                        Lansstyrelsen_number = "Länsstyrelsens nr",
                        Last_saved_by = "Sist lagret av",
                        Last_saved =  "Sist lagret dato",
                        Municipality_number =  "Kommunenummer",
                        Municipality = "Kommune", 
                        County_number = "Fylkenummer",
                        County = "Fylke"))) %>%
  ##-- Add some columns
  dplyr::mutate( 
    ##-- Add "Country" column
    Country_sample = substrRight(County, 3),
    ##-- Change date format
    Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
    ##-- Extract year
    Year = as.numeric(format(Date,"%Y")),
    ##-- Extract month
    Month = as.numeric(format(Date,"%m"))) %>%
  ##-- Filter data
  filter(.,
         ##-- Filter to the focal species
         Species %in% "Bjørn",
         ##-- Filter out samples with 
         Year %in% years) 
  
  ##############################################################################
  numDetDR <- nrow(DR)                ## NUmber of DEAD RECOVERIES in Rovbase
  numIdDR <- length(unique(DR$Id))    ## Number of DEAD RECOVERIES ID in Rovbase
  tabDetDR_sexyear <- table(DR$Year, DR$Sex)
  ##############################################################################
  
  
  
  ##----- 3. MERGE BOTH DATA -----
  
  ##-- Merge DNA and dead recoveries files using all shared names columns
  DATA <- merge( DR, DNA, 
                 by = c("Id","RovbaseID","DNAID","Species","Sex","Date","Year",
                        "East_UTM33","North_UTM33","County","Country_sample","Month"),
                 all = TRUE) 
  
  ##############################################################################
  numDetDATA <- nrow(DATA)                ## NUmber of DEAD RECOVERIES in Rovbase
  numIdDATA <- length(unique(DATA$Id))    ## Number of DEAD RECOVERIES ID in Rovbase
  tabDetDATA_sexyear <- table(DATA$Year, DATA$Sex)
  ##############################################################################
  
  
  
  ##----- 4. CLEAN UP DATA -----
  
  ##-- Set monitoring season for the bear
  if(is.null(samplingMonths)){ samplingMonths <- list(4:11) }
  
  
  ##-- For sampling periods spanning over two calendar years (wolf & wolverine)
  ##-- Set all months in given sampling period to the same year
  index <- DATA$Month < unlist(samplingMonths)[1]
  index[is.na(index)] <- FALSE
  DATA$Year[index] <- DATA$Year[index] - 1
  
  ##-- Fix unknown "Id"
  DATA$Id[DATA$Id == ""] <- NA
  
  ##-- Determine Death and Birth Years
  DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
  DATA$RovbaseID <- as.character(DATA$RovbaseID)
  DATA$Death <- NA
  DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
  DATA$Birth <- DATA$Death - DATA$Age
  
  
  ##############################################################################
  noID <- sum(is.na(DATA$Id))              ## Number of samples without ID
  noDate <- sum(is.na(DATA$Year))          ## Number of samples without Date
  noCoords <- sum(is.na(DATA$East_UTM33))  ## Number of samples without Coords  
  ##############################################################################
  
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
  
  

  ##-- Check sex assignments
  ID <- unique(as.character(DATA$Id))
  DATA$Sex <- as.character(DATA$Sex)
  doubleSexID <- IdDoubleSex <- NULL
  counter <- 1
  for(i in 1:length(ID)){
    ##-- Subset data to individual i
    tmp <- DATA$Sex[DATA$Id == ID[i]]
    
    ##-- Number of times individual i was assigned to each sex
    tab <- table(tmp[tmp %in% c("Hunn","Hann")])
    
    ##-- If conflicting sexes (ID identified as both "Hunn" and "Hann")
    if(length(tab) == 2){
      ##-- If ID assigned the same number of times to the 2 sexes, assign to Ukjent
      if(tab[1] == tab[2]){
        DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"
      } else {
        ##-- Otherwise pick the most common sex
        DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
      }
      IdDoubleSex[counter] <- ID[i]
      counter <- counter + 1
    }
    
    ##-- If only one of "Hunn" or "Hann" registered
    if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
    
    ##-- If anything else registered : "Ukjent"
    if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
    
    doubleSexID[i] <- length(tab)
  }#i
  
  
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
  
  
  ##-- Load most recent "flagged" file from HB
  flagged <- readMostRecent( 
    path = data_dir,
    extension = ".csv",
    pattern = "dna_bear_to_remove", 
    fileEncoding = "Latin1") 
  
  
  ##-- Remove flagged samples 
  keep.alive <- !alive$Barcode_sample %in% flagged$Strekkode
  alive <- alive[keep.alive, ]
  keep.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
  dead.recovery <- dead.recovery[keep.dead, ]
  dead.recovery$Missing <- NA
  dead.recovery$Individ <- NA
  
  
  ##-- Turn into sf points dataframe
  alive <- sf::st_as_sf( x = alive,
                         coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633)) 
  
  
  ##-- Intersect and extract country name
  alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
  
  
  ##-- Turn into sf points dataframe
  dead.recovery <- sf::st_as_sf( x = dead.recovery,
                                 coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633))
  
  
  ##-- Intersect and extract country name
  dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]
  
  
  ##-- Issues
  ##-- Multiple deaths
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
  
  ##-- Detections after death
  id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
  ghosts <- unlist(lapply(id.list, function(id){
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
        }##-- FLAG THOSE FOR HENDRIK
      }
    }, silent = TRUE)
    return(out)
  }))
  samples.to.remove <- unlist(ghosts)
  
  ##-- Remove flagged NGS detections after dead recovery
  alive <- alive[!rownames(alive) %in% samples.to.remove, ]
  
  ##-- Save clean data
  fileName <-  paste0("Data_bear_",DATE,".RData")
  save( alive, 
        dead.recovery,
        IdDoubleSex,
        file = file.path(output_dir, paste0("Data_bear_",DATE,".RData")))
  
  ##-- Return list of metrics for the .rmd report
  params <- list()
  return(params)
  
  
}