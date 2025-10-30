##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE OPSCR analysis 2025 
##
## This R script reproduces the OPSCR analysis of the wolverine data as performed
## by RovQuant in 2024 but with the 'rovquantR' package.
##
## This R script performs:
##  1. Initial cleaning of the wolverine NGS data downloaded from RovBase.3.0
##  2. Data preparation for the OPSCR analysis with the 'nimbleSCR' package
##  3. Model fitting using 'nimble' and 'nimbleSCR'
##  4. Post-processing of the MCMC output
##
## NOTES : Two main updates from last year's script ('53.Cleaned2024.R')
##  1. Added a step to remove duplicated GPS tracks (from Asun)
##  2. Added dead recovery states ('recovered dead legal' and 'recovered dead other') 
##     as in the last wolf ('40.F_2024_sf.R') and bear analyses ('Bear_NORWAY_2015-2024.R').
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 21/10/2025
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2025
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway  
##
##------------------------------------------------------------------------------
##
## Notes: 
## This is based on 'rovquantR' beta version 0.1
##   
##------------------------------------------------------------------------------

## CLEAR-UP ENVIRONMENT ------

rm(list = ls())
gc()


## INSTALL 'rovquantR' FROM GITHUB ------

## Ctrl + Shift + F10 (to restart R session)
#devtools::install_github("PierreDupont/rovquantR")


## LOAD REQUIRED LIBRARIES ------

library(rovquantR)
library(nimbleSCR)
library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(readxl)


##------------------------------------------------------------------------------

## I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test.0.1"

##-- General DropBox directory
dir.dropbox <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant"

##-- SOURCE ADDITIONAL FUNCTIONS
# sourceDirectory(dir.function, modifiedOnly = FALSE)
# sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)
# load(file.path(dir.dropbox, "DATA/MISC DATA/age.lookup.table.RData"))
source("C:/My_documents/RovQuant/Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R")
source("C:/My_documents/RovQuant/Source/DoScale.R")



##------------------------------------------------------------------------------

## 1. GENERAL VARIABLES DECLARATION ------

years <- 2014:2023
n.years <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))
species <- "Jerv"
load("C:/My_documents/rovquantR/R/sysdata.rda")
sampling.months <- list(12,1:6)

two.sex <- TRUE
legal.dead <- NULL
overwrite <- FALSE
plot.check = TRUE


##------------------------------------------------------------------------------

##  cleanRovBaseData() ------

# CleanDataNew2sf( 
#   dna_samples = DNA,
#   dead_recoveries = DEAD,
#   species_id = species,
#   country_polygon = COUNTRIES,
#   threshold_month = unlist(sampling.months)[1],
#   keep_dead = T,
#   age.label.lookup = age.lookup.table)


## ------   1. INITIAL CHECKS -----

##-- Make sure directory structure exists
if(two.sex) {
  makeDirectories( path = working.dir,
                   subFolders = c("female","male"),
                   show.dir = TRUE)
} else {
  makeDirectories( path = working.dir,
                   show.dir = TRUE)
}

##-- Species
if(length(species) > 1) {
  stop('This function can only deal with one species at a time... \nPlease, use one of "bear", "wolf", or "wolverine" for the target species.')
}
if(sum(grep("bear", species, ignore.case = T)) > 0|
   sum(grep("bjørn", species, ignore.case = T)) > 0|
   sum(grep("bjorn", species, ignore.case = T)) > 0) {
  SPECIES <- "Brown bear"
  engSpecies <- "bear"
  norSpecies <- c("Bjørn", "BjÃ¸rn")
} else {
  if(sum(grep("wolf", species, ignore.case = T))>0|
     sum(grep("ulv", species, ignore.case = T))>0) {
    SPECIES <- "Gray wolf"
    engSpecies <- "wolf"
    norSpecies <- "Ulv"
  } else {
    if(sum(grep("wolverine", species, ignore.case = T))>0|
       sum(grep("järv", species, ignore.case = T))>0|
       sum(grep("jerv", species, ignore.case = T))>0) {
      SPECIES <- "Wolverine"
      engSpecies <- "wolverine"
      norSpecies <- "Jerv"
    } else {
      SPECIES <- engSpecies <- norSpecies <- species 
    }
  }
}

##-- Years
if(is.null(years)) { years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }

##-- Sampling months
if(is.null(sampling.months)) {
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

##-- Legal mortality patterns
if(is.null(legal.dead)) {
  if (engSpecies == "bear") {
    legal.dead <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")
  } else {
    if (engSpecies == "wolf") {
      legal.dead <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")
    } else {
      if (engSpecies == "wolverine") {
        legal.dead <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")
      } else {
        legal.dead <- ""
      }
    }
  }
}

##-- Renaming list
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
  Lansstyrelsen_number = "Länsstyrelsens nr",
  Last_saved_by_sample = "Sist lagret av (Prøve)",
  Last_saved_sample = "Sist lagret dato (Prøve)",
  Last_saved_by = "Sist lagret av (Analyse)",
  Last_saved = "Sist lagret dato (Analyse)",
  Last_saved_by = "Sist lagret av",
  Last_saved =  "Sist lagret dato",
  Locality = "Lokalitet",
  Location = "Funnsted",
  Mountain_area = "Fjellområde",
  Method = "Metode",
  Municipality_number = "Kommunenummer",
  Municipality = "Kommune",
  North_original = "Nord (opprinnelig)",
  North_RT90 = "Nord (RT90)",
  North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
  Origin = "Opprinnelse",
  Outcome = "Utfall",
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
  # Sex_analysis = "Kjønn (Analyse)",
  # Sex = "Kjønn (Individ)",
  Sex_analysis = "Kjønn (Individ)",
  Sex = "Kjønn (Analyse)",
  Sex = "Kjønn",
  Sex = "Kön",
  Site_quality = "Stedkvalitet",
  SVAID = "SVAID",
  Uncertain_date = "Usikker dødsdato",
  Weight_slaughter = "Slaktevekt",
  Weight_total =  "Helvekt")

##-- List months
months = c("January","February","March","April","May","June",
           "July","August","September","October","November","December")

##-- Load pre-processed habitat shapefiles
data(COUNTRIES, envir = environment()) 

##-- data info
DATE <- getMostRecent(path = data.dir, pattern = "DNA")

##-- Set file name for clean data
fileName <- paste0("CleanData_", engSpecies, "_", DATE, ".RData")

##-- Check that a file with that name does not already exist to avoid overwriting
if(!overwrite) {
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



## ------    2. CLEAN THE DATA -----

## ------     2.1. RAW NGS DATA -----

##-- Load NGS data (46421 entries)
DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                        extension = ".xls",
                                        pattern = "DNA")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>%
  ##-- Turn potential factors into characters 
  dplyr::mutate(across(where(is.factor), as.character)) %>%
  ##-- Initial filters
  dplyr::filter(
    ##-- Remove DEAD entries from the DNA data [HB] ==> 44349 entries
    ##-- ==> removes 2072 entries !!!! 
    !substr(RovbaseID_sample,1,1) %in% "M",
    ##-- Filter to the focal species ==> 44349 entries
    Species %in% norSpecies,
    ##-- Filter out samples with no ID ==> 37776 entries
    !is.na(Id),
    ##-- Filter out samples with no Coordinates ==> 37776 entries
    !is.na(East_UTM33),
    ##-- Filter out samples with no dates  ==> 37776 entries
    !is.na(Date)) %>%
  ##-- Remove any duplicates left over
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
    Year = ifelse( Month < unlist(sampling.months)[1],
                   Year-1,
                   Year),
    ##-- Fix unknown "Id"
    Id = ifelse(Id %in% "", NA, Id),
    ##-- Fix unknown "Sex"
    Sex = ifelse(Sex %in% "Ukjent", "unknown", Sex),
    Sex = ifelse(is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% "Hunn", "female", Sex),
    Sex = ifelse(Sex %in% "Hann", "male", Sex)) 
dim(DNA)
table(DNA$Year)



## ------     2.2. RAW DEAD RECOVERY DATA -----

##-- Load raw excel file imported from rovbase 
DR <- suppressWarnings(readMostRecent( path = data.dir,
                                       extension = ".xls",
                                       pattern = "dead")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>%
  ##-- Turn potential factors into characters 
  dplyr::mutate(across(where(is.factor), as.character)) %>%
  ##-- Initial filters
  dplyr::filter(
    ##-- Remove un-verified dead recoveries [HB] ==> 2646 entries
    ##-- ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
    ## ==> removes 16 entries !!!!
    !grepl( pattern = "Påskutt", x = as.character(Outcome)),
    ##-- Filter to the focal species ==> 2646 entries
    Species %in% norSpecies,
    ##-- Filter out samples with no ID ==> 2371 entries
    !is.na(Id),
    ##-- Filter out samples with no Coordinates ==> 2371 entries
    !is.na(East_UTM33),
    ##-- Filter out samples with no dates ==> 2371 entries
    !is.na(Date)) %>%
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
    Year = ifelse( Month < unlist(sampling.months)[1],
                   Year-1,
                   Year), 
    ##-- Fix unknown "Id"
    Id = ifelse(Id %in% "", NA, Id),
    ##-- Fix unknown "Sex"
    Sex = ifelse(Sex %in% "Ukjent", "unknown", Sex),
    Sex = ifelse(is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% "Hunn", "female", Sex),
    Sex = ifelse(Sex %in% "Hann", "male", Sex),
    ##-- Identify legal deaths
    Legal = grepl(paste(legal.dead, collapse = "|"),
                  Death_cause))
dim(DR)
table(DR$Year)



## ------     2.3. MERGE -----

##-- Merge (no duplicates in the two files at this point) ==> 40147 samples
DATA <- merge( DNA, DR,
               by = c("Id","RovbaseID","DNAID","Species","Sex",
                      "Date","Year","Month",
                      "East_UTM33","North_UTM33",
                      "County","Country_sample"),
               all = TRUE)
dim(DATA)



## ------     2.4. AGE -----

##-- Determine Death and Birth Years
DATA <- DATA %>%
  mutate(
    Age = suppressWarnings(as.numeric(as.character(Age))),
    Death = ifelse(substr(RovbaseID,1,1) %in% "M", Year, NA),
    Birth = Death - Age)

# [PD] This is useless
# ##-- Reconstruct minimal & maximal ages
# DATA$Age.orig <- DATA$Age
# temp <- as.character(base::levels(DATA$Age.orig)) ## list age levels
# temp <- toupper(temp)                         ## Upper case all
# temp <- gsub("\\s", "", temp)                 ## Remove blank spaces
# DATA$Age.orig2 <- DATA$Age.orig
# levels(DATA$Age.orig2) <- temp
# DATA <- merge( DATA, age.lookup.table[ ,-1],
#                  by.x = "Age.orig2",
#                  by.y = "age.label",
#                  all.x = TRUE)                ## Merge with info from lookup table
# 
# ##-- Fill in the rest of the ages from numeric records
# numeric.age.records <- which(!is.na(as.numeric(as.character(DATA$Age.orig2))) & !is.na(DATA$Age.orig2))
# DATA[numeric.age.records, c("min.age","max.age","age")] <- floor(as.numeric(as.character(DATA$Age.orig2[numeric.age.records])))
# dim(DATA)



## ------     2.5. SEX ASSIGNMENT -----

##-- SET THE SEX OF INDIVIDUALS BASED ON ALL INFORMATION AVAILABLE
##-- list all individual IDs
ID <- unique(as.character(DATA$Id))
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
    IdDoubleSex[counter] <- ID[i]
    counter <- counter + 1
  }
  
  ##-- If only one of "female" or "male" registered
  if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
  
  ##-- If anything else registered : "unknown"
  if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "unknown"}
  
  doubleSexID[i] <- length(tab)
}#i
table(DATA$Sex)

# ##-- Remove individuals with unknwon sex
# DATA <- DATA[which(DATA$Sex %in% c("Hann","Hunn")), ]

##-- Turn into sf points dataframe
DATA <- sf::st_as_sf( x = DATA,
                      coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(., sf::st_crs(32633)) 

##-- Remove all samples outside the polygon of interest (study area)
## [PD] : MIGHT WANT TO KEEP SAMPLES OUTSIDE NOR & SWE IN cleanRovbaseDat
DATA <- DATA[!is.na(as.numeric(st_intersects(DATA, st_union(COUNTRIES)))), ]

##-- Intersect and extract country name
DATA$Country_sf[!is.na(as.numeric(sf::st_intersects(DATA, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
DATA$Country_sf[!is.na(as.numeric(sf::st_intersects(DATA, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"



## ------     2.6. SPLIT DATA -----

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



## ------   3. SPECIES-SPECIFIC CLEANING STEPS ------

## ------     3.1. WOLVERINE -----

## ------       3.1.1. FLAGGED SAMPLES -----

##-- REMOVE SUSPECT NGS SAMPLES ACCORDING TO HENRIK
SUSPECT_NGS_SAMPLES <- readMostRecent(
  path = data.dir,
  extension = ".xlsx",
  pattern = "Remove ngs")

##-- ==> removes 866 entries !!!!
alive <- alive %>%
  filter(!DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB))
dim(alive)

##-- REMOVE SUSPECT DEAD RECOVERIES ACCORDING TO HENRIK
SUSPECT_DeadRecoSAMPLES <- readMostRecent(
  path = data.dir,
  extension = ".xlsx",
  pattern = "Remove dead")

##-- ==> removes 10 entries !!!!
dead.recovery <- dead.recovery %>%
  filter(!RovbaseID %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID))
dim(dead.recovery)

##-- REMOVE ADDITIONAL DEAD RECOVERIES FLAGGED BY HENRIK (email from the 18/12/2024)
dead.recovery <- dead.recovery %>% 
  filter(!RovbaseID %in% c("M495994","M524051","M524052","M524053"))
dim(dead.recovery)



## ------       3.1.2. PUPS & YOUNG INDIVIDUALS -----

## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and 
## recovered dead between March and November
dead.recovery <- dead.recovery %>%
  filter(!(Age_class %in% "Unge" & Month > 2 & Month < 12))
dim(dead.recovery)

## 2) remove individuals that have a weight >0 and <4 between March and November
##-- FORMAT WEIGHTS
dead.recovery <- dead.recovery %>%
  mutate(
    Weight_total = as.numeric(gsub(",", ".", as.character(Weight_total))),
    Weight_slaughter =  as.numeric(gsub(",", ".", as.character(Weight_slaughter))),
    weight = ifelse(!is.na(Weight_total), Weight_total, Weight_slaughter),
    weight = ifelse(is.na(weight), -999, weight))

## Check with Henrik
## This step does not remove dead recoveries on id with weight == 0, should it?
##-- WEIGTH DISTRIBUTION
par(mfrow = c(4,3))
for(t in 1:12){
  hist(dead.recovery$weight[(dead.recovery$weight > -1) &
                              dead.recovery$Month %in% t],
       breaks = 0:20,
       main = months[t],
       xlab = "Dead recovery Weight (kg)")
}#t

##-- AGE DISTRIBUTION
par(mfrow = c(4,3))
for(t in 1:12){
  hist(dead.recovery$Age[(dead.recovery$Age > -1) &
                           dead.recovery$Month %in% t],
       breaks = 0:20,
       main = months[t],
       xlab = "Age at recovery")
}#t

##-- check how many dead reco we remove and remove if more than 0
dead.recovery <- dead.recovery %>%
  filter(!(weight > 0 & weight < 4 & Month < 12 & Month > 2))
dim(dead.recovery)

##-- check how many dead reco with a weight of 0 kg and recovered between March and November
sum(dead.recovery$Age %in% 0 & dead.recovery$Month < 12 & dead.recovery$Month > 2)



## ------   4. DATA ISSUES -----

## ------     4.1. MULTIPLE DEATHS ------

##-- Remove individuals that died twice
IdDoubleDead <- dead.recovery$Id[duplicated(dead.recovery$Id)]
if(length(IdDoubleDead) > 0){
  duplicatedDeath <- NULL
  for(i in IdDoubleDead){
    tmp <- which(dead.recovery$Id == i & is.na(dead.recovery$Death_method))
    if(length(tmp)==0){tmp <- which(dead.recovery$Id == i)[-2]} ##[CM] remove the second record.
    duplicatedDeath <- c(duplicatedDeath, tmp)
  }#i
  dead.recovery <- dead.recovery[-duplicatedDeath, ]
}#if




## ------   5. TURN INTO .sf OBJECTS -----

##-- Turn into sf points dataframe
alive <- sf::st_as_sf( x = alive,
                       coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(., sf::st_crs(32633)) 

##-- Intersect and extract country name
alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"


##-- Turn into sf points dataframe
dead.recovery <- sf::st_as_sf( x = dead.recovery,
                               coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(.,sf::st_crs(32633))

##-- Intersect and extract country name
dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"



## ------   6. SAVE DATA ------
save( alive, 
      dead.recovery,
      file = file.path( working.dir, "data", fileName))


## END OF cleanRovBaseData() ------ 

##-- Checks
dim(alive)
table(alive$Year)
dim(dead.recovery)
table(dead.recovery$Year)
table(dead.recovery$Death_cause)
table(dead.recovery$Death_method)


##------------------------------------------------------------------------------

## makeRovquantData_wolverine() ------

## ------ 0. BASIC SET-UP ------

##-- Set default values for the wolverine model
aug.factor <- 0.8
sampling.months <- list(12,1:6)
habitat.res <- 20000
buffer.size <- 60000
detector.res <- 10000
subdetector.res <- 2000
max.det.dist <- 40000
resize.factor <- 1

##-- Set up list of Habitat characteristics
habitat <- list( resolution = habitat.res,
                 buffer = buffer.size)

##-- Set up list of Detectors characteristics
detectors <- list( resolution = detector.res,
                   resolution.sub = subdetector.res,
                   maxDist = max.det.dist,
                   resize.factor = resize.factor)

##-- Set up list of Data characteristics
data <- list( aug.factor = aug.factor,
              sampling.months = sampling.months)



## ------ I. LOAD & SELECT DATA ------

## ------   1. HABITAT DATA ------

##-- Load pre-defined habitat rasters and shapefiles
data(GLOBALMAP, envir = environment()) 
data(COUNTRIES, envir = environment()) 
data(REGIONS, envir = environment())
data(COUNTIES, envir = environment()) 
data(habitatRasters, envir = environment()) 

##-- Disaggregate habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)

##-- Merge counties for practical reasons
COUNTIES_AGGREGATED <- REGIONS %>%
  mutate(id = case_when(
    county %in% c("Norrbotten") ~ 1,
    county %in% c("Västerbotten") ~ 2,
    county %in% c("Blekinge","Dalarna","Gävleborg","Gotland","Halland","Jämtland",
                  "Jönköping","Kalmar","Kronoberg","Örebro","Östergötland","Skåne",
                  "Södermanland","Stockholm","Uppsala","Värmland","Västernorrland",
                  "Västmanland","Västra Götaland") ~ 3,
    county %in% c("Agder","Akershus","Buskerud","Innlandet","Møre og Romsdal",
                  "Oppland","Oslo","Østfold","Rogaland","Vestland","Telemark",
                  "Vestfold") ~ 4,
    county %in% c("Trøndelag") ~ 5,
    county %in% c("Finnmark") ~ 6,
    county %in% c("Nordland") ~ 7,
    county %in% c("Troms") ~ 8)) %>%
  group_by(id) %>%
  summarize() %>%
  st_simplify( ., preserveTopology = T, dTolerance = 500)
plot(COUNTIES_AGGREGATED)



## ------   2. STUDY AREA ------

##-- CREATE STUDY AREA POLYGON 
myStudyArea <- COUNTRIES %>%
  filter(ISO %in% c("NOR","SWE")) %>%
  mutate(id = 1) %>%
  group_by(id) %>% 
  summarize()

# ##-- CREATE HABITAT POLYGON 
# myBufferedArea <- myStudyArea %>%
#   st_buffer(dist = HABITAT$habBuffer) %>%
#   st_intersection(., GLOBALMAP)
# 
# ##-- PLOT CHECK
# if(plot.check){
#   par(mfrow = c(1,1))
#   plot(st_geometry(COUNTRIES))
#   plot(st_geometry(myStudyArea), add = TRUE, col = "red")
# }



## ------   3. NGS DATA -----

##-- Extract date from the last cleaned data file
DATE <- getMostRecent( 
  path = file.path(working.dir,"data"),
  pattern = "CleanData_wolverine")

##-- Load the most recent Bear data from RovBase
myFullData.sp <- readMostRecent( 
  path = file.path(working.dir,"data"),
  pattern = "CleanData_wolverine",
  extension = ".RData")

##-- Define years
if(is.null(years)){
  years <- sort(unique(c(myFullData.sp$alive$Year,
                         myFullData.sp$dead.recovery$Year)))
}
data$years <- years
n.years <- length(years)



## ------     3.1. FILTER NGS & DEAD RECOVERY DATA FOR DATES ------

myFilteredData.sp <- myFullData.sp

##-- Filter for dates
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  dplyr::filter(
    ##-- Subset to years of interest
    Year %in% years,
    ##-- Subset to monitoring period 
    Month %in% unlist(sampling.months))

##-- Filter for dates
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery %>%
  ##-- Subset to years of interest
  dplyr::filter(Year %in% years)

##-- Checks
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$dead.recovery)
table(myFilteredData.sp$dead.recovery$Year)
table(myFilteredData.sp$dead.recovery$Death_cause)
table(myFilteredData.sp$dead.recovery$Death_method)



## ------     3.2. FILTER OUT DETECTIONS IN NORRBOTTEN EXCEPT IN 2016:18 and 2023 ------ 

##-- list years with or without sampling in Norrbotten
yearsSampledNorrb <- c(2016:2018,2023)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)

##-- Get Norrbotten borders
COMMUNES_NOR <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp")) ## Communal map of Norway
COMMUNES_SWE <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp")) ## Communal map of Sweden
COUNTIESNorrbotten <- rbind(COMMUNES_NOR, COMMUNES_SWE) %>%
  filter(NAME_1 %in% "Norrbotten") %>%
  group_by(NAME_1) %>%
  summarize()

## [PD] : USING REGIONS INSTEAD OF COMMUNES WILL LEAD TO DIFFERENT NUMBERS OF 
## SAMPLES REMOVED ON DIFFERENT YEARS BECAUSE OF SLIGHT DIFFERENCES IN SHAPEFILES
## (APPROX. 10 samples OVERALL)
# COUNTIESNorrbotten <- REGIONS %>%
#   filter(county %in% "Norrbotten") %>%
#   group_by(county) %>%
#   summarize()

##-- Check how many detections have been collected in Norrbotten overall
is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))
sum(is.Norr, na.rm = T)

##-- Check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(is.Norr)), ]$Year) %>% sum()

##-- Subset NGS dataset
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  filter(!(Year %in% yearsSampledNorrb & !is.na(is.Norr)))

dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)

# ## plot check
# for(t in 1:n.years){
#   plot( st_geometry(myStudyArea))
#   plot( st_geometry(COUNTIESNorrbotten), add = T, col = "blue")
#   plot( st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t], ]),
#         col = "red", add = T, pch = 16)
# }



## ------     3.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

## ------       3.3.1. CLEAN GPS TRACKS ------

message("Cleaning GPS tracks... ")

## LOAD GPS SEARCH TRACKS
# TRACKS_SINGLE <- read_sf(file.path(data.dir,
#                                    "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp"))
# TRACKS_MULTI <- read_sf(file.path(data.dir,
#                                   "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp"))
TRACKS_SINGLE <- read_sf(file.path(data.dir,
                                   "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp"))
TRACKS_MULTI <- read_sf(file.path(data.dir,
                                  "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))

## COMBINE ALL TRACKS AND FIX DATES
TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI) %>%
  mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
          Yr = as.numeric(format(Dato,"%Y")),
          Mth = as.numeric(format(Dato,"%m")),
          Dato = as.character(Dato),
          Length = st_length(., byid = T),
          Centroidx = st_coordinates(st_centroid(.))[ ,1]) %>%
  dplyr::filter(
    ## REMOVE HELICOPTER TRACKS
    Helikopter == "0",
    ## KEEP ONLY WOLVERINE TRACKS
    Jerv == "1")

## FIND DUPLICATES BASED ON PERSON, DISTANCE and DATE
df <- data.frame( Dato = TRACKS$Dato,
                  Person = TRACKS$Person,
                  dist = TRACKS$Length,
                  centroidx = TRACKS$Centroidx)
dupIDs <- TRACKS$ID[duplicated(df)]
dupLength <- TRACKS$Length[duplicated(df)]

## STORE CLEAN TRACKS IN A LIST
TRACKS <- TRACKS[-duplicated(df), ]


# ## Check that we only have wolverine tracks
# table(ALL_TRACKS$Jerv) == nrow(ALL_TRACKS)
# ## Check
# hist(ALL_TRACKS$Yr)
# 
# ## GET EXTENT
# myStudyArea.extent <- st_bbox(extent(myStudyArea))
# st_crs(myStudyArea.extent) <- st_crs(COUNTRIES)
# 
# ## [PD] ASP ADDED REMOVAL OF DUPLICATED TRACKS
# dupIDs <- dupDist <- TRACKS_YEAR <- list()
# for(t in 1:nYears){
#   TRACKS <- ALL_TRACKS %>% 
#     ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
#     filter( Yr%in%YEARS[[t]][1] & Mth %in% myVars$DATA$samplingMonths[[1]] |
#               Yr%in%YEARS[[t]][2] & Mth%in%myVars$DATA$samplingMonths[[2]]) %>%
#     ## SUBSET TRACKS TO THE STUDY AREA
#     st_intersection(., st_as_sfc(myStudyArea.extent)) 
#   
#   ## NAME TRACKS
#   TRACKS$ID <- 1:nrow(TRACKS)
#   ## CALCULATE LENGTH OF EACH TRACK TO IDENTIFY DUPLICATES
#   TRACKS$dist <- st_length(TRACKS, byid = T)
#   ## CALCULATE CENTROIDS TO AVOID KEEPING TRACKS WITH THE SAME LENGTHS BUT IN DIFFERENT LOCATIONS
#   TRACKS$centroidx <- st_coordinates(st_centroid(TRACKS))[ ,1]
#   
#   ## FIND DUPLICATES BASED ON PERSON, DISTANCE and DATE
#   df <- data.frame( Dato = TRACKS$Dato,
#                     Person = TRACKS$Person,
#                     dist = TRACKS$dist,
#                     centroidx = TRACKS$centroidx)
#   dupIDs[[t]] <- TRACKS$ID[duplicated(df)]
#   dupDist[[t]] <- TRACKS$dist[duplicated(df)]
#   
#   ## STORE CLEAN TRACKS IN A LIST
#   TRACKS_YEAR[[t]] <- TRACKS[-dupIDs[[t]], ]
#   
#   # try a fast way to identify duplicated tracks
#   # turn to dataframe and identify them
#   # sub_tracks_filter <- TRACKS_YEAR[[t]] %>%
#   #   distinct(Dato, dist, .keep_all = T)
#   # distinct(Person, Dato, dist, .keep_all = T)
# }#t
# 
# ## PLOT CHECK
# if(myVars$plot.check){
#   par(mfrow = c(2,2))
#   ## Length of tracks searched per year
#   lengthPerYear <- unlist(lapply(TRACKS_YEAR,function(x) sum(x$dist)/1000))
#   names(lengthPerYear) <- years
#   barplot(lengthPerYear, ylab = "Track length (km)", main = "Length of tracks searched per year")
#   
#   ## Number of tracks searched per year
#   numPerYear <- unlist(lapply(TRACKS_YEAR,function(x) length(unique(x$ID))))
#   names(numPerYear) <- years
#   barplot(numPerYear, ylab = "Number of tracks", main = "Number of tracks searched per year")
#   
#   ## Length of tracks duplicated per year
#   dupdist <- unlist(lapply(dupDist,function(x) sum(x)/1000))
#   names(dupdist) <- years
#   barplot(dupdist,ylab = "Track length (km)", main = "Length of tracks duplicated per year")
#   
#   ## Number of tracks duplicated per year
#   dup <- unlist(lapply(dupIDs,length))
#   names(dup) <- years
#   barplot(dup, ylab = "Number of tracks", main = "Number of tracks duplicated per year")
# }



## ------       3.3.1. ASSIGN SAMPLES TO TRACKS ------

### [PD] need to check if this has already been done !!! 
### Need to think about the best way to implement it

message("Assigning DNA samples to GPS tracks... ")
message("This can take several minutes... ")

myFilteredData.sp$alive <- assignSearchTracks(
  data = myFilteredData.sp$alive,
  tracks = TRACKS)

##-- SAVE FOR FASTER LOADING
save(myFilteredData.sp, file = file.path(working.dir, "data/myFilteredData.sp.RData"))
load(file.path(working.dir, "data/myFilteredData.sp.RData"))


# ##-- OLD STUFF
# ##-- ASSIGN ROVBASE ID
# myFilteredData.sp$alive$trackID <- NA
# myFilteredData.sp$alive$trackDist <- NA
#
# ##-- ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
# dnatemp <- st_as_sf(myFilteredData.sp$alive)
# ##-- CREATE A BUFFER AROUND EACH DETECTION
# tmp <- st_buffer(dnatemp, dist = 750)
# ##-- Loop over detections
# for(i in 1:1000){#nrow(myFilteredData.sp$alive)){
#   t <- which(years %in% tmp[i, ]$Year)
#   ##-- If no tracks on the same date ==> next sample
#   whichSameDate <- which(as.character(TRACKS_YEAR[[t]]$Dato) == as.character(myFilteredData.sp$alive$Date[i]))
#   if(length(whichSameDate) == 0){next}
#   ##-- If no tracks on the same date within 750m of the NGS sample ==> next sample
#   tmpTRACKS <- st_intersection(TRACKS_YEAR[[t]][whichSameDate, ], tmp[i, ])
#   if(nrow(tmpTRACKS) == 0){next}
#   ##-- Calculate distance to matching tracks
#   dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
#   ##-- Assign sample to the closest matching track
#   myFilteredData.sp$alive$trackID[i] <- tmpTRACKS$RovbaseID[which.min(dist)]
#   myFilteredData.sp$alive$trackDist[i] <- min(dist)
#   print(i)
# }#i



## ------       3.3.2. ASSIGN SAMPLES TO OPPORTUNISTIC OR STRUCTURED ------

distanceThreshold <- 500

# HairTrapSamples <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/hairtrapsNB2024.xlsx"))#, fileEncoding="latin1") ## DNA samples to be removed from Henrik
HairTrapSamples <- readMostRecent( 
  path = data.dir,
  extension = ".xlsx",
  pattern = "hairtrap")

##-- Identify samples from structured and opportunistic sampling
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  mutate( 
    ##-- Collector column was replaced by two columns, merging them now...
    Collector_role <- ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
    ##-- Identify samples collected by hair traps
    hairTrap = DNAID %in% HairTrapSamples$DNAID,
    structured = Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") & 
      !is.na(trackID) &  
      trackDist <= distanceThreshold)
table(myFilteredData.sp$alive$Collector, useNA = "always")
table(myFilteredData.sp$alive$structured, useNA = "always")

##-- PLOT CHECK
if(plot.check){
  plot(REGIONS[REGIONS$county %in% "Norrbotten", ]$geometry)
  tmp <- myFilteredData.sp$alive %>%  filter(hairTrap)
  plot(tmp$geometry, add = T, col = "red", pch = 16)
  if(length(which(myFilteredData.sp$alive$DNAID[myFilteredData.sp$alive$structured] %in% HairTrapSamples$DNAID))>0){
    print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
  }
}



## ------       3.3.3. PLOT CHECKS ------

if(plot.check){
  
  ##-- Barplot of structured vs. opportunistic samples
  pdf(file = file.path(working.dir, "figures/DetectionsStructuredOppBarplot.pdf"))
  par(mfrow = c(2,1), mar = c(4,4,3,2))
  barplot( rbind(table(myFilteredData.sp$alive$Year[myFilteredData.sp$alive$structured]),
                 table(myFilteredData.sp$alive$Year[!myFilteredData.sp$alive$structured])),
           beside = T,
           ylim = c(0,2000),
           col = c(grey(0.2),grey(0.8)),
           ylab = "Number of samples")
  abline(h = seq(0, 2000, by = 500),
         lty = 2, col = grey(0.8))
  title(main = "500m threshold")
  legend("topleft", fill = c(grey(0.2),grey(0.8)), legend = c("Structured","Other"))
  
  ##
  structured2000 <- myFilteredData.sp$alive$Collector %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
    !is.na(myFilteredData.sp$alive$trackID) &
    myFilteredData.sp$alive$trackDist <= 2000
  barplot( rbind(table(myFilteredData.sp$alive$Year[structured2000]),
                 table(myFilteredData.sp$alive$Year[!structured2000])),
           beside = T,
           ylim = c(0,2000),
           col = c(grey(0.2),grey(0.8)),
           ylab = "Number of samples")
  abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
  title(main="2000m threshold")
  legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))
  dev.off()
  
  
  ##-- CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
  tmp <- myFilteredData.sp$alive %>%
    filter(Collector %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"))
  
  # plot check 
  pdf(file = file.path(working.dir, "figures/DetectionsStructuredOpp.pdf"))
  for(t in 1:n.years){
    par(mar = c(0,0,3,0), mfrow = c(1,3))
    
    ##-- samples with tracks
    tmpTracks <- tmp %>% 
      filter( Year %in% years[t],
              !is.na(trackID),
              trackDist <= 750)
    plot( st_geometry(myStudyArea), col = "gray60", main = "Structured with track")
    plot( st_geometry(tmpTracks),
          pch = 21, col = "black",
          cex = 1, bg = "red", add = T)
    
    ##-- samples without tracks
    tmpNoTracks <- tmp %>% 
      filter( Year %in% years[t],
              is.na(trackID) | trackDist >= 750)
    plot( st_geometry(myStudyArea), col = "gray60", main = "Structured without track")
    plot( st_geometry(tmpNoTracks),
          pch = 21, col = "black", 
          cex = 1, bg = "blue", add = T)
    
    ##-- Other samples
    tmpOpp <- myFilteredData.sp$alive %>%
      filter(!Collector %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),
             Year%in% years[t])
    plot( st_geometry(myStudyArea), col = "gray60", main = "Other samples")
    plot( st_geometry(tmpOpp),
          pch = 21, col = "black",
          cex = 1, bg = "green", add = T)
    
    mtext(years[t], adj = -0.8, padj = 1)
  }#t
  
  ##-- Number of samples collected per year
  tab <- table(tmp$Year, tmp$trackID, useNA ="always")
  barplot( tab[ ,which(is.na(colnames(tab)))]/rowSums(tab),
           main = "% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 
  dev.off()
  
  
  
  ##-- plot check 
  pdf( file = file.path(working.dir, "figures/OverallDetectionsDeadRecoveries.pdf"))
  plot( st_geometry(GLOBALMAP))
  plot( st_geometry(myStudyArea), add = T)
  plot( st_geometry(myFullData.sp$alive),
        pch = 16, col = "red", cex = 0.3, add = T)
  plot( st_geometry(myFullData.sp$dead.recovery),
        pch = 16, col = "blue", cex = 0.3, add = T)
  mtext(paste("Live detections", nrow(myFullData.sp$alive),
              "; ID:", nrow(unique(myFullData.sp$alive$Id))),
        line = +1)
  mtext(paste("Dead recovery:", nrow(myFullData.sp$dead.recovery)))
  dev.off()
}



## ------     3.4. SEPARATE MORTALITY CAUSES ------ 

##-- Identify legal mortality causes
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
whichLegalCauses <- unlist(lapply(c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske"), 
                                  function(x)grep(x,MortalityNames)))
legalCauses <- MortalityNames[whichLegalCauses]

##-- Split data based on mortality causes
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery %>% 
  mutate(legal = Death_cause %in% legalCauses)
# legal.death <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]
# Other.death <- myFilteredData.sp$dead.recovery[!myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]


##-- PLOT CHECK
if(plot.check){
  par(mfrow = c(1,3))
  for(t in 1:n.years){
    ## DEAD RECOVERIES TOTAL
    tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
    NGS_TabTotal <- table(tempTotal$Country)
    ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
    ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
    tempIn <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
    NGS_TabIn <- table(tempIn$Country)
    ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
    ## PLOT NGS SAMPLES
    plot(st_geometry(GLOBALMAP), col="gray80")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
    plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
    ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
    ## ADD OVERALL NUMBERS
    mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
    mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
    # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
  }#t
  
  
  ## PLOT TREND DETECTIONS & DEAD RECOVERIES OVER TIME AND SPACE 
  
  ## DETECTIONS
  pdf(file = file.path(working.dir, "figures/TRENDDetections.pdf"))
  
  temp <- unique(myFilteredData.sp$alive[ ,c("Year","Country","DNAID")])
  tab_Country.Year <- table(temp$Year, temp$Country)
  country.colors <- c("goldenrod1","goldenrod3")
  
  par(mfrow = c(1,1), mar = c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
  lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  ## ID DETECTED
  temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country,myFilteredData.sp$alive$Id)
  tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
  country.colors <- c("goldenrod1","goldenrod3")
  par(mfrow = c(1,1), mar = c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
  lines(tab_Country.Year1[,"N"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year1[,"S"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  ## Average number of detection per detected ID  
  tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
  par(mfrow = c(1,1), mar = c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
       ylab="Average Number of detections", xlab="Years")
  lines(tab_Country.Year2[,"N"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year2[,"S"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  ## deadrecovery 
  temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")])
  tab_Country.Year <- table(temp$Year, temp$Country)
  country.colors <- c("goldenrod1","goldenrod3")
  
  par(mfrow = c(1,1), mar = c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
       ylab="N Id Dead recovered", xlab="Years")
  lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("topright",c("N","S"), fill=country.colors)
  dev.off()
}



## ------ II. CREATE OPSCR DATA ------

## ------   1. GENERATE HABITAT ------

message("Preparing habitat characteristics... ")

## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------

##-- Determine study area based on NGS detections
##-- Buffer NGS detections and cut to Swedish and Norwegian borders
studyArea <- myFullData.sp$alive %>%
  sf::st_buffer(., dist = habitat$buffer * 1.4) %>%
  sf::st_union() %>%
  sf::st_intersection(., COUNTRIES) %>%
  sf::st_as_sf()

##-- make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster(
  poly = studyArea,
  habitat.r = habRaster,
  buffer = habitat$buffer,
  plot.check = FALSE) %>%
  append(habitat,.)

habitat.xy <- coordinates(habitat$habitat.r)[habitat$habitat.r[ ]==1, ] 
n.habCells <- nrow(habitat.xy)[1]

# ## [PD] USELESS!
# ##-- Retrieve habitat windows boundaries
# lowerHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1, ] - 0.5*habitat$resolution
# upperHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1, ] + 0.5*habitat$resolution
# nHabCells <- dim(lowerHabCoords)[1]
# 
# # ## [PD] USELESS !
# ##-- CREATE HABITAT GRID
# habIDCells.mx <- habitat$IDCells.mx
# habIDCells.mx[] <- 0
# for(i in 1:nrow(lowerHabCoords)){
#   habIDCells.mx[trunc(lowerHabCoords[i,2])+1,
#                 trunc(lowerHabCoords[i,1])+1] <- i
# }
# # image(habIDCells.mx)
#
# ## [PD] USELESS!
# whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyArea))))
# if(length(whichOut) > 0){
#   myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
# }
# #myFilteredData.sp$alive$Id <- droplevels( myFilteredData.sp$alive$Id)
# ## REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT #[CM] 
# whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, habitat$buffered.habitat.poly))))
# if(length(whichOutBuff) > 0){ myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]}

##-- PLOT CHECK
if(plot.check){
  # par(mfrow = c(1,2))#[CM]
  plot(habitat$habitat.r)
  plot(st_geometry(myStudyArea),
       add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
  plot(st_geometry(GLOBALMAP), add = T)
  plot(st_geometry(habitat$buffered.habitat.poly), add = T)
  plot(st_geometry(myFilteredData.sp$alive),
       pch = 21, bg = "red", cex = 0.5, add = T)
  plot(st_geometry(myFilteredData.sp$dead.recovery),
       pch = 21, bg = "blue", cex = 0.5, add = T)
  
  ##-- Check correlation number of detections ~ between monitoring season
  deadID <- unique(myFilteredData.sp$dead.recovery$Id)
  ndet <- NULL
  timeDiff <- NULL
  for(i in 1:length(deadID)){
    tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i], ]$Year
    
    timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i], ]$Date-
      as.POSIXct(strptime(paste("01-12", tmpYear, sep = "-"), "%d-%m-%Y")) 
    
    ndet[i] <- nrow(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] & 
                                              myFilteredData.sp$alive$Year %in% tmpYear, ])
  }
  
  pdf(file = file.path(working.dir, "figures/Prop id detected_Time available.pdf"))
  plot( ndet ~ timeDiff,
        ylab = "Total number of detections",
        xlab = "Number of days between dec 1 and dead recovery")
  hh <- hist(timeDiff[ndet > 0], breaks = seq(0,400,by=25))
  hh1 <- hist(timeDiff[ndet == 0], breaks = seq(0,400,by=25))
  barplot(rbind(hh$counts/(hh$counts+hh1$counts),
                hh1$counts/(hh$counts+hh1$counts)),
          names.arg = hh$breaks[1:(length(hh$breaks)-1)],
          xlab = "number of days between dead reco and start monitoring",
          ylab = "%")
  legend( "topright",
          fill = c(grey(0.2), grey(0.8)),
          legend = c("detected","notDetected"))
  dev.off()
}



## ------     1.3. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       1.3.1. DEN COUNTS ------

##-- Load the last DEN COUNT data file
#DEN <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/DEN_COUNTS_2009_2024_fromHB.csv"), fileEncoding="latin1")
DEN <- readMostRecent( 
  path = data.dir,
  extension = ".csv",
  pattern = "DEN_COUNTS") %>%
  st_as_sf(., coords = c("UTM33_X", "UTM33_Y")) %>%
  st_set_crs(value = st_crs(myFilteredData.sp$alive))
# colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN))
# DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
# st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
# # DEN.sp$id  <- rep(1, nrow(DEN.sp))
# # DEN.sp <- DEN.sp[ ,("id")]

DEN.r <- raster(
  adehabitatHR::estUDm2spixdf(
    adehabitatHR::kernelUD( as(DEN, "Spatial"),
                            h = 30000,
                            grid = as(habitat$habitat.r, 'SpatialPixels'))))

##-- PLOT CHECK
if(plot.check){
  plot(DEN.r)
  plot(st_geometry(myStudyArea), add = TRUE, border = "black")
}

##-- EXTRACT COVARIATES
denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
denCounts <- round(scale(denCounts), digits = 2)



## ------   6. GENERATE DETECTORS -----

message("Preparing detectors characteristics... ")

## ------     6.1. GENERATE DETECTORS CHARACTERISTICS -----

##-- GENERATE NGS DETECTORS BASED ON THE STUDY AREA
subdetectors.r <- disaggregate(
  habitat$habitat.rWthBuffer,
  fact = res(habitat$habitat.r)[1]/detectors$resolution.sub)

##-- Generate NGS detectors based on the raster of sub-detectors
detectors <- makeSearchGrid( 
  data = subdetectors.r,
  resolution = detectors$detResolution,
  div = (detectors$resolution/detectors$resolution.sub)^2,
  plot = FALSE) %>%
  append(detectors,.)

##-- EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(detectors$main.detector.sp)[1]

##-- FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(detectors$main.detector.sp)
colnames(detector.xy) <- c("x","y")
n.trials <- as.vector(table(detectors$detector.sp$main.cell.id))

##-- IDENTIFY DETECTORS IN NORBOTTEN 
COUNTIESAroundNorrbotten <- REGIONS %>%
  group_by(county) %>%
  summarize() %>%
  filter(county %in% c("Norrbotten","Troms","Västerbotten","Nordland","Finnmark")) %>% 
  st_simplify( dTolerance = 500)

##-- CREATE A NORROBOTTEN DETECTOR GRID
distDetsCounties <- st_distance( detectors$main.detector.sp,
                                 COUNTIESAroundNorrbotten,
                                 byid = T)
detsNorrbotten <- which(apply(distDetsCounties, 1, which.min) == 3)


## PLOT CHECK
if(plot.check){
  plot( st_geometry(COUNTIESAroundNorrbotten))
  plot( st_geometry(detectors$main.detector.sp),
        col = "black", pch = 16, cex = 0.3, add = T)
  plot( st_geometry(detectors$main.detector.sp[detsNorrbotten, ]),
        col = "red", pch = 16, cex = 0.5, add = T)
  
  ## PLOT NGS DETECTORS
  plot( st_geometry(habitat$buffered.habitat.poly),
        main = paste(n.detectors, "Detectors"),
        col = rgb(0.16,0.67,0.16, alpha = 0.3))  
  plot( st_geometry(myStudyArea), add = TRUE,
        col = rgb(0.16,0.67,0.16, alpha = 0.5))
  plot( st_geometry(detectors$main.detector.sp),
        col = "red", pch = 16, cex = 0.1, add = TRUE)
  plot( st_geometry(COUNTRIES), add = TRUE)
}



## ------     6.2. GENERATE DETECTOR-LEVEL COVARIATES -----

## ------       6.4.1. EXTRACT COUNTIES ------

##-- Extract closest county for each detector
detCounties <- detectors$main.detector.sp %>%
  st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
  apply(., 1, function(x) which.min(x))

##-- PLOT CHECK 
if(plot.check){
  myCol <- terrain.colors(nrow(COUNTIES_AGGREGATED))
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(COUNTRIES), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(detectors$main.detector.sp[detCounties %in% 5, ]),
        col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
  plot(st_geometry(detectors$main.detector.sp),
       col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
  plot(st_geometry(COUNTIES_AGGREGATED), add = TRUE)
  plot(st_geometry(detectors$main.detector.sp[detCounties %in% 1, ]),
       col = "red", pch = 16, cex = 0.8, add = T)
}

##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten 
##-- in years without sampling
countyToggle <- matrix(1, nrow = max(detCounties), ncol = n.years)
for(t in yearsNotSampled){
  countyToggle[1,t] <- 0
}



## ------       6.4.2. EXTRACT COUNTRIES ------

##-- Extract closest country for each detector
detCountries <- detectors$main.detector.sp %>%
  st_distance(., COUNTRIES, by_element = F ) %>%
  apply(., 1, function(x) which.min(x)) %>%
  as.factor(.) %>%
  as.numeric(.)

##-- Turn into a matrix
detCountries <- matrix( detCountries,
                        nrow = length(detCountries),
                        ncol = n.years)

##-- Add another category to detCountry if in Norrbotten, to turnoff detection to 0 there. 
for(t in yearsNotSampled){
  detCountries[detCounties %in% 1,t] <- 3
}#t  

##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten in years without sampling
countryToggle <- matrix(1, nrow = max(detCountries), ncol = n.years)
for(t in yearsNotSampled){
  countryToggle[3,t] <- 0
}

##-- PLOT CHECK 
if(plot.check){
  par(mfrow = c(1,1))
  myCol <- c("blue4", "yellow1", "red")
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Countries")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(detectors$main.detector.sp), col = myCol[detCountries[,1]], pch = 16, cex = 0.8, add=T)
  plot(st_geometry(COUNTRIES), add = TRUE)
}



## ------       6.4.3. EXTRACT GPS TRACKS LENGTHS ------

# ## COMBINE ALL TRACKS
# ALL_TRACKS <- rbind(
#   read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp")),
#   read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))) %>%
#   filter( Helikopter=="0", ## REMOVE HELICOPTER TRACKS
#           Jerv == "1")     ## KEEP ONLY WOLVERINE TRACKS
# 
# # ## REMOVE HELICOPTER TRACKS
# # ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helikopter=="0", ]
# #
# # ## KEEP ONLY WOLVERINE TRACKS
# # ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Jerv == "1", ]
# 
# ## SELECT TRACKS YEAR
# TRACKS_YEAR <- list()
# for(t in 1:n.years){
#   ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
#   TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr %in% YEARS[[t]][1] & ALL_TRACKS$Mth %in% sampling.months[[1]], ]
#   TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr %in% YEARS[[t]][2] & ALL_TRACKS$Mth %in% sampling.months[[2]], ]
#   tmpTRACKS <- rbind(TRACKS_1, TRACKS_2)
#   ## SIMPLIFY TRACKS SHAPES
#   # tmpTRACKS <- st_simplify(tmpTRACKS, T, 100)
#   tmpTRACKS <- st_intersection(tmpTRACKS, st_as_sf(myStudyArea))
#   ## NAME TRACKS
#   tmpTRACKS$ID <- 1:nrow(tmpTRACKS)
#   TRACKS_YEAR[[t]] <- tmpTRACKS
#   TRACKS_YEAR[[t]]$RovbaseID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
# }#t
# 
# TRACKS <- do.call(rbind,TRACKS_YEAR)
# 
# save( TRACKS, file = file.path(working.dir, "data/searchTracks.RData"))
# 
# rm(list = c("TRACKS_YEAR", "TRACKS_1", "TRACKS_2", "ALL_TRACKS", "tmpTRACKS"))

load(file = file.path(working.dir, "data/searchTracks.RData"))



##-- INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(detectors$main.detector.sp),
                                      rep(1,nrow(detectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE,
                             merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)


##-- CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
TRACKS.r <- list()
for(t in 1:n.years){
  TRACKSst <- TRACKS %>% filter(Yr == years[t])
  intersection <- st_intersection(detectorGrid, TRACKSst) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN)) ## Get total length searched in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[ ,t]
}#t


##-- PLOT CHECK 
if(plot.check){
  max <- max(unlist(lapply(TRACKS.r, function(x) max(x[], na.rm = T))))
  cuts <- seq(0,max,length.out = 100)   #set breaks
  col <- rev(terrain.colors(100))
  CountriesDetRes <- disaggregate(habitatRasters$Countries, fact = 2)
  CountriesDetRes <- crop(CountriesDetRes, TRACKS.r[[1]])
  rr <- TRACKS.r[[1]]
  rr[CountriesDetRes[]%in% 2] <- 1
  plot(rr)
  
  sum(st_length(TRACKS[TRACKS$Yr == years[t], ]))/1000
  sum(TRACKS.r[[t]][],na.rm=T)/1000
  
  
  pdf(file = file.path(working.dir, "figures/Tracks.pdf"))
  NORTRACKS <- SWETRACKS <- 0
  for(t in 1:n.years){
    plot( TRACKS.r[[t]], main = years[t], breaks = cuts, col = col, legend = FALSE)
    plot(st_geometry(habitat$habitat.poly), main = years[t], add = T)
    plot( TRACKS.r[[t]],
          legend.only = TRUE, breaks = cuts, col = col, legend.width = 2,
          axis.args = list(at = round(seq(0, max, length.out = 5), digits = 1),
                           labels = round(seq(0, max, length.out = 5), digits = 1),
                           cex.axis = 0.6),
          legend.args = list(text = '', side = 4, font = 2, line = 2.5, cex = 0.8))
    ##-- Summary tracks
    NORTRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 2],na.rm = T )/1000
    SWETRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 4],na.rm = T )/1000
  }#t
  years1 <- years + 1
  plot(SWETRACKS ~ years1, col = country.colors[2],
       lwd = 2, pch = 16, type = "b",
       ylim = c(0,300000), ylab = "sum tracks km")
  lines(NORTRACKS ~ years1, col = country.colors[1],
        lwd = 2, pch = 16, type = "b")
  legend("topright",c("N","S"), fill=country.colors)
  dev.off()
}



## ------       6.4.4. EXTRACT DISTANCES TO ROADS ------

##-- Load map of distance to roads (1km resolution)
#DistAllRoads <- raster(file.path(dir.dropbox,"DATA/GISData/Roads/MinDistAllRoads1km.tif"))
DistAllRoads <- raster::raster(file.path(data.dir,"GIS/Roads/MinDistAllRoads1km.tif"))

##-- Fasterize to remove values that fall in the sea
r <- fasterize::fasterize(sf::st_as_sf(GLOBALMAP), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r

rm(list = c("r"))

##-- AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = detectors$resolution/res(DistAllRoads),
                           fun = mean)

##-- EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)

##-- if NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( DistAllRoads,
                        detectors$main.detector.sp[isna, ],
                        buffer = 15000,
                        fun = mean,
                        na.rm = T)
detRoads[isna] <- tmp

##-- PLOT CHECK 
if(plot.check){
  par(mfrow = c(1,1))
  plot( st_geometry(GLOBALMAP),
        col = "gray80", main = "Distance to roads")
  plot( st_geometry(myStudyArea),
        col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot( st_geometry(COUNTRIES),
        col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot( DistAllRoads, add = T)
  plot( st_geometry(detectors$main.detector.sp),
        cex = DoScale(detRoads), pch = 16, add = T)
}



## ------       6.4.5. EXTRACT DAYS OF SNOW ------

#[PD] NEW SNOW FILE FROM ASUN!
#SNOW <- stack(paste0(dir.dropbox,"/DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2014_2025_Wolverine.tif"))

SNOW <- stack(file.path(data.dir,"GIS/AverageSnowCoverModisSeason2008_2024_Wolf.tif"))

##-- RENAME THE LAYERS
names(SNOW) <- paste(2008:2023, (2008:2023) + 1, sep = "_")

##-- SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years + 1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))

##-- EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = n.years)
det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:n.years] <- raster::extract(SNOW, det.sptransf)

##-- if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp

##-- PLOT CHECK
if(plot.check){
  plot( st_geometry(detectors$main.detector.sp),
        cex = DoScale(detSnow[ ,6], l = 0, u = 0.5),
        pch = 16)
}



## ------       6.4.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------         6.4.6.1. SKANDOBS ------

##-- Load the last SkandObs data file
skandObs <- readMostRecent( 
  path = data.dir,
  extension = ".xlsx",
  pattern = "Skandobs")

##-- Replace scandinavian characters
colnames(skandObs) <- translateForeignCharacters(data = colnames(skandObs))

skandObs <- skandObs %>%
  ##-- Extract important info (e.g. month, year)
  dplyr::mutate( date = as.POSIXct(strptime(date, "%Y-%m-%d")),
                 year = as.numeric(format(date,"%Y")),
                 month = as.numeric(format(date,"%m")),
                 species = stringi::stri_trans_general(species, "Latin-ASCII")) %>%
  ##-- Turn into spatial points object
  sf::st_as_sf(., coords = c("longitude","latitude")) %>%
  sf::st_set_crs(. , value = "EPSG:4326") %>%
  sf::st_transform(. ,sf::st_crs(COUNTIES))

##-- SUBSET BASED ON MONITORING SEASON 
subset <- skandObs$month %in% c(unlist(sampling.months))
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) #--- need to change for other species
skandObs <- skandObs[subset, ] 

##-- SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(habitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in% 1, ]
subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace, ]

##-- RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate( subdetectors.r,
                         fact = (detectors$resolution/detectors$resolution.sub))
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(skandObs[skandObs$monitoring.season %in% y, 1], r.detector , fun="count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})
r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))


##-- PLOT CHECK
if(plot.check){
  
  plot(st_geometry(habitat.rWthBufferPol))
  plot(st_geometry(skandObs), col = "red", add = T)
  
  ## SUMMARY SKANDOBS
  pdf(file = file.path(working.dir, "figures/skandObs.pdf"), width = 10)
  barplot(table(skandObs$monitoring.season ))
  barplot(table(skandObs$month ), xlab = "Months")
  barplot(table(skandObs$species))
  
  ## MAPS 
  par(mar = c(0,0,2,0))
  for(t in 1:n.years){
    plot( st_geometry(myStudyArea), main = years[t])
    plot( st_geometry(skandObs[skandObs$monitoring.season %in% years[t], ]),
          pch = 16, col = "red", cex = 0.1, add = T)
  }
  dev.off()
}



## ------         6.4.6.2. ROVBASE ------

##-- Load the last Rovbase data files
#rovbaseObs1 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB2810202415264376.xlsx"))
rovbaseObs1 <- readMostRecent( 
  path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
  extension = ".xlsx",
  pattern = "RIB2810202415264376")
#rovbaseObs2 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152348493.xlsx"))
rovbaseObs2 <- readMostRecent( 
  path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
  extension = ".xlsx",
  pattern = "RIB28102024152348493")
#rovbaseObs3 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152447860.xlsx"))
rovbaseObs3 <- readMostRecent( 
  path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
  extension = ".xlsx",
  pattern = "RIB28102024152447860")
#rovbaseObs4 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152538742.xlsx"))
rovbaseObs4 <- readMostRecent( 
  path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
  extension = ".xlsx",
  pattern = "RIB28102024152538742")
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3,rovbaseObs4)
colnames(rovbaseObs) <- translateForeignCharacters(dat = colnames(rovbaseObs))
rovbaseObs$Proevetype <- translateForeignCharacters(dat = rovbaseObs$Proevetype)

##-- Remove un-necessary objects
rm(list = c("rovbaseObs1","rovbaseObs2","rovbaseObs3","rovbaseObs4"))


##-- GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`), ]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

##-- DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)

# ##-- SUBSET THE DATA 
# filter <- list(
#   species = "Jerv",
#   type = c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
#      "Saliv/Spytt", "Loepeblod", "Vev"),
#   month = unlist(sampling.months))
# 
# ##-- SUBSET MONTH AND TYPE OF SAMPLE
# subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
# rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
# rovbaseObs.sp <- rovbaseObs.sp[subset, ] 
# 
# ##-- SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
# subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$`Art (Proeve)`) 
# rovbaseObs.sp <- rovbaseObs.sp[!subset, ] 
# 
# ##-- SUBSET BASED ON SPACE 
# subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
# rovbaseObs.sp <- rovbaseObs.sp[subsetSpace,] 

##-- GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`), ]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

##-- DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)


## New version [PD]
rovbaseObs.sp2 <- rovbaseObs.sp %>%
  ##-- EXTRACT YEAR AND MONTH
  mutate(
    year = as.numeric(format(Funnetdato,"%Y")),
    month = as.numeric(format(Funnetdato,"%m"))) %>%
  
  filter( 
    ##-- SUBSET MONTH 
    month %in% unlist(sampling.months),
    ##-- SUBSET TYPE OF SAMPLE
    Proevetype %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)",
                       "Sekret (Jerv)","Saliv/Spytt", "Loepeblod", "Vev"),
    ##-- SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
    !(`Art (Analyse)` %in% "Jerv" & !is.na(Individ)),
    ##-- SUBSET BASED ON SPACE 
    !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))) %>%
  mutate(monitoring.season = ifelse(month < 12, year, year+1))


##-- RASTERIZE 
r.detector <- aggregate( subdetectors.r,
                         fact = (detectors$resolution/detectors$resolution.sub))
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(rovbaseObs.sp[rovbaseObs.sp$monitoring.season %in% y, 1], r.detector , fun="count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})

r.OtherSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.OtherSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))

##-- PLOT CHECK
if(plot.check){
  pdf(file = file.path(working.dir, "figures/mapStructuredOthers.pdf"))
  for(t in 1:n.years){
    tmpOthers <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t] &
                                           !myFilteredData.sp$alive$structured, ]
    tmpStruct <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t] &
                                           myFilteredData.sp$alive$structured, ]
    
    par(mfrow=c(2,2),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]], main=paste(years[t],"\n Rovbase Samples Structured"), box=F, axes=F)
    plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    plot(r.OtherSamplesBinary[[t]],main=paste(years[t],"\n Rovbase Samples Opportunistic"), box=F, axes=F)
    plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
    
    plot(r.skandObsSamplesBinary[[t]], main=paste(years[t]), box=F, axes=F)
    plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    plot(r.skandObsSamplesBinary[[t]],main=paste(years[t],"\n SkandObs Opportunistic"), box=F, axes=F)
    plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
  }
  dev.off()
}



## ------         6.4.6.3. COMBINE ROVBASE & SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:n.years){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <- 1
}

##-- PLOT CHECK
if(plot.check){
  for(t in 1:n.years){
    par(mfrow=c(1,3),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]],main=years[t])
    plot(r.skandObsSamplesBinary[[t]])
    plot(r.SkandObsOtherSamplesBinary[[t]])
  }#t  
}



## ------         6.4.6.4. SMOOTH THE BINARY MAP ------

## we tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
habOwin <- spatstat.geom::as.owin(as.vector(extent(r.detector)))
cutoff <- 1
ds.list <- lapply(years,function(y){
  ## ROVBASE DATA 
  pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y,]
  ## SKANDOBS
  pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y,] )
  ## SMOOTH AND RASTERIZE
  p <-  spatstat.geom::ppp(pts[,1], pts[,2], window = habOwin)
  ds <- density(p, adjust=0.02) #---change bandwith (smoothing) with "adjust
  ds <- raster(ds)
  
  ds <- ds1 <- raster::resample(ds, r.detector) #mask(ds,rasterToPolygons(myHabitat.list$habitat.rWthBuffer,function(x) x==1))
  threshold <- 0.1 / prod(res(ds)) #--number per 1 unit of the projected raster (meters)
  ds1[] <- ifelse(ds[]<threshold,0,1)
  ds1 <- mask(ds1, habitat.rWthBufferPol)
  ds <- mask(ds, habitat.rWthBufferPol)
  
  return(list(ds,ds1))
})

ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))
names(ds.brick) <- names(ds.brickCont) <-years

##-- PLOT CHECK
if(plot.check){
  par(mfrow = c(1,3))
  plot(r.SkandObsOtherSamplesBinary[[t]], main = "Raw Binary", axes = F, box = F)
  plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
  plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
}



## ------         6.4.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------
# 
# ## [PD] : USELESS !!!
# 
# ## IDENTIFY HAIR SAMPLES
# tmpHair <- myFilteredData.sp$alive %>% filter(hairTrap)
#   
# ## MANUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL. 
# tmpyr <- unique(tmpHair$Year)
# for( i in 1:length(tmpyr)){
#   t <- which(years %in% tmpyr)
#   whereHair <- raster::extract(r.SkandObsOtherSamplesBinary[[t]],tmpHair,cellnumbers=T)
#   r.SkandObsOtherSamplesBinary[[t]][whereHair[,1]] <- 1
#   plot(r.SkandObsOtherSamplesBinary[[t]])
#   plot(tmpHair$geometry,add=T,col="red")
# }



## ------         6.4.6.6. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
detOtherSamples[ ,1:n.years] <- raster::extract(r.SkandObsOtherSamplesBinary, detectors$main.detector.sp)
colSums(detOtherSamples)



## ------       6.4.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],2))
detCovs[,,1] <- detTracks
detCovs[,,2] <- detSnow

detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],3))
detCovsOth[,,1] <- detSnow
detCovsOth[,,2] <- matrix(detRoads,length(detRoads),n.years)
detCovsOth[,,3] <- detOtherSamples


## CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}


##-- PLOT CHECK
if(plot.check){
  tmp <- detectorGrid.r
  par(mfrow=c(2,5),mar=c(0,0,0,0))
  max <- max(detCovsOth[,,2])
  cuts <- seq(0,max,length.out = 100) 
  col <- rev(terrain.colors(100))
  for(t in 1:n.years){
    plot(detectorGrid.r, col=c(grey(0.2),grey(0.8)),axes=F,legend=F,box=F,)
    tmp[!is.na( detectorGrid.r)] <- detCovsOth[,t,2]
    plot(tmp,axes=F,legend=F,box=F,breaks = cuts, col=col,add=T)
  }
  
  dev.off()
  pdf(file = file.path(working.dir, "figures/detections over space and time.pdf"))
  for(t in 1:n.years){
    ## NGS DETECTIONS TOTAL
    tempTotal <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]
    NGS_TabTotal <- table(tempTotal$Country)
    ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
    ## ALIVE DETECTIONS INSIDE STUDY AREA/SAMPLING PERIOD
    tempIn <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]
    NGS_TabIn <- table(tempIn$Country)
    ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
    ## PLOT NGS SAMPLES
    plot(st_geometry(GLOBALMAP), col="gray80")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
    plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
    # points(tempTotal, pch = 21, bg = "darkred")
    plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
    ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    graphics::text(x = 100000, y = 7200000, labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="N"],"NGS"), cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 100000, y = 7270000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 820000, y = 6780000, labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="S"],"NGS"), cex = 1.1, col = "navyblue", font = 2)
    graphics::text(x = 820000, y = 6850000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
    ## ADD OVERALL NUMBERS
    mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
    mtext(text = paste(sum(NGS_TabIn), "NGS/", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
    mtext(text = paste(sum(NGS_TabTotal)-sum(NGS_TabIn), "NGS/", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
  }#t
  dev.off()
}



## ------   7. RESCALE COORDINATES ------

# ##-- Rescale detector coordinates to the habitat 
# scaledCoords <-  scaleCoordsToHabitatGrid(coordsData = detector.xy,
#                                           coordsHabitatGridCenter = habitat.xy)
# ##-- Scaled habitat windows boundaries
# lowerHabCoords <- scaledCoords$coordsHabitatGridCenterScaled - 0.5
# upperHabCoords <- scaledCoords$coordsHabitatGridCenterScaled + 0.5
# ##-- Scaled detectors coordinates
# scaledDetCoords <- scaledCoords$coordsDataScaled

##-- Rescale coordinates
scaledCoords <- nimbleSCR::scaleCoordsToHabitatGrid(
  coordsData = detectors$detectors.df[ ,c("x","y")],
  coordsHabitatGridCenter = habitat$habitat.df[ ,c("x","y")])

##-- Scaled habitat window coordinates
habitat$scaledCoords <- scaledCoords$coordsHabitatGridCenterScaled
habitat$scaledLowerCoords <- habitat$scaledCoords - 0.5
habitat$scaledUpperCoords <- habitat$scaledCoords + 0.5

##-- Scaled detector coordinates
detectors$scaledCoords <- scaledCoords$coordsDataScaled



## ------   8. CREATE LOCAL OBJECTS -----

## [CM] reduce multiplicator to 3 
maxDistReCalc <- 2.1*detectors$maxDist #+ sqrt(2*(DETECTIONS$resizeFactor*HABITAT$habResolution)^2)

# DetectorIndexLESS <- GetDetectorIndexLESS(
#   habitat.mx = habitat$habitat.mx,
#   detectors.xy = nimData$detector.xy,
#   maxDist = maxDistReCalc/res(habitat$habitat.r)[1],
#   ResizeFactor = 1,
#   plot.check = TRUE)

DetectorIndexLESS <- getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = scaledDetCoords,
  dmax = maxDistReCalc/res(habitat$habitat.r)[1],
  resizeFactor = 1,
  plot.check = TRUE)

##-- Get local detectors
detectors$localObjects <- getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = detectors$scaledCoords,
  dmax = detectors$maxDist/habitat$resolution,
  resizeFactor = detectors$resize.factor,
  plot.check = F)



## -------   5. SAVE STATE-SPACE CHARACTERISTICS -----

save( habitat,
      file = file.path( working.dir, "data",
                        paste0("Habitat_bear_", DATE, ".RData")))

save( detectors,
      file = file.path( working.dir,"data",
                        paste0("Detectors_bear_", DATE, ".RData")))

## -------   9. GENERATE y DETECTION ARRAYS ------

## ------     9.1. ASSIGN SAMPLES TO DETECTORS -----

##-- ALL SAMPLES
myData.alive <- assignDetectors( 
  data = myFilteredData.sp$alive,                
  detectors = detectors$main.detector.sp,
  subDetectors = detectors$detector.sp,
  radius = detectors$resolution)
# ## STRUCTURED
# myData.aliveStruc <- AssignDetectors_v3sf( 
#   myData = myFilteredData.spStructured,                
#   myDetectors = detectors$main.detector.sp,
#   mysubDetectors = detectors$detector.sp,
#   radius = detectors$resolution)
# ## OTHERS
# myData.aliveOthers <- AssignDetectors_v3sf( 
#   myData = myFilteredData.spOthers,                
#   myDetectors = detectors$main.detector.sp,
#   mysubDetectors = detectors$detector.sp,
#   radius = detectors$resolution)
##-- DEAD RECOVERY
myData.dead <- assignDetectors( 
  data = myFilteredData.sp$dead.recovery,
  detectors = detectors$main.detector.sp,
  radius = detectors$resolution)



## ------     9.2. FIX DETECTIONS IN NORRBOTTEN ------

### MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET 
### ASSIGNED TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING.
### FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE
### OF NORRBOTTEN

##-- Identify sub-detectors inside Norrbotten
subDetsNorrbotten <- which( detectors$detector.sp$main.cell.id %in% 
                              detectors$main.detector.sp$main.cell.id[detsNorrbotten])

##-- Identify which detections are assigned to Norrbotten in years not sampled
whichDets <- which(!myData.alive$data.sp$Year %in% yearsSampledNorrb &
                     myData.alive$data.sp$Detector %in% detsNorrbotten)
length(whichDets)

##-- Loop over flagged detections and assign them to closest detector outside Norrbotten
for(i in 1:length(whichDets)){
  tmp <- myData.alive$data.sp[whichDets[i], ]
  ## Calculate distance to all sub-detectors
  dist <- st_distance(tmp, detectors$detector.sp)
  ## Artificially increase distance for detectors in Norrbotten 
  dist[ ,subDetsNorrbotten] <- 500000
  ## Assign detection to closest sub-detector outside Norrbotten
  myData.alive$data.sp$sub.detector[whichDets[i]] <- which.min(dist[1, ])
  ## Assign detection to the corresponding main detector outside Norrbotten
  thisDet <- detectors$detector.sp$main.cell.id[which.min(dist[1, ])]
  myData.alive$data.sp$Detector[whichDets[i]] <- which(detectors$main.detector.sp$main.cell.id == thisDet)
}#i


##-- SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
sum(myData.alive$data.sp$sub.detector[!myData.alive$data.sp$Year %in% yearsSampledNorrb] %in% subDetsNorrbotten)
sum(myData.alive$data.sp$Detector[!myData.alive$data.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)


# ##-- MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET ASSIGNED 
# ##-- TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING
# ##-- FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE OF NORRBOTTEN
# sum(myData.alive$data.sp$Detector[!myData.alive$data.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)
# whichdets <- which(!myData.alive$data.sp$Year %in% yearsSampledNorrb &
#                      myData.alive$data.sp$Detector %in% detsNorrbotten)
# # whichdetsStruc <- which(!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb &
# #                           myData.aliveStruc$myData.sp$Detector %in% detsNorrbotten)
# # whichdetsOther <- which(!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb &
# #                           myData.aliveOthers$myData.sp$Detector %in% detsNorrbotten)
# 
# # ##-- ASSIGN DETECTORS 
# # subdetector.sf <- detectors$detector.sp
# 
# ## ALL
# ## [PD] N
# for(i in 1:length(whichdets)){
#   tmp <- myData.alive$data.sp[whichdets[i], ]
#   ##-- CALCULATE DISTANCE TO MAIN DETECTOR 
#   dist <- st_distance(tmp, detectors$main.detector.sp)
#   ##-- Artificially increase distance for detectors in Norrbotten 
#   dist[ ,detsNorrbotten] <- 500000
#   idMain <- which.min(dist[1, ])
#   myData.alive$data.sp$Detector[whichdets[i]] <- idMain 
#   ##-- SUBDETECTOR
#   dist <- st_distance(st_as_sf(tmp), detectors$detector.sp)
#   ##-- Artificially increase distance for sub-detectors in Norrbotten
#   dist[,detsNorrbotten] <- 500000
#   idSub <- which.min(dist[1,])
#   myData.alive$data.sp$sub.detector[whichdets[i]] <- idSub 
# }
# #STRUCTURED
# for(i in 1:length(whichdetsStruc)){
#   tmp <- myData.aliveStruc$myData.sp[whichdetsStruc[i],]
#   # MAIN DETECTOR 
#   dist <- st_distance(tmp, detectors$main.detector.sp)
#   #Artificially increase distance for detectors in noorrbotten 
#   dist[,detsNorrbotten] <- 500000
#   idMain <- which.min(dist[1,])
#   myData.aliveStruc$myData.sp$Detector[whichdetsStruc[i]] <- idMain 
#   # SUBDETECTOR
#   dist <- st_distance(tmp, subdetector.sf )
#   #Artificially increase distance for detectors in noorrbotten
#   dist[,detsNorrbotten] <- 500000
#   idSub <- which.min(dist[1,])
#   myData.aliveStruc$myData.sp$sub.detector[whichdetsStruc[i]] <- idSub 
# }
# #OTHER
# for(i in 1:length(whichdetsOther)){
#   tmp <- myData.aliveOthers$myData.sp[whichdetsOther[i],]
#   # MAIN DETECTOR 
#   dist <- st_distance(tmp, detectors$main.detector.sp)
#   #Artificially increase distance for detectors in noorrbotten 
#   dist[,detsNorrbotten] <- 500000
#   idMain <- which.min(dist[1,])
#   myData.aliveOthers$myData.sp$Detector[whichdetsOther[i]] <- idMain 
#   # SUBDETECTOR
#   dist <- st_distance(tmp, subdetector.sf )
#   #Artificially increase distance for detectors in noorrbotten
#   dist[,detsNorrbotten] <- 500000
#   idSub <- which.min(dist[1,])
#   myData.aliveOthers$myData.sp$sub.detector[whichdetsOther[i]] <- idSub 
# }
# 
# ## SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
# sum(myData.alive$data.sp$Detector[!myData.alive$data.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
# # sum(myData.aliveOthers$myData.sp$Detector[!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
# # sum(myData.aliveStruc$myData.sp$Detector[!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)


## [PD] : USELESS!!!
# ## EXPORT THE DATA 
# if(DATA$sex=="Hann"){
#   assign("myFilteredData.spM", myFilteredData.sp)
#   assign("myFilteredData.spOthersM", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
#   save(myFilteredData.spM, myFullData.spM,
#        myFilteredData.spOthersM,myFilteredData.spStructuredM,
#        file = file.path(working.dir, "data/NGSData.RData")) 
#   
# } else {
#   assign("myFilteredData.spF", myFilteredData.sp)
#   assign("myFilteredData.spOthersF", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
#   save(myFilteredData.spF, myFullData.spF,
#        myFilteredData.spOthersF,myFilteredData.spStructuredF,
#        file = file.path(working.dir, "data/NGSData.RData")) 
# }



## ------     9.3. GENERATE DETECTION HISTORIES : y.alive[i,j,t] & y.dead[i,t] ------

##-- ALL SAMPLES
y.ar <- makeY( data = myData.alive$data.sp,
               detectors = detectors$main.detector.sp,
               method = "Binomial",
               data2 = myData.dead,
               detectors2 = detectors$main.detector.sp,
               returnIdvector = TRUE)

##-- STRUCTURED
y.arStruc <- makeY( data = myData.alive$data.sp[myData.alive$data.sp$structured, ],
                    detectors = detectors$main.detector.sp,
                    method = "Binomial",
                    data2 = myData.dead,
                    detectors2 = detectors$main.detector.sp,
                    returnIdvector = TRUE)

##-- OTHERS
y.arOth <- makeY( data = myData.alive$data.sp[!myData.alive$data.sp$structured, ],
                  detectors = detectors$main.detector.sp,
                  method = "Binomial",
                  data2 = myData.dead,
                  detectors2 = detectors$main.detector.sp,
                  returnIdvector = TRUE)

##-- MAKE SURE THE Y HAVE THE SAME DIMENSIONS
y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- array( 0, 
                                                   dim = dim(y.ar$y.ar),
                                                   dimnames = dimnames(y.ar$y.ar))

##-- FILL IN THE Y ARRAYS 
y.ar.ALIVEOthers[dimnames(y.arOth$y.ar)[[1]], , ] <- y.arOth$y.ar
y.ar.ALIVEStructured[dimnames(y.arStruc$y.ar)[[1]], , ] <- y.arStruc$y.ar

# all(dimnames(y.arOth$y.ar)[[1]] %in% dimnames(y.ar$y.ar)[[1]])
# sum(y.ar.ALIVEOthers[1307,,8])
# sum(y.ar.ALIVEStructured[1307,,8])


##-- PROJECT THE DEATH TO THE NEXT OCCASION.
y.ar.DEADProjected <- y.ar$y.ar2 
y.ar.DEADProjected[] <- 0
for(t in 2:n.years){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}

##-- CREATE BINARY DEAD RECOVERY HISTORIES (0: not recovered ; 1: recovered dead)
y.ar.DEAD <- apply(y.ar.DEADProjected, c(1,3), function(x){as.numeric(sum(x)>0)})
dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])

# dim(y.ar.DEAD)
# y.ar.DEAD[y.ar.DEAD>0] <- 1
# dim(y.ar$y.ar2)



## ------     9.4. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------

distances <- list()
for(t in 1:n.years){
  
  print(paste("------ ", t ," -------", sep = "" ))
  
  ##-- CALCULATE DISTANCES BETWEEN PAIRS OF DETECTIONS
  distances[[t]] <- CheckDistanceDetections( 
    y = y.ar$y.ar[,,t], 
    detector.xy = detector.xy, 
    max.distance = detectors$maxDist,
    method = "pairwise",
    plot.check = F)
  
  ##-- PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
  if(plot.check){
    
    par(mfrow = c(1,1))
    
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      count <- 1
      for(i in affected.ids){
        ##-- PLOT STUDY AREA & DETECTORS
        plot( st_geometry(myStudyArea),
              main = paste0("t: ",t,"     i: ", names(affected.ids)[count]))
        scalebar( 2 * detectors$maxDist,
                  xy = c(800000,6700000),
                  type = "bar", divs = 2, below = "km",
                  label = c(0, detectors$maxDist/1000, detectors$maxDist/500),
                  cex = 0.8, adj = c(0.5,-0.9))
        plot(st_geometry(COUNTRIES), add = T)
        plot( st_geometry(detectors$main.detector.sp),
              add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        ##-- PLOT LOCATIONS OF DNA SAMPLES
        tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar$y.ar)[[1]][i] &
                                         myFilteredData.sp$alive$Year == years[t], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
        arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
               length = 0.1, lwd = 1)
        
        ##-- PLOT DETECTORS LOCATIONS
        plot( st_geometry(detectors$main.detector.sp[which(y.ar$y.ar[i,,t] > 0), ]),
              pch = 16, col = "red",add=T)
        
        tmp2 <- detectors$main.detector.sp[which(y.ar$y.ar[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
        plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
        
        ##
        count <- count + 1
      }#i
    }#if
  }#if plot.check
  
  ##-- REMOVE DETECTIONS THAT ARE FURTHER THAN  THE THRESHOLD
  y.ar$y.ar[,,t] <- y.ar$y.ar[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  
  ##-- REMOVE DETECTIONS ALSO IN MYDATA TO RUN GETSINITS
  # tmpmyData.sp <- myData.alive$data.sp
  # 
  # for(t in 1:n.years){
  #   affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
    #       # tmp$myData.sp <- tmp$myData.sp[(tmp$myData.sp$Id %in% idd[i] & 
    #       #                           tmp$myData.sp$Detector %in% detIds) &  ,] 
    #       # tmp$myData.sp[tmp$myData.sp$Id %in% idd[i] ,] 
    #       tmpmyData.sp[(tmpmyData.sp$Id %in% idd[i] &
    #                              tmpmyData.sp$Detector %in% detIds &
    #                              tmpmyData.sp$Year %in% years[t]),]
    #       tmpmyData.sp[(tmpmyData.sp$Id %in% idd[i] &
    #                       tmpmyData.sp$Year %in% years[t]),]
    #       
    #       
    #       tmpmyData.sp <- tmpmyData.sp[!(tmpmyData.sp$Id %in% idd[i] &
    #                                                          tmpmyData.sp$Detector %in% detIds &
    #                                                          tmpmyData.sp$Year %in% years[t]),]
    #      
    #       
    #       if(sum(which(!unique(myData.alive$data.sp$Id) %in% unique(tmpmyData.sp$Id))>0)){
    #         print(idd[i])
    #       }
    #   }
    #   }
    
    myData.alive$data.sp <- myData.alive$data.sp[!(myData.alive$data.sp$Id %in% idd[i] &
                                                     myData.alive$data.sp$Detector %in% detIds &
                                                     myData.alive$data.sp$Year %in% years[t]),]
  }
  # tmp <- myData.alive$data.sp[(myData.alive$data.sp$Id %in% idd[i] &
  #                                  myData.alive$data.sp$Year %in% years[t]),]
  # tmp[(tmp$Id %in% idd[i]),]
  # myData.alive$data.sp[(myData.alive$data.sp$Id %in% idd[i]&
  # #                           myData.alive$data.sp$Detector %in% detIds &  myData.alive$data.sp$Year %in% years[t]),]
  #     if(length(which(myData.alive$data.sp$Id %in% idd[i]))==0){
  #       print(paste(idd[i], t, sep="_"))
  #     }
}



## ------     9.5. GENERATE INDIVIDUAL-LEVEL COVARIATES ------

## ------       9.5.1. TRAP-RESPONSE ------

##-- Make matrix of previous capture indicator
already.detected <- makeTrapResponseCov(
  data = myFullData.sp$alive,
  data.dead = myFullData.sp$dead.recovery)

##-- Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar$y.ar)[[3]]]

##-- Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar$y.ar)[[1]], ]

##-- Plot an image of the matrix
if(plot.check){
  par(mfrow = c(1,1))
  barplot(colSums(apply(y.ar$y.ar, c(1,3), function(x) any(x>0))))
  barplot(colSums(already.detected), add = TRUE, col = "gray40")
  legend( x = 0, y = 250, 
          legend = c("newly Det", "already Det"),
          fill = c("gray80", "gray40"))
}



## ------       9.5.2. AGE ------

min.age <- age <- precapture <- matrix( NA, dim(y.ar$y.ar)[1], dim(y.ar$y.ar)[3],
                                        dimnames = list(y.ar$Id.vector, years))

temp <- apply(y.ar$y.ar, c(1,3), sum)
year.first.capture <- apply(temp, 1, function(x)min(years[which(x>0)]))
year.first.capture[is.infinite(year.first.capture)] <- NA
names(year.first.capture) <- y.ar$Id.vector

for(i in y.ar$Id.vector){
  this.set <- myData.dead[myData.dead$Id == i, ]
  year.dead <- myData.dead$Death[myData.dead$Id == i]
  year.first.captured <- year.first.capture[i]
  precapture[i,] <- as.numeric(years < year.first.captured)
  if(all(is.na(precapture[i,])))precapture[i,] <- 1
  latest.recruitment.year <- min(year.dead,year.first.captured, na.rm = TRUE) 
  
  try({
    min.age[i,] <- years-latest.recruitment.year
  },silent = TRUE)
  
  try({
    birth.year <- this.set$Death-this.set$min.age
    if(birth.year < latest.recruitment.year) min.age[i,] <- years-birth.year 
  },silent = TRUE)
  
  try({
    birth.year <- this.set$Death - this.set$age
    age[i,] <- years-birth.year
  }, silent = TRUE)
}

image(t(min.age))
image(t(age))







## ------   10. MAKE AUGMENTATION ------

##-- DATA ARRAYS
y.alive <- makeAugmentation( y = y.ar$y.ar,
                             aug.factor = aug.factor,
                             replace.value = 0)

y.dead <- makeAugmentation( y = y.ar.DEAD,
                            aug.factor = aug.factor,
                            replace.value = 0)

y.aliveOthers <- makeAugmentation( y = y.ar.ALIVEOthers,
                                   aug.factor = aug.factor,
                                   replace.value = 0)

y.aliveStructured <- makeAugmentation( y = y.ar.ALIVEStructured,
                                       aug.factor = aug.factor, 
                                       replace.value = 0)

##-- INDIVIDUAL COVARIATES
already.detected <- makeAugmentation( y = already.detected,
                                      aug.factor = aug.factor,
                                      replace.value = 0)

# age <- makeAugmentation( y = age,
#                          aug.factor = aug.factor, 
#                          replace.value = NA)
# min.age <- makeAugmentation( y = min.age,
#                              aug.factor = aug.factor,
#                              replace.value = NA)
# precapture <- makeAugmentation( y = precapture,
#                                 aug.factor = aug.factor,
#                                 replace.value = 0)



## ------   11. TRANSFORM Y TO SPARSE MATRICES ------

##-- STRUCTURED
#SparseY <- GetSparseY(y.aliveStructured)
y.sparse <- nimbleSCR::getSparseY(y.aliveStructured)

##-- OTHER
#SparseYOth <- GetSparseY(y.aliveOthers)
y.sparseOth <- nimbleSCR::getSparseY(y.aliveOthers)



## ------ III. MODEL SETTING & RUNNING ------- 

## ------   1. NIMBLE MODEL DEFINITION ------

modelCode <- nimbleCode({
  
  ##------ SPATIAL PROCESS ------## 
  
  dmean ~ dunif(0,100)
  lambda <- 1/dmean
  
  betaDens ~ dnorm(0.0,0.01)
  habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i,1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
      upperCoords = upperHabCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    for(t in 2:n.years){
      sxy[i,1:2,t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
        upperCoords = upperHabCoords[1:numHabWindows,1:2],
        s = sxy[i,1:2,t-1],
        lambda = lambda,
        baseIntensities = habIntensity[1:numHabWindows],
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max,
        numWindows = numHabWindows)
    }#i  
  }#t
  
  
  
  ##----- DEMOGRAPHIC PROCESS -----## 
  
  omeg1[1:2] ~ ddirch(alpha[1:2])   
  
  for(t in 1:n.years1){
    # PRIORS 
    gamma[t] ~ dunif(0,1)
    phi[t] ~ dunif(0,1)
    
    # TRANSITION MATRIX
    omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
    omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
    omega[3,1:3,t] <- c(0,0,1)
  }#t
  
  
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:2]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
    }#i 								
  }#t 
  
  
  
  ##----- DETECTION PROCESS -----## 
  
  for(t in 1:n.years){
    
    sigma[t] ~ dunif(0,4)
    
    ## Systematic sampling
    betaResponse[t] ~ dunif(-5,5)
    
    for(c in 1:n.covs){
      betaCovs[c,t] ~ dunif(-5,5)
    }#c 
    
    for(c in 1:n.counties){
      p01[c,t] ~ dunif(0,1)
      p0[c,t] <- p01[c,t] * countyToggle[c,t]## toggle counties
    }#c  
    
    
    ## Opportunistic sampling
    betaResponseOth[t] ~ dunif(-5,5)
    
    for(c in 1:n.covsOth){
      betaCovsOth[c,t] ~ dunif(-5,5)
    }#c 
    
    for(c in 1:n.countries){
      p01Oth[c,t] ~ dunif(0,1)
      p0Oth[c,t] <- p01Oth[c,t] * countryToggle[c,t]## toggle counties
    }#c  
  }#t
  
  
  ## Individual response
  pResponse ~ dunif(0, 1)
  
  for(i in 1:n.individuals){ 
    detResponse[i,1] ~ dbern(pResponse)
  }#i
  
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){
      
      y[i,1:maxDetNums,t] ~ dbinomLocal_normalCovsResponse( 
        detNums = detNums[i,t],
        detIndices = detIndices[i,1:maxDetNums,t],
        size = size[1:n.detectors],
        p0State = p0[1:n.counties,t],
        sigma = sigma[t],
        s = sxy[i,1:2,t],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = resizeFactor,
        lengthYCombined = lengthYCombined,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        trapCountries = detCounties[1:n.detectors],
        trapCovs = detCovs[1:n.detectors,t,1:n.covs],
        trapBetas = betaCovs[1:n.covs,t],
        responseCovs = detResponse[i,t],
        responseBetas = betaResponse[t])
      
      y.Oth[i,1:maxDetNumsOth,t] ~ dbinomLocal_normalCovsResponse( 
        detNums = detNumsOth[i,t],
        detIndices = detIndicesOth[i,1:maxDetNumsOth,t],
        size = size[1:n.detectors],
        p0State = p0Oth[1:n.counties,t],
        sigma = sigma[t],
        s = sxy[i,1:2,t],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = resizeFactor,
        lengthYCombined = lengthYCombined.Oth,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        trapCountries = detCountries[1:n.detectors],
        trapCovs = detCovsOth[1:n.detectors,t,1:n.covsOth],
        trapBetas = betaCovsOth[1:n.covs,t],
        responseCovs = detResponse[i,t],
        responseBetas = betaResponseOth[t])
    }#i
  }#t
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){
      y.dead.legal[i,t] ~ dbern(z[i,t] == 3) 
      y.dead.other[i,t] ~ dbern(z[i,t] == 4) 
    }#i
  }#t
  
  
  ##---------- DERIVED PARAMETERS ----------##
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){ 
      isAlive[i,t] <- (z[i,t] == 2) 
    }#i
    N[t] <- sum(isAlive[1:n.individuals,t])
  }#t
  
})



## ------   2. NIMBLE CONSTANTS ------

nimConstants <- list( 
  n.individuals = dim(y.alive)[1],
  n.detectors = dim(y.alive)[2],
  numHabWindows = nrow(lowerHabCoords)[1],
  n.years = dim(y.alive)[3], 
  n.years1 = dim(y.alive)[3]-1, 
  n.covs = dim(detCovs)[3],
  n.covs.Oth = dim(detCovsOth)[3],
  n.countries = max(detCountries)+1,
  n.counties = max(detCounties),
  y.max = dim(habIDCells.mx)[1],
  x.max = dim(habIDCells.mx)[2],
  countyToggle = countyToggle,
  countryToggle = countryToggle,
  resizeFactor = DetectorIndexLESS$resizeFactor,
  y.maxDet = dim(DetectorIndexLESS$habitatGrid)[1],
  x.maxDet = dim(DetectorIndexLESS$habitatGrid)[2],
  n.cellsSparse = dim(DetectorIndexLESS$detectorIndex)[1],
  maxNBDets = DetectorIndexLESS$maxNBDets,
  maxDetNums = y.sparse$maxDetNums,
  maxDetNums.Oth = y.sparseOth$maxDetNums)



## ------   3. NIMBLE INITS ------

## ------     3.1. GENERATE INITIAL z ------

z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, NA)

z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    zz[range.det[1]:range.det[2]] <- 2
  }
  return(zz)
}))

z.init <- t(apply(z, 1, function(zz){
  out <- zz
  out[] <- 1
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    if(range.det[1]>1)zz[1:(range.det[1]-1)] <- 1
    if(range.det[2]<length(zz))zz[(range.det[2]+1):length(zz)] <- 3
    out[] <- zz
  } 
  return(out)
}))

z.init <- ifelse(!is.na(z), NA, z.init)



## ------     3.2. LATENT VARIABLE DET RESPONSE ------

detResponse <- already.detected 
detResponse[rownames(detResponse) %in% "Augmented", 1]  <- NA
InitsDetResponse <- detResponse
InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1,0.5)
InitsDetResponse[!is.na(detResponse)] <- NA



## ------   4. NIMBLE DATA ------

nimData <- list( 
  z = z,   
  y = y.sparse$y,
  detIndices = y.sparse$detIndices,
  detNums = y.sparse$detNums,
  y.Oth = y.sparseOth$y, 
  detIndicesOth = y.sparseOth$detIndices,
  detNumsOth = y.sparseOth$detNums,
  lowerHabCoords = lowerHabCoords, 
  upperHabCoords = upperHabCoords, 
  detCounties = detCounties,
  detCountries = detCountries,
  detCovs = detCovs,
  detCovsOth = detCovsOth,
  detResponse = detResponse,
  denCounts = denCounts,
  detectorIndex = DetectorIndexLESS$localIndices,
  nDetectorsLESS = DetectorIndexLESS$numLocalIndices,
  habitatIDDet = DetectorIndexLESS$habitatGrid,
  size = n.trials,
  alpha = rep(1,2),
  detector.xy = scaledDetCoords,
  #sxy = sxy.data,
  habitatGrid = habIDCells.mx)



## ------   5. NIMBLE PARAMETERS ------

nimParams <- c( "N", "lambda", "dmean", "betaDens",
                "omeg1", "gamma", "phi",
                "pResponse", "sigma",
                "p0", "betaResponse", "betaCovs",
                "p0Oth", "betaResponseOth", "betaCovsOth")

nimParams2 <- c("z", "sxy")



## ------   6. CONVERT TO CACHED DETECTORS & SPARSE MATRIX ------
#
## ------     6.1. RESCALE COORDINATES ------
# 
# scaledCoords <-  scaleCoordsToHabitatGrid(
#   coordsData = detector.xy,
#   coordsHabitatGridCenter = habitat$habitat.xy)
# 
# lowerHabCoords <- scaledCoords$coordsHabitatGridCenterScaled - 0.5
# upperHabCoords <- scaledCoords$coordsHabitatGridCenterScaled + 0.5
# 
# scaledDetCoords <- scaledCoords$coordsDataScaled
# 
# # # HABITAT
# # ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData = lowerHabCoords,
# #                                               coordsHabitatGridCenter = habitat$habitat.xy,
# #                                               scaleToGrid =T )$coordsDataScaled
# # ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData = upperHabCoords,
# #                                               coordsHabitatGridCenter = habitat$habitat.xy,
# #                                               scaleToGrid =T )$coordsDataScaled
# # 
# # ScaledUpperCoords[ ,2] <- ScaledUpperCoords[ ,2]+1
# # ScaledLowerCoords[ ,2] <- ScaledLowerCoords[ ,2]-1
# # 
# # # DETECTORS
# # ScaledDetectors <- scaleCoordsToHabitatGrid(coordsData = detector.xy,
# #                                             coordsHabitatGridCenter = habitat$habitat.xy,
# #                                             scaleToGrid =T )$coordsDataScaled
# 
# 
# 
#
## ------     6.2. CREATE CACHED DETECTORS OBJECTS ------
#
# ## [CM] reduce multiplicator to 3 
# maxDistReCalc <- 2.1*detectors$maxDist #+ sqrt(2*(DETECTIONS$resizeFactor*HABITAT$habResolution)^2)
# 
# DetectorIndexLESS <- GetDetectorIndexLESS(
#   habitat.mx = habitat$habitat.mx,
#   detectors.xy = nimData$detector.xy,
#   maxDist = maxDistReCalc/res(habitat$habitat.r)[1],
#   ResizeFactor = 1,
#   plot.check = TRUE)
# 
# # ADD TO NIMDATA
# nimConstants$y.maxDet <- dim(DetectorIndexLESS$habitatID)[1]
# nimConstants$x.maxDet <- dim(DetectorIndexLESS$habitatID)[2]
# nimConstants$ResizeFactor <- DetectorIndexLESS$ResizeFactor
# nimConstants$n.cellsSparse <- dim(DetectorIndexLESS$detectorIndex)[1]
# nimConstants$maxNBDets <- DetectorIndexLESS$maxNBDets
#
# nimData$detectorIndex <- DetectorIndexLESS$detectorIndex
# nimData$nDetectorsLESS <- DetectorIndexLESS$nDetectorsLESS
# nimData$habitatIDDet <- DetectorIndexLESS$habitatID
#
#
#
## ------     6.3 TRANSFORM Y TO SPARSE MATRICES ------
#
# # STRUCTURED
# SparseY <- GetSparseY(y.aliveStructured)
# # ADD TO NIMDATA
# nimData$y.alive <- SparseY$y 
# nimData$yDets <- SparseY$yDets
# nimData$nbDetections <- SparseY$nbDetections
# nimConstants$nMaxDetectors <- SparseY$nMaxDetectors
# 
# # OTHER
# SparseYOth <- GetSparseY(y.aliveOthers)
# # ADD TO NIMDATA
# nimData$y.aliveOth <- SparseYOth$y 
# nimData$yDetsOth <- SparseYOth$yDets
# nimData$nbDetectionsOth <- SparseYOth$nbDetections
# nimConstants$nMaxDetectorsOth <- SparseYOth$nMaxDetectors



## ------   7. LOOP THROUGH INITIAL VALUES & SAVE OBJECT ------

# sxy
#create a data.frame with all detection of all Individuals detected
#project death to the next year
myData.deadProj <- myData.dead[ ,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1#project dead reco to the next year
#remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]

AllDets <- rbind(myData.alive$data.sp[ ,c("Id","Year")],
                 myData.deadProj[ ,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  coordsData = AllDetsxy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =T )$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)

idAugmented <- which(rownames(z) %in% "Augmented")
Id.vector <- y.ar$Id.vector
lowerCoords = nimData$lowerHabCoords
upperCoords = nimData$upperHabCoords
habitatGrid = nimData$habitatGrid

sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = lowerCoords,
                       upperCoords = upperCoords,
                       habitatGrid = habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")

# get location of individuals 
sxy.initscaled <- scaleCoordsToHabitatGrid(
  coordsData = sxy.init,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =F )$coordsDataScaled



## ------   8. calculate realized phi ------
#
# #initialize objects 
# recruit <- 0
# z_caculate <- nimData$z#[,1:(n.years-1)]
# z_caculate[is.na(z_caculate)] <- 0
# lev <- levels(habitatRasterResolution$'2.5km'$Countries)
# #EXTRACT LOCATION BASED ON INITIAL AC
# countryId <- list()
# for(t in 1:dim(z_caculate)[2]){
#   tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
#   countryId[[t]] <- raster::extract( habitatRasterResolution$'2.5km'$Countries ,tmp,sparse = F)
# }
# phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=n.years-1,ncol=length(lev[[1]]$ID))
# colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
#   colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'2.5km'$Countries,lev[[1]]$ID)[,1]
# for(c in 1:ncol(phi)){
#   for(t in 2:dim(z_caculate)[2]){
#     #phi
#     alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
#     phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
#     #culled
#     #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
#     #recru
#     notentered <- which(z_caculate[,t-1] == 0)
#     recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
#                               countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
#     recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
#                             countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
#                                                                      countryId[[t-1]] %in% c)
#   }
# }
# phi <- phi[,c(2,4)]        # overall phi
# recruit <- recruit[,c(2,4)]        # overall phi
# recruitnb <- recruitnb[,c(2,4)]        # overall phi
# 
# 
# ###
# pdf(file = file.path( working.dir, "figures/realizedPhiCountry.pdf"))
# par(mfrow = c(1,1))
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z")
# axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)]+1)
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(phi)){
#   points(phi[,c]~yr,pch=16,type="b", col=c)
# }
# legend("bottomright",colnames(phi),col=c(1:4),pch=16)
# 
# 
# par(mfrow = c(1,1))
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruit from z")
# axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)]+1)
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(recruit)){
#   points(recruit[,c]~yr,pch=16,type="b", col=c)
# }
# legend("topright",colnames(recruit),col=c(1:4),pch=16)
# dev.off()
# 
# 
# ###
# lev <- levels(habitatRasterResolution$'10km'$Counties)
# countryId <- list()
# for(t in 1:dim(z)[2]){
#   tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
#   countryId[[t]] <- raster::extract( habitatRasterResolution$'10km'$Counties ,tmp,sparse = F)
# }
# 
# phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=n.years-1,ncol=length(lev[[1]]$ID))
# colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
#   colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'10km'$Counties,lev[[1]]$ID)[,1]
# 
# for(c in 1:ncol(phi)){
#   for(t in 2:dim(z)[2]){
#     #phi
#     alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
#     phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
#     #culled
#     #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
#     #recru
#     notentered <- which(z_caculate[,t-1] == 0)
#     recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
#                               countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
#     recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
#                             countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
#                                                                      countryId[[t-1]] %in% c)
#   }
# }
# phi <- phi[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
# 
# phiNOR <- phi[,c(8,9,10,11,12,13)]
# phiSWE <- phi[,-c(8,9,10,11,12,13,14,15,16)]
# 
# recruitnb <- recruitnb[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
# recruitnbNOR <- recruitnb[,c(8,9,10,11,12,13)]
# recruitnbSWE <- recruitnb[,-c(8,9,10,11,12,13,14,15,16)]
# 
# recruit <- recruit[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
# recruitNOR <- recruit[,c(8,9,10,11,12,13)]
# recruitSWE <- recruit[,-c(8,9,10,11,12,13,14,15,16)]
# 
# 
# ###
# pdf( file = file.path(working.dir, "figures/realizedPhiCounties.pdf"),
#      width = 11, height = 6)
# # PHI
# #NORWAY
# par(mfrow = c(1,2))
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Norway")
# axis(1, at = 1:(n.years-1) , labels = paste( years[1:(n.years-1)]+1, years[1:(n.years-1)]+2,sep="-"))
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(phiNOR)){
#   points(phiNOR[,c]~yr,pch=16,type="b", col=c)
# }
# legend("bottomleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)
# 
# ###SWEDEN
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Sweden")
# axis(1, at = 1:(n.years-1) , labels = paste( years[1:(n.years-1)]+1, years[1:(n.years-1)]+2,sep="-"))
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(phiSWE)){
#   points(phiSWE[,c]~yr,pch=16,type="b", col=c)
# }
# legend("bottomleft",colnames(phiSWE),col=c(1:ncol(phiSWE)),pch=16)
# 
# 
# # RECRUITS
# #NORWAY
# par(mfrow = c(1,2))
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", 
#      ylab = "Realized recruitment from z", main="Norway")
# axis(1, at = 1:(n.years-1) , labels = paste( years[1:(n.years-1)]+1, years[1:(n.years-1)]+2,sep="-"))
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(recruitNOR)){
#   points(recruitNOR[,c]~yr,pch=16,type="b", col=c)
# }
# legend("topleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)
# 
# ###SWEDEN
# plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruitment from z", main="Sweden")
# axis(1, at = 1:(n.years-1) , labels = paste( years[1:(n.years-1)]+1, years[1:(n.years-1)]+2,sep="-"))
# yr <- c(1:(n.years-1))
# for(c in 1:ncol(recruitSWE)){
#   points(recruitSWE[,c]~yr,pch=16,type="b", col=c)
# }
# legend("topleft",colnames(recruitSWE),col=c(1:ncol(phiSWE)),pch=16)
# dev.off()
# 
# ##prop detected vs Alive in z
# propDet <- 0
# for(t in 1:n.years){
#   whichdets <- unique(c(which(nimData$nbDetections[,t]>0),
#                         which(nimData$nbDetectionsOth[,t]>0)))
#   whichAlive <- which(nimData$z[,t]%in%2)
#   propDet[t] <- length(whichdets)/length(whichAlive)
# }



for(c in 1:4){
  
  ##    3.4. LIST NIMBLE INITS ------
  
  nimInits <- list( "sxy" = sxy.init,
                    "dmean" = runif(1,0,10),
                    "z" = z.init,
                    "omeg1" = c(0.5,0.5),
                    "gamma" = runif(dim(y.alive)[3]-1,0,1),
                    "p01" = array(runif(18,0,0.2), c(nimConstants$n.counties,dim(y.alive)[3])),
                    "p01Oth" = array(runif(18,0,0.2), c(nimConstants$n.countries+1,dim(y.alive)[3])),
                    "sigma" = runif(n.years,1,4),
                    "betaDens" = runif(1,-0.1,0.1),
                    "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),c(dim(detCovsOth)[3],n.years)),
                    "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),c(dim(detCovsOth)[3],n.years)),
                    "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
                    "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
                    "detResponse" = InitsDetResponse,
                    "pResponse" = runif(1,0.4,0.5),
                    "phi" = runif(dim(y.alive)[3]-1,0.1,0.3)) 
  
  ### TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  idDEtected <- which(!rownames(z) %in%"Augmented")
  for(i in 1:length(idDEtected)){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i,,t]  
      }else{SXY <- nimData$sxy[i,,t]}
      sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
                                                        1:nimConstants$maxNBDets] 
      DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
      index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
      
      #table(detectorIndex)
      ## GET NECESSARY INFO 
      n.detectors <- length(index)
      #maxDist_squared <- maxDist*maxDist
      
      YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
      YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
      
      ## RECREATE Y
      if(nimData$nbDetections[i, t] > 0){
        for(j in 1:nimData$nbDetections[i, t]){
          ## check if a detection is out of the "detection window"
          if(sum(YDET[j]==index)==0){
            print(paste("id",i,"t",t,"j",j))
          }
        }
      }
    }}
  
  
  plot(nimData$detector.xy[,2]~nimData$detector.xy[,1])
  points(nimData$detector.xy[index,2]~nimData$detector.xy[index,1], col="red")
  points(SXY[2]~SXY[1], col="blue", pch=16)
  points(nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],2]~
           nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],1], col="green", pch=16)
  points(nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],2]~
           nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],1], col="purple", pch=16)
  
  
  plot(st_geometry(COUNTRIES))
  tmp <- myData.alive$data.sp[myData.alive$data.sp$Id %in% row.names(y.ar.ALIVE)[i] & myData.alive$data.sp$Year %in% years[t],]
  plot(st_geometry(tmp),col="red",add=T)
  
  
  nimInits$sxy <- round(nimInits$sxy, 5)#---an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
  
  
  ##CHECK WHERE IS NORRBOTTEN. IT IS ON THE 5TH INDEX
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(detectors$main.detector.sp[detCounties%in% c(1),]), col = myCol[5], pch = 16, cex = 0.8,add=T)
  
  nimConstants$countyToggle <- nimInits$p01
  nimConstants$countyToggle[] <- 1
  
  yearsNotSampled <- which(!years%in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countyToggle[1,t] <- 0
  }
  
  ## add another category to detcountry if in norrbotten, to turnoff detection to 0 there. 
  detCountriesNorb <- matrix(NA, nrow=length(detCountries),ncol=n.years)
  detCountries1 <- detCountries
  detCountries1[detCounties %in% 1] <- 3
  for(t in 1:n.years){
    if(t %in% yearsNotSampled){
      detCountriesNorb[,t] <- detCountries1
    }else{
      detCountriesNorb[,t] <- detCountries
    }
  }  
  
  nimData$detCountries <-  detCountriesNorb
  
  nimConstants$countryToggle <- nimInits$p01Oth
  nimConstants$countryToggle[] <- 1
  yearsNotSampled <- which(!years%in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countryToggle[3,t] <- 0
  }
  
  
  
  ## 7. SAVE NIMBLE INPUT ------
  
  save(nimData,
       nimConstants,
       y.dead,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(working.dir, "nimbleInFiles", 
                        paste0(modelName, "Chain", c, ".RData")))
}#c


## ------   9. SAVE NECESSARY OBJECTS ------
for(c in 1:4){
  
  ## List NIMBLE inits
  nimInits <- list(
    "sxy" = sxy.init,
    "dmean" = runif(1,0,10),
    "z" = z.init,
    "omeg1" = c(0.5,0.5),
    "gamma" = runif(dim(y.alive)[3]-1,0,1),
    "p01" = array(runif(18,0,0.2), c(nimConstants$n.counties,dim(y.alive)[3])),
    "p01Oth" = array(runif(18,0,0.2), c(nimConstants$n.countries+1,dim(y.alive)[3])),
    "sigma" = runif(n.years,1,4),
    "betaDens" = runif(1,-0.1,0.1),
    "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),c(dim(detCovsOth)[3],n.years)),
    "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),c(dim(detCovsOth)[3],n.years)),
    "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
    "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
    "detResponse" = InitsDetResponse,
    "pResponse" = runif(1,0.4,0.5),
    "phi" = runif(dim(y.alive)[3]-1,0.1,0.3)) 
  
  ## TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  idDEtected <- which(!rownames(z) %in%"Augmented")
  for(i in 1:length(idDEtected)){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i,,t]  
      }else{SXY <- nimData$sxy[i,,t]}
      sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
                                                        1:nimConstants$maxNBDets] 
      DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
      index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
      
      #table(detectorIndex)
      ## GET NECESSARY INFO 
      n.detectors <- length(index)
      #maxDist_squared <- maxDist*maxDist
      
      YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
      YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
      
      ## RECREATE Y
      if(nimData$nbDetections[i, t] > 0){
        for(j in 1:nimData$nbDetections[i, t]){
          ## check if a detection is out of the "detection window"
          if(sum(YDET[j]==index)==0){
            print(paste("id",i,"t",t,"j",j))
          }
        }
      }
    }}
  
  nimInits$sxy <- round(nimInits$sxy, 5)#---an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
  
  ##CHECK WHERE IS NORRBOTTEN. IT IS ON THE 5TH INDEX
  # plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  # plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  # plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  # plot(st_geometry(detectors$main.detector.sp[detCounties%in% c(1),]), col = myCol[5], pch = 16, cex = 0.8,add=T)
  # 
  nimConstants$countyToggle <- nimInits$p01
  nimConstants$countyToggle[] <- 1
  
  yearsNotSampled <- which(!years%in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countyToggle[1,t] <- 0
  }
  
  ## add another category to detCountry if in Norrbotten, to turnoff detection to 0 there. 
  detCountriesNorb <- matrix(NA, nrow=length(detCountries),ncol=n.years)
  detCountries1 <- detCountries
  detCountries1[detCounties %in% 1] <- 3
  for(t in 1:n.years){
    if(t %in% yearsNotSampled){
      detCountriesNorb[,t] <- detCountries1
    }else{
      detCountriesNorb[,t] <- detCountries
    }
  }  
  
  nimData$detCountries <- detCountriesNorb
  
  nimConstants$countryToggle <- nimInits$p01Oth
  nimConstants$countryToggle[] <- 1
  yearsNotSampled <- which(!years %in% yearsSampledNorrb)
  for(t in yearsNotSampled){
    nimConstants$countryToggle[3,t] <- 0
  }
  
  
  
  ## 7. SAVE NIMBLE INPUT ------
  
  save(nimData,
       nimConstants,
       y.dead,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(working.dir, "nimbleInFiles", 
                        paste0(modelName, "Chain", c, ".RData")))
}#c
myHabitat.list <- habitat
myDetectors <- myDetectors
myStudyArea.poly <- myStudyArea


save( myHabitat.list,
      myDetectors,
      COUNTRIES,
      myStudyArea.poly,
      COMMUNES,
      COUNTIES_AGGREGATEDSubset,
      myFilteredData.sp,
      myFullData.sp,
      COUNTIES_AGGREGATED,
      file = file.path(working.dir, "data/NecessaryObjects.RData"))




## END OF makeRovquantData_wolverine() ------




##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

## PICK UP HERE NEXT TIME !

##------------------------------------------------------------------------------
##------------------------------------------------------------------------------