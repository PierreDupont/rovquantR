##--------------------------------------------------------------------------- --
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
##--------------------------------------------------------------------------- --
##
## Notes: 
## This is based on 'rovquantR' beta version 0.1
##   
##--------------------------------------------------------------------------- --

## CLEAR-UP ENVIRONMENT ------

rm(list = ls())
gc()


## INSTALL 'rovquantR' FROM GITHUB ------

## Ctrl + Shift + F10 (to restart R session)
devtools::install_github("PierreDupont/rovquantR@devel")


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
#source("C:/My_documents/RovQuant/Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R")
source("C:/My_documents/RovQuant/Source/DoScale.R")



##------------------------------------------------------------------------------

## 1. GENERAL VARIABLES DECLARATION ------

years <- 2014:2023
n.years <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))
species <- "Jerv"
#load("C:/My_documents/rovquantR/R/sysdata.rda")
sampling.months <- list(12,1:6)

two.sex <- TRUE
legal.dead <- NULL
overwrite <- FALSE
plot.check = TRUE


##------------------------------------------------------------------------------

## cleanRovBaseData() ------

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



## ------   2. CLEAN THE DATA -----

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
    Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
  ##-- Filter to the focal years
  dplyr::filter(., Year %in% years)

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
                  Death_cause)) %>%
  ##-- Filter to the focal years
  dplyr::filter(., Year %in% years)

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

# [PD] This is not used in this analysis (and most likely not working)
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
## [PD] : MIGHT WANT TO KEEP SAMPLES OUTSIDE NOR & SWE IN cleanRovbaseData()?
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



## ------       3.1.3. HAIR TRAPS -----

# HairTrapSamples <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/hairtrapsNB2024.xlsx"))#, fileEncoding="latin1") ## DNA samples to be removed from Henrik
HairTrapSamples <- readMostRecent(
  path = data.dir,
  extension = ".xlsx",
  pattern = "hairtrap")

##-- Identify samples collected by hair traps
alive <- alive %>%
  mutate(hairTrap = DNAID %in% HairTrapSamples$DNAID)



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
max.det.dist <- 840000
resize.factor <- 1
sex <- c("female","male")

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
# ##-- Plot check
# if(plot.check){
#   par(mfrow = c(1,1))
#   plot(st_geometry(COUNTRIES))
#   plot(st_geometry(myStudyArea), add = TRUE, col = "red")
# }



## ------   3. NGS DATA -----

##-- Extract date from the last cleaned data file
DATE <- getMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolverine")

##-- Load the most recent wolverine data from RovBase
myFullData.sp <- readMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolverine",
  extension = ".RData")

##-- Define years
if(is.null(years)){
  years <- sort(unique(c(myFullData.sp$alive$Year,
                         myFullData.sp$dead.recovery$Year)))
}
data$years <- years
n.years <- length(years)

##-- list years with or without sampling in Norrbotten
yearsSampledNorrb <- c(2016:2018,2023)
yearsNotSampled <- years[!years %in% yearsSampledNorrb]
whichYearsNotSampled <- which(years %in% yearsNotSampled)

##-- Filter NGS samples for dates
myFullData.sp$alive <- myFullData.sp$alive %>%
  dplyr::filter(
    ##-- Subset to years of interest
    Year %in% years,
    ##-- Subset to monitoring period
    Month %in% unlist(sampling.months))

##-- Filter Dead recoveries for dates
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
  ##-- Subset to years of interest
  dplyr::filter(Year %in% years)



## ------ II. CREATE OPSCR DATA ------

## ------   1. GENERATE HABITAT ------

message("Preparing habitat characteristics... ")

## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------

##-- Determine study area based on NGS detections
##-- Buffer NGS detections and cut to Swedish and Norwegian borders
studyArea <- myFullData.sp$alive %>%
  sf::st_buffer(., dist = habitat$buffer * 1.4) %>%
  mutate(id = 1) %>%
  group_by(id) %>% 
  summarize() %>% 
  sf::st_intersection(., COUNTRIES) %>%
  sf::st_as_sf()

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster(
  poly = studyArea,
  habitat.r = habRaster,
  buffer = habitat$buffer,
  plot.check = FALSE) %>%
  append(habitat,.)

##-- Retrieve number of habitat windows 
isHab <- habitat$habitat.r[] == 1
n.habWindows <- habitat$n.habWindows <- sum(isHab)
habitat$habitat.df <- cbind.data.frame(
  "id" = 1:habitat$n.habWindows,
  "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
  "y" = raster::coordinates(habitat$habitat.r)[isHab,2])

##-- Make a spatial grid from polygon
habitat$grid <- sf::st_as_sf( stars::st_as_stars(habitat$habitat.r), 
                              as_points = FALSE,
                              merge = FALSE) %>%
  filter( Habitat %in% 1) %>%
  st_set_crs( .,value = sf::st_crs(habitat$buffered.habitat.poly)) %>%
  mutate( id = 1:nrow(.),
          x = st_coordinates(st_centroid(.))[ ,1],
          y = st_coordinates(st_centroid(.))[ ,2]) 

##-- Study area grid from habitat raster
habitat.rWthBufferPol <- sf::st_as_sf( 
  stars::st_as_stars(habitat$habitat.rWthBuffer), 
  as_points = FALSE,
  merge = TRUE) %>%
  filter(Habitat %in% 1)

##-- Plot check 
par(mfrow = c(1,2))
plot(habitat$habitat.r,legend = F)
plot(st_geometry(studyArea), add = T)
plot(st_geometry(habitat$grid), add = T)


# ## [PD] DONE LATER
# ##-- Retrieve habitat windows boundaries
# lowerHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1, ] - 0.5*habitat$resolution
# upperHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1, ] + 0.5*habitat$resolution
# nHabCells <- dim(lowerHabCoords)[1]
# 
# ## [PD] DONE LATER
# ##-- CREATE HABITAT GRID
# habIDCells.mx <- habitat$IDCells.mx
# habIDCells.mx[] <- 0
# for(i in 1:nrow(lowerHabCoords)){
#   habIDCells.mx[trunc(lowerHabCoords[i,2])+1,
#                 trunc(lowerHabCoords[i,1])+1] <- i
# }
# # image(habIDCells.mx)
#
# ## [PD] DONE DURING THE INITIAL CLEANING
# whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyArea))))
# if(length(whichOut) > 0){
#   myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
# }
# #myFilteredData.sp$alive$Id <- droplevels( myFilteredData.sp$alive$Id)
# ## REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT #[CM] 
# whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, habitat$buffered.habitat.poly))))
# if(length(whichOutBuff) > 0){ myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]}

##-- Plot check
if(plot.check){
  # par(mfrow = c(1,2))
  plot(habitat$habitat.r)
  plot(st_geometry(myStudyArea),
       add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
  plot(st_geometry(GLOBALMAP), add = T)
  plot(st_geometry(habitat$buffered.habitat.poly), add = T)
  plot(st_geometry(myFullData.sp$alive),
       pch = 21, bg = "red", cex = 0.5, add = T)
  plot(st_geometry(myFullData.sp$dead.recovery),
       pch = 21, bg = "blue", cex = 0.5, add = T)
}



## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       1.2.1. DEN COUNTS ------

##-- Load the last DEN COUNT data file
#DEN <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/DEN_COUNTS_2009_2024_fromHB.csv"), fileEncoding="latin1")
DEN <- readMostRecent( 
  path = data.dir,
  extension = ".csv",
  pattern = "DEN_COUNTS") %>%
  st_as_sf(., coords = c("UTM33_X", "UTM33_Y")) %>%
  st_set_crs(value = st_crs(myFullData.sp$alive)) %>%
  mutate(id = 1)

# colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN))
# DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
# st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
# DEN.sp$id  <- rep(1, nrow(DEN.sp))
# DEN.sp <- DEN.sp[ ,("id")]

DEN.r <- raster(
  adehabitatHR::estUDm2spixdf(
    adehabitatHR::kernelUD( as(DEN[ ,"id"], "Spatial"),
                            h = 30000,
                            grid = as(habitat$habitat.r, 'SpatialPixels'))))

##-- Plot check
if(plot.check){
  plot(DEN.r)
  plot(st_geometry(myStudyArea), add = TRUE, border = "black")
}

##-- EXTRACT COVARIATES
denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
denCounts <- as.vector(round(scale(denCounts), digits = 2))



## ------   2. GENERATE DETECTORS -----

message("Preparing detectors characteristics... ")

## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----

##-- Generate NGS detectors based on the study area 
detectors$subdetectors.r <- disaggregate(
  habitat$habitat.rWthBuffer,
  fact = res(habitat$habitat.r)[1]/detectors$resolution.sub)

##-- Generate NGS detectors based on the raster of sub-detectors
detectors <- makeSearchGrid( 
  data = detectors$subdetectors.r,
  resolution = detectors$detResolution,
  div = (detectors$resolution/detectors$resolution.sub)^2,
  plot = FALSE) %>%
  append(detectors, .)

##-- Format detector locations & number of trials per detector
detectors$detectors.df <- cbind.data.frame(
  "id" = 1:nrow(detectors$main.detector.sp),
  "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1],
  "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2],
  "size" = detectors$main.detector.sp$count)

##-- Generate detector raster 
detectors$raster <- raster::rasterFromXYZ(
  cbind( detectors$detectors.df[ ,c("x","y")],
         rep(1,nrow(detectors$main.detector.sp))))

##-- Make a spatial grid from detector raster
detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
  x = detectors$raster,
  fun = function(x){x>0})) %>%
  st_set_crs(.,value = st_crs(studyArea)) %>%
  mutate( id = 1:nrow(.),
          x = st_coordinates(st_centroid(.))[ ,1],
          y = st_coordinates(st_centroid(.))[ ,2]) 

##-- Extract numbers of detectors
n.detectors <- detectors$n.detectors <- dim(detectors$main.detector.sp)[1]


##-- Identify detectors in Norrbotten 
COUNTIESAroundNorrbotten <- REGIONS %>%
  group_by(county) %>%
  summarize() %>%
  filter(county %in% c("Norrbotten","Troms","Västerbotten","Nordland","Finnmark")) %>% 
  st_simplify( dTolerance = 500)

##-- Create an index of detectors in Norrbotten
distDetsCounties <- st_distance( detectors$main.detector.sp,
                                 COUNTIESAroundNorrbotten,
                                 byid = T)
detsNorrbotten <- which(apply(distDetsCounties, 1, which.min) == 3)


##-- Plot check
if(plot.check){
  ##-- Plot detectors in Norrbotten
  plot( st_geometry(COUNTIESAroundNorrbotten))
  plot( st_geometry(detectors$main.detector.sp),
        col = "black", pch = 16, cex = 0.3, add = T)
  plot( st_geometry(detectors$main.detector.sp[detsNorrbotten, ]),
        col = "red", pch = 16, cex = 0.5, add = T)
  
  ##-- Plot NGS detectors
  plot( st_geometry(habitat$buffered.habitat.poly),
        main = paste(detectors$n.detectors, "Detectors"),
        col = rgb(0.16,0.67,0.16, alpha = 0.3))  
  plot( st_geometry(studyArea), add = TRUE,
        col = rgb(0.16,0.67,0.16, alpha = 0.5))
  plot( st_geometry(detectors$main.detector.sp),
        col = "red", pch = 16, cex = 0.1, add = TRUE)
  plot( st_geometry(COUNTRIES), add = TRUE)
}



## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----

## ------       2.2.1. EXTRACT COUNTIES ------

##-- Extract closest county for each detector
detCounties <- detectors$main.detector.sp %>%
  st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
  apply(., 1, function(x) which.min(x))

##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten 
##-- in years without sampling
countyToggle <- matrix(1, nrow = max(detCounties), ncol = n.years)
for(t in whichYearsNotSampled){
  countyToggle[1,t] <- 0
}

##-- Plot check 
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
  text(st_geometry(COUNTIES_AGGREGATED), labels = COUNTIES_AGGREGATED$id, col = "black")  
}



## ------       2.2.2. EXTRACT COUNTRIES ------

##-- Extract closest country for each detector
detCountries <- detectors$main.detector.sp %>%
  st_distance(., COUNTRIES, by_element = F) %>%
  apply(., 1, function(x) which.min(x)) %>%
  as.factor(.) %>%
  as.numeric(.)

##-- Turn into a matrix with years in columns
detCountries <- matrix( detCountries,
                        nrow = length(detCountries),
                        ncol = n.years)

##-- Add another category to detCountry if in Norrbotten, to turnoff detection to 0 there. 
for(t in whichYearsNotSampled){
  detCountries[detCounties %in% 1,t] <- 3
}#t  

##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten in years without sampling
countryToggle <- matrix(1, nrow = max(detCountries), ncol = n.years)
for(t in whichYearsNotSampled){
  countryToggle[3,t] <- 0
}

##-- Plot check 
if(plot.check){
  par(mfrow = c(1,1))
  myCol <- c("blue4", "yellow1", "red")
  plot( st_geometry(GLOBALMAP), col = "gray80", main = "Countries")
  plot( st_geometry(myStudyArea),
        col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot( st_geometry(detectors$main.detector.sp),
        col = myCol[detCountries[,1]], pch = 16, cex = 0.8, add = T)
  plot( st_geometry(COUNTRIES), add = TRUE)
}



## ------       2.2.3. EXTRACT GPS TRACKS LENGTHS ------

message("Cleaning GPS tracks... ")

## LOAD NEW GPS SEARCH TRACKS !!!
# TRACKS_SINGLE <- read_sf(file.path(data.dir,
#                                    "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp"))
# TRACKS_MULTI <- read_sf(file.path(data.dir,
#                                   "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp"))
##-- Combine all GPS tracks
TRACKS <- rbind(
  read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp")),
  read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))) %>%
  ##-- Process dates
  mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
          Mth = as.numeric(format(Dato,"%m")),
          Yr = as.numeric(format(Dato,"%Y")),
          Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
  ##-- Filter out irrelevant tracks
  filter( Helikopter == "0",      ## Remove helicopter tracks
          Jerv == "1",            ## Keep Wolverine tracks only
          Year %in% years & Mth %in% unlist(sampling.months)) %>% ## Keep tracks during sampling season only
  ##-- Extract track lengths & centroids
  mutate( Length = st_length(., byid = T),
          Centroidx = st_coordinates(st_centroid(.))[ ,1]) 

##-- Find & filter out duplicates based on person, distance and date.
df <- data.frame( Dato = TRACKS$Dato,
                  Year = TRACKS$Year,
                  Person = TRACKS$Person,
                  Length = TRACKS$Length,
                  Centroidx = TRACKS$Centroidx)
dupIDs <- which(duplicated(df))
dupLength <- TRACKS$Length[duplicated(df)]
TRACKS <- TRACKS[-dupIDs, ]

# ## PLOT CHECK
# if(plot.check){
# 
#   par(mfrow = c(2,2))
#   ## Length of tracks searched per year
#   lengthPerYear <- unlist(lapply(years,function(x) sum(TRACKS$Length[TRACKS$Year == x])/1000))
#   names(lengthPerYear) <- years
#   barplot(lengthPerYear, ylab = "Track length (km)", main = "Length of tracks searched per year")
# 
#   ## Number of tracks searched per year
#   numPerYear <- unlist(lapply(years,function(x) sum(TRACKS$Year == x)))
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
# ##-- save 
# save( TRACKS, file = file.path(working.dir, "data/searchTracks.RData"))
# load(file = file.path(working.dir, "data/searchTracks.RData"))


##-- Extract length of GPS search track per detector grid cell
detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
TRACKS.r <- list()
for(t in 1:n.years){
  intersection <- TRACKS %>%
    dplyr::filter(Year == years[t]) %>%
    sf::st_intersection(detectors$grid, .) %>%
    dplyr::mutate(LEN = st_length(.)) %>%
    sf::st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN)) ##-- Get total length searched in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectors$raster
  TRACKS.r[[t]][detectors$raster[] %in% 1] <- detTracks[ ,t]
  print(t)
}#t


##-- Plot check 
if(plot.check){
  max <- max(unlist(lapply(TRACKS.r, function(x) max(x[], na.rm = T))))
  cuts <- seq(0, max,length.out = 100)   ## Set breaks
  col <- rev(terrain.colors(100))
  CountriesDetRes <- disaggregate(habitatRasters$Countries, fact = 2)
  CountriesDetRes <- crop(CountriesDetRes, TRACKS.r[[1]])
  country.colors <- c("goldenrod1","goldenrod3")
  
  pdf(file = file.path(working.dir, "figures/Tracks.pdf"))
  NORTRACKS <- SWETRACKS <- 0
  for(t in 1:n.years){
    plot( TRACKS.r[[t]], main = years[t], breaks = cuts, col = col, legend = FALSE)
    plot( st_geometry(habitat$habitat.poly), main = years[t], add = T)
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



## ------       2.2.4. EXTRACT DISTANCES TO ROADS ------

##-- Load map of distance to roads (1km resolution)
DistAllRoads <- raster::raster(file.path(data.dir,"GIS/Roads/MinDistAllRoads1km.tif"))

##-- Fasterize to remove values that fall in the sea
r <- fasterize::fasterize(sf::st_as_sf(COUNTRIES), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea)
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

##-- Plot check 
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



## ------       2.2.5. EXTRACT DAYS OF SNOW ------

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

##-- if NA returns the average value of the cells within 15km 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp

##-- Plot check
if(plot.check){
  plot( st_geometry(detectors$main.detector.sp),
        cex = DoScale(detSnow[ ,6], l = 0, u = 0.5),
        pch = 16)
}



## ------       2.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------         2.2.6.1. SKANDOBS ------

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
                 species = stringi::stri_trans_general(species, "Latin-ASCII"),
                 monitoring.season = ifelse( month < unlist(sampling.months)[1],
                                             year, year + 1)) %>%
  ##-- Filter based on monitoring season
  dplyr::filter( month %in% unlist(sampling.months)) %>%
  ##-- Turn into spatial points object
  sf::st_as_sf(., coords = c("longitude","latitude")) %>%
  sf::st_set_crs(., value = "EPSG:4326") %>%
  sf::st_transform(., sf::st_crs(COUNTIES)) %>%
  dplyr::filter(!is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))



## ------         2.2.6.2. ROVBASE ------

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


##-- Process Rovbase observations (all species)
rovbaseObs <- rbind(rovbaseObs1, rovbaseObs2, rovbaseObs3, rovbaseObs4) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>%
  ##-- Extract important info (e.g. month, year, country of collection)
  dplyr::mutate(
    ##-- Turn potential factors into characters 
    across(where(is.factor), as.character),
    ##-- Deal with Scandinavian characters
    Species = stringi::stri_trans_general(Species, "Latin-ASCII"),
    Sample_type = translateForeignCharacters(data = Sample_type),
    ##-- Deal with dates
    Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
    year = as.numeric(format(Date,"%Y")),
    month = as.numeric(format(Date,"%m")),
    monitoring.season = ifelse(month < 12, year, year+1)) %>%
  ##-- Filter out unusable samples
  dplyr::filter( 
    ##-- Filter out samples without coordinates,...
    !is.na(East_UTM33),
    ##-- ...based on species
    # Species %in% c("Bjorn","Fjellrev","Gaupe","Hund","Jerv","Rodrev","Ulv"),
    ##-- ...based on sample type
    Sample_type %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
                        "Saliv/Spytt", "Loepeblod", "Vev"),
    ##-- ...based on monitoring season
    month %in% unlist(sampling.months),
    ##-- ... if sample was from the focal species and successfully genotyped 
    !(Species %in% "Jerv" & !is.na(Id))) %>%
  ##-- Turn into spatial points object
  sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(. , sf::st_crs(COUNTIES)) %>%
  filter( 
    ##-- ... based on space 
    !is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))

##-- Remove un-necessary objects
rm(list = c("rovbaseObs1","rovbaseObs2","rovbaseObs3","rovbaseObs4"))



## ------         2.2.6.3. COMBINE ROVBASE & SKANDOBS ------

##-- RASTERIZE AT THE DETECTOR LEVEL
r.list <- lapply(years, function(y){
  ##-- Rasterize Skandobs observations 
  sk.r <- raster::rasterize(
    skandObs[skandObs$monitoring.season %in% y, 1],
    detectors$raster,
    fun = "count")[[1]]
  sk.r[is.na(sk.r[])] <- 0
  ##-- Set cells outside detector area to NA
  sk.r[!detectors$raster[ ]%in% 1] <- NA
  ##-- Turn into binary raster
  sk.r1 <- sk.r
  sk.r1[sk.r1[]>0] <- 1
  
  ##-- Rasterize Rovbase observations 
  rb.r <- raster::rasterize(
    rovbaseObs[rovbaseObs$monitoring.season %in% y, 1],
    detectors$raster,
    fun = "count")[[1]]
  rb.r[is.na(rb.r[])] <- 0
  ##-- Set cells outside detector area to NA
  rb.r[! detectors$raster[ ]%in% 1] <- NA
  ##-- Turn into binary raster
  rb.r1 <- rb.r
  rb.r1[rb.r1[]>0] <- 1
  
  list(sk.r, sk.r1, rb.r, rb.r1)
})

##-- Store in raster bricks
r.skandObsBinary <- brick(lapply(r.list,function(x) x[[2]]))
r.skandObsContinuous <- brick(lapply(r.list,function(x) x[[1]]))
r.rovbaseBinary <- brick(lapply(r.list,function(x) x[[4]]))
r.rovbaseContinuous <- brick(lapply(r.list,function(x) x[[3]]))


##-- Combine both rasters
r.SkandObsRovbaseBinary <- r.rovbaseBinary + r.skandObsBinary
for(t in 1:n.years){
  r.SkandObsRovbaseBinary[[t]][r.SkandObsRovbaseBinary[[t]][]>1 ] <- 1
}


##-- Plot check
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
  
  # pdf(file = file.path(working.dir, "figures/mapStructuredOthers.pdf"))
  # for(t in 1:n.years){
  #   tmpOthers <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t] &
  #                                          !myFilteredData.sp$alive$structured, ]
  #   tmpStruct <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t] &
  #                                          myFilteredData.sp$alive$structured, ]
  # 
  #   par(mfrow=c(2,2),mar=c(0,0,5,0))
  #   plot(r.OtherSamplesBinary[[t]], main=paste(years[t],"\n Rovbase Samples Structured"), box=F, axes=F)
  #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  #   plot(r.OtherSamplesBinary[[t]],main=paste(years[t],"\n Rovbase Samples Opportunistic"), box=F, axes=F)
  #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
  # 
  #   plot(r.skandObsSamplesBinary[[t]], main=paste(years[t]), box=F, axes=F)
  #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  #   plot(r.skandObsSamplesBinary[[t]],main=paste(years[t],"\n SkandObs Opportunistic"), box=F, axes=F)
  #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
  # }
  # dev.off()
  
  for(t in 1:n.years){
    par(mfrow=c(1,3),mar=c(0,0,5,0))
    plot(r.rovbaseBinary[[t]],main=years[t])
    plot(r.skandObsBinary[[t]])
    plot(r.SkandObsRovbaseBinary[[t]])
  }#t
}



# ## ------         2.2.6.4. SMOOTH THE BINARY MAP ------
# 
# ##-- We tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
# habOwin <- spatstat.geom::as.owin(as.vector(extent(detectors$raster)))
# cutoff <- 1
# ds.list <- lapply(years,function(y){
#   ## ROVBASE DATA 
#   pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y,]
#   ## SKANDOBS
#   pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y,] )
#   ## SMOOTH AND RASTERIZE
#   p <-  spatstat.geom::ppp(pts[,1], pts[,2], window = habOwin)
#   ds <- density(p, adjust=0.02) #---change bandwith (smoothing) with "adjust
#   ds <- raster(ds)
#   
#   ds <- ds1 <- raster::resample(ds, detectors$raster) #mask(ds,rasterToPolygons(myHabitat.list$habitat.rWthBuffer,function(x) x==1))
#   threshold <- 0.1 / prod(res(ds)) #--number per 1 unit of the projected raster (meters)
#   ds1[] <- ifelse(ds[]<threshold,0,1)
#   ds1 <- mask(ds1, habitat.rWthBufferPol)
#   ds <- mask(ds, habitat.rWthBufferPol)
#   
#   return(list(ds,ds1))
# })
# 
# ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
# ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))
# names(ds.brick) <- names(ds.brickCont) <-years
# 
# ##-- Plot check
# if(plot.check){
#   par(mfrow = c(1,3))
#   plot(r.SkandObsRovbaseBinary[[t]], main = "Raw Binary", axes = F, box = F)
#   plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
#   plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
# }



## ------         2.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------

##-- IDENTIFY HAIR SAMPLES
tmpHair <- myFullData.sp$alive %>% filter(hairTrap)

##-- MANUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL.
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr)
  whereHair <- raster::extract( r.SkandObsRovbaseBinary[[t]],
                                tmpHair,
                                cellnumbers = T)
  r.SkandObsRovbaseBinary[[t]][whereHair[ ,1]] <- 1
}#t



## ------         2.2.6.6. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
detOtherSamples[ ,1:n.years] <- raster::extract( r.SkandObsRovbaseBinary,
                                                 detectors$main.detector.sp)
colSums(detOtherSamples)



## ------       2.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 2))
detCovs[,,1] <- detTracks
detCovs[,,2] <- detSnow

detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2], 3))
detCovsOth[,,1] <- detSnow
detCovsOth[,,2] <- matrix(detRoads,length(detRoads),n.years)
detCovsOth[,,3] <- detOtherSamples


##-- CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}


##-- Plot check
if(plot.check){
  tmp <- detectors$raster
  par(mfrow=c(2,5),mar=c(0,0,0,0))
  max <- max(detCovsOth[,,2])
  cuts <- seq(0,max,length.out = 100) 
  col <- rev(terrain.colors(100))
  for(t in 1:n.years){
    plot(detectors$raster, col=c(grey(0.2),grey(0.8)),axes=F,legend=F,box=F,)
    tmp[!is.na(detectors$raster[ ])] <- detCovsOth[,t,2]
    plot(tmp,axes=F,legend=F,box=F,breaks = cuts, col=col,add=T)
  }
  dev.off()
  
  pdf(file = file.path(working.dir, "figures/detections over space and time.pdf"))
  for(t in 1:n.years){
    ## NGS DETECTIONS TOTAL
    tempTotal <- myFullData.sp$alive[myFullData.sp$alive$Year == years[t], ]
    NGS_TabTotal <- table(tempTotal$Country_sample)
    ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country_sample), 2, function(x) sum(x>0))
    ## ALIVE DETECTIONS INSIDE STUDY AREA/SAMPLING PERIOD
    tempIn <- myFullData.sp$alive[myFullData.sp$alive$Year == years[t], ]
    NGS_TabIn <- table(tempIn$Country_sample)
    ID_TabIn <- apply(table(tempIn$Id, tempIn$Country_sample), 2, function(x) sum(x>0))
    ## PLOT NGS SAMPLES
    plot(st_geometry(GLOBALMAP), col="gray80")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
    plot(st_geometry(habitat$buffered.habitat.poly), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
    points(tempTotal, pch = 21, bg = "darkred")
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



## ------   3. RESCALE COORDINATES------

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



## ------   4. CREATE LOCAL OBJECTS -----

# DetectorIndexLESS <- GetDetectorIndexLESS(
#   habitat.mx = habitat$habitat.mx,
#   detectors.xy = detectors$scaledCoords,
#   maxDist = detectors$maxDist/habitat$resolution,
#   ResizeFactor = 1,
#   plot.check = TRUE)

##-- Get local detectors
detectors$localObjects <- getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = detectors$scaledCoords,
  dmax = detectors$maxDist/habitat$resolution,
  resizeFactor = detectors$resize.factor,
  plot.check = F)

# DetectorIndexLESS == detectors$localObjects 



## ------   5. SAVE STATE-SPACE CHARACTERISTICS -----

save( habitat,
      file = file.path( working.dir, "data",
                        paste0("Habitat_wolverine_", DATE, ".RData")))

save( detectors,
      file = file.path( working.dir,"data",
                        paste0("Detectors_wolverine_", DATE, ".RData")))



## ------   6. GENERATE y DETECTION ARRAYS ------

for(thisSex in sex){
  
  message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
  
  ## ------     6.1. FILTER DATA BY SEX -----
  
  # load(file.path( working.dir, "data",
  #                 paste0("FilteredData_bear_", DATE, ".RData")))
  myFilteredData.sp <- myFullData.sp
  
  myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
    dplyr::filter(Sex %in% thisSex)
  
  myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery %>%
    dplyr::filter(Sex %in% thisSex)
  
  
  
  # ## ------     6.1. FILTER NGS & DEAD RECOVERY DATA FOR DATES ------
  # 
  # # ##-- Check correlation number of detections ~ between monitoring season
  # # deadID <- unique(myFilteredData.sp$dead.recovery$Id)
  # # ndet <- NULL
  # # timeDiff <- NULL
  # # for(i in 1:length(deadID)){
  # #   tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i], ]$Year
  # # 
  # #   timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i], ]$Date-
  # #     as.POSIXct(strptime(paste("01-12", tmpYear, sep = "-"), "%d-%m-%Y"))
  # # 
  # #   ndet[i] <- nrow(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] &
  # #                                             myFilteredData.sp$alive$Year %in% tmpYear, ])
  # # }
  # # 
  # # pdf(file = file.path(working.dir, "figures/Prop id detected_Time available.pdf"))
  # # plot( ndet ~ timeDiff,
  # #       ylab = "Total number of detections",
  # #       xlab = "Number of days between dec 1 and dead recovery")
  # # hh <- hist(timeDiff[ndet > 0], breaks = seq(0,400,by=25))
  # # hh1 <- hist(timeDiff[ndet == 0], breaks = seq(0,400,by=25))
  # # barplot(rbind(hh$counts/(hh$counts+hh1$counts),
  # #               hh1$counts/(hh$counts+hh1$counts)),
  # #         names.arg = hh$breaks[1:(length(hh$breaks)-1)],
  # #         xlab = "number of days between dead reco and start monitoring",
  # #         ylab = "%")
  # # legend( "topright",
  # #         fill = c(grey(0.2), grey(0.8)),
  # #         legend = c("detected","notDetected"))
  # # dev.off()
  # # 
  # # 
  # # ##-- Filter NGS samples for dates
  # # myFullData.sp$alive <- myFullData.sp$alive %>%
  # #   dplyr::filter(
  # #     ##-- Subset to years of interest
  # #     Year %in% years,
  # #     ##-- Subset to monitoring period
  # #     Month %in% unlist(sampling.months))
  # # 
  # # ##-- Filter Dead recoveries for dates
  # # myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
  # #   ##-- Subset to years of interest
  # #   dplyr::filter(Year %in% years)
  # 
  # 
  ## ------     6.2. FILTER OUT DETECTIONS IN NORRBOTTEN EXCEPT IN 2016:18 and 2023 ------
  
  # ##-- Get Norrbotten borders
  COMMUNES_NOR <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp")) ## Communal map of Norway
  COMMUNES_SWE <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp")) ## Communal map of Sweden
  COUNTIESNorrbotten <- rbind(COMMUNES_NOR, COMMUNES_SWE) %>%
    filter(NAME_1 %in% "Norrbotten") %>%
    group_by(NAME_1) %>%
    summarize()
  
  ## [PD] : USING REGIONS INSTEAD OF COMMUNES WILL LEAD TO DIFFERENT NUMBERS OF
  ## SAMPLES REMOVED IN DIFFERENT YEARS BECAUSE OF SLIGHT DIFFERENCES IN SHAPEFILES
  ## (APPROX. 10 samples OVERALL)
  # COUNTIESNorrbotten <- REGIONS %>%
  #   filter(county %in% "Norrbotten") %>%
  #   group_by(county) %>%
  #   summarize()
  
  ##-- Check how many detections have been collected in Norrbotten overall
  is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))
  sum(is.Norr, na.rm = T)
  
  ##-- Check how many detections are removed per year
  table(myFilteredData.sp$alive$Year[myFilteredData.sp$alive$Year %in% yearsNotSampled & is.Norr %in% 1])
  sum(myFilteredData.sp$alive$Year %in% yearsNotSampled & is.Norr %in% 1)
  
  ##-- Subset NGS dataset
  myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
    filter(!(Year %in% yearsNotSampled & is.Norr %in% 1))
  
  ##-- Checks
  dim(myFilteredData.sp$alive)
  table(myFilteredData.sp$alive$Year)
  
  # ## plot check
  # for(t in 1:n.years){
  #   plot( st_geometry(myStudyArea))
  #   plot( st_geometry(COUNTIESNorrbotten), add = T, col = "blue")
  #   plot( st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t], ]),
  #         col = "red", add = T, pch = 16)
  # }
  
  
  
  ## ------     6.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------
  
  ## ------       6.3.1. ASSIGN SAMPLES TO GPS SEARCH TRACKS ------
  
  message("Assigning DNA samples to GPS tracks... ")
  message("This can take several minutes... ")
  
  myFilteredData.sp$alive <- assignSearchTracks(
    data = myFilteredData.sp$alive,
    tracks = TRACKS)
  
  # ##-- SAVE FOR FASTER LOADING
  # save(myFilteredData.sp, file = file.path(working.dir, "data/myFilteredData.sp.RData"))
  # load(file.path(working.dir, "data/myFilteredData.sp.RData"))
  
  
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
  
  
  
  ## ------       6.3.2. ASSIGN SAMPLES TO OPPORTUNISTIC OR STRUCTURED ------
  
  distanceThreshold <- 500
  
  ##-- Identify samples from structured and opportunistic sampling
  myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
    mutate(
      ##-- Collector column was replaced by two columns, merging them now...
      Collector_role = ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
      ##-- Identify samples collected during structured sampling 
      structured = Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
        !is.na(trackID) &
        trackDist <= distanceThreshold & 
        !hairTrap)
  table(myFilteredData.sp$alive$Collector_role, useNA = "always")
  table(myFilteredData.sp$alive$structured, useNA = "always")
  
  ##-- Plot check
  if(plot.check){
    plot(REGIONS[REGIONS$county %in% "Norrbotten", ]$geometry)
    tmp <- myFilteredData.sp$alive %>%  filter(hairTrap)
    plot(tmp$geometry, add = T, col = "red", pch = 16)
    if(length(which(myFilteredData.sp$alive$DNAID[myFilteredData.sp$alive$structured] %in% HairTrapSamples$DNAID)) > 0){
      print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
    }
  }
  
  
  
  ## ------       6.3.3. PLOT CHECKS ------
  
  if(plot.check){
    
    ##-- Barplot of structured vs. opportunistic samples
    pdf(file = file.path(working.dir, "figures/DetectionsStructuredOppBarplot.pdf"))
    par(mfrow = c(2,1), mar = c(4,4,3,2))
    barplot( rbind(table(myFilteredData.sp$alive$Year[myFilteredData.sp$alive$structured]),
                   table(myFilteredData.sp$alive$Year[!myFilteredData.sp$alive$structured])),
             beside = T,
             ylim = c(0,3000),
             col = c(grey(0.2),grey(0.8)),
             ylab = "Number of samples")
    abline(h = seq(0, 3000, by = 500),
           lty = 2, col = grey(0.8))
    title(main = "500m threshold")
    legend("topleft", fill = c(grey(0.2),grey(0.8)),
           legend = c("Structured","Other"))
    
    ##-- Barplot of structured vs. opportunistic samples (threshold = 2000m)
    structured2000 <- myFilteredData.sp$alive$Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
      !is.na(myFilteredData.sp$alive$trackID) &
      myFilteredData.sp$alive$trackDist <= 2000 & 
      !myFilteredData.sp$alive$hairTrap
    barplot( rbind(table(myFilteredData.sp$alive$Year[structured2000]),
                   table(myFilteredData.sp$alive$Year[!structured2000])),
             beside = T,
             ylim = c(0,3000),
             col = c(grey(0.2),grey(0.8)),
             ylab = "Number of samples")
    abline(h=seq(0,3000,by=500),lty=2,col=grey(0.8))
    title(main="2000m threshold")
    legend("topleft",fill=c(grey(0.2),grey(0.8)),
           legend = c("Structured","Other"))
    dev.off()
    
    
    ##-- CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO"
    tmp <- myFilteredData.sp$alive %>%
      filter(Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"))
    
    ## plot check
    pdf(file = file.path(working.dir, "figures/DetectionsStructuredOpp.pdf"))
    for(t in 1:n.years){
      par(mar = c(0,0,3,0), mfrow = c(1,3))
      
      ##-- samples with tracks
      tmpTracks <- tmp %>%
        filter( Year %in% years[t],
                !is.na(trackID),
                trackDist <= 500)
      plot( st_geometry(myStudyArea), col = "gray60", main = "Structured with track")
      plot( st_geometry(tmpTracks),
            pch = 21, col = "black",
            cex = 1, bg = "red", add = T)
      
      ##-- samples without tracks
      tmpNoTracks <- tmp %>%
        filter( Year %in% years[t],
                is.na(trackID) | trackDist > 500)
      plot( st_geometry(myStudyArea), col = "gray60", main = "Structured without track")
      plot( st_geometry(tmpNoTracks),
            pch = 21, col = "black",
            cex = 1, bg = "blue", add = T)
      
      ##-- Other samples
      tmpOpp <- myFilteredData.sp$alive %>%
        filter(!Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),
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
    
    
    ##-- Plot check
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
  
  
  
  ## ------     6.4. SEPARATE MORTALITY CAUSES ------
  
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
  
  
  ##-- Plot check
  if(plot.check){
    par(mfrow = c(1,3))
    for(t in 1:n.years){
      ## DEAD RECOVERIES TOTAL
      tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
      NGS_TabTotal <- table(tempTotal$Country_sf)
      ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country_sf), 2, function(x) sum(x>0))
      ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
      tempIn <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
      NGS_TabIn <- table(tempIn$Country_sf)
      ID_TabIn <- apply(table(tempIn$Id, tempIn$Country_sf), 2, function(x) sum(x>0))
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
    
    temp <- unique(myFilteredData.sp$alive[ ,c("Year","Country_sf","DNAID")])
    tab_Country.Year <- table(temp$Year, temp$Country_sf)
    country.colors <- c("goldenrod1","goldenrod3")
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
    lines(tab_Country.Year[,"(N)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year[,"(S)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("bottomright",c("N","S"), fill=country.colors)
    
    ## ID DETECTED
    temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country_sf,myFilteredData.sp$alive$Id)
    tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
    country.colors <- c("goldenrod1","goldenrod3")
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
    lines(tab_Country.Year1[,"(N)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year1[,"(S)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("bottomright",c("N","S"), fill=country.colors)
    
    ## Average number of detection per detected ID
    tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
         ylab="Average Number of detections", xlab="Years")
    lines(tab_Country.Year2[,"(N)"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year2[,"(S)"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("bottomright",c("N","S"), fill=country.colors)
    
    ## deadrecovery
    temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country_sf","Id")])
    tab_Country.Year <- table(temp$Year, temp$Country_sf)
    country.colors <- c("goldenrod1","goldenrod3")
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
         ylab="N Id Dead recovered", xlab="Years")
    lines(tab_Country.Year[,"(N)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year[,"(S)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("topright",c("N","S"), fill=country.colors)
    dev.off()
  }
  
  
  
  ## ------     6.5. ASSIGN SAMPLES TO DETECTORS -----
  
  ##-- ALL SAMPLES
  myData.alive <- assignDetectors( 
    data = myFilteredData.sp$alive,                
    detectors = detectors$main.detector.sp,
    subDetectors = detectors$detector.sp,
    radius = detectors$resolution)
  
  ##-- DEAD RECOVERY
  myData.dead <- assignDetectors( 
    data = myFilteredData.sp$dead.recovery,
    detectors = detectors$main.detector.sp,
    radius = detectors$resolution)
  
  
  
  ## ------     6.6. FIX DETECTIONS IN NORRBOTTEN ------
  
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
  
  ##-- Loop over flagged detections & assign them to closest detector outside Norrbotten
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
  
  
  
  ## ------     6.7. GENERATE DETECTION HISTORIES : y.alive[i,j,t] & y.dead[i,t] ------
  
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
  
  ##-- Make sure all detection arrays have the same dimensions
  y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- array( 0, 
                                                     dim = dim(y.ar$y.ar),
                                                     dimnames = dimnames(y.ar$y.ar))
  ##-- Fill in the y arrays
  y.ar.ALIVEOthers[dimnames(y.arOth$y.ar)[[1]], , ] <- y.arOth$y.ar
  y.ar.ALIVEStructured[dimnames(y.arStruc$y.ar)[[1]], , ] <- y.arStruc$y.ar
  
  ##-- Project death to the next occasion.
  y.ar.DEADProjected <- y.ar$y.ar2 
  y.ar.DEADProjected[] <- 0
  for(t in 2:n.years){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}
  
  ##-- Create binary dead recovery histories (0: not recovered ; 1: recovered dead)
  y.ar.DEAD <- apply(y.ar.DEADProjected, c(1,3), function(x){as.numeric(sum(x)>0)})
  dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
  y.ar.DEAD[y.ar.DEAD > 0] <- 1
  
  
  
  
  ## ------     6.8. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------
  
  distances <- list()
  for(t in 1:n.years){
    
    ##[PD] WE NEED TO DISCUSS THE MAX DIST USED HERE 
    ## MUCH SMALLER THAN THE MAX DIST USED IN THE LOCAL EVAL
    ##-- Identify detections further than maxDist
    print(paste("------ ", t ," -------", sep = "" ))
    distances[[t]] <- checkDistanceDetections( 
      y = y.ar$y.ar[,,t], 
      detector.xy = detectors$detectors.df[, c("x","y")], 
      max.distance = 40000, #maxDist
      method = "pairwise",
      plot.check = F)
    
    ##-- Plot individuals that have detections further than the threshold
    if(plot.check){
      par(mfrow = c(1,1))
      if(sum(distances[[t]]$y.flagged) > 0){
        affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
        count <- 0
        for(i in affected.ids){
          count <- count+1
          plot(st_geometry(myStudyArea), main = paste0("t: ",t,"     i: ", names(affected.ids)[count]))
          scalebar( 2 * detectors$maxDist,
                    xy = c(800000,6700000),
                    type = "bar", divs = 2, below = "km",
                    label = c(0, detectors$maxDist/1000, detectors$maxDist/500),
                    cex = 0.8, adj = c(0.5,-0.9))
          plot( st_geometry(COUNTRIES), add = T)
          plot( st_geometry(detectors$main.detector.sp),
                add = T, col = grey(0.8), cex = 0.3, pch = 19)
          
          tmp <- myData.alive$data.sp[myData.alive$data.sp$Id == dimnames(y.ar$y.ar)[[1]][i] &
                                        myData.alive$data.sp$Year == years[t], ]
          tmp <- tmp[order(tmp$Date), ]
          tmp.xy <- st_coordinates(tmp)
          n.det <- nrow(tmp.xy)
          
          plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
          arrows( x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
                  x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
                  length = 0.1, lwd = 1)
          plot(st_geometry(detectors$main.detector.sp[which(y.ar$y.ar[i, ,t] > 0), ]),
               pch = 16, col = "red", add = T)
          
          tmp2 <- detectors$main.detector.sp[which(y.ar$y.ar[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
          plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
        }#i
      }#if
    }#if
    
    ##-- Remove detections that are further then the threshold
    #y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
    
    ##-- Remove detections also in myData.alive to run getSInits later
    idd <- names(affected.ids)
    for(i in 1:length(idd)){
      detIds <- which(distances[[t]]$y.flagged[idd[i], ]>0)
      # myData.alive$data.sp <- myData.alive$data.sp[!(myData.alive$data.sp$Id %in% idd[i] &
      #                                                      myData.alive$data.sp$Detector %in% detIds &
      #                                                      myData.alive$data.sp$Year %in% years[t]), ]
      myData.alive$data.sp <- myData.alive$data.sp %>%
        dplyr::filter(!(Id %in% idd[i] & Detector %in% detIds & Year %in% years[t]))
    }#i
  }#t
  
  
  
  ## ------   7. GENERATE INDIVIDUAL-LEVEL COVARIATES ------
  
  ## ------     7.1. TRAP-RESPONSE ------
  
  ##-- Make matrix of previous capture indicator
  already.detected <- makeTrapResponseCov(
    data = myFullData.sp$alive,
    data.dead = myFullData.sp$dead.recovery)
  
  ##-- Subset to focal years
  already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar$y.ar)[[3]]]
  
  ##-- Subset to focal individuals
  already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar$y.ar)[[1]], ]
  
  ##-- Set first detection for augmented individuals to NA
  already.detected[rownames(already.detected) %in% "Augmented",1]  <- NA
  
  ##-- Plot an image of the matrix
  if(plot.check){
    par(mfrow = c(1,1))
    barplot(colSums(apply(y.ar$y.ar, c(1,3), function(x) any(x>0))))
    barplot(colSums(already.detected), add = TRUE, col = "gray40")
    legend( x = 0, y = 250, 
            legend = c("newly Det", "already Det"),
            fill = c("gray80", "gray40"))
  }
  
  
  
  ## ------     7.2. AGE ------
  
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
  
  
  
  
  
  
  
  ## ------   8. MAKE AUGMENTATION ------
  
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
  
  
  
  ## ------   9. TRANSFORM Y TO SPARSE MATRICES ------
  
  ##-- STRUCTURED
  y.sparse <- nimbleSCR::getSparseY(y.aliveStructured)
  
  ##-- OTHER
  y.sparseOth <- nimbleSCR::getSparseY(y.aliveOthers)
  
  
  
  ## ------ III. MODEL SETTING & RUNNING ------- 
  
  ## ------   1. NIMBLE MODEL DEFINITION ------
  
  modelCode <- nimbleCode({
    
    ##------ SPATIAL PROCESS ------## 
    
    dmean ~ dunif(0,100)
    lambda <- 1/dmean
    
    betaDens ~ dnorm(0.0,0.01)
    habIntensity[1:n.habWindows] <- exp(betaDens * denCounts[1:n.habWindows])
    sumHabIntensity <- sum(habIntensity[1:n.habWindows])
    logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
    logSumHabIntensity <- log(sumHabIntensity)
    
    for(i in 1:n.individuals){
      sxy[i,1:2,1] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
        upperCoords = upperHabCoords[1:n.habWindows,1:2],
        logIntensities = logHabIntensity[1:n.habWindows],
        logSumIntensity = logSumHabIntensity,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
      
      for(t in 2:n.years){
        sxy[i,1:2,t] ~ dbernppACmovement_exp(
          lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
          upperCoords = upperHabCoords[1:n.habWindows,1:2],
          s = sxy[i,1:2,t-1],
          lambda = lambda,
          baseIntensities = habIntensity[1:n.habWindows],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max,
          numWindows = n.habWindows)
      }#i  
    }#t
    
    
    
    ##----- DEMOGRAPHIC PROCESS -----## 
    
    omeg1[1:2] ~ ddirch(alpha[1:2])   
    
    for(t in 1:(n.years-1)){
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
      for(t in 1:(n.years-1)){
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
      for(c in 1:n.covs.Oth){
        betaCovsOth[c,t] ~ dunif(-5,5)
      }#c 
      for(c in 1:n.countries){
        p01Oth[c,t] ~ dunif(0,1)
        p0Oth[c,t] <- p01Oth[c,t] * countryToggle[c,t]## toggle countries
      }#c  
    }#t
    
    
    ## Individual response
    pResponse ~ dunif(0, 1)
    for(i in 1:n.individuals){ 
      detResponse[i,1] ~ dbern(pResponse)
    }#i
    
    ## Individual detection histories
    for(t in 1:n.years){
      for(i in 1:n.individuals){
        
        y[i,1:maxDetNums,t] ~ dbinomLocal_normalCovsResponse2( 
          detNums = detNums[i,t],
          detIndices = detIndices[i,1:maxDetNums,t],
          size = size[1:n.detectors],
          p0State = p0[1:n.counties,t],
          sigma = sigma[t],
          s = sxy[i,1:2,t],
          trapCoords = detector.xy[1:n.detectors,1:2],
          localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
          localTrapsNum = localDetNum[1:n.habWindows],
          resizeFactor = resizeFactor,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          indicator = isAlive[i,t],
          lengthYCombined = lengthYCombined,
          trapCountries = detCounties[1:n.detectors],
          trapCovs = detCovs[1:n.detectors,t,1:n.covs],
          trapBetas = betaCovs[1:n.covs,t],
          responseCovs = detResponse[i,t],
          responseBetas = betaResponse[t])
        
        y.Oth[i,1:maxDetNumsOth,t] ~ dbinomLocal_normalCovsResponse2( 
          detNums = detNumsOth[i,t],
          detIndices = detIndicesOth[i,1:maxDetNumsOth,t],
          size = size[1:n.detectors],
          p0State = p0Oth[1:n.counties,t],
          sigma = sigma[t],
          s = sxy[i,1:2,t],
          trapCoords = detector.xy[1:n.detectors,1:2],
          localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
          localTrapsNum = localDetNum[1:n.habWindows],
          resizeFactor = resizeFactor,
          lengthYCombined = lengthYCombined.Oth,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          indicator = isAlive[i,t],
          trapCountries = detCountries[1:n.detectors,t],
          trapCovs = detCovsOth[1:n.detectors,t,1:n.covs.Oth],
          trapBetas = betaCovsOth[1:n.covs.Oth,t],
          responseCovs = detResponse[i,t],
          responseBetas = betaResponseOth[t])
        
        # y.dead.legal[i,t] ~ dbern(z[i,t] == 3) 
        # y.dead.other[i,t] ~ dbern(z[i,t] == 4) 
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
    n.individuals = dim(y.sparse$y)[1],
    n.detectors = nrow(detectors$scaledCoords),
    n.habWindows = nrow(habitat$scaledLowerCoords),
    n.years =  dim(y.sparse$y)[3], 
    n.covs = dim(detCovs)[3],
    n.covs.Oth = dim(detCovsOth)[3],
    n.countries = max(detCountries),
    n.counties = max(detCounties),
    countyToggle = countyToggle,
    countryToggle = countryToggle,
    resizeFactor = detectors$localObjects$resizeFactor,
    y.max = dim(detectors$localObjects$habitatGrid)[1],
    x.max = dim(detectors$localObjects$habitatGrid)[2],
    numLocalIndicesMax = detectors$localObjects$numLocalIndicesMax,
    maxDetNums = y.sparse$maxDetNums,
    maxDetNumsOth = y.sparseOth$maxDetNums,
    lengthYCombined = y.sparse$lengthYCombined,
    lengthYCombined.Oth = y.sparseOth$lengthYCombined)
  
  
  
  ## ------   3. NIMBLE DATA ------
  
  ## ------     3.1. GENERATE KNOWN z ------
  
  ##-- Set all individuals alive to 2 between first and last detection
  z <- apply(y.alive, c(1,3), function(x)ifelse(any(x>0), 2, NA))
  z <- t(apply(z, 1, function(zz){
    if(any(!is.na(zz))){
      range.det <- range(which(!is.na(zz)))
      zz[range.det[1]:range.det[2]] <- 2
    }
    return(zz)
  }))
  
  
  
  ## ------     3.2. LIST DATA ------
  
  nimData <- list( 
    z = z,   
    y = y.sparse$y,
    detIndices = y.sparse$detIndices,
    detNums = y.sparse$detNums,
    y.Oth = y.sparseOth$y, 
    detIndicesOth = y.sparseOth$detIndices,
    detNumsOth = y.sparseOth$detNums,
    lowerHabCoords = as.matrix(habitat$scaledLowerCoords), 
    upperHabCoords = as.matrix(habitat$scaledUpperCoords), 
    detCounties = detCounties,
    detCountries = detCountries,
    detCovs = detCovs,
    detCovsOth = detCovsOth,
    detResponse = already.detected,
    denCounts = denCounts,
    localDetIndices = detectors$localObjects$localIndices,
    localDetNum = detectors$localObjects$numLocalIndices,
    habitatGrid = detectors$localObjects$habitatGrid,
    size = detectors$detectors.df$size,
    alpha = rep(1,2),
    detector.xy = as.matrix(detectors$scaledCoords))
  #habitatGrid = habIDCells.mx)
  
  
  
  ## ------   4. NIMBLE INITS ------
  
  ## ------     4.1. GENERATE INITIAL z ------
  
  ##-- Set z to 1 before first detection and 3 after last detection
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
  
  ##-- Set initial values to NA when individual state is known
  z.init <- ifelse(!is.na(z), NA, z.init)
  
  
  
  ## ------     4.2. LATENT VARIABLE DET RESPONSE ------
  
  detResponse.inits <- nimData$detResponse
  detResponse.inits[is.na(detResponse.inits)] <- rbinom(sum(is.na(detResponse.inits)),1,0.5)
  detResponse.inits[!is.na(already.detected)] <- NA
  
  
  
  ## ------     4.3. GENERATE INITIAL sxy ------
  
  ##-- sxy
  AllDets <- myData.dead[ ,c("Id","Year")] %>% 
    ##-- Project death to the next year
    mutate(Year = Year + 1) %>%
    ##-- Remove dead reco occuring the last year (not used)
    filter(!Year %in% max(Year)) %>%
    ##-- Combine with detections alive
    rbind(.,myData.alive$data.sp[ ,c("Id","Year")]) %>%
    ##-- Add coordinates
    mutate("x" = st_coordinates(.)[ ,1],
           "y" = st_coordinates(.)[ ,2]) %>%
    as.data.frame()
  
  ##-- Rescale detections
  AllDetsxyscaled <- scaleCoordsToHabitatGrid(
    coordsData = AllDets,
    coordsHabitatGridCenter = habitat$habitat.xy,
    scaleToGrid =T )$coordsDataScaled
  
  
  # ##-- Create a data.frame with all detection of all Individuals detected
  # ##-- Project death to the next year
  # myData.deadProj <- myData.dead[ ,c("Id","Year")]
  # myData.deadProj$Year <- myData.deadProj$Year + 1 
  # 
  # ##-- Remove dead reco occuring the last year (not used)
  # myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]
  # 
  # AllDets <- rbind( myData.deadProj,
  #                   myData.alive$myData.sp[ ,c("Id","Year")])
  # AllDetections <- as.data.frame(AllDets)
  # AllDetsxy <- st_coordinates(AllDets) 
  # colnames(AllDetsxy) <- c("x","y")
  # AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  #   coordsData = AllDetsxy,
  #   coordsHabitatGridCenter = habitat$habitat.xy,
  #   scaleToGrid = T)$coordsDataScaled
  #
  # AllDetections <- cbind(AllDetections, AllDetsxyscaled)
  
  sxy.init <- getSInits( AllDetections = AllDetsxyscaled[,c("Id","Year","x","y")],
                         Id.vector = y.ar$Id.vector,
                         idAugmented = which(rownames(z) %in% "Augmented"),
                         lowerCoords = nimData$lowerHabCoords,
                         upperCoords = nimData$upperHabCoords,
                         habitatGrid = nimData$habitatGrid,
                         intensity = NULL,
                         sd = 4,
                         movementMethod = "dbernppACmovement_normal")
  
  ##-- An extreme number of decimals may cause a number to appear as an integer
  ##-- to Nimble, and then coincide with habitat window boundaries
  sxy.init <- round(sxy.init, 4)
  
  # ##-- Get location of individuals 
  # sxy.initscaled <- scaleCoordsToHabitatGrid(
  #   coordsData = sxy.init,
  #   coordsHabitatGridCenter = habitat$habitat.xy,
  #   scaleToGrid =F )$coordsDataScaled
  
  
  
  ## ------   5. NIMBLE PARAMETERS ------
  
  nimParams <- c( "N", "lambda", "dmean", "betaDens",
                  "omeg1", "gamma", "phi",
                  "pResponse", "sigma",
                  "p0", "betaResponse", "betaCovs",
                  "p0Oth", "betaResponseOth", "betaCovsOth")
  
  nimParams2 <- c("z", "sxy")
  
  
  
  ## ------   6. LOOP THROUGH INITIAL VALUES & SAVE OBJECT ------
  
  for(c in 1:4){
    
    ## ------    6.1. LIST INITIAL VALUES ------
    nimInits <- list(
      "sxy" = sxy.init,
      "z" = z.init,
      "dmean" = runif(1, 0, 10),
      "betaDens" = runif(1, -0.1, 0.1),
      "omeg1" = c(0.5, 0.5),
      "gamma" = runif(dim(y.alive)[3]-1, 0, 1),
      "phi" = runif(dim(y.alive)[3]-1, 0.1, 0.3),
      "pResponse" = runif(1, 0.4, 0.5),
      "detResponse" = detResponse.inits,
      "sigma" = runif(n.years, 1, 4),
      "p01" = array(runif(18, 0, 0.2),
                    c(nimConstants$n.counties, dim(y.alive)[3])),
      "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
      "betaCovs" = array(runif(dim(detCovs)[3], -0.1, 0.1),
                         c(dim(detCovsOth)[3], n.years)),
      "p01Oth" = array(runif(18, 0, 0.2),
                       c(nimConstants$n.countries+1, dim(y.alive)[3])),
      "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
      "betaCovsOth" = array(runif(dim(detCovsOth)[3], -0.1, 0.1),
                            c(dim(detCovsOth)[3], n.years))) 
    
    # ##-- Test if the local evaluation wil work 
    # ##-- Get detector index from the habitat ID matrix
    # idDEtected <- which(!rownames(z) %in%"Augmented")
    # for(i in 1:length(idDEtected)){
    #   for(t in 1:nimConstants$n.years){
    #     if(!is.na(nimInits$sxy[i,1,t])){
    #       SXY <- nimInits$sxy[i,,t]  
    #     }else{SXY <- nimData$sxy[i,,t]}
    #     sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
    #     DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
    #                                                       1:nimConstants$maxNBDets] 
    #     DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
    #     index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
    #     
    #     ## GET NECESSARY INFO 
    #     n.detectors <- length(index)
    #     YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
    #     YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
    #     
    #     ## RECREATE Y
    #     if(nimData$nbDetections[i, t] > 0){
    #       for(j in 1:nimData$nbDetections[i, t]){
    #         ## check if a detection is out of the "detection window"
    #         if(sum(YDET[j]==index)==0){
    #           print(paste("id",i,"t",t,"j",j))
    #         }
    #       }
    #     }
    #   }}

    
    
    ## ------    6.2. SAVE NIMBLE INPUT ------
    
    save(nimData,
         nimConstants,
         y.dead,
         nimParams,
         nimParams2,
         modelCode,
         nimInits,
         file = file.path( working.dir, "nimbleInFiles", thisSex,
                           paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
  }#c
}#thisSex



##------    7. SAVE NECESSARY OBJECTS ------

myHabitat.list <- habitat
myDetectors <- detectors
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
#   colnames(recruitnb) <-  colnames(recruit)  <- 
#     factorValues(habitatRasterResolution$'2.5km'$Countries,lev[[1]]$ID)[,1]
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
# phi <- phi[,c(2,4)]                    # overall phi
# recruit <- recruit[,c(2,4)]            # overall phi
# recruitnb <- recruitnb[,c(2,4)]        # overall phi
# 
# 
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
# lev <- levels(habitatRasterResolution$'10km'$Counties)
# countryId <- list()
# for(t in 1:dim(z)[2]){
#   tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
#   countryId[[t]] <- raster::extract( habitatRasterResolution$'10km'$Counties ,tmp,sparse = F)
# }
# 
# phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=n.years-1,ncol=length(lev[[1]]$ID))
# colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
#   colnames(recruitnb) <-  colnames(recruit)  <- 
#    factorValues(habitatRasterResolution$'10km'$Counties,lev[[1]]$ID)[,1]
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
# phi <- phi[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]               # overall phi
# phiNOR <- phi[,c(8,9,10,11,12,13)]
# phiSWE <- phi[,-c(8,9,10,11,12,13,14,15,16)]
# recruitnb <- recruitnb[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]   # overall phi
# recruitnbNOR <- recruitnb[,c(8,9,10,11,12,13)]
# recruitnbSWE <- recruitnb[,-c(8,9,10,11,12,13,14,15,16)]
# recruit <- recruit[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]       # overall phi
# recruitNOR <- recruit[,c(8,9,10,11,12,13)]
# recruitSWE <- recruit[,-c(8,9,10,11,12,13,14,15,16)]
# 
# 
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




##------------------------------------------------------------------------------


## ------ III. NIMBLE RUN ====

## ------   1. CONFIGURE NIMBLE MODEL ====

load(file.path( working.dir, "nimbleInFiles/male/nimbleInput_2025-11-13_male_1.RData"))
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 
model$calculate()

t <- 1
i <- 1

dbinomLocal_normalCovsResponse( 
  x = nimData$y[i, ,t],
  s = nimInits$sxy[i, ,t],
  sigma = nimInits$sigma[t],
  detNums = nimData$detNums[i,t],
  detIndices = nimData$detIndices[i, ,t],
  trapCoords = nimData$detector.xy,
  size = nimData$size,
  localTrapsIndices = ni$,
  localTrapsNum = nimData$localDetNum
  ResizeFactor = nimConstants$ResizeFactor,
  maxNBDets = nimConstants$maxNBDets,
  habitatID = nimData$habitatIDDet[1:nimConstants$y.maxDet,1:nimConstants$x.maxDet],
  indicator = model$z[i,t]    ,
  p0State = model$p01[1:nimConstants$n.countries,t],
  detCountries = nimData$detCountries[1:nimConstants$n.detectors],
  detCov = nimData$detCovs[1:nimConstants$n.detectors,t,1:nimConstants$n.covs],
  betaCov = nimInits$betaCovs[1:nimConstants$n.covs],
  BetaResponse = nimInits$betaResponse[t],
  detResponse = nimData$detResponse[i,t])




## ------   2. CHECK INITIAL LOG-LIKELIHOODS ====

model$calculate("sxy")
which(is.infinite(model$logProb_sxy), arr.ind = T)
which(is.infinite(model$logProb_y.alive), arr.ind = T)
which(is.infinite(model$logProb_y.aliveOth), arr.ind = T)
which(is.infinite(model$logProb_z), arr.ind = T)

whichlogprop <- names(model)[grep("logProb_",names(model))]
whichlogprop <- whichlogprop[-grep("env_",whichlogprop)]
whichlogprop <- whichlogprop[-grep("row",whichlogprop)]
whichlogprop <- whichlogprop[-grep("name",whichlogprop)]

for(i in 1:length(whichlogprop)){
  tmp <- model[[whichlogprop[i]]] #[[grep("logProb_",names(model))[i]]]
  if(sum(is.infinite(tmp))>0){
    print(whichlogprop[i])
  }
}

model$y.alive[1377, ,9]
model$yDets[1377, ,9]
model$detector.xy[model$yDets[1377,1,9],]
model$sxy[1377, ,9]
model$z[1377, ]
model$detCovs[model$yDets[1377,1,9], ,]
model$detCountries[model$yDets[1377,1,9]]

plot(model$detector.xy[,2]~model$detector.xy[ ,1])
points(model$detector.xy[model$detCountries == 5,2] ~ model$detector.xy[model$detCountries == 5,1], col = "red")
model$trials

which(is.na(model$logProb_sxy),arr.ind = T)
which(is.na(model$logProb_y.alive),arr.ind = T)
which(is.na(model$logProb_y.aliveOth),arr.ind = T)

model$calculate("y.alive")

for( t in 1:3){
  for(i in 1: nimConstants$n.individuals){
    dbin_LESS_Cached_MultipleCovResponse( x=nimData$y.alive[i,1:nimConstants$nMaxDetectors,t]
                                          ,
                                          sxy = nimInits$sxy[i,1:2,t]
                                          ,
                                          
                                          sigma = nimInits$sigma[t]
                                          ,
                                          nbDetections= nimData$nbDetections[i,t]
                                          ,
                                          yDets = nimData$yDets[i,1:nimConstants$nMaxDetectors,t]
                                          ,
                                          detector.xy =  nimData$detector.xy[1:nimConstants$n.detectors,1:2]
                                          ,
                                          trials = nimData$trials[1:nimConstants$n.detectors]
                                          ,
                                          
                                          detectorIndex = nimData$detectorIndex[1:nimConstants$n.cellsSparse,1:nimConstants$maxNBDets]
                                          ,
                                          nDetectorsLESS = nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
                                          ,
                                          ResizeFactor = nimConstants$ResizeFactor
                                          ,
                                          maxNBDets = nimConstants$maxNBDets
                                          ,
                                          habitatID = nimData$habitatIDDet[1:nimConstants$y.maxDet,1:nimConstants$x.maxDet]
                                          ,
                                          indicator = model$z[i,t]    
                                          ,
                                          p0State = model$p01[1:nimConstants$n.countries,t]
                                          ,
                                          detCountries = nimData$detCountries[1:nimConstants$n.detectors]
                                          ,
                                          detCov = nimData$detCovs[1:nimConstants$n.detectors,t,1:nimConstants$n.covs]
                                          ,
                                          betaCov = nimInits$betaCovs[1:nimConstants$n.covs]
                                          ,
                                          BetaResponse = nimInits$betaResponse[t]
                                          ,
                                          detResponse = nimData$detResponse[i,t])
  }
}
dbin_LESS_Cached_MultipleCovResponse




cmodel$calculate()  
cmodel$calculate("y.alive")
which(cmodel$logProb_y.alive == -Inf, arr.ind = TRUE)

cmodel$calculate("sxy")
which(cmodel$sxy == -Inf, arr.ind = TRUE)
which(is.infinite(cmodel$sxy))
sum(is.infinite(cmodel$sxy))
cmodel$calculate("dispSigma")
cmodel$calculate("p0")
cmodel$calculate("mu")
cmodel$calculate("dispSigma")
cmodel$calculate("betaDens")

cmodel$calculate("z")

if(cmodel$calculate("y.alive") == -Inf){
  probs <- which(cmodel$logProb_y.alive == -Inf, arr.ind = TRUE)
  print(dim(probs)[1])
  if(dim(probs)[1] < 11){
    for(d in 1:dim(probs)[1]){
      #   for(d in sample(dim(probs)[1], 10)){
      
      ## RETRIEVE ID & YEAR
      i <- probs[d,1]
      t <- probs[d,3]
      
      ## PLOT HABITAT SET-UP
      plot(myHabitat$habitat.r)
      plot(myHabitat$buffered.habitat.poly, add = TRUE)
      plot(myStudyArea, add = TRUE)
      plot(myDetectors$main.detector.sp, cex = 0.3, col = "gray20", pch = 16, add = TRUE)
      
      ## PLOT SIMULATED ACs & DETECTIONS
      plot(myDetectors$main.detector.sp[which(y.alive[i, ,t]>0), ], pch = 16, col = "navyblue", add = TRUE)
      
      ## RETRIEVE NECESSARY INFOS
      x = cmodel$y.alive[i, ,t]       
      detectionsNum = cmodel$nbDetections[i,t]
      detectionsID = nimData$yDets[i,1:detectionsNum,t]
      pZero = cmodel$p0[t]
      sigma = cmodel$sigma
      sxy = cmodel$sxy[i,1:2,t]
      detectorCoords = nimData$detCoords
      detectorID = nimData$detID
      detectorNum = nimData$detNum
      detectorTrials = nimData$trials
      habitatID = nimData$habID
      habitatFactor = nimConstants$habFactor
      habitatMinX = nimConstants$habMinX
      habitatMaxY = nimConstants$habMaxY
      habitatResolution = nimConstants$habRes
      indicator = cmodel$z[i,t]
      n.detectors <- nimConstants$numDetectors
      
      print(paste("z ==", indicator))
      
      
      dbin_Cached_Sparse(x,
                         detectionsNum,
                         detectionsID,
                         pZero,
                         sigma,
                         sxy,
                         detectorCoords,
                         detectorID,
                         detectorNum,
                         detectorTrials,
                         habitatID,
                         habitatFactor,
                         habitatMinX,
                         habitatMaxY,
                         habitatResolution,
                         indicator,
                         log =  0)
      
      ## GET HABITAT CELL ID FROM THE HABITAT ID MATRIX
      scaledX <- (sxy[1] - habitatMinX) / (habitatResolution * habitatFactor)
      scaledY <- -(sxy[2] - habitatMaxY) / (habitatResolution * habitatFactor)
      sID <- habitatID[trunc(scaledY) + 1, trunc(scaledX) + 1]
      
      ## GET NUMBER OF DETECTORS WITHIN maxDist FROM THE DETECTOR NUMBER MATRIX
      detNum <- detectorNum[sID]
      
      ## GET IDs OF DETECTORS WITHIN maxDist OF THE HABITAT CELL FROM THE DETECTOR ID MATRIX
      detIDs <- detectorID[sID, 1:detNum]
      
      ## PLOT MODEL sxy & ASSOCIATED DETECTORS
      points(sxy[1], sxy[2], pch = 3, col = "red")
      points(myDetectors$main.detector.sp[detIDs,1], myDetectors$main.detector.sp[detIDs,2])
      
      ## PLOT DETECTIONS OUTSIDE THE ALLOWED DETECTORS
      outDets <- detectionsID[which(!detectionsID[1:detectionsNum] %in% detIDs)]
      points(myDetectors$main.detector.sp[outDets,1], myDetectors$main.detector.sp[outDets,2], pch = 16, col = "red")
    }#i
  } else { print("TOO MANY PROBLEMATIC INDIVIDUALS TO BE DISPLAYED (>10)") }
}

if(is.na(cmodel$calculate("y.alive"))){
  
  probs <- which(is.na(cmodel$logProb_y.alive), arr.ind = TRUE)
  print(dim(probs)[1])
  
  if(dim(probs)[1] < 11){
    ## for(d in 1:dim(probs)[1]){
    for(d in sample(dim(probs)[1], 10)){
      
      ## RETRIEVE ID & YEAR
      i <- probs[d,1]
      t <- probs[d,3]
      
      ## PLOT HABITAT SET-UP
      plot( myHabitat$habitat.r)
      plot( myHabitat$buffered.habitat.poly, add = TRUE)
      plot( myStudyArea, add = TRUE)
      plot( myDetectors$main.detector.sp,
            cex = 0.3, col = "gray20", pch = 16, add = TRUE)
      
      ## PLOT SIMULATED ACs & DETECTIONS
      plot( myDetectors$main.detector.sp[which(y.alive[i, ,t]>0), ],
            pch = 16, col = "navyblue", add = TRUE)
      
      dbin_LESS_Cached_OneCov( 
        x = cmodel$y.alive[i, ,t],
        sxy = cmodel$sxy[i, ,t],
        sigma = cmodel$sigma,
        nbDetections = cmodel$nbDetections[i,t],
        yDets = cmodel$yDets[i, ,t],
        detector.xy = cmodel$detector.xy,
        trials = cmodel$trials,
        detectorIndex = cmodel$detectorIndex,
        nDetectorsLESS = cmodel$nDetectorsLESS,
        ResizeFactor = nimConstants$ResizeFactor,
        maxNBDets = nimConstants$maxNBDets,
        habitatID = cmodel$habitatIDDet,
        indicator = cmodel$isAlive[i,t],
        p0State = cmodel$p0[ ,2,t],
        detCountries = cmodel$detCountries,
        detCov = cmodel$detTracks[ ,t],
        betaCov = cmodel$betaTracks,
        log = 0)
      
      ## GET HABITAT CELL ID FROM THE HABITAT ID MATRIX
      scaledX <- (sxy[1] - habitatMinX) / (habitatResolution * habitatFactor)
      scaledY <- -(sxy[2] - habitatMaxY) / (habitatResolution * habitatFactor)
      sID <- habitatID[trunc(scaledY) + 1, trunc(scaledX) + 1]
      
      ## GET NUMBER OF DETECTORS WITHIN maxDist FROM THE DETECTOR NUMBER MATRIX
      detNum <- detectorNum[sID]
      
      ## GET IDs OF DETECTORS WITHIN maxDist OF THE HABITAT CELL FROM THE DETECTOR ID MATRIX
      detIDs <- detectorID[sID, 1:detNum]
      
      ## PLOT MODEL sxy & ASSOCIATED DETECTORS
      points(sxy[1], sxy[2], pch = 3, col = "red")
      points(myDetectors$main.detector.sp[detIDs,1], myDetectors$main.detector.sp[detIDs,2])
      
      ## PLOT DETECTIONS OUTSIDE THE ALLOWED DETECTORS
      outDets <- detectionsID[which(!detectionsID[1:detectionsNum] %in% detIDs)]
      points( myDetectors$main.detector.sp[outDets,1], myDetectors$main.detector.sp[outDets,2],
              pch = 16, col = "red")
    }#i
  } else { print("TOO MANY PROBLEMATIC INDIVIDUALS TO BE DISPLAYED (>10)") }
}

if(cmodel$calculate("sxy") == -Inf){
  probs <- which(cmodel$logProb_sxy == -Inf, arr.ind = TRUE)
  print(dim(probs)[1])
}



## ------   3. RESUME MODEL CONFIGURATION ====

conf <- configureMCMC(model, monitors = nimParams, thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc



## ------   4. RUN NIMBLE MCMC IN SUCCESSIVE BITES ====

## SET NUMBER OF BITES AND NUMBER OF ITERATIONS PER BITE
bite.size <- 100 
bite.number <- 10

## LOOP OVER NUMBER OF BITES
for(nb in 1:bite.number){
  print(nb)
  if(nb == 1){
    ## run initial MCMC
    MCMCRuntime <- system.time(Cmcmc$run(bite.size))
  } else {      
    ## run subsequent MCMCs
    MCMCRuntime <- system.time(Cmcmc$run(bite.size, reset = FALSE))
  }
  
  ## STORE BITE OUTPUT IN A MATRIX
  mcmcSamples <- as.matrix(Cmcmc$mvSamples)
  CumulRunTime <- proc.time() - ptm
  
  ## EXPORT NIMBLE OUTPUT 
  outname <- file.path(path.OUT, paste("NimbleBite", nb, "_FOR", set, sep = ""))
  save(CumulRuntime, MCMCRuntime, mcmcSamples, file = outname)
  
  ## FREE UP MEMORY SPACE 
  rm("mcmcSamples") 
  Cmcmc$mvSamples$resize(0) ## reduce the internal mvSamples object to 0 rows,
  gc() ## run R's garbage collector
}#nb



##------------------------------------------------------------------------------
## ------ IV. TURN OPSCR INTO SCR ------

## ------   1. SCR MODEL ====

modelCode1 <- nimbleCode({
  
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  betaDens  ~ dnorm(0.0,0.01)
  
  habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows =  y.max,
      numGridCols = x.max)
  }#i
  
  
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  
  pResponse ~ dunif(0, 1)
  psi ~ dunif(0, 1)
  for(i in 1:n.individuals){ 
    detResponse[i] ~ dbern(pResponse)
    z[i] ~ dbern(psi)	
  }#t 
  
  
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  sigma ~ dunif(0,4)
  for(c in 1:n.covs){
    betaCovs[c] ~ dunif(-5,5)
  }
  
  for(c in 1:n.covsOth){
    betaCovsOth[c] ~ dunif(-5,5)
  }
  
  betaResponse ~ dunif(-5,5)
  betaResponseOth ~ dunif(-5,5)
  
  for(c in 1:n.counties){
    p01[c] ~ dunif(0,1)
    p0[c] <- p01[c] *countyToggle[c]## toggle counties
  }#c  
  
  for(c in 1:n.countries){
    p01Oth[c] ~ dunif(0,1)
    p0Oth[c] <- p01Oth[c] *countyToggleOth[c]## toggle countries
  }#c  
  
  for(i in 1:n.individuals){
    y.alive[i,1:nMaxDetectors] ~ dbin_LESS_Cached_MultipleCovResponse( 
      sxy = sxy[i,1:2],
      sigma = sigma,
      nbDetections[i],
      yDets = yDets[i,1:nMaxDetectors],
      detector.xy =  detector.xy[1:n.detectors,1:2],
      trials = trials[1:n.detectors],
      detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
      nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
      ResizeFactor = ResizeFactor,
      maxNBDets = maxNBDets,
      habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      p0[1:n.counties],
      detCounties[1:n.detectors],
      detCov = detCovs[1:n.detectors,1:n.covs],
      betaCov = betaCovs[1:n.covs],
      BetaResponse = betaResponse,
      detResponse = detResponse[i])
    
    y.aliveOth[i,1:nMaxDetectorsOth] ~ dbin_LESS_Cached_MultipleCovResponse(
      sxy = sxy[i,1:2],
      sigma = sigma,
      nbDetectionsOth[i],
      yDets = yDetsOth[i,1:nMaxDetectorsOth],
      detector.xy =  detector.xy[1:n.detectors,1:2],
      trials = trials[1:n.detectors],
      detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
      nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
      ResizeFactor = ResizeFactor,
      maxNBDets = maxNBDets,
      habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      p0Oth[1:n.countries],
      detCountries[1:n.detectors],
      detCov = detCovsOth[1:n.detectors,1:n.covsOth],
      betaCov = betaCovsOth[1:n.covsOth],
      BetaResponse = betaResponseOth,
      detResponse = detResponse[i])
  }#i
  
  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  
  N <- sum(z[1:n.individuals])
  
})



## ------   2. LOOP OVER YEARS ====

for(thisSex in c("Hann","Hunn")){
  
  ## CREATE DIRECTORY TO STORE SCR INPUT
  dir.create(file.path(myVars$WD, myVars$modelName, "SCR", thisSex))
  
  for(ch in 1:4){
    for(t in 1:nYears){   
      
      ## ------   2.1. LOAD OPSCR INPUT FILE ====
      
      load( file.path(myVars$WD, myVars$modelName, thisSex,
                      paste0(myVars$modelName, thisSex, "_Chain", c, ".RData")))
      
      ## GET DETECTED INDIVIDUALS       
      detectedStruc <- apply(nimData$nbDetections, 2, function(x) x>0)
      detectedOth <- apply(nimData$nbDetectionsOth, 2, function(x) x>0)
      detected <- detectedOth + detectedStruc
      detected <- detected > 0
      
      ## GET SUM OF INDIVIDUALS DETECTED & DECIDE HOW MUCH YOU WISH TO AUGMENT. 
      ## HERE I CHOSE 3 
      N.det <- sum(detected[ ,t])
      N.aug <- N.det * 3                   ## decide the augmentation factor 
      
      
      
      ## ------   2.2. SUBSET NIMDATA ====
      
      ## SUBSET & AUGMENT ALL OBJECTS BASED ON WHETHER INDIVIDUALS WERE DETECTED OR NOT
      ## y.alive
      nimData$y.alive <- nimData$y.alive[detected[ ,t], ,t]  
      nimData$y.alive <- rbind( nimData$y.alive,
                                matrix( 0,
                                        nrow = N.aug,
                                        ncol = nimConstants$nMaxDetectors))
      
      nimData$y.aliveOth <- nimData$y.aliveOth[detected[ ,t], ,t]  
      nimData$y.aliveOth <- rbind( nimData$y.aliveOth,
                                   matrix( 0,
                                           nrow = N.aug,
                                           ncol = nimConstants$nMaxDetectorsOth))
      
      ## z
      nimData$z <- nimData$z[detected[,t],t]  
      nimData$z[nimData$z %in% c(2)] <- 1 # ALIVE IDS BECOMES 1
      nimData$z <- c(nimData$z, rep(NA,N.aug))
      
      ## SXY 
      nimData$sxy <- NULL
      
      ## nbDetections 
      nimData$nbDetections <- nimData$nbDetections[detected[,t],t]
      nimData$nbDetections  <- c(nimData$nbDetections, rep(0,N.aug))
      nimData$nbDetectionsOth <- nimData$nbDetectionsOth[detected[,t],t]
      nimData$nbDetectionsOth  <- c(nimData$nbDetectionsOth, rep(0,N.aug))
      
      ## yDets 
      nimData$yDets <- nimData$yDets[detected[,t],,t]
      nimData$yDets <- rbind(nimData$yDets, matrix(0,nrow =N.aug, nimConstants$nMaxDetectors))
      nimData$yDetsOth <- nimData$yDetsOth[detected[,t],,t]
      nimData$yDetsOth <- rbind(nimData$yDetsOth, matrix(0,nrow =M, nimConstants$nMaxDetectorsOth))
      
      ## detResponse 
      nimData$detResponse <- nimData$detResponse[detected[,t],t]
      nimData$detResponse  <- c(nimData$detResponse, rep(NA,M))## HERE IT IS ASSUMING IT IS A LATENT INDIVIDUAL COVARIATE
      
      ## detCovs
      nimData$detCovs <- nimData$detCovs[,t,]
      nimData$detCovsOth <- nimData$detCovsOth[,t,]
      
      ## density
      nimData$denCounts <- nimData$denCounts[,1]
      
      
      
      ## ------   2.3. SUBSET NIMCONSTANTS ====
      
      ## countyToggle to toggle off Norrbotten
      nimConstants$countyToggle <- nimConstants$countyToggle[,t]
      nimConstants$countyToggleOth <- nimConstants$countyToggleOth[,t]
      
      
      
      ## ------   2.4. SUBSET NIMINITS ====
      
      ## z
      nimInits$z <- nimInits$z[detected[,t],t]  
      nimInits$z <- c(nimInits$z, rbinom(N.aug,1,0.5))
      
      ## detResponse
      ## HERE IT IS TREATED AS A LATENT COVARIATE
      nimInits$detResponse  <- c(rep(NA, sum(detected[,t])),
                                 rbinom(M,1,0.5))
      
      ## SXY 
      nimInits$sxy <- nimInits$sxy[detected[,t],,t]  
      ## GIVE ACS FROM DETECTED INDIVIDUALS TO AUGMENTED IDS. 
      nimInits$sxy <- rbind(nimInits$sxy,
                            nimInits$sxy[sample(nimInits$sxy, M, replace = T), ])
      
      ## p0
      nimInits$p01 <-  nimInits$p01[,t]
      nimInits$p01Oth <-  nimInits$p01Oth[,t]
      
      ## psi
      nimInits$psi <-  runif(1,0.4,0.6)
      
      ## sigma
      nimInits$sigma <-  runif(1,1,2)
      
      ## trapBetas
      nimInits$betaCovs <- nimInits$betaCovs[,t]
      nimInits$betaCovsOth <- nimInits$betaCovsOth[,t]
      
      ## trapBetas
      nimInits$betaResponse <- nimInits$betaResponse[t]
      nimInits$betaResponseOth <- nimInits$betaResponseOth[t]
      nimData$detCountries <- nimData$detCountries[ ,t]
      
      ## Get the new number of individuals 
      nimConstants$n.individuals <- nrow(nimInits$sxy)
      
      
      
      ## ------   2.5. SUBSET NIMPARAMS ====
      
      nimParams <- c("N", "psi", "pResponse","p0Oth","betaCovsOth","betaResponseOth",
                     "p0", "sigma", "betaDens", "betaCovs","betaResponse","betaResponseOth")
      
      nimParams2 <- c("z", "sxy")
      
      
      
      ## ------   2.6. SAVE SCR INPUT ====
      
      modelCode <- modelCode1
      
      save(nimData,
           nimConstants,
           nimParams,
           nimParams2,
           modelCode,
           nimInits,
           file = file.path(myVars$WD, myVars$modelName, "SCR", thisSex,
                            paste0("Snap",myVars$modelName,thisSex, years[t],"_", ch, ".RData")))
    }#t
  }#c
}#thisSex



## ------   3. CHECK MODEL LIKELIHOODS ====

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()



##------------------------------------------------------------------------------

