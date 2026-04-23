
## ------ IMPORT REQUIRED LIBRARIES ------

rm(list=ls())
gc()

library(raster)
library(coda)
library(nimble)
library(stringr)
library(abind)
library(R.utils)
library(adehabitatHR)
library(sf)
library(fasterize)
library(nimbleSCR)
library(readxl)
library(spatstat)
library(stars)
library(dplyr)
library(ggplot2)
library(rovquantR)


## ------ SOURCE THE REQUIRED FUNCTIONS ------

#source("C:/My_documents/RovQuant/workingDirectories.R")             
sourceDirectory("R")


## ------ SET WORKING DIRECTORIES ------

data.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2025/Data"
working.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2025/Test_F_NEW"

source("C:/My_documents/RovQuant/Temp/CM/functions/Nimble/dbinomLocal_normalWolf.R")

##------------------------------------------------------------------------------
####----------------------------- 
####---- cleanRovBaseData() -----

#cleanRovbaseData <- function(
    ##-- paths
#data.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2025/Data"
#working.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2025/Test_F_NEW"

##-- data
species = "wolf"
years = NULL 
two.sex = TRUE
sampling.months = NULL
rename.list = NULL
legal.dead = NULL 

##-- miscellanious
print.report = TRUE
Rmd.template = NULL
output.dir = NULL
overwrite = FALSE
#) {

plot.check = TRUE


####----------------------------- 

# ## ------ 1. INITIAL CHECKS -----
# 
# ##-- Make sure directory structure exists
# if(two.sex) {
#   makeDirectories( path = working.dir,
#                    subFolders = c("female","male"),
#                    show.dir = TRUE)
# } else {
#   makeDirectories( path = working.dir,
#                    show.dir = TRUE)
# }
# 
# ##-- Species
# if(length(species) > 1) {
#   stop('This function can only deal with one species at a time... \nPlease, use one of "bear", "wolf", or "wolverine" for the target species.')
# }
# if(sum(grep("bear", species, ignore.case = T)) > 0|
#    sum(grep("bjørn", species, ignore.case = T)) > 0|
#    sum(grep("bjorn", species, ignore.case = T)) > 0) {
#   SPECIES <- "Brown bear"
#   engSpecies <- "bear"
#   norSpecies <- c("Bjørn", "BjÃ¸rn")
# } else {
#   if(sum(grep("wolf", species, ignore.case = T))>0|
#      sum(grep("ulv", species, ignore.case = T))>0) {
#     SPECIES <- "Gray wolf"
#     engSpecies <- "wolf"
#     norSpecies <- "Ulv"
#   } else {
#     if(sum(grep("wolverine", species, ignore.case = T))>0|
#        sum(grep("järv", species, ignore.case = T))>0|
#        sum(grep("jerv", species, ignore.case = T))>0) {
#       SPECIES <- "wolverine"
#       engSpecies <- "wolverine"
#       norSpecies <- "Jerv"
#     } else {
#       SPECIES <- engSpecies <- norSpecies <- species
#     }
#   }
# }
# 
# ##-- Years
# if(is.null(years)) { years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }
# 
# ##-- Sampling months
# if(is.null(sampling.months)) {
#   if (engSpecies == "bear") {
#     sampling.months <- list(4:11)
#   } else {
#     if (engSpecies == "wolf") {
#       sampling.months <- list(c(10:12),c(1:3))
#     } else {
#       if (engSpecies == "wolverine") {
#         sampling.months <- list(c(10:12),c(1:4))
#       } else {
#         stop("No default setting available for the monitoring period of this species. \n You must specify the monitoring season months through the 'sampling.months' argument.")
#       }
#     }
#   }
# }
# 
# ##-- Legal mortality patterns
# if(is.null(legal.dead)) {
#   if (engSpecies %in% c("bear","wolf","wolverine")) {
#     legal.dead <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")
#   } else {
#     legal.dead <- ""
#   }
# }
# 
# ##-- Renaming list
# if(is.null(rename.list)) {
#   rename.list = c(
#     Age_estimated = "Alder, vurdert",
#     Age = "Alder, verifisert",
#     Age_verif_by = "Alder, verifisert av",
#     Age_class = "Alder på dødt individ",
#     Age_class_verif = "Aldersklasse verifisert SVA",
#     Analyzed_by = "AnalysertAv",
#     Analysis_priority = "Analyseprioritet",
#     Approved_by = "Godkjent av",
#     Approved_date = "Godkjentdato",
#     Assessment = "Vurdering",
#     Barcode_sample = "Strekkode (Prøve)",
#     Barcode = "Strekkode (Analyse)",
#     Birth_territory = "Født revir",
#     CITES = "CITES-nummer",
#     Collected_by = "Hvem samlet inn",
#     Collector_name = "Samlet selv - Navn",
#     Collector_phone = "Samlet selv - Telefon",
#     Collector_email = "Samlet selv - E-post",
#     Collector_role = "Samlet selv - Rolle",
#     Collector_other_name = "Annen innsamler - Navn" ,
#     Collector_other_phone = "Annen innsamler - Telefon",
#     Collector_other_email = "Annen innsamler - E-post",
#     Collector_other_role = "Annen innsamler - Rolle",
#     Comments_sample = "Merknad (Prøve)",
#     Comments = "Merknad (Analyse)",
#     Control_status = "Kontrollstatus",
#     Coordinate_system = "Koordinatsystem",
#     Counted_off_against_decision = "Regnes av mot vedtak",
#     County_number = "Fylkenummer",
#     County = "Fylke",
#     Date = "Funnetdato",
#     Date = "Dødsdato",
#     Death_cause = "Bakgrunn/årsak",
#     Death_method = "Bakgrunn/årsak metode",
#     Death_purpose = "Bakgrunn/årsak formål",
#     DNAID_sample = "DNAID (Prøve)",
#     DNAID = "DNAID (Analyse)",
#     EventID = "HendelseID",
#     East_Original = "Øst (opprinnelig)",
#     East_RT90 = "Øst (RT90)",
#     East_UTM33 = "Øst (UTM33/SWEREF99 TM)",
#     Felling_site_verif = "Kontroll av fellingsted",
#     Field_personnel = "Feltpersonell",
#     Hunting_date = "Observasjons/Jaktdato",
#     Id = "Individ",
#     Id = "Individ (Rovbase)",
#     IdSimplified = "ROVBASE_IndividID",
#     Juvenile = "Yngling",
#     Mountain_area = "Fjellområde",
#     Method = "Metode",
#     Municipality_number = "Kommunenummer",
#     Municipality = "Kommune",
#     North_original = "Nord (opprinnelig)",
#     North_RT90 = "Nord (RT90)",
#     North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
#     Origin = "Opprinnelse",
#     Outcome = "Utfall",
#     Last_saved_by_sample = "Sist lagret av (Prøve)",
#     Last_saved_sample = "Sist lagret dato (Prøve)",
#     Last_saved_by = "Sist lagret av (Analyse)",
#     Last_saved = "Sist lagret dato (Analyse)",
#     Last_saved_by = "Sist lagret av",
#     Last_saved =  "Sist lagret dato",
#     Locality = "Lokalitet",
#     Location = "Funnsted",
#     Lansstyrelsen_number = "Länsstyrelsens nr",
#     Quality_checked = "Kvalitetssikret av feltpersonell",
#     Quality_check_name = "Kvalitetssikrer - navn",
#     Quality_check_orga = "Kvalitetssikrer - Organisasjon",
#     Release_Date = "Frigivelsesdato",
#     Sample_type = "Prøvetype",
#     Sensitivity = "Følsomhet",
#     Species_sample = "Art (Prøve)",
#     Site_quality = "Stedkvalitet",
#     Time_of_death = "Dødstidspunkt",
#     Tips_name = "Tipser - Navn",
#     Tips_phone = "Tipser - Telefon",
#     Tips_email = "Tipser - E-post",
#     Tips_role = "Tipser - Rolle",
#     Tissue_sample = "Vevsprøve tatt",
#     Release_Date = "Frigivelsesdato",
#     RovbaseID = "RovbaseID (Analyse)",
#     RovbaseID_sample = "RovbaseID (Prøve)",
#     Species = "Art (Analyse)",
#     Species = "Art",
#     Sample_status = "Prøvestatus",
#     Sensitivity = "Følsomhet",
#     Sex_analysis = "Kjønn (Analyse)",
#     Sex = "Kjønn (Individ)",
#     Sex = "Kjønn",
#     Sex = "Kön",
#     Site_quality = "Stedkvalitet",
#     SVAID = "SVAID",
#     Uncertain_date = "Usikker dødsdato",
#     Weight_slaughter = "Slaktevekt",
#     Weight_total =  "Helvekt")
# }
# 
# ##-- Load pre-processed habitat shapefiles
# #data(COUNTRIES, envir = environment()) [PD]: WE SHOULD USE REGIONS FOR EVERYTHING !!!
# data(REGIONS, envir = environment())
# 
# ##-- data info
# DATE <- getMostRecent(path = data.dir, pattern = "DNA")
# 
# ##-- Set file name for clean data
# fileName <- paste0("CleanData_", engSpecies, "_", DATE, "_2.RData")
# 
# ##-- Check that a file with that name does not already exist to avoid overwriting
# if(!overwrite) {
#   existTest <- file.exists(file.path(working.dir, "data", fileName))
#   if (any(existTest)) {
#     message(paste0("A file named '", fileName[existTest], "' already exists in: \n",
#                    file.path(working.dir, "data")))
#     message("Are you sure you want to proceed and overwrite existing clean data? (y/n) ")
#     question1 <- readLines(n = 1)
#     if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
#       message("Not overwriting existing files...")
#       return(invisible(NULL))
#     } else {
#       message(paste0("Now overwriting '", fileName[existTest],"'.\n"))
#     }
#   }
# }
# 
# ## [PD] REMOVE?
# # n.years <- length(years)
# # YEARS <- lapply(years, function(x)c(x,x+1))
# 
# 
# 
# ## ------ 2. CLEAN THE DATA -----
# 
# ## ------   2.1. RAW NGS DATA -----
# 
# ##-- NGS data
# ## [PD]: REMOVE
# # DNA <- read.csv( file.path(data.dir, "RIB22042025133456403_wolfDNA.csv"),
# #                  fileEncoding = "latin1")
# # colnames(DNA) <- translateForeignCharacters( dat = colnames(DNA),
# #                                              dir.translation = dir.analysis)
# 
# DNA <- suppressWarnings(readMostRecent( path = data.dir,
#                                         extension = ".xls",
#                                         pattern = "DNA")) %>%
#   ##-- Rename columns to facilitate manipulation
#   dplyr::rename(., any_of(rename.list)) %>%
#   ##-- Turn potential factors into characters
#   dplyr::mutate(across(where(is.factor), as.character)) %>%
#   ##-- Initial filters
#   dplyr::filter(
#     ##-- Filter to the focal species
#     Species %in% norSpecies,
#     # [CHECK] Should we do this for all species now ???
#     # ##-- Filter dead recoveries (HB for the last wolverine analysis)
#     # !substr(RovbaseID_sample,1,1) %in% "M"
#   ) %>%
#   ##-- Remove any duplicates
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   ##-- Add some columns
#   dplyr::mutate(
#     ##-- Add "Country" column
#     Country_sample = substrRight(County, 3),
#     ##-- Change date format
#     Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
#     ##-- Extract year
#     Year = as.numeric(format(Date,"%Y")),
#     ##-- Extract month
#     Month = as.numeric(format(Date,"%m")),
#     ##-- Extract sampling season
#     ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
#     ##-- Set all months in given sampling period to the same year)
#     Year = ifelse( Month < unlist(sampling.months)[1],
#                    Year-1,
#                    Year),
#     ##-- Fix unknown "Id"
#     Id = ifelse(Id %in% "", NA, Id),
#     ##-- Fix unknown "Sex"
#     Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown", Sex),
#     #Sex = ifelse(is.na(Sex), "unknown", Sex),
#     Sex = ifelse(Sex %in% "Hunn", "female", Sex),
#     Sex = ifelse(Sex %in% "Hann", "male", Sex))
# # [PD] Should we filter for years here ???
# # %>%
# # ##-- Filter to the focal years
# # dplyr::filter(., Year %in% years)
# 
# 
# ##-- Number of NGS samples
# NGS_samples <- table(DNA$Sex, DNA$Year, useNA = "ifany")
# NGS_samples <- rbind(NGS_samples, "Total" = colSums(NGS_samples))
# NGS_samples <- cbind(NGS_samples, "Total" = rowSums(NGS_samples))
# write.csv( NGS_samples,
#            file = file.path( working.dir, "tables",
#                              paste0(engSpecies, "_Raw NGS Samples_",
#                                     years[1]," to ", years[length(years)],
#                                     ".csv")))
# 
# 
# ##-- Number of individuals detected alive
# NGS_ids <- apply(table(DNA$Sex, DNA$Year, DNA$Id, useNA = "ifany"), c(1,2), function(x)sum(x>0))
# NGS_ids <- rbind(NGS_ids, "Total" = apply(table(DNA$Year,DNA$Id, useNA = "ifany"), 1, function(x)sum(x>0)))
# NGS_ids <- cbind(NGS_ids, "Total" = c(apply(table(DNA$Sex,DNA$Id, useNA = "ifany"), 1, function(x)sum(x>0)),length(unique(DNA$Id))))
# write.csv( NGS_ids,
#            file = file.path( working.dir, "tables",
#                              paste0(engSpecies, "_Raw NGS Ids_",
#                                     years[1]," to ", years[length(years)],
#                                     ".csv")))
# 
# 
# 
# ## ------   2.2. RAW DEAD RECOVERY DATA -----
# 
# ##-- Load raw excel file imported from rovbase
# ## [PD]: REMOVE
# # DEAD <- read.csv( file.path(data.dir, "RIB22042025133534832_wolfDEAD.csv"),
# #                   fileEncoding = "latin1")
# # colnames(DEAD) <- translateForeignCharacters(dat = colnames(DEAD), dir.translation = dir.analysis )
# DR <- suppressWarnings(readMostRecent( path = data.dir,
#                                        extension = ".xls",
#                                        pattern = "dead")) %>%
#   ##-- Rename columns to facilitate manipulation
#   dplyr::rename(., any_of(rename.list)) %>%
#   ##-- Initial filters
#   dplyr::filter(
#     ##-- Filter to the focal species
#     Species %in% norSpecies) %>%
#   ##-- Remove any duplicates
#   dplyr::distinct(., .keep_all = TRUE) %>%
#   ##-- Turn potential factors into characters
#   dplyr::mutate(across(where(is.factor), as.character)) %>%
#   ##-- Add some columns
#   dplyr::mutate(
#     ##-- Add "Country" column
#     Country_sample = substrRight(County, 3),
#     ##-- Change date format
#     Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
#     ##-- Extract year
#     Year = as.numeric(format(Date,"%Y")),
#     ##-- Extract month
#     Month = as.numeric(format(Date,"%m")),
#     ##-- Extract sampling season
#     ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
#     ##-- Set all months in given sampling period to the same year)
#     Year = ifelse( Month < unlist(sampling.months)[1],
#                    Year-1,
#                    Year),
#     ##-- Fix unknown "Id"
#     Id = ifelse(Id %in% "", NA, Id),
#     ##-- Fix unknown "Sex"
#     Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown" , Sex),
#     #Sex = ifelse(is.na(Sex), "unknown", Sex),
#     Sex = ifelse(Sex %in% "Hunn", "female", Sex),
#     Sex = ifelse(Sex %in% "Hann", "male", Sex),
#     ##-- Identify legal deaths
#     Legal = grepl(paste(legal.dead, collapse = "|"), Death_cause))
# # [PD] Should we filter for years here ???
# # %>%
# # ##-- Filter to the focal years
# # dplyr::filter(., Year %in% years)
# 
# ##-- Number of DR samples
# DR_samples <- table(DR$Sex, DR$Year, useNA = "ifany")
# DR_samples <- rbind(DR_samples, "Total" = colSums(DR_samples))
# DR_samples <- cbind(DR_samples, "Total" = rowSums(DR_samples))
# write.csv( DR_samples,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Raw DR Samples_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# ##-- Number of individuals recovered dead
# DR_ids <- apply(table(DR$Sex, DR$Year, DR$Id, useNA = "ifany"), c(1,2), function(x)sum(x>0))
# DR_ids <- rbind(DR_ids, "Total" = apply(table(DR$Year,DR$Id, useNA = "ifany"), 1, function(x)sum(x>0)))
# DR_ids <- cbind(DR_ids, "Total" = c(apply(table(DR$Sex,DR$Id, useNA = "ifany"), 1, function(x)sum(x>0)),length(unique(DR$Id))))
# write.csv( DR_ids,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Raw DR Ids_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# 
# 
# ## ------   2.3. CHECKS & FILTERS -----
# 
# ##-- Filter out unusable samples
# numNoID_DNA <- sum(is.na(DNA$Id))              ## number of samples without ID
# numNoDate_DNA <- sum(is.na(DNA$Year))          ## number of samples without Date
# numNoCoords_DNA <- sum(is.na(DNA$East_UTM33))  ## number of samples without Coords
# 
# DNA <- DNA %>%
#   dplyr::filter(##-- Filter out samples with no ID
#     !is.na(Id),
#     ##-- Filter out samples with no Coordinates
#     !is.na(East_UTM33),
#     ##-- Filter out samples with no dates
#     !is.na(Year))
# 
# 
# ##-- Filter out unusable samples
# numNoID_DR <- sum(is.na(DR$Id))                ## number of DR without ID
# numNoDate_DR <- sum(is.na(DR$Year))            ## number of DR without Date
# numNoCoords_DR <- sum(is.na(DR$East_UTM33))    ## number of DR without Coords
# 
# DR <- DR %>%
#   dplyr::filter(##-- Filter out samples with no ID
#     !is.na(Id),
#     ##-- Filter out samples with no Coordinates
#     !is.na(East_UTM33),
#     ##-- Filter out samples with no dates
#     !is.na(Year))
# 
# 
# ##-- Filter out problematic samples
# 
# ##-- Number of 'dead recovery" samples in DNA only: DNAID
# numDNAID_inDNA_notinDR <- sum(!DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$DNAID)
# DNAID_inDNA_notinDR <- NULL
# if(numDNAID_inDNA_notinDR > 0){
#   ##-- Identify dead recoveries only in DNA
#   tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
#   DNAID_inDNA_notinDR <- tmp[!tmp$DNAID %in% DR$DNAID, c("DNAID", "RovbaseID", "Id")]
# }
# 
# ##-- Number of 'dead recovery" samples in DNA only: RovbaseID
# numRovbaseID_inDNA_notinDR <- sum(!DNA$RovbaseID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$RovbaseID)
# RovbaseID_inDNA_notinDR <- NULL
# if(numRovbaseID_inDNA_notinDR > 0){
#   ##-- Identify dead recoveries only in DNA
#   tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
#   RovbaseID_inDNA_notinDR <- tmp[!tmp$RovbaseID %in% DR$RovbaseID, c("DNAID", "RovbaseID", "Id")]
# }
# 
# ##-- Duplicated dead recoveries
# numDupId_DR <- sum(duplicated(DR$Id))
# dupId_DR <- NULL
# if(numDupId_DR > 0){
#   ##-- Identify duplicated dead recoveries
#   dupId_DR <- DR[DR$Id %in% DR$Id[duplicated(DR$Id)], c("DNAID", "RovbaseID", "Id")]
#   ##-- Remove duplicated individuals (keeping the last occurrence)
#   #DR <- dplyr::filter(DR, !duplicated(Id, fromLast = T))
# }
# 
# # inDNA_notinDR <- NULL
# # if(numInDNA_notinDR > 0){
# #   ##-- Identify dead recoveries only in DNA
# #   tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
# #   inDNA_notinDR <- tmp[!tmp$DNAID %in% DR$DNAID, c("DNAID", "RovbaseID", "Id")]
# # }
# 
# ##-- Check which data is duplicated
# duplicateData <- dplyr::inner_join( DNA, DR,
#                                     by = c("Id","RovbaseID","DNAID","Species","Sex",
#                                            "Date","Year","Month",
#                                            "East_UTM33","North_UTM33",
#                                            "County","Country_sample"))
# numDupData <- nrow(duplicateData)
# if(numDupData > 0){
#   write.csv( duplicateData,
#              file = file.path( working.dir, "tables",
#                                paste0( engSpecies, "_DR in DNA_",
#                                        years[1]," to ", years[length(years)],
#                                        ".csv")))
#   ##-- Remove duplicated data in DNA before merging
#   # DNA <- DNA[!DNA$DNAID %in% duplicateData$DNAID, ]
# }
# 
# 
# 
# ## ------   2.4. MERGE -----
# 
# ##-- Merge DNA and dead recoveries files using all shared names columns
# # DATA <- full_join(DNA, DR, by = names(DNA)[names(DNA) %in% names(DR)])
# DATA <- merge( DR, DNA,
#                by = c("Id","RovbaseID","DNAID","Species","Sex",
#                       "Date","Year","Month",
#                       "East_UTM33","North_UTM33",
#                       "County","Country_sample"),
#                all = TRUE)
# 
# 
# 
# ## ------   2.5. AGE -----
# 
# ##-- Determine Death and Birth Years
# DATA <- DATA %>%
#   dplyr::mutate(
#     Age = suppressWarnings(as.numeric(as.character(Age))),
#     Death = ifelse(substr(RovbaseID,1,1) %in% "M", Year, NA),
#     Birth = Death - Age)
# 
# 
# 
# ## ------   2.6. SEX ASSIGNMENT -----
# 
# ##-- If this is the wolf data, we first consolidate all the info we have about
# ##-- individual sex from the different data sources.
# if(engSpecies == "wolf"){
# 
#   ##-- Add simplified ID column
#   DATA <- DATA %>%
#     dplyr::mutate(IdSimplified = unlist(lapply(strsplit(Id, " "), function(x) x[1])))
# 
#   ##-- Load most recent Micke's file
#   INDIVIDUAL_ID <- suppressWarnings(readMostRecent( path = data.dir,
#                                                     extension = ".xls",
#                                                     pattern = "_ID Grouping")) %>%
#     ##-- Rename columns to facilitate manipulation
#     dplyr::rename(., any_of(rename.list),
#                   Year = "ReprodYear (May 1 year y - Apr 30 y+1)") %>%
#     ##-- Turn potential factors into characters
#     dplyr::mutate(across(where(is.factor), as.character)) %>%
#     ##-- Add some columns
#     dplyr::mutate(
#       IdSimplified = unlist(lapply(strsplit(Id, " "), function(x) x[1])),
#       Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
#       Sex = ifelse(is.na(Sex), "unknown", Sex),
#       Sex = ifelse(Sex %in% "Hona", "female", Sex),
#       Sex = ifelse(Sex %in% "Hane", "male", Sex))
# 
# 
#   ##-- THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2022/23.
#   Pack_ID2023 <- suppressWarnings(readMostRecent( path = data.dir,
#                                                   extension = ".xls",
#                                                   pattern = "Genetiskt ID")) %>%
#     ##-- Rename columns to facilitate manipulation
#     dplyr::rename(., any_of(rename.list),
#                   IdSimplified = "Rovbase-ID") %>%
#     ##-- Turn potential factors into characters
#     dplyr::mutate(across(where(is.factor), as.character)) %>%
#     ##-- Add some columns
#     dplyr::mutate(
#       ##-- Add status
#       Status = "Pair",
#       ##-- Fix unknown "Sex"
#       Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
#       Sex = ifelse(is.na(Sex), "unknown", Sex),
#       Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
#       Sex = ifelse(Sex %in% c("Hann","Hane"), "male", Sex),
#       Year = 2023)
# 
# 
#   ##-- THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2023/24.
#   Pack_ID2024 <- suppressWarnings(readMostRecent( path = data.dir,
#                                                   extension = ".xls",
#                                                   pattern = "Bilaga_")) %>%
#     ##-- Rename columns to facilitate manipulation
#     dplyr::rename(.,
#                   IdSimplified = "RovbaseID",
#                   any_of(rename.list)) %>%
#     ##-- Turn potential factors into characters
#     dplyr::mutate(across(where(is.factor), as.character)) %>%
#     ##-- Add some columns
#     dplyr::mutate(
#       ##-- Add status
#       Status = "Pair",
#       ##-- Fix unknown "Sex"
#       Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
#       Sex = ifelse(is.na(Sex), "unknown", Sex),
#       Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
#       Sex = ifelse(Sex %in% c("Hann","Hane"), "male", Sex),
#       Year = 2024)
# 
# 
#   ##-- THIS IS THE PACK ID SENT BY ØYSTEIN FOR THE WINTER 2024/25.
#   Pack_ID2025 <- suppressWarnings(readMostRecent( path = data.dir,
#                                                   extension = ".xls",
#                                                   pattern = "FromOystein")) %>%
#     ##-- Rename columns to facilitate manipulation
#     dplyr::rename(., any_of(rename.list),
#                   IdSimplified = "IndividID") %>%
#     ##-- Turn potential factors into characters
#     dplyr::mutate( across(where(is.factor), as.character),
#                   Year = 2025,
#                   Status = "Pair") %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(
#       ##-- Fix unknown "Sex"
#       Sex = ifelse(any(c_across(Sex1:Sex4) %in% c("Tispe","Tik")),
#                    "female",
#                    ifelse(any(c_across(Sex1:Sex4) %in% c("Hann","Hane")),
#                           "male",
#                           "unknown")))
# 
# 
#   ##-- CONSOLIDATE ALL INFO ON INDIVIDUAL SEX IN ONE DATAFRAME
#   ALL_SEX <- rbind( cbind("IdSimplified" = unlist(lapply(strsplit(DATA$Id, " "), function(x) x[1])),
#                           "Sex" = DATA$Sex),
#                     INDIVIDUAL_ID[ ,c("IdSimplified","Sex")],
#                     Pack_ID2023[ ,c("IdSimplified","Sex")],
#                     Pack_ID2024[ ,c("IdSimplified","Sex")],
#                     Pack_ID2025[ ,c("IdSimplified","Sex")])
# 
# 
#   ##-- CONSOLIDATE ALL INFO ON INDIVIDUAL STATUS IN ONE DATAFRAME
#   ALL_STATUS <-  rbind( INDIVIDUAL_ID[ ,c("IdSimplified","Year","Status")],
#                         Pack_ID2023[ ,c("IdSimplified","Year","Status")],
#                         Pack_ID2024[ ,c("IdSimplified","Year","Status")],
#                         Pack_ID2025[ ,c("IdSimplified","Year","Status")])
# 
#   DATA <- DATA %>%
#     left_join(., ALL_STATUS, by = c("IdSimplified","Year"))
# }#if

#   ############################################################################
#   # ### [PD] CHECK ###
#   # table(INDIVIDUAL_ID$Sex[INDIVIDUAL_ID$Id == "UI414251 G63-20"])
#   # table(DATA$Sex[DATA$Id == "UI414251 G63-20"])
#   #
#   # ##-- Overwrite gender from Micke's data when available
#   # micke.sex <- lapply(DATA$Id,
#   #                     function(i){
#   #                       INDIVIDUAL_ID[INDIVIDUAL_ID$Id %in% i, "Sex"][1]
#   #                     })
#   # DATA$Sex <- ifelse(!is.na(micke.sex), micke.sex, DATA$Sex)
#   #
#   # numOverwiteSex <- sum(unique(INDIVIDUAL_ID$`Individ (Rovbase)`) %in% DATA$Id)
#   #
#   # micke.sex <- as.character(unlist(lapply(myCleanedData.sp$Id, function(i) INDIVIDUAL_ID[as.character(INDIVIDUAL_ID$Individ..Rovbase.)==i,"Sex"][1])))
#   # micke.sex[micke.sex %in% "0"] <- NA
#   # micke.sex[micke.sex %in% names(table(micke.sex))[3]] <- NA
#   # micke.sex[micke.sex %in% "Hona"] <- "Hunn"
#   # micke.sex[micke.sex %in% "Hane"] <- "Hann"
#   # table(!is.na(micke.sex))
#   # new.sex <- ifelse(!is.na(micke.sex), as.character(micke.sex), as.character(myCleanedData.sp$Sex))
#   # table(myCleanedData.sp$Sex, new.sex)
#   # myCleanedData.sp$Sex <- new.sex
#   # table(myCleanedData.sp$Sex, new.sex)
#   #
#   # ## [PD] : REMOVE?
#   # ##-- Make a simplified column to match the rovbase id given by Oystein in Linn's file
#   # DATA$IdSimplified <- unlist(lapply(strsplit(as.character(DATA$Id), " "), function(x) x[1]))
#   #
#   # ##-- OVERWRITE GENDER FROM PACK COMPOSITION (FROM LINN's file 2023-24)
#   # ##-- check the sex in the pair data given by Linn and assign the sex to all detections
#   # ##-- Overwrite sex
#   # for(i in 1:nrow(Pack_ID2023)){
#   #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2023$Rovbase.ID[i]] <- Pack_ID2023$SEx[i]
#   # }#i
#   #
#   # ##-- Overwrite sex
#   # for(i in 1:nrow(Pack_ID2024)){
#   #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2024$Rovbase.ID[i]] <- Pack_ID2024$Sex[i]
#   # }#i
#   #
#   # ##-- Overwrite sex
#   # for(i in 1:nrow(Pack_ID2025)){
#   #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2025$IndividID[i]] <- Pack_ID2025$Sex[i]
#   # }#i
# }
# 
# ##-- Loop over all individuals
# ID <- unique(as.character(DATA$Id))
# doubleSexID <- IdDoubleSex <- NULL
# counter <- 1
# for(i in 1:length(ID)){
#   ##-- Subset data to individual i
#   if(engSpecies == "wolf"){
#     tmp <- ALL_SEX$Sex[ALL_SEX$IdSimplified == unlist(lapply(strsplit(ID[i], " "), function(x) x[1]))]
#   } else{
#     tmp <- DATA$Sex[DATA$Id == ID[i]]
#   }
# 
#   ##-- Number of times individual i was assigned to each sex
#   tab <- table(tmp[tmp %in% c("female","male")])
# 
#   ##-- If conflicting sexes (ID identified as both "female" and "male")
#   if(length(tab) == 2){
#     ##-- If ID assigned the same number of times to the 2 sexes, assign to unknown
#     if(tab[1] == tab[2]){
#       DATA$Sex[DATA$Id == ID[i]] <- "unknown"
#     } else {
#       ##-- Otherwise pick the most common sex
#       DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
#     }
#     # print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))]))
#     IdDoubleSex[counter] <- ID[i]
#     counter <- counter + 1
#   }
# 
#   ##-- If only one of "female" or "male" registered
#   if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
# 
#   ##-- If anything else registered : "unknown"
#   if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "unknown"}
# 
#   ##-- Track number of sexes assigned for this individual
#   ##-- (0 == "unknown", 2 == "both sexes)
#   doubleSexID[i] <- length(tab)
# }#i
# 
# 
# 
# ## ------   2.8. SPLIT DATA -----
# 
# ##-- Split DATA into alive and dead.recovery datasets
# alive <- DATA[is.na(DATA$Death), ]
# dead.recovery <- DATA[!is.na(DATA$Death), ]
# 
# ##-- Add earlier detection index
# alive$detected.earlier <-
#   unlist(lapply(1:nrow(alive),
#                 function(i){
#                   this.id <- alive[i,"Id"]
#                   this.date <- alive[i,"Date"]
#                   any(alive$Id %in% this.id & alive$Date < this.date)
#                 }))
# 
# dead.recovery$detected.earlier <-
#   unlist(lapply(1:nrow(dead.recovery),
#                 function(i){
#                   this.id <- dead.recovery[i,"Id"]
#                   this.date <- dead.recovery[i,"Date"]
#                   any(alive$Id %in% this.id & alive$Date < this.date)
#                 }))
# 
# 
# 
# # ## ------     2.6.2. FILTER DATA FOR SEX -----
# #
# # # myFullData.sp <- FilterDatasf(
# # #   myData = myCleanedData.sp,
# # #   dead.recovery = T,
# # #   sex = DATA$sex,
# # #   setSex = T)
# #
# # ## [PD] part of the cleanRovBaseData function now!!!!
# #
# # ##-- List all individual IDs
# # ID <- unique(as.character(myCleanedData.sp$Id))
# # myCleanedData.sp$Sex <- as.character(myCleanedData.sp$Sex)
# #
# # ##-- Initialize the vector of IDs with conflicting sexes
# # IdDoubleSex <- 0
# # counter <- 1
# #
# # for(i in 1:length(ID)){
# #   ##-- subset data to individual i
# #   tmp <- myCleanedData.sp$Sex[myCleanedData.sp$Id == ID[i]]
# #   ##-- create a table of the number of times individual i was assigned to each sex
# #   tab <- table(tmp[tmp %in% c("female","male")])
# #   ##-- If conflicting sexes (ID identified as both "Hunn" and "Hann")
# #   if(length(tab) == 2){
# #     ##-- If ID assigned the same number of times to the 2 sexes, assign to Ukjent
# #     if(tab[1] == tab[2]){
# #       myCleanedData.sp$Sex[myCleanedData.sp$Id == ID[i]] <- "unknown"
# #     } else {
# #       ##-- Otherwise pick the most common sex
# #       myCleanedData.sp$Sex[myCleanedData.sp$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
# #     }
# #     # print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))]))
# #     IdDoubleSex[counter] <- ID[i]
# #     counter <- counter + 1
# #   }
# #   ##-- If only one of "female" or "male" registered
# #   if(length(tab) == 1){myCleanedData.sp$Sex[myCleanedData.sp$Id == ID[i]] <- names(tab)}
# #
# #   ##-- If anything else registered : "unknown"
# #   if(length(tab) == 0){myCleanedData.sp$Sex[myCleanedData.sp$Id == ID[i]] <- "Ukjent"}
# # }#i
# #
# # myData.dead <- myData[!is.na(myData$Death), ]
# # myData.alive <- myData[is.na(myData$Death), ]
# #
# # myData.dead$Id <- droplevels(myData.dead$Id)
# # myData.alive$Id <- droplevels(myData.alive$Id)
# #
# # IdDoubleDead <- myData.dead$Id[duplicated(myData.dead$Id)]
# #
# # myFullData.sp <- list( alive = myData.alive,
# #                        dead.recovery = myData.dead,
# #                        IdDoubleSex = IdDoubleSex,
# #                        IdDoubleDead = IdDoubleDead)
# #
# #
# #
# # ## ------     2.6.3. REMOVE INDVIDUALS THAT DIED TWICE ------
# #
# # duplicatedDeath <- NULL
# # for(i in myFullData.sp$IdDoubleDead){
# #   tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
# #   if(length(tmp)==0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-1]}
# #   duplicatedDeath <- c(duplicatedDeath, tmp)
# # }#i
# #
# # myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]
# 
# 
# 
# ## ----- 4. DATA ISSUES -----
# 
# ## -----   4.1. MULTIPLE DEATHS ------
# 
# # ##-- Identify and count individuals dead "more than once"
# # ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
# # multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
# # for (i in 1:length(ID)) {
# #   tmp <- dead.recovery[dead.recovery$Id == ID[i], ]
# #   ##-- Multiple death dates
# #   if(length(unique(tmp$Date)) > 1){
# #     multiDeathDate <- c(multiDeathDate, ID[i])
# #   }
# #   ##-- Multiple death years
# #   if(length(unique(tmp$Year)) > 1){
# #     multiDeathYear <- c(multiDeathYear, ID[i])
# #   }
# #   ##-- Multiple death locations
# #   if(length(unique(tmp$East)) > 1 | length(unique(tmp$North)) > 1){
# #     multiDeathLocs <- c(multiDeathLocs, ID[i])
# #   }
# # }#i
# 
# ##-- Remove individuals that died more than once
# IdDoubleDead <- dead.recovery$Id[duplicated(dead.recovery$Id)]
# duplicatedDeath <- NULL
# if(length(IdDoubleDead) > 0){
#   for(i in IdDoubleDead){
#     tmp <- which(dead.recovery$Id == i & is.na(dead.recovery$Death_cause))
#     if(length(tmp)==0){tmp  <- which(dead.recovery$Id == i)[-1]}
#     duplicatedDeath <- c(duplicatedDeath, tmp)
#   }#i
#   dead.recovery <- dead.recovery[-duplicatedDeath, ]
# }#if
# 
# 
# 
# ## -----   4.2. GHOST INDIVIDUALS ------
# 
# id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
# ghosts <- unlist(lapply(id.list, function(id) {
#   out <- NULL
#   try({
#     if(id %in% dead.recovery$Id){
#       mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
#       this.alive <- alive[alive$Id == id, ]
#       ##-- Was it detected alive in any season after death?
#       temp <- this.alive[this.alive$Year > mort.year, ]
#       if(length(temp) > 0){
#         out <- rownames(temp)
#         names(out) <- id
#       }
#     }
#   }, silent = TRUE)
#   return(out)
# }))
# samples.to.remove <- unlist(ghosts)
# 
# ##-- Remove flagged NGS detections after dead recovery
# alive <- alive[!rownames(alive) %in% samples.to.remove, ]
# 
# 
# 
# ## ----- 5. TURN INTO .sf OBJECTS -----
# 
# ##-- Turn into sf points dataframe
# alive <- sf::st_as_sf( x = alive,
#                        coords = c("East_UTM33","North_UTM33")) %>%
#   sf::st_set_crs(., sf::st_crs(32633))
# 
# ##-- Intersect and extract country name
# alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
# alive$Country_sf[!is.na(as.numeric(sf::st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
# 
# ##-- Turn into sf points dataframe
# dead.recovery <- sf::st_as_sf( x = dead.recovery,
#                                coords = c("East_UTM33","North_UTM33")) %>%
#   sf::st_set_crs(.,sf::st_crs(32633))
# 
# ##-- Intersect and extract country name
# dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
# dead.recovery$Country_sf[!is.na(as.numeric(sf::st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
# 
# 
# 
# ## ----- 6. DATA SUMMARY -----
# 
# ## -----   6.1. DATA SUMMARY - TABLES -----
# 
# ##-- Number of NGS samples per year and country (date,rovbase)
# samples <- table(alive$Country_sample, alive$Year)
# samples <- rbind(samples, "Total" = colSums(samples))
# samples <- cbind(samples, "Total" = rowSums(samples))
# write.csv( samples,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Clean NGS Samples_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# ##-- Number of individuals detected alive
# ids <- apply(table(alive$Country_sample, alive$Year, alive$Id), c(1,2), function(x)sum(x>0))
# ids <- rbind(ids, "Total" = apply(table(alive$Year, alive$Id), 1, function(x)sum(x>0)))
# ids <- cbind(ids, "Total" = c(apply(table(alive$Country_sample,alive$Id), 1, function(x)sum(x>0)), length(unique(alive$Id))))
# write.csv( ids,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Clean NGS Ids_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# ##-- Number of DR samples
# deadSamples <- table(dead.recovery$Country_sample, dead.recovery$Year)
# deadSamples <- rbind(deadSamples, "Total" = colSums(deadSamples))
# deadSamples <- cbind(deadSamples, "Total" = rowSums(deadSamples))
# write.csv( deadSamples,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Clean DR Samples_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# ##-- Number of individuals recovered
# deadIds <- apply(table(dead.recovery$Country_sample,dead.recovery$Year,dead.recovery$Id),c(1,2),function(x)sum(x>0))
# deadIds <- rbind(deadIds, "Total" = apply(table(dead.recovery$Year,dead.recovery$Id), 1, function(x)sum(x>0)))
# deadIds <- cbind(deadIds, "Total" = c(apply(table(dead.recovery$Country_sample,dead.recovery$Id), 1, function(x)sum(x>0)), length(unique(dead.recovery$Id))))
# write.csv( deadIds,
#            file = file.path( working.dir, "tables",
#                              paste0( engSpecies, "_Clean DR Ids_",
#                                      years[1]," to ", years[length(years)],
#                                      ".csv")))
# 
# 
# 
# ## -----   6.2. NUMBER OF SAMPLES - FIGURE -----
# 
# ##-- Number of NGS per month
# dat.alive <- alive %>%
#   dplyr::mutate(Date = trunc(Date, "month")) %>%
#   dplyr::group_by(Date) %>%
#   dplyr::summarise(n = dplyr::n())
# dat.alive$type = "NGS"
# 
# ##-- Number of dead recoveries per month
# dat.dead <- dead.recovery %>%
#   dplyr::mutate(Date = trunc(Date, "month")) %>%
#   dplyr::group_by(Date) %>%
#   dplyr::summarise(n = -dplyr::n())
# dat.dead$type = "dead.recovery"
# 
# ##-- Combine NGS and dead recoveries
# dat <- rbind(dat.alive, dat.dead)
# dat$Date <- as.Date(dat$Date)
# 
# ##-- Plot time series of number of samples per month
# NGS_plot <- ggplot(dat) +
#   geom_col(aes(x = Date, y = n, fill = type)) +
#   ylab("Number of samples") +
#   guides(fill = guide_legend(reverse = TRUE)) +
#   theme(legend.title = element_blank(),
#         legend.position.inside = c(0.1,0.9),
#         axis.text.x = element_text(angle = 60,
#                                    hjust = 1)) +
#   scale_x_date( date_breaks = "years",
#                 date_labels = "%Y")
# 
# ##-- Export as .png
# ggsave(filename = file.path(working.dir, "figures",
#                             paste0( engSpecies, "_Clean Rovbase Samples_",
#                                     years[1]," to ", years[length(years)],
#                                     ".png")),
#        plot = NGS_plot,
#        dpi = 300, height = 6, width = 12, device = "png")
# 
# # grDevices::png( filename = file.path(working.dir, "figures",
# #                                      paste0( engSpecies, "_Clean Rovbase Samples_",
# #                                              years[1]," to ", years[length(years)],
# #                                              ".png")),
# #                 width = 8, height = 6,
# #                 units = "in", pointsize = 12,
# #                 res = 300, bg = NA)
# # NGS_plot
# # dev.off()
# 
# 
# 
# ## -----   6.3. NUMBER OF INDIVIDUALS - FIGURE -----
# 
# ##-- Number of IDs
# dat.alive <- alive %>%
#   dplyr::group_by(Year) %>%
#   dplyr::summarise(n = length(unique(Id)))
# dat.alive$type = "NGS"
# 
# ##-- Number of Dead Recoveries
# dat.dead <- dead.recovery %>%
#   dplyr::group_by(Year) %>%
#   dplyr::summarise(n = -length(unique(Id)))
# dat.dead$type = "dead.recovery"
# 
# ##-- Combine NGS and dead recoveries
# dat <- rbind(dat.alive, dat.dead)
# 
# ##-- Plot time series of number of IDs per year
# IDS_plot <- ggplot2::ggplot(dat) +
#   geom_col(aes(x = Year, y = n, fill = type)) +
#   ylab("Number of individuals") +
#   guides(fill = guide_legend(reverse = TRUE)) +
#   theme(legend.title = element_blank(),
#         legend.position.inside = c(0.1,0.9),
#         axis.text.x = element_text(angle = 60,
#                                    hjust = 1)) +
#   scale_x_continuous( breaks = years,
#                       labels = years)
# ##-- Export as .png
# ggsave(filename = file.path(working.dir, "figures",
#                             paste0( engSpecies, "_Clean Rovbase Ids_",
#                                     years[1]," to ", years[length(years)],
#                                     ".png")),
#        plot = IDS_plot,
#        dpi = 300, height = 6, width = 12, device = "png")
# 
# # grDevices::png( filename = file.path(working.dir, "figures",
# #                                      paste0(engSpecies, "_Clean Rovbase Ids_",
# #                                             years[1]," to ", years[length(years)],
# #                                             ".png")),
# #                 width = 8, height = 6,
# #                 units = "in", pointsize = 12,
# #                 res = 300, bg = NA)
# # IDS_plot
# # dev.off()
# 
# 
# 
# ## -----   6.4. SAMPLING MAPS - FIGURE ------
# 
# ##-- Maps layout
# L <- length(years)
# if(L < 6){ nrows <- 1 } else{
#   if(L < 13){ nrows <- 2 } else {
#     if(L < 22){ nrows <- 3 } else {
#       if(L < 33){ nrows <- 4 } else {
#         nrows <- 5
#       }}}}
# ncols <- ceiling(L/nrows)
# 
# ##-- Save maps as .png
# grDevices::png(filename = file.path( working.dir, "figures",
#                                      paste0( engSpecies, "_Clean Rovbase Maps_",
#                                              years[1]," to ", years[length(years)],
#                                              ".png")),
#                width = ncols*2, height = nrows*4,
#                units = "in", pointsize = 12,
#                res = 300, bg = NA)
# 
# ##-- Maps layout
# mx <- matrix(NA, nrow = nrows*2, ncol =  (ncols*2)+1)
# for(r in 1:nrows){
#   mx[r*2-1, ] <- c(1,rep(1:ncols, each = 2)) + (r-1)*ncols
#   mx[r*2, ] <- c(rep(1:ncols, each = 2),ncols) + (r-1)*ncols
# }#r
# nf <- graphics::layout(mx,
#                        widths = c(rep(1,ncol(mx))),
#                        heights = rep(1,2))
# par(mar = c(0,0,0,0))
# 
# for(t in 1:length(years)){
#   ##-- Plot maps
#   plot( sf::st_geometry(COUNTRIES), border = NA, col = c("gray80","gray60"))
#   try(
#     plot( sf::st_geometry(alive[alive$Year == years[t], ]), add = TRUE, col = "orange", pch = 3),
#     silent = TRUE)
#   try(
#     plot( sf::st_geometry(dead.recovery[dead.recovery$Year == years[t], ]), add = TRUE, col = "slateblue", pch = 3),
#     silent = TRUE)
#   plot( sf::st_geometry(COUNTRIES), border = "gray40", col = NA, add = TRUE)
# 
#   ##-- Add year
#   graphics::mtext(text = years[t],
#                   side = 1, line = -18,
#                   adj = 0.18, cex = 1.2)
# }#t
# dev.off()
# 
# 
# 
# ## -----   6.5. PREVIOUSLY DETECTED - FIGURE -----
# 
# ##-- Plot number of individuals with previous NGS detections
# plot1 <- alive %>%
#   dplyr::group_by(Year, detected.earlier) %>%
#   dplyr::summarise(n = length(unique(Id))) %>%
#   ggplot2::ggplot() +
#   geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
#   labs( tag = "A",
#         x = "years",
#         y = "Number of individuals detected through NGS") +
#   theme( legend.position = "none",
#          legend.title = element_blank(),
#          axis.text.x = element_text(angle = 60,
#                                     hjust = 1)) +
#   scale_x_continuous(breaks = years, labels = years) +
#   scale_fill_manual(values = c("gray20", "gray60"))
# 
# ##-- Plot number of dead recoveries with previous detections
# plot2 <- dead.recovery %>%
#   dplyr::group_by(Year, detected.earlier) %>%
#   dplyr::summarise(n = length(unique(Id))) %>%
#   ggplot2::ggplot() +
#   ggplot2::geom_col(aes(x = Year, y = n, fill = detected.earlier)) +
#   labs( tag = "B",
#         x = "years",
#         y = "Number of individuals recovered dead") +
#   theme( legend.title = element_blank(),
#          legend.position.inside = 2,
#          axis.text.x = element_text(angle = 60,
#                                     hjust = 1)) +
#   scale_x_continuous(breaks = years, labels = years) +
#   scale_fill_manual(values = c("gray20", "gray60"))
# plot_total <- plot1 + plot2
# 
# #-- Save figure
# ggsave(filename = file.path(working.dir, "figures",
#                             paste0( engSpecies, "_Previous Detection_",
#                                     years[1]," to ", years[length(years)],
#                                     ".png")),
#        plot = plot_total,
#        dpi = 300, height = 6, width = 12, device = "png")
# # grDevices::png( filename = file.path( working.dir, "figures",
# #                                       paste0( engSpecies, "_Previous Detection_",
# #                                               years[1]," to ", years[length(years)],
# #                                               ".png")),
# #                 width = 12, height = 6,
# #                 units = "in", pointsize = 12,
# #                 res = 300, bg = NA)
# # dev.off()
# 
# 
# 
# ## ----- 7. SAVE DATA ------
# 
# save( alive,
#       dead.recovery,
#       file = file.path( working.dir, "data", fileName))

cleanRovbaseData(
  species = "wolf",
  data.dir = data.dir,
  working.dir = working.dir,
  two.sex = TRUE,
  print.report = TRUE)



##------------------------------------------------------------------------------
# ## ------ 3. FILTER DATA -----
# 
# ## ------   3.1. ALIVE DATA -----
# 
# data.alive <- myFullData.sp$alive %>%
#   dplyr::filter(
#     ##-- Subset to years of interest
#     Year %in% years,
#     ##-- Subset to months of interest
#     Month %in% unlist(sampling.months),
#     ##-- Subset to sex of interest
#     Sex %in% sex) %>%
#   ##-- Filter based on space 
#   sf::st_filter( .,st_as_sfc(st_bbox(studyArea)), .predicate = st_intersects)
# 
# 
# 
# ## ------   3.2. DEAD RECOVERY DATA -----
# 
# data.dead <- myFullData.sp$dead.recovery %>%
#   dplyr::filter(
#     ##-- Subset to years of interest
#     Year %in% years,
#     ##-- Subset to sex of interest
#     Sex %in% sex) %>%
#   ##-- Filter based on space 
#   sf::st_filter( .,st_as_sfc(st_bbox(studyArea)), .predicate = st_intersects)
# 
# # ## Remove all alive detections outside of the study area extent 
# # data.alive <- data.alive[!is.na(as.numeric(st_intersects(data.alive, st_as_sfc(myStudyArea.extent)))), ]
# # ## Remove all dead recoveries outside of the study area polygon 
# # data.dead <- data.dead[!is.na(as.numeric(st_intersects(data.dead, st_as_sfc(myStudyArea.extent)))), ]
# # ## Remove all alive detections outside of the sampling period 
# # data.alive <- data.alive[ data.alive$Month %in% unlist()&  data.alive$Year %in% unlist(DATA$years), ]
# # ## Remove all dead recoveries outside of the sampling period 
# # data.dead <- data.dead[ data.dead$Year %in% unlist(DATA$years),] 
# 
# 
# 
# ## ------   3.3. DATA SUMMARY ------
# 
# ## ------     3.3.1. PLOT CHECKS ------ 
# 
# ## PROPORTION SAMPLES STRUCTURED/OTHERS
# pdf(file = file.path(working.dir, "figures", "ProportionStucturedOther.pdf"))
# par(mfrow = c(1,1), mar = c(4,4,3,2))
# barplot(rbind(table(myFilteredData.spStructured$Year),
#               table(myFilteredData.spOthers$Year)),
#         beside = T, ylim = c(0,2000),
#         col = c(grey(0.2),grey(0.8)),
#         ylab = "Number of samples")
# abline(h = seq(0,2000,by=500), lty = 2, col = grey(0.8))
# title(main = "500m threshold")
# legend("topleft",
#        fill = c(grey(0.2),grey(0.8)),
#        legend = c("Structured","Other"))
# dev.off()
# 
# ## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
# tmp <- data.alive[data.alive$Proevetype %in% 
#                     c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
# tab <- table(tmp$Year, tmp$TrackRovbsID, useNA ="always" )
# 
# ## MAP  SAMPLES STRUCTURED OTHERS
# pdf(file = file.path(working.dir, "figures", "MapStucturedOther.pdf"))
# for(t in 1:n.years){
#   par(mar = c(0,0,3,0), mfrow = c(1,3))
#   tmp1 <- tmp[tmp$Year%in% years[t],]
#   tmpNoTracks <- tmp1[is.na(tmp1$TrackRovbsID), ]
#   tmpTracks <- tmp1[!is.na(tmp1$TrackRovbsID), ]
#   
#   plot(st_geometry(studyArea), main="Structured with track")
#   plot(st_geometry(tmpTracks), pch=21, col="black", cex=1,bg="red",add=T)
#   
#   plot(st_geometry(studyArea), main="Structured without track")
#   plot(st_geometry(tmpNoTracks), pch=21, col="black", cex=1,bg="blue",add=T)
#   
#   tmpOpp <- data.alive[!data.alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
#   tmpOpp <- tmpOpp[tmpOpp$Year%in% years[t],]
#   
#   plot(st_geometry(studyArea), main="Other samples")
#   plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
#   mtext(years[t],adj = -0.8,padj = 1)
# }
# barplot(tab[,which(is.na(colnames(tab)))]/rowSums(tab),main="% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 
# dev.off()
# 
# ## OVERALL MAP DETECTION DEAD RECOVERIES MAP
# pdf(file = file.path(working.dir, "figures", "OverallDetectionsDeadRecoveries.pdf"))
# plot(st_geometry(GLOBALMAP))
# plot(st_geometry(studyArea),add=T)
# plot(st_geometry(myFullData.sp$alive), pch=16, col="red", cex=0.3,add=T)
# plot(st_geometry(myFullData.sp$dead.recovery),pch=16, col="blue", cex=0.3,add=T)
# mtext(paste("Live detections", length(myFullData.sp$alive),
#             "; ID:", length(unique(myFullData.sp$alive$Id))),
#       line = +1)
# mtext(paste("Dead recovery:",length(myFullData.sp$dead.recovery)))
# dev.off()
# 
# 
# ## ------     3.3.1. DETECTIONS ALIVE ------ 
# 
# ## NUMBER OF INDIVIDUALS DETECTED ALIVE
# length(unique(data.alive$Id))
# 
# ## NUMBER OF INDIVIDUALS DETECTED ALIVE & RECOVERED DEAD
# sum(unique(data.alive$Id) %in% unique(data.dead$Id))
# 
# ## NUMBER OF INDIVIDUALS DETECTED/YEAR/COUNTRY
# table.id <- table(data.alive$Id,data.alive$Year, data.alive$Country)
# apply(table.id, c(2,3), function(x) sum(x>0))
# 
# ## NUMBER OF DETECTIONS/YEAR/COUNTRY
# countrytab <- table(data.alive$Year, data.alive$Country)
# countrytab 
# 
# ## NUMBER OF DETECTIONS/YEAR/COUNTRY/SEX
# sex_countrytab <- table(data.alive$Year, data.alive$Country, data.alive$Sex)
# sex_countrytab
# 
# ## PROPORTION OF DETECTIONS PER YEAR/COUNTRY/SEX
# sex_countrytab_prop <- sex_countrytab
# for(i in 1:dim(sex_countrytab)[3]){
#   sex_countrytab_prop[,,i] <- sex_countrytab_prop[,,i]/countrytab
# }
# sex_countrytab_prop
# 
# 
# 
# ## ------     3.3.2. DEAD RECOVERIES ------ 
# 
# ## NUMBER OF INDIVIDUALS RECOVERED
# length(unique(data.dead$Id))
# 
# ## NUMBER OF DEAD RECOVERIES/YEAR/COUNTRY
# table(data.dead$Year, data.dead$Country)
# 
# ## MORTALITY CAUSES
# unique(as.character(data.dead$DeathCause))
# unique(as.character(data.dead$DeathCause_2))  
# MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))
# 
# ## DEFINE LEGAL MORTALITY
# legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
# legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
# 
# 
# 
# ## ------     3.3.3. PLOTS ------ 
# 
# legal.death <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]
# Other.death <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]
# 
# table(legal.death$Year)
# table(Other.death$Year)
# 
# ## PLOT CHECK
# if(plot.check){
#   pdf(file = file.path(working.dir,"figures","DetectionsDeadRecoveries.pdf"))
#   
#   ## PLOT ALL DETECTIONS
#   for(t in 1:n.years){
#     plot(habitat$habitat.r, main = years[t]) 
#     plot(st_geometry(detectors$main.detector.sp), add=T, pch=16, cex=0.1)
#     plot(st_geometry(data.alive[data.alive$Year == years[t], ]), 
#          pch=16,col="red", cex=0.7, add=T)
#     mtext(paste("Live detections",
#                 nrow(data.alive[data.alive$Year == years[t], ]),
#                 "; ID:",
#                 nrow(unique(data.alive[data.alive$Year == years[t], ]$Id))),
#           line = +1)
#     
#     plot(st_geometry(data.dead[data.dead$Year == years[t], ]),
#          pch=16,col="blue", cex=0.7,add=T)
#     mtext(paste("Dead recovery:",nrow(data.dead[data.dead$Year == years[t], ])))
#     plot(st_geometry(GLOBALMAP), add=T) 
#   }#t
#   
#   ## PLOT NUMBER OF MORTALITY EVENTS, LEGAL/OTHERS 
#   barplot(table(legal.death$Country, legal.death$Year), main = "Legal causes")
#   barplot(table(Other.death$Country, Other.death$Year), main = "Other causes")
#   legend("topleft", fill = c(grey(0.3), grey(0.7)), legend = c("NOR","SWE"))
#   dev.off()  
# }
# 
# ## EXPORT NGS DATA 
# if(DATA$sex == "Hann"){
#   assign("myFilteredData.spM", myFilteredData.sp)
#   assign("myFilteredData.spOthersM", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
#   
#   save(myFilteredData.spM, 
#        myFullData.spM,
#        myFilteredData.spOthersM,
#        myFilteredData.spStructuredM,
#        file = file.path(working.dir, "data", "NGSData.RData"))
# } else {
#   assign("myFilteredData.spF", myFilteredData.sp)
#   assign("myFilteredData.spOthersF", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
#   
#   save(myFilteredData.spF,
#        myFullData.spF,
#        myFilteredData.spOthersF,
#        myFilteredData.spStructuredF,
#        file = file.path(working.dir, "data", "NGSData.RData"))
# }
# 
# 

##------------------------------------------------------------------------------
####----------------------- 
####---- makeRovQuantData() -----
# makeRovquantData_wolf <- function(
    # ##-- paths
# data.dir = getwd()
# working.dir = getwd()

##-- data
years = 2015:2024
sex = c("female","male")
aug.factor = 0.8
sampling.months = list(10:12,1:4)

##-- habitat
habitat.res = 20000
x.extent = NULL
y.extent = NULL
buffer.size = 40000
max.move.dist = 250000

##-- detectors
detector.res = 10000
subdetector.res = 1000
max.det.dist = 45000
resize.factor = 1

##-- Miscellanious
rename.list = NULL
#  ){
####----------------------- 

## ------ 0. BASIC SET-UP ------

##-- Set default values for the wolf model
if(is.null(aug.factor)){aug.factor <- 0.8}
if(is.null(sampling.months)){sampling.months <- list(10:12,1:4)}
if(is.null(habitat.res)){habitat.res <- 20000} 
if(is.null(x.extent)){x.extent <- c(210000,760000)}
if(is.null(y.extent)){y.extent <- c(6000000,7050000)}
if(is.null(buffer.size)){buffer.size <- 40000}
if(is.null(detector.res)){detector.res <- 10000}
if(is.null(subdetector.res)){subdetector.res <- 1000}
if(is.null(max.det.dist)){max.det.dist <- 45000}
if(is.null(resize.factor)){resize.factor <- 1}
if(is.null(rename.list)){rename.list = r.list.internal}

##-- Set up list of Habitat characteristics
habitat <- list( resolution = habitat.res,
                 buffer = buffer.size,
                 x.extent = x.extent,
                 y.extent = y.extent)

##-- Set up list of Detectors characteristics
detectors <- list( resolution = detector.res,
                   resolution.sub = subdetector.res,
                   maxDist = max.det.dist,
                   resize.factor = resize.factor)

##-- Set up list of Data characteristics
DATA <- list( sex = sex,
              aug.factor = aug.factor,
              sampling.months = sampling.months)



##------------------------------------------------------------------------------

## ------ I. LOAD AND SELECT DATA ------

## ------   1. HABITAT DATA -----

##-- Load pre-defined habitat rasters and shapefiles
#data(COUNTRIES, envir = environment()) 
#data(COUNTIES, envir = environment()) 
data(habitatRasters, envir = environment()) 
#data(GLOBALMAP, envir = environment()) 
data(REGIONS, envir = environment())

##-- Disaggregate habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)

##-- Merge counties for practical reasons
COUNTIES_AGGREGATED <- REGIONS %>%
  mutate(id = case_when(
    county %in% c("Akershus","Agder","Buskerud","Vestfold","Oslo","Østfold","Telemark") ~ "NO1",
    county %in% c("Innlandet","Møre og Romsdal","Trøndelag") ~ "NO2",
    county %in% c("Jämtland","Västernorrland","Västerbotten","Dalarna","Gävleborg") ~ "SE1",
    county %in% c("Uppsala","Västmanland","Stockholm","Södermanland") ~ "SE2",
    county %in% c("Blekinge","Örebro","Östergötland","Jönköping","Kronoberg","Kalmar","Skåne","Gotlands") ~ "SE3",
    county %in% c("Västra Götaland","Värmland","Halland") ~ "SE4")) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise()

if(plot.check){
  COUNTIES_AGGREGATED %>%
    ggplot(.) +
    geom_sf(aes(fill = id)) +
    geom_sf_label(aes(label = id)) +
    geom_sf(data = studyArea)
}



## ------   2. NGS DATA -----

##-- Extract date from the last cleaned data file
DATE <- getMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolf")

##-- Load the most recent clean wolf data from RovBase
myFullData.sp <- readMostRecent( 
  path = file.path(working.dir,"data"),
  pattern = "CleanData_wolf",
  extension = ".RData") 

##-- List years
if(is.null(years)){
  years <- sort(unique(c(myFullData.sp$alive$Year,
                         myFullData.sp$dead.recovery$Year)))
}
DATA$years <- years
n.years <- length(years)

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



##------------------------------------------------------------------------------

## ------ II. CREATE OPSCR DATA ------

## ------   1. GENERATE HABITAT ------

message("Preparing habitat characteristics... ")

## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------

# [PD]:REMOVE
# ## CREATE STUDY AREA POLYGON BASED x AND y EXTENTS
# myStudyArea.extent  <- st_bbox(extent(x.extent, y.extent))
# st_crs(myStudyArea.extent) <- st_crs(COUNTRIESWaterHumans)
# myStudyArea.poly <- st_crop(COUNTRIESWaterHumans,
#                             extent(x.extent, y.extent))
# # to get only "polygons objects"
# myStudyArea.poly <- st_collection_extract(myStudyArea.poly, "POLYGON")
# 
# myHabitat.list <- MakeHabitatFromRastersf( 
#   poly = myStudyArea.poly,
#   habitat.r = habitatRasters[["Habitat"]],
#   buffer = habitat$buffer,                               
#   plot.check = T)

##-- Determine study area based on predefined extent
studyArea <- st_crop( REGIONS,
                      xmin = x.extent[1], xmax = x.extent[2],
                      ymin = y.extent[1], ymax = y.extent[2]) %>%
  st_collection_extract(., "POLYGON") %>%
  summarise()

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster(
  poly = studyArea,
  habitat.r = habRaster,
  buffer = habitat$buffer,
  plot.check = TRUE) %>%
  append(habitat,.)

##-- Retrieve number of habitat windows 
isHab <- habitat$habitat.r[] == 1
n.habwindows <- habitat$n.habwindows <- sum(isHab)
habitat$habitat.df <- cbind.data.frame(
  "id" = 1:habitat$n.habwindows,
  "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
  "y" = raster::coordinates(habitat$habitat.r)[isHab,2])

##-- Make a spatial grid from polygon
habitat$grid <- sf::st_as_sf( stars::st_as_stars(habitat$habitat.r), 
                              as_points = FALSE,
                              merge = FALSE) %>%
  filter( Habitat %in% 1) %>%
  mutate( id = 1:nrow(.)) %>%
  st_set_crs( ., value = sf::st_crs(habitat$buffered.habitat.poly))

## [PD] CHECK IF NECESSARY
##-- Study area grid from habitat raster
habitat.rWthBufferPol <- sf::st_as_sf( 
  stars::st_as_stars(habitat$habitat.rWthBuffer), 
  as_points = FALSE,
  merge = TRUE) %>%
  dplyr::filter(Habitat %in% 1)



## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       1.2.1. DENSITY OF PACKS/PAIRS ------

##-- Kernel of NGS detections of individuals in pairs 
kern <- list()
habDens <- matrix(NA, nrow = n.habwindows, ncol = n.years)
for(t in 1:n.years){
  ##-- Subset the NGS data to individuals in packs/pairs this year
  data.pairs.t <- myFullData.sp$alive %>%
    dplyr::filter( Year == years[t],
                   Status %in% c("Pair", "Family group"))
  
  ##-- Get mean coordinates of pack coordinates
  IDs <- unique(data.pairs.t$IdSimplified)
  m.xy <- matrix(NA, nrow = length(IDs), ncol = 2)
  colnames(m.xy) <- c("x","y")
  for(i in 1:length(IDs)){
    m.xy[i, ] <- data.pairs.t %>%
      dplyr::filter( IdSimplified == IDs[i]) %>%
      st_coordinates(.) %>%
      colMeans(.)
  }#i
  
  ##-- Check if some coordinates are missing  
  if(sum(is.na(m.xy[ ,1])) > 0){m.xy <- m.xy[!is.na(m.xy[ ,1]), ]}
  
  ##-- Turn into .sf
  locationsFamily <- st_as_sf( as.data.frame(m.xy),
                               coords = c("x","y"),
                               crs = st_crs(habitat$habitat.sp))
  locationsFamily$id <- rep(1, nrow(locationsFamily))
  
  ##-- Calculate kernel of detections
  kern[[t]] <- raster(estUDm2spixdf(kernelUD( 
    as(locationsFamily[ ,"id"], "Spatial"), 
    h = 15000,
    grid = as(habitat$habitat.r, 'SpatialPixels'))))
  
  ##-- Plot check
  plot(kern[[t]], main = years[t])
  plot(habitat$habitat.poly$geometry, add = T, col = NA)
  
  habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ] == 1])
} #t



## ------   2. GENERATE DETECTORS ------ 

message("Preparing detectors characteristics... ")

## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----

##-- Generate NGS detectors based on the study area 
detectors$subdetectors.r <- raster::disaggregate(
  habitat$habitat.rWthBuffer,
  fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)

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
         Detector = rep(1,nrow(detectors$main.detector.sp))))

##-- Make a spatial grid from detector raster
detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
  x = detectors$raster,
  fun = function(x){x>0})) %>%
  sf::st_set_crs(.,value = sf::st_crs(studyArea)) %>%
  dplyr::mutate( id = 1:nrow(.)) 

##-- Extract numbers of detectors
n.detectors <- detectors$n.detectors <- dim(detectors$main.detector.sp)[1]



## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES ------ 

## ------       2.2.1. EXTRACT COUNTRIES ------ 

##-- Extract closest country for each detector
detCountries <- detectors$main.detector.sp %>%
  sf::st_distance(., COUNTRIES, by_element = F) %>%
  apply(., 1, function(x) which.min(x)) %>%
  as.factor(.) %>%
  as.numeric(.)

##-- Put into "nimble2SCR" format
detectors$detectors.df$countries <- detCountries



## ------       2.2.2. EXTRACT COUNTIES ------ 

##-- Extract closest county for each detector
detCounties <- detectors$main.detector.sp %>%
  sf::st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
  apply(., 1, function(x) which.min(x))

##-- Put into "nimble2SCR" shape
detectors$detectors.df$counties <- detCounties



## ------       2.2.3. EXTRACT GPS TRACKS LENGTHS ------ 

message("Cleaning GPS tracks... ")

##-- Combine all GPS tracks
# TRACKS2 <- rbind(
#   sf::read_sf(file.path(data.dir, "Tracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250422.shp")),
#   sf::read_sf(file.path(data.dir, "Tracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422.shp"))) %>%
#   ##-- Process dates
#   dplyr::mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
#                  Mth = as.numeric(format(Dato,"%m")),
#                  Yr = as.numeric(format(Dato,"%Y")),
#                  Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
#   ##-- Filter out irrelevant tracks
#   dplyr::filter( Helikopter == "0",      ## Remove helicopter tracks
#                  # Jerv == "1",          ## [CHECK] should we keep wolf tracks only?
#                  Year %in% years & Mth %in% unlist(sampling.months)) %>% ## Keep tracks during sampling season only
#   ##-- Extract track lengths & centroids
#   dplyr::mutate( Length = sf::st_length(., byid = T),
#                  Centroidx = sf::st_coordinates(sf::st_centroid(.))[ ,1])
# 
# ##-- Find & filter out duplicates based on person, distance and date.
# df <- data.frame( Dato = TRACKS2$Dato,
#                   Year = TRACKS2$Year,
#                   Person = TRACKS2$Person,
#                   Length = TRACKS2$Length,
#                   Centroidx = TRACKS2$Centroidx)
# dupIDs <- which(duplicated(df))
# dupLength <- TRACKS2$Length[duplicated(df)]
# TRACKS2 <- TRACKS2[-dupIDs, ]
TRACKS <- readTracks( data.dir = data.dir,
                      years = years,
                      sampling.months = sampling.months)

##-- Extract length of GPS search track per detector grid cell
detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
##-- Set-up progress bar
pb = utils::txtProgressBar( min = 1, max = n.years, initial = 0, style = 3)
# TRACKS.r <- list()
for(t in 1:n.years){
  intersection <- TRACKS %>%
    dplyr::filter(Year == years[t]) %>%
    sf::st_intersection(detectors$grid, .) %>%
    dplyr::mutate(LEN = st_length(.)) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(transect_L = sum(LEN)) ##-- Get total length searched in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  # TRACKS.r[[t]] <- detectors$raster
  # TRACKS.r[[t]][detectors$raster[] %in% 1] <- detTracks[ ,t]
  ##-- Print progress 
  utils::setTxtProgressBar(pb,t) 
  }#t

##-- Put into "nimble2SCR" format
colnames(detTracks) <- paste0("tracks.", years)
detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detTracks)



## ------       2.2.4. EXTRACT DISTANCES TO ROADS ------ 

##-- Load map of distance to roads (1km resolution)
DistAllRoads <- readMostRecent( path = file.path(data.dir, "Roads"), 
                                extension = ".tif", 
                                stack = FALSE)

##-- Fasterize to remove values that fall in the sea
r <- fasterize::fasterize(sf::st_as_sf(REGIONS), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- raster::crop(DistAllRoads, studyArea)
rm(list = c("r"))

##-- Aggregate to match the detectors resolution 
DistAllRoads <- raster::aggregate( 
  x = DistAllRoads,
  fact = detectors$resolution/raster::res(DistAllRoads),
  fun = mean)

##-- Extract distance to roads for each detector
detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)

##-- If NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( x = DistAllRoads,
                        y = detectors$main.detector.sp[isna, ],
                        buffer = 15000,
                        fun = mean,
                        na.rm = T)
detRoads[isna] <- tmp
detRoads <- round(scale(detRoads), digits = 2)

##-- Put into "nimble2SCR" format
detectors$detectors.df$roads <- detRoads



## ------       2.2.5. EXTRACT DAYS OF SNOW ------ 

##-- Load raster stack of snow cover
SNOW <- readMostRecent( path = file.path(data.dir, "Snow"), 
                        extension = ".tif", 
                        stack = TRUE)

##-- Select snow data corresponding to the monitoring period
SNOW <- SNOW[[paste("X", years, "_", years + 1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))

##-- Extract snow 
detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = n.years)
det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:n.years] <- raster::extract(SNOW, det.sptransf)

##-- if NA returns the average value of the cells within 15km 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp

##-- still some NA... Increase buffer again 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 35000, fun = mean, na.rm = T)
detSnow[isna,1:n.years] <- tmp

##-- Put into "nimble2SCR" format
colnames(detSnow) <- paste0("snow.", years)
detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detSnow)



## ------       2.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------ 

## ------         2.2.6.1. SKANDOBS ------ 

##-- Load the last SkandObs data file
skandObs <- readMostRecent( 
  path = file.path(data.dir,"Skandobs"),
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
                                             year-1, year)) %>%
  ##-- Filter based on monitoring season
  dplyr::filter( month %in% unlist(sampling.months)) %>%
  ##-- Turn into spatial points object
  sf::st_as_sf(., coords = c("longitude","latitude")) %>%
  sf::st_set_crs(., value = "EPSG:4326") %>%
  sf::st_transform(., sf::st_crs(COUNTIES)) %>%
  sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)

# ## RASTERIZE AT THE DETECTOR LEVEL
# r.detector <- aggregate(habitat.subdetectors, fact=(detector.res/subdetector.res))
# r.list <- lapply(years, function(y){
#   rl <- raster::rasterize(skandObs[skandObs$monitoring.season %in% y, 1], r.detector , fun="count")[[1]]
#   rl[is.na(rl[])] <- 0
#   rl[!r.detector[]%in% 1] <- NA
#   rl1 <- rl
#   rl1[rl[]>0] <- 1
#   list(rl1, rl)
# })
# r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
# r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))
# plot(r.skandObsSamplesBinary[[t]])
# 
# ## PLOT CHECK 
# if(plot.check){
#   ## SUMMARY SKANDOBS
#   pdf( file = file.path(working.dir,"figures","skandObs.pdf"),
#        width = 10)
#   barplot(table(skandObs$monitoring.season))
#   barplot(table(skandObs$month), xlab = "Months")
#   barplot(table(skandObs$activity), cex.names = 0.7)
#   barplot(table(skandObs$species))
#   
#   ## MAPS 
#   par(mar = c(0,0,2,0))
#   for(t in 1:n.years){
#     plot(st_geometry(studyArea), main= years[t])
#     plot(st_geometry(skandObs[skandObs$monitoring.season %in% years[t],  ]), pch=16, col="red", cex=0.1)
#   }
#   dev.off()
# }



## ------         2.2.6.2. ROVBASE ------ 

##-- Process Rovbase observations (all species)
rovbaseObs <- readMultiples( 
  path = file.path(data.dir, "AllSamples"),
  extension = ".xlsx") %>%
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
    monitoring.season = ifelse(month < unlist(sampling.months)[1],
                               year-1, year)) %>%
  ##-- Filter out unusable samples
  dplyr::filter( 
    ##-- Filter out samples without coordinates,...
    !is.na(East_UTM33),
    ##-- ...based on species
    Species %in% c("Ulv"),
    ##-- ...based on sample type
    Sample_type %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
                        "Saliv/Spytt", "Loepeblod", "Vev"),
    ##-- ...based on monitoring season
    month %in% unlist(sampling.months)
    ##-- ... if sample was from the focal species and successfully genotyped 
    # !(Species %in% "Ulv" & !is.na(Id))
    ) %>% # NOT ANYMORE
  ##-- Turn into spatial points object
  sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(. , sf::st_crs(COUNTIES)) %>%
  ##-- Filter based on space 
  sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)

# ## RASTERIZE 
# r.detector <- aggregate(habitat.subdetectors, fact=(detector.res/subdetector.res))
# r.list <- lapply(years, function(y){
#   rl <- raster::rasterize(rovbaseObs.sp[rovbaseObs.sp$monitoring.season %in% y, 1], r.detector , fun="count")[[1]]
#   rl[is.na(rl[])] <- 0
#   rl[!r.detector[]%in% 1] <- NA
#   rl1 <- rl
#   rl1[rl[]>0] <- 1
#   list(rl1, rl)
# })
# r.OtherSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
# r.OtherSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))
# 
# ## PLOT CHECK
# if(plot.check){
#   
#   pdf(file = file.path(working.dir,"figures","mapStructuredOthers.pdf"))
#   for(t in 1:n.years){
#     year = years[t]
#     tmpOthers <- myFilteredData.spOthers[myFilteredData.spOthers$Year%in%year, ]
#     tmpStruct <- myFilteredData.spStructured[myFilteredData.spStructured$Year%in%year, ]
#     
#     par(mfrow = c(2,2), mar = c(0,0,5,0))
#     plot( r.OtherSamplesBinary[[t]],
#           main = paste(year,"\n Rovbase Samples Other"), box = F, axes = F)
#     plot( st_geometry(tmpOthers),
#           pch = 16, col = "blue", bg = "blue", cex = 0.6, add = T)
#     plot( r.OtherSamplesBinary[[t]],
#           main = paste (year,"\n Rovbase Samples Structured"), box = F, axes = F)
#     plot( st_geometry(tmpStruct),
#           pch = 16, col = "red", bg = "red", cex = 0.6, add = T)
#     
#     plot( r.skandObsSamplesBinary[[t]],
#           main = paste(year,"\n SkandObs Other"),
#           box = F, axes = F)
#     plot( st_geometry(tmpOthers),
#           pch = 16, col = "blue", bg = "blue", cex = 0.6, add = T)
#     plot( r.skandObsSamplesBinary[[t]],
#           main = paste(year,"\n SkandObs Structured"), box = F, axes = F)
#     plot( st_geometry(tmpStruct),
#           pch = 16, col = "red", bg = "red", cex = 0.5, add = T)
#   }
#   dev.off()
# }



## ------         2.2.6.3. COMBINE ROVBASE & SKANDOBS ------ 

##-- Rasterize at the detector level
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
  ##-- Store in a list
  list(sk.r, sk.r1, rb.r, rb.r1)
})

##-- Store in raster bricks
r.skandObsContinuous <- brick(lapply(r.list,function(x) x[[1]]))
r.skandObsBinary <- brick(lapply(r.list,function(x) x[[2]]))
r.rovbaseContinuous <- brick(lapply(r.list,function(x) x[[3]]))
r.rovbaseBinary <- brick(lapply(r.list,function(x) x[[4]]))

##-- Combine both rasters
r.SkandObsRovbaseBinary <- r.rovbaseBinary + r.skandObsBinary
for(t in 1:n.years){
  r.SkandObsRovbaseBinary[[t]][r.SkandObsRovbaseBinary[[t]][] > 1] <- 1
}#t



## ------         2.2.6.4. ASSIGN THE COVARIATE ------ 

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
detOtherSamples[ ,1:n.years] <- raster::extract( r.SkandObsRovbaseBinary,
                                                 detectors$main.detector.sp)

##-- Put into "nimble2SCR" format
colnames(detOtherSamples) <- paste0("detOtherSamples.", years)
detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detOtherSamples)



## ------       2.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

##-- STRUCTURED 
detCovs <- array(NA, c(dim(detTracks)[1], 2, dim(detTracks)[2]))
detCovs[ ,1, ] <- detTracks
detCovs[ ,2, ] <- detSnow
dimnames(detCovs) <- list( "detectors" = 1:n.detectors,
                           "covariates" = c("tracks", "snow"), 
                           "years" = years)

##-- OTHERS 
detCovsOth <- array(NA, c(dim(detTracks)[1], 3, dim(detTracks)[2]))
detCovsOth[ ,1, ] <- detRoads
detCovsOth[ ,2, ] <- detSnow
detCovsOth[ ,3, ] <- detOtherSamples
dimnames(detCovsOth) <- list( "detectors" = 1:n.detectors,
                              "covariates" = c( "roads","snow","obs"),
                              "years" = years)

##-- Store in the detector list                          
detectors$covariates <- detCovs
detectors$covariates.others <- detCovsOth

##-- Merge with the detector grid
detectors$grid <- dplyr::left_join(
  x = detectors$grid,
  y = detectors$detectors.df,
  by = "id")



## ------   3. RESCALE COORDINATES -----

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

# ##-- Create habitat grid
# habitatGrid <- matrix(0,
#                       nrow = max(habitat$scaledCoords[ ,"y"]) + 1,
#                       ncol = max(habitat$scaledCoords[ ,"x"]) + 1)
# for (c in 1:nrow(habitat$scaledCoords)) {
#   habitatGrid[trunc(habitat$scaledCoords[c,"y"]) + 1,
#               trunc(habitat$scaledCoords[c,"x"]) + 1] <- c
# }#c



## ------   4. CREATE LOCAL OBJECTS -----

# localDets <- getLocalObjects(
#   habitatMask = habitat$habitat.mx,
#   coords = scaledCoords$coordsDataScaled,
#   dmax =  DETECTIONS$maxDist*1.4/res(habitat$habitat.r)[1],
#   resizeFactor = 1,
#   plot.check = TRUE)

##-- Get local detectors
detectors$localObjects <- getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = detectors$scaledCoords,
  dmax = detectors$maxDist*1.4/habitat$resolution,
  resizeFactor = detectors$resize.factor,
  plot.check = F)



## ------   5. SAVE STATE-SPACE CHARACTERISTICS -----

save( habitat,
      file = file.path( working.dir, "data",
                        paste0("Habitat_wolf_", DATE, ".RData")))

save( detectors,
      file = file.path( working.dir,"data",
                        paste0("Detectors_wolf_", DATE, ".RData")))



## ------   6. FILTER DATA -----

## ------     6.1. ALIVE DATA -----

data.alive <- myFullData.sp$alive %>%
  dplyr::filter(
    ##-- Subset to years of interest
    Year %in% years,
    ##-- Subset to months of interest
    Month %in% unlist(sampling.months),
    ##-- Subset to sex of interest
    Sex %in% sex) %>%
  ##-- Filter based on space 
  sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)



## ------     6.2. DEAD RECOVERY DATA -----

data.dead <- myFullData.sp$dead.recovery %>%
  dplyr::filter(
    ##-- Subset to years of interest
    Year %in% years,
    ##-- Subset to sex of interest
    Sex %in% sex) %>%
  ##-- Filter based on space 
  sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)



## ------     6.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

## ------       6.3.1. ASSIGN SAMPLES TO GPS SEARCH TRACKS ------

message("Assigning DNA samples to GPS tracks... ")
message("This can take several minutes... ")

data.alive <- assignSearchTracks(
  data = data.alive,
  tracks = TRACKS)

# ##-- ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
# data.alive$TrackRovbsID <- NA
# data.alive$TrackDist <- NA
# 
# TRACKSSimple_sf <- list()
# for(t in 1:n.years){
#   TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
#   TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbaseID)
#   TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]
# }
# 
# ## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
# dnatemp <- st_as_sf(data.alive)
# ## CREATE A BUFFER AROUND EACH DETECTION
# tmp <-  st_buffer(dnatemp, dist=750)
# 
# for(i in 1:nrow(data.alive)){
#   # INTERSECT POINT WITH TRACKS,
#   t <- which(years %in% tmp[i, ]$Year)
#   whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(data.alive$Date[i]))
#   tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])
#   
#   if(nrow(tmpTRACKS)==0){next}
#   
#   # FIND THE CLOSEST TRACK
#   dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
#   
#   # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
#   # IF NO MATCHING DATE ASSIGN TO NA?
#   if(length(dist)==0){
#     data.alive$TrackRovbsID[i] <- NA
#     data.alive$TrackDist[i] <- NA
#   }
#   # IF MATCHING DATE ASSING TO THAT TRACK
#   if(length(dist)==1){
#     data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
#     data.alive$TrackDist[i] <- dist
#   }
#   # IF SEVERAL MATCHING DATES ASSING TO THE CLOSEST OF THE MATCHING TRACKS
#   if(length(dist)>1){
#     data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
#     data.alive$TrackDist[i] <- min(dist)
#   }
#   print(i)
# }
# 
# ##-- SAVE FOR FASTER LOADING
# save( myFilteredData.sp,
#       file = file.path(working.dir, "data", "_myFilteredData.sp.RData"))
# load(file.path(working.dir, "data", "_myFilteredData.sp.RData"))



## ------       6.3.2. ASSIGN SAMPLES TO OPPORTUNISTIC OR STRUCTURED ------

distanceThreshold <- 500

##-- Identify samples from structured and opportunistic sampling
data.alive <- data.alive %>%
  dplyr::mutate(
    ##-- Collector column was replaced by two columns, merging them now...
    Collector_role = ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
    ##-- Identify samples collected during structured sampling 
    structured = Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
      !is.na(trackID) &
      trackDist <= distanceThreshold)

# distanceThreshold <- 500
# 
# ##-- Proevetype columns was replaced by two columns, merging them now...
# data.alive$Proevetype <-  ifelse(
#   data.alive$Annen.innsamler...Rolle %in% "", 
#   data.alive$Samlet.selv...Rolle,
#   data.alive$Annen.innsamler...Rolle)
# 
# whichStructured <- data.alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
#   !is.na(data.alive$TrackRovbsID) &
#   data.alive$TrackDist <= distanceThreshold
# myFilteredData.spStructured <- data.alive[whichStructured,]
# myFilteredData.spOthers <- data.alive[!whichStructured,]
# 
# ##-- CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
# nrow(myFilteredData.spStructured) + nrow(myFilteredData.spOthers)
# nrow(data.alive)
# 
# ##-- Check number of opp vs. struc each year 
# data.alive$TrackDistCat <- ifelse(data.alive$TrackDist > 500, 0, 1)
# data.alive$TrackDistCat[is.na(data.alive$TrackDistCat)] <- 0
# 
# table( data.alive$Proevetype,
#        data.alive$Year,
#        data.alive$TrackDistCat)



# ## ------     6.4. SEPARATE MORTALITY CAUSES ------
# 
# ##-- Identify legal mortality causes
# MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
# whichLegalCauses <- unlist(lapply(c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske"),
#                                   function(x)grep(x,MortalityNames)))
# legalCauses <- MortalityNames[whichLegalCauses]
# 
# ##-- Identify legal dead recoveries based on mortality causes
# data.dead <- data.dead %>%
#   mutate(legal = Death_cause %in% legalCauses)



## ------     6.4. ASSIGN DETECTORS ------ 

##-- ALL SAMPLES
data.alive <- assignDetectors( 
  data = data.alive,                
  detectors = detectors$main.detector.sp,
  subDetectors = detectors$detector.sp,
  radius = detectors$resolution)

##-- DEAD RECOVERY
data.dead <- assignDetectors( 
  data = data.dead,
  detectors = detectors$main.detector.sp,
  radius = detectors$resolution)

# ##-- ALL SAMPLES
# myData.alive <- AssignDetectors_v3sf( 
#   myData = data.alive,                
#   myDetectors = detectors$main.detector.sp,
#   mysubDetectors = mydetector.sp,
#   radius = detector.res)
# 
# ## STRUCTURED
# myData.aliveStruc <- AssignDetectors_v3sf( 
#   myData = myFilteredData.spStructured,                
#   myDetectors = detectors$main.detector.sp,
#   mysubDetectors = mydetector.sp,
#   radius = detector.res)
# 
# ## OTHERS
# myData.aliveOthers <- AssignDetectors_v3sf( 
#   myData = myFilteredData.spOthers,                
#   myDetectors = detectors$main.detector.sp,
#   mysubDetectors = mydetector.sp,
#   radius = detector.res)
# 
# ## DEAD RECOVERIES
# myData.dead <- AssignDetectors_v3sf(
#   myData = data.dead,
#   myDetectors = detectors$main.detector.sp,
#   radius = detector.res)



## ------     6.5. PLOT NGS and DEAD RECOVERY MAPS ----- 

##-- layout
L <- n.years
if(L < 6){ nrows <- 1 } else{
  if(L < 13){ nrows <- 2 } else {
    if(L < 22){ nrows <- 3 } else {
      if(L < 33){ nrows <- 4 } else {
        nrows <- 5
      }}}}
ncols <- ceiling(L/nrows)


##-- NGS maps
grDevices::png(filename = file.path(working.dir, "figures/NGS_TimeSeries.png"),
               width = ncols*2, height = nrows*4,
               units = "in", pointsize = 12,
               res = 300, bg = NA)
##-- layout
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
    plot( sf::st_geometry(data.alive$data.sp[data.alive$data.sp$Year == years[t], ]), add = TRUE, col = "orange", pch = 3),
    silent = TRUE)
  plot( sf::st_geometry(COUNTRIES), border = "gray20", col = NA, add = TRUE)
  
  ##-- Add year
  graphics::mtext(text = years[t],
                  side = 1, line = -18,
                  adj = 0.18, cex = 1.2)
}#t
dev.off()


##-- Dead recoveries maps
grDevices::png(filename = file.path(working.dir, "figures/DEAD_TimeSeries.png"),
               width = ncols*2, height = nrows*4,
               units = "in", pointsize = 12,
               res = 300, bg = NA)
##-- layout
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
  try( plot( sf::st_geometry(data.dead[data.dead$Year == years[t] & 
                                         data.dead$Legal, ]),
             add = TRUE, 
             col = "slateblue1",
             pch = 3),
       silent = TRUE)
  try( plot( sf::st_geometry(data.dead[data.dead$Year == years[t], ]),
             add = TRUE,
             col = "slateblue4",
             pch = 3),
       silent = TRUE)
  plot( sf::st_geometry(COUNTRIES),
        border = "gray20",
        col = NA,
        add = TRUE)
  
  ##-- Add year
  graphics::mtext(text = years[t],
                  side = 1, line = -18,
                  adj = 0.18, cex = 1.2)
}#t
dev.off()



## ------     6.6. SAVE FILTERED DATA ----- 

save( data.alive, data.dead,
      file = file.path( working.dir, "data",
                        paste0("FilteredData_wolf_", DATE, ".RData")))



## ------   7. GENERATE DETECTION HISTORY ------

for(thisSex in sex){
  
  message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
  
  ## ------     7.1. FILTER DATA BY SEX -----
  
  load(file.path( working.dir, "data",
                  paste0("FilteredData_wolf_", DATE, ".RData")))
  
  data.alive$data.sp <- data.alive$data.sp %>%
    dplyr::filter(Sex %in% thisSex)
  
  data.dead <- data.dead %>%
    dplyr::filter(Sex %in% thisSex)
  
  
  
  ## ------     7.2. GENERATE DETECTION HISTORY ARRAYS -----
  
  ##-- ALL SAMPLES
  y.ar <- makeY( data = data.alive$data.sp,
                 detectors = detectors$main.detector.sp,
                 method = "Binomial",
                 data2 = data.dead,
                 detectors2 = detectors$main.detector.sp,
                 returnIdvector = TRUE)
  
  ##-- STRUCTURED
  y.arStruc <- makeY( data = data.alive$data.sp[data.alive$data.sp$structured, ],
                      detectors = detectors$main.detector.sp,
                      method = "Binomial",
                      data2 = data.dead,
                      detectors2 = detectors$main.detector.sp,
                      returnIdvector = TRUE)
  
  ##-- OTHERS
  y.arOth <- makeY( data = data.alive$data.sp[!data.alive$data.sp$structured, ],
                    detectors = detectors$main.detector.sp,
                    method = "Binomial",
                    data2 = data.dead,
                    detectors2 = detectors$main.detector.sp,
                    returnIdvector = TRUE)
  
  ##-- Make sure all detection arrays have the same dimensions
  y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- array( 0, 
                                                     dim = dim(y.ar$y.ar),
                                                     dimnames = dimnames(y.ar$y.ar))
  ##-- Fill in the y arrays
  y.ar.ALIVEOthers[dimnames(y.arOth$y.ar)[[1]], , ] <- y.arOth$y.ar
  y.ar.ALIVEStructured[dimnames(y.arStruc$y.ar)[[1]], , ] <- y.arStruc$y.ar
  
  ##-- Project death to the next occasion
  y.ar.DEADProjected <- y.ar$y.ar2 
  y.ar.DEADProjected[] <- 0
  for(t in 2:n.years){ y.ar.DEADProjected[ , ,t] <- y.ar$y.ar2[ , ,t-1] }
  
  ##-- Get dead recovery detector index
  y.ar.DEAD <- apply( y.ar.DEADProjected,
                      c(1,3),
                      function(x){
                        if(sum(x)>0){which(x>0)}else{0}
                      })
  dimnames(y.ar.DEAD) <- list( "id" = dimnames(y.ar$y.ar2)[[1]],
                               "year" = dimnames(y.ar$y.ar2)[[3]])
  
  ##-- Create binary dead recovery histories (0: not recovered ; 1: recovered dead)
  y.ar.DEAD <- apply(y.ar.DEADProjected, c(1,3), function(x){as.numeric(sum(x)>0)})
  dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
  y.ar.DEAD[y.ar.DEAD > 0] <- 1
  
  
  
  # ## ------     7.3. GENERATE OBSERVATIONS : y.obs[i,t] ------ 
  #
  # INDIVIDUAL_ID$STATUS_Numeric <- 1
  # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% "Juvenile"] <- 2
  # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Pair" )] <- 3
  # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Family group")] <- 4
  # 
  # indIDRovBase <- unlist(lapply(strsplit(y.ar$Id.vector," "), function(x) x[1]))
  # ALLIDS <- c(unique(indIDRovBase))
  # y.obsALL <- matrix(1, nrow = length(ALLIDS), ncol = dim(y.ar.ALIVE)[3]+1)
  # yrs <- c(years[1]-1, years)
  # dimnames(y.obsALL) <- list(ALLIDS, yrs)
  # for(i in 1:dim(y.obsALL)[1]){
  #   for(t in 1:(n.years+1)){
  #     tmp <- unique(INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$ReprodYear..May.1.year.y...Apr.30.y.1. == yrs[t] & 
  #                                                  INDIVIDUAL_ID$ROVBASE_IndividID == ALLIDS[i]])
  #     if(length(tmp) > 0){
  #       if(length(tmp) > 1){ tmp <- tmp[1] }
  #       y.obsALL[i,t] <- tmp
  #     }
  #   }#t
  #   
  #   ##-- For the last years, get through the pair-based files from Linn.
  #   ## Linn's 2023 file
  #   tmp <-  Pack_ID2023[Pack_ID2023$Rovbase.ID %in% ALLIDS[i],]
  #   if(nrow(tmp) > 0){ y.obsALL[i,n.years-1] <- 3 }
  #   ## Linn's 2024 file
  #   tmp <-  Pack_ID2024[Pack_ID2024$RovbaseID %in% ALLIDS[i],]
  #   if(nrow(tmp) > 0){ y.obsALL[i,n.years] <- 3 }
  #   ## Linn's 2025 file
  #   tmp <-  Pack_ID2025[Pack_ID2025$IndividID %in% ALLIDS[i],]
  #   if(nrow(tmp) > 0){ y.obsALL[i,n.years+1] <- 3 }
  # }#i
  # 
  # ##-- subset y.obs for the individuals present in y.ar
  # indID <- unlist(lapply(strsplit(y.ar$Id.vector, " "), function(x)x[1]))
  # y.obs <- y.obsALL[indID, ]
  # y.obs[y.obs == 4] <- 3
  # y.obs <- y.obs[ ,as.character(years)]
  
  
  
  ## ------     7.3. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------ 
  
  distances <- list()
  for(t in 1:n.years){
    
    print(paste("------ ", t ," -------", sep = "" ))
    distances[[t]] <- checkDistanceDetections(
      y = y.ar$y.ar[ , ,t], 
      detector.xy = detectors$detectors.df[ ,c("x","y")], 
      max.distance = detectors$maxDist,
      method = "pairwise",
      plot.check = F)
    
    ##-- If some detections are flagged
    if(sum(distances[[t]]$y.flagged) > 0){
      
      ##-- Remove detections that are further than the threshold
      #y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
      y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
      y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
      
      ##-- Remove detections also in data.alive$data.sp to run getSInits later
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      idd <- names(affected.ids)
      for(i in 1:length(idd)){
        detIds <- which(distances[[t]]$y.flagged[idd[i], ] > 0)
        data.alive$data.sp <- data.alive$data.sp %>%
          dplyr::filter(!(Id %in% idd[i] & Detector %in% detIds & Year %in% years[t]))
      }#i
    }#if
  }#t
  
  
  
  ## ------     7.4. GENERATE INDIVIDUAL-LEVEL COVARIATES ------ 
  
  ## ------       7.4.1. INDIVIDUAL STATE ------ 
  
  ##-- Turn status into factor for correct order
  data.alive$data.sp$Status <- factor( data.alive$data.sp$Status,
                                       levels = c("Juvenile", "Pair", "Family group"))
  
  ##-- Make a table of individuals states per year
  ##-- 1: no info
  ##-- 2: Juvenile
  ##-- 3: Pair
  ##-- 4: Family group
  tmp <- apply( table(data.alive$data.sp$Id,
                      data.alive$data.sp$Status,
                      data.alive$data.sp$Year),
                c(1,3),
                function(x)ifelse(any(x>0), which(x>0), 0)) + 1
  
  ##-- Resize to match detection array
  y.obs <- y.status <- matrix(1, nrow = nrow(y.ar.DEAD), ncol = ncol(y.ar.DEAD))
  dimnames(y.obs) <- dimnames(y.ar.DEAD)
  y.obs[dimnames(tmp)[[1]], ] <- tmp
  y.obs[y.obs > 3] <- 3
  
  ##-- For the model, set status to 2 from the first time it is "pair" or "family group" to the last occasion
  for(i in 1:dim(y.status)[1]){
    if(any(y.obs[i, ] >= 3)){
      y.status[i, min(which(y.obs[i, ] >= 3)):ncol(y.status)] <- 2
    }
  }#i
  
  
  
  ## ------       7.4.2. TRAP-RESPONSE ------ 
  
  ##-- Make matrix of previous capture indicator
  detResponse <- makeTrapResponseCov(
    data = myFullData.sp$alive,
    data.dead = myFullData.sp$dead.recovery)
  
  ##-- Subset to focal years
  detResponse <- detResponse[ ,dimnames(detResponse)[[2]] %in% dimnames(y.ar$y.ar)[[3]]]
  
  ##-- Subset to focal individuals
  detResponse <- detResponse[dimnames(detResponse)[[1]] %in% dimnames(y.ar$y.ar)[[1]], ]
  
  
  
  ## ------     7.5. AUGMENT DETECTION HISTORIES -----
  
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
  y.status <- makeAugmentation( y = y.status,
                                aug.factor = aug.factor,
                                replace.value = 1)
  
  detResponse <- makeAugmentation( y = detResponse,
                                   aug.factor = aug.factor,
                                   replace.value = 0)
  ##-- Set first detection for augmented individuals to NA
  detResponse[rownames(detResponse) %in% "Augmented",1]  <- NA
  
  
  
  ## ------     7.6. TRANSFORM Y TO SPARSE MATRICES ------
  
  ##-- STRUCTURED
  y.sparse <- nimbleSCR::getSparseY(y.aliveStructured)
  
  ##-- OTHER
  y.sparseOth <- nimbleSCR::getSparseY(y.aliveOthers)
  
  
  
  ## ------ IV. MODEL SETTING ------- 
  
  ## ------   1. NIMBLE MODEL DEFINITION ------
  
  modelCode <- nimbleCode({
    
    ##------ SPATIAL PROCESS ------## 
    for(st in 1:2){
      dmean[st] ~ dunif(0,100)
      lambda[st] <- 1/dmean[st]
    }#st
    
    betaDens ~ dnorm(0.0,0.01)
    
    for(t in 1:n.years){
      habIntensity[1:n.habWindows,t] <- exp(betaDens * habDens[1:n.habWindows,t])
      sumHabIntensity[t] <- sum(habIntensity[1:n.habWindows,t])
      logHabIntensity[1:n.habWindows,t] <- log(habIntensity[1:n.habWindows,t])
      logSumHabIntensity[t] <- log(sumHabIntensity[t])
    }#t
    
    for(i in 1:n.individuals){
      sxy[i,1:2,1] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
        upperCoords = upperHabCoords[1:n.habWindows,1:2],
        logIntensities = logHabIntensity[1:n.habWindows,1],
        logSumIntensity = logSumHabIntensity[1],
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
      
      for(t in 2:n.years){
        sxy[i,1:2,t] ~ dbernppACmovement_exp(
          lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
          upperCoords = upperHabCoords[1:n.habWindows,1:2],
          s = sxy[i,1:2,t-1],
          lambda = lambda[state[i,t-1]+1],
          baseIntensities = habIntensity[1:n.habWindows,t],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max,
          numWindows = n.habWindows)
      }#t
    }#i
    
    
    ##----- DEMOGRAPHIC PROCESS -----##
    ##-- FIRST YEAR
    omeg1[1:3] ~ ddirch(alpha[1:3])  
    
    for(i in 1:n.individuals){
      z[i,1] ~ dcat(omeg1[1:3])
      isAlive[i,1] <- (z[i,1] == 2) + (z[i,1] == 3)
      state[i,1] <- (z[i,1] == 3)
    }#i
    
    ##-- FOLLOWING YEARS 
    for(t in 1:(n.years-1)){
      gamma[t] ~ dunif(0,1)
      psi[t] ~ dunif(0,1)
      
      for(st in 1:2){
        w[st,t] ~ dunif(0,1)
        h[st,t] ~ dunif(0,1)
        rw[st,t] ~ dunif(0,1)
        ones.dead.legal[st,t] ~ dbern(step(1 - (h[st,t] + w[st,t] + rw[st,t])))     
        phi[st,t] <- 1 - h[st,t] - w[st,t] - rw[st,t]
        wAll[st,t] <- w[st,t] + rw[st,t]
      }#st
      
      omega[1,1:6,t] <- c(1-gamma[t], gamma[t]           , 0              , 0     , 0      , 0     ) ## "UNBORN"
      omega[2,1:6,t] <- c(0         , phi[1,t]*(1-psi[t]), phi[1,t]*psi[t], h[1,t], rw[1,t], w[1,t]) ## "NON-PAIRS"
      omega[3,1:6,t] <- c(0         , 0                  , phi[2,t]       , h[2,t], rw[2,t], w[2,t]) ## "PAIRS"
      omega[4,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "NEWLY DEAD LEGAL HUNTING"
      omega[5,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "NEWLY DEAD OTHER SOURCES"
      omega[6,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "DEAD"
      
      for(i in 1:n.individuals){
        z[i,t+1] ~ dcat(omega[z[i,t],1:6,t])
        isAlive[i,t+1] <- (z[i,t+1] == 2) + (z[i,t+1] == 3)
        state[i,t+1] <- (z[i,t+1] == 3)
      }#i                                                                                                        
    }#t
    
    
    ##----- DETECTION PROCESS -----##
    pResponse ~ dunif(0,1)
    for(i in 1:n.individuals){
      detResponse[i,1] ~ dbern(pResponse)
    }
    
    for(t in 1:n.years){
      for(st in 1:2){
        sigma[st,t] ~ dunif(0,50)
      }#st
      
      ##-- Structured
      betaResponse[t] ~ dunif(-5,5)
      for(n in 1:n.covs){
        betaCovs[n,t] ~ dunif(-5,5)
      }#n
      for(c in 1:n.counties){
        p0[c,1,t] ~ dunif(0,1)
        p0[c,2,t] ~ dunif(0,1)
      }#c    
      
      ##-- Opportunistic
      betaResponseOth[t] ~ dunif(-5,5)
      for(n in 1:n.covsOth){
        betaCovsOth[n,t] ~ dunif(-5,5)
      }#n
      for(c in 1:n.countries){
        p0Oth[c,1,t] ~ dunif(0,1)
        p0Oth[c,2,t] ~ dunif(0,1)
      }#c
      
      for(i in 1:n.individuals){
        
        ##-- Structured
        y.alive[i,1:maxDetNums,t] ~ dbinomLocal_normalWolf(
          detNums = detNums[i,t],
          detIndices = detIndices[i,1:maxDetNums,t],
          size = size[1:n.detectors],
          p0 = p0[1:n.counties,1:2,t],
          sigma = sigma[state[i,t]+1,t],
          s = sxy[i,1:2,t],
          trapCoords = detector.xy[1:n.detectors,1:2],
          localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
          localTrapsNum = localDetNum[1:n.habWindows],
          resizeFactor = resizeFactor,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          indicator = isAlive[i,t],
          z = z[i,t]-1,
          trapCovsIntercept = detCounties[1:n.detectors],
          indCov = detResponse[i,t],
          indBeta = betaResponse[t],
          trapCovs = detCovs[1:n.detectors,1:n.covs,t],
          trapBetas = betaCovs[1:n.covs,t],
          lengthYCombined = 1)
        
        ##--- Opportunistic
        y.aliveOth[i,1:maxDetNumsOth,t] ~ dbinomLocal_normalWolf(
          detNums = detNumsOth[i,t],
          detIndices = detIndicesOth[i,1:maxDetNumsOth,t],
          size = size[1:n.detectors],
          p0 = p0Oth[1:n.countries,1:2,t],
          sigma = sigma[state[i,t]+1,t],
          s = sxy[i,1:2,t],
          trapCoords = detector.xy[1:n.detectors,1:2],
          localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
          localTrapsNum = localDetNum[1:n.habWindows],
          resizeFactor = resizeFactor,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          indicator = isAlive[i,t],
          z = z[i,t]-1,
          trapCovsIntercept = detCountries[1:n.detectors],
          indCov = detResponse[i,t],
          indBeta = betaResponseOth[t],
          trapCovs = detCovsOth[1:n.detectors,1:n.covsOth,t],
          trapBetas = betaCovsOth[1:n.covsOth,t],
          lengthYCombined = 1)
        
        ##-- Dead recoveries legal
        x.deadculled[i,t] ~ dbern(z[i,t] == 4)
        
        ##-- Dead recoveries others
        x.deadOther[i,t] ~ dbern(z[i,t] == 5)
      }#i
    }#t
    
    
    ##---------- DERIVED PARAMETERS ----------##
    for(t in 1:n.years){
      N[t] <- sum(isAlive[1:n.individuals,t])
    }#t
  })
  
  
  
  ## ------   2. NIMBLE CONSTANTS ------
  
  nimConstants <- list( 
    n.individuals = dim(y.sparse$y)[1],
    n.habWindows = nrow(habitat$scaledLowerCoords),
    n.detectors = nrow(detectors$scaledCoords),
    n.years = dim(y.sparse$y)[3], 
    n.covs = dim(detCovs)[2],
    n.covsOth = dim(detCovsOth)[2],
    n.countries = max(detCountries),
    n.counties = max(detCounties),
    resizeFactor = detectors$localObjects$resizeFactor,
    y.max = dim(detectors$localObjects$habitatGrid)[1],
    x.max = dim(detectors$localObjects$habitatGrid)[2],
    numLocalIndicesMax = detectors$localObjects$numLocalIndicesMax,
    maxDetNums = y.sparse$maxDetNums,
    maxDetNumsOth = y.sparseOth$maxDetNums)
  
  
  
  ## ------   3. NIMBLE DATA ------
  
  ## ------     3.1. GENERATE KNOWN z ------
  
  ##-- Create known values for z (2 = alive; 3 = recovered dead)
  z.data <- apply(y.alive, c(1,3), function(x) ifelse(any(x>0),2,0))
  z.dead <- apply(y.dead, c(1,2), function(x) ifelse(any(x>0),3,0))
  z.data <- ifelse( z.dead + z.data == 0, NA, z.dead + z.data)
  
  ##-- Fill in the gaps
  z.data <- t(apply(z.data, 1, function(zz){
    if(any(!is.na(zz))){
      ##-- Identify range of detections
      range.det <- range(which(!is.na(zz)))
      ##-- If not dead recovered, set state to 2 from first to last detection
      if(sum(zz == 3, na.rm = T)<1){
        zz[range.det[1]:range.det[2]] <- 2
      }
      ##-- If recovered
      if(sum(zz == 3, na.rm = T)>0){
        ##-- If another detection alive also available
        ##-- Set state to 2 from first detection to occasion before recovered
        if(sum(zz, na.rm = T)>3){
          zz[range.det[1]:(range.det[2]-1)] <- 2
        }
        
        t.recovered <- which(zz==3)
        ##-- If recovered before last occasion
        ##-- Set state to 4 from occ. after recovery to last occ
        if(t.recovered < length(zz)){
          zz[(t.recovered+1):length(zz)] <- 6
        }
        ##-- In any case set state to 2 for occ just before recovery
        if(t.recovered <= length(zz)){
          zz[(t.recovered-1)] <- 2
        }
      }
    }
    return(zz)
  }))
  
  ##-- Identify individuals dead to legal culling and to other reasons
  legal.mx <- other.mx <- matrix( 0,
                                  nrow(y.dead), ncol(y.dead),
                                  dimnames = dimnames(y.dead))
  legal.mx[data.dead$Id[data.dead$Legal], ] <- 1
  other.mx[data.dead$Id[!data.dead$Legal], ] <- 1
  
  ##-- Id culled get state 4, id not culled get 5. 
  z.data <- ifelse(z.data == 3 & other.mx %in% c(1), 5, z.data) 
  z.data <- ifelse(z.data == 3 & legal.mx %in% c(1), 4, z.data)
  
  ##-- Id alive and in pairs get state 3
  z.data <- ifelse(z.data == 2 & y.status %in% c(2), 3, z.data)
  
  
  
  ## ------     3.2. GENERATE x.dead ------
  
  x.deadculled <- x.deadOther <- matrix( 0,
                                         nrow(z.data), ncol(z.data),
                                         dimnames = dimnames(z.data))
  x.deadculled[z.data == 4] <- 1
  x.deadOther[z.data == 5] <- 1
  
  
  
  ## ------     3.3. LIST DATA ------
  
  nimData <- list( 
    z = z.data,   
    #sxy = sxy.data,
    
    y.alive = y.sparse$y,
    detIndices = y.sparse$detIndices,
    detNums = y.sparse$detNums,
    
    y.aliveOth = y.sparseOth$y, 
    detIndicesOth = y.sparseOth$detIndices,
    detNumsOth = y.sparseOth$detNums,
    
    x.deadculled = x.deadculled,
    x.deadOther = x.deadOther,
    
    habitatGrid = detectors$localObjects$habitatGrid,
    habDens = habDens,
    lowerHabCoords = as.matrix(habitat$scaledLowerCoords), 
    upperHabCoords = as.matrix(habitat$scaledUpperCoords), 
    detCounties = detCounties,
    detCountries = detCountries,
    detCovs = detCovs,
    detCovsOth = detCovsOth,
    detResponse = detResponse,
    localDetIndices = detectors$localObjects$localIndices,
    localDetNum = detectors$localObjects$numLocalIndices,
    size = detectors$detectors.df$size,
    detector.xy = as.matrix(detectors$scaledCoords),
    alpha = rep(1,3),
    ones.dead.legal = array(1,c(2,dim(y.alive)[3]-1)))
  
  
  
  ## ------   4. NIMBLE INITS ------
  
  ## ------     4.1. GENERATE INITIAL z ------
  
  ##-- Create initial values for z  
  z.init <- t(apply(z.data, 1, function(zz){
    out <- zz
    out[] <- 1
    
    if(any(!is.na(zz))){
      ##-- If not legally culled, set z to 1 before first detection and to 6 after last detection
      if(sum(zz == 4, na.rm = T) < 1){
        range.det <- range(which(!is.na(zz)))
        if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
        if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 6
      }
      ##-- If alive in a pair, set state to 2 the year before
      if(sum(zz == 3, na.rm = T) > 0){
        reco.3 <- min(which(zz == 3))
        if(reco.3>1){zz[(reco.3-1)] <- 2}
      }
      ##-- If alive not in pair, set state to 1 before first detection
      if(sum(zz == 2, na.rm = T) > 0){
        reco.alive <- min(which(zz == 2))
        if(reco.alive>1){zz[1:(reco.alive-1)] <- 1}
      }
      out[] <- zz
      
      ##-- if still some NAs initialize it with 2
      if(sum(is.na(out) > 0)){
        out[is.na(out)] <- 2
      }
    }
    return(out)
  }))
  z.init[!is.na(z.data)] <- NA
  
  
  
  ## ------     4.2. GENERATE INITIAL sxy ------ 
  
  ##-- Project death to the next year and combine all detections
  AllDetections <- data.dead %>% 
    select(Id, Year) %>%
    mutate(Year = Year + 1) %>%
    filter(Year != max(Year)) %>%
    rbind(., data.alive$data.sp[ ,c("Id", "Year")]) %>%
    mutate( "x" = st_coordinates(.)[,1],
            "y" = st_coordinates(.)[,2]) %>%
    st_drop_geometry() 
  
  ##-- Rescale all detections
  AllDetections <- scaleCoordsToHabitatGrid(
    coordsData = AllDetections,
    coordsHabitatGridCenter = habitat$habitat.xy,
    scaleToGrid = T)$coordsDataScaled
  
  ##-- Identify augmented individuals
  idAugmented <- which(rownames(z.data) %in%"Augmented")
  
  ##-- Generate initial values for sxy
  sxy.init <- getSInits( AllDetections = AllDetections,
                         Id.vector = y.ar$Id.vector,
                         idAugmented = idAugmented,
                         lowerCoords = as.matrix(habitat$scaledLowerCoords),
                         upperCoords = as.matrix(habitat$scaledUpperCoords),
                         habitatGrid = detectors$localObjects$habitatGrid,
                         intensity = NULL,
                         sd = 4,
                         movementMethod = "dbernppACmovement_normal")
  
  ##-- SXY data 
  sxy.data <- sxy.init
  sxy.data[ ] <- NA
  for(i in 1:length(y.ar$Id.vector)){
    for(t in 1:(dim(sxy.data)[3]-1)){
      if(sum(z.data[i,t+1] %in% c(4,5,6)) > 0){
        # print(i)
        # print(t)
        sxy.data[i, ,t+1] <- sxy.init[i, ,t+1] 
        sxy.init[i, ,t+1] <- NA
      }#if
    }#t
  }#i
  
  sxy.init <- round(sxy.init,5) #-- an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
  sxy.data <- round(sxy.data,5)
  
  nimData$sxy <- sxy.data
  
  
  
  ## ------     4.3. LATENT VARIABLE DET RESPONSE ------ 
  
  detResponse.init <- detResponse
  detResponse.init[is.na(detResponse)] <- rbinom(sum(is.na(detResponse)), 1,0.5)
  detResponse.init[!is.na(detResponse)] <- NA
  
  
  
  ## ------   5. NIMBLE PARAMETERS ------ 
  
  nimParams <- c("N","lambda","dmean","betaDens",
                 "omeg1","gamma","psi","phi","h","w","wAll","rw",
                 "pResponse","sigma",
                 "p0","betaResponse","trapBetas",
                 "p0Oth", "betaResponseOth","trapBetasOth")
  
  nimParams2 <- c("z", "sxy")
  
  
  
  ## ------   6. SAVE INPUTS ----- 
  
  for(c in 1:4){
    
    nimInits <- list( 
      "sxy" = sxy.init,
      "dmean" = runif(2,2,4),
      "z" = z.init,
      "omeg1" = c(0.5,0.25,0.25),
      "psi" = runif(dim(y.alive)[3]-1,0.1,0.7),
      "gamma" = runif(dim(y.alive)[3]-1,0,1),
      "h" =  array( runif((dim(y.alive)[3]-1)*2,0.2,0.4),
                    c(2,dim(y.alive)[3]-1)),
      "rw" =  array( runif((dim(y.alive)[3]-1)*2,0.05,0.10),
                     c(2,dim(y.alive)[3]-1)),
      "w" = array( runif((dim(y.alive)[3]-1)*2,0.2,0.4),
                   c(2,dim(y.alive)[3]-1)),
      "pResponse"  = runif(1,0,1),
      "detResponse" = detResponse.init,
      "sigma" = array(runif(2,4,8),
                      c(2,dim(y.alive)[3])),
      "p0" = array( runif(12,0,0.2), 
                    c(nimConstants$n.counties,2,dim(y.alive)[3])),
      "p0Oth" = array( runif(12,0,0.2),
                       c(nimConstants$n.countries,2,dim(y.alive)[3])),
      "betaResponse" = runif(dim(y.alive)[3],-1,1),
      "betaResponseOth" = runif(dim(y.alive)[3],-1,1),
      "betaDens" = runif(1,-1,1),
      "betaCovs" = array( runif(nimConstants$n.covs,-1,1),
                          c(nimConstants$n.covs,dim(y.alive)[3])),
      "betaCovsOth" = array( runif(nimConstants$n.covsOth,-1,1),
                             c(nimConstants$n.covsOth,dim(y.alive)[3])))
    
    ### TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK
    ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
    for(i in 1:nimConstants$n.individuals){
      for(t in 1:nimConstants$n.years){
        if(!is.na(nimInits$sxy[i,1,t])){
          SXY <- nimInits$sxy[i,,t]
        } else {
          SXY <- nimData$sxy[i,,t]
        }
        sID <- nimData$habitatGrid[trunc(SXY[2]/nimConstants$resizeFactor)+1,
                                   trunc(SXY[1]/nimConstants$resizeFactor)+1]
        index <- nimData$localDetIndices[sID, 1:nimData$localDetNum[sID]]
        ## GET NECESSARY INFO
        n.detectors <- length(index)
        YDET <- nimData$detIndices[i,, t]
        ## check if a detection is out of the "detection window"
        if(nimData$detNums[i,t] > 0){
          for(j in 1:nimData$detNums[i,t]){
            if(sum(YDET[j]==index)==0){
              print(paste("id",i,"t",t,"j",j))
            }
          }
        }
      }#t
    }#i
    
    save( modelCode,
          nimData,
          nimConstants,
          nimParams,
          nimParams2,
          nimInits,
          file = file.path( working.dir, "nimbleInFiles", thisSex,
                            paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
  }#c
}#sex



##------------------------------------------------------------------------------
## ------   8. NIMBLE RUN ------ 

load(file.path( working.dir, "nimbleInFiles", thisSex,
                paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
ptm <- proc.time()
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()
cmodel <- compileNimble(model)
cmodel$calculate() r
which(is.infinite(model$logProb_z),arr.ind = T)


probs <- which(is.na(model$logProb_y),arr.ind=T)

model$y



##------------------------------------------------------------------------------

## ------ IV. MAKE THE SINGLE SEASON SCR MODEL ------- 

## ------   1. NIMBLE MODEL DEFINITION ------ 

modelCode1 <- nimbleCode({
  
  ##------ SPATIAL PROCESS ------##  
  beta.dens  ~ dnorm(0.0,0.01)
  habIntensity[1:n.habWindows] <- exp(beta.dens * habDens[1:n.habWindows])
  sumHabIntensity <- sum(habIntensity[1:n.habWindows])
  logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:n.habWindows, 1:2],
      upperCoords = upperHabCoords[1:n.habWindows, 1:2],
      logIntensities = logHabIntensity[1:n.habWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  
  ##----- DEMOGRAPHIC PROCESS -----## 
  psi ~ dunif(0,1)
  probAdult ~ dunif(0,1)
  pResponse ~ dunif(0,1)
  for(i in 1:n.individuals){
    z[i] ~ dbern(psi)
    state[i] ~ dbern(probAdult)
    idResponse[i] ~ dbern(pResponse)
  }
  
  
  
  ##----- DETECTION PROCESS -----## 
  for(st in 1:2){
    sigma[st] ~ dunif(0,20)
  }
  
  betaResponse ~ dunif(-5,5)
  betaResponseOth ~ dunif(-5,5)
  
  ##-- STRUCTURED
  for(c in 1:n.counties){
    p0[c,1] ~ dunif(0,1)
    p0[c,2] ~ dunif(0,1)
  }
  for(c in 1:nTrapCovs){
    trapBetas[c] ~ dunif(-5,5)
  }
  
  ##-- OTHERS
  for(c in 1:n.countries){
    p0Oth[c,1] ~ dunif(0,1)
    p0Oth[c,2] ~ dunif(0,1)
  }
  for(c in 1:nTrapCovsOth){
    trapBetasOth[c] ~ dunif(-5,5)
  }
  
  for(i in 1:n.individuals){
    
    y.alive[i,1:nMaxDetectors] ~ dbinomLocal_normalWolf(
      detNums = nbDetections[i],
      detIndices = yDets[i,1:nMaxDetectors],
      size = trials[1:n.detectors],
      p0 = p0[1:n.counties,1:2],
      sigma = sigma[state[i]+1],
      s = sxy[i,1:2],
      trapCoords = detector.xy[1:n.detectors,1:2],
      localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
      localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      z = state[i]+1,
      trapCovsIntercept =  detCounties[1:n.detectors],
      indCov = idResponse[i],
      indBeta = betaResponse,
      trapCovs =  trapCovs[1:n.detectors,1:nTrapCovs],
      trapBetas = trapBetas[1:nTrapCovs],
      lengthYCombined = 1)
    
    y.aliveOth[i,1:nMaxDetectorsOth] ~ dbinomLocal_normalWolf(
      detNums = nbDetectionsOth[i],
      detIndices = yDetsOth[i,1:nMaxDetectorsOth],
      size = trials[1:n.detectors],
      p0 = p0Oth[1:n.countries,1:2],
      sigma = sigma[state[i]+1],
      s = sxy[i,1:2],
      trapCoords = detector.xy[1:n.detectors,1:2],
      localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
      localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
      resizeFactor = ResizeFactor,
      habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
      indicator = z[i],
      z = state[i]+1,
      trapCovsIntercept =  detCountries[1:n.detectors],
      indCov = idResponse[i],
      indBeta = betaResponseOth,
      trapCovs =  trapCovsOth[1:n.detectors,1:nTrapCovsOth],
      trapBetas = trapBetasOth[1:nTrapCovsOth],
      lengthYCombined = 1)
  }#i
  
  
  ##---------- DERIVED PARAMETERS ----------##
  N <- sum(z[1:n.individuals])
  
})



## ------   2. SAVE SCR NIMBLE INPUTS ------ 

for(ch in 1:4){
  
  for(t in 1:n.years){   
    
    load(file.path( WD,
                    modelName,
                    paste0(modelName, "_INPUTChain", ch, ".RData")))
    
    
    ## ------   2.1. GET WHICH INDIVIDUAL IS DETECTED ------ 
    
    detectedStruc <- apply(nimData$nbDetections,2,function(x) x>0)
    detectedOth <- apply(nimData$nbDetectionsOth,2,function(x) x>0)
    detected <- detectedOth + detectedStruc
    detected <- detected>0
    
    # GET SUM OF INDIVIDUALS DETECTED AND DECIDE HOW MUCH YOU WISH TO AUGMENT. hERE I CHOSE 2
    sumDets <- sum(detected[,t])*3## decide the augmentation factor 
    
    
    
    ## ------   2.2. DATA ------ 
    
    ## y.alive
    nimData$y.alive <- nimData$y.alive[detected[,t],,t]  
    nimData$y.alive <- rbind(nimData$y.alive, matrix(0,nrow =sumDets, nimConstants$nMaxDetectors))
    nimData$y.aliveOth <- nimData$y.aliveOth[detected[,t],,t]  
    nimData$y.aliveOth <- rbind(nimData$y.aliveOth, matrix(0,nrow =sumDets, nimConstants$nMaxDetectorsOth))
    
    ## z, state
    nimData$z <- state <-  nimData$z[detected[,t],t]  
    state[] <- 0
    state[nimData$z %in% c(3)] <- 1
    nimData$z[nimData$z %in% c(2,3)] <- 1 # ALIVE IDS BECOMES 1
    nimData$z <- c(nimData$z, rep(NA,sumDets))
    nimData$state <- c(state, rep(NA,sumDets))
    
    ## sxy 
    nimData$sxy <- nimData$sxy[detected[,t],,t]  
    nimData$sxy <- rbind(nimData$sxy, matrix(NA, nrow =sumDets, ncol=2))
    
    ## nbDetections 
    nimData$nbDetections <- nimData$nbDetections[detected[,t],t]
    nimData$nbDetections  <- c(nimData$nbDetections, rep(0,sumDets))
    nimData$nbDetectionsOth <- nimData$nbDetectionsOth[detected[,t],t]
    nimData$nbDetectionsOth  <- c(nimData$nbDetectionsOth, rep(0,sumDets))
    
    ## yDets 
    nimData$yDets <- nimData$yDets[detected[,t],,t]
    nimData$yDets <- rbind(nimData$yDets, matrix(0,nrow =sumDets, nimConstants$nMaxDetectors))
    nimData$yDetsOth <- nimData$yDetsOth[detected[,t],,t]
    nimData$yDetsOth <- rbind(nimData$yDetsOth, matrix(0,nrow =sumDets, nimConstants$nMaxDetectorsOth))
    
    ## detResponse 
    nimData$idResponse <- nimData$idResponse[detected[,t],t]
    nimData$idResponse  <- c(nimData$idResponse, rep(NA,sumDets))## HERE IT IS ASSUMING IT IS A LATENT INDIVIDUAL COVARIATE
    
    ## detCovs
    nimData$trapCovs <- nimData$trapCovs[,,t]
    nimData$trapCovsOth <- nimData$trapCovsOth[,,t]
    
    ## density
    nimData$habDens <- nimData$habDens[,t]
    
    
    
    ## ------   2.3. INITS ------
    
    ## z
    nimInits$z <- nimInits$z[detected[,t],t]  
    nimInits$z <- c(nimInits$z, rbinom(sumDets,1,0.5))
    
    ## state
    nimInits$state <- nimData$state
    nimInits$state[is.na(nimData$state)] <- rbinom(sum(is.na(nimInits$state)),1,0.5)
    nimInits$state[!is.na(nimData$state)] <- NA
    
    ## detResponse
    ## HERE IT IS TREATED AS A LATENT COVARIATE
    nimInits$idResponse  <- c(rep(NA,sum(detected[,t])), rbinom(sumDets,1,0.5))
    
    ## sxy 
    nimInits$sxy <- nimInits$sxy[detected[,t],,t]  
    # GIVE ACS FROM DETECTED INDIVIDUALS TO AUGMENTED IDS. 
    nimInits$sxy <- rbind(nimInits$sxy, nimInits$sxy[sample(nimInits$sxy, sumDets, replace = T),])
    
    ## p0
    nimInits$p0 <-  nimInits$p0[,,t]
    nimInits$p0Oth <-  nimInits$p0Oth[,,t]
    
    ## psi
    nimInits$psi <-  runif(1,0.4,0.6)
    
    ## sigma
    nimInits$sigma <-  runif(2,4,8)
    
    ## trapBetas
    nimInits$trapBetas <- nimInits$trapBetas[,t]
    nimInits$trapBetasOth <- nimInits$trapBetasOth[,t]
    
    ## trapBetas
    nimInits$betaResponse <- nimInits$betaResponse[t]
    nimInits$betaResponseOth <- nimInits$betaResponseOth[t]
    
    ## beta.dens
    nimInits$beta.dens <- nimInits$beta.dens
    
    ## parameter for the latent detResponse covariate
    nimInits$probDetBefore <-  runif(1,0.4,0.6)
    nimInits$probAdult  <-  runif(1,0.4,0.6)
    
    
    ## ------   2.4. CONSTANTS ------
    
    ## get the new number of individuals 
    nimConstants$n.individuals <- nrow(nimInits$sxy)
    
    
    
    ## ------   2.5. PARAMS ------
    
    nimParams <- c("N", "psi","probAdult", "pResponse","p0Oth","trapBetasOth","betaResponseOth",
                   "p0", "sigma", "beta.dens", "trapBetas","betaResponse")
    
    nimParams2 <- c("state", "z", "sxy")
    
    
    ## ------   2.6. SAVE INPUT ------
    
    modelCode <- modelCode1
    
    save(nimData,
         nimConstants,
         nimParams,
         nimParams2,
         modelCode,
         nimInits,
         file = file.path(
           WD,
           modelName,
           paste0("Snap", modelName, years[t], "_", ch, ".RData")))
  }#t
}#ch



## ------   3. NIMBLE RUN ------ 

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()
model$initializeInfo()
model$nodeFunctionGeneratorNames


cmodel <- compileNimble(model)
cmodel$calculate() 

MCMCconf <- configureMCMC( model = model,
                           monitors = nimParams,
                           control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                           thin = 1) 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC,project = model,resetFunctions = TRUE)
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 0,
                                                      niter = 100,
                                                      nchains = 1,
                                                      inits = nimInits,
                                                      samplesAsCodaMCMC = TRUE))
TotalRuntime <- proc.time()-ptm

save(myNimbleOutput,
     MCMCRuntime,
     TotalRuntime,
     file = file.path(
       WD, 
       modelName,
       paste0(modelName,"_OUTPUT.RData")))



##------------------------------------------------------------------------------
####---------------------------------- 
####---- processRovQuantOutput() -----
#processRovquantOutput_wolf <- function(
    ##-- paths
# data.dir = getwd(),
# working.dir = NULL,
##-- MCMC
nburnin = 0
niter = 100
##-- Density 
extraction.res = 5000
##-- Miscellanious
overwrite = FALSE
#){
####---------------------------------- 


## ------   1. LOAD & PROCESS OPSCR OUTPUTS ------ 

## COMBINE DIFFERENT CHAINS (IF ANY) INTO ONE 
outDirectories <- list.files(file.path(WD, modelName))[grep("NimbleOut", list.files(file.path(WD, modelName)))]
path.list <- file.path(WD, modelName, outDirectories)[c(1:4)]

nthin <- 1
# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)/2
}))
minBites <- min(numBites)
minBites

#Niter to Remove (burn-in) #CM
NSkipBites <- 40
nimOutput <- RUNTIME <- list()
gc()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3] 
    params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
    parmIndex <- which(! params.simple %in% c("sxy","z"))
    nthins <- seq(1,dim(this.sample)[1], by=nthin)
    out[[x]] <- this.sample[nthins,]#[ ,parmIndex] 
    if(sum(is.na(out[[x]]))>0){
      out[[x]] <- out[[x]][-unique(which(is.na(out[[x]]),arr.ind = T)[,1]),]
    }
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  nimOutput[[p]] <- as.mcmc(out.mx)
}#p

lapply(RUNTIME, function(x) x/3600)#[CM]
unlist(lapply(RUNTIME, function(x) x/3600))#[CM]
TIME <- lapply(RUNTIME, function(x) x/3600)

max <- unlist(lapply(RUNTIME, function(x) x/3600))
at=c(1:length(TIME[[1]]))
plot(TIME[[1]]~at, pch=16, col=adjustcolor("red", alpha.f = 0.5), xlim=c(0,length(TIME[[1]])+2),
     ylim=c(0,10), ylab="time hours", xlab="bite number")
for(i in 2:length(TIME)){
  points(TIME[[i]] ~ at, )
}
nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))

t=8
table(myResults$sims.list$z[1,,t]==3 & myResults$sims.list$z[1,,t-1]==2)

table(myResults$sims.list$z[1,,8])
table(myResults$sims.list$z[1,,7])

gc()



## ------   2. PLOT PARAMETERS ESTIMATES ------ 

{#doall}
  pdf(file=file.path(WD, modelName, paste(modelName,"_ESTIMATES.pdf",sep="")))
  
  ## ------  2.1.N ------  
  
  plot(10, xlim = c(0, n.years+1), ylim = c(80,350), type ="n", xaxt="n", xlab = "Years", ylab = "N")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResults$sims.list$N[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("N",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  # ## ------  2.2.rho ------ 
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("gamma",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------  [RB] h[t] ------ 
  
  plot(10, xlim = c(0, n.years+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "h")
  axis(1, c(1:n.years),labels = years)
  myCol <- c("firebrick3","navyblue")
  myDev <- c(-0.2,+0.2)
  for(t in 1:(n.years-1)){
    for(s in 1:2){
      plot.violins(list(myResults$sims.list$h[,s,t]),
                   x = t+ myDev[s],
                   at = t+myDev[s],
                   violin.width = 0.3,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }#s
  }#t
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("h\\[",dimnames(nimOutput[[1]])[[2]])]
  params <- params[1:2]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------  [RB] wall[t] ------  
  
  plot(10, xlim = c(0, n.years+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "wAll")
  axis(1, c(1:n.years),labels = years)
  myDev <- c(-0.2,+0.2)
  for(t in 1:(n.years-1)){
    for(s in 1:2){
      plot.violins(list(myResults$sims.list$wAll[,s,t]),
                   x = t+ myDev[s],
                   at = t+ myDev[s],
                   violin.width = 0.3,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }#s
  }#t
  
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("wAll",dimnames(nimOutput[[1]])[[2]])]
  
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------  [RB] rw[t] ------  
  
  plot(10, xlim = c(0, n.years+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "rw")
  axis(1, c(1:n.years),labels = years)
  myDev <- c(-0.2,+0.2)
  for(t in 1:(n.years-1)){
    for(s in 1:2){
      plot.violins(list(myResults$sims.list$rw[,s,t]),
                   x = t+ myDev[s],
                   at = t+ myDev[s],
                   violin.width = 0.3,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }#s
  }#t
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("rw",dimnames(nimOutput[[1]])[[2]])]
  
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------   2.3.phi ------ 
  
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "phi")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  for(t in 1:(n.years-1)){
    for(s in 1:2){
      plot.violins(list(myResults$sims.list$phi[ ,s, t]),
                   x = t + myDev[s],
                   at = t + myDev[s],
                   violin.width = 0.2,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }
  }
  
  legend("bottomright", fill = c("firebrick3","navyblue"), legend = c("not pair", "pair"))
  params <- dimnames(nimOutput[[1]])[[2]][grep("phi",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------   2.4.p0 ------   
  
  par(mfrow=c(1,2))
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  par(mfrow=c(1,2))
  for(c in 1:6){
    plot(-10, xlim = c(0,n.years+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(n.years) , labels = years[1:(n.years)])
    for(s in 1:2){
      for(t in 1:n.years){
        plot.violins(list(myResults$sims.list$p0[ , c, s, t]),
                     x = t + myDev[s],
                     at = t + myDev[s],
                     violin.width = 0.2,
                     col = myCol[s],
                     add = T,
                     alpha = 0.2,
                     border.col = myCol[s])
      }
    }
    plot(COUNTIESplot$geometry)
    plot(COUNTIES[COUNTIES$id==detCounties.original[c], ],add=T, col="red")
  }
  
  # text(COUNTIESplot
  #      ,labels=COUNTIESplot$id, col="red")
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("p0\\[",dimnames(nimOutput[[1]])[[2]])[-1]]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  ## ------   2.4.p0Oth ------   
  par(mfrow=c(1,2))
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  par(mfrow=c(1,2))
  main= c("Sweden", "Norway")
  for(c in 1:2){
    plot(-10, xlim = c(0,n.years+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=main[c])#paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(n.years) , labels = years[1:(n.years)])
    for(s in 1:2){
      for(t in 1:n.years){
        plot.violins(list(myResults$sims.list$p0Oth[ , c, s, t]),
                     x = t + myDev[s],
                     at = t + myDev[s],
                     violin.width = 0.2,
                     col = myCol[s],
                     add = T,
                     alpha = 0.2,
                     border.col = myCol[s])
      }
    }
    # plot(COUNTIESplot)
    # plot(COUNTIES[COUNTIES$id==detCounties.original[c], ],add=T, col="red")
  }
  
  # text(COUNTIESplot
  # ,labels=COUNTIESplot$id, col="red")
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("p0Oth",dimnames(nimOutput[[1]])[[2]])[-1]]
  
  
  
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  
  
  ## ------   2.5.psi ------   
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "psi")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$psi[ , t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  params <- dimnames(nimOutput[[1]])[[2]][grep("psi",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  ##sigma ##
  par(mfrow=c(1,1))
  offsetVal <- c(-0.2,0.2)
  plot(-10, xlim = c(0,n.years), ylim=c(0,10000), type ="n", xaxt="n", xlab = "Years", ylab = "sigma")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  abline(v=seq(1.5,n.years-0.5,by=1),lty=2)
  for(s in 1:2){
    for(t in 1:(n.years-1)){
      plot.violins(list(myResults$sims.list$sigma[ ,s,t]*habitat$resolution),
                   x = t ,
                   at = t +offsetVal[s] ,
                   violin.width = 0.2,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }
  }
  
  ###
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,5), ylim=c(0,150000), type ="n", xaxt="n", xlab = "", ylab = "km")
  axis(1, at = c(1.5,3.5) , labels = c("sigma","tau"))
  #sigma
  # for(s in 1:2){
  #     plot.violins(list(myResults$sims.list$sigma[ , s]*habitat$resolution),
  #                x = t ,
  #                at = s ,
  #                violin.width = 0.2,
  #                col = myCol[s],
  #                add = T,
  #                alpha = 0.2,
  #                border.col = myCol[s])
  # }
  #tau
  for(s in 1:2){
    plot.violins(list(myResults$sims.list$lambda[ , s]*habitat$resolution),
                 x = t ,
                 at = s+2 ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("sigma",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  ##beta.dens
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta.dens")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$beta.dens),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  params <- dimnames(nimOutput[[1]])[[2]][grep("beta.dens",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  ##trapBetas[1]
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-2,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta Tracks")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$trapBetas[ ,1,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  
  ##trapBetas[2]
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta Snow")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$trapBetas[ ,2,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("trapBetas",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  ##betaOth road 
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta other road")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$trapBetasOth[ ,1,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  # betaOth snow 
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta other snow")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$trapBetasOth[ ,2,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  # betaOth other 
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "beta other location other samples ")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$trapBetasOth[ ,3,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("trapBetasOth",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  ##betaResponse
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponse")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$betaResponse[ ,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  
  par(mfrow=c(1,1))
  plot(-10, xlim = c(0,n.years), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponseOth")
  axis(1, at = 1:(n.years-1) , labels = years[1:(n.years-1)])
  for(t in 1:(n.years-1)){
    plot.violins(list(myResults$sims.list$betaResponseOth[ ,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = myCol[s],
                 add = T,
                 alpha = 0.2,
                 border.col = myCol[s])
  }
  params <- dimnames(nimOutput[[1]])[[2]][grep("betaResponse",dimnames(nimOutput[[1]])[[2]])]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  
  ## ------  2.5.the rest ------ 
  params <- c("lambda[1]",
              "lambda[2]",
              "pResponse")
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  
  dev.off()
}#do all



## ------   3. LOAD & PROCESS SNAP OUTPUTS ------ 

# List the directories containing bite outputs
outDirectories <- list.files(file.path(WD, modelName))[grep("NimbleOut",
                                                            list.files(file.path(WD,
                                                                                 modelName)))]

outDirectories <- outDirectories[grep("Snap",outDirectories)]

#path.list <- file.path(WD, modelName, outDirectories)

t=1
nimOutputList <- myResultsList <-  list()

for(t in 1:n.years){
  path.list <- 0
  for(ch in 1:4){
    path.list[ch] <- file.path(WD, modelName, paste("NimbleOutFORSnap" ,
                                                    modelName,years[t],"_", ch, ".RData", sep = "")) 
  }
  
  # Retrieve the minimum number of bites per chain
  numBites <- unlist(lapply(path.list, function(x){
    files <- list.files(x)
    files <- files[grep(".RData", files)]
    length(files)/2
  }))
  
  path.list <- path.list[!numBites%in%0]
  numBites<- numBites[!numBites%in%0]
  minBites <- min(numBites)
  #Niter to Remove (burn-in) #CM
  NSkipBites <- 55
  nthin <- 1
  nimOutput <- RUNTIME <- list()
  #gc()
  for(p in 1:length(path.list)){
    print(path.list[p])
    outfiles <- list.files(path.list[p])
    out <- runtime <- list()#[CM]
    for(x in NSkipBites:minBites){
      print(x)
      load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
      runtime[[x]] <- RunTime[3] 
      params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
      parmIndex <- which(! params.simple %in% c("sxy","z"))
      nthins <- seq(1,dim(this.sample)[1], by=nthin)
      out[[x]] <- this.sample[nthins,parmIndex]#[ ,parmIndex] 
      
      # out[[x]] <- this.sample[nthins,]#[ ,parmIndex] 
    }#x
    RUNTIME[[p]] <- unlist(runtime)#[CM]
    out.mx <- do.call(rbind, out)
    nimOutput[[p]] <- as.mcmc(out.mx)
  }#p
  
  
  lapply(RUNTIME, function(x) x/3600)
  unlist(lapply(RUNTIME, function(x) x/3600))
  
  
  for(i in 1: length(nimOutput)){
    nimOutput[[i]][,"sigma[2]"] <-  nimOutput[[i]][,"sigma[2]"]*habitat$resolution
    nimOutput[[i]][,"sigma[1]"] <-  nimOutput[[i]][,"sigma[1]"]*habitat$resolution
  }
  
  nimOutputList[[t]] <- nimOutput <- as.mcmc.list(nimOutput)
  
  myResultsList[[t]] <- myResults <- ProcessCodaOutput(nimOutput, params.omit = c("sxy","z"))
  
} 
gc()


### combine all years together 
myResultsListALL <-  myResultsList[[1]] 

myResultsListALL$sims.list$N <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$N )) 
myResultsListALL$sims.list$beta.dens <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$beta.dens )) 
myResultsListALL$sims.list$betaResponse <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$betaResponse )) 

myResultsListALL$sims.list$pResponse <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$pResponse )) 
#myResultsListALL$sims.list$sigma <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$sigma )) 
myResultsListALL$sims.list$psi <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$psi )) 

myResultsListALL$sims.list$trapBetas <- array(NA,c(dim(myResultsListALL$sims.list$trapBetas),n.years) )
myResultsListALL$sims.list$trapBetasOth <- array(NA,c(dim(myResultsListALL$sims.list$trapBetasOth),n.years) )

myResultsListALL$sims.list$p0 <- array(NA,c(dim(myResultsListALL$sims.list$p0),n.years) )
myResultsListALL$sims.list$p0Oth <- array(NA,c(dim(myResultsListALL$sims.list$p0Oth),n.years) )

myResultsListALL$sims.list$sigma <- array(NA,c(dim(myResultsListALL$sims.list$sigma),n.years) )

for(t in 1:n.years){
  myResultsListALL$sims.list$trapBetas[,,t] <- myResultsList[[t]]$sims.list$trapBetas
  myResultsListALL$sims.list$trapBetasOth[,,t] <- myResultsList[[t]]$sims.list$trapBetasOth
  myResultsListALL$sims.list$sigma[,,t] <- myResultsList[[t]]$sims.list$sigma
  myResultsListALL$sims.list$p0[,,,t] <- myResultsList[[t]]$sims.list$p0
  myResultsListALL$sims.list$p0Oth[,,,t] <- myResultsList[[t]]$sims.list$p0Oth
}



## ------   4. PLOT PARAMETERS ESTIMATES ------ 

{#doall}
  pdf(file=file.path(WD, modelName, paste(modelName,"Snap_ESTIMATES.pdf",sep="")))
  
  ## ------  2.1.N ------  
  
  plot(10, xlim = c(0, n.years+1), ylim = c(80,350), type ="n", xaxt="n", xlab = "Years", ylab = "N")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$N[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "N")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.4.p0 ------   
  
  par(mfrow=c(1,2))
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  par(mfrow=c(1,2))
  for(c in 1:6){
    plot(-10, xlim = c(0,n.years+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(n.years) , labels = years[1:(n.years)])
    for(s in 1:2){
      for(t in 1:n.years){
        plot.violins(list(myResultsListALL$sims.list$p0[ , c, s, t]),
                     x = t + myDev[s],
                     at = t + myDev[s],
                     violin.width = 0.2,
                     col = myCol[s],
                     add = T,
                     alpha = 0.2,
                     border.col = myCol[s])
      }
    }
    plot(COUNTIESplot)
    plot(COUNTIES[COUNTIES$id==detCounties.original[c], ],add=T, col="red")
  }
  
  text(COUNTIESplot
       ,labels=COUNTIESplot$id, col="red")
  
  params <- dimnames(nimOutputList[[t]][[1]])[[2]][grep("p0\\[",dimnames(  nimOutputList[[t]][[1]])[[2]])[-1]]
  
  for(t in 1:n.years){
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutputList[[t]], params = params[i])
      title(years[t], line = 0.5)
    }
  }
  
  
  
  ## ------   2.4.p0Oth ------   
  
  par(mfrow=c(1,2))
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  par(mfrow=c(1,2))
  main= c("Sweden", "Norway")
  for(c in 1:2){
    plot(-10, xlim = c(0,n.years+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=main[c])#paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(n.years) , labels = years[1:(n.years)])
    for(s in 1:2){
      for(t in 1:n.years){
        plot.violins(list(myResultsListALL$sims.list$p0Oth[ , c, s, t]),
                     x = t + myDev[s],
                     at = t + myDev[s],
                     violin.width = 0.2,
                     col = myCol[s],
                     add = T,
                     alpha = 0.2,
                     border.col = myCol[s])
      }
    }
  }
  
  params <- dimnames(nimOutputList[[t]][[1]])[[2]][grep("p0Oth",dimnames(  nimOutputList[[t]][[1]])[[2]])[-1]]
  for(t in 1:n.years){
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutputList[[t]], params = params[i])
      title(years[t], line = 0.5)
    }
  }
  
  
  
  ## ------   2.5.sigma ------   
  
  plot(10, xlim = c(0, n.years+1), ylim = c(5000,12000), type ="n", xaxt="n", xlab = "Years", ylab = "sigma")
  axis(1, c(1:n.years),labels = years)
  offset <- c(-0.25,0.25)
  for(t in 1:n.years){
    for(s in 1:2){
      plot.violins(list(myResultsListALL$sims.list$sigma[,s,t]),
                   x = t,
                   at = t+offset[s],
                   violin.width = 0.3,
                   col = myCol[s],
                   add = T,
                   alpha = 0.2,
                   border.col = myCol[s])
    }
  }#t
  
  legend("topleft",legend=c("other","scent-marking"),col=myCol,pch=16)
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "sigma[1]")
    title(years[t],line = 0.5)
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "sigma[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.5.roads ------  
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-3,3), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Tracks")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$trapBetas[,1,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetas[1]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6.tracks  ------ 
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$trapBetas[,2,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetas[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------   
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta road (other)")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,1,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[1]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------ 
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow (other)")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,2,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------ 
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow (opportunistic)")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,3,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[3]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6. pResponse ------   
  
  plot(10, xlim = c(0, n.years+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "pResponse")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$pResponse[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "pResponse")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6. betaResponse ------   
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-1,2), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponse")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$betaResponse[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "betaResponse")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. beta.dens ------   
  
  plot(10, xlim = c(0, n.years+1), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta.dens")
  axis(1, c(1:n.years),labels = years)
  for(t in 1:n.years){
    plot.violins(list(myResultsListALL$sims.list$beta.dens[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:n.years){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "beta.dens")
    title(years[t],line = 0.5)
  }
  
  dev.off()
}#do all

##------------------------------------------------------------------------------