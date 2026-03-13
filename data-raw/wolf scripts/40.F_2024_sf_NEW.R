
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


## ------ SET REQUIRED WORKING DIRECTORIES ------

source("C:/My_documents/RovQuant/workingDirectories.R")             

##-- DATA DIRECTORY
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2024/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/2024/Test_F_NEW"


## ------ SOURCE THE REQUIRED FUNCTIONS ------

source(file.path(dir.git,"Temp/CM/functions/Nimble/dbinomLocal_normalWolf.R"))


## -----------------------------------------------------------------------------

## ------ 0. SET-UP PARAMETERS ------

makeDirectories( path = working.dir,
                 subFolders = c("female","male"),
                 show.dir = TRUE)

plot.check = TRUE

years = 2015:2024
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))
species = "Ulv"           
SPECIES <- "Gray wolf"
engSpecies <- "wolf"
norSpecies <- "Ulv"

sex = c("male","female") 

legal.dead <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")

##-- Set default values for the wolf model
aug.factor <- 0.8
sampling.months <- list(10:12,1:4)
habitat.res <- 20000
x.extent <- c(210000,740000)
y.extent <- c(6000000,7050000)
buffer.size <- 40000
detector.res <- 10000
subdetector.res <- 1000
max.det.dist <- 45000
resize.factor <- 1


##-- Renaming list
#if(is.null(rename.list)) {
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
  Sex = "Kön",
  Site_quality = "Stedkvalitet",
  SVAID = "SVAID",
  Uncertain_date = "Usikker dødsdato",
  Weight_slaughter = "Slaktevekt",
  Weight_total =  "Helvekt")
#}

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



## -----------------------------------------------------------------------------

## ------ I. LOAD & SELECT DATA ------

## ------   1. HABITAT DATA -----

##-- Load pre-defined habitat rasters and shapefiles
data(COUNTRIES, envir = environment()) 
data(COUNTIES, envir = environment()) 
data(habitatRasters, envir = environment()) 
data(GLOBALMAP, envir = environment()) 

##-- Disaggregate habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)

##-- Merge counties for practical reasons
COUNTIES <- COUNTIES %>%
  mutate(id = case_when(
    county %in% c("Akershus","Agder","Buskerud","Vestfold","Oslo","Østfold","Telemark") ~ "NO1",
    county %in% c("Innlandet","Møre og Romsdal","Trøndelag") ~ "NO2",
    county %in% c("Jämtlands","Västernorrlands","Västerbottens","Dalarnas","Gävleborgs") ~ "SE1",
    county %in% c("Uppsala","Västmanlands","Stockholms","Södermanlands") ~ "SE2",
    county %in% c("Blekinge","Örebro","Östergötlands","Jönköpings","Kronobergs","Kalmar","Skåne","Gotlands") ~ "SE3",
    county %in% c("Västra","Värmlands","Hallands") ~ "SE4")) 
# %>%
#   dplyr::group_by(id) %>%
#   dplyr::summarise()

##-- PLOT CHECK 
if(plot.check){
  COUNTIES %>% 
    group_by(id) %>%
    summarize() %>%
    ggplot(.) +
    geom_sf(aes(fill = id)) +
    geom_sf_label(aes(label = id))
}

##-- CREATE STUDY AREA POLYGON BASED x AND y EXTENTS
# myStudyArea.extent <- st_bbox(extent(x.extent, y.extent))
# st_crs(myStudyArea.extent) <- st_crs(COUNTRIES)
studyArea <- st_crop( COUNTIES,
                      xmin = x.extent[1], xmax = x.extent[2],
                      ymin = y.extent[1], ymax = y.extent[2]) %>%
  st_collection_extract(., "POLYGON")  



## ------   2. LOAD ROVBASE FILES ------ 

##-- NGS data
# DNA <- read.csv( file.path(data.dir, "RIB22042025133456403_wolfDNA.csv"), fileEncoding = "latin1")
# colnames(DNA) <- translateForeignCharacters( dat = colnames(DNA), dir.translation = dir.analysis)
DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                        extension = ".xls",
                                        pattern = "DNA")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>%
  ##-- Turn potential factors into characters 
  dplyr::mutate(across(where(is.factor), as.character)) %>%
  ##-- Initial filters
  dplyr::filter(
    ##-- Filter to the focal species
    Species %in% norSpecies,
    # [CHECK] Should we do this for all species now ???
    # ##-- Filter dead recoveries (HB for the last wolverine analysis)
    # !substr(RovbaseID_sample,1,1) %in% "M"
  ) %>%
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
    Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown", Sex),
    #Sex = ifelse(is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% "Hunn", "female", Sex),
    Sex = ifelse(Sex %in% "Hann", "male", Sex))
# [CHECK] Should we filter for years here ???
# %>%
# ##-- Filter to the focal years
# dplyr::filter(., Year %in% years)


## Dead Recoveries from RovBase
# DEAD <- read.csv( file.path(data.dir, "RIB22042025133534832_wolfDEAD.csv"),
#                   fileEncoding = "latin1")
# colnames(DEAD) <- translateForeignCharacters(dat = colnames(DEAD), dir.translation = dir.analysis )
##-- Load raw excel file imported from rovbase 
DR <- suppressWarnings(readMostRecent( path = data.dir,
                                       extension = ".xls",
                                       pattern = "dead")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>%
  ##-- Initial filters
  dplyr::filter(
    ##-- Filter to the focal species
    Species %in% norSpecies) %>%
  ##-- Remove any duplicates
  dplyr::distinct(., .keep_all = TRUE) %>%
  ##-- Turn potential factors into characters 
  dplyr::mutate(across(where(is.factor), as.character)) %>%
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
    Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown" , Sex),
    #Sex = ifelse(is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% "Hunn", "female", Sex),
    Sex = ifelse(Sex %in% "Hann", "male", Sex),
    ##-- Identify legal deaths
    Legal = grepl(paste(legal.dead, collapse = "|"), Death_cause))
# [CHECK] Should we filter for years here ???
# %>%
# ##-- Filter to the focal years
# dplyr::filter(., Year %in% years)


## Wolves infos from Micke
INDIVIDUAL_ID <- suppressWarnings(readMostRecent( path = data.dir,
                                                  extension = ".xls",
                                                  pattern = "_ID Grouping")) %>%
  dplyr::rename(., IdSimplified = "ROVBASE_IndividID")
  
# INDIVIDUAL_ID <- read.csv( file.path(data.dir, "220512_ID Grouping 2006-2021.csv"), fileEncoding = "latin1")  
# colnames(INDIVIDUAL_ID) <- translateForeignCharacters(dat = colnames(INDIVIDUAL_ID), dir.translation = dir.analysis)


## THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2022/23.
Pack_ID2023 <- suppressWarnings(readMostRecent( path = data.dir,
                                                extension = ".xls",
                                                pattern = "Genetiskt ID RM")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>% 
  ##-- Add some columns
  dplyr::mutate( 
    Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
    Sex = ifelse(Sex %in% c("Hane","Hann"), "male", Sex))
# read.csv( file.path(data.dir, "Genetiskt ID RM vargar 2223 Bilaga 4_ØF.csv"), fileEncoding = "latin1")  
# colnames(Pack_ID2023) <- translateForeignCharacters(dat = colnames(Pack_ID2023), dir.translation = dir.analysis)
# Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Tispe"] <- "Hunn"
# Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Hane"] <- "Hann"
# Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Tik"] <- "Hunn"


## THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2023/24.
Pack_ID2024 <- suppressWarnings(readMostRecent( path = data.dir,
                                                extension = ".xls",
                                                pattern = "Bilaga_11")) %>%
  ##-- Rename columns to facilitate manipulation
  dplyr::rename(., any_of(rename.list)) %>% 
  ##-- Add some columns
  dplyr::mutate( 
    Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
    Sex = ifelse(Sex %in% c("Hane","Hann"), "male", Sex))
# read.csv( file.path(data.dir, "Bilaga_11.4_240424_ØF to Cyril.csv"), fileEncoding = "latin1")  
# colnames(Pack_ID2024) <- translateForeignCharacters(dat = colnames(Pack_ID2024), dir.translation = dir.analysis)
# Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Tispe"] <- "Hunn"
# Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Hane"] <- "Hann"
# Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Tik"] <- "Hunn"


## THIS IS THE PACK ID SENT BY Øystein.
Pack_ID2025 <- suppressWarnings(readMostRecent( path = data.dir,
                                                extension = ".xls",
                                                pattern = "estimates2025FromOystein"))
## Here we need to recreate the sex columns as Oystein gave me a list of ids only (losing the sex)
Pack_ID2025$Sex <- apply(Pack_ID2025[ ,c("Sex1","Sex2","Sex3","Sex4")], 1, function(x) x[!is.na(x)][1]) 
Pack_ID2025 <- Pack_ID2025 %>% 
  dplyr::mutate( 
    Sex = ifelse(Sex %in% "Ukjent" | is.na(Sex), "unknown", Sex),
    Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
    Sex = ifelse(Sex %in% c("Hane","Hann"), "male", Sex))
# read.csv( file.path(data.dir, "RovbaseID for Rovquant estimates2025FromOystein.csv"), fileEncoding = "latin1")  
# colnames(Pack_ID2025) <- translateForeignCharacters(dat = colnames(Pack_ID2025), dir.translation = dir.analysis)



## ------   3. SEARCH EFFORT DATA ------ 

# ##-- Combine all GPS tracks
# TRACKS <- rbind(
#   sf::read_sf(file.path(data.dir, "TRACKS/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250422.shp")),
#   sf::read_sf(file.path(data.dir, "TRACKS/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422.shp"))) %>%
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
# df <- data.frame( Dato = TRACKS$Dato,
#                   Year = TRACKS$Year,
#                   Person = TRACKS$Person,
#                   Length = TRACKS$Length,
#                   Centroidx = TRACKS$Centroidx)
# dupIDs <- which(duplicated(df))
# dupLength <- TRACKS$Length[duplicated(df)]
# TRACKS <- TRACKS[-dupIDs, ]
# 
# ## --PLOT CHECK
# if(plot.check){
#   ## Total length of tracks searched
#   TRACKS %>%
#     group_by(Year)%>%
#     summarize(TotLength = sum(Length)) %>%
#     barplot( ., ylab = "Sum length tracks")
#    
#   ## Check number of duplicated tracks removed
#   dup <- unlist(lapply(dupIDs, length))
#   names(dup) <- years
#   barplot(dup, ylab = "Number of duplicated tracks")
#   
#   ## Check length of duplicated tracks removed
#   dupdist <- unlist(lapply(dupDist,sum))
#   names(dupdist) <- years
#   barplot(dupdist, ylab = "Distance of duplicated tracks")
# }
# 
# ##-- Save TRACKS
#save( TRACKS, file = file.path(working.dir, "data", "TRACKSSouthSweden2014202540NotSimplifiedSF.RData"))
load(file.path(working.dir, "data", "TRACKSSouthSweden2014202540NotSimplifiedSF.RData"))



## -----------------------------------------------------------------------------

## ------ II. CREATE OPSCR DATA ------

## ------   1. CLEAN & FILTER NGS DATA ------ 

## ------     1.1. CLEAN NGS & DEAD RECOVERY DATA ------ 

myCleanedData.sp <- CleanDataNew3sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(DATA$sampling.months)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)

##-- make a simplified column to match the rovbase id given by Oystein in Linn's file
myCleanedData.sp$IdSimplified <- unlist(lapply(strsplit(as.character(myCleanedData.sp$Id), " "), function(x) x[1]))

##-- OVERWRITE GENDER FROM MICKE'S DATA WHEN AVAILABLE
micke.sex <- as.character(unlist(lapply(myCleanedData.sp$Id, function(i) INDIVIDUAL_ID[as.character(INDIVIDUAL_ID$IdSimplified)==i,"Sex"][1])))
micke.sex[micke.sex %in% "0"] <- NA
micke.sex[micke.sex %in% names(table(micke.sex))[3]] <- NA
micke.sex[micke.sex %in% "Hona"] <- "female"
micke.sex[micke.sex %in% "Hane"] <- "male"
table(!is.na(micke.sex))
new.sex <- ifelse(!is.na(micke.sex), as.character(micke.sex), as.character(myCleanedData.sp$Sex))
table(myCleanedData.sp$Sex, new.sex)
myCleanedData.sp$Sex <- new.sex
table(myCleanedData.sp$Sex, new.sex)

##-- OVERWRITE GENDER FROM PACK COMPOSITION (FROM LINN's file 2023-24)
tab <- list()
##-- check the sex in the pair data given by Linn and assign the sex to all detections 
for(i in 1:length(Pack_ID2023$Sex)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2023$Rovbase.ID[i]])
  #Overwrite sex 
  # if(length(tab[[i]])>1){print(tab[[i]])}
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2023$Rovbase.ID[i]] <- Pack_ID2023$Sex[i]
}
for(i in 1:length(Pack_ID2024$Sex)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2024$Rovbase.ID[i]])
  #Overwrite sex 
  # if(length(tab[[i]])>1){print(tab[[i]])}
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2024$Rovbase.ID[i]] <- Pack_ID2024$Sex[i]
}
for(i in 1:length(Pack_ID2025$Sex)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2025$IndividID[i]])
  # if(length(tab[[i]])>1){print(tab[[i]])}
  #Overwrite sex 
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2025$IndividID[i]] <- Pack_ID2025$Sex[i]
}



## ------     1.2. FILTER DATA FOR SEX -----

myFullData.sp <- FilterDatasf( 
  myData = myCleanedData.sp,
  dead.recovery = T,
  sex = DATA$sex, 
  setSex = T)



## ------     1.3. REMOVE INDVIDUALS THAT DIED TWICE ------ 

duplicatedDeath <- NULL
for(i in myFullData.sp$IdDoubleDead){
  tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
  if(length(tmp)==0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-1]}
  duplicatedDeath <- c(duplicatedDeath, tmp)
}#i  
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]




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
  sf::st_filter( .,st_as_sfc(st_bbox(studyArea)), .predicate = st_intersects)


## ------     6.2. DEAD RECOVERY DATA -----

data.dead <- myFullData.sp$dead.recovery %>%
  dplyr::filter(
    ##-- Subset to years of interest
    Year %in% years,
    ##-- Subset to sex of interest
    Sex %in% sex) %>%
  ##-- Filter based on space 
  sf::st_filter( .,st_as_sfc(st_bbox(studyArea)), .predicate = st_intersects)

# ## Remove all alive detections outside of the study area extent 
# data.alive <- data.alive[!is.na(as.numeric(st_intersects(data.alive, st_as_sfc(myStudyArea.extent)))), ]
# ## Remove all dead recoveries outside of the study area polygon 
# data.dead <- data.dead[!is.na(as.numeric(st_intersects(data.dead, st_as_sfc(myStudyArea.extent)))), ]
# ## Remove all alive detections outside of the sampling period 
# data.alive <- data.alive[ data.alive$Month %in% unlist()&  data.alive$Year %in% unlist(DATA$years), ]
# ## Remove all dead recoveries outside of the sampling period 
# data.dead <- data.dead[ data.dead$Year %in% unlist(DATA$years),] 



## ------     1.6. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------ 

## ------       1.6.1. ASSIGN SAMPLES TO TRACKS  ------ 

## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
data.alive$TrackRovbsID <- NA
data.alive$TrackDist <- NA

TRACKSSimple_sf <- list()
for(t in 1:nYears){
  TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]
}

## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
dnatemp <- st_as_sf(data.alive)
## CREATE A BUFFER AROUND EACH DETECTION
tmp <-  st_buffer(dnatemp, dist=750)

for(i in 1:nrow(data.alive)){
  # INTERSECT POINT WITH TRACKS,
  t <- which(years %in% tmp[i, ]$Year)
  whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(data.alive$Date[i]))
  tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])
  
  if(nrow(tmpTRACKS)==0){next}
  
  # FIND THE CLOSEST TRACK
  dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
  
  # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
  # IF NO MATCHING DATE ASSIGN TO NA?
  if(length(dist)==0){
    data.alive$TrackRovbsID[i] <- NA
    data.alive$TrackDist[i] <- NA
  }
  # IF MATCHING DATE ASSING TO THAT TRACK
  if(length(dist)==1){
    data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
    data.alive$TrackDist[i] <- dist
  }
  # IF SEVERAL MATCHING DATES ASSING TO THE CLOSEST OF THE MATCHING TRACKS
  if(length(dist)>1){
    data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
    data.alive$TrackDist[i] <- min(dist)
  }
  print(i)
}

##-- SAVE FOR FASTER LOADING
save( myFilteredData.sp,
      file = file.path(working.dir, "data", "_myFilteredData.sp.RData"))
load(file.path(working.dir, "data", "_myFilteredData.sp.RData"))



## ------       1.6.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC AND STRUCTURED ------ 

distanceThreshold <- 500

## Proevetype columns was replaced by two columns, merging them now...
data.alive$Proevetype <-  ifelse(
  data.alive$Annen.innsamler...Rolle %in% "", 
  data.alive$Samlet.selv...Rolle,
  data.alive$Annen.innsamler...Rolle)

whichStructured <- data.alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(data.alive$TrackRovbsID) &
  data.alive$TrackDist <= distanceThreshold
myFilteredData.spStructured <- data.alive[whichStructured,]
myFilteredData.spOthers <- data.alive[!whichStructured,]

## CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
nrow(myFilteredData.spStructured) + nrow(myFilteredData.spOthers)
nrow(data.alive)

## Check number of opp vs. struc each year 
data.alive$TrackDistCat <- ifelse(data.alive$TrackDist > 500, 0, 1)
data.alive$TrackDistCat[is.na(data.alive$TrackDistCat)] <- 0

table( data.alive$Proevetype,
       data.alive$Year,
       data.alive$TrackDistCat)



## ------     1.7. DATA SUMMARY ------
## ------       1.7.1. PLOT CHECKS ------ 

## PROPORTION SAMPLES STRUCTURED/OTHERS
pdf(file = file.path(working.dir, "figures", "ProportionStucturedOther.pdf"))
par(mfrow = c(1,1), mar = c(4,4,3,2))
barplot(rbind(table(myFilteredData.spStructured$Year),
              table(myFilteredData.spOthers$Year)),
        beside = T, ylim = c(0,2000),
        col = c(grey(0.2),grey(0.8)),
        ylab = "Number of samples")
abline(h = seq(0,2000,by=500), lty = 2, col = grey(0.8))
title(main = "500m threshold")
legend("topleft",
       fill = c(grey(0.2),grey(0.8)),
       legend = c("Structured","Other"))
dev.off()

## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
tmp <- data.alive[data.alive$Proevetype %in% 
                                 c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
tab <- table(tmp$Year, tmp$TrackRovbsID, useNA ="always" )

## MAP  SAMPLES STRUCTURED OTHERS
pdf(file = file.path(working.dir, "figures", "MapStucturedOther.pdf"))
for(t in 1:nYears){
  par(mar = c(0,0,3,0), mfrow = c(1,3))
  tmp1 <- tmp[tmp$Year%in% years[t],]
  tmpNoTracks <- tmp1[is.na(tmp1$TrackRovbsID), ]
  tmpTracks <- tmp1[!is.na(tmp1$TrackRovbsID), ]
  
  plot(st_geometry(studyArea), main="Structured with track")
  plot(st_geometry(tmpTracks), pch=21, col="black", cex=1,bg="red",add=T)
  
  plot(st_geometry(studyArea), main="Structured without track")
  plot(st_geometry(tmpNoTracks), pch=21, col="black", cex=1,bg="blue",add=T)
  
  tmpOpp <- data.alive[!data.alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
  tmpOpp <- tmpOpp[tmpOpp$Year%in% years[t],]
  
  plot(st_geometry(studyArea), main="Other samples")
  plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
  mtext(years[t],adj = -0.8,padj = 1)
}
barplot(tab[,which(is.na(colnames(tab)))]/rowSums(tab),main="% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 
dev.off()

## OVERALL MAP DETECTION DEAD RECOVERIES MAP
pdf(file = file.path(working.dir, "figures", "OverallDetectionsDeadRecoveries.pdf"))
plot(st_geometry(GLOBALMAP))
plot(st_geometry(studyArea),add=T)
plot(st_geometry(myFullData.sp$alive), pch=16, col="red", cex=0.3,add=T)
plot(st_geometry(myFullData.sp$dead.recovery),pch=16, col="blue", cex=0.3,add=T)
mtext(paste("Live detections", length(myFullData.sp$alive),
            "; ID:", length(unique(myFullData.sp$alive$Id))),
      line = +1)
mtext(paste("Dead recovery:",length(myFullData.sp$dead.recovery)))
dev.off()


## ------       4.1. DETECTIONS ALIVE ------ 

## NUMBER OF INDIVIDUALS DETECTED ALIVE
length(unique(data.alive$Id))

## NUMBER OF INDIVIDUALS DETECTED ALIVE & RECOVERED DEAD
sum(unique(data.alive$Id) %in% unique(data.dead$Id))

## NUMBER OF INDIVIDUALS DETECTED/YEAR/COUNTRY
table.id <- table(data.alive$Id,data.alive$Year, data.alive$Country)
apply(table.id, c(2,3), function(x) sum(x>0))

## NUMBER OF DETECTIONS/YEAR/COUNTRY
countrytab <- table(data.alive$Year, data.alive$Country)
countrytab 

## NUMBER OF DETECTIONS/YEAR/COUNTRY/SEX
sex_countrytab <- table(data.alive$Year, data.alive$Country, data.alive$Sex)
sex_countrytab

## PROPORTION OF DETECTIONS PER YEAR/COUNTRY/SEX
sex_countrytab_prop <- sex_countrytab
for(i in 1:dim(sex_countrytab)[3]){
  sex_countrytab_prop[,,i] <- sex_countrytab_prop[,,i]/countrytab
}
sex_countrytab_prop



## ------       4.2. DEAD RECOVERIES ------ 

## NUMBER OF INDIVIDUALS RECOVERED
length(unique(data.dead$Id))

## NUMBER OF DEAD RECOVERIES/YEAR/COUNTRY
table(data.dead$Year, data.dead$Country)

## MORTALITY CAUSES
unique(as.character(data.dead$DeathCause))
unique(as.character(data.dead$DeathCause_2))  
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))

## DEFINE LEGAL MORTALITY
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])



## ------       4.3. PLOTS ------ 

legal.death <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]
Other.death <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]

table(legal.death$Year)
table(Other.death$Year)

## PLOT CHECK
if(plot.check){
  pdf(file = file.path(working.dir,"figures","DetectionsDeadRecoveries.pdf"))
  
  ## PLOT ALL DETECTIONS
  for(t in 1:nYears){
    plot(habitat$habitat.r, main = years[t]) 
    plot(st_geometry(detectors$main.detector.sp), add=T, pch=16, cex=0.1)
    plot(st_geometry(data.alive[data.alive$Year == years[t], ]), 
         pch=16,col="red", cex=0.7, add=T)
    mtext(paste("Live detections",
                nrow(data.alive[data.alive$Year == years[t], ]),
                "; ID:",
                nrow(unique(data.alive[data.alive$Year == years[t], ]$Id))),
          line = +1)
    
    plot(st_geometry(data.dead[data.dead$Year == years[t], ]),
         pch=16,col="blue", cex=0.7,add=T)
    mtext(paste("Dead recovery:",nrow(data.dead[data.dead$Year == years[t], ])))
    plot(st_geometry(GLOBALMAP), add=T) 
  }#t
  
  ## PLOT NUMBER OF MORTALITY EVENTS, LEGAL/OTHERS 
  barplot(table(legal.death$Country, legal.death$Year), main = "Legal causes")
  barplot(table(Other.death$Country, Other.death$Year), main = "Other causes")
  legend("topleft", fill = c(grey(0.3), grey(0.7)), legend = c("NOR","SWE"))
  dev.off()  
}

## EXPORT NGS DATA 
if(DATA$sex == "Hann"){
  assign("myFilteredData.spM", myFilteredData.sp)
  assign("myFilteredData.spOthersM", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
  
  save(myFilteredData.spM, 
       myFullData.spM,
       myFilteredData.spOthersM,
       myFilteredData.spStructuredM,
       file = file.path(working.dir, "data", "NGSData.RData"))
} else {
  assign("myFilteredData.spF", myFilteredData.sp)
  assign("myFilteredData.spOthersF", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
  
  save(myFilteredData.spF,
       myFullData.spF,
       myFilteredData.spOthersF,
       myFilteredData.spStructuredF,
       file = file.path(working.dir, "data", "NGSData.RData"))
}



## ------   2. GENERATE HABITAT ------

## ------     2.1. GENERATE HABITAT CHARACTERISTICS ------ 

myHabitat.list <- MakeHabitatFromRastersf( 
  poly = studyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = buffer.size,                               
  plot.check = T)

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- makeHabitatFromRaster(
  poly = studyArea,
  habitat.r = habRaster,
  buffer = habitat$buffer,
  plot.check = FALSE) %>%
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

## PLOT CHECK
if(plot.check){
  plot(myHabitat.list$habitat.r)
  plot(st_geometry(studyArea), add=T)
  plot(st_geometry(COUNTRIES), add=T)
}



## ------     2.2. GENERATE HABITAT-LEVEL COVARIATES ------ 

## KERNEL OF INDIVIDUALS IN PAIRS 
kern <- list()
habDens <- matrix(NA, nrow = numHabWindows, ncol = nYears)
IDS <- unlist(lapply(strsplit(as.character(myFullData.sp$alive$Id), " "), function(x)x[1])) 
for(t in 1:nYears){
  id.fam <- which(y.obsALL[ ,as.character(years[t]-1)] %in% c(3,4), arr.ind = T)
  
  m.xy <- matrix(NA, nrow = length(id.fam), ncol = 2)
  colnames(m.xy) <- c("x","y")
  for(i in 1:length(id.fam)){
    tmp <- myFullData.sp$alive[IDS == row.names(y.obsALL)[i], ]
    m.xy[i, ]<- colMeans(st_coordinates(tmp))
  }#i
  if(sum(is.na(m.xy[ ,1])) > 0){
    m.xy <- m.xy[!is.na(m.xy[ ,1]), ]
  }
  locationsFamily <- st_as_sf( as.data.frame(m.xy),
                               coords = c("x","y"),
                               crs = st_crs(habitat$habitat.sp))
  locationsFamily$id <- rep(1, nrow(locationsFamily))
  kern[[t]] <- raster(estUDm2spixdf(kernelUD( as(locationsFamily[ ,"id"], "Spatial"), 
                                              h = 15000,
                                              grid = as(habitat$habitat.r, 'SpatialPixels'))))
  habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ] == 1])
}

## PLOT CHECK 
for(t in 1:nYears){
  plot(kern[[t]], main = years[t])
  plot(habitat$habitat.poly$geometry, add = T, col = NA)
}#t




## ------   3. GENERATE DETECTORS ------ 

## ------     3.1. GENERATE DETECTORS CHARACTERISTICS ------ 

habitat.subdetectors <- disaggregate( myHabitat.list$habitat.rWthBuffer,
                                      fact = res(myHabitat.list$habitat.r)[1]/subdetector.res)

myDetectors <- MakeSearchGridsf(
  data = habitat.subdetectors,
  resolution = detector.res,
  div = (detector.res/subdetector.res)^2,      
  plot = FALSE,
  fasterize = TRUE)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(mymain.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS and NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(mymain.detector.sp)
n.trials <- as.vector(table(mydetector.sp$main.cell.id))

## PLOT CHECK
if(plot.check){
  par(mfrow = c(1,2))
  plot(st_geometry(studyArea), main = "Detectors Alive")
  plot(st_geometry(mymain.detector.sp), col = "red", pch = 16, cex = 0.1, add = T)
  plot(st_geometry(GLOBALMAP), add = T)
  
  plot(st_geometry(studyArea), main = "Detectors Dead")
  plot(st_geometry(mymain.detector.sp), col = "red", pch = 16, cex = 0.1, add = T)
  plot(st_geometry(GLOBALMAP), add = T)
}



## ------     3.2. GENERATE DETECTOR-LEVEL COVARIATES ------ 

## ------       3.2.1. EXTRACT COUNTRIES ------ 

dist <- st_distance(mymain.detector.sp, COUNTRIES, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))



## ------       3.2.2. EXTRACT COUNTIES ------ 

dist <- st_distance(mymain.detector.sp, COUNTIES, by_element = F )

## MERGE COUNTIES
COUNTIES$id <- 1:nrow(COUNTIES)
#NO1
COUNTIES$NAME_1[c(1,2,4,33,22,39)]
COUNTIES$id[c(1,2,4,33,22,39)] <- 2
#NO2
COUNTIES$NAME_1[c(20,10,16,28,17)]
COUNTIES$id[c(20,10,16,28,17)] <- 10#35
#SE1
COUNTIES$NAME_1[c(12,36,35,5,8)]
COUNTIES$id[c(12,36,35,5,8)] <- 5
#SE2
COUNTIES$NAME_1[c(31,37,26,27)]
COUNTIES$id[c(31,26,27)] <- 4
#SE3
COUNTIES$NAME_1[c(3,21,40,13,15,14,24,7)]
COUNTIES$id[c(3,21,13,15,14,24,7)] <- 3
#SE4
COUNTIES$NAME_1[c(38,34,9)]
COUNTIES$id[c(34,9)] <- 9

## CONVERT TO FACTOR & BACK
COUNTIES$id <- as.numeric(as.factor(COUNTIES$id))

## ASSIGN COUNTIES TO DETECTORS.
detCounties1 <- apply(dist,1, function(x) which.min(x))
detCounties1 <- COUNTIES$id[detCounties1]
detCounties <- as.numeric(as.factor(detCounties1))
table(detCounties)

## CREATE A VECTOR TO KEEP THE ORIGINAL COUNTY ID
detCounties.original <- 0
for(i in 1: max(detCounties)){
  detCounties.original[i] <- detCounties1[which(detCounties==i)][1]
}

## PLOT CHECK 
COUNTIESplot <- st_simplify(COUNTIES, dTolerance = 500) %>%
  #st_intersection(., studyArea) %>%
  group_by(id) %>%
  summarize()

ggplot(st_simplify(COUNTIESplot,dTolerance = 500)) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))

col <- rainbow(length(unique(detCounties)))
plot( st_geometry(COUNTIESplot))
plot( st_geometry(mymain.detector.sp),
      col = col[detCounties], pch = 16,
      cex = 0.8, add = T)



## ------       3.2.3. EXTRACT GPS TRACKS LENGTHS ------ 

## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(mymain.detector.sp),
                                      rep(1,nrow(mymain.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(studyArea)
detectorGrid$id <- 1:nrow(detectorGrid)
plot(st_geometry(detectorGrid))

## CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
for(t in 1:nYears){
  TRACKSst <- TRACKS_YEAR[[t]]
  
  intersection <- st_intersection(detectorGrid, TRACKSst) %>%
    dplyr::mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(transect_L = sum(LEN))    ## Get total length searched in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
}

## Plot check 
par(mar = c(0,0,0,0))
plot(st_geometry(mymain.detector.sp[as.numeric(intersection$id), ]),
     cex = DoScale(detTracks[ ,t], l = 0, u = 2), pch = 16)



## ------       3.2.4. EXTRACT DISTANCES TO ROADS ------ 

##-- Load map of distance to roads (1km resolution)
DistAllRoads <- raster::raster(file.path(data.dir, "Roads/MinDistAllRoads1km.tif"))

##-- Fasterize to remove values that fall in the sea
r <- fasterize::fasterize(sf::st_as_sf(COUNTRIES), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- raster::crop(DistAllRoads, studyArea)
rm(list = c("r"))


##-- AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- raster::aggregate( 
  x = DistAllRoads,
  fact = detectors$resolution/raster::res(DistAllRoads),
  fun = mean)

##-- EXTRACT ROAD DISTANCE FOR EACH DETECTOR
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


## PLOT CHECK
if(plot.check){
  plot(st_geometry(mymain.detector.sp),cex=DoScale(detRoads),pch=16)
}









## ------       3.2.5. EXTRACT DAYS OF SNOW ------ 

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
## UPDATE"!!!!!
SNOW <- stack(file.path(data.dir, "Snow/AverageSnowCoverModisSeason2014_2025_Wolf.tif"))

## RENAME THE LAYERS
names(SNOW) <- paste(2014:2024,(2014:2024)+1, sep = "_")

## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))

## EXTRACT SNOW 
detSnow <- matrix(0, nrow = n.detectors, ncol = nYears)
det.sptransf <- st_transform(mymain.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp
if(plot.check){
  plot(st_geometry(mymain.detector.sp), cex=DoScale(detSnow[,3]),pch=16)
}

## still some NA... Increase buffer again 
isna <- which(is.na(detSnow),arr.ind = T)
isna <- unique(isna[,1])
tmp.list <- raster::extract(SNOW, det.sptransf[isna,], buffer=35000)
detSnow[isna,1:nYears] <- unlist(lapply(tmp.list, function(x) colMeans(x, na.rm=T)))



## ------       3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------ 

## ------         3.2.6.1. SKANDOBS ------ 

skandObs <- read_xlsx(file.path(dir.data, "Richard_Biscof_Skandobs_2012_2025dd.xlsx"))
colnames(skandObs) <- translateForeignCharacters( dat = colnames(skandObs),
                                                  dir.translation = dir.analysis)

## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))

## MAKE IT SPATIAL 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(studyArea))

## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist())
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) #--- need to change for other species
skandObs <- skandObs[subset,] 

## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(myHabitat.list$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1,]

subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace,] 
plot(st_geometry(habitat.rWthBufferPol))
plot(st_geometry(skandObs),col="red",add=T)

## RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate(habitat.subdetectors, fact=(detector.res/subdetector.res))
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
plot(r.skandObsSamplesBinary[[t]])

## PLOT CHECK 
if(plot.check){
  ## SUMMARY SKANDOBS
  pdf( file = file.path(working.dir,"figures","skandObs.pdf"),
       width = 10)
  barplot(table(skandObs$monitoring.season))
  barplot(table(skandObs$month), xlab = "Months")
  barplot(table(skandObs$activity), cex.names = 0.7)
  barplot(table(skandObs$species))
  
  ## MAPS 
  par(mar = c(0,0,2,0))
  for(t in 1:nYears){
    plot(st_geometry(studyArea), main= years[t])
    plot(st_geometry(skandObs[skandObs$monitoring.season %in% years[t],  ]), pch=16, col="red", cex=0.1)
  }
  dev.off()
}



## ------         3.2.6.2. ROVBASE ------ 

## GET ALL SAMPLES COLLECTED
# files need to be loaded separately and binded. 
rovbaseObs1 <- read_xlsx(file.path(data.dir, "RIB15042025111042561.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(data.dir, "RIB15042025112748412.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(data.dir, "RIB15042025112853603.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3)

colnames(rovbaseObs) <- translateForeignCharacters( dat = colnames(rovbaseObs),
                                                    dir.translation = dir.analysis)
rovbaseObs$Proevetype <- translateForeignCharacters( dat = rovbaseObs$Proevetype,
                                                     dir.translation = dir.analysis)

rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Oest (UTM33/SWEREF99 TM)`), ]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

## DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf( rovbaseObs,
                           coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(studyArea)

## SUBSET THE DATA 
filter <- list(
  species = "Ulv",
  type = c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)","Saliv/Spytt"),
  month = unlist())

## SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset,] 

## SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
subset <- rovbaseObs.sp$`Art (Proeve)`%in% filter$species #& !is.na(rovbaseObs.sp$`RovbaseID (Analyse)`) 
rovbaseObs.sp <- rovbaseObs.sp[subset,] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace,] 

## RASTERIZE 
r.detector <- aggregate(habitat.subdetectors, fact=(detector.res/subdetector.res))
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

## PLOT CHECK
if(plot.check){
  
  pdf(file = file.path(working.dir,"figures","mapStructuredOthers.pdf"))
  for(t in 1:nYears){
    year = years[t]
    tmpOthers <- myFilteredData.spOthers[myFilteredData.spOthers$Year%in%year, ]
    tmpStruct <- myFilteredData.spStructured[myFilteredData.spStructured$Year%in%year, ]
    
    par(mfrow = c(2,2), mar = c(0,0,5,0))
    plot( r.OtherSamplesBinary[[t]],
          main = paste(year,"\n Rovbase Samples Other"), box = F, axes = F)
    plot( st_geometry(tmpOthers),
          pch = 16, col = "blue", bg = "blue", cex = 0.6, add = T)
    plot( r.OtherSamplesBinary[[t]],
          main = paste (year,"\n Rovbase Samples Structured"), box = F, axes = F)
    plot( st_geometry(tmpStruct),
          pch = 16, col = "red", bg = "red", cex = 0.6, add = T)
    
    plot( r.skandObsSamplesBinary[[t]],
          main = paste(year,"\n SkandObs Other"),
          box = F, axes = F)
    plot( st_geometry(tmpOthers),
          pch = 16, col = "blue", bg = "blue", cex = 0.6, add = T)
    plot( r.skandObsSamplesBinary[[t]],
          main = paste(year,"\n SkandObs Structured"), box = F, axes = F)
    plot( st_geometry(tmpStruct),
          pch = 16, col = "red", bg = "red", cex = 0.5, add = T)
  }
  dev.off()
}


## ------         3.2.6.3. COMBINE ROVBASE & SKANDOBS ------ 

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <-    1
}

## PLOT CHECK
if(plot.check){
  for(t in 1:nYears){
    par(mfrow = c(1,3), mar = c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]], main = years[t])
    plot(r.skandObsSamplesBinary[[t]])
    plot(r.SkandObsOtherSamplesBinary[[t]])
  }  
}



## ------         3.2.6.4. SMOOTH THE BINARY MAP ------ 

## we tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
habOwin <- as.owin(as.vector(extent(r.detector)))
ds.list <- lapply(years, function(y){
  ## ROVBASE  
  pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y, ]
  ## SKANDOBS
  pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y, ] )
  ## SMOOTH & RASTERIZE
  p <- ppp(pts[,1], pts[,2], window = habOwin)
  ds <- density(p, adjust = 0.02) ##-- change bandwith (smoothing) with "adjust
  ds <- raster(ds)
  ds <- ds1 <- raster::resample(ds, r.detector) 
  threshold <- 0.1 / prod(res(ds)) ##-- number per 1 unit of the projected raster (meters)
  ds1[] <- ifelse(ds[]<threshold,0,1)
  ds1 <- mask(ds1, habitat.rWthBufferPol)
  ds <- mask(ds, habitat.rWthBufferPol)
  return(list(ds,ds1))
})

ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))
names(ds.brick) <- years

## PLOT CHECK
if(plot.check){
  par(mfrow = c(1,3))
  plot(r.SkandObsOtherSamplesBinary[[t]], main = "Raw Binary", axes = F, box = F)
  plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
  plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
}



## ------         3.2.6.5. ASSIGN THE COVARIATE ------ 

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract(r.SkandObsOtherSamplesBinary, mymain.detector.sp)
colSums(detOtherSamples)



## ------       3.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

## CHECK IF CONTAINS NAs
if(sum(is.na(detSnow)) > 0 | sum(is.na(detRoads)) > 0 | sum(is.na(detTracks)) > 0){
  print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")
}


## PLOT COVARIATES 
if(plot.check){
  pdf(file = file.path(working.dir, "figures", "SpatialCovariates.pdf"))
  
  ## Tracks
  max <- max(detTracks)
  cuts <- seq(min(detTracks), max, length.out = 100) ## set breaks
  col <- rev(terrain.colors(100))
  for(t in 1:nYears){
    detectorGrid <- rasterFromXYZ(cbind(st_coordinates(mymain.detector.sp),
                                        rep(1,nrow(mymain.detector.sp))))
    id <- which(detectorGrid[] %in% 1)
    detectorGrid[id] <- detTracks[ ,t]
    
    plot( detectorGrid,
          breaks = cuts, col = col,
          main = paste("Tracks", years[t]),
          legend = F)
    plot( detectorGrid,
          breaks = cuts, col = col,
          legend.only = TRUE,
          legend.width = 2,
          axis.args = list(
            at = round(seq(min(detTracks), max, length.out = 8),digits = 0),
            labels = round(seq(min(detTracks), max, length.out = 8),digits = 0), 
            cex.axis = 0.6),
          legend.args = list(
            text = '',
            side = 4,
            font = 2,
            line = 2.5,
            cex = 0.8))
    plot( st_geometry(data.alive[data.alive$Year == years[t], ]),
          pch = 16, col = "red", cex = 0.2, add = T)
    
    plot( st_geometry(myHabitat.list$buffered.habitat.poly),
          add = T, border = "grey")
  }#t
  
  ## Snow
  for(t in 1:nYears){
    plot(SNOW[[t]], main=paste("Snow", years[t]))
    plot(st_transform(myHabitat.list$buffered.habitat.poly, st_crs(SNOW))$geometry, add=T)
  }#t
  
  ## Roads
  plot(DistAllRoads, main = "Roads")
  plot(st_geometry(myHabitat.list$buffered.habitat.poly), add = T)
  dev.off()
}



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

# scaledCoords <- scaleCoordsToHabitatGrid(
#   coordsData = detectors.xy,
#   coordsHabitatGridCenter = habitat.xy)
#
# scaledLowUpCoords <- getWindowCoords(
#   scaledHabGridCenter = scaledCoords$coordsHabitatGridCenterScaled,
#   scaledObsGridCenter = scaledCoords$coordsDataScaled)



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
  dmax = detectors$maxDist/habitat$resolution,
  resizeFactor = resize.factor,
  plot.check = F)



## ------   4. GENERATE y DETECTION ARRAYS ------ 

## ------     4.1. DATA SUMMARY ------ 

## ------       4.1.1. DETECTIONS ALIVE ------ 

## NUMBER OF INDIVIDUALS DETECTED ALIVE
length(unique(data.alive$Id))

## NUMBER OF INDIVIDUALS DETECTED ALIVE AND RECOVERED DEAD
sum(unique(data.alive$Id) %in% unique(data.dead$Id))

## NUMBER OF INDIVIDUALS DETECTED/YEAR/COUNTRY
table.id <- table(data.alive$Id,data.alive$Year, data.alive$Country)
apply(table.id, c(2,3), function(x) sum(x>0))

## NUMBER OF DETECTIONS/YEAR/COUNTRY
countrytab <- table(data.alive$Year, data.alive$Country)
countrytab 

## NUMBER OF DETECTIONS/YEAR/COUNTRY/SEX
sex_countrytab <- table(data.alive$Year, data.alive$Country, data.alive$Sex)
sex_countrytab

## PROPORTION OF DETECTIONS PER YEAR/COUNTRY/SEX
sex_countrytab_prop <- sex_countrytab
for(i in 1:dim(sex_countrytab)[3]){
  sex_countrytab_prop[,,i] <- sex_countrytab_prop[,,i]/countrytab
}
sex_countrytab_prop



## ------       4.1.2. DEAD RECOVERIES ------ 

## NUMBER OF INDIVIDUALS RECOVERED
length(unique(data.dead$Id))

## NUMBER OF DEAD RECOVERIES/YEAR/COUNTRY
table(data.dead$Year, data.dead$Country)

## MORTALITY CAUSES
unique(as.character(data.dead$DeathCause))
unique(as.character(data.dead$DeathCause_2))  
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))

## DEFINE LEGAL MORTALITY
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])



## ------       4.1.3. PLOTS ------ 

legal.death <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]
Other.death <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$DeathCause %in% legalCauses, ]

table(legal.death$Year)
table(Other.death$Year)

## PLOT CHECK
if(plot.check){
  pdf(file = file.path(working.dir,"figures","DetectionsDeadRecoveries.pdf"))
  
  ## PLOT ALL DETECTIONS
  for(t in 1:nYears){
    plot(myHabitat.list$habitat.r, main = years[t]) 
    plot(st_geometry(mymain.detector.sp), add=T, pch=16, cex=0.1)
    plot(st_geometry(data.alive[data.alive$Year == years[t], ]), 
         pch=16,col="red", cex=0.7, add=T)
    mtext(paste("Live detections",
                nrow(data.alive[data.alive$Year == years[t], ]),
                "; ID:",
                nrow(unique(data.alive[data.alive$Year == years[t], ]$Id))),
          line = +1)
    
    plot(st_geometry(data.dead[data.dead$Year == years[t], ]),
         pch=16,col="blue", cex=0.7,add=T)
    mtext(paste("Dead recovery:",nrow(data.dead[data.dead$Year == years[t], ])))
    plot(st_geometry(GLOBALMAP), add=T) 
  }#t
  
  ## PLOT NUMBER OF MORTALITY EVENTS, LEGAL/OTHERS 
  barplot(table(legal.death$Country, legal.death$Year), main = "Legal causes")
  barplot(table(Other.death$Country, Other.death$Year), main = "Other causes")
  legend("topleft", fill = c(grey(0.3), grey(0.7)), legend = c("NOR","SWE"))
  dev.off()  
}

## EXPORT NGS DATA 
if(DATA$sex == "Hann"){
  assign("myFilteredData.spM", myFilteredData.sp)
  assign("myFilteredData.spOthersM", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
  
  save(myFilteredData.spM, 
       myFullData.spM,
       myFilteredData.spOthersM,
       myFilteredData.spStructuredM,
       file = file.path(working.dir, "data", "NGSData.RData"))
} else {
  assign("myFilteredData.spF", myFilteredData.sp)
  assign("myFilteredData.spOthersF", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
  
  save(myFilteredData.spF,
       myFullData.spF,
       myFilteredData.spOthersF,
       myFilteredData.spStructuredF,
       file = file.path(working.dir, "data", "NGSData.RData"))
}



## ------     4.2. ASSIGN DETECTORS ------ 

## ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = data.alive,                
  myDetectors = mymain.detector.sp,
  mysubDetectors = mydetector.sp,
  radius = detector.res)

## STRUCTURED
myData.aliveStruc <- AssignDetectors_v3sf( 
  myData = myFilteredData.spStructured,                
  myDetectors = mymain.detector.sp,
  mysubDetectors = mydetector.sp,
  radius = detector.res)

## OTHERS
myData.aliveOthers <- AssignDetectors_v3sf( 
  myData = myFilteredData.spOthers,                
  myDetectors = mymain.detector.sp,
  mysubDetectors = mydetector.sp,
  radius = detector.res)

## DEAD RECOVERIES
myData.dead <- AssignDetectors_v3sf(
  myData = data.dead,
  myDetectors = mymain.detector.sp,
  radius = detector.res)



## ------     4.3. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ------ 

## ALL SAMPLES
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 myDetectors = mymain.detector.sp,
                 method = "Binomial",
                 myData2 = myData.dead,
                 myDetectors2 = mymain.detector.sp,
                 returnIdvector = TRUE)
y.ar.ALIVE <- y.ar$y.ar
dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)

## STRUCTURED
y.arStruc <- MakeYsf( myData = myData.aliveStruc$myData.sp,
                      myDetectors = mymain.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      myDetectors2 = mymain.detector.sp,
                      returnIdvector = TRUE)
y.ar.ALIVEStruc <- y.arStruc$y.ar
dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)

## OTHERS
y.arOth <- MakeYsf( myData = myData.aliveOthers$myData.sp,
                    myDetectors = mymain.detector.sp,
                    method = "Binomial",
                    myData2 = myData.dead,
                    myDetectors2 = mymain.detector.sp,
                    returnIdvector = TRUE)
y.ar.ALIVEOth <- y.arOth$y.ar
dimnames(y.ar.ALIVEOth) <- dimnames(y.arOth$y.ar)

## MAKE SURE THE Y HAVE THE SAME DIMENSIONS#
y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0

## FILL IN THE Y ARRAYS 
y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]],,] <-  y.ar.ALIVEOth
y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]],,] <-  y.ar.ALIVEStruc

## PROJECT THE DEATH TO THE NEXT OCCASION.
y.ar.DEADProjected <- y.ar$y.ar2 
y.ar.DEADProjected[] <- 0
for(t in 2:nYears){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}

y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])



## ------     4.4. GENERATE OBSERVATIONS : y.obs[i,t] ------ 

INDIVIDUAL_ID$STATUS_Numeric <- 1
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% "Juvenile"] <- 2
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Pair" )] <- 3
INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Family group")] <- 4

## FOR THE LAST YEAR GET THROUGH THE PAIR BASED-FILE FROM LINN.
indIDRovBase <- unlist(lapply(strsplit(y.ar$Id.vector," "), function(x) x[1]))

ALLIDS <- c(unique(indIDRovBase))

y.obsALL <- matrix(1, nrow = length(ALLIDS), ncol = dim(y.ar.ALIVE)[3]+1)
yrs <- c(years[1]-1, years)
dimnames(y.obsALL) <- list(ALLIDS, yrs)
for(i in 1:dim(y.obsALL)[1]){
  for(t in 1:(nYears+1)){
    tmp <- unique(INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$ReprodYear..May.1.year.y...Apr.30.y.1. == yrs[t] & 
                                                 INDIVIDUAL_ID$ROVBASE_IndividID == ALLIDS[i]])
    if(length(tmp) > 0){
      if(length(tmp) > 1){
        print(tmp)
        tmp <- tmp[1]
      }
      y.obsALL[i,t] <- tmp
    }
  }#t
}#i

## FOR THE LAST YEAR GET TRHOUGH THE PAIR BASED-FILE FROM LINN.
for(i in 1:dim(y.obsALL)[1]){
  # for(t in 1:(nYears+1)){
  t <- nYears
  ## linn's 2023 file
  tmp <-  Pack_ID2023[Pack_ID2023$Rovbase.ID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    y.obsALL[i,t-1] <- 3
  }
  ## linn's 2024 file
  tmp <-  Pack_ID2024[Pack_ID2024$RovbaseID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    y.obsALL[i,t] <- 3
  }
  ## linn's 2025 file
  tmp <-  Pack_ID2025[Pack_ID2025$IndividID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    y.obsALL[i,t+1] <- 3
  }
}#i

## subset y.obs for the individuals present in y.ar
indID <- unlist(lapply(strsplit(y.ar$Id.vector, " "), function(x)x[1]))
y.obs <- y.obsALL[indID, ]
#y.obs.family <- y.obs
y.obs[y.obs == 4] <- 3
y.obs <- y.obs[ ,as.character(years)]

## PLOT CHECK 
if(plot.check){
  pdf(file = file.path(working.dir, "figures", "States.pdf"))
  tabInd <- table(INDIVIDUAL_ID$ReprodYear..May.1.year.y...Apr.30.y.1., INDIVIDUAL_ID$Status)
  SumID <- rowSums(tabInd)
  yearsind <- c(as.numeric(row.names(tabInd)))
  plot(-10, xlim = range(yearsind), ylim=c(0,500), ylab="N individuals")
  polygon(c(yearsind,rev(yearsind)), c(SumID, rep(0,length(SumID))), col="black")
  polygon(c(yearsind,rev(yearsind)), c(SumID-tabInd[,1], rev(SumID)), col="red")
  polygon(c(yearsind,rev(yearsind)), c(SumID-tabInd[,1]-tabInd[,2], rev(SumID-tabInd[,1])), col="blue")
  polygon(c(yearsind,rev(yearsind)), c(SumID-tabInd[,1]-tabInd[,2]-tabInd[,3], rev(SumID-tabInd[,2]-tabInd[,1])), col="green")
  legend("topleft", fill=c("black", "green", "blue", "red"), legend = c("Pair", "Juvenile","Family group", "Unknown" ) )
  dev.off()
}



## ------     4.5. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------ 

distances <- list()
for(t in 1:nYears){
  print(paste("------ ", t ," -------", sep = "" ))
  distances[[t]] <- CheckDistanceDetectionsV2sf(
    y = y.ar.ALIVE[,,t], 
    detector.xy = st_coordinates(mymain.detector.sp), 
    max.distance = maxDist,
    method = "pairwise",
    plot.check = F)
  
  ## PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
  if(plot.check){
    par(mfrow=c(1,1),mar=c(1,1,1,1))
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      for(i in affected.ids){
        plot(st_geometry(studyArea), main = paste("t: ",t,"     i: ", i, sep = ""))
        plot(st_geometry(GLOBALMAP), add = T)
        plot(st_geometry(mymain.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        tmp <- data.alive[data.alive$Id == y.ar$Id.vector[i], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1)
        
        arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2], length = 0.1, lwd = 1)
        
        plot(st_geometry(mymain.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
        
        tmp2 <- mymain.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
        plot(st_geometry(tmp2), col = "blue", pch = 13, cex = 1.5, lwd = 1,add=T)
      }#i
    }#if
  }#if plot.check
  
  ## REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
  y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]),]
  }#i
}#t



## ------   5. GENERATE INDIVIDUAL-LEVEL COVARIATES ------ 

## ------     5.1. INDIVIDUAL STATE ------ 

indSocialState <- matrix(1, nrow = dim(y.ar.ALIVE)[1], ncol = dim(y.ar.ALIVE)[3])
for(i in 1:dim(indSocialState)[1]){
  if(any(y.obs[i, ] >= 3)){
    indSocialState[i, min(which(y.obs[i, ] >= 3)):dim(y.ar.ALIVE)[3]] <- 2
  }
}#i



## ------     5.2. TRAP-RESPONSE ------ 

## Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCovsf(myFullData.sp$alive, myFullData.sp$dead.recovery)

## Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]

## Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]

## Plot an image of the matrix
if(plot.check){image(t(already.detected))}



## ------   6. GENERATE HABITAT-LEVEL COVARIATES ------ 

## KERNEL OF INDIVIDUALS IN PAIRS 
kern <- list()
habDens <- matrix(NA, nrow = nHabCells, ncol = nYears)
IDS <- unlist(lapply(strsplit(as.character(myFullData.sp$alive$Id) , " "), function(x)x[1])) 
for(t in 1:nYears){
  id.fam <- which(y.obsALL[,as.character(years[t]-1)]%in% c(3,4), arr.ind = T)
  m.xy <- matrix(NA, nrow=length(id.fam), ncol=2 )
  colnames(m.xy) <-c("x","y")
  for(i in 1:length(id.fam)){
    tmp <- myFullData.sp$alive[IDS==row.names(y.obsALL)[i],]
    #plot(st_geometry(myHabitat.list$buffered.habitat.poly))
    m.xy[i,]<- colMeans(st_coordinates(tmp))
    #points(m.xy[i,2]~m.xy[i,1], col="red", pch=16)
  }
  if(sum(is.na(m.xy[,1]))>0){
    m.xy <- m.xy[!is.na(m.xy[,1]),]
  }
  locationsFamily <- st_as_sf(as.data.frame(m.xy), coords =c("x","y"),crs=st_crs(myHabitat.list$habitat.sp))
  locationsFamily$id <- rep(1, nrow(locationsFamily))
  kern[[t]] <- raster(estUDm2spixdf(kernelUD(as(locationsFamily[ ,"id"],"Spatial"),h = 15000,
                                             grid = as(myHabitat.list$habitat.r, 'SpatialPixels'))))
  #plot(kern[[1]])
  plot(st_geometry(myHabitat.list$habitat.poly), add=T)
  habDens[,t] <- scale(kern[[t]][myHabitat.list$habitat.r[ ]==1])
}

#check 
for(t in 1:nYears){
  plot(kern[[t]],main=years[t])
  plot(myHabitat.list$habitat.poly$geometry,add=T,col=NA)
}



## ------   7. MAKE AUGMENTATION ------ 

## DATA ARRAYS
y.alive <- MakeAugmentation(y = y.ar.ALIVE, aug.factor = aug.factor, replace.value = 0)
y.aliveOthers <- MakeAugmentation(y = y.ar.ALIVEOthers, aug.factor = aug.factor, replace.value = 0)
y.aliveStructured <- MakeAugmentation(y = y.ar.ALIVEStructured, aug.factor = aug.factor, replace.value = 0)

y.dead <- MakeAugmentation(y = y.ar.DEAD, aug.factor = aug.factor, replace.value = 0)
y.obs <- MakeAugmentation(y = y.obs, aug.factor = aug.factor, replace.value = 1)

## INDIVIDUAL COVARIATES
indSocialState <- MakeAugmentation(y = indSocialState, aug.factor = aug.factor, replace.value = 1)
already.detected <- MakeAugmentation(y = already.detected, aug.factor = aug.factor, replace.value = 0)
age <- MakeAugmentation(y = age, aug.factor = aug.factor, replace.value = NA)
min.age <- MakeAugmentation(y = min.age, aug.factor = aug.factor, replace.value = NA)
precapture <- MakeAugmentation(y = precapture, aug.factor = aug.factor, replace.value = 0)



## -----------------------------------------------------------------------------

## ------ III. MODEL SETTING & RUNNING ------- 

## ------   1. NIMBLE MODEL DEFINITION ------ 

modelCode <- nimbleCode({
  
  ##------ SPATIAL PROCESS ------## 
  
  for(st in 1:2){
    dmean[st] ~ dunif(0,100)
    lambda[st] <- 1/dmean[st]
  }#st
  
  beta.dens ~ dnorm(0.0,0.01)
  
  for(t in 1:n.years){
    habIntensity[1:numHabWindows,t] <- exp(beta.dens * habDens[1:numHabWindows,t])
    sumHabIntensity[t] <- sum(habIntensity[1:numHabWindows,t])
    logHabIntensity[1:numHabWindows,t] <- log(habIntensity[1:numHabWindows,t])
    logSumHabIntensity[t] <- log(sumHabIntensity[t])
  }#t
  
  for(i in 1:n.individuals){
    
    sxy[i,1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
      upperCoords = upperHabCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows,1],
      logSumIntensity = logSumHabIntensity[1],
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    for(t in 2:n.years){
      
      sxy[i,1:2,t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
        upperCoords = upperHabCoords[1:numHabWindows,1:2],
        s = sxy[i,1:2,t-1],
        lambda = lambda[state[i,t-1]+1],
        baseIntensities = habIntensity[1:numHabWindows,t],
        habitatGrid =  habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max,
        numWindows = numHabWindows)
    }#i
  }#t
  
  
  ##----- DEMOGRAPHIC PROCESS -----##
  
  ## FIRST YEAR
  omeg1[1:3] ~ ddirch(alpha[1:3])  
  
  for(i in 1:n.individuals){
    z[i,1] ~ dcat(omeg1[1:3])
    isAlive[i,1] <- (z[i,1] == 2) + (z[i,1] == 3)
    state[i,1] <- (z[i,1] == 3)
  }#i
  
  ## FOLLOWING YEARS 
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
    idResponse[i,1] ~ dbern(pResponse)
  }
  
  for(t in 1:n.years){
    for(st in 1:2){
      sigma[st,t] ~ dunif(0,50)
    }#st
    
    ## Structured
    betaResponse[t] ~ dunif(-5,5)
    for(n in 1:nTrapCovs){
      trapBetas[n,t] ~ dunif(-5,5)
    }#n
    for(c in 1:n.counties){
      p0[c,1,t] ~ dunif(0,1)
      p0[c,2,t] ~ dunif(0,1)
    }#c    
    
    ## Opportunistic
    betaResponseOth[t] ~ dunif(-5,5)
    for(n in 1:nTrapCovsOth){
      trapBetasOth[n,t] ~ dunif(-5,5)
    }#n
    for(c in 1:n.countries){
      p0Oth[c,1,t] ~ dunif(0,1)
      p0Oth[c,2,t] ~ dunif(0,1)
    }#c
    
    for(i in 1:n.individuals){
      
      ## Structured
      y.alive[i,1:nMaxDetectors,t] ~ dbinomLocal_normalWolf(
        detNums = nbDetections[i,t],
        detIndices = yDets[i,1:nMaxDetectors,t],
        size = trials[1:n.detectors],
        p0 = p0[1:n.counties,1:2,t],
        sigma = sigma[state[i,t]+1,t],
        s = sxy[i,1:2,t],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = ResizeFactor,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        z = z[i,t]-1,
        trapCovsIntercept = detCounties[1:n.detectors],
        indCov = idResponse[i,t],
        indBeta = betaResponse[t],
        trapCovs = trapCovs[1:n.detectors,1:nTrapCovs,t],
        trapBetas = trapBetas[1:nTrapCovs,t],
        lengthYCombined = 1)
      
      ## Opportunistic
      y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbinomLocal_normalWolf(
        detNums = nbDetectionsOth[i,t],
        detIndices = yDetsOth[i,1:nMaxDetectorsOth,t],
        size = trials[1:n.detectors],
        p0 = p0Oth[1:n.countries,1:2,t],
        sigma = sigma[state[i,t]+1,t],
        s = sxy[i,1:2,t],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = ResizeFactor,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        z = z[i,t]-1,
        trapCovsIntercept = detCountries[1:n.detectors],
        indCov = idResponse[i,t],
        indBeta = betaResponseOth[t],
        trapCovs = trapCovsOth[1:n.detectors,1:nTrapCovsOth,t],
        trapBetas = trapBetasOth[1:nTrapCovsOth,t],
        lengthYCombined = 1)
      
      ## Dead recoveries legal
      x.deadculled[i,t] ~ dbern(z[i,t] == 4)
      
      ## Dead recoveries others
      x.deadOther[i,t] ~ dbern(z[i,t] == 5)
    }#i
  }#t
  
  ##---------- DERIVED PARAMETERS ----------##
  
  for(t in 1:n.years){
    N[t] <- sum(isAlive[1:n.individuals,t])
  }#t
  
})



## ------   2. Z ------ 

## CREATE Z 
z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, 0)

z.dead <- apply(y.dead, c(1,2), function(x) any(x>0))
z.dead <- ifelse(z.dead, 3, 0)
z <- ifelse(z.dead+z==0 ,NA, z.dead+z)

## IDENTIFY INDIVIDUALS DEAD TO CULLING
aug.id.list <- c(y.ar$Id.vector, rep("AugInd", dim(z)[1]-length(y.ar$Id.vector)))

legal.mx <- do.call(rbind, lapply(aug.id.list, function(x){
  out <- rep(0,dim(z)[2])
  if(x %in% legal.death$Id) out <- rep(1,dim(z)[2])
  return(out)
}))

other.mx <- legal.mx
other.mx <- ifelse(legal.mx==0 & z.dead %in% c(3), 1 ,0)

z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    if(sum(zz ==3,na.rm = T)<1){
      zz[range.det[1]:range.det[2]] <- 2
    }
    if(sum(zz ==3,na.rm = T)>0){
      if(sum(zz,na.rm=T)>3){
        zz[range.det[1]: (range.det[2]-1)] <- 2
      }
      t.recovered <- which(zz==3)
      if(t.recovered< dim(z)[2]){
        zz[(t.recovered+1):(dim(z)[2])] <- 4
      }
      if(t.recovered <= dim(z)[2]){
        zz[(t.recovered-1)] <- 2
        print(1)
      }
    }
  }
  return(zz)
}))


# #if dead last year, make sure it was alive the year before
# if(out[length(out)]==5 & out[length(out)-1]==1 ) {
#   out[length(out)-1] <- 2
# }


z[z %in% 4] <- 6

## id culled get the state 4. id not culled get the 5. 
z <- ifelse(z==3 & other.mx %in% c(1), 5 ,z) # remove other causes of mortality
z <- ifelse(z==3 & legal.mx %in% c(1), 4 ,z)

## IDENTIFY SOCIAL STATE
z <- ifelse(z==2 & indSocialState %in% c(1), 2 ,z)
z <- ifelse(z==2 & indSocialState %in% c(2), 3 ,z)

## CREATE INITIAL Z VALUES  
z.init <- t(apply(z, 1, function(zz){
  out <- zz
  out[] <- 1
  if(any(!is.na(zz))){
    if(sum(zz == 4, na.rm = T)<1){
      range.det <- range(which(!is.na(zz)))
      if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
      if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 6
    }
    
    if(sum(zz == 3, na.rm = T)>0){
      reco.3 <- min(which(zz == 3))
      if(reco.3>1){zz[(reco.3-1)] <- 2}
    }
    if(sum(zz == 2, na.rm = T)>0){
      reco.alive <- min(which(zz == 2))
      if(reco.alive>1){zz[1:(reco.alive-1)] <- 1}
    }
    out[] <- zz
    
    ## if still some NAs initialize it with 2
    if(sum(is.na(out)>0)) {
      out[is.na(out)] <- 2
    }
  }
  return(out)
}))

z.init[!is.na(z)] <- NA

z.age <- z

x.deadculled <- x.deadOther <- z.age
x.deadculled[] <- ifelse(z.age %in% c(4) & legal.mx == 1, 1, 0)
x.deadculled <- t(apply(x.deadculled, 1, function(x){
  out <- x
  out[] <- 0
  if(any(x==1)) out[min(which(x==1))] <- 1
  return(out)
}))

x.deadOther[] <- ifelse(z.age %in% c(5) & other.mx == 1, 1 ,0)
x.deadOther <- t(apply(x.deadOther, 1, function(x){
  out <- x
  out[] <- 0
  if(any(x==1)) out[min(which(x==1))] <- 1
  return(out)
}))



## ------   4. NIMBLE DATA  ------ 

## ------     4.1. TRAP COVARIATES ------ 

#STRUCTURED 
trapCovs <- array(0,c(nrow(detTracks), 2, nYears))
trapCovs[,1,] <- detTracks
trapCovs[,2,] <- detSnow

#OTHERS 
trapCovsOth <- array(0,c(nrow(detTracks), 3, nYears))
trapCovsOth[,1,] <- detRoads
trapCovsOth[,2,] <- detSnow
trapCovsOth[,3,] <- detOtherSamples




## ------     4.2. GENERATE sxy & sxy.init ARAYS  ------ 

## SXY INITS 
#create a data.frame with all detections of all Individuals detected
#project death to the next year
myData.deadProj <- myData.dead[ ,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1#project dead reco to the next year
#remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year ),]
myData.deadProj$Id

AllDets <- rbind(myData.alive$myData.sp[,c("Id","Year")],
                 myData.deadProj[,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  coordsData = AllDetsxy,
  coordsHabitatGridCenter = myHabitat.list$habitat.xy,
  scaleToGrid =T )$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)

idAugmented <- which(rownames(z) %in%"Augmented")

sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = y.ar$Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = ScaledLowUpCoords$lowerHabCoords,
                       upperCoords = ScaledLowUpCoords$upperHabCoords,
                       habitatGrid = ScaledLowUpCoords$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")

## SXY DATA 
sxy.data <- sxy.init
sxy.data[] <- NA
i=783
t=9
for(i in 1:length(y.ar$Id.vector)){
  for(t in 1:(dim(sxy.data)[3]-1)){
    if(sum(z.age[i,t+1] %in% c(4,5,6)) > 0){
      # temp <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Id == y.ar$Id.vector[i] &
      #                                     myFullData.sp$dead.recovery$Year == years[t], ]
      print(i)
      print(t)
      # if(length(temp) > 0){
      #if(!is.na(over(temp,myStudyArea.extent))){
      #if(raster::extract(myHabitat.list$habitat.r, temp)==0){
      
      #id could be dead in the spatial extent but can be outside the habitat because in a cell <49%habitat
      # buff <- gBuffer(temp,width = myHabitat.list$resolution*2)
      # inter <- intersect(buff,habitat.poly)
      sxy.data[i, ,t+1]  <- sxy.init[i, ,t+1] #coordinates(spsample(x = inter,n = 1,type="random"))
      sxy.init[i, ,t+1]  <- NA#coordi
      #}else{sxy.data[i, ,t+1] <- coordinates(temp)}
      #}
      # }
    }
  }
}


# plot check 
# plot(ScaledLowUpCoords$lowerHabCoords[,2]~ScaledLowUpCoords$lowerHabCoords[,1])
# points(sxy.init[i, 2,t+1]~sxy.init[i, 1,t+1],col="red",pch=16)
# plot(myHabitat.list$habitat.r)
# plot(st_geometry(AllDets[AllDets$Id == y.ar$Id.vector[i] &
#                            AllDets$Year == years[t+1], ]),col="red",add=T)
# #
# tmp <- scaleCoordsToHabitatGrid(coordsData = sxy.init,
#                          coordsHabitatGridCenter = habitatxy,scaleToGrid = F)$coordsDataScaled
# tmp[i, ,t+1]
# st_coordinates(myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Id == y.ar$Id.vector[i] &
#                                           myFullData.sp$dead.recovery$Year == years[t], ])

sxy.init <- round(sxy.init, 5)#---an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
sxy.data <- round(sxy.data, 5)
# sxy.init[!is.na(sxy.data)] <- NA




## ------     4.3. TRANSFORM Y TO SPARSE MATRICES  ------ 
#STRUCTURED 
SparseY <- getSparseY(y.aliveStructured)
#OTHER
SparseYOth <- getSparseY(y.aliveOthers)

## ------     4.4. LATENT VARIABLE DET RESPONSE ------ 
detResponse <- already.detected 
detResponse[rownames(detResponse) %in% "Augmented", 1]  <- NA
InitsDetResponse <- detResponse
InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1,0.5)
InitsDetResponse[!is.na(detResponse)] <- NA




## ------   5. NIMBLE DATA ------ 
nimData <- list( z = z.age,   
                 sxy = sxy.data,
                 y.alive = SparseY$y,
                 yDets = SparseY$detIndices,
                 nbDetections = SparseY$detNums,
                 y.aliveOth = SparseYOth$y, 
                 yDetsOth = SparseYOth$detIndices,
                 nbDetectionsOth = SparseYOth$detNums,
                 x.deadculled = x.deadculled,
                 x.deadOther = x.deadOther,
                 
                 ones.dead.legal = array(1,c(2,dim(y.alive)[3]-1)),
                 idResponse = detResponse,
                 alpha = rep(1,3))

## ------   6. NIMBLE CONSTANTS ------ 
nimConstants <- list( n.individuals = dim(y.alive)[1],
                      n.detectors = dim(y.alive)[2],  
                      n.years = dim(y.alive)[3], 
                      n.years1 = dim(y.alive)[3]-1, 
                      numHabWindows = nHabCells,
                      n.counties = max(detCounties),
                      n.countries = max(detCountries),
                      nTrapCovs = 2,
                      nTrapCovsOth = 3,
                      detCountries = detCountries,
                      detCounties = detCounties,
                      y.max = dim(ScaledLowUpCoords$habitatGrid)[1],
                      x.max = dim(ScaledLowUpCoords$habitatGrid)[2],
                      habitatGrid = ScaledLowUpCoords$habitatGrid,
                      lowerHabCoords = ScaledLowUpCoords$lowerHabCoords,
                      upperHabCoords = ScaledLowUpCoords$upperHabCoords,
                      detector.xy = as.matrix(ScaledcoordsDataScaled),
                      lowerHabCoords = as.matrix(ScaledLowUpCoords$lowerHabCoords),
                      upperHabCoords = as.matrix(ScaledLowUpCoords$upperHabCoords),
                      y.maxDet = dim(DetectorIndexLESS$habitatGrid)[1],
                      x.maxDet = dim(DetectorIndexLESS$habitatGrid)[2],
                      ResizeFactor = DetectorIndexLESS$resizeFactor,
                      n.cellsSparse = dim(DetectorIndexLESS$localIndices)[1],
                      maxNBDets = DetectorIndexLESS$numLocalIndicesMax,
                      detectorIndex = DetectorIndexLESS$localIndices,
                      nDetectorsLESS = DetectorIndexLESS$numLocalIndices,
                      habitatIDDet = DetectorIndexLESS$habitatGrid,
                      nMaxDetectors = SparseY$maxDetNums,
                      nMaxDetectorsOth = SparseYOth$maxDetNums,
                      habDens = as.matrix(habDens),
                      trapCovsOth = trapCovsOth,
                      trapCovs = trapCovs,
                      trials = n.trials
)

## ------   7. NIMBLE PARAMETERS ------ 
nimParams <- c("N",
               "omeg1",
               "gamma",
               "p0",
               "phi",
               "h",
               "w",
               "wAll",
               "rw",
               "psi",
               "lambda",
               "p0",
               "p0Oth",
               "sigma",
               "betaResponse",
               "trapBetas",
               "trapBetasOth",
               "betaResponseOth",
               "pResponse",
               "beta.dens")

nimParams2 <- c("z",
                "sxy")

## ------   8. SAVE NECESSARY OBJECTS FOR PLOTTING ------ 
#ONLY SAVE IF IT IS THE FEMALE SCRIPT TO AVOID DUPLICATED SCRIPTS
if(DATA$sex %in% "Hunn"){
  if(!dir.exists(file.path(WD,"/Figures",modelName))){dir.create(file.path(WD,"/Figures", modelName))}
  
  save(myHabitat.list, myDetectors, COUNTRIES,
       studyArea,COMMUNES,habitat.subdetectors,
       myFilteredData.sp, myFullData.sp, COUNTIESplot,
       detCounties.original,
       file = file.path(paste(WD,"/Figures/",modelName,sep=""), "NecessaryObjects.RData" ))
}

## ------   9. SET UP SEVERAL CHAINS WITH DIFFERENT STARTING VALUES ------ 

for(c in 1:4){#----SET UP SEVERAL CHAINS WITH DIFFERENT STARTING VALUES
  
  ## ------    3.4. LIST NIMBLE INITS ------ 
  nimInits <- list( "sxy" = sxy.init,
                    "dmean" = runif(2,2,4),
                    "z" = z.init,
                    "omeg1" = c(0.5,0.25,0.25),
                    "psi" = runif(dim(y.alive)[3]-1,0.1,0.7),
                    "gamma" = runif(dim(y.alive)[3]-1,0,1),
                    "p0" = array(runif(12,0,0.2), c(nimConstants$n.counties,2,dim(y.alive)[3])),
                    "p0Oth" = array(runif(12,0,0.2), c(nimConstants$n.countries,2,dim(y.alive)[3])),
                    "betaResponse" = runif(dim(y.alive)[3],-1,1),
                    "betaResponseOth" = runif(dim(y.alive)[3],-1,1),
                    "beta.dens" = runif(1,-1,1),
                    "trapBetas" = array(runif(nimConstants$nTrapCovs,-1,1),c(nimConstants$nTrapCovs,dim(y.alive)[3])),
                    "trapBetasOth" = array(runif(nimConstants$nTrapCovs,-1,1),c(nimConstants$nTrapCovsOth,dim(y.alive)[3])),
                    "sigma" = array(runif(2,4,8),c(2,dim(y.alive)[3])),
                    "idResponse" = InitsDetResponse,
                    "pResponse"  = runif(1, 0, 1),#[CM]#0,
                    #--[RB]:
                    "h" =  array(runif((dim(y.alive)[3]-1)*2,0.2,0.4), c(2,dim(y.alive)[3]-1)),
                    "rw" =  array(runif((dim(y.alive)[3]-1)*2,0.05,0.10), c(2,dim(y.alive)[3]-1)),
                    "w" = array(runif((dim(y.alive)[3]-1)*2,0.2,0.4), c(2,dim(y.alive)[3]-1)))
  
  
  
  
  
  
  ### TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  i=1
  t=4
  for(i in 1:nimConstants$n.individuals){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i,,t]  
      }else{SXY <- nimData$sxy[i,,t]}
      sxyID <- nimConstants$habitatIDDet[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      index <- nimConstants$detectorIndex[sxyID, 1:nimConstants$nDetectorsLESS[sxyID]]
      ## GET NECESSARY INFO 
      n.detectors <- length(index)
      #maxDist_squared <- maxDist*maxDist
      
      YDET <- nimData$yDets[i,, t]
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
  
  ## ------ 6. SAVE NIMBLE INPUT ------ 
  save(nimData,
       nimConstants,
       nimParams,
       nimParams2,
       
       modelCode,
       nimInits,
       file = paste(WD,"/",modelName,"/",modelName,"_INPUTChain",c,".RData", sep="" ))
}

## ------   7. NIMBLE RUN ------ 
load(file.path(WD, modelName, paste(modelName, "_INPUTChain1.RData", sep="" )))
ptm <- proc.time()
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()
cmodel <- compileNimble(model)
cmodel$calculate() 
which(is.infinite(model$logProb_z),arr.ind = T)



## -----------------------------------------------------------------------------

## ------ IV. MAKE THE SINGLE SEASON SCR MODEL ------- 

## ------   1. NIMBLE MODEL DEFINITION ------ 

modelCode1 <- nimbleCode({
  
  ##------ SPATIAL PROCESS ------##  
  beta.dens  ~ dnorm(0.0,0.01)
  habIntensity[1:numHabWindows] <- exp(beta.dens * habDens[1:numHabWindows])
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
  
  for(t in 1:nYears){   
    
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



## -----------------------------------------------------------------------------

## ------ V. PROCESS RESULTS ------

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
  
  plot(10, xlim = c(0, nYears+1), ylim = c(80,350), type ="n", xaxt="n", xlab = "Years", ylab = "N")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
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
  
  plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "h")
  axis(1, c(1:nYears),labels = years)
  myCol <- c("firebrick3","navyblue")
  myDev <- c(-0.2,+0.2)
  for(t in 1:(nYears-1)){
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
  
  plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "wAll")
  axis(1, c(1:nYears),labels = years)
  myDev <- c(-0.2,+0.2)
  for(t in 1:(nYears-1)){
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
  
  plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "rw")
  axis(1, c(1:nYears),labels = years)
  myDev <- c(-0.2,+0.2)
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "phi")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  for(t in 1:(nYears-1)){
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
    plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
    for(s in 1:2){
      for(t in 1:nYears){
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
    plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=main[c])#paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
    for(s in 1:2){
      for(t in 1:nYears){
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
  plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "psi")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(0,10000), type ="n", xaxt="n", xlab = "Years", ylab = "sigma")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  abline(v=seq(1.5,nYears-0.5,by=1),lty=2)
  for(s in 1:2){
    for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$sigma[ ,s,t]*myHabitat.list$resolution),
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
  #     plot.violins(list(myResults$sims.list$sigma[ , s]*myHabitat.list$resolution),
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
    plot.violins(list(myResults$sims.list$lambda[ , s]*myHabitat.list$resolution),
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
  plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta.dens")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-2,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta Tracks")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta Snow")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta other road")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta other snow")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "beta other location other samples ")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponse")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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
  plot(-10, xlim = c(0,nYears), ylim=c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponseOth")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  for(t in 1:(nYears-1)){
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

for(t in 1:nYears){
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
    nimOutput[[i]][,"sigma[2]"] <-  nimOutput[[i]][,"sigma[2]"]*myHabitat.list$resolution
    nimOutput[[i]][,"sigma[1]"] <-  nimOutput[[i]][,"sigma[1]"]*myHabitat.list$resolution
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

myResultsListALL$sims.list$trapBetas <- array(NA,c(dim(myResultsListALL$sims.list$trapBetas),nYears) )
myResultsListALL$sims.list$trapBetasOth <- array(NA,c(dim(myResultsListALL$sims.list$trapBetasOth),nYears) )

myResultsListALL$sims.list$p0 <- array(NA,c(dim(myResultsListALL$sims.list$p0),nYears) )
myResultsListALL$sims.list$p0Oth <- array(NA,c(dim(myResultsListALL$sims.list$p0Oth),nYears) )

myResultsListALL$sims.list$sigma <- array(NA,c(dim(myResultsListALL$sims.list$sigma),nYears) )

for(t in 1:nYears){
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
  
  plot(10, xlim = c(0, nYears+1), ylim = c(80,350), type ="n", xaxt="n", xlab = "Years", ylab = "N")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$N[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "N")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.4.p0 ------   
  
  par(mfrow=c(1,2))
  myDev <- c(-0.2,+0.2)
  myCol <- c("firebrick3","navyblue")
  par(mfrow=c(1,2))
  for(c in 1:6){
    plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
    for(s in 1:2){
      for(t in 1:nYears){
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
  
  for(t in 1:nYears){
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
    plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.2), type ="n", xaxt="n", xlab = "Years",
         ylab = "p0", main=main[c])#paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
    axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
    for(s in 1:2){
      for(t in 1:nYears){
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
  for(t in 1:nYears){
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutputList[[t]], params = params[i])
      title(years[t], line = 0.5)
    }
  }
  
  
  
  ## ------   2.5.sigma ------   
  
  plot(10, xlim = c(0, nYears+1), ylim = c(5000,12000), type ="n", xaxt="n", xlab = "Years", ylab = "sigma")
  axis(1, c(1:nYears),labels = years)
  offset <- c(-0.25,0.25)
  for(t in 1:nYears){
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
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "sigma[1]")
    title(years[t],line = 0.5)
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "sigma[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.5.roads ------  
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-3,3), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Tracks")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$trapBetas[,1,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetas[1]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6.tracks  ------ 
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$trapBetas[,2,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetas[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------   
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta road (other)")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,1,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[1]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------ 
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow (other)")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,2,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[2]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. Snow ------ 
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-2,2), type ="n", xaxt="n", xlab = "Years", ylab = "Beta Snow (opportunistic)")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$trapBetasOth[,3,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "trapBetasOth[3]")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6. pResponse ------   
  
  plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "pResponse")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$pResponse[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "pResponse")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.6. betaResponse ------   
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-1,2), type ="n", xaxt="n", xlab = "Years", ylab = "betaResponse")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$betaResponse[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "betaResponse")
    title(years[t],line = 0.5)
  }
  
  
  
  ## ------   2.7. beta.dens ------   
  
  plot(10, xlim = c(0, nYears+1), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "beta.dens")
  axis(1, c(1:nYears),labels = years)
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$beta.dens[,t]),
                 x = t,
                 at = t,
                 violin.width = 0.3,
                 col = "firebrick3",
                 add = T,
                 alpha = 0.2,
                 border.col = "firebrick3")
  }#t
  
  for(t in 1:nYears){
    PlotJagsParams(jags.samples = nimOutputList[[t]], params = "beta.dens")
    title(years[t],line = 0.5)
  }
  
  dev.off()
}#do all

## -----------------------------------------------------------------------------