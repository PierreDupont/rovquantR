## ------ IMPORT REQUIRED LIBRARIES ------
rm(list=ls())
gc()
library(raster)
library(coda)
library(nimble)
#library(spdep)
#library(maptools)
library(stringr)
library(abind)
library(R.utils)
library(adehabitatHR)
library(sf)
library(fasterize)
library(nimbleSCR)
library(dplyr)
library(readxl)
library(spatstat)
library(ggplot2)
library(stars)



## ------ SET REQUIRED WORKING DIRECTORIES ------
source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")



## ------ SOURCE THE REQUIRED FUNCTIONS ------
sourceDirectory(dir.function, modifiedOnly = FALSE)
sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)
load(file.path(dir.dropbox,"DATA/MISC DATA/age.lookup.table.RData"))
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



## ------ SET ANALYSIS CHARACTERISTICS -----

myVars <- list(
  ## WORKING DIRECTORY & MODEL NAME
  WD = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/2024",
  
  modelName = "53.aJ_FaCleaned2024_cleaned",
  
  ## HABITAT SPECIFICATIONS
  HABITAT = list( countries = c("SWE","NOR"),
                  habResolution = 20000,
                  habBuffer = 60000),
  
  ## NGS DATA SPECIFICATIONS
  DATA = list( years = 2014:2023,
               species = c("Jerv"),               
               sex = c("Hunn"),                   
               samplingMonths = list(12,1:6)),  
  
  ## DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 2000,
                    detResolution = 10000,
                    detDeadResolution = 15000),
  
  ## DATA GENERATION
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 3,
                     aug.factor = 0.8),
  
  ## OUTPUT PLOTS
  OUTPUT = list(mapResolution = 10000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))
threshold_month <- unlist(myVars$DATA$samplingMonths)[1]


if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(myVars$WD, myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}

## Set data directory
data.dir <- file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023")



################################################################################
##-- Renaming list ----
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
  CITES_number = "CITES-nummer",
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
  Field_personnel = "Feltpersonell",
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
  RovbaseID = "RovbaseID (Analyse)",
  RovbaseID_sample = "RovbaseID (Prøve)",
  Sample_status = "Prøvestatus",
  Sample_type = "Prøvetype",
  Sensitivity = "Følsomhet",
  Sex_analysis = "Kjønn (Analyse)",
  Sex = "Kjønn (Individ)",
  Sex = "Kjønn",
  Site_quality = "Stedkvalitet",
  Species_sample = "Art (Prøve)",
  Species = "Art (Analyse)",
  Species = "Art",
  Sensitivity = "Følsomhet",
  Site_quality = "Stedkvalitet",
  SVAID = "SVAID",
  Time_of_death = "Dødstidspunkt",
  Tips_name = "Tipser - Navn",
  Tips_phone = "Tipser - Telefon",
  Tips_email = "Tipser - E-post",
  Tips_role = "Tipser - Rolle",
  Tissue_sample = "Vevsprøve tatt",
  Uncertain_date = "Usikker dødsdato",
  Weight_slaughter = "Slaktevekt",
  Weight_total =  "Helvekt")



######################## ##
##  cleanRovBaseData()   ##
######################## ##

## ------ I. LOAD & CLEAN DATA ------

## ------   1. LOAD HABITAT DATA ------

##-- Load pre-processed habitat shapefiles
data(GLOBALMAP, envir = environment()) 
#data(COUNTIES, envir = environment()) 
data(COUNTRIES, envir = environment()) 

##-- CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
myStudyArea <- COUNTRIES %>%
  filter(ISO %in% myVars$HABITAT$countries) %>%
  mutate(id = 1) %>%
  group_by(id) %>%
  summarize()



## ------   2. LOAD ROVBASE DATA ------

## ------     2.1. NGS ------

##-- Load .xlsx 
DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                        extension = ".xls",
                                        pattern = "dna")) %>% 
  ##-- Rename columns to facilitate manipulation
  rename(., any_of(rename.list)) %>%
  ##-- Remove DEAD entries from the DNA data 
  filter(., !substr(RovbaseID_sample,1,1) %in% "M")



## ------     2.2. DEAD RECOVERIES ------

##-- Load .xlsx 
DEAD <- suppressWarnings(readMostRecent( path = data.dir,
                                         extension = ".xls",
                                         pattern = "dead")) %>% 
  ##-- Rename columns to facilitate manipulation
  rename(., any_of(rename.list)) %>%
  ##-- Check dead recovery weights
  mutate(.,
         Weight_slaughter = as.numeric(gsub(",", ".", as.character(Weight_slaughter))),
         Weight_total = as.numeric(gsub(",", ".", as.character(Weight_total))),
         Weight = ifelse(is.na(Weight_total), Weight_slaughter, Weight_total),
         Weight = ifelse(is.na(Weight),-999,Weight)) %>%
  ##-- Remove un-verified dead recoveries [HB]
  filter(.,!grepl( pattern = "Påskutt",
                   x = as.character(Outcome)))



## ------   3. CLEAN DATA ------

## ------     3.1. MERGE NGS & DEAD RECOVERIES ------

myData <- merge( DEAD,
                 DNA,
                 by = c("Id","RovbaseID","DNAID","Species",
                        "Sex","Date","East_UTM33","North_UTM33",
                        "County"),
                 all = TRUE)



## ------     3.2. FILTER UNUSABLE SAMPLES ------

myData <- myData %>%
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
    Year = ifelse( Month < threshold_month,
                   Year-1,
                   Year),
    ##-- Fix unknown "Id"
    Id = ifelse(Id %in% "", NA, Id),
    ##-- Fix unknown "Sex"
    Sex = ifelse(Sex %in% "Ukjent", NA, Sex),
    ##-- Determine Death and Birth Years
    Age = suppressWarnings(as.numeric(as.character(Age))),
    RovbaseID = as.character(RovbaseID),
    Death = ifelse(substr(RovbaseID,1,1) %in% "M", Year, NA),
    Birth = Death - Age) %>%
  dplyr::filter(., 
                ##-- Filter out samples with no ID
                !is.na(Id),
                ##-- Filter out samples with no Coordinates
                !is.na(East_UTM33),
                ##-- Filter out samples with no dates  
                !is.na(Year)) %>%
  droplevels(.)



## ------     3.3. AGE ------

## Reconstruct minimal & maximal ages
if(!is.null(age.lookup.table)){
  myData$age.label <- myData$Age
  
  myData <- merge( myData, age.lookup.table[ ,-1],
                   by.x = "age.label",
                   by.y = "age.label",
                   all.x = TRUE)                          ## Merge with info from lookup table
  
  ## FILL IN THE REST OF THE AGES FROM FOR NUMERIC RECORDS
  numeric.age.records <- which(!is.na(as.numeric(as.character(myData$age.label))) & !is.na(myData$age.label))
  myData[numeric.age.records, c("min.age","max.age","age")] <- floor(as.numeric(as.character(myData$Age[numeric.age.records])))
}



## ------     3.4. EXTRACT COUNTRY ------

## Convert samples coordinates to the correct spatial projection
myData <- myData %>%
  sf::st_as_sf(., coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(.,sf::st_crs(32633)) 

myData$Country_sf[!is.na(as.numeric(st_intersects(myData, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
myData$Country_sf[!is.na(as.numeric(st_intersects(myData, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"



## ------     3.5. FACTOR TRICK ------

### [PD] not needed anymore. using dropLevels() instead

# myData$Id <- factor( x = as.character(myData$Id),
#                      levels = unique(as.character(myData$Id)))



## ------     3.6. CHECK SEX ASSIGNMENT ------

##-- List all individual IDs
ID <- unique(as.character(myData$Id))
myData$Sex <- as.character(myData$Sex)

##-- Initialize the vector of IDs with conflicting sexes
IdDoubleSex <- 0
counter <- 1
for(i in 1:length(ID)){
  ## Subset data to individual i
  tmp <- myData$Sex[myData$Id == ID[i]] 
  ## Create a table of the number of times individual i was assigned to each sex
  tab <- table(tmp[tmp %in% c("Hunn","Hann")])
  ## If conflicting sexes (ID identified as both "Hunn" and "Hann")
  if(length(tab) == 2){
    ## If ID assigned the same number of times to the 2 sexes, assign to Ukjent
    if(tab[1] == tab[2]){
      myData$Sex[myData$Id == ID[i]] <- "Ukjent"
    } else {
      ## Otherwise pick the most common sex
      myData$Sex[myData$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
    }
    ## In any case, print a warning
    print(paste("Warnings!!!", "Individuals", ID[i],
                "assigned to both sexes. Now assigned to",
                names(tab)[which(tab == max(tab))])) 
    IdDoubleSex[counter] <- ID[i]
    counter <- counter + 1
  }
  ## If only one of "Hunn" or "Hann" registered
  if(length(tab) == 1){myData$Sex[myData$Id == ID[i]] <- names(tab)}
  ## If anything else registered : "Ukjent"
  if(length(tab) == 0){myData$Sex[myData$Id == ID[i]] <- "Ukjent"}
}#i



## ------     3.7. SPATIAL FILTER ------

## Remove all samples outside the polygon of interest
myData <- myData[!is.na(as.numeric(st_intersects(myData, myStudyArea))), ]



## ------     3.8. SPLIT NGS & DEAD RECOVERY ------

myFullData.sp <- list()
myFullData.sp$dead.recovery <- myData[!is.na(myData$Death), ] 
myFullData.sp$alive <- myData[is.na(myData$Death), ] 

IdDoubleDead <- myFullData.sp$dead.recovery$Id[duplicated(myFullData.sp$dead.recovery$Id)]



## ------     4.9. REMOVE SUSPECT NGS ACCORDING TO HENRIK ------

##-- DNA samples to be removed from Henrik
SUSPECT_NGS_SAMPLES <- read.csv( file.path(data.dir, "Remove ngs samples list wolverine 2024.csv"),
                                 fileEncoding = "latin1") 

myFullData.sp$alive <- myFullData.sp$alive %>%
  filter(.,!as.character(DNAID) %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB))



## ------     4.10. REMOVE SUSPECT DEAD RECOVERIES ACCORDING TO HENRIK ------

##-- DEAD RECOVERIES to be removed from Henrik
SUSPECT_DeadRecoSAMPLES <- read.csv( file.path(data.dir, "Remove dead recoveries list wolverine 2024.csv"),
                                     fileEncoding = "latin1")

myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
  ##-- Dead recoveries in suspect samples list from HB
  filter(.,
         ##-- Remove samples flagged by Henrik
         !(as.character(RovbaseID) %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)),
         ##-- Dead recoveries flagged by Henrik that should always be removed (email from the 18/12/2024)
         !(RovbaseID %in% c("M495994","M524051","M524052","M524053")),
         ##-- Remove dead recoveries of pups killed before recruitment
         !(Age_class %in% "Unge" & Month > 2 & Month < 12),
         ##-- Remove dead recoveries of pups based on weight
         !(Weight > 0 & Weight < 4 & Month < 12 & Month > 2),
         ##-- Remove potential duplicates 
         !(duplicated(Id) & is.na(Death_cause))) 


# ## Remove individuals that died twice
# ## [CM] TO BE CHECKED BECAUSE "length(IdDoubleDead) < 0" and so it was deactivated
# IdDoubleDead <- myFullData.sp$dead.recovery$Id[duplicated(myFullData.sp$dead.recovery$Id)]
# if(length(IdDoubleDead) > 0){
#   duplicatedDeath <- NULL
#   for(i in IdDoubleDead){
#     tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$Death_cause))
#     if(length(tmp) == 0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-2]}##[CM] remove the second record.
#     duplicatedDeath <- c(duplicatedDeath, tmp)
#   }#i
#   myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]
# }#if
#
# ## Remove pups killed before recruitment based on weight (cf. Henrik)
# ## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
# sum(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
#       myFullData.sp$dead.recovery$Month > 2 &
#       myFullData.sp$dead.recovery$Month < 12)
# 
# myFullData.sp$dead.recovery <-  myFullData.sp$dead.recovery %>%
#   filter(.,
#          Age_class %in% "Unge" & Month > 2 & Month < 12)
#
# ## 2) remove individuals that have a weight > 0 and < 4 between March and November
# ## format the weight correctly
# myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
# myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)
# 
# ## convert to decimals
# myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
# myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))
# ## get the two weight columns together.
# myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt),
#                                              myFullData.sp$dead.recovery$Helvekt,
#                                              myFullData.sp$dead.recovery$Slaktevekt)
# ## Assign negative values to NAs to avoid issues
# myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999
#
# ## [PD] was it checked ?
# ## Check with Henrik
# ## this step does not remove dead recoveries on id with weight==0 should it?
# 
# ## Check how many dead reco we remove and remove if more than 0
# if(sum(myFullData.sp$dead.recovery$weight > 0 &
#        myFullData.sp$dead.recovery$weight < 4 &
#        myFullData.sp$dead.recovery$Month < 12 &
#        myFullData.sp$dead.recovery$Month > 2) > 0){
#   myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
#                                                                       myFullData.sp$dead.recovery$weight < 4 &
#                                                                       myFullData.sp$dead.recovery$Month < 12 &
#                                                                       myFullData.sp$dead.recovery$Month > 2), ]
# }
# 
# ## Check how many dead reco with a weight of 0 kg and recovered between march and november
# if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
#        myFullData.sp$dead.recovery$Month < 12 &
#        myFullData.sp$dead.recovery$Month > 2) > 0){
#   myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
#                                 myFullData.sp$dead.recovery$Month < 12 &
#                                 myFullData.sp$dead.recovery$Month > 2,  ]
# }



## ------   5. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

## ------     5.1. LOAD GPS SEARCH TRACKS ------
# 
# TRACKS1 <- read_sf(file.path( dir.dropbox,
#                               "DATA/RovbaseData/TRACK DATA FROM BOUVET 20240830",
#                               "XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp"))
# # whichEquals <- st_equals( TRACKS1, remove_self = TRUE)
# # whichHaveEquals <- which(!unlist(lapply(whichEquals, is.empty)))
# # length(whichHaveEquals)
# # 
# # id <- whichHaveEquals[2]
# # 
# # plot(st_geometry(TRACKS1[id, ]), col = "red")
# # plot(st_geometry(TRACKS1[whichEquals[[id]], ]), col = "blue")
# # 
# # as.data.frame(st_drop_geometry(TRACKS1[id, ]))
# # as.data.frame(st_drop_geometry(TRACKS1[whichEquals[[id]], ]))
# 
# 
# TRACKS2 <- read_sf(file.path(dir.dropbox,
#                              "DATA/RovbaseData/TRACK DATA FROM BOUVET 20240830",
#                              "XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))
# # whichEquals <- st_equals( TRACKS2, remove_self = TRUE)
# # whichHaveEquals <- which(!unlist(lapply(whichEquals, is.empty)))
# # length(whichHaveEquals)
# # 
# # id <- whichHaveEquals[2]
# # 
# # plot(st_geometry(TRACKS1[id, ]), col = "red")
# # plot(st_geometry(TRACKS1[whichEquals[[id]], ]), col = "blue")
# # 
# # as.data.frame(st_drop_geometry(TRACKS1[id, ]))
# # as.data.frame(st_drop_geometry(TRACKS1[whichEquals[[id]], ]))
# 
# 
# ## COMBINE ALL TRACKS
# TRACKS <- rbind(TRACKS1, TRACKS2) %>% 
#   ## Set years according to monitoring season
#   mutate( Year = ifelse( Mth < threshold_month, 
#                          Yr,
#                          Yr - 1),
#           RovbaseID = as.character(RovbaseID)
#           ) %>%
#   filter(.,
#          ## REMOVE HELICOPTER TRACKS
#          Helikopter == "0",
#          ## KEEP ONLY WOLVERINE TRACKS
#          Jerv == "1",
#          ## Keep only sampling months
#          Mth %in% unlist(myVars$DATA$samplingMonths),
#          ## Keep only sampling years
#          Year %in% years) %>%
#   ## Keep tracks in the study area only
#   st_intersection(., st_as_sf(myStudyArea)) %>%
#   ## Count number of search tracks per year
#   group_by(Year) %>%
#   mutate(trackNum = 1:length(Year)) %>%
#   ungroup()
#   
# rm(TRACKS1)
# rm(TRACKS2)
# 
# save(TRACKS, file = file.path(myVars$WD, myVars$modelName, "tracks.RData"))
load(file.path(myVars$WD, myVars$modelName, "tracks.RData"))



## ------     5.2. ASSIGN SAMPLES TO TRACKS ------

# test <- assignSearchtTracks(
#   data = myFilteredData.sp$alive,
#   tracks = TRACKS_YEAR2,
#   dist = 750,
#   progress.bar = T)

myFilteredData.sp <- myFullData.sp
myFilteredData.sp$alive$Track_RovbaseID <- NA
myFilteredData.sp$alive$Track_Dist <- NA

## CREATE A BUFFER AROUND EACH DETECTION
tmp <-  st_buffer( myFilteredData.sp$alive,
                   dist = 750)
## LOOP OVER NGS 
for(i in 1:nrow(myFilteredData.sp$alive)){
  ## Subset to tracks run on the same date
  tmpTRACKS <- TRACKS %>%
    filter(as.character(Dato) %in% as.character(myFilteredData.sp$alive$Date[i]))
  
  ## Extract tracks within the buffered point circle
  isWithinDistance <- st_intersection(tmpTRACKS, tmp[i,])
  
  ## If no tracks within the threshold distance, move to next sample..
  if(nrow(isWithinDistance) < 1){
    next
  } else {
    ## ...else, assign sample to the closest track
    dist <- st_distance( myFilteredData.sp$alive[i, ],
                         isWithinDistance,
                         by_element = F)
    myFilteredData.sp$alive$Track_RovbaseID[i] <- isWithinDistance$RovbaseID[which.min(dist)]
    myFilteredData.sp$alive$Track_Dist[i] <- min(dist)
  }
}#i



## ------     5.3. SPLIT DATA TO OPPORTUNISTIC & STRUCTURED ------

## Check if some hair samples are assigned to structured
## HAIR TRAP SAMPLES from Henrik
HairTrapSamples <- read_xlsx(file.path(data.dir, "hairtrapsNB2024.xlsx")) 

## Distance threshold
distanceThreshold <- 500

## Proeveleverandoer columns was replaced by two columns, merging them now...
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  mutate(.,
         Collector = ifelse(is.na(Collector_other_role), 
                            Collector_role,
                            Collector_other_role),
         Structured = Collector %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
           !is.na(Track_RovbaseID) &
           Track_Dist <= distanceThreshold,
         HairTrap = DNAID %in% HairTrapSamples$DNAID)

if(any(myFilteredData.sp$alive$Structured & myFilteredData.sp$alive$HairTrap)){
  print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
}



## ------   6. SEPARATE MORTALITY CAUSES -----

##-- List mortality causes
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))

##-- List of patterns to match for legal causes of mortality
toMatch <- c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske")

##-- Identify legal death
myFilteredData.sp$dead.recovery <- mutate( 
  myFilteredData.sp$dead.recovery,
  legal = grepl(paste(toMatch, collapse = "|"), Death_cause))

# ## DEFINE LEGAL MORTALITY
# legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
# legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
# 
# ## SPLIT MORTALITY CAUSES
# legal.death <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]
# Other.death <- myFilteredData.sp$dead.recovery[!myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]



## ------   7. SAVE filtered data FOR FASTER LOADING -----

save( myFilteredData.sp,
      file = file.path(myVars$WD, myVars$modelName, "myFilteredData.RData"))
#load(file.path(myVars$WD, myVars$modelName, "myFilteredData.RData"))



################################################################################

######################## ##
##   makeRovQuantData()  ##
######################## ##

##-- Load pre-defined habitat rasters and shapefiles
data(COUNTRIES, envir = environment()) 
data(COUNTIES, envir = environment()) 
data(habitatRasters, envir = environment()) 
data(GLOBALMAP, envir = environment()) 

##-- Disaggregate habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/myVars$HABITAT$habResolution)

##-- Merge Norwegian counties for practical reasons
COUNTIES_AGGREGATED <- COUNTIES
COUNTIES_AGGREGATED$id <- 1:nrow(COUNTIES_AGGREGATED)
COUNTIES_AGGREGATED$id[c(24,3,15,9,14,38,40,21,27,37,31,26,34,5,8,12,36,13,7)] <- 3
COUNTIES_AGGREGATED$id[c(39,33,23,32,29,22,4,11,20,2,10,16,25,1)] <- 4
COUNTIES_AGGREGATED$id[c(19)] <- 1
COUNTIES_AGGREGATED$id[c(35)] <- 2
COUNTIES_AGGREGATED$id[c(17,28)] <- 5
COUNTIES_AGGREGATED$id[c(18)] <- 7
COUNTIES_AGGREGATED$id[c(30)] <- 8
COUNTIES_AGGREGATED <- COUNTIES_AGGREGATED %>%
  dplyr::group_by(id) %>%
  dplyr::summarise() %>%
  st_simplify(. ,
              preserveTopology = T,
              dTolerance = 500)



## ------   1. GENERATE HABITAT ------

## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------

load(file.path(myVars$WD, myVars$modelName, "myFilteredData.RData"))

##-- Determine study area based on NGS detections
##-- Buffer NGS detections and cut to Swedish and Norwegian borders
myStudyArea <- myFilteredData.sp$alive %>%  
  sf::st_buffer(., dist = myVars$HABITAT$habBuffer * 1.4) %>%
  sf::st_union() %>%
  sf::st_intersection(., COUNTRIES) %>%
  sf::st_as_sf() 

##-- Make habitat from predefined Scandinavian raster of suitable habitat
habitat <- MakeHabitatFromRastersf(
  poly = myStudyArea,
  habitat.r = habRaster,
  buffer = myVars$HABITAT$habBuffer,
  plot.check = FALSE) 

##-- Retrieve habitat windows boundaries
lowerHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1,] - 0.5*myVars$HABITAT$habResolution
upperHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1,] + 0.5*myVars$HABITAT$habResolution
numHabCells <- dim(lowerHabCoords)[1]

##-- Create habitat grid
habIDCells.mx <- habitat$IDCells.mx 
habIDCells.mx[] <- 0
scaledHabGridCenters <- scaleCoordsToHabitatGrid(
  coordsData = habitat$habitat.xy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =F )$coordsHabitatGridCenterScaled

scaledHabGridCenters <- scaledHabGridCenters[habitat$habitat.r[]==1, ]
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}#i



## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       1.2.1. DEN COUNTS ------

DEN <- read.csv(file.path(dir.dropbox,
                          "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/DEN_COUNTS_2009_2024_fromHB.csv"),
                fileEncoding = "latin1")

DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[ ,("id")]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD(as(DEN.sp,"Spatial"),
             h = 30000,
             grid = as(habitat$habitat.r, 'SpatialPixels'))))

##-- EXTRACT COVARIATES
denCounts <- DEN.r[habitat$habitat.r[ ]==1] %>%
  scale() %>% 
  round(., digits = 2)



## ------   2. GENERATE DETECTORS ------

## ------    2.1. GENERATE DETECTORS CHARACTERISTICS ------

##-- GENERATE NGS DETECTORS BASED ON THE STUDY AREA
habitat.subdetectors <- disaggregate(
  habitat$habitat.rWthBuffer,
  fact = res(habitat$habitat.r)[1]/myVars$DETECTORS$detSubResolution)

##-- GENERATE NGS DETECTORS BASED ON THE STUDY AREA
detectors <- MakeSearchGridsf(
  data = habitat.subdetectors,
  resolution = myVars$DETECTORS$detResolution,
  div = (myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)^2,
  plot = FALSE,
  fasterize = TRUE)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(detectors$main.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(detectors$main.detector.sp)
n.trials <- as.vector(table(detectors$detector.sp$main.cell.id))


## ------    2.2. IDENTIFY DETECTORS IN NORBOTTEN -----

COUNTIESAroundNorrbotten <- COUNTIES %>%
  filter(.,NAME_1 %in% c("Norrbotten","Troms","Vasterbotten","Nordland","Finnmark")) %>%
  st_simplify(., dTolerance = 500)

distDestsCounties <- st_distance( detectors$main.detector.sp,
                                  COUNTIESAroundNorrbotten,
                                  byid = T)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min) == 3) 


# ## Alternative...but does not account for detectors in Finland
# detsNorrbotten2 <- COUNTIES %>%
#   filter(.,NAME_1 %in% c("Norrbotten")) %>%
#   st_intersects(.,detectors$main.detector.sp) %>%
#   unlist() %>%
#   unique()




## ------    2.3. GENERATE DETECTOR-LEVEL COVARIATES ------

## ------       2.3.1. EXTRACT COUNTRIES ------

detCountries <- st_distance( detectors$main.detector.sp,
                             COUNTRIES,
                             by_element = F) %>%
  apply(., 1, function(x) which.min(x)) %>% 
  as.numeric(as.factor(.))



## ------       2.3.2. EXTRACT COUNTIES ------

dist <- st_distance(detectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F )
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties),]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))



## ------       2.3.3. EXTRACT GPS TRACKS LENGTHS ------

##-- Load pre-processed seach tracks
load(file.path(myVars$WD, myVars$modelName, "tracks.RData"))


##-- INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(detectors$main.detector.sp),
                                      rep(1,nrow(detectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf( stars::st_as_stars(detectorGrid.r), 
                              as_points = FALSE,
                              merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)

intersection2 <- st_intersection(detectorGrid, TRACKS) %>%
  mutate(LEN = st_length(.)) %>%
  st_drop_geometry() %>%
  group_by(id,Year) %>%
  summarise(transect_L = sum(LEN))

detTracks2 <- intersection2 %>%
  tidyr::pivot_wider( names_from = Year,
               values_from = transect_L)

##-- CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
TRACKS.r <- list()
for(t in 1:nYears){
  TRACKSst <- filter(TRACKS, Year == years[t])
  
  intersection <- st_intersection(detectorGrid, TRACKS) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN)) 
  
  ## Get total length searched in each detector grid cell
  # transect_N = length(unique(ID)))## Get total number of visits in each detector grid cell
  # transect_qi = mean(QI))     
  
  ##-- Get mean transects quality index for each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[ ,t]
}#t



## ------       2.3.4. EXTRACT DISTANCES TO ROADS ------

##-- LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(file.path(dir.dropbox,
                                 "DATA/GISData/Roads/MinDistAllRoads1km.tif"))
r <- fasterize(st_as_sf(myStudyArea), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea)

##-- AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

##-- EXTRACT DISTANCE TO ROAD FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)

##-- if NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( DistAllRoads,
                        detectors$main.detector.sp[isna, ],
                        buffer = 15000,
                        fun = mean,
                        na.rm = T)
detRoads[isna] <- tmp



## ------       2.3.5. EXTRACT DAYS OF SNOW ------

##-- SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
SNOW <- stack(file.path(dir.dropbox,
                        "DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2008_2024_Wolf.tif"))
##-- RENAME THE LAYERS
names(SNOW) <- paste(2008:2023,(2008:2023)+1, sep="_")

##-- SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste0("X", years, "_", years + 1)]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))

##-- EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = nYears)

det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

##-- if NA returns the average value of the cells within 15000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW,
                       det.sptransf[isna, ],
                       buffer = 15000,
                       fun = mean,
                       na.rm = T)
detSnow[isna,1:nYears] <- tmp



## ------       2.3.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------          2.3.6.1. SKANDOBS ------

##-- SkandObs
skandObs <- read_xlsx(file.path(dir.dropbox,
                                "DATA/Skandobs/RB_Skandobs_2012_2024/Richard_Bischof_Skandobs_2012_2024dd.xlsx"))
colnames(skandObs) <- translateForeignCharacters(
  dat = colnames(skandObs),
  dir.translation = dir.analysis)

##-- Get time 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))

##-- Make it spatial 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(myStudyArea))

##-- Subste based on season 
subset <- skandObs$month %in% c(unlist(myVars$DATA$samplingMonths))
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) #--- need to change for other species
skandObs <- skandObs[subset,] 

##-- Subset based on space 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(habitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1,]

subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace, ] 

##-- Rasterize at the detector level
r.detector <- aggregate(habitat.subdetectors, fact=(myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution))
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



## ------          2.3.6.2. ROVBASE ------

##-- RovBase Observations
rovbaseObs1 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB2810202415264376.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152348493.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152447860.xlsx"))
rovbaseObs4 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152538742.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3,rovbaseObs4)

##-- Get all samples collected
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`), ]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

colnames(rovbaseObs) <- translateForeignCharacters(
  dat = colnames(rovbaseObs),
  dir.translation = dir.analysis)

rovbaseObs$Proevetype <- translateForeignCharacters(
  dat = rovbaseObs$Proevetype,
  dir.translation = dir.analysis)

##-- Define projections
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)

##-- Subset data 
filter <- list(
  species = "Jerv",
  type = c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)",
           "Sekret (Jerv)","Saliv/Spytt"),
  month = unlist(myVars$DATA$samplingMonths))

##-- Subset month & type of sample
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset, ] 

##-- Subset if sample was successfully genotyped and from the focal speciesS 
subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$`Art (Proeve)`) 
rovbaseObs.sp <- rovbaseObs.sp[-subset, ] 

##-- Subset based on space 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace, ] 

##-- Rasterize
r.detector <- aggregate( habitat.subdetectors,
                         fact=(myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution))
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



## ------          2.3.6.3. COMBINE ROVBASE & SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <-    1
}



## ------          2.3.6.4. SMOOTH THE BINARY MAP ------

## we tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
habOwin <- as.owin(as.vector(extent(r.detector)))
ds.list <- lapply(years,function(y){
  ## ROVBASE DATA 
  pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y,]
  ## SKANDOBS
  pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y,] )
  ## SMOOTH AND RASTERIZE
  p <- ppp(pts[,1], pts[,2], window=habOwin)
  ds <- density(p, adjust=0.02) #---change bandwith (smoothing) with "adjust
  ds <- raster(ds)
  
  ds <- ds1 <- raster::resample(ds, r.detector) #mask(ds,rasterToPolygons(habitat.list$habitat.rWthBuffer,function(x) x==1))
  threshold <- 0.1 / prod(res(ds)) #--number per 1 unit of the projected raster (meters)
  ds1[] <- ifelse(ds[]<threshold,0,1)
  ds1 <- mask(ds1, habitat.rWthBufferPol)
  ds <- mask(ds, habitat.rWthBufferPol)
  
  return(list(ds,ds1))
})

ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))
names(ds.brick) <- years



## ------          2.3.6.5. COLOR CELLS WHERE HAIR TRAPPING OCCURED ------

##-- IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.sp$alive %>%
  filter(HairTrap)

##-- MANUALLY FIND HAIR SMAPLES AND COLOR THE CELL.
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr)
  whereHair <- raster::extract( r.SkandObsOtherSamplesBinary[[t]],
                                tmpHair,
                                cellnumbers = T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[ ,1]] <- 1
}#i


#IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.sp$alive[which(myFilteredData.sp$alive$DNAID %in% HairTrapSamples$DNAID), ]

myFilteredData.sp$alive$Date[which(myFilteredData.sp$alive$DNAID %in% HairTrapSamples$DNAID)]

#MANNUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL. 
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr)
  whereHair <-raster::extract(r.SkandObsOtherSamplesBinary[[t]],tmpHair,cellnumbers=T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[,1]] <- 1
  plot(r.SkandObsOtherSamplesBinary[[t]])
  plot(tmpHair$geometry,add=T,col="red")
}



## ------          2.3.6.6. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract( r.SkandObsOtherSamplesBinary,
                                                detectors$main.detector.sp)



## ------       2.3.7. SCALE AND ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],2))
detCovs[,,1] <- detTracks
detCovs[,,2] <- detSnow

detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],3))
detCovsOth[,,1] <- detSnow
detCovsOth[,,2] <- matrix(detRoads,length(detRoads),nYears)
detCovsOth[,,3] <- detOtherSamples

## CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){
  print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")
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

# ##-- Create habitat grid
# habIDCells.mx <- habitat$IDCells.mx
# habIDCells.mx[] <- 0
# for(i in 1:nrow( habitat$scaledCoords)){
#   habIDCells.mx[trunc( habitat$scaledCoords[i,2])+1,
#                 trunc( habitat$scaledCoords[i,1])+1] <- i
# }



## ------   4. CREATE LOCAL OBJECTS -----

##-- Get local habitat windows
habitat$localObjects <- nimbleSCR::getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = habitat$scaledCoords,
  dmax = habitat$maxDist/habitat$resolution,
  plot.check = F)

##-- Get local detectors
detectors$localObjects <- nimbleSCR::getLocalObjects(
  habitatMask = habitat$habitat.mx,
  coords = detectors$scaledCoords,
  dmax = detectors$maxDist/detectors$resolution,
  plot.check = F)




## ------   5. GENERATE y DETECTION ARRAYS ------

## ------     5.1. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 

##-- Remove samples outside the STUDY AREA 
myStudyArea$idd <- 1
myStudyAreaAggregated <- myStudyArea %>% group_by(idd) %>% summarize()

whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive,
                                                   myStudyAreaAggregated))))
if(length(whichOut)>0){
  myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
}
myFilteredData.sp$alive$Id <- droplevels(myFilteredData.sp$alive$Id)


##-- Remove dead recoveries outside the BUFFERED HABITAT
whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery,
                                                       habitat$buffered.habitat.poly))))
if(length(whichOutBuff)>0){
  myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]
}
myFilteredData.sp$dead.recovery$Id <- droplevels(myFilteredData.sp$dead.recovery$Id)



## ------     5.2. FILTER NGS & DEAD RECOVERY DATA FOR YEARS ------

##-- Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]



## ------     5.3. FILTER NGS DATA FOR MONTHS ------

##-- Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in%
                                                     unlist(myVars$DATA$samplingMonths), ]



## ------     5.4. FILTER DATA IN NORRBOTTEN IN ALL YEARS EXCEPT 2017, 2018 and 2019 ------ 

##-- Norrbotten county
COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018,2023)
is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))

##-- Check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(is.Norr)), ]$Year)

##-- Subset 
myFilteredData.sp$alive <- myFilteredData.sp$alive[- which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                                             !is.na(is.Norr)), ]



## ------     5.5. FILTER NGS & DEAD RECOVERY DATA FOR SEX ------

##-- myFilteredData
myFilteredData.spAllSex <- myFilteredData.sp
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Sex %in% myVars$DATA$sex,]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Sex %in% myVars$DATA$sex,]



## ------     5.6. ASSIGN SAMPLES TO DETECTORS ------

##-- ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$alive,                
  myDetectors = detectors.dead$main.detector.sp,
  mysubDetectors = detectors.dead$detector.sp,
  radius = myVars$DETECTORS$detResolution)

##-- DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$dead.recovery,
  myDetectors = detectors.dead$main.detector.sp,
  radius = myVars$DETECTORS$detResolution)

### MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET ASSIGNED 
### TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING.
### FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM TO THE CLOSEST DETECTOR 
### OUTSIDE OF NORRBOTTEN.
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)
whichdets <- which(!myData.alive$myData.sp$Year %in% yearsSampledNorrb &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)

##-- ASSIGN DETECTORS 
subdetector.sf <- detectors$detector.sp

##-- ALL
for(i in 1:length(whichdets)){
  tmp <- myData.alive$myData.sp[whichdets[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, detectors$main.detector.sp)
  #Artificially increase distance for detectors in Norrbotten 
  dist[,detsNorrbotten] <- 500000
  idMain <- which.min(dist[1,])
  myData.alive$myData.sp$Detector[whichdets[i]] <- idMain 
  # SUBDETECTOR
  dist <- st_distance(st_as_sf(tmp), subdetector.sf)
  #Artificially increase distance for detectors in Norrbotten
  dist[,detsNorrbotten] <- 500000
  idSub <- which.min(dist[1,])
  myData.alive$myData.sp$sub.detector[whichdets[i]] <- idSub 
}

##-- SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)




## ------     5.7. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ------

##-- STRUCTURED
y.arStruc <- MakeYsf( myData = myData.aliveStruc$myData.sp,
                      myDetectors = detectors$main.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      myDetectors2 = detectors.dead$main.detector.sp,
                      returnIdvector = TRUE)
y.ar.ALIVEStruc <- y.arStruc$y.ar
dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)

##-- OTHERS
y.arOth <- MakeYsf( myData = myData.aliveOthers$myData.sp,
                    myDetectors = detectors$main.detector.sp,
                    method = "Binomial",
                    returnIdvector = TRUE)
y.ar.ALIVEOth <- y.arOth$y.ar
dimnames(y.ar.ALIVEOth) <- dimnames(y.arOth$y.ar)

##-- MAKE SURE THE Y HAVE THE SAME DIMENSIONS
y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0

##-- FILL IN THE Y ARRAYS 
y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]],,] <- y.ar.ALIVEOth
y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]],,] <- y.ar.ALIVEStruc

##-- PROJECT DEATHS TO THE NEXT OCCASION.
y.ar.DEADProjected <- y.ar$y.ar2 
y.ar.DEADProjected[] <- 0
for(t in 2:nYears){y.ar.DEADProjected[ , ,t] <- y.ar$y.ar2[ , ,t-1]}#t
y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
y.ar.DEAD[y.ar.DEAD > 0] <- 1



## ------     5.8. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------

distances <- list()
for(t in 1:nYears){
  print(file.path("------ ", t ," -------", sep = "" ))
  distances[[t]] <- CheckDistanceDetectionsV2sf(
    y = y.ar.ALIVE[,,t], 
    detector.xy = detector.xy, 
    max.distance = myVars$DETECTIONS$maxDetDist,
    method = "pairwise",
    plot.check = T)
  
  ## PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
  if(myVars$plot.check){
    par(mfrow = c(1,1))
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      count <- 0
      for(i in affected.ids){
        count <- count+1
        plot(st_geometry(myStudyArea), main = file.path("t: ",t,"     i: ", names(affected.ids)[count]))
        scalebar(2*myVars$DETECTIONS$maxDetDist, xy = c(800000,6700000), type = "bar", divs = 2, below = "km",
                 label = c(0, myVars$DETECTIONS$maxDetDist/1000, myVars$DETECTIONS$maxDetDist/500), cex = 0.8, adj = c(0.5,-0.9))
        plot(st_geometry(COUNTRIES), add = T)
        plot(st_geometry(detectors$main.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar.ALIVE)[[1]][i] &
                                         myFilteredData.sp$alive$Year == years[t], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
        arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
               length = 0.1, lwd = 1)
        plot(st_geometry(detectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
        
        tmp2 <- detectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
        plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
      }#i
    }#if
  }#if plot.check
  
  ## REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
  y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  
  ## REMOVE DETECTIONS ALSO IN MYDATA TO RUN getSinits
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i], ] > 0)
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]),]
  }#i
}#t



## ------     5.9. TRAP-RESPONSE ------

##-- Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCovsf(
  data = myFullData.sp$alive,
  data.dead = myFullData.sp$dead.recovery)

##-- Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]

##-- Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]



## ------   6. MAKE AUGMENTATION ------

##-- DATA ARRAYS
y.alive <- MakeAugmentation( y = y.ar.ALIVE,
                             aug.factor = myVars$DETECTIONS$aug.factor,
                             replace.value = 0)
y.dead <- MakeAugmentation( y = y.ar.DEAD,
                            aug.factor = myVars$DETECTIONS$aug.factor,
                            replace.value = 0)
y.aliveOthers <- MakeAugmentation( y = y.ar.ALIVEOthers,
                                   aug.factor = myVars$DETECTIONS$aug.factor,
                                   replace.value = 0)
y.aliveStructured <- MakeAugmentation( y = y.ar.ALIVEStructured,
                                       aug.factor = myVars$DETECTIONS$aug.factor,
                                       replace.value = 0)

##-- INDIVIDUAL COVARIATES
already.detected <- MakeAugmentation( y = already.detected,
                                      aug.factor = myVars$DETECTIONS$aug.factor,
                                      replace.value = 1)



## -----------------------------------------------------------------------------
## ------ III.MODEL SETTING & RUNNING ------- 

## ------   1. NIMBLE MODEL ------

modelCode <- nimbleCode({
  ##------ SPATIAL PROCESS ------
  
  dmean ~ dunif(0,100)
  lambda <- 1/dmean
  betaDens  ~ dnorm(0.0,0.01)
  
  habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  for(i in 1:n.individuals){
    sxy[i, 1:2, 1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    for(t in 2:n.years){
      sxy[i, 1:2, t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
        upperCoords = upperHabCoords[1:numHabWindows,1:2],
        s = sxy[i,1:2,t-1],
        lambda = lambda,
        baseIntensities = habIntensity[1:numHabWindows],
        habitatGrid =  habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max,
        numWindows= numHabWindows)
    }#i  
  }#t
  
  
  ##----- DEMOGRAPHIC PROCESS ----- 
  
  omeg1[1:2] ~ ddirch(alpha[1:2])   
  
  for(t in 1:n.years1){
    # PRIORS 
    gamma[t] ~ dunif(0,1)
    phi[t] ~ dunif(0,1)
    
    # "UNBORN"
    omega[1,1,t] <- 1-gamma[t]
    omega[1,2,t] <- gamma[t]
    omega[1,3,t] <- 0
    # "Alive"
    omega[2,1,t] <- 0
    omega[2,2,t] <- phi[t]
    omega[2,3,t] <- 1-phi[t]
    # "Dead"
    omega[3,1,t] <- 0
    omega[3,2,t] <- 0
    omega[3,3,t] <- 1
  }#t
  
  
  
  for(i in 1:n.individuals){ 
    detResponse[i,1] ~ dbern(pResponse)
    z[i,1] ~ dcat(omeg1[1:2]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
    }#i 								
  }#t 
  
  
  ##----- DETECTION PROCESS -----## 
  
  pResponse ~ dunif(0,1)
  for(i in 1:n.individuals){ 
    detResponse[i,1] ~ dbern(pResponse)
  }#i
  
  for(t in 1:n.years){
    sigma[t] ~ dunif(0,4)
    
    betaResponse[t] ~ dunif(-5,5)
    betaResponseOth[t] ~ dunif(-5,5)
    
    for(c in 1:n.covs){
      betaCovs[c,t] ~ dunif(-5,5)
    }#c
    for(c in 1:n.covsOth){
      betaCovsOth[c,t] ~ dunif(-5,5)
    }#c
    for(c in 1:n.counties){
      p01[c,t] ~ dunif(0,1)
      p0[c,t] <- p01[c,t] *countyToggle[c,t]## toggle counties
    }#c
    for(c in 1:n.countries){
      p01Oth[c,t] ~ dunif(0,1)
      p0Oth[c,t] <- p01Oth[c,t] *countyToggleOth[c,t]## toggle counties
    }#c
    
    for(i in 1:n.individuals){
      y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCovResponse(
        sxy = sxy[i,1:2,t],
        sigma = sigma[t],
        nbDetections[i,t],
        yDets = yDets[i,1:nMaxDetectors,t],
        detector.xy =  detector.xy[1:n.detectors,1:2],
        trials = trials[1:n.detectors],
        detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
        ResizeFactor = ResizeFactor,
        maxNBDets = maxNBDets,
        habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        p0[1:n.counties,t],
        detCounties[1:n.detectors],
        detCov = detCovs[1:n.detectors,t,1:n.covs],
        betaCov = betaCovs[1:n.covs,t],
        BetaResponse = betaResponse[t],
        detResponse = detResponse[i,t])
      
      y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbin_LESS_Cached_MultipleCovResponse(
        sxy = sxy[i,1:2,t],
        sigma = sigma[t],
        nbDetectionsOth[i,t],
        yDets = yDetsOth[i,1:nMaxDetectorsOth,t],
        detector.xy =  detector.xy[1:n.detectors,1:2],
        trials = trials[1:n.detectors],
        detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
        ResizeFactor = ResizeFactor,
        maxNBDets = maxNBDets,
        habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        p0Oth[1:n.countries,t],
        detCountries[1:n.detectors,t],
        detCov = detCovsOth[1:n.detectors,t,1:n.covsOth],
        betaCov = betaCovsOth[1:n.covsOth,t],
        BetaResponse = betaResponseOth[t],
        detResponse = detResponse[i,t])
    }#i
  }#t
  
  
  
  ##------ DERIVED PARAMETERS ------
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){ 
      isAlive[i,t] <- (z[i,t] == 2) 
    }#i
    N[t] <- sum(isAlive[1:n.individuals,t])
  }#t
  
})



## ------   2. NIMBLE CONSTANTS ------

nimConstants$countyToggle <- nimInits$p01
nimConstants$countyToggle[] <- 1

yearsNotSampled <- which(!years%in% yearsSampledNorrb)
for(t in yearsNotSampled){
  nimConstants$countyToggle[1,t] <- 0
}

nimConstants$countyToggleOth <- nimInits$p01Oth
nimConstants$countyToggleOth[] <- 1
yearsNotSampled <- which(!years%in% yearsSampledNorrb)
for(t in yearsNotSampled){
  nimConstants$countyToggleOth[3,t] <- 0
}

nimConstants <- list( n.individuals = dim(y.alive)[1],
                      n.detectors = dim(y.alive)[2],  
                      n.years = dim(y.alive)[3], 
                      n.years1 = dim(y.alive)[3]-1, 
                      n.covs = dim(detCovs)[3],
                      n.covsOth = dim(detCovsOth)[3],
                      numHabWindows = numHabCells,
                      n.countries = max(detCountries) + 1,
                      n.counties = max(detCounties),
                      y.max = dim(habIDCells.mx)[1],
                      x.max = dim(habIDCells.mx)[2])



## ------   3. NIMBLE DATA ------

## Add another category to detcountry if in Norrbotten, to turn off detection to 0 there. 
detCountriesNorb <- matrix(NA, nrow=length(detCountries),ncol=nYears)
detCountries1 <- detCountries
detCountries1[detCounties %in% 1] <- 3
for(t in 1:nYears){
  if(t %in% yearsNotSampled){
    detCountriesNorb[,t] <- detCountries1
  }else{
    detCountriesNorb[,t] <- detCountries
  }
}  

nimData$detCountries <-  detCountriesNorb


## ------     3.1. GENERATE KNOWN z ------

z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, NA)
z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    zz[range.det[1]:range.det[2]] <- 2
  }
  return(zz)
}))



## ------     3.2. GENERATE detResponse -----

detResponse <- already.detected 
detResponse[rownames(detResponse) %in% "Augmented",1] <- NA



## ------     3.3. LIST DATA ------

nimData <- list( z = z,   
                 y.alive = y.alive,
                 lowerHabCoords = lowerHabCoords/1000, 
                 upperHabCoords = upperHabCoords/1000, 
                 detCounties = detCounties,
                 detCountries = detCountries,
                 detCovs = detCovs,
                 detCovsOth = detCovsOth,
                 detResponse = detResponse,
                 denCounts = denCounts,
                 trials = n.trials,
                 alpha = rep(1,2),
                 detector.xy = detector.xy/1000,
                 habitatGrid = habIDCells.mx )



## ------   4. NIMBLE PARAMETERS ------

nimParams <- c("N", 
               "omeg1", "gamma", "phi", 
               "pResponse","lambda","dmean","p0Oth",
               "p0", "sigma", "betaDens", "betaCovs","betaResponse",
               "betaCovsOth","betaResponseOth")

nimParams2 <- c("z", "sxy")



## ------   5. NIMBLE INITS ------

## ------     5.1. GENERATE z INITIAL VALUES -----

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



## ------     5.2. GENERATE detResponse INITIAL VALUES -----

InitsDetResponse <- detResponse
InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1, 0.5)
InitsDetResponse[!is.na(detResponse)] <- NA



## ------     5.3. GENERATE sxy INITIAL VALUES ------

##-- create a data.frame with all detection of all Individuals detected
##-- project death to the next year
myData.deadProj <- myData.dead[ ,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1 # Project dead reco to the next year

##-- Remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]

AllDets <- rbind(myData.alive$myData.sp[,c("Id","Year")],
                 myData.deadProj[,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  coordsData = AllDetsxy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid = T)$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)
idAugmented <- which(rownames(z) %in%"Augmented")
Id.vector <- y.ar$Id.vector


sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = nimData$lowerHabCoords,
                       upperCoords = nimData$upperHabCoords,
                       habitatGrid = nimData$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")



## ------     5.4. LIST INITS ------

for(c in 1:4){
  
  nimInits <- list( "sxy" = round(sxy.init,5),
                    "dmean" = runif(1,0,10),
                    "z" = z.init,
                    "omeg1" = c(0.5,0.5),
                    "gamma" = runif(dim(y.alive)[3]-1,0,1),
                    "p01" = array( runif(18,0,0.2), 
                                  c(nimConstants$n.counties,dim(y.alive)[3])),
                    "p01Oth" = array( runif(18,0,0.2),
                                      c(nimConstants$n.countries+1,dim(y.alive)[3])),
                    "sigma" = runif(nYears,1,4),
                    "betaDens" = runif(1,-0.1,0.1),
                    "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),
                                        c(dim(detCovsOth)[3],nYears)),
                    "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),
                                           c(dim(detCovsOth)[3],nYears)),
                    "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
                    "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
                    "detResponse" = InitsDetResponse,
                    "pResponse"  = runif(1, 0.4, 0.5),
                    "phi" = runif(dim(y.alive)[3]-1,0.1,0.3)) 
  
  ## TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  idDEtected <- which(!rownames(z) %in%"Augmented")
  for(i in 1:length(idDEtected)){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i,,t]  
      } else {
        SXY <- nimData$sxy[i,,t]
      }
      sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
                                                        1:nimConstants$maxNBDets] 
      DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
      index <- DETECTIndexdetectorIndex[sxyID,1:DETECTLESS[sxyID]]
      
      ## GET NECESSARY INFO 
      n.detectors <- length(index)
      YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
      YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
      
      ## RECREATE Y
      if(nimData$nbDetections[i,t] > 0){
        for(j in 1:nimData$nbDetections[i, t]){
          ## check if a detection is out of the "detection window"
          if(sum(YDET[j]==index) == 0){
            print(file.path("id",i,"t",t,"j",j))
          }
        }#j
      }
    }#t
  }#i
  
  # plot(nimData$detector.xy[ ,2] ~ nimData$detector.xy[ ,1])
  # points(nimData$detector.xy[index,2] ~ nimData$detector.xy[index,1], col = "red")
  # points(SXY[2] ~ SXY[1], col = "blue", pch = 16)
  # points(nimData$detector.xy[YDET[1:nimData$nbDetections[i,t]],2] ~
  #          nimData$detector.xy[YDET[1:nimData$nbDetections[i,t]],1],
  #        col = "green", pch = 16)
  # points(nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i,t]],2] ~
  #          nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i,t]],1],
  #        col = "purple", pch = 16)
  # 
  # plot(st_geometry(COUNTRIES))
  # tmp <- myData.alive$myData.sp[myData.alive$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] & myData.alive$myData.sp$Year %in% years[t], ]
  # plot(st_geometry(tmp), col = "red", add = T)
  # 
  # ##CHECK WHERE IS NORRBOTTEN. IT IS ON THE 5TH INDEX
  # plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  # plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  # plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  # plot(st_geometry(detectors$main.detector.sp[detCounties%in% c(1),]), col = myCol[5], pch = 16, cex = 0.8,add=T)

  
  
  ## ------ 7. SAVE NIMBLE INPUT ------
  save(nimData,
       nimConstants,
       y.dead,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(myVars$WD, myVars$modelName,
                        file.path( myVars$modelName, "Chain", c, ".RData")))
}#c



################################################################################
## ------ 8. SAVE NECESSARY OBJECTS ------

habitat.list <- habitat
myStudyArea.poly <- myStudyArea

save(habitat.list, detectors, COUNTRIES, myStudyArea.poly,COUNTIES_AGGREGATEDSubset,
     myFilteredData.sp, myFullData.sp, COUNTIES_AGGREGATED,
     file = file.path(myVars$WD, myVars$modelName, "NecessaryObjects.RData" ))


