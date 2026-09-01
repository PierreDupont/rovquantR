##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE OPSCR analysis 2025 
##
## This R script reproduces the OPSCR analysis of the wolverine data as performed
## by RovQuant in 2024.
##
## NOTES : Two main updates from last year's script ('53.Cleaned2024.R')
##  1. Added a step to remove duplicated GPS tracks (from Asun)
##  2. Fixed the filtering of Rovbase samples for the detection covariates
##
## For later: add dead recovery states ('recovered dead legal' and 'recovered dead other') 
## as in the last wolf ('40.F_2024_sf.R') and bear analyses ('Bear_NORWAY_2015-2024.R').
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 20/10/2025
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2025
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##   
##------------------------------------------------------------------------------

rm(list=ls())
gc()

## ------ IMPORT REQUIRED LIBRARIES ------

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
library(dplyr)
library(readxl)
library(spatstat)
library(ggplot2)
library(stars)

##-- Identify user and set corresponding DropBox and Git/Rovquant directory
if(Sys.info()['user'] == 'pidu') {
  dir.git <- "C:/My_documents/RovQuant"
  dir.dropbox <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'pierredupont') {
  dir.git <- "/Users/pierredupont/Documents/RovQuant"
  dir.dropbox <- "/Users/pierredupont/Dropbox (AQEG)/AQEG Team Folder/RovQuant/"
}
if(Sys.info()['user'] == 'cymi') {
  dir.git <- "C:/My_documents/rovquant/analyses/Rgit/RovQuant/"
  dir.dropbox <- "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'richbi') {
  dir.git <- "C:/Users/richbi/OneDrive - Norwegian University of Life Sciences/PROJECTS/RovQuant"
  dir.dropbox <- "C:/Users/richbi/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'seasunci') {
  dir.git <- "C:/Users/seasunci/02_RovQuant/RovQuant"
  dir.dropbox <- "C:/Users/seasunci/AQEG Dropbox/AQEG Team Folder/RovQuant"
}


## ------ SET REQUIRED WORKING DIRECTORIES ------

# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/myWorkingDirectories.R")
source(file.path(dir.git,"Temp/PD/myWorkingDirectories.R"))
# source(file.path(dir.git,"Temp/ASP/myWorkingDirectories.R"))


## ------ SOURCE THE REQUIRED FUNCTIONS ------

sourceDirectory(dir.function, modifiedOnly = FALSE)
sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)
load(file.path(dir.dropbox,"DATA/MISC DATA/age.lookup.table.RData"))


## ------ SOURCE THE NIMBLE FUNCTION ------

# source("C:/My_documents/RovQuant/Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R")
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R")
source(file.path(dir.git, "Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R"))


##------------------------------------------------------------------------------

## ------ 0.SET ANALYSIS CHARACTERISTICS -----

myVars <- list(
  ## WORKING DIRECTORY & MODEL NAME
  WD = file.path(dir.dropbox,"wolverine/CM/2025"),
  modelName = "54.Cleaned2025",
  
  ## HABITAT SPECIFICATIONS
  HABITAT = list( countries =  c("SWE","NOR"),
                  habResolution = 20000,
                  habBuffer = 60000),
  
  ## NGS DATA SPECIFICATIONS
  DATA = list( years = 2015:2024, #2014:2023
               species = c("Jerv"),              
               sex = c("Hann","Hunn"),                   
               samplingMonths = list(12,1:6)),   
  
  ## DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 2000,
                    detResolution = 10000,
                    detDeadResolution = 15000),
  
  ## DATA GENERATION
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 1,
                     aug.factor = 0.8),
  
  ## OUTPUT PLOTS
  OUTPUT = list(mapResolution = 10000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(myVars$WD, myVars$modelName))){
  dir.create(file.path(myVars$WD, myVars$modelName))
  dir.create(file.path(myVars$WD, myVars$modelName, "Hunn"))
  dir.create(file.path(myVars$WD, myVars$modelName, "Hann"))
}



##------------------------------------------------------------------------------

## ------ I.LOAD AND SELECT DATA ------

## ------ 1. HABITAT DATA ------

## ------    1.1. LOAD RAW SHAPEFILES ------

##-- POLYGONS OF THE REGION
GLOBALMAP <- file.path( dir.dropbox,
                        "DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp") %>%
  st_read(.) %>%
  filter(area > 80000000) %>%
  st_crop(., st_bbox(extent(c(-70000,1200000,5100000,8080000))))

##-- POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP %>%
  filter(ISO %in% c("SWE","NOR")) %>%
  group_by(ISO) %>%
  summarize()

##-- POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp"))   
COMMUNES_SWE <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp"))   
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)

##-- POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- COMMUNES %>%
  group_by(NAME_1) %>%
  summarize()

##-- AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATE <- COUNTIES
COUNTIES_AGGREGATE$id <- 1:nrow(COUNTIES_AGGREGATE)
#[CM] adjust Counties aggregation
COUNTIES_AGGREGATE$id[c(24,3,15,9,14,38,40,21,27,37,31,26,34,5,8,12,36,13,7)] <- 3
COUNTIES_AGGREGATE$id[c(39,33,23,32,29,22,4,11,20,2,10,16,25,1)] <- 4
COUNTIES_AGGREGATE$id[c(19)] <- 1
COUNTIES_AGGREGATE$id[c(35)] <- 2
COUNTIES_AGGREGATE$id[c(17,28)] <- 5
COUNTIES_AGGREGATE$id[c(18)] <- 7
COUNTIES_AGGREGATE$id[c(30)] <- 8
COUNTIES_AGGREGATE <- COUNTIES_AGGREGATE %>% group_by(id) %>% summarize()
COUNTIES_AGGREGATED <- st_simplify( COUNTIES_AGGREGATE,
                                    preserveTopology = TRUE,
                                    dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id
# ggplot(COUNTIES_AGGREGATED) +
#   geom_sf(aes(fill = id)) +
#   geom_sf_label(aes(label = id))



## ------    1.2. CREATE STUDY AREA POLYGON ------

##-- CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
if(!is.null(myVars$HABITAT$countries)){
  myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% myVars$HABITAT$countries, ]
  
  ##-- CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
  myBufferedArea <- st_buffer(x = st_as_sf(myStudyArea),
                              dist = myVars$HABITAT$habBuffer) %>%
    mutate(id = 1) %>%
    group_by(id) %>% 
    summarize() %>%
    st_intersection(., GLOBALMAP)
}

myStudyArea$id <- myStudyArea %>%
  mutate(id = 1) %>%
  group_by(id) %>% 
  summarize()

##-- PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myBufferedArea), add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
  plot(st_geometry(myStudyArea), add = TRUE, col = "red")
}



## ------ 2. NGS DATA ------

## ------    2.1. LOAD ROVBASE FILES ------

## NGS data from RovBase
DNA <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dna_wolverines.csv"),
                 fileEncoding = "latin1") 

## Dead Recoveries from RovBase
DEAD <- read.csv( file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dead_carnivores.csv"),
                  fileEncoding = "latin1") 

## DNA samples to be removed from Henrik (no samples to be removed from Sweden)
SUSPECT_NGS_SAMPLES <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/Remove ngs samples list wolverine 2024.csv"),
                                 fileEncoding = "latin1") 
SUSPECT_NGS_SAMPLES_2025 <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/Additional samples to add to the Remove ngs samples list wolverine from Norway 2025.csv"),
                                      fileEncoding = "latin1")
SUSPECT_NGS_SAMPLES <- rbind(SUSPECT_NGS_SAMPLES, SUSPECT_NGS_SAMPLES_2025)  # rbind file from last year with the one from 2025

## DEAD RECOVERIES to be removed from Henrik  --> Only file from 2024, NO ADDITIONAL INDIVIDUALS IN 2025.
SUSPECT_DeadRecoSAMPLES <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/Remove dead recoveries list wolverine 2024.csv"),
                                     fileEncoding = "latin1") 

## DNA samples to be removed from Henrik
HairTrapSamples <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/hairtrapsNB2024.xlsx")) 
HairTrapSamples_2025 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/hairtrapsEvaHedmark2025.xlsx")) 
HairTrapSamples_2025$Funnetdato <- as.POSIXct(strptime(HairTrapSamples_2025$Funnetdato, "%Y-%m-%d"))
HairTrapSamples <- rbind(HairTrapSamples, HairTrapSamples_2025)  # rbind file from last year with the one from 2025 (information from Eva Hedmark)
#### !!!!!! QUESTION ---> HairTrapSamples are actually these remove later?
####                      We should also add "HairTrapSamples_2025" in the other observations covariate????


## Wolverine den locations
DEN <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/DEN_COUNTS_2009_2025_fromHB.csv"),
                 fileEncoding = "latin1")

## Skandobs observations
skandObs <- read_xlsx(file.path(dir.dropbox, "DATA/Skandobs/RB_Skandobs_2012_2025/Richard_Bischof_Skandobs_2012_20250915.xlsx"))

## Rovbase observations
rovbaseObs1 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dna_all_batch1.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dna_all_batch2.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dna_all_batch3.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3)
rm(list = c("rovbaseObs1", "rovbaseObs2", "rovbaseObs3"))



## ------    2.2. TRANSLATE SCANDINAVIAN CHARACTERS ------

## Drop a column that makes cleanDataNew to fail
colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA), dir.translation = dir.analysis )
DNA <- DNA[,-which(colnames(DNA)%in% "Kjoenn..Individ.")]

colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD), dir.translation = dir.analysis )
colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN), dir.translation = dir.analysis)
colnames(skandObs) <- translateForeignCharacters(dat=colnames(skandObs), dir.translation = dir.analysis )
colnames(rovbaseObs) <- translateForeignCharacters(dat=colnames(rovbaseObs), dir.translation = dir.analysis )
rovbaseObs$Proevetype <- translateForeignCharacters(dat=rovbaseObs$Proevetype, dir.translation = dir.analysis )



## ------ 3. SEARCH EFFORT DATA ------

## ------    3.1. GPS SEARCH TRACKS ------

## LOAD GPS SEARCH TRACKS
TRACKS_SINGLE <- read_sf(file.path(dir.dropbox,
                                   "DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp", sep = ""))
TRACKS_MULTI <- read_sf(file.path(dir.dropbox,
                                  "DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp", sep = ""))

## COMBINE ALL TRACKS AND FIX DATES
ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI) %>%
  mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
          Yr = as.numeric(format(Dato,"%Y")),
          Mth = as.numeric(format(Dato,"%m")),
          Dato = as.character(Dato)) %>%
  dplyr::filter(
    ## REMOVE HELICOPTER TRACKS
    Helikopter == "0",
    ## KEEP ONLY WOLVERINE TRACKS
    Jerv == "1")

## GET EXTENT
myStudyArea.extent <- st_bbox(extent(myStudyArea))
st_crs(myStudyArea.extent) <- st_crs(COUNTRIES)

## [PD] ASP ADDED REMOVAL OF DUPLICATED TRACKS
dupIDs <- dupDist <- TRACKS_YEAR <- list()
for(t in 1:nYears){
  
  TRACKS <- ALL_TRACKS %>% 
    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
    filter( Yr%in%YEARS[[t]][1] & Mth %in% myVars$DATA$samplingMonths[[1]] |
              Yr%in%YEARS[[t]][2] & Mth%in%myVars$DATA$samplingMonths[[2]]) %>%
    ## SUBSET TRACKS TO THE STUDY AREA
    st_intersection(., st_as_sfc(myStudyArea.extent)) 
  
  ## NAME TRACKS
  TRACKS$ID <- 1:nrow(TRACKS)
  ## CALCULATE LENGTH OF EACH TRACK TO IDENTIFY DUPLICATES
  TRACKS$dist <- st_length(TRACKS, byid = T)
  ## CALCULATE CENTROIDS TO AVOID KEEPING TRACKS WITH THE SAME LENGTHS BUT IN DIFFERENT LOCATIONS
  TRACKS$centroidx <- st_coordinates(st_centroid(TRACKS))[ ,1]
  
  ## FIND DUPLICATES BASED ON PERSON, DISTANCE and DATE
  df <- data.frame( Dato = TRACKS$Dato,
                    Person = TRACKS$Person,
                    dist = TRACKS$dist,
                    centroidx = TRACKS$centroidx)
  dupIDs[[t]] <- TRACKS$ID[duplicated(df)]
  dupDist[[t]] <- TRACKS$dist[duplicated(df)]
  
  ## STORE CLEAN TRACKS IN A LIST
  TRACKS_YEAR[[t]] <- TRACKS[-dupIDs[[t]], ]
}#t



## ------    3.2. DISTANCE TO ROADS ------

## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(file.path( dir.dropbox,
                                  "DATA/GISData/Roads/MinDistAllRoads1km.tif"))

## RASTERIZE DISTANCE TO ROADS
r <- fasterize(st_as_sf(myStudyArea), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea)



## ------    3.3. DAYS OF SNOW ------

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
SNOW <- stack(file.path( dir.dropbox, 
                         "DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2014_2025_Wolverine.tif")) # Average snow from December to June (the official monitoring period for Norway&Sweden)
## RENAME THE LAYERS
names(SNOW) <- paste(2014:2024,(2014:2024)+1, sep="_")

## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))



## ------    3.4. LOAD SCANDINAVIAN 20KM HABITAT ------

load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/Habitat20km.RData"))
load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData"))



##------------------------------------------------------------------------------

## ------ II.CREATE SCR DATA ------

## ------ 1. CLEAN & FILTER NGS DATA ------

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]

## Remove un-verified dead recoveries [HB]
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]


## ------    1.1. CLEAN NGS & DEAD RECOVERY DATA ------

myCleanedData.sp <- CleanDataNew2sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = myVars$DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(myVars$DATA$samplingMonths)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)



## ------    1.2. FILTER DATA ------

myFullData.sp <- FilterDatasf(
  myData = myCleanedData.sp,
  poly = myStudyArea,
  dead.recovery = T ,
  sex = c("Hann","Hunn"), 
  setSex = T)


## REMOVE SUSPECT SAMPLES ACCORDING TO HENRIK (THERE ARE NOT SUSPECT SAMPLES IN SWEDEN)
myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
dim(myFullData.sp$alive)

##EXPORT THE DATA
myFullData.spSaved <- myFullData.sp


## REMOVE SUSPECT DEAD RECOVERIES ACCORDING TO HENRIK
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]

## REMOVE CUBS OF THE YEAR ACCORDING TO EVA HEDMARK (email 23/11/2025, document with explanation in AQEG Dropbox\AQEG Team Folder\RovQuant\DATA\RovbaseData\ROVBASE DOWNLOAD 20251121)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% c("D608410", "D608411", "D605997")),]

## Remove individuals that died twice
myFullData.sp$dead.recovery$Id <- as.character(myFullData.sp$dead.recovery$Id)
IdDoubleDead <- myFullData.sp$dead.recovery$Id[duplicated(myFullData.sp$dead.recovery$Id)]

if(length(IdDoubleDead) > 0){
  duplicatedDeath <- NULL
  for(i in IdDoubleDead){
    tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
    if(length(tmp)==0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-2]}##[CM] remove the second record.
    duplicatedDeath <- c(duplicatedDeath, tmp)
  }#i
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]
}#if

## Dead recoveries flagged by Henrik that should always be removed (email from the 18/12/2024)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$RovBaseId %in% c("M495994", "M524051", "M524052", "M524053"), ]

## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) Remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
sum(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
      myFullData.sp$dead.recovery$Month > 2 &
      myFullData.sp$dead.recovery$Month < 12)

myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
                                                                    myFullData.sp$dead.recovery$Month > 2 &
                                                                    myFullData.sp$dead.recovery$Month < 12),]


## 2) Remove individuals that have a weight >0 and <4 between March and November format the weight correctly
myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)

##-- Convert to decimals
myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))

##-- Get the two weight columns together
myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt),
                                             myFullData.sp$dead.recovery$Helvekt,
                                             myFullData.sp$dead.recovery$Slaktevekt)

##-- Assign negative values to nas to avoid issues
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999

##-- check how many dead reco we remove and remove if more than 0
if(sum(myFullData.sp$dead.recovery$weight > 0 &
       myFullData.sp$dead.recovery$weight < 4 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
                                                                      myFullData.sp$dead.recovery$weight < 4 &
                                                                      myFullData.sp$dead.recovery$Month < 12 &
                                                                      myFullData.sp$dead.recovery$Month > 2),]
}

##-- check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
                                myFullData.sp$dead.recovery$Month < 12 &
                                myFullData.sp$dead.recovery$Month > 2,  ]
}



## ------    1.3. FILTER NGS & DEAD RECOVERY DATA ------

myFilteredData.sp <- myFullData.sp

## Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]

myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]

## Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in% unlist(myVars$DATA$samplingMonths), ]



## ------    1.4. SUBSET DETECTIONS IN NORRBOTTEN IN ALL YEARS EXCEPT 2017, 2018 and 2019 ------ 

COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018,2023,2024)
is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))

## Remove samples in Norrbotten nin years not sampled
myFilteredData.sp$alive <- myFilteredData.sp$alive[- which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                                             !is.na(is.Norr)), ]



## ------    1.5. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

## ------      1.5.1. ASSIGN SAMPLES TO TRACKS ------

## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
myFilteredData.sp$alive$TrackRovbsID <- NA
myFilteredData.sp$alive$TrackDist <- NA
TRACKSSimple_sf <- list()
for(t in 1:nYears){+
    TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
    TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbsID)
    TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]
}

## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
dnatemp <- st_as_sf(myFilteredData.sp$alive)

## CREATE A BUFFER AROUND EACH DETECTION
# tmp <-  st_buffer(dnatemp, dist = 750)
# for(i in 1:nrow(myFilteredData.sp$alive)){
#    ## MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
#    t <- which(years %in% tmp[i,]$Year)
#    whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
#    
#    ## INTERSECT POINT WITH TRACKS
#    tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i, ])
# 
#    ## If not TRACKS that date, move on
#    if(nrow(tmpTRACKS)==0){next}
# 
#    ## Else, find the closest TRACK
#    dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
# 
#    ## IF NO MATCHING DATE TRACK ASSIGN NA
#    if(length(dist)==0){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- NA
#      myFilteredData.sp$alive$TrackDist[i] <- NA
#    }
#    ## IF 1 MATCHING TRACK ASSIGN TO THAT TRACK
#    if(length(dist)==1){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
#      myFilteredData.sp$alive$TrackDist[i] <- dist
#    }
#    ## IF SEVERAL MATCHING DATES ASSIGN TO THE CLOSEST OF THE MATCHING TRACKS
#    if(length(dist)>1){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
#      myFilteredData.sp$alive$TrackDist[i] <- min(dist)
#    }
#    print(i)
#    #if(is.na(myFilteredData.sp$alive$TrackRovbsID[i])){print(i)}
# }#i
#
# ## SAVE FOR FASTER LOADING
# save( myFilteredData.sp,
#       file = file.path(myVars$WD, myVars$modelName, "myFilteredData_original.RData"))
load(file.path(myVars$WD, myVars$modelName, "myFilteredData_original.RData"))



## ------      1.5.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC & STRUCTURED ------

distanceThreshold <- 500

## Proeveleverandoer columns was replaced by two columns, merging them now...
myFilteredData.sp$alive$Proeveleverandoer <- ifelse(myFilteredData.sp$alive$Annen.innsamler...Rolle %in% "", 
                                                    myFilteredData.sp$alive$Samlet.selv...Rolle,
                                                    myFilteredData.sp$alive$Annen.innsamler...Rolle)

## Structured samples are:
##    1. samples from agencies ("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen")  
##    2. samples assigned to a GPS search track
##    3. wihtin distanceThreshold of a GPS search track
whichStructured <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold

myFilteredData.spStructured <- myFilteredData.sp$alive[whichStructured, ]
myFilteredData.spOthers <- myFilteredData.sp$alive[!whichStructured, ]



## ------    1.6. SEPARATE MORTALITY CAUSES ------ 

## DEFINE LEGAL MORTALITY
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])

## SPLIT MORTALITY CAUSES
legal.death <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]
Other.death <- myFilteredData.sp$dead.recovery[!myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]



## ------ 2. GENERATE HABITAT ------

## ------    2.1. REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ------

## DELINEATE A BUFFER AROUND ALL DETECTIONS 
myBufferedArea <- st_buffer( myFilteredData.sp$alive,
                             dist = myVars$HABITAT$habBuffer * 1.4) %>%
  mutate(id = 1) %>%
  group_by(id) %>%
  summarize()

## CUT TO SWEDISH AND NORWEGIAN BORDERS
myStudyArea <- st_intersection(myBufferedArea, myStudyArea)



## ------    2.2. GENERATE HABITAT CHARACTERISTICS FROM THE NEW HABITAT DEFINITION ------

myHabitat <- MakeHabitatFromRastersf( 
  poly = myStudyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = myVars$HABITAT$habBuffer,                               
  plot.check = T)

## RETRIEVE HABITAT WINDOWS BOUNDARIES
lowerHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1, ] - 0.5*myVars$HABITAT$habResolution
upperHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1, ] + 0.5*myVars$HABITAT$habResolution
nHabCells <- dim(lowerHabCoords)[1]

## CREATE HABITAT GRID 
habIDCells.mx <- myHabitat$IDCells.mx 
habIDCells.mx[] <- 0
scaledHabGridCenters <- scaleCoordsToHabitatGrid(
  coordsData = myHabitat$habitat.xy,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid =F )$coordsHabitatGridCenterScaled

scaledHabGridCenters <- scaledHabGridCenters[myHabitat$habitat.r[] == 1, ]
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}



## ------    2.3. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 

## Remove samples outside the STUDY AREA #[CM]
myStudyArea$idd <- 1
myStudyAreaAggregated <- myStudyArea %>% group_by(idd) %>% summarize()
whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyAreaAggregated))))
if(length(whichOut)>0){
  myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
}
myFilteredData.sp$alive$Id <- droplevels( myFilteredData.sp$alive$Id)

## REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT #[CM] 
whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, myHabitat$buffered.habitat.poly))))
if(length(whichOutBuff)>0){
  myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]
}

## check correlation number of detections ~ between monitoring season
myFilteredData.sp$dead.recovery$Id <- as.character(myFilteredData.sp$dead.recovery$Id)
myFilteredData.sp$alive$Id <- as.character(myFilteredData.sp$alive$Id)



## ------    2.4. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       2.4.1. DEN COUNTS ------

DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X","UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id  <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[ ,"id"]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD( as(DEN.sp,"Spatial"),
              h = 30000,
              grid = as(myHabitat$habitat.r, 'SpatialPixels'))))

## EXTRACT COVARIATE
denCounts <- DEN.r[myHabitat$habitat.r[ ] == 1]
denCounts <- round(scale(denCounts), digits = 2)



## ------ 3. GENERATE DETECTORS ------

## ------    3.1. GENERATE DETECTORS CHARACTERISTICS ------

## GENERATE SUB-DETECTORS BASED ON THE STUDY AREA
habitat.subdetectors <- disaggregate(
  myHabitat$habitat.rWthBuffer,
  fact = res(myHabitat$habitat.r)[1]/myVars$DETECTORS$detSubResolution)

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
myDetectors <- myDetectors.dead <- MakeSearchGridsf(
  data = habitat.subdetectors,
  resolution = myVars$DETECTORS$detResolution,
  div = (myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)^2,
  plot = FALSE,
  fasterize = TRUE)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(myDetectors$main.detector.sp)[1]
n.detectors.dead <- dim(myDetectors.dead$main.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(myDetectors$main.detector.sp)
n.trials <- as.vector(table(myDetectors$detector.sp$main.cell.id))
detector.dead.xy <- st_coordinates(myDetectors.dead$main.detector.sp)

## IDENTIFY DETECTORS IN NORBOTTEN 
COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten","Troms","Västerbotten",
                                                            "Nordland","Finnmark"),]
COUNTIESAroundNorrbotten <- st_simplify(COUNTIESAroundNorrbotten, dTolerance = 500)

## CREATE A NORRBOTTEN DETECTOR GRID
distDestsCounties <- st_distance(myDetectors$main.detector.sp, COUNTIESAroundNorrbotten,byid = T)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min)==3)

## RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * myVars$DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * myVars$DETECTORS$detResolution



## ------    3.2. GENERATE DETECTOR-LEVEL COVARIATES ------

## ------      3.2.1. EXTRACT COUNTRIES ------

dist <- st_distance(myDetectors$main.detector.sp, COUNTRIES, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))



## ------      3.2.2. EXTRACT COUNTIES ------

dist <- st_distance(myDetectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F )
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties),]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))



## ------      3.2.3. EXTRACT GPS TRACKS LENGTHS ------

## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp),
                                      rep(1,nrow(myDetectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)

## INITIALIZE MATRIX & RASTERS OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
TRACKS.r <- list()

## CALCULATE THE LENGTH OF THE TRACKS
for(t in 1:nYears){
  intersection <- st_intersection(detectorGrid, TRACKS_YEAR[[t]]) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN))   
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[ ,t]
  print(t)
}#t



## ------      3.2.4. EXTRACT DISTANCES TO ROADS ------

## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)

## if NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( DistAllRoads,
                        myDetectors$main.detector.sp[isna, ],
                        buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp



## ------      3.2.5. EXTRACT DAYS OF SNOW ------

## EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(myDetectors$main.detector.sp)[1], ncol = nYears)
det.sptransf <- st_transform(myDetectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp



## ------      3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------        3.2.6.1. SKANDOBS ------

## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))

## MAKE IT SPATIAL 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(myStudyArea))

## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist(myVars$DATA$samplingMonths))
skandObs$monitoring.season <- ifelse(skandObs$month > 12, skandObs$year, skandObs$year-1)
skandObs <- skandObs[subset, ] 

## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(myHabitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1,]

subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace,] 

## RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate(habitat.subdetectors, fact=(myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution))
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(skandObs[skandObs$monitoring.season %in% y,1], r.detector, fun = "count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})
r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))



## ------        3.2.6.2. ROVBASE ------

## GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`), ]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

## DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)

## SUBSET THE DATA 
filter <- list( 
  species = "Jerv",
  type = c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
            "Saliv/Spytt", "Loepeblod", "Blod", "Vev"), ## !!!! ASP : Blod as well
  month = unlist(myVars$DATA$samplingMonths))

## SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month > 12, rovbaseObs.sp$year, rovbaseObs.sp$year-1) 
rovbaseObs.sp <- rovbaseObs.sp[subset, ] 

## REMOVE SAMPLES THAT WERE SUCCESSFULLY GENOTYPED & FROM THE FOCAL SPECIES 
subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$Individ) 
rovbaseObs.sp <- rovbaseObs.sp[!subset, ] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace, ] 

## RASTERIZE 
r.detector <- aggregate(habitat.subdetectors, fact=(myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution))
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(rovbaseObs.sp[rovbaseObs.sp$monitoring.season %in% y, 1], r.detector , fun="count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})

r.OtherSamplesBinary <- brick(lapply(r.list, function(x) x[[1]]))
r.OtherSamplesContinuous <- brick(lapply(r.list, function(x) x[[2]]))



## ------        3.2.6.3. COMBINE ROVBASE & SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][] > 1] <- 1
}



## ------        3.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------

## IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.sp$alive[which(myFilteredData.sp$alive$DNAID %in% HairTrapSamples$DNAID), ]

## MANUALLY FIND THE HAIR SMAPLES & COLOR THE CELL
tmpyr <- unique(tmpHair$Year)
for(i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr[i])
  whereHair <- raster::extract( r.SkandObsOtherSamplesBinary[[t]],
                                tmpHair[tmpHair$Year %in% tmpyr[i], ],
                                cellnumbers = T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[ ,1]] <- 1
}



## ------        3.2.6.7. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract( r.SkandObsOtherSamplesBinary,
                                                myDetectors$main.detector.sp)



## ------      3.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],2))
detCovs[ , ,1] <- detTracks
detCovs[ , ,2] <- detSnow

detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],3))
detCovsOth[ , ,1] <- detSnow
detCovsOth[ , ,2] <- matrix(detRoads, length(detRoads), nYears)
detCovsOth[ , ,3] <- detOtherSamples



## ------ 4. GENERATE y DETECTION ARRAYS ------

## ------    4.1. ASSIGN SAMPLES TO DETECTORS ------

## ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$alive,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = myVars$DETECTORS$detResolution)

## STRUCTURED
myData.aliveStruc <- AssignDetectors_v3sf(
  myData = myFilteredData.spStructured,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = myVars$DETECTORS$detResolution)

## OTHERS
myData.aliveOthers <- AssignDetectors_v3sf(
  myData = myFilteredData.spOthers,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = myVars$DETECTORS$detResolution)

## DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf(
  myData = myFilteredData.sp$dead.recovery,
  myDetectors = myDetectors.dead$main.detector.sp,
  radius = myVars$DETECTORS$detResolution)


##-- [PD] fixed assignment to sub-detectors in Norrbotten
##-- Identify sub-detectors inside Norrbotten
subDetsNorrbotten <- which( myDetectors$detector.sp$main.cell.id %in% 
                              myDetectors$main.detector.sp$main.cell.id[detsNorrbotten])


## ALL
##-- Loop over flagged detections and assign them to closest sub-detector outside Norrbotten
whichdets <- which(!myData.alive$myData.sp$Year %in% yearsSampledNorrb &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)
for(i in 1:length(whichdets)){
  tmp <- myData.alive$myData.sp[whichdets[i], ]
  ## Calculate distance to all sub-detectors
  dist <- st_distance(tmp, myDetectors$detector.sp)
  ## Artificially increase distance for detectors in Norrbotten 
  dist[ ,subDetsNorrbotten] <- 500000
  ## Assign detection to closest sub-detector outside Norrbotten
  myData.alive$myData.sp$sub.detector[whichdets[i]] <- which.min(dist[1, ])
  ## Assign detection to the corresponding main detector outside Norrbotten
  thisDet <- myDetectors$detector.sp$main.cell.id[which.min(dist[1, ])]
  myData.alive$myData.sp$Detector[whichdets[i]] <- which(myDetectors$main.detector.sp$main.cell.id == thisDet)
}#i


## STRUCTURED
##-- Loop over flagged detections and assign them to closest sub-detector outside Norrbotten
whichdetsStruc <- which(!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveStruc$myData.sp$Detector %in% detsNorrbotten)
for(i in 1:length(whichdetsStruc)){
  tmp <- myData.aliveStruc$myData.sp[whichdetsStruc[i], ]
  ## Calculate distance to all sub-detectors
  dist <- st_distance(tmp, myDetectors$detector.sp)
  ## Artificially increase distance for detectors in Norrbotten 
  dist[ ,subDetsNorrbotten] <- 500000
  ## Assign detection to closest sub-detector outside Norrbotten
  myData.aliveStruc$myData.sp$sub.detector[whichdetsStruc[i]] <- which.min(dist[1, ])
  ## Assign detection to the corresponding main detector outside Norrbotten
  thisDet <- myDetectors$detector.sp$main.cell.id[which.min(dist[1, ])]
  myData.aliveStruc$myData.sp$Detector[whichdetsStruc[i]] <-  which(myDetectors$main.detector.sp$main.cell.id == thisDet)#which(myDetectors$main.detector.sp$main.cell.id == thisDet)
}#i


## OTHER
##-- Loop over flagged detections and assign them to closest sub-detector outside Norrbotten
whichdetsOther <- which(!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveOthers$myData.sp$Detector %in% detsNorrbotten)
for(i in 1:length(whichdetsOther)){
  tmp <- myData.aliveOthers$myData.sp[whichdetsOther[i], ]
  ## Calculate distance to all sub-detectors
  dist <- st_distance(tmp, myDetectors$detector.sp)
  ## Artificially increase distance for detectors in Norrbotten 
  dist[ ,subDetsNorrbotten] <- 500000
  ## Assign detection to closest sub-detector outside Norrbotten
  myData.aliveOthers$myData.sp$sub.detector[whichdetsOther[i]] <- which.min(dist[1, ])
  ## Assign detection to the corresponding main detector outside Norrbotten
  thisDet <- myDetectors$detector.sp$main.cell.id[which.min(dist[1, ])]
  myData.aliveOthers$myData.sp$Detector[whichdetsOther[i]] <-  which(myDetectors$main.detector.sp$main.cell.id == thisDet)#which(myDetectors$main.detector.sp$main.cell.id == thisDet)
}#i



## ------    4.2. SAVE PREPARED DATA ------

# save( myData.alive, myData.aliveStruc, myData.aliveOthers, myData.dead,
#       file = file.path(myVars$WD, myVars$modelName, "myFilteredData.RData"))

## [CM] NEEDS THIS OTHERWISE THE LOOPS OVERWRITES THE SEX
myData.aliveALL <- myData.alive
myData.aliveStrucALL <- myData.aliveStruc
myData.aliveOthersALL <- myData.aliveOthers
myData.deadALL <- myData.dead



## ------    4.3. GENERATE DETECTION HISTORY (FOR BOTH SEXES) ------

for(thisSex in c("Hann","Hunn")){
  
  message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
  
  ## ------    4.4. FILTER DATA BY SEX ------
  
  # load(file.path(myVars$WD, myVars$modelName, "myFilteredData.RData"))
  
  myData.alive$myData.sp <- myData.aliveALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  
  myData.aliveStruc$myData.sp <- myData.aliveStrucALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  
  myData.aliveOthers$myData.sp <- myData.aliveOthersALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  
  myData.dead <- myData.deadALL %>%
    dplyr::filter(Sex %in% thisSex)
  
  # do this for the dead recoveries
  myFullData.spDeadsex <- myFullData.spSaved$dead.recovery %>%
    dplyr::filter(Sex %in% thisSex)
  
  
  
  ## ------    4.5. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ------
  
  ## ALL SAMPLES
  y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                   myDetectors = myDetectors$main.detector.sp,
                   method = "Binomial",
                   myData2 = myData.dead,
                   myDetectors2 = myDetectors.dead$main.detector.sp,
                   returnIdvector = TRUE)
  y.ar.ALIVE <- y.ar$y.ar
  dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)
  
  ## STRUCTURED
  y.arStruc <- MakeYsf( myData = myData.aliveStruc$myData.sp,
                        myDetectors = myDetectors$main.detector.sp,
                        method = "Binomial",
                        myData2 = myData.dead,
                        myDetectors2 = myDetectors.dead$main.detector.sp,
                        returnIdvector = TRUE)
  y.ar.ALIVEStruc <- y.arStruc$y.ar
  dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)
  
  ## OTHERS
  y.arOth <- MakeYsf( myData = myData.aliveOthers$myData.sp,
                      myDetectors = myDetectors$main.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      myDetectors2 = myDetectors.dead$main.detector.sp,
                      returnIdvector = TRUE)
  y.ar.ALIVEOth <- y.arOth$y.ar
  dimnames(y.ar.ALIVEOth) <- dimnames(y.arOth$y.ar)
  
  ## RESIZE DETECTION ARRAYS TO MAKE SURE THEY HAVE THE SAME DIMENSIONS
  y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
  y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0
  
  ## FILL IN THE Y ARRAYS 
  y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]], , ] <- y.ar.ALIVEOth
  y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]], , ] <- y.ar.ALIVEStruc
  
  ## PROJECT THE DEATH TO THE NEXT OCCASION
  y.ar.DEADProjected <- y.ar$y.ar2 
  y.ar.DEADProjected[] <- 0
  for(t in 2:nYears){y.ar.DEADProjected[ , ,t] <- y.ar$y.ar2[ , ,t-1]}
  
  ## TURN INTO ANNUAL DEAD RECOVERY MATRIX
  y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
  y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
  y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
  dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
  y.ar.DEAD[y.ar.DEAD > 0] <- 1
  
  
  
  ## ------    4.6. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------
  
  distances <- list()
  for(t in 1:nYears){
    print(paste0("------ ", t ," -------"))
    distances[[t]] <- CheckDistanceDetectionsV2sf( 
      y = y.ar.ALIVE[ , ,t], 
      detector.xy = detector.xy, 
      max.distance = myVars$DETECTIONS$maxDetDist,
      method = "pairwise",
      plot.check = F)
    
    ## REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
    y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
    
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
    }
    idd <- names(affected.ids)
    for(i in 1:length(idd)){
      detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
      myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                           myData.alive$myData.sp$Detector %in% detIds &
                                                           myData.alive$myData.sp$Year %in% years[t]), ]
    }#i
  }#t
  
  
  
  ## ------    4.7. GENERATE INDIVIDUAL-LEVEL COVARIATES ------
  
  ## ------      4.7.1. TRAP-RESPONSE ------
  
  ## Make matrix of previous capture indicator
  already.detected <- MakeTrapResponseCovsf(
    data = myFullData.sp$alive,
    data.dead = myFullData.sp$dead.recovery)
  
  ## Subset to focal years
  already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]
  
  ## Subset to focal individuals
  already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]
  
  
  
  ## ------    4.8. MAKE AUGMENTATION ------
  
  ## DATA ARRAYS
  y.alive <- MakeAugmentation(
    y = y.ar.ALIVE,
    aug.factor = myVars$DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.aliveStructured <- MakeAugmentation(
    y = y.ar.ALIVEStructured,
    aug.factor = myVars$DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.aliveOthers <- MakeAugmentation(
    y = y.ar.ALIVEOthers,
    aug.factor = myVars$DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.dead <- MakeAugmentation( 
    y = y.ar.DEAD,
    aug.factor = myVars$DETECTIONS$aug.factor,
    replace.value = 0)
  
  ## INDIVIDUAL COVARIATES
  already.detected <- MakeAugmentation( 
    y = already.detected,
    aug.factor = myVars$DETECTIONS$aug.factor,
    replace.value = 0)
  
  
  
  ##------------------------------------------------------------------------------
  
  ## ------ III.MODEL SETTING & RUNNING ------- 
  
  ## ------ 1. NIMBLE MODEL DEFINITION ------
  
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
      sxy[i, 1:2, 1] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
        upperCoords = upperHabCoords[1:numHabWindows, 1:2],
        logIntensities = logHabIntensity[1:numHabWindows],
        logSumIntensity = logSumHabIntensity,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows =  y.max,
        numGridCols = x.max)
      for(t in 2:n.years){
        sxy[i, 1:2, t] ~ dbernppACmovement_exp(
          lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
          upperCoords = upperHabCoords[1:numHabWindows, 1:2],
          s = sxy[i, 1:2, t-1],
          lambda = lambda,
          baseIntensities = habIntensity[1:numHabWindows],
          habitatGrid =  habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max,
          numWindows= numHabWindows)
      }#t 
    }#i
    
    
    
    ##----- DEMOGRAPHIC PROCESS -----## 
    
    omeg1[1:2] ~ ddirch(alpha[1:2])   
    
    for(t in 1:n.years1){
      # PRIORS 
      gamma[t] ~ dunif(0,1)
      phi[t] ~ dunif(0,1)
      
      # "UNBORN"
      omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
      # "Alive"
      omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
      # "Dead"
      omega[3,1:3,t] <- c(0,0,1)
    }#t
    
    pResponse ~ dunif(0, 1)
    
    for(i in 1:n.individuals){ 
      detResponse[i,1] ~ dbern(pResponse)
      z[i,1] ~ dcat(omeg1[1:2]) 
      for(t in 1:n.years1){
        z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
      }#i 								
    }#t 
    
    
    ##----- DETECTION PROCESS -----## 
    
    for(t in 1:n.years){
      
      sigma[t] ~ dunif(0,4)
      betaResponse[t] ~ dunif(-5,5)
      betaResponseOth[t] ~ dunif(-5,5)
      
      for(c in 1:n.covs){
        betaCovs[c,t] ~ dunif(-5,5)
      }
      
      for(c in 1:n.covsOth){
        betaCovsOth[c,t] ~ dunif(-5,5)
      }
      
      for(c in 1:n.counties){
        p01[c,t] ~ dunif(0,1)
        p0[c,t] <- p01[c,t] * countyToggle[c,t]          ## toggle counties
      }#c
      
      for(c in 1:n.countries){
        p01Oth[c,t] ~ dunif(0,1)
        p0Oth[c,t] <- p01Oth[c,t] * countyToggleOth[c,t] ## toggle countries
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
    
    
    ##----- DERIVED PARAMETERS -----##
    
    for(t in 1:n.years){
      for(i in 1:n.individuals){ 
        isAlive[i,t] <- (z[i,t] == 2) 
      }#i
      N[t] <- sum(isAlive[1:n.individuals,t])
    }#t
    
  })
  
  
  
  ## ------ 2. NIMBLE CONSTANTS ------
  
  nimConstants <- list( n.individuals = dim(y.alive)[1],
                        n.detectors = dim(y.alive)[2],  
                        n.years = dim(y.alive)[3], 
                        n.years1 = dim(y.alive)[3]-1, 
                        n.covs = dim(detCovs)[3],
                        n.covsOth = dim(detCovsOth)[3],
                        numHabWindows = nHabCells,
                        n.countries = max(detCountries)+1,# + 1 for Norrbotten
                        n.counties = max(detCounties),
                        y.max = dim(habIDCells.mx)[1],
                        x.max = dim(habIDCells.mx)[2])
  
  
  
  ## ------ 3. NIMBLE INITS ------
  
  ## ------    3.1. GENERATE z KNOWN VALUES ------
  
  z <- apply(y.alive, c(1,3), function(x) any(x>0))
  z <- ifelse(z, 2, NA)
  z <- t(apply(z, 1, function(zz){
    if(any(!is.na(zz))){
      range.det <- range(which(!is.na(zz)))
      zz[range.det[1]:range.det[2]] <- 2
    }
    return(zz)
  }))
  
  
  
  ## ------    3.2. GENERATE z INITIAL values ------
  
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
  
  
  
  ## ------    3.3. GENERATE detResponse INITIAL VALUES ------
  
  detResponse <- already.detected 
  detResponse[rownames(detResponse) %in% "Augmented" ,1]  <- NA
  InitsDetResponse <- detResponse
  InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1,0.5)
  InitsDetResponse[!is.na(detResponse)] <- NA
  
  
  
  ## ------ 4. NIMBLE DATA ------
  
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
                   habitatGrid = habIDCells.mx)
  
  
  
  ## ------ 5. NIMBLE PARAMETERS ------
  
  nimParams <- c("N", "betaDens", "lambda", "dmean",
                 "omeg1", "gamma", "phi",
                 "pResponse", "sigma",
                 "p0", "betaResponse", "betaCovs",
                 "p0Oth", "betaResponseOth", "betaCovsOth")
  
  nimParams2 <- c("z", "sxy")
  
  
  
  ## ------ 6. CONVERT TO CACHED DETECTORS & SPARSE MATRIX ------
  
  ## ------    6.1. RESCALE COORDINATES  ------
  
  ## HABITAT
  ScaledLowerCoords <- scaleCoordsToHabitatGrid(
    coordsData = lowerHabCoords,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid = T)$coordsDataScaled
  
  ScaledUpperCoords <- scaleCoordsToHabitatGrid(
    coordsData = upperHabCoords,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid = T)$coordsDataScaled
  ScaledUpperCoords[ ,2] <- ScaledUpperCoords[ ,2]+1
  ScaledLowerCoords[ ,2] <- ScaledLowerCoords[ ,2]-1
  
  ## DETECTORS
  colnames(detector.xy) <- c("x","y")
  ScaledDetectors <- scaleCoordsToHabitatGrid(
    coordsData = detector.xy,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid =T )$coordsDataScaled
  
  ## ADD TO NIMDATA
  nimData$detector.xy <- as.matrix(ScaledDetectors)          
  nimData$lowerHabCoords <- as.matrix(ScaledLowerCoords)
  nimData$upperHabCoords <- as.matrix(ScaledUpperCoords)
  
  
  
  ## ------    6.2. CREATE CACHED DETECTORS OBJECTS ------
  
  ## [CM] reduce multiplicator to 3 ?????
  maxDistReCalc <- 2.1 * myVars$DETECTIONS$maxDetDist 
  
  DetectorIndexLESS <- GetDetectorIndexLESS(
    habitat.mx = myHabitat$habitat.mx,
    detectors.xy = nimData$detector.xy,
    maxDist = maxDistReCalc/res(myHabitat$habitat.r)[1],
    ResizeFactor = 1,
    plot.check = TRUE)
  
  ## ADD TO nimConstants
  nimConstants$y.maxDet <- dim(DetectorIndexLESS$habitatID)[1]
  nimConstants$x.maxDet <- dim(DetectorIndexLESS$habitatID)[2]
  nimConstants$ResizeFactor <- DetectorIndexLESS$ResizeFactor
  nimConstants$n.cellsSparse <- dim(DetectorIndexLESS$detectorIndex)[1]
  nimConstants$maxNBDets <- DetectorIndexLESS$maxNBDets
  
  ## ADD TO nimData
  nimData$detectorIndex <- DetectorIndexLESS$detectorIndex
  nimData$nDetectorsLESS <- DetectorIndexLESS$nDetectorsLESS
  nimData$habitatIDDet <- DetectorIndexLESS$habitatID
  
  
  
  ## ------    6.3. TRANSFORM Y TO SPARSE MATRICES ------
  
  ## STRUCTURED
  SparseY <- GetSparseY(y.aliveStructured)
  ## ADD TO nimData
  nimData$y.alive <- SparseY$y 
  nimData$yDets <- SparseY$yDets
  nimData$nbDetections <- SparseY$nbDetections
  ## ADD TO nimConstants
  nimConstants$nMaxDetectors <- SparseY$nMaxDetectors
  
  ## OTHER
  SparseYOth <- GetSparseY(y.aliveOthers)
  ## ADD TO nimData
  nimData$y.aliveOth <- SparseYOth$y 
  nimData$yDetsOth <- SparseYOth$yDets
  nimData$nbDetectionsOth <- SparseYOth$nbDetections
  ## ADD TO nimConstants
  nimConstants$nMaxDetectorsOth <- SparseYOth$nMaxDetectors
  
  
  
  ## ------ 7. NIMBLE INITS (CONTINUED) ------
  
  ## ------    7.1.GENERATE sxy INITIAL VALUES ------
  
  ## Project death to the next year
  myData.deadProj <- myData.dead[ ,c("Id","Year")]
  myData.deadProj$Year <- myData.deadProj$Year + 1
  ## Remove dead reco occuring the last year (not used)
  myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]
  
  ## Create a data.frame with all detection of all Individuals detected
  AllDets <- rbind(myData.alive$myData.sp [,c("Id","Year")],
                   myData.deadProj[ ,c("Id","Year")])
  AllDetections <- as.data.frame(AllDets)
  AllDetsxy <- st_coordinates(AllDets) 
  colnames(AllDetsxy) <- c("x","y")
  
  ## Rescale detection coordinates
  AllDetsxyscaled <- scaleCoordsToHabitatGrid(
    coordsData = AllDetsxy,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid =T )$coordsDataScaled
  AllDetections <- cbind(AllDetections, AllDetsxyscaled)
  
  ## GENERATE sxy INITIAL VALUES
  sxy.init <- getSInits(
    AllDetections = AllDetections,
    Id.vector = y.ar$Id.vector,
    idAugmented = which(rownames(z) %in%"Augmented"),
    lowerCoords = nimData$lowerHabCoords,
    upperCoords = nimData$upperHabCoords,
    habitatGrid = nimData$habitatGrid,
    intensity = NULL,
    sd = 4,
    movementMethod = "dbernppACmovement_normal")
  
  ## RESCALE sxy INITIAL VALUES 
  sxy.initscaled <- scaleCoordsToHabitatGrid(
    coordsData = sxy.init,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid = F)$coordsDataScaled
  
  
  
  ## ------ 9. LIST NIMBLE INITS & SAVE NIMBLE INPUTS ------
  
  for(c in 1:4){
    
    ## ------  9.1. LIST NIMBLE INITS ------
    nimInits <- list( "sxy" = sxy.init,
                      "dmean" = runif(1,0,10),
                      "z" = z.init,
                      "omeg1" = c(0.5,0.5),
                      "gamma" = runif(dim(y.alive)[3]-1,0,1),
                      "p01" = array(runif(18,0,0.2),
                                    c(nimConstants$n.counties, dim(y.alive)[3])),
                      "p01Oth" = array(runif(18,0,0.2),
                                       c(nimConstants$n.countries+1, dim(y.alive)[3])),
                      "sigma" = runif(nYears,1,4),
                      "betaDens" = runif(1,-0.1,0.1),
                      "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),
                                          c(dim(detCovsOth)[3], nYears)),
                      "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),
                                             c(dim(detCovsOth)[3], nYears)),
                      "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
                      "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
                      "detResponse" = InitsDetResponse,
                      "pResponse"  = runif(1, 0.4, 0.5),
                      "phi" = runif(dim(y.alive)[3]-1,0.1,0.3))
    
    
    ## An extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
    nimInits$sxy <- round(nimInits$sxy, 5)
    
    nimConstants$countyToggle <- nimInits$p01
    nimConstants$countyToggle[] <- 1
    
    yearsNotSampled <- which(!years%in% yearsSampledNorrb)
    for(t in yearsNotSampled){ nimConstants$countyToggle[1,t] <- 0 }
    
    ## add another category to detcountry if in norrbotten, to turnoff detection to 0 there. 
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
    nimData$detCountries <- detCountriesNorb
    
    nimConstants$countyToggleOth <- nimInits$p01Oth
    nimConstants$countyToggleOth[] <- 1
    yearsNotSampled <- which(!years%in% yearsSampledNorrb)
    for(t in yearsNotSampled){
      nimConstants$countyToggleOth[3,t] <- 0
    }
    
    
    
    ## ------  9.2. SAVE NIMBLE INPUTS ------
    
    save(nimData,
         nimConstants,
         y.dead,
         nimParams,
         nimParams2,
         modelCode,
         nimInits,
         file = file.path(myVars$WD, myVars$modelName, thisSex,
                          paste0(myVars$modelName, thisSex,"_Chain", c, ".RData")))
  }#c
}#thisSex



##------------------------------------------------------------------------------