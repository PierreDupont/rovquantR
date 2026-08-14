##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE analysis 2026 - WIP script
##
## This R script presents the process of packaging the Wolverine analysis of 2025,
## known as "the Chaotic Estimation!" (cf. CM).  OPSCR and SCR models were ran using the model
## 54.Cleaned. ("54.Cleaned2025TestPDScript.R").Results from the OPSCR were extracted 
## using "PlotWolverine54Cleaned2025.R" and "PlotWolverine54Cleaned2025Snap_ASPcm.R" 
## for the SCR model. Figures and tables that combined results from the OPSCR 2024 
## and the SCR model were created using "CombineAndPlotFigureTable.R".
##
## This file combines all four scripts and consists of > 18000 lines before cleaning.
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 10/08/2026
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2026
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
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(xtable)
library(terra)
library(colorspace)

library(rovquantR)

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

data.dir <- file.path(dir.dropbox, "wolverine/2025/Data")
working.dir <- file.path(dir.dropbox, "wolverine/2025/Test_OLD")


## ------ SOURCE THE REQUIRED FUNCTIONS ------

source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")
sourceDirectory(dir.function, modifiedOnly = FALSE)
#sourceDirectory(dir.function.nimble, modifiedOnly = FALSE)
load(file.path(dir.dropbox,"DATA/MISC DATA/age.lookup.table.RData"))


##------------------------------------------------------------------------------

## ------ 0.SET ANALYSIS CHARACTERISTICS -----

## WORKING DIRECTORY & MODEL NAME
modelName = "54.Cleaned2025"

## HABITAT SPECIFICATIONS
HABITAT = list( countries =  c("SWE","NOR"),
                habResolution = 20000,
                habBuffer = 60000)

## NGS DATA SPECIFICATIONS
DATA = list( years = 2015:2024, 
             species = c("Jerv"),              
             sex = c("Hann","Hunn"),                   
             sampling.months = list(12,1:6)) 

## DETECTORS SPECIFICATIONS
DETECTORS = list( detSubResolution = 2000,
                  detResolution = 10000,
                  detDeadResolution = 15000)

## DATA GENERATION
DETECTIONS = list( maxDetDist = 40000,
                   resizeFactor = 1,
                   aug.factor = 0.8)

## OUTPUT PLOTS
OUTPUT = list(mapResolution = 10000)

## MISCELLANEOUS
plot.check = TRUE

years <- DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))



##------------------------------------------------------------------------------

## ------ I. LOAD AND SELECT DATA ------

## ------   1. HABITAT DATA ------

## ------     1.1. LOAD RAW SHAPEFILES ------

## POLYGONS OF THE REGION
GLOBALMAP <- st_read(file.path( dir.dropbox,
                                "DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) 
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP[GLOBALMAP$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- COUNTRIES %>%
  group_by(ISO) %>%
  summarize()

## POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp"))   
COMMUNES_SWE <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp"))   
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)

## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- COMMUNES %>%
  group_by(NAME_1) %>%
  summarize()

## Identify Norrbotten county
COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018,2023,2024)

## AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATE <- COUNTIES
COUNTIES_AGGREGATE$id <- 1:nrow(COUNTIES_AGGREGATE)
COUNTIES_AGGREGATE$id[c(24,3,15,9,14,38,40,21,27,37,31,26,34,5,8,12,36,13,7)] <- 3
COUNTIES_AGGREGATE$id[c(39,33,23,32,29,22,4,11,20,2,10,16,25,1)] <- 4
COUNTIES_AGGREGATE$id[c(19)] <- 1
COUNTIES_AGGREGATE$id[c(35)] <- 2
COUNTIES_AGGREGATE$id[c(17,28)] <- 5
COUNTIES_AGGREGATE$id[c(18)] <- 7
COUNTIES_AGGREGATE$id[c(30)] <- 8
COUNTIES_AGGREGATE <- COUNTIES_AGGREGATE %>% group_by(id) %>% summarize()
COUNTIES_AGGREGATED <- st_simplify(COUNTIES_AGGREGATE,preserveTopology = T,dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id
ggplot(COUNTIES_AGGREGATED) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))



## ------     1.2. CREATE STUDY AREA POLYGON ------

## CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
if(!is.null(HABITAT$countries)){
  myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% HABITAT$countries, ]
  
  ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
  myBufferedArea <- st_buffer(st_as_sf(myStudyArea) ,dist = HABITAT$habBuffer)
  myBufferedArea$id <- 1
  myBufferedArea <- myBufferedArea %>% group_by(id) %>% summarize()
  myBufferedArea <- st_intersection(myBufferedArea, GLOBALMAP)
}

myStudyArea$id <- 1
myStudyArea <- myStudyArea %>% group_by(id) %>% summarize()

## PLOT CHECK
if(plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myBufferedArea), add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
  plot(st_geometry(myStudyArea), add = TRUE, col ="red")
}



## ------   2. NGS DATA ------

## ------     2.1. LOAD ROVBASE FILES ------

## NGS data from RovBase
DNA <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dna_wolverines.csv"),
                 fileEncoding = "latin1") 

## Dead Recoveries from RovBase
DEAD <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20251121/dead_carnivores.csv"),
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
#rovbaseObs4 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/all carnivores NGS/RIB15092025172035833.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3)#,rovbaseObs4)
rm(list = c("rovbaseObs1", "rovbaseObs2", "rovbaseObs3"))#, "rovbaseObs4"))



## ------     2.2. TRANSLATE SCANDINAVIAN CHARACTERS ------

colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA))
## Drop a column that makes cleanDataNew to fail
DNA <- DNA[,-which(colnames(DNA)%in% "Kjoenn..Individ.")]
colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD))
colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN))
colnames(skandObs) <- translateForeignCharacters(dat=colnames(skandObs))
colnames(rovbaseObs) <- translateForeignCharacters(dat=colnames(rovbaseObs))
rovbaseObs$Proevetype <- translateForeignCharacters(dat=rovbaseObs$Proevetype)



## ------   3. SEARCH EFFORT DATA ------

## ------     3.1. GPS SEARCH TRACKS ------

## LOAD GPS SEARCH TRACKS
TRACKS_SINGLE <- read_sf(paste(dir.dropbox,
                               "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp", sep = ""))
TRACKS_MULTI <- read_sf(paste(dir.dropbox,
                              "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp", sep = ""))

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

## Check that we only have wolverine tracks
table(ALL_TRACKS$Jerv) == nrow(ALL_TRACKS)
## Check
hist(ALL_TRACKS$Yr)

## GET EXTENT
myStudyArea.extent <- st_bbox(extent(myStudyArea))
st_crs(myStudyArea.extent) <- st_crs(COUNTRIES)

## [PD] ASP ADDED REMOVAL OF DUPLICATED TRACKS
dupIDs <- dupDist <- TRACKS_YEAR <- list()
for(t in 1:nYears){
  
  TRACKS <- ALL_TRACKS %>% 
    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
    filter( Yr%in%YEARS[[t]][1] & Mth %in% DATA$samplingMonths[[1]] |
              Yr%in%YEARS[[t]][2] & Mth%in%DATA$samplingMonths[[2]]) %>%
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

## PLOT CHECK
if(plot.check){
  
  par(mfrow = c(2,2))
  
  ## Length of tracks searched per year
  lengthPerYear <- unlist(lapply(TRACKS_YEAR,function(x) sum(x$dist)/1000))
  names(lengthPerYear) <- years
  barplot(lengthPerYear, ylab = "Track length (km)", main = "Length of tracks searched per year")
  
  ## Number of tracks searched per year
  numPerYear <- unlist(lapply(TRACKS_YEAR,function(x) length(unique(x$ID))))
  names(numPerYear) <- years
  barplot(numPerYear, ylab = "Number of tracks", main = "Number of tracks searched per year")
  
  ## Length of tracks duplicated per year
  dupdist <- unlist(lapply(dupDist,function(x) sum(x)/1000))
  names(dupdist) <- years
  barplot(dupdist,ylab = "Track length (km)", main = "Length of tracks duplicated per year")
  
  ## Number of tracks duplicated per year
  dup <- unlist(lapply(dupIDs,length))
  names(dup) <- years
  barplot(dup, ylab = "Number of tracks", main = "Number of tracks duplicated per year")
}



## ------     3.2. DISTANCE TO ROADS ------

## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(file.path( dir.dropbox,
                                  "DATA/GISData/Roads/MinDistAllRoads1km.tif"))

## RASTERIZE DISTANCE TO ROADS
r <- fasterize(st_as_sf(myStudyArea), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea)

## PLOT CHECK
if(plot.check){
  par(mfrow = c(1,1))
  plot(DistAllRoads)
  plot(st_geometry(myStudyArea), add = T)
}



## ------     3.3. DAYS OF SNOW ------

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
SNOW <- stack(file.path( dir.dropbox, 
                         "DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2014_2025_Wolverine.tif")) # Average snow from December to June (the official monitoring period for Norway&Sweden)
## RENAME THE LAYERS
#names(SNOW) <- paste(2008:2023,(2008:2023)+1, sep="_")
names(SNOW) <- paste(2014:2024,(2014:2024)+1, sep="_")
## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep = "")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))



##------------------------------------------------------------------------------

## ------ II. CREATE OPSCR DATA ------

## ------   1. CLEAN & FILTER NGS DATA ------

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]
dim(DNA)

## Remove un-verified dead recoveries 
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]
dim(DEAD)



## ------     1.1. CLEAN NGS & DEAD RECOVERY DATA ------

myCleanedData.sp <- CleanDataNew2sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(DATA$sampling.months)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)


##-- Plot check
if(plot.check){
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myStudyArea), add = T, col = "red")
  plot(st_geometry(myCleanedData.sp), add = TRUE, pch = 19, cex = 0.2, col = "blue")
}

myFullData.sp <- FilterDatasf(
  myData = myCleanedData.sp,
  poly = myStudyArea,
  dead.recovery = T ,
  sex = c("Hann","Hunn"), # Do the sex selection at the last moment
  setSex = T)

if(plot.check){
  plot( st_geometry(myFullData.sp$alive), add = TRUE, pch = 19, cex = 0.2, col = "white")
}

## REMOVE SUSPECT SAMPLES ACCORDING TO HENRIK (THERE ARE NOT SUSPECT SAMPLES IN SWEDEN)
myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
dim(myFullData.sp$alive)

##-- EXPORT THE DATA
myFullData.spSaved <- myFullData.sp

##-- REMOVE SUSPECT DEAD RECOVERIES ACCORDING TO HENRIK
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
dim(myFullData.sp$dead.recovery)

##-- REMOVE CUBS OF THE YEAR ACCORDING TO EVA HEDMARK (email 23/11/2025, document with explanation in AQEG Dropbox\AQEG Team Folder\RovQuant\DATA\RovbaseData\ROVBASE DOWNLOAD 20251121)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% c("D608410", "D608411", "D605997")),]

##-- Remove individuals that died twice
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

unique(myFullData.sp$dead.recovery$DeathCause)
unique(myFullData.sp$dead.recovery$DeathCause_2)

## Dead recoveries flagged by Henrik that should always be removed (email from the 18/12/2024)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$RovBaseId %in% c("M495994", "M524051", "M524052", "M524053"), ]

## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) Remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[
  -which(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
           myFullData.sp$dead.recovery$Month > 2 &
           myFullData.sp$dead.recovery$Month < 12), ]

## 2) Remove individuals that have a weight >0 and <4 between March and November format the weight correctly
myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)

##-- Convert to decimals
myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))

##-- Get the two weight columns together.
myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt),
                                             myFullData.sp$dead.recovery$Helvekt,
                                             myFullData.sp$dead.recovery$Slaktevekt)

##-- Assign negative values to NAs to avoid issues
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999
dim(myFullData.sp$dead.recovery)

##-- check with Henrik
# this step does not remove dead recoveries on id with weight==0 should it?
##-- WEIGTH DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$weight[(myFullData.sp$dead.recovery$weight >-1) &
                                            myFullData.sp$dead.recovery$Month%in% t],breaks=c(0:30), main=t,xlab="Weight")
}

##-- AGE DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$Age[(myFullData.sp$dead.recovery$Age >-1) &
                                         myFullData.sp$dead.recovery$Month%in% t],breaks=seq(-0.01,20.99,by=1),
       main=t,
       xlab="Age")
}

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

##-- Checks
dim(myFullData.sp$alive)
table(myFullData.sp$alive$Year)
dim(myFullData.sp$dead.recovery)
table(myFullData.sp$dead.recovery$Year)
table(myFullData.sp$dead.recovery$DeathCause)
table(myFullData.sp$dead.recovery$DeathCause_2)



## ------     1.2. FILTER NGS & DEAD RECOVERY DATA ------

myFilteredData.sp <- myFullData.sp
dim(myFullData.sp$alive)

##-- Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$alive)

myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]
table(myFilteredData.sp$dead.recovery$Year)
dim(myFilteredData.sp$dead.recovery)

##-- Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in% unlist(DATA$sampling.months), ]
table(myFilteredData.sp$alive$Month)

##-- Checks
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$dead.recovery)
table(myFilteredData.sp$dead.recovery$Year)
table(myFilteredData.sp$dead.recovery$DeathCause)
table(myFilteredData.sp$dead.recovery$DeathCause_2)



## ------     1.3. REMOVE DETECTIONS IN NORRBOTTEN IN ALL YEARS EXCEPT 2017, 2018 and 2019 ------ 

is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))
sum(is.Norr, na.rm = T)

## Check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(is.Norr)), ]$Year) %>% sum()

## Remove samples in Norrbotten in years not sampled
myFilteredData.sp$alive <- myFilteredData.sp$alive[- which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                                             !is.na(is.Norr)), ]
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)

## Plot check
for(t in 1:nYears){
  plot( st_geometry(myStudyArea))
  plot( st_geometry(COUNTIESNorrbotten), add = T, col = "blue")
  plot( st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t], ]), col = "red", add = T, pch = 16)
}#t



## ------     1.4. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

## ------       1.4.1. ASSIGN SAMPLES TO TRACKS ------

## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
myFilteredData.sp$alive$TrackRovbsID <- NA
myFilteredData.sp$alive$TrackDist <- NA
TRACKSSimple_sf <- list()
for(t in 1:nYears){
  TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbsID)
}

## CREATE A BUFFER AROUND EACH DETECTION
dnatemp <- st_as_sf(myFilteredData.sp$alive)
tmp <- st_buffer(dnatemp, dist = 750)

## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
for(i in 1:nrow(myFilteredData.sp$alive)){
   ## MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
   t <- which(years %in% tmp[i, ]$Year)
   whichSameDate <- which(as.character(TRACKS_YEAR[[t]]$Dato) == as.character(myFilteredData.sp$alive$Date[i]))

   ## INTERSECT POINT WITH TRACKS
   tmpTRACKS <- st_intersection(TRACKS_YEAR[[t]][whichSameDate, ], tmp[i, ])

   ## If not TRACKS that date, move on
   if(nrow(tmpTRACKS)==0){next}

   ## Else, find the closest TRACK
   dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)

   ## IF NO MATCHING DATE TRACK ASSIGN NA
   if(length(dist)==0){
     myFilteredData.sp$alive$TrackRovbsID[i] <- NA
     myFilteredData.sp$alive$TrackDist[i] <- NA
   }
   ## IF 1 MATCHING TRACK ASSIGN TO THAT TRACK
   if(length(dist)==1){
     myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
     myFilteredData.sp$alive$TrackDist[i] <- dist
   }
   ## IF SEVERAL MATCHING DATES ASSIGN TO THE CLOSEST OF THE MATCHING TRACKS
   if(length(dist)>1){
     myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
     myFilteredData.sp$alive$TrackDist[i] <- min(dist)
   }
}#i

## SAVE FOR FASTER LOADING
save( myFilteredData.sp,
      file = file.path(working.dir, "data/myFilteredData_original.RData"))
load(file.path(working.dir, "data/myFilteredData_original.RData"))
# sum(myFilteredData.sp$alive$TrackDist[myFilteredData.sp$alive$Sex %in% "Hunn"], na.rm = T)
# nrow(myFilteredData.sp$alive[myFilteredData.sp$alive$Sex %in% "Hunn", ] )



## ------       1.4.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC & STRUCTURED ------

distanceThreshold <- 500

## Proeveleverandoer columns was replaced by two columns, merging them now...
myFilteredData.sp$alive$Proeveleverandoer <- ifelse(myFilteredData.sp$alive$Annen.innsamler...Rolle %in% "", 
                                                    myFilteredData.sp$alive$Samlet.selv...Rolle,
                                                    myFilteredData.sp$alive$Annen.innsamler...Rolle)
# table(myFilteredData.sp$alive$Proeveleverandoer, useNA = "always")
# sum(duplicated(myFilteredData.sp$alive$RovbaseID))
# sum(duplicated(myFilteredData.sp$alive$DNAID))

## Structured samples are:
##    1. samples from agencies ("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen")  
##    2. samples assigned to a GPS search track
##    3. wihtin distanceThreshold of a GPS search track
whichStructured <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold
# table(whichStructured, useNA = "always")
myFilteredData.spStructured <- myFilteredData.sp$alive[whichStructured, ]
myFilteredData.spOthers <- myFilteredData.sp$alive[!whichStructured, ]

## CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
nrow(myFilteredData.spStructured) + nrow(myFilteredData.spOthers)
nrow(myFilteredData.sp$alive)

## CHECK THAT ALL RESULTS FROM HAIR TRAPS ARE SET TO OPPORTUNISTIC
whichHair <- which(myFilteredData.sp$alive$DNAID%in% HairTrapSamples$DNAID)

## Plot check 
plot(COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]$geometry)
plot( myFilteredData.sp$alive[whichHair,]$geometry, add = T, col = "red", pch = 16)
if(length(which(myFilteredData.spStructured$alive$DNAID%in% HairTrapSamples$DNAID))>0){
  print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
}



## ------     1.5. SEPARATE MORTALITY CAUSES ------ 

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

# ##-- Plot check
# if(plot.check){
#   par(mfrow = c(1,3))
#   for(t in 1:nYears){
#     ## DEAD RECOVERIES TOTAL
#     tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
#     NGS_TabTotal <- table(tempTotal$Country)
#     ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
#     ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
#     tempIn <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]
#     NGS_TabIn <- table(tempIn$Country)
#     ID_TabIn <- apply(table(tempIn$Id, tempIn$Country), 2, function(x) sum(x>0))
#     ## PLOT NGS SAMPLES
#     plot(st_geometry(GLOBALMAP), col="gray80")
#     plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
#     plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
#     plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
#     ## ADD NUMBER OF NGS samples and IDs per COUNTRY
#     graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
#     graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
#     ## ADD OVERALL NUMBERS
#     mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
#     mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
#   }#t
#   
#   ## PLOT TREND DETECTIONS AND DEAD RECOVERIES OVER TIME AND SPACE
#   ## DETECTIONS
#   #pdf(file=file.path(working.dir, modelName, paste(modelName,"_TRENDDetections.pdf",sep="")))
#   temp <- unique(myFilteredData.sp$alive[ ,c("Year","Country","DNAID")])
#   tab_Country.Year <- table(temp$Year, temp$Country)
#   country.colors <- c("goldenrod1","goldenrod3")
#   par(mfrow=c(1,1), mar=c(5,5,5,5))
#   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
#   lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
#   lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
#   legend("bottomright",c("N","S"), fill=country.colors)
#   
#   ## ID DETECTED
#   temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country,myFilteredData.sp$alive$Id)
#   tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
#   country.colors <- c("goldenrod1","goldenrod3")
#   par(mfrow=c(1,1), mar=c(5,5,5,5))
#   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
#   lines(tab_Country.Year1[,"N"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
#   lines(tab_Country.Year1[,"S"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
#   legend("bottomright",c("N","S"), fill=country.colors)
#   
#   ## Average number of detection per detected ID  
#   tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
#   par(mfrow=c(1,1), mar=c(5,5,5,5))
#   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
#        ylab="Average Number of detections", xlab="Years")
#   lines(tab_Country.Year2[,"N"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
#   lines(tab_Country.Year2[,"S"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
#   legend("bottomright",c("N","S"), fill=country.colors)
#   
#   ## deadrecovery 
#   temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")])
#   tab_Country.Year <- table(temp$Year, temp$Country)
#   country.colors <- c("goldenrod1","goldenrod3")
#   par(mfrow=c(1,1), mar=c(5,5,5,5))
#   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
#        ylab="N Id Dead recovered", xlab="Years")
#   lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
#   lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
#   legend("topright",c("N","S"), fill=country.colors)
#   
#   dev.off()
# }



## ------   2. GENERATE HABITAT ------

## ------     2.1. REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ------

# ## DELINEATE A BUFFER AROUND ALL DEADRECO 
# BuffDead <- st_buffer(myFilteredData.sp$dead.recovery, dist = HABITAT$habBuffer)
# BuffDead$idd <- 1
# BuffDead <- BuffDead %>% group_by(idd) %>% summarize()

##-- DELINEATE A BUFFER AROUND ALL DETECTIONS 
myBufferedArea <- st_buffer( myFilteredData.sp$alive,
                             dist = HABITAT$habBuffer * 1.4) %>%
  mutate(id = 1) %>%
  group_by(id) %>% 
  summarize()

##-- CUT TO SWEDISH & NORWEGIAN BORDERS
myStudyArea <- st_intersection(myBufferedArea, myStudyArea)

myStudyArea$idd <- 1
myStudyAreaAggregated <- myStudyArea %>% group_by(idd) %>% summarize()

##-- Plot check 
# if(plot.check){
#   par(mar=c(0,0,0,0))
#   plot(st_geometry(myStudyArea))
#   plot(st_geometry(myFilteredData.sp$alive), pch=21, bg="red", cex=0.5,add=T)
#   plot(st_geometry(myFilteredData.sp$dead.recovery), pch=21, bg="blue", cex=0.5,add=T)
#   plot(st_geometry(myBufferedArea), border="red", add=T)
#   plot(st_geometry(myStudyArea), border="grey", add=T)
# }



## ------     2.2. GENERATE HABITAT CHARACTERISTICS FROM THE NEW HABITAT DEFINITION ------

myHabitat <- MakeHabitatFromRastersf( 
  poly = myStudyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = HABITAT$habBuffer,                               
  plot.check = T)

##-- Plot check 
# if(plot.check){
# par(mfrow = c(1,2))
# plot(myHabitat$habitat.r,legend = F)
# plot(st_geometry(myStudyArea), add = T)
# }

## RETRIEVE HABITAT WINDOWS BOUNDARIES
lowerHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1, ] - 0.5 * HABITAT$habResolution
upperHabCoords <- coordinates(myHabitat$habitat.r)[myHabitat$habitat.r[]==1, ] + 0.5 * HABITAT$habResolution
nHabCells <- dim(lowerHabCoords)[1]

## CREATE HABITAT GRID 
habIDCells.mx <- myHabitat$IDCells.mx 
habIDCells.mx[] <- 0
scaledHabGridCenters <- scaleCoordsToHabitatGrid(
  coordsData = myHabitat$habitat.xy,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid = F)$coordsHabitatGridCenterScaled

scaledHabGridCenters <- scaledHabGridCenters[myHabitat$habitat.r[] == 1, ]
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}#i



## ------     2.3. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 

##-- Remove samples outside the STUDY AREA #[CM]
whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyArea))))
if(length(whichOut)>0){ myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ] }
myFilteredData.sp$alive$Id <- droplevels(myFilteredData.sp$alive$Id)

##-- Remove dead recoveries outside the HABITAT #[CM] 
whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, myHabitat$buffered.habitat.poly))))
if(length(whichOutBuff)>0){ myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ] }

##-- Plot check 
# if(plot.check){
#   # par(mfrow = c(1,2))
#   plot(myHabitat$habitat.r)
#   plot(st_geometry(myStudyArea), add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
#   plot(st_geometry(GLOBALMAP), add = T)
#   plot(st_geometry(myHabitat$buffered.habitat.poly), add=T)
#   plot(st_geometry(myFilteredData.sp$alive),pch=21, bg="red", cex=0.5,add=T)
#   plot(st_geometry(myFilteredData.sp$dead.recovery),pch=21, bg="blue", cex=0.5,add=T)
# }
# 
# ## check correlation number of detections ~ between monitoring season
# myFilteredData.sp$dead.recovery$Id <- as.character(myFilteredData.sp$dead.recovery$Id)
# myFilteredData.sp$alive$Id <- as.character( myFilteredData.sp$alive$Id)
# deadID <- unique(myFilteredData.sp$dead.recovery$Id)
# ndet <- NULL
# timeDiff <- NULL
# for(i in 1:length(deadID)){
#   tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Year
#   timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i], ]$Date -
#     as.POSIXct(strptime(paste("01-12", tmpYear, sep="-"), "%d-%m-%Y")) 
#   ndet[i] <- length(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] & myFilteredData.sp$alive$Year %in% tmpYear,])
# }
# 
# ##-- Plot check 
# # pdf(file = file.path(working.dir, modelName,"Prop id detected_Time available.pdf"))
# plot(ndet ~ timeDiff,
#      ylab = "Total number of detections",
#      xlab = "Number of days between dec 1 and dead recovery")
# hh <- hist(timeDiff[ndet > 0], breaks = seq(0, 400, by = 25))
# hh1 <- hist(timeDiff[ndet == 0], breaks = seq(0, 400, by = 25))
# barplot(rbind(hh$counts/(hh$counts+hh1$counts),
#               hh1$counts/(hh$counts+hh1$counts)),
#         names.arg = hh$breaks[1:(length(hh$breaks)-1)],
#         xlab = "number of days between dead reco and start monitoring",
#         ylab = "%")
# legend("topright",
#        fill = c(grey(0.2),grey(0.8)),
#        legend = c("detected","notDetected"))
# # dev.off()



## ------     2.4. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       2.4.1. DEN COUNTS ------

DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X","UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[ ,"id"]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD( as(DEN.sp,"Spatial"),
              h = 30000,
              grid = as(myHabitat$habitat.r, 'SpatialPixels'))))

if(plot.check){
  plot(DEN.r)
  plot(st_geometry(myStudyArea), add = TRUE, border = "black")
}

##-- EXTRACT COVARIATE
denCounts <- DEN.r[myHabitat$habitat.r[ ] == 1]
denCounts <- round(scale(denCounts), digits = 2)



## ------   3. GENERATE DETECTORS ------

## ------     3.1. GENERATE DETECTORS CHARACTERISTICS ------

##-- GENERATE SUB-DETECTORS BASED ON THE STUDY AREA
habitat.subdetectors <- disaggregate(
  myHabitat$habitat.rWthBuffer,
  fact = res(myHabitat$habitat.r)[1]/DETECTORS$detSubResolution)

##-- GENERATE NGS DETECTORS BASED ON THE STUDY AREA
myDetectors <- myDetectors.dead <- MakeSearchGridsf(
  data = habitat.subdetectors,
  resolution = DETECTORS$detResolution,
  div = (DETECTORS$detResolution/DETECTORS$detSubResolution)^2,
  plot = FALSE,
  fasterize = TRUE)

##-- EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(myDetectors$main.detector.sp)[1]
n.detectors.dead <- dim(myDetectors.dead$main.detector.sp)[1]

##-- FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(myDetectors$main.detector.sp)
n.trials <- as.vector(table(myDetectors$detector.sp$main.cell.id))
detector.dead.xy <- st_coordinates(myDetectors.dead$main.detector.sp)

##-- IDENTIFY DETECTORS IN NORBOTTEN 
COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten","Troms","Västerbotten","Nordland","Finnmark"), ]
COUNTIESAroundNorrbotten <- st_simplify(COUNTIESAroundNorrbotten, dTolerance = 500)

##-- CREATE A NORRBOTTEN DETECTOR GRID
distDestsCounties <- st_distance(myDetectors$main.detector.sp, COUNTIESAroundNorrbotten, byid = T)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min) == 3)

##-- Plot check 
if(plot.check){
  plot(st_geometry(COUNTIESAroundNorrbotten))
  plot(st_geometry(myDetectors$main.detector.sp), col = "black", pch = 16, cex = 0.3, add = T)
  plot(st_geometry(myDetectors$main.detector.sp[detsNorrbotten, ]), col = "red", pch = 16, cex = 0.3, add = T)
}

##-- RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * DETECTORS$detResolution

# ##-- Plot check 
# if(plot.check){
#   par(mfrow = c(1,2))
#   ## PLOT NGS DETECTORS
#   plot( st_geometry(myHabitat$buffered.habitat.poly),
#         main = paste(n.detectors, "Detectors Alive"),
#         col = rgb(0.16,0.67,0.16, alpha = 0.3))  
#   plot( st_geometry(myStudyArea),
#         col = rgb(0.16,0.67,0.16, alpha = 0.5), add = TRUE)
#   plot( st_geometry(myDetectors$main.detector.sp),
#         col = "red", pch = 16, cex = 0.1, add = TRUE)
#   plot( st_geometry(COUNTRIES), add = TRUE)
#   ## PLOT DEAD DETECTORS
#   plot( st_geometry(myHabitat$buffered.habitat.poly),
#         main = paste(n.detectors.dead, "Detectors Dead"),
#         col = rgb(0.16,0.67,0.16, alpha = 0.3)) 
#   plot( st_geometry(myStudyArea),
#         add = T, col = rgb(0.16,0.67,0.16, alpha = 0.5))
#   plot( st_geometry(myDetectors.dead$main.detector.sp),
#         col = "red", pch = 16, cex = 0.1, add = TRUE)
#   plot( st_geometry(COUNTRIES), add = TRUE)
# }



## ------     3.2. GENERATE DETECTOR-LEVEL COVARIATES ------

## ------       3.2.1. EXTRACT COUNTRIES ------

dist <- st_distance(myDetectors$main.detector.sp, COUNTRIES, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))

##-- Plot check 
if(plot.check){
  par(mfrow = c(1,2))
  myCol <- c("blue4", "yellow1")
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Countries")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(myDetectors$main.detector.sp), col = myCol[detCountries], pch = 16, cex = 0.8, add=T)
  plot(st_geometry(COUNTRIES), add = TRUE)
}



## ------       3.2.2. EXTRACT COUNTIES ------

dist <- st_distance(myDetectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F )
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties), ]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))

##-- Plot check 
if(plot.check){
  myCol <- terrain.colors(nrow(COUNTIES_AGGREGATED))
  plot( st_geometry(GLOBALMAP),
        col = "gray80", main = "Aggregated Counties")
  plot( st_geometry(myStudyArea),
        col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot( st_geometry(myBufferedArea),
        col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot( st_geometry(myDetectors$main.detector.sp[detCounties%in% c(5), ]),
        col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
  plot( st_geometry(myDetectors$main.detector.sp),
        col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
  plot( st_geometry(COUNTIES_AGGREGATED), add = TRUE)
  text(st_geometry(COUNTIES_AGGREGATED), labels = COUNTIES_AGGREGATED$id, col = "black")  
  plot( st_geometry(myDetectors$main.detector.sp[detCounties %in% 3, ]),
        col = "red", pch = 16, cex = 0.8, add = T)
}



## ------       3.2.3. EXTRACT GPS TRACKS LENGTHS ------

## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp),
                                      rep(1,nrow(myDetectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)
#plot(st_geometry(detectorGrid))

## INITIALIZE MATRIX & RASTERS OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
TRACKS.r <- list()

## CALCULATE THE LENGTH OF THE TRACKS
for(t in 1:nYears){
  intersection <- st_intersection(detectorGrid, TRACKS_YEAR[[t]]) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN)) ## Get total length searched in each detector grid cell
  # transect_N = length(unique(ID))) ## Get total number of visits in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[ ,t]
  print(t)
}#t

##-- Plot check 
if(plot.check){
  max <- max(unlist(lapply(TRACKS.r, function(x) max(x[], na.rm = T))))
  cuts <- seq(0,max,length.out = 100) # set breaks
  col <- rev(terrain.colors(100))
  CountriesDetRes <- disaggregate(habitatRasters$Countries, fact = 2)
  CountriesDetRes <- crop(CountriesDetRes,TRACKS.r[[1]])
  rr <- TRACKS.r[[1]]
  rr[CountriesDetRes[] %in% 2] <- 1
  plot(rr)
  sum(st_length(TRACKS_YEAR[[t]]))/1000
  sum(TRACKS.r[[t]][], na.rm = T)/1000
  
  # pdf(file = file.path(working.dir, modelName, "Tracks.pdf"))
  NORTRACKS <- SWETRACKS <- 0
  for(t in 1:nYears){
    plot( TRACKS.r[[t]], main = years[t], breaks = cuts, col = col, legend = FALSE)
    plot( st_geometry(myHabitat$habitat.poly), main = years[t], add = T)
    plot( TRACKS.r[[t]], legend.only = TRUE,
          breaks = cuts, col = col, legend.width = 2,
          axis.args = list(at = round(seq(0, max, length.out = 5), digits = 1),
                           labels = round(seq(0, max, length.out = 5), digits = 1),
                           cex.axis = 0.6),
          legend.args=list(text='', side=4, font=2, line=2.5, cex=0.8))
    # points( myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ],
    #         col="red", pch=16, cex=0.8)
    ## summary tracks
    NORTRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 2],na.rm = T )/1000
    SWETRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 4],na.rm = T )/1000
  }#t
  
  years1 <- years + 1
  plot(SWETRACKS ~ (years1),  col = country.colors[2],
       lwd = 2, pch = 16, type = "b", ylim = c(0,300000), ylab = "sum tracks km")
  lines(NORTRACKS~(years1), col=country.colors[1], lwd=2, pch=16, type="b")
  legend("topright",c("N","S"), fill=country.colors)
  # dev.off()
}



## ------       3.2.4. EXTRACT DISTANCES TO ROADS ------

## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)

## if NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( DistAllRoads,
                        myDetectors$main.detector.sp[isna, ],
                        buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp

##-- Plot check 
# if(plot.check){
#   par(mfrow = c(1,1))
#   plot(st_geometry(GLOBALMAP), col = "gray80", main = "Distance to roads")
#   plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
#   plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
#   plot(DistAllRoads,add=T)
#   plot(st_geometry(myDetectors$main.detector.sp), cex=DoScale(detRoads), pch = 16, add = T)
# }



## ------       3.2.5. EXTRACT DAYS OF SNOW ------

##-- EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(myDetectors$main.detector.sp)[1], ncol = nYears)
det.sptransf <- st_transform(myDetectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

##-- if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp

##-- Plot check 
# if(plot.check){
#   plot( st_geometry(myDetectors$main.detector.sp),
#         cex = DoScale(detSnow[,6],l = 0,u = 0.5),
#         pch = 16)
# }



## ------       3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------         3.2.6.1. SKANDOBS ------

## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))

## MAKE IT SPATIAL 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(myStudyArea))

## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist(DATA$sampling.months))
skandObs$monitoring.season <- ifelse(skandObs$month > 12, skandObs$year, skandObs$year-1) 
skandObs <- skandObs[subset, ] 

## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(myHabitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1, ]
subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace,] 

## RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate(habitat.subdetectors, fact=(DETECTORS$detResolution/DETECTORS$detSubResolution))
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(skandObs[skandObs$monitoring.season %in% y, 1], r.detector, fun="count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})
r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))
plot(r.skandObsSamplesBinary[[t]])

##-- Plot check
if(plot.check){
  # pdf( file = file.path(working.dir, modelName,"skandObs.pdf"), width = 10)
  barplot(table(skandObs$monitoring.season ))
  barplot(table(skandObs$month ), xlab="Months")
  barplot(table(skandObs$species))
  ## MAPS 
  par(mar=c(0,0,2,0))
  for(t in 1:nYears){
    plot(st_geometry(myStudyArea), main = years[t])
    plot( st_geometry(skandObs[skandObs$monitoring.season %in% years[t], ]),
          pch = 16, col = "red", cex = 0.1, add = T)
  }
  # dev.off()
}


## ------         3.2.6.2. ROVBASE ------

## GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`), ]
# rovbaseObs$Funnetdato <- as.POSIXct(strptime(rovbaseObs$Funnetdato, "%d.%m.%Y")) 
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

## DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)

## SUBSET THE DATA 
filter <- list( 
  species = "Jerv",
  type = c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
            "Saliv/Spytt", "Loepeblod", "Blod", "Vev"),   ##### !!!! ASP : Blod as well
  month = unlist(DATA$sampling.months))

## SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month > 12, rovbaseObs.sp$year, rovbaseObs.sp$year-1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset, ] 

## [PD] ASP FIXED FILTERING
## REMOVE SAMPLES THAT WERE SUCCESSFULLY GENOTYPED & FROM THE FOCAL SPECIES 
subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$Individ) 
rovbaseObs.sp <- rovbaseObs.sp[!subset, ] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace, ] 

## Check if the filter is correct  
tmp <- rovbaseObs.sp[!is.na(rovbaseObs.sp$Individ), ]
table(tmp$`Art (Analyse)`) # Correct, all wolverine samples have NAs (not IDs)

## RASTERIZE 
r.detector <- aggregate(habitat.subdetectors, fact=(DETECTORS$detResolution/DETECTORS$detSubResolution))
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

##-- Plot check
if(plot.check){
  # pdf(file = file.path(working.dir, modelName, "mapStructuredOthers.pdf"))
  for(t in 1:nYears){ 
    year = years[t]
    tmpOthers <- myFilteredData.spOthers[myFilteredData.spOthers$Year%in%year, ]
    tmpStruct <- myFilteredData.spStructured[myFilteredData.spStructured$Year%in%year, ]
    
    par(mfrow=c(2,2),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]], main=paste(year,"\n Rovbase Samples Structured"), box=F, axes=F)
    plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    plot(r.OtherSamplesBinary[[t]],main=paste(year,"\n Rovbase Samples Opportunistic"), box=F, axes=F)
    plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
    
    plot(r.skandObsSamplesBinary[[t]], main=paste(year,"\n SkandObs Structured"), box=F, axes=F)
    plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    plot(r.skandObsSamplesBinary[[t]],main=paste(year,"\n SkandObs Opportunistic"), box=F, axes=F)
    plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
  }
  # dev.off()
}



## ------         3.2.6.3. COMBINE ROVBASE & SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][] > 1] <- 1
}

##-- Plot check
if(plot.check){
  for(t in 1:nYears){
    par(mfrow=c(1,3),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]],main=years[t])
    plot(r.skandObsSamplesBinary[[t]])
    plot(r.SkandObsOtherSamplesBinary[[t]])
  }  
}


## ------         3.2.6.4. SMOOTH THE BINARY MAP ------

## we tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
habOwin <- as.owin(as.vector(extent(r.detector)))
y <- 2020
cutoff <- 1
ds.list <- lapply(years,function(y){
  ## ROVBASE DATA 
  pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y,]
  ## SKANDOBS
  pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y,] )
  ## SMOOTH AND RASTERIZE
  p <- ppp(pts[,1], pts[,2], window=habOwin)
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
names(ds.brick) <- years

##-- Plot check
if(plot.check){
  par(mfrow = c(1,3))
  plot(r.SkandObsOtherSamplesBinary[[t]], main = "Raw Binary", axes = F, box = F)
  plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
  plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
}



## ------         3.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------

## IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.sp$alive[which(myFilteredData.sp$alive$DNAID %in% HairTrapSamples$DNAID), ]

## MANUALLY FIND THE HAIR SMAPLES & COLOR THE CELL
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr[i])#[missing an i index
  whereHair <- raster::extract( r.SkandObsOtherSamplesBinary[[t]],
                                tmpHair[tmpHair$Year %in% tmpyr[i],],#issue here as well
                                cellnumbers = T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[ ,1]] <- 1
  plot(r.SkandObsOtherSamplesBinary[[t]])
  plot(tmpHair[tmpHair$Year %in% tmpyr[i],]$geometry, add = T, col = "red")
}



## ------         3.2.6.7. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract( r.SkandObsOtherSamplesBinary,
                                                myDetectors$main.detector.sp)
colSums(detOtherSamples)



## ------       3.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

detCovs <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],2))
detCovs[,,1] <- detTracks
detCovs[,,2] <- detSnow

detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2],3))
detCovsOth[,,1] <- detSnow
detCovsOth[,,2] <- matrix(detRoads, length(detRoads), nYears)
detCovsOth[,,3] <- detOtherSamples

## CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}

##-- Plot check
if(plot.check){
  tmp <- detectorGrid.r
  par(mfrow = c(2,5), mar = c(0,0,0,0))
  max <- max(detCovsOth[ , ,2])
  cuts <- seq(0,max,length.out = 100) # set breaks
  col <- rev(terrain.colors(100))
  for(t in 1:nYears){
    plot(detectorGrid.r, col = c(grey(0.2),grey(0.8)), axes = F, legend = F, box = F)
    tmp[!is.na(detectorGrid.r)] <- detCovsOth[,t,2]
    plot(tmp, axes = F, legend = F, box = F, breaks = cuts, col = col, add = T)
  }
  #dev.off()
  
  # pdf(file = file.path(working.dir, modelName, "Detections over space and time.pdf"))
  for(t in 1:nYears){
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
  # dev.off()
}



## ------   4. GENERATE y DETECTION ARRAYS ------

## ------     4.1. ASSIGN SAMPLES TO DETECTORS ------

## ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$alive,                
  #myDetectors = myDetectors.dead$main.detector.sp, ==> FIXED myDetectors [PD]
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = DETECTORS$detResolution)

## STRUCTURED
myData.aliveStruc <- AssignDetectors_v3sf(
  myData = myFilteredData.spStructured,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = DETECTORS$detResolution)

## OTHERS
myData.aliveOthers <- AssignDetectors_v3sf(
  myData = myFilteredData.spOthers,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = DETECTORS$detResolution)


## MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET ASSIGNED TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING
## FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE OF NORRBOTTEN
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)

whichdets <- which(!myData.alive$myData.sp$Year %in% yearsSampledNorrb &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)
whichdetsStruc <- which(!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveStruc$myData.sp$Detector %in% detsNorrbotten)
whichdetsOther <- which(!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveOthers$myData.sp$Detector %in% detsNorrbotten)


##-- [PD] fixed assignment to sub-detectors in Norrbotten
##-- Identify sub-detectors inside Norrbotten
subDetsNorrbotten <- which( myDetectors$detector.sp$main.cell.id %in% 
                              myDetectors$main.detector.sp$main.cell.id[detsNorrbotten])

##-- Loop over flagged detections and assign them to closest sub-detector outside Norrbotten
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

plot(myDetectors$main.detector.sp$geometry)
plot(tmp$geometry,add=T,col="red")
plot(myDetectors$main.detector.sp[myDetectors$main.detector.sp$main.cell.id%in% myData.alive$myData.sp$Detector[whichdets[i]],]$geometry ,add=T,col="red")
plot(myDetectors$main.detector.sp[myData.alive$myData.sp$Detector[whichdets[i]], ]$geometry ,add=T,col="blue")


## STRUCTURED
##-- Loop over flagged detections and assign them to closest sub-detector outside Norrbotten
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

## SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
sum(myData.aliveOthers$myData.sp$Detector[!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
sum(myData.aliveStruc$myData.sp$Detector[!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)


## DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf(
  myData = myFilteredData.sp$dead.recovery,
  myDetectors = myDetectors.dead$main.detector.sp,
  radius = DETECTORS$detResolution)




## ------     4.2. SAVE PREPARED DATA ------

# save( myData.alive, myData.aliveStruc, myData.aliveOthers, myData.dead,
#       file = file.path(working.dir, modelName, "myFilteredData.RData"))

###[CM] NEEDS THIS OTHERWISE THE LOOPS OVERWRITES THE SEX
myData.aliveALL <- myData.alive
myData.aliveStrucALL <- myData.aliveStruc
myData.aliveOthersALL <- myData.aliveOthers
myData.deadALL <- myData.dead


## ------     4.3. GENERATE DETECTION HISTORY (FOR BOTH SEXES) ------

for(thisSex in c("Hann","Hunn")){
  
  message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
  
  ## ------    4.4. FILTER DATA BY SEX ------
  
  # load(file.path(working.dir, modelName, "myFilteredData.RData"))
  
  dim(myData.alive$myData.sp)
  myData.alive$myData.sp <- myData.aliveALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  dim(myData.alive$myData.sp)
  
  dim(myData.aliveStruc$myData.sp)
  myData.aliveStruc$myData.sp <- myData.aliveStrucALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  dim(myData.aliveStruc$myData.sp)
  
  dim(myData.aliveOthers$myData.sp)
  myData.aliveOthers$myData.sp <- myData.aliveOthersALL$myData.sp %>%
    dplyr::filter(Sex %in% thisSex)
  dim(myData.aliveOthers$myData.sp)
  
  dim(myData.dead)
  myData.dead <- myData.deadALL %>%
    dplyr::filter(Sex %in% thisSex)
  dim(myData.dead)
  
  # do this for the dead recoveries
  myFullData.spDeadsex <-  myFullData.spSaved$dead.recovery %>%
    dplyr::filter(Sex %in% thisSex)
  
  ##-- EXPORT THE DATA 
  if( thisSex%in% "Hann"){
    assign("myData.aliveM", myData.alive)
    assign("myData.aliveOthersM", myData.aliveOthers)
    assign("myData.aliveStrucM", myData.aliveStruc)
    assign("myData.deadM", myData.dead)
    assign("myFullData.spDeadM", myFullData.spDeadsex)
    
    save(myData.aliveM,
         myData.aliveOthersM,
         myData.aliveStrucM,
         myData.deadM,
         myFullData.spDeadM,
         file = file.path(working.dir, modelName,thisSex,
                          paste0(modelName, "_NGSData.RData")))
    
  } else {
    assign("myData.aliveF", myData.alive)
    assign("myData.aliveOthersF", myData.aliveOthers)
    assign("myData.aliveStrucF", myData.aliveStruc)
    assign("myData.deadF", myData.dead)
    assign("myFullData.spDeadF", myFullData.spDeadsex)
    
    save(myData.aliveF,
         myData.aliveOthersF, 
         myData.aliveStrucF,
         myData.deadF,
         myFullData.spDeadF,
         file = file.path(working.dir, modelName,thisSex,
                          paste0(modelName, "_NGSData.RData")))
  }
  
  
  
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
  
  dim(y.ar.ALIVEOth)
  dim(y.ar.ALIVEStruc)
  dim(y.ar.ALIVE)
  ## RESIZE DETECTION ARRAYS TO MAKE SURE THEY HAVE THE SAME DIMENSIONS
  y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
  y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0
  ## FILL IN THE Y ARRAYS 
  y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]],,] <- y.ar.ALIVEOth
  y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]],,] <- y.ar.ALIVEStruc
  
  ## PROJECT THE DEATH TO THE NEXT OCCASION
  y.ar.DEADProjected <- y.ar$y.ar2 
  y.ar.DEADProjected[] <- 0
  for(t in 2:nYears){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}
  
  ## TURN INTO ANNUAL DEAD RECOVERY MATRIX
  y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
  y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
  y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
  dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
  dim(y.ar.DEAD)
  y.ar.DEAD[y.ar.DEAD > 0] <- 1
  dim(y.ar$y.ar2)
  
  
  
  ## ------    4.6. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------
  
  distances <- list()
  for(t in 1:nYears){
    print(paste("------ ", t ," -------", sep = "" ))
    distances[[t]] <- CheckDistanceDetectionsV2sf( 
      y = y.ar.ALIVE[,,t], 
      detector.xy = detector.xy, 
      max.distance = DETECTIONS$maxDetDist,
      method = "pairwise",
      plot.check = F)
    
    # ## PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
    # if(plot.check){
    #   par(mfrow = c(1,1))
    #   if(sum(distances[[t]]$y.flagged) > 0){
    #     affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
    #     count <- 0
    #     for(i in affected.ids){
    #       count <- count+1
    #       plot(st_geometry(myStudyArea), main = paste("t: ",t,"     i: ", names(affected.ids)[count], sep = ""))
    #       scalebar(2*DETECTIONS$maxDetDist, xy = c(800000,6700000), type = "bar", divs = 2, below = "km",
    #                label = c(0, DETECTIONS$maxDetDist/1000, DETECTIONS$maxDetDist/500), cex = 0.8, adj = c(0.5,-0.9))
    #       plot(st_geometry(COUNTRIES), add = T)
    #       plot(st_geometry(myDetectors$main.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
    #       
    #       tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar.ALIVE)[[1]][i] &
    #                                        myFilteredData.sp$alive$Year == years[t], ]
    #       tmp <- tmp[order(tmp$Date), ]
    #       tmp.xy <- st_coordinates(tmp)
    #       n.det <- nrow(tmp.xy)
    #       
    #       plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
    #       arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
    #              x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
    #              length = 0.1, lwd = 1)
    #       plot(st_geometry(myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
    #       
    #       tmp2 <- myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
    #       plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
    #     }#i
    #   }#if
    # }#if
    
    ## REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
    y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
    y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
    
    ## REMOVE DETECTIONS ALSO IN MYDATA TO RUN GETSINITS
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
  
  ## Plot an image of the matrix
  if(plot.check){
    par(mfrow = c(1,1))
    barplot(colSums(apply(y.ar.ALIVE, c(1,3),function(x)any(x>0))))
    barplot(colSums(already.detected), add = TRUE, col = "gray40")
    legend(x = 0, y = 250, legend = c("newly Det", "already Det"),
           fill = c("gray80", "gray40"))
  }
  
  
  
  ## ------      4.7.2. AGE ------
  
  # min.age <- age <- precapture <- matrix(NA, dim(y.ar.ALIVE)[1], dim(y.ar.ALIVE)[3], dimnames = list(y.ar$Id.vector,years))
  # 
  # temp <- apply(y.ar.ALIVE, c(1,3), sum)
  # year.first.capture <- apply(temp, 1, function(x)min(years[which(x>0)]))
  # year.first.capture[is.infinite(year.first.capture)] <- NA
  # names(year.first.capture) <- y.ar$Id.vector
  # 
  # for(i in y.ar$Id.vector){
  #   this.set <- myData.dead[myData.dead$Id == i, ]
  #   year.dead <- myData.dead$Death[myData.dead$Id == i]
  #   year.first.captured <- year.first.capture[i]
  #   precapture[i,] <- as.numeric(years < year.first.captured)
  #   if(all(is.na(precapture[i,])))precapture[i,] <- 1
  #   latest.recruitment.year <- min(year.dead,year.first.captured, na.rm = TRUE) 
  #   
  #   try({
  #     min.age[i,] <- years-latest.recruitment.year
  #   },silent = TRUE)
  #   
  #   try({
  #     birth.year <- this.set$Death-this.set$min.age
  #     if(birth.year<latest.recruitment.year) min.age[i,] <- years-birth.year 
  #   },silent = TRUE)
  #   
  #   try({
  #     birth.year <- this.set$Death - this.set$age
  #     age[i,] <- years-birth.year
  #   }, silent = TRUE)
  # }
  # image(t(min.age))
  # image(t(age))
  
  
  
  ## ------    4.8. MAKE AUGMENTATION ------
  
  ## DATA ARRAYS
  y.alive <- MakeAugmentation(
    y = y.ar.ALIVE,
    aug.factor = DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.aliveStructured <- MakeAugmentation(
    y = y.ar.ALIVEStructured,
    aug.factor = DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.aliveOthers <- MakeAugmentation(
    y = y.ar.ALIVEOthers,
    aug.factor = DETECTIONS$aug.factor,
    replace.value = 0)
  
  y.dead <- MakeAugmentation( 
    y = y.ar.DEAD,
    aug.factor = DETECTIONS$aug.factor,
    replace.value = 0)
  
  ## INDIVIDUAL COVARIATES
  already.detected <- MakeAugmentation( 
    y = already.detected,
    aug.factor = DETECTIONS$aug.factor,
    replace.value = 0)
  
  
  
  ##----------------------------------------------------------------------------
  
  ## ------ III.MODEL SETTING & RUNNING ------- 
  
  ## ------ 1. NIMBLE MODEL DEFINITION ------
  
  modelCode <- nimbleCode({
    
    ##-----------------------------## 
    ##------ SPATIAL PROCESS ------##  
    ##-----------------------------##  
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
    }#i
    
    for(t in 2:n.years){
      for(i in 1:n.individuals){
        sxy[i, 1:2, t] ~ dbernppACmovement_exp(
          lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
          upperCoords = upperHabCoords[1:numHabWindows, 1:2],
          s = sxy[i, 1:2, t-1],
          lambda = lambda,
          baseIntensities = habIntensity[1:numHabWindows],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max,
          numWindows = numHabWindows)
      }#i  
    }#t
    
    
    ##-------------------------------## 
    ##----- DEMOGRAPHIC PROCESS -----## 
    ##-------------------------------##    
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
    
    
    pResponse ~ dunif(0, 1)
    
    for(i in 1:n.individuals){ 
      detResponse[i,1] ~ dbern(pResponse)
      
      z[i,1] ~ dcat(omeg1[1:2]) 
      for(t in 1:n.years1){
        z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
      }#i 								
    }#t 
    
    
    ##-----------------------------##
    ##----- DETECTION PROCESS -----## 
    ##-----------------------------##
    for(t in 1:n.years){
      sigma[t] ~ dunif(0,4)
      for(c in 1:n.covs){
        betaCovs[c,t] ~ dunif(-5,5)
      }
      
      for(c in 1:n.covsOth){
        betaCovsOth[c,t] ~ dunif(-5,5)
      }
      
      betaResponse[t] ~ dunif(-5,5)
      betaResponseOth[t] ~ dunif(-5,5)
    }
    
    for(c in 1:n.counties){
      for(t in 1:n.years){
        p01[c,t] ~ dunif(0,1)
        p0[c,t] <- p01[c,t] *countyToggle[c,t]## toggle counties
      }#t
    }#c  
    
    for(c in 1:n.countries){
      for(t in 1:n.years){
        p01Oth[c,t] ~ dunif(0,1)
        p0Oth[c,t] <- p01Oth[c,t] *countyToggleOth[c,t]## toggle countries
      }#t
    }#c  
    
    for(t in 1:n.years){
      for(i in 1:n.individuals){
        y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCovResponse( 
          sxy = sxy[i,1:2,t],
          sigma = sigma[t],
          nbDetections = nbDetections[i,t],
          yDets = yDets[i,1:nMaxDetectors,t],
          detector.xy = detector.xy[1:n.detectors,1:2],
          trials = trials[1:n.detectors],
          detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
          nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
          ResizeFactor = ResizeFactor,
          maxNBDets = maxNBDets,
          habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
          indicator = isAlive[i,t],
          p0State = p0[1:n.counties,t],
          detCountries = detCounties[1:n.detectors],
          detCov = detCovs[1:n.detectors,t,1:n.covs],
          betaCov = betaCovs[1:n.covs,t],
          BetaResponse = betaResponse[t],
          detResponse = detResponse[i,t])
        
        y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbin_LESS_Cached_MultipleCovResponse(
          sxy = sxy[i,1:2,t],
          sigma = sigma[t],
          nbDetections = nbDetectionsOth[i,t],
          yDets = yDetsOth[i,1:nMaxDetectorsOth,t],
          detector.xy = detector.xy[1:n.detectors,1:2],
          trials = trials[1:n.detectors],
          detectorIndex = detectorIndex[1:n.cellsSparse,1:maxNBDets],
          nDetectorsLESS = nDetectorsLESS[1:n.cellsSparse],
          ResizeFactor = ResizeFactor,
          maxNBDets = maxNBDets,
          habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
          indicator = isAlive[i,t],
          p0State = p0Oth[1:n.countries,t],
          detCountries = detCountries[1:n.detectors,t],
          detCov = detCovsOth[1:n.detectors,t,1:n.covsOth],
          betaCov = betaCovsOth[1:n.covsOth,t],
          BetaResponse = betaResponseOth[t],
          detResponse = detResponse[i,t])
      }#i
    }#t
    
    
    ##----------------------------------------## 
    ##---------- DERIVED PARAMETERS ----------##
    ##----------------------------------------##
    for(i in 1:n.individuals){ 
      isAlive[i,1] <- (z[i,1] == 2) 
      for(t in 1:n.years1){
        isAlive[i,t+1] <- (z[i,t+1] == 2) 
      }
    }
    for(t in 1:n.years){
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
                        n.countries = max(detCountries)+1, # + 1 for Norrbotten
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
  
  ## LATENT VARIABLE DET RESPONSE
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
    scaleToGrid =T )$coordsDataScaled
  ScaledUpperCoords <- scaleCoordsToHabitatGrid(
    coordsData = upperHabCoords,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid =T )$coordsDataScaled
  
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
  maxDistReCalc <- 2.1 * DETECTIONS$maxDetDist 
  
  DetectorIndexLESS <- GetDetectorIndexLESS(
    habitat.mx = myHabitat$habitat.mx,
    detectors.xy = nimData$detector.xy,
    maxDist = maxDistReCalc/res(myHabitat$habitat.r)[1],
    ResizeFactor = 1,
    plot.check = TRUE)
  
  DetectorIndexLESS$nDetectorsLESS
  
  dim(DetectorIndexLESS$detectorIndex)
  
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
  myData.deadProj <- myData.dead[,c("Id","Year")]
  myData.deadProj$Year <- myData.deadProj$Year +1#project dead reco to the next year
  ## Remove dead reco occuring the last year (not used)
  myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]
  
  ## Create a data.frame with all detection of all Individuals detected
  AllDets <- rbind(myData.alive$myData.sp [,c("Id","Year")],
                   myData.deadProj[ ,c("Id","Year")])
  AllDetections <- as.data.frame(AllDets)
  AllDetsxy <- st_coordinates(AllDets) 
  colnames(AllDetsxy) <- c("x","y")
  AllDetsxyscaled <- scaleCoordsToHabitatGrid(
    coordsData = AllDetsxy,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid =T )$coordsDataScaled
  AllDetections <- cbind(AllDetections, AllDetsxyscaled)
  
  idAugmented <- which(rownames(z) %in%"Augmented")
  Id.vector <- y.ar$Id.vector
  lowerCoords = nimData$lowerHabCoords
  upperCoords = nimData$upperHabCoords
  habitatGrid = nimData$habitatGrid
  
  ## GENERATE sxy INITIAL VALUES
  sxy.init <- getSInits(
    AllDetections = AllDetections,
    Id.vector = Id.vector,
    idAugmented = idAugmented,
    lowerCoords = lowerCoords,
    upperCoords = upperCoords,
    habitatGrid = habitatGrid,
    intensity = NULL,
    sd = 4,
    movementMethod = "dbernppACmovement_normal")
  
  ## RESCALE sxy INITIAL VALUES 
  sxy.initscaled <- scaleCoordsToHabitatGrid(
    coordsData = sxy.init,
    coordsHabitatGridCenter = myHabitat$habitat.xy,
    scaleToGrid =F )$coordsDataScaled
  
  
  
  ## ------ 8. calculate realized phi ------
  
  ## Initialize objects 
  recruit <- 0
  z_caculate <- nimData$z#[,1:(nYears-1)]
  z_caculate[is.na(z_caculate)] <- 0
  lev <- levels(habitatRasterResolution$'2.5km'$Countries)
  
  ## EXTRACT LOCATION BASED ON INITIAL AC
  countryId <- list()
  for(t in 1:dim(z_caculate)[2]){
    tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
    countryId[[t]] <- raster::extract( habitatRasterResolution$'2.5km'$Countries ,tmp,sparse = F)
  }
  
  phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=nYears-1,ncol=length(lev[[1]]$ID))
  colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
    colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'2.5km'$Countries,lev[[1]]$ID)[,1]
  for(c in 1:ncol(phi)){
    for(t in 2:dim(z_caculate)[2]){
      #phi
      alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
      phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
      #culled
      #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
      #recru
      notentered <- which(z_caculate[,t-1] == 0)
      recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                                countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
      recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                              countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
                                                                       countryId[[t-1]] %in% c)
    }
  }
  phi <- phi[,c(2,4)]        # overall phi
  recruit <- recruit[,c(2,4)]        # overall phi
  recruitnb <- recruitnb[,c(2,4)]        # overall phi
  
  
  ##-- Plot check
  if(plot.check){
    # pdf(file = file.path(working.dir, modelName, "realizedPhiCountry.pdf"))
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z")
    axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)]+1)
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(phi)){
      points(phi[,c]~yr,pch=16,type="b", col=c)
    }
    legend("bottomright",colnames(phi),col=c(1:4),pch=16)
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruit from z")
    axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)]+1)
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(recruit)){
      points(recruit[,c]~yr,pch=16,type="b", col=c)
    }
    legend("topright",colnames(recruit),col=c(1:4),pch=16)
    dev.off()
  }
  
  
  lev <- levels(habitatRasterResolution$'10km'$Counties)
  countryId <- list()
  for(t in 1:dim(z)[2]){
    tmp <- st_as_sf(data.frame(sxy.initscaled[,,t]), coords = c("x", "y"))
    countryId[[t]] <- raster::extract( habitatRasterResolution$'10km'$Counties ,tmp,sparse = F)
  }
  
  phi <- phiind1 <- culled <- recruit <- recruitnb<- matrix(0,nrow=nYears-1,ncol=length(lev[[1]]$ID))
  colnames(phi) <- colnames(phiind1) <- colnames(culled) <- 
    colnames(recruitnb) <-  colnames(recruit)  <- factorValues(habitatRasterResolution$'10km'$Counties,lev[[1]]$ID)[,1]
  
  for(c in 1:ncol(phi)){
    for(t in 2:dim(z)[2]){
      #phi
      alivet <- which(z_caculate[,t-1] %in% c(2) & countryId[[t-1]] %in% c )
      phi[t-1, c] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
      #culled
      #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
      #recru
      notentered <- which(z_caculate[,t-1] == 0)
      recruitnb[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                                countryId[[t]][notentered] %in% c)#/sum(z_caculate[,t-1] %in% c(2))
      recruit[t-1,c] <- sum(z_caculate[notentered,t] %in% c(2) & 
                              countryId[[t]][notentered] %in% c)/sum(z_caculate[,t-1] %in% c(2) & 
                                                                       countryId[[t-1]] %in% c)
    }
  }
  phi <- phi[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
  
  phiNOR <- phi[,c(8,9,10,11,12,13)]
  phiSWE <- phi[,-c(8,9,10,11,12,13,14,15,16)]
  
  recruitnb <- recruitnb[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
  recruitnbNOR <- recruitnb[,c(8,9,10,11,12,13)]
  recruitnbSWE <- recruitnb[,-c(8,9,10,11,12,13,14,15,16)]
  
  recruit <- recruit[,c(13,16,17,18,19,20,21,27,28,30,31,32,33,35,36,38)]        # overall phi
  recruitNOR <- recruit[,c(8,9,10,11,12,13)]
  recruitSWE <- recruit[,-c(8,9,10,11,12,13,14,15,16)]
  
  
  ##-- Plot check
  if(plot.check){
    # pdf(file = file.path(working.dir, modelName, "realizedPhiCounties.pdf"),
    # width = 11, height = 6)
    # PHI
    ## NORWAY
    par(mfrow = c(1,2))
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Norway")
    axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(phiNOR)){
      points(phiNOR[,c]~yr,pch=16,type="b", col=c)
    }
    legend("bottomleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)
    
    ## SWEDEN
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Sweden")
    axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(phiSWE)){
      points(phiSWE[,c]~yr,pch=16,type="b", col=c)
    }
    legend("bottomleft",colnames(phiSWE),col=c(1:ncol(phiSWE)),pch=16)
    
    
    # RECRUITS
    ## NORWAY
    par(mfrow = c(1,2))
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", 
         ylab = "Realized recruitment from z", main="Norway")
    axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(recruitNOR)){
      points(recruitNOR[,c]~yr,pch=16,type="b", col=c)
    }
    legend("topleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)
    
    ## SWEDEN
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruitment from z", main="Sweden")
    axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
    yr <- c(1:(nYears-1))
    for(c in 1:ncol(recruitSWE)){
      points(recruitSWE[,c]~yr,pch=16,type="b", col=c)
    }
    legend("topleft",colnames(recruitSWE),col=c(1:ncol(phiSWE)),pch=16)
    dev.off()
  }
  
  
  ## prop detected vs Alive in z
  propDet <- 0
  for(t in 1:nYears){
    whichdets <- unique(c(which(nimData$nbDetections[,t]>0),
                          which(nimData$nbDetectionsOth[,t]>0)))
    whichAlive <- which(nimData$z[,t]%in%2)
    propDet[t] <- length(whichdets)/length(whichAlive)
  }
  
  
  
  ## ------ 9. LIST NIMBLE INITS & SAVE NIMBLE INPUTS ------
  
  for(c in 1:4){
    
    ## ------  9.1. LIST NIMBLE INITS ------
    nimInits <- list( "sxy" = sxy.init,
                      "dmean" = runif(1,0,10),
                      "z" = z.init,
                      "omeg1" = c(0.5,0.5),
                      "gamma" = runif(dim(y.alive)[3]-1,0,1),
                      "p01" = array(runif(18,0,0.2), c(nimConstants$n.counties,dim(y.alive)[3])),
                      "p01Oth" = array(runif(18,0,0.2), c(nimConstants$n.countries+1,dim(y.alive)[3])),
                      "sigma" = runif(nYears,1,4),
                      "betaDens" = runif(1,-0.1,0.1),
                      "betaCovs" = array( runif(dim(detCovs)[3],-0.1,0.1),c(dim(detCovsOth)[3],nYears)),
                      "betaCovsOth" = array( runif(dim(detCovsOth)[3],-0.1,0.1),c(dim(detCovsOth)[3],nYears)),
                      "betaResponseOth" = runif(dim(y.alive)[3], -0.1, 0.1),
                      "betaResponse" = runif(dim(y.alive)[3], -0.1, 0.1),
                      "detResponse" = InitsDetResponse,
                      "pResponse"  = runif(1, 0.4, 0.5),
                      "phi" = runif(dim(y.alive)[3]-1,0.1,0.3))
    
    ##TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
    ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
    
    idDEtected <- which(!rownames(z) %in%"Augmented")
    for(i in 1:length(idDEtected)){
      for(t in 1:nimConstants$n.years){
        if(!is.na(nimInits$sxy[i,1,t])){ SXY <- nimInits$sxy[i, ,t] } else {SXY <- nimData$sxy[i, ,t]}
        sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1, trunc(SXY[1]/nimConstants$ResizeFactor)+1]
        DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse,1:nimConstants$maxNBDets] 
        DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
        index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
        
        ## GET NECESSARY INFO 
        n.detectors <- length(index)
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
    
    
    ##-- Plot check
    plot(nimData$detector.xy[,2]~nimData$detector.xy[,1])
    points(nimData$detector.xy[index,2]~nimData$detector.xy[index,1], col="red")
    points(SXY[2]~SXY[1], col="blue", pch=16)
    points(nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],2]~
             nimData$detector.xy[YDET[1:nimData$nbDetections[i, t]],1], col="green", pch=16)
    points(nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],2]~
             nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],1], col="purple", pch=16)
    
    plot(st_geometry(COUNTRIES))
    tmp <- myData.alive$myData.sp[myData.alive$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] & myData.alive$myData.sp$Year %in% years[t],]
    # tmp <- myData.aliveOthers$myData.sp[myData.aliveOthers$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] &
    #                                          myData.aliveOthers$myData.sp$Year %in% years[t],]
    # tmp <- myData.aliveStruc$myData.sp[myData.aliveStruc$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] &
    #                                      myData.aliveStruc$myData.sp$Year %in% years[t],]
    # 
    plot(st_geometry(tmp),col="red",add=T)
    
    
    ## An extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
    nimInits$sxy <- round(nimInits$sxy, 5)
    
    ## CHECK WHERE IS NORRBOTTEN. IT IS ON THE 5TH INDEX
    plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
    plot(st_geometry(myDetectors$main.detector.sp[detCounties%in% c(1),]), col = myCol[5], pch = 16, cex = 0.8,add=T)
    
    nimConstants$countyToggle <- nimInits$p01
    nimConstants$countyToggle[] <- 1
    
    yearsNotSampled <- which(!years%in% yearsSampledNorrb)
    for(t in yearsNotSampled){
      nimConstants$countyToggle[1,t] <- 0
    }
    
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
    
    nimData$detCountries <-  detCountriesNorb
    
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
         file = file.path(working.dir, modelName, thisSex,
                          paste0(modelName, thisSex,"_Chain", c, ".RData")))
  }#c
  
}#thisSex




## ------   10. SAVE NECESSARY OBJECTS ------

# load(file.path(working.dir, modelName, "myFilteredData.RData"))

myHabitat.list <- myHabitat
myDetectors <- myDetectors
COUNTRIES <- COUNTRIES
COUNTIES <- COUNTIES
myStudyArea.poly <- myStudyArea
COMMUNES <- COMMUNES
myFilteredData.sp <- myFilteredData.sp
myFullData.sp <- myFullData.sp
COUNTIES_AGGREGATED <- COUNTIES_AGGREGATED
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATEDSubset

save(myHabitat.list, myDetectors, COUNTRIES, myStudyArea.poly,
     COMMUNES, COUNTIES_AGGREGATEDSubset,
     myFilteredData.sp, myFullData.sp, COUNTIES_AGGREGATED,
     file = file.path(working.dir, modelName, "NecessaryObjects.RData" ))



##------------------------------------------------------------------------------

## ------ III. MAKE THE SINGLE SEASON SCR MODEL ------

## -----  1. SCR MODEL CODE ------

modelCode1 <- nimbleCode({
  

  ##------ SPATIAL PROCESS ------ 
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
  

  ##----- DEMOGRAPHIC PROCESS -----
  pResponse ~ dunif(0, 1)
  psi ~ dunif(0, 1)
  for(i in 1:n.individuals){ 
    detResponse[i] ~ dbern(pResponse)
    z[i] ~ dbern(psi)	
  }#t 
  
  ##----- DETECTION PROCESS -----
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
  

  ##---------- DERIVED PARAMETERS --------
  N <- sum(z[1:n.individuals])
})



## ------- 2. SCRize NIMBLE INPUT DATA ------

for(ch in 1:4){
  for(t in 1:nYears){ 
    for(thisSex in c("Hann","Hunn")){
      
      load( file.path(working.dir, modelName, thisSex,
                      paste0(modelName, thisSex,"_Chain", ch, ".RData")))
      
      ## GET WHICH INDIVIDUAL IS DETECTED       
      detectedStruc <- apply(nimData$nbDetections,2,function(x) x>0)
      detectedOth <- apply(nimData$nbDetectionsOth,2,function(x) x>0)
      detected <- detectedOth + detectedStruc
      detected <- detected>0
      
      ## GET SUM OF INDIVIDUALS DETECTED AND DECIDE HOW MUCH YOU WISH TO AUGMENT. hERE I CHOSE 2
      sumDets <- sum(detected[ ,t])*3## decide the augmentation factor 
      
      ## SUBSET AND AUGMENT ALL OBJECTS BASED ON WHETHER INDIVIDUALS WHERE DETECTED OR NOT
      ###nimData
      ## y.alive
      nimData$y.alive <- nimData$y.alive[detected[ ,t], ,t]  
      nimData$y.alive <- rbind( nimData$y.alive,
                                matrix( 0,
                                        nrow = sumDets,
                                        nimConstants$nMaxDetectors))
      
      nimData$y.aliveOth <- nimData$y.aliveOth[detected[ ,t], ,t]  
      nimData$y.aliveOth <- rbind( nimData$y.aliveOth,
                                   matrix( 0,
                                           nrow = sumDets,
                                           nimConstants$nMaxDetectorsOth))
      
      ## z
      nimData$z  <-  nimData$z[detected[,t],t]  
      nimData$z[nimData$z %in% c(2)] <- 1 # ALIVE IDS BECOMES 1
      nimData$z <- c(nimData$z, rep(NA,sumDets))
      
      ## SXY 
      nimData$sxy <- NULL
      
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
      nimData$detResponse <- nimData$detResponse[detected[,t],t]
      nimData$detResponse  <- c(nimData$detResponse, rep(NA,sumDets))## HERE IT IS ASSUMING IT IS A LATENT INDIVIDUAL COVARIATE
      
      ##detCovs
      nimData$detCovs <- nimData$detCovs[,t,]
      nimData$detCovsOth <- nimData$detCovsOth[,t,]
      
      ## density
      nimData$denCounts <- nimData$denCounts[,1]
      #countyToggle to toggle off norbotten
      nimConstants$countyToggle <- nimConstants$countyToggle[,t]
      nimConstants$countyToggleOth <- nimConstants$countyToggleOth[,t]
      
      
      ### nimInits
      ## z
      nimInits$z <- nimInits$z[detected[,t],t]  
      nimInits$z <- c(nimInits$z, rbinom(sumDets,1,0.5))
      
      ## detResponse
      ## HERE IT IS TREATED AS A LATENT COVARIATE
      nimInits$detResponse  <- c(rep(NA,sum(detected[,t])), rbinom(sumDets,1,0.5))
      
      ## SXY 
      nimInits$sxy <- nimInits$sxy[detected[,t],,t]  
      # GIVE ACS FROM DETECTED INDIVIDUALS TO AUGMENTED IDS. 
      nimInits$sxy <- rbind(nimInits$sxy, nimInits$sxy[sample(nimInits$sxy, sumDets, replace = T),])
      
      ## p0
      nimInits$p01 <-  nimInits$p01[,t]
      nimInits$p01Oth <-  nimInits$p01Oth[,t]
      
      ## psi
      nimInits$psi <-  runif(1,0.4,0.6)
      
      ## sigma
      nimInits$sigma <-  runif(1,1,2)
      
      ## betaCovs
      nimInits$betaCovs <- nimInits$betaCovs[,t]
      nimInits$betaCovsOth <- nimInits$betaCovsOth[,t]
      
      ## betaResponse
      nimInits$betaResponse <- nimInits$betaResponse[t]
      nimInits$betaResponseOth <- nimInits$betaResponseOth[t]

      nimData$detCountries <- nimData$detCountries[ ,t]
      
      ## Nimble constants
      nimConstants$n.individuals <- nrow(nimInits$sxy)
      
      ## Nimble parameters
      nimParams <- c("N", "psi", "pResponse","p0Oth","betaCovsOth","betaResponseOth",
                     "p0", "sigma", "betaDens", "betaCovs","betaResponse","betaResponseOth")
      
      nimParams2 <- c("z", "sxy")
      
      modelCode <- modelCode1
      
      save(nimData,
           nimConstants,
           y.dead,
           nimParams,
           nimParams2,
           modelCode,
           nimInits,
           file = file.path(working.dir, modelName, thisSex, Snapshot,
                            paste0("Snap",years[t],"_",modelName, thisSex,"_Chain", ch, ".RData")))
    }#c
  }
}



##-- Test
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()



##------------------------------------------------------------------------------

## ------ IV. RUN NIMBLE MODEL ------

## ------   1 CONFIGURE NIMBLE MODEL ------
thisSex <- "Hunn"
for(thisSex in c("Hann","Hunn")){
  for(c in c(2,4)){
    # load(file.path( working.dir, modelName, thisSex,
    #                 paste0(modelName, thisSex,"_Chain", c, ".RData")))
    ptm <- proc.time()
    model <- nimbleModel( code = modelCode,
                          constants = nimConstants,
                          inits = nimInits,
                          data = nimData,
                          check = FALSE,
                          calculate = FALSE) 
    system.time(print(model$calculate()))#-740704.4
  }#c
}#thisSex

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



## ------   2. CHECK INITIAL LOG-LIKELIHOODS ------

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



## ------   3. RESUME MODEL CONFIGURATION ------

conf <- configureMCMC(model, monitors = nimParams, thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc



## ------   4. RUN NIMBLE MCMC IN SUCCESSIVE BITES ------

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

## ------ V. PROCESS SCR OUTPUT ------

## ------   1.LOAD AND SELECT DATA ------

## ------   1. HABITAT DATA ------

## ------     1.1.LOAD RAW SHAPEFILES ------

## AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATE <- COUNTIES
#COUNTIES_AGGREGATED <- gSimplify(COUNTIES_AGGREGATED,tol=500, topologyPreserve = TRUE)
COUNTIES_AGGREGATE$id <- 1:nrow(COUNTIES_AGGREGATE)
# ggplot(COUNTIES_AGGREGATE) +
#   geom_sf(aes(fill = id)) +
#   geom_sf_label(aes(label = id))

#[CM] adjust Counties aggregation
COUNTIES_AGGREGATE$id[c(24,3,15,9,14,38,40,21,27,37,31,26,34,5,8,12,36,13,7)] <- 3
COUNTIES_AGGREGATE$id[c(39,33,23,32,29,22,4,11,20,2,10,16,25,1)] <- 4
COUNTIES_AGGREGATE$id[c(19)] <- 1
COUNTIES_AGGREGATE$id[c(35)] <- 2
COUNTIES_AGGREGATE$id[c(17,28)] <- 5
COUNTIES_AGGREGATE$id[c(18)] <- 7
COUNTIES_AGGREGATE$id[c(30)] <- 6

COUNTIES_AGGREGATE <- COUNTIES_AGGREGATE %>% group_by(id) %>% summarize()
#COUNTIES_AGGREGATE <- aggregate(x = gBuffer(COUNTIES_AGGREGATE, width = 0, byid = T), by = "id")
COUNTIES_AGGREGATED <- st_simplify(COUNTIES_AGGREGATE,preserveTopology = T,dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id
ggplot(COUNTIES_AGGREGATED) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))



## ------     1.2.CREATE STUDY AREA POLYGON ------

if(!is.null(HABITAT$countries)){
  myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% HABITAT$countries, ]
  
  myBufferedArea <- st_buffer(st_as_sf(myStudyArea) ,dist = HABITAT$habBuffer)
  myBufferedArea$id <- 1
  myBufferedArea <- myBufferedArea %>% group_by(id) %>% summarize()
  myBufferedArea <- st_intersection(myBufferedArea, GLOBALMAP)
}

##-- Plot check
if(plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myBufferedArea), add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
  plot(st_geometry(myStudyArea), add = TRUE, col ="red")
}



## ------   2. LOAD NECESSARY OBJECTS ------

# LOAD OBJECTS
load(file.path(WD, modelName, "NecessaryObjects.RData" ))

#load the habitat 
load(paste(dir.dropbox,"/DATA/GISData/spatialDomain/Habitat20kmNewSweCounties.RData",sep=""))
load(paste(dir.dropbox,"/DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData",sep=""))

## LOAD NECESSARY OBJECTS
## INFILES
#FEMALES
load(file.path(WD, modelName,"Hunn",paste(modelName,"Hunn","_Chain","1.RData",sep="")))
nimDataF <- nimData

#MALES
load(file.path(WD, modelName,"Hann",paste(modelName,"Hann","_Chain","1.RData",sep="")))
nimDataM <- nimData

WDFigures <- file.path(WD, modelName,"FigureSnap")
WDTables <- file.path(WD, modelName,"TableSnap")

if(!dir.exists(file.path(WDFigures))){dir.create(WDFigures)}
if(!dir.exists(file.path(WDTables))){dir.create(WDTables)}



## ------   8. TABLES OF #NGS SAMPLES, #DEAD RECOVERIES & #IDs DETECTED ------

## ------     8.1 OVERALL NUMBERS ------

load(file.path(WD, modelName, "Hunn",paste(modelName,"_NGSData.RData", sep="")))
load(file.path(WD, modelName, "Hann",paste(modelName,"_NGSData.RData", sep="")))

# load(file.path("C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2022", modelNameM, paste(modelNameM,"_NGSData.RData", sep="")))


## --- SOME TALLIES TO CHECK THINGS
# --- NGS
NGS <- rbind(myData.aliveF$myData.sp, myData.aliveM$myData.sp)
table(NGS$Year)
nrow(NGS)


NGSStructured <- rbind(myData.aliveStrucF$myData.sp, myData.aliveStrucM$myData.sp)
NGSOther <- rbind(myData.aliveOthersF$myData.sp, myData.aliveOthersM$myData.sp)
table(NGSStructured$Year)
nrow(NGSStructured)+
  nrow(NGSOther)


# NGS.all <- rbind(myFullData.spF$alive,myFullData.spM$alive)
# NGS.all <- NGS.all[NGS.all$Year%in%years, ]
# table(NGS.all$Year)
# length(NGS.all)

#### FOR REPORT SUMMARY
length(NGS$Id)
length(NGS$Id[NGS$Sex=="Hunn"])
length(NGS$Id[NGS$Sex=="Hann"])

length(NGS$Id[NGS$Country=="S"])/nrow(NGS)


length(unique(NGS$Id))
length(unique(NGS$Id[NGS$Sex=="Hunn"]))
length(unique(NGS$Id[NGS$Sex=="Hann"]))



#last year
length(NGS$Id[NGS$Year %in%tail(years, n=1)])
length(NGS$Id[NGS$Sex=="Hunn" & NGS$Year %in%tail(years, n=1)])
length(NGS$Id[NGS$Sex=="Hann"& NGS$Year %in%tail(years, n=1)])


###
length(NGSStructured$Id)
length(NGSStructured$Id[NGSStructured$Sex=="Hunn"])
length(NGSStructured$Id[NGSStructured$Sex=="Hann"])

length(NGS$Id[NGS$Country=="S"])/nrow(NGS)

length(NGSOther$Id)
length(NGSOther$Id[NGSOther$Sex=="Hunn"])
length(NGSOther$Id[NGSOther$Sex=="Hann"])



# --- DEAD RECOVERY
dead <- rbind(myFullData.spDeadF, myFullData.spDeadM)
table(dead$Year)
length(dead)
length(unique(dead$Id[dead$Sex=="Hunn"]))
length(unique(dead$Id[dead$Sex=="Hann"]))


###
tmpdead <- dead[dead$Year %in% c(2018:2023),]
tmpdead <- tmpdead[!duplicated(tmpdead$DNAID),]

table(tmpdead$Year,tmpdead$Sex)
table(tmpdead$Year,tmpdead$Month)
table(tmpdead$Year)
nrow(table(tmpdead$Year))
duplicated(tmpdead$DNAID)

# mapview::mapview(tmpdead)

tmpNGS <- NGS[NGS$Year %in% c(2018:2023),]
table(tmpNGS$Year,tmpNGS$Sex)
table(tmpNGS$Year,tmpNGS$Month)
table(tmpNGS$Year)

###HENRIK CHECK WITH PUBLIC SAMPLES
idPublic <-
  c(
    'D555438',
    'D556590',
    'D553322',
    'D555239',
    'D554783',
    'D558440',
    'D555845',
    'D555412',
    'D556343',
    'D556344',
    'D556347',
    'D558235',
    'D557101')

which(NGSStructured$DNAID%in%idPublic)
# 
# # NGSStructured[which(NGSStructured$DNAID%in%idPublic),]@data
# # NGSOther[which(NGSOther$DNAID%in%idPublic),]@data
# 
# 
# idPublic1 <-
#   c(
#     'D555438')
# NGSStructured[which(NGSStructured$DNAID%in%idPublic1),]@data
# NGSOther[which(NGSOther$DNAID%in%idPublic1),]
# # 
# mapview(
# list(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"),
#         NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]),
# layer.name = c("Franconian districts", "Franconian breweries")
# )
#         
# mapview(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"))+
#   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# 
# 
# mapview(as(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),"Spatial"))+
#   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# 
# st_distance(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),
#             st_as_sf(NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]))
# 
# 
# 
# 
# plot(TRACKS_YEAR[[9]][TRACKS_YEAR[[9]]$RovbaseID %in% "T471191",]$geometry)
# plot(TRACKSSimple_sf[[9]][TRACKSSimple_sf[[9]]$RovbaseID %in% "T471191",]$geometry)
# plot(TRACKS[TRACKS$RovbaseID %in% "T471191",]$geometry,col="red",add=T)
# points(NGSStructured[which(NGSStructured$DNAID%in%idPublic),],pch=16)



## ------     8.2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX ------

## ------       8.2.1 ALL------

NGSCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
#colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
colnames(NGSCountrySEX) <- unlist(lapply(YEARS, function(x) c(x[2],x[2])))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

NGSCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- nrow(temp[temp$Country %in% "N", ])
    NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S", ])
    NGSCountrySEX["Total",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))



# print(xtable(NGSCountrySEX, type = "latex",
#              align = paste(c("l",rep("c",ncol(NGSCountrySEX))),collapse = "")),
#       #scalebox = .8,
#       floating = FALSE,include.colnames=F,
#       add.to.row = addtorow,
#       file = file.path(WDTables,paste("NGSCountrySEX.tex",sep="")))
#write.csv(NGSCountrySEX ,file = file.path(WDTables,paste("NGSCountrySEX.csv",sep="")))

sum(as.numeric(NGSCountrySEX["Total",]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))

## ------       8.2.2 PER OBSERVATION PROCESS------
NGSCountrySEXoBS <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBS) <- c("",rep(c("Norway","Sweden","Total"),each=2))
#colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



NGSCountrySEXoBS[1,] <- c("",rep(c("F","M"),nYears))
NGSCountrySEXoBS[,1] <- c("",rep(c("Structured","Unstructured"),3))

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(2,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    ## structured
    tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[1], ye[t] + sex1[s] ] <- nrow(tempStruc[tempStruc$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S" | tempStruc$Country %in% "N" , ])
    
    ## Other
    tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[2], ye[t] + sex1[s] ] <- nrow(tempOther[tempOther$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S" | tempOther$Country %in% "N" , ])
  }#t
}

addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBS)[2:ncol(NGSCountrySEXoBS)])),
                                                              '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBS) <- rep("", ncol(NGSCountrySEXoBS))

#"\\multirow{1}{*}{}",
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
NGSCountrySEXoBS <- data.frame(cbind(multirowadd,NGSCountrySEXoBS))
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



# print(xtable(NGSCountrySEXoBS, type = "latex",
#              align = paste(c("l",rep("c",ncol(NGSCountrySEXoBS))),collapse = "")),
#       #scalebox = .7, 
#       floating = FALSE,
#       add.to.row = addtorow,
#       include.colnames = F,
#       include.rownames = FALSE,
#       sanitize.text.function = function(x){x},
#       file = file.path(WDTables,paste("NGSCountrySEXperObs.tex",sep="")))


sum(as.numeric(NGSCountrySEX["Total",]))

sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))

#PLOT CHECK 
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="red",add=T)
plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="red",add=T)

par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="#E69F00",add=T)
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="#009E73",add=T)

#dev.off()
# 
NGSStructured$Year1 <- NGSStructured$Year+1
NGSOther$Year1 <- NGSOther$Year+1

barplot(rbind(table(NGSStructured$Year1),table(NGSOther$Year1)),
        col=c("#E69F00","#009E73"))
legend("topleft",fill=c("#E69F00","#009E73"),legend=c("Structured","Other") )

#GIVE FILE TO HENRIK
tmp <- NGSOther[NGSOther$Year %in% c(2019,2020,2021),]
tmp

#write.csv(tmp,file= file.path(WDTables,paste("Unstructured2020_2022.csv",sep="")))


## ------     8.4. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------

## ------       8.4.1 ALL ------

NGSidCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
#colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
colnames(NGSidCountrySEX) <- c(unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#


NGSidCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
    NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
    NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
    
  }#t
  
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))

write.csv(NGSidCountrySEX ,file = file.path(WDTables,paste("NGSidCountrySEX.csv",sep="")))

# print(xtable(NGSidCountrySEX, type = "latex",
#              align = paste(c("l",rep("c",ncol(NGSidCountrySEX))),collapse = "")),
#       #scalebox = .8, 
#       floating = FALSE,include.colnames=F,
#       add.to.row = addtorow,
#       file = file.path(WDTables,paste("NGSidCountrySEX.tex",sep="")))


### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
NGSidCountryTotal <- matrix(0, ncol = nYears, nrow = 1)
row.names(NGSidCountryTotal) <- c("Total")
#colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))
colnames(NGSidCountryTotal) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

for(t in 1:nYears){
  temp <- NGS[NGS$Year == years[t] , ]
  NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
}#t

write.csv(NGSidCountryTotal,file = file.path(WDTables,paste("TotalIdDetected.csv",sep="")))



## ------       8.4.2 PER OBSERVATION PROCESS------
NGSCountrySEXoBSid <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBSid) <- c("",rep(c("Norway","Sweden","Total"),each=2))
#colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#


NGSCountrySEXoBSid[1,] <- c("",rep(c("F","M"),nYears))
NGSCountrySEXoBSid[,1] <- c("",rep(c("Structured","Unstructured"),3))

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(2,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    ## structured
    tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
    
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[1], ye[t] + sex1[s] ] <- length(unique(tempStruc$Id[tempStruc$Country %in% "N" ])) 
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id[tempStruc$Country %in% "S" ]))
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id))
    
    ## Other
    tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
    
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[2], ye[t] + sex1[s] ] <- length(unique(tempOther$Id[tempOther$Country %in% "N" ])) 
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id[tempOther$Country %in% "S" ]))
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id))
    
    ###TOTAL 
    
    
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBSid)[2:ncol(NGSCountrySEXoBSid)])),
                                                              '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBSid) <- rep("", ncol(NGSCountrySEXoBSid))

#"\\multirow{1}{*}{}",
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
NGSCountrySEXoBSid <- data.frame(cbind(multirowadd,NGSCountrySEXoBSid))
colnames(NGSCountrySEXoBSid) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



# print(xtable(NGSCountrySEXoBSid, type = "latex",
#              align = paste(c("l",rep("c",ncol(NGSCountrySEXoBSid))),collapse = "")),
#       #scalebox = .7, 
#       floating = FALSE,
#       add.to.row = addtorow,
#       include.colnames = F,
#       include.rownames = FALSE,
#       sanitize.text.function = function(x){x},
#       file = file.path(WDTables,paste("NGSCountrySEXperObsid.tex",sep="")))



## ------     8.5. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------
DeadidCountrySEX <- matrix(0, ncol = nYears*2+1, nrow = 6)
row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
#colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))#c("",unlist(lapply(YEARS,function(x) c(x[2],x[2]))))#
colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

DeadidCountrySEX[1,] <- c("",rep(c("F","M"),nYears))
DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
###
MortalityNames <- unique(as.character(dead$DeathCause))
table(as.character(dead$DeathCause))
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
## SEPARATE MORTALITIES
cause <- c("other","legal culling")

for(t in 1:nYears){
  for(s in 1:2){
    for(d in 1:2){
      if(d==1){temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$DeathCause %in% legalCauses), ]
      }else{
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$DeathCause %in% legalCauses, ]
      }
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Norway" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
      
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Sweden" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
    }#t
    DeadidCountrySEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(DeadidCountrySEX[2:6,ye[t] + sex1[s]+1]))
  }
}


##summary
#Other causes
sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
#legal
sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))

sum(as.numeric(DeadidCountrySEX[c(2,3),2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))


## %of dead reco (legal) in norway
sum(as.numeric(DeadidCountrySEX[4,2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(4,5),2:ncol(DeadidCountrySEX)]))



sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="M")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,] %in% c("F","M"))]))

#write latex
addtorow <- list()

addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))
# REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC


multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))
# addtorow$pos <- list(c(0),0)
# addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
#                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))

# xTableState <- xtable(TableState)

# rownames(TableState)[2:5] <- c("$\\rho$","$\\phi$","h","w")

# rownames(TableState)[2:8] <- c("$\\gamma$","$\\phi$","   ","h","  ","w","  ")


# print(xtable(DeadidCountrySEX, type = "latex",
#              align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
#       #scalebox = .7, 
#       floating = FALSE,
#       add.to.row = addtorow,
#       include.colnames = F,
#       include.rownames = FALSE,
#       sanitize.text.function = function(x){x},
#       file = file.path(WDTables,paste("DeadidCountrySEX.tex", sep="")))


#check 
# tmp <- dead[dead$Year == 2019 & 
#               !dead$DeathCause %in% legalCauses &
#               dead$Country %in% "N"& 
#               dead$Sex %in% "Hunn", ]
# length(unique(tmp$Id))
# unique(tmp$Id)
# 
# DeadidCountrySEX[,"X2022"]
# DeadidCountrySEX[,"X2022.1"]
# 
# dead[dead$Id %in% "JI416817 Ind7303 +",]$Sex
# 
# plot(COUNTRIES$geometry)
# plot(tmp$geometry,add=T,col="red",pch=16)

## ------     8.6. GET THE DETECTED INDIVIDUALS ------
n.detected <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv",sep="")))
n.detected <- n.detected[1,2:ncol(n.detected)]

## ------     8.7. SUMMARY DETECTED INDIVIDUALS PER COUNTIES ------
# myFilteredData.sp$alive$COUNTIES  <- st_intersects(myFilteredData.sp$alive[,1], COUNTIES_AGGREGATED[,1])
# myFilteredData.sp$alive$COUNTIES <- as.numeric(myFilteredData.sp$alive$COUNTIES)
# 
# 
# myFilteredData.sp$alive$COUNTIES <- apply(st_intersects(COUNTIES_AGGREGATED, myFilteredData.sp$alive, sparse = FALSE), 2, 
#       function(col) {which(col)})
# myFilteredData.sp$alive$counties1 <- 0
# for(i in 1:nrow(myFilteredData.sp$alive)){
#   if(length(myFilteredData.sp$alive$COUNTIES[[i]])>0){
#   myFilteredData.sp$alive$counties1[i] <- myFilteredData.sp$alive$COUNTIES[[i]][1]
#   }else{
#     myFilteredData.sp$alive$counties1[[i]] <- 0
#   }
# }
# 
# par(mar=c(0,0,0,0))
# plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
# text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
#      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
#      COUNTIES_AGGREGATED$id,col="red",font=2)
# 
# 
# 
# 
# 
# 
# ## NSAMPLES
# summa <- myFilteredData.sp$alive %>%
#   group_by(counties1,Year) %>%
#   summarise(n=n()) %>%
#   st_drop_geometry()
# summa <- summa[summa$counties1>0,]
# #2023
# tmp <- summa[summa$Year %in% 2023,]
# COUNTIES_AGGREGATED$nSampl2023 <- tmp$n#[1:8,"n"]
# #2022
# tmp1 <- summa[summa$Year %in% 2022,]
# COUNTIES_AGGREGATED$nSampl2022 <- tmp1$n#[1:8,"n"]
# 
# ##
# tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("id","nSampl2022","nSampl2023")])
# bar <- t(as.matrix(tmppp[c(5,7,8),2:3]))
# colnames(bar) <- c(5,7,8)
# barplo <- barplot(bar,beside=T,ylab="N samples")
# legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))
# 
# 
# ## NdetectionsPerID 
# #COUNT NUMBER IDS 
# 
# summa$n1 <- summa$NID <- 0
# dpt <- unique(unlist(summa$counties1))
# yearsss <- c(2022,2023)
# for(t in 1:length(yearsss)){
#   for(i in 1:length(dpt)){
#     tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$counties1 %in% dpt[i] & myFilteredData.sp$alive$Year %in% yearsss[t],]
#     
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$n1 <- nrow(tmp)
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$NID <-  length(unique(tmp$Id))
#   }
# }
# #
# summa$NdetPerIDDet <- summa$n1/summa$NID
# 
# #2023
# tmp <- summa[summa$Year %in% 2023,]
# COUNTIES_AGGREGATED$detPerID2023 <- tmp$NdetPerIDDet#[1:8,"n"]
# #2022
# tmp1 <- summa[summa$Year %in% 2022,]
# COUNTIES_AGGREGATED$detPerID2022 <- tmp1$NdetPerIDDet#[1:8,"n"]
# 
# 
# tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("detPerID2022","detPerID2023")])
# bar <- t(as.matrix(tmppp[c(5,7,8),]))
# colnames(bar) <- c(5,7,8)
# 
# 
# par(mfrow=c(1,2))
# par(mar=c(0,0,0,0))
# plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
# text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
#      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
#      COUNTIES_AGGREGATED$id,col="red",font=2)
# par(mar=c(4,5,1,1))
# barplo <- barplot(bar,beside=T,ylab="average Dets per IDS")
# legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))





## ------   4.GET THE MCMC ESTIMATES  -----

## ------     4.1.GET AND COMPILE BITES ------

# COMPILE CHARACTERISTICS 
bitesize <- 500
burnin <- 20000
NSkipBites <- burnin/bitesize
nthinsxy <- 5 #thinnumber for the sxy AND z values (necessary to save memory) 

## ------       4.1.1 FEMALES ------

#outDirectories <- list.files(file.path(WD, modelNameF))[grep(paste("FORSnap", modelNameF,sep="") ,list.files(file.path(WD, modelNameF)))]
outDirectories <- list.files(file.path(WD, modelName, "Hunn"))[grep(paste("NimbleOutFORSnap2015_", modelName,sep=""),
                                                                                  list.files(file.path(WD, modelName,"Hunn")))]

path.list <- file.path(WD, modelName, "Hunn",outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  #length(files)/2  ## if we have another ifle for sxy and z
  length(files)
}))
minBites <- floor(min(numBites))
#minBites <- 95

## GO TROUGH THE BITES AND GET THEM
nimOutput <- nimOutputSXY <- RUNTIME <- list()
nimOutputList <- nimOutputSXYList <- RUNTIME <- list()
myResultsSXYZ_FList <- myResults_FList <- nidMinFList <- list()
#gc()
for(t in 1:nYears){
  path.list <- 0
  for(ch in 1:4){
    path.list[ch] <- file.path(WD, modelName, "Hunn", paste("NimbleOutFORSnap", years[t],"_",
                                                                          modelName,"Hunn","_Chain", ch, ".RData", sep = "")) 
  }
  
  nidMin <-0
  for(ch in 1:4){
    # load(file.path(WD, modelName, "Hunn",
    #                paste("Snap", modelName,years[t],"_", ch, ".RData", sep = "")))
    load(file.path(WD, modelName, "Hunn/Snapshot",
                   paste("Snap",years[t],"_", modelName,"Hunn_Chain", ch, ".RData", sep = "")))
    
    nidMin[ch] <- nimConstants$n.individuals
  }
  nidMin <- min(nidMin)
  nidMinFList[[t]] <- min(nidMin)
  
  
  for(p in 1:length(path.list)){
    print(path.list[p])
    outfiles <- list.files(path.list[p])
    out <- outSXYZ <- runtime <- list()#[CM]

    for(x in NSkipBites:minBites){
      print(x)
      load(file.path(path.list[p], paste("MCMC_bite_", x, ".RData", sep = "")))
      #runtime[[x]] <- RunTime[3]
      #params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
      # REMOVE SXY AND Z, NO THINING
      #parmIndex <- which(! params.simple %in% c("sxy","z"))
      
      if(sum(is.na(mcmcSamples))>0){
        id <- unique(which(is.na(this.sample),arr.ind = T)[,1])
        mcmcSamples <- mcmcSamples[-id,]
        print(x)
        print("here")
        
      }
      
      out[[x]] <- mcmcSamples#[,parmIndex]#[ ,parmIndex]
      
      # KEEP SXY AND Z, THINING
      outSXYZ[[x]] <- mcmcSamples2#[,subset]#[nthins,parmIndex]#[ ,parmIndex] #this.sampleSxyZ
      
    }#x
    RUNTIME[[p]] <- unlist(runtime)#[CM]
    out.mx <- do.call(rbind, out)
    out.mxSXY <- do.call(rbind, outSXYZ)
    
    nimOutput[[p]] <- as.mcmc(out.mx)
    nimOutputSXY[[p]] <- as.mcmc(out.mxSXY)
  }#p
  
  ## COMPILE THE RESULTS
  nimOutputList[[t]] <- as.mcmc.list(nimOutput)
  nimOutputSXYList[[t]] <- as.mcmc.list(nimOutputSXY)
  
  myResults_FList[[t]] <- ProcessCodaOutput(nimOutputList[[t]], params.omit = c("sxy","z"))
  myResultsSXYZ_FList[[t]] <- ProcessCodaOutput(nimOutputSXYList[[t]], params.omit = c("sxy","z"))
  
}


## ------       1.2 MALES ------
#outDirectories <- list.files(file.path(WD, modelNameF))[grep(paste("FORSnap", modelNameF,sep="") ,list.files(file.path(WD, modelNameF)))]
outDirectories <- list.files(file.path(WD, modelName, "Hann"))[grep(paste("NimbleOutFORSnap2015_", modelName,sep=""),
                                                                                  list.files(file.path(WD, modelName,"Hann")))]

path.list <- file.path(WD, modelName, "Hann",outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  #length(files)/2  ## if we have another ifle for sxy and z
  length(files)
}))
minBites <- floor(min(numBites))
#minBites <- 95

## GO TROUGH THE BITES AND GET THEM
nimOutput <- nimOutputSXY <- RUNTIME <- list()
nimOutputList <- nimOutputSXYList <- RUNTIME <- list()
myResultsSXYZ_MList <- myResults_MList <- nidMinMList <- list()
#gc()
for(t in 1:nYears){
  path.list <- 0
  for(ch in 1:4){
    path.list[ch] <- file.path(WD, modelName, "Hann", paste("NimbleOutFORSnap", years[t],"_",
                                                                          modelName,"Hann","_Chain", ch, ".RData", sep = "")) 
  }
  
  nidMin <-0
  for(ch in 1:4){
    # load(file.path(WD, modelName, "Hann",
    #                paste("Snap", modelName,years[t],"_", ch, ".RData", sep = "")))
    load(file.path(WD, modelName, "Hann/Snapshot",
                   paste("Snap",years[t],"_", modelName,"Hann_Chain", ch, ".RData", sep = "")))
    
    nidMin[ch] <- nimConstants$n.individuals
  }
  nidMin <- min(nidMin)
  nidMinMList[[t]] <- min(nidMin)
  
  
  for(p in 1:length(path.list)){
    print(path.list[p])
    outfiles <- list.files(path.list[p])
    out <- outSXYZ <- runtime <- list()#[CM]

    for(x in NSkipBites:minBites){
      print(x)
      load(file.path(path.list[p], paste("MCMC_bite_", x, ".RData", sep = "")))
      #runtime[[x]] <- RunTime[3]
      #params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
      # REMOVE SXY AND Z, NO THINING
      #parmIndex <- which(! params.simple %in% c("sxy","z"))
      
      if(sum(is.na(mcmcSamples))>0){
        id <- unique(which(is.na(this.sample),arr.ind = T)[,1])
        mcmcSamples <- mcmcSamples[-id,]
        print(x)
        print("here")
        
      }
      
      out[[x]] <- mcmcSamples#[,parmIndex]#[ ,parmIndex]
      
      # KEEP SXY AND Z, THINING
      outSXYZ[[x]] <- mcmcSamples2#[,subset]#[nthins,parmIndex]#[ ,parmIndex] #this.sampleSxyZ
      
    }#x
    RUNTIME[[p]] <- unlist(runtime)#[CM]
    out.mx <- do.call(rbind, out)
    out.mxSXY <- do.call(rbind, outSXYZ)
    
    nimOutput[[p]] <- as.mcmc(out.mx)
    nimOutputSXY[[p]] <- as.mcmc(out.mxSXY)
  }#p
  
  ## COMPILE THE RESULTS
  nimOutputList[[t]] <- as.mcmc.list(nimOutput)
  nimOutputSXYList[[t]] <- as.mcmc.list(nimOutputSXY)
  
  
  myResults_MList[[t]] <- ProcessCodaOutput(nimOutputList[[t]], params.omit = c("sxy","z"))
  myResultsSXYZ_MList[[t]] <- ProcessCodaOutput(nimOutputSXYList[[t]], params.omit = c("sxy","z"))
  
}

### combine all years together 
maxM <- max(unlist(nidMinMList))
maxF <- max(unlist(nidMinFList))
myResultsSXYZ_M <- myResultsSXYZ_MList[[1]] 
myResultsSXYZ_F <- myResultsSXYZ_FList[[1]] 

myResultsSXYZ_M$sims.list$sxy <- array(0,c(dim(myResultsSXYZ_MList[[1]]$sims.list$sxy)[1],maxM ,2,nYears) )
myResultsSXYZ_M$sims.list$z <- array(0,c(dim(myResultsSXYZ_MList[[1]]$sims.list$z)[1],maxM ,nYears) )
myResultsSXYZ_M$sims.list$N <- array(0,c(length(myResults_MList[[1]]$sims.list$N)[1] ,nYears) )
myResultsSXYZ_M$sims.list$sigma <- array(0,c(length(myResults_MList[[1]]$sims.list$sigma)[1] ,nYears) )
myResultsSXYZ_M$sims.list$p0 <- array(0,c(dim(myResults_MList[[1]]$sims.list$p0) ,nYears) )
myResultsSXYZ_M$sims.list$p0Oth <- array(0,c(dim(myResults_MList[[1]]$sims.list$p0Oth) ,nYears) )
myResultsSXYZ_M$sims.list$betaCovs <- array(0,c(dim(myResults_MList[[1]]$sims.list$betaCovs) ,nYears) )
myResultsSXYZ_M$sims.list$betaCovsOth <- array(0,c(dim(myResults_MList[[1]]$sims.list$betaCovsOth) ,nYears) )
myResultsSXYZ_M$sims.list$betaDens <- array(0,c(length(myResults_MList[[1]]$sims.list$betaDens) ,nYears) )
myResultsSXYZ_M$sims.list$betaResponse <- array(0,c(length(myResults_MList[[1]]$sims.list$betaResponse) ,nYears) )
myResultsSXYZ_M$sims.list$betaResponseOth <- array(0,c(length(myResults_MList[[1]]$sims.list$betaResponseOth) ,nYears) )

myResultsSXYZ_F$sims.list$sxy <- array(0,c(dim(myResultsSXYZ_FList[[1]]$sims.list$sxy)[1],maxF ,2,nYears) )
myResultsSXYZ_F$sims.list$z <- array(0,c(dim(myResultsSXYZ_FList[[1]]$sims.list$z)[1],maxF ,nYears) )
myResultsSXYZ_F$sims.list$N <- array(0,c(length(myResults_FList[[1]]$sims.list$N)[1] ,nYears) )
myResultsSXYZ_F$sims.list$sigma <- array(0,c(length(myResults_FList[[1]]$sims.list$sigma)[1] ,nYears) )
myResultsSXYZ_F$sims.list$p0 <- array(0,c(dim(myResults_FList[[1]]$sims.list$p0) ,nYears) )
myResultsSXYZ_F$sims.list$p0Oth <- array(0,c(dim(myResults_FList[[1]]$sims.list$p0Oth) ,nYears) )
myResultsSXYZ_F$sims.list$betaCovs <- array(0,c(dim(myResults_FList[[1]]$sims.list$betaCovs) ,nYears) )
myResultsSXYZ_F$sims.list$betaCovsOth <- array(0,c(dim(myResults_FList[[1]]$sims.list$betaCovsOth) ,nYears) )
myResultsSXYZ_F$sims.list$betaDens <- array(0,c(length(myResults_FList[[1]]$sims.list$betaDens) ,nYears) )
myResultsSXYZ_F$sims.list$betaResponse <- array(0,c(length(myResults_FList[[1]]$sims.list$betaResponse) ,nYears) )
myResultsSXYZ_F$sims.list$betaResponseOth <- array(0,c(length(myResults_FList[[1]]$sims.list$betaResponseOth) ,nYears) )

for(t in 1:nYears){
  myResultsSXYZ_M$sims.list$sxy[,1:nidMinMList[[t]],,t] <-   myResultsSXYZ_MList[[t]]$sims.list$sxy
  myResultsSXYZ_M$sims.list$z[,1:nidMinMList[[t]],t] <-   myResultsSXYZ_MList[[t]]$sims.list$z
  myResultsSXYZ_M$sims.list$N[,t] <-   myResults_MList[[t]]$sims.list$N
  myResultsSXYZ_M$sims.list$sigma[,t] <-   myResults_MList[[t]]$sims.list$sigma
  myResultsSXYZ_M$sims.list$p0[,,t] <-   myResults_MList[[t]]$sims.list$p0
  myResultsSXYZ_M$sims.list$p0Oth[,,t] <-   myResults_MList[[t]]$sims.list$p0Oth
  myResultsSXYZ_M$sims.list$betaCovs[,,t] <-   myResults_MList[[t]]$sims.list$betaCovs
  myResultsSXYZ_M$sims.list$betaCovsOth[,,t] <-   myResults_MList[[t]]$sims.list$betaCovsOth
  myResultsSXYZ_M$sims.list$betaDens[,t] <-   myResults_MList[[t]]$sims.list$betaDens
  myResultsSXYZ_M$sims.list$betaResponse[,t] <-   myResults_MList[[t]]$sims.list$betaResponse
  myResultsSXYZ_M$sims.list$betaResponseOth[,t] <-   myResults_MList[[t]]$sims.list$betaResponseOth

  myResultsSXYZ_F$sims.list$sxy[,1:nidMinFList[[t]],,t] <-   myResultsSXYZ_FList[[t]]$sims.list$sxy
  myResultsSXYZ_F$sims.list$z[,1:nidMinFList[[t]],t] <-   myResultsSXYZ_FList[[t]]$sims.list$z
  myResultsSXYZ_F$sims.list$N[,t] <-   myResults_FList[[t]]$sims.list$N
  myResultsSXYZ_F$sims.list$sigma[,t] <-   myResults_FList[[t]]$sims.list$sigma
  myResultsSXYZ_F$sims.list$p0[,,t] <-   myResults_FList[[t]]$sims.list$p0
  myResultsSXYZ_F$sims.list$p0Oth[,,t] <-   myResults_FList[[t]]$sims.list$p0Oth
  myResultsSXYZ_F$sims.list$betaCovs[,,t] <-   myResults_FList[[t]]$sims.list$betaCovs
  myResultsSXYZ_F$sims.list$betaCovsOth[,,t] <-   myResults_FList[[t]]$sims.list$betaCovsOth
  myResultsSXYZ_F$sims.list$betaDens[,t] <-   myResults_FList[[t]]$sims.list$betaDens
  myResultsSXYZ_F$sims.list$betaResponse[,t] <- myResults_FList[[t]]$sims.list$betaResponse
  myResultsSXYZ_F$sims.list$betaResponseOth[,t] <- myResults_FList[[t]]$sims.list$betaResponseOth
}
myResultsList <-  list(myResultsSXYZ_F,myResultsSXYZ_M)



## ------     1.3 RESCALE SXY SO THEY CAN BE COMPARED ------

dimnames(myResultsSXYZ_F$sims.list$sxy)[[3]] <- c("x", "y")
myResultsSXYZ_F$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ_F$sims.list$sxy,
                                                          coordsHabitatGridCenter = myHabitat.list$habitat.xy,
                                                          scaleToGrid = FALSE)$coordsDataScaled

dimnames(myResultsSXYZ_M$sims.list$sxy)[[3]] <- c("x", "y")
myResultsSXYZ_M$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ_M$sims.list$sxy,
                                                          coordsHabitatGridCenter = myHabitat.list$habitat.xy,
                                                          scaleToGrid = FALSE)$coordsDataScaled


##QUICK CHECK
t=9
#Male
plot(st_geometry(myHabitat.list$buffered.habitat.poly))
points(myResultsSXYZ_M$sims.list$sxy[1,myResultsSXYZ_M$sims.list$z[1,,t] %in% 1, 2, t]~
         myResultsSXYZ_M$sims.list$sxy[1,myResultsSXYZ_M$sims.list$z[1,,t] %in% 1, 1, t], pch=16, col="red")

#Female
plot(st_geometry(myHabitat.list$buffered.habitat.poly))
points(myResultsSXYZ_F$sims.list$sxy[1,myResultsSXYZ_F$sims.list$z[1,,t] %in% 1, 2, t]~
         myResultsSXYZ_F$sims.list$sxy[1,myResultsSXYZ_F$sims.list$z[1,,t] %in% 1, 1, t], pch=16, col="red")



## ------     1.4 IDENTIFY INDIVIDUALS IN THE BUFFER AND ASSIGN THEM A Z NOT ALIVE (Z==5) ------
## ------       1.4.1 CREATE POLYGON WITHOUT BUFFER ------
##MALES
habbRNobuffM <- myHabitat.list$habitat.rWthBuffer
##FEMALES
habbRNobuffF <- myHabitat.list$habitat.rWthBuffer

## ------       1.4.2 IDENTIFY INDIVIDUALS IN THE BUFFER AND GIVE THEM A STATE 5 ------
#MAKE A COPY 
myResultsSXYZ_F$sims.list$z1 <- myResultsSXYZ_F$sims.list$z 
myResultsSXYZ_M$sims.list$z1  <- myResultsSXYZ_M$sims.list$z 


##FEMALES
dim(myResultsSXYZ_F$sims.list$sxy)
dim( myResultsSXYZ_F$sims.list$z)
for(t in 1:nYears){
  for(i in 1:dim( myResultsSXYZ_F$sims.list$z)[1]){
    whichNA <- which(is.na(habbRNobuffF[cellFromXY(habbRNobuffF, myResultsSXYZ_F$sims.list$sxy[i,,1:2,t])]))
    myResultsSXYZ_F$sims.list$z[i,whichNA,t] <- 5
  }
  print(t)
}


###IDENTIFY INDIVIDUALS IN THE BUFFER AND GIVE THEM A STATE 5
##MALES
gc()
dim(myResultsSXYZ_M$sims.list$sxy)
dim( myResultsSXYZ_M$sims.list$z)
for(t in 1:nYears){
  for(i in 1:dim( myResultsSXYZ_M$sims.list$z)[1]){
    whichNA <- which(is.na(habbRNobuffM[cellFromXY(habbRNobuffM,myResultsSXYZ_M$sims.list$sxy[i,,1:2,t])]))
    myResultsSXYZ_M$sims.list$z[i,whichNA,t] <- 5
  }
  # gc()
  print(t)
}

gc()


## ------     1.4 COMBINE MALES AND FEMALES ------
## ASSUMING THE SAME NUMBER OF CHAINS AND ITERATIONS/ WE CAN JUST ADD N ESIMAES OF FEMALES AND MALES
myResultsSXYZ_MF <- myResultsSXYZ_M
## HERE I ONLY HAVE ONE CHAIN FOR THE MALES
dimmF <- dim(myResultsSXYZ_F$sims.list$sxy )[1]
dimmFsxy <- dim(myResultsSXYZ_F$sims.list$sxy)[1]


## ------     1.5 COMBINE SXY AND Z ------
myResultsSXYZ_MF$sims.list$sxy <- abind(myResultsSXYZ_M$sims.list$sxy[1:dimmFsxy[1],,,] , myResultsSXYZ_F$sims.list$sxy, along = 2 )
dimnames(myResultsSXYZ_MF$sims.list$sxy)[[3]] <- c("x", "y")

myResultsSXYZ_MF$sims.list$z <- abind(myResultsSXYZ_M$sims.list$z[1:dimmFsxy[1],,],
                                      myResultsSXYZ_F$sims.list$z, along = 2 )


myResultsSXYZ_MF$sims.list$sigma <- abind(myResultsSXYZ_M$sims.list$sigma[1:dimmF[1],] * res(myHabitat.list$habitat.r)[1],
                                          myResultsSXYZ_F$sims.list$sigma* res(myHabitat.list$habitat.r)[1], along = 1 )


#myResultsSXYZ_MF$sims.list$sigma <- abind(myResultsSXYZ_M$sims.list$sigma[1:dimmF[1]], myResultsSXYZ_F$sims.list$sigma, along = 2 )
myResultsSXYZ_MF$sims.list$sex <- rep(c("M","F"), c(dim(myResultsSXYZ_M$sims.list$sxy)[2], dim(myResultsSXYZ_F$sims.list$sxy)[2]))


#EMPTY USELESS ARRAYS
myResultsSXYZ_F$sims.list$z <-NULL
myResultsSXYZ_M$sims.list$z <-NULL

myResultsSXYZ_M$sims.list$N <-NULL
myResultsSXYZ_F$sims.list$N <-NULL

myResultsSXYZ_F$sims.list$sxy <-NULL
myResultsSXYZ_M$sims.list$sxy <-NULL

gc()

# # #select 100 iterations for RB
# dim(myResultsSXYZ_MF$sims.list$sxy)
# nit <- sample(dim(myResultsSXYZ_MF$sims.list$sxy)[1],100)
# sxy <- round(myResultsSXYZ_MF$sims.list$sxy[nit,,,],digits=5)
# z <- myResultsSXYZ_MF$sims.list$z[nit,,]
# 
# save(sxy, z,
#      file = file.path(WDFigures, "Itera.RData"))
# #contains 100 iterations of "sxy" and "z"
# load("C://Users//cymi//Dropbox (Old)//AQEG Dropbox//AQEG Team Folder//RovQuant//wolverine//CM//2024/plot53Cleaned2024/Figure/Itera.RData")
# 

#MERGE RESULTS IN A LIST
#Results.list <- list(myResults_F,myResults_M)
Results.list <- list(myResultsSXYZ_F,myResultsSXYZ_M)
names(Results.list) <- c("F","M")

## ------     1.6 SAVE AND LOAD DATA ------
# save(Results.list, myResultsSXYZ_MF,
#     file = file.path(paste(dir.dropbox,"/wolverine/CM/2022/plot25Cleaned/Figure/",sep=""), "MCMC.RData" ))
# load(file.path(paste(dir.dropbox,"/wolverine/CM/2022/plot25Cleaned/Figure/",sep=""), "MCMC.RData" ))

Results.list[["F"]]$mean$sigma
Results.list[["M"]]$mean$sigma

myResults_F <- Results.list[["F"]]
myResults_M <- Results.list[["M"]]

## ------     1.7 CHECK RHAT ------
myResults_F$Rhat
myResults_M$Rhat
# basicMCMCplots::chainsPlot(nimOutput,var = "omeg1")
# basicMCMCplots::chainsPlot(nimOutput,var = "N[1]")
# basicMCMCplots::chainsPlot(nimOutput,var = "N")
# basicMCMCplots::chainsPlot(nimOutput,var = "betaDens")

# WHAT STATES ARE CONSIDERED AS ALIVE IN THE MODEL
#alive.states <- c(2) 
alive.states <- c(1) 

#years not sampled in Norrbotten
yearsSampledNorrb <- c(2016:2018,2023,2024)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)


## ------   2. AC BASED DENSITY  (5km) ------
### REMOVE THE BUFFER FROM THE HABITAT ###
### COUNTRIES 
rrCountries <- habitatRasterResolution$`5km`[["Countries"]]
# REMOVE FINLAND AND RUSSIA 
rrCountries[rrCountries%in% c(1,3)] <- NA
plot(rrCountries)

habitat.rWthBuffer <- myHabitat.list$habitat.rWthBuffer
habitat.rWthBuffer[habitat.rWthBuffer %in% 0] <- NA

searchedPolygon <- sf::st_as_sf(stars::st_as_stars(habitat.rWthBuffer), 
                                as_points = FALSE, merge = TRUE)
searchedPolygon <- searchedPolygon[searchedPolygon$Habitat>0,]
# searchedPolygon <- rasterToPolygons(habitat.rWthBuffer, dissolve = T, function(x) x==1 )

rrCountries <- mask(rrCountries, searchedPolygon)
rrCountries <- crop(rrCountries, myHabitat.list$habitat.r)
plot(rrCountries)
### REGIONS AND COUNTIES 
rrRegions <- habitatRasterResolution$`5km`[["Regions"]] 
## deal with the special characters
levels(rrRegions)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <-  c("Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
                                                                    "Värmland","Örebro","Västmanland","Gävleborg",
                                                                    "Västernorrland","Jämtland" ,"Västerbotten")

# REMOVE FINLAND AND RUSSIA 
rrRegions[habitatRasterResolution$`5km`[["Countries"]]%in% c(1, 3)] <- NA
plot(rrRegions)
rrRegions <- mask(rrRegions, searchedPolygon)
rrRegions <- crop(rrRegions, myHabitat.list$habitat.r)

# habitatPolygon <- rasterToPolygons(myHabitat.list$habitat.r, dissolve = T, function(x) x==1 )
habitatPolygon <- sf::st_as_sf(stars::st_as_stars(myHabitat.list$habitat.r), 
                               as_points = FALSE, merge = TRUE)
habitatPolygon <- habitatPolygon[habitatPolygon$Habitat>0,]

habitatPolygon5km <- mask(habitatRasterResolution$`5km`[["Habitat"]], habitatPolygon)
habitatPolygon5km <- crop(habitatRasterResolution$`5km`[["Habitat"]], myHabitat.list$habitat.r)

plot(rrRegions)
plot(habitatPolygon5km)
# points(myDetectors$main.detector.sp,pch=16,cex=0.5)


### SWEDISH REGIONS 
rrRegionsSwe <- habitatRasterResolution$`5km`[["Regions"]]
# REMOVE FINLAND AND RUSSIA 
rrRegionsSwe[habitatRasterResolution$`5km`[["Countries"]]%in% c(1,2, 3)] <- NA
plot(rrRegionsSwe)
rrRegionsSwe <- mask(rrRegionsSwe, searchedPolygon)
rrRegionsSwe <- crop(rrRegionsSwe, myHabitat.list$habitat.r)

plot(rrRegionsSwe)
rrRegionsSwe[]
rrRegionsSwe[rrRegionsSwe[]%in% c(18,19,20,21)] <- 1
rrRegionsSwe[rrRegionsSwe[]%in% c(13,17,16,14,15,12,22,3)] <- 2
rrRegionsSwe[rrRegionsSwe[]%in% c(4,5,10,6,7,9,11,8)] <- 3

rrRegionsSwe <- ratify(rrRegionsSwe)
df <- data.frame("ID"=c(1,2,3), "Regions"= c("Nordre","Midtre","SÃ¸ndre"))
levels(rrRegionsSwe)[[1]] <- df

##NorwegianCounty 
rrCountiesNor <- habitatRasterResolution$`5km`[["Counties"]]
rrCountiesNor[rrCountiesNor[] %in% c(1,2, 3)] <- NA
rrCountiesNor <- mask(rrCountiesNor, searchedPolygon)
rrCountiesNor <- crop(rrCountiesNor, myHabitat.list$habitat.r)
#remove sweden
rrCountiesNor[rrCountries[]%in%4] <- NA
plot(rrCountiesNor)
plot(rrRegions)
plot(rrCountries)


#library(lattice)
#lattice::levelplot(rrRegions, col.regions=rev(terrain.colors(20)), xlab="", ylab="")
##
# Nordre
# Jämtland
# Västernorrland
##
# Midtre
# Värmland
# Gävleborg
# Dalarna
# Örebro
# Västmanland
# Västra Götaland
# Uppsala
# Stockholm
##
# Søndre
# Södermanland
# Östergötland
# Jönköping
# Skåne

gc()
### GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
## COUNTRY
densityInputCountries <- getDensityInput( regions = rrCountries
                                          , 
                                          habitat = habitatPolygon5km
                                          ,
                                          s = myResultsSXYZ_MF$sims.list$sxy
                                          ,
                                          plot.check = TRUE
)
### GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
## REGIONS
densityInputRegions <- getDensityInput( regions = rrRegions
                                        , 
                                        habitat = habitatPolygon5km
                                        ,
                                        s = myResultsSXYZ_MF$sims.list$sxy
                                        ,
                                        plot.check = TRUE
)

#Swedish Regions
densityInputRegionsSwe <- getDensityInput( regions = rrRegionsSwe
                                           , 
                                           habitat = habitatPolygon5km
                                           ,
                                           s = myResultsSXYZ_MF$sims.list$sxy
                                           ,
                                           plot.check = TRUE
)

## MERGE COUNTRY AND REGION MATRICES TO ALLOW SIMULTANEOUS EXTRACTION
regionID <- rbind (densityInputCountries$regions.rgmx,
                   densityInputRegions$regions.rgmx,
                   densityInputRegionsSwe$regions.rgmx)
row.names(regionID) <- c(row.names(densityInputCountries$regions.rgmx),
                         row.names(densityInputRegions$regions.rgmx),
                         row.names(densityInputRegionsSwe$regions.rgmx))


#Norwegian Counties
densityInputRegionsNor <- getDensityInput( regions = rrCountiesNor
                                           , 
                                           habitat = habitatPolygon5km
                                           ,
                                           s = myResultsSXYZ_MF$sims.list$sxy
                                           ,
                                           plot.check = TRUE
)


## ------     2.1.1 GET THE % OF REGIONS COVERED BY THE ANALYSIS ------
## Percentage of each county included in the analysis
## GET THE PERCENTAGE FOR ALL REGIONS 
rrRegions1 <- habitatRasterResolution$`5km`[["Regions"]] 
## deal with the special characters
levels(rrRegions1)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <-  c("Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
                                                                     "Värmland","Örebro","Västmanland","Gävleborg",
                                                                     "Västernorrland","Jämtland" ,"Västerbotten")
## CALCULATE AREA OF EACH COUNTY
AreaStudiedRegion <- table(factorValues(rrRegions, rrRegions[]))*res(rrRegions)[1]*1e-6
TotalArea <- table(factorValues(rrRegions1, rrRegions1[]))*res(rrRegions1)[1]*1e-6
#Percentage of counties included in the analysis
areaRegions <- AreaStudiedRegion/TotalArea[names(AreaStudiedRegion)]

## GET THE PERCENTAGE FOR THE 3 SWEDISH UNITS REGIONS 
rrRegionsSwe1 <- habitatRasterResolution$`5km`[["Regions"]]
# REMOVE FINLAND AND RUSSIA 
rrRegionsSwe1[habitatRasterResolution$`5km`[["Countries"]]%in% c(1,2, 3)] <- NA
plot(rrRegionsSwe1)

rrRegionsSwe1[rrRegionsSwe1[]%in% c(18,19,20,21)] <- 1
rrRegionsSwe1[rrRegionsSwe1[]%in% c(13,17,16,14,15,12,22,3)] <- 2
rrRegionsSwe1[rrRegionsSwe1[]%in% c(4,5,10,6,7,9,11,8)] <- 3

rrRegionsSwe1 <- ratify(rrRegionsSwe1)
df <- data.frame("ID"=c(1,2,3), "Regions"= c("Nordre","Midtre","SÃ¸ndre"))
levels(rrRegionsSwe1)[[1]] <- df
plot(rrRegionsSwe1)
## CALCULATE AREA OF EACH UNIT
AreaStudiedSwe <- table(factorValues(rrRegionsSwe, rrRegionsSwe[]))*res(rrRegionsSwe)[1]*1e-6
TotalAreaSwe <- table(factorValues(rrRegionsSwe1, rrRegionsSwe1[]))*res(rrRegionsSwe1)[1]*1e-6#Percentage of counties included in the analysis
areaRegionsSwe <- AreaStudiedSwe/TotalAreaSwe[names(AreaStudiedSwe)]

## GET THE PERCENTAGE FOR THE COUNTRIES 
rrCountries1 <- habitatRasterResolution$`5km`[["Countries"]]
# REMOVE FINLAND AND RUSSIA 
rrCountries1[rrCountries1%in% c(1,3)] <- NA
TotalAreaCountry <- table(factorValues(rrCountries1, rrCountries1[]))*res(rrCountries1)[1]*1e-6
AreaStudiedCountry <- table(factorValues(rrCountries, rrCountries[]))*res(rrCountries)[1]*1e-6
TotalAreaCountry["Total"] <- sum(TotalAreaCountry)
AreaStudiedCountry["Total"] <- sum(AreaStudiedCountry)
## CALCULATE AREA OF EACH COUNTRY
areaCountry <- AreaStudiedCountry/TotalAreaCountry[names(AreaStudiedCountry)]

### MERGE THE PERCENTAGE 
areaAllRegions <- c(areaCountry, areaRegionsSwe, areaRegions)

## ------     2.1 MALE AND FEMALES (5km) ------
ite <- seq(1,dim(densityInputCountries$sx[,,t])[1],by=25)
gc()
## EXTRACT DENSITY 
DensityCountriesRegions <- list()
for(t in 1:nYears){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,,t],
    sy =  densityInputCountries$sy[ite,,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}

DensityCountriesRegions[[t]]$summary

## SAVE
#save(DensityCountriesRegions, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegions.RData" ))

##save Object to calculate growth rate
posteriorRegions <- list()
for(t in 1:nYears){
  posteriorRegions[[t]] <- DensityCountriesRegions[[t]]$PosteriorRegions
}
save(posteriorRegions, 
     file=file.path(WDFigures, paste("posteriorRegions.RData", sep="")))

gc()



## ------     2.2 MALE (5km) ------

IDMales <- which(myResultsSXYZ_MF$sims.list$sex=="M")

DensityCountriesRegionsM <- list()
for(t in 1:nYears){
  DensityCountriesRegionsM[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDMales,t],
    sy =  densityInputCountries$sy[ite,IDMales,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,IDMales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}

DensityCountriesRegionsM[[t]]$summary

## SAVE
#save(DensityCountriesRegionsM, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegionsM.RData" ))

gc()



## ------     2.3 FEMALE (5km) ------

IDFemales <- which(myResultsSXYZ_MF$sims.list$sex=="F")

DensityCountriesRegionsF <- list()
for(t in 1:nYears){
  DensityCountriesRegionsF[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDFemales,t],
    sy =  densityInputCountries$sy[ite,IDFemales,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,IDFemales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}

DensityCountriesRegionsF[[t]]$summary
## SAVE
#save(DensityCountriesRegionsF, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegionsF.RData" ))


## ------     2.4 COUNTIES NORWAY M AND F (5km)  ------

DensityCountriesRegionsNOR <- list()
for(t in 1:nYears){
  DensityCountriesRegionsNOR[[t]] <- GetDensity_PD(
    sx = densityInputRegionsNor$sx[ite,,t],
    sy =  densityInputRegionsNor$sy[ite,,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputRegionsNor$habitat.id,
    aliveStates = alive.states,
    regionID = densityInputRegionsNor$regions.rgmx,
    returnPosteriorCells = F)
}



## ------     2.4 SUMMARY TABLES ------

## ------       2.4.1 ALL YEARS, BOTH SEX ------

idcounty <- row.names(DensityCountriesRegions[[t]]$summary)
#REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
#GET NORWEGIAN VERSUS SWEDISH COUNTIES 
idcountyNOR <- idcounty[grep("Region",idcounty)]
idcountySWE <- sort(idcounty[-grep("Region",idcounty)])
idcountyTable <- c("Total","Norway", idcountyNOR, "Sweden" ,idcountySWE)

CountyNorth <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 1], layer=1)[,1])
CountyMiddle <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 2], layer=1)[,1])
CountySouth <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 3], layer=1)[,1])

idcountyTable <- c("Total","Norway",
                   idcountyNOR,
                   "Sweden" ,
                   "Nordre",
                   idcountySWE[idcountySWE%in%CountyNorth],
                   "Midtre",
                   idcountySWE[idcountySWE%in%CountyMiddle],
                   "SÃ¸ndre",
                   idcountySWE[idcountySWE%in%CountySouth]
)


#CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimates <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimates) <- c(idcountyTable)
colnames(NCarRegionEstimates) <-  unlist(lapply(YEARS,function(x) c(paste(x, collapse = "/"))))#unlist(lapply(YEARS ,function(x) x[2]))#

#FILL IN THE TABLE 
for(t in 1:nYears){
  for( i in 1:length(idcountyTable)){
    NCarRegionEstimates[idcountyTable[i],t] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                     " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                     round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}

##QUICK CHECK TO MAKE SURE VALUES SUMS UP 
tmp <- DensityCountriesRegions[[t]]$summary[1:(nrow(DensityCountriesRegions[[t]]$summary)),]
# SWE
row.names(DensityCountriesRegions[[t]]$summary)
sum(tmp[idcountySWE,"mean"])
tmp["Sweden","mean"]
#NOR
sum(tmp[idcountyNOR,"mean"])
tmp["Norway","mean"]
#TOTAL
sum(tmp[c(idcountyNOR,idcountySWE),"mean"])
tmp["Total","mean"]



## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")
idcountySWE1 <- idcountySWE

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
# 
# row.names(NCarRegionEstimates) <- idcounty1
# NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste(NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*", sep="")

#print csv
write.csv(NCarRegionEstimates,
          file = file.path(WDTables,paste("NAllYears.csv",sep="")),fileEncoding="latin1")





idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten}"
row.names(NCarRegionEstimates) <- idcounty1
NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*","}", sep="")

NCarRegionEstimates["SWEDEN",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["SWEDEN",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["TOTAL",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["TOTAL",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["Nordre",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["Nordre",yearsNotSampled], "**","}", sep="")


row.names(NCarRegionEstimates) <- c("TOTAL",
                                    paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                    paste("\\hspace{0.5cm} ",
                                          idcountyNOR,sep=""),
                                    paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                    paste("\\hspace{0.5cm}","Norra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                    paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                    paste("\\hspace{0.5cm}","Södra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCarRegionEstimates)[grep("VÃ¤straGÃ¶taland", row.names(NCarRegionEstimates))] <- paste("\\hspace{0.75cm}",
                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

# row.names(NCarRegionEstimates) <- c("TOTAL",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1,sep="")
# )


print(xtable(NCarRegionEstimates, type = "latex",align=paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(NCarRegionEstimates),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("NCountiesCarnivoreRegions.tex",sep="")))

## ------       2.4.2 LAST YEAR N PER SEX PER COUNTY  ------
NCountyEstimatesLastRegions <- matrix("", ncol=3, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste(round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                   " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                   round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste(round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste(round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}



#print csv
write.csv(NCountyEstimatesLastRegions,
          file = file.path(WDTables,paste("NLastYearPerSex.csv",sep="")),fileEncoding="latin1")



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

idcountySWE1 <- idcountySWE
idcountySWE1 <- sort(idcountySWE1)

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"


row.names(NCountyEstimatesLastRegions) <- idcounty1
NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),], "*}", sep="")

NCountyEstimatesLastRegions["SWEDEN",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",], "**}", sep="")
NCountyEstimatesLastRegions["TOTAL",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",], "**}", sep="")
NCountyEstimatesLastRegions["Nordre",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",], "**}", sep="")


row.names(NCountyEstimatesLastRegions) <- c("TOTAL",
                                            paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                            paste("\\hspace{0.5cm} ",
                                                  idcountyNOR,sep=""),
                                            paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                            paste("\\hspace{0.5cm}","Norra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                            paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                            paste("\\hspace{0.5cm}","Södra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCountyEstimatesLastRegions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLastRegions))] <- paste("\\hspace{0.75cm}",
                                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


# WRITE LATEX 
print(xtable(NCountyEstimatesLastRegions, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path(WDTables,paste("NCountiesSexLastYearRegions.tex",sep="")))



## ------       2.4.2 LAST YEAR N PER SEX PER COUNTY WITH PROPORTION OF AREA COVERED  ------
NCountyEstimatesLastRegions <- matrix("", ncol=4, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total","\\% Area")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste(round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                   " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                   round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste(round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste(round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}



NCountyEstimatesLastRegions[names(areaAllRegions),"\\% Area"] <- round(areaAllRegions*100,digits = 0)
NCountyEstimatesLastRegions[NCountyEstimatesLastRegions[,4] %in% c("98","99"),4] <- 100
#NCountyEstimatesLastRegions[names(AreaStudied/TotalArea[names(AreaStudied)]),"Area"] <- round(AreaStudied/TotalArea[names(AreaStudied)]*100,digits = 0)

#print csv
write.csv(NCountyEstimatesLastRegions,
          file = file.path(WDTables,paste("NLastYearPerSexArea.csv",sep="")),fileEncoding="latin1")



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

idcountySWE1 <- idcountySWE
idcountySWE1 <- sort(idcountySWE1)

# idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"


row.names(NCountyEstimatesLastRegions) <- idcounty1
# NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),1:3], "*}", sep="")
# 
# NCountyEstimatesLastRegions["SWEDEN",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["TOTAL",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["Nordre",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",1:3], "**}", sep="")


row.names(NCountyEstimatesLastRegions) <- c("TOTAL",
                                            paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                            paste("\\hspace{0.5cm} ",
                                                  idcountyNOR,sep=""),
                                            paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                            paste("\\hspace{0.5cm}","Norra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                            paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                            paste("\\hspace{0.5cm}","Södra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCountyEstimatesLastRegions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLastRegions))] <- paste("\\hspace{0.75cm}",
                                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


# WRITE LATEX 
print(xtable(NCountyEstimatesLastRegions, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions)-1),"||c"),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path(WDTables,paste("NCountiesSexLastYearRegionsArea.tex",sep="")))




## ------       2.4.3 MAKE A TABLE 2 last years  ------
NCountyEstimatesLast2Regions <- matrix("", ncol=6, nrow=length(idcountyTable))
row.names(NCountyEstimatesLast2Regions) <- c(idcountyTable)
colnames(NCountyEstimatesLast2Regions) <- c(paste("Females", years[nYears-1]),
                                            paste("Males", years[nYears-1]),
                                            paste("Total", years[nYears-1]),
                                            paste("Females", years[nYears]),
                                            paste("Males", years[nYears]),
                                            paste("Total", years[nYears])
                                            
)




## FILL IN TABLE 
for(t in (nYears-1):nYears){
  ## FEMALES
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Females",years[t])] <- paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                      " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Males",years[t])] <- paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Total",years[t])] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


#print csv
write.csv(NCountyEstimatesLast2Regions,
          file = file.path(WDTables,paste("NLast2YearsPerSex.csv",sep="")),fileEncoding="latin1")



##ADD LITLE STAR TO NORRBOTTEN
idcountySWE1 <- idcountySWE

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
row.names(NCountyEstimatesLast2Regions) <- idcounty1
NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),] <- paste(NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),], "*", sep="")





row.names(NCountyEstimatesLast2Regions) <- c("TOTAL",
                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                             paste("\\hspace{0.5cm} ",
                                                   idcountyNOR,sep=""),
                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                             paste("\\hspace{0.5cm}","Norra**",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                             paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                             paste("\\hspace{0.5cm}","Södra",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)
row.names(NCountyEstimatesLast2Regions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLast2Regions))] <- paste("\\hspace{0.75cm}",
                                                                                                                    "VÃ¤stra GÃ¶taland", sep="")


NCountyEstimatesLast2Regions <- rbind(c("F","M","Total","F","M","Total"), NCountyEstimatesLast2Regions)

# WRITE LATEX 
addtorow <- list()

addtorow$pos <- list(c(0),0)
uniqueYEAR <- c(paste(unlist(YEARS[nYears-1]),collapse = "/"),paste(unlist(YEARS[nYears]),collapse = "/"))
addtorow$command <- c(paste0(paste0('& \\multicolumn{3}{c}{', uniqueYEAR,
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

print(xtable(NCountyEstimatesLast2Regions, type = "latex",
             align = paste(c("l",rep("c",3),"|",rep("c",3)),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      file = file.path(WDTables,paste("NCountiesSexLast2YearsRegions.tex",sep="")))

## ------       2.4.4 ALL YEARS N PER SEX PER COUNTY  ------
NCountyEstimatesAllSexRegions <- matrix("", ncol=nYears*3, nrow=length(idcountyTable)+1)
row.names(NCountyEstimatesAllSexRegions) <- c("",idcountyTable)
colnames(NCountyEstimatesAllSexRegions) <- rep(unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/")))),each=3)
NCountyEstimatesAllSexRegions[1,] <- rep(c("Females","Males","Total"),nYears)
## FILL IN TABLE 
for(t in 1:nYears){
  
  ## FEMALES
  cols <- which(colnames(NCountyEstimatesAllSexRegions) %in% unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/"))))[t])
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Females")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## TOTAL 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}

#print csv
write.csv(NCountyEstimatesAllSexRegions,
          file = file.path(WDTables,paste("NAllYearsPerSex.csv",sep="")),fileEncoding="latin1")



# # ADJUST NAMES OF THE TABLE 
# idcounty1 <- idcountyTable
# 
# idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
# idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
# idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"
# 
# idcountySWE1 <- idcountySWE
# idcountySWE1 <- sort(idcountySWE1)
# 
# idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"
# row.names(NCountyEstimatesLastRegions) <- idcounty1
# NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),], "*}", sep="")
# 
# 
# row.names(NCountyEstimatesLastRegions) <- c("TOTAL**",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN**",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1, sep="")
# )
# 
# ## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# # idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# # idcounty1 <- str_remove(idcounty1, "s ")
# # idcounty1 <- str_remove(idcounty1, " ")
# 
# 
# # WRITE LATEX 
# print(xtable(NCountyEstimatesLastRegions, type = "latex",
#              align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))),collapse = "")),
#       sanitize.text.function=function(x){x},
#       # scalebox=.8,
#       floating = FALSE,
#       add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
#       file = file.path(WDTables,paste("NCountiesSexLastYearRegions.tex",sep="")))
## ------       2.4.5 ALL YEARS, BOTH SEX COUNTIES NORWAY ------
idcounty <- row.names(DensityCountriesRegionsNOR[[t]]$summary)
#REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
#GET NORWEGIAN VERSUS SWEDISH COUNTIES 
idcountyNOR <- idcounty[grep("Region",idcounty)]
#idcountySWE <- sort(idcounty[-grep("Region",idcounty)])
idcountyTable <- c("Total","Norway", idcountyNOR, "Sweden" ,idcountySWE)


idcountyTable <- c("Total",
                   idcountyNOR
)


#CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimatesNOR <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimatesNOR) <- c(idcountyTable)
colnames(NCarRegionEstimatesNOR) <- unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/"))))#  unlist(lapply(YEARS ,function(x) x[2]))#

#FILL IN THE TABLE 
for(t in 1:nYears){
  for( i in 1:length(idcountyTable)){
    NCarRegionEstimatesNOR[idcountyTable[i],t] <- paste(round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                        " (",round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                        round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}



##QUICK CHECK TO MAKE SURE VALUES SUMS UP 
tmp <- DensityCountriesRegionsNOR[[t]]$summary#[1:(nrow(DensityCountriesRegions[[t]]$summary)),]
# SWE
row.names(DensityCountriesRegions[[t]]$summary)


# sum(tmp[row.names(tmp) %in% idcountySWE,"mean"])
sum(tmp[row.names(tmp) %in% idcountyNOR,"mean"])

# tmp["Sweden","mean"]
# #NOR
# sum(tmp[idcountyNOR,"mean"])
# tmp["Norway","mean"]
# #TOTAL
# sum(tmp[c(idcountyNOR,idcountySWE),"mean"])
# tmp["Total","mean"]
# 


## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1 <- gsub("Region ", "", idcountyTable)
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")

# row.names(NCarRegionEstimates) <- idcounty1
# NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste(NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*", sep="")
row.names(NCarRegionEstimatesNOR) <- idcounty1

#print csv
# NCarRegionEstimatesNOR <- data.frame(NCarRegionEstimatesNOR)
# NCarRegionEstimatesNOR$name <- row.names(NCarRegionEstimatesNOR)
# Encoding(NCarRegionEstimatesNOR[1,"name"]) <- "UTF-16"#"UTF-16"
#save(NCarRegionEstimatesNOR,file=file.path(WDTables,paste("NAllYearsNorwegianCounties.RData",sep="")))
write.csv(NCarRegionEstimatesNOR,
          file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")),fileEncoding= "latin1")

# Encoding(NCarRegionEstimatesNOR[,"name"])[9] <- "ISO-8859-1"
# mb_convert_encoding($file, 'UTF-8', 'ISO-8859-1')
# write.csv2(NCarRegionEstimatesNOR,
#            file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")),fileEncoding= "UTF-16LE")
# readr::write_excel_csv(NCarRegionEstimatesNOR,
#                         file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")))
## try to join to the Norwegian layer for Richard
# tmp <- data.frame(NCarRegionEstimatesNOR)
# tmp$NAME_1 <- row.names(tmp) 
# COUNTIES_1 <- merge(COUNTIES,tmp[,c("X2024","NAME_1")],by="NAME_1")
# 




row.names(NCarRegionEstimatesNOR) <- c("NORWAY",
                                       # paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                       paste("\\hspace{0.25cm} ",
                                             idcountyNOR,sep="")#
                                       # paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                       # paste("\\hspace{0.5cm}","Norra",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                       # paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                       # paste("\\hspace{0.5cm}","Södra",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCarRegionEstimatesNOR)[grep("VÃ¤straGÃ¶taland", row.names(NCarRegionEstimatesNOR))] <- paste("\\hspace{0.25cm}",
                                                                                                        "VÃ¤stra GÃ¶taland", sep="")

# row.names(NCarRegionEstimates) <- c("TOTAL",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1,sep="")
# )


print(xtable(NCarRegionEstimatesNOR, type = "latex",align=paste(c("l",rep("c",ncol(NCarRegionEstimatesNOR))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(NCarRegionEstimatesNOR),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("NCountiesCarnivoreRegionsNorway.tex",sep="")))









## ------     2.5 PLOT ABUNDANCE TIME SERIES ------
## ------       2.5.1 VIOLINS ------
## ------         2.5.1.1 ALL YEARS ------
SeasonText <- lapply(YEARS,FUN = function(x) paste(x,collapse = "/"))

#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf(file= file.path(WDFigures, paste("NCountriesBars.pdf", sep="")), width = 14, height = 10)
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.2,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 
# n.detected <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv",sep="")))
# n.detected <- as.vector(n.detected[1,2:ncol(n.detected)])
# n.detected[1]
# 
# for(t in 1:nYears){
#   xx <- t
#   yy <- n.detected[1,t]
#   xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
#   yy<-c(0,0,yy,yy)
#   polygon(xx, yy ,border=NA,col=grey(0.9))
#   
# }
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE

quantile95Tot <- quantile50Tot <- list()
quantile95Swe <- quantile50Swe <- list()
quantile95Nor <- quantile50Nor <- list()

#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95Tot[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Tot[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95Tot[[t]][1], quantile95Tot[[t]][1],
                quantile95Tot[[t]][2], quantile95Tot[[t]][2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95Tot[[t]][2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50Tot[[t]][1], quantile50Tot[[t]][1],
                  quantile50Tot[[t]][2], quantile50Tot[[t]][2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
  quantile95Swe[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Swe[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95Swe[[t]] [1], quantile95Swe[[t]] [1],
                quantile95Swe[[t]] [2], quantile95Swe[[t]] [2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95Swe[[t]][2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50Swe[[t]][1], quantile50Swe[[t]][1],
                  quantile50Swe[[t]][2], quantile50Swe[[t]][2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  quantile95Nor[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Nor[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95Nor[[t]][1], quantile95Nor[[t]][1],
                quantile95Nor[[t]][2], quantile95Nor[[t]][2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50Nor[[t]][1], quantile50Nor[[t]][1],
                  quantile50Nor[[t]][2], quantile50Nor[[t]][2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(200, 200, 200)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()

##save it for Post-plotting 
save(quantile95Tot, quantile50Tot,
     quantile95Swe, quantile50Swe,
     quantile95Nor, quantile50Nor, file=file.path(WDFigures, paste("CICounties.RData", sep=""))
)

## ------         2.5.1.2 LAST YEAR ------
pdf(file= file.path(WDFigures, paste("NCountriesViolinsLastYear.pdf", sep="")), width = 12, height = 8)

par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.8, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(nYears-0.1, nYears+0.1), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolves"), xaxt="n")
axis(1, at=c(nYears), labels = SeasonText[nYears], cex.axis=1.2)
at = c(1:nYears)




#for(t in nYears){
# xx <- t
# yy <- n.detected[[nYears]]
# xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
# yy<-c(0,0,yy,yy)
# polygon(xx, yy ,border=NA,col=grey(0.9))

#}

t <- nYears
#TOTAL
plot.violins2(list(colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)),
              x = at[t],
              at= at[t],
              violin.width = 0.02,
              col = TotalColors,
              alpha = violin.alpha,
              border.col = TotalColors,
              add = T
              ,cex=2,median = FALSE)


# text(round(round(DensityCountriesRegions[[t]]$summary["Total","mean"],digits = 1)),x=t+0.1,
#      y=DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]+total.offset, cex=text.cex)
# 
#SWEDEN
plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]),
              x = at[t],
              at=at[t],
              violin.width = 0.02,
              col = country.colors[2],
              alpha = violin.alpha,
              border.col = country.colors[2],
              add = T
              ,cex=2,median = FALSE)
# text(round(round(DensityCountriesRegions[[t]]$summary["Sweden","mean"],digits = 1)),x=t-0.1,
#      y=DensityCountriesRegions[[t]]$summary["Sweden","95%CILow"] + SE.offset, cex=text.cex)
#NORWAY
#print(t)
plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]),
              x = at[t],
              at=at[t],
              violin.width = 0.02,
              col = country.colors[1],
              alpha = violin.alpha,
              border.col = country.colors[1],
              add = T,scale.width = FALSE
              ,cex=2,median = FALSE
              
)
# text(round(round(DensityCountriesRegions[[t]]$summary["Norway","mean"],digits = 1)),x=t+0.1,
#      y=DensityCountriesRegions[[t]]$summary["Norway","95%CIHigh"]+NO.offset,cex=text.cex)
# 


# }
box()
# abline(v=at[1:(nYears-1)]+0.5,lty=2)

#legend
par(xpd=TRUE)
legend(x = nYears+0.01, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(2, 2, 2),
       horiz = T,
       pch=c(16, 16, 16),
       col=c(country.colors, "black"),
       bty = 'n',
       cex = 1.2)
legend(x = nYears+0.01, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(0.5,0.5 , 0.5),
       horiz = T,
       pch=c(16,16,16),
       col=c("white", "white", "white"),
       bty = 'n',
       cex = 1.2)

dev.off()


## ------       2.5.2 BARS ------
## ------         2.5.2.1 ALL YEARS ------
SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

# SeasonText <- lapply(YEARS, FUN = function(x) x[2]) #paste(x[[2]]))


#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf(file= file.path(WDFigures, paste("NCountriesBars.pdf", sep="")), width = 12, height = 8)
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.4,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 
# n.detected <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv",sep="")))
# n.detected <- as.vector(n.detected[1,2:ncol(n.detected)])
# n.detected[1]
# 
# for(t in 1:nYears){
#   xx <- t
#   yy <- n.detected[1,t]
#   xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
#   yy<-c(0,0,yy,yy)
#   polygon(xx, yy ,border=NA,col=grey(0.9))
#   
# }
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(200, 200, 200)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()




## ------         2.5.2.2 ALL YEARS SEX ------
SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

# SeasonText <- lapply(YEARS, FUN = function(x) x[2]) #paste(x[[2]]))


#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf(file= file.path(WDFigures, paste("NCountriesBarsSex.pdf", sep="")), width = 18, height = 8)
par(mfrow=c(1,2), mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,800),
     xlab="", ylab = paste("Estimated number of Females"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.1,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(50, 50, 50)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}

##MALES 
#par(mfrow=c(1,2), mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,800),
     xlab="", ylab = paste("Estimated number of Males"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.1,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a starCOUNTRIES <- COUNTRIES %>%    group_by(ISO) %>%summarize()
  
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegionsM[[t]]$PosteriorRegions["Norway",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(600,600,650,650),
        col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(630, 630, 630)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()





## ------         2.5.1.3 LAST YEAR ------
pdf(file= file.path(WDFigures, paste("NCountriesBarsLastYear.pdf", sep="")), width = 12, height = 8)
plot(-1000, xlim=c(nYears-0.1, nYears+0.1), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolves"), xaxt="n")
axis(1, at=c(nYears), labels = SeasonText[nYears], cex.axis=1.2)
at = c(1:nYears)


# #for(t in nYears){
# xx <- t
# yy <- n.detected[1,nYears]
# xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
# yy<-c(0,0,yy,yy)
# polygon(xx, yy ,border=NA,col=grey(0.9))
# 
# #}
widthPolygon <- 0.01
widthPolygon1 <- 0.01
widthPolygon2 <- 0.01
violin.alpha <- 0.8
displayQuantiles50 <- FALSE
t <- nYears
#for(t in 1:nYears){
#TOTAL
tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t-widthPolygon, t+widthPolygon,
              t+widthPolygon, t-widthPolygon ),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(TotalColors, violin.alpha),
        border= NA)



if(displayQuantiles50){
  polygon(x = c(t-widthPolygon, t+widthPolygon,
                t+widthPolygon, t-widthPolygon ),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(TotalColors, violin.alpha),
          border= NA)
}
#SWEDEN
tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t, t-widthPolygon1*2,
              t-widthPolygon1*2, t),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(country.colors[2], violin.alpha),
        border= NA)
if(displayQuantiles50){
  polygon(x = c(t, t-widthPolygon1*2,
                t-widthPolygon1*2, t),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(country.colors[2], violin.alpha),
          border= NA)
}

#NORWAY
#print(t)
tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t, t+widthPolygon1*2,
              t+widthPolygon1*2, t),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(country.colors[1], violin.alpha),
        border= NA)
if(displayQuantiles50){
  polygon(x = c(t, t+widthPolygon2*2,
                t+widthPolygon2*2, t),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(country.colors[1], violin.alpha),
          border= NA)
}


#}
box()
abline(v=at[1:(nYears-1)]+0.5,lty=2)

#legend
par(xpd=TRUE)
legend(x = 1, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       #pt.cex = c(4, 4, 4),
       horiz = T,
       #pch=c(16, 16, 16),
       fill=c(country.colors, "black"),
       border=NA,
       bty = 'n',
       cex = 1.5)


dev.off()

## ------         2.5.1.4 COMPARE WITH OTHER REPORTS ------
# #INITIALIZE THE OBJECTS
# TOTSWE <- TOTNOR <- TOT <- list()
# TOTNORCIL <- TOTNORCIH <- list()
# TOTSWECIL <- TOTSWECIH <- list()
# TOTCIL <- TOTCIH <- list()
# 
# count <- 1
# # SUMMARY MATRIX THE MODELS
# yearsEnd <- matrix( c(#1,9,"28.F_2021","","WolfRuns20122021","2022",
#   #1,10,"34.F_2022","","WolfRuns20142023","2022",
#   #1,10,"34.F_2022","Snap","WolfRuns20142023","2022",
#   # 1,10,"35.F_2022","","WolfRuns20122022","2022",
#   1,10,"Plot53Cleaned","","2023","2023",
#   2,11,"Plot53Cleaned","","2024","2024"
#   
#   
# ),nrow=2,byrow =T
# 
# )
# 
# 
# 
# xx=1
# # 
# #loop through the model results
# for(xx in 1:nrow(yearsEnd)){
#   modelName <- paste(yearsEnd[xx,3],yearsEnd[xx,4],sep="")
#   
#   # #SET DIRECTORY WHERE WOLF FIGURES WILL BE STORED
#   WDFigures <- file.path("C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM",yearsEnd[xx,5],
#                          modelName)
#   #C:\Users\cymi\Dropbox (Old)\AQEG Dropbox\AQEG Team Folder\RovQuant\wolf\WolfRuns20122021\FIGURES\28.F_2021
#   # "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/WolfRuns20142023/Figures/36.F_2023_sf"
#   # "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/WolfRuns2014023/Figures/36.F_2023_sf"
#   ##TOTAL
#   L <- readLines(file.path(WDFigures,"Table", "NCountiesCarnivoreRegions.tex"))
#   #
#   LTot <- grep("TOTAL", L, value = TRUE)
#   
#   DF <- read.table(text = LTot, sep = "&", header = F,
#                    strip.white = TRUE, check.names = FALSE, comment.char = "\\")
#   
#   TOT[[xx]] <- as.numeric(unlist(lapply(2:ncol(DF), function(x){
#     strsplit(DF[,x],split = " ")[[1]][1]
#   })))
#   #CI
#   CI <- unlist(lapply(2:ncol(DF), function(x){
#     strsplit(DF[,x],split = " ")[[1]][2]
#   }))
#   CI <- gsub("\\(", "", CI)
#   CI <- gsub("\\)", "", CI)
#   
#   TOTCIL[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[1])))
#   TOTCIH[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[2])))
#   
#   
#   ##SWEDEN
#   LSWE <- grep("SWEDEN", L, value = TRUE)
#   LSWE <- gsub("  \\\\hspace", "", LSWE)
#   LSWE <- gsub("  \\\\rowcolor", "", LSWE)
#   LSWE <- gsub("\\\\hspace", "", LSWE)
#   
#   
#   
#   DFSWE <- read.table(text = LSWE, sep = "&", header = F,
#                       strip.white = TRUE, check.names = FALSE, comment.char = "\\")
#   DFSWE[,1]
#   
#   TOTSWE[[xx]] <- as.numeric(unlist(lapply(2:ncol(DFSWE), function(x){
#     strsplit(DFSWE[,x],split = " ")[[1]][1]
#   })))
#   #CI
#   CI <- unlist(lapply(2:ncol(DFSWE), function(x){
#     strsplit(DFSWE[,x],split = " ")[[1]][2]
#   }))
#   CI <- gsub("\\(", "", CI)
#   CI <- gsub("\\)", "", CI)
#   
#   TOTSWECIL[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[1])))
#   TOTSWECIH[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[2])))
#   
#   
#   ##NORWAY
#   LNOR <- grep("NORWAY", L, value = TRUE)
#   
#   LNOR <- gsub("  \\\\rowcolor", "", LNOR)
#   LNOR <- gsub("\\\\hspace", "", LNOR)
#   LSWE <- gsub("\\\\hspace", "", LSWE)
#   
#   DFNOR <- read.table(text = LNOR, sep = "&", header = F,
#                       strip.white = TRUE, check.names = FALSE, comment.char = "\\")
#   DFNOR[,1]
#   
#   
#   TOTNOR[[xx]] <- as.numeric(unlist(lapply(2:ncol(DFNOR), function(x){
#     strsplit(DFNOR[,x],split = " ")[[1]][1]
#   })))
#   
#   #CI
#   CI <- unlist(lapply(2:ncol(DFNOR), function(x){
#     strsplit(DFNOR[,x],split = " ")[[1]][2]
#   }))
#   CI <- gsub("\\(", "", CI)
#   CI <- gsub("\\)", "", CI)
#   
#   TOTNORCIL[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[1])))
#   TOTNORCIH[[xx]] <- unlist(lapply(strsplit(CI ,split = "-"), function(x) as.numeric(x[2])))
#   # xx <- count +1
# }
## ------       2.5.3 MAPS ------
habbdensCropped <- list()
max <- max(unlist(lapply(DensityCountriesRegions, function(x) max(x$MeanCell))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

#PLOT
pdf(file=file.path(WDFigures, paste("DensityMapsAC5kms.pdf",sep="")))
for(t in 1: nYears){
  habbdens <- densityInputRegions$regions.r
  habbdens[] <- NA
  habbdens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  habbdensCropped[[t]] <- habbdens#crop(habbdens, e.sp)
  
  plot(habbdensCropped[[t]], breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
  plot(nngeo::st_remove_holes(myHabitat.list$habitat.poly),add=T, col=NA,border=grey(0.5))
  # points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],],
  #        pch=16, cex=0.4, col=adjustcolor("black",alpha.f = 0.2))
  plot(habbdensCropped[[t]], legend.only=TRUE,breaks=cuts, col=col,
       legend.width = 2,
       axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
                      labels=round(seq(0, max, length.out = 5),digits = 1),
                      cex.axis=0.6),
       legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
  
}
dev.off()

## ------   3. UD BASED DENSITY  (5km) ------
### IDENTIFY PROXIMITY HABITAT CELLS 
habitatMask <- densityInputCountries$habitat.id
habitatMask[!is.na(habitatMask)] <- 1 
# DetIndex <- getLocalObjects(habitatMask =  habitatMask,
#                             coords = densityInputCountries$habitat.xy,
#                             dmax = 15,
#                             resizeFactor = 1)

### COMPUTE THE UD BASED FOR A FEW ITERATIONS
#RESCALE SIGMA TO METERS 
sigma <- myResultsSXYZ_MF$sims.list$sigma#*res(myHabitat.list$habitat.r)[1]
#RESCALE SIGMA TO THE HABITAT SCALE
sigmaRescaled <- sigma/res(rrCountries)[1]

#SELECT xxx ITERATIONS RANDOMLY 
# spaceUSED <- list()
# iter <- sample(1:dim(densityInputCountries$sy)[1], size = 1000)#dim(densityInputCountries$sx)[1])
# for(t in 1:nYears){
#   # spaceUSED[[t]] <- GetSpaceUseLESS( densityInputCountries$sx[iter,,t],
#   #                               densityInputCountries$sy[iter,,t],
#   #                               myResultsSXYZ_MF$sims.list$z[iter,,t],
#   #                               sigmaRescaled[iter],
#   #                               densityInputCountries$habitat.xy,
#   #                               aliveStates = alive.states,
#   #                               regionID = regionID,
#   #                               habitatID = DetIndex$habitatGrid-1,
#   #                               habitatIndex = DetIndex$localIndices-1,
#   #                               nHabitatLESS = DetIndex$numLocalIndices,
#   #                               display_progress = T,
#   #                               returnPosteriorCells = F
#   # )
#   gc()
#   spaceUSED[[t]] <- GetSpaceUse(densityInputCountries$sx[iter,,t],
#                                 densityInputCountries$sy[iter,,t],
#                                 myResultsSXYZ_MF$sims.list$z[iter,,t],
#                                 sigmaRescaled[iter],#sigmaRescaled[iter],
#                                 densityInputCountries$habitat.xy,
#                                 aliveStates = alive.states,
#                                 regionID = densityInputCountries$regions.rgmx,
#                                 display_progress = T,
#                                 returnPosteriorCells = T
#   )
# }
# 
# t=9
# plot(densityInputCountries$habitat.xy[,2]~
#        densityInputCountries$habitat.xy[,1],pch=16,cex=0.1)
# 
# for(iter in 1:100){
#   points(densityInputCountries$sy[iter,myResultsSXYZ_MF$sims.list$z[iter,,t] %in%2,t]~
#            densityInputCountries$sx[iter,myResultsSXYZ_MF$sims.list$z[iter,,t] %in%2,t],col="red",pch=16,cex=0.1)
# }
# 
# 
# 
# 
# 
# plot(densityInputCountries$sy[1,myResultsSXYZ_MF$sims.list$z[1,,t] %in%2,t]~
#        densityInputCountries$sx[1,myResultsSXYZ_MF$sims.list$z[1,,t] %in%2,t])
# 
# spaceUSED1 <- spaceUSED
# spaceUSED <- list()
# for(t in 1:nYears){
#   spaceUSED[[t]] <- list()
#   spaceUSED[[t]][["MeanCell"]] <- spaceUSED1[[t]]$MeanCell
# }
# save(spaceUSED, file = file.path(WDFigures, "spaceUsed5km.RData" ))
load(file = file.path(WDFigures, "spaceUsed5km.RData" ))

## ------          2.1.2.1 PLOT TIME SERIES ------
#### PLOT TIME SERIES 
## PREPARE THE FILES 
#SeasonText <- lapply(YEARS,FUN = function(x) paste(x,collapse = "/"))
SeasonText <- lapply(YEARS,FUN = function(x) paste(x[[2]]))
# habbdensUDCropped[[t]]
COUNTRIESSCA <- COUNTRIES[COUNTRIES$ISO %in% c("NOR","SWE"),]
COUNTRIESsimpFig <- st_simplify(COUNTRIESSCA, preserveTopology = F,dTolerance = 4000)
habbdensFig <- densityInputRegions$regions.r
habbRFig <- densityInputRegions$regions.r
## from 25km2 (5*5raster) to 100km2
spaceUSED100km2 <- lapply(spaceUSED, function(x) x$MeanCell * 4 )


# studDis <- disaggregate(myHabitat.list$habitat.poly)
# studDis$id <- 1:length(studDis)
# plot(studDis)
# text(studDis,studDis$id)
# lakes <- disaggregate(COUNTRIESWaterHumans[COUNTRIESWaterHumans$area>20000000000,])
#lakes <- dropHole(lakes)
# plot(lakes)

#find the lakes 
# which(unlist(lapply(lakes@polygons[[2]]@Polygons, function(x) x@area))>500000000)
# featureNumber=2 ; ringNumber=115
# Lake1 = SpatialPolygons(
#   list(
#     Polygons(
#       list(
#         lakes@polygons[[featureNumber]]@Polygons[[ringNumber]]
#       ),
#       ID=1)))
# featureNumber=2 ; ringNumber=103
# Lake2 = SpatialPolygons(
#   list(
#     Polygons(
#       list(
#         lakes@polygons[[featureNumber]]@Polygons[[ringNumber]]
#       ),
#       ID=1)))

##PLOT
pdf(file=file.path(WDFigures, paste("DensityMapsUD.pdf",sep="")), width = 12, height = 8)
#layout
mx <- rbind(c(1,rep(1:5, each=2)),
            c(rep(1:5, each=2),5))
mx <- rbind(mx, mx+5)
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)
habbdensUDCropped <- list()
for(t in 1:length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  habbdensFig[!is.na(habbRFig[])] <- spaceUSED100km2[[t]]
  habbdensFig[habbRFig[]==0] <- NA
  
  habbdensUDCropped[[t]] <- habbdensFig#mask(habbdensFig, e.sp)
  crs( habbdensUDCropped[[t]]) <- st_crs(myHabitat.list$habitat.poly)
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image(habbdensUDCropped[[t]], add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  # plot(RemoveHolesSp(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  
  mtext(SeasonText[[t]], 1, -2, adj=0.15, cex=1.2)
  #box()
  # spPol <- rasterToPolygons(myHabitat.list$habitat.r,
  #                  fun = function(x) x==1) 
  # spPol <- aggregate(spPol)
  # plot(SpatialPolygons(spPol@polygons[[1]]@Polygons[[1]]))
  # 
  # spPol <- st_as_sf(spPol)
  # state_union <- spPol %>% 
  #   group_by(Habitat) %>%
  #   summarise(geometry = sf::st_union(geometry)) %>%
  #   ungroup() %>% st_as_sf()
  # plot(state_union$geometry)
  # 
  # agg <- aggregate(rasterToPolygons(myHabitat.list$habitat.r,
  #                                   fun = function(x) x==1))
  # agg <- RemoveHolesSp(agg)
  # 
  # plot(agg, add=TRUE, border="black", col=NA)
  agg1 <- aggregate(rasterToPolygons(myHabitat.list$habitat.rWthBuffer,
                                     fun = function(x) x==1))
  #plot(agg1, add=TRUE, border="red", col=NA)
  # plot(RemoveHolesSp(aggregate(myHabitat.list$habitat.poly)), add=TRUE, border=grey(0.5), col=NA)
  
  
  
  # plot(Lake2, add=TRUE, border=grey(0.4), col=NA)
  
  
  if(sum(t%in%yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1320000,x1=1320000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=2, lend=2)  
    text(1280000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    plot(habbdensUDCropped[[t]], legend.only=TRUE,breaks=cuts, col=col,
         legend.width = 2,
         axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        cex.axis=1.6),
         smallplot=c(0.72, 0.75, 0.2, 0.4),
         legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                          side=4, font=1, line=4, cex=1.2))
  }
  
}
dev.off()

## ------          2.1.2.2 PLOT LAST 2 YEARS ------
## PLOT LAST 2 YEARS  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLast2Years.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1:2, each=2)),
            c(rep(1:2, each=2), 3))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in (length(years)-1):length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = NA)
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  mtext(SeasonText[[t]], 1, -2, adj=0.15, cex=1.2)
  #box()
  plot(nngeo::st_remove_holes(myHabitat.list$habitat.poly), add=TRUE, border="black", col=NA)
  #plot(Lake2, add=TRUE, border=grey(0.4), col=NA)
  
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 5), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 5), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.72, 0.75, 0.2, 0.4),
          legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
}
dev.off()


## ------          2.1.2.3 PLOT LAST YEAR ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYear.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)#YEARS[[t]][2]
  #box()
  #plot(RemoveHolesSp(aggregate(myHabitat.list$habitat.poly)), add=TRUE, border="black", col=NA)
  #plot(Lake2, add=TRUE, border=grey(0.4), col=NA)
  
  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 0.5,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.73, 0.75, 0.25, 0.4),
          legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()


## ------          2.1.2.4 PLOT LAST YEAR SUMMARY ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYearSummary.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  # plot(COUNTRIESsimpFig[1], border=grey(0.1), col = NA, add=TRUE, lwd=2)
  
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)
  #box()
  #plot(RemoveHolesSp(aggregate(myHabitat.listF$habitat.poly)), add=TRUE, border="black", col=NA)
  #plot(Lake2, add=TRUE, border="black", col=NA)
  
  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.72, 0.75, 0.2, 0.4),
          legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()



## ------          2.1.2.5 PLOT LAST YEAR SUMMARY NO ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYearSummaryNO.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  # plot(COUNTRIESsimpFig[1], border=grey(0.1), col = NA, add=TRUE, lwd=2)
  
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)
  #box()
  #plot(RemoveHolesSp(aggregate(myHabitat.listF$habitat.poly)), add=TRUE, border="black", col=NA)
  #plot(Lake2, add=TRUE, border="black", col=NA)
  
  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.72, 0.75, 0.2, 0.4),
          legend.args=list(text=expression(paste("Individer/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()




## ------          2.1.2.6 WRITE UD 5km RASTER FOR ROVBASE ------
if(!dir.exists(file.path(WDFigures, "RasterForRovbase"))){dir.create(file.path(WDFigures, "RasterForRovbase"))}

for(t in 1:length(years)){
  raster::crs(habbdensUDCropped[[t]]) <- "EPSG:32633"#st_crs(myHabitat.list$habitat.poly))
  
  
  path <- file.path(WDFigures, "RasterForRovbase",paste("wolverine_5km",paste(YEARS[[t]][1],collapse = "_"),".tif",sep=""))
  writeRaster(habbdensUDCropped[[t]], path, overwrite=TRUE)
}


## ------  4. DERIVED PARAMETERS FROM ABUNDANCE ------ 
## ------    4.1 MAKE A GROWTH RATE TABLE PER COUNTRY  ------
growthRate <- matrix(0, ncol=nYears-1,nrow=3)
row.names(growthRate) <- c("Norway","Sweden","Total")
colnames(growthRate) <- unlist(lapply(YEARS[2:(length(YEARS))],function(x) paste(x,collapse =  "-")))

for(t in 1:(nYears-1)){
  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Norway",] /
    DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  growthRate["Norway",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                  " (",
                                  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                  ")",sep="")
  
  
  
  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Sweden",]/ 
    (DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",])
  growthRate["Sweden",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                  " (",
                                  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                  ")",sep="")
  
  
  growth <- colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"),])/
    colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])#colSums(DensityCountriesRegions[[t+1]]$PosteriorAllRegions) / 
  #colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)
  growthRate["Total",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                 " (",format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                 format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                 ")",sep="")
  
}


### add *** to years and country with OPSCR results 
#row.names(growthRate)[2] <- "Sweden*"
#row.names(growthRate)[3] <- "Total*"
# transitionNotSampled <- colnames(growthRate)
# transitionNotSampled <- unique(unlist(lapply(1:length(years),function(x) grep(as.character((years+1)[yearsNotSampled])[x], transitionNotSampled))))
# transitionNotSampled <- transitionNotSampled[!is.na(transitionNotSampled)] 
# growthRate[2,transitionNotSampled] <- paste(growthRate[2,transitionNotSampled],"*",sep="")
# growthRate[3,transitionNotSampled] <- paste(growthRate[3,transitionNotSampled],"*",sep="")


#print table 
print(xtable(growthRate, type = "latex", align=paste(c("l", rep("c",ncol(growthRate))), collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(growthRate),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("growthRate.tex", sep="")))


colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"),])
## ------    4.2 DERIVE SEX RATIO  ------
PropFemale <- PropFemaleSWE <- PropFemaleNOR <- list()

for(t in 1:nYears){
  PropFemale[[t]] <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])/
    (colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]) + colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]))
  
  PropFemaleSWE[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",])
  
  PropFemaleNOR[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Norway",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",])
}

#OVERALL PROPORTION OF FEMALES
mean(unlist(PropFemale))
round(quantile(unlist(PropFemale), probs=c(0.025,0.975)),digits = 2)

#median(PropFemale[[nYears]])
quantile(PropFemale[[nYears]], probs=c(0.025,0.975))
mean(PropFemale[[nYears]])


# mean(do.call(c,PropFemale))
# quantile(do.call(c,PropFemale), probs=c(0.025,0.975))

propFemale_tab <- matrix(0, ncol=nYears,nrow=3)
row.names(propFemale_tab) <- c("Norway","Sweden","Total")
colnames(propFemale_tab) <- unlist(lapply(YEARS,function(x) c(x[2])))#c(unlist(lapply(YEARS[1:(length(YEARS))],function(x) paste(x,collapse =  "-"))))
for(t in 1:nYears){
  #Sweden
  propFemale_tab["Sweden",t] <- paste(format(round(mean(PropFemaleSWE[[t]]),digits = 2),nsmall = 2),
                                      " (",
                                      format(round(quantile(PropFemaleSWE[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                      format(round(quantile(PropFemaleSWE[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                      ")",sep="")
  #NORWAY
  propFemale_tab["Norway",t] <- paste(format(round(mean(PropFemaleNOR[[t]]),digits = 2),nsmall = 2),
                                      " (",
                                      format(round(quantile(PropFemaleNOR[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                      format(round(quantile(PropFemaleNOR[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                      ")",sep="")
  propFemale_tab["Total",t] <- paste(format(round(mean(PropFemale[[t]]),digits = 2),nsmall = 2),
                                     " (",
                                     format(round(quantile(PropFemale[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                     format(round(quantile(PropFemale[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                     ")",sep="")
}

#print table 
print(xtable(propFemale_tab, type = "latex", align=paste(c("l", rep("c",ncol(propFemale_tab))), collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(propFemale_tab),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("propFemale.tex", sep="")))


## ------    4.3 DERIVE DENSITY  ------
habbRCarRegionsTRY <- rrCountries
habbRCarRegionsTRY[!is.na(habbRCarRegionsTRY[])] <-1 
habbRCarRegionsTRY[] <- as.numeric(habbRCarRegionsTRY[])

habbRCarRegionsTRYpol <- sf::st_as_sf(stars::st_as_stars(habbRCarRegionsTRY), 
                                      as_points = FALSE, merge = F)

plot(rasterToPolygons(habbRCarRegionsTRY, function(x) x>0,dissolve = T))
areaSqKm <- sum(st_area(habbRCarRegionsTRYpol))#*1e-6
units::set_units(areaSqKm, km^2)
#size of the area in km2
areaSqKm

#multiplied by 100 to get per 100km2
DensityCountriesRegions[[t]]$summary["Total","mean"]/areaSqKm*100
DensityCountriesRegions[[t]]$summary["Total","95%CILow"]/areaSqKm*100
DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]/areaSqKm*100

###COULD REMAKE ALL THE TABLES WITH DENSITY INSTEAD OF ABUNDANCE....
## ------    4.4  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED ------
n.detectedTotal <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv", sep="")))
n.detectedTotal <- as.vector(n.detectedTotal[1, 2:ncol(n.detectedTotal)])

n.detectedSex <- read.csv(file.path(WDTables, paste("NGSidCountrySEX.csv", sep="")))
tF <- seq(2, length.out = nYears, by=2)
tM <- seq(3, length.out = nYears, by=2)

propDetected <- matrix("", ncol=nYears,nrow=3)
row.names(propDetected) <- c("M","F","Total")
colnames(propDetected) <- unlist(lapply(YEARS,function(x) c(x[2])))#unlist(lapply(YEARS,function(x) paste(x,collapse = "/")))#

for(t in 1:nYears){
  propDetected["Total",t] <-  paste(format(round(mean(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                    " (",format(round(quantile(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                    format(round(quantile(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
  
  n.detectedSexM <- as.numeric(as.character(n.detectedSex[4,tM[t]]))
  propDetected["M",t] <-  paste(format(round(mean(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                " (",format(round(quantile(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.025)), digits = 2),nsmall = 2), "-",
                                format(round(quantile(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.975)), digits = 2),nsmall = 2),")", sep="")
  
  n.detectedSexF <- as.numeric(as.character(n.detectedSex[4,tF[t]]))
  propDetected["F",t] <-  paste(format(round(mean(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                " (",format(round(quantile(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.025)), digits = 2),nsmall = 2), "-",
                                format(round(quantile(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.975)), digits = 2),nsmall = 2),")", sep="")
  
  
  
}


##PRINT
print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(propDetected), by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("PropDetectedIds.tex",sep="")))




## ------    4.5  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED PER COUNTRIES  ------
n.detectedCountry <- read.csv(file.path(WDTables, paste("NGSidCountrySEX.csv", sep="")))
colnames(n.detectedCountry) <- c("",unlist(lapply(YEARS,function(x) c(x[2],x[2])))) 
#propDetected <- matrix("", ncol=nYears,nrow=3)
#row.names(propDetected) <- c("M","F","Total")
#colnames(propDetected) <- unlist(lapply(YEARS,function(x) c(x[2])))#unlist(lapply(YEARS,function(x) paste(x,collapse = "/")))#
propDetectedCountry <- n.detectedCountry 
propDetectedCountry[2:4, 2:ncol(propDetectedCountry)] <- NA

NCountrySex <- read.csv(file.path(WDTables, paste("NAllYearsPerSex.csv", sep="")))
colnames(NCountrySex) <- c("",unlist(lapply(YEARS,function(x) c(x[2],x[2],x[2]))))

yrs <- unlist(lapply(YEARS, function(x) c(x[2]))) 


lisCountries <- list()

lisCountries[[1]] <- c("Norway")
lisCountries[[2]] <- c("Sweden")
lisCountries[[3]] <- c("Sweden","Norway")

for(t in 1:nYears){
  col <- which(colnames(NCountrySex)  %in% yrs[t])
  NCountrySex[,col]
  
  cols <- which(colnames(propDetectedCountry)  %in% as.character(yrs[t]))
  for(p in 1:2){
    propDetectedCountry[p+1,cols[1]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[1]]) /DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[p]],]), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
    propDetectedCountry[p+1,cols[2]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[2]]) /DensityCountriesRegionsM[[t]]$PosteriorRegions[lisCountries[[p]],]), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
  }
  for(p in 3:3){
    propDetectedCountry[p+1,cols[1]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[1]]) /colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[p]],])), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
    propDetectedCountry[p+1,cols[2]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[2]]) /colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),])), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
  }
  
  
  
}



addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(propDetectedCountry)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                    '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))


# colnames(TableState) <- rep("", ncol(TableState))
# REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC


# multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
# multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
# DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))
# addtorow$pos <- list(c(0),0)
# addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
#                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))

# xTableState <- xtable(TableState)

# rownames(TableState)[2:5] <- c("$\\rho$","$\\phi$","h","w")

# rownames(TableState)[2:8] <- c("$\\gamma$","$\\phi$","   ","h","  ","w","  ")


print(xtable(propDetectedCountry, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry)+1), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry.tex", sep="")))





##SPLIT THE TABLE IN TWO 
propDetectedCountry1 <- propDetectedCountry[,c(1:11)]
propDetectedCountry2 <- propDetectedCountry[,c(1,12:21)]

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[2:11])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[12:21])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))


#SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(propDetectedCountry1, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(propDetectedCountry2, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry2.tex", sep="")))



## ------  5. P0  ------
widthPolygon <- 0.15
## ------    5.1 BARS  ------
pdf(file=file.path(WDFigures, paste("DetectionProbBars.pdf",sep="")),width=8,height=10)
# mx <- matrix(1:24,12,2,byrow = TRUE)
mx <- matrix(c(1,2,5,6,9,10,13,14),4,2,byrow = T)
mx1 <- matrix(c(3,4,7,8,11,12,15,16),4,2,byrow = T)
mx <- cbind(mx,mx1)
# mx <- cbind((max(mx)+1):((max(mx)+1)+dim(mx)[1]-1),mx)
# mx <- rbind((max(mx)+1):((max(mx)+1)+dim(mx)[2]-1),mx)
# 
# nf <- layout(mx,widths=c(0.2,rep(1,dim(mx)[2]-2),0.75),heights=c(0.2,rep(1,dim(mx)[1]-1)))
nf <- layout(mx,widths=c(1,0.5,1,0.5),heights=c(rep(1,dim(mx)[1]-1)))
# layout.show(nf)
country.colors <- c("goldenrod1","goldenrod3")
contry.colors.samples <- c("black",grey(0.4))
names(country.colors) <- c("Norway","Sweden")
COUNTIES_AGGREGATEDSubsetsimp <- st_simplify(COUNTIES_AGGREGATEDSubset, dTolerance = 2000, preserveTopology  = T)
COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique

# ggplot(COUNTIES_AGGREGATEDSubsetsimp) +
#   geom_sf(aes(fill = idunique)) +
#   geom_sf_label(aes(label = idunique))

COUNTIES_AGGREGATEDSubsetsimp$Name <- c("NO5",
                                        "NO4",
                                        "SE3",
                                        "NO3",
                                        "SE2",
                                        "NO2",
                                        "SE1",
                                        "NO1")#,"")

# names(detCounties.original) <- c("NO1","SE1","SE2","NO2","SE3","SE4")
# county.index.by.country <- (1:length(detCounties.original))[rev(order(names(detCounties.original)))]

myCol <-rep(c("lightblue","blue"),15)
myDev <- c(-0.15,+0.15)#seq(-0.3,0.3,length.out = 2)# 

CountyIndex <- COUNTIES_AGGREGATEDSubsetsimp$idunique

# CountyIndex <- unique(unlist(detCounties.original1))
# COUNTIESList <- list(F=COUNTIESF, M=COUNTIESM)



index <-c(6,1,
          4,5,
          7,4,
          2,3)
index <-c(4,6,
          5,3,
          7,2,
          8,1)
#for(c in 1:length(CountyIndex)){
for(c in index){
  
  par(mar=c(4,4,1,1), tck=0)
  # par(mfrow=c(1,2))
  plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,0.04), type ="n", xaxt="n",
       xlab = "Years", ylab = "Detection probability")
  axis(2,tck=-0.02)
  
  axis(1, c(1:nYears),labels = years+1,cex.axis=0.5,padj = -2)
  abline(v=1:(nYears-1)+0.5,lty=2,col=grey(0.5))
  
  # plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.1), type ="n", xaxt="n", xlab = "Years",
  #      ylab = "p0")#, main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
  # axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
  for(s in 1:2){
    myResults <- myResultsList[[s]]
    
    for(t in 1:nYears){
      # coun <- which( detCounties.original1[[t]]  %in%  CountyIndex[c])
      if(length(c)>0){
        tmp <- myResults$sims.list$p0[ , c, t]
        quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
        polygon(x = c(t+ myDev[s] - widthPolygon, t+ myDev[s] + widthPolygon,
                      t+ myDev[s] + widthPolygon, t+ myDev[s] - widthPolygon ),
                y = c(quantile95[1], quantile95[1],
                      quantile95[2], quantile95[2]), 
                col=adjustcolor(myCol[s], 0.6),
                border= adjustcolor(myCol[s], 0.7))
        
        
      }
    }
  }
  
  if(c==6){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.006)
    for(i in 1:2){
      #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
      points(2, 0.035-yoffset[i], pch=15,cex=1.8,col=adjustcolor(cols[i],violin.alpha95+0.1))
      #points(2, 0.035-yoffset[i], pch=15,cex=0.5,col=adjustcolor(pt.col[i],violin.alpha50))
      text(3.5, 0.035-yoffset[i], labels[i],cex=1.3,pos=4)
    }
  }
  
  
  
  ## PLOT REGION MAP    
  par(mar=c(0,0,0,0))
  plot(st_geometry(COUNTIES_AGGREGATEDSubsetsimp),border=grey(0.5),col=grey(0.5),lwd=0.1)
  aggCounty <- COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==c, ]
  plot(st_geometry(aggCounty),
       add=T, col=adjustcolor("red",0.4),border="red")
  # plot(COUNTIESList[["M"]][COUNTIESList[["M"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("red",0.4),border="red")
  
  # plot(COUNTIESList[["F"]][COUNTIESList[["F"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("blue",0.4))#,border=grey(0.5))
  # par(mar=c(4.5,1,1,1),xaxs="i",yaxs="i")
  # Cplot(COUNTIESplot,border=grey(0.5))
  text(st_coordinates((st_centroid(aggCounty))) ,labels=COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==c, ]$Name, col="white")
  
  
}



dev.off()

## ------    5.1 BARS LAST YEAR  ------
#layout
pdf(file=file.path(WDFigures, paste("DetectionProbBarsLastYear.pdf",sep="")),width=10,
    height=7)

mx <- matrix(c(1,2),1,2,byrow = T)
nf <- layout(mx,widths=c(1,0.8),
             heights=c(1))
#colors and plot params
country.colors <- c("goldenrod1","goldenrod3")
contry.colors.samples <- c("black",grey(0.4))
names(country.colors) <- c("Norway","Sweden")

myCol <-rep(c("lightblue","blue"),15)
myDev <- c(-0.15,+0.15)
#deal with the map and index
COUNTIES_AGGREGATEDSubsetsimp <- st_simplify(COUNTIES_AGGREGATEDSubset, 
                                             dTolerance = 2000, preserveTopology  = T)
COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
COUNTIES_AGGREGATEDSubsetsimp$Name <- c("NO5",
                                        "NO4",
                                        "SE3",
                                        "NO3",
                                        "SE2",
                                        "NO2",
                                        "SE1",
                                        "NO1")
COUNTIES_AGGREGATEDSubsetsimp <- COUNTIES_AGGREGATEDSubsetsimp[order(COUNTIES_AGGREGATEDSubsetsimp$Name),]
CountyIndex <- COUNTIES_AGGREGATEDSubsetsimp$idunique


#plot
par(mar=c(4,4,1,1), tck=0)
#plot
plot(10, xlim = c(0.5, length(index)+0.5), ylim = c(0,0.04), type ="n", xaxt="n",
     xlab = "Region", ylab = "Detection probability")
axis(2,tck=-0.02)
axis(1, c(1:length(CountyIndex)),
     labels = COUNTIES_AGGREGATEDSubsetsimp$Name,
     cex.axis=1)#,padj = -2)
abline(v=1:(length(CountyIndex)-1)+0.5,lty=2,col=grey(0.5))

count <- 0
for(c in 1:length(index)){
  count <- count+1
  
  
  # plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.1), type ="n", xaxt="n", xlab = "Years",
  #      ylab = "p0")#, main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
  # axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
  for(s in 1:2){
    myResults <- myResultsList[[s]]
    t <- nYears
    #for(t in 1:nYears){
    # coun <- which( detCounties.original1[[t]]  %in%  CountyIndex[c])
    if(length(c)>0){
      
      
      tmp <- myResults$sims.list$p0[ , CountyIndex[c], t]
      quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
      polygon(x = c(count+ myDev[s] - widthPolygon, count+ myDev[s] + widthPolygon,
                    count+ myDev[s] + widthPolygon, count+ myDev[s] - widthPolygon ),
              y = c(quantile95[1], quantile95[1],
                    quantile95[2], quantile95[2]), 
              col=adjustcolor(myCol[s], 0.6),
              border= adjustcolor(myCol[s], 0.7))
      
      
    }
  }
  
  if(c==6){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.003)
    for(i in 1:2){
      #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
      points(2, 0.035-yoffset[i], pch=15,cex=1.8,col=adjustcolor(cols[i],violin.alpha95+0.1))
      #points(2, 0.035-yoffset[i], pch=15,cex=0.5,col=adjustcolor(pt.col[i],violin.alpha50))
      text(2.5, 0.035-yoffset[i], labels[i],cex=1.3,pos=4)
    }
  }
}

## PLOT REGION MAP    
par(mar=c(0,0,0,0))
plot(st_geometry(COUNTIES_AGGREGATEDSubsetsimp),col=grey(0.5),border=NA)
plot(st_geometry(COUNTIES_AGGREGATEDSubsetsimp),
     add=T, col=NA,border="white",lwd=.5)
for(c in index){
  
  aggCounty <- COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==c, ]
  
  
  # plot(COUNTIESList[["M"]][COUNTIESList[["M"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("red",0.4),border="red")
  
  # plot(COUNTIESList[["F"]][COUNTIESList[["F"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("blue",0.4))#,border=grey(0.5))
  # par(mar=c(4.5,1,1,1),xaxs="i",yaxs="i")
  # Cplot(COUNTIESplot,border=grey(0.5))
  text(st_coordinates((st_centroid(aggCounty))),
       labels=COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==c, ]$Name, 
       col="black")
}

#}


dev.off()



## ------    5.1 BARS Other ------
pdf(file=file.path(WDFigures, paste("DetectionProbBarsOther.pdf",sep="")),
    width=8,height=8)
# mx <- matrix(1:24,12,2,byrow = TRUE)
mx <- matrix(1:4,2,2,byrow = TRUE)

# mx <- cbind((max(mx)+1):((max(mx)+1)+dim(mx)[1]-1),mx)
# mx <- rbind((max(mx)+1):((max(mx)+1)+dim(mx)[2]-1),mx)
# 
# nf <- layout(mx,widths=c(0.2,rep(1,dim(mx)[2]-2),0.75),heights=c(0.2,rep(1,dim(mx)[1]-1)))
nf <- layout(mx,widths=c(1,0.5), heights=c(rep(1,dim(mx)[1]-1)))

# layout.show(nf)
country.colors <- c("goldenrod1","goldenrod3")
contry.colors.samples <- c("black",grey(0.4))
names(country.colors) <- c("Norway","Sweden")
COUNTRIESsimp <- st_simplify(COUNTRIES,dTolerance = 1500,preserveTopology = T)
COUNTRIESsimp$idunique <- COUNTRIES$ISO

# names(detCounties.original) <- c("NO1","SE1","SE2","NO2","SE3","SE4")
# county.index.by.country <- (1:length(detCounties.original))[rev(order(names(detCounties.original)))]

myCol <-rep(c("lightblue","blue"),15)
myDev <- c(-0.15,+0.15)#seq(-0.3,0.3,length.out = 2)# 


# CountyIndex <- unique(unlist(detCounties.original1))
# COUNTIESList <- list(F=COUNTIESF, M=COUNTIESM)
for(c in 1:2){
  
  par(mar=c(4,4,1,1), tck=0)
  # par(mfrow=c(1,2))
  plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,0.015), type ="n", xaxt="n",
       xlab = "Years", ylab = "Detection probability")
  axis(2,tck=-0.02)
  
  axis(1, c(1:nYears),labels = years+1,cex.axis=0.9)
  abline(v=1:(nYears-1)+0.5,lty=2,col=grey(0.5))
  
  # plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.1), type ="n", xaxt="n", xlab = "Years",
  #      ylab = "p0")#, main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
  # axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
  for(s in 1:2){
    myResults <- myResultsList[[s]]
    
    for(t in 1:nYears){
      # coun <- which( detCounties.original1[[t]]  %in%  CountyIndex[c])
      if(length(c)>0){
        tmp <- myResults$sims.list$p0Oth[ , c, t]
        quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
        polygon(x = c(t+ myDev[s] - widthPolygon, t+ myDev[s] + widthPolygon,
                      t+ myDev[s] + widthPolygon, t+ myDev[s] - widthPolygon ),
                y = c(quantile95[1], quantile95[1],
                      quantile95[2], quantile95[2]), 
                col=adjustcolor(myCol[s], 0.6),
                border= adjustcolor(myCol[s], 0.7))
        
        
      }
    }
  }
  
  
  if(c==1){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.002)
    for(i in 1:2){
      #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
      points(2, 0.014-yoffset[i], pch=15,cex=2.1,col=adjustcolor(cols[i],violin.alpha95+0.1))
      #points(2, 0.035-yoffset[i], pch=15,cex=0.5,col=adjustcolor(pt.col[i],violin.alpha50))
      text(2.5, 0.014-yoffset[i], labels[i],cex=1.5,pos=4)
    }
  }
  
  par(mar=c(0,0,0,0))
  plot(st_geometry(COUNTRIESsimp),border=grey(0.5),col=grey(0.5),lwd=0.1)
  aggCounty <- nngeo::st_remove_holes(COUNTRIESsimp[c, ])
  
  #aggCounty <- RemoveHolesSp(aggregate(COUNTRIESsimp[c, ]))
  plot(st_geometry(aggCounty),
       add=T, col=adjustcolor("red",0.4),border="red")
  # plot(COUNTIESList[["M"]][COUNTIESList[["M"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("red",0.4),border="red")
  
  # plot(COUNTIESList[["F"]][COUNTIESList[["F"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("blue",0.4))#,border=grey(0.5))
  # par(mar=c(4.5,1,1,1),xaxs="i",yaxs="i")
  # Cplot(COUNTIESplot,border=grey(0.5))
  if(c%in% 1){
    text(253761.2 ,6775270 ,labels=COUNTRIESsimp[c, ]$idunique, col="white")
  }else{
    text(504651.6,6928887 ,labels=COUNTRIESsimp[c, ]$idunique, col="white")
    
  }
  
  
}


dev.off()





## ------    5.1 BARS Other LAST YEAR ------
#layout
pdf(file=file.path(WDFigures, paste("DetectionProbBarsOtherLastYear.pdf",sep="")),
    width=10,
    height=7)
COUNTRIESsimp <- st_simplify(COUNTRIES,dTolerance = 1500,preserveTopology = T)
COUNTRIESsimp$idunique <- COUNTRIES$ISO
#plot
mx <- matrix(c(1,2),1,2,byrow = T)
nf <- layout(mx,widths=c(1,0.8),
             heights=c(1))
par(mar=c(4,4,1,1), tck=0)
#plot
plot(10, xlim = c(0, length(COUNTRIESsimp)), ylim = c(0,0.01), type ="n", xaxt="n",
     xlab = "Country", ylab = "Detection probability")
axis(2,tck=-0.02)
axis(1, c(1:length(COUNTRIESsimp$idunique)),
     labels = COUNTRIESsimp$ISO,
     cex.axis=1)#,padj = -2)
abline(v=1:(length(COUNTRIESsimp$idunique)-1)+0.5,lty=2,col=grey(0.5))
myDev <- c(-0.2,+0.2)#seq(-0.3,0.3,length.out = 2)# 

count <- 0
for(c in 1:length(COUNTRIESsimp$idunique)){
  count <- count+1
  
  
  # plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.1), type ="n", xaxt="n", xlab = "Years",
  #      ylab = "p0")#, main=paste(COUNTIES[COUNTIES$id==detCounties.original[c], ]$NAME_1, collapse = " "))
  # axis(1, at = 1:(nYears) , labels = years[1:(nYears)])
  for(s in 1:2){
    myResults <- myResultsList[[s]]
    t <- nYears
    #for(t in 1:nYears){
    # coun <- which( detCounties.original1[[t]]  %in%  CountyIndex[c])
    if(length(c)>0){
      
      
      tmp <- myResults$sims.list$p0Oth[ , c, t]
      quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
      polygon(x = c(count+ myDev[s] - widthPolygon, count+ myDev[s] + widthPolygon,
                    count+ myDev[s] + widthPolygon, count+ myDev[s] - widthPolygon ),
              y = c(quantile95[1], quantile95[1],
                    quantile95[2], quantile95[2]), 
              col=adjustcolor(myCol[s], 0.6),
              border= adjustcolor(myCol[s], 0.7))
      
      
    }
  }
  
  if(c==2){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.001)
    for(i in 1:2){
      #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
      points(0.6, 0.009-yoffset[i], pch=15,cex=1.8,col=adjustcolor(cols[i],violin.alpha95+0.1))
      #points(2, 0.035-yoffset[i], pch=15,cex=0.5,col=adjustcolor(pt.col[i],violin.alpha50))
      text(0.7, 0.009-yoffset[i], labels[i],cex=1.3,pos=4)
    }
  }
}

## PLOT REGION MAP    
par(mar=c(0,0,0,0))
aggCounty <- nngeo::st_remove_holes(COUNTRIESsimp)
plot(st_geometry(aggCounty),border="white",col=grey(0.5),lwd=0.1)

for(c in 1:2){
  
  aggCounty <- COUNTRIESsimp[c, ]
  
  
  # plot(COUNTIESList[["M"]][COUNTIESList[["M"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("red",0.4),border="red")
  
  # plot(COUNTIESList[["F"]][COUNTIESList[["F"]]$id==CountyIndex[c], ],
  #      add=T, col=adjustcolor("blue",0.4))#,border=grey(0.5))
  # par(mar=c(4.5,1,1,1),xaxs="i",yaxs="i")
  # Cplot(COUNTIESplot,border=grey(0.5))
  if(c%in% 1){
    text(253761.2 ,6775270 ,labels=COUNTRIESsimp[c, ]$idunique, col="white")
  }else{
    text(504651.6,6928887 ,labels=COUNTRIESsimp[c, ]$idunique, col="white")
    
  }
}

#}


dev.off()
## ------    5.2 MAPS  ------
# t=7
# #plot(myHabitat.listM$habitat.r)
# # pairs<- 2
# # alreadydetected <- 1
# 
# pdf(file=file.path(WDFigures, paste("MapDetectionProb.pdf",sep="")),width=6,height=6)
# for(t in 1:nYears){
#   
#   p0.r <- myHabitat.list$habitat.r
#   P0 <- ilogit(logit(myResults_M$mean$p0[nimDataM$detCountries,t])+
#                  myResults_M$mean$betaResponse*1 +
#                  myResults_M$mean$betaCovs[1]*nimDataM$detCovs[,t,1]+
#                  myResults_M$mean$betaCovs[2]*nimDataM$detCovs[,t,2]+
#                  myResults_M$mean$betaCovs[3]*nimDataM$detCovs[,t,3]
#                
#   )
#   #    
#   cells <- cellFromXY(p0.r,coordinates(myDetectorsM$main.detector.sp))
#   p0.r[cells] <- P0
#   p0.r <- mask(p0.r,myHabitat.list$habitat.poly)
#   # plot(p0.r)
#   # # p0.r <- focal(p0.r, w=matrix(1/25,nrow=5,ncol=5),na.rm=T) 
#   # plot(p0.r,main="p0")
#   # plot(myStudyArea.poly,add=T)
#   
#   # plot(p0.r)
#   logp0.r <- log(p0.r)
#   
#   plot(logp0.r,legend.args=list(text='log(p0)', side=4, font=2, line=2.5, cex=0.8),main=years[t])
# }
# dev.off()
# # nimDataF$detTracks
# 
## ------    5.3 TABLE  ------
CountiesID <- 1:dim(myResultsList[[1]]$sims.list$p0)[2]
state <- c( "Others","Scent-marking adult")
sex <- c("M", "F")
Tablep0 <- matrix(NA, nrow=length(CountiesID)+1, ncol=(nYears)*2)
rownames(Tablep0) <- c("",unlist(lapply(as.list(CountiesID),function(x) rep(x,1))))
colnames(Tablep0) <- c(unlist(lapply(as.list(years),function(x) rep(paste(x,collapse =  "-"),2))))
Tablep0[1,] <- c( rep(sex,(nYears)) )



n.digits = 3
rownamesTablep0 <- ""
for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResultsList[[2]]}else{results <- myResultsList[[1]]}
  for(i in 1:length(CountiesID)){
    rows <- which(rownames(Tablep0)==CountiesID[i])
    col <- which(Tablep0[1,]==sex[s])
    # state 1 
    Tablep0[rows[1],col] <-  paste( apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(median(x),n.digits), nsmall = n.digits)), #median 
                                    " (", apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(quantile(x,probs=0.025),n.digits), nsmall = n.digits)),#UpperCI
                                    "-" , apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(quantile(x,probs=0.975),n.digits), nsmall = n.digits)),#Lower CI
                                    ")", sep="") 
    
    rownamesTablep0[rows[1]] <- COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==i, ]$Name
    
  }
}

rownames(Tablep0) <- rownamesTablep0
Tablep0 <- Tablep0[order(row.names(Tablep0)),]

# rownames(Tablep0) <- c("",unlist(lapply(as.list(c("NO1","SE1","SE2","NO2","SE3","SE4")),function(x) rep(x,2))))
# WRITE THE FILE

write.csv(Tablep0, file=file.path(WDFigures, paste("Tablep0.csv",sep="")))



## ------  6. BETA P0S  ------
## ------    6.1 P0STRUCTURED  ------
## ------      6.1.1 BETATRACKS  ------
pdf(file=file.path(WDFigures, paste("Betap0StructuredTracks.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaCovs[,1,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovs[,1,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()


## ------      6.1.2 BETA SNOW  ------
pdf(file=file.path(WDFigures, paste("Betap0StructuredSnow.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaCovs[,2,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovs[,2,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()



## ------      6.1.2 BETA RESPONSE  ------
pdf(file=file.path(WDFigures, paste("Betap0StructuredResponse.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaResponse[,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaResponse[,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()




## ------    6.2 P0OTHER  ------
## ------      6.2.1 BETASNOW ------
pdf(file=file.path(WDFigures, paste("Betap0OtherSnow.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,1,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,1,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()


## ------      6.2.2 BETA ROADS  ------
pdf(file=file.path(WDFigures, paste("Betap0OtherRoads.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,2,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,2,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()



## ------      6.2.3 BETA SKANDOBS  ------
pdf(file=file.path(WDFigures, paste("Betap0OtherSkandobs.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,3,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,3,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()




## ------      6.2.3 BETA RESPONSE  ------
pdf(file=file.path(WDFigures, paste("Betap0OtherResponse.pdf",sep="")),width=8,height=4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- myResultsList[[s]]
  
  for(t in 1:(nYears)){
    
    quantile95 <- quantile(myResults$sims.list$betaResponseOth[,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaResponseOth[,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
    
    
    
    
  }
  t
}#i


#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}


dev.off()




## ------  6. TABLE OTHERS  ------
## ------    6.1. TABLE DENSITY AND MOVEMENT  OPSCR  ------
parameters <- c("betaDens","sigma")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{dens}^*$","$\\sigma$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")


n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableDensityMovementSCR <- matrix(NA, nrow=length(parameters)+1, ncol=(nYears)*2)
rownames(TableDensityMovementSCR) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableDensityMovementSCR) <-  unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
#c(unlist(lapply(YEARS[1:(length(YEARS))],function(x) rep(paste(x+1,collapse =  "-"),2))))
TableDensityMovementSCR[1,] <- c(rep(sex, (nYears)) )

#ST <- 2:1 ## state 3 had the index 1. state 2 has the index 2. so need to inverse it. 
for(s in 1:2){
  # choose the sex
  
  if(s==1){results <- myResultsList[[2]]
  }else{results <- myResultsList[[1]]}
  
  ### BETA DENSITY 
  # REPEAT ESTIMATES FOR ALL YEARS FOR SUCH PARAMETERS
  param <- "betaDens"
  rows <- which(rownames(TableDensityMovementSCR)==param)[1]
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste( format(round(median(results$sims.list[[param]]), 2), nsmall = 2), #median
                                               " (", format(round(quantile(results$sims.list[[param]], probs=0.025),2), nsmall = 2),#UpperCI
                                               "-" , format(round(quantile(results$sims.list[[param]], probs=0.975),2), nsmall = 2),#Lower CI
                                               ")", sep="")
  
  ### SIGMA
  # for(st in 1:2){
  param <- "sigma"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste( format(round(apply(results$sims.list[[param]] * myHabitat.list$resolution /1000,2,median), 2), nsmall = 2), #median
                                               " (", format( round(apply(results$sims.list[[param]]* myHabitat.list$resolution/1000,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                               "-" , format(round(apply(results$sims.list[[param]]* myHabitat.list$resolution/1000,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                               ")", sep="")
  
  # }
  
  ### LAMBDA
  param <- "dmean"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste( format(round(median(results$sims.list[[param]] * myHabitat.list$resolution /1000), 2), nsmall = 2), #median
                                               " (", format( round(quantile(results$sims.list[[param]]* myHabitat.list$resolution/1000, probs=0.025 ),2), nsmall = 2),#UpperCI
                                               "-" , format(round(quantile(results$sims.list[[param]]* myHabitat.list$resolution/1000, probs=0.975 ),2), nsmall = 2),#Lower CI
                                               ")", sep="")
  
}


##WRITE TABLES 
write.csv(TableDensityMovementSCR, file=file.path(WDFigures, paste("TableDensityMovement.csv",sep="")))
#write latex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command3 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR)[19:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

colnames(TableDensityMovementSCR) <- rep("", ncol(TableDensityMovementSCR))

rownames(TableDensityMovementSCR)[2:nrow(TableDensityMovementSCR)] <- parameters1



print(xtable(TableDensityMovementSCR, type = "latex",align = paste(rep("c",ncol(TableDensityMovementSCR)+1),collapse = "")),
      # scalebox=.7,
      floating = FALSE,
      add.to.row=addtorow,include.colnames=F,sanitize.text.function=function(x){x},
      file = file.path(WDTables,paste("TableDensityMovement.tex", sep="")))


##SPLIT THE TABLE IN TWO 
TableDensityMovementSCR1 <- TableDensityMovementSCR[,c(1:10)]
TableDensityMovementSCR2 <- TableDensityMovementSCR[,c(11:20)]

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
#SAVE TABLE 1

addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableDensityMovementSCR1, type = "latex",
             align = paste(rep("c", ncol(TableDensityMovementSCR1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableDensityMovement1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableDensityMovementSCR2, type = "latex",
             align = paste(rep("c", ncol(TableDensityMovementSCR2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableDensityMovement2.tex", sep="")))


##last year##

addtorow1$command <- command3
TableDensityMovementSCR3 <- TableDensityMovementSCR[,c(19:20)]
print(xtable(TableDensityMovementSCR3 , type = "latex",
             align = paste(rep("c", ncol(TableDensityMovementSCR3)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableDensityMovementLastYear.tex", sep="")))


## ------    6.2. TABLE  PARAMETERS  STRUCTURED ------
parameters <- c("betaResponse","betaCovs","betaCovs")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{1_Structured}$","$\\beta_{2_Structured}$","$\\beta_{3_Structured}$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableStructured <- matrix(NA, nrow=length(parameters1)+1, ncol=(nYears)*2)
rownames(TableStructured) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableStructured) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
TableStructured[1,] <- c(rep(sex, (nYears)) )



for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResultsList[[2]]
  }else{results <- myResultsList[[1]]}
  
  
  
  ### betaResponse
  param <- "betaResponse"
  rows <- which(rownames(TableStructured)==param)
  col <- which(TableStructured[1,]==sex[s])
  
  TableStructured[rows,col] <-  paste( format(round(apply(results$sims.list[[param]] ,2,median), 2), nsmall = 2), #median
                                       " (", format( round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                       "-" , format(round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                       ")", sep="")
  
  
  
  ### trapBetas
  for(st in 1:2){
    param <- "betaCovs"
    rows <- which(rownames(TableStructured)==param)[st]
    col <- which(TableStructured[1,]==sex[s])
    
    TableStructured[rows,col] <-  paste(format(round(apply(results$sims.list[[param]][,st,] ,2,median), 2), nsmall = 2), #median
                                        " (", format( round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                        "-" , format(round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                        ")", sep="")
  }
}

##WRITE TABLES 
write.csv(TableStructured, file=file.path(WDFigures, paste("TableStructured.csv",sep="")))
#write latex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command3 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured)[19:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

colnames(TableStructured) <- rep("", ncol(TableStructured))

rownames(TableStructured)[2:nrow(TableStructured)] <- parameters1



print(xtable(TableStructured, type = "latex",align = paste(rep("c",ncol(TableStructured)+1),collapse = "")),
      # scalebox=.7,
      floating = FALSE,
      add.to.row=addtorow,include.colnames=F,sanitize.text.function=function(x){x},
      file = file.path(WDTables,paste("TableStructured.tex", sep="")))


##SPLIT THE TABLE IN TWO 
TableStructured1 <- TableStructured[,c(1:10)]
TableStructured2 <- TableStructured[,c(11:20)]

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
#SAVE TABLE 1

addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableStructured1, type = "latex",
             align = paste(rep("c", ncol(TableStructured1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableStructured1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableStructured2, type = "latex",
             align = paste(rep("c", ncol(TableStructured2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableStructured2.tex", sep="")))

##last year##
addtorow1$command <- command3
TableStructured3 <- TableStructured[,c(19:20)]
print(xtable(TableStructured3 , type = "latex",
             align = paste(rep("c", ncol(TableStructured3)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableStructuredLastYear.tex", sep="")))




## ------    6.3. TABLE  PARAMETERS  OTHERS ------
parameters <- c("betaResponseOth","betaCovsOth","betaCovsOth","betaCovsOth")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{1_Unstructured}$","$\\beta_{2_Unstructured}$","$\\beta_{3_Unstructured}$","$\\beta_{4_Unstructured}$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableOther <- matrix(NA, nrow=length(parameters1)+1, ncol=(nYears)*2)
rownames(TableOther) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableOther) <- unlist(lapply(YEARS,function(x) c(x[2],x[2])))# unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
TableOther[1,] <- c(rep(sex, (nYears)) )



for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResultsList[[2]]
  }else{results <- myResultsList[[1]]}
  
  
  
  ### betaResponse
  param <- "betaResponseOth"
  rows <- which(rownames(TableOther)==param)
  col <- which(TableOther[1,]==sex[s])
  
  TableOther[rows,col] <-  paste( format(round(apply(results$sims.list[[param]] ,2,median), 2), nsmall = 2), #median
                                  " (", format( round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                  "-" , format(round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                  ")", sep="")
  
  
  
  ### trapBetas
  for(st in 1:3){
    param <- "betaCovsOth"
    rows <- which(rownames(TableOther)==param)[st]
    col <- which(TableOther[1,]==sex[s])
    
    TableOther[rows,col] <-  paste(format(round(apply(results$sims.list[[param]][,st,] ,2,median), 2), nsmall = 2), #median
                                   " (", format( round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                   "-" , format(round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                   ")", sep="")
  }
}

##WRITE TABLES 
write.csv(TableOther, file=file.path(WDFigures, paste("TableOther.csv",sep="")))
#write latex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command3 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther)[19:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(TableOther) <- rep("", ncol(TableOther))

rownames(TableOther)[2:nrow(TableOther)] <- parameters1



print(xtable(TableOther, type = "latex",align = paste(rep("c",ncol(TableOther)+1),collapse = "")),
      # scalebox=.7,
      floating = FALSE,
      add.to.row=addtorow,include.colnames=F,sanitize.text.function=function(x){x},
      file = file.path(WDTables,paste("TableOther.tex", sep="")))


##SPLIT THE TABLE IN TWO 
TableOther1 <- TableOther[,c(1:10)]
TableOther2 <- TableOther[,c(11:20)]

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
#SAVE TABLE 1

addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableOther1, type = "latex",
             align = paste(rep("c", ncol(TableOther1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableOther1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableOther2, type = "latex",
             align = paste(rep("c", ncol(TableOther2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableOther2.tex", sep="")))



##last year##
TableOther3 <- TableOther[,c(19:20)]

addtorow1$command <- command3
TableOther3 <- TableOther[,c(19:20)]
print(xtable(TableOther3 , type = "latex",
             align = paste(rep("c", ncol(TableOther3)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableOtherLastYear.tex", sep="")))



## ------  7. TRANSITION SURFACES  ------
## ------    7.1 SET UP HABITAT  ------
habbR <- raster::disaggregate(myHabitat.list$habitat.r, fact=4)

##
#COUNTRIES <- aggregate(x = COMMUNESM, by = "ISO")
COUNTRIES <- COMMUNES %>%    group_by(ISO) %>%summarize()

COUNTRIESsimp <- COUNTRIES#gSimplify(COUNTIES,tol = 100, topologyPreserve = T)
COUNTRIESsimp$NAME_1 <- COUNTRIES$ISO
plot(st_geometry(myHabitat.list$habitat.poly))
plot(st_geometry(COUNTRIESsimp), add=T)

## IDENTIFY COUNTIES
myStudyArea.polyMnoHoles <- nngeo::st_remove_holes(myHabitat.list$habitat.poly)
habbRCountieswbuff <- habbR <- mask(habbR, myStudyArea.polyMnoHoles)
habbR[habbR==0] <- NA
plot(habbR)
# identify SWEDEN/NORWAY in the raster
SWE <- COUNTRIESsimp[which(COUNTRIESsimp$ISO %in% c("SWE")),]     ## Just take Sweden
NOR <- COUNTRIESsimp[which(COUNTRIESsimp$ISO %in% c("NOR")),]     ## Just take Norway

this.r <- fasterize(SWE,habbR )
habbR[this.r==1] <- 2
this.r <- fasterize(NOR,habbR )
habbR[this.r==1] <- 1

# # NOR <- gSimplify(spgeom = NOR, tol = 500)
# # SWE <- gSimplify(spgeom = SWE, tol = 500)
# #this.r <- RasterizePolygon(poly=SWE, r=habbR,CoverToKeepHabitat=50,fasterize=TRUE)
# habbR[this.r==1] <- 2
# this.r <- RasterizePolygon(poly=NOR, r=habbR,CoverToKeepHabitat=50,fasterize=TRUE)
# habbR[this.r==1] <- 1

plot(habbR)
# plot(COUNTRIES,add=T)
##remove buffer from habitat
# e <- extent(myHabitat.listM$habitat.poly)
# e.sp <- as(e, 'SpatialPolygons')
##GET THE BUFFERLESS AREA
# habitat.r <- myHabitat.listM$habitat.r
# habitat.r[habitat.r!=1] <- NA
plot(habbR)
habbR <- mask(habbR, myStudyArea.polyMnoHoles)

plot(habbR)

# convert to text
habbR[habbR==2] <- "SWE"
habbR[habbR==1] <- "NOR"
habbR[habbR==0] <- NA
gc()

## AC DENSITY BASED
habIDCells.mx <- habbR
habIDCells.mx[] <- 1:ncell(habbR)
habIDCells.mx <- as.matrix(habIDCells.mx)


# SET UP A MATRIX ROW (NUMBER OF REGIONS), COLUMN NUMBER OF CELLS
# FILL IN WITH 1 AND 0 TO ASSIGN EACH CELL TO A REGION
# CELL THAT ARE NOT HABITAT OR WITHIN BUFFER ASSIGNED TO 0
regionID <- habbR[]
regionIDunique <- unique(regionID)
regionIDunique <- regionIDunique[!is.na(regionIDunique)]
regionIDmat <- do.call(rbind,lapply(regionIDunique,function(x)habbR[]== x  ))
regionIDmat[is.na(regionIDmat)] <- 0
row.names(regionIDmat) <- unique(regionIDunique)

#sourceCpp("C:/Personal_Cloud/OneDrive/Work/Coding/Rcpp/GetTransitionSurface1.cpp")
##sourceCpp("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/cpp/GetTransitionSurface1.cpp")
habbRxy <- coordinates(habbR)  
colnames(habbRxy) <- c("x","y")
myResultsSXYZ_MF$sims.list$scaledsxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ_MF$sims.list$sxy,
                                                                 coordsHabitatGridCenter = habbRxy)$coordsDataScaled

## ------    7.2 TRANSITION PROBABILITY OTHER CAUSES  ------
TransitionSurfaceOther <- list()
iter <- sample(1:dim(myResultsSXYZ_MF$sims.list$scaledsxy)[1], size = 100)#dim(densityInputCountries$sx)[1])



#for(t in 6:(nYears-1)){
for(t in 1:(nYears-1)){
  TransitionSurfaceOther[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$scaledsxy[,IDFemales,1,t],
                                                       myResultsSXYZ_MF$sims.list$scaledsxy[,IDFemales,2,t],
                                                       myResultsSXYZ_MF$sims.list$z[,IDFemales,t],
                                                       myResultsSXYZ_MF$sims.list$z[,IDFemales,t+1],
                                                       habIDCells.mx,
                                                       regionID = regionIDmat,
                                                       stateFrom = c(2),
                                                       stateTo = c(3),
                                                       ncell = ncell(habbR),
                                                       probs=c(0.025,0.975),
                                                       returnPosteriors = F)
}

TransitionSurfaceOther[[t]]$PosteriorTransitionRegion
TransitionSurfaceOther[[2]]$SummaryTransitionRegion


### TRY A SPATIAL PLOT ##


SpatialRaster <- habbR
SpatialRaster[] <- TransitionSurfaceOther[[t]]$MeanCell
plot(SpatialRaster)
plot(myHabitat.list$habitat.poly,add=T)



## ------    7.2 TRANSITION PROBABILITY CULLING  ------
TransitionSurfaceCulling <- list()
for(t in 1:(nYears-1)){
  TransitionSurfaceCulling[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$scaledsxy[,,1,t],
                                                         myResultsSXYZ_MF$sims.list$scaledsxy[,,2,t],
                                                         myResultsSXYZ_MF$sims.list$z[,,t],
                                                         myResultsSXYZ_MF$sims.list$z[,,t+1],
                                                         habIDCells.mx,
                                                         regionID = regionIDmat,
                                                         stateFrom = c(2),
                                                         stateTo = c(3),
                                                         ncell = ncell(habbR),
                                                         probs=c(0.025,0.975),
                                                         returnPosteriors = F)
}

TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion
TransitionSurfaceCulling[[t]]$SummaryTransitionRegion

### TRY A SPATIAL PLOT ##


SpatialRaster <- habbR
SpatialRaster[] <- TransitionSurfaceCulling[[t]]$MeanCell
plot(SpatialRaster)
plot(myHabitat.list$habitat.poly,add=T)


## ------    7.3 PLOT  ------
pdf(file=file.path(WDFigures, paste("CountryMortalityRates.pdf",sep="")), width = 9, height = 7)
par(mfrow=c(1,2))
offset <- c(-0.2,0.2)
plot(-10, xlim=c(0,nYears),ylim=c(0,1), ylab="Mortality rate" ,xaxt="n")
axis(1,at=c(1:nYears),labels = years)
col <- c("red","blue")
for(t in 1:(nYears-1)){
  for(r in 1:2){
    plot.violins(list(TransitionSurfaceOther[[t]]$PosteriorTransitionRegion[r,]),
                 at = t +offset[r],
                 x=1, col=col[r], alpha = 0.5,add=T)
  }
}

offset <- c(-0.2,0.2)
plot(-10, xlim=c(0,nYears),ylim=c(0,1), ylab="Mortality rate culling" ,xaxt="n")
axis(1,at=c(1:nYears),labels = years)
col <- c("red","blue")
for(t in 1:(nYears-1)){
  for(r in 1:2){
    plot.violins(list(TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion[r,]),
                 at = t +offset[r],
                 x=1, col=col[r], alpha = 0.5,add=T)
  }
}

legend("topright", fill=col, legend=c("Norway","Sweden"))

dev.off()

# 
# 


## ------  8. OTHER PLOTS  ------
## ------    8.1 SKANDOBS ------
habitat.detectors <- aggregate(rasterToPolygons(disaggregate(myHabitat.list$habitat.rWthBuffer, 
                                                             fact=2),fun=function(x)x==1))
skandobs.r1 <-skandobs.r <- disaggregate(myHabitat.list$habitat.rWthBuffer, 
                                         fact=2)

pdf(file= file.path(WDFigures, paste("SkandobsRovbaseCovariates.pdf", sep="")), width = 12, height = 8)
par(mfrow=c(2,5),mar=c(1,1,3,1))
rrr <- list()
for(t in 1:nYears){
  plot(habitat.detectors,border=grey(0.8),col=grey(0.8))
  skandobs.r1[skandobs.r[]%in%1]<- nimData$detCovsOth[,t,3]
  rrr[[t]] <- skandobs.r1
  pol <- rasterToPolygons(skandobs.r1,fun = function(x) x==1)
  plot(pol,col="darkgreen",border=NA,add=T)
  mtext(paste(YEARS[[t]][2]))
}
dev.off()

save( rrr,file=file.path(WDFigures, paste("Skandobs.RData", sep="")))


## ------    7.4 SAVE  ------

# TransitionPosteriorCulling <- 
#    TransitionPosteriorOther <- list()
# for(t in 1:(nYears-1)){
#    TransitionPosteriorCulling[[t]]<- TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion
#    TransitionPosteriorOther[[t]]<- TransitionSurfaceOther[[t]]$PosteriorTransitionRegion
# }
# 
# save(TransitionPosteriorCulling,TransitionPosteriorOther,
#      file = file.path(paste(dir.dropbox,"/wolverine/CM/2021/plot24/Figure/",sep=""), "CountryTransitionAliveCulledOther.RData" ))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###c++
# t=2
# habbRtransWolverine <-habbRid<-  disaggregate(myHabitat.listM$habitat.r,fact=4)
# habbRid[] <-1:ncell(habbRtransWolverine) 
# myResultsSXYZ_MF$sims.list$sxy <- UTMToGrid(data.sxy = myResultsSXYZ_MF$sims.list$sxy,
#                                             grid.sp = SpatialPoints(coordinates(habbRtransWolverine)) )$data.scaled.xy
# 
# 
# rr <- habbRtransWolverine
# rr[]<- 1:ncell(habbRtransWolverine)
# rr[habbRtransWolverine[]==0] <- 0
# rr[habbRtransWolverine[]==1] <- 1:sum(habbRtransWolverine[]==1)
# image(as.matrix(rr))
# t=1
# 
# ##COMPUTE TRANSITION SURFACES
# TransitionSurfaceWolverine <- list()
# for(t in 1:(nYears-1)){
#    TransitionSurfaceWolverine[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$sxy[,,1,t],
#                                                        myResultsSXYZ_MF$sims.list$sxy[,,2,t],
#                                                        myResultsSXYZ_MF$sims.list$z[,,t],
#                                                        myResultsSXYZ_MF$sims.list$z[,,t+1],
#                                                        as.matrix(rr),
#                                                        stateFrom = 2,
#                                                        stateTo = c(3),
#                                                        ncell = max(rr[]),
#                                                        probs=c(0.025,0.975))
#    
#    TransitionSurfaceWolverine[[t]]$PosteriorTransition <- NULL
#    gc()
#    rrr <-   habbRtransWolverine
#    rrr[habbRtransWolverine==1] <- TransitionSurfaceWolverine[[t]]$MeanCell
#    plot(rrr)
#    # if(t <nYears){
#    # plot(dead[dead$Year==years[t],],pch=16,add=T, cex=0.1, col=adjustcolor("black",alpha.f = 0.2))
#    # }
#    sum(rrr[])
# }
# 
# ##SAVE OBJECTS
# save(TransitionSurfaceWolverine,habbRtransWolverine,
#      file = file.path(paste(dir.dropbox,"/wolverine/CM/2021/plot24/Figure/",sep=""), "TransitionAliveCulled.RData" ))
# 
# ###LARGE RASTER 
# ###c++
# t=2
# habbRtransWolverine <-habbRid<-  aggregate(disaggregate(myHabitat.listM$habitat.r,fact=2),fact=3)
# habbRid[] <-1:ncell(habbRtransWolverine) 
# myResultsSXYZ_MF$sims.list$sxy <- UTMToGrid(data.sxy = myResultsSXYZ_MF$sims.list$sxy,
#                                             grid.sp = SpatialPoints(coordinates(habbRtransWolverine)) )$data.scaled.xy
# 
# 
# rr <- habbRtransWolverine
# rr[]<- 1:ncell(habbRtransWolverine)
# rr[habbRtransWolverine[]==0] <- 0
# rr[habbRtransWolverine[]>0] <- 1:sum(habbRtransWolverine[]>0)
# image(as.matrix(rr))
# t=1
# 
# ##COMPUTE TRANSITION SURFACES
# TransitionSurfaceWolverine <- list()
# for(t in 1:(nYears-1)){
#   TransitionSurfaceWolverine[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$sxy[,,1,t],
#                                                            myResultsSXYZ_MF$sims.list$sxy[,,2,t],
#                                                            myResultsSXYZ_MF$sims.list$z[,,t],
#                                                            myResultsSXYZ_MF$sims.list$z[,,t+1],
#                                                            as.matrix(rr),
#                                                            stateFrom = 2,
#                                                            stateTo = c(3),
#                                                            ncell = max(rr[]),
#                                                            probs=c(0.025,0.975))
#   
#   TransitionSurfaceWolverine[[t]]$PosteriorTransition <- NULL
#   gc()
#   rrr <-   habbRtransWolverine
#   rrr[habbRtransWolverine>0] <- TransitionSurfaceWolverine[[t]]$MeanCell
#   plot(rrr)
#   # if(t <nYears){
#   # plot(dead[dead$Year==years[t],],pch=16,add=T, cex=0.1, col=adjustcolor("black",alpha.f = 0.2))
#   # }
#   sum(rrr[])
# }
# 
# ##SAVE OBJECTS
# save(TransitionSurfaceWolverine,habbRtransWolverine,
#      file = file.path(paste(dir.dropbox,"/wolverine/CM/2021/plot24/Figure/",sep=""), "TransitionAliveCulled30KM.RData" ))
# 
# 
# 
# 
# 
## ------  8. TABLES OF #NGS SAMPLES, #DEAD RECOVERIES & #IDs DETECTED ------
## ------    8.1 OVERALL NUMBERS ------
load(file.path(WD, modelNameM,paste(modelNameM,"_NGSData.RData",sep="")))
load(file.path(WD, modelNameF,paste(modelNameF,"_NGSData.RData",sep="")))

## --- SOME TALLIES TO CHECK THINGS
# --- NGS
NGS <- rbind(myData.aliveF$alive, myData.aliveM$alive)
table(NGS$Year)
length(NGS)

# NGS.all <- rbind(myFullData.spF$alive,myFullData.spM$alive)
# NGS.all <- NGS.all[NGS.all$Year%in%years, ]
# table(NGS.all$Year)
# length(NGS.all)

#### FOR REPORT SUMMARY
length(NGS$Id)
length(NGS$Id[NGS$Sex=="Hunn"])
length(NGS$Id[NGS$Sex=="Hann"])

length(NGS$Id[NGS$Country=="S"])/nrow(NGS)


length(unique(NGS$Id))
length(unique(NGS$Id[NGS$Sex=="Hunn"]))
length(unique(NGS$Id[NGS$Sex=="Hann"]))


# --- DEAD RECOVERY
dead <- rbind(myFullData.spF$dead.recovery,myFullData.spM$dead.recovery)
table(dead$Year)
length(dead)



## ------    8.2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX------
NGSCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
NGSCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- length(temp[temp$Country %in% "N", ])
    NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- length(temp[temp$Country %in% "S", ])
    NGSCountrySEX["Total",ye[t] + sex1[s]] <- length(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))



print(xtable(NGSCountrySEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEX))),collapse = "")),
      #scalebox = .8,
      floating = FALSE,include.colnames=F,
      add.to.row = addtorow,
      file = file.path(WDTables,paste("NGSCountrySEX.tex",sep="")))

## ------    8.2. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------

NGSidCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
NGSidCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
    NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
    NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
    
  }#t
  
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))


write.csv(NGSidCountrySEX ,file = file.path(WDTables,paste("NGSidCountrySEX.csv",sep="")))

print(xtable(NGSidCountrySEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSidCountrySEX))),collapse = "")),
      #scalebox = .8, 
      floating = FALSE,include.colnames=F,
      add.to.row = addtorow,
      file = file.path(WDTables,paste("NGSidCountrySEX.tex",sep="")))


### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
NGSidCountryTotal <- matrix(0, ncol = nYears, nrow = 1)
row.names(NGSidCountryTotal) <- c("Total")
colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))
for(t in 1:nYears){
  temp <- NGS[NGS$Year == years[t] , ]
  NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
}#t

write.csv(NGSidCountryTotal,file = file.path(WDTables,paste("TotalIdDetected.csv",sep="")))

### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR PER SEX



## ------    8.3. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------

DeadidCountrySEX <- matrix(0, ncol = nYears*2+1, nrow = 6)
row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))
DeadidCountrySEX[1,] <- c("",rep(c("F","M"),nYears))
DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
###
MortalityNames <- unique(as.character(dead$DeathCause))
table(as.character(dead$DeathCause))
legalCauses <- MortalityNames[grep("JF", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("9", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("23", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("28", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Rifle", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("18", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("17", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Uspesifisert", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Fellefangst", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Hagle", MortalityNames)])

## SEPARATE MORTALITIES
cause <- c("other","legal culling")

for(t in 1:nYears){
  for(s in 1:2){
    for(d in 1:2){
      if(d==1){temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$DeathCause %in% legalCauses), ]
      }else{
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$DeathCause %in% legalCauses, ]
      }
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Norway" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
      
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Sweden" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
    }#t
    DeadidCountrySEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(DeadidCountrySEX[2:6,ye[t] + sex1[s]+1]))
  }
}


##summary
#Other causes
sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
#legal
sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))

sum(as.numeric(DeadidCountrySEX[c(3,5),2:ncol(DeadidCountrySEX)]))/sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))


#write latex
addtorow <- list()

addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))
# REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC


multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))
# addtorow$pos <- list(c(0),0)
# addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
#                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))

# xTableState <- xtable(TableState)

# rownames(TableState)[2:5] <- c("$\\rho$","$\\phi$","h","w")

# rownames(TableState)[2:8] <- c("$\\gamma$","$\\phi$","   ","h","  ","w","  ")


print(xtable(DeadidCountrySEX, type = "latex",
             align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("DeadidCountrySEX.tex", sep="")))



## ------ 9. plot regions ------

NewCountySwe <- readOGR(file.path(dir.dropbox,"/DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp"))
plot(NewCountySwe, border="red")

COMMUNES_NOR <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp", sep = ""))   ## Communal map of Norway
COMMUNES_SWE <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp", sep = ""))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)
## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- aggregate(x = COMMUNES, by = "NAME_1")

## only select the norwegian counties
COMMUNESNOR <- COMMUNES[COMMUNES$NAME_0=="Norway",]
NORWAY <- aggregate(x = COMMUNESNOR, by = "NAME_1")
plot(NORWAY)
NAME_1 <- as.character(NORWAY$NAME_1)
df.CountiesRegions <- matrix(c(
  "Finnmark",         8,
  "Troms",            8,
  "Nordland",         7,
  NAME_1[15],         6,
  NAME_1[9],         6,
  NAME_1[8],         6,
  "Hedmark",          5,
  "Oppland",          3,
  NAME_1[2],          4,
  "Oslo",             4,
  "Akershus",         4,
  "Sogn og Fjordane", 1,
  "Hordaland",        1,
  "Rogaland" ,        1,
  "Vest-Agder",       1,
  "Aust-Agder",       2,
  "Telemark" ,        2,
  "Buskerud" ,        2,
  "Vestfold",         2), byrow=T,ncol=2)


NORWAY$NAME_1 <- as.character(NORWAY$NAME_1) 
for(i in 1:nrow(df.CountiesRegions)){
  NORWAY$NAME_1[NORWAY$NAME_1 %in% df.CountiesRegions[i,1]] <- df.CountiesRegions[i,2]
}
NORWAY1 <- aggregate(NORWAY,by="NAME_1")
plot(NORWAY1)



# RENAME THE FIELDS SO THEY MATCH BETWEEN NORWEGIAN AND SWEDISH LAYERS
NewCountySwe <- NewCountySwe[,"LANSNAMN"]
colnames(NewCountySwe@data) <- "NAME_1"
NewCountySwe$Country <- "SWE"
NORWAY1$Country <- "NOR"

# MERGE THE 2 LAYERS.
# THERE IS SOME SPACE BETWEEN THE TWO LAYERS, BUT IT DOESNT MATTER. 
# IF A CELL IS NOT ASSIGNED TO ANY COUNTY THEN WE ASSIGN IT THE CLOSEST COUNTY BELOW. 
COUNTIESsimp <- rbind(NORWAY1,NewCountySwe)#, NewCountySwe, makeUniqueIDs = TRUE) 

country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
col <- c("firebrick2", "deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway", "Sweden")
border.col <- NA

pdf(file = file.path(WDTables, "RegionMaps.pdf"),
    width = 9, height = 13, pointsize = 12)
par(mar=c(0,0,0,0))
plot(gSimplify(COUNTIESsimp, tol=5000), col=NA, border=NA)

#NORWAY
CARNIVORE.REGIONS <- COUNTIESsimp[COUNTIESsimp$Country%in% "NOR",]
CARNIVORE.REGIONS1 <- gSimplify(CARNIVORE.REGIONS, tol=200, topologyPreserve = T)
CARNIVORE.REGIONS1$Region <- CARNIVORE.REGIONS$NAME_1
# region <- gUnaryUnion(CARNIVORE.REGIONS1, id = CARNIVORE.REGIONS1$Region)
region <- gSimplify(NORWAY1, tol=500)


#---SPECIFY COLOR PALETTE

this.col <- sequential_hcl(1+length(unique(NORWAY1$NAME_1)), "Reds 3")
this.col <- this.col[-length(this.col)]
set.seed(100)
this.col <- sample(this.col)
plot(region,add=T, col=this.col, border=border.col, lwd=1)


# SWEDEN
#swedenCounties1 <- COUNTIESsimpCarnivoreRegions[COUNTIESsimpCarnivoreRegions$Country%in% "SWE",]
NewCountySwe$NAME_1 <- as.character(NewCountySwe$NAME_1)
#swedenCounties1 <- swedenCounties1[-which(swedenCounties1$NAME_1%in%"Gotlands län"),]

NewCountySwe1 <- gSimplify(NewCountySwe,tol=200,topologyPreserve = T)
swedenCounties2 <- NewCountySwe1
swedenCounties2$NAME_1 <- as.character(NewCountySwe$NAME_1)

swedenCounties2 <- swedenCounties2[!(swedenCounties2$NAME_1%in% "Gotlands lÃ¤n"), ]

#---SPECIFY COLOR PALETTE
this.col <- sequential_hcl(25+length(unique(swedenCounties2$NAME_1)), "Blues 3")
this.col <- this.col[1:length(unique(swedenCounties2$NAME_1))]
set.seed(100)
this.col <- sample(this.col)

plot(RemoveHolesSp(swedenCounties2),add=T, col=this.col, border=border.col, lwd=1)

#---LABELS
CARNIVORE.REGIONS2 <- aggregate(CARNIVORE.REGIONS1, by="Region")
#raster::text(CARNIVORE.REGIONS2,labels = paste("Region ",CARNIVORE.REGIONS2$Region,sep=""),col=grey(0.05),cex=1.3)
raster::text(NORWAY1,labels = NORWAY1$NAME_1, cex=1.3,
             col=ifelse(NORWAY1$NAME_1 %in% c(5), grey(0.7), grey(0)))

swedenCounties2$NAME_1 <- as.character(unlist(lapply( strsplit(swedenCounties2$NAME_1, " "),
                                                      function(x) x[1])))
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="GÃ¤vleborgs"] <- "Gävleborg"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤rmlands"] <- "Värmland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–stergÃ¶tlands"] <- "Östergötland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¤mtlands"] <- "Jämtland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¶nkÃ¶pings"] <- "Jönköping" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SkÃ¥ne"] <- "Skåne" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stmanlands"] <- "Västmanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stra"] <- "Västra Götaland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sternorrlands"] <- "Västernorrland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sterbottens"] <- "Västerbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SÃ¶dermanlands"] <- "Södermanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–rebro"] <- "Örebro" 

swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Norrbottens"] <- "Norrbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Hallands"] <- "Halland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Kronobergs"] <- "Kronoberg" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Stockholms"] <- "Stockholm" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 

####
swedenCounties2$Region <- c("Mellersta",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Norra","Mellersta")
swedenCounties2Regions <- aggregate(swedenCounties2,by="Region")

plot(RemoveHolesSp(swedenCounties2Regions),add=T, border=grey(0.0), lwd=3)

polygonsLabel(swedenCounties2, swedenCounties2$NAME_1,
              method = "buffer", cex=1.1,#gridpoints = 50000,
              col=ifelse(swedenCounties2$NAME_1%in%c("Jämtland"),grey(0.7),grey(0.7)))

dev.off()










##------------------------------------------------------------------------------
## ------ V. PROCESS OPSCR OUTPUT ------

myVars <- list( 
  ## WORKING DIRECTORY & MODEL NAME
  # WD = "C:/My_documents/NIMBLE/WOLVERINE",
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
               sampling.months = list(12,1:6)),   
  ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)
  
  # DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 2000,
                    detResolution = 10000,
                    detDeadResolution = 15000),
  
  # DATA GENERATION 
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 3,
                     aug.factor = 0.8),
  
  ## OUTPUT PLOTS 
  OUTPUT = list(mapResolution = 10000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)


years <- DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

if(is.null(modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(WD, modelName))){dir.create(file.path(WD, modelName))}

## ------ I.LOAD AND SELECT DATA ------
## ------ 1. HABITAT DATA ------
## ------    1.1.LOAD RAW SHAPEFILES ------
## POLYGONS OF THE REGION
GLOBALMAP <- st_read(paste(dir.dropbox,"/DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))
#plot(st_geometry(GLOBALMAP))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP[GLOBALMAP$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- COUNTRIES %>%    group_by(ISO) %>%summarize()


## POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp", sep = ""))   ## Communal map of Norway
COMMUNES_SWE <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp", sep = ""))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)
## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- COMMUNES %>%    group_by(NAME_1) %>%summarize()

#plot(st_geometry(COUNTIES))
## AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATE <- COUNTIES
#COUNTIES_AGGREGATED <- gSimplify(COUNTIES_AGGREGATED,tol=500, topologyPreserve = TRUE)
COUNTIES_AGGREGATE$id <- 1:nrow(COUNTIES_AGGREGATE)
# ggplot(COUNTIES_AGGREGATE) +
#   geom_sf(aes(fill = id)) +
#   geom_sf_label(aes(label = id))

#[CM] adjust Counties aggregation
COUNTIES_AGGREGATE$id[c(24,3,15,9,14,38,40,21,27,37,31,26,34,5,8,12,36,13,7)] <- 3
COUNTIES_AGGREGATE$id[c(39,33,23,32,29,22,4,11,20,2,10,16,25,1)] <- 4
COUNTIES_AGGREGATE$id[c(19)] <- 1
COUNTIES_AGGREGATE$id[c(35)] <- 2
COUNTIES_AGGREGATE$id[c(17,28)] <- 5
COUNTIES_AGGREGATE$id[c(18)] <- 7
COUNTIES_AGGREGATE$id[c(30)] <- 6

COUNTIES_AGGREGATE <- COUNTIES_AGGREGATE %>% group_by(id) %>% summarize()
#COUNTIES_AGGREGATE <- aggregate(x = gBuffer(COUNTIES_AGGREGATE, width = 0, byid = T), by = "id")
COUNTIES_AGGREGATED <- st_simplify(COUNTIES_AGGREGATE,preserveTopology = T,dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id
ggplot(COUNTIES_AGGREGATED) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))

## ------    1.2.CREATE STUDY AREA POLYGON ------
## CREATE STUDY AREA POLYGON BASED ON COUNTY & COMMUNES IDs
# if(!is.null(HABITAT$countyNames)){
#    myStudyArea <- COMMUNES[COMMUNES$NAME_1 %in% HABITAT$countyNames, ]
#    myStudyArea1 <- COMMUNES[COMMUNES$ISO == "NOR" & COMMUNES$ID_2 %in% HABITAT$communeNames, ]
#    myStudyArea <- aggregate(rbind(myStudyArea, myStudyArea1))
#    ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
#    myBufferedArea <- gBuffer(spgeom = myStudyArea, width = HABITAT$habBuffer)
#    myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
# }

## CREATE STUDY AREA POLYGON BASED ON x AND y EXTENTS
# if(!is.null(HABITAT$x.extent)){
#    # myStudyArea <- crop(COUNTRIES, extent(HABITAT$x.extent, HABITAT$y.extent))
#    myStudyArea <- st_crop(COUNTRIES, st_bbox(extent(HABITAT$x.extent, HABITAT$y.extent)))
#    ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
#    myBufferedArea <- gBuffer(spgeom = myStudyArea, width = HABITAT$habBuffer)
#    myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
#    myBufferedArea <- crop(myBufferedArea, extent(HABITAT$x.extent, HABITAT$y.extent))
# }

## CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
if(!is.null(HABITAT$countries)){
  myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% HABITAT$countries, ]
  # myStudyArea$id <- 1
  # myStudyArea <- myStudyArea %>% group_by(id) %>% summarize()
  # 
  ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
  myBufferedArea <- st_buffer(st_as_sf(myStudyArea) ,dist = HABITAT$habBuffer)
  myBufferedArea$id <- 1
  myBufferedArea <- myBufferedArea %>% group_by(id) %>% summarize()
  myBufferedArea <- st_intersection(myBufferedArea, GLOBALMAP)
  
}

## CREATE STUDY AREA POLYGON BASED ON POLYGON COORDINATES
# if(!is.null(HABITAT$coord.x)){
#    myStudyPoly <- MakePolygon( coord.x = HABITAT$coord.x,
#                                coord.y = HABITAT$coord.y)
#    myStudyArea <- gIntersection(myStudyPoly, COUNTRIES)
#    ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
#    myBufferedArea <- gBuffer(spgeom = myStudyArea, width = HABITAT$habBuffer)
#    myBufferedArea <- gIntersection(myBufferedArea, GLOBALMAP)
#    myBufferedArea <- gIntersection(myBufferedArea, myStudyPoly)
# }

##-- Plot check
if(plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myBufferedArea), add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
  plot(st_geometry(myStudyArea), add = TRUE, col ="red")
}

## ------ 2. LOAD NECESSARY OBJECTS ------
# LOAD OBJECTS
load(file.path(WD, modelName, "NecessaryObjects.RData" ))
#load the habitat 
load(paste(dir.dropbox,"/DATA/GISData/spatialDomain/Habitat20kmNewSweCounties.RData",sep=""))
load(paste(dir.dropbox,"/DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData",sep=""))

#load a gis layer
GLOBALMAP <- st_read(paste(dir.dropbox,"/DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))


#Name of the female and male models
years <- DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))


modelNameF <- "54.aJ_FaCleaned2024"
modelNameM <- "54.aJ_MaCleaned2024"


## LOAD NECESSARY OBJECTS
## INFILES
#FEMALES
load(file.path(WD, modelName,"Hunn",paste(modelName,"Hunn","_Chain","1.RData",sep="")))
# load(file.path("C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2022", 
#                modelNameF,paste(modelNameF,"1.RData",sep="")))

nimDataF <- nimData
#MALES
load(file.path(WD, modelName,"Hann",paste(modelName,"Hann","_Chain","1.RData",sep="")))
# load(file.path("C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2022", 
#                modelNameM,paste(modelNameM,"1.RData",sep="")))

nimDataM <- nimData
# 
# 
# #SET DIRECTORY WHERE WOLF FIGURES WILL BE STORED
WDFigures <- file.path(WD, modelName,"Figure")
WDTables <- file.path(WD, modelName,"Table")

if(!dir.exists(file.path(WDFigures))){dir.create(WDFigures)}
if(!dir.exists(file.path(WDTables))){dir.create(WDTables)}



## ------  8. TABLES OF #NGS SAMPLES, #DEAD RECOVERIES & #IDs DETECTED ------
## ------    8.1 OVERALL NUMBERS ------
load(file.path(WD, modelName, "Hunn",paste(modelName,"_NGSData.RData", sep="")))
load(file.path(WD, modelName, "Hann",paste(modelName,"_NGSData.RData", sep="")))

# load(file.path("C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2022", modelNameM, paste(modelNameM,"_NGSData.RData", sep="")))


## --- SOME TALLIES TO CHECK THINGS
# --- NGS
NGS <- rbind(myData.aliveF$myData.sp, myData.aliveM$myData.sp)
table(NGS$Year)
nrow(NGS)


NGSStructured <- rbind(myData.aliveStrucF$myData.sp, myData.aliveStrucM$myData.sp)
NGSOther <- rbind(myData.aliveOthersF$myData.sp, myData.aliveOthersM$myData.sp)
table(NGSStructured$Year)
nrow(NGSStructured)+
  nrow(NGSOther)


# NGS.all <- rbind(myFullData.spF$alive,myFullData.spM$alive)
# NGS.all <- NGS.all[NGS.all$Year%in%years, ]
# table(NGS.all$Year)
# length(NGS.all)

#### FOR REPORT SUMMARY
length(NGS$Id)
length(NGS$Id[NGS$Sex=="Hunn"])
length(NGS$Id[NGS$Sex=="Hann"])

length(NGS$Id[NGS$Country=="S"])/nrow(NGS)


length(unique(NGS$Id))
length(unique(NGS$Id[NGS$Sex=="Hunn"]))
length(unique(NGS$Id[NGS$Sex=="Hann"]))



#last year
length(NGS$Id[NGS$Year %in%tail(years, n=1)])
length(NGS$Id[NGS$Sex=="Hunn" & NGS$Year %in%tail(years, n=1)])
length(NGS$Id[NGS$Sex=="Hann"& NGS$Year %in%tail(years, n=1)])


###
length(NGSStructured$Id)
length(NGSStructured$Id[NGSStructured$Sex=="Hunn"])
length(NGSStructured$Id[NGSStructured$Sex=="Hann"])

length(NGS$Id[NGS$Country=="S"])/nrow(NGS)

length(NGSOther$Id)
length(NGSOther$Id[NGSOther$Sex=="Hunn"])
length(NGSOther$Id[NGSOther$Sex=="Hann"])



# --- DEAD RECOVERY
dead <- rbind(myFullData.spDeadF, myFullData.spDeadM)
table(dead$Year)
length(dead)
length(unique(dead$Id[dead$Sex=="Hunn"]))
length(unique(dead$Id[dead$Sex=="Hann"]))


###
tmpdead <- dead[dead$Year %in% c(2018:2023),]
tmpdead <- tmpdead[!duplicated(tmpdead$DNAID),]

table(tmpdead$Year,tmpdead$Sex)
table(tmpdead$Year,tmpdead$Month)
table(tmpdead$Year)
nrow(table(tmpdead$Year))
duplicated(tmpdead$DNAID)

# mapview::mapview(tmpdead)



tmpNGS <- NGS[NGS$Year %in% c(2018:2023),]
table(tmpNGS$Year,tmpNGS$Sex)
table(tmpNGS$Year,tmpNGS$Month)
table(tmpNGS$Year)

###HENRIK CHECK WITH PUBLIC SAMPLES
idPublic <-
  c(
    'D555438',
    'D556590',
    'D553322',
    'D555239',
    'D554783',
    'D558440',
    'D555845',
    'D555412',
    'D556343',
    'D556344',
    'D556347',
    'D558235',
    'D557101')

which(NGSStructured$DNAID%in%idPublic)
# 
# # NGSStructured[which(NGSStructured$DNAID%in%idPublic),]@data
# # NGSOther[which(NGSOther$DNAID%in%idPublic),]@data
# 
# 
# idPublic1 <-
#   c(
#     'D555438')
# NGSStructured[which(NGSStructured$DNAID%in%idPublic1),]@data
# NGSOther[which(NGSOther$DNAID%in%idPublic1),]
# # 
# mapview(
# list(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"),
#         NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]),
# layer.name = c("Franconian districts", "Franconian breweries")
# )
#         
# mapview(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"))+
#   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# 
# 
# mapview(as(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),"Spatial"))+
#   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# 
# st_distance(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),
#             st_as_sf(NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]))
# 
# 
# 
# 
# plot(TRACKS_YEAR[[9]][TRACKS_YEAR[[9]]$RovbaseID %in% "T471191",]$geometry)
# plot(TRACKSSimple_sf[[9]][TRACKSSimple_sf[[9]]$RovbaseID %in% "T471191",]$geometry)
# plot(TRACKS[TRACKS$RovbaseID %in% "T471191",]$geometry,col="red",add=T)
# points(NGSStructured[which(NGSStructured$DNAID%in%idPublic),],pch=16)


## ------    8.2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX------
## ------      8.2.1 ALL------
NGSCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
# colnames(NGSCountrySEX) <- unlist(lapply(YEARS, function(x) c(x[2],x[2])))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

NGSCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- nrow(temp[temp$Country %in% "N", ])
    NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S", ])
    NGSCountrySEX["Total",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))



print(xtable(NGSCountrySEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEX))),collapse = "")),
      #scalebox = .8,
      floating = FALSE,include.colnames=F,
      add.to.row = addtorow,
      file = file.path(WDTables,paste("NGSCountrySEX.tex",sep="")))

write.csv(NGSCountrySEX ,file = file.path(WDTables,paste("NGSCountrySEX.csv",sep="")))

sum(as.numeric(NGSCountrySEX["Total",]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))

## ------      8.2.2 PER OBSERVATION PROCESS------
NGSCountrySEXoBS <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBS) <- c("",rep(c("Norway","Sweden","Total"),each=2))
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
#colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



NGSCountrySEXoBS[1,] <- c("",rep(c("F","M"),nYears))
NGSCountrySEXoBS[,1] <- c("",rep(c("Structured","Unstructured"),3))

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(2,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    ## structured
    tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
    
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[1], ye[t] + sex1[s] ] <- nrow(tempStruc[tempStruc$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S" | tempStruc$Country %in% "N" , ])
    
    ## Other
    tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
    
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[2], ye[t] + sex1[s] ] <- nrow(tempOther[tempOther$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S" | tempOther$Country %in% "N" , ])
    
    ###TOTAL 
    
    
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBS)[2:ncol(NGSCountrySEXoBS)])),
                                                              '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBS) <- rep("", ncol(NGSCountrySEXoBS))

#"\\multirow{1}{*}{}",
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
NGSCountrySEXoBS <- data.frame(cbind(multirowadd,NGSCountrySEXoBS))
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



print(xtable(NGSCountrySEXoBS, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEXoBS))),collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("NGSCountrySEXperObs.tex",sep="")))


sum(as.numeric(NGSCountrySEX["Total",]))

sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))

#PLOT CHECK 
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="red",add=T)
plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="red",add=T)

par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="#E69F00",add=T)
plot(myHabitat.list$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="#009E73",add=T)

dev.off()
# 
NGSStructured$Year1 <- NGSStructured$Year+1
NGSOther$Year1 <- NGSOther$Year+1

barplot(rbind(table(NGSStructured$Year1),table(NGSOther$Year1)),
        col=c("#E69F00","#009E73"))
legend("topleft",fill=c("#E69F00","#009E73"),legend=c("Structured","Other") )

#GIVE FILE TO HENRIK
tmp <- NGSOther[NGSOther$Year %in% c(2019,2020,2021),]
tmp

write.csv(tmp,file= file.path(WDTables,paste("Unstructured2020_2022.csv",sep="")))


## ------    8.4. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------
## ------      8.4.1 ALL ------

NGSidCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
#colnames(NGSidCountrySEX) <- c(unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#


NGSidCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    
    NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
    NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
    NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
    
  }#t
  
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))


write.csv(NGSidCountrySEX ,file = file.path(WDTables,paste("NGSidCountrySEX.csv",sep="")))

print(xtable(NGSidCountrySEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSidCountrySEX))),collapse = "")),
      #scalebox = .8, 
      floating = FALSE,include.colnames=F,
      add.to.row = addtorow,
      file = file.path(WDTables,paste("NGSidCountrySEX.tex",sep="")))


### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
NGSidCountryTotal <- matrix(0, ncol = nYears, nrow = 1)
row.names(NGSidCountryTotal) <- c("Total")
colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))
#colnames(NGSidCountryTotal) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

for(t in 1:nYears){
  temp <- NGS[NGS$Year == years[t] , ]
  NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
}#t

write.csv(NGSidCountryTotal,file = file.path(WDTables,paste("TotalIdDetected.csv",sep="")))

### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR PER SEX


## ------      8.4.2 PER OBSERVATION PROCESS------
NGSCountrySEXoBSid <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBSid) <- c("",rep(c("Norway","Sweden","Total"),each=2))
colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) )))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#
#colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#


NGSCountrySEXoBSid[1,] <- c("",rep(c("F","M"),nYears))
NGSCountrySEXoBSid[,1] <- c("",rep(c("Structured","Unstructured"),3))

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(2,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    ## structured
    tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
    
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[1], ye[t] + sex1[s] ] <- length(unique(tempStruc$Id[tempStruc$Country %in% "N" ])) 
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id[tempStruc$Country %in% "S" ]))
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id))
    
    ## Other
    tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
    
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[2], ye[t] + sex1[s] ] <- length(unique(tempOther$Id[tempOther$Country %in% "N" ])) 
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id[tempOther$Country %in% "S" ]))
    NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id))
    
    ###TOTAL 
    
    
  }#t
}


addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBSid)[2:ncol(NGSCountrySEXoBSid)])),
                                                              '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBSid) <- rep("", ncol(NGSCountrySEXoBSid))

#"\\multirow{1}{*}{}",
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
NGSCountrySEXoBSid <- data.frame(cbind(multirowadd,NGSCountrySEXoBSid))
colnames(NGSCountrySEXoBSid) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#



print(xtable(NGSCountrySEXoBSid, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEXoBSid))),collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("NGSCountrySEXperObsid.tex",sep="")))



## ------    8.5. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------
DeadidCountrySEX <- matrix(0, ncol = nYears*2+1, nrow = 6)
row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))#c("",unlist(lapply(YEARS,function(x) c(x[2],x[2]))))#
#colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))#

DeadidCountrySEX[1,] <- c("",rep(c("F","M"),nYears))
DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
###
MortalityNames <- unique(as.character(dead$DeathCause))
table(as.character(dead$DeathCause))
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
## SEPARATE MORTALITIES
cause <- c("other","legal culling")

for(t in 1:nYears){
  for(s in 1:2){
    for(d in 1:2){
      if(d==1){temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$DeathCause %in% legalCauses), ]
      }else{
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$DeathCause %in% legalCauses, ]
      }
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Norway" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
      
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[,1]=="Sweden" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
    }#t
    DeadidCountrySEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(DeadidCountrySEX[2:6,ye[t] + sex1[s]+1]))
  }
}


##summary
#Other causes
sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
#legal
sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))

sum(as.numeric(DeadidCountrySEX[c(2,3),2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))


## %of dead reco (legal) in norway
sum(as.numeric(DeadidCountrySEX[4,2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(4,5),2:ncol(DeadidCountrySEX)]))



sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="M")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,] %in% c("F","M"))]))

#write latex
addtorow <- list()

addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))
# REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC


multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))
# addtorow$pos <- list(c(0),0)
# addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
#                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# colnames(TableState) <- rep("", ncol(TableState))

# xTableState <- xtable(TableState)

# rownames(TableState)[2:5] <- c("$\\rho$","$\\phi$","h","w")

# rownames(TableState)[2:8] <- c("$\\gamma$","$\\phi$","   ","h","  ","w","  ")


print(xtable(DeadidCountrySEX, type = "latex",
             align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("DeadidCountrySEX.tex", sep="")))


#check 
# tmp <- dead[dead$Year == 2019 & 
#               !dead$DeathCause %in% legalCauses &
#               dead$Country %in% "N"& 
#               dead$Sex %in% "Hunn", ]
# length(unique(tmp$Id))
# unique(tmp$Id)
# 
# DeadidCountrySEX[,"X2022"]
# DeadidCountrySEX[,"X2022.1"]
# 
# dead[dead$Id %in% "JI416817 Ind7303 +",]$Sex
# 
# plot(COUNTRIES$geometry)
# plot(tmp$geometry,add=T,col="red",pch=16)

## ------    8.6. GET THE DETECTED INDIVIDUALS ------
n.detected <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv",sep="")))
n.detected <- n.detected[1,2:ncol(n.detected)]

## ------    8.7. SUMMARY DETECTED INDIVIDUALS PER COUNTIES ------
# myFilteredData.sp$alive$COUNTIES  <- st_intersects(myFilteredData.sp$alive[,1], COUNTIES_AGGREGATED[,1])
# myFilteredData.sp$alive$COUNTIES <- as.numeric(myFilteredData.sp$alive$COUNTIES)
# 
# 
# myFilteredData.sp$alive$COUNTIES <- apply(st_intersects(COUNTIES_AGGREGATED, myFilteredData.sp$alive, sparse = FALSE), 2, 
#       function(col) {which(col)})
# myFilteredData.sp$alive$counties1 <- 0
# for(i in 1:nrow(myFilteredData.sp$alive)){
#   if(length(myFilteredData.sp$alive$COUNTIES[[i]])>0){
#   myFilteredData.sp$alive$counties1[i] <- myFilteredData.sp$alive$COUNTIES[[i]][1]
#   }else{
#     myFilteredData.sp$alive$counties1[[i]] <- 0
#   }
# }
# 
# par(mar=c(0,0,0,0))
# plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
# text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
#      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
#      COUNTIES_AGGREGATED$id,col="red",font=2)
# 
# 
# 
# 
# 
# 
# ## NSAMPLES
# summa <- myFilteredData.sp$alive %>%
#   group_by(counties1,Year) %>%
#   summarise(n=n()) %>%
#   st_drop_geometry()
# summa <- summa[summa$counties1>0,]
# #2023
# tmp <- summa[summa$Year %in% 2023,]
# COUNTIES_AGGREGATED$nSampl2023 <- tmp$n#[1:8,"n"]
# #2022
# tmp1 <- summa[summa$Year %in% 2022,]
# COUNTIES_AGGREGATED$nSampl2022 <- tmp1$n#[1:8,"n"]
# 
# ##
# tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("id","nSampl2022","nSampl2023")])
# bar <- t(as.matrix(tmppp[c(5,7,8),2:3]))
# colnames(bar) <- c(5,7,8)
# barplo <- barplot(bar,beside=T,ylab="N samples")
# legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))
# 
# 
# ## NdetectionsPerID 
# #COUNT NUMBER IDS 
# 
# summa$n1 <- summa$NID <- 0
# dpt <- unique(unlist(summa$counties1))
# yearsss <- c(2022,2023)
# for(t in 1:length(yearsss)){
#   for(i in 1:length(dpt)){
#     tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$counties1 %in% dpt[i] & myFilteredData.sp$alive$Year %in% yearsss[t],]
#     
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$n1 <- nrow(tmp)
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$NID <-  length(unique(tmp$Id))
#   }
# }
# #
# summa$NdetPerIDDet <- summa$n1/summa$NID
# 
# #2023
# tmp <- summa[summa$Year %in% 2023,]
# COUNTIES_AGGREGATED$detPerID2023 <- tmp$NdetPerIDDet#[1:8,"n"]
# #2022
# tmp1 <- summa[summa$Year %in% 2022,]
# COUNTIES_AGGREGATED$detPerID2022 <- tmp1$NdetPerIDDet#[1:8,"n"]
# 
# 
# tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("detPerID2022","detPerID2023")])
# bar <- t(as.matrix(tmppp[c(5,7,8),]))
# colnames(bar) <- c(5,7,8)
# 
# 
# par(mfrow=c(1,2))
# par(mar=c(0,0,0,0))
# plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
# text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
#      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
#      COUNTIES_AGGREGATED$id,col="red",font=2)
# par(mar=c(4,5,1,1))
# barplo <- barplot(bar,beside=T,ylab="average Dets per IDS")
# legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))




## ------ II.GET THE MCMC ESTIMATES  -----
## ------  1.GET AND COMPILE BITES ------
# COMPILE CHARACTERISTICS 
bitesize <- 250
burnin <- 10000
NSkipBites <- burnin/bitesize
nthinsxy <- 5 #thinnumber for the sxy AND z values (necessary to save memory) 


# Retrieve the minimum number of bites per chain

outDirectoriesF <- list.files(file.path(WD, modelName,"Hunn"))[grep(paste("NimbleOutFOR", modelName,"Hunn",sep=""),
                                                                                  list.files(file.path(WD, modelName,"Hunn")))]
path.listF <- file.path(WD, modelName,"Hunn", outDirectoriesF)
numBitesF <- unlist(lapply(path.listF, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)/2
}))

outDirectoriesM <- list.files(file.path(WD, modelName,"Hann"))[grep(paste("NimbleOutFOR", modelName,"Hann",sep=""),
                                                                                  list.files(file.path(WD, modelName,"Hann")))]
path.listM <- file.path(WD, modelName,"Hann", outDirectoriesM)
numBitesM <- unlist(lapply(path.listM, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)/2
}))

minBites <- floor(min(c(numBitesF,numBitesM)))
#minBites <- 25



## ------    1.1 FEMALES ------
## GO TROUGH THE BITES AND GET THEM
nimOutput <- nimOutputSXY <- RUNTIME <- list()
#gc()
for(p in 1:length(path.listF)){
  print(path.listF[p])
  outfiles <- list.files(path.listF[p])
  out <- outSXYZ <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path(path.listF[p], paste("bite_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3]
    #params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
    # REMOVE SXY AND Z, NO THINING
    #parmIndex <- which(! params.simple %in% c("sxy","z"))
    
    if(sum(is.na(this.sample))>0){
      id <- unique(which(is.na(this.sample),arr.ind = T)[,1])
      this.sample <- this.sample[-id,]
      print(x)
      print("here")
      
    }
    
    out[[x]] <- this.sample#[,parmIndex]#[ ,parmIndex]
    
    # KEEP SXY AND Z, THINING
    #nthins <- seq(1,dim(this.sample)[1], by=nthinsxy)
    load(file.path(path.listF[p], paste("biteSxyZ_", x, ".RData", sep = "")))
    
    #parmIndex <- which(params.simple %in% c("sxy","z","N","sigma"))
    outSXYZ[[x]] <- this.sampleSxyZ#[nthins,parmIndex]#[ ,parmIndex]
    
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  out.mxSXY <- do.call(rbind, outSXYZ)
  
  nimOutput[[p]] <- as.mcmc(out.mx)
  nimOutputSXY[[p]] <- as.mcmc(out.mxSXY)
}#p

## COMPILE THE RESULTS
nimOutput <- as.mcmc.list(nimOutput)
nimOutputSXY <- as.mcmc.list(nimOutputSXY)

myResults_F <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))
myResultsSXYZ_F <- ProcessCodaOutput(nimOutputSXY,params.omit = c("sxy","z"))



## ------    1.2 MALES ------
## GO TROUGH THE BITES AND GET THEM
nimOutput <- nimOutputSXY <- RUNTIME <- list()
#gc()
for(p in 1:length(path.listM)){
  print(path.listM[p])
  outfiles <- list.files(path.listM[p])
  out <- outSXYZ <- runtime <- list()#[CM]
  for(x in NSkipBites:minBites){
    print(x)
    load(file.path(path.listM[p], paste("bite_", x, ".RData", sep = "")))
    runtime[[x]] <- RunTime[3]
    #params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
    # REMOVE SXY AND Z, NO THINING
    #parmIndex <- which(! params.simple %in% c("sxy","z"))
    
    if(sum(is.na(this.sample))>0){
      id <- unique(which(is.na(this.sample),arr.ind = T)[,1])
      this.sample <- this.sample[-id,]
      print(x)
      print("here")
      
    }
    
    out[[x]] <- this.sample#[,parmIndex]#[ ,parmIndex]
    
    # KEEP SXY AND Z, THINING
    #nthins <- seq(1,dim(this.sample)[1], by=nthinsxy)
    load(file.path(path.listM[p], paste("biteSxyZ_", x, ".RData", sep = "")))
    
    #parmIndex <- which(params.simple %in% c("sxy","z","N","sigma"))
    outSXYZ[[x]] <- this.sampleSxyZ#[nthins,parmIndex]#[ ,parmIndex]
    
  }#x
  RUNTIME[[p]] <- unlist(runtime)#[CM]
  out.mx <- do.call(rbind, out)
  out.mxSXY <- do.call(rbind, outSXYZ)
  
  nimOutput[[p]] <- as.mcmc(out.mx)
  nimOutputSXY[[p]] <- as.mcmc(out.mxSXY)
}#p

## COMPILE THE RESULTS
nimOutput <- as.mcmc.list(nimOutput)
nimOutputSXY <- as.mcmc.list(nimOutputSXY)

myResults_M <- ProcessCodaOutput(nimOutput,params.omit = c("sxy","z"))
myResultsSXYZ_M <- ProcessCodaOutput(nimOutputSXY,params.omit = c("sxy","z"))


## check how many ids are left available based on z. 
# nState <- array(NA, c(dim(myResultsSXYZ_M$sims.list$z)[2], 3, nYears))
# for(t in 1:nYears){
#   for(i in 1:dim(myResultsSXYZ_M$sims.list$z)[2]){
#     tmp <- table(myResultsSXYZ_M$sims.list$z[i,,t])
#     nState[i,as.numeric(names(tmp)),t] <- tmp 
#   }
# }
# MIN <- apply(nState,c(2,3),min)
# MAX <- apply(nState,c(2,3),max)
# 
# min(nState[,1,1])
# 
# #ndetected
# Ndete <- 0
# whichDet <- list()
# for(t in 1:nYears){
#   Ndete[t] <- sum(nimDataM$y.alive[,1,t]>0)
#   whichDet[[t]] <- which(nimDataM$y.alive[,1,t]>0)
# }
# 
#  
# whichDet
# plot(Ndete~years, ylim=c(0,1500),pch=16,ylab="N")
# points(MAX[2,]~years,col="red",pch=16)
# points(MIN[1,]~years,col="blue",pch=16)
# points(MIN[3,]~years,col="orange",pch=16)
# legend("topleft",legend=c("NDetected","NAlive Max","N Available","N state 3" ),col=c("black","red","blue","orange"),pch=16)
# 
# 
# whichDetAll <- unique(unlist(whichDet))
# 
# ## check how many ids are left available based on z. 
# nState1 <- array(NA, c(dim(myResultsSXYZ_M$sims.list$z)[2], 3, nYears))
# for(t in 1:nYears){
#   for(i in 1:dim(myResultsSXYZ_M$sims.list$z)[2]){
#     tmp <- table(myResultsSXYZ_M$sims.list$z[i,-whichDetAll,t])
#     nState1[i,as.numeric(names(tmp)),t] <- tmp 
#   }
# }
# MIN1 <- apply(nState1,c(2,3),min)
# MAX1 <- apply(nState1,c(2,3),max)
# 
# plot(Ndete~years, ylim=c(0,1500),pch=16,ylab="N")
# points(MAX1[2,]~years,col="red",pch=16)
# points(MIN1[1,]~years,col="blue",pch=16)
# points(MIN1[3,]~years,col="orange",pch=16)
# legend("topleft",legend=c("NDetected","NAlive Max","N Available","N state 3" ),col=c("black","red","blue","orange"),pch=16)
# 


# plot(myResultsSXYZ_M$sims.list$N[1:5600,1], type="l")
# points(myResultsSXYZ_M$sims.list$N[5601:(5600+5600),1], type="l",col="red")
# 
# plot(myResultsSXYZ_F$sims.list$N[1:5600,1], type="l")
# points(myResultsSXYZ_F$sims.list$N[5601:(5600+5600),1], type="l",col="red")
# 
# 
# plot(myResultsSXYZ_F$sims.list$N[1:5600,1])
## rescale coordinates
## ------    1.3 RESCALE SXY SO THEY CAN BE COMPARED ------
dimnames(myResultsSXYZ_F$sims.list$sxy)[[3]] <- c("x", "y")
myResultsSXYZ_F$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ_F$sims.list$sxy,
                                                          coordsHabitatGridCenter = myHabitat.list$habitat.xy,
                                                          scaleToGrid = FALSE)$coordsDataScaled

dimnames(myResultsSXYZ_M$sims.list$sxy)[[3]] <- c("x", "y")
myResultsSXYZ_M$sims.list$sxy <- scaleCoordsToHabitatGrid(coordsData = myResultsSXYZ_M$sims.list$sxy,
                                                          coordsHabitatGridCenter = myHabitat.list$habitat.xy,
                                                          scaleToGrid = FALSE)$coordsDataScaled


##QUICK CHECK
t=9
#Male
plot(st_geometry(myHabitat.list$buffered.habitat.poly))
points(myResultsSXYZ_M$sims.list$sxy[1,myResultsSXYZ_M$sims.list$z[1,,t] %in% 2, 2, t]~
         myResultsSXYZ_M$sims.list$sxy[1,myResultsSXYZ_M$sims.list$z[1,,t] %in% 2, 1, t], pch=16, col="red")

#Female
plot(st_geometry(myHabitat.list$buffered.habitat.poly))
points(myResultsSXYZ_F$sims.list$sxy[1,myResultsSXYZ_F$sims.list$z[1,,t] %in% 2, 2, t]~
         myResultsSXYZ_F$sims.list$sxy[1,myResultsSXYZ_F$sims.list$z[1,,t] %in% 2, 1, t], pch=16, col="red")

#
# 
# which(myResultsSXYZ_M$sims.list$sxy[1,, 2, t]> 7930000)
# id <- 690
# myResultsSXYZ_M$sims.list$z[1,id,t]
# 
# nimDataM$yDets[id,,t]
# 
# points(myDetectors$main.detector.sp[nimDataM$yDets[id,1,t],],pch=16)
# 
# tmp <- myData.aliveM$alive[ myData.aliveM$alive$Id %in% row.names(nimDataM$z[])[id] &
#                                  myData.aliveM$alive$Year %in%  years[t],]
# 
# 
# points(tmp, pch=16,col="blue")
# 
# 
# 
# points(myResultsSXYZ_F$sims.list$sxy[1,,2,1]~myResultsSXYZ_F$sims.list$sxy[1,,1,1], pch=16, col="blue")


## ------    1.4 IDENTIFY INDIVIDUALS IN THE BUFFER AND ASSIGN THEM A Z NOT ALIVE (Z==5) ------
## ------       1.4.1 CREATE POLYGON WITHOUT BUFFER ------
##MALES
# habitat.r <- myHabitat.listM$habitat.r
# habitat.r[habitat.r!=1] <- NA
# myStudyArea.polyMnoHoles <- RemoveHolesSp(myHabitat.listM$habitat.poly)
# habbRNobuffM <- mask(habitat.r, myStudyArea.polyMnoHoles)
# plot(habbRNobuffM)
# areanobuffM <- aggregate(rasterToPolygons(habbRNobuffM,fun=function(x)x==1))
# 
# plot(areanobuffM, add=T)
# plot(myHabitat.listM$habitat.rWthBuffer,add=T)

habbRNobuffM <- myHabitat.list$habitat.rWthBuffer
##FEMALES
# habitat.r <- myHabitat.listF$habitat.r
# habitat.r[habitat.r!=1] <- NA
# myStudyArea.polyFnoHoles <- RemoveHolesSp(myHabitat.listF$habitat.poly)
# habbRNobuffF <- mask(habitat.r, myStudyArea.polyFnoHoles)
# plot(habbRNobuffF)
# areanobuffF <- aggregate(rasterToPolygons(habbRNobuffF,fun=function(x)x==1))
# plot(areanobuffF, add=T)
# plot(myHabitat.listF$buffered.habitat.poly,add=T)
habbRNobuffF <- myHabitat.list$habitat.rWthBuffer

# points(myResultsSXYZ_F$sims.list$sxy[i,1,2,t]~myResultsSXYZ_F$sims.list$sxy[i,1,1,t], col="red")
# points(myResultsSXYZ_F$sims.list$sxy[i,3126,2,t]~myResultsSXYZ_F$sims.list$sxy[i,3126,1,t], col="red")

## ------       1.4.2 IDENTIFY INDIVIDUALS IN THE BUFFER AND GIVE THEM A STATE 5 ------
#MAKE A COPY 
myResultsSXYZ_F$sims.list$z1 <- myResultsSXYZ_F$sims.list$z 
myResultsSXYZ_M$sims.list$z1  <- myResultsSXYZ_M$sims.list$z 


##FEMALES
dim(myResultsSXYZ_F$sims.list$sxy)
dim( myResultsSXYZ_F$sims.list$z)
for(t in 1:nYears){
  for(i in 1:dim( myResultsSXYZ_F$sims.list$z)[1]){
    whichNA <- which(is.na(habbRNobuffF[cellFromXY(habbRNobuffF, myResultsSXYZ_F$sims.list$sxy[i,,1:2,t])]))
    myResultsSXYZ_F$sims.list$z[i,whichNA,t] <- 5
  }
  print(t)
}


###IDENTIFY INDIVIDUALS IN THE BUFFER AND GIVE THEM A STATE 5
##MALES
gc()
dim(myResultsSXYZ_M$sims.list$sxy)
dim( myResultsSXYZ_M$sims.list$z)
for(t in 1:nYears){
  for(i in 1:dim( myResultsSXYZ_M$sims.list$z)[1]){
    whichNA <- which(is.na(habbRNobuffM[cellFromXY(habbRNobuffM,myResultsSXYZ_M$sims.list$sxy[i,,1:2,t])]))
    myResultsSXYZ_M$sims.list$z[i,whichNA,t] <- 5
  }
  # gc()
  print(t)
}

gc()
## ------   1.4 COMBINE MALES AND FEMALES ------
## ASSUMING THE SAME NUMBER OF CHAINS AND ITERATIONS/ WE CAN JUST ADD N ESIMAES OF FEMALES AND MALES
myResultsSXYZ_MF <- myResultsSXYZ_M
## HERE I ONLY HAVE ONE CHAIN FOR THE MALES
dimmF <- dim(myResults_F$sims.list$N)[1]
dimmFsxy <- dim(myResultsSXYZ_MF$sims.list$sxy)[1]

myResultsSXYZ_MF$sims.list$N <- myResults_M$sims.list$NdimmF[1:dimmF[1],,] + myResults_F$sims.list$N

## ------   1.5 COMBINE SXY AND Z ------
dim(myResultsSXYZ_M$sims.list$sxy)
dim(myResultsSXYZ_F$sims.list$sxy)

myResultsSXYZ_MF$sims.list$sxy <- abind(myResultsSXYZ_M$sims.list$sxy[1:dimmFsxy[1],,,],
                                        myResultsSXYZ_F$sims.list$sxy, along = 2 )
dimnames(myResultsSXYZ_MF$sims.list$sxy)[[3]] <- c("x", "y")

myResultsSXYZ_MF$sims.list$z <- abind(myResultsSXYZ_M$sims.list$z[1:dimmFsxy[1],,],
                                      myResultsSXYZ_F$sims.list$z, along = 2 )
myResultsSXYZ_MF$sims.list$z1 <- abind(myResultsSXYZ_M$sims.list$z1[1:dimmFsxy[1],,],
                                       myResultsSXYZ_F$sims.list$z1, along = 2 )
gc()
myResultsSXYZ_MF$sims.list$sigma <- abind(myResults_M$sims.list$sigma[1:dimmF[1],] * res(myHabitat.list$habitat.r)[1],
                                          myResults_F$sims.list$sigma* res(myHabitat.list$habitat.r)[1], along = 2 )


#myResultsSXYZ_MF$sims.list$sigma <- abind(myResultsSXYZ_M$sims.list$sigma[1:dimmF[1]], myResultsSXYZ_F$sims.list$sigma, along = 2 )
myResultsSXYZ_MF$sims.list$sex <- rep(c("M","F"), c(dim(myResultsSXYZ_M$sims.list$sxy)[2], 
                                                    dim(myResultsSXYZ_F$sims.list$sxy)[2]))

#EMPTY USELESS ARRAYS
myResults_F$sims.list$z <-NULL
myResults_M$sims.list$z <-NULL

myResults_M$sims.list$N <-NULL
myResults_F$sims.list$N <-NULL

myResults_F$sims.list$sxy <-NULL
myResults_M$sims.list$sxy <-NULL

gc()

# # #select 100 iterations for RB
# dim(myResultsSXYZ_MF$sims.list$sxy)
# nit <- sample(dim(myResultsSXYZ_MF$sims.list$sxy)[1],100)
# sxy <- round(myResultsSXYZ_MF$sims.list$sxy[nit,,,],digits=5)
# z <- myResultsSXYZ_MF$sims.list$z[nit,,]
# 
# save(sxy, z,
#      file = file.path(WDFigures, "Itera.RData"))
# #contains 100 iterations of "sxy" and "z"
# load("C://Users//cymi//Dropbox (Old)//AQEG Dropbox//AQEG Team Folder//RovQuant//wolverine//CM//2024/plot53Cleaned2024/Figure/Itera.RData")
# 

#MERGE RESULTS IN A LIST
Results.list <- list(myResults_F,myResults_M)
names(Results.list) <- c("F","M")

## ------   1.6 SAVE AND LOAD DATA ------
# save(Results.list, myResultsSXYZ_MF,
#     file = file.path(paste(dir.dropbox,"/wolverine/CM/2022/plot25Cleaned/Figure/",sep=""), "MCMC.RData" ))
# load(file.path(paste(dir.dropbox,"/wolverine/CM/2022/plot25Cleaned/Figure/",sep=""), "MCMC.RData" ))

Results.list[["F"]]$mean$sigma
Results.list[["M"]]$mean$sigma

myResults_F <- Results.list[["F"]]
myResults_M <- Results.list[["M"]]

## ------   1.7 CHECK RHAT ------
myResults_F$Rhat
myResults_M$Rhat
# basicMCMCplots::chainsPlot(nimOutput,var = "omeg1")
basicMCMCplots::chainsPlot(nimOutput,var = "N[1]")
basicMCMCplots::chainsPlot(nimOutput,var = "betaDens")

# WHAT STATES ARE CONSIDERED AS ALIVE IN THE MODEL
alive.states <- c(2) 

#years not sampled in Norrbotten
yearsSampledNorrb <- c(2016:2018,2023,2024)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)


## ------ 2. AC BASED DENSITY  (5km) ------
### REMOVE THE BUFFER FROM THE HABITAT ###
### COUNTRIES 
rrCountries <- habitatRasterResolution$`5km`[["Countries"]]
# REMOVE FINLAND AND RUSSIA 
rrCountries[rrCountries%in% c(1,3)] <- NA
plot(rrCountries)

habitat.rWthBuffer <- myHabitat.list$habitat.rWthBuffer
habitat.rWthBuffer[habitat.rWthBuffer %in% 0] <- NA

searchedPolygon <- sf::st_as_sf(stars::st_as_stars(habitat.rWthBuffer), 
                                as_points = FALSE, merge = TRUE)
searchedPolygon <- searchedPolygon[searchedPolygon$Habitat>0,]
# searchedPolygon <- rasterToPolygons(habitat.rWthBuffer, dissolve = T, function(x) x==1 )

rrCountries <- mask(rrCountries, searchedPolygon)
rrCountries <- crop(rrCountries, myHabitat.list$habitat.r)
plot(rrCountries)
### REGIONS AND COUNTIES 
rrRegions <- habitatRasterResolution$`5km`[["Regions"]] 
## deal with the special characters
levels(rrRegions)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <-  c("Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
                                                                    "Värmland","Örebro","Västmanland","Gävleborg",
                                                                    "Västernorrland","Jämtland" ,"Västerbotten")

# REMOVE FINLAND AND RUSSIA 
rrRegions[habitatRasterResolution$`5km`[["Countries"]]%in% c(1, 3)] <- NA
plot(rrRegions)
rrRegions <- mask(rrRegions, searchedPolygon)
rrRegions <- crop(rrRegions, myHabitat.list$habitat.r)

# habitatPolygon <- rasterToPolygons(myHabitat.list$habitat.r, dissolve = T, function(x) x==1 )
habitatPolygon <- sf::st_as_sf(stars::st_as_stars(myHabitat.list$habitat.r), 
                               as_points = FALSE, merge = TRUE)
habitatPolygon <- habitatPolygon[habitatPolygon$Habitat>0,]

habitatPolygon5km <- mask(habitatRasterResolution$`5km`[["Habitat"]], habitatPolygon)
habitatPolygon5km <- crop(habitatRasterResolution$`5km`[["Habitat"]], myHabitat.list$habitat.r)

plot(rrRegions)
plot(habitatPolygon5km)
# points(myDetectors$main.detector.sp,pch=16,cex=0.5)


### SWEDISH REGIONS 
rrRegionsSwe <- habitatRasterResolution$`5km`[["Regions"]]
# REMOVE FINLAND AND RUSSIA 
rrRegionsSwe[habitatRasterResolution$`5km`[["Countries"]]%in% c(1,2, 3)] <- NA
plot(rrRegionsSwe)
rrRegionsSwe <- mask(rrRegionsSwe, searchedPolygon)
rrRegionsSwe <- crop(rrRegionsSwe, myHabitat.list$habitat.r)

plot(rrRegionsSwe)
rrRegionsSwe[]

# dfLevel <- levels(rrRegionsSwe)[[1]]
# dfLevel[c(18,19,20,21),"Regions"] <- "1"
# dfLevel[c(13,17,16,14,15,12,22,3),"Regions"] <- "Midtre"
# dfLevel[c(4,5,10,6,7,9,11,8),"Regions"] <- "SÃ¸ndre"
# 

rrRegionsSwe[rrRegionsSwe[]%in% c(18,19,20,21)] <- 1
rrRegionsSwe[rrRegionsSwe[]%in% c(13,17,16,14,15,12,22,3)] <- 2
rrRegionsSwe[rrRegionsSwe[]%in% c(4,5,10,6,7,9,11,8)] <- 3

rrRegionsSwe <- ratify(rrRegionsSwe)
df <- data.frame("ID"=c(1,2,3), "Regions"= c("Nordre","Midtre","SÃ¸ndre"))
levels(rrRegionsSwe)[[1]] <- df

##NorwegianCounty 
rrCountiesNor <- habitatRasterResolution$`5km`[["Counties"]]
rrCountiesNor[rrCountiesNor[] %in% c(1,2, 3)] <- NA
rrCountiesNor <- mask(rrCountiesNor, searchedPolygon)
rrCountiesNor <- crop(rrCountiesNor, myHabitat.list$habitat.r)
#remove sweden
rrCountiesNor[rrCountries[]%in%4] <- NA
plot(rrCountiesNor)
plot(rrRegions)
plot(rrCountries)


#library(lattice)
#lattice::levelplot(rrRegions, col.regions=rev(terrain.colors(20)), xlab="", ylab="")
##
# Nordre
# Jämtland
# Västernorrland
##
# Midtre
# Värmland
# Gävleborg
# Dalarna
# Örebro
# Västmanland
# Västra Götaland
# Uppsala
# Stockholm
##
# Søndre
# Södermanland
# Östergötland
# Jönköping
# Skåne

gc()
### GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
## COUNTRY
densityInputCountries <- getDensityInput( regions = rrCountries
                                          , 
                                          habitat = habitatPolygon5km
                                          ,
                                          s = myResultsSXYZ_MF$sims.list$sxy
                                          ,
                                          plot.check = TRUE
)
### GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
## REGIONS
densityInputRegions <- getDensityInput( regions = rrRegions
                                        , 
                                        habitat = habitatPolygon5km
                                        ,
                                        s = myResultsSXYZ_MF$sims.list$sxy
                                        ,
                                        plot.check = TRUE
)

#Swedish Regions
densityInputRegionsSwe <- getDensityInput( regions = rrRegionsSwe
                                           , 
                                           habitat = habitatPolygon5km
                                           ,
                                           s = myResultsSXYZ_MF$sims.list$sxy
                                           ,
                                           plot.check = TRUE
)

## MERGE COUNTRY AND REGION MATRICES TO ALLOW SIMULTANEOUS EXTRACTION
regionID <- rbind (densityInputCountries$regions.rgmx,
                   densityInputRegions$regions.rgmx,
                   densityInputRegionsSwe$regions.rgmx)
row.names(regionID) <- c(row.names(densityInputCountries$regions.rgmx),
                         row.names(densityInputRegions$regions.rgmx),
                         row.names(densityInputRegionsSwe$regions.rgmx))


#Norwegian Counties
densityInputRegionsNor <- getDensityInput( regions = rrCountiesNor
                                           , 
                                           habitat = habitatPolygon5km
                                           ,
                                           s = myResultsSXYZ_MF$sims.list$sxy
                                           ,
                                           plot.check = TRUE
)


## ------  2.1.1 GET THE % OF REGIONS COVERED BY THE ANALYSIS ------
## Percentage of each county included in the analysis
## GET THE PERCENTAGE FOR ALL REGIONS 
rrRegions1 <- habitatRasterResolution$`5km`[["Regions"]] 
## deal with the special characters
levels(rrRegions1)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <-  c("Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
                                                                     "Värmland","Örebro","Västmanland","Gävleborg",
                                                                     "Västernorrland","Jämtland" ,"Västerbotten")
## CALCULATE AREA OF EACH COUNTY
AreaStudiedRegion <- table(factorValues(rrRegions, rrRegions[]))*res(rrRegions)[1]*1e-6
TotalArea <- table(factorValues(rrRegions1, rrRegions1[]))*res(rrRegions1)[1]*1e-6
#Percentage of counties included in the analysis
areaRegions <- AreaStudiedRegion/TotalArea[names(AreaStudiedRegion)]

## GET THE PERCENTAGE FOR THE 3 SWEDISH UNITS REGIONS 
rrRegionsSwe1 <- habitatRasterResolution$`5km`[["Regions"]]
# REMOVE FINLAND AND RUSSIA 
rrRegionsSwe1[habitatRasterResolution$`5km`[["Countries"]]%in% c(1,2, 3)] <- NA
plot(rrRegionsSwe1)

rrRegionsSwe1[rrRegionsSwe1[]%in% c(18,19,20,21)] <- 1
rrRegionsSwe1[rrRegionsSwe1[]%in% c(13,17,16,14,15,12,22,3)] <- 2
rrRegionsSwe1[rrRegionsSwe1[]%in% c(4,5,10,6,7,9,11,8)] <- 3

rrRegionsSwe1 <- ratify(rrRegionsSwe1)
df <- data.frame("ID"=c(1,2,3), "Regions"= c("Nordre","Midtre","SÃ¸ndre"))
levels(rrRegionsSwe1)[[1]] <- df
plot(rrRegionsSwe1)
## CALCULATE AREA OF EACH UNIT
AreaStudiedSwe <- table(factorValues(rrRegionsSwe, rrRegionsSwe[]))*res(rrRegionsSwe)[1]*1e-6
TotalAreaSwe <- table(factorValues(rrRegionsSwe1, rrRegionsSwe1[]))*res(rrRegionsSwe1)[1]*1e-6#Percentage of counties included in the analysis
areaRegionsSwe <- AreaStudiedSwe/TotalAreaSwe[names(AreaStudiedSwe)]

## GET THE PERCENTAGE FOR THE COUNTRIES 
rrCountries1 <- habitatRasterResolution$`5km`[["Countries"]]
# REMOVE FINLAND AND RUSSIA 
rrCountries1[rrCountries1%in% c(1,3)] <- NA
TotalAreaCountry <- table(factorValues(rrCountries1, rrCountries1[]))*res(rrCountries1)[1]*1e-6
AreaStudiedCountry <- table(factorValues(rrCountries, rrCountries[]))*res(rrCountries)[1]*1e-6
TotalAreaCountry["Total"] <- sum(TotalAreaCountry)
AreaStudiedCountry["Total"] <- sum(AreaStudiedCountry)
## CALCULATE AREA OF EACH COUNTRY
areaCountry <- AreaStudiedCountry/TotalAreaCountry[names(AreaStudiedCountry)]

### MERGE THE PERCENTAGE 
areaAllRegions <- c(areaCountry, areaRegionsSwe, areaRegions)

## ------   2.1 MALE AND FEMALES (5km) ------
ite <- seq(1,dim(densityInputCountries$sx[,,t])[1],by=1)
gc()
## EXTRACT DENSITY 
DensityCountriesRegions <- list()
for(t in 1:nYears){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,,t],
    sy =  densityInputCountries$sy[ite,,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
  gc()
}

DensityCountriesRegions[[t]]$summary

## SAVE
#save(DensityCountriesRegions, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegions.RData" ))

## ------   2.2 MALE (5km) ------
IDMales <- which(myResultsSXYZ_MF$sims.list$sex=="M")

DensityCountriesRegionsM <- list()
for(t in 1:nYears){
  DensityCountriesRegionsM[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDMales,t],
    sy =  densityInputCountries$sy[ite,IDMales,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,IDMales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
  gc()
}

DensityCountriesRegionsM[[t]]$summary

## SAVE
#save(DensityCountriesRegionsM, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegionsM.RData" ))

## ------   2.3 FEMALE (5km) ------
IDFemales <- which(myResultsSXYZ_MF$sims.list$sex=="F")

DensityCountriesRegionsF <- list()
for(t in 1:nYears){
  DensityCountriesRegionsF[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDFemales,t],
    sy =  densityInputCountries$sy[ite,IDFemales,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,IDFemales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
  gc()
}

DensityCountriesRegionsF[[t]]$summary
## SAVE
#save(DensityCountriesRegionsF, file = file.path(paste(dir.dropbox, "/wolverine/CM/2021/plot25Cleaned/Figure/",sep=""), "DensityCountriesRegionsF.RData" ))


## ------   2.4 COUNTIES NORWAY M AND F (5km)  ------
DensityCountriesRegionsNOR <- list()
for(t in 1:nYears){
  DensityCountriesRegionsNOR[[t]] <- GetDensity_PD(
    sx = densityInputRegionsNor$sx[ite,,t],
    sy =  densityInputRegionsNor$sy[ite,,t],
    z = myResultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputRegionsNor$habitat.id,
    aliveStates = alive.states,
    regionID = densityInputRegionsNor$regions.rgmx,
    returnPosteriorCells = F)
  gc()
}

# sum(DensityCountriesRegionsNOR[[t]]$summary[14:31,"mean"])
# sum(DensityCountriesRegionsNOR[[t]]$summary[14:31,"mean"])
# 



##



## ------   2.4 SUMMARY TABLES ------
## ------    2.4.1 ALL YEARS, BOTH SEX ------
idcounty <- row.names(DensityCountriesRegions[[t]]$summary)
#REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
#GET NORWEGIAN VERSUS SWEDISH COUNTIES 
idcountyNOR <- idcounty[grep("Region",idcounty)]
idcountySWE <- sort(idcounty[-grep("Region",idcounty)])
idcountyTable <- c("Total","Norway", idcountyNOR, "Sweden" ,idcountySWE)

CountyNorth <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 1], layer=1)[,1])
CountyMiddle <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 2], layer=1)[,1])
CountySouth <- unique(factorValues(rrRegions, rrRegions[rrRegionsSwe[] %in% 3], layer=1)[,1])

idcountyTable <- c("Total","Norway",
                   idcountyNOR,
                   "Sweden" ,
                   "Nordre",
                   idcountySWE[idcountySWE%in%CountyNorth],
                   "Midtre",
                   idcountySWE[idcountySWE%in%CountyMiddle],
                   "SÃ¸ndre",
                   idcountySWE[idcountySWE%in%CountySouth]
)


#CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimates <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimates) <- c(idcountyTable)
colnames(NCarRegionEstimates) <-  unlist(lapply(YEARS,function(x) c(paste(x, collapse = "/"))))#unlist(lapply(YEARS ,function(x) x[2]))#

#FILL IN THE TABLE 
for(t in 1:nYears){
  for( i in 1:length(idcountyTable)){
    NCarRegionEstimates[idcountyTable[i],t] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                     " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                     round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}

##QUICK CHECK TO MAKE SURE VALUES SUMS UP 
tmp <- DensityCountriesRegions[[t]]$summary[1:(nrow(DensityCountriesRegions[[t]]$summary)),]
# SWE
row.names(DensityCountriesRegions[[t]]$summary)
sum(tmp[idcountySWE,"mean"])
tmp["Sweden","mean"]
#NOR
sum(tmp[idcountyNOR,"mean"])
tmp["Norway","mean"]
#TOTAL
sum(tmp[c(idcountyNOR,idcountySWE),"mean"])
tmp["Total","mean"]



## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")
idcountySWE1 <- idcountySWE

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
# 
# row.names(NCarRegionEstimates) <- idcounty1
# NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste(NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*", sep="")

#print csv
write.csv(NCarRegionEstimates,
          file = file.path(WDTables,paste("NAllYears.csv",sep="")),fileEncoding="latin1")





idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten}"
row.names(NCarRegionEstimates) <- idcounty1
NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*","}", sep="")

NCarRegionEstimates["SWEDEN",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["SWEDEN",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["TOTAL",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["TOTAL",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["Nordre",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["Nordre",yearsNotSampled], "**","}", sep="")


row.names(NCarRegionEstimates) <- c("TOTAL",
                                    paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                    paste("\\hspace{0.5cm} ",
                                          idcountyNOR,sep=""),
                                    paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                    paste("\\hspace{0.5cm}","Norra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                    paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                    paste("\\hspace{0.5cm}","Södra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCarRegionEstimates)[grep("VÃ¤straGÃ¶taland", row.names(NCarRegionEstimates))] <- paste("\\hspace{0.75cm}",
                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

# row.names(NCarRegionEstimates) <- c("TOTAL",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1,sep="")
# )


print(xtable(NCarRegionEstimates, type = "latex",align=paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(NCarRegionEstimates),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("NCountiesCarnivoreRegions.tex",sep="")))


#####
NCarRegionEstimatesOPSCR <- read.csv(file.path(WDTables,paste("NAllYears.csv",sep="")))
NCarRegionEstimates1 <- NCarRegionEstimates
NCarRegionEstimates1[,]

NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*","}", sep="")

NCarRegionEstimates["SWEDEN",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["SWEDEN",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["TOTAL",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["TOTAL",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["Nordre",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["Nordre",yearsNotSampled], "**","}", sep="")


row.names(NCarRegionEstimates) <- c("TOTAL",
                                    paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                    paste("\\hspace{0.5cm} ",
                                          idcountyNOR,sep=""),
                                    paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                    paste("\\hspace{0.5cm}","Norra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                    paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                    paste("\\hspace{0.5cm}","Södra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCarRegionEstimates)[grep("VÃ¤straGÃ¶taland", row.names(NCarRegionEstimates))] <- paste("\\hspace{0.75cm}",
                                                                                                  "VÃ¤stra GÃ¶taland", sep="")




## ------    2.4.2 LAST YEAR N PER SEX PER COUNTY  ------
NCountyEstimatesLastRegions <- matrix("", ncol=3, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste(round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                   " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                   round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste(round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste(round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}



#print csv
write.csv(NCountyEstimatesLastRegions,
          file = file.path(WDTables,paste("NLastYearPerSex.csv",sep="")),fileEncoding="latin1")



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

idcountySWE1 <- idcountySWE
idcountySWE1 <- sort(idcountySWE1)

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"


row.names(NCountyEstimatesLastRegions) <- idcounty1
NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),], "*}", sep="")

NCountyEstimatesLastRegions["SWEDEN",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",], "**}", sep="")
NCountyEstimatesLastRegions["TOTAL",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",], "**}", sep="")
NCountyEstimatesLastRegions["Nordre",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",], "**}", sep="")


row.names(NCountyEstimatesLastRegions) <- c("TOTAL",
                                            paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                            paste("\\hspace{0.5cm} ",
                                                  idcountyNOR,sep=""),
                                            paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                            paste("\\hspace{0.5cm}","Norra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                            paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                            paste("\\hspace{0.5cm}","Södra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCountyEstimatesLastRegions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLastRegions))] <- paste("\\hspace{0.75cm}",
                                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


# WRITE LATEX 
print(xtable(NCountyEstimatesLastRegions, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path(WDTables,paste("NCountiesSexLastYearRegions.tex",sep="")))



## ------    2.4.2 LAST YEAR N PER SEX PER COUNTY WITH PROPORTION OF AREA COVERED  ------
NCountyEstimatesLastRegions <- matrix("", ncol=4, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total","\\% Area")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste(round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                   " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                   round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste(round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste(round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                 " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                 round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
}



NCountyEstimatesLastRegions[names(areaAllRegions),"\\% Area"] <- round(areaAllRegions*100,digits = 0)
NCountyEstimatesLastRegions[NCountyEstimatesLastRegions[,4] %in% c("98","99"),4] <- 100
#NCountyEstimatesLastRegions[names(AreaStudied/TotalArea[names(AreaStudied)]),"Area"] <- round(AreaStudied/TotalArea[names(AreaStudied)]*100,digits = 0)

#print csv
write.csv(NCountyEstimatesLastRegions,
          file = file.path(WDTables,paste("NLastYearPerSexArea.csv",sep="")),fileEncoding="latin1")



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable

idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

idcountySWE1 <- idcountySWE
idcountySWE1 <- sort(idcountySWE1)

# idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"


row.names(NCountyEstimatesLastRegions) <- idcounty1
# NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),1:3], "*}", sep="")
# 
# NCountyEstimatesLastRegions["SWEDEN",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["TOTAL",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["Nordre",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",1:3], "**}", sep="")


row.names(NCountyEstimatesLastRegions) <- c("TOTAL",
                                            paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                            paste("\\hspace{0.5cm} ",
                                                  idcountyNOR,sep=""),
                                            paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                            paste("\\hspace{0.5cm}","Norra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                            paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                            paste("\\hspace{0.5cm}","Södra",sep=""),
                                            paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCountyEstimatesLastRegions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLastRegions))] <- paste("\\hspace{0.75cm}",
                                                                                                                  "VÃ¤stra GÃ¶taland", sep="")

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


# WRITE LATEX 
print(xtable(NCountyEstimatesLastRegions, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions)-1),"||c"),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path(WDTables,paste("NCountiesSexLastYearRegionsArea.tex",sep="")))




## ------    2.4.3 MAKE A TABLE 2 last years  ------
NCountyEstimatesLast2Regions <- matrix("", ncol=6, nrow=length(idcountyTable))
row.names(NCountyEstimatesLast2Regions) <- c(idcountyTable)
colnames(NCountyEstimatesLast2Regions) <- c(paste("Females", years[nYears-1]),
                                            paste("Males", years[nYears-1]),
                                            paste("Total", years[nYears-1]),
                                            paste("Females", years[nYears]),
                                            paste("Males", years[nYears]),
                                            paste("Total", years[nYears])
                                            
)




## FILL IN TABLE 
for(t in (nYears-1):nYears){
  ## FEMALES
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Females",years[t])] <- paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                      " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Males",years[t])] <- paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Total",years[t])] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}



# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")


#print csv
write.csv(NCountyEstimatesLast2Regions,
          file = file.path(WDTables,paste("NLast2YearsPerSex.csv",sep="")),fileEncoding="latin1")



##ADD LITLE STAR TO NORRBOTTEN
idcountySWE1 <- idcountySWE

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
row.names(NCountyEstimatesLast2Regions) <- idcounty1
NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),] <- paste(NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),], "*", sep="")





row.names(NCountyEstimatesLast2Regions) <- c("TOTAL",
                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                             paste("\\hspace{0.5cm} ",
                                                   idcountyNOR,sep=""),
                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                             paste("\\hspace{0.5cm}","Norra**",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                             paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                             paste("\\hspace{0.5cm}","Södra",sep=""),
                                             paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)
row.names(NCountyEstimatesLast2Regions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLast2Regions))] <- paste("\\hspace{0.75cm}",
                                                                                                                    "VÃ¤stra GÃ¶taland", sep="")


NCountyEstimatesLast2Regions <- rbind(c("F","M","Total","F","M","Total"), NCountyEstimatesLast2Regions)

# WRITE LATEX 
addtorow <- list()

addtorow$pos <- list(c(0),0)
uniqueYEAR <- c(paste(unlist(YEARS[nYears-1]),collapse = "/"),paste(unlist(YEARS[nYears]),collapse = "/"))
addtorow$command <- c(paste0(paste0('& \\multicolumn{3}{c}{', uniqueYEAR,
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

print(xtable(NCountyEstimatesLast2Regions, type = "latex",
             align = paste(c("l",rep("c",3),"|",rep("c",3)),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      file = file.path(WDTables,paste("NCountiesSexLast2YearsRegions.tex",sep="")))

## ------    2.4.4 ALL YEARS N PER SEX PER COUNTY  ------
NCountyEstimatesAllSexRegions <- matrix("", ncol=nYears*3, nrow=length(idcountyTable)+1)
row.names(NCountyEstimatesAllSexRegions) <- c("",idcountyTable)
colnames(NCountyEstimatesAllSexRegions) <- rep(unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/")))),each=3)
NCountyEstimatesAllSexRegions[1,] <- rep(c("Females","Males","Total"),nYears)
## FILL IN TABLE 
for(t in 1:nYears){
  
  ## FEMALES
  cols <- which(colnames(NCountyEstimatesAllSexRegions) %in% unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/"))))[t])
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Females")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## MALES 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
  
  ## TOTAL 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                         " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                         round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}

#print csv
write.csv(NCountyEstimatesAllSexRegions,
          file = file.path(WDTables,paste("NAllYearsPerSex.csv",sep="")),fileEncoding="latin1")



# # ADJUST NAMES OF THE TABLE 
# idcounty1 <- idcountyTable
# 
# idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
# idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
# idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"
# 
# idcountySWE1 <- idcountySWE
# idcountySWE1 <- sort(idcountySWE1)
# 
# idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"
# row.names(NCountyEstimatesLastRegions) <- idcounty1
# NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),], "*}", sep="")
# 
# 
# row.names(NCountyEstimatesLastRegions) <- c("TOTAL**",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN**",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1, sep="")
# )
# 
# ## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# # idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# # idcounty1 <- str_remove(idcounty1, "s ")
# # idcounty1 <- str_remove(idcounty1, " ")
# 
# 
# # WRITE LATEX 
# print(xtable(NCountyEstimatesLastRegions, type = "latex",
#              align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))),collapse = "")),
#       sanitize.text.function=function(x){x},
#       # scalebox=.8,
#       floating = FALSE,
#       add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
#       file = file.path(WDTables,paste("NCountiesSexLastYearRegions.tex",sep="")))
## ------    2.4.5 ALL YEARS, BOTH SEX COUNTIES NORWAY ------
idcounty <- row.names(DensityCountriesRegionsNOR[[t]]$summary)
#REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
#GET NORWEGIAN VERSUS SWEDISH COUNTIES 
idcountyNOR <- idcounty[grep("Region",idcounty)]
#idcountySWE <- sort(idcounty[-grep("Region",idcounty)])
idcountyTable <- c("Total","Norway", idcountyNOR, "Sweden" ,idcountySWE)


idcountyTable <- c("Total",
                   idcountyNOR
)


#CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimatesNOR <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimatesNOR) <- c(idcountyTable)
colnames(NCarRegionEstimatesNOR) <- unlist(lapply(YEARS ,function(x) c(paste(x, collapse = "/"))))#  unlist(lapply(YEARS ,function(x) x[2]))#

#FILL IN THE TABLE 
for(t in 1:nYears){
  for( i in 1:length(idcountyTable)){
    NCarRegionEstimatesNOR[idcountyTable[i],t] <- paste(round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                        " (",round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                        round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }
}



##QUICK CHECK TO MAKE SURE VALUES SUMS UP 
tmp <- DensityCountriesRegionsNOR[[t]]$summary#[1:(nrow(DensityCountriesRegions[[t]]$summary)),]
# SWE
row.names(DensityCountriesRegions[[t]]$summary)


# sum(tmp[row.names(tmp) %in% idcountySWE,"mean"])
sum(tmp[row.names(tmp) %in% idcountyNOR,"mean"])

# tmp["Sweden","mean"]
# #NOR
# sum(tmp[idcountyNOR,"mean"])
# tmp["Norway","mean"]
# #TOTAL
# sum(tmp[c(idcountyNOR,idcountySWE),"mean"])
# tmp["Total","mean"]
# 


## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1 <- gsub("Region ", "", idcountyTable)
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"

## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")

# row.names(NCarRegionEstimates) <- idcounty1
# NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste(NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*", sep="")
row.names(NCarRegionEstimatesNOR) <- idcounty1

#print csv
# NCarRegionEstimatesNOR <- data.frame(NCarRegionEstimatesNOR)
# NCarRegionEstimatesNOR$name <- row.names(NCarRegionEstimatesNOR)
# Encoding(NCarRegionEstimatesNOR[1,"name"]) <- "UTF-16"#"UTF-16"
#save(NCarRegionEstimatesNOR,file=file.path(WDTables,paste("NAllYearsNorwegianCounties.RData",sep="")))
write.csv(NCarRegionEstimatesNOR,
          file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")),fileEncoding= "latin1")

# Encoding(NCarRegionEstimatesNOR[,"name"])[9] <- "ISO-8859-1"
# mb_convert_encoding($file, 'UTF-8', 'ISO-8859-1')
# write.csv2(NCarRegionEstimatesNOR,
#            file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")),fileEncoding= "UTF-16LE")
# readr::write_excel_csv(NCarRegionEstimatesNOR,
#                         file = file.path(WDTables,paste("NAllYearsNorwegianCounties.csv",sep="")))
## try to join to the Norwegian layer for Richard
# tmp <- data.frame(NCarRegionEstimatesNOR)
# tmp$NAME_1 <- row.names(tmp) 
# COUNTIES_1 <- merge(COUNTIES,tmp[,c("X2024","NAME_1")],by="NAME_1")
# 




row.names(NCarRegionEstimatesNOR) <- c("NORWAY",
                                       # paste("\\hspace{0.25cm}","NORWAY",sep=""),
                                       paste("\\hspace{0.25cm} ",
                                             idcountyNOR,sep="")#
                                       # paste("\\hspace{0.25cm}","SWEDEN",sep=""),
                                       # paste("\\hspace{0.5cm}","Norra",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                       # paste("\\hspace{0.5cm}","Mellersta",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                       # paste("\\hspace{0.5cm}","Södra",sep=""),
                                       # paste("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
)

row.names(NCarRegionEstimatesNOR)[grep("VÃ¤straGÃ¶taland", row.names(NCarRegionEstimatesNOR))] <- paste("\\hspace{0.25cm}",
                                                                                                        "VÃ¤stra GÃ¶taland", sep="")

# row.names(NCarRegionEstimates) <- c("TOTAL",
#                                             paste("\\hspace{0.25cm}","NORWAY",sep=""),
#                                             paste("\\hspace{0.5cm} ",
#                                                   idcountyNOR,sep=""),
#                                             paste("\\hspace{0.25cm}","SWEDEN",sep=""),
#                                             paste("\\hspace{0.5cm}", idcountySWE1,sep="")
# )


print(xtable(NCarRegionEstimatesNOR, type = "latex",align=paste(c("l",rep("c",ncol(NCarRegionEstimatesNOR))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(NCarRegionEstimatesNOR),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("NCountiesCarnivoreRegionsNorway.tex",sep="")))



## ------   2.5 PLOT ABUNDANCE TIME SERIES ------
## ------    2.5.1 VIOLINS ------
## ------      2.5.1.1 ALL YEARS ------
SeasonText <- lapply(YEARS,FUN = function(x) paste(x,collapse = "/"))

#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha <- 1

pdf(file= file.path(WDFigures, paste("NCountriesViolins.pdf", sep="")), width = 12, height = 8)
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.8, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+0.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolves"), xaxt="n")
axis(1, at=c(1:nYears), labels = SeasonText, cex.axis=1.2)
at = c(1:nYears)




for(t in 1:nYears){
  xx <- t
  yy <- n.detected[[t]]
  xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
  yy<-c(0,0,yy,yy)
  polygon(xx, yy ,border=NA,col=grey(0.9))
  
}


for(t in 1:nYears){
  #TOTAL
  plot.violins2(list(colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)),
                x = at[t]+0.05,
                at= at[t]+0.05,
                violin.width = 0.3,
                col = TotalColors,
                alpha = violin.alpha,
                border.col = TotalColors,
                add = T
                ,cex=2,median = FALSE)
  
  
  text(round(round(DensityCountriesRegions[[t]]$summary["Total","mean"],digits = 1)),x=t+0.1,
       y=DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]+total.offset, cex=text.cex)
  
  #SWEDEN
  plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]),
                x = at[t],
                at=at[t]-0.05,
                violin.width = 0.3,
                col = country.colors[2],
                alpha = violin.alpha,
                border.col = country.colors[2],
                add = T
                ,cex=2,median = FALSE)
  text(round(round(DensityCountriesRegions[[t]]$summary["Sweden","mean"],digits = 1)),x=t-0.1,
       y=DensityCountriesRegions[[t]]$summary["Sweden","95%CILow"] + SE.offset, cex=text.cex)
  #NORWAY
  #print(t)
  plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]),
                x = at[t],
                at=at[t],
                violin.width = 0.3,
                col = country.colors[1],
                alpha = violin.alpha,
                border.col = country.colors[1],
                add = T,scale.width = FALSE
                ,cex=2,median = FALSE
                
  )
  text(round(round(DensityCountriesRegions[[t]]$summary["Norway","mean"],digits = 1)),x=t+0.1,
       y=DensityCountriesRegions[[t]]$summary["Norway","95%CIHigh"]+NO.offset,cex=text.cex)
  
  
  
}
box()
abline(v=at[1:(nYears-1)]+0.5,lty=2)

#legend
par(xpd=TRUE)
legend(x = 1, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(4, 4, 4),
       horiz = T,
       pch=c(16, 16, 16),
       col=c(country.colors, "black"),
       bty = 'n',
       cex = 1.5)
legend(x = 1, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(1.3,1.3 , 1.3),
       horiz = T,
       pch=c(16,16,16),
       col=c("white", "white", "white"),
       bty = 'n',
       cex = 1.5)

dev.off()

## ------      2.5.1.2 LAST YEAR ------
pdf(file= file.path(WDFigures, paste("NCountriesViolinsLastYear.pdf", sep="")), width = 12, height = 8)

par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.8, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(nYears-0.1, nYears+0.1), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolves"), xaxt="n")
axis(1, at=c(nYears), labels = SeasonText[nYears], cex.axis=1.2)
at = c(1:nYears)




#for(t in nYears){
xx <- t
yy <- n.detected[[nYears]]
xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
yy<-c(0,0,yy,yy)
polygon(xx, yy ,border=NA,col=grey(0.9))

#}

t <- nYears
#TOTAL
plot.violins2(list(colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)),
              x = at[t],
              at= at[t],
              violin.width = 0.02,
              col = TotalColors,
              alpha = violin.alpha,
              border.col = TotalColors,
              add = T
              ,cex=2,median = FALSE)


# text(round(round(DensityCountriesRegions[[t]]$summary["Total","mean"],digits = 1)),x=t+0.1,
#      y=DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]+total.offset, cex=text.cex)
# 
#SWEDEN
plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]),
              x = at[t],
              at=at[t],
              violin.width = 0.02,
              col = country.colors[2],
              alpha = violin.alpha,
              border.col = country.colors[2],
              add = T
              ,cex=2,median = FALSE)
# text(round(round(DensityCountriesRegions[[t]]$summary["Sweden","mean"],digits = 1)),x=t-0.1,
#      y=DensityCountriesRegions[[t]]$summary["Sweden","95%CILow"] + SE.offset, cex=text.cex)
#NORWAY
#print(t)
plot.violins2(list(DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]),
              x = at[t],
              at=at[t],
              violin.width = 0.02,
              col = country.colors[1],
              alpha = violin.alpha,
              border.col = country.colors[1],
              add = T,scale.width = FALSE
              ,cex=2,median = FALSE
              
)
# text(round(round(DensityCountriesRegions[[t]]$summary["Norway","mean"],digits = 1)),x=t+0.1,
#      y=DensityCountriesRegions[[t]]$summary["Norway","95%CIHigh"]+NO.offset,cex=text.cex)
# 


# }
box()
# abline(v=at[1:(nYears-1)]+0.5,lty=2)

#legend
par(xpd=TRUE)
legend(x = nYears+0.01, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(2, 2, 2),
       horiz = T,
       pch=c(16, 16, 16),
       col=c(country.colors, "black"),
       bty = 'n',
       cex = 1.2)
legend(x = nYears+0.01, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       pt.cex = c(0.5,0.5 , 0.5),
       horiz = T,
       pch=c(16,16,16),
       col=c("white", "white", "white"),
       bty = 'n',
       cex = 1.2)

dev.off()


## ------    2.5.2 BARS ------
## ------      2.5.2.1 ALL YEARS ------
SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

# SeasonText <- lapply(YEARS, FUN = function(x) x[2]) #paste(x[[2]]))


#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf(file= file.path(WDFigures, paste("NCountriesBarsSeason.pdf", sep="")),   width = 14, height = 10)
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.2,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 
# n.detected <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv",sep="")))
# n.detected <- as.vector(n.detected[1,2:ncol(n.detected)])
# n.detected[1]
# 
# for(t in 1:nYears){
#   xx <- t
#   yy <- n.detected[1,t]
#   xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
#   yy<-c(0,0,yy,yy)
#   polygon(xx, yy ,border=NA,col=grey(0.9))
#   
# }
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE

quantile95Tot <- quantile50Tot <- list()
quantile95Swe <- quantile50Swe <- list()
quantile95Nor <- quantile50Nor <- list()

#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95Tot[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Tot[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95Tot[[t]][1], quantile95Tot[[t]][1],
                quantile95Tot[[t]][2], quantile95Tot[[t]][2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95Tot[[t]][2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50Tot[[t]][1], quantile50Tot[[t]][1],
                  quantile50Tot[[t]][2], quantile50Tot[[t]][2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
  quantile95Swe[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Swe[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95Swe[[t]] [1], quantile95Swe[[t]] [1],
                quantile95Swe[[t]] [2], quantile95Swe[[t]] [2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95Swe[[t]][2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50Swe[[t]][1], quantile50Swe[[t]][1],
                  quantile50Swe[[t]][2], quantile50Swe[[t]][2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  quantile95Nor[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Nor[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95Nor[[t]][1], quantile95Nor[[t]][1],
                quantile95Nor[[t]][2], quantile95Nor[[t]][2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50Nor[[t]][1], quantile50Nor[[t]][1],
                  quantile50Nor[[t]][2], quantile50Nor[[t]][2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(200, 200, 200)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()

##save it for Post-plotting 
save(quantile95Tot, quantile50Tot,
     quantile95Swe, quantile50Swe,
     quantile95Nor, quantile50Nor, file=file.path(WDFigures, paste("CICounties.RData", sep=""))
)




## ------      2.5.2.2 ALL YEARS SEX ------
SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

#SeasonText <- lapply(YEARS, FUN = function(x) x[2]) #paste(x[[2]]))


#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf(file= file.path(WDFigures, paste("NCountriesBarsSexSeason.pdf", sep="")), width = 18, height = 8)
par(mfrow=c(1,2), mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,800),
     xlab="", ylab = paste("Estimated number of Females"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.1,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(50, 50, 50)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}

##MALES 
#par(mfrow=c(1,2), mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,800),
     xlab="", ylab = paste("Estimated number of Males"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.1,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA
for(t in 1:nYears){
  #TOTAL
  tmp <- colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a starCOUNTRIES <- COUNTRIES %>%    group_by(ISO) %>%summarize()
  
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  tmp <- DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95[2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  #print(t)
  tmp <- DensityCountriesRegionsM[[t]]$PosteriorRegions["Norway",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  # if(sum(t %in% yearsNotSampled)){
  #    text(x=t+widthPolygon1*2+offsetstar ,y= quantile95[2], "*",cex=cexStar)
  #    
  # }
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(600,600,650,650),
        col=adjustcolor("white",alpha.f = 0.9),border="white")

# legend(x = 1, y = 600,
#        legend= c(" Norway  ", " Sweden  ", " Total"),
#        #pt.cex = c(4, 4, 4),
#        horiz = T,
#        #pch=c(16, 16, 16),
#        fill=c(country.colors, "black"),
#        border=NA,
#        bty = 'n',
#        cex = 1.5)



labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(630, 630, 630)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()





## ------      2.5.1.3 LAST YEAR ------
pdf(file= file.path(WDFigures, paste("NCountriesBarsLastYear.pdf", sep="")), width = 12, height = 8)
plot(-1000, xlim=c(nYears-0.1, nYears+0.1), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolves"), xaxt="n")
axis(1, at=c(nYears), labels = SeasonText[nYears], cex.axis=1.2)
at = c(1:nYears)


# #for(t in nYears){
# xx <- t
# yy <- n.detected[1,nYears]
# xx <- c(xx-0.5,xx+0.5,xx+0.5,xx-0.5)
# yy<-c(0,0,yy,yy)
# polygon(xx, yy ,border=NA,col=grey(0.9))
# 
# #}
widthPolygon <- 0.01
widthPolygon1 <- 0.01
widthPolygon2 <- 0.01
violin.alpha <- 0.8
displayQuantiles50 <- FALSE
t <- nYears
#for(t in 1:nYears){
#TOTAL
tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t-widthPolygon, t+widthPolygon,
              t+widthPolygon, t-widthPolygon ),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(TotalColors, violin.alpha),
        border= NA)



if(displayQuantiles50){
  polygon(x = c(t-widthPolygon, t+widthPolygon,
                t+widthPolygon, t-widthPolygon ),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(TotalColors, violin.alpha),
          border= NA)
}
#SWEDEN
tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t, t-widthPolygon1*2,
              t-widthPolygon1*2, t),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(country.colors[2], violin.alpha),
        border= NA)
if(displayQuantiles50){
  polygon(x = c(t, t-widthPolygon1*2,
                t-widthPolygon1*2, t),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(country.colors[2], violin.alpha),
          border= NA)
}

#NORWAY
#print(t)
tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
quantile50 <- quantile(tmp, prob=c(0.25, 0.75))

polygon(x = c(t, t+widthPolygon1*2,
              t+widthPolygon1*2, t),
        y = c(quantile95[1], quantile95[1],
              quantile95[2], quantile95[2]), 
        col=adjustcolor(country.colors[1], violin.alpha),
        border= NA)
if(displayQuantiles50){
  polygon(x = c(t, t+widthPolygon2*2,
                t+widthPolygon2*2, t),
          y = c(quantile50[1], quantile50[1],
                quantile50[2], quantile50[2]), 
          col=adjustcolor(country.colors[1], violin.alpha),
          border= NA)
}


#}
box()
abline(v=at[1:(nYears-1)]+0.5,lty=2)

#legend
par(xpd=TRUE)
legend(x = 1, y = 200,
       legend= c(" Norway  ", " Sweden  ", " Total"),
       #pt.cex = c(4, 4, 4),
       horiz = T,
       #pch=c(16, 16, 16),
       fill=c(country.colors, "black"),
       border=NA,
       bty = 'n',
       cex = 1.5)


dev.off()

## ------      2.5.1.4 UPDATE WITH LAST YEAR RESULTS 10years ------
### LOAD LAST YEAR RESULTS 
load(file.path("C:/Users/cymi/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2024/plot53Cleaned2024/Figure" ,
               "CICounties.RData"))
quantile95Tot2024 <- quantile95Tot
quantile50Tot2024 <- quantile50Tot
quantile95Swe2024 <- quantile95Swe
quantile50Swe2024 <- quantile50Swe
quantile95Nor2024 <- quantile95Nor
quantile50Nor2024 <- quantile50Nor



#define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7
#update season text 
SeasonText1 <- c(list(c("2014/2015")), SeasonText)

pdf(file= file.path(WDFigures,
                    paste("NCountriesBarsOPSCR_SCR.pdf", sep="")),
    width = 14, height = 10)

par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+0.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(1:(length(SeasonText))), labels = SeasonText, cex.axis=1.2,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE

#get the CI from the SCR model 
quantile95Tot <- quantile50Tot <- list()
quantile95Swe <- quantile50Swe <- list()
quantile95Nor <- quantile50Nor <- list()

for(t in 1:nYears){
  tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95Tot[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Tot[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
  quantile95Swe[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Swe[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
  tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  quantile95Nor[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50Nor[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
}

#USE CI from OPSCR model published last year, except for last year. 
quantile95Tot2024[[nYears+1]] <- quantile95Tot[[nYears]]
quantile50Tot2024[[nYears+1]] <- quantile50Tot[[nYears]]
quantile95Swe2024[[nYears+1]] <- quantile95Swe[[nYears]]
quantile50Swe2024[[nYears+1]] <- quantile50Swe[[nYears]]
quantile95Nor2024[[nYears+1]] <- quantile95Nor[[nYears]]
quantile50Nor2024[[nYears+1]] <- quantile50Nor[[nYears]]


for(t in 1:(nYears)){
  #TOTAL
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95Tot2024[[t+1]][1], quantile95Tot2024[[t+1]][1],
                quantile95Tot2024[[t+1]][2], quantile95Tot2024[[t+1]][2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum(t %in% yearsNotSampled)){
    text(x=t+widthPolygon+offsetstar ,y= quantile95Tot2024[[t+1]][2], "*",cex=cexStar)
    
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50Tot2024[[t+1]][1], quantile50Tot2024[[t+1]][1],
                  quantile50Tot2024[[t+1]][2], quantile50Tot2024[[t+1]][2]), 
            col=adjustcolor(TotalColors, violin.alpha50),
            border= NA)
  }
  #SWEDEN
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95Swe2024[[t+1]] [1], quantile95Swe2024[[t+1]] [1],
                quantile95Swe2024[[t+1]] [2], quantile95Swe2024[[t+1]] [2]), 
          col=adjustcolor(country.colors[2], violin.alpha95),
          border= NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95Swe2024[[t+1]][2], "*",cex=cexStar)
    
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50Swe2024[[t+1]][1], quantile50Swe2024[[t+1]][1],
                  quantile50Swe2024[[t+1]][2], quantile50Swe2024[[t+1]][2]), 
            col=adjustcolor(country.colors[2], violin.alpha50),
            border= NA)
  }
  
  #NORWAY
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95Nor2024[[t+1]][1], quantile95Nor2024[[t+1]][1],
                quantile95Nor2024[[t+1]][2], quantile95Nor2024[[t+1]][2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50Nor2024[[t+1]][1], quantile50Nor2024[[t+1]][1],
                  quantile50Nor2024[[t+1]][2], quantile50Nor2024[[t+1]][2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
  
  
}
box()
abline(v=at[1:(nYears)]+0.5,lty=2)

#legend
par(xpd=TRUE)
polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")


labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(200, 200, 200)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}



dev.off()
## ------      2.5.1.4 UPDATE WITH LAST YEAR RESULTS 11yeara ------
# 
# 
# #define colors
# text.cex <- 1.5
# total.offset <- 37
# NO.offset <- -37
# SE.offset <- +37
# xlim <- c(0.5, nYears + 0.5)
# 
# TotalColors <- "black"
# country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
# names(country.colors) <- c("Norway","Sweden")
# violin.alpha95 <- 0.3
# violin.alpha50 <- 0.7
# #update season text 
# SeasonText1 <- c(list(c("2014/2015")), SeasonText)
# 
# pdf(file= file.path(WDFigures,
#                     paste("NCountriesBarsSeasonCombined202411yrs.pdf", sep="")),
#     width = 14, height = 10)
# 
# par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
# plot(-1000, xlim=c(0.5, nYears+1.5), ylim=c(0,1300),
#      xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
# axis(1, at=c(1:(length(SeasonText1))), labels = SeasonText1, cex.axis=1.2,padj = -1)
# at = c(1:nYears)
# abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))
# 
# widthPolygon <- 0.15
# widthPolygon1 <- 0.15
# widthPolygon2 <- 0.15
# offsetstar <- 0.05
# cexStar <- 1.5
# displayQuantiles50 <- TRUE
# 
# quantile95Tot <- quantile50Tot <- list()
# quantile95Swe <- quantile50Swe <- list()
# quantile95Nor <- quantile50Nor <- list()
# 
# for(t in 1:nYears){
# tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
#   quantile95Tot[[t+1]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Tot[[t+1]] <- quantile(tmp, prob=c(0.25, 0.75))
# tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
#   quantile95Swe[[t+1]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Swe[[t+1]] <- quantile(tmp, prob=c(0.25, 0.75))
# tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
#   quantile95Nor[[t+1]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Nor[[t+1]] <- quantile(tmp, prob=c(0.25, 0.75))
# }
# 
# quantile95Tot[[1]] <- quantile95Tot2024[[1]]
# quantile50Tot[[1]] <- quantile50Tot2024[[1]]
# quantile95Swe[[1]] <- quantile95Swe2024[[1]]
# quantile50Swe[[1]] <- quantile50Swe2024[[1]]
# quantile95Nor[[1]] <- quantile95Nor2024[[1]]
# quantile50Nor[[1]] <- quantile50Nor2024[[1]]
#   
# 
# 
#   
# for(t in 1:(nYears+1)){
#   #TOTAL
#   polygon(x = c(t - widthPolygon, t + widthPolygon,
#                 t + widthPolygon, t - widthPolygon ),
#           y = c(quantile95Tot[[t]][1], quantile95Tot[[t]][1],
#                 quantile95Tot[[t]][2], quantile95Tot[[t]][2]), 
#           col=adjustcolor(TotalColors, violin.alpha95),
#           border= NA)
#   
#   # add a star
#   if(sum(t %in% yearsNotSampled)){
#     text(x=t+widthPolygon+offsetstar ,y= quantile95Tot[[t]][2], "*",cex=cexStar)
#     
#   }
#   
#   if(displayQuantiles50){
#     polygon(x = c(t-widthPolygon, t+widthPolygon,
#                   t+widthPolygon, t-widthPolygon ),
#             y = c(quantile50Tot[[t]][1], quantile50Tot[[t]][1],
#                   quantile50Tot[[t]][2], quantile50Tot[[t]][2]), 
#             col=adjustcolor(TotalColors, violin.alpha50),
#             border= NA)
#   }
#   #SWEDEN
#   polygon(x = c(t, t -widthPolygon1*2,
#                 t - widthPolygon1*2, t),
#           y = c(quantile95Swe[[t]] [1], quantile95Swe[[t]] [1],
#                 quantile95Swe[[t]] [2], quantile95Swe[[t]] [2]), 
#           col=adjustcolor(country.colors[2], violin.alpha95),
#           border= NA)
#   
#   if(sum(t %in% yearsNotSampled)){
#     text(x= t+offsetstar ,y= quantile95Swe[[t]][2], "*",cex=cexStar)
#     
#   }
#   if(displayQuantiles50){
#     polygon(x = c(t, t-widthPolygon1*2,
#                   t-widthPolygon1*2, t),
#             y = c(quantile50Swe[[t]][1], quantile50Swe[[t]][1],
#                   quantile50Swe[[t]][2], quantile50Swe[[t]][2]), 
#             col=adjustcolor(country.colors[2], violin.alpha50),
#             border= NA)
#   }
#   
#   #NORWAY
#   polygon(x = c(t, t+widthPolygon1*2,
#                 t+widthPolygon1*2, t),
#           y = c(quantile95Nor[[t]][1], quantile95Nor[[t]][1],
#                 quantile95Nor[[t]][2], quantile95Nor[[t]][2]), 
#           col=adjustcolor(country.colors[1], violin.alpha95),
#           border= NA)
# 
#   if(displayQuantiles50){
#     polygon(x = c(t, t+widthPolygon2*2,
#                   t+widthPolygon2*2, t),
#             y = c(quantile50Nor[[t]][1], quantile50Nor[[t]][1],
#                   quantile50Nor[[t]][2], quantile50Nor[[t]][2]), 
#             col=adjustcolor(country.colors[1], violin.alpha50),
#             border= NA)
#   }
#   
#   
# }
# box()
# abline(v=at[1:(nYears)]+0.5,lty=2)
# 
# #legend
# par(xpd=TRUE)
# polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")
# 
# 
# labels <- c(" Norway  ", " Sweden  ", " Total")
# pch <- rep(19,4)
# cex <- rep(4,4)
# y <-  c(200, 200, 200)
# x <- c(1,3,5)
# mycol1 <- c(country.colors, "black")
# # add transparent background polygon
# #polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
# for(i in 1:3){
#   #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
#   points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
#   points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
#   text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
# }
# 
# 
# 
# dev.off()
# 
## ------      2.5.1.5 UPDATE ALL YEARS WITH LAST YEAR RESULTS ------

# pdf(file= file.path(WDFigures,
#                     paste("NCountriesBarsSeason2024withLastYear2025.pdf", sep="")),
#     width = 14, height = 10)
# 
# par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
# plot(-1000, xlim=c(0.5, nYears+1.5), ylim=c(0,1300),
#      xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
# axis(1, at=c(1:(length(SeasonText1))), labels = SeasonText1, cex.axis=1.2,padj = -1)
# at = c(1:nYears)
# abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))
# 
# widthPolygon <- 0.15
# widthPolygon1 <- 0.15
# widthPolygon2 <- 0.15
# offsetstar <- 0.05
# cexStar <- 1.5
# displayQuantiles50 <- TRUE
# 
# quantile95Tot <- quantile50Tot <- list()
# quantile95Swe <- quantile50Swe <- list()
# quantile95Nor <- quantile50Nor <- list()
# 
# for(t in nYears){
#   tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
#   quantile95Tot[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Tot[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
#   tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
#   quantile95Swe[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Swe[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
#   tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
#   quantile95Nor[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Nor[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
# }
# 
# quantile95Tot2024[[nYears+1]] <- quantile95Tot[[nYears]]
# quantile50Tot2024[[nYears+1]] <- quantile50Tot[[nYears]]
# quantile95Swe2024[[nYears+1]] <- quantile95Swe[[nYears]]
# quantile50Swe2024[[nYears+1]] <- quantile50Swe[[nYears]]
# quantile95Nor2024[[nYears+1]] <- quantile95Nor[[nYears]]
# quantile50Nor2024[[nYears+1]] <- quantile50Nor[[nYears]]
# 
# 
# 
# 
# for(t in 1:(nYears+1)){
#   #TOTAL
#   polygon(x = c(t - widthPolygon, t + widthPolygon,
#                 t + widthPolygon, t - widthPolygon ),
#           y = c(quantile95Tot2024[[t]][1], quantile95Tot2024[[t]][1],
#                 quantile95Tot2024[[t]][2], quantile95Tot2024[[t]][2]), 
#           col=adjustcolor(TotalColors, violin.alpha95),
#           border= NA)
#   
#   # add a star
#   if(sum(t %in% yearsNotSampled)){
#     text(x=t+widthPolygon+offsetstar ,y= quantile95Tot2024[[t]][2], "*",cex=cexStar)
#     
#   }
#   
#   if(displayQuantiles50){
#     polygon(x = c(t-widthPolygon, t+widthPolygon,
#                   t+widthPolygon, t-widthPolygon ),
#             y = c(quantile50Tot2024[[t]][1], quantile50Tot2024[[t]][1],
#                   quantile50Tot2024[[t]][2], quantile50Tot2024[[t]][2]), 
#             col=adjustcolor(TotalColors, violin.alpha50),
#             border= NA)
#   }
#   #SWEDEN
#   polygon(x = c(t, t -widthPolygon1*2,
#                 t - widthPolygon1*2, t),
#           y = c(quantile95Swe2024[[t]] [1], quantile95Swe2024[[t]] [1],
#                 quantile95Swe2024[[t]] [2], quantile95Swe2024[[t]] [2]), 
#           col=adjustcolor(country.colors[2], violin.alpha95),
#           border= NA)
#   
#   if(sum(t %in% yearsNotSampled)){
#     text(x= t+offsetstar ,y= quantile95Swe2024[[t]][2], "*",cex=cexStar)
#     
#   }
#   if(displayQuantiles50){
#     polygon(x = c(t, t-widthPolygon1*2,
#                   t-widthPolygon1*2, t),
#             y = c(quantile50Swe2024[[t]][1], quantile50Swe2024[[t]][1],
#                   quantile50Swe2024[[t]][2], quantile50Swe2024[[t]][2]), 
#             col=adjustcolor(country.colors[2], violin.alpha50),
#             border= NA)
#   }
#   
#   #NORWAY
#   polygon(x = c(t, t+widthPolygon1*2,
#                 t+widthPolygon1*2, t),
#           y = c(quantile95Nor2024[[t]][1], quantile95Nor2024[[t]][1],
#                 quantile95Nor2024[[t]][2], quantile95Nor2024[[t]][2]), 
#           col=adjustcolor(country.colors[1], violin.alpha95),
#           border= NA)
#   
#   if(displayQuantiles50){
#     polygon(x = c(t, t+widthPolygon2*2,
#                   t+widthPolygon2*2, t),
#             y = c(quantile50Nor2024[[t]][1], quantile50Nor2024[[t]][1],
#                   quantile50Nor2024[[t]][2], quantile50Nor2024[[t]][2]), 
#             col=adjustcolor(country.colors[1], violin.alpha50),
#             border= NA)
#   }
#   
#   
# }
# box()
# abline(v=at[1:(nYears)]+0.5,lty=2)
# 
# #legend
# par(xpd=TRUE)
# polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")
# 
# 
# labels <- c(" Norway  ", " Sweden  ", " Total")
# pch <- rep(19,4)
# cex <- rep(4,4)
# y <-  c(200, 200, 200)
# x <- c(1,3,5)
# mycol1 <- c(country.colors, "black")
# # add transparent background polygon
# #polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
# for(i in 1:3){
#   #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
#   points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
#   points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
#   text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
# }
# 
# 
# 
# dev.off()
# 
# 
## ------      2.5.1.5 PLOT 2024 and 2025 FOR COMPARISONS ------

# pdf(file= file.path(WDFigures,
#                     paste("NCountriesBarsSeason2024_2025Comp.pdf", sep="")),
#     width = 14, height = 10)
# 
# par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
# plot(-1000, xlim=c(0.5, nYears+1.5), ylim=c(0,1300),
#      xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
# axis(1, at=c(1:(length(SeasonText1))), labels = SeasonText1, cex.axis=1.2,padj = -1)
# at = c(1:nYears)
# abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))
# 
# for(t in 1:nYears){
#   tmp <- colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])
#   quantile95Tot[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Tot[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
#   tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",]
#   quantile95Swe[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Swe[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
#   tmp <- DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
#   quantile95Nor[[t]] <- quantile(tmp, prob=c(0.0275, 0.975))
#   quantile50Nor[[t]] <- quantile(tmp, prob=c(0.25, 0.75))
# }
# 
# quantile95Tot2024[[nYears+1]] <- quantile95Tot[[nYears]]
# quantile50Tot2024[[nYears+1]] <- quantile50Tot[[nYears]]
# quantile95Swe2024[[nYears+1]] <- quantile95Swe[[nYears]]
# quantile50Swe2024[[nYears+1]] <- quantile50Swe[[nYears]]
# quantile95Nor2024[[nYears+1]] <- quantile95Nor[[nYears]]
# quantile50Nor2024[[nYears+1]] <- quantile50Nor[[nYears]]
# 
# for(t in 1:(nYears)){
#   #TOTAL
#   polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                 t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#           y = c(quantile95Tot2024[[t]][1], quantile95Tot2024[[t]][1],
#                 quantile95Tot2024[[t]][2], quantile95Tot2024[[t]][2]), 
#           col=adjustcolor(TotalColors, violin.alpha95),
#           border= NA)
#   
#   # add a star
#   # if(sum(t %in% yearsNotSampled)){
#   #   text(x=t+widthPolygon+offsetstar ,y= quantile95Tot2024[[t]][2], "*",cex=cexStar)
#   #   
#   # }
#   
#   if(displayQuantiles50){
#     polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                   t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#             y = c(quantile50Tot2024[[t]][1], quantile50Tot2024[[t]][1],
#                   quantile50Tot2024[[t]][2], quantile50Tot2024[[t]][2]), 
#             col=adjustcolor(TotalColors, violin.alpha50),
#             border= NA)
#   }
#   #SWEDEN
#   polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                 t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#           y = c(quantile95Swe2024[[t]] [1], quantile95Swe2024[[t]] [1],
#                 quantile95Swe2024[[t]] [2], quantile95Swe2024[[t]] [2]), 
#           col=adjustcolor(country.colors[2], violin.alpha95),
#           border= NA)
#   
#   # if(sum(t %in% yearsNotSampled)){
#   #   text(x= t+offsetstar ,y= quantile95Swe2024[[t]][2], "*",cex=cexStar)
#   #   
#   # }
#   if(displayQuantiles50){
#     polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                   t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#             y = c(quantile50Swe2024[[t]][1], quantile50Swe2024[[t]][1],
#                   quantile50Swe2024[[t]][2], quantile50Swe2024[[t]][2]), 
#             col=adjustcolor(country.colors[2], violin.alpha50),
#             border= NA)
#   }
#   
#   #NORWAY
#   polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                 t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#           y = c(quantile95Nor2024[[t]][1], quantile95Nor2024[[t]][1],
#                 quantile95Nor2024[[t]][2], quantile95Nor2024[[t]][2]), 
#           col=adjustcolor(country.colors[1], violin.alpha95),
#           border= NA)
#   
#   if(displayQuantiles50){
#     polygon(x = c(t - widthPolygon-0.20, t + widthPolygon-0.20,
#                   t + widthPolygon-0.20, t - widthPolygon-0.20 ),
#             y = c(quantile50Nor2024[[t]][1], quantile50Nor2024[[t]][1],
#                   quantile50Nor2024[[t]][2], quantile50Nor2024[[t]][2]), 
#             col=adjustcolor(country.colors[1], violin.alpha50),
#             border= NA)
#   }
#   
#   
# }
# 
# ###
# 
# for(t in 1:(nYears)){
#   #TOTAL
#   polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                 t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#           y = c(quantile95Tot[[t]][1], quantile95Tot[[t]][1],
#                 quantile95Tot[[t]][2], quantile95Tot[[t]][2]), 
#           col=adjustcolor(TotalColors, violin.alpha95),
#           border= NA)
#   
#   # add a star
#   # if(sum(t %in% yearsNotSampled)){
#   #   text(x=t+widthPolygon+offsetstar ,y= quantile95Tot2024[[t]][2], "*",cex=cexStar)
#   #   
#   # }
#   
#   if(displayQuantiles50){
#     polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                   t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#             y = c(quantile50Tot[[t]][1], quantile50Tot[[t]][1],
#                   quantile50Tot[[t]][2], quantile50Tot[[t]][2]), 
#             col=adjustcolor(TotalColors, violin.alpha50),
#             border= NA)
#   }
#   #SWEDEN
#   polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                 t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#           y = c(quantile95Swe[[t]] [1], quantile95Swe[[t]] [1],
#                 quantile95Swe[[t]] [2], quantile95Swe[[t]] [2]), 
#           col=adjustcolor(country.colors[2], violin.alpha95),
#           border= NA)
#   
#   # if(sum(t %in% yearsNotSampled)){
#   #   text(x= t+offsetstar ,y= quantile95Swe2024[[t]][2], "*",cex=cexStar)
#   #   
#   # }
#   if(displayQuantiles50){
#     polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                   t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#             y = c(quantile50Swe[[t]][1], quantile50Swe[[t]][1],
#                   quantile50Swe[[t]][2], quantile50Swe[[t]][2]), 
#             col=adjustcolor(country.colors[2], violin.alpha50),
#             border= NA)
#   }
#   
#   #NORWAY
#   polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                 t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#           y = c(quantile95Nor[[t]][1], quantile95Nor[[t]][1],
#                 quantile95Nor[[t]][2], quantile95Nor[[t]][2]), 
#           col=adjustcolor(country.colors[1], violin.alpha95),
#           border= NA)
#   
#   if(displayQuantiles50){
#     polygon(x = c(t+1 - widthPolygon+0.20, t+1 + widthPolygon+0.20,
#                   t+1 + widthPolygon+0.20, t+1 - widthPolygon+0.20 ),
#             y = c(quantile50Nor[[t]][1], quantile50Nor[[t]][1],
#                   quantile50Nor[[t]][2], quantile50Nor[[t]][2]), 
#             col=adjustcolor(country.colors[1], violin.alpha50),
#             border= NA)
#   }
#   
#   
# }
# 
# box()
# abline(v=at[1:(nYears)]+0.5,lty=2)
# 
# #legend
# par(xpd=TRUE)
# polygon(x=c(0.8,7.2,7.2,0.8),y=c(170,170,230,230),col=adjustcolor("white",alpha.f = 0.9),border="white")
# 
# 
# labels <- c(" Norway  ", " Sweden  ", " Total")
# pch <- rep(19,4)
# cex <- rep(4,4)
# y <-  c(200, 200, 200)
# x <- c(1,3,5)
# mycol1 <- c(country.colors, "black")
# # add transparent background polygon
# #polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
# for(i in 1:3){
#   #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
#   points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
#   points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
#   text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
# }
# 
# 
# 
# dev.off()





## ------    2.5.3 MAPS ------

habbdensCropped <- list()
max <- max(unlist(lapply(DensityCountriesRegions, function(x) max(x$MeanCell))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

#PLOT
pdf(file=file.path(WDFigures, paste("DensityMapsAC5kms.pdf",sep="")))
for(t in 1: nYears){
  habbdens <- densityInputRegions$regions.r
  habbdens[] <- NA
  habbdens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  habbdensCropped[[t]] <- habbdens#crop(habbdens, e.sp)
  
  plot(habbdensCropped[[t]], breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
  plot(nngeo::st_remove_holes(myHabitat.list$habitat.poly),add=T, col=NA,border=grey(0.5))
  # points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],],
  #        pch=16, cex=0.4, col=adjustcolor("black",alpha.f = 0.2))
  plot(habbdensCropped[[t]], legend.only=TRUE,breaks=cuts, col=col,
       legend.width = 2,
       axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
                      labels=round(seq(0, max, length.out = 5),digits = 1),
                      cex.axis=0.6),
       legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
  
}
dev.off()




## ------  3. UD BASED DENSITY (5km) ------

### IDENTIFY PROXIMITY HABITAT CELLS 
habitatMask <- densityInputCountries$habitat.id
habitatMask[!is.na(habitatMask)] <- 1 
# DetIndex <- getLocalObjects(habitatMask =  habitatMask,
#                             coords = densityInputCountries$habitat.xy,
#                             dmax = 15,
#                             resizeFactor = 1)

### COMPUTE THE UD BASED FOR A FEW ITERATIONS
#RESCALE SIGMA TO METERS 
sigma <- myResultsSXYZ_MF$sims.list$sigma#*res(myHabitat.list$habitat.r)[1]
#RESCALE SIGMA TO THE HABITAT SCALE
sigmaRescaled <- sigma/res(rrCountries)[1]

#SELECT xxx ITERATIONS RANDOMLY 
spaceUSED <- list()
iter <- sample(1:dim(densityInputCountries$sy)[1], size = 1000)#dim(densityInputCountries$sx)[1])
for(t in 1:nYears){
  ## spaceUSED[[t]] <- GetSpaceUseLESS( densityInputCountries$sx[iter,,t],
  ##                               densityInputCountries$sy[iter,,t],
  ##                               myResultsSXYZ_MF$sims.list$z[iter,,t],
  ##                               sigmaRescaled[iter],
  ##                               densityInputCountries$habitat.xy,
  ##                               aliveStates = alive.states,
  ##                               regionID = regionID,
  ##                               habitatID = DetIndex$habitatGrid-1,
  ##                               habitatIndex = DetIndex$localIndices-1,
  ##                               nHabitatLESS = DetIndex$numLocalIndices,
  ##                               display_progress = T,
  ##                               returnPosteriorCells = F
  ## )
  gc()
  # spaceUSED[[t]] <- GetSpaceUse(densityInputCountries$sx[iter,,t],
  #                               densityInputCountries$sy[iter,,t],
  #                               myResultsSXYZ_MF$sims.list$z[iter,,t],
  #                               sigmaRescaled[iter],#sigmaRescaled[iter],
  #                               densityInputCountries$habitat.xy,
  #                               aliveStates = alive.states,
  #                               regionID = densityInputCountries$regions.rgmx,
  #                               display_progress = T,
  #                               returnPosteriorCells = T
  # )
}

t=9
plot(densityInputCountries$habitat.xy[,2]~
       densityInputCountries$habitat.xy[,1],pch=16,cex=0.1)

for(iter in 1:100){
  points(densityInputCountries$sy[iter,myResultsSXYZ_MF$sims.list$z[iter,,t] %in%2,t]~
           densityInputCountries$sx[iter,myResultsSXYZ_MF$sims.list$z[iter,,t] %in%2,t],col="red",pch=16,cex=0.1)
}





plot(densityInputCountries$sy[1,myResultsSXYZ_MF$sims.list$z[1,,t] %in%2,t]~
       densityInputCountries$sx[1,myResultsSXYZ_MF$sims.list$z[1,,t] %in%2,t])

spaceUSED1 <- spaceUSED
spaceUSED <- list()
for(t in 1:nYears){
  spaceUSED[[t]] <- list()
  spaceUSED[[t]][["MeanCell"]] <- spaceUSED1[[t]]$MeanCell
}
# save(spaceUSED, file = file.path(WDFigures, "spaceUsed5km.RData" ))
load(file = file.path(WDFigures, "spaceUsed5km.RData" ))

## ------          2.1.2.1 PLOT TIME SERIES ------
#### PLOT TIME SERIES 
## PREPARE THE FILES 
SeasonText <- lapply(YEARS,FUN = function(x) paste(x,collapse = "/"))
#SeasonText <- lapply(YEARS,FUN = function(x) paste(x[[2]]))
# habbdensUDCropped[[t]]
COUNTRIESSCA <- COUNTRIES[COUNTRIES$ISO %in% c("NOR","SWE"),]
COUNTRIESsimpFig <- st_simplify(COUNTRIESSCA, preserveTopology = F,dTolerance = 4000)
habbdensFig <- densityInputRegions$regions.r
habbRFig <- densityInputRegions$regions.r
## from 25km2 (5*5raster) to 100km2
spaceUSED100km2 <- lapply(spaceUSED, function(x) x$MeanCell * 4 )


# studDis <- disaggregate(myHabitat.list$habitat.poly)
# studDis$id <- 1:length(studDis)
# plot(studDis)
# text(studDis,studDis$id)
# lakes <- disaggregate(COUNTRIESWaterHumans[COUNTRIESWaterHumans$area>20000000000,])
#lakes <- dropHole(lakes)
# plot(lakes)

#find the lakes 
# which(unlist(lapply(lakes@polygons[[2]]@Polygons, function(x) x@area))>500000000)
# featureNumber=2 ; ringNumber=115
# Lake1 = SpatialPolygons(
#   list(
#     Polygons(
#       list(
#         lakes@polygons[[featureNumber]]@Polygons[[ringNumber]]
#       ),
#       ID=1)))
# featureNumber=2 ; ringNumber=103
# Lake2 = SpatialPolygons(
#   list(
#     Polygons(
#       list(
#         lakes@polygons[[featureNumber]]@Polygons[[ringNumber]]
#       ),
#       ID=1)))

##PLOT
pdf(file=file.path(WDFigures, paste("DensityMapsUD.pdf",sep="")), width = 12, height = 8)
#layout
mx <- rbind(c(1,rep(1:5, each=2)),
            c(rep(1:5, each=2),5))
mx <- rbind(mx, mx+5)
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)
habbdensUDCropped <- list()
for(t in 1:length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  habbdensFig[!is.na(habbRFig[])] <- spaceUSED100km2[[t]]
  habbdensFig[habbRFig[]==0] <- NA
  
  habbdensUDCropped[[t]] <- habbdensFig#mask(habbdensFig, e.sp)
  crs( habbdensUDCropped[[t]]) <- st_crs(myHabitat.list$habitat.poly)
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image(habbdensUDCropped[[t]], add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  # plot(RemoveHolesSp(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  
  mtext(SeasonText[[t]], 1, -2, adj=0.15, cex=1.2)
  #box()
  # spPol <- rasterToPolygons(myHabitat.list$habitat.r,
  #                  fun = function(x) x==1) 
  # spPol <- aggregate(spPol)
  # plot(SpatialPolygons(spPol@polygons[[1]]@Polygons[[1]]))
  # 
  # spPol <- st_as_sf(spPol)
  # state_union <- spPol %>% 
  #   group_by(Habitat) %>%
  #   summarise(geometry = sf::st_union(geometry)) %>%
  #   ungroup() %>% st_as_sf()
  # plot(state_union$geometry)
  # 
  # agg <- aggregate(rasterToPolygons(myHabitat.list$habitat.r,
  #                                   fun = function(x) x==1))
  # agg <- RemoveHolesSp(agg)
  # 
  # plot(agg, add=TRUE, border="black", col=NA)
  agg1 <- aggregate(rasterToPolygons(myHabitat.list$habitat.rWthBuffer,
                                     fun = function(x) x==1))

  if(sum(t%in%yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1320000,x1=1320000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=2, lend=2)  
    text(1280000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    plot(habbdensUDCropped[[t]], legend.only=TRUE,breaks=cuts, col=col,
         legend.width = 2,
         axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        cex.axis=1.6),
         smallplot=c(0.72, 0.75, 0.2, 0.4),
         legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                          side=4, font=1, line=4, cex=1.2))
  }
  
}
dev.off()

## ------          2.1.2.3 PLOT LAST YEAR ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYear.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)#YEARS[[t]][2]

  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 0.5,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.73, 0.75, 0.25, 0.4),
          legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()


## ------          2.1.2.4 PLOT LAST YEAR SUMMARY ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYearSummary.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  # plot(COUNTRIESsimpFig[1], border=grey(0.1), col = NA, add=TRUE, lwd=2)
  
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)
  
  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.72, 0.75, 0.2, 0.4),
          legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()



## ------          2.1.2.5 PLOT LAST YEAR SUMMARY NO ------
## PLOT LAST  YEAR  
pdf(file=file.path(WDFigures, paste("DensityMapsUDLastYearSummaryNO.pdf",sep="")), 
    width = 8, height = 8)

#layout
mx <- rbind(c(1,rep(1, each=)),
            c(rep(1, each=), 2))
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(spaceUSED100km2, function(x) max(x))))
cuts <- seq(0, max, length.out = 100)   #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

for(t in length(years)){
  par(mar=c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  
  plot(st_geometry(COUNTRIESsimpFig), border=NA,col = grey(0.85))
  
  # habbdensFig[habbRFig[]>0] <- spaceUSED100km2[[t]]
  # habbdensCropped <- mask(habbdensFig, e.sp)
  
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image( habbdensUDCropped[[t]] , add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  # plot(COUNTRIESsimpFig[1], border=grey(0.1), col = NA, add=TRUE, lwd=2)
  
  mtext(SeasonText[[t]], 1, -4, adj=0.25, cex=1.2)
  
  #PLOT COUNTIES 
  if(sum(t %in% yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1130000,x1=1130000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=4, lend=2)  
    text(1100000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    length(years)
    
    plot( habbdensUDCropped[[t]] , legend.only=TRUE,breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                         cex.axis=1.6),
          smallplot=c(0.72, 0.75, 0.2, 0.4),
          legend.args=list(text=expression(paste("Individer/100 km"^ 2, "", sep="")),
                           side=4, font=1, line=4.5, cex=1.2))
    
    
  }
  
  
  
  
}
dev.off()




## ------          2.1.2.6 WRITE UD 5km RASTER FOR ROVBASE ------

if(!dir.exists(file.path(WDFigures, "RasterForRovbase"))){dir.create(file.path(WDFigures, "RasterForRovbase"))}

for(t in 1:length(years)){
  raster::crs(habbdensUDCropped[[t]]) <- "EPSG:32633"#st_crs(myHabitat.list$habitat.poly))
  path <- file.path(WDFigures, "RasterForRovbase",paste("wolverine_5km",paste(YEARS[[t]][1],collapse = "_"),".tif",sep=""))
  writeRaster(habbdensUDCropped[[t]], path, overwrite=TRUE)
}



## ------  4. DERIVED PARAMETERS FROM ABUNDANCE ------ 

## ------    4.1 MAKE A GROWTH RATE TABLE PER COUNTRY  ------

growthRate <- matrix(0, ncol=nYears-1,nrow=3)
row.names(growthRate) <- c("Norway","Sweden","Total")
colnames(growthRate) <- unlist(lapply(YEARS[2:(length(YEARS))],function(x) paste(x,collapse =  "-")))

for(t in 1:(nYears-1)){
  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Norway",] /
    DensityCountriesRegions[[t]]$PosteriorRegions["Norway",]
  growthRate["Norway",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                  " (",
                                  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                  ")",sep="")
  
  
  
  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Sweden",]/ 
    (DensityCountriesRegions[[t]]$PosteriorRegions["Sweden",])
  growthRate["Sweden",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                  " (",
                                  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                  ")",sep="")
  
  
  growth <- colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"),])/
    colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])#colSums(DensityCountriesRegions[[t+1]]$PosteriorAllRegions) / 
  #colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)
  growthRate["Total",t] <- paste(format(round(mean(growth),digits = 2),nsmall = 2),
                                 " (",format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                 format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),
                                 ")",sep="")
  
}


### add *** to years and country with OPSCR results 
#row.names(growthRate)[2] <- "Sweden*"
#row.names(growthRate)[3] <- "Total*"
# transitionNotSampled <- colnames(growthRate)
# transitionNotSampled <- unique(unlist(lapply(1:length(years),function(x) grep(as.character((years+1)[yearsNotSampled])[x], transitionNotSampled))))
# transitionNotSampled <- transitionNotSampled[!is.na(transitionNotSampled)] 
# growthRate[2,transitionNotSampled] <- paste(growthRate[2,transitionNotSampled],"*",sep="")
# growthRate[3,transitionNotSampled] <- paste(growthRate[3,transitionNotSampled],"*",sep="")


#print table 
print(xtable(growthRate, type = "latex", align=paste(c("l", rep("c",ncol(growthRate))), collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(growthRate),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("growthRate.tex", sep="")))


colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"),])
## ------    4.2 DERIVE SEX RATIO  ------

PropFemale <- PropFemaleSWE <- PropFemaleNOR <- list()

for(t in 1:nYears){
  PropFemale[[t]] <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])/
    (colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]) + colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]))
  
  PropFemaleSWE[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",])
  
  PropFemaleNOR[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Norway",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",])
}

#OVERALL PROPORTION OF FEMALES
mean(unlist(PropFemale))
round(quantile(unlist(PropFemale), probs=c(0.025,0.975)),digits = 2)

#median(PropFemale[[nYears]])
round(quantile(PropFemale[[nYears]], probs=c(0.025,0.975)),digits=2)
round(mean(PropFemale[[nYears]]),digits=2)

lapply(PropFemale,function(x) round(quantile(x,probs=c(0.025,0.975)),digits=2) )
# mean(do.call(c,PropFemale))
# quantile(do.call(c,PropFemale), probs=c(0.025,0.975))

propFemale_tab <- matrix(0, ncol=nYears,nrow=3)
row.names(propFemale_tab) <- c("Norway","Sweden","Total")
colnames(propFemale_tab) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/")))# unlist(lapply(YEARS,function(x) c(x[2])))#c(unlist(lapply(YEARS[1:(length(YEARS))],function(x) paste(x,collapse =  "-"))))
for(t in 1:nYears){
  #Sweden
  propFemale_tab["Sweden",t] <- paste(format(round(mean(PropFemaleSWE[[t]]),digits = 2),nsmall = 2),
                                      " (",
                                      format(round(quantile(PropFemaleSWE[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                      format(round(quantile(PropFemaleSWE[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                      ")",sep="")
  #NORWAY
  propFemale_tab["Norway",t] <- paste(format(round(mean(PropFemaleNOR[[t]]),digits = 2),nsmall = 2),
                                      " (",
                                      format(round(quantile(PropFemaleNOR[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                      format(round(quantile(PropFemaleNOR[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                      ")",sep="")
  propFemale_tab["Total",t] <- paste(format(round(mean(PropFemale[[t]]),digits = 2),nsmall = 2),
                                     " (",
                                     format(round(quantile(PropFemale[[t]], probs=c(0.025)), digits = 2),nsmall = 2),"-",
                                     format(round(quantile(PropFemale[[t]], probs=c(0.975)), digits = 2),nsmall = 2),
                                     ")",sep="")
}

#print table 
print(xtable(propFemale_tab, type = "latex", align=paste(c("l", rep("c",ncol(propFemale_tab))), collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(propFemale_tab),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("propFemale.tex", sep="")))



## ------    4.3 DERIVE DENSITY  ------

habbRCarRegionsTRY <- rrCountries
habbRCarRegionsTRY[!is.na(habbRCarRegionsTRY[])] <-1 
habbRCarRegionsTRY[] <- as.numeric(habbRCarRegionsTRY[])

habbRCarRegionsTRYpol <- sf::st_as_sf(stars::st_as_stars(habbRCarRegionsTRY), 
                                      as_points = FALSE, merge = F)

plot(rasterToPolygons(habbRCarRegionsTRY, function(x) x>0,dissolve = T))
areaSqKm <- sum(st_area(habbRCarRegionsTRYpol))#*1e-6
units::set_units(areaSqKm, km^2)
#size of the area in km2
areaSqKm

#multiplied by 100 to get per 100km2
DensityCountriesRegions[[t]]$summary["Total","mean"]/areaSqKm*100
DensityCountriesRegions[[t]]$summary["Total","95%CILow"]/areaSqKm*100
DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]/areaSqKm*100

###COULD REMAKE ALL THE TABLES WITH DENSITY INSTEAD OF ABUNDANCE....


## ------    4.4  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED ------

n.detectedTotal <- read.csv(file.path(WDTables, paste("TotalIdDetected.csv", sep="")))
n.detectedTotal <- as.vector(n.detectedTotal[1, 2:ncol(n.detectedTotal)])

n.detectedSex <- read.csv(file.path(WDTables, paste("NGSidCountrySEX.csv", sep="")))
tF <- seq(2, length.out = nYears, by=2)
tM <- seq(3, length.out = nYears, by=2)

propDetected <- matrix("", ncol=nYears,nrow=3)
row.names(propDetected) <- c("M","F","Total")

colnames(propDetected) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/"))) #unlist(lapply(YEARS,function(x) c(x[2])))#unlist(lapply(YEARS,function(x) paste(x,collapse = "/")))#

for(t in 1:nYears){
  propDetected["Total",t] <-  paste(format(round(mean(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                    " (",format(round(quantile(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                    format(round(quantile(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
  
  n.detectedSexM <- as.numeric(as.character(n.detectedSex[4,tM[t]]))
  propDetected["M",t] <-  paste(format(round(mean(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                " (",format(round(quantile(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.025)), digits = 2),nsmall = 2), "-",
                                format(round(quantile(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.975)), digits = 2),nsmall = 2),")", sep="")
  
  n.detectedSexF <- as.numeric(as.character(n.detectedSex[4,tF[t]]))
  propDetected["F",t] <-  paste(format(round(mean(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])), digits = 2),nsmall = 2),
                                " (",format(round(quantile(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.025)), digits = 2),nsmall = 2), "-",
                                format(round(quantile(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]), probs=c(0.975)), digits = 2),nsmall = 2),")", sep="")
}

##PRINT
print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row=list(list(seq(1,nrow(propDetected), by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(WDTables, paste("PropDetectedIds.tex",sep="")))
gc()


## ------    4.5  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED PER COUNTRIES ------

n.detectedCountry <- read.csv(file.path(WDTables, paste("NGSidCountrySEX.csv", sep="")))
paste(x,collapse = "/")
colnames(n.detectedCountry) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/"))))) 
#propDetected <- matrix("", ncol=nYears,nrow=3)
#row.names(propDetected) <- c("M","F","Total")
#colnames(propDetected) <- unlist(lapply(YEARS,function(x) c(x[2])))#unlist(lapply(YEARS,function(x) paste(x,collapse = "/")))#
propDetectedCountry <- n.detectedCountry 
propDetectedCountry[2:4, 2:ncol(propDetectedCountry)] <- NA

NCountrySex <- read.csv(file.path(WDTables, paste("NAllYearsPerSex.csv", sep="")))
colnames(NCountrySex) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))#c(x[2],x[2],x[2]))))

yrs <- unlist(lapply(YEARS, function(x) c(paste(x,collapse = "/"))))#c(x[2]))) 


lisCountries <- list()
lisCountries[[1]] <- c("Norway")
lisCountries[[2]] <- c("Sweden")
lisCountries[[3]] <- c("Sweden","Norway")

for(t in 1:nYears){
  col <- which(colnames(NCountrySex)  %in% yrs[t])
  NCountrySex[,col]
  
  cols <- which(colnames(propDetectedCountry)  %in% as.character(yrs[t]))
  for(p in 1:2){
    propDetectedCountry[p+1,cols[1]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[1]]) /DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[p]],]), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
    propDetectedCountry[p+1,cols[2]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[2]]) /DensityCountriesRegionsM[[t]]$PosteriorRegions[lisCountries[[p]],]), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),],probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
  }
  for(p in 3:3){
    propDetectedCountry[p+1,cols[1]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[1]]) /colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[p]],])), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[1]])/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
    propDetectedCountry[p+1,cols[2]] <-  paste(format(round(mean(as.numeric(n.detectedCountry[p+1,cols[2]]) /colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),])), digits = 2),nsmall = 2),
                                               " (",format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.025)),digits = 2),nsmall = 2),"-",
                                               format(round(quantile(as.numeric(n.detectedCountry[p+1,cols[2]])/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c(lisCountries[[p]]),]),probs=c(0.975)) ,digits = 2),nsmall = 2),")",sep="")
    
  }
  
  
  
}



addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(propDetectedCountry)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                    '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))

print(xtable(propDetectedCountry, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry)+1), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry.tex", sep="")))

##SPLIT THE TABLE IN TWO 
propDetectedCountry1 <- propDetectedCountry[,c(1:11)]
propDetectedCountry2 <- propDetectedCountry[,c(1,12:21)]

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[2:11])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[12:21])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))


#SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(propDetectedCountry1, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(propDetectedCountry2, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("propDetectedCountry2.tex", sep="")))



## ------  5. VITAL RATES  ------

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15



## ------    5.1 SURVIVAL BARS  ------

pdf(file=file.path(WDFigures, paste("SurvivalBars.pdf",sep="")),width=10,height=6)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears-1+0.5), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Survival")
axis(2,tck=-0.02)
abline(v=1:(nYears-1)+0.5,lty=2)
SeasonTextvec <- unlist(lapply(YEARS,FUN = function(x) paste(substring(x, nchar(x) - 1),collapse = "/")))

#axis(1, c(1:nYears), labels = paste(years,years+1,sep=" to\n "),  cex.axis=1.1,padj  = 0.5)
#axis(1, c(1:nYears), labels = paste(years+1,years+2,sep=" to\n "),cex.axis=0.5,padj = -2)
axis(1, c(1:(nYears-1)),labels = paste(SeasonTextvec[1:(nYears-1)], SeasonTextvec[2:(nYears)],sep=" to\n"),
     cex.axis=1.2,padj = 0.2,tick = F)

myCol <- c("#E69F00","#009E73")

#myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
widthPolygon <- 0.15
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  
  for(t in 1:(nYears-1)){
    
    quantile95 <- quantile(myResults$sims.list$phi[,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$phi[,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}#i

#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)

labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
# add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------    5.2 MORTALITY BARS  ------

## DERIVE MORTALITY FROM POSTERIOR AND NUMBER OF DEAD RECOVERIES 
IDF <- row.names(nimDataF$z)
IDM <- row.names(nimDataM$z)

IDMF <- c(IDM, IDF)
load(file.path(dir.dropbox,"wolverine/CM/2025/54.Cleaned_DR_2015-2025/Hann/54.Cleaned_DR_2015-2025Hann_Chain1.RData"))
y.deadM <- y.dead
load(file.path(dir.dropbox,"wolverine/CM/2025/54.Cleaned_DR_2015-2025/Hunn/54.Cleaned_DR_2015-2025Hunn_Chain1.RData"))
y.deadF <- y.dead
y.deadMF <- rbind(y.deadM, y.deadF)

## ARRAY 
itera <- seq(1,dim(myResultsSXYZ_MF$sims.list$z1)[1],by=10)# 1:dim(myResultsSXYZ_MF$sims.list$z1)[1]#sample(1:dim(myResultsSXYZ_MF$sims.list$z1)[1], size = 1000)

MortalityAll  <- array(NA, c(length(itera),nYears-1,2))
MortalityCulled <- array(NA, c(length(itera),nYears-1,2))
MortalityOther <- array(NA, c(length(itera),nYears-1,2))

dimnames(MortalityAll)[[3]] <- dimnames(MortalityCulled)[[3]] <- dimnames(MortalityOther)[[3]] <- c("M","F")

for(iter in 1:length(itera) ){
  for(t in 1:(nYears-1)){
    #MALE
    sumDead <-  sum(myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2 &
                      myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "M",t+1] %in% 3) 
    sumAlive <- sum(myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2)
    
    MortalityAll[iter,t,1] <- sumDead / sumAlive
    
    MortalityOther[iter,t,1] <- (sumDead  -  sum(y.deadM[,t+1])) /sumAlive
    MortalityCulled[iter,t,1] <- MortalityAll[iter,t,1]- MortalityOther[iter,t,1] 
    
    #FEMALE
    sumDead <-  sum(myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "F",t] %in% 2 &
                      myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "F",t+1] %in% 3) 
    sumAlive <- sum(myResultsSXYZ_MF$sims.list$z1[itera[iter], myResultsSXYZ_MF$sims.list$sex %in% "F",t] %in% 2)
    
    MortalityAll[iter,t,2] <- sumDead / sumAlive
    
    MortalityOther[iter,t,2] <- (sumDead  -  sum(y.deadF[,t+1])) /sumAlive
    MortalityCulled[iter,t,2] <- MortalityAll[iter,t,2]- MortalityOther[iter,t,2] 
  }
}

gc()

Results.list[["M"]]$sims.list$w <- MortalityOther[,,1]
Results.list[["M"]]$sims.list$h<- MortalityCulled[,,1]

Results.list[["F"]]$sims.list$w <- MortalityOther[,,2]
Results.list[["F"]]$sims.list$h <- MortalityCulled[,,2]


pdf(file=file.path(WDFigures, paste("MortalityBars.pdf",sep="")),width=10,height=8)

nf <- layout(cbind(c(3,4,5),c(7,1,2),c(9,8,6)),widths=c(0.15,1,0.35),heights=c(0.15,1,1))
# layout.show(nf)
for(i in c("F","M")){
  
  par(mar=c(4,4.5,0.5,1),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
  plot(10, xlim = c(0.5, nYears-1+0.5), ylim = c(0,0.8), type ="n", xaxt="n", xlab = "Years", ylab = "Mortality")
  axis(2,tck=-0.02)
  abline(v=1:(nYears-1)+0.5,lty=2)
  
  axis(1, c(1:nYears), labels = paste(years+1,years+2,sep=" to\n "),  cex.axis=1.1,padj  = 0.5)
  myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
  myDev <- c(-0.15,+0.15)
  
  
  
  for(t in 1:(nYears-1)){
    #culled
    quantile95 <- quantile( MortalityCulled[,t,i], prob=c(0.0275, 0.975))
    quantile50 <- quantile( MortalityCulled[,t,i], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[1] - widthPolygon, t+myDev[1] + widthPolygon,
                  t+myDev[1] + widthPolygon, t+myDev[1] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[1], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[1]-widthPolygon, t+myDev[1]+widthPolygon,
                  t+myDev[1]+widthPolygon, t+myDev[1]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[1], violin.alpha50),
            border= NA)
    
    
    #other
    quantile95 <- quantile( MortalityOther[,t,i], prob=c(0.0275, 0.975))
    quantile50 <- quantile( MortalityOther[,t,i], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[2] - widthPolygon, t+myDev[2] + widthPolygon,
                  t+myDev[2] + widthPolygon, t+myDev[2] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[2], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[2]-widthPolygon, t+myDev[2]+widthPolygon,
                  t+myDev[2]+widthPolygon, t+myDev[2]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[2], violin.alpha50),
            border= NA)
    
  }
  t
}#i

par(mar=c(0,0,0,0))
plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")

my.labels=c("Female","Male")
for(i in my.labels){
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
  text(0,0,labels=i,srt=90,cex=1,font=2)
}

#----LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
cols <- myCol#c("orange","darkorange2")
pt.col <- myCol#c("white","white")
labels<-c("Legal\nculling","Other\nmortality")
for(i in 1:2){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(2,7-i*2,pch=15,cex=3.5,col=adjustcolor(cols[i],violin.alpha95))
  points(2,7-i*2,pch=15,cex=1.5,col=adjustcolor(pt.col[i],violin.alpha50))
  text(3.5,7-i*2,labels[i],cex=2,pos=4)
}
dev.off()


## ------    5.2 PER CAPITA RECRUITMENT  ------

## ------       5.2.1 PLOT NUMBER OF RECRUITS + PER CAPITA RECRUITMENT ------

pdf(file=file.path(WDFigures, paste("NbRecruitsPerCapita.pdf",sep="")),width=10,height=8)

nf <- layout(rbind(c(3,5,6,7),
                   c(3,1,2,4),
                   c(8,1,2,4)),
             widths=c(0.15,1,1,0.1), heights=c(0.15,0.5,0.5))
##PER CAPITA RECRUITMENT 
par(mar=c(4.5,4.5,1,1) ,xaxs="i", cex.axis=1.3, cex.lab=1.6)#
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,1.2), type ="n", xaxt="n", xlab = "Years", ylab = "Per-capita recruitment")
axis(2, tck=-0.02)
# axis(1, c(2:(nYears+1)), labels = paste(years, years+1, sep= " to "), 
#      cex.axis = 0.8)
SeasonTextvec <- unlist(lapply(YEARS,FUN = function(x) paste(substring(x, nchar(x) - 1),collapse = "/")))

#axis(1, c(1:nYears), labels = paste(years,years+1,sep=" to\n "),  cex.axis=1.1,padj  = 0.5)
#axis(1, c(1:nYears), labels = paste(years+1,years+2,sep=" to\n "),cex.axis=0.5,padj = -2)
axis(1, c(1:(nYears-1)),labels = paste(SeasonTextvec[1:(nYears-1)], SeasonTextvec[2:(nYears)],sep=" to\n"),
     cex.axis=1.2,padj = 0.2,tick = F)

abline(v=1:(nYears-1)+0.5,lty=2)
myCol <- grey(0.2)#c("red","blue")

myResults <- myResultsSXYZ_MF

for(t in 1:(nYears-1)){
  #available <-apply(myResults$sims.list$z[,,t+1],1,function(x)sum(x%in%c(1)))
  # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
  n.recruit <- apply(myResults$sims.list$z[itera,,c(t,t+1)], 1, function(x) sum( x[,1] %in% c(1) & x[,2] %in% c(2) ))
  # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
  alivetminus1 <- apply(myResults$sims.list$z[itera,,t], 1, function(x)sum(x %in% c(2,3)))
  temp <- n.recruit/alivetminus1
  plot.violins2(list(temp),
                x = t,
                at = t+1,
                violin.width = 0.15,
                col = myCol,
                add = T,
                alpha = 0.6,
                plot.ci=0.95,
                border.col = myCol,
                cex=1)
  gc()
}#t

##NUMBER OF RECRUITS 
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,270), type ="n", xaxt="n", xlab = "Years", ylab = "Number of recruits")
axis(1, c(1:(nYears)),labels = paste(years, years+1,sep=" to "), cex.axis=0.8)
axis(2,tck=-0.02)
abline(v=1:(nYears-1)+0.5,lty=2)

for(t in 1:(nYears-1)){
  n.recruit <- apply(myResults$sims.list$z[itera,,c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
  plot.violins2(list(n.recruit),
                x = t,
                at = t+1,
                violin.width = 0.15,
                col = myCol,
                add = T,
                alpha = 0.6,
                plot.ci=0.95,
                border.col = myCol,
                cex=1)
  gc()
}#t

par(mar=c(0,0,0,0))
plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
dev.off()



## ------       5.2.1 PLOT NUMBER OF RECRUITS  ------

pdf( file = file.path(WDFigures, "Recruitment.pdf"),
     width = 10, height = 6)
widthPolygon <- 0.15

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
             widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
myResults <- myResultsSXYZ_MF

##NUMBER OF RECRUITS 
plot(10, xlim = c(0.5, nYears-0.5), ylim = c(0,300), type ="n", xaxt="n", xlab = "Years", ylab = "Number of recruits")
SeasonTextvec <- unlist(lapply(YEARS,FUN = function(x) paste(substring(x, nchar(x) - 1),collapse = "/")))

axis(1, c(1:(nYears-1)),labels = paste(SeasonTextvec[1:(nYears-1)], SeasonTextvec[2:(nYears)],sep=" to\n"),
     cex.axis=1.2,padj = 0.2,tick = F)

axis(2,tck=-0.02)
abline(v=1:(nYears-1)+0.5,lty=2)
myCol <- c("#E69F00","#009E73")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

myDev <- c(-0.16,+0.16)

gc()
for(s in 1:2){
  if(s==1){
    IDSex <- which(myResultsSXYZ_MF$sims.list$sex=="F")
  }else{
    IDSex <- which(myResultsSXYZ_MF$sims.list$sex=="M")
  }
  
  for(t in 1:(nYears-1)){
    n.recruit <- apply(myResults$sims.list$z[itera,IDSex,c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
    gc()
    
    quantile95 <- quantile(n.recruit, prob=c(0.0275, 0.975))
    quantile50 <- quantile(n.recruit, prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------    5.3 SUMMARY DEMOGRAPHIC RATES ------

parameters <- c("gamma","phi","h","w")
sex <- c("M", "F")
TableState <- matrix(NA, nrow=length(parameters)+1, ncol=(nYears)*2-2)
rownames(TableState) <- c("",unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableState) <- c(unlist(lapply(YEARS[1:(length(YEARS)-1)],function(x) rep(paste(x+1,collapse =  "-"),2))))
TableState[ ] <- c(rep(sex,nYears-1))

for(s in 1:2){
  # choose the sex
  if(s==1){ results <- Results.list[["M"]] } else { results <- Results.list[["F"]] }
  for(i in 1:length(parameters)){
    # for psi and gamma that are not sate spefici 
    if(length(dim(results$sims.list[[parameters[i]]])) == 2){   
      rows <- which(rownames(TableState) == parameters[i])[1]
      col <- which(TableState[1, ] == sex[s])
      gc()
      if(length(results$mean[parameters[i]][[1]]) == 2){
        TableState[rows,col[c(2,5)]] <- paste0(
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(median(x),2), nsmall = 2)), " (",
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.025),2), nsmall = 2)), "-",
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.975),2), nsmall = 2)), ")")  
      } else {
        TableState[rows,col] <-  paste0(
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(mean(x),2), nsmall = 2))," (",
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.025),2), nsmall = 2)), "-",
          apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.975),2), nsmall = 2)), ")")
      }
    }
  }
}

## ADD DERIVED RECRUITMENT 
for(s in 1:2){
  # choose the sex
  if(s==1){
    results <- myResults_M
    results$sims.list$z <- myResultsSXYZ_MF$sims.list$z[,myResultsSXYZ_MF$sims.list$sex=="M",]
  } else {
    results <- myResults_F
    results$sims.list$z <- myResultsSXYZ_MF$sims.list$z[,myResultsSXYZ_MF$sims.list$sex=="F",]
  }
  for(t in 1:(nYears-1)){
    n.recruit <- apply(results$sims.list$z[itera,,c(t,t+1)], 1, function(x) sum( x[,1] %in% c(1) & x[,2] %in% c(2) ))
    # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
    alivetminus1 <- apply(results$sims.list$z[itera,,t], 1, function(x)sum(x %in% c(2)))
    temp <- n.recruit/alivetminus1
    col <- which(TableState[1,]==sex[s])
    gc()
    TableState["gamma",col[t]] <- paste0(
      format(round(median(temp),2), nsmall = 2), " (",
      format(round(quantile(temp,probs=0.025),2), nsmall = 2), "-",
      format(round(quantile(temp,probs=0.975),2), nsmall = 2), ")")
  }
}

##
write.csv(TableState, file=file.path(WDFigures, paste("TableParametersState.csv",sep="")))

#write latex
colnames(TableState) <- c(unlist(lapply(YEARS[1:(length(YEARS)-1)],function(x) rep(paste(x+1,collapse =  "-"),2))))
rownames(TableState)[2:3] <- c("$\\rho$","$\\phi$")

addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

print(xtable( TableState, 
              type = "latex",
              align = paste(rep("c",ncol(TableState)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableParametersState.tex"))


##SPLIT THE TABLE IN TWO 
TableState1 <- TableState[,c(1:10)]
TableState2 <- TableState[,c(11:18)]
command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState1))),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState2))),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

#SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableState1, type = "latex",
             align = paste(rep("c", ncol(TableState1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableParametersState1.tex"))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableState2, type = "latex",
             align = paste(rep("c", ncol(TableState2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableParametersState2.tex"))



## ------  5. P0  ------

widthPolygon <- 0.15

## ------    5.1 BARS  ------

pdf(file = file.path(WDFigures, "DetectionProbBars.pdf"),
    width = 9.5, height = 10)

mx <- matrix(c(1,2,5,6,9,10,13,14),4,2,byrow = T)
mx1 <- matrix(c(3,4,7,8,11,12,15,16),4,2,byrow = T)
mx <- cbind(mx,mx1)
nf <- layout(mx,widths=c(1,0.5,1,0.5),heights=c(rep(1,dim(mx)[1]-1)))
# layout.show(nf)

country.colors <- c("goldenrod1","goldenrod3")
contry.colors.samples <- c("black",grey(0.4))
names(country.colors) <- c("Norway","Sweden")
COUNTIES_AGGREGATEDSubsetsimp <- st_simplify(COUNTIES_AGGREGATEDSubset, dTolerance = 2000, preserveTopology  = T)
COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
COUNTIES_AGGREGATEDSubsetsimp$Name <- c("NO5",
                                        "NO4",
                                        "SE3",
                                        "NO3",
                                        "SE2",
                                        "NO2",
                                        "SE1",
                                        "NO1")

myCol <-rep(c("lightblue","blue"),15)
myDev <- c(-0.15,+0.15)#seq(-0.3,0.3,length.out = 2)# 

CountyIndex <- COUNTIES_AGGREGATEDSubsetsimp$idunique

index <-c(6,1,
          4,5,
          7,4,
          2,3)
index <-c(4,6,
          5,3,
          7,2,
          8,1)

for(c in index){
  par(mar = c(4,4,1,1), tck=0)
  plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,0.04), type ="n", xaxt="n",
       xlab = "Years", ylab = "Detection probability")
  axis(2,tck=-0.02)
  yrs1 <- unlist(lapply(YEARS, function(x) {
    paste0(x[1], "/"," \n", x[2])
  }))
  axis(1, c(1:nYears),labels = yrs1 ,cex.axis=0.55,padj = -0.5)
  abline(v=1:(nYears-1)+0.5,lty=2,col=grey(0.5))
  
  for(s in 1:2){
    myResults <- Results.list[[s]]
    
    for(t in 1:nYears){
      if(length(c)>0){
        tmp <- myResults$sims.list$p0[ , c, t]
        quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
        polygon(x = c(t+ myDev[s] - widthPolygon, t+ myDev[s] + widthPolygon,
                      t+ myDev[s] + widthPolygon, t+ myDev[s] - widthPolygon ),
                y = c(quantile95[1], quantile95[1],
                      quantile95[2], quantile95[2]), 
                col=adjustcolor(myCol[s], 0.6),
                border= adjustcolor(myCol[s], 0.7))
      }
    }
  }
  
  if(c==6){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.006)
    for(i in 1:2){
      points(2, 0.035-yoffset[i], pch=15,cex=1.8,col=adjustcolor(cols[i],violin.alpha95+0.1))
      text(3.5, 0.035-yoffset[i], labels[i],cex=1.3,pos=4)
    }
  }
  
  ## PLOT REGION MAP    
  par(mar=c(0,0,0,0))
  plot(st_geometry(COUNTIES_AGGREGATEDSubsetsimp),border=grey(0.5),col=grey(0.5),lwd=0.1)
  aggCounty <- COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==c, ]
  plot(st_geometry(aggCounty),
       add = T, col = adjustcolor("red",0.4), border = "red")
  text( st_coordinates(st_centroid(aggCounty)),
        labels = COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique == c, ]$Name, col="white")
}
dev.off()



## ------    5.1 BARS Other ------

pdf(file = file.path(WDFigures, "DetectionProbBarsOther.pdf"),
    width = 9, height = 8)

mx <- matrix(1:4,2,2,byrow = TRUE)
nf <- layout(mx, widths=c(1,0.5), heights=c(rep(1,dim(mx)[1]-1)))
# layout.show(nf)

country.colors <- c("goldenrod1","goldenrod3")
contry.colors.samples <- c("black",grey(0.4))
names(country.colors) <- c("Norway","Sweden")
COUNTRIESsimp <- st_simplify(COUNTRIES,dTolerance = 1500,preserveTopology = T)
COUNTRIESsimp$idunique <- COUNTRIES$ISO

myCol <-rep(c("lightblue","blue"),15)
myDev <- c(-0.15,+0.15)#seq(-0.3,0.3,length.out = 2)# 

yrs1 <- unlist(lapply(YEARS, function(x){ paste0(x[1], "/"," \n", x[2]) }))

for(c in 1:2){
  par(mar=c(4,4,1,1), tck=0)
  # par(mfrow=c(1,2))
  plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,0.015), type ="n", xaxt="n",
       xlab = "Years", ylab = "Detection probability")
  axis(2,tck=-0.02)
  axis(1, c(1:nYears),labels = yrs1,cex.axis=0.9)
  abline(v=1:(nYears-1)+0.5,lty=2,col=grey(0.5))
  
  for(s in 1:2){
    myResults <- Results.list[[s]]
    for(t in 1:nYears){
      # coun <- which( detCounties.original1[[t]]  %in%  CountyIndex[c])
      if(length(c)>0){
        tmp <- myResults$sims.list$p0Oth[ , c, t]
        quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
        polygon(x = c(t+ myDev[s] - widthPolygon, t+ myDev[s] + widthPolygon,
                      t+ myDev[s] + widthPolygon, t+ myDev[s] - widthPolygon ),
                y = c(quantile95[1], quantile95[1],
                      quantile95[2], quantile95[2]), 
                col=adjustcolor(myCol[s], 0.6),
                border= adjustcolor(myCol[s], 0.7))
      }
    }#t
  }#s
  
  if(c==1){
    cols <- myCol#c("orange","darkorange2")
    pt.col <- myCol#c("white","white")
    labels <- c("Female","Male")
    yoffset <- c(0, 0.002)
    for(i in 1:2){
      #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
      points(2, 0.014-yoffset[i], pch=15,cex=2.1,col=adjustcolor(cols[i],violin.alpha95+0.1))
      #points(2, 0.035-yoffset[i], pch=15,cex=0.5,col=adjustcolor(pt.col[i],violin.alpha50))
      text(2.5, 0.014-yoffset[i], labels[i],cex=1.5,pos=4)
    }
  }
  
  par(mar=c(0,0,0,0))
  plot(st_geometry(COUNTRIESsimp),border=grey(0.5),col=grey(0.5),lwd=0.1)
  aggCounty <- nngeo::st_remove_holes(COUNTRIESsimp[c, ])
  plot(st_geometry(aggCounty),
       add=T, col=adjustcolor("red",0.4),border="red")
  if(c%in% 1){
    text(253761.2, 6775270, labels = COUNTRIESsimp[c, ]$idunique, col="white")
  }else{
    text(504651.6, 6928887, labels = COUNTRIESsimp[c, ]$idunique, col="white")
  }
}
dev.off()



## ------    5.3 TABLE  ------

CountiesID <- 1:dim(myResults_F$sims.list$p0)[2]
state <- c( "Others","Scent-marking adult")
sex <- c("M", "F")
Tablep0 <- matrix(NA, nrow=length(CountiesID)+1, ncol=(nYears)*2)
rownames(Tablep0) <- c("",unlist(lapply(as.list(CountiesID),function(x) rep(x,1))))
colnames(Tablep0) <- c(unlist(lapply(as.list(years),function(x) rep(paste(x,collapse =  "-"),2))))
Tablep0[1,] <- c( rep(sex,(nYears)) )

n.digits = 3
rownamesTablep0 <- ""
for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResults_M}else{results <- myResults_F}
  for(i in 1:length(CountiesID)){
    rows <- which(rownames(Tablep0)==CountiesID[i])
    col <- which(Tablep0[1,]==sex[s])
    # state 1 
    Tablep0[rows[1],col] <- paste0( 
      apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(median(x),n.digits), nsmall = n.digits)), " (",
      apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(quantile(x,probs=0.025),n.digits), nsmall = n.digits)), "-",
      apply(results$sims.list$p0[,CountiesID[i],],2, function(x) format(round(quantile(x,probs=0.975),n.digits), nsmall = n.digits)), ")") 
    
    rownamesTablep0[rows[1]] <- COUNTIES_AGGREGATEDSubsetsimp[COUNTIES_AGGREGATEDSubsetsimp$idunique==i, ]$Name
  }
}

rownames(Tablep0) <- rownamesTablep0
Tablep0 <- Tablep0[order(row.names(Tablep0)), ]

# WRITE THE FILE
write.csv(Tablep0, file = file.path(WDFigures, "Tablep0.csv"))



## ------  6. BETA P0S  ------

## ------    6.1 P0STRUCTURED  ------

## ------      6.1.1 BETATRACKS  ------

pdf(file = file.path(WDFigures, "Betap0StructuredTracks.pdf"),
    width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)
axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:(nYears)){
    quantile95 <- quantile(myResults$sims.list$betaCovs[,1,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovs[,1,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}#i

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------      6.1.2 BETA SNOW  ------

pdf(file = file.path(WDFigures, "Betap0StructuredSnow.pdf"),
    width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)
axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:(nYears)){
    quantile95 <- quantile(myResults$sims.list$betaCovs[,2,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovs[,2,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}#i

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------      6.1.2 BETA RESPONSE  ------

pdf(file = file.path(WDFigures, "Betap0StructuredResponse.pdf"),
    width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)
axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:(nYears)){
    quantile95 <- quantile(myResults$sims.list$betaResponse[,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaResponse[,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}#i

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------    6.2 P0OTHER  ------

## ------      6.2.1 BETASNOW ------

pdf( file = file.path(WDFigures, "Betap0OtherSnow.pdf"),
     width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:nYears){
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,1,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,1,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col = adjustcolor(myCol[s], violin.alpha95),
            border = NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col = adjustcolor(myCol[s], violin.alpha50),
            border = NA)
  }#t
}#i

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------      6.2.2 BETA ROADS  ------

pdf(file = file.path(WDFigures, "Betap0OtherRoads.pdf"),
    width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:nYears){
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,2,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,2,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col=adjustcolor(myCol[s], violin.alpha95),
            border= NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col=adjustcolor(myCol[s], violin.alpha50),
            border= NA)
  }#t
}#i

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}
dev.off()



## ------      6.2.3 BETA SKANDOBS  ------

pdf( file = file.path(WDFigures, "Betap0OtherSkandobs.pdf"),
     width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:(nYears)){
    quantile95 <- quantile(myResults$sims.list$betaCovsOth[,3,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaCovsOth[,3,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col = adjustcolor(myCol[s], violin.alpha95),
            border = NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col = adjustcolor(myCol[s], violin.alpha50),
            border = NA)
  }#t
}#i

##-- LEGEND
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(2, 3)
for(i in 1:2){
  points(4,y[i],pch=15,cex=5.5,col=adjustcolor(myCol[i],violin.alpha95))
  points(4,y[i],pch=15,cex=3,col=adjustcolor(myCol[i],violin.alpha50))
  text(5.3,y[i],labels[i],cex=1.2,pos=4)
}#i
dev.off()



## ------      6.2.3 BETA RESPONSE  ------

pdf( file = file.path(WDFigures, "Betap0OtherResponse.pdf"), 
     width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
             widths = c(0.05,1,0.30),
             heights = c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(-1,1), type ="n", xaxt="n", xlab = "Years", ylab = "Beta")
axis(2,tck=-0.02)
abline(v=1:(nYears)+0.5,lty=2)
abline(h=0,lty=1)

axis(1, c(1:nYears), labels = years+1,  cex.axis=0.8)
myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
myDev <- c(-0.15,+0.15)
ss <- c("F","M")
for(s in 1:2){
  myResults <- Results.list[[ss[s]]]
  for(t in 1:(nYears)){
    quantile95 <- quantile(myResults$sims.list$betaResponseOth[,t], prob=c(0.0275, 0.975))
    quantile50 <- quantile(myResults$sims.list$betaResponseOth[,t], prob=c(0.25, 0.75))
    polygon(x = c(t+myDev[s] - widthPolygon, t+myDev[s] + widthPolygon,
                  t+myDev[s] + widthPolygon, t+myDev[s] - widthPolygon ),
            y = c(quantile95[1], quantile95[1],
                  quantile95[2], quantile95[2]), 
            col = adjustcolor(myCol[s], violin.alpha95),
            border = NA)
    polygon(x = c(t+myDev[s]-widthPolygon, t+myDev[s]+widthPolygon,
                  t+myDev[s]+widthPolygon, t+myDev[s]-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col = adjustcolor(myCol[s], violin.alpha50),
            border = NA)
  }#t
}#i

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
labels <- c("Females", "Males")
pch <- rep(19,4)
cex <- rep(4,4)
y <- c(2, 3)
for(i in 1:2){
  points(4, y[i], pch = 15, cex = 5.5, col = adjustcolor(myCol[i],violin.alpha95))
  points(4, y[i], pch = 15, cex = 3, col = adjustcolor(myCol[i],violin.alpha50))
  text(5.3, y[i], labels[i], cex = 1.2, pos = 4)
}#i
dev.off()



## ------  6. TABLE OTHERS  ------

## ------    6.1. TABLE DENSITY & MOVEMENT  OPSCR  ------

parameters <- c("betaDens","sigma","dmean")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{dens}^*$","$\\sigma$","$\\lambda^*$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableDensityMovementSCR <- matrix(NA, nrow = length(parameters)+1, ncol = (nYears)*2)
rownames(TableDensityMovementSCR) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableDensityMovementSCR) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))) 
TableDensityMovementSCR[1,] <- c(rep(sex, (nYears)) )

for(s in 1:2){
  ##-- choose the sex
  if(s == 1){ results <- myResults_M }else { results <- myResults_F }
  
  ##-- BETA DENSITY 
  # REPEAT ESTIMATES FOR ALL YEARS FOR SUCH PARAMETERS
  param <- "betaDens"
  rows <- which(rownames(TableDensityMovementSCR)==param)[1]
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  TableDensityMovementSCR[rows,col] <- paste0(
    format(round(median(results$sims.list[[param]]), 2), nsmall = 2), " (",
    format(round(quantile(results$sims.list[[param]], probs=0.025),2), nsmall = 2), "-",
    format(round(quantile(results$sims.list[[param]], probs=0.975),2), nsmall = 2), ")")
  
  ##-- SIGMA
  param <- "sigma"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  TableDensityMovementSCR[rows,col] <- paste0( 
    format(round(apply(results$sims.list[[param]] * myHabitat.list$resolution /1000,2,median), 2), nsmall = 2), " (",
    format( round(apply(results$sims.list[[param]]* myHabitat.list$resolution/1000,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2), "-",
    format(round(apply(results$sims.list[[param]]* myHabitat.list$resolution/1000,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2), ")")
  
  ##-- LAMBDA
  param <- "dmean"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  TableDensityMovementSCR[rows,col] <- paste0( 
    format(round(median(results$sims.list[[param]] * myHabitat.list$resolution /1000), 2), nsmall = 2), " (",
    format( round(quantile(results$sims.list[[param]]* myHabitat.list$resolution/1000, probs=0.025 ),2), nsmall = 2), "-",
    format(round(quantile(results$sims.list[[param]]* myHabitat.list$resolution/1000, probs=0.975 ),2), nsmall = 2),  ")")
}#s

##-- Export .csv
write.csv(TableDensityMovementSCR, file = file.path(WDFigures, "TableDensityMovement.csv"))

##-- Export .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableDensityMovementSCR)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

colnames(TableDensityMovementSCR) <- rep("", ncol(TableDensityMovementSCR))
rownames(TableDensityMovementSCR)[2:nrow(TableDensityMovementSCR)] <- parameters1

print(xtable( TableDensityMovementSCR,
              type = "latex",
              align = paste(rep("c",ncol(TableDensityMovementSCR)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableDensityMovement.tex"))

##-- SPLIT THE TABLE IN TWO 
TableDensityMovementSCR1 <- TableDensityMovementSCR[,c(1:10)]
TableDensityMovementSCR2 <- TableDensityMovementSCR[,c(11:20)]

##-- SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable( TableDensityMovementSCR1,
              type = "latex",
              align = paste(rep("c", ncol(TableDensityMovementSCR1)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableDensityMovement1.tex", sep="")))

##-- SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableDensityMovementSCR2,
             type = "latex",
             align = paste(rep("c", ncol(TableDensityMovementSCR2)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableDensityMovement2.tex"))



## ------    6.2. TABLE PARAMETERS STRUCTURED ------

parameters <- c("betaResponse","betaCovs","betaCovs")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{1_Structured}$","$\\beta_{2_Structured}$","$\\beta_{3_Structured}$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableStructured <- matrix(NA, nrow=length(parameters1)+1, ncol=(nYears)*2)
rownames(TableStructured) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableStructured) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/"))))#unlist(lapply(YEARS,function(x) c(x[2],x[2])))# unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
TableStructured[1, ] <- c(rep(sex, nYears))

for(s in 1:2){
  ##-- choose the sex
  if(s==1){results <- myResults_M } else { results <- myResults_F }
  
  ##-- betaResponse
  param <- "betaResponse"
  rows <- which(rownames(TableStructured)==param)
  col <- which(TableStructured[1,]==sex[s])
  TableStructured[rows,col] <- paste0( 
    format(round(apply(results$sims.list[[param]] ,2,median), 2), nsmall = 2), " (",
    format( round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2), "-",
    format(round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),")")
  
  ##-- trapBetas
  for(st in 1:2){
    param <- "betaCovs"
    rows <- which(rownames(TableStructured)==param)[st]
    col <- which(TableStructured[1,]==sex[s])
    TableStructured[rows,col] <- paste0(
      format(round(apply(results$sims.list[[param]][,st,] ,2,median), 2), nsmall = 2), " (",
      format( round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2), "-",
      format(round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),")")
  }#st
}#s

##-- Export .csv
write.csv(TableStructured, file = file.path(WDFigures, "TableStructured.csv"))

##-- Export .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableStructured)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

colnames(TableStructured) <- rep("", ncol(TableStructured))
rownames(TableStructured)[2:nrow(TableStructured)] <- parameters1
print(xtable( TableStructured,
              type = "latex",
              align = paste(rep("c",ncol(TableStructured)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableStructured.tex"))

##-- SPLIT THE TABLE IN TWO 
TableStructured1 <- TableStructured[,c(1:10)]
TableStructured2 <- TableStructured[,c(11:20)]

##-- SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable( TableStructured1, 
              type = "latex",
              align = paste(rep("c", ncol(TableStructured1)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableStructured1.tex"))

##-- SAVE TABLE 2
addtorow1$command <- command2
print(xtable( TableStructured2,
              type = "latex",
              align = paste(rep("c", ncol(TableStructured2)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "TableStructured2.tex"))



## ------    6.3. TABLE  PARAMETERS  OTHERS ------

parameters <- c("betaResponseOth","betaCovsOth","betaCovsOth","betaCovsOth")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{1_Unstructured}$","$\\beta_{2_Unstructured}$","$\\beta_{3_Unstructured}$","$\\beta_{4_Unstructured}$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableOther <- matrix(NA, nrow=length(parameters1)+1, ncol=(nYears)*2)
rownames(TableOther) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableOther) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/"))))# unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
TableOther[1, ] <- c(rep(sex, nYears))

for(s in 1:2){
  ##-- choose the sex
  if(s == 1){ results <- myResults_M } else { results <- myResults_F }
  
  ##-- betaResponse
  param <- "betaResponseOth"
  rows <- which(rownames(TableOther) == param)
  col <- which(TableOther[1, ] == sex[s])
  TableOther[rows,col] <- paste0(
    format(round(apply(results$sims.list[[param]] ,2,median), 2), nsmall = 2), " (",
    format( round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2), "-",
    format(round(apply(results$sims.list[[param]] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),")")
  
  ##-- trapBetas
  for(st in 1:3){
    param <- "betaCovsOth"
    rows <- which(rownames(TableOther) == param)[st]
    col <- which(TableOther[1, ] == sex[s])
    TableOther[rows,col] <- paste0(
      format(round(apply(results$sims.list[[param]][,st,] ,2,median), 2), nsmall = 2), " (",
      format( round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2), "-" ,
      format(round(apply(results$sims.list[[param]][,st,] ,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),")")
  }#st
}#s

##-- Export .csv 
write.csv( TableOther, 
           file = file.path(WDFigures, "TableOther.csv"))

##-- Export .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther)[1:10])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableOther)[11:20])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

colnames(TableOther) <- rep("", ncol(TableOther))
rownames(TableOther)[2:nrow(TableOther)] <- parameters1

print(xtable(TableOther, type = "latex",align = paste(rep("c",ncol(TableOther)+1),collapse = "")),
      # scalebox=.7,
      floating = FALSE,
      add.to.row=addtorow,include.colnames=F,sanitize.text.function=function(x){x},
      file = file.path(WDTables,paste("TableOther.tex", sep="")))

##-- SPLIT THE TABLE IN TWO 
TableOther1 <- TableOther[,c(1:10)]
TableOther2 <- TableOther[,c(11:20)]

##-- SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableOther1, type = "latex",
             align = paste(rep("c", ncol(TableOther1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableOther1.tex", sep="")))

##-- SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableOther2, type = "latex",
             align = paste(rep("c", ncol(TableOther2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables,paste("TableOther2.tex", sep="")))



## ------  7. TRANSITION SURFACES  ------

## ------    7.1 SET UP HABITAT ------

habbR <- raster::disaggregate(myHabitat.list$habitat.r, fact=4)

#COUNTRIES <- aggregate(x = COMMUNESM, by = "ISO")
COUNTRIES <- COMMUNES %>% group_by(ISO) %>% summarize()
COUNTRIESsimp <- COUNTRIES
COUNTRIESsimp$NAME_1 <- COUNTRIES$ISO
plot(st_geometry(myHabitat.list$habitat.poly))
plot(st_geometry(COUNTRIESsimp), add = T)

## IDENTIFY COUNTIES
myStudyArea.polyMnoHoles <- nngeo::st_remove_holes(myHabitat.list$habitat.poly)
habbRCountieswbuff <- habbR <- mask(habbR, myStudyArea.polyMnoHoles)
habbR[habbR==0] <- NA
plot(habbR)
# identify SWEDEN/NORWAY in the raster
SWE <- COUNTRIESsimp[which(COUNTRIESsimp$ISO %in% c("SWE")),]     ## Just take Sweden
NOR <- COUNTRIESsimp[which(COUNTRIESsimp$ISO %in% c("NOR")),]     ## Just take Norway

this.r <- fasterize(SWE,habbR )
habbR[this.r==1] <- 2
this.r <- fasterize(NOR,habbR )
habbR[this.r==1] <- 1
plot(habbR)

##-- GET THE BUFFERLESS AREA
plot(habbR)
habbR <- mask(habbR, myStudyArea.polyMnoHoles)
plot(habbR)

# convert to text
habbR[habbR==2] <- "SWE"
habbR[habbR==1] <- "NOR"
habbR[habbR==0] <- NA
gc()

## AC DENSITY BASED
habIDCells.mx <- habbR
habIDCells.mx[] <- 1:ncell(habbR)
habIDCells.mx <- as.matrix(habIDCells.mx)

# SET UP A MATRIX : ROW = NUMBER OF REGIONS, COLUMN = NUMBER OF CELLS
# FILL IN WITH 1 AND 0 TO ASSIGN EACH CELL TO A REGION
# CELL THAT ARE NOT HABITAT OR WITHIN BUFFER ASSIGNED TO 0
regionID <- habbR[]
regionIDunique <- unique(regionID)
regionIDunique <- regionIDunique[!is.na(regionIDunique)]
regionIDmat <- do.call(rbind,lapply(regionIDunique,function(x)habbR[]== x  ))
regionIDmat[is.na(regionIDmat)] <- 0
row.names(regionIDmat) <- unique(regionIDunique)

#sourceCpp("C:/Personal_Cloud/OneDrive/Work/Coding/Rcpp/GetTransitionSurface1.cpp")
##sourceCpp("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/cpp/GetTransitionSurface1.cpp")
habbRxy <- coordinates(habbR)  
colnames(habbRxy) <- c("x","y")
myResultsSXYZ_MF$sims.list$scaledsxy <- scaleCoordsToHabitatGrid(
  coordsData = myResultsSXYZ_MF$sims.list$sxy,
  coordsHabitatGridCenter = habbRxy)$coordsDataScaled


## ------    7.2 TRANSITION PROBABILITY OTHER CAUSES  ------

TransitionSurfaceOther <- list()
iter <- sample(1:dim(myResultsSXYZ_MF$sims.list$scaledsxy)[1], size = 100)#dim(densityInputCountries$sx)[1])
for(t in 1:(nYears-1)){
  TransitionSurfaceOther[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$scaledsxy[ ,IDFemales,1,t],
                                                       myResultsSXYZ_MF$sims.list$scaledsxy[ ,IDFemales,2,t],
                                                       myResultsSXYZ_MF$sims.list$z[ ,IDFemales,t],
                                                       myResultsSXYZ_MF$sims.list$z[ ,IDFemales,t+1],
                                                       habIDCells.mx,
                                                       regionID = regionIDmat,
                                                       stateFrom = c(2),
                                                       stateTo = c(3),
                                                       ncell = ncell(habbR),
                                                       probs = c(0.025,0.975),
                                                       returnPosteriors = F)
}#t
TransitionSurfaceOther[[t]]$PosteriorTransitionRegion
TransitionSurfaceOther[[2]]$SummaryTransitionRegion

##-- TRY A SPATIAL PLOT 
SpatialRaster <- habbR
SpatialRaster[] <- TransitionSurfaceOther[[t]]$MeanCell
plot(SpatialRaster)
plot(myHabitat.list$habitat.poly,add=T)



## ------    7.2 TRANSITION PROBABILITY CULLING  ------

TransitionSurfaceCulling <- list()
for(t in 1:(nYears-1)){
  TransitionSurfaceCulling[[t]] <- GetTransitionSurface( myResultsSXYZ_MF$sims.list$scaledsxy[,,1,t],
                                                         myResultsSXYZ_MF$sims.list$scaledsxy[,,2,t],
                                                         myResultsSXYZ_MF$sims.list$z[,,t],
                                                         myResultsSXYZ_MF$sims.list$z[,,t+1],
                                                         habIDCells.mx,
                                                         regionID = regionIDmat,
                                                         stateFrom = c(2),
                                                         stateTo = c(3),
                                                         ncell = ncell(habbR),
                                                         probs = c(0.025,0.975),
                                                         returnPosteriors = F)
}#t
TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion
TransitionSurfaceCulling[[t]]$SummaryTransitionRegion

##-- TRY A SPATIAL PLOT 
SpatialRaster <- habbR
SpatialRaster[] <- TransitionSurfaceCulling[[t]]$MeanCell
plot(SpatialRaster)
plot(myHabitat.list$habitat.poly, add = T)



## ------    7.3 PLOT  ------

pdf( file = file.path(WDFigures, "CountryMortalityRates.pdf"),
     width = 9, height = 7)
par(mfrow = c(1,2))
offset <- c(-0.2,0.2)
plot(-10, xlim = c(0,nYears), ylim = c(0,1), ylab = "Mortality rate", xaxt = "n")
axis(1, at = c(1:nYears), labels = years)
col <- c("red","blue")
for(t in 1:(nYears-1)){
  for(r in 1:2){
    plot.violins(list(TransitionSurfaceOther[[t]]$PosteriorTransitionRegion[r,]),
                 at = t +offset[r],
                 x = 1, col = col[r], alpha = 0.5,add=T)
  }
}

offset <- c(-0.2,0.2)
plot(-10, xlim=c(0,nYears), ylim=c(0,1), ylab="Mortality rate culling" ,xaxt="n")
axis(1, at = c(1:nYears), labels = years)
col <- c("red","blue")
for(t in 1:(nYears-1)){
  for(r in 1:2){
    plot.violins(list(TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion[r, ]),
                 at = t + offset[r],
                 x = 1, col = col[r], alpha = 0.5, add = T)
  }#r
}#t

legend("topright", fill=col, legend=c("Norway","Sweden"))

dev.off()



## ------  8. OTHER PLOTS  ------

## ------    8.1 SKANDOBS ------

habitat.detectors <- aggregate(rasterToPolygons(disaggregate(myHabitat.list$habitat.rWthBuffer, 
                                                             fact = 2),
                                                fun = function(x)x==1))
skandobs.r1 <- skandobs.r <- disaggregate(myHabitat.list$habitat.rWthBuffer, 
                                          fact = 2)

pdf( file= file.path(WDFigures, "SkandobsRovbaseCovariates.pdf"),
     width = 12, height = 8)
par(mfrow = c(2,5), mar = c(1,1,3,1))
rrr <- list()
for(t in 1:nYears){
  plot(habitat.detectors,border=grey(0.8),col=grey(0.8))
  skandobs.r1[skandobs.r[]%in%1]<- nimData$detCovsOth[,t,3]
  rrr[[t]] <- skandobs.r1
  pol <- rasterToPolygons(skandobs.r1,fun = function(x) x==1)
  plot(pol,col="darkgreen",border=NA,add=T)
  mtext(paste(YEARS[[t]][2]))
}
dev.off()

save( rrr, file = file.path(WDFigures, "Skandobs.RData"))



## ------  8. TABLES OF #NGS SAMPLES, #DEAD RECOVERIES & #IDs DETECTED ------

## ------    8.1 OVERALL NUMBERS ------

load(file.path(WD, modelNameM, paste0(modelNameM,"_NGSData.RData")))
load(file.path(WD, modelNameF, paste0(modelNameF,"_NGSData.RData")))

##-- SOME TALLIES TO CHECK THINGS
##-- NGS
NGS <- rbind(myData.aliveF$alive, myData.aliveM$alive)
table(NGS$Year)
length(NGS)

# NGS.all <- rbind(myFullData.spF$alive,myFullData.spM$alive)
# NGS.all <- NGS.all[NGS.all$Year%in%years, ]
# table(NGS.all$Year)
# length(NGS.all)

##-- FOR REPORT SUMMARY
length(NGS$Id)
length(NGS$Id[NGS$Sex=="Hunn"])
length(NGS$Id[NGS$Sex=="Hann"])
length(NGS$Id[NGS$Country=="S"])/nrow(NGS)
length(unique(NGS$Id))
length(unique(NGS$Id[NGS$Sex=="Hunn"]))
length(unique(NGS$Id[NGS$Sex=="Hann"]))

##-- DEAD RECOVERY
dead <- rbind(myFullData.spF$dead.recovery,myFullData.spM$dead.recovery)
table(dead$Year)
length(dead)



## ------    8.2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX ------

NGSCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
NGSCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex == sex[s], ]
    NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- length(temp[temp$Country %in% "N", ])
    NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- length(temp[temp$Country %in% "S", ])
    NGSCountrySEX["Total",ye[t] + sex1[s]] <- length(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  }#t
}"s"

addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))

print(xtable( NGSCountrySEX, 
              type = "latex",
              align = paste(c("l",rep("c",ncol(NGSCountrySEX))),collapse = "")),
      floating = FALSE,
      include.colnames=F,
      add.to.row = addtorow,
      file = file.path(WDTables, "NGSCountrySEX.tex"))



## ------    8.2. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------

NGSidCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
NGSidCountrySEX[1,] <- rep(c("F","M"),nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N" ]))
    NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
    NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
  }#t
}#s

colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))

write.csv( NGSidCountrySEX, 
           file = file.path(WDTables, "NGSidCountrySEX.csv"))

addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
print( xtable(NGSidCountrySEX,
              type = "latex",
              align = paste(c("l",rep("c",ncol(NGSidCountrySEX))),collapse = "")),
       floating = FALSE,
       include.colnames = FALSE,
       add.to.row = addtorow,
       file = file.path(WDTables, "NGSidCountrySEX.tex"))

##-- PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
NGSidCountryTotal <- matrix(0, ncol = nYears, nrow = 1)
row.names(NGSidCountryTotal) <- c("Total")
colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))
for(t in 1:nYears){
  temp <- NGS[NGS$Year == years[t] , ]
  NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
}#t

##-- PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR PER SEX
write.csv( NGSidCountryTotal,
           file = file.path(WDTables, "TotalIdDetected.csv"))




## ------    8.3. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------

DeadidCountrySEX <- matrix(0, ncol = nYears*2+1, nrow = 6)
row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))
DeadidCountrySEX[1,] <- c("",rep(c("F","M"),nYears))
DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,nYears*2,by=2)

MortalityNames <- unique(as.character(dead$DeathCause))
table(as.character(dead$DeathCause))
legalCauses <- MortalityNames[grep("JF", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("9", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("23", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("28", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Rifle", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("18", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("17", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Uspesifisert", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Fellefangst", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Hagle", MortalityNames)])

## SEPARATE MORTALITIES
cause <- c("other","legal culling")
for(t in 1:nYears){
  for(s in 1:2){
    for(d in 1:2){
      if(d==1){temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$DeathCause %in% legalCauses), ]
      }else{
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$DeathCause %in% legalCauses, ]
      }
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[ ,1] == "Norway" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country %in% "N"]))
      row <- which(rownames(DeadidCountrySEX)==cause[d] & DeadidCountrySEX[ ,1] == "Sweden" )
      DeadidCountrySEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
    }#t
    DeadidCountrySEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(DeadidCountrySEX[2:6,ye[t] + sex1[s]+1]))
  }
}

##summary
#Other causes
sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
#legal
sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))
sum(as.numeric(DeadidCountrySEX[c(3,5),2:ncol(DeadidCountrySEX)]))/sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))

#write latex
addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
# REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC
multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))

print(xtable( DeadidCountrySEX,
              type = "latex",
              align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(WDTables, "DeadidCountrySEX.tex"))



## ------ 9. plot regions ------

NewCountySwe <- readOGR(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp"))
plot(NewCountySwe, border="red")

COMMUNES_NOR <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp", sep = ""))   ## Communal map of Norway
COMMUNES_SWE <- readOGR(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp", sep = ""))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)

## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- aggregate(x = COMMUNES, by = "NAME_1")

## only select the norwegian counties
COMMUNESNOR <- COMMUNES[COMMUNES$NAME_0=="Norway",]
NORWAY <- aggregate(x = COMMUNESNOR, by = "NAME_1")
plot(NORWAY)
NAME_1 <- as.character(NORWAY$NAME_1)
df.CountiesRegions <- matrix(c(
  "Finnmark",         8,
  "Troms",            8,
  "Nordland",         7,
  NAME_1[15],         6,
  NAME_1[9],         6,
  NAME_1[8],         6,
  "Hedmark",          5,
  "Oppland",          3,
  NAME_1[2],          4,
  "Oslo",             4,
  "Akershus",         4,
  "Sogn og Fjordane", 1,
  "Hordaland",        1,
  "Rogaland" ,        1,
  "Vest-Agder",       1,
  "Aust-Agder",       2,
  "Telemark" ,        2,
  "Buskerud" ,        2,
  "Vestfold",         2), byrow=T,ncol=2)

NORWAY$NAME_1 <- as.character(NORWAY$NAME_1) 
for(i in 1:nrow(df.CountiesRegions)){
  NORWAY$NAME_1[NORWAY$NAME_1 %in% df.CountiesRegions[i,1]] <- df.CountiesRegions[i,2]
}
NORWAY1 <- aggregate(NORWAY,by="NAME_1")
plot(NORWAY1)

##-- RENAME THE FIELDS SO THEY MATCH BETWEEN NORWEGIAN AND SWEDISH LAYERS
NewCountySwe <- NewCountySwe[,"LANSNAMN"]
colnames(NewCountySwe@data) <- "NAME_1"
NewCountySwe$Country <- "SWE"
NORWAY1$Country <- "NOR"

# MERGE THE 2 LAYERS.
# THERE IS SOME SPACE BETWEEN THE TWO LAYERS, BUT IT DOESNT MATTER. 
# IF A CELL IS NOT ASSIGNED TO ANY COUNTY THEN WE ASSIGN IT THE CLOSEST COUNTY BELOW. 
COUNTIESsimp <- rbind(NORWAY1,NewCountySwe)#, NewCountySwe, makeUniqueIDs = TRUE) 

country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
col <- c("firebrick2", "deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway", "Sweden")
border.col <- NA

pdf(file = file.path(WDTables, "RegionMaps.pdf"),
    width = 9, height = 13, pointsize = 12)
par(mar=c(0,0,0,0))
plot(gSimplify(COUNTIESsimp, tol=5000), col=NA, border=NA)

##-- NORWAY
CARNIVORE.REGIONS <- COUNTIESsimp[COUNTIESsimp$Country%in% "NOR",]
CARNIVORE.REGIONS1 <- gSimplify(CARNIVORE.REGIONS, tol=200, topologyPreserve = T)
CARNIVORE.REGIONS1$Region <- CARNIVORE.REGIONS$NAME_1
region <- gSimplify(NORWAY1, tol=500)

##-- SPECIFY COLOR PALETTE
this.col <- sequential_hcl(1+length(unique(NORWAY1$NAME_1)), "Reds 3")
this.col <- this.col[-length(this.col)]
set.seed(100)
this.col <- sample(this.col)
plot(region,
     add = T,
     col = this.col,
     border = border.col,
     lwd = 1)

##-- SWEDEN
NewCountySwe$NAME_1 <- as.character(NewCountySwe$NAME_1)
NewCountySwe1 <- gSimplify(NewCountySwe,tol=200,topologyPreserve = T)
swedenCounties2 <- NewCountySwe1
swedenCounties2$NAME_1 <- as.character(NewCountySwe$NAME_1)
swedenCounties2 <- swedenCounties2[!(swedenCounties2$NAME_1%in% "Gotlands lÃ¤n"), ]

##-- SPECIFY COLOR PALETTE
this.col <- sequential_hcl(25+length(unique(swedenCounties2$NAME_1)), "Blues 3")
this.col <- this.col[1:length(unique(swedenCounties2$NAME_1))]
set.seed(100)
this.col <- sample(this.col)

plot( RemoveHolesSp(swedenCounties2),
      add = T,
      col = this.col,
      border = border.col,
      lwd = 1)

##-- LABELS
CARNIVORE.REGIONS2 <- aggregate(CARNIVORE.REGIONS1, by = "Region")
raster::text(NORWAY1, labels = NORWAY1$NAME_1, cex=1.3,
             col=ifelse(NORWAY1$NAME_1 %in% c(5), grey(0.7), grey(0)))

swedenCounties2$NAME_1 <- as.character(unlist(lapply( strsplit(swedenCounties2$NAME_1, " "),
                                                      function(x) x[1])))
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="GÃ¤vleborgs"] <- "Gävleborg"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤rmlands"] <- "Värmland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–stergÃ¶tlands"] <- "Östergötland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¤mtlands"] <- "Jämtland"
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="JÃ¶nkÃ¶pings"] <- "Jönköping" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SkÃ¥ne"] <- "Skåne" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stmanlands"] <- "Västmanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤stra"] <- "Västra Götaland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sternorrlands"] <- "Västernorrland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="VÃ¤sterbottens"] <- "Västerbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="SÃ¶dermanlands"] <- "Södermanland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Ã–rebro"] <- "Örebro" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Norrbottens"] <- "Norrbotten" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Hallands"] <- "Halland" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Kronobergs"] <- "Kronoberg" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Stockholms"] <- "Stockholm" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 
swedenCounties2$NAME_1[swedenCounties2$NAME_1=="Dalarnas"] <- "Dalarna" 

swedenCounties2$Region <- c("Mellersta",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Södra",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Mellersta",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Norra",
                            "Norra","Mellersta")
swedenCounties2Regions <- aggregate(swedenCounties2,by="Region")

plot( RemoveHolesSp(swedenCounties2Regions),
      add = T,
      border = grey(0.0),
      lwd = 3)

polygonsLabel( swedenCounties2,
               swedenCounties2$NAME_1,
               method = "buffer",
               cex = 1.1,
               col = ifelse(swedenCounties2$NAME_1%in%c("Jämtland"),grey(0.7),grey(0.7)))
dev.off()



##------------------------------------------------------------------------------

## ------ VI. COMBINE OUTPUT(S) ------

## ------   1. LOAD & PROCESS NIMBLE OPSCR OUTPUTS ------

# List the directories containing bite outputs
for(thisSex in c("Hann","Hunn")){
  modelName <- modelName#"23.J_Fa"
  outDirectories <- list.files(file.path(working.dir, modelName,thisSex))[grep(paste("NimbleOutFOR",modelName,sep=""), 
                                                                               list.files(file.path(working.dir, modelName,thisSex)))]
  path.list <- file.path(working.dir, modelName,thisSex, outDirectories)
  
  # Retrieve the minimum number of bites per chain
  numBites <- unlist(lapply(path.list, function(x){
    files <- list.files(x)
    files <- files[grep(".RData", files)]
    length(files)
  }))
  numBites
  minBites <- min(numBites)/2
  
  
  # Niter to Remove (burn-in) #CM
  NSkipBites <- 3
  nthin <- 1
  nimOutput <- RUNTIME <- list()
  
  for(p in 1:length(path.list)){
    print(path.list[p])
    outfiles <- list.files(path.list[p])
    out <- runtime <- list()#[CM]
    for(x in NSkipBites:minBites){
      print(x)
      load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
      runtime[[x]] <- RunTime[3] 
      #params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
      #parmIndex <- which(! params.simple %in% c("sxy","z"))
      nthins <- seq(1,dim(this.sample)[1], by=nthin)
      out[[x]] <- this.sample[nthins,]#[ ,parmIndex] 
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
  
  
  {#doall[CM]
    
    ## ------ 2. PLOT PARAMETERS ESTIMATES ------
    
    pdf(file = file.path( working.dir,
                          modelName,
                          thisSex,
                          paste0(modelName,thisSex,".pdf")))
    
    
    
    ## ------    2.1.PLOT DETECTIONS ------
    
    myLayOut.mx <- cbind(c(2,1), c(1,1))
    myLayOut <- layout(myLayOut.mx, width = c(1,2), heights = c(1,2))
    #layout.show(myLayOut)
    
    ## PLOT STUDY AREA
    par(mar = c(0,0,0,0))
    #plot(myHabitat$buffered.habitat.poly,  col = rgb(34/250, 139/250, 34/250, alpha = 0))
    plot(st_geometry(GLOBALMAP), col = "gray80")#, add = TRUE)
    plot(st_geometry(myHabitat$habitat.poly), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
    plot(st_geometry(myHabitat$buffered.habitat.poly), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
    
    ## PLOT DETECTIONS
    plot(st_geometry(myFilteredData.sp$alive), pch = 3, col = "darkred", cex = 0.5,add=T)
    plot(st_geometry(myFilteredData.sp$dead.recovery), pch = 3, col = "darkblue", cex = 0.5,add=T)
    
    ## ADD NUMBER OF DETECTIONS
    graphics::text(x = 190000, y = 7820000, cex = 1.5,
                   labels = paste(dim(myFilteredData.sp$alive)[1], "NGS samples"))
    
    ## ADD NUMBER OF INDIVIDUALS DETECTED
    graphics::text(x = 190000, y = 7880000, cex = 1.5,
                   labels = paste(length(unique(myFilteredData.sp$alive$Id)), "Individuals"))
    
    
    
    ## ------    2.1.N ------ 
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(10, xlim = c(0, nYears+1), ylim = c(200,1200), type ="n", xaxt="n", xlab = "Years", ylab = "N")
    axis(1, c(1:nYears),labels = years)
    for(t in 1:nYears){
      plot.violins(list(myResults$sims.list$N[,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = "firebrick3",
                   add = T,
                   alpha = 0.3,
                   border.col = "firebrick3")
    }#t
    
    params <- dimnames(nimOutput[[1]])[[2]][grep("N",dimnames(nimOutput[[1]])[[2]])]
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
    }
    
    
    
    ## ------    2.2.h ------
    
    # par(mfrow = c(1,1))
    # plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "h")
    # axis(1, c(1:nYears),labels = years)
    # myCol <- "firebrick3"
    # for(t in 1:(nYears-1)){
    #   plot.violins(list(myResults$sims.list$h[ ,t]),
    #                x = t,
    #                at = t,
    #                violin.width = 0.3,
    #                col = myCol,
    #                add = T,
    #                alpha = 0.3,
    #                border.col = myCol)
    # }#t
    # 
    # params <- dimnames(nimOutput[[1]])[[2]][grep("h\\[",dimnames(nimOutput[[1]])[[2]])]
    # for(i in 1:length(params)){
    #   PlotJagsParams(jags.samples = nimOutput, params = params[i])
    # }
    
    
    
    ## ------    2.3.gamma ------ 
    
    par(mfrow = c(1,1))
    plot(10, xlim = c(0, nYears), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "gamma")
    axis(1, c(1:nYears),labels = years)
    myCol <- "firebrick3"
    for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$gamma[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.3,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
    }#t
    
    params <- dimnames(nimOutput[[1]])[[2]][grep("gamma",dimnames(nimOutput[[1]])[[2]])]
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
    }
    
    
    
    ## ------    2.4.w ------ 
    
    # par(mfrow = c(1,1))
    # plot(10, xlim = c(0, nYears+1), ylim = c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "w")
    # axis(1, c(1:nYears),labels = years)
    # for(t in 1:(nYears-1)){
    #   plot.violins(list(myResults$sims.list$w[ ,t]),
    #                x = t,
    #                at = t,
    #                violin.width = 0.3,
    #                col = myCol,
    #                add = T,
    #                alpha = 0.3,
    #                border.col = myCol)
    # }#t
    # 
    # params <- dimnames(nimOutput[[1]])[[2]][grep("w",dimnames(nimOutput[[1]])[[2]])]
    # for(i in 1:length(params)){
    #   PlotJagsParams(jags.samples = nimOutput, params = params[i])
    # }
    
    ## ------    2.5.phi ------  
    
    load(file.path(working.dir,
                   modelName,
                   thisSex,
                   paste0(modelName,thisSex,"_Chain",1, ".RData")))
    
    phi <- phiind1 <- culled <- recruit <- 0
    
    recruit <- 0
    z_caculate <- nimData$z
    z_caculate[is.na(z_caculate)] <- 0
    for(t in 2:dim(nimData$z)[2]){
      #phi
      alivet <- which(z_caculate[,t-1] %in% c(2))
      phi[t-1] <- sum(z_caculate[alivet,t] %in% c(2))/length(alivet)
      #culled
      #culled[t-1] <- sum(z_caculate[alivet,t] %in% c(3))/length(alivet)
      
      #recru
      notentered <- which(z_caculate[,t-1] == 0)
      recruit[t-1] <- sum(z_caculate[notentered,t] %in% c(2))/sum(z_caculate[,t-1] %in% c(2))
    }
    phi        # overall phi
    recruit    # per capita recruitment based on all individuals present at t-1
    #culled
    par(mfrow = c(1,1))
    plot(-10, type = "n", xaxt = "n",
         xlim = c(0,nYears), ylim=c(0,1),
         xlab = "Years", ylab = "Realized phi from z")
    axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
    yr <- c(1:(nYears-1))
    points(phi ~ yr, pch = 16)
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "phi")
    axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
    myCol <- c("firebrick3")
    for(t in 1:(nYears-1)){
      plot.violins(list(myResults$sims.list$phi[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.2,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
    }#t
    
    params <- dimnames(nimOutput[[1]])[[2]][grep("phi",dimnames(nimOutput[[1]])[[2]])]
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
    }
    
    
    
    ## ------    2.6.p0 ------
    
    ## by country and trap-response
    par(mfrow = c(1,2))
    myDev <- c(-0.3, 0.3)
    myCol <- c("blue4", "yellow3")
    
    COUNTIES_AGGREGATEDSubsetsimp <- st_simplify(COUNTIES_AGGREGATEDSubset, dTolerance = 500,preserveTopology =  T)
    COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
    
    for(c in 1:dim(myResults$sims.list$p0)[2]){
      plot(st_geometry(myStudyArea))
      plot(st_geometry(COUNTIES_AGGREGATEDSubset[COUNTIES_AGGREGATEDSubsetsimp$idunique %in% c, ]), add=T, col="red")
      
      plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
      axis(1, at = 1:(nYears), labels = years[1:(nYears)])
      
      for(t in 1:nYears){
        plot.violins(list(myResults$sims.list$p0[ ,c,t]),
                     x = t ,
                     at = t ,
                     violin.width = 0.2,
                     col = "red",
                     add = T,
                     alpha = 0.3,
                     border.col = "red")
      }#t
    }#c
    
    par(mfrow = c(1,2))
    myDev <- c(-0.3, 0.3)
    myCol <- c("blue4", "yellow3")
    
    # COUNTIES_AGGREGATEDSubsetsimp <- gSimplify(COUNTIES_AGGREGATEDSubset,tol = 500,topologyPreserve = T)
    # COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
    for(c in 1:(dim(myResults$sims.list$p0Oth)[2])){
      plot(st_geometry(myStudyArea))
      if(c %in% 3){ 
        plot(st_geometry(myDetectors$main.detector.sp[detCounties %in% 1,]), add=T, col="red")
      } else {
        plot(st_geometry(COUNTRIES[c, ]), add=T, col="red")
      }
      
      plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
      axis(1, at = 1:(nYears), labels = years[1:(nYears)])
      
      for(t in 1:nYears){
        plot.violins(list(myResults$sims.list$p0Oth[ ,c,t]),
                     x = t ,
                     at = t ,
                     violin.width = 0.2,
                     col = "red",
                     add = T,
                     alpha = 0.3,
                     border.col = "red")
      }#t
    }#c
    
    params <- dimnames(nimOutput[[1]])[[2]][grep("p0", dimnames(nimOutput[[1]])[[2]])[-1]]
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
    }
    
    ## ------    2.7.the rest ------
    #structured
    par(mfrow = c(1,1))
    trapCovs <- c("Tracks","Snow")
    for(b in 1:dim(myResults$sims.list$betaCovs)[2]){
      plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n", xlab = "beta Structured", ylab = "beta",main=trapCovs[b])
      abline(h=0)
      axis(1, at = 1:nYears , labels = years)
      myCol <- c("firebrick3")
      for(t in 1:dim(myResults$sims.list$betaCovs)[3]){
        plot.violins(list(myResults$sims.list$betaCovs[ ,b,t]),
                     x = t,
                     at = t,
                     violin.width = 0.2,
                     col = myCol,
                     add = T,
                     alpha = 0.3,
                     border.col = myCol)
      }#t
    }
    #other
    par(mfrow = c(1,1))
    trapCovs <- c("Snow","Road","Skandobs")
    for(b in 1:dim(myResults$sims.list$betaCovsOth)[2]){
      plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n", xlab = "beta Other", ylab = "beta",main=trapCovs[b])
      abline(h=0)
      axis(1, at = 1:nYears , labels = years)
      myCol <- c("firebrick3")
      for(t in 1:dim(myResults$sims.list$betaCovsOth)[3]){
        plot.violins(list(myResults$sims.list$betaCovsOth[ ,b,t]),
                     x = t,
                     at = t,
                     violin.width = 0.2,
                     col = myCol,
                     add = T,
                     alpha = 0.3,
                     border.col = myCol)
      }#t
    }
    
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears+0.5), ylim=c(0,10000), type ="n", xaxt="n",
         xlab = "beta Other", ylab = "sigma")
    #abline(h=0)
    axis(1, at = 1:nYears , labels = years)
    myCol <- c("firebrick3")
    for(t in 1:dim(myResults$sims.list$sigma)[2]){
      plot.violins(list(myResults$sims.list$sigma[ ,t]*myHabitat$resolution),
                   x = t,
                   at = t,
                   violin.width = 0.2,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
    }#t
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n",
         xlab = "beta Other", ylab = "betaResponse")
    abline(h=0)
    axis(1, at = 1:nYears , labels = years)
    myCol <- c("firebrick3")
    for(t in 1:dim(myResults$sims.list$betaResponse)[2]){
      plot.violins(list(myResults$sims.list$betaResponse[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.2,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
    }#t
    
    
    par(mfrow = c(1,1))
    plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n",
         xlab = "beta Other", ylab = "betaResponseOth")
    abline(h=0)
    axis(1, at = 1:nYears , labels = years)
    myCol <- c("firebrick3")
    for(t in 1:dim(myResults$sims.list$betaResponseOth)[2]){
      plot.violins(list(myResults$sims.list$betaResponseOth[ ,t]),
                   x = t,
                   at = t,
                   violin.width = 0.2,
                   col = myCol,
                   add = T,
                   alpha = 0.3,
                   border.col = myCol)
    }#t
    
    # par(mfrow = c(1,1))
    # plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n",
    #      xlab = "beta Other", ylab = "beta Dens")
    # abline(h=0)
    # axis(1, at = 1:nYears , labels = years)
    # myCol <- c("firebrick3")
    # for(t in 1:dim(myResults$sims.list$betaDens)[2]){
    #   plot.violins(list(myResults$sims.list$betaDens[ ,t]),
    #                x = t,
    #                at = t,
    #                violin.width = 0.2,
    #                col = myCol,
    #                add = T,
    #                alpha = 0.3,
    #                border.col = myCol)
    # }#t
    
    params <- c( "lambda",  "betaDens")
    for(i in 1:length(params)){
      PlotJagsParams(jags.samples = nimOutput, params = params[i])
    }#i
    dev.off()
  }#doall[CM]
}#thisSex



## ------   0.SET ANALYSIS CHARACTERISTICS -----

myVars <- list( 
  ## WORKING DIRECTORY & MODEL NAME
  working.dirSCR2025 = file.path(dir.dropbox,"wolverine/CM/2025"),
  working.dirSCROPSCR2024 = file.path(dir.dropbox,"wolverine/CM/2024"),
  
  modelNameSCR2025 = "54.Cleaned2025",
  modelNameOPSCR2024 = "plot53Cleaned2024",
  
  ## HABITAT SPECIFICATIONS
  HABITAT = list( countries =  c("SWE","NOR"),
                  habResolution = 20000,
                  habBuffer = 60000),
  
  ## NGS DATA SPECIFICATIONS
  DATA = list( years = 2015:2024, #2014:2023
               species = c("Jerv"),              
               sex = c("Hann","Hunn"),                   
               sampling.months = list(12,1:6)),   
  ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)
  
  # DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 2000,
                    detResolution = 10000,
                    detDeadResolution = 15000),
  
  # DATA GENERATION 
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 3,
                     aug.factor = 0.8),
  
  ## OUTPUT PLOTS 
  OUTPUT = list(mapResolution = 10000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)


years <- DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

yearsSampledNorrb <- c(2016:2018,2023,2024)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)

if(!dir.exists(file.path(working.dirSCR2025,"FigureTableReport"))){
  dir.create(file.path(working.dirSCR2025, "FigureTableReport"))
}



## ------   1. ABUNDANCE TIME SERIES ------

## ------     1.1 LOAD DATA OPSCR2024 ------

load(file.path(working.dirSCROPSCR2024,modelNameOPSCR2024,"Figure","CICounties.RData"))
#store the data          
quantile95Tot2024 <- quantile95Tot
quantile50Tot2024 <- quantile50Tot
quantile95Swe2024 <- quantile95Swe
quantile50Swe2024 <- quantile50Swe
quantile95Nor2024 <- quantile95Nor
quantile50Nor2024 <- quantile50Nor



## ------     1.2 LOAD DATA SCR2025 ------

load(file.path(working.dirSCR2025,
               modelNameSCR2025,
               "FigureSnap/CICounties.RData"))

## add SCR results to the list
quantile95Tot2024[[nYears+1]] <- quantile95Tot[[nYears]]
quantile50Tot2024[[nYears+1]] <- quantile50Tot[[nYears]]
quantile95Swe2024[[nYears+1]] <- quantile95Swe[[nYears]]
quantile50Swe2024[[nYears+1]] <- quantile50Swe[[nYears]]
quantile95Nor2024[[nYears+1]] <- quantile95Nor[[nYears]]
quantile50Nor2024[[nYears+1]] <- quantile50Nor[[nYears]]



## ------     1.3 PLOT ------

SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

#define colors
TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")#c("turquoise","darkmagenta")# c("goldenrod1","goldenrod3")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE

##pdf
pdf(file= file.path(working.dirSCR2025, 
                    "FigureTableReport/NCountriesBarsSeason.pdf"),
    width = 14, height = 10)

#plot
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(1.5, nYears+.5+1), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(2:(nYears+1)), labels = SeasonText, cex.axis=1.2,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

for(t in 2:(nYears+1)){
  ##-- TOTAL
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95Tot2024[[t]][1], quantile95Tot2024[[t]][1],
                quantile95Tot2024[[t]][2], quantile95Tot2024[[t]][2]), 
          col=adjustcolor(TotalColors, violin.alpha95),
          border= NA)
  
  # add a star
  if(sum((t-1) %in% yearsNotSampled)){
    text(x = t+widthPolygon+offsetstar,
         y = quantile95Tot2024[[t]][2],
         "*", cex = cexStar)
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50Tot2024[[t]][1], quantile50Tot2024[[t]][1],
                  quantile50Tot2024[[t]][2], quantile50Tot2024[[t]][2]), 
            col = adjustcolor(TotalColors, violin.alpha50),
            border = NA)
  }
  
  ##-- SWEDEN
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95Swe2024[[t]] [1], quantile95Swe2024[[t]] [1],
                quantile95Swe2024[[t]] [2], quantile95Swe2024[[t]] [2]), 
          col = adjustcolor(country.colors[2], violin.alpha95),
          border = NA)
  
  if(sum((t-1) %in% yearsNotSampled)){
    text(x= t+offsetstar ,y= quantile95Swe2024[[t]][2], "*",cex=cexStar)
  }
  if(displayQuantiles50){
    polygon(x = c(t, t-widthPolygon1*2,
                  t-widthPolygon1*2, t),
            y = c(quantile50Swe2024[[t]][1], quantile50Swe2024[[t]][1],
                  quantile50Swe2024[[t]][2], quantile50Swe2024[[t]][2]), 
            col = adjustcolor(country.colors[2], violin.alpha50),
            border = NA)
  }
  
  ##-- NORWAY
  polygon(x = c(t, t+widthPolygon1*2,
                t+widthPolygon1*2, t),
          y = c(quantile95Nor2024[[t]][1], quantile95Nor2024[[t]][1],
                quantile95Nor2024[[t]][2], quantile95Nor2024[[t]][2]), 
          col=adjustcolor(country.colors[1], violin.alpha95),
          border= NA)
  
  if(displayQuantiles50){
    polygon(x = c(t, t+widthPolygon2*2,
                  t+widthPolygon2*2, t),
            y = c(quantile50Nor2024[[t]][1], quantile50Nor2024[[t]][1],
                  quantile50Nor2024[[t]][2], quantile50Nor2024[[t]][2]), 
            col=adjustcolor(country.colors[1], violin.alpha50),
            border= NA)
  }
}#t
box()
abline(v = at[1:(nYears)]+0.5, lty = 2)

##-- legend
par(xpd = TRUE)
polygon(x = c(0.8,7.2,7.2,0.8)+1,
        y = c(170,170,230,230),
        col = adjustcolor("white",alpha.f = 0.9),
        border = "white")
labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <- c(200, 200, 200)
x <- c(1,3,5)+1
mycol1 <- c(country.colors, "black")
# add transparent background polygon
for(i in 1:3){
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}#i

dev.off()



## ------   2. UD TIME SERIES ------

## ------     2.1. LOAD OBJECTS ------

##-- LOAD UD MAPS (OPSCR2024 + SCR2025 for the last year)
UDmap <- list()
for(t in 1:(nYears-1)){
  path <- file.path( working.dirSCROPSCR2024,
                     modelNameOPSCR2024,
                     "Figure/RasterForRovbase",
                     paste0("wolverine_5km",paste(YEARS[[t]][1],collapse = "_"),".tif"))
  UDmap[[t]] <- raster(path)
}#t

t <- nYears
path <- file.path( working.dirSCR2025,
                   modelNameSCR2025,
                   "FigureSnap/RasterForRovbase",
                   paste0("wolverine_5km",paste(YEARS[[t]][1],collapse = "_"),".tif"))
UDmap[[nYears]] <- raster(path)

##-- Load necessary objects SCR 2025
load(file.path(working.dirSCR2025, modelNameSCR2025, "NecessaryObjects.RData"))
myHabitat.listOPSCR2024 <- myHabitat.list

##-- Load necessary objects OPSCR 2024
load(file.path(working.dirSCROPSCR2024, modelNameOPSCR2024, "NecessaryObjects.RData"))
myHabitat.listSCR2025 <- myHabitat.list



## ------     2.2. UD-DENSITY TIME SERIES ------

pdf(file = file.path(working.dirSCR2025, "FigureTableReport/DensityMapsUD.pdf"),
    width = 12, height = 8)

##-- layout
mx <- rbind(c(1,rep(1:5, each=2)),
            c(rep(1:5, each=2),5))
mx <- rbind(mx, mx+5)
nf <- layout(mx, widths = c(rep(1,ncol(mx))), heights=rep(1,2))
#layout.show(nf)

max <- max(unlist(lapply(UDmap, function(x) max(x[], na.rm = T))))
cuts <- seq(0, max, length.out = 100) #set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)
habbdensUDCropped <- list()
for(t in 1:length(years)){
  par(mar = c(0,0,0,0))#, bg="transparent")#-country polygons will only show on dark background
  plot(st_geometry(COUNTRIESsimpFig), border = NA, col = grey(0.85))
  ##-- BECAUSE raster::plot MESSES UP THE LAYOUT
  image(UDmap[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  # plot(RemoveHolesSp(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  plot(st_geometry(COUNTRIESsimpFig), border = grey(0.4), col = NA, add = TRUE)
  
  mtext(SeasonText[[t]], 1, -2, adj=0.15, cex=1.2)
  ##-- if OPSCR results used, use habitat from OPSCR model otherwise use the SCR ones
  if(years[t] < 2024){
    agg1 <- aggregate(rasterToPolygons(myHabitat.listOPSCR2024$habitat.rWthBuffer,
                                       fun = function(x) x == 1))
  } else {
    agg1 <- aggregate(rasterToPolygons(myHabitat.listSCR2025$habitat.rWthBuffer,
                                       fun = function(x) x == 1))
  }
  
  if(sum(t%in%yearsNotSampled)){
    plot(st_geometry(COUNTIESNorrbotten),add=T, lwd=2)
  }
  
  if(t==nYears){
    segments(x0= 1320000,x1=1320000,
             y0= 6900000,y1=6900000 + 1000000, col=grey(0.3), lwd=2, lend=2)  
    text(1280000, 6900000+1000000/2,labels="1000 km", srt=90 )
    
    plot(UDmap[[t]], legend.only=TRUE,breaks=cuts, col=col,
         legend.width = 2,
         axis.args=list(at=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        labels=round(seq(0, max-0.05, length.out = 4), digits = 1),
                        cex.axis=1.6),
         smallplot=c(0.72, 0.75, 0.2, 0.4),
         legend.args=list(text=expression(paste("Individuals/100 km"^ 2, "", sep="")),
                          side=4, font=1, line=4, cex=1.2))
  }
}#t
dev.off()



## ------     2.3. WRITE RASTER & COG ------

#run all but last year. the last year comes from the
for(t in 1:(length(years))){
  setworking.dir(file.path(working.dirSCR2025, modelNameSCR2025, "FigureSnap/RasterForRovbase"))
  r <- UDmap[[t]]
  
  ## write raster
  writeRaster(x = r,
              filename = file.path( working.dirSCR2025, "FigureTableReport/RasterForRovbase",
                                    paste0("wolverine_5km",years[t],".tif")),
              overwrite=TRUE)
  
  ## reproject as asked by rovbase
  r_3006 <- terra::project(rast(r), "EPSG:3006")
  
  ## write Cloud Optimized GeoTIFF should be 256*256
  writeRaster(
    r_3006,
    file.path(working.dirSCR2025, "FigureTableReport/CogRasterForRovbase",
              paste0("wolverine_5km",years[t],"_cog.tif")),
    filetype = "COG",
    overwrite = TRUE,
    gdal = c( "BLOCKSIZE=256",
              "COMPRESS=DEFLATE",
              "PREDICTOR=2"))
}#t



## ------   3. TABLE TIME SERIES ------

## ------     3.1 NORWEGIAN REGIONS AND SWEDISH ------

## ------     3.1. 2 LOAD DATA OPSCR2024 AND SCR 2025 ------

##-- OPSCR2024
NAllYears2024 <- read.csv(file.path(working.dirSCROPSCR2024,modelNameOPSCR2024,"Table/NAllYears.csv"))

##-- SCR2025
NAllYears2025 <- read.csv(file.path(working.dirSCR2025,modelNameSCR2025,"TableSnap/NAllYears.csv"))


## TABLE
idcountyTable <- NAllYears2024$X
idcountyTable <- iconv(idcountyTable, from = "latin1", to = "UTF-8")

##-- CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimates <- NAllYears2024[,2:(ncol(NAllYears2024))]
NCarRegionEstimates[,1:(ncol(NCarRegionEstimates)-1)] <- NAllYears2024[,3:(ncol(NAllYears2024))]#matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimates) <- c(idcountyTable)
colnames(NCarRegionEstimates) <-  unlist(SeasonText)#unlist(lapply(YEARS ,function(x) x[2]))#
NCarRegionEstimates[,SeasonText[[nYears]]] <- NAllYears2025[,ncol(NAllYears2025)]


## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

##-- print csv
write.csv(NCarRegionEstimates,
          file = file.path(working.dirSCR2025,"FigureTableReport",paste("NAllYears.csv",sep="")),
          fileEncoding="latin1")

idcounty1[which(idcounty1 %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten}"
row.names(NCarRegionEstimates) <- idcounty1
NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*","}", sep="")

NCarRegionEstimates["SWEDEN",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["SWEDEN",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["TOTAL",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["TOTAL",yearsNotSampled], "**","}", sep="")
NCarRegionEstimates["Nordre",yearsNotSampled] <- paste("\\textcolor[gray]{.5}{",NCarRegionEstimates["Nordre",yearsNotSampled], "**","}", sep="")

idcountyNOR <- idcounty1[3:10]
idcountySWENorth <- idcounty1[13:16]
idcountySWEMiddle <- idcounty1[18:25]
idcountySWESouth <- idcounty1[27:28]

row.names(NCarRegionEstimates) <- c("TOTAL",
                                    paste("\\hspace{0.25cm}", "NORWAY",sep=""),
                                    paste("\\hspace{0.5cm} ", idcountyNOR,sep=""),
                                    paste("\\hspace{0.25cm}", "SWEDEN",sep=""),
                                    paste("\\hspace{0.5cm}", "Norra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWENorth, sep=""),
                                    paste("\\hspace{0.5cm}", "Mellersta",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWEMiddle, sep=""),
                                    paste("\\hspace{0.5cm}", "Södra",sep=""),
                                    paste("\\hspace{0.75cm}", idcountySWESouth, sep=""))
row.names(NCarRegionEstimates)[grep("VästraGötaland", row.names(NCarRegionEstimates))] <- paste("\\hspace{0.75cm}",
                                                                                                "Västra Götaland", sep="")
##-- Export .tex
print(xtable(NCarRegionEstimates,
             type = "latex",
             align = paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
      # scalebox=.8,
      floating = FALSE,
      sanitize.text.function = function(x){x},
      add.to.row = list(list(seq(1,nrow(NCarRegionEstimates), by = 2)),"\\rowcolor[gray]{.96} "),
      file = file.path( working.dirSCR2025,
                        "FigureTableReport/NCountiesCarnivoreRegions.tex"))



## ------     3.2. NORWEGIAN COUNTIES AND SWEDISH ------

## ------       3.2.1 LOAD DATA OPSCR2024 ------

NAllYearsNorwegianCounties2024 <- read.csv(file.path(working.dirSCROPSCR2024,
                                                     modelNameOPSCR2024,
                                                     "Table/NAllYearsNorwegianCounties.csv"))


## ------       3.2.2. LOAD DATA SCR2025 ------

NAllYearsNorwegianCounties2025 <- read.csv(file.path( working.dirSCR2025,
                                                      modelNameSCR2025,
                                                      "TableSnap/NAllYearsNorwegianCounties.csv"))



## ------       3.2.3. COMBINE IN ONE TABLE ------

idcountyTable <- NAllYearsNorwegianCounties2024$X
idcountyTable[idcountyTable %in% "\xc3fstfold"] <- "Østfold"
idcountyTable <- iconv(
  iconv(idcountyTable, from = "UTF-8", to = "latin1"),
  from = "latin1",
  to   = "UTF-8")

##-- CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimatesNOR <- NAllYearsNorwegianCounties2024[,2:(ncol(NAllYearsNorwegianCounties2024))]
NCarRegionEstimatesNOR[ ,1:(ncol(NCarRegionEstimatesNOR)-1)] <- NAllYearsNorwegianCounties2024[,3:(ncol(NAllYearsNorwegianCounties2024))]#matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimatesNOR) <- c(idcountyTable)
colnames(NCarRegionEstimatesNOR) <- unlist(SeasonText)#unlist(lapply(YEARS ,function(x) x[2]))#
NCarRegionEstimatesNOR[,SeasonText[[nYears]]] <- NAllYearsNorwegianCounties2025[,ncol(NAllYearsNorwegianCounties2025)]


## WRITE LATEX TABLE 
# ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "TOTAL")] <- "NORWAY"

##-- Export .csv
write.csv(NCarRegionEstimatesNOR,
          file = file.path(working.dirSCR2025,"FigureTableReport",paste("NCountiesCarnivoreRegionsNorway.csv",sep="")),
          fileEncoding="latin1")

row.names(NCarRegionEstimatesNOR) <- c("NORWAY",
                                       paste("\\hspace{0.25cm} ",
                                             idcounty1[2:length(idcounty1)],sep=""))

print(xtable(NCarRegionEstimatesNOR,
             type = "latex",
             align = paste(c("l",rep("c",ncol(NCarRegionEstimatesNOR))),collapse = "")),
      # scalebox=.8,
      floating = FALSE,
      sanitize.text.function = function(x){x},
      add.to.row = list(list(seq(1,nrow(NCarRegionEstimatesNOR),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dirSCR2025,
                       "FigureTableReport/NCountiesCarnivoreRegionsNorway.tex"))


## ------   4. GROWTH RATE ------

## ------     4.1 GET AND TRANSFORM .TEX FILE ------

tex <- readLines(file.path(working.dirSCROPSCR2024,
                           modelNameOPSCR2024,
                           "Table/growthRate.tex"),
                 encoding = "UTF-8")
tex <- gsub("\\\\rowcolor\\[[^]]*\\]\\{[^}]*\\}", "", tex)

# This strips \rowcolor[gray]{.96} but leaves "Sweden & ..." intact.
# 3. Keep only table rows
rows <- tex[grepl("&", tex)]
rows <- rows[!grepl("\\\\hline|\\\\begin|\\\\end", rows)]

# 4. Clean row endings
rows <- gsub("\\\\\\\\", "", rows)  # remove \\
rows <- trimws(rows)

# 5. Split into columns
split_rows <- strsplit(rows, "\\s*&\\s*")
growthRate <- do.call(rbind, split_rows)

# Step 6: Set row and column names
# The first row contains column names, first column contains row names.
colnames(growthRate) <- growthRate[1, ]
rownames(growthRate) <- growthRate[, 1]
growthRate <- growthRate[-1, -1]

# skip first year of OPSCR 2024
growthRate[,1:(ncol(growthRate)-1)] <- growthRate[,2:ncol(growthRate)]

# UPDATE LAST YEAR 
load(file.path( working.dirSCROPSCR2024,
                modelNameOPSCR2024,
                "Figure/posteriorRegions.RData"))
posteriorRegionsOPSCR2024 <- posteriorRegions

load(file.path( working.dirSCR2025,
                modelNameSCR2025,
                "FigureSnap/posteriorRegions.RData"))
posteriorRegionsSCR2025 <- posteriorRegions



## ------     4.2. CALCULATE GROWTH RATE OPSCR2024 TO SCR2025 ------

t = 10

##-- MAKE SURE WE HAVE THE SAME NB OF ITERATIONS BETWEEN SCR AND OPSCR
ite <- seq(1, length(posteriorRegionsOPSCR2024[[t]]["Norway", ]),
           length.out = length(posteriorRegionsSCR2025[[t]]["Norway", ]))
# UPDATE LAST GROWTH RATE ESTIMTATE
##-- NORWAY
growth <- posteriorRegionsSCR2025[[t]]["Norway", ] / posteriorRegionsOPSCR2024[[t]]["Norway",ite]
growthRate["Norway",t-1] <- paste0(
  format(round(mean(growth),digits = 2),nsmall = 2)," (",
  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2), ")")

##-- SWEDEN
growth <- posteriorRegionsSCR2025[[t]]["Sweden", ] / posteriorRegionsOPSCR2024[[t]]["Sweden",ite]
growthRate["Sweden",t-1] <- paste0(
  format(round(mean(growth),digits = 2),nsmall = 2)," (",
  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),")")

##--TOTAL
growth <- posteriorRegionsSCR2025[[t]][c("Sweden","Norway"), ] / posteriorRegionsOPSCR2024[[t]][c("Sweden","Norway"),ite]
growthRate["Total",t-1] <- paste0(
  format(round(mean(growth),digits = 2),nsmall = 2)," (",
  format(round(quantile(growth, probs=c(0.025)), digits = 2),nsmall = 2),"-",
  format(round(quantile(growth, probs=c(0.975)), digits = 2),nsmall = 2),")")

##-- UPDATE YEARS COLUMN
colnames(growthRate) <- unlist(lapply(YEARS, FUN = function(x) paste(x, collapse ="-")))[2:length(YEARS)] #paste(x[[2]]))

##-- Export .tex
print(xtable( growthRate,
              type = "latex",
              align = paste(c("l", rep("c",ncol(growthRate))), collapse = "")),
      # scalebox=.8,
      floating = FALSE,
      sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(growthRate),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dirSCR2025, "FigureTableReport/growthRate.tex"))



## ------   5. PROPORTION OF IDs DETECTED ------

## ------     5.1 GET AND TRANSFORM .TEX FILE ------

#prop detected OPSCR 2024
texDetected2024 <- readLines(file.path(working.dirSCROPSCR2024,modelNameOPSCR2024,"Table",
                                       "propDetectedCountry.tex"), encoding = "UTF-8")


tex <- gsub("\\\\rowcolor\\[[^]]*\\]\\{[^}]*\\}", "", texDetected2024)

# This strips \rowcolor[gray]{.96} but leaves "Sweden & ..." intact.
# 3. Keep only table rows
rows <- tex[grepl("&", tex)]
rows <- rows[!grepl("\\\\hline|\\\\begin|\\\\end", rows)]

# 4. Clean row endings
rows <- gsub("\\\\\\\\", "", rows)  # remove \\
rows <- trimws(rows)

# 5. Split into columns
split_rows <- strsplit(rows, "\\s*&\\s*")
Detected2024 <- do.call(rbind, split_rows)

#prop detected SCR 2025
texDetected2025 <- readLines(file.path(working.dirSCR2025,modelNameSCR2025,"Table",
                                       "propDetectedCountry.tex"), encoding = "UTF-8")


tex <- gsub("\\\\rowcolor\\[[^]]*\\]\\{[^}]*\\}", "", texDetected2025)

# This strips \rowcolor[gray]{.96} but leaves "Sweden & ..." intact.
# 3. Keep only table rows
rows <- tex[grepl("&", tex)]
rows <- rows[!grepl("\\\\hline|\\\\begin|\\\\end", rows)]

# 4. Clean row endings
rows <- gsub("\\\\\\\\", "", rows)  # remove \\
rows <- trimws(rows)

# 5. Split into columns
split_rows <- strsplit(rows, "\\s*&\\s*")
Detected2025 <- do.call(rbind, split_rows)



## ------     5.2 UPDATE LAST YEAR ------

Detected <- Detected2024
##-- skip first year
Detected[ ,2:(ncol(Detected)-2)] <- Detected2024[ ,4:(ncol(Detected))]
##-- update last year
Detected[ ,c(20,21)] <- Detected2025[ ,c(20,21)]



## ------     5.3. SHAPE TABLE ------

row.names(Detected) <- c("","Norway","Sweden","Total")
colnames(Detected) <- c("",unlist(lapply(YEARS,function(x) c(paste(x[1],x[2],sep="/"),
                                                             paste(x[1],x[2],sep="/")))))
addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(Detected)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', uniqueYEAR, '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))

print( xtable(Detected,
              type = "latex",
              align = paste(rep("c", ncol(Detected)+1), collapse = "")),
       #scalebox = .7, 
       floating = FALSE,
       add.to.row = addtorow,
       include.colnames = FALSE,
       include.rownames = FALSE,
       sanitize.text.function = function(x){x},
       file = file.path(working.dirSCR2025, "FigureTableReport/propDetectedCountry.tex"))

##-- SPLIT THE TABLE IN TWO 
propDetectedCountry1 <- Detected[,c(1:11)]
propDetectedCountry2 <- Detected[,c(1,12:21)]

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(Detected)[2:11])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(Detected)[12:21])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))


##-- SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(propDetectedCountry1,
             type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry1)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dirSCR2025, "FigureTableReport/propDetectedCountry1.tex"))

##-- SAVE TABLE 2
addtorow1$command <- command2
print(xtable(propDetectedCountry2,
             type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dirSCR2025, "FigureTableReport/propDetectedCountry2.tex"))



## ------   6. SAVE RASTER MAPS ------

years <- c(2016:2025)
t =  1

##-- Run all but last year. the last year comes from the
for(t in 1:(length(years)-1)){
  
  ##-- Load raster
  setworking.dir(file.path(working.dirSCR2025, modelNameSCR2025, "FigureSnap/RasterForRovbase"))
  r <- rast(paste0("wolverine_5km",years[t],".tif"))
  
  ##-- Export .tif
  writeRaster( x = r,
               filename = file.path(working.dirSCR2025,
                                    "FigureTableReport/RasterForRovbase",
                                    paste0("wolverine_5km",years[t],".tif")),
               overwrite = TRUE)
  
  ##-- Reproject as asked by rovbase
  r_3006 <- project(r, "EPSG:3006")
  
  ##-- Write Cloud Optimized GeoTIFF should be 256*256
  writeRaster( r_3006,
               file.path( working.dirSCR2025,
                          "FigureTableReport/CogRasterForRovbase",
                          paste0("wolverine_5km", years[t], "_cog.tif")),
               filetype = "COG",
               overwrite = TRUE,
               gdal = c( "BLOCKSIZE=256",
                         "COMPRESS=DEFLATE",
                         "PREDICTOR=2"))
}#t



##------------------------------------------------------------------------------