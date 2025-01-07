## ------ IMPORT REQUIRED LIBRARIES ------

library(raster)
library(coda)
library(nimble)
library(spdep)
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

library(rovquantR)
library(dplyr)



## -----------------------------------------------------------------------------

## ------ 0. SET ANALYSIS CHARACTERISTICS ------

# ## WORKING DIRECTORY & MODEL NAME
# WD = "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/CM/2024"
# modelName = "53.aJ_FaCleaned2024"

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/2024"

##-- HABITAT SPECIFICATIONS
HABITAT = list( countries = c("SWE","NOR"),
                habResolution = 20000,
                habBuffer = 60000)

##-- NGS DATA SPECIFICATIONS
DATA = list( years = 2014:2023,
             species = c("Jerv"),           
             sex = c("Hunn"),                 
             samplingMonths = list(12,1:6))

##-- DETECTORS SPECIFICATIONS
DETECTORS = list( detSubResolution = 2000,
                  detResolution = 10000,
                  detDeadResolution = 15000)

##-- DATA GENERATION
DETECTIONS = list( maxDetDist = 40000,
                   resizeFactor = 3,
                   aug.factor = 0.8)

##-- OUTPUT PLOTS
OUTPUT = list(mapResolution = 10000)

##-- MISCELLANEOUS
plot.check = TRUE

years <- DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))



## -----------------------------------------------------------------------------

## ------ I. LOAD AND SELECT DATA ------

## ------   1. HABITAT DATA ------

## ------     1.1. LOAD RAW SHAPEFILES ------

##-- Load map of Scandinavia (including Finland & parts of Russia)
data(GLOBALMAP, envir = environment()) 

##-- Load polygons of Sweden and Norway
data(COUNTRIES, envir = environment())
COUNTRIES <- COUNTRIES %>% 
  filter(ISO %in% c("SWE","NOR")) %>%
  group_by(ISO) %>%
  summarize()

##-- Load polygons of counties in Sweden and Norway
data(COUNTIES, envir = environment()) 

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
  dplyr::summarise() 



## ------     1.2. CREATE STUDY AREA POLYGON ------

studyArea <- COUNTRIES %>%
  filter(ISO %in% HABITAT$countries) %>%
  mutate(id = 1) %>%
  group_by(id) %>% summarize()



## ------   2. NGS DATA ------

## ------     2.1. LOAD CLEANED ROVBASE FILES ------

# ## NGS data from RovBase#[CM update to 20190627]
# DNA <- read.csv(file.path(data.dir,
#                           "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dna_wolverines.csv"),
#                 fileEncoding = "latin1")
# 
# ## Dead Recoveries from RovBase#[CM update to 20190627]
# DEAD <- read.csv(file.path(data.dir,
#                            "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dead_carnivores.csv"),
#                  fileEncoding = "latin1")

##-- Extract date from the last cleaned data file
DATE <- getMostRecent(
  path = file.path(working.dir,"data"),
  pattern = "CleanData_wolverine")

##-- Load the most recent Bear data from RovBase
myFullData.sp <- readMostRecent(
  path = file.path(working.dir,"data"),
  pattern = "CleanData_wolverine",
  extension = ".RData")
  


## ------     2.2. LOAD OTHER FILES ------

## DNA samples to be removed from Henrik
SUSPECT_NGS_SAMPLES <- read.csv(file.path(data.dir,
                                          "Remove ngs samples list wolverine 2024.csv"),
                                fileEncoding = "latin1") 

## DNA samples to be removed from Henrik
SUSPECT_DeadRecoSAMPLES <- read.csv(file.path(data.dir,
                                              "Remove dead recoveries list wolverine 2024.csv"),
                                    fileEncoding = "latin1") 
## ??
HairTrapSamples <- read_xlsx(file.path(data.dir,
                                       "hairtrapsNB2024.xlsx"))

## Wolverine den locations
DEN <- read.csv(file.path(data.dir,
                      "DEN_COUNTS_2009_2024_fromHB.csv"),
                fileEncoding = "latin1")

## Skandobs
skandObs <- read_xlsx(file.path(data.dir, "Skandobs/Richard_Bischof_Skandobs_2012_2024dd.xlsx"))

## RovBase detections
rovbaseObs1 <- read_xlsx(file.path(data.dir, "ALL SPECIES IN SEPERATE YEARS/RIB2810202415264376.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(data.dir, "ALL SPECIES IN SEPERATE YEARS/RIB28102024152348493.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(data.dir, "ALL SPECIES IN SEPERATE YEARS/RIB28102024152447860.xlsx"))
rovbaseObs4 <- read_xlsx(file.path(data.dir, "ALL SPECIES IN SEPERATE YEARS/RIB28102024152538742.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3,rovbaseObs4)



## ------     2.3. TRANSLATE SCANDINAVIAN CHARACTERS ------

# colnames(DNA) <- translateForeignCharacters(
#   dat = colnames(DNA),
#   dir.translation = dir.analysis)
# 
# # drop a column that makes cleanDataNew to fail
# DNA <- DNA[ ,-which(colnames(DNA)%in% "Kjoenn..Individ.")]
# 
# colnames(DEAD) <- translateForeignCharacters(
#   dat = colnames(DEAD),
#   dir.translation = dir.analysis)

colnames(DEN) <- translateForeignCharacters(
  dat = colnames(DEN),
  dir.translation = dir.analysis)

colnames(skandObs) <- translateForeignCharacters(
  dat = colnames(skandObs),
  dir.translation = dir.analysis)

colnames(rovbaseObs) <- translateForeignCharacters(
  dat = colnames(rovbaseObs),
  dir.translation = dir.analysis)

rovbaseObs$Proevetype <- translateForeignCharacters(
  dat = rovbaseObs$Proevetype,
  dir.translation = dir.analysis)



## ------   3. SEARCH EFFORT DATA ------

## ------     3.1. GPS SEARCH TRACKS ------

## LOAD TRACKS
TRACKS_SINGLE <- read_sf(file.path(
  data.dir,
  "GIS/SearchTracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp"))
TRACKS_MULTI <- read_sf(file.path(
  data.dir,
  "GIS/SearchTracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))

## COMBINE ALL TRACKS
ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI) %>%
  filter(Helikopter == "0",
         Jerv == "1")

## SELECT TRACKS YEAR
TRACKS_YEAR <- list()
for(t in 1:nYears){
  ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
  TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%DATA$samplingMonths[[1]], ]
  TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%DATA$samplingMonths[[2]], ]
  TRACKS <- rbind(TRACKS_1, TRACKS_2)
  TRACKS <- st_intersection(TRACKS, st_as_sf(studyArea))
  ## NAME TRACKS
  TRACKS$ID <- 1:nrow(TRACKS)
  TRACKS_YEAR[[t]] <- TRACKS
}#t



## ------     3.2. DISTANCE TO ROADS ------

## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(file.path(data.dir, "GIS/Roads/MinDistAllRoads1km.tif"))
r <- fasterize(st_as_sf(studyArea), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, studyArea)



## ------     3.3. DAYS OF SNOW ------

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
SNOW <- stack(file.path(data.dir,"GIS/Snow/AverageSnowCoverModisSeason2008_2024_Wolf.tif"))

## RENAME THE LAYERS
names(SNOW) <- paste(2008:2023,(2008:2023)+1, sep = "_")

## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep="")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))



## ------     3.4. LOAD SCANDINAVIAN 20KM HABIAT ------

##-- Load pre-defined habitat rasters and shapefiles
data(habitatRasters, envir = environment()) 

##-- Disaggregate habitat raster to the desired resolution
habRaster <- raster::disaggregate(
  x = habitatRasters[["Habitat"]],
  fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)




## -----------------------------------------------------------------------------

## ------ II. CREATE SCR DATA ------

## ------   1. CLEAN AND FILTER NGS DATA ------

## ------     1.1. CLEAN NGS & DEAD RECOVERY DATA ------

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]

## Remove un-verified dead recoveries [HB]
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]

## Clean data 
myCleanedData.sp <- CleanDataNew2sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(DATA$samplingMonths)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)



## ------     1.2. FILTER DATA ------

myFullData.sp <- FilterDatasf( 
  myData = myCleanedData.sp,
  poly = studyArea,
  dead.recovery = T,
  sex = c("Hann","Hunn"),
  setSex = T)



## ------       1.2.1. FILTER SUSPECT SAMPLES ACCORDING TO HENRIK ------

myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]



## ------       1.2.2. FILTER INDIVIDUALS THAT DIES TWICE ------

myFullData.sp$dead.recovery$Id <- as.character(myFullData.sp$dead.recovery$Id)
IdDoubleDead <- myFullData.sp$dead.recovery$Id[duplicated(myFullData.sp$dead.recovery$Id)]

if(length(IdDoubleDead) > 0){
  duplicatedDeath <- NULL
  for(i in IdDoubleDead){
    tmp <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
    if(length(tmp)==0){tmp <- which(myFullData.sp$dead.recovery$Id == i)[-2]}##[CM] remove the second record.
    duplicatedDeath <- c(duplicatedDeath, tmp)
  }#i
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]
}#if



## ------       1.2.3. FILTER OUT PUPS ------

## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and 
## recovered dead between March and November
sum(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
      myFullData.sp$dead.recovery$Month > 2 &
      myFullData.sp$dead.recovery$Month < 12)

myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
                                                                    myFullData.sp$dead.recovery$Month > 2 &
                                                                    myFullData.sp$dead.recovery$Month < 12),]


## 2) remove individuals that have a weight > 0 and < 4 between March and November
## format the weight correctly
myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)

## Convert to decimals
myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))
## Get the two weight columns together.
myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt),
                                             myFullData.sp$dead.recovery$Helvekt,
                                             myFullData.sp$dead.recovery$Slaktevekt)
## Assign negative values to NAs to avoid issues
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999

## Check with Henrik
## this step does not remove dead recoveries on id with weight==0 should it?

## check how many dead reco we remove and remove if more than 0
if(sum(myFullData.sp$dead.recovery$weight > 0 &
       myFullData.sp$dead.recovery$weight < 4 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
                                                                      myFullData.sp$dead.recovery$weight < 4 &
                                                                      myFullData.sp$dead.recovery$Month < 12 &
                                                                      myFullData.sp$dead.recovery$Month > 2),]
} 

## check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
                                myFullData.sp$dead.recovery$Month < 12 &
                                myFullData.sp$dead.recovery$Month > 2,  ]
}


## ------       1.2.4. FILTER FOR DATES ------

myFilteredData.sp <- list()# myFullData.sp

## Subset to years and months of interest
myFilteredData.sp$alive <- myFullData.sp$alive %>%
  filter(Year %in% years,
         Month %in% unlist(DATA$samplingMonths))
  
## Subset to years of interest
myFilteredData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
  filter(Year %in% years)



## ------       1.2.5. FILTER SAMPLES IN NORRBOTTEN IN ALL YEARS EXCEPT 2017, 2018 and 2019 ------ 

COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018,2023)

##-- Identify samples in Norrbotten
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  mutate(is.Norr = as.numeric(st_intersects(., COUNTIESNorrbotten)))

## check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(myFilteredData.sp$is.Norr)), ]$Year)

# subset
myFilteredData.sp$alive <- myFilteredData.sp$alive %>%
  filter(!Year %in% yearsSampledNorrb,
         !is.na(is.Norr))



## ------       1.2.6. FILTER FOR SEX ------

## SELECT THE SEX
# MYFULLDATA
myFullData.sp <- myFullData.sp
myFullData.sp$alive <- myFullData.sp$alive[myFullData.sp$alive$Sex %in% DATA$sex,]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Sex %in% DATA$sex,]

# myFilteredData
myFilteredData.spAllSex <- myFilteredData.sp
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Sex %in% DATA$sex,]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Sex %in% DATA$sex,]




## ------     1.4. SEPARATE STRUCTURED AND OPPORTUNISTIC SAMPLING ------
## ------       1.4.1. ASSIGN SAMPLES TO TRACKS  ------
## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
myFilteredData.sp$alive$TrackRovbsID <- NA
myFilteredData.sp$alive$TrackDist <- NA

TRACKSSimple_sf <- list()
for(t in 1:nYears){
  TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbsID)
  TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]#st_simplify(TRACKS_YEAR[[t]],preserveTopology = T,dTolerance = 1000)
}

## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
dnatemp <- st_as_sf(myFilteredData.sp$alive)
## CREATE A BUFFER AROUND EACH DETECTION
tmp <-  st_buffer(dnatemp, dist = 750)
for(i in 1:nrow(myFilteredData.sp$alive)){
  # INTERSECT POINT WITH TRACKS,
  t <- which(years %in% tmp[i,]$Year)
  # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
  whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
  
  tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])
  
  if(nrow(tmpTRACKS)==0){next}
  # FIND THE CLOSEST TRACK
  dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
  # IF NO MATCHING DATE ASSIGN TO NA?
  if(length(dist)==0){
    myFilteredData.sp$alive$TrackRovbsID[i] <- NA
    myFilteredData.sp$alive$TrackDist[i] <- NA
  }
  # IF  MATCHING DATE ASSING TO THAT TRACK
  if(length(dist)==1){
    myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
    myFilteredData.sp$alive$TrackDist[i] <- dist
  }
  # IF SEVERAL MATCHING DATES ASSING TO THE CLOSEST OF THE MATCHING TRACKS
  if(length(dist)>1){
    myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
    myFilteredData.sp$alive$TrackDist[i] <- min(dist)
  }
  print(i)
  #if(is.na(myFilteredData.sp$alive$TrackRovbsID[i])){print(i)}
}



## ------       1.4.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC AND STRUCTURED ------
distanceThreshold <- 500

## Proeveleverandoer columns was replaced by two columns, merging them now...
myFilteredData.sp$alive$Proeveleverandoer <-  ifelse(myFilteredData.sp$alive$Annen.innsamler...Rolle %in% "" , 
                                                     myFilteredData.sp$alive$Samlet.selv...Rolle, myFilteredData.sp$alive$Annen.innsamler...Rolle)

whichStructured <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold

myFilteredData.spStructured <- myFilteredData.sp$alive[whichStructured, ]
myFilteredData.spOthers <- myFilteredData.sp$alive[!whichStructured, ]



## ------     1.5. SEPARATE MORTALITY CAUSES ------ 

## MORTALITY CAUSES
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))

## DEFINE LEGAL MORTALITY
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])

## SPLIT MORTALITY CAUSES
legal.death <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]
Other.death <- myFilteredData.sp$dead.recovery[!myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]



## ------   2. GENERATE HABITAT ------

## ------     2.1. REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ------

## DELINEATE A BUFFER AROUND ALL DETECTIONS 
myBufferedArea <- myFilteredData.spAllSex$alive %>%
  st_buffer(., dist = HABITAT$habBuffer * 1.4) %>%
  mutate(idd = 1) %>%
  group_by(idd) %>%
  summarize()

## CUT TO SWEDISH AND NORWEGIAN BORDERS
studyArea <- st_intersection(myBufferedArea, studyArea)



## ------     2.2. GENERATE HABITAT CHARACTERISTICS FROM THE NEW HABITAT DEFINITION ------

habitat <- MakeHabitatFromRastersf(
  poly = studyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = HABITAT$habBuffer,                               
  plot.check = T)

## RETRIEVE HABITAT WINDOWS BOUNDARIES
lowerHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1,] - 0.5*HABITAT$habResolution
upperHabCoords <- coordinates(habitat$habitat.r)[habitat$habitat.r[]==1,] + 0.5*HABITAT$habResolution
nHabCells <- dim(lowerHabCoords)[1]

## CREATE HABITAT GRID 
habIDCells.mx <- habitat$IDCells.mx 
habIDCells.mx[] <- 0
scaledHabGridCenters <- scaleCoordsToHabitatGrid(
  coordsData = habitat$habitat.xy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid = F)$coordsHabitatGridCenterScaled

scaledHabGridCenters <- scaledHabGridCenters[habitat$habitat.r[] == 1, ]
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}



## ------       2.2.1. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 

## Remove samples outside the STUDY AREA 
studyArea$idd <- 1
studyAreaAggregated <- studyArea %>% group_by(idd) %>% summarize()

whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, studyAreaAggregated))))
if(length(whichOut) > 0){
  myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
}
myFilteredData.sp$alive$Id <- droplevels( myFilteredData.sp$alive$Id)

## REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT 
whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, habitat$buffered.habitat.poly))))
if(length(whichOutBuff)>0){
  myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]
}

### check correlation number of detections ~ between monitoring season
myFilteredData.sp$dead.recovery$Id <- as.character(  myFilteredData.sp$dead.recovery$Id)
myFilteredData.sp$alive$Id <- as.character( myFilteredData.sp$alive$Id)
deadID <- unique(myFilteredData.sp$dead.recovery$Id)
ndet <- NULL
timeDiff <- NULL
for(i in 1:length(deadID)){
  tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Year
  
  timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Date-
    as.POSIXct(strptime(paste("01-12",tmpYear,sep="-"), "%d-%m-%Y")) 
  
  ndet[i] <- length(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] & 
                                              myFilteredData.sp$alive$Year %in% tmpYear,])
}#i



## ------     2.3. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       2.3.1. DEN COUNTS ------

DEN.sp <-  st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id  <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[,("id")]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD(as(DEN.sp,"Spatial"),
             h = 30000,
             grid = as(habitat$habitat.r,'SpatialPixels'))
  ))

if(plot.check){
  plot(DEN.r)#,main=paste("Sum Den-2013:2018Kernel H=", H[h], sep=""))
  plot(st_geometry(studyArea), add = TRUE, border = "black")
}

## EXTRACT COVARIATEs
denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
denCounts <- round(scale(denCounts), digits = 2)



## ------   3. GENERATE DETECTORS ------

## ------     3.1. GENERATE DETECTORS CHARACTERISTICS ------

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
habitat.subdetectors <- disaggregate(
  habitat$habitat.rWthBuffer,
  fact = res(habitat$habitat.r)[1]/DETECTORS$detSubResolution)

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
detectors <- MakeSearchGridsf(
  data = habitat.subdetectors,
  resolution = DETECTORS$detResolution,
  div = (DETECTORS$detResolution/DETECTORS$detSubResolution)^2,
  plot = FALSE,
  fasterize = TRUE)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(detectors$main.detector.sp)[1]
n.detectors.dead <- dim(detectors.dead$main.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(detectors$main.detector.sp)
n.trials <- as.vector(table(detectors$detector.sp$main.cell.id))
detector.dead.xy <- st_coordinates(detectors.dead$main.detector.sp)

## IDENTIFY DETECTORS IN NORBOTTEN 
COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten","Troms","Västerbotten",
                                                            "Nordland","Finnmark"),]
COUNTIESAroundNorrbotten <- st_simplify(COUNTIESAroundNorrbotten, dTolerance = 500)

## CREATE A NORROBOTTEN DETECTOR GRID
distDestsCounties <- st_distance(detectors$main.detector.sp, COUNTIESAroundNorrbotten,byid = T)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min) == 3)

## PLOT CHECK 
plot(st_geometry(COUNTIESAroundNorrbotten))
plot(st_geometry(detectors$main.detector.sp), col="black",pch=16,cex=0.3,add=T)
plot(st_geometry(detectors$main.detector.sp[detsNorrbotten,]), col="red",pch=16,cex=0.3,add=T)

## RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * DETECTORS$detResolution
lowerDetCoords.dead <- detector.dead.xy - 0.5 * DETECTORS$detDeadResolution
upperDetCoords.dead <- detector.dead.xy + 0.5 * DETECTORS$detDeadResolution



## ------     3.2. GENERATE DETECTOR-LEVEL COVARIATES ------
## ------       3.2.1. EXTRACT COUNTRIES ------
dist <- st_distance(detectors$main.detector.sp, COUNTRIES, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))



## ------       3.2.2. EXTRACT COUNTIES ------
## ASSIGN COUNTIES TO DETECTORS
dist <- st_distance(detectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F )
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties),]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))



## ------       3.2.3. EXTRACT GPS TRACKS LENGTHS ------
## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(detectors$main.detector.sp),
                                      rep(1,nrow(detectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(studyArea)
detectorGrid$id <- 1:nrow(detectorGrid)

##CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
TRACKS.r <- list()
for(t in 1:nYears){
  TRACKSst <- TRACKS_YEAR[[t]]
  intersection <- st_intersection(detectorGrid, TRACKSst) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN))## Get total length searched in each detector grid cell
  # transect_N = length(unique(ID)))## Get total number of visits in each detector grid cell
  # transect_qi = mean(QI))          ## Get mean transects quality index for each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[,t]
}

## PLOT CHECK 
max <- max(unlist(lapply(TRACKS.r, function(x) max(x[],na.rm=T))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))
CountriesDetRes <- disaggregate(habitatRasters$Countries,fact=2)
CountriesDetRes <- crop(CountriesDetRes,TRACKS.r[[1]])
rr <- TRACKS.r[[1]]
rr[CountriesDetRes[]%in% 2] <- 1

sum(st_length(TRACKS_YEAR[[t]]))/1000
sum(TRACKS.r[[t]][],na.rm=T)/1000


pdf(file=file.path(WD, modelName, paste(modelName,"Tracks.pdf",sep="")))
if(plot.check){
  NORTRACKS <- SWETRACKS <- 0
  # par(mfrow = c(2,2))#[CM]
  for(t in 1:nYears){
    plot( TRACKS.r[[t]],main=years[t], breaks=cuts, col = col,legend=FALSE)#[CM]
    plot(st_geometry(habitat$habitat.poly), main = years[t],add=T)
    plot( TRACKS.r[[t]], legend.only=TRUE, breaks=cuts, col=col,
          legend.width = 2,
          axis.args=list(at=round(seq(0, max, length.out = 5), digits = 1),
                         labels=round(seq(0, max, length.out = 5), digits = 1),
                         cex.axis=0.6),
          legend.args=list(text='', side=4, font=2, line=2.5, cex=0.8))
    # points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year==years[t], ], col="red", pch=16, cex=0.8)
    ## summary tracks
    NORTRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 2],na.rm = T )/1000
    SWETRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 4],na.rm = T )/1000
  }#t
  years1 <- years+1
  plot(SWETRACKS~(years1),  col=country.colors[2],
       lwd=2, pch=16, type="b", ylim=c(0,300000), ylab="sum tracks km")
  lines(NORTRACKS~(years1), col=country.colors[1], lwd=2, pch=16, type="b")
  legend("topright",c("N","S"), fill=country.colors)
  
}
dev.off()





## ------       3.2.4. EXTRACT DISTANCES TO ROADS ------
## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)

## if NA returns the average value of the cells within 20000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract(DistAllRoads, detectors$main.detector.sp[isna,], buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp



## ------       3.2.5. EXTRACT DAYS OF SNOW ------
## EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = nYears)

det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp



## ------       3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------

## ------          3.2.6.1. SKANDOBS ------

## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))

## MAKE IT SPATIAL 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(studyArea))

## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist(DATA$samplingMonths))
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) 
skandObs <- skandObs[subset,] 

## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(habitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in% 1, ]
subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace,] 

## RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate(habitat.subdetectors, fact=(DETECTORS$detResolution/DETECTORS$detSubResolution))
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



## ------          3.2.6.2. ROVBASE ------

## GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`),]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

## DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(studyArea)

## SUBSET THE DATA 
filter <- list(
  species = "Jerv",
  type = c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)","Saliv/Spytt"),
  month = unlist(DATA$samplingMonths))

## SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset,] 

## SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$`Art (Proeve)`) 
rovbaseObs.sp <- rovbaseObs.sp[-subset,] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace,] 

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



## ------          3.2.6.3. COMBINE ROVBASE AND SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <-    1
}



## ------          3.2.6.4. SMOOTH THE BINARY MAP ------

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



## ------          3.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------

## IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.spAllSex$alive[which(myFilteredData.spAllSex$alive$DNAID%in% HairTrapSamples$DNAID),]

## MANUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL. 
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr)
  whereHair <-raster::extract(r.SkandObsOtherSamplesBinary[[t]],tmpHair,cellnumbers=T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[,1]] <- 1
  plot(r.SkandObsOtherSamplesBinary[[t]])
  plot(tmpHair$geometry,add=T,col="red")
}


## ------          3.2.6.6. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract(r.SkandObsOtherSamplesBinary, detectors$main.detector.sp)


## ------       3.2.7. SCALE AND ROUND DETECTOR-LEVEL COVARIATES ------

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

# CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}



## ------   4. RESCALE COORDINATES ------
# HABITAT
ScaledLowerCoords <- scaleCoordsToHabitatGrid(
  coordsData = lowerHabCoords,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =T )$coordsDataScaled
ScaledUpperCoords <- scaleCoordsToHabitatGrid(
  coordsData = upperHabCoords,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =T )$coordsDataScaled

ScaledUpperCoords[ ,2] <- ScaledUpperCoords[ ,2]+1
ScaledLowerCoords[ ,2] <- ScaledLowerCoords[ ,2]-1

# DETECTORS
colnames(detector.xy) <- c("x","y")
ScaledDetectors <- scaleCoordsToHabitatGrid(coordsData = detector.xy,
                                            coordsHabitatGridCenter = habitat$habitat.xy,
                                            scaleToGrid =T )$coordsDataScaled






## ------   5. GENERATE y DETECTION ARRAYS ------

## ------     5.1. GENERATE NGS DETECTIONS : y.alive[i,j,t] ------

## ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$alive,                
  detectors = detectors.dead$main.detector.sp,
  mysubDetectors = detectors.dead$detector.sp,
  radius = DETECTORS$detResolution)

## STRUCTURED
myData.aliveStruc <- AssignDetectors_v3sf( 
  myData = myFilteredData.spStructured,                
  detectors = detectors$main.detector.sp,
  mysubDetectors = detectors$detector.sp,
  radius = DETECTORS$detResolution)

## OTHERS
myData.aliveOthers <- AssignDetectors_v3sf( 
  myData = myFilteredData.spOthers,                
  detectors = detectors$main.detector.sp,
  mysubDetectors = detectors$detector.sp,
  radius = DETECTORS$detResolution)

## DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$dead.recovery,
  detectors = detectors.dead$main.detector.sp,
  radius = DETECTORS$detResolution)

### MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET ASSIGNED 
### TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING FIND THE CASES
### WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE OF NORRBOTTEN
whichdets <- which(!myData.alive$myData.sp$Year %in% yearsSampledNorrb &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)
whichdetsStruc <- which(!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveStruc$myData.sp$Detector %in% detsNorrbotten)
whichdetsOther <- which(!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveOthers$myData.sp$Detector %in% detsNorrbotten)

## ASSIGN DETECTORS 
subdetector.sf <- detectors$detector.sp

## ALL
for(i in 1:length(whichdets)){
  tmp <- myData.alive$myData.sp[whichdets[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, detectors$main.detector.sp)
  #Artificially increase distance for detectors in noorrbotten 
  dist[,detsNorrbotten] <- 500000
  idMain <- which.min(dist[1,])
  myData.alive$myData.sp$Detector[whichdets[i]] <- idMain 
  # SUBDETECTOR
  dist <- st_distance(st_as_sf(tmp), subdetector.sf )
  #Artificially increase distance for detectors in noorrbotten
  dist[,detsNorrbotten] <- 500000
  idSub <- which.min(dist[1,])
  myData.alive$myData.sp$sub.detector[whichdets[i]] <- idSub 
}

## STRUCTURED
for(i in 1:length(whichdetsStruc)){
  tmp <- myData.aliveStruc$myData.sp[whichdetsStruc[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, detectors$main.detector.sp)
  #Artificially increase distance for detectors in noorrbotten 
  dist[,detsNorrbotten] <- 500000
  idMain <- which.min(dist[1,])
  myData.aliveStruc$myData.sp$Detector[whichdetsStruc[i]] <- idMain 
  # SUBDETECTOR
  dist <- st_distance(tmp, subdetector.sf )
  #Artificially increase distance for detectors in noorrbotten
  dist[,detsNorrbotten] <- 500000
  idSub <- which.min(dist[1,])
  myData.aliveStruc$myData.sp$sub.detector[whichdetsStruc[i]] <- idSub 
}

## OTHER
for(i in 1:length(whichdetsOther)){
  tmp <- myData.aliveOthers$myData.sp[whichdetsOther[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, detectors$main.detector.sp)
  #Artificially increase distance for detectors in noorrbotten 
  dist[,detsNorrbotten] <- 500000
  idMain <- which.min(dist[1,])
  myData.aliveOthers$myData.sp$Detector[whichdetsOther[i]] <- idMain 
  # SUBDETECTOR
  dist <- st_distance(tmp, subdetector.sf )
  #Artificially increase distance for detectors in noorrbotten
  dist[,detsNorrbotten] <- 500000
  idSub <- which.min(dist[1,])
  myData.aliveOthers$myData.sp$sub.detector[whichdetsOther[i]] <- idSub 
}



## ------     5.2. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ------

## ALL SAMPLES
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 detectors = detectors$main.detector.sp,
                 method = "Binomial",
                 myData2 = myData.dead,
                 detectors2 = detectors.dead$main.detector.sp,
                 returnIdvector = TRUE)
y.ar.ALIVE <- y.ar$y.ar
dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)

## STRUCTURED
y.arStruc <- MakeYsf( myData = myData.aliveStruc$myData.sp,
                      detectors = detectors$main.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      detectors2 = detectors.dead$main.detector.sp,
                      returnIdvector = TRUE)
y.ar.ALIVEStruc <- y.arStruc$y.ar
dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)

## OTHERS
y.arOth <- MakeYsf( myData = myData.aliveOthers$myData.sp,
                    detectors = detectors$main.detector.sp,
                    method = "Binomial",
                    myData2 = myData.dead,
                    detectors2 = detectors.dead$main.detector.sp,
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
y.ar.DEAD[y.ar.DEAD>0] <- 1



## ------     5.3. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------

distances <- list()
for(t in 1:nYears){
  print(paste("------ ", t ," -------", sep = "" ))
  distances[[t]] <- CheckDistanceDetectionsV2sf( 
    y = y.ar.ALIVE[ , ,t], 
    detector.xy = detector.xy, 
    max.distance = DETECTIONS$maxDetDist,
    method = "pairwise",
    plot.check = T)
  
  ## REMOVE DETECTIONS THAT ARE FURTHER THAN  THE THRESHOLD
  y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  
  ## REMOVE DETECTIONS ALSO IN MYDATA TO RUN GETSINITS
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
    
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]),]
  }
}#t



## ------     5.4. GENERATE INDIVIDUAL-LEVEL COVARIATES ------

## ------       5.4.1. TRAP-RESPONSE ------

## Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCovsf( data = myFullData.sp$alive,
                                           data.dead = myFullData.sp$dead.recovery)
## Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]

## Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]



## ------   6. MAKE AUGMENTATION ------

## DATA ARRAYS
y.alive <- MakeAugmentation( y = y.ar.ALIVE,
                             aug.factor = DETECTIONS$aug.factor,
                             replace.value = 0)
y.dead <- MakeAugmentation( y = y.ar.DEAD,
                            aug.factor = DETECTIONS$aug.factor,
                            replace.value = 0)
y.aliveOthers <- MakeAugmentation( y = y.ar.ALIVEOthers, 
                                   aug.factor = DETECTIONS$aug.factor,
                                   replace.value = 0)
y.aliveStructured <- MakeAugmentation( y = y.ar.ALIVEStructured,
                                       aug.factor = DETECTIONS$aug.factor,
                                       replace.value = 0)

## INDIVIDUAL COVARIATES
already.detected <- MakeAugmentation( y = already.detected,
                                      aug.factor = DETECTIONS$aug.factor,
                                      replace.value = 0)



## ------   7. CONVERT TO CACHED DETECTORS AND SPARSE MATRIX ------

## ------     7.1. CREATE CACHED DETECTORS OBJECTS ------

maxDistReCalc <- 2.1*DETECTIONS$maxDetDist 

DetectorIndexLESS <- GetDetectorIndexLESS( 
  habitat.mx = habitat$habitat.mx,
  detectors.xy = nimData$detector.xy,
  maxDist = maxDistReCalc/res(habitat$habitat.r)[1],
  resizeFactor = 1,
  plot.check = TRUE)



## ------     7.2. TRANSFORM Y TO SPARSE MATRICES  ------

# STRUCTURED
SparseY <- GetSparseY(y.aliveStructured)

# OTHER
SparseYOth <- GetSparseY(y.aliveOthers)


## -----------------------------------------------------------------------------
## ------ III. MODEL SETTING & RUNNING ------- 

## ------   1. NIMBLE MODEL DEFINITION ------

modelCode <- nimbleCode({
  
  ##----- SPATIAL PROCESS ------ 
  dmean ~ dunif(0,100)
  lambda <- 1/dmean
  betaDens  ~ dnorm(0.0,0.01)
  ##
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
  }#i
  
  for(t in 2:n.years){
    for(i in 1:n.individuals){
      sxy[i, 1:2, t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
        upperCoords = upperHabCoords[1:numHabWindows, 1:2],
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
    gamma[t] ~ dunif(0,1)
    phi[t] ~ dunif(0,1)
    
    omega[1,1,t] <- 1-gamma[t]
    omega[1,2,t] <- gamma[t]
    omega[1,3,t] <- 0
    omega[2,1,t] <- 0
    omega[2,2,t] <- phi[t]
    omega[2,3,t] <- 1-phi[t]
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
  
  
  ##----- DETECTION PROCESS -----
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
      y.alive[i,1:nMaxDetectors,t] ~ dbinomLocal_normalCovsResponse( 
        detNums = nbDetections[i,t],
        detIndices = yDets[i,1:nMaxDetectors,t],
        size = trials[1:n.detectors],
        s = sxy[i,1:2,t],
        sigma = sigma[t],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = resizeFactor,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        lengthYCombined = maxNBDets,
        indicator = isAlive[i,t],
        p0State = p0[1:n.counties,t],
        trapCountries = detCounties[1:n.detectors],
        trapCovs = detCovs[1:n.detectors,t,1:n.covs],
        trapBetas = betaCovs[1:n.covs,t],
        responseCovs = detResponse[i,t],
        responseBetas = betaResponse[t])
      
      y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbinomLocal_normalCovsResponse(
        detNums = nbDetectionsOth[i,t],
        detIndices = yDetsOth[i,1:nMaxDetectorsOth,t],
        size = trials[1:n.detectors],
        s = sxy[i,1:2,t],
        sigma = sigma[t],
        trapCoords =  detector.xy[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cellsSparse,1:maxNBDets],
        localTrapsNum = nDetectorsLESS[1:n.cellsSparse],
        resizeFactor = resizeFactor,
        lengthYCombined = maxNBDets,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = isAlive[i,t],
        p0State = p0Oth[1:n.countries,t],
        trapCountries = detCountries[1:n.detectors,t],
        trapCovs = detCovsOth[1:n.detectors,t,1:n.covsOth],
        betaCov = betaCovsOth[1:n.covsOth,t],
        responseCovs = detResponse[i,t],
        responseCovs = betaResponseOth[t])
    }#i
  }#t
  
  ##----- DERIVED PARAMETERS ------
  for(t in 1:n.years){
    for(i in 1:n.individuals){ 
      isAlive[i,t] <- (z[i,t] == 2) 
    }#i
    N[t] <- sum(isAlive[1:n.individuals,t])
  }#t
  
})



## ------   2. NIMBLE CONSTANTS ------

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

# ADD TO NIMDATA
nimConstants$y.maxDet <- dim(DetectorIndexLESS$habitatID)[1]
nimConstants$x.maxDet <- dim(DetectorIndexLESS$habitatID)[2]
nimConstants$resizeFactor <- DetectorIndexLESS$resizeFactor
nimConstants$n.cellsSparse <- dim(DetectorIndexLESS$detectorIndex)[1]
nimConstants$maxNBDets <- DetectorIndexLESS$maxNBDets

nimConstants$nMaxDetectors <- SparseY$nMaxDetectors

nimConstants$nMaxDetectorsOth <- SparseYOth$nMaxDetectors


yearsNotSampled <- which(!years%in% yearsSampledNorrb)

nimConstants$countyToggle <- nimInits$p01
nimConstants$countyToggle[] <- 1
for(t in yearsNotSampled){
  nimConstants$countyToggle[1,t] <- 0
}#t

nimConstants$countryToggle <- nimInits$p01Oth
nimConstants$countryToggle[] <- 1
yearsNotSampled <- which(!years%in% yearsSampledNorrb)
for(t in yearsNotSampled){
  nimConstants$countryToggle[3,t] <- 0
}#t



## ------   3. NIMBLE INITS ------

## ------     3.1. GENERATE z DATA VALUES ------

z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, NA)

z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    zz[range.det[1]:range.det[2]] <- 2
  }
  return(zz)
}))



## ------     3.2. GENERATE z INITIAL VALUES ------

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

## LATENT VARIABLE DET RESPONSE
detResponse <- already.detected 
detResponse[rownames(detResponse) %in% "Augmented" ,1]  <- NA
InitsDetResponse <- detResponse
InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1,0.5)
InitsDetResponse[!is.na(detResponse)] <- NA



## ------     3.3. GENERATE sxy INITIAL VALUES ------

## sxy
## Create a data.frame with all detection of all Individuals detected
## project death to the next year
myData.deadProj <- myData.dead[,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1

## Remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year ), ]

AllDets <- rbind(myData.alive$myData.sp[,c("Id","Year")],
                 myData.deadProj[,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  coordsData = AllDetsxy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid =T )$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)

idAugmented <- which(rownames(z) %in%"Augmented")
Id.vector <-   y.ar$Id.vector

sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = nimData$lowerHabCoords,
                       upperCoords = nimData$upperHabCoords,
                       habitatGrid = nimData$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")

## ------   4. NIMBLE DATA ------

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

## ADD TO NIMDATA
nimData$detector.xy <- as.matrix(ScaledDetectors)          
nimData$lowerHabCoords <- as.matrix(ScaledLowerCoords)
nimData$upperHabCoords <- as.matrix(ScaledUpperCoords)

nimData$detectorIndex <- DetectorIndexLESS$detectorIndex
nimData$nDetectorsLESS <- DetectorIndexLESS$nDetectorsLESS
nimData$habitatIDDet <- DetectorIndexLESS$habitatID

## ADD TO NIMDATAx
nimData$y.alive <- SparseY$y 
nimData$yDets <- SparseY$yDets
nimData$nbDetections <- SparseY$nbDetections

# ADD TO NIMDATA
nimData$y.aliveOth <- SparseYOth$y 
nimData$yDetsOth <- SparseYOth$yDets
nimData$nbDetectionsOth <- SparseYOth$nbDetections

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

nimData$detCountries <- detCountriesNorb



## ------   5. NIMBLE PARAMETERS ------

nimParams <- c("N", 
               "lambda","dmean","betaDens",
               "omeg1","gamma","phi",
               "sigma","pResponse",
               "p0","betaCovs","betaResponse",
               "p0Oth","betaCovsOth","betaResponseOth")

nimParams2 <- c("z", "sxy")



## ------   6. LOOP THROUGH INITIAL VALUES & SAVE OBJECT ------
for(c in 1:4){

    ## LIST INITIAL VALUES
    nimInits <- list(
    "sxy" = sxy.init,
    "dmean" = runif(1,0,10),
    "z" = z.init,
    "omeg1" = c(0.5,0.5),
    "gamma" = runif(dim(y.alive)[3]-1,0,1),
    "p01" = array(runif(18,0,0.2),
                  c(nimConstants$n.counties,dim(y.alive)[3])),
    "p01Oth" = array(runif(18,0,0.2),
                     c(nimConstants$n.countries+1,dim(y.alive)[3])),
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
  
    
    nimInits$sxy <- round(nimInits$sxy, 5)
    
    
  # ## TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  # ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  # idDEtected <- which(!rownames(z) %in%"Augmented")
  # for(i in 1:length(idDEtected)){
  #   for(t in 1:nimConstants$n.years){
  #     if(!is.na(nimInits$sxy[i,1,t])){
  #       SXY <- nimInits$sxy[i, ,t]  
  #     } else {
  #       SXY <- nimData$sxy[i, ,t]
  #       }
  #     sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$resizeFactor)+1, trunc(SXY[1]/nimConstants$resizeFactor)+1]
  #     DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse, 
  #                                                       1:nimConstants$maxNBDets] 
  #     DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
  #     index <- DETECTIndexdetectorIndex [sxyID,1:DETECTLESS[sxyID]]
  #     
  #     ## GET NECESSARY INFO 
  #     n.detectors <- length(index)
  #     YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors,t]
  #     YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth,t]
  #     
  #     ## RECREATE Y
  #     if(nimData$nbDetections[i, t] > 0){
  #       for(j in 1:nimData$nbDetections[i, t]){
  #         ## check if a detection is out of the "detection window"
  #         if(sum(YDET[j] == index)==0){ 
  #           print(paste("id",i,"t",t,"j",j))
  #           }
  #       }
  #     }
  #   }
  #   }



  
  ## SAVE NIMBLE INPUT 
  save( nimData,
        nimConstants,
        y.dead,
        nimParams,
        nimParams2,
        modelCode,
        nimInits,
        file = file.path(WD,
                         modelName,
                         paste0(modelName,"Chain",c,".RData")))
}#c



## ------ 8. SAVE NECESSARY OBJECTS ------

save(habitat,
     detectors,
     studyArea,
     COUNTIES_AGGREGATEDSubset,
     myFilteredData.sp,
     myFullData.sp,
     COUNTIES_AGGREGATED,
     file = file.path(WD, modelName, "NecessaryObjects.RData"))

