##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE OPSCR analysis 2025 
##
## This R script reproduces the OPSCR analysis of the wolverine data as performed
## by RovQuant in 2024.
##
## NOTES : Two main updates from last year's script ('53.Cleaned2024.R')
##  1. Added a step to remove duplicated GPS tracks (from Asun)
##  2. Added dead recovery states ('recovered dead legal' and 'recovered dead other') 
##     as in the last wolf ('40.F_2024_sf.R') and bear analyses ('Bear_NORWAY_2015-2024.R').
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
library(spdep)
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
load(paste(dir.dropbox,"/DATA/MISC DATA/age.lookup.table.RData",sep=""))


## ------ SOURCE THE NIMBLE FUNCTION ------

source("C:/My_documents/RovQuant/Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R")


##------------------------------------------------------------------------------

## ------ 0.SET ANALYSIS CHARACTERISTICS -----

myVars <- list(
  ## WORKING DIRECTORY & MODEL NAME
  WD = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025",
  modelName = "Test2025",
  
  ## HABITAT SPECIFICATIONS
  HABITAT = list( countries =  c("SWE","NOR"),
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
if(!dir.exists(file.path(myVars$WD, myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}



##------------------------------------------------------------------------------
## ------ I.LOAD AND SELECT DATA ------

### ==== 1. HABITAT DATA ====

### ====    1.1.LOAD RAW SHAPEFILES ====

## POLYGONS OF THE REGION
GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) 
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP[GLOBALMAP$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- COUNTRIES %>%
  group_by(ISO) %>%
  summarize()

## POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp"))   
COMMUNES_SWE <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp"))   
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)

## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- COMMUNES %>%
  group_by(NAME_1) %>%
  summarize()

## AGGREGATE COUNTIES (OPTIONAL)
COUNTIES_AGGREGATE <- COUNTIES
#COUNTIES_AGGREGATED <- gSimplify(COUNTIES_AGGREGATED,tol=500, topologyPreserve = TRUE)
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
COUNTIES_AGGREGATED <- st_simplify(COUNTIES_AGGREGATE,preserveTopology = T,dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id
ggplot(COUNTIES_AGGREGATED) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))



### ====    1.2.CREATE STUDY AREA POLYGON ====

## CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
if(!is.null(myVars$HABITAT$countries)){
  myStudyArea <- COUNTRIES[COUNTRIES$ISO %in% myVars$HABITAT$countries, ]
  
  ## CREATE A POLYGON OF THE ACTUAL HABITAT POLYGON CONSIDERED (different from buffered.habitat.poly)
  myBufferedArea <- st_buffer(st_as_sf(myStudyArea) ,dist = myVars$HABITAT$habBuffer)
  myBufferedArea$id <- 1
  myBufferedArea <- myBufferedArea %>% group_by(id) %>% summarize()
  myBufferedArea <- st_intersection(myBufferedArea, GLOBALMAP)
}

myStudyArea$id <- 1
myStudyArea <- myStudyArea %>% group_by(id) %>% summarize()

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myBufferedArea), add = TRUE, col = rgb(0.72,0.14,0.14,0.3))
  plot(st_geometry(myStudyArea), add = TRUE, col ="red")
}



### ==== 2. NGS DATA ====

### ====    2.1.LOAD ROVBASE FILES ====

## [[ PD : UPDATE WITH NEW FILE NAMES!!!! ]]

DNA <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dna_wolverines.csv"),
                 fileEncoding = "latin1") ## NGS data from RovBase#[CM update to 20190627]

DEAD <- read.csv( file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/dead_carnivores.csv"),
                  fileEncoding = "latin1") ## Dead Recoveries from RovBase#[CM update to 20190627]

SUSPECT_NGS_SAMPLES <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/Remove ngs samples list wolverine 2024.csv"),
                                 fileEncoding = "latin1") ## DNA samples to be removed from Henrik

SUSPECT_DeadRecoSAMPLES <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/Remove dead recoveries list wolverine 2024.csv"),
                                     fileEncoding = "latin1") ## DNA samples to be removed from Henrik

HairTrapSamples <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/hairtrapsNB2024.xlsx")) ## DNA samples to be removed from Henrik

DEN <- read.csv( file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/DEN_COUNTS_2009_2024_fromHB.csv"),
                 fileEncoding = "latin1")
skandObs <- read_xlsx(file.path(dir.dropbox, "DATA/Skandobs/RB_Skandobs_2012_2024/Richard_Bischof_Skandobs_2012_2024dd.xlsx"))

rovbaseObs1 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB2810202415264376.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152348493.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152447860.xlsx"))
rovbaseObs4 <- read_xlsx(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/ALL SPECIES IN SEPERATE YEARS/RIB28102024152538742.xlsx"))

rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3,rovbaseObs4)



### ====    2.2.TRANSLATE SCANDINAVIAN CHARACTERS ====

colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA), dir.translation = dir.analysis )
#drop a column that makes cleanDataNew to fail
DNA <- DNA[,-which(colnames(DNA)%in% "Kjoenn..Individ.")]
colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD), dir.translation = dir.analysis )
colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN), dir.translation = dir.analysis)
colnames(skandObs) <- translateForeignCharacters(dat=colnames(skandObs), dir.translation = dir.analysis )
colnames(rovbaseObs) <- translateForeignCharacters(dat=colnames(rovbaseObs), dir.translation = dir.analysis )
rovbaseObs$Proevetype <- translateForeignCharacters(dat=rovbaseObs$Proevetype, dir.translation = dir.analysis )



## ==== 3. SEARCH EFFORT DATA ====

### ====    3.1.GPS SEARCH TRACKS ====

## LOAD GPS SEARCH TRACKS
TRACKS_SINGLE <- read_sf(file.path(dir.dropbox,
                                   "DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp"))
TRACKS_MULTI <- read_sf(file.path(dir.dropbox,
                                  "DATA/RovbaseData/ROVBASE DOWNLOAD 20250915/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp"))

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
## check
hist(ALL_TRACKS$Yr)

## SELECT TRACKS YEAR
myStudyArea.extent <- st_bbox(extent(myStudyArea))
st_crs(myStudyArea.extent) <- st_crs(COUNTRIES)

## REMOVE DUPLICATED TRACKS
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
  
  # try a fast way to identify duplicated tracks
  # turn to dataframe and identify them
  # sub_tracks_filter <- TRACKS_YEAR[[t]] %>%
  #   distinct(Dato, dist, .keep_all = T)
  # distinct(Person, Dato, dist, .keep_all = T)

}#t


## PLOT CHECK
if(myVars$plot.check){
  
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



### ====    3.2.DISTANCE TO ROADS ====

## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(file.path( dir.dropbox,
                                  "DATA/GISData/Roads/MinDistAllRoads1km.tif"))

## RASTERIZE DISTANCE TO ROADS
r <- fasterize(st_as_sf(myStudyArea), DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea)

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow= c(1,1))
  plot(DistAllRoads)
  plot(st_geometry(myStudyArea), add = T)
}



### ====    3.3.DAYS OF SNOW ====

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
SNOW <- stack(file.path( dir.dropbox, 
                         "DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2008_2024_Wolf.tif"))
## RENAME THE LAYERS
names(SNOW) <- paste(2008:2023,(2008:2023)+1, sep="_")
## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep="")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))



### ====    3.4. LOAD SCANDINAVIAN 20KM HABITAT ====

load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/Habitat20km.RData"))
load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData"))

plot(habitatRasters,"Habitat")
plot(st_geometry(myStudyArea),add=T)



##------------------------------------------------------------------------------
## ------ II.CREATE SCR DATA ------

### ==== 1. CLEAN AND FILTER NGS DATA ====

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ]
dim(DNA)

## Remove un-verified dead recoveries [HB]
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]



### ====    1.1. CLEAN NGS & DEAD RECOVERY DATA ====

myCleanedData.sp <- CleanDataNew2sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = myVars$DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(myVars$DATA$samplingMonths)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)

dim(myCleanedData.sp)

## PLOT CHECK
if(myVars$plot.check){
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myStudyArea),
       add = T, col ="red")
  plot( st_geometry(myCleanedData.sp),
        add = TRUE, pch = 19, cex = 0.2, col = "blue")
}



### ====    1.2. FILTER DATA ====

myFullData.sp <- FilterDatasf(
  myData = myCleanedData.sp,
  poly = myStudyArea,
  dead.recovery = T ,
  sex = c("Hann","Hunn"), # do the sex selection at the last moment
  setSex = T)

if(myVars$plot.check){
  plot(st_geometry(myFullData.sp$alive), add=TRUE,pch=19,cex=0.2, col="white")
}

## REMOVE SUSPECT SAMPLES ACCORDING TO HENRIK
myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)
dim(myFullData.sp$alive)
dim(myFullData.sp$dead.recovery)

myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
dim(myFullData.sp$alive)
dim(myFullData.sp$dead.recovery)

## Remove individuals that died twice# [CM] TO BE CHECKED BECAUSE "length(IdDoubleDead) < 0" and so it was desactivated
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

##dead recoveries flagged by Henrik that should always be removed (email from the 18/12/2024)
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!myFullData.sp$dead.recovery$RovBaseId %in% c("M495994", "M524051", "M524052", "M524053"),]


## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
sum(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
      myFullData.sp$dead.recovery$Month > 2 &
      myFullData.sp$dead.recovery$Month < 12)

myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
                                                                    myFullData.sp$dead.recovery$Month > 2 &
                                                                    myFullData.sp$dead.recovery$Month < 12),]


# 2) remove individuals that have a weight >0 and <4 between March and November
# format the weight correctly
myFullData.sp$dead.recovery$Helvekt <- as.character(myFullData.sp$dead.recovery$Helvekt)
myFullData.sp$dead.recovery$Slaktevekt <- as.character(myFullData.sp$dead.recovery$Slaktevekt)

#convert to decimals
myFullData.sp$dead.recovery$Helvekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Helvekt))
myFullData.sp$dead.recovery$Slaktevekt <- as.numeric(gsub(",", ".", myFullData.sp$dead.recovery$Slaktevekt))
# get the two weight columns together.
myFullData.sp$dead.recovery$weight <- ifelse(!is.na(myFullData.sp$dead.recovery$Helvekt),
                                             myFullData.sp$dead.recovery$Helvekt,
                                             myFullData.sp$dead.recovery$Slaktevekt)
#assign negative values to nas to avoid issues
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999
dim(myFullData.sp$dead.recovery)

# check with Henrik
# this step does not remove dead recoveries on id with weight==0 should it?
# WEIGTH DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$weight[(myFullData.sp$dead.recovery$weight >-1) &
                                            myFullData.sp$dead.recovery$Month%in% t],breaks=c(0:30), main=t,xlab="Weight")
}
# AGE DISTRIBUTION
par(mfrow=c(4,3))
for(t in 1:12){
  hist(myFullData.sp$dead.recovery$Age[(myFullData.sp$dead.recovery$Age >-1) &
                                         myFullData.sp$dead.recovery$Month%in% t],breaks=seq(-0.01,20.99,by=1),
       main=t,xlab="Age")
}

# check how many dead reco we remove and remove if more than 0
if(sum(myFullData.sp$dead.recovery$weight > 0 &
       myFullData.sp$dead.recovery$weight < 4 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
                                                                      myFullData.sp$dead.recovery$weight < 4 &
                                                                      myFullData.sp$dead.recovery$Month < 12 &
                                                                      myFullData.sp$dead.recovery$Month > 2),]
}
# check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2)>0){
  myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
                                myFullData.sp$dead.recovery$Month < 12 &
                                myFullData.sp$dead.recovery$Month > 2,  ]
}

table(myFullData.sp$alive$Year)



### ====    1.3. FILTER NGS & DEAD RECOVERY DATA ====
myFilteredData.sp <- myFullData.sp
dim(myFullData.sp$alive)

## Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]
table(myFilteredData.sp$alive$Year)
dim(myFilteredData.sp$alive)

myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]
table(myFilteredData.sp$dead.recovery$Year)
dim(myFilteredData.sp$dead.recovery)

## Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in% unlist(myVars$DATA$samplingMonths), ]
table(myFilteredData.sp$alive$Year)



### ====    1.4. SUBSET DETECTIONS IN NORRBOTTEN IN ALL YEARS EXCEPT 2017, 2018 and 2019 ==== 

COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018,2023)
is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))
sum(is.Norr, na.rm = T)

# check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(is.Norr)), ]$Year) %>% sum()

# subset
myFilteredData.sp$alive <- myFilteredData.sp$alive[- which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                                             !is.na(is.Norr)), ]
dim(myFilteredData.sp$alive)
table(myFilteredData.sp$alive$Year)

#plot check
for(t in 1:nYears){
  plot(st_geometry(myStudyArea))
  plot(st_geometry(COUNTIESNorrbotten),add=T,col="blue")
  plot(st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],]), col="red",add=T,pch=16)
}

# ## SELECT THE SEX
# #MYFULLDATA
# myFullData.sp <- myFullData.sp
# myFullData.sp$alive <- myFullData.sp$alive[myFullData.sp$alive$Sex %in% myVars$DATA$sex,]
# myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Sex %in% myVars$DATA$sex,]
# #myFilteredData
myFilteredData.spAllSex <- myFilteredData.sp

# myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Sex %in% myVars$DATA$sex,]
# myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Sex %in% myVars$DATA$sex,]
# 
# ## PLOT CHECK
# if(myVars$plot.check){
#    par(mfrow = c(1,3))
#    for(t in 1:nYears){
#       ## DEAD RECOVERIES TOTAL
#       tempTotal <-  myFilteredData.sp$dead.recovery[ myFilteredData.sp$dead.recovery$Year == years[t], ]
#       NGS_TabTotal <- table(tempTotal$Country)
#       ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x) sum(x>0))
#       ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
#       ## PLOT NGS SAMPLES
#       plot(st_geometry(GLOBALMAP), col="gray80")
#       plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
#       plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
#       plot(st_geometry(tempTotal), pch = 21, bg = "darkred",add=T)
#       # ## ADD NUMBER OF NGS samples and IDs per COUNTRY
#       graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
#       graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
#       ## ADD OVERALL NUMBERS
#       mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
#     }#t
# }



### ====    1.5. SEPARATE STRUCTURED AND OPPORTUNISTIC SAMPLING ====

### ====      1.5.1. ASSIGN SAMPLES TO TRACKS ====

# ## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
# myFilteredData.sp$alive$TrackRovbsID <- NA
# myFilteredData.sp$alive$TrackDist <- NA
# 
# TRACKSSimple_sf <- list()
# for(t in 1:nYears){
#   TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
#   TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbsID)
#   TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]#st_simplify(TRACKS_YEAR[[t]],preserveTopology = T,dTolerance = 1000)
# }
# 
# ## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
# dnatemp <- st_as_sf(myFilteredData.sp$alive)
# ## CREATE A BUFFER AROUND EACH DETECTION
# tmp <-  st_buffer(dnatemp, dist=750)
# for(i in 1:nrow(myFilteredData.sp$alive)){
#    # INTERSECT POINT WITH TRACKS,
#    t <- which(years %in% tmp[i,]$Year)
# 
#    whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
# 
#    tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])
# 
#    if(nrow(tmpTRACKS)==0){next}
#    # FIND THE CLOSEST TRACK
#    dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
#    #whichSameDate <- numeric()
#    # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
#    #whichSameDate <- which(as.character(tmpTRACKS$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
#    # IF NO MATCHING DATE ASSIGN TO NA?
#    if(length(dist)==0){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- NA
#      myFilteredData.sp$alive$TrackDist[i] <- NA
#    }
#    # IF  MATCHING DATE ASSING TO THAT TRACK
#    if(length(dist)==1){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
#      myFilteredData.sp$alive$TrackDist[i] <- dist
#    }
#    # IF SEVERAL MATCHING DATES ASSING TO THE CLOSEST OF THE MATCHING TRACKS
#    if(length(dist)>1){
#      myFilteredData.sp$alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
#      myFilteredData.sp$alive$TrackDist[i] <- min(dist)
#    }
#    print(i)
#    #if(is.na(myFilteredData.sp$alive$TrackRovbsID[i])){print(i)}
# }
# 
# 
# # SAVE THE FOR FASTER LOADING
# save(myFilteredData.sp, file = file.path(myVars$WD, myVars$modelName,
#                                          paste("_myFilteredData.sp", ".RData", sep = "")))
# load(file.path(myVars$WD, myVars$modelName,paste("_myFilteredData.sp", ".RData", sep = "")))
# save( myFilteredData.sp,
#       file ="C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test/myFilteredData_original.RData")
load("C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test/myFilteredData_original.RData")



### ====      1.5.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC & STRUCTURED ====

distanceThreshold <- 500

#Proeveleverandoer columns was replaced by two columns, merging them now...
myFilteredData.sp$alive$Proeveleverandoer <- ifelse(myFilteredData.sp$alive$Annen.innsamler...Rolle %in% "" , 
                                                    myFilteredData.sp$alive$Samlet.selv...Rolle,
                                                    myFilteredData.sp$alive$Annen.innsamler...Rolle)
table(myFilteredData.sp$alive$Proeveleverandoer, useNA = "always")
sum(duplicated(myFilteredData.sp$alive$RovbaseID))
sum(duplicated(myFilteredData.sp$alive$DNAID))


whichStructured <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold
table(whichStructured, useNA = "always")

myFilteredData.spStructured <- myFilteredData.sp$alive[whichStructured, ]
myFilteredData.spOthers <- myFilteredData.sp$alive[!whichStructured, ]

## CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
nrow(myFilteredData.spStructured) + nrow(myFilteredData.spOthers)
nrow(myFilteredData.sp$alive)

## PUSH ALL RESULTS FROM HAIR TRAPS TO OPPORTUNISTIC
whichHair <- which(myFilteredData.sp$alive$DNAID%in% HairTrapSamples$DNAID)

## Plot check 
plot(COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]$geometry)
plot( myFilteredData.sp$alive[whichHair,]$geometry,
      add = T, col = "red", pch = 16)
if(length(which(myFilteredData.spStructured$alive$DNAID%in% HairTrapSamples$DNAID))>0){
  print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
}



### ====      1.5.3. PLOT CHECKS ====

##pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"DetectionsStructuredOppBarplot",".pdf", sep="" ))
par(mfrow=c(2,1),mar=c(4,4,3,2))
barplot(rbind(table(myFilteredData.spStructured$Year),
              table(myFilteredData.spOthers$Year)),beside=T,ylim=c(0,2000),col=c(grey(0.2),grey(0.8)),ylab="Number of samples")
abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
title(main="500m threshold")
legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))

##
distanceThreshold1 <- 2000
whichStructured2000 <- myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold1
myFilteredData.spStructured2000 <- myFilteredData.sp$alive[whichStructured2000,]
myFilteredData.spOthers2000 <- myFilteredData.sp$alive[!whichStructured2000,]

barplot(rbind(table(myFilteredData.spStructured2000$Year),
              table(myFilteredData.spOthers2000$Year)),beside=T,ylim=c(0,2000),col=c(grey(0.2),grey(0.8)),ylab="Number of samples")
abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
title(main="2000m threshold")
legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))
#dev.off()

## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Proeveleverandoer %in% 
                                 c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
tab <- table(tmp$Year,tmp$TrackRovbsID,useNA ="always" )


## plot check 
# pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"DetectionsStructuredOpp",".pdf", sep="" ))
for(t in 1:nYears){
  par(mar = c(0,0,3,0), mfrow = c(1,3))
  tmp1 <- tmp[tmp$Year%in% years[t],]
  tmpNoTracks <- tmp1[is.na(tmp1$TrackRovbsID), ]
  tmpTracks <- tmp1[!is.na(tmp1$TrackRovbsID), ]
  plot( st_geometry(myStudyArea), col = "gray60",  main="Structured with track")
  plot(st_geometry(tmpTracks), pch=21, col="black", cex=1,bg="red",add=T)
  plot( st_geometry(myStudyArea), col = "gray60",  main="Structured without track")
  plot(st_geometry(tmpNoTracks), pch=21, col="black", cex=1,bg="blue",add=T)
  
  tmpOpp <- myFilteredData.sp$alive[!myFilteredData.sp$alive$Proeveleverandoer %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
  tmpOpp <- tmpOpp[tmpOpp$Year%in% years[t],]
  plot( st_geometry(myStudyArea), col = "gray60",  main="Other samples")
  plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
  mtext(years[t],adj = -0.8,padj = 1)
}
barplot(tab[ ,which(is.na(colnames(tab)))]/rowSums(tab),
        main = "% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 
# dev.off()

### plot check 
# pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"OverallDetectionsDeadRecoveries",".pdf", sep="" ))
plot(st_geometry(GLOBALMAP))
plot(st_geometry(myStudyArea),add=T)
plot(st_geometry(myFullData.sp$alive),pch=16, col="red", cex=0.3,add=T)
plot(st_geometry(myFullData.sp$dead.recovery),pch=16, col="blue", cex=0.3,add=T)
mtext(paste("Live detections", nrow(myFullData.sp$alive),
            "; ID:", nrow(unique(myFullData.sp$alive$Id))),
      line = +1)
mtext(paste("Dead recovery:", nrow(myFullData.sp$dead.recovery)))
# dev.off()



### ====    1.6. SEPARATE MORTALITY CAUSES ==== 

## MORTALITY CAUSES
length(myFilteredData.sp$dead.recovery)
length(myFullData.sp$dead.recovery)

MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))

## DEFINE LEGAL MORTALITY
# legalCauses <- MortalityNames[grep("JF", MortalityNames)]
# legalCauses <- c(legalCauses, MortalityNames[grep("Rifle", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("18", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("17", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Uspesifisert", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Fellefangst", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Hagle", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Skyddsjakt", MortalityNames)])
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])

## SPLIT MORTALITY CAUSES
legal.death <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]
Other.death <- myFilteredData.sp$dead.recovery[!myFilteredData.sp$dead.recovery$DeathCause %in% legalCauses, ]

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,3))
  for(t in 1:nYears){
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
    plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
    # points(tempTotal, pch = 21, bg = "darkred")
    plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
    ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
    ## ADD OVERALL NUMBERS
    mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
    mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
    # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
  }#t
}


## PLOT TREND DETECTIONS AND DEAD RECOVERIES OVER TIME AND SPACE #[CM]
## DETECTIONS
#pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_TRENDDetections.pdf",sep="")))
temp <- unique(myFilteredData.sp$alive[,c("Year","Country","DNAID")])
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## ID DETECTED
temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country,myFilteredData.sp$alive$Id)
tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
country.colors <- c("goldenrod1","goldenrod3")
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
lines(tab_Country.Year1[,"N"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year1[,"S"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## Average number of detection per detected ID  
tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
     ylab="Average Number of detections", xlab="Years")
lines(tab_Country.Year2[,"N"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year2[,"S"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("bottomright",c("N","S"), fill=country.colors)

## deadrecovery 
temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")])
tab_Country.Year <- table(temp$Year, temp$Country)
country.colors <- c("goldenrod1","goldenrod3")

par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
     ylab="N Id Dead recovered", xlab="Years")
lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
legend("topright",c("N","S"), fill=country.colors)
dev.off()



### ==== 2. GENERATE HABITAT ====

### ====    2.1. REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ====

## DELINEATE A BUFFER AROUND ALL DEADRECO 
BuffDead <- st_buffer(myFilteredData.spAllSex$dead.recovery, dist = myVars$HABITAT$habBuffer)
BuffDead$idd <- 1
BuffDead <- BuffDead %>% group_by(idd) %>% summarize()

## DELINEATE A BUFFER AROUND ALL DETECTIONS 
BuffAlive <- st_buffer(myFilteredData.spAllSex$alive, dist = myVars$HABITAT$habBuffer*1.4)
BuffAlive$idd <- 1
myBufferedArea <- BuffAlive %>% group_by(idd) %>% summarize()

## CUT TO SWEDISH AND NORWEGIAN BORDERS
myStudyArea <- st_intersection(myBufferedArea, myStudyArea)#[myStudyArea$ISO %in% c("SWE","NOR"), ])

# plot check 
par(mar=c(0,0,0,0))
plot(st_geometry(myStudyArea))
plot(st_geometry(myFilteredData.spAllSex$alive), pch=21, bg="red", cex=0.5,add=T)
plot(st_geometry(myFilteredData.spAllSex$dead.recovery), pch=21, bg="blue", cex=0.5,add=T)
plot(st_geometry(BuffDead),add=T, border="blue")
plot(st_geometry(myBufferedArea), border="red", add=T)
plot(st_geometry(myStudyArea), border="grey", add=T)



### ====    2.2. GENERATE HABITAT CHARACTERISTICS FROM THE NEW HABITAT DEFINITION ====

myHabitat <- MakeHabitatFromRastersf( 
  poly = myStudyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = myVars$HABITAT$habBuffer,                               
  plot.check = T)

## PLOT CHECK 
par(mfrow = c(1,2))
plot(myHabitat$habitat.r,legend = F)
plot(st_geometry(myStudyArea), add = T)

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

scaledHabGridCenters <- scaledHabGridCenters[myHabitat$habitat.r[]==1, ]
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}
image(habIDCells.mx)



### ====    2.3. SUBSET DETECTIONS BASED ON HABITAT EXTENT ==== 

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

## PLOT CHECK
if(myVars$plot.check){
  # par(mfrow = c(1,2))#[CM]
  plot(myHabitat$habitat.r)
  plot(st_geometry(myStudyArea), add = T, col = rgb(150/250,150/250,150/250, alpha = 0.75))
  plot(st_geometry(GLOBALMAP), add = T)
  plot(st_geometry(myHabitat$buffered.habitat.poly), add=T)
  plot(st_geometry(myFilteredData.sp$alive),pch=21, bg="red", cex=0.5,add=T)
  plot(st_geometry(myFilteredData.sp$dead.recovery),pch=21, bg="blue", cex=0.5,add=T)
}

## check correlation number of detetions ~ between monitoring season
myFilteredData.sp$dead.recovery$Id <- as.character(  myFilteredData.sp$dead.recovery$Id)
myFilteredData.sp$alive$Id <- as.character( myFilteredData.sp$alive$Id)
deadID <- unique(myFilteredData.sp$dead.recovery$Id)
ndet <- NULL
timeDiff <- NULL
i=36
for(i in 1:length(deadID)){
  tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Year
  
  timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Date-
    as.POSIXct(strptime(paste("01-12",tmpYear,sep="-"), "%d-%m-%Y")) 
  
  ndet[i] <- length(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] & 
                                              myFilteredData.sp$alive$Year %in% tmpYear,])
}

## PLOT CHECK
##pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Prop id deteced_Time available.pdf",sep="")))
plot(ndet ~ timeDiff,
     ylab = "Total number of detections",
     xlab = "Number of days between dec 1 and dead recovery")
hh <- hist(timeDiff[ndet>0], breaks = seq(0,400, by = 25))
hh1 <- hist(timeDiff[ndet == 0], breaks = seq(0,400, by = 25))
barplot(rbind(hh$counts/(hh$counts+hh1$counts),
              hh1$counts/(hh$counts+hh1$counts)),
        names.arg = hh$breaks[1:(length(hh$breaks)-1)],
        xlab = "number of days between dead reco and start monitoring",
        ylab = "%")
legend("topright",
       fill = c(grey(0.2),grey(0.8)),
       legend = c("detected","notDetected"))
dev.off()



### ====    2.4. GENERATE HABITAT-LEVEL COVARIATES ====

### ====       2.4.1. DEN COUNTS ====

DEN.sp <-  st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id  <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[,("id")]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD(as(DEN.sp,"Spatial"), h = 30000,#H[h],#[CM]
             grid = as(myHabitat$habitat.r, 'SpatialPixels'))
  )
)

if(myVars$plot.check){
  plot(DEN.r)#,main=paste("Sum Den-2013:2018Kernel H=", H[h], sep=""))
  plot(st_geometry(myStudyArea), add = TRUE, border = "black")
}

##EXTRACT COVARIATEs
denCounts <- DEN.r[myHabitat$habitat.r[ ]==1]
denCounts <- round(scale(denCounts),digits = 2)



### ==== 3. GENERATE DETECTORS ====

### ====    3.1. GENERATE DETECTORS CHARACTERISTICS ====

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
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

plot(st_geometry(myDetectors$main.detector.sp))
plot(st_geometry(myDetectors$main.detector.sp))

## IDENTIFY DETECTORS IN NORBOTTEN 
COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten","Troms","Västerbotten",
                                                            "Nordland","Finnmark"),]
COUNTIESAroundNorrbotten <- st_simplify(COUNTIESAroundNorrbotten, dTolerance = 500)

## CREATE A NORROBOTTEN DETECTOR GRID
distDestsCounties <- st_distance(myDetectors$main.detector.sp, COUNTIESAroundNorrbotten,byid = T)

# distDestsCounties <- gDistance(myDetectors$main.detector.sp, COUNTIESAroundNorrbotten,byid = T)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min)==3)
## PLOT CHECK 
plot(st_geometry(COUNTIESAroundNorrbotten))
plot(st_geometry(myDetectors$main.detector.sp),
     col = "black", pch = 16, cex = 0.3, add = T)
plot(st_geometry(myDetectors$main.detector.sp[detsNorrbotten,]),
     col = "red", pch = 16, cex = 0.3, add = T)

## RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * myVars$DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * myVars$DETECTORS$detResolution
lowerDetCoords.dead <- detector.dead.xy - 0.5 * myVars$DETECTORS$detDeadResolution
upperDetCoords.dead <- detector.dead.xy + 0.5 * myVars$DETECTORS$detDeadResolution

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,2))
  ## PLOT NGS DETECTORS
  plot( st_geometry(myHabitat$buffered.habitat.poly),
        main = paste(n.detectors, "Detectors Alive"),
        col = rgb(0.16,0.67,0.16, alpha = 0.3))  
  plot( st_geometry(myStudyArea),
        col = rgb(0.16,0.67,0.16, alpha = 0.5), add = TRUE)
  plot( st_geometry(myDetectors$main.detector.sp),
        col = "red", pch = 16, cex = 0.1, add = TRUE)
  plot( st_geometry(COUNTRIES), add = TRUE)
  ## PLOT DEAD DETECTORS
  plot( st_geometry(myHabitat$buffered.habitat.poly),
        main = paste(n.detectors.dead, "Detectors Dead"),
        col = rgb(0.16,0.67,0.16, alpha = 0.3)) 
  plot( st_geometry(myStudyArea),
        add = T, col = rgb(0.16,0.67,0.16, alpha = 0.5))
  plot( st_geometry(myDetectors.dead$main.detector.sp),
        col = "red", pch = 16, cex = 0.1, add = TRUE)
  plot( st_geometry(COUNTRIES), add = TRUE)
}



### ====    3.2. GENERATE DETECTOR-LEVEL COVARIATES ====

### ====      3.2.1. EXTRACT COUNTRIES ====

dist <- st_distance(myDetectors$main.detector.sp, COUNTRIES, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))

## PLOT CHECK 
if(myVars$plot.check){
  par(mfrow = c(1,2))
  myCol <- c("blue4", "yellow1")
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Countries")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(myDetectors$main.detector.sp), col = myCol[detCountries], pch = 16, cex = 0.8, add=T)
  plot(st_geometry(COUNTRIES), add = TRUE)
}



### ====      3.2.2. EXTRACT COUNTIES ====

## ASSIGN COUNTIES TO DETECTORS
dist <- st_distance(myDetectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F )
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties),]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))

## PLOT CHECK 
if(myVars$plot.check){
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
  # text(st_geometry(COUNTIES_AGGREGATED), labels = COUNTIES_AGGREGATED$id, col = "black")  
  plot( st_geometry(myDetectors$main.detector.sp[detCounties %in% 3, ]),
        col = "red", pch = 16, cex = 0.8, add = T)
}



### ====      3.2.3. EXTRACT GPS TRACKS LENGTHS ====

## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp),
                                      rep(1,nrow(myDetectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)
plot(st_geometry(detectorGrid))


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
  #transect_qi = mean(QI))          ## Get mean transects quality index for each detector grid cell
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
plot(rr)

sum(st_length(TRACKS_YEAR[[t]]))/1000
sum(TRACKS.r[[t]][],na.rm=T)/1000


#pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Tracks.pdf",sep="")))
if(myVars$plot.check){
  NORTRACKS <- SWETRACKS <- 0
  # par(mfrow = c(2,2))#[CM]
  for(t in 1:nYears){
    plot( TRACKS.r[[t]],main=years[t], breaks=cuts, col = col,legend=FALSE)#[CM]
    plot(st_geometry(myHabitat$habitat.poly), main = years[t],add=T)
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
#dev.off()



### ====      3.2.4. EXTRACT DISTANCES TO ROADS ====

## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)

## if NA returns the average value of the cells within 20000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract(DistAllRoads, myDetectors$main.detector.sp[isna,], buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp

if(myVars$plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Distance to roads")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(DistAllRoads,add=T)
  plot(st_geometry(myDetectors$main.detector.sp), cex=DoScale(detRoads), pch = 16, add = T)
}



### ====      3.2.5. EXTRACT DAYS OF SNOW ====

## EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(myDetectors$main.detector.sp)[1], ncol = nYears)

det.sptransf <- st_transform(myDetectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp
if(myVars$plot.check){
  plot(st_geometry(myDetectors$main.detector.sp),cex=DoScale(detSnow[,6],l = 0,u = 0.5),pch=16)
}



### ====      3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ====

### ====        3.2.6.1. SKANDOBS ====

## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))
## MAKE IT SPATIAL 
# coordinates(skandObs) <- skandObs[,c("longitude","latitude")]
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")

skandObs <- st_transform(skandObs, st_crs(myStudyArea))
## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist(myVars$DATA$samplingMonths))
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) #--- need to change for other species
skandObs <- skandObs[subset,] 
## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(myHabitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1,]

subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace,] 
plot(st_geometry(habitat.rWthBufferPol))
plot(st_geometry(skandObs),col="red",add=T)
## PLOT CHECK 
## SUMMARY SKANDOBS
##pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"skandObs",".pdf", sep="" ), width = 10)
barplot(table(skandObs$monitoring.season ))
barplot(table(skandObs$month ), xlab="Months")
# barplot(table(skandObs$activity),cex.names=0.7)
barplot(table(skandObs$species))
## MAPS 
par(mar=c(0,0,2,0))
for(t in 1:nYears){
  plot(st_geometry(myStudyArea), main= years[t])
  plot(st_geometry(skandObs[skandObs$monitoring.season %in% years[t],  ]), pch=16, col="red", cex=0.1,add=T)
}
dev.off()

## RASTERIZE AT THE DETECTOR LEVEL
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
plot(r.skandObsSamplesBinary[[t]])



### ====        3.2.6.2. ROVBASE ====

## GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Nord (UTM33/SWEREF99 TM)`),]
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
            "Saliv/Spytt", "Loepeblod", "Vev"),
  month = unlist(myVars$DATA$samplingMonths))

## SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset, ] 

## SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
subset <- (rovbaseObs.sp$`Art (Analyse)` %in% filter$species) & !is.na(rovbaseObs.sp$Individ) 
rovbaseObs.sp <- rovbaseObs.sp[!subset, ] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace, ] 

## Check if the filter is correct  
tmp <- rovbaseObs.sp[!is.na(rovbaseObs.sp$Individ), ]
table(tmp$`Art (Analyse)`) # Correct, all wolverine samples have NAs (not IDs)

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

r.OtherSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.OtherSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))

## PLOT
#pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"mapStructuredOthers",".pdf", sep="" ))
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
#dev.off()



### ====        3.2.6.3. COMBINE ROVBASE AND SKANDOBS ====

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <-    1
}

for(t in 1:nYears){
  par(mfrow=c(1,3),mar=c(0,0,5,0))
  plot(r.OtherSamplesBinary[[t]],main=years[t])
  plot(r.skandObsSamplesBinary[[t]])
  plot(r.SkandObsOtherSamplesBinary[[t]])
}  



### ====        3.2.6.4. SMOOTH THE BINARY MAP ====

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
par(mfrow=c(1,3))
plot(r.SkandObsOtherSamplesBinary[[t]], main="Raw Binary",axes=F,box=F)
plot(ds.brick[[t]], main="Smoothed",axes=F,box=F)
plot(ds.brickCont[[t]], main="Binary after smoothing",axes=F,box=F)



### ====        3.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ====

#IDENTIFY HAIR SAMPLES
tmpHair <- myFilteredData.spAllSex$alive[which(myFilteredData.spAllSex$alive$DNAID%in% HairTrapSamples$DNAID),]
#MANNUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL. 
tmpyr <- unique(tmpHair$Year)
for( i in 1:length(tmpyr)){
  t <- which(years %in% tmpyr)
  whereHair <-raster::extract(r.SkandObsOtherSamplesBinary[[t]],tmpHair,cellnumbers=T)
  r.SkandObsOtherSamplesBinary[[t]][whereHair[,1]] <- 1
  plot(r.SkandObsOtherSamplesBinary[[t]])
  plot(tmpHair$geometry,add=T,col="red")
}



### ====        3.2.6.7. ASSIGN THE COVARIATE ====

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract(r.SkandObsOtherSamplesBinary, myDetectors$main.detector.sp)
colSums(detOtherSamples)



### ====      3.2.7. SCALE AND ROUND DETECTOR-LEVEL COVARIATES ====

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
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}

## PLOT CHECK 
tmp <- detectorGrid.r
par(mfrow=c(2,5),mar=c(0,0,0,0))
max <- max(detCovsOth[,,2])
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))
for(t in 1:nYears){
  plot(detectorGrid.r, col=c(grey(0.2),grey(0.8)),axes=F,legend=F,box=F,)
  
  tmp[!is.na( detectorGrid.r)] <- detCovsOth[,t,2]
  plot(tmp,axes=F,legend=F,box=F,breaks = cuts, col=col,add=T)
}



## PLOT CHECK
dev.off()

#pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_detections over space and time.pdf",sep="")))
if(myVars$plot.check){
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
}
#dev.off()



### ==== 4. GENERATE y DETECTION ARRAYS ====

### ====    4.1. GENERATE NGS DETECTIONS : y.alive[i,j,t] ====

## ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$alive,                
  myDetectors = myDetectors.dead$main.detector.sp,
  mysubDetectors = myDetectors.dead$detector.sp,
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

## MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET ASSIGNED TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING
## FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE OF NORRBOTTEN
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)
whichdets <- which(!myData.alive$myData.sp$Year %in% yearsSampledNorrb &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)
whichdetsStruc <- which(!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveStruc$myData.sp$Detector %in% detsNorrbotten)
whichdetsOther <- which(!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb &
                          myData.aliveOthers$myData.sp$Detector %in% detsNorrbotten)

## ASSIGN DETECTORS 
subdetector.sf <- myDetectors$detector.sp

i=11
## ALL
for(i in 1:length(whichdets)){
  tmp <- myData.alive$myData.sp[whichdets[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, myDetectors$main.detector.sp)
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
  dist <- st_distance(tmp, myDetectors$main.detector.sp)
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

##OTHER
for(i in 1:length(whichdetsOther)){
  tmp <- myData.aliveOthers$myData.sp[whichdetsOther[i],]
  # MAIN DETECTOR 
  dist <- st_distance(tmp, myDetectors$main.detector.sp)
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

## SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
sum(myData.aliveOthers$myData.sp$Detector[!myData.aliveOthers$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)
sum(myData.aliveStruc$myData.sp$Detector[!myData.aliveStruc$myData.sp$Year %in% yearsSampledNorrb] %in%detsNorrbotten)

## DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf(
  myData = myFilteredData.sp$dead.recovery,
  myDetectors = myDetectors.dead$main.detector.sp,
  radius = myVars$DETECTORS$detResolution)

# ##EXPORT THE DATA 
# if(myVars$DATA$sex=="Hann"){
#   assign("myFilteredData.spM", myFilteredData.sp)
#   assign("myFilteredData.spOthersM", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
#   
#   save(myFilteredData.spM, myFullData.spM,
#        myFilteredData.spOthersM,myFilteredData.spStructuredM,
#        file = file.path(myVars$WD, myVars$modelName,
#                         paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
#   
# }else{
#   assign("myFilteredData.spF", myFilteredData.sp)
#   assign("myFilteredData.spOthersF", myFilteredData.spOthers)
#   assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
#   
#   save(myFilteredData.spF, myFullData.spF,
#        myFilteredData.spOthersF,myFilteredData.spStructuredF,
#        file = file.path(myVars$WD, myVars$modelName,
#                         paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
#   
# }


### ====    4.3. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ====

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

## MAKE SURE THE Y HAVE THE SAME DIMENSIONS#
y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0

## FILL IN THE Y ARRAYS 
y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]],,] <-  y.ar.ALIVEOth
y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]],,] <-  y.ar.ALIVEStruc

sum(y.ar.ALIVEOthers[1307,,8])
sum(y.ar.ALIVEStructured[1307,,8])

## PROJECT THE DEATH TO THE NEXT OCCASION.
y.ar.DEADProjected <- y.ar$y.ar2 
y.ar.DEADProjected[] <- 0
for(t in 2:nYears){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}

y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])

dim(y.ar.DEAD)
y.ar.DEAD[y.ar.DEAD>0] <- 1
dim(y.ar$y.ar2)



### ====    4.3. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ====

distances <- list()
for(t in 1:nYears){
  print(paste("------ ", t ," -------", sep = "" ))
  distances[[t]] <- CheckDistanceDetectionsV2sf( y = y.ar.ALIVE[,,t], 
                                                 detector.xy = detector.xy, 
                                                 max.distance = myVars$DETECTIONS$maxDetDist,
                                                 method = "pairwise",
                                                 plot.check = F)
  
  # PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
  if(myVars$plot.check){
    par(mfrow = c(1,1))
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      count <- 0
      for(i in affected.ids){
        count <- count+1
        plot(st_geometry(myStudyArea), main = paste("t: ",t,"     i: ", names(affected.ids)[count], sep = ""))
        scalebar(2*myVars$DETECTIONS$maxDetDist, xy = c(800000,6700000), type = "bar", divs = 2, below = "km",
                 label = c(0, myVars$DETECTIONS$maxDetDist/1000, myVars$DETECTIONS$maxDetDist/500), cex = 0.8, adj = c(0.5,-0.9))
        plot(st_geometry(COUNTRIES), add = T)
        plot(st_geometry(myDetectors$main.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar.ALIVE)[[1]][i] &
                                         myFilteredData.sp$alive$Year == years[t], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
        arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2],
               length = 0.1, lwd = 1)
        plot(st_geometry(myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
        
        tmp2 <- myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
        plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
      }#i
    }#if
  }#if plot.check
  ## REMOVE DETECTIONS THAT ARE FURTHER THAN  THE THRESHOLD
  y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  ## REMOVE DETECTIONS ALSO IN MYDATA TO RUN GETSINITS
  
  # tmpmyData.sp <- myData.alive$myData.sp
  # 
  # for(t in 1:nYears){
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
    #       if(sum(which(!unique(myData.alive$myData.sp$Id) %in% unique(tmpmyData.sp$Id))>0)){
    #         print(idd[i])
    #       }
    #   }
    #   }
    
    
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]),]
  }
  # tmp <- myData.alive$myData.sp[(myData.alive$myData.sp$Id %in% idd[i] &
  #                                  myData.alive$myData.sp$Year %in% years[t]),]
  # tmp[(tmp$Id %in% idd[i]),]
  # myData.alive$myData.sp[(myData.alive$myData.sp$Id %in% idd[i]&
  # #                           myData.alive$myData.sp$Detector %in% detIds &  myData.alive$myData.sp$Year %in% years[t]),]
  # [1] "JI000075 Ind75"  
  # [2] "JI000544 Ind544" 
  # [3] "JI001005 Ind1006"
  # [4] "JI001093 Ind1106"
  # [5] "JI001095 Ind1108"
  # [6] "JI001097 Ind1110"
  # [7] "JI001163 Ind768" 
  # [8] "JI001194 Ind804" 
  # [9] "JI001202 Ind812"
  
  
  # [1] "id 586 t 10 j 1"
  # [1] "id 945 t 10 j 1"
  # [1] "id 1340 t 10 j 1"
  # [1] "id 1601 t 10 j 1"
  # [1] "id 1620 t 10 j 1"
  # [1] "id 1631 t 10 j 1"
  # [1] "id 1631 t 10 j 2"
  # [1] "id 1649 t 10 j 1"
  # [1] "id 1649 t 10 j 2"
  # [1] "id 1650 t 10 j 1"
  # [1] "id 1674 t 10 j 1"
  # [1] "id 1682 t 10 j 1"
  # [1] "id 1692 t 10 j 1"
  # [1] "id 586 t 10 j 1"
  # [1] "id 945 t 10 j 1"
  # [1] "id 1340 t 10 j 1"
  # [1] "id 1601 t 10 j 1"
  # [1] "id 1620 t 10 j 1"
  # [1] "id 1631 t 10 j 1"
  # [1] "id 1631 t 10 j 2"
  # [1] "id 1649 t 10 j 1"
  # [1] "id 1649 t 10 j 2"
  # [1] "id 1650 t 10 j 1"
  # [1] "id 1674 t 10 j 1"
  # [1] "id 1682 t 10 j 1"
  # [1] "id 1692 t 10 j 1
  
  #     if(length(which(myData.alive$myData.sp$Id %in% idd[i]))==0){
  #       print(paste(idd[i], t, sep="_"))
  #     }
  #   
}



### ====    4.4. GENERATE INDIVIDUAL-LEVEL COVARIATES ====

### ====      4.4.1. TRAP-RESPONSE ====

## Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCovsf(data=myFullData.sp$alive,
                                          data.dead=myFullData.sp$dead.recovery)
## Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]
## Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]

## Plot an image of the matrix
if(myVars$plot.check){
  par(mfrow = c(1,1))
  barplot(colSums(apply(y.ar.ALIVE, c(1,3),function(x)any(x>0))))
  barplot(colSums(already.detected), add = TRUE, col = "gray40")
  legend(x = 0, y = 250, legend = c("newly Det", "already Det"),
         fill = c("gray80", "gray40"))
}



### ====      4.4.2. AGE ====

min.age <- age <- precapture <- matrix(NA, dim(y.ar.ALIVE)[1], dim(y.ar.ALIVE)[3], dimnames = list(y.ar$Id.vector,years))

temp <- apply(y.ar.ALIVE, c(1,3), sum)
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
    if(birth.year<latest.recruitment.year) min.age[i,] <- years-birth.year 
  },silent = TRUE)
  
  try({
    birth.year <- this.set$Death - this.set$age
    age[i,] <- years-birth.year
  }, silent = TRUE)
}
image(t(min.age))
image(t(age))



### ==== 5. MAKE AUGMENTATION ====

## DATA ARRAYS
y.alive <- MakeAugmentation(y = y.ar.ALIVE, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.dead <- MakeAugmentation(y = y.ar.DEAD, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveOthers <- MakeAugmentation(y = y.ar.ALIVEOthers, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveStructured <- MakeAugmentation(y = y.ar.ALIVEStructured, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

## INDIVIDUAL COVARIATES
already.detected <- MakeAugmentation(y = already.detected, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
# dispersalToggle <- MakeAugmentation(y = dispersalToggle, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
age <- MakeAugmentation(y = age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
min.age <- MakeAugmentation(y = min.age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
precapture <- MakeAugmentation(y = precapture, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)



##------------------------------------------------------------------------------

## ------ III.MODEL SETTING AND RUNNING ------- 

### ==== 1. NIMBLE MODEL DEFINITION ====

modelCode <- nimbleCode({
  
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
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
      numGridRows =  y.max,
      numGridCols = x.max
      
    )
  }#i
  
  for(t in 2:n.years){
    for(i in 1:n.individuals){
      sxy[i, 1:2, t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows, 1:2]
        ,
        upperCoords = upperHabCoords[1:numHabWindows, 1:2]
        ,
        s = sxy[i, 1:2, t - 1]
        ,
        lambda = lambda
        ,
        baseIntensities = habIntensity[1:numHabWindows]
        ,
        habitatGrid =  habitatGrid[1:y.max,1:x.max]
        ,
        numGridRows = y.max
        ,
        numGridCols = x.max
        ,
        numWindows= numHabWindows
      )
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



### ==== 2. NIMBLE CONSTANTS ====

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
                      x.max = dim(habIDCells.mx)[2])#max(detCountries))



### ==== 3. NIMBLE INITS ====

### ====    3.1. GENERATE INITIAL z ====

z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, NA)
z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    zz[range.det[1]:range.det[2]] <- 2
  }
  return(zz)
}))



### ====    3.2. GENERATE z INITIAL values ====

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



### ==== 4. NIMBLE DATA ====

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



### ==== 5. NIMBLE PARAMETERS ====

nimParams <- c("N", 
               "omeg1", "gamma", "phi", "pResponse","lambda","dmean","p0Oth",
               "p0", "sigma", "betaDens", "betaCovs","betaResponse","betaCovsOth","betaResponseOth")

nimParams2 <- c("z", "sxy")



### ==== 6. CONVERT TO CACHED DETECTORS AND SPARSE MATRIX ====

### ====   6.1. RESCALE COORDINATES  ====

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



### ====   6.2. CREATE CACHED DETECTORS OBJECTS ====

## [CM] reduce multiplicator to 3 
maxDistReCalc <- 2.1 * myVars$DETECTIONS$maxDetDist 

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



### ====   6.3. TRANSFORM Y TO SPARSE MATRICES ====

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



### ==== 7. LOOP THROUGH INITIAL VALUES AND SAVE OBJECT ====

# sxy
#create a data.frame with all detection of all Individuals detected
#project death to the next year
myData.deadProj <- myData.dead[,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year +1#project dead reco to the next year
#remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year ),]

AllDets <- rbind(myData.alive$myData.sp[,c("Id","Year")],
                 myData.deadProj[,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(coordsData = AllDetsxy,
                                            coordsHabitatGridCenter = myHabitat$habitat.xy,
                                            scaleToGrid =T )$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)

idAugmented <- which(rownames(z) %in%"Augmented")
Id.vector <-   y.ar$Id.vector
lowerCoords = nimData$lowerHabCoords
upperCoords = nimData$upperHabCoords
habitatGrid = nimData$habitatGrid

# AllDetections[AllDetections$Id %in% Id.vector[586],]
# 
# i=945
# dimnames(y.aliveStructured)[[1]][i]#[1] "JI413571 Ind5343"
# Id.vector[i]
# tmp <- AllDetections[AllDetections$Id %in% Id.vector[945] & AllDetections$Year %in%c(2023),]
# plot(COUNTRIES$geometry)
# plot(tmp$geometry,add=T,col="red",pc=16)
# 
# par(mfrow=c(2,5),mar=c(0,0,0,0))
# for(t in 1:nYears){
#   plot(COUNTRIES$geometry)
#   
# points(sxy.initscaled[which(dimnames(sxy.init)[[1]]%in% Id.vector[i]),2,t]~
#          sxy.initscaled[which(dimnames(sxy.init)[[1]]%in% Id.vector[i]),1,t],col="red",pch=16)
# 
#   tmp <- AllDetections[AllDetections$Id %in% Id.vector[i] & AllDetections$Year %in%years[t],]
#   plot(tmp$geometry,add=T,col="blue",pch=16)
#   
#   
#   YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
#   YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
#   
#   points(detector.xy[YDET[1:nimData$nbDetections[i, t]],2]~
#            detector.xy[YDET[1:nimData$nbDetections[i, t]],1], col="green", pch=16)
#   points(detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],2]~
#            detector.xy[YDETOth[1:nimData$nbDetectionsOth[i, t]],1], col="purple", pch=16)
# }

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



### ==== 8. calculate realized phi ====

# get location of individuals 
sxy.initscaled <- scaleCoordsToHabitatGrid(coordsData = sxy.init,
                                           coordsHabitatGridCenter = myHabitat$habitat.xy,
                                           scaleToGrid =F )$coordsDataScaled

## initialize objects 
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

###
pdf(file = file.path( myVars$WD,
                      myVars$modelName,
                      paste0(myVars$modelName,"realizedPhiCountry.pdf")))

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

###
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


###
pdf(file = file.path( myVars$WD,
                      myVars$modelName,
                      paste0(myVars$modelName,"realizedPhiCounties.pdf")),
    width = 11, height = 6)
# PHI
#NORWAY
par(mfrow = c(1,2))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Norway")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(phiNOR)){
  points(phiNOR[,c]~yr,pch=16,type="b", col=c)
}
legend("bottomleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)

###SWEDEN
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z", main="Sweden")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(phiSWE)){
  points(phiSWE[,c]~yr,pch=16,type="b", col=c)
}
legend("bottomleft",colnames(phiSWE),col=c(1:ncol(phiSWE)),pch=16)


# RECRUITS
#NORWAY
par(mfrow = c(1,2))
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", 
     ylab = "Realized recruitment from z", main="Norway")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(recruitNOR)){
  points(recruitNOR[,c]~yr,pch=16,type="b", col=c)
}
legend("topleft",colnames(phiNOR),col=c(1:ncol(phiNOR)),pch=16)

###SWEDEN
plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized recruitment from z", main="Sweden")
axis(1, at = 1:(nYears-1) , labels = paste( years[1:(nYears-1)]+1, years[1:(nYears-1)]+2,sep="-"))
yr <- c(1:(nYears-1))
for(c in 1:ncol(recruitSWE)){
  points(recruitSWE[,c]~yr,pch=16,type="b", col=c)
}
legend("topleft",colnames(recruitSWE),col=c(1:ncol(phiSWE)),pch=16)
dev.off()

## prop detected vs Alive in z
propDet <- 0
for(t in 1:nYears){
  whichdets <- unique(c(which(nimData$nbDetections[,t]>0),
                        which(nimData$nbDetectionsOth[,t]>0)))
  whichAlive <- which(nimData$z[,t]%in%2)
  propDet[t] <- length(whichdets)/length(whichAlive)
}



### ==== 9. LIST NIMBLE INITS ====

for(c in 1:4){
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
  
  # SXY
  # nimInits$sxy <- UTMToGrid(data.sxy = sxy.init,
  #                           grid.sp = myHabitat$habitat.sp)$data.scaled.xy
  # dimnames(sxy.init)[[2]] <- c("x","y")
  # nimInits$sxy <- scaleCoordsToHabitatGrid(coordsData = sxy.init,
  #                                             coordsHabitatGridCenter = myHabitat$habitat.xy,
  #                                             scaleToGrid =F )$coordsDataScaled
  # 
  # nimData$sxy <- UTMToGrid(data.sxy = nimData$sxy,
  #                          grid.sp = myHabitat$habitat.sp)$data.scaled.xy
  
  ### TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK 
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  i=586
  t=10
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
  
  ## PLOT CHECK 
  nimInits$sxy <- round(nimInits$sxy, 5)#---an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
  
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
  
  
  
  ### ==== 7. SAVE NIMBLE INPUT ====
  
  save(nimData,
       nimConstants,
       y.dead,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(myVars$WD, myVars$modelName,
                        paste(myVars$modelName,"Chain", c, ".RData", sep = "")))
}#c



### ==== 10. SAVE NECESSARY OBJECTS ====

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
save(myHabitat.list, myDetectors, COUNTRIES, myStudyArea.poly, COMMUNES,COUNTIES_AGGREGATEDSubset,
     myFilteredData.sp, myFullData.sp, COUNTIES_AGGREGATED,
     file = file.path(myVars$WD, myVars$modelName, "NecessaryObjects.RData" ))



### ==== 11. NIMBLE RUN ====

### ====    11.1 CONFIGURE NIMBLE MODEL ====

for(c in 2:4){
  load(file.path(myVars$WD, myVars$modelName,
                 paste(myVars$modelName,"Chain", c, ".RData", sep = "")))
  ptm <- proc.time()
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        inits = nimInits,
                        data = nimData,
                        check = FALSE,
                        calculate = FALSE) 
  system.time(print(model$calculate()))#-740704.4
}#c

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
                                          detResponse = nimData$detResponse[i,t]
    )
    
    
  }
}
dbin_LESS_Cached_MultipleCovResponse


## MAKE THE SINGLE SEASON SCR MODEL ###

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
      numGridCols = x.max
      
    )
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



###
for(ch in 1:4){
  
  for(t in 1:nYears){   
    
    load(file.path(myVars$WD, myVars$modelName,
                   paste0(myVars$modelName,"Chain", ch, ".RData")))
    
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
    
    
    ##################################  
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
    # sigma
    nimInits$sigma <-  runif(1,1,2)
    # trapBetas
    nimInits$betaCovs <- nimInits$betaCovs[,t]
    nimInits$betaCovsOth <- nimInits$betaCovsOth[,t]
    
    # trapBetas
    #nimInits$detResponse <- nimInits$detResponse[t]
    nimInits$betaResponse <- nimInits$betaResponse[t]
    nimInits$betaResponseOth <- nimInits$betaResponseOth[t]
    
    # beta.dens
    #nimInits$betaDens <- nimInits$betaDens#[t]
    
    
    ## parameter for the latent detRespons covariate
    #nimInits$detRespone <-  runif(1,0.4,0.6)
    # nimInits$probAdult  <-  runif(1,0.4,0.6)
    nimData$detCountries <- nimData$detCountries[,t]
    
    
    # get the new number of indivudals 
    nimConstants$n.individuals <- nrow(nimInits$sxy)
    
    
    ##Nimble parameters
    nimParams <- c("N", "psi", "pResponse","p0Oth","betaCovsOth","betaResponseOth",
                   "p0", "sigma", "betaDens", "betaCovs","betaResponse","betaResponseOth")
    
    modelCode <- modelCode1
    nimParams2 <- c(
      "z", "sxy")
    
    # save(nimData,
    #      nimConstants,
    #      nimParams,
    #      nimParams2,
    #      modelCode,
    #      nimInits,
    #      file = file.path(myVars$WD, myVars$modelName,
    #                       paste("Snap",myVars$modelName,years[t],"_", ch, ".RData", sep = "")))
    #####
    
  }
  
  
}

length(model$z)


model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  
model$calculate()




model$calculate("y.alive")
which(is.infinite(model$logProb_y.alive), arr.ind = TRUE)
which(is.infinite(model$logProb_z), arr.ind = TRUE)
which(is.infinite(model$logProb_z), arr.ind = TRUE)

model$calculate("sxy")
which(model$sxy == -Inf, arr.ind = TRUE)
model$calculate("z")

model$calculate("dispSigma")
model$calculate("p0")
model$calculate("y.dead")

cmodel <- compileNimble(model)

### ====    11.2. CHECK INITIAL LOG-LIKELIHOODS ====
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



### ====    11.3. RESUME MODEL CONFIGURATION ====

conf <- configureMCMC(model, monitors = nimParams, thin = 1)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model, mcmc = Rmcmc))
Cmodel <- compiledList$model
Cmcmc <- compiledList$mcmc



### ====    11.4. RUN NIMBLE MCMC IN SUCCESSIVE BITES ====

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
## ------ IV.PROCESS RESULTS ------

### ===== 1. LOAD & PROCESS NIMBLE OPSCR OUTPUTS ====

# List the directories containing bite outputs
myVars$modelName <- myVars$modelName#"23.J_Fa"
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep(paste("NimbleOutFOR",myVars$modelName,sep=""), list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

# Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x){
  files <- list.files(x)
  files <- files[grep(".RData", files)]
  length(files)
}))
minBites <- min(numBites)/2


#Niter to Remove (burn-in) #CM
NSkipBites <- 15
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
    params.simple <- sapply(strsplit(colnames(this.sample), "\\["), "[", 1)
    parmIndex <- which(! params.simple %in% c("sxy","z"))
    nthins <- seq(1,dim(this.sample)[1], by=nthin)
    out[[x]] <- this.sample[nthins,parmIndex]#[ ,parmIndex] 
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



# for(i in 1: length(nimOutput)){
#   
#   nimOutput[[i]][,"sigma"] <-  nimOutput[[i]][,"sigma"]*myHabitat$resolution
#   nimOutput[[i]][,"lambda"] <-  nimOutput[[i]][,"lambda"]*myHabitat$resolution
# }

{#doall[CM]
  
  ### ==== 2. PLOT PARAMETERS ESTIMATES ====
  
  # pdf(file = file.path( myVars$WD,
  # myVars$modelName,paste(myVars$modelName,".pdf",sep="")
  # ))
  
  ### ====    2.1.PLOT DETECTIONS ====
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
  
  # ## PLOT GLOBAL MAP
  # par(mar = c(0,0,0,0))
  # plot(GLOBALMAP, col = "gray80")
  # plot(myStudyArea, col = "red", add = TRUE)
  
  ### ====    2.1.N ==== 
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
  
  ### ====    2.2.h ==== 
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
  
  ### ====    2.3.gamma ==== 
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
  
  ### ====    2.4.w ==== 
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
  
  ### ====    2.5.phi ====  
  # ### ==== 8.CALCULATE TRANSITION RATES FROM Z ====
  load(file.path(myVars$WD, myVars$modelName,
                 paste(myVars$modelName,"Chain",1, ".RData", sep = "")))
  
  phi <- phiind1 <- culled<- recruit  <- 0
  
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
  plot(-10, xlim = c(0,nYears), ylim=c(0,1), type ="n", xaxt="n", xlab = "Years", ylab = "Realized phi from z")
  axis(1, at = 1:(nYears-1) , labels = years[1:(nYears-1)])
  yr <- c(1:(nYears-1))
  points(phi~yr,pch=16)
  
  
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
  
  ### ====    2.6.p0 ====  
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
      # }#g
    }#t
  }#c
  # legend(x = 0, y = 0.25, legend = c("Norway", "Sweden"), fill= myCol)
  
  par(mfrow = c(1,2))
  myDev <- c(-0.3, 0.3)
  myCol <- c("blue4", "yellow3")
  
  # COUNTIES_AGGREGATEDSubsetsimp <- gSimplify(COUNTIES_AGGREGATEDSubset,tol = 500,topologyPreserve = T)
  # COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique
  # 
  for(c in 1:(dim(myResults$sims.list$p0Oth)[2])){
    plot(st_geometry(myStudyArea))
    if(c %in% 3){ 
      plot(st_geometry(myDetectors$main.detector.sp[detCounties %in% 1,]), add=T, col="red")
      
    }else{
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
      # }#g
    }#t
  }#c
  
  
  params <- dimnames(nimOutput[[1]])[[2]][grep("p0", dimnames(nimOutput[[1]])[[2]])[-1]]
  for(i in 1:length(params)){
    PlotJagsParams(jags.samples = nimOutput, params = params[i])
  }
  
  ### ====    2.7.the rest ====
  #structured<
  
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



### ===== 2. LOAD & PROCESS SNAP OUTPUTS ====

# List the directories containing bite outputs
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut",
                                                                          list.files(file.path(myVars$WD,
                                                                                               myVars$modelName)))]

outDirectories <- outDirectories[grep("Snap",outDirectories)]

#path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

t=1
nimOutputList <- myResultsList <-  list()

for(t in 1:nYears){
  path.list <- 0
  for(ch in 1:4){
    path.list[ch] <- file.path(myVars$WD, myVars$modelName, paste("NimbleOutFORSnap" ,
                                                                  myVars$modelName,years[t],"_", ch, ".RData", sep = "")) 
  }
  
  # Retrieve the minimum number of bites per chain
  numBites <- unlist(lapply(path.list, function(x){
    files <- list.files(x)
    files <- files[grep(".RData", files)]
    length(files)/2
  }))
  
  path.list <- path.list[!numBites%in%0]
  numBites <- numBites[!numBites%in%0]
  minBites <- min(numBites)
  #Niter to Remove (burn-in) #CM
  NSkipBites <- 50
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
  
  
  lapply(RUNTIME, function(x) x/3600)#[CM]
  unlist(lapply(RUNTIME, function(x) x/3600))#[CM]
  # TIME <- lapply(RUNTIME, function(x) x/3600)
  # 
  # max <- unlist(lapply(RUNTIME, function(x) x/3600))
  # at=c(1:length(TIME[[1]]))
  # plot(TIME[[1]]~at, pch=16, col=adjustcolor("red", alpha.f = 0.5), xlim=c(0,length(TIME[[1]])+2),
  #      ylim=c(0,10), ylab="time hours", xlab="bite number")
  # for(i in 2:length(TIME)){
  #    points(TIME[[i]] ~ at, )
  # }
  
  
  
  for(i in 1: length(nimOutput)){
    
    #nimOutput[[i]][,"sigma[2]"] <-  nimOutput[[i]][,"sigma[2]"]*myHabitat.list$resolution
    nimOutput[[i]][,"sigma"] <-  nimOutput[[i]][,"sigma"]*myHabitat$resolution
    
  }
  
  nimOutputList[[t]] <- nimOutput <- as.mcmc.list(nimOutput)
  
  #nimOutput <- as.mcmc.list(nimOutput)
  myResultsList[[t]] <- myResults <- ProcessCodaOutput(nimOutput, params.omit = c("sxy","z"))
  
} 
gc()

### combine all years together 
myResultsListALL <-  myResultsList[[1]] 

myResultsListALL$sims.list$N <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$N )) 
myResultsListALL$sims.list$betaDens <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$betaDens )) 
myResultsListALL$sims.list$betaResponse <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$betaResponse )) 
myResultsListALL$sims.list$betaResponseOth <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$betaResponseOth )) 

myResultsListALL$sims.list$pResponse <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$pResponse )) 
myResultsListALL$sims.list$sigma <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$sigma )) 
myResultsListALL$sims.list$psi <- do.call(cbind,lapply(myResultsList, function(x) x$sims.list$psi )) 

myResultsListALL$sims.list$betaCovs <- array(NA,c(dim(myResultsListALL$sims.list$betaCovs),nYears) )
myResultsListALL$sims.list$betaCovsOth <- array(NA,c(dim(myResultsListALL$sims.list$betaCovsOth),nYears) )

#myResultsListALL$sims.list$sxy <- array(0,c(dim(myResultsListALL$sims.list$sxy),nYears) )
# myResultsListALL$sims.list$z <- array(0,c(dim(myResultsListALL$sims.list$z),nYears) )
myResultsListALL$sims.list$p0 <- array(NA,c(dim(myResultsListALL$sims.list$p0),nYears) )
myResultsListALL$sims.list$p0Oth <- array(NA,c(dim(myResultsListALL$sims.list$p0Oth),nYears) )

#myResultsListALL$sims.list$sigma <- array(NA,c(dim(myResultsListALL$sims.list$sigma),nYears) )

for(t in 1:nYears){
  myResultsListALL$sims.list$betaCovs[,,t] <- myResultsList[[t]]$sims.list$betaCovs
  myResultsListALL$sims.list$betaCovsOth[,,t] <- myResultsList[[t]]$sims.list$betaCovsOth
  
  # myResultsListALL$sims.list$sigma[,,t] <- myResultsList[[t]]$sims.list$sigma
  
  myResultsListALL$sims.list$p0[,,t] <- myResultsList[[t]]$sims.list$p0
  myResultsListALL$sims.list$p0Oth[,,t] <- myResultsList[[t]]$sims.list$p0Oth
  
  #myResultsListALL$sims.list$sxy[,1:dim(myResultsList[[t]]$sims.list$sxy)[2],,t] <- myResultsList[[t]]$sims.list$sxy
  #myResultsListALL$sims.list$z[,1:dim(myResultsList[[t]]$sims.list$sxy)[2],t] <- myResultsList[[t]]$sims.list$z
}



{#doall[CM]
  
  ### ==== 2. PLOT PARAMETERS ESTIMATES ====
  ##pdf(file = file.path( myVars$WD,
  myVars$modelName,paste(myVars$modelName,"Snap.pdf",sep="")
  ))


### ====    2.1.N ==== 
par(mfrow = c(1,1), mar = c(5,5,5,5))
plot(10, xlim = c(0, nYears+1), ylim = c(200,1200), type ="n", xaxt="n", xlab = "Years", ylab = "N")
axis(1, c(1:nYears),labels = years)
for(t in 1:nYears){
  plot.violins(list(myResultsListALL$sims.list$N[,t]),
               x = t,
               at = t,
               violin.width = 0.3,
               col = "firebrick3",
               add = T,
               alpha = 0.3,
               border.col = "firebrick3")
}#t


### ====    2.6.p0 ====  
## by country and trap-response
par(mfrow = c(1,2))
myDev <- c(-0.3, 0.3)
myCol <- c("blue4", "yellow3")

COUNTIES_AGGREGATEDSubsetsimp <- st_simplify(COUNTIES_AGGREGATEDSubset,dTolerance = 500,preserveTopology = T)
COUNTIES_AGGREGATEDSubsetsimp$idunique <- COUNTIES_AGGREGATEDSubset$idunique

for(c in 1:dim(myResultsListALL$sims.list$p0)[2]){
  plot(st_geometry(myStudyArea))
  plot(st_geometry(COUNTIES_AGGREGATEDSubset[COUNTIES_AGGREGATEDSubsetsimp$idunique %in% c, ]), add=T, col="red")
  
  
  plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
  axis(1, at = 1:(nYears), labels = years[1:(nYears)])
  
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$p0[ ,c,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = "red",
                 add = T,
                 alpha = 0.3,
                 border.col = "red")
    # }#g
  }#t
}#c
# legend(x = 0, y = 0.25, legend = c("Norway", "Sweden"), fill= myCol)

par(mfrow = c(1,2))
myDev <- c(-0.3, 0.3)
myCol <- c("blue4", "yellow3")


for(c in 1:(dim(myResultsListALL$sims.list$p0Oth)[2])){
  plot(st_geometry(myStudyArea))
  if(c %in% 3){ 
    plot(st_geometry(myDetectors$main.detector.sp[detCounties %in% 1,]), add=T, col="red")
    
  }else{
    plot(st_geometry(COUNTRIES[c, ]), add=T, col="red")
  }
  
  plot(-10, xlim = c(0,nYears+1), ylim=c(0,0.06), type ="n", xaxt="n", xlab = "Years", ylab = "p0")
  axis(1, at = 1:(nYears), labels = years[1:(nYears)])
  
  for(t in 1:nYears){
    plot.violins(list(myResultsListALL$sims.list$p0Oth[ ,c,t]),
                 x = t ,
                 at = t ,
                 violin.width = 0.2,
                 col = "red",
                 add = T,
                 alpha = 0.3,
                 border.col = "red")
    # }#g
  }#t
}#c



### ====    2.7.the rest ====
#structured
par(mfrow = c(1,1))
trapCovs <- c("Tracks","Snow")
for(b in 1:dim(myResultsListALL$sims.list$betaCovs)[2]){
  plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n", xlab = "beta Structured", ylab = "beta",main=trapCovs[b])
  abline(h=0)
  axis(1, at = 1:nYears , labels = years)
  myCol <- c("firebrick3")
  for(t in 1:dim(myResultsListALL$sims.list$betaCovs)[3]){
    plot.violins(list(myResultsListALL$sims.list$betaCovs[ ,b,t]),
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
for(b in 1:dim(myResultsListALL$sims.list$betaCovsOth)[2]){
  plot(-10, xlim = c(0,nYears+0.5), ylim=c(-2,2), type ="n", xaxt="n", xlab = "beta Other", ylab = "beta",main=trapCovs[b])
  abline(h=0)
  axis(1, at = 1:nYears , labels = years)
  myCol <- c("firebrick3")
  for(t in 1:dim(myResultsListALL$sims.list$betaCovsOth)[3]){
    plot.violins(list(myResultsListALL$sims.list$betaCovsOth[ ,b,t]),
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
for(t in 1:dim(myResultsListALL$sims.list$sigma)[2]){
  plot.violins(list(myResultsListALL$sims.list$sigma[ ,t]),
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
for(t in 1:dim(myResultsListALL$sims.list$betaResponse)[2]){
  plot.violins(list(myResultsListALL$sims.list$betaResponse[ ,t]),
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
for(t in 1:dim(myResultsListALL$sims.list$betaResponseOth)[2]){
  plot.violins(list(myResultsListALL$sims.list$betaResponseOth[ ,t]),
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
     xlab = "beta Other", ylab = "beta Dens")
abline(h=0)
axis(1, at = 1:nYears , labels = years)
myCol <- c("firebrick3")
for(t in 1:dim(myResultsListALL$sims.list$betaDens)[2]){
  plot.violins(list(myResultsListALL$sims.list$betaDens[ ,t]),
               x = t,
               at = t,
               violin.width = 0.2,
               col = myCol,
               add = T,
               alpha = 0.3,
               border.col = myCol)
}#t




dev.off()
}#doall[CM]





## TRY DENSITY ESTIMATES 
gc()
# ORIGIN DENSITY
habbR <- myHabitat$habitat.r#raster::aggregate(myHabitat.list$habitat.r,fact=2)
habbR <- raster::disaggregate(myHabitat$habitat.r, fact=1)
habbR[habbR>0] <- 1
habbdens <- habbR
plot(habbdens)
dim(myResults$sims.list$sxy)
gc()

# RESCALE THE COORDINATES 
sxyScaled <- GridToUTM(data.sxy = myResults$sims.list$sxy,
                       grid.sp = myHabitat$habitat.sp )
gc()

# CHE
plot(habbR)
dim(sxyScaled$data.scaled.xy)
points(sxyScaled$data.scaled.xy[1500,,2,1]~sxyScaled$data.scaled.xy[1500,,1,1], pch=16)
ite <- (dim(sxyScaled$data.scaled.xy)[1]-10000):dim(sxyScaled$data.scaled.xy)[1]
dim(myResults$sims.list$z)


Densi <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
  Densi[[t]] <- EstimateN_v3(habRaster = habbR, # RASTER  FOR WHICH REGION ESTIMATES SHOULD BE ESTIMATED. 
                             posterior.sxy = sxyScaled$data.scaled.xy[,,,t],#sxyScaled$data.scaled.xy[ite,,,t],
                             posterior.z =  myResults$sims.list$z[,,t],#myResults$sims.list$z[ite,,t],
                             alive.states =  c(2),
                             return.all = FALSE,
                             regionEstimates = FALSE)
  print(t)
}
Densi[[t]]$PosteriorsCellsMean

max <- max(unlist(lapply(Densi, function(x) max(x$PosteriorsCellsMean))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

##pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_DENSITY.pdf",sep="")))

for(t in 1: dim(sxyScaled$data.scaled.xy)[4]){
  habbdens[] <- Densi[[t]]$PosteriorsCellsMean
  
  plot(habbdens, breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
  plot(myStudyArea.poly,add=T, border=grey(0.5))
  points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t],],
         pch=16, cex=0.4, col=adjustcolor("black",alpha.f = 0.2))
  plot(habbdens, legend.only=TRUE,breaks=cuts, col=col,
       legend.width = 2,
       axis.args=list(at=round(seq(0, max, length.out = 5),digits = 1),
                      labels=round(seq(0, max, length.out = 5),digits = 1), 
                      cex.axis=0.6),
       legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
  
}
dev.off()


### ===== COUNTRY-SPECIFIC DENSITY ESTIMATES ====

habbR <- myHabitat$habitat.r#raster::aggregate(myHabitat.list$habitat.r,fact=2)
plot(myHabitat$habitat.r)


plot(myHabitat$habitat.poly,add=T)
habbR <- raster::disaggregate(myHabitat$habitat.r, fact=1)
# identify SWEDEN in the raster
SWE <- COUNTRIES[which(COUNTRIES$ISO %in% c("SWE")),]     ## Just take Sweden
NOR <- COUNTRIES[which(COUNTRIES$ISO %in% c("NOR")),]     ## Just take Norway
NOR <- gSimplify(spgeom = NOR, tol = 500)
SWE <- gSimplify(spgeom = SWE, tol = 500)
this.r <- RasterizePolygon(poly=SWE, r=habbR,CoverToKeepHabitat=50,fasterize=TRUE)
habbR[this.r==1] <- 2
this.r <- RasterizePolygon(poly=NOR, r=habbR,CoverToKeepHabitat=50,fasterize=TRUE)
habbR[this.r==1] <- 1
# 
# habbR.sp <- SpatialPoints(coordinates(habbR),proj4string = CRS(proj4string(myHabitat.list$habitat.r)))
# spIN <- sp::over(habbR.sp, SWE)
# habbR[!is.na(spIN[])] <- 2
plot(habbR)
##remove buffer from habitat
e <- extent(myStudyArea.poly)
e.sp <- as(e, 'SpatialPolygons')  
habitat.r <- myHabitat$habitat.r
habitat.r[habitat.r!=1] <- NA
habbR <- mask(habbR, e.sp)
# habbR.sp <- SpatialPoints(coordinates(habbR),proj4string = CRS(proj4string(myHabitat.list$habitat.poly)))
# spOUTBUFF <-sp::over(habbR.sp, habitat.r.wobuffer)
# spOUTBUFF <- sp::over(habbR.sp, myHabitat.list$habitat.poly)
# habbR[is.na(spOUTBUFF[,1])] <- 0
# plot(habbR)

# convert to text 
habbR[habbR==2] <- "SWE"
habbR[habbR==1] <- "NOR"
habbR[habbR==0] <- NA
gc()
## APPLY THE FUNCTION 
DensiCountries <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
  DensiCountries[[t]] <- EstimateN_v3(habRaster = habbR, # RASTER  FOR WHICH REGION ESTIMATES SHOULD BE ESTIMATED. 
                                      posterior.sxy = sxyScaled$data.scaled.xy[,,,t],#sxyScaled$data.scaled.xy[ite,,,t],
                                      posterior.z =  myResults$sims.list$z[,,t],#myResults$sims.list$z[ite,,t],
                                      alive.states =  c(2),
                                      return.all = FALSE,
                                      regionEstimates = TRUE)
  print(t)
}
DensiCountries[[1]]$PosteriorsRegion$`1`
DensiCountries[[1]]$summary

### ===== COUNTY SPECIFIC DENSITY ESTIMATES ====

habbRCounties <- myHabitat$habitat.r#raster::aggregate(myHabitat.list$habitat.r,fact=2)
plot(myHabitat$habitat.r)

## get the counties 
COUNTIES <- aggregate(x = COMMUNES, by = "NAME_1")
COUNTIESsimp <- gSimplify(COUNTIES,tol = 500, topologyPreserve = T)
COUNTIESsimp$NAME_1 <- COUNTIES$NAME_1
plot(myHabitat$habitat.poly,add=T)
plot(COUNTIESsimp, add=T)
# habbRCounties.sp <- SpatialPoints(coordinates(habbRCounties),proj4string = CRS(proj4string(myHabitat.list$habitat.r)))

##remove buffer from habitat
habbRCountieswbuff <- habbRCounties <- mask(habbRCounties, e.sp)

for(i in 1:length(COUNTIESsimp)){
  # identify SWEDEN in the raster
  this.r <- RasterizePolygon(poly = COUNTIESsimp[i,], r = habbRCountieswbuff, CoverToKeepHabitat = 50, fasterize = TRUE)
  habbRCounties[this.r==1] <- as.character(COUNTIESsimp$NAME_1[i])
}
plot(habbRCounties)

# habbRCounties.sp <- SpatialPoints(coordinates(habbRCounties),proj4string = CRS(proj4string(myHabitat.list$habitat.poly)))
# habbRCounties.sp$id <- rep(1,length(habbRCounties.sp))
# spOUTBUFF <- sp::over(habbRCounties.sp, habitat.r.wobuffer)
# habbRCounties[spOUTBUFF[,1]==0] <- NA
# habbRCounties[habbRCounties==0] <- NA
# habbRCounties[habbRCounties==1] <- NA
# habbRCounties[] <- as.numeric(as.factor(habbRCounties[]))+4

plot(habbRCounties)
plot(myHabitat$habitat.poly,add=T)



gc()

DensiCounties <- list()
for(t in 1:dim(sxyScaled$data.scaled.xy)[4]){
  DensiCounties[[t]] <- EstimateN_v3(habRaster = habbRCounties
                                     , # RASTER  FOR WHICH REGION ESTIMATES SHOULD BE ESTIMATED. 
                                     posterior.sxy = sxyScaled$data.scaled.xy[,,,t]
                                     ,#sxyScaled$data.scaled.xy[ite,,,t],
                                     posterior.z =  myResults$sims.list$z[,,t]
                                     ,#myResults$sims.list$z[ite,,t],
                                     alive.states = c(2)
                                     ,
                                     return.all = FALSE
                                     ,
                                     regionEstimates = TRUE)
  print(t)
}



### ===== PLOT COUNTIES/COUNTRIES SPECIFIC DENSITY ESTIMATES ====
## TOTAL 
{#doall
  ##pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_NRegions.pdf",sep="")))
  plot(-1000, xlim=c(0.5,nYears+0.5), ylim=c(0,1200), ylab="N Total",xaxt="n")
  axis(1, at=c(1:nYears), labels = years)
  at = c(1:nYears)
  # points(unlist(lapply(DensiCountries, function(x) x$summary["1","mean"]))~at)
  # points(unlist(lapply(DensiCountries, function(x) x$summary["2","mean"]))~at)
  country.colors <- c("goldenrod1","goldenrod3")
  names(country.colors) <- c("Norway","Sweden")
  for(t in 1:nYears){
    plot.violins2(DensiCountries[[t]]$PosteriorsRegion["ALL"],
                  x = at[t],
                  at=at[t],
                  violin.width = 0.3,
                  col = "red",
                  alpha = 0.2,
                  border.col = "red",
                  add = T)
    text(DensiCountries[[t]]$summary["Total","median"],x=t,y=DensiCountries[[t]]$summary["Total","median"]+50)
    
  }
  
  
  ## COUNTRIES 
  plot(-1000, xlim=c(0.5,nYears+0.5), ylim=c(0,700), ylab="N",xaxt="n")
  axis(1, at=c(1:nYears), labels = years)
  at = c(1:nYears)
  # points(unlist(lapply(DensiCountries, function(x) x$summary["1","mean"]))~at)
  # points(unlist(lapply(DensiCountries, function(x) x$summary["2","mean"]))~at)
  country.colors <- c("goldenrod1","goldenrod3")
  names(country.colors) <- c("Norway","Sweden")
  for(t in 1:nYears){
    plot.violins2(DensiCountries[[t]]$PosteriorsRegion["NOR"],
                  x = at[t],
                  at=at[t],
                  violin.width = 0.3,
                  col = country.colors[1],
                  alpha = 0.2,
                  border.col = country.colors[1],
                  add = T)
    text(DensiCountries[[t]]$summary["NOR","median"],x=t,y=DensiCountries[[t]]$summary["NOR","median"]+25)
    plot.violins2(DensiCountries[[t]]$PosteriorsRegion["SWE"],
                  x = at[t],
                  at=at[t],
                  violin.width = 0.3,
                  col = country.colors[2],
                  alpha = 0.2,
                  border.col =  country.colors[2],
                  add = T)
    text(DensiCountries[[t]]$summary["SWE","median"],x=t,y=DensiCountries[[t]]$summary["SWE","median"]+25)
    
  }
  
  legend("topright", legend = names(country.colors), fill=country.colors, bty="n" )
  
  ## COUNTY
  # points(unlist(lapply(DensiCountries, function(x) x$summary["1","mean"]))~at)
  # points(unlist(lapply(DensiCountries, function(x) x$summary["2","mean"]))~at)
  country.colors <- c("goldenrod1","goldenrod3")
  names(country.colors) <- c("Norway","Sweden")
  idcounty <- rownames(DensiCounties[[1]]$summary)
  par(mfrow=c(2,2))
  
  for(co in 2:(length(idcounty)-1)){
    plot(-1000, xlim=c(0.5,nYears+0.5), ylim=c(0,300), ylab=paste("N",idcounty[co],sep=""),xaxt="n")
    axis(1, at=c(1:nYears), labels = years)
    at = c(1:nYears)
    
    for(t in 1:nYears){
      
      plot.violins2(DensiCounties[[t]]$PosteriorsRegion[idcounty[co]],
                    x = at[t],
                    at=at[t],
                    violin.width = 0.3,
                    col = country.colors[1],
                    alpha = 0.2,
                    border.col = country.colors[1],
                    add = T)
      text(DensiCounties[[t]]$summary[idcounty[co],"median"],x=t,y=DensiCounties[[t]]$summary[idcounty[co],"median"]+15)
      
    }      
    
    plot(myStudyArea.poly)
    plot(COUNTIESsimp,add=T)
    plot(COUNTIESsimp[COUNTIESsimp$NAME_1==idcounty[co],], col="red",add=T)
    
  }
  
  dev.off()      
  
}#doall





### ===== 3. PLOT DENSITY MAPS ====
### =====    3.1. PROCESS SXY OUTPUT ====
## List the directories containing bite outputs
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut", list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)

## Retrieve the minimum number of bites per chain
numBites <- unlist(lapply(path.list, function(x)length(list.files(x))))
minBites <- min(numBites)

nimOutput <- list()
for(p in 1:length(path.list)){
  print(path.list[p])
  outfiles <- list.files(path.list[p])
  out <- list()
  for(x in 61:minBites){
    print(x)
    load(file.path(path.list[p], paste("bite_", x, ".RData", sep = "")))
    params.simple <- sapply(strsplit(colnames(mcmcSamples), "\\["), "[", 1)
    parmIndex <- which(params.simple %in% c("sxy","z"))
    out[[x]] <- mcmcSamples[ ,parmIndex] 
  }#x
  out.mx <- do.call(rbind, out)
  nimOutput[[p]] <- as.mcmc(out.mx)
}#p

nimOutput <- as.mcmc.list(nimOutput)
myResults <- ProcessCodaOutput(x = nimOutput)

## RESCALE POSTERIOR sxy
postSxy <- GridToUTM( data.sxy = myResults$sims.list$sxy,
                      grid.sp = myHabitat$habitat.sp)$data.scaled.xy

### =====    3.2. PLOT DENSITY MAPS ====
## NEW GetDensity FUNCTION
GetPopDensityFromPolygon <- function(habPolygon = NULL,
                                     numDiv.x = 1,
                                     numDiv.y = 1,
                                     habResolution = NULL,
                                     habRaster = NULL,
                                     posterior.sxy,
                                     posterior.z,
                                     alive.states,
                                     plot.check = FALSE){
  
  ## Filter out posterior sxy for dead individuals
  deadId <- !which(posterior.z %in% alive.states, arr.ind = TRUE)
  post.x <- posterior.sxy[ , ,1]
  post.x[deadId] <- -999
  posterior.sxy[ , ,1] <- post.x
  
  ## Divide the original polygon using SF package
  habPolygons <- st_make_grid(x = habPolygon, n = c(numDiv.x, numDiv.y))
  habPolygons <- as(habPolygons, "Spatial")
  habPolygons <- habPolygons[!is.na(habPolygons %over% habPolygon)]
  habPolygons <- intersect(habPolygons, as(habPolygon, "SpatialPolygons"))
  
  ## Initialize list of outputs
  thisMean <- thisLowerCI  <- thisUpperCI  <- thisMedian <- list()
  
  for(n in 1:length(habPolygons)){
    
    print(paste("Processing raster bite #", n, " of ", length(habPolygons), sep = ""))
    
    ## Create a raster of the zone of interest
    habRaster <- raster( habPolygons[n],
                         resolution = habResolution,
                         origin = c(extent(habPolygon)[c(1,3)]))
    # res(habRaster) <- habResolution
    habRaster <- rasterize(habPolygons[n], habRaster)
    # origin(habRaster) <- c(150000, 6400000) # c(extent(habPolygon)[c(1,3)])
    
    # Initialize the list of output
    thisMean[[n]] <- thisLowerCI[[n]] <- thisUpperCI[[n]] <- thisMedian[[n]] <- habRaster
    
    if(all(is.na(habRaster[]))){next}
    
    ## Convert habRaster to SpatialGrid object 
    spRaster <- as(habRaster, "SpatialGridDataFrame")
    
    ## Calculate cell-specific density for each iteration
    Density <- apply(X = posterior.sxy, FUN = function(curCoords, inRaster){
      ## Get the cell index for each coordinate
      curCoords <- curCoords[curCoords[ ,1] != -999, ]
      curGridIndeces <- getGridIndex( curCoords, getGridTopology(inRaster),
                                      all.inside = FALSE)
      ## Get rid of coordinate indeces that do not fall in a cell
      curGridIndeces <- curGridIndeces[!is.na(curGridIndeces)]
      ## Create a frequency table of grid indeces
      freqTable <- table(as.character(curGridIndeces))
      ## Initialise an output vector
      outVector <- rep(0, nrow(coordinates(inRaster)))
      outVector[as.integer(names(freqTable))] <- freqTable
      outVector[is.na(spRaster@data$layer)] <- NA
      outVector
    }, inRaster = spRaster, MARGIN = 1)
    
    ## EXTRACT SUMMARY STATISTICS PER HABITAT CELL
    meanDensity <- apply(Density, 1, mean)
    CI <- apply(Density, 1, function(x)quantile(x, c(0.025,0.5,0.975), na.rm = TRUE))
    
    ## FEED IN RASTER
    thisMean[[n]][] <- meanDensity
    thisLowerCI[[n]][] <- CI[1, ]
    thisMedian[[n]][] <- CI[2, ]
    thisUpperCI[[n]][] <- CI[3, ]
  }#n
  
  if(plot.check){
    par(mfrow = c(1,2))
    ## Plot mean density 
    max.cut <- max(meanDensity.r[], na.rm = TRUE)
    cuts <- seq(0, max.cut, length.out = 101) 
    pal <- rev(terrain.colors(length(cuts)))
    plot(meanDensity.r, breaks = cuts, col = pal, main ="Mean density", legend=FALSE, axes=FALSE)
    plot(meanDensity.r, legend.only=TRUE, col=pal,
         legend.width=1, legend.shrink=1,
         axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                        labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                        cex.axis=0.6),
         legend.args = list(text=paste("Density (ind.",(habResolution/1000)^2,"km-2)",sep=""), side=4, font=2, line=2.5, cex=0.8))
    plot(habPolygon, add = TRUE)
    
    ## Plot CI width
    max.cut <- max(upCIDensity.r[]-lowCIDensity.r[], na.rm = TRUE)
    cuts <- seq(0, max.cut, length.out = 101) 
    pal <- rev(heat.colors(length(cuts)))
    plot(upCIDensity.r-lowCIDensity.r, breaks = cuts, col = pal, main ="95% CI width", legend = FALSE, axes=FALSE)
    plot(upCIDensity.r, legend.only=TRUE, col=pal,
         legend.width=1, legend.shrink=1,
         axis.args=list(at=round(seq(0, max.cut, length.out = 11),digits = 2),
                        labels=round(seq(0, max.cut, length.out = 11),digits = 2), 
                        cex.axis=0.6))
    plot(habPolygon, add = TRUE)
  }
  return(list("mean.Density" = test <- do.call(merge, thisMean),
              "median.Density" = do.call(merge,thisMedian),
              "upperCI.Density" = do.call(merge,thisUpperCI),
              "lowerCI.Density" = do.call(merge,thisLowerCI)))
}
## CALCULATE DENSITY
# DENSITY <- list()
# for(t in 1:dim(postSxy)[4]){
#    DENSITY[[t]] <- GetPopDensityFromPolygon(habPolygon = myStudyArea,
#                                             numDiv.x = 10,
#                                             numDiv.y = 20,
#                                             habResolution = 20000,
#                                             posterior.sxy = postSxy[ , , ,t],
#                                             posterior.z = myResults$sims.list$z[,,t],
#                                             alive.states = 2,
#                                             plot.check = F)
#    plot(DENSITY[[t]]$mean.Density)
#    plot(myStudyArea, add = TRUE)
# }#t
# 
# save(DENSITY, file = file.path( myVars$WD,
#                                 myVars$modelName,
#                                 paste("9JM_DENSITY_PD.RData", sep = "")))
## PLOT DENSITY
#pdf(file = file.path( myVars$WD,
# myVars$modelName,
# paste("9JM_DENSITY_PD.pdf", sep = "")))

max.cut <- max(unlist(lapply(DENSITY, function(x) x$mean.Density[])), na.rm = TRUE)
cuts <- seq(0, max.cut, length.out = 11) 
pal <- rev(terrain.colors(length(cuts)))
# max.cut1 <- max(unlist(lapply(myPopDensity, function(x) x$upperCI.Density[]-x$lowerCI.Density[])), na.rm = TRUE)
# cuts1 <- seq(0, max.cut1, length.out = 101) 
# pal1 <- rev(heat.colors(length(cuts1)))

for(t in 1:nYears){
  par(mfrow = c(1,1))
  
  ## Plot mean density 
  plot(DENSITY[[t]]$mean.Density, breaks = cuts, col = pal,
       main = paste( "N =", round(sum(DENSITY[[t]]$mean.Density[], na.rm = T))),
       legend = F, axes = F)  
  plot(myStudyArea, add = T, border = grey(0.5))
  plot(myHabitat$buffered.habitat.poly, border = grey(0.5), add = TRUE)
  plot(DENSITY[[t]]$mean.Density, legend.only = T, col = pal,
       legend.width = 1, legend.shrink = 1,
       axis.args = list( at = round(seq(0, max.cut, length.out = 11), digits = 2),
                         labels = round(seq(0, max.cut, length.out = 11), digits = 2), 
                         cex.axis = 0.6),
       legend.args = list( text = paste("Density (ind.",(myVars$OUTPUT$mapResolution/1000)^2,"km-2)",sep=""),
                           side = 4, font = 2, line = 2.5, cex = 0.8))
  # ## Plot CI width
  # plot(myPopDensity[[t]]$upperCI.Density-myPopDensity[[t]]$lowerCI.Density, breaks = cuts1, col = pal1, main ="95% CI width", legend = FALSE, axes=FALSE)
  # plot(myStudyArea, add = T, border = grey(0.5))
  # plot(myHabitat$buffered.habitat.poly, border = grey(0.5), add = TRUE)
  # plot(myPopDensity[[t]]$upperCI.Density-myPopDensity[[t]]$lowerCI.Density, legend.only=TRUE, col=pal1,
  #      legend.width=1, legend.shrink=1,
  #      axis.args=list(at=round(seq(0, max.cut1, length.out = 11),digits = 2),
  #                     labels=round(seq(0, max.cut1, length.out = 11),digits = 2), 
  #                     cex.axis=0.6))
  # plot(myHabitat$buffered.habitat.poly, add = TRUE)
}#t
dev.off()

### =====    3.3. PLOT DENSITY MAPS (CM) ====
## TRY EstimateN_v3 (CM)
gc()
habbR <- raster::disaggregate(myHabitat$habitat.r, fact = 2)
habbR[habbR > 0] <- 1
plot(habbR)

## CALCULATE DENSITY RASTERS
Densi <- list()
for(t in 1:dim(postSxy)[4]){
  Densi[[t]] <- EstimateN_v3( habRaster = habbR, 
                              posterior.sxy = postSxy[,,,t],
                              posterior.z = myResults$sims.list$z[,,t],
                              alive.states = 2,
                              return.all = FALSE,
                              regionEstimates = FALSE)
  print(t)
}#t

## PLOT
max.cuts <- max(unlist(lapply(Densi, function(x) max(x$PosteriorsCellsMean))))
cuts <- round(seq(0, max.cuts, length.out = 101) , digits = 2)
pal <- rev(terrain.colors(length(cuts)))


for(t in 1:dim(postSxy)[4]){
  habdens[] <- Densi[[t]]$PosteriorsCellsMean
  
  plot(habdens, breaks = cuts, col = pal,
       main = paste( "N =", round(sum(Densi[[t]]$PosteriorsCellsMean[], na.rm = T))),
       legend = F, axes = F)   
  plot(myStudyArea, add = T, border = grey(0.5))
  plot(habdens, legend.only=TRUE, breaks=cuts, col=pal,
       legend.width = 2,
       axis.args = list(at=round(seq(0, max.cuts, length.out = 6), digits = 1),
                        labels=round(seq(0, max.cuts, length.out = 6), digits = 1),
                        cex.axis=0.6),
       legend.args=list(text='Density', side=4, font=2, line=2.5, cex=0.8))
}
dev.off()
##------------------------------------------------------------------------------
