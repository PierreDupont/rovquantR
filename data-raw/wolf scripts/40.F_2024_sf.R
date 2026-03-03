
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


## ------ SET REQUIRED WORKING DIRECTORIES ------
source("C:/My_documents/RovQuant/workingDirectories.R")             


## ------ SOURCE THE REQUIRED FUNCTIONS ------
source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/Nimble/dbin_LESSCachedAllSparseWolf.R")
source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/Nimble/dbinomLocal_normalWolf.R")


## -----------------------------------------------------------------------------

## ------   0. SET ANALYSIS CHARACTERISTICS -----

## ------   1. GENERAL VARIABLES DECLARATION ------ 

myVars <- list( 
  ## WORKING DIRECTORY & MODEL NAME
  WD = "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant/wolf/WolfRuns20152024",
  modelName = "40.F_2024_sf",
  
  # HABITAT SPECIFICATIONS
  HABITAT = list( x.extent = c(210000,740000),
                  y.extent = c(6000000,7050000),
                  habBuffer = 40000),
  
  # NGS DATA SPECIFICATIONS
  DATA = list( years = 2015:2024,   
               species = c("Ulv"),                ## "Ulv","Jerv","Bjorn"
               sex = c("Hunn"),                   ## "Hann","Hunn","Ukjent" 
               samplingMonths = list(10:12,1:4)), ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)
  
  # DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 1000,
                    detResolution = 10000,
                    detDeadResolution = 10000),
  
  # DATA GENERATION 
  DETECTIONS = list( maxDist = 45000,
                     aug.factor = 0.8),
  
  ## OUTPUT PLOTS 
  OUTPUT = list(mapResolution = 5000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

if(is.null(myVars$modelName))stop("YOU SHOULD PROBABLY CHOOSE A NAME FOR THIS ANALYSIS/MODEL")
if(is.null(myVars$WD))stop("YOU SHOULD PROBABLY CHOOSE A WORKING DIRECTORY FOR THIS ANALYSIS/MODEL")
if(!dir.exists(file.path(myVars$WD,myVars$modelName))){dir.create(file.path(myVars$WD, myVars$modelName))}



## -----------------------------------------------------------------------------

## ------ I. LOAD & SELECT DATA ------

## ------   1. HABITAT DATA ------ 

## ------     1.1. LOAD RAW SHAPEFILES ------ 

COUNTRIES <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/countries_multipart.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
COUNTRIES <- COUNTRIES[which(COUNTRIES$ISO %in% c("NOR","SWE","FIN")), ]                                 ## Just take Sweden and Norway
GLOBALMAP <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/Scandinavia_border_33NNoLakes.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
GLOBALMAP <- st_simplify(GLOBALMAP, dTolerance =  500)
COMMUNES_NOR <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp",sep=""))   ## Communal map of Norway
COMMUNES_SWE <- st_read(paste(dir.dropbox,"/DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp",sep=""))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)
COUNTIES <- COMMUNES %>% group_by(NAME_1) %>% summarize()

## LOAD POLYGONS OF WATER + HUMANS WITH AREAS >80000M2 (CREATED IN TEMP/CM/GIS/buildingsWaterPolygons.R)
COUNTRIESWaterHumans <- st_read(paste(dir.dropbox,"/DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
## SELECT POLGYONS WITH AREAS SIZE WITH WaterHumans >25000000 m2 (5*5km)
COUNTRIESWaterHumans <- COUNTRIESWaterHumans[COUNTRIESWaterHumans$ISO %in% c("SWE","NOR"), ]
## Remove small polygons (islands and things)
COUNTRIESWaterHumans <- COUNTRIESWaterHumans[COUNTRIESWaterHumans$area>80000000,]
plot(st_geometry(COUNTRIESWaterHumans))

## ------     1.2. CLEAN SHAPEFILES ------ 
# GLOBALMAP <- st_read(paste(dir.dropbox,"/DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp",sep="")) ## Map of Scandinavia (including Finland & parts of Russia)
# GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
# GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))
#plot(st_geometry(GLOBALMAP))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- COUNTRIESWaterHumans[COUNTRIESWaterHumans$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- COUNTRIES %>%    group_by(ISO) %>%summarize()

SWE <- COUNTRIES[which(COUNTRIES$ISO %in% c("SWE")),]     ## Just take Sweden
NOR <- COUNTRIES[which(COUNTRIES$ISO %in% c("NOR")),]     ## Just take Norway

NOR1 <- st_simplify(NOR, dTolerance = 500)
SWE1 <- st_simplify(SWE, dTolerance = 500)
# NOR_buffer <- st_buffer(NOR1, dist = 20000)
# NOR_buffer <- NOR_buffer %>%  group_by(ISO) %>%summarize()
# SWE_buffer <- st_buffer(SWE1, dist = 20000)
# SWE_buffer <- SWE_buffer %>%  group_by(ISO) %>%summarize()
# NOR_buffer <- st_difference(NOR_buffer,SWE)
# plot(st_geometry(NOR_buffer))
# 
# NOR_buffer <- st_intersection(NOR_buffer, SWE)
# 
# SWE_buffer <- st_difference(SWE_buffer, NOR)
NOR1$ID <- "S"
SWE1$ID <- "N"
country <- rbind(NOR1, SWE1)
plot(st_geometry(country))
## ------     1.3. SAVE SHAPEFILES OBJECTS FOR FASTER RUNS ------ 
save( GLOBALMAP,
      COUNTRIES,
      COUNTIES,
      COUNTRIESWaterHumans,
      COMMUNES,
      country, file = file.path(myVars$WD,"HABITATsf.RData"))
load(file.path(myVars$WD,"HABITATsf.RData"))

## ------     1.3. LOAD SCANDINAVIAN 20KM HABIAT  ------ 
load(paste(dir.dropbox,"/DATA/GISData/spatialDomain/Habitat20km.RData",sep=""))

## ------     1.4. CREATE STUDY AREA POLYGON ------ 
## CREATE STUDY AREA POLYGON BASED x AND y EXTENTS
myStudyArea.extent  <- st_bbox(extent(myVars$HABITAT$x.extent, myVars$HABITAT$y.extent))
st_crs(myStudyArea.extent) <- st_crs(COUNTRIESWaterHumans)
myStudyArea.poly <- st_crop(COUNTRIESWaterHumans, extent(myVars$HABITAT$x.extent, myVars$HABITAT$y.extent))
# to get only "polygons objects"
myStudyArea.poly <- st_collection_extract(myStudyArea.poly, "POLYGON")



## ------   2. NGS DATA ------ 
## ------     2.1.LOAD ROVBASE FILES ------ 
DNA <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/RIB22042025133456403_wolfDNA.csv",sep=""), fileEncoding="latin1")## NGS data from RovBase
DEAD <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/RIB22042025133534832_wolfDEAD.csv",sep=""), fileEncoding="latin1") ## Dead Recoveries from RovBase
INDIVIDUAL_ID <- read.csv(paste(dir.dropbox,"/DATA/RovbaseData/ROVBASE DOWNLOAD 20220513/220512_ID Grouping 2006-2021.csv",sep=""), fileEncoding="latin1")  ## Wolves infos from Micke
## THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2022/23.
Pack_ID2023 <- read.csv(paste(dir.dropbox,
                              "/DATA/MickeIndividual/Genetiskt ID RM vargar 2223 Bilaga 4_ØF.csv",sep=""),
                        fileEncoding="latin1")  
Pack_ID2024 <- read.csv(paste(dir.dropbox,
                              "/DATA/MickeIndividual/Bilaga_11.4_240424_ØF to Cyril.csv",sep=""),
                        fileEncoding="latin1")  

Pack_ID2025 <- read.csv(paste(dir.dropbox,
                              "/DATA/MickeIndividual/RovbaseID for Rovquant estimates2025FromOystein.csv",sep=""),
                        fileEncoding="latin1")  
## Here we need to recreate the sex columns as Oystein gave me a list of ids only (losing the sex)
Pack_ID2025$Sex <- apply(Pack_ID2025[,c("Sex1","Sex2","Sex3","Sex4")],1, function(x) x[which(!x%in% "")][1] )

## THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2023/24.

## ------     2.2.TRANSLATE SCANDINAVIAN CHARACTERS ------ 
colnames(DNA) <- translateForeignCharacters(dat=colnames(DNA), dir.translation = dir.analysis )
colnames(DEAD) <- translateForeignCharacters(dat=colnames(DEAD), dir.translation = dir.analysis )
colnames(INDIVIDUAL_ID) <- translateForeignCharacters(dat=colnames(INDIVIDUAL_ID), dir.translation = dir.analysis )
colnames(Pack_ID2023) <- translateForeignCharacters(dat=colnames(Pack_ID2023), dir.translation = dir.analysis )
colnames(Pack_ID2024) <- translateForeignCharacters(dat=colnames(Pack_ID2024), dir.translation = dir.analysis )
colnames(Pack_ID2025) <- translateForeignCharacters(dat=colnames(Pack_ID2025), dir.translation = dir.analysis )

## ------   3. SEARCH EFFORT DATA ------ 

## ------     3.1. GPS SEARCH TRACKS ------ 

# #LOAD GPS SEARCH TRACKS FROM ROVBASE
# TRACKS_SINGLE <- read_sf(paste(dir.dropbox, "/DATA/RovbaseData/TRACK DATA FROM BOUVET 20240426/Aktivitetslogg_20240426/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240426_dateSf.shp", sep = ""))
# TRACKS_MULTI <- read_sf(paste(dir.dropbox, "/DATA/RovbaseData/TRACK DATA FROM BOUVET 20240426/Aktivitetslogg_20240426/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240426_dateSf.shp", sep = ""))
# 
# TRACKS_SINGLE <- read_sf(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250422.shp", sep = ""))
# TRACKS_SINGLE$Dato <- as.POSIXct(strptime(TRACKS_SINGLE$Dato, "%Y-%m-%d"))
# TRACKS_SINGLE$Yr <- as.numeric(format(TRACKS_SINGLE$Dato,"%Y"))
# TRACKS_SINGLE$Mth <- as.numeric(format(TRACKS_SINGLE$Dato,"%m"))
# TRACKS_SINGLE$Dato <- as.character(TRACKS_SINGLE$Dato)
# 
# #
# TRACKS_MULTI <- read_sf(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422.shp", sep = ""))
# TRACKS_MULTI$Dato <- as.POSIXct(strptime(TRACKS_MULTI$Dato, "%Y-%m-%d"))
# TRACKS_MULTI$Yr <- as.numeric(format(TRACKS_MULTI$Dato,"%Y"))
# TRACKS_MULTI$Mth <- as.numeric(format(TRACKS_MULTI$Dato,"%m"))
# TRACKS_MULTI$Dato <- as.character(TRACKS_MULTI$Dato)
# 
# # ALL_TRACKS <- st_read(paste(dir.dropbox,
# #                               "//DATA//RovbaseData//ROVBASE DOWNLOAD 20230508//eksport_rovquant//Shapefile/spor.shp",sep=""))
# 
# # COMBINE ALL TRACKS
# ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI)
# # REMOVE HELICOPTER TRACKS
# ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helikopter=="0", ]
# 
# #TRACKS_MULTI <- readOGR(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20220513/AktivitetsloggRovQuant20220518/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20220512_date.shp", sep = ""))
# #check
# hist(ALL_TRACKS$Yr )
# 
# 
# # SELECT TRACKS YEAR
# 
# dupIDs <- dupDist <- length <- TRACKS_YEAR <- TRACKS_YEAR.sp <- list()
# for(t in 1:nYears){
#    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
#    TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[1]], ]
#    TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[2]], ]
#    TRACKS <- rbind(TRACKS_1, TRACKS_2)
#    ## SIMPLIFY TRACKS SHAPES
#    ## SUBSET TRACKS TO THE STUDY AREA
#    TRACKS <- st_intersection(TRACKS, st_as_sfc(myStudyArea.extent))
#    ## NAME TRACKS
#    TRACKS$ID <- 1:nrow(TRACKS)
#    TRACKS_YEAR[[t]] <- TRACKS
#    #calculate length to identify duplicates
#    TRACKS_YEAR[[t]]$dist <- st_length(TRACKS_YEAR[[t]],byid = T)
#    #calculate centroids to also identify tracks that could have the same length but in different location
#    TRACKS_YEAR[[t]]$centroidx <-  st_coordinates(st_centroid(TRACKS_YEAR[[t]]))[,1]
#    # try a fast way to identify duplicated tracks
#    # turn to dataframe and identify them
#    # sub_tracks_filter <- TRACKS_YEAR[[t]] %>%
#    #   distinct(Dato, dist, .keep_all = T)
#     # distinct(Person, Dato, dist, .keep_all = T)
#    df <- data.frame(ID = TRACKS_YEAR[[t]]$ID,
#                     Dato=TRACKS_YEAR[[t]]$Dato,
#                     Person =TRACKS_YEAR[[t]]$Person,
#                     dist=TRACKS_YEAR[[t]]$dist,
#                     centroidx = TRACKS_YEAR[[t]]$centroidx)
#    dupIDs[[t]] <- duplicated(df[,2:5])# find duplicates based on person and distance and date
#    dupIDs[[t]] <- df$ID[dupIDs[[t]]]
#    dupDist[[t]] <- TRACKS_YEAR[[t]][dupIDs[[t]],]$dist
#   #  plot(TRACKS_YEAR[[t]][dupIDs[[t]][1:4],]$geometry)
#   #
#   #  i=7
#   # whicdup <- which(TRACKS_YEAR[[t]]$Dato %in% df[dupIDs[[t]][i],]$Dato &
#   #          TRACKS_YEAR[[t]]$Person %in% df[dupIDs[[t]][i],]$Person &
#   #          TRACKS_YEAR[[t]]$dist %in% df[dupIDs[[t]][i],]$dist  &
#   #          TRACKS_YEAR[[t]]$centroidx %in% df[dupIDs[[t]][i],]$centroidx
#   #          )
#   # TRACKS_YEAR[[t]][whicdup,]
#   # st_coordinates(st_centroid(TRACKS_YEAR[[t]][whicdup,]))
#   # mapview::mapview(TRACKS_YEAR[[t]][whicdup,])
#   #
#   # mapview::mapview(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T451513",])
#    TRACKS_YEAR[[t]] <-  TRACKS_YEAR[[t]][-dupIDs[[t]],]
#    # tmp <- st_length(TRACKS_YEAR[[t]],byid = T)
#    # tmp1 <- st_length(TRACKS_YEAR[[t]], by_element = FALSE)
#    #
# 
# }#t
# 
# #check
# #number of tracks
# barplot(unlist(lapply(TRACKS_YEAR,function(x) sum(x$dist))),ylab="sum length tracks")
# #chekc number of duplicated tracks removed
# dup <- (unlist(lapply(dupIDs,length)))
# names(dup) <- years
# barplot(dup,ylab="Number of duplicated tracks")
# ##distance
# dupdist <- (unlist(lapply(dupDist,sum)))
# names(dupdist) <- years
# barplot(dupdist,ylab="Distance")



## ------     3.2. DISTANCE TO ROADS ------ 

## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
DistAllRoads <- raster(paste(dir.dropbox,"/DATA/GISData/Roads/MinDistAllRoads1km.tif", sep=""))
r <- fasterize(myStudyArea.poly, DistAllRoads)
r[!is.na(r)] <- DistAllRoads[!is.na(r)]
DistAllRoads <- r
DistAllRoads <- crop(DistAllRoads, myStudyArea.poly)

## PLOT CHECK
if(myVars$plot.check){
  plot((DistAllRoads))
  plot(st_geometry(myStudyArea.poly),add=T)
}



## ------     3.3. DAYS OF SNOW ------ 

## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
## UPDATE"!!!!!
SNOW <- stack(file.path(dir.dropbox,
                        "DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2014_2025_Wolf.tif"))

## RENAME THE LAYERS
names(SNOW) <- paste(2014:2024,(2014:2024)+1, sep="_")

## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
SNOW <- SNOW[[paste("X", years, "_", years+1, sep="")]]
SNOW <- raster::crop(SNOW, c(0,40,55,75))



## ------     3.4. SAVE SEARCH EFFORT OBJECTS FOR FASTER RUNS ------ 

save(TRACKS_YEAR,
     SNOW,
     DistAllRoads,
     file = file.path(myVars$WD, "TRACKSSouthSweden2014202540NotSimplifiedSF.RData"))
load(file.path(myVars$WD, "TRACKSSouthSweden2014202540NotSimplifiedSF.RData"))



## -----------------------------------------------------------------------------

## ------ II. CREATE OPSCR DATA ------

## ------   1. CLEAN & FILTER NGS DATA ------ 

## ------     1.1. CLEAN NGS & DEAD RECOVERY DATA ------ 

myCleanedData.sp <- CleanDataNew3sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = myVars$DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(myVars$DATA$samplingMonths)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)

##---- OVERWRITE GENDER FROM MICKE'S DATA WHEN AVAILABLE
micke.sex <- as.character(unlist(lapply(myCleanedData.sp$Id, function(i) INDIVIDUAL_ID[as.character(INDIVIDUAL_ID$Individ..Rovbase.)==i,"Sex"][1])))
micke.sex[micke.sex %in% "0"] <- NA
micke.sex[micke.sex %in% names(table(micke.sex))[3]] <- NA
micke.sex[micke.sex %in% "Hona"] <- "Hunn"
micke.sex[micke.sex %in% "Hane"] <- "Hann"
table(!is.na(micke.sex))
new.sex <- ifelse(!is.na(micke.sex), as.character(micke.sex), as.character(myCleanedData.sp$Sex))
table(myCleanedData.sp$Sex, new.sex)
myCleanedData.sp$Sex <- new.sex
table(myCleanedData.sp$Sex, new.sex)

##---- OVERWRITE GENDER FROM PACK COMPOSITION (FROM LINN's file 2022-23)
Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Tispe"] <- "Hunn"
Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Hane"] <- "Hann"
Pack_ID2023$Kon[Pack_ID2023$Kon %in% "Tik"] <- "Hunn"

Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Tispe"] <- "Hunn"
Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Hane"] <- "Hann"
Pack_ID2024$Kon[Pack_ID2024$Kon %in% "Tik"] <- "Hunn"

Pack_ID2025$Sex[Pack_ID2025$Sex %in% "Tispe"] <- "Hunn"
Pack_ID2025$Sex[Pack_ID2025$Sex %in% "Hane"] <- "Hann"
Pack_ID2025$Sex[Pack_ID2025$Sex %in% "Tik"] <- "Hunn"

##---- make a simplified column to match the rovbase id given by Oystein in Linn's file
myCleanedData.sp$IdSimplified <- unlist(lapply(strsplit(as.character(myCleanedData.sp$Id), " "), function(x) x[1]))
#myCleanedData.sp$Id <- as.character(myCleanedData.sp$Id)
##---- OVERWRITE GENDER FROM PACK COMPOSITION (FROM LINN's file 2023-24)



table(myCleanedData.sp$Sex)
tab <- list()

##check the sex in the pair data given by Linn and assign the sex to all detections 
for(i in 1:length(Pack_ID2023$Kon)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2023$Rovbase.ID[i]])
  #Overwrite sex 
  # if(length(tab[[i]])>1){print(tab[[i]])}
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2023$Rovbase.ID[i]  ] <- Pack_ID2023$Kon[i]
}
##
for(i in 1:length(Pack_ID2024$Kon)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2024$Rovbase.ID[i]])
  #Overwrite sex 
  # if(length(tab[[i]])>1){print(tab[[i]])}
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2024$Rovbase.ID[i]  ] <- Pack_ID2024$Kon[i]
}

##
i=130
for(i in 1:length(Pack_ID2025$Sex)){
  tab[[i]] <- table(myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2025$IndividID[i]])
  # if(length(tab[[i]])>1){print(tab[[i]])}
  #Overwrite sex 
  myCleanedData.sp$Sex[myCleanedData.sp$IdSimplified %in% Pack_ID2025$IndividID[i]  ] <- Pack_ID2025$Sex[i]
}


## ------     1.2. FILTER DATA FOR SEX myFullData.sp====
## ONLY FILTER FOR ONE SEX
myFullData.sp <- FilterDatasf( myData = myCleanedData.sp
                               ,
                               dead.recovery = T
                               ,
                               sex = myVars$DATA$sex
                               , 
                               setSex = T
)


## ------     1.3. REMOVE INDVIDUALS THAT DIED TWICE ------ 
duplicatedDeath <- NULL
for(i in myFullData.sp$IdDoubleDead){
  tmp  <- which(myFullData.sp$dead.recovery$Id == i & is.na(myFullData.sp$dead.recovery$DeathCause_2))
  if(length(tmp)==0){tmp  <- which(myFullData.sp$dead.recovery$Id == i)[-1]}
  duplicatedDeath <- c(duplicatedDeath, tmp)
}#i  
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-duplicatedDeath, ]

##EXPORT THE DATA 
if(myVars$DATA$sex=="Hann"){
  assign("myFullData.spM", myFullData.sp)
}else{
  assign("myFullData.spF", myFullData.sp)
}
## ------     1.4. FILTER DATA FOR SPACE ------ 
myFilteredData.sp <- myFullData.sp

## Remove all alive detections outside of the study area extent 
myFilteredData.sp$alive <- myFilteredData.sp$alive[!is.na(as.numeric(st_intersects(myFilteredData.sp$alive, st_as_sfc(myStudyArea.extent)))), ]


## Remove all dead recoveries outside of the study area polygon 
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[!is.na(as.numeric(st_intersects(myFilteredData.sp$dead.recovery, st_as_sfc(myStudyArea.extent)))), ]

## ------     1.5. FILTER DATA FOR DATES ------ 
## Remove all alive detections outside of the sampling period 
myFilteredData.sp$alive <- myFilteredData.sp$alive[ myFilteredData.sp$alive$Month %in% unlist(myVars$DATA$samplingMonths) 
                                                    & myFilteredData.sp$alive$Year %in% unlist(myVars$DATA$years), ]

## Remove all dead recoveries outside of the sampling period 
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[ myFilteredData.sp$dead.recovery$Year %in% unlist(myVars$DATA$years),] 

## ------     1.6. SEPARATE STRUCTURED AND OPPORTUNISTIC SAMPLING ------ 
## ------       1.6.1. ASSIGN SAMPLES TO TRACKS  ------ 
## ASSIGN ROVBASE ID AND SIMPLIFY TRACKS
myFilteredData.sp$alive$TrackRovbsID <- NA
myFilteredData.sp$alive$TrackDist <- NA


#TRACKS_YEAR <- tmpT1[c(2:nYears,1)]#tmp1[c(2:nYears,1)]

TRACKSSimple_sf <- list()
for(t in 1:nYears){
  TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbaseID)
  TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]
}

## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
dnatemp <- st_as_sf(myFilteredData.sp$alive)
#CREATE A BUFFER AROUND EACH DETECTION
tmp <-  st_buffer(dnatemp, dist=750)


## i=349
##tmpT[[5]][tmpT[[5]]$RovbaseID %in% "T442459",]$Dato
##tmpT1[[6]][tmpT1[[6]]$RovbaseID %in% "T442459",]$Dato



#
for(i in 1:nrow(myFilteredData.sp$alive)){
  # INTERSECT POINT WITH TRACKS,
  t <- which(years %in% tmp[i,]$Year)
  
  
  whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
  
  tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])
  
  # plot(COUNTIES$geometry, border=NA,col=grey(0.8))
  # plot(TRACKSSimple_sf[[t]][whichSameDate,],add=T,col="red")
  # plot(tmp[i,]$geometry,add=T,col="green", cex=10)
  # dnatemp[i,]$Annen.innsamler...Rolle
  
  if(nrow(tmpTRACKS)==0){next}
  # FIND THE CLOSEST TRACK
  dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)
  
  #whichSameDate <- numeric()
  # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
  #whichSameDate <- which(as.character(tmpTRACKS$Dato)==as.character(myFilteredData.sp$alive$Date[i]))
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

# #SAVE THE FOR FASTER LOADING
save(myFilteredData.sp, file = file.path(myVars$WD, myVars$modelName,
                                         paste("_myFilteredData.sp", ".RData", sep = "")))
load(file.path(myVars$WD, myVars$modelName,paste("_myFilteredData.sp", ".RData", sep = "")))


## ------       1.6.2. SPLIT MYFILTERED DATA TO OPPORTUNISTIC AND STRUCTURED ------ 
distanceThreshold <- 500
#Proevetype columns was replaced by two columns, merging them now...
myFilteredData.sp$alive$Proevetype <-  ifelse(myFilteredData.sp$alive$Annen.innsamler...Rolle %in% "" , 
                                              myFilteredData.sp$alive$Samlet.selv...Rolle, myFilteredData.sp$alive$Annen.innsamler...Rolle)

table(myFilteredData.sp$alive$Proevetype)
whichStructured <- myFilteredData.sp$alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  !is.na(myFilteredData.sp$alive$TrackRovbsID) &
  myFilteredData.sp$alive$TrackDist <= distanceThreshold
myFilteredData.spStructured <- myFilteredData.sp$alive[whichStructured,]
myFilteredData.spOthers <- myFilteredData.sp$alive[!whichStructured,]

## CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
nrow(myFilteredData.spStructured) +
  nrow(myFilteredData.spOthers)
nrow(myFilteredData.sp$alive)


myFilteredData.sp$alive$TrackDistCat <- ifelse(myFilteredData.sp$alive$TrackDist>500,0,1)
myFilteredData.sp$alive$TrackDistCat[is.na(myFilteredData.sp$alive$TrackDistCat)] <- 0

table(myFilteredData.sp$alive$Proevetype,myFilteredData.sp$alive$Year,myFilteredData.sp$alive$TrackDistCat)


# tmp <- myFilteredData.spStructured[myFilteredData.spStructured$Year %in% 2013,]
# tmp <- tmp[tmp$Month %in% 1,]
# tmp$Id
# 
# tmp[tmp$Id%in%"UI040578 G85-13",]
# tmp[tmp$Id%in%"UI040754 G83-13",]$DNAID
# 
# hist(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% 2023, ]$Annen.innsamler...Rolle)
# 
# 
# table(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% 2023, ]$Annen.innsamler...Rolle)

## ------       1.6.3. PLOT CHECKS   ------ 
## PROPORTION SAMPLES STRUCTURED/OTHERS
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"ProportionStucturedOther",".pdf", sep="" ))
par(mfrow=c(1,1),mar=c(4,4,3,2))
barplot(rbind(table(myFilteredData.spStructured$Year),
              table(myFilteredData.spOthers$Year)),beside=T,ylim=c(0,2000),col=c(grey(0.2),grey(0.8)),ylab="Number of samples")
abline(h=seq(0,2000,by=500),lty=2,col=grey(0.8))
title(main="500m threshold")
legend("topleft",fill=c(grey(0.2),grey(0.8)),legend=c("Structured","Other"))
dev.off()

## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO" 
tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Proevetype %in% 
                                 c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
tab <- table(tmp$Year, tmp$TrackRovbsID, useNA ="always" )


## MAP  SAMPLES STRUCTURED OTHERS
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"MapStucturedOther",".pdf", sep="" ))
for(t in 1:nYears){
  par(mar=c(0,0,3,0),mfrow=c(1,3))
  tmp1 <- tmp[tmp$Year%in% years[t],]
  tmpNoTracks <-  tmp1[is.na(tmp1$TrackRovbsID), ]
  tmpTracks <-  tmp1[!is.na(tmp1$TrackRovbsID), ]
  
  plot(st_geometry(myStudyArea.poly), main="Structured with track")
  plot(st_geometry(tmpTracks), pch=21, col="black", cex=1,bg="red",add=T)
  
  plot(st_geometry(myStudyArea.poly), main="Structured without track")
  plot(st_geometry(tmpNoTracks), pch=21, col="black", cex=1,bg="blue",add=T)
  
  tmpOpp <- myFilteredData.sp$alive[!myFilteredData.sp$alive$Proevetype %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),]
  tmpOpp <- tmpOpp[tmpOpp$Year%in% years[t],]
  
  plot(st_geometry(myStudyArea.poly), main="Other samples")
  plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
  
  mtext(years[t],adj = -0.8,padj = 1)
  
}
barplot(tab[,which(is.na(colnames(tab)))]/rowSums(tab),main="% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track") 

dev.off()

### OVERALL MAP DETECTION DEAD RECOVERIES MAP
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"OverallDetectionsDeadRecoveries",".pdf", sep="" ))
plot(st_geometry(GLOBALMAP))
plot(st_geometry(myStudyArea.poly),add=T)
plot(st_geometry(myFullData.sp$alive), pch=16, col="red", cex=0.3,add=T)
plot(st_geometry(myFullData.sp$dead.recovery),pch=16, col="blue", cex=0.3,add=T)
mtext(paste("Live detections", length(myFullData.sp$alive),
            "; ID:", length(unique(myFullData.sp$alive$Id))
),line = +1)
mtext(paste("Dead recovery:",length(myFullData.sp$dead.recovery)))
dev.off()



## ------   2. GENERATE HABITAT ------ 
## ------     2.1. GENERATE HABITAT CHARACTERISTICS ------ 
myHabitat.list <- MakeHabitatFromRastersf( poly = myStudyArea.poly
                                           ,
                                           habitat.r = habitatRasters[["Habitat"]]
                                           ,
                                           buffer = myVars$HABITAT$habBuffer
                                           ,                               
                                           plot.check = T
)

nHabCells <- sum(myHabitat.list$habitat.r[ ]==1)

## PLOT CHECK
if(myVars$plot.check){
  plot(myHabitat.list$habitat.r)
  plot(st_geometry(myStudyArea.poly), add=T)
  plot(st_geometry(COUNTRIESWaterHumans), add=T)
}

## ------   3. GENERATE DETECTORS ------ 
## ------     3.1. GENERATE DETECTORS CHARACTERISTICS ------ 
habitat.subdetectors <- disaggregate(myHabitat.list$habitat.rWthBuffer, fact=res(myHabitat.list$habitat.r)[1]/myVars$DETECTORS$detSubResolution)

myDetectors <- MakeSearchGridsf(data = habitat.subdetectors,
                                resolution = myVars$DETECTORS$detResolution,
                                div = (myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)^2,      
                                plot = FALSE,
                                fasterize = TRUE)
# proj4string(myDetectors$main.detector.sp) <- CRS(proj4string(myHabitat.list$habitat.sp))
plot(st_geometry(myHabitat.list$buffered.habitat.poly))
plot(st_geometry(myDetectors$main.detector.sp), pch=16, cex=0.1,add=T)

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(myDetectors$main.detector.sp)[1]

## FORMAT DETECTOR LOCATIONS and NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(myDetectors$main.detector.sp)
n.trials <- as.vector(table(myDetectors$detector.sp$main.cell.id))

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,2))
  plot(st_geometry(myStudyArea.poly), main = "Detectors Alive")
  plot(st_geometry(myDetectors$main.detector.sp), col = "red", pch = 16, cex = 0.1, add = T)
  plot(st_geometry(GLOBALMAP), add = T)
  
  plot(st_geometry(myStudyArea.poly), main = "Detectors Dead")
  plot(st_geometry(myDetectors$main.detector.sp), col = "red", pch = 16, cex = 0.1, add = T)
  plot(st_geometry(GLOBALMAP), add = T)
}

## ------     3.2. GENERATE DETECTOR-LEVEL COVARIATES ------ 
## ------       3.2.1. EXTRACT COUNTRIES ------ 
dist <- st_distance(myDetectors$main.detector.sp, country, by_element = F )
detCountries <- apply(dist,1, function(x) which.min(x))
detCountries <- as.numeric(as.factor(detCountries))

## ------       3.2.2. EXTRACT COUNTIES ------ 
dist <- st_distance(myDetectors$main.detector.sp, COUNTIES, by_element = F )

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
COUNTIES$id[c(31,37,26,27)] <- 4
#SE3
COUNTIES$NAME_1[c(3,21,40,13,15,14,24,7)]
COUNTIES$id[c(3,21,40,13,15,14,24,7)] <- 3
#SE4
COUNTIES$NAME_1[c(38,34,9)]
COUNTIES$id[c(38,34,9)] <- 9

#CONVERT TO FACTOR AND BACK
COUNTIES$id <- as.numeric(as.factor(COUNTIES$id))

## ASSIGN COUNTIES TO DETECTORS.
detCounties1 <- apply(dist,1, function(x) which.min(x))
detCounties1 <- COUNTIES$id[detCounties1]
detCounties <- as.numeric(as.factor(detCounties1))
table(detCounties)
# CREATE A VECTOR TO KEEP THE ORIGINAL COUNTY ID
detCounties.original <- 0
for(i in 1: max(detCounties)){
  detCounties.original[i] <-   detCounties1[which(detCounties==i)][1]
}

##PLOT CHECK 
COUNTIESplot <- st_intersection(st_simplify(COUNTIES,dTolerance = 500),myStudyArea.poly)
COUNTIESplot <- COUNTIESplot %>% group_by(id) %>% summarize()

ggplot(st_simplify(COUNTIESplot,dTolerance = 500)) +
  geom_sf(aes(fill = id)) +
  geom_sf_label(aes(label = id))

col <- rainbow(length(unique(detCounties)))
plot(st_geometry(COUNTIESplot))
plot(st_geometry(myDetectors$main.detector.sp), col=col[detCounties], pch=16, cex=0.8,add=T)


## ------       3.2.3. EXTRACT GPS TRACKS LENGTHS ------ 
## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
# detectorGrid <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp), rep(1,nrow(myDetectors$main.detector.sp))))
# detectorGrid <- rasterToPolygons(detectorGrid)
# detectorGrid$id <- 1:nrow(detectorGrid)
# plot(detectorGrid)
# proj4string(detectorGrid) <- CRS(proj4string(myDetectors$main.detector.sp))
# plot(detectorGrid)
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp),
                                      rep(1,nrow(myDetectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf(stars::st_as_stars(detectorGrid.r), 
                             as_points = FALSE, merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea.poly)
detectorGrid$id <- 1:nrow(detectorGrid)
plot(st_geometry(detectorGrid))



##CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
for(t in 1:nYears){
  TRACKSst <- TRACKS_YEAR[[t]]
  
  intersection <- st_intersection(detectorGrid, TRACKSst) %>%
    dplyr::mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(transect_L = sum(LEN))#,               ## Get total length searched in each detector grid cell
  # transect_N = length(unique(ID)))#,   ## Get total number of visits in each detector grid cell
  #transect_qi = mean(QI))              ## Get mean transects quality index for each detector grid cell
  
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
}

#check 
t=1
par(mar=c(0,0,0,0))
plot(st_geometry(myDetectors$main.detector.sp[as.numeric(intersection$id),]),
     cex=DoScale(detTracks[,t],l = 0,u = 2),pch=16)


## ------       3.2.4. EXTRACT DISTANCES TO ROADS ------ 
## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)
if(myVars$plot.check){
  plot(st_geometry(myDetectors$main.detector.sp),cex=DoScale(detRoads),pch=16)
}
isna <- which(is.na(detRoads))
tmp <- raster::extract(DistAllRoads, myDetectors$main.detector.sp[isna, ], buffer = 15000, fun = mean, na.rm = T)
detRoads[isna] <- tmp
if(myVars$plot.check){
  plot(st_geometry(myDetectors$main.detector.sp),cex=DoScale(detRoads),pch=16)
}

## ------       3.2.5. EXTRACT DAYS OF SNOW ------ 
## EXTRACT SNOW 
detSnow <- matrix(0, nrow = n.detectors, ncol = nYears)
det.sptransf <- st_transform(myDetectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 20000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract(SNOW, det.sptransf[isna, ], buffer = 15000, fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp
if(myVars$plot.check){
  plot(st_geometry(myDetectors$main.detector.sp), cex=DoScale(detSnow[,3]),pch=16)
}
# still some NA... Increase buffer again 
isna <- which(is.na(detSnow),arr.ind = T)
isna <- unique(isna[,1])
tmp.list <- raster::extract(SNOW, det.sptransf[isna,], buffer=35000)
detSnow[isna,1:nYears] <- unlist(lapply(tmp.list, function(x) colMeans(x, na.rm=T)))

## ------       3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------ 
## ------         3.2.6.1. SKANDOBS ------ 
skandObs <- read_xlsx(file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/Richard_Biscof_Skandobs_2012_2025dd.xlsx"))
colnames(skandObs) <- translateForeignCharacters(dat=colnames(skandObs), dir.translation = dir.analysis )
## GET TIME 
skandObs$date1 <- as.POSIXct(strptime(skandObs$date, "%Y-%m-%d"))
skandObs$year <- as.numeric(format(skandObs$date1,"%Y"))
skandObs$month <- as.numeric(format(skandObs$date1,"%m"))
## MAKE IT SPATIAL 
skandObs <- st_as_sf(skandObs, coords = c("longitude", "latitude"))
st_crs(skandObs) <- st_crs("EPSG:4326")
skandObs <- st_transform(skandObs, st_crs(myStudyArea.poly))
## SUBSET BASED ON SEASON 
subset <- skandObs$month %in% c(unlist(myVars$DATA$samplingMonths))
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

## PLOT CHECK 
## SUMMARY SKANDOBS
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"skandObs",".pdf", sep="" ), width = 10)
barplot(table(skandObs$monitoring.season ))
barplot(table(skandObs$month ), xlab="Months")
barplot(table(skandObs$activity),cex.names=0.7)
barplot(table(skandObs$species))
## MAPS 
par(mar=c(0,0,2,0))
for(t in 1:nYears){
  plot(st_geometry(myStudyArea.poly), main= years[t])
  plot(st_geometry(skandObs[skandObs$monitoring.season %in% years[t],  ]), pch=16, col="red", cex=0.1)
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


## ------         3.2.6.2. ROVBASE ------ 
## GET ALL SAMPLES COLLECTED
#files need to be loaded separately and binded. 
rovbaseObs1 <- read_xlsx(file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/RIB15042025111042561.xlsx"))
rovbaseObs2 <- read_xlsx(file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/RIB15042025112748412.xlsx"))
rovbaseObs3 <- read_xlsx(file.path(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20250415/RIB15042025112853603.xlsx"))
rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3)

colnames(rovbaseObs) <- translateForeignCharacters(dat=colnames(rovbaseObs), dir.translation = dir.analysis )
rovbaseObs$Proevetype <- translateForeignCharacters(dat=rovbaseObs$Proevetype, dir.translation = dir.analysis )

rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$`Oest (UTM33/SWEREF99 TM)`),]
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

rovbaseObs.sp <- rovbaseObs
#---DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf(rovbaseObs, coords = c("Oest (UTM33/SWEREF99 TM)","Nord (UTM33/SWEREF99 TM)"))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea.poly)

# SUBSET THE DATA 
filter <- list(
  species = "Ulv"
  , type = c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)","Saliv/Spytt")
  , month = unlist(myVars$DATA$samplingMonths)
)

### SUBSET MONTH AND TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset,] 
### SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
subset <- rovbaseObs.sp$`Art (Proeve)`%in% filter$species #& !is.na(rovbaseObs.sp$`RovbaseID (Analyse)`) 
rovbaseObs.sp <- rovbaseObs.sp[subset,] 
## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace,] 


# RASTERIZE 
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


#PLOT
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"mapStructuredOthers",".pdf", sep="" ))
for(t in 1:nYears){
  year= years[t]
  tmpOthers <- myFilteredData.spOthers[myFilteredData.spOthers$Year%in%year, ]
  tmpStruct <- myFilteredData.spStructured[myFilteredData.spStructured$Year%in%year, ]
  
  par(mfrow=c(2,2),mar=c(0,0,5,0))
  plot(r.OtherSamplesBinary[[t]], main=paste(year,"\n Rovbase Samples Other"), box=F, axes=F)
  plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  plot(r.OtherSamplesBinary[[t]],main=paste(year,"\n Rovbase Samples Structured"), box=F, axes=F)
  plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
  
  plot(r.skandObsSamplesBinary[[t]], main=paste(year,"\n SkandObs Other"), box=F, axes=F)
  plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  plot(r.skandObsSamplesBinary[[t]],main=paste(year,"\n SkandObs Structured"), box=F, axes=F)
  plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
}
dev.off()


## ------         3.2.6.3. COMBINE ROVBASE AND SKANDOBS ------ 
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
par(mfrow=c(1,3))
plot(r.SkandObsOtherSamplesBinary[[t]], main="Raw Binary",axes=F,box=F)
plot(ds.brick[[t]], main="Smoothed",axes=F,box=F)
plot(ds.brickCont[[t]], main="Binary after smoothing",axes=F,box=F)

## ------         3.2.6.5. ASSIGN THE COVARIATE ------ 
detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract(r.SkandObsOtherSamplesBinary, myDetectors$main.detector.sp)
colSums(detOtherSamples)



## ------       3.2.7. SCALE AND ROUND DETECTOR-LEVEL COVARIATES ------ 
detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

# CHECK IF CONTAINS NAs
if(sum(is.na(detSnow))>0 | sum(is.na(detRoads))>0| sum(is.na(detTracks))>0 ){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}


## PLOT COVARIATES ##
#Tracks
pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"SpatialCovariates",".pdf", sep="" ))

max <- max(detTracks)
cuts <- seq(min(detTracks),max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

for(t in 1:nYears){
  detectorGrid <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp), rep(1,nrow(myDetectors$main.detector.sp))))
  id <- which(detectorGrid[]%in%1)
  detectorGrid[id] <- detTracks[,t]
  #plot(TRACKS.r[[t]], main=paste("Tracks", years[t]))
  plot(detectorGrid, breaks=cuts, col = col, main=paste("Tracks", years[t]),legend=F) #p
  plot(detectorGrid, legend.only=TRUE,breaks=cuts, col=col,
       legend.width = 2,
       axis.args=list(at=round(seq(min(detTracks), max, length.out = 8),digits = 0),
                      labels=round(seq(min(detTracks), max, length.out = 8),digits = 0), 
                      cex.axis=0.6),
       legend.args=list(text='', side=4, font=2, line=2.5, cex=0.8))
  plot(st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t],]), pch=16, col="red", cex=0.2,add=T)
  
  plot(st_geometry(myHabitat.list$buffered.habitat.poly), add=T,border="grey")
}
#snow
for(t in 1:nYears){
  plot(SNOW[[t]], main=paste("Snow", years[t]))
  plot(st_transform(myHabitat.list$buffered.habitat.poly, st_crs(SNOW))$geometry, add=T)
}
#roads
plot(DistAllRoads, main="Roads")
plot(st_geometry(myHabitat.list$buffered.habitat.poly), add=T)
dev.off()



## ------   4. GENERATE y DETECTION ARRAYS ------ 
## ------     4.1. DATA SUMMARY ------ 
## ------       4.1.1. DETECTIONS ALIVE ------ 
## NUMBER OF INDIVIDUALS DETECTED ALIVE
length(unique(myFilteredData.sp$alive$Id))

## NUMBER OF INDIVIDUALS DETECTED ALIVE AND RECOVERED DEAD
sum(unique(myFilteredData.sp$alive$Id) %in% unique(myFilteredData.sp$dead.recovery$Id))

## NUMBER OF INDIVIDUALS DETECTED/YEAR/COUNTRY
table.id <- table(myFilteredData.sp$alive$Id,myFilteredData.sp$alive$Year, myFilteredData.sp$alive$Country)
apply(table.id, c(2,3), function(x) sum(x>0))

## NUMBER OF DETECTIONS/YEAR/COUNTRY
countrytab <- table(myFilteredData.sp$alive$Year, myFilteredData.sp$alive$Country)
countrytab 


## NUMBER OF DETECTIONS/YEAR/COUNTRY/SEX
sex_countrytab <- table(myFilteredData.sp$alive$Year, myFilteredData.sp$alive$Country, myFilteredData.sp$alive$Sex)
sex_countrytab

## PROPORTION OF DETECTIONS PER YEAR/COUNTRY/SEX
sex_countrytab_prop <- sex_countrytab
for(i in 1:dim(sex_countrytab)[3]){
  sex_countrytab_prop[,,i] <- sex_countrytab_prop[,,i]/countrytab
}
sex_countrytab_prop

## ------       4.1.2. DEAD RECOVERIES ------ 
## NUMBER OF INDIVIDUALS RECOVERED
length(unique(myFilteredData.sp$dead.recovery$Id))

## NUMBER OF DEAD RECOVERIES/YEAR/COUNTRY
table(myFilteredData.sp$dead.recovery$Year, myFilteredData.sp$dead.recovery$Country)

## MORTALITY CAUSES
unique(as.character(myFilteredData.sp$dead.recovery$DeathCause))
unique(as.character(myFilteredData.sp$dead.recovery$DeathCause_2))  

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


if(myVars$plot.check){
  pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"DetectionsDeadRecoveries",".pdf", sep="" ))
  
  ## PLOT ALL DETECTIONS
  for(t in 1:nYears){
    plot(myHabitat.list$habitat.r, main = years[t]) 
    plot(st_geometry(myDetectors$main.detector.sp), add=T, pch=16, cex=0.1)
    plot(st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]), 
         pch=16,col="red", cex=0.7, add=T)
    mtext(paste("Live detections", nrow(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]),
                "; ID:", nrow(unique(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ]$Id))
    ),line = +1)
    
    plot(st_geometry(myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ]),
         pch=16,col="blue", cex=0.7,add=T)
    mtext(paste("Dead recovery:",nrow(myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ])))
    plot(st_geometry(GLOBALMAP), add=T) 
  }#t
  
  ## PLOT NUMBER OF MORTALITY EVENTS, LEGAL/OTHERS 
  #par(mfrow=c(1,2))
  barplot(table(legal.death$Country,legal.death$Year ), main="Legal causes")
  barplot(table(Other.death$Country,Other.death$Year ), main="Other causes")
  legend("topleft", fill=c(grey(0.3), grey(0.7)), legend=c("NOR","SWE"))
  dev.off()  
}#if plot.check

##### EXPORT NGS DATA 
if(myVars$DATA$sex=="Hann"){
  assign("myFilteredData.spM", myFilteredData.sp)
  assign("myFilteredData.spOthersM", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredM", myFilteredData.spStructured)
  
  save(myFilteredData.spM, myFullData.spM,
       myFilteredData.spOthersM,myFilteredData.spStructuredM,
       file = file.path(myVars$WD, myVars$modelName,
                        paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
  
}else{
  assign("myFilteredData.spF", myFilteredData.sp)
  assign("myFilteredData.spOthersF", myFilteredData.spOthers)
  assign("myFilteredData.spStructuredF", myFilteredData.spStructured)
  
  save(myFilteredData.spF, myFullData.spF,
       myFilteredData.spOthersF,myFilteredData.spStructuredF,
       file = file.path(myVars$WD, myVars$modelName,
                        paste(myVars$modelName, "_NGSData", ".RData", sep = "")) )
  
}
## ------     4.2. ASSIGN DETECTORS ------ 
#ALL SAMPLES
myData.alive <- AssignDetectors_v3sf( myData = myFilteredData.sp$alive,                
                                      myDetectors = myDetectors$main.detector.sp,
                                      mysubDetectors = myDetectors$detector.sp,
                                      radius = myVars$DETECTORS$detResolution)
#STRUCTURED
myData.aliveStruc <- AssignDetectors_v3sf( myData = myFilteredData.spStructured,                
                                           myDetectors = myDetectors$main.detector.sp,
                                           mysubDetectors = myDetectors$detector.sp,
                                           radius = myVars$DETECTORS$detResolution)
#OTHERS
myData.aliveOthers <- AssignDetectors_v3sf( myData = myFilteredData.spOthers,                
                                            myDetectors = myDetectors$main.detector.sp,
                                            mysubDetectors = myDetectors$detector.sp,
                                            radius = myVars$DETECTORS$detResolution)


myData.dead <- AssignDetectors_v3sf( myData = myFilteredData.sp$dead.recovery,
                                     myDetectors = myDetectors$main.detector.sp,
                                     radius = myVars$DETECTORS$detResolution)

## ------     4.3. GENERATE NGS & DEAD RECOVERIES : y.alive[i,j,t] & y.dead[i,t] ------ 
#ALL SAMPLES
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 myDetectors = myDetectors$main.detector.sp,
                 method = "Binomial",
                 myData2 = myData.dead,
                 myDetectors2 = myDetectors$main.detector.sp,
                 returnIdvector = TRUE)
y.ar.ALIVE <- y.ar$y.ar
dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)

#STRUCTURED
y.arStruc <- MakeYsf( myData = myData.aliveStruc$myData.sp,
                      myDetectors = myDetectors$main.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      myDetectors2 = myDetectors$main.detector.sp,
                      returnIdvector = TRUE)
y.ar.ALIVEStruc <- y.arStruc$y.ar
dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)
#OTHERS
y.arOth <- MakeYsf( myData = myData.aliveOthers$myData.sp,
                    myDetectors = myDetectors$main.detector.sp,
                    method = "Binomial",
                    myData2 = myData.dead,
                    myDetectors2 = myDetectors$main.detector.sp,
                    returnIdvector = TRUE)
y.ar.ALIVEOth <- y.arOth$y.ar
dimnames(y.ar.ALIVEOth) <- dimnames(y.arOth$y.ar)

### MAKE SURE THE Y HAVE THE SAME DIMENSIONS#
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

##FOR THE LAST YEAR GET THROUGH THE PAIR BASED-FILE FROM LINN.
indIDRovBase <- unlist(lapply(strsplit(y.ar$Id.vector," "), function(x) x[1]))


ALLIDS <- c(unique(indIDRovBase))#, unique(INDIVIDUAL_ID$Individ..Rovbase.))




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


#sum(unlist(lapply(strsplit(as.character(myFilteredData.sp$alive$Id)," "), function(x) x[1])) %in% "UI418704")

##FOR THE LAST YEAR GET TRHOUGH THE PAIR BASED-FILE FROM LINN.
for(i in 1:dim(y.obsALL)[1]){
  # for(t in 1:(nYears+1)){
  t <- nYears
  #linn's 2023 file
  tmp <-  Pack_ID2023[Pack_ID2023$Rovbase.ID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    y.obsALL[i,t-1] <- 3
  }
  #linn's 2024 file
  
  tmp <-  Pack_ID2024[Pack_ID2024$RovbaseID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    #print(i)
    
    y.obsALL[i,t] <- 3
  }
  
  tmp <-  Pack_ID2025[Pack_ID2025$IndividID %in% ALLIDS[i],]
  if(nrow(tmp) > 0){
    #print(i)
    
    y.obsALL[i,t+1] <- 3
  }
  
}#i

# sum(ALLIDS %in% "UI418401")
# sum(ALLIDS %in% "UI418834")
apply(y.obsALL,2, table)


##subset y.obs for the individuals present in y.ar
indID <- unlist(lapply(strsplit(y.ar$Id.vector, " "), function(x)x[1]))
y.obs <- y.obsALL[indID,]


y.obs.family <- y.obs
y.obs[y.obs==4] <- 3
y.obs <- y.obs[,as.character(years)]

# ## PLOT CHECK 
if(myVars$plot.check){
  pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"States",".pdf", sep="" ))
  
  #apply(y.obs,MARGIN = 2,table)
  
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
  distances[[t]] <- CheckDistanceDetectionsV2sf( y = y.ar.ALIVE[,,t], 
                                                 detector.xy = st_coordinates(myDetectors$main.detector.sp), 
                                                 max.distance = myVars$DETECTIONS$maxDist,
                                                 method = "pairwise",
                                                 plot.check = F)
  
  # PLOT INDIVIDUALS THAT DO HAVE DETECTIONS FURTHER AWAY THAN THRESHOLD DISTANCE
  if(myVars$plot.check){
    par(mfrow=c(1,1),mar=c(1,1,1,1))
    if(sum(distances[[t]]$y.flagged) > 0){
      affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      for(i in affected.ids){
        plot(st_geometry(myStudyArea.poly), main = paste("t: ",t,"     i: ", i, sep = ""))
        plot(st_geometry(GLOBALMAP), add = T)
        plot(st_geometry(myDetectors$main.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == y.ar$Id.vector[i], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1)
        
        arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2], length = 0.1, lwd = 1)
        
        plot(st_geometry(myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
        
        tmp2 <- myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
        plot(st_geometry(tmp2), col = "blue", pch = 13, cex = 1.5, lwd = 1,add=T)
      }#i
    }#if
  }#if plot.check
  
  ## REMOVE DETECTIONS THAT ARE FURTHER THAN  THE THRESHOLD
  y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEOthers[,,t] <- y.ar.ALIVEOthers[,,t] * (1-distances[[t]]$y.flagged)
  y.ar.ALIVEStructured[,,t] <- y.ar.ALIVEStructured[,,t] * (1-distances[[t]]$y.flagged)
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]),]
  }
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
if(myVars$plot.check){image(t(already.detected))}

## ------     5.3. TELEPORTATION COVARIATE ------ 
## CHECK DISTANCES BETWEEN DETECTIONS BETWEEN YEARS 
# distancesACs <- CheckDistanceACS( y = y.ar.ALIVE + y.ar.DEADProjected,
#                                   detector.sp = myDetectors$main.detector.sp,
#                                   myStudyArea.sp = myStudyArea.poly,
#                                   plot.check = myVars$plot.check)
# 
# long.dispersal.id <- which(distancesACs$distances.not.consecutive > myVars$DETECTIONS$maxDispersalDist, arr.ind = T)
# dispersalToggle <- matrix(0, nrow = dim(y.ar.ALIVE)[1], ncol = nYears)
# for(i in 1:nrow(long.dispersal.id)){
#   dispersalToggle[long.dispersal.id[i,1],long.dispersal.id[i,3]] <- 1
# }   


## ------     5.4. AGE ------ 
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
y.alive <- MakeAugmentation(y = y.ar.ALIVE, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveOthers <- MakeAugmentation(y = y.ar.ALIVEOthers, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveStructured <- MakeAugmentation(y = y.ar.ALIVEStructured, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

y.dead <- MakeAugmentation(y = y.ar.DEAD, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.obs <- MakeAugmentation(y = y.obs, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 1)

## INDIVIDUAL COVARIATES
indSocialState <- MakeAugmentation(y = indSocialState, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 1)
already.detected <- MakeAugmentation(y = already.detected, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
age <- MakeAugmentation(y = age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
min.age <- MakeAugmentation(y = min.age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
precapture <- MakeAugmentation(y = precapture, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

## -----------------------------------------------------------------------------

## ------ III. MODEL SETTING & RUNNING ------- 
## ------   1. NIMBLE MODEL DEFINITION ------ 
modelCode <- nimbleCode({
  ##--------------------------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------##  
  for(st in 1:2){
    dmean[st] ~ dunif(0,100)
    lambda[st] <- 1/dmean[st]
  }
  beta.dens ~ dnorm(0.0,0.01)
  
  for(t in 1:n.years){
    habIntensity[1:numHabWindows,t] <- exp(beta.dens * habDens[1:numHabWindows,t])
    sumHabIntensity[t] <- sum(habIntensity[1:numHabWindows,t])
    logHabIntensity[1:numHabWindows,t] <- log(habIntensity[1:numHabWindows,t])
    logSumHabIntensity[t] <- log(sumHabIntensity[t])
  }
  
  for(i in 1:n.individuals){
    sxy[i, 1:2, 1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows,1],
      logSumIntensity = logSumHabIntensity[1],
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
        lambda = lambda[state[i,t-1]+1]#dispSigma[state[i,t-1]+1]
        ,
        baseIntensities = habIntensity[1:numHabWindows,t]
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
  
  ##--------------------------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##    
  omeg1[1:3] ~ ddirch(alpha[1:3])   
  
  # LEGAL HUNTING (TIME-DEPENDENT)[RB]
  for(t in 1:n.years1){
    w[1,t] ~ dunif(0,1)
    w[2,t] ~ dunif(0,1)
    h[1,t] ~ dunif(0,1)
    h[2,t] ~ dunif(0,1)
    rw[1,t] ~ dunif(0,1)
    rw[2,t] ~ dunif(0,1)
    
    ones.dead.legal[1,t] ~ dbern(step(1 - (h[1,t] + w[1,t] + rw[1,t])))      
    ones.dead.legal[2,t] ~ dbern(step(1 - (h[2,t] + w[2,t] + rw[2,t])))  
  } 
  
  # FIRST YEAR
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:3]) 
    
  }#i
  
  for(t in 1:n.years1){
    # PRIORS 
    gamma[t] ~ dunif(0,1)
    phi[1,t] <- 1-h[1,t]-w[1,t]-rw[1,t]
    phi[2,t] <- 1-h[2,t]-w[2,t]-rw[2,t]
    
    wAll[1,t] <- w[1,t]+rw[1,t]
    wAll[2,t] <- w[2,t]+rw[2,t]
    
    psi[t] ~ dunif(0,1)
    
    # "UNBORN"
    omega[1,1,t] <- 1-gamma[t]
    omega[1,2,t] <- gamma[t]
    omega[1,3,t] <- 0
    omega[1,4,t] <- 0
    omega[1,5,t] <- 0
    omega[1,6,t] <- 0
    
    # "NON-PAIRS"
    omega[2,1,t] <- 0
    omega[2,2,t] <- phi[1,t]*(1-psi[t])
    omega[2,3,t] <- phi[1,t]*psi[t]
    omega[2,4,t] <- h[1,t]
    omega[2,5,t] <- rw[1,t]
    omega[2,6,t] <- w[1,t]
    
    # "PAIRS"
    omega[3,1,t] <- 0
    omega[3,2,t] <- 0
    omega[3,3,t] <- phi[2,t]
    omega[3,4,t] <- h[2,t]
    omega[3,5,t] <- rw[2,t]
    omega[3,6,t] <- w[2,t]
    
    # "NEWLY DEAD LEGAL HUNTING"[RB]
    omega[4,1,t] <- 0
    omega[4,2,t] <- 0
    omega[4,3,t] <- 0
    omega[4,4,t] <- 0
    omega[4,5,t] <- 0
    omega[4,6,t] <- 1
    
    # "NEWLY DEAD OTHER SOURCES AND DEAD"[RB]
    omega[5,1,t] <- 0
    omega[5,2,t] <- 0
    omega[5,3,t] <- 0
    omega[5,4,t] <- 0
    omega[5,5,t] <- 0
    omega[5,6,t] <- 1
    
    omega[6,1,t] <- 0
    omega[6,2,t] <- 0
    omega[6,3,t] <- 0
    omega[6,4,t] <- 0
    omega[6,5,t] <- 0
    omega[6,6,t] <- 1
    
    for(i in 1:n.individuals){ 
      z[i,t+1] ~ dcat(omega[z[i,t],1:6,t]) 
    }#i 								
    
  }#t 
  
  
  ##---------------------------------------------------------------------------------------------   
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  # PRIORS
  for(t in 1:n.years){
    betaResponse[t] ~ dunif(-5,5)
    betaResponseOth[t] ~ dunif(-5,5)
    
    for(st in 1:2){
      sigma[st,t] ~ dunif(0, 50)
    }
    
    for(n in 1:nTrapCovs){
      trapBetas[n,t] ~ dunif(-5,5)
    }
    
    for(n in 1:nTrapCovsOth){
      trapBetasOth[n,t] ~ dunif(-5,5)
    }
  }
  
  for(c in 1:n.counties){
    for(t in 1:n.years){
      p0[c,1,t] ~ dunif(0,1)
      p0[c,2,t] ~ dunif(0,1)
    }#t
  }#c     
  
  pResponse ~ dunif(0, 1)
  
  for(i in 1:n.individuals){ 
    idResponse[i,1] ~ dbern(pResponse)
  }
  
  
  for(c in 1:n.countries){
    for(t in 1:n.years){
      p0Oth[c,1,t] ~ dunif(0,1)
      p0Oth[c,2,t] ~ dunif(0,1)
    }#t
  }#
  
  
  for(t in 1:n.years){
    for(i in 1:n.individuals){
      # STRUCTURED 
      y.alive[i,1:nMaxDetectors,t] ~ dbinomLocal_normalWolf(detNums = nbDetections[i,t],
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
                                                            trapCovsIntercept =  detCounties[1:n.detectors],
                                                            indCov = idResponse[i,t],
                                                            indBeta = betaResponse[t],
                                                            trapCovs =  trapCovs[1:n.detectors,1:nTrapCovs,t],
                                                            trapBetas = trapBetas[1:nTrapCovs,t],
                                                            lengthYCombined = 1
      )
      
      ##OTHERS
      y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbinomLocal_normalWolf(detNums = nbDetectionsOth[i,t],
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
                                                                  trapCovsIntercept =  detCountries[1:n.detectors],
                                                                  indCov = idResponse[i,t],
                                                                  indBeta = betaResponseOth[t],
                                                                  trapCovs =  trapCovsOth[1:n.detectors,1:nTrapCovsOth,t],
                                                                  trapBetas = trapBetasOth[1:nTrapCovsOth,t],
                                                                  lengthYCombined = 1
      )
      
      x.deadculled[i,t] ~ dbern(z[i,t]==4) 
      x.deadOther[i,t] ~ dbern(z[i,t]==5) 
      
    }#i
  }#t
  
  
  
  ##---------------------------------------------------------------------------------------------										  
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(i in 1:n.individuals){ 
    isAlive[i,1] <- (z[i,1] == 2) + (z[i,1] == 3)
    state[i,1] <- (z[i,1] == 3)
    for(t in 1:n.years1){
      isAlive[i,t+1] <- (z[i,t+1] == 2) + (z[i,t+1] == 3)
      state[i,t+1] <- (z[i,t+1] == 3)
    }
  }
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

### IDENTIFY INDIVIDUALS DEAD TO CULLING
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
      if(t.recovered<= dim(z)[2]){
        #if(t.recovered< dim(z)[2]){
        
        zz[(t.recovered-1)] <- 2
        print(1)
      }
      
      # if(t.recovered< dim(z)[2]){
      #   zz[(t.recovered-1)] <- 2
      # }
      
    }
  }
  return(zz)
}))


# #if dead last year, make sure it was alive the year before
# if(out[length(out)]==5 & out[length(out)-1]==1 ) {
#   out[length(out)-1] <- 2
# }


z[z%in%4] <- 6
z[50,]

## id culled get the state 4. id not culled get the 5. 
z <- ifelse(z==3 & other.mx %in% c(1), 5 ,z)#remove other causes of mortality
z <- ifelse(z==3 & legal.mx %in% c(1), 4 ,z)

### IDENTIFY SOCIAL STATE
z <- ifelse(z==2 & indSocialState %in% c(1), 2 ,z)
z <- ifelse(z==2 & indSocialState %in% c(2), 3 ,z)

# CREATE INITIAL Z VALUES  
z.init <- t(apply(z, 1, function(zz){
  out <- zz
  out[] <- 1
  if(any(!is.na(zz))){
    
    # if(sum(zz == 5, na.rm = T)<1){
    #   range.det <- range(which(!is.na(zz)))
    #   if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
    #   if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 6
    # }
    # 
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
x.deadculled[] <- ifelse(z.age%in%c(4) & legal.mx==1,1,0)
x.deadculled <- t(apply(x.deadculled, 1, function(x){
  out <- x
  out[] <- 0
  if(any(x==1)) out[min(which(x==1))] <- 1
  return(out)
}))

x.deadOther[] <- ifelse(z.age%in%c(5) & other.mx==1,1,0)
x.deadOther <- t(apply(x.deadOther, 1, function(x){
  out <- x
  out[] <- 0
  if(any(x==1)) out[min(which(x==1))] <- 1
  return(out)
}))


#z.init["UI417471 G155-22 V1047 +","2021"] <- 2
# z.init["UI418546 G172-22 +",9] <- 2

## ------   3. GENERATE sxy & sxy.init ARAYS  ------ 


## ------   4. NIMBLE DATA  ------ 
## ------     4.1. TRAP COVARIATES  ------ 
#STRUCTURED 
trapCovs <- array(0,c(nrow(detTracks), 2, nYears))
trapCovs[,1,] <- detTracks
trapCovs[,2,] <- detSnow
#OTHERS 
trapCovsOth <- array(0,c(nrow(detTracks), 3, nYears))
trapCovsOth[,1,] <- detRoads
trapCovsOth[,2,] <- detSnow
trapCovsOth[,3,] <- detOtherSamples

## ------     4.2. RESCALE COORDINATES  ------ 
detectorsxy <- st_coordinates(myDetectors$main.detector.sp)
habitatxy <- st_coordinates(myHabitat.list$habitat.sp[myHabitat.list$habitat.r[]==1,])
colnames(habitatxy) <- colnames(detectorsxy) <- c("x","y")

ScaledDetectors <- scaleCoordsToHabitatGrid(coordsData = detectorsxy,
                                            coordsHabitatGridCenter = habitatxy)

ScaledLowUpCoords <- getWindowCoords(scaledHabGridCenter = ScaledDetectors$coordsHabitatGridCenterScaled,
                                     scaledObsGridCenter = ScaledDetectors$coordsDataScaled)


#SXY INITS 
#create a data.frame with all detections of all Individuals detected
#project death to the next year
myData.deadProj <- myData.dead[,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1#project dead reco to the next year
#remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year ),]
myData.deadProj$Id

AllDets <- rbind(myData.alive$myData.sp[,c("Id","Year")],
                 myData.deadProj[,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")
AllDetsxyscaled <- scaleCoordsToHabitatGrid(coordsData = AllDetsxy,
                                            coordsHabitatGridCenter = myHabitat.list$habitat.xy,
                                            scaleToGrid =T )$coordsDataScaled

AllDetections <- cbind(AllDetections, AllDetsxyscaled)

idAugmented <- which(rownames(z) %in%"Augmented")




which(!y.ar$Id.vector%in% unique(AllDets$Id))

y.ar$Id.vector[788]

myData.dead[myData.dead$Id %in%y.ar$Id.vector[788], ]
y.dead[y.ar$Id.vector[788],]
sum(y.alive[y.ar$Id.vector[788],,])



sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = y.ar$Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = ScaledLowUpCoords$lowerHabCoords,
                       upperCoords = ScaledLowUpCoords$upperHabCoords,
                       habitatGrid = ScaledLowUpCoords$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal"
                       
)

sxy.init[880,,]
# sxy.init <- MakeInitsXY( y = y.alive,
#                          detector.xy = detector.xy,
#                          habitat.r = myHabitat.list$habitat.r,
#                          ydead = y.dead,
#                          detector.xyDead = detector.xy,    
#                          dist.move = 5000)
# 
# 
# habitat.poly <- aggregate(rasterToPolygons(myHabitat.list$habitat.r,fun = function(x){x==1}))
# colnames(sxy.init) <- c("x","y")

# SXY DATA 
sxy.data <- sxy.init
sxy.data[] <- NA
i=783
t=9
for(i in 1:length(y.ar$Id.vector)){
  for(t in 1:(dim(sxy.data)[3]-1)){
    if(sum(z.age[i,t+1] %in% c(4,5,6))>0 ){
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

## ------     4.3. CREATE CACHE DETECTORS OBJECTS  ------ 
DetectorIndexLESS <- getLocalObjects(habitatMask = myHabitat.list$habitat.mx,
                                     coords = ScaledDetectors$coordsDataScaled,
                                     dmax =  myVars$DETECTIONS$maxDist*1.4/res(myHabitat.list$habitat.r)[1],
                                     resizeFactor = 1,
                                     plot.check = TRUE
)


## ------     4.4. TRANSFORM Y TO SPARSE MATRICES  ------ 
#STRUCTURED 
SparseY <- getSparseY(y.aliveStructured)
#OTHER
SparseYOth <- getSparseY(y.aliveOthers)

## ------     4.5. LATENT VARIABLE DET RESPONSE ------ 
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
                      detector.xy = as.matrix(ScaledDetectors$coordsDataScaled),
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
if(myVars$DATA$sex %in% "Hunn"){
  if(!dir.exists(file.path(myVars$WD,"/Figures",myVars$modelName))){dir.create(file.path(myVars$WD,"/Figures", myVars$modelName))}
  
  save(myHabitat.list, myDetectors, COUNTRIES,
       myStudyArea.poly,COMMUNES,habitat.subdetectors,
       myFilteredData.sp, myFullData.sp, COUNTIESplot,
       detCounties.original,
       file = file.path(paste(myVars$WD,"/Figures/",myVars$modelName,sep=""), "NecessaryObjects.RData" ))
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
       file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"_INPUTChain",c,".RData", sep="" ))
}

## ------   7. NIMBLE RUN ------ 
load(file.path(myVars$WD, myVars$modelName, paste(myVars$modelName, "_INPUTChain1.RData", sep="" )))
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
    
    load(file.path( myVars$WD,
                    myVars$modelName,
                    paste0(myVars$modelName, "_INPUTChain", ch, ".RData")))
    
    
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
           myVars$WD,
           myVars$modelName,
           paste0("Snap", myVars$modelName, years[t], "_", ch, ".RData")))
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
       myVars$WD, 
       myVars$modelName,
       paste0(myVars$modelName,"_OUTPUT.RData")))



## -----------------------------------------------------------------------------

## ------ V. PROCESS RESULTS ------

## ------   1. LOAD & PROCESS OPSCR OUTPUTS ------ 

## COMBINE DIFFERENT CHAINS (IF ANY) INTO ONE 
outDirectories <- list.files(file.path(myVars$WD, myVars$modelName))[grep("NimbleOut", list.files(file.path(myVars$WD, myVars$modelName)))]
path.list <- file.path(myVars$WD, myVars$modelName, outDirectories)[c(1:4)]

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
  pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_ESTIMATES.pdf",sep="")))
  
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
  pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Snap_ESTIMATES.pdf",sep="")))
  
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