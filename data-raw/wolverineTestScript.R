##------------------------------------------------------------------------------
##
## Script name: RovQuant BEAR OPSCR analysis 
##
## Purpose of script: 
## This R script performs:
## 1. the initial cleaning of the Brown bear NGS data downloaded from RovBase.3.0
## 2. the data preparation for the RovQuant OPSCR analysis with the 'nimbleSCR' package
## 3. the model fitting using 'nimble' and 'nimbleSCR'
## 4. the post-processing of the MCMC output
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), `r paste(format(Sys.Date(), "%Y"))`
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
##------------------------------------------------------------------------------
##
## Notes: 
## This is based on 'rovquantR' beta version 0.2
##   
##------------------------------------------------------------------------------
rm(list = ls())
gc()


## ------ IMPORT REQUIRED LIBRARIES ------

devtools::install_github("PierreDupont/rovquantR")
## Ctrl + Shift + F10 (to restart R session)
library(rovquantR)
library(nimbleSCR)



##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/2024"



##------------------------------------------------------------------------------
## ----- II. CLEAN NGS DATA -----

cleanRovbaseData( 
  species = "wolverine",
  years = 2015:2024,
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  species = "wolverine",
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- IV. FIT ROVQUANT MODELS ------

## -----   1. Females ------

##-- List all prepared input files
inputFiles <- list.files(file.path( working.dir, "nimbleInFiles/Hunn"),
                         full.names = T)

##-- Load the first one
load(inputFiles[1]) 

##-- Build nimble model object
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 
model$calculate()
cmodel <- compileNimble(model)
conf <- configureMCMC( model,
                       monitors = nimParams,
                       thin = 1,
                       monitors2 = nimParams2,
                       thin2 = 5)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble( list(model = model,
                                    mcmc = Rmcmc),
                               showCompilerOutput = F)
Cmcmc <- compiledList$mcmc

##-- RUN NIMBLE MCMC IN SUCCESSIVE BITES
system.time(runMCMCbites( mcmc = Cmcmc,
                          bite.size = 100,
                          bite.number = 5,
                          path = file.path(working.dir,"nimbleOutfiles/Hunn")))



## -----   2. Males ------

##-- List all prepared input files
inputFiles <- list.files(file.path(working.dir, "nimbleInFiles/Hann"),
                         full.names = T)

##-- Load the first one
load(inputFiles[1]) 

##-- Build nimble model object
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 
model$calculate()
cmodel <- compileNimble(model)
conf <- configureMCMC(model,
                      monitors = nimParams,
                      thin = 1,
                      monitors2 = nimParams2,
                      thin2 = 5)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model,
                                   mcmc = Rmcmc),
                              showCompilerOutput = F)
Cmcmc <- compiledList$mcmc

##-- RUN NIMBLE MCMC IN SUCCESSIVE BITES
system.time(runMCMCbites( mcmc = Cmcmc,
                          bite.size = 100,
                          bite.number = 5,
                          path = file.path(working.dir,"nimbleOutfiles/Hann")))



##------------------------------------------------------------------------------
## ----- V. PROCESS ROVQUANT OUTPUT ------

processRovquantOutput(   
  ##-- paths
  species = "Wolverines",
  data.dir = data.dir,
  working.dir = working.dir)


##------------------------------------------------------------------------------

source(file.path(dir.git, "Temp/CM/functions/Nimble/dbin_LESS_Cached_MultipleCovResponse.R"))



## -----------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----

myVars <- list(
  ## WORKING DIRECTORY & MODEL NAME
  WD = file.path(dir.dropbox,"wolverine/CM/2023"),
  WD_ASP = file.path(dir.dropbox,"wolverine/ASP"),
  modelName = "53.aJ_MaCleaned2023",
  modelName_ASP = "Multistate_noYear",
  
  ## HABITAT SPECIFICATIONS
  HABITAT = list( countries = c("SWE","NOR"),
                  habResolution = 20000,
                  habBuffer = 60000),
  
  ## NGS DATA SPECIFICATIONS
  DATA = list( years = 2013:2022,
               species = c("Jerv"),            
               sex = c("Hann"),                
               samplingMonths = list(12,1:6)), 
  
  ## DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 2000,
                    detResolution = 10000),
  
  ## DATA GENERATION
  DETECTIONS = list( maxDetDist = 84000,
                     resizeFactor = 2,
                     aug.factor = 0.8),
  
  # ## [PD] PLEASE USE THE ACTUAL VALUES IN THE LIST !
  # ## (i.e. USE myVars$DETECTIONS$maxDetDist = 2.1*40000 = 84000)
  # maxDistReCalc <- 2.1*myVars$DETECTIONS$maxDetDist '
  # ## (i.e. USE myVars$DETECTIONS$resizeFactor = 2 as stated in the code not 3 as it was stated here...)
  
  ## OUTPUT PLOTS
  OUTPUT = list(mapResolution = 10000),
  
  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

if(!dir.exists(file.path(myVars$WD_ASP, myVars$modelName))){dir.create(file.path(myVars$WD_ASP, myVars$modelName))}



## -----------------------------------------------------------------------------
## ------ I. CLEAN ROVBASE DATA ------
## ------ II. CLEAN ROVBASE TRACKS ------
## ------   3. SEARCH EFFORT DATA ------

## ------     3.1.GPS SEARCH TRACKS ------

### STEPS:
### 1 - Check if there is a "tracks directory" in the "data" folder
### 2 - Load all tracks in said folder
### 3 - Combine all tracks
### 4 - Check for any duplication
### 5 - Print report


dir.exists(file.path(dir.data, "GIS/tracks"))

TRACKS_SINGLE <- st_read(file.path(dir.dropbox,
                                   "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Aktivitetslogg_20231020/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20231020_date.shp")
TRACKS_MULTI <- st_read(file.path(dir.dropbox, 
                                  "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Aktivitetslogg_20231020/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20231020_date.shp")

## COMBINE ALL TRACKS
ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI)
## REMOVE HELICOPTER TRACKS
ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helikopter=="0", ]
# KEEP ONLY WOLVERINE TRACKS
ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Jerv == "1", ]

## SELECT TRACKS YEAR
TRACKS_YEAR <- TRACKS_YEAR.sp <- list()
for(t in 1:nYears){
   ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
   TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[1]], ]
   TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[2]], ]
   TRACKS <- rbind(TRACKS_1, TRACKS_2)
   ## SIMPLIFY TRACKS SHAPES
  # TRACKS <- st_simplify(TRACKS, T, 100)
   TRACKS <- st_intersection(TRACKS, st_as_sf(myStudyArea))
   ## NAME TRACKS
   TRACKS$ID <- 1:nrow(TRACKS)
   TRACKS_YEAR[[t]] <- TRACKS
}#t






## ------ III. MAKE ROVQUANT DATA ------
## ------ IV. PROCESS ROVQUANT OUTPUT ------
## -----------------------------------------------------------------------------
## ------   1. HABITAT DATA ------

## ------     1.1.LOAD RAW SHAPEFILES ------

## POLYGONS OF THE REGION
GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) ## Map of Scandinavia (including Finland & parts of Russia)
GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))

## POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP[GLOBALMAP$ISO %in% c("SWE","NOR"), ]
COUNTRIES <- COUNTRIES %>%    group_by(ISO) %>% summarize()

## POLYGONS OF COMMUNES IN SWEDEN & NORWAY
COMMUNES_NOR <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp"))   ## Communal map of Norway
COMMUNES_SWE <- st_read(file.path(dir.dropbox, "DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp"))    ## Communal map of Sweden
COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)

## POLYGONS OF COUNTIES IN SWEDEN & NORWAY
COUNTIES <- COMMUNES %>% group_by(NAME_1) %>% summarize()

## AGGREGATE COUNTIES 
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
COUNTIES_AGGREGATED <- st_simplify( COUNTIES_AGGREGATE,
                                    preserveTopology = T,
                                    dTolerance = 500)
COUNTIES_AGGREGATED$id <- COUNTIES_AGGREGATE$id



## ------     1.2.CREATE STUDY AREA POLYGON ------

## CREATE STUDY AREA POLYGON BASED ON COUNTRY NAMES
myStudyArea <- COUNTRIES %>%
  filter(ISO %in% myVars$HABITAT$countries) %>%
  mutate(id = 1) %>%
  group_by(id) %>%
  summarize()

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myStudyArea), add = TRUE, col ="red")
}



## ------   2. NGS DATA ------

## ------     2.1.LOAD ROVBASE FILES ------

DNA <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/DNA_carnivores.csv"), fileEncoding = "latin1")## NGS data from RovBase#[CM update to 20190627]
DEAD <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/dead_carnivores.csv"), fileEncoding = "latin1") ## Dead Recoveries from RovBase#[CM update to 20190627]
SUSPECT_NGS_SAMPLES <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Remove ngs samples list wolverine 2023.csv"), fileEncoding = "latin1") ## DNA samples to be removed from Henrik
SUSPECT_DeadRecoSAMPLES <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Remove dead recoveries list wolverine 2023.csv"), fileEncoding = "latin1") ## DNA samples to be removed from Henrik
DEN <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/DEN_COUNTS_2009_2023_fromHB.csv"), fileEncoding = "latin1")
skandObs <- read_xlsx(file.path(dir.dropbox, "DATA/Skandobs/20231024 - Richard_Bischof_Skandobs_2012_2023/20231024 - Richard_Bischof_Skandobs_2012_2023.xlsx"))
rovbaseObs <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/DNA_carnivores.csv"), fileEncoding = "latin1")



## ------     2.2.TRANSLATE SCANDINAVIAN CHARACTERS ------

colnames(DNA) <- translateForeignCharacters(dat = colnames(DNA), dir.translation = dir.analysis)
colnames(DEAD) <- translateForeignCharacters(dat = colnames(DEAD), dir.translation = dir.analysis)
colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN), dir.translation = dir.analysis)
colnames(skandObs) <- translateForeignCharacters(dat = colnames(skandObs), dir.translation = dir.analysis)
colnames(rovbaseObs) <- translateForeignCharacters(dat = colnames(rovbaseObs), dir.translation = dir.analysis)
rovbaseObs$Proevetype <- translateForeignCharacters(dat = rovbaseObs$Proevetype, dir.translation = dir.analysis)



## ------   3. SEARCH EFFORT DATA ------

# ## ------     3.1.GPS SEARCH TRACKS ------
#
# TRACKS_SINGLE <- st_read(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Aktivitetslogg_20231020/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20231020_date.shp", sep = ""))
# TRACKS_MULTI <- st_read(paste(dir.dropbox, "/DATA/RovbaseData/ROVBASE DOWNLOAD 20231020/Aktivitetslogg_20231020/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20231020_date.shp", sep = ""))
# 
# ## COMBINE ALL TRACKS
# ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI)
# ## REMOVE HELICOPTER TRACKS
# ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helikopter=="0", ]
# # KEEP ONLY WOLVERINE TRACKS
# ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Jerv == "1", ]
# 
# ## SELECT TRACKS YEAR
# TRACKS_YEAR <- TRACKS_YEAR.sp <- list()
# for(t in 1:nYears){
#    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
#    TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[1]], ]
#    TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[2]], ]
#    TRACKS <- rbind(TRACKS_1, TRACKS_2)
#    ## SIMPLIFY TRACKS SHAPES
#   # TRACKS <- st_simplify(TRACKS, T, 100)
#    TRACKS <- st_intersection(TRACKS, st_as_sf(myStudyArea))
#    ## NAME TRACKS
#    TRACKS$ID <- 1:nrow(TRACKS)
#    TRACKS_YEAR[[t]] <- TRACKS
# }#t
# 
#
#
# ## ------     3.2.DISTANCE TO ROADS ------
#
# ## LOAD MAP OF DISTANCES TO ROADS (1km resolution)
# DistAllRoads <- raster(paste(dir.dropbox,"/DATA/GISData/Roads/MinDistAllRoads1km.tif", sep=""))
# r   <- fasterize(st_as_sf(myStudyArea), DistAllRoads)
# r[!is.na(r)] <- DistAllRoads[!is.na(r)]
# DistAllRoads <- r
# DistAllRoads <- crop(DistAllRoads, myStudyArea)
# ## PLOT CHECK
# if(myVars$plot.check){
#    plot(DistAllRoads)
#    plot(myStudyArea,add=T)
# }
# 
#
#
# ## ------     3.3.DAYS OF SNOW ------
#
# ## SEASONAL MAPS (CREATED IN TEMP/CM/GIS/snowMODIS)
# SNOW <- stack(paste(dir.dropbox,"/DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2008_2023_09.tif", sep=""))
# ## RENAME THE LAYERS
# names(SNOW) <- paste(2008:2022,(2008:2022)+1, sep="_")
# ## SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
# SNOW <- SNOW[[paste("X", years, "_", years+1, sep="")]]
# SNOW <- raster::crop(SNOW, c(0,40,55,75))
#
#
#
# ## ------     3.4.SAVE SEARCH EFFORT OBJECTS FOR FASTER RUNS ------
#
# save(TRACKS_YEAR, SNOW, DistAllRoads, file = file.path(myVars$WD, "TRACKS0mSouthSweden20132023Cleaned.RData"))
load(file.path(myVars$WD, "TRACKS0mSouthSweden20132023Cleaned.RData"))



## ------   4. LOAD SCANDINAVIAN 20KM HABITAT ------

load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/Habitat20km.RData"))
load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData"))



## ------   1. CLEAN & FILTER NGS DATA ------

## ------     1.1. CLEAN NGS & DEAD RECOVERY DATA ------

## Remove DEAD entries from the DNA data [HB]
DNA <- DNA[substr(DNA$RovbaseID..Proeve.,1,1) != "M", ] 

## Remove un-verified dead recoveries [HB]
## ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
DEAD <- DEAD[!grepl(pattern = "Påskutt", x = as.character(DEAD$Utfall)), ]

## CLEAN DATA
myCleanedData.sp <- CleanDataNew2sf( 
  dna_samples = DNA,
  dead_recoveries = DEAD,
  species_id = myVars$DATA$species,
  country_polygon = COUNTRIES,
  threshold_month = unlist(myVars$DATA$samplingMonths)[1],
  keep_dead = T,
  age.label.lookup = age.lookup.table)

## PLOT CHECK
if(myVars$plot.check){
  plot(st_geometry(COUNTRIES))
  plot(st_geometry(myStudyArea), add = T, col ="red")
  plot(st_geometry(myCleanedData.sp), add = TRUE, pch = 19, cex = 0.2, col = "purple")
} 



## ------     1.2. CHECK SEX ASSIGNMENT ------

myFullData.sp <- FilterDatasf( 
  myData = myCleanedData.sp,
  poly = myStudyArea,
  dead.recovery = T,
  sex = c("Hann","Hunn"),
  setSex = T)

if(myVars$plot.check){
  plot( st_geometry(myFullData.sp$alive), add = TRUE, pch = 19, cex = 0.2, col = "lightblue")
}



## ------     1.3. FILTER OUT SUSPECT SAMPLES ACCORDING TO HENRIK ------

myFullData.sp$alive$DNAID <- as.character(myFullData.sp$alive$DNAID)
myFullData.sp$dead.recovery$DNAID <- as.character(myFullData.sp$dead.recovery$DNAID)
myFullData.sp$alive <- myFullData.sp$alive[!(myFullData.sp$alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[!(myFullData.sp$dead.recovery$RovBaseId %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]

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

## Remove pups killed before recruitment based on weight (cf. Henrik)
## 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and 
## recovered dead between March and November
sum(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
      myFullData.sp$dead.recovery$Month > 2 &
      myFullData.sp$dead.recovery$Month < 12)

myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$Alder.pa.doedt.individ %in% "Unge" &
                                                                    myFullData.sp$dead.recovery$Month > 2 &
                                                                    myFullData.sp$dead.recovery$Month < 12),]

## 2) remove individuals that have a weight >0 and <4 between March and November
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

## Assign negative values to nas to avoid issues
myFullData.sp$dead.recovery$weight[is.na(myFullData.sp$dead.recovery$weight)] <- -999

# WEIGHT DISTRIBUTION
# par(mfrow=c(4,3))
# for(t in 1:12){
# hist(myFullData.sp$dead.recovery$weight[(myFullData.sp$dead.recovery$weight >-1) &
#                                           myFullData.sp$dead.recovery$Month%in% t],breaks=c(0:30), main=t,xlab="Weight")
# }
# # AGE DISTRIBUTION
# par(mfrow=c(4,3))
# for(t in 1:12){
#     hist(myFullData.sp$dead.recovery$Age[(myFullData.sp$dead.recovery$Age >-1) &
#                                          myFullData.sp$dead.recovery$Month%in% t],breaks=seq(-0.01,20.99,by=1),
#          main=t,xlab="Age")
# }

## Check how many dead recoveries we remove and remove if more than 0
if(sum(myFullData.sp$dead.recovery$weight > 0 &
       myFullData.sp$dead.recovery$weight < 4 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2) > 0){
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[-which(myFullData.sp$dead.recovery$weight > 0 &
                                                                      myFullData.sp$dead.recovery$weight < 4 &
                                                                      myFullData.sp$dead.recovery$Month < 12 &
                                                                      myFullData.sp$dead.recovery$Month > 2), ]
}

## Check how many dead reco with a weight of 0 kg and recovered between march and november
if(sum(myFullData.sp$dead.recovery$Age %in% 0 &
       myFullData.sp$dead.recovery$Month < 12 &
       myFullData.sp$dead.recovery$Month > 2) > 0){
  myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Age %in% 0 &
                                myFullData.sp$dead.recovery$Month < 12 &
                                myFullData.sp$dead.recovery$Month > 2, ]
}



## ------     1.4. FILTER NGS & DEAD RECOVERY DATA FOR DATES ------

myFilteredData.sp <- myFullData.sp

## Subset to years of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years, ]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year %in% years, ]

## Subset to months of interest
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Month %in% unlist(myVars$DATA$samplingMonths), ]

## [PD] ADD CHECKS OF THE NUMBER OF SAMPLES/INDIVIDUALS REMOVED FOR rovquantR REPORT



## ------     1.5. FILTER OUT NGS IN NORRBOTTEN EXCEPT IN 2017, 2018 and 2019 ------ 

### [CM] We subset data in Norrbotten in specific years because there was no 
### comprehensive sampling (few detections that we can't use), so we remove 
### detections in all years except 2017/18- 2018/2019 --> 2016:2018 
### (for the wolverine we call those winters 2017-2019 because there are many 
### detections collected in December. mostly January-April)

COUNTIESNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten", ]
yearsSampledNorrb <- c(2016:2018)

## Identify DNA samples collected in Norrbotten
is.Norr <- as.numeric(st_intersects(myFilteredData.sp$alive, COUNTIESNorrbotten))

## Check how many detections are removed.
table(myFilteredData.sp$alive[which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                      !is.na(is.Norr)), ]$Year)

## Remove DNA samples collected inside Norrbotten in those years
myFilteredData.sp$alive <- myFilteredData.sp$alive[- which(!myFilteredData.sp$alive$Year %in% yearsSampledNorrb &
                                                             !is.na(is.Norr)), ]
## plot check
for(t in 1:nYears){
  plot(st_geometry(myStudyArea))
  plot(st_geometry(COUNTIESNorrbotten), add = T, col = "blue")
  plot(st_geometry(myFilteredData.sp$alive[myFilteredData.sp$alive$Year %in% years[t], ]),
       col = "red", add = T, pch = 16)
}#t



## ------     1.6. FILTER NGS & DEAD RECOVERY DATA FOR SEX ------

## myFullData.sp
myFullData.sp$alive <- myFullData.sp$alive[myFullData.sp$alive$Sex %in% myVars$DATA$sex, ]
myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Sex %in% myVars$DATA$sex, ]

## myFilteredData
myFilteredData.spAllSex <- myFilteredData.sp
myFilteredData.sp$alive <- myFilteredData.sp$alive[myFilteredData.sp$alive$Sex %in% myVars$DATA$sex, ]
myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Sex %in% myVars$DATA$sex, ]

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,3))
  for(t in 1:nYears){
    ## DEAD RECOVERIES TOTAL
    tempTotal <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t] &
                                                   myFilteredData.sp$dead.recovery$Sex %in% myVars$DATA$sex, ]
    NGS_TabTotal <- table(tempTotal$Country)
    ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country), 2, function(x)sum(x > 0))
    ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD 
    ## PLOT NGS SAMPLES
    plot(st_geometry(GLOBALMAP), col = "gray80")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    #plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
    plot(st_geometry(tempTotal), pch = 21, bg = "darkred", add = T)
    ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    graphics::text( x = 100000, y = 7250000,
                    labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"),
                    cex = 1.1, col = "firebrick3", font = 2)
    graphics::text( x = 820000, y = 6820000,
                    labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"),
                    cex = 1.1, col = "navyblue", font = 2)
    ## ADD OVERALL NUMBERS
    mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
  }#t
}



## ------     1.7. ASSIGN SAMPLES TO GPS TRACKS ------

### It is not possible to distinguish between structured and unstructured samples.
### We assign each sample to structured or unstructured based on whether a given
### sample matches in time and space with recorded search tracks:
### - a sample is assigned to "structured" if it was collected by the authorities
### ("Statsforvalteren" or "SNO").
### -  a sample is assigned to "structured" if it was located within 500 m from
### a GPS search track recorded the same day.

## ASSIGN GPS TRACK ROVBASE ID
myFilteredData.sp$alive$TrackRovbsID <- NA
myFilteredData.sp$alive$TrackDist <- NA

for(i in 1:nrow(myFilteredData.sp$alive)){
  ## SUBSET TO THIS DETECTION
  thisDet <- myFilteredData.sp$alive[i, ]
  t <- which(years == thisDet$Year)
  
  ## IDENTIFY TRACKS WITH THE SAME DATE AS THE DETECTION
  whichSameDate <- which(as.character(TRACKS_YEAR[[t]]$Dato) == as.character(thisDet$Date))
  ## IF NO TRACKS WITH THE SAME DATE ==> MOVE TO NEXT SAMPLE
  if(length(whichSameDate) == 0){next}
  
  ## CHECK IF SOME TRACKS ARE WITHIN 750m OF THE DETECTION
  theseTracks <- TRACKS_YEAR[[t]][whichSameDate, ]
  tmpTRACKS <- st_intersection( x = theseTracks,
                                y = st_buffer( x = thisDet, dist = 750))
  if(nrow(tmpTRACKS)==0){next}
  dist <- st_distance(x = thisDet, y = tmpTRACKS, by_element = F)
  myFilteredData.sp$alive$TrackRovbsID[i] <- as.character(tmpTRACKS$RovbaseID[which.min(dist)])
  myFilteredData.sp$alive$TrackDist[i] <- min(dist)
}#i
# ## SAVE THE FOR FASTER LOADING
# save( myFilteredData.sp,
#       file = file.path(myVars$WD, myVars$modelName,"_myFilteredData.sp.RData"))
# load(file.path(myVars$WD, myVars$modelName,"_myFilteredData.sp.RData"))



## ------     1.8. SEPARATE MORTALITY CAUSES ------ 

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
    #plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
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
  
  
  ## PLOT TREND DETECTIONS AND DEAD RECOVERIES OVER TIME AND SPACE #[CM]
  #DETECTIONS
  #pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_TRENDDetections.pdf",sep="")))
  temp <- unique(myFilteredData.sp$alive[,c("Year","Country","DNAID")])
  tab_Country.Year <- table(temp$Year, temp$Country)
  country.colors <- c("goldenrod1","goldenrod3")
  par(mfrow=c(1,1), mar=c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
  lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  #ID DETECTED
  temp <- table(myFilteredData.sp$alive$Year,myFilteredData.sp$alive$Country,myFilteredData.sp$alive$Id)
  tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
  country.colors <- c("goldenrod1","goldenrod3")
  par(mfrow=c(1,1), mar=c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
  lines(tab_Country.Year1[,"N"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year1[,"S"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  ## Average number of detection per detected ID  #[CM]
  tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
  par(mfrow=c(1,1), mar=c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year2))), ylim=c(0,max(tab_Country.Year2)),
       ylab="Average Number of detections", xlab="Years")
  lines(tab_Country.Year2[,"N"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year2[,"S"]~as.numeric(row.names(tab_Country.Year2)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("bottomright",c("N","S"), fill=country.colors)
  
  ## deadrecovery #[CM]
  temp <- unique(myFilteredData.sp$dead.recovery[,c("Year","Country","Id")])
  tab_Country.Year <- table(temp$Year, temp$Country)
  country.colors <- c("goldenrod1","goldenrod3")
  par(mfrow=c(1,1), mar=c(5,5,5,5))
  plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)),
       ylab="N Id Dead recovered", xlab="Years")
  lines(tab_Country.Year[,"N"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
  lines(tab_Country.Year[,"S"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
  legend("topright",c("N","S"), fill=country.colors)
  #dev.off()
}




## ------   2. GENERATE HABITAT ------

## ------     2.1. REDUCE THE AREA OF THE STATE-SPACE BASED ON DETECTIONS ------

## DELINEATE A BUFFER AROUND ALL DEAD RECOVERIES 
BuffDead <-  myFilteredData.spAllSex$dead.recovery %>%
  st_buffer(.,dist = myVars$HABITAT$habBuffer) %>%
  summarize() 

## DELINEATE A BUFFER AROUND ALL DETECTIONS
### This is the buffer used to delineate the study area (dead buffer is not 
### included because there were some dead recoveries found very far)
myBufferedArea <- myFilteredData.spAllSex$alive %>%
  st_buffer(.,dist = myVars$HABITAT$habBuffer * 1.4) %>%
  summarize()  

## CUT TO SWEDISH AND NORWEGIAN BORDERS
myStudyArea <- st_intersection(myBufferedArea, myStudyArea) %>%
  summarize()

## PLOT CHECK
if(myVars$plot.check){
  par(mar = c(0,0,0,0))
  plot(st_geometry(myStudyArea))
  plot( st_geometry(myFilteredData.spAllSex$alive),
        pch = 21, bg = "red", cex = 0.5, add = T)
  plot( st_geometry(myFilteredData.spAllSex$dead.recovery),
        pch = 21, bg = "blue", cex = 0.5, add = T)
  plot(st_geometry(BuffDead), add = T, border = "blue")
  plot(st_geometry(myBufferedArea), border = "red", add = T)
  plot(st_geometry(myStudyArea), border = "grey", add = T)
}



## ------     2.2. GENERATE HABITAT CHARACTERISTICS ------

## MAKE HABITAT FROM RASTER
myHabitat <- MakeHabitatFromRastersf( 
  poly = myStudyArea,
  habitat.r = habitatRasters[["Habitat"]],
  buffer = myVars$HABITAT$habBuffer,                               
  plot.check = T)

## RESCALE HABITAT COORDINATES
scaledHabGridCenters <- scaleCoordsToHabitatGrid(
  coordsData = myHabitat$habitat.xy,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid = F)$coordsHabitatGridCenterScaled
scaledHabGridCenters <- scaledHabGridCenters[myHabitat$habitat.r[] == 1, ]

## CREATE HABITAT GRID 
habIDCells.mx <- myHabitat$IDCells.mx 
habIDCells.mx[] <- 0
for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,2))
  plot(myHabitat$habitat.r, legend = F)
  plot(st_geometry(myStudyArea), add = T)
  
  image(habIDCells.mx)
}



## ------     2.3. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 

## Remove samples outside the STUDY AREA #[CM]
# myStudyArea$idd <- 1
# myStudyAreaAggregated <- myStudyArea %>% group_by(idd) %>% summarize()
whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyArea))))
if(length(whichOut)>0){
  myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
}
myFilteredData.sp$alive$Id <- droplevels(myFilteredData.sp$alive$Id)

## REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT #[CM] 
whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, myHabitat$buffered.habitat.poly))))
if(length(whichOutBuff)>0){
  myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]
}

## PLOT CHECK
if(myVars$plot.check){
  # par(mfrow = c(1,2))#[CM]
  plot( myHabitat$habitat.r)
  plot( st_geometry(myStudyArea), add = T,
        col = rgb(150/250,150/250,150/250, alpha = 0.75))
  plot( st_geometry(GLOBALMAP), add = T)
  plot( st_geometry(myHabitat$buffered.habitat.poly), add = T)
  plot( st_geometry(myFilteredData.sp$alive),
        pch = 21, bg = "red", cex = 0.5, add = T)
  plot( st_geometry(myFilteredData.sp$dead.recovery),
        pch = 21, bg = "blue", cex=0.5,add=T)
}

## Check correlation number of detections ~ between monitoring season
myFilteredData.sp$dead.recovery$Id <- as.character(myFilteredData.sp$dead.recovery$Id)
myFilteredData.sp$alive$Id <- as.character(myFilteredData.sp$alive$Id)
deadID <- unique(myFilteredData.sp$dead.recovery$Id)
ndet <- NULL
timeDiff <- NULL
for(i in 1:length(deadID)){
  tmpYear <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Year
  
  timeDiff[i] <- myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Id %in% deadID[i],]$Date-
    as.POSIXct(strptime(paste("01-12", tmpYear, sep = "-"), "%d-%m-%Y")) 
  
  ndet[i] <- length(myFilteredData.sp$alive[myFilteredData.sp$alive$Id %in% deadID[i] & 
                                              myFilteredData.sp$alive$Year %in% tmpYear, ])
}#i

## PLOT CHECK
if(myVars$plot.check){
  #pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Prop id deteced_Time available.pdf",sep="")))
  plot( ndet ~ timeDiff,
        ylab = "Total number of detections",
        xlab = "Number of days between dec 1 and dead recovery")
  hh <- hist(timeDiff[ndet > 0], breaks = seq(0, 400, by = 25))
  hh1 <- hist(timeDiff[ndet == 0], breaks = seq(0, 400, by = 25))
  barplot( rbind(hh$counts/(hh$counts+hh1$counts),
                 hh1$counts/(hh$counts+hh1$counts)),
           names.arg = hh$breaks[1:(length(hh$breaks)-1)],
           xlab = "number of days between dead reco and start monitoring",
           ylab = "%")
  legend( "topright", fill = c(grey(0.2),grey(0.8)),
          legend = c("detected","notDetected"))
  #dev.off()
}



## ------     2.4. GENERATE HABITAT-LEVEL COVARIATES ------

## ------       2.4.1. DEN COUNTS ------

DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
DEN.sp$id  <- rep(1,nrow(DEN.sp))
DEN.sp <- DEN.sp[ ,"id"]

DEN.r <- raster(
  estUDm2spixdf(
    kernelUD(as(DEN.sp,"Spatial"), h = 30000,
             grid = as(myHabitat$habitat.r, 'SpatialPixels'))))

if(myVars$plot.check){
  plot(DEN.r)
  plot(st_geometry(myStudyArea), add = TRUE, border = "black")
}

## EXTRACT COVARIATES
denCounts <- DEN.r[myHabitat$habitat.r[ ] == 1]
denCounts <- round(scale(denCounts), digits = 2)



## ------     2.5. RESCALE HABITAT COORDINATES ------ 

### [PD] NO NEED TO RERUN THE FUNCTION TWICE HERE
### COULD SIMPLY USE scaledHabGridCenters - 0.5 and scaledHabGridCenters + 0.5
# ScaledLowerCoords <- scaleCoordsToHabitatGrid(
#   coordsData = lowerHabCoords,
#   coordsHabitatGridCenter = myHabitat$habitat.xy,
#   scaleToGrid = T)$coordsDataScaled
# ScaledLowerCoords[ ,2] <- ScaledLowerCoords[ ,2] - 1
# ScaledUpperCoords <- scaleCoordsToHabitatGrid(
#   coordsData = upperHabCoords,
#   coordsHabitatGridCenter = myHabitat$habitat.xy,
#   scaleToGrid = T)$coordsDataScaled
# ScaledUpperCoords[ ,2] <- ScaledUpperCoords[ ,2] + 1

ScaledLowerCoords <- scaledHabGridCenters - 0.5
ScaledUpperCoords <- scaledHabGridCenters + 0.5



## ------   3. GENERATE DETECTORS ------
## ------     3.1. GENERATE DETECTORS CHARACTERISTICS ------

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
habitat.subdetectors <- disaggregate(
  myHabitat$habitat.rWthBuffer,
  fact = res(myHabitat$habitat.r)[1]/myVars$DETECTORS$detSubResolution)

## GENERATE NGS DETECTORS BASED ON THE STUDY AREA
myDetectors <- MakeSearchGridsf( 
  data = habitat.subdetectors,
  resolution = myVars$DETECTORS$detResolution,
  div = (myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)^2,
  plot = FALSE,
  fasterize = TRUE)

## FORMAT DETECTOR LOCATIONS & NUMBER OF TRIALS PER DETECTOR IN ARRAYS/MATRICES
detector.xy <- st_coordinates(myDetectors$main.detector.sp)
colnames(detector.xy) <- c("x","y")
n.trials <- as.vector(table(myDetectors$detector.sp$main.cell.id))

## RETRIEVE DETECTION WINDOWS BOUNDARIES
lowerDetCoords <- detector.xy - 0.5 * myVars$DETECTORS$detResolution
upperDetCoords <- detector.xy + 0.5 * myVars$DETECTORS$detResolution

## EXTRACT NUMBERS OF DETECTORS
n.detectors <- dim(myDetectors$main.detector.sp)[1]


# ## IDENTIFY DETECTORS IN NORBOTTEN 
# COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten","Troms","Västerbotten",
#                                                             "Nordland","Finnmark"), ]
# COUNTIESAroundNorrbotten <- st_simplify(COUNTIESAroundNorrbotten, dTolerance = 500)
# 
# ## CREATE A NORROBOTTEN DETECTOR GRID
# distDestsCounties <- st_distance( myDetectors$main.detector.sp,
#                                   COUNTIESAroundNorrbotten,
#                                   byid = TRUE)
# detsNorrbotten <- which(apply(distDestsCounties, 1, which.min)==3)
#
# ## PLOT CHECK 
# plot(st_geometry(COUNTIESAroundNorrbotten))
# plot(st_geometry(myDetectors$main.detector.sp), col="black",pch=16,cex=0.3,add=T)
# plot(st_geometry(myDetectors$main.detector.sp[detsNorrbotten,]), col="red",pch=16,cex=0.3,add=T)
#
# ## PLOT CHECK
# if(myVars$plot.check){
#   par(mfrow = c(1,2))
#   ## PLOT NGS DETECTORS
#   plot(st_geometry(myHabitat$buffered.habitat.poly), main = paste(n.detectors, "Detectors Alive"), col = rgb(0.16,0.67,0.16, alpha = 0.3))  
#   plot(st_geometry(myStudyArea), add = TRUE, col = rgb(0.16,0.67,0.16,alpha = 0.5))
#   plot(st_geometry(myDetectors$main.detector.sp), col = "red", pch = 16, cex = 0.1, add = TRUE)
#   plot(st_geometry(COUNTRIES), add = TRUE)
#   ## PLOT DEAD DETECTORS
#   plot(st_geometry(myHabitat$buffered.habitat.poly), main = paste(n.detectors.dead, "Detectors Dead"), col = rgb(0.16,0.67,0.16, alpha = 0.3)) 
#   plot(st_geometry(myStudyArea), add = T, col = rgb(0.16,0.67,0.16,alpha = 0.5))
#   plot(st_geometry(myDetectors.dead$main.detector.sp), col = "red", pch = 16, cex = 0.1, add = TRUE)
#   plot(st_geometry(COUNTRIES), add = TRUE)
# }



## ------     3.2. GENERATE DETECTOR-LEVEL COVARIATES ------

## ------       3.2.1. EXTRACT COUNTIES ------

## ASSIGN COUNTIES TO DETECTORS BASED ON DISTANCE
dist <- st_distance(myDetectors$main.detector.sp, COUNTIES_AGGREGATED, by_element = F)
detCounties <- apply(dist, 1, function(x) which.min(x))
COUNTIES_AGGREGATEDSubset <- COUNTIES_AGGREGATED[unique(detCounties), ]
COUNTIES_AGGREGATEDSubset$idunique <- as.numeric(as.factor(unique(detCounties)))
detCounties <- as.numeric(as.factor(detCounties))

## PLOT CHECK 
if(myVars$plot.check){
  myCol <- terrain.colors(nrow(COUNTIES_AGGREGATED))
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(st_geometry(myDetectors$main.detector.sp[detCounties%in% c(5),]), col = myCol[detCounties], pch = 16, cex = 0.8,add=T)
  plot(st_geometry(myDetectors$main.detector.sp), col = myCol[detCounties], pch = 16, cex = 0.8, add=T)
  plot(st_geometry(COUNTIES_AGGREGATED), add = TRUE)
  plot(st_geometry(myDetectors$main.detector.sp[detCounties %in%3 ,]), col = "red", pch = 16, cex = 0.8, add=T)
}



## ------       3.2.2. EXTRACT COUNTRIES ------

## ASSIGN COUNTRIES TO DETECTORS BASED ON DISTANCE
dist <- st_distance(myDetectors$main.detector.sp, COUNTRIES, by_element = F)
detCountries <- apply(dist, 1, function(x)which.min(x))

## ADD ANOTHER CATEGORY to detcountry if in Norrbotten, to turnoff detection to 0 there. 
detCountriesNorb <- matrix(NA, nrow = length(detCountries), ncol = nYears)
detCountries1 <- detCountries
detCountries1[detCounties %in% 1] <- 3
for(t in 1:nYears){
  if(t %in% which(!years %in% yearsSampledNorrb)){
    detCountriesNorb[ ,t] <- detCountries1
  } else {
    detCountriesNorb[ ,t] <- detCountries
  }
}#t 

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



## ------       3.2.3. EXTRACT GPS TRACKS LENGTHS ------

## INITIALIZE MATRIX OF GPS TRACKS LENGTH FOR EACH DETECTOR & YEAR
detectorGrid.r <- rasterFromXYZ(cbind(st_coordinates(myDetectors$main.detector.sp),
                                      rep(1,nrow(myDetectors$main.detector.sp))))
detectorGrid <- sf::st_as_sf( stars::st_as_stars(detectorGrid.r), 
                              as_points = FALSE, 
                              merge = F)
st_crs(detectorGrid) <- st_crs(myStudyArea)
detectorGrid$id <- 1:nrow(detectorGrid)

## CALCULATE THE LENGTH OF THE TRACKS
detTracks <- matrix(0, nrow = n.detectors, ncol = nYears)
TRACKS.r <- list()
for(t in 1:nYears){
  intersection <- st_intersection(detectorGrid, TRACKS_YEAR[[t]]) %>%
    mutate(LEN = st_length(.)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(transect_L = sum(LEN))## Get total length searched in each detector grid cell
  detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
  TRACKS.r[[t]] <- detectorGrid.r
  TRACKS.r[[t]][detectorGrid.r[] %in% 1] <- detTracks[ ,t]
}#t

## PLOT CHECK 
if(myVars$plot.check){
  
  max <- max(unlist(lapply(TRACKS.r, function(x) max(x[],na.rm=T))))
  cuts <- seq(0,max,length.out = 100)   #set breaks
  col <- rev(terrain.colors(100))
  CountriesDetRes <- disaggregate(habitatRasters$Countries,fact=2)
  CountriesDetRes <- crop(CountriesDetRes,TRACKS.r[[1]])
  rr <- TRACKS.r[[1]]
  rr[CountriesDetRes[]%in% 2] <- 1
  
  # pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"Tracks.pdf",sep="")))
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
  # dev.off()
}



## ------       3.2.4. EXTRACT DISTANCES TO ROADS ------

## AGGREGATE TO MATCH THE DETECTORS RESOLUTION
DistAllRoads <- aggregate( DistAllRoads,
                           fact = myVars$DETECTORS$detResolution/res(DistAllRoads),
                           fun = mean)

## EXTRACT ROAD DISTANCE FOR EACH DETECTOR
detRoads <- raster::extract(DistAllRoads, myDetectors$main.detector.sp)

## If NA returns the average value of the cells within 15000m 
isna <- which(is.na(detRoads))
tmp <- raster::extract( DistAllRoads,
                        myDetectors$main.detector.sp[isna, ],
                        buffer = 15000,
                        fun = mean,
                        na.rm = T)
detRoads[isna] <- tmp

## PLOT CHECK 
if(myVars$plot.check){
  par(mfrow = c(1,1))
  plot(st_geometry(GLOBALMAP), col = "gray80", main = "Distance to roads")
  plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
  plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
  plot(DistAllRoads,add=T)
  plot(st_geometry(myDetectors$main.detector.sp), cex=DoScale(detRoads), pch = 16, add = T)
}



## ------       3.2.5. EXTRACT DAYS OF SNOW ------

## EXTRACT SNOW 
detSnow <- matrix(0, nrow = dim(myDetectors$main.detector.sp)[1], ncol = nYears)
det.sptransf <- st_transform(myDetectors$main.detector.sp, st_crs(SNOW))
detSnow[ ,1:nYears] <- raster::extract(SNOW, det.sptransf)

## if NA returns the average value of the cells within 15000m 
isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                        buffer = 15000,
                        fun = mean, na.rm = T)
detSnow[isna,1:nYears] <- tmp

## PLOT CHECK
if(myVars$plot.check){
  plot( st_geometry(myDetectors$main.detector.sp),
        cex = DoScale(detSnow[,6],l = 0,u = 0.5),
        pch = 16)
}



## ------       3.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------
## ------          3.2.6.1. SKANDOBS ------

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
skandObs$monitoring.season <- ifelse(skandObs$month < 12, skandObs$year, skandObs$year+1) #--- need to change for other species
skandObs <- skandObs[subset, ] 

## SUBSET BASED ON SPACE 
habitat.rWthBufferPol <- sf::st_as_sf(stars::st_as_stars(myHabitat$habitat.rWthBuffer), 
                                      as_points = FALSE, merge = TRUE)
habitat.rWthBufferPol <- habitat.rWthBufferPol[habitat.rWthBufferPol$Habitat %in%1, ]
subsetSpace <- !is.na(as.numeric(st_intersects(skandObs, habitat.rWthBufferPol)))
skandObs <- skandObs[subsetSpace, ] 

## RASTERIZE AT THE DETECTOR LEVEL
r.detector <- aggregate( habitat.subdetectors,
                         fact = myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)
r.list <- lapply(years, function(y){
  rl <- raster::rasterize(skandObs[skandObs$monitoring.season %in% y, 1],
                          r.detector,
                          fun = "count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})
r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))

## PLOT CHECK 
if(myVars$plot.check){
  
  # plot(st_geometry(habitat.rWthBufferPol))
  # plot(st_geometry(skandObs),col="red",add=T)
  # plot(r.skandObsSamplesBinary[[t]])
  
  # pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"skandObs",".pdf", sep="" ), width = 10)
  ## SUMMARY SKANDOBS
  barplot(table(skandObs$monitoring.season))
  barplot(table(skandObs$month), xlab = "Months")
  barplot(table(skandObs$species))
  
  ## MAPS 
  par(mar=c(0,0,2,0))
  for(t in 1:nYears){
    plot(st_geometry(myStudyArea), main = years[t])
    plot(st_geometry(skandObs[skandObs$monitoring.season %in% years[t],  ]),
         pch = 16, col = "red", cex = 0.1, add = T)
  }#t
  # dev.off()
}



## ------          3.2.6.2. ROVBASE ------

## GET ALL SAMPLES COLLECTED
rovbaseObs <- rovbaseObs[!is.na(rovbaseObs$Oest..UTM33.SWEREF99.TM.), ]
rovbaseObs$Funnetdato <- as.POSIXct(strptime(rovbaseObs$Funnetdato,"%d.%m.%Y")) 
rovbaseObs$year <- as.numeric(format(rovbaseObs$Funnetdato,"%Y"))
rovbaseObs$month <- as.numeric(format(rovbaseObs$Funnetdato,"%m"))

## DEFINE PROJECTIONS
rovbaseObs.sp <- st_as_sf( rovbaseObs,
                           coords = c( "Oest..UTM33.SWEREF99.TM.",
                                       "Nord..UTM33.SWEREF99.TM."))
st_crs(rovbaseObs.sp) <- st_crs(myStudyArea)

## SUBSET THE DATA TO SAMPLES OF THE SAME SPECIES COLLECTED DURING SAMPLING PERIOD
filter <- list(
  species = "Jerv",
  type = c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)","Saliv/Spytt"),
  month = unlist(myVars$DATA$samplingMonths))

## SUBSET MONTH & TYPE OF SAMPLE
subset <- rovbaseObs.sp$month %in% filter$month & rovbaseObs.sp$Proevetype %in% filter$type
rovbaseObs.sp$monitoring.season <- ifelse(rovbaseObs.sp$month < 12, rovbaseObs.sp$year, rovbaseObs.sp$year+1) #--- need to change for other species
rovbaseObs.sp <- rovbaseObs.sp[subset, ] 

## SUBSET IF SAMPLE WAS SUCCESSFULLY GENOTYPED AND FROM THE FOCAL SPECIES 
subset <- (rovbaseObs.sp$Art..Analyse. %in% filter$species) & !is.na(rovbaseObs.sp$RovbaseID..Analyse.) 
rovbaseObs.sp <- rovbaseObs.sp[-subset, ] 

## SUBSET BASED ON SPACE 
subsetSpace <- !is.na(as.numeric(st_intersects(rovbaseObs.sp, habitat.rWthBufferPol)))
rovbaseObs.sp <- rovbaseObs.sp[subsetSpace, ] 

## RASTERIZE 
r.detector <- aggregate( habitat.subdetectors,
                         fact = myVars$DETECTORS$detResolution/myVars$DETECTORS$detSubResolution)
r.list <- lapply(years, function(y){
  rl <- raster::rasterize( rovbaseObs.sp[rovbaseObs.sp$monitoring.season %in% y,1],
                           r.detector,
                           fun = "count")[[1]]
  rl[is.na(rl[])] <- 0
  rl[!r.detector[]%in% 1] <- NA
  rl1 <- rl
  rl1[rl[]>0] <- 1
  list(rl1, rl)
})
r.OtherSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
r.OtherSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))

## PLOT CHECK
if(myVars$plot.check){
  # pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"mapStructuredOthers",".pdf", sep="" ))
  for(t in 1:nYears){
    year = years[t]
    tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Year%in%year, ]
    #tmpStruct <- myFilteredData.spStructured[myFilteredData.spStructured$Year%in%year, ]
    
    par(mfrow=c(2,2),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]], main=paste(year,"\n Rovbase Samples"), box=F, axes=F)
    plot(st_geometry(tmp), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    #plot(r.OtherSamplesBinary[[t]],main=paste(year,"\n Rovbase Samples Opportunistic"), box=F, axes=F)
    #plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
    
    plot(r.skandObsSamplesBinary[[t]], main=paste(year,"\n SkandObs"), box=F, axes=F)
    plot(st_geometry(tmp), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    #plot(r.skandObsSamplesBinary[[t]],main=paste(year,"\n SkandObs Opportunistic"), box=F, axes=F)
    #plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
  }#t
  # dev.off()
}



## ------          3.2.6.3. COMBINE ROVBASE & SKANDOBS ------

r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
for(t in 1:nYears){
  r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][ ] > 1] <- 1
}#t

## PLOT CHECK
if(myVars$plot.check){
  for(t in 1:nYears){
    par(mfrow=c(1,3),mar=c(0,0,5,0))
    plot(r.OtherSamplesBinary[[t]],main=years[t])
    plot(r.skandObsSamplesBinary[[t]])
    plot(r.SkandObsOtherSamplesBinary[[t]])
  }  
}



## ------          3.2.6.4. SMOOTH THE BINARY MAP ------

## We tried adjust = 0.05, 0.037, 0.02 and decided to go for 0.02 
## ([PD] 0.02 apparently ; not 0.037 as mentionned in previous comment) 
habOwin <- as.owin(as.vector(extent(r.detector)))
ds.list <- lapply(years,function(y){
  ## ROVBASE DATA 
  pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y, ]
  ## SKANDOBS
  pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y, ])
  ## SMOOTH AND RASTERIZE
  p <- ppp(pts[,1], pts[,2], window = habOwin)
  ds <- density(p, adjust = 0.02) #---change bandwith (smoothing) with "adjust
  ds <- raster(ds)
  ds <- ds1 <- raster::resample(ds, r.detector) #mask(ds,rasterToPolygons(myHabitat.list$habitat.rWthBuffer,function(x) x==1))
  threshold <- 0.1 / prod(res(ds)) #--number per 1 unit of the projected raster (meters)
  ds1[] <- ifelse(ds[] < threshold, 0, 1)
  ds1 <- mask(ds1, habitat.rWthBufferPol)
  ds <- mask(ds, habitat.rWthBufferPol)
  return(list(ds,ds1))
})

ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))

names(ds.brick) <- names(ds.brickCont) <- years

## PLOT CHECK
if(myVars$plot.check){
  par(mfrow = c(1,3))
  plot(r.SkandObsOtherSamplesBinary[[t]], main = "Raw Binary", axes = F, box = F)
  plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
  plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
}



## ------          3.2.6.5. ASSIGN THE COVARIATE ------

detOtherSamples <- matrix(0, nrow = n.detectors, ncol = nYears)
detOtherSamples[ ,1:nYears] <- raster::extract( r.SkandObsOtherSamplesBinary,
                                                myDetectors$main.detector.sp)
#colSums(detOtherSamples)



## ------       3.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------

## SCALE COVARIATES
detSnow <- round(scale(detSnow), digits = 2)
detRoads <- round(scale(detRoads), digits = 2)
detTracks <- round(scale(detTracks), digits = 2)

## FORMAT DETECTOR COVARIATES (SYSTEMATIC SAMPLING)
detCovs <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 2))
detCovs[,,1] <- detTracks
detCovs[,,2] <- detSnow

## FORMAT DETECTOR COVARIATES (OPPORTUNISTIC SAMPLING)
detCovsOth <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 3))
detCovsOth[,,1] <- detSnow
detCovsOth[,,2] <- matrix(detRoads,length(detRoads),nYears)
detCovsOth[,,3] <- detOtherSamples

## CHECK IF CONTAINS NAs
if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}

## PLOT CHECK
if(myVars$plot.check){
  tmp <- detectorGrid.r
  par(mfrow = c(2,5), mar = c(0,0,0,0))
  max <- max(detCovsOth[ , ,2])
  cuts <- seq(0,max,length.out = 100)   
  col <- rev(terrain.colors(100))
  
  for(t in 1:nYears){
    plot(detectorGrid.r, col = c(grey(0.2), grey(0.8)), axes = F, legend = F, box = F)
    tmp[!is.na(detectorGrid.r)] <- detCovsOth[ ,t,2]
    plot(tmp, axes = F, legend = F, box = F, breaks = cuts, col = col, add = T)
  }#t
  
  #dev.off()
  #pdf(file=file.path(myVars$WD, myVars$modelName, paste(myVars$modelName,"_detections over space and time.pdf",sep="")))
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
    plot(st_geometry(GLOBALMAP), col = "gray80")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    plot(st_geometry(myBufferedArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
    plot(st_geometry(tempIn), pch = 21, bg = "blue", add = T)
    
    ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    graphics::text(x = 100000, y = 7200000,
                   labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="N"],"NGS"), 
                   cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 100000, y = 7270000,
                   labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"),
                   cex = 1.1, col = "firebrick3", font = 2)
    graphics::text(x = 820000, y = 6780000,
                   labels = paste(NGS_TabTotal[names(NGS_TabTotal)=="S"],"NGS"),
                   cex = 1.1, col = "navyblue", font = 2)
    graphics::text(x = 820000, y = 6850000,
                   labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"),
                   cex = 1.1, col = "navyblue", font = 2)
    ## ADD OVERALL NUMBERS
    mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
    mtext(text = paste(sum(NGS_TabIn), "NGS/", sum(ID_TabIn), "IDs IN"),
          side = 3, line = 0)
    mtext(text = paste(sum(NGS_TabTotal)-sum(NGS_TabIn), "NGS/", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"),
          side = 3, line = -1)
  }#t
}
#dev.off()



## ------     3.3. RESCALE DETECTOR COORDINATES ------

ScaledDetectors <- scaleCoordsToHabitatGrid(
  coordsData = detector.xy,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid = T)$coordsDataScaled



## ------     3.4. CREATE CACHED DETECTORS OBJECTS ------

DetectorIndexLESS <- GetDetectorIndexLESS(
  habitat.mx = myHabitat$habitat.mx,
  detectors.xy = as.matrix(ScaledDetectors),
  maxDist = myVars$DETECTIONS$maxDetDist/res(myHabitat$habitat.r)[1],
  ResizeFactor = myVars$DETECTIONS$resizeFactor,
  plot.check = TRUE)



## ------   4. GENERATE y DETECTION ARRAYS ------
## ------     4.1. ASSIGN DETECTORS ------

## ALIVE
myData.alive <- AssignDetectors_v3sf(
  myData = myFilteredData.sp$alive,                
  myDetectors = myDetectors$main.detector.sp,
  mysubDetectors = myDetectors$detector.sp,
  radius = myVars$DETECTORS$detResolution)

## DEAD RECOVERY
myData.dead <- AssignDetectors_v3sf( 
  myData = myFilteredData.sp$dead.recovery,
  myDetectors = myDetectors$main.detector.sp,
  radius = myVars$DETECTORS$detResolution)



## ------     4.2. FIX DETECTIONS IN NORRBOTTEN ------

# plot(st_geometry(myDetectors$main.detector.sp))
# plot(st_geometry(myDetectors$main.detector.sp[1500, ]),pch=19,add=T,col="red")
# myDetectors$detector.sp %>%
#   filter(main.cell.id == myDetectors$main.detector.sp$main.cell.id[1500]) %>%
#   st_geometry(.) %>%
#   plot(., pch = 19, add = T, col = "green")


### MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET 
### ASSIGNED TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING.
### FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE
### OF NORRBOTTEN

# sum(myData.alive$myData.sp$sub.detector[!myData.alive$myData.sp$Year %in% c(2016:2018)] %in% subdetsNorrbotten)

## Identify detectors inside Norrbotten
COUNTIESAroundNorrbotten <- COUNTIES[COUNTIES$NAME_1 %in% c("Norrbotten",
                                                            "Troms",
                                                            "Västerbotten",
                                                            "Nordland",
                                                            "Finnmark"), ]
COUNTIESAroundNorrbotten <- st_simplify( COUNTIESAroundNorrbotten,
                                         dTolerance = 500)
# distSubdetsCounties <- st_distance( myDetectors$detector.sp,
#                                     COUNTIESAroundNorrbotten,
#                                     byid = TRUE)
# subdetsNorrbotten <- which(apply(distSubdetsCounties, 1, which.min) %in% 3)
# detsNorrbotten <- unique(myDetectors$detector.sp$main.cell.id[subdetsNorrbotten])
distDestsCounties <- st_distance( myDetectors$main.detector.sp,
                                  COUNTIESAroundNorrbotten,
                                  byid = TRUE)
detsNorrbotten <- which(apply(distDestsCounties, 1, which.min) == 3)

## Identify corresponding sub-detectors
subDetsNorrbotten <- which(myDetectors$detector.sp$main.cell.id %in% 
                             myDetectors$main.detector.sp$main.cell.id[detsNorrbotten])

## Identify which detections are assigned to those main detectors in Norrbotten 
## in years where the county was not sampled.
whichDets <- which(!myData.alive$myData.sp$Year %in% c(2016:2018) &
                     myData.alive$myData.sp$Detector %in% detsNorrbotten)

## Loop over flagged detections and assign them to closest detector outside Norrbotten
for(i in 1:length(whichDets)){
  tmp <- myData.alive$myData.sp[whichDets[i], ]
  ## Calculate distance to all sub-detectors
  dist <- st_distance(tmp, myDetectors$detector.sp)
  ## Artificially increase distance for sub-detectors in Norrbotten 
  dist[ ,subDetsNorrbotten] <- 500000
  ## Assign detection to closest sub-detector outside Norrbotten
  myData.alive$myData.sp$sub.detector[whichDets[i]] <- which.min(dist[1, ])
  ## Assign detection to the corresponding main detector outside Norrbotten
  thisDet <- which(myDetectors$main.detector.sp$main.cell.id %in% 
                     myDetectors$detector.sp$main.cell.id[which.min(dist[1, ])])
  myData.alive$myData.sp$Detector[whichDets[i]] <- thisDet
}#i

# plot(st_geometry(myDetectors$main.detector.sp))
# plot(st_geometry(myDetectors$detector.sp[subDetsNorrbotten, ]), col = "green", add = T)
# mainDets <- myDetectors$detector.sp$main.cell.id[subDetsNorrbotten]
# plot(st_geometry(myDetectors$main.detector.sp[myDetectors$main.detector.sp$main.cell.id %in% mainDets, ]),
#      col = "red", add = T)

## SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
sum(myData.alive$myData.sp$Detector[!myData.alive$myData.sp$Year %in% c(2016:2018)] %in% detsNorrbotten)
sum(myData.alive$myData.sp$sub.detector[!myData.alive$myData.sp$Year %in% c(2016:2018)] %in% subDetsNorrbotten)



## ------     4.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------

### It is not possible to distinguish between structured and unstructured samples.
### We assign each sample to structured or unstructured based on whether a given 
### sample matches in time and space with recorded search tracks:
### - a sample is assigned to "structured" if it was collected by the authorities 
### ("Statsforvalteren" or "SNO").
### - a sample is assigned to "structured" if it was located within 500 m from 
### a GPS search track recorded the same day.

distanceThreshold <- 500
whichStructured <- myData.alive$myData.sp$Proeveleverandoer %in% c("Statsforvalteren","SNO","Fylkesmannen") &
  !is.na(myData.alive$myData.sp$TrackRovbsID) &
  myData.alive$myData.sp$TrackDist <= distanceThreshold
myData.aliveStruct <- myData.alive$myData.sp[whichStructured, ]
myData.aliveOthers <- myData.alive$myData.sp[!whichStructured, ]

## CHECK IF A SAMPLE IS NOT MISSING SOMEWHERE
nrow(myData.aliveStruct) + nrow(myData.aliveOthers)
nrow(myData.alive$myData.sp)

## PLOT CHECK
if(myVars$plot.check){
  
  # pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"DetectionsStructuredOppBarplot",".pdf", sep="" ))
  ## USE 500m THRESHOLD
  par(mfrow = c(2,1), mar = c(4,4,3,2))
  barplot(rbind(table(myData.aliveStruct$Year),
                table(myData.aliveOthers$Year)),
          beside = T, ylim = c(0,2000),
          col = c(grey(0.2), grey(0.8)), ylab = "Number of samples")
  abline(h = seq(0,2000, by = 500), lty = 2, col = grey(0.8))
  title(main = "500m threshold")
  legend("topleft", fill = c(grey(0.2), grey(0.8)), legend = c("Structured","Other"))
  
  ## USE 2000m THRESHOLD
  distanceThreshold1 <- 2000
  whichStructured2000 <- myData.alive$myData.sp$Proeveleverandoer %in% c("Statsforvalteren","SNO","Fylkesmannen") &
    !is.na(myData.alive$myData.sp$TrackRovbsID) &
    myData.alive$myData.sp$TrackDist <= distanceThreshold1
  myData.aliveStruc2000 <- myData.alive$myData.sp[whichStructured2000,]
  myData.aliveOthers2000 <- myData.alive$myData.sp[!whichStructured2000,]
  barplot(rbind(table(myData.aliveStruc2000$Year),
                table(myData.aliveOthers2000$Year)),
          beside = T, ylim = c(0,2000),
          col = c(grey(0.2), grey(0.8)), ylab = "Number of samples")
  abline(h = seq(0,2000, by = 500), lty = 2, col = grey(0.8))
  title(main = "2000m threshold")
  legend("topleft", fill = c(grey(0.2), grey(0.8)), legend = c("Structured","Other"))
  # dev.off()
  
  ## CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO"
  tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Proeveleverandoer %in%
                                   c("Statsforvalteren","SNO","Fylkesmannen"),]
  tab <- table(tmp$Year,tmp$TrackRovbsID,useNA ="always" )
  
  ## plot check
  # pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"DetectionsStructuredOpp",".pdf", sep="" ))
  for(t in 1:nYears){
    par(mar = c(0,0,3,0), mfrow = c(1,3))
    tmp1 <- tmp[tmp$Year%in% years[t],]
    tmpNoTracks <- tmp1[is.na(tmp1$TrackRovbsID), ]
    tmpTracks <- tmp1[!is.na(tmp1$TrackRovbsID), ]
    
    plot( myStudyArea, main = "Structured with track")
    plot( st_geometry(tmpTracks),
          pch = 21, col = "black", cex = 1, bg = "red", add = T)
    
    plot( myStudyArea, main = "Structured without track")
    plot( st_geometry(tmpNoTracks),
          pch = 21, col = "black", cex = 1, bg = "blue", add = T)
    
    tmpOpp <- myData.alive$myData.sp[!myData.alive$myData.sp$Proeveleverandoer %in% c("Statsforvalteren","SNO"), ]
    tmpOpp <- tmpOpp[tmpOpp$Year %in% years[t], ]
    
    plot(myStudyArea, main="Other samples")
    plot(st_geometry(tmpOpp), pch=21, col="black", cex=1,bg="green",add=T)
    mtext(years[t],adj = -0.8,padj = 1)
  }#t
  
  barplot(tab[ ,which(is.na(colnames(tab)))]/rowSums(tab),
          main = "% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track")
  # dev.off()
  
  ### plot check
  # pdf(file = paste(myVars$WD,"/",myVars$modelName,"/",myVars$modelName,"OverallDetectionsDeadRecoveries",".pdf", sep="" ))
  plot(st_geometry(GLOBALMAP))
  plot(st_geometry(myStudyArea), add = T)
  plot(st_geometry(myFullData.sp$alive),
       pch = 16, col = "red", cex = 0.3, add = T)
  plot(st_geometry(myFullData.sp$dead.recovery),
       pch = 16, col = "blue", cex = 0.3, add = T)
  mtext(paste("Live detections", nrow(myFullData.sp$alive),
              "; ID:", length(unique(myFullData.sp$alive$Id))
  ),line = +1)
  mtext(paste("Dead recovery:",nrow(myFullData.sp$dead.recovery)))
  # dev.off()
}



## ------     4.4. GENERATE DETECTION HISTORIES : y.alive[i,j,t] & y.dead[i,t] ------

## ALL SAMPLES
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 myDetectors = myDetectors$main.detector.sp,
                 method = "Binomial",
                 myData2 = myData.dead,
                 myDetectors2 = myDetectors$main.detector.sp,
                 returnIdvector = TRUE)
y.ar.ALIVE <- y.ar$y.ar
dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)

## STRUCTURED
y.arStruc <- MakeYsf( myData = myData.aliveStruct,
                      myDetectors = myDetectors$main.detector.sp,
                      method = "Binomial",
                      myData2 = myData.dead,
                      myDetectors2 = myDetectors$main.detector.sp,
                      returnIdvector = TRUE)
y.ar.ALIVEStruc <- y.arStruc$y.ar
dimnames(y.ar.ALIVEStruc) <- dimnames(y.arStruc$y.ar)

## OTHERS
y.arOth <- MakeYsf( myData = myData.aliveOthers,
                    myDetectors = myDetectors$main.detector.sp,
                    method = "Binomial",
                    myData2 = myData.dead,
                    myDetectors2 = myDetectors$main.detector.sp,
                    returnIdvector = TRUE)
y.ar.ALIVEOth <- y.arOth$y.ar
dimnames(y.ar.ALIVEOth) <- dimnames(y.arOth$y.ar)

## MAKE SURE THE Y HAVE THE SAME DIMENSIONS
y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- y.ar.ALIVE
y.ar.ALIVEOthers[] <- y.ar.ALIVEStructured[] <- 0

## FILL IN THE Y ARRAYS 
y.ar.ALIVEOthers[dimnames(y.ar.ALIVEOth)[[1]], , ] <- y.ar.ALIVEOth
y.ar.ALIVEStructured[dimnames(y.ar.ALIVEStruc)[[1]], , ] <- y.ar.ALIVEStruc

## PROJECT THE DEATH TO THE NEXT OCCASION.
# y.ar.DEADProjected <- y.ar$y.ar2 
# y.ar.DEADProjected[] <- 0
# for(t in 2:nYears){y.ar.DEADProjected[,,t] <- y.ar$y.ar2[,,t-1]}
y.ar.DEAD <- apply(y.ar$y.ar2, c(1,3), function(x){if(sum(x)>0){which(x>0)}else{0}})
y.ar.DEAD <- cbind(rep(0, dim(y.ar.DEAD)[1]), y.ar.DEAD)
y.ar.DEAD <- y.ar.DEAD[ ,1:nYears]
y.ar.DEAD[y.ar.DEAD>0] <- 1
dimnames(y.ar.DEAD) <- list( dimnames(y.ar$y.ar2)[[1]],
                             dimnames(y.ar$y.ar2)[[3]])



## ------     4.5. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------

distances <- list()
for(t in 1:nYears){
  print(paste("------ ", t ," -------", sep = "" ))
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
        plot( st_geometry(myStudyArea), 
              main = paste("t: ",t,"     i: ", names(affected.ids)[count], sep = ""))
        scalebar( 2*myVars$DETECTIONS$maxDetDist,
                  xy = c(800000,6700000),
                  type = "bar", divs = 2, below = "km",
                  label = c(0,
                            myVars$DETECTIONS$maxDetDist/1000,
                            myVars$DETECTIONS$maxDetDist/500),
                  cex = 0.8, adj = c(0.5,-0.9))
        plot( st_geometry(COUNTRIES), add = T)
        plot( st_geometry(myDetectors$main.detector.sp),
              add = T, col = grey(0.8), cex = 0.3, pch = 19)
        
        tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$Id == dimnames(y.ar.ALIVE)[[1]][i] &
                                         myFilteredData.sp$alive$Year == years[t], ]
        tmp <- tmp[order(tmp$Date), ]
        tmp.xy <- st_coordinates(tmp)
        n.det <- nrow(tmp.xy)
        
        plot(st_geometry(tmp), col = "pink", pch = 16, cex = 1,add=T)
        arrows(x0 = tmp.xy[1:(n.det-1),1],
               y0 = tmp.xy[1:(n.det-1),2],
               x1 = tmp.xy[2:n.det,1],
               y1 = tmp.xy[2:n.det,2],
               length = 0.1, lwd = 1)
        plot(st_geometry(myDetectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]),
             pch = 16, col = "red",add=T)
        
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
  idd <- names(affected.ids)
  for(i in 1:length(idd)){
    detIds <- which(distances[[t]]$y.flagged[idd[i],]>0)
    myData.alive$myData.sp <- myData.alive$myData.sp[!(myData.alive$myData.sp$Id %in% idd[i] &
                                                         myData.alive$myData.sp$Detector %in% detIds &
                                                         myData.alive$myData.sp$Year %in% years[t]), ]
  }#i
}#t



## ------     4.6. GENERATE INDIVIDUAL-LEVEL COVARIATES ------

## ------       4.6.1. PREVIOUS DETECTION ------

## Make matrix of previous capture indicator
already.detected <- MakeTrapResponseCovsf( data = myFullData.sp$alive,
                                           data.dead = myFullData.sp$dead.recovery)

## Subset to focal years
already.detected <- already.detected[ ,dimnames(already.detected)[[2]] %in% dimnames(y.ar.ALIVE)[[3]]]

## Subset to focal individuals
already.detected <- already.detected[dimnames(already.detected)[[1]] %in% dimnames(y.ar.ALIVE)[[1]], ]

## PLOT CHECK
if(myVars$plot.check){
  ## Plot an image of the matrix
  par(mfrow = c(1,1))
  barplot(colSums(apply(y.ar.ALIVE, c(1,3),function(x)any(x>0))))
  barplot(colSums(already.detected), add = TRUE, col = "gray40")
  legend(x = 0, y = 250, legend = c("newly Det", "already Det"),
         fill = c("gray80", "gray40"))
}



# ## ------       4.6.2. AGE ------
# 
# min.age <- age <- precapture <- matrix( NA,
#                                         nrow = dim(y.ar.ALIVE)[1],
#                                         ncol = dim(y.ar.ALIVE)[3],
#                                         dimnames = list(y.ar$Id.vector,years))
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
# }#i
# 
# ## PLOT CHECK
# if(myVars$plot.check){
#   image(t(min.age))
#   image(t(age))
# }
# 
# 
## ------     4.7. MAKE AUGMENTATION ------

## DATA ARRAYS
y.alive <- MakeAugmentation(y = y.ar.ALIVE, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.dead <- MakeAugmentation(y = y.ar.DEAD, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveOthers <- MakeAugmentation(y = y.ar.ALIVEOthers, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
y.aliveStructured <- MakeAugmentation(y = y.ar.ALIVEStructured, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)

## INDIVIDUAL COVARIATES
already.detected <- MakeAugmentation(y = already.detected, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
# dispersalToggle <- MakeAugmentation(y = dispersalToggle, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
# age <- MakeAugmentation(y = age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
# min.age <- MakeAugmentation(y = min.age, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
# precapture <- MakeAugmentation(y = precapture, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)



## ------     4.8. TRANSFORM Y TO SPARSE MATRICES ------
## STRUCTURED
SparseY <- GetSparseY(y.aliveStructured)

## OTHER
SparseYOth <- GetSparseY(y.aliveOthers)



## ------   1. NIMBLE MODEL DEFINITION ------

modelCode <- nimbleCode({
  ##----------------------------------------------------------------------------
  ##-----------------------------## 
  ##------ SPATIAL PROCESS ------##  
  ##-----------------------------## 
  
  ## One dmean/Lambda per state
  dmean0 ~ dunif(0,10) 
  delta.dmean ~ dunif(0,20)
  
  dmean[1] <- dmean0
  dmean[2] <- dmean0
  dmean[3] <- dmean0 + delta.dmean # dispersers have larger lambda
  dmean[4] <- dmean0
  dmean[5] <- dmean0
  
  lambda[1] <- 1/dmean[1]
  lambda[2] <- 1/dmean[2]
  lambda[3] <- 1/dmean[3]
  lambda[4] <- 1/dmean[4]
  lambda[5] <- 1/dmean[5]
  
  ## Spatial covariate effects
  betaDens ~ dnorm(0.0,0.01)
  habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ## [PD] for later: use different 'betaDens' for the different dispersal states
  ## FIRST YEAR
  for(i in 1:n.individuals){
    sxy[i,1:2,1] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
      upperCoords = upperHabCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
    
    ## T > 1 --> location of ACs depends on :
    ##        1) the distance to the AC location in the previous year and 
    ##        2) Spatial covariates
    for(t in 2:n.years){
      sxy[i,1:2,t] ~ dbernppACmovement_exp(
        lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
        upperCoords = upperHabCoords[1:numHabWindows,1:2],
        s = sxy[i,1:2,t-1],
        lambda = lambda[z[i,t-1]],
        baseIntensities = habIntensity[1:numHabWindows],
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max,
        numWindows= numHabWindows)
    }#t
  }#i
  
  
  
  ##----------------------------------------------------------------------------
  ##-------------------------------## 
  ##----- DEMOGRAPHIC PROCESS -----## 
  ##-------------------------------##  
  
  omeg1[1:4] ~ ddirch(alpha[1:4])   ## The Dirichlet distribution is parameterized by the vector α
  
  disp ~ dunif(0,1) 
  set ~ dunif(0,1)
  phi.juv ~ dunif(0,1)  
  phi.disp ~ dunif(0,1)
  phi.ad ~ dunif(0,1)
  
  for(t in 1:n.years1){
    gamma[t] ~ dunif(0,1)
    omega[1,1:5,t] <- c(1-gamma[t],gamma[t],0,0,0)
    omega[2,1:5,t] <- c(0,phi.juv*(1-disp),phi.juv*disp,0,1-phi.juv)  ## Juvenile
    omega[3,1:5,t] <- c(0,0,phi.disp*(1-set),phi.disp*set,1-phi.disp) ## Disperser
    omega[4,1:5,t] <- c(0,0,0,phi.ad,1-phi.ad)                        ## Adult
    omega[5,1:5,t] <- c(0,0,0,0,1)
  }#t
  
  for(i in 1:n.individuals){ 
    z[i,1] ~ dcat(omeg1[1:4]) 
    for(t in 1:n.years1){
      z[i,t+1] ~ dcat(omega[z[i,t],1:5,t]) 
    }#i 								
  }#t 
  
  
  
  ##----------------------------------------------------------------------------
  ##-----------------------------##
  ##----- DETECTION PROCESS -----## 
  ##-----------------------------##
  
  ## Scale parameter
  sigma0 ~ dunif(0,5)                  ## sigma for juveniles
  delta.sigma.disp ~ dunif(0,5)        ## difference in sigma for dispersers
  delta.sigma.ad ~ dunif(0,5)          ## difference in sigma for adults
  
  sigma[1] <- sigma0
  sigma[2] <- sigma0
  sigma[3] <- sigma0 + delta.sigma.disp ## dispersers have larger sigma
  sigma[4] <- sigma0 + delta.sigma.ad   ## adults have larger sigma?
  sigma[5] <- sigma0
  
  ## Covariates coefficients
  for(t in 1:n.years){
    for(c in 1:n.covs){
      betaCovs[c,t] ~ dunif(-5,5)
    }
    for(c in 1:n.covsOth){
      betaCovsOth[c,t] ~ dunif(-5,5)
    }
    betaResponse[t] ~ dunif(-5,5)
    betaResponseOth[t] ~ dunif(-5,5)
  }#t
  
  ## One p0 per state and county
  for(c in 1:n.counties){
    for(t in 1:n.years){
      for(st in 1:5){
        p01[c,t,st] ~ dunif(0,1)
        p0[c,t,st] <- p01[c,t,st] * countyToggle[c,t]  ## toggle counties
      }
    }#t
  }#c
  
  for(c in 1:n.countries){
    for(t in 1:n.years){
      for(st in 1:5){
        p01Oth[c,t,st] ~ dunif(0,1)
        p0Oth[c,t,st] <- p01Oth[c,t,st] * countryToggle[c,t] ## toggle countries
      }
    }#t
  }#c
  
  pResponse ~ dunif(0, 1)
  
  for(i in 1:n.individuals){
    detResponse[i,1] ~ dbern(pResponse)
    for(t in 1:n.years){
      y.alive[i,1:nMaxDetectors,t] ~ dbin_LESS_Cached_MultipleCovResponse(
        sxy = sxy[i,1:2,t],
        sigma = sigma[z[i,t]], ## One sigma per state 
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
        p0State = p0[1:n.counties,t,z[i,t]], 
        detCountries = detCounties[1:n.detectors],
        detCov = detCovs[1:n.detectors,t,1:n.covs],
        betaCov = betaCovs[1:n.covs,t],
        BetaResponse = betaResponse[t],
        detResponse = detResponse[i,t])
      
      y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbin_LESS_Cached_MultipleCovResponse( 
        sxy = sxy[i,1:2,t],
        sigma = sigma[z[i,t]], ## One sigma per state 
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
        p0State = p0Oth[1:n.countries,t,z[i,t]], ## One p0Oth for county & state
        detCountries = detCountries[1:n.detectors,t],
        detCov = detCovsOth[1:n.detectors,t,1:n.covsOth],
        betaCov = betaCovsOth[1:n.covsOth,t],
        BetaResponse = betaResponseOth[t],
        detResponse = detResponse[i,t])
    }#i
  }#t
  
  
  
  ##----------------------------------------------------------------------------
  ##----------------------------------------## 
  ##---------- DERIVED PARAMETERS ----------##
  ##----------------------------------------##
  for(i in 1:n.individuals){
    for(t in 1:n.years){
      isAlive[i,t] <- (z[i,t] == 2) + (z[i,t] == 3) + (z[i,t] == 4)
    }#t
  }#i
  
  for(t in 1:n.years){
    N.juv[t] <- sum(z[1:n.individuals,t] == 2)
    N.disp[t] <- sum(z[1:n.individuals,t] == 3)
    N.ad[t] <- sum(z[1:n.individuals,t] == 4)
    N[t] <- N.juv[t] + N.disp[t] + N.ad[t]
  }#t
})



## ------   2. NIMBLE PARAMETERS ------

nimParams <- c( "dmean0","delta.dmean","dmean","betaDens","lambda",
                "omeg1", "gamma", "phi.juv","phi.disp", "phi.ad","disp","set",
                "sigma0", "delta.sigma.disp", "delta.sigma.ad","sigma", 
                "pResponse","betaResponse","betaResponseOth",
                "p0", "p0Oth","betaCovs","betaCovsOth",
                "N","N.juv","N.disp","N.ad")

nimParams2 <- c("z","sxy")



## ------   3. NIMBLE CONSTANTS ------
countyToggle <- matrix(1, nrow = max(detCounties), ncol = nYears)
countryToggle <- matrix(1, nrow = max(detCountries)+1, ncol = nYears)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)
for(t in yearsNotSampled){
  countyToggle[1,t] <- 0
  countryToggle[3,t] <- 0
}#t

nimConstants <- list( n.individuals = dim(y.alive)[1],
                      n.detectors = dim(y.alive)[2],  
                      n.years = dim(y.alive)[3], 
                      n.years1 = dim(y.alive)[3]-1, 
                      n.covs = dim(detCovs)[3],
                      n.covsOth = dim(detCovsOth)[3],
                      numHabWindows = nrow(ScaledLowerCoords),
                      n.countries = max(detCountries)+1,
                      n.counties = max(detCounties),
                      countyToggle = countyToggle,
                      countryToggle = countryToggle,
                      y.max = dim(habIDCells.mx)[1],
                      x.max = dim(habIDCells.mx)[2],
                      y.maxDet = dim(DetectorIndexLESS$habitatID)[1],
                      x.maxDet = dim(DetectorIndexLESS$habitatID)[2],
                      ResizeFactor = DetectorIndexLESS$ResizeFactor,
                      n.cellsSparse = dim(DetectorIndexLESS$detectorIndex)[1],
                      maxNBDets = DetectorIndexLESS$maxNBDets,
                      nMaxDetectors = SparseY$nMaxDetectors,
                      nMaxDetectorsOth = SparseYOth$nMaxDetectors)



## ------   4. NIMBLE DATA ------

## ------     4.1. RECONSTRUCT KNOWN z ------

z <- apply(y.alive, c(1,3), function(x) any(x>0))
z <- ifelse(z, 2, NA)
z <- t(apply(z, 1, function(zz){
  if(any(!is.na(zz))){
    range.det <- range(which(!is.na(zz)))
    zz[range.det[1]:range.det[2]] <- 2
  }
  return(zz)
}))



## ------     4.2. GENERATE detResponse DATA ------

detResponse <- already.detected 
detResponse[rownames(detResponse) %in% "Augmented",1] <- NA



## ------     4.3. LIST OF NIMBLE DATA ------

nimData <- list( z = z,   
                 y.alive = SparseY$y,
                 yDets = SparseY$yDets,
                 nbDetections = SparseY$nbDetections,
                 y.aliveOth = SparseYOth$y,
                 yDetsOth = SparseYOth$yDets,
                 nbDetectionsOth = SparseYOth$nbDetections,
                 lowerHabCoords = as.matrix(ScaledLowerCoords),
                 upperHabCoords = as.matrix(ScaledUpperCoords), 
                 detCounties = detCounties,
                 detCountries = detCountriesNorb,
                 detCovs = detCovs,
                 detCovsOth = detCovsOth,
                 detResponse = detResponse,
                 denCounts = denCounts,
                 trials = n.trials,
                 alpha = rep(1,4),
                 detector.xy = as.matrix(ScaledDetectors),
                 habitatGrid = habIDCells.mx,
                 detectorIndex = DetectorIndexLESS$detectorIndex,
                 nDetectorsLESS = DetectorIndexLESS$nDetectorsLESS,
                 habitatIDDet = DetectorIndexLESS$habitatID)



## ------   5. NIMBLE INITS ------

## ------     5.1. GENERATE z INITIAL VALUES ------

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



## ------     5.2. GENERATE detResponse INITIAL VALUES ------

## LATENT VARIABLE DET RESPONSE
InitsDetResponse <- nimData$detResponse
InitsDetResponse[is.na(InitsDetResponse)] <- rbinom(sum(is.na(InitsDetResponse)), 1,0.5)
InitsDetResponse[!is.na(detResponse)] <- NA



## ------     5.3. GENERATE sxy INITIAL VALUES ------

## Create a data.frame with all detections of all Individuals detected
## Project death to the next year
myData.deadProj <- myData.dead[ ,c("Id","Year")]
myData.deadProj$Year <- myData.deadProj$Year + 1

## Remove dead reco occuring the last year (not used)
myData.deadProj <- myData.deadProj[!myData.deadProj$Year %in% max(myData.deadProj$Year), ]

## Combine dead recoveries and NGS detections in one file
AllDets <- rbind(myData.alive$myData.sp[ ,c("Id","Year")],
                 myData.deadProj[ ,c("Id","Year")])
AllDetections <- as.data.frame(AllDets)
AllDetsxy <- st_coordinates(AllDets) 
colnames(AllDetsxy) <- c("x","y")

## Rescale detection coordinates
AllDetsxyscaled <- scaleCoordsToHabitatGrid(
  coordsData = AllDetsxy,
  coordsHabitatGridCenter = myHabitat$habitat.xy,
  scaleToGrid = T)$coordsDataScaled

## Generate sxy initial values
AllDetections <- cbind(AllDetections, AllDetsxyscaled)

sxy.init <- getSInits( AllDetections = AllDetections,
                       Id.vector = y.ar$Id.vector,
                       idAugmented = which(rownames(z) %in% "Augmented"),
                       lowerCoords = nimData$lowerHabCoords,
                       upperCoords = nimData$upperHabCoords,
                       habitatGrid = nimData$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal")  



## ------     5.4. LOOP THROUGH INITIAL VALUES & SAVE OBJECT ------

for(c in 1:4){
  
  ## List initial values
  nimInits <- list( "sxy" = round(sxy.init,5),
                    "dmean0" = runif(1,0,10), 
                    "delta.dmean" = runif(1,0,10),
                    "betaDens" = runif(1,-0.1,0.1),
                    
                    "z" = z.init,
                    "omeg1" = c(0.55,0.1,0.05,0.3),
                    "gamma" = runif(dim(y.alive)[3]-1,0,1),
                    "phi.juv" = runif(1,0.3,0.9),
                    "phi.disp" = runif(1,0.3,0.9),
                    "phi.ad" = runif(1,0.3,0.9),
                    "disp" = runif(1,0,1),
                    "set" = runif(1,0,1),
                    
                    "sigma0" = runif(1,0,2),
                    "delta.sigma.disp" = runif(1,0,4),
                    "delta.sigma.ad" = runif(1,0,2),
                    "pResponse"  = runif(1,0.4,0.5),
                    "detResponse" = InitsDetResponse,
                    "p01" = array(runif(nimConstants$n.counties*dim(y.alive)[3]*5,0,0.2),
                                  c(nimConstants$n.counties,dim(y.alive)[3],5)),
                    "p01Oth" = array(runif(nimConstants$n.counties*dim(y.alive)[3]*5,0,0.2),
                                     c(nimConstants$n.counties,dim(y.alive)[3],5)),
                    "betaCovs" = array(runif(dim(detCovs)[3],-0.1,0.1),
                                       c(dim(detCovs)[3],nYears)),
                    "betaCovsOth" = array(runif(dim(detCovsOth)[3],-0.1,0.1),
                                          c(dim(detCovsOth)[3],nYears)),
                    "betaResponseOth" = runif(dim(y.alive)[3],-0.1,0.1),
                    "betaResponse" = runif(dim(y.alive)[3],-0.1,0.1))
  
  ## TEST IF THE LESS RESTRICTION ON DETECTORS WILL WORK
  ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
  idDetected <- which(!rownames(z) %in%"Augmented")
  for(i in 1:length(idDetected)){
    for(t in 1:nimConstants$n.years){
      if(!is.na(nimInits$sxy[i,1,t])){
        SXY <- nimInits$sxy[i, ,t]
      } else { 
        SXY <- nimData$sxy[i, ,t] 
      }
      sxyID <- nimData$habitatID[trunc(SXY[2]/nimConstants$ResizeFactor)+1,
                                 trunc(SXY[1]/nimConstants$ResizeFactor)+1]
      DETECTIndexdetectorIndex <- nimData$detectorIndex[1:nimConstants$n.cellsSparse,
                                                        1:nimConstants$maxNBDets]
      DETECTLESS <- nimData$nDetectorsLESS[1:nimConstants$n.cellsSparse]
      index <- DETECTIndexdetectorIndex[sxyID,1:DETECTLESS[sxyID]]
      
      ## GET NECESSARY INFO
      n.detectors <- length(index)
      YDET <- nimData$yDets[i,1:nimConstants$nMaxDetectors, t]
      YDETOth <- nimData$yDetsOth[i,1:nimConstants$nMaxDetectorsOth, t]
      
      ## RECREATE Y
      if(nimData$nbDetections[i, t] > 0){
        for(j in 1:nimData$nbDetections[i,t]){
          ## check if a detection is out of the "detection window"
          if(sum(YDET[j]==index)==0){
            print(paste("id",i,"t",t,"j",j))
          }
        }#j
      }#if
    }#t
  }#i
  
  # plot(nimData$detector.xy[ ,2] ~ nimData$detector.xy[ ,1])
  # points(nimData$detector.xy[index,2] ~ nimData$detector.xy[index,1], col = "red")
  # points(SXY[2] ~ SXY[1], col = "blue", pch = 16)
  # points(nimData$detector.xy[YDET[1:nimData$nbDetections[i,t]],2] ~
  #          nimData$detector.xy[YDET[1:nimData$nbDetections[i,t]],1],
  #        col = "green", pch = 16)
  # points(nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i,t]],2]~
  #          nimData$detector.xy[YDETOth[1:nimData$nbDetectionsOth[i,t]],1],
  #        col = "purple", pch = 16)
  # plot(st_geometry(COUNTRIES))
  # tmp <- myData.alive$myData.sp[myData.alive$myData.sp$Id %in% row.names(y.ar.ALIVE)[i] &
  #                                 myData.alive$myData.sp$Year %in% years[t], ]
  # plot(st_geometry(tmp), col = "red", add = T)
  # 
  # ## PLOT CHECK 
  # nimInits$sxy <- round(nimInits$sxy, 5)
  
  
  
  ## ------ 8. SAVE NIMBLE INPUT ------
  # save( modelCode,
  #       nimData,
  #       nimConstants,
  #       nimParams,
  #       nimParams2,
  #       nimInits,
  #       y.dead,
  #       file = file.path( myVars$WD_ASP, myVars$modelName_ASP,
  #                         paste0(myVars$modelName_ASP,"Chain",c,".RData")))
}#c



## ------   6. SAVE NECESSARY OBJECTS ------

# save( myHabitat.list = myHabitat,
#       myDetectors = myDetectors,
#       COUNTRIES = COUNTRIES,
#       myStudyArea.poly = myStudyArea,
#       COMMUNES = COMMUNES,
#       COUNTIES_AGGREGATEDSubset = COUNTIES_AGGREGATEDSubset,
#       myFilteredData.sp = myFilteredData.sp,
#       myFullData.sp = myFullData.sp,
#       COUNTIES_AGGREGATED = COUNTIES_AGGREGATED,
#       file = file.path( myVars$WD_ASP,
#                         myVars$modelName_ASP,
#                         "NecessaryObjects.RData"))



##------------------------------------------------------------------------------