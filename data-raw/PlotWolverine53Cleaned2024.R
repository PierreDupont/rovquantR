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
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
library(xtable)
library(nimbleSCR)
library(readxl)
library(dplyr)
library(ggplot2)
library(colorspace)


## ------ SET REQUIRED WORKING DIRECTORIES ------

# source("C:/PROJECTS/RovQuant/Temp/RB/myWorkingDirectories.R")   
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/myWorkingDirectories.R")
source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")             


## ------ SOURCE THE REQUIRED FUNCTIONS ------

sourceDirectory(dir.function, modifiedOnly = FALSE)
# load(file.path(dir.dropbox,"DATA/MISC DATA/age.lookup.table.RData"))
# source("C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/dbin_LESS_Cached_MultipleCov.R")
source("C:/My_documents/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/dbin_LESS_Cached_MultipleCov.R")

##-- LOAD CPP FUNCTIONS 
for(i in list.files(dir.function.cpp)){sourceCpp(filePath(dir.function.cpp,i))}



##------------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----

myVars <- list( 
  WD = file.path(dir.dropbox, "wolverine/2025"),
  modelName = "Test2025",
  years = 2014:2023,
  plot.check = TRUE)

working.dir <- file.path(myVars$WD, myVars$modelName)


##-- Number of years
years <- myVars$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

##-- Name of the female and male models
modelNameF <- "Hunn"
modelNameM <- "Hann"


nburnin <- 10

##-- Extract date from the last cleaned data file
DATE <- getMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolverine")


##------------------------------------------------------------------------------

## ------ I. LOAD & SELECT DATA ------

## ------   1. LOAD SHAPEFILES ------

# ##-- LOAD GLOBAL MAP
# GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) %>%
#   dplyr::filter(area > 80000000) %>%
#   st_crop(., st_bbox(extent(c(-70000,1200000,5100000,8080000))))

# ##-- POLYGONS OF SWEDEN & NORWAY
# COUNTRIES <- GLOBALMAP %>%
#   dplyr::filter(ISO %in% c("SWE","NOR")) %>%
#   group_by(ISO) %>%
#   summarize()

# ##-- POLYGONS OF COMMUNES IN SWEDEN & NORWAY
# COMMUNES_NOR <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp")) ## Communal map of Norway
# COMMUNES_SWE <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp")) ## Communal map of Sweden
# COMMUNES <- rbind(COMMUNES_NOR, COMMUNES_SWE)
# 
# ##-- POLYGONS OF COUNTIES IN SWEDEN & NORWAY
# COUNTIES <- COMMUNES %>%   
#  group_by(NAME_1) %>%
#   summarize()

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


##-- SHAPEFILE OF NEW COUNTIES IN SWEDEN
NewCountySwe <- readOGR(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp"))
plot(NewCountySwe, border = "red")

##-- HABITAT RASTERS
# load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/Habitat20kmNewSweCounties.RData"))
# load(file.path(dir.dropbox,"DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData"))
# 


## ------   2. LOAD NECESSARY OBJECTS ------

##-- Females
load(file.path( working.dir, "nimbleInFiles/Hunn",
                paste0(myVars$modelName, modelNameF,"_Chain1.RData")))
nimDataF <- nimData
nimInitsF <- nimInits

##-- Males
load(file.path( working.dir, "nimbleInFiles/Hann",
                paste0(myVars$modelName, modelNameM,"_Chain1.RData")))
nimDataM <- nimData
nimInitsM <- nimInits

##-- Remove unnecessary objects from memory
rm(list = c("nimInits", "nimData"))
gc(verbose = FALSE)

##-- Load Necessary objects
load(file.path( myVars$WD, myVars$modelName, "NecessaryObjects.RData"))

habitat <- myHabitat.list
detectors <- myDetectors
rm(list = c("myHabitat.list","myDetectors"))

##-- Load filtered data
load(file.path( myVars$WD, myVars$modelName, "myFilteredData.RData"))

   

##------------------------------------------------------------------------------
## ------ II. INPUT SUMMARY ------

## ------   1. OVERALL NUMBERS ------

##-- SOME TALLIES TO CHECK THINGS
##-- NGS
NGS <- myFilteredData.sp$alive #rbind(myFilteredData.spF$alive, myFilteredData.spM$alive)
NGSStructured <- myData.aliveStruc$myData.sp #rbind( myFilteredData.spStructuredF, myFilteredData.spStructuredM)
NGSOther <- myData.aliveOthers$myData.sp #rbind( myFilteredData.spOthersF, myFilteredData.spOthersM)

##-- FOR REPORT SUMMARY
length(NGS$Id)                                ## Number of NGS samples
length(NGS$Id[NGS$Sex == "Hunn"])             ## Number of Female NGS samples
length(NGS$Id[NGS$Sex == "Hann"])             ## Number of Male NGS samples
length(NGS$Id[NGS$Country == "S"])/nrow(NGS)  ## Proportion of samples in Sweden

length(unique(NGS$Id))                        ## Number of individuals
length(unique(NGS$Id[NGS$Sex == "Hunn"]))     ## Number of Female individuals
length(unique(NGS$Id[NGS$Sex == "Hann"]))     ## Number of Male individuals

## Last year
length(NGS$Id[NGS$Year %in% tail(years, n = 1)])                     ## Number of NGS in the last year
length(NGS$Id[NGS$Sex == "Hunn" & NGS$Year %in% tail(years, n = 1)]) ## Number of female NGS in the last year
length(NGS$Id[NGS$Sex == "Hann" & NGS$Year %in% tail(years, n = 1)]) ## Number of Male NGS in the last year

## NGS structured
length(NGSStructured$Id)                            ## Number of structured NGS samples
length(NGSStructured$Id[NGSStructured$Sex=="Hunn"]) ## Number of structured Female NGS samples
length(NGSStructured$Id[NGSStructured$Sex=="Hann"]) ## Number of structured Male NGS samples

## NGS Other
length(NGSOther$Id)                           ## Number of other NGS samples
length(NGSOther$Id[NGSOther$Sex=="Hunn"])     ## Number of other Female NGS samples
length(NGSOther$Id[NGSOther$Sex=="Hann"])     ## Number of other Male NGS samples


##-- DEAD RECOVERY
dead <- myFullData.sp$dead.recovery #rbind(myFullData.sp$dead.recovery, myFullData.spM$dead.recovery)
table(dead$Year)
length(dead$Id)                               ## Number of individuals
length(unique(dead$Id[dead$Sex=="Hunn"]))     ## Number of Female individuals
length(unique(dead$Id[dead$Sex=="Hann"]))     ## Number of Male individuals

# tmpdead <- dead[dead$Year %in% c(2018:2023),]
# tmpdead <- tmpdead[!duplicated(tmpdead$DNAID), ]
# table(tmpdead$Year,tmpdead$Sex)
# table(tmpdead$Year,tmpdead$Month)
# table(tmpdead$Year)
# nrow(table(tmpdead$Year))
# duplicated(tmpdead$DNAID)
# mapview::mapview(tmpdead)
# 
# tmpNGS <- NGS[NGS$Year %in% c(2018:2023), ]
# table(tmpNGS$Year,tmpNGS$Sex)
# table(tmpNGS$Year,tmpNGS$Month)
# table(tmpNGS$Year)
# 
# ###HENRIK CHECK WITH PUBLIC SAMPLES
# idPublic <- c(
# 'D555438',
# 'D556590',
# 'D553322',
# 'D555239',
# 'D554783',
# 'D558440',
# 'D555845',
# 'D555412',
# 'D556343',
# 'D556344',
# 'D556347',
# 'D558235',
# 'D557101')
# 
# which(NGSStructured$DNAID%in%idPublic)
# 
# # # NGSStructured[which(NGSStructured$DNAID%in%idPublic),]@data
# # # NGSOther[which(NGSOther$DNAID%in%idPublic),]@data
# # 
# # idPublic1 <- c('D555438')
# # NGSStructured[which(NGSStructured$DNAID%in%idPublic1),]@data
# # NGSOther[which(NGSOther$DNAID%in%idPublic1),]
# # 
# # mapview(
# # list(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"),
# #         NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]),
# # layer.name = c("Franconian districts", "Franconian breweries")
# # )
# #         
# # mapview(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"))+
# #   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# # 
# # mapview(as(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),"Spatial"))+
# #   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
# # 
# # st_distance(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),
# #             st_as_sf(NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]))
# # 
# # plot(TRACKS_YEAR[[9]][TRACKS_YEAR[[9]]$RovbaseID %in% "T471191",]$geometry)
# # plot(TRACKSSimple_sf[[9]][TRACKSSimple_sf[[9]]$RovbaseID %in% "T471191",]$geometry)
# # plot(TRACKS[TRACKS$RovbaseID %in% "T471191",]$geometry,col="red",add=T)
# # points(NGSStructured[which(NGSStructured$DNAID%in%idPublic),],pch=16)



## ------   2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX ------

## ------     2.1. ALL ------

NGSCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSCountrySEX) <- unlist(lapply(YEARS, function(x) c(x[2],x[2])))
#colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))

NGSCountrySEX[1, ] <- rep(c("F","M"), nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1, nYears*2, by = 2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
    NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- nrow(temp[temp$Country %in% "N", ])
    NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S", ])
    NGSCountrySEX["Total",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  }#t
}

##-- Export .tex table
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
                                    '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))
print(xtable(NGSCountrySEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEX))), collapse = "")),
      #scalebox = .8,
      floating = FALSE,include.colnames=F,
      add.to.row = addtorow,
      file = file.path(working.dir, "tables","NGSCountrySEX.tex"))

##-- Export .csv table
write.csv( NGSCountrySEX, 
           file = file.path(working.dir, "tables", "NGSCountrySEX.csv"))

# ##-- Checks
# sum(as.numeric(NGSCountrySEX["Total",]))
# sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
# sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))



## ------     2.2. PER OBSERVATION PROCESS ------

NGSCountrySEXoBS <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBS) <- c("", rep(c("Norway","Sweden","Total"), each = 2))
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))

NGSCountrySEXoBS[1, ] <- c("", rep(c("F","M"), nYears))
NGSCountrySEXoBS[ ,1] <- c("", rep(c("Structured","Unstructured"), 3))

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(2,nYears*2,by=2)
for(s in 1:2){
  for(t in 1:nYears){
    ## Structured
    tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[1], ye[t] + sex1[s] ] <- nrow(tempStruc[tempStruc$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[1], ye[t] + sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S" | tempStruc$Country %in% "N", ])
    
    ## Other
    tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[2], ye[t] + sex1[s] ] <- nrow(tempOther[tempOther$Country %in% "N", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S", ])
    NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[2], ye[t] + sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S" | tempOther$Country %in% "N", ])
    }#t
}#s

##-- Export .tex table
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0('& \\multicolumn{1}{c}{} & \\multicolumn{2}{c}{',
                             sort(unique(colnames(NGSCountrySEXoBS)[2:ncol(NGSCountrySEXoBS)])),
                             '}\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBS) <- rep("", ncol(NGSCountrySEXoBS))
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""), ncol = 1)
NGSCountrySEXoBS <- data.frame(cbind(multirowadd,NGSCountrySEXoBS))
colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),
                                                                     paste(x,collapse = "/")))))
print(xtable( NGSCountrySEXoBS,
              type = "latex",
              align = paste(c("l",rep("c",ncol(NGSCountrySEXoBS))), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","NGSCountrySEXperObs.tex"))

# ##-- Checks
# sum(as.numeric(NGSCountrySEX["Total",]))
# sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
# sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))
# 
# #PLOT CHECK 
# plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
# plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="red",add=T)
# plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="red",add=T)
# 
# par(mfrow=c(1,2),mar=c(0,0,0,0))
# plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
# plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="#E69F00",add=T)
# plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
# plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="#009E73",add=T)
# 
# dev.off()
# 
# NGSStructured$Year1 <- NGSStructured$Year+1
# NGSOther$Year1 <- NGSOther$Year+1
# 
# barplot(rbind(table(NGSStructured$Year1),table(NGSOther$Year1)),
#         col=c("#E69F00","#009E73"))
# legend("topleft",fill=c("#E69F00","#009E73"),legend=c("Structured","Other") )

##-- GIVE FILE TO HENRIK ([PD] for what????)
tmp <- NGSOther[NGSOther$Year %in% c(2019,2020,2021), ]
tmp
write.csv( tmp,
           file = file.path(working.dir, "tables", "Unstructured2020_2022.csv"))



## ------   3. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------

## ------     3.1. ALL ------

NGSidCountrySEX <- matrix("", ncol = nYears*2, nrow = 4)
row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
colnames(NGSidCountrySEX) <- c(unlist(lapply(YEARS, function(x)c(x[2],x[2]))))
#colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))

NGSidCountrySEX[1,] <- rep(c("F","M"), nYears)
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1, nYears*2, by = 2)
for(s in 1:2){
  for(t in 1:nYears){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex == sex[s], ]
    NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N"]))
    NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
    NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
  }#t
}#s

##-- Export .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))

print(xtable( NGSidCountrySEX,
              type = "latex",
             align = paste(c("l",rep("c",ncol(NGSidCountrySEX))), collapse = "")),
      #scalebox = .8, 
      floating = FALSE,
      include.colnames = F,
      add.to.row = addtorow,
      file = file.path(working.dir, "tables","NGSidCountrySEX.tex"))

##-- Export .csv
write.csv( NGSidCountrySEX,
           file = file.path(working.dir, "tables", "NGSidCountrySEX.csv"))


##-- Export .csv TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
NGSidCountryTotal <- matrix(0, ncol = nYears, nrow = 1)
row.names(NGSidCountryTotal) <- c("Total")
colnames(NGSidCountryTotal) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))
#colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))

for(t in 1:nYears){
  temp <- NGS[NGS$Year == years[t] , ]
  NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
}#t
write.csv( NGSidCountryTotal,
           file = file.path(working.dir, "tables","TotalIdDetected.csv"))

### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR PER SEX



## ------     3.2. PER OBSERVATION PROCESS ------

NGSCountrySEXoBSid <- matrix("", ncol = nYears*2+1, nrow = 7)
row.names(NGSCountrySEXoBSid) <- c("",rep(c("Norway","Sweden","Total"),each=2))
colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))

NGSCountrySEXoBSid[1,] <- c("", rep(c("F","M"), nYears))
NGSCountrySEXoBSid[,1] <- c("", rep(c("Structured","Unstructured"), 3))

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
}#s


##-- Export .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBSid)[2:ncol(NGSCountrySEXoBSid)])),
                                                              '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
colnames(NGSCountrySEXoBSid) <- rep("", ncol(NGSCountrySEXoBSid))
multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
NGSCountrySEXoBSid <- data.frame(cbind(multirowadd, NGSCountrySEXoBSid))
colnames(NGSCountrySEXoBSid) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))

print(xtable(NGSCountrySEXoBSid, type = "latex",
             align = paste(c("l",rep("c",ncol(NGSCountrySEXoBSid))),collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "NGSCountrySEXperObsid.tex"))



## ------   4. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------

DeadidCountrySEX <- matrix(0, ncol = nYears*2+1, nrow = 6)
row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))
#colnames(DeadidCountrySEX) <- c("", unlist(lapply(YEARS, function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))
DeadidCountrySEX[1,] <- c("",rep(c("F","M"),nYears))
DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1, nYears*2, by = 2)

##-- Identify legal mortality causes
MortalityNames <- unique(as.character(dead$DeathCause))
table(as.character(dead$DeathCause))
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("9", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("23", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("28", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Rifle", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("18", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("17", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Uspesifisert", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Fellefangst", MortalityNames)])
# legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Hagle", MortalityNames)])


##-- SEPARATE MORTALITIES
cause <- c("other","legal culling")
for(t in 1:nYears){
  for(s in 1:2){
    for(d in 1:2){
      if(d == 1){
        temp <- dead[dead$Year == years[t] & dead$Sex == sex[s] & !(dead$DeathCause %in% legalCauses), ]
      } else {
        temp <- dead[dead$Year == years[t] & dead$Sex == sex[s] & dead$DeathCause %in% legalCauses, ]
      }
      row <- which(rownames(DeadidCountrySEX) == cause[d] & DeadidCountrySEX[,1] == "Norway" )
      DeadidCountrySEX[row,ye[t]+sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "N"]))
      
      row <- which(rownames(DeadidCountrySEX) == cause[d] & DeadidCountrySEX[,1] == "Sweden" )
      DeadidCountrySEX[row,ye[t]+sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
    }#d
    DeadidCountrySEX[6,ye[t]+sex1[s]+1] <- sum(as.numeric(DeadidCountrySEX[2:6,ye[t]+sex1[s]+1]))
  }#s
}#t


##-- summary
###-- Other causes
sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
##-- legal
sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))
sum(as.numeric(DeadidCountrySEX[c(2,3),2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))

##-- %of dead reco (legal) in Norway
sum(as.numeric(DeadidCountrySEX[4,2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(4,5),2:ncol(DeadidCountrySEX)]))

##-- ?? 
sum(as.numeric(DeadidCountrySEX[c(3,5),2:ncol(DeadidCountrySEX)]))/
  sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))

##-- ??
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="M")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="F")]))
sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,] %in% c("F","M"))]))

##-- Export. tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
                                                                     '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
##-- REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC
multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))

print(xtable(DeadidCountrySEX, type = "latex",
             align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
      #scalebox = .7, 
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "DeadidCountrySEX.tex"))

# check 
# tmp <- dead[dead$Year == 2019 & 
#               !dead$DeathCause %in% legalCauses &
#               dead$Country %in% "N" & 
#               dead$Sex %in% "Hunn", ]
# length(unique(tmp$Id))
# unique(tmp$Id)
# DeadidCountrySEX[,"X2022"]
# DeadidCountrySEX[,"X2022.1"]
# dead[dead$Id %in% "JI416817 Ind7303 +",]$Sex
# plot(COUNTRIES$geometry)
# plot(tmp$geometry,add=T,col="red",pch=16)



## ------   5. GET THE DETECTED INDIVIDUALS ------

n.detected <- read.csv(file.path(working.dir, "tables", "TotalIdDetected.csv"))
n.detected <- n.detected[1,2:ncol(n.detected)]



## ------   6. SUMMARY DETECTED INDIVIDUALS PER COUNTIES ------

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
# ## NSAMPLES
# summa <- myFilteredData.sp$alive %>%
#   group_by(counties1,Year) %>%
#   summarise(n=n()) %>%
#   st_drop_geometry()
# summa <- summa[summa$counties1>0,]
# ##-- 2023
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
# ## NdetectionsPerID 
# #COUNT NUMBER IDS 
# summa$n1 <- summa$NID <- 0
# dpt <- unique(unlist(summa$counties1))
# yearsss <- c(2022,2023)
# for(t in 1:length(yearsss)){
#   for(i in 1:length(dpt)){
#     tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$counties1 %in% dpt[i] & myFilteredData.sp$alive$Year %in% yearsss[t],]
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$n1 <- nrow(tmp)
#     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$NID <-  length(unique(tmp$Id))
#   }
# }
# 
# summa$NdetPerIDDet <- summa$n1/summa$NID
# 
# ##-- 2023
# tmp <- summa[summa$Year %in% 2023,]
# COUNTIES_AGGREGATED$detPerID2023 <- tmp$NdetPerIDDet#[1:8,"n"]
# #2022
# tmp1 <- summa[summa$Year %in% 2022,]
# COUNTIES_AGGREGATED$detPerID2022 <- tmp1$NdetPerIDDet#[1:8,"n"]
# 
# tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("detPerID2022","detPerID2023")])
# bar <- t(as.matrix(tmppp[c(5,7,8),]))
# colnames(bar) <- c(5,7,8)
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



## ------   7. EXTRACTION REGIONS MAP ------

##-- only select Norwegian counties
COMMUNESNOR <- COMMUNES[COMMUNES$NAME_0 == "Norway", ]
NORWAY <- aggregate(x = COMMUNESNOR, by = "NAME_1")
plot(NORWAY)
NAME_1 <- as.character(NORWAY$NAME_1)
df.CountiesRegions <- matrix(c(
  "Finnmark",         8,
  "Troms",            8,
  "Nordland",         7,
  NAME_1[15],         6,
  NAME_1[9],          6,
  NAME_1[8],          6,
  "Hedmark",          5,
  "Oppland",          3,
  NAME_1[2],          4,
  "Oslo",             4,
  "Akershus",         4,
  "Sogn og Fjordane", 1,
  "Hordaland",        1,
  "Rogaland",         1,
  "Vest-Agder",       1,
  "Aust-Agder",       2,
  "Telemark",         2,
  "Buskerud",         2,
  "Vestfold",         2),
  byrow = T, ncol = 2)

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

pdf(file = file.path(working.dir, "figures", "RegionMaps.pdf"),
    width = 9, height = 13, pointsize = 12)
par(mar=c(0,0,0,0))
plot(gSimplify(COUNTIESsimp, tol=5000), col=NA, border=NA)

##-- NORWAY
CARNIVORE.REGIONS <- COUNTIESsimp[COUNTIESsimp$Country%in% "NOR",]
CARNIVORE.REGIONS1 <- gSimplify(CARNIVORE.REGIONS, tol=200, topologyPreserve = T)
CARNIVORE.REGIONS1$Region <- CARNIVORE.REGIONS$NAME_1
# region <- gUnaryUnion(CARNIVORE.REGIONS1, id = CARNIVORE.REGIONS1$Region)
region <- gSimplify(NORWAY1, tol=500)

##-- SPECIFY COLOR PALETTE
this.col <- sequential_hcl(1+length(unique(NORWAY1$NAME_1)), "Reds 3")
this.col <- this.col[-length(this.col)]
set.seed(100)
this.col <- sample(this.col)
plot(region, add = T, col=this.col, border=border.col, lwd=1)

##-- SWEDEN
NewCountySwe$NAME_1 <- as.character(NewCountySwe$NAME_1)
NewCountySwe1 <- gSimplify(NewCountySwe, tol = 200, topologyPreserve = T)
swedenCounties2 <- NewCountySwe1
swedenCounties2$NAME_1 <- as.character(NewCountySwe$NAME_1)

##-- Remove Gotland
swedenCounties2 <- swedenCounties2[!(swedenCounties2$NAME_1%in% "Gotlands lÃ¤n"), ]

##-- SPECIFY COLOR PALETTE
this.col <- sequential_hcl(25+length(unique(swedenCounties2$NAME_1)), "Blues 3")
this.col <- this.col[1:length(unique(swedenCounties2$NAME_1))]
set.seed(100)
this.col <- sample(this.col)

##-- Plot
plot( RemoveHolesSp(swedenCounties2),
      add = T, col = this.col, border = border.col, lwd = 1)

##-- Add Norwegian county names labels
CARNIVORE.REGIONS2 <- aggregate(CARNIVORE.REGIONS1, by = "Region")
raster::text( NORWAY1,
              labels = NORWAY1$NAME_1,
              cex = 1.3,
              col = ifelse(NORWAY1$NAME_1 %in% c(5), grey(0.7), grey(0)))

##-- Fix county names in Sweden
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

##-- Fix management region names in Sweden
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
                            "Norra",
                            "Mellersta")
swedenCounties2Regions <- aggregate(swedenCounties2, by = "Region")

##-- Plot Swedish management region borders
plot( RemoveHolesSp(swedenCounties2Regions),
      add = T, border = grey(0.0), lwd = 3)

##-- Add Swedish county names labels
polygonsLabel( swedenCounties2,
               swedenCounties2$NAME_1,
               method = "buffer",
               cex = 1.1,
               col = ifelse( swedenCounties2$NAME_1%in%c("Jämtland"),
                             grey(0.7), grey(0.7)))
dev.off()





##------------------------------------------------------------------------------

## ------ III. PROCESS MODEL OUTPUTS -----

## ------ 1. GET AND COMPILE BITES ------
message("## Processing model MCMC outputs...")

##-- Check that a file with that name does not already exist to avoid overwriting
mcmcTest <- TRUE
if (!overwrite) {
  fileName <- paste0("MCMC_wolverine_", DATE, ".RData")
  if (file.exists(file.path(working.dir, "data", fileName))) {
    message(paste0("A processed MCMC output file named '", fileName, "' already exists in: \n",
                   file.path(working.dir, "data")))
    message("Do you want to proceed and overwrite the existing processed MCMC output file? (y/n) ")
    question1 <- readLines(n = 1)
    if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
      message("Not overwriting existing files...")
      message(paste0("Loading '", fileName, "' instead..."))
      load(file.path(working.dir, "data", fileName))
      mcmcTest <- FALSE
    } else {
      message(paste0("Now overwriting '", fileName,"'.\n"))
    }
  }
}

if(mcmcTest){
  
## ------   2.1. FEMALES ------

##-- Compile MCMC bites
nimOutput_F <- collectMCMCbites( path = file.path(working.dir, "nimbleOutFiles/female"),
                                 burnin = nburnin)

##-- Traceplots
grDevices::pdf(file.path(working.dir, "figures/traceplots_F.pdf"))
plot(nimOutput_F$samples[ ,!is.na(nimOutput_F$samples[[1]][1, ])])
grDevices::dev.off()

##-- Process MCMC output
results_F <- processCodaOutput( nimOutput_F$samples,
                                params.omit = c("sxy","z"))
resultsSXYZ_F <- processCodaOutput(nimOutput_F$samples2)

##-- Remove unnecessary objects from memory
rm(list = c("nimOutput_F"))
gc(verbose = FALSE)

##-- Rescale sxy to the original coordinate system
dimnames(resultsSXYZ_F$sims.list$sxy)[[3]] <- c("x","y")
resultsSXYZ_F$sims.list$sxy <- nimbleSCR::scaleCoordsToHabitatGrid(
  coordsData = resultsSXYZ_F$sims.list$sxy,
  coordsHabitatGridCenter = habitat$habitat.xy,
  scaleToGrid = FALSE)$coordsDataScaled

##-- RESCALE sigma AND dmean TO THE ORIGINAL COORDINATE SYSTEM
results_F$sims.list$sigma <- results_F$sims.list$sigma * raster::res(habitat$habitat.r)[1]
results_F$sims.list$dmean <- results_F$sims.list$dmean * raster::res(habitat$habitat.r)[1]



## ------   2.2. MALES -----

##-- Compile MCMC bites
nimOutput_M <- collectMCMCbites( path = file.path(working.dir, "NimbleOutFiles/male"),
                                 burnin = nburnin)

##-- Traceplots
grDevices::pdf(file.path(working.dir, "figures/traceplots_M.pdf"))
plot(nimOutput_M$samples[ ,!is.na(nimOutput_M$samples[[1]][1, ])])
dev.off()

##-- Process MCMC output
results_M <- processCodaOutput( nimOutput_M$samples,
                                params.omit = c("sxy","z"))
resultsSXYZ_M <- processCodaOutput(nimOutput_M$samples2)

##-- Remove unnecessary objects from memory
rm(list = c("nimOutput_M"))
gc(verbose = FALSE)


##-- RESCALE SXY TO THE ORIGINAL COORDINATE SYSTEM
dimnames(resultsSXYZ_M$sims.list$sxy)[[3]] <- c("x","y")
resultsSXYZ_M$sims.list$sxy <- nimbleSCR::scaleCoordsToHabitatGrid(
  coordsData = resultsSXYZ_M$sims.list$sxy,
  coordsHabitatGridCenter = habitat$habitat.df,
  scaleToGrid = FALSE)$coordsDataScaled

##-- RESCALE sigma AND dmean TO THE ORIGINAL COORDINATE SYSTEM
results_M$sims.list$sigma <- results_M$sims.list$sigma * raster::res(habitat$habitat.r)[1]
results_M$sims.list$dmean <- results_M$sims.list$dmean * raster::res(habitat$habitat.r)[1]



## ------   2.3. COMBINE MALES & FEMALES -----

resultsSXYZ_MF <- resultsSXYZ_M

##-- Get minimum number of iterations between model F and M
minIter <- min(dim(resultsSXYZ_F$sims.list$sxy)[1],
               dim(resultsSXYZ_M$sims.list$sxy)[1])

##-- sxy
resultsSXYZ_MF$sims.list$sxy <- abind::abind(resultsSXYZ_M$sims.list$sxy[1:minIter, , , ],
                                             resultsSXYZ_F$sims.list$sxy[1:minIter, , , ],
                                             along = 2)
dimnames(resultsSXYZ_MF$sims.list$sxy)[[3]] <- c("x","y")

##-- z
resultsSXYZ_MF$sims.list$z <- abind::abind(resultsSXYZ_M$sims.list$z[1:minIter, , ],
                                           resultsSXYZ_F$sims.list$z[1:minIter, , ],
                                           along = 2)

##-- sigma
resultsSXYZ_MF$sims.list$sigma <- cbind(results_M$sims.list$sigma[1:minIter, ],
                                        results_F$sims.list$sigma[1:minIter, ])

##-- sex
resultsSXYZ_MF$sims.list$sex <- rep(c("M","F"),
                                    c(dim(resultsSXYZ_M$sims.list$sxy)[2],
                                      dim(resultsSXYZ_F$sims.list$sxy)[2]))

##-- SAVE AND LOAD DATA
save( results_F, results_M, resultsSXYZ_MF,
      file = file.path( working.dir, "data", paste0("MCMC_wolverine.RData")))

}


## ------       1.4.1. CREATE POLYGON WITHOUT BUFFER ------

##-- MALES
habbRNobuffM <- habitat$habitat.rWthBuffer

##-- FEMALES
habbRNobuffF <- habitat$habitat.rWthBuffer





## ------     1.8. CHECK RHAT ------

results_F$Rhat
results_M$Rhat



##------------------------------------------------------------------------------
## ------   2. EXTRACT DENSITY (5km) ------

##-- WHAT STATES ARE CONSIDERED AS ALIVE IN THE MODEL
alive.states <- c(2) 

##-- years not sampled in Norrbotten
yearsSampledNorrb <- c(2016:2018,2023)
yearsNotSampled <- which(!years %in% yearsSampledNorrb)



## ------     2.1. AC BASED DENSITY (5km) ------

habitat.rWthBuffer <- habitat$habitat.rWthBuffer
habitat.rWthBuffer[habitat.rWthBuffer %in% 0] <- NA

searchedPolygon <- sf::st_as_sf(stars::st_as_stars(habitat.rWthBuffer), 
                                as_points = FALSE, merge = TRUE)
searchedPolygon <- searchedPolygon[searchedPolygon$Habitat > 0, ]
# searchedPolygon <- rasterToPolygons(habitat.rWthBuffer, dissolve = T, function(x) x==1 )

##-- REMOVE THE BUFFER FROM THE HABITAT 
##-- COUNTRIES 
rrCountries <- habitatRasterResolution$`5km`[["Countries"]]
##-- REMOVE FINLAND AND RUSSIA 
rrCountries[rrCountries%in% c(1,3)] <- NA
plot(rrCountries)
rrCountries <- mask(rrCountries, searchedPolygon)
rrCountries <- crop(rrCountries, habitat$habitat.r)
plot(rrCountries)

##-- REGIONS & COUNTIES 
rrRegions <- habitatRasterResolution$`5km`[["Regions"]] 
##-- Deal with the special characters
levels(rrRegions)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <- c(
  "Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
  "Värmland","Örebro","Västmanland","Gävleborg",
  "Västernorrland","Jämtland" ,"Västerbotten")
##-- REMOVE FINLAND & RUSSIA 
rrRegions[habitatRasterResolution$`5km`[["Countries"]]%in% c(1, 3)] <- NA
plot(rrRegions)
rrRegions <- mask(rrRegions, searchedPolygon)
rrRegions <- crop(rrRegions, habitat$habitat.r)

habitatPolygon <- sf::st_as_sf(stars::st_as_stars(habitat$habitat.r), 
                                as_points = FALSE, merge = TRUE)
habitatPolygon <- habitatPolygon[habitatPolygon$Habitat > 0, ]
habitatPolygon5km <- mask(habitatRasterResolution$`5km`[["Habitat"]], habitatPolygon)
habitatPolygon5km <- crop(habitatRasterResolution$`5km`[["Habitat"]], habitat$habitat.r)
plot(rrRegions)
plot(habitatPolygon5km)

##-- SWEDISH REGIONS 
rrRegionsSwe <- habitatRasterResolution$`5km`[["Regions"]]

##-- REMOVE FINLAND AND RUSSIA 
rrRegionsSwe[habitatRasterResolution$`5km`[["Countries"]]%in% c(1,2,3)] <- NA
plot(rrRegionsSwe)
rrRegionsSwe <- mask(rrRegionsSwe, searchedPolygon)
rrRegionsSwe <- crop(rrRegionsSwe, habitat$habitat.r)
plot(rrRegionsSwe)
rrRegionsSwe[]

##-- TRANSFORM INTO MANAGEMENT REGIONS 
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
rrCountiesNor <- crop(rrCountiesNor, habitat$habitat.r)
#remove sweden
rrCountiesNor[rrCountries[]%in%4] <- NA
plot(rrCountiesNor)
plot(rrRegions)
plot(rrCountries)


gc()
##-- GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
##-- COUNTRY
densityInputCountries <- getDensityInput(
  regions = rrCountries, 
  habitat = habitatPolygon5km,
  s = resultsSXYZ_MF$sims.list$sxy,
  plot.check = TRUE)

##-- GET THE OBJECTS TO RUN THE DENSITY FUNCTION 
##-- REGIONS
densityInputRegions <- getDensityInput( 
  regions = rrRegions, 
  habitat = habitatPolygon5km,
  s = resultsSXYZ_MF$sims.list$sxy,
  plot.check = TRUE)

##-- Swedish Regions
densityInputRegionsSwe <- getDensityInput( 
  regions = rrRegionsSwe, 
  habitat = habitatPolygon5km,
  s = resultsSXYZ_MF$sims.list$sxy,
  plot.check = TRUE)

##-- MERGE COUNTRY AND REGION MATRICES TO ALLOW SIMULTANEOUS EXTRACTION
regionID <- rbind( densityInputCountries$regions.rgmx,
                   densityInputRegions$regions.rgmx,
                   densityInputRegionsSwe$regions.rgmx)
row.names(regionID) <- c(row.names(densityInputCountries$regions.rgmx),
                         row.names(densityInputRegions$regions.rgmx),
                         row.names(densityInputRegionsSwe$regions.rgmx))

##-- Norwegian Counties
densityInputRegionsNor <- getDensityInput( 
  regions = rrCountiesNor, 
  habitat = habitatPolygon5km,
  s = resultsSXYZ_MF$sims.list$sxy,
  plot.check = TRUE)



## ------       2.1.1. GET THE % OF REGIONS COVERED BY THE ANALYSIS ------

## Percentage of each county included in the analysis
## GET THE PERCENTAGE FOR ALL REGIONS 
rrRegions1 <- habitatRasterResolution$`5km`[["Regions"]] 
##-- Deal with the special characters
levels(rrRegions1)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <- c(
  "Södermanland", "Östergötland","Jönköping", "Skåne", "VästraGötaland",
  "Värmland","Örebro","Västmanland","Gävleborg",
  "Västernorrland","Jämtland" ,"Västerbotten")
##-- CALCULATE AREA OF EACH COUNTY
AreaStudiedRegion <- table(factorValues(rrRegions, rrRegions[]))*res(rrRegions)[1]*1e-6
TotalArea <- table(factorValues(rrRegions1, rrRegions1[]))*res(rrRegions1)[1]*1e-6

##-- Percentage of counties included in the analysis
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
levels(rrRegionsSwe1)[[1]] <- data.frame( "ID" = c(1,2,3),
                                          "Regions" = c("Nordre","Midtre","SÃ¸ndre"))
plot(rrRegionsSwe1)

##-- CALCULATE AREA OF EACH UNIT
AreaStudiedSwe <- table(factorValues(rrRegionsSwe, rrRegionsSwe[]))*res(rrRegionsSwe)[1]*1e-6
TotalAreaSwe <- table(factorValues(rrRegionsSwe1, rrRegionsSwe1[]))*res(rrRegionsSwe1)[1]*1e-6#Percentage of counties included in the analysis
areaRegionsSwe <- AreaStudiedSwe/TotalAreaSwe[names(AreaStudiedSwe)]

##-- GET THE PERCENTAGE FOR THE COUNTRIES 
rrCountries1 <- habitatRasterResolution$`5km`[["Countries"]]
##-- REMOVE FINLAND AND RUSSIA 
rrCountries1[rrCountries1%in% c(1,3)] <- NA
TotalAreaCountry <- table(factorValues(rrCountries1, rrCountries1[]))*res(rrCountries1)[1]*1e-6
AreaStudiedCountry <- table(factorValues(rrCountries, rrCountries[]))*res(rrCountries)[1]*1e-6
TotalAreaCountry["Total"] <- sum(TotalAreaCountry)
AreaStudiedCountry["Total"] <- sum(AreaStudiedCountry)

##-- CALCULATE AREA OF EACH COUNTRY
areaCountry <- AreaStudiedCountry/TotalAreaCountry[names(AreaStudiedCountry)]

##-- MERGE THE PERCENTAGE 
areaAllRegions <- c(areaCountry, areaRegionsSwe, areaRegions)



## ------       2.1.2 MALE & FEMALES (5km) ------

ite <- seq(1,dim(densityInputCountries$sx[ , ,t])[1], by = 1)
 
##-- EXTRACT DENSITY 
DensityCountriesRegions <- list()
for(t in 1:nYears){
  DensityCountriesRegions[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,,t],
    sy =  densityInputCountries$sy[ite,,t],
    z = resultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}#t
DensityCountriesRegions[[t]]$summary

##-- SAVE .RData
save( DensityCountriesRegions,
      file = file.path( dir.dropbox,
                        "wolverine/CM/2021/plot25Cleaned/Figure/DensityCountriesRegions.RData"))



## ------       2.1.3. MALE (5km) ------

##-- Identify Males
IDMales <- which(resultsSXYZ_MF$sims.list$sex == "M")

##-- EXTRACT DENSITY 
DensityCountriesRegionsM <- list()
for(t in 1:nYears){
  DensityCountriesRegionsM[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDMales,t],
    sy =  densityInputCountries$sy[ite,IDMales,t],
    z = resultsSXYZ_MF$sims.list$z[ite,IDMales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}#t
DensityCountriesRegionsM[[t]]$summary

##-- SAVE .RData
save( DensityCountriesRegionsM,
      file = file.path( dir.dropbox,
                       "wolverine/CM/2021/plot25Cleaned/Figure/DensityCountriesRegionsM.RData"))



## ------       2.1.4. FEMALE (5km) ------

##-- Identify Females
IDFemales <- which(resultsSXYZ_MF$sims.list$sex == "F")

##-- EXTRACT DENSITY 
DensityCountriesRegionsF <- list()
for(t in 1:nYears){
  DensityCountriesRegionsF[[t]] <- GetDensity_PD(
    sx = densityInputCountries$sx[ite,IDFemales,t],
    sy =  densityInputCountries$sy[ite,IDFemales,t],
    z = resultsSXYZ_MF$sims.list$z[ite,IDFemales,t],
    IDmx = densityInputCountries$habitat.id,
    aliveStates = alive.states,
    regionID = regionID,
    returnPosteriorCells = F)
}
DensityCountriesRegionsF[[t]]$summary

##-- SAVE .RData
save( DensityCountriesRegionsF,
      file = file.path( dir.dropbox,
                       "wolverine/CM/2021/plot25Cleaned/Figure/DensityCountriesRegionsF.RData" ))



## ------       2.1.5. COUNTIES NORWAY M & F (5km) ------

DensityCountriesRegionsNOR <- list()
for(t in 1:nYears){
  DensityCountriesRegionsNOR[[t]] <- GetDensity_PD(
    sx = densityInputRegionsNor$sx[ite,,t],
    sy =  densityInputRegionsNor$sy[ite,,t],
    z = resultsSXYZ_MF$sims.list$z[ite,,t],
    IDmx = densityInputRegionsNor$habitat.id,
    aliveStates = alive.states,
    regionID = densityInputRegionsNor$regions.rgmx,
    returnPosteriorCells = F)
}#t



## ------       2.1.6. SUMMARY TABLES ------

## ------         2.1.6.1. ALL YEARS, BOTH SEX ------

idcounty <- row.names(DensityCountriesRegions[[t]]$summary)

##-- REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
##-- GET NORWEGIAN VERSUS SWEDISH COUNTIES 
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
                   idcountySWE[idcountySWE%in%CountySouth])

##-- CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimates <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimates) <- c(idcountyTable)
colnames(NCarRegionEstimates) <- unlist(lapply(YEARS ,function(x) x[2]))#

##-- FILL IN THE TABLE 
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

##-- NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")
idcountySWE1 <- idcountySWE
idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
# row.names(NCarRegionEstimates) <- idcounty1
# NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"),yearsNotSampled] <- paste(NCarRegionEstimates[which(idcounty1 %in% "Norrbotten"), yearsNotSampled], "*", sep="")

##-- Export .csv
write.csv( NCarRegionEstimates,
           file = file.path(working.dir, "tables","NAllYears.csv"),
           fileEncoding = "latin1")


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
#                                             paste("\\hspace{0.5cm}", idcountySWE1,sep=""))

print( xtable( NCarRegionEstimates,
              type = "latex",
              align = paste(c("l",rep("c",ncol(NCarRegionEstimates))), collapse = "")),
      # scalebox=.8,
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(NCarRegionEstimates),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "NCountiesCarnivoreRegions.tex"))



## ------         2.1.6.2. LAST YEAR N PER SEX PER COUNTY ------

NCountyEstimatesLastRegions <- matrix("", ncol=3, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste0(
    round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste0(
    round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

## TOTAL 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste0(
    round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

##-- Export .csv
write.csv( NCountyEstimatesLastRegions,
          file = file.path(working.dir, "tables", "NLastYearPerSex.csv"),
          fileEncoding = "latin1")

##-- ADJUST NAMES OF THE TABLE 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
idcounty1[which(idcounty1 %in% "Norway")] <- "NORWAY"
idcounty1[which(idcounty1 %in% "Sweden")] <- "SWEDEN"

idcountySWE1 <- idcountySWE
idcountySWE1 <- sort(idcountySWE1)
idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "\\textcolor[gray]{.5}{Norrbotten*}"

row.names(NCountyEstimatesLastRegions) <- idcounty1
NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"), ] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions[which(idcounty1 %in% "Norrbotten"),], "*}", sep="")
NCountyEstimatesLastRegions["SWEDEN",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",], "**}", sep="")
NCountyEstimatesLastRegions["TOTAL",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",], "**}", sep="")
NCountyEstimatesLastRegions["Nordre",] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",], "**}", sep="")

row.names(NCountyEstimatesLastRegions) <- c(
  "TOTAL",
  paste0("\\hspace{0.25cm}", "NORWAY"),
  paste0("\\hspace{0.50cm}", idcountyNOR),
  paste0("\\hspace{0.25cm}", "SWEDEN"),
  paste0("\\hspace{0.50cm}", "Norra"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth]),
  paste0("\\hspace{0.50cm}", "Mellersta"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle]),
  paste0("\\hspace{0.50cm}", "Södra"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth]))

row.names(NCountyEstimatesLastRegions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLastRegions))] <- "\\hspace{0.75cm}VÃ¤stra GÃ¶taland"
## NECESSARY WITH THE NEW COUNTY DEFINITION IN SWEDEN
# idcounty1 <- str_remove(idcounty1, "lÃ¤n")
# idcounty1 <- str_remove(idcounty1, "s ")
# idcounty1 <- str_remove(idcounty1, " ")

##-- Export .tex 
print(xtable( NCountyEstimatesLastRegions,
              type = "latex",
              align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
      sanitize.text.function = function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions), by = 2)),
                        "\\rowcolor[gray]{.95} "),
      file = file.path(working.dir, "tables","NCountiesSexLastYearRegions.tex"))



## ------         2.1.6.2. LAST YEAR N PER SEX PER COUNTY WITH PROPORTION OF AREA COVERED ------

NCountyEstimatesLastRegions <- matrix("", ncol=4, nrow=length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total","\\% Area")

## FILL IN TABLE 
## FEMALES
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- paste0(
    round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegionsF[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

## MALES 
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- paste0(
    round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegionsM[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

## TOTAL
for( i in 1:length(idcountyTable)){
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- paste0(
    round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"mean"],digits = 1),
    " (",round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
    round(DensityCountriesRegions[[nYears]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}

NCountyEstimatesLastRegions[names(areaAllRegions),"\\% Area"] <- round(areaAllRegions*100,digits = 0)
NCountyEstimatesLastRegions[NCountyEstimatesLastRegions[,4] %in% c("98","99"),4] <- 100
#NCountyEstimatesLastRegions[names(AreaStudied/TotalArea[names(AreaStudied)]),"Area"] <- round(AreaStudied/TotalArea[names(AreaStudied)]*100,digits = 0)

##--  Export .tex
write.csv( NCountyEstimatesLastRegions,
          file = file.path(working.dir, "tables", "NLastYearPerSexArea.csv"),
          fileEncoding = "latin1")

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
# NCountyEstimatesLastRegions["SWEDEN",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["SWEDEN",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["TOTAL",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["TOTAL",1:3], "**}", sep="")
# NCountyEstimatesLastRegions["Nordre",1:3] <- paste("\\textcolor[gray]{.5}{",NCountyEstimatesLastRegions["Nordre",1:3], "**}", sep="")


row.names(NCountyEstimatesLastRegions) <- c(
  "TOTAL",
                                            paste0("\\hspace{0.25cm}","NORWAY"),
                                            paste0("\\hspace{0.5cm} ",
                                                  idcountyNOR,sep=""),
                                            paste0("\\hspace{0.25cm}","SWEDEN"),
                                            paste0("\\hspace{0.5cm}","Norra"),
                                            paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth], sep=""),
                                            paste0("\\hspace{0.5cm}","Mellersta"),
                                            paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle], sep=""),
                                            paste0("\\hspace{0.5cm}","Södra"),
                                            paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth], sep="")
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
      file = file.path(working.dir, "tables",paste("NCountiesSexLastYearRegionsArea.tex")))




## ------         2.1.6.3. MAKE A TABLE 2 last years  ------
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
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Females",years[t])] <- paste0(
      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                      " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Males",years[t])] <- paste0(
      round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }
  
  ## MALES 
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesLast2Regions[idcountyTable[i],paste("Total",years[t])] <- paste0(
      round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
                                                                                    " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
                                                                                    round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
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
write.csv( NCountyEstimatesLast2Regions,
          file = file.path(working.dir, "tables","NLast2YearsPerSex.csv"),
          fileEncoding="latin1")



##ADD LITLE STAR TO NORRBOTTEN
idcountySWE1 <- idcountySWE

idcountySWE1[which(idcountySWE %in% "Norrbotten")] <- "Norrbotten*"
idcountySWE1 <- sort(idcountySWE1)
row.names(NCountyEstimatesLast2Regions) <- idcounty1
NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),] <- paste0(
  NCountyEstimatesLast2Regions[which(idcounty1 %in% "Norrbotten"),], "*")

row.names(NCountyEstimatesLast2Regions) <- c(
  "TOTAL",
  paste0("\\hspace{0.25cm}","NORWAY"),
  paste0("\\hspace{0.5cm} ",idcountyNOR),
  paste0("\\hspace{0.25cm}","SWEDEN"),
  paste0("\\hspace{0.5cm}","Norra**"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyNorth]),
  paste0("\\hspace{0.5cm}","Mellersta"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountyMiddle]),
  paste0("\\hspace{0.5cm}","Södra"),
  paste0("\\hspace{0.75cm}", idcountySWE1[idcountySWE%in%CountySouth])
)
row.names(NCountyEstimatesLast2Regions)[grep("VÃ¤straGÃ¶taland", row.names(NCountyEstimatesLast2Regions))] <- paste("\\hspace{0.75cm}",
                                                                                                                  "VÃ¤stra GÃ¶taland", sep="")
NCountyEstimatesLast2Regions <- rbind(c("F","M","Total","F","M","Total"), NCountyEstimatesLast2Regions)

##-- Export .tex 
addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- c(paste(unlist(YEARS[nYears-1]),collapse = "/"),paste(unlist(YEARS[nYears]),collapse = "/"))
addtorow$command <- c(paste0(paste0('& \\multicolumn{3}{c}{', uniqueYEAR,
                                    '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
print(xtable(NCountyEstimatesLast2Regions, type = "latex",
             align = paste(c("l",rep("c",3),"|",rep("c",3)),collapse = "")),
      sanitize.text.function=function(x){x},
      # scalebox=.8,
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      file = file.path(working.dir, "tables","NCountiesSexLast2YearsRegions.tex"))



## ------         2.1.6.4. ALL YEARS N PER SEX PER COUNTY  ------

NCountyEstimatesAllSexRegions <- matrix("", ncol = nYears*3, nrow = length(idcountyTable)+1)
row.names(NCountyEstimatesAllSexRegions) <- c("",idcountyTable)
colnames(NCountyEstimatesAllSexRegions) <- rep(unlist(lapply(YEARS ,function(x) x[2])),each=3)
NCountyEstimatesAllSexRegions[1,] <- rep(c("Females","Males","Total"),nYears)

## FILL IN TABLE 
for(t in 1:nYears){
  ## FEMALES
  cols <- which(colnames(NCountyEstimatesAllSexRegions) %in% unlist(lapply(YEARS ,function(x) x[2]))[t])
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Females")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste0(
      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
      " (",round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
      round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }
  
  ## MALES 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste0(
      round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
      " (",round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
      round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }
  
  ## TOTAL 
  colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
  for( i in 1:length(idcountyTable)){
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- paste0(
      round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
      " (",round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
      round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }
}

##--  Export .csv
write.csv( NCountyEstimatesAllSexRegions,
          file = file.path(working.dir, "tables", "NAllYearsPerSex.csv"),
          fileEncoding="latin1")

# # ADJUST NAMES OF THE TABLE 
# idcounty1 <- idcountyTable
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
# # WRITE LATEX 
# print(xtable(NCountyEstimatesLastRegions, type = "latex",
#              align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))),collapse = "")),
#       sanitize.text.function=function(x){x},
#       # scalebox=.8,
#       floating = FALSE,
#       add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
#       file = file.path(working.dir, "tables","NCountiesSexLastYearRegions.tex"))



## ------         2.1.6.5. ALL YEARS, BOTH SEX COUNTIES NORWAY ------

idcounty <- row.names(DensityCountriesRegionsNOR[[t]]$summary)
#REMOVE Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
#GET NORWEGIAN VERSUS SWEDISH COUNTIES 
idcountyNOR <- idcounty[grep("Region",idcounty)]
idcountyTable <- c("Total", idcountyNOR)

#CREATE TABLE TO STORE ABUNDANCE AND CI
NCarRegionEstimatesNOR <- matrix("", ncol=nYears, nrow=length(idcountyTable))
row.names(NCarRegionEstimatesNOR) <- c(idcountyTable)
colnames(NCarRegionEstimatesNOR) <- unlist(lapply(YEARS ,function(x) x[2]))#

#FILL IN THE TABLE 
for(t in 1:nYears){
  for( i in 1:length(idcountyTable)){
    NCarRegionEstimatesNOR[idcountyTable[i],t] <- paste0(
      round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"mean"],digits = 1),
      " (",round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
      round(DensityCountriesRegionsNOR[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
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
#save(NCarRegionEstimatesNOR,file=file.path(working.dir, "tables",paste("NAllYearsNorwegianCounties.RData",sep="")))
write.csv( NCarRegionEstimatesNOR,
          file = file.path(working.dir, "tables","NAllYearsNorwegianCounties.csv"),
          fileEncoding = "latin1")

# Encoding(NCarRegionEstimatesNOR[,"name"])[9] <- "ISO-8859-1"
# mb_convert_encoding($file, 'UTF-8', 'ISO-8859-1')
# write.csv2(NCarRegionEstimatesNOR,
#            file = file.path(working.dir, "tables",paste("NAllYearsNorwegianCounties.csv",sep="")),fileEncoding= "UTF-16LE")
# readr::write_excel_csv(NCarRegionEstimatesNOR,
#                         file = file.path(working.dir, "tables",paste("NAllYearsNorwegianCounties.csv",sep="")))
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
      add.to.row = list(list(seq(1,nrow(NCarRegionEstimatesNOR),by=2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "NCountiesCarnivoreRegionsNorway.tex"))



## ------     2.2. PLOT ABUNDANCE TIME SERIES ------

## ------       2.2.2. BARS ------

## ------         2.2.2.1. ALL YEARS ------

SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) #paste(x[[2]]))

SeasonText <- lapply(YEARS, FUN = function(x) x[2]) 

##-- Define colors
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

pdf(file= file.path(working.dir, "figures" , "NCountriesBars.pdf"),
    width = 12, height = 8)
par(mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,1300),
     xlab="", ylab = paste("Estimated number of wolverines"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.4,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

## GET THE DETECTED INDIVIDUALS 
# n.detected <- read.csv(file.path(working.dir, "tables", paste("TotalIdDetected.csv",sep="")))
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




## ------         2.2.2.2. ALL YEARS SEX ------

#SeasonText <- lapply(YEARS, FUN = function(x) paste(x, collapse ="/")) 
SeasonText <- lapply(YEARS, FUN = function(x) x[2]) 

##-- Define colors
text.cex <- 1.5
total.offset <- 37
NO.offset <- -37
SE.offset <- +37
xlim <- c(0.5, nYears + 0.5)

TotalColors <- "black"
country.colors <- c("firebrick2","deepskyblue2")
names(country.colors) <- c("Norway","Sweden")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

pdf( file = file.path(working.dir, "figures", "NCountriesBarsSex.pdf"), width = 18, height = 8)
par( mfrow = c(1,2), mar = c(5,8,3,1),
     las = 1, cex.lab = 2, cex.axis = 1.3, mgp = c(6,2,0),
     xaxs = "i", yaxs = "i")
plot(-1000, xlim = c(0.5, nYears+.5), ylim = c(0,800),
     xlab = "", ylab = "Estimated number of Females", xaxt = "n")
axis(1, at = c(1:(nYears)), labels = SeasonText, cex.axis = 1.1, padj = -1)
at = c(1:nYears)
abline(h = seq(100,1200,by = 100), lty = 2, col = grey(0.90))

##-- GET THE DETECTED INDIVIDUALS 
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE

for(t in 1:nYears){
  ##-- TOTAL
  tmp <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon ),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col = adjustcolor(TotalColors, violin.alpha95),
          border = NA)

  ##-- add a star
  if(sum(t %in% yearsNotSampled)){
    text(x = t + widthPolygon + offsetstar,
         y = quantile95[2],
         "*",
         cex = cexStar)
  }
  ##-- add 50% quantil
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col = adjustcolor(TotalColors, violin.alpha50),
            border = NA)
  }
  
  ##-- SWEDEN
  tmp <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t -widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col = adjustcolor(country.colors[2], violin.alpha95),
          border = NA)
  
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

##-- legend
par(xpd = TRUE)
polygon( x = c(0.8,7.2,7.2,0.8),
         y = c(170,170,230,230),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "white")
labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(50, 50, 50)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
## add transparent background polygon
#polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}

##-- MALES 
#par(mfrow=c(1,2), mar = c(5,8,3,1),las=1, cex.lab=2, cex.axis=1.3, mgp=c(6, 2, 0), xaxs="i", yaxs="i")
plot(-1000, xlim=c(0.5, nYears+.5), ylim=c(0,800),
     xlab="", ylab = paste("Estimated number of Males"), xaxt="n")
axis(1, at=c(1:(nYears)), labels = SeasonText, cex.axis=1.1,padj = -1)
at = c(1:nYears)
abline(h=seq(100,1200,by=100), lty=2, col=grey(0.90))

##-- GET THE DETECTED INDIVIDUALS 
widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15
offsetstar <- 0.05
cexStar <- 1.5
displayQuantiles50 <- TRUE
#yearsNotSampled <- NA

for(t in 1:nYears){
  ## TOTAL
  tmp <- colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"), ])
  quantile95 <- quantile(tmp, prob = c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob = c(0.25, 0.75))
  
  polygon(x = c(t - widthPolygon, t + widthPolygon,
                t + widthPolygon, t - widthPolygon),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col = adjustcolor(TotalColors, violin.alpha95),
          border = NA)
  
  ## add a star
  if(sum(t %in% yearsNotSampled)){
    text( x = t + widthPolygon + offsetstar,
          y = quantile95[2],
          "*", cex = cexStar)
  }
  
  if(displayQuantiles50){
    polygon(x = c(t-widthPolygon, t+widthPolygon,
                  t+widthPolygon, t-widthPolygon ),
            y = c(quantile50[1], quantile50[1],
                  quantile50[2], quantile50[2]), 
            col = adjustcolor(TotalColors, violin.alpha50),
            border = NA)
  }
  
  ## SWEDEN
  tmp <- DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",]
  quantile95 <- quantile(tmp, prob=c(0.0275, 0.975))
  quantile50 <- quantile(tmp, prob=c(0.25, 0.75))
  
  polygon(x = c(t, t - widthPolygon1*2,
                t - widthPolygon1*2, t),
          y = c(quantile95[1], quantile95[1],
                quantile95[2], quantile95[2]), 
          col = adjustcolor(country.colors[2], violin.alpha95),
          border = NA)
  
  if(sum(t %in% yearsNotSampled)){
    text(x = t + offsetstar,
         y = quantile95[2],
         "*", cex = cexStar)
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

## legend
par(xpd=TRUE)
polygon( x = c(0.8,7.2,7.2,0.8),
         y = c(600,600,650,650),
        col = adjustcolor("white",alpha.f = 0.9),
        border="white")
labels <- c(" Norway  ", " Sweden  ", " Total")
pch <- rep(19,4)
cex <- rep(4,4)
y <-  c(630, 630, 630)
x <- c(1,3,5)
mycol1 <- c(country.colors, "black")
# add transparent background polygon
# polygon(c(6.7,8,8,6.7),c(10,10,30,30), col=adjustcolor("white",alpha.f = 0.9), border=NA)
for(i in 1:3){
  #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
  points(x[i],y[i],pch=15,cex=3.5,col=adjustcolor(mycol1[i],violin.alpha95))
  points(x[i],y[i],pch=15,cex=1.5,col=adjustcolor(mycol1[i],violin.alpha50))
  text(x[i]+0.1,y[i],labels[i],cex=1.6,pos=4)
}

dev.off()



## ------         2.2.1.3. LAST YEAR ------

pdf(file= file.path(working.dir, "figures" , "NCountriesBarsLastYear.pdf"), width = 12, height = 8)
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

## ------         2.2.1.4. COMPARE WITH OTHER REPORTS ------

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
#   2,11,"Plot53Cleaned","","2024","2024"),nrow=2,byrow =T)
# 
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
#   ## SWEDEN
#   LSWE <- grep("SWEDEN", L, value = TRUE)
#   LSWE <- gsub("  \\\\hspace", "", LSWE)
#   LSWE <- gsub("  \\\\rowcolor", "", LSWE)
#   LSWE <- gsub("\\\\hspace", "", LSWE)
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
#   ## NORWAY
#   LNOR <- grep("NORWAY", L, value = TRUE)
#   LNOR <- gsub("  \\\\rowcolor", "", LNOR)
#   LNOR <- gsub("\\\\hspace", "", LNOR)
#   LSWE <- gsub("\\\\hspace", "", LSWE)
#   
#   DFNOR <- read.table(text = LNOR, sep = "&", header = F,
#                       strip.white = TRUE, check.names = FALSE, comment.char = "\\")
#   DFNOR[,1]
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



## ------       2.2.3. MAPS ------

habbdensCropped <- list()
max <- max(unlist(lapply(DensityCountriesRegions, function(x) max(x$MeanCell))))
cuts <- seq(0,max,length.out = 100)   #set breaks
col <- rev(terrain.colors(100))

##-- PLOT
pdf(file=file.path(working.dir, "figures", "DensityMapsAC5kms.pdf"))
for(t in 1: nYears){
  habbdens <- densityInputRegions$regions.r
  habbdens[] <- NA
  habbdens[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegions[[t]]$MeanCell
  habbdensCropped[[t]] <- habbdens#crop(habbdens, e.sp)
  
  plot(habbdensCropped[[t]], breaks=cuts, col = col,legend=FALSE, main=years[t]) #p
  plot(nngeo::st_remove_holes(habitat$habitat.poly),add=T, col=NA,border=grey(0.5))
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



## ------     2.3. UD BASED DENSITY (5km) ------

### IDENTIFY PROXIMITY HABITAT CELLS 
habitatMask <- densityInputCountries$habitat.id
habitatMask[!is.na(habitatMask)] <- 1 

## COMPUTE THE UD BASED FOR A FEW ITERATIONS
## RESCALE SIGMA TO METERS 
sigma <- resultsSXYZ_MF$sims.list$sigma#*res(habitat$habitat.r)[1]
##-- RESCALE SIGMA TO THE HABITAT SCALE
sigmaRescaled <- sigma/res(rrCountries)[1]

##-- SELECT xxx ITERATIONS RANDOMLY 
spaceUSED <- list()
iter <- sample(1:dim(densityInputCountries$sy)[1], size = 1000)
# for(t in 1:nYears){
#   spaceUSED[[t]] <- GetSpaceUse(densityInputCountries$sx[iter,,t],
#                                 densityInputCountries$sy[iter,,t],
#                                 resultsSXYZ_MF$sims.list$z[iter,,t],
#                                 sigmaRescaled[iter],#sigmaRescaled[iter],
#                                 densityInputCountries$habitat.xy,
#                                 aliveStates = alive.states,
#                                 regionID = densityInputCountries$regions.rgmx,
#                                 display_progress = T,
#                                 returnPosteriorCells = T)
# }
# 
# spaceUSED1 <- spaceUSED
# spaceUSED <- list()
# for(t in 1:nYears){
#   spaceUSED[[t]] <- list()
#   spaceUSED[[t]][["MeanCell"]] <- spaceUSED1[[t]]$MeanCell
# }
# save(spaceUSED, file = file.path(working.dir, "figures" , "spaceUsed5km.RData" ))
load(file = file.path(working.dir, "figures" , "spaceUsed5km.RData" ))



## ------       2.3.1. PLOT TIME SERIES ------
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


# studDis <- disaggregate(habitat$habitat.poly)
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

#get Norrbotten Counties
#COUNTIESNorrbotten <- gSimplify(COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten",],tol=4000,topologyPreserve = F)
COUNTIESNorrbotten <- st_simplify(COUNTIES[COUNTIES$NAME_1 %in% "Norrbotten",], preserveTopology = F,dTolerance = 4000)


##PLOT
pdf(file=file.path(working.dir, "figures" , "DensityMapsUD.pdf"), width = 12, height = 8)
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
  crs( habbdensUDCropped[[t]]) <- st_crs(habitat$habitat.poly)
  #---BECAUSE raster::plot MESSES UP THE LAYOUT
  image(habbdensUDCropped[[t]], add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
 # plot(RemoveHolesSp(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  plot(st_geometry(COUNTRIESsimpFig), border=grey(0.4), col = NA, add=TRUE)
  
  mtext(SeasonText[[t]], 1, -2, adj=0.15, cex=1.2)
  #box()
  # spPol <- rasterToPolygons(habitat$habitat.r,
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
  # agg <- aggregate(rasterToPolygons(habitat$habitat.r,
  #                                   fun = function(x) x==1))
  # agg <- RemoveHolesSp(agg)
  # 
  # plot(agg, add=TRUE, border="black", col=NA)
   agg1 <- aggregate(rasterToPolygons(habitat$habitat.rWthBuffer,
                                     fun = function(x) x==1))
   #plot(agg1, add=TRUE, border="red", col=NA)
  # plot(RemoveHolesSp(aggregate(habitat$habitat.poly)), add=TRUE, border=grey(0.5), col=NA)
   
  
  
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

## ------       2.3.2. PLOT LAST 2 YEARS ------

## PLOT LAST 2 YEARS  
pdf(file=file.path(working.dir, "figures" , "DensityMapsUDLast2Years.pdf"), 
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
  plot(nngeo::st_remove_holes(habitat$habitat.poly), add=TRUE, border="black", col=NA)
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


## ------       2.3.3. PLOT LAST YEAR ------

## PLOT LAST  YEAR  
pdf(file=file.path(working.dir, "figures", "DensityMapsUDLastYear.pdf"), 
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
  #plot(RemoveHolesSp(aggregate(habitat$habitat.poly)), add=TRUE, border="black", col=NA)
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


## ------       2.3.4. PLOT LAST YEAR SUMMARY ------

## PLOT LAST  YEAR  
pdf(file=file.path(working.dir, "figures" , "DensityMapsUDLastYearSummary.pdf"), 
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
  #plot(RemoveHolesSp(aggregate(habitatF$habitat.poly)), add=TRUE, border="black", col=NA)
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



## ------       2.3.5. PLOT LAST YEAR SUMMARY NO ------
## PLOT LAST  YEAR  
pdf(file=file.path(working.dir, "figures" , paste("DensityMapsUDLastYearSummaryNO.pdf",sep="")), 
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
  #plot(RemoveHolesSp(aggregate(habitatF$habitat.poly)), add=TRUE, border="black", col=NA)
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
    
    plot( habbdensUDCropped[[t]], legend.only=TRUE,breaks=cuts, col=col,
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




## ------       2.3.6. WRITE UD 5km RASTER FOR ROVBASE ------

if(!dir.exists(file.path(working.dir, "rasters"))){
  dir.create(file.path(working.dir, "rasters"))
}

for(t in 1:length(years)){
  raster::crs(habbdensUDCropped[[t]]) <- "EPSG:32633"#st_crs(habitat$habitat.poly))
  path <- file.path(working.dir, "rasters", paste0("wolverine_5km",YEARS[[t]][1],".tif"))
  writeRaster(habbdensUDCropped[[t]], path, overwrite=TRUE)
}



## ------   3. DERIVED PARAMETERS FROM ABUNDANCE ------ 

## ------     3.1. MAKE A GROWTH RATE TABLE PER COUNTRY  ------

growthRate <- matrix(0, ncol = nYears-1, nrow = 3)
row.names(growthRate) <- c("Norway","Sweden","Total")
colnames(growthRate) <- unlist(lapply( YEARS[2:(length(YEARS))],
                                       function(x) paste(x, collapse =  "-")))
for(t in 1:(nYears-1)){
  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Norway", ] /
    DensityCountriesRegions[[t]]$PosteriorRegions["Norway", ]
  growthRate["Norway",t] <- getCleanEstimates(growth)

  growth <- DensityCountriesRegions[[t+1]]$PosteriorRegions["Sweden",] / 
    (DensityCountriesRegions[[t]]$PosteriorRegions["Sweden", ])
  growthRate["Sweden",t] <- getCleanEstimates(growth)
  
  growth <- colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"), ])/
    colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"), ])
  growthRate["Total",t] <- getCleanEstimates(growth)
  
}


### add *** to years and country with OPSCR results 
#row.names(growthRate)[2] <- "Sweden*"
#row.names(growthRate)[3] <- "Total*"
# transitionNotSampled <- colnames(growthRate)
# transitionNotSampled <- unique(unlist(lapply(1:length(years),function(x) grep(as.character((years+1)[yearsNotSampled])[x], transitionNotSampled))))
# transitionNotSampled <- transitionNotSampled[!is.na(transitionNotSampled)] 
# growthRate[2,transitionNotSampled] <- paste(growthRate[2,transitionNotSampled],"*",sep="")
# growthRate[3,transitionNotSampled] <- paste(growthRate[3,transitionNotSampled],"*",sep="")


##-- Print table 
print(xtable( growthRate, type = "latex",
             align = paste(c("l", rep("c",ncol(growthRate))), collapse = "")),
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(growthRate), by = 2)), "\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "growthRate.tex"))

#colSums(DensityCountriesRegions[[t+1]]$PosteriorRegions[c("Sweden","Norway"), ]) 


## ------     3.2. DERIVE SEX RATIO  ------

PropFemale <- PropFemaleSWE <- PropFemaleNOR <- list()
for(t in 1:nYears){
  PropFemale[[t]] <- colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),])/
    (colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"),]) + colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"),]))
  
  PropFemaleSWE[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Sweden",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Sweden",])
  
  PropFemaleNOR[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",]/
    (DensityCountriesRegionsM[[t]]$PosteriorRegions["Norway",] + DensityCountriesRegionsF[[t]]$PosteriorRegions["Norway",])
}

## OVERALL PROPORTION OF FEMALES
mean(unlist(PropFemale))
round(quantile(unlist(PropFemale), probs=c(0.025,0.975)),digits = 2)

## median(PropFemale[[nYears]])
quantile(PropFemale[[nYears]], probs=c(0.025,0.975))
mean(PropFemale[[nYears]])

propFemale_tab <- matrix(0, ncol = nYears, nrow = 3)
row.names(propFemale_tab) <- c("Norway","Sweden","Total")
colnames(propFemale_tab) <- unlist(lapply(YEARS,function(x) c(x[2])))
for(t in 1:nYears){
  ##-- Sweden
  propFemale_tab["Sweden",t] <- getCleanEstimates(PropFemaleSWE[[t]])
  ##-- NORWAY
  propFemale_tab["Norway",t] <- getCleanEstimates(PropFemaleNOR[[t]])
  ##-- Total
  propFemale_tab["Total",t] <- getCleanEstimates(PropFemale[[t]])
}#t

##-- Print table 
print(xtable( propFemale_tab, type = "latex",
              align = paste(c("l", rep("c",ncol(propFemale_tab))), collapse = "")),
      floating = FALSE, sanitize.text.function = function(x){x},
      add.to.row=list(list(seq(1,nrow(propFemale_tab), by = 2)), "\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "propFemale.tex"))



## ------     3.3. DERIVE DENSITY  ------

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

##-- COULD REMAKE ALL THE TABLES WITH DENSITY INSTEAD OF ABUNDANCE....



## ------     3.4.  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED ------

n.detectedTotal <- read.csv(file.path(working.dir, "tables", "TotalIdDetected.csv"))
n.detectedTotal <- as.vector(n.detectedTotal[1, 2:ncol(n.detectedTotal)])

n.detectedSex <- read.csv(file.path(working.dir, "tables", "NGSidCountrySEX.csv"))
tF <- seq(2, length.out = nYears, by=2)
tM <- seq(3, length.out = nYears, by=2)

propDetected <- matrix("", ncol=nYears,nrow=3)
row.names(propDetected) <- c("M","F","Total")
colnames(propDetected) <- unlist(lapply(YEARS,function(x) c(x[2])))

for(t in 1:nYears){
  propDetected["Total",t] <- getCleanEstimates(n.detected[1,t]/colSums(DensityCountriesRegions[[t]]$PosteriorRegions[c("Sweden","Norway"), ]))

  n.detectedSexM <- as.numeric(as.character(n.detectedSex[4,tM[t]]))
  propDetected["M",t] <- getCleanEstimates(n.detectedSexM/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[c("Sweden","Norway"), ]))
  
  n.detectedSexF <- as.numeric(as.character(n.detectedSex[4,tF[t]]))
  propDetected["F",t] <- getCleanEstimates(n.detectedSexF/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[c("Sweden","Norway"), ]))
}#t

##-- PRINT .tex table
print( xtable( propDetected, type = "latex",
               align = paste(c("l",rep("c",ncol(propDetected))), collapse = "")),
       floating = FALSE, sanitize.text.function=function(x){x},
       add.to.row = list(list(seq(1,nrow(propDetected), by = 2)), "\\rowcolor[gray]{.96} "),
       file = file.path(working.dir, "tables", "PropDetectedIds.tex"))



## ------     3.5.  MAKE A TABLE PROPORTION OF INDIVIDUALS DETECTED PER COUNTRIES  ------

n.detectedCountry <- read.csv(file.path(working.dir, "tables", "NGSidCountrySEX.csv"))
colnames(n.detectedCountry) <- c("", unlist(lapply(YEARS,function(x) c(x[2],x[2])))) 
propDetectedCountry <- n.detectedCountry 
propDetectedCountry[2:4, 2:ncol(propDetectedCountry)] <- NA

NCountrySex <- read.csv(file.path(working.dir, "tables", "NAllYearsPerSex.csv"))
colnames(NCountrySex) <- c("", unlist(lapply(YEARS,function(x) c(x[2],x[2],x[2]))))

yrs <- unlist(lapply(YEARS, function(x) c(x[2]))) 

lisCountries <- list()
lisCountries[[1]] <- c("Norway")
lisCountries[[2]] <- c("Sweden")
lisCountries[[3]] <- c("Sweden","Norway")

for(t in 1:nYears){
  cols <- which(colnames(propDetectedCountry) %in% as.character(yrs[t]))
  for(p in 1:2){
    ##-- Female
    propDetectedCountry[p+1,cols[1]] <- getCleanEstimates(
      as.numeric(n.detectedCountry[p+1,cols[1]])/
        DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[p]], ])
    ##-- Male
    propDetectedCountry[p+1,cols[2]] <- getCleanEstimates(
      as.numeric(n.detectedCountry[p+1,cols[2]])/
        DensityCountriesRegionsM[[t]]$PosteriorRegions[lisCountries[[p]], ])
  }#p
  ##-- Female
  propDetectedCountry[4,cols[1]] <- getCleanEstimates(
    as.numeric(n.detectedCountry[4,cols[1]])/
      colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions[lisCountries[[3]], ]))
  ##-- Male
  propDetectedCountry[4,cols[2]] <- getCleanEstimates(
    as.numeric(n.detectedCountry[4,cols[1]])/
      colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions[lisCountries[[3]], ]))
}#t

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
      file = file.path(working.dir, "tables", "propDetectedCountry.tex"))


##-- SPLIT THE TABLE IN TWO 
propDetectedCountry1 <- propDetectedCountry[,c(1:11)]
propDetectedCountry2 <- propDetectedCountry[,c(1,12:21)]
command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[2:11])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(propDetectedCountry)[12:21])),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

##-- SAVE TABLE 1
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
      file = file.path(working.dir, "tables","propDetectedCountry1.tex"))

##-- SAVE TABLE 2
addtorow1$command <- command2
print(xtable(propDetectedCountry2, type = "latex",
             align = paste(rep("c", ncol(propDetectedCountry2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","propDetectedCountry2.tex"))



## ------   4. VITAL RATES  ------

widthPolygon <- 0.15
widthPolygon1 <- 0.15
widthPolygon2 <- 0.15

## ------     4.1. SURVIVAL BARS  ------

pdf( file = file.path(working.dir, "figures", "SurvivalBars.pdf"),
     width = 8, height = 4)

nf <- layout( cbind(c(6,3), c(4,1), c(5,2)),
              widths = c(0.05,1,0.30), heights = c(0.15,1))

par(mar = c(5,4.5,0.5,0.5), tck = 0, xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
plot( 10, xlim = c(0.5, nYears-1+0.5),
      ylim = c(0,1), type ="n", xaxt = "n", xlab = "Years", ylab = "Survival")
axis(2, tck = -0.02)
abline(v = 1:(nYears-1)+0.5, lty = 2)
#axis(1, c(1:nYears), labels = paste(years,years+1,sep=" to\n "),  cex.axis=1.1,padj  = 0.5)
#axis(1, c(1:nYears), labels = paste(years+1,years+2,sep=" to\n "),cex.axis=0.5,padj = -2)
axis( 1, c(1:(nYears)),
      labels = paste(years+1, years+2, sep = " to\n"),
     cex.axis = 1.2, padj = 0.2,tick = F)

myCol <- c("#E69F00","#009E73")

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

## ------     4.2. MORTALITY BARS  ------
# pdf(file=file.path(working.dir, "figures" , paste("Mortality.pdf",sep="")),width=10,height=8)
# 
# nf <- layout(cbind(c(3,4,5),c(7,1,2),c(9,8,6)),widths=c(0.15,1,0.35),heights=c(0.15,1,1))
# # layout.show(nf)
# for(i in c("F","M")){
#   myResults <- Results.list[[i]]
#   
#   par(mar=c(4.5,4.5,1,1),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
#   plot(10, xlim = c(0.5, nYears-1+0.5), ylim = c(0,0.8), type ="n", xaxt="n", xlab = "Years", ylab = "Mortality")
#   axis(2,tck=-0.02)
#   abline(v=1:(nYears-1)+0.5,lty=2)
#   axis(1, c(1:nYears), labels = paste(years,years+1,sep=" to "),  cex.axis=0.8)
#   myCol <- c("green","darkgreen")#c("orange","darkorange2")#c("bisque3","burlywood4")#c("lightgreen","darkgreen")
#   myDev <- c(-0.15,+0.15)
#   
#   for(t in 1:(nYears-1)){
#     plot.violins2(list(myResults$sims.list$h[,t]),
#                   x = t,#+ myDev[s],
#                   at = t+myDev[1],
#                   violin.width = 0.15,
#                   col = myCol[1],#[s],
#                   add = T,
#                   alpha = 0.6,
#                   plot.ci=0.95,
#                   border.col = myCol[1],#[s],
#                   cex=0.8)
#     plot.violins2(list(myResults$sims.list$w[,t]),
#                   x = t,#+ myDev[s],
#                   at = t+myDev[2],
#                   violin.width = 0.15,
#                   col = myCol[2],#[s],
#                   add = T,
#                   alpha = 0.6,
#                   plot.ci=0.95,
#                   border.col = myCol[2],#[s],
#                   cex=0.8)
#   }
#   t
# }#i
# 
# par(mar=c(0,0,0,0))
# plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
# 
# my.labels=c("Female","Male")
# for(i in my.labels){
#   par(mar=c(0.5,0.5,0.5,0.5))
#   plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
#   text(0,0,labels=i,srt=90,cex=2.2,font=2)
# }
# 
# # my.labels=c("Female","Male")
# # for(i in my.labels){
# #    par(mar=c(0.5,0.5,0.5,0.5))
# #    plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
# #    text(0,0,labels=i,srt=0,cex=2.2,font=2)
# # }
# 
# #----LEGEND
# 
# par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
# plot(1,ylim=c(-1,7),xlim=c(0,15),type="n",axes=FALSE)
# cols<-myCol#c("orange","darkorange2")
# pt.col=c("white","white")
# labels<-c("Legal\nculling","Other\nmortality")
# for(i in 1:2){
#   #segments(1,i,2,i,col=cols[i],lwd=10)#,lend=4)
#   points(2,7-i*2,pch=15,cex=3.5,col=adjustcolor(cols[i],0.6))
#   points(2,7-i*2,pch=19,cex=1,col=pt.col[i])
#   text(3.5,7-i*2,labels[i],cex=2,pos=4)
# }
# 
# dev.off()



## DERIVE MORTALITY FROM POSTERIOR AND NUMBER OF DEAD RECOVERIES #
IDF <- row.names(nimDataF$z)
IDM <- row.names(nimDataM$z)

IDMF <- c(IDM, IDF)
load(file.path(dir.dropbox,"wolverine/CM/2024/54.aJ_MaCleaned2024/54.aJ_MaCleaned2024Chain1.RData"))
y.deadM <- nimData$y.dead
load(file.path(dir.dropbox,"wolverine/CM/2024/54.aJ_FaCleaned2024/54.aJ_FaCleaned2024Chain1.RData"))
y.deadF <- nimData$y.dead


y.deadMF <- rbind(y.deadM, y.deadF)


##ARRAY 
itera <- 1:dim(resultsSXYZ_MF$sims.list$z1)[1]#sample(1:dim(resultsSXYZ_MF$sims.list$z1)[1], size = 1000)

MortalityAll  <- array(NA, c(length(itera),nYears-1,2))
MortalityCulled <- array(NA, c(length(itera),nYears-1,2))
MortalityOther <- array(NA, c(length(itera),nYears-1,2))

dimnames(MortalityAll)[[3]] <- dimnames(MortalityCulled)[[3]] <- 
  dimnames(MortalityOther)[[3]] <- c("M","F")

t=1
# tmp <-lapply(1:dim(resultsSXYZ_MF$sims.list$z1)[1],FUN = function(x) sum(resultsSXYZ_MF$sims.list$z1[x,resultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2))
# hist(sum(y.deadM[,t+1])/ unlist(tmp))

for(iter in 1:length(itera) ){
  for(t in 1:(nYears-1)){
    #MALE
    
    sumDead <-  sum(resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2 &
                      resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "M",t+1] %in% 3) 
    sumAlive <- sum(resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2)
    
    MortalityAll[iter,t,1] <- sumDead / sumAlive
    
    MortalityOther[iter,t,1] <- (sumDead  -  sum(y.deadM[,t+1])) /sumAlive
    MortalityCulled[iter,t,1] <- MortalityAll[iter,t,1]- MortalityOther[iter,t,1] 
    
    
    #
    # MortalityCulled[iter,t,1] <- sum(y.deadM[,t+1]) / sum(resultsSXYZ_MF$sims.list$z1[iter,resultsSXYZ_MF$sims.list$sex %in% "M",t] %in% 2)
    # 
    # MortalityOther[iter,t,1] <-  MortalityAll[iter,t,1] - MortalityCulled[iter,t,1]
    # 
    #FEMALE
    sumDead <-  sum(resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "F",t] %in% 2 &
                      resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "F",t+1] %in% 3) 
    sumAlive <- sum(resultsSXYZ_MF$sims.list$z1[itera[iter], resultsSXYZ_MF$sims.list$sex %in% "F",t] %in% 2)
    
    MortalityAll[iter,t,2] <- sumDead / sumAlive
    
    MortalityOther[iter,t,2] <- (sumDead  -  sum(y.deadF[,t+1])) /sumAlive
    MortalityCulled[iter,t,2] <- MortalityAll[iter,t,2]- MortalityOther[iter,t,2] 
    
    
  }
}



# Results.list[["M"]]$sims.list$w  <-  Results.list[["M"]]$sims.list$h  <-   Results.list[["M"]]$sims.list$phi  
# Results.list[["F"]]$sims.list$w  <-  Results.list[["F"]]$sims.list$h  <-   Results.list[["F"]]$sims.list$phi  
#
Results.list[["M"]]$sims.list$w <- MortalityOther[,,1]
Results.list[["M"]]$sims.list$h<- MortalityCulled[,,1]

Results.list[["F"]]$sims.list$w <- MortalityOther[,,2]
Results.list[["F"]]$sims.list$h <- MortalityCulled[,,2]


pdf(file=file.path(working.dir, "figures", "MortalityBars.pdf"),
    width=10,height=8)

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

## ------     4.3. PER CAPITA RECRUITMENT  ------
## ------       4.3.1. PLOT NUMBER OF RECRUITS + PER CAPITA RECRUITMENT ------
pdf(file = file.path(working.dir, "figures", "NbRecruitsPerCapita.pdf"),
    width=10,height=8)

nf <- layout(rbind(c(3,5,6,7),
                   c(3,1,2,4),
                   c(8,1,2,4)),
             widths=c(0.15,1,1,0.1), heights=c(0.15,0.5,0.5))
##PER CAPITA RECRUITMENT 
par(mar=c(4.5,4.5,1,1) ,xaxs="i", cex.axis=1.3, cex.lab=1.6)#
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,1.2), type ="n", xaxt="n", xlab = "Years", ylab = "Per-capita recruitment")
axis(2, tck=-0.02)
axis(1, c(2:(nYears+1)), labels = paste(years, years+1, sep= " to "), 
     cex.axis = 0.8)
abline(v=1:(nYears-1)+0.5,lty=2)
myCol <- grey(0.2)#c("red","blue")

myResults <- resultsSXYZ_MF

for(t in 1:(nYears-1)){
  #available <-apply(myResults$sims.list$z[,,t+1],1,function(x)sum(x%in%c(1)))
  # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
  n.recruit <- apply(myResults$sims.list$z[,,c(t,t+1)], 1, function(x) sum( x[,1] %in% c(1) & x[,2] %in% c(2) ))
  # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
  alivetminus1 <- apply(myResults$sims.list$z[,,t], 1, function(x)sum(x %in% c(2,3)))
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
}#t


##NUMBER OF RECRUITS 
plot(10, xlim = c(0.5, nYears+0.5), ylim = c(0,270), type ="n", xaxt="n", xlab = "Years", ylab = "Number of recruits")
axis(1, c(1:(nYears)),labels = paste(years, years+1,sep=" to "), cex.axis=0.8)
axis(2,tck=-0.02)
abline(v=1:(nYears-1)+0.5,lty=2)


for(t in 1:(nYears-1)){
  n.recruit <- apply(myResults$sims.list$z[,,c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
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
}#t



par(mar=c(0,0,0,0))
plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")

dev.off()

## ------       4.3.2. PLOT NUMBER OF RECRUITS  ------
pdf(file = file.path(working.dir, "figures", "Recruitment.pdf"),
    width=8,height=4)
widthPolygon <- 0.15

# nf <- layout(rbind(c(3,5,6,7),
#                    c(3,1,2,4),
#                    c(8,1,2,4)),
#              widths=c(0.15,1,1,0.1), heights=c(0.15,0.5,0.5))
##PER CAPITA RECRUITMENT 

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),widths=c(0.05,1,0.30),heights=c(0.15,1))

par(mar=c(5,4.5,0.5,0.5),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
#par(mar=c(4.5,4.5,1,1) ,xaxs="i", cex.axis=1.3, cex.lab=1.6)#
myResults <- resultsSXYZ_MF

##NUMBER OF RECRUITS 
plot(10, xlim = c(0.5, nYears-0.5), ylim = c(0,300), type ="n", xaxt="n", xlab = "Years", ylab = "Number of recruits")
axis(1, c(1:(nYears)),labels = paste(years+1, years+2,sep=" to\n"),
     cex.axis=1.2,padj = 0.2,tick = F)
axis(2,tck=-0.02)
abline(v=1:(nYears-1)+0.5,lty=2)
myCol <- c("#E69F00","#009E73")
violin.alpha95 <- 0.3
violin.alpha50 <- 0.7

myDev <- c(-0.16,+0.16)

t=8
n.recruitF <- apply(resultsSXYZ_MF$sims.list$z[,which(resultsSXYZ_MF$sims.list$sex=="F"),c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
n.recruitM <- apply(resultsSXYZ_MF$sims.list$z[,which(resultsSXYZ_MF$sims.list$sex=="M"),c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
 quantile(n.recruitF+n.recruitM, prob=c(0.0275, 0.975))
 quantile(n.recruitF, prob=c(0.0275, 0.975))+quantile(n.recruitM, prob=c(0.0275, 0.975))


for(s in 1:2){
  if(s==1){
    IDSex <- which(resultsSXYZ_MF$sims.list$sex=="F")
  }else{
    IDSex <- which(resultsSXYZ_MF$sims.list$sex=="M")
  }
  
  for(t in 1:(nYears-1)){
    
    n.recruit <- apply(myResults$sims.list$z[,IDSex,c(t,t+1)], 1, function(x)sum( x[,1]%in%c(1) & x[,2]%in%c(2) ))
    # plot.violins2(list(n.recruit),
    #               x = t,
    #               at = t+1,
    #               violin.width = 0.15,
    #               col = myCol,
    #               add = T,
    #               alpha = 0.6,
    #               plot.ci=0.95,
    #               border.col = myCol,
    #               cex=1)
    
    
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

# par(mar=c(0,0,0,0))
# plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")

dev.off()


## ------     4.4. SUMMARY DEMOGRAPHIC RATES ------
parameters <- c("gamma","phi","h","w")
sex <- c("M", "F")
TableState <- matrix(NA, nrow=length(parameters)+1, ncol=(nYears)*2-2)
rownames(TableState) <- c("",unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableState) <- c(unlist(lapply(YEARS[1:(length(YEARS)-1)],function(x) rep(paste(x+1,collapse =  "-"),2))))
TableState[1,] <- c(rep(sex,(nYears-1)) )

for(s in 1:2){
  # choose the sex
  if(s==1){results <- Results.list[["M"]]
  }else{results <- Results.list[["F"]]}
  for(i in 1:length(parameters)){
    # for psi and gamma that are not sate spefici 
    if(length(dim(results$sims.list[[parameters[i]]]))==2){   
      rows <- which(rownames(TableState)==parameters[i])[1]
      col <- which(TableState[1,]==sex[s])
      
      if(length(results$mean[parameters[i]][[1]])==2){
        TableState[rows,col[c(2,5)]] <-  paste( apply(results$sims.list[[parameters[i]]],2, function(x) format(round(median(x),2), nsmall = 2)), #median 
                                                " (", apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.025),2), nsmall = 2)),#UpperCI
                                                "-" , apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.975),2), nsmall = 2)),#Lower CI
                                                ")", sep="")  
        
      }else{
        TableState[rows,col] <-  paste( apply(results$sims.list[[parameters[i]]],2, function(x) format(round(mean(x),2), nsmall = 2)), #median 
                                        " (", apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.025),2), nsmall = 2)),#UpperCI
                                        "-" , apply(results$sims.list[[parameters[i]]],2, function(x) format(round(quantile(x,probs=0.975),2), nsmall = 2)),#Lower CI
                                        ")", sep="")  
        
      }
      
      
      # 3 dimensions vector# 
    }
  }
}

## ADD DERIVED RECRUITMENT 
for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResults_M
  results$sims.list$z <- resultsSXYZ_MF$sims.list$z[,resultsSXYZ_MF$sims.list$sex=="M",]
  
  }else{results <- myResults_F
  results$sims.list$z <- resultsSXYZ_MF$sims.list$z[,resultsSXYZ_MF$sims.list$sex=="F",]
  
  }
  for(t in 1:(nYears-1)){
    n.recruit <- apply(results$sims.list$z[,,c(t,t+1)], 1, function(x) sum( x[,1] %in% c(1) & x[,2] %in% c(2) ))
    # NUMBER of individuals with state 1 at t-1 and 2 at t==> number of recruits
    alivetminus1 <- apply(results$sims.list$z[,,t], 1, function(x)sum(x %in% c(2)))
    temp <- n.recruit/alivetminus1
    col <- which(TableState[1,]==sex[s])
    
    TableState["gamma",col[t]]  <- paste( format(round(median(temp),2), nsmall = 2), #median 
                                          " (",format(round(quantile(temp,probs=0.025),2), nsmall = 2),#UpperCI
                                          "-" ,  format(round(quantile(temp,probs=0.975),2), nsmall = 2),#Lower CI
                                          ")", sep="") 
    
    
  }
}


##
write.csv( TableState,
          file = file.path(working.dir, "figures", "TableParametersState.csv"))
#write latex
addtorow <- list()
addtorow$pos <- list(c(0),0)

colnames(TableState) <- c(unlist(lapply(YEARS[1:(length(YEARS)-1)],function(x) rep(paste(x+1,collapse =  "-"),2))))


addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
#colnames(TableState) <- rep("", ncol(TableState))

rownames(TableState)[2:3] <- c("$\\rho$","$\\phi$")


print(xtable(TableState,
             type = "latex",
             align = paste(rep("c",ncol(TableState)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableParametersState.tex"))


##-- SPLIT THE TABLE IN TWO 
TableState1 <- TableState[,c(1:10)]
TableState2 <- TableState[,c(11:18)]

command1 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState1))),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
command2 <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(TableState2))),
                            '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
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
      file = file.path(working.dir, "tables", "TableParametersState1.tex"))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableState2, type = "latex",
             align = paste(rep("c", ncol(TableState2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "TableParametersState2.tex"))




## ------   5. P0.  ------

widthPolygon <- 0.15

## ------     5.1. BARS ------

pdf(file = file.path(working.dir, "figures", "DetectionProbBars.pdf"),
    width=8,height=10)
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
    myResults <- Results.list[[s]]
    
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



## ------     5.1. BARS Other ------

pdf(file=file.path(working.dir, "figures" , "DetectionProbBarsOther.pdf"),width=8,height=8)
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




## ------     5.2. MAPS ------
# t=7
# #plot(habitatM$habitat.r)
# # pairs<- 2
# # alreadydetected <- 1
# 
# pdf(file=file.path(working.dir, "figures" , paste("MapDetectionProb.pdf",sep="")),width=6,height=6)
# for(t in 1:nYears){
#   
#   p0.r <- habitat$habitat.r
#   P0 <- ilogit(logit(myResults_M$mean$p0[nimDataM$detCountries,t])+
#                  myResults_M$mean$betaResponse*1 +
#                  myResults_M$mean$betaCovs[1]*nimDataM$detCovs[,t,1]+
#                  myResults_M$mean$betaCovs[2]*nimDataM$detCovs[,t,2]+
#                  myResults_M$mean$betaCovs[3]*nimDataM$detCovs[,t,3]
#                
#   )
#   #    
#   cells <- cellFromXY(p0.r,coordinates(detectorsM$main.detector.sp))
#   p0.r[cells] <- P0
#   p0.r <- mask(p0.r,habitat$habitat.poly)
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
## ------     5.3. TABLE ------
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

write.csv( Tablep0,
           file = file.path(working.dir, "figures","Tablep0.csv"))



## ------   6. BETA P0S ------
## ------     6.1. P0STRUCTURED ------
## ------       6.1.1. BETA TRACKS ------
pdf(fil e =file.path(working.dir, "figures" , "Betap0StructuredTracks.pdf"),
    width=8,height=4)

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


## ------       6.1.2. BETA SNOW ------

pdf(file = file.path(working.dir, "figures" , "Betap0StructuredSnow.pdf"),
    width=8,height=4)

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



## ------       6.1.3. BETA RESPONSE ------

pdf(file=file.path(working.dir, "figures" ,"Betap0StructuredResponse.pdf"),width=8,height=4)

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




## ------     6.2. P0_OTHER ------
## ------       6.2.1. BETA SNOW ------

pdf(file=file.path(working.dir, "figures" ,"Betap0OtherSnow.pdf"),width=8,height=4)

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


## ------       6.2.2. BETA ROADS ------

pdf(file=file.path(working.dir, "figures" , "Betap0OtherRoads.pdf",),width=8,height=4)

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



## ------       6.2.3. BETA SKANDOBS ------

pdf(file=file.path(working.dir, "figures" , "Betap0OtherSkandobs.pdf"),width=8,height=4)

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




## ------       6.2.4. BETA RESPONSE ------

pdf(file=file.path(working.dir, "figures" , "Betap0OtherResponse.pdf"),width=8,height=4)

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




## ------   7. TABLE OTHERS  ------

## ------     7.1. TABLE DENSITY & MOVEMENT OPSCR ------

parameters <- c("betaDens","sigma","dmean")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{dens}^*$","$\\sigma$","$\\lambda^*$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableDensityMovementSCR <- matrix(NA, nrow=length(parameters)+1, ncol=(nYears)*2)
rownames(TableDensityMovementSCR) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableDensityMovementSCR) <- unlist(lapply(YEARS,function(x) c(x[2],x[2])))# unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
#c(unlist(lapply(YEARS[1:(length(YEARS))],function(x) rep(paste(x+1,collapse =  "-"),2))))
TableDensityMovementSCR[1,] <- c(rep(sex, (nYears)) )

#ST <- 2:1 ## state 3 had the index 1. state 2 has the index 2. so need to inverse it. 
for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResults_M
  }else{results <- myResults_F}
  
  ### BETA DENSITY 
  # REPEAT ESTIMATES FOR ALL YEARS FOR SUCH PARAMETERS
  param <- "betaDens"
  rows <- which(rownames(TableDensityMovementSCR)==param)[1]
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste0( format(round(median(results$sims.list[[param]]), 2), nsmall = 2), #median
                                               " (", format(round(quantile(results$sims.list[[param]], probs=0.025),2), nsmall = 2),#UpperCI
                                               "-" , format(round(quantile(results$sims.list[[param]], probs=0.975),2), nsmall = 2),#Lower CI
                                               ")")
  ## SIGMA
  # for(st in 1:2){
  param <- "sigma"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste0( format(round(apply(results$sims.list[[param]]/1000,2,median), 2), nsmall = 2), #median
                                               " (", format( round(apply(results$sims.list[[param]]/1000,2,function(x) quantile(x, probs=0.025 )),2), nsmall = 2),#UpperCI
                                               "-" , format(round(apply(results$sims.list[[param]]/1000,2,function(x) quantile(x, probs=0.975 )),2), nsmall = 2),#Lower CI
                                               ")")
  
  ### LAMBDA
  param <- "dmean"
  rows <- which(rownames(TableDensityMovementSCR)==param)
  col <- which(TableDensityMovementSCR[1,]==sex[s])
  
  TableDensityMovementSCR[rows,col] <-  paste0( format(round(median(results$sims.list[[param]]/1000), 2), nsmall = 2), #median
                                               " (", format( round(quantile(results$sims.list[[param]]/1000, probs=0.025 ),2), nsmall = 2),#UpperCI
                                               "-" , format(round(quantile(results$sims.list[[param]]/1000, probs=0.975 ),2), nsmall = 2),#Lower CI
                                               ")")
}

##WRITE TABLES 
write.csv( TableDensityMovementSCR,
          file = file.path(working.dir, "figures", "TableDensityMovement.csv"))
#write latex
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
      sanitize.text.function=function(x){x},
      file = file.path(working.dir, "tables","TableDensityMovement.tex"))


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
      file = file.path(working.dir, "tables",paste("TableDensityMovement1.tex", sep="")))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableDensityMovementSCR2, type = "latex",
             align = paste(rep("c", ncol(TableDensityMovementSCR2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableDensityMovement2.tex"))



## ------     7.2. TABLE PARAMETERS STRUCTURED ------
parameters <- c("betaResponse","betaCovs","betaCovs")#,"betaResponse", "betaTracks","betaRoads", "betaSnow")
parameters1 <- c("$\\beta_{1_Structured}$","$\\beta_{2_Structured}$","$\\beta_{3_Structured}$")#,"$\\beta_1$", "$\\beta_2$", "$\\beta_3$","$\\beta_4$")

n.digits = 2
n.digitsigma = 0

sex <- c("M", "F")
TableStructured <- matrix(NA, nrow=length(parameters1)+1, ncol=(nYears)*2)
rownames(TableStructured) <- c("", unlist(lapply(as.list(parameters),function(x) rep(x,1))))
colnames(TableStructured) <- unlist(lapply(YEARS,function(x) c(x[2],x[2])))# unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
TableStructured[1,] <- c(rep(sex, (nYears)) )



for(s in 1:2){
  # choose the sex
  if(s==1){results <- myResults_M
  }else{results <- myResults_F}
  
  
  
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
write.csv( TableStructured,
           file = file.path(working.dir, "figures" ,"TableStructured.csv"))
#write latex
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

print(xtable(TableStructured,
             type = "latex",
             align = paste(rep("c",ncol(TableStructured)+1),collapse = "")),
      # scalebox=.7,
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableStructured.tex"))


##-- SPLIT THE TABLE IN TWO 
TableStructured1 <- TableStructured[,c(1:10)]
TableStructured2 <- TableStructured[,c(11:20)]

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
#SAVE TABLE 1

addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableStructured1,
             type = "latex",
             align = paste(rep("c", ncol(TableStructured1)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableStructured1.tex"))

#SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableStructured2,
             type = "latex",
             align = paste(rep("c", ncol(TableStructured2)+1), collapse = "")),
      # scalebox = .7,
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "TableStructured2.tex"))



## ------     7.3. TABLE PARAMETERS OTHERS ------
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
  if(s==1){results <- myResults_M
  }else{results <- myResults_F}
  
  
  
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
write.csv( TableOther,
          file = file.path(working.dir, "figures", "TableOther.csv"))
#write latex
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



print(xtable(TableOther,
             type = "latex",
             align = paste(rep("c",ncol(TableOther)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableOther.tex"))


##SPLIT THE TABLE IN TWO 
TableOther1 <- TableOther[ ,c(1:10)]
TableOther2 <- TableOther[ ,c(11:20)]

# uniqueYEAR1 <- uniqueYEAR[1: ((ncol(TableState1)-2)/2)]
# uniqueYEAR2 <- uniqueYEAR[((ncol(TableState2)-2)/2+1): length(uniqueYEAR) ]
#SAVE TABLE 1
addtorow1 <- addtorow
addtorow1$command <- command1
print(xtable(TableOther1,
             type = "latex",
             align = paste(rep("c", ncol(TableOther1)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables","TableOther1.tex"))

# SAVE TABLE 2
addtorow1$command <- command2
print(xtable(TableOther2,
             type = "latex",
             align = paste(rep("c", ncol(TableOther2)+1), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow1,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "TableOther2.tex"))



## ------   8. TRANSITION SURFACES  ------
## ------     8.1. SET UP HABITAT ------
habbR <- raster::disaggregate(habitat$habitat.r, fact=4)

##
#COUNTRIES <- aggregate(x = COMMUNESM, by = "ISO")
COUNTRIES <- COMMUNES %>%    group_by(ISO) %>%summarize()

COUNTRIESsimp <- COUNTRIES#gSimplify(COUNTIES,tol = 100, topologyPreserve = T)
COUNTRIESsimp$NAME_1 <- COUNTRIES$ISO
plot(st_geometry(habitat$habitat.poly))
plot(st_geometry(COUNTRIESsimp), add=T)

## IDENTIFY COUNTIES
myStudyArea.polyMnoHoles <- nngeo::st_remove_holes(habitat$habitat.poly)
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
# e <- extent(habitatM$habitat.poly)
# e.sp <- as(e, 'SpatialPolygons')
##GET THE BUFFERLESS AREA
# habitat.r <- habitatM$habitat.r
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
resultsSXYZ_MF$sims.list$scaledsxy <- scaleCoordsToHabitatGrid(coordsData = resultsSXYZ_MF$sims.list$sxy,
                                                                   coordsHabitatGridCenter = habbRxy)$coordsDataScaled

## ------     8.2. TRANSITION PROBABILITY OTHER CAUSES ------
TransitionSurfaceOther <- list()
iter <- sample(1:dim(resultsSXYZ_MF$sims.list$scaledsxy)[1], size = 100)#dim(densityInputCountries$sx)[1])



#for(t in 6:(nYears-1)){
for(t in 1:(nYears-1)){
   TransitionSurfaceOther[[t]] <- GetTransitionSurface( resultsSXYZ_MF$sims.list$scaledsxy[,IDFemales,1,t],
                                                        resultsSXYZ_MF$sims.list$scaledsxy[,IDFemales,2,t],
                                                        resultsSXYZ_MF$sims.list$z[,IDFemales,t],
                                                        resultsSXYZ_MF$sims.list$z[,IDFemales,t+1],
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
plot(habitat$habitat.poly,add=T)



## ------     8.3. TRANSITION PROBABILITY CULLING  ------
TransitionSurfaceCulling <- list()
for(t in 1:(nYears-1)){
   TransitionSurfaceCulling[[t]] <- GetTransitionSurface( resultsSXYZ_MF$sims.list$scaledsxy[,,1,t],
                                                          resultsSXYZ_MF$sims.list$scaledsxy[,,2,t],
                                                          resultsSXYZ_MF$sims.list$z[,,t],
                                                          resultsSXYZ_MF$sims.list$z[,,t+1],
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
plot(habitat$habitat.poly,add=T)


## ------     8.4. PLOT  ------

pdf( file = file.path(working.dir, "figures", "CountryMortalityRates.pdf"),
     width = 9,
     height = 7)

par(mfrow = c(1,2))
offset <- c(-0.2,0.2)
col <- c("red","blue")

plot(-10, xlim = c(0,nYears), ylim = c(0,1), ylab = "Mortality rate", xaxt = "n")
axis(1, at = c(1:nYears), labels = years)
for(t in 1:(nYears-1)){
   for(r in 1:2){
      plot.violins(list(TransitionSurfaceOther[[t]]$PosteriorTransitionRegion[r, ]),
                   at = t + offset[r],
                   x = 1,
                   col = col[r],
                   alpha = 0.5,
                   add = T)
   }#r
}#t

plot(-10, xlim = c(0,nYears), ylim = c(0,1), ylab = "Mortality rate culling", xaxt = "n")
axis(1, at = c(1:nYears), labels = years)
for(t in 1:(nYears-1)){
   for(r in 1:2){
      plot.violins(list(TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion[r,]),
                   at = t + offset[r],
                   x = 1,
                   col = col[r],
                   alpha = 0.5,
                   add = T)
   }#r
}#t

legend("topright", fill=col, legend=c("Norway","Sweden"))
dev.off()



## ------   9. OTHER PLOTS ------

## ------     9.1. SKANDOBS ------

habitat.detectors <- aggregate(rasterToPolygons(disaggregate(habitat$habitat.rWthBuffer, 
                                                             fact=2),fun=function(x)x==1))
skandobs.r1 <- skandobs.r <- disaggregate(habitat$habitat.rWthBuffer, 
                                         fact=2)

pdf( file= file.path(working.dir, "figures", "SkandobsRovbaseCovariates.pdf"),
     width = 12, height = 8)
par(mfrow = c(2,5), mar = c(1,1,3,1))
rrr <- list()
for(t in 1:nYears){
  plot(habitat.detectors, border = grey(0.8), col = grey(0.8))
  skandobs.r1[skandobs.r[] %in% 1]<- nimData$detCovsOth[ ,t,3]
  rrr[[t]] <- skandobs.r1
  pol <- rasterToPolygons(skandobs.r1,fun = function(x) x==1)
  plot(pol, col = "darkgreen", border = NA, add = T)
  mtext(paste(YEARS[[t]][2]))
}
dev.off()

save(rrr, file = file.path(working.dir, "figures", "Skandobs.RData"))



## ------   10. SAVE ------

# TransitionPosteriorCulling <- 
#    TransitionPosteriorOther <- list()
# for(t in 1:(nYears-1)){
#    TransitionPosteriorCulling[[t]]<- TransitionSurfaceCulling[[t]]$PosteriorTransitionRegion
#    TransitionPosteriorOther[[t]]<- TransitionSurfaceOther[[t]]$PosteriorTransitionRegion
# }
# 
# save( TransitionPosteriorCulling,
#       TransitionPosteriorOther,
#       file = file.path(dir.dropbox, "wolverine/CM/2021/plot24/Figure/CountryTransitionAliveCulledOther.RData"))
# 
# ###c++
# t=2
# habbRtransWolverine <- habbRid <- disaggregate(habitatM$habitat.r,fact=4)
# habbRid[] <-1:ncell(habbRtransWolverine) 
# resultsSXYZ_MF$sims.list$sxy <- UTMToGrid(data.sxy = resultsSXYZ_MF$sims.list$sxy,
#                                             grid.sp = SpatialPoints(coordinates(habbRtransWolverine)))$data.scaled.xy
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
# TransitionSurfaceWolverine[[t]] <- GetTransitionSurface( 
#   resultsSXYZ_MF$sims.list$sxy[,,1,t],
#   resultsSXYZ_MF$sims.list$sxy[,,2,t],
#   resultsSXYZ_MF$sims.list$z[,,t],
#   resultsSXYZ_MF$sims.list$z[,,t+1],
#   as.matrix(rr),
#   stateFrom = 2,
#   stateTo = c(3),
#   ncell = max(rr[]),
#   probs=c(0.025,0.975))
#    
#    TransitionSurfaceWolverine[[t]]$PosteriorTransition <- NULL
#    gc()
#    rrr <- habbRtransWolverine
#    rrr[habbRtransWolverine==1] <- TransitionSurfaceWolverine[[t]]$MeanCell
#    plot(rrr)
#    # if(t <nYears){
#    # plot(dead[dead$Year==years[t],],pch=16,add=T, cex=0.1, col=adjustcolor("black",alpha.f = 0.2))
#    # }
#    sum(rrr[])
# }
# 
# ##-- SAVE OBJECTS
# save( TransitionSurfaceWolverine,
#       habbRtransWolverine,
#       file = file.path(dir.dropbox,"wolverine/CM/2021/plot24/Figure/TransitionAliveCulled.RData"))
# 
# ###LARGE RASTER 
# ###c++
# t=2
# habbRtransWolverine <-habbRid<-  aggregate(disaggregate(habitatM$habitat.r,fact=2),fact=3)
# habbRid[] <-1:ncell(habbRtransWolverine) 
# resultsSXYZ_MF$sims.list$sxy <- UTMToGrid(data.sxy = resultsSXYZ_MF$sims.list$sxy,
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
#   TransitionSurfaceWolverine[[t]] <- GetTransitionSurface( 
#     resultsSXYZ_MF$sims.list$sxy[,,1,t],
#     resultsSXYZ_MF$sims.list$sxy[,,2,t],
#     resultsSXYZ_MF$sims.list$z[,,t],
#     resultsSXYZ_MF$sims.list$z[,,t+1],
#     as.matrix(rr),
#     stateFrom = 2,
#     stateTo = c(3),
#     ncell = max(rr[]),
#     probs=c(0.025,0.975))
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
# save( TransitionSurfaceWolverine,
#       habbRtransWolverine,
#      file = file.path(dir.dropbox,"wolverine/CM/2021/plot24/Figure/TransitionAliveCulled30KM.RData"))
#
# 
# 
##------------------------------------------------------------------------------