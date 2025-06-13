library(rovquantR)
library(nimbleSCR)
library(dplyr)
library(sf)

##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

source("C:/My_documents/RovQuant/Temp/PD/myWorkingDirectories.R")             

data.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2024/Data"
working.dir = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2024/Analysis_report"

species = "bbear"
nburnin = 5
niter = 100
extraction.res = 5000
print.report = TRUE
Rmd.template = NULL
output.dir = NULL
overwrite = FALSE

sex <- c("female","male")
years <- 2015:2024 
n.years <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))


## -----------------------------------------------------------------------------
## ------ 0. BASIC SET-UP ------

if(is.null(working.dir)){working.dir <- getwd()}

##-- Extract date from the last cleaned data file
DATE <- getMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_bear")

##-- Initialize output list
out <- list( SPECIES = "Brown bear",
             engSpecies = "bear",
             DATE = DATE)



## ------ 1. LOAD NECESSARY INPUTS -----

##-- Females
load(list.files(file.path(working.dir, "nimbleInFiles/female"), full.names = T)[1])
nimDataF <- nimData
nimInitsF <- nimInits
y.deadF <- nimData$y.dead

##-- Males
load(list.files(file.path(working.dir, "nimbleInFiles/male"), full.names = T)[1])
nimDataM <- nimData
nimInitsM <- nimInits
y.deadM <- nimData$y.dead
y.deadMF <- rbind(y.deadM, y.deadF)

##-- Habitat
load(file.path( working.dir, "data",
                paste0("Habitat_bear_", DATE, ".RData")))

##-- Detectors
load(file.path( working.dir, "data",
                paste0("Detectors_bear_", DATE, ".RData")))

##-- Habitat Rasters
if(extraction.res <= 1000){
  data(habitatRasterResolution, envir = environment()) 
  extraction.raster <- habitatRasterResolution$'1km'
  extraction.res <- 1000
} else {
  if(extraction.res <= 2000){
    data(habitatRasterResolution, envir = environment()) 
    extraction.raster <- habitatRasterResolution$'2km'
    extraction.res <- 2000
  } else {
    if(extraction.res <= 5000){
      data(habitatRasterResolution, envir = environment()) 
      extraction.raster <- habitatRasterResolution$'5km'
      extraction.res <- 5000
    } else {
      if(extraction.res <= 10000){
        data(habitatRasterResolution, envir = environment()) 
        extraction.raster <- habitatRasterResolution$'10km'
        extraction.res <- 10000
      } else {
        data(habitatRasters, envir = environment()) 
        extraction.raster <- habitatRasters
        extraction.res <- 20000
      }}}}






## ------ II. DATA CHARACTERISTICS -----

myFullData.sp <- readMostRecent( path = file.path(working.dir,"data"),
                                 pattern = "Clean",
                                 extension = ".RData")

##-- MERGE SOME NORWEGIAN COUNTIES
COUNTIES$NAME_1[COUNTIES$NAME_1 %in% c("Sør-Trøndelag",
                                       "Nord-Trøndelag",
                                       "Nordland")] <- "Nord-Trøndelag"
COUNTIES$NAME_1[COUNTIES$NAME_1 %in% c("Troms",
                                       "Finnmark")] <- "Finnmark"
COUNTIES$NAME_1[!COUNTIES$NAME_1 %in% c("Finnmark",
                                        "Troms",
                                        "Nordland", 
                                        "Nord-Trøndelag")] <- "Hedmark"
COUNTIES <- COUNTIES %>%
  group_by(NAME_1) %>%
  summarize()


##-- SIMPLIFY COUNTIES
COUNTIES_s <- sf::st_simplify(sf::st_as_sf(COUNTIES), preserveTopology = T, dTolerance = 500)
COUNTIES_s$index <- c(1,3,2)
COUNTIES_s$Name <- c("NO1","NO3","NO2")

##-- POLYGONS OF COUNTIES IN NORWAY
##-- POLYGONS OF COMMUNES IN NORWAY
COMMUNES <- read_sf(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp")) 
REGIONS <- raster::aggregate(x = as_Spatial(COMMUNES), by = "NAME_1")

##-- LARGE CARNIVORE MANAGEMENT REGIONS
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Troms","Finnmark")] <- "Region 8"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Nordland")] <- "Region 7"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Sør-Trøndelag","Nord-Trøndelag","Møre og Romsdal")] <- "Region 6"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Hedmark")] <- "Region 5"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Akershus","Ãstfold","Oslo")] <- "Region 4"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Oppland")] <- "Region 3"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Vestfold","Telemark","Buskerud","Aust-Agder")] <- "Region 2"
REGIONS$NAME_1[REGIONS$NAME_1 %in% c("Hordaland","Sogn og Fjordane","Rogaland","Vest-Agder")] <- "Region 1"
REGIONS <- raster::aggregate(x = REGIONS, by = "NAME_1")

##-- SIMPLIFY COUNTIES
REGIONS <- st_simplify(st_as_sf(REGIONS), preserveTopology = T, dTolerance = 500)
REGIONS$region <- 1:8
REGIONS$name <- c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6","Region 7", "Region 8")

##-- PREPARE MAPS BACKGROUND
COUNTRIES_BG <- read_sf(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) %>%
  filter(.,
         ISO %in% c("SWE","NOR"),
         area > 80000000)



## -----------------------------------------------------------------------------

## ------ III. REPORT PLOTS & TABLES -----

## ------   0. GENERAL SET-UP -----

##-- load processed density
load(file.path(working.dir, "data", "Density5km.RData")))

##-- load processed MCMC outputs
load(file.path(working.dir, "data", "_MCMC.RData")))

##-- Create 5km raster for plotting
rrNorway <- habitatRasterResolution$`5km`[["Countries"]]
rrNorway[!rrNorway[] %in% 2] <- NA
rrNorway[rrNorway[] %in% 2] <- 1
rrNorway <- crop(rrNorway, myHabitat$habitat.r)

##-- Extract number of individuals detected
isDetected <- rbind(nimDataM$y.alive[ ,1, ],nimDataF$y.alive[ ,1, ]) > 0
n.detected <- apply(isDetected, 2, sum)

##-- Identify individual sex
isFemale <- myResultsSXYZ_MF$sims.list$sex == "F"
isMale <- myResultsSXYZ_MF$sims.list$sex == "M"

##-- Identify individual status
isAvail <- myResultsSXYZ_MF$sims.list$z == 1 
isAlive <- myResultsSXYZ_MF$sims.list$z == 2 
isDead <- myResultsSXYZ_MF$sims.list$z >= 3

##-- Get number of iterations
n.iter <- dim(myResultsSXYZ_MF$sims.list$z)[1]

##-- Plot parameters
diffSex <- 0.2
colSex <- adjustcolor(c("firebrick3", "deepskyblue2", "black"),0.5)
colCause  <- adjustcolor(c("#E69F00","#009E73"), 0.5)



## ------   1. PLOTS -----

## ------     1.1. DENSITY MAPS -----

print("## plotting density maps...")

##-- Names of regions to extract density for
regions.NOR <- c("Region 1","Region 2","Region 3","Region 4",
                 "Region 5","Region 6","Region 7","Region 8")

##-- Load habitat
load(file.path(working.dir, "MODELS", modelNameF, "DATA/Habitat.RData"))

##-- Remove buffer from the habitat
habitat.rWthBuffer <- myHabitat$habitat.rWthBuffer
habitat.rWthBuffer[habitat.rWthBuffer %in% 0] <- NA
searchedPolygon <- rasterToPolygons( habitat.rWthBuffer,
                                     dissolve = T,
                                     function(x) x == 1)

##-- Habitat raster with extent used in the model
habitatPolygon5km <- crop(habitatRasterResolution$`5km`[["Habitat"]],
                          myHabitat$habitat.r)

##-- Create 5km raster for extraction
rrRegions <- habitatRasterResolution$`5km`[["Regions"]]
rrRegions <- mask(rrRegions, myHabitat$habitat.poly)
rrRegions <- crop(rrRegions, myHabitat$habitat.r)
plot(rrRegions)

##-- Load model input
load(file = file.path( working.dir, "MODELS", modelNameF, "NimbleInFiles",
                       paste0(modelNameF, "_input_1.RData")))

##-- Load processed results
load(file.path(working.dir, "RESULTS", paste0(modelName,"_MCMC.RData")))

##-- Get the objects to run the density function
densityInputRegions <- getDensityInput( 
  regions = rrRegions, 
  habitat = habitatPolygon5km,
  s = myResultsSXYZ_MF$sims.list$sxy,
  plot.check = F)

##-- Subset to regions of interest 
regionID <- densityInputRegions$regions.rgmx
row.names(regionID) <- row.names(densityInputRegions$regions.rgmx) 
regionID <- as.matrix(regionID[row.names(regionID) %in% regions.NOR, ])

##-- Select n.iter iterations randomly 
n.iter <- 1000
iter <- sample(1:dim(densityInputRegions$sy)[1], size = n.iter)



## ------       1.1.0. DENSITY PLOT SET-UP ------

##-- Convert densities from 25km2 (5*5raster) to 100km2
UD100km2 <- lapply(spaceUSED[3:12], function(x) x$MeanCell * 4)
AC100km2 <- lapply(ACdensity[3:12], function(x) x$MeanCell * 4)
rrCombined <- rrRegions + rrNorway

##-- Crop density maps to Norway
ACCropped <- UDCropped <- list()
for(t in 1:length(years)){
  ACCropped[[t]] <- densityInputRegions$regions.r
  ACCropped[[t]][] <- NA
  ACCropped[[t]][!is.na(densityInputRegions$regions.r[])] <- AC100km2[[t]]
  ACCropped[[t]][is.na(rrCombined[])] <- NA
  
  UDCropped[[t]] <- densityInputRegions$regions.r
  UDCropped[[t]][] <- NA
  UDCropped[[t]][!is.na(densityInputRegions$regions.r[])] <- UD100km2[[t]]
  UDCropped[[t]][is.na(rrCombined[])] <- NA
  
  proj4string(ACCropped[[t]]) <- proj4string(ACCropped[[t]]) <- st_crs(myHabitat$habitat.poly)
}#t



## ------       1.1.2. AC-Density ------

pdf(file = file.path(working.dir, "figures", paste0("AC_DensityMaps.pdf")),
    width = 12, height = 8)

##-- Set color scale
max <- max(unlist(lapply(ACCropped, function(x) max(x[], na.rm = T))))
cuts <- seq(0, max, length.out = 100) ##-- set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

##-- layout
mx <- rbind(c(1,rep(1:5, each = 2)),
            c(rep(1:5, each = 2), 5))
mx <- rbind(mx, mx + 5)
nf <- layout(mx,
             widths = c(rep(1,ncol(mx))),
             heights = rep(1,2))
# layout.show(nf)
par(mar = c(0,0,0,0))

##-- Plot AC maps
for(t in 1:length(years)){
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
  image(ACCropped[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
  mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)
  
  if(t == n.years){
    segments(x0 = 750000, x1 = 750000,
             y0 = 6750000, y1 = 6750000 + 500000,
             col = grey(0.3), lworking.dir = 4, lend = 2)
    text(700000, 6750000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( UDCropped[[t]],
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           cex.axis = 1.2),
          smallplot = c(0.81, 0.85, 0.2, 0.6),
          legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                             side = 2, font = 1, line = 0.5, cex = 1))
  }#if
  
  ##-- Export rasters
  writeRaster(ACCropped[[t]],
              file.path(working.dirRasters,
                        paste0("AC5kmRaster100km2", years[t],".tif")),
              overwrite = TRUE)
}#t
dev.off()


ACCropped[[t]]



## ------       1.1.3. UD-Density ------

pdf(file = file.path(working.dir, "figures", paste0("UD_DensityMaps.pdf")),
    width = 12, height = 8)

##-- Set color scale
max <- max(unlist(lapply(UDCropped, function(x) max(x[], na.rm = T))))
cuts <- seq(0, max, length.out = 100) ##-- set breaks
colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
col <- colfunc(100)

##-- layout
mx <- rbind(c(1,rep(1:5, each = 2)),
            c(rep(1:5, each = 2), 5))
mx <- rbind(mx, mx + 5)
nf <- layout(mx,
             widths = c(rep(1,ncol(mx))),
             heights = rep(1,2))
#layout.show(nf)
par(mar = c(0,0,0,0))

##-- Plot UD maps
for(t in 1:length(years)){
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  image(UDCropped[[t]], add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = grey(0.4), col = NA, add=TRUE)
  mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)
  
  if(t == n.years){
    segments(x0 = 750000, x1 = 750000,
             y0 = 6750000, y1 = 6750000 + 500000,
             col = grey(0.3), lworking.dir = 4, lend = 2)
    text(700000, 6750000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( UDCropped[[t]],
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           cex.axis = 1.2),
          smallplot = c(0.81, 0.85, 0.2, 0.6),
          legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                             side = 2, font = 1, line = 0.5, cex = 1))
  }#if
  
  ##-- Export rasters
  writeRaster(UDCropped[[t]],
              file.path(working.dirRasters,
                        paste0("UD5kmRaster100km2", years[t],".tif")),
              overwrite = TRUE)
}#t
dev.off()



## ------       1.1.5. LAST YEAR ------

pdf(file = file.path(working.dir, "figures", paste0("UD_DensityMaps_LastYear.pdf")), 
    width = 6, height = 8)

t <- length(years)

par(mar = c(0,0,0,0))
plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
image( UDCropped[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
plot(COUNTRIESsimpFig[1, ], border = grey(0.4), col = NA, add = TRUE)
mtext(text = years[t], side = 1, -25, adj = 0.25, cex = 2, font = 2)
segments(x0 = 760000, x1 = 760000,
         y0 = 6650000, y1 = 6650000 + 500000,
         col = grey(0.3), lworking.dir = 4, lend = 2)  
text(720000, 6650000+500000/2, labels = "500 km", srt = 90, cex = 1.2)
plot( UDCropped[[t]],
      legend.only = T,
      breaks = cuts,
      col = col,
      legend.width = 2,
      axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       cex.axis = 1.2),
      smallplot = c(0.75, 0.78, 0.2, 0.4),
      legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                         side = 2, font = 1, line = 0, cex = 1.2))
dev.off()



## ------       1.1.6. LAST YEAR SUMMARY ------

pdf(file = file.path(working.dir, "figures", paste0("UD_DensityMaps_LastYearSummary.pdf")), 
    width = 8, height = 8)

t <- length(years)

par(mar = c(0,0,0,0))
plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
image( UDCropped[[t]], add = T, breaks = c(cuts, max(cuts)+1000), col = col, legend = F)
plot(COUNTRIESsimpFig[1, ], border = grey(0.4), col = NA, add = TRUE)
mtext(text = years[t], side = 1, -25, adj = 0.25, cex = 3, font = 2)
segments(x0 = 800000, x1 = 800000,
         y0 = 6650000, y1 = 6650000 + 500000,
         col = grey(0.3), lworking.dir = 4, lend = 2)  
text(760000, 6650000+500000/2, labels = "500 km", srt = 90, cex = 1.2)
plot( UDCropped[[t]],
      legend.only = T,
      breaks = cuts,
      col = col,
      legend.width = 2,
      axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       cex.axis = 1.2),
      smallplot = c(0.75, 0.78, 0.2, 0.4),
      legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                         side = 2, font = 1, line = 0, cex = 1.2))
dev.off()




## ------       1.1.7. LAST YEAR SUMMARY NO ------
pdf(file = file.path(working.dir, "figures", paste0("UD_DensityMaps_LastYearSummary_NO.pdf")), 
    width = 8, height = 8)

t <- length(years)

par(mar = c(0,0,0,0))
plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
image( UDCropped[[t]], add = T, breaks = c(cuts, max(cuts)+1000), col = col, legend = F)
plot(COUNTRIESsimpFig[1, ], border=grey(0.1), col = NA, add = TRUE)
mtext(text = years[t], side = 1, -25, adj = 0.25, cex = 3, font = 2)
segments(x0 = 800000, x1 = 800000,
         y0 = 6650000, y1 = 6650000 + 500000,
         col = grey(0.3), lworking.dir = 4, lend = 2)  
text(760000, 6650000+500000/2, labels = "500 km", srt = 90, cex = 1.2)
plot( UDCropped[[t]],
      legend.only = T,
      breaks = cuts,
      col = col,
      legend.width = 2,
      axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                       cex.axis = 1.2),
      smallplot = c(0.75, 0.78, 0.2, 0.4),
      legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                         side = 2, font = 1, line = 0, cex = 1.2))
dev.off()



## -----------------------------------------------------------------------------
## ------     1.2. ABUNDANCE ------

print("## plotting abundance...")

## ------       1.2.2. BARS ------

pdf(file = file.path(working.dir, "figures", paste0("N_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, n.years + 0.5), ylim = c(0,180),
     xlab = "", ylab = paste("Number of bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
axis(2, at = seq(0,180,20), labels = seq(0,180,20), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,180, by = 10), lty = 2, col = "gray60")

for(t in 1:n.years){
  ##-- FEMALES 
  plotQuantiles(x = colSums(ACdensityF[[t+2]]$PosteriorAllRegions),
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = colSums(ACdensityM[[t+2]]$PosteriorAllRegions),
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  ##-- TOTAL 
  plotQuantiles(x = colSums(ACdensity[[t+2]]$PosteriorAllRegions),
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- ADD NUMBER OF INDIVIDUALS DETECTED
  xx <- c(t-0.25,t+0.25,t+0.25,t-0.25)
  yy <- c(n.detected[t+2]-1,n.detected[t+2]-1,n.detected[t+2]+1,n.detected[t+2]+1)
  polygon(xx, yy, border = NA, col = "goldenrod1")
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(5.1,6.6,7.8,9)
yy <- c(5,5,5,5)
labs <- c("Females", "Males", "Total", "Detected")
polygon(x = c(4.6,10.6,10.6,4.6),
        y = c(-2.5,-2.5,12.5,12.5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col = colSex)
points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

polygon(x = c(xx[4]-0.3,xx[4]+0.1,xx[4]+0.1,xx[4]-0.3),
        y = c(yy[4]-1,yy[4]-1,yy[4]+1,yy[4]+1),
        col = "goldenrod1",
        border = F)

dev.off()



## ------       1.2.3. BARS per REGION ------
pdf(file = file.path(working.dir, "figures", paste0("N_bars_regions.pdf")),
    width = 8, height = 15)

nf <- layout(rbind(c(1,2),
                   c(3,4),
                   c(5,6),
                   c(7,8),
                   c(9,10)),
             widths = c(1,0.5),
             heights = 1)

for(cc in c(8,7,6,5,3)){
  par(mar = c(2,5,2,0), tck = 0, mgp = c(2,0.2,0))
  plot(-1000,
       xlim = c(0.5, n.years + 0.5), ylim = c(0,60),
       xlab = "", ylab = paste("Number of bears"),
       xaxt = "n", axes = F, cex.lab = 1.6)
  axis(1, at = c(1:n.years), labels = years, cex.axis = 1.15)
  axis(2, at = seq(0,60,10), labels = seq(0,60,10),
       cex.axis = 1.2, las = 1)
  abline(v = (0:n.years)+0.5, lty = 2)
  abline(h = seq(0,60, by = 10), lty = 2, col = "gray60")
  mtext(text = REGIONS$name[REGIONS$region == cc],
        side = 3, font = 2, line = 0.5)
  
  for(t in 1:n.years){
    ##-- TOTAL 
    plotQuantiles(x = ACdensity[[t+2]]$PosteriorRegions[cc-1, ],
                  at = t,
                  width = 0.4,
                  col = colSex[3])
    
    ##-- FEMALES 
    plotQuantiles(x = ACdensityF[[t+2]]$PosteriorRegions[cc-1, ],
                  at = t - diffSex,
                  width = 0.18,
                  col = colSex[1])
    
    ##-- MALES 
    plotQuantiles(x = ACdensityM[[t+2]]$PosteriorRegions[cc-1, ],
                  at = t + diffSex,
                  width = 0.18,
                  col = colSex[2])
    
  }#t
  ##-- legend
  if(cc %in% 3){
    par(xpd = TRUE)
    xx <- c(1.1,3.2,5)
    yy <- c(50,50,50)
    labs <- c("Females", "Males", "Total")
    polygon(x = c(0.6,6.4,6.4,0.6),
            y = c(50-7,50-7,50+7,50+7),
            col = adjustcolor("white", alpha.f = 0.9),
            border = "gray60")
    
    points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col = colSex)
    points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col = colSex)
    text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  }
  box()
  
  ##-- Map insert
  par(mar = c(5,0,4,5))
  plot(st_geometry(REGIONS), border = grey(0.5), col = grey(0.5), lworking.dir = 0.1)
  plot(st_geometry(REGIONS[REGIONS$region == cc, ]),
       add = T, col = adjustcolor("red",0.5), border = "red")
}#c
dev.off()




## -----------------------------------------------------------------------------
## ------     1.3. VITAL RATES ------

print("## plotting vital rates")

n.iter <- dim(myResultsSXYZ_MF$sims.list$z)[1]

##-- Calculate mortality from estimated hazard rates (mhH and mhW) if necessary
mhH1 <- exp(myResults_F$sims.list$mhH[ ,3:11])
mhW1 <- exp(myResults_F$sims.list$mhW[ ,3:11])
myResults_F$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
myResults_F$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
myResults_F$sims.list$phi <- 1 - myResults_F$sims.list$h - myResults_F$sims.list$w

mhH1 <- exp(myResults_M$sims.list$mhH[ ,3:11])
mhW1 <- exp(myResults_M$sims.list$mhW[ ,3:11])
myResults_M$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
myResults_M$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
myResults_M$sims.list$phi <- 1 - myResults_M$sims.list$h - myResults_M$sims.list$w



## ------       1.3.1. SURVIVAL ------

pdf(file = file.path(working.dir, "figures", paste0("Survival.pdf")),
    width = 10, height = 6)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
             widths = c(0.05,1,0.30),
             heights = c(0.15,1))

par(mar = c(10,4.5,0.5,0.5), tck = 0, xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
plot(10, xlim = c(0.5, n.years-1 + 0.5), ylim = c(0,1),
     type ="n", xaxt = "n", xlab = "", ylab = "Survival")
axis(1, c(1:n.years),
     labels = paste(years, years+1, sep = "\n to \n"),
     padj = 0.6, cex = 0.9)
mtext("Years",side = 1,line = 5)
abline(v = 1:(n.years-1) + 0.5, lty = 2)
axis(2, tck = -0.02)

for(t in 1:(n.years-1)){
  ##-- FEMALES
  plotQuantiles( x = myResults_F$sims.list$phi[ ,t],
                 at = t - diffSex,
                 col = colSex[1])
  
  ##-- MALES
  plotQuantiles( x = myResults_M$sims.list$phi[ ,t],
                 at = t + diffSex,
                 col = colSex[2])
}#t

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
points(c(4,4), c(4,3), pch = 15, cex = 5.5, col = colSex)
points(c(4,4), c(4,3), pch = 15, cex = 3, col = colSex)
text(c(5.3,5.3), c(4,3),  c("Females", "Males"), cex = 2, pos = 4)
dev.off()




## ------       1.3.2. MORTALITY ------

pdf(file = file.path(working.dir, "figures", paste0("Mortality.pdf")),
    width = 10, height = 9)

nf <- layout(rbind(c(3,7,6),
                   c(4,1,6),
                   c(5,2,6)),
             widths = c(0.15,1,0.35),
             heights = c(0.15,1,1))

for(s in 1:2){
  if(s == 1){results <- myResults_F$sims.list} else {results <- myResults_M$sims.list}
  par(mar=c(8,4.5,0.5,1),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
  plot(10, xlim = c(0.5, n.years-1+0.5), ylim = c(0,0.5),
       type = "n", xaxt = "n", xlab = "", ylab = "Mortality")
  axis(1, c(1:n.years),
       labels = paste(years, years+1, sep = "\n to \n"),
       padj = 0.6, cex = 0.9)
  mtext("Years",side = 1,line = 5)
  axis(2, tck = -0.02)
  abline(v = 1:(n.years-1)+0.5, lty = 2)
  myDev <- c(-0.15,+0.15)
  
  for(t in 1:(n.years-1)){
    plotQuantiles(x = results$h[ ,t],
                  at = t - diffSex,
                  col = colCause[1])
    plotQuantiles(x = results$w[ ,t],
                  at = t + diffSex,
                  col = colCause[2])
  }#t
}#s

##-- LABELS
par(mar = c(0,0,0,0))
plot(1, axes = FALSE, ylim = c(-1,1), xlim = c(-1,1), type = "n")
my.labels=c("Female","Male")
for(i in my.labels){
  par(mar=c(0.5,0.5,0.5,0.5))
  plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
  text(0,0.25,labels=i,srt=90,cex=3,font=2)
}

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = F)
points(c(2,2),c(3.5,2.7),pch=15,cex=3.5,col = colCause)
points(c(2,2),c(3.5,2.7),pch=15,cex=1.5,col = colCause)
text(c(3.5,3.5),c(3.5,2.7),c("Legal\nculling","Other\nmortality"),cex=2,pos=4)
dev.off()



## ------       1.3.3. RECRUITMENT ------
## ------         1.3.3.1. PLOT PER CAPITA RECRUITMENT -----

pdf(file = file.path( working.dir, "figures",
                      paste0("PerCapita.pdf")), width = 10, height = 8)

##-- PER CAPITA RECRUITMENT 
par(mar = c(4.5,4.5,1,1), xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
plot(10,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,0.5),
     type ="n", xaxt = "n", xlab = "Years", ylab = "Per-capita recruitment")
axis(2, tck = -0.02)
axis(1, 1:(n.years-1),
     labels = paste(years[1:(n.years-1)], years[2:n.years], sep = " to "), 
     cex.axis = 0.8)
abline(v = 1:(n.years-1)+0.5,lty=2)

for(t in 1:(n.years-1)){
  ##-- Number of recruits at t ==> ids with state 1 at t and 2 at t+1 
  n.recruit <- apply(myResultsSXYZ_MF$sims.list$z[ , ,c(t+2,t+3)], 1, function(x) sum(x[ ,1]%in%1 & x[ ,2]%in%2))
  ##-- Number of reproducing ids at t ==> ids with state 2 at t
  alivetminus1 <- apply(myResultsSXYZ_MF$sims.list$z[ , ,t+2], 1, function(x)sum(x %in% 2))
  ##-- Plot quantiles
  plotQuantiles(x = n.recruit/alivetminus1,
                at = t, 
                width = 0.3,
                col = colCause[1])
}#t
dev.off()



## ------         1.3.3.2. PLOT NUMBER OF RECRUITS ----
N_recruit <- N_recruit_F <- N_recruit_M <- matrix(NA, n.iter, n.years-1)
for(t in 1:(n.years-1)){
  for(iter in 1:n.iter){
    N_recruit_F[iter,t] <- sum(isAvail[iter, ,t+2] & isAlive[iter, ,t+3] & isFemale)
    N_recruit_M[iter,t] <- sum(isAvail[iter, ,t+2] & isAlive[iter, ,t+3] & isMale)
    N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]
  }#iter
  print(t)
}#t


pdf(file = file.path(working.dir, "figures", paste0("NumRecruitTotal.pdf")),
    width = 10, height = 6)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
             widths = c(0.05,1,0.30),
             heights = c(0.15,1))

par(mar = c(4.5,4.5,1,1), xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)

plot(10,
     xlim = c(0.5, n.years-0.5), ylim = c(0,100),
     type ="n", xaxt = "n",
     xlab = "Years", ylab = "Number of recruits")
axis(2, tck = -0.02)
axis(1, 1:(n.years-1),
     labels = years[2:n.years], cex.axis=0.8)
abline(v = 1:(n.years-1)+0.5,lty=2)
for(t in 1:(n.years-1)){
  plotQuantiles(x = N_recruit[ ,t],
                at = t,
                width = 0.3,
                col = colSex[3])
  
  plotQuantiles(x = N_recruit_F[ ,t],
                at = t - diffSex,
                col = colSex[1])
  
  plotQuantiles(x = N_recruit_M[ ,t],
                at = t + diffSex,
                col = colSex[2])
}#t

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
points(c(4,4,4), c(4,3,2), pch = 15, cex = 5.5, col = colSex)
points(c(4,4,4), c(4,3,2), pch = 15, cex = 3, col = colSex)
text(c(5.3,5.3,5.3), c(4,3,2),
     c("Females", "Males", "Total"), cex = 2, pos = 4)

dev.off()



## ------     1.4. POPULATION FLUXES IN/OUT NORWAY -----

norRaster <- habitatRasterResolution$`5km`[["Countries"]]
norRaster[norRaster[] != 2] <- NA

##-- Calculate number of individuals 
N_surv <- N_surv_F <- N_surv_M <- matrix(NA,n.iter,n.years-1)
N_recruit <- N_recruit_F <- N_recruit_M <- matrix(NA,n.iter,n.years-1)
N_immig <- N_immig_F <- N_immig_M <- matrix(NA,n.iter,n.years-1)
N_emig <- N_emig_F <- N_emig_M <- matrix(NA,n.iter,n.years-1)

for(t in 1:(n.years-1)){
  for(iter in 1:n.iter){
    isOut <- is.na(norRaster[cellFromXY(norRaster,myResultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+3])])
    wasOut <- is.na(norRaster[cellFromXY(norRaster,myResultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+2])])
    
    N_recruit_F[iter,t] <- sum(isAvail[iter, ,t+2] & isAlive[iter, ,t+3] & (!isOut) & isFemale)
    N_recruit_M[iter,t] <- sum(isAvail[iter, ,t+2] & isAlive[iter, ,t+3] & (!isOut) & isMale)
    N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]
    
    N_immig_F[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & wasOut & (!isOut) & isFemale)
    N_immig_M[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & wasOut & (!isOut) & isMale)
    N_immig[iter,t] <- N_immig_F[iter,t] + N_immig_M[iter,t]
    
    N_emig_F[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & (!wasOut) & isOut & isFemale)
    N_emig_M[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & (!wasOut) & isOut & isMale)
    N_emig[iter,t] <- N_emig_F[iter,t] + N_emig_M[iter,t]
    
    N_surv_F[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & (!wasOut) & (!isOut) & isFemale)
    N_surv_M[iter,t] <- sum(isAlive[iter, ,t+2] & isAlive[iter, ,t+3] & (!wasOut) & (!isOut) & isMale)
    N_surv[iter,t] <- N_surv_F[iter,t] + N_surv_M[iter,t]
  }#iter
  print(t)
}#t



## ------       1.4.1. PLOT RECRUITS -----
pdf(file = file.path(working.dir, "figures", paste0("NumRecruits_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of bears recruited"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")

for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_recruit[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_recruit_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_recruit_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(1,2.5,4)
yy <- c(55,55,55)
labs <- c("Females", "Males", "Total")
polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
        y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx, y = yy,  pch = 15, cex = 3.5, col = colSex)
points(x = xx, y = yy,  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

dev.off()






## ------       1.4.2. PLOT SURVIVORS -----
pdf(file = file.path(working.dir, "figures", paste0("NumSurvival_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,120),
     xlab = "", ylab = paste("Number of surviving bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,120,20), labels = seq(0,120,20), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,120, by = 10), lty = 2, col = "gray60")

for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_surv[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_surv_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_surv_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(1,2.5,4)
yy <- c(55,55,55)
labs <- c("Females", "Males", "Total")
polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
        y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx, y = yy,  pch = 15, cex = 3.5, col = colSex)
points(x = xx, y = yy,  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

dev.off()







## ------       1.4.3. PLOT IMMIGRANTS -----
pdf(file = file.path(working.dir, "figures", paste0("NumImmigrants_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of immigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")

for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_immig[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_immig_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_immig_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(1,2.5,4)
yy <- c(55,55,55)
labs <- c("Females", "Males", "Total")
polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
        y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx, y = yy,  pch = 15, cex = 3.5, col = colSex)
points(x = xx, y = yy,  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

dev.off()




## ------       1.4.4. PLOT EMIGRANTS -----
pdf(file = file.path(working.dir, "figures", paste0("NumEmigrants_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of emigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")

for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_emig[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_emig_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_emig_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(1,2.5,4)
yy <- c(55,55,55)
labs <- c("Females", "Males", "Total")
polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
        y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx, y = yy,  pch = 15, cex = 3.5, col = colSex)
points(x = xx, y = yy,  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

dev.off()







## ------       1.4.5. PLOT ALL -----
pdf(file = file.path(working.dir, "figures", paste0("NumFluxes_bars.pdf")),
    width = 12, height = 8.5)
par(mfrow = c(2,2))
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of bears recruited"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")
for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_recruit[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_recruit_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_recruit_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()


plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,120),
     xlab = "", ylab = paste("Number of surviving bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,120,20), labels = seq(0,120,20), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,120, by = 10), lty = 2, col = "gray60")
for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_surv[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_surv_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_surv_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()


plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of immigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")
for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_immig[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_immig_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_immig_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
}#t
box()


plot(-1000,
     xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,60),
     xlab = "", ylab = paste("Number of emigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
axis(2, at = seq(0,60,10), labels = seq(0,60,10), cex.axis = 1.6)
abline(v = (0:n.years)+0.5, lty = 2)
abline(h = seq(0,60, by = 5), lty = 2, col = "gray60")
for(t in 1:(n.years-1)){
  ##-- TOTAL 
  plotQuantiles(x = N_emig[ ,t],
                at = t,
                width = 0.4,
                col = colSex[3])
  
  ##-- FEMALES 
  plotQuantiles(x = N_emig_F[ ,t],
                at = t - diffSex,
                width = 0.18,
                col = colSex[1])
  
  ##-- MALES 
  plotQuantiles(x = N_emig_M[ ,t],
                at = t + diffSex,
                width = 0.18,
                col = colSex[2])
  
}#t
box()

##-- legend
par(xpd = TRUE)
xx <- c(1,2.5,4)
yy <- c(55,55,55)
labs <- c("Females", "Males", "Total")
polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
        y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
        col = adjustcolor("white", alpha.f = 0.9),
        border = "gray60")

points(x = xx, y = yy,  pch = 15, cex = 3.5, col = colSex)
points(x = xx, y = yy,  pch = 15, cex = 1.5, col = colSex)
text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)

dev.off()









## ------     1.5. p0 ------

print("## plotting p0...")

## ------       1.5.1. p0 bars ------

pdf(file = file.path(working.dir, "figures", paste0("p0_mod.pdf")),
    width = 8, height = 12)

nf <- layout(rbind(c(1,2),
                   c(3,4),
                   c(5,6)),
             widths = c(1,0.5),
             heights = 1)

for(c in 1:nrow(COUNTIES_s)){
  par(mar=c(4,4,1,1), tck=0)
  
  plot(10, xlim = c(0.5, n.years+0.5), ylim = c(0,0.01), type ="n", xaxt="n",
       xlab = "Years", ylab = "Baseline detection probability") 
  
  axis(1, 1:n.years, labels = years)
  axis(2, tck = -0.02)
  abline(v = 1:(n.years-1) + 0.5, lty = 2)
  
  for(t in 1:n.years){
    plotQuantiles( x = myResults_F$sims.list$p0[ ,COUNTIES_s$index == c,t+2],
                   at = t - diffSex,
                   col = colSex[1])
    
    plotQuantiles( x = myResults_M$sims.list$p0[ ,COUNTIES_s$index == c,t+2],
                   at = t + diffSex,
                   col = colSex[2])
  }#t
  
  ##-- LEGEND
  if(c == 1){
    polygon(x = c(0.5,4.2,4.2,0.5),
            y = c(0.006,0.006,0.01,0.01),
            col = adjustcolor("white", alpha.f = 0.9),
            border = NA)
    points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 5.5, col = colSex)
    points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 3, col = colSex)
    text(c(1.2,1.2),c(0.0077,0.0092),  c("Females", "Males"), cex = 2, pos = 4)
  }
  
  
  par(mar = c(0,0,0,0))
  plot(st_geometry(COUNTIES_s), border = grey(0.5), col = grey(0.5), lworking.dir = 0.1)
  plot(st_geometry(COUNTIES_s[COUNTIES_s$index == c, ]),
       add = T, col = adjustcolor("red",0.5), border = "red")
  text(as_Spatial(COUNTIES_s[COUNTIES_s$index == c, ]),
       labels = COUNTIES_s$Name[COUNTIES_s$index == c],
       col = "white")
}#c
dev.off()




## ------       1.5.2. p0 maps ------

pdf(file = file.path(working.dir, "figures", paste0("p0_map_mod.pdf")),
    width = 10, height = 6)
for(t in 1:n.years){
  par(mfrow = c(1,2))
  
  ##-- FEMALE
  myDetectors$main.detector.sp$p0_F <- 
    ilogit(logit(myResults_F$mean$p0[nimConstants$county,t+2]) +
             myResults_F$mean$betaDet[1] * nimDataF$detCovs[ ,1,t+2] +
             myResults_F$mean$betaDet[2] * nimDataF$detCovs[ ,2,t+2])
  
  p0_F.R <- rasterFromXYZ(cbind(myDetectors$main.detector.sp$main.cell.x,
                                myDetectors$main.detector.sp$main.cell.y,
                                myDetectors$main.detector.sp$p0_F))
  
  plot(p0_F.R,
       main =  paste0("Females ", years[t]),
       legend.args = list(text = 'p0',
                          side = 4, font = 2, line = 2.5, cex = 0.8))
  
  ##-- MALE
  myDetectors$main.detector.sp$p0_M <- 
    ilogit(logit(myResults_M$mean$p0[nimConstants$county,t+2]) +
             myResults_M$mean$betaDet[1] * nimDataM$detCovs[ ,1,t+2] +
             myResults_M$mean$betaDet[2] * nimDataM$detCovs[ ,2,t+2])
  
  p0_M.R <- rasterFromXYZ(cbind(myDetectors$main.detector.sp$main.cell.x,
                                myDetectors$main.detector.sp$main.cell.y,
                                myDetectors$main.detector.sp$p0_M))
  plot(p0_M.R,
       main = paste0("Males ", years[t]),
       legend.args = list(text = 'p0',
                          side = 4, font = 2, line = 2.5, cex = 0.8))
}#t 
dev.off()




## ------       1.5.3. p0 betas ------

pdf(file = file.path(working.dir, "figures", paste0("p0_beta.pdf")),
    width = 8, height = 4)

nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
             widths = c(0.05,1,0.30),
             heights = c(0.15,1))

##-- PLOT BETAS
par(mar = c(5,4.5,0.5,0.5), tck = 0, xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
plot( 10, xlim = c(0.5, 2.5), ylim = c(-5,5),
      type = "n", xaxt = "n", xlab = "Years", ylab = "beta")
axis(1, c(1,2), labels =  c("distance \nto roads","presence of \nother obs."))
axis(2, tck = -0.02)
abline(v = 1.5, lty = 2)
abline(h = 0, lty = 1)
for(b in 1:2){
  plotQuantiles( myResults_F$sims.list$betaDet[ ,b],
                 at = b - diffSex,
                 col = colSex[1])
  plotQuantiles( myResults_M$sims.list$betaDet[ ,b],
                 at = b + diffSex,
                 col = colSex[2])
}#b

##-- LEGEND
par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
points(c(4,4), c(4,3), pch = 15, cex = 5.5, col = colSex)
points(c(4,4), c(4,3), pch = 15, cex = 3, col = colSex)
text(c(5.3,5.3), c(4,3),  c("Females", "Males"), cex = 3, pos = 4)
dev.off()




## ------     1.6. DETECTABILITY ------

## ------       1.6.1. SET-UP ------

##-- Create 5km raster for extraction
regions.r <- habitatRasterResolution$`5km`[["Regions"]]
regions.r <- crop(regions.r, myHabitat$habitat.r)

##-- Set all habitat cells outside Norway to NA
regions.r[regions.r[] < 23] <- NA

##-- Identify habitat cells
isHabitat <- which(!is.na(regions.r[]))

##-- Get habitat cell coordinates
regions.xy <- coordinates(regions.r)[isHabitat, ]

##-- Get region name for each habitat cell (format to run the C++ function)
regionsNames <- sort(unique(na.omit(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regionsNames, function(x)regions.r[] == x))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- factorValues(regions.r, regionsNames)[,1]
regions.rgmx <- regions.rgmx[ ,isHabitat]

##-- Get detectors coordinates (original and scaled)
detectors.xy <- as.matrix(st_coordinates(myDetectors$main.detector.sp))

##-- Select n.iter iterations randomly
iter <- sample(1:dim(myResultsSXYZ_MF$sims.list$sigma)[1], size = 250)

##-- Load detectability calculation
load(file.path(working.dir, "data/Detectability5km.RData"))



## ------       1.6.2. PLOT DETECTABILITY MAPS -----

pdf(file = file.path(working.dir, "figures/Detectability_maps.pdf")),
    width = 10, height = 6)

##-- Set color scale
max <- max(c(unlist(lapply( DetectabilityRegionsM[[t]]$MeanCell,
                            function(x) max(x[], na.rm = T))),
             unlist(lapply( DetectabilityRegionsF[[t]]$MeanCell,
                            function(x) max(x[], na.rm = T)))))
cuts <- seq(0, max, length.out = 100) ##-- set breaks
col <- rev(terrain.colors(100))

##-- layout
mx <- rbind(c(1,rep(1:5, each = 2)),
            c(rep(1:5, each = 2), 5))
mx <- rbind(mx, mx + 5)
nf <- layout(mx,
             widths = c(rep(1,ncol(mx))),
             heights = rep(1,2))
par(mar = c(0,0,0,0))

##-- FEMALES
for(t in 1:length(years)){
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  detectab.r <- regions.r
  detectab.r[isHabitat] <- DetectabilityRegionsF[[t]]$MeanCell
  image(detectab.r, add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])),
       border = grey(0.4), col = NA, add = TRUE)
  mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)
  
  if(t == n.years){
    segments(x0 = 830000, x1 = 830000,
             y0 = 6730000, y1 = 6730000 + 500000,
             col = grey(0.3), lworking.dir = 4, lend = 2)
    text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( detectab.r,
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           cex.axis = 1.2),
          smallplot = c(0.95, 1.00, 0.2, 0.6),
          legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                             side = 2, font = 1, line = 1, cex = 1))
  }#if
}#t

##-- MALES
for(t in 1:length(years)){
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  detectab.r <- regions.r
  detectab.r[isHabitat] <- DetectabilityRegionsM[[t]]$MeanCell
  image(detectab.r, add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])),
       border = grey(0.4), col = NA, add = TRUE)
  mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)
  
  if(t == n.years){
    segments(x0 = 830000, x1 = 830000,
             y0 = 6730000, y1 = 6730000 + 500000,
             col = grey(0.3), lworking.dir = 4, lend = 2)
    text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( detectab.r,
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           cex.axis = 1.2),
          smallplot = c(0.95, 1.00, 0.2, 0.6),
          legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
                             side = 2, font = 1, line = 1, cex = 1))
  }#if
}#t
dev.off()



## ------     1.7. SEX-RATIO ------

## ------       1.7.1. SEX-RATIO BARS ------

pdf(file = file.path(working.dir, "figures/SexRatio_bars.pdf")),
    width = 12, height = 8.5)
par(mar = c(5,5,1,1))
plot(-1000,
     xlim = c(0.5, n.years + 0.5), ylim = c(0.2,0.8),
     xlab = "", ylab = paste("% Female bears"), xaxt = "n", axes = F, cex.lab = 1.6)
axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
axis(2, at = seq(0.2,0.8,0.1), labels = seq(0.2,0.8,0.1), cex.axis = 1.6)
# abline(v = (0:n.years)+0.5, lty = 2)
abline(h = 0.5, lty = 2, col = "gray60")
for(t in 1:n.years){
  plotQuantiles(x = colSums(ACdensityF[[t+2]]$PosteriorAllRegions)/
                  colSums(ACdensity[[t+2]]$PosteriorAllRegions),
                at = t - diffSex,
                width = 0.18,
                col = colSex[3])
}#t
#box()
dev.off()



## ------       1.7.2. SEX-RATIO DISTANCE ------
###--- TBD



## ------       1.7.3. SEX-RATIO MAPS ------
##-- Convert densities from 25km2 (5*5 raster) to 100km2
SexRatio <- list()
for(t in 1:length(years)){
  SexRatio[[t]] <- ACdensityF[[t+2]]$MeanCell/ACdensity[[t+2]]$MeanCell
}

##-- Crop density maps to Norway
rrCombined <- rrRegions + rrNorway
SexRatioMap <- list()
for(t in 1:length(years)){
  SexRatioMap[[t]] <- densityInputRegions$regions.r
  SexRatioMap[[t]][] <- NA
  SexRatioMap[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
  SexRatioMap[[t]][is.na(rrCombined[])] <- NA
  
  crs(SexRatioMap[[t]]) <- st_crs(myHabitat$habitat.poly)
}#t


pdf(file = file.path(working.dir, "figures/SexRatio_Maps.pdf")), 
    width = 12, height = 8)

##-- Set color scale
max <- 1
cuts <- seq(0, max, length.out = 100) ##-- set breaks
colfunc <- colorRampPalette(c(colSex[2],"white", colSex[1]))
col <- colfunc(100)

##-- layout
mx <- rbind(c(1,rep(1:5, each = 2)),
            c(rep(1:5, each = 2), 5))
mx <- rbind(mx, mx + 5)
nf <- layout(mx,
             widths = c(rep(1,ncol(mx))),
             heights = rep(1,2))
#layout.show(nf)
par(mar = c(0,0,0,0))

##-- Plot Sex-ratio maps
for(t in 1:length(years)){
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
  image(SexRatioMap[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
  mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)
  
  if(t == n.years){
    segments(x0 = 830000, x1 = 830000,
             y0 = 6730000, y1 = 6730000 + 500000,
             col = grey(0.3), lworking.dir = 4, lend = 2)
    text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( SexRatioMap[[t]],
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                           cex.axis = 1.2),
          smallplot = c(0.95, 1.00, 0.2, 0.6),
          legend.args = list(text = "% Female bears",
                             side = 2, font = 1, line = 1, cex = 1))
  }#if
}#t
dev.off()




## ------       1.7.4. SEX-RATIO MAPS (SMOOTHED) ------

##-- Set smoothing factor
smoothingFactor <- c(3,7,11,21)
for(SF in smoothingFactor){
  
  SexRatioMap.smooth <- list()
  for(t in 1:length(years)){
    SexRatioMap.smooth[[t]] <- densityInputRegions$regions.r
    SexRatioMap.smooth[[t]][] <- NA
    SexRatioMap.smooth[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
    SexRatioMap.smooth[[t]] <- terra::focal(terra::rast(SexRatioMap[[t]]), SF, "mean", na.rm=TRUE)
    SexRatioMap.smooth[[t]][is.na(rrCombined[])] <- NA
  }#t
  
  
  pdf(file = file.path(working.dir, "figures", paste0("SexRatio_Maps_smooth_",SF,".pdf")), 
      width = 12, height = 8)
  
  ##-- Set color scale
  max <- 1
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  colfunc <- colorRampPalette(c(inferno(100)[1], inferno(100)[25], inferno(100)[50], inferno(100)[75], inferno(100)[100]))
  col <- colfunc(100)
  
  
  ##-- layout
  mx <- rbind(c(1,rep(1:5, each = 2)),
              c(rep(1:5, each = 2), 5))
  mx <- rbind(mx, mx + 5)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)
  par(mar = c(0,0,0,0))
  
  ##-- Plot AC maps
  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
    image( SexRatioMap.smooth[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
    mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)
    
    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lworking.dir = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
      plot( raster(SexRatioMap.smooth[[t]]),
            legend.only = T,
            breaks = cuts,
            col = col,
            legend.width = 2,
            axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                             labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                             cex.axis = 1.2),
            smallplot = c(0.82, 0.84, 0.2, 0.6),
            legend.args = list(text = "% Female bears",
                               side = 2, font = 1, line = 0.2, cex = 1))
    }#if
  }#t
  dev.off()
}#SF



## -----------------------------------------------------------------------------
## ------   2. TABLES -----
gc()

## ------     2.1. ABUNDANCE -----

## ------       2.1.1. ALL YEARS, BOTH SEX -----

idcounty <- row.names(ACdensity[[1]]$summary)

##-- Remove Finland, Norway, Russia, Sweden 
idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
idcounty <- unique(idcounty)

##-- Get names of Norwegian carnivore regions
idcountyNOR <- idcounty[grep("Region",idcounty)]
idcountyTable <- c(idcountyNOR, "Total")

##-- Create table to store N estimates (and CI)
NCarRegionEstimates <- NCarRegionMean <- matrix("", ncol = n.years, nrow = length(idcountyTable))
row.names(NCarRegionEstimates) <- row.names(NCarRegionMean) <- c(idcountyTable)
colnames(NCarRegionEstimates) <- colnames(NCarRegionMean) <- years

##-- Fill-in the table
for(t in 1:n.years){
  for(i in 1:length(idcountyTable)){
    NCarRegionEstimates[idcountyTable[i],t] <- 
      paste0(round(ACdensity[[t+2]]$summary[idcountyTable[i],"mean"], digits = 1), " (",
             round(ACdensity[[t+2]]$summary[idcountyTable[i],"95%CILow"], digits = 0), "-",
             round(ACdensity[[t+2]]$summary[idcountyTable[i],"95%CIHigh"], digits = 0), ")")
    
    NCarRegionMean[idcountyTable[i],t] <- round(ACdensity[[t+2]]$summary[idcountyTable[i],"mean"])
  }#i
}#t

##-- Adjust names
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
row.names(NCarRegionEstimates) <- idcounty1

##-- print .csv
write.csv(NCarRegionEstimates,
          file = file.path(working.dir, "tables",
                           paste0(modelName, "","_NAllYears.csv")))
write.csv(NCarRegionMean,
          file = file.path(working.dir, "tables",
                           paste0(modelName, "","_NAllYears_mean.csv")))

##-- print .tex
row.names(NCarRegionEstimates) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
print(xtable( NCarRegionEstimates,
              type = "latex",
              align = paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
      floating = FALSE,
      sanitize.text.function = function(x){x},
      add.to.row = list(list(seq(1, nrow(NCarRegionEstimates), by = 2)), "\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables",
                       paste0(modelName, "","_NAllYears.tex")))




## ------       2.1.2. LAST YEAR N PER SEX PER COUNTY -----

NCountyEstimatesLastRegions <- matrix("", ncol = 3, nrow = length(idcountyTable))
row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")

for(i in 1:length(idcountyTable)){
  ##-- FEMALES
  NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <- 
    paste0(round(ACdensityF[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(ACdensityF[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(ACdensityF[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  
  ##-- MALES
  NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <- 
    paste0(round(ACdensityM[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(ACdensityM[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(ACdensityM[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  
  ##-- BOTH SEXES
  NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <- 
    paste0(round(ACdensity[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(ACdensity[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(ACdensity[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}#i

##-- ADJUST NAMES 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
row.names(NCarRegionEstimates) <- idcounty1

##-- print .csv
write.csv( NCountyEstimatesLastRegions,
           file = file.path(working.dir, "tables",
                            paste0(modelName, "","_NLastYearPerSex.csv")))

##-- print .tex
row.names(NCountyEstimatesLastRegions) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
print(xtable(NCountyEstimatesLastRegions, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
      sanitize.text.function=function(x){x},
      floating = FALSE,
      add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path( working.dir, "tables",
                        paste0(modelName, "","_NLastYearPerSex.tex")))



##-- UD-Density 
NCountyEstimatesLastRegions_UD <- matrix("", ncol = 3, nrow = length(idcountyTable))
row.names(NCountyEstimatesLastRegions_UD) <- c(idcountyTable)
colnames(NCountyEstimatesLastRegions_UD) <- c("Females","Males","Total")

for(i in 1:length(idcountyTable)){
  ##-- FEMALES
  NCountyEstimatesLastRegions_UD[idcountyTable[i],"Females"] <- 
    paste0(round(spaceUSEDF[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(spaceUSEDF[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(spaceUSEDF[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  
  ##-- MALES
  NCountyEstimatesLastRegions_UD[idcountyTable[i],"Males"] <- 
    paste0(round(spaceUSEDM[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(spaceUSEDM[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(spaceUSEDM[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  
  ##-- BOTH SEXES
  NCountyEstimatesLastRegions_UD[idcountyTable[i],"Total"] <- 
    paste0(round(spaceUSED[[n.years+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
           round(spaceUSED[[n.years+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
           round(spaceUSED[[n.years+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
}#i


##-- ADJUST NAMES 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
row.names(NCarRegionEstimates) <- idcounty1

##-- print .csv
write.csv( NCountyEstimatesLastRegions_UD,
           file = file.path(working.dir, "tables", paste0(modelName, "","_NLastYearPerSex_UD.csv")))

##-- print .tex
row.names(NCountyEstimatesLastRegions_UD) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
print(xtable(NCountyEstimatesLastRegions_UD, type = "latex",
             align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions_UD))), collapse = "")),
      sanitize.text.function=function(x){x},
      floating = FALSE,
      add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions_UD),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path( working.dir, "tables",
                        paste0(modelName, "","_NLastYearPerSex_UD.tex")))



## ------       2.1.4. ALL YEARS N PER SEX PER COUNTY ------

NCountyEstimatesAllSexRegions <- matrix("", ncol = n.years*3, nrow = length(idcountyTable)+1)
row.names(NCountyEstimatesAllSexRegions) <- c("", idcountyTable)
colnames(NCountyEstimatesAllSexRegions) <- rep(years, each = 3)
NCountyEstimatesAllSexRegions[1, ] <- rep(c("Females","Males","Total"), n.years)

##-- Fill-in table
for(t in 1:n.years){
  cols <- which(colnames(NCountyEstimatesAllSexRegions) %in% years[t])
  for( i in 1:length(idcountyTable)){
    ##-- FEMALES
    colss <- which(NCountyEstimatesAllSexRegions[1,cols] %in% "Females")
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- 
      paste(round(ACdensityF[[t+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
            round(ACdensityF[[t+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
            round(ACdensityF[[t+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
    
    ##-- MALES
    colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- 
      paste(round(ACdensityM[[t+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
            round(ACdensityM[[t+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
            round(ACdensityM[[t+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
    
    ##-- TOTAL
    colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
    NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <- 
      paste(round(ACdensity[[t+2]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
            round(ACdensity[[t+2]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
            round(ACdensity[[t+2]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
  }#i
}#t

##-- ADJUST NAMES 
idcounty1 <- idcountyTable
idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
row.names(NCountyEstimatesAllSexRegions) <- c("", idcounty1)

##-- print .csv
write.csv( NCountyEstimatesAllSexRegions,
           file = file.path(working.dir, "tables", paste0(modelName, "","_NAllYearsPerSex.csv")))


##-- print .tex
row.names(NCountyEstimatesAllSexRegions) <- c("", paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")

print(xtable( NCountyEstimatesAllSexRegions, type = "latex",
              align = paste(c("l",rep("c",ncol(NCountyEstimatesAllSexRegions))),collapse = "")),
      sanitize.text.function = function(x){x},
      floating = FALSE,
      add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
      file = file.path(working.dir, "tables", paste0(modelName, "","_NAllYearsPerSex.tex")))




## ------     2.2. NUMBER OF NGS SAMPLES, IDs & DEAD RECOVERIES ------

##-- Load filtered datasets
load(file.path(working.dir, "data", paste0(modelNameM,"_NGSData.RData")))
myData.alive_M <- myData.alive$myData.sp
myData.dead_M <- myData.dead

load(file.path(working.dir, "data", paste0(modelNameF,"_NGSData.RData")))
myData.alive_F <- myData.alive$myData.sp
myData.dead_F <- myData.dead

rm(myData.alive, myData.dead)

##-- SOME TALLIES TO CHECK THINGS
##-- NGS
NGS <- rbind(myData.alive_F, myData.alive_M)
table(NGS$Year)
nrow(NGS)

##-- FOR REPORT SUMMARY
length(NGS$Id)
length(NGS$Id[NGS$Sex=="Hunn"])
length(NGS$Id[NGS$Sex=="Hann"])
length(NGS$Id[NGS$Country=="(S)"])/nrow(NGS)
length(unique(NGS$Id))
length(unique(NGS$Id[NGS$Sex=="Hunn"]))
length(unique(NGS$Id[NGS$Sex=="Hann"]))

##-- DEAD RECOVERY
dead <- rbind(myData.dead_F, myData.dead_M)
table(dead$Year)
nrow(dead)
dead <- dead[dead$Sex %in% c("Hunn","Hann"), ]
table(dead$Year)
length(dead)
length(unique(dead$Id[dead$Sex=="Hunn"]))
length(unique(dead$Id[dead$Sex=="Hann"]))




## ------       2.2.1. NGS SAMPLES & IDs ------

NGS_SEX <- matrix("", ncol = n.years*2, nrow = 3)
row.names(NGS_SEX) <- c( "",
                         "number of NGS samples",
                         "number of NGS individuals")
colnames(NGS_SEX) <- rep(years, each = 2)
NGS_SEX[1, ] <- rep(c("F","M"), n.years)

sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1, n.years*2, by = 2)
for(s in 1:2){
  for(t in 1:n.years){
    temp <- NGS[NGS$Year == years[t] & NGS$Sex == sex[s], ]
    NGS_SEX["number of NGS samples", ye[t] + sex1[s]] <- nrow(temp)
    NGS_SEX["number of NGS individuals", ye[t] + sex1[s]] <- length(unique(temp$Id))
  }#t
}#s

##-- print .csv
write.csv( NGS_SEX, file = file.path(working.dir, "tables", paste0("NGS_SEX.csv")))

##-- print .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGS_SEX))),
                                    '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
colnames(NGS_SEX) <- rep("", ncol(NGS_SEX))
print(xtable(NGS_SEX, type = "latex",
             align = paste(c("l",rep("c",ncol(NGS_SEX))), collapse = "")),
      floating = FALSE, include.colnames = FALSE,
      add.to.row = addtorow,
      file = file.path(working.dir, "tables", paste0("NGS_SEX.tex")))




## ------       2.2.2. DEAD RECOVERIES by CAUSE ------

Dead_SEX <- matrix(0, ncol = n.years*2+1, nrow = 6)
row.names(Dead_SEX) <- c("","other","other","legal culling","legal culling","")
colnames(Dead_SEX) <- c("",unlist(lapply(years, function(x) c(x,x))))
Dead_SEX[1,] <- c("",rep(c("F","M"),n.years))
Dead_SEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
sex <- c("Hunn","Hann")
sex1 <- c(0,1)
ye <- seq(1,n.years*2,by=2)

##-- Define legal mortality causes
MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
MortalityNames[!MortalityNames %in% legalCauses]

##-- Separate mortalities
cause <- c("other","legal culling")
for(t in 1:n.years){
  for(s in 1:2){
    for(d in 1:2){
      if(d==1){
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$Death_cause %in% legalCauses), ]
      } else {
        temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$Death_cause %in% legalCauses, ]
      }
      row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[,1]=="Norway")
      Dead_SEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country_sample %in% "(N)" ]))
      
      row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[,1]=="Sweden" )
      Dead_SEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country_sample %in% "(S)"]))
    }#t
    Dead_SEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(Dead_SEX[2:6,ye[t] + sex1[s]+1]))
  }
}


##-- summary
##-- Other causes
sum(as.numeric(Dead_SEX[2:3,2:ncol(Dead_SEX)]))
sum(as.numeric(Dead_SEX[2:3,which(Dead_SEX[1,]=="F")]))
sum(as.numeric(Dead_SEX[2:3,which(Dead_SEX[1,]=="M")]))
##-- legal
sum(as.numeric(Dead_SEX[4:5,2:ncol(Dead_SEX)]))
sum(as.numeric(Dead_SEX[4:5,which(Dead_SEX[1,]=="F")]))
sum(as.numeric(Dead_SEX[4:5,which(Dead_SEX[1,]=="M")]))

sum(as.numeric(Dead_SEX[c(2,3),2:ncol(Dead_SEX)]))/sum(as.numeric(Dead_SEX[c(2:5),2:ncol(Dead_SEX)]))

##-- %of dead reco (legal) in norway
sum(as.numeric(Dead_SEX[4,2:ncol(Dead_SEX)]))/sum(as.numeric(Dead_SEX[c(4,5),2:ncol(Dead_SEX)]))
sum(as.numeric(Dead_SEX[6,which(Dead_SEX[1,]=="M")]))
sum(as.numeric(Dead_SEX[6,which(Dead_SEX[1,]=="F")]))
sum(as.numeric(Dead_SEX[6,which(Dead_SEX[1,] %in% c("F","M"))]))

##-- print .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
uniqueYEAR <- sort(unique(colnames(Dead_SEX)))
uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",
                             paste0('& \\multicolumn{2}{c}{',
                                    uniqueYEAR,
                                    '}', collapse=''), '\\\\'),
                      rep("\\rowcolor[gray]{.95}",1))
multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{TOTAL}}"), ncol = 1)
Dead_SEX <- data.frame(cbind(multirowadd,Dead_SEX))

print(xtable(Dead_SEX, type = "latex",
             align = rep("c", ncol(Dead_SEX)+1)),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = FALSE,
      include.rownames = FALSE,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "DeadidCountrySEX.tex"))



## ------       2.2.3. PROPORTION OF INDIVIDUALS DETECTED ------

##-- Get the number of individuals detected each year
n.detected_F <- apply(nimDataF$y.alive[ ,1, ], 2, function(x)sum(x>0))
colSums(nimDataF$y.alive[ ,1, ]>0)
n.detected_F

n.detected_M <- apply(nimDataM$y.alive[ ,1, ], 2, function(x)sum(x>0))
colSums(nimDataM$y.alive[ ,1, ]>0)
n.detected_M

propDetected <- matrix("", ncol = n.years, nrow = 3)
row.names(propDetected) <- c("F","M","Total")
colnames(propDetected) <- years
for(t in 1:n.years){
  propDetected["F",t] <- getCleanEstimates(n.detected_F[t+2]/colSums(ACdensityF[[t+2]]$PosteriorRegions))
  propDetected["M",t] <- getCleanEstimates(n.detected_M[t+2]/colSums(ACdensityM[[t+2]]$PosteriorRegions))
  propDetected["Total",t] <- getCleanEstimates((n.detected_F[t+2]+n.detected_M[t+2])/
                                                 (colSums(ACdensityF[[t+2]]$PosteriorRegions)+
                                                    colSums(ACdensityM[[t+2]]$PosteriorRegions)))
}#t

##-- print .csv
write.csv(propDetected,
          file = file.path(working.dir, "tables", "PropDetectedIds.csv"))

##-- print .tex
print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(propDetected), by = 2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "PropDetectedIds.tex"))




## ------       2.2.4. NUMBER OF IDs w/ ACs OUTSIDE NORWAY ------

##-- Prepare raster of countries
countryRaster <- habitatRasterResolution$`5km`[["Countries"]]

##-- Calculate number of individuals alive with their AC in each country each year
N_det_by_country <- matrix(NA,5,n.years)
dimnames(N_det_by_country) <- list("Countries" = c("Norway","Sweden","Finland","Russia","Out"),
                                   "Years" = c(years))
for(t in 1:n.years){
  N_fin_F <- N_fin_M <- N_fin <- rep(NA,n.iter)
  N_nor_F <- N_nor_M <- N_nor <- rep(NA,n.iter)
  N_rus_F <- N_rus_M <- N_rus <- rep(NA,n.iter)
  N_swe_F <- N_swe_M <- N_swe <- rep(NA,n.iter)
  N_out_F <- N_out_M <- N_out <- rep(NA,n.iter)
  for(iter in 1:n.iter){
    
    country <- countryRaster[cellFromXY(norRaster,myResultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+2])]
    isFin <- country %in% 1
    isNor <- country %in% 2
    isRus <- country %in% 3
    isSwe <- country %in% 4
    
    ##-- Detected individuals
    N_fin_F[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isFin & isFemale)
    N_fin_M[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isFin & isMale)
    N_fin[iter] <- N_fin_F[iter] + N_fin_M[iter]
    
    N_nor_F[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor & isFemale)
    N_nor_M[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor & isMale)
    N_nor[iter] <- N_nor_F[iter] + N_nor_M[iter]
    
    N_rus_F[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t] & isRus & isFemale)
    N_rus_M[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t] & isRus & isMale)
    N_rus[iter] <- N_rus_F[iter] + N_rus_M[iter]
    
    N_swe_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t+2] & isSwe & isFemale)
    N_swe_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t+2] & isSwe & isMale)
    N_swe[iter] <- N_swe_F[iter] + N_swe_M[iter]
    
    N_out_F[iter] <- N_fin_F[iter] + N_rus_F[iter] + N_swe_F[iter]
    N_out_M[iter] <- N_fin_M[iter] + N_rus_M[iter] + N_swe_M[iter]
    N_out[iter] <- N_out_F[iter] + N_out_M[iter]
  }#iter
  
  N_det_by_country["Norway",t] <- getCleanEstimates(N_nor)
  N_det_by_country["Sweden",t] <- getCleanEstimates(N_swe)
  N_det_by_country["Finland",t] <- getCleanEstimates(N_fin)
  N_det_by_country["Russia",t] <- getCleanEstimates(N_rus)
  N_det_by_country["Out",t] <- getCleanEstimates(N_out)
  
  print(t)
}#t
##-- Print number of individuals detected per country
print(N_det_by_country)


##-- Calculate number of individuals detected/undetected in Norway
N_NOR <- matrix(NA,4,n.years)
dimnames(N_NOR) <- list("#individuals" = c("Detected","Undetected","Total","%"),
                        "Years" = c(years))
for(t in 1:n.years){
  N_det <- N_undet <- N_tot <- rep(NA,n.iter)
  for(iter in 1:n.iter){
    country <- countryRaster[cellFromXY(norRaster,myResultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+2])]
    isNor <- country %in% 2
    
    ##-- Detected individuals
    N_det[iter] <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor)
    
    ##-- Undetected individuals
    N_undet[iter] <- sum(!isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor)
    
    ##-- Total Norway
    N_tot[iter] <- N_det[iter] + N_undet[iter]
  }#iter
  
  N_NOR["Detected",t] <- getCleanEstimates(N_det)
  N_NOR["Undetected",t] <- getCleanEstimates(N_undet)
  N_NOR["Total",t] <- getCleanEstimates(N_tot)
  
  print(t)
}#t
##-- Print number of individuals with AC in Norway
print(N_NOR)



##-- Calculate % of the Norwegian bear population detected 
prop <- matrix(NA,3,n.years)
dimnames(prop) <- list("% individuals" = c("F","M","Total"),
                       "Years" = c(years))
for(t in 1:n.years){
  prop_F <- prop_M <- prop_tot <- rep(NA,n.iter)
  for(iter in 1:n.iter){
    country <- countryRaster[cellFromXY(norRaster,myResultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+2])]
    isNor <- country %in% 2
    
    ##-- Detected female
    N_F <- sum(isAlive[iter, ,t+2] & isNor & isFemale)
    N_det_F <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor & isFemale)
    prop_F[iter] <- N_det_F / N_F 
    
    ##-- Detected male
    N_M <- sum(isAlive[iter, ,t+2] & isNor & isMale)
    N_det_M <- sum(isDetected[ ,t+2] & isAlive[iter, ,t+2] & isNor & isMale)
    prop_M[iter] <- N_det_M / N_M 
    
    ##-- Detected total
    prop_tot[iter] <- (N_det_F + N_det_M) / (N_F + N_M) 
  }#iter
  
  prop["F",t] <- getCleanEstimates(prop_F)
  prop["M",t] <- getCleanEstimates(prop_M)
  prop["Total",t] <- getCleanEstimates(prop_tot)
  
  print(t)
}#t
##-- Print number of individuals with AC in Norway
print(prop)

##-- print .csv
write.csv(prop,
          file = file.path(working.dir, "tables", "PropDetectedNorway.csv"))

##-- print .tex
print(xtable(prop, type = "latex", align=paste(c("l",rep("c",ncol(prop))),collapse = "")),
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(prop), by = 2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "PropDetectedNorway.tex"))



## ------     2.3. VITAL RATES ------

parameters <- c("rho","phi", "h", "w")
sex <- c("F", "M")
vitalRate <- matrix(NA, nrow = length(parameters)+1, ncol = (n.years)*2-2)
rownames(vitalRate) <- c("", unlist(lapply(as.list(parameters), function(x)rep(x,1))))
colnames(vitalRate) <- c(unlist(lapply(years[1:(length(years)-1)], function(x)rep(paste(x,x+1,sep=" to "),2))))
vitalRate[1, ] <- c(rep(sex,(n.years-1)) )

for(s in 1:2){
  if(s == 1){results <- myResults_F} else {results <- myResults_M}
  
  col <- which(vitalRate[1, ] == sex[s])
  
  ##-- Per capita recruitment
  if(any(grep("rho",names(results$sims.list)))){
    vitalRate["rho",col] <- getCleanEstimates(results$sims.list$rho, moment = "median")
  } else {
    if(s == 1){
      for(t in 1:(n.years-1)){
        n.recruits <- rowSums(isAvail[ ,isFemale,t+2] * isAlive[ ,isFemale,t+3])
        alivetminus1 <- rowSums(isAlive[ ,isFemale,t+2])
        vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
      }#t
    } else {
      for(t in 1:(n.years-1)){
        n.recruits <- rowSums(isAvail[ ,isMale,t+2] * isAlive[ ,isMale,t+3])
        alivetminus1 <- rowSums(isAlive[ ,isMale,t+2])
        vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
      }#t
    }
  }             
  
  
  ##-- Mortality & Survival
  if(any(grep("mhH",names(results$sims.list)))){
    ##-- Calculate mortality from estimated hazard rates (mhH and mhW)
    mhH1 <- exp(results$sims.list$mhH[ ,3:11])
    mhW1 <- exp(results$sims.list$mhW[ ,3:11])
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    vitalRate["phi",col] <- apply(phi, 2, function(x) getCleanEstimates(x, moment = "median"))
    vitalRate["h",col] <- apply(h, 2, function(x) getCleanEstimates(x, moment = "median"))
    vitalRate["w",col] <- apply(w, 2, function(x) getCleanEstimates(x, moment = "median"))
  } else {
    if(s == 1){
      y.dead <- y.deadF[ ,3:12]
      z <- myResultsSXYZ_MF$sims.list$z[ ,isFemale,3:12] 
    } else {
      y.dead <- y.deadM[ ,3:12]
      z <- myResultsSXYZ_MF$sims.list$z[ ,isMale,3:12] 
    }
    
    ##-- Extract survival from posteriors
    vitalRate["phi",col] <- apply(results$sims.list$phi, 2, 
                                  function(x) getCleanEstimates(x, moment = "median"))
    
    ##-- Derive mortality from posterior z and dead recoveries 
    isDead <- apply((z[ , ,1:(n.years-1)] == 2)*(z[ , ,2:n.years] == 3), c(1,3), sum)
    wasAlive <- apply(z[ , ,1:(n.years-1)] == 2, c(1,3), sum)
    mortality <- isDead / wasAlive
    h <- sapply(1:(n.years-1), function(t)sum(y.dead[ ,t+1])/wasAlive[ ,t])
    w <- mortality - h
    vitalRate["h",col] <- apply(h, 2, function(x) getCleanEstimates(x, moment = "median"))
    vitalRate["w",col] <- apply(w, 2, function(x) getCleanEstimates(x, moment = "median"))
  }#else
}#s

##-- Print .csv
write.csv( vitalRate,
           file = file.path(working.dir, "tables", "VitalRates.csv"))

##-- Print .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(vitalRate))),
                                    '}', collapse = ''), '\\\\'), rep("\\rowcolor[gray]{.95}",1))
colnames(vitalRate) <- rep("", ncol(vitalRate))
rownames(vitalRate)[2:5] <- c("$\\rho$","$\\phi$","h","w")

print(xtable(vitalRate, type = "latex",
             align = paste(rep("c", ncol(vitalRate)+1),collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "VitalRates.tex"))




## ------     2.4. DERIVED PARAMETERS FROM ABUNDANCE ------

## ------       2.4.1. DERIVE SEX-RATIO ------

##-- REGION-SPECIFIC PROPORTION OF FEMALES
PropFemale_regions <- list()
for(t in 1:n.years){
  PropFemale_regions[[t]] <- ACdensityF[[t+2]]$PosteriorRegions/
    (ACdensityM[[t+2]]$PosteriorRegions +
       ACdensityF[[t+2]]$PosteriorRegions)
  print(rowMeans(PropFemale_regions[[t]],na.rm = T))
}#t


##-- OVERALL PROPORTION OF FEMALES
PropFemale <- list()
for(t in 1:n.years){
  PropFemale[[t]] <- colSums(ACdensityF[[t+2]]$PosteriorAllRegions)/
    (colSums(ACdensityM[[t+2]]$PosteriorAllRegions) +
       colSums(ACdensityF[[t+2]]$PosteriorAllRegions))
}#t


##-- Format table
propFemale_tab <- matrix(0, ncol = n.years, nrow = 8)
row.names(propFemale_tab) <- idcountyTable
colnames(propFemale_tab) <- years
for(t in 1:n.years){
  for(c in 1:7){
    propFemale_tab[c,t] <- getCleanEstimates(na.omit(PropFemale_regions[[t]][c, ])) 
  }#c
  propFemale_tab[8,t] <- getCleanEstimates(PropFemale[[t]] ) 
}#t


##-- print .tex
row.names(propFemale_tab) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
print(xtable( propFemale_tab,
              type = "latex",
              align = paste(c("l",rep("c",ncol(propFemale_tab))),collapse = "")),
      floating = FALSE,
      sanitize.text.function = function(x){x},
      add.to.row = list(list(seq(1, nrow(propFemale_tab), by = 2)), "\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "propFemale.tex"))




## ------       2.4.2. DERIVE DENSITY ------

habbRCarRegionsTRY <- rrRegions
habbRCarRegionsTRY[!is.na(habbRCarRegionsTRY[])] <- 1
habbRCarRegionsTRY[] <- as.numeric(habbRCarRegionsTRY[])
areaSqKm <- sum(na.omit(habbRCarRegionsTRY[ ] > 0)) * 25

##-- Multiplied by 100 to get per 100km2
ACdensity[[12]]$summary["Total","mean"]/areaSqKm*100
ACdensity[[12]]$summary["Total","95%CILow"]/areaSqKm*100
ACdensity[[12]]$summary["Total","95%CIHigh"]/areaSqKm*100




## ------       2.4.3. PROPORTION OF INDIVIDUALS DETECTED ------

##-- Get the number of individuals detected each year
n.detected_F <- apply(nimDataF$y.alive[ ,1, ], 2, function(x)sum(x>0))
n.detected_M <- apply(nimDataM$y.alive[ ,1, ], 2, function(x)sum(x>0))

propDetected <- matrix("", ncol = n.years, nrow = 3)
row.names(propDetected) <- c("F","M","Total")
colnames(propDetected) <- years
for(t in 1:n.years){
  propDetected["F",t] <- getCleanEstimates(n.detected_F[t+2]/colSums(ACdensityF[[t+2]]$PosteriorRegions))
  propDetected["M",t] <- getCleanEstimates(n.detected_M[t+2]/colSums(ACdensityM[[t+2]]$PosteriorRegions))
  propDetected["Total",t] <- getCleanEstimates((n.detected_F[t+2]+n.detected_M[t+2])/
                                                 (colSums(ACdensityF[[t+2]]$PosteriorRegions)+
                                                    colSums(ACdensityM[[t+2]]$PosteriorRegions)))
}#t

##-- print .csv
write.csv(propDetected,
          file = file.path(working.dir, "tables", "PropDetectedIds.csv"))

##-- print .tex
print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
      floating = FALSE, sanitize.text.function=function(x){x},
      add.to.row = list(list(seq(1,nrow(propDetected), by = 2)),"\\rowcolor[gray]{.96} "),
      file = file.path(working.dir, "tables", "PropDetectedIds.tex"))




## ------       2.4.4. GROWTH RATE ------

growthRate <- list()
for(t in 1:(n.years-1)){
  growthRate[[t]] <- colSums(ACdensity[[t+3]]$PosteriorAllRegions)/
    colSums(ACdensity[[t+2]]$PosteriorAllRegions) 
}#t

##-- Put in a table format
growthRate_tab <- matrix(0, ncol = (n.years-1), nrow = 1)
colnames(growthRate_tab) <- paste(years[-n.years], years[-1], sep = " to ")
for(t in 1:(n.years-1)){
  growthRate_tab[1,t] <- getCleanEstimates(growthRate[[t]]) 
}#t

##-- Print .tex
addtorow <- list()
addtorow$pos <- list(c(0),0)
addtorow$command <- c(paste0(paste('& {', sort(unique(colnames(growthRate_tab))),
                                   '}', collapse = ''), '\\\\'), rep("\\rowcolor[gray]{.95}",1))
colnames(growthRate_tab) <- rep("", ncol(growthRate_tab))
rownames(growthRate_tab) <- c("$\\lambda$")

print(xtable(growthRate_tab, type = "latex",
             align = paste(c("l", rep("c",ncol(growthRate_tab))), collapse = "")),
      floating = FALSE,
      add.to.row = addtorow,
      include.colnames = F,
      sanitize.text.function = function(x){x},
      file = file.path(working.dir, "tables", "GrowthRates.tex"))




## ------     2.5. TABLE OTHERS ------

parameters <- c("tau",
                "betaDead","betaDens",
                "betaDead","betaDens",
                "sigma",
                "betaDet","betaDet")
sex <- c("F","M")
TableOthers <- matrix(NA, nrow = length(parameters), ncol = 3)
rownames(TableOthers) <- parameters
colnames(TableOthers) <- c("", sex)
TableOthers[ ,1] <- c("$\\tau$",
                      "$\\beta_{dead}_{1}$","$\\beta_{skandobs}_{1}$",
                      "$\\beta_{dead}_{2}$","$\\beta_{skandobs}_{2}$",
                      "$\\sigma$",
                      "$\\beta_{roads}$","$\\beta_{obs}$")

for(s in 1:2){
  if(s == 1){results <- myResults_F} else {results <- myResults_M}
  TableOthers["tau",sex[s]] <- getCleanEstimates(results$sims.list$tau/1000, moment = "median")
  TableOthers[which(parameters == "betaDead"),sex[s]] <- apply(results$sims.list$betaDens[,1,], 2,
                                                               function(x) getCleanEstimates(x,moment = "median"))
  TableOthers[which(parameters == "betaDens"),sex[s]] <- apply(results$sims.list$betaDens[,2,],
                                                               2, function(x) getCleanEstimates(x,moment = "median"))
  TableOthers["sigma",sex[s]] <- getCleanEstimates(results$sims.list$sigma/1000, moment = "median")
  TableOthers[which(parameters == "betaDet"),sex[s]] <- apply(results$sims.list$betaDet, 2,
                                                              function(x) getCleanEstimates(x,moment = "median"))
}#s
##-- Deal with negative values 
TableOthers <- gsub("--", "-(-)", TableOthers)

##-- Change row names
row.names(TableOthers) <- c("Spatial process","","","","",
                            "Detection Process","","")

##-- Print .tex
multirow <- c("\\multirow{5}{*}{\\textbf{Spatial process}}","\\multirow{3}{*}{\\textbf{Detection process}}")
multirowadd <- matrix(c(multirow[1],"","","","",multirow[2],"",""), ncol = 1)
TableOthers <- data.frame(cbind(multirowadd, TableOthers))

addtorow <- list()
addtorow$pos <- list(0,5)
addtorow$command <- c(paste0("& {\\textbf{Parameters}} ",
                             paste0('& {\\textbf{',  sex, '}}', collapse = ''), '\\\\'),
                      "\\hline")

print(xtable(TableOthers, type = "latex",
             align = paste(c("ll", rep("c",ncol(TableOthers)-1)), collapse = "")),
      sanitize.text.function = function(x){x},
      floating = FALSE,
      include.rownames = FALSE,
      include.colnames = FALSE,
      add.to.row = addtorow,
      file = file.path(working.dir, "tables", "TableParametersOthers.tex"))




## -----------------------------------------------------------------------------



## ------       1.1.1. NGS, Dead recoveries & Carnivore obs ------
### NEEDS TO GO INTO makeRovquantData()
{
  myFullData.sp <- readMostRecent( 
    path = file.path(dir.dropbox,"DATA/RovbaseData_clean/bear"),
    extension = ".RData")
  
  myFilteredData.sp <- myFullData.sp
  detArea <- st_as_sf(aggregate(rasterToPolygons(myHabitat$habitat.rWthBuffer,function(x) x==1)))
  myFilteredData.sp$alive <- myFullData.sp$alive[!is.na(as.numeric(st_intersects(myFullData.sp$alive,detArea))), ]
  myFilteredData.sp$dead.recovery <- myFullData.sp$dead.recovery[!is.na(as.numeric(st_intersects(myFullData.sp$dead.recovery, detArea))), ]
  
  ##-- Plot NGS & Dead recovery maps
  pdf(file = file.path(working.dir, "figures", paste0("NGS_DR_maps.pdf")),
      width = 18, height = 12)
  
  ##-- layout
  mx <- rbind(c(1,rep(1:5, each = 2)),
              c(rep(1:5, each = 2), 5))
  mx <- rbind(mx, mx + 5)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  par(mar = c(0,0,0,0))
  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
    points(myFilteredData.sp$alive[myFilteredData.sp$alive$Year == years[t], ],
           pch = 3, col = "orange", lworking.dir = 0.7)
    points(myFilteredData.sp$dead.recovery[myFilteredData.sp$dead.recovery$Year == years[t], ],
           pch = 3, col = "slateblue", lworking.dir = 0.7)
    mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
    
    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lworking.dir = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
      
      ##-- LEGEND
      par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
      plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
      points(c(3,3), c(2.6,1.8), pch = 3, lworking.dir = 1.5, cex = 3, col = c("orange","slateblue"))
      text(c(3.5,3.5), c(2.55,1.75),  c("NGS samples", "Dead recoveries"), cex = 2, pos = 4)
    }#if
  }#t
  dev.off()
  
  
  
  ##-- Plot Carnivore observations maps
  pdf(file = file.path(working.dir, "figures", paste0("CarnivoreObs_maps.pdf")),
      width = 18, height = 12)
  
  ##-- layout
  mx <- rbind(c(1,rep(1:5, each = 2)),
              c(rep(1:5, each = 2), 5))
  mx <- rbind(mx, mx + 5)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  par(mar = c(0,0,0,0))
  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
    image(mask(ds.brickCont[[t]],COUNTRIESsimpFig[1,]), add = TRUE, col = c("white","forestgreen"), legend = FALSE)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = "gray80", col = NA, add = TRUE)
    
    mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
    
    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lworking.dir = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
      
      ##-- LEGEND
      par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
      plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
    }#if
  }#t
  dev.off()
  
}




