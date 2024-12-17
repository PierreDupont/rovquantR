#' @title RovQuant OPSCR wolverine data preparation.
#'
#' @description
#' \code{makeRovquantData_wolverine} formats the available bear data for the OPSCR analysis using nimble and nimbleSCR.
#' The data preparation process is composed of three main steps:
#'  - defining and formatting habitat characteristics
#'  - defining and formatting detectors characteristics
#'  - defining and formatting individual detection histories
#'
#' @name makeRovquantData_wolverine
#'
#' @param data.dir A \code{path}.
#' @param working.dir A \code{path}.
#' @param years A \code{list}.  
#' @param sex A \code{character}.
#' @param aug.factor A \code{Numeric}.
#' @param sampling.months A \code{list}.
#' @param habitat.res A \code{Numeric}.  
#' @param buffer.size A \code{Numeric}.
#' @param max.move.dist A \code{Numeric}.
#' @param detector.res A \code{Numeric}.
#' @param subdetector.res A \code{Numeric}.
#' @param max.det.dist A \code{Numeric}.  
#' @param resize.factor A \code{Numeric}.
#' @param plot.check A \code{Logical}.
#' @param print.report A \code{Logical}.
#' 
#' 
#' @return 
#' A \code{html} report summarizing the data preparation process
#' Additional \code{.png} images that can be reused somewhere else.
#'
#' @author Pierre Dupont
#' 
#' @import sf 
#' @import raster
#' @import dplyr
#' @importFrom adehabitatHR estUDm2spixdf kernelUD
#' @importFrom fasterize fasterize
#' @importFrom nimbleSCR getSparseY scaleCoordsToHabitatGrid
#' @importFrom sp SpatialPoints CRS
#' @importFrom spatstat.geom as.owin ppp
#' @importFrom spatstat.explore density.ppp
#' @importFrom stars st_as_stars
#' @importFrom stats runif
#' @importFrom stringi stri_trans_general 
#' @importFrom utils data
#' 
NULL
#' @rdname makeRovquantData_wolverine
#' @export
makeRovquantData_wolverine <- function(
    ##-- paths
  data.dir = getwd(),
  working.dir = getwd(),
  
  ##-- data
  years = NULL,
  sex = c("Hann","Hunn"),
  aug.factor = 0.8,
  samplingMonths = list(12,1:6),
  
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 60000,
  max.move.dist = 250000,
  
  ##-- detectors
  detector.res = 10000,
  subdetector.res = 2000,
  max.det.dist = 84000,
  resize.factor = 2 
){
  ## ---------------------------------------------------------------------------
  
  ## ------ 0. BASIC SET-UP ------
  
  ##-- Set default values for the wolverine model
  if(is.null(aug.factor)){aug.factor <- 0.8}
  if(is.null(sampling.months)){sampling.months <- list(12,1:6)}
  if(is.null(habitat.res)){habitat.res <- 20000} 
  if(is.null(buffer.size)){buffer.size <- 60000}
  if(is.null(max.move.dist)){max.move.dist <- 250000}
  if(is.null(detector.res)){detector.res <- 10000}
  if(is.null(subdetector.res)){subdetector.res <- 2000}
  if(is.null(max.det.dist)){max.det.dist <- 84000}
  if(is.null(resize.factor)){resize.factor <- 2}
  
  
  ##-- Set up list of Habitat characteristics
  habitat <- list( resolution = habitat.res,
                   buffer = buffer.size,
                   maxDist = max.move.dist)
  
  ##-- Set up list of Detectors characteristics
  detectors <- list( resolution = detector.res,
                     resolution.sub = subdetector.res,
                     maxDist = max.det.dist)
  
  ##-- Set up list of Data characteristics
  data <- list( sex = sex,
                aug.factor = aug.factor,
                sampling.months = sampling.months)
  
  
  
  ## ---------------------------------------------------------------------------
  
  ## ------ I. LOAD AND SELECT DATA ------
  
  ## ------   1. NGS DATA -----
  ##-- Extract date from the last cleaned data file
  DATE <- getMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_wolverine")
  
  ##-- Load the most recent Bear data from RovBase
  myFullData.sp <- readMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_wolverine",
    extension = ".RData")
  
  # ##-- Define legal mortality causes
  # MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
  # legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
  # legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
  # 
  # myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
  #   dplyr::mutate( legal = ifelse(Death_cause %in% legalCauses, "yes", "no"))
  # 
  # legal.death <- dplyr::filter(myFullData.sp$dead.recovery, legal %in% "yes")
  # Other.death <- dplyr::filter(myFullData.sp$dead.recovery, legal %in% "no")
  
  
  if(is.null(years)){
    years <- sort(unique(c(myFullData.sp$alive$Year,
                           myFullData.sp$dead.recovery$Year)))
  }
  data$years <- years
  n.years <- length(years)
  
  
  
  ## ------   2. HABITAT DATA -----
  
  ##-- Load pre-defined habitat rasters and shapefiles
  data(COUNTRIES, envir = environment()) 
  data(COUNTIES, envir = environment()) 
  data(habitatRasters, envir = environment()) 
  data(GLOBALMAP, envir = environment()) 
  
  ##-- Disaggregate habitat raster to the desired resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)
  
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
  
  
  
  ## ------     2.1. DEN COUNTS ------
  
  ##-- Load the last DEN_COUNTS data file
  DEN <- readMostRecent( 
    path = data.dir,
    extension = ".csv",
    pattern = "DEN_COUNTS",
    fileEncoding = "latin1")
  
  
  
  ## ------   3. DETECTORS DATA ----- 
  
  ## ------     3.1. GPS SEARCH TRACKS ------
  
  TRACKS_SINGLE <- sf::st_read(file.path(data.dir, "GIS/SearchTracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20231020_date.shp"))
  TRACKS_MULTI <- sf::st_read(file.path(data.dir, "GIS/SearchTracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20231020_date.shp"))
  
  ##-- Combine all tracks
  ALL_TRACKS <- rbind(TRACKS_SINGLE, TRACKS_MULTI)
  ##-- Remove helicopter tracks
  ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Helikopter == "0", ]
  ##-- Keep only wolverine tracks
  ALL_TRACKS <- ALL_TRACKS[ALL_TRACKS$Jerv == "1", ]
  
  ##-- Split by year
  TRACKS_YEAR <- list()
  for(t in 1:n.years){
    ##-- SUBSET GPS TRACKS TO THE SAMPLING PERIOD
    TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][1] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[1]], ]
    TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in%YEARS[[t]][2] & ALL_TRACKS$Mth%in%myVars$DATA$samplingMonths[[2]], ]
    TRACKS <- rbind(TRACKS_1, TRACKS_2)
    ##-- SIMPLIFY TRACKS SHAPES
    TRACKS <- st_intersection(TRACKS, st_as_sf(myStudyArea))
    ##-- NAME TRACKS
    TRACKS$ID <- 1:nrow(TRACKS)
    TRACKS_YEAR[[t]] <- TRACKS
  }#t
  
  
  
  ## ------     3.2. DAYS OF SNOW -----
  
  ##-- Load raster stack of yearly snow cover
  SNOW <- stack(file.path(data.dir,"GIS/Snow/AverageSnowCoverModisSeason2008_2023_09.tif"))
  
  ##-- Rename layers
  names(SNOW) <- paste(2008:2022, (2008:2022)+1, sep = "_")
  
  ##-- Select years corresponding to the monitoring period
  SNOW <- SNOW[[paste("X", years, "_", years+1, sep = "")]]
  
  ##-- Crop raster stack
  SNOW <- raster::crop(SNOW, c(0,40,55,75))
  
  
  
  
  ## ------     3.3. DISTANCE TO ROADS -----
  
  ##-- Load map of distance to roads (1km resolution)
  DistAllRoads <- raster::raster(file.path(data.dir,"GIS/Roads/MinDistAllRoads1km.tif"))
  
  ##-- Fasterize to remove values that fall in the sea
  r <- fasterize::fasterize(sf::st_as_sf(GLOBALMAP), DistAllRoads)
  r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  DistAllRoads <- r
  
  
  
  ## ------     3.2. SKANDOBS ------
  
  ##-- Load the last SkandObs data file
  skandObs <- readMostRecent( 
    path = file.path(data.dir, "Skandobs"),
    extension = ".xlsx",
    pattern = "Skandobs")
  
  ##-- Replace scandinavian characters
  colnames(skandObs) <- translateForeignCharacters(data = colnames(skandObs))
  
  skandObs <- skandObs %>%
    ##-- Extract important info (e.g. month, year)
    dplyr::mutate( date = as.POSIXct(strptime(date, "%Y-%m-%d")),
                   year = as.numeric(format(date,"%Y")),
                   month = as.numeric(format(date,"%m")),
                   species = stringi::stri_trans_general(species, "Latin-ASCII")) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf(., coords = c("longitude","latitude")) %>%
    sf::st_set_crs(. , value = "EPSG:4326") %>%
    sf::st_transform(. ,sf::st_crs(COUNTIES))
  
  
  
  ## ------     3.3. ROVBASE OBS ------
  
  ##-- GET ALL SAMPLES COLLECTED (all species)
  rovbaseObs <- readMostRecent( path = data.dir,
                                extension = ".csv",
                                pattern = "all_samples") %>%
    ##-- Deal with Scandinavian characters
    dplyr::mutate(Species = stringi::stri_trans_general(Species, "Latin-ASCII")) %>%
    ##-- Filter out samples without coordinates
    dplyr::filter( !is.na(East_UTM33),
                   Species %in% c("Bjorn","Fjellrev","Gaupe","Hund","Jerv","Rodrev","Ulv"),
                   Sample_type %in% c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)",
                                      "Sekret (Jerv)","Saliv/Spytt")) %>%
    ##-- Extract important info (e.g. month, year, country of collection)
    dplyr::mutate( Sample_type = translateForeignCharacters(data = Sample_type),
                   Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
                   year = as.numeric(format(Date,"%Y")),
                   month = as.numeric(format(Date,"%m")),
                   country = substrRight(County,3)) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(. , sf::st_crs(COUNTIES))
  
  
  
  ## ---------------------------------------------------------------------------
  
  ## ------ II. CREATE OPSCR DATA ------
  
  ## ------   1. GENERATE HABITAT ------
  
  message("Preparing habitat characteristics... ")
  
  ## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------
  
  ##-- Determine study area based on NGS detections
  ##-- Buffer NGS detections and cut to Swedish and Norwegian borders
  studyArea <- myFullData.sp$alive %>%  
    sf::st_buffer(., dist = habitat$buffer * 1.4) %>%
    sf::st_union() %>%
    sf::st_intersection(., COUNTRIES) %>%
    sf::st_as_sf() 
  
  ##-- Make habitat from predefined scandinavian raster of suitable habitat
  habitat <- MakeHabitatFromRaster(
    poly = studyArea,
    habitat.r = habRaster,
    buffer = habitat$buffer,
    plot.check = FALSE) %>%
    append(habitat,.)
  
  ##-- ???
  habitat$buffered.habitat.poly <- sf::st_simplify( habitat$buffered.habitat.poly,
                                                    dTolerance = 0)
  sf::st_crs(habitat$buffered.habitat.poly) <- sf::st_crs(habitat$habitat.sp)
  
  ##-- Retrieve number of habitat windows 
  isHab <- habitat$habitat.r[] == 1
  n.habwindows <- habitat$n.habwindows <- sum(isHab)
  habitat$habitat.df <- cbind.data.frame(
    "id" = 1:habitat$n.habwindows,
    "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
    "y" = raster::coordinates(habitat$habitat.r)[isHab,2])
  
  ##-- Make a spatial grid from polygon
  habitat$grid <- sf::st_as_sf(raster::rasterToPolygons(habitat$habitat.r,
                                                        fun = function(x){x>0}))
  habitat$grid$id <- 1:nrow(habitat$grid)
  sf::st_crs(habitat$grid) <- sf::st_crs(habitat$buffered.habitat.poly)
  
  
  ##-- Rescale habitat coordinates
  scaledHabGridCenters <- scaleCoordsToHabitatGrid(
    coordsData = habitat$habitat.xy,
    coordsHabitatGridCenter = habitat$habitat.xy,
    scaleToGrid = F)$coordsHabitatGridCenterScaled
  scaledHabGridCenters <- scaledHabGridCenters[habitat$habitat.r[] == 1, ]
  
  ##-- Create habitat matrix 
  habIDCells.mx <- habitat$IDCells.mx 
  habIDCells.mx[] <- 0
  for(i in 1:nrow(scaledHabGridCenters)){
    habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                  trunc(scaledHabGridCenters[i,1])+1] <- i
  }#i
  
  
  
  ## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES -----
  
  ## ------       1.2.1. DEAD RECOVERIES (ALL YEARS) -----
  
  kern <- sf::st_coordinates(myFullData.sp$dead.recovery) %>%
    sp::SpatialPoints(.,  
                      proj4string = sp::CRS(raster::proj4string(habitat$habitat.r)),
                      bbox = NULL) %>%
    adehabitatHR::kernelUD( .,
                            h = 60000,
                            grid = as(habitat$habitat.r, 'SpatialPixels')) %>%
    raster::raster(.)
  
  ##-- Truncate to max value observed in Norway
  kern.truncated <- kern
  maxKern <- max(kern[habitat$habitat.rWthBuffer[] %in% 1])
  kern.truncated[kern.truncated[] > maxKern] <- maxKern
  
  ##-- Scale covariates
  habDens1 <- base::scale(kern[isHab])
  habDens2 <- base::scale(kern.truncated[isHab])
  
  ##-- Put into "nimble2SCR" format
  habitat$habitat.df <- cbind.data.frame( habitat$habitat.df,
                                          "dead.reco.trunc" = habDens2)
  
  
  
  ## ------       1.2.2. NUMBER OF OBSERVATIONS ------
  
  ##-- Rasterize SkandObs bear observations at the habitat level
  rl <- raster::rasterize( x = skandObs[skandObs$species %in% "Bjorn",1],
                           y = habitat$habitat.r,
                           fun = "count")[[1]]
  rl[is.na(rl[ ])] <- 0
  r.skandObsContinuous <- rl
  
  ##-- Binary version  
  rl1 <- rl
  rl1[rl[ ]>0] <- 1
  r.skandObsBinary <- rl1
  
  ##-- Smooth the covariates based on a 5*5 moving window
  mw <- matrix(1,5,5)
  r.skandObsContinuous.smooth <- raster::focal(
    x = r.skandObsContinuous,
    w = mw,
    fun = function(x)mean(x,na.rm = T),
    pad = T)
  
  r.skandObsBinary.smooth <- raster::focal(
    x = r.skandObsBinary,
    w = mw,
    fun = function(x)mean(x,na.rm = T),
    pad = T)
  
  ##-- Set values outside the habitat to NAs
  r.skandObsContinuous[!isHab] <- NA
  r.skandObsBinary[!isHab] <- NA
  r.skandObsContinuous.smooth[!isHab] <- NA
  r.skandObsBinary.smooth[!isHab] <- NA
  
  ##-- Put into "nimble2SCR" format
  habitat$habitat.df <- cbind.data.frame(
    habitat$habitat.df,
    "skandObs.smooth" = base::scale(r.skandObsBinary.smooth[isHab]))
  
  ##-- Merge with the habitat grid
  habitat$grid <- left_join(
    x = habitat$grid,
    y = habitat$habitat.df,
    by = "id")
  
  
  
  
  
  
  ## ------     2.3. SUBSET DETECTIONS BASED ON HABITAT EXTENT ------ 
  
  ##-- Remove samples outside the STUDY AREA 
  whichOut <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$alive, myStudyArea))))
  if(length(whichOut)>0){
    myFilteredData.sp$alive <- myFilteredData.sp$alive[whichOut, ]
  }
  myFilteredData.sp$alive$Id <- droplevels(myFilteredData.sp$alive$Id)
  
  ##-- REMOVE DEAD RECOVERIES OUTSIDE THE HABITAT 
  whichOutBuff <- which(!as.numeric(unlist(st_intersects(myFilteredData.sp$dead.recovery, habitat$buffered.habitat.poly))))
  if(length(whichOutBuff)>0){
    myFilteredData.sp$dead.recovery <- myFilteredData.sp$dead.recovery[whichOutBuff, ]
  }
  
  ##-- Check correlation number of detections ~ between monitoring season
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
  
  
  
  ## ------     2.4. GENERATE HABITAT-LEVEL COVARIATES ------
  
  ## ------       2.4.1. DEN COUNTS ------
  
  DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
  st_crs(DEN.sp) <- st_crs(myFilteredData.sp$alive)
  DEN.sp$id  <- rep(1,nrow(DEN.sp))
  DEN.sp <- DEN.sp[ ,"id"]
  
  DEN.r <- raster(
    estUDm2spixdf(
      kernelUD(as(DEN.sp,"Spatial"), h = 30000,
               grid = as(habitat$habitat.r, 'SpatialPixels'))))
  
  
  ##-- EXTRACT COVARIATES
  denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
  denCounts <- round(scale(denCounts), digits = 2)
  
  
  
  ## ------     2.5. RESCALE HABITAT COORDINATES ------ 
  
  scaledLowerCoords <- scaledHabGridCenters - 0.5
  scaledUpperCoords <- scaledHabGridCenters + 0.5
  
  
  
  
  
  ## ------   2. GENERATE DETECTORS -----
  
  message("Preparing detectors characteristics... ")
  
  ## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----
  
  ##-- Generate raster of sub-detectors based on the study area
  subdetectors.r <- raster::disaggregate(
    x = habitat$habitat.rWthBuffer,
    fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)
  
  ##-- Generate NGS detectors based on the raster of sub-detectors
  detectors <- MakeSearchGrid( 
    data = subdetectors.r,
    resolution = detectors$detResolution,
    div = (detectors$resolution/detectors$resolution.sub)^2,
    plot = FALSE) %>%
    append(detectors,.)
  
  ##-- Extract numbers of detectors
  n.detectors <- detectors$n.detectors <- dim(detectors$main.detector.sp)[1]
  
  ##-- Format detector locations & number of trials per detector
  n.trials <- as.vector(table(detectors$detector.sp$main.cell.id))
  detectors$detectors.df <- cbind.data.frame(
    "id" = 1:detectors$n.detectors,
    "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1],
    "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2],
    "size" = n.trials)
  
  ##-- Make a spatial grid from polygon
  detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
    x = raster::aggregate( x = subdetectors.r,
                           fact = detectors$resolution/detectors$resolution.sub),
    fun = function(x){x>0}))
  detectors$grid$id <- 1:nrow(detectors$grid)
  
  
  
  ## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----
  
  ## ------       2.2.1. EXTRACT COUNTIES -----
  
  ##-- Assign counties to detectors
  dist <- sf::st_distance(detectors$main.detector.sp, COUNTIES)
  detCounties1 <- apply(dist, 1, function(x) COUNTIES$NAME_1[which.min(x)])
  
  ##-- Re-order to account for some counties being never sampled
  detCounties <- as.numeric(as.factor(detCounties1))
  
  ##-- Create a vector of original county names
  detCounties.original <- 0
  for(i in 1: max(detCounties)){
    detCounties.original[i] <- detCounties1[which(detCounties==i)][1]
  }#i
  
  ##-- Put into "nimble2SCR" shape
  detectors$detectors.df$counties <- detCounties
  detectors$detectors.df$counties1 <- detCounties1
  
  
  
  ## ------       2.2.2. EXTRACT DISTANCES TO ROADS -----
  
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
  detectors$detectors.df$roads <- c(detRoads)
  
  # ##-- CHECK IF CONTAINS NAs
  # if(any(is.na(detRoads)))print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")
  
  
  
  ## ------       2.2.3. EXTRACT PRESENCE OF OTHER SAMPLES ------
  
  habitat.rWthBufferPol <- stars::st_as_stars(habitat$habitat.rWthBuffer) %>%
    sf::st_as_sf(., 
                 as_points = FALSE,
                 merge = TRUE) %>%
    dplyr::filter(Habitat %in% 1)
  
  r.detector <- raster::aggregate( 
    subdetectors.r,
    fact = (detectors$resolution/detectors$resolution.sub))
  
  ##-- Subset SkandObs 
  skandObs <- skandObs %>%
    dplyr::filter( 
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months),
      ##-- ... based on space 
      !is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))
  
  # ##-- Rasterize at the detector level
  # r.list <- lapply(data$years, function(y){
  #   if(y %in% skandObs$year){
  #     rl <- raster::rasterize(skandObs[skandObs$year %in% y, 1],
  #                         r.detector,
  #                         fun = "count")[[1]]
  #   } else {
  #     rl <- r.detector
  #     rl[ ] <- 0
  #   }
  #   rl[is.na(rl[ ])] <- 0
  #   rl[!r.detector[ ]%in% 1] <- NA
  #   rl1 <- rl
  #   rl1[rl[ ]>0] <- 1
  #   list(rl1, rl)
  # })
  # r.skandObsSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
  # r.skandObsSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))
  
  
  ##-- Subset RovbaseObs 
  rovbaseObs <- rovbaseObs %>%
    filter( 
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months),
      ##-- ... based on space 
      !is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))
  
  # ##-- Rasterize at the detector level
  # r.list <- lapply(years, function(y){
  #   rl <- raster::rasterize(rovbaseObs[rovbaseObs$year %in% y, 1], 
  #                           r.detector,
  #                           fun = "count")[[1]]
  #   rl[is.na(rl[])] <- 0
  #   rl[!r.detector[]%in% 1] <- NA
  #   rl1 <- rl
  #   rl1[rl[]>0] <- 1
  #   list(rl1, rl)
  # })
  # r.OtherSamplesBinary <- brick(lapply(r.list,function(x) x[[1]]))
  # r.OtherSamplesContinuous <- brick(lapply(r.list,function(x) x[[2]]))
  # 
  # ##-- Combine RovbaseObs and SkandObs  
  # r.SkandObsOtherSamplesBinary <- r.OtherSamplesBinary + r.skandObsSamplesBinary
  # for(t in 1:n.years){
  #   r.SkandObsOtherSamplesBinary[[t]][r.SkandObsOtherSamplesBinary[[t]][]>1 ] <- 1
  # }
  
  
  ##-- Smooth binary map
  ##-- We tried adjust = 0.05, 0.037,0.02 and decided to go for 0.02 
  habOwin <- spatstat.geom::as.owin(as.vector(raster::extent(r.detector)))
  ds.list <- lapply( data$years, function(y){
    ##-- ROVBASE DATA 
    pts <- sf::st_coordinates(rovbaseObs)[rovbaseObs$year %in% y, ]
    ##-- SKANDOBS
    pts <- rbind(pts, sf::st_coordinates(skandObs)[skandObs$year %in% y, ])
    ##-- SMOOTH AND RASTERIZE
    ds <- spatstat.geom::ppp(pts[ ,1], pts[ ,2], window = habOwin) %>%
      spatstat.explore::density.ppp(., adjust = 0.02) %>%    
      raster::raster(.)
    
    ds <- ds1 <- raster::resample(ds, r.detector) #-- mask(ds,rasterToPolygons(habitat$habitat.rWthBuffer,function(x) x==1))
    threshold <- 0.1 / prod(raster::res(ds))      #-- number per 1 unit of the projected raster (meters)
    ds1[] <- ifelse(ds[] < threshold,0,1)
    ds1 <- raster::mask(ds1, raster::rasterToPolygons(habitat$habitat.rWthBuffer, function(x) x==1))
    ds <- raster::mask(ds, raster::rasterToPolygons(habitat$habitat.rWthBuffer, function(x) x==1))
    
    return(list(ds,ds1))
  })
  ds.brick <- raster::brick(lapply(ds.list, function(x) x[[1]]))
  ds.brickCont <- raster::brick(lapply(ds.list, function(x) x[[2]]))
  names(ds.brick) <- years
  
  ##-- Assign covariates to detectors     
  detOtherSamples <- matrix(0, nrow = detectors$n.detectors, ncol = n.years)
  detOtherSamples[ ,1:n.years] <- raster::extract(ds.brickCont,
                                                  detectors$main.detector.sp)
  detectors$detectors.df$detOtherSamples <- detOtherSamples
  
  
  
  ## ------       2.2.4. FORMAT detCovs ------
  
  detCovs <- array(NA, c(n.detectors, 2, n.years))
  for(t in 1:n.years){
    detCovs[ ,1,t] <- detectors$detectors.df[ ,"roads"]
  }
  detCovs[ ,2, ] <- detectors$detectors.df$detOtherSamples
  dimnames(detCovs) <- list( "detectors" = 1:n.detectors,
                             "covariates" = c("roads", "obs"),
                             "years" = years)
  detectors$covariates <- detCovs
  
  
  ##-- Merge with the habitat grid
  detectors$grid <- left_join(
    x = detectors$grid,
    y = detectors$detectors.df,
    by = "id")
  
  
  
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
  habitat$localObjects <- getLocalObjects(
    habitatMask = habitat$habitat.mx,
    coords = habitat$scaledCoords,
    dmax = habitat$maxDist/habitat$resolution,
    plot.check = F)
  
  ##-- Get local detectors
  detectors$localObjects <- getLocalObjects(
    habitatMask = habitat$habitat.mx,
    coords = detectors$scaledCoords,
    dmax = detectors$maxDist/detectors$resolution,
    plot.check = F)
  
  
  
  ## ------   5. SAVE STATE-SPACE CHARACTERISTICS -----
  
  save( habitat,
        file = file.path( working.dir, "data",
                          paste0("Habitat_bear_", DATE, ".RData")))
  
  save( detectors,
        file = file.path( working.dir,"data",
                          paste0("Detectors_bear_", DATE, ".RData")))
  
  
  
  
  ## ------   6. FILTER DATA -----
  ## ------     6.1. ALIVE DATA -----
  
  data.alive <- myFullData.sp$alive %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to months of interest
      Month %in% unlist(sampling.months),
      ##-- Subset to sex of interest
      Sex %in% sex,
      ##-- Filter data for space
      !is.na(as.numeric(sf::st_intersects(.,habitat.rWthBufferPol)))
    ) %>%
    AssignDetectors(
      myData = .,                
      myDetectors = detectors$main.detector.sp,
      mysubDetectors = detectors$detector.sp,
      radius = detectors$resolution)
  
  
  
  ## ------     6.2. DEAD RECOVERY DATA -----
  
  data.dead <- myFullData.sp$dead.recovery %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to sex of interest
      Sex %in% sex,
      ##-- Filter data for space
      !is.na(as.numeric(sf::st_intersects(.,habitat.rWthBufferPol)))
    ) %>% 
    AssignDetectors(
      myData = .,
      myDetectors = detectors$main.detector.sp,
      radius = detectors$resolution)
  
  
  
  ## ------     6.3. SAVE FILTERED DATA ----- 
  
  save( data.alive, data.dead,
        file = file.path( working.dir, "data",
                          paste0("FilteredData_bear_", DATE, ".RData")))
  
  
  
  ## ------   7. GENERATE DETECTION HISTORY ------
  
  for(thisSex in sex){
    
    message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
    
    ## ------     7.1. FILTER DATA BY SEX -----
    
    load(file.path( working.dir, "data",
                    paste0("FilteredData_bear_", DATE, ".RData")))
    
    data.alive$myData.sp <- data.alive$myData.sp %>%
      dplyr::filter(Sex %in% thisSex)
    
    data.dead <-  data.dead %>%
      dplyr::filter(Sex %in% thisSex)
    
    
    
    ## ------     7.2. GENERATE DETECTION HISTORY ARRAYS -----
    
    y.ar <- MakeY( myData = data.alive$myData.sp,
                   myDetectors = detectors$main.detector.sp,
                   method = "Binomial",
                   myData2 = data.dead,
                   myDetectors2 = detectors$main.detector.sp,
                   returnIdvector = TRUE)
    y.ar.ALIVE <- y.ar$y.ar
    dimnames(y.ar.ALIVE) <- dimnames(y.ar$y.ar)
    
    ##-- Project death to the next occasion
    y.ar.DEADProjected <- y.ar$y.ar2 
    y.ar.DEADProjected[] <- 0
    for(t in 2:n.years){
      y.ar.DEADProjected[ , ,t] <- y.ar$y.ar2[ , ,t-1]
    }
    
    ##-- Get dead recovery detector index
    y.ar.DEAD <- apply( y.ar.DEADProjected,
                        c(1,3),
                        function(x){
                          if(sum(x)>0){which(x>0)}else{0}
                        })
    dimnames(y.ar.DEAD) <- list( "id" = dimnames(y.ar$y.ar2)[[1]],
                                 "year" = dimnames(y.ar$y.ar2)[[3]])
    
    
    
    ## ------     7.3. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR -----
    
    distances <- list()
    for(t in 1:n.years){
      # print(paste("------ ", t ," -------", sep = "" ))
      distances[[t]] <- CheckDistanceDetections(
        y = y.ar.ALIVE[,,t], 
        detector.xy = detectors$detectors.df[ ,c("x","y")], 
        max.distance = detectors$maxDist,
        method = "pairwise",
        plot.check = F,
        verbose = F)
      
      # ##-- Plot individuals with detections further than the threshold distance
      # if(plot.check){
      #   par(mfrow = c(1,1))
      #   if(sum(distances[[t]]$y.flagged) > 0){
      #     affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
      #     count <- 0
      #     for(i in affected.ids){
      #       count <- count+1
      #       plot(st_geometry(studyArea), main = paste("t: ",t,"     i: ", names(affected.ids)[count], sep = ""))
      #       scalebar(2*myVars$DETECTIONS$maxDist, xy = c(800000,6700000), type = "bar", divs = 2, below = "km",
      #                label = c(0, myVars$DETECTIONS$maxDist/1000, myVars$DETECTIONS$maxDist/500), cex = 0.8, adj = c(0.5,-0.9))
      #       plot(st_geometry(COUNTRIES), add = T)
      #       plot(st_geometry(detectors$main.detector.sp), add = T, col = grey(0.8), cex = 0.3, pch = 19)
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
      #       plot(st_geometry(detectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0), ]), pch = 16, col = "red",add=T)
      #       
      #       tmp2 <- detectors$main.detector.sp[which(y.ar.ALIVE[i,,t] > 0 & distances[[t]]$y.flagged[i,] == 1), ]
      #       plot(st_geometry(tmp2), add = T, col = "blue", pch = 13, cex = 1.5, lwd = 1)
      #     }#i
      #   }#if
      # }#if plot.check
      # 
      ##-- REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
      y.ar.ALIVE[ , ,t] <- y.ar.ALIVE[ , ,t] * (1-distances[[t]]$y.flagged)
    }#t
    
    
    
    ## ------     7.4. AUGMENT DETECTION HISTORIES -----
    
    y.alive <- MakeAugmentation( 
      y = y.ar.ALIVE,
      aug.factor = data$aug.factor,
      replace.value = 0)
    
    y.dead.ar <- MakeAugmentation( 
      y = y.ar.DEAD,
      aug.factor = data$aug.factor,
      replace.value = 0)
    
    
    
    ## ------     7.5. TRANSFORM Y TO SPARSE MATRICES -----
    
    y.sparse <- nimbleSCR::getSparseY(y.alive)
    
    
    
    ## ------ IV. MODEL SETTING ------- 
    
    ## -----    1. NIMBLE CODE ------
  modelCode <- nimbleCode({
    
      ##----- SPATIAL PROCESS -----
    
      dmean ~ dunif(0,100)
      lambda <- 1/dmean
      
      betaDens  ~ dnorm(0.0,0.01)
      habIntensity[1:numHabWindows] <- exp(betaDens * denCounts[1:numHabWindows,1])
      sumHabIntensity <- sum(habIntensity[1:numHabWindows])
      logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
      logSumHabIntensity <- log(sumHabIntensity)
      
      for(i in 1:n.individuals){
        
        sxy[i,1:2,1] ~ dbernppAC(
          lowerCoords = lowerHabCoords[1:numHabWindows,1:2],
          upperCoords = upperHabCoords[1:numHabWindows,1:2],
          logIntensities = logHabIntensity[1:numHabWindows],
          logSumIntensity = logSumHabIntensity,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max)
        
        for(t in 2:n.years){
          
          sxy[i,1:2,t] ~ dbernppACmovement_exp(
            lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
            upperCoords = upperHabCoords[1:n.habwindows,1:2],
            s = sxy[i,1:2,t-1],
            lambda = lambda,
            baseIntensities = habIntensity[1:n.habwindows],
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            numGridRows = y.max,
            numGridCols = x.max,
            numWindows = n.habwindows)
        }#i
      }#i
      
      
      
      ##----- DEMOGRAPHIC PROCESS ----- 
      
      omeg1[1:2] ~ ddirch(alpha[1:2])   
      
      for(t in 1:n.years1){
        gamma[t] ~ dunif(0,1)
        phi[t] ~ dunif(0,1)
        
        omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)
        omega[2,1:3,t] <- c(0,phi[t],1-phi[t])
        omega[3,1:3,t] <- c(0,0,1)
      }#t
      
      for(i in 1:n.individuals){ 
        z[i,1] ~ dcat(omeg1[1:2]) 
        for(t in 1:n.years1){
          z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
        }#t								
      }#i
      
      
      
      ##----- DETECTION PROCESS -----
      
      for(t in 1:n.years){
        sigma[t] ~ dunif(0,4)
        
        ## Trap-response
        betaResponse[t] ~ dunif(-5,5)
        betaResponseOth[t] ~ dunif(-5,5)
        
        ## Covariates
        for(c in 1:n.covs){
          betaCovs[c,t] ~ dunif(-5,5)
        }#c
        
        for(c in 1:n.covsOth){
          betaCovsOth[c,t] ~ dunif(-5,5)
        }#c
        
        ## County/country-specific baselines
        for(c in 1:n.counties){
          p01[c,t] ~ dunif(0,1)
          p0[c,t] <- p01[c,t] * countyToggle[c,t]         ## toggle counties
        }#c
        
        for(c in 1:n.countries){
          p01Oth[c,t] ~ dunif(0,1)
          p0Oth[c,t] <- p01Oth[c,t] * countryToggle[c,t]  ## toggle counties
        }#c
      }#t 
      
      pResponse ~ dunif(0, 1)
      
      for(i in 1:n.individuals){
        
        detResponse[i,1] ~ dbern(pResponse)
        
        for(t in 1:n.years){
          
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
            p0 = p0[1:n.counties,t],
            detCounties[1:n.detectors],
            detCov = detCovs[1:n.detectors,t,1:n.covs],
            betaCov = betaCovs[1:n.covs,t],
            BetaResponse = betaResponse[t],
            detResponse = detResponse[i,t])
          
          y.aliveOth[i,1:nMaxDetectorsOth,t] ~ dbin_LESS_Cached_MultipleCovResponse(
            sxy = sxy[i,1:2,t],
            sigma = sigma[t],
            nbDetections = nbDetectionsOth[i,t],
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
      
      
      
      ##------ DERIVED PARAMETERS -----
      
      for(t in 1:n.years){
        for(i in 1:n.individuals){ 
          isAlive[i,t] <- (z[i,t] == 2) 
        }#i
        N[t] <- sum(isAlive[1:n.individuals,t])
      }#t
    })
    
    
    
    ## ------   2. NIMBLE CONSTANTS -----
    
    nimConstants <- list( n.individuals = dim(y.alive)[1], 
                          n.detectors = dim(detectors$scaledCoords)[1], 
                          n.detCovs = dim(detCovs)[2],
                          n.years = dim(y.alive)[3],
                          n.habwindows = nrow(habitat$habitat.df),
                          n.counties = max(detectors$detectors.df$counties),
                          y.max = nrow(habitat$localObjects$habitatGrid),
                          x.max = ncol(habitat$localObjects$habitatGrid),
                          lengthYCombined = y.sparse$lengthYCombined,
                          localDetNumMax = detectors$localObjects$numLocalIndicesMax,
                          localHabNumMax = habitat$localObjects$numLocalIndicesMax,
                          resizeFactor = habitat$localObjects$resizeFactor,
                          county = detectors$detectors.df$counties)
    
    
    
    ## ------   3. NIMBLE INITS -----
    
    ## ------     3.1. RECONSTRUCT z -----
    
    ##-- Reconstruct monthly z based on ALL detections and dead recoveries
    zMonths <- MakeZfromScratch( 
      data.alive = myFullData.sp$alive,
      data.dead = myFullData.sp$dead.recovery,
      samplingMonths = unlist(sampling.months))
    
    ##-- Subset to focal years
    zMonths <- zMonths[ , ,dimnames(zMonths)[[3]] %in% dimnames(y.alive)[[3]]]
    
    ##-- Subset to focal individuals
    zMonths <- zMonths[dimnames(zMonths)[[1]] %in% dimnames(y.alive)[[1]], , ]
    
    ##-- Augment zMonths
    zMonths <- MakeAugmentation( y = zMonths,
                                 aug.factor = data$aug.factor,
                                 replace.value = NA)
    
    ##-- Compress back to yearly z
    zYears <- apply(zMonths, c(1,3), function(x){
      if(any(x[1:length(unlist(sampling.months))] == 1, na.rm = T)){
        2
      } else {
        if(any(x[1:length(unlist(sampling.months))] >= 2, na.rm = T)){
          4 }
        else {NA}}})
    #table(zYears)
    
    z.data <- zYears
    allDead <- apply(z.data, 1, function(x)all(x==4))
    z.data[allDead, ] <- 1
    
    
    
    ## ------     3.2. GENERATE INITIAL z -----
    
    z.init <- t(apply(z.data, 1, function(zz){
      out <- zz
      out[] <- 1
      if(any(!is.na(zz))){
        ##-- Set all occasions after the last detection to 4
        if(sum(zz == 4, na.rm = T)<1){
          range.det <- range(which(!is.na(zz)))
          if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
          if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 4
        }
        ##-- Set occasion before recovery to 2  
        if(sum(zz == 3, na.rm = T)>0){
          reco.3 <- min(which(zz == 3))
          if(reco.3>1){zz[(reco.3-1)] <- 2}
        }
        ##-- Set occasions before first detection alive to 1
        if(sum(zz == 2, na.rm = T)>0){
          reco.alive <- min(which(zz == 2))
          if(reco.alive>1){zz[1:(reco.alive-1)] <- 1}
        }
        out[] <- zz
      }
      return(out)
    }))
    
    ##-- Set all known states to NA in the inits
    z.init[!is.na(z.data)] <- NA
    #table(z.init, useNA = "always")
    
    
    
    ## ------     3.3. GENERATE y.dead -----
    legal.mx <- do.call(rbind, lapply(dimnames(y.alive)[[1]], function(x){
      out <- rep(0,dim(z.data)[2])
      if(x %in% myFullData.sp$dead.recovery$Id[myFullData.sp$dead.recovery$legal == "yes"]) out <- rep(1,dim(z.data)[2])
      return(out)
    }))
    
    y.dead <- z.data
    y.dead[] <- ifelse(z.data %in% c(4) & legal.mx == 1,1,0)
    y.dead <- t(apply(y.dead, 1, function(x){
      out <- x
      out[] <- 0
      if(any(x==1)) out[min(which(x==1))] <- 1
      return(out)
    }))
    
    ##-- DISTINGUISH MORTALITY SOURCE IN z.data and z.init
    z.data[] <- ifelse(y.dead == 1, 3, z.data)
    # z.staggered[] <- ifelse(y.dead == 1, 3, z.staggered)
    
    
    
    ## ------     3.4. GENERATE sxy & sxy.init ARRAYS -----
    
    ##-- Provide sxy as data for recovered individuals
    s.data <- array(NA, c(dim(y.alive)[1], 2, n.years))
    for(i in 1:length(dimnames(y.ar.ALIVE)[[1]])){
      for(t in 2:n.years){
        if(z.data[i,t] %in% c(3,4)){
          temp <- myFullData.sp$dead.recovery[myFullData.sp$dead.recovery$Id %in% dimnames(y.ar.ALIVE)[[1]][i] &
                                                myFullData.sp$dead.recovery$Year %in% years[t-1], ]
          if(nrow(temp) > 0){
            if(!is.na(as.numeric(sf::st_intersects(temp, sf::st_as_sf(habitat$buffered.habitat.poly))))){
              if(raster::extract(habitat$habitat.r, temp)==0){
                ##-- Id could be dead in the spatial extent but can be outside the habitat 
                ##-- because in a cell < 49% habitat
                buff <- sf::st_buffer(temp, dist = habitat$resolution)
                inter <-sf::st_intersection(buff, habitat$buffered.habitat.poly)
                s.data[i, ,t] <- sf::st_coordinates(sf::st_sample( x = inter,
                                                                   size = 1,
                                                                   type = "random"))
              } else {
                s.data[i, ,t] <- sf::st_coordinates(temp)
              }
            }
          }
        }
      }
    }
    
    ##-- Rescale s.data
    dimnames(s.data) <- list("id" = dimnames(y.alive)[[1]],
                             "coord" = c("x","y"),
                             "year" = years)
    s.data <- nimbleSCR::scaleCoordsToHabitatGrid(
      coordsData = s.data,
      coordsHabitatGridCenter = habitat$habitat.df[ ,c("x","y")])$coordsDataScaled
    
    ##-- Generate s.init (taking into account s.data)
    s.init <- getInits.s( y = y.sparse$yCombined,
                          known.s = s.data,
                          trapCoords = detectors$scaledCoords,
                          lowerCoords = habitat$scaledLowerCoords,                    
                          upperCoords = habitat$scaledUpperCoords,
                          habitatGrid = habitat$localObjects$habitatGrid,
                          baseIntensities = rep(2,nimConstants$n.habwindows),
                          sd = 1)
    s.init[!is.na(s.data)] <- NA
    
    ##-- Check AC-movement distances
    test <- s.init
    test[is.na(test)] <- s.data[is.na(test)]
    dist <- sapply(2:dim(test)[3],
                   function(t){
                     sqrt((test[ ,1,t]-test[ ,1,t-1])^2 + (test[ ,2,t]-test[ ,2,t-1])^2)
                   })
    # hist(dist)
    # max(dist)
    
    ##-- Fix some annoying individuals
    # which(dist > habitat$maxDist/habitat$resolution, arr.ind = T)
    # dist[16, ]
    # y.sparse$yCombined[16,1, ]
    # s.init[16, ,3] <-  c(18.9,46.34)
    # tmp <- trunc(s.init[16, ,3]) + 1
    # habitat$localObjects$habitatGrid[tmp[2],tmp[1]]
    #
    # dist[222, ]
    # y.sparse$yCombined[222,1, ]
    # s.init[222, ,7] <- (s.init[222, ,8] + s.init[222, ,6])/2
    # tmp <- trunc(s.init[222, ,7]) + 1
    # habitat$localObjects$habitatGrid[tmp[2],tmp[1]]
    # 
    # ##-- Check AC-movement distances (again)
    # test <- s.init
    # test[is.na(test)] <- s.data[is.na(test)]
    # dist <- sapply(2:dim(test)[3],
    #                function(t){
    #                  sqrt((test[ ,1,t]-test[ ,1,t-1])^2 + (test[ ,2,t]-test[ ,2,t-1])^2)
    #                })
    # hist(dist)
    # max(dist)
    
    
    
    ## ------   4. NIMBLE DATA -----
    
    nimData <- list( z = z.data,   
                     sxy = s.data,
                     y.alive = y.sparse$yCombined,
                     y.dead = y.dead,
                     lowerHabCoords = habitat$scaledLowerCoords, 
                     upperHabCoords = habitat$scaledUpperCoords, 
                     habDens = cbind(habitat$habitat.df$dead.reco.trunc,
                                     habitat$habitat.df$skandObs.smooth),
                     habitatGrid = habitat$localObjects$habitatGrid,
                     localHabIndices = habitat$localObjects$localIndices,
                     localHabNum = habitat$localObjects$numLocalIndices,
                     alpha = c(1,1),
                     detCoords = detectors$scaledCoords,
                     size = detectors$detectors.df$size,
                     detCovs = detCovs,
                     localDetIndices = detectors$localObjects$localIndices,
                     localDetNum = detectors$localObjects$numLocalIndices)
    
    
    
    ## ------   5. NIMBLE PARAMETERS -----
    
    nimParams <- c("N", 
                   "lambda","dmean","betaDens",
                   "omeg1","gamma","phi",
                   "sigma","pResponse",
                   "p0","betaCovs","betaResponse",
                   "p0Oth","betaCovsOth","betaResponseOth")
    
    nimParams2 <- c("z", "sxy")
    
    
    
    ## ------   6. SAVE INPUTS ----- 
    
    for(c in 1:4){
      nimInits <- list( "sxy" = s.init,
                        "z" = z.init,
                        "tau" = 0.4,
                        "betaDens" = matrix(stats::runif(4,0,1),nrow = 2),
                        "omeg1" = c(0.7,0.3),
                        "gamma" = stats::runif(dim(y.alive)[3]-1,0.02,0.1),
                        "mhW" = stats::runif(dim(y.alive)[3]-1,0.1,0.3),
                        "mhH" = stats::runif(dim(y.alive)[3]-1,0.1,0.2),
                        "p0" = array(stats::runif(nimConstants$n.counties*n.years,0,0.1),
                                     c(nimConstants$n.counties,n.years)),
                        "betaDet" = stats::runif(2,-0.5,0.5),
                        "sigma" = stats::runif(1,0.1,0.4))
      
      save( modelCode,
            nimData,
            nimConstants,
            nimParams,
            nimParams2,
            nimInits,
            detCounties.original,
            file = file.path( working.dir, "nimbleInFiles", thisSex,
                              paste0("nimbleInput_", DATE, "_", c, ".RData")))
    }#c
    
  }#thisSex
  
  
  
  
  ## ------   8. RETURN IMPORTANT INFOS FOR REPORT ------
  return(list( SPECIES = "Brown bear",
               engSpecies = "bear",
               YEARS = years,
               SEX = sex,
               DATE = DATE))
}

