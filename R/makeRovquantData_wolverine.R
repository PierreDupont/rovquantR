#' @title RovQuant OPSCR wolverine data preparation.
#'
#' @description
#' \code{makeRovquantData_wolverine} formats the available wolverine data for the OPSCR analysis using nimble and nimbleSCR.
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
#' @param rename.list A \code{list}.
#' 
#' 
#' @return 
#' A \code{html} report summarizing the data preparation process
#' Additional \code{.png} images that can be reused somewhere else.
#'
#' @author Pierre Dupont
#' 
#' @import dplyr
#' @import raster
#' @import sf 
#' @importFrom adehabitatHR estUDm2spixdf kernelUD
#' @importFrom fasterize fasterize
#' @importFrom grDevices grey
#' @importFrom nimbleSCR getSparseY scaleCoordsToHabitatGrid getLocalObjects
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
  sex = c("female","male"),
  aug.factor = 0.8,
  sampling.months = list(12,1:6),
  
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 60000,
  max.move.dist = 250000,
  
  ##-- detectors
  detector.res = 10000,
  subdetector.res = 2000,
  max.det.dist = 84000,
  resize.factor = 1,
  
  ##-- Miscellanious
  rename.list = NULL
){
  
  ## ------ 0. BASIC SET-UP ------
  
  ##-- Set default values for the wolverine model
  if(is.null(aug.factor)){aug.factor <- 0.8}
  if(is.null(sampling.months)){sampling.months <- list(12,1:6)}
  if(is.null(habitat.res)){habitat.res <- 20000} 
  if(is.null(buffer.size)){buffer.size <- 60000}
  if(is.null(max.move.dist)){max.move.dist <- 250000}
  if(is.null(detector.res)){detector.res <- 10000}
  if(is.null(subdetector.res)){subdetector.res <- 2000}
  if(is.null(max.det.dist)){max.det.dist <- 40000}
  if(is.null(resize.factor)){resize.factor <- 1}
  if(is.null(rename.list)) {
    if(!exists("r.list.internalWolf")) stop("Default 'rename.list' not available")
    rename.list <- r.list.internalWolf
  }  
  
  ##-- Set up list of Habitat characteristics
  habitat <- list( resolution = habitat.res,
                   buffer = buffer.size,
                   maxDist = max.move.dist)
  
  ##-- Set up list of Detectors characteristics
  detectors <- list( resolution = detector.res,
                     resolution.sub = subdetector.res,
                     maxDist = max.det.dist,
                     resize.factor = resize.factor)
  
  ##-- Set up list of Data characteristics
  DATA <- list( sex = sex,
                aug.factor = aug.factor,
                sampling.months = sampling.months)
  
  
  
  ##----------------------------------------------------------------------------
  
  ## ------ I. LOAD AND SELECT DATA ------
  
  ## ------   1. HABITAT DATA -----
  
  ##-- Load pre-defined habitat rasters and shapefiles
  data("habitatRasters", envir = environment()) 
  data("REGIONS", envir = environment())
  data("COUNTRIES", envir = environment())
  data("GLOBALMAP", envir = environment())
  
  ##-- Disaggregate habitat raster to the desired resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)
  
  ##-- Merge counties for practical reasons
  COUNTIES_AGGREGATED <- REGIONS %>%
    dplyr::mutate(id = case_when(
      ## Sweden
      county %in% c("Blekinge","Dalarna","Gävleborg","Gotland","Halland","Jämtland",
                    "Jönköping","Kalmar","Kronoberg","Örebro","Östergötland","Skåne",
                    "Södermanland","Stockholm","Uppsala","Värmland","Västernorrland",
                    "Västmanland","Västra Götaland") ~ 3,
      county %in% c("Norrbotten") ~ 1,
      county %in% c("Västerbotten") ~ 2,
      ## Norway
      county %in% c("Agder","Akershus","Buskerud","Innlandet","Møre og Romsdal",
                    "Oppland","Oslo","Østfold","Rogaland","Vestland","Telemark",
                    "Vestfold") ~ 4,
      county %in% c("Trøndelag") ~ 5,
      county %in% c("Nordland") ~ 7,
      county %in% c("Troms") ~ 8,
      county %in% c("Finnmark") ~ 6)) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize()
  
  ##-- Get Norrbotten borders
  COUNTIESNorrbotten <- REGIONS %>%
    dplyr::filter(county %in% "Norrbotten") %>%
    dplyr::group_by(county) %>%
    dplyr::summarize()
  
  
  
  ## ------   2. NGS DATA -----
  
  ##-- Extract date from the last cleaned data file
  DATE <- getMostRecent( 
    path = file.path(working.dir, "data"),
    pattern = "CleanData_wolverine")
  
  
  ##-- Load the most recent clean wolverine data from RovBase
  myFullData.sp <- readMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_wolverine",
    extension = ".RData")
  
  ##-- List years
  if(is.null(years)){
    years <- sort(unique(c(myFullData.sp$alive$Year,
                           myFullData.sp$dead.recovery$Year)))
  }
  DATA$years <- years
  n.years <- length(years)
  
  ##-- Filter NGS samples for dates
  myFullData.sp$alive <- myFullData.sp$alive %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to monitoring period
      Month %in% unlist(sampling.months), 
      ##-- Subset to samples collected in Norway and Sweden
      myFullData.sp$alive$Country_sf %in% c("(N)","(S)")) 
    
  ##-- Filter NGS samples outside the GLOBAL MAP
  ##-- [PD]: should be removed!!
  ##-- This removes 125 samples falling into lakes
  myStudyArea <- GLOBALMAP %>% 
    dplyr::filter(ISO %in% c("SWE","NOR")) %>%
    mutate(id = 1) %>% 
    group_by(id) %>% 
    summarize() 
  
  myFullData.sp$alive <- myFullData.sp$alive %>% 
    dplyr::filter(!is.na(as.numeric(st_intersects(., myStudyArea))))
  
  
  ##-- Filter out detections in Norrbotten except in 2016:18 and after 2023
  ##-- list years with or without sampling in Norrbotten
  yearsSampledNorrb <- c(2016:2018,2023:as.numeric(format(Sys.Date(), "%Y")))
  yearsNotSampled <- years[!years %in% yearsSampledNorrb]
  whichYearsNotSampled <- which(years %in% yearsNotSampled)
  ##-- Identify detections collected in Norrbotten 
  is.Norr <- as.numeric(st_intersects(myFullData.sp$alive, COUNTIESNorrbotten))
  ##-- Filter out detections in Norrbotten in years without sampling
  myFullData.sp$alive <- myFullData.sp$alive %>%
    dplyr::filter(!(Year %in% yearsNotSampled & is.Norr %in% 1))
  
  
  ##-- Filter Dead recoveries for dates
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
    ##-- Subset to years of interest
    dplyr::filter(Year %in% years)
  
  
  
  
  ##----------------------------------------------------------------------------
  
  ## ------ II. CREATE OPSCR DATA ------
  
  ## ------   1. GENERATE HABITAT ------
  
  message("Preparing habitat characteristics... ")
  
  ## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------
  
  ##-- Determine study area based on NGS detections
  ##-- Buffer NGS detections and cut to Swedish and Norwegian borders
  studyArea <- myFullData.sp$alive %>%
    sf::st_buffer(., dist = habitat$buffer * 1.4) %>%
    dplyr::mutate(id = 1) %>%
    dplyr::group_by(id) %>% 
    dplyr::summarize() %>% 
    sf::st_intersection(., COUNTRIES) %>%
    sf::st_as_sf()
  
  ##-- Get study area extent
  studyArea.extent <- st_bbox(extent(studyArea))
  st_crs(studyArea.extent) <- st_crs(studyArea)
  
  ##-- Make habitat from predefined Scandinavian raster of suitable habitat
  habitat <- makeHabitatFromRaster(
    poly = studyArea,
    habitat.r = habRaster,
    buffer = habitat$buffer,
    plot.check = FALSE) %>%
    append(habitat,.)
  
  ##-- Retrieve number of habitat windows 
  isHab <- habitat$habitat.r[] == 1
  n.habWindows <- habitat$n.habWindows <- sum(isHab)
  habitat$habitat.df <- cbind.data.frame(
    "id" = 1:n.habWindows,
    "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
    "y" = raster::coordinates(habitat$habitat.r)[isHab,2])
  
  ##-- Make a spatial grid from polygon
  habitat$grid <- sf::st_as_sf( stars::st_as_stars(habitat$habitat.r), 
                                as_points = FALSE,
                                merge = FALSE) %>%
    dplyr::filter( Habitat %in% 1) %>%
    dplyr::mutate( id = 1:nrow(.)) %>%
    sf::st_set_crs( .,value = sf::st_crs(habitat$buffered.habitat.poly))
  
  ##-- Study area grid from habitat raster
  habitat.rWthBufferPol <- sf::st_as_sf( 
    stars::st_as_stars(habitat$habitat.rWthBuffer), 
    as_points = FALSE,
    merge = TRUE) %>%
    dplyr::filter(Habitat %in% 1)
  
  
  
  ## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------
  
  ## ------       1.2.1. DEN COUNTS ------
  
  ##-- Load the last DEN COUNT data file
  DEN <- readMostRecent( 
    path = data.dir,
    extension = ".csv",
    pattern = "DEN_COUNTS") %>%
    sf::st_as_sf(., coords = c("UTM33_X", "UTM33_Y")) %>%
    sf::st_set_crs(value = sf::st_crs(myFullData.sp$alive)) %>%
    dplyr::mutate(id = 1)
  
  DEN.r <- raster::raster(
    adehabitatHR::estUDm2spixdf(
      adehabitatHR::kernelUD( as(DEN[ ,"id"], "Spatial"),
                              h = 30000,
                              grid = as(habitat$habitat.r, 'SpatialPixels'))))
  
  ##-- EXTRACT COVARIATES
  denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
  denCounts <- as.vector(round(scale(denCounts), digits = 2))
  
  ##-- Put into "nimble2SCR" format
  habitat$habitat.df <- cbind.data.frame( habitat$habitat.df,
                                          "den.counts" = denCounts)
  
  ##-- Merge with the habitat grid
  habitat$grid <- dplyr::left_join(
    x = habitat$grid,
    y = habitat$habitat.df,
    by = "id")
  
  
  
  ## ------   2. GENERATE DETECTORS -----
  
  message("Preparing detectors characteristics... ")
  
  ## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----
  
  ##-- Generate NGS detectors based on the study area 
  detectors$subdetectors.r <- raster::disaggregate(
    habitat$habitat.rWthBuffer,
    fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)
  
  ##-- Generate NGS detectors based on the raster of sub-detectors
  detectors <- makeSearchGrid( 
    data = detectors$subdetectors.r,
    resolution = detectors$detResolution,
    div = (detectors$resolution/detectors$resolution.sub)^2,
    plot = FALSE) %>%
    append(detectors, .)
  
  ##-- Format detector locations & number of trials per detector
  detectors$detectors.df <- cbind.data.frame(
    "id" = 1:nrow(detectors$main.detector.sp),
    "x" = sf::st_coordinates(detectors$main.detector.sp)[ ,1],
    "y" = sf::st_coordinates(detectors$main.detector.sp)[ ,2],
    "size" = detectors$main.detector.sp$count)
  
  ##-- Generate detector raster 
  detectors$raster <- raster::rasterFromXYZ(
    cbind( detectors$detectors.df[ ,c("x","y")],
           Detector = rep(1,nrow(detectors$main.detector.sp))))
  
  ##-- Make a spatial grid from detector raster
  detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
    x = detectors$raster,
    fun = function(x){x>0})) %>%
    sf::st_set_crs(.,value = sf::st_crs(studyArea)) %>%
    dplyr::mutate( id = 1:nrow(.)) 
  
  ##-- Extract numbers of detectors
  n.detectors <- detectors$n.detectors <- dim(detectors$main.detector.sp)[1]
  
  ##-- Identify detectors in Norrbotten 
  COUNTIESAroundNorrbotten <- REGIONS %>%
    dplyr::group_by(county) %>%
    dplyr::summarize() %>%
    dplyr::filter(county %in% c("Norrbotten","Troms","Västerbotten","Nordland","Finnmark")) %>% 
    sf::st_simplify( dTolerance = 500)
  
  ##-- Create an index of detectors in Norrbotten
  distDetsCounties <- sf::st_distance(
    detectors$main.detector.sp,
    COUNTIESAroundNorrbotten,
    byid = T)
  detsNorrbotten <- which(apply(distDetsCounties, 1, which.min) == 3)
  
  
  
  ## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----
  
  ## ------       2.2.1. EXTRACT COUNTIES ------
  
  ##-- Extract closest county for each detector
  detCounties <- detectors$main.detector.sp %>%
    sf::st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
    apply(., 1, function(x) which.min(x))
  
  ##-- Put into "nimble2SCR" shape
  detectors$detectors.df$counties <- detCounties
  
  ##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten 
  ##-- in years without sampling
  countyToggle <- matrix(1, nrow = max(detCounties), ncol = n.years)
  for(t in whichYearsNotSampled){
    countyToggle[1,t] <- 0
  }
  
  
  
  ## ------       2.2.2. EXTRACT COUNTRIES ------
  
  ##-- Extract closest country for each detector
  detCountries <- detectors$main.detector.sp %>%
    sf::st_distance(., COUNTRIES, by_element = F) %>%
    apply(., 1, function(x) which.min(x)) %>%
    as.factor(.) %>%
    as.numeric(.)
  
  ##-- Put into "nimble2SCR" format
  detectors$detectors.df$countries <- detCountries
  
  ##-- Turn into a matrix with years in columns
  detCountries <- matrix( detCountries,
                          nrow = length(detCountries),
                          ncol = n.years)
  
  ##-- Add another category to detCountry if in Norrbotten, to turnoff detection to 0 there. 
  for(t in whichYearsNotSampled){
    detCountries[detCounties %in% 1,t] <- 3
  }#t  
  
  ##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten in years without sampling
  countryToggle <- matrix(1, nrow = max(detCountries), ncol = n.years)
  for(t in whichYearsNotSampled){
    countryToggle[3,t] <- 0
  }
  
  
  
  ## ------       2.2.3. EXTRACT GPS TRACKS LENGTHS ------
  
  message("Cleaning GPS tracks... ")
  
  ##-- Combine all GPS tracks
  TRACKS <- readTracks( data.dir = data.dir,
                        years = years,
                        sampling.months = sampling.months) %>%
    ##-- Filter out irrelevant tracks
    dplyr::filter(Jerv == "1") %>%
    ##-- Subset to the study area
    st_intersection(., st_as_sfc(studyArea.extent)) 
  
  ##-- Extract length of GPS search track per detector grid cell
  detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
  for(t in 1:n.years){
    intersection <- TRACKS %>%
      dplyr::filter(Year == years[t]) %>%
      sf::st_intersection(detectors$grid, .) %>%
      dplyr::mutate(LEN = st_length(.)) %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(transect_L = sum(LEN)) 
    detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
    print(t)
  }#t
  
  ##-- Put into "nimble2SCR" format
  colnames(detTracks) <- paste0("tracks.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detTracks)
  
  
  
  ## ------       2.2.4. EXTRACT DISTANCES TO ROADS ------
  
  ##-- Load map of distance to roads (1km resolution)
  DistAllRoads <- raster::raster(file.path(data.dir,"Roads/MinDistAllRoads1km.tif"))
  
  ##-- Fasterize to remove values that fall in the sea
  r <- fasterize::fasterize(sf::st_as_sf(REGIONS), DistAllRoads)
  r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  DistAllRoads <- r
  DistAllRoads <- raster::crop(DistAllRoads, studyArea)
  rm(list = c("r"))
  
  ##-- AGGREGATE TO MATCH THE DETECTORS RESOLUTION
  DistAllRoads <- raster::aggregate( 
    x = DistAllRoads,
    fact = detectors$resolution/raster::res(DistAllRoads),
    fun = mean)
  
  ##-- EXTRACT ROAD DISTANCE FOR EACH DETECTOR
  detRoads <- raster::extract(DistAllRoads, detectors$main.detector.sp)
  
  ##-- if NA returns the average value of the cells within 15000m 
  isna <- which(is.na(detRoads))
  tmp <- raster::extract( DistAllRoads,
                          detectors$main.detector.sp[isna, ],
                          buffer = 15000,
                          fun = mean,
                          na.rm = T)
  detRoads[isna] <- tmp
  detRoads <- round(scale(detRoads), digits = 2)
  
  ##-- Put into "nimble2SCR" format
  detectors$detectors.df$roads <- detRoads
  
  
  
  ## ------       2.2.5. EXTRACT DAYS OF SNOW ------
  
  ##-- Average snow from December to June (the official monitoring period for Norway&Sweden)
  SNOW <- stack(file.path(data.dir,"Snow/AverageSnowCoverModisSeason2014_2025_Wolverine.tif"))
  
  # ##-- RENAME THE LAYERS
  # names(SNOW) <- paste(years, (years) + 1, sep = "_")
  
  ##-- SELECT SNOW DATA CORRESPONDING TO THE MONITORING PERIOD
  SNOW <- SNOW[[paste("X", years, "_", years + 1, sep = "")]]
  SNOW <- raster::crop(SNOW, c(0,40,55,75))
  
  ##-- EXTRACT SNOW 
  detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = n.years)
  det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
  detSnow[ ,1:n.years] <- raster::extract(SNOW, det.sptransf)
  
  ##-- if NA returns the average value of the cells within 15km 
  isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
  tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                          buffer = 15000, fun = mean, na.rm = T)
  detSnow[isna,1:n.years] <- tmp
  
  ##-- Put into "nimble2SCR" format
  colnames(detSnow) <- paste0("snow.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detSnow)
  
  
  
  ## ------       2.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------
  
  ## ------         2.2.6.1. SKANDOBS ------
  
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
                   species = stringi::stri_trans_general(species, "Latin-ASCII"),
                   monitoring.season = ifelse( month > unlist(sampling.months)[1],
                                               year, year-1)) %>%
    ## [PD] the version above corresponds to last year's analysis ("54.Cleaned2025TestPDScript.R")
    ## It is wrong because it does nothing. Only samples collected AFTER December get their year changed (i.e. no samples at all).
    ## Below is the correct version, similar to what we do in cleanRovbaseData():
    ## monitoring.season = ifelse( month < unlist(sampling.months)[1],
    ##                             year - 1, year)) %>%
    ##-- Filter based on monitoring season
    dplyr::filter( month %in% unlist(sampling.months)) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf(., coords = c("longitude","latitude")) %>%
    sf::st_set_crs(., value = "EPSG:4326") %>%
    sf::st_transform(., sf::st_crs(REGIONS)) %>%
    sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------         2.2.6.2. ROVBASE ------
  
  ##-- Process Rovbase observations (all species)
  rovbaseObs <- readMultiples( 
    path = file.path(data.dir, "AllSamples"),
    extension = ".xlsx") %>%
    ##-- Rename columns to facilitate manipulation
    dplyr::rename(., any_of(rename.list)) %>%
    ##-- Extract important info (e.g. month, year, country of collection)
    dplyr::mutate(
      ##-- Turn potential factors into characters 
      across(where(is.factor), as.character),
      ##-- Deal with Scandinavian characters
      Species = stringi::stri_trans_general(Species, "Latin-ASCII"),
      Sample_type = translateForeignCharacters(data = Sample_type),
      ##-- Deal with dates
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      year = as.numeric(format(Date,"%Y")),
      month = as.numeric(format(Date,"%m")),
      monitoring.season = ifelse( month > unlist(sampling.months)[1],
                                  year, year-1)) %>%
    ## [PD] the version above corresponds to last year's analysis ("54.Cleaned2025TestPDScript.R")
    ## It is wrong because it does nothing. Only samples collected AFTER December get their year changed (i.e. no samples at all).
    ## Below is the correct version, similar to what we do in cleanRovbaseData():
    ## monitoring.season = ifelse( month < unlist(sampling.months)[1],
    ##                             year - 1, year)) %>%
    ##-- Filter out unusable samples
    dplyr::filter( 
      ##-- Filter out samples without coordinates,...
      !is.na(East_UTM33),
      ##-- ...based on species
      # Species %in% c("Bjorn","Fjellrev","Gaupe","Hund","Jerv","Rodrev","Ulv"),
      ##-- ...based on sample type
      Sample_type %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
                          "Saliv/Spytt","Loepeblod","Blod","Vev"),
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months),
      ##-- ... if sample was from the focal species and successfully genotyped 
      !(Species %in% "Jerv" & !is.na(Id))) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(. , sf::st_crs(REGIONS)) %>%
    ##-- Filter based on space 
    sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------         2.2.6.3. COMBINE ROVBASE & SKANDOBS ------
  
  ##-- Rasterize at the detector level
  r.list <- lapply(years, function(y){
    ##-- Rasterize Skandobs observations 
    sk.r <- raster::rasterize(
      skandObs[skandObs$monitoring.season %in% y, 1],
      detectors$raster,
      fun = "count")[[1]]
    sk.r[is.na(sk.r[])] <- 0
    ##-- Set cells outside detector area to NA
    sk.r[!detectors$raster[ ]%in% 1] <- NA
    ##-- Turn into binary raster
    sk.r1 <- sk.r
    sk.r1[sk.r1[]>0] <- 1
    
    ##-- Rasterize Rovbase observations 
    rb.r <- raster::rasterize(
      rovbaseObs[rovbaseObs$monitoring.season %in% y, 1],
      detectors$raster,
      fun = "count")[[1]]
    rb.r[is.na(rb.r[])] <- 0
    ##-- Set cells outside detector area to NA
    rb.r[! detectors$raster[ ]%in% 1] <- NA
    ##-- Turn into binary raster
    rb.r1 <- rb.r
    rb.r1[rb.r1[]>0] <- 1
    ##-- Store in a list
    list(sk.r, sk.r1, rb.r, rb.r1)
  })
  
  ##-- Store in raster bricks
  r.skandObsBinary <- brick(lapply(r.list,function(x) x[[2]]))
  r.skandObsContinuous <- brick(lapply(r.list,function(x) x[[1]]))
  r.rovbaseBinary <- brick(lapply(r.list,function(x) x[[4]]))
  r.rovbaseContinuous <- brick(lapply(r.list,function(x) x[[3]]))
  
  ##-- Combine both rasters
  r.SkandObsRovbaseBinary <- r.rovbaseBinary + r.skandObsBinary
  for(t in 1:n.years){
    r.SkandObsRovbaseBinary[[t]][r.SkandObsRovbaseBinary[[t]][] > 1] <- 1
  }
  
  
  ##-- Plot check
  # if(plot.check){
  #   plot(st_geometry(habitat.rWthBufferPol))
  #   plot(st_geometry(skandObs), col = "red", add = T)
  #   
  #   ##-- Summary SkanbObs
  #   pdf(file = file.path(working.dir, "figures/skandObs.pdf"), width = 10)
  #   barplot(table(skandObs$monitoring.season ))
  #   barplot(table(skandObs$month ), xlab = "Months")
  #   barplot(table(skandObs$species))
  #   
  #   ##-- Maps
  #   par(mar = c(0,0,2,0))
  #   for(t in 1:n.years){
  #     plot( st_geometry(studyArea), main = years[t])
  #     plot( st_geometry(skandObs[skandObs$monitoring.season %in% years[t], ]),
  #           pch = 16, col = "red", cex = 0.1, add = T)
  #   }
  #   dev.off()
  #   
  #   # pdf(file = file.path(working.dir, "figures/mapStructuredOthers.pdf"))
  #   # for(t in 1:n.years){
  #   #   tmpOthers <- data.alive[data.alive$Year %in% years[t] &
  #   #                                          !data.alive$structured, ]
  #   #   tmpStruct <- data.alive[data.alive$Year %in% years[t] &
  #   #                                          data.alive$structured, ]
  #   # 
  #   #   par(mfrow=c(2,2),mar=c(0,0,5,0))
  #   #   plot(r.OtherSamplesBinary[[t]], main=paste(years[t],"\n Rovbase Samples Structured"), box=F, axes=F)
  #   #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  #   #   plot(r.OtherSamplesBinary[[t]],main=paste(years[t],"\n Rovbase Samples Opportunistic"), box=F, axes=F)
  #   #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
  #   # 
  #   #   plot(r.skandObsSamplesBinary[[t]], main=paste(years[t]), box=F, axes=F)
  #   #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
  #   #   plot(r.skandObsSamplesBinary[[t]],main=paste(years[t],"\n SkandObs Opportunistic"), box=F, axes=F)
  #   #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
  #   # }
  #   # dev.off()
  #   
  #   for(t in 1:n.years){
  #     par(mfrow=c(1,3),mar=c(0,0,5,0))
  #     plot(r.rovbaseBinary[[t]],main=years[t])
  #     plot(r.skandObsBinary[[t]])
  #     plot(r.SkandObsRovbaseBinary[[t]])
  #   }#t
  # }
  
  
  
  ## ------         2.2.6.4. IDENTIFY CELLS WITH HAIR TRAPS AS OPPORTUNISTIC ------
  
  ##-- IDENTIFY HAIR SAMPLES
  tmpHair <- myFullData.sp$alive %>% dplyr::filter(hairTrap)
  
  ##-- MANUALLY FIND THE HAIR SAMPLES & COLOR THE CELL.
  tmpyr <- unique(tmpHair$Year)
  for(i in 1:length(tmpyr)){
    t <- which(years %in% tmpyr[i])
    whereHair <- raster::extract( r.SkandObsRovbaseBinary[[t]],
                                  tmpHair[tmpHair$Year %in% tmpyr[i], ],
                                  cellnumbers = T)
    r.SkandObsRovbaseBinary[[t]][whereHair[ ,1]] <- 1
  }#t
  
  
  
  ## ------         2.2.6.5. ASSIGN THE COVARIATE ------
  
  detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
  detOtherSamples[ ,1:n.years] <- raster::extract( r.SkandObsRovbaseBinary,
                                                   detectors$main.detector.sp)
  colnames(detOtherSamples) <- paste0("detOtherSamples.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detOtherSamples)
  
  
  
  ## ------       2.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------
  
  detSnow <- round(scale(detSnow), digits = 2)
  detRoads <- round(scale(detRoads), digits = 2)
  detTracks <- round(scale(detTracks), digits = 2)
  
  detCovs <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 2))
  detCovs[,,1] <- detTracks
  detCovs[,,2] <- detSnow
  dimnames(detCovs) <- list( "detectors" = 1:n.detectors,
                             "years" = years,
                             "covariates" = c("tracks", "snow"))
  
  detCovsOth <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 3))
  detCovsOth[,,1] <- detSnow
  detCovsOth[,,2] <- matrix(detRoads,length(detRoads),n.years)
  detCovsOth[,,3] <- detOtherSamples
  dimnames(detCovsOth) <- list( "detectors" = 1:n.detectors,
                                "years" = years,
                                "covariates" = c("snow", "roads", "obs"))
  
  ##-- Store in the detector list                          
  detectors$covariates <- detCovs
  detectors$covariates.others <- detCovsOth
  
  ##-- Merge with the detector grid
  detectors$grid <- dplyr::left_join(
    x = detectors$grid,
    y = detectors$detectors.df,
    by = "id")
  
  
  
  ## ------   3. RESCALE COORDINATES ------
  
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
  
  
  
  ## ------   4. CREATE LOCAL OBJECTS -----
  
  ##-- Get local detectors
  detectors$localObjects <- getLocalObjects(
    habitatMask = habitat$habitat.mx,
    coords = detectors$scaledCoords,
    dmax = detectors$maxDist/habitat$resolution,
    resizeFactor = detectors$resize.factor,
    plot.check = F)
  
  
  
  ## ------   5. SAVE STATE-SPACE CHARACTERISTICS -----
  
  save( habitat,
        file = file.path( working.dir, "data",
                          paste0("Habitat_wolverine_", DATE, ".RData")))
  
  save( detectors,
        file = file.path( working.dir,"data",
                          paste0("Detectors_wolverine_", DATE, ".RData")))
  
  
  
  ## ------   6. FILTER DATA -----
  
  ## ------     6.1. ALIVE DATA -----
  
  data.alive <- myFullData.sp$alive %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to months of interest
      Month %in% unlist(sampling.months),
      ##-- Subset to sex of interest
      Sex %in% sex) %>%
    ##-- Filter based on space 
    sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------     6.2. DEAD RECOVERY DATA -----
  
  data.dead <- myFullData.sp$dead.recovery %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to sex of interest
      Sex %in% sex) %>%
    ##-- Filter based on space 
    sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------     6.3. FILTER OUT DETECTIONS IN NORRBOTTEN EXCEPT IN 2016:18 and after 2023 ------
  
  # ##-- Get Norrbotten borders
  # COUNTIESNorrbotten <- REGIONS %>%
  #   dplyr::filter(county %in% "Norrbotten") %>%
  #   dplyr::group_by(county) %>%
  #   dplyr::summarize()
  # 
  # ##-- Identify detections collected in Norrbotten 
  # is.Norr <- as.numeric(st_intersects(data.alive, COUNTIESNorrbotten))
  # 
  # # ##-- Check how many detections are removed per year
  # # table(data.alive$Year[data.alive$Year %in% yearsNotSampled & is.Norr %in% 1])
  # # sum(data.alive$Year %in% yearsNotSampled & is.Norr %in% 1)
  # 
  # ##-- Filter out detections in Norrbotten in years without sampling
  # data.alive <- data.alive %>%
  #   dplyr::filter(!(Year %in% yearsNotSampled & is.Norr %in% 1))
  
  
  
  ## ------     6.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------
  
  ## ------       6.3.1. ASSIGN SAMPLES TO GPS SEARCH TRACKS ------
  
  message("Assigning DNA samples to GPS tracks... ")
  message("This can take several minutes... ")
  
  data.alive <- assignSearchTracks(
    data = data.alive,
    tracks = TRACKS)
  
  # ##-- SAVE FOR FASTER LOADING
  # save(myFilteredData.sp, file = file.path(working.dir, "data/myFilteredData.sp.RData"))
  # load(file.path(working.dir, "data/myFilteredData.sp.RData"))
  
  
  
  ## ------       6.3.2. ASSIGN SAMPLES TO OPPORTUNISTIC OR STRUCTURED ------
  
  distanceThreshold <- 500
  
  ##-- Identify samples from structured and opportunistic sampling
  data.alive <- data.alive %>%
    dplyr::mutate(
      ##-- Collector column was replaced by two columns, merging them now...
      Collector_role = ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
      ##-- Identify samples collected during structured sampling 
      structured = Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
        !is.na(trackID) &
        trackDist <= distanceThreshold & 
        !hairTrap)
  
  
  
  ## ------       6.3.3. PLOT CHECKS ------
  
  # if(plot.check){
  #   ##-- Barplot of structured vs. opportunistic samples
  #   pdf(file = file.path(working.dir, "figures/DetectionsStructuredOppBarplot.pdf"))
  #   par(mfrow = c(2,1), mar = c(4,4,3,2))
  #   barplot( rbind(table(data.alive$Year[data.alive$structured]),
  #                  table(data.alive$Year[!data.alive$structured])),
  #            beside = T,
  #            ylim = c(0,3000),
  #            col = c(grey(0.2),grey(0.8)),
  #            ylab = "Number of samples")
  #   abline(h = seq(0, 3000, by = 500),
  #          lty = 2, col = grey(0.8))
  #   title(main = "500m threshold")
  #   legend("topleft", fill = c(grey(0.2),grey(0.8)),
  #          legend = c("Structured","Other"))
  #   
  #   ##-- Barplot of structured vs. opportunistic samples (threshold = 2000m)
  #   structured2000 <- data.alive$Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
  #     !is.na(data.alive$trackID) &
  #     data.alive$trackDist <= 2000 & 
  #     !data.alive$hairTrap
  #   barplot( rbind(table(data.alive$Year[structured2000]),
  #                  table(data.alive$Year[!structured2000])),
  #            beside = T,
  #            ylim = c(0,3000),
  #            col = c(grey(0.2),grey(0.8)),
  #            ylab = "Number of samples")
  #   abline(h=seq(0,3000,by=500),lty=2,col=grey(0.8))
  #   title(main="2000m threshold")
  #   legend("topleft",fill=c(grey(0.2),grey(0.8)),
  #          legend = c("Structured","Other"))
  #   dev.off()
  #   
  #   ##-- CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO"
  #   tmp <- data.alive %>%
  #     filter(Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"))
  #   
  #   ## plot check
  #   pdf(file = file.path(working.dir, "figures/DetectionsStructuredOpp.pdf"))
  #   for(t in 1:n.years){
  #     par(mar = c(0,0,3,0), mfrow = c(1,3))
  #     
  #     ##-- samples with tracks
  #     tmpTracks <- tmp %>%
  #       filter( Year %in% years[t],
  #               !is.na(trackID),
  #               trackDist <= 500)
  #     plot( st_geometry(studyArea), col = "gray60", main = "Structured with track")
  #     plot( st_geometry(tmpTracks),
  #           pch = 21, col = "black",
  #           cex = 1, bg = "red", add = T)
  #     
  #     ##-- samples without tracks
  #     tmpNoTracks <- tmp %>%
  #       filter( Year %in% years[t],
  #               is.na(trackID) | trackDist > 500)
  #     plot( st_geometry(studyArea), col = "gray60", main = "Structured without track")
  #     plot( st_geometry(tmpNoTracks),
  #           pch = 21, col = "black",
  #           cex = 1, bg = "blue", add = T)
  #     
  #     ##-- Other samples
  #     tmpOpp <- data.alive %>%
  #       filter(!Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),
  #              Year%in% years[t])
  #     plot( st_geometry(studyArea), col = "gray60", main = "Other samples")
  #     plot( st_geometry(tmpOpp),
  #           pch = 21, col = "black",
  #           cex = 1, bg = "green", add = T)
  #     mtext(years[t], adj = -0.8, padj = 1)
  #   }#t
  #   
  #   ##-- Number of samples collected per year
  #   tab <- table(tmp$Year, tmp$trackID, useNA ="always")
  #   barplot( tab[ ,which(is.na(colnames(tab)))]/rowSums(tab),
  #            main = "% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track")
  #   dev.off()
  #   
  #   ##-- plot check
  #   pdf( file = file.path(working.dir, "figures/OverallDetectionsDeadRecoveries.pdf"))
  #   plot( st_geometry(GLOBALMAP))
  #   plot( st_geometry(studyArea), add = T)
  #   plot( st_geometry(myFullData.sp$alive),
  #         pch = 16, col = "red", cex = 0.3, add = T)
  #   plot( st_geometry(myFullData.sp$dead.recovery),
  #         pch = 16, col = "blue", cex = 0.3, add = T)
  #   mtext(paste("Live detections", nrow(myFullData.sp$alive),
  #               "; ID:", nrow(unique(myFullData.sp$alive$Id))),
  #         line = +1)
  #   mtext(paste("Dead recovery:", nrow(myFullData.sp$dead.recovery)))
  #   dev.off()
  # }
  
  
  
  ## ------     6.4. SEPARATE MORTALITY CAUSES ------
  
  # ##-- Identify legal mortality causes
  # MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
  # whichLegalCauses <- unlist(lapply(c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske"),
  #                                   function(x)grep(x,MortalityNames)))
  # legalCauses <- MortalityNames[whichLegalCauses]
  # 
  # ##-- Identify legal dead recoveries based on mortality causes
  # data.dead <- data.dead %>%
  #   mutate(legal = Death_cause %in% legalCauses)
  # 
  # 
  # ##-- Plot check
  # if(plot.check){
  #   # par(mfrow = c(1,3))
  #   # for(t in 1:n.years){
  #   #   ## DEAD RECOVERIES TOTAL
  #   #   tempTotal <- data.dead[data.dead$Year == years[t], ]
  #   #   NGS_TabTotal <- table(tempTotal$Country_sf)
  #   #   ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country_sf), 2, function(x) sum(x>0))
  #   #   ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
  #   #   tempIn <- data.dead[data.dead$Year == years[t], ]
  #   #   NGS_TabIn <- table(tempIn$Country_sf)
  #   #   ID_TabIn <- apply(table(tempIn$Id, tempIn$Country_sf), 2, function(x) sum(x>0))
  #   #   ## PLOT NGS SAMPLES
  #   #   plot(st_geometry(GLOBALMAP), col="gray80")
  #   #   plot(st_geometry(studyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
  #   #   plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
  #   #   ## ADD NUMBER OF NGS samples and IDs per COUNTRY
  #   #   graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
  #   #   graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
  #   #   ## ADD OVERALL NUMBERS
  #   #   mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
  #   #   mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
  #   #   # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
  #   # }#t
  #    
  #   
  #   ##-- Plot temporal trends
  #   
  #   ##-- Number of detections
  #   pdf(file = file.path(working.dir, "figures/TRENDDetections.pdf"))
  #   
  #   temp <- unique(data.alive[ ,c("Year","Country_sf","DNAID")])
  #   tab_Country.Year <- table(temp$Year, temp$Country_sf)
  #   country.colors <- c("goldenrod1","goldenrod3")
  #   
  #   par(mfrow = c(1,1), mar = c(5,5,5,5))
  #   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
  #   lines(tab_Country.Year[,"(N)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
  #   lines(tab_Country.Year[,"(S)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
  #   legend("bottomright",c("N","S"), fill=country.colors)
  #   
  #   ##-- Number of IDs detected
  #   temp <- table(data.alive$Year,data.alive$Country_sf,data.alive$Id)
  #   tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
  #   country.colors <- c("goldenrod1","goldenrod3")
  #   par(mfrow = c(1,1), mar = c(5,5,5,5))
  #   plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
  #   lines(tab_Country.Year1[,"(N)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
  #   lines(tab_Country.Year1[,"(S)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
  #   legend("bottomright",c("N","S"), fill=country.colors)
  #   
  #   ##-- Average number of detection per detected ID
  #   tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
  #   par(mfrow = c(1,1), mar = c(5,5,5,5))
  #   plot(-10,
  #        xlim = range(as.numeric(row.names(tab_Country.Year2))),
  #        ylim = c(0,max(tab_Country.Year2)),
  #        ylab = "Average Number of detections", xlab="Years")
  #   lines(tab_Country.Year2[ ,"(N)"] ~ as.numeric(row.names(tab_Country.Year2)),
  #         col = country.colors[1], lwd = 2, pch = 16, type = "b")
  #   lines(tab_Country.Year2[ ,"(S)"] ~ as.numeric(row.names(tab_Country.Year2)),
  #         col = country.colors[2], lwd = 2, pch = 16, type = "b")
  #   legend("bottomright", c("N","S"), fill = country.colors)
  #   
  #   ##-- Dead recoveries
  #   temp <- unique(data.dead[,c("Year","Country_sf","Id")])
  #   tab_Country.Year <- table(temp$Year, temp$Country_sf)
  #   country.colors <- c("goldenrod1","goldenrod3")
  #   
  #   par(mfrow = c(1,1), mar = c(5,5,5,5))
  #   plot(-10,
  #        xlim = range(as.numeric(row.names(tab_Country.Year))),
  #        ylim = c(0,max(tab_Country.Year)),
  #        ylab = "N Id Dead recovered", xlab="Years")
  #   lines(tab_Country.Year[ ,"(N)"] ~ as.numeric(row.names(tab_Country.Year)),
  #         col = country.colors[1], lwd = 2, pch = 16, type = "b")
  #   lines(tab_Country.Year[ ,"(S)"] ~ as.numeric(row.names(tab_Country.Year)),
  #         col = country.colors[2], lwd = 2, pch = 16, type = "b")
  #   legend("topright", c("N","S"), fill = country.colors)
  #   dev.off()
  # }
  
  
  
  ## ------     6.5. ASSIGN SAMPLES TO DETECTORS -----
  
  ##-- ALL SAMPLES
  data.alive <- assignDetectors( 
    data = data.alive,                
    detectors = detectors$main.detector.sp,
    subDetectors = detectors$detector.sp,
    radius = detectors$resolution)
  
  ##-- DEAD RECOVERY
  data.dead <- assignDetectors( 
    data = data.dead,
    detectors = detectors$main.detector.sp,
    radius = detectors$resolution)
  
  
  
  ## ------     6.6. FIX DETECTIONS IN NORRBOTTEN ------
  
  ### MAKE SURE THAT INDIVIDUALS DETECTED OUTSIDE OF NORRBOTTEN DO NOT GET 
  ### ASSIGNED TO A DETECTOR IN NORRBOTTEN IN YEARS WERE THERE IS NO SAMPLING.
  ### FIND THE CASES WHERE IT HAPPENS AND ASSIGN THEM THE CLOSEST DETECTOR OUTSIDE
  ### OF NORRBOTTEN
  
  ##-- Identify sub-detectors inside Norrbotten
  subDetsNorrbotten <- which( detectors$detector.sp$main.cell.id %in% 
                                detectors$main.detector.sp$main.cell.id[detsNorrbotten])
  
  ##-- Identify which detections are assigned to Norrbotten in years not sampled
  whichDets <- which(!data.alive$data.sp$Year %in% yearsSampledNorrb &
                       data.alive$data.sp$Detector %in% detsNorrbotten)
  
  ##-- Loop over flagged detections & assign them to closest detector outside Norrbotten
  for(i in 1:length(whichDets)){
    tmp <- data.alive$data.sp[whichDets[i], ]
    ## Calculate distance to all sub-detectors
    dist <- st_distance(tmp, detectors$detector.sp)
    ## Artificially increase distance for detectors in Norrbotten 
    dist[ ,subDetsNorrbotten] <- 500000
    ## Assign detection to closest sub-detector outside Norrbotten
    data.alive$data.sp$sub.detector[whichDets[i]] <- which.min(dist[1, ])
    ## Assign detection to the corresponding main detector outside Norrbotten
    thisDet <- detectors$detector.sp$main.cell.id[which.min(dist[1, ])]
    data.alive$data.sp$Detector[whichDets[i]] <- which(detectors$main.detector.sp$main.cell.id == thisDet)
  }#i
  
  # ##-- SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
  # sum(data.alive$data.sp$sub.detector[!data.alive$data.sp$Year %in% yearsSampledNorrb] %in% subDetsNorrbotten)
  # sum(data.alive$data.sp$Detector[!data.alive$data.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)
  
  
  
  ## ------     6.7. PLOT NGS & DEAD RECOVERY MAPS ----- 
  
  ##-- layout
  L <- n.years
  if(L < 6){ nrows <- 1 } else{
    if(L < 13){ nrows <- 2 } else {
      if(L < 22){ nrows <- 3 } else {
        if(L < 33){ nrows <- 4 } else { nrows <- 5 }}}}
  ncols <- ceiling(L/nrows)
  
  
  ##-- NGS maps
  grDevices::png(filename = file.path(working.dir, "figures/NGS_TimeSeries.png"),
                 width = ncols*2, height = nrows*4,
                 units = "in", pointsize = 12,
                 res = 300, bg = NA)
  ##-- layout
  mx <- matrix(NA, nrow = nrows*2, ncol =  (ncols*2)+1)
  for(r in 1:nrows){
    mx[r*2-1, ] <- c(1,rep(1:ncols, each = 2)) + (r-1)*ncols
    mx[r*2, ] <- c(rep(1:ncols, each = 2),ncols) + (r-1)*ncols
  }#r
  nf <- graphics::layout(mx,
                         widths = c(rep(1,ncol(mx))),
                         heights = rep(1,2))
  par(mar = c(0,0,0,0))
  
  for(t in 1:length(years)){
    ##-- Plot maps
    plot( sf::st_geometry(COUNTRIES), border = NA, col = c("gray80","gray60"))
    try(
      plot( sf::st_geometry(data.alive$data.sp[data.alive$data.sp$Year == years[t], ]), add = TRUE, col = "orange", pch = 3),
      silent = TRUE)
    plot( sf::st_geometry(COUNTRIES), border = "gray20", col = NA, add = TRUE)
    
    ##-- Add year
    graphics::mtext(text = years[t],
                    side = 1, line = -18,
                    adj = 0.18, cex = 1.2)
  }#t
  dev.off()
  
  
  ##-- Dead recoveries maps
  grDevices::png(filename = file.path(working.dir, "figures/DEAD_TimeSeries.png"),
                 width = ncols*2, height = nrows*4,
                 units = "in", pointsize = 12,
                 res = 300, bg = NA)
  ##-- layout
  mx <- matrix(NA, nrow = nrows*2, ncol =  (ncols*2)+1)
  for(r in 1:nrows){
    mx[r*2-1, ] <- c(1,rep(1:ncols, each = 2)) + (r-1)*ncols
    mx[r*2, ] <- c(rep(1:ncols, each = 2), ncols) + (r-1)*ncols
  }#r
  nf <- graphics::layout(mx,
                         widths = c(rep(1,ncol(mx))),
                         heights = rep(1,2))
  par(mar = c(0,0,0,0))
  
  for(t in 1:length(years)){
    ##-- Plot maps
    plot( sf::st_geometry(COUNTRIES), border = NA, col = c("gray80","gray60"))
    try( plot( sf::st_geometry(data.dead[data.dead$Year == years[t] & 
                                           data.dead$Legal, ]),
               add = TRUE, 
               col = "slateblue1",
               pch = 3),
         silent = TRUE)
    try( plot( sf::st_geometry(data.dead[data.dead$Year == years[t], ]),
               add = TRUE,
               col = "slateblue4",
               pch = 3),
         silent = TRUE)
    plot( sf::st_geometry(COUNTRIES),
          border = "gray20",
          col = NA,
          add = TRUE)
    
    ##-- Add year
    graphics::mtext(text = years[t],
                    side = 1, line = -18,
                    adj = 0.18, cex = 1.2)
  }#t
  dev.off()
  
  
  
  ## ------     6.8. SAVE FILTERED DATA ----- 
  
  save( data.alive,
        data.dead,
        file = file.path( working.dir, "data",
                          paste0("FilteredData_wolverine_", DATE, ".RData")))
  
  
  
  ## ------   7. GENERATE DETECTION HISTORY ------
  
  for(thisSex in sex){
    
    message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
    
    ## ------     7.1. FILTER DATA BY SEX -----
    
    load(file.path( working.dir, "data",
                    paste0("FilteredData_wolverine_", DATE, ".RData")))
    
    data.alive$data.sp <- data.alive$data.sp %>%
      dplyr::filter(Sex %in% thisSex)
    
    data.dead <- data.dead %>%
      dplyr::filter(Sex %in% thisSex)
    
    
    
    ## ------     7.2. GENERATE DETECTION HISTORY ARRAYS -----
    
    ##-- ALL SAMPLES
    y.ar <- makeY( data = data.alive$data.sp,
                   detectors = detectors$main.detector.sp,
                   method = "Binomial",
                   data2 = data.dead,
                   detectors2 = detectors$main.detector.sp,
                   returnIdvector = TRUE)
    
    ##-- STRUCTURED
    y.arStruc <- makeY( data = data.alive$data.sp[data.alive$data.sp$structured, ],
                        detectors = detectors$main.detector.sp,
                        method = "Binomial",
                        data2 = data.dead,
                        detectors2 = detectors$main.detector.sp,
                        returnIdvector = TRUE)
    
    ##-- OTHERS
    y.arOth <- makeY( data = data.alive$data.sp[!data.alive$data.sp$structured, ],
                      detectors = detectors$main.detector.sp,
                      method = "Binomial",
                      data2 = data.dead,
                      detectors2 = detectors$main.detector.sp,
                      returnIdvector = TRUE)
    
    ##-- Make sure all detection arrays have the same dimensions
    y.ar.ALIVEOthers <- y.ar.ALIVEStructured <- array( 0, 
                                                       dim = dim(y.ar$y.ar),
                                                       dimnames = dimnames(y.ar$y.ar))
    ##-- Fill in the y arrays
    y.ar.ALIVEOthers[dimnames(y.arOth$y.ar)[[1]], , ] <- y.arOth$y.ar
    y.ar.ALIVEStructured[dimnames(y.arStruc$y.ar)[[1]], , ] <- y.arStruc$y.ar
    
    ##-- Project death to the next occasion
    y.ar.DEADProjected <- y.ar$y.ar2 
    y.ar.DEADProjected[] <- 0
    for(t in 2:n.years){ y.ar.DEADProjected[ , ,t] <- y.ar$y.ar2[ , ,t-1] }
    
    ##-- Get dead recovery detector index
    y.ar.DEAD <- apply( y.ar.DEADProjected,
                        c(1,3),
                        function(x){
                          if(sum(x)>0){which(x>0)}else{0}
                        })
    dimnames(y.ar.DEAD) <- list( "id" = dimnames(y.ar$y.ar2)[[1]],
                                 "year" = dimnames(y.ar$y.ar2)[[3]])
    
    ##-- Create binary dead recovery histories (0: not recovered ; 1: recovered dead)
    y.ar.DEAD <- apply(y.ar.DEADProjected, c(1,3), function(x){as.numeric(sum(x)>0)})
    dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
    y.ar.DEAD[y.ar.DEAD > 0] <- 1
    
    
    
    ## ------     7.3. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------
    
    distances <- list()
    for(t in 1:n.years){
      
      ##-- Identify detections further than maxDist
      print(paste0("------ ", t ," -------"))
      distances[[t]] <- checkDistanceDetections( 
        y = y.ar$y.ar[ , ,t], 
        detector.xy = detectors$detectors.df[ ,c("x","y")], 
        max.distance = detectors$maxDist,
        method = "pairwise",
        plot.check = F)
      
      ##-- If any detection flagged
      if(sum(distances[[t]]$y.flagged) > 0){
        
        ##-- Remove detections that are further than the threshold
        #y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
        y.ar.ALIVEOthers[ , ,t] <- y.ar.ALIVEOthers[ , ,t] * (1-distances[[t]]$y.flagged)
        y.ar.ALIVEStructured[ , ,t] <- y.ar.ALIVEStructured[ , ,t] * (1-distances[[t]]$y.flagged)
        
        ##-- Remove detections also in data.alive$data.sp to run getSInits later
        affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
        idd <- names(affected.ids)
        for(i in 1:length(idd)){
          detIds <- which(distances[[t]]$y.flagged[idd[i], ] > 0)
          data.alive$data.sp <- data.alive$data.sp %>%
            dplyr::filter(!(Id %in% idd[i] & Detector %in% detIds & Year %in% years[t]))
        }#i
      }#if
      
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
      #       tmp <- data.alive[data.alive$Id == dimnames(y.ar.ALIVE)[[1]][i] &
      #                                        data.alive$Year == years[t], ]
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
    }#t
    
    
    
    ## ------     7.4. GENERATE INDIVIDUAL-LEVEL COVARIATES ------
    
    ##-- Make matrix of previous capture indicator
    detResponse <- makeTrapResponseCov(
      data = myFullData.sp$alive,
      data.dead = myFullData.sp$dead.recovery)
    
    ##-- Subset to focal years
    detResponse <- detResponse[ ,dimnames(detResponse)[[2]] %in% dimnames(y.ar$y.ar)[[3]]]
    
    ##-- Subset to focal individuals
    detResponse <- detResponse[dimnames(detResponse)[[1]] %in% dimnames(y.ar$y.ar)[[1]], ]
    
    
    
    ## ------     7.5. AUGMENT DETECTION HISTORIES -----
    
    ##-- DATA ARRAYS
    y.alive <- makeAugmentation( y = y.ar$y.ar,
                                 aug.factor = aug.factor,
                                 replace.value = 0)
    
    y.dead <- makeAugmentation( y = y.ar.DEAD,
                                aug.factor = aug.factor,
                                replace.value = 0)
    
    y.aliveOthers <- makeAugmentation( y = y.ar.ALIVEOthers,
                                       aug.factor = aug.factor,
                                       replace.value = 0)
    
    y.aliveStructured <- makeAugmentation( y = y.ar.ALIVEStructured,
                                           aug.factor = aug.factor, 
                                           replace.value = 0)
    
    ##-- INDIVIDUAL COVARIATES
    detResponse <- makeAugmentation( y = detResponse,
                                     aug.factor = aug.factor,
                                     replace.value = 0)
    
    ##-- Set first detection for augmented individuals to NA
    detResponse[rownames(detResponse) %in% "Augmented",1]  <- NA
    
    
    
    ## ------     7.6. TRANSFORM Y TO SPARSE MATRICES ------
    
    ##-- STRUCTURED
    y.sparse <- nimbleSCR::getSparseY(y.aliveStructured)
    
    ##-- OTHER
    y.sparseOth <- nimbleSCR::getSparseY(y.aliveOthers)
    
    
    
    ## ------ IV. MODEL SETTING ------- 
    
    ## ------   1. NIMBLE MODEL DEFINITION ------
    
    modelCode <- nimbleCode({
      
      ## ----- SPATIAL PROCESS ------ 
      
      dmean ~ dunif(0,100)
      lambda <- 1/dmean
      
      betaDens ~ dnorm(0.0,0.01)
      habIntensity[1:n.habWindows] <- exp(betaDens * denCounts[1:n.habWindows])
      sumHabIntensity <- sum(habIntensity[1:n.habWindows])
      logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
      logSumHabIntensity <- log(sumHabIntensity)
      
      for(i in 1:n.individuals){
        sxy[i, 1:2, 1] ~ dbernppAC(
          lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
          upperCoords = upperHabCoords[1:n.habWindows,1:2],
          logIntensities = logHabIntensity[1:n.habWindows],
          logSumIntensity = logSumHabIntensity,
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max)
        
        for(t in 2:n.years){
          sxy[i, 1:2, t] ~ dbernppACmovement_exp(
            lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
            upperCoords = upperHabCoords[1:n.habWindows,1:2],
            s = sxy[i,1:2,t-1],
            lambda = lambda,
            baseIntensities = habIntensity[1:n.habWindows],
            habitatGrid =  habitatGrid[1:y.max,1:x.max],
            numGridRows = y.max,
            numGridCols = x.max,
            numWindows = n.habWindows)
        }#i  
      }#t
      
      
      ## ----- DEMOGRAPHIC PROCESS -----
      
      omeg1[1:2] ~ ddirch(alpha[1:2])   
      
      for(t in 1:(n.years-1)){
        gamma[t] ~ dunif(0,1)
        phi[t] ~ dunif(0,1)
        
        omega[1,1:3,t] <- c(1-gamma[t],gamma[t],0)  # "Unborn"
        omega[2,1:3,t] <- c(0,phi[t],1-phi[t])      # "Alive"
        omega[3,1:3,t] <- c(0,0,1)                  # "Dead"
      }#t
      
      for(i in 1:n.individuals){ 
        z[i,1] ~ dcat(omeg1[1:2]) 
        for(t in 1:(n.years-1)){
          z[i,t+1] ~ dcat(omega[z[i,t],1:3,t]) 
        }#t 								
      }#i
      
      
      ## ----- DETECTION PROCESS ----- 
      
      pResponse ~ dunif(0, 1)
      
      for(t in 1:n.years){
        
        sigma[t] ~ dunif(0,4)
        betaResponse[t] ~ dunif(-5,5)
        betaResponseOth[t] ~ dunif(-5,5)
        
        for(c in 1:n.covs){
          betaCovs[c,t] ~ dunif(-5,5)
        }#c
        
        for(c in 1:n.covs.Oth){
          betaCovsOth[c,t] ~ dunif(-5,5)
        }#c
        
        for(c in 1:n.counties){
          p01[c,t] ~ dunif(0,1)
          p0[c,t] <- p01[c,t] * countyToggle[c,t]         ## toggle counties
        }#c  
        
        for(c in 1:n.countries){
          p01Oth[c,t] ~ dunif(0,1)
          p0Oth[c,t] <- p01Oth[c,t] * countryToggle[c,t] ## toggle countries
        }#c  
      }#t
      
      
      for(i in 1:n.individuals){
        
        detResponse[i,1] ~ dbern(pResponse)
        
        for(t in 1:n.years){
          
          # y.alive[i,1:maxDetNums,t] ~ dbin_LESS_Cached_MultipleCovResponse(  
          #   sxy = sxy[i,1:2,t],
          #   sigma = sigma[t],
          #   nbDetections = nbDetections[i,t],
          #   yDets = detIndices[i,1:maxDetNums,t],
          #   detector.xy = detector.xy[1:n.detectors,1:2],
          #   trials = trials[1:n.detectors],
          #   detectorIndex = detectorIndex[1:n.habWindows,1:numLocalIndicesMax],
          #   nDetectorsLESS = nDetectorsLESS[1:n.habWindows],
          #   ResizeFactor = resizeFactor,
          #   maxNBDets = maxDetNums,
          #   habitatID = habitatIDDet[1:y.max,1:x.max],
          #   indicator = isAlive[i,t],
          #   p0State = p0[1:n.counties,t],
          #   detCountries = detCounties[1:n.detectors],
          #   detCov = detCovs[1:n.detectors,t,1:n.covs],
          #   betaCov = betaCovs[1:n.covs,t],
          #   BetaResponse = betaResponse[t],
          #   detResponse = detResponse[i,t])
          
          y[i,1:maxDetNums,t] ~ dbinomLocal_normalWolverine(
            detNums = detNums[i,t],
            detIndices = detIndices[i,1:maxDetNums,t],
            size = size[1:n.detectors],
            p0 = p0[1:n.counties,t],
            sigma = sigma[t],
            s = sxy[i,1:2,t],
            trapCoords = detector.xy[1:n.detectors,1:2],
            localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
            localTrapsNum = localDetNum[1:n.habWindows],
            resizeFactor = resizeFactor,
            lengthYCombined = maxDetNums,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            indicator = isAlive[i,t],
            trapCovsIntercept = detCounties[1:n.detectors],
            trapCovs = detCovs[1:n.detectors,t,1:n.covs],
            trapBetas = betaCovs[1:n.covs,t],
            indBeta = betaResponse[t],
            indCov = detResponse[i,t])
          
          y.Oth[i,1:maxDetNumsOth,t] ~ dbinomLocal_normalWolverine(
            detNums = detNumsOth[i,t],
            detIndices = detIndicesOth[i,1:maxDetNumsOth,t],
            size = size[1:n.detectors],
            p0 = p0Oth[1:n.countries,t],
            sigma = sigma[t],
            s = sxy[i,1:2,t],
            trapCoords = detector.xy[1:n.detectors,1:2],
            localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
            localTrapsNum = localDetNum[1:n.habWindows],
            resizeFactor = resizeFactor,
            lengthYCombined = maxDetNumsOth,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            indicator = isAlive[i,t],
            trapCovsIntercept = detCountries[1:n.detectors,t],
            trapCovs = detCovsOth[1:n.detectors,t,1:n.covs.Oth],
            trapBetas = betaCovsOth[1:n.covs.Oth,t],
            indBeta = betaResponseOth[t],
            indCov = detResponse[i,t])
        }#t
      }#i
      
      
      
      ## ----- DERIVED PARAMETERS -----
      
      for(i in 1:n.individuals){ 
        for(t in 1:n.years){
          isAlive[i,t] <- (z[i,t] == 2) 
        }#t
      }#i
      
      for(t in 1:n.years){
        N[t] <- sum(isAlive[1:n.individuals,t])
      }#t
    })
    
    
    
    ## ------   2. NIMBLE CONSTANTS ------
    
    nimConstants <- list( 
      n.individuals = dim(y.sparse$y)[1],
      n.detectors = nrow(detectors$scaledCoords),
      n.habWindows = nrow(habitat$scaledLowerCoords),
      n.years = dim(y.sparse$y)[3], 
      n.covs = dim(detCovs)[3],
      n.covs.Oth = dim(detCovsOth)[3],
      n.countries = max(detCountries),
      n.counties = max(detCounties),
      countyToggle = countyToggle,
      countryToggle = countryToggle,
      resizeFactor = detectors$localObjects$resizeFactor,
      y.max = dim(detectors$localObjects$habitatGrid)[1],
      x.max = dim(detectors$localObjects$habitatGrid)[2],
      numLocalIndicesMax = detectors$localObjects$numLocalIndicesMax,
      maxDetNums = y.sparse$maxDetNums,
      maxDetNumsOth = y.sparseOth$maxDetNums)
    #,
    #lengthYCombined = y.sparse$lengthYCombined,
    #lengthYCombined.Oth = y.sparseOth$lengthYCombined)
    
    
    
    ## ------   3. NIMBLE DATA ------
    
    ## ------     3.1. GENERATE KNOWN z ------
    
    ##-- Set all individuals alive to 2 between first and last detection
    z <- apply(y.alive, c(1,3), function(x)ifelse(any(x>0), 2, NA))
    z <- t(apply(z, 1, function(zz){
      if(any(!is.na(zz))){
        range.det <- range(which(!is.na(zz)))
        zz[range.det[1]:range.det[2]] <- 2
      }
      return(zz)
    }))
    
    
    
    ## ------     3.2. LIST DATA ------
    
    nimData <- list( 
      z = z,   
      y = y.sparse$y,
      detIndices = y.sparse$detIndices,
      detNums = y.sparse$detNums,
      y.Oth = y.sparseOth$y, 
      detIndicesOth = y.sparseOth$detIndices,
      detNumsOth = y.sparseOth$detNums,
      lowerHabCoords = as.matrix(habitat$scaledLowerCoords), 
      upperHabCoords = as.matrix(habitat$scaledUpperCoords), 
      detCounties = detCounties,
      detCountries = detCountries,
      detCovs = detCovs,
      detCovsOth = detCovsOth,
      detResponse = detResponse,
      denCounts = denCounts,
      localDetIndices = detectors$localObjects$localIndices,
      localDetNum = detectors$localObjects$numLocalIndices,
      habitatGrid = detectors$localObjects$habitatGrid,
      size = detectors$detectors.df$size,
      alpha = rep(1,2),
      detector.xy = as.matrix(detectors$scaledCoords))
    #habitatGrid = habIDCells.mx)
    
    
    
    ## ------   4. NIMBLE INITS ------
    
    ## ------     4.1. GENERATE INITIAL z ------
    
    ##-- Set z to 1 before first detection and 3 after last detection
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
    
    ##-- Set initial values to NA when individual state is known
    z.init <- ifelse(!is.na(z), NA, z.init)
    
    
    
    ## ------     4.2. LATENT VARIABLE DET RESPONSE ------
    
    detResponse.inits <- nimData$detResponse
    detResponse.inits[is.na(detResponse.inits)] <- rbinom(sum(is.na(detResponse.inits)),1,0.5)
    detResponse.inits[!is.na(detResponse)] <- NA
    
    
    
    ## ------     4.3. GENERATE INITIAL sxy ------
    
    ##-- sxy
    AllDets <- data.dead[ ,c("Id","Year")] %>% 
      ##-- Project death to the next year
      mutate(Year = Year + 1) %>%
      ##-- Remove dead reco occuring the last year (not used)
      filter(!Year %in% max(Year)) %>%
      ##-- Combine with detections alive
      rbind(.,data.alive$data.sp[ ,c("Id","Year")]) %>%
      ##-- Add coordinates
      mutate("x" = st_coordinates(.)[ ,1],
             "y" = st_coordinates(.)[ ,2]) %>%
      as.data.frame()
    
    ##-- Rescale detections
    AllDets <- scaleCoordsToHabitatGrid(
      coordsData = AllDets,
      coordsHabitatGridCenter = habitat$habitat.xy,
      scaleToGrid =T )$coordsDataScaled
    
    ##-- Generate initial sxy values
    sxy.init <- getSInits( AllDetections = AllDets[,c("Id","Year","x","y")],
                           Id.vector = y.ar$Id.vector,
                           idAugmented = which(rownames(z) %in% "Augmented"),
                           lowerCoords = nimData$lowerHabCoords,
                           upperCoords = nimData$upperHabCoords,
                           habitatGrid = nimData$habitatGrid,
                           intensity = NULL,
                           sd = 4,
                           movementMethod = "dbernppACmovement_normal")
    
    ##-- An extreme number of decimals may cause a number to appear as an integer
    ##-- to Nimble, and then coincide with habitat window boundaries
    sxy.init <- round(sxy.init, 4)
    
    
    
    ## ------   5. NIMBLE PARAMETERS ------
    
    nimParams <- c( "N", "lambda", "dmean", "betaDens",
                    "omeg1", "gamma", "phi",
                    "pResponse", "sigma",
                    "p0", "betaResponse", "betaCovs",
                    "p0Oth", "betaResponseOth", "betaCovsOth")
    
    nimParams2 <- c("z", "sxy")
    
    
    
    ## ------   6. SAVE INPUTS ----- 
    
    for(c in 1:4){
      
      nimInits <- list(
        "sxy" = sxy.init,
        "z" = z.init,
        "dmean" = stats::runif(1, 0, 10),
        "betaDens" = stats::runif(1, -0.1, 0.1),
        "omeg1" = c(0.5, 0.5),
        "gamma" = stats::runif(dim(y.alive)[3]-1, 0, 1),
        "phi" = stats::runif(dim(y.alive)[3]-1, 0.1, 0.3),
        "pResponse" = stats::runif(1, 0.4, 0.5),
        "detResponse" = detResponse.inits,
        "sigma" = stats::runif(n.years, 1, 4),
        "p01" = array(stats::runif(18, 0, 0.2),
                      c(nimConstants$n.counties, dim(y.alive)[3])),
        "betaResponse" = stats::runif(dim(y.alive)[3], -0.1, 0.1),
        "betaCovs" = array(stats::runif(dim(detCovs)[3], -0.1, 0.1),
                           c(dim(detCovs)[3], n.years)),
        "p01Oth" = array(stats::runif(18, 0, 0.2),
                         c(nimConstants$n.countries, dim(y.alive)[3])),
        "betaResponseOth" = stats::runif(dim(y.alive)[3], -0.1, 0.1),
        "betaCovsOth" = array(stats::runif(dim(detCovsOth)[3], -0.1, 0.1),
                              c(dim(detCovsOth)[3], n.years))) 
      
      save( modelCode,
            nimData,
            nimConstants,
            nimParams,
            nimParams2,
            nimInits,
            file = file.path( working.dir, "nimbleInFiles", thisSex,
                              paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
    }#c
  }#thisSex
  
  
  
  ##----------------------------------------------------------------------------
  
  ## ------ III. CREATE SCR DATA ------
  
  ## ------   1. SCR NIMBLE MODEL DEFINITION ------
  
  modelCode1 <- nimbleCode({
    
    ##------ SPATIAL PROCESS ------ 
    
    betaDens  ~ dnorm(0.0,0.01)
    habIntensity[1:n.habWindows] <- exp(betaDens * denCounts[1:n.habWindows])
    sumHabIntensity <- sum(habIntensity[1:n.habWindows])
    logHabIntensity[1:n.habWindows] <- log(habIntensity[1:n.habWindows])
    logSumHabIntensity <- log(sumHabIntensity)
    
    for(i in 1:n.individuals){
      sxy[i,1:2] ~ dbernppAC(
        lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
        upperCoords = upperHabCoords[1:n.habWindows,1:2],
        logIntensities = logHabIntensity[1:n.habWindows],
        logSumIntensity = logSumHabIntensity,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        numGridRows = y.max,
        numGridCols = x.max)
    }#i
    
    
    ##----- DEMOGRAPHIC PROCESS -----
    
    psi ~ dunif(0, 1)
    
    for(i in 1:n.individuals){ 
      z[i] ~ dbern(psi)	
    }#i
    
    
    ##----- DETECTION PROCESS -----
    
    sigma ~ dunif(0,4)
    
    pResponse ~ dunif(0, 1)
    betaResponse ~ dunif(-5,5)
    betaResponseOth ~ dunif(-5,5)
    
    for(c in 1:n.covs){
      betaCovs[c] ~ dunif(-5,5)
    }
    
    for(c in 1:n.covs.Oth){
      betaCovsOth[c] ~ dunif(-5,5)
    }
    
    for(c in 1:n.counties){
      p01[c] ~ dunif(0,1)
      p0[c] <- p01[c] * countyToggle[c]         ## toggle counties
    }#c  
    
    for(c in 1:n.countries){
      p01Oth[c] ~ dunif(0,1)
      p0Oth[c] <- p01Oth[c] * countryToggle[c]  ## toggle countries
    }#c  
    
    for(i in 1:n.individuals){
      
      detResponse[i] ~ dbern(pResponse)
      
      y[i,1:maxDetNums] ~ dbinomLocal_normalWolverine(
        detNums = detNums[i],
        detIndices = detIndices[i,1:maxDetNums],
        size = size[1:n.detectors],
        p0 = p0[1:n.counties],
        sigma = sigma,
        s = sxy[i,1:2],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
        localTrapsNum = localDetNum[1:n.habWindows],
        resizeFactor = resizeFactor,
        lengthYCombined = maxDetNums,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        indicator = z[i],
        trapCovsIntercept = detCounties[1:n.detectors],
        trapCovs = detCovs[1:n.detectors,1:n.covs],
        trapBetas = betaCovs[1:n.covs],
        indBeta = betaResponse,
        indCov = detResponse[i])
      
      y.Oth[i,1:maxDetNumsOth] ~ dbinomLocal_normalWolverine(
        detNums = detNumsOth[i],
        detIndices = detIndicesOth[i,1:maxDetNumsOth],
        size = size[1:n.detectors],
        p0 = p0Oth[1:n.countries],
        sigma = sigma,
        s = sxy[i,1:2],
        trapCoords = detector.xy[1:n.detectors,1:2],
        localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
        localTrapsNum = localDetNum[1:n.habWindows],
        resizeFactor = resizeFactor,
        lengthYCombined = maxDetNumsOth,
        habitatGrid = habitatGrid[1:y.max,1:x.max],
        indicator = z[i],
        trapCovsIntercept = detCountries[1:n.detectors],
        trapCovs = detCovsOth[1:n.detectors,1:n.covs.Oth],
        trapBetas = betaCovsOth[1:n.covs.Oth],
        indBeta = betaResponseOth,
        indCov = detResponse[i])
    }#i
    
    
    ##----- DERIVED PARAMETERS -----
    
    N <- sum(z[1:n.individuals])
    
  })
  
  
  
  ## ------   2. SCRize NIMBLE INPUT DATA ------
  
  for(thisSex in c("female","male")){
    
    ## [PD] : For now, I am only creating the last year's SCR dataset (faster)
    #for(t in 1:n.years){ 
    t <- n.years
    
    message(paste0("Preparing SCR input for sex: ", thisSex, "... "))
    
    for(ch in 1:4){
      
      ##-- Load OPSCR input
      load( file.path( working.dir, "nimbleInFiles", thisSex,
                       paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
      
      ##-- Identify detected individuals     
      detectedStruc <- apply(nimData$detNums,2,function(x) x > 0)
      detectedOth <- apply(nimData$detNumsOth,2,function(x) x > 0)
      detected <- detectedOth + detectedStruc
      detected <- detected > 0
      n.detected <- sum(detected[ ,t])    ## Number of detected ids that year
      n.augmented <- n.detected * 2       ## Decide the augmentation factor here
      
      
      ## ------     2.1. NIMBLE DATA ------
      
      ##-- y.alive
      nimData$y <- nimData$y[detected[ ,t], ,t]  
      nimData$y <- rbind( nimData$y,
                          matrix( 0,
                                  nrow = n.augmented,
                                  ncol = nimConstants$maxDetNums))
      
      nimData$y.Oth <- nimData$y.Oth[detected[ ,t], ,t]  
      nimData$y.Oth <- rbind( nimData$y.Oth,
                              matrix( 0,
                                      nrow = n.augmented,
                                      ncol = nimConstants$maxDetNumsOth))
      
      ##-- detNums 
      nimData$detNums <- nimData$detNums[detected[ ,t],t]
      nimData$detNums  <- c(nimData$detNums, rep(0,n.augmented))
      nimData$detNumsOth <- nimData$detNumsOth[detected[ ,t],t]
      nimData$detNumsOth  <- c(nimData$detNumsOth, rep(0,n.augmented))
      
      ##-- detIndices 
      nimData$detIndices <- nimData$detIndices[detected[ ,t], ,t]
      nimData$detIndices <- rbind( nimData$detIndices,
                                   matrix(0, 
                                          nrow = n.augmented,
                                          ncol = nimConstants$maxDetNums))
      nimData$detIndicesOth <- nimData$detIndicesOth[detected[ ,t], ,t]
      nimData$detIndicesOth <- rbind( nimData$detIndicesOth, 
                                      matrix( 0,
                                              nrow = n.augmented,
                                              ncol = nimConstants$maxDetNumsOth))
      
      ##-- z
      nimData$z <- nimData$z[detected[ ,t],t]  
      nimData$z[nimData$z %in% c(2)] <- 1 # ALIVE IDS BECOMES 1
      nimData$z <- c(nimData$z, rep(NA,n.augmented))
      
      ## [PD]: alternative
      nimData$z <- c(rep(1,n.detected), rep(NA,n.augmented))
      
      
      ##-- sxy 
      nimData$sxy <- NULL
      
      ##-- detResponse 
      nimData$detResponse <- nimData$detResponse[detected[ ,t],t]
      nimData$detResponse <- c(nimData$detResponse, rep(NA,n.augmented))
      
      ##-- detCovs
      nimData$detCovs <- nimData$detCovs[ ,t, ]
      nimData$detCovsOth <- nimData$detCovsOth[ ,t, ]
      
      ##-- density covariate
      nimData$denCounts <- nimData$denCounts
      
      ##-- detCountries
      nimData$detCountries <- nimData$detCountries[ ,t]
      
      
      
      ## ------     2.2. NIMBLE CONSTANTS ------
      
      ##-- countyToggle to toggle off Norbotten
      nimConstants$countyToggle <- nimConstants$countyToggle[ ,t]
      nimConstants$countyToggleOth <- nimConstants$countyToggleOth[ ,t]
      
      ##-- number of individuals
      nimConstants$n.individuals <- n.detected + n.augmented
      
      
      
      ## ------     2.3. NIMBLE INITS ------
      
      ##-- z
      nimInits$z <- nimInits$z[detected[,t],t]  
      nimInits$z <- c(nimInits$z, rbinom(n.augmented,1,0.5))
      ## [PD]: alternative
      nimInits$z <- c(rep(NA,n.detected), rbinom(n.augmented,1,0.5))
      
      ##-- detResponse
      nimInits$detResponse <- c( rep(NA, n.detected),
                                 rbinom(n.augmented,1,0.5))   ## HERE IT IS TREATED AS A LATENT COVARIATE
      
      ##-- sxy 
      nimInits$sxy <- nimInits$sxy[detected[ ,t], ,t]  
      nimInits$sxy <- rbind( nimInits$sxy,              ## GIVE ACS FROM DETECTED INDIVIDUALS TO AUGMENTED IDS. 
                             nimInits$sxy[sample(nimInits$sxy, n.augmented, replace = T), ])
      
      ##-- p0
      nimInits$p01 <- nimInits$p01[ ,t]
      nimInits$p01Oth <- nimInits$p01Oth[ ,t]
      
      ##-- psi
      nimInits$psi <- runif(1,0.4,0.6)
      
      ##-- sigma
      nimInits$sigma <- runif(1,1,2)
      
      ##-- betaCovs
      nimInits$betaCovs <- nimInits$betaCovs[ ,t]
      nimInits$betaCovsOth <- nimInits$betaCovsOth[ ,t]
      
      ##-- betaResponse
      nimInits$betaResponse <- nimInits$betaResponse[t]
      nimInits$betaResponseOth <- nimInits$betaResponseOth[t]
      
      
      
      ## ------     2.4. NIMBLE INITS ------
      
      nimParams <- c("N", "betaDens", "psi", 
                     "pResponse", "betaResponse", "sigma", 
                     "p0", "betaCovs", "betaResponse",
                     "p0Oth", "betaCovsOth","betaResponseOth")
      
      nimParams2 <- c("z", "sxy")
      
      
      
      ## ------     2.5. SAVE INPUTS ------
      
      modelCode <- modelCode1
      
      save( modelCode,
            nimData,
            nimConstants,
            nimInits,
            nimParams,
            nimParams2,
            file = file.path( working.dir, "nimbleInFiles", thisSex,
                              paste0("SCRinput", years[t], "_", DATE, "_", thisSex, "_", c, ".RData")))
    }#c
  }#thisSex
  
  
  
  ##----------------------------------------------------------------------------
  
  ## ------ IV. RETURN IMPORTANT INFOS FOR REPORT ------
  
  return(list( SPECIES = "Wolverine",
               engSpecies = "wolverine",
               YEARS = years,
               SEX = sex,
               DATE = DATE))
  
  
  
  ##----------------------------------------------------------------------------
}

