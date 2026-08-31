#' @title RovQuant OPSCR wolf data preparation.
#'
#' @description
#' \code{makeRovquantData_wolf} formats the available wolf data for the OPSCR analysis using nimble and nimbleSCR.
#' The data preparation process is composed of three main steps:
#'  - defining and formatting habitat characteristics
#'  - defining and formatting detectors characteristics
#'  - defining and formatting individual detection histories
#'
#' @name makeRovquantData_wolf
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
#' @param x.extent (Optional) A \code{Numeric}.
#' @param y.extent (Optional) A \code{Numeric}.
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
#' @importFrom readxl read_excel
#' 
NULL
#' @rdname makeRovquantData_wolf
#' @export
makeRovquantData_wolf <- function(
  ##-- paths
  data.dir = getwd(),
  working.dir = getwd(),
  
  ##-- data
  years = NULL,
  sex = c("female","male"),
  aug.factor = 0.8,
  sampling.months = list(10:12,1:3),
  
  ##-- habitat
  habitat.res = 20000,
  x.extent = NULL,
  y.extent = NULL,
  buffer.size = 40000,
  max.move.dist = 250000,
  
  ##-- detectors
  detector.res = 10000,
  subdetector.res = 1000,
  max.det.dist = 45000,
  resize.factor = 1,
  
  ##-- Miscellanious
  rename.list = NULL
){
  
  ## ------ 0. BASIC SET-UP ------
  
  ##-- Set default values for the wolf model
  if(is.null(aug.factor)){aug.factor <- 0.8}
  if(is.null(sampling.months)){sampling.months <- list(10:12,1:3)}
  if(is.null(habitat.res)){habitat.res <- 20000} 
  if(is.null(x.extent)){x.extent <- c(210000,760000)}
  if(is.null(y.extent)){y.extent <- c(6000000,7050000)}
  if(is.null(buffer.size)){buffer.size <- 40000}
  if(is.null(max.move.dist)){max.move.dist <- 250000}
  if(is.null(detector.res)){detector.res <- 10000}
  if(is.null(subdetector.res)){subdetector.res <- 1000}
  if(is.null(max.det.dist)){max.det.dist <- 45000}
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
  
  
  ## ---------------------------------------------------------------------------
  
  ## ------ I. LOAD AND SELECT DATA ------
  
  ## ------   1. HABITAT DATA -----
  
  ##-- Load pre-defined habitat rasters and shapefiles
  data(habitatRasters, envir = environment()) 
  data(REGIONS, envir = environment())
  data(studyAreaWolf, envir = environment())
  
  
  ##-- Disaggregate habitat raster to the desired resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)
  
  ##-- Merge counties for practical reasons
  COUNTIES_AGGREGATED <- REGIONS %>%
    mutate(id = case_when(
      county %in% c("Akershus","Agder","Buskerud","Vestfold","Oslo","Østfold","Telemark") ~ "NO1",
      county %in% c("Innlandet","Møre og Romsdal","Trøndelag") ~ "NO2",
      county %in% c("Jämtland","Västernorrland","Västerbotten","Dalarna","Gävleborg") ~ "SE1",
      county %in% c("Uppsala","Västmanland","Stockholm","Södermanland") ~ "SE2",
      county %in% c("Blekinge","Örebro","Östergötland","Jönköping","Kronoberg","Kalmar","Skåne","Gotlands") ~ "SE3",
      county %in% c("Västra Götaland","Värmland","Halland") ~ "SE4")) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise()
  

  ## ------   2. NGS DATA -----
  
  ##-- Extract date from the last cleaned data file
  DATE <- getMostRecent( 
    path = file.path(working.dir, "data"),
    pattern = "CleanData_wolf")
  
  ##-- Load the most recent clean wolf data from RovBase
  myFullData.sp <- readMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_wolf",
    extension = ".RData")
  
  ##-- List years
  if(is.null(years)){
    years <- sort(unique(c(myFullData.sp$alive$Year,
                           myFullData.sp$dead.recovery$Year)))
  }
  DATA$years <- years
  n.years <- length(years)
  
  
  
  ## ---------------------------------------------------------------------------
  
  ## ------ II. CREATE OPSCR DATA ------
  
  ## ------   1. GENERATE HABITAT ------
  
  message("Preparing habitat characteristics... ")
  
  ## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------
  ##-- Determine study area based on predefined extent
  #[CM] use the same study area to align with previous analyses
  #[CM] new dataset added to the data package
  studyArea <- myStudyArea.poly
  studyAreaExtent <- myStudyArea.extent
  
  #[CM] Commented out
  # studyArea <- sf::st_crop( REGIONS,
  #                           xmin = x.extent[1], xmax = x.extent[2],
  #                           ymin = y.extent[1], ymax = y.extent[2]) %>%
  #   sf::st_collection_extract(., "POLYGON") %>%
  #   summarise()  

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
    sf::st_set_crs( ., value = sf::st_crs(habitat$buffered.habitat.poly))
  
  ##-- Study area grid from habitat raster
  habitat.rWthBufferPol <- sf::st_as_sf( 
    stars::st_as_stars(habitat$habitat.rWthBuffer), 
    as_points = FALSE,
    merge = TRUE) %>%
    dplyr::filter(Habitat %in% 1)
  
  
  
  ## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------
  
  ## ------       1.2.1. DENSITY OF PACKS/PAIRS ------
  
  ##-- Kernel of NGS detections of individuals in pairs
  #[CM] move this to be sex-specific as in previous analyses 
  # kern <- list()
  # habDens <- matrix(NA, nrow = n.habWindows, ncol = n.years)
  # for(t in 1:n.years){
  #   ##-- Subset the NGS data to individuals in packs/pairs this year
  #   data.pairs.t <- myFullData.sp$alive %>%
  #     dplyr::filter( Year == years[t],
  #                    STATUS %in% c(3,4),
  #                    Sex=="male")
  #   
  #   ##-- Get mean coordinates of packs
  #   IDs <- unique(data.pairs.t$IdSimplified)
  #   m.xy <- matrix(NA, nrow = length(IDs), ncol = 2)
  #   colnames(m.xy) <- c("x","y")
  #   for(i in 1:length(IDs)){
  #     m.xy[i, ] <- data.pairs.t %>%
  #       dplyr::filter( IdSimplified == IDs[i]) %>%
  #       st_coordinates(.) %>%
  #       colMeans(.)
  #   }#i
  #   
  #   ##-- Check if some coordinates are missing  
  #   if(sum(is.na(m.xy[ ,1])) > 0){m.xy <- m.xy[!is.na(m.xy[ ,1]), ]}
  #   
  #   ##-- Turn into .sf
  #   locationsFamily <- st_as_sf( as.data.frame(m.xy),
  #                                coords = c("x","y"),
  #                                crs = st_crs(habitat$habitat.sp))
  #   locationsFamily$id <- rep(1, nrow(locationsFamily))
  #   
  #   ##-- Calculate kernel of detections
  #   kern[[t]] <- raster(estUDm2spixdf(kernelUD( 
  #     as(locationsFamily[ ,"id"], "Spatial"), 
  #     h = 15000,
  #     grid = as(habitat$habitat.r, 'SpatialPixels'))))
  #   
  #   ##-- Plot check
  #   plot(kern[[t]], main = years[t])
  #   plot(habitat$habitat.poly$geometry, add = T, col = NA)
  #   
  #   ##-- Scale covariate
  #   habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ] == 1])
  # } #t
  
  
  
  ## ------   2. GENERATE DETECTORS -----
  
  message("Preparing detectors characteristics... ")
  
  ## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----
  
  ##-- Generate raster of sub-detectors based on the study area 
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
  
  
  
  ## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----
  
  ## ------       2.2.1. EXTRACT COUNTRIES ------
  
  ##-- Extract closest country for each detector
  detCountries <- detectors$main.detector.sp %>%
    sf::st_distance(., COUNTRIES, by_element = F) %>%
    apply(., 1, function(x) which.min(x)) %>%
    as.factor(.) %>%
    as.numeric(.)
  
  ##-- Put into "nimble2SCR" format
  detectors$detectors.df$countries <- detCountries
  
  
  
  ## ------       2.2.2. EXTRACT COUNTIES ------
  
  ##-- Extract closest county for each detector
  detCounties <- detectors$main.detector.sp %>%
    sf::st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
    apply(., 1, function(x) which.min(x))
  
  ##-- Put into "nimble2SCR" shape
  detectors$detectors.df$counties <- detCounties
  
  
  
  ## ------       2.2.3. EXTRACT GPS TRACKS LENGTHS ------
  
  message("Cleaning GPS tracks... ")
  TRACKS_MULTI <- sf::read_sf(file.path(data.dir, "Tracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250422.shp")) %>%
    ##-- Process dates
    dplyr::mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
                   Mth = as.numeric(format(Dato,"%m")),
                   Yr = as.numeric(format(Dato,"%Y")),
                   Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
    ##-- Filter out irrelevant tracks
    dplyr::filter( Helikopter == "0",      ## Remove helicopter tracks
                   Year %in% years & Mth %in% unlist(sampling.months))

  ##-- LOAD GPS SEARCH TRACKS FROM BOUVET (EXPORT 2025)
  TRACKS_SINGLE_25 <- sf::read_sf(file.path(data.dir, "Tracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250422.shp")) %>%
    ##-- Process dates
    dplyr::mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
                   Mth = as.numeric(format(Dato,"%m")),
                   Yr = as.numeric(format(Dato,"%Y")),
                   Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
    ##-- Filter out irrelevant tracks
    dplyr::filter( Yr < 2020,
                   Helikopter == "0",      ## Remove helicopter tracks
                   Year %in% years & Mth %in% unlist(sampling.months))

  ##-- LOAD GPS SEARCH TRACKS FROM BOUVET (EXPORT 2026)
  TRACKS_SINGLE_26 <- sf::read_sf(file.path(data.dir, "Tracks/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20260420.shp")) %>%
    ##-- Process dates
    dplyr::mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
                   Mth = as.numeric(format(Dato,"%m")),
                   Yr = as.numeric(format(Dato,"%Y")),
                   Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
    ##-- Filter out irrelevant tracks
    dplyr::filter( Helikopter == "0",      ## Remove helicopter tracks
                   Year %in% years & Mth %in% unlist(sampling.months))

  ## COMBINE ALL TRACKS
  ALL_TRACKS <- rbind(TRACKS_SINGLE_25, TRACKS_SINGLE_26, TRACKS_MULTI)
  rm(list = c("TRACKS_SINGLE_25", "TRACKS_SINGLE_26", "TRACKS_MULTI"))

  ## Check
  ## SELECT TRACKS YEAR
  dupIDs <- dupDist <- length <- TRACKS_YEAR <- TRACKS_YEAR.sp <- list()
  for(t in 1:length(years)){
    ## SUBSET GPS TRACKS TO THE SAMPLING PERIOD
    TRACKS_1 <- ALL_TRACKS[ALL_TRACKS$Yr%in% years[t] &
                             ALL_TRACKS$Mth%in%sampling.months[[1]], ]
    TRACKS_2 <- ALL_TRACKS[ALL_TRACKS$Yr%in% (years[t]+1) &
                             ALL_TRACKS$Mth%in%sampling.months[[2]], ]
    TRACKS <- rbind(TRACKS_1, TRACKS_2)
    ## SIMPLIFY TRACKS SHAPES
    ## SUBSET TRACKS TO THE STUDY AREA
    TRACKS <- st_intersection(TRACKS, st_as_sfc(studyAreaExtent))
    ## NAME TRACKS
    TRACKS$ID <- 1:nrow(TRACKS)
    TRACKS_YEAR[[t]] <- TRACKS
    # calculate length to identify duplicates
    TRACKS_YEAR[[t]]$dist <- st_length(TRACKS_YEAR[[t]], byid = T)
    # calculate centroids to also identify tracks that could have the same length but in different location
    TRACKS_YEAR[[t]]$centroidx <-  st_coordinates(st_centroid(TRACKS_YEAR[[t]]))[,1]
    df <- data.frame(ID = TRACKS_YEAR[[t]]$ID,
                     Dato = TRACKS_YEAR[[t]]$Dato,
                     Person = TRACKS_YEAR[[t]]$Person,
                     dist = TRACKS_YEAR[[t]]$dist,
                     centroidx = TRACKS_YEAR[[t]]$centroidx)
    dupIDs[[t]] <- duplicated(df[ ,2:5])# find duplicates based on person, distance and date
    dupIDs[[t]] <- df$ID[dupIDs[[t]]]
    dupDist[[t]] <- TRACKS_YEAR[[t]][dupIDs[[t]], ]$dist
    TRACKS_YEAR[[t]] <-  TRACKS_YEAR[[t]][-dupIDs[[t]], ]
  }#t
  
  ##-- Combine all GPS tracks
  # [CM] comment out this function. Use the same script than in previous analysis instead
  # TRACKS <- readTracks( data.dir = data.dir,
  #                       years = years,
  #                       sampling.months = sampling.months)
  
  ##-- Extract length of GPS search track per detector grid cell
  detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
  ##-- Set-up progress bar
  pb = utils::txtProgressBar( min = 1, max = n.years, initial = 0, style = 3)
  for(t in 1:n.years){
    intersection <- TRACKS_YEAR[[t]] %>%
      #dplyr::filter(Year == years[t]) %>%
      sf::st_intersection(detectors$grid, .) %>%
      dplyr::mutate(LEN = st_length(.)) %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(transect_L = sum(LEN)) 
    detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
    ##-- Print progress 
    utils::setTxtProgressBar(pb,t) 
  }#t
  
  ##-- Put into "nimble2SCR" format
  colnames(detTracks) <- paste0("tracks.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detTracks)
  
  
  
  ## ------       2.2.4. EXTRACT DISTANCES TO ROADS ------
  
  ##-- Load map of distance to roads (1km resolution)
  DistAllRoads <- raster(file.path(data.dir, "Roads/MinDistAllRoads1km.tif"))
  
  r <- fasterize(studyArea, DistAllRoads)
  r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  DistAllRoads <- r
  DistAllRoads <- crop(DistAllRoads, studyArea)
  
  # [CM] comment out this function. Use the same script than in previous analysis instead
  # DistAllRoads <- readMostRecent( path = file.path(data.dir, "Roads"), 
  #                                 extension = ".tif", 
  #                                 stack = FALSE)
  
  ##-- Fasterize to remove values that fall in the sea
  ##-- Aggregate to match the detectors resolution
  DistAllRoads <- raster::aggregate( 
    x = DistAllRoads,
    fact = detectors$resolution/raster::res(DistAllRoads),
    fun = mean)
  
  ##-- Extract distance to road for each detector
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
  
  SNOW <- stack(file.path(data.dir, "Snow/AverageSnowCoverModisSeason2016_2026_Wolf.tif"))

  ##-- RENAME THE LAYERS
  names(SNOW) <- paste(years,(years)+1, sep = "_")
  
  ##-- Load raster stack of snow cover
  # [CM] comment out this function. Use the same script than in previous analysis instead
  # SNOW <- readMostRecent( path = file.path(data.dir, "Snow"), 
  #                         extension = ".tif", 
  #                         stack = TRUE)
  
  ##-- Select snow data corresponding to the monitoring period
  SNOW <- SNOW[[paste("X", years, "_", years + 1, sep = "")]]
  SNOW <- raster::crop(SNOW, c(0,40,55,75))
  
  ##-- Extract snow 
  detSnow <- matrix(0, nrow = dim(detectors$main.detector.sp)[1], ncol = n.years)
  det.sptransf <- st_transform(detectors$main.detector.sp, st_crs(SNOW))
  detSnow[ ,1:n.years] <- raster::extract(SNOW, det.sptransf)
  
  ##-- if NA returns the average value of the cells within 15km 
  isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
  tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                          buffer = 15000, fun = mean, na.rm = T)
  detSnow[isna,1:n.years] <- tmp
  
  ##-- still some NA... Increase buffer again 
  isna <- which(apply(detSnow, 1, function(x)any(is.na(x))))
  tmp <- raster::extract( SNOW, det.sptransf[isna, ],
                          buffer = 35000, fun = mean, na.rm = T)
  detSnow[isna,1:n.years] <- tmp
  
  ##-- Put into "nimble2SCR" format
  colnames(detSnow) <- paste0("snow.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detSnow)
  
  
  
  ## ------       2.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------
  
  ## ------         2.2.6.1. SKANDOBS ------
  
  ##-- Load the last SkandObs data file
  skandObs <- readMostRecent( path = file.path(data.dir, "Skandobs"),
                              extension = ".xlsx",
                              pattern = "Skandobs")
  
  ##-- Replace scandinavian characters
  colnames(skandObs) <- rovquantR::translateForeignCharacters(dat = colnames(skandObs))#,[CM]
                                                   #dir.translation = dir.analysis)
  
  skandObs <- skandObs %>%
    ##-- Extract important info (e.g. month, year)
    dplyr::mutate( date = as.POSIXct(strptime(date, "%Y-%m-%d")),
                   year = as.numeric(format(date,"%Y")),
                   month = as.numeric(format(date,"%m")),
                   species = stringi::stri_trans_general(species, "Latin-ASCII"),
                   monitoring.season = ifelse(month < 12,#[CM]UPDATE to 10 #unlist(sampling.months)[1],
                                              year-1, year)) %>%
    ##-- Filter based on monitoring season
    dplyr::filter( month %in% unlist(sampling.months)) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf(., coords = c("longitude","latitude")) %>%
    sf::st_set_crs(., value = "EPSG:4326") %>%
    sf::st_transform(., sf::st_crs(REGIONS)) %>%
    sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------         2.2.6.2. ROVBASE ------
  
  rovbaseObs1 <- read_excel(file.path(data.dir, "AllSamples/RIB1804202607011165.xlsx"))
  rovbaseObs2 <- read_excel(file.path(data.dir, "AllSamples/RIB1804202607015508.xlsx"))
  rovbaseObs3 <- read_excel(file.path(data.dir, "AllSamples/RIB18042026065900566.xlsx"))
  rovbaseObs4 <- read_excel(file.path(data.dir, "AllSamples/RIB18042026070006976.xlsx"))
  
  rovbaseObs <- rbind(rovbaseObs1,rovbaseObs2,rovbaseObs3,rovbaseObs4)
  rm(list = c("rovbaseObs1", "rovbaseObs2", "rovbaseObs3", "rovbaseObs4"))
  
  colnames(rovbaseObs) <- translateForeignCharacters( dat = colnames(rovbaseObs))
  rovbaseObs$Sample_type <- translateForeignCharacters( dat = rovbaseObs$Proevetype)
  
  ##-- Process Rovbase observations (all species)
  rovbaseObs <- rovbaseObs %>%
                # [CM] comment out readMultiples. Use the same script than in previous analysis instead
                # readMultiples( path = file.path(data.dir, "AllSamples"),
                               #extension = ".xlsx") %>%
    ##-- Rename columns to facilitate manipulation
    dplyr::rename(., any_of(rename.list)) %>%
    ##-- Extract important info (e.g. month, year, country of collection)
    dplyr::mutate(
      ##-- Turn potential factors into characters 
      across(where(is.factor), as.character),
      ##-- Deal with Scandinavian characters
      Species = stringi::stri_trans_general(Species, "Latin-ASCII"),
      # [CM]
      # Sample_type = translateForeignCharacters(dat=Sample_type,dir.analysis),
      ##-- Deal with dates
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      year = as.numeric(format(Date,"%Y")),
      month = as.numeric(format(Date,"%m")),
      monitoring.season = ifelse(month < 12,#[CM] UPDATE to 10 #unlist(sampling.months)[1],
                                 year-1, year)) %>%
    ##-- Filter out unusable samples
    dplyr::filter( 
      ##-- Filter out samples without coordinates,...
      !is.na(North_UTM33),
      ##-- ...based on species
      Species %in% c("Ulv"),
      ##-- ...based on sample type
      Sample_type %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
                          "Saliv/Spytt"),#[CM] comment out , "Loepeblod", "Vev"),
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months)
      # ##-- ... if sample was from the focal species and successfully genotyped 
      # !(Species %in% "Ulv" & !is.na(Id))
    ) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf( ., coords = c("Oest (UTM33/SWEREF99 TM)","North_UTM33")) %>%
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
  r.skandObsContinuous <- brick(lapply(r.list,function(x) x[[1]]))
  r.skandObsBinary <- brick(lapply(r.list,function(x) x[[2]]))
  r.rovbaseContinuous <- brick(lapply(r.list,function(x) x[[3]]))
  r.rovbaseBinary <- brick(lapply(r.list,function(x) x[[4]]))
  
  ##-- Combine both rasters
  r.SkandObsRovbaseBinary <- r.rovbaseBinary + r.skandObsBinary
  for(t in 1:n.years){
    r.SkandObsRovbaseBinary[[t]][r.SkandObsRovbaseBinary[[t]][] > 1] <- 1
  }#t
  
  
  
  ## ------         2.2.6.4. ASSIGN THE COVARIATE ------
  
  detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
  detOtherSamples[ ,1:n.years] <- raster::extract( r.SkandObsRovbaseBinary,
                                                   detectors$main.detector.sp)
  
  ##-- Put into "nimble2SCR" format
  colnames(detOtherSamples) <- paste0("detOtherSamples.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detOtherSamples)
  
  
  
  ## ------       2.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------
  
  detSnow <- round(scale(detSnow), digits = 2)
  detRoads <- round(scale(detRoads), digits = 2)
  detTracks <- round(scale(detTracks), digits = 2)

  ##-- STRUCTURED 
  detCovs <- array(NA, c(dim(detTracks)[1], 2, dim(detTracks)[2]))
  detCovs[ ,1, ] <- detTracks
  detCovs[ ,2, ] <- detSnow
  dimnames(detCovs) <- list( "detectors" = 1:n.detectors,
                             "covariates" = c("tracks", "snow"), 
                             "years" = years)
  
  ##-- OTHERS 
  detCovsOth <- array(NA, c(dim(detTracks)[1], 3, dim(detTracks)[2]))
  detCovsOth[ ,1, ] <- detRoads
  detCovsOth[ ,2, ] <- detSnow
  detCovsOth[ ,3, ] <- detOtherSamples
  dimnames(detCovsOth) <- list( "detectors" = 1:n.detectors,
                                "covariates" = c( "roads","snow","obs"),
                                "years" = years)
  
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
  detectors$localObjects <- getLocalObjects(#[CM] why not using nimbleSCR version?
    habitatMask = habitat$habitat.mx,
    coords = detectors$scaledCoords,
    dmax = detectors$maxDist*1.4/habitat$resolution,
    resizeFactor = detectors$resize.factor,
    plot.check = F)

  
  
  ## ------   5. SAVE STATE-SPACE CHARACTERISTICS -----
  
  save( habitat,
        file = file.path( working.dir, "data",
                          paste0("Habitat_wolf_", DATE, ".RData")))
  
  save( detectors,
        file = file.path( working.dir,"data",
                          paste0("Detectors_wolf_", DATE, ".RData")))
  
  
  
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
    ##-- Filter based on space #[CM] use extent instead of the grid based as in previous analyses
    sf::st_filter( .,sf::st_as_sfc(myStudyArea.extent), .predicate = st_intersects)
    # sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------     6.2. DEAD RECOVERY DATA -----
  
  data.dead <- myFullData.sp$dead.recovery %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to sex of interest
      Sex %in% sex) %>%
    ##-- Filter based on space #[CM] use extent instead of the grid based as in previous analyses
    sf::st_filter( .,sf::st_as_sfc(myStudyArea.extent), .predicate = st_intersects)
    # sf::st_filter( .,habitat.rWthBufferPol, .predicate = st_intersects)
  
  
  
  ## ------     6.3. SEPARATE STRUCTURED & OPPORTUNISTIC SAMPLING ------
  
  ## ------       6.3.1. ASSIGN SAMPLES TO GPS SEARCH TRACKS ------
  
  message("Assigning DNA samples to GPS tracks... ")
  message("This can take several minutes... ")
  
  #[CM] Couldnt get this to reproduce the results commented out 
  # data.alive <- assignSearchTracks(
  #   data = data.alive,
  #   tracks = TRACKS)
  # rm(list = c("TRACKS"))
  
  TRACKSSimple_sf <- list()
  for(t in 1:length(years)){
    TRACKS_YEAR[[t]]$RovbsID <- as.character(TRACKS_YEAR[[t]]$RovbaseID)
    TRACKS_YEAR[[t]]$RovbasID <- 1:length(TRACKS_YEAR[[t]]$RovbaseID)
    TRACKSSimple_sf[[t]] <- TRACKS_YEAR[[t]]
  }#t
  
  data.alive$TrackRovbsID <- NA
  data.alive$trackDist <- NA
  
  ## ASSIGN EACH SAMPLE TO THE CLOSEST TRACK
  dnatemp <- st_as_sf(data.alive)
  
  ## CREATE A BUFFER AROUND EACH DETECTION
  tmp <- st_buffer(dnatemp, dist = 750)
   
  for(i in 1:nrow(data.alive)){
    # INTERSECT POINT WITH TRACKS,
    t <- which(years %in% tmp[i, ]$Year)
    whichSameDate <- which(as.character(TRACKSSimple_sf[[t]]$Dato)==as.character(data.alive$Date[i]))
    tmpTRACKS <- st_intersection(TRACKSSimple_sf[[t]][whichSameDate,], tmp[i,])

    if(nrow(tmpTRACKS)==0){next}

    # FIND THE CLOSEST TRACK
    dist <- st_distance(dnatemp[i,], tmpTRACKS, by_element = F)

    # MAKE SURE THE SAMPLE WAS COLLECTED AT THE SAME TIME THAN THE TRACK
    # IF NO MATCHING DATE ASSIGN TO NA?
    if(length(dist)==0){
      data.alive$TrackRovbsID[i] <- NA
      data.alive$trackDist[i] <- NA
    }
    # IF MATCHING DATE ASSING TO THAT TRACK
    if(length(dist)==1){
      data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID
      data.alive$trackDist[i] <- dist
    }
    # IF SEVERAL MATCHING DATES ASSING TO THE CLOSEST OF THE MATCHING TRACKS
    if(length(dist)>1){
      data.alive$TrackRovbsID[i] <- tmpTRACKS$RovbsID[which.min(dist)]
      data.alive$trackDist[i] <- min(dist)
    }
    #print(i)
  }
  
  
  
  ## ------       6.3.2. ASSIGN SAMPLES TO OPPORTUNISTIC OR STRUCTURED ------
  
  distanceThreshold <- 500
  
  ##-- Identify samples from structured and opportunistic sampling
  data.alive <- data.alive %>%
    dplyr::mutate(
      ##-- Collector column was replaced by two columns, merging them now...
      Collector_role1 = ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
      ##-- Identify samples collected during structured sampling 
      structured = Collector_role1 %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
        !is.na(TrackRovbsID) &
        trackDist <= distanceThreshold)
  
  

  ## ------     6.4. ASSIGN SAMPLES TO DETECTORS -----
  
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
  
  
  
  ## ------     6.5. PLOT NGS & DEAD RECOVERY MAPS ----- 
  
  ##-- layout
  L <- n.years
  if(L < 6){ nrows <- 1 } else{
    if(L < 13){ nrows <- 2 } else {
      if(L < 22){ nrows <- 3 } else {
        if(L < 33){ nrows <- 4 } else {
          nrows <- 5
        }}}}
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
      plot( sf::st_geometry(data.alive$data.sp[data.alive$data.sp$Year == years[t], ]),
            add = TRUE, col = "orange", pch = 3),
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
    mx[r*2, ] <- c(rep(1:ncols, each = 2),ncols) + (r-1)*ncols
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
  
  
  
  ## ------     6.6. SAVE FILTERED DATA ----- 
  
  save( data.alive, data.dead,
        file = file.path( working.dir, "data",
                          paste0("FilteredData_wolf_", DATE, ".RData")))
  
  
  
  ## ------   7. GENERATE DETECTION HISTORY ------
  
  for(thisSex in sex){
    
    message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
    
    ## ------     7.1. FILTER DATA BY SEX -----
    
    load(file.path( working.dir, "data",
                    paste0("FilteredData_wolf_", DATE, ".RData")))
    
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
      distances[[t]] <- checkDistanceDetections( 
        y = y.ar$y.ar[ , ,t], 
        detector.xy = detectors$detectors.df[ ,c("x","y")], 
        max.distance = detectors$maxDist,
        method = "pairwise",
        plot.check = F)
      
      ##-- If some detections are flagged
      if(sum(distances[[t]]$y.flagged) > 0){
        ##-- Remove detections that are further then the threshold
        y.ar.ALIVEOthers[ , ,t] <- y.ar.ALIVEOthers[ , ,t] * (1-distances[[t]]$y.flagged)
        y.ar.ALIVEStructured[ , ,t] <- y.ar.ALIVEStructured[ , ,t] * (1-distances[[t]]$y.flagged)
        
        ##-- Remove detections also in data.alive$data.sp to run getSInits later
        affected.ids <- which(apply(distances[[t]]$y.flagged,1,sum)>0)
        idd <- names(affected.ids)
        for(i in 1:length(idd)){
          detIds <- which(distances[[t]]$y.flagged[idd[i], ] > 0)
          tmp <- data.alive$data.sp %>%#
            dplyr::filter(!(Id %in% idd[i] &
                              Detector %in% detIds &
                              Year %in% years[t]))
          data.alive$data.sp <- tmp#[CM] for whatever reasons i need this extra step
        }#i
      }#if
    }#t
    
    
    
    ## ------     7.4. GENERATE INDIVIDUAL-LEVEL COVARIATES ------ 
    
    ## ------       7.4.1. INDIVIDUAL STATE ------ 
    ##[CM] Do the state assignement here instead of within cleanData function
    ##[CM] State assignment based on NGS alive is challenging as we also need to do it for dead wolves. 
    ##[CM] It is much easier to do it once we create the entire y dataframe
    ##-- Load most recent Micke's file
    INDIVIDUAL_ID <- suppressWarnings(readMostRecent( path = data.dir,
                                                      extension = ".xls",
                                                      pattern = "Grouping")) %>%
      ##-- Rename columns to facilitate manipulation
      dplyr::rename(., any_of(rename.list),
                    IdSimplified = "ROVBASE_IndividID") %>%
      ##-- Turn potential factors into characters 
      dplyr::mutate(across(where(is.factor), as.character)) %>%
      ##-- Add some columns
      dplyr::mutate( 
        Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
        Sex = ifelse(is.na(Sex), "unknown", Sex),
        Sex = ifelse(Sex %in% "Hona", "female", Sex),
        Sex = ifelse(Sex %in% "Hane", "male", Sex),
        Year = 2021)  
    
    
    ##-- THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2022/23.
    Pack_ID2023 <- suppressWarnings(readMostRecent( path = data.dir,
                                                    extension = ".xls",
                                                    pattern = "Genetiskt ID")) %>%
      ##-- Rename columns to facilitate manipulation
      dplyr::rename(.,
                    any_of(rename.list),
                    IdSimplified = "Rovbase-ID") %>%
      ##-- Turn potential factors into characters
      dplyr::mutate(across(where(is.factor), as.character)) %>%
      ##-- Add some columns
      dplyr::mutate(
        ##-- Add status 
        Status = "Pair",
        ##-- Fix unknown "Sex"
        Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
        Sex = ifelse(is.na(Sex), "unknown", Sex),
        Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
        Sex = ifelse(Sex %in% c("Hann","Hane"), "male", Sex),
        Year = 2022)
    
    
    ##-- THIS IS THE PACK ID SENT BY LINN FOR THE WINTER 2023/24.
    Pack_ID2024 <- suppressWarnings(readMostRecent( path = data.dir,
                                                    extension = ".xls",
                                                    pattern = "Bilaga_")) %>%
      ##-- Rename columns to facilitate manipulation
      dplyr::rename(.,
                    IdSimplified = "RovbaseID",
                    any_of(rename.list)) %>%
      ##-- Turn potential factors into characters
      dplyr::mutate(across(where(is.factor), as.character)) %>%
      ##-- Add some columns
      dplyr::mutate(
        ##-- Add status 
        Status = "Pair",
        ##-- Fix unknown "Sex"
        Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
        Sex = ifelse(is.na(Sex), "unknown", Sex),
        Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
        Sex = ifelse(Sex %in% c("Hann","Hane"), "male", Sex),
        Year = 2023)
    
    
    ##-- THIS IS THE PACK ID SENT BY ØYSTEIN FOR THE WINTER 2024/25.
    Pack_ID2025 <- suppressWarnings(readMostRecent( path = data.dir,
                                                    extension = ".xls",
                                                    pattern = "FromOystein")) %>%
      ##-- Rename columns to facilitate manipulation
      dplyr::rename(.,
                    any_of(rename.list),
                    IdSimplified = "IndividID") %>%
      ##-- Turn potential factors into characters
      dplyr::mutate( across(where(is.factor), as.character),
                     Year = 2024,
                     Status = "Pair") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        ##-- Fix unknown "Sex"
        Sex = ifelse(any(c_across(Sex1:Sex4) %in% c("Tispe","Tik")),
                     "female",
                     ifelse(any(c_across(Sex1:Sex4) %in% c("Hann","Hane")),
                            "male",
                            "unknown")))
    
    
    ##-- THIS IS THE PACK ID SENT BY ØYSTEIN FOR THE WINTER 2024/25.
    Pack_ID2026 <- suppressWarnings(readMostRecent( path = data.dir,
                                                    extension = ".csv",
                                                    pattern = "RovbaseID ØF")) %>%
      ##-- Rename columns to facilitate manipulation
      dplyr::rename(.,
                    IdSimplified = "RovbaseID",
                    any_of(rename.list)) %>%
      ##-- Turn potential factors into characters
      dplyr::mutate(across(where(is.factor), as.character)) %>%
      ##-- Add some columns
      dplyr::mutate(
        ##-- Add status 
        Status = "Pair",
        ##-- Fix unknown "Sex"
        Sex = ifelse(Sex %in% "Okänt", "unknown", Sex),
        Sex = ifelse(is.na(Sex), "unknown", Sex),
        Sex = ifelse(Sex %in% c("Tispe","Tik"), "female", Sex),
        Sex = ifelse(Sex %in% c("Hann","Hane"), "male", Sex),
        Year = 2025)
    
    
    ##-- Consolidate all info on individual sex in one dataframe
    # [CM] Commented out
    # ALL_SEX <- rbind( DATA[ ,c("IdSimplified","Sex")],
    #                   INDIVIDUAL_ID[ ,c("IdSimplified","Sex")],
    #                   Pack_ID2023[ ,c("IdSimplified","Sex")],
    #                   Pack_ID2024[ ,c("IdSimplified","Sex")],
    #                   Pack_ID2025[ ,c("IdSimplified","Sex")],
    #                   Pack_ID2026[ ,c("IdSimplified","Sex")])
    # 
    # 
    ##-- Consolidate all info on individual status in one dataframe
    # ALL_STATUS <- rbind( INDIVIDUAL_ID[ ,c("IdSimplified","Year","Status")],
    #                      Pack_ID2023[ ,c("IdSimplified","Year","Status")],
    #                      Pack_ID2024[ ,c("IdSimplified","Year","Status")],
    #                      Pack_ID2025[ ,c("IdSimplified","Year","Status")],
    #                      Pack_ID2026[ ,c("IdSimplified","Year","Status")]) 
    ##-- Merge with detection data 
    # DATA <- DATA %>%
    #   left_join(., ALL_STATUS, by = c("IdSimplified","Year"))
    #overwrite Sex with Micke info
    # micke.sex <- unlist(lapply(DATA$Id,
    #                            function(i){ 
    #                              INDIVIDUAL_ID[INDIVIDUAL_ID$Id %in% i, "Sex"][1,]
    #                            }))
    # DATA$Sex <- ifelse(!is.na(micke.sex), micke.sex, DATA$Sex)
    #
    ##statut 
    # DATA$STATUS <- NA 
    # 
    # INDIVIDUAL_ID$STATUS_Numeric <- 1
    # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% "Juvenile"] <- 2
    # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Pair" )] <- 3
    # INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Family group")] <- 4
    # 
    # # indIDRovBase <- unlist(lapply(strsplit(DATA$IdSimplified," "), function(x) x[1]))
    # ALLIDS <- c(unique(DATA$IdSimplified))
    # # ALLIDS <- c(unique(DATA$Id))
    # 
    # # y.obsALL <- matrix(1, nrow = length(ALLIDS), ncol = dim(y.ar.ALIVE)[3]+1)
    # # yrs <- c(years[1]-1, years)
    # # dimnames(y.obsALL) <- list(ALLIDS, yrs)
    # nyears <- years
    # yrs <- c(years[1]-1, years)
    # t=1
    # for(i in 1:length(ALLIDS)[1]){
    #   for(t in 1:(length(years+1))){
    #     tmp <- unique(INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$`ReprodYear (May 1 year y - Apr 30 y+1)` == yrs[t] & 
    #                                                  INDIVIDUAL_ID$IdSimplified == ALLIDS[i]])
    #     if(length(tmp) > 0){
    #       if(length(tmp) > 1){
    #         # print(tmp)
    #         tmp <- tmp[1]
    #       }
    #       DATA$STATUS[DATA$IdSimplified %in% ALLIDS[i] & DATA$Year %in% yrs[t]] <- tmp
    #     }
    #   }#t
    # }#i
    # # 
    # # tmp <- DATA[DATA$IdSimplified %in% ALLIDS[i] & DATA$Year %in% years[t],]
    # # tmp1 <- INDIVIDUAL_ID[INDIVIDUAL_ID$`ReprodYear (May 1 year y - Apr 30 y+1)` == yrs[t] & 
    # #                                INDIVIDUAL_ID$IdSimplified == ALLIDS[i],]
    # # tmp1$`ReprodYear (May 1 year y - Apr 30 y+1)`
    # # 
    # numOverwiteSex <- sum(unique(INDIVIDUAL_ID$Id) %in% DATA$Id)
    # #overwrite Sex with pack id info
    # for(i in 1:nrow(Pack_ID2023)){
    #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2023$IdSimplified[i]] <- Pack_ID2023$Sex[i]
    #   DATA$STATUS[DATA$IdSimplified %in% Pack_ID2023$IdSimplified[i] & 
    #                 DATA$Year %in% 2022 ] <- 3
    #   
    # }#i
    # 
    # ##-- Overwrite sex 
    # for(i in 1:nrow(Pack_ID2024)){
    #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2024$IdSimplified[i]] <- Pack_ID2024$Sex[i]
    #   DATA$STATUS[DATA$IdSimplified %in% Pack_ID2024$IdSimplified[i] & 
    #                 DATA$Year %in% 2023 ] <- 3
    # }#i
    # 
    # ##-- Overwrite sex 
    # for(i in 1:nrow(Pack_ID2025)){
    #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2025$IdSimplified[i]] <- Pack_ID2025$Sex[i]
    #   DATA$STATUS[DATA$IdSimplified %in% Pack_ID2025$IdSimplified[i] & 
    #                 DATA$Year %in% 2024 ] <- 3
    # }#i
    # 
    # for(i in 1:nrow(Pack_ID2026)){
    #   DATA$Sex[DATA$IdSimplified %in% Pack_ID2026$IdSimplified[i]] <- Pack_ID2026$Sex[i]
    #   DATA$STATUS[DATA$IdSimplified %in% Pack_ID2026$IdSimplified[i] & 
    #                 DATA$Year %in% 2025 ] <- 3
    # }#i
    
    ##statut[CM] assignation of state to the y array
    INDIVIDUAL_ID$STATUS_Numeric <- 1
    INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% "Juvenile"] <- 2
    INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Pair" )] <- 3
    INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$Status %in% c("Family group")] <- 4
    
    ## FOR THE LAST YEAR GET THROUGH THE PAIR BASED-FILE FROM LINN.
    indIDRovBase <- unlist(lapply(strsplit(y.ar$Id.vector," "), function(x) x[1]))
    
    ALLIDS <- c(unique(indIDRovBase))
    
    y.obsALL <- matrix(1, nrow = length(ALLIDS), ncol = dim(y.ar$y.ar)[3]+1)
    yrs <- c(years[1]-1, years)
    dimnames(y.obsALL) <- list(ALLIDS, yrs)
    for(i in 1:dim(y.obsALL)[1]){
      for(t in 1:(length(years)+1)){
        tmp <- unique(INDIVIDUAL_ID$STATUS_Numeric[INDIVIDUAL_ID$`ReprodYear (May 1 year y - Apr 30 y+1)` == yrs[t] & 
                                                     INDIVIDUAL_ID$IdSimplified == ALLIDS[i]])
        if(length(tmp) > 0){
          if(length(tmp) > 1){
            print(tmp)
            tmp <- tmp[1]
          }
          y.obsALL[i,t] <- tmp
        }
      }#t
    }#i

    
    ## FOR THE LAST YEAR GET TRHOUGH THE PAIR BASED-FILE FROM LINN.
    for(i in 1:dim(y.obsALL)[1]){
      t <- length(years)
      
      #linn's 2023 file
      tmp <- Pack_ID2023[Pack_ID2023$IdSimplified %in% ALLIDS[i],]
      if(nrow(tmp) > 0){
        y.obsALL[i,t-2] <- 3
      }
      
      #linn's 2024 file
      tmp <-  Pack_ID2024[Pack_ID2024$IdSimplified %in% ALLIDS[i],]
      if(nrow(tmp) > 0){
        y.obsALL[i,t-1] <- 3
      }
      
      #linn's 2025 file
      tmp <-  Pack_ID2025[Pack_ID2025$IdSimplified %in% ALLIDS[i],]
      if(nrow(tmp) > 0){
        y.obsALL[i,t] <- 3
      }
      
      #linn's 2026 file
      tmp <- Pack_ID2026[Pack_ID2026$IdSimplified %in% ALLIDS[i],]
      if(nrow(tmp) > 0){
        y.obsALL[i,t+1] <- 3
      }
    }#i
    
    ## subset y.obs for the individuals present in y.ar
    indID <- unlist(lapply(strsplit(y.ar$Id.vector, " "),
                           function(x)x[1]))
    y.obs <- y.obsALL[indID, ]
    #y.obs.family <- y.obs
    y.obs[y.obs == 4] <- 3
    y.obs <- y.obs[ ,as.character(years)]
    
    indSocialState <- matrix(1, nrow = dim(y.ar$y.ar)[1], 
                             ncol = dim(y.ar$y.ar)[3])
    dimnames(indSocialState) <- list(dimnames(y.ar$y.ar)[[1]],
                                     years) 
    for(i in 1:dim(indSocialState)[1]){
      if(any(y.obs[i, ] >= 3)){
        indSocialState[i, min(which(y.obs[i, ] >= 3)):dim(y.ar$y.ar)[3]] <- 2
      }
    }#i
    
    y.status <- indSocialState
    
  # [CM] Commented out 
  #   # data.alive$data.sp$Status <- data.alive$data.sp$STATUS
  #   ##-- Turn status into factor for correct order
  #   # data.alive$data.sp$Status <- factor( data.alive$data.sp$Status,
  #   #                                      levels = c("Juvenile", "Pair", "Family group"))
  #   # 
  #   ##-- Make a table of individuals states per year
  #   ##-- 1: no info
  #   ##-- 2: Juvenile
  #   ##-- 3: Pair
  #   ##-- 4: Family group
  #   # tmp <- apply( table(data.alive$data.sp$Id,
  #   #                     data.alive$data.sp$Status,
  #   #                     data.alive$data.sp$Year),
  #   #               c(1,3),
  #   #               function(x)ifelse(any(x>0), which(x>0), 0)) #+ 1
  #   
  #   ##-- Resize to match detection array
  #   y.obs <- y.status <- matrix(1, nrow = nrow(y.ar.DEAD), ncol = ncol(y.ar.DEAD))
  #   dimnames(y.obs) <- dimnames(y.ar.DEAD)
  #   y.obs[dimnames(tmp)[[1]], ] <- tmp
  #   y.obs[y.obs > 2] <- 3
  #   
  #   ##-- For the model, set status to 2 from the first time it is "pair" or "family group" to the last occasion
  #   dimnames(y.status)[1] <- list(dimnames(y.ar$y.ar)[[1]]) 
  #   
  #   for(i in 1:dim(y.status)[1]){
  #     if(any(y.obs[i, ] >= 3)){
  #       y.status[i, min(which(y.obs[i, ] >= 3)):ncol(y.status)] <- 2
  #     }
  #   }#i

  
    
    
    ## ------       7.4.2. TRAP-RESPONSE ------ 
    
    ##-- Make matrix of previous capture indicator
    detResponse <- makeTrapResponseCov(
      data = myFullData.sp$alive,
      data.dead = myFullData.sp$dead.recovery)
    
    ##-- Subset to focal years
    detResponse <- detResponse[ ,dimnames(detResponse)[[2]] %in% dimnames(y.ar$y.ar)[[3]]]
    
    ##-- Subset to focal individuals
    detResponse <- detResponse[dimnames(y.ar$y.ar)[[1]], ]
    
    
    
    ## ------     7.5. HAB DENSITY ------ 
    
    ##-- KERNEL OF INDIVIDUALS IN PAIRS
    #[CM] Commented out 
    #[CM] use the what we had in the previous script
    # kern <- list()
    # habDens <- matrix(NA, nrow = n.habWindows, ncol = n.years)
    # for(t in 1:n.years){
    #   ##-- Subset the NGS data to individuals in packs/pairs this year
    #   data.pairs.t <- myFullData.sp$alive %>%
    #     dplyr::filter( Year == years[t],
    #                    STATUS %in% c(3,4),
    #                    Sex=="male")
    #   
    #   ##-- Get mean coordinates of packs
    #   IDs <- unique(data.pairs.t$IdSimplified)
    #   m.xy <- matrix(NA, nrow = length(IDs), ncol = 2)
    #   colnames(m.xy) <- c("x","y")
    #   for(i in 1:length(IDs)){
    #     m.xy[i, ] <- data.pairs.t %>%
    #       dplyr::filter( IdSimplified == IDs[i]) %>%
    #       st_coordinates(.) %>%
    #       colMeans(.)
    #   }#i
    #   
    #   ##-- Check if some coordinates are missing  
    #   if(sum(is.na(m.xy[ ,1])) > 0){m.xy <- m.xy[!is.na(m.xy[ ,1]), ]}
    #   
    #   ##-- Turn into .sf
    #   locationsFamily <- st_as_sf( as.data.frame(m.xy),
    #                                coords = c("x","y"),
    #                                crs = st_crs(habitat$habitat.sp))
    #   locationsFamily$id <- rep(1, nrow(locationsFamily))
    #   
    #   ##-- Calculate kernel of detections
    #   kern[[t]] <- raster(estUDm2spixdf(kernelUD( 
    #     as(locationsFamily[ ,"id"], "Spatial"), 
    #     h = 15000,
    #     grid = as(habitat$habitat.r, 'SpatialPixels'))))
    #   
    #   ##-- Plot check
    #   plot(kern[[t]], main = years[t])
    #   plot(habitat$habitat.poly$geometry, add = T, col = NA)
    #   
    #   ##-- Scale covariate
    #   habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ] == 1])
    # } #t
    
    
    ##-- KERNEL OF INDIVIDUALS IN PAIRS
    #[CM] Future improvements => use both sex to construct the map.
    #[CM] few diffs because maps are constructed using all data since 2012 in "original" script
    kern <- list()
    habDens <- matrix(NA, nrow = n.habWindows, ncol = n.years)
    IDS <- unlist(lapply(strsplit(as.character(myFullData.sp$alive$Id) , " "), function(x)x[1])) 
    for(t in 1:n.years){
      id.fam <- which(y.obsALL[ ,as.character(years[t]-1)] 
                      %in%
                        c(3,4), arr.ind = T)
      
      ## [PD] ADDED THE IF STATEMENT HERE AS A BANDAID UNTIL WE HAVE THE FINAL FILE FROM LINN
      if(length(id.fam) > 0){
        m.xy <- matrix(NA, nrow = length(id.fam), ncol = 2)
        colnames(m.xy) <- c("x","y")
        for(i in 1:length(id.fam)){
          tmp <- myFullData.sp$alive[IDS == row.names(y.obsALL)[i] & 
                                       myFullData.sp$alive$Sex %in% thisSex, ]
          m.xy[i, ] <- colMeans(st_coordinates(tmp))
        }
        if(sum(is.na(m.xy[,1]))>0){
          m.xy <- m.xy[!is.na(m.xy[ ,1]), ]
        }
      }
      locationsFamily <- st_as_sf( as.data.frame(m.xy),
                                   coords = c("x","y"),
                                   crs = st_crs(habitat$habitat.sp))
      locationsFamily$id <- rep(1, nrow(locationsFamily))
      kern[[t]] <- raster(estUDm2spixdf(kernelUD( as(locationsFamily[ ,"id"],"Spatial"),
                                                  h = 15000,
                                                  grid = as(habitat$habitat.r, 'SpatialPixels'))))
      habDens[ ,t] <- scale(kern[[t]][habitat$habitat.r[ ]==1])
    }#t
    
    # ##-- Check 
    # for(t in 1:n.years){
    #   plot(kern[[t]], main = years[t])
    #   plot(habitat$habitat.poly$geometry, add = T, col = NA)
    # }#t
     
    
    
    ## ------     7.6. AUGMENT DETECTION HISTORIES -----
    
    ##-- Data arrays
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
    
    ##-- Individual covariates
    y.status <- makeAugmentation( y = y.status,
                                  aug.factor = aug.factor,
                                  replace.value = 1)
    
    detResponse <- makeAugmentation( y = detResponse,
                                     aug.factor = aug.factor,
                                     replace.value = 0)
    ##-- Set first detection for augmented individuals to NA
    detResponse[rownames(detResponse) %in% "Augmented",1]  <- NA
    
    
    
    ## ------     7.7. TRANSFORM Y TO SPARSE MATRICES ------
    
    ##-- STRUCTURED
    y.sparse <- nimbleSCR::getSparseY(y.aliveStructured)
    
    ##-- OTHER
    y.sparseOth <- nimbleSCR::getSparseY(y.aliveOthers)
    
    
    
    ## ------ IV. MODEL SETTING ------- 
    
    ## ------   1. NIMBLE MODEL DEFINITION ------
    
    modelCode <- nimbleCode({
      
      ##------ SPATIAL PROCESS ------## 
      for(st in 1:2){
        dmean[st] ~ dunif(0,100)
        lambda[st] <- 1/dmean[st]
      }#st
      
      betaDens ~ dnorm(0.0,0.01)
      
      for(t in 1:n.years){
        habIntensity[1:n.habWindows,t] <- exp(betaDens * habDens[1:n.habWindows,t])
        sumHabIntensity[t] <- sum(habIntensity[1:n.habWindows,t])
        logHabIntensity[1:n.habWindows,t] <- log(habIntensity[1:n.habWindows,t])
        logSumHabIntensity[t] <- log(sumHabIntensity[t])
      }#t
      
      for(i in 1:n.individuals){
        sxy[i,1:2,1] ~ dbernppAC(
          lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
          upperCoords = upperHabCoords[1:n.habWindows,1:2],
          logIntensities = logHabIntensity[1:n.habWindows,1],
          logSumIntensity = logSumHabIntensity[1],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max)
        
        for(t in 2:n.years){
          sxy[i,1:2,t] ~ dbernppACmovement_exp(
            lowerCoords = lowerHabCoords[1:n.habWindows,1:2],
            upperCoords = upperHabCoords[1:n.habWindows,1:2],
            s = sxy[i,1:2,t-1],
            lambda = lambda[state[i,t-1]+1],
            baseIntensities = habIntensity[1:n.habWindows,t],
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            numGridRows = y.max,
            numGridCols = x.max,
            numWindows = n.habWindows)
        }#t
      }#i
      
      
      ##----- DEMOGRAPHIC PROCESS -----##
      ##-- FIRST YEAR
      omeg1[1:3] ~ ddirch(alpha[1:3])  
      
      for(i in 1:n.individuals){
        z[i,1] ~ dcat(omeg1[1:3])
        isAlive[i,1] <- (z[i,1] == 2) + (z[i,1] == 3)
        state[i,1] <- (z[i,1] == 3)
      }#i
      
      ##-- FOLLOWING YEARS 
      for(t in 1:(n.years-1)){
        gamma[t] ~ dunif(0,1)
        psi[t] ~ dunif(0,1)
        
        for(st in 1:2){
          w[st,t] ~ dunif(0,1)
          h[st,t] ~ dunif(0,1)
          rw[st,t] ~ dunif(0,1)
          ones.dead.legal[st,t] ~ dbern(step(1 - (h[st,t] + w[st,t] + rw[st,t])))     
          phi[st,t] <- 1 - h[st,t] - w[st,t] - rw[st,t]
          wAll[st,t] <- w[st,t] + rw[st,t]
        }#st
        
        omega[1,1:6,t] <- c(1-gamma[t], gamma[t]           , 0              , 0     , 0      , 0     ) ## "UNBORN"
        omega[2,1:6,t] <- c(0         , phi[1,t]*(1-psi[t]), phi[1,t]*psi[t], h[1,t], rw[1,t], w[1,t]) ## "NON-PAIRS"
        omega[3,1:6,t] <- c(0         , 0                  , phi[2,t]       , h[2,t], rw[2,t], w[2,t]) ## "PAIRS"
        omega[4,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "NEWLY DEAD LEGAL HUNTING"
        omega[5,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "NEWLY DEAD OTHER SOURCES"
        omega[6,1:6,t] <- c(0         , 0                  , 0              , 0     , 0      , 1     ) ## "DEAD"
        
        for(i in 1:n.individuals){
          z[i,t+1] ~ dcat(omega[z[i,t],1:6,t])
          isAlive[i,t+1] <- (z[i,t+1] == 2) + (z[i,t+1] == 3)
          state[i,t+1] <- (z[i,t+1] == 3)
        }#i                                                                                                        
      }#t
      
      
      ##----- DETECTION PROCESS -----##
      pResponse ~ dunif(0,1)
      for(i in 1:n.individuals){
        detResponse[i,1] ~ dbern(pResponse)
      }
      
      for(t in 1:n.years){
        for(st in 1:2){
          sigma[st,t] ~ dunif(0,50)
        }#st
        
        ##-- Structured
        betaResponse[t] ~ dunif(-5,5)
        for(n in 1:n.covs){
          betaCovs[n,t] ~ dunif(-5,5)
        }#n
        for(c in 1:n.counties){
          p0[c,1,t] ~ dunif(0,1)
          p0[c,2,t] ~ dunif(0,1)
        }#c    
        
        ##-- Opportunistic
        betaResponseOth[t] ~ dunif(-5,5)
        for(n in 1:n.covsOth){
          betaCovsOth[n,t] ~ dunif(-5,5)
        }#n
        for(c in 1:n.countries){
          p0Oth[c,1,t] ~ dunif(0,1)
          p0Oth[c,2,t] ~ dunif(0,1)
        }#c
        
        for(i in 1:n.individuals){
          
          ##-- Structured
          y.alive[i,1:maxDetNums,t] ~ dbinomLocal_normalWolf(
            detNums = detNums[i,t],
            detIndices = detIndices[i,1:maxDetNums,t],
            size = size[1:n.detectors],
            p0 = p0[1:n.counties,1:2,t],
            sigma = sigma[state[i,t]+1,t],
            s = sxy[i,1:2,t],
            trapCoords = detector.xy[1:n.detectors,1:2],
            localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
            localTrapsNum = localDetNum[1:n.habWindows],
            resizeFactor = resizeFactor,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            indicator = isAlive[i,t],
            z = z[i,t]-1,
            trapCovsIntercept = detCounties[1:n.detectors],
            indCov = detResponse[i,t],
            indBeta = betaResponse[t],
            trapCovs = detCovs[1:n.detectors,1:n.covs,t],
            trapBetas = betaCovs[1:n.covs,t],
            lengthYCombined = 1)
          
          ##--- Opportunistic
          y.aliveOth[i,1:maxDetNumsOth,t] ~ dbinomLocal_normalWolf(
            detNums = detNumsOth[i,t],
            detIndices = detIndicesOth[i,1:maxDetNumsOth,t],
            size = size[1:n.detectors],
            p0 = p0Oth[1:n.countries,1:2,t],
            sigma = sigma[state[i,t]+1,t],
            s = sxy[i,1:2,t],
            trapCoords = detector.xy[1:n.detectors,1:2],
            localTrapsIndices = localDetIndices[1:n.habWindows,1:numLocalIndicesMax],
            localTrapsNum = localDetNum[1:n.habWindows],
            resizeFactor = resizeFactor,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            indicator = isAlive[i,t],
            z = z[i,t]-1,
            trapCovsIntercept = detCountries[1:n.detectors],
            indCov = detResponse[i,t],
            indBeta = betaResponseOth[t],
            trapCovs = detCovsOth[1:n.detectors,1:n.covsOth,t],
            trapBetas = betaCovsOth[1:n.covsOth,t],
            lengthYCombined = 1)
          
          ##-- Dead recoveries legal
          x.deadculled[i,t] ~ dbern(z[i,t] == 4)
          
          ##-- Dead recoveries others
          x.deadOther[i,t] ~ dbern(z[i,t] == 5)
        }#i
      }#t
      
      
      ##---------- DERIVED PARAMETERS ----------##
      for(t in 1:n.years){
        N[t] <- sum(isAlive[1:n.individuals,t])
      }#t
    })
    
    
    
    ## ------   2. NIMBLE CONSTANTS ------
    
    nimConstants <- list( 
      n.individuals = dim(y.sparse$y)[1],
      n.habWindows = nrow(habitat$scaledLowerCoords),
      n.detectors = nrow(detectors$scaledCoords),
      n.years = dim(y.sparse$y)[3], 
      n.covs = dim(detCovs)[2],
      n.covsOth = dim(detCovsOth)[2],
      n.countries = max(detCountries),
      n.counties = max(detCounties),
      resizeFactor = detectors$localObjects$resizeFactor,
      y.max = dim(detectors$localObjects$habitatGrid)[1],
      x.max = dim(detectors$localObjects$habitatGrid)[2],
      numLocalIndicesMax = detectors$localObjects$numLocalIndicesMax,
      maxDetNums = y.sparse$maxDetNums,
      maxDetNumsOth = y.sparseOth$maxDetNums)
    
    
    
    ## ------   3. NIMBLE DATA ------
    
    ## ------     3.1. GENERATE KNOWN z ------
    
    ##-- Create known values for z (2 = alive; 3 = recovered dead)
    z.data <- apply(y.alive, c(1,3), function(x) ifelse(any(x>0),2,0))
    z.dead <- apply(y.dead, c(1,2), function(x) ifelse(any(x>0),3,0))
    z.data <- ifelse( z.dead + z.data == 0, NA, z.dead + z.data)
    
    ##-- Fill in the gaps
    z.data <- t(apply(z.data, 1, function(zz){
      if(any(!is.na(zz))){
        ##-- Identify range of detections
        range.det <- range(which(!is.na(zz)))
        ##-- If not dead recovered, set state to 2 from first to last detection
        if(sum(zz == 3, na.rm = T)<1){
          zz[range.det[1]:range.det[2]] <- 2
        }
        ##-- If recovered
        if(sum(zz == 3, na.rm = T)>0){
          ##-- If another detection alive also available
          ##-- Set state to 2 from first detection to occasion before recovered
          if(sum(zz, na.rm = T)>3){
            zz[range.det[1]:(range.det[2]-1)] <- 2
          }
          
          t.recovered <- which(zz==3)
          ##-- If recovered before last occasion
          ##-- Set state to 4 from occ. after recovery to last occ
          if(t.recovered < length(zz)){
            zz[(t.recovered+1):length(zz)] <- 6
          }
          ##-- In any case set state to 2 for occ just before recovery
          if(t.recovered <= length(zz)){
            zz[(t.recovered-1)] <- 2
          }
        }
      }
      return(zz)
    }))
    
    ##-- Identify individuals dead to legal culling and to other reasons
    legal.mx <- other.mx <- matrix( 0,
                                    nrow(y.dead), ncol(y.dead),
                                    dimnames = dimnames(y.dead))
    legal.mx[data.dead$Id[data.dead$Legal], ] <- 1
    other.mx[data.dead$Id[!data.dead$Legal], ] <- 1
    
    ##-- Id culled get state 4, id not culled get 5. 
    z.data <- ifelse(z.data == 3 & other.mx %in% c(1), 5, z.data) 
    z.data <- ifelse(z.data == 3 & legal.mx %in% c(1), 4, z.data)
    
    ##-- Id alive and in pairs get state 3
    z.data <- ifelse(z.data == 2 & y.status %in% c(2), 3, z.data)
    
    
    
    ## ------     3.2. GENERATE x.dead ------
    
    x.deadculled <- x.deadOther <- matrix( 0,
                                           nrow(z.data), ncol(z.data),
                                           dimnames = dimnames(z.data))
    x.deadculled[z.data == 4] <- 1
    x.deadOther[z.data == 5] <- 1
    
    
    
    ## ------     3.3. LIST DATA ------
    
    nimData <- list( 
      z = z.data,   
      y.alive = y.sparse$y,
      detIndices = y.sparse$detIndices,
      detNums = y.sparse$detNums,
      y.aliveOth = y.sparseOth$y, 
      detIndicesOth = y.sparseOth$detIndices,
      detNumsOth = y.sparseOth$detNums,
      x.deadculled = x.deadculled,
      x.deadOther = x.deadOther,
      habitatGrid = detectors$localObjects$habitatGrid,
      habDens = habDens,
      lowerHabCoords = as.matrix(habitat$scaledLowerCoords), 
      upperHabCoords = as.matrix(habitat$scaledUpperCoords), 
      detCounties = detCounties,
      detCountries = detCountries,
      detCovs = detCovs,
      detCovsOth = detCovsOth,
      detResponse = detResponse,
      localDetIndices = detectors$localObjects$localIndices,
      localDetNum = detectors$localObjects$numLocalIndices,
      size = detectors$detectors.df$size,
      detector.xy = as.matrix(detectors$scaledCoords),
      alpha = rep(1,3),
      ones.dead.legal = array(1,c(2,dim(y.alive)[3]-1)))
    
    
    
    ## ------   4. NIMBLE INITS ------
    
    ## ------     4.1. GENERATE INITIAL z ------
    
    ##-- Create initial values for z  
    z.init <- t(apply(z.data, 1, function(zz){
      out <- zz
      out[] <- 1
      
      if(any(!is.na(zz))){
        ##-- If not legally culled, set z to 1 before first detection and to 6 after last detection
        if(sum(zz == 4, na.rm = T) < 1){
          range.det <- range(which(!is.na(zz)))
          if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
          if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 6
        }
        ##-- If alive in a pair, set state to 2 the year before
        if(sum(zz == 3, na.rm = T) > 0){
          reco.3 <- min(which(zz == 3))
          if(reco.3>1){zz[(reco.3-1)] <- 2}
        }
        ##-- If alive not in pair, set state to 1 before first detection
        if(sum(zz == 2, na.rm = T) > 0){
          reco.alive <- min(which(zz == 2))
          if(reco.alive>1){zz[1:(reco.alive-1)] <- 1}
        }
        out[] <- zz
        
        ##-- if still some NAs initialize it with 2
        if(sum(is.na(out) > 0)){
          out[is.na(out)] <- 2
        }
      }
      return(out)
    }))
    z.init[!is.na(z.data)] <- NA
    
    
    
    ## ------     4.2. GENERATE INITIAL sxy ------ 
    
    ##-- Project death to the next year and combine all detections
    AllDetections <- data.dead[ ,c("Id","Year")] %>% 
      mutate(Year = Year + 1) %>%
      filter(Year != max(Year)) %>%
      rbind(., data.alive$data.sp[ ,c("Id","Year")]) %>%
      mutate( "x" = st_coordinates(.)[,1],
              "y" = st_coordinates(.)[,2]) %>%
      st_drop_geometry() 
    
    ##-- Rescale all detections
    AllDetections <- scaleCoordsToHabitatGrid(
      coordsData = AllDetections,
      coordsHabitatGridCenter = habitat$habitat.xy,
      scaleToGrid = T)$coordsDataScaled
    
    ##-- Identify augmented individuals
    idAugmented <- which(rownames(z.data) %in%"Augmented")
    
    ##-- Generate initial values for sxy
    sxy.init <- getSInits( AllDetections = AllDetections,
                           Id.vector = y.ar$Id.vector,
                           idAugmented = idAugmented,
                           lowerCoords = as.matrix(habitat$scaledLowerCoords),
                           upperCoords = as.matrix(habitat$scaledUpperCoords),
                           habitatGrid = detectors$localObjects$habitatGrid,
                           intensity = NULL,
                           sd = 4,
                           movementMethod = "dbernppACmovement_normal")
    
    ##-- SXY data 
    sxy.data <- sxy.init
    sxy.data[ ] <- NA
    for(i in 1:length(y.ar$Id.vector)){
      for(t in 1:(dim(sxy.data)[3]-1)){
        if(z.data[i,t+1] %in% c(4,5,6)){
          sxy.data[i, ,t+1] <- sxy.init[i, ,t+1] 
          sxy.init[i, ,t+1] <- NA
        }#if
      }#t
    }#i
    
    sxy.init <- round(sxy.init,5) #-- an extreme number of decimals may cause a number to appear as an integer to Nimble, and then coincide with habitat window boundaries
    sxy.data <- round(sxy.data,5)
    
    nimData$sxy <- sxy.data
    
    
    
    ## ------     4.3. LATENT VARIABLE DET RESPONSE ------ 
    
    detResponse.init <- detResponse
    detResponse.init[is.na(detResponse)] <- rbinom(sum(is.na(detResponse)), 1,0.5)
    detResponse.init[!is.na(detResponse)] <- NA
    
    
    
    ## ------   5. NIMBLE PARAMETERS ------ 
    
    nimParams <- c("N", "lambda", "dmean", "betaDens",
                   "omeg1", "gamma", "psi", "phi", "h", "w", "wAll", "rw",
                   "pResponse", "sigma",
                   "p0", "betaResponse", "betaCovs",
                   "p0Oth", "betaResponseOth", "betaCovsOth")
    
    nimParams2 <- c("z", "sxy")
    
    
    
    ## ------   6. SAVE INPUTS ----- 
    
    for(c in 1:4){
      
      nimInits <- list( 
        "sxy" = sxy.init,
        "dmean" = runif(2,2,4),
        "z" = z.init,
        "omeg1" = c(0.5,0.25,0.25),
        "psi" = runif(dim(y.alive)[3]-1,0.1,0.7),
        "gamma" = runif(dim(y.alive)[3]-1,0,1),
        "h" = array(runif((dim(y.alive)[3]-1)*2,0.2,0.4), c(2,dim(y.alive)[3]-1)),
        "rw" = array(runif((dim(y.alive)[3]-1)*2,0.05,0.10), c(2,dim(y.alive)[3]-1)),
        "w" = array(runif((dim(y.alive)[3]-1)*2,0.2,0.4), c(2,dim(y.alive)[3]-1)),
        "p0" = array(runif(12,0,0.2), c(nimConstants$n.counties,2,dim(y.alive)[3])),
        "p0Oth" = array(runif(12,0,0.2), c(nimConstants$n.countries,2,dim(y.alive)[3])),
        "betaResponse" = runif(dim(y.alive)[3],-1,1),
        "betaResponseOth" = runif(dim(y.alive)[3],-1,1),
        "betaDens" = runif(1,-1,1),
        "betaCovs" = array(runif(nimConstants$n.covs,-1,1),c(nimConstants$n.covs,dim(y.alive)[3])),
        "betaCovsOth" = array(runif(nimConstants$n.covsOth,-1,1),c(nimConstants$n.covsOth,dim(y.alive)[3])),
        "sigma" = array(runif(2,4,8),c(2,dim(y.alive)[3])),
        "detResponse" = detResponse.init,
        "pResponse"  = runif(1,0,1))
      
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
  
  #[CM] Sum of objects obtained from the Original script
  #Male
  # lapply(nimData, function(x) sum(x,na.rm=T))
  # $detCountries
  # [1] 5727
  # $detCounties
  # [1] 12461
  # $y.max
  # [1] 48
  # $x.max
  # [1] 29¨
  #$detCountries
  # [1] 5727
  # $detCounties
  # [1] 12461
  # $y.max
  # [1] 48
  # $x.max
  # [1] 29
  # $detector.xy
  # [1] 112192
  # $lowerHabCoords
  # [1] 28813
  # $upperHabCoords
  # [1] 30623
  # $habitatGrid
  # [1] 409965
  # $habDens
  # [1] 1.085243e-14
  # $trapCovsOth
  # [1] 14801.96
  # $trapCovs
  # [1] 3.42
  # $trials
  # [1] 322400
  # $n.individuals
  # [1] 1957
  # $n.detectors
  # [1] 3224
  # $n.years
  # [1] 10
  # $numHabWindows
  # [1] 905
  # $idResponse
  # [1] 2006
  # $y.aliveOth
  # [1] -187766
  # $z
  # [1] 14791
  # $x.deadculled
  # [1] 316
  # 
  # $x.deadOther
  # [1] 57
  ## ------   8. RETURN IMPORTANT INFOS FOR REPORT ------
  
  return(list( SPECIES = "Wolf",
               engSpecies = "wolf",
               YEARS = years,
               SEX = sex,
               DATE = DATE))
}




