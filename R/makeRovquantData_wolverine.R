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
  sex = c("female","male"),
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
  resize.factor = 1)
{

  ## ------ 0. BASIC SET-UP ------
  
  ##-- Set default values for the wolverine model
  if(is.null(aug.factor)){aug.factor <- 0.8}
  if(is.null(sampling.months)){sampling.months <- list(12,1:6)}
  if(is.null(habitat.res)){habitat.res <- 20000} 
  if(is.null(buffer.size)){buffer.size <- 60000}
  if(is.null(detector.res)){detector.res <- 10000}
  if(is.null(subdetector.res)){subdetector.res <- 2000}
  if(is.null(max.det.dist)){max.det.dist <- 84000}
  if(is.null(resize.factor)){resize.factor <- 1}
  
  ##-- Set up list of Habitat characteristics
  habitat <- list( resolution = habitat.res,
                   buffer = buffer.size)
  
  ##-- Set up list of Detectors characteristics
  detectors <- list( resolution = detector.res,
                     resolution.sub = subdetector.res,
                     maxDist = max.det.dist,
                     resize.factor = resize.factor)
  
  ##-- Set up list of Data characteristics
  data <- list( sex = sex,
                aug.factor = aug.factor,
                sampling.months = sampling.months)
  
  
  
  ## ---------------------------------------------------------------------------
  
  ## ------ I. LOAD AND SELECT DATA ------
  
  ## ------   1. HABITAT DATA -----
  
  ##-- Load pre-defined habitat rasters and shapefiles
  data(COUNTRIES, envir = environment()) 
  data(COUNTIES, envir = environment()) 
  data(habitatRasters, envir = environment()) 
  data(GLOBALMAP, envir = environment()) 
  data(REGIONS, envir = environment())
  
  
  ##-- Disaggregate habitat raster to the desired resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)
  
  
  ##-- Merge counties for practical reasons
  COUNTIES_AGGREGATED <- REGIONS %>%
    mutate(id = case_when(
      county %in% c("Norrbotten") ~ 1,
      county %in% c("Västerbotten") ~ 2,
      county %in% c("Blekinge","Dalarna","Gävleborg","Gotland","Halland","Jämtland",
                    "Jönköping","Kalmar","Kronoberg","Örebro","Östergötland","Skåne",
                    "Södermanland","Stockholm","Uppsala","Värmland","Västernorrland",
                    "Västmanland","Västra Götaland") ~ 3,
      county %in% c("Agder","Akershus","Buskerud","Innlandet","Møre og Romsdal",
                    "Oppland","Oslo","Østfold","Rogaland","Vestland","Telemark",
                    "Vestfold") ~ 4,
      county %in% c("Trøndelag") ~ 5,
      county %in% c("Finnmark") ~ 6,
      county %in% c("Nordland") ~ 7,
      county %in% c("Troms") ~ 8)) %>%
    group_by(id) %>%
    summarize() %>%
    st_simplify( ., preserveTopology = T, dTolerance = 500)
  
  COUNTIES_AGGREGATED$id <- as.character(1:nrow(COUNTIES_AGGREGATED))
  
  
  
  ## ------   2. STUDY AREA ------
  
  ##-- CREATE STUDY AREA POLYGON 
  myStudyArea <- COUNTRIES %>%
    filter(ISO %in% c("NOR","SWE")) %>%
    mutate(id = 1) %>%
    group_by(id) %>% 
    summarize()
  
  
  
  # ## ------   2. DETECTORS DATA ----- 
  # 
  # ## ------     2.1. DISTANCE TO ROADS -----
  # 
  # ##-- Load map of distance to roads (1km resolution)
  # DistAllRoads <- raster::raster(file.path(data.dir,"GIS/Roads/MinDistAllRoads1km.tif"))
  # 
  # ##-- Fasterize to remove values that fall in the sea
  # r <- fasterize::fasterize(sf::st_as_sf(GLOBALMAP), DistAllRoads)
  # r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  # DistAllRoads <- r
  # 
  # 
  # 
  # ## ------     2.2. SKANDOBS ------
  # 
  # ##-- Load the last SkandObs data file
  # skandObs <- readMostRecent( 
  #   path = file.path(data.dir, "Skandobs"),
  #   extension = ".xlsx",
  #   pattern = "Skandobs")
  # 
  # ##-- Replace scandinavian characters
  # colnames(skandObs) <- translateForeignCharacters(data = colnames(skandObs))
  # 
  # skandObs <- skandObs %>%
  #   ##-- Extract important info (e.g. month, year)
  #   dplyr::mutate( date = as.POSIXct(strptime(date, "%Y-%m-%d")),
  #                  year = as.numeric(format(date,"%Y")),
  #                  month = as.numeric(format(date,"%m")),
  #                  species = stringi::stri_trans_general(species, "Latin-ASCII")) %>%
  #   ##-- Turn into spatial points object
  #   sf::st_as_sf(., coords = c("longitude","latitude")) %>%
  #   sf::st_set_crs(. , value = "EPSG:4326") %>%
  #   sf::st_transform(. ,sf::st_crs(COUNTIES))
  # 
  # 
  # 
  # ## ------     2.3. ROVBASE OBS ------
  # 
  # ##-- GET ALL SAMPLES COLLECTED (all species)
  # rovbaseObs <- readMostRecent( path = data.dir,
  #                               extension = ".xls",
  #                               pattern = "all_samples") %>%
  #   ##-- Deal with Scandinavian characters
  #   dplyr::mutate(Species = stringi::stri_trans_general(Species, "Latin-ASCII")) %>%
  #   ##-- Filter out samples without coordinates
  #   dplyr::filter( !is.na(East_UTM33),
  #                  Species %in% c("Bjorn","Fjellrev","Gaupe","Hund","Jerv","Rodrev","Ulv"),
  #                  Sample_type %in% c("Ekskrement","Har","Urin","Valpeekskrement (Ulv)",
  #                                     "Sekret (Jerv)","Saliv/Spytt")) %>%
  #   ##-- Extract important info (e.g. month, year, country of collection)
  #   dplyr::mutate( Sample_type = translateForeignCharacters(data = Sample_type),
  #                  Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
  #                  year = as.numeric(format(Date,"%Y")),
  #                  month = as.numeric(format(Date,"%m")),
  #                  country = substrRight(County,3)) %>%
  #   ##-- Turn into spatial points object
  #   sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
  #   sf::st_set_crs(. , sf::st_crs(COUNTIES))
  # 
  # 
  # 
  ## ------   3. NGS DATA -----
  
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
  data$years <- years
  n.years <- length(years)

  ##-- list years with or without sampling in Norrbotten
  yearsSampledNorrb <- c(2016:2018,2023)
  yearsNotSampled <- years[!years %in% yearsSampledNorrb]
  whichYearsNotSampled <- which(years %in% yearsNotSampled)
  
  ##-- Filter NGS samples for dates
  myFullData.sp$alive <- myFullData.sp$alive %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to monitoring period
      Month %in% unlist(sampling.months))
  
  ##-- Filter Dead recoveries for dates
  myFullData.sp$dead.recovery <- myFullData.sp$dead.recovery %>%
    ##-- Subset to years of interest
    dplyr::filter(Year %in% years)
  
  
  
  ## ---------------------------------------------------------------------------
  ## ------ II. CREATE OPSCR DATA ------
  
  ## ------   1. GENERATE HABITAT ------
  
  message("Preparing habitat characteristics... ")
  
  ## ------     1.1. GENERATE HABITAT CHARACTERISTICS ------
  
  ##-- Determine study area based on NGS detections
  ##-- Buffer NGS detections and cut to Swedish and Norwegian borders
  studyArea <- myFullData.sp$alive %>%
    sf::st_buffer(., dist = habitat$buffer * 1.4) %>%
    mutate(id = 1) %>%
    group_by(id) %>% 
    summarize() %>% 
    sf::st_intersection(., COUNTRIES) %>%
    sf::st_as_sf()
  
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
    "id" = 1:habitat$n.habWindows,
    "x" = raster::coordinates(habitat$habitat.r)[isHab,1],
    "y" = raster::coordinates(habitat$habitat.r)[isHab,2])
  
  ##-- Make a spatial grid from polygon
  habitat$grid <- sf::st_as_sf( stars::st_as_stars(habitat$habitat.r), 
                                as_points = FALSE,
                                merge = FALSE) %>%
    filter( Habitat %in% 1) %>%
    st_set_crs( .,value = sf::st_crs(habitat$buffered.habitat.poly)) %>%
    mutate( id = 1:nrow(.),
            x = st_coordinates(st_centroid(.))[ ,1],
            y = st_coordinates(st_centroid(.))[ ,2]) 
  
  ##-- Study area grid from habitat raster
  habitat.rWthBufferPol <- sf::st_as_sf( 
    stars::st_as_stars(habitat$habitat.rWthBuffer), 
    as_points = FALSE,
    merge = TRUE) %>%
    filter(Habitat %in% 1)
  
  
  
  ## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------
  
  ## ------       1.2.1. DEN COUNTS ------
  
  ##-- Load the last DEN COUNT data file
  #DEN <- read.csv(file.path(dir.dropbox, "DATA/RovbaseData/ROVBASE DOWNLOAD 20241023/DEN_COUNTS_2009_2024_fromHB.csv"), fileEncoding="latin1")
  DEN <- readMostRecent( 
    path = data.dir,
    extension = ".csv",
    pattern = "DEN_COUNTS") %>%
    st_as_sf(., coords = c("UTM33_X", "UTM33_Y")) %>%
    st_set_crs(value = st_crs(myFullData.sp$alive)) %>%
    mutate(id = 1)
  
  # colnames(DEN) <- translateForeignCharacters(dat = colnames(DEN))
  # DEN.sp <- st_as_sf(DEN, coords = c("UTM33_X", "UTM33_Y"))
  # st_crs(DEN.sp) <- st_crs(data.alive)
  # DEN.sp$id  <- rep(1, nrow(DEN.sp))
  # DEN.sp <- DEN.sp[ ,("id")]
  
  DEN.r <- raster(
    adehabitatHR::estUDm2spixdf(
      adehabitatHR::kernelUD( as(DEN[ ,"id"], "Spatial"),
                              h = 30000,
                              grid = as(habitat$habitat.r, 'SpatialPixels'))))
  
  ##-- Plot check
  if(plot.check){
    plot(DEN.r)
    plot(st_geometry(myStudyArea), add = TRUE, border = "black")
  }
  
  ##-- EXTRACT COVARIATES
  denCounts <- DEN.r[habitat$habitat.r[ ] == 1]
  denCounts <- as.vector(round(scale(denCounts), digits = 2))
  
  
  
  ## ------   2. GENERATE DETECTORS -----
  
  message("Preparing detectors characteristics... ")
  
  ## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----
  
  ##-- Generate NGS detectors based on the study area 
  detectors$subdetectors.r <- disaggregate(
    habitat$habitat.rWthBuffer,
    fact = res(habitat$habitat.r)[1]/detectors$resolution.sub)
  
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
           rep(1,nrow(detectors$main.detector.sp))))
  
  ##-- Make a spatial grid from detector raster
  detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
    x = detectors$raster,
    fun = function(x){x>0})) %>%
    st_set_crs(.,value = st_crs(studyArea)) %>%
    mutate( id = 1:nrow(.),
            x = st_coordinates(st_centroid(.))[ ,1],
            y = st_coordinates(st_centroid(.))[ ,2]) 
  
  ##-- Extract numbers of detectors
  n.detectors <- detectors$n.detectors <- dim(detectors$main.detector.sp)[1]
  
  
  ##-- Identify detectors in Norrbotten 
  COUNTIESAroundNorrbotten <- REGIONS %>%
    group_by(county) %>%
    summarize() %>%
    filter(county %in% c("Norrbotten","Troms","Västerbotten","Nordland","Finnmark")) %>% 
    st_simplify( dTolerance = 500)
  
  ##-- Create an index of detectors in Norrbotten
  distDetsCounties <- st_distance( detectors$main.detector.sp,
                                   COUNTIESAroundNorrbotten,
                                   byid = T)
  detsNorrbotten <- which(apply(distDetsCounties, 1, which.min) == 3)
  
  
  ##-- Plot check
  if(plot.check){
    ##-- Plot detectors in Norrbotten
    plot( st_geometry(COUNTIESAroundNorrbotten))
    plot( st_geometry(detectors$main.detector.sp),
          col = "black", pch = 16, cex = 0.3, add = T)
    plot( st_geometry(detectors$main.detector.sp[detsNorrbotten, ]),
          col = "red", pch = 16, cex = 0.5, add = T)
    
    ##-- Plot NGS detectors
    plot( st_geometry(habitat$buffered.habitat.poly),
          main = paste(detectors$n.detectors, "Detectors"),
          col = rgb(0.16,0.67,0.16, alpha = 0.3))  
    plot( st_geometry(studyArea), add = TRUE,
          col = rgb(0.16,0.67,0.16, alpha = 0.5))
    plot( st_geometry(detectors$main.detector.sp),
          col = "red", pch = 16, cex = 0.1, add = TRUE)
    plot( st_geometry(COUNTRIES), add = TRUE)
  }
  
  
  
  ## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----
  
  ## ------       2.2.1. EXTRACT COUNTIES ------
  
  ##-- Extract closest county for each detector
  detCounties <- detectors$main.detector.sp %>%
    st_distance(., COUNTIES_AGGREGATED, by_element = F) %>%
    apply(., 1, function(x) which.min(x))
  
  ##-- Create a toggle matrix to turn detection probability to 0 in Norrbotten 
  ##-- in years without sampling
  countyToggle <- matrix(1, nrow = max(detCounties), ncol = n.years)
  for(t in whichYearsNotSampled){
    countyToggle[1,t] <- 0
  }
  
  ##-- Plot check 
  if(plot.check){
    myCol <- terrain.colors(nrow(COUNTIES_AGGREGATED))
    plot(st_geometry(GLOBALMAP), col = "gray80", main = "Aggregated Counties")
    plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    plot(st_geometry(COUNTRIES), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
    plot(st_geometry(detectors$main.detector.sp[detCounties %in% 5, ]),
         col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
    plot(st_geometry(detectors$main.detector.sp),
         col = myCol[detCounties], pch = 16, cex = 0.8, add = T)
    plot(st_geometry(COUNTIES_AGGREGATED), add = TRUE)
    plot(st_geometry(detectors$main.detector.sp[detCounties %in% 1, ]),
         col = "red", pch = 16, cex = 0.8, add = T)
    text(st_geometry(COUNTIES_AGGREGATED), labels = COUNTIES_AGGREGATED$id, col = "black")  
  }
  
  
  
  ## ------       2.2.2. EXTRACT COUNTRIES ------
  
  ##-- Extract closest country for each detector
  detCountries <- detectors$main.detector.sp %>%
    st_distance(., COUNTRIES, by_element = F) %>%
    apply(., 1, function(x) which.min(x)) %>%
    as.factor(.) %>%
    as.numeric(.)
  
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
  
  ##-- Plot check 
  if(plot.check){
    par(mfrow = c(1,1))
    myCol <- c("blue4", "yellow1", "red")
    plot( st_geometry(GLOBALMAP), col = "gray80", main = "Countries")
    plot( st_geometry(myStudyArea),
          col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    plot( st_geometry(detectors$main.detector.sp),
          col = myCol[detCountries[,1]], pch = 16, cex = 0.8, add = T)
    plot( st_geometry(COUNTRIES), add = TRUE)
  }
  
  
  
  ## ------       2.2.3. EXTRACT GPS TRACKS LENGTHS ------
  
  message("Cleaning GPS tracks... ")
  
  ## LOAD NEW GPS SEARCH TRACKS !!!
  # TRACKS_SINGLE <- read_sf(file.path(data.dir,
  #                                    "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20250908.shp"))
  # TRACKS_MULTI <- read_sf(file.path(data.dir,
  #                                   "GPS/eksport_rovquant_aktivitetslogg_20250908/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20250908.shp"))
  ##-- Combine all GPS tracks
  TRACKS <- rbind(
    read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_multilinestring_20240829_dateSfAll.shp")),
    read_sf(file.path(data.dir, "GPS/XX_eksport_rovquant_aktivitetslogg_alle_spor_linestring_20240829_dateSfAll.shp"))) %>%
    ##-- Process dates
    mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
            Mth = as.numeric(format(Dato,"%m")),
            Yr = as.numeric(format(Dato,"%Y")),
            Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
    ##-- Filter out irrelevant tracks
    filter( Helikopter == "0",      ## Remove helicopter tracks
            Jerv == "1",            ## Keep Wolverine tracks only
            Year %in% years & Mth %in% unlist(sampling.months)) %>% ## Keep tracks during sampling season only
    ##-- Extract track lengths & centroids
    mutate( Length = st_length(., byid = T),
            Centroidx = st_coordinates(st_centroid(.))[ ,1]) 
  
  ##-- Find & filter out duplicates based on person, distance and date.
  df <- data.frame( Dato = TRACKS$Dato,
                    Year = TRACKS$Year,
                    Person = TRACKS$Person,
                    Length = TRACKS$Length,
                    Centroidx = TRACKS$Centroidx)
  dupIDs <- which(duplicated(df))
  dupLength <- TRACKS$Length[duplicated(df)]
  TRACKS <- TRACKS[-dupIDs, ]
  
  # ## PLOT CHECK
  # if(plot.check){
  # 
  #   par(mfrow = c(2,2))
  #   ## Length of tracks searched per year
  #   lengthPerYear <- unlist(lapply(years,function(x) sum(TRACKS$Length[TRACKS$Year == x])/1000))
  #   names(lengthPerYear) <- years
  #   barplot(lengthPerYear, ylab = "Track length (km)", main = "Length of tracks searched per year")
  # 
  #   ## Number of tracks searched per year
  #   numPerYear <- unlist(lapply(years,function(x) sum(TRACKS$Year == x)))
  #   names(numPerYear) <- years
  #   barplot(numPerYear, ylab = "Number of tracks", main = "Number of tracks searched per year")
  # 
  #   ## Length of tracks duplicated per year
  #   dupdist <- unlist(lapply(dupDist,function(x) sum(x)/1000))
  #   names(dupdist) <- years
  #   barplot(dupdist,ylab = "Track length (km)", main = "Length of tracks duplicated per year")
  # 
  #   ## Number of tracks duplicated per year
  #   dup <- unlist(lapply(dupIDs,length))
  #   names(dup) <- years
  #   barplot(dup, ylab = "Number of tracks", main = "Number of tracks duplicated per year")
  # }
  # ##-- save 
  # save( TRACKS, file = file.path(working.dir, "data/searchTracks.RData"))
  # load(file = file.path(working.dir, "data/searchTracks.RData"))
  
  
  ##-- Extract length of GPS search track per detector grid cell
  detTracks <- matrix(0, nrow = n.detectors, ncol = n.years)
  TRACKS.r <- list()
  for(t in 1:n.years){
    intersection <- TRACKS %>%
      dplyr::filter(Year == years[t]) %>%
      sf::st_intersection(detectors$grid, .) %>%
      dplyr::mutate(LEN = st_length(.)) %>%
      sf::st_drop_geometry() %>%
      group_by(id) %>%
      summarise(transect_L = sum(LEN)) ##-- Get total length searched in each detector grid cell
    detTracks[intersection$id,t] <- as.numeric(intersection$transect_L)
    TRACKS.r[[t]] <- detectors$raster
    TRACKS.r[[t]][detectors$raster[] %in% 1] <- detTracks[ ,t]
    print(t)
  }#t
  
  
  ##-- Plot check 
  if(plot.check){
    max <- max(unlist(lapply(TRACKS.r, function(x) max(x[], na.rm = T))))
    cuts <- seq(0, max,length.out = 100)   ## Set breaks
    col <- rev(terrain.colors(100))
    CountriesDetRes <- disaggregate(habitatRasters$Countries, fact = 2)
    CountriesDetRes <- crop(CountriesDetRes, TRACKS.r[[1]])
    country.colors <- c("goldenrod1","goldenrod3")
    
    pdf(file = file.path(working.dir, "figures/Tracks.pdf"))
    NORTRACKS <- SWETRACKS <- 0
    for(t in 1:n.years){
      plot( TRACKS.r[[t]], main = years[t], breaks = cuts, col = col, legend = FALSE)
      plot( st_geometry(habitat$habitat.poly), main = years[t], add = T)
      plot( TRACKS.r[[t]],
            legend.only = TRUE, breaks = cuts, col = col, legend.width = 2,
            axis.args = list(at = round(seq(0, max, length.out = 5), digits = 1),
                             labels = round(seq(0, max, length.out = 5), digits = 1),
                             cex.axis = 0.6),
            legend.args = list(text = '', side = 4, font = 2, line = 2.5, cex = 0.8))
      ##-- Summary tracks
      NORTRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 2],na.rm = T )/1000
      SWETRACKS[t] <- sum(TRACKS.r[[t]][CountriesDetRes[]%in% 4],na.rm = T )/1000
    }#t
    
    years1 <- years + 1
    plot(SWETRACKS ~ years1, col = country.colors[2],
         lwd = 2, pch = 16, type = "b",
         ylim = c(0,300000), ylab = "sum tracks km")
    lines(NORTRACKS ~ years1, col = country.colors[1],
          lwd = 2, pch = 16, type = "b")
    legend("topright",c("N","S"), fill=country.colors)
    dev.off()
  }
  
  
  
  ## ------       2.2.4. EXTRACT DISTANCES TO ROADS ------
  
  ##-- Load map of distance to roads (1km resolution)
  DistAllRoads <- raster::raster(file.path(data.dir,"GIS/Roads/MinDistAllRoads1km.tif"))
  
  ##-- Fasterize to remove values that fall in the sea
  r <- fasterize::fasterize(sf::st_as_sf(COUNTRIES), DistAllRoads)
  r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  DistAllRoads <- r
  DistAllRoads <- crop(DistAllRoads, myStudyArea)
  rm(list = c("r"))
  
  ##-- AGGREGATE TO MATCH THE DETECTORS RESOLUTION
  DistAllRoads <- aggregate( DistAllRoads,
                             fact = detectors$resolution/res(DistAllRoads),
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
  
  ##-- Plot check 
  if(plot.check){
    par(mfrow = c(1,1))
    plot( st_geometry(GLOBALMAP),
          col = "gray80", main = "Distance to roads")
    plot( st_geometry(myStudyArea),
          col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add = T)
    plot( st_geometry(COUNTRIES),
          col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add = T)
    plot( DistAllRoads, add = T)
    plot( st_geometry(detectors$main.detector.sp),
          cex = DoScale(detRoads), pch = 16, add = T)
  }
  
  
  
  ## ------       2.2.5. EXTRACT DAYS OF SNOW ------
  
  #[PD] NEW SNOW FILE FROM ASUN!
  #SNOW <- stack(paste0(dir.dropbox,"/DATA/GISData/SNOW/ModisSnowCover0.1degrees/AverageSnowCoverModisSeason2014_2025_Wolverine.tif"))
  SNOW <- stack(file.path(data.dir,"GIS/AverageSnowCoverModisSeason2008_2024_Wolf.tif"))
  
  ##-- RENAME THE LAYERS
  names(SNOW) <- paste(2008:2023, (2008:2023) + 1, sep = "_")
  
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
  
  ##-- Plot check
  if(plot.check){
    plot( st_geometry(detectors$main.detector.sp),
          cex = DoScale(detSnow[ ,6], l = 0, u = 0.5),
          pch = 16)
  }
  
  
  
  ## ------       2.2.6. EXTRACT PRESENCE OF OTHER SAMPLES ------
  
  ## ------         2.2.6.1. SKANDOBS ------
  
  ##-- Load the last SkandObs data file
  skandObs <- readMostRecent( 
    path = data.dir,
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
                   monitoring.season = ifelse( month < unlist(sampling.months)[1],
                                               year, year + 1)) %>%
    ##-- Filter based on monitoring season
    dplyr::filter( month %in% unlist(sampling.months)) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf(., coords = c("longitude","latitude")) %>%
    sf::st_set_crs(., value = "EPSG:4326") %>%
    sf::st_transform(., sf::st_crs(COUNTIES)) %>%
    dplyr::filter(!is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))
  
  
  
  ## ------         2.2.6.2. ROVBASE ------
  
  ##-- Load the last Rovbase data files
  rovbaseObs1 <- readMostRecent( 
    path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
    extension = ".xlsx",
    pattern = "RIB2810202415264376")
  rovbaseObs2 <- readMostRecent( 
    path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
    extension = ".xlsx",
    pattern = "RIB28102024152348493")
  rovbaseObs3 <- readMostRecent( 
    path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
    extension = ".xlsx",
    pattern = "RIB28102024152447860")
  rovbaseObs4 <- readMostRecent( 
    path = file.path(data.dir,"ALL SPECIES IN SEPERATE YEARS"),
    extension = ".xlsx",
    pattern = "RIB28102024152538742")
  
  
  ##-- Process Rovbase observations (all species)
  rovbaseObs <- rbind(rovbaseObs1, rovbaseObs2, rovbaseObs3, rovbaseObs4) %>%
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
      monitoring.season = ifelse(month < 12, year, year+1)) %>%
    ##-- Filter out unusable samples
    dplyr::filter( 
      ##-- Filter out samples without coordinates,...
      !is.na(East_UTM33),
      ##-- ...based on species
      # Species %in% c("Bjorn","Fjellrev","Gaupe","Hund","Jerv","Rodrev","Ulv"),
      ##-- ...based on sample type
      Sample_type %in% c( "Ekskrement","Har","Urin","Valpeekskrement (Ulv)","Sekret (Jerv)",
                          "Saliv/Spytt", "Loepeblod", "Vev"),
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months),
      ##-- ... if sample was from the focal species and successfully genotyped 
      !(Species %in% "Jerv" & !is.na(Id))) %>%
    ##-- Turn into spatial points object
    sf::st_as_sf( ., coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(. , sf::st_crs(COUNTIES)) %>%
    filter( 
      ##-- ... based on space 
      !is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))
  
  ##-- Remove un-necessary objects
  rm(list = c("rovbaseObs1","rovbaseObs2","rovbaseObs3","rovbaseObs4"))
  
  
  
  ## ------         2.2.6.3. COMBINE ROVBASE & SKANDOBS ------
  
  ##-- Rastertize at the detector level
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
    r.SkandObsRovbaseBinary[[t]][r.SkandObsRovbaseBinary[[t]][]>1 ] <- 1
  }
  
  
  ##-- Plot check
  if(plot.check){
    plot(st_geometry(habitat.rWthBufferPol))
    plot(st_geometry(skandObs), col = "red", add = T)
    
    ##-- Summary SkanbObs
    pdf(file = file.path(working.dir, "figures/skandObs.pdf"), width = 10)
    barplot(table(skandObs$monitoring.season ))
    barplot(table(skandObs$month ), xlab = "Months")
    barplot(table(skandObs$species))
    
    ##-- Maps
    par(mar = c(0,0,2,0))
    for(t in 1:n.years){
      plot( st_geometry(myStudyArea), main = years[t])
      plot( st_geometry(skandObs[skandObs$monitoring.season %in% years[t], ]),
            pch = 16, col = "red", cex = 0.1, add = T)
    }
    dev.off()
    
    # pdf(file = file.path(working.dir, "figures/mapStructuredOthers.pdf"))
    # for(t in 1:n.years){
    #   tmpOthers <- data.alive[data.alive$Year %in% years[t] &
    #                                          !data.alive$structured, ]
    #   tmpStruct <- data.alive[data.alive$Year %in% years[t] &
    #                                          data.alive$structured, ]
    # 
    #   par(mfrow=c(2,2),mar=c(0,0,5,0))
    #   plot(r.OtherSamplesBinary[[t]], main=paste(years[t],"\n Rovbase Samples Structured"), box=F, axes=F)
    #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    #   plot(r.OtherSamplesBinary[[t]],main=paste(years[t],"\n Rovbase Samples Opportunistic"), box=F, axes=F)
    #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.6,add=T)
    # 
    #   plot(r.skandObsSamplesBinary[[t]], main=paste(years[t]), box=F, axes=F)
    #   plot(st_geometry(tmpOthers), pch=16, col="blue",bg="blue", cex=0.6,add=T)
    #   plot(r.skandObsSamplesBinary[[t]],main=paste(years[t],"\n SkandObs Opportunistic"), box=F, axes=F)
    #   plot(st_geometry(tmpStruct), pch=16, col="red",bg="red", cex=0.5,add=T)
    # }
    # dev.off()
    
    for(t in 1:n.years){
      par(mfrow=c(1,3),mar=c(0,0,5,0))
      plot(r.rovbaseBinary[[t]],main=years[t])
      plot(r.skandObsBinary[[t]])
      plot(r.SkandObsRovbaseBinary[[t]])
    }#t
  }
  
  
  
  # ## ------         2.2.6.4. SMOOTH THE BINARY MAP ------
  # 
  # ##-- We tried adjust = 0.05, 0.037,0.02 and decided to go for 0.037 
  # habOwin <- spatstat.geom::as.owin(as.vector(extent(detectors$raster)))
  # cutoff <- 1
  # ds.list <- lapply(years,function(y){
  #   ## ROVBASE DATA 
  #   pts <- st_coordinates(rovbaseObs.sp)[rovbaseObs.sp$monitoring.season %in% y,]
  #   ## SKANDOBS
  #   pts <- rbind(pts, st_coordinates(skandObs)[skandObs$monitoring.season %in% y,] )
  #   ## SMOOTH AND RASTERIZE
  #   p <-  spatstat.geom::ppp(pts[,1], pts[,2], window = habOwin)
  #   ds <- density(p, adjust=0.02) #---change bandwith (smoothing) with "adjust
  #   ds <- raster(ds)
  #   
  #   ds <- ds1 <- raster::resample(ds, detectors$raster) #mask(ds,rasterToPolygons(myHabitat.list$habitat.rWthBuffer,function(x) x==1))
  #   threshold <- 0.1 / prod(res(ds)) #--number per 1 unit of the projected raster (meters)
  #   ds1[] <- ifelse(ds[]<threshold,0,1)
  #   ds1 <- mask(ds1, habitat.rWthBufferPol)
  #   ds <- mask(ds, habitat.rWthBufferPol)
  #   
  #   return(list(ds,ds1))
  # })
  # 
  # ds.brick <- brick(lapply(ds.list, function(x) x[[1]]))
  # ds.brickCont <- brick(lapply(ds.list, function(x) x[[2]]))
  # names(ds.brick) <- names(ds.brickCont) <-years
  # 
  # ##-- Plot check
  # if(plot.check){
  #   par(mfrow = c(1,3))
  #   plot(r.SkandObsRovbaseBinary[[t]], main = "Raw Binary", axes = F, box = F)
  #   plot(ds.brick[[t]], main = "Smoothed", axes = F, box = F)
  #   plot(ds.brickCont[[t]], main = "Binary after smoothing", axes = F, box = F)
  # }
  
  
  
  ## ------         2.2.6.5. COLOR CELLS WHERE HAIR TRAP COLLECTED ------
  
  ##-- IDENTIFY HAIR SAMPLES
  tmpHair <- myFullData.sp$alive %>% filter(hairTrap)
  
  ##-- MANUALLY FIND THE HAIR SMAPLES AND COLOR THE CELL.
  tmpyr <- unique(tmpHair$Year)
  for( i in 1:length(tmpyr)){
    t <- which(years %in% tmpyr)
    whereHair <- raster::extract( r.SkandObsRovbaseBinary[[t]],
                                  tmpHair,
                                  cellnumbers = T)
    r.SkandObsRovbaseBinary[[t]][whereHair[ ,1]] <- 1
  }#t
  
  
  
  ## ------         2.2.6.6. ASSIGN THE COVARIATE ------
  
  detOtherSamples <- matrix(0, nrow = n.detectors, ncol = n.years)
  detOtherSamples[ ,1:n.years] <- raster::extract( r.SkandObsRovbaseBinary,
                                                   detectors$main.detector.sp)
  colSums(detOtherSamples)
  
  
  
  ## ------       2.2.7. SCALE & ROUND DETECTOR-LEVEL COVARIATES ------
  
  detSnow <- round(scale(detSnow), digits = 2)
  detRoads <- round(scale(detRoads), digits = 2)
  detTracks <- round(scale(detTracks), digits = 2)
  
  detCovs <- array(NA, c(dim(detTracks)[1], dim(detTracks)[2], 2))
  detCovs[,,1] <- detTracks
  detCovs[,,2] <- detSnow
  
  detCovsOth <- array(NA, c(dim(detTracks)[1],dim(detTracks)[2], 3))
  detCovsOth[,,1] <- detSnow
  detCovsOth[,,2] <- matrix(detRoads,length(detRoads),n.years)
  detCovsOth[,,3] <- detOtherSamples
  
  
  ##-- CHECK IF CONTAINS NAs
  if(any(is.na(detCovs))){print("WARNINGS!!!!!!! ONE OF THE DETECTOR MATRIX CONTAINS NA")}
  
  
  ##-- Plot check
  if(plot.check){
    tmp <- detectors$raster
    par(mfrow=c(2,5),mar=c(0,0,0,0))
    max <- max(detCovsOth[,,2])
    cuts <- seq(0,max,length.out = 100) 
    col <- rev(terrain.colors(100))
    for(t in 1:n.years){
      plot(detectors$raster, col=c(grey(0.2),grey(0.8)),axes=F,legend=F,box=F,)
      tmp[!is.na(detectors$raster[ ])] <- detCovsOth[,t,2]
      plot(tmp,axes=F,legend=F,box=F,breaks = cuts, col=col,add=T)
    }
    dev.off()
    
    pdf(file = file.path(working.dir, "figures/detections over space and time.pdf"))
    for(t in 1:n.years){
      ## NGS DETECTIONS TOTAL
      tempTotal <- myFullData.sp$alive[myFullData.sp$alive$Year == years[t], ]
      NGS_TabTotal <- table(tempTotal$Country_sample)
      ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country_sample), 2, function(x) sum(x>0))
      ## ALIVE DETECTIONS INSIDE STUDY AREA/SAMPLING PERIOD
      tempIn <- myFullData.sp$alive[myFullData.sp$alive$Year == years[t], ]
      NGS_TabIn <- table(tempIn$Country_sample)
      ID_TabIn <- apply(table(tempIn$Id, tempIn$Country_sample), 2, function(x) sum(x>0))
      ## PLOT NGS SAMPLES
      plot(st_geometry(GLOBALMAP), col="gray80")
      plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
      plot(st_geometry(habitat$buffered.habitat.poly), col = rgb(34/250, 139/250, 34/250, alpha = 0.2), add=T)
      points(tempTotal, pch = 21, bg = "darkred")
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
    dev.off()
  }
  
  
  
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
      Sex %in% sex,
      ##-- Filter data for space
      !is.na(as.numeric(sf::st_intersects(.,habitat.rWthBufferPol)))) 
  
  
  ## ------     6.2. DEAD RECOVERY DATA -----
  
  data.dead <- myFullData.sp$dead.recovery %>%
    dplyr::filter(
      ##-- Subset to years of interest
      Year %in% years,
      ##-- Subset to sex of interest
      Sex %in% sex,
      ##-- Filter data for space
      !is.na(as.numeric(sf::st_intersects(.,habitat.rWthBufferPol)))) 
  
  

  ## ------     6.3. FILTER OUT DETECTIONS IN NORRBOTTEN EXCEPT IN 2016:18 and 2023 ------
  
  ##-- Get Norrbotten borders
  COUNTIESNorrbotten <- COUNTIES %>%
    dplyr::filter(county %in% "Norrbotten") %>%
    dplyr::group_by(county) %>%
    dplyr::summarize()
  
  ##-- Identify detections collected in Norrbotten 
  is.Norr <- as.numeric(st_intersects(data.alive, COUNTIESNorrbotten))

  # ##-- Check how many detections are removed per year
  # table(data.alive$Year[data.alive$Year %in% yearsNotSampled & is.Norr %in% 1])
  # sum(data.alive$Year %in% yearsNotSampled & is.Norr %in% 1)
  
  ##-- Filter out detections in Norrbotten in years without sampling
  data.alive <- data.alive %>%
    filter(!(Year %in% yearsNotSampled & is.Norr %in% 1))

  
  
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
    mutate(
      ##-- Collector column was replaced by two columns, merging them now...
      Collector_role = ifelse(is.na(Collector_other_role), Collector_role, Collector_other_role),
      ##-- Identify samples collected during structured sampling 
      structured = Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
        !is.na(trackID) &
        trackDist <= distanceThreshold & 
        !hairTrap)
  
  # ##-- Plot check
  # if(plot.check){
  #   plot(REGIONS[REGIONS$county %in% "Norrbotten", ]$geometry)
  #   tmp <- data.alive %>%  filter(hairTrap)
  #   plot(tmp$geometry, add = T, col = "red", pch = 16)
  #   if(length(which(data.alive$DNAID[data.alive$structured] %in% HairTrapSamples$DNAID)) > 0){
  #     print("WARNING SAMPLES FROM HAIR TRAPS ASSIGNED TO STRUCTURED")
  #   }
  # }
  
  
  
  ## ------       6.3.3. PLOT CHECKS ------
  
  if(plot.check){
    
    ##-- Barplot of structured vs. opportunistic samples
    pdf(file = file.path(working.dir, "figures/DetectionsStructuredOppBarplot.pdf"))
    par(mfrow = c(2,1), mar = c(4,4,3,2))
    barplot( rbind(table(data.alive$Year[data.alive$structured]),
                   table(data.alive$Year[!data.alive$structured])),
             beside = T,
             ylim = c(0,3000),
             col = c(grey(0.2),grey(0.8)),
             ylab = "Number of samples")
    abline(h = seq(0, 3000, by = 500),
           lty = 2, col = grey(0.8))
    title(main = "500m threshold")
    legend("topleft", fill = c(grey(0.2),grey(0.8)),
           legend = c("Structured","Other"))
    
    ##-- Barplot of structured vs. opportunistic samples (threshold = 2000m)
    structured2000 <- data.alive$Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen") &
      !is.na(data.alive$trackID) &
      data.alive$trackDist <= 2000 & 
      !data.alive$hairTrap
    barplot( rbind(table(data.alive$Year[structured2000]),
                   table(data.alive$Year[!structured2000])),
             beside = T,
             ylim = c(0,3000),
             col = c(grey(0.2),grey(0.8)),
             ylab = "Number of samples")
    abline(h=seq(0,3000,by=500),lty=2,col=grey(0.8))
    title(main="2000m threshold")
    legend("topleft",fill=c(grey(0.2),grey(0.8)),
           legend = c("Structured","Other"))
    dev.off()
    
    
    ##-- CONSTRAIN TO SAMPLES COLLECTED "Fylkesmannen","SNO"
    tmp <- data.alive %>%
      filter(Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"))
    
    ## plot check
    pdf(file = file.path(working.dir, "figures/DetectionsStructuredOpp.pdf"))
    for(t in 1:n.years){
      par(mar = c(0,0,3,0), mfrow = c(1,3))
      
      ##-- samples with tracks
      tmpTracks <- tmp %>%
        filter( Year %in% years[t],
                !is.na(trackID),
                trackDist <= 500)
      plot( st_geometry(myStudyArea), col = "gray60", main = "Structured with track")
      plot( st_geometry(tmpTracks),
            pch = 21, col = "black",
            cex = 1, bg = "red", add = T)
      
      ##-- samples without tracks
      tmpNoTracks <- tmp %>%
        filter( Year %in% years[t],
                is.na(trackID) | trackDist > 500)
      plot( st_geometry(myStudyArea), col = "gray60", main = "Structured without track")
      plot( st_geometry(tmpNoTracks),
            pch = 21, col = "black",
            cex = 1, bg = "blue", add = T)
      
      ##-- Other samples
      tmpOpp <- data.alive %>%
        filter(!Collector_role %in% c("Statsforvalteren","Länsstyrelsen","SNO","Fylkesmannen"),
               Year%in% years[t])
      plot( st_geometry(myStudyArea), col = "gray60", main = "Other samples")
      plot( st_geometry(tmpOpp),
            pch = 21, col = "black",
            cex = 1, bg = "green", add = T)
      
      mtext(years[t], adj = -0.8, padj = 1)
    }#t
    
    ##-- Number of samples collected per year
    tab <- table(tmp$Year, tmp$trackID, useNA ="always")
    barplot( tab[ ,which(is.na(colnames(tab)))]/rowSums(tab),
             main = "% of samples from Statsforvalteren and \nSNO that cannot be assigned to a track")
    dev.off()
    
    
    ##-- plot check
    pdf( file = file.path(working.dir, "figures/OverallDetectionsDeadRecoveries.pdf"))
    plot( st_geometry(GLOBALMAP))
    plot( st_geometry(myStudyArea), add = T)
    plot( st_geometry(myFullData.sp$alive),
          pch = 16, col = "red", cex = 0.3, add = T)
    plot( st_geometry(myFullData.sp$dead.recovery),
          pch = 16, col = "blue", cex = 0.3, add = T)
    mtext(paste("Live detections", nrow(myFullData.sp$alive),
                "; ID:", nrow(unique(myFullData.sp$alive$Id))),
          line = +1)
    mtext(paste("Dead recovery:", nrow(myFullData.sp$dead.recovery)))
    dev.off()
  }
  
  
  
  ## ------     6.4. SEPARATE MORTALITY CAUSES ------
  
  ##-- Identify legal mortality causes
  MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$Death_cause))
  whichLegalCauses <- unlist(lapply(c("Lisensfelling","tamdyr","SNO","Skadefelling","Politibeslutning","menneske"),
                                    function(x)grep(x,MortalityNames)))
  legalCauses <- MortalityNames[whichLegalCauses]
  
  ##-- Identify legal dead recoveries based on mortality causes
  data.dead <- data.dead %>%
    mutate(legal = Death_cause %in% legalCauses)

  
  ##-- Plot check
  if(plot.check){
    # par(mfrow = c(1,3))
    # for(t in 1:n.years){
    #   ## DEAD RECOVERIES TOTAL
    #   tempTotal <- data.dead[data.dead$Year == years[t], ]
    #   NGS_TabTotal <- table(tempTotal$Country_sf)
    #   ID_TabTotal <- apply(table(tempTotal$Id, tempTotal$Country_sf), 2, function(x) sum(x>0))
    #   ## DEAD RECOVERIES INSIDE STUDY AREA/SAMPLING PERIOD
    #   tempIn <- data.dead[data.dead$Year == years[t], ]
    #   NGS_TabIn <- table(tempIn$Country_sf)
    #   ID_TabIn <- apply(table(tempIn$Id, tempIn$Country_sf), 2, function(x) sum(x>0))
    #   ## PLOT NGS SAMPLES
    #   plot(st_geometry(GLOBALMAP), col="gray80")
    #   plot(st_geometry(myStudyArea), col = rgb(34/250, 139/250, 34/250, alpha = 0.5), add=T)
    #   plot(st_geometry(tempIn), pch = 21, bg = "blue",add=T)
    #   ## ADD NUMBER OF NGS samples and IDs per COUNTRY
    #   graphics::text(x = 100000, y = 7250000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="N"], "IDs"), cex = 1.1, col = "firebrick3", font = 2)
    #   graphics::text(x = 820000, y = 6820000, labels = paste(ID_TabTotal[names(NGS_TabTotal)=="S"], "IDs"), cex = 1.1, col = "navyblue", font = 2)
    #   ## ADD OVERALL NUMBERS
    #   mtext(text = years[t], side = 3, line = 1, cex = 1.5, font = 2)
    #   mtext(text = paste(sum(NGS_TabIn), "Dead Recoveries /", sum(ID_TabIn), "IDs IN"), side = 3, line = 0)
    #   # mtext(text = paste(sum(NGS_TabTotal), "Recoveries /", sum(ID_TabTotal)-sum(ID_TabIn), "IDs OUT"), side = 3, line = -1)
    # }#t
     
    
    ##-- Plot dtemporal trends
    
    ##-- Number of detections
    pdf(file = file.path(working.dir, "figures/TRENDDetections.pdf"))
    
    temp <- unique(data.alive[ ,c("Year","Country_sf","DNAID")])
    tab_Country.Year <- table(temp$Year, temp$Country_sf)
    country.colors <- c("goldenrod1","goldenrod3")
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year))), ylim=c(0,max(tab_Country.Year)), ylab="N Detections", xlab="Years")
    lines(tab_Country.Year[,"(N)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year[,"(S)"]~as.numeric(row.names(tab_Country.Year)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("bottomright",c("N","S"), fill=country.colors)
    
    ##-- Number of IDs detected
    temp <- table(data.alive$Year,data.alive$Country_sf,data.alive$Id)
    tab_Country.Year1 <- apply(temp,c(1,2),function(x) sum(x>0))
    country.colors <- c("goldenrod1","goldenrod3")
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10, xlim=range(as.numeric(row.names(tab_Country.Year1))), ylim=c(0,max(tab_Country.Year1)), ylab="N Id detected", xlab="Years")
    lines(tab_Country.Year1[,"(N)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[1], lwd=2, pch=16, type="b")
    lines(tab_Country.Year1[,"(S)"]~as.numeric(row.names(tab_Country.Year1)), col=country.colors[2], lwd=2, pch=16, type="b")
    legend("bottomright",c("N","S"), fill=country.colors)
    
    ##-- Average number of detection per detected ID
    tab_Country.Year2 <- tab_Country.Year/tab_Country.Year1
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10,
         xlim = range(as.numeric(row.names(tab_Country.Year2))),
         ylim = c(0,max(tab_Country.Year2)),
         ylab = "Average Number of detections", xlab="Years")
    lines(tab_Country.Year2[ ,"(N)"] ~ as.numeric(row.names(tab_Country.Year2)),
          col = country.colors[1], lwd = 2, pch = 16, type = "b")
    lines(tab_Country.Year2[ ,"(S)"] ~ as.numeric(row.names(tab_Country.Year2)),
          col = country.colors[2], lwd = 2, pch = 16, type = "b")
    legend("bottomright", c("N","S"), fill = country.colors)
    
    ##-- Dead recoveries
    temp <- unique(data.dead[,c("Year","Country_sf","Id")])
    tab_Country.Year <- table(temp$Year, temp$Country_sf)
    country.colors <- c("goldenrod1","goldenrod3")
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(-10,
         xlim = range(as.numeric(row.names(tab_Country.Year))),
         ylim = c(0,max(tab_Country.Year)),
         ylab = "N Id Dead recovered", xlab="Years")
    lines(tab_Country.Year[ ,"(N)"] ~ as.numeric(row.names(tab_Country.Year)),
          col = country.colors[1], lwd = 2, pch = 16, type = "b")
    lines(tab_Country.Year[ ,"(S)"] ~ as.numeric(row.names(tab_Country.Year)),
          col = country.colors[2], lwd = 2, pch = 16, type = "b")
    legend("topright", c("N","S"), fill = country.colors)
    dev.off()
  }
  
  
  
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
  
  ##-- SHOULD NOT BE ANY INDIVIDUAL DETECTED IN NORRBOTTEN NOW 
  sum(data.alive$data.sp$sub.detector[!data.alive$data.sp$Year %in% yearsSampledNorrb] %in% subDetsNorrbotten)
  sum(data.alive$data.sp$Detector[!data.alive$data.sp$Year %in% yearsSampledNorrb] %in% detsNorrbotten)
  
  

  ## ------     6.9. SAVE FILTERED DATA ----- 

  save( data.alive, data.dead,
        file = file.path( working.dir, "data",
                          paste0("FilteredData_wolverine_", DATE, ".RData")))
  
  
  
  ## ------   7. GENERATE DETECTION HISTORY ------
  
  for(thisSex in sex){
    
    message(paste0("Preparing individual detection histories for sex: ", thisSex, "... "))
    
    ## ------     7.1. FILTER DATA BY SEX -----
    
    load(file.path( working.dir, "data",
                    paste0("FilteredData_wolverine_", DATE, ".RData")))
    
    data.alive$myData.sp <- data.alive$myData.sp %>%
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
    
    ##-- Create binary dead recovery histories (0: not recovered ; 1: recovered dead)
    y.ar.DEAD <- apply(y.ar.DEADProjected, c(1,3), function(x){as.numeric(sum(x)>0)})
    dimnames(y.ar.DEAD) <- list(dimnames(y.ar$y.ar2)[[1]], dimnames(y.ar$y.ar2)[[3]])
    y.ar.DEAD[y.ar.DEAD > 0] <- 1

    
    
    ## ------     6.8. CHECK DISTANCES BETWEEN DETECTIONS WITHIN A YEAR ------
      
    distances <- list()
    for(t in 1:n.years){
      
      ##[PD] WE NEED TO DISCUSS THE maxDist USED HERE 
      ## MUCH SMALLER THAN THE MAX DIST USED IN THE LOCAL EVAL
      
      ##-- Identify detections further than maxDist
      print(paste0("------ ", t ," -------"))
      distances[[t]] <- checkDistanceDetections( 
        y = y.ar$y.ar[ , ,t], 
        detector.xy = detectors$detectors.df[ ,c("x","y")], 
        max.distance = 40000,
        method = "pairwise",
        plot.check = F)
      
      ##-- REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
      y.ar.ALIVE[ , ,t] <- y.ar.ALIVE[ , ,t] * (1-distances[[t]]$y.flagged)
    }#t
      
      ##-- Remove detections that are further then the threshold
      #y.ar.ALIVE[,,t] <- y.ar.ALIVE[,,t] * (1-distances[[t]]$y.flagged)
      y.ar.ALIVEOthers[ , ,t] <- y.ar.ALIVEOthers[ , ,t] * (1-distances[[t]]$y.flagged)
      y.ar.ALIVEStructured[ , ,t] <- y.ar.ALIVEStructured[ , ,t] * (1-distances[[t]]$y.flagged)
      
      ##-- Remove detections also in data.alive$data.sp to run getSInits later
      idd <- names(affected.ids)
      for(i in 1:length(idd)){
        detIds <- which(distances[[t]]$y.flagged[idd[i], ] > 0)
        data.alive <- data.alive %>%
          dplyr::filter(!(Id %in% idd[i] & Detector %in% detIds & Year %in% years[t]))
      }#i
    }#t
    
    
  
      
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
      
      ##-- REMOVE DETECTIONS THAT ARE FURTHER THAN THE THRESHOLD
      y.ar.ALIVE[ , ,t] <- y.ar.ALIVE[ , ,t] * (1-distances[[t]]$y.flagged)
    }#t
    
    
    
    ## ------     7.4. AUGMENT DETECTION HISTORIES -----
    
    y.alive <- makeAugmentation( 
      y = y.ar.ALIVE,
      aug.factor = data$aug.factor,
      replace.value = 0)
    
    y.dead.ar <- makeAugmentation( 
      y = y.ar.DEAD,
      aug.factor = data$aug.factor,
      replace.value = 0)
    
    
    
    ## ------     7.5. TRANSFORM Y TO SPARSE MATRICES -----
    
    y.sparse <- nimbleSCR::getSparseY(y.alive)
    
    
    
    ## ------ IV. MODEL SETTING ------- 
    
    ## -----    1. NIMBLE CODE ------
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
    zMonths <- makeZfromScratch( 
      data.alive = myFullData.sp$alive,
      data.dead = myFullData.sp$dead.recovery,
      samplingMonths = unlist(sampling.months))
    
    ##-- Subset to focal years
    zMonths <- zMonths[ , ,dimnames(zMonths)[[3]] %in% dimnames(y.alive)[[3]]]
    
    ##-- Subset to focal individuals
    zMonths <- zMonths[dimnames(zMonths)[[1]] %in% dimnames(y.alive)[[1]], , ]
    
    ##-- Augment zMonths
    zMonths <- makeAugmentation( y = zMonths,
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
                              paste0("nimbleInput_", DATE, "_", thisSex, "_", c, ".RData")))
    }#c
    
  }#thisSex
  
  
  
  
  ## ------   8. RETURN IMPORTANT INFOS FOR REPORT ------
  return(list( SPECIES = "Brown bear",
               engSpecies = "bear",
               YEARS = years,
               SEX = sex,
               DATE = DATE))
}

