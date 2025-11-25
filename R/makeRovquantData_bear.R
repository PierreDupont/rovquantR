#' @title RovQuant OPSCR brown bear data preparation.
#'
#' @description
#' \code{makeRovquantData_bear} formats the available bear data for the OPSCR analysis using nimble and nimbleSCR.
#' The data preparation process is composed of three main steps:
#'  - defining and formatting habitat characteristics
#'  - defining and formatting detectors characteristics
#'  - defining and formatting individual detection histories
#'
#' @name makeRovquantData_bear
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
#' @importFrom graphics mtext layout 
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
#' @importFrom grDevices grey
#' 
NULL
#' @rdname makeRovquantData_bear
#' @export
makeRovquantData_bear <- function(
    ##-- paths
  data.dir = getwd(),
  working.dir = getwd(),
  
  ##-- data
  years = NULL,
  sex = c("female","male"),
  aug.factor = 2,
  sampling.months = list(c(4,5,6,7,8,9,10,11)),
  
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 50000,
  max.move.dist = 300000,
  
  ##-- detectors
  detector.res = 5000,
  subdetector.res = 1000,
  max.det.dist = 60000,
  resize.factor = 1 
){
  ## ---------------------------------------------------------------------------
  
  ## ------ 0. BASIC SET-UP ------
  
  ##-- Set default values for the brown bear model
  if(is.null(aug.factor)){aug.factor <- 2}
  if(is.null(sampling.months)){sampling.months <- list(c(4,5,6,7,8,9,10,11))}
  if(is.null(habitat.res)){habitat.res <- 20000} 
  if(is.null(buffer.size)){buffer.size <- 50000}
  if(is.null(max.move.dist)){max.move.dist <- 250000}
  if(is.null(detector.res)){detector.res <- 5000}
  if(is.null(subdetector.res)){subdetector.res <- 1000}
  if(is.null(max.det.dist)){max.det.dist <- 60000}
  if(is.null(resize.factor)){resize.factor <- 1}
  
  
  ##-- Set up list of Habitat characteristics
  habitat <- list( resolution = habitat.res,
                   buffer = buffer.size,
                   maxDist = max.move.dist)
  
  ##-- Set up list of Detectors characteristics
  detectors <- list( resolution = detector.res,
                     resolution.sub = subdetector.res,
                     maxDist = max.det.dist)
  
  ##-- Set up list of Data characteristics
  DATA <- list( sex = sex,
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
  
  ##-- Disaggregate habitat raster to the desired resolution
  habRaster <- raster::disaggregate(
    x = habitatRasters[["Habitat"]],
    fact = raster::res(habitatRasters[["Habitat"]])/habitat.res)
  
  ##-- Merge Norwegian counties for practical reasons
  COUNTIES <- COUNTIES %>%
    mutate(county = case_when(
      county %in% c("Trøndelag", "Nordland") ~ "Trøndelag",
      county %in% c("Troms", "Finnmark") ~ "Finnmark",
      county %in% c("Akershus","Agder", "Buskerud",
                    "Innlandet", "Hordaland",
                    "Møre og Romsdal","Oslo", "Oppland",
                    "Rogaland", "Vestland","Telemark",
                    "Vestfold","Østfold") ~ "Innlandet")) %>%
    dplyr::filter( , county %in% c("Trøndelag","Innlandet","Finnmark")) %>%
    dplyr::group_by(county) %>%
    dplyr::summarise() 
  
  COUNTIES$id <- as.character(1:nrow(COUNTIES))
  
  
  
  ## ------   2. NGS DATA -----
  
  ##-- Extract date from the last cleaned DATA file
  DATE <- getMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_bear")
  
  ##-- Load the most recent Bear data from RovBase
  myFullData.sp <- readMostRecent( 
    path = file.path(working.dir,"data"),
    pattern = "CleanData_bear",
    extension = ".RData")
  
  ##-- Define years
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
  
  ##-- Determine study area based on NGS detections
  ##-- Buffer NGS detections and cut to Swedish and Norwegian borders
  studyArea <- myFullData.sp$alive %>%
    sf::st_buffer(., dist = habitat$buffer) %>%
    sf::st_union() %>%
    sf::st_intersection(., COUNTIES) %>%
    sf::st_as_sf()
  
  ##-- Make habitat from predefined Scandinavian raster of suitable habitat
  habitat <- makeHabitatFromRaster(
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
  habitat$grid <- sf::st_as_sf( stars::st_as_stars(habitat$habitat.r), 
                                as_points = FALSE,
                                merge = FALSE) %>%
    filter( Habitat %in% 1) %>%
    mutate( id = 1:nrow(.)) %>%
    st_set_crs( ., value = sf::st_crs(habitat$buffered.habitat.poly))
  
  
  
  ## ------     1.2. GENERATE HABITAT-LEVEL COVARIATES ------
  
  ## ------       1.2.1. DEAD RECOVERIES (ALL YEARS) ------
  
  kern <- sf::st_coordinates(myFullData.sp$dead.recovery) %>%
    sp::SpatialPoints( .,  
                       proj4string = sp::CRS(raster::proj4string(habitat$habitat.r)),
                       bbox = NULL) %>%
    adehabitatHR::kernelUD( .,
                            h = 60000,
                            grid = as( habitat$habitat.r,
                                       'SpatialPixels')) %>%
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
                                          #"dead.reco" = habDens1,
                                          "dead.reco.trunc" = habDens2)
  
  
  
  ## ------       1.2.2. NUMBER OF OBSERVATIONS ------
  
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
  habitat$grid <- dplyr::left_join(
    x = habitat$grid,
    y = habitat$habitat.df,
    by = "id")
  
  
  
  ## ------       1.2.3. PLOTS -----
  
  ##-- [PD]: NEED TO ADD PLOTS OF HABITAT COVARIATES
  
  
  
  ## ------   2. GENERATE DETECTORS -----
  
  message("Preparing detectors characteristics... ")
  
  ## ------     2.1. GENERATE DETECTORS CHARACTERISTICS -----
  
  ##-- Generate raster of sub-detectors based on the study area
  subdetectors.r <- raster::disaggregate(
    x = habitat$habitat.rWthBuffer,
    fact = raster::res(habitat$habitat.r)[1]/detectors$resolution.sub)
  
  ##-- Generate NGS detectors based on the raster of sub-detectors
  detectors <- makeSearchGrid( 
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
  
  ##-- make a spatial grid from polygon
  detectors$grid <- sf::st_as_sf(raster::rasterToPolygons(
    x = raster::aggregate( x = subdetectors.r,
                           fact = detectors$resolution/detectors$resolution.sub),
    fun = function(x){x>0})) %>%
    mutate(id = 1:nrow(.))
  
  
  
  ## ------     2.2. GENERATE DETECTOR-LEVEL COVARIATES -----
  
  ## ------       2.2.1. EXTRACT COUNTIES -----
  
  ##-- Assign counties to detectors
  dist <- sf::st_distance(detectors$main.detector.sp, COUNTIES)
  detCounties1 <- apply(dist, 1, function(x) COUNTIES$county[which.min(x)])
  
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
  
  ##-- Load map of distance to roads (1km resolution)
  DistAllRoads <- raster::raster(file.path(data.dir,"Roads/MinDistAllRoads1km.tif"))
  
  ##-- Fasterize to remove values that fall in the sea
  r <- fasterize::fasterize(sf::st_as_sf(GLOBALMAP), DistAllRoads)
  r[!is.na(r)] <- DistAllRoads[!is.na(r)]
  DistAllRoads <- r
  
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
  
  
  
  ## ------       2.2.3. EXTRACT PRESENCE OF OTHER SAMPLES ------
  
  habitat.rWthBufferPol <- stars::st_as_stars(habitat$habitat.rWthBuffer) %>%
    sf::st_as_sf(., 
                 as_points = FALSE,
                 merge = TRUE) %>%
    dplyr::filter(Habitat %in% 1)

  r.detector <- raster::aggregate( 
    subdetectors.r,
    fact = (detectors$resolution/detectors$resolution.sub))
  
  
  
  ## ------         2.2.3.1. SKANDOBS ------
  
  ##-- Subset SkandObs 
  skandObs <- skandObs %>%
    dplyr::filter( 
      ##-- ...based on monitoring season
      month %in% unlist(sampling.months),
      ##-- ... based on space 
      !is.na(as.numeric(sf::st_intersects(., habitat.rWthBufferPol))))
  
  # ##-- Rasterize at the detector level
  # r.list <- lapply(DATA$years, function(y){
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
  
  
  
  ## ------         2.2.3.2. ROVBASE ------
  
  ##-- Get all samples collected (all species)
  rovbaseObs <- readMostRecent( path = data.dir,
                                extension = ".xls",
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
  
  
  
  ## ------         2.2.3.3. COMBINE ROVBASE & SKANDOBS ------
  
  ##-- Smooth binary map
  ##-- We tried adjust = 0.05, 0.037,0.02 and decided to go for 0.02 
  habOwin <- spatstat.geom::as.owin(as.vector(raster::extent(r.detector)))
  ds.list <- lapply( DATA$years, function(y){
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
  colnames(detOtherSamples) <- paste0("detOtherSamples.", years)
  detectors$detectors.df <- cbind.data.frame(detectors$detectors.df, detOtherSamples)
  
  
  
  ## ------       2.2.4. FORMAT detCovs ------
  
  detCovs <- array(NA, c(n.detectors, 2, n.years))
  for(t in 1:n.years){
    detCovs[ ,1,t] <- detectors$detectors.df[ ,"roads"]
  }
  detCovs[ ,2, ] <- detOtherSamples
  dimnames(detCovs) <- list( "detectors" = 1:n.detectors,
                             "covariates" = c("roads", "obs"),
                             "years" = years)
  detectors$covariates <- detCovs
  
  
  ##-- Merge with the habitat grid
  detectors$grid <- dplyr::left_join(
    x = detectors$grid,
    y = detectors$detectors.df,
    by = "id")
  
  
  
  ## ------       2.2.5. PLOTS -----
  
  ##-- Plot Carnivore observations maps
  pdf(file = file.path(working.dir, "figures", "CarnivoreObs.pdf"),
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
    plot(st_geometry(COUNTRIES[1,]), border = NA, col = "gray80")
    
    thisColumn <- which(names(st_drop_geometry(detectors$grid)) == paste0("detOtherSamples.",years[t]))
    thisCol <- ifelse(st_drop_geometry(detectors$grid)[,thisColumn] == 0, "white", "forestgreen")
    plot(detectors$grid[, paste0("detOtherSamples.",years[t])],
         border = NA, add = TRUE, legend = FALSE, col = thisCol)
    mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
    
    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
    }#if
  }#t
  dev.off()
  
  
  ##-- Plot distance to roads maps
  pdf(file = file.path(working.dir, "figures", "DistanceToRoads.pdf"),
      width = 18, height = 12)
  par(mar = c(0,0,0,0))
  plot(st_geometry(COUNTRIES[1, ]), border = NA, col = "gray80")
  plot(detectors$grid[, "roads"], add = T, border = NA,)
  segments(x0 = 830000, x1 = 830000,
           y0 = 6730000, y1 = 6730000 + 500000,
           col = grey(0.3), lwd = 4, lend = 2)
  text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
  dev.off()  
  
  
  
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
  
  
  
  ## ------   4. CREATE LOCAL OBJECTS -----
  
  # ##-- Get local habitat windows
  # habitat$localObjects <- getLocalObjects(
  #   habitatMask = habitat$habitat.mx,
  #   coords = habitat$scaledCoords,
  #   dmax = habitat$maxDist/habitat$resolution,
  #   plot.check = F)
  
  ##-- Get local detectors
  detectors$localObjects <- getLocalObjects(
    habitatMask = habitat$habitat.mx,
    coords = detectors$scaledCoords,
    dmax = detectors$maxDist/habitat$resolution,
    resizeFactor = resize.factor,
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
    ##-- Assign detector based on distance
    assignDetectors(
      data = .,                
      detectors = detectors$main.detector.sp,
      subDetectors = detectors$detector.sp,
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
      # ##-- Filter data for space
      # Country_sample == "(N)"
    ) %>% 
    ##-- Assign detector based on distance
    assignDetectors(
      data = .,
      detectors = detectors$main.detector.sp,
      radius = detectors$resolution)
  
  
  ## ------     6.3. PLOT NGS and DEAD RECOVERY MAPS ----- 
  
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
      plot( sf::st_geometry(data.dead[data.dead$Year == years[t] & 
                                        data.dead$Legal, ]), add = TRUE, col = "slateblue1", pch = 3),
      silent = TRUE)
    try(
      plot( sf::st_geometry(data.dead[data.dead$Year == years[t], ]), add = TRUE, col = "slateblue4", pch = 3),
      silent = TRUE)
    plot( sf::st_geometry(COUNTRIES), border = "gray20", col = NA, add = TRUE)
    
    ##-- Add year
    graphics::mtext(text = years[t],
                    side = 1, line = -18,
                    adj = 0.18, cex = 1.2)
  }#t
  dev.off()
  
  
  
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
    
    data.alive$data.sp <- data.alive$data.sp %>%
      dplyr::filter(Sex %in% thisSex)
    
    data.dead <- data.dead %>%
      dplyr::filter(Sex %in% thisSex)
    
    
    
    ## ------     7.2. GENERATE DETECTION HISTORY ARRAYS -----
    
    y.ar <- makeY( data = data.alive$data.sp,
                   detectors = detectors$main.detector.sp,
                   method = "Binomial",
                   data2 = data.dead,
                   detectors2 = detectors$main.detector.sp,
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
      
      ##-- Identify detections further then maxDist
      distances[[t]] <- checkDistanceDetections(
        y = y.ar.ALIVE[,,t], 
        detector.xy = detectors$detectors.df[ ,c("x","y")], 
        max.distance = detectors$maxDist,
        method = "pairwise",
        plot.check = F,
        verbose = F)
      
      ##-- Remove detections that are further than the threshold
      y.ar.ALIVE[ , ,t] <- y.ar.ALIVE[ , ,t] * (1-distances[[t]]$y.flagged)
    }#t
    
    
    
    ## ------     7.4. AUGMENT DETECTION HISTORIES -----
    
    y.alive <- makeAugmentation( 
      y = y.ar.ALIVE,
      aug.factor = DATA$aug.factor,
      replace.value = 0)
    
    y.dead.ar <- makeAugmentation( 
      y = y.ar.DEAD,
      aug.factor = DATA$aug.factor,
      replace.value = 0)
    
    
    
    ## ------     7.5. TRANSFORM Y TO SPARSE MATRICES -----
    
    y.sparse <- nimbleSCR::getSparseY(y.alive)
    
    
    
    ## ------ IV. MODEL SETTING ------- 
    
    ## -----    1. NIMBLE CODE ------
    
    modelCode <- nimbleCode({
      ##----- SPATIAL PROCESS -----
      tau ~ dunif(0,4)
      
      for(tp in 1:2){
        betaDens[1,tp] ~ dnorm(0,0.01)
        betaDens[2,tp] ~ dnorm(0,0.01)
        
        habIntensity[1:n.habwindows,tp] <- exp(betaDens[1,tp] * habCovs[1:n.habwindows,1] +
                                                 betaDens[2,tp] * habCovs[1:n.habwindows,2])
        sumHabIntensity[tp] <- sum(habIntensity[1:n.habwindows,tp])
        logHabIntensity[1:n.habwindows,tp] <- log(habIntensity[1:n.habwindows,tp])
        logSumHabIntensity[tp] <- log(sumHabIntensity[tp])
      }#tp
      
      
      for(i in 1:n.individuals){
        sxy[i,1:2,1] ~ dbernppAC(
          lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
          upperCoords = upperHabCoords[1:n.habwindows,1:2],
          logIntensities = logHabIntensity[1:n.habwindows,1],
          logSumIntensity = logSumHabIntensity[1],
          habitatGrid = habitatGrid[1:y.max,1:x.max],
          numGridRows = y.max,
          numGridCols = x.max)
        
        for(t in 2:n.years){
          sxy[i,1:2,t] ~ dbernppACmovement_normal(
            lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
            upperCoords = upperHabCoords[1:n.habwindows,1:2],
            s = sxy[i,1:2,t-1],
            sd = tau,
            baseIntensities = habIntensity[1:n.habwindows,2],
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            numGridRows = y.max,
            numGridCols = x.max,
            numWindows = n.habwindows)
        }#i
      }#t
      
      
      
      ##----- DEMOGRAPHIC PROCESS -----
      omeg1[1:2] ~ ddirch(alpha[1:2])   
      
      for(i in 1:n.individuals){ 
        z[i,1] ~ dcat(omeg1[1:2]) 
      }#i
      
      for(t in 1:n.years){
        n.available[t] <- sum(z[1:n.individuals,t] == 1)
        N[t] <- sum(z[1:n.individuals,t] == 2)
      }#t
      
      # for(t in 1:(n.years-1)){
      #   gamma[t] ~ dunif(0,1)
      #   mhW[t] ~ dunif(-10,10)
      #   mhH[t] ~ dunif(-10,10)
      #   
      #   for(i in 1:n.individuals){
      #     z[i,t+1] ~ dcatHR( z = z[i,t],
      #                        gamma = gamma[t],
      #                        mhH = mhH[t],
      #                        mhW = mhW[t])
      #   }#i
      # }#t
      
      for(t in 1:(n.years-1)){
        
        gamma[t] ~ dunif(0,1)
        w[t] ~ dunif(0,1)
        h[t] ~ dunif(0,1)
        r[t] ~ dunif(0,1)
        phi[t] <- 1-h[t]-w[t]
        
        ones.dead[t] ~ dbern(step(1 - (h[t] + w[t])))  
        
        omega[1,1:5,t] <- c(1-gamma[t], gamma[t], 0, 0, 0)
        omega[2,1:5,t] <- c(0,phi[t],h[t],w[t]*r[t],w[t]*(1-r[t]))
        omega[3,1:5,t] <- c(0,0,0,0,1)
        omega[4,1:5,t] <- c(0,0,0,0,1)
        omega[5,1:5,t] <- c(0,0,0,0,1)
        
        
        for(i in 1:n.individuals){ 
          z[i,t+1] ~ dcat(omega[z[i,t],1:5,t]) 
        }#i 								
      }#t 
      
      
      ##----- DETECTION PROCESS -----
      sigma ~ dunif(0,5)
      
      for(d in 1:n.detCovs){
        betaDet[d] ~ dunif(-5,5)
      }#d
      
      for(t in 1:n.years){
        for(c in 1:n.counties){
          p0[c,t] ~ dunif(0,1)
        }#c  
        
        for(i in 1:n.individuals){
          y.alive[i,1:lengthYCombined,t] ~ dbinomLocal_normalCovs(
            s = sxy[i,1:2,t],
            size = size[1:n.detectors],
            p0Traps = p0[1:n.counties,t], 
            sigma = sigma,
            trapCoords = detCoords[1:n.detectors,1:2],
            localTrapsIndices = localDetIndices[1:n.habwindows,1:localDetNumMax],
            localTrapsNum = localDetNum[1:n.habwindows],
            resizeFactor = 1,
            habitatGrid = habitatGrid[1:y.max,1:x.max],
            indicator = (z[i,t] == 2),
            lengthYCombined = lengthYCombined,
            allowNoLocal = 1,
            trapCovsIntercept = county[1:n.detectors],
            trapCovs = detCovs[1:n.detectors,1:n.detCovs,t],
            trapBetas = betaDet[1:n.detCovs])
        }#i
      }#t
      
      for(t in 1:n.years){
        for(i in 1:n.individuals){
          y.dead.legal[i,t] ~ dbern(z[i,t] == 3) 
          y.dead.other[i,t] ~ dbern(z[i,t] == 4) 
        }#i
      }#t
    })
    
    
    
    ## ------   2. NIMBLE CONSTANTS -----
    
    nimConstants <- list( n.individuals = dim(y.alive)[1], 
                          n.detectors = dim(detectors$scaledCoords)[1], 
                          n.detCovs = dim(detCovs)[2],
                          n.years = dim(y.alive)[3],
                          n.habwindows = nrow(habitat$habitat.df),
                          n.counties = max(detectors$detectors.df$counties),
                          y.max = nrow(detectors$localObjects$habitatGrid),
                          x.max = ncol(detectors$localObjects$habitatGrid),
                          lengthYCombined = y.sparse$lengthYCombined,
                          localDetNumMax = detectors$localObjects$numLocalIndicesMax,
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
                                 aug.factor = DATA$aug.factor,
                                 replace.value = NA)
    
    ##-- Compress back to yearly z
    zYears <- apply(zMonths, c(1,3), function(x){
      if(any(x[1:length(unlist(sampling.months))] == 1, na.rm = T)){
        2
      } else {
        if(any(x[1:length(unlist(sampling.months))] >= 2, na.rm = T)){
          5 }
        else {NA}}})
    
    z.data <- zYears
    allDead <- apply(z.data, 1, function(x)all(x==5))
    z.data[allDead, ] <- 1
    
    
    
    ## ------     3.2. GENERATE INITIAL z -----
    
    z.init <- t(apply(z.data, 1, function(zz){
      out <- zz
      out[] <- 1
      if(any(!is.na(zz))){
        ##-- Set all occasions after the last detection to 5
        if(sum(zz == 5, na.rm = T)<1){
          range.det <- range(which(!is.na(zz)))
          if(range.det[1]>1) zz[1:(range.det[1]-1)] <- 1
          if(range.det[2]<length(zz)) zz[(range.det[2]+1):length(zz)] <- 5
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
    
    
    
    ## ------     3.3. GENERATE y.dead -----
    
    legal.mx <- do.call(rbind, lapply(dimnames(y.alive)[[1]], function(x){
      out <- rep(0,dim(z.data)[2])
      if(x %in% data.dead$Id[data.dead$Legal]) out <- rep(1,dim(z.data)[2])
      return(out)
    }))
    
    y.dead.legal <- z.data
    y.dead.legal[] <- ifelse(z.data %in% c(5) & legal.mx == 1,1,0)
    y.dead.legal <- t(apply(y.dead.legal, 1, function(x){
      out <- x
      out[] <- 0
      if(any(x==1)) out[min(which(x==1))] <- 1
      return(out)
    }))
    
    other.mx <- do.call(rbind, lapply(dimnames(y.alive)[[1]], function(x){
      out <- rep(0,dim(z.data)[2])
      if(x %in% data.dead$Id[!data.dead$Legal]) out <- rep(1,dim(z.data)[2])
      return(out)
    }))
    
    y.dead.other <- z.data
    y.dead.other[] <- ifelse(z.data %in% c(5) & other.mx == 1,1,0)
    y.dead.other <- t(apply(y.dead.other, 1, function(x){
      out <- x
      out[] <- 0
      if(any(x==1)) out[min(which(x==1))] <- 1
      return(out)
    }))
    
    
    ##-- DISTINGUISH MORTALITY SOURCE IN z.data and z.init
    z.data[] <- ifelse(y.dead.legal == 1, 3, z.data)
    z.data[] <- ifelse(y.dead.other == 1, 4, z.data)
    
    
    
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
                          habitatGrid = detectors$localObjects$habitatGrid,
                          baseIntensities = rep(2,nimConstants$n.habwindows),
                          sd = 1,
                          radius = habitat$maxDist)
    s.init[!is.na(s.data)] <- NA
    
    ##-- Check AC-movement distances
    # test <- s.init
    # test[is.na(test)] <- s.data[is.na(test)]
    # dist <- sapply(2:dim(test)[3],
    #                function(t){
    #                  sqrt((test[ ,1,t]-test[ ,1,t-1])^2 + (test[ ,2,t]-test[ ,2,t-1])^2)
    #                })
    # hist(dist)
    # max(dist)
    
    ##-- Fix some annoying individuals
    # which(dist > habitat$maxDist/habitat$resolution, arr.ind = T)
    # dist[16, ]
    # y.sparse$yCombined[16,1, ]
    # s.init[16, ,3] <-  c(18.9,46.34)
    # tmp <- trunc(s.init[16, ,3]) + 1
    # detectors$localObjects$habitatGrid[tmp[2],tmp[1]]
    #
    # dist[222, ]
    # y.sparse$yCombined[222,1, ]
    # s.init[222, ,7] <- (s.init[222, ,8] + s.init[222, ,6])/2
    # tmp <- trunc(s.init[222, ,7]) + 1
    # detectors$localObjects$habitatGrid[tmp[2],tmp[1]]
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
                     y.dead.legal = y.dead.legal,
                     y.dead.other = y.dead.other,
                     lowerHabCoords = habitat$scaledLowerCoords, 
                     upperHabCoords = habitat$scaledUpperCoords, 
                     habCovs = cbind(habitat$habitat.df$dead.reco.trunc,
                                     habitat$habitat.df$skandObs.smooth),
                     habitatGrid = detectors$localObjects$habitatGrid,
                     alpha = c(1,1),
                     ones.dead = rep(1,dim(y.alive)[3]-1),
                     detCoords = detectors$scaledCoords,
                     size = detectors$detectors.df$size,
                     detCovs = detCovs,
                     localDetIndices = detectors$localObjects$localIndices,
                     localDetNum = detectors$localObjects$numLocalIndices)
    
    
    
    ## ------   5. NIMBLE PARAMETERS -----
    
    nimParams <- c("N", "n.available", "omeg1",
                   "gamma", "p0", "h", "w", "phi", "r",
                   "tau",
                   "p0", "sigma",
                   "betaDet", "betaDens")
    
    nimParams2 <- c("z", "sxy")
    
    
    
    ## ------   6. CHECK INPUT VALIDITY ------
    
    
    
    ## ------   7. SAVE INPUTS ----- 
    
    for(c in 1:4){
      nimInits <- list( "sxy" = s.init,
                        "z" = z.init,
                        "tau" = 0.4,
                        "betaDens" = matrix(stats::runif(4,0,1),nrow = 2),
                        "omeg1" = c(0.7,0.3),
                        "gamma" = stats::runif(dim(y.alive)[3]-1,0.02,0.1),
                        "h" = stats::runif(dim(y.alive)[3]-1,0.05,0.1),
                        "w" = stats::runif(dim(y.alive)[3]-1,0.1,0.4),
                        "r" = stats::runif(dim(y.alive)[3]-1,0,0.5),
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

