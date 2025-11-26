#' @title RovQuant OPSCR wolverine output processing
#' 
#' @description
#' \code{processRovquantOutput_wolverine} calls a custom Rmarkdown template that combines 
#' and processes MCMC outputs from NIMBLE models and produces figures,
#' tables and rasters of interest (e.g. population density maps)
#' 
#' @param data.dir A \code{path}
#' @param working.dir A \code{path}
#' @param nburnin An \code{integer} denoting the number of iterations to be removed from each MCMC as burnin.
#' @param niter An \code{integer} denoting the number of MCMC iterations to be used for density extraction.
#' @param extraction.res A \code{integer} denoting the raster resolution for density extraction.
#' 
#' @return 
#' A \code{.RData} file with the clean NGS and dead recovery data objects
#' for the species and period specified.
#' A \code{html} report summarizing the data cleaning process
#' Additional \code{.png} images that can be reused somewhere else.
#'
#' @author Pierre Dupont
#' 
#' @import sf 
#' @import raster
#' @import dplyr
#' @importFrom fasterize fasterize
#' @importFrom adehabitatHR estUDm2spixdf kernelUD
#' @importFrom stats density
#' @importFrom grDevices adjustcolor dev.off pdf png grey
#' @importFrom graphics axis abline par
#' @importFrom stars st_as_stars
#' @importFrom nimbleSCR scaleCoordsToHabitatGrid
#' @importFrom abind abind
#' @importFrom utils data
#' @importFrom xtable xtable
#' 
#' @rdname processRovquantOutput_wolverine
#' @export
processRovquantOutput_wolverine <- function(
  ##-- paths
  data.dir = getwd(),
  working.dir = NULL,
  ##-- MCMC
  nburnin = 0,
  niter = 100,
  ##-- Density 
  extraction.res = 5000,
  ##-- Miscellanious
  overwrite = FALSE
){

  ## ------ 0. BASIC SET-UP ------
  
  if(is.null(working.dir)){working.dir <- getwd()}
  
  ##-- Extract date from the last cleaned data file
  DATE <- getMostRecent( 
    path = file.path(working.dir, "data"),
    pattern = "CleanData_wolverine")
  
  ##-- Initialize output list
  out <- list( SPECIES = "Wolverine",
               engSpecies = "wolverine",
               DATE = DATE)
  
  ##-- States alive
  alive.states <- 2
  
  
  
  ## ------ 1. LOAD NECESSARY INPUTS -----
 
  ##-- Females
  load(list.files(file.path(working.dir, "nimbleInFiles/female"), full.names = T)[1])
  nimDataF <- nimData
  nimInitsF <- nimInits
  
  ##-- Males
  load(list.files(file.path(working.dir, "nimbleInFiles/male"), full.names = T)[1])
  nimDataM <- nimData
  nimInitsM <- nimInits
  
  ##-- Remove unnecessary objects from memory
  rm(list = c("nimInits", "nimData"))
  gc(verbose = FALSE)
  
  ##-- Habitat
  load(file.path( working.dir, "data",
                  paste0("Habitat_wolverine_", DATE, ".RData")))
  
  ##-- Detectors
  load(file.path( working.dir, "data",
                  paste0("Detectors_wolverine_", DATE, ".RData")))
  
  ##-- Load filtered data
  load(file.path( working.dir, "data",
                  paste0("FilteredData_wolverine_", DATE, ".RData")))
  
  ##-- Habitat Rasters
  if(extraction.res <= 1000) {
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
  
  ##-- Extract years
  years <- as.numeric(dimnames(nimDataF$z)[[2]])
  n.years <- length(years) 
  
  ##-- years not sampled in Norrbotten
  yearsSampledNorrb <- c(2016:2018,2023)
  yearsNotSampled <- which(!years %in% yearsSampledNorrb)
  
  ##-- Polygons of Sweden & Norway
  COUNTRIES <- REGIONS %>%
    dplyr::filter(country %in% c("SWE","NOR")) %>%
    dplyr::group_by(country) %>%
    dplyr::summarize()
  
  ##-- Polygons of counties in Sweden & Norway
  COUNTIES <- REGIONS %>%
    group_by(county) %>%
    summarize()
  
  ##-- Merge counties for practical reasons
  COUNTIES_AGGREGATED <- REGIONS %>%
    dplyr::mutate(id = case_when(
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
    dplyr::group_by(id) %>%
    dplyr::summarize() %>%
    sf::st_simplify( ., preserveTopology = T, dTolerance = 500)

  ##-- Prepare raster of countries
  countryRaster <- habitatRasterResolution$`5km`[["Countries"]]
  
  
  ## ------ 2. PROCESS MCMC SAMPLES -----
  
  message("## Processing model MCMC outputs...")
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  mcmcTest <- TRUE
  if(!overwrite){
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
    ## ------   2.1. FEMALES -----
    
    ##-- Compile MCMC bites
    gc(verbose = FALSE)
    nimOutput_F <- collectMCMCbites( path = file.path(working.dir, "nimbleOutFiles/female"),
                                     burnin = nburnin)
    
    ##-- Traceplots
    gc(verbose = FALSE)
    grDevices::pdf(file.path(working.dir, "figures/traceplots_F.pdf"))
    plot(nimOutput_F$samples[ ,!is.na(nimOutput_F$samples[[1]][1, ])])
    grDevices::dev.off()

    ##-- Process MCMC output
    gc(verbose = FALSE)
    results_F <- processCodaOutput( nimOutput_F$samples,
                                    params.omit = c("sxy","z"))
    
    gc(verbose = FALSE)
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
    
    ##-- Rescale sigma & dmean to the original coordinate system
    results_F$sims.list$sigma <- results_F$sims.list$sigma * raster::res(habitat$habitat.r)[1]
    results_F$sims.list$dmean <- results_F$sims.list$dmean * raster::res(habitat$habitat.r)[1]
    
    
    
    ## ------   2.2. MALES -----
    
    ##-- Compile MCMC bites
    gc(verbose = FALSE)
    nimOutput_M <- collectMCMCbites( path = file.path(working.dir, "nimbleOutFiles/male"),
                                     burnin = nburnin)
    
    ##-- Traceplots
    gc(verbose = FALSE)
    grDevices::pdf(file.path(working.dir, "figures/traceplots_M.pdf"))
    plot(nimOutput_M$samples[ ,!is.na(nimOutput_M$samples[[1]][1, ])])
    dev.off()
    
    ##-- Process MCMC output
    gc(verbose = FALSE)
    results_M <- processCodaOutput( nimOutput_M$samples,
                                    params.omit = c("sxy","z"))
    gc(verbose = FALSE)
    resultsSXYZ_M <- processCodaOutput(nimOutput_M$samples2)
    
    ##-- Remove unnecessary objects from memory
    rm(list = c("nimOutput_M"))
    gc(verbose = FALSE)
    
    ##-- Rescale sxy to the original coordinate system
    dimnames(resultsSXYZ_M$sims.list$sxy)[[3]] <- c("x","y")
    resultsSXYZ_M$sims.list$sxy <- nimbleSCR::scaleCoordsToHabitatGrid(
      coordsData = resultsSXYZ_M$sims.list$sxy,
      coordsHabitatGridCenter = habitat$habitat.df,
      scaleToGrid = FALSE)$coordsDataScaled
    
    ##-- Rescale sigma & dmean to the original coordinate system
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
    minIterSigma <- min(dim(results_F$sims.list$sigma)[1],dim(results_M$sims.list$sigma)[1])
    iterSigma <- seq(1, minIterSigma, by = minIterSigma/minIter)
    resultsSXYZ_MF$sims.list$sigma <- abind::abind(results_M$sims.list$sigma[iterSigma, ],
                                                   results_F$sims.list$sigma[iterSigma, ],
                                                   along = 3)
    dimnames(resultsSXYZ_MF$sims.list$sigma)[[3]] <- c("M","F")
    
    ##-- sex
    resultsSXYZ_MF$sims.list$sex <- rep(c("M","F"),
                                        c(dim(resultsSXYZ_M$sims.list$sxy)[2],
                                          dim(resultsSXYZ_F$sims.list$sxy)[2]))
    
    ##-- SAVE & LOAD DATA
    save( results_F, results_M, resultsSXYZ_MF,
          file = file.path( working.dir, "data",
                            paste0("MCMC_wolverine_", DATE, ".RData")))
  }

  ##-- Number of activity center posterior samples
  n.mcmc <- dim(resultsSXYZ_MF$sims.list$z)[1]
  gc(verbose = FALSE)
  
  
  
  ## ------ 3. EXTRACT DENSITY -----
  
  message("## Processing density outputs...")
  
  ##-- Reduce the size of the output to 'niter' MCMC samples
  if(n.mcmc >= niter) {
    iter <- round(seq(1, n.mcmc, length.out = niter))
  } else {
    message(paste0( "The number of MCMC samples available (", n.mcmc,
                    ") is less than niter = ", niter,
                    ".\nusing niter = ", n.mcmc, " instead."))
    iter <- 1:n.mcmc
  }
  
  ##-- Remove buffer from the habitat
  habitat.rWthBuffer <- habitat$habitat.rWthBuffer
  habitat.rWthBuffer[habitat.rWthBuffer[] %in% 0] <- NA
  searchedPolygon <- sf::st_as_sf(stars::st_as_stars(habitat.rWthBuffer), 
                                  as_points = FALSE, merge = TRUE)
  searchedPolygon <- searchedPolygon[searchedPolygon$Habitat > 0, ]
  
  ##-- Habitat raster with extent used in the model
  habitatPolygon5km <- raster::crop(extraction.raster$Habitat, habitat$habitat.r)
  
  ##-- Create raster of countries for extraction
  rrCountries <- extraction.raster$Countries
  rrCountries[extraction.raster$Countries[] %in% c(1,3)] <- NA
  areaCountriesTotal <- table(raster::factorValues(rrCountries, rrCountries[]))*raster::res(rrCountries)[1]*1e-6
  rrCountries <- raster::mask(rrCountries, searchedPolygon)
  rrCountries <- raster::crop(rrCountries, habitat$habitat.r)
  #plot(rrCountries)
  ##-- Calculate studied area of each county
  areaCountries <- table(raster::factorValues(rrCountries, rrCountries[]))*raster::res(rrCountries)[1]*1e-6 
  percCountries <- round(areaCountries/areaCountriesTotal, 2)
  percTotal <- round(sum(areaCountries)/sum(areaCountriesTotal),2)
  
  ##-- Create raster of counties for extraction
  levels(extraction.raster$Counties)[[1]][c(4,5,6,10,12,13,14,15,17,18,19,20),2] <- c(
    "Södermanland", "Östergötland","Jönköping", "Skåne", "Västra Götaland",
    "Värmland","Örebro","Västmanland","Gävleborg",
    "Västernorrland","Jämtland","Västerbotten")
  rrCounties <- extraction.raster$Counties
  rrCounties[extraction.raster$Countries[] %in% c(1,3)] <- NA
  areaCountiesTotal <- table(raster::factorValues(rrCounties, rrCounties[]))*res(rrCounties)[1]*1e-6
  rrCounties <- raster::mask(rrCounties, searchedPolygon)
  rrCounties <- raster::crop(rrCounties, habitat$habitat.r)
  #plot(rrCounties)
  ##-- Calculate studied area of each county
  areaCounties <- table(raster::factorValues(rrCounties, rrCounties[]))*res(rrCounties)[1]*1e-6
  areaCountiesTotal <- areaCountiesTotal[names(areaCountiesTotal) %in% names(areaCounties)]
  percCounties <- round(areaCounties/areaCountiesTotal, 2)
  
  ##-- Create raster of carnivore regions for extraction
  rrRegions <- extraction.raster$Regions
  rrRegions[extraction.raster$Countries %in% c(1,3)] <- NA
  rrRegions[rrRegions[ ] %in% c(18,19,20,21)] <- 1
  rrRegions[rrRegions[ ] %in% c(13,17,16,14,15,12,22,3)] <- 2
  rrRegions[rrRegions[ ] %in% c(4,5,10,6,7,9,11,8)] <- 3
  rrRegions <- ratify(rrRegions)
  levels(rrRegions)[[1]] <- data.frame(
    "ID" = c(1,2,3,23:30),
    "Regions"= c( "Nordre","Midtre","Söndre",
                  "Region 3","Region 1","Region 2","Region 4",
                  "Region 7","Region 6","Region 8","Region 5"))
  areaRegionsTotal <- table(factorValues(rrRegions, rrRegions[]))*res(rrRegions)[1]*1e-6
  rrRegions <- mask(rrRegions, searchedPolygon)
  rrRegions <- crop(rrRegions, habitat$habitat.r)
  #plot(rrRegions)
  ##-- Calculate studied area of each county
  areaRegions <- table(factorValues(rrRegions,rrRegions[]))*res(rrRegions)[1]*1e-6
  areaRegionsTotal <- areaRegionsTotal[names(areaRegionsTotal) %in% names(areaRegions)]
  percRegions <- round(areaRegions/areaRegionsTotal, 2)
  
  ##-- Merge the percentages
  percAllRegions <- c(percTotal, percCountries, percRegions, percCounties)
  names(percAllRegions)[1] <- "Total"
  
  
  ##-- Calculate density only if necessary
  ##-- Check that a file with that name does not already exist to avoid overwriting
  densTest <- TRUE
  if(!overwrite){
    fileName <- paste0("Density_wolverine_", DATE, ".RData")
    if (file.exists(file.path(working.dir, "data", fileName))) {
      message(paste0("A density output file named '", fileName, "' already exists in: \n",
                     file.path(working.dir, "data")))
      message("Do you want to proceed and overwrite the existing density output file? (y/n) ")
      question1 <- readLines(n = 1)
      if (regexpr(question1, 'y', ignore.case = TRUE) != 1) {
        message("Not overwriting existing files...")
        message(paste0("Loading '", fileName, "' instead..."))
        load(file.path(working.dir, "data", fileName))
        densTest <- FALSE
      } else {
        message(paste0("Now overwriting '", fileName,"'.\n"))
      }
    }
  } 

  
  if(densTest){
    
    message("## Extracting population density... \n## This might take a while...")
    
    
    ## ------   1. PREPARE DENSITY EXTRACTION ------
    
    ##-- Get the objects to run the density function
    ##-- COUNTRY
    densityInputCountries <- suppressWarnings(getDensityInput(
      regions = rrCountries, 
      habitat = habitatPolygon5km,
      s = resultsSXYZ_MF$sims.list$sxy,
      plot.check = FALSE))
   # rownames(densityInputCountries$regions.rgmx)
    
    ##-- COUNTIES
    densityInputCounties <- suppressWarnings(getDensityInput( 
      regions = rrCounties, 
      habitat = habitatPolygon5km,
      s = resultsSXYZ_MF$sims.list$sxy,
      plot.check = FALSE))
    # rownames(densityInputCounties$regions.rgmx)
    
    ##-- REGIONS
    densityInputRegions <- suppressWarnings(getDensityInput( 
      regions = rrRegions, 
      habitat = habitatPolygon5km,
      s = resultsSXYZ_MF$sims.list$sxy,
      plot.check = FALSE))
   # rownames(densityInputRegions$regions.rgmx)
    
    ##-- Merge country, county & region matrices to allow simultaneous estimation
    regionID <- rbind( densityInputCountries$regions.rgmx,
                       densityInputRegions$regions.rgmx,
                       densityInputCounties$regions.rgmx)
    row.names(regionID) <- c(row.names(densityInputCountries$regions.rgmx),
                             row.names(densityInputRegions$regions.rgmx),
                             row.names(densityInputCounties$regions.rgmx))
    
    ##-- Free up space
    rm(list = c("densityInputCounties","densityInputCountries"))
    sx_extract <- densityInputRegions$sx
    sy_extract <- densityInputRegions$sy
    z_extract <- resultsSXYZ_MF$sims.list$z
    habitat.id_extract <- densityInputRegions$habitat.id
    habitat.xy_extract <- densityInputRegions$habitat.xy
    inputRaster <- densityInputRegions$regions.r
    rm(list = c("densityInputRegions"))
    
    
    
    ## ------   2. AC-BASED DENSITY ------
    ## ------     2.1. MALE & FEMALES ------
    
    ##-- EXTRACT DENSITY 
    ACdensity <- list()
    for(t in 1:n.years){
      ACdensity[[t]] <- GetDensity(
        sx = as.matrix(sx_extract[ , ,t]),
        sy = as.matrix(sy_extract[ , ,t]),
        z = as.matrix(z_extract[ , ,t]),
        IDmx = habitat.id_extract,
        aliveStates = alive.states,
        regionID = regionID,
        returnPosteriorCells = F)
    }#t
    names(ACdensity) <- years+1
    
    
    
    ## ------     2.2. MALE -----
    
    IDMales <- which(resultsSXYZ_MF$sims.list$sex == "M")
    
    ACdensityM <- list()
    for(t in 1:n.years){
      ACdensityM[[t]] <- GetDensity(
        sx = sx_extract[ ,IDMales,t],
        sy = sy_extract[ ,IDMales,t],
        z = z_extract[ ,IDMales,t],
        IDmx = habitat.id_extract,
        aliveStates = alive.states,
        regionID = regionID,
        returnPosteriorCells = F)
    }#t
    names(ACdensityM) <- years+1
    
    
    
    ## ------     2.3. FEMALE -----
    
    IDFemales <- which(resultsSXYZ_MF$sims.list$sex == "F")
    
    ACdensityF <- list()
    for(t in 1:n.years){
      ACdensityF[[t]] <- GetDensity(
        sx = sx_extract[ ,IDFemales,t],
        sy = sy_extract[ ,IDFemales,t],
        z = z_extract[ ,IDFemales,t],
        IDmx = habitat.id_extract,
        aliveStates = alive.states,
        regionID = regionID,
        returnPosteriorCells = F)
    }
    names(ACdensityF) <- years+1
    

    
    ## ------   3. UD-BASED DENSITY ------
    
    ##-- Combine male and female sigma
    sigma <- array(NA, c( dim(sx_extract)[1],
                          length(resultsSXYZ_MF$sims.list$sex),
                          dim(sx_extract)[3]))
    for(i in 1:length(resultsSXYZ_MF$sims.list$sex)){
      if(resultsSXYZ_MF$sims.list$sex[i] == "M"){
        sigma[ ,i, ] <- resultsSXYZ_MF$sims.list$sigma[ , ,"M"]
      } else {
        sigma[ ,i, ] <- resultsSXYZ_MF$sims.list$sigma[ , ,"F"]
      }
    }#i

    ##-- Rescale sigma to the raster resolution
    sigma <- sigma/raster::res(rrRegions)[1]
    
    UDdensity <- list()
    for(t in 1:n.years){
      UDdensity[[t]] <- rovquantR::GetSpaceUse(
        sx = sx_extract[iter, ,t],
        sy = sy_extract[iter, ,t],
        z = z_extract[iter, ,t],
        sigma = sigma[iter,,t],
        habitatxy = habitat.xy_extract,
        aliveStates = alive.states,
        regionID = regionID[c("Norway","Sweden"), ], 
        display_progress = TRUE,
        returnPosteriorCells = FALSE)
      
      ##-- Free up space
      UDdensity[[t]]$MedianCell <- NULL
      UDdensity[[t]]$CVCell <- NULL
      UDdensity[[t]]$CILCell <- NULL
      UDdensity[[t]]$CIHCell <- NULL
    }#t
    names(UDdensity) <- years+1
    
    
    
    ## ------   4. SAVE DENSITY OBJECTS ------
    
    save( inputRaster,
          ACdensity,
          ACdensityF,
          ACdensityM,
          UDdensity,
          file = file.path( working.dir, "data",
                            paste0("Density_wolverine_", DATE, ".RData")))
  } 
  
  
  
  ## ------ 4. FIGURES -----
  
  ##-- Plot parameters
  diffSex <- 0.2
  colCountries <- c("firebrick2", "deepskyblue2", "black")
  names(colCountries) <- c("Norway","Sweden", "Total")
  colCause  <- adjustcolor( c("#E69F00","#009E73"), 0.5)
  
  
  
  ## ------   4.1. DENSITY MAPS -----
  
  message("## Plotting population density maps...") 
  
  ##-- Create 5km raster for plotting
  rrNorway <- extraction.raster[["Countries"]]
  rrNorway[!rrNorway[] %in% c(2,4)] <- NA
  rrNorway[rrNorway[] %in% c(2,4)] <- 1
  rrNorway <- raster::crop(rrNorway, habitat$habitat.r)
  rrCombined <- rrRegions + rrNorway
  
  ##-- AC-density maps
  plotDensityMaps(
    input = inputRaster,
    estimates = ACdensity,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES,
    type = c("time.series", "last.year"),
    path = working.dir,
    name = "AC_Density")
  
  ##-- UD-density maps
  plotDensityMaps( 
    input = inputRaster,
    estimates = UDdensity,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES,
    type = c("time.series", "last.year"),#,"summary","summary_NOR"),
    species = "wolverine",
    labels = list("nor" = ACdensity[[n.years]]$summary["Norway",c("95%CILow","95%CIHigh")],
                  "swe" = ACdensity[[n.years]]$summary["Sweden",c("95%CILow","95%CIHigh")],
                  "both" = ACdensity[[n.years]]$summary["Total",c("95%CILow","95%CIHigh")]),
    x.labels = c(0.3,0.75,0.7),
    y.labels = c(0.8,0.7,0.05),
    path = working.dir,
    name = "UD_Density")
  dev.off()
  
  
  
  ## ------   4.2. ABUNDANCE TIME SERIES ------
  
  message("## Plotting abundance...")
  
  ##-- Plot N  
  # pdf(file = file.path(working.dir, "figures/Abundance_TimeSeries.pdf"),
  #     width = 12, height = 8.5)
  grDevices::png(filename = file.path(working.dir,"figures/Abundance_TimeSeries.png"),
                 width = 12, height = 8.5, units = "in", pointsize = 12,
                 res = 300, bg = NA)
  
  graphics::par(mar = c(5,8,3,1),
                las = 1,
                cex.lab = 2,
                cex.axis = 1.3,
                mgp = c(6, 2, 0),
                xaxs = "i",
                yaxs = "i")
  
  ymax <- 100*(trunc(max(unlist(lapply(ACdensity, function(x)max(colSums(x$PosteriorAllRegions)))))/100)+1)
  
  plot(-1000,
       xlim = c(0.5, n.years+0.5),
       ylim = c(0,ymax),
       xlab = "", ylab = paste("Estimated number of wolverines"),
       xaxt = "n", axes = F, cex.lab = 1.6)
  graphics::axis(1, at = c(1:(n.years)), labels = years+1, cex.axis = 1.5, padj = -1)
  graphics::axis(2, at = seq(0,ymax,200), labels = seq(0,ymax,200), cex.axis = 1.5, hadj = 0.5)
  graphics::abline(v = (1:n.years)+0.5, lty = 2)
  graphics::abline(h = seq(0,ymax, by = 100), lty = 2, col = "gray90")
  
  for(t in 1:n.years){
    ##-- Norway
    plotQuantiles(x = ACdensity[[t]]$PosteriorRegions["Norway", ],
                  at = t + diffSex,
                  width = 0.15,
                  col = colCountries[1])
    
    ##-- Sweden 
    add.star <- t %in% yearsNotSampled
    plotQuantiles(x = ACdensity[[t]]$PosteriorRegions["Sweden", ],
                  at = t - diffSex,
                  width = 0.15,
                  col = colCountries[2],
                  add.star = add.star)
    
    ##-- TOTAL
    plotQuantiles(x = colSums(ACdensity[[t]]$PosteriorAllRegions),
                  at = t,
                  width = 0.15,
                  col = colCountries[3],
                  add.star = add.star)
    
    # ##-- ADD NUMBER OF INDIVIDUALS DETECTED
    # xx <- c(t-0.25,t+0.25,t+0.25,t-0.25)
    # yy <- c(n.detected[t]-1,n.detected[t]-1,n.detected[t]+1,n.detected[t]+1)
    # polygon(xx, yy, border = NA, col = "goldenrod1")
  }#t
  box()
  
  ##-- legend
  par(xpd = TRUE)
  xx <- c(0.11*n.years, 0.24*n.years, 0.37*n.years) 
  yy <- c(200,200,200)
  labs <- c("Norway", "Sweden", "Total")
  polygon(x = c(0.08*n.years,0.46*n.years,0.46*n.years,0.08*n.years),
          y = c(150,150,250,250),
          col = adjustcolor("white", alpha.f = 0.9),
          border = "gray90")
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col = adjustcolor(colCountries,0.3))
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col = adjustcolor(colCountries,0.7))
  text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.2, pos = 4)
  
  dev.off()
  
  ##-- Remove unnecessary objects from memory
  gc(verbose = FALSE)
  
  
  
  
  ## ------   4.3. ABUNDANCE TIME SERIES BY SEX ------
  
  grDevices::png(filename = file.path(working.dir,"figures/Abundance_TimeSeries_bySex.png"),
                 width = 18, height = 8, units = "in", pointsize = 12,
                 res = 300, bg = NA)
  
  graphics::par( mfrow = c(1,2),
                 mar = c(5,8,3,1),
                 las = 1,
                 cex.lab = 2,
                 cex.axis = 1.3,
                 mgp = c(6, 2, 0),
                 xaxs = "i",
                 yaxs = "i")
  
  ymax <- 800
  
  ##-- FEMALES 
  plot(-1000,
       xlim = c(0.5, n.years+0.5), ylim = c(0,ymax),
       xlab = "", ylab = "Estimated number of females",
       xaxt = "n", axes = F, cex.lab = 1.6)
  graphics::axis(1, at = c(1:(n.years)), labels = years+1, cex.axis = 1.5, padj = -1)
  graphics::axis(2, at = seq(0,ymax,200), labels = seq(0,ymax,200), cex.axis = 1.5, hadj = 0.5)
  graphics::abline(v = (1:n.years)+0.5, lty = 2)
  graphics::abline(h = seq(0,ymax, by = 100), lty = 2, col = "gray90")
  
  for(t in 1:n.years){
    ##-- Norway
    plotQuantiles(x = ACdensityF[[t]]$PosteriorRegions["Norway", ],
                  at = t + diffSex,
                  width = 0.15,
                  col = colCountries[1])
    
    ##-- Sweden 
    add.star <- t %in% yearsNotSampled
    plotQuantiles(x = ACdensityF[[t]]$PosteriorRegions["Sweden", ],
                  at = t - diffSex,
                  width = 0.15,
                  col = colCountries[2],
                  add.star = add.star)
    
    ##-- TOTAL
    plotQuantiles(x = colSums(ACdensityF[[t]]$PosteriorAllRegions),
                  at = t,
                  width = 0.15,
                  col = colCountries[3],
                  add.star = add.star)
  }#t
  box()
  
  
  ##-- LEGEND
  par(xpd = TRUE)
  xx <- c(0.11*n.years, 0.3*n.years, 0.48*n.years) 
  yy <- c(100,100,100)
  labs <- c("Norway", "Sweden", "Total")
  polygon(x = c(0.08*n.years,0.6*n.years,0.6*n.years,0.08*n.years),
          y = c(50,50,150,150),
          col = adjustcolor("white", alpha.f = 0.9),
          border = "gray90")
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col = adjustcolor(colCountries,0.3))
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col = adjustcolor(colCountries,0.7))
  text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.2, pos = 4)
  
  
  ##-- MALES 
  plot(-1000,
       xlim = c(0.5, n.years+0.5), ylim = c(0,ymax),
       xlab = "", ylab = "Estimated number of males",
       xaxt = "n", axes = F, cex.lab = 1.6)
  graphics::axis(1, at = c(1:(n.years)), labels = years+1, cex.axis = 1.5, padj = -1)
  graphics::axis(2, at = seq(0,ymax,200), labels = seq(0,ymax,200), cex.axis = 1.5, hadj = 0.5)
  graphics::abline(v = (1:n.years)+0.5, lty = 2)
  graphics::abline(h = seq(0,ymax, by = 100), lty = 2, col = "gray90")
  
  for(t in 1:n.years){
    ##-- Norway
    plotQuantiles(x = ACdensityM[[t]]$PosteriorRegions["Norway", ],
                  at = t + diffSex,
                  width = 0.15,
                  col = colCountries[1])
    
    ##-- Sweden 
    add.star <- t %in% yearsNotSampled
    plotQuantiles(x = ACdensityM[[t]]$PosteriorRegions["Sweden", ],
                  at = t - diffSex,
                  width = 0.15,
                  col = colCountries[2],
                  add.star = add.star)
    
    ##-- TOTAL
    plotQuantiles(x = colSums(ACdensityM[[t]]$PosteriorAllRegions),
                  at = t,
                  width = 0.15,
                  col = colCountries[3],
                  add.star = add.star)
  }#t
  box()
  dev.off()
  
  
  
  # ## ------   4.4. VITAL RATES ------
  # 
  # message("## Plotting vital rates...")
  # 
  # ##-- Calculate mortality from estimated hazard rates (mhH and mhW) if necessary
  # if("mhH" %in% names(results_F$sims.list)){
  #   mhH1 <- exp(results_F$sims.list$mhH)
  #   mhW1 <- exp(results_F$sims.list$mhW)
  #   results_F$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
  #   results_F$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
  #   results_F$sims.list$phi <- 1 - results_F$sims.list$h - results_F$sims.list$w
  # }
  # 
  # if("mhH" %in% names(results_M$sims.list)){
  #   mhH1 <- exp(results_M$sims.list$mhH)
  #   mhW1 <- exp(results_M$sims.list$mhW)
  #   results_M$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
  #   results_M$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
  #   results_M$sims.list$phi <- 1 - results_M$sims.list$h - results_M$sims.list$w
  # }
  # 
  # # if(!"rho" %in% names(results_M$sims.list)){
  # #   results_M$sims.list$rho <- array(NA, c())
  # #   for(t in 1:(n.years-1)){
  # #     ##-- Number of recruits at t ==> ids with state 1 at t and 2 at t+1
  # #     n.recruit <- apply(resultsSXYZ_MF$sims.list$z[ , ,c(t,t+1)], 1, function(x) sum(x[ ,1]%in%1 & x[ ,2]%in%2))
  # #     ##-- Number of reproducing ids at t ==> ids with state 2 at t
  # #     alivetminus1 <- apply(resultsSXYZ_MF$sims.list$z[ , ,t], 1, function(x)sum(x %in% 2))
  # #     ##-- Store per-capita recruitment rate in the results
  # #   }#t
  # # }
  # 
  # 
  # 
  # ## ------     4.4.1. SURVIVAL ------
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("Survival.pdf")),
  # #     width = 10, height = 6)
  # grDevices::png(filename = file.path(working.dir,"figures/Survival.png"),
  #                width = 10, height = 6, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
  #              widths = c(0.05,1,0.30),
  #              heights = c(0.15,1))
  # 
  # par(mar = c(10,4.5,0.5,0.5), tck = 0, xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
  # plot(10, xlim = c(0.5, n.years-1 + 0.5), ylim = c(0,1),
  #      type ="n", xaxt = "n", xlab = "", ylab = "Survival")
  # axis(1, c(1:n.years),
  #      labels = paste(years, years+1, sep = "\n to \n"),
  #      padj = 0.6, cex = 0.9)
  # mtext("Years",side = 1,line = 5)
  # abline(v = 1:(n.years-1) + 0.5, lty = 2)
  # axis(2, tck = -0.02)
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- FEMALES
  #   plotQuantiles( x = results_F$sims.list$phi[ ,t],
  #                  at = t - diffSex,
  #                  col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles( x = results_M$sims.list$phi[ ,t],
  #                  at = t + diffSex,
  #                  col = colSex[2])
  # }#t
  # 
  # ##-- LEGEND
  # par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  # plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
  # points(c(4,4), c(4,3), pch = 15, cex = 5.5, col = colSex)
  # points(c(4,4), c(4,3), pch = 15, cex = 3, col = colSex)
  # text(c(5.3,5.3), c(4,3), c("Females", "Males"), cex = 2, pos = 4)
  # dev.off()
  # 
  # 
  # 
  # ## ------     4.4.2. MORTALITY ------
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("Mortality.pdf")),
  # #     width = 10, height = 9)
  # grDevices::png(filename = file.path(working.dir, "figures/Mortality.png"),
  #                width = 10, height = 9, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # nf <- layout(rbind(c(3,7,6),
  #                    c(4,1,6),
  #                    c(5,2,6)),
  #              widths = c(0.15,1,0.35),
  #              heights = c(0.15,1,1))
  # 
  # for(s in 1:2){
  #   if(s == 1){results <- results_F$sims.list} else {results <- results_M$sims.list}
  #   par(mar=c(8,4.5,0.5,1),tck=0,xaxs="i",cex.axis=1.3,cex.lab=1.6)
  #   plot(10, xlim = c(0.5, n.years-1+0.5), ylim = c(0,0.5),
  #        type = "n", xaxt = "n", xlab = "", ylab = "Mortality")
  #   axis(1, c(1:n.years),
  #        labels = paste(years, years+1, sep = "\n to \n"),
  #        padj = 0.6, cex = 0.9)
  #   mtext("Years",side = 1,line = 5)
  #   axis(2, tck = -0.02)
  #   abline(v = 1:(n.years-1)+0.5, lty = 2)
  #   myDev <- c(-0.15,+0.15)
  # 
  #   for(t in 1:(n.years-1)){
  #     plotQuantiles(x = results$h[ ,t],
  #                   at = t - diffSex,
  #                   col = colCause[1])
  #     plotQuantiles(x = results$w[ ,t],
  #                   at = t + diffSex,
  #                   col = colCause[2])
  #   }#t
  # }#s
  # 
  # ##-- LABELS
  # par(mar = c(0,0,0,0))
  # plot(1, axes = FALSE, ylim = c(-1,1), xlim = c(-1,1), type = "n")
  # my.labels=c("Female","Male")
  # for(i in my.labels){
  #   par(mar=c(0.5,0.5,0.5,0.5))
  #   plot(1,axes=FALSE,ylim=c(-1,1),xlim=c(-1,1),type="n")
  #   text(0,0.25,labels=i,srt=90,cex=3,font=2)
  # }
  # 
  # ##-- LEGEND
  # par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  # plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = F)
  # points(c(2,2),c(3.5,2.7),pch=15,cex=3.5,col = colCause)
  # points(c(2,2),c(3.5,2.7),pch=15,cex=1.5,col = colCause)
  # text(c(3.5,3.5),c(3.5,2.7),c("Legal\nculling","Other\nmortality"),cex=2,pos=4)
  # dev.off()
  # 
  # 
  # 
  # ## ------     4.4.3. RECRUITMENT ------
  # 
  # ## ------       4.4.3.1. PLOT PER CAPITA RECRUITMENT -----
  # 
  # # pdf(file = file.path(working.dir,"figures/erCapita.pdf"),width = 10,height = 8)
  # grDevices::png(filename = file.path( working.dir, "figures/PerCapita.png"),
  #                width = 10, height = 8, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # ##-- PER CAPITA RECRUITMENT
  # par(mar = c(4.5,4.5,1,1), xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
  # plot(10,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,0.5),
  #      type ="n", xaxt = "n", xlab = "Years", ylab = "Per-capita recruitment")
  # axis(2, tck = -0.02)
  # axis(1, 1:(n.years-1),
  #      labels = paste(years[1:(n.years-1)], years[2:n.years], sep = " to "),
  #      cex.axis = 0.8)
  # abline(v = 1:(n.years-1)+0.5,lty=2)
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- Number of recruits at t ==> ids with state 1 at t and 2 at t+1
  #   n.recruit <- apply(resultsSXYZ_MF$sims.list$z[ , ,c(t,t+1)], 1, function(x) sum(x[ ,1]%in%1 & x[ ,2]%in%2))
  #   ##-- Number of reproducing ids at t ==> ids with state 2 at t
  #   alivetminus1 <- apply(resultsSXYZ_MF$sims.list$z[ , ,t], 1, function(x)sum(x %in% 2))
  #   ##-- Plot quantiles
  #   plotQuantiles(x = n.recruit/alivetminus1,
  #                 at = t,
  #                 width = 0.3,
  #                 col = colCause[1])
  # }#t
  # dev.off()
  # 
  # 
  # 
  # ## ------       4.4.3.2. PLOT NUMBER OF RECRUITS ----
  # 
  # ##-- Identify individual status
  # isAvail <- resultsSXYZ_MF$sims.list$z == 1
  # isAlive <- resultsSXYZ_MF$sims.list$z == 2
  # 
  # ##-- Identify individual sex
  # isFemale <- resultsSXYZ_MF$sims.list$sex == "F"
  # isMale <- resultsSXYZ_MF$sims.list$sex == "M"
  # 
  # N_recruit <- N_recruit_F <- N_recruit_M <- matrix(NA, n.mcmc, n.years-1)
  # for(t in 1:(n.years-1)){
  #   for(iter in 1:n.mcmc){
  #     N_recruit_F[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & isFemale)
  #     N_recruit_M[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & isMale)
  #     N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]
  #   }#iter
  # }#t
  # 
  # 
  # # pdf(file = file.path(working.dir,"figures/NumRecruitTotal.pdf"),width = 10,height = 6)
  # grDevices::png( filename = file.path(working.dir, "figures/NumRecruitTotal.png"),
  #                 width = 10, height = 6, units = "in", pointsize = 12,
  #                 res = 300, bg = NA)
  # 
  # nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
  #              widths = c(0.05,1,0.30),
  #              heights = c(0.15,1))
  # 
  # par(mar = c(4.5,4.5,1,1), xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
  # 
  # ymax <- 10*(trunc(max(N_recruit)/10)+1)
  # 
  # plot(10,
  #      xlim = c(0.5, n.years-0.5), ylim = c(0,ymax),
  #      type ="n", xaxt = "n",
  #      xlab = "Years", ylab = "Number of recruits")
  # axis(2, tck = -0.02)
  # axis(1, 1:(n.years-1),
  #      labels = years[2:n.years], cex.axis=0.8)
  # abline(v = 1:(n.years-1)+0.5,lty=2)
  # for(t in 1:(n.years-1)){
  #   plotQuantiles(x = N_recruit[ ,t],
  #                 at = t,
  #                 width = 0.3,
  #                 col = colSex[3])
  # 
  #   plotQuantiles(x = N_recruit_F[ ,t],
  #                 at = t - diffSex,
  #                 col = colSex[1])
  # 
  #   plotQuantiles(x = N_recruit_M[ ,t],
  #                 at = t + diffSex,
  #                 col = colSex[2])
  # }#t
  # 
  # ##-- LEGEND
  # par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  # plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
  # points(c(4,4,4), c(4,3,2), pch = 15, cex = 5.5, col = colSex)
  # points(c(4,4,4), c(4,3,2), pch = 15, cex = 3, col = colSex)
  # text(c(5.3,5.3,5.3), c(4,3,2),
  #      c("Females", "Males", "Total"), cex = 2, pos = 4)
  # 
  # dev.off()
  # 
  # 

  # ## ------   4.6. p0 ------
  # 
  # message("## Plotting detectability...")
  # 
  # ## ------     4.6.1. p0 bars ------
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("p0_bars.pdf")),
  # #     width = 8, height = 12)
  # grDevices::png(filename = file.path(working.dir, "figures/p0_bars.png"),
  #                width = 8, height = 12, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # nf <- layout(rbind(c(1,2),
  #                    c(3,4),
  #                    c(5,6)),
  #              widths = c(1,0.5),
  #              heights = 1)
  # 
  # for(c in 1:nrow(COUNTIES)){
  #   par(mar=c(4,4,1,1), tck=0)
  # 
  #   plot(10, xlim = c(0.5, n.years+0.5), ylim = c(0,0.01), type ="n", xaxt="n",
  #        xlab = "Years", ylab = "Baseline detection probability")
  # 
  #   axis(1, 1:n.years, labels = years)
  #   axis(2, tck = -0.02)
  #   abline(v = 1:(n.years-1) + 0.5, lty = 2)
  # 
  #   for(t in 1:n.years){
  #     plotQuantiles( x = results_F$sims.list$p0[ ,COUNTIES$index == c,t],
  #                    at = t - diffSex,
  #                    col = colSex[1])
  # 
  #     plotQuantiles( x = results_M$sims.list$p0[ ,COUNTIES$index == c,t],
  #                    at = t + diffSex,
  #                    col = colSex[2])
  #   }#t
  # 
  #   ##-- LEGEND
  #   if(c == 1){
  #     polygon(x = c(0.5,4.2,4.2,0.5),
  #             y = c(0.006,0.006,0.01,0.01),
  #             col = adjustcolor("white", alpha.f = 0.9),
  #             border = NA)
  #     points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 5.5, col = colSex)
  #     points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 3, col = colSex)
  #     text(c(1.2,1.2),c(0.0077,0.0092),  c("Females", "Males"), cex = 2, pos = 4)
  #   }
  # 
  # 
  #   par(mar = c(0,0,0,0))
  #   plot(st_geometry(COUNTIES), border = grey(0.5), col = grey(0.5), lwd = 0.1)
  #   plot(st_geometry(COUNTIES[COUNTIES$index == c, ]),
  #        add = T, col = adjustcolor("red",0.5), border = "red")
  #   text(COUNTIES[COUNTIES$index == c, ],
  #        labels = COUNTIES$Name[COUNTIES$index == c],
  #        col = "white")
  # }#c
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------     4.6.2. p0 maps ------
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("p0_maps.pdf")),
  # #     width = 10, height = 6)
  # grDevices::png(filename = file.path(working.dir, "figures/p0_maps.png"),
  #                width = 10, height = 6, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # for(t in 1:n.years){
  #   par(mfrow = c(1,2))
  # 
  #   ##-- FEMALE
  #   detectors$main.detector.sp$p0_F <-
  #     ilogit(logit(results_F$mean$p0[nimConstants$county,t]) +
  #              results_F$mean$betaDet[1] * nimDataF$detCovs[ ,1,t] +
  #              results_F$mean$betaDet[2] * nimDataF$detCovs[ ,2,t])
  # 
  #   p0_F.R <- raster::rasterFromXYZ(cbind(detectors$main.detector.sp$main.cell.x,
  #                                         detectors$main.detector.sp$main.cell.y,
  #                                         detectors$main.detector.sp$p0_F))
  # 
  #   raster::plot(p0_F.R,
  #                main =  paste0("Females ", years[t]),
  #                legend.args = list(text = 'p0',
  #                                   side = 4, font = 2, line = 2.5, cex = 0.8))
  # 
  #   ##-- MALE
  #   detectors$main.detector.sp$p0_M <-
  #     ilogit(logit(results_M$mean$p0[nimConstants$county,t]) +
  #              results_M$mean$betaDet[1] * nimDataM$detCovs[ ,1,t] +
  #              results_M$mean$betaDet[2] * nimDataM$detCovs[ ,2,t])
  # 
  #   p0_M.R <- raster::rasterFromXYZ(cbind(detectors$main.detector.sp$main.cell.x,
  #                                         detectors$main.detector.sp$main.cell.y,
  #                                         detectors$main.detector.sp$p0_M))
  #   raster::plot(p0_M.R,
  #                main = paste0("Males ", years[t]),
  #                legend.args = list(text = 'p0',
  #                                   side = 4, font = 2, line = 2.5, cex = 0.8))
  # }#t
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------     4.6.3. p0 betas ------
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("p0_beta.pdf")),
  # #     width = 6, height = 4)
  # grDevices::png(filename = file.path(working.dir, "figures/p0_beta.png"),
  #                width = 6, height = 4, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # nf <- layout(cbind(c(6,3),c(4,1),c(5,2)),
  #              widths = c(0.05,1,0.30),
  #              heights = c(0.15,1))
  # 
  # ##-- PLOT BETAS
  # par(mar = c(5,4.5,0.5,0.5), tck = 0, xaxs = "i", cex.axis = 1.3, cex.lab = 1.6)
  # plot( 10, xlim = c(0.5, 2.5), ylim = c(-5,5),
  #       type = "n", xaxt = "n", xlab = "Years", ylab = "beta")
  # axis(1, c(1,2), labels =  c("distance \nto roads","presence of \nother obs."), padj = 0.5)
  # axis(2, tck = -0.02)
  # abline(v = 1.5, lty = 2)
  # abline(h = 0, lty = 1)
  # for(b in 1:2){
  #   plotQuantiles( results_F$sims.list$betaDet[ ,b],
  #                  at = b - diffSex,
  #                  col = colSex[1])
  #   plotQuantiles( results_M$sims.list$betaDet[ ,b],
  #                  at = b + diffSex,
  #                  col = colSex[2])
  # }#b
  # 
  # ##-- LEGEND
  # par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  # plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
  # points(c(4,4), c(4,3), pch = 15, cex = 5.5, col = colSex)
  # points(c(4,4), c(4,3), pch = 15, cex = 3, col = colSex)
  # text(c(5.3,5.3), c(4,3),  c("Females", "Males"), cex = 1.5, pos = 4)
  # dev.off()
  # 
  # 
  # 

  ## ------   4.7. NGS, Dead recoveries & Carnivore obs ------
  
  ##-- Plot NGS & Dead recovery maps
  # pdf(file = file.path(working.dir, "figures", "NGS_DR_maps.pdf"),
  #     width = 18, height = 12)
  grDevices::png(filename = file.path(working.dir, "figures/NGS_DR_maps.png"),
                 width = 18, height = 12, units = "in", pointsize = 12,
                 res = 300, bg = NA)
  
  ##-- layout
  mx <- rbind(c(1,rep(1:5, each = 2)),
              c(rep(1:5, each = 2), 5))
  mx <- rbind(mx, mx + 5)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  par(mar = c(0,0,0,0))
  for(t in 1:length(years)){
    plot(sf::st_geometry(COUNTIES), border = NA, col = "gray80")
    points(data.alive$data.sp[data.alive$data.sp$Year == years[t], ],
           pch = 3, col = "orange", lwd = 0.7)
    points(data.dead[data.dead$Year == years[t], ],
           pch = 3, col = "slateblue", lwd = 0.7)
    mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
    
    if(t == n.years){
      ##-- LEGEND
      xLeg <- 830000
      yLeg <- 6730000
      segments(x0 = xLeg, x1 = xLeg,
               y0 = yLeg, y1 = yLeg + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
      text(xLeg-80000, yLeg+500000/2, labels = "500 km", srt = 90, cex = 2)
      
      points(x = c(xLeg-200000,xLeg-200000),
             y = c(yLeg-100000,yLeg-180000),
             pch = 3, lwd = 1.5, cex = 3,
             col = c("orange","slateblue"))
      text(x = c(xLeg-150000,xLeg-150000),
           y = c(yLeg-100000,yLeg-180000),
           c("NGS samples", "Dead recoveries"), cex = 2, pos = 4)
    }#if
  }#t
  dev.off()
  
  # ##-- Plot Carnivore observations maps
  # pdf(file = file.path(working.dir, "figures", paste0("CarnivoreObs_maps_classic.pdf")),
  #     width = 18, height = 12)
  # grDevices::png(filename = file.path(working.dir, "figures/CarnivoreObs_maps_classic.png"),
  #     width = 18, height = 12, units = "in", pointsize = 12,
  #     res = 300, bg = NA)
  #
  # ##-- layout
  # mx <- rbind(c(1,rep(1:5, each = 2)),
  #             c(rep(1:5, each = 2), 5))
  # mx <- rbind(mx, mx + 5)
  # nf <- layout(mx,
  #              widths = c(rep(1,ncol(mx))),
  #              heights = rep(1,2))
  # par(mar = c(0,0,0,0))
  # for(t in 1:length(years)){
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
  #   image(mask(ds.brickCont[[t]],COUNTRIESsimpFig[1,]), add = TRUE, col = c("white","forestgreen"), legend = FALSE)
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = "gray80", col = NA, add = TRUE)
  #   
  #   mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
  #   
  #   if(t == n.years){
  #     segments(x0 = 830000, x1 = 830000,
  #              y0 = 6730000, y1 = 6730000 + 500000,
  #              col = grey(0.3), lwd = 4, lend = 2)
  #     text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
  #     
  #     ##-- LEGEND
  #     par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  #     plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
  #   }#if
  # }#t
  # dev.off()
  
  
  
  ## ------ 5. TABLES -----
  
  gc(verbose = FALSE)
  
  ## ------   5.1. ABUNDANCE ------
  
  ##-- Set-up names for tables
  
  regionNames_NOR <- sort(row.names(ACdensity[[t]]$summary)[grep("Region",row.names(ACdensity[[t]]$summary))])
  countyNames_NOR <- sort(unique(factorValues(rrCounties, rrCounties[rrRegions[] > 3], layer = 1)[,1]))
  countyNames_North <- sort(unique(factorValues(rrCounties, rrCounties[rrRegions[] %in% 1], layer=1)[,1]))
  countyNames_Middle <- sort(unique(factorValues(rrCounties, rrCounties[rrRegions[] %in% 2], layer=1)[,1]))
  countyNames_South <- sort(unique(factorValues(rrCounties, rrCounties[rrRegions[] %in% 3], layer=1)[,1]))
  
  rownames_Table <- c("Total",
                      "Norway", regionNames_NOR,
                      "Sweden",
                      "Nordre", countyNames_North,
                      "Midtre", countyNames_Middle,
                      "Söndre", countyNames_South)
  
  
  ##-- Fix names for .tex table and add grey color for Norrbotten
  countyNames_North_tex <- countyNames_North
  countyNames_North_tex[countyNames_North_tex %in% "Norrbotten"] <- "\\textcolor[gray]{.5}{Norrbotten}"
  
  rownames_Table_tex1 <- c("TOTAL",
                           "\\hspace{0.25cm}NORWAY",
                           paste0("\\hspace{0.5cm}", regionNames_NOR),
                           "\\hspace{0.25cm}SWEDEN",
                           "\\hspace{0.5cm}Norra",
                           paste0("\\hspace{0.75cm}", countyNames_North),
                           "\\hspace{0.5cm}Mellersta",
                           paste0("\\hspace{0.75cm}", countyNames_Middle),
                           "\\hspace{0.5cm}Södra",
                           paste0("\\hspace{0.75cm}", countyNames_South))
  
  rownames_Table_tex2 <- c("TOTAL",
                           "\\hspace{0.25cm}NORWAY",
                           paste0("\\hspace{0.5cm}", regionNames_NOR),
                           "\\hspace{0.25cm}SWEDEN",
                           "\\hspace{0.5cm}Norra",
                           paste0("\\hspace{0.75cm}", countyNames_North_tex),
                           "\\hspace{0.5cm}Mellersta",
                           paste0("\\hspace{0.75cm}", countyNames_Middle),
                           "\\hspace{0.5cm}Södra",
                           paste0("\\hspace{0.75cm}", countyNames_South))
  
  
  
  ## ------     5.1.1. ALL YEARS, BOTH SEX COMBINED ------
  
  ##-- Create table to store abundance & CI
  NCarRegionEstimates <- matrix("", ncol = n.years, nrow = length(rownames_Table))
  row.names(NCarRegionEstimates) <- rownames_Table
  colnames(NCarRegionEstimates) <- years+1
  
  ##-- Fill in the table 
  for(t in 1:n.years){
    NCarRegionEstimates[rownames_Table,t] <- paste0(
      round(ACdensity[[t]]$summary[rownames_Table,"mean"],digits = 1)," (",
      round(ACdensity[[t]]$summary[rownames_Table,"95%CILow"],digits = 0),"-",
      round(ACdensity[[t]]$summary[rownames_Table,"95%CIHigh"],digits = 0),")")
  }#t
  
  # ##-- Quick check to make sure values sums up 
  # tmp <- ACdensity[[t]]$summary[1:(nrow(ACdensity[[t]]$summary)),]
  # # SWE
  # row.names(ACdensity[[t]]$summary)
  # sum(tmp[countyNames_North,"mean"])+
  #   sum(tmp[countyNames_Middle,"mean"])+
  #   sum(tmp[countyNames_South,"mean"])
  # sum(tmp[c("Nordre","Midtre","Söndre"),"mean"])
  # tmp["Sweden","mean"]
  # #NOR
  # sum(tmp[regionNames_NOR,"mean"])
  # tmp["Norway","mean"]
  # #TOTAL
  # tmp["Sweden","mean"]+tmp["Norway","mean"]
  # tmp["Total","mean"]
  
  ##-- Export .csv
  write.csv( NCarRegionEstimates,
             file = file.path(working.dir, "tables/NAllYears.csv"),
             fileEncoding = "latin1")
  
  ##-- Add grey color for years without sampling in Norrbotten
  NCarRegionEstimates["Norrbotten",yearsNotSampled] <- paste0("\\textcolor[gray]{.5}{",NCarRegionEstimates["Norrbotten", yearsNotSampled], "*}")
  NCarRegionEstimates["Nordre",yearsNotSampled] <- paste0("\\textcolor[gray]{.5}{",NCarRegionEstimates["Nordre",yearsNotSampled], "**}")
  NCarRegionEstimates["Sweden",yearsNotSampled] <- paste0("\\textcolor[gray]{.5}{",NCarRegionEstimates["Sweden",yearsNotSampled], "**}")
  NCarRegionEstimates["Total",yearsNotSampled] <- paste0("\\textcolor[gray]{.5}{",NCarRegionEstimates["Total",yearsNotSampled], "**}")
  
  ##-- Fix row names
  row.names(NCarRegionEstimates) <- rownames_Table_tex2
  
  ##-- Export .tex
  print( xtable( NCarRegionEstimates,
                 type = "latex",
                 align = paste(c("l",rep("c",ncol(NCarRegionEstimates))), collapse = "")),
         # scalebox=.8,
         floating = FALSE,
         sanitize.text.function = function(x){x},
         add.to.row = list( list(seq(1,nrow(NCarRegionEstimates), by = 2)),
                            "\\rowcolor[gray]{.96}"),
         file = file.path(working.dir, "tables/NCountiesCarnivoreRegions.tex"))
  
  
  
  ## ------     5.1.2. LAST YEAR, N PER SEX PER REGION ------
  
  NCountyEstimatesLastRegions <- matrix("", ncol = 3, nrow = length(rownames_Table))
  row.names(NCountyEstimatesLastRegions) <- rownames_Table
  colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")
  
  ##-- Fill in table 
  ##-- FEMALES
  NCountyEstimatesLastRegions[ ,"Females"] <- paste0(
    round(ACdensityF[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensityF[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensityF[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  ##-- MALES 
  NCountyEstimatesLastRegions[ ,"Males"] <- paste0(
    round(ACdensityM[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensityM[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensityM[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  ##-- TOTAL 
  NCountyEstimatesLastRegions[ ,"Total"] <- paste0(
    round(ACdensity[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensity[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensity[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  
  ##-- Export .csv
  write.csv( NCountyEstimatesLastRegions,
             file = file.path(working.dir, "tables/N_LastYearPerSex_region.csv"),
             fileEncoding = "latin1")
  
  ##-- Fix row names
  row.names(NCountyEstimatesLastRegions) <- rownames_Table_tex1
  
  ##-- Export .tex 
  print(xtable( NCountyEstimatesLastRegions,
                type = "latex",
                align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
        sanitize.text.function = function(x){x},
        # scalebox=.8,
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions), by = 2)),
                          "\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_LastYearPerSex_region.tex"))
  
  
  
  ## ------     5.1.3. LAST YEAR, N PER SEX PER REGION WITH PROPORTION OF AREA COVERED ------
  
  NCountyEstimatesLastRegions <- matrix("", ncol = 4, nrow = length(rownames_Table))
  row.names(NCountyEstimatesLastRegions) <- rownames_Table
  colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total","\\% Area")
  
  ##-- Fill in table 
  ##-- FEMALES
  NCountyEstimatesLastRegions[ ,"Females"] <- paste0(
    round(ACdensityF[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensityF[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensityF[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  ##-- MALES 
  NCountyEstimatesLastRegions[ ,"Males"] <- paste0(
    round(ACdensityM[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensityM[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensityM[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  ##-- TOTAL 
  NCountyEstimatesLastRegions[ ,"Total"] <- paste0(
    round(ACdensity[[n.years]]$summary[rownames_Table,"mean"], digits = 1)," (",
    round(ACdensity[[n.years]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
    round(ACdensity[[n.years]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  
  ##-- AREA
  NCountyEstimatesLastRegions[ ,"\\% Area"] <- round(percAllRegions[rownames_Table]*100, digits = 0)
  ##-- Round up
  NCountyEstimatesLastRegions[NCountyEstimatesLastRegions[,4] %in% c("98","99"),4] <- 100
  
  ##--  Export .csv
  write.csv( NCountyEstimatesLastRegions,
             file = file.path(working.dir, "tables", "NLastYearPerSexArea.csv"),
             fileEncoding = "latin1")
  
  
  ##--  Export .tex
  row.names(NCountyEstimatesLastRegions) <- rownames_Table_tex1
  
  print(xtable(NCountyEstimatesLastRegions,
               type = "latex",
               align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions)-1),"||c"),collapse = "")),
        sanitize.text.function=function(x){x},
        # scalebox=.8,
        floating = FALSE,
        add.to.row = list( list(seq(1,nrow(NCountyEstimatesLastRegions), by = 2)),
                           "\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/NCountiesSexLastYearRegionsArea.tex"))
  
  
  
  # ## ------     5.1.4. 2 LAST YEARS, N PER SEX PER REGION ------
  # 
  # NCountyEstimatesLast2Regions <- matrix("", ncol = 6, nrow = length(rownames_Table))
  # row.names(NCountyEstimatesLast2Regions) <- rownames_Table
  # colnames(NCountyEstimatesLast2Regions) <- c(paste("Females", years[n.years-1]),
  #                                             paste("Males", years[n.years-1]),
  #                                             paste("Total", years[n.years-1]),
  #                                             paste("Females", years[n.years]),
  #                                             paste("Males", years[n.years]),
  #                                             paste("Total", years[n.years]))
  # 
  # ##-- FILL IN TABLE 
  # for(t in (n.years-1):n.years){
  #   ##-- FEMALES
  #   NCountyEstimatesLast2Regions[ ,paste("Females", years[t])] <- paste0(
  #     round(ACdensityF[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
  #     round(ACdensityF[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
  #     round(ACdensityF[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  #   
  #   ##-- MALES 
  #   NCountyEstimatesLast2Regions[ ,paste("Males", years[t])] <- paste0(
  #     round(ACdensityM[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
  #     round(ACdensityM[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
  #     round(ACdensityM[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  #   
  #   ##-- TOTAL 
  #   NCountyEstimatesLast2Regions[ ,paste("Total", years[t])] <- paste0(
  #     round(ACdensity[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
  #     round(ACdensity[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
  #     round(ACdensity[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  # }#t
  # 
  # ##-- Export .csv
  # write.csv( NCountyEstimatesLast2Regions,
  #            file = file.path(working.dir, "tables/NLast2YearsPerSex.csv"),
  #            fileEncoding = "latin1")
  # 
  # 
  # ##-- Fix row names
  # row.names(NCountyEstimatesLast2Regions) <- rownames_Table_tex2
  # 
  # ##-- Fix col names
  # NCountyEstimatesLast2Regions <- rbind( c("F","M","Total","F","M","Total"),
  #                                        NCountyEstimatesLast2Regions)
  # 
  # ##-- Export .tex 
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # uniqueYEAR <- c( paste(unlist(YEARS[n.years-1]),collapse = "/"),
  #                  paste(unlist(YEARS[n.years]),collapse = "/"))
  # addtorow$command <- c(paste0(paste0('& \\multicolumn{3}{c}{', uniqueYEAR,
  #                                     '}', collapse=''), '\\\\'),
  #                       rep("\\rowcolor[gray]{.95}",1))
  # print(xtable( NCountyEstimatesLast2Regions, 
  #               type = "latex",
  #               align = paste(c("l",rep("c",3),"|", rep("c",3)),collapse = "")),
  #       sanitize.text.function = function(x){x},
  #       # scalebox=.8,
  #       floating = FALSE,
  #       add.to.row = addtorow,
  #       include.colnames = F,
  #       file = file.path(working.dir, "tables/NCountiesSexLast2YearsRegions.tex"))
  # 
  # 
  # 
  ## ------     5.1.5. ALL YEARS, N PER SEX PER COUNTY ------
  
  NCountyEstimatesAllSexRegions <- matrix("", ncol = n.years*3, nrow = length(rownames_Table)+1)
  row.names(NCountyEstimatesAllSexRegions) <- c("",rownames_Table)
  colnames(NCountyEstimatesAllSexRegions) <- rep((years+1), each = 3)
  NCountyEstimatesAllSexRegions[1, ] <- rep(c("Females","Males","Total"), n.years)
  
  ## FILL IN TABLE 
  for(t in 1:n.years){
    ##-- Identify columns for that year
    cols <- which(colnames(NCountyEstimatesAllSexRegions) %in% (years+1)[t])
    
    #-- Identify the column for each sex that year
    colsF <- which(NCountyEstimatesAllSexRegions[1,cols] %in% "Females")
    colsM <- which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
    colsT <- which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
    
    ##-- FEMALES
    NCountyEstimatesAllSexRegions[rownames_Table,cols[colsF]] <- paste0(
      round(ACdensityF[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
      round(ACdensityF[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
      round(ACdensityF[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
    
    ##-- MALES 
    NCountyEstimatesAllSexRegions[rownames_Table,cols[colsM]] <- paste0(
      round(ACdensityM[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
      round(ACdensityM[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
      round(ACdensityM[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
    
    ##-- TOTAL 
    NCountyEstimatesAllSexRegions[rownames_Table,cols[colsT]] <- paste0(
      round(ACdensity[[t]]$summary[rownames_Table,"mean"], digits = 1)," (",
      round(ACdensity[[t]]$summary[rownames_Table,"95%CILow"], digits = 0),"-",
      round(ACdensity[[t]]$summary[rownames_Table,"95%CIHigh"], digits = 0),")")
  }#t
  
  ##-- Export .csv
  write.csv( NCountyEstimatesAllSexRegions,
             file = file.path(working.dir, "tables/NAllYearsPerSex.csv"),
             fileEncoding = "latin1")
  
  ##-- Add grey color for years without sampling in Norrbotten
  index <- which(colnames(NCountyEstimatesAllSexRegions) %in% (years[yearsNotSampled]+1))
  NCountyEstimatesAllSexRegions["Norrbotten",index] <- paste0("\\textcolor[gray]{.5}{",NCountyEstimatesAllSexRegions["Norrbotten", index], "*}")
  NCountyEstimatesAllSexRegions["Nordre",index] <- paste0("\\textcolor[gray]{.5}{",NCountyEstimatesAllSexRegions["Nordre",index], "**}")
  NCountyEstimatesAllSexRegions["Sweden",index] <- paste0("\\textcolor[gray]{.5}{",NCountyEstimatesAllSexRegions["Sweden",index], "**}")
  NCountyEstimatesAllSexRegions["Total",index] <- paste0("\\textcolor[gray]{.5}{",NCountyEstimatesAllSexRegions["Total",index], "**}")
  
  ##-- Fix row names
  row.names(NCountyEstimatesAllSexRegions) <- c("", rownames_Table_tex2)
  
  ##-- Export .tex WRITE LATEX
  print(xtable(NCountyEstimatesAllSexRegions,
               type = "latex",
               align = paste(c("l",rep("c",ncol(NCountyEstimatesAllSexRegions))),collapse = "")),
        sanitize.text.function=function(x){x},
        # scalebox=.8,
        floating = FALSE,
        add.to.row=list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),
                        "\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/NAllYearsPerSex.tex"))
  
  
  
  ## ------     5.1.6. ALL YEARS, BOTH SEX COMBINED NORWAY ------
  
  rownames_Table_NOR <- c("Norway", countyNames_NOR)
  
  ##-- Create table to store abundance and CIs
  NCarRegionEstimatesNOR <- matrix("", ncol = n.years, nrow = length(rownames_Table_NOR))
  row.names(NCarRegionEstimatesNOR) <- rownames_Table_NOR
  colnames(NCarRegionEstimatesNOR) <- years+1
  
  ##-- Fill in the table 
  for(t in 1:n.years){
    NCarRegionEstimatesNOR[rownames_Table_NOR,t] <- paste0(
      round(ACdensity[[t]]$summary[rownames_Table_NOR,"mean"],digits = 1)," (",
      round(ACdensity[[t]]$summary[rownames_Table_NOR,"95%CILow"],digits = 0),"-",
      round(ACdensity[[t]]$summary[rownames_Table_NOR,"95%CIHigh"],digits = 0),")")
  }#t
  
  ##-- Export .csv
  row.names(NCarRegionEstimatesNOR) <- c("TOTAL", countyNames_NOR)
  write.csv( NCarRegionEstimatesNOR,
             file = file.path(working.dir, "tables/NAllYearsNorwegianCounties.csv"),
             fileEncoding = "latin1")
  
  ##-- Export .tex
  row.names(NCarRegionEstimatesNOR) <- c("TOTAL", paste0("\\hspace{0.25cm} ", countyNames_NOR))
  print(xtable( NCarRegionEstimatesNOR, 
                type = "latex",
                align = paste(c("l",rep("c",ncol(NCarRegionEstimatesNOR))),collapse = "")),
        # scalebox=.8,
        floating = FALSE,
        sanitize.text.function=function(x){x},
        add.to.row = list( list(seq(1,nrow(NCarRegionEstimatesNOR), by = 2)),
                           "\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/NAllYearsNorwegianCounties.tex"))
  
  
  
  
  ## ------   5.2. DATA SUMMARY ------

  ##-- SOME TALLIES TO CHECK THINGS
  ##-- NGS
  NGS <- data.alive$data.sp
  dead <- data.dead
  
  ##-- FOR REPORT SUMMARY
  dataSummary <- cbind.data.frame(
    "female" = c(length(NGS$Id[NGS$Sex=="female"]),
                 length(unique(dead$Id[dead$Sex=="female"])),
                 length(unique(c(NGS$Id[NGS$Sex=="female"],dead$Id[dead$Sex=="female"])))),
    "male" = c(length(NGS$Id[NGS$Sex=="male"]),
               length(unique(dead$Id[dead$Sex=="male"])),
               length(unique(c(NGS$Id[NGS$Sex=="male"],dead$Id[dead$Sex=="male"])))),
    "Total" = c(length(NGS$Id),
                length(dead$Id),
                length(unique(c(NGS$Id,dead$Id)))))
  row.names(dataSummary) <- c("N_NGS", "N_DR", "N_IDs")
  ##-- print .csv
  write.csv(dataSummary, file = file.path(working.dir, "tables/dataSummary.csv"))

  
  ##-- Prepare raster of countries
  countryRaster <- habitatRasterResolution$`5km`[["Countries"]]
  
  
  
  ## ------     5.2.1. NGS SAMPLES & IDs ------
  
  NGS_SEX <- matrix("", ncol = n.years*2, nrow = 3)
  row.names(NGS_SEX) <- c( "",
                           "number of NGS samples",
                           "number of NGS individuals")
  colnames(NGS_SEX) <- rep(years, each = 2)
  NGS_SEX[1, ] <- rep(c("F","M"), n.years)
  
  sex <- c("female","male")
  sex1 <- c(0,1)
  ye <- seq(1, n.years*2, by = 2)
  for(s in 1:2){
    for(t in 1:n.years){
      temp <- NGS[NGS$Year == years[t] & NGS$Sex == sex[s], ]
      NGS_SEX["number of NGS samples", ye[t] + sex1[s]] <- nrow(temp)
      NGS_SEX["number of NGS individuals", ye[t] + sex1[s]] <- length(unique(temp$Id))
    }#t
  }#s
  
  ##-- Remove unnecessary objects from memory
  rm(list = c( "temp"))
  gc(verbose = FALSE)    
  
  ##-- print .csv
  write.csv( NGS_SEX, file = file.path(working.dir, "tables", paste0("NGS_SEX.csv")))
  
  ##-- print .tex
  addtorow <- list()
  addtorow$pos <- list(c(0),0)
  addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGS_SEX))),
                                      '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
  colnames(NGS_SEX) <- rep("", ncol(NGS_SEX))
  print(xtable::xtable(NGS_SEX, type = "latex",
                       align = paste(c("l",rep("c",ncol(NGS_SEX))), collapse = "")),
        floating = FALSE, include.colnames = FALSE,
        add.to.row = addtorow,
        file = file.path(working.dir, "tables", paste0("NGS_SEX.tex")))
  
  
  
  
  ## ------     5.2.2. DEAD RECOVERIES by CAUSE ------
  
  Dead_SEX <- matrix(0, ncol = n.years*2+1, nrow = 6)
  row.names(Dead_SEX) <- c("","other","other","legal culling","legal culling","")
  colnames(Dead_SEX) <- c("",unlist(lapply(years, function(x) c(x,x))))
  Dead_SEX[1,] <- c("",rep(c("F","M"),n.years))
  Dead_SEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
  sex <- c("female","male")
  sex1 <- c(0,1)
  ye <- seq(1,n.years*2,by=2)
  
  ##-- Separate mortalities
  cause <- c("other","legal culling")
  for(t in 1:n.years){
    for(s in 1:2){
      for(d in 1:2){
        if(d==1){
          temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !dead$Legal, ]
        } else {
          temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$Legal, ]
        }
        row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[ ,1]=="Norway")
        Dead_SEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country_sample %in% "(N)" ]))
        
        row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[,1]=="Sweden" )
        Dead_SEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country_sample %in% "(S)"]))
      }#t
      Dead_SEX[6, ye[t] + sex1[s]+1] <-  sum(as.numeric(Dead_SEX[2:6,ye[t] + sex1[s]+1]))
    }
  }
  ##-- Remove unnecessary objects from memory
  rm(list = c( "temp"))
  gc(verbose = FALSE)    
  
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
  
  print(xtable::xtable(Dead_SEX, type = "latex",
                       align = rep("c", ncol(Dead_SEX)+1)),
        floating = FALSE,
        add.to.row = addtorow,
        include.colnames = FALSE,
        include.rownames = FALSE,
        sanitize.text.function = function(x){x},
        file = file.path(working.dir, "tables/DeadidCountrySEX.tex"))
  
  
  
  ## ------     5.2.3. PROPORTION OF INDIVIDUALS DETECTED OVERALL ------
  
  ##-- Get the number of individuals detected each year
  n.detected_F <- apply(nimDataF$detNums+nimDataF$detNumsOth, 2, function(x)sum(x>0))
  n.detected_M <- apply(nimDataM$detNums+nimDataM$detNumsOth, 2, function(x)sum(x>0))

  propDetected <- matrix("", ncol = n.years, nrow = 3)
  row.names(propDetected) <- c("F","M","Total")
  colnames(propDetected) <- years
  for(t in 1:n.years){
    propDetected["F",t] <- getCleanEstimates(n.detected_F[t]/colSums(ACdensityF[[t]]$PosteriorAllRegions))
    propDetected["M",t] <- getCleanEstimates(n.detected_M[t]/colSums(ACdensityM[[t]]$PosteriorAllRegions))
    propDetected["Total",t] <- getCleanEstimates((n.detected_F[t]+n.detected_M[t])/
                                                   (colSums(ACdensityF[[t]]$PosteriorAllRegions)+
                                                      colSums(ACdensityM[[t]]$PosteriorAllRegions)))
  }#t
  
  ##-- Remove unnecessary objects from memory
  rm(list = c( "n.detected_F", "n.detected_M"))
  gc(verbose = FALSE)    
  
  ##-- print .csv
  write.csv(propDetected,
            file = file.path(working.dir, "tables/PropDetectedIds.csv"))
  
  ##-- print .tex
  print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
        floating = FALSE, sanitize.text.function=function(x){x},
        add.to.row = list(list(seq(1,nrow(propDetected), by = 2)),"\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/PropDetectedIds.tex"))
  
  
  
  ## ------     5.2.4. PROPORTION OF THE POPULATION DETECTED ------
  
  ##-- Extract number of individuals detected
  isDetected <- rbind(nimDataM$detNums+nimDataM$detNumsOth,
                      nimDataF$detNums+nimDataF$detNumsOth) > 0
  
  ##-- Identify individual sex
  isFemale <- resultsSXYZ_MF$sims.list$sex == "F"
  isMale <- resultsSXYZ_MF$sims.list$sex == "M"
  
  ##-- Identify individual status
  #isAvail <- resultsSXYZ_MF$sims.list$z == 1 
  isAlive <- resultsSXYZ_MF$sims.list$z == 2 
  
  ##-- Calculate % of the wolverine population detected 
  prop <- matrix(NA,3,n.years)
  dimnames(prop) <- list("% individuals" = c("F","M","Total"),
                         "Years" = c(years))
  for(t in 1:n.years){
    prop_F <- prop_M <- prop_tot <- rep(NA,n.mcmc)
    for(iter in 1:n.mcmc){
      
      country <- countryRaster[raster::cellFromXY(countryRaster, resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
      isIn <- country %in% c(2,4)
      
      ##-- Detected female
      N_F <- sum(isAlive[iter, ,t] & isIn & isFemale)
      N_det_F <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isIn & isFemale)
      prop_F[iter] <- N_det_F / N_F 
      
      ##-- Detected male
      N_M <- sum(isAlive[iter, ,t] & isIn & isMale)
      N_det_M <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isIn & isMale)
      prop_M[iter] <- N_det_M / N_M 
      
      ##-- Detected total
      prop_tot[iter] <- (N_det_F + N_det_M) / (N_F + N_M) 
    }#iter
    
    prop["F",t] <- getCleanEstimates(prop_F)
    prop["M",t] <- getCleanEstimates(prop_M)
    prop["Total",t] <- getCleanEstimates(prop_tot)
  }#t
  
  ##-- Remove unnecessary objects from memory
  rm(list = c( "N_F", "N_M",
               "N_det_F", "N_det_M",
               "prop_F", "prop_M", "prop_tot"))
  gc(verbose = FALSE)  
  

  ##-- print .csv
  write.csv(prop,
            file = file.path(working.dir, "tables/PropDetected.csv"))
  
  if(n.years > 8){
    ##-- Print .tex (split in two tables to print in the overleaf document)
    splitYear <- ceiling(n.years/2)
    
    tab1 <- rbind(colnames(prop[ ,1:splitYear]), prop[ ,1:splitYear])
    
    if(length(1:splitYear) == length((splitYear+1):n.years)){
      tab2 <-  rbind(colnames(prop[ ,(splitYear+1):n.years]),
                     prop[ ,(splitYear+1):n.years])
    } else {
      diffYears <- length((splitYear+1):n.years) - length(1:splitYear) 
      tab2 <-  rbind(c(colnames(prop[ ,(splitYear+1):n.years]), rep("", diffYears)),
                     cbind(prop[ ,(splitYear+1):n.years], matrix("", nrow = nrow(prop), ncol = diffYears)))
    }
  
  splitProp <- rbind(tab1, tab2)
  splitProp <- cbind(c("","F","M","Total","","F","M","Total"), splitProp)
  addtorow <- list()
  addtorow$pos <- list(2,4,6)
  addtorow$command <- c("\\rowcolor[gray]{.95}", 
                        "\\hline \\\\",
                        "\\rowcolor[gray]{.95}")
  print(xtable( splitProp,
                type = "latex",
                align = paste(c("l",rep("c", ncol(splitProp))), collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        include.colnames = FALSE,
        include.rownames = FALSE,
        hline.after = c(0,1,4,5, nrow(splitProp)),
        add.to.row = addtorow,
        file = file.path(working.dir, "tables/PropDetected.tex"))
  } else {
    
    ##-- Print  a single .tex 
    addtorow <- list()
    addtorow$pos <- list(2)
    addtorow$command <- c("\\rowcolor[gray]{.95}")
    print(xtable( prop,
                  type = "latex",
                  align = paste(c("l",rep("c", ncol(prop))), collapse = "")),
          floating = FALSE,
          sanitize.text.function = function(x){x},
          hline.after = c(-1,0, nrow(prop)),
          add.to.row = addtorow,
          file = file.path(working.dir, "tables/PropDetected.tex"))
  }
    

  
  # ## ------     5.2.5. NUMBER OF IDs w/ ACs OUTSIDE NORWAY ------
  # 
  # ##-- Calculate number of individuals alive with their AC in each country each year
  # N_det_by_country <- matrix(NA,5,n.years)
  # dimnames(N_det_by_country) <- list("Countries" = c("Norway","Sweden","Finland","Russia","Out"),
  #                                    "Years" = c(years))
  # for(t in 1:n.years){
  #   
  #   N_fin_F <- N_fin_M <- N_fin <- rep(NA,n.mcmc)
  #   N_nor_F <- N_nor_M <- N_nor <- rep(NA,n.mcmc)
  #   N_rus_F <- N_rus_M <- N_rus <- rep(NA,n.mcmc)
  #   N_swe_F <- N_swe_M <- N_swe <- rep(NA,n.mcmc)
  #   N_out_F <- N_out_M <- N_out <- rep(NA,n.mcmc)
  #   
  #   for(iter in 1:n.mcmc){
  #     
  #     country <- countryRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
  #     isFin <- country %in% 1
  #     isNor <- country %in% 2
  #     isRus <- country %in% 3
  #     isSwe <- country %in% 4
  #     
  #     ##-- Detected individuals
  #     N_fin_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isFin & isFemale)
  #     N_fin_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isFin & isMale)
  #     N_fin[iter] <- N_fin_F[iter] + N_fin_M[iter]
  #     
  #     N_nor_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isFemale)
  #     N_nor_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isMale)
  #     N_nor[iter] <- N_nor_F[iter] + N_nor_M[iter]
  #     
  #     N_rus_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isRus & isFemale)
  #     N_rus_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isRus & isMale)
  #     N_rus[iter] <- N_rus_F[iter] + N_rus_M[iter]
  #     
  #     N_swe_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isSwe & isFemale)
  #     N_swe_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isSwe & isMale)
  #     N_swe[iter] <- N_swe_F[iter] + N_swe_M[iter]
  #     
  #     N_out_F[iter] <- N_fin_F[iter] + N_rus_F[iter] + N_swe_F[iter]
  #     N_out_M[iter] <- N_fin_M[iter] + N_rus_M[iter] + N_swe_M[iter]
  #     N_out[iter] <- N_out_F[iter] + N_out_M[iter]
  #   }#iter
  #   
  #   N_det_by_country["Norway",t] <- getCleanEstimates(N_nor)
  #   N_det_by_country["Sweden",t] <- getCleanEstimates(N_swe)
  #   N_det_by_country["Finland",t] <- getCleanEstimates(N_fin)
  #   N_det_by_country["Russia",t] <- getCleanEstimates(N_rus)
  #   N_det_by_country["Out",t] <- getCleanEstimates(N_out)
  # 
  # }#t
  # 
  # ##-- print .csv
  # write.csv(N_det_by_country,
  #           file = file.path(working.dir, "tables/NumDetectedIds_country.csv"))
  #
  # ##-- Calculate number of individuals detected/undetected in Norway
  # N_NOR <- matrix(NA,4,n.years)
  # dimnames(N_NOR) <- list("#individuals" = c("Detected","Undetected","Total","%"),
  #                         "Years" = c(years))
  # for(t in 1:n.years){
  #   N_det <- N_undet <- N_tot <- rep(NA,n.mcmc)
  #   for(iter in 1:n.mcmc){
  #     country <- countryRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
  #     isNor <- country %in% 2
  #     
  #     ##-- Detected individuals
  #     N_det[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor)
  #     
  #     ##-- Undetected individuals
  #     N_undet[iter] <- sum(!isDetected[ ,t] & isAlive[iter, ,t] & isNor)
  #     
  #     ##-- Total Norway
  #     N_tot[iter] <- N_det[iter] + N_undet[iter]
  #   }#iter
  #   
  #   N_NOR["Detected",t] <- getCleanEstimates(N_det)
  #   N_NOR["Undetected",t] <- getCleanEstimates(N_undet)
  #   N_NOR["Total",t] <- getCleanEstimates(N_tot)
  #   
  #   print(t)
  # }#t
  # ##-- Print number of individuals with AC in Norway
  # print(N_NOR)
  # 
  # ##-- Remove unnecessary objects from memory
  # rm(list = c("N_fin_F", "N_fin_M", "N_fin",
  #             "N_nor_F", "N_nor_M", "N_nor",
  #             "N_rus_F", "N_rus_M", "N_rus",
  #             "N_swe_F", "N_swe_M", "N_swe",
  #             "N_out_F", "N_out_M", "N_out",
  #             "isDetected"))
  # gc(verbose = FALSE)  
  # 
  # 
  # 
  # 
  # ##----------------------------------------------------------------------------
  # ## ------   2.1. OVERALL NUMBERS ------
  # 
  # ##-- SOME TALLIES TO CHECK THINGS
  # ##-- NGS
  # NGS <- myFilteredData.sp$alive #rbind(myFilteredData.spF$alive, myFilteredData.spM$alive)
  # NGSStructured <- myData.aliveStruc$myData.sp #rbind( myFilteredData.spStructuredF, myFilteredData.spStructuredM)
  # NGSOther <- myData.aliveOthers$myData.sp #rbind( myFilteredData.spOthersF, myFilteredData.spOthersM)
  # 
  # ##-- FOR REPORT SUMMARY
  # length(NGS$Id)                                ## Number of NGS samples
  # length(NGS$Id[NGS$Sex == "Hunn"])             ## Number of Female NGS samples
  # length(NGS$Id[NGS$Sex == "Hann"])             ## Number of Male NGS samples
  # length(NGS$Id[NGS$Country == "S"])/nrow(NGS)  ## Proportion of samples in Sweden
  # 
  # length(unique(NGS$Id))                        ## Number of individuals
  # length(unique(NGS$Id[NGS$Sex == "Hunn"]))     ## Number of Female individuals
  # length(unique(NGS$Id[NGS$Sex == "Hann"]))     ## Number of Male individuals
  # 
  # ## Last year
  # length(NGS$Id[NGS$Year %in% tail(years, n = 1)])                     ## Number of NGS in the last year
  # length(NGS$Id[NGS$Sex == "Hunn" & NGS$Year %in% tail(years, n = 1)]) ## Number of female NGS in the last year
  # length(NGS$Id[NGS$Sex == "Hann" & NGS$Year %in% tail(years, n = 1)]) ## Number of Male NGS in the last year
  # 
  # ## NGS structured
  # length(NGSStructured$Id)                            ## Number of structured NGS samples
  # length(NGSStructured$Id[NGSStructured$Sex=="Hunn"]) ## Number of structured Female NGS samples
  # length(NGSStructured$Id[NGSStructured$Sex=="Hann"]) ## Number of structured Male NGS samples
  # 
  # ## NGS Other
  # length(NGSOther$Id)                           ## Number of other NGS samples
  # length(NGSOther$Id[NGSOther$Sex=="Hunn"])     ## Number of other Female NGS samples
  # length(NGSOther$Id[NGSOther$Sex=="Hann"])     ## Number of other Male NGS samples
  # 
  # 
  # ##-- DEAD RECOVERY
  # dead <- myFullData.sp$dead.recovery #rbind(myFullData.sp$dead.recovery, myFullData.spM$dead.recovery)
  # table(dead$Year)
  # length(dead$Id)                               ## Number of individuals
  # length(unique(dead$Id[dead$Sex=="Hunn"]))     ## Number of Female individuals
  # length(unique(dead$Id[dead$Sex=="Hann"]))     ## Number of Male individuals
  # 
  # # tmpdead <- dead[dead$Year %in% c(2018:2023),]
  # # tmpdead <- tmpdead[!duplicated(tmpdead$DNAID), ]
  # # table(tmpdead$Year,tmpdead$Sex)
  # # table(tmpdead$Year,tmpdead$Month)
  # # table(tmpdead$Year)
  # # nrow(table(tmpdead$Year))
  # # duplicated(tmpdead$DNAID)
  # # mapview::mapview(tmpdead)
  # # 
  # # tmpNGS <- NGS[NGS$Year %in% c(2018:2023), ]
  # # table(tmpNGS$Year,tmpNGS$Sex)
  # # table(tmpNGS$Year,tmpNGS$Month)
  # # table(tmpNGS$Year)
  # # 
  # # ###HENRIK CHECK WITH PUBLIC SAMPLES
  # # idPublic <- c(
  # # 'D555438',
  # # 'D556590',
  # # 'D553322',
  # # 'D555239',
  # # 'D554783',
  # # 'D558440',
  # # 'D555845',
  # # 'D555412',
  # # 'D556343',
  # # 'D556344',
  # # 'D556347',
  # # 'D558235',
  # # 'D557101')
  # # 
  # # which(NGSStructured$DNAID%in%idPublic)
  # # 
  # # # # NGSStructured[which(NGSStructured$DNAID%in%idPublic),]@data
  # # # # NGSOther[which(NGSOther$DNAID%in%idPublic),]@data
  # # # 
  # # # idPublic1 <- c('D555438')
  # # # NGSStructured[which(NGSStructured$DNAID%in%idPublic1),]@data
  # # # NGSOther[which(NGSOther$DNAID%in%idPublic1),]
  # # # 
  # # # mapview(
  # # # list(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"),
  # # #         NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]),
  # # # layer.name = c("Franconian districts", "Franconian breweries")
  # # # )
  # # #         
  # # # mapview(as(st_geometry(TRACKS_YEAR[[t]][TRACKS_YEAR[[t]]$RovbaseID %in% "T477952",]),"Spatial"))+
  # # #   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
  # # # 
  # # # mapview(as(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),"Spatial"))+
  # # #   NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]
  # # # 
  # # # st_distance(st_geometry(TRACKS[TRACKS$RovbaseID %in% "T477952",]),
  # # #             st_as_sf(NGSOther[which(NGSOther$DNAID%in%idPublic),][1,]))
  # # # 
  # # # plot(TRACKS_YEAR[[9]][TRACKS_YEAR[[9]]$RovbaseID %in% "T471191",]$geometry)
  # # # plot(TRACKSSimple_sf[[9]][TRACKSSimple_sf[[9]]$RovbaseID %in% "T471191",]$geometry)
  # # # plot(TRACKS[TRACKS$RovbaseID %in% "T471191",]$geometry,col="red",add=T)
  # # # points(NGSStructured[which(NGSStructured$DNAID%in%idPublic),],pch=16)
  # 
  # 
  # 
  # ## ------   2.2. TABLE 1 NGS SAMPLES YEAR/COUNTRIES/SEX ------
  # 
  # ## ------     2.2.1. ALL ------
  # 
  # NGSCountrySEX <- matrix("", ncol = n.years*2, nrow = 4)
  # row.names(NGSCountrySEX) <- c("","Norway","Sweden","Total")
  # colnames(NGSCountrySEX) <- unlist(lapply(YEARS, function(x) c(x[2],x[2])))
  # #colnames(NGSCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
  # 
  # NGSCountrySEX[1, ] <- rep(c("F","M"), n.years)
  # sex <- c("Hunn","Hann")
  # sex1 <- c(0,1)
  # ye <- seq(1, n.years*2, by = 2)
  # for(s in 1:2){
  #   for(t in 1:n.years){
  #     temp <- NGS[NGS$Year == years[t] & NGS$Sex==sex[s], ]
  #     NGSCountrySEX["Norway",ye[t] + sex1[s] ] <- nrow(temp[temp$Country %in% "N", ])
  #     NGSCountrySEX["Sweden",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S", ])
  #     NGSCountrySEX["Total",ye[t] + sex1[s]] <- nrow(temp[temp$Country %in% "S" | temp$Country %in% "N" , ])
  #   }#t
  # }
  # 
  # ##-- Export .tex table
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEX))),
  #                                     '}', collapse=''), '\\\\'),
  #                       rep("\\rowcolor[gray]{.95}",1))
  # colnames(NGSCountrySEX) <- rep("", ncol(NGSCountrySEX))
  # print(xtable(NGSCountrySEX, type = "latex",
  #              align = paste(c("l",rep("c",ncol(NGSCountrySEX))), collapse = "")),
  #       #scalebox = .8,
  #       floating = FALSE,include.colnames=F,
  #       add.to.row = addtorow,
  #       file = file.path(working.dir, "tables","NGSCountrySEX.tex"))
  # 
  # ##-- Export .csv table
  # write.csv( NGSCountrySEX, 
  #            file = file.path(working.dir, "tables", "NGSCountrySEX.csv"))
  # 
  # # ##-- Checks
  # # sum(as.numeric(NGSCountrySEX["Total",]))
  # # sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
  # # sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))
  # 
  # 
  # 
  # ## ------     2.2.2. PER OBSERVATION PROCESS ------
  # 
  # NGSCountrySEXoBS <- matrix("", ncol = n.years*2+1, nrow = 7)
  # NGSCountrySEXoBS[1, ] <- c("", rep(c("F","M"), n.years))
  # NGSCountrySEXoBS[ ,1] <- c("", rep(c("Structured","Unstructured"), 3))
  # row.names(NGSCountrySEXoBS) <- c("", rep(c("Norway","Sweden","Total"), each = 2))
  # colnames(NGSCountrySEXoBS) <- c("", unlist(lapply(YEARS, function(x) c(x[2],x[2]))))
  # 
  # sex <- c("Hunn","Hann")
  # sex1 <- c(0,1)
  # ye <- seq(2, n.years*2, by = 2)
  # for(s in 1:2){
  #   for(t in 1:n.years){
  #     ##-- Structured
  #     tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex == sex[s], ]
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[1],ye[t]+sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "N", ])
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[1],ye[t]+sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S", ])
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[1], ye[t]+sex1[s]] <- nrow(tempStruc[tempStruc$Country %in% "S" | tempStruc$Country %in% "N", ])
  #     
  #     ##-- Other
  #     tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex == sex[s], ]
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Norway")[2],ye[t]+sex1[s]] <- nrow(tempOther[tempOther$Country %in% "N", ])
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Sweden")[2],ye[t]+sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S", ])
  #     NGSCountrySEXoBS[which(row.names(NGSCountrySEXoBS) %in% "Total")[2], ye[t]+sex1[s]] <- nrow(tempOther[tempOther$Country %in% "S" | tempOther$Country %in% "N", ])
  #   }#t
  # }#s
  # 
  # ##-- Export .tex table
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # addtorow$command <- c(paste0('& \\multicolumn{2}{c}{}',
  #                              paste0(' & \\multicolumn{2}{c}{',
  #                                     sort(unique(colnames(NGSCountrySEXoBS)[2:ncol(NGSCountrySEXoBS)])),
  #                                     '}', collapse = ""),
  #                              '\\\\'),
  #                       rep("\\rowcolor[gray]{.95}",1))
  # colnames(NGSCountrySEXoBS) <- rep("", ncol(NGSCountrySEXoBS))
  # multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
  # multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""), ncol = 1)
  # NGSCountrySEXoBS <- data.frame(cbind(multirowadd,NGSCountrySEXoBS))
  # colnames(NGSCountrySEXoBS) <- c("",unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),
  #                                                                      paste(x,collapse = "/")))))
  # print(xtable( NGSCountrySEXoBS,
  #               type = "latex",
  #               align = paste(c("l",rep("c",ncol(NGSCountrySEXoBS))), collapse = "")),
  #       #scalebox = .7, 
  #       floating = FALSE,
  #       add.to.row = addtorow,
  #       include.colnames = F,
  #       include.rownames = FALSE,
  #       sanitize.text.function = function(x){x},
  #       file = file.path(working.dir, "tables","NGSCountrySEXperObs.tex"))
  # 
  # # ##-- Checks
  # # sum(as.numeric(NGSCountrySEX["Total",]))
  # # sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "F"]))
  # # sum(as.numeric(NGSCountrySEX["Total",NGSCountrySEX[1,]%in% "M"]))
  # # 
  # # #PLOT CHECK 
  # # plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
  # # plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="red",add=T)
  # # plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="red",add=T)
  # # 
  # # par(mfrow=c(1,2),mar=c(0,0,0,0))
  # # plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
  # # plot(st_geometry(NGSStructured),pch=21,col="black",cex=0.5,bg="#E69F00",add=T)
  # # plot(habitat$habitat.r,axes=F,legend=F,box=F,col=c(grey(0.99),grey(0.8)))
  # # plot(st_geometry(NGSOther),pch=21,col="black",cex=0.5,bg="#009E73",add=T)
  # # 
  # # dev.off()
  # # 
  # # NGSStructured$Year1 <- NGSStructured$Year+1
  # # NGSOther$Year1 <- NGSOther$Year+1
  # # 
  # # barplot(rbind(table(NGSStructured$Year1),table(NGSOther$Year1)),
  # #         col=c("#E69F00","#009E73"))
  # # legend("topleft",fill=c("#E69F00","#009E73"),legend=c("Structured","Other") )
  # 
  # ##-- GIVE FILE TO HENRIK ([PD] for what????)
  # tmp <- NGSOther[NGSOther$Year %in% c(2019,2020,2021), ]
  # #tmp
  # write.csv(tmp, file = file.path(working.dir, "tables", "Unstructured2020_2022.csv"))
  # 
  # 
  # 
  # ## ------   2.3. TABLE 2 NGS ID YEAR/COUNTRIES/SEX ------
  # 
  # ## ------     2.3.1. ALL ------
  # 
  # NGSidCountrySEX <- matrix("", ncol = n.years*2, nrow = 4)
  # row.names(NGSidCountrySEX) <- c("","Norway","Sweden","Total")
  # colnames(NGSidCountrySEX) <- c(unlist(lapply(YEARS, function(x)c(x[2],x[2]))))
  # #colnames(NGSidCountrySEX) <- unlist(lapply(YEARS,function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")) ))
  # 
  # NGSidCountrySEX[1,] <- rep(c("F","M"), n.years)
  # sex <- c("Hunn","Hann")
  # sex1 <- c(0,1)
  # ye <- seq(1, n.years*2, by = 2)
  # for(s in 1:2){
  #   for(t in 1:n.years){
  #     temp <- NGS[NGS$Year == years[t] & NGS$Sex == sex[s], ]
  #     NGSidCountrySEX["Norway",ye[t] + sex1[s] ] <- length(unique(temp$Id[temp$Country %in% "N"]))
  #     NGSidCountrySEX["Sweden",ye[t] + sex1[s]] <- length(unique(temp$Id[temp$Country %in% "S"]))
  #     NGSidCountrySEX["Total",ye[t] + sex1[s]] <- length(unique(temp$Id))
  #   }#t
  # }#s
  # 
  # ##-- Export .tex
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSidCountrySEX))),
  #                                     '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))
  # colnames(NGSidCountrySEX) <- rep("", ncol(NGSidCountrySEX))
  # 
  # print(xtable( NGSidCountrySEX,
  #               type = "latex",
  #               align = paste(c("l",rep("c",ncol(NGSidCountrySEX))), collapse = "")),
  #       #scalebox = .8, 
  #       floating = FALSE,
  #       include.colnames = F,
  #       add.to.row = addtorow,
  #       file = file.path(working.dir, "tables","NGSidCountrySEX.tex"))
  # 
  # ##-- Export .csv
  # write.csv( NGSidCountrySEX,
  #            file = file.path(working.dir, "tables", "NGSidCountrySEX.csv"))
  # 
  # 
  # ##-- Export .csv TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR
  # NGSidCountryTotal <- matrix(0, ncol = n.years, nrow = 1)
  # row.names(NGSidCountryTotal) <- c("Total")
  # colnames(NGSidCountryTotal) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))
  # #colnames(NGSidCountryTotal) <- unlist(lapply(YEARS,function(x) paste(x,collapse = "/") ))
  # 
  # for(t in 1:n.years){
  #   temp <- NGS[NGS$Year == years[t] , ]
  #   NGSidCountryTotal["Total", t] <- length(unique(temp$Id))
  # }#t
  # write.csv( NGSidCountryTotal,
  #            file = file.path(working.dir, "tables","TotalIdDetected.csv"))
  # 
  # ### PRINT A CSV TABLE WITH THE NUMBER OF TOTAL IDS PER YEAR PER SEX
  # 
  # 
  # 
  # ## ------     2.3.2. PER OBSERVATION PROCESS ------
  # 
  # NGSCountrySEXoBSid <- matrix("", ncol = n.years*2+1, nrow = 7)
  # row.names(NGSCountrySEXoBSid) <- c("",rep(c("Norway","Sweden","Total"),each=2))
  # colnames(NGSCountrySEXoBSid) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))
  # 
  # NGSCountrySEXoBSid[1,] <- c("", rep(c("F","M"), n.years))
  # NGSCountrySEXoBSid[,1] <- c("", rep(c("Structured","Unstructured"), 3))
  # 
  # sex <- c("Hunn","Hann")
  # sex1 <- c(0,1)
  # ye <- seq(2,n.years*2,by=2)
  # for(s in 1:2){
  #   for(t in 1:n.years){
  #     ## structured
  #     tempStruc <- NGSStructured[NGSStructured$Year == years[t] & NGSStructured$Sex==sex[s], ]
  #     
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[1], ye[t] + sex1[s] ] <- length(unique(tempStruc$Id[tempStruc$Country %in% "N" ])) 
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id[tempStruc$Country %in% "S" ]))
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[1], ye[t] + sex1[s]] <- length(unique(tempStruc$Id))
  #     
  #     ## Other
  #     tempOther <- NGSOther[NGSOther$Year == years[t] & NGSOther$Sex==sex[s], ]
  #     
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Norway")[2], ye[t] + sex1[s] ] <- length(unique(tempOther$Id[tempOther$Country %in% "N" ])) 
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Sweden")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id[tempOther$Country %in% "S" ]))
  #     NGSCountrySEXoBSid[which(row.names(NGSCountrySEXoBSid) %in% "Total")[2], ye[t] + sex1[s]] <- length(unique(tempOther$Id))
  #     
  #     ###TOTAL 
  #     
  #     
  #   }#t
  # }#s
  # 
  # 
  # ##-- Export .tex
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # addtorow$command <- c(paste0("& \\multicolumn{1}{c}{}",paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(NGSCountrySEXoBSid)[2:ncol(NGSCountrySEXoBSid)])),
  #                                                               '}', collapse=''), '\\\\'),
  #                       rep("\\rowcolor[gray]{.95}",1))
  # colnames(NGSCountrySEXoBSid) <- rep("", ncol(NGSCountrySEXoBSid))
  # multirow <- paste0( paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Norway","Sweden","Total"), "}}"))
  # multirowadd <- matrix(c("",multirow[1],"",multirow[2],"",multirow[3],""),ncol=1)
  # NGSCountrySEXoBSid <- data.frame(cbind(multirowadd, NGSCountrySEXoBSid))
  # colnames(NGSCountrySEXoBSid) <- c(unlist(lapply(YEARS, function(x) c(x[2]))))
  # 
  # print(xtable(NGSCountrySEXoBSid, type = "latex",
  #              align = paste(c("l",rep("c",ncol(NGSCountrySEXoBSid))),collapse = "")),
  #       #scalebox = .7, 
  #       floating = FALSE,
  #       add.to.row = addtorow,
  #       include.colnames = F,
  #       include.rownames = FALSE,
  #       sanitize.text.function = function(x){x},
  #       file = file.path(working.dir, "tables", "NGSCountrySEXperObsid.tex"))
  # 
  # 
  # 
  # ## ------   2.4. TABLE 3 DEAD CAUSE ID YEAR/COUNTRIES/SEX ------
  # 
  # DeadidCountrySEX <- matrix(0, ncol = n.years*2+1, nrow = 6)
  # row.names(DeadidCountrySEX) <- c("","other","other","legal culling","legal culling","")
  # colnames(DeadidCountrySEX) <- c("",unlist(lapply(YEARS, function(x) c(x[2],x[2]))))
  # #colnames(DeadidCountrySEX) <- c("", unlist(lapply(YEARS, function(x) c(paste(x,collapse = "/"),paste(x,collapse = "/")))))
  # DeadidCountrySEX[1,] <- c("",rep(c("F","M"),n.years))
  # DeadidCountrySEX[,1] <- c("","Norway","Sweden","Norway","Sweden","Total")
  # sex <- c("Hunn","Hann")
  # sex1 <- c(0,1)
  # ye <- seq(1, n.years*2, by = 2)
  # 
  # ##-- Identify legal mortality causes
  # MortalityNames <- unique(as.character(dead$DeathCause))
  # table(as.character(dead$DeathCause))
  # legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
  # legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
  # legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("9", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("23", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("28", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("Rifle", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("18", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("17", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Uspesifisert", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Fellefangst", MortalityNames)])
  # # legalCauses <- c(legalCauses, MortalityNames[grep("Jakt - Hagle", MortalityNames)])
  # 
  # 
  # ##-- SEPARATE MORTALITIES
  # cause <- c("other","legal culling")
  # for(t in 1:n.years){
  #   for(s in 1:2){
  #     for(d in 1:2){
  #       if(d == 1){
  #         temp <- dead[dead$Year == years[t] & dead$Sex == sex[s] & !(dead$DeathCause %in% legalCauses), ]
  #       } else {
  #         temp <- dead[dead$Year == years[t] & dead$Sex == sex[s] & dead$DeathCause %in% legalCauses, ]
  #       }
  #       row <- which(rownames(DeadidCountrySEX) == cause[d] & DeadidCountrySEX[,1] == "Norway" )
  #       DeadidCountrySEX[row,ye[t]+sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "N"]))
  #       
  #       row <- which(rownames(DeadidCountrySEX) == cause[d] & DeadidCountrySEX[,1] == "Sweden" )
  #       DeadidCountrySEX[row,ye[t]+sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "S"]))
  #     }#d
  #     DeadidCountrySEX[6,ye[t]+sex1[s]+1] <- sum(as.numeric(DeadidCountrySEX[2:6,ye[t]+sex1[s]+1]))
  #   }#s
  # }#t
  # 
  # 
  # ##-- summary
  # ###-- Other causes
  # sum(as.numeric(DeadidCountrySEX[2:3,2:ncol(DeadidCountrySEX)]))
  # sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="F")]))
  # sum(as.numeric(DeadidCountrySEX[2:3,which(DeadidCountrySEX[1,]=="M")]))
  # ##-- legal
  # sum(as.numeric(DeadidCountrySEX[4:5,2:ncol(DeadidCountrySEX)]))
  # sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="F")]))
  # sum(as.numeric(DeadidCountrySEX[4:5,which(DeadidCountrySEX[1,]=="M")]))
  # sum(as.numeric(DeadidCountrySEX[c(2,3),2:ncol(DeadidCountrySEX)]))/
  #   sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))
  # 
  # ##-- %of dead reco (legal) in Norway
  # sum(as.numeric(DeadidCountrySEX[4,2:ncol(DeadidCountrySEX)]))/
  #   sum(as.numeric(DeadidCountrySEX[c(4,5),2:ncol(DeadidCountrySEX)]))
  # 
  # ##-- ?? 
  # sum(as.numeric(DeadidCountrySEX[c(3,5),2:ncol(DeadidCountrySEX)]))/
  #   sum(as.numeric(DeadidCountrySEX[c(2:5),2:ncol(DeadidCountrySEX)]))
  # 
  # ##-- ??
  # sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="M")]))
  # sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,]=="F")]))
  # sum(as.numeric(DeadidCountrySEX[6,which(DeadidCountrySEX[1,] %in% c("F","M"))]))
  # 
  # ##-- Export. tex
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # uniqueYEAR <- sort(unique(colnames(DeadidCountrySEX)))
  # uniqueYEAR <- uniqueYEAR[2:length(uniqueYEAR)]
  # addtorow$command <- c(paste0("& \\multicolumn{1}{c}{Country}",paste0('& \\multicolumn{2}{c}{', uniqueYEAR,
  #                                                                      '}', collapse=''), '\\\\'),
  #                       rep("\\rowcolor[gray]{.95}",1))
  # ##-- REMOVE ROWS WHERE PARAMETERS ARE NOT STATE SPECIFIC
  # multirow <- paste0("\\multirow{", 2, "}{*}{\\textbf{", c("Other","Legal culling"), "}}")
  # multirowadd <- matrix(c("",multirow[1],"",multirow[2],"","{\\textbf{Total}}"),ncol=1)
  # DeadidCountrySEX <- data.frame(cbind(multirowadd,DeadidCountrySEX))
  # 
  # print(xtable(DeadidCountrySEX, type = "latex",
  #              align = paste(rep("c", ncol(DeadidCountrySEX)+1), collapse = "")),
  #       #scalebox = .7, 
  #       floating = FALSE,
  #       add.to.row = addtorow,
  #       include.colnames = F,
  #       include.rownames = FALSE,
  #       sanitize.text.function = function(x){x},
  #       file = file.path(working.dir, "tables", "DeadidCountrySEX.tex"))
  # 
  # # check 
  # # tmp <- dead[dead$Year == 2019 & 
  # #               !dead$DeathCause %in% legalCauses &
  # #               dead$Country %in% "N" & 
  # #               dead$Sex %in% "Hunn", ]
  # # length(unique(tmp$Id))
  # # unique(tmp$Id)
  # # DeadidCountrySEX[,"X2022"]
  # # DeadidCountrySEX[,"X2022.1"]
  # # dead[dead$Id %in% "JI416817 Ind7303 +",]$Sex
  # # plot(COUNTRIES$geometry)
  # # plot(tmp$geometry,add=T,col="red",pch=16)
  # 
  # 
  # 
  # ## ------   2.5. GET THE DETECTED INDIVIDUALS ------
  # 
  # n.detected <- read.csv(file.path(working.dir, "tables", "TotalIdDetected.csv"))
  # n.detected <- n.detected[1,2:ncol(n.detected)]
  # 
  # 
  # 
  # ## ------   2.6. SUMMARY DETECTED INDIVIDUALS PER COUNTIES ------
  # 
  # # myFilteredData.sp$alive$COUNTIES  <- st_intersects(myFilteredData.sp$alive[,1], COUNTIES_AGGREGATED[,1])
  # # myFilteredData.sp$alive$COUNTIES <- as.numeric(myFilteredData.sp$alive$COUNTIES)
  # # 
  # # 
  # # myFilteredData.sp$alive$COUNTIES <- apply(st_intersects(COUNTIES_AGGREGATED, myFilteredData.sp$alive, sparse = FALSE), 2, 
  # #       function(col) {which(col)})
  # # myFilteredData.sp$alive$counties1 <- 0
  # # for(i in 1:nrow(myFilteredData.sp$alive)){
  # #   if(length(myFilteredData.sp$alive$COUNTIES[[i]])>0){
  # #   myFilteredData.sp$alive$counties1[i] <- myFilteredData.sp$alive$COUNTIES[[i]][1]
  # #   }else{
  # #     myFilteredData.sp$alive$counties1[[i]] <- 0
  # #   }
  # # }
  # # 
  # # par(mar=c(0,0,0,0))
  # # plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
  # # text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
  # #      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
  # #      COUNTIES_AGGREGATED$id,col="red",font=2)
  # #
  # # ## NSAMPLES
  # # summa <- myFilteredData.sp$alive %>%
  # #   group_by(counties1,Year) %>%
  # #   summarise(n=n()) %>%
  # #   st_drop_geometry()
  # # summa <- summa[summa$counties1>0,]
  # # ##-- 2023
  # # tmp <- summa[summa$Year %in% 2023,]
  # # COUNTIES_AGGREGATED$nSampl2023 <- tmp$n#[1:8,"n"]
  # # #2022
  # # tmp1 <- summa[summa$Year %in% 2022,]
  # # COUNTIES_AGGREGATED$nSampl2022 <- tmp1$n#[1:8,"n"]
  # # 
  # # ##
  # # tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("id","nSampl2022","nSampl2023")])
  # # bar <- t(as.matrix(tmppp[c(5,7,8),2:3]))
  # # colnames(bar) <- c(5,7,8)
  # # barplo <- barplot(bar,beside=T,ylab="N samples")
  # # legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))
  # # 
  # # ## NdetectionsPerID 
  # # #COUNT NUMBER IDS 
  # # summa$n1 <- summa$NID <- 0
  # # dpt <- unique(unlist(summa$counties1))
  # # yearsss <- c(2022,2023)
  # # for(t in 1:length(yearsss)){
  # #   for(i in 1:length(dpt)){
  # #     tmp <- myFilteredData.sp$alive[myFilteredData.sp$alive$counties1 %in% dpt[i] & myFilteredData.sp$alive$Year %in% yearsss[t],]
  # #     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$n1 <- nrow(tmp)
  # #     summa[summa$counties1 %in% dpt[i] & summa$Year %in% yearsss[t], ]$NID <-  length(unique(tmp$Id))
  # #   }
  # # }
  # # 
  # # summa$NdetPerIDDet <- summa$n1/summa$NID
  # # 
  # # ##-- 2023
  # # tmp <- summa[summa$Year %in% 2023,]
  # # COUNTIES_AGGREGATED$detPerID2023 <- tmp$NdetPerIDDet#[1:8,"n"]
  # # #2022
  # # tmp1 <- summa[summa$Year %in% 2022,]
  # # COUNTIES_AGGREGATED$detPerID2022 <- tmp1$NdetPerIDDet#[1:8,"n"]
  # # 
  # # tmppp<- st_drop_geometry(COUNTIES_AGGREGATED[,c("detPerID2022","detPerID2023")])
  # # bar <- t(as.matrix(tmppp[c(5,7,8),]))
  # # colnames(bar) <- c(5,7,8)
  # # 
  # # par(mfrow=c(1,2))
  # # par(mar=c(0,0,0,0))
  # # plot(COUNTIES_AGGREGATED$geometry,border="white",col=grey(0.5))
  # # text(st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,1],
  # #      st_coordinates(st_centroid(COUNTIES_AGGREGATED$geometry))[,2],
  # #      COUNTIES_AGGREGATED$id,col="red",font=2)
  # # par(mar=c(4,5,1,1))
  # # barplo <- barplot(bar,beside=T,ylab="average Dets per IDS")
  # # legend("topright",fill=c(grey(0.3),grey(0.6)),legend=c(2023,2024))
  # 
  # 
  # ##--------------------------------------------------------------------------
  # 
  # 
  # ## ------   5.3. VITAL RATES ------
  # 
  # parameters <- c("rho","phi", "h", "w", "r")
  # sex <- c("F", "M")
  # vitalRate <- matrix(NA, nrow = length(parameters)+1, ncol = (n.years)*2-2)
  # rownames(vitalRate) <- c("", unlist(lapply(as.list(parameters), function(x)rep(x,1))))
  # colnames(vitalRate) <- c(unlist(lapply(years[1:(length(years)-1)], function(x)rep(paste(x,x+1,sep=" to "),2))))
  # vitalRate[1, ] <- c(rep(sex,(n.years-1)) )
  # 
  # for(s in 1:2){
  #   if(s == 1){results <- results_F} else {results <- results_M}
  # 
  #   col <- which(vitalRate[1, ] == sex[s])
  # 
  #   ##-- Per capita recruitment
  #   if(any(grep("rho",names(results$sims.list)))){
  #     vitalRate["rho",col] <- getCleanEstimates(results$sims.list$rho, moment = "median")
  #   } else {
  #     if(s == 1){
  #       for(t in 1:(n.years-1)){
  #         n.recruits <- rowSums(isAvail[ ,isFemale,t] * isAlive[ ,isFemale,t+1])
  #         alivetminus1 <- rowSums(isAlive[ ,isFemale,t])
  #         vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
  #       }#t
  #     } else {
  #       for(t in 1:(n.years-1)){
  #         n.recruits <- rowSums(isAvail[ ,isMale,t] * isAlive[ ,isMale,t+1])
  #         alivetminus1 <- rowSums(isAlive[ ,isMale,t])
  #         vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
  #       }#t
  #     }
  #   }
  # 
  # 
  #   ##-- Mortality & Survival
  #   if(any(grep("mhH",names(results$sims.list)))){
  #     ##-- Calculate mortality from estimated hazard rates (mhH and mhW)
  #     mhH1 <- exp(results$sims.list$mhH[ ,3:11])
  #     mhW1 <- exp(results$sims.list$mhW[ ,3:11])
  #     h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
  #     w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
  #     phi <- 1-h-w
  #     vitalRate["phi",col] <- apply(phi, 2, function(x) getCleanEstimates(x, moment = "median"))
  #     vitalRate["h",col] <- apply(h, 2, function(x) getCleanEstimates(x, moment = "median"))
  #     vitalRate["w",col] <- apply(w, 2, function(x) getCleanEstimates(x, moment = "median"))
  #   } else {
  # 
  #     ##-- Extract survival from posteriors
  #     vitalRate["phi",col] <- apply(results$sims.list$phi, 2, function(x) getCleanEstimates(x, moment = "median"))
  #     vitalRate["r",col] <- apply(results$sims.list$r, 2, function(x) getCleanEstimates(x, moment = "median"))
  # 
  #     if("h" %in% names(results$sims.list)) {
  #       vitalRate["h",col] <- apply(results$sims.list$h, 2, function(x) getCleanEstimates(x, moment = "median"))
  #       vitalRate["w",col] <- apply(results$sims.list$w, 2, function(x) getCleanEstimates(x, moment = "median"))
  #     } else {
  #       if(s == 1){
  #         y.dead <- nimDataF$y.dead
  #         z <- resultsSXYZ_MF$sims.list$z[ ,isFemale, ]
  #       } else {
  #         y.dead <- nimDataM$y.dead
  #         z <- resultsSXYZ_MF$sims.list$z[ ,isMale, ]
  #       }
  #       ##-- Derive mortality from posterior z and dead recoveries
  #       isDead <- apply((z[ , ,1:(n.years-1)] == 2)*(z[ , ,2:n.years] == 3), c(1,3), sum)
  #       wasAlive <- apply(z[ , ,1:(n.years-1)] == 2, c(1,3), sum)
  #       mortality <- isDead / wasAlive
  #       h <- sapply(1:(n.years-1), function(t)sum(y.dead[ ,t+1])/wasAlive[ ,t])
  #       w <- mortality - h
  #       vitalRate["h",col] <- apply(h, 2, function(x) getCleanEstimates(x, moment = "median"))
  #       vitalRate["w",col] <- apply(w, 2, function(x) getCleanEstimates(x, moment = "median"))
  #     }#else
  #   }#else
  # }#s
  # 
  # ##-- Print .csv
  # write.csv( vitalRate,
  #            file = file.path(working.dir, "tables/VitalRates.csv"))
  # 
  # ##-- Print .tex (print two separate tables to put in the overleaf document)
  # #colnames(vitalRate) <- rep("", ncol)
  # rownames(vitalRate)[2:6] <- c("$\\rho$","$\\phi$","h","w", "r")
  # 
  # 
  # if((n.years-1) > 8){
  #   ##-- Print .tex (split in two tables to print in the overleaf document)
  #   splitYear <- ceiling(n.years/2)
  # 
  #   ##-- Table Vital Rates #1
  #   addtorow <- list()
  #   addtorow$pos <- list(0,0)
  #   addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(vitalRate)))[1:splitYear],'}', collapse = ''), '\\\\'),
  #                         "\\rowcolor[gray]{.95}")
  #   print(xtable(vitalRate[,1:(splitYear*2)], type = "latex",
  #                align = paste(rep("c", (splitYear*2)+1), collapse = "")),
  #         floating = FALSE,
  #         add.to.row = addtorow,
  #         include.colnames = F,
  #         sanitize.text.function = function(x){x},
  #         file = file.path(working.dir, "tables/VitalRates_1.tex"))
  # 
  #   ##-- Table Vital Rates #2
  #   addtorow <- list()
  #   addtorow$pos <- list(0,0)
  #   addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(vitalRate)))[(splitYear+1):(n.years-1)],'}', collapse = ''), '\\\\'),
  #                         "\\rowcolor[gray]{.95}")
  # 
  #   print(xtable(vitalRate[ ,(splitYear*2+1):(2*(n.years-1))], type = "latex",
  #                align = paste(rep("c", length((splitYear*2+1):(2*(n.years-1)))+1), collapse = "")),
  #         floating = FALSE,
  #         add.to.row = addtorow,
  #         include.colnames = F,
  #         sanitize.text.function = function(x){x},
  #         file = file.path(working.dir, "tables/VitalRates_2.tex"))
  # } else {
  #   ##-- Table Vital Rates
  #   addtorow <- list()
  #   addtorow$pos <- list(0,0)
  #   addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{', sort(unique(colnames(vitalRate))),'}', collapse = ''), '\\\\'),
  #                         "\\rowcolor[gray]{.95}")
  # 
  #   print(xtable(vitalRate, type = "latex",
  #                align = paste(rep("c", ncol(vitalRate)+1), collapse = "")),
  #         floating = FALSE,
  #         add.to.row = addtorow,
  #         include.colnames = F,
  #         sanitize.text.function = function(x){x},
  #         file = file.path(working.dir, "tables/VitalRates.tex"))
  # 
  # }
  # 
  # ##-- Remove unnecessary objects from memory
  # rm(list = c("isMale", "isFemale","isAlive", "isAvail"))
  # gc(verbose = FALSE)
  # 
  # 
  # 
  # ## ------   5.4. DERIVED PARAMETERS FROM ABUNDANCE ------
  # 
  # ## ------     5.4.1. DERIVE SEX-RATIO ------
  # 
  # ##-- REGION-SPECIFIC PROPORTION OF FEMALES
  # PropFemale_regions <- list()
  # for(t in 1:n.years){
  #   PropFemale_regions[[t]] <- ACdensityF[[t]]$PosteriorRegions/
  #     (ACdensityM[[t]]$PosteriorRegions +
  #        ACdensityF[[t]]$PosteriorRegions)
  #   #print(rowMeans(PropFemale_regions[[t]], na.rm = T))
  # }#t
  # 
  # ##-- OVERALL PROPORTION OF FEMALES
  # PropFemale <- list()
  # for(t in 1:n.years){
  #   PropFemale[[t]] <- colSums(ACdensityF[[t]]$PosteriorAllRegions)/
  #     (colSums(ACdensityM[[t]]$PosteriorAllRegions) +
  #        colSums(ACdensityF[[t]]$PosteriorAllRegions))
  # }#t
  # 
  # ##-- Format table
  # propFemale_tab <- matrix(0, ncol = n.years, nrow = length(idregionTable))
  # row.names(propFemale_tab) <- idregionTable
  # colnames(propFemale_tab) <- years
  # for(t in 1:n.years){
  #   for(c in 1:7){
  #     propFemale_tab[idregionTable[c],t] <- getCleanEstimates(na.omit(PropFemale_regions[[t]][idregionTable[c], ]))
  #   }#c
  #   propFemale_tab[8,t] <- getCleanEstimates(PropFemale[[t]])
  # }#t
  # 
  # ##-- print .tex
  # row.names(propFemale_tab) <- c(paste0("\\hspace{0.1cm} ", idregionNOR), "TOTAL")
  # print(xtable( propFemale_tab,
  #               type = "latex",
  #               align = paste(c("l",rep("c",ncol(propFemale_tab))),collapse = "")),
  #       floating = FALSE,
  #       sanitize.text.function = function(x){x},
  #       add.to.row = list(list(seq(1, nrow(propFemale_tab), by = 2)), "\\rowcolor[gray]{.96} "),
  #       file = file.path(working.dir, "tables/propFemale.tex"))
  # 
  # 
  # 
  # 
  # ## ------     5.4.2. DERIVE AVERAGE DENSITY ------
  # 
  # ##-- Format table
  # averageDensity <- matrix(0, ncol = n.years, nrow = length(idregionTable))
  # row.names(averageDensity) <- idregionTable
  # colnames(averageDensity) <- years
  # 
  # ##-- Extract region ID levels
  # regionsLevels <- as.data.frame(raster::levels(rrRegions$Regions[[1]]))
  # 
  # ##-- Fill in table
  # for(c in 1:7){
  #   thisRegion <- regionsLevels$ID[regionsLevels$Regions == idregionTable[c]]
  #   thisArea <- sum(na.omit(rrRegions[ ] == thisRegion)) * 25
  #   for(t in 1:n.years){
  #     tmp <- format(round(100*(ACdensity[[t]]$summary[idregionTable[c],c("mean","95%CILow","95%CIHigh")]/thisArea), digits = 3), digits = 3)
  #     averageDensity[idregionTable[c],t] <- paste0(tmp[1], " (", tmp[2], "-", tmp[3], ")")
  #   }#t
  # }#c
  # 
  # ##-- Add total
  # areaSqKm2 <- sum(na.omit(rrRegions[ ] > 0)) * 25
  # for(t in 1:n.years){
  #   tmp <- format(round(100*(ACdensity[[t]]$summary["Total",c("mean","95%CILow","95%CIHigh")]/areaSqKm2), digits = 3), digits = 3)
  #   averageDensity["Total",t] <- paste0(tmp[1], " (", tmp[2], "-", tmp[3], ")")
  # }#t
  # 
  # ##-- Print .csv
  # write.csv( averageDensity,
  #            file = file.path(working.dir, "tables/AverageDensity.csv"))
  # 
  # 
  # 
  # ## ------     5.4.3. GROWTH RATE ------
  # 
  # growthRate <- list()
  # for(t in 1:(n.years-1)){
  #   growthRate[[t]] <- colSums(ACdensity[[t+1]]$PosteriorAllRegions)/
  #     colSums(ACdensity[[t]]$PosteriorAllRegions)
  # }#t
  # 
  # ##-- Put in a table format
  # growthRate_tab <- matrix(0, ncol = (n.years-1), nrow = 1)
  # colnames(growthRate_tab) <- paste(years[-n.years], years[-1], sep = " to ")
  # for(t in 1:(n.years-1)){
  #   growthRate_tab[1,t] <- getCleanEstimates(growthRate[[t]])
  # }#t
  # 
  # ##-- Print .tex
  # addtorow <- list()
  # addtorow$pos <- list(c(0),0)
  # addtorow$command <- c(paste0(paste('& {', sort(unique(colnames(growthRate_tab))),
  #                                    '}', collapse = ''), '\\\\'), rep("\\rowcolor[gray]{.95}",1))
  # colnames(growthRate_tab) <- rep("", ncol(growthRate_tab))
  # rownames(growthRate_tab) <- c("$\\lambda$")
  # 
  # print(xtable(growthRate_tab, type = "latex",
  #              align = paste(c("l", rep("c",ncol(growthRate_tab))), collapse = "")),
  #       floating = FALSE,
  #       add.to.row = addtorow,
  #       include.colnames = F,
  #       sanitize.text.function = function(x){x},
  #       file = file.path(working.dir, "tables/GrowthRates.tex"))
  # 
  # 
  # 
  # 
  # ## ------   5.5. TABLE OTHERS ------
  # 
  # parameters <- c("tau",
  #                 "betaDead","betaDens",
  #                 "betaDead","betaDens",
  #                 "sigma",
  #                 "betaDet","betaDet")
  # sex <- c("F","M")
  # TableOthers <- matrix(NA, nrow = length(parameters), ncol = 3)
  # rownames(TableOthers) <- parameters
  # colnames(TableOthers) <- c("", sex)
  # TableOthers[ ,1] <- c("$\\tau$",
  #                       "$\\beta_{dead_1}$","$\\beta_{skandobs_1}$",
  #                       "$\\beta_{dead_2}$","$\\beta_{skandobs_2}$",
  #                       "$\\sigma$",
  #                       "$\\beta_{roads}$","$\\beta_{obs}$")
  # 
  # for(s in 1:2){
  #   if(s == 1){results <- results_F} else {results <- results_M}
  #   TableOthers["tau",sex[s]] <- getCleanEstimates(results$sims.list$tau/1000, moment = "median")
  #   TableOthers[which(parameters == "betaDead"),sex[s]] <- apply(results$sims.list$betaDens[,1,], 2,
  #                                                                function(x) getCleanEstimates(x,moment = "median"))
  #   TableOthers[which(parameters == "betaDens"),sex[s]] <- apply(results$sims.list$betaDens[,2,],
  #                                                                2, function(x) getCleanEstimates(x,moment = "median"))
  #   TableOthers["sigma",sex[s]] <- getCleanEstimates(results$sims.list$sigma/1000, moment = "median")
  #   TableOthers[which(parameters == "betaDet"),sex[s]] <- apply(results$sims.list$betaDet, 2,
  #                                                               function(x) getCleanEstimates(x,moment = "median"))
  # }#s
  # ##-- Deal with negative values
  # TableOthers <- gsub("--", "-(-)", TableOthers)
  # 
  # ##-- Change row names
  # row.names(TableOthers) <- c("Spatial process","","","","",
  #                             "Detection Process","","")
  # 
  # ##-- Print .tex
  # multirow <- c("\\multirow{5}{*}{\\textbf{Spatial process}}","\\multirow{3}{*}{\\textbf{Detection process}}")
  # multirowadd <- matrix(c(multirow[1],"","","","",multirow[2],"",""), ncol = 1)
  # TableOthers <- data.frame(cbind(multirowadd, TableOthers))
  # 
  # addtorow <- list()
  # addtorow$pos <- list(0,5)
  # addtorow$command <- c(paste0("& {\\textbf{Parameters}} ",
  #                              paste0('& {\\textbf{',  sex, '}}', collapse = ''), '\\\\'),
  #                       "\\hline")
  # 
  # print(xtable(TableOthers, type = "latex",
  #              align = paste(c("ll", rep("c",ncol(TableOthers)-1)), collapse = "")),
  #       sanitize.text.function = function(x){x},
  #       floating = FALSE,
  #       include.rownames = FALSE,
  #       include.colnames = FALSE,
  #       add.to.row = addtorow,
  #       file = file.path(working.dir, "tables/TableParametersOthers.tex"))
  # 
  # 
  # 

  ## ------ 6. OUTPUT -----
  out$YEARS <- years
  
  return(out)
}
