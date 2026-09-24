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
  ##-- Years
  years = NULL,
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
    dplyr::mutate(id = dplyr::case_when(
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
    paste(years, years+1, sep = "-")
    
    
    
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
    paste(years, years+1, sep = "-")
    
    
    
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
    paste(years, years+1, sep = "-")
    

    
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
    names(UDdensity) <- paste(years, years+1, sep = "-")
    
    
    
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
  
  seasons <- paste(years, "/", substr(years+1, 3, 4), sep = "")
  intervals <- paste(years[-length(years)], years[-1], sep = "\n to \n")
  
  
  
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
  graphics::axis(1, at = c(1:(n.years)), labels = seasons, cex.axis = 1.5, padj = -1)
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
  graphics::axis(1, at = c(1:(n.years)), labels = seasons, cex.axis = 1.5, padj = -1)
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
  graphics::axis(1, at = c(1:(n.years)), labels = seasons, cex.axis = 1.5, padj = -1)
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
  
  
  
  ## ------   4.4. NGS, Dead recoveries & Carnivore obs ------
  
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
    mtext(text = seasons[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
    
    if(t == n.years){
      ##-- LEGEND
      xLeg <- 1000000
      yLeg <- 6350000
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
  colnames(NCarRegionEstimates) <- seasons
  
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
             file = file.path(working.dir, "tables/NAllYears.csv"))
             #fileEncoding = "latin1")
  
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
             file = file.path(working.dir, "tables/N_LastYearPerSex_region.csv"))
             #,fileEncoding = "latin1")
  
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
             file = file.path(working.dir, "tables", "NLastYearPerSexArea.csv"))
             #,fileEncoding = "latin1")
  
  
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
             file = file.path(working.dir, "tables/NAllYearsPerSex.csv"))
             #fileEncoding = "latin1")
  
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
             file = file.path(working.dir, "tables/NAllYearsNorwegianCounties.csv"))
             #fileEncoding = "latin1")
  
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
  n.detected_F <- apply(nimDataF$detNums + nimDataF$detNumsOth, 2, function(x)sum(x>0))
  n.detected_M <- apply(nimDataM$detNums + nimDataM$detNumsOth, 2, function(x)sum(x>0))

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
  write.csv( propDetected,
             file = file.path(working.dir, "tables/PropDetectedIds.csv"))
  
  ##-- print .tex
  print(xtable( propDetected,
                type = "latex",
                align = paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list( list(seq(1,nrow(propDetected), by = 2)),
                           "\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/PropDetectedIds.tex"))
  
  
  
  ## ------     5.2.4. PROPORTION OF THE POPULATION DETECTED ------
  
  ##-- Extract number of individuals detected
  isDetected <- rbind(nimDataM$detNums + nimDataM$detNumsOth,
                      nimDataF$detNums + nimDataF$detNumsOth) > 0
  
  ##-- Identify individual sex
  isFemale <- resultsSXYZ_MF$sims.list$sex == "F"
  isMale <- resultsSXYZ_MF$sims.list$sex == "M"
  
  ##-- Identify individual status
  #isAvail <- resultsSXYZ_MF$sims.list$z == 1 
  isAlive <- resultsSXYZ_MF$sims.list$z == 2 
  
  ##-- Calculate % of the wolverine population detected 
  prop <- matrix(NA,3,n.years)
  dimnames(prop) <- list("% individuals" = c("F","M","Total"),
                         "Years" = c(years+1))
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
    


  ## ------ 6. OUTPUT -----
  out$YEARS <- years+1
  
  return(out)
  
}
