#' @title RovQuant OPSCR bear output processing
#' 
#' @description
#' \code{processRovquantOutput_bear} calls a custom Rmarkdown template that combines 
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
#' @rdname processRovquantOutput_bear
#' @export
processRovquantOutput_bear <- function(
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
  
  ##-- Males
  load(list.files(file.path(working.dir, "nimbleInFiles/male"), full.names = T)[1])
  nimDataM <- nimData
  nimInitsM <- nimInits
  
  ##-- Remove unnecessary objects from memory
  rm(list = c("nimInits", "nimData"))
  gc(verbose = FALSE)
  
  ##-- Habitat
  load(file.path( working.dir, "data", paste0("Habitat_bear_", DATE, ".RData")))
  
  ##-- Detectors
  load(file.path( working.dir, "data", paste0("Detectors_bear_", DATE, ".RData")))
  
  ##-- Load filtered data
  load(file.path( working.dir, "data", paste0("FilteredData_bear_", DATE, ".RData")))
  
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
  years <- as.numeric(dimnames(detectors$covariates)[[3]])
  n.years <- length(years) 
  
  ##-- MERGE & SIMPLIFY SOME NORWEGIAN COUNTIES
  COUNTIES_s <- COUNTIES[COUNTIES$country %in% "NOR", ] %>% sf::st_intersection(COUNTRIES)
  COUNTIES_s$county[COUNTIES_s$county %in% c("Trøndelag","Nordland")] <- "Trøndelag"
  COUNTIES_s$county[COUNTIES_s$county %in% c("Troms","Finnmark")] <- "Finnmark"
  COUNTIES_s$county[!COUNTIES_s$county %in% c("Finnmark","Trøndelag")] <- "Innlandet"
  COUNTIES_s <- COUNTIES_s %>%
    group_by(county) %>%
    dplyr::summarize() 
 
  COUNTIES_s <- sf::st_simplify(sf::st_as_sf(COUNTIES_s), preserveTopology = T, dTolerance = 500)
  COUNTIES_s$index <- c(1,3,2)
  COUNTIES_s$Name <- c("NO1","NO3","NO2")
  
  
  
  ## ------ 2. PROCESS MCMC SAMPLES -----
  
  message("## Processing model MCMC outputs...")
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  mcmcTest <- TRUE
  if (!overwrite) {
    fileName <- paste0("MCMC_bear_", DATE, ".RData")
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
    
    ##-- RESCALE sigma AND tau TO THE ORIGINAL COORDINATE SYSTEM
    results_F$sims.list$sigma <- results_F$sims.list$sigma * raster::res(habitat$habitat.r)[1]
    results_F$sims.list$tau <- results_F$sims.list$tau * raster::res(habitat$habitat.r)[1]
    
    
    
    ## ------   2.2. MALES -----
    
    ##-- Compile MCMC bites
    nimOutput_M <- collectMCMCbites( path = file.path(working.dir, "nimbleOutFiles/male"),
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
    
    ##-- RESCALE sigma AND tau TO THE ORIGINAL COORDINATE SYSTEM
    results_M$sims.list$sigma <- results_M$sims.list$sigma * raster::res(habitat$habitat.r)[1]
    results_M$sims.list$tau <- results_M$sims.list$tau * raster::res(habitat$habitat.r)[1]
    
    
    
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
    resultsSXYZ_MF$sims.list$sigma <- cbind(results_M$sims.list$sigma[1:minIter],
                                            results_F$sims.list$sigma[1:minIter])
    dimnames(resultsSXYZ_MF$sims.list$sigma)[[2]] <- c("M","F")
    
    ##-- sex
    resultsSXYZ_MF$sims.list$sex <- rep(c("M","F"),
                                        c(dim(resultsSXYZ_M$sims.list$sxy)[2],
                                          dim(resultsSXYZ_F$sims.list$sxy)[2]))
    
    ##-- SAVE AND LOAD DATA
    save( results_F, results_M, resultsSXYZ_MF,
          file = file.path( working.dir, "data", paste0("MCMC_bear_", DATE, ".RData")))
  }

  ##-- Number of activity center posterior samples
  n.mcmc <- dim(resultsSXYZ_MF$sims.list$z)[1]
  
  
  
  ## ------ 3. EXTRACT DENSITY -----
  
  message("## Processing density outputs...")
  
  ##-- Select niter iterations randomly
  if(n.mcmc >= niter){
    iter <- round(seq(1, n.mcmc, length.out = niter))
  } else {
    message(paste0( "The number of MCMC samples available (", n.mcmc,
                    ") is less than niter = ", niter,
                    ".\nusing niter = ", n.mcmc, " instead."))
    iter <- 1:n.mcmc
  }

  ##-- Remove buffer from the habitat
  ## [PD] maybe remove?
  # habitat.rWthBuffer <- habitat$habitat.rWthBuffer
  # habitat.rWthBuffer[habitat.rWthBuffer[] %in% 0] <- NA
  # searchedPolygon <- raster::rasterToPolygons( habitat.rWthBuffer,
  #                                              dissolve = T,
  #                                              function(x) x == 1)
  
  ##-- Habitat raster with extent used in the model
  habitatPolygon5km <- raster::crop( extraction.raster$Habitat,
                                     habitat$habitat.r)
  
  ##-- Create 5km raster of carnivore regions for extraction
  rrRegions <- extraction.raster$Regions
  rrRegions <- raster::mask(rrRegions, habitat$habitat.poly)
  rrRegions <- raster::crop(rrRegions, habitat$habitat.r)
  
  ##-- Create 5km raster of counties for extraction
  rrCounties <- extraction.raster$Counties
  rrCounties <- raster::mask(rrCounties, habitat$habitat.poly)
  rrCounties <- raster::crop(rrCounties, habitat$habitat.r)
  
  ##-- Calculate density only if necessary
  ##-- Check that a file with that name does not already exist to avoid overwriting
  densTest <- TRUE
  if(!overwrite){
    fileName <- paste0("Density_bear_", DATE, ".RData")
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
    densityInputRegions <- suppressWarnings(getDensityInput(
      regions = rrRegions,
      habitat = habitatPolygon5km,
      s = resultsSXYZ_MF$sims.list$sxy[iter, , , ],
      plot.check = F))
    
    inputRaster <- densityInputRegions$regions.r
    
    ##-- Subset to regions of interest
    regions.names <- c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6","Region 7","Region 8")
    regionID <- densityInputRegions$regions.rgmx
    row.names(regionID) <- row.names(densityInputRegions$regions.rgmx)
    regionID <- as.matrix(regionID[row.names(regionID) %in% regions.names, ])
    
    
    ##-- Get the objects to run the density function
    densityInputCounties <- suppressWarnings(getDensityInput(
      regions = rrCounties,
      habitat = habitatPolygon5km,
      s = resultsSXYZ_MF$sims.list$sxy[iter, , , ],
      plot.check = F))
    
    ##-- Subset to Counties of interest
    county.names <- COUNTIES$county[COUNTIES$country == "NOR"]
    countyID <- densityInputCounties$regions.rgmx
    row.names(countyID) <- row.names(densityInputCounties$regions.rgmx)
    countyID <- as.matrix(countyID[row.names(countyID) %in% county.names, ])
    
    
    
    ## ------   2. AC-BASED DENSITY (5km) ------
    
    ## ------     2.1. MALE & FEMALES -----
    
    ACdensity <- list()
    for(t in 1:n.years){
      ACdensity[[t]] <- GetDensity(
        sx = densityInputRegions$sx[ , ,t],
        sy = densityInputRegions$sy[ , ,t],
        z = resultsSXYZ_MF$sims.list$z[ , ,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = rbind(regionID,countyID),
        returnPosteriorCells = F)
    }#t
    names(ACdensity) <- years
    
    
    
    ## ------     2.2. MALE -----
    
    IDMales <- which(resultsSXYZ_MF$sims.list$sex == "M")
    
    ACdensityM <- list()
    for(t in 1:n.years){
      ACdensityM[[t]] <- GetDensity(
        sx = densityInputRegions$sx[ ,IDMales,t],
        sy =  densityInputRegions$sy[ ,IDMales,t],
        z = resultsSXYZ_MF$sims.list$z[ ,IDMales,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = rbind(regionID,countyID),
        returnPosteriorCells = F)
    }#t
    names(ACdensityM) <- years
    
    
    
    ## ------     2.3. FEMALE -----
    
    IDFemales <- which(resultsSXYZ_MF$sims.list$sex == "F")
    
    ACdensityF <- list()
    for(t in 1:n.years){
      ACdensityF[[t]] <- GetDensity(
        sx = densityInputRegions$sx[ ,IDFemales,t],
        sy = densityInputRegions$sy[ ,IDFemales,t],
        z = resultsSXYZ_MF$sims.list$z[ ,IDFemales,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = rbind(regionID,countyID),
        returnPosteriorCells = F)
    }
    names(ACdensityF) <- years
    
    
    
    ## ------   3. UD-BASED DENSITY (5km) ------
    
    ##-- Combine male and female sigma
    sigma <- do.call(cbind,lapply(resultsSXYZ_MF$sims.list$sex,
                                  function(x){
                                    if(x == "M"){
                                      resultsSXYZ_MF$sims.list$sigma[iter,"M"]
                                    } else {
                                      resultsSXYZ_MF$sims.list$sigma[iter,"F"]
                                    }
                                  }))
    
    ##-- Rescale sigma to the raster resolution
    sigma <- sigma/raster::res(rrRegions)[1]
    
    
    
    ## ------     3.1. MALE -----
    
    IDMales <- which(resultsSXYZ_MF$sims.list$sex=="M")
    
    UDdensityM <- list()
    for(t in 1:n.years){
      UDdensityM[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[ ,IDMales,t],
        sy = densityInputRegions$sy[ ,IDMales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDMales,t],
        sigma = sigma[ ,IDMales],
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = F)
    }#t
    names(UDdensityM) <- years
    
    
    
    ## ------     3.2. FEMALE -----
    
    IDFemales <- which(resultsSXYZ_MF$sims.list$sex=="F")
    
    UDdensityF <- list()
    for(t in 1:n.years){
      UDdensityF[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[ ,IDFemales,t],
        sy = densityInputRegions$sy[ ,IDFemales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDFemales,t],
        sigma = sigma[ ,IDFemales],
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = F)
    }#t
    names(UDdensityF) <- years
    
    
    
    ## ------     3.3. MALE & FEMALES -----
    
    UDdensity <- list()
    for(t in 1:n.years){
      UDdensity[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[ , ,t],
        sy = densityInputRegions$sy[ , ,t],
        z = resultsSXYZ_MF$sims.list$z[iter, ,t],
        sigma = sigma,
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = F)
    }#t
    names(UDdensity) <- years
    
    
    
    ## ------   4. SAVE DENSITY OBJECTS ------
    
    save( inputRaster,
          ACdensity,
          ACdensityF,
          ACdensityM,
          UDdensity,
          UDdensityF,
          UDdensityM,
          file = file.path( working.dir, "data",
                            paste0("Density_bear_",DATE,".RData")))
  }
  
  
  
  ## ------ 4. FIGURES -----
  
  ##-- Plot parameters
  diffSex <- 0.2
  colSex <- c("firebrick2", "deepskyblue2", "black")
  colCause  <- adjustcolor( c("#E69F00","#009E73"), 0.5)
  
  
  
  ## ------   4.1. DENSITY MAPS -----
  
  message("## Plotting population density maps...") 
  
  ##-- Create 5km raster for plotting
  rrNorway <- extraction.raster[["Countries"]]
  rrNorway[!rrNorway[] %in% 2] <- NA
  rrNorway[rrNorway[] %in% 2] <- 1
  rrNorway <- raster::crop(rrNorway, habitat$habitat.r)
  rrCombined <- rrRegions + rrNorway
  
  ##-- AC-density maps
  plotDensityMaps(
    input = inputRaster,
    estimates = ACdensity,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES[1, ],
    type = c("time.series", "last.year"),
    path = working.dir,
    name = "AC_Density")
  
  ##-- UD-density maps
  plotDensityMaps( 
    input = inputRaster,
    estimates = UDdensity,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES[1, ],
    type = c("time.series", "last.year","summary","summary_NOR"),
    species = "bear",
    labels = list("nor" = ACdensity[[n.years]]$summary["Total",c("95%CILow","95%CIHigh")]),
    x.labels = 0.3,
    y.labels = 0.8,
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
  
  graphics::par(mar = c(5,5,1,1))
  
  n.detected <- apply(rbind(nimDataM$y.alive[ ,1, ],nimDataF$y.alive[ ,1, ]) > 0, 2, sum)
  
  ymax <- 10*(trunc(max(c(unlist(lapply(ACdensity, function(x)max(colSums(x$PosteriorAllRegions)))), n.detected))/10)+1)
  
  plot(-1000,
       xlim = c(0.5, n.years + 0.5),
       ylim = c(0,ymax),
       xlab = "", ylab = paste("Number of bears"),
       xaxt = "n", axes = F, cex.lab = 1.6)
  graphics::axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
  graphics::axis(2, at = seq(0,ymax,20), labels = seq(0,ymax,20), cex.axis = 1.6)
  graphics::abline(v = (0:n.years)+0.5, lty = 2)
  graphics::abline(h = seq(0,ymax, by = 10), lty = 2, col = "gray90")
  
  for(t in 1:n.years){
    ##-- FEMALES
    plotQuantiles(x = colSums(ACdensityF[[t]]$PosteriorAllRegions),
                  at = t - diffSex,
                  width = 0.18,
                  col = colSex[1])
    
    ##-- MALES
    plotQuantiles(x = colSums(ACdensityM[[t]]$PosteriorAllRegions),
                  at = t + diffSex,
                  width = 0.18,
                  col = colSex[2])
    ##-- TOTAL
    plotQuantiles(x = colSums(ACdensity[[t]]$PosteriorAllRegions),
                  at = t,
                  width = 0.4,
                  col = colSex[3])
    
    ##-- ADD NUMBER OF INDIVIDUALS DETECTED
    xx <- c(t-0.25,t+0.25,t+0.25,t-0.25)
    yy <- c(n.detected[t]-1,n.detected[t]-1,n.detected[t]+1,n.detected[t]+1)
    polygon(xx, yy, border = NA, col = "goldenrod1")
  }#t
  box()
  
  ##-- legend
  par(xpd = TRUE)
  xx <- c(0.6*n.years,0.72*n.years,0.82*n.years,0.92*n.years)#c(6.1,7.6,8.8,10) 
  yy <- c(5,5,5,5)
  labs <- c("Females", "Males", "Total", "Detected")
  polygon(x = c(0.58*n.years,1.05*n.years,1.05*n.years,0.58*n.years),
          y = c(-2.5,-2.5,12.5,12.5),
          col = adjustcolor("white", alpha.f = 0.9),
          border = "gray60")
  
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col =  adjustcolor(colSex,0.3))
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col =  adjustcolor(colSex,0.7))
  
  text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.2, pos = 4)
  
  polygon(x = c(xx[4]-0.3,xx[4]+0.1,xx[4]+0.1,xx[4]-0.3),
          y = c(yy[4]-1,yy[4]-1,yy[4]+1,yy[4]+1),
          col = "goldenrod1",
          border = F)
  
  dev.off()

  ##-- Remove unnecessary objects from memory
  rm(list = c("n.detected"))
  gc(verbose = FALSE)
  
  
  
  ## ------   4.3. ABUNDANCE per REGION ------

  regPlotNames <- c("Region 8","Region 7","Region 6","Region 5","Region 3")

  grDevices::png(filename = file.path(working.dir,"figures/Abundance_Regions.png"),
                 width = 8, height = 15, units = "in", pointsize = 12,
                 res = 300, bg = NA)

  nf <- layout(rbind(c(1,2),
                     c(3,4),
                     c(5,6),
                     c(7,8),
                     c(9,10)),
               widths = c(1,0.5),
               heights = 1)

  for(cc in regPlotNames){

    ymax <- max(c(10*(trunc(max(ACdensity[[n.years]]$PosteriorRegions[cc, ])/10)+1),
                60))

    par(mar = c(2,5,2,0), tck = 0, mgp = c(2,0.2,0))
    plot(-1000,
         xlim = c(0.5, n.years + 0.5), ylim = c(0,ymax),
         xlab = "", ylab = paste("Number of bears"),
         xaxt = "n", axes = F, cex.lab = 1.6)
    axis(1, at = c(1:n.years), labels = years, cex.axis = 1.15)
    axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10),
         cex.axis = 1.2, las = 1)
    abline(v = (0:n.years)+0.5, lty = 2)
    abline(h = seq(0,ymax, by = 10), lty = 2, col = "gray60")
    mtext(text = cc, side = 3, font = 2, line = 0.5)

    for(t in 1:n.years){
      ##-- TOTAL
      plotQuantiles(x = ACdensity[[t]]$PosteriorRegions[cc, ],
                    at = t,
                    width = 0.4,
                    col = colSex[3])

      ##-- FEMALES
      plotQuantiles(x = ACdensityF[[t]]$PosteriorRegions[cc, ],
                    at = t - diffSex,
                    width = 0.18,
                    col = colSex[1])

      ##-- MALES
      plotQuantiles(x = ACdensityM[[t]]$PosteriorRegions[cc, ],
                    at = t + diffSex,
                    width = 0.18,
                    col = colSex[2])
    }#t

    ##-- legend
    if(cc == "Region 3"){
      par(xpd = TRUE)
      xx <- c(1.1,3.2,5)
      yy <- 0.9*ymax
      labs <- c("Females", "Males", "Total")
      polygon(x = c(0.6,6.4,6.4,0.6),
              y = c(yy-5,yy-5,yy+5,yy+5),
              col = adjustcolor("white", alpha.f = 0.9),
              border = "gray60")

      points(x = xx, y = rep(yy,3),  pch = 15, cex = 3.5, col =  adjustcolor(colSex,0.3))
      points(x = xx, y = rep(yy,3),  pch = 15, cex = 1.5, col =  adjustcolor(colSex,0.7))
      text(x = xx + 0.1, y = rep(yy,3), labels = labs, cex = 1.4, pos = 4)
    }
    box()

    ##-- Map insert
    par(mar = c(5,0,4,5))
    plot(st_geometry(REGIONS), border = grey(0.5), col = grey(0.5), lwd = 0.1)
    plot(st_geometry(REGIONS[REGIONS$region == cc, ]),
         add = T, col = adjustcolor("red",0.5), border = "red")
  }#c
  dev.off()



  ##----------------------------------------------------------------------------
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
  # points(c(4,4), c(4,3), pch = 15, cex = 5.5, col =  adjustcolor(colSex,0.3))
  # points(c(4,4), c(4,3), pch = 15, cex = 3, col =  adjustcolor(colSex,0.7))
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
  # points(c(4,4,4), c(4,3,2), pch = 15, cex = 5.5, col =  adjustcolor(colSex,0.3))
  # points(c(4,4,4), c(4,3,2), pch = 15, cex = 3, col =  adjustcolor(colSex,0.7))
  # text(c(5.3,5.3,5.3), c(4,3,2),
  #      c("Females", "Males", "Total"), cex = 2, pos = 4)
  # 
  # dev.off()
  # 
  # 
  # 
  # ## ------   4.5. POPULATION FLUXES IN/OUT NORWAY -----
  # 
  # message("## Plotting population fluxes...")
  # 
  # norRaster <- habitatRasterResolution$`5km`[["Countries"]]
  # norRaster[norRaster[] != 2] <- NA
  # 
  # ##-- Calculate number of individuals
  # N_surv <- N_surv_F <- N_surv_M <- matrix(NA,n.mcmc,n.years-1)
  # N_recruit <- N_recruit_F <- N_recruit_M <- matrix(NA,n.mcmc,n.years-1)
  # N_immig <- N_immig_F <- N_immig_M <- matrix(NA,n.mcmc,n.years-1)
  # N_emig <- N_emig_F <- N_emig_M <- matrix(NA,n.mcmc,n.years-1)
  # 
  # for(t in 1:(n.years-1)){
  #   for(iter in 1:n.mcmc){
  #     isOut <- is.na(norRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+1])])
  #     wasOut <- is.na(norRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])])
  # 
  #     N_recruit_F[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & (!isOut) & isFemale)
  #     N_recruit_M[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & (!isOut) & isMale)
  #     N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]
  # 
  #     N_immig_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & wasOut & (!isOut) & isFemale)
  #     N_immig_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & wasOut & (!isOut) & isMale)
  #     N_immig[iter,t] <- N_immig_F[iter,t] + N_immig_M[iter,t]
  # 
  #     N_emig_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & isOut & isFemale)
  #     N_emig_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & isOut & isMale)
  #     N_emig[iter,t] <- N_emig_F[iter,t] + N_emig_M[iter,t]
  # 
  #     N_surv_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & (!isOut) & isFemale)
  #     N_surv_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & (!isOut) & isMale)
  #     N_surv[iter,t] <- N_surv_F[iter,t] + N_surv_M[iter,t]
  #   }#iter
  #   #print(t)
  # }#t
  # 
  # ##-- Remove unnecessary objects from memory
  # rm(list = c("isMale", "isFemale","isAlive", "isAvail"))
  # gc(verbose = FALSE)
  # 
  # 
  # 
  # ## ------     4.5.1. PLOT RECRUITS -----
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("NumRecruits_bars.pdf")),
  # #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path(working.dir, "figures/NumRecruits_bars.png"),
  #                width = 12, height = 8.5, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # ymax <- 10*(trunc(max(N_recruit)/10)+1)
  # par(mar = c(5,5,1,1))
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of bears recruited"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_recruit[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_recruit_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_recruit_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # 
  # }#t
  # box()
  # 
  # ##-- legend
  # par(xpd = TRUE)
  # xx <- c(1,2.5,4)
  # yy <- c(55,55,55)
  # labs <- c("Females", "Males", "Total")
  # polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
  #         y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
  #         col = adjustcolor("white", alpha.f = 0.9),
  #         border = "gray60")
  # 
  # points(x = xx, y = yy,  pch = 15, cex = 3.5, col =  adjustcolor(colSex,0.3))
  # points(x = xx, y = yy,  pch = 15, cex = 1.5, col =  adjustcolor(colSex,0.7))
  # text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  # 
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------     4.5.2. PLOT SURVIVORS -----
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("NumSurvival_bars.pdf")),
  # #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path(working.dir, "figures/NumSurvival_bars.png"),
  #                width = 12, height = 8.5, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # ymax <- 10*(trunc(max(N_surv)/10)+1)
  # par(mar = c(5,5,1,1))
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of surviving bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,20), labels = seq(0,ymax,20), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 10), lty = 2, col = "gray60")
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_surv[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_surv_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_surv_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # 
  # }#t
  # box()
  # 
  # ##-- legend
  # par(xpd = TRUE)
  # xx <- c(1,2.5,4)
  # yy <- c(55,55,55)
  # labs <- c("Females", "Males", "Total")
  # polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
  #         y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
  #         col = adjustcolor("white", alpha.f = 0.9),
  #         border = "gray60")
  # 
  # points(x = xx, y = yy,  pch = 15, cex = 3.5, col =  adjustcolor(colSex,0.3))
  # points(x = xx, y = yy,  pch = 15, cex = 1.5, col =  adjustcolor(colSex,0.7))
  # text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  # 
  # dev.off()
  # 
  # 
  # 
  # ## ------     4.5.3. PLOT IMMIGRANTS -----
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("NumImmigrants_bars.pdf")),
  # #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path( working.dir, "figures/NumImmigrants_bars.png"),
  #                width = 12, height = 8.5, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # par(mar = c(5,5,1,1))
  # ymax <- 10*(trunc(max(N_immig)/10)+1)
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of immigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_immig[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_immig_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_immig_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # }#t
  # box()
  # 
  # ##-- legend
  # par(xpd = TRUE)
  # xx <- c(1,2.5,4)
  # yy <- c(55,55,55)
  # labs <- c("Females", "Males", "Total")
  # polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
  #         y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
  #         col = adjustcolor("white", alpha.f = 0.9),
  #         border = "gray60")
  # 
  # points(x = xx, y = yy,  pch = 15, cex = 3.5, col = adjustcolor(colSex,0.3))
  # points(x = xx, y = yy,  pch = 15, cex = 1.5, col = adjustcolor(colSex,0.7))
  # text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  # 
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------     4.5.4. PLOT EMIGRANTS -----
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("NumEmigrants_bars.pdf")),
  # #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path(working.dir, "figures/NumEmigrants_bars.png"),
  #                width = 12, height = 8.5, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # ymax <- 10*(trunc(max(N_emig)/10)+1)
  # par(mar = c(5,5,1,1))
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of emigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # 
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_emig[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_emig_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_emig_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # 
  # }#t
  # box()
  # 
  # ##-- legend
  # par(xpd = TRUE)
  # xx <- c(1,2.5,4)
  # yy <- c(55,55,55)
  # labs <- c("Females", "Males", "Total")
  # polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
  #         y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
  #         col = adjustcolor("white", alpha.f = 0.9),
  #         border = "gray60")
  # 
  # points(x = xx, y = yy,  pch = 15, cex = 3.5, col = adjustcolor(colSex,0.3))
  # points(x = xx, y = yy,  pch = 15, cex = 1.5, col = adjustcolor(colSex,0.7))
  # text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  # 
  # dev.off()
  # 
  # 
  # 
  # ## ------     4.5.5. PLOT ALL -----
  # 
  # # pdf(file = file.path(working.dir, "figures", paste0("NumFluxes_bars.pdf")),
  # #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path(working.dir, "figures/NumFluxes_bars.png"),
  #                width = 12, height = 8.5, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
  # 
  # par(mfrow = c(2,2), mar = c(5,5,1,1))
  # 
  # ymax <- 10*(trunc(max(N_recruit)/10)+1)
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of bears recruited"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_recruit[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_recruit_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_recruit_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # }#t
  # box()
  # 
  # 
  # ymax <- 10*(trunc(max(N_surv)/10)+1)
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of surviving bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,20), labels = seq(0,ymax,20), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 10), lty = 2, col = "gray60")
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_surv[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_surv_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_surv_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # 
  # }#t
  # box()
  # 
  # 
  # ymax <- 10*(trunc(max(N_immig)/10)+1)
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of immigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_immig[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_immig_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_immig_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # }#t
  # box()
  # 
  # ymax <- 10*(trunc(max(N_emig)/10)+1)
  # plot(-1000,
  #      xlim = c(0.5, (n.years-1) + 0.5), ylim = c(0,ymax),
  #      xlab = "", ylab = paste("Number of emigrant bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:(n.years-1)), labels = years[2:n.years], cex.axis = 1.6)
  # axis(2, at = seq(0,ymax,10), labels = seq(0,ymax,10), cex.axis = 1.6)
  # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = seq(0,ymax, by = 5), lty = 2, col = "gray60")
  # for(t in 1:(n.years-1)){
  #   ##-- TOTAL
  #   plotQuantiles(x = N_emig[ ,t],
  #                 at = t,
  #                 width = 0.4,
  #                 col = colSex[3])
  # 
  #   ##-- FEMALES
  #   plotQuantiles(x = N_emig_F[ ,t],
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[1])
  # 
  #   ##-- MALES
  #   plotQuantiles(x = N_emig_M[ ,t],
  #                 at = t + diffSex,
  #                 width = 0.18,
  #                 col = colSex[2])
  # 
  # }#t
  # box()
  # 
  # ##-- legend
  # par(xpd = TRUE)
  # xx <- c(1,2.5,4)
  # yy <- c(0.8*ymax,0.8*ymax,0.8*ymax)
  # labs <- c("Females", "Males", "Total")
  # polygon(x = c(min(xx)-0.5,max(xx)+1.5,max(xx)+1.5,min(xx)-0.5),
  #         y = c(min(yy)-5,min(yy)-5,min(yy)+5,min(yy)+5),
  #         col = adjustcolor("white", alpha.f = 0.9),
  #         border = "gray60")
  # 
  # points(x = xx, y = yy,  pch = 15, cex = 3.5, col = adjustcolor(colSex,0.3))
  # points(x = xx, y = yy,  pch = 15, cex = 1.5, col = adjustcolor(colSex,0.7))
  # text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.4, pos = 4)
  # 
  # dev.off()
  # 
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
  # for(c in 1:nrow(COUNTIES_s)){
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
  #     plotQuantiles( x = results_F$sims.list$p0[ ,COUNTIES_s$index == c,t],
  #                    at = t - diffSex,
  #                    col = colSex[1])
  # 
  #     plotQuantiles( x = results_M$sims.list$p0[ ,COUNTIES_s$index == c,t],
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
  #     points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 5.5, col = adjustcolor(colSex,0.3))
  #     points(c(0.8,0.8), c(0.0077,0.0092), pch = 15, cex = 3, col = adjustcolor(colSex,0.7))
  #     text(c(1.2,1.2),c(0.0077,0.0092),  c("Females", "Males"), cex = 2, pos = 4)
  #   }
  # 
  # 
  #   par(mar = c(0,0,0,0))
  #   plot(st_geometry(COUNTIES_s), border = grey(0.5), col = grey(0.5), lwd = 0.1)
  #   plot(st_geometry(COUNTIES_s[COUNTIES_s$index == c, ]),
  #        add = T, col = adjustcolor("red",0.5), border = "red")
  #   text(COUNTIES_s[COUNTIES_s$index == c, ],
  #        labels = COUNTIES_s$Name[COUNTIES_s$index == c],
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
  # points(c(4,4), c(4,3), pch = 15, cex = 5.5, col = adjustcolor(colSex,0.3))
  # points(c(4,4), c(4,3), pch = 15, cex = 3, col = adjustcolor(colSex,0.7))
  # text(c(5.3,5.3), c(4,3),  c("Females", "Males"), cex = 1.5, pos = 4)
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------   4.7. NGS, Dead recoveries & Carnivore obs ------
  # 
  # ##-- Plot NGS & Dead recovery maps
  # # pdf(file = file.path(working.dir, "figures", "NGS_DR_maps.pdf"),
  # #     width = 18, height = 12)
  # grDevices::png(filename = file.path(working.dir, "figures/NGS_DR_maps.png"),
  #                width = 18, height = 12, units = "in", pointsize = 12,
  #                res = 300, bg = NA)
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
  #   plot(sf::st_geometry(COUNTIES_s), border = NA, col = "gray80")
  #   points(data.alive$data.sp[data.alive$data.sp$Year == years[t], ],
  #          pch = 3, col = "orange", lwd = 0.7)
  #   points(data.dead[data.dead$Year == years[t], ],
  #          pch = 3, col = "slateblue", lwd = 0.7)
  #   mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
  #   
  #   if(t == n.years){
  #     ##-- LEGEND
  #     xLeg <- 830000
  #     yLeg <- 6730000
  #     segments(x0 = xLeg, x1 = xLeg,
  #              y0 = yLeg, y1 = yLeg + 500000,
  #              col = grey(0.3), lwd = 4, lend = 2)
  #     text(xLeg-80000, yLeg+500000/2, labels = "500 km", srt = 90, cex = 2)
  #     
  #     points(x = c(xLeg-200000,xLeg-200000),
  #            y = c(yLeg-100000,yLeg-180000),
  #            pch = 3, lwd = 1.5, cex = 3,
  #            col = c("orange","slateblue"))
  #     text(x = c(xLeg-150000,xLeg-150000),
  #          y = c(yLeg-100000,yLeg-180000),
  #          c("NGS samples", "Dead recoveries"), cex = 2, pos = 4)
  #   }#if
  # }#t
  # dev.off()
  # 
  # 
  # 
  # # ##-- Plot Carnivore observations maps
  # # pdf(file = file.path(working.dir, "figures", paste0("CarnivoreObs_maps_classic.pdf")),
  # #     width = 18, height = 12)
  # # grDevices::png(filename = file.path(working.dir, "figures/CarnivoreObs_maps_classic.png"),
  # #     width = 18, height = 12, units = "in", pointsize = 12,
  # #     res = 300, bg = NA)
  # #
  # # ##-- layout
  # # mx <- rbind(c(1,rep(1:5, each = 2)),
  # #             c(rep(1:5, each = 2), 5))
  # # mx <- rbind(mx, mx + 5)
  # # nf <- layout(mx,
  # #              widths = c(rep(1,ncol(mx))),
  # #              heights = rep(1,2))
  # # par(mar = c(0,0,0,0))
  # # for(t in 1:length(years)){
  # #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
  # #   image(mask(ds.brickCont[[t]],COUNTRIESsimpFig[1,]), add = TRUE, col = c("white","forestgreen"), legend = FALSE)
  # #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = "gray80", col = NA, add = TRUE)
  # #   
  # #   mtext(text = years[t], side = 1, -25, adj=0.2, cex=1.8, font = 2)
  # #   
  # #   if(t == n.years){
  # #     segments(x0 = 830000, x1 = 830000,
  # #              y0 = 6730000, y1 = 6730000 + 500000,
  # #              col = grey(0.3), lwd = 4, lend = 2)
  # #     text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 2)
  # #     
  # #     ##-- LEGEND
  # #     par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  # #     plot(1, ylim = c(-1,7), xlim = c(0,15), type = "n", axes = FALSE)
  # #   }#if
  # # }#t
  # # dev.off()
  # 
  # 
  # 
  # ## ------   4.7. DETECTABILITY ------
  # 
  # ## ------     4.7.1. SET-UP ------
  # 
  # ##-- Create 5km raster for extraction
  # regions.r <- habitatRasterResolution$`5km`[["Regions"]]
  # regions.r <- raster::crop(regions.r, habitat$habitat.r)
  # 
  # ##-- Set all habitat cells outside Norway to NA
  # regions.r[regions.r[] < 23] <- NA
  # 
  # ##-- Identify habitat cells
  # isHabitat <- which(!is.na(regions.r[]))
  # 
  # ##-- Get habitat cell coordinates
  # regions.xy <- raster::coordinates(regions.r)[isHabitat, ]
  # 
  # ##-- Get region name for each habitat cell (format to run the C++ function)
  # regionsNames <- sort(unique(na.omit(regions.r[])))
  # regions.rgmx <- do.call(rbind, lapply(regionsNames, function(x)regions.r[] == x))
  # regions.rgmx[is.na(regions.rgmx)] <- 0
  # row.names(regions.rgmx) <- raster::factorValues(regions.r, regionsNames)[,1]
  # regions.rgmx <- regions.rgmx[ ,isHabitat]
  # 
  # ##-- Get detectors coordinates (original and scaled)
  # detectors.xy <- as.matrix(st_coordinates(detectors$main.detector.sp))
  # 
  # ##-- Load detectability calculation
  # load(file.path(working.dir, "data/Detectability5km.RData"))
  # 
  # 
  # 
  # ## ------     4.7.2. PLOT DETECTABILITY MAPS -----
  # 
  # pdf(file = file.path(working.dir, "figures/Detectability_maps.pdf"),
  #     width = 10, height = 6)
  # grDevices::png(filename = file.path(working.dir, "figures/Detectability_maps.png"),
  #     width = 10, height = 6, units = "in", pointsize = 12,
  #     res = 300, bg = NA)
  #
  # ##-- Set color scale
  # max <- max(c(unlist(lapply( DetectabilityRegionsM[[t]]$MeanCell,
  #                             function(x) max(x[], na.rm = T))),
  #              unlist(lapply( DetectabilityRegionsF[[t]]$MeanCell,
  #                             function(x) max(x[], na.rm = T)))))
  # cuts <- seq(0, max, length.out = 100) ##-- set breaks
  # col <- rev(terrain.colors(100))
  # 
  # ##-- layout
  # mx <- rbind(c(1,rep(1:5, each = 2)),
  #             c(rep(1:5, each = 2), 5))
  # mx <- rbind(mx, mx + 5)
  # nf <- layout(mx,
  #              widths = c(rep(1,ncol(mx))),
  #              heights = rep(1,2))
  # par(mar = c(0,0,0,0))
  # 
  # ##-- FEMALES
  # for(t in 1:length(years)){
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  #   detectab.r <- regions.r
  #   detectab.r[isHabitat] <- DetectabilityRegionsF[[t]]$MeanCell
  #   image(detectab.r, add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])),
  #        border = grey(0.4), col = NA, add = TRUE)
  #   mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)
  #   
  #   if(t == n.years){
  #     segments(x0 = 830000, x1 = 830000,
  #              y0 = 6730000, y1 = 6730000 + 500000,
  #              col = grey(0.3), lwd = 4, lend = 2)
  #     text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  #     plot( detectab.r,
  #           legend.only = T,
  #           breaks = cuts,
  #           col = col,
  #           legend.width = 2,
  #           axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            cex.axis = 1.2),
  #           smallplot = c(0.95, 1.00, 0.2, 0.6),
  #           legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
  #                              side = 2, font = 1, line = 1, cex = 1))
  #   }#if
  # }#t
  # 
  # ##-- MALES
  # for(t in 1:length(years)){
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  #   detectab.r <- regions.r
  #   detectab.r[isHabitat] <- DetectabilityRegionsM[[t]]$MeanCell
  #   image(detectab.r, add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])),
  #        border = grey(0.4), col = NA, add = TRUE)
  #   mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)
  #   
  #   if(t == n.years){
  #     segments(x0 = 830000, x1 = 830000,
  #              y0 = 6730000, y1 = 6730000 + 500000,
  #              col = grey(0.3), lwd = 4, lend = 2)
  #     text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  #     plot( detectab.r,
  #           legend.only = T,
  #           breaks = cuts,
  #           col = col,
  #           legend.width = 2,
  #           axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            cex.axis = 1.2),
  #           smallplot = c(0.95, 1.00, 0.2, 0.6),
  #           legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
  #                              side = 2, font = 1, line = 1, cex = 1))
  #   }#if
  # }#t
  # dev.off()
  # 
  # 
  # 
  # ## ------   4.8. SEX-RATIO ------
  # 
  # ## ------     4.8.1. SEX-RATIO BARS ------
  # 
  # pdf(file = file.path(working.dir, "figures/SexRatio_bars.pdf"),
  #     width = 12, height = 8.5)
  # grDevices::png(filename = file.path(working.dir, "figures/SexRatio_bars.png"),
  #     width = 12, height = 8.5, units = "in", pointsize = 12,
  #     res = 300, bg = NA)
  #
  # par(mar = c(5,5,1,1))
  # plot(-1000,
  #      xlim = c(0.5, n.years + 0.5), ylim = c(0.2,0.8),
  #      xlab = "", ylab = paste("% Female bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  # axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
  # axis(2, at = seq(0.2,0.8,0.1), labels = seq(0.2,0.8,0.1), cex.axis = 1.6)
  # # abline(v = (0:n.years)+0.5, lty = 2)
  # abline(h = 0.5, lty = 2, col = "gray60")
  # for(t in 1:n.years){
  #   plotQuantiles(x = colSums(ACdensityF[[t]]$PosteriorAllRegions)/
  #                   colSums(ACdensity[[t]]$PosteriorAllRegions),
  #                 at = t - diffSex,
  #                 width = 0.18,
  #                 col = colSex[3])
  # }#t
  # #box()
  # dev.off()
  # 
  # 
  # 
  # ## ------     4.8.2. SEX-RATIO DISTANCE ------
  # ###--- TBD
  # 
  # 
  # 
  # ## ------     4.8.3. SEX-RATIO MAPS ------
  # 
  # ##-- Convert densities from 25km2 (5*5 raster) to 100km2
  # SexRatio <- list()
  # for(t in 1:length(years)){
  #   SexRatio[[t]] <- ACdensityF[[t]]$MeanCell/ACdensity[[t]]$MeanCell
  # }
  # 
  # ##-- Crop density maps to Norway
  # rrCombined <- rrRegions + rrNorway
  # SexRatioMap <- list()
  # for(t in 1:length(years)){
  #   SexRatioMap[[t]] <- densityInputRegions$regions.r
  #   SexRatioMap[[t]][] <- NA
  #   SexRatioMap[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
  #   SexRatioMap[[t]][is.na(rrCombined[])] <- NA
  #   
  #   crs(SexRatioMap[[t]]) <- st_crs(habitat$habitat.poly)
  # }#t
  # 
  # 
  # # pdf(file = file.path(working.dir, "figures/SexRatio_Maps.pdf"),
  # #    width = 12, height = 8)
  # grDevices::png(filename = file.path(working.dir, "figures/SexRatio_Maps.png"), 
  #     width = 12, height = 8, units = "in", pointsize = 12,
  #     res = 300, bg = NA)
  # 
  # ##-- Set color scale
  # max <- 1
  # cuts <- seq(0, max, length.out = 100) ##-- set breaks
  # colfunc <- colorRampPalette(c(colSex[2],"white", colSex[1]))
  # col <- colfunc(100)
  # 
  # ##-- layout
  # mx <- rbind(c(1,rep(1:5, each = 2)),
  #             c(rep(1:5, each = 2), 5))
  # mx <- rbind(mx, mx + 5)
  # nf <- layout(mx,
  #              widths = c(rep(1,ncol(mx))),
  #              heights = rep(1,2))
  # #layout.show(nf)
  # par(mar = c(0,0,0,0))
  # 
  # ##-- Plot Sex-ratio maps
  # for(t in 1:length(years)){
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
  #   image(SexRatioMap[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  #   plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
  #   mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)
  #   
  #   if(t == n.years){
  #     segments(x0 = 830000, x1 = 830000,
  #              y0 = 6730000, y1 = 6730000 + 500000,
  #              col = grey(0.3), lwd = 4, lend = 2)
  #     text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  #     plot( SexRatioMap[[t]],
  #           legend.only = T,
  #           breaks = cuts,
  #           col = col,
  #           legend.width = 2,
  #           axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                            cex.axis = 1.2),
  #           smallplot = c(0.95, 1.00, 0.2, 0.6),
  #           legend.args = list(text = "% Female bears",
  #                              side = 2, font = 1, line = 1, cex = 1))
  #   }#if
  # }#t
  # dev.off()
  # 
  # 
  # 
  # 
  # ## ------     4.8.4. SEX-RATIO MAPS (SMOOTHED) ------
  # 
  # ##-- Set smoothing factor
  # smoothingFactor <- c(3,7,11,21)
  # for(SF in smoothingFactor){
  #   
  #   SexRatioMap.smooth <- list()
  #   for(t in 1:length(years)){
  #     SexRatioMap.smooth[[t]] <- densityInputRegions$regions.r
  #     SexRatioMap.smooth[[t]][] <- NA
  #     SexRatioMap.smooth[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
  #     SexRatioMap.smooth[[t]] <- terra::focal(terra::rast(SexRatioMap[[t]]), SF, "mean", na.rm=TRUE)
  #     SexRatioMap.smooth[[t]][is.na(rrCombined[])] <- NA
  #   }#t
  #   
  #   
  # # pdf(file = file.path(working.dir, "figures", paste0("SexRatio_Maps_smooth_",SF,".pdf")), 
  # #      width = 12, height = 8)
  #   grDevices::png(filename = file.path(working.dir, "figures", paste0("SexRatio_Maps_smooth_",SF,".png")), 
  #       width = 12, height = 8, units = "in", pointsize = 12,
  #       res = 300, bg = NA)
  #   
  #   ##-- Set color scale
  #   max <- 1
  #   cuts <- seq(0, max, length.out = 100) ##-- set breaks
  #   colfunc <- colorRampPalette(c(inferno(100)[1], inferno(100)[25], inferno(100)[50], inferno(100)[75], inferno(100)[100]))
  #   col <- colfunc(100)
  #   
  #   
  #   ##-- layout
  #   mx <- rbind(c(1,rep(1:5, each = 2)),
  #               c(rep(1:5, each = 2), 5))
  #   mx <- rbind(mx, mx + 5)
  #   nf <- layout(mx,
  #                widths = c(rep(1,ncol(mx))),
  #                heights = rep(1,2))
  #   #layout.show(nf)
  #   par(mar = c(0,0,0,0))
  #   
  #   ##-- Plot AC maps
  #   for(t in 1:length(years)){
  #     plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  #     image( SexRatioMap.smooth[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
  #     plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
  #     mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)
  #     
  #     if(t == n.years){
  #       segments(x0 = 830000, x1 = 830000,
  #                y0 = 6730000, y1 = 6730000 + 500000,
  #                col = grey(0.3), lwd = 4, lend = 2)
  #       text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  #       plot( raster(SexRatioMap.smooth[[t]]),
  #             legend.only = T,
  #             breaks = cuts,
  #             col = col,
  #             legend.width = 2,
  #             axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                              labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
  #                              cex.axis = 1.2),
  #             smallplot = c(0.82, 0.84, 0.2, 0.6),
  #             legend.args = list(text = "% Female bears",
  #                                side = 2, font = 1, line = 0.2, cex = 1))
  #     }#if
  #   }#t
  #   dev.off()
  # }#SF
  # 
  # 
  # 
  # 
  ##----------------------------------------------------------------------------
  ## ------ 5. TABLES -----
  
  gc()
  
  ## ------   5.1. ABUNDANCE -----
  
  ## ------     5.1.1. ALL YEARS, BOTH SEX, PER REGION -----
  
  idregion <- row.names(ACdensity[[1]]$summary)
  
  ##-- Remove Finland, Norway, Russia, Sweden
  idregion <- idregion[-which(idregion %in% c("Finland","Norway","Russia","Sweden","Total"))]
  idregion <- sort(unique(idregion))
  
  ##-- Get names of Norwegian carnivore regions
  idregionNOR <- idregion[grep("Region",idregion)]
  idregionTable <- c(idregionNOR, "Total")
  
  ##-- Create table to store N estimates (and CI)
  NCarRegionEstimates <- NCountyMean <- matrix("", ncol = n.years, nrow = length(idregionTable))
  row.names(NCarRegionEstimates) <- row.names(NCountyMean) <- c(idregionTable)
  colnames(NCarRegionEstimates) <- colnames(NCountyMean) <- years
  
  ##-- Fill-in the table
  for(t in 1:n.years){
    for(i in 1:length(idregionTable)){
      NCarRegionEstimates[idregionTable[i],t] <-
        paste0(round(ACdensity[[t]]$summary[idregionTable[i],"mean"],digits = 1), " (",
               round(ACdensity[[t]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
               round(ACdensity[[t]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
      
      NCountyMean[idregionTable[i],t] <- round(ACdensity[[t]]$summary[idregionTable[i],"mean"])
    }#i
  }#t
  
  ##-- Adjust names
  idregion1 <- idregionTable
  idregion1[which(idregion1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idregion1
  
  ##-- print .csv
  write.csv(NCarRegionEstimates,
            file = file.path(working.dir, "tables/N_AllYears_region.csv"))
  
  ##-- print .tex
  row.names(NCarRegionEstimates) <- c(paste0("\\hspace{0.1cm} ", idregionNOR), "TOTAL")
  print(xtable::xtable( NCarRegionEstimates,
                        type = "latex",
                        align = paste(c("l",rep("c",ncol(NCarRegionEstimates))), collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list(list(seq(1, nrow(NCarRegionEstimates), by = 2)), "\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/N_AllYears_region.tex"))
  
  
  
  
  ## ------     5.1.2. LAST YEAR, PER SEX, PER REGION -----
  
  NCarRegionEstimatesLast <- matrix("", ncol = 3, nrow = length(idregionTable))
  row.names(NCarRegionEstimatesLast) <- c(idregionTable)
  colnames(NCarRegionEstimatesLast) <- c("Females","Males","Total")
  
  for(i in 1:length(idregionTable)){
    ##-- FEMALES
    NCarRegionEstimatesLast[idregionTable[i],"Females"] <-
      paste0(round(ACdensityF[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(ACdensityF[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensityF[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- MALES
    NCarRegionEstimatesLast[idregionTable[i],"Males"] <-
      paste0(round(ACdensityM[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(ACdensityM[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensityM[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- BOTH SEXES
    NCarRegionEstimatesLast[idregionTable[i],"Total"] <-
      paste0(round(ACdensity[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(ACdensity[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensity[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
  }#i
  
  ##-- ADJUST NAMES
  idregion1 <- idregionTable
  idregion1[which(idregion1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idregion1
  
  ##-- print .csv
  write.csv( NCarRegionEstimatesLast,
             file = file.path(working.dir, "tables/N_LastYearPerSex_region.csv"))
  
  ##-- print .tex
  row.names(NCarRegionEstimatesLast) <- c(paste0("\\hspace{0.1cm} ",idregionNOR),"TOTAL")
  print(xtable::xtable(NCarRegionEstimatesLast, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCarRegionEstimatesLast))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCarRegionEstimatesLast),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_LastYearPerSex_region.tex"))
  
  
  
  ##-- UD-Density
  NCarRegionEstimatesLast_UD <- matrix("", ncol = 3, nrow = length(idregionTable))
  row.names(NCarRegionEstimatesLast_UD) <- c(idregionTable)
  colnames(NCarRegionEstimatesLast_UD) <- c("Females","Males","Total")
  
  for(i in 1:length(idregionTable)){
    ##-- FEMALES
    NCarRegionEstimatesLast_UD[idregionTable[i],"Females"] <-
      paste0(round(UDdensityF[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(UDdensityF[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(UDdensityF[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- MALES
    NCarRegionEstimatesLast_UD[idregionTable[i],"Males"] <-
      paste0(round(UDdensityM[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(UDdensityM[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(UDdensityM[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- BOTH SEXES
    NCarRegionEstimatesLast_UD[idregionTable[i],"Total"] <-
      paste0(round(UDdensity[[n.years]]$summary[idregionTable[i],"mean"],digits = 1)," (",
             round(UDdensity[[n.years]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
             round(UDdensity[[n.years]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")")
  }#i
  
  ##-- ADJUST NAMES
  idregion1 <- idregionTable
  idregion1[which(idregion1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idregion1
  
  ##-- print .csv
  write.csv( NCarRegionEstimatesLast_UD,
             file = file.path(working.dir, "tables/N_LastYearPerSex_region_UD.csv"))
  
  ##-- print .tex
  row.names(NCarRegionEstimatesLast_UD) <- c(paste0("\\hspace{0.1cm} ",idregionNOR),"TOTAL")
  print(xtable::xtable(NCarRegionEstimatesLast_UD, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCarRegionEstimatesLast_UD))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCarRegionEstimatesLast_UD),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_LastYearPerSex_region_UD.tex"))
  
  
  
  ## ------     5.1.3. ALL YEARS, PER SEX, PER REGION ------
  
  NCarRegionEstimatesAllSex <- matrix("", ncol = n.years*3, nrow = length(idregionTable)+1)
  row.names(NCarRegionEstimatesAllSex) <- c("", idregionTable)
  colnames(NCarRegionEstimatesAllSex) <- rep(years, each = 3)
  NCarRegionEstimatesAllSex[1, ] <- rep(c("Females","Males","Total"), n.years)
  
  ##-- Fill-in table
  for(t in 1:n.years){
    cols <- which(colnames(NCarRegionEstimatesAllSex) %in% years[t])
    for( i in 1:length(idregionTable)){
      ##-- FEMALES
      colss <- which(NCarRegionEstimatesAllSex[1,cols] %in% "Females")
      NCarRegionEstimatesAllSex[idregionTable[i],cols[colss]] <-
        paste(round(ACdensityF[[t]]$summary[idregionTable[i],"mean"],digits = 1)," (",
              round(ACdensityF[[t]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensityF[[t]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- MALES
      colss <-  which(NCarRegionEstimatesAllSex[1,cols] %in% "Males")
      NCarRegionEstimatesAllSex[idregionTable[i],cols[colss]] <-
        paste(round(ACdensityM[[t]]$summary[idregionTable[i],"mean"],digits = 1)," (",
              round(ACdensityM[[t]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensityM[[t]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- TOTAL
      colss <-  which(NCarRegionEstimatesAllSex[1,cols] %in% "Total")
      NCarRegionEstimatesAllSex[idregionTable[i],cols[colss]] <-
        paste(round(ACdensity[[t]]$summary[idregionTable[i],"mean"],digits = 1)," (",
              round(ACdensity[[t]]$summary[idregionTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensity[[t]]$summary[idregionTable[i],"95%CIHigh"],digits = 0),")",sep="")
    }#i
  }#t
  
  ##-- ADJUST NAMES
  idregion1 <- idregionTable
  idregion1[which(idregion1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimatesAllSex) <- c("", idregion1)
  
  ##-- print .csv
  write.csv(NCarRegionEstimatesAllSex,
            file = file.path(working.dir, "tables/N_AllYearsPerSex_region.csv"))
  
  
  ##-- print .tex
  row.names(NCarRegionEstimatesAllSex) <- c("", paste0("\\hspace{0.1cm} ", idregionNOR), "TOTAL")
  
  print(xtable::xtable(NCarRegionEstimatesAllSex, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCarRegionEstimatesAllSex))),collapse = "")),
        sanitize.text.function = function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCarRegionEstimatesLast),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_AllYearsPerSex_region.tex"))
  
  
  
  ## ------     5.1.4. ALL YEARS, BOTH SEX, PER COUNTY -----
  
  idcounty <- row.names(ACdensity[[1]]$summary)
  
  ##-- Remove Finland, Norway, Russia, Sweden
  idcounty <- idcounty[-which(idcounty %in% c("Finland","Norway","Russia","Sweden","Total"))]
  idcounty <- sort(unique(idcounty))
  
  ##-- Get names of Norwegian counties
  idcountyNOR <- idcounty[-grep("Region",idcounty)]
  idcountyTable <- c(idcountyNOR, "Total")
  
  ##-- Create table to store N estimates (and CI)
  NCountyEstimates <- NCountyMean <- matrix("", ncol = n.years, nrow = length(idcountyTable))
  row.names(NCountyEstimates) <- row.names(NCountyMean) <- c(idcountyTable)
  colnames(NCountyEstimates) <- colnames(NCountyMean) <- years
  
  ##-- Fill-in the table
  for(t in 1:n.years){
    for(i in 1:length(idcountyTable)){
      NCountyEstimates[idcountyTable[i],t] <-
        paste0(round(ACdensity[[t]]$summary[idcountyTable[i],"mean"],digits = 1), " (",
               round(ACdensity[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
               round(ACdensity[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
      
      NCountyMean[idcountyTable[i],t] <- round(ACdensity[[t]]$summary[idcountyTable[i],"mean"])
    }#i
  }#t
  
  ##-- Adjust names
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCountyEstimates) <- idcounty1
  
  ##-- print .csv
  write.csv(NCountyEstimates,
            file = file.path(working.dir, "tables/N_AllYears_county.csv"))
  
  ##-- print .tex
  row.names(NCountyEstimates) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  print(xtable::xtable( NCountyEstimates,
                        type = "latex",
                        align = paste(c("l",rep("c",ncol(NCountyEstimates))),collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list(list(seq(1, nrow(NCountyEstimates), by = 2)), "\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/N_AllYears_county.tex"))
  
  
  
  ## ------     5.1.5. LAST YEAR, PER SEX, PER COUNTY -----
  
  ##-- AC-Density
  NCountyEstimatesLastRegions <- matrix("", ncol = 3, nrow = length(idcountyTable))
  row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
  colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")
  
  for(i in 1:length(idcountyTable)){
    ##-- FEMALES
    NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <-
      paste0(round(ACdensityF[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(ACdensityF[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensityF[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- MALES
    NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <-
      paste0(round(ACdensityM[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(ACdensityM[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensityM[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- BOTH SEXES
    NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <-
      paste0(round(ACdensity[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(ACdensity[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(ACdensity[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }#i
  
  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCountyEstimates) <- idcounty1
  
  ##-- print .csv
  write.csv( NCountyEstimatesLastRegions,
             file = file.path(working.dir, "tables/N_LastYearPerSex_county.csv"))
  
  ##-- print .tex
  row.names(NCountyEstimatesLastRegions) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  print(xtable::xtable(NCountyEstimatesLastRegions, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_LastYearPerSex_county.tex"))
  
  
  # ##-- UD-Density
  # NCountyEstimatesLastRegions_UD <- matrix("", ncol = 3, nrow = length(idcountyTable))
  # row.names(NCountyEstimatesLastRegions_UD) <- c(idcountyTable)
  # colnames(NCountyEstimatesLastRegions_UD) <- c("Females","Males","Total")
  # 
  # for(i in 1:length(idcountyTable)){
  #   ##-- FEMALES
  #   NCountyEstimatesLastRegions_UD[idcountyTable[i],"Females"] <-
  #     paste0(round(UDdensityF[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
  #            round(UDdensityF[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
  #            round(UDdensityF[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  #   
  #   ##-- MALES
  #   NCountyEstimatesLastRegions_UD[idcountyTable[i],"Males"] <-
  #     paste0(round(UDdensityM[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
  #            round(UDdensityM[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
  #            round(UDdensityM[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  #   
  #   ##-- BOTH SEXES
  #   NCountyEstimatesLastRegions_UD[idcountyTable[i],"Total"] <-
  #     paste0(round(UDdensity[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
  #            round(UDdensity[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
  #            round(UDdensity[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  # }#i
  # 
  # ##-- ADJUST NAMES
  # idcounty1 <- idcountyTable
  # idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  # row.names(NCountyEstimates) <- idcounty1
  # 
  # ##-- print .csv
  # write.csv( NCountyEstimatesLastRegions_UD,
  #            file = file.path(working.dir, "tables/N_LastYearPerSex_county_UD.csv"))
  # 
  # ##-- print .tex
  # row.names(NCountyEstimatesLastRegions_UD) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  # print(xtable::xtable(NCountyEstimatesLastRegions_UD, type = "latex",
  #                      align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions_UD))), collapse = "")),
  #       sanitize.text.function=function(x){x},
  #       floating = FALSE,
  #       add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions_UD),by=2)),"\\rowcolor[gray]{.95} "),
  #       file = file.path(working.dir, "tables/N_LastYearPerSex_county_UD.tex"))
  
  
  
  ## ------     5.1.6. ALL YEARS, PER SEX, PER COUNTY ------
  
  NCountyEstimatesAllSex <- matrix("", ncol = n.years*3, nrow = length(idcountyTable)+1)
  row.names(NCountyEstimatesAllSex) <- c("", idcountyTable)
  colnames(NCountyEstimatesAllSex) <- rep(years, each = 3)
  NCountyEstimatesAllSex[1, ] <- rep(c("Females","Males","Total"), n.years)
  
  ##-- Fill-in table
  for(t in 1:n.years){
    cols <- which(colnames(NCountyEstimatesAllSex) %in% years[t])
    for( i in 1:length(idcountyTable)){
      ##-- FEMALES
      colss <- which(NCountyEstimatesAllSex[1,cols] %in% "Females")
      NCountyEstimatesAllSex[idcountyTable[i],cols[colss]] <-
        paste(round(ACdensityF[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(ACdensityF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensityF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- MALES
      colss <- which(NCountyEstimatesAllSex[1,cols] %in% "Males")
      NCountyEstimatesAllSex[idcountyTable[i],cols[colss]] <-
        paste(round(ACdensityM[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(ACdensityM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensityM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- TOTAL
      colss <- which(NCountyEstimatesAllSex[1,cols] %in% "Total")
      NCountyEstimatesAllSex[idcountyTable[i],cols[colss]] <-
        paste(round(ACdensity[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(ACdensity[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(ACdensity[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
    }#i
  }#t
  
  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCountyEstimatesAllSex) <- c("", idcounty1)
  
  ##-- print .csv
  write.csv(NCountyEstimatesAllSex,
            file = file.path(working.dir, "tables/N_AllYearsPerSex_county.csv"))
  
  ##-- print .tex
  row.names(NCountyEstimatesAllSex) <- c("", paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  
  print(xtable::xtable(NCountyEstimatesAllSex, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCountyEstimatesAllSex))),collapse = "")),
        sanitize.text.function = function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/N_AllYearsPerSex_county.tex"))
  
  
  
  ## ------   5.2. NUMBER OF NGS SAMPLES, IDs & DEAD RECOVERIES ------

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
  n.detected_F <- apply(nimDataF$y.alive[ ,1, ], 2, function(x)sum(x>0))
  n.detected_M <- apply(nimDataM$y.alive[ ,1, ], 2, function(x)sum(x>0))

  propDetected <- matrix("", ncol = n.years, nrow = 3)
  row.names(propDetected) <- c("F","M","Total")
  colnames(propDetected) <- years
  for(t in 1:n.years){
    propDetected["F",t] <- getCleanEstimates(n.detected_F[t]/colSums(ACdensityF[[t]]$PosteriorRegions))
    propDetected["M",t] <- getCleanEstimates(n.detected_M[t]/colSums(ACdensityM[[t]]$PosteriorRegions))
    propDetected["Total",t] <- getCleanEstimates((n.detected_F[t]+n.detected_M[t])/
                                                   (colSums(ACdensityF[[t]]$PosteriorRegions)+
                                                      colSums(ACdensityM[[t]]$PosteriorRegions)))
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
  
  
  
  ## ------     5.2.4. PROPORTION OF THE NORWEGIAN POPULATION DETECTED ------
  
  ##-- Extract number of individuals detected
  isDetected <- rbind(nimDataM$y.alive[ ,1, ],nimDataF$y.alive[ ,1, ]) > 0
  
  ##-- Identify individual sex
  isFemale <- resultsSXYZ_MF$sims.list$sex == "F"
  isMale <- resultsSXYZ_MF$sims.list$sex == "M"
  
  ##-- Identify individual status
  isAvail <- resultsSXYZ_MF$sims.list$z == 1 
  isAlive <- resultsSXYZ_MF$sims.list$z == 2 
  
  
  ##-- Calculate % of the Norwegian bear population detected 
  prop <- matrix(NA,3,n.years)
  dimnames(prop) <- list("% individuals" = c("F","M","Total"),
                         "Years" = c(years))
  for(t in 1:n.years){
    prop_F <- prop_M <- prop_tot <- rep(NA,n.mcmc)
    for(iter in 1:n.mcmc){
      country <- countryRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
      isNor <- country %in% 2
      
      ##-- Detected female
      N_F <- sum(isAlive[iter, ,t] & isNor & isFemale)
      N_det_F <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isFemale)
      prop_F[iter] <- N_det_F / N_F 
      
      ##-- Detected male
      N_M <- sum(isAlive[iter, ,t] & isNor & isMale)
      N_det_M <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isMale)
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
            file = file.path(working.dir, "tables/PropDetectedNorway.csv"))
  
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
        file = file.path(working.dir, "tables/PropDetectedNorway.tex"))
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
          file = file.path(working.dir, "tables/PropDetectedNorway.tex"))
  }
    

  
  ## ------     5.2.5. NUMBER OF IDs w/ ACs OUTSIDE NORWAY ------
  
  ##-- Calculate number of individuals alive with their AC in each country each year
  N_det_by_country <- matrix(NA,5,n.years)
  dimnames(N_det_by_country) <- list("Countries" = c("Norway","Sweden","Finland","Russia","Out"),
                                     "Years" = c(years))
  for(t in 1:n.years){
    
    N_fin_F <- N_fin_M <- N_fin <- rep(NA,n.mcmc)
    N_nor_F <- N_nor_M <- N_nor <- rep(NA,n.mcmc)
    N_rus_F <- N_rus_M <- N_rus <- rep(NA,n.mcmc)
    N_swe_F <- N_swe_M <- N_swe <- rep(NA,n.mcmc)
    N_out_F <- N_out_M <- N_out <- rep(NA,n.mcmc)
    
    for(iter in 1:n.mcmc){
      
      country <- countryRaster[raster::cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
      isFin <- country %in% 1
      isNor <- country %in% 2
      isRus <- country %in% 3
      isSwe <- country %in% 4
      
      ##-- Detected individuals
      N_fin_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isFin & isFemale)
      N_fin_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isFin & isMale)
      N_fin[iter] <- N_fin_F[iter] + N_fin_M[iter]
      
      N_nor_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isFemale)
      N_nor_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor & isMale)
      N_nor[iter] <- N_nor_F[iter] + N_nor_M[iter]
      
      N_rus_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isRus & isFemale)
      N_rus_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isRus & isMale)
      N_rus[iter] <- N_rus_F[iter] + N_rus_M[iter]
      
      N_swe_F[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isSwe & isFemale)
      N_swe_M[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isSwe & isMale)
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
  
  }#t
  
  ##-- print .csv
  write.csv(N_det_by_country,
            file = file.path(working.dir, "tables/NumDetectedIds_country.csv"))
  
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
  
  
  ##-- Remove unnecessary objects from memory
  rm(list = c("N_fin_F", "N_fin_M", "N_fin",
              "N_nor_F", "N_nor_M", "N_nor",
              "N_rus_F", "N_rus_M", "N_rus",
              "N_swe_F", "N_swe_M", "N_swe",
              "N_out_F", "N_out_M", "N_out",
              "isDetected"))
  gc(verbose = FALSE)  
  
  
  
  ##----------------------------------------------------------------------------
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
