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
#' @importFrom grDevices adjustcolor dev.off pdf
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
  extraction.res = 5000
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
  load(list.files(file.path(working.dir, "NimbleInFiles/Hunn"), full.names = T)[1])
  nimDataF <- nimData
  nimInitsF <- nimInits
  
  ##-- Males
  load(list.files(file.path(working.dir, "NimbleInFiles/Hann"), full.names = T)[1])
  nimDataM <- nimData
  nimInitsM <- nimInits
  
  ##-- Habitat
  ##-- WE MAY WANT TO ADD WARNINGS HERE ???
  load(file.path( working.dir, "data",
                  paste0("Habitat_bear_", DATE, ".RData")))
  
  ##-- Detectors
  ##-- WE MAY WANT TO ADD WARNINGS HERE ???
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
  
  years <- as.numeric(dimnames(detectors$covariates)[[3]])
  n.years <- length(years) 
  
  
  
  
  ## ------ 2. PROCESS MCMC SAMPLES -----
  message("## Processing model MCMC outputs...")
  
  if(file.exists(file.path( working.dir, "data", paste0("MCMC_bear_", DATE, ".RData")))){
    load(file.path( working.dir, "data", paste0("MCMC_bear_", DATE, ".RData")))
  } else {
    ## ------   2.1. FEMALES -----
    ##-- Compile MCMC bites
    nimOutput_F <- collectMCMCbites( path = file.path(working.dir, "NimbleOutFiles/Hunn"),
                                     burnin = nburnin)
    
    ##-- Traceplots
    grDevices::pdf(file.path(working.dir, "figures/traceplots_F.pdf"))
    plot(nimOutput_F$samples[ ,!is.na(nimOutput_F$samples[[1]][1, ])])
    grDevices::dev.off()
    
    ##-- Process MCMC output
    results_F <- ProcessCodaOutput( nimOutput_F$samples,
                                    params.omit = c("sxy","z"))
    resultsSXYZ_F <- ProcessCodaOutput(nimOutput_F$samples2)
    
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
    nimOutput_M <- collectMCMCbites( path = file.path(working.dir, "NimbleOutFiles/Hann"),
                                     burnin = nburnin)
    
    ##-- Traceplots
    grDevices::pdf(file.path(working.dir, "figures/traceplots_M.pdf"))
    plot(nimOutput_M$samples[ ,!is.na(nimOutput_M$samples[[1]][1, ])])
    dev.off()
    
    ##-- Process MCMC output
    results_M <- ProcessCodaOutput( nimOutput_M$samples,
                                    params.omit = c("sxy","z"))
    resultsSXYZ_M <- ProcessCodaOutput(nimOutput_M$samples2)
    
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
  
  
  
  ## ------ 3. EXTRACT DENSITY -----
  ##-- Names of regions to extract density for
  regions.names <- c("Region 1","Region 2","Region 3","Region 4",
                     "Region 5","Region 6","Region 7","Region 8")
  
  ##-- Remove buffer from the habitat
  habitat.rWthBuffer <- habitat$habitat.rWthBuffer
  habitat.rWthBuffer[habitat.rWthBuffer[] %in% 0] <- NA
  searchedPolygon <- raster::rasterToPolygons( habitat.rWthBuffer,
                                               dissolve = T,
                                               function(x) x == 1)
  
  ##-- Habitat raster with extent used in the model
  habitatPolygon5km <- raster::crop( extraction.raster$Habitat,
                                     habitat$habitat.r)
  
  ##-- Create 5km raster for extraction
  rrRegions <- extraction.raster$Regions
  rrRegions <- raster::mask(rrRegions, habitat$habitat.poly)
  rrRegions <- raster::crop(rrRegions, habitat$habitat.r)
  
  ##-- Get the objects to run the density function
  densityInputRegions <- getDensityInput(
    regions = rrRegions,
    habitat = habitatPolygon5km,
    s = resultsSXYZ_MF$sims.list$sxy,
    plot.check = F)
  
  ##-- Subset to regions of interest
  regionID <- densityInputRegions$regions.rgmx
  row.names(regionID) <- row.names(densityInputRegions$regions.rgmx)
  regionID <- as.matrix(regionID[row.names(regionID) %in% regions.names, ])
  
  ##-- Calculate area of extraction
  #sum(colSums(regionID)>0)
  
  ##-- Select n.iter iterations randomly
  if(dim(densityInputRegions$sy)[1] >= niter){
    iter <- seq(1, dim(densityInputRegions$sy)[1], length.out = niter)
  } else {
    newN <- dim(densityInputRegions$sy)[1]
    warning(paste0( "The number of MCMC samples available (", newN,
                    ") is less than niter = ",  dim(densityInputRegions$sy)[1],
                    ".\nusing niter = ", newN, " instead."))
    iter <- 1:newN
  }
  
  ##-- Calculate density only if necessary
  if(file.exists(file.path(working.dir, "data", paste0("Density_bear_", DATE, ".RData")))){
    
    message("## Loading pre-processed population density...")
    
    load(file.path( working.dir, "data", paste0("Density_bear_", DATE, ".RData")))
    
  } else {
    
    message("## Extracting population density... \n## This might take a while...")
    
    ## ------   1. AC-BASED DENSITY (5km) ------
    ## ------     1.1. MALE & FEMALES -----
    DensityCountriesRegions <- list()
    for(t in 1:n.years){
      DensityCountriesRegions[[t]] <- GetDensity(
        sx = densityInputRegions$sx[iter, ,t],
        sy = densityInputRegions$sy[iter, ,t],
        z = resultsSXYZ_MF$sims.list$z[iter, ,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = regionID,
        returnPosteriorCells = F)
    }#t
    
    
    
    ## ------     1.2. MALE -----
    IDMales <- which(resultsSXYZ_MF$sims.list$sex == "M")
    
    DensityCountriesRegionsM <- list()
    for(t in 1:n.years){
      DensityCountriesRegionsM[[t]] <- GetDensity(
        sx = densityInputRegions$sx[iter,IDMales,t],
        sy =  densityInputRegions$sy[iter,IDMales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDMales,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = regionID,
        returnPosteriorCells = F)
    }#t
    
    
    
    ## ------     1.3. FEMALE -----
    IDFemales <- which(resultsSXYZ_MF$sims.list$sex == "F")
    
    DensityCountriesRegionsF <- list()
    for(t in 1:n.years){
      DensityCountriesRegionsF[[t]] <- GetDensity(
        sx = densityInputRegions$sx[iter,IDFemales,t],
        sy = densityInputRegions$sy[iter,IDFemales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDFemales,t],
        IDmx = densityInputRegions$habitat.id,
        aliveStates = 2,
        regionID = regionID,
        returnPosteriorCells = F)
    }
    
    
    
    ## ------   2. UD-BASED DENSITY (5km) ------
    ##--  Combine male and female sigma
    sigma <- do.call(cbind,lapply(resultsSXYZ_MF$sims.list$sex,
                                  function(x){
                                    if(x == "M"){
                                      resultsSXYZ_MF$sims.list$sigma[ ,"M"]
                                    } else {
                                      resultsSXYZ_MF$sims.list$sigma[ ,"F"]
                                    }
                                  }))
    
    ##-- Rescale sigma to the raster resolution
    sigma <- sigma/raster::res(rrRegions)[1]
    
    
    
    ## ------     2.1. MALE -----
    IDMales <- which(resultsSXYZ_MF$sims.list$sex=="M")
    
    spaceUSEDM <- list()
    for(t in 1:n.years){
      spaceUSEDM[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[iter,IDMales,t],
        sy = densityInputRegions$sy[iter,IDMales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDMales,t],
        sigma = sigma[iter,IDMales],
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = T)
    }#t
    
    
    
    ## ------     2.2. FEMALE -----
    IDFemales <- which(resultsSXYZ_MF$sims.list$sex=="F")
    
    spaceUSEDF <- list()
    for(t in 1:n.years){
      spaceUSEDF[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[iter,IDFemales,t],
        sy = densityInputRegions$sy[iter,IDFemales,t],
        z = resultsSXYZ_MF$sims.list$z[iter,IDFemales,t],
        sigma = sigma[iter,IDFemales],
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = T)
    }#t
    
    
    
    ## ------     2.3. MALE & FEMALES -----
    spaceUSED <- list()
    for(t in 1:n.years){
      spaceUSED[[t]] <- GetSpaceUse(
        sx = densityInputRegions$sx[iter, ,t],
        sy = densityInputRegions$sy[iter, ,t],
        z = resultsSXYZ_MF$sims.list$z[iter, ,t],
        sigma = sigma[iter, ],
        habitatxy = densityInputRegions$habitat.xy,
        aliveStates = 2,
        regionID = regionID,
        display_progress = T,
        returnPosteriorCells = T)
    }#t
    
    
    
    ## ------   3. SAVE DENSITY OBJECTS ------
    save( DensityCountriesRegions,
          DensityCountriesRegionsF,
          DensityCountriesRegionsM,
          spaceUSED,
          spaceUSEDF,
          spaceUSEDM,
          file = file.path( working.dir, "data", paste0("Density_bear_", DATE, ".RData")))
  }
  
  
  
  ## ------ 4. FIGURES -----
  
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
    input = densityInputRegions,
    estimates = DensityCountriesRegions,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES[1,],
    type = c("all"),# "last.year", "time.series"),
    path = file.path(working.dir, "figures"),
    name = "AC_Density")
  
  ##-- UD-density maps
  plotDensityMaps( 
    input = densityInputRegions,
    estimates = spaceUSED,
    unit = 100,
    mask = rrCombined,
    background = COUNTRIES[1,],
    type = c("all"),
    path = file.path(working.dir, "figures"),
    name = "UD_Density")
  
  
  
  
  ## ------   4.2. ABUNDANCE TIME SERIES ------
  message("## Plotting abundance...")
  
  ##-- Extract number of individuals detected
  isDetected <- rbind(nimDataM$y.alive[ ,1, ],nimDataF$y.alive[ ,1, ]) > 0
  n.detected <- apply(isDetected, 2, sum)
  
  ##-- Plot parameters
  diffSex <- 0.2
  colSex <- adjustcolor(c("firebrick3", "deepskyblue2", "black"),0.5)
  colCause  <- adjustcolor( c("#E69F00","#009E73"), 0.5)
  
  ##-- Plot N  
  # pdf(file = file.path(working.dir, "figures/Abundance_TimeSeries.pdf"),
  #     width = 12, height = 8.5)
  grDevices::png(filename = file.path(working.dir, "figures/Abundance_TimeSeries.png"),
                 width = 12, height = 8.5, units = "in", pointsize = 12,
                 res = 300, bg = NA)
  
  par(mar = c(5,5,1,1))
  plot(-1000,
       xlim = c(0.5, n.years + 0.5),
       ylim = c(0,180),
       xlab = "", ylab = paste("Number of bears"),
       xaxt = "n", axes = F, cex.lab = 1.6)
  axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
  axis(2, at = seq(0,170,20), labels = seq(0,170,20), cex.axis = 1.6)
  abline(v = (0:n.years)+0.5, lty = 2)
  abline(h = seq(0,170, by = 10), lty = 2, col = "gray90")
  
  for(t in 1:n.years){
    ##-- FEMALES
    plotQuantiles(x = colSums(DensityCountriesRegionsF[[t]]$PosteriorAllRegions),
                  at = t - diffSex,
                  width = 0.18,
                  col = colSex[1])
    
    ##-- MALES
    plotQuantiles(x = colSums(DensityCountriesRegionsM[[t]]$PosteriorAllRegions),
                  at = t + diffSex,
                  width = 0.18,
                  col = colSex[2])
    ##-- TOTAL
    plotQuantiles(x = colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions),
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
  
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 3.5, col = colSex)
  points(x = xx[1:3], y = yy[1:3],  pch = 15, cex = 1.5, col = colSex)
  text(x = xx + 0.1, y = yy-1, labels = labs, cex = 1.2, pos = 4)
  
  polygon(x = c(xx[4]-0.3,xx[4]+0.1,xx[4]+0.1,xx[4]-0.3),
          y = c(yy[4]-1,yy[4]-1,yy[4]+1,yy[4]+1),
          col = "goldenrod1",
          border = F)
  
  dev.off()
  
  
  
  ## ------ 5. TABLES -----
  gc()
  ## ------   5.1. ABUNDANCE -----
  ## ------     5.1.1. ALL YEARS, BOTH SEX -----
  idcounty <- row.names(DensityCountriesRegions[[1]]$summary)
  
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
        paste0(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1), " (",
               round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
               round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
      
      NCarRegionMean[idcountyTable[i],t] <- round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"])
    }#i
  }#t
  
  ##-- Adjust names
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idcounty1
  
  ##-- print .csv
  write.csv(NCarRegionEstimates,
            file = file.path(working.dir, "tables/Abundance_AllYears.csv"))
  
  ##-- print .tex
  row.names(NCarRegionEstimates) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  print(xtable::xtable( NCarRegionEstimates,
                        type = "latex",
                        align = paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list(list(seq(1, nrow(NCarRegionEstimates), by = 2)), "\\rowcolor[gray]{.96} "),
        file = file.path(working.dir, "tables/Abundance_AllYears.tex"))
  
  
  
  
  ## ------     5.1.2. LAST YEAR N PER SEX PER COUNTY -----
  NCountyEstimatesLastRegions <- matrix("", ncol = 3, nrow = length(idcountyTable))
  row.names(NCountyEstimatesLastRegions) <- c(idcountyTable)
  colnames(NCountyEstimatesLastRegions) <- c("Females","Males","Total")
  
  for(i in 1:length(idcountyTable)){
    ##-- FEMALES
    NCountyEstimatesLastRegions[idcountyTable[i],"Females"] <-
      paste0(round(DensityCountriesRegionsF[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(DensityCountriesRegionsF[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(DensityCountriesRegionsF[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- MALES
    NCountyEstimatesLastRegions[idcountyTable[i],"Males"] <-
      paste0(round(DensityCountriesRegionsM[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(DensityCountriesRegionsM[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(DensityCountriesRegionsM[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- BOTH SEXES
    NCountyEstimatesLastRegions[idcountyTable[i],"Total"] <-
      paste0(round(DensityCountriesRegions[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(DensityCountriesRegions[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(DensityCountriesRegions[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }#i
  
  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idcounty1
  
  ##-- print .csv
  write.csv( NCountyEstimatesLastRegions,
             file = file.path(working.dir, "tables/Abundance_LastYearPerSex.csv"))
  
  ##-- print .tex
  row.names(NCountyEstimatesLastRegions) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  print(xtable::xtable(NCountyEstimatesLastRegions, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/Abundance_LastYearPerSex.tex"))
  
  
  
  ##-- UD-Density
  NCountyEstimatesLastRegions_UD <- matrix("", ncol = 3, nrow = length(idcountyTable))
  row.names(NCountyEstimatesLastRegions_UD) <- c(idcountyTable)
  colnames(NCountyEstimatesLastRegions_UD) <- c("Females","Males","Total")
  
  for(i in 1:length(idcountyTable)){
    ##-- FEMALES
    NCountyEstimatesLastRegions_UD[idcountyTable[i],"Females"] <-
      paste0(round(spaceUSEDF[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(spaceUSEDF[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(spaceUSEDF[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- MALES
    NCountyEstimatesLastRegions_UD[idcountyTable[i],"Males"] <-
      paste0(round(spaceUSEDM[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(spaceUSEDM[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(spaceUSEDM[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
    
    ##-- BOTH SEXES
    NCountyEstimatesLastRegions_UD[idcountyTable[i],"Total"] <-
      paste0(round(spaceUSED[[n.years]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
             round(spaceUSED[[n.years]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
             round(spaceUSED[[n.years]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")")
  }#i
  
  
  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCarRegionEstimates) <- idcounty1
  
  ##-- print .csv
  write.csv( NCountyEstimatesLastRegions_UD,
             file = file.path(working.dir, "tables/Abundance_LastYearPerSex_UD.csv"))
  
  ##-- print .tex
  row.names(NCountyEstimatesLastRegions_UD) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  print(xtable::xtable(NCountyEstimatesLastRegions_UD, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions_UD))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions_UD),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/Abundance_LastYearPerSex_UD.tex"))
  
  
  
  
  ## ------     5.1.3. ALL YEARS N PER SEX PER COUNTY ------
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
        paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- MALES
      colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Males")
      NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <-
        paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
      
      ##-- TOTAL
      colss <-  which(NCountyEstimatesAllSexRegions[1,cols] %in% "Total")
      NCountyEstimatesAllSexRegions[idcountyTable[i],cols[colss]] <-
        paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
    }#i
  }#t
  
  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCountyEstimatesAllSexRegions) <- c("", idcounty1)
  
  ##-- print .csv
  write.csv(NCountyEstimatesAllSexRegions,
            file = file.path(working.dir, "tables/Abundance_AllYearsPerSex.csv"))
  
  
  ##-- print .tex
  row.names(NCountyEstimatesAllSexRegions) <- c("", paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  
  print(xtable::xtable(NCountyEstimatesAllSexRegions, type = "latex",
                       align = paste(c("l",rep("c",ncol(NCountyEstimatesAllSexRegions))),collapse = "")),
        sanitize.text.function = function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(working.dir, "tables/Abundance_AllYearsPerSex.tex"))
  
  
  
  
  ## ------ 6. OUTPUT -----
  out$YEARS <- years
  
  return(out)
}
