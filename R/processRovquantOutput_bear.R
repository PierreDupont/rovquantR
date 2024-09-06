#' @title RovQuant OPSCR bear output processing
#' 
#' @description
#' \code{processRovquantOutput_bear} calls a custom Rmarkdown template that combines 
#' and processes MCMC outputs from NIMBLE models and produces figures,
#' tables and rasters of interest (e.g. population density maps)
#' 
#' @param data_dir A \code{path}
#' @param working_dir A \code{path}
#' @param nburnin An \code{integer} denoting the number of iterations to be removed from each MCMC as burnin.
#' @param niter An \code{integer} denoting the number of MCMC iterations to be used for density extraction.
#' @param extraction.res A \code{integer} denoting the raster resolution for density extraction.
#' @param years A \code{path}
#' @param sex An \code{integer} denoting the number of iterations to be removed from each MCMC as burnin.
#' @param plot.check An \code{integer} denoting the number of MCMC iterations to be used for density extraction.
#' @param print.report A \code{integer} denoting the raster resolution for density extraction.
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
#' @importFrom grDevices adjustcolor 
#' @importFrom stars st_as_stars
#' 
#' @rdname processRovquantOutput_bear
#' @export
processRovquantOutput_bear <- function(
  ##-- paths
  data_dir = "./Data",
  working_dir = NULL,

  ##-- MCMC processing
  nburnin = 0,

  ##-- Density extraction
  niter = 100,
  extraction.res = 5000,

  ##-- miscellanious
  plot.check = FALSE,
  print.report = TRUE
)
{
  ## ------ 0. BASIC SET-UP ------
  if(is.null(working_dir)){working_dir <- getwd()}


  
  ## ------ 1. LOAD NECESSARY INPUTS -----
  ##-- Females
  load(list.files(file.path(working_dir, "NimbleInFiles/Hunn"), full.names = T)[1])
  nimDataF <- nimData
  nimInitsF <- nimInits

  ##-- Males
  load(list.files(file.path(working_dir, "NimbleInFiles/Hann"), full.names = T)[1])
  nimDataM <- nimData
  nimInitsM <- nimInits

  ##-- Habitat
  load(file.path(working_dir, "data/Habitat.RData"))

  ##-- Detectors
  load(file.path(working_dir, "data/Detectors.RData"))

  ##-- Habitat Rasters
  if(extraction.res <= 1000){
    extraction.raster <- habitatRasterResolution$'1km'
  } else {
    if(extraction.res <= 2000){
      extraction.raster <- habitatRasterResolution$'2km'
    } else {
      if(extraction.res <= 5000){
        extraction.raster <- habitatRasterResolution$'5km'
      } else {
        if(extraction.res <= 10000){
          extraction.raster <- habitatRasterResolution$'10km'
        } else {
          extraction.raster <- habitatRasters
        }}}}



  ## ------ 2. PROCESS MCMC SAMPLES -----
  message("## Processing model MCMC outputs...")

  if(file.exists(file.path(working_dir, "data/MCMC.RData"))){
    load(file.path(working_dir, "data/MCMC.RData"))
  } else {
    ## ------   2.1. FEMALES -----
    ##-- Compile MCMC bites
    nimOutput_F <- collectMCMCbites( path = file.path(working_dir, "NimbleOutFiles/Hunn"),
                                     burnin = nburnin)

    ##-- Traceplots
    pdf(file.path(working_dir, "figures/traceplots_F.pdf"))
    plot(nimOutput_F$samples[ ,!is.na(nimOutput_F$samples[[1]][1, ])])
    dev.off()

    ##-- Process MCMC output
    results_F <- ProcessCodaOutput( nimOutput_F$samples,
                                       params.omit = c("sxy","z"))
    resultsSXYZ_F <- ProcessCodaOutput(nimOutput_F$samples2)

    ##-- Rescale sxy to the original coordinate system
    dimnames(resultsSXYZ_F$sims.list$sxy)[[3]] <- c("x","y")
    resultsSXYZ_F$sims.list$sxy <- scaleCoordsToHabitatGrid(
      coordsData = resultsSXYZ_F$sims.list$sxy,
      coordsHabitatGridCenter = habitat$habitat.xy,
      scaleToGrid = FALSE)$coordsDataScaled

    ##-- RESCALE sigma AND tau TO THE ORIGINAL COORDINATE SYSTEM
    results_F$sims.list$sigma <- results_F$sims.list$sigma * res(habitat$habitat.r)[1]
    results_F$sims.list$tau <- results_F$sims.list$tau * res(habitat$habitat.r)[1]



    ## ------   2.2. MALES -----
    ##-- Compile MCMC bites
    nimOutput_M <- collectMCMCbites( path = file.path(working_dir, "NimbleOutFiles/Hann"),
                                     burnin = nburnin)

    ##-- Traceplots
    pdf(file.path(working_dir, "figures/traceplots_M.pdf"))
    plot(nimOutput_F$samples[ ,!is.na(nimOutput_F$samples[[1]][1, ])])
    dev.off()

    ##-- Process MCMC output
    results_M <- ProcessCodaOutput( nimOutput_M$samples,
                                       params.omit = c("sxy","z"))
    resultsSXYZ_M <- ProcessCodaOutput(nimOutput_M$samples2)

    ##-- RESCALE SXY TO THE ORIGINAL COORDINATE SYSTEM
    dimnames(resultsSXYZ_M$sims.list$sxy)[[3]] <- c("x","y")
    resultsSXYZ_M$sims.list$sxy <- scaleCoordsToHabitatGrid(
      coordsData = resultsSXYZ_M$sims.list$sxy,
      coordsHabitatGridCenter = habitat$habitat.df,
      scaleToGrid = FALSE)$coordsDataScaled

    ##-- RESCALE sigma AND tau TO THE ORIGINAL COORDINATE SYSTEM
    results_M$sims.list$sigma <- results_M$sims.list$sigma * res(habitat$habitat.r)[1]
    results_M$sims.list$tau <- results_M$sims.list$tau * res(habitat$habitat.r)[1]



    ## ------   2.3. COMBINE MALES & FEMALES -----
    resultsSXYZ_MF <- resultsSXYZ_M

    ##-- Get minimum number of iterations between model F and M
    minIter <- min(dim(resultsSXYZ_F$sims.list$sxy)[1],
                   dim(resultsSXYZ_M$sims.list$sxy)[1])

    ##-- sxy
    resultsSXYZ_MF$sims.list$sxy <- abind(resultsSXYZ_M$sims.list$sxy[1:minIter, , , ],
                                          resultsSXYZ_F$sims.list$sxy[1:minIter, , , ],
                                          along = 2)
    dimnames(resultsSXYZ_MF$sims.list$sxy)[[3]] <- c("x","y")

    ##-- z
    resultsSXYZ_MF$sims.list$z <- abind(resultsSXYZ_M$sims.list$z[1:minIter, , ],
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
          file = file.path(working_dir, "data/MCMC.RData"))
  }



  ## ------ 3. EXTRACT DENSITY -----
  message("## Extracting population density...")

  ##-- Names of regions to extract density for
  regions.names <- c("Region 1","Region 2","Region 3","Region 4",
                     "Region 5","Region 6","Region 7","Region 8")

  ##-- Remove buffer from the habitat
  habitat.rWthBuffer <- habitat$habitat.rWthBuffer
  habitat.rWthBuffer[habitat.rWthBuffer %in% 0] <- NA
  searchedPolygon <- raster::rasterToPolygons( habitat.rWthBuffer,
                                       dissolve = T,
                                       function(x) x == 1)

  ##-- Habitat raster with extent used in the model
  habitatPolygon5km <- raster::crop( habitatRasterResolution$`5km`[["Habitat"]],
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
  if(dim(densityInputRegions$sy)[1] > niter){
    iter <- seq(1, dim(densityInputRegions$sy)[1], length.out = niter)
  } else {
    iter <- 1:dim(densityInputRegions$sy)[1]
  }


  if(file.exists(file.path(working_dir, "results/Density.RData"))){
    load(file.path(working_dir, "results/Density.RData"))
  } else {
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
      spaceUSEDF[[t]] <- GetSpaceUse_PD(
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
      spaceUSEDF[[t]] <- GetSpaceUse_PD(
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
      spaceUSED[[t]] <- GetSpaceUse_PD(
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
    save(DensityCountriesRegions,
         DensityCountriesRegionsF,
         DensityCountriesRegionsM,
         spaceUSED,
         spaceUSEDF,
         spaceUSEDM,
         file = file.path( working_dir, "results/Density.RData"))
  }



  ##----------------------------------------------------------------------------
  ## ------   1. PLOTS -----

  y.deadF <- nimDataF$y.dead
  y.deadM <- nimDataM$y.dead
  y.deadMF <- rbind(y.deadM, y.deadF)

  ##-- Create 5km raster for plotting
  rrNorway <- habitatRasterResolution$`5km`[["Countries"]]
  rrNorway[!rrNorway[] %in% 2] <- NA
  rrNorway[rrNorway[] %in% 2] <- 1
  rrNorway <- raster::crop(rrNorway, habitat$habitat.r)

  ##-- Extract number of individuals detected
  isDetected <- rbind(nimDataM$y.alive[ ,1, ],nimDataF$y.alive[ ,1, ]) > 0
  n.detected <- apply(isDetected, 2, sum)

  ##-- Identify individual sex
  isFemale <- resultsSXYZ_MF$sims.list$sex == "F"
  isMale <- resultsSXYZ_MF$sims.list$sex == "M"

  ##-- Identify individual status
  isAvail <- resultsSXYZ_MF$sims.list$z == 1
  isAlive <- resultsSXYZ_MF$sims.list$z == 2
  isDead <- resultsSXYZ_MF$sims.list$z >= 3

  ##-- Get number of iterations
  n.iter <- dim(resultsSXYZ_MF$sims.list$z)[1]

  ##-- Plot parameters
  diffSex <- 0.2
  colSex <- adjustcolor(c("firebrick3", "deepskyblue2", "black"),0.5)
  colCause  <- adjustcolor( c("#E69F00","#009E73"), 0.5)




  ## ------     1.1. DENSITY MAPS -----
  message("## Plotting population density maps...")

  ## ------       1.1.0. DENSITY PLOT SET-UP ------
  ##-- Convert densities from 25km2 (5*5raster) to 100km2
  UD100km2 <- lapply(spaceUSED, function(x) x$MeanCell * 4)

  AC100km2 <- lapply(DensityCountriesRegions, function(x) x$MeanCell * 4)
  AC_F100km2 <- lapply(DensityCountriesRegionsF, function(x) x$MeanCell * 4)
  AC_M100km2 <- lapply(DensityCountriesRegionsM, function(x) x$MeanCell * 4)

  rrCombined <- rrRegions + rrNorway

  ##-- Crop density maps to Norway
  ACCropped <- UDCropped <- UDCropped_F <- UDCropped_M <- list()
  for(t in 1:length(years)){
    ACCropped[[t]] <- densityInputRegions$regions.r
    ACCropped[[t]][] <- NA
    ACCropped[[t]][!is.na(densityInputRegions$regions.r[])] <- AC100km2[[t]]
    ACCropped[[t]][is.na(rrCombined[])] <- NA

    UDCropped[[t]] <- densityInputRegions$regions.r
    UDCropped[[t]][] <- NA
    UDCropped[[t]][!is.na(densityInputRegions$regions.r[])] <- UD100km2[[t]]
    UDCropped[[t]][is.na(rrCombined[])] <- NA

    UDCropped_F[[t]] <- densityInputRegions$regions.r
    UDCropped_F[[t]][] <- NA
    UDCropped_F[[t]][!is.na(densityInputRegions$regions.r[])] <- UD_F100km2[[t]]
    UDCropped_F[[t]][is.na(rrCombined[])] <- NA

    UDCropped_M[[t]] <- densityInputRegions$regions.r
    UDCropped_M[[t]][] <- NA
    UDCropped_M[[t]][!is.na(densityInputRegions$regions.r[])] <- UD_M100km2[[t]]
    UDCropped_M[[t]][is.na(rrCombined[])] <- NA
    raster::proj4string(ACCropped[[t]]) <- proj4string(ACCropped[[t]]) <- st_crs(habitat$habitat.poly)
  }#t



  ## ------       1.1.1. AC-Density ------
  pdf(file = file.path(working_dir, "figures/AC_DensityMaps.pdf"),
      width = 12, height = 8)

  ##-- Set color scale
  max <- max(unlist(lapply(ACCropped, function(x) max(x[], na.rm = T))))
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
  col <- colfunc(100)

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
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
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
      plot( UDCropped[[t]],
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
    ##-- Export rasters
    writeRaster(ACCropped[[t]],
                file.path(WDRasters,
                          paste0("AC5kmRaster100km2_classic", years[t],".tif")),
                overwrite = TRUE)
  }#t
  dev.off()



  ## ------       1.1.2. UD-Density ------
  pdf(file = file.path(working_dir, "figures/UD_DensityMaps.pdf"),
      width = 12, height = 8)

  ##-- Set color scale
  max <- max(unlist(lapply(UDCropped, function(x) max(x[], na.rm = T))))
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  colfunc <- colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
  col <- colfunc(100)

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
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
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
      plot( UDCropped[[t]],
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
    ##-- Export rasters
    writeRaster(UDCropped[[t]],
                file.path(WDRasters,
                          paste0("UD5kmRaster100km2_classic", years[t],".tif")),
                overwrite = TRUE)
  }#t
  dev.off()




  ## ------       1.1.3. LAST YEAR SUMMARY ------
  pdf(file = file.path(working_dir, "figures/UD_DensityMaps_LastYear.pdf"),
      width = 8, height = 8)

  t <- length(years)

  par(mar = c(0,0,0,0))
  plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
  image( UDCropped[[t]], add = T, breaks = c(cuts, max(cuts)+1000), col = col, legend = F)
  plot(COUNTRIESsimpFig[1, ], border = grey(0.4), col = NA, add = TRUE)
  mtext(text = years[t], side = 1, -25, adj = 0.25, cex = 3, font = 2)
  segments(x0 = 800000, x1 = 800000,
           y0 = 6650000, y1 = 6650000 + 500000,
           col = grey(0.3), lwd = 4, lend = 2)
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




  ##----------------------------------------------------------------------------
  ## ------     1.2. ABUNDANCE ------
  message("## Plotting abundance...")

  ## ------       1.2.2 BARS ------
  pdf(file = file.path(working_dir, "figures/Abundance_TimeSeries.pdf"),
      width = 12, height = 8.5)
  par(mar = c(5,5,1,1))
  plot(-1000,
       xlim = c(0.5, n.years + 0.5), ylim = c(0,180),
       xlab = "", ylab = paste("Number of bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
  axis(2, at = seq(0,170,20), labels = seq(0,170,20), cex.axis = 1.6)
  abline(v = (0:n.years)+0.5, lty = 2)
  abline(h = seq(0,170, by = 10), lty = 2, col = "gray60")

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

    # ##-- Add mean estimates
    # text(round(round(DensityCountriesRegions[[t]]$summary["Total","mean"], digits = 1)),
    #      x = t,
    #      y = DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"] + 6,
    #      cex = 1.5)
  }#t
  box()

  ##-- legend
  par(xpd = TRUE)
  xx <- c(6.1,7.6,8.8,10)
  yy <- c(5,5,5,5)
  labs <- c("Females", "Males", "Total", "Detected")
  polygon(x = c(5.6,11.6,11.6,5.6),
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



  ## ------       1.2.2 BARS per REGION ------
  pdf(file = file.path(working_dir, "figures/Abundance_TimeSeries_regions.pdf"),
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
      plotQuantiles(x = DensityCountriesRegions[[t]]$PosteriorRegions[cc-1, ],
                    at = t,
                    width = 0.4,
                    col = colSex[3])

      ##-- FEMALES
      plotQuantiles(x = DensityCountriesRegionsF[[t]]$PosteriorRegions[cc-1, ],
                    at = t - diffSex,
                    width = 0.18,
                    col = colSex[1])

      ##-- MALES
      plotQuantiles(x = DensityCountriesRegionsM[[t]]$PosteriorRegions[cc-1, ],
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
    plot(st_geometry(REGIONS), border = grey(0.5), col = grey(0.5), lwd = 0.1)
    plot(st_geometry(REGIONS[REGIONS$region == cc, ]),
         add = T, col = adjustcolor("red",0.5), border = "red")
  }#c
  dev.off()




  ##----------------------------------------------------------------------------
  ## ------     1.3. VITAL RATES ------
  message("## Plotting vital rates...")

  ## ------       1.3.1. SURVIVAL ------
  pdf(file = file.path(working_dir, "figures/SurvivalBars_classic.pdf"),
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

  ##-- Calculate mortality from estimated hazard rates (mhH and mhW) if necessary
  if(any(grep("mhH",names(results_F$sims.list)))){
    mhH1 <- exp(results_F$sims.list$mhH)
    mhW1 <- exp(results_F$sims.list$mhW)
    results_F$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
    results_F$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
    results_F$sims.list$phi <- 1 - results_F$sims.list$h - results_F$sims.list$w

    mhH1 <- exp(results_M$sims.list$mhH)
    mhW1 <- exp(results_M$sims.list$mhW)
    results_M$sims.list$h <- (1-exp(-(mhH1+mhW1))) * (mhH1/(mhH1+mhW1))
    results_M$sims.list$w <- (1-exp(-(mhH1+mhW1))) * (mhW1/(mhH1+mhW1))
    results_M$sims.list$phi <- 1 - results_M$sims.list$h - results_M$sims.list$w
  }

  for(t in 1:(n.years-1)){
    ##-- FEMALES
    plotQuantiles(x = results_F$sims.list$phi[ ,t],
                  at = t - diffSex,
                  col = colSex[1])

    ##-- MALES
    plotQuantiles(x = results_M$sims.list$phi[ ,t],
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
  ##-- Get posterior z and sex
  z_F <- resultsSXYZ_MF$sims.list$z[ ,isFemale, ]
  z_M <- resultsSXYZ_MF$sims.list$z[ ,isMale, ]

  ##-- Calculate mortality from estimated hazard rates (mhH and mhW) if necessary
  if(!any(grep("mhH",names(results_F$sims.list)))){
    ##-- FEMALE
    isDead_F <- apply((z_F[ , ,1:(n.years-1)] == 2)*(z_F[ , ,2:n.years]%in%c(3,4)), c(1,3), sum)
    wasAlive_F <- apply(z_F[ , ,1:(n.years-1)] == 2, c(1,3), sum)
    MortalityAll <- isDead_F / wasAlive_F
    results_F$sims.list$h <- sapply(1:(n.years-1), function(t)sum(y.deadF[ ,t+1]) / wasAlive_F[ ,t])
    results_F$sims.list$w <- MortalityAll - results_F$sims.list$h

    ##-- MALE
    isDead_M <- apply((z_M[ , ,1:(n.years-1)] == 2)*(z_M[ , ,2:n.years]%in%c(3,4)), c(1,3), sum)
    wasAlive_M <- apply(z_M[ , ,1:(n.years-1)] == 2, c(1,3), sum)
    MortalityAll <- isDead_M / wasAlive_M
    results_M$sims.list$h <- sapply(1:(n.years-1), function(t)sum(y.deadM[ ,t+1]) / wasAlive_M[ ,t])
    results_M$sims.list$w  <- MortalityAll - results_M$sims.list$h
  }


  ##-- Plot bars
  pdf(file = file.path(working_dir, "figures", paste0("MortalityBars_classic.pdf")),
      width = 10, height = 9)

  nf <- layout(rbind(c(3,7,6),
                     c(4,1,6),
                     c(5,2,6)),
               widths = c(0.15,1,0.35),
               heights = c(0.15,1,1))

  for(s in 1:2){
    if(s == 1){results <- results_F$sims.list} else {results <- results_M$sims.list}
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
  pdf(file = file.path(working_dir, "figures", paste0("NbRecruitsPerCapita_classic.pdf")),
      width = 10, height = 8)

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
    n.recruit <- apply(resultsSXYZ_MF$sims.list$z[ , ,c(t,t+1)], 1, function(x) sum(x[ ,1]%in%1 & x[ ,2]%in%2))
    ##-- Number of reproducing ids at t ==> ids with state 2 at t
    alivetminus1 <- apply(resultsSXYZ_MF$sims.list$z[ , ,t], 1, function(x)sum(x %in% 2))
    ##-- Plot quantiles
    plotQuantiles(x = n.recruit/alivetminus1,
                  at = t,
                  width = 0.3,
                  col = colCause[1])
  }#t

  dev.off()




  ## ------         1.3.3.2. PLOT NUMBER OF RECRUITS (GLOBAL) ----
  N_recruit <- N_recruit_F <- N_recruit_M <- matrix(NA,n.iter,n.years-1)
  for(t in 1:(n.years-1)){
    for(iter in 1:n.iter){
      N_recruit_F[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & isFemale)
      N_recruit_M[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & isMale)
      N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]
    }#iter
    print(t)
  }#t

  pdf(file = file.path(working_dir, "figures", paste0("NumRecruitTotal_classic.pdf")),
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
    for(iter in 1:dim(resultsSXYZ_MF$sims.list$z)[1]){
      isOut <- is.na(norRaster[cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t+1])])
      wasOut <- is.na(norRaster[cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])])

      N_recruit_F[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & (!isOut) & isFemale)
      N_recruit_M[iter,t] <- sum(isAvail[iter, ,t] & isAlive[iter, ,t+1] & (!isOut) & isMale)
      N_recruit[iter,t] <- N_recruit_F[iter,t] + N_recruit_M[iter,t]

      N_immig_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & wasOut & (!isOut) & isFemale)
      N_immig_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & wasOut & (!isOut) & isMale)
      N_immig[iter,t] <- N_immig_F[iter,t] + N_immig_M[iter,t]

      N_emig_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & isOut & isFemale)
      N_emig_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & isOut & isMale)
      N_emig[iter,t] <- N_emig_F[iter,t] + N_emig_M[iter,t]

      N_surv_F[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & (!isOut) & isFemale)
      N_surv_M[iter,t] <- sum(isAlive[iter, ,t] & isAlive[iter, ,t+1] & (!wasOut) & (!isOut) & isMale)
      N_surv[iter,t] <- N_surv_F[iter,t] + N_surv_M[iter,t]
    }#iter
    print(t)
  }#t



  ## ------       1.4.1. PLOT RECRUITS -----
  pdf(file = file.path(working_dir, "figures", paste0("NumRecruits_bars_classic.pdf")),
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
  pdf(file = file.path(working_dir, "figures", paste0("NumSurvival_bars_classic.pdf")),
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
  pdf(file = file.path(working_dir, "figures", paste0("NumImmigrants_bars_classic.pdf")),
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
  pdf(file = file.path(working_dir, "figures", paste0("NumEmigrants_bars_classic.pdf")),
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
  pdf(file = file.path(working_dir, "figures", paste0("NumFluxes_bars_classic.pdf")),
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
  pdf(file = file.path(working_dir, "figures", paste0("p0_mod_classic.pdf")),
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
      plotQuantiles( x = results_F$sims.list$p0[ ,COUNTIES_s$index == c,t],
                     at = t - diffSex,
                     col = colSex[1])

      plotQuantiles( x = results_M$sims.list$p0[ ,COUNTIES_s$index == c,t],
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
    plot(st_geometry(COUNTIES_s), border = grey(0.5), col = grey(0.5), lwd = 0.1)
    plot(st_geometry(COUNTIES_s[COUNTIES_s$index == c, ]),
         add = T, col = adjustcolor("red",0.5), border = "red")
    text(as_Spatial(COUNTIES_s[COUNTIES_s$index == c, ]),
         labels = COUNTIES_s$Name[COUNTIES_s$index == c],
         col = "white")
  }#c
  dev.off()




  ## ------       1.5.2. p0 maps ------
  pdf(file = file.path(working_dir, "figures", paste0("p0_map_mod_classic.pdf")),
      width = 10, height = 6)
  for(t in 1:n.years){
    par(mfrow = c(1,2))

    ##-- FEMALE
    myDetectors$main.detector.sp$p0_F <-
      ilogit(logit(results_F$mean$p0[nimConstants$county,t]) +
               results_F$mean$betaDet[1] * nimDataF$detCovs[ ,1,t] +
               results_F$mean$betaDet[2] * nimDataF$detCovs[ ,2,t])

    p0_F.R <- rasterFromXYZ(cbind(myDetectors$main.detector.sp$main.cell.x,
                                  myDetectors$main.detector.sp$main.cell.y,
                                  myDetectors$main.detector.sp$p0_F))

    plot(p0_F.R,
         main =  paste0("Females ", years[t]),
         legend.args = list(text = 'p0',
                            side = 4, font = 2, line = 2.5, cex = 0.8))

    ##-- MALE
    myDetectors$main.detector.sp$p0_M <-
      ilogit(logit(results_M$mean$p0[nimConstants$county,t]) +
               results_M$mean$betaDet[1] * nimDataM$detCovs[ ,1,t] +
               results_M$mean$betaDet[2] * nimDataM$detCovs[ ,2,t])

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
  pdf(file = file.path(working_dir, "figures", paste0("p0_beta_classic.pdf")),
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
    plotQuantiles( results_F$sims.list$betaDet[ ,b],
                   at = b - diffSex,
                   col = colSex[1])
    plotQuantiles( results_M$sims.list$betaDet[ ,b],
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

  #}#m




  ## ------     1.6. DETECTABILITY ------
  ## ------       1.6.1. SET-UP ------
  ##-- Create 5km raster for extraction
  regions.r <- habitatRasterResolution$`5km`[["Regions"]]
  regions.r <- crop(regions.r, habitat$habitat.r)

  ##-- Set all habitat cells outside Norway to NA
  regions.r[regions.r[] < 23] <- NA
  plot(regions.r)

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
  n.iter <- 100 # dim(results_M$sims.list$p0)[1]
  iter <- sample(1:dim(resultsSXYZ_MF$sims.list$sigma)[1], size = n.iter)



  ## ------       1.6.2. MALE -----
  DetectabilityRegionsM <- list()
  for(t in 1:n.years){
    p0_M <- matrix(NA,n.iter,nimConstants$n.detectors)
    for(j in 1:nimConstants$n.detectors){
      p0_M[ ,j] <-
        ilogit(logit(results_M$sims.list$p0[iter,nimConstants$county[j],t]) +
                 results_M$sims.list$betaDet[iter,1] * nimDataM$detCovs[j,1,t] +
                 results_M$sims.list$betaDet[iter,2] * nimDataM$detCovs[j,2,t])
    }#j

    system.time(
      DetectabilityRegionsM[[t]] <- GetDetectability_normal(
        p0 = p0_M,
        sigma = results_M$sims.list$sigma[iter],
        habitatxy = regions.xy,
        detectorxy = detectors.xy,
        localDist = 50000,
        regionID = regions.rgmx)
    )
    print(DetectabilityRegionsM[[t]]$summary)
  }#t

  par(mfrow = c(2,6))
  for(t in 1:n.years){
    detectab.r <- regions.r
    detectab.r[isHabitat] <- DetectabilityRegionsM[[t]]$MeanCell
    plot(detectab.r)
  }




  ## ------       1.6.3. FEMALE -----
  DetectabilityRegionsF <- list()
  for(t in 1:n.years){
    p0_F <- matrix(NA,n.iter,nimConstants$n.detectors)
    for(j in 1:nimConstants$n.detectors){
      p0_F[ ,j] <-
        ilogit(logit(results_F$sims.list$p0[iter,nimConstants$county[j],t]) +
                 results_F$sims.list$betaDet[iter,1] * nimDataF$detCovs[j,1,t] +
                 results_F$sims.list$betaDet[iter,2] * nimDataF$detCovs[j,2,t])
    }#j

    system.time(
      DetectabilityRegionsF[[t]] <- GetDetectability_normal(
        p0 = p0_F,
        sigma = results_F$sims.list$sigma[iter],
        habitatxy = regions.xy,
        detectorxy = detectors.xy,
        localDist = 200000,
        regionID = regions.rgmx)
    )
    print(DetectabilityRegionsF[[t]]$summary)
  }#t



  ## ------       1.6.4. SAVE DETECTABILITY OBJECTS ------
  save(DetectabilityRegionsF,
       DetectabilityRegionsM,
       file = file.path( working_dir, "figures",
                         paste0("Detectability5km_classic.RData")))

  # load(file.path( working_dir, "figures",
  #                 paste0("Density5km_classic.RData")))




  ## ------       1.6.3. PLOT DETECTABILITY MAPS -----
  pdf(file = file.path(working_dir, "figures", paste0("Detectability_maps.pdf")),
      width = 10, height = 6)

  ##-- Set color scale
  max <- max(c(unlist(lapply( DetectabilityRegionsM[[t]]$MeanCell,
                              function(x) max(x[], na.rm = T))),
               unlist(lapply( DetectabilityRegionsF[[t]]$MeanCell,
                              function(x) max(x[], na.rm = T)))))
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  col <- rev(terrain.colors(100))

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)
  par(mar = c(0,0,0,0))

  ##-- FEMALES
  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
    detectab.r <- regions.r
    detectab.r[isHabitat] <- DetectabilityRegionsF[[t]]$MeanCell
    image(detectab.r, add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = grey(0.4), col = NA, add=TRUE)
    mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)

    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
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
    image(detectab.r, add=TRUE, breaks=c(cuts, max(cuts)+1000), col = col, legend=FALSE,)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = grey(0.4), col = NA, add=TRUE)
    mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)

    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
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
  pdf(file = file.path(working_dir, "figures", paste0("SexRatio_bars_classic.pdf")),
      width = 12, height = 8.5)
  par(mar = c(5,5,1,1))
  plot(-1000,
       xlim = c(0.5, n.years + 0.5), ylim = c(0.2,0.8),
       xlab = "", ylab = paste("% Female bears"), xaxt = "n", axes = F, cex.lab = 1.6)
  axis(1, at = c(1:n.years), labels = years, cex.axis = 1.6)
  axis(2, at = seq(0.2,0.8,0.1), labels = seq(0.2,0.8,0.1), cex.axis = 1.6)
  #abline(v = (0:n.years)+0.5, lty = 2)
  abline(h = 0.5, lty = 2, col = "gray60")
  for(t in 1:n.years){
    plotQuantiles(x = colSums(DensityCountriesRegionsF[[t]]$PosteriorAllRegions)/
                    colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions),
                  at = t - diffSex,
                  width = 0.18,
                  col = colSex[3])
  }#t
  #box()
  dev.off()




  ## ------       1.7.2. SEX-RATIO RASTER MAPS ------
  ##-- Calculate sex-ratio per 100km2
  SexRatio <- list()
  for(t in 1:length(years)){
    SexRatio[[t]] <- DensityCountriesRegionsF[[t]]$MeanCell/DensityCountriesRegions[[t]]$MeanCell
  }

  ##-- Crop density maps to Norway
  rrCombined <- rrRegions + rrNorway
  SexRatioMap <- list()
  for(t in 1:length(years)){
    SexRatioMap[[t]] <- densityInputRegions$regions.r
    SexRatioMap[[t]][] <- NA
    SexRatioMap[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
    SexRatioMap[[t]][is.na(rrCombined[])] <- NA

    crs(SexRatioMap[[t]]) <- st_crs(habitat$habitat.poly)
  }#t


  pdf(file = file.path(working_dir, "figures", paste0("SexRatio_Maps_classic.pdf")),
      width = 12, height = 8)

  ##-- Set color scale
  max <- 1
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  colfunc <- colorRampPalette(c( colSex[2],"white", colSex[1]))
  col <- colfunc(100)


  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)
  par(mar = c(0,0,0,0))

  ##-- Plot AC maps
  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")
    image(SexRatioMap[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
    mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)

    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
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








  ## ------       1.7.3. SEX-RATIO RASTER MAPS (SMOOTHED) ------
  ##-- Set smoothing factor
  smoothingFactor <- c(3,7,11,21)
  for(SF in smoothingFactor){

    SexRatioMap <- list()
    for(t in 1:length(years)){
      SexRatioMap[[t]] <- densityInputRegions$regions.r
      SexRatioMap[[t]][] <- NA
      SexRatioMap[[t]][!is.na(densityInputRegions$regions.r[])] <- SexRatio[[t]]
      SexRatioMap[[t]] <- terra::focal(terra::rast(SexRatioMap[[t]]), SF, "mean", na.rm=TRUE)
      SexRatioMap[[t]][is.na(rrCombined[])] <- NA
      # proj4string(SexRatioMap[[t]]) <- CRS(proj4string(habitat$habitat.poly))
    }#t



    pdf(file = file.path(working_dir, "figures", paste0("SexRatio_Maps_smooth_",SF,".pdf")),
        width = 12, height = 8)

    ##-- Set color scale
    max <- 1
    cuts <- seq(0, max, length.out = 100) ##-- set breaks
    # colfunc <- colorRampPalette(c( colSex[2],"white", colSex[1]))
    # col <- colfunc(100)

    col <- inferno(100)

    ##-- layout
    mx <- rbind(c(1,rep(1:6, each = 2)),
                c(rep(1:6, each = 2), 6))
    mx <- rbind(mx, mx + 6)
    nf <- layout(mx,
                 widths = c(rep(1,ncol(mx))),
                 heights = rep(1,2))
    #layout.show(nf)
    par(mar = c(0,0,0,0))

    ##-- Plot AC maps
    for(t in 1:length(years)){
      plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
      image( SexRatioMap[[t]], add = TRUE, breaks = c(cuts, max(cuts)+1000), col = col, legend = FALSE)
      plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = grey(0.4), col = NA, add = TRUE)
      mtext(text = years[t], side = 1, -20, adj=0.2, cex=1.2)

      if(t == n.years){
        segments(x0 = 830000, x1 = 830000,
                 y0 = 6730000, y1 = 6730000 + 500000,
                 col = grey(0.3), lwd = 4, lend = 2)
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
  }#SF



  ## ------       1.7.4. SEX-RATIO CONCENTRIC MAPS ------

  ##-- CREATE ARBITRARY DENSITY HOTSPOT LOCATIONS (five regions in Norway)
  threshold <- 0.27
  whichOverThreshold <- which(UDCropped[[1]][] > threshold)
  hotspots <- coordinates(UDCropped[[1]])[whichOverThreshold, ] %>%
    as.data.frame(.) %>%
    st_as_sf(., coords = c("x","y"))
  st_crs(hotspots) <- st_crs(COUNTRIESsimpFig)
  plot( RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])),
        border = NA,
        col = "gray80")
  plot(hotspots, add = T, col = "red", pch = 3)

  ##-- Set sex-ratio color scale
  cuts <- seq(0, 1, length.out = 100)
  colfunc <- colorRampPalette(c("blue4","cornflowerblue","gold","firebrick3","firebrick4"))
  col <- colfunc(100)



  ## ------         1.7.4.1. CUMULATIVE SEX-RATIO MAPS ------

  ##-- CREATE CONCENTRIC AREAS FROM THESE HOTSPOTS
  bandWidth <- 10000
  n.bands <- 50
  bands <- list()
  bands[[1]] <- st_buffer(hotspots, dist = bandWidth) %>%
    st_intersection(.,COUNTRIESsimpFig[1, ]) %>%
    st_union() %>%
    st_as_sf()
  for(b in 2:n.bands){
    bands[[b]] <- st_buffer(bands[[b-1]], dist = bandWidth) %>%
      st_intersection(.,COUNTRIESsimpFig[1, ]) %>%
      st_union() %>%
      st_as_sf()
  }#b


  pdf(file = file.path(working_dir, "figures", paste0("SexRatio_Bands_1.pdf")),
      width = 12, height = 8)

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)

  par(mar = c(0,0,0,0))
  for(t in 1:length(years)){

    N_M.r <- densityInputRegions$regions.r
    N_M.r[] <- NA
    N_M.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsM[[t]]$MeanCell

    N_F.r <- densityInputRegions$regions.r
    N_F.r[] <- NA
    N_F.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsF[[t]]$MeanCell

    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")

    for(b in n.bands:1){
      thisBand <- bands[[b]]

      thisN_F <- mask(N_F.r, as_Spatial(thisBand))
      thisN_F <- sum(thisN_F[], na.rm = T)

      thisN_M <- mask(N_M.r, as_Spatial(thisBand))
      thisN_M <- sum(thisN_M[], na.rm = T)

      this_N <- thisN_F + thisN_M

      thisSR <- round(thisN_F/this_N,2)*100

      plot( bands[[b]], border = NA, add = T, col = col[thisSR])
      print(b)
    }#b
    mtext(text = years[t], side = 1, -20, adj = 0.2, cex = 1.2)
  }#t

  segments(x0 = 830000, x1 = 830000,
           y0 = 6730000, y1 = 6730000 + 500000,
           col = grey(0.3), lwd = 4, lend = 2)
  text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  plot( SexRatioMap[[t]],
        legend.only = T,
        breaks = cuts,
        col = col,
        legend.width = 2,
        axis.args = list(at = c(0.0,0.25,0.50,0.75,1.0),
                         labels = c(0.0,0.25,0.50,0.75,1.0),
                         cex.axis = 1.2),
        smallplot = c(0.85, 0.88, 0.2, 0.6),
        legend.args = list(text = "% Male bears",
                           side = 2, font = 1, line = 1, cex = 1))
  dev.off()



  ## ------         1.7.4.2. BAND-SPECIFIC SEX-RATIO MAPS ------

  n.bands <- 15
  bandWidth <- 25000

  bands <- list()
  bands[[1]] <- st_buffer(hotspots, dist = 10000) %>%
    st_union()
  for(b in 2:n.bands){
    bandwidth <- alpha + beta * b
    bands[[b]] <- st_buffer(bands[[b-1]], dist = bandWidth) %>%
      st_union()
  }#b

  for(b in n.bands:2){
    bands[[b]] <- bands[[b]] %>%
      st_difference(., bands[[b-1]]) %>%
      st_intersection(.,COUNTRIESsimpFig[1, ]) %>%
      st_union() %>%
      st_as_sf()
  }#b
  bands[[1]] <- bands[[1]] %>%
    st_intersection(.,COUNTRIESsimpFig[1, ]) %>%
    st_union() %>%
    st_as_sf()


  ##-- Calculate sex-ratio per band
  pdf(file = file.path(working_dir, "figures", paste0("SexRatio_Bands_2.pdf")),
      width = 12, height = 8)

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)
  par(mar = c(0,0,0,0))
  for(t in 1:length(years)){

    N_M.r <- densityInputRegions$regions.r
    N_M.r[] <- NA
    N_M.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsM[[t]]$MeanCell

    N_F.r <- densityInputRegions$regions.r
    N_F.r[] <- NA
    N_F.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsF[[t]]$MeanCell

    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1,])), border = NA, col = "gray80")

    for(b in 1:n.bands){
      thisBand <- bands[[b]]

      thisN_F <- mask(N_F.r, as_Spatial(thisBand))
      thisN_F <- sum(thisN_F[], na.rm = T)

      thisN_M <- mask(N_M.r, as_Spatial(thisBand))
      thisN_M <- sum(thisN_M[], na.rm = T)

      this_N <- thisN_F + thisN_M
      print(this_N)

      thisSR <- round(thisN_F/this_N,2)*100

      thisAlpha <- round(this_N/20,1)
      thisAlpha <- ifelse(thisAlpha > 1, 1, thisAlpha)

      #print(thisAlpha)

      if(this_N > 7) {
        plot( bands[[b]], border = NA, add = T,
              col = adjustcolor(col[thisSR], alpha.f = 1))
      }
    }#b
    mtext(text = years[t], side = 1, -20, adj = 0.2, cex = 1.2)
  }#t

  segments(x0 = 830000, x1 = 830000,
           y0 = 6730000, y1 = 6730000 + 500000,
           col = grey(0.3), lwd = 4, lend = 2)
  text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
  plot( SexRatioMap[[t]],
        legend.only = T,
        breaks = cuts,
        col = col,
        legend.width = 2,
        axis.args = list(at = c(0.0,0.25,0.50,0.75,1.0),
                         labels = c(0.0,0.25,0.50,0.75,1.0),
                         cex.axis = 1.2),
        smallplot = c(0.85, 0.88, 0.2, 0.6),
        legend.args = list(text = "% Female bears",
                           side = 2, font = 1, line = 1, cex = 1))
  dev.off()





  ## ------         1.7.4.3. SEX-RATIO Histograms WITH DISTANCE ------
  ##-- Calculate sex-ratio per band
  # pdf(file = file.path(working_dir, "figures", paste0("SexRatio_Histograms.pdf")),
  #     width = 12, height = 8)


  N <- list()
  for(t in 1:length(years)){

    N_M.r <- densityInputRegions$regions.r
    N_M.r[] <- NA
    N_M.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsM[[t]]$MeanCell

    N_F.r <- densityInputRegions$regions.r
    N_F.r[] <- NA
    N_F.r[!is.na(densityInputRegions$regions.r[])] <- DensityCountriesRegionsF[[t]]$MeanCell

    N_f <- cbind.data.frame(
      "N" = rep(NA, n.bands),
      "distance" = seq( from = bandWidth/2,
                        to = bandWidth*n.bands,
                        by = bandWidth)/1000,
      "sex" = "Female",
      "year" = years[t])

    N_m <- cbind.data.frame(
      "N" = rep(NA, n.bands),
      "distance" = seq( from = bandWidth/2,
                        to = bandWidth*n.bands,
                        by = bandWidth)/1000,
      "sex" = "Male",
      "year" = years[t])

    for(b in 1:n.bands){
      thisBand <- bands[[b]]

      thisN_F <- mask(N_F.r, as_Spatial(thisBand))
      N_f$N[b] <- sum(thisN_F[], na.rm = T)

      thisN_M <- mask(N_M.r, as_Spatial(thisBand))
      N_m$N[b] <- sum(thisN_M[], na.rm = T)
    }#b

    N[[t]] <- rbind(N_m, N_f)
  }#t

  do.call(rbind, N) %>%
    filter(., distance < 200) %>%
    ggplot(., aes(x = distance, y = N, fill= sex)) +
    geom_bar(stat = "identity", position = 'dodge') +
    labs(title = "Evolution of sex-ratio with distance",
         subtitle = "Norwegian bear population",
         y = "number of bears", x = "Distance from core areas (km)") +
    facet_wrap(~year, nrow = 2)


  #dev.off()




  ## ------         1.7.4.4. MAPS OF DIFFERENCES IN UD-DENSITIES ------

  ACDiff <- list()
  for(t in 1:length(years)){
    ACDiff[[t]] <- densityInputRegions$regions.r
    ACDiff[[t]][] <- NA
    ACDiff[[t]][!is.na(densityInputRegions$regions.r[])] <- AC_M100km2[[t]]- AC_F100km2[[t]]
    ACDiff[[t]][is.na(rrCombined[])] <- NA
  }#t


  ##-- Set color scale
  min <- min(unlist(lapply(ACDiff, function(x) min(x[], na.rm = T))))
  max <- max(unlist(lapply(ACDiff, function(x) max(x[], na.rm = T))))
  cuts <- seq(min, max, length.out = 100) ##-- set breaks

  pdf(file = file.path(working_dir, "figures", paste0("AC_DifferenceMaps.pdf")),
      width = 12, height = 8)

  ##-- layout
  mx <- rbind(c(1,rep(1:6, each = 2)),
              c(rep(1:6, each = 2), 6))
  mx <- rbind(mx, mx + 6)
  nf <- layout(mx,
               widths = c(rep(1,ncol(mx))),
               heights = rep(1,2))
  #layout.show(nf)
  par(mar = c(0,0,0,0))

  for(t in 1:length(years)){
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = NA, col = "gray80")
    image(ACDiff[[t]], add=TRUE, breaks=c(cuts, max(cuts)+0.1), col = col, legend=FALSE,)
    plot(RemoveHolesSp(as_Spatial(COUNTRIESsimpFig[1, ])), border = grey(0.4), col = NA, add=TRUE)
    mtext(text = years[t], side = 1, -17, adj = 0.2, cex = 1.2, font = 2)

    if(t == n.years){
      segments(x0 = 830000, x1 = 830000,
               y0 = 6730000, y1 = 6730000 + 500000,
               col = grey(0.3), lwd = 4, lend = 2)
      text(750000, 6730000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
      plot( UDCropped_M[[t]],
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

  ## ------   2. TABLES -----
  gc()
  ## ------     2.1. ABUNDANCE -----
  ## ------       2.1.1. ALL YEARS, BOTH SEX -----
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
            file = file.path(WDTables,
                             paste0(myVars$modelName, "_classic","_NAllYears.csv")))
  write.csv(NCarRegionMean,
            file = file.path(WDTables,
                             paste0(myVars$modelName, "_classic","_NAllYears_mean.csv")))

  ##-- print .tex
  row.names(NCarRegionEstimates) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  print(xtable( NCarRegionEstimates,
                type = "latex",
                align = paste(c("l",rep("c",ncol(NCarRegionEstimates))),collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list(list(seq(1, nrow(NCarRegionEstimates), by = 2)), "\\rowcolor[gray]{.96} "),
        file = file.path(WDTables,
                         paste0(myVars$modelName, "_classic","_NAllYears.tex")))




  ## ------       2.1.2. LAST YEAR N PER SEX PER COUNTY -----
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
             file = file.path(WDTables, paste0(myVars$modelName, "_classic","_NLastYearPerSex.csv")))

  ##-- print .tex
  row.names(NCountyEstimatesLastRegions) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  print(xtable(NCountyEstimatesLastRegions, type = "latex",
               align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path( WDTables,
                          paste0(myVars$modelName, "_classic","_NLastYearPerSex.tex")))



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
             file = file.path(WDTables, paste0(myVars$modelName, "_classic","_NLastYearPerSex_UD.csv")))

  ##-- print .tex
  row.names(NCountyEstimatesLastRegions_UD) <- c(paste0("\\hspace{0.1cm} ",idcountyNOR),"TOTAL")
  print(xtable(NCountyEstimatesLastRegions_UD, type = "latex",
               align = paste(c("l",rep("c",ncol(NCountyEstimatesLastRegions_UD))), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions_UD),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path( WDTables,
                          paste0(myVars$modelName, "_classic","_NLastYearPerSex_UD.tex")))




  ## ------       2.1.3. MAKE A TABLE 2 last years -----
  NCountyEstimatesLast2Regions <- matrix("", ncol = 6, nrow = length(idcountyTable))
  row.names(NCountyEstimatesLast2Regions) <- c(idcountyTable)
  colnames(NCountyEstimatesLast2Regions) <- c(paste("Females", years[n.years-1]),
                                              paste("Males", years[n.years-1]),
                                              paste("Total", years[n.years-1]),
                                              paste("Females", years[n.years]),
                                              paste("Males", years[n.years]),
                                              paste("Total", years[n.years]))

  ##-- Fill-in table
  for(t in (n.years-1):n.years){
    for(i in 1:length(idcountyTable)){
      ##-- FEMALES
      NCountyEstimatesLast2Regions[idcountyTable[i],paste("Females",years[t])] <-
        paste(round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegionsF[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")

      ##-- MALES
      NCountyEstimatesLast2Regions[idcountyTable[i],paste("Males",years[t])] <-
        paste(round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegionsM[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")

      ##-- TOTAL
      NCountyEstimatesLast2Regions[idcountyTable[i],paste("Total",years[t])] <-
        paste(round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"mean"],digits = 1)," (",
              round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CILow"],digits = 0),"-",
              round(DensityCountriesRegions[[t]]$summary[idcountyTable[i],"95%CIHigh"],digits = 0),")",sep="")
    }#i
  }#t

  ##-- ADJUST NAMES
  idcounty1 <- idcountyTable
  idcounty1[which(idcounty1 %in% "Total")] <- "TOTAL"
  row.names(NCountyEstimatesLast2Regions) <- idcounty1

  ##-- print .csv
  write.csv(NCountyEstimatesLast2Regions,
            file = file.path(WDTables, paste0(myVars$modelName, "_classic","_NLast2YearsPerSex.csv")))


  ##-- print .tex
  row.names(NCountyEstimatesLast2Regions) <- c(paste0("\\hspace{0.5cm} ", idcountyNOR), "TOTAL")

  NCountyEstimatesLast2Regions <- rbind(c("F","M","Total","F","M","Total"),
                                        NCountyEstimatesLast2Regions)

  addtorow <- list()
  addtorow$pos <- list(c(0),0)
  addtorow$command <- c(paste0(paste0('& \\multicolumn{3}{c}{', years[(n.years-1):n.years],
                                      '}', collapse=''), '\\\\'),rep("\\rowcolor[gray]{.95}",1))

  print(xtable(NCountyEstimatesLast2Regions, type = "latex",
               align = paste(c("l",rep("c",3),"|",rep("c",3)), collapse = "")),
        sanitize.text.function=function(x){x},
        floating = FALSE,
        add.to.row = addtorow,
        include.colnames = F,
        file = file.path(WDTables,
                         paste0(myVars$modelName, "_classic","_NLast2YearsPerSex.tex")))




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
            file = file.path(WDTables, paste0(myVars$modelName, "_classic","_NAllYearsPerSex.csv")))


  ##-- print .tex
  row.names(NCountyEstimatesAllSexRegions) <- c("", paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")

  print(xtable(NCountyEstimatesAllSexRegions, type = "latex",
               align = paste(c("l",rep("c",ncol(NCountyEstimatesAllSexRegions))),collapse = "")),
        sanitize.text.function = function(x){x},
        floating = FALSE,
        add.to.row = list(list(seq(1,nrow(NCountyEstimatesLastRegions),by=2)),"\\rowcolor[gray]{.95} "),
        file = file.path(WDTables, paste0(myVars$modelName, "_classic","_NAllYearsPerSex.tex")))




  ## ------     2.2. NUMBER OF NGS SAMPLES, IDs & DEAD RECOVERIES ------
  ##-- Load filtered datasets
  load(file.path(myVars$WD, "MODELS", modelNameM, "DATA", paste0(modelNameM,"_NGSData.RData")))
  myData.aliveM <- myData.alive$myData.sp
  myData.deadM <- myData.dead

  load(file.path(myVars$WD, "MODELS", modelNameF, "DATA", paste0(modelNameF,"_NGSData.RData")))
  myData.aliveF <- myData.alive$myData.sp
  myData.deadF <- myData.dead

  ##-- SOME TALLIES TO CHECK THINGS
  ##-- NGS
  NGS <- rbind(myData.aliveM, myData.aliveF)
  table(NGS$Year)
  length(NGS)

  ##-- FOR REPORT SUMMARY
  length(NGS$Id)
  length(NGS$Id[NGS$Sex=="Hunn"])
  length(NGS$Id[NGS$Sex=="Hann"])
  length(NGS$Id[NGS$Country=="(S)"])/nrow(NGS)
  length(unique(NGS$Id))
  length(unique(NGS$Id[NGS$Sex=="Hunn"]))
  length(unique(NGS$Id[NGS$Sex=="Hann"]))

  ##-- DEAD RECOVERY
  dead <- rbind(myData.deadF, myData.deadM)
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
      NGS_SEX["number of NGS samples", ye[t] + sex1[s]] <-  nrow(temp)
      NGS_SEX["number of NGS individuals", ye[t] + sex1[s]] <- length(unique(temp$Id))
    }#t
  }#s

  ##-- print .csv
  write.csv(NGS_SEX,
            file = file.path(WDTables, paste0("NGS_SEX_classic.csv")))

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
        file = file.path(WDTables, paste0("NGS_SEX_classic.tex")))





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
  MortalityNames <- unique(as.character(myFullData.sp$dead.recovery$DeathCause))
  legalCauses <- MortalityNames[grep("Lisensfelling", MortalityNames)]
  legalCauses <- c(legalCauses, MortalityNames[grep("tamdyr", MortalityNames)])
  legalCauses <- c(legalCauses, MortalityNames[grep("SNO", MortalityNames)])
  legalCauses <- c(legalCauses, MortalityNames[grep("Skadefelling", MortalityNames)])
  legalCauses <- c(legalCauses, MortalityNames[grep("Politibeslutning", MortalityNames)])
  legalCauses <- c(legalCauses, MortalityNames[grep("menneske", MortalityNames)])
  MortalityNames[!MortalityNames %in% legalCauses]

  ##-- SEPARATE MORTALITIES
  cause <- c("other","legal culling")
  for(t in 1:n.years){
    for(s in 1:2){
      for(d in 1:2){
        if(d==1){
          temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & !(dead$DeathCause %in% legalCauses), ]
        } else {
          temp <- dead[dead$Year == years[t] & dead$Sex==sex[s] & dead$DeathCause %in% legalCauses, ]
        }
        row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[,1]=="Norway" )
        Dead_SEX[row,ye[t] + sex1[s]+1 ] <- length(unique(temp$Id[temp$Country %in% "(N)" ]))

        row <- which(rownames(Dead_SEX)==cause[d] & Dead_SEX[,1]=="Sweden" )
        Dead_SEX[row,ye[t] + sex1[s]+1] <- length(unique(temp$Id[temp$Country %in% "(S)"]))
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
        file = file.path(WDTables, paste0("DeadidCountrySEX_classic.tex")))




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
    propDetected["F",t] <- getCleanEstimates(
      n.detected_F[t]/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions))

    propDetected["M",t] <- getCleanEstimates(
      n.detected_M[t]/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions))

    propDetected["Total",t] <- getCleanEstimates(
      (n.detected_F[t]+n.detected_M[t])/
        (colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions)+
           colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions)))
  }#t

  ##-- print .csv
  write.csv(propDetected,
            file = file.path(WDTables, "PropDetectedIds_classic.csv"))

  ##-- print .tex
  print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
        floating = FALSE, sanitize.text.function=function(x){x},
        add.to.row = list(list(seq(1,nrow(propDetected), by = 2)),"\\rowcolor[gray]{.96} "),
        file = file.path(WDTables, "PropDetectedIds_classic.tex"))




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

      country <- countryRaster[cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
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

    print(t)
  }#t
  ##-- Print number of individuals detected per country
  print(N_det_by_country)


  ##-- Calculate number of individuals detected/undetected in Norway
  N_NOR <- matrix(NA,3,n.years)
  dimnames(N_NOR) <- list("#individuals" = c("Detected","Undetected","Total"),
                          "Years" = c(years))
  for(t in 1:n.years){
    N_det <- N_undet <- N_tot <- rep(NA,n.iter)
    for(iter in 1:n.iter){
      country <- countryRaster[cellFromXY(norRaster,resultsSXYZ_MF$sims.list$sxy[iter, ,1:2,t])]
      isNor <- country %in% 2

      ##-- Detected individuals
      N_det[iter] <- sum(isDetected[ ,t] & isAlive[iter, ,t] & isNor)

      ##-- Undetected individuals
      N_undet[iter] <- sum(!isDetected[ ,t] & isAlive[iter, ,t] & isNor)

      ##-- Total Norway
      N_tot[iter] <- N_det[iter] + N_undet[iter]
    }#iter

    N_NOR["Detected",t] <- getCleanEstimates(N_det)
    N_NOR["Undetected",t] <- getCleanEstimates(N_det)
    N_NOR["Total",t] <- getCleanEstimates(N_tot)
    print(t)
  }#t
  ##-- Print number of individuals with AC in Norway
  print(N_NOR)





  ## ------     2.3. VITAL RATES ------
  parameters <- c("rho","phi", "h", "w")
  sex <- c("F", "M")
  vitalRate <- matrix(NA, nrow = length(parameters)+1, ncol = (n.years)*2-2)
  rownames(vitalRate) <- c("", unlist(lapply(as.list(parameters), function(x)rep(x,1))))
  colnames(vitalRate) <- c(unlist(lapply(years[1:(length(years)-1)], function(x)rep(paste(x,x+1,sep=" to "),2))))
  vitalRate[1, ] <- c(rep(sex,(n.years-1)) )

  for(s in 1:2){
    if(s == 1){results <- results_F} else {results <- results_M}

    col <- which(vitalRate[1, ] == sex[s])

    ##-- Per capita recruitment
    if(any(grep("rho",names(results$sims.list)))){
      vitalRate["rho",col] <- getCleanEstimates(results$sims.list$rho, moment = "median")
    } else {
      if(s == 1){
        for(t in 1:(n.years-1)){
          n.recruits <- rowSums(isAvail[ ,isFemale,t] * isAlive[ ,isFemale,t+1])
          alivetminus1 <- rowSums(isAlive[ ,isFemale,t])
          vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
        }#t
      } else {
        for(t in 1:(n.years-1)){
          n.recruits <- rowSums(isAvail[ ,isMale,t] * isAlive[ ,isMale,t+1])
          alivetminus1 <- rowSums(isAlive[ ,isMale,t])
          vitalRate["rho",col[t]] <- getCleanEstimates(n.recruits/alivetminus1, moment = "median")
        }#t
      }
    }


    ##-- Mortality & Survival
    if(any(grep("mhH",names(results$sims.list)))){
      ##-- Calculate mortality from estimated hazard rates (mhH and mhW)
      mhH1 <- exp(results$sims.list$mhH)
      mhW1 <- exp(results$sims.list$mhW)
      h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
      w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
      phi <- 1-h-w
      vitalRate["phi",col] <- apply(phi, 2, function(x) getCleanEstimates(x, moment = "median"))
      vitalRate["h",col] <- apply(h, 2, function(x) getCleanEstimates(x, moment = "median"))
      vitalRate["w",col] <- apply(w, 2, function(x) getCleanEstimates(x, moment = "median"))
    } else {
      if(s == 1){
        y.dead <- y.deadF
        z <- resultsSXYZ_MF$sims.list$z[ ,isFemale, ]
      } else {
        y.dead <- y.deadM
        z <- resultsSXYZ_MF$sims.list$z[ ,isMale, ]
      }

      ##-- Extract survival from posteriors
      vitalRate["phi",col] <- apply(results$sims.list$phi, 2, function(x) getCleanEstimates(x, moment = "median"))

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
             file = file.path(WDTables, "VitalRates_classic.csv"))

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
        file = file.path(WDTables, "VitalRates_classic.tex"))




  ## ------     2.4. DERIVED PARAMETERS FROM ABUNDANCE ------
  ## ------       2.4.1. DERIVE SEX-RATIO ------

  ##-- REGION-SPECIFIC PROPORTION OF FEMALES
  PropFemale_regions <- list()
  for(t in 1:n.years){
    PropFemale_regions[[t]] <- DensityCountriesRegionsF[[t]]$PosteriorRegions/
      (DensityCountriesRegionsM[[t]]$PosteriorRegions +
         DensityCountriesRegionsF[[t]]$PosteriorRegions)
    print(rowMeans(PropFemale_regions[[t]],na.rm = T))
  }#t


  ##-- OVERALL PROPORTION OF FEMALES
  dim(DensityCountriesRegionsF[[t]]$PosteriorAllRegions)
  PropFemale <- list()
  for(t in 1:n.years){
    PropFemale[[t]] <- colSums(DensityCountriesRegionsF[[t]]$PosteriorAllRegions)/
      (colSums(DensityCountriesRegionsM[[t]]$PosteriorAllRegions) +
         colSums(DensityCountriesRegionsF[[t]]$PosteriorAllRegions))
  }#t


  ##-- Format table
  propFemale_tab <- matrix(0, ncol = n.years, nrow = 8)
  row.names(propFemale_tab) <- idcountyTable
  colnames(propFemale_tab) <- years
  for(t in 1:n.years){
    for(c in 1:7){
      propFemale_tab[c,t] <- getCleanEstimates(na.omit(PropFemale_regions[[t]][c, ]))
    }#c
    propFemale_tab[8,t] <- getCleanEstimates( PropFemale[[t]] )
  }#t


  ##-- print .tex
  row.names(propFemale_tab) <- c(paste0("\\hspace{0.1cm} ", idcountyNOR), "TOTAL")
  print(xtable( propFemale_tab,
                type = "latex",
                align = paste(c("l",rep("c",ncol(propFemale_tab))),collapse = "")),
        floating = FALSE,
        sanitize.text.function = function(x){x},
        add.to.row = list(list(seq(1, nrow(propFemale_tab), by = 2)), "\\rowcolor[gray]{.96} "),
        file = file.path(WDTables, "propFemale_classic.tex"))





  # ## ------       2.4.2. DERIVE DENSITY ------
  # habbRCarRegionsTRY <- rrRegions
  # habbRCarRegionsTRY[!is.na(habbRCarRegionsTRY[])] <-1
  # habbRCarRegionsTRY[] <- as.numeric(habbRCarRegionsTRY[])
  # areaSqKm <- gArea(rasterToPolygons(habbRCarRegionsTRY,function(x) x > 0, dissolve = T))*1e-6
  #
  # ##-- Multiplied by 100 to get per 100km2
  # DensityCountriesRegions[[t]]$summary["Total","mean"]/areaSqKm*100
  # DensityCountriesRegions[[t]]$summary["Total","95%CILow"]/areaSqKm*100
  # DensityCountriesRegions[[t]]$summary["Total","95%CIHigh"]/areaSqKm*100



  ## ------       2.4.3. PROPORTION OF INDIVIDUALS DETECTED ------
  ##-- Get the number of individuals detected each year
  n.detected_F <-  apply(nimDataF$y.alive[ ,1, ], 2, function(x)sum(x>0))
  n.detected_M <-  apply(nimDataM$y.alive[ ,1, ], 2, function(x)sum(x>0))

  propDetected <- matrix("", ncol = n.years, nrow = 3)
  row.names(propDetected) <- c("F","M","Total")
  colnames(propDetected) <- years
  for(t in 1:n.years){
    propDetected["F",t] <- getCleanEstimates(n.detected_F[t]/colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions))
    propDetected["M",t] <- getCleanEstimates(n.detected_M[t]/colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions))
    propDetected["Total",t] <- getCleanEstimates((n.detected_F[t]+n.detected_M[t])/
                                                   (colSums(DensityCountriesRegionsF[[t]]$PosteriorRegions)+
                                                      colSums(DensityCountriesRegionsM[[t]]$PosteriorRegions)))
  }#t

  ##-- print .csv
  write.csv(propDetected,
            file = file.path(WDTables, "PropDetectedIds_classic.csv"))

  ##-- print .tex
  print(xtable(propDetected, type = "latex", align=paste(c("l",rep("c",ncol(propDetected))),collapse = "")),
        floating = FALSE, sanitize.text.function=function(x){x},
        add.to.row = list(list(seq(1,nrow(propDetected), by = 2)),"\\rowcolor[gray]{.96} "),
        file = file.path(WDTables, "PropDetectedIds_classic.tex"))




  ## ------       2.4.4. GROWTH RATE ------
  growthRate <- list()
  for(t in 1:(n.years-1)){
    growthRate[[t]] <- colSums(DensityCountriesRegions[[t+1]]$PosteriorAllRegions)/
      colSums(DensityCountriesRegions[[t]]$PosteriorAllRegions)
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
        file = file.path(WDTables, "GrowthRates_classic.tex"))




  ## ------     2.5. TABLE OTHERS ------
  if(is.null(dim(results$sims.list$betaDens))){
    parameters <- c("tau","betaDead","sigma","betaDet","betaDet")
    sex <- c("F","M")
    TableOthers <- matrix(NA, nrow = length(parameters), ncol = 3)
    rownames(TableOthers) <- parameters
    colnames(TableOthers) <- c("", sex)
    TableOthers[ ,1] <- c("$\\tau$","$\\beta_{dead}$","$\\sigma$","$\\beta_{roads}$","$\\beta_{obs}$")

    for(s in 1:2){
      if(s == 1){results <- results_F} else {results <- results_M}
      TableOthers["tau",sex[s]] <- getCleanEstimates(results$sims.list$tau/res(habitat$habitat.r)[1], moment = "median")
      TableOthers[which(parameters == "betaDead"),sex[s]] <- getCleanEstimates(results$sims.list$betaDens,moment = "median")
      TableOthers["sigma",sex[s]] <- getCleanEstimates(results$sims.list$sigma/res(habitat$habitat.r)[1],moment = "median")
      TableOthers[which(parameters == "betaDet"),sex[s]] <- apply(results$sims.list$betaDet, 2,function(x) getCleanEstimates(x,moment = "median"))
    }#s

    ##-- Deal with negative values
    TableOthers <- gsub("--", "-(-)", TableOthers)

    ##-- Change row names
    row.names(TableOthers) <- c("Spatial process","","Detection Process","","")


    ##-- Print .tex
    multirow <- c("\\multirow{2}{*}{\\textbf{Spatial process}}",
                  "\\multirow{3}{*}{\\textbf{Detection process}}")
    multirowadd <- matrix(c(multirow[1],"",multirow[2],"",""), ncol = 1)
    TableOthers <- data.frame(cbind(multirowadd, TableOthers))

    addtorow <- list()
    addtorow$pos <- list(0,0,3)
    addtorow$command <- c(paste0("& {\\textbf{Parameters}} ",
                                 paste0('& {\\textbf{',  sex, '}}', collapse = ''), '\\\\'),
                          rep("\\rowcolor[gray]{.95}",1),
                          "\\hline")

    print(xtable(TableOthers, type = "latex",
                 align = paste(c("ll", rep("c",ncol(TableOthers)-1)), collapse = "")),
          sanitize.text.function = function(x){x},
          floating = FALSE,
          include.rownames = FALSE,
          include.colnames = FALSE,
          add.to.row = addtorow,
          file = file.path(WDTables, paste0("TableParametersOthers_classic.tex")))
  } else {
    parameters <- c("tau","betaDens1","betaDens1",
                    "betaDead2","betaDead2",
                    "sigma","betaDet","betaDet")
    sex <- c("F","M")
    TableOthers <- matrix(NA, nrow = length(parameters), ncol = 3)
    rownames(TableOthers) <- parameters
    colnames(TableOthers) <- c("", sex)
    TableOthers[ ,1] <- c("$\\tau$",
                          "$\\beta_{dead}_{1}$","$\\beta_{obs}_{1}$",
                          "$\\beta_{dead}_{2}$","$\\beta_{obs}_{2}$",
                          "$\\sigma$","$\\beta_{roads}$","$\\beta_{skandobs}$")

    for(s in 1:2){
      if(s == 1){results <- results_F} else {results <- results_M}
      TableOthers["tau",sex[s]] <- getCleanEstimates(results$sims.list$tau/1000, moment = "median")
      TableOthers[which(parameters == "betaDead1"),sex[s]] <- apply(results$sims.list$betaDens[,,1], 2,function(x) getCleanEstimates(x,moment = "median"))
      TableOthers[which(parameters == "betaDead2"),sex[s]] <- apply(results$sims.list$betaDens[,,2], 2,function(x) getCleanEstimates(x,moment = "median"))
      TableOthers["sigma",sex[s]] <- getCleanEstimates(results$sims.list$sigma/1000, moment = "median")
      TableOthers[which(parameters == "betaDet"),sex[s]] <- apply(results$sims.list$betaDet, 2,function(x) getCleanEstimates(x,moment = "median"))
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
          file = file.path(WDTables, "TableParametersOthers_classic.tex"))
  }

}

