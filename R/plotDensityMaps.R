#' Plot population density maps
#'
#' Function to plot (and print as \code{.pdf} files) population density maps from SCR or OPSCR models.
#'
#' The \code{plotDensityMaps} function uses the outputs from the \code{getDensityInput} and \code{getDensity} functions as inputs.
#'
#' @param input Raster used for extraction, as prepared by the \code{getDensityInput} function.
#' @param estimates Density estimates as provided by the \code{getDensity} or \code{getSpaceUse} functions.
#' @param unit Unit (in km-2) to be used for plotting density.
#' @param mask A \code{raster} object denoting the region for which density should be plotted.
#' All grid cells with values different from NA will be plotted.
#' Note that the extent and resolution of this raster must match the extent and resolution of the raster in \code{input}.
#' @param background A \code{sf} object of the overall region for which density should be plotted.
#' @param type A \code{character} string denoting the type of plot to be printed. Can be one of \code{"last.year"}, \code{"time.series"} or \code{"all"}.
#' @param path A \code{path} denoting where to save the density maps.
#' @param name A \code{character} string with the name to be used when saving .png and .tif files.
#' 
#' @return This function prints out \code{.png} image of density maps in the specified folder.
#' 
#' @author Pierre Dupont
#'
#' @importFrom raster res image extent writeRaster
#' @importFrom graphics plot layout par segments mtext text rasterImage
#' @importFrom sf st_geometry 
#' @importFrom grDevices png colorRampPalette
#'
#' @rdname plotDensityMaps
#' @export
plotDensityMaps <- function( 
    input,
    estimates,
    unit = 100,
    mask,
    background = NULL,
    type = c("time.series", "last.year", "summary"),
    path = getwd(),
    species = NULL,
    labels = NULL,
    x.labels = NULL,
    y.labels = NULL,
    caption = NULL,
    export.raster = TRUE,
    name = "UD_Density")
{
  
  ##-- Initial checks
  if(is.null(mask)){ mask <- input }
  if(!inherits(estimates,"list")){ estimates <- list(estimates) }
  if(is.null(background)){background <- COUNTRIES}
  if(sum(grep("bear", species, ignore.case = T)) > 0|
     sum(grep("bjørn", species, ignore.case = T)) > 0|
     sum(grep("bjorn", species, ignore.case = T)) > 0) { 
    species <- "bear" 
    engSpecies <- "brown bear"
    norSpecies <- "brunbjørn"
    }
  if(sum(grep("wolf", species, ignore.case = T)) > 0|
     sum(grep("wolves", species, ignore.case = T)) > 0|
     sum(grep("ulv", species, ignore.case = T)) > 0) {
    species <- "wolf" 
    engSpecies <- "wolf" 
    norSpecies <- "ulven"
    }
  if(sum(grep("wolverine", species, ignore.case = T))>0|
     sum(grep("jerv", species, ignore.case = T))>0|
     sum(grep("järv", species, ignore.case = T))>0) {
    species <- "wolverine" 
    engSpecies <- "wolverine" 
    norSpecies <- "jerv"
    }
  
  
  ##-- Convert densities to the desired density unit (usually inds.100km-2)
  conversionFactor <- unit/( raster::res(input)[1]/1000)^2
  
  
  ##-- Rasterize and mask density maps 
  density <- list()
  for(t in 1:length(estimates)){
    density[[t]] <- input
    density[[t]][] <- NA
    density[[t]][!is.na(input[])] <- estimates[[t]]$MeanCell * conversionFactor
    density[[t]][is.na(mask[])] <- NA
  }#t
  names(density) <- names(estimates)
  
  
  ##-- Set color scale
  max <- max(unlist(lapply(density, function(x) max(x[], na.rm = T))))
  cuts <- seq(0, max, length.out = 100) ##-- set breaks
  colfunc <- grDevices::colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
  col <- colfunc(100)
  
  
  ##-- Set x- and y-limits
  xLims <- raster::extent(background)[1:2]
  xRange <- diff(xLims)
  yLims <- raster::extent(background)[3:4]
  yRange <- diff(yLims)
  
  
  ##-- Density maps time series
  if("time.series" %in% type){
    
    ##-- layout
    L <- length(density)
    if(L < 6){ nrows <- 1 } else {
      if(L < 13){ nrows <- 2 } else {
        if(L < 22){ nrows <- 3 } else {
          if(L < 33){ nrows <- 4 } else {
            nrows <- 5
          }}}}
    ncols <- ceiling(L/nrows)
    
    grDevices::png(filename = file.path(path, "figures", paste0(name,"_TimeSeries.png")),
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
    
    ##-- legend coordinates
    legend.x <- xLims[1] + 0.7 * xRange
    legend.y <- yLims[1] + 0.4 * yRange
    
    ##-- Plot density maps
    graphics::par(mar = c(0,0,0,0))
    
    for(t in 1:length(density)){
      ##-- Plot density
      plot(sf::st_geometry(background), border = NA, col = "gray80")
      raster::image( density[[t]], add = TRUE,
                     breaks = c(cuts, max(cuts) + 1000),
                     col = col, legend = FALSE)
      plot( sf::st_geometry(background),
            border = "gray40", col = NA, add = TRUE)
      
      ##-- Add year if available
      if(!is.null(names(estimates))){
        graphics::mtext(text = names(estimates)[t],
                        side = 1, line =  -20,
                        adj = 0.2, cex = 1.2)
      }
      
      ##-- Add legend
      if(t == length(density)){
        graphics::segments(
          x0 = legend.x, x1 = legend.x,
          y0 = legend.y-250000, y1 = legend.y + 250000,
          col = "grey30", lwd = 4, lend = 2)
        graphics::text(
          x = legend.x - 0.05 * xRange,
          y = legend.y,
          labels = "500 km", srt = 90, cex = 1.4)
        raster::plot( density[[t]],
                      legend.only = T, breaks = cuts,
                      col = col, legend.width = 2,
                      axis.args = list( at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                        labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                        cex.axis = 1.2),
                      smallplot = c(0.8, 0.85, 0.2, 0.6), 
                      legend.args = list(text = paste0("Individuals/", unit, " km2"),
                                         side = 2, font = 1, line = 0, cex = 1))
        ######----- NEED TO FIX LEGEND TEXT 
        ######----- expression("Individuals/100 km"^ 2)
      }#if
      
      ##-- Export rasters
      if(export.raster){
        writeRaster( density[[t]],
                     file.path( path, "rasters",
                                paste0( species, "_",
                                        raster::res(input)[1]/1000, "km", 
                                        names(estimates)[t], ".tif")),
                     overwrite = TRUE)
      }
      
    }#t
    dev.off()
  }
  
  
  ##-- Last year's density map
  if("last.year" %in% type){
    
    grDevices::png(filename = file.path(path, paste0(name,"_LastYear.png")),
                   width = 8, height = 8, units = "in", pointsize = 12,
                   res = 300, bg = NA)
    
    graphics::par(mar = c(0,0,0,0))
    plot(sf::st_geometry(background), border = NA, col = "gray80")
    raster::image( density[[length(density)]], add = TRUE,
                   breaks = c(cuts, max(cuts)+1000),
                   col = col, legend = FALSE)
    plot( sf::st_geometry(background),
          border = "gray40", col = NA, add = TRUE)
    
    ##-- Add caption if available
    if(!is.null(names(estimates))){
      mtext(text = names(estimates)[t],
            side = 1, line =  -25,
            adj = 0.25, cex = 3, font = 2)
    }
    
    ##-- Add legend
    legend.x <- xLims[1] + 0.9 * xRange
    legend.y <- yLims[1] + 0.3 * yRange
    
    graphics::segments(
      x0 = legend.x, x1 = legend.x,
      y0 = legend.y-250000, y1 = legend.y + 250000,
      col = "gray30", lwd = 4, lend = 2)
    graphics::text(
      x = legend.x - 0.05 * xRange,
      y = legend.y,
      labels = "500 km", srt = 90, cex = 1.4)
    
    raster::plot( density[[t]],
                  legend.only = T, breaks = cuts,
                  col = col, legend.width = 2,
                  axis.args = list( at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                    labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                    cex.axis = 1.2),
                  smallplot = c(0.85, 0.88, 0.2, 0.4),
                  legend.args = list(text = paste0("Individuals/", unit, " km2"),
                                     side = 2, font = 1, line = 0, cex = 1))
    dev.off()
  }
  
  
  ##-- Summary density map
  if("summary" %in% type){
    
    if(is.null(species)){
      message("You must provide one of 'bear', 'wolf' or 'wolverine' as the 'species' argument to be able plot the density map summary figure!")
    } else {
      ##-- Plot last year's density map
      # grDevices::pdf(file = file.path(path, paste0(name,"_Summary.pdf")),
      #                width = 8, height = 8, pointsize = 12)
      grDevices::png(filename = file.path(path, paste0(name,"_Summary.png")),
                     width = 8, height = 8, units = "in", pointsize = 12,
                     res = 300, bg = NA)
      
      graphics::par(mar = c(5,0,0,0))
      plot(sf::st_geometry(background), border = NA, col = "gray80")
      raster::image( density[[length(density)]], add = TRUE,
                     breaks = c(cuts, max(cuts)+1000),
                     col = col, legend = FALSE)
      plot( sf::st_geometry(background),
            border = "gray40", col = NA, add = TRUE)
      
      ##-- Add colour scale 
      raster::plot( density[[t]],
                    legend.only = T,
                    breaks = cuts,
                    col = col,
                    legend.width = 2,
                    axis.args = list( at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                      labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                      cex.axis = 1.2),
                    smallplot = c(0.80, 0.83, 0.25, 0.45),
                    legend.args = list(text = paste0("Individuals/", unit, " km2"),
                                       side = 2, font = 1, line = 0, cex = 1))
      
      ##-- Add km scale 
      addScale(x = 0.75, y = 0.25, size = 500000)

      ##-- Add species silhouette 
      addPNG( x = 0.8, y = 0.5, name = species, size = 0.2)
      
      ##-- Add flags 
      if(!is.null(labels)){addPopSize( x = x.labels, y = y.labels, labels = labels)}
      
      ##-- Add caption
      if(is.null(caption)){
      mtext(text = paste0("Density map and estimated ", engSpecies,
                          "\nabundance range estimated in ",
                          names(estimates)[length(density)]),
            side = 1,line = 2, adj = 0.5, cex = 1.2, font = 2)
      } else {
        mtext(text = caption, side = 1,line = 2, adj = 0.5, cex = 1.2, font = 2)
      }
      dev.off()
    }
  }
  
  
  ##-- Summary density map (Norwegian)
  if("summary_NOR" %in% type){
    
    if(is.null(species)){
      message("You must provide one of 'bear', 'wolf' or 'wolverine' as the 'species' argument to be able plot the density map summary figure!")
    } else {
      ##-- Plot last year's density map
      # grDevices::pdf(file = file.path(path, paste0(name,"_Summary.pdf")),
      #                width = 8, height = 8, pointsize = 12)
      grDevices::png(filename = file.path(path, paste0(name,"_Summary_NOR.png")),
                     width = 8, height = 8, units = "in", pointsize = 12,
                     res = 300, bg = NA)
      
      graphics::par(mar = c(5,0,0,0))
      plot(sf::st_geometry(background), border = NA, col = "gray80")
      raster::image( density[[length(density)]], add = TRUE,
                     breaks = c(cuts, max(cuts)+1000),
                     col = col, legend = FALSE)
      plot( sf::st_geometry(background),
            border = "gray40", col = NA, add = TRUE)
      

      ##-- Add colour scale 
      raster::plot( density[[t]],
                    legend.only = T,
                    breaks = cuts,
                    col = col,
                    legend.width = 2,
                    axis.args = list( at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                      labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                      cex.axis = 1.2),
                    smallplot = c(0.80, 0.83, 0.25, 0.45),
                    legend.args = list(text = paste0("Individer/", unit, " km2"),
                                       side = 2, font = 1, line = 0, cex = 1))
      
      ##-- Add km scale 
      addScale(x = 0.75, y = 0.25, size = 500000)
      
      ##-- Add species silhouette 
      addPNG( x = 0.8, y = 0.5, name = species, size = 0.2)
      
      ##-- Add flags 
      if(!is.null(labels)){addPopSize( x = x.labels, y = y.labels, labels = labels)}
      
      ##-- Add caption
      if(is.null(caption)){
      mtext(text = paste0("Kart som viser tetthet av ", norSpecies,
                          " med \nintervall for estimert antall ", norSpecies,
                          " i ", names(estimates)[length(density)]),
            side = 1,line = 2, adj = 0.5, cex = 1.2, font = 2)
        } else {
          mtext(text = caption, side = 1,line = 2, adj = 0.5, cex = 1.2, font = 2)
        }
      dev.off()
    }
  }
  
}
