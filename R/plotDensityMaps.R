#' Plot population density maps
#'
#' Function to plot (and print as \code{.pdf} files) population density maps from SCR or OPSCR models.
#'
#' The \code{plotDensityMaps} function uses the outputs from the \code{getDensityInput} and \code{getDensity} functions as inputs.
#'
#' @param input Density inuts, as prepared by the \code{getDensityInput} function.
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
#' @importFrom raster res image extent
#' @importFrom graphics plot layout par segments mtext text 
#' @importFrom sf st_geometry 
#' @importFrom grDevices pdf colorRampPalette
#'
#' @rdname plotDensityMaps
#' @export
plotDensityMaps <- function( 
    input,
    estimates,
    unit = 100,
    mask,
    background = NULL,
    type = c("all"),# "last.year", "time.series"),
    path = file.path(getwd(), "figures"),
    image = NULL,
    name = "UD_Density")
{
  
  ##-- Initial checks
  if(is.null(mask)){ mask <- input$regions.r }
  if(!inherits(estimates,"list")){ estimates <- list(estimates) }
  if(is.null(background)){background <- COUNTRIES}
  
  ##-- Convert densities to the desired density unit (usually inds.100km-2)
  conversionFactor <- unit/( raster::res(input$regions.r)[1]/1000)^2
  estimates <- lapply(estimates, function(x) x$MeanCell * conversionFactor)
  
  ##-- Rasterize and mask density maps 
  density <- list()
  for(t in 1:length(estimates)){
    density[[t]] <- input$regions.r
    density[[t]][] <- NA
    density[[t]][!is.na(input$regions.r[])] <- estimates[[t]]
    density[[t]][is.na(mask[])] <- NA
  }#t
  names(density) <- names(estimates)
  
  ##-- Density maps time series
  if(type %in% c("time.series","all")){
    
    ##-- layout
    L <- length(density)
    if(L < 6){ nrows <- 1 } else{
      if(L < 13){ nrows <- 2 } else {
        if(L < 22){ nrows <- 3 } else {
          if(L < 33){ nrows <- 4 } else {
            nrows <- 5
          }}}}
    ncols <- ceiling(L/nrows)
    
    # grDevices::pdf(file = file.path(path, paste0(name,"_TimeSeries.pdf")),
    #                width = ncols*2, height = nrows*4)
    
    grDevices::png(filename = file.path(path, paste0(name,"_TimeSeries.png")),
                   width = ncols*2, height = nrows*4,
                   units = "in", pointsize = 12,
                   res = 300, bg = NA)
    
    ##-- Set color scale
    max <- max(unlist(lapply(density, function(x) max(x[], na.rm = T))))
    cuts <- seq(0, max, length.out = 100) ##-- set breaks
    colfunc <- grDevices::colorRampPalette(c("white", "slateblue", "yellow", "orange", "red", "red"))
    col <- colfunc(100)
    
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
    legend.x <- raster::extent(background)[1] + 0.77*diff(raster::extent(background)[1:2])
    legend.y <- raster::extent(background)[3] + 0.4*diff(raster::extent(background)[3:4])
    
    ##-- Plot density maps
    graphics::par(mar = c(0,0,0,0))
    
    for(t in 1:length(density)){
      ##-- Plot density
      plot(sf::st_geometry(background), border = NA, col = "gray80")
      raster::image( density[[t]], add = TRUE,
                     breaks = c(cuts, max(cuts)+1000),
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
          x = legend.x-0.05*diff(raster::extent(background)[1:2]),
          y = legend.y,
          labels = "500 km", srt = 90, cex = 1.4)
        raster::plot( density[[t]],
                      legend.only = T, breaks = cuts,
                      col = col, legend.width = 2,
                      axis.args = list( at = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                        labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
                                        cex.axis = 1.2),
                      smallplot = c(0.85, 0.9, 0.2, 0.6), 
                      legend.args = list(text = paste0("Individuals/", unit, " km2"),
                                         side = 2, font = 1, line = 0, cex = 1))
        ######----- NEED TO FIX LEGEND TEXT 
        ######----- expression("Individuals/100 km"^ 2)
      }#if
      # ##-- Export rasters
      # writeRaster(density[[t]],
      #             file.path(path,
      #                       paste0("DensityRaster100km2", years[t],".tif")),
      #             overwrite = TRUE)
    }#t
    dev.off()
  }
  
  
  ##-- Last year's density map
  if(type %in% c("last.year","all")){
    
    # grDevices::pdf(file = file.path(path, paste0(name,"_LastYear.pdf")),
    #                width = 8, height = 8)
    
    grDevices::png(filename = file.path(path, paste0(name,"_LastYear.png")),
        width = 8, height = 8, units = "in", pointsize = 12,
        res = 300, bg = NA)
    
    t <- length(density)
    
    graphics::par(mar = c(0,0,0,0))
    plot(sf::st_geometry(background), border = NA, col = "gray80")
    raster::image( density[[t]], add = TRUE,
                   breaks = c(cuts, max(cuts)+1000),
                   col = col, legend = FALSE)
    plot( sf::st_geometry(background),
          border = "gray40", col = NA, add = TRUE)
    
    ##-- Add year if available
    if(!is.null(names(estimates))){
      mtext(text = names(estimates)[t],
            side = 1, line =  -25,
            adj = 0.25, cex = 3, font = 2)
    }
    
    ##-- Add legend
    legend.x <- raster::extent(background)[1] + 0.8*diff(raster::extent(background)[1:2])
    legend.y <- raster::extent(background)[3] + 0.3* diff(raster::extent(background)[3:4])
    
    graphics::segments(
      x0 = legend.x, x1 = legend.x,
      y0 = legend.y-250000, y1 = legend.y + 250000,
      col = "gray30", lwd = 4, lend = 2)
    graphics::text(
      x = legend.x-0.05*diff(raster::extent(background)[1:2]),
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
    ######----- NEED TO FIX LEGEND TEXT 
    # graphics::segments(x0 = 800000, x1 = 800000,
    #          y0 = 6650000, y1 = 6650000 + 500000,
    #          col = "gray30", lwd = 4, lend = 2)
    # graphics::text(760000, 6650000+500000/2, labels = "500 km", srt = 90, cex = 1.2)
    # raster::plot( density[[t]],
    #       legend.only = T,
    #       breaks = cuts,
    #       col = col,
    #       legend.width = 2,
    #       axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
    #                        labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
    #                        cex.axis = 1.2),
    #       smallplot = c(0.75, 0.78, 0.2, 0.4),
    #       legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
    #                          side = 2, font = 1, line = 0, cex = 1.2))
    dev.off()
  }
  
  
  ##-- Summary density map
  if(type %in% c("summary","all")){
    
     grDevices::png(filename = file.path(path, paste0(name,"_Summary.png")),
                   width = 8, height = 8, units = "in", pointsize = 12,
                   res = 300, bg = NA)
    
    t <- length(density)
    
    graphics::par(mar = c(0,0,0,0))
    plot(sf::st_geometry(background), border = NA, col = "gray80")
    raster::image( density[[t]], add = TRUE,
                   breaks = c(cuts, max(cuts)+1000),
                   col = col, legend = FALSE)
    plot( sf::st_geometry(background),
          border = "gray40", col = NA, add = TRUE)
    

    
    ##-- If species is specified; insert the corresponding image file
    if(!is.null(species)){
      
    ##-- Add year if available
    if(!is.null(names(estimates))){
      mtext(text = paste0("Density map and ranges of abundance /nestimated for ", species," in " names(estimates)[t]), 
            side = 1, line =  -25,
            adj = 0.25, cex = 3, font = 2)
    }
    
    
      myImage <- system.file("extdata", "file.csv", package = "mypackage")
    }
    
    
    ##-- Add legend
    legend.x <- raster::extent(background)[1] + 0.8*diff(raster::extent(background)[1:2])
    legend.y <- raster::extent(background)[3] + 0.3* diff(raster::extent(background)[3:4])
    
    graphics::segments(
      x0 = legend.x, x1 = legend.x,
      y0 = legend.y-250000, y1 = legend.y + 250000,
      col = "gray30", lwd = 4, lend = 2)
    graphics::text(
      x = legend.x-0.05*diff(raster::extent(background)[1:2]),
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
    ######----- NEED TO FIX LEGEND TEXT 
    # graphics::segments(x0 = 800000, x1 = 800000,
    #          y0 = 6650000, y1 = 6650000 + 500000,
    #          col = "gray30", lwd = 4, lend = 2)
    # graphics::text(760000, 6650000+500000/2, labels = "500 km", srt = 90, cex = 1.2)
    # raster::plot( density[[t]],
    #       legend.only = T,
    #       breaks = cuts,
    #       col = col,
    #       legend.width = 2,
    #       axis.args = list(at = round(seq(0, max-0.05, length.out = 4), digits = 1),
    #                        labels = round(seq(0, max-0.05, length.out = 4), digits = 1),
    #                        cex.axis = 1.2),
    #       smallplot = c(0.75, 0.78, 0.2, 0.4),
    #       legend.args = list(text = expression(paste("Individuals/100 km"^ 2)),
    #                          side = 2, font = 1, line = 0, cex = 1.2))
    dev.off()
  }
  
}
