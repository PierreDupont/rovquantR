#' @title Wolverine results processing
#'
#' @description
#' \code{processRovquantOutput_wolverine_report} processes both SCR and OPSCR outputs 
#' before combining them (together with previously published results) to create 
#' the figures and tables presented in the RovQuant wolverine report.
#' 
#' @param working.dir \code{Path} to the directory containing the \code{nimbleOutFiles} folder where MCMC outputs are stored. 
#' @param nburnin Number of initial, pre-thinning, MCMC bites to discard. Default value is 0.
#' @param niter Number of MCMC iterations to be used for density extraction. Default value is 100.
#' @param thin Thinning interval for collecting MCMC samples, corresponding to monitors. Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' @param thin2 Thinning interval for collecting MCMC samples, corresponding to the second, optional set of monitors2. Thinning occurs after the initial nburnin samples are discarded. Default value is 1.
#' @param extraction.res Resolution (in meters) to use for density extraction. Default value is 5000m.
#' @param Rmd.template (optional) \code{Path} to a custom .Rmd template to use instead of the default one provided in 'rovquantR'.
#' @param output.dir (optional) \code{Path} to the location of the final .html report. Default is to print in the 'reports' folder of the working directory.
#' @param overwrite \code{logical} Whether to silently overwrite (TRUE) or ask before overwriting potentially existing .html reports (FALSE).
#'
#'
#' @return This function creates the  following figures and tables:
#' \enumerate{
#' \item Summary figure in English and Norwegian (SCR)
#' \item Figure 1. Abundance time series (previous results + SCR)
#' \item Figure 2. Last year density map (SCR)
#' \item Figure 3. Survival (OPSCR)
#' \item Figure 4. Number of recruits (OPSCR)
#' \item Figure A.1. SkandObs + RovBase obs covariate (?)
#' \item Figure A.2. Regions & Counties map (?)
#' \item Figure A.3. Density maps time series (previous results + SCR)
#' \item Figure A.4. p0 structured (SCR)
#' \item Figure A.5. P0 other (SCR)
#' \item Table 1. last year N estimates based (SCR)
#' \item Table A.1. Number of NGS (OPSCR)
#' \item Table A.2. Number of individuals (OPSCR)
#' \item Table A.3. Number of dead recoveries (OPSCR)
#' \item Table A.4. Annual abundance (previous results + SCR)
#' \item Table A.5. Annual estimates in Norway (previous results + SCR)
#' \item Table A.6. Growth rate (previous results + SCR)
#' \item Table A.7. Demographic parameters (OPSCR)
#' \item Table A.8. Prop detected (previous results + SCR)
#' \item Table A.9. Density coefficient + sigma (SCR)
#' \item Table A.10. p0 structured coefficients (SCR)
#' \item Table A.11. p0 other coefficients (SCR)
#' }
#'
#' @author Pierre Dupont
#' 
#' @import sf 
#' @import raster
#' @import dplyr
#' @importFrom grDevices adjustcolor dev.off pdf png grey
#' @importFrom graphics axis abline par
#' @importFrom nimbleSCR scaleCoordsToHabitatGrid
#' @importFrom abind abind
#' @importFrom utils data
#' @importFrom xtable xtable
#' 
#' @rdname processRovquantOutput_wolverine_report
#' @export
processRovquantOutput_wolverine_report <- function(
    working.dir = NULL,
    nburnin = 0,
    niter = 100,
    thin = 1,
    thin2 = 1,
    extraction.res = 5000,
    years = NULL,
    overwrite = FALSE
){

  
  ##---- 2. OPSCR RESULTS PROCESSING -----
  
  ##-- Process the model output
  out2 <- processRovquantOutput_wolverine(
    working.dir = working.dir,
    nburnin = nburnin,
    niter = niter,
    thin = thin,
    thin2 = thin2,
    extraction.res = extraction.res,
    years = years,
    overwrite = overwrite)
  
  
  
  ##---- 3. COMBINE OPSCR & SCR RESULTS -----
  
  ## ------ 3.1. LOAD PREVIOUS RESULTS ------
  
  
  ## ------ Figure 1. Abundance time series (previous results + SCR) ------

  message("## Plotting combined abundance time series...")
  
  ##-- Plot N  
  grDevices::png(filename = file.path(working.dir,"figures/Abundance_TimeSeries_Report.png"),
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
  
  
  for(t in 1:(n.years-1)){
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
  
  
  
  
  
  
  
  
  ## ------ Figure A.3. Density maps time series (previous results + SCR) ------
  ## ------ Table A.4. Annual abundance (previous results + SCR) ------
  ## ------ Table A.5. Annual estimates in Norway (previous results + SCR) ------
  ## ------ Table A.6. Growth rate (previous results + SCR) -------
  ## ------ Table A.8. Prop detected (previous results + SCR) ------

  
  ##-- Combine model results
  out3 <- list()
    
  
  ##---- 4. PRINT REPORT -----
  
  if(print.report){
    
    ##-- Find the correct .rmd template for the report.
    if(is.null(Rmd.template)){
      Rmd.template <- system.file("rmd", "RovQuant_FullReport.Rmd", package = "rovquantR")
      if(!file.exists(Rmd.template)) {
        stop('Can not find the default .rmd template called "RovQuant_FullReport.Rmd". \n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
      } 
    }
    
    ##-- Find the directory to print the report.
    if(is.null(output.dir)){ output.dir <- file.path(working.dir, "reports") }
    
    ##-- Clean the data and print report
    rmarkdown::render(
      input = Rmd.template,
      params = list( species = out$SPECIES,
                     years = out$YEARS,
                     seasons = out$SEASONS,
                     date = out$DATE,
                     working.dir = working.dir),
      output_dir = output.dir,
      output_file = paste0("FullResults_", out$engSpecies, "_", out$DATE,".html"))
  }
}
