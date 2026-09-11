#' Wolverine Package Analysis 2026 - Problems & Fixes
#'
#' Still uses 'sex(analysis)' for the sex assignment to match previous analysis (same as wolves).
#'  ==> Switch to 'sex (individ)' for next analysis.
#'  
#' Still uses wrong assignment of monitoring seasons for skandobs and rovbase.
#' It should be from the start of the monitoring season and not from December
#'  ==> fix already in the script (commented out)
#'   
#' Could not get readTracks() to reproduce the previous analysis script
#' Could not get readMostRecent() to read .tif files (SNOW AND ROADS)
#' Could not get readMultiples() to work all samples from Rovbase
#' 'Sample_type' used from rovbase samples do not use "Loepeblod", "Vev" to match previous analysis. 
#' Could not get assignSearchTracks() to reproduce the previous analysis script 
#' 
#' 
#' FURTHER IMPROVEMENTS:
#' 
#' 1.use 'nimbleSCR' format for everything:
#'  1.1. use a sf spatial grid dataframe for all habitat and detector characteristics (i.e. get rid of rasters and spatial points)
#'  
#'  2. distance to roads : no need to remove values where habitat is NA before assigning to each detector; this is what creates NAs when assigning to detectors.