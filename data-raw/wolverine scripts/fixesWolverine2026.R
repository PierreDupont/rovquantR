#' Wolverine Package Analysis 2026 - Problems & Fixes
#'
#'
#' FIXES:
#' 1. Still uses wrong assignment of monitoring seasons for skandobs and rovbase.
#' It should be from the start of the monitoring season and not from December
#'   ==> fix already in the script 
#'   ==> DONE in branch 'wolverine2'
#'  
#' 2. Problem generating 'RovQuant_DataReport.Rmd'
#'   ==> "Quitting from RovQuant_DataReport.Rmd:414-569 [nimble data]" 
#'  
#' 3. Fix monitoring season names throughout the markdown reports
#'   ==> use 'seasons' instead of years for all titles, plot labels, tables, etc... 
#'   
#' 4. Fix the set-up for SCR models directories 
#'   ==> remove the unused folder with the year as a name 
#'   ==> DONE in branch 'wolverine2'
#' 
#' 
#' IMPROVEMENT
#' 1. Use fixed extent instead of habitat based on detections 
#'   ==> DONE in branch 'wolverine2'
#' 
#' 2. Use 'nimbleSCR' format for everything:
#'   ==> Use a sf spatial grid dataframe for all habitat and detector characteristics 
#'  (i.e. get rid of rasters and spatial points)
#'  
#' 3. Distance to roads : 
#' - no need to remove values where habitat is NA before assigning to each detector; 
#' - this is what creates NAs when assigning to detectors.
#' 
#' 
#' OTHER TASKS:
#' 1. Check w/ CM if I can delete branches "wolf" and "wolfdevel"
#' 2. Continue uniformization/cleaning of the wolf functions
#' 