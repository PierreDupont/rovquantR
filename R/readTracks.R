#' Read and combine multiple .shp files containing GPS tracks of recorded search effort
#'
#' The \code{readTracks} function identifies and loads multiple .shp files, before combining them using rbind.
#' 
#' @name readTracks
#'
#' @param data.dir A character string denoting the data directory path.
#' 
#' @return The data loaded
#'
#' @author Pierre Dupont
#'
#' @importFrom sf read_sf st_length st_coordinates st_centroid
#' @importFrom dplyr mutate filter
#' 
NULL
#' @rdname readTracks
#' @export
readTracks <- function(data.dir, 
                       years = 1950:format(Sys.Date(), "%Y"),
                       sampling.months = 1:12){
  
  ##-- Check that a "Tracks" folder exists
  if(!dir.exists(file.path(data.dir,"Tracks"))){
    stop(paste0("Did not find GPS search tracks in '", data.dir,"'.\nPlease store search effort shapefiles in a folder named 'Tracks'."))
  }
  
  ##-- List all files with the requested extension (including in sub-directories)
  infiles <- list.files( path = file.path(data.dir, "Tracks"),
                         pattern = ".shp",
                         recursive = TRUE,
                         ignore.case = TRUE)
  if(length(infiles) < 1){stop(paste0("No file was found with a '.shp' extension in the requested directory."))}
  
  ##-- Load all GPS tracks in a list
  data <- list()
  for(f in 1:length(infiles)){
    message(paste0('Loading file ', infiles[f], '...\n'))
    
    ##-- Load and filter GPS search tracks
    data[[f]] <- sf::read_sf(file.path(data.dir, "Tracks", infiles[f])) %>%
      ##-- Process dates
      dplyr::mutate( Dato = as.POSIXct(strptime(Dato, "%Y-%m-%d")),
                     Mth = as.numeric(format(Dato,"%m")),
                     Yr = as.numeric(format(Dato,"%Y")),
                     Year = ifelse( Mth < unlist(sampling.months)[1], Yr-1,Yr)) %>%
      ##-- Filter out irrelevant tracks
      dplyr::filter( Helikopter == "0",      ## Remove helicopter tracks
                     # Jerv == "1",          ## [CHECK] should we keep wolf tracks only?
                     Year %in% years & Mth %in% unlist(sampling.months)) %>% ## Keep tracks during sampling season only
      ##-- Extract track lengths & centroids
      dplyr::mutate( Length = sf::st_length(., byid = T),
                     Centroidx = sf::st_coordinates(sf::st_centroid(.))[ ,1])
  }#f
  
  ##-- combine data files
  data <- do.call(rbind, data)
  
  ##-- Remove duplicates
  ##-- Find & filter out duplicates based on person, distance and date.
  df <- data.frame( Dato = data$Dato,
                    Year = data$Year,
                    Person = data$Person,
                    Length = data$Length,
                    Centroidx = data$Centroidx)
  data <- data[!duplicated(df), ]  
  
  ##-- Output
  return(data)
}
