#' @title Brown bear dataset clean-up 
#'
#' @description
#' \code{cleanRovbaseData_bear} is an internal function called by the \code{cleanRovbaseData} function.
#' It identifies and loads the most recent Brown bear data downloaded from Rovbase
#' and conducts a set of cleaning steps that include:
#'  - removing un-identified samples
#'  - checking sex-assignment
#'  - removing problematic samples flagged by RovData.
#'
#' @name cleanRovbaseData_bear
#' 
#' @param species A \code{character} string with the name of the focal species
#'  ("bear", "wolf", or "wolverine").
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param sampling.months A \code{list} containing the months of interest, e.g. 
#' list(c(10:12),c(1:3)) for a sampling period extending from October to March. 
#' @param rename.list (Optional) A named \code{character} vector used to rename raw Rovbase files (.
#' @param data.dir the \code{path} pointing to the directory containing the raw 
#' data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. 
#' The cleaned data will be stored in the 'data' sub-folder. 
#' The associated \code{.html} report describing the content of the clean data will 
#' be placed in the 'report' sub-folder.
#' @param Rmd.template the \code{path} to the \code{.rmd} template to be used for
#'  cleaning the data. By default, the \code{.rmd} template provided with the 
#'  \code{rovquantR} package is used.  
#' @param overwrite A \code{logical} Whether previously existing clean data should
#' be overwritten or not (default = FALSE).
#' 
#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects
#'  for the species and period specified.
#' \item A \code{.html} report summarizing the data cleaning process. 
#' \item Additional \code{.png} images that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @importFrom dplyr distinct rename filter
#' @importFrom sf st_as_sf st_crs st_set_crs st_intersects

NULL
#' @rdname cleanRovbaseData_bear
#' @export
cleanRovbaseData_bear <- function(
    species = c("bear","wolf","wolverine"),
    years = NULL, 
    sex = c("Hunn","Hann"),
    sampling.months = NULL,
    rename.list = NULL,
    data.dir,
    working.dir,
    #print.report = TRUE,
    Rmd.template = NULL,
    overwrite = FALSE){
  
  ##----- 0. INITIAL CHECKS -----
  
  ##-- Make sure subfolders exist
  makeDirectories( path = working.dir,
                   two.sex = length(sex)>1,
                   show.dir = TRUE)
  
  
  ##-- Species
  if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjørn", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
    engSpecies <- "bear"
    norSpecies <- "Bjørn"
  } else {
    if(sum(grep("wolf", species, ignore.case = T))>0|sum(grep("ulv", species, ignore.case = T))>0){
      engSpecies <- "wolf"
      norSpecies <- "Ulv"
    } else {
      if(sum(grep("wolverine", species, ignore.case = T))>0|sum(grep("jerv", species, ignore.case = T))>0){
        engSpecies <- "wolverine"
        norSpecies <- "Jerv"
      } else {
       engSpecies <- norSpecies <- species 
      }
    }
  }
  
  
  ##-- Years
  if(is.null(years)){ years <- 2012:as.numeric(format(Sys.Date(), "%Y")) }
  
  
  ##-- Sampling months
  if(is.null(sampling.months)){
    if(engSpecies == "bear"){
      sampling.months <- list(4:11)
    } else {
      if(engSpecies == "wolf"){
        sampling.months <- list(c(10:12),c(1:4))
      } else {
        if(engSpecies == "wolverine"){
          sampling.months <- list(c(10:12),c(1:4))
        } else {
          stop("No default setting available for the monitoring period of this species. \n You must specify the monitoring season months through the 'sampling.months' argument.")
        }
      }
    }
  }
  

  ##-- Extract the date from the last .xlsx data file
  DATE <- getMostRecent( path = data.dir, pattern = "DNA.xlsx")
  
  
  ##-- Set file name for clean data
  fileName <- paste0("Data_", engSpecies, "_", DATE, ".RData")
  
  
  ##-- Check that a file with that name does not already exist to avoid overwriting
  if(!overwrite){
    if(any(file.exists(file.path(working.dir, "data", fileName)))){
      message(paste0("A file named '", fileName, "' already exists in: \n",
                     file.path(working.dir, "data")))
      message("Are you sure you want to proceed and overwrite existing clean data? (y/n) ")
      question1 <- readLines(n = 1)
      if(regexpr(question1, 'y', ignore.case = TRUE) != 1){
        message("Not overwriting existing files...")
        return(invisible(NULL))
      } else {
        message(paste0("Now overwriting '", fileName,"'.\n"))
      }
    }
  }
  
  
  ##-- Render .Rmd report
  if(is.null(Rmd.template)) {
    Rmd.template <- system.file("rmd", "RovQuant_CleaningReport.Rmd", package = "rovquantR")
    if (!file.exists(Rmd.template)) {
      stop("Can not find the Rmarkdown document to use for cleaning Rovbase.3.0 data.\n You must provide the path to the Rmarkdown template through the \"Rmd_template\" argument.")
    }
  }
  
  rmarkdown::render(input = Rmd.template,
                    params = list( species = species, 
                                   years = years,
                                   sampling.months = sampling.months,
                                   rename.list = rename.list,
                                   dir.in = data.dir,
                                   dir.out = working.dir,
                                   modDate = DATE),
                    output_dir = file.path(working.dir, "reports"), 
                    output_file = paste0("Data_", species, "_", DATE,".html"))
}
#   
#   ##-- Set-up list of parameters for the .Rmd report
#   params <- list()
#   
#   
#   ##----- 1. LOAD RAW NGS DATA -----
#   DNA <- suppressWarnings(readMostRecent( path = data.dir,
#                                           extension = ".xls",
#                                           pattern = "DNA")) %>%
#     ##-- Rename columns to facilitate manipulation
#     rename(., any_of(rename.list)) %>%
#     ##-- Filter to the focal species
#     filter(., Species %in% norSpecies) %>%
#     ##-- Remove any duplicates
#     distinct(., .keep_all = TRUE) %>%
#     ##-- Add some columns
#     dplyr::mutate( 
#       ##-- Add "Country" column
#       Country_sample = substrRight(County, 3),
#       ##-- Change date format
#       Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
#       ##-- Extract year
#       Year = as.numeric(format(Date,"%Y")),
#       ##-- Extract month
#       Month = as.numeric(format(Date,"%m")),
#       ##-- Extract sampling season
#       ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
#       ##-- Set all months in given sampling period to the same year)
#       Season = ifelse( Month < unlist(sampling.months)[1],
#                        Year,
#                        Year-1),
#       ##-- Fix unknown "Id"
#       Id = ifelse(Id %in% "", NA, Id),
#       ##-- Fix unknown "Sex"
#       Sex = ifelse(Sex %in% "Ukjent", NA, Sex)) %>%
#     ##-- Filter to the focal years
#     filter(., Year %in% years) 
#   
#   ##############################################################################
#   params$numDetDNA <- nrow(DNA)                     ## Number of NGS samples in Rovbase
#   params$numIdDNA <- length(unique(DNA$Id))        ## Number of NGS individuals in Rovbase 
#   params$tabDetDNA_sexyear <- table(DNA$Year, DNA$Sex, useNA = "always")
#   params$numUnknownSexDNA <- sum(is.na(DNA$Sex))
#   params$numUnknownCoordsDNA <- sum(is.na(DNA$East_UTM33))
#   params$numUnknownIDDNA <- sum(is.na(DNA$Id))
#   ##############################################################################
#   
#   
#   
#   ##----- 2. LOAD RAW DEAD RECOVERY DATA -----
#   
#   ##-- Load raw excel file imported from rovbase 
#   DR <- suppressWarnings(readMostRecent( path = data.dir,
#                                          extension = ".xls",
#                                          pattern = "dead")) %>%
#     ##-- Rename columns to facilitate manipulation
#     rename(., any_of(rename.list)) %>%
#     ##-- Filter to the focal species
#     filter(., Species %in% norSpecies) %>%
#     ##-- Remove any duplicates
#     distinct(., .keep_all = TRUE) %>%
#     ##-- Add some columns
#     dplyr::mutate( 
#       ##-- Add "Country" column
#       Country_sample = substrRight(County, 3),
#       ##-- Change date format
#       Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
#       ##-- Extract year
#       Year = as.numeric(format(Date,"%Y")),
#       ##-- Extract month
#       Month = as.numeric(format(Date,"%m")),
#       ##-- Extract sampling season
#       ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
#       ##-- Set all months in given sampling period to the same year)
#       Season = ifelse( Month < unlist(sampling.months)[1],
#                        Year,
#                        Year-1),
#       ##-- Fix unknown "Id"
#       Id = ifelse(Id %in% "", NA, Id),
#       ##-- Fix unknown "Sex"
#       Sex = ifelse(Sex %in% "Ukjent", NA, Sex)) %>%
#     ##-- Filter to the focal years
#     filter(., Year %in% years) 
#   
#   ##############################################################################
#   params$numDetDR <- nrow(DR)                ## Number of DEAD RECOVERIES in Rovbase
#   params$numIdDR <- length(unique(DR$Id))    ## Number of DEAD RECOVERIES ID in Rovbase
#   params$tabDetDR_sexyear <- table(DR$Year, DR$Sex)
#   params$numUnknownSexDR <- sum(is.na(DR$Sex))
#   params$numUnknownCoordsDR <- sum(is.na(DR$East_UTM33))
#   params$numUnknownIDDR <- sum(is.na(DR$Id))
#   ##############################################################################
#   
#   
#   
#   ##----- 3. CHECK DEAD RECOVERIES -----
# 
#   ##-- Make sure all dead recoveries in DNA are in DR
#   check1 <- all(DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$DNAID) 
#   params$probs_DR_in_DNA <- NULL
#   if(!check1){
#     tmp <- DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"]
#     params$probs_DR_in_DNA <- tmp[!tmp %in% DR$DNAID]
#   } 
#   
#   tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
#   test <- anti_join(tmp,DR, by = names(tmp)[names(tmp) %in% names(DR)])
#   
#   
#   ##-- Make sure that only "Dead recoveries" are in DR 
#   check2 <- all(substr(DR$RovbaseID,1,1) %in% "M")
#   params$probs_DNA_in_DR <- NULL
#   if(!check2){
#     params$probs_DNA_in_DR <- DR$DNAID[!substr(DR$RovbaseID,1,1) %in% "M"]
#   }
#   
# 
#   
#   ##----- 4. MERGE DATASETS -----
#   DATA <- full_join(DNA, DR, by = names(DNA)[names(DNA) %in% names(DR)]) 
#   # DATA <- merge( DR, DNA, 
#   #                by = c("Id","RovbaseID","DNAID","Species", "Sex","Date","East_UTM33","North_UTM33", "County"),
#   #                all = TRUE)
#   
#   
#    
#   ##----- 5. CLEAN UP DATA -----
#   # index <- DATA$Month < unlist(sampling.months)[1]
#   # index[is.na(index)] <- FALSE
#   # DATA$Year[index] <- DATA$Year[index] - 1
#    
#   
#   ##-- Determine Death and Birth Years
#   DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
#   DATA$RovbaseID <- as.character(DATA$RovbaseID)
#   DATA$Death <- NA
#   DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
#   DATA$Birth <- DATA$Death - DATA$Age
#   
#   
#   ##############################################################################
#   params$noID <- sum(is.na(DATA$Id))              ## Number of samples without ID
#   params$noDate <- sum(is.na(DATA$Year))          ## Number of samples without Date
#   params$noCoords <- sum(is.na(DATA$East_UTM33))  ## Number of samples without Coords  
#   ##############################################################################
#   
#   ##-- Filter out unusable samples
#   DATA <- DATA %>%
#     dplyr::filter(., 
#                   ##-- Filter out samples with no ID
#                   !is.na(Id),
#                   ##-- Filter out samples with no Coordinates
#                   !is.na(East_UTM33),
#                   ##-- Filter out samples with no dates  
#                   !is.na(Year)) %>%
#     droplevels(.)
#   
#   
#   
#   ##----- 6. CHECK SEX ASSIGNMENT -----
#   
#   ID <- unique(as.character(DATA$Id))
#   DATA$Sex <- as.character(DATA$Sex)
#   doubleSexID <- IdDoubleSex <- NULL
#   counter <- 1
#   for(i in 1:length(ID)){
#     ##-- Subset data to individual i
#     tmp <- DATA$Sex[DATA$Id == ID[i]]
#     
#     ##-- Number of times individual i was assigned to each sex
#     tab <- table(tmp[tmp %in% c("Hunn","Hann")])
#     
#     ##-- If conflicting sexes (ID identified as both "Hunn" and "Hann")
#     if(length(tab) == 2){
#       ##-- If ID assigned the same number of times to the 2 sexes, assign to Ukjent
#       if(tab[1] == tab[2]){
#         DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"
#       } else {
#         ##-- Otherwise pick the most common sex
#         DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
#       }
#       IdDoubleSex[counter] <- ID[i]
#       counter <- counter + 1
#     }
#     
#     ##-- If only one of "Hunn" or "Hann" registered
#     if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
#     
#     ##-- If anything else registered : "Ukjent"
#     if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
#     
#     doubleSexID[i] <- length(tab)
#   }#i
#   
#   
#   
#   ##----- 7. SPLIT DATA -----
#   
#   ##-- Split DATA into alive and dead.recovery datasets
#   alive <- DATA[is.na(DATA$Death), ]
#   dead.recovery <- DATA[!is.na(DATA$Death), ]
#   
#   
#   ##-- Add earlier detection index
#   alive$detected.earlier <-
#     unlist(lapply(1:nrow(alive),
#                   function(i){
#                     this.id <- alive[i,"Id"]
#                     this.date <- alive[i,"Date"]
#                     any(alive$Id %in% this.id & alive$Date < this.date)
#                   }))
#   
#   dead.recovery$detected.earlier <-
#     unlist(lapply(1:nrow(dead.recovery),
#                   function(i){
#                     this.id <- dead.recovery[i,"Id"]
#                     this.date <- dead.recovery[i,"Date"]
#                     any(alive$Id %in% this.id & alive$Date < this.date)
#                   }))
#   
#   
#   
#   ##----- 8. REMOVE FLAGGED SAMPLES -----
#   
#   ##-- Load most recent "flagged" file from HB
#   flagged <- readMostRecent( 
#     path = data.dir,
#     extension = ".csv",
#     pattern = "dna_bear_to_remove", 
#     fileEncoding = "Latin1") 
#   
#   
#   ##-- Remove flagged samples 
#   keep.alive <- !alive$Barcode_sample %in% flagged$Strekkode
#   alive <- alive[keep.alive, ]
#   keep.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
#   dead.recovery <- dead.recovery[keep.dead, ]
#   dead.recovery$Missing <- NA
#   dead.recovery$Individ <- NA
#   
#   
#   
#   ##----- 9. TURN INTO .sf OBJECTS -----
#   
#   ##-- Turn into sf points dataframe
#   alive <- sf::st_as_sf( x = alive,
#                          coords = c("East_UTM33","North_UTM33")) %>%
#     sf::st_set_crs(.,sf::st_crs(32633)) 
#   
#   
#   ##-- Intersect and extract country name
#   alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
#   
#   
#   ##-- Turn into sf points dataframe
#   dead.recovery <- sf::st_as_sf( x = dead.recovery,
#                                  coords = c("East_UTM33","North_UTM33")) %>%
#     sf::st_set_crs(.,sf::st_crs(32633))
#   
#   
#   ##-- Intersect and extract country name
#   dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]
#   
#   
#   
#   ##----- 10. IDENTIFY ISSUES -----
#   
#   ##-- Multiple deaths
#   ##-- Identify and count individuals dead "more than once"
#   ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
#   multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
#   for(i in 1:length(ID)){
#     tmp <- dead.recovery[dead.recovery$Id == ID[i], ] 
#     ##-- Multiple death dates
#     if(length(unique(tmp$Date)) > 1){
#       multiDeathDate <- c(multiDeathDate, ID[i])
#     }
#     ##-- Multiple death years
#     if(length(unique(tmp$Year)) > 1){
#       multiDeathYear <- c(multiDeathYear, ID[i])
#     }
#     ##-- Multiple death locations
#     if(length(unique(tmp$East)) > 1 | length(unique(tmp$North)) > 1){
#       multiDeathLocs <- c(multiDeathLocs, ID[i])
#     }
#   }#i
#   
#   ##-- Remove individuals that died more than once
#   dead.recovery$Id <- as.character(dead.recovery$Id)
#   IdDoubleDead <- names(table(dead.recovery$Id))[table(dead.recovery$Id) > 1]
#   if(length(IdDoubleDead) > 0){
#     for(i in IdDoubleDead){
#       ##-- Identify repeated deaths
#       tmp <- which(dead.recovery$Id %in% i) 
#       ##-- Try to keep death with known death cause 
#       tmp2 <- which(!is.na(dead.recovery$DeathCause_2[tmp]))[1]
#       if(length(tmp2) == 0){tmp <- tmp[-1]} else {tmp <- tmp[!tmp %in% tmp2]}
#       ##-- Remove repeated deaths
#       dead.recovery <- dead.recovery[-tmp, ]
#     }#i
#   }#if
#   
#   ##-- Detections after death
#   id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
#   ghosts <- unlist(lapply(id.list, function(id){
#     out <- NULL
#     try({
#       if(id %in% dead.recovery$Id){
#         mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
#         this.alive <- alive[alive$Id == id, ]
#         ##-- Was it detected alive in any season after death?
#         temp <- this.alive[this.alive$Year > mort.year, ]
#         if(length(temp) > 0){
#           out <- rownames(temp)
#           names(out) <- id
#         }##-- FLAG THOSE FOR HENDRIK
#       }
#     }, silent = TRUE)
#     return(out)
#   }))
#   samples.to.remove <- unlist(ghosts)
#   
#   ##-- Remove flagged NGS detections after dead recovery
#   alive <- alive[!rownames(alive) %in% samples.to.remove, ]
#   
#   
#   
#   ##----- 11. RENDER .html REPORT -----
#   
#   if(print.report){
#     
#     ##-- Find the .rmd template
#     if(is.null(Rmd.template)){
#       Rmd.template <- system.file("rmd", "RovQuant_CleaningReport.Rmd", package = "rovquantR")
#       if(!file.exists(Rmd.template)) {
#         stop('Can not find the Rmarkdown document to use for cleaning Rovbase.3.0 data.\n You must provide the path to the Rmarkdown template through the "Rmd.template" argument.')
#       } 
#     }
#     
#     ##-- Render the .html report
#     reportName <- paste0("Data_", engSpecies, "_", DATE,".html")
#     
#     message(paste0("\nRendering report '", reportName, "' in: \n",
#                    file.path(working.dir, "data")))
#     
#     rmarkdown::render(
#       input = Rmd.template,
#       params = params,
#       output_dir = file.path(working.dir, "reports"),
#       output_file = reportName)
#   }
#   
#   
#   
#   ##----- 12. SAVE CLEAN DATA -----
#   
#   ##-- Save .RData file
#   message(paste0("\nSaving file '", fileName, "' in: \n",
#                  file.path(working.dir, "data")))
#   save( alive, 
#         dead.recovery,
#         IdDoubleSex,
#         file = file.path(working.dir, "data", fileName))
#   
#   
#   ##-- If explicit return also data as a list
#   # if(silent){
#   #   
#   #   dat.list <- list(
#   #     "alive" = alive, 
#   #     "dead.recovery" = dead.recovery,
#   #     "IdDoubleSex" = IdDoubleSex)
#   #   
#   #   return(dat.list)
#   # }