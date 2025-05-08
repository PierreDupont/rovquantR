#' @title Brown bear Data set clean-up.
#'
#' @description
#' \code{cleanRovbaseData_bear} identifies and loads the most recent Brown bear data
#'  extracted from Rovbase.no in the specified \code{data.dir} location, and 
#'  conducts a set of data cleaning steps.  @seealso [cleanRovbaseData()] for more details on the data cleaning process.
#'
#' @name cleanRovbaseData_bear
#' 
#' @param data.dir the \code{path} pointing to the directory containing the raw 
#' data from Rovbase.
#' @param working.dir the \code{path} pointing to the working directory. By default,
#'  the cleaned data will be stored in a subfolder of this working directory called 'data'.
#' @param years A \code{numeric} vector containing the years of interest. 
#' Only data for those years will be cleaned and returned.
#' @param sex A \code{character} vector containing the sex of interest. 
#' Can be "Hunn" for females or "Hann" for males. Default is both sexes (c("Hunn","Hann")).
#' @param sampling.months (Optional) A \code{list} containing the sampling period months.
#' If the sampling period overlaps two calendar years, the list should contain one element per year.
#' (e.g. samplingMonths <- list(c(11,12), c(1,2,3,4))) for a sampling period extending from November to April of the following year.
#' @param rename.list (Optional) A named \code{character} vector used to rename columns in the raw Rovbase files.


#' @return This function returns:
#' \enumerate{
#' \item A \code{.RData} file with the clean NGS and dead recovery data objects 
#' for the species and period specified. The clean data file is saved as an \code{.RData} 
#' file named using the species name and the date of extraction of the raw Rovbase data 
#' to facilitate replicability (e.g. 'CleanData_bear_2024-08-10.RData').
#' \item A \code{.html} report summarizing the data cleaning process. 
#' The \code{.RData} report is using the same naming convention as the clean \code{.RData} (e.g. 'CleanData_bear_2024-08-10.html').
#' \item Additional \code{.png} images that can be reused somewhere else.
#' }
#'
#' @author Pierre Dupont
#' 
#' @import ggplot2 

NULL
#' @rdname cleanRovbaseData_bear
#' @export
cleanRovbaseData_bear <- function(
    ##-- paths
  data.dir,
  working.dir,
  
  ##-- data
  years, 
  sex = c("female","male"),
  sampling.months = NULL,
  rename.list,
) {
  
  ##----- 1. INITIAL CHECKS -----
  
  ##-- Species
  species <- "bear"
  engSpecies <- "Brown bear"
  norSpecies <- "Bjørn"
  
  ##-- Sampling months
  if (is.null(sampling.months)) { sampling.months <- list(4:11) }
  
  
  
  ##----- 2. CLEAN THE DATA -----
  
  ##-- Load pre-processed habitat shapefiles
  data(COUNTRIES, envir = environment()) 
  # ##-- List months
  # months = c("January","February","March","April","May","June",
  #            "July","August","September","October","November","December")
  
  ##-- data info
  DATE <- getMostRecent(path = data.dir, pattern = "DNA")
  
  
  ##-- NGS data
  DNA <- suppressWarnings(readMostRecent( path = data.dir,
                                          extension = ".xls",
                                          pattern = "DNA")) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    filter(., Species %in% norSpecies) %>%
    ##-- Remove any duplicates
    distinct(., .keep_all = TRUE) %>%
    ##-- Add some columns
    dplyr::mutate( 
      ##-- Add "Country" column
      Country_sample = substrRight(County, 3),
      ##-- Change date format
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      ##-- Extract year
      Year = as.numeric(format(Date,"%Y")),
      ##-- Extract month
      Month = as.numeric(format(Date,"%m")),
      ##-- Extract sampling season
      ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
      ##-- Set all months in given sampling period to the same year)
      Season = ifelse( Month < unlist(sampling.months)[1],
                       Year,
                       Year-1),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", NA, Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    filter(., Year %in% years) 
  
  ##-- Number of NGS samples
  NGS_samples <- table(DNA$Sex, DNA$Year, useNA = "ifany")
  NGS_samples <- rbind(NGS_samples, "Total" = colSums(NGS_samples))
  NGS_samples <- cbind(NGS_samples, "Total" = rowSums(NGS_samples))
  write.csv(NGS_samples,
            file = file.path(working.dir, "tables",
                             paste0(species, "_NGS samples_", years[1]," to ", years[length(years)], ".csv")))
  
  ##-- Number of individuals detected alive
  NGS_ids <- apply(table(DNA$Sex, DNA$Year, DNA$Id, useNA = "ifany"),
               c(1,2), function(x)sum(x>0))
  NGS_ids <- rbind(NGS_ids,
               "Total" = apply(table(DNA$Year,DNA$Id, useNA = "ifany"),
                               1, function(x)sum(x>0)))
  NGS_ids <- cbind(NGS_ids,
               "Total" = c(apply(table(DNA$Sex,DNA$Id, useNA = "ifany"),
                                 1, function(x)sum(x>0)),length(unique(DNA$Id))))
  write.csv(NGS_ids,
            file = file.path(working.dir, "tables",
                             paste0(species, "_NGS ids_", years[1]," to ", years[length(years)], ".csv")))

  
  
  
  ## ---------------------------------------------------------------------------
  ## -- MOVE TO .RMD -----------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ##-- NGS sample table
  NGS_samples <- read.csv( file.path(working.dir, "tables", paste0(species, "_NGS samples_", years[1]," to ", years[length(years)], ".csv")),
                           check.names = FALSE)
  knitr::kable(NGS_samples, align = "lc",
        caption = "Number of NGS samples collected by year and sex") %>%
    kableExtra::kable_styling(full_width = F)
  
  
  ##-- NGS id table data
  NGS_ids <- read.csv( file.path(working.dir, "tables", paste0(species, "_NGS ids_", years[1]," to ", years[length(years)], ".csv")),
                           check.names = FALSE)
  knitr::kable(NGS_ids, align = "lc",
         caption = "Number of individuals detected through NGS by year and sex") %>%
    kableExtra::kable_styling(full_width = F)
  ## ---------------------------------------------------------------------------
  ## -- MOVE TO .RMD -----------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  
  
  
  
  ## ----DR data ---------------------------------------------------------------
  ##-- Load raw excel file imported from rovbase 
  DR <- suppressWarnings(readMostRecent( path = data.dir,
                                         extension = ".xls",
                                         pattern = "dead")) %>%
    ##-- Rename columns to facilitate manipulation
    rename(., any_of(rename.list)) %>%
    ##-- Filter to the focal species
    filter(., Species %in% norSpecies) %>%
    ##-- Remove any duplicates
    distinct(., .keep_all = TRUE) %>%
    ##-- Add some columns
    dplyr::mutate( 
      ##-- Add "Country" column
      Country_sample = substrRight(County, 3),
      ##-- Change date format
      Date = as.POSIXct(strptime(Date, "%Y-%m-%d")),
      ##-- Extract year
      Year = as.numeric(format(Date,"%Y")),
      ##-- Extract month
      Month = as.numeric(format(Date,"%m")),
      ##-- Extract sampling season
      ##-- (for sampling periods spanning over two calendar years (wolf & wolverine)
      ##-- Set all months in given sampling period to the same year)
      Season = ifelse( Month < unlist(sampling.months)[1],
                       Year,
                       Year-1),
      ##-- Fix unknown "Id"
      Id = ifelse(Id %in% "", NA, Id),
      ##-- Fix unknown "Sex"
      Sex = ifelse(Sex %in% "Ukjent", NA, Sex),
      Sex = ifelse(Sex %in% "Hunn", "female", Sex),
      Sex = ifelse(Sex %in% "Hann", "male", Sex)) %>%
    ##-- Filter to the focal years
    filter(., Year %in% years) 
  
  
  ##-- Number of NGS samples
  DR_samples <- table(DR$Sex, DR$Year, useNA = "ifany")
  DR_samples <- rbind(DR_samples, "Total" = colSums(DR_samples))
  DR_samples <- cbind(DR_samples, "Total" = rowSums(DR_samples))
  write.csv(DR_samples,
            file = file.path(working.dir, "tables",
                             paste0(species, "_DR samples_", years[1]," to ", years[length(years)], ".csv")))
  
  ##-- Number of individuals detected alive
  DR_ids <- apply(table(DR$Sex, DR$Year, DR$Id, useNA = "ifany"),
               c(1,2),
               function(x)sum(x>0))
  DR_ids <- rbind(DR_ids,
               "Total" = apply(table(DR$Year,DR$Id, useNA = "ifany"),
                               1,
                               function(x)sum(x>0)))
  DR_ids <- cbind(DR_ids,
               "Total" = c(apply(table(DR$Sex,DR$Id, useNA = "ifany"),
                                 1,
                                 function(x)sum(x>0)),length(unique(DR$Id))))
  write.csv(DR_ids,
            file = file.path(working.dir, "tables",
                             paste0(species, "_DR ids_", years[1]," to ", years[length(years)], ".csv")))
  
  
  
  
  ## ---------------------------------------------------------------------------
  ## -- MOVE TO .RMD -----------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## ----DR sample table, echo = F, collapse = TRUE----------------------------------------------------------------
  kable( DR_samples,
         align = "lc",
         caption = "Number of dead recoveries by year and sex") %>%
    kable_styling(full_width = F)
  
  
  
  ## ----DR id table data, echo = F, collapse = TRUE---------------------------------------------------------------
  kable( ids,
         align = "lc",
         caption = "Number of individuals identified from dead recoveries by year and sex")%>%
    kable_styling(full_width = F)
  # samples_kable <- kable(samples, align = "lc",
  #       caption = "Number of individuals detected through NGS by year and sex",
  #       format = "latex")
  # 
  # ids_kable <- kable(ids, align = "lc",
  #       caption = "Number of individuals detected through NGS by year and sex",
  #       format = "latex")
  # 
  # cat(c("\\begin{table}[h] \\centering ", 
  #       ids_kable,
  #     "\\hspace{1cm} \\centering ",
  #       ids_kable,
  #     "\\caption{My tables} \\end{table}"))  
  ## ---------------------------------------------------------------------------
  ## -- MOVE TO .RMD -----------------------------------------------------------
  ## ---------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
                  ## START AGAIN FROM HERE ON 06.05.2025 !!!!!!
  
  
  
  
  
  
  
  
  
  ## ----checks, echo = F, collapse = TRUE-------------------------------------------------------------------------
  ##-- Make sure all dead recoveries in DNA are in DR
  check1 <- all(DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"] %in% DR$DNAID) 
  probs_DR_in_DNA <- NULL
  if(!check1){
    tmp <- DNA$DNAID[substr(DNA$RovbaseID,1,1) %in% "M"]
    probs_DR_in_DNA <- tmp[!tmp %in% DR$DNAID]
  } 
  
  tmp <- DNA[substr(DNA$RovbaseID,1,1) %in% "M", ]
  test <- anti_join(tmp,DR, by = names(tmp)[names(tmp) %in% names(DR)])
  
  
  ##-- Make sure that only "Dead recoveries" are in DR 
  check2 <- all(substr(DR$RovbaseID,1,1) %in% "M")
  probs_DNA_in_DR <- NULL
  if(!check2){
    probs_DNA_in_DR <- DR$DNAID[!substr(DR$RovbaseID,1,1) %in% "M"]
  }
  
  
  ## ----merge, echo = F, collapse = TRUE--------------------------------------------------------------------------
  ##-- Merge DNA and dead recoveries files using all shared names columns
  #DATA <- full_join(DNA, DR, by = names(DNA)[names(DNA) %in% names(DR)]) 
  DATA <- merge( DR, DNA,
                 by = c("Id","RovbaseID","DNAID","Species","Sex",
                        "Date","Year","Month","Season",
                        "East_UTM33","North_UTM33",
                        "County","Country_sample"),
                 all = TRUE)
  
  ##-- Determine Death and Birth Years
  DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
  DATA$RovbaseID <- as.character(DATA$RovbaseID)
  DATA$Death <- NA
  DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
  DATA$Birth <- DATA$Death - DATA$Age
  
  ##-- Extract useful numbers
  noID <- sum(is.na(DATA$Id))              ## number of samples without ID
  noDate <- sum(is.na(DATA$Year))          ## number of samples without Date
  noCoords <- sum(is.na(DATA$East_UTM33))  ## number of samples without Coords  
  # notInDR <- sum(as.numeric(substr(DATA$RovbaseID,1,1) %in% "M")
  #                * as.numeric(!(DATA$DNAID %in% DR$DNAID)))
  
  ##-- Filter out unusable samples
  DATA <- DATA %>%
    dplyr::filter(., 
                  ##-- Filter out samples with no ID
                  !is.na(Id),
                  ##-- Filter out samples with no Coordinates
                  !is.na(East_UTM33),
                  ##-- Filter out samples with no dates  
                  !is.na(Year),
                  ##-- Filter out samples with 
                  Year %in% years) %>%
    droplevels(.)
  
  
  ## ----sex assignment, echo = F----------------------------------------------------------------------------------
  ID <- unique(as.character(DATA$Id))
  DATA$Sex <- as.character(DATA$Sex)
  doubleSexID <- IdDoubleSex <- NULL
  counter <- 1
  for(i in 1:length(ID)){
    ##-- Subset data to individual i
    tmp <- DATA$Sex[DATA$Id == ID[i]]
    
    ##-- Number of times individual i was assigned to each sex
    tab <- table(tmp[tmp %in% c("female","male")])
    
    ##-- If conflicting sexes (ID identified as both "female" and "male")
    if(length(tab) == 2){
      ##-- If ID assigned the same number of times to the 2 sexes, assign to Ukjent
      if(tab[1] == tab[2]){
        DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"
      } else {
        ##-- Otherwise pick the most common sex
        DATA$Sex[DATA$Id == ID[i]] <- names(tab)[which(tab == max(tab))]
      }
      # print(paste("Warnings!!!", "Individuals", ID[i], "assigned to both sexes. Now assigned to", names(tab)[which(tab == max(tab))])) 
      IdDoubleSex[counter] <- ID[i]
      counter <- counter + 1
    }
    
    ##-- If only one of "female" or "male" registered
    if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
    
    ##-- If anything else registered : "Ukjent"
    if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
    
    doubleSexID[i] <- length(tab)
  }#i
  
  
  ## ----split DATA, echo = F--------------------------------------------------------------------------------------
  ##-- Split DATA into alive and dead.recovery datasets
  alive <- DATA[is.na(DATA$Death), ]
  dead.recovery <- DATA[!is.na(DATA$Death), ]
  
  ##-- Add earlier detection index
  alive$detected.earlier <-
    unlist(lapply(1:nrow(alive),
                  function(i){
                    this.id <- alive[i,"Id"]
                    this.date <- alive[i,"Date"]
                    any(alive$Id %in% this.id & alive$Date < this.date)
                  }))
  
  dead.recovery$detected.earlier <-
    unlist(lapply(1:nrow(dead.recovery),
                  function(i){
                    this.id <- dead.recovery[i,"Id"]
                    this.date <- dead.recovery[i,"Date"]
                    any(alive$Id %in% this.id & alive$Date < this.date)
                  }))
  
  
  ## ----wolverine, echo = F, collapse = TRUE----------------------------------------------------------------------
  if(engSpecies == "wolverine"){
    ##-- Remove un-verified dead recoveries [HB] 
    ##-- ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
    dead.recovery <- dead.recovery[!grepl(pattern = "Påskutt",
                                          x = as.character(dead.recovery$Outcome)), ]
    
    ##-- Remove suspect NGS samples according to Henrik
    SUSPECT_NGS_SAMPLES <- readMostRecent(
      path = dir.in,
      extension = ".xls",
      pattern = "Remove ngs samples list wolverine")
    alive$DNAID <- as.character(alive$DNAID)
    alive <- alive[!(alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
    
    ##-- Remove suspect dead recoveries according to Henrik
    SUSPECT_DeadRecoSAMPLES <- readMostRecent(
      path = dir.in,
      extension = ".xls",
      pattern = "Remove dead recoveries list wolverine")
    dead.recovery$DNAID <- as.character(dead.recovery$DNAID)
    dead.recovery <- dead.recovery[!(dead.recovery$RovbaseID %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
    
    
    ##-- Remove pups killed before recruitment based on weight (cf. Henrik)
    ##-- 1) remove individuals that are "Ja" in column "Doedt.individ..Unge" and recovered dead between March and November
    youngDeads <- which(dead.recovery$Age_class %in% "Unge" &
                          dead.recovery$Month > 2 &
                          dead.recovery$Month < 12)
    if(length(youngDeads) > 0){
      dead.recovery <- dead.recovery[-youngDeads, ]
    }
    
    
    ##-- 2) remove individuals with 0 <= weight < 4kg between March and November 
    ##-- Format the weight correctly 
    dead.recovery$Weight_total <- as.character(dead.recovery$Weight_total)
    dead.recovery$Weight_slaughter <- as.character(dead.recovery$Weight_slaughter)
    ##-- Convert to decimals
    dead.recovery$Weight_total <- as.numeric(gsub(",", ".", dead.recovery$Weight_total))
    dead.recovery$Weight_slaughter <- as.numeric(gsub(",", ".", dead.recovery$Weight_slaughter))
    ##-- Get the two weight columns together.
    dead.recovery$weight <- ifelse(!is.na(dead.recovery$Weight_total),
                                   dead.recovery$Weight_total,
                                   dead.recovery$Weight_slaughter)
    ##-- Assign negative values to nas to avoid issues
    dead.recovery$weight[is.na(dead.recovery$weight)] <- -999
    
    
    ##-- Check with Henrik (this step does not remove dead recoveries on id with weight==0 should it?)
    ##-- Check how many dead reco we remove and remove if more than 0
    lowWeightDeads <- which(dead.recovery$weight > 0 & dead.recovery$weight < 4 &
                              dead.recovery$Month > 2 & dead.recovery$Month < 12)
    if(length(lowWeightDeads) > 0){
      dead.recovery <- dead.recovery[-lowWeightDeads, ]
    }
    
    ##-- Check how many dead reco with a weight of 0 kg and recovered between March and November
    zeroWeightDeads <- which(dead.recovery$Age %in% 0 &
                               dead.recovery$Month > 2 &
                               dead.recovery$Month < 12)
  }
  
  
  ## ----wolf, echo = F, collapse = TRUE---------------------------------------------------------------------------
  if(engSpecies == "wolf"){
    ##-- Load most recent Micke's file
    INDIVIDUAL_ID <- readMostRecent.csv(
      path = dir.in,
      pattern = "_ID Grouping ",
      fileEncoding = "latin1")  
    
    ##-- Translate Scandinavian characters
    colnames(INDIVIDUAL_ID) <- translateForeignCharacters( data = colnames(INDIVIDUAL_ID))
    
    ##-- Overwrite gender from Micke's data when available
    micke.sex <- as.character(unlist(lapply(DATA$Id,
                                            function(i){ 
                                              INDIVIDUAL_ID[as.character(INDIVIDUAL_ID$Individ..Rovbase.) %in% i,"Sex"][1]
                                            })))
    micke.sex[micke.sex %in% "0"] <- NA
    micke.sex[micke.sex %in% names(table(micke.sex))[3]] <- NA
    micke.sex[micke.sex %in% "Hona"] <- "female"
    micke.sex[micke.sex %in% "Hane"] <- "male"
    new.sex <- ifelse(!is.na(micke.sex), as.character(micke.sex), as.character(DATA$Sex))
    DATA$Sex <- new.sex
    
    numOverwiteSex <- sum(unique(as.character(INDIVIDUAL_ID$Individ..Rovbase.)) %in% DATA$Id)
  }
  
  
  ## ----bear, echo = F, collapse = TRUE---------------------------------------------------------------------------
  if(engSpecies == "bear"){
    ##-- Load most recent "flagged" file from HB
    flagged <- readMostRecent( 
      path = dir.in,
      extension = ".csv",
      pattern = "dna_bear_to_remove", 
      fileEncoding = "Latin1") 
    
    ##-- Remove flagged samples 
    remove.alive <- !alive$Barcode_sample %in% flagged$Strekkode
    alive <- alive[remove.alive, ]
    remove.dead <- !dead.recovery$Barcode_sample %in% flagged$Strekkode
    dead.recovery <- dead.recovery[remove.dead, ]
    dead.recovery$Missing <- NA
    dead.recovery$Individ <- NA
  }
  
  
  ## ----turn into sf, echo = F, collapse = TRUE-------------------------------------------------------------------
  ##-- Turn into sf points dataframe
  alive <- sf::st_as_sf( x = alive,
                         coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633)) 
  
  ##-- Intersect and extract country name
  #alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
  alive$Country_sf[!is.na(as.numeric(st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  alive$Country_sf[!is.na(as.numeric(st_intersects(alive, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  ##-- Turn into sf points dataframe
  dead.recovery <- sf::st_as_sf( x = dead.recovery,
                                 coords = c("East_UTM33","North_UTM33")) %>%
    sf::st_set_crs(.,sf::st_crs(32633))
  
  ##-- Intersect and extract country name
  #dead.recovery$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(dead.recovery, COUNTRIES))]
  dead.recovery$Country_sf[!is.na(as.numeric(st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
  dead.recovery$Country_sf[!is.na(as.numeric(st_intersects(dead.recovery, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"
  
  
  
  ## ---- Rovbase data summary ---------------------------------------------------
  ##-- Number of NGS samples
  samples <- table(alive$Country_sample, alive$Year)
  samples2 <- table(alive$Country_sf, alive$Year)
  samples <- rbind(samples, "Total" = colSums(samples))
  samples <- cbind(samples, "Total" = rowSums(samples))
  
  
  ##-- Number of individuals detected alive
  ids <- apply(table(alive$Country_sample,
                     alive$Year,
                     alive$Id),
               c(1,2),
               function(x)sum(x>0))
  
  ids <- rbind(ids,
               "Total" = apply(table(alive$Year,
                                     alive$Id),
                               1,
                               function(x)sum(x>0)))
  
  ids <- cbind(ids,
               "Total" = c(apply(table(alive$Country_sample,alive$Id),
                                 1,
                                 function(x)sum(x>0)),
                           length(unique(alive$Id))))
  
  
  ##-- Number of DR samples
  deadSamples <- table(dead.recovery$Country_sample, dead.recovery$Year)
  deadSamples <- rbind(deadSamples, "Total" = colSums(deadSamples))
  deadSamples <- cbind(deadSamples, "Total" = rowSums(deadSamples))
  
  
  ##-- Number of individuals recovered
  deadIds <- apply(table(dead.recovery$Country_sample,
                         dead.recovery$Year,
                         dead.recovery$Id),
                   c(1,2),
                   function(x)sum(x>0))
  
  deadIds <- rbind(deadIds,
                   "Total" = apply(table(dead.recovery$Year,
                                         dead.recovery$Id),
                                   1,
                                   function(x)sum(x>0)))
  
  deadIds <- cbind(deadIds,
                   "Total" = c(apply(table(dead.recovery$Country_sample,
                                           dead.recovery$Id),
                                     1,
                                     function(x)sum(x>0)),
                               length(unique(dead.recovery$Id))))
  
  
  
  
  ## ---- num samples - FIGURE ---------------------------------------------------
  
  ##-- Number of NGS per month
  dat.alive <- alive %>%
    mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = n())
  dat.alive$type = "NGS"
  
  ##-- Number of dead recoveries per month
  dat.dead <- dead.recovery %>%
    mutate(Date = trunc(Date, "month")) %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(n = -n())
  dat.dead$type = "dead.recovery"
  
  ##-- Combine NGS and dead recoveries
  dat <- rbind(dat.alive, dat.dead)
  dat$Date <- as.Date(dat$Date)
  
  ##-- Plot time series of number of samples per month
  samplesTimeSeries <- ggplot(dat) +
    geom_col(aes(x = Date, y = n, fill = type)) +
    ylab("Number of samples") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position.inside = c(0.1,0.9),
          axis.text.x=element_text(angle = 60,
                                   hjust = 1)) +
    scale_x_date( date_breaks = "years",
                  date_labels = "%Y") 
  
  ##-- Save plot as .png
  grDevices::png( filename = file.path(dir.out, "figures",
                                       paste0(species, "_monitoring_", years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  samplesTimeSeries
  graphics.off()
  
  
  
  ## ---- num samples - TABLES ---------------------------------------------------
  ##-- Number of NGS samples
  kable( samples,
         align = "lc",
         caption = "Number of NGS samples per year and country") %>%
    kable_styling(full_width = F)
  
  ##-- Number of dead recoveries
  kable( deadSamples,
         align = "lc",
         caption = "Number of DNA samples from dead animals per year and country") %>% 
    kable_styling(full_width = F)
  
  
  
  ## ---- num ids ---------------------------------
  ##-- Number of IDs
  dat.alive <- alive %>% 
    dplyr::group_by(Year) %>% 
    dplyr::summarise(n = length(unique(Id)))
  dat.alive$type = "NGS"
  
  ##-- Number of Dead Recoveries
  dat.dead <- dead.recovery %>% 
    dplyr::group_by(Year) %>% 
    dplyr::summarise(n = -length(unique(Id)))
  dat.dead$type = "dead.recovery"
  
  ##-- Combine NGS and dead recoveries
  dat <- rbind(dat.alive, dat.dead)
  
  ##-- Plot time series of number of IDs per year
  idsTimeSeries <- ggplot(dat) +
    geom_col(aes(x = Year, y = n, fill = type)) +
    ylab("Number of individuals") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.title = element_blank(),
          legend.position = c(0.1,0.9)) +
    scale_x_continuous( breaks = years,
                        labels = years)
  
  ##-- Save plot as .png
  grDevices::png( filename = file.path(dir.out, "figures",
                                       paste0(species, "_monitoring_", years[1]," to ", years[length(years)], ".png")),
                  width = 8, height = 6,
                  units = "in", res = 300)
  samplesTimeSeries
  graphics.off()
  
  
  ## ---- Kables ----------------------------------------------------------------
  ##-- Number of ID detected alive
  kable(ids, align = "lc",
        caption = "Number of individuals detected through NGS per year and country") %>%
    kable_styling(full_width = F)
  
  ##-- Number of ID detected alive
  kable(deadIds, align = "lc",
        caption = "Number of identified dead animals per year and country") %>% 
    kable_styling(full_width = F)
  
  
  
  ## ---- previous det ----------------------------
  dat.alive <- alive %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id)))
  dat.alive$type = "NGS"
  
  dat.dead <- dead.recovery %>%
    dplyr::group_by(Year, detected.earlier) %>%
    dplyr::summarise(n = length(unique(Id)))
  dat.dead$type = "dead.recovery"
  
  ggplot() +
    geom_col(data = dat.alive,
             aes(x = Year, y = n, fill = detected.earlier)) +
    ylab("Number of individuals detected through NGS") +
    theme(legend.position = c(0.1,0.9)) +
    scale_x_continuous(breaks = years, labels = years)
  
  ggplot() +
    geom_col(data = dat.dead,
             aes(x = Year, y = n, fill = detected.earlier)) +
    ylab("Number of dead individuals") +
    theme(legend.position = c(0.1,0.9)) +
    scale_x_continuous(breaks = years, labels = years)
  
  
  
  ## ---- multiple sex -----------------------------------------------------------
  sexTab <- cbind.data.frame(
    "problems" = c("Unknown sex", "both 'female' and 'male'"),
    "number of individuals" = as.numeric(table(doubleSexID)[c(1,3)]))
  
  kable(sexTab, align = "lc") %>%
    kable_styling(full_width = F)
  
  
  
  ## ---- multiple deaths --------------------------------------------------------
  ##-- Identify and count individuals dead "more than once"
  ID <- names(table(dead.recovery$Id))[table(dead.recovery$Id)>1]
  multiDeathDate <- multiDeathYear <- multiDeathLocs <-  NULL
  for(i in 1:length(ID)){
    tmp <- dead.recovery[dead.recovery$Id == ID[i], ] 
    ##-- Multiple death dates
    if(length(unique(tmp$Date)) > 1){
      multiDeathDate <- c(multiDeathDate, ID[i])
    }
    ##-- Multiple death years
    if(length(unique(tmp$Year)) > 1){
      multiDeathYear <- c(multiDeathYear, ID[i])
    }
    ##-- Multiple death locations
    if(length(unique(tmp$East)) > 1 | length(unique(tmp$North)) > 1){
      multiDeathLocs <- c(multiDeathLocs, ID[i])
    }
  }#i
  
  ##-- Remove individuals that died more than once
  dead.recovery$Id <- as.character(dead.recovery$Id)
  IdDoubleDead <- names(table(dead.recovery$Id))[table(dead.recovery$Id) > 1]
  if(length(IdDoubleDead) > 0){
    for(i in IdDoubleDead){
      ##-- Identify repeated deaths
      tmp <- which(dead.recovery$Id %in% i) 
      ##-- Try to keep death with known death cause 
      tmp2 <- which(!is.na(dead.recovery$DeathCause_2[tmp]))[1]
      if(length(tmp2) == 0){tmp <- tmp[-1]} else {tmp <- tmp[!tmp %in% tmp2]}
      ##-- Remove repeated deaths
      dead.recovery <- dead.recovery[-tmp, ]
    }#i
  }#if
  
  
  
  ## ---- ghost individuals ------------------------------------------------------
  id.list <- unique(c(as.character(dead.recovery$Id), as.character(alive$Id)))
  ghosts <- unlist(lapply(id.list, function(id){
    out <- NULL
    try({
      if(id %in% dead.recovery$Id){
        mort.year <- min(dead.recovery$Year[dead.recovery$Id == id])
        this.alive <- alive[alive$Id == id, ]
        ## Was it detected alive in any season after death?
        temp <- this.alive[this.alive$Year > mort.year, ]
        if(length(temp) > 0){
          out <- rownames(temp)
          names(out) <- id
        }## FLAG THOSE FOR HENDRIK
      }
    }, silent = TRUE)
    return(out)
  }))
  samples.to.remove <- unlist(ghosts)
  
  ##-- Remove flagged NGS detections after dead recovery
  alive <- alive[!rownames(alive) %in% samples.to.remove, ]
  
  
  
  ## ---- NGS maps ----------------------------------
  ##-- NGS map
  numRows <- ceiling(length(years)/5)
  numCols <- 5
  NGS_map <- ggplot(data = alive) +
    geom_sf(data = COUNTRIES,
            aes(fill = ISO),
            alpha = 0.3,
            color = NA) +
    geom_sf(color = "black",
            alpha = 0.3, size = 0.8, pch = 3) +
    facet_wrap(~Year, nrow = numRows, ncol = numCols) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  
  ##-- Save maps as .png
  grDevices::png(filename = file.path(dir.out, "figures",
                                      paste0(species, "_NGS_", years[1]," to ", years[length(years)], ".png")),
                 width = 8, height = 6, units = "in", res = 300)
  NGS_map
  graphics.off()
  
  ## ---- Dead recovery maps ------------------------
  dead_map <- ggplot(data = dead.recovery) +
    geom_sf(data = COUNTRIES, 
            aes(fill = ISO),
            alpha = 0.3,
            color = NA) + 
    geom_sf(color = "black", alpha = 0.5, size = 0.8, pch = 3) +
    facet_wrap(~Year, nrow = numRows, ncol = numCols) +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  
  ##-- Save maps as .png
  grDevices::png(filename = file.path(dir.out, "figures",
                                      paste0(engSpecies, "_DEAD_", years[1]," to ", years[length(years)], ".png")),
                 width = 8, height = 6, units = "in", res = 300)
  dead_map
  graphics.off()
  
  
  ## ---- save data ---------------------------------------------------------------------------------------
  fileName <-  paste0("CleanData_", engSpecies, "_",DATE,".RData")
  
  save( alive, 
        dead.recovery,
        IdDoubleSex,
        file = file.path(dir.out, "data", fileName))
  
  
  
  
}