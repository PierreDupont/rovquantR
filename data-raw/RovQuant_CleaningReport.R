params <- list(
  species = "Jerv",
  years = 2014:2023,
  sex = c("Hunn", "Hann"),
  sampling.months = list(12,1:6),
  dir.in = "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/rovquantR/wolverine/2024",
  dir.out = "C:/My_documents/rovquantR",
  date = structure(20115, class = "Date"))

##-- load libraries
library(kableExtra)
library(ggplot2)

##-- Rendering parameters
species <- params$species
years <- params$years
sex <- params$sex
sampling.months <- params$sampling.months
rename.list <- params$rename.list
dir.in <- params$dir.in
dir.out <- params$dir.out

##-- Species
if(sum(grep("bear", species, ignore.case = T))>0|sum(grep("bjørn", species, ignore.case = T))>0|sum(grep("bjorn", species, ignore.case = T))>0){
  engSpecies <- "bear"
  norSpecies <- c("Bjørn", "BjÃ¸rn")
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

##-- Sampling months
if(is.null(sampling.months)){
  if(engSpecies == "bear"){
    sampling.months <- list(4:11)
  } else {
    if(engSpecies == "wolf"){
      sampling.months <- list(c(10:12),c(1:4))
    } else {
      if(engSpecies == "wolverine"){
        sampling.months <- list(c(12),c(1:6))
      } else {
        stop("No default setting available for the monitoring period of this species. \n You must specify the monitoring season months through the 'sampling.months' argument.")
      }
    }
  }
}

##-- Renaming list
if(is.null(rename.list)){
rename.list = c(
  Age_estimated = "Alder, vurdert",
  Age = "Alder, verifisert",
  Age_verif_by = "Alder, verifisert av",
  Age_class = "Alder på dødt individ",
  Age_class_verif = "Aldersklasse verifisert SVA",
  Analyzed_by = "AnalysertAv",
  Analysis_priority = "Analyseprioritet",
  Approved_by = "Godkjent av",
  Approved_date = "Godkjentdato",
  Assessment = "Vurdering",
  Barcode_sample = "Strekkode (Prøve)",
  Barcode = "Strekkode (Analyse)",
  Birth_territory = "Født revir",
  CITES = "CITES-nummer",
  Collected_by = "Hvem samlet inn",
  Collector_name = "Samlet selv - Navn",
  Collector_phone = "Samlet selv - Telefon",
  Collector_email = "Samlet selv - E-post",
  Collector_role = "Samlet selv - Rolle",
  Collector_other_name = "Annen innsamler - Navn" ,
  Collector_other_phone = "Annen innsamler - Telefon",
  Collector_other_email = "Annen innsamler - E-post",
  Collector_other_role = "Annen innsamler - Rolle",
  Comments_sample = "Merknad (Prøve)",
  Comments = "Merknad (Analyse)",
  Control_status = "Kontrollstatus",
  Coordinate_system = "Koordinatsystem",
  Counted_off_against_decision = "Regnes av mot vedtak",
  County_number = "Fylkenummer",
  County = "Fylke",
  Date = "Funnetdato",
  Date = "Dødsdato",
  Death_cause = "Bakgrunn/årsak",
  Death_method = "Bakgrunn/årsak metode",
  Death_purpose = "Bakgrunn/årsak formål",
  DNAID_sample = "DNAID (Prøve)",
  DNAID = "DNAID (Analyse)",
  EventID = "HendelseID",
  East_Original = "Øst (opprinnelig)",
  East_RT90 = "Øst (RT90)",
  East_UTM33 = "Øst (UTM33/SWEREF99 TM)",
  Felling_site_verif = "Kontroll av fellingsted",
  Field_personnel ="Feltpersonell",
  Hunting_date = "Observasjons/Jaktdato",
  Id = "Individ",
  Juvenile = "Yngling",
  Mountain_area = "Fjellområde",
  Method = "Metode",
  Municipality_number = "Kommunenummer",
  Municipality = "Kommune",
  North_original = "Nord (opprinnelig)",
  North_RT90 = "Nord (RT90)",
  North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
  Origin = "Opprinnelse",
  Outcome = "Utfall",
  Last_saved_by_sample = "Sist lagret av (Prøve)",
  Last_saved_sample = "Sist lagret dato (Prøve)",
  Last_saved_by = "Sist lagret av (Analyse)",
  Last_saved = "Sist lagret dato (Analyse)",
  Last_saved_by = "Sist lagret av",
  Last_saved =  "Sist lagret dato",
  Locality = "Lokalitet",
  Location = "Funnsted",
  Lansstyrelsen_number = "Länsstyrelsens nr",
  Quality_checked = "Kvalitetssikret av feltpersonell",
  Quality_check_name = "Kvalitetssikrer - navn",
  Quality_check_orga = "Kvalitetssikrer - Organisasjon",
  Release_Date = "Frigivelsesdato",
  Sample_type = "Prøvetype",
  Sensitivity = "Følsomhet",
  Species_sample = "Art (Prøve)",
  Site_quality = "Stedkvalitet",
  Time_of_death = "Dødstidspunkt",
  Tips_name = "Tipser - Navn",
  Tips_phone = "Tipser - Telefon",
  Tips_email = "Tipser - E-post",
  Tips_role = "Tipser - Rolle",
  Tissue_sample = "Vevsprøve tatt",
  Release_Date = "Frigivelsesdato",
  RovbaseID = "RovbaseID (Analyse)",
  RovbaseID_sample = "RovbaseID (Prøve)",
  Species = "Art (Analyse)",
  Species = "Art",
  Sample_status = "Prøvestatus",
  Sensitivity = "Følsomhet",
  Sex_analysis = "Kjønn (Analyse)",
  Sex = "Kjønn (Individ)",
  Sex = "Kjønn",
  Site_quality = "Stedkvalitet",
  SVAID = "SVAID",
  Uncertain_date = "Usikker dødsdato",
  Weight_slaughter = "Slaktevekt",
  Weight_total =  "Helvekt")
}


## ------   1. LOAD HABITAT DATA ------
##-- Load pre-processed habitat shapefiles
data(COUNTRIES, envir = environment()) 



## ------   2. NGS DATA ------
## ------     2.1. LOAD NGS ------
## ------     2.2. RENAME ------
## ------   7. FORMAT DATES -----
DNA <- suppressWarnings(readMostRecent( path = dir.in,
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
    MonSeason = ifelse( Month < unlist(sampling.months)[1],
                        Year,
                        Year-1),
    ##-- Fix unknown "Id"
    Id = ifelse(Id %in% "", NA, Id),
    ##-- Fix unknown "Sex"
    Sex = ifelse(Sex %in% "Ukjent", NA, Sex)) %>%
  ##-- Filter to the focal years
  filter(., Year %in% years) 



## ------   3. DEAD RECOVERIES ------
## ------     3.1. LOAD DEAD RECOVERIES ------
## ------     3.2. RENAME ------
## ------   7. FORMAT DATES -----
##-- Load raw excel file imported from rovbase 
DR <- suppressWarnings(readMostRecent( path = dir.in,
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
    Sex = ifelse(Sex %in% "Ukjent", NA, Sex)) %>%
  ##-- Filter to the focal years
  filter(., Year %in% years) 



## ------   6. MERGE NGS & DEAD RECOVERIES ------
##-- Merge DNA and dead recoveries files using all shared names columns
#DATA <- full_join(DNA, DR, by = names(DNA)[names(DNA) %in% names(DR)]) 
DATA <- merge( DR, DNA,
               by = c("Id","RovbaseID","DNAID","Species","Sex",
                      "Date","Year","Month","Season",
                      "East_UTM33","North_UTM33",
                      "County","Country_sample"),
               all = TRUE)



## ------     8.2. DETERMINE AGE ------
##-- Determine Death and Birth Years
DATA$Age <- suppressWarnings(as.numeric(as.character(DATA$Age))) 
DATA$RovbaseID <- as.character(DATA$RovbaseID)
DATA$Death <- NA
DATA$Death[substr(DATA$RovbaseID,1,1) %in% "M"] <- DATA$Year[substr(DATA$RovbaseID,1,1) %in% "M"]
DATA$Birth <- DATA$Death - DATA$Age



## ------     8.1. REMOVE UNUSABLE SAMPLES ------
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



## ------     8.4. CHECK SEX ASSIGNMENT ------
##-- Check sex assignment
ID <- unique(as.character(DATA$Id))
DATA$Sex <- as.character(DATA$Sex)
doubleSexID <- IdDoubleSex <- NULL
counter <- 1
for(i in 1:length(ID)){
  ##-- Subset data to individual i
  tmp <- DATA$Sex[DATA$Id == ID[i]]
  
  ##-- Number of times individual i was assigned to each sex
  tab <- table(tmp[tmp %in% c("Hunn","Hann")])
  
  ##-- If conflicting sexes (ID identified as both "Hunn" and "Hann")
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
  
  ##-- If only one of "Hunn" or "Hann" registered
  if(length(tab) == 1){DATA$Sex[DATA$Id == ID[i]] <- names(tab)}
  
  ##-- If anything else registered : "Ukjent"
  if(length(tab) == 0){DATA$Sex[DATA$Id == ID[i]] <- "Ukjent"}
  
  doubleSexID[i] <- length(tab)
}#i



## ------     8.3. EXTRACT COUNTRY ------
##-- Turn into sf points dataframe
DATA <- sf::st_as_sf( x = DATA,
                       coords = c("East_UTM33","North_UTM33")) %>%
  sf::st_set_crs(.,sf::st_crs(32633)) 

##-- Intersect and extract country name
#alive$Country_sf <- COUNTRIES$ISO[as.numeric(sf::st_intersects(alive, COUNTRIES))]
DATA$Country_sf[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[COUNTRIES$ISO %in% "NOR", ])))] <- "(N)"
DATA$Country_sf[!is.na(as.numeric(st_intersects(DATA, COUNTRIES[COUNTRIES$ISO %in% "SWE", ])))] <- "(S)"





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

if(engSpecies == "wolverine"){
  ## ------     3.3. REMOVE UNIVERIFIED DR ------
  ##-- Remove un-verified dead recoveries [HB] 
  ##-- ("Påskutt ikke belastet kvote" & "Påskutt belastet kvote")
  dead.recovery <- dead.recovery[!grepl(pattern = "Påskutt",
                                        x = as.character(dead.recovery$Outcome)), ]
  
  
  ## ------     8.6. REMOVE SUSPECT NGS ACCORDING TO HENRIK ------
  ##-- Remove suspect NGS samples according to Henrik
  SUSPECT_NGS_SAMPLES <- readMostRecent(
    path = dir.in,
    extension = ".xls",
    pattern = "Remove ngs samples list wolverine")
  alive$DNAID <- as.character(alive$DNAID)
  alive <- alive[!(alive$DNAID %in% as.character(SUSPECT_NGS_SAMPLES$DNAID_RB)), ]
  
  
  ## ------     8.7. REMOVE SUSPECT DEAD RECOVERIES ACCORDING TO HENRIK ------
  ##-- Remove suspect dead recoveries according to Henrik
  SUSPECT_DeadRecoSAMPLES <- readMostRecent(
    path = dir.in,
    extension = ".xls",
    pattern = "Remove dead recoveries list wolverine")
  dead.recovery$DNAID <- as.character(dead.recovery$DNAID)
  dead.recovery <- dead.recovery[!(dead.recovery$RovbaseID %in% as.character(SUSPECT_DeadRecoSAMPLES$Rovbase_ID)), ]
  
  
  ## ------     8.9. Remove pups killed before recruitment based on weight (cf. Henrik) ------
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



## ------     8.8. Remove individuals that died twice ------

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



##-- Identify individuals detected after their supposed death
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




