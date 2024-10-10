##------------------------------------------------------------------------------
##
## Script name: RovQuant BEAR OPSCR analysis 
##
## Purpose of script: 
## This R script performs:
## 1. the initial cleaning of the large carnivore NGS data downloaded from RovBase.3.0
## 2. the data preparation for analysis with the 'nimbleSCR' package
## 3. the model fitting using 'nimble' and 'nimbleSCR'
## 4. the post-processing of the MCMC output
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: `r paste(Sys.Date())`
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), `r paste(format(Sys.Date(), "%Y"))`
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
##------------------------------------------------------------------------------
##
## Notes: 
##   
##
##------------------------------------------------------------------------------
rm(list = ls())
gc()


## ------ IMPORT REQUIRED LIBRARIES ------

library(devtools)
library(dplyr)
library(readxl)
library(nimbleSCR)
library(rovquantR)



##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data_dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/Data"


##-- WORKING DIRECTORY (= main folder for the analysis)
working_dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/test2"


##-- Create folder structure for the analysis
makeDirectories(path = working_dir)



##------------------------------------------------------------------------------
## ----- II. LOAD AND FIX RAW ROVBASE.3.0 DATA -----

## -----    1. NGS SAMPLES ------

##-- Load raw excel file imported from rovbase 
DNA <- readxl::read_xlsx(file.path(data_dir,"DNA.xlsx")) %>%
  ##-- Remove any duplicates
  distinct(., .keep_all = TRUE) %>%
  ##-- Rename columns to facilitate manipulation
  rename(., any_of(c( DNAID_sample = "DNAID (Prøve)",
                      Barcode_sample = "Strekkode (Prøve)",
                      RovbaseID_sample = "RovbaseID (Prøve)",
                      EventID = "HendelseID",
                      Analysis_priority = "Analyseprioritet",
                      Species_sample = "Art (Prøve)",
                      Sample_type = "Prøvetype", 
                      Date = "Funnetdato",
                      Mountain_area = "Fjellområde",
                      Sensitivity = "Følsomhet",
                      Release_Date = "Frigivelsesdato", 
                      Locality = "Lokalitet",
                      Location = "Funnsted",
                      North_original = "Nord (opprinnelig)",
                      East_Original = "Øst (opprinnelig)",
                      North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
                      East_UTM33 = "Øst (UTM33/SWEREF99 TM)" ,
                      North_RT90 = "Nord (RT90)",
                      East_RT90 = "Øst (RT90)" ,
                      Coordinate_system = "Koordinatsystem" ,
                      Site_quality = "Stedkvalitet",
                      Collected_by = "Hvem samlet inn",
                      Collector_name = "Samlet selv - Navn",
                      Collector_phone = "Samlet selv - Telefon", 
                      Collector_email = "Samlet selv - E-post",
                      Collector_role = "Samlet selv - Rolle",
                      Collector_other_name = "Annen innsamler - Navn" ,
                      Collector_other_phone = "Annen innsamler - Telefon", 
                      Collector_other_email = "Annen innsamler - E-post",
                      Collector_other_role = "Annen innsamler - Rolle",
                      Tips_name = "Tipser - Navn",
                      Tips_phone = "Tipser - Telefon",
                      Tips_email = "Tipser - E-post",
                      Tips_role = "Tipser - Rolle",
                      Quality_checked = "Kvalitetssikret av feltpersonell",
                      Quality_check_name = "Kvalitetssikrer - navn",
                      Quality_check_orga = "Kvalitetssikrer - Organisasjon",
                      Comments_sample = "Merknad (Prøve)",
                      Last_saved_by_sample = "Sist lagret av (Prøve)",
                      Last_saved_sample = "Sist lagret dato (Prøve)",
                      DNAID = "DNAID (Analyse)",
                      Barcode = "Strekkode (Analyse)",
                      RovbaseID = "RovbaseID (Analyse)",
                      Species = "Art (Analyse)",
                      Sample_status = "Prøvestatus",
                      Id = "Individ",
                      Sex_analysis = "Kjønn (Analyse)",
                      Sex = "Kjønn (Individ)",
                      Method = "Metode",
                      Analyzed_by = "AnalysertAv",
                      Comments = "Merknad (Analyse)",
                      Last_saved_by = "Sist lagret av (Analyse)" ,
                      Last_saved = "Sist lagret dato (Analyse)" ,
                      Municipality_number = "Kommunenummer",
                      Municipality = "Kommune", 
                      County_number = "Fylkenummer",
                      County = "Fylke"))) %>%
  ##-- Filter to the focal species
  filter(., Species == "Bjørn") %>%
  ##-- Save as.csv 
  writeMostRecent.csv( ., file = file.path(data_dir, "dna_bear.csv"))



## -----    2. DEAD RECOVERY DATA ------

##-- Load raw excel file imported from rovbase 
DR <- readxl::read_xlsx(file.path(data_dir,"dead carnivores.xlsx")) %>%
  ##-- Remove any duplicates
  distinct(., .keep_all = TRUE) %>%
  ##-- Rename columns to facilitate manipulation
  rename(., any_of(c( Species = "Art",
                      Death_cause = "Bakgrunn/årsak",
                      Death_method = "Bakgrunn/årsak metode",
                      Death_purpose = "Bakgrunn/årsak formål", 
                      Date = "Dødsdato",
                      Uncertain_date = "Usikker dødsdato",
                      Time_of_death = "Dødstidspunkt",
                      Age_class = "Alder på dødt individ",
                      Age_class_verif = "Aldersklasse verifisert SVA",
                      Juvenile = "Yngling",
                      Sex = "Kjønn", 
                      Age_estimated = "Alder, vurdert",
                      Age = "Alder, verifisert",
                      Age_verif_by = "Alder, verifisert av",
                      Full_weight =  "Helvekt",
                      Slaughter_weight =  "Slaktevekt",
                      CITES = "CITES-nummer",
                      Felling_site_verif = "Kontroll av fellingsted",
                      Tissue_sample = "Vevsprøve tatt",
                      Hunting_date = "Observasjons/Jaktdato",
                      Field_personnel ="Feltpersonell",
                      Control_status = "Kontrollstatus",
                      Assessment = "Vurdering",
                      Location = "Funnsted",
                      North_original = "Nord (opprinnelig)",
                      East_Original = "Øst (opprinnelig)",
                      North_UTM33 = "Nord (UTM33/SWEREF99 TM)",
                      East_UTM33 = "Øst (UTM33/SWEREF99 TM)",
                      North_RT90 = "Nord (RT90)",
                      East_RT90 = "Øst (RT90)",
                      Coordinate_system = "Koordinatsystem",
                      Site_quality = "Stedkvalitet",
                      Approved_date = "Godkjentdato",
                      Outcome = "Utfall", 
                      Counted_off_against_decision = "Regnes av mot vedtak", 
                      Approved_by = "Godkjent av",
                      Sensitivity = "Følsomhet",
                      Release_Date = "Frigivelsesdato", 
                      Id = "Individ", 
                      SVAID = "SVAID",
                      Lansstyrelsen_number = "Länsstyrelsens nr",
                      Last_saved_by = "Sist lagret av",
                      Last_saved =  "Sist lagret dato",
                      Municipality_number =  "Kommunenummer",
                      Municipality = "Kommune", 
                      County_number = "Fylkenummer",
                      County = "Fylke"))) %>%
  ##-- Filter to the focal species
  filter(., Species == "Bjørn") %>%
  ##-- Save as.csv 
  writeMostRecent.csv( ., file = file.path(data_dir, "dead_bear.csv"))



##------------------------------------------------------------------------------
## ----- III. CLEAN NGS DATA -----

cleanRovbaseData( species = "bear",
                  years = 2020:2024,
                  data_dir = data_dir,
                  output_dir = file.path(working_dir, "data"))



##------------------------------------------------------------------------------
## ----- IV. PREPARE OPSCR DATA ------

makeRovquantData(    
  ##-- paths
  species = "bear",
  data_dir = data_dir,
  working_dir = working_dir,
  ##-- data
  #years = 2020:2023,
  sex = c("Hann","Hunn"),
  aug.factor = 2,
  sampling.months = list(4,5,6,7,8,9,10,11),
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 50000,
  max.move.dist = 250000,
  ##-- detectors
  detector.res = 5000,
  subdetector.res = 1000,
  max.det.dist = 70000,
  resize.factor = 1,
  ##-- miscellanious
  print.report = F)



##------------------------------------------------------------------------------
## ----- V. FIT ROVQUANT MODELS ------
## -----   1. Females ------
##-- List all prepared input files
inputFiles <- list.files(file.path(working_dir, "nimbleInFiles/Hunn"),
                         full.names = T)

##-- Load the first one
load(inputFiles[1]) 

##-- Build nimble model object
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 
model$calculate()
cmodel <- compileNimble(model)
conf <- configureMCMC( model,
                       monitors = nimParams,
                       thin = 1,
                       monitors2 = nimParams2,
                       thin2 = 5)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble( list(model = model,
                                    mcmc = Rmcmc),
                               showCompilerOutput = F)
Cmcmc <- compiledList$mcmc

##-- RUN NIMBLE MCMC IN SUCCESSIVE BITES
system.time(runMCMCbites( mcmc = Cmcmc,
                          bite.size = 100,
                          bite.number = 5,
                          path = file.path(working_dir,"nimbleOutfiles/Hunn")))



## -----   2. Males ------
##-- List all prepared input files
inputFiles <- list.files(file.path(working_dir, "nimbleInFiles/Hann"),
                         full.names = T)

##-- Load the first one
load(inputFiles[1]) 

##-- Build nimble model object
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      inits = nimInits,
                      data = nimData,
                      check = FALSE,
                      calculate = FALSE) 
model$calculate()
cmodel <- compileNimble(model)
conf <- configureMCMC(model,
                      monitors = nimParams,
                      thin = 1,
                      monitors2 = nimParams2,
                      thin2 = 5)
Rmcmc <- buildMCMC(conf)
compiledList <- compileNimble(list(model = model,
                                   mcmc = Rmcmc),
                              showCompilerOutput = F)
Cmcmc <- compiledList$mcmc

##-- RUN NIMBLE MCMC IN SUCCESSIVE BITES
system.time(runMCMCbites( mcmc = Cmcmc,
                          bite.size = 100,
                          bite.number = 5,
                          path = file.path(working_dir,"nimbleOutfiles/Hann")))



##------------------------------------------------------------------------------
## ----- VI. PROCESS ROVQUANT OUTPUT ------
processRovquantOutput(   
  ##-- paths
  species = "bear",
  data_dir = data_dir,
  working_dir = working_dir,
  ##-- MCMC processing
  nburnin = 0,
  ##-- Density extraction
  niter = 100,
  extraction.res = 5000,
  ##-- miscellanious
  plot.check = FALSE,
  print.report = TRUE)



##------------------------------------------------------------------------------