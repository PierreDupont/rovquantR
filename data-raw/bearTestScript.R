##------------------------------------------------------------------------------
##
## Script name: RovQuant BEAR OPSCR analysis 
##
## Purpose of script: 
## This R script performs:
## 1. the initial cleaning of the Brown bear NGS data downloaded from RovBase.3.0
## 2. the data preparation for the RovQuant OPSCR analysis with the 'nimbleSCR' package
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

# devtools::install_github("PierreDupont/rovquantR")
## Ctrl + Shift + F10 (to restart R session)
library(rovquantR)
library(nimbleSCR)



##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/test.2"



##------------------------------------------------------------------------------
## ----- II. CLEAN NGS DATA -----

cleanRovbaseData_bear( 
  species = "bear",
  years = 2020:2024,
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  ##-- paths
  species = "bear",
  data.dir = data.dir,
  working.dir = working.dir,
  ##-- data
  sex = c("Hann","Hunn"),
  aug.factor = 2,
  ##-- habitat
  habitat.res = 20000, 
  buffer.size = 50000,
  max.move.dist = 250000,
  ##-- detectors
  detector.res = 5000,
  subdetector.res = 1000,
  max.det.dist = 70000)



##------------------------------------------------------------------------------
## ----- IV. FIT ROVQUANT MODELS ------

## -----   1. Females ------

##-- List all prepared input files
inputFiles <- list.files(file.path(working.dir, "nimbleInFiles/Hunn"),
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
                          path = file.path(working.dir,"nimbleOutfiles/Hunn")))



## -----   2. Males ------

##-- List all prepared input files
inputFiles <- list.files(file.path(working.dir, "nimbleInFiles/Hann"),
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
                          path = file.path(working.dir,"nimbleOutfiles/Hann")))



##------------------------------------------------------------------------------
## ----- VI. PROCESS ROVQUANT OUTPUT ------

processRovquantOutput(   
  ##-- paths
  species = "bear",
  data.dir = data.dir,
  working.dir = working.dir,
  ##-- MCMC processing
  nburnin = 0,
  ##-- Density extraction
  niter = 100,
  extraction.res = 5000,
  ##-- miscellanious
  print.report = TRUE)


##------------------------------------------------------------------------------