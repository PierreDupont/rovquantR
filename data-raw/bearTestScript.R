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
## This is based on 'rovquantR' beta version 0.2
##   
##------------------------------------------------------------------------------
rm(list = ls())
gc()


## ------ IMPORT REQUIRED LIBRARIES ------

devtools::install_github("PierreDupont/rovquantR")
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

cleanRovbaseData( 
  species = "bear",
  years = 2012:2024,
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  ##-- paths
  species = "bear",
  #years = 2020:2024,
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- IV. FIT ROVQUANT MODELS ------

## -----   1. Females ------

##-- List all prepared input files
inputFiles <- list.files(file.path( working.dir, "nimbleInFiles/Hunn"),
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
## ----- V. PROCESS ROVQUANT OUTPUT ------

processRovquantOutput(   
  ##-- paths
  species = "bear",
  data.dir = data.dir,
  working.dir = working.dir)


##------------------------------------------------------------------------------
