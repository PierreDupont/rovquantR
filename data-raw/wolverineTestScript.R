##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE OPSCR analysis 
##
## Purpose of script: 
## This R script performs:
## 1. the initial cleaning of the wolverine NGS data downloaded from RovBase.3.0
## 2. the data preparation for the RovQuant OPSCR analysis with the 'nimbleSCR' package
## 3. the model fitting using 'nimble' and 'nimbleSCR'
## 4. the post-processing of the MCMC output
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 26/09/2025
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2025
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
##------------------------------------------------------------------------------
##
## Notes: 
## This is based on 'rovquantR' beta version 0.1
##   
##------------------------------------------------------------------------------
rm(list = ls())
gc()


## ------ IMPORT REQUIRED LIBRARIES ------

## Ctrl + Shift + F10 (to restart R session)
devtools::install_github("PierreDupont/rovquantR@devel")


## ------ LOAD REQUIRED LIBRARIES ------

library(rovquantR)
library(nimbleSCR)



##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY 
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test.0.2"



##------------------------------------------------------------------------------
## ----- II. CLEAN NGS DATA -----

cleanRovbaseData( 
  species = "wolverine",
  years = 2014:2023,
  data.dir = data.dir,
  working.dir = working.dir)


##-- Load the most recent clean wolverine data from RovBase
myFullData.sp <- readMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolverine",
  extension = ".RData")

##-- Checks
dim(myFullData.sp$alive[myFullData.sp$alive$Country_sf %in% c("(N)","(S)"), ])
table(myFullData.sp$alive$Season[myFullData.sp$alive$Country_sf %in% c("(N)","(S)")])
dim(dead.recovery)
table(dead.recovery$Season)
table(dead.recovery$Death_cause)
table(dead.recovery$Death_method)



##------------------------------------------------------------------------------
## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  species = "wolverine",
  data.dir = data.dir,
  working.dir = working.dir)



##------------------------------------------------------------------------------
## ----- IV. FIT ROVQUANT MODELS ------

## -----   1. Females ------

##-- List all prepared input files
inputFiles <- list.files(file.path( working.dir, "nimbleInFiles/female"),
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
                          path = file.path(working.dir,"nimbleOutfiles/female")))



## -----   2. Males ------

##-- List all prepared input files
inputFiles <- list.files(file.path(working.dir, "nimbleInFiles/male"),
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
                          path = file.path(working.dir,"nimbleOutfiles/male")))



##------------------------------------------------------------------------------
## ----- V. PROCESS ROVQUANT OUTPUT ------

processRovquantOutput(   
  species = "Wolverines",
  data.dir = data.dir,
  working.dir = working.dir)


## -----------------------------------------------------------------------------
