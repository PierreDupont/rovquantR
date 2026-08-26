##------------------------------------------------------------------------------
##
## Script name: RovQuant WOLVERINE analysis 2026 - test script
##
## This R script tests the packaged version of the Wolverine analysis of 2025,
## known as "the Chaotic Estimation!" (cf. CM).  OPSCR and SCR models were ran using the model
## 54.Cleaned. ("54.Cleaned2025TestPDScript.R").Results from the OPSCR were extracted 
## using "PlotWolverine54Cleaned2025.R" and "PlotWolverine54Cleaned2025Snap_ASPcm.R" 
## for the SCR model. Figures and tables that combined results from the OPSCR 2024 
## and the SCR model were created using "CombineAndPlotFigureTable.R".
##
## This file tests if the rovquantR version of those functions can produce the 
## same nimble inputs and outputs than the script "wolverineWIPscript.R".
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 11/08/2026
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2026
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
#pak::pak("PierreDupont/rovquantR@devel")
devtools::install_github("PierreDupont/rovquantR@devel")



## ------ LOAD REQUIRED LIBRARIES ------

library(rovquantR)
library(nimbleSCR)

library(raster)
library(sf)
library(dplyr)
library(ggplot2)


##------------------------------------------------------------------------------

## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY 
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/wolverine/2025/Test_rovquantR"



##------------------------------------------------------------------------------

## ----- II. CLEAN NGS DATA -----

cleanRovbaseData( 
  species = "wolverine",
  years = 2012:2025,
  data.dir = data.dir,
  working.dir = working.dir)


##-- Load the most recent clean wolverine data from RovBase
myFullData.sp <- readMostRecent( 
  path = file.path(working.dir, "data"),
  pattern = "CleanData_wolverine",
  extension = ".RData")

##-- Checks
dim(myFullData.sp$alive[myFullData.sp$alive$Country_sf %in% c("(N)","(S)"), ])
table(myFullData.sp$alive$Year[myFullData.sp$alive$Country_sf %in% c("(N)","(S)")])
dim(myFullData.sp$dead.recovery)
table(myFullData.sp$dead.recovery$Year)
table(myFullData.sp$dead.recovery$Death_cause)
table(myFullData.sp$dead.recovery$Death_method)



##------------------------------------------------------------------------------

## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  species = "wolverine",
  years = 2015:2024,
  data.dir = data.dir,
  working.dir = working.dir)

## Load and explore input data
load(file.path(working.dir,"nimbleInFiles/female/nimbleInput_2026-08-20_female_1.RData")) 
nimData_NEW <- nimData
nimConstants_NEW <- nimConstants
lapply(nimData_NEW,dim)
lapply(nimData_NEW,sum)



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


##------------------------------------------------------------------------------