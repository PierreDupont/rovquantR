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
## Date Created: 05/04/2025
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2025
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
##------------------------------------------------------------------------------
##
## Notes: 
## This is based on 'rovquantR' beta version 0.2
##   
##------------------------------------------------------------------------------

## ------ CLEAR-UP ENVIRONMENT ------

rm(list = ls())
gc()

 
## ------ INSTALL 'rovquantR' FROM GITHUB ------

## Ctrl + Shift + F10 (to restart R session)
devtools::install_github("PierreDupont/rovquantR")



## ------ LOAD REQUIRED LIBRARIES ------

library(rovquantR)
?GetDensity
library(nimbleSCR)


##------------------------------------------------------------------------------
## ------ I. SET-UP WORKING ENVIRONMENT ------

##-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2024/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working.dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2024/Analysis_test"


##------------------------------------------------------------------------------
## ----- II. CLEAN NGS DATA -----

cleanRovbaseData( 
  species = "bear",
  years = 2020:2025,
  data.dir = data.dir,
  working.dir = working.dir,
  two.sex = TRUE,
  print.report = TRUE)


##------------------------------------------------------------------------------
## ----- III. PREPARE OPSCR DATA ------

makeRovquantData(    
  species = "bear",
  data.dir = data.dir,
  working.dir = working.dir)


##------------------------------------------------------------------------------
## ----- IV. CHECK OPSCR DATA ------

for (s in c("female", "male")) {

  ##-- Load first input
  input <- list.files( file.path( working.dir, "nimbleInFiles", s ),
                       full.names = TRUE)
  load(input[1])

  ##-- Table of #of detections 
  detsPerYear <- apply(nimData$y.alive,c(1,3), function(x){
    ifelse(x[1] > 0, sum(x[2:(x[1]+1)]), 0)
  })
  numSamplesPerYear <- colSums(detsPerYear)
  numIdsPerYear <- colSums(detsPerYear > 0)
  
  ##-- Build nimble model object
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        inits = nimInits,
                        data = nimData,
                        check = FALSE,
                        calculate = FALSE)
  model$calculate()
  # which( is.infinite(model$logProb_sxy), arr.ind = TRUE)
  # model$sxy[15, , ]
  # colSums(nimData$y.alive[15,,]>0)
  # model$sxy[15, ,4] <- (model$sxy[15, ,3]*0.2 +  model$sxy[15, ,5]*0.8)
}#




##------------------------------------------------------------------------------
## ----- V. FIT ROVQUANT MODELS ------

## -----   1. Females ------

##-- List all prepared input files
inputFiles <- list.files( file.path(working.dir, "nimbleInFiles/female"),
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
inputFiles <- list.files( file.path(working.dir, "nimbleInFiles/male"),
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

åcmodel <- compileNimble(model)
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

# Rmd.template <- system.file("rmd", "RovQuant_FullReport.Rmd", package = "rovquantR")

processRovquantOutput(   
  species = "Brown bear"
  ,
  data.dir = data.dir
  ,
  working.dir = working.dir
  ,
  nburnin = 0
  ,
  niter = 500
  ,
  extraction.res = 5000
  ,
  print.report = TRUE
  # ,
  # Rmd.template = Rmd.template
  ,
  #output.dir = working.dir
  ,
  overwrite = FALSE
)


## TODAY, I NEED TO:
## - Add the norwegian summary to the full report?
## - Check numbers in the overleaf report

##------------------------------------------------------------------------------
