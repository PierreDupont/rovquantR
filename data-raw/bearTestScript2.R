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



##------------------------------------------------------------------------------
## ----- IV. PREPARE OPSCR DATA ------

##-- paths
species = "bear"
data_dir = data_dir
working_dir = working_dir
##-- data
years = 2020:2023
sex = c("Hann","Hunn")
aug.factor = 2
sampling.months = list(4,5,6,7,8,9,10,11)
##-- habitat
habitat.res = 20000
buffer.size = 50000
max.move.dist = 250000
##-- detectors
detector.res = 5000
subdetector.res = 1000
max.det.dist = 70000
resize.factor = 1
##-- miscellanious
print.report = F


