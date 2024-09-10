## -----------------------------------------------------------------------------
##
## Script name: 'rovquantR' dataset preparation.
##
## Purpose of script: This R script loads and clean data for internal use in 'rovquantR'.
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 2024-09-04
##
## Copyright (c) Applied Quantitative Ecology Group (AQEG), 2024
## Faculty of Environmental Sciences and Natural Resource Management (MINA)
## Norwegian University of Life Sciences (NMBU), Ås, Norway 
##
##
## -----------------------------------------------------------------------------
##
## Notes: this script is for internal use only!
##   
## -----------------------------------------------------------------------------
rm(list = ls())
gc()



##-- Identify user and set corresponding DropBox directory
if(Sys.info()['user'] == 'pidu') {
  dir.git <- "C:/My_documents/RovQuant"
  dir.dropbox <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'pierredupont') {
  dir.git <- "/Users/pierredupont/Documents/RovQuant"
  dir.dropbox <- "/Users/pierredupont/Dropbox (AQEG)/AQEG Team Folder/RovQuant/"
}
if(Sys.info()['user'] == 'cymi') {
  dir.git <- "C:/My_documents/rovquant/"
  dir.dropbox <- "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant" 
}
if(Sys.info()['user'] == 'richbi') {
  dir.git <- "C:/Users/richbi/OneDrive - Norwegian University of Life Sciences/PROJECTS/RovQuant"
  dir.dropbox <- "C:/Users/richbi/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'seasunci') {
  dir.git <- "C:/Users/seasunci/02_RovQuant/RovQuant"
  dir.dropbox <- "C:/Users/seasunci/AQEG Dropbox/AQEG Team Folder/RovQuant"
}



##-- Create a folder to contain this script, as well as other R scripts used during package development
library(usethis)
#usethis::use_data_raw()



##-- Load and prepare translation data 
##-- This is saved in R/sysdata.rda as it is only used inside the 'translateForeignCharacters' function.
load(file.path(dir.git, "Analyses/CharacterTranslation.RData"))
head(fromto)
usethis::use_data(fromto, internal = TRUE, overwrite = TRUE)



##-- Load and prepare other necessary data (COUNTRIES and COUNTIES maps in our case)
##-- These are saved in ./data 
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/Habitat_shp.RData"))
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/Habitat20kmNewSweCounties.RData"))
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/HabitatAllResolutionsNewSweCounties.RData"))


##-- Check for non-ASCII characters
##-- Use encoding() to learn the current encoding of the elements in a 
##-- character vector and functions such as enc2utf8() or iconv() to convert 
##-- between encodings.
COUNTIES$NAME_1 <- stri_trans_general(COUNTIES$NAME_1, "Latin-ASCII")
stri_enc_isutf8(COUNTIES$NAME_1)


##-- Save necessary data in the right folder (./data)
use_data(COUNTIES, overwrite = TRUE)
use_data(COUNTRIES, overwrite = TRUE)
use_data(COUNTRIESWaterHumans, overwrite = TRUE)
use_data(GLOBALMAP, overwrite = TRUE)
use_data(habitatRasterResolution, overwrite = TRUE)
use_data(habitatRasters, overwrite = TRUE)

