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

##-- Libraries
library(usethis)
library(stringi)
library(dplyr)
library(sf)


##-- Cresfheaders##-- Create a folder to contain this script, as well as other R scripts used 
##-- during package development ; Only done once when creating the package
#usethis::use_data_raw()


##-- Load and prepare translation data 
##-- This is saved in R/sysdata.rda as it is only used inside the 'translateForeignCharacters' function.
load(file.path(dir.git, "Analyses/CharacterTranslation.RData"))
head(fromto)


##-- Load and prepare spatial data (COUNTRIES and COUNTIES maps in our case)
##-- These are saved in ./data 
#load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/Habitat_shp.RData"))
GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) %>%
  filter(.,
         area > 80000000,
         ISO %in% c("SWE","NOR")) %>%
  st_crop(., st_bbox(extent(c(-70000,1200000,5100000,8080000)))) 
# GLOBALMAP <- GLOBALMAP[GLOBALMAP$area > 80000000, ]
# GLOBALMAP <- st_crop(GLOBALMAP, st_bbox(extent(c(-70000,1200000,5100000,8080000))))


##-- POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- GLOBALMAP %>%
  filter(ISO %in% c("SWE","NOR")) %>%
  group_by(ISO) %>%
  summarize()


##-- POLYGONS OF COUNTIES IN SWEDEN & NORWAY
##-- Simplification at the end is needed because of large size otherwise
COUNTIES <- rbind(
  st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp")),
  st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/SWE_adm2_UTM33.shp"))) %>%
  group_by(NAME_1) %>%
  summarize() %>%
  ## Simplification because of large size otherwise (from 13.4 MB to 400 KB)
  st_simplify(., dTolerance = 500) %>%
  ##-- Check for non-ASCII characters
  ##-- Use encoding() to learn the current encoding of the elements in a 
  ##-- character vector and functions such as enc2utf8() or iconv() to convert 
  ##-- between encodings.
  mutate(.,NAME_1 = stri_trans_general(NAME_1, "Latin-ASCII"))
## Check encoding of county names
stri_enc_isutf8(COUNTIES$NAME_1)


##-- HABITAT RASTERS AT DIFFERENT RESOLUTIONS (REFERENCE RASTERS)
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/Habitat20kmNewNorCounties.RData"))
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/HabitatAllResolutionsNewNorCounties.RData"))


##-- Save necessary data in the right folder (./data)
use_data(fromto, internal = TRUE, overwrite = TRUE)
use_data(COUNTIES, overwrite = TRUE)
use_data(COUNTRIES, overwrite = TRUE)
#use_data(COUNTRIESWaterHumans, overwrite = TRUE)
use_data(GLOBALMAP, overwrite = TRUE)
use_data(habitatRasterResolution, overwrite = TRUE)
use_data(habitatRasters, overwrite = TRUE)

