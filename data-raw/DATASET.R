## -----------------------------------------------------------------------------
##
## Script name: 'rovquantR' dataset preparation.
##
## Purpose of script: This R script loads and cleans data for internal use in 'rovquantR'.
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
library(raster)
library(rmapshaper)



## -----------------------------------------------------------------------------

##-- Create a folder to contain this script, as well as other R scripts used 
##-- during package development ; Only done once when creating the package
#usethis::use_data_raw()

## -----------------------------------------------------------------------------

##-- Load and prepare translation data 
##-- This is saved in R/sysdata.rda as it is only used inside the 'translateForeignCharacters' function.
load(file.path(dir.git, "Analyses/CharacterTranslation.RData"))
head(fromto)

## -----------------------------------------------------------------------------

##-- Load and prepare spatial data (COUNTRIES and COUNTIES maps in our case)
##-- These are saved in ./data 
# GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) %>%
#   filter(.,
#          area > 80000000) %>%
#   st_crop(., st_bbox(raster::extent(c(-70000,1200000,5100000,8080000)))) 
GLOBALMAP <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/Scandinavia_border_33NNoLakes.shp")) %>% 
  st_simplify(., dTolerance =  500)



##-- POLYGONS OF SWEDEN & NORWAY
COUNTRIES <- st_read(file.path(dir.dropbox,"DATA/GISData/vegetation/Countries_waterHumans25000000m2_multimulti.shp")) %>%
  filter(ISO %in% c("SWE","NOR"), 
         area > 80000000) %>%
  group_by(ISO) %>%
  summarize()



##-- POLYGONS OF COUNTIES IN NORWAY 
COUNTIES_NOR <- st_read(file.path(dir.dropbox,"DATA/GISData/new_scandinavian_border/fylker-2024.shp")) %>%
  st_transform(crs = st_crs(COUNTRIES)) %>%
  rename(county = fylkesnavn,
         area = SHAPE_Area) %>%
  mutate(county = unlist(lapply(strsplit(county, " - "), function(x)x[1])),
         country = "NOR")
COUNTIES_NOR <- COUNTIES_NOR[ ,c("county","country","area")]

##-- POLYGONS OF COUNTIES IN SWEDEN 
COUNTIES_SWE <- st_read(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/rk_lan_07_WGS84.shp")) %>%
  rename(county = LANSNAMN,
         area = AREA_METER) %>%
  mutate(county = unlist(lapply(strsplit(county, " "), function(x)x[1])),
         country = "SWE")
COUNTIES_SWE <- COUNTIES_SWE[ ,c("county","country","area")]

##-- JOIN COUNTY POLYGONS
##-- Simplification is needed because of size for the package
COUNTIES <- rbind(COUNTIES_NOR, COUNTIES_SWE) %>% 
  ms_simplify( .,
               keep = 0.1,
               keep_shapes = FALSE)

##-- POLYGONS OF CARNIVORE MANAGEMENT REGIONS IN SWEDEN & NORWAY 
##-- Simplification at the end is needed because of large size otherwise
REGIONS <- sf::read_sf(file.path(dir.dropbox,"DATA/GISData/scandinavian_border/NOR_adm2_UTM33.shp"))  %>%
  group_by(NAME_1) %>%
  summarise() %>%
  mutate(
    ISO = "NOR",
    name = case_when(
    NAME_1 %in% c("Troms","Finnmark") ~ "Region 8",
    NAME_1 %in% c("Nordland") ~ "Region 7",
    NAME_1 %in% c("Sør-Trøndelag","Nord-Trøndelag","Møre og Romsdal") ~ "Region 6",
    NAME_1 %in% c("Hedmark") ~ "Region 5",
    NAME_1 %in% c("Akershus", "Ãstfold", "Oslo") ~"Region 4",
    NAME_1 %in% c("Oppland") ~ "Region 3",
    NAME_1 %in% c("Vestfold","Telemark","Buskerud","Aust-Agder") ~ "Region 2",
    NAME_1 %in% c("Hordaland","Sogn og Fjordane","Rogaland","Vest-Agder") ~ "Region 1",
    TRUE ~ NAME_1)) %>%
  group_by(name) %>%
  summarise(country = unique(ISO, na.rm = TRUE)) 
##-- SIMPLIFY COUNTIES
REGIONS <- sf::st_simplify(REGIONS, preserveTopology = F, dTolerance = 500)
REGIONS$region <- 1:8

# plot(st_geometry(REGIONS2))
# text(REGIONS2,labels = REGIONS2$region)
# REGIONS_NOR <- sf::st_read(file.path( dir.dropbox,
#                                   "DATA/GISData/NorwegianManagementRegions/rovviltregioner2024.shp"))

##-- HABITAT RASTERS AT DIFFERENT RESOLUTIONS (REFERENCE RASTERS)
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/Habitat20kmNewNorCounties.RData"))
load(file.path(dir.dropbox, "DATA/GISData/spatialDomain/HabitatAllResolutionsNewNorCounties.RData"))


##-- Save necessary data in the right folder (./data)
use_data(fromto, internal = TRUE, overwrite = TRUE)
use_data(COUNTIES, overwrite = TRUE)
use_data(REGIONS, overwrite = TRUE)
use_data(COUNTRIES, overwrite = TRUE)
use_data(GLOBALMAP, overwrite = TRUE)
use_data(habitatRasterResolution, overwrite = TRUE)
use_data(habitatRasters, overwrite = TRUE)
# use_data(COUNTRIESWaterHumans, overwrite = TRUE)



##------------------------------------------------------------------------------