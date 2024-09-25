## -----------------------------------------------------------------------------
##
## Script name: 'rovquantR' vignettes preparation.
##
## Purpose of script: This R script sets up the 'rovquantR' package's vignettes
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 2024-09-25
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


library(usethis)

usethis::use_vignette("Bear_Vignette")

devtools::build_rmd("vignettes/Bear_Vignette.Rmd")