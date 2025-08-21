## -----------------------------------------------------------------------------
##
## Script name: 'rovquantR' package building.
##
## Purpose of script: This R script builds and checks the 'rovquantR' package.
##
## Author: 2024-09-04
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


## ------ IMPORT REQUIRED LIBRARIES ------
library(devtools)
library(roxygen2)


## ------ SET PACKAGE DIRECTORY ------
if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/rovquantR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/My_documents/rovquantR/'                   ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/rovquantR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/rovquantR/'                       ## Richard
} else stop('unknown user')


## ------ SET-UP README.rmd FILE ------
#usethis::use_readme_rmd()


## ------ DOCUMENT PACKAGE ------
document(baseDir)
Rcpp::compileAttributes()           # this updates the Rcpp layer from C++ to R
roxygen2::roxygenize(roclets="rd") 

#useDynLib(rovquantR, .registration = TRUE)


## ------ BUILD PACKAGE ------
system(paste0('R CMD build ', baseDir))


## ------ CHECK PACKAGE ------
check(baseDir)


## ------ INSTALL PACKAGE FROM .tar ------
suppressMessages(try(remove.packages('rovquantR'), silent = TRUE))
tarFiles <- grep( pattern = '\\.tar\\.gz',
                  x = list.files(baseDir, include.dirs = TRUE),
                  value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
system(paste0('R CMD INSTALL --build ', lastTarFile))


## ------ (ALTERNATIVE) INSTALL PACKAGE FROM GITHUB ------
devtools::install_github("PierreDupont/rovquantR")


## ------ NOW QUIT R -----
q()


## ------ NOW RESTART R ------


## ------ TRY TO LOAD PACKAGE ------
library(rovquantR)

rovquantR::GetDensity()



## ------ FURTHER TESTING ------
# # Run testthat tests
# if(Sys.info()['user'] == 'dturek') {
#     baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
# } else if(Sys.info()['user'] == 'pidu') {
#     baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
# } else if(Sys.info()['user'] == 'cymi') {
#     baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
# } else if(Sys.info()['user'] == 'arichbi') {
#     baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
# } else stop('unknown user')
# 
# # Utility function (fast!)
# source(file.path(baseDir,'nimbleSCR/test/testthat/testUtilityFunctions.R'))
# # distributions (slow)
# source(file.path(baseDir,'nimbleSCR/test/testthat/testDistributionFunctions.R'))
# library(testthat)
# ##test_package('nimbleSCR')
# devtools::test('rovquantR')
# 
# # dDispersal_exp
# # ?dDispersal_exp
# 
# 
# 
# ## Inspect package vignettes
# (vignettes <- vignette(package = 'rovquantR'))
# for(v in vignettes$results[, 'Item'])   print(vignette(v))
# 
#
# ## time building of all vignettes:
# library(rovquantR)
# library(knitr)
# library(rmarkdown)
# 
# setwd('~/github/nimble/nimbleSCR/nimbleSCR/vignettes')
# f <- list.files()
# f2 <- grep('\\.[Rr]md$', f, value = TRUE)
# times <- numeric(length(f2))
# names(times) <- f2
# 
# for(iii in 1:length(f2)) {
#     message(f2[iii], ':')
#     tm <- system.time(render(f2[iii]))
#     times[iii] <- tm[3]
# }
# 
# times
# sum(times) / 60   ## total in minutes