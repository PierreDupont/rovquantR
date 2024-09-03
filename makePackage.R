rm(list=ls())
library(devtools)
library(roxygen2)
#library(nimble, warn.conflicts = FALSE)
#library(basicMCMCplots)

if(Sys.info()['user'] == 'dturek') {
    baseDir <- '~/github/nimble/rovquantR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
    baseDir <- 'C:/My_documents/rovquantR/'          ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
    baseDir <- 'C:/Personal_Cloud/OneDrive/Work/rovquantR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
    baseDir <- 'C:/PROJECTS/rovquantR/'                       ## Richard
} else stop('unknown user')


if(!('makePackage.R' %in% list.files(baseDir))) stop('')


document(paste0(baseDir, 'rovquantR'))


if(.Platform$OS.type == 'windows') {
    system(paste0('R CMD build ', baseDir, 'rovquantR'))
} else {
    system(paste0('R CMD BUILD ', baseDir, 'rovquantR'))
}

check(paste0(baseDir, 'rovquantR'))


suppressMessages(try(remove.packages('rovquantR'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
if(.Platform$OS.type == 'windows') {
    system(paste0('R CMD INSTALL --build ', lastTarFile))
} else {
    system(paste0('R CMD install ', lastTarFile))
}

## now quit R

## now restart R

library(rovquantR)


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
library(testthat)
##test_package('nimbleSCR')
devtools::test('rovquantR')

# dDispersal_exp
# ?dDispersal_exp



## inspect package vignettes
(vignettes <- vignette(package = 'rovquantR'))
for(v in vignettes$results[, 'Item'])   print(vignette(v))



# ## time building of all vignettes:
# 
# library(nimbleSCR)
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
