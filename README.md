
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rovquantR <a href="https://www.nmbu.no/en/research/projects/rovquant"><img src="inst/images/RovQuant.png" align="right" height="180"/></a>

<h4 align="center">

A user-friendly interface to reproduce RovQuant’s analyses.

</h4>

<!-- badges: start -->

<p align="center">

<a href="https://github.com/PierreDupont/rovquantR/commits/master">
<img src="https://img.shields.io/github/last-commit/PierreDupont/rovquantR.svg?style=flat-square&amp;logo=github&amp;logoColor=white" alt="GitHub last commit"/>
<a href="https://github.com/PierreDupont/rovquantR/issues">
<img src="https://img.shields.io/github/issues-raw/PierreDupont/rovquantR.svg?style=flat-square&amp;logo=github&amp;logoColor=white" alt="GitHub issues"/>
<a href="https://github.com/PierreDupont/rovquantR/pulls">
<img src="https://img.shields.io/github/issues-pr-raw/PierreDupont/rovquantR.svg?style=flat-square&amp;logo=github&amp;logoColor=white" alt="GitHub pull requests"/>

</p>

<!-- badges: end -->

<!-- <p align="center"> -->

<!--   <a href="#installation">Installation</a> • -->

<!--   <a href="#example">Example</a> • -->

<!--   <a href="#goodrm">Good Readme</a>  -->

<!-- </p> -->

------------------------------------------------------------------------

The goal of **rovquantR** is to provide user-friendly `R` and `nimble`
functions to facilitate the reproduction of SCR and OPSCR analyses of
the Scandinavian large carnivore monitoring data originally performed by
the [Applied Quantitative Ecology Group
(AQEG)](https://www.nmbu.no/en/research/groups/applied-quantitative-ecology-group-aqeg)
during project
[RovQuant](https://www.nmbu.no/forside/en/projects/rovquant).

This package represents a collaborative effort between the
[nimble](https://r-nimble.org/) development team and Project
[RovQuant](https://www.nmbu.no/forside/en/projects/rovquant).

## Installation

You can install and load the development version of **rovquantR** from
[GitHub](https://github.com/) with:

``` r
##-- Install "rovquantR" package
# install.packages("devtools")
devtools::install_github("PierreDupont/rovquantR")

##-- Load "rovquantR" package
library(rovquantR)
```

## Example

The `rovquantR` package can then be used to set-up the directory
structure for a new analysis:

``` r
#-- DATA DIRECTORY
##-- Directory containing the raw data necessary for the analysis
##-- (NB: This is NOT the working directory; NOTHING SHOULD BE SAVED/WRITTEN IN THIS DIRECTORY)
data_dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/Data"

##-- WORKING DIRECTORY (= main folder for the analysis)
working_dir <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant/bear/2023/RovQuant_test/test2"

##-- Create folder structure for the analysis
makeDirectories(path = working_dir)
```

<!-- ## Good Readme -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
