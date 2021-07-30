
<!-- README.md is generated from README.Rmd. Please edit that file -->

# expressalyzr

<!-- badges: start -->

<!-- badges: end -->

The goal of expressalyzr is to automate the analysis of flow cytometry
data.

## Installation

You can install the expressalyzr from
[GitHub](https://github.com/freitim/expressalyzr) with:

``` r
devtools::install_github("freitim/expressalyzr", ref = "main")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(expressalyzr)
library(ggcyto)

data_dt <- run_pipeline(data_path = "path to your data", view_config = TRUE)
```
