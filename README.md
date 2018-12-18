
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lungct

[![Travis-CI Build
Status](https://travis-ci.org/muschellij2/lungct.svg?branch=master)](https://travis-ci.org/muschellij2/lungct)

The goal of lungct is to process Lung CT Scans.

## Installation

You can install lungct from github with:

``` r
# install.packages("devtools")
devtools::install_github("muschellij2/lungct")
```

## Example

This is a basic example of lung segmentation:

``` r
library(lungct)
filename = "example.nii.gz"
reg = segment_lung(filename)
```
