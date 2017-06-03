# lungct

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
