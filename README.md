
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialCCC

<!-- badges: start -->
<!-- badges: end -->

The goal of spatialCCC is to investigate cell-cell signaling, by
analyzing ligand-receptor interactions in spatial transcriptomic data.

## Installation

You can install the development version of spatialCCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dolchan/spatialCCC")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SpatialExperiment)
library(scater)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(patchwork)

library(ggraph)
library(tidygraph)
library(RColorBrewer)
```

``` r
library(spatialCCC)
## basic example code
```

## Example: Visium Spatial Transcriptomic Data

``` r
data_dir <- file.path("example", "visium_tutorial")
spe <- SpatialExperiment::read10xVisium(samples = data_dir, type = "HDF5", data = "filtered")
```

## Ligand-Receptor Database: built-in

``` r
LRdb_mouse <- get_LRdb("mouse")
```

## Log-Normalize

``` r
spe <- logNormCounts(spe)
```

## Compute Cell-Cell Communications over ligand-receptor pairs

``` r
future::plan(future::multisession, workers = 6)

ccc_tbl <- compute_spatial_ccc(spe = spe, 
                               assay_name = "logcounts",
                               LRdb = LRdb_mouse)

future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 2,089,252 × 10
#>    src          dst       d norm.d LR    ligand receptor LRscore weight WLRscore
#>    <chr>        <chr> <dbl>  <dbl> <chr> <chr>  <chr>      <dbl>  <dbl>    <dbl>
#>  1 ACCGCGGTGGA… TCTT…  138.   1.00 Apoe… Apoe   Lrp1       0.880  0.992    0.873
#>  2 GCAACAGCAGT… TCTT…  138.   1.01 Apoe… Apoe   Lrp1       0.879  0.980    0.861
#>  3 TGCTGGTTGGA… TTGT…  138.   1.01 Apoe… Apoe   Lrp1       0.878  0.980    0.860
#>  4 TCGTCCGCTGG… TTGT…  137    1    Apoe… Apoe   Lrp1       0.877  1        0.877
#>  5 CCATGCTCTGC… TTGT…  138.   1.01 Apoe… Apoe   Lrp1       0.876  0.980    0.858
#>  6 ATCGCCAGTCA… TCTT…  138    1.01 Apoe… Apoe   Lrp1       0.875  0.986    0.863
#>  7 TTCGCTATCTG… GGCT…  138.   1.01 Apoe… Apoe   Lrp1       0.875  0.980    0.857
#>  8 TTCTAGGCCAA… TTGT…  138    1.01 Apoe… Apoe   Lrp1       0.873  0.986    0.861
#>  9 TAACAAAGGGA… GCTT…  138.   1.01 Apoe… Apoe   Lrp1       0.873  0.987    0.861
#> 10 TTCTTGTAACC… GCTT…  138.   1.01 Apoe… Apoe   Lrp1       0.873  0.987    0.861
#> # ℹ 2,089,242 more rows
```
