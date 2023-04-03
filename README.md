
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
LRdb_mouse <- 
  get_LRdb("mouse") %>%
  dplyr::slice_sample(n = 50)
```

## Log-Normalize

``` r
spe <- logNormCounts(spe)
```

## Compute Cell-Cell Communications over ligand-receptor pairs

``` r
# For full LRdb analysis, future::plan can be used 
#   for parallelization
# future::plan(future::multisession, workers = 4)

tictoc::tic()
ccc_tbl <- compute_spatial_ccc(spe = spe, 
                               assay_name = "logcounts",
                               LRdb = LRdb_mouse)
tictoc::toc()
#> 16.888 sec elapsed

# future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 43,276 × 10
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
#> # ℹ 43,266 more rows
```

``` r
tictoc::tic("to_spatical_ccc_graph ...")

sp_col_data <- get_spatial_coords(spe)

ccc_graph_list <-
  to_spatial_ccc_graph_list(ccc_tbl, sp_col_data, workers = 4)

tictoc::toc()
#> to_spatical_ccc_graph ...: 13.731 sec elapsed
```

``` r
tictoc::tic()

ccc_graph_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list)

tictoc::toc()
#> 0.062 sec elapsed
```

``` r
ccc_graph_metrics_summary_df 
#> # A tibble: 16 × 12
#>    LR        graph_n_nodes graph_n_edges graph_component_count graph_motif_count
#>    <chr>             <int>         <dbl>                 <dbl>             <int>
#>  1 Trf_Tfr2           1227          1520                    83              3481
#>  2 Jag1_Not…           346           263                   108               210
#>  3 Igf2_Igf…          1074          1265                   134              1941
#>  4 Col4a1_C…           583           516                   138               570
#>  5 Nlgn1_Nr…          2580          9213                     3             20396
#>  6 Icam2_It…           151            99                    56                45
#>  7 Inhba_Ac…           314           349                    58               521
#>  8 Apoe_Lrp1          2701         14893                     3             27822
#>  9 Mdk_Itgb1          2124          5200                    17             11114
#> 10 Fgf13_Fg…          2363          5649                    16             13251
#> 11 Adam12_I…            83            52                    33                21
#> 12 Wnt5b_Lr…           345           278                    91               295
#> 13 Fgf2_Not…           424           340                   108               328
#> 14 Vip_Vipr1           609           684                    77              1092
#> 15 Sema4d_M…           647           787                    65              1393
#> 16 Efnb1_Ep…          1457          2168                    75              4431
#> # ℹ 7 more variables: graph_diameter <dbl>, graph_un_diameter <dbl>,
#> #   graph_mean_dist <dbl>, graph_circuit_rank <dbl>, graph_reciprocity <dbl>,
#> #   graph_clique_num <int>, graph_clique_count <int>
```

``` r
tictoc::tic()

ccc_graph_group_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list, level = "group")

tictoc::toc()
#> 0.053 sec elapsed
```

``` r
ccc_graph_group_metrics_summary_df 
#> # A tibble: 1,065 × 12
#>    LR       group group_n_nodes group_n_edges group_adhesion group_motif_count
#>    <chr>    <int>         <int>         <dbl>          <dbl>             <int>
#>  1 Trf_Tfr2     5            45            58              0               109
#>  2 Trf_Tfr2    31            13            12              0                31
#>  3 Trf_Tfr2    45             7             6              0                15
#>  4 Trf_Tfr2    79             4             3              0                 3
#>  5 Trf_Tfr2    23            15            18              0                45
#>  6 Trf_Tfr2    46             7             6              0                15
#>  7 Trf_Tfr2    16            22            42              0                95
#>  8 Trf_Tfr2    18            20            29              0                64
#>  9 Trf_Tfr2     4            50            59              0               154
#> 10 Trf_Tfr2     3            57            79              0               173
#> # ℹ 1,055 more rows
#> # ℹ 6 more variables: group_diameter <dbl>, group_un_diameter <dbl>,
#> #   group_mean_dist <dbl>, group_girth <dbl>, group_circuit_rank <dbl>,
#> #   group_reciprocity <dbl>
```
