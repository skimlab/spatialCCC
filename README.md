
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
LRdb_m <- 
  get_LRdb("mouse") 

LR_pair_list <- c(
  "Apoe_Sdc1",
  "App_Dcc",
  "Cxcl12_Ackr3",
  "Jam3_Itgam",
  "Cx3cl1_Itga4")

LRdb_m <-
  dplyr::filter(LRdb_m, LR %in% LR_pair_list)

LRdb_m %>% 
  arrange(ligand_gene_symbol, receptor_gene_symbol)
#> # A tibble: 5 × 10
#>   LR     ligand_gene_symbol receptor_gene_symbol ligand_gene_id receptor_gene_id
#>   <chr>  <chr>              <chr>                         <dbl>            <dbl>
#> 1 Apoe_… Apoe               Sdc1                          11816            20969
#> 2 App_D… App                Dcc                           11820            13176
#> 3 Cx3cl… Cx3cl1             Itga4                         20312            16401
#> 4 Cxcl1… Cxcl12             Ackr3                         20315            12778
#> 5 Jam3_… Jam3               Itgam                         83964            16409
#> # ℹ 5 more variables: ligand_ensembl_protein_id <chr>,
#> #   receptor_ensembl_protein_id <chr>, ligand_ensembl_gene_id <chr>,
#> #   receptor_ensembl_gene_id <chr>, evidence <chr>
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
                               LRdb = LRdb_m)
tictoc::toc()
#> 18.077 sec elapsed

# future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 8,083 × 10
#>    src          dst       d norm.d LR    ligand receptor LRscore weight WLRscore
#>    <chr>        <chr> <dbl>  <dbl> <chr> <chr>  <chr>      <dbl>  <dbl>    <dbl>
#>  1 GCTATCGCGGC… CCCG…  138.   1.01 App_… App    Dcc        0.821  0.980    0.804
#>  2 CTGGCGGGAAT… CCCG…  138.   1.01 App_… App    Dcc        0.819  0.980    0.803
#>  3 CTGGAAATGGA… TAGA…  138    1.01 App_… App    Dcc        0.818  0.986    0.806
#>  4 GAGACTGATGG… TAGA…  138.   1.01 App_… App    Dcc        0.816  0.980    0.800
#>  5 TAGCTAAGTCC… TAGA…  137    1    App_… App    Dcc        0.816  1        0.816
#>  6 TGATCGGTTTG… TAGA…  138.   1.01 App_… App    Dcc        0.815  0.980    0.799
#>  7 TTCGGCTAGAG… CCCG…  138.   1.01 App_… App    Dcc        0.815  0.980    0.799
#>  8 CATCCTCTCAA… CATG…  138    1.01 App_… App    Dcc        0.814  0.986    0.802
#>  9 GACACGAGTTA… CCCG…  138    1.01 App_… App    Dcc        0.813  0.986    0.801
#> 10 GATCATTCCAA… TAGA…  138.   1.01 App_… App    Dcc        0.812  0.980    0.796
#> # ℹ 8,073 more rows
```

``` r
tictoc::tic("to_spatical_ccc_graph ...")

sp_col_data <- get_spatial_coords(spe)

ccc_graph_list <-
  to_spatial_ccc_graph_list(ccc_tbl, sp_col_data, workers = 4)

tictoc::toc()
#> to_spatical_ccc_graph ...: 6.397 sec elapsed
```

``` r
tictoc::tic()

ccc_graph_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list)

tictoc::toc()
#> 0.023 sec elapsed
```

``` r
ccc_graph_metrics_summary_df %>%
  arrange(graph_component_count)
#> # A tibble: 4 × 12
#>   LR         graph_n_nodes graph_n_edges graph_component_count graph_motif_count
#>   <chr>              <int>         <dbl>                 <dbl>             <int>
#> 1 App_Dcc             1989          3267                    42              8152
#> 2 Cx3cl1_It…          1400          1994                    65              4790
#> 3 Cxcl12_Ac…          1171          1418                   112              2436
#> 4 Jam3_Itgam          1191          1404                   118              2317
#> # ℹ 7 more variables: graph_diameter <dbl>, graph_un_diameter <dbl>,
#> #   graph_mean_dist <dbl>, graph_circuit_rank <dbl>, graph_reciprocity <dbl>,
#> #   graph_clique_num <int>, graph_clique_count <int>
```

``` r
tictoc::tic()

ccc_graph_group_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list, level = "group")

tictoc::toc()
#> 0.014 sec elapsed
```

``` r
ccc_graph_group_metrics_summary_df 
#> # A tibble: 337 × 12
#>    LR      group group_n_nodes group_n_edges group_adhesion group_motif_count
#>    <chr>   <int>         <int>         <dbl>          <dbl>             <int>
#>  1 App_Dcc     1          1095          1957              0              4947
#>  2 App_Dcc     2           183           310              0               807
#>  3 App_Dcc    20             8             9              0                15
#>  4 App_Dcc     5            52            92              0               217
#>  5 App_Dcc     4            70           115              0               292
#>  6 App_Dcc     6            45            59              0               151
#>  7 App_Dcc     9            27            30              0                78
#>  8 App_Dcc     3           136           275              0               645
#>  9 App_Dcc     7            30            48              0               111
#> 10 App_Dcc    17            12            11              0                26
#> # ℹ 327 more rows
#> # ℹ 6 more variables: group_diameter <dbl>, group_un_diameter <dbl>,
#> #   group_mean_dist <dbl>, group_girth <dbl>, group_circuit_rank <dbl>,
#> #   group_reciprocity <dbl>
```

``` r
LR_of_interest <- "App_Dcc"
```

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       tissue_img = imgRaster(spe),
                       node_color = "group_diameter",
                       node_size = 1,
                       edge_color = "group_diameter",
                       which_on_top = "edge")
#> Warning: The `guide` argument in `scale_*()` cannot be `FALSE`. This was deprecated in
#> ggplot2 3.3.4.
#> ℹ Please use "none" instead.
#> ℹ The deprecated feature was likely used in the spatialCCC package.
#>   Please report the issue at <https://github.com/dolchan/spatialCCC/issues>.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       # tissue_img = imgRaster(spe),
                       node_color = "group",
                       node_size = 0.1,
                       edge_color = "group_diameter",
                       edge_width = 0.1,
                       which_on_top = "edge")
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       # tissue_img = imgRaster(spe),
                       graph_layout = "stress",
                       node_color = "group",
                       edge_color = "group_diameter",
                       edge_width = 0.25,
                       which_on_top = "edge")
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />
