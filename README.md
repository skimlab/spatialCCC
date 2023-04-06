
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialCCC <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/dolchan/spatialCCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dolchan/spatialCCC/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goal of **spatialCCC** package is to investigate cell-cell
signaling, by analyzing ligand-receptor interactions in spatial
transcriptomic data.

## Installation

You can install the development version of spatialCCC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dolchan/spatialCCC")
```

# Example

This is a basic example which shows you a basic workflow of the package:

### Load necessary packages

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

library(ggtree)
```

### Load spatialCCC package

``` r
library(spatialCCC)
## basic example code
```

## Ligand-Receptor Database

### Load built-in LRdb

``` r
LRdb_mm <- 
  get_LRdb("mouse")

LRdb_m <-
  LRdb_mm %>%
  slice_sample(n = 100)

LR_pair_list <- c(
  "Apoe_Sdc1",
  "App_Dcc",
  "Cxcl12_Ackr3",
  "Jam3_Itgam",
  "Cx3cl1_Itga4")

LRdb_m <-
  rbind(LRdb_m,
        dplyr::filter(LRdb_mm, LR %in% LR_pair_list)) %>%
  distinct()

LRdb_m %>% 
  arrange(ligand_gene_symbol, receptor_gene_symbol)
#> # A tibble: 104 × 10
#>    LR    ligand_gene_symbol receptor_gene_symbol ligand_gene_id receptor_gene_id
#>    <chr> <chr>              <chr>                         <dbl>            <dbl>
#>  1 Adam… Adam10             Fcer2a                        11487            14128
#>  2 Adam… Adam10             Gpnmb                         11487            93695
#>  3 Adam… Adam17             Cd9                           11491            12527
#>  4 Adam… Adam2              Itgb7                         11495            16421
#>  5 Agt_… Agt                Agtr1a                        11606            11607
#>  6 Agtr… Agtrap             Rack1                         11610            14694
#>  7 Angp… Angpt1             Tek                           11600            21687
#>  8 Angp… Angpt2             Itgb2                         11601            16414
#>  9 Angp… Angptl1            Tek                           72713            21687
#> 10 Angp… Angptl8            Gpihbp1                      624219            68453
#> # ℹ 94 more rows
#> # ℹ 5 more variables: ligand_ensembl_protein_id <chr>,
#> #   receptor_ensembl_protein_id <chr>, ligand_ensembl_gene_id <chr>,
#> #   receptor_ensembl_gene_id <chr>, evidence <chr>
```

## Load Visium spatial transcriptomic data

``` r
data_dir <- file.path("example", "visium_tutorial")
spe <- SpatialExperiment::read10xVisium(samples = data_dir, type = "HDF5", data = "filtered")
```

### Log-Normalize

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
#> 7.017 sec elapsed

# future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 101,757 × 10
#>    src          dst       d norm.d LR    ligand receptor LRscore weight WLRscore
#>    <chr>        <chr> <dbl>  <dbl> <chr> <chr>  <chr>      <dbl>  <dbl>    <dbl>
#>  1 AGAAGTGATTC… AATG…  138    1.01 Calm… Calm1  Kcnq3      0.858  0.986    0.846
#>  2 TTGCATGCTGA… GTGC…  138.   1.01 Calm… Calm1  Kcnq3      0.858  0.980    0.840
#>  3 CCCAACATACG… GTGC…  138.   1.00 Calm… Calm1  Kcnq3      0.857  0.992    0.850
#>  4 CGGTGCAGATA… GTGC…  138    1.01 Calm… Calm1  Kcnq3      0.857  0.986    0.845
#>  5 GTGGTTTCCGC… GTGC…  138    1.01 Calm… Calm1  Kcnq3      0.857  0.986    0.845
#>  6 TCTGGGTAGCG… ACAT…  138.   1.01 Calm… Calm1  Kcnq3      0.857  0.987    0.845
#>  7 GGGAGGATGCC… CCAA…  138    1.01 Calm… Calm1  Kcnq3      0.855  0.986    0.843
#>  8 TTGCATGCTGA… ACAT…  137    1    Calm… Calm1  Kcnq3      0.855  1        0.855
#>  9 AGATACCGGTG… ACAT…  138.   1.01 Calm… Calm1  Kcnq3      0.855  0.980    0.837
#> 10 TTCTTATCCGC… AGTT…  138.   1.00 Calm… Calm1  Kcnq3      0.855  0.992    0.848
#> # ℹ 101,747 more rows
```

### Convert CCC table to CCC graph

The conversion also adds various graph metrics to each CCC graph.

``` r
tictoc::tic("to_spatical_ccc_graph ...")

sp_col_data <- get_spatial_coords(spe)

ccc_graph_list <-
  to_spatial_ccc_graph_list(ccc_tbl, sp_col_data, workers = 4)

tictoc::toc()
#> to_spatical_ccc_graph ...: 12.265 sec elapsed
```

summarize_ccc_graph_metrics() summarize those graph metrics for each LR
pair.

``` r
tictoc::tic()

ccc_graph_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list)

tictoc::toc()
#> 0.066 sec elapsed
```

``` r
ccc_graph_metrics_summary_df %>%
  arrange(graph_component_count)
#> # A tibble: 48 × 12
#>    LR        graph_n_nodes graph_n_edges graph_component_count graph_motif_count
#>    <chr>             <int>         <dbl>                 <dbl>             <int>
#>  1 Calm1_Kc…          2697         11629                     3             24463
#>  2 Mfge8_It…          2568          8167                     3             18869
#>  3 Mfge8_It…          2556          9374                     3             20536
#>  4 Tln1_Itg…          2359          6971                     3             15470
#>  5 Fgf12_Fg…          2512          9821                     4             21180
#>  6 L1cam_It…          2345          7629                     5             16831
#>  7 Flt3l_Fl…            78            46                    32                14
#>  8 Pltp_Abc…          1907          3555                    35              7087
#>  9 Cxcl12_G…          1763          2760                    42              5310
#> 10 App_Dcc            1989          3267                    42              8152
#> # ℹ 38 more rows
#> # ℹ 7 more variables: graph_diameter <dbl>, graph_un_diameter <dbl>,
#> #   graph_mean_dist <dbl>, graph_circuit_rank <dbl>, graph_reciprocity <dbl>,
#> #   graph_clique_num <int>, graph_clique_count <int>
```

summarize_ccc_graph_metrics(…, level = “group”) summarizes the metrics
for each subgraph (group) in CCC graph (LR)

``` r
tictoc::tic()

ccc_graph_group_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list, level = "group")

tictoc::toc()
#> 0.06 sec elapsed
```

``` r
ccc_graph_group_metrics_summary_df 
#> # A tibble: 3,330 × 12
#>    LR         group group_n_nodes group_n_edges group_adhesion group_motif_count
#>    <chr>      <int>         <int>         <dbl>          <dbl>             <int>
#>  1 Agtrap_Ra…     4           158           294              0               691
#>  2 Agtrap_Ra…    30            11            16              0                33
#>  3 Agtrap_Ra…    15            24            28              0                76
#>  4 Agtrap_Ra…     8            61            84              0               211
#>  5 Agtrap_Ra…    26            13            14              0                25
#>  6 Agtrap_Ra…     2           172           308              0               742
#>  7 Agtrap_Ra…    22            16            18              0                44
#>  8 Agtrap_Ra…    17            23            36              0                91
#>  9 Agtrap_Ra…     1           184           344              0               846
#> 10 Agtrap_Ra…    36             7             6              0                15
#> # ℹ 3,320 more rows
#> # ℹ 6 more variables: group_diameter <dbl>, group_un_diameter <dbl>,
#> #   group_mean_dist <dbl>, group_girth <dbl>, group_circuit_rank <dbl>,
#> #   group_reciprocity <dbl>
```

## Visualization

``` r
LR_of_interest <- "App_Dcc"
```

### spatial CCC graph plot with tissue image

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       tissue_img = imgRaster(spe),
                       node_color = "group_diameter",
                       node_size = 1,
                       edge_color = "group_diameter",
                       which_on_top = "edge")
```

<img src="man/figures/README-plot-spatial-ccc-graph-tissue-image-1.png" width="100%" />

### spatial CCC graph plot without tissue image

In this case, graph layout can be “manual” which keeps the original
spatial locations, or other graph layout algorithm supported by igraph
package. Below uses “auto” layout (“kk” spring layout)

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       # tissue_img = imgRaster(spe),
                       node_color = "group",
                       node_size = 0.1,
                       edge_color = "group_diameter",
                       edge_width = 0.1,
                       which_on_top = "edge")
```

<img src="man/figures/README-plot-spatial-ccc-graph-kk-layout-1.png" width="100%" />

In this case, below is “stress” layout.

``` r
plot_spatial_ccc_graph(ccc_graph = ccc_graph_list[[LR_of_interest]],
                       # tissue_img = imgRaster(spe),
                       graph_layout = "stress",
                       node_color = "group",
                       edge_color = "group_diameter",
                       edge_width = 0.25,
                       which_on_top = "edge")
```

<img src="man/figures/README-plot-spatial-ccc-graph-stress-layout-1.png" width="100%" />

``` r
tictoc::tic()

cell_overlap_dist <-
  dist_cell_overlap_ccc_tbl(ccc_tbl)

tictoc::toc()
#> 0.342 sec elapsed
```

``` r
tictoc::tic()

cell_overlap_lf <-
  lf_cell_overlap_ccc_tbl(ccc_tbl)

tictoc::toc()
#> 0.306 sec elapsed
```

``` r
tictoc::tic()

LRs_high_cell_overlap <-
  cell_overlap_lf %>% 
  dplyr::filter(d < 1) %>%
  dplyr::select(lr1, lr2) %>%
  unlist() %>% unique()

tictoc::toc()
#> 0.002 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist <-
  cell_overlap_dist[LRs_high_cell_overlap, LRs_high_cell_overlap]

tictoc::toc()
#> 0.001 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist2 <-
  cell_overlap_lf %>%
  dplyr::filter(d < 1) %>%
  dplyr::select(lr1, lr2, d) %>%
  lf_to_dist()

tictoc::toc()
#> 0.006 sec elapsed
```

``` r
LR_ccc_summary_tbl <-
  ccc_tbl %>% 
  pull(LR) %>%
  table() %>%
  as_tibble() %>%
  rename("LR" = ".") %>%
  arrange(desc(n)) %>%
  left_join(
    ccc_tbl %>% 
      select(LR, ligand, receptor) %>%
      distinct(),
    by = "LR"
  )
```

``` r
hclus.res <- fastcluster::hclust(as.dist(high_cell_overlap_dist),
                                 method = "complete")
```

``` r
ape::as.phylo(hclus.res) %>%
  ggtree(layout="dendrogram") %<+% LR_ccc_summary_tbl +
  aes(color=receptor) +
  theme(legend.position = "none") +
  geom_tippoint(aes(size=n), alpha=0.5) +
  geom_tiplab(size=2, offset=-0.15) +
  xlim_tree(3)
```

<img src="man/figures/README-LR-similarity-flat-dendrogram-1.png" width="100%" />

``` r
ape::as.phylo(hclus.res) %>%
  ggtree(layout = "circular") %<+% LR_ccc_summary_tbl +
  aes(color=receptor) +
  # aes(color=ligand) +
  theme(legend.position = "none") +
  geom_tippoint(aes(size=n), alpha=0.5) +
  geom_tiplab(size=2, offset=0.05)
```

<img src="man/figures/README-LR-similarity-circular-dendrogram-1.png" width="100%" />
