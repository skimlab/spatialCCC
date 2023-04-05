
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
#> # A tibble: 103 × 10
#>    LR    ligand_gene_symbol receptor_gene_symbol ligand_gene_id receptor_gene_id
#>    <chr> <chr>              <chr>                         <dbl>            <dbl>
#>  1 Adam… Adam17             Itgb1                         11491            16412
#>  2 Adam… Adam2              Itgb1                         11495            16412
#>  3 Afdn… Afdn               Nectin4                       17356            71740
#>  4 Agt_… Agt                Agtr1a                        11606            11607
#>  5 Ahsg… Ahsg               Insr                          11625            16337
#>  6 Apoc… Apoc1              Vldlr                         11812            22359
#>  7 Apoe… Apoe               Lrp6                          11816            16974
#>  8 Apoe… Apoe               Sdc1                          11816            20969
#>  9 App_… App                Dcc                           11820            13176
#> 10 Bdnf… Bdnf               Sorcs2                        12064            81840
#> # ℹ 93 more rows
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
#> 19.214 sec elapsed

# future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 92,210 × 10
#>    src          dst       d norm.d LR    ligand receptor LRscore weight WLRscore
#>    <chr>        <chr> <dbl>  <dbl> <chr> <chr>  <chr>      <dbl>  <dbl>    <dbl>
#>  1 GTTCATCGTTT… GCGA…  138    1.01 Apoe… Apoe   Lrp6       0.869  0.986    0.857
#>  2 ACTTTACCCTC… GCGA…  138.   1.01 Apoe… Apoe   Lrp6       0.868  0.980    0.851
#>  3 AGATGATGGAG… GCGA…  138.   1.01 Apoe… Apoe   Lrp6       0.867  0.980    0.849
#>  4 CTCTGCAGGCA… CTTC…  137    1    Npy_… Npy    Npy2r      0.858  1        0.858
#>  5 ACACATGATCA… GCTG…  138    1.01 Apoe… Apoe   Lrp6       0.852  0.986    0.840
#>  6 TATTCAATTCT… GGGA…  138.   1.01 Apoe… Apoe   Lrp6       0.851  0.980    0.834
#>  7 GGAGCGAGGCC… GAGG…  137    1    Apoe… Apoe   Lrp6       0.851  1        0.851
#>  8 TGACGAATATT… TCGT…  138.   1.01 Apoe… Apoe   Lrp6       0.850  0.980    0.833
#>  9 GTCCCAACGTA… GGGA…  138    1.01 Apoe… Apoe   Lrp6       0.850  0.986    0.837
#> 10 GGTAAATGTGC… TAGT…  138.   1.01 Vip_… Vip    Ramp1      0.850  0.987    0.838
#> # ℹ 92,200 more rows
```

### Convert CCC table to CCC graph

The conversion also adds various graph metrics to each CCC graph.

``` r
tictoc::tic("to_spatical_ccc_graph ...")

sp_col_data <- get_spatial_coords(spe)

ccc_graph_list <-
  to_spatial_ccc_graph_list(ccc_tbl, sp_col_data, workers = 4)

tictoc::toc()
#> to_spatical_ccc_graph ...: 38.355 sec elapsed
```

summarize_ccc_graph_metrics() summarize those graph metrics for each LR
pair.

``` r
tictoc::tic()

ccc_graph_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list)

tictoc::toc()
#> 0.183 sec elapsed
```

``` r
ccc_graph_metrics_summary_df %>%
  arrange(graph_component_count)
#> # A tibble: 49 × 12
#>    LR        graph_n_nodes graph_n_edges graph_component_count graph_motif_count
#>    <chr>             <int>         <dbl>                 <dbl>             <int>
#>  1 S100b_Fg…          2578          8316                     3             18789
#>  2 Fgf13_Fg…          2423          7295                     6             16329
#>  3 Gas6_Tyr…          2372          7658                     7             16934
#>  4 Apoe_Lrp6          2497          5640                     9             13761
#>  5 Hsp90b1_…          2461          5543                     9             13394
#>  6 Ptn_Ptprb          2381          4685                    14             11391
#>  7 Vip_Ramp1          1761          4214                    29              9757
#>  8 Cd14_C3a…            83            57                    33                22
#>  9 Efna3_Ep…          1591          2637                    37              5314
#> 10 App_Dcc            1989          3267                    42              8152
#> # ℹ 39 more rows
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
#> 0.151 sec elapsed
```

``` r
ccc_graph_group_metrics_summary_df 
#> # A tibble: 3,720 × 12
#>    LR       group group_n_nodes group_n_edges group_adhesion group_motif_count
#>    <chr>    <int>         <int>         <dbl>          <dbl>             <int>
#>  1 Vim_Cd44     8            21            26              0                55
#>  2 Vim_Cd44    79             2             1              0                 0
#>  3 Vim_Cd44     1            45            80              0               150
#>  4 Vim_Cd44     6            27            37              0                65
#>  5 Vim_Cd44     9            16            32              0                66
#>  6 Vim_Cd44    27             6             5              0                10
#>  7 Vim_Cd44    45             4             3              0                 3
#>  8 Vim_Cd44    16             9             8              0                13
#>  9 Vim_Cd44     4            30            44              0               109
#> 10 Vim_Cd44    17             9            14              0                23
#> # ℹ 3,710 more rows
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

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

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

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

``` r
tictoc::tic()

cell_overlap_dist <-
  dist_cell_overlap_ccc_tbl(ccc_tbl)

tictoc::toc()
#> 1.107 sec elapsed
```

``` r
tictoc::tic()

cell_overlap_lf <-
  lf_cell_overlap_ccc_tbl(ccc_tbl)

tictoc::toc()
#> 1.149 sec elapsed
```

``` r
tictoc::tic()

LRs_high_cell_overlap <-
  cell_overlap_lf %>% 
  filter(d < 1) %>%
  select(lr1, lr2) %>%
  unlist() %>% unique()

tictoc::toc()
#> 0.007 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist <-
  cell_overlap_dist[LRs_high_cell_overlap, LRs_high_cell_overlap]

tictoc::toc()
#> 0.003 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist2 <-
  cell_overlap_lf %>%
  filter(d < 1) %>%
  select(lr1, lr2, d) %>%
  lf_to_dist()

tictoc::toc()
#> 0.02 sec elapsed
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

<img src="man/figures/README-unnamed-chunk-24-1.png" width="100%" />

``` r
ape::as.phylo(hclus.res) %>%
  ggtree(layout = "circular") %<+% LR_ccc_summary_tbl +
  aes(color=receptor) +
  # aes(color=ligand) +
  theme(legend.position = "none") +
  geom_tippoint(aes(size=n), alpha=0.5) +
  geom_tiplab(size=2, offset=0.05)
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />

``` r
gp_hex <-
  ccc_graph_list[[LR_of_interest]] %>%
  plot_spatial_ccc_graph(
    tissue_img = imgRaster(spe),
    clip = TRUE,
    node_color = "group",
    node_size = 1,
    edge_color = "group_n_edges",
    edge_width = 0.5,
    ghost_img = TRUE,
    which_on_top = "edge"
  ) +
  theme(legend.position = "none")

gp_hex
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="100%" />

``` r
library(hexSticker)
gp_hex_sticker <- sticker(
  gp_hex,
  package = "spatialCCC",
  # p_color = "#1F1F1F",
  p_color = "snow",
  p_family = "sans",
  p_fontface = "bold",
  p_y = 1.1,
  p_size = 40,
  
  s_x = 1.0,
  s_y = 1.2,
  s_width = 5,
  s_height = 5*3/4,
  
  # h_fill = "purple",
  h_color = "snow",
  h_size = 1.5,
  
  # spotlight = TRUE,
  # l_x = 1.0,
  # l_y = 0.4,
  # l_alpha = 0.25,
  
  filename = "spatialCCC.png",
  dpi = 600
)
```
