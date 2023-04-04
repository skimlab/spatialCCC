
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialCCC <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/dolchan/spatialCCC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dolchan/spatialCCC/actions/workflows/R-CMD-check.yaml)
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
#> # A tibble: 105 × 10
#>    LR    ligand_gene_symbol receptor_gene_symbol ligand_gene_id receptor_gene_id
#>    <chr> <chr>              <chr>                         <dbl>            <dbl>
#>  1 Adam… Adam10             Epha3                         11487            13837
#>  2 Adam… Adam10             Tspan10                       11487           208634
#>  3 Adam… Adam15             Itgb1                         11490            16412
#>  4 Adam… Adam17             Ghr                           11491            14600
#>  5 Adam… Adam17             Itgb6                         11491            16420
#>  6 Adam… Adam2              Cd9                           11495            12527
#>  7 Adcy… Adcyap1            Vipr2                         11516            22355
#>  8 Anxa… Anxa1              Dysf                          16952            26903
#>  9 Apel… Apela              Aplnr                     100038489            23796
#> 10 Apob… Apob               Apobr                        238055           171504
#> # ℹ 95 more rows
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
#> 18.206 sec elapsed

# future::plan(future::sequential)
```

``` r
ccc_tbl %>% 
  dplyr::arrange(desc(LRscore))
#> # A tibble: 138,939 × 10
#>    src          dst       d norm.d LR    ligand receptor LRscore weight WLRscore
#>    <chr>        <chr> <dbl>  <dbl> <chr> <chr>  <chr>      <dbl>  <dbl>    <dbl>
#>  1 GTTATTAACGT… TTGC…  138    1.01 Apoe… Apoe   Abca1      0.862  0.986    0.850
#>  2 AGAAGTGATTC… AATG…  138    1.01 Calm… Calm1  Kcnq3      0.858  0.986    0.846
#>  3 CGTGCTGGCCT… GGGT…  137    1    Hsp9… Hsp90… Fgfr3      0.858  1        0.858
#>  4 TTGCATGCTGA… GTGC…  138.   1.01 Calm… Calm1  Kcnq3      0.858  0.980    0.840
#>  5 CCCAACATACG… GTGC…  138.   1.00 Calm… Calm1  Kcnq3      0.857  0.992    0.850
#>  6 CGGTGCAGATA… GTGC…  138    1.01 Calm… Calm1  Kcnq3      0.857  0.986    0.845
#>  7 GTGGTTTCCGC… GTGC…  138    1.01 Calm… Calm1  Kcnq3      0.857  0.986    0.845
#>  8 TCTGGGTAGCG… ACAT…  138.   1.01 Calm… Calm1  Kcnq3      0.857  0.987    0.845
#>  9 TGTCGGCATGG… GGGT…  138    1.01 Hsp9… Hsp90… Fgfr3      0.856  0.986    0.844
#> 10 AGTTTGGCACG… GCAT…  137.   1.00 Hsp9… Hsp90… Fgfr3      0.856  0.999    0.856
#> # ℹ 138,929 more rows
```

### Convert CCC table to CCC graph

The conversion also adds various graph metrics to each CCC graph.

``` r
tictoc::tic("to_spatical_ccc_graph ...")

sp_col_data <- get_spatial_coords(spe)

ccc_graph_list <-
  to_spatial_ccc_graph_list(ccc_tbl, sp_col_data, workers = 4)

tictoc::toc()
#> to_spatical_ccc_graph ...: 40.503 sec elapsed
```

summarize_ccc_graph_metrics() summarize those graph metrics for each LR
pair.

``` r
tictoc::tic()

ccc_graph_metrics_summary_df <-
  summarize_ccc_graph_metrics(ccc_graph_list)

tictoc::toc()
#> 0.191 sec elapsed
```

``` r
ccc_graph_metrics_summary_df %>%
  arrange(graph_component_count)
#> # A tibble: 53 × 12
#>    LR        graph_n_nodes graph_n_edges graph_component_count graph_motif_count
#>    <chr>             <int>         <dbl>                 <dbl>             <int>
#>  1 Mfge8_It…          2568          8167                     3             18869
#>  2 Calm1_Kc…          2697         11629                     3             24463
#>  3 Hsp90aa1…          2695         12051                     3             25363
#>  4 Nlgn2_Nr…          2636         11910                     4             24385
#>  5 Adam15_I…          2344          6849                     4             15221
#>  6 Apoe_Abc…          2548          5733                     5             14085
#>  7 App_Lrp8           2524          6406                     6             15469
#>  8 Calm3_Ad…          2396          5839                     7             14043
#>  9 Efna3_Ep…          2361          8104                     8             17563
#> 10 Vegfa_Fl…          1972          3998                    13              8407
#> # ℹ 43 more rows
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
#> 0.169 sec elapsed
```

``` r
ccc_graph_group_metrics_summary_df 
#> # A tibble: 3,341 × 12
#>    LR       group group_n_nodes group_n_edges group_adhesion group_motif_count
#>    <chr>    <int>         <int>         <dbl>          <dbl>             <int>
#>  1 Bsg_Cav1     1          1913          3392              0              8432
#>  2 Bsg_Cav1    19             5             4              0                 6
#>  3 Bsg_Cav1     2            63           114              0               276
#>  4 Bsg_Cav1    14             7             6              0                15
#>  5 Bsg_Cav1     5            23            24              0                65
#>  6 Bsg_Cav1     7            17            18              0                43
#>  7 Bsg_Cav1    12             8             9              0                15
#>  8 Bsg_Cav1     4            24            36              0                87
#>  9 Bsg_Cav1    11            11            13              0                24
#> 10 Bsg_Cav1     9            12            19              0                35
#> # ℹ 3,331 more rows
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
#> 1.035 sec elapsed
```

``` r
tictoc::tic()

cell_overlap_lf <-
  lf_cell_overlap_ccc_tbl(ccc_tbl)

tictoc::toc()
#> 0.924 sec elapsed
```

``` r
tictoc::tic()

LRs_high_cell_overlap <-
  cell_overlap_lf %>% 
  filter(d < 1) %>%
  select(lr1, lr2) %>%
  unlist() %>% unique()

tictoc::toc()
#> 0.005 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist <-
  cell_overlap_dist[LRs_high_cell_overlap, LRs_high_cell_overlap]

tictoc::toc()
#> 0.002 sec elapsed
```

``` r
tictoc::tic()

high_cell_overlap_dist2 <-
  cell_overlap_lf %>%
  filter(d < 1) %>%
  select(lr1, lr2, d) %>%
  lf_to_dist()

tictoc::toc()
#> 0.017 sec elapsed
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
    edge_color = "group",
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
  p_size = 20,
  s_x = 1.5,
  s_y = 0.5,
  s_width = 4,
  s_height = 4,
  filename = "spatialCCC.png"
)
```
