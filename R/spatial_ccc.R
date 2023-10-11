#' Compute spatial ligand-receptor interactions
#'
#' Compute ligand-receptor interactions (LRscore) between neighboring spots.
#'   Neighboring spots are those with norm.d < spot_dist_cutoff.
#'
#'   \deqn{LRscore = \frac{\sqrt{x_L \cdot x_R}}{\sqrt{x_L \cdot x_R} + C}}
#'
#'   where \eqn{x_L} and \eqn{x_R} are expression values of
#'   ligand (L) and receptor (R), and \eqn{C} is the average gene expression
#'   across all genes and all samples.
#'
#'   source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7261168/
#'
#' @param spe SpatialExperiment object
#' @param assay_name assay name in string
#' @param ligand ligand
#' @param receptor receptor
#' @param expression_min_prop minimum proportion of samples
#'    with non-zero expression value (default: 0.05)
#' @param spot_dist distances between spots, an output of [calc_spot_dist()]
#' @param spot_dist_cutoff cutoff value for norm.d in spot
#'    distances. Default cutoff is 1.5 and
#'    it comes from sqrt(3) ~ 1.73, based on how Visium spot-array is
#'    arranged, so by default, only computing nearest neighbors.
#' @param LRscore_cutoff minimum LRscore to keep
#'
#'
#' @returns LRscore table in data.frame with the following columns:
#' \describe{
#'   \item{src, dst}{source and destination cells/spots}
#'   \item{d, norm.d}{cell-cell/spot-spot distance and normalized distance}
#'   \item{ligand, receptor}{ligand & receptor genes}
#'   \item{LRscore}{LRscore -- see above}
#'   \item{weight}{weight -- reciprocal of norm.d}
#'   \item{WLRscore}{weighted LRscore = \eqn{weight \cdot LRscore}}
#' }
#'
#' @export
compute_spatial_ccc_tbl <-
  function(spe,
           assay_name = "logcounts",
           ligand, receptor,
           expression_min_prop = 0.05,
           spot_dist = calc_spot_dist(spe),
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5) {
    #
    # spot_dist <- calc_spot_dist(spe)

    gexp <-
      as.matrix(SummarizedExperiment::assays(spe)[[assay_name]])
    rowAnnots <- SingleCellExperiment::rowData(spe)

    # first remove genes with no expression at all
    #   maybe, we need to remove those with expression in only 5% (cutoff) or less
    n_cutoff <- ncol(gexp) * expression_min_prop
    subset_idx <- rowSums(gexp > 0) > n_cutoff
    gexp <- gexp[subset_idx, ]
    rowAnnots <- subset(rowAnnots, subset = subset_idx)

    # background/baseline expression
    c_mean <- mean(gexp)

    compute_LRscore <- function(lig, rec) {

      # ligand and receptor expression values
      lig_ids <- which(rowAnnots$symbol == lig)
      rec_ids <- which(rowAnnots$symbol == rec)

      if (length(lig_ids) > 1) {
        lig_exp <- colMeans(gexp[lig_ids, ])
      } else {
        lig_exp <- gexp[lig_ids, ]
      }

      if (length(rec_ids) > 1) {
        rec_exp <- colMeans(gexp[rec_ids, ])
      } else {
        rec_exp <- gexp[rec_ids, ]
      }

      sqrt.prod <-
        sqrt(lig_exp[spot_dist$src] * rec_exp[spot_dist$dst])

      tibble(spot_dist, LRscore = sqrt.prod / (sqrt.prod + c_mean)) %>%
        # Weighted LR score to account for
        #   the attenuation of cell signaling due to
        #   traveling distance.
        #
        # when norm.d = 0, we will use 1 to compute weight because
        #   1. this is not still a single cell level, hence, no guarantee it's
        #      autocrine, but still very close to each other.
        dplyr::mutate(weight = ifelse(norm.d == 0, 1, 1 / norm.d ^ 2)) %>%
        dplyr::mutate(WLRscore = LRscore * weight) %>%

        # temporarily disable the line below to calculate "p-value"?
        # now back to LRscore filtering
        dplyr::filter(LRscore > LRscore_cutoff) %>%

        dplyr::mutate(LRscore_perc_rank = percent_rank(LRscore)) %>%

        # keep as many as possible, except LRscore = 0
        # dplyr::filter(LRscore > 0) %>%

        dplyr::mutate(ligand = lig, receptor = rec)
    }

    compute_LRscore(lig = ligand, rec = receptor) %>%
      dplyr::relocate(ligand, receptor, .after = "norm.d")
  }


#' Compute a spatial CCC graph
#'
#' @inheritParams compute_spatial_ccc_tbl
#' @inheritParams to_spatial_ccc_graph
#'
#' @returns a spatial CCC graphs in `tidygraph` format;
#'   spatial CCC graph also include the following graph metrics
#'
#' @section graph metrics:
#'   1. For overall graph,
#'      * graph_n_nodes,
#'      * graph_n_edges
#'      * graph_component_count
#'      * graph_motif_count
#'      * graph_diameter
#'      * graph_un_diameter
#'      * graph_mean_dist
#'      * graph_circuit_rank = graph_n_edges - graph_n_nodes + graph_component_count
#'      * graph_reciprocity
#'      * graph_clique_num (sp_ccc_graph assumed as undirected)
#'      * graph_clique_count (sp_ccc_graph assumed as undirected)
#'
#'  2. For each sub-graph (after [group_components()])
#'      * group_n_nodes
#'      * group_n_edges
#'      * group_adhesion
#'      * group_motif_count
#'      * group_diameter
#'      * group_un_diameter
#'      * group_mean_dist
#'      * group_girth
#'      * group_circuit_rank
#'      * group_reciprocity
#'
#' @export
compute_spatial_ccc_graph <-
  function(spe,
           assay_name = "logcounts",
           ligand,
           receptor,
           expression_min_prop = 0.05,
           spot_dist = calc_spot_dist(spe),
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5) {
    compute_spatial_ccc_tbl(
      spe,
      assay_name,
      ligand,
      receptor,
      expression_min_prop,
      spot_dist,
      spot_dist_cutoff,
      LRscore_cutoff
    ) %>%
      to_spatial_ccc_graph(get_spatial_data(spe))
  }


#' Calculate distance between spots
#'
#' @inheritParams compute_spatial_ccc_tbl
#' @param long_format if TRUE (default), return long format data.frame.
#'   Otherwise, matrix format.
#'
#' @returns A data fame with computed distance including both `d` and `norm.d`.
#'   `norm.d` is a multiple of the shortest distance between spots.
#'
#' @export
calc_spot_dist <-
  function(spe,
           long_format = TRUE,
           spot_dist_cutoff = 1.5) {
    spe_cd <- get_spatial_data(spe)

    spot_dist <-
      calc_spot_dist0(spe_cd[c("cell_id", "pxl_col_in_fullres", "pxl_row_in_fullres")],
                      long_format = long_format)

    spot_dist %>%
      dplyr::filter(norm.d < spot_dist_cutoff)
  }


#' Compute spatial CCC graphs over a list of LR pairs
#'
#' @inheritParams compute_spatial_ccc_graph
#' @param LRdb LRdb in data.frame (see [get_LRdb()])
#' @param workers the number of processes to be used for parallel
#'   processing
#'
#' @returns a list of spatial CCC graphs
#'
#' @inheritSection compute_spatial_ccc_graph graph metrics
#'
#' @export
compute_spatial_ccc_graph_list <-
  function(spe,
           assay_name = "logcounts",
           LRdb,
           expression_min_prop = 0.05,
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5,
           workers = 1) {

    # will save some overhead in compute_spatial_ccc_graph
    spot_dist <- calc_spot_dist(spe)

    # sp_col_data <-
    #   get_spatial_data(spe)

    # probably redundant, but needed to avoid unnecessary calculation later
    gexp <-
      as.matrix(SummarizedExperiment::assays(spe)[[assay_name]])
    rowAnnots <- SingleCellExperiment::rowData(spe)

    # first remove genes with no expression at all
    #   maybe, we need to remove those with expression in only 5% (cutoff) or less
    n_cutoff <- ncol(gexp) * expression_min_prop
    subset_idx <- rowSums(gexp > 0) > n_cutoff
    gexp <- gexp[subset_idx,]
    rowAnnots <- subset(rowAnnots, subset = subset_idx)

    # make sure to process only LR pairs present in gene expression data
    LRdb <-
      LRdb %>%
      dplyr::filter(
        ligand_gene_symbol %in% rowAnnots$symbol,
        receptor_gene_symbol %in% rowAnnots$symbol
      )

    Ks <- setNames(1:nrow(LRdb), nm = LRdb$LR)

    if (workers > 1) {
      # set up future::plan and its release on exit
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential), add = TRUE)

      furrr::future_map(Ks,
                        function(k) {
                          compute_spatial_ccc_graph(
                            spe,
                            assay_name,
                            ligand = LRdb$ligand_gene_symbol[k],
                            receptor = LRdb$receptor_gene_symbol[k],
                            expression_min_prop,
                            spot_dist,
                            spot_dist_cutoff,
                            LRscore_cutoff
                          )
                        },
                        # making sure "dplyr" package is attached to future parallel processes
                        .options = furrr::furrr_options(seed = TRUE,
                                                        packages = c("dplyr")))
    } else {
      purrr::map(Ks,
                 function(k) {
                   compute_spatial_ccc_graph(
                     spe,
                     assay_name,
                     ligand = LRdb$ligand_gene_symbol[k],
                     receptor = LRdb$receptor_gene_symbol[k],
                     expression_min_prop,
                     spot_dist,
                     spot_dist_cutoff,
                     LRscore_cutoff
                   )
                 })
    }
  }


#' Flatten multiple CCC graphs into aggregate one
#'
#' @param ccc_graph_ls a list of spatial CCC graphs
#'
#' @returns a spatial CCC graph
#'
#' @section flattening_CCC_graph:
#'
#' The following node parameters will be kept:
#'   * name
#'   * spot_x, spot_y
#'   * in_tissue
#'   * array_row, array_col
#'   * sizeFactor
#'
#' The following node parameters will be aggregated over all LR pairs:
#'   * n, n.src, n.dst, LRscore.sum.src, and LRscore.sum.dst will be summed.
#'   * node_is_src and node_is_dst will be OR'd
#'   * n.inflow and LRscore.sum.inflow will be summed.
#'
#' @export
flatten_ccc_graph_list <- function(ccc_graph_ls) {
  ccog_nodes <-
    ccc_graph_ls %>%
    purrr::map(function(ccog) {
      ccog %>%
        tidygraph::activate("nodes") %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(name,
             spot_x,
             spot_y,
             in_tissue,
             array_row,
             array_col,
             sizeFactor) %>%
    summarise(
      n = n(),
      n.src = sum(n.src),
      LRscore.sum.src = sum(LRscore.sum.src),
      # WLRscore.sum.src = sum(WLRscore.sum.src),
      n.dst = sum(n.dst),
      LRscore.sum.dst = sum(LRscore.sum.dst),
      # WLRscore.sum.dst = sum(WLRscore.sum.dst),
      node_is_src = any(node_is_src),
      node_is_dst = any(node_is_dst),
      n.inflow = sum(n.inflow),
      LRscore.sum.inflow = sum(LRscore.sum.inflow),
      # WLRscore.sum.inflow = sum(WLRscore.sum.inflow),
      .groups = "drop"
    ) %>%
    rename(cell_id = name)

  ccog_edges <-
    ccc_graph_ls %>%
    purrr::map(function(ccog) {
      ccog %>%
        tidygraph::activate("edges") %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(src, dst, d, norm.d) %>%
    summarise(n = n(),
              LRscore = max(LRscore),
              # WLRscore = max(WLRscore),
              LRscore.sum = sum(LRscore),
              # WLRscore.sum = sum(WLRscore),
              LRscore.sd = sd(LRscore),
              # WLRscore.sd = sd(WLRscore),
              .groups = "drop")

  ccog_edges %>%
    to_barebone_spatial_ccc_graph() %>%
    tidygraph::activate(nodes) %>%
    left_join(ccog_nodes,
              by = c("name" = "cell_id"))
}



#' Summarize spatial CCC graph for LR
#'
#' @inheritParams summarize_ccc_graph_metrics
#'
#' @returns spatial CCC graph summarized for LR,
#'    aggregating LRscores over the graph and computing percentile ranks.
#'
#' @inheritSection calculate_LR_percent_rank CCC_graph_percent_rank
#' @export
to_LR_summary_tbl <- function(ccc_graph_ls) {
  left_join(
    calculate_LR_percent_rank(ccc_graph_ls),
    summarize_ccc_graph_metrics(ccc_graph_ls, level = "graph"),
    by = "LR"
  ) %>%
    split_ccc_tbl_LR()
}


#' Convert spatial CCC graphs to spatial CCC table
#'
#' @param ccc_graph_ls a list of spatial CCC graphs
#'
#' @returns spatial CCC table
#'
#' @export
extract_spatial_ccc_graph_edges <-
  function(ccc_graph_ls) {
    lapply(ccc_graph_ls,
           function(ccog) {
             ccog %>%
               tidygraph::activate("edges") %>%
               tibble::as_tibble() %>%
               dplyr::select(-c("from", "to"))
           }) %>%
      dplyr::bind_rows(.id = "LR")
  }





#' Summarize spatial CCC graph list to graph/group metrics table
#'
#' @param ccc_graph_ls list of spatial CCC graph,
#'   each of which is an output of to_spatial_ccc_graph.
#' @param level extract graph metrics from either overall graph ("graph") or
#'   subgraph ("group")
#'
#' @returns graph metrics table summarized for each LR pair
#'
#' @export
summarize_ccc_graph_metrics <- function(ccc_graph_ls,
                                        level = c("graph", "group")) {
  ccc_graph_ls %>%
    purrr::map(function(ccc_graph) {
      ccc_graph %>%
        extract_ccc_graph_metrics(level)
    }) %>%
    dplyr::bind_rows(.id = "LR")
}


#' Amend CCC table with annotations for source and destination cells
#'
#' @param ccc_tbl CCC table, an output of [compute_spatial_ccc_tbl()]
#' @param spe SpatialExperiment object
#' @param annot_cols column names in colData(spe), which will be matched
#'   against source and destination cells and added separately:
#'   as in annot_cols{.src, .dst}
#' @param overwrite If TRUE, the previous columns colliding with
#'   the new columns (annot_cols.{src, dst}) will be reset.
#'
#' @export
amend_ccc_tbl_with_cell_annots <-
  function(ccc_tbl, spe, annot_cols, overwrite = TRUE) {

  # temp internal functions
  add_dot_src <- function(s) {
    paste0(s, ".src")
  }

  add_dot_dst <- function(s) {
    paste0(s, ".dst")
  }

  if (overwrite) {
    annot_cols.after <-
      lapply(annot_cols,
             function(s) {
               c(add_dot_src(s), add_dot_dst(s))
             }) %>% unlist()

    ccc_tbl <-
      ccc_tbl %>%
      dplyr::select(-dplyr::any_of(annot_cols.after))
  }

  annot_df <-
    SingleCellExperiment::colData(spe) %>%
    tibble::as_tibble(rownames = "cell_id") %>%
    dplyr::select(cell_id, dplyr::all_of(annot_cols))

  ccc_tbl %>%
    # add annotations to source cells
    dplyr::left_join(annot_df, by = c("src" = "cell_id")) %>%
    dplyr::rename_with(add_dot_src, dplyr::all_of(annot_cols)) %>%

    # add annotations to destination cells
    dplyr::left_join(annot_df, by = c("dst" = "cell_id")) %>%
    dplyr::rename_with(add_dot_dst, dplyr::all_of(annot_cols))
}


#' Set default cluster in CCC table
#'
#' @param ccc_tbl CCC table, output of [compute_spatial_ccc_tbl()]
#' @param cluster_name string, the name of default cluster.
#'   {cluster_name}.src and {cluster_name}.dst columns will be copied to
#'   "cluster.src" and "cluster.dst".
#'
#' @returns CCC table with amended columns
#'
#' @export
set_ccc_tbl_default_cluster <-
  function(ccc_tbl, cluster_name) {
    ccc_tbl %>%
      dplyr::mutate(cluster.src = get(paste0(cluster_name, ".src")),
                    cluster.dst = get(paste0(cluster_name, ".dst")))
  }

#' Tidy up CCC graph by removing isolated nodes
#'
#' @param ccog CCC graph
#'
#' @returns CCC graph with all isolated nodes removed.
#'
#' @export
tidy_up_ccc_graph <-
  function(ccog) {
    connected_nodes <-
      ccog %>%
      tidygraph::activate("edges") %>%
      tibble::as_tibble() %>%
      dplyr::select(src, dst) %>%
      unlist()

    ccog %>%
      tidygraph::activate("nodes") %>%
      tidygraph::filter(name %in% connected_nodes)
  }


#' Tag nodes if included in COIs
#'
#' Add additional columns to all nodes in CCC graph:
#'   1. InCOI = TRUE if a node is included in COIs or is connected
#'                  to a node in COIs
#'   2. tagged = TRUE if a node is included in a group that includes
#'                   any of COIs, when edges_expanded_to_group = TRUE
#'
#' @param ccog CCC graph
#' @param COIs array of nodes to highlight
#' @param edges_expanded_to_group If TRUE,
#'
#' @returns CCC graph with InCOI and tagged columns (see above)
tag_cells_in_ccc_graph <- function(ccog, COIs, edges_expanded_to_group = FALSE) {
  ccog <-
    ccog %>%
    tidygraph::activate("nodes") %>%
    tidygraph::mutate(InCOI = name %in% COIs) %>%
    tidygraph::mutate(tagged = InCOI) %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(InCOI = src %in% COIs | dst %in% COIs) %>%
    tidygraph::mutate(tagged = InCOI) # by default

  # now including all cells that belong to
  #   the groups that cells_of_interest belong to
  if (edges_expanded_to_group) {
    GOIs <-
      ccog %>%
      tidygraph::activate("nodes") %>%
      tidygraph::filter(InCOI) %>%
      tidygraph::pull(group) %>%
      unique()

    cells_in_GOI <-
      ccog %>%
      tidygraph::activate("nodes") %>%
      tidygraph::filter(group %in% GOIs) %>%
      tidygraph::pull(name) %>%
      unique()

    ccog <-
      ccog %>%
      tidygraph::activate("nodes") %>%
      tidygraph::mutate(tagged = name %in% cells_in_GOI) %>%
      tidygraph::activate("edges") %>%
      tidygraph::mutate(tagged = src %in% cells_in_GOI |
                          dst %in% cells_in_GOI)
  }

  ccog
}

#
# Internal functions =====
#





#' Convert CCC table to CCC graph
#'
#'
#' @param ccc_tbl an output of [compute_spatial_ccc_tbl()]
#' @param spatial_data a table with spatial coordinates of
#'   spots/cells in spatial transcriptomic data.  if spatial_data
#'   includes cell_id, spot_x, spot_y, they will be used as is;
#'   Otherwise, use the first three columns as if they are
#'   cell_id, spot_x, spot_y
#'
#' @returns spatial CCC graph -- a `tidygraph`
#'
#' @inheritSection compute_spatial_ccc_graph_list graph metrics
#'
#' @export
to_spatial_ccc_graph <-
  function(ccc_tbl,
           spatial_data) {

    # check if we have these three columns
    cell_spot_columns <-
      match(c("cell_id", "spot_x", "spot_y"),
            colnames(spatial_data))

    if (any(is.na(cell_spot_columns))) {
      cli::cli_abort(message = "{.var spatial_data} should include: 'cell_id', 'spot_x', 'spot_y'")
    } else {
      spatial_data <-
        dplyr::relocate(spatial_data,
                        cell_id, spot_x, spot_y)
    }

    ccc_by_cell <-
      summarise_ccc_by_cell(ccc_tbl) %>%
      add_spatial_data(spatial_data)

    ccc_tbl %>%
      to_barebone_spatial_ccc_graph() %>%
      tidygraph::activate("nodes") %>%
      dplyr::left_join(ccc_by_cell,
                       by = c("name" = "cell_id")) %>%
      #
      # add various graph metrics to ccc_graph
      #
      add_spatial_ccc_graph_metrics()
  }


#' Calculate (combined) percent rank of n and LRscore
#'
#' @param ccc_graph_ls list of spatial CCC graph
#'
#' @returns data frame (tibble)
#'
#' @section CCC_graph_percent_rank:
#'
#' In the table returned,
#'   * n:  number of edges in CCC graph (LR),
#'   * LRscore: median of all LRscores in CCC graph
#'   * n_perc_rank:  percent rank of n
#'   * LRscore_perc_rank: percent rank of LRscore (median)
#'   * perc_rank: combined percent rank (geometric mean of n_perc_rank and LRscore_perc_rank)
#'   * perc_group: binned perc_rank
#'
calculate_LR_percent_rank <- function(ccc_graph_ls) {
  ccc_graph_el_tbl <-
    ccc_graph_ls %>%
    extract_spatial_ccc_graph_edges()

  ccc_graph_el_tbl %>%
    group_by(LR) %>%
    summarise(n = n(),
              LRscore = median(LRscore)) %>%
    mutate(n_perc_rank = percent_rank(n),
           LRscore_perc_rank = percent_rank(LRscore)) %>%
    # geometric mean of percent ranks of n and LRscore
    mutate(perc_rank = sqrt(n_perc_rank * LRscore_perc_rank)) %>%
    mutate(perc_group = sprintf("perc_%02d", floor(perc_rank * 10) * 10))
}


#' Calculate distance between spots
#'
#' @inheritParams calc_spot_dist
#' @param coords spatial coordinates in data.frame
#'   with `cell_id`, `x`, and `y` in that order.  The column headings
#'   do not matter, just the order.
#'
#' @returns A data fame with computed distance including both `d` and `norm.d`.
#'   `norm.d` is a multiple of the shortest distance between spots.
#'
#' @noRd
calc_spot_dist0 <- function(coords, long_format = TRUE) {
  cell_ids <- coords[[1]]

  x <-
    coords[-1] %>%
    dist()

  x <- as.matrix(x)
  x_min <- min(x[x > 0])

  colnames(x) <- cell_ids
  rownames(x) <- cell_ids

  if (long_format) {
    x %>%
      tibble::as_tibble(rownames = "src") %>%
      tidyr::pivot_longer(cols = -src,
                          names_to = "dst",
                          values_to = "d") %>%
      # initially not included,
      #   but now included because spots are not single cells
      #   maybe, later when spatial transcriptomics become
      #   single cell resolution, might need to separate this
      #   into auocrine and paracrine
      # dplyr::filter(d > 0) %>%
      dplyr::mutate(norm.d = d / x_min)
  } else {
    x
  }
}


#' Get spatial coordinates
#'
#' Returns default spatial coordinates and scaled coordinates,
#'   `spot_x` and `spot_y`.
#'
#' @param spe SpatialExperiment object
#'
#' @returns spatial coordinates in data.frame with
#'   `cell_id`, ... (from `spatialCoords(spe)`),
#'   `spot_x`, `spot_y`.
#'
#' @export
get_spatial_data <- function(spe) {
  spe_col_data <- SingleCellExperiment::colData(spe)
  spe_spatial_coords <- SpatialExperiment::spatialCoords(spe)

  # for some reasons, as_tibble() convert column names with checknames = TRUE
  #   so, need to put them back
  cnames_spe_col_data <- colnames(spe_col_data)
  spe_col_data <- tibble::as_tibble(spe_col_data, rownames = "cell_id")
  colnames(spe_col_data) <- c("cell_id", cnames_spe_col_data)

  cnames_spe_spatial_coords <- colnames(spe_spatial_coords)
  spe_spatial_coords <- tibble::as_tibble(spe_spatial_coords, rownames = "cell_id")
  colnames(spe_spatial_coords) <- c("cell_id", cnames_spe_spatial_coords)

  dplyr::left_join(
    spe_col_data,
    spe_spatial_coords,
    by = "cell_id"
  ) %>%
    dplyr::mutate(
      spot_x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe),
      spot_y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe)
    )
}

#' Split LR into L and R
#'
#' @param ccc_tbl an output of [extract_spatial_ccc_graph_edges()] or
#'                any table with "LR" column which should be in "$L_$R" format
split_ccc_tbl_LR <- function(ccc_tbl) {
  ccc_tbl %>%
    mutate(L = sub("_.+", "", LR),
           R = sub(".+_", "", LR)) %>%
    relocate(L, R, .after = "LR")
}

#' Create barebone CCC graph
#'
#' @inheritParams compute_spatial_ccc
#'
#' @returns spatial CCC graph w/o LRscore
#'
#' @noRd
barebone_spatial_ccc_graph <-
  function(spe,
           spot_dist_cutoff = 1.5) {
    calc_spot_dist(spe) %>%
      to_barebone_spatial_ccc_graph()
  }


#' Add spatial data to CCC table
#'
#' internal function
#' @param df CCC table. see [compute_spatial_ccc()].
#' @param spatial_data spatial data: cell_id, spot_x, spot_y, ...;
#'   they should contain cell_id, spot_x, spot_y, and the rest will
#'   be added
#'
#' @noRd
add_spatial_data <-
  function(df,
           spatial_data) {

    # adjust colnames to make them consistent
    spatial_data <-
      spatial_data %>%
      dplyr::relocate(cell_id, spot_x, spot_y)

    dplyr::left_join(df,
                     spatial_data,
                     by = "cell_id")
  }



#' Add various graph metrics to spatial CCC graph
#'
#' See below and [tidygraph::graph_measures] for how these metrics are computed.
#'
#' @param sp_ccc_graph an output of [to_spatial_ccc_graph()]
#' @param from_scratch if TRUE, the existing metrics are wiped clean.
#'
#' @noRd
#'
add_spatial_ccc_graph_metrics <-
  function(sp_ccc_graph, from_scratch = TRUE) {

    #
    # clean up previous metrics
    #
    if (from_scratch) {
      sp_ccc_graph %<>%
        # clean up nodes
        tidygraph::activate("nodes") %>%
        dplyr::select(-dplyr::starts_with("graph"), -dplyr::starts_with("group")) %>%

        # clean up edges
        tidygraph::activate("edges") %>%
        dplyr::select(-dplyr::starts_with("graph"), -dplyr::starts_with("group"))
    }

    #
    # first compute graph metrics for nodes
    #   then, later copy those to edges
    sp_ccc_graph %<>%

      # assign graph metrics to node
      tidygraph::activate("nodes") %>%
      ## overall graph
      tidygraph::mutate(
        graph_n_nodes = tidygraph::graph_order(),
        graph_n_edges = tidygraph::graph_size(),

        graph_component_count = tidygraph::graph_component_count(),
        graph_motif_count = tidygraph::graph_motif_count(),
        graph_diameter = tidygraph::graph_diameter(directed = TRUE),
        graph_un_diameter = tidygraph::graph_diameter(directed = FALSE),

        # maybe not useful for differentiating
        graph_mean_dist = tidygraph::graph_mean_dist(),

        graph_circuit_rank = graph_n_edges - graph_n_nodes + graph_component_count,
        graph_reciprocity = tidygraph::graph_reciprocity()
      ) %>%
      tidygraph::morph(tidygraph::to_undirected) %>%
      tidygraph::mutate(graph_clique_num = tidygraph::graph_clique_num(),
                        graph_clique_count = tidygraph::graph_clique_count()) %>%
      tidygraph::unmorph() %>%

      ## find subgraphs
      tidygraph::mutate(group = tidygraph::group_components()) %>%
      tidygraph::morph(tidygraph::to_components) %>%

      ### for each subgraph
      tidygraph::mutate(
        group_n_nodes = tidygraph::graph_order(),
        group_n_edges = tidygraph::graph_size(),

        group_adhesion = tidygraph::graph_adhesion(),

        group_motif_count = tidygraph::graph_motif_count(),
        group_diameter = tidygraph::graph_diameter(directed = TRUE),
        group_un_diameter = tidygraph::graph_diameter(directed = FALSE),

        # maybe not useful for differentiating
        group_mean_dist = tidygraph::graph_mean_dist(),

        # probably not informative
        group_girth = tidygraph::graph_girth(),

        # group_n_edges - group_n_nodes + group_componet_count (=1)
        group_circuit_rank = group_n_edges - group_n_nodes + 1,
        group_reciprocity = tidygraph::graph_reciprocity()

      ) %>%
      tidygraph::unmorph()

    sp_ccc_graph %>%
      add_spatial_ccc_graph_metrics_to_edges()
  }


#' Copy graph metrics from nodes to edges
#'
#' Internal function
#' @inheritParams add_spatial_ccc_graph_metrics
#'
#' @noRd
add_spatial_ccc_graph_metrics_to_edges <-
  function(sp_ccc_graph, from_scratch = TRUE) {
    sp_ccc_graph_nodes_df <-
      sp_ccc_graph %>%
      tidygraph::activate("nodes") %>%
      dplyr::as_tibble()

    if (from_scratch) {
      sp_ccc_graph %<>%
        tidygraph::activate("edges") %>%
        # clean up previous metrics
        dplyr::select(-dplyr::starts_with("graph"),-dplyr::starts_with("group"))
    }

    sp_ccc_graph %>%
      dplyr::left_join(
        sp_ccc_graph_nodes_df %>%
          dplyr::select(name, dplyr::starts_with("graph_"), group, dplyr::starts_with("group_")),
        by = c("src" = "name")
      )
  }


#' Extract spatial CCC graph metrics
#'
#' @param ccc_graph spatial CCC graph, in `tidygraph` format,
#'                  an element in the output of [to_spatial_ccc_graph_list()]
#'                  or [compute_spatial_ccc_graph_list()]
#' @inheritParams summarize_ccc_graph_metrics
#'
#' @returns one row data.frame with graph metrics
#'
#' @noRd
extract_ccc_graph_metrics <- function(ccc_graph,
                                      level = c("graph", "group")) {

  if (length(level) > 1) {
    level <- level[1]
  }

  level_option <- c("graph", "group")

  if (length(level) != 1 ||
      is.na(pmatch(level, level_option))) {
    cli::cli_abort(message = "{.var level} should be a string, either 'graph' or 'group'")
  }

  level <- level_option[pmatch(level, level_option)]

  if (level == "graph") {
    ccc_graph %>%
      tidygraph::activate("nodes") %>%
      tibble::as_tibble() %>%
      dplyr::select(dplyr::starts_with("graph_")) %>%
      dplyr::slice_head(n = 1)
  } else {
    ccc_graph %>%
      tidygraph::activate("nodes") %>%
      tibble::as_tibble() %>%
      dplyr::select(group, dplyr::starts_with("group_")) %>%
      dplyr::distinct()
  }
}


#' Summarize CCC table by (cell, LR)
#'
#' internal function
#'
#' @param ccc_tbl CCC graph table.  see [compute_spatial_ccc_tbl()].
#'
#' @export
summarise_ccc_by_cell <-
  function(ccc_tbl) {
    src_summary <-
      ccc_tbl %>%
      dplyr::group_by(src) %>%
      dplyr::summarise(
        n.src = dplyr::n(),  # to avoid confusion between tidygraph::n() and dplyr::n()
        LRscore.sum.src = sum(LRscore),
        # WLRscore.sum.src = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = src) %>%
      dplyr::mutate(node_is_src = TRUE)

    dst_summary <-
      ccc_tbl %>%
      dplyr::group_by(dst) %>%
      dplyr::summarise(
        n.dst = dplyr::n(),  # to avoid confusion between tidygraph::n() and dplyr::n()
        LRscore.sum.dst = sum(LRscore),
        # WLRscore.sum.dst = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = dst) %>%
      dplyr::mutate(node_is_dst = TRUE)

    cell_summary <-
      dplyr::full_join(src_summary,
                       dst_summary,
                       by = c("cell_id"))

    # replace all NAs with 0
    cell_summary[is.na(cell_summary)] <- 0

    cell_summary %>%
      dplyr::relocate(node_is_src, node_is_dst, .after = dplyr::last_col()) %>%
      dplyr::mutate(n.inflow = n.dst - n.src,
                    # WLRscore.sum.inflow = WLRscore.sum.dst - WLRscore.sum.src,
                    LRscore.sum.inflow = LRscore.sum.dst - LRscore.sum.src)
  }





#' Convert CCC table to barebone spatial CCC graph
#'
#' @param ccc_tbl an output of [compute_spatial_ccc()]
#'
#' @returns CCC graph (barebone)
#'
#' @noRd
to_barebone_spatial_ccc_graph <-
  function(ccc_tbl) {
    ccc_tbl %>%
      dplyr::mutate(from = src,
                    to = dst) %>%
      tidygraph::as_tbl_graph(directed = TRUE)
  }



