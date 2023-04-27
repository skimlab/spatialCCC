#' Compute ligand-receptor interactions
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
#' @param LRdb LRdb in data.frame (see [get_LRdb()])
#' @param expression_min_prop minimum proportion of samples
#'    with non-zero expression value (default: 0.05)
#' @param spot_dist_cutoff cutoff value for norm.d in spot
#'    distances. Default cutoff is 1.5 and
#'    it comes from sqrt(3) ~ 1.73, based on how Visium spot-array is
#'    arranged, so by default, only computing nearest neighbors.
#' @param LRscore_cutoff minimum LRscore to keep
#'
#' @return LRscore table in data.frame with the following columns:
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
#'
#' @examples
#' \dontrun{
#' # read Visium spatial transcriptomic data
#' #   and log-normalize data
#' spe <-
#'   SpatialExperiment::read10xVisium(
#'     samples = data_dir,
#'     type = "HDF5",
#'     data = "filtered") %>%
#'     logNormCounts()
#'
#' # get built-in LR datbase for mouse
#' LRdb_mouse <- get_LRdb("mouse")
#'
#' # calculate LRscores
#' ccc_table <-
#'   compute_spatial_ccc(
#'     spe = spe,
#'     assay_name = "logcounts",
#'     LRdb = LRdb_mouse)#'
#' }
#'
compute_spatial_ccc_tbl <-
  function(spe,
           assay_name = "logcounts",
           LRdb,
           expression_min_prop = 0.05,
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5) {

    spot_dist <- calc_spot_dist(spe)

    gexp <- as.matrix(SummarizedExperiment::assays(spe)[[assay_name]])
    rowAnnots <- SingleCellExperiment::rowData(spe)

    # first remove genes with no expression at all
    #   maybe, we need to remove those with expression in only 5% (cutoff) or less
    n_cutoff <- ncol(gexp) * expression_min_prop
    subset_idx <- rowSums(gexp > 0) > n_cutoff
    gexp <- gexp[subset_idx,]
    rowAnnots <- subset(rowAnnots, subset = subset_idx)

    # background/baseline expression
    c_mean <- mean(gexp)

    # make sure to process only LR pairs present in gene expression data
    LRdb %<>%
      dplyr::filter(
        ligand_gene_symbol %in% rowAnnots$symbol,
        receptor_gene_symbol %in% rowAnnots$symbol
      )

    # this seems to be faster, than using apply(.., MARGIN, ...)
    k <- 1:nrow(LRdb)
    names(k) <- LRdb$LR

    ccc_list <-
      furrr::future_map(k,
                        function(k) {
                          # ligand and receptor names
                          lig <- LRdb$ligand_gene_symbol[k]
                          rec <- LRdb$receptor_gene_symbol[k]

                          # ligand and receptor expression values
                          lig_ids <- which(rowAnnots$symbol == lig)
                          rec_ids <- which(rowAnnots$symbol == rec)

                          if (length(lig_ids) > 1) {
                            lig_exp <- colMeans(gexp[lig_ids,])
                          } else {
                            lig_exp <- gexp[lig_ids,]
                          }

                          if (length(rec_ids) > 1) {
                            rec_exp <- colMeans(gexp[rec_ids,])
                          } else {
                            rec_exp <- gexp[rec_ids,]
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
                            dplyr::mutate(weight = ifelse(norm.d == 0, 1, 1 / norm.d^2)) %>%
                            dplyr::mutate(WLRscore = LRscore * weight) %>%

                            dplyr::filter(LRscore > LRscore_cutoff) %>%
                            dplyr::mutate(ligand = lig, receptor = rec)
                        })

    dplyr::bind_rows(ccc_list,
                     .id = "LR") %>%
      dplyr::relocate(LR, ligand, receptor, .after = "norm.d")
  }



#' Compute cell-cell communication graphs
#'
#' @inheritParams compute_spatial_ccc_tbl
#' @inheritParams to_spatial_ccc_graph_list
#'
#' @return a list of spatial CCC graphs, each of which is in `tidygraph` format;
#'   spatial CCC graph also include the following graph metrics to characterize
#'   each CCC graph.
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
compute_spatial_ccc_graph_list <-
  function(spe,
           assay_name = "logcounts",
           LRdb,
           expression_min_prop = 0.05,
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5,
           workers = 1) {

    sp_col_data <-
      get_spatial_data(spe)

    compute_spatial_ccc_tbl(spe,
                            assay_name,
                            LRdb,
                            expression_min_prop,
                            spot_dist_cutoff,
                            LRscore_cutoff) %>%
      to_spatial_ccc_graph_list(
        spatial_data = sp_col_data,
        workers = workers
      )
  }


#' Convert CCC table to a list of spatial CCC graphs
#'
#' @param ccc_tbl an output of [compute_spatial_ccc_tbl()]
#' @param spatial_data a table with spatial coordinates of
#'   spots/cells in spatial transcriptomic data.  if spatial_data
#'   includes cell_id, spot_x, spot_y, they will be used as is;
#'   Otherwise, use the first three columns as if they are
#'   cell_id, spot_x, spot_y
#' @param workers the number of processes to be used for parallel
#'   processing
#'
#' @inheritSection compute_spatial_ccc_graph_list graph metrics
#'
#' @return a list of spatial CCC graph
#'
#' @export
to_spatial_ccc_graph_list <-
  function(ccc_tbl,
           spatial_data,
           workers = 1) {
    # staging list of LR pairs to process
    LR_list <-
      ccc_tbl$LR %>%
      unique()

    names(LR_list) <- LR_list

    future::plan(future::multisession, workers = workers)

    on.exit(future::plan(future::sequential), add = TRUE)

    LR_list %>%
      furrr::future_map(function(LR_of_interest) {
        to_spatial_ccc_graph(ccc_tbl,
                             spatial_data,
                             LR_of_interest)
      },
      .options = furrr::furrr_options(seed = TRUE,
                                      packages = c("tidygraph", "dplyr")))
  }



#' Convert spatial CCC graphs to spatial CCC table
#'
#' @param ccc_graph_list a list of spatial CCC graphs
#'
#' @return spatial CCC table
#'
#' @export
to_spatial_ccc_tbl <-
  function(ccc_graph_list) {
    lapply(ccc_graph_list,
           function(ccog) {
             ccog %>%
               tidygraph::activate("edges") %>%
               tibble::as_tibble() %>%
               dplyr::select(-c("from", "to"))
           }) %>%
      dplyr::bind_rows(.id = "LR")
  }







#' Summarize spatial CCC graph list to table
#'
#' @param ccc_graph_list list of spatial CCC graph,
#'   each of which is an output of to_spatial_ccc_graph.
#' @param level extract graph metrics from either overall graph ("graph") or
#'   subgraph ("group")
#'
#' @return graph metrics table summarized for each LR pair
#'
#' @export
summarize_ccc_graph_metrics <- function(ccc_graph_list,
                                        level = c("graph", "group")) {
  ccc_graph_list %>%
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
      dplyr::select(-any_of(annot_cols.after))
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
#' @return CCC table with amended columns
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
#' @return CCC graph with all isolated nodes removed.
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
#' @return CCC graph with InCOI and tagged columns (see above)
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
#' @inheritParams to_spatial_ccc_graph_list
#' @param LR_of_interest LR pair name in string.  If NULL (default),
#'   CCC graph will be created using all LR pairs.
#'
#' @return spatial CCC graph -- a `tidygraph`
#'
#' @inheritSection compute_spatial_ccc_graph_list graph metrics
#'
#' @noRd
to_spatial_ccc_graph <-
  function(ccc_tbl,
           spatial_data,
           LR_of_interest = NULL) {
    if (!is.null(LR_of_interest)) {
      ccc_tbl <-
        dplyr::filter(ccc_tbl, LR == LR_of_interest)
    }

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
      summarise_ccc_by_cell_lr(ccc_tbl) %>%
      collapse_to_ccc_by_cell() %>%
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



#' Calculate distance between spots
#'
#' @param coords spatial coordinates in data.frame
#'   with `cell_id`, `x`, and `y` in that order.  The column headings
#'   do not matter, just the order.
#' @param long_format if TRUE (default), return long format data.frame.
#'   Otherwise, matrix format.
#'
#' @return A data fame with computed distance including both `d` and `norm.d`.
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

#' Calculate distance between spots
#'
#' @inheritParams compute_spatial_ccc
#' @inheritParams calc_spot_dist0
#'
#' @return A data fame with computed distance including both `d` and `norm.d`.
#'   `norm.d` is a multiple of the shortest distance between spots.
#'
#' @noRd
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

#' Get spatial coordinates
#'
#' Returns default spatial coordinates and scaled coordinates,
#'   `spot_x` and `spot_y`.
#'
#' @param spe SpatialExperiment object
#'
#' @return spatial coordinates in data.frame with
#'   `cell_id`, ... (from `spatialCoords(spe)`),
#'   `spot_x`, `spot_y`.
#'
#' @export
get_spatial_data <- function(spe) {
  spe_col_data <- SingleCellExperiment::colData(spe)
  spe_spatial_coords <- SpatialExperiment::spatialCoords(spe)

  dplyr::left_join(
    tibble::as_tibble(spe_col_data, rownames = "cell_id"),
    tibble::as_tibble(spe_spatial_coords, rownames = "cell_id"),
    by = "cell_id"
  ) %>%
    dplyr::mutate(
      spot_x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe),
      spot_y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe)
    )
}

#' Create barebone CCC graph
#'
#' @inheritParams compute_spatial_ccc
#'
#' @return spatial CCC graph w/o LRscore
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

        group_circuit_rank = graph_n_edges - graph_n_nodes + graph_component_count,
        group_reciprocity = tidygraph::graph_reciprocity()

      ) %>%
      tidygraph::unmorph()

    sp_ccc_graph %>%
      add_spatial_ccc_graph_metrics_to_edges()
  }


#' Extract spatial CCC graph metrics
#'
#' @param ccc_graph spatial CCC graph, in `tidygraph` format,
#'                  an element in the output of [to_spatial_ccc_graph_list()]
#'                  or [compute_spatial_ccc_graph_list()]
#' @inheritParams summarize_ccc_graph_metrics
#'
#' @return one row data.frame with graph metrics
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
#' @param ccc_tbl CCC graph table.  see [compute_spatial_ccc()].
#'
#' @noRd
summarise_ccc_by_cell_lr <-
  function(ccc_tbl) {
    src_summary <-
      ccc_tbl %>%
      dplyr::group_by(src, LR) %>%
      dplyr::summarise(
        src.n = n(),
        src.WLR_total = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = src) %>%
      dplyr::mutate(src = TRUE)

    dst_summary <-
      ccc_tbl %>%
      dplyr::group_by(dst, LR) %>%
      dplyr::summarise(
        dst.n = n(),
        dst.WLR_total = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = dst) %>%
      dplyr::mutate(dst = TRUE)

    cell_lr_summary <-
      dplyr::full_join(src_summary,
                       dst_summary,
                       by = c("cell_id", "LR"))

    # replace all NAs with 0
    cell_lr_summary[is.na(cell_lr_summary)] <- 0

    cell_lr_summary %>%
      dplyr::relocate(src, dst, .after = dplyr::last_col()) %>%
      dplyr::mutate(inflow.n = dst.n - src.n,
                    inflow.WLR_total = dst.WLR_total - src.WLR_total)
  }

#' Aggregate ccc_by_cell_lr to cells
#'
#' Internal function
#' @param ccc_by_cell_lr see [summarise_ccc_by_cell_lr()]
#'
#' @noRd
collapse_to_ccc_by_cell <- function(ccc_by_cell_lr) {
  ccc_by_cell_lr %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      n = n(),
      LR = paste(LR, collapse = ";"),
      src.n = mean(src.n),
      src.WLR_total = mean(src.WLR_total),
      dst.n = mean(dst.n),
      dst.WLR_total = mean(dst.WLR_total),
      src = any(src),
      dst = any(dst),
      inflow.n = mean(inflow.n),
      inflow.WLR_total = mean(inflow.WLR_total)
    )
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

#' Convert CCC table to barebone spatial CCC graph
#'
#' @param ccc_tbl an output of [compute_spatial_ccc()]
#'
#' @return CCC graph (barebone)
#'
#' @noRd
to_barebone_spatial_ccc_graph <-
  function(ccc_tbl) {
    ccc_tbl %>%
      dplyr::mutate(from = src,
                    to = dst) %>%
      tidygraph::as_tbl_graph(directed = TRUE)
  }
