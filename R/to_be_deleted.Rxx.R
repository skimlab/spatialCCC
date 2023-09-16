#
# TO-BE-DELETED =====
#

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
#' @param workers the number of processes to be used for parallel
#'   processing
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
#' @noRd
compute_spatial_ccc_tbl_list <-
  function(spe,
           assay_name = "logcounts",
           LRdb,
           expression_min_prop = 0.05,
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5,
           workers = 1) {
    spot_dist <- calc_spot_dist(spe)

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

    # make sure to process only LR pairs present in gene expression data
    LRdb %<>%
      dplyr::filter(
        ligand_gene_symbol %in% rowAnnots$symbol,
        receptor_gene_symbol %in% rowAnnots$symbol
      )

    # this seems to be faster, than using apply(.., MARGIN, ...)

    compute_lr_score <- function(k) {
      # ligand and receptor names
      lig <- LRdb$ligand_gene_symbol[k]
      rec <- LRdb$receptor_gene_symbol[k]

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

    multiworkers <- function(Ks,
                             workers) {
      # set up future::plan and its release on exit
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential), add = TRUE)

      furrr::future_map(Ks,
                        compute_lr_score,
                        # making sure "dplyr" package is attached to future parallel processes
                        .options = furrr::furrr_options(seed = TRUE,
                                                        packages = c("dplyr")))
    }


    Ks <- setNames(1:nrow(LRdb), nm = LRdb$LR)

    if (workers > 1) {
      # ccc_list
      ccc_ls <-
        multiworkers(Ks, workers)
    } else {
      ccc_ls <-
        # single worker
        purrr::map(Ks, compute_lr_score)
    }

    dplyr::bind_rows(ccc_ls,
                     .id = "LR") %>%
      dplyr::relocate(LR, ligand, receptor, .after = "norm.d")
  }


#' Convert CCC table to CCC graph
#'
#' @inheritParams to_spatial_ccc_graph_list
#' @param LR_of_interest LR pair name in string.  If NULL (default),
#'   CCC graph will be created using all LR pairs.
#'
#' @returns spatial CCC graph -- a `tidygraph`
#'
#' @inheritSection compute_spatial_ccc_graph_list graph metrics
#'
#' @noRd
to_spatial_ccc_graph_old <-
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


#' Convert CCC table to a list of spatial CCC graphs
#'
#' @inheritParams to_spatial_ccc_graph
#' @param workers the number of processes to be used for parallel
#'   processing
#'
#' @inheritSection compute_spatial_ccc_graph_list graph metrics
#'
#' @returns a list of spatial CCC graph
#'
#'
#' @noRd
to_spatial_ccc_graph_list <-
  function(ccc_tbl,
           spatial_data,
           workers = 1) {
    # staging list of LR pairs to process
    LR_ls <-
      ccc_tbl$LR %>%
      unique()

    names(LR_ls) <- LR_ls

    if (workers > 1) {
      # set up future::plan and its release on exit
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential), add = TRUE)

      ccc_graph_ls <-
        furrr::future_map(LR_ls,
                          function(LR_of_interest) {
                            ccc_tbl %>%
                              dplyr::filter(LR == LR_of_interest) %>%
                              to_spatial_ccc_graph(ccc_tbl,
                                                   spatial_data)
                          },
                          # making sure c("tidygraph", "dplyr") packages are attached to future parallel processes
                          .options = furrr::furrr_options(seed = TRUE,
                                                          packages = c("tidygraph", "dplyr")))
    } else {
      ccc_graph_ls <-
        purrr::map(LR_ls,
                   function(LR_of_interest) {
                     ccc_tbl %>%
                       dplyr::filter(LR == LR_of_interest) %>%
                       to_spatial_ccc_graph(ccc_tbl,
                                            spatial_data)
                   })

    }

    ccc_graph_ls

    # purrr::map(LR_ls,
    #            function(LR_of_interest) {
    #              to_spatial_ccc_graph(ccc_tbl,
    #                                   spatial_data,
    #                                   LR_of_interest)
    #            })
  }


#' Compute cell-cell communication graphs
#'
#' @inheritParams compute_spatial_ccc_graph
#' @param LRdb LRdb in data.frame (see [get_LRdb()])
#' @param workers the number of processes to be used for parallel
#'   processing
#'
#' @returns a list of spatial CCC graphs, each of which is in `tidygraph` format;
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
#'
#' @noRd
compute_spatial_ccc_graph_list_obsolete <-
  function(spe,
           assay_name = "logcounts",
           LRdb,
           expression_min_prop = 0.05,
           spot_dist_cutoff = 1.5,
           LRscore_cutoff = 0.5,
           workers = 1) {

    sp_col_data <-
      get_spatial_data(spe)

    compute_spatial_ccc_tbl_list(
      spe,
      assay_name,
      LRdb,
      expression_min_prop,
      spot_dist_cutoff,
      LRscore_cutoff,
      workers = workers
    ) %>%
      to_spatial_ccc_graph_list(spatial_data = sp_col_data,
                                workers = workers)
  }


#' Summarize CCC table by (cell, LR)
#'
#' internal function
#'
#' @param ccc_tbl CCC graph table.  see [compute_spatial_ccc_tbl()].
#'
#' @noRd
summarise_ccc_by_cell_lr <-
  function(ccc_tbl) {
    src_summary <-
      ccc_tbl %>%
      dplyr::group_by(src, LR) %>%
      dplyr::summarise(
        src.n = dplyr::n(),  # to avoid confusion between tidygraph::n() and dplyr::n()
        src.WLR_total = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = src) %>%
      dplyr::mutate(node_is_src = TRUE)

    dst_summary <-
      ccc_tbl %>%
      dplyr::group_by(dst, LR) %>%
      dplyr::summarise(
        dst.n = dplyr::n(),  # to avoid confusion between tidygraph::n() and dplyr::n()
        dst.WLR_total = sum(WLRscore),
        .groups = "drop"
      ) %>%
      dplyr::rename(cell_id = dst) %>%
      dplyr::mutate(node_is_dst = TRUE)

    cell_lr_summary <-
      dplyr::full_join(src_summary,
                       dst_summary,
                       by = c("cell_id", "LR"))

    # replace all NAs with 0
    cell_lr_summary[is.na(cell_lr_summary)] <- 0

    cell_lr_summary %>%
      dplyr::relocate(node_is_src, node_is_dst, .after = dplyr::last_col()) %>%
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
      n = dplyr::n(),
      LR = paste(LR, collapse = ";"),
      src.n = mean(src.n),
      src.WLR_total = mean(src.WLR_total),
      dst.n = mean(dst.n),
      dst.WLR_total = mean(dst.WLR_total),
      node_is_src = any(node_is_src),
      node_is_dst = any(node_is_dst),
      inflow.n = mean(inflow.n),
      inflow.WLR_total = mean(inflow.WLR_total)
    )
}
