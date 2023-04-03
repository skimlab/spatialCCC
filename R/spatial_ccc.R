#' Calculate distance between spots
#'
#' @param coords spatial coordinates in data.frame
#'   with `cell_id`, `x`, and `y` in that order.  The column headings
#'   do not matter, just the order.
#' @param long_format if TRUE (default), return long format data.frame.
#'   Otherwise, matrix format.
#'
#' @return computed distance including both `d` and `norm.d`.
#'   `norm.d` is a multiple of the shortest distance between spots.
#'
#' @export
calc_spot_dist <- function(coords, long_format = TRUE) {
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
      dplyr::filter(d > 0) %>%
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
#' @return spatial coordinates in data.frame with
#'   `cell_id`, ... (from `spatialCoords(spe)`),
#'   `spot_x`, `spot_y`.
#'
#' @export
get_spatial_coords <- function(spe) {
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
#'    distances ([see calc_spot_dist()])
#' @param LRscore_cutoff minimum LRscore to keep
#'
#' @return LRscore table in data.frame
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
#' ccc_tbl <-
#'   compute_spatial_ccc(
#'     spe = spe,
#'     assay_name = "logcounts",
#'     LRdb = LRdb_mouse)#'
#' }
#'
compute_spatial_ccc <-
  function(spe,
           assay_name = "logcounts",
           LRdb,

           # expression filtering
           #   default is 5%
           expression_min_prop = 0.05,

           #
           # spot distance cutoff
           #   the cutoff is in terms of norm.d (normalized distance)
           #   norm.d is multiple of
           #   minimum of spot-to-spot distances
           #
           # cutoff 1.5 comes from sqrt(3) ~ 1.73
           # so by default, only computing nearest neighbors
           spot_dist_cutoff = 1.5,

           #
           # source:
           #   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7261168/
           #
           LRscore_cutoff = 0.5) {
    # calculate relative distance between cells/spots
    # spe_cd <- spatialCoords(spe) %>%
    #   as_tibble(rownames = "cell_id")
    spe_cd <- get_spatial_coords(spe)

    spot_dist <-
      calc_spot_dist(spe_cd[c("cell_id", "pxl_col_in_fullres", "pxl_row_in_fullres")])

    spot_dist %<>%
      dplyr::filter(norm.d < spot_dist_cutoff)

    gexp <- as.matrix(assays(spe)[[assay_name]])
    rowAnnots <- rowData(spe)

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
      filter(
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
                            #   traveling distance
                            mutate(weight = 1 / norm.d ^ 2) %>%
                            mutate(WLRscore = LRscore * weight) %>%

                            filter(LRscore > LRscore_cutoff) %>%
                            mutate(ligand = lig, receptor = rec)
                        })

    bind_rows(ccc_list,
              .id = "LR") %>%
      relocate(LR, ligand, receptor, .after = "norm.d")
  }
