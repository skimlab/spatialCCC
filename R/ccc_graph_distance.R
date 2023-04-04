#' Compute cell overlap distance
#'
#' @param ccc_tbl CCC graph in tabular format
#'
#' @return Jaccard distance in matrix format
#'
#' @export
dist_cell_overlap_ccc_tbl <- function(ccc_tbl) {
  ccc_tbl %>%
    prep_cell_overlap_ccc_tbl() %>%
    dist_cell_overlap_edgesets()
}

#' Compute cell overlap distance
#'
#' @param ccc_tbl CCC graph in tabular format
#'
#' @return Jaccard distance in long format
#'
#' @export
lf_cell_overlap_ccc_tbl <- function(ccc_tbl) {
  ccc_tbl %>%
    prep_cell_overlap_ccc_tbl() %>%
    lf_cell_overlap_edgesets()
}

#' Convert distance in long format to matrix
#'
#' It ensures columns and rows have the same features
#'
#' @param lf distance in long format
#'
#' @return distance in matrix format
#'
#' @export
lf_to_dist <- function(lf) {
  col_name <- colnames(lf)[2]
  row_name <- colnames(lf)[1]
  value_name <- colnames(lf)[3]

  wf <-
    lf %>% pivot_wider(
      id_cols = row_name,
      names_from = {{ col_name }},
      values_from = {{ value_name }},
      values_fill = 1
    )

  full_features <- union(lf[[1]], lf[[2]])
  dd <-
    matrix(
      NA,
      nrow = length(full_features),
      ncol = length(full_features),
      dimnames = list(full_features, full_features)
    )

  dd[wf[[1]], colnames(wf)[-1]] <- as.matrix(wf[-1])

  diag(dd) <- 0
  dd[is.na(dd)] <- 0

  dd + t(dd)
}

# ****************************
# ==== Internal functions ====
# ****************************

#' Prepare CCC table for distance calcuation
#'
#' Reorganize CCC table into a list where each element is a CCC graph,
#' an edge list where each edge, (src, dst) pair, is not concatenated
#' into src_dst string for easy enumeration later.
#'
#' @param ccc_tbl CCC table
#'
#' @return a list of CCC graph which is a list of concatenated edges
prep_cell_overlap_ccc_tbl <- function(ccc_tbl) {
  LR_list <- unique(ccc_tbl$LR)
  names(LR_list) <- LR_list

  LR_list %>%
    map(function(lr) {
      ccc_tbl %>%
        dplyr::filter(LR == lr) %>%
        dplyr::mutate(src_dst = paste(src, dst, sep = "_")) %>%
        dplyr::pull(src_dst)
    })
}

#' Compute distance between LRs
#'
#' @description
#' `dist_cell_overlap_edgesets()` returns Jaccard distance matrix.
#' `lf_cell_overlap_edgesets()` returns Jaccard distance in long format.
#'
#' Calculate Jaccard distance between LRs where each LR has
#' a set of concatenated edges
#'
#' @param edgeset_list a set of concatenated edges
#'   (see [prep_cell_overlap_ccc_tbl()])
#'
dist_cell_overlap_edgesets <- function(edgeset_list) {
  edgeset_names <- names(edgeset_list)

  #
  # prepare empty distance matrix
  #   again distance matrix is pre-made for speed up
  cell_overlap_dist <-
    matrix(0,
           nrow = length(edgeset_names),
           ncol = length(edgeset_names))

  # calculate (Jaccard) indices
  #   only for upper triangle
  for (i in 1:(length(edgeset_names) - 1)) {
    for (j in (i + 1):length(edgeset_names)) {
      cell_overlap_dist[i, j] <-
        length(intersect(edgeset_list[[i]], edgeset_list[[j]])) /
        length(union(edgeset_list[[i]], edgeset_list[[j]]))
    }
  }

  cell_overlap_dist <- cell_overlap_dist + t(cell_overlap_dist)

  # set diagonal to 1 so that d(x,x) = 0 (see below)
  diag(cell_overlap_dist) <- 1

  colnames(cell_overlap_dist) <-
    rownames(cell_overlap_dist) <- edgeset_names

  # convert to distance
  1 - cell_overlap_dist
}

#' @rdname dist_cell_overlap_edgesets
lf_cell_overlap_edgesets <- function(edgeset_list) {
  edgeset_names <- names(edgeset_list)

  #
  # prepare distance matrix
  #   again distance matrix is pre-made for speed up
  cell_overlap_lf <-
    combn(edgeset_names, 2) %>%
    t() %>%
    data.frame()

  colnames(cell_overlap_lf) <- c("lr1", "lr2")

  cell_overlap_lf %<>%
    mutate(
      n1 = furrr::future_map_int(.$lr1,
                                 function(lr_1) {
                                   length(edgeset_list[[lr_1]])
                                 }),
      n2 = furrr::future_map_int(.$lr2,
                                 function(lr_2) {
                                   length(edgeset_list[[lr_2]])
                                 }),
      intersect =
        furrr::future_map2_int(.$lr1, .$lr2,
                               function(lr_1, lr_2) {
                                 length(intersect(edgeset_list[[lr_1]], edgeset_list[[lr_2]]))
                               }),
      union =
        furrr::future_map2_int(.$lr1, .$lr2,
                               function(lr_1, lr_2) {
                                 length(union(edgeset_list[[lr_1]], edgeset_list[[lr_2]]))
                               })
    )

  cell_overlap_lf %>%
    mutate(d = 1 - intersect / union)
}

