#' Combine two CCC graphs
#'
#' @param ccog1,ccog2 spatial CCC graphs
#'
#' @return flattened spatial CCC graph
#' @export
ccc_graph_union <- function(ccog1, ccog2) {
  flatten_ccc_graph_list(list(ccog1, ccog2))
}

#' Create a new spatial CCC graph with common edges
#'
#' @inheritParams ccc_graph_union
#'
#' @return a spatial CCC graph with common edges
#' @export
ccc_graph_intersect <- function(ccog1, ccog2) {
  ccog2_src_dst <-
    ccog2 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(src, dst, sep = "_")) %>%
    tidygraph::pull(src_dst)

  ccog1 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(src, dst, sep = "_")) %>%
    tidygraph::filter(src_dst %in% ccog2_src_dst) %>%
    tidygraph::select(-src_dst) %>%
    tidy_up_ccc_graph()
}

#' Subtract a spatial CCC graph (ccog2) from another (ccog1)
#'
#' @inheritParams ccc_graph_union
#'
#' @return a spatial CCC graph with common edges
#' @export
ccc_graph_diff <- function(ccog1, ccog2) {
  ccog2_src_dst <-
    ccog2 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(src, dst, sep = "_")) %>%
    tidygraph::pull(src_dst)

  ccog1 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(src, dst, sep = "_")) %>%
    tidygraph::filter(!(src_dst %in% ccog2_src_dst)) %>%
    tidygraph::select(-src_dst) %>%
    tidy_up_ccc_graph()
}
