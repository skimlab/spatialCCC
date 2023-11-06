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
    tidygraph::mutate(src_dst = paste(.data$src, .data$dst, sep = "_")) %>%
    tidygraph::pull(.data$src_dst)

  ccog1 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(.data$src, .data$dst, sep = "_")) %>%
    tidygraph::filter(.data$src_dst %in% ccog2_src_dst) %>%
    tidygraph::select(-.data$src_dst) %>%
    tidy_up_ccc_graph()
}

#' Subtract a spatial CCC graph (ccog2) from another (ccog1)
#'
#' @inheritParams ccc_graph_union
#'
#' @return a spatial CCC graph with common edges
#' @export
ccc_graph_subtract <- function(ccog1, ccog2) {
  ccog2_src_dst <-
    ccog2 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(.data$src, .data$dst, sep = "_")) %>%
    tidygraph::pull(.data$src_dst)

  ccog1 %>%
    tidygraph::activate("edges") %>%
    tidygraph::mutate(src_dst = paste(.data$src, .data$dst, sep = "_")) %>%
    tidygraph::filter(!(.data$src_dst %in% ccog2_src_dst)) %>%
    tidygraph::select(-.data$src_dst) %>%
    tidy_up_ccc_graph()
}


#' Calculate difference between spatial CCC graphs (ccog1, ccog1)
#'
#' @inheritParams ccc_graph_union
#'
#' @return a spatial CCC graph with common edges
#' @export
ccc_graph_diff <- function(ccog1, ccog2) {
  ccc_graph_union(ccog1, ccog2) %>%
    ccc_graph_subtract(ccc_graph_intersect(ccog1, ccog2)) %>%
    tidy_up_ccc_graph()
}
