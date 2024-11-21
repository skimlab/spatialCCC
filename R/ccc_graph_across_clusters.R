#' Cross-reference ccc_graph across clusters
#'
#' @description
#' Cross-reference *ccc_graph* across clusters defined for each spot and
#'   test the frequency of inter-cluster communications and
#'   calculate the likelihood of the frequency by chance.
#'
#' @details
#'   ccc_graph is summarized in *cog_tbl* and
#'   spatial configuration is stored in *sp_dist* which
#'   also includes (cluster_src, cluster_dst).
#'
#' **Step 1** Create base distribution by cross-referencing
#'   *sp_dist* against cluster configuration
#'
#' **Step 2** Map ccc_graph over clusters and take the average across
#'   LRs to create "null" distribution of inter/intra-cluster
#'   interaction frequency.
#'
#' **Step 3** Calculate the likelihood of observed inter/intra-cluster
#'   interaction frequency for each LR by chance using the "null"
#'   distribution estimated in **Step 2**.
#'
#' @param cog_tbl table (data.frame) converted from ccc_graph_list.
#'                Node and edge statistics are calculated for each
#'                (LR, cluster_src, cluster_dst) tuple.
#' @param sp_dist an output of calc_spot_dist() amended by cluster IDs
#'                (cluster_src, cluster_dst) should be used for cluster IDs.
#'                The easiest way to make sure this is satisfied is to add
#'                "cluster" column in SpatialExperiment object with cluster/group
#'                information.  This will be automatically used in
#'                calc_spot_dist() function when it creates spot_dist data frame.
#'
#' @returns data frame with estimated p-values
#'
#' @export
#'
map_ccc_graph_across_clusters <- function(cog_tbl, sp_dist) {
  sp_dist_base_by_clust <-
    sp_dist %>%
    # filter(cluster_src != cluster_dst) %>%
    group_by(cluster_src, cluster_dst) %>%
    count_interactions_across_clusters_ccc_graph_tbl() %>%
    select(-n_clust, -n_edges_between, -n_edges_within) %>%
    mutate(ccc_type = ifelse(cluster_src != cluster_dst,
                             "between", "within"))

  # summary of graph for # of edges and nodes
  #   by LR, cluster_src, cluster_dst
  z <-
    cog_tbl %>%
    group_by(LR, cluster_src, cluster_dst) %>%
    count_interactions_across_clusters_ccc_graph_tbl()

  # summary of graph by LR
  z1 <-
    cog_tbl %>%
    group_by(LR) %>%
    count_interactions_across_clusters_ccc_graph_tbl()

  # averaging # of edges and nodes over LRs
  #
  # outcome is average number of nodes and edges for
  #   (cluster_src, cluster_dst) pair
  z_avg <-
    z %>%
    group_by(cluster_src, cluster_dst) %>%
    average_ccc_graph_tbl_summary()

  # this will work as null distribution
  #   to compute p-value for (LR, cluster_src, cluster_dst) tuple
  #
  # outcome is also by (cluster_src, cluster_dst) pair
  z_base <-
    sp_dist_base_by_clust %>%
    left_join(
      z_avg,
      by = c("cluster_src", "cluster_dst"),
      suffix = c(".base", ".avg")
    )

  # combining parameters for actual obs (z) and
  #   null distribution (z_base)
  #   to compute p-value for (LR, cluster_src, cluster_dst) tuple
  z %>%
    left_join(z_base,
              by = c("cluster_src", "cluster_dst")) %>%
    left_join(z1,
              by = "LR",
              suffix = c("", ".LR")) %>%
    # now actually computing p-values
    calc_ccc_graph_summary_pval
}


#' Summarize p-value across clusters
#'
#' Summarize p-value (LR, cluster_src, cluster_dst)
#'   across (cluster_src, cluster_dst).
#'
#' @param LR_pval_across data.frame, an output of map_ccc_graph_across_clusters.
#'    the columns include cluster_src, cluster_dst, n_nodes, p_val_node,
#'    n_edges, p_val_edge
#'
#' @export
summarize_by_LR_pval_across_clusters <- function(LR_pval_across) {
  # compute the average p_val_nodes by
  #   taking weighted sum of p_val_nodes/p_val_edges with n_nodes/n_edges as weight
  #          divided by the sum of n_nodes/n_edges
  #          over (cluster_src, cluster_dst) pairs
  #
  # **within** means nodes are within the same cluster
  # **between** means nodes are between two different clusters
  #
  LR_pval_within <-
    LR_pval_across %>%
    filter(cluster_src == cluster_dst) %>%
    group_by(LR) %>%
    summarise(
      # is taking weighted average a reasonable method to aggregate p-values?
      p_val_node_mean = sum(n_nodes * p_val_node) / sum(n_nodes),  # weighted average
      p_val_edge_mean = sum(n_edges * p_val_edge) / sum(n_edges)   # weighted average
    ) %>%
    arrange(p_val_edge_mean)

  LR_pval_between <-
    LR_pval_across %>%
    filter(cluster_src != cluster_dst) %>%
    group_by(LR) %>%
    summarise(
      p_val_node_mean = sum(n_nodes * p_val_node) / sum(n_nodes),
      p_val_edge_mean = sum(n_edges * p_val_edge) / sum(n_edges)
    ) %>%
    arrange(p_val_edge_mean)

  # combine cellualar interactions both within and between clusters
  LR_pval <-
    full_join(
      LR_pval_between,
      LR_pval_within,
      by = "LR",
      suffix = c("_between", "_within")
    )

  # combine p-values from both "within" and "between"
  LR_pval[["p_val_node"]] <-
    LR_pval %>%
    select(starts_with("p_val_node")) %>%
    apply(MARGIN = 1,
          function(x) {
            if (any(is.na(x))) {
              NA
            } else {
              r <- metap::sumlog(x)
              r$p
            }
          })

  # combine p-values from both "within" and "between"
  LR_pval[["p_val_edge"]] <-
    LR_pval %>%
    select(starts_with("p_val_edge")) %>%
    apply(MARGIN = 1,
          function(x) {
            if (any(is.na(x))) {
              NA
            } else {
              r <- metap::sumlog(x)
              r$p
            }
          })

  LR_pval %>%
    # adding nodes/edges statistics for LR
    left_join(LR_pval_across %>%
                select(LR, ends_with(".LR")) %>%
                distinct(),
              by = "LR") %>%
    relocate(LR, p_val_node, p_val_edge)
}

#' Calculate p-value for LR in ccc_graph
#'
#' @rdname map_ccc_graph_across_clusters
#'
#' @details
#' *map_ccc_graph_across_clusters_shortcut()* maps and summarize
#' the counts by LR only, not by (LR, cluster_src, cluster_dst) like
#' map_ccc_graph_across_clusters().
#'
#' @export
map_ccc_graph_across_clusters_shortcut <- function(cog_tbl, sp_dist) {
  sp_dist_base <-
    sp_dist %>%
    count_interactions_across_clusters_ccc_graph_tbl() %>%
    select(-n_clust, -n_edges_between, -n_edges_within)

  # summary by LR
  z1 <-
    cog_tbl %>%
    group_by(LR) %>%
    count_interactions_across_clusters_ccc_graph_tbl()

  # averaging # of edges and nodes over LRs
  #
  # outcome is average number of nodes and edges for
  #   the whole ccc_graph
  z1_avg <-
    z1 %>%
    average_ccc_graph_tbl_summary()

  z1_merged <-
    z1 %>%
    bind_cols(sp_dist_base %>%
                rename_with(~ paste0(.x, ".base"))) %>%
    bind_cols(z1_avg %>%
                rename_with(~ paste0(.x, ".avg")))

  # now actually computing p-values
  z1_merged_pval <-
    calc_ccc_graph_summary_pval(z1_merged)

  z1_merged_pval
}



# Don't need this any more
#
# merge_ccc_graph_pval <- function(cog_pval_LR, cog_pval) {
#   left_join(
#     cog_pval_LR,
#     cog_pval %>%
#       select(LR, starts_with("adj_p")),
#     by = "LR"
#   )
# }





#' Count interactions across clusters
#'
#' @inheritParams map_ccc_graph_across_clusters
#' @param include_edges If TRUE (default), include edge counts.
#'
#' @returns data frame with the counts of nodes and edges interacting
#'   across clusters
count_interactions_across_clusters_ccc_graph_tbl <-
  function(cog_tbl, include_edges = TRUE) {
    res <-
      cog_tbl %>%
      summarise(
        n_clust = length(unique(c(
          cluster_src, cluster_dst
        ))),

        n_nodes = length(unique(c(src, dst))),
        n_nodes_src = length(unique(src)),
        n_nodes_dst = length(unique(dst)),

        # summarizing edge stats overall/between/within clusters
        n_edges = n(),
        n_edges_between = sum(cluster_src != cluster_dst),
        n_edges_within = sum(cluster_src == cluster_dst),

        .groups = "drop"
      )

    if (!include_edges) {
      res <-
        res %>% select(-n_edges,-n_edges_between,-n_edges_within)
    }

    res
  }

#'
#' @rdname count_interactions_across_clusters_ccc_graph_tbl
#'
summarize_ccc_graph_tbl <- function(cog_tbl, include_edges = TRUE) {
  count_interactions_across_clusters_ccc_graph_tbl(cog_tbl, include_edges)
}

#'
#' @rdname count_interactions_across_clusters_ccc_graph_tbl
#' @param cog spatial CCC graph
#'
count_interactions_across_clusters_ccc_graph <- function(cog, include_edges = TRUE) {
  cog %>%
    tidygraph::activate("edges") %>%
    tibble::as_tibble() %>%
    count_interactions_across_clusters_ccc_graph_tbl(include_edges)
}

#'
#' @rdname count_interactions_across_clusters_ccc_graph_tbl
summarize_ccc_graph <- function(cog, include_edges = TRUE) {
  count_interactions_across_clusters_ccc_graph(
    cog, include_edges
  )
}


# In this, the means are actually being over-estimated
#   because we are not including those w/ no results
#   i.e. n_nodes_* == 0 and/or n_edges == 0
# To be more precise, it should be sum(.) / # of LRs in LRdb
# However, this will produce more "conservative" results
average_ccc_graph_tbl_summary <-
  function(cog_tbl_summary) {
    cog_tbl_summary %>%
      summarize(
        n_nodes = mean(n_nodes),
        n_nodes_src = mean(n_nodes_src),
        n_nodes_dst = mean(n_nodes_dst),
        n_edges = mean(n_edges),
        .groups = "drop"
      )
  }

# now actually computing p-values
calc_ccc_graph_summary_pval <- function(z_m) {
  z_m %>%
    mutate(
      p_val_edge = pbinom(
        q = n_edges,
        size = n_edges.base,
        prob = n_edges.avg / n_edges.base,
        lower.tail = FALSE
      ),
      p_val_node = pbinom(
        q = n_nodes,
        size = n_nodes.base,
        prob = n_nodes.avg / n_nodes.base,
        lower.tail = FALSE
      )
    ) %>%
    mutate(
      adj_p_val_edge = p.adjust(p_val_edge),
      adj_p_val_node = p.adjust(p_val_node)
    ) %>%
    relocate(p_val_edge,
             adj_p_val_edge,
             .before = "n_edges") %>%
    relocate(p_val_node,
             adj_p_val_node,
             .before = "n_nodes")
}
