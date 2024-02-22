

addMetadata <- function(spe, mdata, replace = FALSE) {
  cData <- SummarizedExperiment::colData(spe)
  colnames_cData <- colnames(cData)
  colnames_mdata <- colnames(mdata)

  if (all(length(intersect(colnames_cData, colnames_mdata)) > 0, !replace)) {
    cli_abort(paste("column names conflicts with the existing column names:",
                    paste(intersect(colnames_cData, colnames_mdata), collapse = ", ")))
  }

  mdata <- mdata[rownames(cData), ]

  cData[, colnames_mdata] <- mdata

  SummarizedExperiment::colData(spe) <- cData

  spe
}


summarize_ccc_graph <- function(cog, include_edges = TRUE) {
  cog %>%
    tidygraph::activate("edges") %>%
    tibble::as_tibble() %>%
    summarize_ccc_graph_tbl(include_edges)
}

summarize_ccc_graph_tbl <-
  function(cog_tbl, include_edges = TRUE) {
    res <-
      cog_tbl %>%
      summarise(
        n_clust = length(unique(c(mclust_src, mclust_dst))),

        n_nodes = length(unique(c(src, dst))),
        n_nodes_src = length(unique(src)),
        n_nodes_dst = length(unique(dst)),

        # summarizing edge stats overall/between/within clusters
        n_edges = n(),
        n_edges_between = sum(mclust_src != mclust_dst),
        n_edges_within = sum(mclust_src == mclust_dst),

        .groups = "drop"
      )

    if (!include_edges) {
      res <-
        res %>% select(-n_edges, -n_edges_between, -n_edges_within)
    }

    res
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

# not used in any case
summarize_ccc_graph_tbl_by_LR <- function(cog_tbl) {
  x1 <-
    cog_tbl %>%
    group_by(LR) %>%
    summarize_ccc_graph_tbl() %>%
    mutate(r = n_edges_between / (n_edges_within + n_edges_between))

  x2 <-
    ccc_graph_tbl %>%
    filter(mclust_src != mclust_dst) %>%
    group_by(LR) %>%
    summarize_ccc_graph_tbl(include_edges = FALSE)

  left_join(x1, x2, by = "LR", suffix = c("", "_between")) %>%
    relocate("LR",
             starts_with("n_nodes"),
             starts_with("n_edges"),
             starts_with("n_clust"))
}


# now actually computing p-values
calc_z_merged_pval <- function(z_m) {
  z_m %>%
    mutate(
      p_val_edges = pbinom(
        q = n_edges,
        size = n_edges.base,
        prob = n_edges.avg / n_edges.base,
        lower.tail = FALSE
      ),
      p_val_nodes = pbinom(
        q = n_nodes,
        size = n_nodes.base,
        prob = n_nodes.avg / n_nodes.base,
        lower.tail = FALSE
      )
    ) %>%
    mutate(
      adj_p_val_edges = p.adjust(p_val_edges),
      adj_p_val_nodes = p.adjust(p_val_nodes)
    ) %>%
    relocate(p_val_edges,
             adj_p_val_edges,
             .before = "n_edges") %>%
    relocate(p_val_nodes,
             adj_p_val_nodes,
             .before = "n_nodes")
}

calc_ccc_graph_pval_LR <- function(cog_tbl, sp_dist) {
  sp_dist_base_by_clust <-
    sp_dist %>%
    # filter(mclust_src != mclust_dst) %>%
    group_by(mclust_src, mclust_dst) %>%
    summarize_ccc_graph_tbl() %>%
    select(-n_clust,-n_edges_between,-n_edges_within)

  # summary of graph for # of edges and nodes
  #   by LR, mclust_src, mclust_dst
  z <-
    cog_tbl %>%
    group_by(LR, mclust_src, mclust_dst) %>%
    summarize_ccc_graph_tbl()

  # averaging # of edges and nodes over LRs
  #
  # outcome is average number of nodes and edges for
  #   (mclust_src, mclust_dst) pair
  z_avg <-
    z %>%
    group_by(mclust_src, mclust_dst) %>%
    average_ccc_graph_tbl_summary()

  # this will work as null distribution
  #   to compute p-value for (LR, mclust_src, mclust_dst) tuple
  #
  # outcome is also by (mclust_src, mclust_dst) pair
  z_base <-
    sp_dist_base_by_clust %>%
    left_join(z_avg,
              by = c("mclust_src", "mclust_dst"),
              suffix = c(".base", ".avg"))

  # combining parameters for actual obs (z) and
  #   null distribution (z_base)
  #   to compute p-value for (LR, mclust_src, mclust_dst) tuple
  z_merged <-
    z %>%
    left_join(z_base,
              by = c("mclust_src", "mclust_dst"))

  # now actually computing p-values
  z_merged_pval <-
    calc_z_merged_pval(z_merged)

  # compute the average p_val_nodes by
  #   taking weighted sum of p_val_nodes/p_val_edges with n_nodes/n_edges as weight
  #          divided by the sum of n_nodes/n_edges
  #          over (mclust_src, mclust_dst) pairs
  #
  # **within** means nodes are within the same cluster
  # **between** means nodes are between two different clusters
  z_merged_pval_LR_within <-
    z_merged_pval %>%
    filter(mclust_src == mclust_dst) %>%
    group_by(LR) %>%
    summarise(
      # is taking weighted average a reasonable method to aggregate p-values?
      p_val_nodes_mean = sum(n_nodes * p_val_nodes) / sum(n_nodes),  # weighted average
      p_val_edges_mean = sum(n_edges * p_val_edges) / sum(n_edges)   # weighted average
    ) %>%
    arrange(p_val_edges_mean)

  z_merged_pval_LR_between <-
    z_merged_pval %>%
    filter(mclust_src != mclust_dst) %>%
    group_by(LR) %>%
    summarise(
      p_val_nodes_mean = sum(n_nodes * p_val_nodes) / sum(n_nodes),
      p_val_edges_mean = sum(n_edges * p_val_edges) / sum(n_edges)
    ) %>%
    arrange(p_val_edges_mean)

  z_merged_pval_LR <-
    left_join(
      z_merged_pval_LR_between,
      z_merged_pval_LR_within,
      by = "LR",
      suffix = c("_between", "_within")
    )

  z_merged_pval_LR$p_value_nodes <-
    z_merged_pval_LR %>%
    select(starts_with("p_val_nodes")) %>%
    apply(MARGIN = 1,
          function(x) {
            if (any(is.na(x))) {
              NA
            } else {
              r <- metap::sumlog(x)
              r$p
            }
          })

  z_merged_pval_LR$p_value_edges <-
    z_merged_pval_LR %>%
    select(starts_with("p_val_edges")) %>%
    apply(MARGIN = 1,
          function(x) {
            if (any(is.na(x))) {
              NA
            } else {
              r <- metap::sumlog(x)
              r$p
            }
          })

  z_merged_pval_LR %>%
    relocate(LR, p_value_nodes, p_value_edges)
}


calc_ccc_graph_pval <- function(cog_tbl, sp_dist) {
  sp_dist_base <-
    sp_dist %>%
    summarize_ccc_graph_tbl() %>%
    select(-n_clust, -n_edges_between, -n_edges_within)

  # summary by LR
  z1 <-
    ccc_graph_tbl %>%
    group_by(LR) %>%
    summarize_ccc_graph_tbl()

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
    calc_z_merged_pval(z1_merged)

  z1_merged_pval
}

merge_ccc_graph_pval <- function(z_m_pv_LR, z1_m_pv) {
  pv_m <-
    left_join(
      z_m_pv_LR,
      z1_m_pv %>%
        select(LR, starts_with("adj_p")),
      by = "LR"
    )

  pv_m
}

