% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_ccc.R
\name{compute_spatial_ccc_graph}
\alias{compute_spatial_ccc_graph}
\title{Compute a spatial CCC graph}
\usage{
compute_spatial_ccc_graph(
  spe,
  assay_name = "logcounts",
  ligand,
  receptor,
  expression_min_prop = 0.05,
  spot_dist = calc_spot_dist(spe),
  spot_dist_cutoff = 1.5,
  LRscore_type = "sqrt.additive",
  LRpvalue_cutoff = 0.05,
  LRscore_cutoff = 0.5
)
}
\arguments{
\item{spe}{SpatialExperiment object}

\item{assay_name}{assay name in string}

\item{ligand}{ligand}

\item{receptor}{receptor}

\item{expression_min_prop}{minimum proportion of samples
with non-zero expression value (default: 0.05)}

\item{spot_dist}{distances between spots, an output of \code{\link[=calc_spot_dist]{calc_spot_dist()}}}

\item{spot_dist_cutoff}{cutoff value for norm.d in spot
distances. Default cutoff is 1.5 and
it comes from sqrt(3) ~ 1.73, based on how Visium spot-array is
arranged, so by default, only computing nearest neighbors.}

\item{LRscore_type}{LRscore type: "sqrt.prod" or "sqrt.additive"}

\item{LRpvalue_cutoff}{LRpvalue cutoff}

\item{LRscore_cutoff}{minimum LRscore to keep}
}
\value{
a spatial CCC graphs in \code{tidygraph} format;
spatial CCC graph also include the following graph metrics
}
\description{
Compute a spatial CCC graph
}
\section{graph metrics}{

\enumerate{
\item For overall graph,
\itemize{
\item graph_n_nodes,
\item graph_n_edges
\item graph_component_count
\item graph_motif_count
\item graph_diameter
\item graph_un_diameter
\item graph_mean_dist
\item graph_circuit_rank = graph_n_edges - graph_n_nodes + graph_component_count
\item graph_reciprocity
\item graph_clique_num (sp_ccc_graph assumed as undirected)
\item graph_clique_count (sp_ccc_graph assumed as undirected)
}
\item For each sub-graph (after \code{\link[=group_components]{group_components()}})
\itemize{
\item group_n_nodes
\item group_n_edges
\item group_adhesion
\item group_motif_count
\item group_diameter
\item group_un_diameter
\item group_mean_dist
\item group_girth
\item group_circuit_rank
\item group_reciprocity
}
}
}

