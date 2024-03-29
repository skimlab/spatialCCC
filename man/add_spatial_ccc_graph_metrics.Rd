% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_ccc.R
\name{add_spatial_ccc_graph_metrics}
\alias{add_spatial_ccc_graph_metrics}
\title{Add various graph metrics to spatial CCC graph}
\usage{
add_spatial_ccc_graph_metrics(sp_ccc_graph, from_scratch = TRUE)
}
\arguments{
\item{sp_ccc_graph}{an output of \code{\link[=to_spatial_ccc_graph]{to_spatial_ccc_graph()}}}

\item{from_scratch}{if TRUE, the existing metrics are wiped clean.}
}
\description{
See below and \link[tidygraph:graph_measures]{tidygraph::graph_measures} for how these metrics are computed.
}
