% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_ccc.R
\name{add_metadata}
\alias{add_metadata}
\title{Add metadata to an object}
\usage{
add_metadata(spe, mdata, replace = FALSE)
}
\arguments{
\item{spe}{SpatialExperiment object}

\item{mdata}{additional metadata to be added \link{data.frame}}

\item{replace}{existing columns are replaced if TRUE, and
the function stops if FALSE, when there is a conflict in column headings.}
}
\value{
SpatialExperiment object
}
\description{
Add metadata to an object
}
