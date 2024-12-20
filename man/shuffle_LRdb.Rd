% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_LRdb.R
\name{shuffle_LRdb}
\alias{shuffle_LRdb}
\title{Shuffle LRdb data frame}
\usage{
shuffle_LRdb(lrdb, reverse = TRUE)
}
\arguments{
\item{lrdb}{LRdb table in data frame}

\item{reverse}{LR pairing is reversed if TRUE.
If FALSE, LR pairing is shuffled.}
}
\value{
LRdb table in data frame
}
\description{
Shuffle ligand and receptor pairing, utilized mainly to generate
NULL distribution of LRscore, which later can be used to compute
p-value.
}
