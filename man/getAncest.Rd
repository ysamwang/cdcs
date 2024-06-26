% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAncest.R
\name{getAncest}
\alias{getAncest}
\title{Confidence sets for Causal Discovery}
\usage{
getAncest(tab)
}
\arguments{
\item{tab}{a matrix with p columns (corresponding to variables) and each row is an ordering}
}
\value{
a p x p matrix where A[i,j] indicates the proportion of orderings where j precedes i
}
\description{
Take a set of orderings and see which ancestral relations hold
}
