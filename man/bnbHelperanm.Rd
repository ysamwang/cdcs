% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bnbHelperanm}
\alias{bnbHelperanm}
\title{bnbHelperanm}
\usage{
bnbHelperanm(ancest, children, G, withinAgg, aggType, bs, intercept)
}
\arguments{
\item{ancest}{matrix of ancestor basis functions}

\item{children}{matrix of observed children}

\item{G}{an n x k x |ancest| array of where each variable (corresponding to slice) has k test functions}

\item{withinAgg}{indicates which norm to take when combining test statistics within
a variable but across test functions}

\item{aggType}{indicates which norm to take when combining test statistics across variables}

\item{bs}{number of bootstrap draws}

\item{intercept}{should an additional intercept be included?}
}
\value{
pval: p-value for each child in children
}
\description{
This function is a helper for the branchAndBound_anm function
}
