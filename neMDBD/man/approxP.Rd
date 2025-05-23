% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neMDBD.R
\name{approxP}
\alias{approxP}
\title{approxP}
\usage{
approxP(x_new, X, beta, P_unnormalized, perplexity)
}
\arguments{
\item{x_new}{The new x.}

\item{X}{The rest of X matrix.}

\item{beta}{beta for corresponding X.}

\item{P_unnormalized}{The unnormalized P matrix rearranged to put x_new as the first one.}
}
\value{
approximated P matrix when changing x_new.
}
\description{
R function to calculate approximated P matrix of t-SNE
}
\examples{
approxP(x_new, X, beta, P_unnormalized, perplexity)
}
