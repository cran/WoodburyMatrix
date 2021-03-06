% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahalanobis.R
\name{mahalanobis}
\alias{mahalanobis}
\alias{mahalanobis,ANY,ANY,SWoodburyMatrix-method}
\title{Mahalanobis distance}
\usage{
mahalanobis(x, center, cov, inverted = FALSE, ...)

\S4method{mahalanobis}{ANY,ANY,SWoodburyMatrix}(x, center, cov, inverted = FALSE, ...)
}
\arguments{
\item{x}{Numeric vector or matrix}

\item{center}{Numeric vector representing the mean; if omitted, defaults to
zero mean}

\item{cov}{Covariance matrix}

\item{inverted}{Whether to treat \code{cov} as a precision matrix; must be
\code{FALSE} for \code{SWoodburyMatrix} objects.}

\item{...}{Passed to the \code{\link[Matrix]{Cholesky}} function.}
}
\description{
Generic for computing the squared Mahalanobis distance.
}
\section{Methods (by class)}{
\itemize{
\item \code{x = ANY,center = ANY,cov = SWoodburyMatrix}: Use the Woodbury matrix identity to compute the
squared Mahalanobis distance with the implicit matrix as the covariance.
}}

\seealso{
\link[stats]{mahalanobis}
}
