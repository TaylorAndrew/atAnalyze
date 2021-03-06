% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrStars.R
\name{corrStars}
\alias{corrStars}
\title{corrStars}
\usage{
corrStars(x, method = "spearman", dec = 2, includeN = FALSE)
}
\arguments{
\item{x}{Numeric matrix/data.frame}

\item{method}{correlation method: spearman, pearson, kendall}

\item{dec}{number of decimal places for correlation coefficients}

\item{N}{Logical, whether or not to include the N for each correlation}
}
\value{
data.frame containing correlation table
}
\description{
corrStars computes a correlation matrix, along with asterisks to indicate significance.
}
\examples{
#NULL
}

