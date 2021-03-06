% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multKM.R
\name{multKM}
\alias{multKM}
\title{multKM}
\usage{
multKM(data, time, event, catList)
}
\arguments{
\item{data}{data.frame containing all data}

\item{time}{Variable name contaning time data}

\item{event}{Variable name containing event data (event = 1, censor = 0)}

\item{catList}{Vector of categorical variable name containing grouping variables}
}
\value{
R data.frame containing all km info
}
\description{
multKM provides iterative Kaplan-Meier survival analysis for a vector of categorical variables, including median survival, event counts, and pairwise comparisons.
}
\examples{
#None
}

