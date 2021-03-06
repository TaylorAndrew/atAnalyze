% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/descriptiveConfInts.R
\name{ciProp}
\alias{ciProp}
\title{ciProp}
\usage{
ciProp(x, conf.int = 0.95)
}
\arguments{
\item{x}{Vector of categorical data}

\item{conf.int}{confidence interval, between .01 and .99.}
}
\value{
Vector containing the lower and upper confidence interval for each level of x
}
\description{
ciProp provides confidence intervals for proportions.
}
\examples{
#x <- sample(c(1,2,3), 100, replace = T)
#ciProp(x)
}

