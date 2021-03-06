% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multGroup.R
\name{multGroup}
\alias{multGroup}
\title{multGroup}
\usage{
multGroup(data, PcontinVars = NULL, PcatVars = NULL, NPcontinVars = NULL,
  NPcatVars = NULL, SortVars = NULL, varLabelTable = NULL,
  grouping = NULL, pdec = 3, dec = 2, mu = 0, ChiProbabilities = NULL,
  labels = NULL, percent = "column", NPdescriptives = NULL,
  NPvariance = "range", padjust = F, provideP = T, include = "none",
  verbose = T)
}
\arguments{
\item{data}{Input data.frame}

\item{PcontinVars}{Vector of continous variables names from data that should be tested parametrically}

\item{PcatVars}{Vector of categorical variables names from data that should be tested parametrically}

\item{NPcontinVars}{Vector of continous variables names from data that should be tested non-parametrically}

\item{NPcatVars}{Vector of categorical variables names from data that should be tested non-parametrically}

\item{SortVars}{Vector of variable names indicating in what order the variables should be presented in the output table}

\item{varLabelTable}{A varLabelTable object. Converts variable names to variable labels in output table.}

\item{grouping}{Grouping variable of any number of levels that will be used to stratify each variable}

\item{pdec}{Number of decimals for p-values}

\item{dec}{Number of decimals for descriptive statistics}

\item{mu}{For a one sample t-test (no grouping), the reference value to compare population to. Default of 0.}

\item{ChiProbabilities}{For a one sample chi-square (no grouping), the underlying probability distribution to compare the population to. If NULL, will compare to an even distribution across all levels.}

\item{labels}{Labels for the grouping variable}

\item{percent}{Indicates if the output table should present 'row', 'column', or 'overall' percents for categorical variables.}

\item{NPdescriptives}{Vector of continuous variable names that should present median and range in the output table.}

\item{padjust}{If TRUE and grouping has 3 or more levels, only pairwise comparisons that are significant at 0.05 after using a multiple comparisons adjustment are presented}

\item{provideP}{If TRUE, in the case of no grouping variable, the p-values for the 1-sample tests will be included in the output, otherwise only the descriptive statistics are returned}

\item{include}{If 'none', no additional output is included. If 'range', range is included. If '95ci' then 95% Confidence intervals are included.}

\item{verbose}{If TRUE, additional information is printed to the console}
}
\value{
A data.frame including all univariate results.
}
\description{
multGroup provides summary statistics for continous and categorical variables either overall or stratified by a grouping variable. It also provides parametric or nonparametric tests for each variable.
}
\examples{

#df <- data.frame(time = runif(100, 1, 100),
#                  event = sample(c(0,1), 100, replace = T, prob = c(.1, .9)),
#                  x1 = rnorm(100),
#                  x2 = sample(c(1,2), 100, replace = T))
#multGroup(data = df,
#          NPcontinVars = "time",
#          grouping = "event",
#          PcontinVars = "x1",
#          NPcatVars = "x2",
#          keepOutput=T,
#          printOutput = F)
}

