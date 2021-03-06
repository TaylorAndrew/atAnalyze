% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multAovInteraction.R
\name{multAovInteraction}
\alias{multAovInteraction}
\title{multAovInteraction}
\usage{
multAovInteraction(data, grp_var1, grp_var2, varlist, round = 2, pround = 3,
  includeN = F)
}
\arguments{
\item{data}{Input dataset}

\item{grp_var1}{Variable name of variable 1}

\item{grp_var2}{Variable name of variable 2}

\item{varlist}{Vector of variable names to be used as the DV(s)}

\item{round}{Number of decimals for descriptive statistics}

\item{pround}{Number of decimals for p-values}

\item{includeN}{Logical, whether to include the N for descriptives}
}
\value{
data.frame containing the output summary
}
\description{
multAovInteraction provides descriptive summary of multiple depedent variables, stratified by two grouping variables, along with p-values from aov() models testing the main effects and interaction between the two grouping variables.
}
\examples{
#df<-data.frame(grp1=sample(c(1,2),100,replace=TRUE),
#                grp2=sample(c(1,2,3),100,replace=TRUE),
#                var1=rnorm(100,32,12),
#                var2=runif(100),
#                var3=rnorm(100))
#multAovInteraction(df, grp_var1 = "grp1", grp_var1 = "grp2", varlist = c("var2","var3","var1"))
}

