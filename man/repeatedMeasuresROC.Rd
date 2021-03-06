% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/repeatedMeasuresROC.R
\name{repeatedMeasuresROC}
\alias{repeatedMeasuresROC}
\title{repeatedMeasuresROC}
\usage{
repeatedMeasuresROC(glmerModel)
}
\arguments{
\item{glmerModel}{an lme4::glmer() model}
}
\value{
a list containing: ROC_table and AUC
}
\description{
repeatedMeasuresROC provides a table of sensitivities and specificities, as well as the AUC using Wilcoxon non-parametric approach, both in line with Liu and Wu (2003)
}
\examples{
hdp <- read.csv("http://www.ats.ucla.edu/stat/data/hdp.csv")
hdp <- within(hdp, {
  Married <- factor(Married, levels = 0:1, labels = c("no", "yes"))
  DID <- factor(DID)
  HID <- factor(HID)
})
library(lme4)
m <- glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience + (1 | DID),
           data = hdp,
           family = binomial,
           control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 10)
out <- repeatedMeasuresROC(m)
plot(x = 1-out$ROC_table$spec, y = out$ROC_table$sens)
abline(0,1)
(AUC = out$AUC)
}

