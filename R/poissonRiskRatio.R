#' poissonRiskRatio
#'
#' poissonRiskRatio take a poisson regression with a log link and will output the risk ratio. This function is to be used when an investigator wants the risk ratio instead of the odds ratio, and a glm model with family = binomial(link = 'log') does not converge.
#'
#' @param model a glm object with family = binomial(link = 'log') set
#' @param dec number of  decimals to report parameter estimates
#' @param pdec number of decimals to report p-value
#'
#' @return data.frame containing risk ratio, 95% confidence intervals, and p-values
#' @export
#'
#' @examples
#' #none
poissonRiskRatio <- function(model, dec = 2, pdec = 3) {
    sum1 <- summary(model)
    robustSE <- sqrt(diag(vcovHC(model)))
    sum1$coefficients[,2] <- robustSE
    RR <- round(exp(sum1$coefficients[, 1]), dec)
    pval <- round(sum1$coefficients[, 4], pdec)
    lb <- round(exp(sum1$coefficients[, 1]) - (1.95*sum1$coefficients[, 2]), dec)
    ub <- round(exp(sum1$coefficients[, 1]) + (1.95*sum1$coefficients[, 2]), dec)
    return(data.frame(RR = RR, lb = lb, ub = ub, pval = pval))
}