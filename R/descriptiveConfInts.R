#' ciMean
#'
#' ciMean provides confidence intervals for a mean.
#'
#' @param x Vector of numeric data
#' @param conf.int confidence interval, between .01 and .99.
#'
#' @return Vector containing the lower and upper confidence interval
#' @export
#'
#' @examples
#' #x <- rnorm(100)
#' #ciMean(x)
ciMean <- function(x, conf.int = .95) {
  mean <- mean(x, na.rm =T)
  sd <- sd(x, na.rm = T)
  n  <- length(x[!is.na(x)])
  SE = sd/sqrt(n)
  E = qt(conf.int +(1-conf.int)/2, df = n -1)*SE
  out <- mean + c(-E, E)
  names(out) <- c("ci.lower", "ci.upper")
  return(out)
}

#' ciProp
#'
#' ciProp provides confidence intervals for proportions.
#'
#' @param x Vector of categorical data
#' @param conf.int confidence interval, between .01 and .99.
#'
#' @return Vector containing the lower and upper confidence interval for each level of x
#' @export
#'
#' @examples
#' #x <- sample(c(1,2,3), 100, replace = T)
#' #ciProp(x)
ciProp <- function(x, conf.int = .95) {
  x <- x[!is.na(x)]
  n <- length(x)
  pbar =  unlist(table(x))/n
  SE = sqrt(pbar*(1-pbar)/n)
  E = qnorm(conf.int +(1-conf.int)/2)*SE
  out = do.call(rbind, lapply(1:length(pbar), function(i) pbar[i] + c(-E[i], E[i])))
  row.names(out) <- names(pbar)
  colnames(out) <- c("ci.lower", "ci.upper")
  return(out)
}