#' corrStars
#'
#' corrStars computes a correlation matrix, along with asterisks to indicate significance.
#'
#' @param x Numeric matrix/data.frame
#' @param method correlation method: spearman, pearson, kendall
#' @param dec number of decimal places for correlation coefficients
#'
#' @return data.frame containing correlation table
#' @export
#'
#' @examples
#' #NULL
corrStars <- function(x, method = "spearman", dec = 2) {
  x <- as.matrix(x)
  R <- rcorr(x,type = method)$r
  p <- rcorr(x,type = method)$P
  ## define notions for significance levels; spacing is important.
  mystars <-
    ifelse(p <= .001, "***", ifelse(p <= .01, "** ", ifelse(p <= .05, "*", " ")))
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(
    x
  )), R), dec))[,-1]
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep = "")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep = "")
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew) - 1])
  print("Note: ***: p \u2264 0.001; ** : p \u2264 0.01; *: p \u2264 0.05, no asterisk: p > 0.05")
  return(Rnew)
}
