#' corrStars
#'
#' corrStars computes a correlation matrix, along with asterisks to indicate significance.
#'
#' @param x Numeric matrix/data.frame
#' @param method correlation method: spearman, pearson, kendall
#' @param dec number of decimal places for correlation coefficients
#' @param N Logical, whether or not to include the N for each correlation
#'
#' @return data.frame containing correlation table
#' @export
#'
#' @examples
#' #NULL
corrStars <- function(x, method = "spearman", dec = 2, N = FALSE) {
  x <- as.matrix(x)
  R <- rcorr(x,type = method)$r
  p <- rcorr(x,type = method)$P
  N <- rcorr(x, type = method)$n
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
  #Add pairwise Ns if N = TRUE
  if(N == TRUE) {
    Rnew <- matrix(paste(N, Rnew, sep = ", "), ncol = ncol(Rnew))
  }
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- Rnew[, -length(colnames(x))]
  Rnew <- Rnew[-1,]
  colnames(Rnew) <- colnames(x)[1:(length(colnames(x))-1)]
  rownames(Rnew) <- colnames(x)[2:length(colnames(x))]
  print("Note: ***: p \u2264 0.001; ** : p \u2264 0.01; *: p \u2264 0.05, no asterisk: p > 0.05")
  return(Rnew)
}
