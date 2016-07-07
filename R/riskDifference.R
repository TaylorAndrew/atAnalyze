#' riskDifference
#'
#' riskDifference provides the risk difference for a categorical predictor and binary outcome. The reference/control group is the first level of the predictor variable.
#'
#' @param data data.frame containing x and y
#' @param y outcome variable name
#' @param x categorical predictor variable name
#'
#' @return A data.frame object containing relative risks with 95% CI
#' @export
#'
#' @examples
#' #test_df <- data.frame(x=sample(c(1,2,3), 50, replace = T),
#' #                      y=sample(c(0,1), 50, replace = T))
#' #riskDifference(test_df, x = 'x', y = 'y')
riskDifference <- function(data, y, x) {
  tb <- table(data[, y], data[, x])
  doOne <- function(i) {
    CE = tb[2,1]
    CN = tb[1,1]
    EE = tb[2, i+1]
    EN = tb[1, i+1]
    ES = EE + EN
    EER = EE/ES
    CS = CE + CN
    CER = CE/CS

    RD = EER - CER
    Int = 1.96 * sqrt(((EER*(1-EER))/EN) + ((CER*(1-CER))/CN))
    RDci = round(100*(RD + c(-Int, Int)))
    RD = round(100*RD)
    return(paste0(RD, "% (", RDci[1], "%, ", RDci[2], "%)"))
  }
  out <- do.call(rbind, lapply(1:(length(levels(factor(data[, x])))-1), doOne))
  row.names(out) <- levels(factor(data[, x]))[-1]
  colnames(out) <- "RD (95%)"
  return(out)
}