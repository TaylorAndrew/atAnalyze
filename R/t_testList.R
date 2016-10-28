#' t_testList
#'
#' t_testList performs paired and unpaired parametric and non-parametric t.tests between columns within a data.frame, for multiple pairs of varaibles given in list form.
#'
#' @param data data.frame containing all columns provided in `list_of_t.tests`
#' @param list_of_t.tests A list, in which each element contains the names of two variables from `data`
#' @param dec Number of decimals to round descriptive statistics to
#' @param pdec Number of decimals to round p-values to
#' @param var.equal If TRUE, and parametric==TRUE, student's t.test will be used, else Welch's t.test methods are used
#' @param parametric If TRUE, t.test() methods are used, else wilcox.test() methods are used
#' @param paired If TRUE, the two variables will be treated as paired.
#' @param NPDescriptives If TRUE, Median and either IQR or range are provided, else Mean and Std are provided
#' @param NPtype Can be 'IQR' or 'range'. If NPDescriptives==TRUE, then the specified variance estimate(s) will be provided
#'
#' @return A data.frame containing the results from the t.tests along with descriptive statistics.
#' @export
#'
#' @examples
#' #example_df <- data.frame(x1 = rnorm(100),
#' #                         x2 = rnorm(100),
#' #                         x3 = rnorm(100),
#' #                         x4 = rnorm(100),
#' #                         x5 = rnorm(100))
#' #list <- list(c("x1", "x2"), c("x1", "x3"), c("x4", "x5"))
#' #t_testList(example_df, list)
#' #t_testList(example_df, list, parametric = F)
#' #t_testList(example_df, list, paired = T)
#' #t_testList(example_df, list, NPDescriptives = T)
#' #t_testList(example_df, list, NPDescriptives = T, NPtype = "range")
t_testList <- function(data,
                       list_of_t.tests,
                       dec = 2,
                       pdec = 3,
                       var.equal = TRUE,
                       parametric = TRUE,
                       paired = FALSE,
                       NPDescriptives = F,
                       NPtype = "IQR",
                       verbose = FALSE) {
  data <- as.data.frame(data)
  do_one <- function(i) {
    varnames <- unlist(list_of_t.tests[i])
    if(verbose==TRUE) print(paste0('Now analyzing ', varnames[1], ' and ', varnames[2]))
    subdat <- data[, unlist(list_of_t.tests[i])]
    if (paired == TRUE)
      subdat <- subdat[complete.cases(subdat), ]
    subdat[, 1] <- as.numeric(as.character(subdat[, 1]))
    subdat[, 2] <- as.numeric(as.character(subdat[, 2]))
    if (NPDescriptives == T) {
      m1 <- round(median(subdat[, 1], na.rm = TRUE), digits = dec)
      m2 <- round(median(subdat[, 2], na.rm = TRUE), digits = dec)
      if (NPtype == "range") {
        sd1 <-
          paste0(round(range(subdat[, 1], na.rm = TRUE), digits = dec)[1],
                 ", ",
                 round(range(subdat[, 1], na.rm = TRUE), digits = dec)[2])
        sd2 <-
          paste0(round(range(subdat[, 2], na.rm = TRUE), digits = dec)[1],
                 ", ",
                 round(range(subdat[, 2], na.rm = TRUE), digits = dec)[2])
      }
      if (NPtype == "IQR") {
        sd1 <- round(IQR(subdat[, 1], na.rm = TRUE), digits = dec)
        sd2 <- round(IQR(subdat[, 2], na.rm = TRUE), digits = dec)
      }
    }
    if (NPDescriptives == F) {
      m1 <- round(mean(subdat[, 1], na.rm = TRUE), digits = dec)
      m2 <- round(mean(subdat[, 2], na.rm = TRUE), digits = dec)
      sd1 <- round(sd(subdat[, 1], na.rm = TRUE), digits = dec)
      sd2 <- round(sd(subdat[, 2], na.rm = TRUE), digits = dec)
    }
    n1 <- length(subdat[, 1][!is.na(subdat[, 1])])
    n2 <- length(subdat[, 2][!is.na(subdat[, 2])])
    des1 <- paste0(n1, "  ", m1, " (", sd1, ")")
    des2 <- paste0(n2, "  ", m2, " (", sd2, ")")
    if (parametric == TRUE) {
      pvalue <- sprintf(paste0("%.", pdec, "f"),
                        round(
                          t.test(
                            subdat[, 1],
                            subdat[, 2],
                            equal.var = equal.var,
                            na.omit = TRUE,
                            paired = paired
                          )$p.value,
                          digits = pdec
                        ))
    }
    if (parametric == FALSE) {
      pvalue <- sprintf(paste0("%.", pdec, "f"),
                        round(
                          wilcox.test(subdat[, 1], subdat[, 2],
                                      na.omit = TRUE,
                                      paired = paired)$p.value,
                          digits = pdec
                        ))
    }
    pvalue <-
      gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
    pvalue <- gsub("^1\\.", "> 0.", pvalue)
    pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
    output <- data.frame(
      Group1 = varnames[1],
      Group1_desc = des1,
      Group2 = varnames[2],
      Group2_desc = des2,
      p_value = pvalue
    )
    return(output)
  }
  output_full <-
    do.call(rbind, lapply(1:length(list_of_t.tests), do_one))
  return(output_full)
}
