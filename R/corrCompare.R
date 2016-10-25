#' corrCompare
#'
#' corrCompare computes a correlation table for each level of a binary grouping variable separately, and then tests differences in correlation coefficients between grouping variables
#'
#' @param data data.frame containing 1 grouping variable and all variables to be used to compute a correlation table
#' @param group variable name of the grouping variable
#' @param method correlation method: pearson or spearman
#' @param rdec number of decimals for correlation coefficients
#' @param pdec number of decimals for p-values
#' @param tdec number of decimals for t-values
#'
#' @return A data.frame containing group level correlation coefficients as well as relevant tests between groups.
#' @export
#'
#' @examples
#' #dat<-data.frame(g<-rep(c(1,2),50),
#' #                v2<-rnorm(100,5,1),
#' #                v1<-rnorm(100,2,1),
#' #                v3<-rnorm(100,499,12))
#' #names(dat) <- c("g","v1","v2","v3")
#' #jf <- corrCompare(data=dat, group="g")
#' #View(jf)
corrCompare <-
  function(data = NULL,
           group = NULL,
           method = "pearson",
           rdec = 3,
           pdec = 3,
           tdec = 3) {
    data <- as.data.frame(data)
    glevels <- sort(unique(data[, group]))
    if (length(glevels) != 2)
      print("Sorry, Group must be two levels")
    else {
      g1 <- subset(data, data[, group] == glevels[1])
      g2 <- subset(data, data[, group] == glevels[2])
      g1[, group] <- NULL
      g2[, group] <- NULL
      g1cor <-
        corr.test(g1,
                  use = "pairwise",
                  method = method,
                  adjust = "none")
      g2cor <-
        corr.test(g2,
                  use = "pairwise",
                  method = method,
                  adjust = "none")
      colNames <- names(g1)
      name <- t(combn(colNames, 2))
      rg1 <-
        round(as.vector(g1cor$r[lower.tri(g1cor$r, diag = FALSE) == TRUE]),
              rdec)
      ng1 <-
        as.vector(g1cor$n[lower.tri(g1cor$n, diag = FALSE) == TRUE])
      pg1 <-
        round(as.vector(g1cor$p[lower.tri(g1cor$p, diag = FALSE) == TRUE]),
              pdec)
      rg2 <-
        round(as.vector(g2cor$r[lower.tri(g2cor$r, diag = FALSE) == TRUE]),
              rdec)
      ng2 <-
        as.vector(g2cor$n[lower.tri(g2cor$n, diag = FALSE) == TRUE])
      pg2 <-
        round(as.vector(g2cor$p[lower.tri(g2cor$p, diag = FALSE) == TRUE]),
              pdec)
      if (length(g2cor$n) == 1)
        ng2 <- rep(g2cor$n, length(name[, 1]))
      if (length(g1cor$n) == 1)
        ng1 <- rep(g1cor$n, length(name[, 1]))
      merged <- data.frame(name, rg1, ng1, pg1,
                           rg2, ng2, pg2)
      merged2 <- subset(merged, merged$ng1 > 3 & merged$ng2 > 3)
      merged2$zg1 <- fisherz(merged2$rg1)
      merged2$zg2 <- fisherz(merged2$rg2)
      getp <- function(i) {
        d <- merged2[i, ]$zg1 - merged2[i, ]$zg2
        sigma <-
          sqrt((1 / (merged2[i, ]$ng1 - 3)) + (1 / (merged2[i, ]$ng2 - 3)))
        teststat <- round(d / sigma, digits = tdec)
        lowerlimit <- round(d - ((1.96) * (sigma)), digits = tdec)
        upperlimit <- round(d + ((1.96) * (sigma)), digits = tdec)
        if (lowerlimit > 0 & upperlimit > 0)
          test = "different"
        else if (lowerlimit < 0 & upperlimit < 0)
          test = "different"
        else
          test = "not different"
        pvalue = 2 * pnorm(-abs(teststat))
        j <-
          c(
            as.character(merged2$X1[i]),
            as.character(merged2$X2[i]),
            merged2[i, 3:8],
            teststat,
            lowerlimit,
            upperlimit,
            test,
            round(pvalue, digits = pdec)
          )
        names(j) <-
          (
            c(
              "Comparison Var1",
              "Comparison Var2",
              "Group1 rcoef",
              "Group1 N",
              "Group1 Pvalue",
              "Group2 rcoef",
              "Group2 N",
              "Group2 Pvalue",
              "teststat",
              "lowerlimit",
              "upperlimit",
              "crosszero.test",
              "pvalue"
            )
          )
        return(j)
      }
      leng <- 1:length(merged2[, 1])
      stats <- data.frame(do.call(rbind, lapply(leng, getp)))
      return(stats)
    }
  }
