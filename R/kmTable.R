#' kmTable
#'
#' kmTable provides summary information from a kaplan-meier survival model, including median survival, event counts, and pairwise comparisons.
#'
#' @param data data.frame containing all data
#' @param time Variable name contaning time data
#' @param event Variable name containing event data (event = 1, censor = 0)
#' @param group Variable name containing groupiong variable
#'
#' @return R data.frame containing all km info
#' @export
#'
#' @examples
#' #example_df <- data.frame(a = sample(c(0, 1), 100, replace = T),
#' #                         b = sample(letters[24:26], 100, replace = T),
#' #                         c = rnorm(100))
#' #kmTable(data=example_df, time = "c", event = "a", group = "b")
kmTable <- function(data, time, event, group) {
  fit <- survfit(Surv(data[, time], data[, event]) ~ data[, group])
  diffP <- round(1 - pchisq(survdiff(Surv(data[, time], data[, 
        event]) ~ data[, group])$chisq, length(survdiff(Surv(data[, time], data[, 
        event]) ~ data[, group])$obs) - 1), digits = 3)


  descript <- data.frame(
    Group = gsub("=","",gsub(
      "^[^\\=]+","",
      names(fit$strata), perl =
        T
    )),
    summary(fit)$table[,c(3,4)],
    `Median Survival` = paste0(
      round(summary(fit)$table[,c(7)], digits = 2), " (",
      round(summary(fit)$table[,c(8)], digits = 2), ", ",
      round(summary(fit)$table[,c(9)], digits = 2), ")"
    )
  )

  SurvPairwise <-
    function(data,time,event,group) {
      levs <- levels(factor(data[,group]))
      len <- length(levs)
      chisq <- matrix(0., len,len)
      for (i in 1:len) {
        for (j in (1:len)[-i]) {
          temp <- survdiff(Surv(data[, time], data[, event]) ~ data[, group],
                           subset = (data[, group] %in% (levs)[c(i,j)]))
          chisq[i,j] <- 1 - pchisq(temp$chisq, length(temp$n) - 1)
        }
      }

      pvalList <-
        round(chisq[lower.tri(chisq, diag = FALSE) == TRUE],digits = 3)


      l <- t(combn(levs,2))
      lab <- paste(l[,1], " vs. ", l[,2])
      out <- cbind(lab,pvalList)
      return(out)
    }
  pairwise <-
    SurvPairwise(
      data = data,time = time,event = event,group = group)
  mLen <- max(length(descript[,1]), length(pairwise[,1]))

  outputAsMatrix <- matrix(nrow = mLen, ncol = 7)
  outputAsMatrix[c(1:length(descript[,1])),1] <- descript[,1]
  outputAsMatrix[c(1:length(descript[,1])),2] <- descript[,2]
  outputAsMatrix[c(1:length(descript[,1])),3] <- descript[,3]
  outputAsMatrix[c(1:length(descript[,1])),4] <-
    as.character(descript[,4])
  outputAsMatrix[c(1:length(pairwise[,1])),5] <- pairwise[,1]
  outputAsMatrix[c(1:length(pairwise[,1])),6] <- pairwise[,2]
  outputAsMatrix[1, 7] <- diffP
  output <- as.data.frame(outputAsMatrix)
  names(output) <-
    c(
      "Group", "N Total", "N Events" ,"Median Survival (95% CI)",
      "Pairwise Comparison", "Pairwise P-value", "Overall Log-Rank Test"
    )
  return(output)
}
