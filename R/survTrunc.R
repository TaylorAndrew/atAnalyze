#' survTrunc
#'
#' survTrunc iteratively truncates patient observations at specified time intervals and provides Kaplan-Meier survival curves analysis at each.
#'
#' @param data data.frame containing all variables
#' @param time variable name for the time variable
#' @param event variable name for the event variable
#' @param group variable name for the group variable
#' @param times vector of times to truncate and test survival curves at
#' @param grpLabels labels for the group variable
#' @param timeLabels labels for the time variable
#'
#' @return data.frame containing results
#' @export
#'
#' @examples
#' #NULL
survTrunc <- function(data, time, event, group, times, grpLabels=NULL, timeLabels=NULL) {
  doOne <- function(OneTime) {
    if(!is.null(timeLabels)) {
      if(length(timeLabels)!=length(times)){
        return(print(paste0("The length of the time labels: [",
                            length(timeLabels),
                            "] does not match the length of times: [",
                            length(times),
                            "]")))
      }
    }
    df <- data
    df[, event] <- ifelse(df[, time]>= OneTime, 0, df[, event])
    df[, time] <- ifelse(df[, time]>= OneTime, OneTime, df[, time])
    mod <- survfit(Surv(df[, time], df[, event])~factor(df[, group]))
    pmod <- survdiff(Surv(df[, time], df[, event])~factor(df[, group]))
    out <- summary(mod, times = OneTime)
    if(!is.null(grpLabels)) {
       if(length(grpLabels)!=length(rownames(out$table))){
         return(print(paste0("The length of the group labels: [",
                             length(grpLabels),
                             "] does not match the length of times: [",
                             length(rownames(out$table)),
                             "]")))
       }
     }
    mat <- matrix(
      paste0(
        round(100 * out$surv, 2),
        "% (",
        round(100 * out$lower, 2),
        "%, ",
        round(100 * out$upper, 2),
        "%)"
      ),
    ncol = if (is.vector(out$table)) {
      1
    } else {
      length(out$table[,1])
    },
    nrow = 1,
    byrow = F
  )
    modOverall <- survfit(Surv(df[, time], df[, event])~1)
    outOv <- summary(modOverall, times = OneTime)
    matOv <- matrix(
      paste0(
        round(100 * outOv$surv, 2),
        "% (",
        round(100 * outOv$lower, 2),
        "%, ",
        round(100 * outOv$upper, 2),
        "%)"
      ),
    ncol = if (is.vector(outOv$table)) {
      1
    } else {
      length(outOv$table[,1])
    },
    nrow = 1,
    byrow = F
  )
    Line <- data.frame(OneTime, matOv, mat, kmPval(pmod))
    names(Line) <- c("Time", "Overall", if(!is.null(grpLabels)){grpLabels} else {
      gsub("=",
             "",
             gsub("^[^`=]+",
                  "",
                  rownames(out$table),
                  perl = T))}, "Pval")
  return(Line)
  }
  Output <- do.call(rbind, lapply(times, doOne))
  if(!is.null(timeLabels)) Output[,1] <- timeLabels
  return(Output)
}