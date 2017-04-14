#' survPercent
#'
#' survPercent provides survival percentages for a survfit(Surv()) model for specified times.
#'
#' @param model A survfit(Surv()) model from the survival package
#' @param times A vector of times to provide survival percentages for
#' @param labels Labels for the grouping variable, if applicable
#' @param variance Must be either 'confint' or 'stderr'. Decides whether to provide survival as a percent with 95% confidence intervals or as an estimate with a standard error.
#'
#' @return A data.frame containing survival percentages.
#' @export
#'
#' @examples
#' #example_dat <- data.frame(time = runif(100, 2, 100),
#' #                          event = sample(c(0, 1), 100, prob = c(.2, .8), replace = T),
#' #                          grp = sample(letters[1:3], 100, replace = T))
#' #library(survival)
#' #sMod <- survfit(Surv(time, event) ~ grp, data = example_dat)
#' #survPercent(model = sMod, times = c(10, 20, 30, 40, 50))
#' #survPercent(model = sMod, times = 1:5*10, labels=c("Group 1", "Group 2", "Group 3"))
survPercent <- function(model,
                        times,
                        labels = NULL,
                        variance = 'confint') {
  if(!variance%in%c('confint', 'stderr')) {
      return(print("variance must be either 'confint' or 'stderr'"))
  }
  out <- summary(model, times = times)
  if(variance == 'confint'){
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
    nrow = length(times),
    byrow = F
  )
  } else {
  mat <- matrix(
    paste0(
      round(out$surv, 3),
      " (",
      round(out$std.err, 3),
      ")"
    ),
    ncol = if (is.vector(out$table)) {
      1
    } else {
      length(out$table[,1])
    },
    nrow = length(times),
    byrow = F
  )
  }
  m <- data.frame(Time = times, data.frame(mat))
  if (is.null(labels)) {
    names(m) <-
      c("Time", if (is.null(rownames(out$table))) {
        "Overall"
      } else {
        gsub("=",
             "",
             gsub("^[^`=]+",
                  "",
                  rownames(out$table),
                  perl = T))
      })
  } else {
    names(m) <- c("Time", labels)
  }
    return(m)
}
