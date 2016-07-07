#' multCox
#'
#' multCox will compute univariate Cox proportional hazards model models for multiple predictor variables.
#'
#' @param data Input data.frame
#' @param timevar Time variable
#' @param event Event variable, where event=1, censor=0
#' @param vars Vector of variable names to be modeled
#' @param classvars Vector of variable names that should be treated as categorical
#' @param varLabelTable A varLabelTable object that will convert variable names to nicer labels in the output table
#' @param dec Number of decimals for descriptive statistics
#' @param pdec Number of decimals for p-values
#'
#' @return Output table including the results from each univariate model
#' @export
#'
#' @examples
#' #NULL
multCox <- function(data,
                    timevar,
                    event,
                    vars,
                    classvars = NULL,
                    varLabelTable = NULL,
                    dec = 2,
                    pdec = 3) {
  data <- as.data.frame(data)

  if (is.null(classvars)) contvars <- vars
  if (!is.null(classvars))
  contvars <- vars[!(vars %in% classvars)]
  get <- function(x) {
    getline <- function(x) {
      data[,timevar] <- as.numeric(as.character(data[,timevar]))
      data[,event] <- as.numeric(as.character(data[,event]))
      data[,x] <- as.numeric(as.character(data[,x]))
      tC <-
        tryCatch(
          coxph(Surv(data[,timevar], data[,event]) ~ data[,x],
                data = data),
          warning = function(w)
            w,
          error = function(e)
            e
        )
      tC2 <-
        tryCatch(
          summary(coxph(Surv(data[,timevar], data[,event]) ~ data[,x],
                        data = data)),
          warning = function(w)
            w,
          error = function(e)
            e
        )
      if (is(tC, "warning") | is(tC2, "warning")) {
        print(paste0(x, ": Nonconvergence of likelihood function"))
      }
      if (is(tC, "warning") | is(tC2, "warning")) {
        data$timeVariable <- as.numeric(as.character(data[, timevar]))
        data$censorVariable <-
          as.numeric(as.character(data[, event]))
        data$xxxxx <- as.numeric(as.character(data[,x]))

        dataTemp <- data[complete.cases(data[, c("timeVariable",
                                                 "censorVariable",
                                                 "xxxxx")]) ,]
        print(paste0(x, ": Using Firth Penalized Maximum Likelihood Method"))
        coxmodel <- coxphf(Surv(timeVariable, censorVariable) ~
                             xxxxx,
                           data = dataTemp)
        LRmethod = "pML"
      } else {
        print(paste0(x, ": Using ML Method"))
        coxmodel <-
          coxph(Surv(
            as.numeric(as.character(data[,timevar])),
            as.numeric(as.character(data[,event]))
          ) ~
            as.numeric(as.character(data[,x])))
        LRmethod = "ML"
      }
      MeanDV <- mean(as.numeric(as.character(data[,x])), na.rm = TRUE)
      SDDV <- sd(as.numeric(as.character(data[,x])), na.rm = TRUE)
      if (is(tC, "warning") | is(tC2, "warning")) {
        sumcox <- coxmodel
      } else {
        sumcox <- summary(coxmodel)
      }
      Line <- c(
        x,
        "N  Mean (SD)",
        paste0(
          coxmodel$n,"  ",
          round(MeanDV,digits = dec),
          " (",
          round(SDDV,digits = dec),
          ")"
        ),
        if (is(tC, "warning") | is(tC2, "warning")) {
          length(data$censorVariable[data$censorVariable == 1])
        } else {
          coxmodel$nevent
        },
        paste(
          round(if (is(tC, "warning") | is(tC2, "warning")) {
            exp(sumcox$coefficients)
          } else {
            exp(coxmodel$coef)
          },digits = dec),
          paste0(
            "(",
            round(if (is(tC, "warning") |
                      is(tC2, "warning")) {
              sumcox$ci.lower
            } else {
              sumcox$conf.int[,3]
            },digits = 2),
            " , ",
            round(if (is(tC, "warning") |
                      is(tC2, "warning")) {
              sumcox$ci.upper
            } else {
              sumcox$conf.int[,4]
            },digits = 2),
            ")"
          )
        ),
        round(if (is(tC, "warning") |
                  is(tC2, "warning")) {
          sumcox$prob
        } else {
          sumcox$coefficients[5]
        },digits = pdec),
        LRmethod
      )
      names(Line) <- c(
        "Variable",
        "Response",
        "Descriptives",
        "N Events",
        "Hazard Ratio",
        "p-value",
        "Method"
      )
      return(Line)
    }
    getchunk <- function(x) {
      data[,x] <- factor(data[,x])
      mat <- matrix(nrow = length(levels(data[,x])),ncol = 7,NA)
      mat[,1] <- c(x,rep("",length(levels(data[,x])) - 1))
      mat[,4] <- t(table(data[,event],data[,x]))[,2]
      mat[,3] <- paste0(table(data[,x]),
                        " (",
                        round(100 * prop.table(table(data[,x])), digits = 2),
                        "%)"
                        )
      mat[,2] <- levels(factor(data[,x]))
      tC <- tryCatch(
        coxph(Surv(as.numeric(
          as.character(data[,timevar])
        ),
        data[,event]) ~ data[,x],
        data = data),
        warning = function(w)
          w,
        error = function(e)
          e
      )
      tC2 <-
        tryCatch(
          summary(coxph(
            Surv(as.numeric(as.character(data[,timevar])),
                 data[,event]) ~ data[,x],
            data = data
          )),
          warning = function(w)
            w,
          error = function(e)
            e
        )
      if (is(tC, "warning") | is(tC2, "warning")) {
        print(paste0(x, ": Nonconvergence of likelihood function"))
      }
      if (is(tC, "warning") | is(tC2, "warning")) {
        data$timeVariable <- as.numeric(as.character(data[, timevar]))
        data$censorVariable <-
          as.numeric(as.character(data[, event]))
        data$xxxxx <- factor(data[,x])
        dataTemp <- data[complete.cases(data[, c("timeVariable",
                                                 "censorVariable",
                                                 "xxxxx")]) ,
                         c("timeVariable",
                           "censorVariable",
                           "xxxxx")]
        print(paste0(x, ": Using Firth Penalized Maximum Likelihood Method"))
        coxmodel <-
          coxphf(data = dataTemp, Surv(timeVariable, censorVariable) ~ xxxxx)
        LRmethod = "pML"
      } else {
        print(paste0(x, ": Using ML Method"))
        coxmodel <-
          coxph(Surv(
            as.numeric(as.character(data[,timevar])),
            as.numeric(as.character(data[,event]))
          ) ~ data[,x])
        LRmethod = "ML"
      }
      if (is(tC, "warning") | is(tC2, "warning")) {
        sumcox <- coxmodel
      } else {
        sumcox <- summary(coxmodel)
      }
      mat[,5] <- paste(c("",
                         round(if (is(tC, "warning") |
                                   is(tC2, "warning")) {
                           exp(sumcox$coefficients)
                         } else {
                           sumcox$coef[,2]
                         },
                         digits = dec)),
                       c("",
                         paste0(
                           "(",
                           round(if (is(tC, "warning") |
                                     is(tC2, "warning")) {
                             sumcox$ci.lower
                           } else {
                             sumcox$conf.int[,3]
                           },
                           digits = 2),
                           " , " ,
                           round(if (is(tC, "warning") |
                                     is(tC2, "warning")) {
                             sumcox$ci.upper
                           } else {
                             sumcox$conf.int[,4]
                           },
                           digits = 2),
                           ")"
                         )))
        logp <- if (is(tC, "warning") | is(tC2, "warning")) {
          ""
        } else {
          round(sumcox$logtest[3], digits = pdec)
        }
        otherps <- round(if (is(tC, "warning") | is(tC2, "warning")) {
          sumcox$prob
        } else {
          sumcox$coefficients[,5]
        },digits = pdec)

      pvalues <- c(logp,otherps)
      pvalues <- ifelse(pvalues == 0, "< 0.001",
                        ifelse(pvalues == 1, "> 0.999",
                               pvalues))
      mat[, 6] <- pvalues
      mat[1, 7] <- LRmethod
      names(mat) <- c(
        "Variable",
        "Response",
        "Descriptives",
        "N Events",
        "Hazard Ratio",
        "p-value",
        "Method"
      )
      return(mat)
    }
    if (!x %in% classvars) {
      return(getline(x))
    } else {
      return(getchunk(x))
    }
  }
  finaltable <- do.call(rbind,lapply(vars, get))
  allvars <- vars
  if (!is.null(varLabelTable)) {
    finaltable[,1] <- as.character(finaltable[,1])
    for (i in 1:length(finaltable[,1])) {
      for (j in 1:length(varLabelTable[,1])) {
        if (as.character(finaltable[i,1]) == as.character(varLabelTable[j,1])) {
          finaltable[i,1] <- as.character(varLabelTable[j,2])
        }
      }
    }
  }
finaltable <- as.data.frame(finaltable)
return(finaltable)
}
