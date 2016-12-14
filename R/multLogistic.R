#' multLogistic
#'
#' multLogistic performs univariate, multivariate, and/or conditional logistic regression models over a list of predictor variables.
#'
#' @param data data.frame containing all data for analysis
#' @param y Outcome variable name
#' @param predlist Vector of predictor variable names
#' @param catlist Vector of categorical variable names
#' @param controlCatList Vector of categorical variable names to control for in each model
#' @param controlContinList Vetor of continuous variable names to control for in each model
#' @param conditionalLogistic If TRUE, perform conditional logisitic regression
#' @param matchIDs If conditionalLogistic==TRUE, variable name that contains matchIDs
#' @param varLabelTable varLabelTable object containing variable names in column 1 and variable labels in column 2
#' @param minObsForML The minimum number of events (y==1) observed in the dataset before the method should be switched to Firth's Penalized Maximum Likelihood. Default is 20, which is commonly used as a good rule of thumb.
#' @param dec Number of decimals to round descriptive statistics and odds ratios to
#' @param pdec Number of decimals to round p-values to
#' @param rowwisePercent If TRUE, present row-wise percents for descriptive statistics
#' @param includeCI If TRUE, confidence intervals for the means/proportions are included in the output.
#' @param verbose If TRUE, print additional information to the console
#' @param riskRatio If TRUE, output the risk ratio instead of odds ratios
#'
#' @return An R data.frame object contianing all results.
#' @export
#'
#' @examples
#' #NULL
multLogistic <- function(data,
                         y,
                         predlist,
                         catlist = NULL,
                         controlCatList = NULL,
                         controlContinList = NULL,
                         conditionalLogistic = FALSE,
                         matchIDs = NULL,
                         varLabelTable = NULL,
                         minObsForML = 20,
                         dec = 2,
                         pdec = 3,
                         rowwisePercent = T,
                         includeCI = F,
                         verbose = TRUE,
                         relativeRisk = FALSE) {
  yvar <- y
  data <- as.data.frame(data)
  data <- subset(data, !is.na(data[, yvar]))
  if (conditionalLogistic == FALSE) {
    if (is.null(controlCatList) & is.null(controlContinList)) {
      title <-
        paste0("Table X. Univariate binary logistic regression predicting ",
               yvar)
    }
    if (!is.null(controlCatList) & is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary logistic regression predicting ",
          yvar,
          " controlling for ",
          paste(controlCatList,collapse = ', ')
        )
    }
    if (is.null(controlCatList) & !is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary logistic regression predicting ",
          yvar," controlling for ",
          paste(controlContinList,collapse = ', ')
        )
    }
    if (!is.null(controlCatList) & !is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary logistic regression predicting ",
          yvar,
          " controlling for ",
          paste(c(controlCatList,controlContinList),collapse = ', ')
        )
    }
  }
  if (conditionalLogistic == TRUE) {
    if (is.null(controlCatList) & is.null(controlContinList)) {
      title <-
        paste0("Table X. Univariate binary conditional logistic regression predicting ",
               y)
    }
    if (!is.null(controlCatList) & is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary conditional logistic regression predicting ",
          y,
          "controlling for ",
          paste(controlCatList, collapse = ', ')
        )
    }
    if (is.null(controlCatList) & !is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary conditional logistic regression predicting ",
          y,
          "controlling for ",
          paste(controlContinList, collapse = ', ')
        )
    }
    if (!is.null(controlCatList) & !is.null(controlContinList)) {
      title <-
        paste0(
          "Table X. Multivariable binary conditional logistic regression predicting ",
          y,
          "controlling for ",
          paste(c(controlCatList,controlCatList), collapse = ', ')
        )
    }
  }
  getall <- function(x) {
    if (verbose == T)
      print(paste0("Current Variable is: ",x))
    if (!is.null(controlCatList) & !is.null(controlContinList)) {
      newdat <- data[,c(x,yvar,controlCatList,controlContinList)]
    }
    if (is.null(controlCatList) & !is.null(controlContinList)) {
      newdat <- data[,c(x,yvar,controlContinList)]
    }
    if (!is.null(controlCatList) & is.null(controlContinList)) {
      newdat <- data[,c(x,yvar,controlCatList)]
    }
    if (is.null(controlCatList) & is.null(controlContinList)) {
      newdat <- data[,c(x,yvar)]
    }
    newdat[,yvar] <- as.numeric(newdat[,yvar])
    if (x %in% catlist)
      newdat[,x] <- factor(newdat[,x])
    if (!(x %in% catlist))
      newdat[,x] <- as.numeric(as.character(newdat[,x]))
    #
    newdat <- newdat[complete.cases(newdat),]
    numEvents <- length(newdat[, yvar][newdat[,yvar]==1 & !is.na(newdat[,yvar])])<minObsForML

    if (conditionalLogistic == FALSE) {
      if (!is.null(controlCatList) & !is.null(controlContinList)) {
        f <- paste("newdat[,yvar] ~ newdat[,x] +",
                   paste(sprintf(
                     "`%s`", c(controlContinList,controlCatList)
                   ),
                   collapse = "+"))
      }
      if (is.null(controlCatList) & !is.null(controlContinList)) {
        f <- paste("newdat[,yvar] ~ newdat[,x] +",
                   paste(sprintf("`%s`", controlContinList), collapse = "+"))
      }
      if (!is.null(controlCatList) & is.null(controlContinList)) {
        f <- paste("newdat[,yvar] ~ newdat[,x] +",
                   paste(sprintf("`%s`", controlCatList), collapse = "+"))
      }
    }
    if (conditionalLogistic == TRUE) {
      if (!is.null(controlCatList) & !is.null(controlContinList)) {
        f <- paste(
          "newdat[,yvar] ~ newdat[,x] +",
          paste(sprintf(
            "`%s`",
            c(controlContinList,controlCatList)
          ),
          collapse = "+"),"+ strata(newdat[,matchIDs])"
        )
      }
      if (is.null(controlCatList) & !is.null(controlContinList)) {
        f <- paste(
          "newdat[,yvar] ~ newdat[,x] +",
          paste(sprintf("`%s`", controlContinList),
                collapse = "+"),
          "+ strata(newdat[,matchIDs])"
        )
      }
      if (!is.null(controlCatList) & is.null(controlContinList)) {
        f <- paste(
          "newdat[,yvar] ~ newdat[,x] +",
          paste(sprintf("`%s`", controlCatList), collapse = "+"),
          "+ strata(newdat[,matchIDs])"
        )
      }
    }
    if (!(x %in% catlist)) {
      ciMean <- function(x, conf.int = .95) {
        mean <- mean(x, na.rm = T)
        sd <- sd(x, na.rm = T)
        n  <- length(x[!is.na(x)])
        SE = sd / sqrt(n)
        E = qt(conf.int + (1 - conf.int) / 2, df = n - 1) * SE
        return(mean + c(-E, E))
      }
      newdat[,x] <- as.numeric(as.character(newdat[,x]))
      N <- length(newdat[!is.na(x),x])
      mean <-
        round(mean(as.numeric(as.character(newdat[,x]))),digits = dec)
      stdev <-
        round(sd(as.numeric(as.character(newdat[,x]))),digits = dec)
      ci = ciMean(newdat[,x])
      d0 <- subset(newdat,newdat[,yvar] == 0)
      N0 <- length(d0[,x])
      mean0 <- round(mean(d0[,x]),digits = dec)
      stdev0 <- round(sd(d0[,x]),digits = dec)
      ci0 <- ciMean(d0[,x])
      d1 <- subset(newdat,newdat[,yvar] == 1)
      N1 <- length(d1[,x])
      mean1 <- round(mean(d1[,x]),digits = dec)
      stdev1 <- round(sd(d1[,x]),digits = dec)
      ci1 <- ciMean(d1[,x])
      if (conditionalLogistic == FALSE) {
        if (is.null(controlCatList) & is.null(controlContinList)) {
          tC <- tryCatch(
            if(relativeRisk==T) {
            glm(
              newdat[,yvar] ~ newdat[,x],
              family = poisson(link = 'log'),
              data = newdat
            )
            } else {
            glm(
              newdat[,yvar] ~ newdat[,x],
              family = binomial(logit),
              data = newdat
            )
            }
            ,
            warning = function(w)
              w,
            error = function(e)
              e
          )
          tC2 <- tryCatch(
            confint(
              if(relativeRisk==T) {
            glm(
              newdat[,yvar] ~ newdat[,x],
              family = poisson(link = 'log'),
              data = newdat
            )
            } else {
            glm(
              newdat[,yvar] ~ newdat[,x],
              family = binomial(logit),
              data = newdat
            )
            }
            ),
            warning = function(w)
              w,
            error = function(e)
              e
          )
          if (is(tC, "warning") |
              is(tC2, "warning"))
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if(numEvents) {
              print(paste0(x, ": Less than 20 events"))
          }
          if (is(tC, "warning") | is(tC2, "warning") | numEvents) {
            print(paste0(
              x, ": Using Firth Penalized Maximum Likelihood Method"
            ))
            mod <- logistf(newdat[,yvar] ~ newdat[,x], data = newdat)
            LRmethod = "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            if(relativeRisk==T) {
            mod <-
              glm(newdat[,yvar] ~ newdat[,x], family = poisson(link = 'log'))
            } else {
            mod <-
              glm(newdat[,yvar] ~ newdat[,x], family = binomial(logit))
            }
              LRmethod = 'ML'
          }
        }
        if (!is.null(controlCatList) & !is.null(controlContinList)) {
          tC <- tryCatch(
            glm(
              newdat[,yvar] ~ factor(newdat[,x]),
              family = binomial(logit)
            ),
            warning = function(w)
              w,
            error = function(e)
              e
          )
          tC2 <- tryCatch(
            confint(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]),
                family = binomial(logit)
              )
            ),
            warning = function(w)
              w,
            error = function(e)
              e
          )
          if (is(tC, "warning") |
              is(tC2, "warning"))
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if(numEvents) {
              print(paste0(x, ": Less than 20 events"))
          }
          if (is(tC, "warning") | is(tC2, "warning") | numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <- logistf(formula(f), family = binomial, data = newdat)
            LRmethod  <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit))
            LRmethod <- "ML"
          }
        }
        if (!is.null(controlCatList) & is.null(controlContinList)) {
          tC <-
            tryCatch(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit), data = newdat
              ),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          tC2 <-
            tryCatch(
              confint(glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit)
              ), data = newdat),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          if (is(tC, "warning") |
              is(tC2, "warning"))
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if(numEvents) {
              print(paste0(x, ": Less than 20 events"))
          }
          if (is(tC, "warning") | is(tC2, "warning") | numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <- logistf(formula(f), family = binomial, data = newdat)
            LRmethod <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit), data = newdat)
            LRmethod <- "ML"
          }
        }
        if (is.null(controlCatList) & !is.null(controlContinList)) {
          tC <- tryCatch(
            glm(
              newdat[,yvar] ~ factor(newdat[,x]),
              family = binomial(logit),
              data = newdat
            ),
            warning = function(w)
              w,
            error = function(e)
              e
          )
          tC2 <- tryCatch(
            confint(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]),
                family = binomial(logit),
                data = newdat
              )
            ),
            warning = function(w)
              w,
            error = function(e)
              e
          )
          if (is(tC, "warning") |
              is(tC2, "warning"))
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if(numEvents) {
             print(paste0(x, ": Less than 20 total events"))
          }
          if (is(tC, "warning") | is(tC2, "warning") |
              numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <- logistf(formula(f), family = binomial, data = newdat)
            LRmethod <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit), data = newdat)
            LRmethod <- "ML"
          }
        }

        if (is(tC, "warning") | is(tC2, "warning") | numEvents)
          model = mod
        if (!(is(tC, "warning") |
              is(tC2, "warning")| numEvents))
          model <- summary(mod)
        LogisticOR <- function(model) {
          k <-
            as.data.frame(cbind(paste0(
              round(exp(cbind(
                coef(model)
              )),digits = 2),
              " (",
              round(exp(cbind(
                confint(model)[,1]
              )),digits = 2),
              ", ",
              ifelse(exp(cbind(
                confint(model)[,2]
              )) > 100,"Inf",
              round(exp(
                cbind(confint(model)[,2])
              ),digits = 2)),")"
            )))
          rownames(k) <- rownames(cbind(coef(model)))
          return(k)
        }
        OR <- LogisticOR(mod)[2,1]
        if (!(is(tC, "warning") |
              is(tC2, "warning") | numEvents))
          pvalue <- round(model$coefficients[-1,4],digits = pdec)
        if (is(tC, "warning") |
            is(tC2, "warning")|
              numEvents){
          print(model)
          pvalue <- round(model$prob[-1],digits = pdec)
        }

        pvalue <-
          ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
        pvalue <- pvalue[1]
      }
      if (conditionalLogistic == TRUE) {
        if (is.null(controlCatList) & is.null(controlContinList)) {
          mod <- clogit(newdat[,yvar] ~ newdat[,x] + strata(newdat[,matchIDs]),
                        method = "exact", data = newdat)
          LRmethod = "Cond.Logit"
        }
        if (!is.null(controlCatList) & !is.null(controlContinList)) {
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod = "Cond.Logit"
        }
        if (!is.null(controlCatList) & is.null(controlContinList)) {
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod = "Cond.Logit"
        }
        if (is.null(controlCatList) & !is.null(controlContinList)) {
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod = "Cond.Logit"
        }
        model <- summary(mod)
        LogisticOR <- function(model) {
          k <-
            as.data.frame(cbind(paste0(
              round(exp(cbind(
                coef(model)
              )),digits = 2),
              " (",
              round(exp(cbind(
                confint(model)[,1]
              )),digits = 2),
              ", ",
              ifelse(exp(cbind(
                confint(model)[,2]
              )) > 100,"Inf",
              round(exp(
                cbind(confint(model)[,2])
              ),digits = 2)),")"
            )))
          rownames(k) <- rownames(cbind(coef(model)))
          return(k)
        }
        OR <- LogisticOR(mod)[1,1]
        pvalue <- round(model$coefficients[1,5],digits = pdec)
        pvalue <-
          ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
        pvalue <- pvalue[1]
      }
    if(includeCI==F) {
      line <- c(
        x,"N  Mean (Std.Dev)",
        paste0(N,"  ",mean," (",stdev, ")"),
        paste0(N0,"  ",mean0," (",stdev0, ")"),
        paste0(N1,"  ",mean1," (",stdev1, ")"),
        as.character(OR),pvalue, LRmethod
      )
      names(line) <- c(
        "Variable",
        "Response",
        "Overall",
        paste0(yvar," = 0"),
        paste0(yvar," = 1"),
        if(relativeRisk==T) "RR (95% CI)" else "OR (95% CI)",
        "p-value",
        "Method"
      )
    }
       if(includeCI==T) {
      line <- c(
        x,"N  Mean (Std.Dev)\n (95% CI)",
        paste0(N,"  ",mean," (",stdev, ")\n", paste0("(", round(ci[1], dec), ", ",  round(ci[2], dec), ")")),
        paste0(N0,"  ",mean0," (",stdev0, ")\n", paste0("(",  round(ci0[1], dec), ", ",  round(ci0[2], dec), ")")),
        paste0(N1,"  ",mean1," (",stdev1, ")\n", paste0("(",  round(ci1[1], dec), ", ",  round(ci1[2], dec), ")")),
        as.character(OR),pvalue, LRmethod
      )
      names(line) <- c(
        "Variable",
        "Response",
        "Overall",
        paste0(yvar," = 0"),
        paste0(yvar," = 1"),
        if(relativeRisk==T) "RR (95% CI)" else "OR (95% CI)",
        "p-value",
        "Method"
      )
    }
      return(line)
    }

    if (x %in% catlist) {
      ciProp <- function(x, conf.int = .95) {
        x <- x[!is.na(x)]
        n <- length(x)
        pbar =  unlist(table(x)) / n
        SE = sqrt(pbar * (1 - pbar) / n)
        E = qnorm(conf.int + (1 - conf.int) / 2) * SE
        out = do.call(rbind, lapply(1:length(pbar), function(i)
          pbar[i] + c(-E[i], E[i])))
        row.names(out) <- names(pbar)
        colnames(out) <- c("ci.lower", "ci.upper")
        return(out)
      }

      if(rowwisePercent==T) {
        levsci <- levels(factor(newdat[,x]))
        getCorrect <- function(lev) {
          ll <- newdat[factor(newdat[,x])==lev,]
        if(length(ciProp(ll[, yvar])[,1])==2){
          return(data.frame(yvar0 = paste0("(",
                                              round(100*ciProp(ll[, yvar])[1,][1]),
                                              "%, ",
                                              round(100*ciProp(ll[, yvar])[1,][2]),
                                              "%)"),
                            yvar1 = paste0("(",
                                              round(100*ciProp(ll[, yvar])[2,][1]),
                                              "%, ",
                                              round(100*ciProp(ll[, yvar])[2,][2]),
                                              "%)")
          )) } else {return(data.frame(yvar0= "(NA, NA)", yvar1="(NA, NA)"))}
        }
        ciDF <- do.call(rbind, lapply(levsci, getCorrect))
      }
      if(rowwisePercent==F) {
         levsci <- levels(factor(newdat[,yvar]))
        getCorrect <- function(lev) {
          ll <- newdat[factor(newdat[,yvar])==lev,]
          return(data.frame(yvar0 = paste0("(",
                                              round(100*ciProp(ll[, x])[,1]),
                                              "%, ",
                                              round(100*ciProp(ll[, x])[,2]),
                                              "%)")

          ))
        }
        ciDF <- do.call(cbind, lapply(levsci, getCorrect))
      }
      len <- length(levels(factor(newdat[,x])))
      tab <- t(table(newdat[,x]))
      proptab <- round(t(100 * prop.table(tab)),digits = dec)

      mat <- matrix(ncol = 8,nrow = len)
      mat[1,1] <- x
      mat[,2] <- levels(factor(newdat[,x]))
      mat[,3] <-
        paste0(tab,"  (",if (rowwisePercent == T) {
          100
        }else{
          proptab
        },"%)")
      fulltab <- table(newdat[,x],newdat[,yvar])
      proptab <-
        round(100 * prop.table(fulltab,if (rowwisePercent == T) {
          1
        }else{
          2
        }),digits = dec)
      if(includeCI == F) {
        mat[,4] <- paste0(fulltab[,1],"  (",proptab[,1],"%)")
        mat[,5] <- paste0(fulltab[,2],"  (",proptab[,2],"%)")
      }
      if(includeCI == T) {
        mat[,4] <- paste0(fulltab[,1],"  (",proptab[,1],"%)\n", ciDF[,1])
        mat[,5] <- paste0(fulltab[,2],"  (",proptab[,2],"%)\n", ciDF[,2])
      }
       numEvents <- length(newdat[, yvar][newdat[,yvar]==1 & !is.na(newdat[,yvar])])<minObsForML
      if (conditionalLogistic == FALSE) {
        if (is.null(controlCatList) & is.null(controlContinList)) {
          tC <-
            tryCatch(
            if(relativeRisk==T) {
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = poisson(link='log'), data = newdat
              )
            } else {
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit), data = newdat
              )
            },
              warning = function(w)
                w,
              error = function(e)
                e
            )
          tC2 <-
            tryCatch(
              confint(
                if(relativeRisk==T) {
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = poisson(link = 'log')
              , data = newdat)} else {
                glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit)
              , data = newdat)
                  }
              )
              ,
              warning = function(w)
                w,
              error = function(e)
                e
            )

          if (is(tC, "warning") |
              is(tC2, "warning")|
              numEvents)
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if (is(tC, "warning") | is(tC2, "warning")|
              numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <- logistf(newdat[,yvar] ~ newdat[,x], data = newdat)
            LRmethod <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            if(relativeRisk==T) {
                mod <-
              glm(newdat[,yvar] ~ newdat[,x], family = poisson(link = 'log'))
            } else {
                mod <-
              glm(newdat[,yvar] ~ newdat[,x], family = binomial(logit))
            }

            LRmethod <- "ML"
          }
        }
        if (!is.null(controlCatList) & !is.null(controlContinList)) {
          tC <-
            tryCatch(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit), data = newdat
              ),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          tC2 <-
            tryCatch(
              confint(glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit)
              ), data = newdat),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          if (is(tC, "warning") |
              is(tC2, "warning")|
              numEvents)
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if (is(tC, "warning") | is(tC2, "warning")|
              numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <- logistf(formula(f), family = binomial, data = newdat)
            LRmethod = "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit), data = newdat)
            LRmethod = "ML"
          }
        }
        if (!is.null(controlCatList) & is.null(controlContinList)) {
          tC <-
            tryCatch(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit), data = newdat
              ),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          tC2 <-
            tryCatch(
              confint(glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit)
              ), data = newdat),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          if (is(tC, "warning") |
              is(tC2, "warning")|
              numEvents)
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if (is(tC, "warning") | is(tC2, "warning")|
              numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <-
              logistf(formula(f), family = binomial, data = newdat)
            LRmethod <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit), data = newdat)
            LRmethod <- "ML"
          }
        }
        if (is.null(controlCatList) & !is.null(controlContinList)) {
          tC <-
            tryCatch(
              glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit), data = newdat
              ),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          tC2 <-
            tryCatch(
              confint(glm(
                newdat[,yvar] ~ factor(newdat[,x]), family = binomial(logit)
              ), data = newdat),
              warning = function(w)
                w,
              error = function(e)
                e
            )
          if (is(tC, "warning") |
              is(tC2, "warning")|
              numEvents)
            print(paste0(x, ": Quasi- or Complete-Separation"))
          if (is(tC, "warning") | is(tC2, "warning")|
              numEvents) {
            print(paste0(x, ": Using Firth Method"))
            mod <-
              logistf(formula(f), family = binomial, data = newdat)
            LRmethod <- "pML"
          } else {
            print(paste0(x, ": Using ML Method"))
            mod <-
              glm(formula(f), family = binomial(logit), data = newdat)
            LRmethod <- "ML"
          }
        }
      }
      if (conditionalLogistic == TRUE) {
        if (is.null(controlCatList) &
            is.null(controlContinList))
          mod <-
            clogit(newdat[,yvar] ~ newdat[,x], method = "exact", data = newdat)
          LRmethod <- "Conditional Logistic"
        if (!is.null(controlCatList) &
            !is.null(controlContinList))
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod <- "Conditional Logistic"
        if (!is.null(controlCatList) &
            is.null(controlContinList))
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod <- "Conditional Logistic"
        if (is.null(controlCatList) &
            !is.null(controlContinList))
          mod <- clogit(formula(f), method = "exact", data = newdat)
          LRmethod <- "Conditional Logistic"
      }
      if (conditionalLogistic == FALSE) {
        if (is(tC, "warning") | is(tC2, "warning")|
              numEvents)
          model = mod
        if (!(is(tC, "warning") |
              is(tC2, "warning")))
          model <- summary(mod)
        if (!(is(tC, "warning")) & !(is(tC2, "warning")) &
            !numEvents) {
          LogisticOR <- function(mod) {
            model <- summary(mod)
            k <- as.data.frame(cbind(paste0(
              ifelse(
                round(exp(coef(model)[,1]),digits = 2) > 10000,
                "NA",
                round(exp(coef(model)[,1]),digits = 2)
              ),
              " (",
              round(exp(confint(mod)[,1]),digits = 2),
              ", ",
              ifelse(exp(confint(mod)[,2]) > 10000,"Inf",
                     round(exp(
                       confint(mod)[,2]
                     ),digits = 2)),
              ")"
            )))
            rownames(k) <- rownames(cbind(coef(model)))
            return(k)
          }
        } else {
          LogisticOR <- function(mod) {
            model <- mod
            k <- as.data.frame(cbind(paste0(
              ifelse(
                round(exp(mod$coefficients),digits = 2) > 10000,
                "NA",
                round(exp(mod$coefficients),digits = 2)
              ),
              " (",
              round(exp(confint(mod)[,1]),digits = 2),
              ", ",
              ifelse(exp(confint(mod)[,2]) > 10000,"Inf",
                     round(exp(
                       confint(mod)[,2]
                     ),digits = 2)),
              ")"
            )))
            rownames(k) <- rownames(cbind(coef(model)))
            return(k)
          }
        }
        if (!(is(tC, "warning") |
              is(tC2, "warning") |
            numEvents))
          pvalue <- round(model$coefficients[-1,4],digits = pdec)
        if (is(tC, "warning") |
            is(tC2, "warning")|
              numEvents)
          pvalue <- round(model$prob[-1],digits = pdec)
        pvalue <-
          ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
        OR <- LogisticOR(mod = mod)[2:len,1]

        pvalue <- pvalue[1:(len - 1)]

        mat[,6] <- c("",as.character(OR))

        pvalue <- c("Ref",pvalue)

        mat[, 7] <- pvalue
        mat[1, 8] <- LRmethod
        mat[is.na(mat)] <- " "
      }
      if (conditionalLogistic == TRUE) {
        LogisticOR <- function(model) {
          k <- as.data.frame(paste0(
            round(model$coefficients[,2],dec),
            " (",round(model$conf.int[,3],dec),
            ", ",
            round(model$conf.int[,4],dec),")"
          ),
          stringsAsFactors = FALSE)
          rownames(k) <- rownames(cbind(coef(model)))
          names(k) <- "OddsRatios"
          return(k)
        }
        pvalue <- round(model$coefficients[-1,5],digits = pdec)
        pvalue <-
          ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
        pvalue <- pvalue[1:(len - 1)]
        OR <- LogisticOR(model)[2:len,1]
        mat[,6] <- c("",as.character(OR))
        pvalue <- c("Ref",pvalue)
        mat[,7] <- pvalue
        mat[is.na(mat)] <- " "
      }
      line <- mat
      return(line)
    }
    return(line)
  }
    results <- do.call(rbind, lapply(predlist,getall))

  if (!is.null(varLabelTable)) {
    if (length(varLabelTable[1,]) != 2)
      print("varLabelTable does not have exactly 2 columns")
    if (length(varLabelTable[1,]) == 2) {
      names(varLabelTable) <- c("Variable","Variable Name")
      results <- data.frame(id = 1:length(results[,1]),results)
      results <- merge(varLabelTable,results,by = "Variable", all.y = TRUE)
      results <- results[order(results$id),]
      results$Variable <-
        ifelse(!is.na(results$`Variable Name`),as.character(results$`Variable Name`),"")
      results$`Variable Name` <- NULL
      results$id <- NULL
    }
  }
  results <- data.frame(results)
  names(results) <- c(
    "Variable",
    "Response",
    "Overall",
    paste0(yvar, c(" = 0", " = 1")),
    if(relativeRisk==T) "RR (95% CI)" else "OR (95% CI)",
    "p-value",
    "Method"
  )
 return(results)
}
