#' multReg
#'
#' multReg performs univariate simple linear regression models for a single outcome variable over a list of predictor variables.
#'
#' @param data data.frame containing all variables for analysis
#' @param yvar Variable name of the outcome variable
#' @param predlist Vector of predictor variables
#' @param catlist Vector of variable names found in `predlist` that should be treated as categorical variables
#' @param varLabelTable A varLabelTable object containing variable names in column 1 and variable labels in column 2
#' @param sortVars Vector of variable names in the order to be presented in the output
#' @param dec Number of decimals to round descriptive statistics
#' @param pdec Number of decimals to round p-values
#' @param rsqdec Number of decimals to round r-squared values
#'
#' @return A data.frame object conatining all results.
#' @export
#'
#' @examples
#' #NULL
multReg <- function(data,
                    yvar,
                    predlist,
                    catlist = NULL,
                    varLabelTable = FALSE,
                    sortVars = FALSE,
                    dec = 2,
                    pdec = 3,
                    rsqdec = 4) {

  data <- as.data.frame(data)
  if (!is.null(catlist))
    predlist <- unique(c(predlist,catlist))

  getall <- function(x) {
    newdat <- data[,c(x,yvar)]
    newdat[,yvar] <- as.numeric(newdat[,yvar])
    if (x %in% catlist)
      newdat[,x] <- factor(newdat[,x])
    if (!(x %in% catlist))
      newdat[,x] <- as.numeric(as.character(newdat[,x]))
    newdat <-
      subset(newdat,!is.na(newdat[,yvar]) & !is.na(newdat[,x]))

    if (!(x %in% catlist)) {
      N <- length(newdat[,x])
      mean <- round(mean(newdat[,x]),digits = dec)
      stdev <- round(sd(newdat[,x]),digits = dec)
      model <- summary(lm(newdat[,yvar] ~ newdat[,x]))
      est <- paste0(
        round(model$coefficients[2,1],digits = dec),
        " (",round(model$coefficients[2,2],digits = dec),")"
      )
      rsq <- round(model$r.squared,digits = rsqdec)
      pvalue <- round(model$coefficients[2,4],digits = pdec)
      pvalue <-
        ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
      line <-
        c(x,"N  Mean (Std.Dev.)",paste0(N,"  ",mean," (",stdev, ")"),est,rsq,pvalue)
      names(line) <-
        c(
          "Variable","Response","Descriptive Stats", "Coeff. Est.","R Squared","p-value"
        )
      return(line)
    }

    if (x %in% catlist) {
      len <- length(levels(factor(newdat[,x])))
      tab <- t(table(newdat[,x]))
      proptab <- round(t(100 * prop.table(tab)),digits = dec)
      mat <- data.frame(matrix(ncol = 6,nrow = len))
      names(mat) <-
        c(
          "Variable","Response","Descriptive Stats", "Coeff. Est.","R Squared","p-value"
        )
      mat[1,1] <- x
      mat[,2] <- levels(factor(newdat[,x]))
      mat[,3] <- paste0(tab,"  (",proptab,"%)")
      model <- summary(lm(newdat[,yvar] ~ newdat[,x]))
      coeff <-
        c("Reference",paste0(
          round(model$coefficients[-1,1],digits = dec),
          " (",round(model$coefficients[-1,2],digits =
                       dec),")"
        ))
      mat[,4] <- coeff
      mat[1,5] <- round(model$r.squared,digits = rsqdec)
      pvalue1 <- round(
        pf(
          q = model$fstatistic[1],
          df1 = model$fstatistic[2],
          df2 = model$fstatistic[3], lower.tail = FALSE
        ), digits = pdec
      )
      pvalueRest <-
        round(model$coefficients[2:length(model$coefficients[,1]),4],digits = pdec)
      pvalue <- c(pvalue1, pvalueRest)
      pvalue <-
        ifelse(pvalue == 1,"> 0.999",ifelse(pvalue == 0,"< 0.001",pvalue))
      mat[,6] <- pvalue
      mat[is.na(mat)] <- " "
      line <- mat
      return(line)
    }
    return(line)
  }
  results <- do.call(rbind, lapply(predlist,getall))
  #Sort the vars
  if (sortVars != FALSE) {
    list <- data.frame(sortnums = 1:(length(sortVars)),Variable = sortVars)
    total <- merge(list,results,by = "Variable")
    finaldat <- total[with(total, order(sortnums)),]
    results <- finaldat[,-(2)]
  }




  if (varLabelTable != FALSE) {
    results[,1] <- as.character(results[,1])
    for (i in 1:length(results[,1])) {
      for (j in 1:length(varLabelTable[,1])) {
        if (as.character(results[i,1]) == as.character(varLabelTable[j,1])) {
          results[i,1] <- as.character(varLabelTable[j,2])
        }
      }
    }
  }
return(results)
}
