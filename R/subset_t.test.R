#Perform a t.test on a viable across a binary grouping variable for each level
#of a subset variable
subset_t.test <- function(data, #data.frame containing all variables
                          var, #variable name for continuous variable of interest
                          grouping, #variable name for the grouping variable to test differences in 'var' between
                          subsetVar, #variable name for the variable to subset the data by to test at each unique level.
                          var.equal = F, #Logical, if TRUE, variances will be treated as equal.
                          paired = F,#Logical, if TRUE, data will be treated as paired.
                          verbose = F,
                          aovPairedVar = NULL,
                          includeTukey = FALSE) { #Logical, if TRUE, print additional information to console
  lenlevsgrp <- length(levels(factor(data[, grouping])))
  levsgrp <- levels(factor(data[, grouping]))
  # if(lenlevsgrp!=2) return(paste0("Grouping variable ", grouping, " does not have two levels. Number of levels = ", lenlevsgrp))
  subsetLevels <- levels(factor(data[, subsetVar]))
  oneSubset <- function(subsetLev) {
    if(verbose==T) print(paste0("Now doing: ", subsetLev))
    d <- subset(data, data[,subsetVar]==subsetLev)
    descript <- function(grplev) {
      x = d[factor(d[,grouping]) == grplev, var]
      meanx <- round(mean(x, na.rm = T),2)
      sdx <- round(sd(x, na.rm = T),2)
      nx <- length(x[!is.na(x)])
      desc <- data.frame(Level1 = paste0(nx, "  ", meanx, " (", sdx, ")"))
      names(desc) <- grplev
      return(desc)
    }
    descriptives <- do.call(cbind, lapply(levsgrp, descript))
    if(lenlevsgrp!=2) {
    if(paired==F) {
      if(length(d[,var][!is.na(d[,var])])>9 &
       length(d[,grouping][!is.na(d[,grouping])])>9) {
       ttestResult <- round(summary(aov(d[, var]~
                                          factor(d[, grouping]))
                                    )[[1]][["Pr(>F)"]][1], 3)
    } else { ttestResult = "NA"}
    }
    if(paired==T) {
    if(length(d[,var][!is.na(d[,var])])>9 &
       length(d[,grouping][!is.na(d[,grouping])])>9) {
      print(aovPairedVar)
      print(factor(d[, aovPairedVar]))
       ttestResult <- round(summary(aov(d[, var]~
                                          factor(d[, grouping]) +
                                          factor(d[, aovPairedVar]))
                                    )[[1]][["Pr(>F)"]][1], 3)
    } else { ttestResult = "NA"}
    }
    }

    if(lenlevsgrp==2) {
    if(length(d[,var][!is.na(d[,var])])>9 &
       length(d[,grouping][!is.na(d[,grouping])])>9) {
       ttestResult <- round(t.test(d[, var]~ factor(d[, grouping]),
                          var.equal = var.equal, paired = paired)$p.value,
                         3)
    } else { ttestResult = "NA"}
    }
    line <- data.frame(Subset = subsetLev, descriptives, pvalue = ttestResult)
    return(line)
  }
  output <- do.call(rbind, lapply(subsetLevels, oneSubset))
  if(includeTukey==T) {
    tukeyResultsOne <- function(subsetLev) {
    if(verbose==T) print(paste0("Now doing paired comparisons for: ", subsetLev))
    d <- subset(data, data[,subsetVar]==subsetLev)
    if(lenlevsgrp!=2) {
    if(paired==F) {
      if(length(d[,var][!is.na(d[,var])])>9 &
       length(d[,grouping][!is.na(d[,grouping])])>9) {
       mod <- aov(d[, var]~ factor(d[, grouping]))
       tukResult <- TukeyHSD(mod)[[1]]
       tukResult <- data.frame(1:length(tukResult[,1]), tukResult[,4])
       comps <- tukResult[tukResult[,2]<= 0.05, 1]
       comps <- paste0(comps, collapse = ', ')
    } else { comps = "NA"}
    }
    if(paired==T) {
    if(length(d[,var][!is.na(d[,var])])>9 &
       length(d[,grouping][!is.na(d[,grouping])])>9) {
      print(aovPairedVar)
      print(factor(d[, aovPairedVar]))
       mod <- aov(d[, var]~
                          factor(d[, grouping]) +
                          factor(d[, aovPairedVar]))
       tukResult <- TukeyHSD(mod)[[1]]
       tukResult <- data.frame(1:length(tukResult[,1]), tukResult[,4])
       comps <- tukResult[tukResult[,2]<= 0.05, 1]
       comps <- paste0(comps, collapse = ', ')
    } else { comps = "NA"}
    }
    }

    if(lenlevsgrp==2) {comps = "NA"}
    return(comps)
    }
    allcomps <- do.call(rbind, lapply(subsetLevels, tukeyResultsOne))
    allcomps <- data.frame(subset = subsetLevels,
                           SigComps = allcomps)
  }
  if(includeTukey==T) {
    output <- list(output = output, comps = allcomps)
  }
  return(output)
}