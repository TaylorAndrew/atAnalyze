#' multGroup
#'
#' multGroup provides summary statistics for continous and categorical variables either overall or stratified by a grouping variable. It also provides parametric or nonparametric tests for each variable.
#'
#' @param data Input data.frame
#' @param PcontinVars Vector of continous variables names from data that should be tested parametrically
#' @param PcatVars Vector of categorical variables names from data that should be tested parametrically
#' @param NPcontinVars Vector of continous variables names from data that should be tested non-parametrically
#' @param NPcatVars Vector of categorical variables names from data that should be tested non-parametrically
#' @param SortVars Vector of variable names indicating in what order the variables should be presented in the output table
#' @param varLabelTable A varLabelTable object. Converts variable names to variable labels in output table.
#' @param grouping Grouping variable of any number of levels that will be used to stratify each variable
#' @param pdec Number of decimals for p-values
#' @param dec Number of decimals for descriptive statistics
#' @param mu For a one sample t-test (no grouping), the reference value to compare population to. Default of 0.
#' @param ChiProbabilities For a one sample chi-square (no grouping), the underlying probability distribution to compare the population to. If NULL, will compare to an even distribution across all levels.
#' @param labels Labels for the grouping variable
#' @param verbose If TRUE, additional information is printed to the console
#' @param percent Indicates if the output table should present 'row', 'column', or 'overall' percents for categorical variables.
#' @param NPdescriptives Vector of continuous variable names that should present median and range in the output table.
#' @param padjust If TRUE and grouping has 3 or more levels, only pairwise comparisons that are significant at 0.05 after using a multiple comparisons adjustment are presented
#' @param provideP If TRUE, in the case of no grouping variable, the p-values for the 1-sample tests will be included in the output, otherwise only the descriptive statistics are returned
#' @param include If 'none', no additional output is included. If 'range', range is included. If '95ci' then 95% Confidence intervals are included.
#'
#' @return A data.frame including all univariate results.
#' @export
#'
#' @examples
#'
#' #df <- data.frame(time = runif(100, 1, 100),
#' #                  event = sample(c(0,1), 100, replace = T, prob = c(.1, .9)),
#' #                  x1 = rnorm(100),
#' #                  x2 = sample(c(1,2), 100, replace = T))
#' #multGroup(data = df,
#' #          NPcontinVars = "time",
#' #          grouping = "event",
#' #          PcontinVars = "x1",
#' #          NPcatVars = "x2",
#' #          keepOutput=T,
#' #          printOutput = F)
multGroup <- function(data,
                      PcontinVars = NULL,
                      PcatVars = NULL,
                      NPcontinVars = NULL,
                      NPcatVars = NULL,
                      SortVars = NULL,
                      varLabelTable = NULL,
                      grouping = NULL,
                      pdec = 3,
                      dec = 2,
                      mu = 0,
                      ChiProbabilities = NULL,
                      labels = NULL,
                      percent = "column",
                      NPdescriptives = NULL,
                      NPvariance = 'range',
                      padjust = F,
                      provideP = T,
                      include = 'none',
                      verbose = T) {
  deleteExtraneous <- F
  data <- as.data.frame(data)
  if(!is.null(grouping)) {
    data <- subset(data, !is.na(data[, grouping]))
  }
  if (is.null(grouping)) {
    keepNames <- names(data)
    data1 <- data.frame(data, fakeVar = "All Patients")
    data2 <- data.frame(data, fakeVar = "zzzzzzz")
    data <- rbind(data1, data2)
    names(data) <- c(keepNames, "fakeVar")
    data$fakeVar <- factor(data$fakeVar)
    grouping <- "fakeVar"
    deleteExtraneous <- T
  }
  if (!is.null(grouping) & !is.null(labels)) {
    if (length(labels) != length(levels(factor(data[, grouping])))) {
      return(paste0(
        "Number of labels does not match number of levels of ",
        "grouping variable"
      ))
    }
  }
  ciProp <- function(x, conf.int = .95) {
  x <- x[!is.na(x)]
  n <- length(x)
  pbar =  unlist(table(x))/n
  SE = sqrt(pbar*(1-pbar)/n)
  E = qnorm(conf.int +(1-conf.int)/2)*SE
  out = do.call(rbind, lapply(1:length(pbar), function(i) pbar[i] + c(-E[i], E[i])))
  row.names(out) <- names(pbar)
  colnames(out) <- c("ci.lower", "ci.upper")
  return(out)
}
  #For parametric continuous variables
  if (!is.null(PcontinVars) | !is.null(NPcontinVars)) {
    if (verbose == T) {
      nnv <- function(x) {
        is.na1 <- is.na(x)
        x2 <- suppressWarnings(as.numeric(as.character(x)))
        is.na2 <- is.na(x2)
        x[is.na2 & !is.na1]
      }
    }
  }
  if (!is.null(PcontinVars)) {
    #Get a list of the groups being compared
    flist <- sort(unique(data[, grouping]))
    #Function to get descriptive statistics per group
    descriptive <- function(x) {
      if (verbose == T)
        print(paste0(x, " is now being analyzed"))
      if (verbose == T) {
        toNA <- nnv(data[, x])
        if (length(toNA) > 0) {
          print(paste0(
            x,
            ": The following entries will be converted to NA prior to analysis"
          ))
          print(toNA)
        }
      }
      #data for a specified variable
      xlist <- suppressWarnings(as.numeric(as.character(data[, x])))
      #Number of groups being compared
      numlevs <- length(flist)
      #making a list of 1 to number of groups
      numdos <- 1:numlevs
      #finally the meat and potatoes of getting the desc. stats.
      overallxs <- xlist[!is.na(xlist)]
      overallN <- length(overallxs)
      overallmean <-
        sprintf(paste0("%.", dec, "f"), round(mean(overallxs), digits = dec))
      overallsd <-
        sprintf(paste0("%.", dec, "f"), round(sd(overallxs), digits = dec))
      overallmed <- round(median(overallxs), digits = dec)
        overallrange <-
          paste0(round(range(overallxs)[1], digits = dec),
                 ", ",
                 round(range(overallxs)[2], digits = dec))
        overallIQR <- round(IQR(overallxs), digits = dec)
        overallQuartiles <-  paste0(round(quantile(overallxs)[2], digits = dec),
                 ", ",
                 round(quantile(overallxs)[4], digits = dec))
      if (include == '95ci') {
        overallSE <- sd(overallxs) / sqrt(overallN)
        E <- qt(.975, overallN - 1) * overallSE
        overall95ci <- paste0(" [",
                              sprintf(paste0("%.", dec, "f"), round(mean(overallxs) - E, digits = dec)),
                              ", ",
                              sprintf(paste0("%.", dec, "f"), round(mean(overallxs) + E, digits = dec)),
                              "]")
      }

      if (include == 'none') {
        overallstats <-
          paste0(
            overallN,
            "  ",
            if(!is.null(NPdescriptives) &
                (x %in% NPdescriptives)) {
              overallmed
            } else {overallmean} ,
            " ",
            "(",
            if(!is.null(NPdescriptives) &
                (x %in% NPdescriptives)) {
              if(NPvariance=='range'){
                overallrange } else if(NPvariance=='IQR') {
                overallIQR } else if(NPvariance=='innerQuartiles') {
                overallQuartiles
                }
                } else {
              overallsd
            },
             ")"
            )
      }
      if (include == 'range') {
        overallstats <- paste0(overallN,
                               "  ",
                               overallmean,
                               " ",
                               "(",
                               overallsd,
                               ") ",
                               overallrange)
      }
      if (include == '95ci') {
        overallstats <- paste0(overallN,
                               "  ",
                               overallmean,
                               " ",
                               "(",
                               overallsd,
                               ") ",
                               overall95ci)
      }
      getPstats <- function(level) {
        #what group level are we looking at?
        grouplevel <- level
        #get only the data for individuals in that group
        xs <- as.numeric(xlist[data[, grouping] == flist[level]])
        #get rid of all NAs
        xclean <- xs[!is.na(xs)]
        #get the mean, rounded
        mean <-
          sprintf(paste0("%.", dec, "f"), round(mean(xclean), digits = dec))
        #get the std. dev, rounded
        sd1 <-
          sprintf(paste0("%.", dec, "f"), round(sd(xclean), digits = dec))
        #get the n
        n <- length(xclean)

        if (include == '95ci') {
          oneSE <- sd(xclean) / sqrt(n)
          E <- qt(.975, n - 1) * oneSE
          one95ci <- paste0(" [",
                            sprintf(paste0("%.", dec, "f"), round(mean(xclean) - E, digits = dec)),
                            ", ",
                            sprintf(paste0("%.", dec, "f"), round(mean(xclean) + E, digits = dec)),
                            "]")
        }


        #bind it all together
        if (include == 'none') {
          stats <- paste0(n, "  ", mean, " ", "(", sd1, ")")
        }
        if (include == 'range') {
          stats <- paste0(n,
                          "  ",
                          mean,
                          " ",
                          "(",
                          sd1,
                          ") ",
                          paste0(
                            round(range(xclean)[1], digits = dec),
                            ", ",
                            round(range(xclean)[2], digits =
                                    dec)
                          ))
        }
        if (include == '95ci') {
          stats <- paste0(n,
                          "  ",
                          mean,
                          " ",
                          "(",
                          sd1,
                          ") ",
                          one95ci)
        }
        #return the result 'up' a level to be used in the next function
        return(stats)
      }
      getNPstats <- function(level) {
        #what group level are we looking at?
        grouplevel <- level
        #get only the data for individuals in that group
        xs <- as.numeric(xlist[data[, grouping] == flist[level]])
        #get rid of all NAs
        xclean <- xs[!is.na(xs)]
        med <-
          sprintf(paste0("%.", dec, "f"), round(median(xclean), digits = dec))
        range <-
          paste0(round(range(xclean)[1], digits = dec),
                 ", ",
                 round(range(xclean)[2], digits =
                         dec))
        IQR <- round(IQR(xclean), digits = dec)
        Quartiles <-  paste0(round(quantile(xclean)[2], digits = dec),
                 ", ",
                 round(quantile(xclean)[4], digits = dec))
        #get the n
        n <- length(xclean)
        #bind it all together
        if(NPvariance=='range') stats <- paste0(n, "  ", med, " ", "(", range, ")")
        if(NPvariance=='IQR') stats <- paste0(n, "  ", med, " ", "(", IQR, ")")
        if(NPvariance=='innerQuartiles') stats <- paste0(n, "  ", med, " ", "(", Quartiles, ")")
        #return the result 'up' a level to be used in the next function
        return(stats)
      }
      #we column bind the descriptive results so that we have them for each
      #level of group
      firstpart <-
        do.call(cbind, lapply(numdos, if (x %in% NPdescriptives) {
          getNPstats
        } else {
          getPstats
        }))
      #now we decide what test to run, based off the number of groups and we
      #do that test
      #also decides based off parametric or non-parametric
      if (numlevs == 2)
        Test <-
        "t.test"
      else if (numlevs == 1)
        Test <- "Single Sample"
      else
        Test <- "ANOVA"
      if (Test == "Single Sample") {
        pvalue <-
          sprintf(paste0("%.", pdec, "f"),
                  round(
                    t.test(data[, x],
                           alternative = "two.sided",
                           mu = mu)$p.value,
                    digits = pdec
                  ))
        if (pvalue == 0)
          pvalue <- "<.0001"
        testname <- "1 sample ttest"
        sigpcol = "None"
      }
      if (Test == "t.test") {
        x1 <-
          suppressWarnings(as.numeric(as.character(data[, x][data[, grouping] == flist[1]])))
        x1clean <- x1[!is.na(x1)]
        x2 <-
          suppressWarnings(as.numeric(as.character(data[, x][data[, grouping] == flist[2]])))
        x2clean <- x2[!is.na(x2)]
        tc <- tryCatch(t.test(x1clean, x2clean),
                        warning = function(w)
                  w,
                error = function (e)
                  e)
        if (is(tc, "error")) {
        testname <- "None"
        sigpcol = "None"
        pvalue <- 'None'
        sigpcol <- "None"
        } else {
        t <- t.test(x1clean, x2clean)
        testname <- "2-sample t-test"
        sigpcol = "None"
        pvalue <-
          sprintf(paste0("%.", pdec, "f"),
                  round(t$p.value, digits = pdec))
        pvalue <-
          gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
        pvalue <- gsub("^1\\.", "> 0.", pvalue)
        pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
        sigpcol <- "None"
        }

      }
      if (Test == "ANOVA") {
        anovasum <-
          summary(aov(data[, x] ~ factor(data[, grouping])))[[1]][["Pr(>F)"]]
        testname <- "ANOVA"
        modelfortuk <- aov(data[, x] ~ factor(data[, grouping]))
        #extra bit for ANOVA is doing a TukeyHSD multiple comparisons
        if (anovasum[1] <= .05) {
          pairps <- TukeyHSD(modelfortuk)[[1]][, 4]
          sigpsTF <- pairps <= .05
          print(pairps)
          print(sigpsTF)
          sigps <- which(sigpsTF %in% TRUE)
          print(sigps)
          if (length(sigps) == 0) {
            sigpcol = "None"
          } else {
            sigpcol = paste(unlist(sigps), collapse = " ")
          }
        } else {
          sigpcol = "None"
        }
        pvalue <-
          sprintf(paste0("%.", pdec, "f"),
                  round(anovasum[1], digits = pdec))
        pvalue <-
          gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
        pvalue <- gsub("^1\\.", "> 0.", pvalue)
        pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
      }
      if (x %in% NPdescriptives == 1) {
        resp <-
          paste0("N  Median (", NPvariance, ")")
      } else if (include == 'none') {
        resp <- "N  Mean (Std.Dev.)"
      } else if (include == 'range') {
        resp <- "N  Mean (Std.Dev.) Range"
      } else {
        resp <- "N  Mean (Std.Dev.) (95% CI)"
      }
      #pull together the descriptives, and the test results for a single DV
      Poneline <-
        data.frame(resp,
                   overallstats,
                   firstpart,
                   as.character(pvalue),
                   testname,
                   sigpcol)
    }
    #then pull toegether the results for all DVs, and label the columns
    #and rows and print it
    Pmorepre <- do.call(rbind, lapply(PcontinVars, descriptive))
    Pmore <- data.frame(PcontinVars, PcontinVars, Pmorepre)
    colnames(Pmore) <- c(
      "ForSortVars",
      "Variable",
      "Response",
      "Overall",
      if (!is.null(labels)) {
        labels
      } else {
        flist
      },
      "p value",
      "Test",
      "Pairwise"
    )
  }
  #For nonparametric continuous variables
  if (!is.null(NPcontinVars)) {
    #Get a list of the groups being compared
    flist <- sort(unique(data[, grouping]))
    #Function to get descriptive statistics per group
    descriptive <- function(x) {
      if (verbose == T) print(paste0(x, " is now being analyzed"))
      if (verbose == T) {
        toNA <- nnv(data[, x])
        if (length(toNA) > 0) {
          print(paste0(
            x,
            ": The following entries will be converted to NA prior to analysis"
          ))
          print(toNA)
        }
      }
      #data for a specified variable
      xlist <- suppressWarnings(as.numeric(as.character(data[, x])))
      #Number of groups being compared
      numlevs <- length(flist)
      #making a list of 1 to number of groups
      numdos <- 1:numlevs
      #finally the meat and potatoes of getting the desc. stats.
      overallxs <- xlist[!is.na(xlist)]
      overallN <- length(overallxs)
      overallmean <-
        sprintf(paste0("%.", dec, "f"), round(mean(overallxs), digits = dec))
      overallsd <-
        sprintf(paste0("%.", dec, "f"), round(sd(overallxs), digits = dec))
      overallmed <- round(median(overallxs), digits = dec)
      overallrange <-
        paste0(round(range(overallxs)[1], digits = dec),
               ", ",
               round(range(overallxs)[2], digits = dec))

      if (include == '95ci') {
        overallSE <- sd(overallxs) / sqrt(overallN)
        E <- qt(.975, overallN - 1) * overallSE
        overall95ci <- paste0(" [",
                              sprintf(paste0("%.", dec, "f"), round(mean(overallxs) - E, digits = dec)),
                              ", ",
                              sprintf(paste0("%.", dec, "f"), round(mean(overallxs) + E, digits = dec)),
                              "]")
      }



      if (include == 'none') {
        overallstats <- paste0(
          overallN,
          "  ",
          ifelse(
            !is.null(NPdescriptives) & (x %in% NPdescriptives),
            overallmed,
            overallmean
          ),
          " ",
          "(",
          if(!is.null(NPdescriptives) & (x %in% NPdescriptives)) {
            overallrange
          } else {overallsd}
            ,
            ")"
          )

      }
      if (include == 'range') {
        overallstats <- paste0(overallN,
                               "  ",
                               overallmean,
                               " ",
                               "(",
                               overallsd,
                               ") ",
                               overallrange)
      }
      if (include == '95ci') {
        overallstats <- paste0(overallN,
                               "  ",
                               overallmean,
                               " ",
                               "(",
                               overallsd,
                               ") ",
                               overall95ci)
      }
      getPstats <- function(level) {
        #what group level are we looking at?
        grouplevel <- level
        #get only the data for individuals in that group
        xs <-
          suppressWarnings(as.numeric(as.character(xlist[data[, grouping] == flist[level]])))
        #get rid of all NAs
        xclean <- xs[!is.na(xs)]
        #get the mean, rounded
        mean <-
          sprintf(paste0("%.", dec, "f"), round(mean(xclean), digits = dec))
        #get the std. dev, rounded
        sd1 <-
          sprintf(paste0("%.", dec, "f"), round(sd(xclean), digits = dec))
        #get the n
        n <- length(xclean)


        if (include == '95ci') {
          oneSE <- sd(xclean) / sqrt(n)
          E <- qt(.975, n - 1) * oneSE
          one95ci <- paste0(" [",
                            sprintf(paste0("%.", dec, "f"), round(mean(xclean) - E, digits = dec)),
                            ", ",
                            sprintf(paste0("%.", dec, "f"), round(mean(xclean) + E, digits = dec)),
                            "]")
        }


        #bind it all together
        if (include == 'none') {
          stats <- paste0(n, "  ", mean, " ", "(", sd1, ")")
        }
        if (include == 'range') {
          stats <- paste0(n,
                          "  ",
                          mean,
                          " ",
                          "( ",
                          sd1,
                          ") ",
                          paste0(
                            round(range(xclean)[1], digits = dec),
                            ", ",
                            round(range(xclean)[2], digits = dec)
                          ))
        }
        if (include == '95ci') {
          stats <- paste0(n,
                          "  ",
                          mean,
                          " ",
                          "(",
                          sd1,
                          ") ",
                          one95ci)
        }
        #return the result 'up' a level to be used in the next function
        return(stats)
      }
      getNPstats <- function(level) {
        #what group level are we looking at?
        grouplevel <- level
        #get only the data for individuals in that group
        xs <-
          suppressWarnings(as.numeric(as.character(xlist[data[, grouping] == flist[level]])))
        #get rid of all NAs
        xclean <- xs[!is.na(xs)]
        med <-
          sprintf(paste0("%.", dec, "f"), round(median(xclean), digits = dec))
        IQRstat <- paste0(round(range(xclean)[1], digits = dec),
                          ", ",
                          round(range(xclean)[2], digits = dec))
        #get the n
        n <- length(xclean)
        #bind it all together
        stats <- paste0(n, "  ", med, " ", "(", IQRstat, ")")
        #return the result 'up' a level to be used in the next function
        return(stats)
      }
      #we column bind the descriptive results so that we have them for each
      #level of group
      firstpart <- do.call(cbind,
                           lapply(numdos,
                                  if (x %in% NPdescriptives) {
                                    getNPstats
                                  } else {
                                    getPstats
                                  }))
      #now we decide what test to run, based off the number of groups and
      #we do that test also decides based off parametric or non-parametric
      if (numlevs == 2)
        Test <-
        "t.test"
      else if (numlevs == 1)
        Test <- "Single Sample"
      else
        Test <- "ANOVA"
      if (Test == "Single Sample") {
        pvalue <- sprintf(paste0("%.",
                                 pdec,
                                 "f"),
                          round(wilcox.exact(data[, x])$p.value, digits =
                                  pdec))
        pvalue <-
          gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
        pvalue <- gsub("^1\\.", "> 0.", pvalue)
        pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
        testname <- "Rank Sum"
        sigpcol = "None"
      }
      if (Test == "t.test") {
        x1 <-
          suppressWarnings(as.numeric(as.character(data[, x][data[, grouping] == flist[1]])))
        x1clean <- x1[!is.na(x1)]
        x2 <-
          suppressWarnings(as.numeric(as.character(data[, x][data[, grouping] == flist[2]])))
        x2clean <- x2[!is.na(x2)]
        tC <- tryCatch(
          wilcox.test(x1clean, x2clean, exact = T),
          warning = function(w)
            w,
          error = function(e)
            e
        )
        if (is(tC, "warning")) {
          ddd <-
            data.frame(grp = factor(c(
              rep(0, length(x1clean)), rep(1, length(x2clean))
            )),
            val = c(x1clean, x2clean))
          t <- wilcox_test(val ~ grp, ddd)
          pvalue <- sprintf(paste0("%.", pdec, "f"),
                            round(pvalue(t, distribution = exact),
                                  digits = pdec))
          if (verbose == T) {
            print(paste0(
              x,
              ": Ties found, exact p-value is computed using mid-ranks method"
            ))
          }
        } else {
          t <- wilcox.test(x1clean, x2clean, exact = T)
          pvalue <-
            sprintf(paste0("%.", pdec, "f"),
                    round(t$p.value, digits = pdec))
        }
        testname <- "Mann-Whitney U"
        sigpcol = "None"
        pvalue <-
          gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
        pvalue <- gsub("^1\\.", "> 0.", pvalue)
        pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
        sigpcol <- "None"
      }
      if (Test == "ANOVA") {
        anovasum <-
          kruskal.test(data[, x] ~ data[, grouping])$p.value
        testname <- "Kruskal-Wallis"
        #Again we have to do multiple comparisons since more than 2 groups
        sigpsTF <-
          kruskalmc(data[, x] ~ data[, grouping])$dif.com[, 3]
        sigps <- which(sigpsTF %in% TRUE)
        if (length(sigps) == 0) {
          sigpcol = "None"
        } else {
          sigpcol = paste(unlist(sigps), collapse = " ")
        }
        pvalue <-
          sprintf(paste0("%.", pdec, "f"),
                  round(anovasum[1], digits = pdec))
        if (pvalue == 0)
          pvalue <- "<.0001"
        if (pvalue == 1)
          pvalue <- ">.999"
      }
      if (x %in% NPdescriptives == 1) {
        resp <-
          "N  Median (Range)"
      } else if (include == 'none') {
        resp <- "N  Mean (Std.Dev.)"
      } else if (include == 'range') {
        resp <- "N  Mean (Std.Dev.) Range"
      } else {
        resp <- "N  Mean (Std.Dev.) (95% CI)"
      }
      #pull together the descriptives, and the test results for a single DV
      NPoneline <-
        data.frame(resp,
                   overallstats,
                   firstpart,
                   as.character(pvalue),
                   testname,
                   sigpcol)
    }
    #then pull toegether the results for all DVs, and label the columns and
    #rows and print it
    NPmorepre <- do.call(rbind, lapply(NPcontinVars, descriptive))
    NPmore <- data.frame(NPcontinVars, NPcontinVars, NPmorepre)
    colnames(NPmore) <-
      c(
        "ForSortVars",
        "Variable",
        "Response",
        "Overall",
        if (!is.null(labels)) {
          labels
        } else {
          flist
        },
        "p value",
        "Test",
        "Pairwise"
      )
  }
  if (!is.null(PcontinVars) & !is.null(NPcontinVars)) {
    more <- rbind(Pmore, NPmore)
  }
  if (!is.null(PcontinVars) & is.null(NPcontinVars)) {
    more <- Pmore
  }
  if (is.null(PcontinVars) & !is.null(NPcontinVars)) {
    more <- NPmore
  }
  if (!is.null(PcatVars)) {
    posthoc.CHI <-
      function (chi,
                popsInRows = FALSE,
                control = c("fdr",
                            "BH",
                            "BY",
                            "bonferroni",
                            "holm",
                            "hochberg",
                            "hommel"),
                digits = 4)
      {
        control <- match.arg(control)
        tbl <- chi$observed
        if (!popsInRows)
          tbl <- t(tbl)
        popsNames <- rownames(tbl)
        prs <- combn(1:nrow(tbl), 2)
        tests <- ncol(prs)
        pvals <- numeric(tests)
        lbls <- character(tests)
        for (i in 1:tests) {
          pvals[i] <- suppressWarnings(chisq.test(tbl[prs[, i],])$p.value)
          lbls[i] <- paste(popsNames[prs[, i]], collapse = " vs. ")
        }
        adj.pvals <- p.adjust(pvals, method = control)
        #       cat("Adjusted p-values used the", control, "method.\n\n")
        data.frame(
          comparison = lbls,
          raw.p = round(pvals, digits),
          adj.p = round(adj.pvals, digits)
        )
      }
    #Location on catVars for which variable name to pull
    indexedListLocation = 0
    get.cat.stats <-
      function(catVars,
               group = data[, grouping],
               dec = dec,
               ...) {
        get.chi.stuff <- function(var) {
          #Don't use <<, unless you have to...
          indexedListLocation <<- indexedListLocation + 1
          if (verbose == T)
            print(paste0(catVars[indexedListLocation], " is now being analyzed"))
          #Unique levels of group
          flist <- sort(unique(group))
          #is our table a contingency or not?
          if (length(flist) == 1) {
            long <- table(var)
          }
          else {
            long <- table(var, group)
          }
          overalllong <- table(var)
          overallprop <-
            t(t(sprintf(
              paste0("%.", dec, "f"),
              round(100 * prop.table(overalllong), digits = dec)
            )))
           if(include == "95ci") {
            ciPropTable <- round(100*ciProp(var), 1)
            ciPropTable <- paste0("[", ciPropTable[,1], "%, ", ciPropTable[,2], "%]")
          }
          #get proportions for the table
          longprop <- sprintf(paste0("%.", dec, "f"),
                              round(
                                100 * prop.table(long, if (percent == "column" &
                                                           length(flist) != 1) {
                                  2
                                }
                                else if (percent == "column") {
                                  NULL
                                }
                                else if (percent == "row") {
                                  1
                                }
                                else if (percent == "overall") {
                                  NULL
                                }),
                                digits = dec
                              ))

          #pick a test to use
          if (provideP == TRUE) {
            if (length(flist) == 1) {
              if (mean(ChiProbabilities == FALSE) == 1 |
                  mean(c(
                    j <- mean(ChiProbabilities != FALSE),
                    k <-
                    (length(levels(
                      factor(var)
                    )) != length(ChiProbabilities))
                  )) == 1) {
                if (ChiProbabilities != FALSE) {
                  print(
                    paste0(
                      "Number of Probabilites listed in ChiProbabilities",
                      "\ndoes not equal the number of levels of ",
                      catVars[indexlocation]
                    )
                  )
                }
                test <-
                  chisq.test(table(as.character(var)))
              } else {
                test <- chisq.test(table(as.character(var)), p = ChiProbabilities)
              }
              pvalue <-
                sprintf(paste0("%.", pdec, "f"),
                        round(test$p.value, digits = pdec))
              pvalue <-
                gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
              pvalue <- gsub("^1\\.", "> 0.", pvalue)
              pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
              sigpcol <- "N/A"
              method <- "Chi-Square Test"
            }
          }
          if (provideP == FALSE)
          pvalue = "NA"
          sigpcol = "NA"
          method = "NA"
          if (length(flist) > 1) {
            if (length(levels(factor(var))) > 1) {
              if (length(levels(factor(group))) > 1) {
                tryC <- tryCatch(
                  test <- chisq.test(var, group),
                  warning = function(w)
                    w,
                  error = function (e)
                    e
                )
                if (grepl("'x' and 'y' must have at least", tryC[1]) |
                    grepl("'x' and 'y' must have at least", tryC[2])) {
                  pvalue <- "N/A"
                  sigps <- NULL
                  sigpcol <- "N/A"
                } else {
                  test <- chisq.test(as.character(as.factor(var)), group)
                  pvalue <-
                    sprintf(paste0("%.", pdec, "f"),
                            round(test$p.value, digits = pdec))
                  ps <-
                    posthoc.CHI(test)$raw.p <= (ifelse(padjust == TRUE, .05 /
                                                         sum(1:(
                                                           length(flist)
                                                         ) - 1), .05))
                  sigps <- which(ps %in% TRUE)
                }
              } else {
                pvalue <- "NA"
                sigps <- NULL
              }
              if (length(sigps) == 0) {
                sigpcol = "None"
              } else {
                sigpcol = paste(unlist(sigps), collapse = " ")
              }
            }

            if (length(levels(factor(var))) == 1) {
              pvalue <- "N/A"
              sigpcol <- "N/A"
            }
            pvalue <-
              gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
            pvalue <- gsub("^1\\.", "> 0.", pvalue)
            pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
            method <- "Chi-Square Test"
          }
          #get test stats
          teststats <- c(pvalue, method)
          #make a matrix to put all of the results into
          res <-
            data.frame(matrix(NA, nrow(long), 7 + if (length(flist) == 1) {
              1
            } else {
              ncol(long)
            }))
          #make the column names for the matrix
          names(res) <-
            c(
              "ForSortVars",
              "Variable",
              "Response",
              "Overall",
              if (!is.null(labels)) {
                labels
              } else {
                flist
              },
              "p value",
              "Test",
              "Pairwise"
            )
          #make the first column, first row, the name of the variable,
          #from the list you provided
          res[, 1] <-
            rep(catVars[indexedListLocation], length(res[, 1]))
          res[1, 2] <- catVars[indexedListLocation]
          #here the row.names are the possible responses to the categorical
          #variable
          res[, 3] <- row.names(long)
          #now we are putting in our table of raw counts (and percentages)
          #overall
          if (percent == "row") {
            res[, 4] <-
            paste(overalllong, paste("(", "100", "%", ")", sep = ''))
          }
          if (percent != "row") {
            if(include != "95ci") {
              res[, 4] <-
            paste(overalllong, paste("(", overallprop, "%", ")", sep = ''))
            }
            if(include == "95ci") {
              res[, 4] <-
             paste(overalllong, paste("(", overallprop, "%", ") ", ciPropTable, sep = ''))
            }
          }

          #number of columns being used is contengent on the size of the
          #grouping variable
          res[, 5:(4 + length(flist))] <-
            paste(long, paste("(", longprop, "%", ")", sep = ''))
          #add in test stats to the appropriate columns
          res[1, (5 + length(flist)):(6 + length(flist))] <-
            teststats
          res[1, (7 + length(flist))] <- sigpcol
          #make all NA cells into blank cells
          res[is.na(res)] <- " "
          #and return it for later use in the higher up function
          return(res)
        }
        #do get.chi.stuff for all categorical variables in the list
        if (length(PcatVars) > 1) {
          chitable <- do.call(rbind, lapply(data[, PcatVars], get.chi.stuff))
        }
        if (length(PcatVars) == 1) {
          chitable <- get.chi.stuff(data[, PcatVars])
        }
        #and return that up a level
        return(chitable)
      }
    Pcattable <-
      get.cat.stats(
        catVars = PcatVars,
        group = data[, grouping],
        data = data,
        dec = dec,
        pdec =
          pdec,
        testtype = testtype
      )
  }
  if (!is.null(NPcatVars)) {
    need <- c("0", "1")
    posthoc.CHI <-
      function (chi,
                popsInRows = FALSE,
                control = c("fdr",
                            "BH",
                            "BY",
                            "bonferroni",
                            "holm",
                            "hochberg",
                            "hommel"),
                digits = 4)
      {
        control <- match.arg(control)
        tbl <- chi$observed
        if (!popsInRows)
          tbl <- t(tbl)
        popsNames <- rownames(tbl)
        prs <- combn(1:nrow(tbl), 2)
        tests <- ncol(prs)
        pvals <- numeric(tests)
        lbls <- character(tests)
        for (i in 1:tests) {
          pvals[i] <- suppressWarnings(chisq.test(tbl[prs[, i],])$p.value)
          lbls[i] <- paste(popsNames[prs[, i]], collapse = " vs. ")
        }
        adj.pvals <- p.adjust(pvals, method = control)
        #       cat("Adjusted p-values used the", control, "method.\n\n")
        data.frame(
          comparison = lbls,
          raw.p = round(pvals, digits),
          adj.p = round(adj.pvals, digits)
        )
      }
    #Location on catVars for which variable name to pull
    indexedListLocation = 0
    get.cat.stats <-
      function(catVars,
               group = data[, grouping],
               dec = dec,
               ...) {
        get.chi.stuff <- function(var) {
          #Don't use <<, unless you have to...
          indexedListLocation <<- indexedListLocation + 1
          if (verbose == T)
            print(paste0(catVars[indexedListLocation], " is now being analyzed"))
          #Unique levels of group
          flist <- sort(unique(group))
          #is our table a contingency or not?
          if (length(flist) == 1) {
            long <- table(var)
          } else {
            long <- table(var, group)
          }
          #get proportions for the table
          overalllong <- table(var)
          overallprop <-
            t(t(sprintf(
              paste0("%.", dec, "f"),
              round(100 * prop.table(overalllong), digits = dec)
            )))
          longprop <-
            sprintf(paste0("%.", dec, "f"),
                    round(100 * prop.table(long, if (percent == "column") {
                      2
                    }
                    else if (percent == "row") {
                      1
                    }
                    else if (percent == "overall") {
                      NULL
                    }), digits = dec))
          #pick a test to use
          if (length(flist) == 1) {
            if (mean(ChiProbabilities == FALSE) == 1 |
                mean(c(
                  j <-
                  mean(ChiProbabilities != FALSE),
                  k <-
                  (length(levels(
                    factor(var)
                  )) != length(ChiProbabilities))
                )) == 1) {
              if (ChiProbabilities != FALSE) {
                print(
                  paste0(
                    "Number of Probabilites listed in ChiProbabilities",
                    "\ndoes not equal the number of levels of ",
                    catVars[indexlocation]
                  )
                )
              }
              test <-
                chisq.test(table(as.factor(as.character(var))))
            } else {
              test <- chisq.test(table(as.factor(as.character(var))), p = ChiProbabilities)
            }
            if (pvalue != "N/A")
              pvalue <-
                sprintf(paste0("%.", pdec, "f"),
                        round(test$p.value, digits = pdec))
            if (pvalue != "N/A")
              pvalue <-
                gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
            if (pvalue != "N/A")
              pvalue <- gsub("^1\\.", "> 0.", pvalue)
            if (pvalue != "N/A")
              pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
            sigpcol <- "N/A"
            method <- "Fisher's Exact Test"
          }
          if (length(flist) > 1) {
            if (length(levels(factor(var))) > 1) {
              tryC <- tryCatch(
                fisher.test(var, group),
                warning = function(w)
                  w,
                error = function (e)
                  e
              )
              if (grepl("'x' and 'y' must have at least", tryC[1])) {
                pvalue <- "N/A"
                sigpcol <- "N/A"
                sigps <- NULL
              } else {
                test <- fisher.test(as.factor(as.character(var)), group)
                test2 <- suppressWarnings(chisq.test(as.factor(as.character(var)), group))
                pvalue <-
                  sprintf(paste0("%.", pdec, "f"),
                          round(test$p.value, digits = pdec))
                if (pvalue != "N/A")
                  pvalue <-
                  gsub("(?:^1\\.|\\G)\\K0(?=0*$)",
                       "9",
                       pvalue,
                       perl = T)
                if (pvalue != "N/A")
                  pvalue <- gsub("^1\\.", "> 0.", pvalue)
                if (pvalue != "N/A")
                  pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
                ps <-
                  posthoc.CHI(test2)$raw.p <= (ifelse(padjust == TRUE, .05 /
                                                        sum(1:(
                                                          length(flist)
                                                        ) - 1), .05))
                sigps <- which(ps %in% TRUE)
                if (length(sigps) == 0) {
                  sigpcol = "None"
                } else {
                  sigpcol = paste(unlist(sigps), collapse = " ")
                }
              }
            } else if (length(levels(factor(var))) == 1) {
              pvalue <- "N/A"
              sigpcol <- "N/A"
            }
            if (pvalue != "N/A")
              pvalue <-
                gsub("(?:^1\\.|\\G)\\K0(?=0*$)", "9", pvalue, perl = T)
            if (pvalue != "N/A")
              pvalue <- gsub("^1\\.", "> 0.", pvalue)
            if (pvalue != "N/A")
              pvalue <- gsub("^(0\\.0*)0$", "< \\11", pvalue)
            method <- "Fisher's Exact"
          }

          # }
          #get test stats
          teststats <- c(pvalue, method)
          #make a matrix to put all of the results into
          res <-
            data.frame(matrix(NA, nrow(long), 7 + if (length(flist) == 1) {
              1
            } else {
              ncol(long)
            }))
          #make the column names for the matrix
          names(res) <-
            c(
              "ForSortVars",
              "Variable",
              "Response",
              "Overall",
              if (!is.null(labels)) {
                labels
              } else {
                flist
              },
              "p value",
              "Test",
              "Pairwise"
            )
          #make the first column, first row, the name of the variable,
          #from the list you provided
          res[, 1] <-
            rep(catVars[indexedListLocation], length(res[, 1]))
          res[1, 2] <- catVars[indexedListLocation]
          #here the row.names are the possible responses to the categorical
          #variable
          res[, 3] <- row.names(long)
          #now we are putting in our table of raw counts (and percentages)
          #overall
          if (percent == "row")
            res[, 4] <-
            paste(overalllong, paste("(", "100", "%", ")", sep = ''))
          if (percent != "row")
            res[, 4] <-
            paste(overalllong, paste("(", overallprop, "%", ")", sep = ''))
          #number of columns being used is contengent on the size of the
          #grouping variable
          res[, 5:(4 + length(flist))] <-
            paste(long, paste("(", longprop, "%", ")", sep = ''))
          #add in test stats to the appropriate columns
          res[1, (5 + length(flist)):(6 + length(flist))] <-
            teststats
          res[1, (7 + length(flist))] <- sigpcol
          #make all NA cells into blank cells
          res[is.na(res)] <- " "
          #and return it for later use in the higher up function
          return(res)
        }
        #do get.chi.stuff for all categorical variables in the lsit
        if (length(NPcatVars) > 1) {
          chitable <- do.call(rbind, lapply(data[, NPcatVars], get.chi.stuff))
        }
        if (length(NPcatVars) == 1) {
          chitable <- get.chi.stuff(data[, NPcatVars])
        }
        #and return that up a level
        return(chitable)
      }
    NPcattable <-
      get.cat.stats(
        catVars = NPcatVars,
        group = data[, grouping],
        data = data,
        dec = dec,
        pdec =
          pdec,
        testtype = testtype
      )
  }
  if (!is.null(PcatVars) & !is.null(NPcatVars)) {
    cattable <- rbind(Pcattable, NPcattable)
  }
  if (!is.null(PcatVars) & is.null(NPcatVars)) {
    cattable <- Pcattable
  }
  if (is.null(PcatVars) & !is.null(NPcatVars)) {
    cattable <- NPcattable
  }
  catVars <- NULL
  continVars <- NULL
  if (!is.null(PcatVars) &
      !is.null(NPcatVars))
    catVars <- c(PcatVars, NPcatVars)
  if (!is.null(PcatVars) & is.null(NPcatVars))
    catVars <- PcatVars
  if (is.null(PcatVars) & !is.null(NPcatVars))
    catVars <- NPcatVars
  if (!is.null(PcontinVars) &
      !is.null(NPcontinVars))
    continVars <- c(PcontinVars, NPcontinVars)
  if (!is.null(PcontinVars) &
      is.null(NPcontinVars))
    continVars <- PcontinVars
  if (is.null(PcontinVars) &
      !is.null(NPcontinVars))
    continVars <- NPcontinVars
  if (!is.null(catVars) & !is.null(continVars)) {
    fullresults <- rbind(more, cattable)
    newdat <-
      data.frame(x = 1:length(fullresults[, 1]), fullresults)
  }
  else if (!is.null(catVars)) {
    fullresults <- cattable
    newdat <-
      data.frame(x = 1:length(fullresults[, 1]), fullresults)
  }
  else if (!is.null(continVars)) {
    fullresults <- more
    newdat <-
      data.frame(x = 1:length(fullresults[, 1]), fullresults)
  }
  if (!is.null(SortVars)) {
    list <-
      data.frame(sortnums = 1:(length(SortVars)), ForSortVars = SortVars)
    total <- merge(list, newdat, by = "ForSortVars")
    finaldat <- total[with(total, order(sortnums, x)),]
    final <- finaldat[,-(1:3)]
  }
  else if (is.null(SortVars)) {
    final <- newdat[,-(1:2)]
  }
  if (!is.null(varLabelTable)) {
    final[, 1] <- as.character(final[, 1])
    for (i in 1:length(final[, 1])) {
      for (j in 1:length(varLabelTable[, 1])) {
        if (as.character(final[i, 1]) == as.character(varLabelTable[j, 1])) {
          final[i, 1] <- as.character(varLabelTable[j, 2])
        }
      }
    }
  }
  names(final) <- c("Variable",
                    "Response",
                    "Overall",
                    if (!is.null(labels)) {
                      labels
                    } else {
                      as.character(sort(unique(data[, grouping])))
                    },
                    "p value",
                    "Test",
                    "Pairwise")
  if (deleteExtraneous == T) {
    final <- final[,-c(3, 5:8)]
  }
  if (deleteExtraneous == F &
      length(levels(factor(data[, grouping]))) < 3) {
    final <- final[,-c(length(final[1,]))]
  }
  return(final)
}
