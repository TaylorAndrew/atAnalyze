#' multOutcomesLogistic
#'
#' multOutcomesLogistic provides results for univariate binary logistic regression models for a single predictor and multiple outcomes.
#'
#' @param data data.frame containing all data for analyses
#' @param predictor Variable name for the predictor variable
#' @param outcomes Vector of variable names for the outcome variables
#' @param predictorAsFactor If TRUE, the predictor variable will be treated as a factor
#' @param rowwise If TRUE, row-wise percents are provided in the case of categorical predictors
#' @param dec Number of decimals to round descriptive
#' @param pdec Number of decimals to round p-vales
#' @param verbose Logical: should extra information be printed to the console
#'
#' @return A data.frame R object is returned
#' @export
#'
#' @examples
#' #NULL
multOutcomesLogistic <-
  function(data,
           predictor,
           outcomes,
           predictorAsFactor = FALSE,
           rowwise = T,
           dec = 2,
           pdec = 3,
           verbose = TRUE) {
    data <- as.data.frame(data)

    makeone <- function(out) {
    if(verbose==T) {
        print(paste0(out, 'is now being analyzed'))
    }
      nlvs <- length(levels(factor(data[, out])))
      if (nlvs != 2) {
        print(paste0(
          out,
          " does not have 2 levels. ",
          "Only first two levels will be used"
        ))
        data[, out][!factor(data[, out]) %in% levels(factor(data[, out]))[c(1:2)]] <- NA
      }
      levelsNum <- levels(factor(data[, out]))
      if(predictorAsFactor == F) {
        x_y1 <- data[data[, out] == levelsNum[1], predictor]
        x_y1 <- x_y1[!is.na(x_y1)]
        x_y2 <- data[data[, out] == levelsNum[2], predictor]
        x_y2 <- x_y2[!is.na(x_y2)]
        desc_y1 <- paste0(round(mean(x_y1), dec),
                          " (",
                          round(sd(x_y1), dec),
                          ")")
        desc_y2 <- paste0(round(mean(x_y2), dec),
                          " (",
                          round(sd(x_y2), dec),
                          ")")
      }
      if(predictorAsFactor == T) {
        counts <- table(data[, predictor], data[, out])
        percts <- round(prop.table(counts, if(rowwise==T) {1} else {2}), dec)
        desc_y1 <- paste0(counts[,1], " (", percts[,1], ")")
        desc_y2 <- paste0(counts[,2], " (", percts[,2], ")")
      }
      if (max(as.numeric(levelsNum) != c(0, 1))==1) {
        print(paste0(out, " is not coded c(0, 1). ",
                     "Will be re-coded to c(0, 1) for analysis purposes, but will be presented with the original levels."))
        print(paste0("Models will model ", out, "==", levels(factor(data[, out]))[2]))
        data[, out] <- as.numeric(factor(data[, out])) -  1
      }
      if (predictorAsFactor == F)
        mod <-
        glm(data[, out] ~ data[, predictor], data = data, family = "binomial")
      if (predictorAsFactor == T)
        mod <-
        glm(data[, out] ~ factor(data[, predictor]), data = data, family = "binomial")
      output <- broom::tidy(mod, exponentiate = T, conf.int = T)
      output[, 1] <-
        gsub("data\\[, predictor\\]", paste0(predictor, " "), output[, 1])
      output[, 1] <-
        gsub("factor\\(", "", output[, 1])
      output[, 1] <-
        gsub("\\ \\)", ":", output[, 1])

      output$estimate <- paste0(
        round(output$estimate, dec),
        " (",
        round(output$conf.low, dec),
        ", ",
        round(output$conf.high, dec),
        ")"
      )
      output <- output[-1, -c(3:4, 6:7)]
      output$p.value <- round(output$p.value, pdec)
      if(predictorAsFactor == F) {
        output2 <- data.frame(Outcome = paste0(out, ":", levelsNum[2]),
                              output[, 1],
                              desc_y1,
                              desc_y2,
                              output[, 2],
                              output[, 3])
        names(output2) <- c("Outcome: Modeling Level",
                           "Predictor",
                           paste0('0:', " Mean (Std.Dev)"),
                           paste0('1:', " Mean (Std.Dev)"),
                           "Odds Ratio (95% CI)",
                           "p.value")
      }
      if(predictorAsFactor == T) {
        outputNullRow <- c(paste0(predictor, ":", levels(factor(data[, predictor]))[1]), NA, NA)
        names(outputNullRow) <- c("term", "estimate", "p.value")
        output <- rbind(outputNullRow, output)
        output2 <- data.frame(Outcome = paste0(out, ":", levelsNum[2]),
                              output[, 1],
                              desc_y1,
                              desc_y2,
                              output[, 2],
                              output[, 3])
        names(output2) <- c("Output",
                           "Predictor",
                           paste0(0, " N (%)"),
                           paste0(1, " N (%)"),
                           "Odds Ratio (95% CI)",
                           "p.value")
      }
      print(names(output2))
      return(output2)
    }
 fullOut <- do.call(rbind, lapply(outcomes, makeone))
 return(fullOut)
  }
