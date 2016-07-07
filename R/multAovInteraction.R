#' multAovInteraction
#'
#' multAovInteraction provides descriptive summary of multiple depedent variables, stratified by two grouping variables, along with p-values from aov() models testing the main effects and interaction between the two grouping variables.
#'
#' @param data Input dataset
#' @param grp_var1 Variable name of variable 1
#' @param grp_var2 Variable name of variable 2
#' @param varlist Vector of variable names to be used as the DV(s)
#' @param round Number of decimals for descriptive statistics
#' @param pround Number of decimals for p-values
#' @return data.frame containing the output summary
#' @export
#'
#' @examples
#' #df<-data.frame(grp1=sample(c(1,2),100,replace=TRUE),
#' #                grp2=sample(c(1,2,3),100,replace=TRUE),
#' #                var1=rnorm(100,32,12),
#' #                var2=runif(100),
#' #                var3=rnorm(100))
#' #multAovInteraction(df, grp_var1 = "grp1", grp_var1 = "grp2", varlist = c("var2","var3","var1"))
multAovInteraction <- function(data,
                             grp_var1,
                             grp_var2,
                             varlist,
                             round = 2,
                             pround = 3) {
    do_one <- function(var) {
    rowlevs <- ifelse(length(levels(factor(data[,grp_var1]))) < 3,3,
                      length(levels(factor(data[,grp_var1]))))
    mat <- matrix(nrow = rowlevs,
                  ncol = (length(levels(factor(
                    data[,grp_var2]
                  ))) + 4))
    levs1 <- sort(levels(factor(data[,grp_var1])))
    levs2 <- sort(levels(factor(data[,grp_var2])))
    do_cell <- function(i,j) {
      data_small <-
        as.numeric(unlist(subset(
          data,data[,grp_var1] == levs1[i] &
            data[,grp_var2] == levs2[j],select =
            var
        )))
      mat[i,j + 2] <<-
        paste0(round(mean(data_small,na.rm = TRUE),digits = round),
               " ? ",
               round(sd(data_small,na.rm = TRUE),digits = round))
    }
    for (i in 1:length(levs1)) {
      for (j in 1:length(levs2)) {
        do_cell(i,j)
      }
    }
    mat[1:length(levs1),2] <- levs1
    mat[1,1] <- var

    pvals <- round(summary(aov(data[,var] ~ factor(data[,grp_var1]) * factor(data[,grp_var2])))[[1]][["Pr(>F)"]]
                   ,pround)

    mat[1,(length(mat[1,]) - 1):length(mat[1,])] <- c(grp_var1,pvals[1])
    mat[2,(length(mat[1,]) - 1):length(mat[1,])] <- c(grp_var2,pvals[2])
    mat[3,(length(mat[1,]) - 1):length(mat[1,])] <-
      c("Interaction",pvals[3])


    mat[is.na(mat)] <- " "
    output <- data.frame(mat)
    names(output) <-
      c("Variable",grp_var1,paste(grp_var2,levs2),"Pval Effect","Pvalue")
    return(output)
  }
  if (length(varlist) > 1)
    output_combined <- do.call(rbind,lapply(varlist,do_one))
  if (length(varlist) == 1)
    output_combined <- do_one(varlist)
    return(output_combined)
}
