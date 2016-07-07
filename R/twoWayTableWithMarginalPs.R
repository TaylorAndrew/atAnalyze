#Get the means and sd for a continuous variable 'var', as it's stratified by both
#columnGroup and rowGroup.
#Get the p-values compared var across colGroup at each level of rowGroup and vice
#versa
#' twoWayTableWithMarginalPs
#'
#' @param data data.frame contianing all data
#' @param columnGroup Grouping variable to stratify across columns
#' @param rowGroup Grouping variable to stratify across rows
#' @param var Outcome variable, numeric
#' @param columnPaired Logical, whether or not the column grouping variable should be treated as paired
#' @param rowPaired Logical, whether or not the row grouping variable should be treated as paired
#' @param aovPairedID If paired, id variable name.
#'
#' @return data.frame with all output
#' @export
#'
#' @examples
#' # Needs and example
twoWayTableWithMarginalPs <- function(data, columnGroup, rowGroup, var,
                                      columnPaired = FALSE, rowPaired = FALSE,
                                      aovPairedID) {
  Table_Plus_rowPsList <- atAnalyze::subset_t.test(data = data,
                         var = var,
                         grouping = columnGroup,
                         subsetVar = rowGroup,
                         paired = columnPaired,
                         aovPairedVar = aovPairedID,
                         includeTukey = TRUE)
  Table_Plus_rowPs <- Table_Plus_rowPsList$output
  TableComps <- Table_Plus_rowPsList$comps
  col <- atAnalyze::subset_t.test(data = data,
                         var = var,
                         grouping = rowGroup,
                         subsetVar = columnGroup,
                         paired = rowPaired,
                         aovPairedVar = aovPairedID,
                       includeTukey = TRUE)
  colPs <- col$output$pvalue
  colTuks <- col$comps
  colP_line <- data.frame(NA,
                          t(colPs),
                          NA)
  makeChar <- function(i) {
    col <- as.character(Table_Plus_rowPs[,i])
    return(col)
  }
  nms <- names(Table_Plus_rowPs)
  Table_Plus_rowPs <- data.frame(do.call(cbind,
                                         lapply(1:length(Table_Plus_rowPs[1,]),
                                            makeChar)), stringsAsFactors = F)
  # Table_Plus_rowPs <- as.character(Table_Plus_rowPs)
  names(colP_line) <- nms
  names(Table_Plus_rowPs) <- nms
  out <- rbind(Table_Plus_rowPs, colP_line)
  out <- list(Output = out, Comparisons = rbind(TableComps, colTuks))
  return(out)
}
