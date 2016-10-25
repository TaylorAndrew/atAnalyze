#' repeatedMeasuresROC
#'
#' repeatedMeasuresROC provides a table of sensitivities and specificities, as well as the AUC using Wilcoxon non-parametric approach, both in line with Liu and Wu (2003)
#'
#' @param glmerModel an lme4::glmer() model
#'
#' @return a list containing: ROC_table and AUC
#' @export
#'
#' @examples
#' hdp <- read.csv("http://www.ats.ucla.edu/stat/data/hdp.csv")
#' hdp <- within(hdp, {
#'   Married <- factor(Married, levels = 0:1, labels = c("no", "yes"))
#'   DID <- factor(DID)
#'   HID <- factor(HID)
#' })
#' m <- glmer(remission ~ IL6 + CRP + CancerStage + LengthofStay + Experience + (1 | DID),
#'            data = hdp,
#'            family = binomial,
#'            control = glmerControl(optimizer = "bobyqa"),
#'            nAGQ = 10)
#'            out <- repeatedMeasuresROC(m)
#'plot(x = 1-out$ROC_table$spec, y = out$ROC_table$sens)
#'abline(0,1)
#'(AUC = out$AUC)
repeatedMeasuresROC <-  function(glmerModel) {
## Get predicted probabilities from glmer model
preds <- predict(glmerModel, type = "response")

## Create a funciton that will get sensitivity and specificity for a single cut point:
getSensSpec <- function(cutpoint) {
    predBasedCut <- ifelse(preds <= cutpoint, 0, 1) # Get predicted probs based on cutpoint
    predBasedCut <-  factor(predBasedCut, levels = c(0, 1)) # Just in case all 1 or 0, make both levels
    freqs <- data.frame(table(predBasedCut, glmerModel@frame[, names(glmerModel@frame)[1]])) # Get frequecies
    sens <- freqs[4, 3] / (freqs[4, 3] + freqs[3, 3]) # Compute sensitivity (TP / TP + FN)
    spec <-  freqs[1, 3] / (freqs[1, 3] + freqs[2, 3]) # Compute specificity (TN / TN + FP)
    return(data.frame(cut = cutpoint,
                      sens = sens,
                      spec = spec))
}
## Run 'getSensSpec' for every cut point from 0 to 1, interating by 0.005
out <- do.call(rbind, lapply(seq(0, 1, by = 0.005), getSensSpec))
## Analysis for AUC:
# Run wilcox.test of rank differences in predicted probabilities between actual positive and negative observations
wi = wilcox.test(preds[glmerModel@frame[, names(glmerModel@frame)[1]]==1],preds[glmerModel@frame[, names(glmerModel@frame)[1]]==0])
# Get the wilcoxon statistic
w = wi$statistic
# Compute the AUC as the wilcoxon statistic divided by the number of positive pairs
AUC = w/(length(preds[glmerModel@frame[, names(glmerModel@frame)[1]]==1])*length(preds[glmerModel@frame[, names(glmerModel@frame)[1]]==0]))
names(AUC) <- 'AUC'
# Put both outputs together:
output = list(ROC_table = out, AUC = AUC)
# And return it
return(output)
}





