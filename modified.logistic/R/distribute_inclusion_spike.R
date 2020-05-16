#' Get Rid of an Inculsion Spike in the Data
#'
#' @description distribute.inclusion.spike is used when the epidemiological data has a spike in the number of cases at one point, which was caused by widening the diagnostics criteria or inclusion of a large number of individuals who have been known to be infected, but were not officially reported, as it was the case, for example, with France during the COVID-19 epidemic in 2020, which on one day included all the infections from elderly care facilities which were not counted previously in the reported cases.
#'
#' @param data Data frame with predictor variables and outcome
#' @param spikePointIndex Integer value of the index of the t.vec vector at which the inclusion spike occurs in I.vec.
#' @param additionKnown Boolean value. If TRUE, then the program takes the additionNum as the number of cases added at the inclusion spike. If FALSE, the program estimates the additionNum and uses a smoothing procedoure to distribute the inclusion spike proportionally over the previous entries of I.vec.
#' @param additionNum Number of cases added which form the inclusion spike.They will be distributed proportionally over the previous values of I.vec.
#'
#' @return Data frame with I.vec and t.vec containing the smoothed data.
#' @export
#'
#' @examples distribute.inclusion.spike(I.vec = data$Cases, t.vec = data$Days, spikePointIndex = 22, additionKnown = FALSE)
distribute.inclusion.spike <- function(I.vec, t.vec, spikePointIndex, additionKnown = TRUE, additionNum = 0) {
  data <- data.frame(I.vec, t.vec)
  if (additionKnown == FALSE) {
    x <- sqrt(data$I.vec[spikePointIndex - 1] * data$I.vec[spikePointIndex + 1])
    additionNum <- data$I.vec[spikePointIndex] - x
  }
  data$I.vec[spikePointIndex] <- data$I.vec[spikePointIndex] - additionNum
  for (i in 1:spikePointIndex) {
    data$I.vec[i] <- data$I.vec[i] + (additionNum * data$I.vec[i] / data$I.vec[spikePointIndex])
  }

  data
}
