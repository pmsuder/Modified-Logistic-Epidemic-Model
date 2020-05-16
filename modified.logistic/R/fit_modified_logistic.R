#' Fitting Modified Logistic
#'
#' @description fit.modified.logistic is used to fit the modified logistic model as specified in the description of the package.
#'
#' @param I.vec Vector of values of the output variable - number of infected at given point in time in case of epidemic data.
#' @param t.vec Vector of values of the input variable - value of time in case of time series.
#' @param K.start Starting value of parameter K given to the optimization algorithm which adjusts the coefficients.
#' @param theta.start Starting value of parameter theta given to the optimization algorithm which adjusts the coefficients.
#' @param omega.start Starting value of parameter omega given to the optimization algorithm which adjusts the coefficients.
#' @param gamma.start Starting value of parameter gamma given to the optimization algorithm which adjusts the coefficients.
#' @param autoAdjustStartingVals IF TRUE the algorithm will first fit the model without the theta term and then use the obtained values as starting values for fitting of the full model. It is recommended to leave this TRUE.
#' @param autoAdjustedTheta.start Starting value of parameter theta in case the model is fitted using the other parameters from the reduced model fitting. It is recommended to keep this one close to 0.
#' @param adjustmentOmega.start Starting value of parameter omega given to the optimization algorithm which adjusts the coefficients in the reduced model.
#' @param adjustmentGamma.start Starting value of parameter gamma given to the optimization algorithm which adjusts the coefficients in the reduced model.
#' @param adjustmentK.start Starting value of parameter K given to the optimization algorithm which adjusts the coefficients in the reduced model.
#' @param maxiter A positive integer specifying the maximum number of iterations allowed for both reduced and full model fitting.
#' @param tol 	A positive numeric value specifying the tolerance level for the relative offset convergence criterion for both reduced and full model fitting.
#' @param minFactor A positive numeric value specifying the minimum step-size factor allowed on any step in the iteration.
#' @param adjustForMinCutoff If TRUE the function will delete the data points for which the value of I.vec is less than the minCutoff.
#' @param minCutoff If adjustForMinCutoff is TRUE the function will delete the data points for which the value of I.vec is less than the minCutoff.
#' @param adjust_t_ForMinCutoff If the function adjusted for the minCutoff and this parameter is TRUE, the valule of t will be adjusted so that t starts running from the first data point which was not deleted.
#' @param constFor_t_Adjustment This value will be added after the adjustment to t is made post miCutoff adjustment. It is recommended to keep this parameter at value greater than 0, so that t is never 0.
#'
#' @return Nonlinear regression model.
#' @export
#'
#' @examples fit.modified.logistic(I.vec = data$Cases, t.vec = data$Days, K.start = 200000, omega.start = 0.2, gamma.start = 2, autoAdjustedTheta.start = 1e-7, minCutoff = 500)
fit.modified.logistic <- function(I.vec, t.vec, K.start = 100000,
                                  theta.start = 0.1, omega.start = 0.1, gamma.start = 1,
                                  autoAdjustStartingVals = TRUE, autoAdjustedTheta.start = 1e-6,
                                  adjustmentOmega.start = -0.1, adjustmentGamma.start = 1, adjustmentK.start = 100000,
                                  maxiter = 10000, tol = 1e-4, minFactor = 1e-15,
                                  adjustForMinCutoff = TRUE, minCutoff = 0,
                                  adjust_t_ForMinCutoff = TRUE, constFor_t_Adjustment = 1) {


  data <- data.frame(I.vec, t.vec)

  if (adjustForMinCutoff == TRUE) {
    data <- data[(data$I.vec > minCutoff),]
    if (adjust_t_ForMinCutoff == TRUE) {
      data$t.vec <- data$t.vec - data$t.vec[1] + constFor_t_Adjustment
    }
  }

  control <- nls.control(maxiter = maxiter, tol = tol, minFactor = minFactor, warnOnly = TRUE)

  if (autoAdjustStartingVals == TRUE) {
    adjustmentK.start <- 2 * max(data$I.vec)
    model_0 <- nls(I.vec ~ K / (1 + exp((omega * t.vec) + gamma)),
                   start = list(omega = adjustmentOmega.start, gamma = adjustmentGamma.start,
                                K = adjustmentK.start), control = control)
    summary_0 <- summary(model_0)
    omega_0_hat <- summary_0$coefficients[1,1]
    gamma_0_hat <- summary_0$coefficients[2,1]
    K_0_hat <- summary_0$coefficients[3,1]

    omega.start <- omega_0_hat
    gamma.start <- gamma_0_hat
    K.start <- K_0_hat
    theta.start <- autoAdjustedTheta.start
  }


  model <- nls(I.vec ~ K / (1 + exp(theta * (t.vec ^ -1) + (omega * t.vec) + gamma)),
               start = list(omega = omega.start, gamma = gamma.start, theta = theta.start,
                            K = K.start), control = control)
  model
}
