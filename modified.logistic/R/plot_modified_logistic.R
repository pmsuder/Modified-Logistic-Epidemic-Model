#' Plot the Fitted Model Anlong with Data Points
#'
#' @description plot.modified.logistic is used to plot the modified logistic model as specified in the description of the package together with the data points which were used to fit it.
#' @param model Model to be plotted.
#' @param I.vec Vector of y-coordinates of data points to be plotted.
#' @param t.vec Vector of x-coordinates of data points to be plotted.
#' @param numIntervals Number of intervals into which we divide the scope of the plot to provide for a smoothly-looking curve. It is recommended to keep this one large.
#' @param xlim Vector with two values - lower and upper limit for the x-coordinate in the scope of the plot.
#' @param ylimVector Vector with two values - lower and upper limit for the y-coordinate in the scope of the plot.
#'
#' @return
#' @export
#'
#' @examples plot.modified.logistic(model, data$Cases, data$Days, xlim = c(1, 200), ylim = c(0, 1000000))
plot.modified.logistic <- function(model, I.vec, t.vec, numIntervals = 100000, xlim, ylim) {

  x_vals <- seq(xlim[1] + ((xlim[2] - xlim[1]) / numIntervals), xlim[2], ((xlim[2] - xlim[1]) / numIntervals))
  y_vals <- predict(model, list(I.vec = x_vals), type = "response")
  plot(I.vec ~ t.vec, xlim = xlim, ylim = ylim)
  lines(x_vals, y_vals)
}
