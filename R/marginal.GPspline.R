
#' Predict the the marginal response surface wrt Z using the Gaussian process fit.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param newX A matrix with new data. If not presented, using the training data.
#' @param newZ A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param marginal A logical

#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.

marginal <- function(x) UseMethod("marginal", x)

marginal.GPspline <- function(object,newX,newZ,causal=FALSE){
  predict(object,newX,newZ, marginal = TRUE,causal = causal)
}
