#' Predict the the response surface and the marginal response wrt Z using the Gaussian process fit.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param newX A matrix with new data. If not presented, using the training data.
#' @param newZ A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param marginal A logical
#' @param causal logical

#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.


predict.GPspline <- function(object,newX,newZ, marginal = FALSE, causal = FALSE){
  if( length(unique(object$train_data$Z))==2 ) { isbinary <- TRUE }
    else { isbinary <- FALSE }
  #this function returns the prediction for the fitted Gaussian process
  if(missing(newX) || missing(newZ)){
    if(missing(newX)){
      cat("No newX given, using training X.\n")
      newX <- object$train_data$X #normalized
    }
    if(missing(newZ)){
      cat("No newZ given, using training Z.\n")
      newZ <- object$train_data$Z #normalized
    }
  }  else if( !missing(newX) && !missing(newZ)) {
    newX <- matrix(newX); newZ <- matrix(newZ)
    if(ncol(newX)!=ncol(object$train_data$X)){ stop("Error: Dimension mismatch of X with newX", call. = FALSE) }
    if((ncol(newZ)!=ncol(object$train_data$Z)) && (length(newZ)!=1) ){ stop("Error: Dimension mismatch of Z with newZ", call. = FALSE) }
    if(length(newZ)==1){ newZ <- rep(newZ,nrow(newX)) }

    #normalize the non-binary variables
    normalize_test(newX,newZ,object$moments)
  }

  #get appropriate basis
  if(marginal==FALSE){
    pred_list <- object$Kernel$predict(object$train_data$y,
                                       object$train_data$X,
                                       object$Spline$B,
                                       newX, object$Spline$testbasis(newZ)$B,
                                       object$moments[1,1],object$moments[1,2])
  } else {
    pred_list <- object$Kernel$predict_treat(object$train_data$y,
                                             object$train_data$X,
                                             object$Spline$B,
                                             newX, object$Spline$testbasis(newZ)$dB,
                                             object$moments[1,1],object$moments[1,2],(isbinary && causal))
  }
  pred_list
}








