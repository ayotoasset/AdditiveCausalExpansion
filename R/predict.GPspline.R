#' Predict the the response surface and the marginal response wrt Z using the Gaussian process fit.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param newX A matrix with new data. If not presented, using the training data.
#' @param newZ A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param marginal A logical
#' @param causal logical

#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.


predict.GPspline <- function(object,newX,newZ, marginal = FALSE, causal = FALSE){
  isbinary <- object$train_data$Zbinary
  #this function returns the prediction for the fitted Gaussian process
  if(missing(newX) && missing(newZ)){
    cat("No newX given, using training X.\n")
    newX.intern <- object$train_data$X #normalized

    cat("No newZ given, using training Z.\n")
    newZ.intern <- object$train_data$Z #normalized

  } else if(!missing(newX) && missing(newZ)){
    newX.intern <- as.matrix(data.table::copy(newX));
    cat("No newZ given, using 0.\n")
    newZ.intern <- matrix(rep(0,nrow(newX.intern)),nrow(newX.intern),1)

    #normalize the non-binary variables
    normalize_test(newX.intern,newZ.intern,object$moments)

  }  else if( !missing(newX) && !missing(newZ)) {
    newX.intern <- as.matrix(data.table::copy(newX)); newZ.intern <- as.matrix(data.table::copy(newZ))

    if(ncol(newX.intern)!=ncol(object$train_data$X)){ stop("Error: Dimension mismatch of X with newX", call. = FALSE) }
    if((ncol(newZ.intern)!=ncol(object$train_data$Z)) && (length(newZ)!=1) ){ stop("Error: Dimension mismatch of Z with newZ", call. = FALSE) }
    if(length(newZ)==1){ newZ.intern <- matrix(rep(newZ.intern,nrow(newX.intern)),nrow(newX.intern),1) }
    #normalize the non-binary variables
    normalize_test(newX.intern,newZ.intern,object$moments)

  }

  if(any(abs(newX.intern)>1)){
    cat("Some values of X outside the unit-circle of the training data.\n")
    apply(abs(newX.intern),2,max)

  }


  #get appropriate basis
  if(marginal==FALSE){
    pred_list <- object$Kernel$predict(object$train_data$y,
                                       object$train_data$X,
                                       object$Spline$B,
                                       newX.intern, object$Spline$testbasis(newZ.intern)$B,
                                       object$moments[1,1],object$moments[1,2])
  } else {
    if(class(object$Kernel)=="SqExpNonAddKernel"){
      pred_list <- object$Kernel$predict_marginal(object$train_data$y,
                                                  object$train_data$X,
                                                  object$Spline$B,
                                                  newX.intern, newZ.intern,
                                                  object$moments[1,1],object$moments[1,2],(isbinary && causal))
    } else {
      pred_list <- object$Kernel$predict_marginal(y = object$train_data$y,
                                                  X = object$train_data$X,
                                                  Z = object$Spline$B,
                                                  X2 = newX.intern,
                                                  dZ2 = object$Spline$testbasis(newZ.intern)$dB,
                                                  mean_y = object$moments[1,1],
                                                  std_y = object$moments[1,2],
                                                  std_Z = object$moments[dim(object$moments)[1],2],
                                                  calculate_ate = (isbinary && causal))
    }
  }
  pred_list
}








