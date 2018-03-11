#' Predict the the response surface and the marginal response wrt Z using the Gaussian process fit.
#'
#' @param object An "ace" object from the ace
#' @param newX A matrix with new data. If not presented, using the training data.
#' @param newZ A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param marginal A logical
#' @param causal logical

#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.


predict.ace <- function(object, newX, newZ, marginal = FALSE, causal = FALSE, normalize = TRUE){
  isbinary <- object$train_data$Zbinary
  px <- ncol(object$train_data$X)
  pz <- ncol(object$train_data$Z)
  #this function returns the prediction for the fitted Gaussian process

  if (missing(newX) && missing(newZ)){
    cat("No newX given, using training X.\n")
    newX.intern <- object$train_data$X #normalized
    cat("No newZ given, using training Z.\n")
    newZ.intern <- object$train_data$Z #normalized

  } else if (missing(newX) && !missing(newZ)){
    newX.intern <- as.matrix(rlang::duplicate(object$train_data$X)) #normalized

    if (length(newZ) == 1) newZ.intern <- matrix(rep(as.numeric(newZ.intern), nrow(newX.intern)), nrow(newX.intern), 1)
    else {
      #normalize the non-binary variables
      newZ.intern <- as.matrix(as.numeric(rlang::duplicate(newZ)))
      if(normalize) normalize_test(newX.intern, newZ.intern, object$moments)
      newX.intern <- object$train_data$X
    }
  } else if (!missing(newX) && missing(newZ)) {
    newX.intern <- as.matrix(rlang::duplicate(newX))
    if (!marginal) cat("No newZ given, using 0.\n")

    newZ.intern <- matrix(rep(0, nrow(newX.intern)), nrow(newX.intern), 1)

    #normalize the non-binary variables
    if(normalize) normalize_test(newX.intern, newZ.intern, object$moments)

    newZ.intern <- matrix(rep(0, nrow(newX.intern)), nrow(newX.intern), 1)

  }  else if (!missing(newX) && !missing(newZ)) {
    newX.intern <- as.matrix(rlang::duplicate(newX))
    newZ.intern <- as.matrix(rlang::duplicate(newZ))

    if (ncol(newX.intern) != ncol(object$train_data$X)) stop("Dimension mismatch of X with newX", call. = FALSE)
    #if ((ncol(newZ.intern) != ncol(object$train_data$Z)) && (length(newZ) != 1)) stop("Dimension mismatch of Z with newZ", call. = FALSE)
    if (length(newZ) == 1) newZ.intern <- matrix(rep(newZ.intern, nrow(newX.intern)), nrow(newX.intern), 1)

    #normalize the non-binary variables
    if(normalize) normalize_test(newX.intern, newZ.intern, object$moments)
  }

  if (any(abs(newX.intern) > 1)) cat("Some values of X outside the unit-circle of the training data.\n")

  #get appropriate basis
  if(!marginal) pred_list <- object$Kernel$predict(object$train_data$y,
                                                   object$train_data$X,
                                                   object$Basis$B,
                                                   newX.intern,
                                                   object$Basis$testbasis(newZ.intern)$B,
                                                   object$moments[1, 1],
                                                   object$moments[1, 2])
  else {
    test.basis <- object$Basis$testbasis(newZ.intern)
    pred_list <- object$Kernel$predict_marginal(y = object$train_data$y,
                                                X = object$train_data$X,
                                                Z = object$Basis$B,
                                                X2 = newX.intern,
                                                Z2 = test.basis$B,
                                                dZ2 = test.basis$dB,
                                                mean_y = object$moments[1, 1],
                                                std_y = object$moments[1, 2],
                                                std_Z = object$moments[(2 + px):(1 + px + pz), 2],
                                                calculate_ate = (isbinary && causal))
  }
  invisible(pred_list)
}

