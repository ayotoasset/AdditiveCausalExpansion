#' Predict the response surfaces using a Gaussian process or Student-t process fit.
#'
#' @param object An "GPspline" object from the GPspline function
#' @param X A matrix with new data. If not presented, using the training data.
#' @param z A vector, matrix or scalar with new treatment data. If it is a scalar, it predicts using the same value for all observations. If missing, it uses the training data.
#' @param pscore A vector with the propensity score. Throws an error if use is different to the CausalStump fit.
#' @param nsampling (optional) A number that overwrites the number of samples used for the TP-kernel sampling.
#' @return Returns the MAP and the 95 percent credible interval of the fitted process as a list.


predict.GPspline <- function(object,newX,newZ, marginal = FALSE){
  #this function returns the prediction for the fitted Gaussian process
  if(missing(newX) || missing(Z)){
    if(missing(newX)){
      cat("No newX given, using training X.\n")
      newX <- object$train_data$X #normalized
    }
    if(missing(newZ)){
      cat("No newZ given, using training Z.\n")
      newZ <- object$train_data$Z #normalized
    }
  }  else if( !missing(newX) && !missing(Z)) {

    if(ncol(newX)!=ncol(object$train_data$X)){ stop("Error: Dimension mismatch of X with newX", call. = FALSE) }
    if((ncol(newZ)!=ncol(object$train_data$Z)) && (length(z)!=1)){ stop("Error: Dimension mismatch of Z with newZ", call. = FALSE) }
    if(length(newZ)==1){ newZ <- rep(newZ,nrow(newX)) }

    #normalize the non-binary variables
    normalize_test(newX,newZ,object$moments)
  }

  #get appropriate basis
  if(marginal==FALSE){ Zbasis <- object$Spline$testbasis(newZ)$B  }
  else { Zbasis <- object$Spline$testbasis(newZ)$dB }

  #remaining kernel calculations using the kernel class method

  pred_list = object$Kernel$predict(object$train_data$y,
                                    object$train_data$X,
                                    object$Spline$B,
                                    newX, Zbasis,
                                    object$moments[1,1],object$moments[1,2])

  list(map = pred_list$map,ci = pred_list$ci)
}





