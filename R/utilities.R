plot_train_stats <- function(stats) {
  graphics::par(mfrow=c(1, 2))
  graphics::plot(stats[2, ], type="l", ylab="log Evidence", xlab="Iteration")
  graphics::plot(stats[1, ], type="l", ylab="training RMSE", xlab="Iteration")
  graphics::par(mfrow=c(1, 1))
}

set_optimizer <- function(optimizer, myKernel, learning_rate, momentum,
                          beta1, beta2,
                          norm.clip, clip.at) {
  if (optimizer=="Adam") myOptimizer <- optAdam$new(myKernel, lr = learning_rate, beta1 = beta1, beta2 = beta2,
                                                    norm.clip =  norm.clip, clip.at = clip.at)
  else if (optimizer=="Nadam") myOptimizer <- optNadam$new(myKernel, lr = learning_rate, beta1 = beta1, beta2 = beta2,
                                                           norm.clip =  norm.clip, clip.at = clip.at)
  else if (optimizer == "GD" || optimizer=="NAG") {
    if (optimizer == "GD") momentum <- 0.0
    myOptimizer <- optNesterov$new(myKernel, lr = learning_rate, momentum = momentum,
                                   norm.clip =  norm.clip, clip.at = clip.at)
  }
  invisible(myOptimizer)
}

set_basis <-function(basis, isuniv) {
  if ( (basis == "binary") || (basis == "linear")) myBasis <- linear_spline$new()
  else if (isuniv && (basis == "B")) myBasis <- B_spline$new()
  else if (isuniv && (basis == "square")) myBasis <- square_spline$new()
  else if (isuniv && (basis == "cubic")) myBasis <- ns_spline$new()
  else if (isuniv) myBasis <- ns_spline$new()
  invisible(myBasis)
}



