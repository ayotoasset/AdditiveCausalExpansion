NewtonStep <- function(object,maxiter=1,tol,lr=object$OptimSettings$lr,momentum=0.0){
  #create backup of parameters
  backup_parameters <- object$Kernel$parameters

  #initialize Newton optimizer
  myOptimizer = optNewton$new(lr = lr, momentum = momentum);
  myOptimizer$initOpt(object$Kernel);

  #old stats:
  stats = matrix(0,2,1); oldstats=stats
  oldstats = object$Kernel$get_train_stats(object$train_data$y,object$train_data$X,object$Spline$B)
  originalstats=oldstats

  #update and new stats
  if(!missing(tol)){
    mytol <- tol+1
    iter <- object$OptimSettings$iter
    while(mytol>tol){
      oldstats <- object$Kernel$para_update(100,object$train_data$y,object$train_data$X,object$Spline$B,myOptimizer,printevery=101)
      stats <- object$Kernel$get_train_stats(object$train_data$y,object$train_data$X,object$Spline$B)
      mytol <- abs(stats[2]- oldstats[2]); #oldstats = stats
      iter <- iter+1
      cat("\r ",iter,"| Tolerance: ",mytol, "")
    }

    cat("\n")
  } else {
    for(iter in seq(maxiter)+object$OptimSettings$iter){
      stats <- object$Kernel$para_update(iter,object$train_data$y,object$train_data$X,object$Spline$B,myOptimizer,printevery=1)
      mytol <- abs(stats[2]- oldstats[2]); oldstats=stats
    }
    stats <- object$Kernel$get_train_stats(object$train_data$y,object$train_data$X,object$Spline$B)
    mytol <- stats[2]-oldstats[2]
    object$OptimSettings$iter <- iter
  }
  cat("After",iter-object$OptimSettings$iter <- iter,"steps:\n")
  cat("Final evidence: ", stats[2], "| Error:", stats[1],"| Tolerance: ", mytol,"\n")

  object$OptimSettings$iter <- iter

  if(originalstats[2]>stats[2]){
    cat("No improvement, maybe try a smaller learning rate\nKeeping original parameter values\n")
    object$Kernel$parameters <- backup_parameters
  }

  object
}
