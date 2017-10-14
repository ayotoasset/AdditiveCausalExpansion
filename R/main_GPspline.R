#' Fit a spline model with Gaussian process coefficients
#'
#' @param y A vector
#' @param X A data.frame
#' @param Z A vector (multivariate / tensor splines might be added in the future)
#' @param kernelfun A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and Matern 5/2 kernel
#' @param spline A string (default: "ns" natural cubic spline for continuous or discrete Z and "binary" if Z is binary (factor)
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z
#' @param myoptim A string (default: Gradient Descent -- GD)
#' @param maxiter (default: 5000) Maximum number of iterations of the empirical Bayes optimization
#' @param tol (default: 1e-4) Stopping tolerance for the empirical Bayes optimization
#' @param learning_rate (default: 0.01) Learning rate for the empirical Bayes optimization
#' @param beta1 (default: 0.9) Learning parameter ("first moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param beta2 (default: 0.999) Learning parameter ("second moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param momentum (default: 0.0) Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' @return The function returns the fitted process as a GPspline class object. Predictions can be obtained using the corresponding S3 methods "prediction" and "marginal".
#' The latter is the predicted curve with a differentiated basis of Z.
#' @examples
#' #Example
#' #Generate data
#' n = 120
#' Z = rbinom(n, 1, 0.3)


GPspline <- function(y,X,Z,kernel = "SE",spline="ns",n.knots=3,myoptim = "Nadam",maxiter=5000,tol=1e-4,learning_rate=0.01,beta1=0.9,beta2=0.999,momentum=0.0){
  #check_inputs(y,X,z); needs revision

  n = length(y); px = ncol(X); pz = ncol(Z);

  #normalize variables
  #########  #########  #########  #########  #########  ######### WRITE IN CPP
  norm_ret = norm_variables(y,X,Z)
  moments = norm_ret$moments; y = norm_ret$y; X = norm_ret$X; Z = norm_ret$Z;
  # make X an appropriately sized matrix

  #check whether Z is univariate
  if( pz > 1 ) isuniv = FALSE else isuniv = TRUE
  #check whether Z is binary
  if( (unique(c(Z)) > 2) ) {cat("Non-Binary Z detected"); isbinary = FALSE}
  else {cat("Binary Z detected");isbinary = TRUE; spline = "binary"}

  #### select chosen spline or appropriate based on data ###
  if( isuniv && isbinary && (spline == "binary")) {
    cat("Binary spline selected\n");  mySpline <- binary_spline$new()  }
  else if( isuniv && (spline == "B") ) {
    cat("B-spline selected\n");  mySpline <- B_spline$new()  }
  else if( isuniv ){
    cat("NC-splineselected\n");   mySpline <- ns_spline$new()  }

  #generate basis
  mySpline$trainbasis(Z,n_knots) #binary "spline" discards n_knots

  #Gaussian Process kernel (only SE implemented)
  if (kernel = "Matern") {
    #myKernel <- KernelClass_Matern$new(px = px,pz = pz)
  } else if(kernel = "Polynomial") {
    #myKernel <- KernelClass_Poly$new(px = px,pz = pz)
  } else {
    myKernel <- KernelClass_SE$new(px = px,pz = pz)
  }

  #initialize Kernel parameters given the spline basis dimension (binary:2, ncs: n_knots+3)
  myKernel$parainit(y,mySPline$dim)

  #initialize optimizer
  if((myoptim=="Adam") || (myoptim=="Nadam")){
    if(myoptim=="Adam") {
      myOptimizer = optAdam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    } else {
      myOptimizer = optNadam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    }
  } else if(myoptim=="GD" || myoptim=="Nesterov"){
    if(myoptim=="GD"){ momentum=0.0 }
    myOptimizer = optNestorov$new(lr = learning_rate, momentum = momentum)
  }

  #set optimization variables
  myOptimizer$initOpt(myKernel);

  stats = matrix(0,2,maxiter+2) #Evidence and RMSE

  ### write the loop in C++ at one point together with the optimizer initialization ?
  for(iter in 1:maxiter){
    stats[,iter+1] = myKernel$para_update(iter,y,X,mySpline$B,myOptimizer)

    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol) && (iter > 3)){ cat( sprintf("Stopped: change smaller than tolerance after %d iterations",iter)); break; }
  }

  if(iter == maxiter) cat("Stopped: maximum iterations reached")

  stats[,iter+2] = myKernel$get_train_stats(y,X,Z)

  graphics::par(mfrow=c(1,2))
  graphics::plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iteration")
  graphics::plot(stats[1,3:(iter+2)],type="l",ylab="training RMSE",xlab="Iteration")

  structure(list(Kernel = myKernel, spline=mySpline,
                 moments=moments,
                 train_data=list(y = y,X = X,Z = Z)),
                 class = "GPspline")

}
