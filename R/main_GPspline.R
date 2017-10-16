#' Fit a spline model with Gaussian process coefficients
#'
#' @param y A numeric vector
#' @param X A numeric vector or matrix
#' @param Z A vector or matrix (multivariate / tensor splines might be added in the future) Factors will be transformed to numeric using as.numeric. Hence, non-binary actors are discouraged especially when they are not ordinal.
#' @param kernel A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and Matern 5/2 kernel
#' @param spline A string (default: "ns" natural cubic spline for continuous or discrete Z and "binary" if Z is binary (factor)
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z
#' @param myoptim A string (default: "GD" gradient descent). Other options are "Nesterov" (accelerated gradient/momentum), "Adam", and "Nadam" (Nesterov-Adam).
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


GPspline <- function(y,X,Z,kernel = "SE",spline="ns",n.knots=3,myoptim = "GD",maxiter=1000,tol=1e-4,learning_rate=0.01,beta1=0.9,beta2=0.999,momentum=0.0){
  if(class(y)=="factor") stop("y is not numeric. This package does not support classification tasks.")
  n = length(y); px = ncol(X);

  if((class(Z) == "matrix") || (class(Z)== "data.frame")) { pz = ncol(Z); }
  else if(length(c(Z))==n ){ pz=1; }
  else { stop("Dimension/filetype of Z invalid.\n") }

  X <- as.matrix(as.numeric(X))
  if(class(Z)=="factor") {Z <- (as.numeric(Z)-1) }
  Z <- as.matrix(as.numeric(Z))
  if( !all(dim(X)==c(n,px)) ) stop("Dimension of X not correct. Use the observations as rows and variables as columns and check the number of observations with respect to y.\n")
  if( !all(dim(Z)==c(n,pz)) ) stop("Dimension of Z not correct. Use the observations as rows and variables as columns and check the number of observations with respect to y.\n")

  #normalize variables
  moments <- normalize_train(y,X,Z)

  #check whether Z is univariate
  if( pz > 1 ) isuniv = FALSE else isuniv = TRUE
  #check whether Z is binary
  if( (length(unique(c(Z))) > 2) ) {cat("Non-Binary Z detected"); isbinary = FALSE}
  else {cat("Binary Z detected\n"); isbinary = TRUE; spline = "binary"}

  #### select chosen spline or appropriate based on data ###
  if( isuniv && isbinary && (spline == "binary")) {
    cat("Using binary spline\n");  mySpline <- binary_spline$new()  }
  else if( isuniv && (spline == "B") ) {
    cat("Using B-spline\n"); mySpline <- B_spline$new()  }
  else if( isuniv ){
    cat("Using NC-spline\n"); mySpline <- ns_spline$new()  }

  #generate basis
  mySpline$trainbasis(Z,n.knots) #binary "spline" discards n_knots

  #Gaussian Process kernel (only SE implemented)
  if (kernel == "Matern") {
    #myKernel <- KernelClass_Matern$new(px = px,pz = pz)
  } else if(kernel == "Polynomial") {
    #myKernel <- KernelClass_Poly$new(px = px,pz = pz)
  } else {
    myKernel <- KernelClass_SE$new()
  }

  #initialize Kernel parameters given the spline basis dimension (e.g.: binary:2, ncs: n_knots+3)
  myKernel$parainit(y,p=px,mySpline$dim())

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

  if(iter == maxiter) cat("Optimization stopped: maximum iterations reached\n")

  stats[,iter+2] = myKernel$get_train_stats(y,X,Z)

  graphics::par(mfrow=c(1,2))
  graphics::plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iteration")
  graphics::plot(stats[1,3:(iter+2)],type="l",ylab="training RMSE",xlab="Iteration")

  structure(list(Kernel = myKernel, spline=mySpline,
                 moments=moments,
                 train_data=list(y = y,X = X,Z = Z)),
                 class = "GPspline")

}

