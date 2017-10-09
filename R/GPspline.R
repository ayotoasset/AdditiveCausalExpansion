#' Fit a spline model with Gaussian process coefficients
#'
#' @param y A vector
#' @param X A data.frame
#' @param Z A data.frame
#' @param kernelfun A string (default: SE) -- no other kernels implemented
#' @param myoptim A string (default: Gradient Descent -- GD)
#' @param maxiter (default: 5000)
#' @param tol (default: 1e-4)
#' @param learning_rate (default: 0.01)
#' @param beta1 (default: 0.9)
#' @param beta2 (default: 0.999)
#' @param momentum (default: 0.0)
#' @return The function returns the fitted process as a GPspline class object
#' @examples
#' #Generate data
#' n = 120
#' Z = rbinom(n, 1, 0.3)
# 'X1 = runif(sum(Z), min = 20, max = 40)
#' X0 = runif(n-sum(Z), min = 20, max = 40)
#' X = data.frame(matrix(NaN,n,1))
#' X[Z==1,] = X1; X[Z==0,] = X0
#' y0_true = as.matrix(72 + 3 * sqrt(X))
#' y1_true = as.matrix(90 + exp(0.06 * X))
#' Y0 = rnorm(n, mean = y0_true, sd = 1)
#' Y1 = rnorm(n, mean = y1_true, sd = 1)
#' Y = Y0*(1-Z) + Y1*Z
#' mystump <- CausalStump(Y,X,Z)

GPspline <- function(y,X,Z,myoptim = "Nadam",maxiter=5000,tol=1e-4,learning_rate=0.01,beta1=0.9,beta2=0.999,momentum=0.0){
  #check_inputs(y,X,z);

  n = length(y); px = ncol(X); pz = ncol(Z);

  #normalize variables
  norm_ret = norm_variables(y,X,Z)
  moments = norm_ret$moments; y = norm_ret$y; X = norm_ret$X; Z = norm_ret$Z;

  #cat("\nFitting the Gaussian process:\n")
  myKernel = KernelClass_SE$new(px = px,pz = pz)
  myKernel$parainit(y);

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

  ### write the loop in C++ at one point together with the optimizer initialization
  for(iter in 1:maxiter){
    stats[,iter+1] = myKernel$para_update(iter,y,X,Z,myOptimizer)

    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol) && (iter > 3)){ cat( sprintf("Stopped: change smaller than tolerance after %d iterations",iter)); break; }
  }

  if(iter == maxiter){ cat("Stopped: maximum iterations reached") }

  stats[,iter+2] = myKernel$get_train_stats(y,X,Z)

  graphics::par(mfrow=c(1,2))
  graphics::plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iteration")
  graphics::plot(stats[1,3:(iter+2)],type="l",ylab="training RMSE",xlab="Iteration")

  structure(list(Kernel = myKernel,moments=moments,train_data=list(y=y,X=X,Z=Z)), class = "GPspline")

}
