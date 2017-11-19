#' Fit a VCM spline model with Gaussian process coefficients
#'
#' @param y A numeric vector
#' @param X A numeric vector or matrix
#' @param Z A vector or matrix (multivariate / tensor splines might be added in the future) Factors will be transformed to numeric using as.numeric. Hence, non-binary actors are discouraged especially when they are not ordinal.
#' @param kernel A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and Matern 5/2 kernel
#' @param spline A string (default: "ns" natural cubic spline for continuous or discrete Z and "binary" if Z is binary (factor)
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z
#' @param myoptim A string (default: "GD" gradient descent). Other options are "Nesterov" (accelerated gradient/momentum), "Adam", and "Nadam" (Nesterov-Adam).
#' @param maxiter  (default: 5000) Maximum number of iterations of the empirical Bayes optimization
#' @param tol (default: 1e-4) Stopping tolerance for the empirical Bayes optimization
#' @param learning_rate (default: 0.001) Learning rate for the empirical Bayes optimization
#' @param beta1 (default: 0.9) Learning parameter ("first moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param beta2 (default: 0.999) Learning parameter ("second moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param momentum (default: 0.0) Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' @return The function returns the fitted process as a GPspline class object. Predictions can be obtained using the corresponding S3 methods "prediction" and "marginal".
#' The latter is the predicted curve with a differentiated basis of Z.
#' @examples
#' #Example replicating CausalStump with binary uni-variate Z
#' #Generate data
#' library(GPspline)
#' set.seed(1231)
#' n <- 120
#' Z <- rbinom(n, 1, 0.3)
#' X1 <- rnorm(sum(Z), mean = 30,sd = 10)
#' X0 <- rnorm(n-sum(Z), mean = 20, sd = 10)
#' X <- matrix(NaN,n,1)
#' X[Z==1,] <- X1; X[Z==0,] <- X0
#' sort.idx <- sort(X,index.return=TRUE)$ix
#' y_truefun <- function(x,z) {mat0 <- matrix(72 + 3 * (x[z==0,]>0)*sqrt(abs(x[z==0,])),sum(z==0),1); mat1 <- matrix(90 + exp(0.06 * x[z==1,]),sum(z==1),1) ; mat <- matrix(NaN,length(z),1); mat[z==0,1] <- mat0; mat[z==1,1] <- mat1; c(mat) }
#' y0_true <- y_truefun(X,rep(0,n))
#' y1_true <- y_truefun(X,rep(1,n))
#' Y0 <- rnorm(n, mean = y0_true, sd = 1)
#' Y1 <- rnorm(n, mean = y1_true, sd = 1)
#' Y <- Y0*(1-Z) + Y1*Z
#' my.GPS <- GPspline.train(Y,X,Z,myoptim="GD",learning_rate=0.0001)
#' #print (sample) average treatment effect (ATE)
#' predict(my.GPS,marginal=TRUE,causal=TRUE)$ate_map
#' #true ATE
#' mean(y1_true-y0_true)
#' #plot response curves
#' plot(my.GPS,1,marginal=FALSE,truefun=y_truefun)
#' #plot treatment curve
#' treat_truefun <- function(x) {y_truefun(x,rep(1,nrow(x))) - y_truefun(x,rep(0,nrow(x)))}
#' plot(my.GPS,1,marginal=TRUE,truefun=treat_truefun)
#'
#' #continuous Z
#' set.seed(1234)
#' n2 <- 200
#' X2 <- matrix(runif(n2, min = 1, max = 2))
#' X3 <- matrix(runif(n2, min = 1, max = 2))
#' X <- data.frame(X2,X3)
#' Z2 <- rnorm(n2, exp(X2)-10, 1)
#' y_truefun <- function(x,z) {as.matrix(10*x + (x-1.5) * ((z+9)^2 - 2*z))}
#' marg_truefun <- function(x,z) {as.matrix( (x-1.5) * (2*(z+9) - 2))}
#' y2_true <- y_truefun(X2,Z2)
#' Y2 <- rnorm(n2, mean = y2_true, sd = 2)
#' marg_true <- marg_truefun(X2,Z2)
#' my.GPS <- GPspline.train(Y2,X2,Z2,myoptim="GD",learning_rate = 0.0001,spline="cubic",n.knots=1)
#' plot(my.GPS,1,marginal=F,plot3D=FALSE,truefun=y_truefun)
#' my.pred <- predict(my.GPS)
#' #quality of prediction
#' plot(Y2,my.pred$map); abline(0,1,lty=2)
#' #comparison with the true curve:
#' plot(my.GPS,1,marginal=FALSE,truefun = y_truefun)
#' #plotting of the marginal curve:
#' plot(my.GPS,1,marginal=TRUE)
#' #plot of the 2D curve with only Z
#' plot(my.GPS,marginal=FALSE,truefun=y_truefun)
#'

GPspline.train <- function(y,X,Z,kernel = "SE",spline="ns",n.knots=1,myoptim = "GD",maxiter=1000,tol=1e-4,learning_rate=0.001,beta1=0.9,beta2=0.999,momentum=0.0){

  if(class(y)=="factor") stop("y is not numeric. This package does not support classification tasks.")
  n <- length(y); px <- ncol(X);
  y <- matrix(data.table::copy(y));
  X.intern <- as.matrix(data.table::copy(X))
  if(class(Z)=="factor") {Z <- (as.numeric(Z)-1) }
  Z.intern <- matrix(as.numeric(data.table::copy(Z)))

  if((class(Z) == "matrix") || (class(Z)== "data.frame")) { pz <- ncol(Z); }
  else if(length(c(Z))==n ) { pz <- 1; }
  else { stop("Dimension/filetype of Z invalid.\n") }


  if( !all(dim(X.intern)==c(n,px)) ) stop("Dimension of X not correct. Use the observations as rows and variables as columns and check the number of observations with respect to y.\n")
  if( !all(dim(Z.intern)==c(n,pz)) ) stop("Dimension of Z not correct. Use the observations as rows and variables as columns and check the number of observations with respect to y.\n")

  #normalize variables
  moments <- normalize_train(y,X.intern,Z.intern)

  #check whether Z is univariate
  if( pz > 1 ) isuniv = FALSE else isuniv = TRUE
  #check whether Z is binary
  if( isbinary <- (length(unique(c(Z))) == 2) ) {cat("Binary Z detected\n"); spline = "binary"}
  else {cat("Non-Binary Z detected\n"); isbinary = FALSE}

  #Gaussian Process kernel (only SE and Matern so far)
  if (kernel == "Matern12") {
    cat("Using Matern 1/2 kernel\n"); myKernel <- KernelClass_Matern12$new()
  } else if (kernel == "Matern32") {
    cat("Using Matern 3/2 kernel\n"); myKernel <- KernelClass_Matern32$new()
  } else if (kernel == "Matern52") {
    cat("Using Matern 5/2 kernel\n"); myKernel <- KernelClass_Matern52$new()
  } else if(kernel == "SEiso"){
    cat("Using isotropic SE kernel (treats each basis dimension equally)\n"); myKernel <- KernelClass_SE_iso$new()
  } else if(kernel == "SEnonadd"){
    cat("Using non-additive SE kernel (regular GP)\n"); myKernel <- KernelClass_SE_nonadd$new()
    cat("Only linear spline supported for regular GP features\n"); spline <- "linear"
  } else {
    cat("Using SE kernel\n"); myKernel <- KernelClass_SE$new()
  }

  #### select chosen spline or appropriate based on data ###
  if( isuniv && ((isbinary && (spline == "binary")) || (spline=="linear"))) {
    cat("Using binary/linear-basis\n");  mySpline <- linear_spline$new()  }
  else if( isuniv && (spline == "B") ) {
    cat("Using B-spline\n"); mySpline <- B_spline$new()  }
  else if( isuniv && (spline == "square") ) {
    cat("Using square-basis\n"); mySpline <- square_spline$new() }
  else if( isuniv && (spline == "cubic") ) {
    cat("Using cubic-basis\n"); mySpline <- ns_spline$new(); n.knots=0;}
  else if( isuniv ){
    cat("Using NC-spline\n"); mySpline <- ns_spline$new()  }

  #generate basis
  mySpline$trainbasis(Z.intern,n.knots) #binary "spline" discards n_knots

  #initialize Kernel parameters given the spline basis dimension (e.g.: binary:2, ncs: n_knots+3)
  myKernel$parainit(y,p=px,mySpline$dim(),mySpline$B)

  #initialize optimizer
  if((myoptim=="Adam") || (myoptim=="Nadam")){
    if(myoptim=="Adam") {
      myOptimizer = optAdam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    } else {
      myOptimizer = optNadam$new(lr = learning_rate, beta1 = beta1, beta2 = beta2)
    }
  } else if(myoptim=="Newton"){
    myOptimizer = optNewton$new(lr = learning_rate, momentum = momentum)
  } else if(myoptim=="GD" || myoptim=="Nesterov"){
    if(myoptim=="GD"){ momentum=0.0 }
    myOptimizer = optNesterov$new(lr = learning_rate, momentum = momentum)
  }

  #set optimization variables
  myOptimizer$initOpt(myKernel);

  stats = matrix(0,2,maxiter+2) #Evidence and RMSE

  ### write the loop in C++ at one point together with the optimizer initialization ?
  for(iter in 1:maxiter){
    stats[,iter+1] = myKernel$para_update(iter,y,X.intern,mySpline$B,myOptimizer)

    change = abs(stats[2,iter+1] - stats[2,iter])
    if((change < tol) && (iter > 3)){ cat( sprintf("Stopped: change smaller than tolerance after %d iterations\n",iter)); break; }
  }

  if(iter == maxiter) cat("Optimization stopped: maximum iterations reached\n")

  stats[,iter+2] = myKernel$get_train_stats(y,X.intern,mySpline$B)

  graphics::par(mfrow=c(1,2))
  graphics::plot(stats[2,3:(iter+2)],type="l",ylab="log Evidence",xlab="Iteration")
  graphics::plot(stats[1,3:(iter+2)],type="l",ylab="training RMSE",xlab="Iteration")
  graphics::par(mfrow=c(1,1))

  structure(list(Kernel = myKernel, Spline=mySpline,OptimSettings = list(optim=myoptim,
                                                                         lr=learning_rate,
                                                                         momentum = momentum,
                                                                         beta1 = beta1,
                                                                         beta2 = beta2,
                                                                         iter = iter),
                 moments=moments,
                 train_data=list(y = y,X = X.intern,Z = Z.intern, Zbinary = isbinary)),
                 class = "GPspline")
}

#' Fit GPspline, data.frame wrapper of GPspline.train
#' @param formula Formula of the form \code{y ~ x | z} where the \code{x} variables determine the coefficients of the (spline) basis expansion of \code{z}.
#' @param data A data.frame containing the variables in the formula.
#' @param kernel A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and Matern 5/2 kernel
#' @param spline A string (default: "ns" natural cubic spline for continuous or discrete Z and "binary" if Z is binary (factor)
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z
#' @param myoptim A string (default: "GD" gradient descent). Other options are "Nesterov" (accelerated gradient/momentum), "Adam", and "Nadam" (Nesterov-Adam).
#' @param maxiter  (default: 5000) Maximum number of iterations of the empirical Bayes optimization
#' @param tol (default: 1e-4) Stopping tolerance for the empirical Bayes optimization
#' @param learning_rate (default: 0.001) Learning rate for the empirical Bayes optimization
#' @param beta1 (default: 0.9) Learning parameter ("first moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param beta2 (default: 0.999) Learning parameter ("second moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param momentum (default: 0.0) Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' @return GPspline object that can be used for prediction and plotting
#' @examples
#' #continuous Z
#' set.seed(1234)
#' n2 <- 300
#' df <- data.frame(x = runif(n2, min = 1, max = 2))
#' df$x2 <- runif(n2, min = -1, max = 1)
#' df$z = rnorm(n2, exp(df$x)-14, 1)
#' y_truefun <- function(x,z) {as.matrix(sqrt(x[,1]) + x[,2] *3 * ((z+8)^2 - 2*z))}
#' y2_true <- y_truefun(df[,c("x","x2")],df$z)
#' df$y <- rnorm(n2, mean = y2_true, sd = 1)
#' my.GPS <- GPspline(y ~ x + x2 | z,data=df,myoptim="GD",learning_rate = 0.0001,spline="ns",n.knots=2)
#' my.pred <- predict(my.GPS)
#' plot(df$y,my.pred$map); abline(0,1,lty=2)
#' #prediction of the curve
#' plot(my.GPS,"x2",marginal=FALSE,plot3D=TRUE)
#' #difference to the true marginal curve:
#' marg_truefun <- function(x,z) {as.matrix(sqrt(x[,1]) + x[,2] *3 * (2*(z+8) - 2))}
#' plot(my.GPS,"x2",marginal=TRUE,show.observations=TRUE,truefun=marg_truefun)

GPspline <- function(formula,data,kernel = "SE",spline="ns",n.knots=1,myoptim = "GD",maxiter=1000,tol=1e-4,learning_rate=0.001,beta1=0.9,beta2=0.999,momentum=0.0){
  myformula <- Formula(formula)
  data <- model.frame(myformula,data)

  y <- data[[attr(attr(data, "terms"),"response")]]
  X <- model.matrix(update(formula(myformula,lhs=0,rhs=1), ~ . + 0),data) #no intercept
  Z <- model.matrix(update(formula(myformula,lhs=0,rhs=2), ~ . + 0),data) #no intercept

  out_object <- GPspline.train(y,X,Z,kernel = kernel,spline=spline,n.knots=n.knots,myoptim = myoptim,maxiter=maxiter,tol=tol,
           learning_rate=learning_rate,beta1=beta1,beta2=beta2,momentum=momentum)

  colnames(out_object$train_data$y) <- all.vars(formula(myformula,lhs=1,rhs=0))
  colnames(out_object$train_data$X) <- all.vars(formula(myformula,lhs=0,rhs=1))
  colnames(out_object$train_data$Z) <- all.vars(formula(myformula,lhs=0,rhs=2))

  out_object
}
