#' Fit a VCM basis model with Gaussian process priors on the coefficients
#'
#' @param y A numeric vector
#' @param X A numeric vector or matrix
#' @param Z A vector or matrix (multivariate / tensor splines might be added in the future) Factors will be transformed to numeric using as.numeric. Hence, non-binary actors are discouraged especially when they are not ordinal.
#' @param kernel A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and "Matern32", the Matern 3/2 kernel with ARD
#' @param basis A string (default: "binary" if Z is binary (factor) and "ns", natural cubic spline, for continuous or discrete Z
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z. Ignored if basis not a spline.
#' @param optimizer A string (default: "Nadam" Nesterov-accelerated Adam). Other options are "GD" (gradient descent), "Nesterov" (accelerated gradient/momentum), "Adam".
#' @param maxiter  (default: 1000) Maximum number of iterations of the empirical Bayes optimization
#' @param tol (default: 1e-4) Stopping tolerance for the empirical Bayes optimization
#' @param learning_rate (default: 0.001) Learning rate for the empirical Bayes optimization
#' @param beta1 (default: 0.9) Learning parameter ("first moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param beta2 (default: 0.999) Learning parameter ("second moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param momentum (default: 0.0) Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' @return The function returns the fitted process as a GPspline class object. Predictions can be obtained using the corresponding S3 methods "prediction" and "marginal".
#' The latter is the predicted curve with a differentiated basis of Z.
#' @examples
#' ## Example replicating CausalStump with binary uni-variate Z
#' # Generate data
#' library(ace)
#' set.seed(1231)
#' n <- 300
#' Z <- rbinom(n, 1, 0.3)
#' X1 <- rnorm(sum(Z), mean = 30,sd = 10)
#' X0 <- rnorm(n-sum(Z), mean = 20, sd = 10)
#' X2 <- runif(n)
#' X <- matrix(NaN, n, 1)
#' X[Z==1, ] <- X1
#' X[Z==0, ] <- X0
#' X <- data.frame(X, X2)
#' sort.idx <- sort(X[,1], index.return = TRUE)$ix
#' y_truefun <- function(x,z) {
#'     mat0 <- matrix(72 + 3 * (x[z == 0,1] > 0) * sqrt(abs(x[z == 0,1])), sum(z == 0), 1)
#'     mat1 <- matrix(90 + exp(0.06 * x[z == 1,2]), sum(z == 1), 1)
#'     mat <- matrix(NaN, length(z), 1)
#'     mat[z==0, 1] <- mat0
#'     mat[z==1, 1] <- mat1
#'     c(mat)
#'     }
#' y0_true <- y_truefun(X, rep(0, n))
#' y1_true <- y_truefun(X, rep(1, n))
#' Y0 <- rnorm(n, mean = y0_true, sd = 1)
#' Y1 <- rnorm(n, mean = y1_true, sd = 1)
#' Y <- Y0 * (1-Z) + Y1 * Z
#' pihatz = nnet(X, Z, size=10)$fitted.values
#' # train model:
#' my.GPS <- ace.train(Y, X, Z, pi = pihatz, kernel="SE", optimizer = "GD", learning_rate = 0.002, maxiter=2000)
#' # print (sample) average treatment effect (ATE)
#' predict(my.GPS, marginal = TRUE, causal = TRUE)$ate
#' predict(my.GPS, marginal = TRUE, causal = TRUE)$att
#' robust_treatment(my.GPS)
#' abline(a=mean(y1_true - y0_true), b=0)
#' #true ATE
#' mean(y1_true - y0_true)
#' # plot response curves
#' plot(my.GPS, 1, marginal = TRUE)
#' # plot treatment curve
#' treat_truefun <- function(x) {y_truefun(x, rep(1, nrow(x))) - y_truefun(x, rep(0, nrow(x)))}
#' plot(my.GPS, 1, marginal = TRUE, truefun = treat_truefun)
#'
#' ## Example with continuous Z
#' set.seed(1234)
#' n2 <- 200
#' X2 <- matrix(runif(n2, min = 1, max = 2))
#' X3 <- matrix(runif(n2, min = 1, max = 2))
#' X <- data.frame(X2, X3)
#' Z2 <- rnorm(n2, exp(X2) - 10, 1)
#' y_truefun <- function(x, z) {as.matrix(10 * x[,1] + (x[,1] - 1.5) * ((z + 9)^2 - 2 * z))}
#' marg_truefun <- function(x, z) {as.matrix( (x[,1] - 1.5) * (2 * (z + 9) - 2))}
#' y2_true <- y_truefun(X2, Z2)
#' Y2 <- rnorm(n2, mean = y2_true, sd = 2)
#' marg_true <- marg_truefun(X2, Z2)
#' my.GPS <- ace.train(Y2, X, Z2,
#'                     optimizer = "Nadam",
#'                     learning_rate = 0.005,
#'                     basis = "ncs")
#' my.pred <- predict(my.GPS)
#' # plot quality of prediction:
#' plot(Y2, my.pred$map)
#' abline(0, 1, lty = 2, col = "blue", lwd = 2)
#' # comparison with the true curve:
#' plot(my.GPS, 1, truefun = y_truefun, plot3D = TRUE)
#' # plotting of the marginal curve:
#' plot(my.GPS, 1, marginal=TRUE)
#' # plot of the 2D curve with only Z
#' plot(my.GPS, truefun = y_truefun)
#'

ace.train <- function(y, X, Z, pi,
                      kernel          = "SE",
                      basis           = "linear",
                      n.knots         = 1,
                      optimizer       = "Nadam",
                      maxiter         = 1000,
                      tol             = 1e-4,
                      learning_rate   = 0.01,
                      beta1           = 0.9,
                      beta2           = 0.999,
                      momentum        = 0.0,
                      norm.clip       = (optimizer!="Nadam") | (optimizer!="Adam")  ,
                      clip.at         = 1,
                      init_parameters = set_initial_parameters(px, myBasis$dim(), n, y, Z)) {



  if (is.factor(y)) stop("y is not numeric. This package does not support classification tasks.")

  n  <- length(y)
  px <- ncol(X)
  y  <- matrix(y)

  if (is.factor(Z)) Z <- (as.numeric(Z)-1)

  # create manual copy of variable since calling cpp functions with reference for normalization
  Z.intern <- as.matrix(rlang::duplicate(Z))
  X.intern <- as.matrix(rlang::duplicate(X))

  if (is.matrix(Z) || is.data.frame(Z)) pz <- ncol(Z)
  else if (length(c(Z)) == n ) pz <- 1
  else stop("Dimension/filetype of Z invalid.\n")

  if (!all(dim(X.intern) == c(n, px))) stop("Dimension of X not correct. Use the observations
                                            as rows and variables as columns and check the number
                                            of observations with respect to y.\n")
  if (!all(dim(Z.intern) == c(n, pz))) stop("Dimension of Z not correct. Use the observations as
                                            rows and variables as columns and check the number of
                                            observations with respect to y.\n")

  #normalize variables
  moments <- normalize_train(y, X.intern, Z.intern)

  #check whether Z is univariate
  isuniv <- (pz == 1)

  #check whether Z is binary
  isbinary <- (moments[(2 + px):(1 + pz + px), 3] == 1)
  if (all(isbinary)) {
    cat("Assuming binary Z\n")
    if (missing(pi)) basis <- "binary"
  } else cat("Non-Binary Z detected\n")

  if (!missing(pi) && isuniv) {
    pi.intern <- as.matrix(rlang::duplicate(pi))
    if (!isbinary) pi.intern <- (pi.intern - moments[2 + px, 1]) / moments[2 + px, 2]
    Z.intern <- Z.intern - pi.intern
    #isbinary <- FALSE
  }

  #### select chosen spline or appropriate based on data ###
  if ( (basis == "binary") || (basis == "linear")) myBasis <- linear_spline$new()
  else if (isuniv && (basis == "B")) myBasis <- B_spline$new()
  else if (isuniv && (basis == "square")) myBasis <- square_spline$new()
  else if (isuniv && (basis == "cubic")) myBasis <- ns_spline$new()
  else if (isuniv) myBasis <- ns_spline$new()

  #generate basis
  myBasis$trainbasis(Z.intern, n.knots) #binary "spline" discards n_knots

  #Gaussian Process kernel (only SE and Matern so far)
  if (kernel == "Matern32") myKernel <- KernelClass_Matern32_R6$new(px, myBasis$dim(), init_parameters)
  else                      myKernel <- KernelClass_SE_R6$new(px, myBasis$dim(), init_parameters)

  #initialize optimizer
  if (optimizer=="Adam") myOptimizer = optAdam$new(myKernel, lr = learning_rate, beta1 = beta1, beta2 = beta2,
                                                   norm.clip =  norm.clip, clip.at = clip.at)
  else if (optimizer=="Nadam") myOptimizer = optNadam$new(myKernel, lr = learning_rate, beta1 = beta1, beta2 = beta2,
                                                          norm.clip =  norm.clip, clip.at = clip.at)
  else if (optimizer == "GD" || optimizer=="Nesterov") {
    if(optimizer == "GD") momentum = 0.0
    myOptimizer = optNesterov$new(myKernel, lr = learning_rate, momentum = momentum,
                                  norm.clip =  norm.clip, clip.at = clip.at)
  }


  ### run
  stats = matrix(0, 2, maxiter + 2) #Evidence and RMSE

  ### write the loop in C++ at one point together with the optimizer initialization ?
  for(iter in 1:maxiter){
    stats[, iter+1] = myKernel$para_update(iter, y, X.intern, myBasis$B, myOptimizer)

    change = abs(stats[2, iter + 1] - stats[2, iter])
    if ((change < tol) && (iter > 3)) {
      cat( sprintf("Stopped: change smaller than tolerance after %d iterations\n",
                   iter))
      break
      }
  }
  convergence.flag = (iter < maxiter)
  if(!convergence.flag) cat("Optimization stopped: maximum iterations reached\n")

  stats[, iter+2] = myKernel$get_train_stats(y, X.intern, myBasis$B)

  graphics::par(mfrow=c(1, 2))
  graphics::plot(stats[2, 3:(iter + 2)], type="l", ylab="log Evidence", xlab="Iteration")
  graphics::plot(stats[1, 3:(iter + 2)], type="l", ylab="training RMSE", xlab="Iteration")
  graphics::par(mfrow=c(1, 1))

  #if(!missing(pi)) Z.intern <- (as.matrix(rlang::duplicate(Z)) - moments[2 + px, 1]) / moments[2 + px, 2]

  invisible(structure(list(Kernel = myKernel, Basis = myBasis, OptimSettings = list(optim = optimizer,
                                                                            lr = learning_rate,
                                                                            momentum = momentum,
                                                                            beta1 = beta1,
                                                                            beta2 = beta2,
                                                                            iter = iter),
                 moments = moments,
                 train_data = list(y = y, X = X.intern, Z = Z.intern, Zbinary = isbinary),
                 RIC_bias_corrected = !missing(pi),
                 convergence = convergence.flag),
                 class = "ace"))
}

#' Fit GPspline, data.frame wrapper of GPspline.train
#' @param formula Formula of the form \code{y ~ x | z} where the \code{x} variables determine the coefficients of the (spline) basis expansion of \code{z}.
#' @param data A data.frame containing the variables in the formula.
#' @param kernel A string (default: "SE" Squared exponential with ARD) -- has no effect, might include (ARD) polynomial and Matern 5/2 kernel
#' @param basis A string (default: "ns" natural cubic spline for continuous or discrete Z and "binary" if Z is binary (factor)
#' @param n.knots An integer denoting the  umber of internal knots of the spline of Z
#' @param optimizer A string (default: "GD" gradient descent). Other options are "Nesterov" (accelerated gradient/momentum), "Adam", and "Nadam" (Nesterov-Adam).
#' @param maxiter  (default: 5000) Maximum number of iterations of the empirical Bayes optimization
#' @param tol (default: 1e-4) Stopping tolerance for the empirical Bayes optimization
#' @param learning_rate (default: 0.001) Learning rate for the empirical Bayes optimization
#' @param beta1 (default: 0.9) Learning parameter ("first moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param beta2 (default: 0.999) Learning parameter ("second moment") for the empirical Bayes optimization when using Adam or Nadam optimizers
#' @param momentum (default: 0.0) Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' @return GPspline object that can be used for prediction and plotting
#' @examples
#' ## Example continuous treatment Z
#' set.seed(1234)
#' n2 <- 300
#' df <- data.frame(x = runif(n2, min = 1, max = 2))
#' df$x2 <- runif(n2, min = -1, max = 1)
#' df$z = rnorm(n2, exp(df$x) - 14, 1)
#' y_truefun <- function(x, z) {as.matrix(sqrt(x[, 1]) + x[, 2] *3 * ((z + 8)^2 - 2 * z))}
#' y2_true <- y_truefun(df[, c("x", "x2")], df$z)
#' df$y <- rnorm(n2, mean = y2_true, sd = 2)
#' my.GPS <- ace(y ~ x + x2 | z,
#'               data = df,
#'               optimizer = "GD",
#'               learning_rate = 0.0001,
#'               basis = "ns",
#'               n.knots = 2)
#' my.pred <- predict(my.GPS)
#' plot(df$y, my.pred$map)
#' abline(0, 1, lty=2)
#' # prediction of the curve
#' plot(my.GPS, "x2", plot3D = TRUE)
#' # difference to the true marginal curve:
#' # marg_truefun <- function(x,z) {as.matrix(sqrt(x[,1]) + x[,2] *3 * (2*(z+8) - 2))}
#' # plot(my.GPS,"x2",marginal=TRUE,show.observations=TRUE,truefun=marg_truefun)

ace <- function(formula, data,
                     kernel = "SE",
                     basis = "linear",
                     n.knots = 1,
                     optimizer = "Nadam",
                     maxiter = 1000,
                     tol = 1e-4,
                     learning_rate = 0.001,
                     beta1 = 0.9,
                     beta2 = 0.999,
                     momentum = 0.0,
                     norm.clip     = TRUE,
                     clip.at       = 1) {
  myformula <- Formula::Formula(formula)
  data <- stats::model.frame(myformula, data)

  y <- data[[attr(attr(data, "terms"), "response")]]
  X <- stats::model.matrix(update(formula(myformula, lhs=0, rhs=1), ~ . + 0), data) #no intercept
  Z <- stats::model.matrix(update(formula(myformula, lhs=0, rhs=2), ~ . + 0), data) #no intercept

  out_object <- ace.train(y, X, Z,
                             kernel = kernel,
                             basis = basis,
                             n.knots = n.knots,
                             optimizer = optimizer,
                             maxiter = maxiter,
                             tol = tol,
                             learning_rate = learning_rate,
                             beta1 = beta1,
                             beta2 = beta2,
                             momentum = momentum)

  colnames(out_object$train_data$y) <- all.vars(formula(myformula, lhs=1, rhs=0))
  colnames(out_object$train_data$X) <- all.vars(formula(myformula, lhs=0, rhs=1))
  colnames(out_object$train_data$Z) <- all.vars(formula(myformula, lhs=0, rhs=2))

  invisible(out_object)
}
