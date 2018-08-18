#' Fit an Additive Causal Expansion model
#'
#' This function fits a varying coefficient model with Gaussian process priors on the coefficients. The number of additive elements is determined by the expansion of the treatment variable.
#' \itemize{
#' \item For binary treatments this corresponds to the linear model form:\cr
#'     ```y = m(x) + tau(x)  z + eps```
#' \item A continuous treatment with e.g. 3-rd order polynomial is fit using:\cr
#'     ```y = m(x) + g_1(x) * z + g_2(x) * z^2 + g_3(x) * z^3 + eps```
#' }
#'
#' @author Philip Pilgerstorfer <philip at pilgerstorfer.org>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib ace, .registration=TRUE
#'
#' @param y A numeric vector of the outcomes.
#' @param X A numeric matrix of confounder variables.
#' @param Z A vector or matrix (multivariate / tensor splines might be added in the future) Factors will be transformed to numeric using as.numeric. Hence, non-binary actors are discouraged especially when they are not ordinal.
#' @param pi (Optional) A vector of propenbsity score that can be used for binary Z when the propensity score function is assumed to be in a different RKHS than the potential outcome function. NOTE that prediction needs to be done using `newZ = (Z - pi)`. For refrence see Hahn et al. (2018).
#' @param kernel A string denoting the Mercer/Covariance Kernel (default: "SE"). Options:  "SE" - Squared exponential kernelwith ARD, and "Matern32" - the Matern 3/2 kernel with ARD.
#' @param basis A string (default: "binary"/"linear" if Z is binary and "ns", natural cubic spline, for continuous or discrete Z). Options: "linear" - 1rd order polynomial, "square" - 2rd order polynomial, "cubic" - 3rd order polynomial, "B" - B-spline, and "ns" - Natural cubic spline.
#' @param n.knots An integer denoting the number of internal knots of the spline of Z (default: 1). Ignored if the basis is not a spline.
#' @param optimizer A string selecting the gradient optimizer (default: "NAG"). Options: "GD" - gradient descent, "NAG" - Nesterov accelerated gradient descent, "Adam" - Adaptive moment estimation (Kingma & Ba, 2015), "Nadam" - Nesterov-accelerated Adam (Dozat, 2016).
#' @param maxiter  An integer defining the maximum number of iterations  of the gradient-based empirical Bayes optimization (default: 1000).
#' @param tol A numeric scalar defining the stopping tolerance for the gradient-based empirical Bayes optimization (default: 1e-4).
#' @param learning_rate An integer defining the learning rate for the gradient-based empirical Bayes optimization (default: 0.01).
#' @param beta1 A numeric scalar defining the ("first moment") learning parameter for the Adam and Nadam optimizers (default: 0.5, in NN often 0.9).
#' @param beta2 A numeric scalar defining the ("second moment") learning parameter for the Adam and Nadam optimizers (default: 0.0, in NN often 0.999).
#' @param momentum Momentum for the empirical Bayes optimization when using Nesterov. Equivalent to gradient descent ("GD") if momentum is 0.
#' (default: 0.0). Ignored for Nadam and Adam optimizers.
#' @param norm.clip A boolean scalar defining whether the gradients should be clipped based on their norm (default: TRUE if optimizer is GD or NAG, FALSE if optimizer is Adam or Nadam)
#' @param clip.at A numeric scalar defining the maximum L2-length of the gradients (default: 1). If the Euclidian length of the gradients is above this value, they will be proportionally shrunk to this maximum. Ignored if norm.clip is FALSE.
#' @param init.sigma A numeric scaler defining the initial value for the noise variance (default: NULL). If NULL, uses the estimated noise variance of a linear regression.
#' @param init.length_scale A numeric scalar defining the initial value of the length sclaes of the ARD kernel (default: 20). For large values the Kernel approximates a linear model, which can be a good initial guess, especially for discrete confounders.
#' @param plot_stats A boolean scalar to turn off the plotting of the likelihood and RMSE learning curves (default: TRUE).
#' @param init_iter An integer scalar determining how many pre-rpocessing iterations for the noise variance should be done. (Default: 0)
#' @param verbose A boolean scalar to turn off any text output of the function (default: FALSE).
#'
#' @return The function returns the fitted process, together with relevant training settings, as an "ace" class object.
#' Predictions of the outcome curves as well as the marginal response (over Z) can be obtained using the accompanying "prediction" S3 method; see examples.
#'
#' @references
#' * Dozat, T. (2016). "Incoproating Nesterov Momentum into Adam." Proceedings of the International Conference on Learning Representations
#' * Hahn, P. R., C. M. Carvalho, D. Puelz, and J. He (2018). "Regularization and Confounding in Linear Regression for Treatment Effect Estimation." Bayesian Analysis, Vol 13, No 1, pp. 163-182.
#' * Hill, J. (2011). "Bayesian Nonparametric Modeling for Causal Inference." Journal of Computational and Graphical Statistics, Vol. 20, No. 1, pp. 217–240.
#' * Kingma, D. and Ba, J. (2015). "Adam: A method for stochastic optimization." Proceedings of the International Conference on Learning Representations
#'
#' @examples
#' library(ace)
#' ## Example with binary treatment similar to Hill (2011)'s
#'
#' set.seed(1231)
#' n <- 300
#'
#' # generate treatment
#' Z <- rbinom(n, 1, 0.3)
#'
#' # generate confounder and exogenous variable
#' X <- matrix(NaN, n, 1)
#' X[Z==1, ] <- rnorm(sum(Z), mean = 30,sd = 10)
#' X[Z==0, ] <- rnorm(n - sum(Z), mean = 20, sd = 10)
#' E <- runif(n) # exogenous variable
#' X <- data.frame(X, E)
#'
#' # sort Confounder for visualizations
#' sort.idx <- sort(X[, 1], index.return = TRUE)$ix
#'
#' # define and draw the reponse function
#' y_truefun <- function(x, z) {
#'     mat <- matrix(NaN, length(z), 1)
#'     mat[z==0, 1] <- matrix(72 + 3 * (x[z == 0,1] > 0) * sqrt(abs(x[z == 0, 1])), sum(z == 0), 1)
#'     mat[z==1, 1] <- matrix(90 + exp(0.06 * x[z == 1, 2]), sum(z == 1), 1)
#'     c(mat)}
#' y0_true <- y_truefun(X, rep(0, n))
#' y1_true <- y_truefun(X, rep(1, n))
#' Y0 <- rnorm(n, mean = y0_true, sd = 1)
#' Y1 <- rnorm(n, mean = y1_true, sd = 1)
#' Y <- Y0 * (1 - Z) + Y1 * Z
#'
#' # run model
#' my.ace <- ace.train(Y, X, Z,
#'                     kernel = "SE", optimizer = "Nadam",
#'                     learning_rate = 0.005, maxiter = 1000)
#'
#' # print (sample) average treatment effect (ATE)
#' my.pred <- predict(my.ace, marginal = TRUE, return_average_treatments = TRUE)
#'
#' # predicted vs actual ATE
#' my.pred$ate
#' mean(y1_true - y0_true)
#'
#' # plot outcome/response curves
#' plot_ace(my.ace, 1, marginal = FALSE)
#'
#' # plot treatment curve
#' plot_ace(my.ace, 1, marginal = TRUE)
#'
#' # Check for robustness of ATE:
#' robust_treatment(my.ace, n.steps=5)
#'
#' ## Example with continuous Z
#'
#' # generate confounder and treatment (dosage)
#' set.seed(1234)
#' n2 <- 200
#' Xc <- data.frame(matrix(runif(2 * n2, min = 1, max = 2), n2, 2))
#' Zc <- rnorm(n2, Xc[,1] + exp(Xc[,1]), exp(Xc[,1]))
#' yc_truefun <- function(x, z) {as.matrix(20 * sin(10 * x[, 2] - 1.5)
#'                                         + (x[, 1] - 1.5) * abs(z^3 - 2 * z)/100)
#'                              }
#'
#' # generate response curve
#' yc_true <- yc_truefun(Xc, Zc)
#' Yc <- rnorm(n2, mean = yc_true, sd = 2)
#'
#' # train and predict model (in-sample)
#' my.ace <- ace.train(Yc, Xc, Zc,
#'                     optimizer = "Nadam",
#'                     learning_rate = 0.01,
#'                     basis = "ncs")
#' my.pred <- predict(my.ace)
#'
#' # plot a 3D surface of the outcome response with respect to the confounder (column 1)
#' plot_ace(my.ace, 1, plot3D = TRUE)
#'
#' # plot a 3D surface of the marginal response
#' plot_ace(my.ace, 1, marginal = TRUE, plot3D = TRUE)
#'
#' # plot of the 2D curve with only Z using the mean for each dimension of X
#' plot_ace(my.ace)
#'
#'@export
ace.train <- function(y, X, Z, pi,
                      kernel            = "SE",
                      basis             = "linear",
                      n.knots           = 1L,
                      optimizer         = "NAG",
                      maxiter           = 1000L,
                      tol               = 1e-4,
                      learning_rate     = 0.01,
                      beta1             = 0.9,
                      beta2             = 0.999,
                      momentum          = 0.0,
                      norm.clip         = (optimizer %in% c("GD", "NAG")),
                      clip.at           = 1,
                      init.sigma        = NULL, #Default:
                      init.length_scale = 20, #gives a roughly linear model (-> infty)
                      plot_stats        = TRUE,
                      init_iter         = 0L,
                      verbose=TRUE) {

  if (is.factor(y)) stop("y is not numeric. This package does not support classification tasks.")

  n  <- length(y)
  px <- ncol(X)
  y  <- matrix(y)

  if (is.factor(Z)) Z <- (as.numeric(Z) - 1)

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
    if (verbose) cat("Assuming binary Z\n")
    if (missing(pi)) basis <- "binary"
  } else if (verbose) cat("Non-Binary Z detected\n")

  if (!missing(pi) && isuniv) {
    pi.intern <- as.matrix(rlang::duplicate(pi))
    if (!isbinary) pi.intern <- (pi.intern - moments[2 + px, 1]) / moments[2 + px, 2]
    Z.intern <- Z.intern - pi.intern
    #isbinary <- FALSE
  }

  #### select chosen spline or appropriate based on data ###
  myBasis <- set_basis(basis, isuniv)
  # generate basis
  myBasis$trainbasis(Z.intern, n.knots, verbose) #binary "spline" discards n_knots

  #Set kernel parameters
  init_parameters = set_initial_parameters(px, myBasis$dim(), n,
                                           y, X.intern, Z.intern,
                                           init.sigma = init.sigma, init.length_scale = init.length_scale,
                                           verbose = verbose)
  #Gaussian Process kernel (only SE and Matern so far)
  if (kernel == "Matern32") myKernel <- KernelClass_Matern32_R6$new(px, myBasis$dim(), init_parameters, moments[1,2], verbose)
  else                      myKernel <- KernelClass_SE_R6$new(px, myBasis$dim(), init_parameters, moments[1,2], verbose)

  #initialize optimizer
  myOptimizer <- set_optimizer(optimizer, myKernel, learning_rate = learning_rate, momentum,
                               beta1 = beta1, beta2 = beta2,
                               norm.clip =  norm.clip, clip.at = clip.at)

  ### run iterations:
  stats <- matrix(0, 2, maxiter + 2) #Evidence and RMSE

  if ((init_iter > 0) & (class(myKernel) == "SqExpKernel")) { # not implemented for Matern
      myKernel$init_variance(init_iter, y, X.intern, myBasis$B, myOptimizer, verbose=verbose)
  }


  for(iter in 1:maxiter){
    stats[, iter + 1] <- myKernel$para_update(iter, y, X.intern, myBasis$B, myOptimizer, verbose=verbose)

    # prints it correctly:
    #stats[, iter + 1] <- myKernel$get_train_stats(y, X.intern, myBasis$B)

    change = abs(stats[2, iter + 1] - stats[2, iter])
    if ((change < tol) && (iter > 3)) {
      if (verbose) cat(sprintf("Stopped: change smaller than tolerance after %d iterations\n",
                               iter))
      break
    }
  }
  convergence.flag <- (iter < maxiter)

  stats[, iter + 2] <- myKernel$get_train_stats(y, X.intern, myBasis$B)

  if (verbose) if(!convergence.flag) cat("WARNING NO CONVERGENCE - Optimization stopped: maximum iterations reached\n")
  if (verbose) cat("Final training log Evidence: ", stats[2, iter + 2], "\n")

  stats <- stats[, 3:(iter + 2)]

  if (plot_stats) plot_train_stats(stats)

  #if(!missing(pi)) Z.intern <- (as.matrix(rlang::duplicate(Z)) - moments[2 + px, 1]) / moments[2 + px, 2]
  # reset graphics setting
  graphics::par(mfrow=c(1,1))
  invisible(structure(list(Kernel = myKernel, Basis = myBasis, OptimSettings = list(optim = optimizer,
                                                                            lr = learning_rate,
                                                                            momentum = momentum,
                                                                            beta1 = beta1,
                                                                            beta2 = beta2),
                           moments = moments,
                           train_data = list(y = y, X = X.intern, Z = as.matrix(Z.intern), Zbinary = isbinary),
                           train_stats = list(RIC_bias_corrected = !missing(pi),
                                              init.length_scale = init.length_scale,
                                              convergence = convergence.flag,
                                              final_evidence = utils::tail(stats[2,], n=1) )),
                 class = "ace"))
}
