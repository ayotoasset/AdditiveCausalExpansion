set_initial_parameters <- function(p, B, n, y, X, Z, init.sigma=NULL, init.length_scale=20, verbose=FALSE) {
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  if (is.null(init.sigma)) {
    #least square
    QR <- qr.default(cbind(X,Z,1))
    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
    init.sigma <- log( t(y) %*% (diag(n) - tcrossprod(Q)) %*% y / (n - 1) )
  }
  else {
    init.sigma <- log(init.sigma)
  }
  # RF:
  #rf_fit <- ranger::ranger(y ~ ., data=data.frame(y=y, X=X, Z=Z))
  #init_sigma <- log( sum( (y - rf_fit$predictions)^2) / (n-1)  )

  if (verbose) cat("Initial noise variance: ", exp(init.sigma), "\n")
  parameters <- matrix(c(init.sigma, # "log( sigma)"
                         0, # mu
                         -log(rep(1, B)),# kernel-scales
                         log(rep(init.length_scale, B * p)) # ARD length-scales
                         )
                       )

  invisible(parameters)
}
