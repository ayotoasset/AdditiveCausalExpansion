set_initial_parameters <- function(p, B, n, y, X, Z, init.sigma, init.length_scale=20, verbose=FALSE) {
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  if (!missing(init.sigma)) {
    #least square
    QR <- qr.default(cbind(X,Z,1))
    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
    init.sigma <- log( t(y) %*% (diag(n) - tcrossprod(Q)) %*% y / (n - 1) )
  }
  # RF:
  #rf_fit <- ranger::ranger(y ~ ., data=data.frame(y=y, X=X, Z=Z))
  #init_sigma <- log( sum( (y - rf_fit$predictions)^2) / (n-1)  )

  if (verbose) cat("Initial noise variance: ", exp(init.sigma), "\n")
  parameters <- matrix(c(init.sigma, #sigma
                         0, # mu
                         -log(rep(1, B)),# kernel-scales
                         log(rep(init.length_scale, B * p)) # ARD length-scales
                         )
                       )

  invisible(parameters)
}
