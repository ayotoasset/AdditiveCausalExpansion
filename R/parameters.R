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
                         -log(rep(1, ncol(Z) + 1)),# kernel-scales
                         log(rep(init.length_scale, B*p)) # ARD length-scales
                         )
                       )

  invisible(parameters)
}

#generate_random_parameters <- function(p, B) {
#  parameters <- matrix(c(runif(1, min=-2, max=2), #sigma
#                         rnorm(1), # mu
#                         runif(B, min=-2, max=2),
#                         rep(runif(1, min=-1, max=1), B*p) # ARD parameters
#                         )
#                       )
#  invisible(parameters)
#}
#
#gen_random_para_list <- function(len, p, B){
#  output <- list(generate_random_parameters(p, B))
#  if (len >= 2) {
#    for (i in 2:len) {
#      output[[i]] <- generate_random_parameters(p, B)
#    }
#  }
#
#  output
#}


#n <- length(y)
#parameters <<- matrix(c(2 * B,0, # sigma, mu
#                         (c(1, rep(1, B - 1))), # use inner product to scale
#                        rep(c(-0.1, rep(0.1, B - 1)), p) # ARD parameters
#                        )
#                     )
#
#
#
#
#
#
