set_initial_parameters <- function(p, B, n, y, Z) {
  Z <- as.matrix(Z)
  parameters <- matrix(c(0,0, # sigma, mu
                          -log(c(1, diag(t(Z) %*% Z) / n)), # use inner product to scale
                          rep(c(-0.1, rep(-0.1, B - 1)), p) # ARD parameters
                         )
                       )

  invisible(parameters)
}




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
