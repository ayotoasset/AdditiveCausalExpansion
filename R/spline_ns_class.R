ns_spline <- setRefClass("ns_spline",
                              fields = list(n_knots = "numeric", #internal knots
                                            m = "numeric", #cubic
                                            knots = "matrix", #knots placements for test data
                                            B = "matrix",   #design matrix
                                            dB = "matrix"), #derivative off design matrix
                              methods = list(
                                init = function(Z,nknots) {
                                  m = 4 # with intercept?
                                  n_knots <- nknots #equal to (n_knots + m) spline elements

                                  if(ncol(Z)>1){cat("More than one dimension of Z not supported (yet)")}
                                  #knots based on quantiles easier in one-dimension
                                  # equally spaced?
                                  knots <- sort(c(rep(min(Z),4),rep(max(Z),4),quantile(Z,probs = seq(n_knots)/(n_knots+1))))

                                  #generate matrix
                                  #require(splines)
                                  B <- splines::splineDesign(knots=knots, x=Z, ord=m, derivs=0)
                                  dB <- splines::splineDesign(knots=knots, x=Z, ord=m, derivs=1)


                                  #basis <- splineDesign(Aknots, x, 4)
                                  #const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2)) ?
                                  #if (!intercept) {
                                  #  const <- const[, -1, drop = FALSE]
                                  #  basis <- basis[, -1, drop = FALSE]
                                  #}
                                  #qr.const <- qr(t(const))
                                  #basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),
                                  #                                                   drop = FALSE])

                                },
                                predict = function(Znew){
                                  Bnew <- splines::splineDesign(knots=knots, x=Znew, ord=m, derivs=0)
                                  dBnew <- splines::splineDesign(knots=knots, x=Znew, ord=m, derivs=1)

                                  list(B = Bnew, dB = dBnew)
                                }


                              )
)
