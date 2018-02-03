B_spline <- R6::R6Class("B_spline",
                        cloneable = FALSE,
                        class = FALSE,
                        portable = FALSE,
                        public = list(
                          B = NULL,
                          dB = NULL,
                          dim = function(){
                            (ncol(B)+1); #basis plus nuisance term m()
                          },
                          trainbasis = function(Z, n_knots, m = 4) {
                            cat("Using B-spline\n")
                            #internal knots (within unit circle)
                            IKnots <- quantile(Z,probs = seq(n_knots)/(n_knots+1))
                            B <<- splines2::bSpline(Z,
                                                    knots = IKnots,
                                                    Boundary.knots = c(-1, 1), #unit circle is boundary
                                                    degree = m - 1, #cubic
                                                    intercept = FALSE)
                            dB <<- splines2::dbs(Z,
                                                 derivs = 1L,
                                                 knots = IKnots,
                                                 Boundary.knots = c(-1, 1), #unit circle is boundary
                                                 degree = m - 1, #cubic
                                                 intercept = FALSE)
                            invisible(B)
                          },
                          testbasis = function(Znew){
                            #get new basis matrices
                            if(!missing(Znew)){
                              a <- c(list(x = Znew), attributes(B)[c("degree", "knots","Boundary.knots", "intercept")])
                              Bnew <- do.call(splines2::bSpline, a)
                              a <- c(list(x = Znew), attributes(dB)[c("degree", "knots","Boundary.knots", "intercept", "derivs")])
                              dBnew <- do.call(splines2::dbs, a)
                            } else {
                              #return stored matrices
                              Bnew <- B
                              dBnew <- dB
                            }
                            invisible(list(B = Bnew, dB = dBnew))
                          }
                        )
)


