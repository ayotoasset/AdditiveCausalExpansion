B_spline <- setRefClass("B_spline",
                              fields = list(#n_knots = "numeric", #number of internal knots (input)
                                            #m = "numeric", #cubic = 4
                                            #myknots = "numeric", #knots placements for test data
                                            B = "matrix",   #design matrix and bSpline2 class
                                            dB = "matrix"), #derivative off design matrix
                              methods = list(
                                basisdim = function(){
                                  (ncol(B)+1); #basis plus nuisance term m()
                                },
                                trainbasis = function(Z,n_knots) {
                                  m <- 4 # with intercept

                                  #internal knots (within unit circle)
                                  IKnots <- quantile(Z,probs = seq(n_knots)/(n_knots+1))
                                  B <<- splines2::bSpline(Z,
                                                          knots = IKnots,
                                                          Boundary.knots = c(-1,1), #unit circle is boundary
                                                          degree = m-1, #cubic
                                                          intercept = FALSE)
                                  dB <<- splines2::dbs(Z,
                                                       derivs = 1L,
                                                       knots = IKnots,
                                                       Boundary.knots = c(-1,1), #unit circle is boundary
                                                       degree = m-1, #cubic
                                                       intercept = FALSE)
                                  B
                                },
                                getbasis = function(Znew){
                                  #get new basis matrices
                                  Bnew <- splines2::predict.bSpline2(B,ewx=Znew)
                                  dBnew <- splines2::predict.dbs(dB,newx=Znew)

                                  list(B = Bnew, dB = dBnew)
                                }
                              )
)
