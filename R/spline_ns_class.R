ns_spline <- setRefClass("ns_spline",
                        fields = list(myknots = "vector", #store knot positions
                                      B = "matrix",   #design matrix
                                      dB = "matrix"), #derivative of design matrix
                        methods = list(
                          dim = function(){
                            (ncol(B)+1); #basis plus nuisance term m()
                          },
                          trainbasis = function(Z,n_knots) {
                            #internal knots (within unit circle)
                            IKnots <- quantile(Z,probs = seq(n_knots)/(n_knots+1))
                            Boundary.knots <- c(-1,1)
                            myknots <<- c(IKnots,Boundary.knots);

                            B <<- ncs_basis(Z,myknots)
                            #not really needed:
                            dB <<- ncs_basis_deriv(Z,myknots)

                            B
                          },
                          getbasis = function(Znew){
                            #get new basis matrices
                            if(!missing(Znew)){
                              Bnew <- ncs_basis(Znew,myknots)
                              dBnew <- ncs_basis_deriv(Znew,myknots)
                            } else {
                              #return stored matrices
                              Bnew <- B
                              dBnew <- dB
                            }

                            list(B = Bnew, dB = dBnew)
                          }
                        )
)
