ns_spline <- R6::R6Class("ns_spline",
                         cloneable = FALSE,
                         class = FALSE,
                         portable = FALSE,
                         public = list(
                           B = NULL,
                           dB = NULL,
                           myknots = NULL,
                           dim = function(){
                             (ncol(B) + 1) #basis plus nuisance term m()
                           },
                           trainbasis = function(Z, n_knots, verbose=FALSE) {
                             if (verbose) cat("Using natural cubic-spline\n")
                             #internal knots (within unit circle)
                             if (n_knots > 0) IKnots <- quantile(Z, probs = seq(n_knots) / (n_knots + 1))
                             else IKnots <- NULL
                             Boundary.knots <- c(-1, 1)
                             myknots <<- c(IKnots, Boundary.knots)
                             B <<- ncs_basis(Z, myknots)
                             #not really needed:
                             dB <<- ncs_basis_deriv(Z, myknots)
                             invisible(B)
                           },
                           testbasis = function(Znew){
                             #get new basis matrices
                             if(!missing(Znew)){
                               Bnew <- ncs_basis(Znew, myknots)
                               dBnew <- ncs_basis_deriv(Znew, myknots)
                             } else {
                               #return stored matrices
                               Bnew <- B
                               dBnew <- dB
                             }
                             invisible(list(B = Bnew, dB = dBnew))
                           }
                         )
)
