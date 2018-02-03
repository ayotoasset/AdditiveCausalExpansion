square_spline <- R6::R6Class("square_spline",
                             cloneable = FALSE,
                             class = FALSE,
                             portable = FALSE,
                             public = list(
                               B = NULL,
                               dB = NULL,
                               dim = function(){
                                 (ncol(B) + 1) #basis plus nuisance term m()
                               },
                               trainbasis = function(Z, n_knots) {
                                 cat("Using square-basis\n")
                                 n <- length(Z)

                                 #generate matrices without intercept!
                                 B <<- matrix(c(Z, Z ^ 2), n, 2)
                                 dB <<- matrix(c(rep(1, n), 2 * Z), n, 2)
                                 invisible(B)
                               },
                               testbasis = function(Znew){
                                 #get new basis matrices
                                 if(!missing(Znew)){
                                   n <- length(Znew)
                                   Bnew <- matrix(c(Znew, Znew ^ 2), n, 2)
                                   dBnew <- matrix(c(rep(1, n), 2 * Znew), n, 2)
                                 } else {
                                   #return stored matrices
                                   Bnew <- matrix(B)
                                   dBnew <- matrix(dB)
                                 }
                                 invisible(list(B = Bnew, dB = dBnew))
                               }
                             )
)

