linear_spline <- R6::R6Class("linear_spline",
                             cloneable = FALSE,
                             class = FALSE,
                             portable = FALSE,
                             public = list(
                               B = NULL,
                               dB = NULL,
                               dim = function(){
                                 (ncol(B) + 1) #basis plus nuisance term m()
                               },
                               trainbasis = function(Z, n_knots, verbose=FALSE) {
                                 if (verbose) cat("Using binary/linear-basis\n")
                                 n <- length(Z)
                                 #generate matrices without intercept!
                                 B <<- matrix(Z,n,1)
                                 dB <<- matrix(1,n,1) # 0 -> 1
                                 invisible(B)
                               },
                               testbasis = function(Znew){
                                 #get new basis matrices
                                 if(!missing(Znew)){
                                   Bnew <- matrix(Znew)
                                   dBnew <- matrix(1, length(Znew), 1) # 0 -> 1
                                 } else {
                                   #return stored matrices
                                   Bnew <- matrix(B)
                                   dBnew <- matrix(dB)
                                 }
                                 invisible(list(B = Bnew, dB = dBnew))
                               }
                             )
)

