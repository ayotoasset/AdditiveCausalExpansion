square_spline <- setRefClass("square_spline",
                             fields = list(B = "matrix",   #design matrix
                                           dB = "matrix"), #derivative off design matrix
                             methods = list(
                               dim = function(){
                                 (ncol(B)+1); #basis plus nuisance term m()
                               },
                               trainbasis = function(Z,n_knots) {
                                 n <- length(Z)

                                 #generate matrices without intercept!
                                 B <<- matrix(c(Z,Z^2),n,2)
                                 dB <<- matrix(c(rep(1,n),Z) ,n,2) # 0 -> 1

                                 B
                               },
                               testbasis = function(Znew){
                                 #get new basis matrices
                                 if(!missing(Znew)){
                                   n <- length(Znew)
                                   Bnew <- matrix(c(Znew,Znew^2),n,2)
                                   dBnew <- matrix(c(rep(1,n),Znew) ,n,2)
                                 } else {
                                   #return stored matrices
                                   Bnew <- matrix(B)
                                   dBnew <- matrix(dB)
                                 }
                                 list(B = Bnew, dB = dBnew)
                               }
                             )
)
