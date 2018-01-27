linear_spline <- setRefClass("linear_spline",
                         fields = list(B = "matrix",   #design matrix
                                       dB = "matrix"), #derivative off design matrix
                         methods = list(
                           dim = function(){
                             (ncol(B)+1); #basis plus nuisance term m()
                           },
                           trainbasis = function(Z,n_knots) {
                             n <- length(Z)

                             #generate matrices without intercept!
                             B <<- matrix(Z,n,1)
                             dB <<- matrix(1,n,1) # 0 -> 1

                             B
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
                             list(B = Bnew, dB = dBnew)
                           }
                         )
)
