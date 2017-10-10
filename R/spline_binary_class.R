binary_spline <- setRefClass("binary_spline",
                         fields = list(B = "matrix",   #design matrix
                                       dB = "matrix"), #derivative off design matrix
                         methods = list(
                           dim = function(){
                             (ncol(B)+1); #basis plus nuisance term m()
                           },
                           trainbasis = function(Z) {
                             n <- length(Z)

                             #generate matrices without intercept!
                             B <<- as.matrix(Z,n,1)
                             dB <<- matrix(1,n,1) # 0 -> 1

                             B
                           },
                           testbasis = function(Znew){
                             #get new basis matrices
                             if(!missing(Znew)){
                               n <- length(Znew)
                               Bnew <- Znew
                               dBnew <- matrix(1,n,1) # 0 -> 1
                             } else {
                               #return stored matrices
                               Bnew <- B
                               dBnew <- dB
                             }
                             list(B = Bnew, dB = dBnew)
                           }


                         )
)
