binary_spline <- setRefClass("binary_spline",
                         fields = list(B = "matrix",   #design matrix
                                       dB = "matrix"), #derivative off design matrix
                         methods = list(
                           trainbasis = function(Z) {
                             m <<- 4 # with intercept
                             n <- length(Z)

                             #generate matrices without intercept!
                             B <<- as.matrix(Z,n,1)
                             dB <<- matrix(1,n,1) # 0 -> 1

                             B
                           },
                           testbasis = function(Znew){
                             n <- length(Znew)

                             Bnew <- Znew
                             dBnew <- matrix(1,n,1) # 0 -> 1

                             list(B = Bnew, dB = dBnew)
                           }


                         )
)
