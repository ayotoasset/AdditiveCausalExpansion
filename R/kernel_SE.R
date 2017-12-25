KernelClass_SE <- setRefClass("SqExpKernel",
                                 fields = list(parameters = "matrix",
                                               invKmatn = "matrix", #training inverse of the Kernel matrix with noise
                                               Kmat = "matrix", #kernel matrix
                                               Karray = "array", #training Kernel matrices
                                               B = "numeric",
                                               p = "numeric"), #basis dimension
                                 methods = list(
                                   parainit = function(y,p,B,Z) {
                                     Z <- as.matrix(Z)
                                     B <<- B
                                     p <<- p
                                     n <- length(y)
                                     parameters <<- matrix(c(0, #sigma
                                                      0, #mu
                                                      -log(c(1, diag(t(Z) %*% Z) / n)),#(c(0.3,0.3)), #lambda as B(Z) B(Z) becomes very small
                                                      rep(c(-0.1, rep(-0.1, B - 1)), p)  #L: for nuisance term -0.1 for others 0.1
                                                      ))
                                     #old:
                                     #parameters <<- list(sigma=log(var(y)),
                                                         #sigma_z=log(runif(1,min=1,max=2)),
                                                         #lambdam=log(runif(1,min=0.2,max=1)),
                                                         #lambdaa=log(runif(1,min=0.2,max=1)),
                                                         #Lm=rep(-0.1,p),La=rep(0.1,p),
                                                         #mu = mean(y))

                                     # [ 1 | 1 | B | B | B | B ]
                                     # B is repeated p+1 times with the first block being lambda[1:B]

                                     #selectors: C++               | R
                                     #- sigma:   [0,0]             | [1]
                                     #- mu:      [1,1]             | [2]
                                     #- lambda:  [2,1+B]           | [3:(2+B)]
                                     #- L:       [2+B,1+B*(p+1)] | [(3+B):(2+B*(p+1))]
                                     #in Cpp:
                                     # L(i,b) = parameters[2+(B-1)+b+B*i] = parameters[1+B+b+B*i] = parameters[1+b+B*(i+1)]
                                   },
                                   kernel_mat = function(X1, X2, Z1, Z2) {
                                     #intended use for prediction
                                     Klist <- kernmat_SE_cpp(X1, X2, Z1, Z2, parameters) #lambda and L

                                   },
                                   kernel_mat_sym = function(X, Z) {
                                     #intended use for the gradient step
                                     Klist <- kernmat_SE_symmetric_cpp(X, Z, parameters)
                                     #Klist <- kernmat_SE_cpp(X,X,Z,Z,parameters)
                                     Kmat <<- Klist$full
                                     Karray <<- Klist$elements
                                     Klist
                                   },
                                   getinv_kernel = function(X, Z, noeigen=FALSE) {
                                     #get matrices and return inverse for prediction

                                     #Klist <- kernel_mat(X,X,Z,Z)
                                     #Kmat <<- Klist$full
                                     #Karray <<- Klist$elements

                                     kernel_mat_sym(X, Z)
                                     if(noeigen){
                                       invKmatList <- invkernel_no_eigen_cpp(Kmat, c(parameters[1]) )
                                     } else {
                                       invKmatList <- invkernel_cpp(Kmat, c(parameters[1])) #no error handling, return eigenvalues!
                                     }
                                     invKmatn <<- invKmatList$inv
                                     invKmatList
                                   },
                                   para_update = function(iter, y, X, Z, Optim, printevery=100) {
                                     #update Kmat and invKmat in the class environment
                                     stats <- c(0,0)
                                     invKmatList <- getinv_kernel(X, Z);

                                     gradients <- grad_SE_cpp(y, X, Z,
                                                              Kmat, Karray,
                                                              invKmatn, invKmatList$eigenval,
                                                              parameters, stats, B)
                                     parameters <<- Optim$update(iter, parameters, gradients)

                                     mean_solution(y) #overwrites mu gradient update

                                     if (iter%%printevery == 0) {
                                       cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f \n", iter, stats[2], stats[1]))
                                      }

                                     stats
                                   },
                                   get_train_stats = function(y, X, Z, invKmatList){
                                     if (missing(invKmatList)) {
                                       #do not update the inverse
                                       Klist = kernel_mat_sym(X, Z)
                                       invKmatList <- invkernel_cpp(Klist$full, c(parameters[1]))
                                     }

                                     stats <- stats_cpp(y, Kmat, invKmatList$inv,
                                                        invKmatList$eigenval, parameters[2])

                                   },
                                   mean_solution = function(y) {
                                     #using analytic solution
                                     parameters[2] <<- mu_solution_cpp(y, invKmatn)
                                   },
                                   predict = function(y,X,Z,X2,Z2,mean_y,std_y){
                                     n2 <- nrow(X2)

                                     K_xX <- kernmat_SE_cpp(X2,X,Z2,Z,parameters)$full
                                     K_xx <- kernmat_SE_symmetric_cpp(X2,Z2,parameters)$full

                                     outlist <- pred_cpp(y, parameters[1], parameters[2],
                                                         invKmatn, K_xX, K_xx, mean_y, std_y)
                                   },
                                   predict_marginal = function(y, X, Z, X2, dZ2,
                                                               mean_y, std_y, std_Z,
                                                               calculate_ate) {
                                     n <- length(y)
                                     n2 <- nrow(X2)

                                     Kmarginal_xX <- kernmat_SE_cpp(X2, X, dZ2, Z, parameters)$elements
                                     Kmarginal_xx <- kernmat_SE_symmetric_cpp(X2, dZ2, parameters)$elements

                                     outlist <- pred_marginal_cpp(y, parameters[1], parameters[2],
                                                                  invKmatn,
                                                                  Kmarginal_xX, Kmarginal_xx,
                                                                  mean_y, std_y, std_Z, calculate_ate)
                                   }
                                 )
)
