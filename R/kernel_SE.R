KernelClass_SE <- setRefClass("SqExpKernel",
                                 fields = list(parameters = "matrix",
                                               invKmatn = "matrix", #training inverse of the Kernel matrix with noise
                                               Kmat = "matrix", #kernel matrix
                                               Karray = "array", #training Kernel matrices
                                               B = "numeric",
                                               p = "numeric"), #basis dimension
                                 methods = list(
                                   parainit = function(y,p,B) {
                                     B <<- B
                                     p <<- p
                                     parameters <<- matrix(c(log(1), #sigma
                                                      mean(y), #mu
                                                      log(seq(B,1)),
                                                      rep(0.1,(p) * B) # lambda and L
                                                      ))
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
                                   kernel_mat = function(X1,X2,Z1,Z2) {
                                     #intended use for prediction
                                     Klist <- kernmat_SE_cpp(X1,X2,Z1,Z2,parameters) #lambda and L

                                   },
                                   kernel_mat_sym = function(X,Z) {
                                     #intended use for the gradient step
                                     Klist <- kernmat_SE_symmetric_cpp(X,Z,parameters)
                                     Kmat <<- Klist$full
                                     Karray <<- Klist$elements
                                     Klist
                                   },
                                   getinv_kernel = function(X,Z,noeigen=FALSE) {
                                     #get matrices and return inverse for prediction
                                     kernel_mat_sym(X,Z)
                                     if(noeigen){
                                       invKmatList <- invkernel_no_eigen_cpp(Kmat,c(parameters[1]))
                                     } else {
                                       invKmatList <- invkernel_cpp(Kmat,c(parameters[1])) #no error handling, return eigenvalues!
                                     }
                                     invKmatn <<- invKmatList$inv
                                     invKmatList
                                   },
                                   para_update = function(iter,y,X,Z,Optim) {
                                     #update Kmat and invKmat in the class environment
                                     stats <- c(0,0)
                                     invKmatList <- getinv_kernel(X,Z);
                                     gradients <- grad_SE_cpp(y,X,Z,Kmat,Karray,invKmatn,invKmatList$eigenval,
                                                              parameters,
                                                              stats,B)

                                     parameters <<- Optim$update(iter,parameters,gradients)
                                     mean_solution(y) #overwrites mu gradient update

                                     if(iter%%100 == 0){ cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f \n", iter, stats[2], stats[1])) }

                                     stats
                                   },
                                   get_train_stats = function(y,X,Z,invKmatList){
                                     if(missing(invKmatList)){
                                       #do not update the inverse
                                       Klist = kernel_mat_sym(X,Z)
                                       invKmatList <- invkernel_cpp(Klist$full,c(parameters[1]))
                                     }

                                     stats <- stats_SE_cpp(y,Kmat, invKmatList$inv,invKmatList$eigenval, parameters[2])

                                   },
                                   mean_solution = function(y){
                                     #using analytic solution
                                     parameters[2] <<- mu_solution_cpp(y, invKmatn)
                                   },
                                   predict = function(y,X,Z,X2,Z2,mean_y,var_y){
                                     n2 <- nrow(X2)

                                     K_xX <- kernel_mat(X2,X,Z2,Z)$full
                                     K_xx <- kernel_mat_sym(X2,Z2)$full

                                     outlist <- pred_cpp(y,parameters[1],parameters[2],invKmatn, K_xX, K_xx, mean_y,var_y)
                                   },
                                   predict_marginal = function(y,X,Z,X2,dZ2,mean_y,std_y,std_Z,calculate_ate){

                                     n <- length(y);
                                     n2 <- nrow(X2);

                                     Kmarginal_xX <- kernel_mat(X2,X,dZ2,Z)$elements
                                     Kmarginal_xx <- kernel_mat_sym(X2,dZ2)$elements

                                     outlist <- pred_marginal_cpp(y,parameters[1],parameters[2],
                                                                  invKmatn,
                                                                  Kmarginal_xX, Kmarginal_xx,
                                                                  mean_y,std_y,std_Z,calculate_ate)
                                   }
                                 )
)
