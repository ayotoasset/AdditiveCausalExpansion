KernelClass_SE <- setRefClass("SqExpKernel",
                                 fields = list(parameters = "list",
                                               invKmatn = "matrix", #training inverse of the Kernel matrix with noise
                                               Kmat = "matrix", #kernel matrix
                                               Karray = "array", #training Kernel matrices
                                               B = "numeric",
                                               p = "numeric"), #basis dimension
                                 methods = list(
                                   parainit = function(y,p,B) {
                                     B <<- B
                                     p <<- p

                                     #for momentum and adam
                                     parameters <<- list(sigma=log(var(y)),
                                                         lambda=log(runif(B,min=0.2,max=1)), #vector of lambdas
                                                         L = matrix(0.1,p,B), #list of vectors with ARD weights for SE kernel
                                                         mu = mean(y)
                                     )
                                   },
                                   kernel_mat = function(X1,X2,Z1,Z2) {
                                     #intended use for prediction
                                     Klist <- kernmat_SE_cpp(X1,X2,Z1,Z2,parameters$lambda,matrix(parameters$L,p,B))

                                   },
                                   kernel_mat_sym = function(X,Z) {
                                     #intended use for the gradient step
                                     Klist = kernmat_SE_symmetric_cpp(X,Z,parameters$lambda,matrix(parameters$L,p,B))
                                     Kmat <<- Klist$full
                                     Karray <<- Klist$elements
                                   },
                                   getinv_kernel = function(X,Z) {
                                     #get matrices and return inverse for prediction

                                     #get new kernel and invert with noise
                                     Klist = kernel_mat_sym(X,Z)

                                     invKmatList <- invkernel_cpp(Kmat,parameters$sigma) #no error handling, return eigenvalues!
                                     invKmatn <<- invKmatList$inv

                                     invKmatList
                                   },
                                   para_update = function(iter,y,X,Z,Optim) {
                                     #update Kmat and invKmat in the class environment
                                     stats <- c(0,0)

                                     invKmatList <- getinv_kernel(X,Z);

                                     gradients <- grad_SE_cpp(y,X,Z,Kmat,Karray,invKmatn,invKmatList$eigenval,
                                                              parameters$sigma,parameters$lambda,
                                                              matrix(parameters$L,p,B), # otherwise issues with p=1
                                                              parameters$mu,
                                                              stats)

                                     parameters <<- Optim$update(iter,parameters,gradients)
                                     mean_solution(y) #overwrites mu gradient update

                                     if(iter%%100 == 0){ cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f \n", iter, stats[2], stats[1])) }
                                     iter
                                   },
                                   get_train_stats = function(y,X,z,invKmatList){
                                     if(missing(invKmatList)){
                                       invKmatList <- getinv_kernel(X,z)
                                     }

                                     stats <- stats_SE_cpp(y,Kmat, invKmatList$inv,invKmatList$eigenval, parameters$mu)

                                   },
                                   mean_solution = function(y){
                                     #using analytic solution
                                     parameters$mu <<- mu_solution_cpp(y, invKmatn)
                                     #parameters <<- parameters #??????
                                   },
                                   predict = function(y,X,Z,X2,Z2){
                                     n2 = nrow(X2)

                                     K_xX = kernel_mat(X2,X,Z2,Z)$full
                                     K_xx = kernel_mat_sym(X2,Z2)$full

                                     outlist <- pred_cpp(y, parameters$sigma,parameters$mu,invKmatn, K_xX, K_xx)
                                   },
                                   predict_treat = function(y,X,Z,X2,dZ2,isbinary){

                                     n = length(y);
                                     n2 = nrow(X2)

                                     Kmarginal_xX = kernel_mat(X2,X,dZ2,Z)$elements
                                     Kmarginal_xx = kernel_mat_sym(X2,dZ2)$elements

                                     outlist <- pred_cpp(y,
                                                         parameters$sigma,parameters$mu,
                                                         invKmatn,
                                                         Kmarginal_xX, Kmarginal_xx,
                                                         isbinary)
                                   }
                                 )
)
