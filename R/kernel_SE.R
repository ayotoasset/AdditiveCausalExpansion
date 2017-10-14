KernelClass_SE <- setRefClass("SqExpKernel",
                                 fields = list(parameters = "list",
                                               invKmatn = "matrix", #training inverse of the Kernel matrix with noise
                                               Kmat = "matrix", #kernel matrix
                                               Karray = "array", #training Kernel matrices
                                               K = "list", #list of kernels
                                               B = "numeric", #basis dimension
                                               px = "numeric"),
                                 methods = list(
                                   parainit = function(y,B) {

                                     #for momentum and adam
                                     parameters <<- list(sigma=log(var(y)),
                                                         lambda=log(runif(B,min=0.2,max=1)), #vector of lambdas
                                                         L = matrix(0.1,p,B), #list of vectors with ARD weights for SE kernel
                                                         mu = mean(y)
                                     )
                                   },
                                   kernel_mat = function(X1,X2,Z1,Z2) {
                                     #intended use for prediction
                                     X1 = as.matrix(X1);
                                     X2 = as.matrix(X2);
                                     Z1 = as.matrix(Z1);
                                     Z2 = as.matrix(Z2);

                                     Klist = kernmat_SE_cpp(X1,X2,Z1,Z2,parameters$lambda,parameters$L)
                                     #Kmat <<- Klist$full
                                     #Karray <<- Klist$elements

                                   },
                                   kernel_mat_sym = function(X1,Z2) {
                                     #intended use for the gardaient step
                                     X = as.matrix(X1);
                                     Z = as.matrix(Z2);

                                     Klist = kernmat_SE_sym_cpp(X,Z,parameters$lambda,parameters$L)
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

                                     gradients <- grad_SE_cpp(y,as.matrix(X),as.matrix(Z),Kmat,Karray,invKmatList$inv,invKmatList$eigenval,
                                                              parameters$sigma,parameters$mu,parameters$lambda,parameters$L,stats)

                                     parameters <<- Optim$update(iter,parameters,gradients)
                                     mean_solution(y,z) #overwrites mu gradient update

                                     if(iter%%100 == 0){ cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f \n", iter, stats[2], stats[1])) }
                                     iter
                                   },
                                   get_train_stats = function(y,X,z,invKmatList){
                                     if(missing(invKmatList)){
                                       invKmatList <- getinv_kernel(X,z)
                                     }

                                     stats <- stats_SE(y,Kmat, invKmatList$inv,invKmatList$eigenval, parameters$mu)

                                   },
                                   mean_solution = function(y,z){
                                     #using analytic solution
                                     parameters$mu <<- mu_solution_cpp(y, invKmatn)
                                     parameters <<- parameters #??????
                                   },
                                   predict = function(y,X,z,X2,z2){
                                     n2 = nrow(X2)

                                     K_xX = kernel_mat(X2,X,z2,z)$full
                                     K_xx = kernel_mat_sym(X2,z2)$full

                                     #map: put in cpp
                                     map = K_xX %*% invKmatn %*% (y - parameters$mu) + parameters$mu
                                     cov = K_xx - K_xX %*% invKmatn %*% t(K_xX) + exp(parameters$sigma) * diag(n2)
                                     uncentered_ci = cbind(-1.96*diag(cov),1.96 * diag(cov))
                                     list(map=map,ci=uncentered_ci)
                                   },
                                   predict_treat = function(y,X,z,X2){
                                     n = length(y);
                                     n2 = nrow(X2)
                                     z2 = rep(1,n2)

                                     K_xX = kernel_mat(X2,X,z2,z)$elements
                                     K_xx = kernel_mat_sym(X2,z2)$elements

                                     #map
                                     map = K_xX %*% invKmatn %*% (y - parameters$mu)
                                     ate = mean(map);

                                     #ci
                                     cov = K_xx - K_xX %*% invKmatn %*% t(K_xX)
                                     uncentered_ci = cbind(-1.96*diag(cov),1.96 * diag(cov))
                                     ate_cov = sum(sum(cov))/(n^2);
                                     ate_uncentered_ci = cbind( -1.96 * sqrt(ate_cov), 1.96 * sqrt(ate_cov) )

                                     list(map=map,ci=uncentered_ci,ate_map = ate , ate_ci =  ate_uncentered_ci)
                                   }
                                 )
)
