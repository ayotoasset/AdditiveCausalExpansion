KernelClass_SE_nonadd <- setRefClass("SqExpNonAddKernel",
                                  fields = list(parameters = "matrix",
                                                invKmatn = "matrix", #training inverse of the Kernel matrix with noise
                                                Kmat = "matrix", #kernel matrix
                                                Karray = "array", #training Kernel matrices
                                                B = "numeric",
                                                p = "numeric"), #basis dimension
                                  methods = list(
                                    parainit = function(y,p,voidB) {
                                      #if(B!=2){cat("Isotropic Kernel only allowed for uni-variate basis (linear spline)\n")}
                                      B <<- 1 #only for univariate Z
                                      p <<- p+1 #attach Z to X
                                      parameters <<- matrix(c(log(1), #sigma
                                                              0, #mu
                                                              log(seq(1,1)),
                                                              rep(0.1,(p) * 1) # lambda and L
                                      ))
                                    },
                                    kernel_mat = function(X1,X2,Z1,Z2) {
                                      #intended use for prediction

                                      Klist <- kernmat_SE_noadd_cpp(cbind(X1,Z1),cbind(X2,Z2),parameters) #lambda and L

                                    },
                                    kernel_mat_sym = function(X,Z) {
                                      #intended use for the gradient step
                                      Klist <- kernmat_SE_symmetric_noadd_cpp(cbind(X,Z),parameters)
                                      Kmat <<- Klist$full
                                      Karray <<- Klist$elements
                                      Klist
                                    },
                                    getinv_kernel = function(X,Z) {
                                      #get matrices and return inverse for prediction
                                      kernel_mat_sym(X,Z)
                                      invKmatList <- invkernel_cpp(Kmat,c(parameters[1])) #no error handling, return eigenvalues!
                                      invKmatn <<- invKmatList$inv
                                      invKmatList
                                    },
                                    para_update = function(iter,y,X,Z,Optim) {
                                      #update Kmat and invKmat in the class environment
                                      stats <- c(0,0)
                                      invKmatList <- getinv_kernel(X,Z);
                                      gradients <- grad_SE_cpp(y,X,Z,Kmat,Karray,invKmatn,invKmatList$eigenval,
                                                               parameters,
                                                               stats,1)

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

                                      stats <- stats_cpp(y,Kmat, invKmatList$inv,invKmatList$eigenval, parameters[2])

                                    },
                                    mean_solution = function(y){
                                      #using analytic solution
                                      parameters[2] <<- mu_solution_cpp(y, invKmatn)
                                    },
                                    predict = function(y,X,Z,X2,Z2,mean_y,std_y){
                                      n2 = nrow(X2)

                                      K_xX <- kernel_mat(X2,X,Z2,Z)$full
                                      K_xx <- kernel_mat_sym(X2,Z2)$full

                                      outlist <- pred_cpp(y,parameters[1],parameters[2],invKmatn, K_xX, K_xx, mean_y,std_y)
                                    },
                                    predict_marginal = function(y,X,Z,X2,Z2,mean_y,std_y,calculate_ate){

                                      n = length(y);
                                      n2 = nrow(X2);

                                      Kmarginal_xX <- kernel_mat(X2,X,Z2,Z)$full
                                      Kmarginal_xX <- kernmat_marginal_noadd_cpp(Kmarginal_xX,Z2,Z,parameters)
                                      Kmarginal_xX <- array(Kmarginal_xX,dim=c(n2,n,1))

                                      Kmarginal_xx <- kernel_mat_sym(X2,Z2)$full
                                      Kmarginal_xx <- kernmat_marginal_noadd_cpp(Kmarginal_xx,Z2,Z2,parameters)
                                      Kmarginal_xx <- array(Kmarginal_xx,dim=c(n2,n2,1))


                                      outlist <- pred_marginal_cpp(y,parameters[1],parameters[2],
                                                                   invKmatn,
                                                                   Kmarginal_xX, Kmarginal_xx,
                                                                   mean_y,std_y,calculate_ate)
                                    }
                                  )
)
