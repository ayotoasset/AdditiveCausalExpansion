## R6 object #####
KernelClass_SE_R6 <- R6::R6Class("SqExpKernel",
                             cloneable = FALSE,
                             class = FALSE,
                             portable = FALSE,
                             public = list(
                                 parameters = NULL,
                                 invKmatn = NULL,
                                 Kmat = NULL,
                                 Karray = NULL,
                                 B = NULL,
                                 p = NULL,
                                 stdy = 1,
                                 initialize = function(p_arg, B_arg, ext_init_parameters, std_y_arg=1, verbose=FALSE) {
                                   if (verbose) cat("Using SE kernel\n")
                                   B <<- B_arg
                                   p <<- p_arg
                                   parameters <<- ext_init_parameters
                                   stdy <<- std_y_arg# hand over moment for correct RMSE
                                 },
                                 kernel_mat = function(X1, X2, Z1, Z2) {
                                   #intended use for prediction
                                   Klist <- kernmat_SE_cpp(X1, X2, Z1, Z2, parameters) #lambda and L

                                 },
                                 kernel_mat_sym = function(X, Z) {
                                   #intended use for the gradient step
                                   Klist <- kernmat_SE_symmetric_cpp(X, Z, parameters)
                                   Kmat <<- Klist$full
                                   Karray <<- Klist$elements
                                   invisible(Klist)
                                 },
                                 getinv_kernel = function(X, Z) {
                                   # get matrices and return inverse for prediction
                                   kernel_mat_sym(X, Z)
                                   invKmatList <- invkernel_cpp(Kmat, c(parameters[1])) #no error handling
                                   invKmatn <<- invKmatList$inv
                                   invisible(invKmatList)
                                 },
                                 init_variance = function(init_iter, y, X, Z, Optim, verbose=TRUE, lr_multiplier=10) {
                                   # TODO: make efficient by calculating only the variance gradient instead of all gradients (in grad_SE_cpp)

                                   stats <- get_train_stats(y, X, Z)
                                   if (verbose) cat(sprintf("PRE |  0 | log Evidence %9.4f | RMSE %9.4f | Norm. noise var: %3.4f\n",
                                               stats[2], stats[1], exp(parameters[1])))

                                   Optim$lr <- Optim$lr * lr_multiplier
                                   for (iter in 1:init_iter) {
                                     eigenval <- getinv_kernel(X, Z)$eigenval;
                                     gradients <- grad_SE_cpp(y, X, Z,
                                                              Kmat, Karray,
                                                              invKmatn, eigenval,
                                                              parameters, stats, B, stdy)

                                     gradients[2:length(gradients)] <- 0
                                     parameters <<- Optim$update(iter, parameters, gradients)
                                     private$mean_solution(y)
                                   }
                                   #if ((class(Optim)=="AdamOpt") | (class(Optim)=="NadamOpt")) {
                                   #  Optim$m[0] <- 0
                                   #  Optim$v[0] <- 0
                                   #}

                                   if (verbose) cat(sprintf("PRE | %d | log Evidence %9.4f | RMSE %9.4f | Norm. noise var: %3.4f\n",
                                               init_iter, stats[2], stats[1], exp(parameters[1])))
                                   Optim$lr <- Optim$lr / lr_multiplier

                                 },
                                 para_update = function(iter, y, X, Z, Optim, printevery=100, verbose=TRUE) {
                                   #update Kmat and invKmat in the class environment
                                   stats <- c(0, 0)
                                   eigenval <- getinv_kernel(X, Z)$eigenval;

                                   if (iter==1) private$mean_solution(y) #overwrites mu gradient update

                                   gradients <- grad_SE_cpp(y, X, Z,
                                                            Kmat, Karray,
                                                            invKmatn, eigenval,
                                                            parameters, stats, B, stdy)
                                   #print(gradients)

                                   parameters <<- Optim$update(iter, parameters, gradients)

                                   private$mean_solution(y)

                                   if ((iter %% printevery == 0) && verbose) {
                                     cat(sprintf("%5d | log Evidence %9.4f | RMSE %9.4f | Norm. noise var: %3.4f | Gradient L2: %3.4f\n",
                                                 iter, stats[2], stats[1], exp(parameters[1]), norm(gradients)))
                                   }

                                   invisible(stats)
                                 },
                                 get_train_stats = function(y, X, Z, invKmatList){
                                   if (missing(invKmatList)) {
                                     #do not update the inverse
                                     Klist = kernel_mat_sym(X, Z)
                                     invKmatList <- invkernel_cpp(Klist$full, c(parameters[1]))
                                   }

                                   #private$mean_solution(y)
                                   stats <- stats_cpp(y, Kmat, invKmatList$inv,
                                                      invKmatList$eigenval, c(parameters[2]), stdy)

                                 },
                                 predict = function(y, X, Z, X2, Z2, mean_y, std_y) {
                                   n2 <- nrow(X2)

                                   K_xX <- kernmat_SE_cpp(X2,X,Z2,Z,parameters)$full
                                   K_xx <- kernmat_SE_symmetric_cpp(X2,Z2,parameters)$full

                                   outlist <- pred_cpp(y, parameters[1], parameters[2],
                                                       invKmatn, K_xX, K_xx, mean_y, std_y)
                                 },
                                 predict_marginal = function(y, X, Z, X2, Z2, dZ2,
                                                             mean_y, std_y, std_Z,
                                                             calculate_ate) {
                                   n <- length(y)
                                   n2 <- nrow(X2)

                                   Kmarginal_xX <- kernmat_SE_cpp(X2, X, dZ2, Z, parameters)$elements
                                   Kmarginal_xx <- kernmat_SE_symmetric_cpp(X2, dZ2, parameters)$elements

                                   outlist <- pred_marginal_cpp(y, Z2, parameters[1], parameters[2],
                                                                invKmatn,
                                                                Kmarginal_xX, Kmarginal_xx,
                                                                mean_y, std_y, std_Z, calculate_ate)
                                 }), # end public
                             private = list(
                                   mean_solution = function(y) {
                                     #using analytic solution
                                     parameters[2] <<- mu_solution_cpp(y, invKmatn)
                                   })
                             ) # end R6 class
