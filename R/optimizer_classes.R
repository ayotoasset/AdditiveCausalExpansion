## define optimizer classes
optAdam <- setRefClass("AdamOpt",
                       fields = list(m = "vector",
                                     v = "vector",
                                     lr = "numeric",
                                     beta1 = "numeric",
                                     beta2 = "numeric",
                                     useHessian="logical"),
                       methods = list(
                         update = function(iter,parameters,gradients) {
                           if(Adam_cpp(iter,lr,beta1,beta2,1e-8,m,v,gradients,parameters)==FALSE) { stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.") }
                           parameters
                         },
                         initOpt = function(KernelObj){
                           m <<- KernelObj$parameters*0
                           v <<- KernelObj$parameters*0
                           useHessian <<- FALSE
                         })
)

optNadam <- setRefClass("NadamOpt",
                        fields = list(m = "vector",
                                      v = "vector",
                                      lr = "numeric",
                                      beta1 = "numeric",
                                      beta2 = "numeric",
                                      useHessian="logical"),
                        methods = list(
                          update = function(iter,parameters,gradients) {
                            if(Nadam_cpp(iter,lr,beta1,beta2,1e-8,m,v,gradients,parameters)==FALSE) { stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.") }
                            parameters
                          },
                          initOpt = function(KernelObj){
                            m <<- KernelObj$parameters*0
                            v <<- KernelObj$parameters*0
                            useHessian <<- FALSE
                          })
)

optNesterov <- setRefClass("NesterovOpt",
                           fields = list(nu = "vector",
                                         lr = "numeric",
                                         momentum = "numeric",
                                         useHessian="logical"),
                           methods = list(
                             update = function(iter,parameters,gradients) {
                               if(Nesterov_cpp(lr,momentum,nu,gradients,parameters)==FALSE) { stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.") }
                               parameters
                             },
                             initOpt = function(KernelObj){
                               nu <<- KernelObj$parameters*0; #same structure
                               useHessian <<- FALSE
                             })
)

optNewton <- setRefClass("NewtonOpt",
                         fields = list(nu = "vector",
                                       lr = "numeric",
                                       momentum = "numeric",
                                       inititer="numeric",
                                       useHessian="logical"),
                         methods = list(
                           update = function(iter,parameters,gradients,Hessian=diag(length(momentum)),B) {
                             if(iter>inititer){
                               if(Newton_cpp(iter,lr,momentum,nu,gradients,Hessian,parameters,B)==FALSE) { stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.") }
                             } else {
                               if(Nesterov_cpp(lr,momentum,nu,gradients,parameters)==FALSE) { stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.") }
                             }
                             parameters

                           },
                           initOpt = function(KernelObj,inititer = 5){
                             nu <<- KernelObj$parameters*0
                             inititer <<- inititer
                             useHessian <<- TRUE
                           })
)
