## define optimizer classes
optAdam <- setRefClass("AdamOpt",
                       fields = list(m = "vector",
                                     v = "vector",
                                     lr = "numeric",
                                     beta1 = "numeric",
                                     beta2 = "numeric",
                                     norm.clip = "logical",
                                     clip.at = "numeric"),
                       methods = list(
                         update = function(iter, parameters, gradients) {

                           if (norm.clip) {
                             gradients <- norm.clip(gradients, max.length = clip.at, type = "2")
                           }

                           if (Adam_cpp(iter, lr, beta1, beta2, 1e-8,
                                        m, v, gradients, parameters) == FALSE) {
                             stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                           }
                           parameters
                         },
                         initOpt = function(KernelObj){
                           m <<- KernelObj$parameters*0
                           v <<- KernelObj$parameters*0
                         }
                       )
)
optNadam <- setRefClass("NadamOpt",
                        fields = list(m = "vector",
                                      v = "vector",
                                      lr = "numeric",
                                      beta1 = "numeric",
                                      beta2 = "numeric",
                                      norm.clip = "logical",
                                      clip.at = "numeric"),
                        methods = list(
                          update = function(iter, parameters, gradients) {

                            if (norm.clip) {
                              gradients <- norm.clip(gradients, max.length = clip.at, type = "2")
                            }

                            if (Nadam_cpp(iter, lr, beta1, beta2, 1e-8,
                                          m, v, gradients, parameters) == FALSE) {
                              stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                            }
                            parameters
                          },
                          initOpt = function(KernelObj){
                            m <<- KernelObj$parameters*0
                            v <<- KernelObj$parameters*0
                          }
                        )
)

optNesterov <- setRefClass("NesterovOpt",
                           fields = list(nu = "vector",
                                         lr = "numeric",
                                         momentum = "numeric",
                                         norm.clip = "logical",
                                         clip.at = "numeric"),
                           methods = list(
                             update = function(iter, parameters, gradients) {

                               if (norm.clip) {
                                 gradients <- norm.clip(gradients, max.length = clip.at, type = "2")
                               }

                               if(Nesterov_cpp(lr, momentum, nu, gradients, parameters) == FALSE) {
                                 stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                               }
                               parameters
                             },
                             initOpt = function(KernelObj){
                               nu <<- KernelObj$parameters*0;
                             }
                           )
)


norm.clip <- function(gradients, max.length = 1, type = "2") {
  # normalizes gradient if too large
  # arg in:
  # - gradients: vector of the gradients
  # - max.length: normalize above this value
  # - type: type of norm as in base::norm

  mynorm <- base::norm(gradients, type = type)
  if (mynorm > max.length) {
    gradients <- gradients / mynorm
  }

  gradients

}
