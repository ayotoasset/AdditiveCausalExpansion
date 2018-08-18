## define optimizer classes - R6 objects #####
optAdam <- R6::R6Class("AdamOpt",
                       cloneable = FALSE,
                       class = FALSE,
                       portable = FALSE,
                       public = list(
                           m = NULL, v = NULL,
                           lr = 0.0001,
                           beta1 = 0.9,
                           beta2 = 0.999,
                           norm.clip = FALSE,
                           clip.at = 1,
                           initialize = function(KernelObj, lr, beta1, beta2, norm.clip, clip.at){
                             m <<- KernelObj$parameters*0
                             v <<- KernelObj$parameters*0
                             lr <<- lr
                             beta1 <<- beta1
                             beta2 <<- beta2
                             norm.clip <<- norm.clip
                             clip.at <<- clip.at

                             if (beta2 <= 0.02) {
                               cat("Adam:Beta2 is very small. This can lead to numeric instability.\n")
                             }
                           },
                           update = function(iter, parameters, gradients) {

                             norm_clip_cpp(norm.clip, gradients, clip.at)

                             if (Adam_cpp(iter, lr, beta1, beta2, 1e-8,
                                          m, v, gradients, parameters) == FALSE) {
                               stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                             }
                             parameters
                           }) # end public
                       ) # end Adam class

optNadam <- R6::R6Class("NadamOpt",
                        cloneable = FALSE,
                        class = FALSE,
                        portable = FALSE,
                        public = list(
                          m = NULL, v = NULL,
                          lr = 0.0001,
                          beta1 = 0.9,
                          beta2 = 0.999,
                          norm.clip = FALSE,
                          clip.at = 1,
                          initialize = function(KernelObj, lr, beta1, beta2, norm.clip, clip.at){
                            m <<- KernelObj$parameters*0
                            v <<- KernelObj$parameters*0
                            lr <<- lr
                            beta1 <<- beta1
                            beta2 <<- beta2
                            norm.clip <<- norm.clip
                            clip.at <<- clip.at

                            if (beta2 <= 0.02) {
                              cat("Nadam:Beta2 is very small. This can lead to numeric instability.\n")
                            }
                          },
                          update = function(iter, parameters, gradients) {

                            norm_clip_cpp(norm.clip, gradients, clip.at)

                            if (Nadam_cpp(iter, lr, beta1, beta2, 1e-8,
                                          m, v, gradients, parameters) == FALSE) {
                              stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                            }
                            parameters
                          }) # end public
) # end Nadam class

optNesterov <- R6::R6Class("NesterovOpt",
                           cloneable = FALSE,
                           class = FALSE,
                           portable = FALSE,
                           public = list(
                             nu = NULL,
                             lr = 0.0001,
                             momentum = 0.0,
                             norm.clip = TRUE,
                             clip.at = 1,
                             initialize = function(KernelObj, lr,  momentum, norm.clip, clip.at){
                               nu <<- KernelObj$parameters*0;
                               lr <<- lr
                               momentum <<- momentum
                               norm.clip <<- norm.clip
                               clip.at <<- clip.at
                             },
                             update = function(iter, parameters, gradients) {

                               norm_clip_cpp(norm.clip, gradients, clip.at)

                               if (Nesterov_cpp(lr, momentum, nu, gradients, parameters) == FALSE) {
                                 stop("Some gradients are not finite, NaN, or NA. Often this is due to too large learning rates.")
                               }
                               parameters
                             }) # end public
) # end NAG class

