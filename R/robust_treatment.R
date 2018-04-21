#' Plot a sequence of average treatment effects with increasing posterior certainty
#'
#' @description {The method plots a sequence of average treatment effects by sequentially discarding points with the largest treatment uncertainty.
#' The shape of the curve indicates the robustness of the average treatment effect with respect to the treatment balance in the covariate space.
#' A substantial change of the treatment effect over the levels of removal is a good indicator for limited support but the reverse need not be true.
#'
#' Note that the average treatment effect on the treated/untreated (ATT and ATU) is only supported for training set evaluation (if `newX` missing).
#' }
#'
#'
#' @param object ace object
#' @param newX (optional) The covariate matrix for the treatment effect predictions. If omitted, uses training samples.
#' @param n.steps (default: 10) The number of steps to be evaluated. The steps are equally spaced from 0\% (1 observation) to 100\% (full dataset).
#'
#' @return A matrix with the index for each discarding point along the x-axis.
#'
#' @examples
#' library(ace)
#' ## Example with binary treatment similar to Hill (2011)'s
#'
#' set.seed(1231)
#' n <- 300
#'
#' # generate treatment
#' Z <- rbinom(n, 1, 0.3)
#'
#' # generate confounder and exogenous variable
#' X <- matrix(NaN, n, 1)
#' X[Z==1, ] <- rnorm(sum(Z), mean = 30,sd = 10)
#' X[Z==0, ] <- rnorm(n - sum(Z), mean = 20, sd = 10)
#' E <- runif(n) # exogenous variable
#' X <- data.frame(X, E)
#'
#' # sort Confounder for visualizations
#' sort.idx <- sort(X[, 1], index.return = TRUE)$ix
#'
#' # define and draw the reponse function
#' y_truefun <- function(x, z) {
#'     mat <- matrix(NaN, length(z), 1)
#'     mat[z==0, 1] <- matrix(72 + 3 * (x[z == 0,1] > 0) * sqrt(abs(x[z == 0, 1])), sum(z == 0), 1)
#'     mat[z==1, 1] <- matrix(90 + exp(0.06 * x[z == 1, 2]), sum(z == 1), 1)
#'     c(mat)}
#' y0_true <- y_truefun(X, rep(0, n))
#' y1_true <- y_truefun(X, rep(1, n))
#' Y0 <- rnorm(n, mean = y0_true, sd = 1)
#' Y1 <- rnorm(n, mean = y1_true, sd = 1)
#' Y <- Y0 * (1 - Z) + Y1 * Z
#'
#' # run model
#' my.ace <- ace.train(Y, X, Z,
#'                     kernel = "SE", optimizer = "Nadam",
#'                     learning_rate = 0.005, maxiter = 1000)
#'
#' # plot treatment curve
#' plot_ace(my.ace, 1, marginal = TRUE)
#'
#' # Check for robustness of ATE:
#' idx = robust_treatment(my.ace, n.steps=5)
#'
#' # Let's see which points were discarded (X[, 1] is the sole confounder here):
#' plot(X[, 1], idx[, 4], xlab="Retained Observations Flag", ylab="Confounder")
#' # When comparing to the treatment plot above, we keep the
#' # observations for which we have local overlap.
#'
#'
#' @export
robust_treatment <- function(object, newX, n.steps = 5) {
  if(!object$train_data$Zbinary) {
    stop("Average treatment effect requires bianry treatment indicator.")
  }

  insample = missing(newX)
  if (insample) {
    invisible(prediction <- predict.ace(object, marginal = TRUE))

    att_results <- matrix(NaN, n.steps+1, 4)
    atu_results <- matrix(NaN, n.steps+1, 4)
  } else {
    # out-of-sample
    invisible(prediction <- predict.ace(object, newX = newX, marginal = TRUE))
  }

  quantile.list <- seq(0, 1, 1 / (n.steps))
  discard.steps <- stats::quantile(prediction$var, probs = quantile.list)

  ate_results <- matrix(NaN, n.steps+1, 5 + 2 * insample) #includes other baseline results

  # define return matrix with selected index
  if (missing(newX)) n.pred <- nrow(object$train_data$X)
  else n.pred <- nrow(newX)
  idx = matrix(FALSE, n.pred, n.steps+1)

  for (i in 1:(n.steps+1)) {
    idx[, i] <- c(prediction$var <= discard.steps[i])
    ate_results[i, 5] <- sum(idx[, i])

    if(missing(newX)) {
      # within-of-sample
      invisible(predace_robust <- predict.ace(object,
                                              newX = matrix(object$train_data$X[idx[, i], ], sum(idx[, i]), ncol(object$train_data$X)),
                                              newZ = object$train_data$Z[c(idx[, i])],
                                              marginal=TRUE, return_average_treatments=TRUE, normalize=FALSE))
    } else {
      # out-of-sample
      predace_robust <- predict.ace(object, newX = newX[idx[, i], ], newZ = 1, marginal=TRUE, return_average_treatments=TRUE)
    }
    ate_results[i, 1:4] <- unlist(predace_robust$ate)
    if(insample) {
      #remaining treated:
      ate_results[i, 6] <- t(idx[, i]) %*% object$train_data$Z
      #remaining untreated:
      ate_results[i, 7] <- t(idx[, i]) %*% (1 - object$train_data$Z)
      #ATT
      if(ate_results[i, 6] > 0) {
        att_results[i, 1:4] <- unlist(predace_robust$att)
      } else {
        att_results[i, 1:4] <- rep(NA, 4)
      }

      #ATU
      if(ate_results[i, 7] > 0) {
        atu_results[i, 1:4] <- unlist(predace_robust$atu)
      } else {
        atu_results[i, 1:4] <- rep(NA, 4)
      }

    }
  }


  #For the colors see http://jfly.iam.u-tokyo.ac.jp/color/#what
  Blue100         <- grDevices::rgb(0,   114, 178, max = 255, alpha = (100 -  0) * 255 / 100)
  Blue50          <- grDevices::rgb(0,   114, 178, max = 255, alpha = (100 - 50) * 255 / 100)
  Vermillion100   <- grDevices::rgb(213,  94,   0, max = 255, alpha = (100 -  0) * 255 / 100)

  graphics::par(mfrow=c(1, 1 + insample), xpd=TRUE, cex.main=0.8)
  if (insample) {
    graphics::plot(ate_results[, 5],  ate_results[, 6]/ate_results[, 5],
                   lty=1, ylim = c(0, 1), type="l",
                   xlab="Remaining obs.", ylab="% Treated", main="Balance of Groups")
  }
  if (insample) main_text <- "Average Treatment Effect\nblue: ATE (95% CI)\n black-dash: ATT\n red-dash: ATU"
  else main_text <- "Average Treatment Effect\nblue: ATE (95% CI)"
  graphics::plot(ate_results[, 5], ate_results[, 1], type="l", lty = 1, lwd = 2,
                 ylim=c(min(ate_results[, 2]),
                        max(ate_results[, 3])),
                 xlim=c(0, length(prediction$var)),
                 col=Blue100,
                 xlab="Remaining obs.",
                 ylab="Treatment effect",
                 main=main_text)
  #axis(1, quantile.list, n.list)
  #graphics::lines(ate_results[, 5], rep(0, n.steps+1), lty=2, col="gray")
  #confidence interval:
  graphics::polygon(c(ate_results[, 5], rev(ate_results[, 5])),
                    c(ate_results[, 2], rev(ate_results[, 3])),
                    lty = 1, border = FALSE, density = 30, col = Blue50, angle = 45)

  if (insample) {
    graphics::lines(ate_results[, 5], att_results[, 1], lty=2)
    graphics::lines(ate_results[, 5], atu_results[, 1], lty=2, col=Vermillion100)
    #graphics::legend("bottomleft",
    #                 inset=c(-0.36, -0.3),
    #                 c("ATE", "ATT", "ATU"),
    #                 lty=c(1, 2, 2),
    #                 col=c(Blue100, "black", Vermillion100),
    #                 lwd=2,
    #                 cex=0.8)
  }

  #stats for within-sample:
  # % of treated
  # % of untreated

  # reset graphics setting
  graphics::par(mfrow=c(1,1))
  # return matrix with selection indices
  invisible(idx)
}
