#' Calculate robust average treatment effect using posterior uncertainty
#'
#' @param X The covariates used in the causal regression.
#'
#' @return A trained causal forest object.
#'
#' @examples \dontrun{
#'
#' }
#'
#' @export

robust_treatment <- function(object, newX, n.steps = 10) {
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
  discard.steps <- quantile(prediction$var, probs = quantile.list)

  ate_results <- matrix(NaN, n.steps+1, 5 + 2 * insample) #includes other baseline results

  for (i in 1:(n.steps+1)) {
    idx <- c(prediction$var <= discard.steps[i])
    ate_results[i, 5] <- sum(idx)

    if(missing(newX)) {
      # within-of-sample
      invisible(predace_robust <- predict.ace(object,
                                              newX = matrix(object$train_data$X[idx, ], sum(idx), ncol((object$train_data$X))),
                                              newZ = object$train_data$Z[idx],
                                              marginal=TRUE, causal=TRUE, normalize=FALSE))
    } else {
      # out-of-sample
      predace_robust <- predict.ace(object, newX = newX[idx, ], newZ = 1, marginal=TRUE, causal=TRUE)
    }
    ate_results[i, 1:4] <- unlist(predace_robust$ate)
    if(insample) {
      #remaining treated:
      ate_results[i, 6] <- t(idx) %*% object$train_data$Z
      #remaining untreated:
      ate_results[i, 7] <- t(idx) %*% (1 - object$train_data$Z)
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

  graphics::par(mfrow=c(1, 1 + insample), xpd=TRUE)
  if (insample) {
    plot(ate_results[, 5],  ate_results[, 6]/ate_results[, 5],
         lty=1, ylim = c(0, 1), type="l",
         xlab="Remaining obs.", ylab="Ratio of treated")
  }

  graphics::plot(ate_results[, 5], ate_results[, 1], type="l", lty = 1, lwd = 2,
                 ylim=c(min(ate_results[, 2]),
                        max(ate_results[, 3])),
                 xlim=c(0, length(prediction$var)),
                 col=Blue100,
                 xlab="Remaining obs.",
                 ylab="Avg Treatment Effect (95% CI)")
  #axis(1, quantile.list, n.list)
  #graphics::lines(ate_results[, 5], rep(0, n.steps+1), lty=2, col="gray")
  #confidence interval:
  graphics::polygon(c(ate_results[, 5], rev(ate_results[, 5])),
                    c(ate_results[, 2], rev(ate_results[, 3])),
                    lty = 1, border = FALSE, density = 30, col = Blue50, angle = 45)

  if (insample) {
    lines(ate_results[, 5], att_results[, 1], lty=2)
    lines(ate_results[, 5], atu_results[, 1], lty=2, col=Vermillion100)
    legend("bottomleft",
           inset=c(-0.36, -0.3),
           c("ATE", "ATT", "ATU"),
           lty=c(1, 2, 2),
           col=c(Blue100, "black", Vermillion100),
           lwd=2,
           cex=0.8)
  }

  #stats for within-sample:
  # % of treated
  # % of untreated

}
