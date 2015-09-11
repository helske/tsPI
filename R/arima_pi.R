#' Compute Bayesian prediction intervals for ARIMA processes with exogenous variables using importance sampling
#'
#' Function \code{arima_pi} computes prediction intervals for ARIMA processes
#' with exogenous variables using importance sampling. For regression coefficients,
#' diffuse (uninformative) prior is used, whereas multiple options for
#' prior distributions for ARMA coefficients are supported.
#'
#' @export
#' @import KFAS
#' @name arima_pi
#
#' @param y vector containing the time series
#' @param xreg matrix or data frame containing the exogenous variables
#' (not including the intercept which is always included for non-differenced series)
#' @param order vector of length 3 with values p,d,q
#' corresponding to the number of AR parameters, degree of differencing and number of MA parameters.
#' @param nsim number of simulations used in importance sampling.
#' @param n.ahead length of the forecast horizon.
#' @param level desired frequentist coverage probability of the prediction intervals.
#' @param median compute the median of the prediction interval.
#' @param se_limits compute the standard errors of the prediction interval limits.
#' @param priors priors to be used in importance sampling. Multiple choices are allowed. See \code{\link{jeffreys}} for details.
#' @param custom_prior function for computing custom prior.
#' First argument must be a vector containing the AR and MA parameters (in that order).
#' @param custom_prior_args list containing additional arguments to \code{custom_prior}.
#' @param last_only compute the prediction intervals only for the last prediction step.
#' @param return_weights Return (scaled) weights used in importance sampling.
#' @param ... Additional arguments for \code{arima}.
#' @return a list containing the prediction intervals.
#'  @references
#' \enumerate{
#' \item{Helske, J. and Nyblom, J. (2013). Improved frequentist prediction
#' intervals for autoregressive models by simulation.
#' In Siem Jan Koopman and Neil Shephard, editors,
#' Unobserved Components and Time Series Econometrics. Oxford University Press. In press.}
#'  \item{Helske, J. and Nyblom, J. (2014). Improved frequentist prediction intervals for
#'  ARMA models by simulation.
#'  In Johan Knif and Bernd Pape, editors,
#' Contributions to Mathematics, Statistics, Econometrics, and Finance:
#' essays in honour of professor Seppo Pynnonen,
#' number 296 in Acta Wasaensia, pages 71–86. University of Vaasa.}
#' }
#' @examples
#'
#' pred_arima <- predict(arima(lh, order = c(3,0,0)), n.ahead = 12, se.fit = TRUE)
#' pred_arima <- cbind(pred = pred_arima$pred,
#'  lwr = pred_arima$pred - qnorm(0.975)*pred_arima$se^2,
#'  upr = pred_arima$pred + qnorm(0.975)*pred_arima$se^2)
#'
#' pred <- arima_pi(lh, order = c(3,0,0), n.ahead = 12, prior = "uniform")
#'
#' ts.plot(ts.union(lh,pred_arima, pred$uniform[,1:3]), col = c(1,2,2,2,3,3,3),
#' lty = c(1,1,2,2,1,2,2))
arima_pi <- function(y, order, xreg = NULL, n.ahead = 1, level = 0.95, median = TRUE, se_limits = TRUE,
  priors = c("uniform", "approx_joint_jeffreys", "approx_marginal_jeffreys", "exact_joint_jeffreys", "exact_marginal_jeffreys"),
  custom_prior, custom_prior_args, nsim = 1000, last_only = FALSE, return_weights = FALSE, ...){

  distfkt <- function(a, prob, ey, sdy, w){
    sum(w * pnorm(q = a, mean = ey, sd = sdy)) - prob
  }

  priors <- match.arg(priors, several.ok = TRUE)
  if (!missing(custom_prior))
    priors <- c(priors, "custom")
  n <- length(y)
  fit <- arima(y, order, xreg = xreg[1:n, ], ...)
  if (fit$code != 0 || sum(diag(fit$var.coef) < 1e-7) != 0)
  {
    stop("arima function returned non-convergence or coefficient variances smaller than 1e-7.")
  }

  p <- order[1]
  d <- order[2]
  q <- order[3]
  npar <- p + q
  psihat <- as.numeric(fit$coef[1:npar])
  psivar <- matrix(fit$var.coef[1:npar, 1:npar], npar, npar)
  psivarchol <- try(t(chol(psivar)), TRUE)
  if (class(psivarchol) == "try-error" || any(diag(psivarchol) < 1e-12))
    stop("Covariance matrix fit$var.coef[1:npar, 1:npar] obtained from arima function is not positive definite.")

  psisim <- array(rnorm(n = nsim * npar),c(nsim, npar))
  for (i in 1:nsim)
    psisim[i,] <- psihat + psivarchol %*% psisim[i, ]

  sigma2hat <- fit$sigma2
  k <- length(fit$coef) - npar
  sigmasim <- 1/rchisq(nsim, df = n - k)

  #stationarity
  if (p > 0) {
    w <- apply(matrix(abs(apply(cbind(rep(1, nsim),-psisim[, 1:p]), 1, polyroot)) > 1, p, nsim), 2, sum) == p
  } else w <- rep(1, nsim)
  #invertibility
  if (q > 0) {
    w <- w * apply(matrix(abs(apply(cbind(rep(1, nsim), psisim[, (p + 1):(p + q)]),1, polyroot)) > 1, q, nsim), 2, sum) == q
  }

  y <- c(y, rep(NA, n.ahead))
  ey <- sdy <- matrix(0,n.ahead, nsim)
  psivarinv <- crossprod(solve(t(psivarchol)))

  valid <- which(w > 0)[1]
  valid_ar <- if (p > 0) psisim[valid, 1:p] else NULL
  valid_ma <- if ( q > 0) psisim[valid, (p + 1):(p + q)] else NULL
  if (is.null(xreg)) {
    model <- SSModel(y ~ SSMarima(ar = valid_ar, ma = valid_ma, d = order[2], Q = 1), H = 0)
  } else {
    model <- SSModel(y ~ xreg + SSMarima(ar = valid_ar, ma = valid_ma, d = order[2], Q = 1), H = 0)
  }
  kd <- k + d
  m <- max(p, q + 1)
  nd <- (kd + 1):(kd + m)

  w <- matrix(w, nsim, length(priors))
  colnames(w) <- priors
  for (i in which(w[, 1] > 0)) {
    if (q > 0)
      model$R[(kd + 2):(kd + 1 + q)] <- psisim[i, (p + 1):npar]
    if (p > 0)
      model$T[(kd + 1):(kd + p),kd + 1,] = psisim[i, 1:p]
    model$P1[nd, nd] <- solve(a = diag(m^2) - as.matrix(kronecker(model$T[nd, nd, ], model$T[nd, nd, ])),
      b = c(model$R[nd] %*% t(model$R[nd])))

    out <- KFS(model, filtering = "mean", smoothing = "none")
    if (out$d > 0) {
      s2 <- sum(out$v[1:out$d][out$Finf == 0]^2/out$F[1:out$d][out$Finf == 0]) +
        sum(out$v[(out$d + 1):n]^2/out$F[(out$d + 1):n])
    } else s2 <- sum(c(out$v[1:n])^2/out$F[1:n])

    sigmasim[i] <-  sqrt(s2 * sigmasim[i])
    ey[1:n.ahead,i] <- out$m[(n + 1):(n + n.ahead)]
    sdy[1:n.ahead,i] <- sqrt(out$P_mu[(n + 1):(n + n.ahead)]) * sigmasim[i]

    #s2 is replaced by (s2/n)/sigma2hat so that w>>0
    detVXinvX <- prod(c(out$Finf[out$Finf > 0], out$F[1:out$d][out$Finf == 0], out$F[(out$d + 1):n]))
    weight <-
      (((s2 / n) / sigma2hat)^(-0.5 * (n - k)) / sqrt(detVXinvX)) /
      exp(-0.5 * t(psisim[i,] - psihat) %*% psivarinv %*% (psisim[i, ] - psihat))

    if("uniform" %in% priors)
      w[i, "uniform"] <- weight
    if("exact_joint_jeffreys" %in% priors)
      w[i, "exact_joint_jeffreys"] <- weight * exact_joint_jeffreys(psisim[i, ], xreg, p, q, n)
    if("exact_marginal_jeffreys" %in% priors)
      w[i, "exact_marginal_jeffreys"] <- weight * exact_marginal_jeffreys(psisim[i, ], p, q, n)
    if("approx_joint_jeffreys" %in% priors)
      w[i, "approx_joint_jeffreys"] <- weight * approx_joint_jeffreys(psisim[i, ], xreg, p, q, n)
    if("approx_marginal_jeffreys" %in% priors)
      w[i, "approx_marginal_jeffreys"] <- weight * approx_marginal_jeffreys(psisim[i, ], p, q)
    if("custom" %in% priors)
      w[i, "custom"] <- weight * do.call(custom_prior, list(psisim[i, ], custom_prior_args))
  }

  out <- vector("list",length(priors))
  names(out) <- priors
  for (i in 1:length(priors)) {
    w[, i] <- w[, i]/sum(w[, i])
    if (!last_only) {
      out[[i]] <- ts(matrix(0, n.ahead, 2 + median + 2 * se_limits), end = end(model$y), frequency = frequency(model$y))
      colnames(out[[i]]) <- c(if (median) "median", "lwr", "upr", if (se_limits) c("se_lwr", "se_upr"))
      if (sum(is.na(w[, i])) == 0)
      {
        nz_w <- w[w[, i] != 0, i]
        for (j in 1:n.ahead)
        {
          nz_ey <- ey[j, w[, i] != 0]
          nz_sdy <- sdy[j, w[, i] != 0]
          interval <- c(mean(nz_ey) + c(-1, 1) * 8 * max(nz_sdy))
          out[[i]][j, "lwr"] <- uniroot(distfkt, interval = interval, prob = (1 - level) / 2,
            ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
          out[[i]][j, "upr"] <- uniroot(distfkt,  interval = interval, prob = 1 - (1 - level) / 2,
            ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
          if (median) {
            out[[i]][j, "median"] <- uniroot(distfkt,  interval = interval, prob = 0.5,
              ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
          }
          if (se_limits) {
            out[[i]][j, "se_lwr"] <-
              sqrt(sum((nz_w * (((1 - level) / 2) - pnorm(q = out[[i]][j, 1], nz_ey, nz_sdy)))^2) / (nsim - 1)) /
              (sum(nz_w * dnorm(x = out[[i]][j, 1], nz_ey, nz_sdy)/sqrt(nsim)))
            out[[i]][j, "se_upr"] <-
              sqrt(sum((nz_w * ((1 - (1 - level) / 2) - pnorm(q = out[[i]][j, 2], nz_ey, nz_sdy)))^2) / (nsim - 1)) /
              (sum(nz_w * dnorm(x = out[[i]][j, 2], nz_ey, nz_sdy)/sqrt(nsim)))
          }
        }
      }

    } else {
      out[[i]] <- ts(numeric(2 + median + 2 * se_limits), end = end(model$y), frequency = frequency(model$y))
      names(out[[i]]) <- c(if (median) "median", "lwr", "upr", if (se_limits) c("se_lwr", "se_upr"))
      if (sum(is.na(w[, i])) == 0)
      {
        nz_w <- w[w[, i] != 0, i]
        nz_ey <- ey[w[, i] != 0]
        nz_sdy <- sdy[w[, i] != 0]
        interval <- c(mean(nz_ey) + c(-1, 1) * 8 * max(nz_sdy))
        out[[i]]["lwr"] <- uniroot(distfkt, interval = interval, prob = (1 - level) / 2,
          ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
        out[[i]]["upr"] <- uniroot(distfkt,  interval = interval, prob = 1 - (1 - level) / 2,
          ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
        if (median) {
          out[[i]]["median"] <- uniroot(distfkt,  interval = interval, prob = 0.5,
            ey = nz_ey, sdy = nz_sdy,w = nz_w, tol = 1e-12)$root
        }
        if (se_limits) {
          out[[i]]["se_lwr"] <-
            sqrt(sum((nz_w * (((1 - level) / 2) - pnorm(q = out[[i]][1], nz_ey, nz_sdy)))^2) / (nsim - 1)) /
            (sum(nz_w * dnorm(x = out[[i]][1], nz_ey, nz_sdy)/sqrt(nsim)))
          out[[i]]["se_upr"] <-
            sqrt(sum((nz_w * ((1 - (1 - level) / 2) - pnorm(q = out[[i]][2], nz_ey, nz_sdy)))^2) / (nsim - 1)) /
            (sum(nz_w * dnorm(x = out[[i]][2], nz_ey, nz_sdy)/sqrt(nsim)))
        }

      }
    }
  }
  if (return_weights)
    out <- list(pred = out, weights = w)
  out$fit <- fit
  out
}