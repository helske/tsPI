#' Compute Bayesian prediction intervals for Structural Time Series with exogenous variables using importance sampling
#'
#' Function \code{struct_pi} computes prediction intervals for structural time series
#' with exogenous variables using importance sampling.
#'
#' @export
#' @import KFAS
#' @name struct_pi
#
#' @param x vector containing the time series
#' @param xreg matrix or data frame containing the exogenous variables
#' (not including the intercept which is always included for non-differenced series)
#' @param type type of model. Possible options are \code{"level"}, \code{"trend"} and \code{"BSM"},
#' corresponding to local level, local linear trend, and local linear trend model with seasonal component.
#' @param nsim number of simulations used in importance sampling.
#' @param n_ahead length of the forecast horizon.
#' @param level desired frequentist coverage probability of the prediction intervals.
#' @param median compute the median of the prediction interval.
#' @param se_limits compute the standard errors of the prediction interval limits.
#' @param prior prior to be used in importance sampling for log-sd parameters.
#' Defaults to uniform prior on logarithm of standard deviations.
#' If "custom", a user-defined custom prior is used (see next arguments).
#' @param custom_prior function for computing custom prior.
#' First argument must be a vector containing the log-variance parameters.
#' @param custom_prior_args list containing additional arguments to \code{custom_prior}.
#' @param inits initial values for log-sds
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
#'
#'
#' pred_StructTS <- predict(StructTS(Nile, type ="level"), n.ahead = 10, se.fit = TRUE)
#' pred_StructTS <- cbind(pred = pred_StructTS$pred,
#'   lwr = pred_StructTS$pred - qnorm(0.975)*pred_StructTS$se,
#'  upr = pred_StructTS$pred + qnorm(0.975)*pred_StructTS$se)
#'
#' pred <- struct_pi(Nile, type = "level", n_ahead = 10)
#'
#' ts.plot(ts.union(Nile,pred_StructTS, pred[,1:3]), col = c(1,2,2,2,3,3,3),
#'   lty = c(1,1,2,2,1,2,2))
#'
#'
struct_pi <- function(x, type = c("level", "trend", "BSM"), xreg = NULL,
  n_ahead = 1, level = 0.95, median = TRUE, se_limits = TRUE,
  prior = "uniform", custom_prior, custom_prior_args, nsim = 1000, inits,
  last_only = FALSE, return_weights = FALSE, ...){

  distfkt <- function(a, prob, ex, sdx, w){
    sum(w * pnorm(q = a, mean = ex, sd = sdx)) - prob
  }

  type <- match.arg(type)
  prior <- match.arg(prior, c("uniform", "custom"))
  if (prior == "custom" && missing(custom_prior))
    stop("Missing custom prior.")

  n <- length(x)
  x <- window(x, end = end(x) + c(0, n_ahead), extend = TRUE)

  if (!is.null(xreg)) {
    model <- switch(type,
      level = SSModel(x ~ xreg + SSMtrend(1, NA), H = NA),
      trend = SSModel(x ~ xreg + SSMtrend(2, list(NA, NA)), H = NA),
      BSM = SSModel(x ~ xreg + SSMtrend(2, list(NA, NA)) + SSMseasonal(frequency(x), Q = NA), H = NA))
  } else {
    model <- switch(type,
      level = SSModel(x ~ SSMtrend(1, NA), H = NA),
      trend = SSModel(x ~ SSMtrend(2, list(NA, NA)), H = NA),
      BSM = SSModel(x ~ SSMtrend(2, list(NA, NA)) + SSMseasonal(frequency(x), Q = NA), H = NA))
  }
  npar <- switch(type,
    level = 2,
    trend = 3,
    BSM = 4)

  dx <- 1 + 0:(npar - 2) * npar

  likfn <- function(pars, model){
    # parameters are log(standard deviation)
    model$Q[dx] <- exp(2 * pars[-npar])
    model$H[1] <- exp(2 * pars[npar])
    - logLik(model)
  }
  fit <- optim(fn = likfn, par = if (missing(inits)) log(rep(sd(x, na.rm = TRUE), npar)) else inits,
    method = "BFGS", hessian = TRUE, model = model)

  if (any(exp(2*fit$par) < 1e-7))
    warning("Some of the variance parameters were estimated as smaller than 1e-7.
      Boundaries of parameter space can cause problems in importance sampling. Consider fixing variances to zero.")

  psihat <- as.numeric(fit$par)
  psivarchol <- try(t(chol(solve(fit$hessian))), TRUE)
  if (inherits(psivarchol, "try-error"))
    stop("Hessian obtained from optim is not positive definite.")
  psivarinv <- fit$hessian
  psisim <- array(rnorm(n = nsim * npar),c(nsim, npar))
  for (i in 1:nsim)
    psisim[i,] <- psihat + psivarchol %*% psisim[i, ]

  ex <- sdx <- matrix(0,n_ahead, nsim)

  w <- apply(psisim, 1, function(x) all(x < log(sqrt(1e7))))


  for (i in which(w)) {
    model$Q[dx] <- exp(2 * psisim[i, -npar])
    model$H[1] <- exp(2 * psisim[i, npar])
    out <- KFS(model, filtering = "mean", smoothing = "none")

    ex[1:n_ahead,i] <- out$m[(n + 1):(n + n_ahead)]
    sdx[1:n_ahead,i] <- sqrt(out$P_mu[(n + 1):(n + n_ahead)] + model$H[1])

    weight <- exp(out$logLik + fit$value)  /
      exp(-0.5 * t(psisim[i,] - psihat) %*% psivarinv %*% (psisim[i, ] - psihat))

    w[i] <- weight * switch(prior,
      uniform = 1,
      custom = do.call(custom_prior, list(psisim[i, ], custom_prior_args)))
  }
  w <- w/sum(w)
  out <- ts(matrix(NA, n_ahead, 2 + median + 2 * se_limits), end = end(model$y), frequency = frequency(model$y))
  colnames(out) <- c(if (median) "median", "lwr", "upr", if (se_limits) c("se_lwr", "se_upr"))

  if (sum(is.na(w)) == 0) {
    nz_w <- w[w != 0]
    for (j in 1:n_ahead) {
      nz_ex <- ex[j, w != 0]
      nz_sdx <- sdx[j, w != 0]
      interval <- c(mean(nz_ex) + c(-1, 1) * 8 * max(nz_sdx))
      out[j, "lwr"] <- uniroot(distfkt, interval = interval, prob = (1 - level) / 2,
        ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      out[j, "upr"] <- uniroot(distfkt,  interval = interval, prob = 1 - (1 - level) / 2,
        ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      if (median) {
        out[j, "median"] <- uniroot(distfkt,  interval = interval, prob = 0.5,
          ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      }
      if (se_limits) {
        out[j, "se_lwr"] <-
          sqrt(sum((nz_w * ((1 - level) / 2 - pnorm(q = out[j, 1], nz_ex, nz_sdx)))^2) / (nsim - 1)) /
          (sum(nz_w * dnorm(x = out[j, 1], nz_ex, nz_sdx)/sqrt(nsim)))
        out[j, "se_upr"] <-
          sqrt(sum((nz_w * ((1 - level) / 2 - pnorm(q = out[j, 2], nz_ex, nz_sdx)))^2) / (nsim - 1)) /
          (sum(nz_w * dnorm(x = out[j, 2], nz_ex, nz_sdx)/sqrt(nsim)))
      }
    }
  } else stop("NA values in weights.")


  if (return_weights)
    out <- list(pred = out, weights = w)
  out
}