#' Compute the average coverage of the prediction intervals computed by
#' naive plug-in method and \code{\link{struct_pi}}
#'
#' Computes expected coverage probabilities of the prediction intervals of structural time series model.
#'
#' @export
#' @seealso \code{\link{struct_pi}}.
#
#' @param type Type of model. See \code{\link{struct_pi}}.
#' @param sds vector containing the standard deviations of the model.
#' @param n length of the time series
#' @param n_ahead length of the forecast horizon
#' @param nsim number of simulations used in importance sampling
#' @param nsim2 number of simulations used in computing the expected coverage
#' @param level desired coverage probability of the prediction intervals
#' @param prior prior to be used in importance sampling.
#' @param return_all_coverages return raw results i.e. coverages for each simulations. When \code{FALSE} (default), summary statistics are returned.
#' @param ... additional arguments to \code{\link{struct_pi}}.
#' @return a list containing the coverage probabilities
avg_coverage_struct <- function(type = c("level", "trend", "BSM"), sds, n, n_ahead = 1,
  nsim2, nsim = 100, level = 0.95, prior = "uniform", return_all_coverages = FALSE, ...){



  prior <- match.arg(prior, c("uniform", "custom"))

  type <- match.arg(type)

  npar <- switch(type,
    level = 2,
    trend = 3,
    BSM = 4)
  if (npar != length(sds))
    stop("Incorrect number of standard deviations provided.")

  model <- true_model <- switch(type,
    level = SSModel(rep(NA,n) ~ SSMtrend(1, NA), H = NA),
    trend = SSModel(rep(NA,n) ~ SSMtrend(2, list(NA, NA)), H = NA),
    BSM = SSModel(rep(NA,n) ~ SSMtrend(2, list(NA, NA)) + SSMseasonal(frequency(x), Q = NA), H = NA))
  model$y[1] <- 0

  covprobs <- array(0, c(n_ahead, nsim2, 2))
  dimnames(covprobs)[[3]] <- c("plug-in", prior)

  dx <- 1 + 0:(npar - 2) * npar

  likfn <- function(pars, model){
    # parameters are log(standard deviation)
    model$Q[dx] <- exp(2 * pars[-1])
    model$H[1] <- exp(2 * pars[1])
    - logLik(model)
  }

  true_model$Q[dx] <- exp(2 * log(sds[-1]))
  true_model$H[1] <- exp(2 * log(sds[1]))

  for (i in 1:nsim2) {
    true_model$y[-1] <- NA
    model$y[] <- true_model$y[] <- simulateSSM(true_model, "obs")
    # use true parameters as initial values, in real application multiple initial values would be used
    fit <- optim(fn = likfn, par = log(sds), method = "BFGS", model = model)

    model$Q[dx] <- exp(2 * fit$par[-1])
    model$H[1] <- exp(2 * fit$par[1])

    pred <- predict(model, n.ahead = n_ahead, level = level, interval = "prediction")

    true_pred <- predict(true_model, n.ahead = n_ahead, se.fit = TRUE)


    ipi <- try(struct_pi(model$y, type = type, inits = fit$par, nsim = nsim, level = level, n_ahead = n_ahead,
      prior = prior, median = FALSE, se_limits = FALSE, ...), TRUE)
    if (!inherits(ipi, "try-error")) {

      covprobs[, i, 2] <- pnorm(q = ipi[, "upr"], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2)) -
        pnorm(q = ipi[, "lwr"], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2))
    }
    covprobs[, i, 1] <- pnorm(q = pred[, 3],mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2)) -
      pnorm(q = pred[, 2], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2))

  }
  if(!return_all_coverages){
    out <- vector("list", 2)
    names(out) <- c("plugin", prior)
    for(i in 1:2){
      out[[i]] <- matrix(0, n_ahead, 7)
      colnames(out[[i]]) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.","se")
      rownames(out[[i]]) <- paste("n_ahead =", 1:n_ahead)
      for(j in 1:n_ahead)
        out[[i]][j,] <- c(summary(covprobs[j,,i]), sd(covprobs[j,,i])/sqrt(nsim2))
    }
    return(out)
  } else return(covprobs)
}
