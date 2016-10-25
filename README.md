[![Build Status](https://travis-ci.org/helske/tsPI.png?branch=master)](https://travis-ci.org/helske/tsPI)
[![codecov.io](http://codecov.io/github/helske/tsPI/coverage.svg?branch=master)](http://codecov.io/github/helske/tsPI?branch=master)
[![downloads](http://cranlogs.r-pkg.org/badges/tsPI)](http://cranlogs.r-pkg.org/badges/tsPI)
[![cran version](http://www.r-pkg.org/badges/version/tsPI)](http://cran.r-project.org/package=tsPI)

# tsPI: R package for improved prediction intervals for ARIMA and structural time series models with exogenous variables

Package tsPI computes prediction intervals for ARIMA and Gaussian structural time series models by using importance sampling approach with uninformative priors for model parameters, leading to more accurate coverage probabilities in frequentist sense. Instead of sampling the future observations and hidden states of the state space representation of the model, only model parameters are sampled, and the method is based solving the equations corresponding to the conditional coverage probability of the prediction intervals. This makes method relatively fast compared to for example MCMC methods, and standard errors of prediction limits can also be computed straightforwardly.

ARIMA case is based on articles 

- Jouni Helske and Jukka Nyblom. Improved frequentist prediction intervals for autoregressive models by simulation. In Siem Jan Koopman and Neil Shephard, editors, Unobserved Components and Time Series Econometrics. Oxford University Press, 2015.
- Jouni Helske and Jukka Nyblom. [Improved frequentist prediction intervals for ARMA models by simulation](http://www.uva.fi/materiaali/pdf/isbn_978-952-476-523-7.pdf). In Johan Knif and Bernd Pape, editors, Contributions to Mathematics, Statistics, Econometrics, and Finance: essays in honour of professor Seppo Pynnönen, number 296 in Acta Wasaensia, pages 71–86. University of Vaasa, 2014. 

Structural time series model case is based on a straightforward generalization presented in 
- Helske, J. (2015). Prediction and interpolation of time series by state space models. University of Jyväskylä.  [PhD thesis.](https://jyx.jyu.fi/dspace/handle/123456789/49043)

### Example: 95 % prediction intervals for AR(1) process ###

```{r, fig.height = 4, fig.width = 8}
library(tsPI)
library(KFAS) #for plug-in intervals

set.seed(12345)
x <- arima.sim(n = 30, model = list(ar = 0.9))
fit <- arima(x, c(1, 0, 0))
model <- SSModel(x ~ SSMarima(ar = fit$coef[1], Q = fit$sigma), H = 0)
model$P1inf[1,1] <- 0
model$a1[1] <- fit$coef[2]
pred_plugin <- predict(model, n.ahead = 10, interval = "prediction")
pred_tspi  <- arima_pi(x, c(1, 0, 0), n_ahead = 10, se_limits = FALSE, nsim = 1000)

ylim <- round(range(c(pred_plugin, pred_tspi, x)) + c(-1, 1))

plot(ts.union(x, pred_plugin, pred_tspi), plot.type = "single",
  col = c(1, 2, 2, 2, 4, 4, 4), pch = c(19, 15, 15, 15, 15, 15, 15),
  lty = c(1, 1, 1, 1, 1), type = "b",
  ylim = ylim, xlab = "time", ylab= "value", axes = FALSE)
axis(1, at = 1:40, labels = 1:40)
axis(2, at = ylim[1]:ylim[2])
legend("topleft", c("observations", "plug-in",  "tsPI"), 
  lty = 1, pch= c(19, 15, 15), col = c(1, 2, 4))
```  
![imfs](https://github.com/helske/tsPI/blob/master/ar1.png)

### Installing tsPI ###

Package is now available at CRAN. The latest development version can be installed from the github using the `devtools` package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/tsPI")
```
