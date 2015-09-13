[![Build Status](https://travis-ci.org/helske/tsPI.png?branch=master)](https://travis-ci.org/helske/tsPI)
[![cran version](http://www.r-pkg.org/badges/version/tsPI)](http://cran.r-project.org/package=tsPI)

# tsPI: R package for improved prediction intervals for ARIMA and structural time series models with exogenous variables

Prediction intervals for ARIMA and structural time series models by using importance sampling approach with uninformative priors for model parameters, leading to more accurate coverage probabilities in frequentist sense. Instead of sampling the future observations and hidden states of the state space representation of the model, only model parameters are sampled, and the method is based solving the equations corresponding to the conditional coverage probability of the prediction intervals. This makes method relatively fast compared to for example MCMC methods, and standard errors of prediction limits can also be computed straightforwardly.

ARIMA case is based on articles 

- Jouni Helske and Jukka Nyblom. Improved frequentist prediction intervals for autoregressive models by simulation. In Siem Jan Koopman and Neil Shephard, editors, Unobserved Components and Time Series Econometrics. Oxford University Press, 2015. In press.
- Jouni Helske and Jukka Nyblom. Improved frequentist prediction intervals for arma models by simulation. In Johan Knif and Bernd Pape, editors, Contributions to Mathematics, Statistics, Econometrics, and Finance: essays in honour of professor Seppo Pynnönen, number 296 in Acta Wasaensia, pages 71–86. University of Vaasa, 2014. http://www.uva.fi/materiaali/pdf/isbn_978-952-476-523-7.pdf

Structural time series model case is based on a straightforward generalization presented in 
- Helske, J. (2015). Prediction and interpolation of time series by state space models. University of Jyväskylä. PhD thesis. In press.

