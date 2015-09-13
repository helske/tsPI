context("Testing arima_pi")


test_that("bogus arguments throw error",{
  x <- rnorm(10)
  expect_error(arima_pi(x, phi = -2))
  expect_error(arima_pi(x, order = "ar"))
  expect_error(arima_pi(x, xreg = 1))
  expect_error(arima_pi(x, n_ahead = -1))
  expect_error(arima_pi(x, level = 2))
  expect_error(arima_pi(x, median = "true"))
  expect_error(arima_pi(x, prior = "custom"))
  expect_error(arima_pi(x, prior = "custom", custom_prior = "f"))
  expect_error(arima_pi(x, nsim = 0))
})

test_that("output of arima_pi is of correct size and form",{
  x <- ts(rnorm(10), start = 2000, frequency = 12)
  pred <- arima_pi(x, order = c(1, 0, 0), n_ahead = 5, nsim = 50)
  expect_identical(dim(pred), c(5L, 5L))
  expect_identical(class(pred), c("mts", "ts", "matrix"))
  expect_identical(frequency(pred), frequency(x))
  expect_identical(start(pred), end(x)+c(0, 1))
})

test_that("same seeds give same results",{
    x <- rnorm(10)
    set.seed(1)
    pred1 <- arima_pi(x, c(1, 0, 0), nsim = 50)
    set.seed(1)
    pred2 <- arima_pi(x, c(1, 0, 0), nsim = 50)
    expect_identical(pred1, pred2)
})

test_that("larger nsim gives smaller se",{
  x <- rnorm(10)
  set.seed(1)
  pred1 <- arima_pi(x, c(1, 0, 0), nsim = 50)
  set.seed(1)
  pred2 <- arima_pi(x, c(1, 0, 0), nsim = 100)
  expect_more_than(pred1[, "se_upr"], pred2[, "se_upr"])
})

test_that("arima_pi gives same results each time",{
  set.seed(1)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50)
  expect_equal(pred[1,"median"], 2.707872101, tol = 1e-4, check.attributes = FALSE)
  expect_equal(pred[1,"lwr"], 1.809644512, tol = 1e-4, check.attributes = FALSE)
  expect_equal(pred[1,"upr"], 3.606841626, tol = 1e-4, check.attributes = FALSE)
  expect_equal(pred[1,"se_lwr"], 0.07736541451, tol = 1e-4, check.attributes = FALSE)
  expect_equal(pred[1,"se_upr"], 0.01355272753, tol = 1e-4, check.attributes = FALSE)
})
