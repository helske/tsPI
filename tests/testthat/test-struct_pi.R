context("Testing struct_pi")

set.seed(1)

test_that("bogus arguments throw error",{
  x <- rnorm(10)
  expect_error(struct_pi(x, type = "arima"))
  expect_error(struct_pi(x, xreg = 1))
  expect_error(struct_pi(x, n_ahead = -1))
  expect_error(struct_pi(x, level = 2))
  expect_error(struct_pi(x, median = "true"))
  expect_error(struct_pi(x, prior = "custom"))
  expect_error(struct_pi(x, prior = "custom", custom_prior = "f"))
  expect_error(struct_pi(x, nsim = 0))
})

test_that("output of struct_pi is of correct size and form",{
  x <- ts(rnorm(10), start = 2000, frequency = 12)
  pred <- struct_pi(x, n_ahead = 5, nsim = 50)
  expect_identical(dim(pred), c(5L, 5L))
  expect_identical(class(pred), c("mts", "ts", "matrix"))
  expect_identical(frequency(pred), frequency(x))
  expect_identical(start(pred), end(x)+c(0, 1))
})

test_that("same seeds give same results",{
    x <- rnorm(10)
    set.seed(1)
    pred1 <- struct_pi(x, nsim = 50)
    set.seed(1)
    pred2 <- struct_pi(x, nsim = 50)
    expect_identical(pred1, pred2)
})

test_that("larger nsim gives smaller se",{
  x <- rnorm(10)
  set.seed(1)
  pred1 <- struct_pi(x, nsim = 50)
  set.seed(1)
  pred2 <- struct_pi(x, nsim = 100)
  expect_more_than(pred1[, "se_upr"], pred2[, "se_upr"])
})

test_that("struct_pi gives same results each time",{
  set.seed(1)
  pred <- struct_pi(Nile, nsim = 50)
  expect_equivalent(pred[1,"median"], 800.20719)
  expect_equivalent(pred[1,"lwr"], 505.08147)
  expect_equivalent(pred[1,"upr"], 1092.89858)
  expect_equivalent(pred[1,"se_lwr"], 26.176364)
  expect_equivalent(pred[1,"se_upr"], 5.3530965)

})
