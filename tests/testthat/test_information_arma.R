context("Testing information_arma")


test_that("bogus arguments throw error",{
  expect_error(information_arma(NA, Inf))
})

test_that("output of information_arma is of correct size and form",{
  mat <- information_arma(0.9, NULL)

})
