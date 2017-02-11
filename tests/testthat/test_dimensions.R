source("main.R")

context("Checking dimension of output objects.")

test_that("Dimension of gam_vb is compatible with the data", {
  expect_equal(dim(vb$gam_vb), c(p, d))
  expect_equal(dim(vb_z$gam_vb), c(p, d))
  expect_equal(dim(vb_logit$gam_vb), c(p, d))
  expect_equal(dim(vb_logit_z$gam_vb), c(p, d))
})


test_that("Length of om_vb is compatible with the data", {
  expect_equal(length(vb$om_vb), p)
  expect_equal(length(vb_z$om_vb), p)
  expect_equal(length(vb_logit$om_vb), p)
  expect_equal(length(vb_logit_z$om_vb), p)
})
