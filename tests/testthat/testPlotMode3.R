test_that("Error if FactorsX$Mode3 is missing", {
  fake <- list(FactorsX = list())
  expect_error(plot_nplsda_blockX_mode3(fake))
})

test_that("Error if pcs request more components than available", {
  fake <- list(FactorsX = list(Mode3 = matrix(1:6, nrow = 3, ncol = 2)))
  expect_error(plot_nplsda_blockX_mode3(fake, pcs = c(1, 3)))
})

test_that("Output (invisible) is a proper data.frame", {
  M3 <- matrix(1:12, nrow = 4, ncol = 3)
  rownames(M3) <- paste0("tp", 1:4)
  fake <- list(FactorsX = list(Mode3 = M3))
  out <- plot_nplsda_blockX_mode3(fake, pcs = c(1, 2), edge = c(0.1, 0.1))
  
  expect_s3_class(out, "data.frame")
  expect_equal(ncol(out), 3)
  expect_named(out, c("time", "x", "y"))
  expect_equal(nrow(out), 4)
})

