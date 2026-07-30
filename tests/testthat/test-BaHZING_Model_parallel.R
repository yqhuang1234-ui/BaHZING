test_that("parallel chain execution matches serial execution given the same seed", {
  skip_on_cran()

  data("iHMP_Reduced")
  formatted_data <- Format_BaHZING(iHMP_Reduced)
  x <- c("soft_drinks_dietnum", "diet_soft_drinks_dietnum")

  shared_args <- list(
    formatted_data = formatted_data,
    x = x,
    exposure_standardization = "standard_normal",
    n.chains = 2,
    n.adapt = 60,
    n.iter.burnin = 2,
    n.iter.sample = 2,
    seed = 123,
    verbose = FALSE
  )

  serial_results <- do.call(BaHZING_Model, c(shared_args, list(parallel = FALSE)))
  parallel_results <- do.call(BaHZING_Model, c(shared_args, list(parallel = TRUE, n.cores = 2)))

  testthat::expect_equal(serial_results, parallel_results)
})
