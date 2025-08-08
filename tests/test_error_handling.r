# Test di gestione errori per prepare_omics()

test_that("prepare_omics handles invalid file paths", {
  expect_error(
    prepare_omics(data = "nonexistent_file.csv", id_col = "ID"),
    class = "error"
  )
})

test_that("prepare_omics handles missing id_col", {
  test_data <- data.frame(x = 1:3, y = 4:6)
  
  expect_error(
    prepare_omics(data = test_data, id_col = "missing_column"),
    class = "error"
  )
})

test_that("prepare_omics handles invalid cohort_filter", {
  test_data <- data.frame(
    Individual.Id = c(1, 2),
    Feature1 = c(1.1, 2.1)
  )
  
  cohort_data <- data.frame(
    Individual.Id = c(1, 2),
    Group = c("A", "B")
  )
  
  expect_error(
    prepare_omics(
      data = test_data,
      id_col = "Individual.Id",
      cohort = cohort_data,
      cohort_filter = "invalid R syntax here"
    ),
    class = "error"
  )
})

test_that("prepare_omics handles unsupported input types", {
  expect_error(
    prepare_omics(data = "not a dataframe", id_col = "ID"),
    "Unsupported input type"
  )
})

test_that("prepare_omics handles invalid transpose mode", {
  test_data <- data.frame(Individual.Id = 1, Feature1 = 1.1)
  
  expect_error(
    prepare_omics(data = test_data, id_col = "Individual.Id", transpose = "invalid"),
    class = "error"
  )
})

test_that("prepare_omics handles invalid deduplicate mode", {
  test_data <- data.frame(
    Individual.Id = c(1, 1),
    Time.to.IA = c(-6, -6),
    Feature1 = c(1.1, 1.2)
  )
  
  expect_error(
    prepare_omics(
      data = test_data, 
      id_col = "Individual.Id", 
      time_col = "Time.to.IA",
      deduplicate = "invalid"
    ),
    class = "error"
  )
})