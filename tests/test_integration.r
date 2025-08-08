# Test di integrazione per prepare_omics() 

test_that("prepare_omics integration: basic 3D workflow", {
  # Dati sintetici controllati
  test_data <- data.frame(
    Individual.Id = c(1, 1, 2, 2, 3, 3),
    Time.to.IA = c(-6, -3, -6, -3, -6, -3),
    Feature1 = c(1.1, 1.2, 2.1, 2.2, 3.1, 3.2),
    Feature2 = c(4.1, 4.2, 5.1, 5.2, 6.1, 6.2)
  )
  
  result <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id",
    time_col = "Time.to.IA",
    transpose = "never"
  )
  
  # Verifica output
  expect_true(validate_omics_output(result))
  expect_equal(dim(result), c(3, 2, 2))  # 3 soggetti, 2 feature, 2 tempi
  expect_equal(dimnames(result)$Subject, c("1", "2", "3"))
  expect_equal(dimnames(result)$Feature, c("Feature1", "Feature2"))
  expect_equal(dimnames(result)$Time, c("-6", "-3"))
})

test_that("prepare_omics integration: workflow without time returns NULL (current behavior)", {
  test_data <- data.frame(
    Individual.Id = c(1, 2, 3),
    Feature1 = c(1.1, 2.1, 3.1),
    Feature2 = c(1.2, 2.2, 3.2)
  )
  
  result <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id",
    transpose = "never"
  )
  
  # Attualmente .build_tensor restituisce NULL per dati senza time_col
  expect_null(result)
})

test_that("prepare_omics integration: returns data.frame when tensor=FALSE", {
  test_data <- data.frame(
    Individual.Id = c(1, 2, 3),
    Feature1 = c(1.1, 2.1, 3.1),
    Feature2 = c(1.2, 2.2, 3.2)
  )
  
  result <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id",
    tensor = FALSE,  # Bypassa .build_tensor
    transpose = "never"
  )
  
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
})

test_that("prepare_omics integration: cohort filtering with 3D data", {
  test_data <- data.frame(
    Individual.Id = c(1, 1, 2, 2, 3, 3),
    Time.to.IA = c(-6, -3, -6, -3, -6, -3),
    Feature1 = c(1.1, 1.2, 2.1, 2.2, 3.1, 3.2)
  )
  
  cohort_data <- data.frame(
    Individual.Id = c(1, 2, 3),
    Group = c("A", "B", "A"),
    Model.or.Validation = c("Model", "Validation", "Model")
  )
  
  result <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id",
    time_col = "Time.to.IA",
    cohort = cohort_data,
    cohort_filter = "Model.or.Validation == 'Model'",
    transpose = "never"
  )
  
  expect_equal(dim(result)[1], 2)  # Solo soggetti 1 e 3
  expect_equal(dimnames(result)$Subject, c("1", "3"))
})

test_that("prepare_omics integration: transpose workflow", {
  # Dati in formato "wide" da trasporre
  wide_data <- data.frame(
    row_names = c("Individual.Id", "Time.to.IA", "Feature1", "Feature2"),
    Subj_1_T1 = c("1", "-6", "1.1", "4.1"),
    Subj_1_T2 = c("1", "-3", "1.2", "4.2"),
    Subj_2_T1 = c("2", "-6", "2.1", "5.1")
  )
  
  result <- prepare_omics(
    data = wide_data,
    id_col = "Individual.Id",
    time_col = "Time.to.IA",
    transpose = "always",
    legacy_na = TRUE
  )
  
  expect_true(is.array(result))
  expect_equal(length(dim(result)), 3)
})

# Test per documentare il comportamento attuale di .build_tensor
test_that("prepare_omics documents current .build_tensor limitations", {
  test_data <- data.frame(
    Individual.Id = c(1, 2, 3),
    Feature1 = c(1.1, 2.1, 3.1)
  )
  
  # Con tensor = TRUE ma senza time_col, dovrebbe restituire NULL
  result_null <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id",
    tensor = TRUE
  )
  expect_null(result_null)
  
  # Con tensor = FALSE, dovrebbe restituire data.frame processato
  result_df <- prepare_omics(
    data = test_data,
    id_col = "Individual.Id", 
    tensor = FALSE
  )
  expect_true(is.data.frame(result_df))
})

test_that("numeric_coercion non converte id_col", {
  df <- data.frame(Individual.Id=c("ID001","ID002"),
                   Time.to.IA=c(0,0), F1=c("1","x"))
  out <- prepare_omics(df, "Individual.Id", "Time.to.IA",
                       tensor=FALSE, numeric_coercion=TRUE)
  expect_type(out$Individual.Id, "character")
  expect_true(is.na(out$F1[2]))
})