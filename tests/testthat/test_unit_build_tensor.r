# Test unitari per .build_tensor() - Array3D

test_that(".build_tensor creates 3D array when time_col is present", {
  test_data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(-6, -3, -6, -3),
    feat = c(1.1, 1.2, 2.1, 2.2)
  )
  
  result <- .build_tensor(test_data, id_col = "id", time_col = "time")
  
  expect_equal(length(dim(result)), 3)
  expect_equal(dim(result), c(2, 1, 2))  # 2 soggetti, 1 feature, 2 tempi
})

test_that(".build_tensor returns NULL when time_col is NULL", {
  test_data <- data.frame(
    id = c(1, 2, 3),
    feat1 = c(1.1, 2.1, 3.1),
    feat2 = c(1.2, 2.2, 3.2)
  )
  
  result <- .build_tensor(test_data, id_col = "id", time_col = NULL)
  
  # Attualmente .build_tensor non gestisce il caso 2D, quindi dovrebbe restituire NULL
  expect_null(result)
})

test_that(".build_tensor returns NULL when time_col not in dataframe", {
  test_data <- data.frame(
    id = c(1, 2, 3),
    feat1 = c(1.1, 2.1, 3.1),
    feat2 = c(1.2, 2.2, 3.2)
  )
  
  # time_col specificato ma non esiste nel dataframe
  result <- .build_tensor(test_data, id_col = "id", time_col = "nonexistent")
  
  # Dovrebbe restituire NULL perché la colonna non esiste
  expect_null(result)
})

test_that(".build_tensor passes arguments correctly to .create_array3d", {
  test_data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(-6, -3, -6, -3),
    feat1 = c(1.1, 1.2, 2.1, 2.2),
    feat2 = c(3.1, 3.2, 4.1, 4.2),
    exclude_me = c("a", "b", "c", "d")
  )
  
  result <- .build_tensor(test_data,
                          id_col = "id",
                          time_col = "time", 
                          exclude_cols = c("id", "time", "exclude_me"),
                          subjects = c(2, 1),  # Ordine specifico
                          time_points = c(-3, -6))  # Ordine specifico
  
  # Verifica che l'ordine sia rispettato
  expect_equal(dimnames(result)$Subject, c("2", "1"))
  expect_equal(dimnames(result)$Time, c("-3", "-6"))
  expect_equal(dimnames(result)$Feature, c("feat1", "feat2"))
})

test_that(".build_tensor handles legacy_na parameter", {
  test_data <- data.frame(
    id = c(1, 1),
    time = c(-6, -3),
    feat = c(1.1, 1.2)
  )
  
  # Test con legacy_na = TRUE
  result_legacy <- .build_tensor(test_data, "id", "time", legacy_na = TRUE)
  expect_true(is.numeric(result_legacy))
  expect_equal(length(dim(result_legacy)), 3)
  
  # Test con legacy_na = FALSE
  result_modern <- .build_tensor(test_data, "id", "time", legacy_na = FALSE)
  expect_true(is.numeric(result_modern))
  expect_equal(length(dim(result_modern)), 3)
})

test_that(".build_tensor handles edge cases with time columns", {
  # Data.frame con una sola riga
  single_row <- data.frame(id = 1, time = -6, feat = 1.1)
  result_single <- .build_tensor(single_row, "id", "time")
  expect_equal(dim(result_single), c(1, 1, 1))
  
  # Data.frame con multiple features
  multi_feat <- data.frame(
    id = c(1, 1),
    time = c(-6, -3),
    feat1 = c(1.1, 1.2),
    feat2 = c(2.1, 2.2),
    feat3 = c(3.1, 3.2)
  )
  result_multi <- .build_tensor(multi_feat, "id", "time")
  expect_equal(dim(result_multi), c(1, 3, 2))  # 1 soggetto, 3 feature, 2 tempi
})

test_that(".build_tensor handles duplicates through .create_array3d", {
  # Test che i duplicati vengano gestiti correttamente
  test_data <- data.frame(
    id = c(1, 1, 1),  # Stesso soggetto, stesso tempo
    time = c(-6, -6, -6),
    feat = c(1.0, 1.1, 1.2)
  )
  
  result_first <- .build_tensor(test_data, "id", "time", deduplicate = "first")
  expect_equal(result_first[1, 1, 1], 1.0)
  
  result_mean <- .build_tensor(test_data, "id", "time", deduplicate = "mean")
  expect_equal(result_mean[1, 1, 1], 1.1)  # Media di 1.0, 1.1, 1.2
})