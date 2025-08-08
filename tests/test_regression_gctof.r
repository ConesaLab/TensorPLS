# Test di regressione con dati GCTOF reali

test_that("GCTOF data processing regression test", {
  # Skip se i file non esistono
  raw_path <- system.file("extdata", "GCTOF_Data_Processed.csv",
                         package = "NPLS.MultiomiX")
  cohort_path <- system.file("extdata", "CohortData.csv",
                            package = "NPLS.MultiomiX")
  
  skip_if(raw_path == "", "GCTOF test data not available")
  skip_if(cohort_path == "", "Cohort test data not available")
  
  # Leggi dati per verifica
  raw_df <- read.csv(raw_path, check.names = FALSE)
  tmp <- t(raw_df[-1, ])
  colnames(tmp) <- tmp[1, ]
  feature_cols_vec <- colnames(tmp)[4:367]
  
  CohortData <- read.csv(cohort_path, header = TRUE)
  subjects_vec <- CohortData$Group.Id[CohortData$Model.or.Validation == "Model"]
  
  # Test della tua funzione
  arr <- prepare_omics(
    data          = raw_path,
    id_col        = "Individual.Id",
    time_col      = "Time.to.IA",
    transpose     = "always",
    legacy_na     = TRUE,
    subjects      = subjects_vec,
    feature_cols  = feature_cols_vec,
    cohort        = cohort_path,
    cohort_id_col = "Group.Id",
    cohort_filter = "Model.or.Validation=='Model'"
  )
  
  # Verifica dimensioni
  expect_equal(dim(arr)[1], length(subjects_vec))
  expect_equal(dim(arr)[2], length(feature_cols_vec))
  expect_true(is.array(arr))
  expect_equal(length(dim(arr)), 3)
  
  # Verifica nomi
  expect_equal(dimnames(arr)$Feature, feature_cols_vec)
  expect_true(all(dimnames(arr)$Subject %in% as.character(subjects_vec)))
  
  # Verifica sanità dati
  expect_true(all(is.finite(arr) | is.na(arr)))
  
  # Salva reference la prima volta (decommentare per creare reference)
  # ref_dir <- "tests/testdata"
  # if (!dir.exists(ref_dir)) dir.create(ref_dir, recursive = TRUE)
  # saveRDS(arr, file.path(ref_dir, "gctof_reference.rds"))
  
  # Confronta con reference (se esiste)
  ref_file <- system.file("testdata", "gctof_reference.rds", 
                         package = "NPLS.MultiomiX")
  if (file.exists(ref_file)) {
    reference <- readRDS(ref_file)
    expect_equal(arr, reference, tolerance = 1e-10)
  }
})