setwd("/Users/alessandrogiordano/Desktop/AnaConesa/NPLSDAPackage/tests/manual/")


raw_path    <- system.file("extdata", "GCTOF_Data_Processed.csv",
                           package = "NPLS-MultiomiX")
cohort_path <- system.file("extdata", "CohortData.csv",
                           package = "NPLS-MultiomiX")

#I have to extract colnames from dataframe
raw_df         <- read.csv(raw_path, check.names = FALSE)
head(raw_df)
tmp            <- t(raw_df[-1, ])
head(tmp)
colnames(tmp)  <- tmp[1, ]
colnames(tmp)
feature_cols_vec <- colnames(tmp)[4:367]

feature_cols_vec

CohortData   <- read.csv(cohort_path, header = TRUE)

subjects_vec <- CohortData$Group.Id[CohortData$Model.or.Validation == "Model"]

subjects_vec

arr <- prepare_omics(
  data          = raw_path,
  id_col        = "Individual.Id",
  time_col      = "Time.to.IA",
  transpose     = "always",
  legacy_na     = FALSE,
  subjects      = subjects_vec,
  feature_cols  = feature_cols_vec,
  cohort        = cohort_path,
  cohort_id_col = "Group.Id",
  cohort_filter = "Model.or.Validation=='Model'"
)

head(arr[,3,5])

head(arr)
dim(arr)
dimnames(arr)
dimnames(arr)$Subject
identical(arr,arr1)
arr1 <- prepare_omics(
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
identical(arr1,arr)
dim(arr)
arr[,,1]
head(arr)
dimnames(arr)$Subject
summary(arr)
identical(arr,)
load("/Users/alessandrogiordano/Desktop/TEDDY/Metabolomics/RData/fullarrayGCTOF136Proc.RData")
load("/Users/alessandrogiordano/Desktop/TEDDY/Metabolomics/RData/arrayGCTOFXprocessed.RData")
identical(arrayGCTOFX,arr)
head(arr)
summary(arr)
summary(arrayGCTOFX)
fullarrayGCTOF136
arr

arr[,1,1]


raw_df <- read.csv("/Users/alessandrogiordano/Desktop/TEDDY/GeneExpression/GeneExpressionDataProcessed.csv", check.names = FALSE)

tmp    <- t(raw_df[-1, ])
colnames(tmp) <- tmp[1, ]
tmp
feature_cols_vec <- colnames(tmp)[4:length(colnames(tmp))]
length(feature_cols_vec)
CohortData<-read.csv ("/Users/alessandrogiordano/Desktop/TEDDY/SupplementaryData/GeneExpression/CohortData.csv",header = TRUE)
subjects_vec <- CohortData$Group.Id[CohortData$Model.or.Validation == "Model"]
subjects_vec
length(subjects_vec)









library(tibble)
library(tidyr)
library(dplyr)
counts = read.csv("/Users/alessandrogiordano/Downloads/RNA_seq/STATegra.RNAseq.allSamples.counts.csv")
# 1) counts: righe=gene (rownames), colonne=sample (colnames)
genes   <- counts$GeneName
genes
samples <- setdiff(colnames(counts), "GeneName")
samples
# 2) Metadati dai nomi colonna: ExpBatch_<rep>_<cond>_<H>H
meta <- tibble(sample_id = samples) %>%
  extract(sample_id, into = c("replicate_id","condition","time_h"),
          regex = "^ExpBatch_(\\d+)_(Ctr|Ik)_(\\d+)H$", remove = FALSE) %>%
  mutate(
    replicate_id = as.integer(replicate_id),
    time         = as.integer(time_h),
    subject_id   = paste0(replicate_id, "_", condition),
    condition    = factor(condition, levels = c("Ctr","Ik"))
  )
counts_mat <- counts %>% 
  column_to_rownames(var = "GeneName")  # ora rownames(counts_mat) = ENSMUSG...
# Verifica:
stopifnot(all(feature_cols %in% rownames(counts_mat)))# --- transpose counts in "samples × genes" conservando i nomi dei geni ---
sxg <- as.data.frame(t(counts_mat), check.names = FALSE)
dim(sxg)    

names(sxg)[1:5]

sxg$sample_id <- rownames(sxg)

# 4) Unisci con i tuoi metadati
long_df <- meta %>% 
  left_join(sxg, by = "sample_id")

# controllo: i nomi colonna devono includere i geneID
stopifnot(all(feature_cols %in% names(long_df)))
# --- unione metadati + features ---


# 7) Costruisci il tensor Subject × Feature × Time con le TUE funzioni
arr <- prepare_omics(
  data            = long_df %>% select(subject_id, time, all_of(feature_cols)),
  id_col          = "subject_id",
  time_col        = "time",
  tensor          = TRUE,
  transpose       = "never",
  numeric_coercion= FALSE,   # già numerico, ma fa comodo
  legacy_na       = FALSE,
  min_timepoints  = 2,      # coerente con il filtro
  feature_cols    = feature_cols
)
dim(arr)
arr[,1,1]
dim(arr)



## ====== 0) Utility comuni ======
subj <- dimnames(arr)$Subject
stopifnot(length(subj) == 6)

## Etichette: 0=Ctr, 1=Ik
Y_vec <- ifelse(grepl("_Ik$", subj), 1L, 0L)
Y_mat <- matrix(Y_vec, ncol = 1,
                dimnames = list(subj, "class"))
Y_mat
## Matricizza il 3D (n × p*k)
flatten_tensor <- function(X) {
  n <- dim(X)[1]; p <- dim(X)[2]; k <- dim(X)[3]
  matrix(aperm(X, c(1,2,3)), nrow = n, ncol = p*k,
         dimnames = list(dimnames(X)$Subject, as.vector(outer(dimnames(X)$Feature, dimnames(X)$Time, paste, sep="_"))))
}

## Pieghe leave-pair-out per replicato: (1_Ctr,1_Ik) | (2_Ctr,2_Ik) | (3_Ctr,3_Ik)
make_pairs_folds <- function(subjects) {
  rep_id <- sub("^([0-9]+)_.*$", "\\1", subjects)
  lapply(sort(unique(rep_id)), function(r) which(rep_id == r))
}

## Dato un fit di compute_npls_factors, costruisci la matrice K (colonne = kron(Wk_f, Wj_f))
build_K <- function(Wj, Wk, ncomp) {
  do.call(cbind, lapply(seq_len(ncomp), function(f) kronecker(Wk[, f], Wj[, f])))
}

## Proietta X (matricizzato) sui fattori (scores Mode1)
project_scores <- function(Xmat, Wj, Wk, ncomp) {
  K <- build_K(Wj, Wk, ncomp)
  Xmat %*% K
}
## ====== 1) Valutazione CV con Factors (B) ======
eval_cv_factors <- function(X, Y, ncomp = 2, center_opt = 2) {
  subjects <- dimnames(X)$Subject
  folds <- make_pairs_folds(subjects)
  acc <- numeric(length(folds))
  
  for (i in seq_along(folds)) {
    te_idx <- folds[[i]]
    tr_idx <- setdiff(seq_along(subjects), te_idx)
    
    Xtr <- X[tr_idx, , , drop = FALSE]
    Xte <- X[te_idx, , , drop = FALSE]
    Ytr <- matrix(Y[tr_idx], ncol = 1, dimnames = list(dimnames(X)$Subject[tr_idx], "class"))
    Yte <- matrix(Y[te_idx], ncol = 1)
    
    ## 1) Fattori NPLS dal TRAIN
    fit_tr <- compute_npls_factors(Xtr, Ytr, ncomp = ncomp, center_opt = center_opt)
    
    ## 2) Scores TRAIN (già disponibili) e TEST (proiezione)
    T_train <- fit_tr$FactorsX$Mode1                   # n_tr × ncomp
    Xte_mat <- flatten_tensor(Xte)
    Wj <- fit_tr$FactorsX$Mode2                         # p × ncomp
    Wk <- fit_tr$FactorsX$Mode3                         # k × ncomp
    T_test  <- project_scores(Xte_mat, Wj, Wk, ncomp)   # n_te × ncomp
    
    ## 3) Classificatore semplice sui T: glm binario
    df_tr <- data.frame(y = as.factor(drop(Ytr)), T_train)
    mdl   <- glm(y ~ ., data = df_tr, family = binomial())
    pr    <- predict(mdl, newdata = data.frame(T_test), type = "response")
    ypred <- ifelse(pr >= 0.5, 1L, 0L)
    acc[i] <- mean(ypred == drop(Yte))
  }
  list(mean_acc = mean(acc), acc_per_fold = acc)
}

res_B <- eval_cv_factors(arr, Y_vec, ncomp = 2, center_opt = 2)
res_B

res_C <- eval_cv_factors(arr_C, Y_vec, ncomp = 2, center_opt = 0)
res_C
#' Compute NPLS Factors with Core Arrays
#' 
#' Calcola i fattori NPLS e i core arrays G e Gu seguendo  
#' l'algoritmo di NPLSDAmod originale.
#'
#' @param X Array 3D (n × p × k) dove n = individui, p = variabili, k = tempo
#' @param Y Matrice (n × q) per classic NPLS-DA o Array 3D (n × q × k) per regression class 2
#' @param ncomp Numero di componenti da estrarre
#' @param center_opt Opzione di centratura: 0=no, 1=Mode1, 2=Mode2, 3=Mode3. Default: 2
#' @param tol Tolleranza per convergenza algoritmo iterativo. Default: 1e-16
#' @param max_iter Numero massimo di iterazioni per convergenza. Default: 10000
#' 
#' @return Lista contenente:
#'   \item{FactorsX}{Lista con Mode1 (scores), Mode2 (loadings variabili), Mode3 (loadings tempo)}
#'   \item{FactorsY}{Lista con fattori Y (solo se Y è 3D)}
#'   \item{B}{Matrice coefficienti regressione}
#'   \item{G}{Core array 3D (ncomp × ncomp × ncomp)}
#'   \item{Gu}{Lista di core arrays per ogni componente}
#' 
#' @export
compute_npls_factors <- function(X, Y, ncomp, 
                                 center_opt = 2,
                                 tol = 1e-16, 
                                 max_iter = 10000) {
  
  cat("Calcolo Factors NPLS con core arrays (logica esatta NPLSDAmod) in corso...\n")
  
  # Stesso pre-processing della funzione principale
  if (center_opt != 0) {
    X <- center_mode(X, mode = center_opt)
    if (length(dim(Y)) == 3)
      Y <- center_mode(Y, mode = center_opt)
    else
      Y <- scale(Y, center = TRUE, scale = FALSE)
  }
  
  # Calcola Factors E core arrays 
  npls_result <- .npls_factors_with_core_arrays(X, Y, ncomp, tol, max_iter)
  
  Y_is_3way <- (length(dim(Y)) == 3L)
  
  FactorsX <- list(
    Mode1 = npls_result$Tt,
    Mode2 = npls_result$WsupraJ,
    Mode3 = npls_result$WsupraK
  )
  
  FactorsY <- NULL
  if (Y_is_3way) {
    FactorsY <- list(
      Mode1 = npls_result$U,
      Mode2 = npls_result$QsupraJ,
      Mode3 = npls_result$QsupraK
    )
  }
  
  cat("Factors NPLS e core arrays calcolati con successo!\n")
  
  return(list(
    FactorsX = FactorsX, 
    FactorsY = FactorsY,
    B = npls_result$B,
    G = npls_result$G,        # Core array 3D
    Gu = npls_result$Gu       # Lista core arrays
  ))
}

#' NPLS Core Algorithm con G e Gu (Internal Function)
#' 
#' Implementa l'algoritmo core NPLS includendo calcolo di G e Gu durante 
#' il loop principale, esattamente come NPLSDAmod.
#'
#' @param X Array 3D centrato
#' @param Y Matrice o Array 3D centrato  
#' @param ncomp Numero di componenti
#' @param tol Tolleranza convergenza
#' @param max_iter Iterazioni massime
#' 
#' @return Lista con matrici di loadings, scores, coefficienti e core arrays
#' 
#' @noRd
.npls_factors_with_core_arrays <- function(X, Y, ncomp, tol = 1e-16, max_iter = 10000) {
  
  n <- dim(X)[1]; p <- dim(X)[2]; k <- dim(X)[3]
  Y_is_3way <- (length(dim(Y)) == 3L)
  q <- if (Y_is_3way) dim(Y)[2] else ncol(Y)
  
  # Inizializzazione ESATTA come NPLSDAmod
  Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
  B <- G <- matrix(0, ncol = ncomp, nrow = ncomp)  # ← G inizializzato come NPLSDAmod
  Gu <- vector("list", ncomp)                      # ← Gu inizializzato come NPLSDAmod
  
  # Unfold
  Xmat <- matrix(aperm(X, c(1, 2, 3)), n, p * k)
  if (Y_is_3way) {
    Ymat <- matrix(aperm(Y, c(1, 2, 3)), n, q * k)
  } else {
    # Per Classic NPLS-DA, duplica Y (come originale)
    Y_ext <- array(0, dim = c(n, q, k))
    for (i in 1:k) Y_ext[, , i] <- as.matrix(Y)
    Ymat <- matrix(aperm(Y_ext, c(1, 2, 3)), n, q * k)
  }
  
  # Copie per deflazione (ESATTO come originale)
  X_work <- Xmat  # X NON viene deflazionato nel loop!
  Y_work <- Ymat  # Solo Y viene deflazionato
  
  # Loop principale 
  for (f in 1:ncomp) {
    
    # Inizializzazione u 
    if (ncol(Y_work) > 0) {
      Uf <- svd(Y_work)$u
      u <- Uf[, 1]
    } else {
      u <- rnorm(n)
    }
    
    # Iterazioni convergenza 
    for (it in 1:max_iter) {
      
      # STEP 1: Calcola Z e SVD (ESATTO)
      tX <- t(X_work)
      Zrow <- tX %*% u
      Z <- matrix(Zrow, nrow = p, ncol = k)
      svd_z <- svd(Z)
      wsupraj <- svd_z$u[, 1]
      wsuprak <- svd_z$v[, 1]
      
      # STEP 2: Calcola scores X (ESATTO)
      tf <- X_work %*% kronecker(wsuprak, wsupraj)
      
      # STEP 3: Calcola V e SVD (ESATTO)
      Vrow <- t(Y_work) %*% tf
      V <- matrix(Vrow, nrow = q, ncol = k)
      svd_v <- svd(V)
      qsupraj <- svd_v$u[, 1]
      qsuprak <- svd_v$v[, 1]
      
      # STEP 4: Calcola scores Y (ESATTO)
      uf <- Y_work %*% kronecker(qsuprak, qsupraj)
      
      # Test convergenza (ESATTO come originale)
      if (sum((uf - u)^2) < tol) {
        cat("Componente", f, "convergenza in", it, "iterazioni\n")
        break
      }
      u <- uf
    }
    
    # Salva risultati 
    Tt <- cbind(Tt, tf)
    WsupraJ <- cbind(WsupraJ, wsupraj)
    WsupraK <- cbind(WsupraK, wsuprak)
    U <- cbind(U, uf)
    QsupraJ <- cbind(QsupraJ, qsupraj)
    QsupraK <- cbind(QsupraK, qsuprak)
    
    # Calcola coefficiente regressione (ESATTO come originale)
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package MASS richiesto per ginv()")
    }
    
    bf <- MASS::ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
    B[1:length(bf), f] <- bf
    
    # ========================================================================
    # CALCOLA Gu per questa componente (ESATTO come NPLSDAmod)
    # ========================================================================
    
    # Pseudo-inverse matrices 
    TM <- MASS::ginv(t(Tt) %*% Tt) %*% t(Tt)
    WkM <- MASS::ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
    WjM <- MASS::ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
    
    # Core per questa componente 
    Gu[[f]] <- TM %*% X_work %*% kronecker(t(WkM), t(WjM))
    
    # Deflazione SOLO di Y 
    Y_work <- Y_work - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
    
    # Aggiorna u per prossima componente (come originale)
    if (f < ncomp) {
      Uf <- svd(Y_work)$u
      u <- Uf[, 1]
    }
  }
  
  # ========================================================================
  # CALCOLA G: Core array 3D 
  # ========================================================================
  
  #  usa la seconda componente per G, o la prima se ncomp=1
  if (ncomp >= 2) {
    G <- array(as.vector(Gu[[2]]), c(ncomp, ncomp, ncomp))
  } else {
    G <- array(as.vector(Gu[[1]]), c(ncomp, ncomp, ncomp))
  }
  
  # Nomi ESATTI come originale
  rownames(Tt) <- dimnames(X)[[1]]
  rownames(WsupraJ) <- dimnames(X)[[2]]
  rownames(WsupraK) <- dimnames(X)[[3]]
  rownames(U) <- dimnames(X)[[1]]
  
  if (Y_is_3way) {
    rownames(QsupraJ) <- dimnames(Y)[[2]]
    rownames(QsupraK) <- dimnames(Y)[[3]]
  } else {
    rownames(QsupraJ) <- colnames(Y)
    rownames(QsupraK) <- dimnames(X)[[3]]
  }
  
  return(list(
    WsupraJ = WsupraJ, WsupraK = WsupraK,
    QsupraJ = QsupraJ, QsupraK = QsupraK,
    Tt = Tt, U = U, B = B,
    G = G,   # ← Core array 3D
    Gu = Gu  # ← Lista core arrays per componente
  ))
}

center_mode <- function(X, mode = 2) {
  # X  : array 3-D
  # mode = 0 → nessuna centratura
  #        1 → centra per "soggetto"      (dimensione 1)
  #        2 → centra per "variabile"     (dimensione 2)
  #        3 → centra per "time-point"    (dimensione 3)
  
  stopifnot(length(dim(X)) == 3L, mode %in% 0:3)
  if (mode == 0) return(X)               # niente da fare
  
  d  <- dim(X)                            # c(I, J, K)
  dn <- dimnames(X)                       # salviamo le etichette
  
  ## 1. Portiamo l’asse da centrare al primo posto ---------------
  perm <- switch(as.character(mode),
                 "1" = c(1, 2, 3),        # soggetti già primo
                 "2" = c(2, 1, 3),        # variabili → primo
                 "3" = c(3, 1, 2))        # tempo → primo
  
  Xp  <- aperm(X, perm)                   # ruota l’array
  
  ## 2. Appiattiamo in matrice (righe = asse da centrare) --------
  mat <- matrix(Xp,
                nrow = d[perm[1]],        # lunghezza dell’asse scelto
                ncol = prod(d[-perm[1]]))# tutte le altre combinate
  
  ## 3. Trasformiamo eventuali fattori / char in numerico --------
  mat <- apply(mat, 2, as.numeric)
  
  ## 4. Sottraiamo la media di OGNI riga -------------------------
  mat <- mat - matrix(rowMeans(mat, na.rm = TRUE),
                      nrow(mat),
                      ncol(mat))
  
  ## 5. Torniamo alla forma 3-D originale ------------------------
  Xc  <- aperm(array(mat, d[perm]), order(perm))  # srotola e ruota indietro
  dimnames(Xc) <- dn                              # rimettiamo le etichette
  Xc
}
