###############################################################
## 00_setup_utils.R
## Common utilities for NMF stability and cross-sex analysis
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppML)   
  library(clue)     
  library(dplyr)
  library(ggplot2)
  library(pracma)   
})

set.seed(123)
options(stringsAsFactors = FALSE)

## ---- Helper: Extract counts and metadata ----
get_counts_and_meta <- function(seurat_obj, assay = "RNA") {
  DefaultAssay(seurat_obj) <- assay
  counts <- GetAssayData(seurat_obj, slot = "counts")
  meta   <- seurat_obj@meta.data
  return(list(counts = counts, meta = meta))
}

## ---- Helper: Run RcppML NMF ----
run_nmf_rcppml <- function(mat, k, max_epochs = 200L, tol = 1e-4, 
                           n_threads = NULL, verbose = FALSE, seed = NULL) {
  if (!is.null(n_threads)) {
    RcppML::setRcppMLthreads(as.integer(n_threads))
  }
  
  fit <- RcppML::nmf(
    A       = mat,
    k       = k,
    tol     = tol,
    maxit   = as.integer(max_epochs),   
    verbose = verbose,
    L1      = c(0, 0),
    seed    = seed
  )
  return(list(W = fit$w, H = fit$h))
}

## ---- Helper: Topic similarity (Pearson correlation) ----
topic_similarity_cor <- function(W1, W2) {
  common_genes <- intersect(rownames(W1), rownames(W2))
  W1c <- W1[common_genes, , drop = FALSE]
  W2c <- W2[common_genes, , drop = FALSE]
  sim <- cor(as.matrix(W1c), as.matrix(W2c))
  return(sim)
}

## ---- Helper: Hungarian matching algorithm ----
hungarian_match <- function(sim_mat, maximize = TRUE) {
  assignment <- clue::solve_LSAP(sim_mat, maximum = maximize)
  matched_scores <- sim_mat[cbind(seq_len(nrow(sim_mat)), assignment)]
  return(list(assignment = assignment, scores = matched_scores))
}

## ---- Helper: Donor-level bootstrap stability ----
nmf_stability_one_k <- function(counts, meta, donor_col = "donor_id", k, 
                                n_boot = 30, donor_frac = 0.8, 
                                max_epochs = 200L, n_threads = 1L) {
  donor_ids <- unique(meta[[donor_col]])
  W_list <- list()
  
  for (b in seq_len(n_boot)) {
    sampled_donors <- sample(donor_ids, size = max(2, floor(length(donor_ids) * donor_frac)), replace = TRUE)
    cells_sub <- rownames(meta)[meta[[donor_col]] %in% sampled_donors]
    mat_sub   <- counts[, cells_sub, drop = FALSE]
    mat_sub   <- mat_sub[rowSums(mat_sub) > 0, , drop = FALSE]
    
    nmf_out <- run_nmf_rcppml(mat_sub, k = k, max_epochs = max_epochs, n_threads = n_threads)
    W_tmp <- nmf_out$W
    rownames(W_tmp) <- rownames(mat_sub)
    colnames(W_tmp) <- paste0("k", k, "_b", b, "_t", seq_len(k))
    W_list[[b]] <- W_tmp
  }
  
  all_scores <- c()
  if (length(W_list) > 1) {
    for (i in 1:(length(W_list) - 1)) {
      for (j in (i + 1):length(W_list)) {
        sim_mat <- topic_similarity_cor(W_list[[i]], W_list[[j]])
        match   <- hungarian_match(sim_mat, maximize = TRUE)
        all_scores <- c(all_scores, match$scores)
      }
    }
  }
  
  return(list(stability_mean = mean(all_scores, na.rm = TRUE), stability_sd = sd(all_scores, na.rm = TRUE)))
}

## ---- Helper: Stability grid search ----
nmf_stability_grid <- function(counts, meta, donor_col = "donor_id", k_values = c(4,6,8,10), 
                               n_boot = 30, donor_frac = 0.8, max_epochs = 200L, n_threads = 1L) {
  res <- data.frame()
  for (k in k_values) {
    message("Running stability for k = ", k)
    stab <- nmf_stability_one_k(counts, meta, donor_col, k, n_boot, donor_frac, max_epochs, n_threads)
    res <- bind_rows(res, data.frame(k = k, stability_mean = stab$stability_mean, stability_sd = stab$stability_sd))
  }
  return(res)
}

## ---- Helper: Run final NMF and save ----
run_final_nmf_and_save <- function(counts, k, prefix, out_dir = "output/cross_sex", 
                                   max_epochs = 400L, n_threads = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  mat <- counts[rowSums(counts) > 0, , drop = FALSE]
  
  nmf_out <- run_nmf_rcppml(mat, k = k, max_epochs = max_epochs, n_threads = n_threads)
  W <- nmf_out$W; H <- nmf_out$H
  
  rownames(W) <- rownames(mat)
  colnames(W) <- paste0(prefix, "_t", seq_len(k))
  rownames(H) <- colnames(W)
  colnames(H) <- colnames(mat)
  
  saveRDS(list(W = W, H = H), file = file.path(out_dir, paste0(prefix, "_NMF_k", k, ".rds")))
  return(list(W = W, H = H))
}

## ---- Helper: Project new data onto fixed basis (NNLS) ----
project_to_basis_pracma <- function(X_new, W_basis) {
  common_genes <- intersect(rownames(X_new), rownames(W_basis))
  Xc <- as.matrix(X_new[common_genes, , drop = FALSE])
  Wc <- as.matrix(W_basis[common_genes, , drop = FALSE])
  
  k <- ncol(Wc)
  n_cells <- ncol(Xc)
  H_new <- matrix(0, nrow = k, ncol = n_cells, dimnames = list(colnames(Wc), colnames(Xc)))
  
  for (i in seq_len(n_cells)) {
    H_new[, i] <- pracma::lsqnonneg(Wc, Xc[, i])$x
  }
  return(H_new)
}

## ---- Helper: Reconstruction Error ----
reconstruction_error <- function(X, W, H) {
  common_genes <- intersect(rownames(X), rownames(W))
  Xc <- X[common_genes, , drop = FALSE]
  Wc <- W[common_genes, , drop = FALSE]
  X_hat <- Wc %*% H
  return(norm(Xc - X_hat, type = "F")^2)
}