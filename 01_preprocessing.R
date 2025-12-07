###############################################################
## 01_preprocessing.R
## Preprocessing pipeline for female and male astrocytes
###############################################################

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

female_path <- "data/Astrocytes_73F.rds"
male_path   <- "data/Astrocytes_70M.rds"

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# Load and update objects
astro_F <- UpdateSeuratObject(readRDS(female_path))
astro_M <- UpdateSeuratObject(readRDS(male_path))

# QC Metrics Function
qc_add <- function(obj) {
  assay <- DefaultAssay(obj)
  mat <- GetAssayData(obj, assay = assay, slot = "counts")
  obj$nFeature_RNA <- Matrix::colSums(mat > 0)
  obj$nCount_RNA <- Matrix::colSums(mat)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = assay)
  return(obj)
}

# QC Filter Function
qc_filter <- function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
}

# Seurat Workflow
preprocess_seurat <- function(obj, npcs = 30, res = 0.6) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = npcs)
  obj <- RunUMAP(obj, dims = 1:20)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = res)
  return(obj)
}

# Wrapper
process_one_group <- function(obj, group_label) {
  message("Processing group: ", group_label)
  obj <- qc_add(obj)
  obj <- qc_filter(obj)
  obj <- preprocess_seurat(obj)
  
  p <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0("Astrocytes (", group_label, ")"))
  ggsave(file.path("results/plots", paste0("UMAP_astro_", group_label, ".pdf")), p, width = 6, height = 5)
  return(obj)
}

# Execution
astro_F_clean <- process_one_group(astro_F, "F")
astro_M_clean <- process_one_group(astro_M, "M")

saveRDS(astro_F_clean, "data/astro_F_clean.rds")
saveRDS(astro_M_clean, "data/astro_M_clean.rds")
message("Preprocessing complete.")