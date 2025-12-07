###############################################################
## 02_nmf_stability_female.R
## NMF rank stability analysis for Female astrocytes
###############################################################

rm(list = ls())
set.seed(123)
source("scripts/00_setup_utils.R")

female_rds   <- "data/astro_F_clean.rds"
out_dir      <- "output/nmf_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
k_values     <- c(4, 6, 8, 10, 12)
n_boot       <- 30
donor_frac   <- 0.8
n_threads    <- 8L

# Load data
astro_F <- readRDS(female_rds)
cm <- get_counts_and_meta(astro_F)

# Run stability analysis
res_F <- nmf_stability_grid(
  counts     = cm$counts,
  meta       = cm$meta,
  donor_col  = "projid",
  k_values   = k_values,
  n_boot     = n_boot,
  donor_frac = donor_frac,
  n_threads  = n_threads
)

# Save results
write.csv(res_F, file.path(out_dir, "Female_NMF_Stability.csv"), row.names = FALSE)
saveRDS(res_F, file.path(out_dir, "Female_NMF_Stability.rds"))

# Plot
p_F <- plot_k_vs_stability(res_F, title = "Female astrocytes")
ggsave(file.path(out_dir, "Female_NMF_Stability_k_plot.pdf"), p_F, width = 5, height = 4)