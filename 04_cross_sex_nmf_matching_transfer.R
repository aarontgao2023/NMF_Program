###############################################################
## 04_cross_sex_nmf_matching_transfer.R
## Cross-sex NMF topic matching and transferability assessment
###############################################################

rm(list = ls())
set.seed(123)
source("scripts/00_setup_utils.R")

# Configuration
k_use      <- 6       
out_dir    <- "output/cross_sex"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n_threads  <- max(1L, floor(parallel::detectCores() / 2L))

# Load data
astro_F <- readRDS("data/astro_F_clean.rds")
astro_M <- readRDS("data/astro_M_clean.rds")
cmF <- get_counts_and_meta(astro_F)
cmM <- get_counts_and_meta(astro_M)

# 1. Run Final NMF
message("Running final NMF (k=", k_use, ")...")
resF <- run_final_nmf_and_save(cmF$counts, k_use, "Female", out_dir, n_threads = n_threads)
resM <- run_final_nmf_and_save(cmM$counts, k_use, "Male", out_dir, n_threads = n_threads)

# 2. Topic Matching (Hungarian Algorithm)
sim_FM <- topic_similarity_cor(resF$W, resM$W)
match_FM <- hungarian_match(sim_FM, maximize = TRUE)

matched_pairs <- data.frame(
  Female_topic = colnames(resF$W),
  Male_topic   = colnames(resM$W)[match_FM$assignment],
  similarity   = match_FM$scores
)
write.csv(matched_pairs, file.path(out_dir, paste0("CrossSex_topic_matching_k", k_use, ".csv")), row.names = FALSE)

# 3. Transferability Metrics
message("Computing transferability metrics...")
H_M_fromF <- project_to_basis_pracma(cmM$counts, resF$W)
RE_M_self <- reconstruction_error(cmM$counts, resM$W, resM$H)
RE_M_fromF <- reconstruction_error(cmM$counts, resF$W, H_M_fromF)

H_F_fromM <- project_to_basis_pracma(cmF$counts, resM$W)
RE_F_self <- reconstruction_error(cmF$counts, resF$W, resF$H)
RE_F_fromM <- reconstruction_error(cmF$counts, resM$W, H_F_fromM)

transfer_df <- data.frame(
  direction      = c("Female->Male", "Male->Female"),
  transfer_ratio = c(RE_M_fromF / RE_M_self, RE_F_fromM / RE_F_self)
)
write.csv(transfer_df, file.path(out_dir, paste0("CrossSex_transferability_k", k_use, ".csv")), row.names = FALSE)

# 4. Specific Analysis: Female_t6 vs All Male Topics
female_topic <- "Female_t6"
sim_f6_allM <- sapply(colnames(resM$W), function(tp_m) {
  topic_similarity_cor(resF$W[, female_topic, drop = FALSE], resM$W[, tp_m, drop = FALSE])[1, 1]
})

sim_f6_df <- data.frame(Male_topic = colnames(resM$W), similarity = sim_f6_allM)

# Visualization: Dot plot
p_dot <- ggplot(sim_f6_df, aes(x = Male_topic, y = "Female_t6", size = similarity, color = similarity)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "steelblue") +
  scale_size_continuous(range = c(3, 12)) +
  theme_bw() +
  labs(title = "Female_t6 vs Male topics Similarity", x = "Male topics", y = NULL)

ggsave("output/sex_specific/Female_t6_vs_Male_topics_dotplot.pdf", p_dot, width = 5, height = 2)