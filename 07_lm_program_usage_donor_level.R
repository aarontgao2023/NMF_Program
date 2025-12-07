###############################################################
## 07_lm_program_usage_donor_level.R
## Linear modeling of donor-level topic usage against covariates
###############################################################

rm(list = ls())
set.seed(123)
library(Seurat); library(dplyr); library(readr); library(broom); library(ggplot2); library(tidyr)

out_dir <- "output/lm_donor_level"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load Data
nmfF <- readRDS("output/cross_sex/Female_NMF_k6.rds")
nmfM <- readRDS("output/cross_sex/Male_NMF_k6.rds")
meta_donor <- read_csv("data/meta_multiomics_78F_71M.csv", show_col_types = FALSE)
astro_F <- readRDS("data/astro_F_clean.rds")
astro_M <- readRDS("data/astro_M_clean.rds")

# 2. Compute Mean Usage per Donor
compute_donor_usage <- function(H, meta, sex) {
  as.data.frame(t(H)) %>%
    mutate(cell = rownames(.)) %>%
    inner_join(meta %>% select(cell = rownames(.), projid), by = "cell") %>%
    mutate(sex = sex)
}

usage_F <- compute_donor_usage(nmfF$H, astro_F@meta.data, "F")
usage_M <- compute_donor_usage(nmfM$H, astro_M@meta.data, "M")

usage_donor <- bind_rows(usage_F, usage_M) %>%
  group_by(projid, sex) %>%
  summarise(across(starts_with(c("Female_", "Male_")), mean), .groups = "drop") %>%
  mutate(projid = as.character(projid)) %>%
  inner_join(meta_donor %>% mutate(projid = as.character(projid)), by = "projid")

# 3. Normalization (Rank-based Inverse Normal Transformation)
rankNorm <- function(x) {
  r <- rank(x, na.last = "keep")
  n <- sum(!is.na(x))
  qnorm((r - 0.5) / n)
}

data_norm <- usage_donor %>%
  mutate(across(starts_with(c("Female_", "Male_")), rankNorm))

# 4. Fit Linear Models
prog_cols <- grep("^Female_|^Male_", colnames(usage_donor), value = TRUE)
covariates <- c("age_death_limited", "pmi", "apoe4_pos")

results <- bind_rows(lapply(prog_cols, function(topic) {
  form <- as.formula(paste(topic, "~", paste(covariates, collapse = "+")))
  lm(form, data = data_norm) %>% tidy() %>% mutate(topic = topic)
}))

write_csv(results, file.path(out_dir, "donor_LM_results.csv"))

# 5. Visualization (Raw Usage vs Age)
plot_age_assoc <- function(topic, label) {
  ggplot(usage_donor, aes_string(x = "age_death_limited", y = topic)) +
    geom_point(alpha = 0.7) + geom_smooth(method = "lm", color = "black") +
    labs(title = label, x = "Age", y = "Mean Usage (Raw)") + theme_bw()
}

p1 <- plot_age_assoc("Male_t5", "Positive Ctrl: Male_t5")
p2 <- plot_age_assoc("Female_t6", "Female_t6")
ggsave("figs/Scatter_Age_Male_t5_Female_t6.pdf", gridExtra::grid.arrange(p1, p2, ncol = 2), width = 8, height = 4)