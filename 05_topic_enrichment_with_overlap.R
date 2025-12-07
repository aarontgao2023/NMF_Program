###############################################################
## 05_topic_enrichment_with_overlap.R
## Topic gene ranking, overlap calculation, and GO enrichment
###############################################################

rm(list = ls())
set.seed(123)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(clusterProfiler); library(org.Hs.eg.db); library(tibble)
})

# Parameters
out_dir     <- "output/topic_enrichment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
top_n_genes <- 100

nmfF <- readRDS("output/cross_sex/Female_NMF_k6.rds")
nmfM <- readRDS("output/cross_sex/Male_NMF_k6.rds")
topic_match <- read_csv("output/cross_sex/CrossSex_topic_matching_k6.csv", show_col_types = FALSE)

# Helper: Rank genes
rank_genes_by_loading <- function(W, sex_label) {
  as.data.frame(W) %>%
    mutate(gene = rownames(W)) %>%
    pivot_longer(-gene, names_to = "topic", values_to = "loading") %>%
    group_by(topic) %>%
    arrange(desc(loading)) %>%
    mutate(rank = row_number(), sex = sex_label) %>%
    ungroup()
}

genes_F <- rank_genes_by_loading(nmfF$W, "F")
genes_M <- rank_genes_by_loading(nmfM$W, "M")

get_top_genes <- function(df, topic_name, n = 100) {
  df %>% filter(topic == topic_name) %>% arrange(rank) %>% slice_head(n = n) %>% pull(gene)
}

# 1. Matched Pair Overlap
overlap_list <- lapply(seq_len(nrow(topic_match)), function(i) {
  f_tp <- topic_match$Female_topic[i]
  m_tp <- topic_match$Male_topic[i]
  g_f <- get_top_genes(genes_F, f_tp, top_n_genes)
  g_m <- get_top_genes(genes_M, m_tp, top_n_genes)
  
  tibble(
    Female_Topic = f_tp, Male_Topic = m_tp,
    Correlation = topic_match$similarity[i],
    Jaccard_Index = length(intersect(g_f, g_m)) / length(union(g_f, g_m)),
    Intersection_Count = length(intersect(g_f, g_m))
  )
})
write_csv(bind_rows(overlap_list), file.path(out_dir, "Matched_Topic_Gene_Overlap_Stats.csv"))

# 2. Within-Sex Overlap Calculation
calc_within_sex_overlap <- function(genes_df, sex_label, top_n) {
  topics <- unique(genes_df$topic)
  n_t <- length(topics)
  mat_jaccard <- matrix(0, n_t, n_t, dimnames = list(topics, topics))
  
  for (i in 1:n_t) {
    for (j in 1:n_t) {
      g1 <- get_top_genes(genes_df, topics[i], top_n)
      g2 <- get_top_genes(genes_df, topics[j], top_n)
      mat_jaccard[topics[i], topics[j]] <- length(intersect(g1, g2)) / length(union(g1, g2))
    }
  }
  return(mat_jaccard)
}

write.csv(calc_within_sex_overlap(genes_F, "F", top_n_genes), file.path(out_dir, "Female_Within_Overlap_Jaccard.csv"))
write.csv(calc_within_sex_overlap(genes_M, "M", top_n_genes), file.path(out_dir, "Male_Within_Overlap_Jaccard.csv"))

# 3. GO Enrichment
run_enrich <- function(genes, background) {
  enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = background)
}

background <- union(rownames(nmfF$W), rownames(nmfM$W))

# Enrich all topics and save
process_enrichment <- function(topic_list, genes_df, prefix) {
  for (tp in topic_list) {
    ego <- run_enrich(get_top_genes(genes_df, tp, top_n_genes), background)
    if (!is.null(ego) && nrow(ego) > 0) {
      write_csv(as.data.frame(ego), file.path(out_dir, paste0(prefix, "_", tp, "_GO_BP.csv")))
    }
  }
}

process_enrichment(unique(genes_F$topic), genes_F, "Female")
process_enrichment(unique(genes_M$topic), genes_M, "Male")