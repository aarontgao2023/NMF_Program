###############################################################
## 06_plot_enrichment.R
## Visualization of GO enrichment results
###############################################################

rm(list = ls())
set.seed(123)
library(dplyr); library(readr); library(ggplot2); library(forcats)

enrich_dir  <- "output/topic_enrichment"
fig_dir     <- "figs/enrichment"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

topic_match <- read_csv("output/cross_sex/CrossSex_topic_matching_k6.csv", show_col_types = FALSE)

# Helper: Categorize GO terms
assign_pathway_category <- function(description) {
  d <- tolower(description)
  case_when(
    grepl("synap|dendrit|axon|neuron", d) ~ "Synaptic / Structural",
    grepl("membrane potential|ion channel", d) ~ "Electrophysiology",
    grepl("calcium", d) ~ "Calcium handling",
    grepl("junction|adhesion", d) ~ "Adhesion",
    grepl("immune|inflam", d) ~ "Immune",
    grepl("lipid|metabolism", d) ~ "Metabolism",
    TRUE ~ "Other"
  )
}

# Helper: Plot Bar
plot_go_bar <- function(csv_file, title) {
  if (!file.exists(csv_file)) return(NULL)
  df <- read_csv(csv_file, show_col_types = FALSE) %>%
    arrange(p.adjust) %>% slice_head(n = 12) %>%
    mutate(Category = assign_pathway_category(Description),
           Description = fct_reorder(Description, -log10(p.adjust)))
  
  ggplot(df, aes(x = -log10(p.adjust), y = Description, fill = Category)) +
    geom_col() + theme_bw() + labs(title = title, x = "-log10(adj P)")
}

# 1. Plot Matched Pairs
for (i in seq_len(nrow(topic_match))) {
  f_tp <- topic_match$Female_topic[i]
  m_tp <- topic_match$Male_topic[i]
  
  p1 <- plot_go_bar(file.path(enrich_dir, paste0("Female_", f_tp, "_GO_BP.csv")), f_tp)
  if(!is.null(p1)) ggsave(file.path(fig_dir, paste0("Bar_", f_tp, ".pdf")), p1, width = 6, height = 4)
}

# 2. Plot Female_t6
p_f6 <- plot_go_bar(file.path(enrich_dir, "Female_Female_t6_GO_BP.csv"), "Female_t6")
if(!is.null(p_f6)) ggsave(file.path(fig_dir, "Bar_Female_t6.pdf"), p_f6, width = 7, height = 5)