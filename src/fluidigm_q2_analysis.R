#Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

#Data
dat <- readr::read_tsv("data/FluidigmPlate.tab", show_col_types = FALSE)

gene_cols <- grep("^Gene", names(dat), value = TRUE)
meta <- dat %>% dplyr::select(sample, group)
X <- as.matrix(dat[, gene_cols])

##PCA

pca <- prcomp(X, center = TRUE, scale. = FALSE)

#Percent variance explained
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- as.data.frame(pca$x[, 1:3])
colnames(scores) <- c("PC1", "PC2", "PC3")
plot_df <- dplyr::bind_cols(meta, scores)

#PCA: PC1 vs PC2
p12 <- ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  labs(
    title = "PCA of 96-gene expression across 96 biopsies",
    x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
    color = "Group"
  ) +
  theme_minimal(base_size = 12)

#PCA: PC2 vs PC3
p23 <- ggplot(plot_df, aes(x = PC2, y = PC3, color = group)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(level = 0.95, linewidth = 0.8) +
  labs(
    title = "PCA of 96-gene expression (PC2 vs PC3)",
    x = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
    y = sprintf("PC3 (%.1f%%)", 100 * var_expl[3]),
    color = "Group"
  ) +
  theme_minimal(base_size = 12)

##K-means clustering on PCs

set.seed(42)
km <- kmeans(scores[, 1:3], centers = 3, nstart = 50)

plot_df$cluster <- factor(km$cluster)

p_km <- ggplot(plot_df, aes(x = PC1, y = PC2, color = cluster, shape = group)) +
  geom_point(size = 2.2, alpha = 0.9) +
  stat_ellipse(aes(color = cluster), level = 0.95, linewidth = 0.8) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "K-means clusters (k = 3) overlaid on PCA",
    x = sprintf("PC1 (%.1f%%)", 100 * var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_expl[2]),
    color = "Cluster",
    shape = "True group"
  ) +
  theme_minimal(base_size = 12)

##Mean expression per gene and group

df_long <- cbind(group = meta$group, as.data.frame(X)) %>%
  as_tibble() %>%
  pivot_longer(
    cols = starts_with("Gene"),
    names_to = "Gene",
    values_to = "Expression"
  )

gene_means <- df_long %>%
  group_by(Gene, group) %>%
  summarise(mean_expression = mean(Expression), .groups = "drop")

##Correlations of group means

gene_means_wide <- gene_means %>%
  pivot_wider(names_from = group, values_from = mean_expression)

plot_correlation <- function(df, group_x, group_y) {
  cor_test <- cor.test(df[[group_x]], df[[group_y]], method = "pearson")
  r_val <- cor_test$estimate
  p_val <- cor_test$p.value
  
  ggplot(df, aes_string(x = group_x, y = group_y)) +
    geom_point(color = "steelblue", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(
      title = sprintf("Mean expression: %s vs %s", group_x, group_y),
      subtitle = sprintf("Pearson r = %.3f, p = %.2e", r_val, p_val),
      x = sprintf("%s mean expression", group_x),
      y = sprintf("%s mean expression", group_y)
    ) +
    theme_minimal(base_size = 12)
}

p_uc_cd  <- plot_correlation(gene_means_wide, "UC",  "CD")
p_cd_con <- plot_correlation(gene_means_wide, "CD",  "CON")
p_uc_con <- plot_correlation(gene_means_wide, "UC",  "CON")

##Top-10 genes by |log2 fold change|

compute_log2fc <- function(df, g1, g2, prefix) {
  df %>%
    mutate(!!paste0("log2FC_", prefix) := log2((.data[[g1]] + 1e-6) / (.data[[g2]] + 1e-6)))
}

gene_means_fc <- gene_means_wide %>%
  compute_log2fc("CD",  "CON", "CD_CON") %>%
  compute_log2fc("CD",  "UC",  "CD_UC")  %>%
  compute_log2fc("UC",  "CON", "UC_CON")

top_cd_con <- gene_means_fc %>%
  arrange(desc(abs(log2FC_CD_CON))) %>%
  slice_head(n = 10) %>%
  select(Gene, log2FC_CD_CON)

top_cd_uc <- gene_means_fc %>%
  arrange(desc(abs(log2FC_CD_UC))) %>%
  slice_head(n = 10) %>%
  select(Gene, log2FC_CD_UC)

top_uc_con <- gene_means_fc %>%
  arrange(desc(abs(log2FC_UC_CON))) %>%
  slice_head(n = 10) %>%
  select(Gene, log2FC_UC_CON)

##Expression distributions for top genes

genes_cd_con <- top_cd_con %>% arrange(desc(abs(log2FC_CD_CON))) %>% pull(Gene)
genes_cd_uc  <- top_cd_uc  %>% arrange(desc(abs(log2FC_CD_UC)))  %>% pull(Gene)
genes_uc_con <- top_uc_con %>% arrange(desc(abs(log2FC_UC_CON))) %>% pull(Gene)

make_dist_plot <- function(genes_vec, title_text){
  df_plot <- df_long %>%
    filter(Gene %in% genes_vec) %>%
    mutate(Gene = factor(Gene, levels = genes_vec))
  
  ggplot(df_plot, aes(x = Gene, y = Expression, color = group,
                      group = interaction(group, Gene))) +
    geom_boxplot(position = position_dodge(width = 0.75),
                 outlier.shape = NA, width = 0.65) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width = 0.75),
                alpha = 0.6, size = 1.6) +
    labs(title = title_text,
         x = "Gene (top 10 by |log2FC|)",
         y = "Expression") +
    theme_minimal(base_size = 12)
}

p_dist_cd_con <- make_dist_plot(genes_cd_con, "CD vs CON: distribution of top-10 genes")
p_dist_cd_uc  <- make_dist_plot(genes_cd_uc,  "CD vs UC: distribution of top-10 genes")
p_dist_uc_con <- make_dist_plot(genes_uc_con, "UC vs CON: distribution of top-10 genes")

##Plots

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

ggsave("results/pca_pc1_pc2.png",      p12,             width = 7, height = 5, dpi = 300)
ggsave("results/pca_pc2_pc3.png",      p23,             width = 7, height = 5, dpi = 300)
ggsave("results/kmeans_pca.png",       p_km,            width = 7, height = 5, dpi = 300)
ggsave("results/uc_vs_cd_corr.png",    p_uc_cd,         width = 7, height = 5, dpi = 300)
ggsave("results/cd_vs_con_corr.png",   p_cd_con,        width = 7, height = 5, dpi = 300)
ggsave("results/uc_vs_con_corr.png",   p_uc_con,        width = 7, height = 5, dpi = 300)
ggsave("results/topgenes_cd_con.png",  p_dist_cd_con,   width = 7, height = 5, dpi = 300)
ggsave("results/topgenes_cd_uc.png",   p_dist_cd_uc,    width = 7, height = 5, dpi = 300)
ggsave("results/topgenes_uc_con.png",  p_dist_uc_con,   width = 7, height = 5, dpi = 300)
