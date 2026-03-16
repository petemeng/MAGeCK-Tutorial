source('analysis/_common.R')
suppressPackageStartupMessages({
  library(MAGeCKFlute)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
})

gdata <- ReadRRA('results/mageck_test.gene_summary.txt')
gene_raw <- read_tsv('results/mageck_test.gene_summary.txt', show_col_types = FALSE)
cat('RRA genes:', nrow(gdata), '\n')

rank_df <- as_tibble(gdata) %>% mutate(rank = row_number(), label = if_else(row_number() <= 8, id, NA_character_))
rank_plot <- ggplot(rank_df, aes(rank, Score)) +
  geom_line(color = screen_colors$essential, linewidth = 0.8) +
  geom_point(color = screen_colors$essential, size = 0.8) +
  geom_text_repel(aes(label = label), size = 3, fontface = 'italic', max.overlaps = 12) +
  labs(x = 'Gene rank', y = 'Flute score', title = 'MAGeCKFlute-style gene ranking') +
  theme_screen(10)
save_plot('results/figures/pub_flute_rankview.png', rank_plot, 10, 5.6)

essential_ids <- gene_raw %>%
  filter(`neg|fdr` < 0.05) %>%
  pull(id)
gene_map <- bitr(essential_ids, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
cat('mapped essential genes:', nrow(gene_map), '\n')

go_bp <- enrichGO(gene = gene_map$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', qvalueCutoff = 0.2, readable = TRUE)
cat('significant GO BP terms:', nrow(go_bp@result), '\n')
write_tsv(as_tibble(go_bp@result), 'results/flute_go_bp.tsv')

p_kegg <- dotplot(go_bp, showCategory = 10) +
  labs(title = 'Essential genes: GO BP enrichment') +
  theme_screen(10)
save_plot('results/figures/pub_flute_kegg.png', p_kegg, 8, 6.5)

mle_raw <- read_tsv('results/mageck_mle.gene_summary.txt', show_col_types = FALSE)
cat('MLE rows:', nrow(mle_raw), '\n')

square_df <- mle_raw %>%
  transmute(
    Gene,
    beta = `KBM7|beta`,
    neglogfdr = -log10(`KBM7|wald-fdr` + 1e-30),
    class = case_when(
      `KBM7|wald-fdr` < 0.05 & `KBM7|beta` < -0.5 ~ 'Essential',
      `KBM7|wald-fdr` < 0.05 & `KBM7|beta` > 0.5 ~ 'Enriched',
      TRUE ~ 'NS'
    ),
    label = if_else(Gene %in% c('AATF', 'BCL2', 'C1orf109', 'ANAPC4', 'BUB3'), Gene, NA_character_)
  )

p_square <- ggplot(square_df, aes(beta, neglogfdr, color = class)) +
  geom_point(size = 0.7, alpha = 0.65) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  geom_text_repel(aes(label = label), size = 3, fontface = 'italic', max.overlaps = 12) +
  scale_color_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = screen_colors$ns)) +
  labs(x = 'KBM7 beta', y = '-log10(wald FDR)', title = 'SquareView-style beta classification') +
  theme_screen()
save_plot('results/figures/pub_flute_squareview.png', p_square, 7, 6)

cn_data <- read_tsv('data/mle/cnv_data.txt', show_col_types = FALSE)
cn_merged <- gene_raw %>%
  transmute(id, neg.lfc = `neg|lfc`) %>%
  left_join(cn_data, by = c('id' = 'SYMBOL')) %>%
  mutate(log_copy_number = .[[ncol(.)]]) %>%
  dplyr::select(id, neg.lfc, log_copy_number) %>%
  filter(!is.na(log_copy_number))
cat('CN-merged genes:', nrow(cn_merged), '\n')
cat('CN vs LFC correlation:', round(cor(cn_merged$log_copy_number, cn_merged$neg.lfc), 3), '\n')

p_cn <- ggplot(cn_merged, aes(log_copy_number, neg.lfc)) +
  geom_point(size = 0.5, alpha = 0.35, color = 'grey55') +
  geom_smooth(method = 'lm', color = screen_colors$essential, linewidth = 0.8, se = TRUE) +
  annotate('text', x = max(cn_merged$log_copy_number), y = max(cn_merged$neg.lfc), hjust = 1, vjust = 1, label = sprintf('r = %.3f', cor(cn_merged$log_copy_number, cn_merged$neg.lfc)), size = 4) +
  labs(x = 'log2 copy number', y = 'negative-selection LFC', title = 'Copy number bias check') +
  theme_screen()
save_plot('results/figures/pub_cn_bias.png', p_cn, 7, 6)

depmap_demo <- read_tsv('refs/depmap_dependency_demo.tsv', show_col_types = FALSE)
heatmap_mat <- depmap_demo %>%
  tidyr::pivot_wider(names_from = cell_line, values_from = dependency) %>%
  arrange(gene_name)
mat <- as.matrix(heatmap_mat[, -1])
rownames(mat) <- heatmap_mat$gene_name
pheatmap(
  mat,
  color = colorRampPalette(c('#08306B', 'white', '#B2182B'))(60),
  clustering_method = 'complete',
  filename = 'results/figures/pub_depmap_heatmap.png',
  width = 7,
  height = 6
)
cat('DepMap-demo heatmap genes:', nrow(mat), '\n')
