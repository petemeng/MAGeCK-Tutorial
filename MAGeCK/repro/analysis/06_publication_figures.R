source('analysis/_common.R')
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
})

sgrna <- read_tsv('results/mageck_test.sgrna_summary.txt', show_col_types = FALSE) %>%
  arrange(LFC) %>% mutate(rank = row_number())
gene_res <- read_tsv('results/mageck_test.gene_summary.txt', show_col_types = FALSE)
mle <- read_tsv('results/mageck_mle.gene_summary.txt', show_col_types = FALSE)
depmap_ref <- read_tsv('refs/depmap_common_essential.tsv', show_col_types = FALSE)

# Panel A
label_genes <- gene_res %>% arrange(`neg|fdr`, `neg|lfc`) %>% slice_head(n = 4) %>% pull(id)
label_df <- sgrna %>% filter(Gene %in% label_genes) %>% group_by(Gene) %>% slice_min(LFC, n = 1, with_ties = FALSE) %>% ungroup()
p_a <- ggplot(sgrna, aes(rank, LFC)) +
  geom_point(size = 0.25, alpha = 0.3, color = 'grey70') +
  geom_point(data = filter(sgrna, Gene %in% label_genes), color = screen_colors$essential, size = 0.8) +
  geom_text_repel(data = label_df, aes(label = Gene), size = 2.7, fontface = 'italic', max.overlaps = 10) +
  labs(x = 'sgRNA rank', y = 'log2 FC', title = 'sgRNA rank plot') +
  theme_screen(8) + theme(legend.position = 'none')

# Panel B
gene_plot <- gene_res %>% mutate(
  neg_log_fdr = -log10(`neg|fdr` + 1e-30),
  class = case_when(`neg|fdr` < 0.05 & `neg|lfc` < -0.5 ~ 'Essential', `pos|fdr` < 0.05 & `pos|lfc` > 0.5 ~ 'Enriched', TRUE ~ 'NS'),
  label = if_else(id %in% c(label_genes, gene_res %>% arrange(`pos|fdr`) %>% slice_head(n = 3) %>% pull(id)), id, NA_character_)
)
p_b <- ggplot(gene_plot, aes(`neg|lfc`, neg_log_fdr, color = class)) +
  geom_point(size = 0.45, alpha = 0.6) +
  geom_text_repel(aes(label = label), size = 2.5, fontface = 'italic', max.overlaps = 10) +
  scale_color_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = 'grey82')) +
  labs(x = 'Negative-selection LFC', y = '-log10(FDR)', title = 'Gene volcano') +
  theme_screen(8)

# Panel C
beta_df <- mle %>% transmute(beta = `KBM7|beta`, fdr = `KBM7|fdr`, class = case_when(fdr < 0.05 & beta < -0.5 ~ 'Essential', fdr < 0.05 & beta > 0.5 ~ 'Enriched', TRUE ~ 'NS'))
p_c <- ggplot(beta_df, aes(beta, fill = class)) +
  geom_histogram(bins = 50, alpha = 0.9) +
  scale_fill_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = 'grey82')) +
  labs(x = 'KBM7 beta', y = 'Gene count', title = 'MLE beta distribution') +
  theme_screen(8)

# Panel D
essential_ids <- gene_res %>% filter(`neg|fdr` < 0.05) %>% pull(id)
gene_map <- bitr(essential_ids, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
go_bp <- enrichGO(gene = gene_map$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', qvalueCutoff = 0.2, readable = TRUE)
go_df <- as_tibble(go_bp@result) %>% slice_head(n = 8) %>% mutate(GeneCount = Count, Description = factor(Description, levels = rev(Description)))
p_d <- ggplot(go_df, aes(GeneCount, Description)) +
  geom_col(fill = screen_colors$enriched, width = 0.7) +
  labs(x = 'Mapped genes in term', y = NULL, title = 'Top GO BP terms') +
  theme_screen(8)

# Panel E
top6 <- gene_res %>% arrange(`neg|fdr`, `neg|lfc`) %>% slice_head(n = 6) %>% pull(id)
sgrna_top <- sgrna %>% filter(Gene %in% top6) %>% mutate(Gene = factor(Gene, levels = top6))
p_e <- ggplot(sgrna_top, aes(sgrna, LFC, fill = Gene)) +
  geom_col(width = 0.7) +
  facet_wrap(~Gene, scales = 'free_x', nrow = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = rep(screen_colors$essential, length(top6))) +
  labs(x = NULL, y = 'log2 FC', title = 'sgRNA consistency') +
  theme_screen(7) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none', strip.text = element_text(face = 'italic'))

# Panel F
our_essential <- gene_res %>% filter(`neg|fdr` < 0.05) %>% pull(id)
overlap_n <- length(intersect(our_essential, depmap_ref$Gene))
depmap_bar <- tibble(
  Category = c('Our essential', 'DepMap-demo essential', 'Overlap'),
  Count = c(length(our_essential), nrow(depmap_ref), overlap_n),
  Fill = c('Our', 'Ref', 'Overlap')
)
p_f <- ggplot(depmap_bar, aes(Category, Count, fill = Fill)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Our = screen_colors$essential, Ref = 'grey60', Overlap = '#8B4513')) +
  labs(x = NULL, y = 'Gene count', title = 'Reference overlap') +
  theme_screen(8) + theme(legend.position = 'none')

combined <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = 'bold', size = 11)))

ggsave('results/figures/Figure_main.png', combined, width = 180, height = 240, units = 'mm', dpi = 300)
ggsave('results/figures/Figure_main.pdf', combined, width = 180, height = 240, units = 'mm')
cat('Main figure exported.\n')

count_summary <- read_tsv('results/mageck_count.countsummary.txt', show_col_types = FALSE)
crispri <- read_tsv('results/crispri_test.gene_summary.txt', show_col_types = FALSE)
drug_mle <- read_tsv('results/drug_mle.gene_summary.txt', show_col_types = FALSE)

stats_table <- tribble(
  ~Analysis, ~Method, ~Key_Result,
  'MAGeCK count', 'Official demo FASTQ', sprintf('mapped %.1f%%-%.1f%%', min(count_summary$Percentage) * 100, max(count_summary$Percentage) * 100),
  'RRA essential genes', 'MAGeCK test', sprintf('%d genes (FDR<0.05)', sum(gene_res$`neg|fdr` < 0.05)),
  'MLE essential genes', 'MAGeCK mle', sprintf('%d genes (FDR<0.05)', sum(mle$`KBM7|fdr` < 0.05)),
  'MLE vs RRA', 'Correlation', sprintf('r = %.3f', cor(mle$`KBM7|beta`, gene_res$`neg|lfc`, use = 'complete.obs')),
  'CRISPRi essential genes', 'MAGeCK test', sprintf('%d genes (FDR<0.05)', sum(crispri$`neg|fdr` < 0.05)),
  'CRISPRi ∩ CRISPRko', 'Overlap', sprintf('%d genes', length(intersect(crispri$id[crispri$`neg|fdr` < 0.05], gene_res$id[gene_res$`neg|fdr` < 0.05]))),
  'Drug synthetic lethal', 'MLE interaction', sprintf('%d genes', sum(drug_mle$`drug|fdr` < 0.05 & drug_mle$`drug|beta` < -0.5)),
  'Drug resistance', 'MLE interaction', sprintf('%d genes', sum(drug_mle$`drug|fdr` < 0.05 & drug_mle$`drug|beta` > 0.5)),
  'DepMap-demo overlap', 'Cross-reference', sprintf('%d genes', overlap_n)
)
write_tsv(stats_table, 'results/statistics_summary.tsv')
cat('Statistics summary rows:', nrow(stats_table), '\n')
