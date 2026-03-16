script_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub('^--file=', '', script_args[grep('^--file=', script_args)]))
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..', '..'))

source(file.path(repo_root, 'MAGeCK/repro/analysis/_common.R'))
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

results_dir <- file.path(repo_root, 'MAGeCK/full/results/article1_basic_full_raw')
fig_dir <- file.path(repo_root, 'MAGeCK/full/reports/figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_res <- read_tsv(file.path(results_dir, 'mageck_test.gene_summary.txt'), show_col_types = FALSE)
sgrna <- read_tsv(file.path(results_dir, 'mageck_test.sgrna_summary.txt'), show_col_types = FALSE)
depmap_ref <- read_tsv(file.path(repo_root, 'MAGeCK/repro/refs/depmap_common_essential.tsv'), show_col_types = FALSE)

cat('sgRNA rows:', nrow(sgrna), '\n')
cat('gene rows:', nrow(gene_res), '\n')

neg_hits <- gene_res %>%
  filter(`neg|fdr` < 0.05) %>%
  arrange(`neg|fdr`, `neg|lfc`)
pos_hits <- gene_res %>%
  filter(`pos|fdr` < 0.05) %>%
  arrange(`pos|fdr`, desc(`pos|lfc`))

top_neg <- neg_hits %>% slice_head(n = 6) %>% pull(id)
top_pos <- pos_hits %>% slice_head(n = 4) %>% pull(id)
label_genes <- unique(c(top_neg[1:min(4, length(top_neg))], top_pos[1:min(3, length(top_pos))]))

sgrna_rank <- sgrna %>%
  arrange(LFC) %>%
  mutate(rank = row_number(), direction = if_else(LFC < 0, 'Negative', 'Positive'))
label_df <- sgrna_rank %>%
  filter(Gene %in% label_genes) %>%
  group_by(Gene) %>%
  slice_min(abs(LFC), n = 1, with_ties = FALSE) %>%
  ungroup()

p_rank <- ggplot(sgrna_rank, aes(rank, LFC)) +
  geom_point(color = 'grey75', size = 0.3, alpha = 0.45) +
  geom_point(data = filter(sgrna_rank, Gene %in% top_neg), color = screen_colors$essential, size = 0.6, alpha = 0.7) +
  geom_point(data = filter(sgrna_rank, Gene %in% top_pos), color = screen_colors$enriched, size = 0.6, alpha = 0.7) +
  geom_text_repel(data = label_df, aes(label = Gene), size = 3, fontface = 'italic', max.overlaps = 20) +
  geom_hline(yintercept = 0, color = 'grey40', linewidth = 0.3) +
  labs(x = 'sgRNA rank', y = 'log2 fold change', title = 'MAGeCK sgRNA ranking (full raw cohort)') +
  theme_screen()
save_plot(file.path(fig_dir, 'article1_pub_sgrna_rank_full.png'), p_rank, 8, 5)

gene_plot <- gene_res %>%
  mutate(
    neg_log_fdr = -log10(`neg|fdr` + 1e-30),
    class = case_when(
      `neg|fdr` < 0.05 & `neg|lfc` < -0.5 ~ 'Essential',
      `pos|fdr` < 0.05 & `pos|lfc` > 0.5 ~ 'Enriched',
      TRUE ~ 'NS'
    ),
    label = if_else(id %in% unique(c(top_neg[1:min(5, length(top_neg))], top_pos[1:min(5, length(top_pos))])), id, NA_character_)
  )

p_volcano <- ggplot(gene_plot, aes(`neg|lfc`, neg_log_fdr, color = class)) +
  geom_point(size = 0.8, alpha = 0.65) +
  geom_text_repel(aes(label = label), size = 3, fontface = 'italic', max.overlaps = 15) +
  scale_color_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = screen_colors$ns)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  labs(x = 'Negative-selection LFC', y = '-log10(FDR)', title = 'Gene-level volcano plot (full raw cohort)') +
  theme_screen()
save_plot(file.path(fig_dir, 'article1_pub_gene_volcano_full.png'), p_volcano, 7, 6)

bar_df <- sgrna %>%
  filter(Gene %in% top_neg[1:min(6, length(top_neg))]) %>%
  mutate(Gene = factor(Gene, levels = top_neg[1:min(6, length(top_neg))]))

p_bar <- ggplot(bar_df, aes(sgrna, LFC, fill = Gene)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~Gene, scales = 'free_x', nrow = 1) +
  scale_fill_manual(values = rep(screen_colors$essential, length(unique(bar_df$Gene)))) +
  labs(x = NULL, y = 'log2 fold change', title = 'Top essential genes: sgRNA consistency (full raw cohort)') +
  theme_screen(9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none', strip.text = element_text(face = 'italic'))
save_plot(file.path(fig_dir, 'article1_pub_sgrna_barplot_full.png'), p_bar, 11, 3.6)

essential_ids <- neg_hits %>% pull(id)
gene_map <- bitr(essential_ids, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
cat('mapped essential genes:', nrow(gene_map), '\n')

go_bp <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENTREZID',
  ont = 'BP',
  pAdjustMethod = 'BH',
  qvalueCutoff = 0.2,
  readable = TRUE
)
cat('significant GO BP terms:', nrow(go_bp@result), '\n')
write_tsv(as_tibble(go_bp@result), file.path(results_dir, 'go_bp.tsv'))

p_go <- dotplot(go_bp, showCategory = 10) +
  labs(title = 'Essential genes: GO BP enrichment (full raw cohort)') +
  theme_screen(10)
save_plot(file.path(fig_dir, 'article1_pub_go_enrichment_full.png'), p_go, 8, 6)

depmap_overlap <- tibble(Gene = essential_ids) %>% inner_join(depmap_ref, by = 'Gene')
cat('essential genes (FDR<0.05):', length(essential_ids), '\n')
cat('DepMap overlap:', nrow(depmap_overlap), '\n')
