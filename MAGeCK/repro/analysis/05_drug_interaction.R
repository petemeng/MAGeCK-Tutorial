source('analysis/_common.R')
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

mle_res <- read_tsv('results/drug_mle.gene_summary.txt', show_col_types = FALSE)
dose_res <- read_tsv('results/drug_mle_dose.gene_summary.txt', show_col_types = FALSE)

cat('Drug MLE genes:', nrow(mle_res), '\n')

gene_classes <- mle_res %>%
  transmute(
    Gene,
    `time|beta`,
    `time|fdr`,
    `drug|beta`,
    `drug|fdr`,
    class = case_when(
      `drug|fdr` < 0.05 & `drug|beta` < -0.5 ~ 'Synthetic lethal',
      `drug|fdr` < 0.05 & `drug|beta` > 0.5 ~ 'Resistance',
      `time|fdr` < 0.05 & `time|beta` < -0.5 ~ 'Common essential',
      TRUE ~ 'NS'
    )
  )

cat('Class counts:\n')
print(table(gene_classes$class))

top_sl <- gene_classes %>% filter(class == 'Synthetic lethal') %>% arrange(`drug|beta`) %>% slice_head(n = 15)
top_res <- gene_classes %>% filter(class == 'Resistance') %>% arrange(desc(`drug|beta`)) %>% slice_head(n = 15)
cat('\nTop synthetic lethal:\n')
print(top_sl, n = 15)
cat('\nTop resistance:\n')
print(top_res, n = 15)

highlight_genes <- c('BRCA1', 'BRCA2', 'PALB2', 'PARP1', 'TP53BP1', 'ATM', 'RAD51C', 'RPS19', 'RPL11', 'PCNA')
plot_df <- gene_classes %>% filter(Gene %in% highlight_genes)

p_diff <- ggplot(gene_classes, aes(`time|beta`, `drug|beta`, color = class)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_point(data = plot_df, size = 2) +
  geom_text_repel(data = plot_df, aes(label = Gene), size = 3, fontface = 'italic', max.overlaps = 20) +
  scale_color_manual(values = c('Synthetic lethal' = screen_colors$synthetic_lethal, 'Resistance' = screen_colors$resistance, 'Common essential' = screen_colors$common_essential, 'NS' = 'grey84')) +
  geom_hline(yintercept = 0, color = 'grey35', linewidth = 0.3) +
  geom_vline(xintercept = 0, color = 'grey35', linewidth = 0.3) +
  labs(x = 'time|beta', y = 'drug|beta', title = 'Differential beta plot for drug interaction') +
  theme_screen()
save_plot('results/figures/pub_diff_beta.png', p_diff, 8.5, 7)

waterfall_df <- gene_classes %>%
  arrange(`drug|beta`) %>%
  mutate(rank = row_number(), label = if_else(Gene %in% c('BRCA1', 'BRCA2', 'PALB2', 'PARP1', 'TP53BP1', 'ATM'), Gene, NA_character_))

p_water <- ggplot(waterfall_df, aes(rank, `drug|beta`, color = class)) +
  geom_point(size = 0.5, alpha = 0.55) +
  geom_text_repel(data = filter(waterfall_df, !is.na(label)), aes(label = label), size = 3, fontface = 'italic', max.overlaps = 15) +
  scale_color_manual(values = c('Synthetic lethal' = screen_colors$synthetic_lethal, 'Resistance' = screen_colors$resistance, 'Common essential' = screen_colors$common_essential, 'NS' = 'grey84')) +
  geom_hline(yintercept = 0, color = 'grey35', linewidth = 0.3) +
  labs(x = 'Gene rank', y = 'drug|beta', title = 'Waterfall plot of drug-specific effects') +
  theme_screen()
save_plot('results/figures/pub_waterfall.png', p_water, 10, 5.5)

sl_all <- gene_classes %>% filter(class == 'Synthetic lethal') %>% arrange(`drug|beta`)
res_all <- gene_classes %>% filter(class == 'Resistance') %>% arrange(desc(`drug|beta`))
sl_ids <- bitr(sl_all$Gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
cat('Synthetic lethal gene mapping:', nrow(sl_ids), '/', nrow(sl_all), '\n')

go_sl <- enrichGO(gene = sl_ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', qvalueCutoff = 0.2, readable = TRUE)
cat('Significant GO BP terms:', nrow(go_sl@result), '\n')
write_tsv(as_tibble(go_sl@result), 'results/sl_go_bp.tsv')
write_tsv(sl_all, 'results/synthetic_lethal_genes.tsv')
write_tsv(res_all, 'results/resistance_genes.tsv')

p_sl <- dotplot(go_sl, showCategory = 10) +
  labs(title = 'Synthetic-lethal genes: GO BP enrichment') +
  theme_screen(10)
save_plot('results/figures/pub_sl_go.png', p_sl, 8, 6)

dose_df <- dose_res %>%
  transmute(Gene, lo_beta = `drug_lo|beta`, hi_beta = `drug_hi|beta`) %>%
  mutate(label = if_else(Gene %in% top_sl$Gene[1:10], Gene, NA_character_))

p_dose <- ggplot(dose_df, aes(lo_beta, hi_beta)) +
  geom_point(size = 0.7, alpha = 0.45, color = 'grey70') +
  geom_point(data = filter(dose_df, !is.na(label)), color = screen_colors$synthetic_lethal, size = 1.8) +
  geom_text_repel(data = filter(dose_df, !is.na(label)), aes(label = label), size = 3, fontface = 'italic', max.overlaps = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
  labs(x = 'Low-dose drug beta', y = 'High-dose drug beta', title = 'Dose-response of drug interaction beta') +
  theme_screen()
save_plot('results/figures/pub_dose_response.png', p_dose, 7, 6.5)
