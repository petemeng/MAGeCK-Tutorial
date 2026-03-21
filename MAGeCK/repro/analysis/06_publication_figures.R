source('analysis/_common.R')
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

full_root <- normalizePath('../full')

# Article 1 full raw data
rra_gene <- read_tsv(file.path(full_root, 'results/article1_basic_full_raw/mageck_test.gene_summary.txt'), show_col_types = FALSE)
rra_sgrna <- read_tsv(file.path(full_root, 'results/article1_basic_full_raw/mageck_test.sgrna_summary.txt'), show_col_types = FALSE)
count_summary <- read_tsv(file.path(full_root, 'counts/article1_basic_full_raw/mageck_count.countsummary.txt'), show_col_types = FALSE)
go_bp <- read_tsv(file.path(full_root, 'results/article1_basic_full_raw/go_bp.tsv'), show_col_types = FALSE)

# Article 2 full raw data
mle_gene <- read_tsv(file.path(full_root, 'results/article2_full_raw/mageck_mle.gene_summary.txt'), show_col_types = FALSE)

# Article 3 summary
article3_overlap <- read_tsv(file.path(full_root, 'results/article3_integrative/depmap_overlap_summary.tsv'), show_col_types = FALSE)

# Article 4 / 5 repro data
crispri <- read_tsv('results/crispri_test.gene_summary.txt', show_col_types = FALSE)
crisprko_repro <- read_tsv('results/mageck_test.gene_summary.txt', show_col_types = FALSE)
drug_mle <- read_tsv('results/drug_mle.gene_summary.txt', show_col_types = FALSE)

# Panel A: article 1 sgRNA rank
rank_sgrna <- rra_sgrna %>% arrange(LFC) %>% mutate(rank = row_number())
label_genes_a <- rra_gene %>% arrange(`neg|fdr`, `neg|lfc`) %>% slice_head(n = 4) %>% pull(id)
label_df_a <- rank_sgrna %>%
  filter(Gene %in% label_genes_a) %>%
  group_by(Gene) %>%
  slice_min(LFC, n = 1, with_ties = FALSE) %>%
  ungroup()

p_a <- ggplot(rank_sgrna, aes(rank, LFC)) +
  geom_point(size = 0.22, alpha = 0.3, color = 'grey72') +
  geom_point(data = filter(rank_sgrna, Gene %in% label_genes_a), color = screen_colors$essential, size = 0.85) +
  geom_text_repel(data = label_df_a, aes(label = Gene), size = 2.6, fontface = 'italic', max.overlaps = 10) +
  labs(x = 'sgRNA rank', y = 'log2 FC', title = 'Article 1 · sgRNA rank') +
  theme_screen(8)

# Panel B: article 1 volcano
label_genes_b <- unique(c(
  rra_gene %>% arrange(`neg|fdr`, `neg|lfc`) %>% slice_head(n = 4) %>% pull(id),
  rra_gene %>% arrange(`pos|fdr`, desc(`pos|lfc`)) %>% slice_head(n = 3) %>% pull(id)
))
volcano_df <- rra_gene %>% mutate(
  neg_log_fdr = -log10(`neg|fdr` + 1e-30),
  class = case_when(
    `neg|fdr` < 0.05 & `neg|lfc` < -0.5 ~ 'Essential',
    `pos|fdr` < 0.05 & `pos|lfc` > 0.5 ~ 'Enriched',
    TRUE ~ 'NS'
  ),
  label = if_else(id %in% label_genes_b, id, NA_character_)
)

p_b <- ggplot(volcano_df, aes(`neg|lfc`, neg_log_fdr, color = class)) +
  geom_point(size = 0.45, alpha = 0.62) +
  geom_text_repel(data = filter(volcano_df, !is.na(label)), aes(label = label), size = 2.5, fontface = 'italic', max.overlaps = 12) +
  scale_color_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = 'grey82')) +
  labs(x = 'Negative-selection LFC', y = '-log10(FDR)', title = 'Article 1 · gene volcano') +
  theme_screen(8)

# Panel C: article 2 MLE vs RRA
joined <- rra_gene %>%
  transmute(Gene = id, rra_neg_lfc = `neg|lfc`, rra_neg_fdr = `neg|fdr`, rra_pos_fdr = `pos|fdr`) %>%
  inner_join(mle_gene %>% transmute(Gene, mle_beta = `treatment|beta`, mle_fdr = `treatment|fdr`), by = 'Gene') %>%
  mutate(
    agreement = case_when(
      rra_neg_fdr < 0.05 & mle_fdr < 0.05 & mle_beta < -0.5 ~ 'Both essential',
      rra_neg_fdr < 0.05 & !(mle_fdr < 0.05 & mle_beta < -0.5) ~ 'RRA only essential',
      !(rra_neg_fdr < 0.05) & mle_fdr < 0.05 & mle_beta < -0.5 ~ 'MLE only essential',
      rra_pos_fdr < 0.05 & mle_fdr < 0.05 & mle_beta > 0.5 ~ 'Both enriched',
      TRUE ~ 'Other / NS'
    ),
    mle_beta_plot = pmax(pmin(mle_beta, 6), -6)
  )
label_genes_c <- unique(c(
  joined %>% filter(agreement == 'Both essential', abs(mle_beta) <= 10) %>% arrange(mle_fdr, mle_beta) %>% slice_head(n = 5) %>% pull(Gene),
  joined %>% filter(agreement == 'Both enriched') %>% arrange(mle_fdr, desc(mle_beta)) %>% slice_head(n = 3) %>% pull(Gene)
))
label_df_c <- joined %>% filter(Gene %in% label_genes_c)

p_c <- ggplot(joined, aes(rra_neg_lfc, mle_beta_plot, color = agreement)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 'dashed', linewidth = 0.3, color = 'grey55') +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', linewidth = 0.3, color = 'grey55') +
  geom_point(size = 0.55, alpha = 0.68) +
  geom_text_repel(data = label_df_c, aes(label = Gene), size = 2.4, fontface = 'italic', max.overlaps = 12) +
  scale_color_manual(values = c(
    'Both essential' = screen_colors$essential,
    'Both enriched' = screen_colors$enriched,
    'RRA only essential' = '#8bbf99',
    'MLE only essential' = '#3a7d6d',
    'Other / NS' = 'grey82'
  )) +
  labs(x = 'RRA negative-selection LFC', y = 'MLE beta (clipped to ±6)', title = 'Article 2 · MLE vs RRA') +
  theme_screen(8)

# Panel D: article 3 GO summary
p_d_df <- go_bp %>%
  slice_head(n = 8) %>%
  mutate(Description = factor(Description, levels = rev(Description)), neglog10 = -log10(p.adjust))

p_d <- ggplot(p_d_df, aes(Count, Description, size = Count, color = neglog10)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = '#8bbf99', high = '#7b2d26') +
  labs(x = 'Gene count', y = NULL, color = '-log10(adj.P)', title = 'Article 3 · GO summary') +
  theme_screen(8)

# Panel E: article 4 CRISPRi vs KO
comp <- crispri %>%
  transmute(Gene = id, i_lfc = `neg|lfc`, i_fdr = `neg|fdr`) %>%
  inner_join(crisprko_repro %>% transmute(Gene = id, ko_lfc = `neg|lfc`, ko_fdr = `neg|fdr`), by = 'Gene') %>%
  mutate(category = case_when(
    i_fdr < 0.05 & ko_fdr < 0.05 ~ 'Both',
    i_fdr < 0.05 & ko_fdr >= 0.05 ~ 'CRISPRi only',
    i_fdr >= 0.05 & ko_fdr < 0.05 ~ 'CRISPRko only',
    TRUE ~ 'NS'
  ))

p_e <- ggplot(comp, aes(ko_lfc, i_lfc, color = category)) +
  geom_point(size = 0.65, alpha = 0.55) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  scale_color_manual(values = c('Both' = screen_colors$essential, 'CRISPRi only' = screen_colors$crispri, 'CRISPRko only' = '#F39B7F', 'NS' = 'grey82')) +
  labs(x = 'CRISPRko LFC', y = 'CRISPRi LFC', title = 'Article 4 · CRISPRi vs CRISPRko') +
  theme_screen(8)

# Panel F: article 5 drug interaction
plot_df <- drug_mle %>%
  transmute(
    Gene,
    time_beta = `time|beta`,
    time_fdr = `time|fdr`,
    drug_beta = `drug|beta`,
    drug_fdr = `drug|fdr`,
    class = case_when(
      drug_fdr < 0.05 & drug_beta < -0.5 ~ 'Synthetic lethal',
      drug_fdr < 0.05 & drug_beta > 0.5 ~ 'Resistance',
      time_fdr < 0.05 & time_beta < -0.5 ~ 'Common essential',
      TRUE ~ 'NS'
    ),
    label = if_else(Gene %in% c('BRCA2', 'RAD51', 'FANCD2', 'PARP1', 'TP53BP1', 'RNF8'), Gene, NA_character_)
  )

p_f <- ggplot(plot_df, aes(time_beta, drug_beta, color = class)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_text_repel(data = filter(plot_df, !is.na(label)), aes(label = label), size = 2.4, fontface = 'italic', max.overlaps = 12) +
  scale_color_manual(values = c('Synthetic lethal' = screen_colors$synthetic_lethal, 'Resistance' = screen_colors$resistance, 'Common essential' = screen_colors$common_essential, 'NS' = 'grey84')) +
  geom_hline(yintercept = 0, color = 'grey35', linewidth = 0.3) +
  geom_vline(xintercept = 0, color = 'grey35', linewidth = 0.3) +
  labs(x = 'time|beta', y = 'drug|beta', title = 'Article 5 · differential beta') +
  theme_screen(8)

combined <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = 'bold', size = 11)))

ggsave('results/figures/Figure_main.png', combined, width = 180, height = 240, units = 'mm', dpi = 300)
ggsave('results/figures/Figure_main.pdf', combined, width = 180, height = 240, units = 'mm')
cat('Main figure exported.\n')

# Summary table
rra_essential <- sum(rra_gene$`neg|fdr` < 0.05)
rra_enriched <- sum(rra_gene$`pos|fdr` < 0.05)
mle_essential <- sum(mle_gene$`treatment|fdr` < 0.05 & mle_gene$`treatment|beta` < -0.5)
mle_enriched <- sum(mle_gene$`treatment|fdr` < 0.05 & mle_gene$`treatment|beta` > 0.5)
rra_mle_overlap <- sum(joined$agreement == 'Both essential')
crispri_essential <- sum(crispri$`neg|fdr` < 0.05)
crispri_overlap <- sum(comp$i_fdr < 0.05 & comp$ko_fdr < 0.05)
drug_sl <- sum(drug_mle$`drug|fdr` < 0.05 & drug_mle$`drug|beta` < -0.5)
drug_res <- sum(drug_mle$`drug|fdr` < 0.05 & drug_mle$`drug|beta` > 0.5)

overlap_lookup <- setNames(article3_overlap$Count, article3_overlap$Category)

stats_table <- tribble(
  ~Analysis, ~Method, ~Key_Result,
  'Article 1 count QC', 'full raw MAGeCK count', sprintf('mapped %.1f%%-%.1f%%', min(count_summary$Percentage) * 100, max(count_summary$Percentage) * 100),
  'Article 1 essential', 'MAGeCK test (RRA)', sprintf('%d essential / %d enriched', rra_essential, rra_enriched),
  'Article 2 essential', 'MAGeCK mle', sprintf('%d essential / %d enriched', mle_essential, mle_enriched),
  'Article 2 overlap', 'RRA ∩ MLE', sprintf('%d shared essential genes', rra_mle_overlap),
  'Article 3 DepMap', 'local reference overlap', sprintf('RRA %d / MLE %d / shared %d', overlap_lookup['RRA ∩ DepMap'], overlap_lookup['MLE ∩ DepMap'], overlap_lookup['RRA ∩ MLE ∩ DepMap']),
  'Article 4 CRISPRi', 'MAGeCK test', sprintf('%d essential, overlap with CRISPRko = %d', crispri_essential, crispri_overlap),
  'Article 5 synthetic lethal', 'MLE interaction', sprintf('%d genes', drug_sl),
  'Article 5 resistance', 'MLE interaction', sprintf('%d genes', drug_res)
)
write_tsv(stats_table, 'results/statistics_summary.tsv')
cat('Statistics summary rows:', nrow(stats_table), '\n')
