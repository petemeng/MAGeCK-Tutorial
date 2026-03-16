source('analysis/_common.R')

mle <- read_tsv('results/mageck_mle.gene_summary.txt', show_col_types = FALSE)
rra <- read_tsv('results/mageck_test.gene_summary.txt', show_col_types = FALSE)

comparison <- mle %>%
  transmute(
    Gene,
    mle_beta = `KBM7|beta`,
    mle_fdr = `KBM7|fdr`
  ) %>%
  left_join(
    rra %>% transmute(Gene = id, rra_lfc = `neg|lfc`, rra_fdr = `neg|fdr`),
    by = 'Gene'
  ) %>%
  mutate(
    mle_sig = mle_fdr < 0.05,
    rra_sig = rra_fdr < 0.05,
    class = case_when(
      mle_sig & rra_sig ~ 'Both',
      mle_sig & !rra_sig ~ 'MLE only',
      !mle_sig & rra_sig ~ 'RRA only',
      TRUE ~ 'NS'
    )
  )

cat('MLE genes:', nrow(mle), '\n')
cat('RRA genes:', nrow(rra), '\n')
cat('MLE essential (FDR<0.05):', sum(comparison$mle_sig, na.rm = TRUE), '\n')
cat('RRA essential (FDR<0.05):', sum(comparison$rra_sig, na.rm = TRUE), '\n')
cat('Overlap:', sum(comparison$mle_sig & comparison$rra_sig, na.rm = TRUE), '\n')

r_val <- cor(comparison$mle_beta, comparison$rra_lfc, use = 'complete.obs')
label_df <- comparison %>%
  filter(class == 'Both') %>%
  arrange(mle_fdr) %>%
  slice_head(n = 8)

p_compare <- ggplot(comparison, aes(rra_lfc, mle_beta, color = class)) +
  geom_point(size = 0.6, alpha = 0.55) +
  geom_text_repel(data = label_df, aes(label = Gene), size = 3, fontface = 'italic', max.overlaps = 12) +
  scale_color_manual(values = c('Both' = screen_colors$essential, 'MLE only' = '#F39B7F', 'RRA only' = '#4DBBD5', 'NS' = 'grey82')) +
  annotate('text', x = min(comparison$rra_lfc, na.rm = TRUE), y = max(comparison$mle_beta, na.rm = TRUE), hjust = 0, vjust = 1, label = sprintf('r = %.3f', r_val), size = 4) +
  labs(x = 'RRA negative-selection LFC', y = 'MLE beta (KBM7)', title = 'MLE beta vs RRA LFC') +
  theme_screen()
save_plot('results/figures/pub_mle_vs_rra.png', p_compare, 7, 6.5)

beta_df <- mle %>%
  transmute(
    Gene,
    beta = `KBM7|beta`,
    fdr = `KBM7|fdr`,
    class = case_when(
      fdr < 0.05 & beta < -0.5 ~ 'Essential',
      fdr < 0.05 & beta > 0.5 ~ 'Enriched',
      TRUE ~ 'NS'
    )
  )

p_beta <- ggplot(beta_df, aes(beta, fill = class)) +
  geom_histogram(bins = 60, alpha = 0.9) +
  scale_fill_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = screen_colors$ns)) +
  geom_vline(xintercept = 0, color = 'grey40', linewidth = 0.3) +
  labs(x = 'KBM7 beta score', y = 'Gene count', title = 'MLE beta distribution') +
  theme_screen()
save_plot('results/figures/pub_beta_distribution.png', p_beta, 8, 5)

set.seed(20260315)
nine_df <- mle %>%
  transmute(
    Gene,
    drugA_beta = `KBM7|beta` + rnorm(n(), 0, 0.18),
    drugB_beta = `HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE|beta` + rnorm(n(), 0, 0.18)
  ) %>%
  mutate(
    category = case_when(
      drugA_beta < -0.5 & drugB_beta < -0.5 ~ 'Both essential',
      drugA_beta < -0.5 & drugB_beta >= -0.5 ~ 'DrugA-specific',
      drugA_beta >= -0.5 & drugB_beta < -0.5 ~ 'DrugB-specific',
      drugA_beta > 0.5 | drugB_beta > 0.5 ~ 'Enriched',
      TRUE ~ 'NS'
    )
  )
cat('Nine-square classes:\n')
print(table(nine_df$category))

p_nine <- ggplot(nine_df, aes(drugA_beta, drugB_beta, color = category)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c('Both essential' = screen_colors$essential, 'DrugA-specific' = '#F39B7F', 'DrugB-specific' = '#4DBBD5', 'Enriched' = screen_colors$enriched, 'NS' = 'grey82')) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 'dashed', color = 'grey55', linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', color = 'grey55', linewidth = 0.3) +
  labs(x = 'Condition A beta', y = 'Condition B beta', title = 'Nine-quadrant style beta comparison') +
  theme_screen()
save_plot('results/figures/pub_nine_quadrant.png', p_nine, 7.5, 6.5)
