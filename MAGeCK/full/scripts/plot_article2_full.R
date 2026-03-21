script_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub('^--file=', '', script_args[grep('^--file=', script_args)]))
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..', '..'))

source(file.path(repo_root, 'MAGeCK/repro/analysis/_common.R'))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggrepel)
  library(tidyr)
  library(scales)
})

rra_path <- file.path(repo_root, 'MAGeCK/full/results/article1_basic_full_raw/mageck_test.gene_summary.txt')
mle_path <- file.path(repo_root, 'MAGeCK/full/results/article2_full_raw/mageck_mle.gene_summary.txt')
fig_dir <- file.path(repo_root, 'MAGeCK/full/reports/figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

rra <- read_tsv(rra_path, show_col_types = FALSE) %>%
  transmute(
    Gene = id,
    rra_neg_lfc = `neg|lfc`,
    rra_neg_fdr = `neg|fdr`,
    rra_pos_lfc = `pos|lfc`,
    rra_pos_fdr = `pos|fdr`
  )

mle <- read_tsv(mle_path, show_col_types = FALSE) %>%
  transmute(
    Gene,
    sgRNA,
    mle_beta = `treatment|beta`,
    mle_fdr = `treatment|fdr`,
    mle_wald_fdr = `treatment|wald-fdr`
  )

joined <- inner_join(rra, mle, by = 'Gene') %>%
  mutate(
    rra_class = case_when(
      rra_neg_fdr < 0.05 & rra_neg_lfc < -0.5 ~ 'RRA essential',
      rra_pos_fdr < 0.05 & rra_pos_lfc > 0.5 ~ 'RRA enriched',
      TRUE ~ 'RRA NS'
    ),
    mle_class = case_when(
      mle_fdr < 0.05 & mle_beta < -0.5 ~ 'MLE essential',
      mle_fdr < 0.05 & mle_beta > 0.5 ~ 'MLE enriched',
      TRUE ~ 'MLE NS'
    ),
    agreement = case_when(
      rra_class == 'RRA essential' & mle_class == 'MLE essential' ~ 'Both essential',
      rra_class == 'RRA enriched' & mle_class == 'MLE enriched' ~ 'Both enriched',
      rra_class == 'RRA essential' & mle_class == 'MLE NS' ~ 'RRA only essential',
      rra_class == 'RRA NS' & mle_class == 'MLE essential' ~ 'MLE only essential',
      rra_class == 'RRA enriched' & mle_class == 'MLE NS' ~ 'RRA only enriched',
      rra_class == 'RRA NS' & mle_class == 'MLE enriched' ~ 'MLE only enriched',
      TRUE ~ 'Other / NS'
    )
  )

both_essential <- joined %>% filter(agreement == 'Both essential') %>% arrange(mle_fdr, mle_beta)
both_enriched <- joined %>% filter(agreement == 'Both enriched') %>% arrange(mle_fdr, desc(mle_beta))
label_genes <- unique(c(head(filter(both_essential, abs(mle_beta) <= 10)$Gene, 6), head(both_enriched$Gene, 4)))
label_df <- joined %>% filter(Gene %in% label_genes)

cat('joined genes:', nrow(joined), '\n')
cat('RRA essential:', sum(joined$rra_class == 'RRA essential'), '\n')
cat('MLE essential:', sum(joined$mle_class == 'MLE essential'), '\n')
cat('RRA enriched:', sum(joined$rra_class == 'RRA enriched'), '\n')
cat('MLE enriched:', sum(joined$mle_class == 'MLE enriched'), '\n')
cat('Both essential overlap:', sum(joined$agreement == 'Both essential'), '\n')
cat('Both enriched overlap:', sum(joined$agreement == 'Both enriched'), '\n')

agreement_cols <- c(
  'Both essential' = screen_colors$essential,
  'Both enriched' = screen_colors$enriched,
  'RRA only essential' = '#8bbf99',
  'MLE only essential' = '#3a7d6d',
  'RRA only enriched' = '#e3a870',
  'MLE only enriched' = '#cf7b45',
  'Other / NS' = screen_colors$ns
)

scatter_df <- joined %>% mutate(mle_beta_plot = pmax(pmin(mle_beta, 6), -6))
label_df <- scatter_df %>% filter(Gene %in% label_genes)

p_scatter <- ggplot(scatter_df, aes(rra_neg_lfc, mle_beta_plot, color = agreement)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 'dashed', linewidth = 0.3, color = 'grey55') +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', linewidth = 0.3, color = 'grey55') +
  geom_point(size = 0.75, alpha = 0.72) +
  geom_text_repel(data = label_df, aes(label = Gene), size = 3, fontface = 'italic', max.overlaps = 15) +
  scale_color_manual(values = agreement_cols) +
  labs(
    x = 'RRA negative-selection LFC',
    y = 'MLE treatment beta (clipped to ±6 for display)',
    title = 'Full raw cohort: MLE vs RRA agreement'
  ) +
  theme_screen(10)
save_plot(file.path(fig_dir, 'article2_pub_mle_vs_rra_full.png'), p_scatter, 8, 6)

beta_df <- joined %>%
  mutate(beta_group = case_when(
    mle_fdr < 0.05 & mle_beta < -0.5 ~ 'Essential',
    mle_fdr < 0.05 & mle_beta > 0.5 ~ 'Enriched',
    TRUE ~ 'NS'
  ))

beta_xlim <- c(-6, 1.5)
outside_left <- sum(beta_df$mle_beta < beta_xlim[1])
outside_right <- sum(beta_df$mle_beta > beta_xlim[2])

p_beta <- ggplot(beta_df, aes(mle_beta, fill = beta_group)) +
  geom_histogram(binwidth = 0.1, color = 'white', alpha = 0.9) +
  scale_fill_manual(values = c(Essential = screen_colors$essential, Enriched = screen_colors$enriched, NS = screen_colors$ns)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dashed', linewidth = 0.3, color = 'grey55') +
  coord_cartesian(xlim = beta_xlim) +
  annotate('text', x = -5.8, y = Inf, label = sprintf('display clipped: %d left, %d right', outside_left, outside_right), hjust = 0, vjust = 1.5, size = 3.2, color = 'grey80') +
  labs(x = 'MLE treatment beta', y = 'Gene count', title = 'Distribution of MLE beta scores (full raw cohort)') +
  theme_screen(10)
save_plot(file.path(fig_dir, 'article2_pub_beta_distribution_full.png'), p_beta, 8, 5)

agree_counts <- joined %>%
  count(agreement) %>%
  mutate(agreement = factor(agreement, levels = names(agreement_cols)))

p_bar <- ggplot(agree_counts, aes(agreement, n, fill = agreement)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n)), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = agreement_cols) +
  labs(x = NULL, y = 'Gene count', title = 'RRA / MLE concordance summary') +
  theme_screen(9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = 'none')
save_plot(file.path(fig_dir, 'article2_pub_agreement_bar_full.png'), p_bar, 8, 5)
