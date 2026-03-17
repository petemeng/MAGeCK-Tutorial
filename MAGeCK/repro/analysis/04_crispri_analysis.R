source('analysis/_common.R')

counts <- read_tsv('results/crispri_test.normalized.txt', show_col_types = FALSE) %>%
  mutate(
    tss_dist = as.numeric(str_extract(sgRNA, '[+-]\\d+')),
    lfc = log2((T12_rep1 + T12_rep2 + 1) / (T0_rep1 + T0_rep2 + 1))
  )
cat('CRISPRi sgRNAs:', nrow(counts), '\n')
cat('With TSS distance:', sum(!is.na(counts$tss_dist)), '\n')

crispri <- read_tsv('results/crispri_test.gene_summary.txt', show_col_types = FALSE)
crisprko <- read_tsv('results/mageck_test.gene_summary.txt', show_col_types = FALSE)

known_essential <- crispri %>% arrange(`neg|fdr`, `neg|lfc`) %>% slice_head(n = 8) %>% pull(id)
ess_sgrna <- counts %>% filter(Gene %in% known_essential, !is.na(tss_dist))
cat('Essential-gene sgRNAs:', nrow(ess_sgrna), '\n')
cat('TSS range:', paste(range(ess_sgrna$tss_dist), collapse = ' to '), '\n')

p_tss <- ggplot(ess_sgrna, aes(tss_dist, lfc)) +
  geom_point(size = 1.2, alpha = 0.55, color = screen_colors$essential) +
  geom_smooth(method = 'loess', span = 0.5, color = '#3C5488', fill = '#3C5488', alpha = 0.18) +
  annotate('rect', xmin = -50, xmax = 300, ymin = -Inf, ymax = Inf, fill = screen_colors$essential, alpha = 0.06) +
  geom_vline(xintercept = c(-50, 300), linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  labs(x = 'Distance to TSS (bp)', y = 'log2 FC (T12 / T0)', title = 'CRISPRi efficiency depends on TSS distance') +
  theme_screen()
save_plot('results/figures/pub_tss_distance.png', p_tss, 8, 5.8)

comp <- crispri %>%
  transmute(Gene = id, i_lfc = `neg|lfc`, i_fdr = `neg|fdr`) %>%
  inner_join(
    crisprko %>% transmute(Gene = id, ko_lfc = `neg|lfc`, ko_fdr = `neg|fdr`),
    by = 'Gene'
  ) %>%
  mutate(
    category = case_when(
      i_fdr < 0.05 & ko_fdr < 0.05 ~ 'Both',
      i_fdr < 0.05 & ko_fdr >= 0.05 ~ 'CRISPRi only',
      i_fdr >= 0.05 & ko_fdr < 0.05 ~ 'CRISPRko only',
      TRUE ~ 'NS'
    )
  )
cat('Shared genes:', nrow(comp), '\n')
cat('CRISPRi essential:', sum(comp$i_fdr < 0.05), '\n')
cat('CRISPRko essential:', sum(comp$ko_fdr < 0.05), '\n')
cat('Overlap:', sum(comp$i_fdr < 0.05 & comp$ko_fdr < 0.05), '\n')

p_comp <- ggplot(comp, aes(ko_lfc, i_lfc, color = category)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey50', linewidth = 0.3) +
  scale_color_manual(values = c('Both' = screen_colors$essential, 'CRISPRi only' = screen_colors$crispri, 'CRISPRko only' = '#F39B7F', 'NS' = 'grey82')) +
  labs(x = 'CRISPRko LFC', y = 'CRISPRi LFC', title = 'CRISPRi vs CRISPRko effect size') +
  theme_screen()
save_plot('results/figures/pub_crispri_vs_ko.png', p_comp, 7, 6.5)

ko_only <- comp %>% filter(ko_fdr < 0.05 & i_fdr >= 0.05) %>% pull(Gene)
cn_tbl <- read_tsv('data/mle/cnv_data.txt', show_col_types = FALSE)
cn_check <- tibble(Gene = ko_only) %>% left_join(cn_tbl, by = c('Gene' = 'SYMBOL'))
cat('CRISPRko-only genes:', length(ko_only), '\n')
cat('Mean log2 CN of CRISPRko-only genes:', round(mean(cn_check$HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE, na.rm = TRUE), 3), '\n')
cat('Genome-wide mean log2 CN:', round(mean(cn_tbl$HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE, na.rm = TRUE), 3), '\n')

crispra <- read_tsv('results/crispra_test.normalized.txt', show_col_types = FALSE) %>%
  mutate(pos_lfc = log2((Drug_rep1 + Drug_rep2 + 1) / (DMSO_rep1 + DMSO_rep2 + 1)))

p_crispra <- ggplot(crispra, aes(pos_lfc)) +
  geom_histogram(bins = 80, fill = 'grey72', color = 'white') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = '#3C5488', linewidth = 0.5) +
  labs(x = 'log2 FC (Drug / DMSO)', y = 'sgRNA count', title = 'CRISPRa positive-selection signal is sparse') +
  theme_screen()
save_plot('results/figures/pub_crispra_dist.png', p_crispra, 8, 5)
