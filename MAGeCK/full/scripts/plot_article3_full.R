script_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub('^--file=', '', script_args[grep('^--file=', script_args)]))
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..', '..'))

source(file.path(repo_root, 'MAGeCK/repro/analysis/_common.R'))
suppressPackageStartupMessages({
  library(MAGeCKFlute)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

rra_path <- file.path(repo_root, 'MAGeCK/full/results/article1_basic_full_raw/mageck_test.gene_summary.txt')
mle_path <- file.path(repo_root, 'MAGeCK/full/results/article2_full_raw/mageck_mle.gene_summary.txt')
go_path <- file.path(repo_root, 'MAGeCK/full/results/article1_basic_full_raw/go_bp.tsv')
depmap_path <- file.path(repo_root, 'MAGeCK/repro/refs/depmap_common_essential.tsv')
fig_dir <- file.path(repo_root, 'MAGeCK/full/reports/figures')
out_dir <- file.path(repo_root, 'MAGeCK/full/results/article3_integrative')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rra_flute <- ReadRRA(rra_path)
rank_vec <- rra_flute$Score
names(rank_vec) <- rra_flute$id
png(file.path(fig_dir, 'article3_pub_flute_rankview_full.png'), width = 1800, height = 1100, res = 200)
suppressWarnings(RankView(rank_vec, top = 6, bottom = 0, cutoff = 2, main = 'MAGeCKFlute RankView (full raw cohort)'))
dev.off()

beta <- ReadBeta(mle_path)
square_df <- data.frame(Control = rep(0, nrow(beta)), Treatment = beta$treatment, Gene = beta$Gene)
png(file.path(fig_dir, 'article3_pub_flute_squareview_full.png'), width = 1600, height = 1100, res = 200)
SquareView(square_df, ctrlname = 'Control', treatname = 'Treatment', top = 6, main = 'MAGeCKFlute SquareView (full raw cohort)')
dev.off()

go <- read_tsv(go_path, show_col_types = FALSE) %>%
  slice_head(n = 10) %>%
  mutate(Description = factor(Description, levels = rev(Description)), neglog10 = -log10(p.adjust))

p_go <- ggplot(go, aes(Count, Description, size = Count, color = neglog10)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = '#8bbf99', high = '#7b2d26') +
  labs(x = 'Hit genes', y = NULL, color = '-log10(adj.P)', title = 'GO BP summary from full RRA essentials') +
  theme_screen(10)
save_plot(file.path(fig_dir, 'article3_pub_flute_go_full.png'), p_go, 8, 5.8)

depmap <- read_tsv(depmap_path, show_col_types = FALSE)
rra_raw <- read_tsv(rra_path, show_col_types = FALSE)
mle_raw <- read_tsv(mle_path, show_col_types = FALSE)
rra_ess <- rra_raw %>% filter(`neg|fdr` < 0.05) %>% pull(id)
mle_ess <- mle_raw %>% filter(`treatment|fdr` < 0.05, `treatment|beta` < -0.5) %>% pull(Gene)
dep_genes <- depmap$Gene

overlap_df <- tibble(
  Category = c('RRA ∩ DepMap', 'MLE ∩ DepMap', 'RRA ∩ MLE ∩ DepMap'),
  Count = c(
    length(intersect(rra_ess, dep_genes)),
    length(intersect(mle_ess, dep_genes)),
    length(intersect(intersect(rra_ess, mle_ess), dep_genes))
  )
) %>% mutate(Category = factor(Category, levels = Category))
write_tsv(overlap_df, file.path(out_dir, 'depmap_overlap_summary.tsv'))

p_depmap <- ggplot(overlap_df, aes(Category, Count, fill = Category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = Count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c('RRA ∩ DepMap' = '#7c9970', 'MLE ∩ DepMap' = '#b66a3c', 'RRA ∩ MLE ∩ DepMap' = '#3b6c8c')) +
  labs(x = NULL, y = 'Gene count', title = 'Overlap with local DepMap essential reference') +
  theme_screen(10) +
  theme(legend.position = 'none')
save_plot(file.path(fig_dir, 'article3_pub_depmap_overlap_full.png'), p_depmap, 7, 5)

cat('rra_flute_rows:', nrow(rra_flute), '\n')
cat('beta_rows:', nrow(beta), '\n')
cat('go_terms_plotted:', nrow(go), '\n')
cat('rra_depmap_overlap:', overlap_df$Count[overlap_df$Category == 'RRA ∩ DepMap'], '\n')
cat('mle_depmap_overlap:', overlap_df$Count[overlap_df$Category == 'MLE ∩ DepMap'], '\n')
cat('shared_depmap_overlap:', overlap_df$Count[overlap_df$Category == 'RRA ∩ MLE ∩ DepMap'], '\n')
