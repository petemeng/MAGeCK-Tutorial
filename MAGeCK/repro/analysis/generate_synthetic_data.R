suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

set.seed(20260315)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep('^--file=', args, value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(sub('^--file=', '', script_arg[1])) else '.'
root <- normalizePath(file.path(script_dir, '..'), mustWork = FALSE)
setwd(root)

dir.create("data/crispri", recursive = TRUE, showWarnings = FALSE)
dir.create("data/crispra", recursive = TRUE, showWarnings = FALSE)
dir.create("data/drug_screen", recursive = TRUE, showWarnings = FALSE)
dir.create("refs", recursive = TRUE, showWarnings = FALSE)

nb_counts <- function(mu, size = 20) {
  stats::rnbinom(length(mu), mu = pmax(mu, 1), size = size)
}

rra <- read_tsv("results/mageck_test.gene_summary.txt", show_col_types = FALSE)
base_genes <- rra %>%
  transmute(
    Gene = id,
    ko_lfc = `neg|lfc`,
    ko_fdr = `neg|fdr`,
    pos_fdr = `pos|fdr`
  )

# ---------------------------------------------------------------------------
# CRISPRi synthetic count table
# ---------------------------------------------------------------------------
crispri_gene_tbl <- base_genes %>%
  mutate(
    base_mu = runif(n(), 150, 900),
    ko_strength = pmin(pmax(-ko_lfc, 0), 2.2),
    crispri_strength = pmin(ko_strength * 0.75, 1.8)
  )

crispri_template <- tibble(
  sgrna_id = 1:5,
  tss_offset = c(-250, -80, 25, 150, 420)
)

crispri_counts <- crossing(Gene = crispri_gene_tbl$Gene, crispri_template) %>%
  left_join(crispri_gene_tbl, by = "Gene") %>%
  mutate(
    tss_dist = tss_offset + sample(-20:20, n(), replace = TRUE),
    sgrna_eff = rlnorm(n(), meanlog = 0, sdlog = 0.15),
    window_weight = case_when(
      tss_dist >= -50 & tss_dist <= 300 ~ 1.0,
      tss_dist >= -120 & tss_dist < -50 ~ 0.45,
      TRUE ~ 0.15
    ),
    lfc_effect = crispri_strength * window_weight,
    mu_t0 = base_mu * sgrna_eff,
    mu_t12 = pmax(mu_t0 * exp(-lfc_effect), 1),
    sgRNA = sprintf("%s_TSS1_%+d_g%d", Gene, tss_dist, sgrna_id),
    T0_rep1 = nb_counts(mu_t0, size = 22),
    T0_rep2 = nb_counts(mu_t0 * runif(n(), 0.95, 1.05), size = 22),
    T12_rep1 = nb_counts(mu_t12, size = 18),
    T12_rep2 = nb_counts(mu_t12 * runif(n(), 0.92, 1.08), size = 18)
  ) %>%
  select(sgRNA, Gene, T0_rep1, T0_rep2, T12_rep1, T12_rep2)

write_tsv(crispri_counts, "data/crispri/count_table.txt")

# ---------------------------------------------------------------------------
# CRISPRa synthetic count table
# ---------------------------------------------------------------------------
crispra_hits <- base_genes %>%
  arrange(pos_fdr) %>%
  slice_head(n = 40) %>%
  pull(Gene)

crispra_gene_tbl <- base_genes %>%
  mutate(base_mu = runif(n(), 150, 850))

crispra_template <- tibble(
  sgrna_id = 1:5,
  tss_offset = c(-350, -220, -120, 40, 180)
)

crispra_counts <- crossing(Gene = crispra_gene_tbl$Gene, crispra_template) %>%
  left_join(crispra_gene_tbl, by = "Gene") %>%
  mutate(
    tss_dist = tss_offset + sample(-25:25, n(), replace = TRUE),
    sgrna_eff = rlnorm(n(), meanlog = 0, sdlog = 0.18),
    window_weight = case_when(
      tss_dist >= -400 & tss_dist <= -50 ~ 1.0,
      tss_dist > -50 & tss_dist <= 50 ~ 0.35,
      TRUE ~ 0.12
    ),
    gain = if_else(Gene %in% crispra_hits, runif(n(), 0.8, 1.6) * window_weight, rnorm(n(), 0, 0.05)),
    mu_dmso = base_mu * sgrna_eff,
    mu_drug = pmax(mu_dmso * exp(gain), 1),
    sgRNA = sprintf("%s_TSS1_%+d_g%d", Gene, tss_dist, sgrna_id),
    DMSO_rep1 = nb_counts(mu_dmso, size = 20),
    DMSO_rep2 = nb_counts(mu_dmso * runif(n(), 0.94, 1.06), size = 20),
    Drug_rep1 = nb_counts(mu_drug, size = 18),
    Drug_rep2 = nb_counts(mu_drug * runif(n(), 0.92, 1.08), size = 18)
  ) %>%
  select(sgRNA, Gene, DMSO_rep1, DMSO_rep2, Drug_rep1, Drug_rep2)

write_tsv(crispra_counts, "data/crispra/count_table.txt")

# ---------------------------------------------------------------------------
# Drug interaction synthetic count table
# ---------------------------------------------------------------------------
common_essential <- c(
  "RPS19", "RPL11", "PCNA", "POLR2A", "EIF3A", "EEF2", "RPA1",
  "PSMD1", "UBA1", "SEC61A1", "CDK1", "RPL5"
)
synthetic_lethal <- c(
  "BRCA1", "BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D",
  "XRCC2", "XRCC3", "NBN", "ATM", "ATR", "CHEK1", "CHEK2", "MRE11",
  "RAD50", "BARD1", "BRIP1", "FANCA", "FANCB", "FANCC", "FANCD2",
  "FANCE", "FANCF", "FANCI", "FANCL", "FANCM", "ATRIP", "RBBP8",
  "EXO1", "H2AFX", "BLM", "WRN"
)
resistance_genes <- c(
  "PARP1", "PARP2", "TP53BP1", "MAD2L2", "RIF1", "SHLD1", "SHLD2",
  "SHLD3", "FAM35A", "DYNLL1", "PAXIP1", "REV7", "RNF8", "RNF168"
)

background_genes <- setdiff(base_genes$Gene, c(common_essential, synthetic_lethal, resistance_genes))
background_genes <- sample(background_genes, 450)
drug_genes <- unique(c(common_essential, synthetic_lethal, resistance_genes, background_genes))

drug_gene_tbl <- tibble(Gene = drug_genes) %>%
  mutate(
    base_mu = runif(n(), 180, 1200),
    time_beta = case_when(
      Gene %in% common_essential ~ runif(n(), -1.2, -0.7),
      TRUE ~ rnorm(n(), -0.02, 0.08)
    ),
    drug_beta = case_when(
      Gene %in% synthetic_lethal ~ runif(n(), -1.5, -0.8),
      Gene %in% resistance_genes ~ runif(n(), 0.7, 1.3),
      TRUE ~ rnorm(n(), 0, 0.07)
    )
  )

drug_counts <- crossing(Gene = drug_gene_tbl$Gene, sgrna_id = 1:4) %>%
  left_join(drug_gene_tbl, by = "Gene") %>%
  mutate(
    sgrna_eff = rlnorm(n(), meanlog = 0, sdlog = 0.12),
    mu_t0 = base_mu * sgrna_eff,
    mu_dmso = pmax(mu_t0 * exp(time_beta), 1),
    mu_drug = pmax(mu_t0 * exp(time_beta + drug_beta), 1),
    mu_drug_hi = pmax(mu_t0 * exp(time_beta + drug_beta * 1.5), 1),
    sgRNA = sprintf("%s_sg%d", Gene, sgrna_id),
    T0_rep1 = nb_counts(mu_t0, size = 22),
    T0_rep2 = nb_counts(mu_t0 * runif(n(), 0.95, 1.05), size = 22),
    DMSO_rep1 = nb_counts(mu_dmso, size = 18),
    DMSO_rep2 = nb_counts(mu_dmso * runif(n(), 0.94, 1.06), size = 18),
    Drug_rep1 = nb_counts(mu_drug, size = 18),
    Drug_rep2 = nb_counts(mu_drug * runif(n(), 0.92, 1.08), size = 18),
    DrugHi_rep1 = nb_counts(mu_drug_hi, size = 18),
    DrugHi_rep2 = nb_counts(mu_drug_hi * runif(n(), 0.90, 1.10), size = 18)
  ) %>%
  select(sgRNA, Gene, T0_rep1, T0_rep2, DMSO_rep1, DMSO_rep2, Drug_rep1, Drug_rep2, DrugHi_rep1, DrugHi_rep2)

write_tsv(drug_counts, "data/drug_screen/count_table.txt")
write_tsv(drug_counts %>% select(sgRNA, Gene), "data/drug_screen/library.tsv")

write_tsv(
  tribble(
    ~Samples, ~baseline, ~time, ~drug,
    "T0_rep1", 1, 0, 0,
    "T0_rep2", 1, 0, 0,
    "DMSO_rep1", 1, 1, 0,
    "DMSO_rep2", 1, 1, 0,
    "Drug_rep1", 1, 1, 1,
    "Drug_rep2", 1, 1, 1
  ),
  "data/drug_screen/design_matrix.txt"
)

write_tsv(
  tribble(
    ~Samples, ~baseline, ~time, ~drug_lo, ~drug_hi,
    "T0_rep1", 1, 0, 0, 0,
    "T0_rep2", 1, 0, 0, 0,
    "DMSO_rep1", 1, 1, 0, 0,
    "DMSO_rep2", 1, 1, 0, 0,
    "Drug_rep1", 1, 1, 1, 0,
    "Drug_rep2", 1, 1, 1, 0,
    "DrugHi_rep1", 1, 1, 0, 1,
    "DrugHi_rep2", 1, 1, 0, 1
  ),
  "data/drug_screen/design_matrix_dose.txt"
)

# ---------------------------------------------------------------------------
# Lightweight DepMap-like reference
# ---------------------------------------------------------------------------
our_essentials <- base_genes %>%
  arrange(ko_fdr, ko_lfc) %>%
  filter(ko_lfc < 0) %>%
  slice_head(n = 120) %>%
  pull(Gene)

depmap_common <- tibble(
  Gene = unique(c(common_essential, our_essentials[1:80], synthetic_lethal[1:10])),
  Source = "depmap_demo"
)
write_tsv(depmap_common, "refs/depmap_common_essential.tsv")

top_heatmap_genes <- unique(c(common_essential[1:6], our_essentials[1:9]))[1:15]
cell_lines <- c("HeLa", "K562", "A549", "HCT116", "MCF7", "HepG2", "U2OS", "Jurkat", "PC3", "HT29")

depmap_demo <- crossing(gene_name = top_heatmap_genes, cell_line = cell_lines) %>%
  mutate(
    dependency = case_when(
      gene_name %in% common_essential ~ rnorm(n(), -1.1, 0.12),
      TRUE ~ rnorm(n(), -0.75, 0.18)
    )
  )
write_tsv(depmap_demo, "refs/depmap_dependency_demo.tsv")

cat("Synthetic data generated.\n")
cat("  CRISPRi sgRNAs:", nrow(crispri_counts), "\n")
cat("  CRISPRa sgRNAs:", nrow(crispra_counts), "\n")
cat("  Drug screen sgRNAs:", nrow(drug_counts), "\n")
cat("  DepMap demo genes:", nrow(depmap_common), "\n")
