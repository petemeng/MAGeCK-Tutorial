#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
export PATH="$HOME/.local/bin:$PATH"

mkdir -p results/figures logs

mageck count \
  -l data/basic/library.txt \
  --fastq data/basic/CTRL_rep1.fastq data/basic/CTRL_rep2.fastq data/basic/TREAT_rep1.fastq data/basic/TREAT_rep2.fastq \
  --sample-label CTRL_rep1,CTRL_rep2,TREAT_rep1,TREAT_rep2 \
  -n results/mageck_count

mageck test \
  -k data/mle/leukemia.new.csv \
  -t HL60.final,KBM7.final \
  -c HL60.initial,KBM7.initial \
  -n results/mageck_test \
  --normcounts-to-file

mageck mle \
  -k data/mle/leukemia.new.csv \
  -d data/mle/designmat.txt \
  -n results/mageck_mle \
  --cnv-norm data/mle/cnv_data.txt \
  --permutation-round 2

Rscript analysis/generate_synthetic_data.R

mageck test \
  -k data/crispri/count_table.txt \
  -t T12_rep1,T12_rep2 \
  -c T0_rep1,T0_rep2 \
  -n results/crispri_test \
  --gene-lfc-method alphamedian \
  --normcounts-to-file

mageck test \
  -k data/crispra/count_table.txt \
  -t Drug_rep1,Drug_rep2 \
  -c DMSO_rep1,DMSO_rep2 \
  -n results/crispra_test \
  --normcounts-to-file

mageck mle \
  -k data/drug_screen/count_table.txt \
  -d data/drug_screen/design_matrix.txt \
  -n results/drug_mle \
  --permutation-round 2

mageck mle \
  -k data/drug_screen/count_table.txt \
  -d data/drug_screen/design_matrix_dose.txt \
  -n results/drug_mle_dose \
  --permutation-round 2

for script in \
  analysis/01_basic_visualization.R \
  analysis/02_mle_analysis.R \
  analysis/03_mageckflute_analysis.R \
  analysis/04_crispri_analysis.R \
  analysis/05_drug_interaction.R \
  analysis/06_publication_figures.R; do
  Rscript "$script"
done

echo "Done. Results are in $ROOT/results"
