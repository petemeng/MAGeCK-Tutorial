from __future__ import annotations

from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
MAG = ROOT / 'MAGeCK'
DOCS = ROOT / 'docs'
FIG_DST = DOCS / 'assets' / 'figures'

PAGES = [
    ('1.md', '01-basic-mageck.md', '1. 基础流程：从计数到必需基因'),
    ('2.md', '02-mle-vispr.md', '2. MLE 与 VISPR'),
    ('3.md', '03-flute-integrative.md', '3. MAGeCKFlute 整合分析'),
    ('4.md', '04-crispri-crispra.md', '4. CRISPRi / CRISPRa'),
    ('5.md', '05-drug-interaction.md', '5. 药物互作与合成致死'),
    ('6.md', '06-publication-figures.md', '6. 投稿图表与审稿问答'),
]

FIGURES = [
    'Figure_main.png',
    'pub_beta_distribution.png',
    'pub_cn_bias.png',
    'pub_crispra_dist.png',
    'pub_crispri_vs_ko.png',
    'pub_depmap_heatmap.png',
    'pub_diff_beta.png',
    'pub_dose_response.png',
    'pub_flute_kegg.png',
    'pub_flute_rankview.png',
    'pub_flute_squareview.png',
    'pub_gene_volcano.png',
    'pub_go_enrichment.png',
    'pub_mle_vs_rra.png',
    'pub_nine_quadrant.png',
    'pub_sgrna_barplot.png',
    'pub_sgrna_rank.png',
    'pub_sl_kegg.png',
    'pub_tss_distance.png',
    'pub_waterfall.png',
]


def transform_markdown(text: str) -> str:
    text = text.replace('(repro/results/figures/', '(assets/figures/')
    text = text.replace('`MAGeCK/repro/', '`repro/')
    text = text.replace('`MAGeCK/', '`')
    text = text.replace('MAGeCK/repro/results/', 'repro/results/')
    # remove huge manual truncation artifacts if present
    text = text.replace('…10689 chars truncated…', '')
    return text


def main() -> None:
    FIG_DST.mkdir(parents=True, exist_ok=True)
    for fig in FIGURES:
        src = MAG / 'repro' / 'results' / 'figures' / fig
        if src.exists():
            shutil.copy2(src, FIG_DST / fig)

    for src_name, dst_name, nav_title in PAGES:
        text = (MAG / src_name).read_text()
        text = transform_markdown(text)
        front_matter = f'---\ntitle: {nav_title}\n---\n\n'
        (DOCS / dst_name).write_text(front_matter + text)


if __name__ == '__main__':
    main()
