# MAGeCK-Tutorial

A clickable documentation website for the MAGeCK / CRISPR screen tutorial series.

## Website

After GitHub Pages is enabled, the site will be available at:

- `https://petemeng.github.io/MAGeCK-Tutorial/`

## Contents

- `docs/`: MkDocs website pages
- `MAGeCK/`: original tutorial markdown and reproducible analysis assets
- `scripts/build_docs.py`: syncs tutorial markdown into the docs site
- `mkdocs.yml`: MkDocs configuration

## Local preview

```bash
python -m venv .site-venv
source .site-venv/bin/activate
pip install -r requirements-docs.txt
python scripts/build_docs.py
mkdocs serve
```
