from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import markdown
from bs4 import BeautifulSoup, NavigableString, Tag

ROOT = Path(__file__).resolve().parents[2]
WORK_DIR = ROOT / 'MAGeCK' / 'wechat_drafts'
RENDER_DIR = WORK_DIR / 'rendered'
OUTPUT_DIR = WORK_DIR / 'output'

ARTICLE_CONFIG: list[dict[str, Any]] = [
    {
        'source': ROOT / 'MAGeCK' / '1.md',
        'slug': '01-mageck-basic',
        'series': 'CRISPR 筛选最佳实践（一）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_gene_volcano.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_sgrna_rank.png',
        'digest': 'sgRNA计数、RRA与必需基因',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK①基础流程',
    },
    {
        'source': ROOT / 'MAGeCK' / '2.md',
        'slug': '02-mageck-mle-vispr',
        'series': 'CRISPR 筛选最佳实践（二）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_nine_quadrant.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_mle_vs_rra.png',
        'digest': '设计矩阵、MLE与可视化',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK②MLE可视化',
    },
    {
        'source': ROOT / 'MAGeCK' / '3.md',
        'slug': '03-mageckflute-integrative',
        'series': 'CRISPR 筛选最佳实践（三）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_flute_squareview.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_depmap_heatmap.png',
        'digest': '通路、CN偏差与DepMap',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK③Flute整合',
    },
    {
        'source': ROOT / 'MAGeCK' / '4.md',
        'slug': '04-crispri-analysis',
        'series': 'CRISPR 筛选最佳实践（四）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_crispri_vs_ko.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_tss_distance.png',
        'digest': 'CRISPRi窗口和ko对照',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK④CRISPRi/a',
    },
    {
        'source': ROOT / 'MAGeCK' / '5.md',
        'slug': '05-drug-interaction-screen',
        'series': 'CRISPR 筛选最佳实践（五）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_diff_beta.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_dose_response.png',
        'digest': '合成致死、耐药与剂量效应',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK⑤药物互作',
    },
    {
        'source': ROOT / 'MAGeCK' / '6.md',
        'slug': '06-publication-figures-review',
        'series': 'CRISPR 筛选最佳实践（六）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'Figure_main.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'Figure_main.png',
        'digest': '投稿主图与审稿人问答',
        'author': 'SongLab-Cal',
        'wechat_title': 'MAGeCK⑥投稿图表',
    },
]

BADGE_STYLE = 'display:inline-block;padding:4px 10px;margin:0 8px 8px 0;border-radius:999px;background:#EEF2FF;color:#4338CA;font-size:12px;line-height:1.6;'


def extract_title(md_text: str) -> str:
    for line in md_text.splitlines():
        if line.startswith('# '):
            return line[2:].strip()
    raise ValueError('Missing H1 title')


def parse_meta_block(lines: list[str]) -> dict[str, str]:
    meta: dict[str, str] = {}
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('> - '):
            text = stripped[4:]
            if '：' in text:
                key, value = text.split('：', 1)
                meta[key.strip()] = value.strip()
    return meta


def build_hero(title: str, cfg: dict[str, Any], meta: dict[str, str], output_path: Path) -> str:
    badges = []
    for key in ['预计阅读', '难度', '前置知识']:
        value = meta.get(key)
        if value:
            badges.append(f'<span style="{BADGE_STYLE}">{key}：{value}</span>')
    from os.path import relpath
    hero_rel = cfg['hero_image'].relative_to(ROOT / 'MAGeCK')
    hero_src = Path(relpath(cfg['hero_image'].resolve(), output_path.parent)).as_posix()
    return f'''
    <section style="margin:0 0 28px 0;padding:24px 20px;background:linear-gradient(180deg,#F8FAFC 0%,#FFFFFF 100%);border:1px solid #E5E7EB;border-radius:20px;">
      <div style="font-size:12px;color:#2563EB;font-weight:700;letter-spacing:0.08em;text-transform:uppercase;margin-bottom:10px;">{cfg['series']}</div>
      <h1 style="font-size:28px;line-height:1.35;margin:0 0 14px 0;color:#111827;">{title}</h1>
      <p style="font-size:15px;line-height:1.8;color:#374151;margin:0 0 12px 0;">{cfg['digest']}</p>
      <div style="margin:0 0 18px 0;">{''.join(badges)}</div>
      <img src="{hero_src}" alt="hero" style="width:100%;display:block;border-radius:16px;box-shadow:0 12px 28px rgba(15,23,42,0.10);" />
      <p style="margin:10px 2px 0 2px;font-size:12px;line-height:1.7;color:#6B7280;">配图来自实跑结果：<code style="font-size:12px;background:#F3F4F6;padding:2px 6px;border-radius:6px;">{hero_rel.as_posix()}</code></p>
    </section>
    '''


def rel_for_output(target: Path, output_path: Path) -> str:
    return target.relative_to(ROOT / 'MAGeCK').as_posix() if target.is_absolute() else target.as_posix()


def local_src_for_render(src: str, source_path: Path, output_path: Path) -> str:
    if src.startswith('http://') or src.startswith('https://') or src.startswith('data:'):
        return src
    resolved = (source_path.parent / src).resolve()
    from os.path import relpath
    return Path(relpath(resolved, output_path.parent)).as_posix()


def style_tag(tag: Tag, style: str) -> None:
    existing = tag.get('style', '')
    tag['style'] = (existing + (';' if existing and not existing.endswith(';') else '') + style).strip(';')


def beautify_html(html_body: str, cfg: dict[str, Any], source_path: Path, output_path: Path) -> str:
    soup = BeautifulSoup(html_body, 'html.parser')

    if soup.h1:
        soup.h1.decompose()

    for tag in soup.find_all(['h2', 'h3', 'h4']):
        level = tag.name
        if level == 'h2':
            style_tag(tag, 'font-size:22px;line-height:1.5;color:#111827;margin:34px 0 16px;padding-left:12px;border-left:4px solid #2563EB;')
        elif level == 'h3':
            style_tag(tag, 'font-size:18px;line-height:1.6;color:#1F2937;margin:24px 0 12px;')
        else:
            style_tag(tag, 'font-size:16px;line-height:1.6;color:#1F2937;margin:18px 0 10px;')

    for p in soup.find_all('p'):
        style_tag(p, 'font-size:16px;line-height:1.9;color:#374151;margin:14px 0;word-break:break-word;')

    for ul in soup.find_all('ul'):
        style_tag(ul, 'padding-left:1.4em;margin:12px 0 16px 0;color:#374151;')
    for ol in soup.find_all('ol'):
        style_tag(ol, 'padding-left:1.4em;margin:12px 0 16px 0;color:#374151;')
    for li in soup.find_all('li'):
        style_tag(li, 'margin:8px 0;font-size:16px;line-height:1.9;')

    for hr in soup.find_all('hr'):
        style_tag(hr, 'border:none;border-top:1px solid #E5E7EB;margin:28px 0;')

    for quote in soup.find_all('blockquote'):
        style_tag(quote, 'margin:18px 0;padding:14px 16px;background:#F8FAFC;border-left:4px solid #60A5FA;border-radius:8px;color:#475569;')
        for inner in quote.find_all('p'):
            style_tag(inner, 'margin:8px 0;font-size:15px;line-height:1.9;color:#475569;')

    for code in soup.find_all('code'):
        if code.parent.name != 'pre':
            style_tag(code, 'font-family:Menlo,Consolas,monospace;font-size:14px;background:#F3F4F6;color:#BE123C;padding:2px 6px;border-radius:6px;')

    output_markers: list[Tag] = []
    for p in soup.find_all('p'):
        text = p.get_text(strip=True)
        if text == '📊 输出：':
            output_markers.append(p)
            p.clear()
            badge = soup.new_tag('span')
            badge.string = '📊 输出'
            badge['style'] = 'display:inline-block;font-size:13px;font-weight:700;line-height:1;padding:8px 12px;background:#DBEAFE;color:#1D4ED8;border-radius:999px;'
            p.append(badge)
            p['style'] = 'margin:18px 0 10px 0;'

    for pre in soup.find_all('pre'):
        code = pre.code
        lang = None
        if code and code.get('class'):
            for cls in code.get('class', []):
                if cls.startswith('language-'):
                    lang = cls.split('-', 1)[1]
                    break
                if cls and cls != 'codehilite':
                    lang = cls
        raw_text = pre.get_text()
        is_output = raw_text.lstrip().startswith('📊 输出：')
        if is_output:
            cleaned = raw_text.lstrip()
            cleaned = cleaned.replace('📊 输出：', '', 1).lstrip('\n')
            if code is not None:
                code.string = cleaned
            else:
                pre.string = cleaned
            badge = soup.new_tag('div')
            badge.string = '📊 输出'
            badge['style'] = 'display:inline-block;margin:0 0 8px 0;padding:6px 10px;border-radius:999px;background:#DBEAFE;color:#1D4ED8;font-size:12px;font-weight:700;'
            pre.insert_before(badge)
            pre['style'] = 'white-space:pre-wrap;word-break:break-word;overflow:auto;background:#F8FAFC;border:1px solid #DBEAFE;border-radius:14px;padding:16px 18px;color:#1F2937;font-size:13px;line-height:1.75;margin:10px 0 20px 0;'
        else:
            pre['style'] = 'white-space:pre-wrap;word-break:break-word;overflow:auto;background:#0F172A;border-radius:14px;padding:16px 18px;color:#E5E7EB;font-size:13px;line-height:1.75;margin:10px 0 20px 0;'
            if lang:
                label = soup.new_tag('div')
                label.string = lang.upper()
                label['style'] = 'display:inline-block;margin:0 0 8px 0;padding:4px 8px;border-radius:999px;background:#1E293B;color:#93C5FD;font-size:11px;font-weight:700;letter-spacing:0.06em;'
                pre.insert_before(label)
        if code:
            code['style'] = 'font-family:Menlo,Consolas,monospace;background:transparent;padding:0;color:inherit;'

    for img in soup.find_all('img'):
        src = img.get('src', '')
        img['src'] = local_src_for_render(src, source_path, output_path)
        img['style'] = 'width:100%;display:block;border-radius:14px;box-shadow:0 10px 24px rgba(15,23,42,0.10);margin:18px auto;'

    for table in soup.find_all('table'):
        wrapper = soup.new_tag('div')
        wrapper['style'] = 'overflow-x:auto;margin:16px 0 22px 0;border:1px solid #E5E7EB;border-radius:14px;'
        table.wrap(wrapper)
        table['style'] = 'width:100%;border-collapse:collapse;font-size:14px;line-height:1.7;background:#FFFFFF;'
        for th in table.find_all('th'):
            th['style'] = 'background:#F8FAFC;border-bottom:1px solid #E5E7EB;padding:10px 12px;text-align:left;color:#111827;'
        for td in table.find_all('td'):
            td['style'] = 'border-top:1px solid #F1F5F9;padding:10px 12px;color:#374151;vertical-align:top;'

    wrapper = BeautifulSoup('', 'html.parser')
    meta = parse_meta_block(source_path.read_text().splitlines())
    hero = BeautifulSoup(build_hero(extract_title(source_path.read_text()), cfg, meta, output_path), 'html.parser')
    container = wrapper.new_tag('section')
    container['style'] = 'max-width:820px;margin:0 auto;padding:24px 16px 40px 16px;background:#FFFFFF;'
    container.append(hero)
    for node in list(soup.contents):
        container.append(node)
    footer = BeautifulSoup(f'''<section style="margin-top:36px;padding:18px 16px;background:#F8FAFC;border-radius:16px;border:1px solid #E5E7EB;">
      <p style="margin:0 0 8px 0;font-size:15px;line-height:1.8;color:#374151;"><strong>实跑说明：</strong>本文图表与结果均来自 <code style="background:#F3F4F6;padding:2px 6px;border-radius:6px;">MAGeCK/repro</code>，公众号版草稿由本地脚本自动生成。</p>
      <p style="margin:0;font-size:13px;line-height:1.8;color:#6B7280;">建议发布前在公众号后台再做一次封面与摘要微调。</p>
    </section>''', 'html.parser')
    container.append(footer)
    wrapper.append(container)
    return str(wrapper)


def render_article(cfg: dict[str, Any]) -> dict[str, Any]:
    source_path: Path = cfg['source']
    md_text = source_path.read_text()
    title = extract_title(md_text)
    output_path = RENDER_DIR / f"{cfg['slug']}.html"

    html_body = markdown.markdown(
        md_text,
        extensions=['fenced_code', 'tables', 'sane_lists', 'toc'],
        output_format='html5',
    )
    polished = beautify_html(html_body, cfg, source_path, output_path)
    full_html = f'''<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
</head>
<body style="margin:0;background:#F3F6FB;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,'PingFang SC','Hiragino Sans GB','Microsoft YaHei',sans-serif;">{polished}</body>
</html>'''
    output_path.write_text(full_html)

    return {
        'slug': cfg['slug'],
        'title': title,
        'wechat_title': cfg.get('wechat_title', title),
        'series': cfg['series'],
        'author': cfg['author'],
        'digest': cfg['digest'],
        'source_markdown': str(source_path.relative_to(ROOT)),
        'html_file': str(output_path.relative_to(ROOT)),
        'cover_image': str(cfg['cover_image'].relative_to(ROOT)),
        'hero_image': str(cfg['hero_image'].relative_to(ROOT)),
    }


def main() -> None:
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    articles = [render_article(cfg) for cfg in ARTICLE_CONFIG]
    manifest = {
        'account': 'wechat-official-account',
        'draft_mode': 'separate_single_articles',
        'articles': articles,
    }
    (OUTPUT_DIR / 'articles.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(f'Generated {len(articles)} rendered HTML files.')
    print(OUTPUT_DIR / 'articles.json')


if __name__ == '__main__':
    main()
