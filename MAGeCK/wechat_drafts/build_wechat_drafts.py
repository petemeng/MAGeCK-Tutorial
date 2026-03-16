from __future__ import annotations

import json
from pathlib import Path
from os.path import relpath
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
        'cover_image': ROOT / 'MAGeCK' / 'full' / 'reports' / 'figures' / 'article1_pub_gene_volcano_full.png',
        'hero_image': ROOT / 'MAGeCK' / 'full' / 'reports' / 'figures' / 'article1_pub_sgrna_rank_full.png',
        'digest': '基于 Sanson 2018 / SRP172473 全量原始数据，从 Brunello 文库计数、RRA 排序到 essential gene 与 GO 富集，完整走通 MAGeCK 基础分析。',
        'author': 'Songlab',
    },
    {
        'source': ROOT / 'MAGeCK' / '2.md',
        'slug': '02-mageck-mle-vispr',
        'series': 'CRISPR 筛选最佳实践（二）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_nine_quadrant.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_mle_vs_rra.png',
        'digest': '讲透 MAGeCK MLE 的设计矩阵、多条件建模与 beta score 解释，并把当前版本可运行的 mageck-vispr 工作流一起整理清楚。',
        'author': 'Songlab',
    },
    {
        'source': ROOT / 'MAGeCK' / '3.md',
        'slug': '03-mageckflute-integrative',
        'series': 'CRISPR 筛选最佳实践（三）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_flute_squareview.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_depmap_heatmap.png',
        'digest': '把 MAGeCK 结果接入 MAGeCKFlute，补上通路富集、拷贝数偏差检查和 DepMap 参考背景，做成更完整的整合分析。',
        'author': 'Songlab',
    },
    {
        'source': ROOT / 'MAGeCK' / '4.md',
        'slug': '04-crispri-analysis',
        'series': 'CRISPR 筛选最佳实践（四）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_crispri_vs_ko.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_tss_distance.png',
        'digest': '聚焦 CRISPRi：从 TSS 距离窗口、alphamedian 汇总，到与 CRISPRko 的结果对照和拷贝数假阳性解释。',
        'author': 'Songlab',
    },
    {
        'source': ROOT / 'MAGeCK' / '5.md',
        'slug': '05-drug-interaction-screen',
        'series': 'CRISPR 筛选最佳实践（五）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_diff_beta.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'pub_dose_response.png',
        'digest': '用 MLE 互作项解析 drug-gene interaction screen，区分合成致死、共同 essential 与耐药基因，并加入剂量效应验证。',
        'author': 'Songlab',
    },
    {
        'source': ROOT / 'MAGeCK' / '6.md',
        'slug': '06-publication-figures-review',
        'series': 'CRISPR 筛选最佳实践（六）',
        'cover_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'Figure_main.png',
        'hero_image': ROOT / 'MAGeCK' / 'repro' / 'results' / 'figures' / 'Figure_main.png',
        'digest': '把前 5 篇的核心结果重排成投稿级图表，同时整理 CRISPR screening 论文最常见的审稿问题与回答框架。',
        'author': 'Songlab',
    },
]

ROOT_WRAP_STYLE = "font-family:-apple-system,BlinkMacSystemFont,'PingFang SC','Hiragino Sans GB','Noto Sans CJK SC',sans-serif;color:#2d2d2d;"
TITLE_CARD_STYLE = 'margin:0 0 26px;padding:28px 24px 24px;border-radius:22px;background:#f5f0e6;border:1px solid #e7dcc7;'
TITLE_BAR_STYLE = 'width:92px;height:4px;background:#7c9970;border-radius:999px;margin:0 0 18px;'
INFO_CARD_STYLE = 'margin:22px 0;padding:18px 18px 10px;background:#f6f1e7;border:1px solid #eadfca;border-radius:14px;'
NOTE_CARD_STYLE = 'margin:18px 0 24px;padding:18px 18px 10px;background:#eef5ee;border:1px solid #d8e6d8;border-radius:14px;'
DIVIDER_STYLE = 'margin:28px auto;width:72px;height:1px;background:#d8ccb7;'
H2_STYLE = 'margin:38px 0 18px;padding:10px 16px;border-left:4px solid #7c9970;background:#eef3ea;color:#203124;font-size:22px;line-height:1.45;'
H3_STYLE = 'margin:28px 0 14px;color:#203124;font-size:19px;line-height:1.55;'
H4_STYLE = 'margin:20px 0 10px;color:#203124;font-size:17px;line-height:1.55;'
P_STYLE = 'margin:0 0 18px;color:#2d2d2d;font-size:16px;line-height:1.85;text-align:left;word-break:break-word;'
LIST_STYLE = 'margin:0 0 18px;padding-left:1.35em;color:#2d2d2d;line-height:1.85;'
ITEM_STYLE = 'margin:0 0 10px;'
QUOTE_STYLE = 'margin:18px 0;padding:16px 16px 2px;background:#eef3ea;border-left:4px solid #7c9970;border-radius:12px;'
INLINE_CODE_STYLE = 'font-family:Menlo,Consolas,monospace;font-size:14px;background:#f3efe6;color:#8a3b12;padding:2px 6px;border-radius:6px;'
CODE_STYLE = 'white-space:pre-wrap;word-break:break-word;overflow:auto;background:#1f2937;border-radius:14px;padding:16px 18px;color:#f8fafc;font-size:13px;line-height:1.75;margin:10px 0 20px 0;'
OUTPUT_STYLE = 'white-space:pre-wrap;word-break:break-word;overflow:auto;background:#f8fafc;border:1px solid #dbeafe;border-radius:14px;padding:16px 18px;color:#1f2937;font-size:13px;line-height:1.75;margin:10px 0 20px 0;'
OUTPUT_BADGE_STYLE = 'display:inline-block;margin:0 0 8px 0;padding:6px 10px;border-radius:999px;background:#dbeafe;color:#1d4ed8;font-size:12px;font-weight:700;'
LANG_BADGE_STYLE = 'display:inline-block;margin:0 0 8px 0;padding:4px 8px;border-radius:999px;background:#334155;color:#bfdbfe;font-size:11px;font-weight:700;letter-spacing:0.06em;'
IMG_STYLE = 'width:100%;display:block;border-radius:14px;box-shadow:0 10px 24px rgba(15,23,42,0.10);margin:18px auto;'
TABLE_WRAP_STYLE = 'overflow-x:auto;margin:16px 0 22px 0;border:1px solid #e5e7eb;border-radius:14px;'
TABLE_STYLE = 'width:100%;border-collapse:collapse;font-size:14px;line-height:1.7;background:#ffffff;'
TH_STYLE = 'background:#f6f1e7;border-bottom:1px solid #e5e7eb;padding:10px 12px;text-align:left;color:#203124;'
TD_STYLE = 'border-top:1px solid #f1f5f9;padding:10px 12px;color:#374151;vertical-align:top;'
LINK_STYLE = 'color:#3b6b54;text-decoration:none;border-bottom:1px solid #b7ccb7;'


def extract_title(md_text: str) -> str:
    for line in md_text.splitlines():
        if line.startswith('# '):
            return line[2:].strip()
    raise ValueError('Missing H1 title')


def normalize_display_title(title: str) -> str:
    return title.replace('最佳实践系列（', '最佳实践（')


def style_tag(tag: Tag, style: str) -> None:
    current = tag.get('style', '').strip()
    if style in current:
        return
    tag['style'] = f'{current};{style}' if current else style



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


def extract_leading_quote_blocks(md_text: str) -> list[list[str]]:
    blocks: list[list[str]] = []
    current: list[str] = []
    seen_quote = False
    for line in md_text.splitlines():
        stripped = line.rstrip()
        if stripped.startswith('>'):
            seen_quote = True
            current.append(stripped.lstrip('>').strip())
            continue
        if current:
            blocks.append(current)
            current = []
        if seen_quote and stripped.strip() and not stripped.startswith('#'):
            break
    if current:
        blocks.append(current)
    return [block for block in blocks if any(item.strip() for item in block)]


def append_markdown_fragment(parent: Tag, text: str) -> None:
    fragment = BeautifulSoup(markdown.markdown(text, extensions=['fenced_code', 'tables', 'sane_lists']), 'html.parser')
    children = list(fragment.contents)
    if not children:
        parent.append(NavigableString(text))
        return
    if len(children) == 1 and getattr(children[0], 'name', None) == 'p':
        children = list(children[0].contents)
    for child in children:
        parent.append(child)



def build_title_card(soup: BeautifulSoup, title: str) -> Tag:
    section = soup.new_tag('section')
    section['style'] = TITLE_CARD_STYLE

    bar = soup.new_tag('section')
    bar['style'] = TITLE_BAR_STYLE
    section.append(bar)

    h1 = soup.new_tag('h1')
    h1.string = title
    h1['style'] = 'margin:0;color:#203124;font-size:28px;line-height:1.45;'
    section.append(h1)
    return section



def build_divider(soup: BeautifulSoup) -> Tag:
    divider = soup.new_tag('section')
    divider['style'] = DIVIDER_STYLE
    return divider



def build_meta_card(_: Tag | BeautifulSoup, lines: list[str], note: bool = False) -> Tag:
    factory = BeautifulSoup('', 'html.parser')
    section = factory.new_tag('section')
    section['style'] = NOTE_CARD_STYLE if note else INFO_CARD_STYLE

    title = lines[0] if lines else ('✅ 实跑修订' if note else '📋 教程信息')
    p = factory.new_tag('p')
    p['style'] = P_STYLE
    append_markdown_fragment(p, title)
    section.append(p)

    items = [line[2:].strip() if line.startswith('- ') else line.strip() for line in lines[1:] if line.strip()]
    if items:
        ul = factory.new_tag('ul')
        ul['style'] = LIST_STYLE
        for item in items:
            li = factory.new_tag('li')
            li['style'] = ITEM_STYLE
            append_markdown_fragment(li, item)
            ul.append(li)
        section.append(ul)

    for code in section.find_all('code'):
        code['style'] = INLINE_CODE_STYLE
    for link in section.find_all('a'):
        style_tag(link, LINK_STYLE)
    return section



def local_src_for_render(src: str, source_path: Path, output_path: Path) -> str:
    if src.startswith('http://') or src.startswith('https://') or src.startswith('data:'):
        return src
    resolved = (source_path.parent / src).resolve()
    return relpath(resolved, output_path.parent.resolve()).replace('\\', '/')



def is_meta_quote(tag: Tag) -> bool:
    text = tag.get_text(' ', strip=True)
    return any(key in text for key in ['教程信息', '实跑修订', '可执行版代码', '一键复现'])



def beautify_html(
    html_body: str,
    cfg: dict[str, Any],
    source_path: Path,
    output_path: Path,
    quote_blocks: list[list[str]],
    display_title: str,
) -> str:
    soup = BeautifulSoup(html_body, 'html.parser')

    if soup.h1:
        soup.h1.decompose()

    for blockquote in soup.find_all('blockquote'):
        if is_meta_quote(blockquote):
            blockquote.decompose()
            continue
        style_tag(blockquote, QUOTE_STYLE)
        for inner in blockquote.find_all(['ul', 'ol']):
            style_tag(inner, LIST_STYLE)
        for inner in blockquote.find_all('li'):
            style_tag(inner, ITEM_STYLE)

    for tag in soup.find_all(['h2', 'h3', 'h4']):
        if tag.name == 'h2':
            style_tag(tag, H2_STYLE)
        elif tag.name == 'h3':
            style_tag(tag, H3_STYLE)
        else:
            style_tag(tag, H4_STYLE)

    for p in soup.find_all('p'):
        style_tag(p, P_STYLE)

    for ul in soup.find_all('ul'):
        style_tag(ul, LIST_STYLE)
    for ol in soup.find_all('ol'):
        style_tag(ol, LIST_STYLE)
    for li in soup.find_all('li'):
        style_tag(li, ITEM_STYLE)

    for hr in soup.find_all('hr'):
        divider = soup.new_tag('section')
        divider['style'] = DIVIDER_STYLE
        hr.replace_with(divider)

    for code in soup.find_all('code'):
        if code.parent.name != 'pre':
            style_tag(code, INLINE_CODE_STYLE)

    for link in soup.find_all('a'):
        style_tag(link, LINK_STYLE)

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
        if raw_text.lstrip().startswith('📊 输出：'):
            cleaned = raw_text.lstrip().replace('📊 输出：', '', 1).lstrip('\n')
            if code is not None:
                code.string = cleaned
            else:
                pre.string = cleaned
            badge = soup.new_tag('div')
            badge.string = '📊 输出'
            badge['style'] = OUTPUT_BADGE_STYLE
            pre.insert_before(badge)
            pre['style'] = OUTPUT_STYLE
        else:
            pre['style'] = CODE_STYLE
            if lang:
                badge = soup.new_tag('div')
                badge.string = lang.upper()
                badge['style'] = LANG_BADGE_STYLE
                pre.insert_before(badge)
        if code:
            code['style'] = 'font-family:Menlo,Consolas,monospace;background:transparent;padding:0;color:inherit;'

    for img in soup.find_all('img'):
        src = img.get('src', '')
        img['src'] = local_src_for_render(src, source_path, output_path)
        img['style'] = IMG_STYLE

    for table in soup.find_all('table'):
        wrapper = soup.new_tag('div')
        wrapper['style'] = TABLE_WRAP_STYLE
        table.wrap(wrapper)
        table['style'] = TABLE_STYLE
        for th in table.find_all('th'):
            th['style'] = TH_STYLE
        for td in table.find_all('td'):
            td['style'] = TD_STYLE

    outer = BeautifulSoup('', 'html.parser')
    root = outer.new_tag('section')
    root['style'] = ROOT_WRAP_STYLE

    root.append(build_title_card(outer, display_title))
    display_blocks = [
        block for block in quote_blocks
        if not any(key in ' '.join(block) for key in ['实跑修订', '可执行版代码', '真实结果目录', '一键复现'])
    ]
    for block in display_blocks:
        root.append(build_meta_card(outer, block, note=False))
    if display_blocks:
        root.append(build_divider(outer))

    previous_divider = bool(display_blocks)
    for node in list(soup.contents):
        if isinstance(node, NavigableString) and not node.strip():
            continue
        is_divider = isinstance(node, Tag) and node.name == 'section' and node.get('style') == DIVIDER_STYLE
        if is_divider and previous_divider:
            continue
        root.append(node)
        previous_divider = is_divider

    outer.append(root)
    return str(outer)



def render_article(cfg: dict[str, Any]) -> dict[str, Any]:
    source_path: Path = cfg['source']
    md_text = source_path.read_text()
    title = extract_title(md_text)
    display_title = cfg.get('wechat_title', normalize_display_title(title))
    output_path = RENDER_DIR / f"{cfg['slug']}.html"

    html_body = markdown.markdown(
        md_text,
        extensions=['fenced_code', 'tables', 'sane_lists', 'toc'],
        output_format='html5',
    )
    quote_blocks = extract_leading_quote_blocks(md_text)
    polished = beautify_html(html_body, cfg, source_path, output_path, quote_blocks, display_title)
    full_html = f'''<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{display_title}</title>
</head>
<body style="margin:0;background:#ffffff;">{polished}</body>
</html>'''
    output_path.write_text(full_html)

    return {
        'slug': cfg['slug'],
        'title': display_title,
        'wechat_title': display_title,
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
