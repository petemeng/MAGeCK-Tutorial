from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
manifest = json.loads((ROOT / 'MAGeCK' / 'wechat_drafts' / 'output' / 'articles.json').read_text())
out = ROOT / 'MAGeCK' / 'wechat_drafts' / 'rendered' / 'index.html'
items = []
for art in manifest['articles']:
    rel = Path(art['html_file']).name
    items.append(f'''<a href="{rel}" style="display:block;text-decoration:none;color:inherit;border:1px solid #E5E7EB;border-radius:16px;padding:18px 18px 14px 18px;background:#fff;box-shadow:0 10px 20px rgba(15,23,42,.05);margin:0 0 14px 0;">\n<div style="font-size:12px;font-weight:700;color:#2563EB;margin-bottom:8px;">{art['series']}</div>\n<div style="font-size:20px;line-height:1.5;font-weight:700;color:#111827;margin-bottom:8px;">{art['title']}</div>\n<div style="font-size:15px;line-height:1.8;color:#4B5563;">{art['digest']}</div>\n</a>''')
html = f'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>MAGeCK WeChat Draft Preview</title></head>
<body style="margin:0;background:#F3F6FB;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,'PingFang SC','Microsoft YaHei',sans-serif;">
<section style="max-width:860px;margin:0 auto;padding:28px 16px 40px;">
<h1 style="font-size:30px;line-height:1.35;color:#111827;margin:0 0 12px 0;">MAGeCK 系列公众号草稿预览</h1>
<p style="font-size:16px;line-height:1.9;color:#4B5563;margin:0 0 24px 0;">共 {len(manifest['articles'])} 篇。若微信接口白名单配置完成，可直接运行 <code style="background:#fff;padding:2px 6px;border-radius:6px;">python MAGeCK/wechat_drafts/publish_to_wechat.py</code> 上传到公众号草稿箱。</p>
{''.join(items)}
</section></body></html>'''
out.write_text(html)
print(out)
