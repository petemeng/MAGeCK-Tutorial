from __future__ import annotations

import argparse
import json
import mimetypes
import os
from pathlib import Path
from typing import Any

import requests
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / 'MAGeCK' / 'wechat_drafts' / 'output' / 'articles.json'
DEFAULT_RESULT = ROOT / 'MAGeCK' / 'wechat_drafts' / 'output' / 'publish_results.json'
DEFAULT_PREVIOUS = ROOT / 'MAGeCK' / 'wechat_drafts' / 'output' / 'publish_results.json'
TOKEN_URL = 'https://api.weixin.qq.com/cgi-bin/token'
UPLOAD_IMG_URL = 'https://api.weixin.qq.com/cgi-bin/media/uploadimg'
ADD_MATERIAL_URL = 'https://api.weixin.qq.com/cgi-bin/material/add_material'
ADD_DRAFT_URL = 'https://api.weixin.qq.com/cgi-bin/draft/add'
UPDATE_DRAFT_URL = 'https://api.weixin.qq.com/cgi-bin/draft/update'
DELETE_DRAFT_URL = 'https://api.weixin.qq.com/cgi-bin/draft/delete'


class WeChatError(RuntimeError):
    pass


def parse_json_response(resp: requests.Response, action: str) -> dict[str, Any]:
    resp.raise_for_status()
    try:
        data = json.loads(resp.content.decode('utf-8'))
    except Exception as exc:
        raise WeChatError(f'{action} returned non-UTF8/invalid JSON: {exc}') from exc
    return ensure_ok(data, action)


def post_json_utf8(url: str, payload: dict[str, Any], *, params: dict[str, Any], timeout: int, action: str) -> dict[str, Any]:
    resp = requests.post(
        url,
        params=params,
        data=json.dumps(payload, ensure_ascii=False).encode('utf-8'),
        headers={'Content-Type': 'application/json; charset=utf-8'},
        timeout=timeout,
    )
    return parse_json_response(resp, action)


def ensure_ok(data: dict[str, Any], action: str) -> dict[str, Any]:
    if 'errcode' in data and data.get('errcode', 0) not in (0,):
        raise WeChatError(f"{action} failed: errcode={data.get('errcode')} errmsg={data.get('errmsg')}")
    return data


def get_access_token(appid: str, secret: str) -> str:
    response = requests.get(
        TOKEN_URL,
        params={'grant_type': 'client_credential', 'appid': appid, 'secret': secret},
        timeout=30,
    )
    data = parse_json_response(response, 'get_access_token')
    token = data.get('access_token')
    if not token:
        raise WeChatError('Missing access_token in response')
    return token


def upload_inline_image(access_token: str, image_path: Path) -> str:
    mime = mimetypes.guess_type(str(image_path))[0] or 'application/octet-stream'
    with image_path.open('rb') as handle:
        resp = requests.post(
            UPLOAD_IMG_URL,
            params={'access_token': access_token},
            files={'media': (image_path.name, handle, mime)},
            timeout=60,
        )
    data = parse_json_response(resp, f'upload_inline_image({image_path.name})')
    url = data.get('url')
    if not url:
        raise WeChatError(f'upload_inline_image({image_path.name}) returned no url')
    return url


def upload_cover_material(access_token: str, image_path: Path) -> str:
    mime = mimetypes.guess_type(str(image_path))[0] or 'application/octet-stream'
    with image_path.open('rb') as handle:
        resp = requests.post(
            ADD_MATERIAL_URL,
            params={'access_token': access_token, 'type': 'image'},
            files={'media': (image_path.name, handle, mime)},
            timeout=60,
        )
    data = parse_json_response(resp, f'upload_cover_material({image_path.name})')
    media_id = data.get('media_id')
    if not media_id:
        raise WeChatError(f'upload_cover_material({image_path.name}) returned no media_id')
    return media_id


def convert_sections_to_divs(soup: BeautifulSoup) -> None:
    for tag_name in ['section', 'article', 'main']:
        for tag in soup.find_all(tag_name):
            tag.name = 'div'


def clean_for_wechat(html: str) -> str:
    soup = BeautifulSoup(html, 'html.parser')
    body = soup.body
    frag = BeautifulSoup(body.decode_contents() if body else html, 'html.parser')
    convert_sections_to_divs(frag)

    for p in frag.find_all('p'):
        txt = p.get_text(' ', strip=True)
        if '配图来自实跑结果' in txt:
            p.decompose()

    for blockquote in frag.find_all('blockquote'):
        txt = blockquote.get_text(' ', strip=True)
        if '教程信息' in txt or '实跑修订' in txt or '可执行版代码' in txt or '真实结果目录' in txt or '一键复现' in txt:
            blockquote.decompose()

    for div in frag.find_all('div'):
        txt = div.get_text(' ', strip=True)
        if ('公众号版草稿由本地脚本自动生成' in txt or '建议发布前在公众号后台再做一次封面与摘要微调' in txt) and len(txt) < 220:
            div.decompose()

    hrs = frag.find_all('hr')
    for hr in hrs[:1]:
        hr.decompose()

    # remove empty nodes
    for tag in frag.find_all():
        if tag.name in {'div', 'p'} and not tag.get_text(strip=True) and not tag.find('img'):
            tag.decompose()

    return ''.join(str(node) for node in frag.contents).strip()


def replace_local_images(access_token: str, html: str, html_path: Path) -> tuple[str, dict[str, str]]:
    soup = BeautifulSoup(html, 'html.parser')
    uploaded: dict[str, str] = {}
    for img in soup.find_all('img'):
        src = img.get('src', '')
        if not src or src.startswith('http://') or src.startswith('https://') or src.startswith('data:'):
            continue
        candidates = [
            (html_path.parent / src).resolve(),
            (ROOT / src).resolve(),
            (ROOT / 'MAGeCK' / src).resolve(),
        ]
        image_path = next((candidate for candidate in candidates if candidate.exists()), None)
        if image_path is None:
            raise FileNotFoundError(f'Image not found for src={src}; tried: ' + ', '.join(str(p) for p in candidates))
        if str(image_path) not in uploaded:
            uploaded[str(image_path)] = upload_inline_image(access_token, image_path)
        img['src'] = uploaded[str(image_path)]
    return str(soup), uploaded


def add_draft(access_token: str, article: dict[str, Any]) -> dict[str, Any]:
    return post_json_utf8(
        ADD_DRAFT_URL,
        {'articles': [article]},
        params={'access_token': access_token},
        timeout=60,
        action=f'add_draft({article.get("title")})',
    )


def update_draft(access_token: str, media_id: str, article: dict[str, Any], index: int = 0) -> dict[str, Any]:
    return post_json_utf8(
        UPDATE_DRAFT_URL,
        {'media_id': media_id, 'index': index, 'articles': article},
        params={'access_token': access_token},
        timeout=60,
        action=f'update_draft({article.get("title")})',
    )


def delete_draft(access_token: str, media_id: str) -> dict[str, Any]:
    return post_json_utf8(
        DELETE_DRAFT_URL,
        {'media_id': media_id},
        params={'access_token': access_token},
        timeout=60,
        action=f'delete_draft({media_id})',
    )


def build_article_payload(item: dict[str, Any], access_token: str, args: argparse.Namespace) -> tuple[dict[str, Any], dict[str, str], str]:
    html_path = ROOT / item['html_file']
    cover_path = ROOT / item['cover_image']
    html = html_path.read_text()
    cleaned = clean_for_wechat(html)
    cleaned, image_map = replace_local_images(access_token, cleaned, html_path)
    thumb_media_id = upload_cover_material(access_token, cover_path)
    article_payload = {
        'title': item.get('wechat_title') or item['title'],
        'author': item.get('author', ''),
        'digest': item.get('digest', ''),
        'content': cleaned,
        'content_source_url': '',
        'thumb_media_id': thumb_media_id,
        'need_open_comment': 1 if args.open_comment else 0,
        'only_fans_can_comment': 1 if args.fans_only_comment else 0,
    }
    return article_payload, image_map, thumb_media_id


def main() -> None:
    parser = argparse.ArgumentParser(description='Upload or update WeChat official account drafts.')
    parser.add_argument('--manifest', default=str(DEFAULT_MANIFEST), help='Path to article manifest JSON')
    parser.add_argument('--appid', default=os.environ.get('WECHAT_APP_ID'), help='WeChat appid')
    parser.add_argument('--secret', default=os.environ.get('WECHAT_APP_SECRET'), help='WeChat app secret')
    parser.add_argument('--result', default=str(DEFAULT_RESULT), help='Where to save publish results JSON')
    parser.add_argument('--previous-results', default=str(DEFAULT_PREVIOUS), help='Existing publish_results.json used for update mode')
    parser.add_argument('--mode', choices=['add', 'update'], default='add', help='Add new drafts or update existing ones')
    parser.add_argument('--open-comment', action='store_true', help='Enable comments on drafts')
    parser.add_argument('--fans-only-comment', action='store_true', help='Allow only followers to comment')
    args = parser.parse_args()

    if not args.appid or not args.secret:
        raise SystemExit('Please provide --appid/--secret or set WECHAT_APP_ID/WECHAT_APP_SECRET')

    manifest = json.loads(Path(args.manifest).read_text())
    previous = None
    if args.mode == 'update':
        previous = json.loads(Path(args.previous_results).read_text())
        if previous.get('status') != 'ok':
            raise SystemExit('previous results file does not contain successful draft uploads')
        if len(previous.get('articles', [])) != len(manifest.get('articles', [])):
            raise SystemExit('previous results and manifest article counts do not match')

    try:
        access_token = get_access_token(args.appid, args.secret)
    except Exception as exc:
        result = {
            'status': 'failed_before_publish',
            'error': str(exc),
            'hint': 'If errmsg mentions invalid ip, add this machine IP to the WeChat official account API whitelist and rerun.',
        }
        Path(args.result).write_text(json.dumps(result, ensure_ascii=False, indent=2))
        raise

    publish_results: list[dict[str, Any]] = []
    try:
        for idx, item in enumerate(manifest['articles']):
            article_payload, image_map, thumb_media_id = build_article_payload(item, access_token, args)
            if args.mode == 'update':
                old = previous['articles'][idx]
                media_id = old['media_id']
                update_draft(access_token, media_id, article_payload)
            else:
                resp = add_draft(access_token, article_payload)
                media_id = resp.get('media_id')
            publish_results.append({
                'title': item['title'],
                'wechat_title': article_payload['title'],
                'media_id': media_id,
                'thumb_media_id': thumb_media_id,
                'image_map': image_map,
            })
            print(f'{args.mode.title()} draft: {item["title"]}')
    except Exception as exc:
        Path(args.result).write_text(json.dumps({
            'status': 'failed_during_publish',
            'mode': args.mode,
            'error': str(exc),
            'articles_done': publish_results,
        }, ensure_ascii=False, indent=2))
        raise

    Path(args.result).write_text(json.dumps({'status': 'ok', 'mode': args.mode, 'articles': publish_results}, ensure_ascii=False, indent=2))
    print(f'Saved publish results to {args.result}')


if __name__ == '__main__':
    main()
