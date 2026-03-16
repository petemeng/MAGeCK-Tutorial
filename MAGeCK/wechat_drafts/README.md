# WeChat Draft Package

这个目录保存了把 `MAGeCK/1.md` 到 `MAGeCK/6.md` 渲染成微信公众号草稿所需的全部文件。

## 目录说明

- `build_wechat_drafts.py`：把 Markdown 渲染成公众号友好的 HTML。
- `publish_to_wechat.py`：调用微信公众号草稿 API 上传正文图片、封面并创建草稿。
- `rendered/`：生成的本地预览 HTML。
- `output/articles.json`：文章清单与元数据。
- `output/publish_results.json`：实际上传后的返回结果。

## 生成草稿 HTML

```bash
python MAGeCK/wechat_drafts/build_wechat_drafts.py
```

## 上传到微信公众号草稿箱

```bash
export WECHAT_APP_ID='你的 appid'
export WECHAT_APP_SECRET='你的 secret'
python MAGeCK/wechat_drafts/publish_to_wechat.py
```

## 当前已知限制

如果报错包含 `invalid ip`，说明当前机器 IP 没加到公众号后台的接口白名单里。需要先在微信公众平台后台添加当前出口 IP，再重新执行上传脚本。
