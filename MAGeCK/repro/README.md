# MAGeCK 可复现实验目录

这个目录是对 `MAGeCK/1.md` 到 `MAGeCK/6.md` 的实跑修订版。

## 我实际修掉的关键问题
- `MAGeCK` 官方 Python 包不在当前镜像索引里，改为使用仓库内已下载的官方源码安装。
- `mageck count --trim-5 ACCG` 是错误写法；`--trim-5` 不是碱基字符串，实跑版改为让 MAGeCK 自动判定 trim。
- 复用同一个 FASTQ 路径会导致 `mageck count` 标签混乱；实跑版改成独立文件名。
- `pip install vispr` / `vispr --version` 不对；当前可用的是 `mageck-vispr`。
- `mageck mle` 在 NumPy 新版本上会因 `np.float` 报错；已在本地源码中修补。
- `MAGeCKFlute::FluteRRA()` 文稿中的参数签名已过时；实跑版改为 `ReadRRA` + 手动画图/富集。
- `enrichKEGG()` 在当前环境触发 KEGG SSL 错误；实跑版统一切到本地可跑的 `GO BP` 富集。

## 目录说明
- `analysis/`：每篇教程对应的可执行脚本。
- `data/`：官方 demo 数据和本地生成的 CRISPRi/CRISPRa/药物互作模拟数据。
- `results/`：真实运行得到的表格、统计结果和图片。
- `refs/`：轻量级 DepMap-demo 参考表。
- `run_all.sh`：一键重跑整个流程。

## 运行方式
先确保 `mageck` 和 `mageck-vispr` 已安装到 `$HOME/.local/bin`，然后在本目录运行：

```bash
./run_all.sh
```

## 主要产物
- `results/figures/Figure_main.png`
- `results/figures/Figure_main.pdf`
- `results/statistics_summary.tsv`
- `results/figures/pub_*.png`
