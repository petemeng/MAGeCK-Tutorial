---
title: 1. 基础流程：从计数到必需基因
---

# CRISPR 筛选最佳实践（一）：MAGeCK 分析——从 sgRNA 计数到必需基因

> 📋 教程信息
> - GitHub 仓库：[petemeng/MAGeCK-Tutorial](https://github.com/petemeng/MAGeCK-Tutorial)（完整代码、结果与更新记录）
> - 在线网页：[petemeng.github.io/MAGeCK-Tutorial](https://petemeng.github.io/MAGeCK-Tutorial/)（可点击阅读的网页版教程）
> - 数据来源：Sanson et al., 2018, *Nature Communications*；原始筛选测序数据：SRA `SRP172473`
> - 分析对象：Brunello modified tracrRNA 文库在 A375 细胞中的 dropout screen
> - 实验设计：1 个 pDNA 基线样本 + 3 个终点筛选重复（A/B/C）
> - 预计阅读：45 分钟 | 实操：60–120 分钟（取决于网络和磁盘）
> - 难度：⭐⭐⭐⭐（5 星制）
> - 前置知识：Linux 命令行基础，了解 CRISPR-Cas9 pooled screen 的基本概念

---

## 本篇目标

这篇不再用几百 KB 的教学玩具 FASTQ，而是直接从公开原始数据出发，完整走一遍真正的 MAGeCK 基础流程：

1. 下载 Brunello 文库注释和原始 FASTQ
2. 用 `mageck count` 从 FASTQ 生成 sgRNA count table
3. 用 `mageck test` 做负筛选 / 正筛选分析
4. 输出可直接用于论文和公众号的 rank plot、volcano plot、sgRNA 一致性图和 GO 富集图
5. 结合真实 QC 指标判断这次 screen 的质量是否过关

如果你之前已经会跑 `mageck test`，这一篇最有价值的地方不是命令本身，而是：**真实数据该怎么看、QC 该怎么看、结果该怎么解释。**

---

## 生物学背景：dropout screen 在找什么

在 pooled CRISPR knockout screen 里，每个细胞通常只携带一条 sgRNA。感染完成后，所有携带不同 sgRNA 的细胞会一起培养；如果某个基因对细胞存活或增殖至关重要，那么敲掉它的细胞就会逐渐减少，对应 sgRNA 的测序计数也会下降。

这类实验最常见的两类命中是：

- **negative selection / dropout hits**：计数显著下降，通常对应 essential genes
- **positive selection / enriched hits**：计数显著上升，通常对应抑癌、细胞周期刹车或特定压力下的耐受相关基因

MAGeCK 的核心思想是两步：

- `mageck count`：把 FASTQ 读段精确匹配到 sgRNA 文库，得到每个样本的 sgRNA 计数
- `mageck test`：比较 treatment 与 control 的计数变化，再用 RRA（Robust Rank Aggregation）把 sgRNA 层面的信号聚合到基因层面

---

## 本教程使用的数据集

本篇使用 Sanson 等人在 `SRP172473` 中公开的人类全基因组 CRISPRko screen。我们选择其中最适合作为基础教学案例的一组：

- **文库**：Brunello modified tracrRNA
- **细胞系**：A375
- **对照**：pDNA（质粒池基线）
- **筛选终点**：A375 dropout 终点样本 3 个生物学重复

### 样本对应关系

| 样本标签 | 条件 | SRA run | 说明 |
|---|---|---|---|
| `BrunelloMod_pDNA` | pDNA | `SRR8297997` | 文库初始分布 |
| `BrunelloMod_A375_repA` | dropout | `SRR8297836` + `SRR8297837` | 重复 A，原始数据分成 2 个 run |
| `BrunelloMod_A375_repB` | dropout | `SRR8297838` + `SRR8297839` | 重复 B，原始数据分成 2 个 run |
| `BrunelloMod_A375_repC` | dropout | `SRR8297840` + `SRR8297841` | 重复 C，原始数据分成 2 个 run |

本篇要回答的问题很直接：**在 A375 细胞里，哪些基因在这个筛选条件下表现出稳定而显著的 dropout？**

---

## 环境准备

```bash
# 建议新建独立环境
conda create -n mageck_full python=3.9 -y
conda activate mageck_full

pip install mageck

# 如果要重新生成图，还需要 R 包
# install.packages(c('ggplot2', 'dplyr', 'readr', 'ggrepel', 'patchwork'))
# BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))

mageck -v
```

```
📊 输出：
0.5.9.5
```

---

## Step 1：下载 Brunello 文库和原始 FASTQ

仓库里已经把全量工作流整理在 `full/`。从仓库根目录执行：

```bash
# 1) 下载并转换 Brunello 文库注释
python full/scripts/fetch_brunello_library.py

# 2) 下载 article 1 需要的全部 FASTQ，并自动把 lane-split run 合并成样本级文件
bash full/scripts/download_sra.sh --cohort article1_basic_full_raw

# 3) 快速检查文库和合并后的 FASTQ
head -5 full/external/brunello/library.tsv
wc -l full/external/brunello/library.tsv
ls -lh full/raw/article1_basic_full_raw/merged
```

```
📊 输出：
A1BG_sg1  CATCTTCTTTCACCTGAACG  A1BG
A1BG_sg2  CTCCGGGGAGAACTCCGGCG  A1BG
A1BG_sg3  TCTCCATGGTGCATCAGCAC  A1BG
A1BG_sg4  TGGAAGTCCACTCCACTCAG  A1BG
A2M_sg1   ACTGCATCTGTGCAAACGGG  A2M

77441 full/external/brunello/library.tsv

total 6.2G
-rw-rw-r-- 1 t060551 t060551 1.9G Mar 17 01:34 BrunelloMod_A375_repA.fastq.gz
-rw-rw-r-- 1 t060551 t060551 2.1G Mar 17 01:35 BrunelloMod_A375_repB.fastq.gz
-rw-rw-r-- 1 t060551 t060551 1.9G Mar 17 01:36 BrunelloMod_A375_repC.fastq.gz
-rw-rw-r-- 1 t060551 t060551 335M Mar 17 01:00 BrunelloMod_pDNA.fastq.gz
```

这一步有两个关键信息：

- 文库是**完整 Brunello**，共有 **77,441** 条 sgRNA，覆盖 **19,115** 个基因 / 分组（含 `Non-Targeting Control`）
- 原始数据不是教学玩具，而是真正的全量 FASTQ：4 个合并后样本合计约 **6.2 GB**

### 一个很容易忽略的细节：`trim-5`

CRISPR screen 的 read 里并不是从 sgRNA 第 1 个碱基直接开始读。不同载体、不同建库方式，5' 端都会有一段固定前缀。

这一批数据我先在 pDNA 样本上用 `mageck count --test-run` 扫了 `20–35` 的候选偏移窗口，主要比较每组参数的比对率、零计数 sgRNA 数和 Gini 指数。更宽的窗口在 test-run 里确实能吃到更多 reads，但也会放宽 sgRNA 起始位候选范围；最终正文保留的是已经在正式全量流程里验证稳定的保守窗口：

- **方向**：正向匹配即可，不需要 reverse-complement
- **正式分析使用的 `--trim-5`**：`23,24,25,26,27,28,29,30`

也就是说，这一批 reads 里 sgRNA 起始位并不完全固定在单一偏移量上，而是落在一个小窗口里。`--test-run` 的价值就在于先把这个窗口探出来，再进入正式统计；这个设置对比对率影响非常大。

---

## Step 2：从 FASTQ 生成 sgRNA count table

`mageck count` 的输入是文库文件和 FASTQ，输出是每条 sgRNA 在每个样本中的原始计数和归一化计数。这里直接用已经验证过的参数：

```bash
mageck count \
  -l full/external/brunello/library.tsv \
  --fastq \
    full/raw/article1_basic_full_raw/merged/BrunelloMod_pDNA.fastq.gz \
    full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repA.fastq.gz \
    full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repB.fastq.gz \
    full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repC.fastq.gz \
  --sample-label BrunelloMod_pDNA,BrunelloMod_A375_repA,BrunelloMod_A375_repB,BrunelloMod_A375_repC \
  --trim-5 23,24,25,26,27,28,29,30 \
  -n full/counts/article1_basic_full_raw/mageck_count
```

跑完以后，先看 MAGeCK 自动输出的 QC 摘要：

```bash
cat full/counts/article1_basic_full_raw/mageck_count.countsummary.txt
```

```
📊 输出：
File  Label  Reads  Mapped  Percentage  TotalsgRNAs  Zerocounts  GiniIndex
full/raw/article1_basic_full_raw/merged/BrunelloMod_pDNA.fastq.gz        BrunelloMod_pDNA        9821128   7719675   0.7860  77441  13    0.07916
full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repA.fastq.gz   BrunelloMod_A375_repA   76471324  56113866  0.7338  77441  1115  0.10780
full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repB.fastq.gz   BrunelloMod_A375_repB   85301059  61564564  0.7217  77441  1031  0.10680
full/raw/article1_basic_full_raw/merged/BrunelloMod_A375_repC.fastq.gz   BrunelloMod_A375_repC   75356900  55558003  0.7373  77441  1114  0.11380
```

### 这些 QC 指标说明了什么

先看好消息：

- pDNA 样本有 **982 万** reads，其中 **771 万**成功匹配到文库，比对率 **78.6%**
- 三个 A375 终点样本分别有 **7536 万–8530 万** reads，比对率稳定在 **72.2%–73.7%**
- 对 Brunello 这类 **20 nt guide** 文库来说，**70%+** 的比对率已经属于比较健康的范围；如果长期低于 **50%**，最先该怀疑的通常不是生物学，而是 `--trim-5` 没调对，或者文库版本与注释文件不匹配
- pDNA 的 Gini index 只有 **0.079**，说明起始文库分布非常均匀

再看真正体现 dropout 的地方：

- pDNA 中只有 **13** 条 sgRNA 为 0，几乎可以看作文库覆盖完整
- 三个终点样本各有约 **1031–1115** 条 sgRNA 归零，说明经过筛选后确实出现了大规模 depletion
- 三个终点样本 Gini 都升到 **0.106–0.114**，但仍在一个很健康的范围内，说明筛选产生了富集/耗竭，却没有把整个文库拖到极端失衡

这正是一个合格 dropout screen 应该呈现的样子：**基线文库均匀，终点样本出现系统性耗竭。**

### 看一眼 count table 的真实内容

```bash
head -5 full/counts/article1_basic_full_raw/mageck_count.count.txt | sed 's/\t/    /g'
```

```
📊 输出：
sgRNA    Gene    BrunelloMod_pDNA    BrunelloMod_A375_repA    BrunelloMod_A375_repB    BrunelloMod_A375_repC
CCDC69_sg4    CCDC69    131    1572    1322    1201
IDUA_sg3      IDUA      45     444     326     214
IFNAR2_sg3    IFNAR2    70     1221    987     766
HELT_sg3      HELT      158    720     1059    1200
```

这里看到的是原始计数，不同样本测序深度差异还没有完全折叠进统计结果里。真正的显著性判断在下一步 `mageck test` 完成。

---

## Step 3：基因层面差异分析——`mageck test`

对于这种最经典的“一个 baseline vs 多个终点重复”的设计，RRA 是最直接也最稳妥的起点：

```bash
mageck test \
  -k full/counts/article1_basic_full_raw/mageck_count.count.txt \
  -t BrunelloMod_A375_repA,BrunelloMod_A375_repB,BrunelloMod_A375_repC \
  -c BrunelloMod_pDNA \
  -n full/results/article1_basic_full_raw/mageck_test \
  --normcounts-to-file
```

跑完以后，我们先用一个很简单的脚本统计显著命中数和 top hits：

```bash
python - <<'PY'
import csv
path = 'full/results/article1_basic_full_raw/mageck_test.gene_summary.txt'
with open(path) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

neg = sorted([r for r in rows if float(r['neg|fdr']) < 0.05], key=lambda x: float(x['neg|fdr']))
pos = sorted([r for r in rows if float(r['pos|fdr']) < 0.05], key=lambda x: float(x['pos|fdr']))

print('negative hits:', len(neg))
print('positive hits:', len(pos))
print('\nTop negative hits:')
for r in neg[:10]:
    print(r['id'], r['num'], r['neg|lfc'], r['neg|fdr'], sep='\t')
print('\nTop positive hits:')
for r in pos[:10]:
    print(r['id'], r['num'], r['pos|lfc'], r['pos|fdr'], sep='\t')
PY
```

```
📊 输出：
negative hits: 1181
positive hits: 5

Top negative hits:
BUB3    4    -7.1050   0.000381
RPS3A   4    -9.4150   0.000381
RPL19   4    -6.0180   0.000381
RPP21   4    -9.0775   0.000381
NUBP1   4    -8.3540   0.000381
HSPA9   4    -8.3440   0.000381
RPL7    4    -5.8605   0.000381
ITPA    4    -6.3977   0.000381
ALDOA   4    -5.1118   0.000381
PRMT5   4    -4.7908   0.000381

Top positive hits:
TADA1   4    0.50482   0.008251
DLC1    4    1.61100   0.008251
CDKN2A  4    1.33170   0.008251
TSPY8   4    1.33360   0.028465
TAF5L   4    1.48080   0.032673
```

### 怎么读这份结果

先看 negative hits。前面几位非常符合生物学预期：

- `RPS3A`、`RPL19`、`RPL7`：核糖体相关，典型 housekeeping / translation essential genes
- `PRMT5`：RNA processing 和细胞存活相关的经典 essential hit
- `HSPA9`：线粒体伴侣蛋白，增殖细胞常见依赖因子
- `BUB3`：有丝分裂检查点核心基因，dropout 非常合理

更重要的是，这批 top hit 大多都是 **4 条 sgRNA 同方向下降**。这意味着信号不是单条 sgRNA 偶然偏低，而是一个基因内部多条 guide 的一致支持。

你还会注意到，这些 top negative genes 的 `neg|fdr` 都是同一个 `0.000381`。这不是排序出错，而是 RRA 依赖排列检验时显著性分辨率有限：多个达到最小显著性档位的基因会共享同一个最小 FDR，这在强信号 dropout screen 里很常见。

再看 positive hits：

- 数量只有 **5 个**，远少于 negative hits
- 这很符合 dropout screen 的经验：**必需基因更容易在终点样本中“掉下去”，而真正能给细胞明显生长优势的敲除位点往往少得多**
- 像 `CDKN2A`、`DLC1` 这类更偏“细胞周期刹车 / 抑癌方向”的 enriched hit，并不是完全没有生物学解释，但必须结合 **A375 的细胞背景** 和 **sgRNA-level 一致性** 再判断
- 当 positive hits 这么少时，更稳妥的做法是把它们先当作 **待复核线索**：回看单条 sgRNA、文献背景、拷贝数状态和后续验证，而不是直接把它们当成确定的增殖优势基因

另外两个结果文件也值得记一下规模：

```bash
wc -l full/results/article1_basic_full_raw/mageck_test.gene_summary.txt
wc -l full/results/article1_basic_full_raw/mageck_test.sgrna_summary.txt
```

```
📊 输出：
20115 full/results/article1_basic_full_raw/mageck_test.gene_summary.txt
77442 full/results/article1_basic_full_raw/mageck_test.sgrna_summary.txt
```

也就是说，这一轮分析最终给出了：

- **20,114** 行 gene-level 结果
- **77,441** 条 sgRNA 的 summary 结果

---

## Step 4：生成发表 / 推文可直接使用的图

我把本篇的 4 张主图整理成了一个独立脚本：

```bash
Rscript full/scripts/plot_article1_full.R
```

```
📊 输出：
sgRNA rows: 77441
gene rows: 20114
mapped essential genes: 1110
significant GO BP terms: 4039
essential genes (FDR<0.05): 1181
DepMap overlap: 38
```

这里顺手输出的 `DepMap overlap: 38`，只是拿本次 essential gene 名单去对本地 common-essential 参考集做一个快速 sanity check：方向对不对、有没有抓住一批公认依赖基因。DepMap 参考集本身怎么定义、这个数字该怎么系统解释，我放到第 3 篇再展开。

这一步会生成 4 张图：

1. sgRNA 全局 rank plot
2. gene-level volcano plot
3. top essential genes 的 sgRNA 一致性条形图
4. essential genes 的 GO Biological Process 富集图

### 图 1：sgRNA rank plot

![图 1：sgRNA rank plot](assets/figures/article1_pub_sgrna_rank_full.png)

这张图是整个 screen 最直观的“全局像”。左侧大幅下沉的 sgRNA 对应耗竭信号，右侧少量抬升的 sgRNA 对应 enriched signal。对这套数据来说，左侧尾部远比右侧明显，和前面看到的“negative hits 1181 vs positive hits 5”完全一致。

### 图 2：基因层面 volcano plot

![图 2：基因火山图](assets/figures/article1_pub_gene_volcano_full.png)

这里横轴是 negative-selection LFC，纵轴是 `-log10(FDR)`。大部分显著点都集中在左侧，说明主信号来自 dropout；标出来的 `RPS3A`、`RPL19`、`PRMT5`、`BUB3` 等，都是很有说服力的 essential genes。

### 图 3：top essential genes 的 sgRNA 一致性

![图 3：sgRNA 一致性](assets/figures/article1_pub_sgrna_barplot_full.png)

这张图在审稿和答辩里非常有用。你真正想给人看的不是“某个基因显著”，而是“**这个基因的多条 sgRNA 都在往同一个方向走**”。对 `BUB3`、`RPS3A`、`RPL19`、`RPP21`、`NUBP1`、`HSPA9` 来说，这一点非常清楚。

### 图 4：GO Biological Process 富集

![图 4：GO 富集](assets/figures/article1_pub_go_enrichment_full.png)

为了不只停留在“基因名单”，再看一下富集到的过程。我们把 `neg|fdr < 0.05` 的 essential genes 做 GO BP enrichment，前五个条目如下：

```bash
python - <<'PY'
import pandas as pd
path = 'full/results/article1_basic_full_raw/go_bp.tsv'
df = pd.read_csv(path, sep='\t')
print(df[['Description', 'Count', 'p.adjust']].head(5).to_string(index=False))
PY
```

```
📊 输出：
                         Description  Count      p.adjust
ribonucleoprotein complex biogenesis    206 1.222617e-124
                 ribosome biogenesis    146  9.181014e-94
              rRNA metabolic process    119  1.912439e-73
             cytoplasmic translation     93  4.138190e-69
                     rRNA processing    105  2.324772e-67
```

方法上，这一步用的是 `clusterProfiler::enrichGO` + `org.Hs.eg.db`，本篇选择 `ont = 'BP'`、`pAdjustMethod = 'BH'`，输入基因集是 `neg|fdr < 0.05` 的 essential genes，经 `bitr()` 映射成 Entrez ID。脚本没有额外指定 `universe`，所以背景并不是本文这 **20,114** 个 MAGeCK 结果基因，而是 `org.Hs.eg.db` 中可映射的人类基因。

这套结果非常“像样”：核糖体生物发生、rRNA 代谢、翻译过程都排在最前面，说明这次 screen 抓到的是细胞最基础、最核心的生长依赖网络。

---

## Step 5：这次筛选质量到底怎么样？

如果只看 top genes，很多坏数据也能“看起来挺像回事”。真正判断这次 screen 是否靠谱，要同时看以下几层：

### 1）文库起始状态是否均匀

- pDNA `GiniIndex = 0.07916`
- pDNA 只有 `13` 条 sgRNA 零计数

这说明文库初始分布非常健康，几乎没有明显掉队的 guide。

### 2）终点样本是否出现了真实 dropout

- 三个终点样本各有约 `1031–1115` 条 sgRNA 归零
- Gini 从 `0.079` 升到 `0.107–0.114`

这正是筛选压力生效后的典型表现。

### 3）生物学上是否讲得通

- top negative genes 以核糖体、RNA processing、细胞分裂相关基因为主
- GO 富集集中在 ribosome biogenesis、rRNA processing、translation 等核心通路
- 与 `DepMap` common essential 参考集有 **38** 个重叠，方向完全正确

也就是说，这不是“统计上显著但生物学莫名其妙”的结果，而是统计和生物学都对得上的结果。

---

## 本篇关键输出文件

```bash
du -h \
  full/counts/article1_basic_full_raw/mageck_count.count.txt \
  full/counts/article1_basic_full_raw/mageck_count.count_normalized.txt \
  full/results/article1_basic_full_raw/mageck_test.gene_summary.txt \
  full/results/article1_basic_full_raw/mageck_test.sgrna_summary.txt \
  full/results/article1_basic_full_raw/go_bp.tsv \
  full/reports/figures/article1_pub_sgrna_rank_full.png \
  full/reports/figures/article1_pub_gene_volcano_full.png \
  full/reports/figures/article1_pub_sgrna_barplot_full.png \
  full/reports/figures/article1_pub_go_enrichment_full.png
```

```
📊 输出：
2.5M  full/counts/article1_basic_full_raw/mageck_count.count.txt
6.7M  full/counts/article1_basic_full_raw/mageck_count.count_normalized.txt
1.6M  full/results/article1_basic_full_raw/mageck_test.gene_summary.txt
9.3M  full/results/article1_basic_full_raw/mageck_test.sgrna_summary.txt
872K  full/results/article1_basic_full_raw/go_bp.tsv
112K  full/reports/figures/article1_pub_sgrna_rank_full.png
764K  full/reports/figures/article1_pub_gene_volcano_full.png
48K   full/reports/figures/article1_pub_sgrna_barplot_full.png
180K  full/reports/figures/article1_pub_go_enrichment_full.png
```

如果你要继续做下游分析，第一个优先读的文件通常是：

- `mageck_test.gene_summary.txt`：基因层面的主结果
- `mageck_test.sgrna_summary.txt`：查看 guide 一致性
- `mageck_test.normalized.txt`：后续做热图、散点、聚类时最方便

---

## 本篇小结

这一篇我们用的不是 toy 数据，而是公开原始数据里真正的全量 cohort：

- **文库**：77,441 条 sgRNA，19,115 个基因 / 分组
- **原始 FASTQ**：4 个样本合计约 6.2 GB
- **count QC**：pDNA Gini 0.079，终点样本 0.107–0.114
- **RRA 结果**：negative hits 1181，positive hits 5
- **代表性 essential genes**：`BUB3`、`RPS3A`、`RPL19`、`RPP21`、`NUBP1`、`HSPA9`、`PRMT5`
- **代表性 GO 过程**：ribosome biogenesis、rRNA processing、cytoplasmic translation

如果你只想记住一句话，那就是：**一套靠谱的 MAGeCK 基础分析，首先要有健康的 pDNA 分布，其次要有一致的 sgRNA dropout，最后才是漂亮的 top hit 名单。**

---

## FAQ：常见问题

**Q1：为什么这里用 pDNA，而不是 T0 细胞样本做对照？**

很多公开 screen 会同时提供 pDNA 和细胞内基线样本。`pDNA` 更接近“文库初始分布”，适合做最基础的入口分析，也是在缺少早期细胞样本时最常见的 fallback；但如果实验里有质量可靠的 `T0/T1` 细胞基线，通常更推荐拿它做效应估计，因为它能吸收感染、药筛和早期培养阶段已经发生的一轮漂移。换句话说：`pDNA` 更像文库参照，`T0` 往往更像真正的生物学起点。

**Q2：为什么 `--trim-5` 不是一个固定数字，而是一串 `23,24,25,26,27,28,29,30`？**

因为这批 reads 里 sgRNA 起始位并不完全固定。MAGeCK 允许给一个候选窗口，它会在这些偏移量上尝试匹配。对 pooled CRISPR 数据来说，这个参数往往比很多人想象得更重要。

**Q3：negative hits 特别多，正常吗？**

对一个高质量 dropout screen 来说完全正常。真正的 essential network 往往包含上千个基因，尤其在增殖活跃的肿瘤细胞系里更常见。要警惕的不是“命中太多”，而是 top hits 和 GO 过程完全不合生物学常识。

**Q4：positive hits 只有 5 个，会不会太少？**

不会。dropout screen 天然更擅长检测 loss-of-fitness signal。能让细胞显著获得优势的 knockout 位点通常远少于造成生长缺陷的 essential genes。

**Q5：下一步该做什么？**

最自然的下一步是：

1. 对复杂设计改用 `mageck mle`
2. 用 `MAGeCKFlute` 做整合分析
3. 结合 DepMap、拷贝数和通路注释做二次筛选
4. 选 top hit 做单基因验证

---

## 本系列导航

| 篇目 | 主题 | 定位 |
|---|---|---|
| **第 1 篇** | **MAGeCK 分析——从 sgRNA 计数到必需基因** | **📍 当前阅读** |
| 第 2 篇 | MAGeCK MLE + VISPR——复杂实验设计与交互可视化 | 进阶篇 |
| 第 3 篇 | MAGeCKFlute 整合分析——基因筛选的全景图 | 整合篇 |
| 第 4 篇 | CRISPRi/CRISPRa 筛选分析策略——不切 DNA 的基因扰动 | 扩展篇 |
| 第 5 篇 | 药物-基因互作筛选与合成致死分析——一加一大于二 | 应用篇 |
| 第 6 篇 | 发表级图表与审稿人常见问题——最后一公里 | 收官篇 |
