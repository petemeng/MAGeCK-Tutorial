---
title: 4. CRISPRi / CRISPRa
---

# CRISPR 筛选最佳实践（四）：CRISPRi/CRISPRa 筛选分析策略——不切 DNA 的基因扰动

> 📋 教程信息
> - GitHub 仓库：[petemeng/MAGeCK-Tutorial](https://github.com/petemeng/MAGeCK-Tutorial)（完整代码、结果与更新记录）
> - 在线网页：[petemeng.github.io/MAGeCK-Tutorial](https://petemeng.github.io/MAGeCK-Tutorial/)（可点击阅读的网页版教程）
> - 数据来源：仓库内可复现教学数据与实跑结果：`repro/data/crispri/`、`repro/data/crispra/`
> - 分析对象：CRISPRi 的 TSS 距离效应、CRISPRi vs CRISPRko 对照、CRISPRa 的稀疏富集信号
> - 本篇重点：`analysis/04_crispri_analysis.R` 的实跑输出与图表解释
> - 预计阅读：35 分钟 | 实操：10–20 分钟
> - 难度：⭐⭐⭐⭐（5 星制）
> - 前置知识：完成第 1 篇，理解 dropout screen 的基本判读逻辑

---

## 本篇目标

CRISPRko 的直觉最强：Cas9 把 DNA 切断，基因失活，sgRNA 计数跟着变。但很多真实问题并不适合直接切 DNA：

- 你想抑制表达，而不是永久敲除
- 你想激活一个基因，而不是让它失活
- 你担心 copy-number 放大区域带来的 DSB 假阳性

这就是 CRISPRi / CRISPRa 的价值。

这篇不再沿用错误的外部 accession 叙事，而是直接基于仓库里已经跑通的教学结果，回答三个实用问题：

1. CRISPRi 的效果为什么强依赖 TSS 距离？
2. CRISPRi 和 CRISPRko 的 essential signal 有多少重合？
3. CRISPRa 的正筛选为什么通常比 knockout dropout 稀疏得多？

---

## 本篇使用的真实入口

分析脚本：`repro/analysis/04_crispri_analysis.R`

直接运行：

```bash
cd MAGeCK/repro
Rscript analysis/04_crispri_analysis.R
```

```
📊 输出：
CRISPRi sgRNAs: 5000
With TSS distance: 5000
Essential-gene sgRNAs: 40
TSS range: -262 to 436
Shared genes: 1000
CRISPRi essential: 9
CRISPRko essential: 36
Overlap: 8
CRISPRko-only genes: 28
Mean log2 CN of CRISPRko-only genes: 0.044
Genome-wide mean log2 CN: 0.002
```

这几个数字已经把本篇最核心的事实交代清楚了：

- CRISPRi 数据里所有 5000 条 sgRNA 都带有 TSS 距离信息
- 我们抽出的 essential-gene sgRNA 一共有 40 条，足够看出一个稳定趋势
- CRISPRi 与 CRISPRko 在这套教学结果里有 **8 个共同 essential**，但也并不完全重合
- CRISPRko-only 命中的平均 copy number 略高于全局平均，这正是 DSB 偏差要警惕的方向

---

## Step 1：CRISPRi 的第一原则——TSS 距离很重要

CRISPRi 不切 DNA，它依赖 dCas9-KRAB 一类抑制复合物压制转录。因此 guide 的位置比 CRISPRko 更敏感：

- 离 TSS 太远，抑制效率会迅速下降
- 过了最佳窗口，很多 sgRNA 即使设计得很好也打不出强表型

这次脚本里直接从 sgRNA 名称中解析了 TSS 偏移，并计算：

- `lfc = log2((T12_rep1 + T12_rep2 + 1) / (T0_rep1 + T0_rep2 + 1))`

### 图 1：TSS 距离与 CRISPRi 效率

![图 1：TSS 距离](assets/figures/pub_tss_distance.png)

这张图是本篇最重要的一张。阴影区域标出的是常用经验窗口（约 `-50` 到 `+300 bp`）。真实结果能看到：

- 越接近这个窗口，sgRNA 的负向 LFC 越明显
- 偏离窗口后，信号会明显变弱

这也是为什么 CRISPRi 文库设计和 CRISPRko 完全不是一个思路：**CRISPRi 的 guide 不是“只要能切就行”，而是“必须打在正确的调控窗口”。**

---

## Step 2：CRISPRi vs CRISPRko——重合，但不完全相同

同一个基因被判为 essential，并不意味着 CRISPRi 和 CRISPRko 一定会给出相同效应量。脚本这一步把两份结果并到一起：

- `results/crispri_test.gene_summary.txt`
- `results/mageck_test.gene_summary.txt`

并按 `neg|fdr < 0.05` 划分：

- Both
- CRISPRi only
- CRISPRko only
- NS

### 图 2：CRISPRi vs CRISPRko effect size

![图 2：CRISPRi vs CRISPRko](assets/figures/pub_crispri_vs_ko.png)

从这张图看，有三条结论最实用：

1. **大方向是一致的。** 真正的核心 essential genes 往往在两种系统里都往负方向走。
2. **CRISPRko 的动态范围通常更大。** 敲除是“硬失活”，CRISPRi 更像“压低表达”。
3. **CRISPRko-only 命中要特别小心。** 有些基因在 ko 里显著，但在 i 里不显著，可能是生物学差异，也可能部分掺了 copy-number / DSB 偏差。

这也是脚本最后专门去看 copy number 的原因。

---

## Step 3：为什么要盯着 CRISPRko-only hit 的 copy number

脚本把 `ko_fdr < 0.05 且 i_fdr >= 0.05` 的基因单独拿出来，并和本地 CN 表做了合并：

- `CRISPRko-only genes: 28`
- `Mean log2 CN of CRISPRko-only genes: 0.044`
- `Genome-wide mean log2 CN: 0.002`

这不是一个“惊天动地”的差异，但方向上很合理：

- ko-only hits 稍微偏向更高 CN 区域
- 说明 CRISPRko 的部分特异命中，确实可能受到额外 DNA damage 影响

所以在真实项目里，如果一个 hit：

- 只在 CRISPRko 显著
- 在 CRISPRi 完全不见
- 又落在放大区域

那就非常值得多看一眼，而不是立刻把它当作生物学强发现。

---

## Step 4：CRISPRa 的信号通常更稀疏

脚本最后还读了：

- `data/crispra/count_table.txt`

并计算：

- `pos_lfc = log2((Drug_rep1 + Drug_rep2 + 1) / (DMSO_rep1 + DMSO_rep2 + 1))`

### 图 3：CRISPRa positive-selection signal

![图 3：CRISPRa 分布](assets/figures/pub_crispra_dist.png)

这张分布图会让人很快建立一个正确直觉：

- 大多数 sgRNA 都堆在 0 附近
- 真正向右偏得很厉害的 guide 只是少数

也就是说，**CRISPRa 的“有效阳性”往往比 CRISPRko 的 dropout 更稀疏。**

原因很简单：

- 不是所有基因被激活都会产生强表型
- 即使会，也不是所有 promoter/TSS 区域都同样容易被激活
- 激活系统本身的上限和背景噪声也更复杂

所以做 CRISPRa 时，一个非常常见的误区就是：用 CRISPRko 的预期去要求 CRISPRa。实际上这两类筛选的信号密度本来就不一样。

---

## 看一眼 CRISPRi 的 gene-level 结果

```bash
head -5 repro/results/crispri_test.gene_summary.txt | sed 's/\t/    /g'
```

```
📊 输出：
id    num    neg|score    neg|p-value    neg|fdr    neg|rank    neg|goodsgrna    neg|lfc
C1orf109    5    1.5566e-07    4.9505e-06    0.00165    1    3    -2.0368
C7orf26     5    4.8931e-06    4.9505e-06    0.00165    2    3    -2.2519
ANAPC7      5    5.8429e-06    4.9505e-06    0.00165    3    4    -0.98229
ANAPC4      5    1.5466e-05    1.4851e-05    0.003713   4    4    -1.3130
```

前几位像 `ANAPC7`、`ANAPC4` 这类细胞周期 / 蛋白降解复合物相关基因，方向是完全合理的。它们也帮助我们确认：

- CRISPRi 这套教学结果并不是随机噪声
- 它仍然能抓到非常经典的细胞存活依赖因子

---

## 本篇关键输出文件

```bash
du -h \
  repro/results/crispri_test.gene_summary.txt \
  repro/results/crispra_test.gene_summary.txt \
  repro/results/figures/pub_tss_distance.png \
  repro/results/figures/pub_crispri_vs_ko.png \
  repro/results/figures/pub_crispra_dist.png
```

```
📊 输出：
56K   repro/results/figures/pub_crispra_dist.png
80K   repro/results/crispri_test.gene_summary.txt
88K   repro/results/crispra_test.gene_summary.txt
176K  repro/results/figures/pub_crispri_vs_ko.png
188K  repro/results/figures/pub_tss_distance.png
```

---

## 本篇小结

这一篇真正想让你建立的是三个判断标准：

1. **CRISPRi 先看 TSS 距离。** 没有这个窗口概念，很多 guide 好坏根本解释不清。
2. **CRISPRi 和 CRISPRko 要对照着看。** 重合越多，说明你抓到的是真正稳的 dependency；分歧越大，越要看机制和偏差来源。
3. **CRISPRa 的强阳性天然稀疏。** 不要拿 CRISPRko dropout 的信号密度去要求它。

如果第 1 篇解决的是“怎么找 essential gene”，那第 4 篇解决的是：

> “当你不再直接切 DNA 时，信号应该怎么重新解读。”

---

## FAQ：常见问题

**Q1：CRISPRi 最佳窗口一定是 `-50` 到 `+300 bp` 吗？**

不是绝对值，但这是一个非常常见、很有用的经验窗口。不同文库和不同启动子结构会有偏差。

**Q2：CRISPRko-only 命中一定是假阳性吗？**

不一定。有些确实是机制差异，有些则可能掺了 copy-number / DSB 偏差。要结合 CN、sgRNA 一致性和功能背景一起看。

**Q3：为什么本篇不再写外部 GEO accession？**

因为这一版内容明确基于仓库内已经实跑通过的教学数据和结果。写清楚本地来源，比继续挂错号更重要。

---

## 本系列导航

| 篇目 | 主题 | 定位 |
|---|---|---|
| 第 1 篇 | MAGeCK 分析——从 sgRNA 计数到必需基因 | 基础篇 |
| 第 2 篇 | MAGeCK MLE + VISPR——复杂实验设计与交互可视化 | 进阶篇 |
| 第 3 篇 | MAGeCKFlute 整合分析——基因筛选的全景图 | 整合篇 |
| **第 4 篇** | **CRISPRi/CRISPRa 筛选分析策略——不切 DNA 的基因扰动** | **📍 当前阅读** |
| 第 5 篇 | 药物-基因互作筛选与合成致死分析——一加一大于二 | 应用篇 |
| 第 6 篇 | 发表级图表与审稿人常见问题——最后一公里 | 收官篇 |
