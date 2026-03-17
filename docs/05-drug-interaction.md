---
title: 5. 药物互作与合成致死
---

# CRISPR 筛选最佳实践（五）：药物-基因互作筛选与合成致死分析——一加一大于二

> 📋 教程信息
> - GitHub 仓库：[petemeng/MAGeCK-Tutorial](https://github.com/petemeng/MAGeCK-Tutorial)（完整代码、结果与更新记录）
> - 在线网页：[petemeng.github.io/MAGeCK-Tutorial](https://petemeng.github.io/MAGeCK-Tutorial/)（可点击阅读的网页版教程）
> - 数据性质：仓库内可复现的教学 MLE 结果集，信号结构参考典型 PARP 抑制剂 / DNA repair interaction screen，用于演示 interaction 分析逻辑，不是某个公开 FASTQ 项目的 raw 重跑
> - 教学库规模：`509` 个基因、`2036` 条 sgRNA，属于教学规模，不是全基因组文库
> - 分析对象：drug-gene interaction、synthetic lethal、resistance 与 dose response
> - 本篇重点：`analysis/05_drug_interaction.R` 的实跑输出与图表解释
> - 预计阅读：40 分钟 | 实操：10–20 分钟
> - 难度：⭐⭐⭐⭐⭐（5 星制）
> - 前置知识：完成第 2 篇，理解 design matrix 与 `beta` 的含义

---

## 本篇目标

前面的篇目大多还是“单条件问题”：

- 哪些基因在筛选中 dropout？
- 哪些基因更像 essential？

但药物筛选真正关心的是一个双因素问题：

> 某个基因本身也许不是通用 essential，可它会不会在药物存在时，变成让细胞特别脆弱的“合成致死”靶点？

这就不是简单的 `Drug vs T0`，也不是简单的 `Drug vs DMSO`。你真正要看的，是 **interaction effect**。

本篇直接基于仓库里已经实跑通过的 `drug_mle` 教学结果，回答四件事：

1. synthetic lethal / resistance / common essential 怎么分
2. `drug|beta` 图应该怎么看
3. top hits 是否符合 DNA repair / PARP biology 的常识
4. 剂量升高后，interaction beta 会不会同步增强

---

## 本篇使用的真实入口

分析脚本：`repro/analysis/05_drug_interaction.R`

直接运行：

```bash
cd MAGeCK/repro
Rscript analysis/05_drug_interaction.R
```

```
📊 输出：
Drug MLE genes: 509

Class counts:
Common essential 12
NS 452
Resistance 13
Synthetic lethal 32

Synthetic lethal gene mapping: 31 / 32
Significant GO BP terms: 625
```

这几行已经足够说明这套教学结果是“有结构的”：

- 这套教学 library 本身只有 `509` 个基因，所以这里的 `509` 不是后处理过滤，而是文库规模
- 最终分出 **32 个 synthetic lethal**、**13 个 resistance**、**12 个 common essential**
- 合成致死基因里有 `31 / 32` 能映射到注释库并做 GO BP 富集

换句话说，这一篇不是在假装“全基因组重分析”，而是在一个小而完整的教学数据集上，把 interaction 的判断逻辑跑清楚。

---

## Step 1：interaction MLE 结果长什么样

```bash
head -5 repro/results/drug_mle.gene_summary.txt | sed 's/\t/    /g'
```

```
📊 输出：
Gene    sgRNA    time|beta    time|z    time|p-value    time|fdr    drug|beta    drug|z    drug|p-value    drug|fdr
AAAS    4    0.19113    1.1809    0.1277    0.75581    -0.27164    -1.6561    0.14538    0.86667
AAK1    4   -0.07659   -0.52105   0.57171    0.99607    -0.11693    -0.78578   0.54617    0.99803
AATF    4    0.12261    0.76877   0.35560    0.95378     0.27014     1.6725    0.12967    0.86667
AATK    4    0.06987    0.57618   0.61690    0.99607     0.02014     0.16401   0.92534    0.99803
```

这里最关键的不是 `time|beta`，而是 `drug|beta`：

- `drug|beta < 0`：在药物存在时更容易 dropout，偏向 synthetic lethal
- `drug|beta > 0`：在药物存在时更占优势，偏向 resistance

而这套模型之所以能拆出 `drug|beta`，是因为 `design_matrix.txt` 用的是一个三列设计：`baseline + time + drug`。

- `T0_rep1/2` 只有基线项
- `DMSO_rep1/2` 和 `Drug_rep1/2` 都共享 `time = 1`
- 只有 `Drug_rep1/2` 额外带 `drug = 1`

因此：

- `time|beta` 更像“终点筛选共同带来的变化”
- `drug|beta` 才是“在同样终点条件下，相对 DMSO 的药物特异效应”

这就是为什么 interaction 分析必须用 MLE，而不能只看普通两组差异。

---

## Step 2：先看分类结果

脚本把每个基因按互斥规则分成四类：

- `Synthetic lethal`
- `Resistance`
- `Common essential`
- `NS`

分类是**按顺序分配**的：先判 `Synthetic lethal`，再判 `Resistance`，再从剩余基因里判 `Common essential`，最后才归到 `NS`。

具体规则是：

- `drug|fdr < 0.05` 且 `drug|beta < -0.5` → `Synthetic lethal`
- `drug|fdr < 0.05` 且 `drug|beta > 0.5` → `Resistance`
- `time|fdr < 0.05` 且 `time|beta < -0.5` → `Common essential`
- 其余 → `NS`

这里给 `Common essential` 加上 `time|fdr < 0.05` 很重要：

- 它避免把“方向像 essential、但统计上并不稳”的基因也归进去
- 在当前这套结果里，加不加这一条，`Common essential` 仍然都是 `12` 个，说明这 12 个本来也都统计显著

这和真实项目里常做的第一轮优先级划分非常接近。

---

## Step 3：合成致死和耐药 top hits 是否像真的

脚本打印出的 top hits 非常值得看，因为它们不是“算法上显著但毫无生物学关系”的杂乱基因。

### Top synthetic lethal

```
📊 输出：
XRCC3   -1.55
FANCC   -1.46
FANCD2  -1.44
BRCA2   -1.38
RAD51   -1.36
FANCB   -1.34
RAD51D  -1.28
FANCA   -1.27
RAD51B  -1.26
RBBP8   -1.25
```

### Top resistance

```
📊 输出：
RNF8     1.24
SHLD1    1.21
PARP1    1.20
MAD2L2   1.18
DYNLL1   1.18
FAM35A   1.04
REV7     1.04
RIF1     1.03
SHLD2    0.991
TP53BP1  0.950
```

如果你对 DNA damage repair / PARP biology 熟一点，会立刻发现：

- synthetic lethal 侧聚了很多 **HR / Fanconi anemia / RAD51** 相关基因
- resistance 侧则出现了 **TP53BP1 / REV7 / SHLD / RIF1** 这类和修复路径切换高度相关的因子

这正说明这份教学结果虽然规模不大，但方向是“像真的”：top hits 的信号结构符合典型 PARP 抑制剂筛选里常见的 repair-pathway 逻辑。

---

## Step 4：Differential beta plot 是药物互作里最关键的一张图

### 图 1：Differential beta plot

![图 1：Differential beta](assets/figures/pub_diff_beta.png)

横轴是 `time|beta`，纵轴是 `drug|beta`。脚本在图上额外标了 `BRCA1`、`BRCA2`、`PALB2`、`PARP1`、`TP53BP1`、`ATM`、`RAD51C`、`RPS19`、`RPL11`、`PCNA` 这 10 个代表性基因，方便读图。

这张图的四个象限都值得解释：

- 左下：本身也在掉，在药物里掉得更厉害 → synthetic lethal / strong dependency 最常聚在这里
- 右上：药物下更占优势 → resistance 候选
- 左上：在普通终点里不明显掉，甚至略偏正，但一到药物条件就转为脆弱 → 更像药物特异 sensitizer
- 右下：本身已经在掉，但药物条件下反而没那么掉了 → 更像 drug-buffering / partial rescue，而不是 synthetic lethal

尤其要记住最后一种情况：

> `time|beta` 很负，不等于它就是 synthetic lethal。

如果 `drug|beta` 接近 0，甚至偏正，那它更可能只是 common essential，或者是被药物部分缓冲的基础依赖基因。

这也是为什么 interaction 模型比简单差异更有信息量：

> 它把“本来就 essential”与“药物特异性脆弱”拆开了。

---

## Step 5：Waterfall plot 最适合做 top hit 浏览

### 图 2：Waterfall plot

![图 2：Waterfall](assets/figures/pub_waterfall.png)

waterfall 把所有基因按 `drug|beta` 从最负排到最正：

- 最左端：最强 synthetic lethal
- 最右端：最强 resistance

脚本在这张图上只保留了 `BRCA1`、`BRCA2`、`PALB2`、`PARP1`、`TP53BP1`、`ATM` 这 6 个标签，目的是让图保持干净，同时把最有代表性的 repair / resistance 基因钉在图上。

这种图非常适合：

- 给老板 / 合作者快速看 top list
- 选验证对象
- 看某些经典基因有没有出现在你预期的位置

在这套结果里，`BRCA2`、`RAD51`、`FANCD2` 这些 repair genes 都落在左端，非常合理。

---

## Step 6：合成致死不是 gene list，还是 pathway story

脚本用 `clusterProfiler::enrichGO` 对 synthetic lethal genes 做了 GO Biological Process 富集：

- 输入基因：`32` 个 synthetic lethal genes
- 成功映射：`31 / 32`
- 多重检验：`BH`
- 显著条目：`625` 个 GO BP terms

### 图 3：Synthetic lethal GO BP enrichment

![图 3：合成致死富集](assets/figures/pub_sl_go.png)

这里的 `625` 不能直接理解成“625 条完全不同的通路”，因为 GO 本身有很强的 parent-child 层级冗余。对于这 32 个高度集中在 DNA repair 方向的输入基因来说，大量相近条目会一起显著。

所以这一步真正要看的不是“条目数很多”，而是：

- 富集是不是收束到 DNA repair / homologous recombination / replication stress 相关过程
- 最显著的条目是不是和 top hits 的 biology 对得上

图里展示的是按显著性排序后最有代表性的前 10 个 GO BP 条目；完整结果见 `repro/results/sl_go_bp.tsv`。

也就是说，**通路层面的自洽性，是 synthetic lethal screen 可信度的重要加分项。**

---

## Step 7：剂量效应让 interaction 更有说服力

脚本还额外读了：

- `results/drug_mle_dose.gene_summary.txt`

它对应的是 `design_matrix_dose.txt` 这套四列设计：`baseline + time + drug_lo + drug_hi`。

- `Drug_rep1/2` 对应 `drug_lo`
- `DrugHi_rep1/2` 对应 `drug_hi`
- 仓库里没有把它们绑定到具体 µM 数值，所以这里更适合把它理解成两个药物强度层级，而不是固定浓度

### 图 4：Dose-response of interaction beta

![图 4：剂量效应](assets/figures/pub_dose_response.png)

这张图是 `drug_lo|beta` 对 `drug_hi|beta` 的散点图，脚本额外标了 top 10 synthetic lethal genes。

它最有价值的地方在于，它不只告诉你“某基因在 drug arm 下显著”，还告诉你：

- 随着药物强度变高，这个 interaction effect 是不是同步增强

对于 synthetic lethal 这类 `beta < 0` 的 hit 来说，如果低强度和高强度都保持同方向，而且高强度更负，那它通常会更像真实的药物依赖，而不是噪声。

---

## 本篇关键输出文件

```bash
du -h \
  repro/results/drug_mle.gene_summary.txt \
  repro/results/drug_mle_dose.gene_summary.txt \
  repro/results/synthetic_lethal_genes.tsv \
  repro/results/resistance_genes.tsv \
  repro/results/figures/pub_diff_beta.png \
  repro/results/figures/pub_waterfall.png \
  repro/results/figures/pub_sl_go.png \
  repro/results/figures/pub_dose_response.png
```

```
📊 输出：
4.0K   repro/results/resistance_genes.tsv
4.0K   repro/results/synthetic_lethal_genes.tsv
56K    repro/results/drug_mle.gene_summary.txt
76K    repro/results/drug_mle_dose.gene_summary.txt
144K   repro/results/figures/pub_waterfall.png
172K   repro/results/figures/pub_sl_go.png
220K   repro/results/figures/pub_dose_response.png
264K   repro/results/figures/pub_diff_beta.png
```

---

## 本篇小结

这篇最重要的升级，不是“又多跑了一个 MLE”，而是学会把 interaction effect 单独拿出来看：

1. **`drug|beta` 才是药物特异效应。**
2. **`time|beta` 很负 ≠ synthetic lethal。** 那可能只是 common essential。
3. **top hits 必须和已知 repair biology 对得上。**
4. **剂量效应能显著增强你对 interaction hit 的信心。**

如果第 2 篇讲的是“怎么理解 beta”，那第 5 篇讲的就是：

> “怎么把 beta 解释成真正的药物互作生物学。”

---

## FAQ：常见问题

**Q1：为什么不能直接拿 `Drug vs T0` 做 synthetic lethal？**

因为那会把本来就 essential 的基因和 drug-specific effect 混在一起。

**Q2：为什么 `Drug vs DMSO` 也还不够？**

因为你仍然需要在模型里明确拆开时间效应和药物效应。一个很典型的反例是：某个基因在 `Drug` 和 `DMSO` 两个终点里都各掉了 50%，那么它是 common essential，但 `Drug vs DMSO` 的差异接近 0；这时它不是“对药物没关系”，而是“对两边都重要”。interaction MLE 才能把这类情况和真正的 drug-specific hit 区分开。

**Q3：本篇为什么明确写成“教学 MLE 矩阵”？**

因为这版内容就是基于仓库内本地可复现结果，重点是把逻辑讲对、代码跑通，并在一个小规模教学库里把 synthetic lethal / resistance / common essential 的判断边界讲清楚，而不是冒充论文级 raw reanalysis。

---

## 本系列导航

| 篇目 | 主题 |
|---|---|
| 第 1 篇 | MAGeCK 分析——从 sgRNA 计数到必需基因 |
| 第 2 篇 | MAGeCK MLE + VISPR——复杂实验设计与交互可视化 |
| 第 3 篇 | MAGeCKFlute 整合分析——基因筛选的全景图 |
| 第 4 篇 | CRISPRi/CRISPRa 筛选分析策略——不切 DNA 的基因扰动 |
| **第 5 篇** | **药物-基因互作筛选与合成致死分析——一加一大于二** |
