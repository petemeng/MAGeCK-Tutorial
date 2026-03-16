# MAGeCK 全量版方案

## 1. 目标

当前仓库里的 `MAGeCK/repro/` 是**教程可复现版**，重点是让文章中的代码、图和结果都能在本地快速跑通；它不是论文级的全量原始数据仓库。

全量版的目标不是把几十 GB 到几百 GB 的原始数据直接塞进 Git，而是建立一套：

1. **数据源可追溯**：每篇教程都绑定到经过核对的官方公开数据源
2. **下载可自动化**：原始 FASTQ、supplementary count table、source data 可脚本化获取
3. **分析可重跑**：从下载、计数、RRA/MLE、富集、作图到网页/推文产物都能重建
4. **教程可双轨并存**：保留当前 `repro/` 作为轻量 demo，新增 `full/` 作为全量版工作流

---

## 2. 当前仓库的真实状态

我已经核对过当前仓库，现状是：

- `MAGeCK/repro/data` 只有约 **3.5 MB**
- `MAGeCK/repro/results` 只有约 **13 MB**
- 第 1 篇的 FASTQ 样本只有约 **298 KB/样本**，正文也明确写了“每个样本只有 2,500 reads”
- 第 1、2、3 篇主分析结果只有 **1000 个基因**量级，不是全基因组完整筛选
- 第 4、5 篇包含本地模拟或教学子集，不是完整公开原始 screen

所以“数据很小”不是你看错了，而是这套仓库本来就是**教学版**。

---

## 3. 先修正两个数据源错误

当前文稿里有两条公开数据编号不能继续沿用到全量版：

1. `GSE178354` 不是当前教程声称的人类 CRISPR Brunello/HeLa screen；我核对官方 GEO 页面时，它对应的是**拟南芥 male meiosis/leaf variegation 相关数据**。
   - 官方页：<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178354>
2. `GSE152916` 也不是当前教程声称的 Horlbeck/Wessels CRISPRi K562 screen；我核对官方 GEO 页面时，它对应的是**budding yeast heat shock / translation / phosphoproteome** 研究。
   - 官方页：<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152916>

结论：**全量版必须重新绑定正确的一手公开来源，不能在现有编号上打补丁。**

---

## 4. 我建议的“全量版”定义：分两层做

为了让项目既真实又可维护，我建议把“全量版”拆成两层：

### A. `full-count` 层（优先落地）

特点：
- 用**全基因组 count table / phenotype table / source data**
- 数据已经是论文级规模，不再是 1000-gene demo
- 对 GitHub 友好，不需要先处理超大原始 FASTQ
- 最适合先把第 1–6 篇教程的逻辑改成“真正的大数据版本”

### B. `full-raw` 层（第二阶段）

特点：
- 从 SRA/ENA 下载原始 FASTQ
- 从 `mageck count` 开始重建
- 最接近论文复现，但磁盘和时间成本都明显更高

我的建议是：

- **第一期先做 `full-count`**，把教程主体从教学子集升级成真实全量结果
- **第二期再补 `full-raw`**，至少把第 1 篇的 `count` 流程和第 4 篇的 CRISPRi/CRISPRa 入口做成原始数据可重跑版

---

## 5. 推荐的一手公开数据骨架

### 5.1 统一主骨架：Sanson et al., Nature Communications 2018

这篇文章非常适合做全量版骨架，因为它同时给出了：

- CRISPRko（Brunello / Avana / GeCKOv2）
- CRISPRi（Dolcetto / hCRISPRi-v2）
- CRISPRa（Calabrese / hCRISPRa-v2）
- 多个细胞系（A375、HT29、MelJuSo 等）
- Supplementary Data 1–10
- 原始筛选测序数据的 SRA accession：`SRP172473`

官方来源：
- 文章与补充材料：<https://pmc.ncbi.nlm.nih.gov/articles/PMC7063732/>
- 文中 data availability 明确给出：raw screen sequencing data 在 `SRP172473`，且 genome-wide screen data 在 Supplementary Data 1–10

这篇文章最大的优势是：**它能同时覆盖第 1、2、4 篇的“全量版”需求，并且第 3 篇可以直接建立在其 RRA/MLE 结果之上。**

### 5.2 第 5 篇候选骨架：Tsujino et al., Nature Communications 2023

当前第 5 篇是药物-基因互作 / olaparib 合成致死主题。就“公开论文级 source data”而言，一个可行候选是：

- Tsujino et al. 2023，PARPi / olaparib sensitivity-resistance 相关 CRISPR screen
- 文章提供 supplementary data / source data，可作为**全量 gene-level/source-data 版**候选

官方来源：
- 文章：<https://pmc.ncbi.nlm.nih.gov/articles/PMC10635828/>

但需要强调：

- 我这轮已经核对到它适合作为**full-count / full-source-data 候选**
- 但**raw sgRNA count / FASTQ accession 我还没有在这轮验证到足够可靠的一手入口**
- 所以第 5 篇在第一期应做成“**全量公开结果版**”，第二期再单独找 raw/count-level accession

也就是说：**第 5 篇不建议现在硬改成伪“全量原始重跑版”**，否则很容易再写错 accession。

---

## 6. 每篇教程怎么升级成全量版

| 教程 | 当前状态 | 全量版建议 | 第一阶段可行性 |
|---|---|---|---|
| 1. 基础流程 | 2,500-read toy FASTQ + 1000-gene demo | 换成 `SRP172473` 中的 Brunello negative-selection screen；`full-count` 先用 Supplementary Data 1，`full-raw` 再从 SRA 重跑 `mageck count` | 高 |
| 2. MLE + VISPR | `leukemia.new.csv` 仅 10,000 行，仍非全基因组 | 用 Sanson 2018 的全量 screen count table 组合成多条件 design matrix（如 A375 vs HT29 / library factor） | 高 |
| 3. MAGeCKFlute | 建立在 demo RRA/MLE 上 | 直接接第 1/2 篇的 full-count/full-MLE 输出，并切换到 full DepMap 参考 | 高 |
| 4. CRISPRi/a | 目前混合了教学数据和模拟数据，且 GEO 号错误 | 用 Sanson 2018 Supplementary Data 7–10 做 full CRISPRi/full CRISPRa；若坚持 K562 叙事，再补 Horlbeck/Wessels 正确来源 | 中高 |
| 5. 药物互作 | 当前是模拟数据 | 第一阶段切到 Tsujino 2023 的 full source data / supplementary tables；第二阶段再补 raw/count-level 入口 | 中 |
| 6. 投稿图表 | 由 demo 结果重排而来 | 用 1–5 篇 full outputs 重新生成 Figure_main / summary tables | 高 |

---

## 7. 推荐的新目录结构

建议不要把当前 `MAGeCK/repro/` 直接改烂，而是新建：

```text
MAGeCK/
├── repro/                      # 现有教学版，保留
└── full/
    ├── README.md               # 本方案文档
    ├── manifests/
    │   ├── datasets.tsv        # 官方数据源清单
    │   ├── samples_basic.tsv   # 第1篇 sample sheet
    │   ├── samples_mle.tsv     # 第2篇 sample sheet
    │   ├── samples_crispri.tsv # 第4篇 sample sheet
    │   └── samples_drug.tsv    # 第5篇 sample sheet
    ├── config/
    │   ├── paths.yaml
    │   ├── resources.yaml
    │   └── depmap_release.yaml
    ├── scripts/
    │   ├── download_sra.sh
    │   ├── fetch_supplementary.py
    │   ├── build_design_matrix.py
    │   └── prepare_depmap_reference.R
    ├── workflow/
    │   └── Snakefile
    ├── envs/
    │   ├── mageck.yml
    │   ├── r-flute.yml
    │   └── sra-tools.yml
    ├── raw/                    # 原始 FASTQ，gitignore
    ├── external/               # supplementary/source data，gitignore
    ├── counts/                 # full count table，gitignore
    ├── results/                # full outputs，gitignore 或只保留摘要
    └── reports/
        ├── figures/
        └── summary_tables/
```

关键原则：

- **原始数据不进 Git**
- **下载脚本、sample sheet、design matrix、图表脚本进 Git**
- 如需长期缓存大文件，优先考虑对象存储、NAS、或 DVC/Git LFS，而不是直接推 GitHub

---

## 8. 推荐工作流

### Phase 0：先把 source-of-truth 纠正

目标：把现有教程中错误或模糊的数据源声明全部摘出来，统一改成：

- `demo` 版：保留当前轻量流程
- `full` 版：引用已验证的官方来源

这一步只改文案和 source mapping，不跑大数据。

### Phase 1：第 1 篇 full-count + full-raw 打通

目标：先打通最有教学价值的一篇。

交付物：

- `MAGeCK/full/manifests/samples_basic.tsv`
- `MAGeCK/full/scripts/download_sra.sh`
- `MAGeCK/full/workflow/Snakefile` 中的 `download -> count -> test -> plot`
- 新版 `pub_sgrna_rank_full.png`、`pub_gene_volcano_full.png`

建议先做：

1. 用 Supplementary Data 1 跑 **full-count** 版 `mageck test`
2. 再用 `SRP172473` 选择一组 Brunello triplicate screen 跑 **full-raw** 版 `mageck count`
3. 对比 `count-table 入口` 和 `raw-fastq 入口` 是否能得出结构一致的 top hits

### Phase 2：第 2/3 篇切到 full MLE + Flute

目标：把 MLE 和 Flute 从 toy 表切到真正的全量矩阵。

交付物：

- 多条件 full count table
- 新 design matrix
- `mageck mle` full outputs
- full VISPR demo
- full Flute rankview / squareview / enrichment

这里最重要的是：

- 设计矩阵不再用 `leukemia.new.csv` 这种 10k 行 toy 数据
- 改成一个**真实全基因组多条件 count matrix**

### Phase 3：第 4 篇切到 full CRISPRi/CRISPRa

目标：把现在的“教学 mix 版”改成真正的论文级 CRISPRi/a 数据版。

建议做法：

- 第一版先用 Sanson 2018 的 Supplementary Data 7–10
- 如果后续一定要保留 K562/TSS-window 叙事，再单独加第二套 K562 来源

这一步的重点不是 raw FASTQ，而是：

- full sgRNA phenotype table
- full gene-level phenotype table
- 更真实的 TSS window / phenotype 分布

### Phase 4：第 5 篇改成 full published-result 版

目标：先把“模拟数据”升级成“真实公开药物筛选结果”。

建议分两步：

1. **v1**：用 Tsujino 2023 的 supplementary / source data，做 full published-result 版
2. **v2**：等 raw/count accession 核对完，再升级成 raw/count-level 重跑版

原因很简单：

- 现在第 5 篇是最容易因为 accession 写错而误导读者的一篇
- 先做成真实公开结果版，比继续保留模拟数据更有价值
- 但不建议现在硬编一个 FASTQ/GEO 号进去

### Phase 5：第 6 篇重新出 full publication figures

目标：用真正的全量结果重新组织投稿级图表。

交付物：

- `Figure_main_full.png`
- `Figure_main_full.pdf`
- `statistics_summary_full.tsv`
- full supplementary tables

---

## 9. 资源预算（估算）

以下是我基于公开 screen 规模做的**估算**，不是下载后精确测量值：

- `full-count` 层：通常几十 MB 到几百 MB，完全可控
- `full-raw` 层：建议预留 **100–300 GB** 的本地 scratch 空间
- `fasterq-dump` 解压临时空间通常会显著放大，建议 raw 下载盘和临时盘分开
- CPU 建议至少 8 线程；正式全量 `mageck count` / `mle` 时间会从分钟级上升到小时级

所以全量版的关键不是“能不能写脚本”，而是：

- 数据下载
- 临时空间
- 失败重试
- sample sheet / checksum 管理

这些都要一开始就设计好。

---

## 10. 我建议的实施顺序

如果现在开始做，我建议按下面顺序推进：

1. **先搭 `MAGeCK/full/` 骨架**
2. **优先落第 1 篇 full-count + full-raw**
3. **再切第 2/3 篇到 full MLE/Flute**
4. **第 4 篇切到 full CRISPRi/a supplementary data 版**
5. **第 5 篇先做 full published-result 版，再补 raw/count-level**
6. **最后重出第 6 篇投稿图**

原因：

- 第 1 篇最能证明“这不是 toy data 了”
- 第 2/3 篇最能把分析深度拉起来
- 第 5 篇最容易踩 source 问题，应该晚一点做

---

## 11. 验收标准

我认为全量版至少要满足以下标准：

- 第 1 篇：不再使用 2,500-read toy FASTQ，且 `mageck count` 可从公开 raw data 重建
- 第 2 篇：不再使用 10k 行 `leukemia.new.csv` toy 表
- 第 3 篇：不再使用 `DepMap-demo`
- 第 4 篇：不再使用错误 GEO 号和本地教学混合数据作为主来源
- 第 5 篇：不再以纯模拟数据作为正文主分析来源
- 第 6 篇：所有主图来自 full outputs，而不是 demo outputs

---

## 12. 最终建议

**最稳妥的路线不是“把当前 demo 仓库强行放大”，而是“保留 demo，新增 full 工作流”。**

也就是说：

- `MAGeCK/repro/`：继续服务教程演示、网页发布、微信文章
- `MAGeCK/full/`：服务真正的大数据复现、全量结果、投稿级图表

这是成本最低、也最不容易把现有可运行版本搞坏的方案。

如果你让我继续做，我建议下一步就不是再写文案了，而是直接开始：

1. 建 `datasets.tsv`
2. 建 `samples_basic.tsv`
3. 写 `download_sra.sh`
4. 先把第 1 篇 full raw/count 跑通


---

## 13. 当前已落地的 Phase 1 骨架

本轮已经在仓库里补好了以下文件：

- `MAGeCK/full/manifests/samples_basic.tsv`
  - 已根据 `SRP172473` 的 ENA run metadata 填好第 1 篇推荐原始样本集
  - 当前默认选择的是 `Brunello_mod_tracr` 在 `A375` 的 `pDNA + 3 个 dropout replicates`
  - lane-split 的 run 已在 manifest 中按 `;` 分隔，供下载后合并
- `MAGeCK/full/scripts/download_sra.sh`
  - 可直接读取 manifest 下载 FASTQ，并自动把同一样本的多 run 拼接成一个 merged FASTQ
- `MAGeCK/full/config/paths.yaml`
  - 定义了第 1 篇 full raw/count 的路径、control label、treatment labels
- `MAGeCK/full/workflow/Snakefile`
  - 已写好 `merge_fastq -> mageck count -> mageck test` 的最小工作流骨架
- `MAGeCK/full/envs/mageck.yml`
- `MAGeCK/full/envs/sra-tools.yml`
- `MAGeCK/full/scripts/fetch_brunello_library.py`

### 当前唯一的硬阻塞

第 1 篇的 full-raw 工作流现在还差一个关键输入：

- **完整的 Brunello 文库注释文件**（全量 sgRNA sequence -> gene 映射）

也就是说，`download_sra.sh` 现在已经能用，sample sheet 也已经能用，但 `mageck count` 正式开跑前，还需要把 full Brunello library annotation 放到：

- `MAGeCK/full/external/brunello/library.tsv`

这一步现在已经接入了抓取脚本，但 Addgene 媒体源在当前环境中响应较慢，正式下载时建议单独执行脚本并保留重试。

## 9. 当前已验证到的关键执行细节

- Brunello 全库注释已经成功组装为 `MAGeCK/full/external/brunello/library.tsv`，共 `77,441` 条 sgRNA，覆盖 `19,115` 个基因/分组（包含 `Non-Targeting Control`）。
- 由于当前环境访问 Addgene CDN 经常超时，`MAGeCK/full/scripts/fetch_brunello_library.py` 已改成：优先尝试官方 Addgene，失败后自动切换到 GitHub 镜像，并使用分块 `Range` 下载保证可恢复。
- `SRR8297997` 的已下载 pDNA 快照上，前 `50,000` 条 reads 的最佳匹配出现在正向 `trim-5=23–30`；反向互补匹配几乎为零。
- 基于同一快照的小样本验证，执行 `mageck count` 后共 `50,000` reads 中成功映射 `38,897` 条，因此 article 1 的全量工作流已固定使用 `--trim-5 23,24,25,26,27,28,29,30`。
- 在完整 `BrunelloMod_pDNA` FASTQ 上执行 `mageck count --test-run`，前 `1,000,001` reads 中成功映射 `786,514` 条，未命中的库位点仅 `782` 个，说明 Brunello 注释和 trim 参数都已经与真实原始数据对齐。
- `MAGeCK/full/scripts/download_sra.sh` 已改为优先使用 `aria2c` 断点续传下载 ENA FASTQ；在当前网络条件下明显快于原来的单线程 `curl`。
