# Reproduction Guide / 复现指南

This document explains how to (1) reproduce the bundled minimal demo end to
end, and (2) reproduce the manuscript's *Chlamydomonas reinhardtii* results
from the public raw data. Other evaluation datasets (E. coli NeuCode, etc.)
are listed in the manuscript's Data availability section.

---

## 1. Bundled demo (engine-free, ~2 minutes)

See [Example/README.md](../Example/README.md). Everything needed ships in the
repository/release archive: one pH5 fraction (3,784 MS/MS spectra), the
original pFind search plus the open/restrict arms, and the reference proteome
(NCBI GCF_000002595.2). No search engine or GPU is required - the demo
exercises decoy removal, cross-strategy competition, Stage-A arbitration and
the two-part report on real data.

Expected values (deterministic):

| item | value |
|---|---|
| spectra processed | 3,784 |
| consistent / candidate misidentification | 3,698 (97.7%) / 86 (2.3%) |
| Part 1 annotation-correction report | 3,784 rows, 86 judged wrong |
| Part 2 novel-peptide report | 0 rows (reference-only database) |
| master table | 3,784 rows x 58 columns |

## 2. Full *C. reinhardtii* results from public data

**Raw data**: PRIDE **PXD003236** (*C. reinhardtii*, 15N labeling, SCX
fractionation; strain cw15). Download the raw files and convert to MGF with
pParse/msconvert.

**Reference proteome**: NCBI GCF_000002595.2 (v5.5), bundled with the demo at
`Example/data/`.

**Pipeline** (parameters exactly as in the manuscript Methods; software side
defaults, see `Software_Manual.md` and `docs/parameters.md`):

1. **Original search** - pFind closed search against the reference proteome
   (the `psm_file` input).
2. **Supplementary arms** - pFind open search (`pfind_open`) and restricted
   multi-mod search (`pfind_restrict`); enable with `auto_search=1` +
   `pfind_exe`, or provide external-engine results (`Comet`, `MSFragger`,
   `Sage`, `MaxQuant`, `MetaMorpheus`) via `search_file`.
3. **isoPSV steps 1-3** - unlabeled-spectrum extraction, decoy removal +
   multi-strategy merge, Stage-A competitive arbitration; for labeled
   channels set `label_type=N15` (or C13/SILAC/NeuCode) and fill `all_mgf`
   to enable Stage-B paired-spectrum verification.
4. **Novel events** (manuscript Fig. 6) - search the labeled/unlabeled
   fractions against the six-frame translated genome plus predicted ORFs
   (AUGUSTUS + joingenes, see Methods), then run the Part-2 novel-peptide
   competition; representative manuscript-level targets: 179 novel spectra
   supporting 39 novel events (16 gene-level / 12 isoform-level / 11
   N-terminal extension).

Compute: the demo runs on a laptop; the full dataset (dozens of fractions)
takes hours on a desktop with pFind installed (the open/restrict arms are
the long poles; results are fingerprint-cached across reruns).

---

## 中文版

本文说明（1）端到端复现内置最小演示，（2）从公开原始数据复现正文衣藻
结果。其他评测数据集（E. coli NeuCode 等）见正文 Data availability。

### 1. 内置演示（免引擎，约 2 分钟）

见 [Example/README.md](../Example/README.md)。仓库/发布包自带一个 pH5 馏分
（3,784 张谱）、原始搜库＋open/restrict 三臂结果与参考蛋白库，无需搜库引
擎与 GPU，即可在真实数据上走完去诱饵、跨策略竞争、阶段 A 仲裁与两分法
报告。预期值（确定性）：3,784 谱 → 3,698 consistent（97.7%）/ 86 判错
（2.3%）；Part1 3,784 行（判错 86）；Part2 0 行（纯参考库口径）；主表
3,784 行 × 58 列。

### 2. 衣藻全量复现

- **原始数据**：PRIDE **PXD003236**（衣藻 15N 标记、SCX 分级；cw15 株）。
  下载 raw 后经 pParse/msconvert 转 MGF。
- **参考库**：NCBI GCF_000002595.2（v5.5），演示目录已随包附带。
- **流程**（参数与正文 Methods 一致；软件侧见 `Software_Manual.md` 与
  `docs/parameters.md`）：
  1. 原始搜库——pFind 闭合式搜参考库（作为 `psm_file` 输入）。
  2. 补充臂——pFind 开放式（`pfind_open`）与限定式多参数（`pfind_restrict`），
     `auto_search=1`＋`pfind_exe` 自动执行；外部引擎结果（Comet/MSFragger/
     Sage/MaxQuant/MetaMorpheus）经 `search_file` 登记亦可。
  3. isoPSV 三步——无标记谱提取、去诱饵＋多策略合并、阶段 A 竞争仲裁；
     标记通道设 `label_type=N15`并填 `all_mgf`，
     即启用阶段 B 成对谱图验证。
  4. 新事件（正文图 6）——六框翻译库＋预测 ORF（AUGUSTUS＋joingenes，见
     Methods）搜库后走 Part2 新肽段竞争；正文级参考值：179 张新谱支撑
     39 个新事件（基因层 16 / 异构体层 12 / N 端延伸 11）。
- **算力**：演示级笔记本即可；全量数十个馏分在装有 pFind 的台式机上需
  数小时（open/restrict 臂为长杆，重跑按指纹缓存复用）。
