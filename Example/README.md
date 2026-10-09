# Minimal Demo / 最小示例

A complete walk-through on one public dataset: download -> config -> run ->
expected outputs.

---

# 中文版

## 1. 数据

从 PRIDE/iPROX 下载一个小型 Thermo raw 数据集（例如本项目评估所用衣藻
数据集 PRIDE: PXD003236，或任何小规模 raw）。放到 `D:\demo\raw\`。

## 2. 软件准备

- 解压 `isoPSV_v1.0.zip`（Release 下载），双击 `isoPSV.exe` 打开控制台
- pFind Studio 已安装（pParse/Searcher 随之可用；`pfind_exe` 一键检测回填）
- pDeep2/pDeep3 按需安装（[docs/installation.md](installation.md)；不用预测器可跳过）

## 3. 网页填写（对照本目录 `configs/` 示例）

- result_dir＝`D:\demo`
- psm_file＝先用 pFind 对该 raw 搜一次的结果（或任意带 File_Name 列的结果表）
- fasta＝对应物种蛋白质组 FASTA
- **mgf＝直接填 `D:\demo\raw`**（raw 自动经 pParse 转 MGF——无需手工转换）
- 其余全默认 → 参数汇总页生成 → RUN

## 4. 命令行等价（`configs/` 已给最小可用示例）

```bat
set CONFIG=D:\demo\config_step1.txt
Rscript code\user_run.R
```

## 5. 预期产出（result_时间戳\ 下）

- `unlabeled\`：逐 raw 无标记 MGF＋step1_done.txt
- `Competitive_search_result\`＋`merge\`：去诱饵清洗件与合并候选表
- `1.unlabeled_result\`：report_A.tsv、classification_A.tsv、两分法报告
- 标记样本另有 `2.labeled_result\`：成对/冲突报告、master_table.tsv

---

# English

## 1. Data

Download a small Thermo raw dataset from PRIDE/iPROX (e.g. the Chlamydomonas
dataset used in our evaluation, PRIDE: PXD003236, or any small raw set) into
`D:\demo\raw\`.

## 2. Software

- Extract `isoPSV_v1.0.zip` (from Releases), double-click `isoPSV.exe`
- pFind Studio installed (pParse/Searcher come with it; one-click detect
  fills `pfind_exe`)
- pDeep2/pDeep3 optionally ([docs/installation.md](installation.md); skip if
  predictors are not used)

## 3. Web entry (mirrored by the `configs/` examples)

- result_dir = `D:\demo`
- psm_file = a pFind search of that raw (or any table with a File_Name column)
- fasta = the species proteome FASTA
- **mgf = just fill `D:\demo\raw`** (raw auto-converts to MGF via pParse)
- keep everything else default -> generate on the Summary page -> RUN

## 4. CLI equivalent (minimal working example in `configs/`)

```bat
set CONFIG=D:\demo\config_step1.txt
Rscript code\user_run.R
```

## 5. Expected outputs (under result_<timestamp>\)

- `unlabeled\`: per-raw unlabeled MGFs + step1_done.txt
- `Competitive_search_result\` + `merge\`: decoy-cleaned copies and the merged candidate table
- `1.unlabeled_result\`: report_A.tsv, classification_A.tsv, two-part reports
- labeled samples add `2.labeled_result\`: paired/conflict reports, master_table.tsv
