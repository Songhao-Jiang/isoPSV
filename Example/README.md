# Minimal Demo / 最小示例

A self-contained demo that runs the full isoPSV pipeline (steps 1-3) with **no
search engine and no GPU required**: one real LC-MS/MS fraction, three real
pFind result tables, and the reference proteome are all bundled under
`Example/data`. Expected outputs are in `Example/expected_outputs/`.

- Spectra: *Chlamydomonas reinhardtii*, pH5 strong-cation-exchange fraction,
  3,784 MS/MS spectra (unlabeled channel). Source raw files: **PRIDE
  PXD003236** (not bundled - download from PRIDE if you need the raws).
- Search results: the original pFind search plus the open-search and
  restricted-search arms from the same project (pH5 rows, real production
  outputs), so the demo exercises genuine cross-strategy competition.
- Reference proteome: NCBI GCF_000002595.2 (*C. reinhardtii* v5.5).

## Quick run

1. Extract the release zip (or clone the repo and use your own R 4.x), copy
   `Example/data` to `D:\isoPSV_demo\data` (any folder works; keep the paths
   below in sync).
2. Copy `configs/demo_config_step1.txt` to `D:\isoPSV_demo\config_step1.txt`.
3. Run step 1 (web console: load the config and Run; or CLI):
   ```bat
   set CONFIG=D:\isoPSV_demo\config_step1.txt
   Rscript code\user_run.R
   ```
   Step 1 generates `config_step2.txt`, then stops at the manual-search point.
4. Append the four lines from `configs/demo_config_step2_additions.txt` to
   `config_step2.txt` (they register the two bundled supplementary
   strategies), and run again.
5. Steps 2-3 complete in ~2 minutes. Compare `result_*\1.unlabeled_result\`
   with `expected_outputs\`.

## Expected results (predictor-off mode)

| item | value |
|---|---|
| spectra processed | 3,784 |
| Stage A: consistent | 3,698 (97.7%) |
| Stage A: candidate misidentification | 86 (2.3%) |
| Part 1 annotation-correction report | 3,784 rows (86 judged wrong) |
| Part 2 novel-peptide report | 0 rows (reference-only database - expected) |
| master table | 3,784 rows x 58 columns |

Row counts are deterministic; per-row scores can differ in the last digits
across R builds. The demo runs without pDeep2 (prediction-based metrics use
defaults); see `Software_Manual.md` for predictor-enabled workflows.

---

# 中文版

自带全部输入、**无需任何搜库引擎与 GPU** 即可跑通 isoPSV 三步流程的最小演示：
一个真实分级馏分的 MGF、三份真实 pFind 结果表与参考蛋白库都在
`Example/data`，预期输出在 `Example/expected_outputs/`。

- 谱图：莱氏衣藻 pH5 离子交换馏分，3,784 张 MS/MS（无标记通道）。原始 raw
  文件＝PRIDE **PXD003236**（不随库分发，需要时去 PRIDE 下载）。
- 搜库结果：同项目的原始搜库＋开放式＋限定式三臂真实产物（pH5 行），演示
  的是真实跨策略竞争。
- 参考库：NCBI GCF_000002595.2（衣藻 v5.5）。

## 快速运行

1. 解压 Release 包（或克隆本仓库并自备 R 4.x），把 `Example/data` 复制到
   `D:\isoPSV_demo\data`（任意目录均可，同步改下述路径）。
2. 把 `configs/demo_config_step1.txt` 复制为 `D:\isoPSV_demo\config_step1.txt`。
3. 跑第一步（网页控制台载入配置点 Run；或命令行）：
   ```bat
   set CONFIG=D:\isoPSV_demo\config_step1.txt
   Rscript code\user_run.R
   ```
   第一步会生成 `config_step2.txt` 并停在手动搜库停点。
4. 把 `configs/demo_config_step2_additions.txt` 的四行追加进
   `config_step2.txt`（登记两份补充策略），再次运行。
5. 第二、三步约 2 分钟完成，将 `result_*\1.unlabeled_result\` 与
   `expected_outputs\` 对照。

## 预期结果（预测器关闭模式）

| 项目 | 数值 |
|---|---|
| 处理谱图 | 3,784 |
| 阶段 A：consistent | 3,698（97.7%） |
| 阶段 A：candidate misidentification | 86（2.3%） |
| Part1 注释校正报告 | 3,784 行（判错 86） |
| Part2 新肽段报告 | 0 行（纯参考库口径，属预期） |
| 主表 | 3,784 行 × 58 列 |

行数判定完全确定；逐行打分在不同 R 版本间末位可能有浮点差异。演示未启用
pDeep2（预测类指标用缺省值），预测器全功能见 `Software_Manual.md`。
