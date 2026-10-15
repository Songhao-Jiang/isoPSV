# Function Skills / 函数技能（场景锁定调用）

"help only tells you the parameters; a skill constrains how to set them per
scenario." Each block below locks ONE tool call for ONE purpose, with the
exact arguments per scenario.

---

# 中文版

## Skill：跑完整流程

- **唯一工具**：`Rscript code/user_run.R`（环境变量 `CONFIG` 指 config_step1.txt）
- **场景约束**：
  - 常规标记样本：`label_type=N15`＋`all_mgf=<标记 MGF 文件夹>`＋`unlabeled_channel=1`
  - 全标记样本：`unlabeled_channel=0`（不填 all_mgf）
  - 无标记样本：`label_type` 留空；输出＝阶段 A 报告＋无标记主表
  - 小样本（FDR 不可信）：`candidate_scope=top5`
  - 只看新肽段：`proteogenomics_mode=1`；只做注释矫正：`=2`
- **禁止**：直接调用 user_step1/2/3.R（内部驱动，停点语义由总入口管理）

## Skill：raw 转 MGF

- **唯一工具**：第一步 mgf 键直接填 raw（内部 `pParse.exe <cfg>`；配置 GBK＋CRLF）
- **参数**：`input_format=raw`、`output_mgf=1`、`co-elute=1`、`outputpath=<result_dir>/pparse_mgf`
- **注意**：pParse 失败时退出码仍为 0——以 MGF 产物落盘为成功判据

## Skill：非 pFind 结果接入第一步

- **唯一工具**：psm_file 填结果文件＋`search_engine` 指认引擎
- **场景约束**：
  - MaxQuant msms.txt：`search_engine=maxquant`（自动转换，无需其他参数）
  - X!Tandem XML / Mascot / PEAKS CSV：自动识别（Mascot/PEAKS 需 `psm_type=mascot|peaks`）
  - Sage/Comet/MSFragger 结果：**不经第一步**——放 config_step2 的 search_file 槽位
- **禁止**：手改 File_Name 列名（由适配器生成）

## Skill：外部引擎自动搜库

- **唯一工具**：config_step2 写 `<engine>_exe`＋`<engine>_mode=1`＋`search_file=(auto:engine)`
- **参数**：MaxQuant 另需 `search_file_raw`；`auto_search=1` 必开（网页自动置位）
- **场景**：MaxQuant（需 raw）/MSFragger/Comet/Sage/MetaMorpheus（用第一步 MGF）

## Skill：W8 阶段 B 权重

- **唯一工具**：网页 W8 页点 15N/13C/SILAC 模板（写 result_dir 的 config_step3.txt：
  `weights_B`＋`threshold_B`）；管线重建 config_step3 时保留该二键

## Skill：单谱导出

- **唯一工具**：`Rscript code/single_spectrum_export_standalone.R`（自包含单文件）
- **参数**：对话框选 pFind-Filtered.spectra＋MGF；产出 scan<编号>.pdf（16:9）＋300-dpi PNG

---

# English

## Skill: run the full pipeline

- **Only tool**: `Rscript code/user_run.R` (env `CONFIG` -> config_step1.txt)
- **Scenario constraints**:
  - Labeled sample: `label_type=N15` (or C13 / custom rows) + `all_mgf=<labeled MGF folder>` + `unlabeled_channel=1`
  - Fully labeled sample: `unlabeled_channel=0` (no all_mgf)
  - Unlabeled sample: leave `label_type` empty; outputs = stage-A reports + unlabeled master table
  - Small sample (FDR unreliable): `candidate_scope=top5`
  - Novel peptides only: `proteogenomics_mode=1`; annotated correction only: `=2`
- **Forbidden**: calling user_step1/2/3.R directly (internal drivers; stop-point semantics belong to the master entry)

## Skill: raw to MGF

- **Only tool**: fill raw in the step-1 mgf key (internally `pParse.exe <cfg>`; GBK + CRLF config)
- **Parameters**: `input_format=raw`, `output_mgf=1`, `co-elute=1`, `outputpath=<result_dir>/pparse_mgf`
- **Note**: pParse exits 0 even on failure - judge success by the output MGfs

## Skill: feed a non-pFind result at step 1

- **Only tool**: psm_file + `search_engine` naming the engine
- **Scenario constraints**:
  - MaxQuant msms.txt: `search_engine=maxquant` (auto-converted, nothing else needed)
  - X!Tandem XML / Mascot / PEAKS CSV: auto-detected (Mascot/PEAKS need `psm_type=mascot|peaks`)
  - Sage/Comet/MSFragger results: NOT via step 1 - put them into config_step2 search_file slots
- **Forbidden**: hand-editing File_Name column names (adapters generate them)

## Skill: external-engine auto search (STEP 2)

- **Only tool**: config_step2 `<engine>_exe` + `<engine>_mode=1` + `search_file=(auto:engine)`
- **Parameters**: MaxQuant additionally needs `search_file_raw`; `auto_search=1` must be on (the web shell sets it)
- **Scenarios**: MaxQuant (raw) / MSFragger / Comet / Sage / MetaMorpheus (reuse the step-1 MGF)

## Skill: W8 stage-B weights

- **Only tool**: the web W8 page templates (write `weights_B` + `threshold_B` into result_dir/config_step3.txt); pipeline regeneration preserves both keys

## Skill: single-spectrum export

- **Only tool**: `Rscript code/single_spectrum_export_standalone.R` (self-contained single file)
- **Parameters**: pick pFind-Filtered.spectra + the MGF files in the dialogs; outputs scan<number>.pdf (16:9) + 300-dpi PNG
