# Pipeline Reference / 流程参考

## Strategy slots (single-letter system, up to 26 strategies)

Every strategy type is hardcoded to a fixed single-letter slot; the merge
orders the slots alphabetically and only the enabled slots carry columns.

| letter | strategy type | post-processing | denovo look-back |
|---|---|---|---|
| a | original (any engine's converted result) | - | |
| b | pfind_open | full (13C drop + terminal loss + SAP exchange) | |
| c | pfind_restrict | none | |
| d | punifind_rescore | terminal (13C drop + terminal loss) | |
| e | punifind_denovo | none | yes |
| f | maxquant | none | |
| g | comet | none | |
| h | msfragger | none | |
| i | sage | none | |
| j | msfragger_boost (MSBooster rescored MSFragger result; same pepXML/TSV format, same adapter) | none | |
| k | xtandem | none | |
| l | metamorpheus | none | |
| m | peaks | none | |
| n | mascot | none | |
| o | pnovo | none | |
| p | msfragger_open | none | |
| p-z | other (generic slots x12, reusable for any engine) | none | |

- Slot letters appear in the merged table as column names `a.Sequence`,
  `b.Sequence`, ... (letter-leading names are valid R names - no X prefix is
  needed, unlike the digit-leading names of the classic edition).
- The strategy state strings (Peptide_state/Modification_state) use the same
  letters, e.g. "abd;c".
- The master table keeps the candidate-1..5 columns identical to the classic
  96-column layout and appends cand6..candN (15 columns per extra candidate).

## Structure

```
isoPSV/
├── code/                            all R code (base R only)
│   ├── user_run.R                   master driver (the only entry point; checkpoint-resumable)
│   ├── unlabeled_extract.R          step 1: unlabeled-spectrum extraction
│   ├── user_step2.R                 step 2: decoy removal + multi-strategy search +
│   │                                per-strategy post-processing + merge + SAP list +
│   │                                denovo database look-back
│   ├── user_step3.R                 step 3: stage-A per-spectrum evaluation
│   ├── extend_report_A.R            stage-A report extension (SAP annotation)
│   ├── stageA_classify.R            stage-A conclusion classification (auto)
│   ├── stageB_paired_match.R        stage B: paired labeled-spectrum matching (auto)
│   ├── stageB_conflict_report.R     conflict report (auto)
│   ├── master_table_export.R        automatic master-table export at the end of the run:
│   │                                master_table.tsv in the same 96-column format as the
│   │                                master table of the Fig.7 panel-f data; for unlabeled
│   │                                samples the export is skipped
│   ├── spectrum_pdf_export.R        automatic labeled-side PDF export at the end of the
│   │                                run (PDF_ON=1 gate; only the confident
│   │                                original-wrong spectra; per-raw split files)
│   ├── spectrum_pdf_export_A.R      automatic unlabeled-side PDF export (stage-A judgment)
│   ├── single_spectrum_export.R     single-spectrum exporter (library version)
│   ├── single_spectrum_export_standalone.R
│   │                                single-spectrum exporter, fully self-contained single file
│   │                                (embedded mass table, choose.files dialogs, base R only)
│   ├── ReadMGF_script.R             core library (MGF parsing / scoring / merge / look-back /
│   │                                decoy removal / result-folder resolution)
│   └── pFind_templates/             pFind search templates
├── tests/                            unit tests (run from tests/ - expect all PASS)
├── isoPSV.exe                       one-click web-shell launcher (double-click; prefers the
│                                    bundled portable runtimes, falls back to PATH)
├── start_isoPSV.bat                 same launcher as a batch script
├── runtime/                         portable runtimes: runtime\python (official Python 3.12.10
│                                    embeddable) and runtime\R (R 4.5.1, base + recommended
│                                    packages only) - no installation needed
├── web/                             web shell (server.py / index.html / spectrum.js;
│                                    Python stdlib only)
├── docs/                            manual PDFs shipped with the package
├── config_step1_template.txt        step-1 parameter file template
└── Software_Manual.md               software manual (pDeep2 install, normalization, web shell)
```

## Run the pipeline

```
Rscript code/user_run.R          (CONFIG env var points to config_step1.txt)
```

- Step 1 extracts the unlabeled spectra; step 2 removes decoy hits from every
  search result (rows whose protein list is all-REV are dropped; mixed lists
  keep only the target proteins - the cleaned copies go into
  Competitive_search_result/ and are the sole input downstream), then runs the
  multi-strategy competitive search (built-in automation for pFind
  restricted/open + pUniFind rescore/denovo; results of other engines such as
  MaxQuant plug in through the unified interchange format - a minimal
  File_Name/Sequence/Modification table, with optional Proteins/Positions) with
  the per-strategy post-processing (restricted: none; open: 13C drop +
  terminal-loss normalization + SAP exchange; rescore: terminal-loss
  normalization only; denovo: database look-back - the peptide is matched
  against the transformed database after I->L and N->GG substitution, matches
  are replaced by the original database sequence with Proteins/Positions
  annotation, non-matches are dropped); step 3 runs stage A and stage B,
  writes the conflict report, exports the master table (labeled samples) and
  automatically exports the spectrum PDFs (both the unlabeled-side and the
  labeled-side two-PDF sets). Engine scores are never used: all judgment
  metrics are recomputed from the spectra and pDeep2 predictions.
- Result folder: every pipeline output goes into a timestamped result folder
  `result_YYYYMMDDHHMMSS\` created next to the parameter file, with four
  subfolders: Competitive_search_result\ (decoy-cleaned search results + SAP
  index), merge\ (revised copies, 1.combine_result.spectra, peptide.txt,
  predict.txt, SAP index), 1.unlabeled_result\ (stage A) and
  2.labeled_result\ (stage B, master table, PDFs). An existing result_* folder
  means the newest one is reused (checkpoint resume in place); removing all
  result_* folders and rerunning creates a fresh timestamped one. A stale
  config_step3.txt (pointing into a removed result folder) is rebuilt
  automatically.
- Strategy typing: the optional parallel key `search_file_type` in
  config_step2.txt names the strategy type of each search_file entry (see the
  slot table above); when omitted, the file/folder naming rules infer it.
  Duplicate fixed types are an error; use "other" repeatedly to fill the
  p..z generic slots.
- Engine adapters: code/search_adapter_maxquant.R converts MaxQuant msms.txt
  into the interchange format (mod-name mapping, raw+scan -> spectrum-title
  mapping, decoy/contaminant/PEP filtering); further engines follow the same
  pattern. An engine's own scores are never used by the judgment.
- Dependencies: base R only (pFind / pUniFind / pDeep2 are external search
  tools; install them per Software_Manual.md).
- Web shell: double-click `isoPSV.exe` (or `start_isoPSV.bat`) at the package
  root. The launcher prefers the bundled portable runtimes (`runtime\python`,
  `runtime\R`) and falls back to PATH, so machines without Python/R installed
  can still open and run the software.
- Self-checks: run tests/run_all.bat (or Rscript tests/run_all.R) - the
  unified entry auto-discovers every `*_test.R` under tests/ and writes
  PASS/FAIL to tests/last_run.log. The test scripts locate code/ relative to
  themselves; worker clusters bootstrap from the probed code directory and the
  driver scripts load the core libraries via boot.R, so the pipeline does not
  depend on being launched from code/.
- Master table scope: the full-spectrum master table (master_table.tsv) is
  available for labeled samples only; for unlabeled samples use the stage-A
  report_A.tsv and classification_A.tsv.
- Version: the root VERSION file records the build version; user_run.R prints
  it at startup for log self-attestation.
- Windows note: with R 4.5 on paths containing Chinese characters, launch
  from a shell without a UTF-8 locale (e.g. cmd.exe), or via
  `env -u LANG -u LC_ALL -u LC_CTYPE` under Git Bash/MSYS; user_run.R has a
  built-in normalizePath fallback, but other tools may still be affected.

## Engine environment prerequisites

- **MetaMorpheus**: needs the **.NET 10 runtime** (`winget install
  Microsoft.DotNet.Runtime.10`) AND the **official full command-line package**
  (the zip containing `Data/`, `Mods/`; a partial extraction with only
  exe/dll does NOT work). Point `metamorpheus_exe` at the `CMD.exe` INSIDE the
  full package.
- **Engine workspace and project directories must be ASCII paths**; clear
  `LC_CTYPE`/`LANG`/`LC_ALL` before starting R (otherwise `list.files` and the
  pUniFind argument assembly fail under Chinese paths).
- **Built-in MGF→mzML converter**: every MS2 spectrum carries
  `MS:1000127 centroid spectrum` (required by Sage - without it Sage panics
  with "contains profile data", exit 1033 - and by MetaMorpheus' MS2
  counting). When `msconvert_exe` is not configured, Sage/MSFragger/
  MetaMorpheus run on the built-in converter output. The converter is a
  vectorized implementation (bulk spectrum-boundary parsing plus a table-lookup
  base64 encoder), and its output is cached at project level (`mzml_cache\`):
  when several engines consume the same MGFs, the file converts once and later
  slots copy; an output newer than its source skips at zero cost.
- **Decoy cleaning is a line-level filter**: non-decoy rows pass through
  verbatim and only decoy-carrying rows are parsed, so million-row open-search
  tables clean in minutes with byte-identical output to the reference
  implementation (see tests/decoy_clean_fast_test.R).
- **label_type class-2 labels are user-defined modification rows** in the
  `name|position|aa-set|composition` grammar (positions ANY/PEP_N/PEP_C/
  PRO_N/PRO_C; isotope elements 2H/13C/15N/18O/34S with signed counts;
  modification.ini style, e.g. `NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)`),
  alongside the N15/C13 full-atom class and the legacy preset names - see
  Software_Manual.md.
- **pUniFind rescore/denovo run per single MGF (per raw)**: multi-raw folders
  are chunked automatically (`chunks\<raw>\` with per-chunk done markers that
  also store an input fingerprint, so a re-search invalidates stale chunks);
  qry.res is copied into each chunk (in-chunk sanitize/renumber rewrites in
  place and must not touch the project originals), the read-only MGF stays
  hardlinked; memory scales with one raw's spectra, a failed chunk resumes on
  rerun, denovo merges with a row-count conservation check and rescore merges
  only column-identical tables. The chunked run exposes three GPU-supply knobs
  in the config: `punifind_batch_size` (inference batch, default 64, try
  128-256 on large sets), `punifind_chunk_parallel` (concurrent chunks,
  default 1; 2-3 fits a 16 GB card - each chunk is a separate ~3 GB process
  with a staggered MASTER_PORT, failures are collected across all chunks
  before stopping) and `punifind_chunk_cleanup` (drop per-chunk data after a
  successful merge, keeping done markers and merged products). Chunk failures
  embed the child log tail directly in the error message.

## Master table

- Pipeline end (automatic): code/master_table_export.R writes
  master_table.tsv in the same 96-column format as the Fig.7 panel-f master
  table (candidate-1..5 parameters, judgment conclusion, protein mapping,
  difference type, per-candidate specific-ion ratios). Labeled samples only:
  for unlabeled samples the export is skipped; the unlabeled needs are fully
  covered by report_A.tsv and classification_A.tsv.

## Spectrum PDF export

- Gated off by default: the PDF stages run only when the environment variable
  `PDF_ON=1` is set (the pages are review material, not a pipeline product);
  `PDF_LIMIT=N` caps the page count (self-checks).
- Labeled side: code/spectrum_pdf_export.R writes the category
  error_match_spectra into 2.labeled_result/ (one PDF per raw,
  `<category>__<raw>.pdf`, sliced into __part2/3... past
  pdf_pages_per_file): only the confident "original result wrong" spectra -
  and novel_peptide_spectra__*.pdf (independent of labeling): every credible
  novel peptide of the two-part report (credible_novel = TRUE), top panel =
  the novel peptide (corrected_*), bottom panel = its closest annotated
  competitor. Each page: candidate-1 spectrum + its ppm panel, winner
  spectrum + its ppm panel, info lines (Type / spectrum id / scores /
  modification change / judgment + conclusion), ion ladder, Base peak and
  Mass error annotations; protein and flanking residues come from the master
  table (cand1_*/corrected_*).
- Unlabeled side: code/spectrum_pdf_export_A.R writes the stage-A
  error_match_spectra__*.pdf (candidate misidentifications), same gate, plus
  novel_peptide_spectra__*.pdf: credible novel peptides are independent of
  labeling - the unlabeled two-part report computes the same six doubt rules
  from the stage-A per-spectrum candidates and the step-1 unlabeled spectra,
  and the PDF renders each credible row (top = the novel peptide, bottom =
  its strongest annotated competitor).
- Storage split: all three PDF categories are written per raw
  (`<category>__<raw>.pdf`), sliced into `__part2/3...` whenever a group
  exceeds `pdf_pages_per_file` (config_step1 key, default 1000 pages); the
  per-spectrum index `spectrum_pdf_index.tsv` (File_Name -> category / file /
  page) lands next to the PDFs in 1.unlabeled_result/ and 2.labeled_result/.
- Result browsing: the web Results panel can load a whole project folder
  copied from another computer (pick the folder, Load): reports, spectra and
  PDF cards are then read from that folder only - no rerun needed. Clear
  switches back to this machine's current run.
- Single spectra: code/single_spectrum_export.R (library) or the standalone
  single_spectrum_export_standalone.R (self-contained, portable) - choose
  pFind-Filtered.spectra and the MGF files in the dialogs; outputs one
  scan<number>.pdf (16:9) + 300-dpi PNG per spectrum, split per MGF into
  subfolders named by the MGF stem.

---

# isoPSV 多引擎竞争验证流程

## 策略槽位（单字母制，至多 26 个策略）

每种策略类型固定占用一个单字母槽位；合并表按字母序排列，仅启用的槽位带列。

| 字母 | 策略类型 | 后处理 | 从头测序回帖 |
|---|---|---|---|
| a | 原始结果（任意引擎转换结果） | - | |
| b | pfind_open | 完整（13C 剔除＋末端丢失归一＋SAP 交换） | |
| c | pfind_restrict | 无 | |
| d | punifind_rescore | 末端（13C 剔除＋末端丢失归一） | |
| e | punifind_denovo | 无 | 是 |
| f | maxquant | 无 | |
| g | comet | 无 | |
| h | msfragger | 无 | |
| i | sage | 无 | |
| j | msfragger_boost（MSBooster 重打分的 MSFragger 结果；pepXML/TSV 格式与适配器同 h） | 无 | |
| k | xtandem | 无 | |
| l | metamorpheus | 无 | |
| m | peaks | 无 | |
| n | mascot | 无 | |
| o | pnovo | 无 | |
| p | msfragger_open | 无 | |
| p-z | other（通用槽位 ×12，可复用于任意引擎） | 无 | |

- 槽位字母即合并表列名 `a.Sequence`、`b.Sequence`……（字母开头的列名是合法
  R 名称，不像经典版的数字开头列名需要 X 前缀）。
- 策略状态串（Peptide_state/Modification_state）使用相同字母，如 "abd;c"。
- 主表的候选 1..5 列与经典 96 列布局完全一致，其后追加 cand6..candN
  （每个额外候选 15 列）。

## 目录结构

```
isoPSV/
├── code/                            全部 R 代码（仅基础包）
│   ├── user_run.R                   总入口（唯一入口；断点续跑）
│   ├── unlabeled_extract.R          第一步：无标记谱图提取
│   ├── user_step2.R                 第二步：去诱饵＋多策略搜库＋逐策略后处理＋
│   │                                合并＋SAP 清单＋从头测序回帖
│   ├── user_step3.R                 第三步：阶段 A 逐谱评估
│   ├── extend_report_A.R            阶段 A 报告扩充（SAP 注释）
│   ├── stageA_classify.R            阶段 A 结论分类（自动）
│   ├── stageB_paired_match.R        阶段 B：成对标记谱匹配（自动）
│   ├── stageB_conflict_report.R     冲突报告（自动）
│   ├── master_table_export.R        运行结束自动导出主表：
│   │                                master_table.tsv 与图 7f 面板主表同 96 列格式；
│   │                                无标记样本跳过导出
│   ├── spectrum_pdf_export.R        运行结束自动导出标记侧 PDF
│   │                                （PDF_ON=1 门控；仅可信的原始判错谱图；按 raw 拆分）
│   ├── spectrum_pdf_export_A.R      自动导出无标记侧 PDF（阶段 A 判定）
│   ├── single_spectrum_export.R     单谱导出器（库版本）
│   ├── single_spectrum_export_standalone.R
│   │                                单谱导出器，完全自包含单文件
│   │                                （内嵌质量表、choose.files 对话框、仅基础 R）
│   ├── ReadMGF_script.R             核心库（MGF 解析/打分/合并/回帖/去诱饵/
│   │                                结果目录解析）
│   └── pFind_templates/             pFind 搜库模板
├── tests/                           单元测试（在 tests/ 内运行——应全部 PASS）
├── isoPSV.exe                       一键网页壳启动器（双击；优先自带绿色运行时，
│                                    回退 PATH）
├── start_isoPSV.bat                 同一启动器的批处理版本
├── runtime/                         绿色运行时：runtime\python（官方 Python 3.12.10
│                                    embeddable）与 runtime\R（R 4.5.1，基础＋推荐包）
│                                    ——免安装
├── web/                             网页壳（server.py / index.html / spectrum.js；
│                                    仅 Python 标准库）
├── docs/                            随包手册 PDF
├── config_step1_template.txt        第一步参数文件模板
└── Software_Manual.md               软件手册（pDeep2 安装、归一化、网页壳）
```

## 运行流程

```
Rscript code/user_run.R          （CONFIG 环境变量指向 config_step1.txt）
```

- 第一步提取无标记谱图；第二步对各搜库结果去诱饵（蛋白列表全 REV 的整行
  丢弃；混合列表只保留正库蛋白——清洗件落 Competitive_search_result/，是下游
  唯一输入），随后执行多策略竞争搜库（内置 pFind 限定式/开放式＋pUniFind
  rescore/denovo 自动化；MaxQuant 等其他引擎结果经统一交换格式接入——最小
  File_Name/Sequence/Modification 表，可选 Proteins/Positions），并做逐策略
  后处理（限定式：无；开放式：13C 剔除＋末端丢失归一＋SAP 交换；rescore：
  仅末端丢失归一；denovo：数据库回帖——肽段在 I->L 与 N->GG 替换后的转换库
  中匹配，命中替换回原库序列并带 Proteins/Positions 注释，未命中丢弃）；
  第三步跑阶段 A 与阶段 B、写冲突报告、导出主表（标记样本）并自动导出谱图
  PDF（无标记侧与标记侧两套）。引擎打分从不参与判定：全部判定指标由谱图
  与 pDeep2 预测重新计算。
- 结果目录：管线全部输出写入参数文件旁的时间戳目录
  `result_YYYYMMDDHHMMSS\`，含四个子目录：Competitive_search_result\（去诱饵
  清洗件＋SAP 索引）、merge\（修订件、1.combine_result.spectra、peptide.txt、
  predict.txt、SAP 索引）、1.unlabeled_result\（阶段 A）与 2.labeled_result\
  （阶段 B、主表、PDF）。已有 result_* 目录时复用最新一个（断点续跑）；删除
  全部 result_* 再跑会新建时间戳目录。指向已删目录的过期 config_step3.txt
  会自动重建。
- 策略类型：config_step2.txt 的可选并行键 `search_file_type` 逐条指明
  search_file 的策略类型（见上表）；缺省时按文件/文件夹命名规则推断。固定
  类型重复为错误；可反复用 "other" 填充 p..z 通用槽位。
- 引擎适配器：code/search_adapter_maxquant.R 把 MaxQuant msms.txt 转成交换
  格式（修饰名映射、raw+scan → 谱图标题映射、诱饵/污染物/PEP 过滤）；其他
  引擎同模式扩展。引擎自身打分从不参与判定。
- 依赖：仅基础 R（pFind / pUniFind / pDeep2 为外部搜库工具；按
  Software_Manual.md 安装）。
- 网页壳：双击包根 `isoPSV.exe`（或 `start_isoPSV.bat`）。启动器优先使用
  自带绿色运行时（`runtime\python`、`runtime\R`）并回退 PATH，未装 Python/R
  的电脑也能直接打开运行。
- 自检：运行 tests/run_all.bat（或 Rscript tests/run_all.R）——统一入口自动
  发现 tests/ 下全部 `*_test.R` 并把 PASS/FAIL 写入 tests/last_run.log。测试
  脚本相对自身定位 code/；worker 集群从探测到的代码目录自举，驱动脚本经
  boot.R 加载核心库，管线不依赖从 code/ 启动。
- 主表范围：全谱主表（master_table.tsv）仅标记样本输出；无标记样本用阶段 A
  的 report_A.tsv 与 classification_A.tsv。
- 版本：包根 VERSION 文件记录构建版本；user_run.R 启动时打印用于日志自证。
- Windows 注意：R 4.5 在含中文字符的路径下，请从无 UTF-8 locale 的终端
  （如 cmd.exe）启动，或在 Git Bash/MSYS 下用
  `env -u LANG -u LC_ALL -u LC_CTYPE`；user_run.R 有内建 normalizePath
  兜底，但其他工具仍可能受影响。

## 引擎环境前提（工作站实测归纳）

- **MetaMorpheus**：需 **.NET 10 运行时**（`winget install
  Microsoft.DotNet.Runtime.10`）且须**官方完整命令行包**（含 `Data/`、
  `Mods/` 的 zip；仅 exe/dll 的部分解压不可用）。`metamorpheus_exe` 应指向
  完整包内的 `CMD.exe`。
- **引擎工作区与项目目录必须 ASCII 路径**；启动 R 前清除
  `LC_CTYPE`/`LANG`/`LC_ALL`（否则中文路径下 `list.files` 与 pUniFind 参数
  组装会失败）。
- **内置 MGF→mzML 转换器**：每张 MS2 谱带 `MS:1000127 centroid spectrum`
  （Sage 必需——缺它报 "contains profile data"、退出码 1033；
  MetaMorpheus 的 MS2 计数同样需要）。未配置 `msconvert_exe` 时 Sage/
  MSFragger/MetaMorpheus 走内置转换器产物。转换器为向量化实现（整块谱边界
  解析＋查表式 base64 编码器），产物按项目级缓存（`mzml_cache\`）：多引擎
  消费同一批 MGF 时只转换一次、后续槽位复制；产物比源新则零成本跳过。
- **诱饵清洗是行级过滤**：非诱饵行原样直通、仅含诱饵的行才解析，百万行开放
  式搜库表数分钟清洗完成，与参考实现逐字节一致（见
  tests/decoy_clean_fast_test.R）。
- **label_type 第二类标记＝用户自定义修饰行**，语法
  `名称|位置|氨基酸集合|元素组成`（位置 ANY/PEP_N/PEP_C/PRO_N/PRO_C；同位素
  元素 2H/13C/15N/18O/34S 可带符号计数；modification.ini 风格，如
  `NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)`），与 N15/C13 全原子类及旧预设名
  并存——详见 Software_Manual.md。
- **pUniFind rescore/denovo 逐单 MGF（逐 raw）运行**：多 raw 文件夹自动分块
  （`chunks\<raw>\`，每块带完成标记并保存输入指纹——重新搜库会使过期分块
  失效）；qry.res 拷入各分块（块内 sanitize/重编号原地改写、不碰项目原件），
  只读 MGF 保持硬链接；内存按单 raw 谱量伸缩，失败分块重跑即续，denovo 合并
  带行数守恒校验、rescore 仅合并列一致表。分块运行在 config 暴露三个 GPU
  供给旋钮：`punifind_batch_size`（推理批，默认 64，大数据集可试 128-256）、
  `punifind_chunk_parallel`（并发分块数，默认 1；16 GB 显卡可 2-3——每块是
  约 3 GB 的独立进程、MASTER_PORT 错峰，全部分块失败收集齐才停）与
  `punifind_chunk_cleanup`（合并成功后丢弃分块数据，保留完成标记与合并产物）。
  分块失败把子进程日志尾部直接嵌进报错信息。

## 主表

- 管线收尾（自动）：code/master_table_export.R 按图 7f 面板主表同款 96 列格式
  写 master_table.tsv（候选 1..5 参数、判定结论、蛋白回贴、差异类型、逐候选
  特异子离子比例）。仅标记样本；无标记样本跳过导出——report_A.tsv 与
  classification_A.tsv 已完整覆盖其需求。

## 谱图 PDF 导出

- 默认关闭：仅当环境变量 `PDF_ON=1` 时 PDF 阶段运行（页面属复核材料而非管线
  产物）；`PDF_LIMIT=N` 限定页数（自检用）。
- 标记侧：code/spectrum_pdf_export.R 把 error_match_spectra 类目写入
  2.labeled_result/（每 raw 一个 PDF，`<类目>__<raw>.pdf`，超
  pdf_pages_per_file 拆 __part2/3…）：仅可信的"原始判错"谱图——以及
  novel_peptide_spectra__*.pdf（与标记无关）：两分法报告的每条可信新肽段
  （credible_novel = TRUE），上面板＝新肽段（corrected_*），下面板＝打分最
  相近的已注释竞争者。每页含候选 1 谱图＋ppm 面板、胜选谱图＋ppm 面板、信息
  行（类型/谱图号/打分/修饰变化/判定与结论）、离子梯、Base peak 与 Mass
  error 标注；蛋白与侧翼残基取自主表（cand1_*/corrected_*）。
- 无标记侧：code/spectrum_pdf_export_A.R 写阶段 A 的
  error_match_spectra__*.pdf（候选误鉴定），同一门控，另有
  novel_peptide_spectra__*.pdf：可信新肽段与标记无关——无标记两分法报告从
  阶段 A 逐谱候选与第一步无标记谱计算同一套六重存疑规则，PDF 逐条渲染
  （上＝新肽段，下＝最强已注释竞争者）。
- 存储拆分：三类 PDF 均按 raw 落盘（`<类目>__<raw>.pdf`），组内超过
  `pdf_pages_per_file`（config_step1 键，默认 1000 页）即拆 `__part2/3…`；
  逐谱索引 `spectrum_pdf_index.tsv`（File_Name → 类目/文件/页码）落在
  1.unlabeled_result/ 与 2.labeled_result/ 的 PDF 旁。
- 结果浏览：网页结果页可整夹加载从别的电脑拷来的项目文件夹（选目录→Load）：
  报告、谱图与 PDF 卡片只读该文件夹——无需重跑。Clear 切回本地当前结果。
- 单谱导出：code/single_spectrum_export.R（库版）或独立的
  single_spectrum_export_standalone.R（自包含、可移植）——在对话框里选
  pFind-Filtered.spectra 与 MGF 文件；每谱输出一个 scan<编号>.pdf（16:9）＋
  300-dpi PNG，按 MGF 词干分子文件夹存放。
