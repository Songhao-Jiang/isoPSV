# Parameters Reference / 参数参考

Every configuration key of isoPSV, its allowed values, default and effect.
Parameter files are tab-separated `key<TAB>value` lines with `#` comments;
they live in the result directory (`result_dir`).

---

# 中文版

isoPSV 全部配置键、取值、默认值与作用。参数文件为制表符 `键<TAB>值` 行＋
`#` 注释；位于结果目录（`result_dir`）。

## config_step1.txt（第一步）

| 键 | 取值 | 默认 | 作用 |
|---|---|---|---|
| result_dir | 绝对路径 | 必填 | 项目总目录；一切产出落此处。建议本地硬盘，勿用 U 盘 |
| search_engine | pfind / maxquant / msfragger / sage / comet | pfind | 第一步原始结果 psm_file 来自哪个引擎——决定读取口径与非 pFind 引擎的自动适配（maxquant 走两参适配器自动转换；sage/comet/msfragger 结果建议经 STEP2 槽位） |
| psm_file | 绝对路径 | 必填 | 原始搜库结果表；任何带 File_Name 列的制表符表可直接用；X!Tandem XML/Mascot/PEAKS CSV 自动识别转换 |
| fasta | 绝对路径 | 必填 | 参考蛋白质组数据库（FASTA） |
| label_type | 多行键：N15 / C13 / 自定义行 | unlabeled | 标记类型。自定义行语法 `名称\|位置\|氨基酸\|组成`，位置 ANY/PEP_N/PEP_C/PRO_N/PRO_C；如 `NeuCode8\|ANY\|K\|C(-6)13C(6)N(-2)15N(2)` |
| unlabeled_channel | 1 / 0 | 1 | 1＝样本含无标记通道（按 Label 列提取）；0＝全标记样本（提取全部鉴定谱，标记按固定修饰搜库） |
| enzyme | trypsin / lysarginase / lysc / lysn / argc / gluc / aspn / chymotrypsin / none | trypsin | 酶切类型；写入各引擎模板并参与分类规则。Trypsin/LysargiNase 均 P 位不切 |
| mgf | 多行/文件夹/raw | 必填 | MGF 列表：每行一个文件、或文件夹（展开 .mgf）、或 **raw 文件/文件夹（第一步自动经 pParse 转 MGF）** |
| all_mgf | 文件夹 | — | 标记样本必填：原始 MGF 文件夹（含标记谱） |
| auto_search | 1 / 0 | 1 | pFind 自动搜库开关；任何外部引擎选自动模式时自动置 1 |
| pfind_exe | 安装文件夹 | 检测 | pFind 安装目录（含 bin\Searcher.exe 那一级） |
| use_pdeep2 | 1 / 2 / 0 | 1 | pDeep2：1＝自动运行；2＝手动（predict_file）；0＝关闭 |
| pdeep2_exe | 安装文件夹 | 检测 | pDeep2 安装目录 |
| instrument | Lumos / QE / Velos / Elite | Lumos | 质谱仪类型；预测器 params 留空时取此键 |
| pdf_on | 1 / 0 | 0 | 逐谱 PDF 导出开关（复核材料） |
| denovo_lookback | keep / all | keep | 从头测序回帖模式（保留原始/全部候选） |
| candidate_scope | filtered / top5 | filtered | 候选摄入口径：filtered＝依赖引擎 1% FDR；top5＝忽略 FDR 每谱前 5 候选进竞争池（小样本备用） |
| proteogenomics_mode | 0 / 1 / 2 | 0 | 0＝全量；1＝纯蛋白质基因组学；2＝已注释矫正 |
| mem_guard_mb | 数值 | 自适应 | 内存守卫上限 MB |
| max_cores | 数值 | 自动评估 | 手动覆盖读取并行度 |

## config_step2.txt（第二步：搜库策略）

| 键 | 取值 | 作用 |
|---|---|---|
| search_file | 多行：文件/文件夹/`(auto:引擎名)` | 各策略结果来源；auto 标记由自动搜库段处理 |
| search_file_type | 与 search_file 平行 | 策略槽位类型（maxquant/comet/msfragger/sage/…）；缺省按命名推断 |
| search_file_raw | raw 列表（; 连接） | MaxQuant 自动搜库的 raw 输入（MaxQuant 不能读 MGF） |
| `<engine>_exe` | 路径 | 外部引擎可执行文件；多版本共存时检测自动选版本号最新 |
| `<engine>_mode` | 1 | 自动运行标记 |
| msconvert_exe | 路径 | ProteoWizard msconvert（Sage/MSFragger 的 MGF→mzML） |
| msfragger_xmx | 如 8G | MSFragger JVM 堆上限 |
| punifind_batch_size / _chunk_parallel / _chunk_cleanup | 64 / 1 / 0 | pUniFind GPU 供给三旋钮 |

## config_step3.txt（第三步：仲裁）

| 键 | 取值 | 默认 | 作用 |
|---|---|---|---|
| weights_A | 四数逗号分隔 | 0.1306,0.4418,0.2381,0.1895 | 阶段 A 权重（ul_by,ul_all,ul_cov,predict） |
| threshold_A | 数值 | 0.0204 | 阶段 A 阈值（ACS 差值门） |
| weights_B | 6 或 8 数逗号分隔 | — | 阶段 B 分标记权重；缺省用统一锁定值。W8 模板一键写入 |
| threshold_B | 数值 | — | 阶段 B 阈值；W8 模板写入 |
| predict_file | 路径 | 自动 | pDeep2 predict.txt（手动模式填） |
| use_pdeep2 | 0/1/2 | 继承 | 同第一步 |

## config_predictor.txt（预测器）

| 键 | 作用 |
|---|---|
| predictor_primary | 主预测器（stage A 的 predict 列来源）；none＝关闭 |
| predictorN_type | pdeep2 / pdeep3 / ms2pip / alphapeptdeep |
| predictorN_python | 安装文件夹（pdeep2/pdeep3）或 Python 环境文件夹（ms2pip/alphapeptdeep） |
| predictorN_params | 留空＝默认（instrument 键 / 内置 HCD 模型） |
| predictorN_mode | 1＝启用 0＝停用 |

---

# English

## config_step1.txt (step 1)

| Key | Values | Default | Effect |
|---|---|---|---|
| result_dir | absolute path | required | Project root; all outputs go here. Local disk recommended |
| search_engine | pfind / maxquant / msfragger / sage / comet | pfind | Which engine produced psm_file - sets the reading convention; maxquant auto-converts via its two-arg adapter; sage/comet/msfragger results are best fed through STEP2 slots |
| psm_file | absolute path | required | Original search result; any tab-separated table with a File_Name column works directly; X!Tandem XML / Mascot / PEAKS CSV auto-detected and converted |
| fasta | absolute path | required | Reference proteome database (FASTA) |
| label_type | multi-row: N15 / C13 / custom rows | unlabeled | Label types. Custom grammar `name\|position\|aa\|composition`, positions ANY/PEP_N/PEP_C/PRO_N/PRO_C, e.g. `NeuCode8\|ANY\|K\|C(-6)13C(6)N(-2)15N(2)` |
| unlabeled_channel | 1 / 0 | 1 | 1 = sample has an unlabeled channel (extract by Label column); 0 = fully labeled (extract all identified spectra; labels searched as fixed mods) |
| enzyme | trypsin / lysarginase / lysc / lysn / argc / gluc / aspn / chymotrypsin / none | trypsin | Enzyme; written into engine templates, drives classification rules. Trypsin/LysargiNase both blocked at P |
| mgf | multi-row / folder / raw | required | MGF list: one file per line, a folder (expanded to .mgf), or raw files/folders (auto-converted to MGF via pParse at step 1) |
| all_mgf | folder | - | Required for labeled samples: original MGF folder (labeled spectra) |
| auto_search | 1 / 0 | 1 | pFind auto-search switch; set to 1 automatically when any external engine is in Auto mode |
| pfind_exe | install folder | detected | pFind install directory (the level containing bin\Searcher.exe) |
| use_pdeep2 | 1 / 2 / 0 | 1 | pDeep2: 1 = auto-run; 2 = manual (predict_file); 0 = off |
| pdeep2_exe | install folder | detected | pDeep2 install directory |
| instrument | Lumos / QE / Velos / Elite | Lumos | Mass spectrometer type; used by predictors when params empty |
| pdf_on | 1 / 0 | 0 | Per-spectrum PDF export switch |
| denovo_lookback | keep / all | keep | De novo look-back mode (keep originals / all candidates) |
| candidate_scope | filtered / top5 | filtered | Candidate intake: filtered = rely on engine 1% FDR; top5 = ignore FDR, top-5 candidates per spectrum enter the pool (small-sample fallback) |
| proteogenomics_mode | 0 / 1 / 2 | 0 | 0 = full; 1 = pure proteogenomics; 2 = annotated-correction only |
| mem_guard_mb | number | adaptive | Memory-guard cap in MB |
| max_cores | number | auto-assessed | Manual override of read parallelism |

## config_step2.txt (step 2: search strategies)

| Key | Values | Effect |
|---|---|---|
| search_file | multi-row: file/folder/`(auto:engine)` | Per-strategy result sources; auto markers handled by the auto-search stage |
| search_file_type | parallel to search_file | Strategy slot type (maxquant/comet/msfragger/sage/...); inferred from naming when omitted |
| search_file_raw | raw list (";-joined) | Raw inputs for the MaxQuant auto search (MaxQuant cannot read MGF) |
| `<engine>_exe` | path | External engine executable; with multiple versions installed, detection picks the newest by parsed version tuples |
| `<engine>_mode` | 1 | Auto-run marker |
| msconvert_exe | path | ProteoWizard msconvert (MGF->mzML for Sage/MSFragger) |
| msfragger_xmx | e.g. 8G | MSFragger JVM heap cap |
| punifind_batch_size / _chunk_parallel / _chunk_cleanup | 64 / 1 / 0 | pUniFind GPU supply knobs |

## config_step3.txt (step 3: arbitration)

| Key | Values | Default | Effect |
|---|---|---|---|
| weights_A | four comma-separated numbers | 0.1306,0.4418,0.2381,0.1895 | Stage-A weights (ul_by, ul_all, ul_cov, predict) |
| threshold_A | number | 0.0204 | Stage-A threshold (ACS gap gate) |
| weights_B | 6 or 8 comma-separated numbers | - | Stage-B per-label weights; unified locked defaults when absent; one-click via the W8 templates |
| threshold_B | number | - | Stage-B threshold; written by the W8 templates |
| predict_file | path | auto | pDeep2 predict.txt (fill for manual mode) |
| use_pdeep2 | 0/1/2 | inherited | Same as step 1 |

## config_predictor.txt (predictors)

| Key | Effect |
|---|---|
| predictor_primary | Primary predictor (source of the stage-A predict column); none = off |
| predictorN_type | pdeep2 / pdeep3 / ms2pip / alphapeptdeep |
| predictorN_python | Installation folder (pdeep2/pdeep3) or Python env folder (ms2pip/alphapeptdeep) |
| predictorN_params | Empty = defaults (instrument key / built-in HCD model) |
| predictorN_mode | 1 = enabled, 0 = disabled |
