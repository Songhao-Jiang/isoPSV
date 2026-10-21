# isoPSV Multi-Engine Competitive Validation

Multi-engine edition of the isoPSV pipeline: identification results of
multiple search engines (pFind, pUniFind, MaxQuant, MSFragger, Sage, Comet,
...) over the same spectra set are unified, decoy-removed and arbitrated per
spectrum - weighted scoring, paired isotope-label evidence and six doubt
rules produce error-spectrum correction and novel-peptide validation, with
per-spectrum PDFs. All judgment metrics are recomputed from the spectra;
engine scores are never used.

## Documentation / 文档

| Topic / 主题 | File / 文件 |
|---|---|
| **Installation, all engines & predictors (links, env isolation, operation) / 安装指南（全部搜库软件与预测器：下载链接、环境隔离、操作）** | [docs/installation.md](docs/installation.md) |
| GUI (web console) / 网页界面 | [docs/gui.md](docs/gui.md) |
| Command line / 命令行 | [docs/cli.md](docs/cli.md) |
| **All parameters (complete key reference) / 全部参数（键值速查表）** | [docs/parameters.md](docs/parameters.md) |
| Pipeline internals (strategies, arbitration, outputs) / 流程原理 | [docs/pipeline.md](docs/pipeline.md) |
| Function skills (scenario-locked tool calls) / 函数技能（场景锁定调用） | [docs/skills.md](docs/skills.md) |
| Minimal demo (data -> config -> run) / 最小示例 | [Example/README.md](Example/README.md) |
| Core R package / 核心 R 包 | [Songhao-Jiang/ReadMGF](https://github.com/Songhao-Jiang/ReadMGF) |

Quick start: extract the release zip, double-click `isoPSV.exe` (or the equivalent `start_isoPSV.bat`; bundled
R/Python runtimes; console opens at http://127.0.0.1:8780), fill STEP 1-3,
Run. Details: [docs/installation.md](docs/installation.md) /
[Example/README.md](Example/README.md).

Search engines (pFind, pUniFind, Comet, MSFragger, Sage, MaxQuant,
MetaMorpheus, X!Tandem) and predictors (pDeep2, pDeep3, AlphaPeptDeep,
MS2PIP) are separate software; every official download link, the per-predictor
Python environment isolation and the STEP-by-STEP operation guide are in
[docs/installation.md](docs/installation.md).

### Predictor setup (concrete commands)

```bat
:: 1) Python 3.7 environment (pick one)
conda create -n pdeep python=3.7 -y && conda activate pdeep
::    or install Python 3.7.9 from python.org and use its pip

:: 2) dependencies (TensorFlow 1.15 requires Python <= 3.7)
pip install tensorflow==1.15.0 scipy pandas matplotlib

:: 3) software
git clone https://github.com/pFindStudio/pDeep.git  D:\Bioinformatics\pDeep2
git clone https://github.com/pFindStudio/pDeep3.git D:\Bioinformatics\pDeep3

:: 4) verify (optional)
python -c "import tensorflow; print(tensorflow.__version__)"
```

In STEP 3 fill the installation folders (`D:\Bioinformatics\pDeep2`,
`D:\Bioinformatics\pDeep3`); parameters stay empty (defaults apply).
MS2PIP: `pip install ms2pip` (needs Python >= 3.11); AlphaPeptDeep: `pip install peptdeep`
(Python 3.10 + torch). Each predictor keeps its own environment; fill the
environment folder (auto-resolved to `python.exe` or `Scripts\python.exe`) and
the built-in one-click detection fills it for you.

### Pure prediction mode

Submit STEP 2 with **no strategy row** (all competitive slots off) to run the
pure prediction mode: the base identification table (STEP 1 `psm_file`) feeds
the predictors directly, stage A/B arbitration and the master-table export are
skipped, and the output `pure_prediction.tsv` is the input table unchanged plus
one `predict_cosine_<predictor>` column per enabled predictor (observed vs
predicted theoretical-spectrum cosine; missing predictions stay empty, never
0). Short peptides below `min_pep_len` do not enter prediction. Any predictor
fails individually without interrupting the run (see
`merge/predictor_run_manifest.txt`).


### PTM-enrichment correction

For PTM-enriched samples (phosphorylation, ubiquitin, ...): set
`unlabeled_channel=0` (no unlabeled channel) and declare the enrichment PTM
as a custom modification label (Class-2 row), parallel to the dimethyl/NeuCode
presets - e.g. phosphorylation `Phospho|ANY|STY|H(1)O(3)P(1)`, ubiquitin K-GG
remnant `GG|ANY|K|H(4)C(2)N(2)O(1)`; stage A then performs competitive
arbitration and error-match correction on the enriched spectra.

---

# isoPSV 多引擎竞争验证流程

isoPSV 多引擎版：汇集多个搜库引擎（pFind、pUniFind、MaxQuant、MSFragger、
Sage、Comet 等）对同一谱图集的鉴定结果，经交换格式统一、诱饵去除与修订后
逐谱竞争仲裁——加权评分、成对同位素标记证据与六重存疑规则输出错误谱图校正
与新肽段验证，并生成逐谱 PDF。全部判定指标由谱图重新计算，引擎打分从不参与。

## 文档

| 主题 | 文件 |
|---|---|
| 安装与预测器 | [docs/installation.md](docs/installation.md) |
| 网页界面 | [docs/gui.md](docs/gui.md) |
| 命令行 | [docs/cli.md](docs/cli.md) |
| **全部参数（键值速查表）** | [docs/parameters.md](docs/parameters.md) |
| 流程原理（策略槽位/仲裁/输出） | [docs/pipeline.md](docs/pipeline.md) |
| 函数技能（场景锁定调用） | [docs/skills.md](docs/skills.md) |
| 最小示例（数据→配置→运行） | [Example/README.md](Example/README.md) |
| 核心 R 包 | [Songhao-Jiang/ReadMGF](https://github.com/Songhao-Jiang/ReadMGF) |

快速开始：解压 Release 的 zip，双击 `isoPSV.exe`（或等价的 `start_isoPSV.bat`；自带 R/Python 运行时；
控制台自动开 http://127.0.0.1:8780），填 STEP 1–3，Run。详见
[docs/installation.md](docs/installation.md) 与 [Example/README.md](Example/README.md)。
### PTM 富集数据的鉴定校正

PTM 富集样本（磷酸化/泛素化等）：`unlabeled_channel=0`（不含无标记通道），把富集 PTM 写成第二类自定义修饰行（与二甲基/NeuCode 预设并列）——如磷酸化 `Phospho|ANY|STY|H(1)O(3)P(1)`、泛素化 K-GG 残基 `GG|ANY|K|H(4)C(2)N(2)O(1)`；阶段 A 对富集谱图做竞争仲裁与错误匹配校正。


### 预测器安装（具体命令）

```bat
:: 1) Python 3.7 环境（任选其一）
conda create -n pdeep python=3.7 -y && conda activate pdeep
::    或从 python.org 安装 3.7.9 后用其 pip

:: 2) 依赖（TensorFlow 1.15 要求 Python <= 3.7）
pip install tensorflow==1.15.0 scipy pandas matplotlib

:: 3) 软件
git clone https://github.com/pFindStudio/pDeep.git  D:\Bioinformatics\pDeep2
git clone https://github.com/pFindStudio/pDeep3.git D:\Bioinformatics\pDeep3

:: 4) 验证（可选）
python -c "import tensorflow; print(tensorflow.__version__)"
```

STEP 3 里填安装文件夹（`D:\Bioinformatics\pDeep2`、
`D:\Bioinformatics\pDeep3`）即可；参数留空＝全部默认。
MS2PIP：`pip install ms2pip`（需 Python ≥ 3.11）；AlphaPeptDeep：`pip install peptdeep`
（Python 3.10 ＋ torch）。各预测器各用独立环境，只填环境文件夹即可（自动解析到
`python.exe` 或 `Scripts\python.exe`），也可用一键检测自动回填。

### 纯预测模式

第二步**不加任何策略行**（竞争槽位全关）即进入纯预测模式：第一步的基础鉴定表
（`psm_file`）直接作为预测输入，跳过阶段 A/B 仲裁与主表导出，输出
`pure_prediction.tsv` ＝ 输入表原样 ＋ 每个启用预测器一列
`predict_cosine_<预测器>`（观测谱对理论谱余弦；无预测留空、不填 0）。低于
`min_pep_len` 的短肽不进入预测。任一预测器失败不打断整轮（见
`merge/predictor_run_manifest.txt`）。

