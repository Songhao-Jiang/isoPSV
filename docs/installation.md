# Installation & Operation / 安装与操作

Everything below is optional except the isoPSV package itself: isoPSV runs
standalone (import result tables, skip every engine and predictor), and each
external engine or predictor is enabled only when you fill its path in STEP 2
or STEP 3. All official download links are listed so that nothing has to be
guessed.

除 isoPSV 主包本身外，下列软件全部可选：isoPSV 可独立运行（直接导入结果表、
不装任何引擎与预测器），装哪个就在第二步/第三步填哪个的路径。所有官方下载
链接一并列出，无需自行查找。

---

# English

## 1. At a glance

| Component | Role | Runtime / isolation | Official link |
|---|---|---|---|
| **isoPSV package** | web console + pipeline | bundled portable R 4.5.1 + Python 3.12.10 (no install) | [Releases](https://github.com/Songhao-Jiang/isoPSV/releases) |
| pFind 3 | database search (a/b/c/q slots) | native Windows GUI/CLI | https://pfind.ict.ac.cn/se/pfind/ |
| pUniFind | re-score / de novo (d/e slots) | its own Python + torch (CUDA optional) | https://github.com/pFindStudio/pUniFind |
| pNovo 3 | de novo, GUI export only | native Windows GUI | https://pfind.ict.ac.cn/se/pnovo/ |
| Comet | external search engine | native exe | https://comet-ms.sourceforge.net/ |
| MSFragger | external search engine (Java) | Java 11+ + msconvert | https://msfragger.nesvilab.org/ |
| Sage | external search engine (Rust) | native exe + msconvert | https://github.com/lazear/sage |
| MaxQuant | external search engine | .NET, needs `.raw` files | https://maxquant.org/ |
| MetaMorpheus | external search engine | .NET, native exe | https://github.com/smith-chem-wisc/MetaMorpheus |
| X!Tandem | external search engine | native exe | https://www.thegpm.org/tandem/ |
| PEAKS Studio | GUI-export engine (Manual only) | Windows GUI | https://www.bioinfor.com/peaks-studio/ |
| Mascot | GUI-export engine (Manual only) | local/remote server | https://www.matrixscience.com/ |
| MSConvert (ProteoWizard) | MGF->mzML for Sage / MSFragger | native exe | https://proteowizard.sourceforge.io/download.html |
| **pDeep2** | predictor | Python 3.7 + TensorFlow 1.15 | https://github.com/pFindStudio/pDeep |
| **pDeep3** | predictor | Python 3.7 + TF 1.15 (reuses pDeep2 env) | https://github.com/pFindStudio/pDeep3 |
| **AlphaPeptDeep** | predictor | Python 3.10 + torch | https://github.com/MannLabs/alphapeptdeep |
| **MS2PIP** | predictor | Python **>= 3.11** + xgboost | https://github.com/compomics/ms2pip |

Only the last four feed the theoretical-spectrum cosine columns; the engines
feed candidate competition. A run without any engine (pure prediction mode)
needs no engine at all.

## 2. isoPSV itself

1. Download `isoPSV_v1.1.zip` from the Releases page and extract it to any
   folder (an ASCII path is recommended; avoid spaces if possible).
2. Double-click `isoPSV.exe` (or run `start_isoPSV.bat` - both are equivalent).
   The console opens automatically at http://127.0.0.1:8780.
3. No administrator rights, no internet connection and no separate R/Python
   installation are needed: the portable runtimes ship inside `bin/runtime/`.

Optional verification: the console footer shows the version; `GET
/api/version` returns it as JSON.

## 3. Search engines

All engines are separate software. isoPSV never installs or downloads them for
you; it reads their result files. Two modes are supported per row in STEP 2:

- **Manual result** - you run the engine yourself and fill the result file
  (works even if the engine is installed on another machine).
- **Auto search** - you fill the engine's executable/installation path and
  isoPSV launches it on the MGF produced by STEP 1.

### pFind 3 (a, b, c, q slots)
- Download from the official pFind Studio site: https://pfind.ict.ac.cn/se/pfind/
  (site mirror: http://pfind.net/)
- Install by extracting the archive; no separate runtime needed.
- STEP 2: fill the pFind **installation folder** (the folder containing
  `bin/Searcher.exe`); a version folder such as `pFind3.2` is accepted and
  `bin` is appended automatically. The three arms are rows:
  - `b` pFind open (open search)
  - `c` pFind restrict (restricted, relaxed modification/cleavage)
  - `q` pFind standard (conventional strict parameters)
- The pFind parameter templates live in `code/pFind_templates/` and can be
  edited; isoPSV fills the database and spectra paths at run time.

### pUniFind (d, e slots)
- Source and models: https://github.com/pFindStudio/pUniFind
- Requires its own Python environment with PyTorch. Install PyTorch first
  (CUDA build if you have a GPU), then the package requirements.
- STEP 2: fill `punifind_exe` (the folder containing the `official_*.exe`
  binaries, `ckpts/`, `consts/`) and `punifind_python` (the environment
  folder; resolved to `python.exe` or `Scripts\python.exe`).
- Re-score (`d`) requires an existing open-pFind result: isoPSV blocks the
  re-score row when the open arm is absent.
- GPU supply knobs: `punifind_batch_size`, `punifind_chunk_parallel`,
  `punifind_chunk_cleanup` (in config_step2.txt).

### Comet, X!Tandem, MetaMorpheus, Sage, MSFragger, MaxQuant
| Engine | Install | STEP 2 field | Notes |
|---|---|---|---|
| Comet | https://comet-ms.sourceforge.net/ (binary) | `comet_exe` | uses STEP-1 MGF directly |
| X!Tandem | https://www.thegpm.org/tandem/ | `xtandem_exe` | uses STEP-1 MGF |
| MetaMorpheus | https://github.com/smith-chem-wisc/MetaMorpheus | `metamorpheus_exe` | uses STEP-1 MGF |
| Sage | https://github.com/lazear/sage (release binary) | `sage_exe` + `msconvert_exe` | needs MGF->mzML; Sage cannot read MGF directly |
| MSFragger | https://msfragger.nesvilab.org/ (`MSFragger.jar`, Java 11+) | `msfragger_exe` + `msconvert_exe` | JVM heap via `msfragger_xmx` |
| MaxQuant | https://maxquant.org/ (`MaxQuantCmd.exe`) | `maxquant_exe` + `search_file_raw` | MaxQuant reads `.raw`, not MGF |

- MSConvert (ProteoWizard): https://proteowizard.sourceforge.io/download.html -
  fill `msconvert_exe` once; isoPSV converts MGF to mzML for Sage/MSFragger.
- PEAKS Studio (https://www.bioinfor.com/peaks-studio/) and Mascot
  (https://www.matrixscience.com/) are GUI-only: export the result tables
  yourself and use **Manual result** rows (engine score formats: see
  `docs/parameters.md`).

## 4. Theoretical-spectrum predictors

Each predictor writes a `predict.txt` in the same shape, and isoPSV adds one
`predict_cosine_<predictor>` column per successful predictor. In STEP 3 fill:

- pDeep2 / pDeep3: the **installation folder**
- AlphaPeptDeep / MS2PIP: the **environment folder** (the folder alone is
  enough - see section 5)

### 4.1 pDeep2 (Python 3.7 + TensorFlow 1.15)
```bat
:: environment (TensorFlow 1.15 requires Python <= 3.7)
conda create -n pdeep python=3.7 -y && conda activate pdeep
pip install tensorflow==1.15.0 scipy pandas matplotlib

:: software
git clone https://github.com/pFindStudio/pDeep.git D:\Bioinformatics\pDeep2
:: verify
python -c "import tensorflow; print(tensorflow.__version__)"
```
STEP 3: fill `D:\Bioinformatics\pDeep2`. isoPSV resolves
`<root>\py37\python.exe` (or `python37`) and runs
`<root>\pDeep\pDeep2\predict.py`.

### 4.2 pDeep3 (reuses the pDeep2 environment)
```bat
git clone https://github.com/pFindStudio/pDeep3.git D:\Bioinformatics\pDeep3
```
STEP 3: fill `D:\Bioinformatics\pDeep3`. isoPSV sets `PDEEP3_ROOT` to that
folder and resolves the interpreter through the same heal chain, falling back
to the sibling `pDeep2\py37\python.exe`. pDeep3 supports QE/Lumos; other
instruments map to Lumos.

### 4.3 AlphaPeptDeep (Python 3.10 + torch)
```bat
conda create -n peptdeep python=3.10 -y && conda activate peptdeep
pip install torch:: CPU build is enough; use the CUDA index URL for GPU
pip install "transformers==4.39.3" "alphabase==1.1.0"
pip install peptdeep:: or: pip install -e <alphapeptdeep source folder>
:: verify
python -c "from peptdeep.pretrained_models import ModelManager; ModelManager().load_installed_models(); print('ok')"
```
STEP 3: fill the environment folder (e.g. `D:\Bioinformatics\Miniconda` or a
venv root). Models download on first use (network required once).

### 4.4 MS2PIP (Python >= 3.11)
```bat
:: MS2PIP 4.x requires Python >= 3.11 - a separate environment is mandatory
py -3.11 -m venv D:\Bioinformatics\ms2pip\venv
D:\Bioinformatics\ms2pip\venv\Scripts\python.exe -m pip install ms2pip
:: verify
D:\Bioinformatics\ms2pip\venv\Scripts\python.exe -m ms2pip --version
```
STEP 3: fill `D:\Bioinformatics\ms2pip\venv` (environment folder). The HCD
model (~1 GB: `model_20210416_HCD2021_B.xgboost` + `..._Y.xgboost`) downloads
on first use into `~/.ms2pip`; if the download is slow, fetch the two files
from https://genesis.ugent.be/uvpublicdata/ms2pip/ (mirror used by MS2PIP;
Zenodo record 13270668 also serves them) and place them in `%USERPROFILE%\.ms2pip`.
Model name stays default `HCD` unless you pass another name in the params field.

## 5. Environment isolation and Python resolution

**Why isolation is mandatory.** The four predictors cannot share one
interpreter: pDeep2/pDeep3 need Python 3.7 + TensorFlow 1.15, AlphaPeptDeep
needs 3.10 + torch, MS2PIP needs >= 3.11. Installing them into one environment
breaks at least two of them. Keep one environment per predictor (conda env or
venv); isoPSV launches each predictor with its own interpreter.

**What you fill.** Only the folder. isoPSV resolves it in this order:

1. `<folder>\python.exe` (conda layout, Miniconda/Anaconda root)
2. `<folder>\Scripts\python.exe` (Windows venv layout)
3. `<folder>\python37\python.exe`, `<folder>\py37\python.exe` (portable addons)
4. `dirname(<folder>)\python37\python.exe`, `dirname(<folder>)\pDeep2\py37\python.exe`
   (the pDeep3 -> pDeep2 sibling layout)

A full path to `python.exe` is also accepted and used as-is. Every resolved
path is printed in the run log as `[predictor] python path healed: ...`.

**Automatic detection.** The console's one-click detection scans common install
roots, finds `site-packages\peptdeep` and `site-packages\ms2pip` (and pDeep2 /
pDeep3 / pFind) and fills the empty predictor paths for you; user-entered
values are never overwritten.

**Subprocess hygiene.** isoPSV strips `PYTHONHOME` / `PYTHONPATH` before
launching any Python subprocess, so a machine-wide Python configuration cannot
hijack a predictor environment. It also clears `LC_CTYPE` / `LC_ALL` / `LANG`
for R subprocesses so non-ASCII paths behave.

**GPU.** pUniFind and AlphaPeptDeep optionally use CUDA; everything else is
CPU. Without a CUDA build, install the CPU wheels (as in section 4.3).

## 6. Operation guide

1. **STEP 1** - result table (`psm_file`, any tab-separated table with a
   `File_Name` column), MGF folder, sampling mode (labeled / unlabeled),
   `result_dir`. Click *Generate config*.
2. **STEP 2** - add one row per search strategy (Auto = engine path, Manual =
   result file). Default seeding is pFind open + restrict.
   **Pure prediction mode**: add **no row at all** and click *Generate* - the
   base table goes straight to the predictors.
3. **STEP 3** - predictor rows (pDeep2/pDeep3 install folders;
   AlphaPeptDeep/MS2PIP environment folders) and the arbitration weights.
4. **RUN** - live log, stop points (missing inputs are listed with the exact
   file to fill), resume by re-running: completed stages are reused via
   fingerprints.

Outputs: `<result_dir>/result_*/` with `1.unlabeled_result/`,
`2.labeled_result/`, `merge/predictors/<type>/predict.txt`,
`merge/predictor_run_manifest.txt` (one line per predictor: OK/FAIL),
and - in pure prediction mode - `pure_prediction.tsv`.

A worked example (data -> config -> run -> expected outputs) is in
[Example/README.md](../Example/README.md); the full parameter reference is in
[parameters.md](parameters.md).

## 7. Quick troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `MS2PIP` column empty, manifest FAIL | interpreter inside a venv not found (older builds) or models missing | fill the venv **root** (v1.1 resolves `Scripts\python.exe`); check `~/.ms2pip` |
| pDeep2 `py37 not found` | wrong folder level | fill the pDeep2 root that contains `py37\` and `pDeep\` |
| pDeep3 import error | TF version mismatch | use the pDeep2 `py37` environment (TF 1.15) |
| Sage / MSFragger fail | no `msconvert_exe` | install ProteoWizard and fill `msconvert_exe` |
| MaxQuant row skipped | MGF given instead of `.raw` | fill `search_file_raw` with the raw files |
| pUniFind re-score blocked | open-pFind arm missing | enable the `b` row first |
| Chinese path errors | external tools read GBK | extract isoPSV and put data under an ASCII path |

---

# 中文

## 1. 一览表

| 组件 | 作用 | 运行时／隔离 | 官方链接 |
|---|---|---|---|
| **isoPSV 主包** | 网页控制台＋流程 | 自带 R 4.5.1＋Python 3.12.10（免安装） | [Releases](https://github.com/Songhao-Jiang/isoPSV/releases) |
| pFind 3 | 数据库搜库（a/b/c/q 槽） | Windows 原生 | https://pfind.ict.ac.cn/se/pfind/ |
| pUniFind | re-score／从头测序（d/e 槽） | 独立 Python＋torch（可选 CUDA） | https://github.com/pFindStudio/pUniFind |
| pNovo 3 | 从头测序，仅 GUI 导出 | Windows 原生 GUI | https://pfind.ict.ac.cn/se/pnovo/ |
| Comet | 外部搜库引擎 | 原生 exe | https://comet-ms.sourceforge.net/ |
| MSFragger | 外部搜库引擎（Java） | Java 11+ ＋ msconvert | https://msfragger.nesvilab.org/ |
| Sage | 外部搜库引擎（Rust） | 原生 exe ＋ msconvert | https://github.com/lazear/sage |
| MaxQuant | 外部搜库引擎 | .NET，需 `.raw` | https://maxquant.org/ |
| MetaMorpheus | 外部搜库引擎 | .NET 原生 exe | https://github.com/smith-chem-wisc/MetaMorpheus |
| X!Tandem | 外部搜库引擎 | 原生 exe | https://www.thegpm.org/tandem/ |
| PEAKS Studio | GUI 导出引擎（仅 Manual） | Windows GUI | https://www.bioinfor.com/peaks-studio/ |
| Mascot | GUI 导出引擎（仅 Manual） | 本地／远程服务器 | https://www.matrixscience.com/ |
| MSConvert（ProteoWizard） | Sage／MSFragger 的 MGF→mzML | 原生 exe | https://proteowizard.sourceforge.io/download.html |
| **pDeep2** | 理论谱预测 | Python 3.7＋TensorFlow 1.15 | https://github.com/pFindStudio/pDeep |
| **pDeep3** | 理论谱预测 | Python 3.7＋TF 1.15（复用 pDeep2 环境） | https://github.com/pFindStudio/pDeep3 |
| **AlphaPeptDeep** | 理论谱预测 | Python 3.10＋torch | https://github.com/MannLabs/alphapeptdeep |
| **MS2PIP** | 理论谱预测 | Python **≥3.11**＋xgboost | https://github.com/compomics/ms2pip |

只有后四个预测器产生理论谱余弦列；搜库引擎用于候选竞争。纯预测模式不需要
任何搜库引擎。

## 2. isoPSV 主软件

1. 从 Releases 下载 `isoPSV_v1.1.zip`，解压到任意文件夹（建议 ASCII 路径，
   尽量不含空格）。
2. 双击 `isoPSV.exe`（或 `start_isoPSV.bat`，两者等价），控制台自动打开
   http://127.0.0.1:8780。
3. 无需管理员权限、无需联网、无需另装 R/Python——便携运行时都在
   `bin/runtime/` 内。

自检：控制台页脚显示版本；`GET /api/version` 返回 JSON 版本信息。

## 3. 搜库引擎

全部为独立软件，isoPSV 不代装、不代下，只读其结果。第二步每行两种方式：

- **Manual result**：自行运行引擎后填结果文件（引擎装在别的机器上也可）。
- **Auto search**：填引擎可执行文件／安装夹，isoPSV 用第一步切分出的 MGF
  自动调用。

### pFind 3（a、b、c、q 槽）
- 官方下载：https://pfind.ict.ac.cn/se/pfind/（镜像 http://pfind.net/）
- 解压即用，无需额外运行时。
- 第二步填 pFind **安装文件夹**（含 `bin/Searcher.exe`）；填 `pFind3.2`
  这类版本目录也行，软件自动补 `bin`。三臂各占一行：
  `b` 开放式、`c` 限定式（放宽修饰与酶切）、`q` 常规限定式。
- 参数模板在 `code/pFind_templates/` 可改；数据库与谱图路径运行时自动填。

### pUniFind（d、e 槽）
- 源码与模型：https://github.com/pFindStudio/pUniFind
- 需独立 Python 环境＋PyTorch（有 GPU 就装 CUDA 版）。
- 第二步填 `punifind_exe`（含 `official_*.exe`、`ckpts/`、`consts/`）与
  `punifind_python`（环境文件夹，自动解析到 `python.exe` 或
  `Scripts\python.exe`）。
- `d` 槽（re-score）依赖开放式 pFind 结果：缺 `b` 行时该行被前置守卫挡下。
- GPU 供给旋钮：`punifind_batch_size`、`punifind_chunk_parallel`、
  `punifind_chunk_cleanup`（config_step2.txt）。

### Comet、X!Tandem、MetaMorpheus、Sage、MSFragger、MaxQuant
| 引擎 | 安装 | 第二步字段 | 备注 |
|---|---|---|---|
| Comet | https://comet-ms.sourceforge.net/ | `comet_exe` | 直接用第一步 MGF |
| X!Tandem | https://www.thegpm.org/tandem/ | `xtandem_exe` | 直接用 MGF |
| MetaMorpheus | https://github.com/smith-chem-wisc/MetaMorpheus | `metamorpheus_exe` | 直接用 MGF |
| Sage | https://github.com/lazear/sage | `sage_exe`＋`msconvert_exe` | 需 MGF→mzML |
| MSFragger | https://msfragger.nesvilab.org/（`MSFragger.jar`，Java 11+） | `msfragger_exe`＋`msconvert_exe` | JVM 堆用 `msfragger_xmx` |
| MaxQuant | https://maxquant.org/（`MaxQuantCmd.exe`） | `maxquant_exe`＋`search_file_raw` | 只吃 `.raw`，不吃 MGF |

- MSConvert（ProteoWizard）：https://proteowizard.sourceforge.io/download.html
  ——填一次 `msconvert_exe`，Sage／MSFragger 的 MGF→mzML 由软件自动完成。
- PEAKS Studio（https://www.bioinfor.com/peaks-studio/）与 Mascot
  （https://www.matrixscience.com/）只能 GUI 导出，用 **Manual result** 行
  导入（分数格式见 `docs/parameters.md`）。

## 4. 理论谱图预测器

四个预测器产物同构（`predict.txt`），每个成功的预测器在输出里加一列
`predict_cosine_<预测器>`。第三步填法：pDeep2／pDeep3 填**安装文件夹**；
AlphaPeptDeep／MS2PIP 填**环境文件夹**（只填文件夹即可，见第 5 节）。

### 4.1 pDeep2（Python 3.7＋TensorFlow 1.15）
```bat
conda create -n pdeep python=3.7 -y && conda activate pdeep
pip install tensorflow==1.15.0 scipy pandas matplotlib
git clone https://github.com/pFindStudio/pDeep.git D:\Bioinformatics\pDeep2
python -c "import tensorflow; print(tensorflow.__version__)"
```
第三步填 `D:\Bioinformatics\pDeep2`；软件解析 `<root>\py37\python.exe`
（或 `python37`）并调用 `<root>\pDeep\pDeep2\predict.py`。

### 4.2 pDeep3（复用 pDeep2 环境）
```bat
git clone https://github.com/pFindStudio/pDeep3.git D:\Bioinformatics\pDeep3
```
第三步填 `D:\Bioinformatics\pDeep3`；软件把 `PDEEP3_ROOT` 指过去，解释器经
同一自愈链解析，必要时回落到同级 `pDeep2\py37\python.exe`。pDeep3 支持
QE/Lumos，其他仪器映射为 Lumos。

### 4.3 AlphaPeptDeep（Python 3.10＋torch）
```bat
conda create -n peptdeep python=3.10 -y && conda activate peptdeep
pip install torch
pip install "transformers==4.39.3" "alphabase==1.1.0"
pip install peptdeep
python -c "from peptdeep.pretrained_models import ModelManager; ModelManager().load_installed_models(); print('ok')"
```
第三步填环境文件夹（如 `D:\Bioinformatics\Miniconda` 或 venv 根）；模型首次
使用需联网下载一次。

### 4.4 MS2PIP（Python ≥3.11）
```bat
py -3.11 -m venv D:\Bioinformatics\ms2pip\venv
D:\Bioinformatics\ms2pip\venv\Scripts\python.exe -m pip install ms2pip
D:\Bioinformatics\ms2pip\venv\Scripts\python.exe -m ms2pip --version
```
第三步填 `D:\Bioinformatics\ms2pip\venv`。HCD 模型（约 1 GB）首次使用时下载到
`~/.ms2pip`；下载慢可从 https://genesis.ugent.be/uvpublicdata/ms2pip/
（MS2PIP 自带镜像；Zenodo record 13270668 同文件）取这两个 xgboost 文件放入
`%USERPROFILE%\.ms2pip`。模型名默认 `HCD`。

## 5. 环境隔离与 Python 解析

**为什么要隔离**：四个预测器无法共用一个解释器——pDeep2/pDeep3 要
Python 3.7＋TF 1.15，AlphaPeptDeep 要 3.10＋torch，MS2PIP 要 ≥3.11。混装至少
有两个会坏。每个预测器一个环境（conda env 或 venv），isoPSV 各用各的。

**填什么**：只填文件夹。解析顺序：

1. `<文件夹>\python.exe`（conda 布局）
2. `<文件夹>\Scripts\python.exe`（Windows venv 布局）
3. `<文件夹>\python37\python.exe`、`<文件夹>\py37\python.exe`（便携附件包）
4. `上级\python37\python.exe`、`上级\pDeep2\py37\python.exe`（pDeep3→pDeep2 同级布局）

直接填 `python.exe` 全路径也接受，原样使用。每次解析都会在运行日志打印
`[预测器] python 路径自愈：… -> …`。

**自动检测**：网页“一键检测”扫描常见安装根，按 `site-packages\peptdeep`、
`site-packages\ms2pip`（以及 pFind／pDeep2／pDeep3）自动回填空的预测器路径，
不覆盖手填值。

**子进程净化**：启动任何 Python 子进程前剥离 `PYTHONHOME`／`PYTHONPATH`，
机器级配置不会劫持预测器环境；R 子进程清 `LC_CTYPE`／`LC_ALL`／`LANG`，
避免非 ASCII 路径异常。

**GPU**：pUniFind 与 AlphaPeptDeep 可选 CUDA，其余 CPU；无 CUDA 时装 CPU 轮子
（见 4.3）。

## 6. 操作流程

1. **第一步**：结果表（`psm_file`，任何带 `File_Name` 列的制表符表）、MGF
   文件夹、采样模式（标记／无标记）、`result_dir`，点“生成参数文件”。
2. **第二步**：每个策略一行（Auto＝填引擎路径，Manual＝填结果文件），默认
   预置 pFind open＋restrict。**纯预测模式**：一行都不加直接生成——基础表直送
   预测器。
3. **第三步**：预测器行（pDeep2/pDeep3 填安装夹；AlphaPeptDeep/MS2PIP 填环境
   夹）与仲裁权重。
4. **运行**：实时日志、停点提示（缺什么会列出要填的文件）、重跑续接（已完成
   阶段按指纹复用）。

输出：`<result_dir>/result_*/` 下 `1.unlabeled_result/`、
`2.labeled_result/`、`merge/predictors/<类型>/predict.txt`、
`merge/predictor_run_manifest.txt`（每个预测器一行 OK/FAIL）；纯预测模式另有
`pure_prediction.tsv`。

完整示例见 [Example/README.md](../Example/README.md)，全部参数见
[parameters.md](parameters.md)。

## 7. 常见问题

| 现象 | 原因 | 处理 |
|---|---|---|
| MS2PIP 列为空、manifest FAIL | venv 解释器未解析到（早期实现）或缺模型 | 填 venv 根（v1.1 会解析 `Scripts\python.exe`）；检查 `~/.ms2pip` |
| pDeep2 报 `py37 not found` | 文件夹层级不对 | 填含 `py37\` 与 `pDeep\` 的 pDeep2 根目录 |
| pDeep3 导入报错 | TF 版本不符 | 用 pDeep2 的 `py37`（TF 1.15） |
| Sage／MSFragger 失败 | 缺 `msconvert_exe` | 装 ProteoWizard 后填 `msconvert_exe` |
| MaxQuant 行被跳过 | 给的是 MGF 不是 `.raw` | 用 `search_file_raw` 填 raw 文件 |
| pUniFind re-score 被挡 | 缺开放式 pFind 行 | 先启用 `b` 行 |
| 中文路径报错 | 外部工具按 GBK 读路径 | isoPSV 与数据放 ASCII 路径 |
