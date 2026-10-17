# Installation / 安装

## Installation

No installation is required. Download the two packages from the
[Releases](https://github.com/Songhao-Jiang/isoPSV/releases) page:

1. **isoPSV_v1.0.zip** (~100 MB) - the main package: the complete software
   with bundled portable runtimes (R 4.5.1 + Python 3.12.10; no admin
   rights, no internet), the local web console, and the bilingual user
   manuals (PDF).
2. **Theoretical-spectrum predictors (optional)** - pDeep2 and pDeep3 are
   separate software from pFind Studio; download them from the official
   sources and set up the environment:
   - pDeep2: https://github.com/pFindStudio/pDeep
   - pDeep3: https://github.com/pFindStudio/pDeep3 (models included)
   - Environment (both): Python 3.7 + TensorFlow 1.15, then
     `pip install tensorflow==1.15.0 scipy pandas matplotlib`
   - In STEP 3 fill each predictor's installation folder; parameters
     default automatically. MS2PIP (`pip install ms2pip`) and AlphaPeptDeep
     (`pip install peptdeep`) follow the same pattern with their own
     Python environments.

Steps: extract the main zip to any folder (an ASCII path is recommended),
double-click `isoPSV.exe` (or the equivalent `start_isoPSV.bat`), and the console opens automatically at
http://127.0.0.1:8780. Search engines (pFind, pUniFind, MaxQuant,
MSFragger, ...) are separate software (see Software_Manual.md); their
result files can also be imported without a local installation.


---

## 安装

软件免安装。从 [Releases](https://github.com/Songhao-Jiang/isoPSV/releases)
页面下载两个包：

1. **isoPSV_v1.0.zip**（约 100 MB）——主软件包：完整软件＋自带绿色运行时
   （R 4.5.1＋Python 3.12.10；无需管理员权限、无需联网）＋本地网页控制台＋
   双语使用说明书（PDF）。
2. **理论谱图预测器（可选）**——pDeep2 与 pDeep3 为 pFind Studio 的独立
   软件，请从官方来源下载并配置环境：
   - pDeep2：https://github.com/pFindStudio/pDeep
   - pDeep3：https://github.com/pFindStudio/pDeep3（含模型）
   - 运行环境（两者相同）：Python 3.7＋TensorFlow 1.15，然后
     `pip install tensorflow==1.15.0 scipy pandas matplotlib`
   - 在 STEP 3 填各预测器的安装文件夹即可，参数全部默认自动。MS2PIP
     （`pip install ms2pip`）与 AlphaPeptDeep（`pip install peptdeep`）
     同模式、各自配 Python 环境。

步骤：主包解压到任意文件夹（建议 ASCII 路径），双击 `isoPSV.exe`（或等价的 `start_isoPSV.bat`），控制台
自动打开 http://127.0.0.1:8780。搜库引擎（pFind、pUniFind、MaxQuant、
MSFragger 等）为独立软件（见 Software_Manual.md）；也可不装引擎、直接导入
其结果文件。

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
MS2PIP: `pip install ms2pip`; AlphaPeptDeep: `pip install peptdeep`.

---

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
MS2PIP：`pip install ms2pip`；AlphaPeptDeep：`pip install peptdeep`。

