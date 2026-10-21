# isoPSV v1.1 Release Notes

isoPSV（isotope-labeled paired-spectrum verification）multi-engine competitive validation software v1.1 (2026-10-21 update).

- Multi-engine competitive searching (pFind/pUniFind family + Comet/MaxQuant/MSFragger/Sage/MetaMorpheus and more)
- Two-stage arbitration (spectrum-level competition + paired-isotope verification), eight weighted indicators with locked weights
- Web console (English by default, switchable to Chinese): parameter generation, live log, stop-point resume, result browser
- **Pure prediction mode (new)**: submit STEP 2 with no strategy row and the base identification table feeds the predictors directly; the output `pure_prediction.tsv` is the input table plus one `predict_cosine_<predictor>` column per enabled predictor (observed vs predicted theoretical-spectrum cosine)
- **Four bundled predictor integrations (new)**: pDeep2, pDeep3, AlphaPeptDeep and MS2PIP in the same run; each keeps its own Python environment (pDeep2 3.7/TF1.15, AlphaPeptDeep 3.10/torch, MS2PIP >= 3.11) and only the environment folder needs to be filled (auto-resolved to `python.exe` or `Scripts\python.exe`); one-click detection fills it automatically
- Unified strategy slots (pFind three arms and pUniFind as peer rows, slot letters alphabetical); pFind standard arm; re-score precondition guard
- Bundled demo dataset (PRIDE PXD003236 Chlamydomonas fraction) + reproduction guide; crosstalk-safe Python subprocess environment; pre-check cache

Usage: extract the zip and double-click `isoPSV.exe` (or the equivalent `start_isoPSV.bat`); the browser opens the console automatically.
See `Software_Manual.md` and `README.md` for details.

---

# isoPSV v1.1 发布说明

isoPSV（同位素标记成对谱图验证）多引擎竞争验证软件 v1.1（2026-10-21 更新）。

- 多引擎竞争搜库（pFind/pUniFind 生态＋Comet/MaxQuant/MSFragger/Sage/MetaMorpheus 等）
- 两阶段竞争仲裁（谱图级竞争＋成对同位素验证），八指标加权评分与锁定权重
- 网页控制台（默认英文，可切换中文）：参数生成、实时日志、停点续接、结果浏览
- **纯预测模式（新增）**：第二步不加任何策略行，基础鉴定表直接作为预测输入；输出 `pure_prediction.tsv` ＝ 输入表原样 ＋ 每个启用预测器一列 `predict_cosine_<预测器>`（观测谱对理论谱余弦）
- **四预测器同轮集成（新增）**：pDeep2、pDeep3、AlphaPeptDeep、MS2PIP 可同轮运行；各用独立 Python 环境（pDeep2 3.7/TF1.15、AlphaPeptDeep 3.10/torch、MS2PIP ≥ 3.11），只需填环境文件夹（自动解析到 `python.exe` 或 `Scripts\python.exe`），一键检测可自动回填
- 策略槽位行统一（pFind 三臂与 pUniFind 同构行、槽位字母字典序）；pFind standard 新臂；re-score 前置守卫
- 内置演示数据集（PRIDE PXD003236 衣藻馏分）＋复现指南；python 子进程环境净化；预检缓存

使用方法：解压后双击 `isoPSV.exe`（或 `start_isoPSV.bat`，两者等价），浏览器自动打开控制台。
详见 `Software_Manual.md` 与 `README.md`。
