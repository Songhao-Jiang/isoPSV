# GUI Guide / 网页界面指南

The web console (double-click `isoPSV.exe`, opens http://127.0.0.1:8780;
default English, 🌐 toggles Chinese).

---

# 中文版

## STEP 1 参数填写与配置文件生成

1. **自动识别软件路径**：点"One-click detect & fill"，按"扫描范围/Extra
   scan roots"（≤10 个目录，分号分隔）搜索本地安装并自动回填；多版本共存
   时自动选版本号最新（按版本数字推断）。GPU 行经注册表探测（滤虚拟显卡）。
2. **输入数据**：result_dir（结果目录）、psm_file（原始搜库结果——占位提示
   随 search_engine 联动）、fasta、mgf（MGF/raw 均可，raw 自动经 pParse 转）。
3. **标记类型**：N15/C13 单选；第二类自定义修饰行（名称｜位置｜氨基酸｜组成
   四格＋常用预设一键填充）。
4. **其他补充选项**：pdf_on、denovo_lookback、instrument、candidate_scope、
   proteogenomics_mode。

## STEP 2 多引擎策略

所有策略（含 pFind 三臂与 pUniFind）统一为**槽位行**：下拉按槽位字母排序
（b pFind open / c restrict / d pUniFind re-score / e de novo / f–p 外部
引擎按软件名字典序 / q pFind standard）。每行 **Auto search / Manual
result**：Auto＝填软件路径（pFind 臂填 pFind 安装目录，pUniFind 填 punifind
目录——自动模式另需下方 punifind_python 字段；仅 MaxQuant 另需 raw 输入，
其余自动用第一步 MGF）；Manual＝填结果文件。GUI 类引擎（PEAKS/Mascot/
pNovo/X!Tandem 等）仅 Manual。默认预置 open＋restrict 两行 Auto。有
pUniFind 矫正行而无 pFind open 行时显示前置提示（矫正输入＝open 结果）。

**纯预测模式**：策略行一个都不加并点生成，即进入纯预测模式——第二步不再竞争
搜库，第一步的基础鉴定表直接作为预测输入，跳过阶段 A/B 仲裁与主表导出；预测器
跑完后输出 `pure_prediction.tsv`（输入表原样＋每个预测器一列
`predict_cosine_<预测器>`，无预测留空）。

## STEP 3 理论谱图预测器

pDeep2/pDeep3 行填安装文件夹，MS2PIP/AlphaPeptDeep 行填 Python 环境文件夹
（只填环境文件夹即可，自动解析到 `python.exe` 或 `Scripts\python.exe`；MS2PIP 需 Python ≥ 3.11）；
参数全部默认自动；"⚡ 一键填写路径与参数"从检测缓存补空路径（不覆盖手填）。

## 参数汇总与生成（Parameters summary & generate）

统一"生成全部参数文件（config_step1＋step2＋predictor）/保存参数/加载参数"；
切页自动生成；参数概览表实时渲染。

## RUN

启动前自动检查并按需生成参数；运行时状态条＋时间线（各模块墙钟/CPU）＋
日志面板；停点时按提示补填后再点运行（断点续跑）。

## 结果报告

浏览阶段 A/B 报告、两分法报告与逐谱 PDF 卡片；谱图查看器同标注规则；
结果浏览模式可只读加载另一台电脑的项目文件夹。

---

# English

## STEP 1 Parameter entry

1. **Auto-detect software**: "One-click detect & fill" scans the "Extra scan
   roots" (up to 10 directories, semicolon-separated) and fills the paths;
   with multiple versions installed the newest version wins (parsed version
   tuples). GPU shown via a registry probe (virtual adapters filtered).
2. **Input data**: result_dir, psm_file (placeholder follows search_engine),
   fasta, mgf (MGF or raw - raw auto-converts via pParse).
3. **Label types**: N15/C13 radio; Class-2 custom modification rows
   (name | position | amino acids | composition) with one-click presets.
4. **Other options**: pdf_on, denovo_lookback, instrument, candidate_scope,
   proteogenomics_mode.

## STEP 2 Multi-engine strategies

Each engine row has an **Auto search / Manual result** dropdown: Auto = fill
the software path (only MaxQuant additionally needs raw input; other engines
reuse the step-1 MGF); Manual = fill the result file. GUI-export engines
(PEAKS/Mascot/pNovo/X!Tandem ...) are Manual-only.

**Pure prediction mode**: add no strategy row at all and click generate. STEP 2
then runs no competitive search; the STEP 1 base identification table feeds the
predictors directly, stage A/B arbitration and the master-table export are
skipped, and the run writes `pure_prediction.tsv` (the input table unchanged
plus one `predict_cosine_<predictor>` column per predictor; missing predictions
stay empty).

## STEP 3 Predictors

pDeep2/pDeep3 rows take installation folders; MS2PIP/AlphaPeptDeep rows take
Python environment folders (the folder alone is enough - it is resolved to
`python.exe` or `Scripts\python.exe`; MS2PIP needs Python >= 3.11); parameters
all default automatically;
"Auto-fill paths & params" fills empty paths from the detection cache
(never overwrites user-entered values).

## Parameters summary & generate

One place to generate all config files (config_step1 + step2 + predictor),
save and load form values; auto-generates on page entry; live overview table.

## RUN

Checks (and if needed generates) the configs before starting; live status,
per-module timeline (wall/CPU) and a log panel; at stop points fill the
prompted keys and click Run again (checkpoint resume).

## Results

Browse stage-A/B reports, the two-part reports and per-spectrum PDF cards;
the spectrum viewer mirrors the PDF annotations; a read-only browse mode can
load a project folder copied from another computer.
