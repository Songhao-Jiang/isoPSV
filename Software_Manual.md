# 无标记 isoPSV 自动化流程 — 客户端代码使用说明

本目录为无标记 isoPSV 三步自动化流程的客户端代码（全部 R，R 4.3.1；**不依赖任何外部 R 包**——stringr/reshape2 的调用由 `code\compat.R` 用 base R 实现并保留原始函数名；parallel 为 base 包）。核心函数将并入 ReadMGF 包更新版。

## 开发背景

蛋白质组学质谱数据的肽段鉴定，靠搜库引擎把实验二级谱图与数据库理论谱图（或预测谱图）匹配打分，再以靶标—诱饵（Target-Decoy）竞争策略在 1% 假发现率（FDR）水平质控；鉴定结果的可信程度，基本取决于所选引擎及其打分模型。FDR 却是数据集层面的统计量：设定 1%，不等于这次实验里恰好有 1% 的鉴定必然出错，也指不出通过阈值的具体肽段中哪些 PSM 属于错误匹配。谱图质量过滤、酶切一致性检查、诱饵库 FDR 控制、重复性检验、污染与混合谱剔除这些常规质控，同样只在整体层面压错误率，回答不了单张谱图上哪一个候选才是正确结果。

单谱层面的问题不挑谱图质量：高质量谱图上，不同引擎常给出不同却同样高分的候选，单一引擎的分数裁决不了；中等质量谱图受共洗脱干扰，混合谱常见；低质量谱图证据稀疏，靠任何一个引擎都难下结论。不同引擎、不同搜库策略对同一张谱图给出矛盾鉴定的情况并不少见；到了六框翻译等扩展数据库的蛋白质基因组学，搜索空间更大、高度相似的候选更多，新肽段子集的鉴定歧义随之放大。要在这个层面得到可靠结论，需要两类互补的证据：让各引擎候选进入同一竞争池、以物理证据裁决胜负的仲裁机制，以及独立于常规匹配打分的实验约束。

isoPSV 方法针对的正是这个问题。以 *Mycolicibacterium smegmatis* 为模型的系统比较表明，15N/13C 代谢全标记把肽段及其碎片所含的 N 或 C 原子数，转化为可直接测量的非标记/标记质量差；母离子总原子数、碎片原子数、相邻碎片间的原子数增量，共同构成随序列而定的"原子组成编码"（atomic-composition encoding）。isoPSV（isotope-labeled paired-spectrum verification）框架据此用成对非标记/标记谱图，对数据库候选序列做正交物理核查。验证效果上，15N 与 13C 的可验证错误事件数和可验证比例都高于 K/R-SILAC；用于六框蛋白质基因组学时，未注释 GSSP 谱图可被已注释肽段重新解释的比例从 41–53% 降到 19–37%；再结合基因组比对、转录组证据和 PCR/Sanger 测序，确认了 37 个蛋白质的 40 个参考基因组插入/缺失异常，以及 9 个候选单氨基酸变异肽段。

本软件（isoPSV_multiengine）把这套方法做成三步流程：多引擎候选同池竞争与谱图级注释矫正、两分法新肽段报告，成对标记谱图核查作为终审。数据集层面的统计质控，由此延伸到单张谱图的实验物理核查。

## 目录结构

```
isoPSV\                        ← 独立软件发布目录
├── code\                      ← 全部 R 代码（总程序会自动切换到本目录）
│   ├── compat.R               base R 兼容层（str_* / melt，原始函数名）
│   ├── ReadMGF_script.R       核心库（MGF 读取/解析/打分/诱饵去除/结果文件夹/CLI 函数）
│   ├── acs_core.R             ACS 核心
│   ├── paired_params.R        谱参数提取（spectrum_params 等）
│   ├── AA_mass_dataframe.txt  氨基酸与修饰质量表（修饰等价知识库来源）
│   ├── pFind_templates\       pFind 搜库模板（open_pFind.cfg / restrict_pFind.cfg；约定目录不存在时自动生成）
│   ├── user_run.R             ★ 总程序（用户唯一入口，断点续跑）
│   ├── unlabeled_extract.R    用户第一步驱动
│   ├── user_step2.R           用户第二步驱动（诱饵去除＋五策略＋后处理＋合并＋SAP）
│   ├── user_step3.R           用户第三步驱动
│   ├── extend_report_A.R      第三步后处理（参考格式扩充）
│   ├── stageA_classify.R      阶段 A 结论分类（总程序自动调用）
│   ├── spectrum_pdf_export.R / spectrum_pdf_export_A.R
│   │                          谱图 PDF 导出（标记侧两分法／无标记侧，总程序自动调用）
│   ├── master_table_export.R  全谱主表导出（96 列，总程序自动调用）
│   ├── single_spectrum_export(_standalone).R  单谱输出（库版／独立版）
│   ├── denovo_db_lookback_test.R / normalize_terminal_loss_test.R  单元测试
│   └── punifind_rescore_convert.py  pUniFind rescore pkl→.spectra 转换器
├── web\                       网页版壳子（server.py＋index.html；双击 start_web.bat 打开）
├── config_step1_template.txt  第一步参数文件模板（用户复制填写）
├── README.md / Software_Manual.md  本说明
├── software_audit_report.txt / comment_audit.txt  软件审计报告／注释审计表
```

## pDeep2 安装（Windows，一次性）

第三步的理论谱预测使用 pDeep2（use_pdeep2=1 时需要）。本流程不捆绑 pDeep2，请按以下步骤一次性安装（约 15 分钟，约 600 MB 磁盘；安装根目录可自定，下文以 D:\Bioinformatics\pDeep2 为例，换成您的路径即可）。安装完成并接入 config_step1.txt 后，pDeep2 由软件自动执行，第三步不再停点。

1. **安装独立 Python 3.7.9**（不污染系统 Python；TensorFlow 1.15 只提供 py37 的 Windows 版）：
```bat
curl -sL -o python37.exe https://www.python.org/ftp/python/3.7.9/python-3.7.9-amd64.exe
python37.exe /quiet InstallAllUsers=0 TargetDir=D:\Bioinformatics\pDeep2\py37 Include_pip=1 Include_launcher=0
```
2. **下载并解压 pDeep 代码**（预训练模型已包含在仓库内，无需另外下载）：
```bat
curl -sL -o pDeep.zip https://github.com/pFindStudio/pDeep/archive/refs/heads/master.zip
powershell -Command "Expand-Archive -Path D:\Bioinformatics\pDeep2\pDeep.zip -DestinationPath D:\Bioinformatics\pDeep2 -Force"
ren D:\Bioinformatics\pDeep2\pDeep-master pDeep
```
3. **安装依赖**（清华镜像；protobuf 必须为 3.20.3，否则 TF 1.15 导入报错）：
```bat
D:\Bioinformatics\pDeep2\py37\python.exe -m pip install --upgrade pip -i https://pypi.tuna.tsinghua.edu.cn/simple
D:\Bioinformatics\pDeep2\py37\python.exe -m pip install "tensorflow==1.15.0" "numpy<1.19" "protobuf==3.20.3" -i https://pypi.tuna.tsinghua.edu.cn/simple
```
4. **验证安装**：
```bat
cd D:\Bioinformatics\pDeep2\pDeep\pDeep2
D:\Bioinformatics\pDeep2\py37\python.exe predict.py -e 0.3 -i Lumos -in peptide.txt -out predict_test.txt
```
成功标志：生成 predict_test.txt（MGF 格式）。运行时的 cudart64_100.dll 警告在无 GPU 电脑上属正常（自动用 CPU）。

5. **接入本流程**：config_step1.txt 填 `pdeep2_exe=D:\Bioinformatics\pDeep2`、`pdeep2_instrument=仪器类型`（QE/Lumos/Velos/Elite，默认 Lumos）。第二步完成后软件自动执行 pDeep2 生成 `result_*\merge\predict.txt` 并写入 config_step3，第三步不再停点。

**目录结构要求**：安装根下必须保持 `<安装根>\py37\python.exe` 与 `<安装根>\pDeep\pDeep2\predict.py`。若 python 3.7 装在非标准位置，可设置环境变量 `PDEEP2_PYTHON` 指向该 python.exe（软件优先读它）。实测：D:\Bioinformatics\pDeep2 安装成功，示例预测约 5 秒（CPU）。

## 网页版壳子（可选；命令行方式完全保留）

`web\` 提供一个本地网页界面（双击 `start_web.bat` 即打开浏览器）：

- **参数表单（四色分块）**：表单按大类分成四块浅色圆角卡片——①输入数据（浅蓝：项目文件夹 result_dir、psm_file、fasta、mgf）②样本标记类型（浅绿：label_type）③竞争性搜库（浅黄：auto_search＝pFind 竞争性搜库方法、pfind_exe、use_punifind；auto_search 选"0-用户手动搜库"时 pfind_exe 自动隐藏；use_punifind 选"1-启用"时显示 pUniFind 的四个选项——punifind_exe（pUniFind 安装目录选择）、punifind_rescore 与 punifind_denovo 两个三选项下拉框：0 不运行（默认）/ 1 自动运行（需填 punifind_exe，软件调官方 CLI 自动完成）/ 2 手动运行（自行提供结果，软件打印拷贝指引并停点）、punifind_python（运行 CLI 与转换器的 python，需 torch/numpy；自动 rescore 必填，留空退回手动模式））④pDeep2 预测（浅紫，每项独占一行纵向排布：use_pdeep2、pdeep2_exe、pdeep2_instrument）；所有路径均须填写绝对路径，可经系统文件对话框选择或直接输入；
- **mgf 只填一次**：阶段 B 的原始 MGF 文件夹 all_mgf 由服务器自动取自 ① 的 mgf——mgf 填单个文件夹时直接用该文件夹，填多个文件且同在一个文件夹时用该文件夹，否则留空（阶段 B 停点时用"停点续接"补 config_step1.txt 的 all_mgf）；
- **use_pdeep2 三选项**：1＝使用 pDeep2，按 pdeep2_exe 安装目录自动运行（必填安装路径，第三步不再停点）；2＝使用 pDeep2，用户对第二步产出的 `result_*\merge\peptide.txt` 手动运行 pDeep2 生成 predict.txt 后填入（pdeep2_exe 留空）；0＝不使用 pDeep2（预测指标按缺省处理）；
- **结果保存文件夹（result_dir＝项目文件夹）**：表单第一项 result_dir 指定项目文件夹——点"生成参数文件"会把 config_step1.txt 写进该文件夹；流程的全部产出写进该文件夹下的结果文件夹 `result_YYYYMMDDHHMMSS\`（仿 pGlycoQuant 时间戳命名，内含 `Competitive_search_result\`（诱饵清洗件＋SAP 索引）、`merge\`、`1.unlabeled_result\`（阶段 A）、`2.labeled_result\`（阶段 B＋主表＋PDF）），项目根另生成可复用件 unlabeled\、config_step2/3.txt（留空则退回客户端目录）。已存在 `result_*` 时取最新就地复用（断点续跑）；删除全部 `result_*` 后重跑即生成新的时间戳文件夹；
- **一键运行**：点"运行流程"＝执行 `Rscript code/user_run.R`（与命令行完全同一命令），日志实时显示，可随时"停止"；运行日志中文一行、英文一行的双语两行格式输出；
- **停点续接**：流程提示"等待输入"时，在页面选择参数文件（config_step1/2/3.txt）＋键＋值（如 config_step1.txt 的 all_mgf 或 config_step3.txt 的 predict_file）点"写入并继续"，再点"运行流程"即续接；
- "保存当前参数"把表单存入浏览器（localStorage）；"载入参数文件"弹出文件对话框选择已有参数文件（config_step1.txt 格式，UTF-8/GBK 均识别），解析其中的键值填回表单；运行历史（开始/结束时间、退出码）表格展示，并附"最近运行时间线"：每个模块的开始/结束时间、用时（墙钟）与运算（CPU 时间，两者之差即读数据与等待）、模块内各阶段耗时；pUniFind 自动 rescore/denovo 各占一行独立模块（开始/结束墙钟与用时，失败标"失败"）；运行中状态徽章显示当前所处模块；
- **结果报告（新窗口/页内嵌）**：侧栏"结果报告"页实时生成双语 HTML 报告——运行时间线、阶段 A 判定与逐谱分类分布、阶段 B 成对与冲突结论分布、错误原因统计与输出文件清单；运行结束自动刷新，"在新窗口打开报告"可整页独立查看（/api/report）；
- **停点续接**：选择参数文件后，页面自动列出该文件中的全部键（含待填写的占位键）并生成对应的数值输入控件——列表键（search_file、mgf）显示多行输入框并附"选择多个文件/选择文件夹"按钮，路径键附"选择文件/文件夹"按钮；填写后写入并重新运行流程即续接；
- 侧栏模块顺序：初始参数填写 → 运行与日志 → 停点续接 → 结果报告 → 其他信息（引用文献格式与等价命令行）。
- 服务器仅监听 127.0.0.1:8780（本地），Python 标准库实现、无第三方依赖；关闭启动窗口即停止服务器。修改 shell 或 R 代码后需重启 `start_web.bat`（浏览器刷新只重载界面）。
- 命令行方式不受影响（见下），两者可交替使用。

## 环境要求

- R ≥ 4.3（Windows 为 UCRT 版）；无外部包依赖。
- 输入数据：pFind 常规搜库结果（含 File_Name、Label 列）、pFind 衍生 MGF、参考蛋白库 FASTA、四个补充搜库结果（pFind 限定式/开放式、pUniFind 矫正/从头）、pDeep2 输出（**predict.txt**，pDeep2 仅输出此文件；预测谱的 MGF 转换由软件内部完成，不再需要 predict_HCDFT.mgf/predict_parameters.txt）。
- **引擎环境前提**：MetaMorpheus 需 .NET 10
  运行时（`winget install Microsoft.DotNet.Runtime.10`）且 `metamorpheus_exe`
  须指向官方**完整命令行包**内的 CMD.exe（含 `Data/`、`Mods/`；仅 exe/dll
  的部分解压不可用），该槽位暂冻结；引擎工作区与项目目录必须 ASCII 路径，
  启动 R 前清除 LC_CTYPE/LANG/LC_ALL；未配置 msconvert 时 Sage/MSFragger/
  MetaMorpheus 走内置 MGF→mzML 转换器（v1.0 起每谱带 centroid spectrum，
  Sage 已实测出结果）。

## 三步流程（用户视角：只调用总程序 user_run.R，断点续跑）

用户全程**只与总程序 `user_run.R` 交互**：把新生成的参数文件复制修改好后，通过 `CONFIG` 环境变量告知其（绝对）路径，然后运行总程序；总程序自动判断哪些步骤已完成、哪些需要用户补参数，并在需要用户输入时打印双语提示后正常退出——用户填好参数文件后**重跑同一条命令**即可续接，已完成的步骤自动跳过。核数与并行数由软件自动计算；SAP 结果由软件在第二步自动从补充搜库结果（开放式搜库＋pUniFind）中提取，均无需用户提供选项。

| 步骤 | 跳过条件 | 停点（打印提示后退出，填完重跑即续接） |
|---|---|---|
| 第一步 | `unlabeled\step1_done.txt`（第一步完成标记）已存在 | 第一步完成后：请填写自动生成的 `config_step2.txt` 的 search_file×4 |
| 第二步 | `result_*\merge\peptide.txt`（pDeep2 输入）已存在 | ① use_punifind=1 且启用的 pUniFind 项（punifind_rescore/denovo＝1 自动或 2 手动）产出未齐：第二步停在待办提示（**不生成 peptide.txt**）——自动模式提示"检查日志后重跑"，手动模式打印拷贝指引（GUI 产出拷到指定目录）；完成后直接重跑总程序续接（无需删除文件）；② 第二步完成后：请填写自动生成的 `config_step3.txt` 的 predict_file |
| 第三步＋后处理 | 不跳过（每次运行） | use_pdeep2=1 且 config_step1.txt 填了 pdeep2_exe 时：软件自动执行 pDeep2 生成 predict.txt、无停点；use_pdeep2=1 但未填 pdeep2_exe 时：提示补 pdeep2_exe（或改 use_pdeep2=2 手填 predict_file）；use_pdeep2=2 时：请填写 config_step3.txt 的 predict_file（对 `result_*\merge\peptide.txt` 手动运行 pDeep2 生成的 predict.txt）；use_pdeep2=0 时无此停点（预测指标按缺省处理） |
| 阶段 A 分类 | 不跳过（每次运行，纯表格运算） | 输出 `result_*\1.unlabeled_result\classification_A.tsv`：全部谱图逐谱结论状态（英文：no competition / competition but original correct / original result wrong＋error_reason＋modification_diff / chimeric spectrum / ambiguous: poor spectrum quality / ambiguous: identical fragment ions） |
| 阶段 B（仅标记样本） | `result_*\2.labeled_result\paired_report_<标记类型>.tsv` 已存在；无标记样本自动跳过 | 标记样本且 config_step1.txt 未填 all_mgf 时：请填写 all_mgf（原始 MGF 文件夹，含无标记＋标记谱） |
| 冲突报告（仅标记样本） | 阶段 B 之后每次执行（纯表格运算，无停点） | 输出 `result_*\2.labeled_result\conflict_report_<标记类型>.tsv`：以 report_A.tsv 的冲突行为基础合并阶段 B 新增列，逐谱结论状态（英文：competition but original correct／original result wrong／chimeric spectrum／ambiguous: poor spectrum quality／ambiguous: vote tie／ambiguous: identical fragment ions／not pair-checked），错误行附图 4 口径错误原因与修饰差异细节（如"original extra: 12,Oxidation[M] ｜ winner extra: 14,Oxidation[M]"）；阶段 B 列为 paired_label_cosine_B/choosen_paired_label_cosine_B/paired_state_B（Paired/Single）并附 winner_sequence/winner_modification（投票胜选候选），同步生成 `conflict_report_<标记类型>_summary.txt`：谱图级按 Single/Paired 分组统计（结论状态、错误原因）＋肽段层面三个 panel（Single 谱贡献肽段、Paired 谱贡献肽段、合并去重肽段，含每类结论的去重肽段数）；网页端结果报告页含同样的分组与肽段三 panel 表格 |

### 命令行操作手册（Windows cmd）

1. **打开命令行**：按 `Win + R`，输入 `cmd`，回车。
2. **进入工作目录**（本目录＝含 `code\` 与 `config_step1_template.txt` 的文件夹，复制粘贴时换成实际路径）：
   ```bat
   cd /d "isoPSV"
   ```
3. **准备参数文件**：把 `config_step1_template.txt` 复制一份并改名（如 `config_step1.txt`），用记事本打开，把 `<绝对路径>` 等占位符全部替换为实际值（**路径一律写绝对路径**），保存。
4. **复制粘贴运行**（命令中的 `set` 只需在窗口打开后执行一次，之后每次运行只粘贴最后一行即可；已 `cd` 到参数文件所在目录时，CONFIG 直接写文件名即可，也可写绝对路径）：
   ```bat
   set RSCRIPT=D:\Bioinformatics\RStudio\R-4.3.1\bin\x64\Rscript.exe
:: CONFIG 可省略（默认用本目录的 config_step1.txt；用其他文件名时才需要 set CONFIG=...）
   "%RSCRIPT%" code\user_run.R
   ```
   若中文显示为乱码，先运行 `chcp 65001` 再运行上面的命令（不影响功能，仅影响显示）。
5. **程序停点续接**：
   - 屏幕提示"请填写 `config_step2.txt` 的 search_file"→ 用记事本打开该文件，把 4 个 `search_file` 占位行改为四个补充搜库结果的**绝对路径**（每个一行），保存 → 重新粘贴运行第 4 步最后一行命令。
   - 屏幕提示"请填写 `config_step3.txt` 的 predict_file"→ 打开该文件，把 `predict_file` 占位行改为 pDeep2 的 predict.txt **绝对路径**，保存 → 重新粘贴运行。
   - 标记样本还会提示"请填写 `config_step1.txt` 的 all_mgf"→ 把 `all_mgf` 行改为原始 MGF 文件夹（含无标记＋标记谱）的**绝对路径**，保存 → 重新粘贴运行（阶段 B 成对匹配）。无标记样本自动跳过阶段 B。
6. **完成与断点续跑**：全部完成后横幅打印"全部完成 all done"，最终结果在结果文件夹 `result_YYYYMMDDHHMMSS\` 下——`1.unlabeled_result\`（阶段 A：report_A.tsv、classification_A.tsv、谱图 PDF）与 `2.labeled_result\`（标记样本：paired/conflict 报告、master_table.tsv、两本谱图 PDF）。任何时候中断后重跑同一命令即可续接——第一步/第二步只要中间产物存在就自动跳过，不会重复计算；第三步每次重新计算最终报告。

各步驱动也可单独运行（调试用），参数文件链式自动生成：

| 步骤 | 参数文件 | 命令 | 输出 |
|---|---|---|---|
| 第一步 | 用户填写 `config_step1.txt`（或复制模板） | `Rscript unlabeled_extract.R` | `unlabeled\` 逐 raw 纯无标记 MGF（`<raw>_HCDFT.mgf`，按标题 raw 前缀重切；**不产出合并文件**）＋`step1_done.txt` 完成标记；**自动生成第二步参数文件** |
| 第二步 | 用户填写第一步生成的 `config_step2.txt` 的 search_file×4 | `Rscript user_step2.R` | `result_*\Competitive_search_result\`（诱饵清洗件＋SAP 索引）与 `result_*\merge\` 的 1.combine_result.spectra＋peptide.txt；**自动生成第三步参数文件** |
| 第三步 | 用户填写第二步生成的 `config_step3.txt` 的 predict_file | `Rscript user_step3.R` | `result_*\1.unlabeled_result\` 的 report_A.tsv / summary_A.txt / scan_details_A.tsv |
| 后处理 | 同第三步参数文件 | `Rscript extend_report_A.R` | report_A.tsv 扩充参考格式列 |

## 参数文件（键值制表符＋# 注释，pFind.cfg 风格）

- **用户填写的参数文件内所有路径一律写为绝对路径**；软件自动生成的参数文件同样写绝对路径（同一口径）。`mgf`/`search_file` 每行一个，文件路径直接收、**文件夹路径展开为其中全部匹配文件**；软件对参数文件只读（重新生成时保留用户已填条目，绝不覆盖成占位符）。

## 标记类型（label_type）与酶切类型（enzyme）

**label_type** 是多行键，每行写一个标记；空着或整键不填，即为无标记样本。两类标记可单独配置，也可同时配置。

第一类：全原子稳定同位素替换。整行填 N15 或 对应 15N / 13C 全标记。旧单值写法在任何版本都兼容。

第二类：修饰形式标记。每行由四个竖线分隔的字段自定义：

　　`名称|位置码|氨基酸集合|元素组成`

- 位置码五种：ANY（任意位置）、PEP_N（肽段N端）、PEP_C（肽段C端）、PRO_N（蛋白N端）、PRO_C（蛋白C端）；中文"任意位置 / 肽段N端 / 肽段C端 / 蛋白N端 / 蛋白C端"同样可写。
- 氨基酸集合：`*` 表示全部氨基酸，或写单字母组合（如 STY、KR）。
- 元素组成：modification.ini 风格，如 H(1)O(3)P(1)；计数可带负号，支持同位素元素 2H、13C、15N、18O、34S。标记质量由软件按内置元素表计算，不需要手工查数填质量。

示例：

　　`Phospho|ANY|STY|H(1)O(3)P(1)`　　　　　磷酸化，任意位置的 S/T/Y
　　`NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)`　NeuCode 赖氨酸标记（8.014199 Da，等价旧预设 Label_13C(6)15N(2)[K]）
　　`Dimethyl|PEP_N|*|H(4)C(2)`　　　　　　　二甲基，肽段N端、任意残基

旧预设名（Label_13C(6)[K]、Dimethyl[AnyN-term] 等按 AA_mass_dataframe 查表的写法）继续可用，新旧写法可混在同一组里。修饰名称中不要使用竖线、方括号、分号、逗号。

网页表单里第二类每行四个格子（名称、位置下拉、氨基酸"全部/自定义"、组成），另备六个常用预设一键填充；混合样本每行添加一个标记即可。

与标记配合使用的还有 unlabeled_channel 键：1＝样本含无标记通道（默认，第一步按 Label 提取无标记谱）；0＝不含（NeuCode 等全标记样本），此时第一步提取全部鉴定谱，阶段 A 把配置的标记当作固定修饰搜库。

引擎侧的行为：comet 的端部标记折算为裸 n/c（任意残基端）；MSFragger 折算为 `[^` / `]$`（任意残基端）或 n/c 前缀（特定残基端）；蛋白端两种引擎没有原生的可变修饰语法，按肽段端近似注入并在控制台提示；sage 的 variable_mods 只认残基位点，端部标记不注入 sage，其余引擎照常携带。

**enzyme** 填酶切类型，默认 trypsin。可选：trypsin / lysarginase / lysc / lysn / argc / gluc / aspn / chymotrypsin（糜蛋白酶）/ none（无酶切）。Trypsin 与 LysargiNase 的切点均受脯氨酸阻断（P 位切不开）。酶切类型会改写各引擎的酶参数：comet 按其参数文件自带的 COMET_ENZYME_INFO 表按名取编号，表内没有的酶（如 LysargiNase、GluC）按同格式追加编号；MSFragger 改写酶 1 的四行定义（none → nonspecific）；sage 改写 cleave_at / c_terminal / restrict。仅在该键显式填写时改写，键留空则引擎模板保持原样（模板即真相）。酶参数带侧车指纹，重跑时参数未变即复用引擎结果、不再重搜。

- `config_step1.txt`：psm_file、fasta、label_type（多行键，口径见上节）、mgf、all_mgf（原始 MGF 文件夹，含无标记＋标记谱，阶段 B 成对匹配用；无标记样本留空；网页壳子自动取自 mgf）、auto_search（第二步搜库方式：1=软件自动调 pFind 搜库（config_step1.txt 填 pfind_exe＝pFind 安装路径，版本目录如 pFind3.1 自动补 bin）；0=用户手动搜库后填 search_file）、use_punifind（1=第二步使用 pUniFind rescore/denovo；0=不用，全 pFind 流程，忽略 pUniFind 输入）、punifind_exe（pUniFind 安装的 punifind 目录，含 official_score_workflow.exe / official_denovo_workflow.exe、ckpts、consts；自动模式下必填，留空即手动模式）、punifind_rescore（开放式搜库结果重新打分校正：0=不运行（默认）；1=自动运行——需 punifind_exe，软件调官方 CLI 完成，路径含中文且直接失败时自动改 C:\temp 英文临时目录重跑（结束后删除临时目录）；2=手动运行——用户在 GUI 完成 rescore 后把产出拷入指定目录，软件打印拷贝指引并停点）、punifind_denovo（从头测序：0/1/2 同上，默认 0）、punifind_python（运行 pUniFind 官方 CLI 与 pkl→.spectra 转换器的 python 可执行文件，需 torch/numpy——自动 rescore 必填，未配置时退回 GUI 导出待办口径）、use_pdeep2（是否使用 pDeep2 理论预测谱：1=用，按 pdeep2_exe 安装目录自动运行；2=用，用户手动运行 pDeep2 生成 predict.txt 后填 predict_file；0=不用，预测指标按缺省处理——缺失不填 0、predict_missing 标注、投票按谱降级）、pdeep2_exe（pDeep2 安装路径，use_pdeep2=1 时必填：软件对第二步 peptide.txt 自动执行 pDeep2 生成 `result_*\merge\predict.txt` 并写入 config_step3——第三步不再停点；use_pdeep2=2 时留空）。python 3.7 路径两级解析：环境变量 `PDEEP2_PYTHON`（不同电脑可覆盖）→ 安装文档标准布局 `<pdeep2_exe>\py37\python.exe`；均不存在时报错提示、pdeep2_instrument（仪器类型，QE/Lumos/Velos/Elite，默认 Lumos）。
- **mgf 只在第一步输出一次，第二步起只读不写**：第一步按标题 raw 前缀把无标记谱图重切为 `unlabeled\<raw>_HCDFT.mgf` 纯逐 raw 文件（pParse 共洗脱拆分的输入文件可能 raw 混杂；pUniFind rescore 按 qry.res 的 raw 名加载 `<raw>_HCDFT.mgf` 建 TITLE 字典，混杂文件会报 KeyError），不产出合并 unlabeled.mgf；pFind 自动搜库直接用这批文件（msmsnum=N），pUniFind GUI 的 mgf 目录也选 `unlabeled\`。搜库结果已存在时自动复用（删除项目 `result\` 可强制重新搜库）。
- **pFind 搜库线程写死为总核数-1**：自动搜库时软件把运行 cfg 的 thread 改写为总核数-1（覆盖模板 cfg 里的 thread 值），留一核给系统与交互（GUI/网页壳子/其他程序），搜库期间机器不会被占满卡死；实际线程数打印在"运行 pFind"日志行。
- **阶段 B 轻谱定位（raw+scan 兜底）**：搜库谱图名不会变，title 格式为 raw.scan.scan.价态.拆分序号.dta——阶段 B 先按标题精确匹配，失败时按标题前两段（raw 文件名＋scan）在原始 MGF 中查找，共洗脱重切导致的价态/拆分序号差异不影响定位。MGF 并行加载器的分块索引错误（大文件与小文件混载时中小文件丢份、曾致阶段 B 丢失约半数谱块）已修复，重跑后冲突谱全部进入成对核查。
- **混合谱判定**：对一条谱的全部不同竞争序列两两比较，各自特异子离子比例（只被该序列匹配到的 b/y 槽位数÷该序列匹配到的槽位数）的**最小值 ≥ 10%** 即判为混合谱（分类报告与冲突报告的 chimeric 列按此口径、经离子匹配计算）。
- **pUniFind 三项运行方式**（use_punifind=1 时，punifind_rescore/punifind_denovo 三值）：第二步先备好两个项目——rescore 项目＝pFind_open（`.pac` 在项目根、`mgfs\` 为纯逐 raw 的 `<raw>_HCDFT.mgf`、pFind 逐 raw 搜库结果在 `result\`）；denovo 项目＝`pUniFind_denovo\`（mgf 目录选 `unlabeled\` 纯逐 raw）。随后按三值执行：**1 自动运行**——软件调官方 CLI（official_score_workflow.exe / official_denovo_workflow.exe，含完整参数集；路径含中文且直接失败时自动改 C:\temp 英文临时目录重跑，结束后删除）。rescore 与 denovo 都按**单 MGF（逐 raw）分片**运行：mgf 目录多于一个文件时自动逐 raw 建片（`chunks\<raw>\`，片级完成标记），每片独立推理后合并产物（denovo 拼接 CSV、rescore 合并谱表），内存需求只与单个 raw 的谱数挂钩；某片失败时日志给出该片身份、日志路径并附该片子日志尾（真实报错直接可见），直接重跑总程序即从失败片续跑（已成功片跳过）。分片运行配三个 GPU 供给旋钮（v1.0）：`punifind_batch_size`（推理批大小，默认 64，大数据集可试 128–256）、`punifind_chunk_parallel`（片级并行路数，默认 1＝顺序即停；16 GB 显存建议 2–3，每片独立进程约 3 GB、各路 MASTER_PORT 自动错开；并行路径跑完全部片再汇总报错，成功片标记保留、失败片可续跑）、`punifind_chunk_cleanup`（1＝合并成功后清理片内 mgf/qry/data_tmp，保留完成标记与合并产物，默认 0）。完成后自动登记为策略 4/5，运行时间显示在"最近运行时间线"（独立模块行，失败标"失败"并打印日志提示后停点）；产物已存在时跳过重复运行（与 pFind 搜库复用同口径，删除对应 result\ 下产物可强制重跑）；**2 手动运行**——软件打印拷贝指引并停点：rescore 输出 `<项目名>fdr0.01_pUniFind.spectra` 拷入 `pFind_open\result\`，denovo 输出 `pUniFind_result\test_5_0.csv` 拷入 `pUniFind_denovo\pUniFind_result\`；**0 不运行**——该项完全跳过（对应策略列留空）。**启用的项产出未齐时第二步不合并、不生成 peptide.txt**（停点收尾），完成后**直接重跑总程序**即续接：软件自动发现产出登记为策略 4/5、pFind 搜库结果自动复用（不重复搜库）、重新五策略合并——无需删除任何文件。
- - **逐 raw 并行搜库键**：`pfind_workers`
  （config_step1，默认 1=逐 raw 串行；如 6 个 raw 设 6 则 pFind open/restrict
  墙钟约缩到 1/6）与 `<engine>_workers`（config_step2：msfragger_workers、
  comet_workers、sage_workers、metamorpheus_workers，默认 1=串行）——并行时
  线程与内存（MSFragger -Xmx）按 worker 数自动摊薄，产物按 raw 顺序合并、
  与串行结果集合等价。
- **python_exe 键**：Python 解释器绝对路径（config_step1，
  可留空）。第二步重编号步骤 renumber_mgf.py 需要 Python；留空＝按环境目录（PATH）
  顺序探测 python / py -3 / python3。用户的 Python 若没有写入环境目录，请在此填
  完整路径（如 D:\Python38\python.exe），软件优先使用本键；所填路径不存在时报
  提示并回退 PATH 探测。网页参数表单的 pFind / pUniFind 生态区块有对应输入框
  （可经文件对话框选择，生成 config_step1 时写入本键）。
- **网页多引擎策略下拉**：策略配置区下拉仅列
  外部引擎（MaxQuant / MSFragger / Sage / Comet / X!Tandem / MetaMorpheus /
  PEAKS / Mascot / pNovo）。pFind 家族策略 b–e（open / restrict / rescore /
  denovo）由第一步 config_step1 的 auto_search 与 use_punifind 键管理，在
  pFind / pUniFind 生态区块配置，不再经此下拉添加——避免与自动注册槽位重复
  冲突（此前显式 pFind 槽与 auto_search 同启曾致衣藻全引擎臂报错）；浏览器
  localStorage 里存留的旧 pFind 行会自动被过滤。
- **MaxQuant 与 FAIMS 拆电压数据**：FAIMS
  480 采集通常把数据按补偿电压拆分为多个 raw（文件名 `_45`/`_60`/`_65`
  后缀）。这类拆分 raw 直接交给 MaxQuant 自动搜库会在启动时崩溃（逐文件
  FAIMS 组参数越界，`GetGroupParamsForFile` IndexOutOfRange）。**MaxQuant
  必须使用未拆分电压的原始 raw**（通常在原始数据目录的 `original\` 子文件
  夹内、无电压后缀），由 MaxQuant 自行处理 FAIMS。建议流程：用户在
  MaxQuant 中手动搜库原始 raw（`maxquant_mode` 留 0，不启用自动搜库），
  再把产出的 `msms.txt` 作为策略结果挂入（config_step2 增加一行
  `search_file\t<msms.txt 路径>` 与 `search_file_type\tmaxquant`，软件的
  msms.txt 适配器会转为交换格式进入合并）。其余引擎（Comet / MSFragger /
  Sage / MetaMorpheus）不受影响，可直接吃 unlabeled MGF 自动搜库。
- `config_step2.txt`（第一步自动生成）：base_config、search_file×4（占位行由用户填写原始版搜库结果；I→L/N→GG 与 4 种修饰校正由软件执行）。**自动搜库（可选）**：把 pFind 开放式/限定式模板参数文件按约定目录放在工作目录下（`pFind_open\param\pFind.cfg` 与 `pFind_restrict*\param\pFind.cfg`，即模板参数文件所在目录名为 pFind_open / pFind_restrict*，输出目录与目录名一致），软件自动发现并写入 pfind_cfg_open/pfind_cfg_restrict 键；第二步由软件自动调 pFind 完成开放式与限定式搜库（模板中的 msmspath1/fastapath/outputpath/outputname 自动替换为本流程的 unlabeled.mgf 与输出目录），无需用户手工搜库。 **搜库侧**：open 搜库的运行参数自动把 modpath 指向随客户端分发的 `code\pUniFind_mods.ini`（pFind modification.ini 格式的兼容修饰库，由生产 open 结果提取、经 pUniFind rescore 验证的修饰名构建），open 搜库只在已知修饰集合内解析（如 Methyl[D] 类问题在源头消除，不动 pFind 原始 modification.ini）；搜库后仍按兼容集合过滤一遍作为保险丝（pFind 动态生成的库外名）。
- `config_step3.txt`（第二步自动生成）：base_config、merge_result、weights_A、threshold_A、predict_file（占位行由用户填写 predict.txt 路径）。
- **SAP 清单无需用户提供**：第二步自动从四个补充搜库结果（开放式搜库与 pUniFind 结果）中提取 SAP 替换行（生产 1.SAP_extract.R 移植），写入 `merge\SAP_spectra_index.txt`，第三步后处理自动读取用于 Error_reason 的 SAP 判定。

## 运行约定

- 命令行环境变量传路径：`CONFIG`、`CONFIG2`、`CONFIG3`。**Windows 中文路径不进命令行参数**，一律走环境变量。
- **核数与并行数由软件自动计算**（第一步按内存/最大 MGF 评估读取核数；第三步扫描并行数＝核数−1），用户无需提供任何并行选项。
- 总程序自动切换到 `code\` 目录执行各步（单步驱动需手动 cd）；单步运行前 `unset LC_CTYPE LC_ALL`（zcode 终端导出的 C.UTF-8 在 Windows 无效，会导致中文路径读取失败；普通 cmd 窗口无此问题）。R 4.5 下 `normalizePath` 在该环境中会腐蚀路径中的中文——user_run.R 已加字符串回退兜底但仍推荐清 locale 或经 cmd.exe 启动。
- 第一步运行时自动资源评估（systeminfo 可用内存＋核数＋最大 MGF 大小 → max_cores，仅打印与使用）；第三步并行逐 scan 评估，控制台显示进度条（已处理谱图数/总数百分比）；每步主要阶段双语横幅打印到控制台。
- 全部输出为 **tab 分隔 txt**（无 xlsx 交互）；蛋白位置列写为 `'174-190`（前导撇号防止 Excel 将 "XX-XX" 识别为日期）。
- **大规模数据的两处内建优化，无需任何开关**：①诱饵清洗为行级过滤实现（未命中行原样透写），百万行级 open 结果表的清洗为分钟级；②MGF→mzML 转换产物缓存于项目级 `mzml_cache\`——msfragger/sage 等多个引擎消费同一批 MGF 时只转换一次，第二处起直接复制；产物已新于源 MGF 时零成本跳过。另可选：安装 ProteoWizard 并在 config 填 msconvert_exe 键，转换即走外部工具原生速度。
- **绿色运行时与一键启动（v1.0）**：包根目录双击 `isoPSV.exe`（或
  `start_isoPSV.bat`）即启动本地网页界面并自动打开浏览器，关闭窗口即停止
  服务。启动器优先使用包内 `runtime\python`（官方 Python 3.12.10 embeddable，
  含 LICENSE）与 `runtime\R`（R 4.5.1 基础包＋推荐包裁剪版，约 158 MB），
  未安装 Python/R 的电脑可直接运行；检测到内置 R 时自动加入 PATH 并设置
  RSCRIPT 环境变量，网页壳子的"运行"按钮经内置 R 执行三步流程。原
  `web\start_web.bat` 保留为转发入口。命令行方式：
  `runtime\R\bin\x64\Rscript.exe code\user_run.R`（CONFIG 指向
  config_step1.txt）。pFind/pUniFind/pDeep2 与各外部引擎仍为独立软件，需
  用户按需自备。
- **两个新模式键（默认关闭＝现状）**：
  - `candidate_scope`（filtered 默认 / top5）：top5＝忽略各搜库软件 1% FDR
    过滤，每个引擎每张谱取引擎自身排序前 5 个肽段+修饰组合进竞争池；判定
    指标仍由软件从谱图独立计算、去诱饵清洗照常，小样本不依赖软件 FDR 也可
    用。合并表变为逐谱 1..K 行（K=各槽候选数最大者、封顶 5），下游报告与
    判定口径不变。
  - `proteogenomics_mode`（0 默认 / 1）：1＝纯蛋白质基因组学模式——只对含
    非参考库候选的谱做竞争验证（新肽段/新基因鉴定），已注释谱标记
    annotated_skipped 并跳过错误校正（不参与判对/判错与五类错误统计），
    计算与输出聚焦新肽段证据；定侧口径与两分法报告一致（Proteins 任一蛋白
    命中用户 FASTA 即"已注释"）。

## 输出文件说明（格式、内容与表头缩写含义）

全部表格输出为 tab 分隔 txt（UTF-8；纯 ASCII 摘要文件按网页渲染中英标题），
MGF 谱峰数值为完整十进制（不使用科学计数法，pFind 不识别 "E" 写法）。
蛋白位置列写为 `'174-190`（前导撇号防 Excel 把 "XX-XX" 识别为日期）。

### 第一步（项目根 unlabeled\）

- `<raw>_HCDFT.mgf`：逐 raw 的无标记（轻标）谱图，pParse 标题命名
  `raw.scan.scan.charge.fraction.dta`；每谱块 BEGIN IONS/TITLE/CHARGE/RTINSECONDS/
  PEPMASS/峰列表（m/z 强度）。第二步起只读，不再改写。
- `step1_done.txt`：第一步完成标记（总程序断点续跑依据）。

### 第二步·自动搜库（pFind_open\result、pFind_restrict\result）

- `pFind-Filtered.spectra`：1% FDR 过滤后的鉴定表。表头：File_Name（谱图标题）、
  Scan_No、Exp.MH+、Charge、Q-value、Sequence、Calc.MH+、Mass_Shift(Exp.-Calc.)、
  Raw_Score、Final_Score、Modification（"位点,修饰名[残基];"）、Specificity、
  Proteins（"ACC/" 直接相邻拼接）、Positions（"位置,前残基,后残基/"）、
  Label（1|＝无标记、2|＝标记）、Target/Decoy、Miss.Clv.Sites、
  Avg.Frag.Mass.Shift、Others。
- `pFind.spectra`：未过滤的全部 PSM；`pFind.summary`：搜索统计（ID rate、
  修饰/长度/电荷分布）；`<outputname>.F*.L*.qry.res`：pFind 二进制结果
  （供 pUniFind rescore 输入）；`1.aa/1.mod/2.aa/2.mod`：修饰辅助文件；
  `backup_pFind.cfg`：本次搜索实际使用的配置备份。

### 第二步·合并（result_*\merge\）

- `revised\`：各策略经诱饵清洗与 revise 修正后的中间件（命名
  `<策略名>_<源名>_revised.spectra`）。
- `1.combine_result.spectra`：多策略合并表（逐谱一行）。策略列按单字母槽位：
  `a.Sequence/a.Modification`、`b.Sequence/b.Modification`……（每策略两列）；
  `Peptide_state`＝按序列分组的策略字母串（如 "abd;c"），`Modification_state`
  ＝按修饰分组的字母串；其余列与 pFind-Filtered.spectra 同义（File_Name、
  Scan_No、Exp.MH.、Charge、Q.value、Calc.MH.、Mass_Shift.Exp..Calc..、
  Raw_Score、Final_Score、Specificity、Proteins、Positions、Label、
  Target.Decoy、Miss.Clv.Sites、Avg.Frag.Mass.Shift、Others）。
- `peptide.txt`：pDeep2 输入。表头 `peptide/modification/charge`（tab 分隔），
  SAP 型替换修饰（X->Y）已代入序列后再交给 pDeep2；越界修饰位点已过滤。
- `predict.txt`：pDeep2 预测输出（谱块格式：BEGIN IONS/TITLE=序列|修饰|价态/
  CHARGE/pepinfo/b+1,b+2,y+1,y+2 概率数组/PEPMASS/峰行 "m/z 强度 离子名"）。
- `SAP_spectra_index.txt`：SAP 单独记录（原始序列＋"pos,XXX->YYY[site]" 修饰
  与 revised_peptide/ref_protein 映射），pDeep2 只用替换后序列，下游可还原。

### 第三步·阶段 A（result_*\1.unlabeled_result\）

- `report_A.tsv`：逐谱一行——candidates（肽段＋修饰组合，"|" 拼接）、逐候选
  覆盖度/PSM_by/PSM_all/预测余弦/特异子离子比例（"|" 拼接，第 1 个候选留空）、
  best_group/best_sequence、spec_ratio_2..5（策略 1 vs 2–5 的混合谱特异子离子
  比例＝两肽各自特异子离子÷两肽总子离子取最小值）、decision_A
  （consistent/candidate_misidentification）、candidate_ACS_A/best_ACS_A
  （权重待成对标记部分完成后标定）。
- `classification_A.tsv`：阶段 A 逐谱判定结论。
- `scan_details_A.tsv`：每 scan 逐候选参数明细（阈值敏感性分析免重读谱图）。
- `summary_A.txt`：阶段 A 统计。

### 标记样本·阶段 B 与主表（result_*\2.labeled_result\）

- `paired_report_<标记类型>.tsv`：不一致谱的成对标记匹配结果——逐候选
  match_state（Matched/Single）、轻/标记谱覆盖度与 PSM_by/PSM_all、
  label_File_names、paired_label_cosine、choosen_paired_label_cosine、
  paired_state（Paired/Single，余弦阈值 0.7）、best_group。
- `conflict_report_<标记类型>.tsv` 与 `conflict_report_<标记类型>_summary.txt`：
  冲突谱全参数与统计。
- `master_table.tsv`：**全谱主表**（判校正行全列填写，无差异行结论列留空、
  corrected_* 复制候选 1）。表头缩写含义（96 列，多引擎版随策略数扩展
  cand6..candN 每组 15 列）：

| 列名 | 含义 |
|---|---|
| File_Name | 谱图标题（.dta 命名 raw.scan.scan.charge.fraction） |
| cand1_protein / cand1_location | 原始鉴定（候选 1）所在蛋白 accession；蛋白定位（前残基[0 起位置]后残基，如 R[215]T） |
| cand1_sequence / cand1_modification | 原始鉴定肽段序列；修饰串（"位点,修饰名[残基];"） |
| cand1_ul_by / cand1_lb_by | 候选 1 无标记 / 标记 PSM 打分（b/y 离子） |
| cand1_ul_all / cand1_lb_all | 候选 1 无标记 / 标记 PSM 打分（全部离子） |
| cand1_pcos | 候选 1 成对标记谱余弦（paired cosine） |
| cand1_predict | 候选 1 pDeep2 预测谱余弦（predict cosine） |
| cand1_paired_state | 候选 1 成对标记状态（paired/single） |
| cand1_unlabeled_coverage | 候选 1 无标记谱子离子覆盖率（%） |
| cand1_labeled_scan | 候选 1 对应的标记谱标题（paired 时） |
| cand1_score | 候选 1 加权总得分（Score=Σwⱼ·pⱼ；权重定稿值 w=(ul_by 0.6200, lb_by 0.0340, ul_all 0.1186, lb_all 0.1701, pcos 0.0573)，predict 按 θ=0.30 救援、不计入主打分） |
| candN_*（N=2..5） | 候选 N 的同组 14 列（sequence/modification/ul_by/lb_by/ul_all/lb_all/pcos/predict/paired_state/unlabeled_coverage/labeled_scan/score） |
| candN_delta_score | 候选 N 与候选 1 的得分差（ΔScore） |
| candN_pred_rescue | 候选 N 是否触发预测救援（δ_pred≥θ＝0.30 时 TRUE） |
| candN_specific_ratio_vs_cand1 | 候选 N 相对候选 1 的特异子离子比例（混合谱判定） |
| rev_idx | 被判定选中（替换原始鉴定）的候选编号 |
| max_delta_score | 候选间最大得分差 |
| tau | 判定阈值 τ |
| gap | τ−ΔScore（难区分带判据 |gap|≤0.015） |
| call_base_on_score | 按得分差判定（ΔScore≥τ）的标记 |
| delta_pred | 原始鉴定与对手的预测得分差（δ_pred） |
| pred_rescue | 预测救援标记 |
| call_base_on_score&pred | 最终判需校正标记（得分差或预测救援任一触发） |
| ions_identical | 子离子完全一致标记 |
| ions_small_diff | 子离子微小差异标记 |
| coverage_below_20 | 覆盖率低于 20% 标记 |
| min_specific_ratio | 混合谱特异子离子比例最小值（chimeric 判定） |
| corrected_protein / corrected_location | 校正后（判定为正确一方）的蛋白 / 定位 |
| corrected_sequence / corrected_modification | 校正后的肽段序列 / 修饰串 |
| corrected_labeled_scan | 校正后对应的标记谱标题 |
| diff_type | 六类错误类型：Chimera (same charge)、Chimera (different charge)、N-terminal non specific cleavage、Different modification、Modification location changed、Same weight AA pairs、Peptide length less than 7 AAs、SAP |
| merged_error_type | 五类合并口径（混合谱含同价/异价、序列替换含等重对与长度类、修饰改变、SAP、N 端非特异）——用于 Fig.7d 柱状图与网页报告 |
| mod_changed | 修饰变化对比（原始修饰串＝校正后修饰串） |
| mod_type_changed | 修饰类型改变 |
| aa_changed | 氨基酸改变 |

- `correction_paired_summary.txt`：判校正行的成对标记证据统计（总数、paired
  数、余弦>0.7 数，按五类错误分列）。
- `1.annotation_correction_report.tsv`：**两分法第 1 部分（注释校正报告）**——
  原始鉴定蛋白在用户 FASTA 内的行；判校正行（diff_type 非空）＝错误匹配清单。
- `2.novel_peptide_competition_report.tsv`：**两分法第 2 部分（新肽段竞争报告）**——
  原始鉴定蛋白不在用户 FASTA 内的行。除主表列外追加（**master_table.tsv 本身
  同样携带这五列**，已注释行留空）：
  - `credible_novel`（TRUE＝可信：竞争保留且非混合谱且无任何存疑。**不要求
    存在竞争候选**——真实新肽段可能在已注释库中无对应物，无竞争对象不构成
    证据不足；设计约定）；
  - `n_frag`：新肽段匹配子离子数（b/y 位点 1..n−1 中有命中的位数）；
  - `closest_annotated_competitor`：打分最相近（ΔScore 最大＝最强）的已注释
    竞争候选，格式 `序列|修饰|ΔScore`（Δ＝竞争者分数−新肽分数，正＝竞争者
    分数更高）；
  - `novel_status`：状态分类——空＝可信；**判错行写具体校正类型（diff_type，
    如 "Chimera (same charge)"＝与已注释肽段混合谱、 "Different modification"
    ＝修饰改变、"N-terminal non specific cleavage"＝N 端非特异）**；
    `close_annotated`＝无法判定但有已注释有力竞争（最强已注释竞争候选
    |ΔScore|≤τ，无需判错即认为新肽段高错误概率、不能保留）；
    `chimeric`/`ions_identical`/`consecutive_missing`＝其余存疑；
  - `novel_doubt`（存疑详情：判错行＝具体校正类型；`chimeric`、
    `close_annotated:<序列>|<修饰>|<Δ>`、`ions_identical:<序列>`＝与已注释
    竞争候选匹配离子集完全一致、`consecutive_missing:<n>`＝连续 ≥n 个位点
    无子离子命中；`low_score:<分>`＝PSM 分数低于 0.2 的低质量谱卡值
）。
- `novel_spectra_pBuild.spectra`：**pBuild 谱图清单**（pFind-Filtered.spectra
  同格式）——可信新肽段的无标记行（Label=1，字段取自原始搜库结果）＋
  成对标记行（Label=2，File_Name/电荷/前体质量换成 15N 标记对应体），
  供 pBuild 画谱。
- `two_part_report_summary.txt`：两分法摘要（part/rows/headline），网页报告
  按此渲染两个面板。
- `raw_loaded_<标记类型>.rds/.sig`：阶段 B 中间缓存（可删除）；`worker_error_log.txt`
：并行 worker 错误日志（正常情况下为空）。
- **谱图 PDF 导出键 `pdf_on`**：config_step1 键
  pdf_on（1=生成错误匹配/存疑谱图 PDF、0=不生成，默认 0），环境变量
  PDF_ON 兜底（任一为 1 即开启）；网页参数填写区有下拉选项。
- **外部引擎自动运行**：
  config_step2 键 `<engine>_exe`／`<engine>_params`／`<engine>_mode`（0 手动、
  1 自动）——msfragger（java -jar）、comet（exe -Pparams）、sage（msconvert
  转 mzML＋sage.exe，需另配 `msconvert_exe`）逐 unlabeled MGF 自动重编号→运行
  →适配器转换→拼接 `<engine>_auto.spectra` 进对应槽位（产物新于 MGF 时复用；
  参数模板在 code/engine_templates/，数据库路径运行时按项目 FASTA 重写）。
  msfragger_boost（FragPipe headless 工作流）与 maxquant（mqpar.xml）提供
  params 后运行，未提供打印双语待办。
- **竞争投票不含 predict**：阶段 A 多数投票仅用
  by-ion/all-ion 两指标——任何报告结果都不如 pDeep2 原生生成谱的相似性高
  （理论注入谱即 pDeep2 产物、余弦理论 1.00），predict 会系统性偏向生成谱
  而失实；predict 仍计算并写入候选列供报告（predict_missing 照常标注）。
  旧数字（六轮与理论 2 各轮）在"predict 参与投票且带修饰候选键尾分号假
  缺失（4,758 行）"口径下取得，已记录；修正口径以 *_c 重跑为准。
- **结果文件夹定向 `RESULT_DIR` 环境变量**：master 级生效——
  同项目多套口径并行重跑时指定各自新文件夹，避免复用旧 merge 的假重跑；
  未设时按最新 result_* 定位（断点续跑）。
- **嵌合谱图错误谱图标注**：嵌合判定（chimeric /
  Chimera same charge）的新肽段谱图定性为**错误谱图**——Part2 报告
  （2.novel_peptide_competition_report.tsv）的 `spectrum_verdict` 列标
  "错误谱图（混合谱判定）"，行归属与可信判定不变。
- **denovo 回帖模式键 `denovo_lookback`**：
  `keep`（默认）＝仅保留能回帖到搜索库的 denovo 肽段（历史行为）；
  `all`＝放开全部——no_match 的 denovo 序列带 `e.DB_match=no_match` 标记
  保留进竞争（SAP 类库外变异肽段由此进入仲裁；终点残留的可解释错误案例
  ALAEFQLR 即属此类）。config_step1/2 均可设。
- **pFind 结果后处理**：自动搜库（开放式/限定式）完成
  或复用后，软件对两个 pFind 项目执行 `pfind_res_postprocess`：
  ①读取 result 目录全部 `*.qry.res`（`pFind_res_reader`，候选行含 rank/序列/
  分数/特异性/修饰（"位点,名称;"）/蛋白命中（`id;诱饵旗标;位置;前后残基`，
  旗标 1＝decoy）/is_decoy；`res_field_map` 给出固定列名表）；
  ②以**未过滤 `pFind.spectra`（全量 PSM）替换** Filtered 参与竞争匹配；
  ③其中 `Target/Decoy＝decoy` 的行用该谱 res 中 **Raw_Score 最高的非 decoy
  候选替换**（序列/修饰/特异性/分数重算；蛋白名与位置按目标库序列回贴解析，
  不依赖 res 蛋白 id；res 无非 decoy 候选的行保留原样交下游诱饵去除）；
  ④结果写回 `pFind-Filtered.spectra`（官方原版备份 `pFind-Filtered.origin.spectra`），
  供后续结果合并、理论谱图生成、竞争打分消费。幂等（origin 不重复备份）。

### 结果文件夹根与项目根

- `code_fingerprint.txt`：核心代码指纹（断点续跑时与当前代码比对，不一致
  会醒目告警——建议删除全部 result_* 重跑）。
- `pdf_on_flag.txt`：PDF_ON 状态记录（网页"谱图 PDF"卡片显示"已开启/默认关闭"）。
- `config_step2.txt / config_step3.txt`：第二步/第三步参数文件（软件自动生成，
  保留用户已填条目；predict_file 由 pDeep2 自动运行后回写）。

### 谱图 PDF（PDF_ON=1 时；默认关闭）

- 三类谱图 PDF 按 **raw 分组存储**（v1.0 起）：每 raw 一个文件
  `<类别>__<raw>.pdf`，组内页数超过 `pdf_pages_per_file`（config_step1 键，
  默认 1000）时再切 `__part2/__part3…`；页内版式与页眉不变。
- `2.labeled_result\error_match_spectra__*.pdf`：确信原始结果错误。
- `2.labeled_result\novel_peptide_spectra__*.pdf`：可信新肽段谱图——两分法
  报告 `credible_novel=TRUE` 的每条新肽段一页四面板（上=新肽段本身 corrected_*
  口径，下=最强已注释竞争者 closest_annotated_competitor）。
- `1.unlabeled_result\`：阶段 A 侧同名三类 PDF（error/ambiguous/novel，
  v1.0 起可信新肽段与标记无关）。
- **索引**：`spectrum_pdf_index.tsv`（各结果目录内）记录每谱
  File_Name → 类别 / PDF 文件 / 页码，按谱标题定位文件，无需整本加载。
- `1.unlabeled_result\`：阶段 A 侧 PDF（candidate misidentification 谱）。

## 多引擎版（isoPSV_multiengine）策略槽位

- 每种策略类型硬编码为固定单字母槽位（a=原始结果，b=pfind_open，c=pfind_restrict，
  d=punifind_rescore，e=punifind_denovo，f=maxquant，g=comet，h=msfragger，
  i=sage，j=msfragger_boost，k=xtandem，l=metamorpheus，m=peaks，n=mascot，o=pnovo，
  p=msfragger_open；msfragger_boost（MSBooster 重打分结果，格式与常规
  相同、复用同一适配器）在 j 槽；other 通用槽顺延 q..z（共 10 个），
  上限 26 个策略；合并按字母序自然排序。
- 合并表与 report_A 的槽列名直接为 `a.Sequence`、`b.Sequence` 等（字母开头即合法
  R 列名，无需 X 前缀）；Peptide_state/Modification_state 同样用单字母（如 "abd;c"）。
- config_step2 新增并行键 `search_file_type`（与 search_file 一一对应，取值即上表
  类型名；未填时按文件名/目录规则推断，旧配置零改动可用；固定类型重复报错，
  other 可重复使用依次占用 q..z）。
- 主表候选列随启用策略数扩展：cand1..cand5 与经典版 96 列完全一致，追加
  cand6..candN（每组 15 列）。
- 引擎适配器：`search_adapter_maxquant.R` 把 MaxQuant msms.txt 转为交换格式
  （修饰全名映射、raw+scan→谱图标题映射、Reverse/污染/PEP 过滤）；其他引擎照此
  模式添加适配器即可挂到对应槽位。

## 已确认的口径（与新数据集核对用）

1. 候选＝（肽段序列＋修饰）组合；同肽不同修饰/位点＝不同候选（修饰等价归一后合并）。
2. 修饰等价：`1,Carbamidomethyl[C]` ≡ `1,Carbamidomethyl[AnyN-term]`，N 端修饰位置 0≡1、AnyC-term 对应末位——按修饰类型逐一判断（AA_mass_dataframe.txt 位点形式清单），首位残基可修饰才替换。
3. 特异子离子比例仅用于混合谱判定（SI 1.8(3)：两肽各自特异子离子÷两肽总子离子，取最小；>10% 且特异离子数>1 判混合谱）。
4. Error_reason 按 SI 1.8(3)–(5) 逐级判断；SAP 严格命中（File_Name＋revise_peptide＝SAP 清单 revised_peptide）优先叠加报告。
5. ref_Protein 空值＝修正肽段不在参考蛋白组（保持空值，作为归因证据）。
6. **搜库引擎无关**：本流程的全部判定指标（无标记/标记离子匹配得分、覆盖度、
   预测余弦、成对余弦）均由软件从谱图独立计算，搜库结果只提供"谱图号＋肽段
   序列＋修饰"候选；因此任何搜库引擎（pFind、pUniFind、MaxQuant 等）的鉴定
   结果经统一交换格式转换后都可接入——作为原始结果（策略 1）或补充竞争策略
   （search_file 槽位，最少三列 File_Name/Sequence/Modification 即可），引擎
   打分不参与判定。
7. **酶切位点口径**：trypsin 含 **P 位阻断**（K/R 后遇 P 不切），等价于
   pFind 的 `Trypsin_P`；LysargiNase 同口径（K/R 后遇 P 不切，N 端酶）。
8. **第一步输入口径**：第一步 `psm_file` 接受**任何带
   `File_Name` 列的制表符结果表**——pFind 家族表（File_Name+Label）、适配器
   交换格式（`.spectra`）与其他搜库软件的表格输出均可用。非表格格式自动转换：
   X!Tandem XML（`<bioml>`）自动识别并经适配器转换（同目录 `<同名>_renum_map.tsv`
   存在时自动按 map 解析谱标题）；Mascot/PEAKS CSV 经 `psm_type` 键声明后自动
   转换；需重编号映射/参数文件的引擎（msfragger/pnovo/mzIdentML）给定向指引
   （走其表格输出或 config_step2 search_file 槽位）。无标记/标记谱区分：有
   `Label` 列按 pFind 语义（首字符 1＝无标记）；无 `Label` 列时无标记样本全部
   谱按无标记处理，标记样本须设 `psm_all_unlabeled=1` 显式确认。
9. **谱图 PDF 第三类**：可信新肽段谱图
   纳入 PDF，**与标记与否无关**——标记侧与无标记侧两分法报告都以同一套六重
   存疑规则（判错 > 嵌合 > 峰集全同 > 已注释有力竞争 > 连续缺失≥2 > 低分
   卡值）计算 `novel_doubt/novel_status/credible_novel`，全部通过者
   `credible_novel=TRUE`；无标记侧数据源＝阶段 A 逐谱候选（candidates 列）＋
   第一步无标记谱（τ 与主表同口径）。`novel_peptide_spectra__<raw>.pdf`
   （v1.0 起按 raw 拆分，超 pdf_pages_per_file 再切 part）落在
   `2.labeled_result\`（标记）与 `1.unlabeled_result\`（无标记），每页四面板：
   上=新肽段本身、下=最强已注释竞争者。
   **证据列口径**：`n_frag`（匹配子离子数）与
   `max_missing_run`（最大连续缺失位点数）仅当该谱在第一步无标记 MGF
   （`unlabeled\<raw>_HCDFT.mgf`）中定位成功且存在相应匹配证据时填写；
   谱图未能定位（标题未命中/谱缺失）或全部位点均有匹配时留空——谱图依赖的
   存疑规则（峰集全同/连续缺失）在该情形下跳过，其余规则照常判定。标记侧
   主表（master_table_export.R）同口径。

## 结果与归一化口径

软件端输出已统一以下口径（历史结果文件仍可直接被网页版读取，读取端自动归一）：

- **末端丢失归一化**：Met-loss 修饰转为 N 端 M 去除（末端类修饰位置 1→0、锚定被去
  残基的修饰丢弃、后续修饰位置前移）；Arg-loss 转为 C 端 R 去除；conflict 报告的
  X1–X5 序列同口径归一。
- **SAP 处理**：候选带 X→Y 替换记法时，pDeep2 预测键先做序列代入（如
  `DASAVILVSDAR` 的 `1,Asp->Asn[D]` 按 `NASAVILVSDAR` 预测）；Xle 歧义位生成 I/L
  双条目，匹配取最大余弦；SAP 清单仍由第二步自动提取。
- **状态二值小写**：`match_state`/`paired_state` 输出 `single`/`paired`
  （历史文件的 Paired/Matched/Double 等取值由消费端归一）。
- **逐谱 PDF 谱图呈现**（`谱图可视化.R`，v4）：双面板（原始/校正候选），
  a/b/c/x/y/z（1+/2+）＋M 系理论梯、20 ppm 窗匹配；颜色表 a #004000/b #008000/
  c #CA1D52/x #650E29/y #CA1D52/z #943AA4/M #808080/I #FF8C00/O #9ACD32；
  电荷号不省略（y10+/y10++）；同峰多离子合并单标签（b/y 为准、同现以 y 为准）；
  相邻峰（≤50 Da）标注窜开。

## 网页版（web）

双击 `start_web.bat` 打开（Python 标准库本地服务器，无第三方依赖）：

- 表单生成 config_step1 → 一键运行 user_run.R（实时日志、停点续接）；
- **结果浏览模式（v1.0）**：把另一台电脑跑完的整个项目文件夹
  （`unlabeled\` ＋ `result_*`）拷到本地，在结果报告面板选择该文件夹并"加载"——
  报告、谱图与 PDF 即直接从该目录只读展示，无需本地重跑；清除即切回本地当前结果；
- `/api/report` 双语结果报告：运行时间线＋阶段 A/B 统计＋输出文件清单；阶段 B 状态
  列自动兼容历史口径（paired_state/double_state 列名与 Paired/Matched/Double 取值
  均归一为 Single/Paired）；
- **谱图浏览卡片**按 File_Name 加载任意谱图（MGF 按标题索引
  定位），双面板候选对比（上＝候选 1，下＝下拉任选竞争/校正候选，默认生产胜选或
  候选 2），理论离子梯/匹配/颜色/标注规则与 PDF v4 完全一致（canvas 本地绘制，
  无外部 JS 依赖）；谱图 PDF 卡片把错误/存疑匹配（error/ambiguous matches）PDF 嵌入窗口内查看打印
  （自动扫描结果文件夹与 web\pdf\，也可输入任意文件夹路径）。

## 附录：可选引擎清单（用途与下载链接）／Appendix: optional engines, roles & download links

（v1.0 起"环境报告"页面并入本附录；本地检测与路径回填由第一步"One-click
detect & fill"完成，此处只列各软件在管线中的角色与官方下载地址。）
(Since v1.0 the "Engines" report page is merged into this appendix; live
detection and path filling are done by "One-click detect & fill" in STEP 1 -
this appendix only lists each engine's role and its official download page.)

| 软件 Software | 管线角色 Role | 下载 Download |
|---|---|---|
| pFind | 主搜库引擎（策略 a/b/c）primary search engine (strategies a/b/c) | http://www.pfind.net |
| pDeep2 | 理论谱图预测（predict 指标）theoretical spectrum prediction (predict metric) | https://github.com/pFindStudio/pDeep |
| pUniFind | rescore / de novo 补充策略 supplementary strategies | http://www.pfind.net |
| Python + torch venv | pUniFind 运行环境 pUniFind runtime | https://www.python.org |
| MSFragger | 搜库引擎策略 f search engine strategy f | https://msfragger.nesvilab.org |
| Comet | 搜库引擎策略 g search engine strategy g | https://github.com/UWPR/Comet |
| Sage | 搜库引擎策略 i search engine strategy i | https://github.com/lazear/sage |
| MaxQuant | 搜库引擎（策略 f 槽，自动搜库需 raw）search engine (auto search needs raw) | https://www.maxquant.org |
| MetaMorpheus | 搜库引擎策略 l search engine strategy l | https://github.com/smith-chem-wisc/MetaMorpheus |
| ProteoWizard msconvert | RAW→MGF 转换 RAW-to-MGF conversion | https://proteowizard.sourceforge.io |
| MSBooster | MSFragger 重打分（策略 j，可选）MSFragger rescoring (strategy j, optional) | https://github.com/Nesvilab/MSBooster |
| R runtime | 管线主体运行时（必需；≥4.1，推荐 4.5；软件自带可免装）pipeline runtime (required; bundled runtime may suffice) | https://cran.r-project.org |
| X!Tandem | 搜库引擎策略 k（可选，手动结果）search engine strategy k (optional, manual) | https://www.thegpm.org/TANDEM/ |
| Mascot | 搜库引擎策略 n（可选，手动结果）search engine strategy n (optional, manual) | https://www.matrixscience.com |

环境前提（实测归纳／prerequisites from live-fire notes）：

1. MetaMorpheus 需 .NET 10 运行时（`winget install Microsoft.DotNet.Runtime.10`），
   且 metamorpheus_exe 须指向官方**完整命令行包**内的 CMD.exe（含 Data/、Mods/；
   仅 exe/dll 的部分解压不可用）。MetaMorpheus needs the .NET 10 runtime and the
   full command-line package (Data/, Mods/); a partial exe/dll extract will not work.
2. 引擎工作区与项目目录必须 ASCII 路径；启动 R 前清除 LC_CTYPE / LANG / LC_ALL。
   Engine workspaces and project directories must use ASCII paths; clear
   LC_CTYPE / LANG / LC_ALL before starting R.
3. 未配置 msconvert 时，Sage / MSFragger / MetaMorpheus 走内置 MGF→mzML 转换器
   （每谱带 centroid spectrum，Sage 已实测出结果）。Without msconvert,
   Sage/MSFragger/MetaMorpheus use the built-in MGF→mzML converter (centroid
   spectrum per scan; verified with Sage).

## 附录：理论谱图预测器安装（pDeep2 / pDeep3）／Appendix: predictor installation

（口径：预测器软件与环境由用户自行从官方来源安装配置，本软件
不再分发其安装包。）
(Since the predictor software and environments are installed by
the user from the official sources; this software no longer redistributes
the packages.)

| 软件 | 官方来源 | 运行环境 |
|---|---|---|
| pDeep2 | https://github.com/pFindStudio/pDeep | Python 3.7 ＋ TensorFlow 1.15 |
| pDeep3 | https://github.com/pFindStudio/pDeep3（含模型） | Python 3.7 ＋ TensorFlow 1.15 |
| MS2PIP | `pip install ms2pip` | Python 3.8+ |
| AlphaPeptDeep | `pip install peptdeep` | Python 3.8+ |

pDeep2/pDeep3 环境配置（两者相同）：

1. 安装 Python 3.7（或复用 pFind Studio 安装自带的 py37 环境，位于
   `<pFind 安装目录>\bin` 同级或 pDeep2 目录下的 py37）；
2. `pip install tensorflow==1.15.0 scipy pandas matplotlib`；
3. `git clone` 上述仓库（pDeep3 仓库含全部模型权重）；
4. isoPSV 网页 STEP 3：pDeep2/pDeep3 行填**安装文件夹**，MS2PIP/
   AlphaPeptDeep 行填**装入该软件的 Python 环境文件夹**；参数留空＝全部
   默认（pDeep3/AlphaPeptDeep 取第一步质谱仪键，MS2PIP 用内置 HCD 模型）；
   可点"⚡ 一键填写路径与参数"从检测结果自动补路径。

pDeep2 vs pUniFind Python 不可互换：pDeep2/3 需 Python 3.7＋TF 1.15；
pUniFind 需任意 Python＋PyTorch。
