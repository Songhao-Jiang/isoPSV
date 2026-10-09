# build_user_guide.py - 生成《isoPSV 用户指南》完整中文 Word 版（图片占位框）
# 章节：1 软件介绍 / 2 软件原理 / 3 软件包含文件 / 4 操作流程（网页+命令行）/
# 5 结果文件格式与内容 / 6 列名解释表格（解析 结果文件列名对照.md）/
# 7 常见问题与配置键落点。图片占位框供用户自行插图。
import io, os, re
from docx import Document
from docx.shared import Pt, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = r"G:\JKY_work\Double_label\26.isoPSV正式流程\isoPSV_multiengine"
OUT = os.path.join(ROOT, "docs", "isoPSV_User_Guide_CN.docx")

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(10.5)
style._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
style.paragraph_format.line_spacing = 1.5
style.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

def run(p, text, bold=False, color=None, size=None, mono=False):
    r = p.add_run(text)
    r.bold = bold
    r.font.name = "Consolas" if mono else "Times New Roman"
    r._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    if color: r.font.color.rgb = color
    if size: r.font.size = Pt(size)
    return r

def h1(text):
    p = doc.add_paragraph()
    run(p, text, bold=True, color=RGBColor(0x0D, 0x3B, 0x66), size=15)
    return p

def h2(text):
    p = doc.add_paragraph()
    run(p, text, bold=True, color=RGBColor(0x14, 0x50, 0x7C), size=12)
    return p

def para(label, text):
    p = doc.add_paragraph()
    if label:
        run(p, label + "：", bold=True, color=RGBColor(0x14, 0x50, 0x7C))
    run(p, text)
    return p

def bullet(text, indent=0.5):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Cm(indent)
    run(p, "• " + text)
    return p

def codeblock(lines):
    for ln in lines:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(0.5)
        p.paragraph_format.line_spacing = 1.0
        run(p, ln, mono=True, size=9.5)
    return

def table(headers, rows, widths=None):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Table Grid"
    for i, htxt in enumerate(headers):
        c = t.rows[0].cells[i]
        c.paragraphs[0].paragraph_format.line_spacing = 1.2
        run(c.paragraphs[0], htxt, bold=True, size=9)
    for row in rows:
        cells = t.add_row().cells
        rv = list(row) + [""] * len(headers)
        for i in range(len(headers)):
            cells[i].paragraphs[0].paragraph_format.line_spacing = 1.2
            run(cells[i].paragraphs[0], str(rv[i]), size=9)
    return t

def placeholder(caption, height_cm=5.5):
    t = doc.add_table(rows=1, cols=1)
    t.style = "Table Grid"
    cell = t.rows[0].cells[0]
    cell.paragraphs[0].paragraph_format.line_spacing = 1.5
    cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
    run(cell.paragraphs[0], "【图片占位：%s】\n（在此处插入截图，插入后删除本行文字）" % caption,
        color=RGBColor(0x88, 0x88, 0x88), size=9)
    trPr = t.rows[0]._tr.get_or_add_trPr()
    trH = OxmlElement("w:trHeight")
    trH.set(qn("w:val"), str(int(height_cm * 567)))
    trH.set(qn("w:hRule"), "atLeast")
    trPr.append(trH)

# ================= 封面 =================
p = doc.add_paragraph()
run(p, "isoPSV 多引擎竞争验证软件 · 用户指南", bold=True, color=RGBColor(0x0D, 0x3B, 0x66), size=18)
para("", "版本 v1.0（2026-09-28）。全文的【图片占位】框请自行插入对应界面的截图，插入后删除占位文字。")

# ================= 1 软件介绍 =================
h1("1 软件介绍")
h2("1.1 定位与核心能力")
para("", "isoPSV 是一款纯 R 实现的肽段鉴定竞争验证软件：把多个搜库引擎与搜库策略对同一批谱图的鉴定结果汇集为竞争候选集，以谱图本身的可测指标、理论预测谱图相似性与稳定同位素代谢标记的成对谱图证据做竞争仲裁，逐谱判定并输出错误谱图校正与新肽段验证结果。")
bullet("多引擎竞争：pFind（常规/开放式/限定式）、pUniFind（重打分/从头测序）与 MaxQuant、Comet、MSFragger、Sage、X!Tandem、MetaMorpheus、PEAKS、Mascot、pNovo 等引擎结果经统一交换格式接入，按固定槽位字母合并竞争。")
bullet("引擎打分无关：全部判定指标（子离子打分、b/y 覆盖度、理论预测余弦、成对标记余弦）由软件从谱图独立计算，任何引擎的鉴定结果只提供“谱图＋肽段序列＋修饰”候选。")
bullet("新肽段验证：与样本是否标记无关的两分法报告（已注释校正 / 新肽段竞争），六重存疑规则判定可信新肽段，附逐谱 PDF 证据。")
bullet("纯蛋白质基因组学模式与全候选竞争模式（v1.0+）：只看新基因/新肽段鉴定（跳过已注释错误校正），或忽略引擎 1% FDR 过滤、每引擎每谱前 5 候选进竞争池（小样本可用）。")
h2("1.2 名称来源")
para("", "isoPSV 由三部分组成：iso 取自 isotope（同位素），代表软件所依赖的稳定同位素代谢标记（15N/13C/SILAC）成对谱图这一核心物理证据；P 取自 Peptide（肽段）；SV 取自 Spectral Validation（谱图验证），即以单张谱图为单位的逐谱验证口径。")
h2("1.3 应用场景")
bullet("常规蛋白质组：多引擎相互矛盾的鉴定结果仲裁，输出错误谱图校正清单与证据。")
bullet("蛋白质基因组学：新肽段/新基因鉴定证据（纯蛋白质基因组学模式）。")
bullet("稳定同位素标记定量研究：成对标记谱图物理证据复核（阶段 B）。")
h2("1.4 运行环境")
bullet("操作系统：Windows 10 / 11；R ≥ 4.1（推荐 4.5.x，仅基础包，无第三方 R 包依赖）；Python 3.8+（重编号与网页界面；未写入 PATH 时在 config_step1 填 python_exe）。")
bullet("绿色运行时：包内自带 runtime\\python（官方 Python 3.12.10 embeddable）与 runtime\\R（R 4.5.1 基础+推荐包裁剪版）——未安装 Python/R 的电脑双击根目录 isoPSV.exe 即可使用。")
bullet("外部工具（独立软件，按需安装）：pFind 3.1、pUniFind（GPU 推荐）、pDeep2、Comet、MSFragger、Sage、MetaMorpheus、MaxQuant 等；未安装的引擎可提供其人工产出结果文件。")

# ================= 2 软件原理 =================
h1("2 软件原理")
h2("2.1 三步流程")
table(["步骤", "内容", "输出"],
      [["第一步", "无标记谱图提取（按 raw 重切），多引擎结果表输入（任何带 File_Name 列的制表符表；XML/CSV 自动转换）", "unlabeled\\ 逐 raw 纯无标记 MGF + step1_done.txt + 自动生成 config_step2.txt"],
       ["第二步", "各搜库策略执行（自动调 pFind/pUniFind/外部引擎或读用户产物）→ 去诱饵 → 修订（I/L 归一、N=GG、四类修饰校正）→ 合并为逐谱一行候选表 → pDeep2 理论预测", "result_*\\Competitive_search_result\\ 与 merge\\（1.combine_result.spectra、peptide.txt、predict.txt）+ 自动生成 config_step3.txt"],
       ["第三步", "阶段 A 无标记逐谱仲裁 + 阶段 B 成对标记复核 + 主表 + 两分法报告 + 逐谱 PDF（可开关）", "1.unlabeled_result\\ 与 2.labeled_result\\ 全套报告"]])
h2("2.2 竞争仲裁机制")
bullet("阶段 A（无标记侧）：对每谱全部候选计算 PSM 子离子打分（by-ion 与 all-ion）、b/y 覆盖度与 pDeep2 预测余弦，可用指标多数投票并加权定胜负；候选两两特异子离子比例最小值 >10% 且特异离子数 >1 判混合谱。")
bullet("阶段 B（标记侧）：无标记谱与 15N/13C/SILAC 标记谱按保留时间配对（容差 20 秒），按标记态重算理论质量后计算成对谱图余弦——独立于搜库打分的物理证据；成对余弦 >0.7 判成对。")
bullet("预测救援：分数差接近判定阈值时按 0.3 阈值找回校正结果（pred_rescue）。")
bullet("新肽段验证：判错 > 嵌合 > 峰集全同 > 已注释有力竞争 > 连续缺失≥2 > 低分卡值六重存疑规则，全部通过判可信新肽段。")
h2("2.3 成对稳定同位素标记与金标准")
para("", "同一细胞/菌株分别在普通与重同位素培养基（15N/13C 代谢标记或 SILAC）中培养，轻重样品混合后一同酶切、色谱分离与质谱采集：同一肽段的轻重谱图成对出现，两谱质量差由肽段原子组成编码唯一确定。对三组标定样本数万张谱逐谱复核标注判对/判错/存疑真值构造金标准数据集，用于质控模型评估（约登指数）与判定权重、阈值的拟合，贯穿开发与质控。")
h2("2.4 判定口径（锁定值）")
bullet("候选＝（肽段序列＋修饰）组合；同肽不同修饰/位点＝不同候选。")
bullet("阶段 A 权重 ul_by/ul_all/ul_cov/predict＝0.1306/0.4418/0.2381/0.1895，阈值 threshold_A＝0.0204；阶段 B 六指标权重（ul_by,lb_by,ul_all,lb_all,predict,pcos），阈值 threshold_B＝0.0217。")
bullet("判定阈值 τ：n15/c13/unlabeled＝0.0083，silac＝0.0183；pred_rescue θ＝0.30；成对余弦阈值 0.7；混合谱特异子离子比例 >10% 且特异离子数 >1。")
bullet("定侧口径：Proteins 任一蛋白命中用户 FASTA 即“已注释”，否则为新肽候选侧（与两分法报告一致）。")
bullet("Trypsin 与 LysargiNase 均为 K/R 后遇脯氨酸不切；跨版本可比性注意 Release_Notes 各版本条目。")

# ================= 3 软件包含文件 =================
h1("3 软件包含文件")
h2("3.1 目录结构")
codeblock([
 "isoPSV_multiengine/",
 "├─ isoPSV.exe / start_isoPSV.bat   一键启动器（双击开网页界面；优先内置运行时）",
 "├─ runtime/                        绿色运行时：python（官方 3.12.10 embeddable）、R（4.5.1 基础+推荐包）",
 "├─ code/                           全部 R 代码（仅基础包）",
 "│   ├─ user_run.R                  总入口（三步调度、断点续跑）",
 "│   ├─ unlabeled_extract.R / user_step2.R / user_step3.R   三步驱动",
 "│   ├─ stageA_classify.R / master_table_export.R / extend_report_A.R  阶段 A 分类/主表/报告扩充",
 "│   ├─ ReadMGF_script.R / compat.R / paired_params.R / acs_core.R     核心库",
 "│   ├─ search_adapter_*.R          各引擎结果→交换格式适配器",
 "│   ├─ predictors\\pred_*.R         理论谱图预测器适配器",
 "│   ├─ engine_templates\\ / pFind_templates\\    引擎与 pFind 参数模板",
 "│   ├─ AA_mass_dataframe.txt       氨基酸质量与修饰表",
 "│   └─ renumber_mgf.py             重编号脚本（需 Python）",
 "├─ web/                            本地网页界面（server.py / index.html / spectrum.js）",
 "├─ tests/                          自动化回归测试（run_all.R 30 套 + web_smoke.py 25 项）",
 "├─ docs/                           随包手册 PDF 与选项详解 Word",
 "├─ config_step1_template.txt       第一步参数模板",
 "├─ VERSION / Release_Notes.md / Software_Manual.md / README.md / 结果文件列名对照.md",
 "└─ 项目目录（运行时生成）：unlabeled\\、result_时间戳\\"])
h2("3.2 各组成部分说明")
table(["文件/目录", "说明"],
      [["isoPSV.exe / start_isoPSV.bat", "一键启动本地网页界面（127.0.0.1:8780）；优先使用 runtime\\ 绿色运行时，其次回退 PATH；关闭窗口即停止服务"],
       ["runtime\\", "绿色 Python 与 R 运行时（无安装电脑可直接运行）"],
       ["code\\", "全部 R 代码与资源；用户不直接修改"],
       ["web\\", "网页壳子；server.py 仅监听本机端口"],
       ["tests\\", "自检套件；run_all.bat 一键运行"],
       ["config_step1_template.txt", "第一步参数模板（复制改名为 config_step1.txt 后填写）"]])

# ================= 4 操作流程 =================
h1("4 操作流程")
h2("4.1 网页流程")
para("", "双击包根目录 isoPSV.exe（或 start_isoPSV.bat）启动本地服务，浏览器自动打开 http://127.0.0.1:8780。左侧导航按顺序完成以下各页。")
para("4.1.1 参数填写（STEP 1）", "填写项目文件夹（result_dir）、原始搜库结果（psm_file）、参考蛋白库（fasta）、标记类型、酶切类型、MGF 文件；可选“一键检测并填写路径”自动扫描本机引擎安装位置回填；其他补充选项（PDF 开关、回帖模式、质谱仪类型、候选摄入口径、工作流程）按需设置；点击“生成参数文件”写入 config_step1.txt。")
placeholder("图 1：参数填写界面（STEP 1）")
para("4.1.2 多引擎策略（STEP 2）", "按策略槽位添加外部引擎行（类型下拉仅列外部引擎；pFind 家族由第一步 auto_search / use_punifind 管理），填写结果文件/文件夹与引擎类型；可一键检测填充自动搜库引擎的软件路径；完成后点击“生成 config_step2”。")
placeholder("图 2：多引擎策略配置界面（STEP 2）")
para("4.1.3 理论谱图预测（STEP 3）", "预测器行配置各预测软件（pDeep2 行为安装路径+自动/手动/关闭；其他行为 Python 路径+模型目录+启用/关闭），选择主预测器（可选“不进行理论谱图预测”），点击“写入 config_step1”。")
placeholder("图 3：理论谱图预测器配置界面（STEP 3）")
para("4.1.4 运行与日志", "点击“运行”启动三步流程；状态条显示当前阶段，时间线列出各模块墙钟与 CPU 时间；步骤停点时按提示补填后再次点击运行（断点续跑）。")
placeholder("图 4：运行与日志界面")
para("4.1.5 结果报告", "完成后在结果报告页选择结果目录，浏览阶段 A/阶段 B 报告、两分法报告与逐谱 PDF 卡片；谱图查看器按 PDF 卡片相同标注规则在浏览器内画谱。结果浏览模式可只读加载另一台电脑的项目/结果文件夹。")
placeholder("图 5：结果报告浏览界面")
para("4.1.6 环境报告", "列出本机已安装的搜库与依赖软件及其目录，未安装的给出官网下载链接。")
placeholder("图 6：环境报告界面")
h2("4.2 命令行流程")
para("", "用户全程只与总程序 user_run.R 交互：准备好参数文件后经 CONFIG 环境变量告知其绝对路径并运行；总程序自动判断已完成步骤，需要补参数时打印双语提示并退出——填好后重跑同一条命令即续接。")
codeblock([
 "cd /d \"G:\\isoPSV_multiengine\"",
 "set CONFIG=G:\\项目文件夹\\config_step1.txt",
 "\"%RSCRIPT%\" code\\user_run.R        （RSCRIPT 指向 Rscript.exe；中文乱码时先 chcp 65001）"])
para("三轮运行", "第一轮：第一步提取无标记谱、第二步执行搜库策略，停点提示填 config_step2.txt 的 search_file（手动搜库时）；第二轮：合并候选表并自动生成 config_step3.txt，停点提示填 predict_file（手动 pDeep2 时）；第三轮：阶段 A/B 仲裁与全部报告输出。标记样本另需在 config_step1.txt 填 all_mgf（原始 MGF 文件夹）。")
table(["步骤", "参数文件", "输出", "跳过条件"],
      [["第一步", "config_step1.txt（用户填）", "unlabeled\\ 逐 raw MGF + step1_done.txt；自动生成 config_step2.txt", "step1_done.txt 已存在"],
       ["第二步", "config_step2.txt（填 search_file）", "Competitive_search_result\\ 与 merge\\（合并表 + peptide.txt）；自动生成 config_step3.txt", "merge\\peptide.txt 已存在"],
       ["第三步", "config_step3.txt（填 predict_file）", "1.unlabeled_result\\ 阶段 A 全套报告", "每次重算"],
       ["阶段 B（仅标记）", "同第三步", "2.labeled_result\\ 成对/冲突报告与主表", "paired_report 已存在"]])
para("运行约定", "环境变量传路径（CONFIG/CONFIG2/CONFIG3），中文路径一律走环境变量、不进命令行参数；核数与并行数自动评估；全部输出为制表符分隔文本（UTF-8）；引擎工作区与项目目录须 ASCII 路径；任何一步中断后重跑同一命令即续接，参数/数据变化按指纹强制重算相关步骤。")

# ================= 5 结果文件格式与内容 =================
h1("5 结果文件格式与内容")
h2("5.1 结果目录结构")
codeblock([
 "result_YYYYMMDDHHMMSS\\",
 "├─ Competitive_search_result\\   各策略去诱饵清洗件 + SAP 索引",
 "├─ merge\\                       合并候选表、peptide.txt、predict.txt、SAP 索引、修订件",
 "├─ 1.unlabeled_result\\          阶段 A：report_A.tsv、classification_A.tsv、两分法报告、逐谱 PDF",
 "└─ 2.labeled_result\\（标记样本）阶段 B：成对/冲突报告、master_table.tsv、两分法报告、逐谱 PDF"])
h2("5.2 各文件内容")
table(["文件", "内容"],
      [["unlabeled\\<raw>_HCDFT.mgf", "逐 raw 纯无标记谱图（第一步产物，第二步起只读）"],
       ["Competitive_search_result\\*.spectra", "各策略结果去诱饵清洗件（Proteins 全诱饵行删除、混合行只留正库）"],
       ["merge\\1.combine_result.spectra", "多策略合并候选表：逐谱一行，槽位列 a.Sequence/a.Modification…，含 Peptide_state/Modification_state"],
       ["merge\\peptide.txt", "pDeep2 理论预测输入（谱图、肽段、修饰、电荷）"],
       ["merge\\predict.txt", "pDeep2 预测产物"],
       ["1.unlabeled_result\\report_A.tsv", "阶段 A 逐谱报告（全部谱图，含逐候选指标与投票）"],
       ["1.unlabeled_result\\classification_A.tsv", "阶段 A 逐谱判定结论（含 error_reason 与修饰差异细节）"],
       ["1.unlabeled_result\\summary_A.txt", "阶段 A 统计"],
       ["1.unlabeled_result\\1.annotation_correction_report.tsv / 2.novel_peptide_competition_report.tsv", "两分法报告（已注释校正 / 新肽段竞争，无标记样本同构输出）"],
       ["2.labeled_result\\paired_report_<标记类型>.tsv", "成对标记匹配报告"],
       ["2.labeled_result\\conflict_report_<标记类型>.tsv + _summary.txt", "冲突谱全参数与分组统计"],
       ["2.labeled_result\\master_table.tsv", "全谱主表（判校正行全列填写）"],
       ["两类 PDF 目录", "判错/存疑/新肽段三类逐谱 PDF（pdf_on=1 时），按 raw 拆分 + spectrum_pdf_index.tsv 索引"]])

# ================= 6 列名解释表格 =================
h1("6 列名解释表格")
md = io.open(os.path.join(ROOT, "结果文件列名对照.md"), encoding="utf-8").read()
cur_rows = None
def flush_rows():
    global cur_rows
    if cur_rows:
        table(cur_rows[0], cur_rows[1:])
        cur_rows = None
for ln in md.splitlines():
    if ln.startswith("### "):
        flush_rows()
        h2("6." + ln[4:].strip())
    elif ln.startswith("## "):
        flush_rows()
        h2("6." + ln[3:].strip())
    elif ln.startswith("|"):
        cells = [c.strip() for c in ln.strip("|").split("|")]
        if cells and all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
            continue
        if cur_rows is None:
            cur_rows = []
        cur_rows.append(cells)
    elif ln.strip() and not ln.startswith(">"):
        flush_rows()
        para("", ln)
flush_rows()

# ================= 7 常见问题与配置键落点 =================
h1("7 常见问题与配置键落点")
h2("7.1 停点与断点续跑")
bullet("任何步骤中断后重跑同一条命令即从断点续接；已完成步骤按完成标记/指纹自动跳过。")
bullet("自动搜库 pUniFind 项产出未齐时第二步停在待办提示（不生成 peptide.txt）；手动模式按拷贝指引把 GUI 产出拷到指定目录后直接重跑。")
bullet("删除结果目录 result_* 可整体强制重算；参数/数据变化按指纹只重算相关步骤。")
h2("7.2 引擎与运行排查")
bullet("引擎失败日志区分“引擎非零退出”与“产物未发现”两种情形，分别提示排查方向。")
bullet("引擎工作区与项目目录须 ASCII 路径；FASTA 路径含空格/中文时 Comet/MSFragger 无法读取——请用 ASCII 路径副本。")
bullet("外部引擎自动运行需 config_step1 的 auto_search=1 且 pfind_exe 非空；step1 会按模板重生成 config_step2（保留 search_file 行、丢弃引擎键），引擎键须在 step1 完成后写入 config_step2 再重跑。")
bullet("纯蛋白质基因组学模式（proteogenomics_mode=1/2）：搜库前预检输入 psm_file 的 Proteins 列，无目标侧行时直接结束（退出码 0、写 proteogenomics_summary.txt）；被跳过的谱在 classification_A.tsv 中补记为 annotated_skipped（纯蛋白质基因组学）或 novel_skipped（已注释矫正），指标列留空。")
h2("7.3 配置键落点速查")
table(["键", "读取位置（生效来源）", "说明"],
      [["proteogenomics_mode", "第一步（config_step1）", "第三步经 cfg3 的 base_config（直指第一步）读取；搜库前早退预检同源；v1.0 起 config_step3 重建时保留该键"],
       ["candidate_scope", "第一步（config_step1）", "第二步与第一步转换路径均读"],
       ["punifind_exe / punifind_python / punifind_rescore / punifind_denovo", "第一步优先（step1→step2 取第一个非空合法值）", "两文件同键不同值会打印双语告警（v1.0），按口径取 step1"],
       ["use_pdeep2 / pdeep2_exe / pdeep2_instrument", "第一步（网页表单写入第一步）", "第三步按 cfg3 复制的值执行"],
       ["<engine>_exe / <engine>_mode / <engine>_workers", "第二步（config_step2）", "外部引擎自动运行与逐 raw 并行"],
       ["pfind_workers", "第一步（config_step1）", "pFind 自动搜库逐 raw 并行"],
       ["python_exe", "第一步（config_step1）", "Python 解释器显式路径（未入 PATH 时）"]])

doc.save(OUT)
print("written:", OUT)
