# gen_options_figs.py - 为说明手册生成每个可选项的"下拉示意小图"（界面示意图）
# 环境无显示通道（IAB 截图 surface 超时），按真实界面值与样式用 PIL 确定性绘制：
# 每图 = 卡片（键名+中文名标签）+ 下拉框（当前值+▾）+ 展开的选项清单。
# 每项生成 中/英 两版，供 CN/EN 手册嵌入。
from PIL import Image, ImageDraw, ImageFont
import os, json

OUT = r"G:\JKY_work\Double_label\26.isoPSV正式流程\isoPSV_multiengine\docs\manual_figures"
os.makedirs(OUT, exist_ok=True)
FONT = "C:/Windows/Fonts/msyh.ttc"
CARD_W = 760
PAD = 16
LBL_H = 34
SEL_H = 40
OPT_H = 34
OPT_IND = 24
COL_BORDER = (186, 196, 206)
COL_CARD = (250, 252, 254)
COL_SEL = (232, 238, 244)
COL_OPT = (255, 255, 255)
COL_TXT = (26, 26, 26)
COL_SUB = (90, 100, 110)
COL_ACC = (13, 59, 102)

def draw_fig(path, key, name, default_v, options, desc_lines):
    n = len(options)
    H = PAD * 2 + LBL_H + SEL_H + OPT_H * n + (12 if desc_lines else 0) + sum(OPT_H for _ in (desc_lines or []))
    img = Image.new("RGB", (CARD_W, H), COL_CARD)
    d = ImageDraw.Draw(img)
    d.rounded_rectangle([0, 0, CARD_W - 1, H - 1], radius=10, fill=COL_CARD, outline=COL_BORDER, width=2)
    f_lbl = ImageFont.truetype(FONT, 22)
    f_sel = ImageFont.truetype(FONT, 20)
    f_opt = ImageFont.truetype(FONT, 18)
    f_sub = ImageFont.truetype(FONT, 16)
    y = PAD
    d.text((PAD, y), key + "  ·  " + name, font=f_lbl, fill=COL_ACC)
    y += LBL_H
    d.rounded_rectangle([PAD, y, CARD_W - PAD, y + SEL_H], radius=8, fill=COL_SEL, outline=COL_BORDER, width=1)
    d.text((PAD + 12, y + 7), default_v, font=f_sel, fill=COL_TXT)
    tw = d.textlength("\u25be", font=f_sel)
    d.text((CARD_W - PAD - 26, y + 5), "\u25be", font=f_sel, fill=COL_SUB)
    y += SEL_H
    for val, txt in options:
        d.rounded_rectangle([PAD, y, CARD_W - PAD, y + OPT_H - 1], radius=4, fill=COL_OPT, outline=(226, 232, 238), width=1)
        d.text((PAD + OPT_IND, y + 5), val, font=f_opt, fill=COL_TXT)
        if txt:
            d.text((PAD + OPT_IND + d.textlength(val, font=f_opt) + 18, y + 7), txt, font=f_sub, fill=COL_SUB)
        y += OPT_H
    for ln in (desc_lines or []):
        d.text((PAD + 4, y + 4), ln, font=f_sub, fill=COL_SUB)
        y += OPT_H
    img.save(path)
    return img.size

# 每项定义：key, 中文名, 英文名, 默认值显示, [(值, 中文说明), ...], 图注行
ITEMS = [
 ("search_engine", "搜库引擎", "Search engine", "pFind（默认）",
  [("pFind", "常规搜库结果表（默认，含 File_Name 与 Label 列）"),
   ("MaxQuant", "msms.txt 结果，经内置适配器转交换格式"),
   ("MSFragger", "tsv/pepXML 结果，经内置适配器转换"),
   ("Sage", "results.sage.tsv 结果，经内置适配器转换"),
   ("Comet", "pep.xml/txt 结果，经内置适配器转换")], None),
 ("auto_search", "pFind 搜库方式", "pFind search mode", "自动 Auto（默认）",
  [("自动 Auto", "软件自动调 pFind 完成开放式与限定式搜库（需填 pfind_exe）"),
   ("手动 Manual", "用户自行跑 pFind，结果填入 config_step2 的 search_file")], None),
 ("use_punifind", "启用 pUniFind", "Enable pUniFind", "启用 Enable（默认）",
  [("启用 Enable", "第二步执行 pUniFind 重打分与从头测序（策略 d/e）"),
   ("禁用 Disable", "纯 pFind 流程，d/e 槽不参与")], None),
 ("punifind_rescore", "重打分运行方式", "Open-search re-score mode", "不运行 Off（默认）",
  [("不运行 Off", "策略 d 不产生"),
   ("自动 Auto", "软件调 pUniFind 官方 CLI 自动完成（需 punifind_exe 与 punifind_python）"),
   ("手动 Manual", "用户 GUI 完成并拷入结果，软件打印拷贝指引并停点")], None),
 ("punifind_denovo", "从头测序运行方式", "De novo sequencing mode", "不运行 Off（默认）",
  [("不运行 Off", "策略 e 不产生"),
   ("自动 Auto", "软件调 pUniFind 官方 CLI 自动完成（GPU，约 1.5 小时量级）"),
   ("手动 Manual", "用户 GUI 完成并拷入 pUniFind_result，软件打印指引并停点")], None),
 ("label_type", "标记类型", "Labeling type", "无标记 Unlabeled（默认）",
  [("无标记 Unlabeled", "只跑阶段 A 竞争验证"),
   ("15N 全标记", "氮-15 代谢标记，阶段 B 成对匹配用 15N 质量重算"),
   ("13C 全标记", "碳-13 代谢标记，同上"),
   ("SILAC（通用）", "细胞培养氨基酸标记，通用档"),
   ("SILAC NeuCode K602/K080/R004/R040/R604/R640", "NeuCode 双标记精细档")], None),
 ("enzyme", "酶切类型", "Cleavage enzyme", "Trypsin（默认，P 切不开）",
  [("Trypsin", "K/R 后切，遇脯氨酸不切（默认）"),
   ("LysargiNase", "K/R 前切（N 端酶），遇脯氨酸不切"),
   ("LysC / LysN", "K 后切 / K 前切"),
   ("ArgC / GluC / AspN", "R 后切 / E、D 后切 / D 前切"),
   ("Chymotrypsin", "F/Y/W/L 后切（糜蛋白酶）"),
   ("无酶切 Unspecific", "不做酶切限制")], None),
 ("pdf_on", "谱图 PDF 导出", "Spectrum PDF export", "不生成 Off（默认）",
  [("不生成 Off", "不导出逐谱 PDF（默认）"),
   ("生成 On", "导出判错/存疑/新肽段三类逐谱 PDF（按 raw 拆分+索引）")], None),
 ("denovo_lookback", "从头测序回帖", "De novo lookback mode", "keep（默认）",
  [("keep", "仅保留能回帖到搜索数据库的从头序列（默认）"),
   ("all", "放开全部从头序列参与竞争（库外序列带标记进入仲裁）")], None),
 ("pdeep2_instrument", "质谱仪类型", "Mass spectrometer type", "Lumos（默认）",
  [("Lumos", "默认仪器"),
   ("Q Exactive", "QE 系列"),
   ("Velos / Elite", "老机型")], None),
 ("candidate_scope", "候选摄入口径", "Candidate scope", "依赖搜库软件 1% FDR（默认）",
  [("依赖搜库软件 1% FDR", "各引擎自身过滤后的结果集进竞争（默认，现状）"),
   ("全候选竞争模式（每谱前 5 候选）", "忽略引擎 FDR 过滤，每个引擎每谱取前 5 个肽段+修饰组合进竞争池（小样本可用）")], None),
 ("proteogenomics_mode", "工作流程", "Workflow scope", "全量模式（默认）",
  [("全量模式", "所有谱竞争验证与错误校正（默认）"),
   ("纯蛋白质基因组学模式", "只验证含非参考库候选的谱（新肽段/新基因鉴定），已注释谱跳过错误校正"),
   ("已注释肽段矫正模式", "只矫正已注释谱，新肽候选谱跳过")], None),
 ("engine_type", "引擎类型（策略行）", "Engine type (strategy row)", "Comet（示例）",
  [("Comet / MaxQuant / MSFragger / Sage / MetaMorpheus", "具备自动搜库能力，行内可填软件路径"),
   ("X!Tandem / Mascot / PEAKS / pNovo", "GUI 导出类，结果文件/文件夹填左侧框"),
   ("MSFragger boost / MSFragger open", "MSFragger 的 boost（MSBooster）与开放搜索变体")], None),
 ("predictor_primary", "主预测器", "Primary predictor", "pDeep2（默认）",
  [("pDeep2", "内置主预测器（独立软件，需用户安装并填路径）"),
   ("MS2PIP / Prosit / 自定义脚本", "外部预测器（按行添加后可选）"),
   ("不进行理论谱图预测", "关闭预测：predict 列按缺失处理，use_pdeep2 自动置 0")], None),
 ("predictor_row", "预测器行（类型+开关）", "Predictor row (type + switch)", "pDeep2 · 自动 Auto",
  [("pDeep2 行", "安装路径 + 自动/手动/关闭"),
   ("MS2PIP/Prosit/自定义行", "Python 路径 + 模型目录 + 启用/关闭"),
   ("每行一个软件，可随时改换类型；至多一行 pDeep2")], None),
]

# 每项的中英图注（负责什么/内容/影响）
NOTES = {
 "search_engine": ("声明第一步 psm_file 的常规搜库结果来自哪个引擎（仅影响第一步输入口径说明与配置文件该键取值；多引擎竞争在第二步策略槽独立配置）。",
                   "Declares which engine produced the step-1 psm_file (affects the step-1 input declaration only; multi-engine competition is configured independently in the step-2 strategy slots)."),
 "auto_search": ("负责第二步 pFind 开放式/限定式搜库的执行方式。自动＝软件调 pFind 真跑搜库（耗时、需 pfind_exe）；手动＝只读用户提供的搜库产物（填 config_step2 search_file）。",
                 "How step-2 pFind open/restricted search runs. Auto = the software runs pFind (slow, needs pfind_exe); manual = user-provided search results are read (search_file in config_step2)."),
 "use_punifind": ("总开关：决定 pUniFind 重打分与从头测序（策略 d/e）是否参与第二步。禁用时这两个策略槽不产生。",
                  "Master switch for pUniFind re-scoring and de novo (strategies d/e) in step 2. Disabled = these slots produce nothing."),
 "punifind_rescore": ("策略 d（开放式搜库结果重打分校正）的运行方式。自动需 punifind_exe 与 punifind_python（GPU 约十分钟）；手动由用户 GUI 完成并拷入结果。",
                      "How strategy d (open-search re-scoring) runs. Auto needs punifind_exe and punifind_python (GPU, ~10 min); manual = user runs the GUI and copies the result."),
 "punifind_denovo": ("策略 e（从头测序）的运行方式。自动走 GPU 官方 CLI（约 1.5 小时）；手动由用户拷入 pUniFind_result。从头序列经数据库回帖后才参与竞争。",
                     "How strategy e (de novo) runs. Auto uses the GPU CLI (~1.5 h); manual = user copies pUniFind_result. De novo sequences enter competition only after database look-back."),
 "label_type": ("样本的代谢标记类型。决定阶段 B 成对匹配时标记谱的理论质量重算（TAU_MAP）与标记侧报告；无标记样本只跑阶段 A。",
                "Metabolic labeling type of the sample. Determines the labeled-spectrum mass recalculation (TAU_MAP) and the labeled-side reports in stage B; unlabeled samples run stage A only."),
 "enzyme": ("酶切类型。影响：① 自动搜库时改写各引擎（pFind/MSFragger/Sage/Comet/MaxQuant/MetaMorpheus）的酶切参数；② 管线内部漏切计数与五类分类的 N 端非特异判据。",
            "Cleavage enzyme. Affects (1) the enzyme parameters rewritten into each engine template during auto search, and (2) the in-pipeline missed-cleavage count and the N-terminal non-specific rule of the five-class classification."),
 "pdf_on": ("是否导出逐谱谱图 PDF（评审材料）。开启后按 raw 拆分输出判错/存疑/新肽段三类 PDF 与 spectrum_pdf_index.tsv 索引，生成与查看两端内存有界。",
            "Whether to export per-spectrum PDFs (review material). On = three classes (error/ambiguous/novel) split per raw plus the spectrum_pdf_index.tsv index; bounded memory on both generation and viewing."),
 "denovo_lookback": ("从头测序序列的回帖口径。keep＝仅保留能回帖到搜索数据库的序列（默认，噪声低）；all＝放开全部序列参与竞争（库外变异可进仲裁，候选更多）。",
                     "De novo look-back mode. keep = only sequences matching the search database survive (default, low noise); all = every de novo sequence competes (out-of-database variants enter arbitration, more candidates)."),
 "pdeep2_instrument": ("pDeep2 理论谱图预测的仪器参数。不同仪器的碎裂特征影响预测谱形态，进而影响阶段 A 的 predict 余弦一票。",
                       "Instrument parameter for pDeep2 theoretical-spectrum prediction. Fragmentation characteristics differ per instrument, shaping the predicted spectra and thus the stage-A predict-cosine vote."),
 "candidate_scope": ("候选摄入口径。依赖 FDR＝各引擎自身过滤后的结果集（默认）；全候选＝忽略引擎 1% FDR 过滤，每个引擎每谱取前 5 个肽段+修饰组合进竞争池——判定指标仍由软件从谱图独立计算，小样本不依赖引擎 FDR 也可用。",
                     "Candidate scope. FDR-reliant = each engine's own filtered result set (default); full-candidate = engine FDR filtering ignored and the top-5 peptide+modification candidates per spectrum per engine enter the pool - metrics are still recomputed from spectra, so small samples work without engine FDR."),
 "proteogenomics_mode": ("工作流程门控。全量＝所有谱竞争验证与错误校正；纯蛋白质基因组学＝只验证含非参考库候选的谱（新肽段/新基因鉴定），已注释谱标记 annotated_skipped 跳过；已注释矫正＝只矫正已注释谱。",
                         "Workflow gate. Full = competition validation and error correction on all spectra; pure proteogenomics = only spectra with non-reference candidates are validated (novel peptide/gene identification), annotated spectra are marked annotated_skipped; annotated-correction = only annotated spectra are corrected."),
 "engine_type": ("每个策略行的引擎类型。决定：槽位字母（g=Comet/h=MSFragger/i=Sage 等）、所用适配器与 revise 模式、是否具备自动搜库能力（行内出现软件路径输入框）。",
                 "Engine type of each strategy row. Determines the slot letter (g=Comet, h=MSFragger, i=Sage, ...), the adapter and revise mode used, and whether the row is auto-search capable (a software-path input appears in the row)."),
 "predictor_primary": ("阶段 A「预测余弦」一票的来源。pDeep2 为默认；外部预测器按行添加后可改选；选「不进行理论谱图预测」即关闭预测（predict 列按缺失处理、use_pdeep2 自动置 0）。",
                       "Source of the stage-A predict-cosine vote. pDeep2 is the default; external predictors become selectable after rows are added; 'no theoretical prediction' turns prediction off (predict column treated as missing, use_pdeep2 set to 0)."),
 "predictor_row": ("每个预测器行的类型与开关。pDeep2 行＝安装路径+自动/手动/关闭（use_pdeep2）；其他行＝Python 路径+模型目录+启用/关闭。每行一个软件、类型可随时改换、至多一行 pDeep2。",
                    "Type and switch of each predictor row. The pDeep2 row = install path + auto/manual/off (use_pdeep2); other rows = Python path + model dir + enable/disable. One software per row, type changeable anytime, at most one pDeep2 row."),
}

def option_list(key, lang):
    # 从 ITEMS 取选项，按语言渲染说明（简化：值固定、说明用双语短文本）
    return None

sizes = {}
for key, cn_name, en_name, default_v, opts, _ in ITEMS:
    vals = [(o[0], "") if isinstance(o, tuple) else (o, "") for o in opts]
    draw_fig(os.path.join(OUT, key + "_cn.png"), key, cn_name, default_v, vals, None)
    draw_fig(os.path.join(OUT, key + "_en.png"), key, en_name, default_v, vals, None)
print("figures:", len(ITEMS) * 2, "->", OUT)
