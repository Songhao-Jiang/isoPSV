# build_options_guide.py - 生成《选项详解》手册 HTML（中/英）→ 交给 LibreOffice 转 PDF
import io, os

OUT = r"G:\JKY_work\Double_label\26.isoPSV正式流程\isoPSV_multiengine\docs\manual_figures"
DOC = r"G:\JKY_work\Double_label\26.isoPSV正式流程\isoPSV_multiengine\docs"
WORK = r"C:\Users\MSI-NB\AppData\Local\Temp\optguide"
if not os.path.isdir(WORK): os.makedirs(WORK)

# key -> (中文名, 英文名, 默认值, 负责什么[cn,en], 内容[(值, cn, en)], 影响[cn,en])
ITEMS = [
 ("search_engine", "搜库引擎", "Search engine", "pFind（默认）",
  "声明第一步 psm_file（常规搜库结果）来自哪个引擎。",
  "Declares which engine produced the step-1 psm_file.",
  [("pFind", "常规搜库结果表（含 File_Name、Label 列）", "regular result table (File_Name + Label columns)"),
   ("MaxQuant", "msms.txt，经内置适配器转交换格式", "msms.txt, converted by the built-in adapter"),
   ("MSFragger", "tsv / pepXML，经内置适配器转换", "tsv / pepXML, converted by the built-in adapter"),
   ("Sage", "results.sage.tsv，经内置适配器转换", "results.sage.tsv, converted by the built-in adapter"),
   ("Comet", "pep.xml / txt，经内置适配器转换", "pep.xml / txt, converted by the built-in adapter")],
  "只影响第一步输入口径与 config_step1 该键取值；多引擎竞争在第二步策略槽独立配置，与本选项无关。",
  "Affects only the step-1 input declaration; the multi-engine competition is configured independently in the step-2 strategy slots."),
 ("auto_search", "pFind 搜库方式", "pFind search mode", "自动 Auto（默认）",
  "负责第二步 pFind 开放式 / 限定式搜库的执行方式。",
  "How the step-2 pFind open/restricted search is performed.",
  [("自动 Auto", "软件自动调 pFind 完成搜库（需填 pfind_exe）", "the software runs pFind (needs pfind_exe)"),
   ("手动 Manual", "用户自行跑 pFind，产物填 config_step2 的 search_file", "user runs pFind; results go into config_step2 search_file")],
  "自动模式真正执行搜库（耗时，open 约 14 分钟级 / restrict 约 149 分钟级于 7 组分示例）；手动模式只读用户产物，搜索参数由用户完全自控。",
  "Auto truly runs the search (open ~14 min / restrict ~149 min on the 7-fraction example); manual only reads user outputs, giving the user full control of search parameters."),
 ("use_punifind", "启用 pUniFind", "Enable pUniFind", "启用 Enable（默认）",
  "pUniFind 重打分与从头测序（策略 d / e）的总开关。",
  "Master switch for pUniFind re-scoring and de novo (strategies d/e).",
  [("启用 Enable", "第二步执行 pUniFind 流程", "pUniFind workflow runs in step 2"),
   ("禁用 Disable", "纯 pFind 流程，d/e 槽不产生", "pure-pFind workflow; slots d/e produce nothing")],
  "禁用时 rescore / denovo 全部跳过，合并表不含 d / e 槽，相关停点也不再出现。",
  "Disabled = rescore/denovo are skipped entirely; the merged table has no d/e slots and their stop points disappear."),
 ("punifind_rescore", "重打分运行方式", "Open-search re-score mode", "不运行 Off（默认）",
  "策略 d（开放式搜库结果重打分校正）的运行方式。",
  "How strategy d (open-search re-scoring) runs.",
  [("不运行 Off", "策略 d 不产生", "strategy d produces nothing"),
   ("自动 Auto", "软件调官方 CLI 自动完成（需 punifind_exe 与 punifind_python）", "the software runs the official CLI (needs punifind_exe and punifind_python)"),
   ("手动 Manual", "用户 GUI 完成并拷入结果，软件打印指引并停点", "user runs the GUI and copies the result; the software prints instructions and stops")],
  "决定策略 d 的来源与第二步是否调 GPU CLI（自动约 10 分钟）；手动模式在结果未到位时停点等待。",
  "Determines where strategy d comes from and whether step 2 calls the GPU CLI (~10 min); manual mode stops until the result is in place."),
 ("punifind_denovo", "从头测序运行方式", "De novo sequencing mode", "不运行 Off（默认）",
  "策略 e（从头测序）的运行方式。",
  "How strategy e (de novo) runs.",
  [("不运行 Off", "策略 e 不产生", "strategy e produces nothing"),
   ("自动 Auto", "软件调官方 CLI 自动完成（GPU，约 1.5 小时量级）", "the software runs the official CLI (GPU, ~1.5 h)"),
   ("手动 Manual", "用户 GUI 完成并拷入 pUniFind_result，软件打印指引并停点", "user runs the GUI and copies pUniFind_result; the software prints instructions and stops")],
  "决定策略 e 的来源；从头序列一律经数据库回帖（denovo_lookback）后才参与竞争。",
  "Determines where strategy e comes from; de novo sequences always pass the database look-back (denovo_lookback) before competing."),
 ("label_type", "标记类型", "Labeling type", "无标记 Unlabeled（默认）",
  "样本的代谢标记类型。",
  "The metabolic labeling type of the sample.",
  [("无标记 Unlabeled", "只跑阶段 A", "stage A only"),
   ("15N 全标记", "氮-15 代谢标记", "nitrogen-15 metabolic labeling"),
   ("13C 全标记", "碳-13 代谢标记", "carbon-13 metabolic labeling"),
   ("SILAC（通用）", "细胞培养氨基酸标记通用档", "generic SILAC"),
   ("SILAC NeuCode K602 / K080 / R004 / R040 / R604 / R640", "NeuCode 双标记精细档", "fine NeuCode dual-labeling tiers")],
  "决定阶段 B 成对匹配时标记谱的理论质量重算（TAU_MAP 与阈值）与标记侧报告口径；无标记样本只跑阶段 A、无阶段 B。",
  "Determines the labeled-spectrum mass recalculation (TAU_MAP and thresholds) and the labeled-side reports in stage B; unlabeled samples run stage A only."),
 ("enzyme", "酶切类型", "Cleavage enzyme", "Trypsin（默认，P 切不开）",
  "样品的酶切方式。",
  "The cleavage enzyme of the sample.",
  [("Trypsin", "K/R 后切，遇脯氨酸不切（默认）", "cleaves after K/R, blocked by P (default)"),
   ("LysargiNase", "K/R 前切（N 端酶），遇脯氨酸不切", "cleaves before K/R (N-terminal), blocked by P"),
   ("LysC / LysN", "K 后切 / K 前切", "after K / before K"),
   ("ArgC / GluC / AspN", "R 后切 / E、D 后切 / D 前切", "after R / after E,D / before D"),
   ("Chymotrypsin", "F/Y/W/L 后切（糜蛋白酶）", "after F/Y/W/L"),
   ("无酶切 Unspecific", "不做酶切限制", "no cleavage restriction")],
  "影响两处：自动搜库时改写各引擎（pFind / MSFragger / Sage / Comet / MaxQuant / MetaMorpheus）模板的酶切参数；管线内部漏切计数与五类分类的 N 端非特异判据。",
  "Two effects: the enzyme parameters rewritten into each engine template during auto search, and the in-pipeline missed-cleavage count plus the N-terminal non-specific rule of the five-class classification."),
 ("pdf_on", "谱图 PDF 导出", "Spectrum PDF export", "不生成 Off（默认）",
  "是否导出逐谱谱图 PDF（评审材料）。",
  "Whether to export per-spectrum PDFs (review material).",
  [("不生成 Off", "不导出（默认）", "no export (default)"),
   ("生成 On", "导出逐谱 PDF", "export per-spectrum PDFs")],
  "开启后输出判错 / 存疑 / 新肽段三类 PDF，按 raw 拆分（pdf_pages_per_file 再切 part）并带 spectrum_pdf_index.tsv 索引；生成与查看两端内存有界。",
  "On = three classes (error/ambiguous/novel) exported per raw (sliced by pdf_pages_per_file) with the spectrum_pdf_index.tsv index; bounded memory on both ends."),
 ("denovo_lookback", "从头测序回帖", "De novo lookback mode", "keep（默认）",
  "从头测序序列的回帖口径。",
  "The look-back mode for de novo sequences.",
  [("keep", "仅保留能回帖到搜索数据库的序列（默认）", "only database-matching sequences survive (default)"),
   ("all", "放开全部从头序列参与竞争", "all de novo sequences compete")],
  "keep 噪声低、候选精；all 让库外序列（如 SAP 型变异）带标记进入仲裁，候选更多、计算更大。",
  "keep = low noise, precise candidates; all = out-of-database sequences (e.g. SAP-type variants) enter arbitration, with more candidates and more compute."),
 ("pdeep2_instrument", "质谱仪类型", "Mass spectrometer type", "Lumos（默认）",
  "理论谱图预测所用的仪器参数。",
  "The instrument parameter for theoretical-spectrum prediction.",
  [("Lumos", "默认仪器", "default instrument"),
   ("Q Exactive", "QE 系列", "QE series"),
   ("Velos / Elite", "老机型", "older instruments")],
  "不同仪器的碎裂特征影响预测谱形态，进而影响阶段 A 的 predict 余弦一票（预测仅作辅助证据）。",
  "Instrument fragmentation characteristics shape the predicted spectra and thus the stage-A predict-cosine vote (prediction is auxiliary evidence only)."),
 ("candidate_scope", "候选摄入口径", "Candidate scope", "依赖搜库软件 1% FDR（默认）",
  "竞争候选池的摄入口径。",
  "How the competition candidate pool is ingested.",
  [("依赖搜库软件 1% FDR", "各引擎自身过滤后的结果集（默认，现状）", "each engine's own filtered result set (default)"),
   ("全候选竞争模式（每谱前 5 候选）", "忽略引擎 FDR，每引擎每谱取前 5 个肽段+修饰组合", "engine FDR ignored; top-5 peptide+modification candidates per spectrum per engine")],
  "全候选模式下候选池约膨胀 1.4～2.2 倍、更多谱进入竞争；判定指标仍由软件从谱图独立计算（判定口径不变），小样本不依赖引擎 FDR 也可用；Comet / MSFragger 会自动输出每谱前 5（v1.0），Sage 仍每谱 1 候选（引擎输出所限）。",
  "Full-candidate mode inflates the pool ~1.4-2.2x and brings more spectra into competition; metrics are still recomputed from spectra (judgment rules unchanged) so small samples work without engine FDR. Comet/MSFragger auto-emit top-5 per spectrum (v1.0); Sage still contributes 1 candidate per spectrum (engine output limit)."),
 ("proteogenomics_mode", "工作流程", "Workflow scope", "全量模式（默认）",
  "哪些谱进入竞争验证与错误校正的工作流程门控。",
  "The workflow gate deciding which spectra enter validation and error correction.",
  [("全量模式", "所有谱竞争验证与错误校正（默认）", "all spectra validated and corrected (default)"),
   ("纯蛋白质基因组学模式", "只验证含非参考库候选的谱（新肽段/新基因鉴定），已注释谱跳过", "only spectra with non-reference candidates (novel peptide/gene); annotated spectra skipped"),
   ("已注释肽段矫正模式", "只矫正已注释谱，新肽候选谱跳过", "only annotated spectra corrected; novel-candidate spectra skipped")],
  "定侧口径与两分法报告一致（Proteins 任一命中用户 FASTA 即已注释）；被跳过的谱标记 annotated_skipped / novel_skipped 补记进 classification_A；搜库前预检无目标侧行时直接结束（退出码 0）。",
  "Side assignment matches the two-part report (any Proteins hit in the user FASTA = annotated); skipped spectra are appended to classification_A as annotated_skipped / novel_skipped; a pre-search check exits early (code 0) when no target-side rows exist."),
 ("engine_type", "引擎类型（策略行）", "Engine type (strategy row)", "Comet（示例）",
  "第二步策略行内所接引擎的类型。",
  "The engine type of a step-2 strategy row.",
  [("Comet / MaxQuant / MSFragger / Sage / MetaMorpheus", "具备自动搜库能力，行内出现软件路径输入框", "auto-search capable; a software-path input appears in the row"),
   ("X!Tandem / Mascot / PEAKS / pNovo", "GUI 导出类，结果文件/文件夹填行内结果框", "GUI exports; the result file/folder goes into the row's result box"),
   ("MSFragger boost / MSFragger open", "MSBooster 重打分与开放搜索变体", "MSBooster re-scoring and open-search variants")],
  "决定槽位字母（f=MaxQuant、g=Comet、h=MSFragger、i=Sage…）、所用适配器与 revise 模式；pFind 家族（b–e）不在此下拉，由第一步 auto_search / use_punifind 管理。",
  "Determines the slot letter (f=MaxQuant, g=Comet, h=MSFragger, i=Sage, ...), the adapter and the revise mode; the pFind family (b-e) is managed by step-1 auto_search / use_punifind, not this dropdown."),
 ("predictor_primary", "主预测器", "Primary predictor", "pDeep2（默认）",
  "阶段 A「预测余弦」一票的来源。",
  "The source of the stage-A predict-cosine vote.",
  [("pDeep2", "默认主预测器（独立软件，需安装并填路径）", "default primary (separate software; path required)"),
   ("MS2PIP / Prosit / 自定义脚本", "外部预测器，按行添加后可改选", "external predictors, selectable after rows are added"),
   ("不进行理论谱图预测", "关闭预测：predict 列按缺失处理", "prediction off: predict column treated as missing")],
  "选「不进行理论谱图预测」时自动把 use_pdeep2 置 0（predict 列按缺失处理、投票降级）；改回预测器且 use_pdeep2 仍为 0 时自动恢复为 1（自动模式）。",
  "Selecting 'no theoretical prediction' sets use_pdeep2 to 0 (predict column missing, votes degrade); selecting a predictor again restores use_pdeep2 to 1 (auto mode) when it is still 0."),
 ("predictor_row", "预测器行（类型+开关）", "Predictor row (type + switch)", "pDeep2 · 自动 Auto",
  "每个预测器行的类型与开关。",
  "The type and switch of each predictor row.",
  [("pDeep2 行", "安装路径 + 自动/手动/关闭（use_pdeep2）", "install path + auto/manual/off (use_pdeep2)"),
   ("MS2PIP / Prosit / 自定义行", "Python 路径 + 模型目录 + 启用/关闭", "Python path + model dir + enable/disable"),
   ("行规则", "每行一个软件、类型可随时改换、至多一行 pDeep2", "one software per row, type changeable, at most one pDeep2 row")],
  "各行开关决定对应预测软件是否运行；产物 merge/predict_<type>.txt 与 pDeep2 同构，运行失败只记 manifest、该列按缺失处理；删除 pDeep2 行后「添加预测器」先补回。",
  "Each row's switch decides whether that predictor runs; outputs (merge/predict_<type>.txt) are pDeep2-shaped, failures only mark the manifest; deleting the pDeep2 row makes 'Add predictor' restore it first."),
]

def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def build(lang):
    cn = lang == "cn"
    title = "isoPSV 网页版选项详解" if cn else "isoPSV Web Options Guide"
    sub = "每个可选项：负责什么 · 选项内容 · 具体影响（v1.0）" if cn else \
          "Every option: what it does, what it contains, and its impact (v1.0)"
    note = ("小图为界面示意图（按实际界面的键名、选项与样式绘制）。" if cn else
            "The small figures are UI sketches drawn from the real key names, options and styling.")
    parts = ['<!DOCTYPE html><html><head><meta charset="utf-8"><title>%s</title><style>' % title,
             '@page { size: A4; margin: 2cm 1.9cm; }',
             'body { font-family: "%s", sans-serif; font-size: 10.5pt; color:#1a1a1a; line-height:1.5; }' % ("Microsoft YaHei" if cn else "Segoe UI"),
             'h1 { font-size: 18pt; color:#0d3b66; margin: 0 0 2pt 0; }',
             '.sub { color:#444; margin-bottom: 10pt; }',
             'h2 { font-size: 13pt; color:#0d3b66; border-bottom: 1.5pt solid #0d3b66; padding-bottom:2pt; margin:14pt 0 6pt 0; page-break-after: avoid; }',
             'p { margin: 3pt 0; }',
             'img { max-width: 92%; border: 0.6pt solid #c9d4de; border-radius: 6pt; margin: 4pt 0; page-break-inside: avoid; }',
             '.kv { font-weight:600; color:#14507c; }',
             'ul { margin: 3pt 0; padding-left: 16pt; }',
             'li { margin-bottom: 1pt; }',
             '</style></head><body>',
             '<h1>%s</h1><div class="sub">%s</div><p class="sub">%s</p>' % (title, sub, note)]
    for i, (key, nm_cn, nm_en, default_v, what_cn, what_en, opts, imp_cn, imp_en) in enumerate(ITEMS, 1):
        nm = nm_cn if cn else nm_en
        what = what_cn if cn else what_en
        imp = imp_cn if cn else imp_en
        parts.append('<h2>%d. %s · %s</h2>' % (i, key, nm))
        parts.append('<p><span class="kv">%s</span>：%s</p>' % ("负责什么" if cn else "What it does", esc(what)))
        parts.append('<p><span class="kv">%s</span>（当前值：%s）</p>' % ("选项内容" if cn else "Options", esc(default_v)))
        parts.append('<ul>')
        for v, d_cn, d_en in opts:
            parts.append('<li><span class="kv">%s</span> — %s</li>' % (esc(v), esc(d_cn if cn else d_en)))
        parts.append('</ul>')
        parts.append('<p><span class="kv">%s</span>：%s</p>' % ("具体影响" if cn else "Impact", esc(imp)))
        parts.append('<img src="figs/%s_%s.png" style="width:86%%">' % (key, "cn" if cn else "en"))
    parts.append('</body></html>')
    out = os.path.join(WORK, "_options_guide_%s.html" % lang)
    with io.open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(parts))
    print("html:", out)

build("cn")
build("en")
