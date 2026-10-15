# =============================================================================
# user_run.R - master driver: the ONLY program the user calls.
#
# The user tells the program the path of the (user-edited) step-1 parameter
# file through the CONFIG environment variable (all paths inside that file
# must be absolute). The program walks the three-step pipeline and skips
# already-done work:
#   step 1  skipped when the completion marker exists
#           (unlabeled/step1_done.txt); otherwise run
#           unlabeled_extract.R, which auto-generates config_step2.txt
#   step 2  skipped when merge/peptide.txt (the pDeep2 input) exists;
#           otherwise run user_step2.R (needs the search_file entries the
#           user filled into config_step2.txt), which auto-generates
#           config_step3.txt and the SAP index
#   step 3  always run (user_step3.R), then extend_report_A.R
#   stage B after step 3 for labeled samples (label_type N15/C13/SILAC):
#   stageB_paired_match.R needs the all_mgf key (original MGF folder with
#   unlabeled + labeled spectra) from config_step1.txt; skipped for unlabeled
#   samples, skipped when stageB/paired_report_<label>.tsv already exists
#
# The SAP index is extracted automatically by step 2 from the supplementary
# search results; core counts are assessed automatically at every step (no
# user options for either).
#
# Whenever a parameter file still needs the user's input (search_file /
# predict_file), the program stops with a bilingual message; the user fills
# the file and calls this program again.
#
# Usage (single-sample workflow: no group concept):
#   CONFIG=config_step1.txt Rscript user_run.R
# =============================================================================
#
# user_run.R —— 总程序主控：用户唯一需要调用的程序。
# 用户通过 CONFIG 环境变量把（用户编辑过的）第一步参数文件路径告诉程序
# （该文件内的所有路径必须写为绝对路径）。程序依次执行三步流程并跳过已完成的工作：
#   第一步 当完成标记（unlabeled/step1_done.txt）已存在时跳过；
#          否则运行 unlabeled_extract.R，它会自动生成 config_step2.txt
#   第二步 当 merge/peptide.txt（pDeep2 输入）已存在时跳过；
#          否则运行 user_step2.R（需要用户填进 config_step2.txt 的 search_file 条目），
#          它会自动生成 config_step3.txt 和 SAP 索引
#   第三步 始终运行（user_step3.R），随后运行 extend_report_A.R
#   阶段 B 对标记样本在第三步之后执行：
#          stageB_paired_match.R 需要 config_step1.txt 里的 all_mgf 键
#          （含无标记＋标记谱的原始 MGF 文件夹）；无标记样本跳过，
#          当 stageB/paired_report_<label>.tsv 已存在时也跳过
#
# SAP 索引由第二步从补充搜库结果中自动提取；每一步的核数均自动评估
# （这两项都不需要用户选项）。
#
# 每当参数文件仍缺用户输入（search_file / predict_file）时，程序以双语消息停下；
# 用户填好文件后再次调用本程序即可。
#
# 用法（单样本流程：没有分组概念）：
#   CONFIG=config_step1.txt Rscript user_run.R

# ---- 先按用户启动目录解析 CONFIG（必须在工作目录切到脚本所在目录之前）。
# ---- 参数文件不存在（或未设置 CONFIG）时，在启动目录生成双语模板 config_step1.txt 并提示用户填写后重跑。
# Resolve CONFIG relative to the user's launch directory BEFORE switching the working directory
# to the script's folder. When the parameter file is missing (or CONFIG is unset), a bilingual
# template config_step1.txt is generated in the launch directory and the program exits telling
# the user to fill it and rerun.
launch_dir <- getwd()
cfg_arg <- Sys.getenv("CONFIG")
cfg_path <- if (nzchar(cfg_arg)) cfg_arg else file.path(launch_dir, "config_step1.txt")
if (!file.exists(cfg_path)) {
  # 仅用于提示信息；R 4.5 + MSYS/Git Bash（LC_ALL=C.UTF-8）下 normalizePath
  # 会腐蚀路径中的中文，此处只做分隔符清理
  cfg_path <- gsub("\\\\", "/", cfg_path)
  # 组装双语参数文件模板的每一行（键值制表符 + # 注释，pFind.cfg 风格）
  # Assemble the bilingual parameter-file template lines (key<TAB>value + # comments, pFind.cfg style)
  template_lines <- c(
    "# 用户第一步参数文件模板",
    "# user step-1 parameter file template",
    "# 键值制表符＋# 注释（pFind.cfg 风格）",
    "# key<TAB>value with # comments (pFind.cfg style)",
    "# 用户填写的所有路径一律写为绝对路径",
    "# all user-written paths must be ABSOLUTE",
    "# psm_file: pFind 常规搜库结果（含 File_Name、Label 列，Label 首字符 1＝无标记）",
    "#           pFind regular search result (File_Name + Label columns; Label first char 1 = unlabeled)",
    "psm_file\t<绝对路径 absolute path>",
    "# fasta: 参考蛋白库",
    "#        reference proteome FASTA",
    "fasta\t<绝对路径 absolute path>",
    "# label_type: C13 / N15 / SILAC / unlabeled",
    "label_type\t<标记类型 label type>",
    "# mgf: 每行一个；文件路径直接收，文件夹路径展开为其中全部 .mgf",
    "#      one entry per line; a file path is taken as-is, a folder expands to all .mgf inside",
    "mgf\t<绝对路径 absolute path>",
    "# all_mgf: 原始 MGF 文件夹（含无标记＋标记谱；阶段 B 成对匹配用；无标记样本留空）",
    "#          original MGF folder (unlabeled + labeled spectra; used by stage B paired matching; leave empty for unlabeled samples)",
    "all_mgf\t<绝对路径 absolute path，无标记样本可留空 may be empty for unlabeled samples>",
    "# auto_search: 第二步搜库方式（1=软件自动调 pFind 搜库；0=用户手动搜库）",
    "#              step-2 search mode (1=automatic pFind search; 0=manual, fill search_file)",
    "auto_search\t1",
    "# pfind_exe: pFind 安装路径（以 bin 结尾；填版本目录如 pFind3.1 时软件自动补 bin）",
    "#            pFind install path (must end with bin; a version folder like pFind3.1 gets bin appended automatically)",
    "pfind_exe\t<绝对路径 absolute path>",
    "# use_punifind: 第二步是否使用 pUniFind 两项（1=用 pUniFind rescore/denovo；0=不用，全 pFind 流程）",
    "#               use pUniFind rescore/denovo in step 2 (1=yes; 0=no, pure-pFind workflow)",
    "use_punifind\t1",
    "# use_pdeep2: 是否使用 pDeep2 理论预测谱（1=用，按 pdeep2_exe 安装目录自动运行；2=用，用户手动运行 pDeep2 生成 predict.txt 后填 predict_file；0=不用，预测指标按缺省处理）",
    "#             use pDeep2 predicted spectra (1=yes, auto-run via the pdeep2_exe install folder; 2=yes, the user runs pDeep2 manually and fills predict_file; 0=no, the predict metric is treated as missing)",
    "use_pdeep2\t1",
    "# pdeep2_exe: pDeep2 安装路径（use_pdeep2=1 时必填，软件自动执行 pDeep2 生成 predict.txt，第三步不再停点；use_pdeep2=2 时留空）",
    "#             pDeep2 install path (required with use_pdeep2=1: the software runs pDeep2 automatically and step 3 does not stop; leave empty with use_pdeep2=2)",
    "pdeep2_exe\t<绝对路径 absolute path>",
    "# instrument: 仪器类型（QE / Lumos / Velos / Elite，默认 Lumos）",
    "#                   instrument type (QE / Lumos / Velos / Elite; default Lumos)",
    "instrument\tLumos"
  )
  writeLines(template_lines, cfg_path, useBytes = TRUE)
  # 写出模板文件并用双语消息提示用户填写后重跑
  # Write the template file and prompt the user (bilingually) to fill it and rerun
  cat("==== [总程序 / Master] 参数文件不存在，已生成模板：", cfg_path, " ====\n")
  cat("==== [Master] parameter file missing, template generated: ", cfg_path, " ====\n")
  cat("请填写：", cfg_path,
      "（路径写绝对路径），然后运行：",
      paste0("CONFIG=", if (nzchar(cfg_arg)) cfg_arg else "config_step1.txt"),
      "，并重跑本程序。\n")
  cat("fill in: ", cfg_path,
      " (absolute paths), then run: ",
      paste0("CONFIG=", if (nzchar(cfg_arg)) cfg_arg else "config_step1.txt"),
      ", and rerun this program.\n")
  quit(status = 0)
}
cfg_path0 <- cfg_path
cfg_path <- normalizePath(cfg_path, winslash = "/")
# R 4.5 + MSYS/Git Bash（LC_ALL=C.UTF-8）下 normalizePath 会腐蚀路径中的中文
# （宽字符经 ANSI 代码页往返损坏）。结果不可用时退化为纯字符串解析：绝对路径
# 直接清理分隔符，相对路径按启动目录拼接；文件存在性已在上方校验过。
if (!file.exists(cfg_path)) {
  cfg_path <- if (grepl("^[A-Za-z]:", cfg_path0)) cfg_path0 else file.path(launch_dir, cfg_path0)
  cfg_path <- gsub("\\\\", "/", cfg_path)
}

# 所有运行都切到本脚本所在目录（库文件按相对路径 source，必须先于 source）。
# Rscript 前端以 --file=<路径> 传脚本名（另有 --no-echo），故按模式定位而非位置。
# Run everything from this script's own folder (the library sources are relative to it, so this
# must happen BEFORE source()). The Rscript front-end passes the script as --file=<path> (plus
# --no-echo), so locate it by pattern, not position.
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
if (is.na(script_arg)) script_arg <- args0[grepl("\\.R$", args0)][1]
if (is.na(script_arg)) {
  stop("无法定位脚本路径：请用 Rscript 运行（Rscript user_run.R）\ncannot locate the script path: run via Rscript (Rscript user_run.R)")
}
script_arg <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_arg, winslash = "/"))
# 同上：normalizePath 腐蚀中文时按字符串解析（绝对路径直接取目录，相对路径按当前目录拼接）
if (!dir.exists(script_dir)) {
  script_dir <- dirname(if (grepl("^[A-Za-z]:", script_arg)) script_arg else file.path(getwd(), script_arg))
  script_dir <- gsub("\\\\", "/", script_dir)
}
# setwd 容错控制台源文件编码受限时上述路径仍可能腐蚀致
# setwd 失败；若当前目录已含 compat.R（即从代码目录启动）则直接继续
# tolerant setwd under a restricted console encoding the
# path above can still be corrupted and setwd fails; when the working directory
# already holds compat.R (i.e. launched from the code folder) just continue
ok_wd <- tryCatch({ setwd(script_dir); TRUE }, error = function(e) FALSE)
if (!ok_wd) {
  if (file.exists(file.path(getwd(), "compat.R"))) {
    cat("note: already in the code directory, continuing without setwd\n")
  } else {
    stop("cannot change into the code directory: ", script_dir,
         "\n请从软件代码目录启动，或检查路径编码\nrun from the code folder or check the path encoding")
  }
}
# 定位 Rscript 可执行文件（先试 x64 子目录，再退回 bin 根目录）
# Locate the Rscript executable (try the x64 subfolder first, then the plain bin folder)
rscript <- file.path(R.home("bin"), "x64", "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript.exe")

# 载入库文件（兼容层 + MGF 读取模块），随后解析第一步参数文件
# Load the library files (compatibility layer + MGF reader), then parse the step-1 parameter file
source("compat.R")
source("ReadMGF_script.R")

cfg1 <- read_param_file(cfg_path)
cfg_dir <- cfg1$config_dir

# 用户/用户反馈①b理论谱图预测器独立参数文件 config_predictor.txt
# 多路径查找：①config_step1 同目录 ②config_step1 的 result_dir 键 ③包根
# ⚠ v1.0 修复：不能用 R 变量 result_dir（第 290 行才定义），改用 cfg1 键
cfg.pred.path <- ""
.rd.from.cfg <- if (!is.null(cfg1$result_dir) && length(cfg1$result_dir) > 0)
  tail(cfg1$result_dir, 1) else ""
for (.cand in c(
  file.path(dirname(cfg_path), "config_predictor.txt"),
  if (nzchar(.rd.from.cfg)) file.path(.rd.from.cfg, "config_predictor.txt") else "",
  file.path(getwd(), "..", "config_predictor.txt"))) {
  if (nzchar(.cand) && file.exists(.cand)) {
    cfg.pred.path <- .cand
    break
  }
}
if (nzchar(cfg.pred.path)) {
  cfg.pred <- read_param_file(cfg.pred.path)
  if (!is.null(cfg.pred) && length(cfg.pred) > 0) {
    for (k.pred in names(cfg.pred)) {
      if (grepl("^predictor", k.pred))
        cfg1[[k.pred]] <- cfg.pred[[k.pred]]
    }
    cat("[config_predictor] 已读取并合并 predictor 键：", cfg.pred.path, "\n")
    cat("[config_predictor] merged predictor keys from: ", cfg.pred.path, "\n")
  }
} else {
  cat("[config_predictor] 未找到（查找了 config 目录、result_dir、包根）——跳过\n")
  cat("[config_predictor] not found (config dir, result_dir, package root) - skipped\n")
}

# PDF 导出选项config_step1 键 pdf_on（1=生成、0=不生成，
# 默认 0）优先，环境变量 PDF_ON 兜底（两者取"1"即开启）
# PDF export option the config_step1 key pdf_on
# (1=on, 0=off, default 0) wins, the PDF_ON environment variable is the
# fallback (either one "1" enables the export)
# P0-1 修复旧行同时要求 pdf_on != "1" 与 == "1"，恒假死代码，
# PDF_ON 环境变量从未生效；改为 trimws 读取
pdf_on <- "0"
pv <- cfg1$pdf_on
if (!is.null(pv) && length(pv) > 0 && pv[1] %in% c("0", "1")) pdf_on <- pv[1]
if (identical(trimws(Sys.getenv("PDF_ON", unset = "0")), "1")) pdf_on <- "1"

# 内存守卫上限config_step1 键 ISOMEM_LIMIT_MB（MB 数，
# 网页高级区可填）注入环境变量，iso_mem_guard 由此读取；留空/非法值时
# 沿用默认（物理内存 80% 自适应，上限 40000）。用户手填超过物理内存的值
# iso_mem_guard 侧只警告不拦截（显式覆盖，后果自负）。
# Memory guard ceiling the config_step1 key ISOMEM_LIMIT_MB
# (MB, fillable in the web advanced panel) is injected into the environment
# for iso_mem_guard; empty/invalid values keep the default (80% of physical
# RAM, capped at 40000). A user value above physical RAM only warns inside
# iso_mem_guard (explicit override, at the user's own risk).
mg_mb <- cfg1$ISOMEM_LIMIT_MB
if (!is.null(mg_mb) && length(mg_mb) > 0 && nzchar(trimws(mg_mb[1]))) {
  mg_n <- suppressWarnings(as.numeric(trimws(mg_mb[1])))
  if (!is.na(mg_n) && mg_n > 0) Sys.setenv(ISOMEM_LIMIT_MB = format(mg_n, scientific = FALSE))
}

# 酶切类型config_step1 键 enzyme（trypsin/trypsin_p/
# lysarginase/lysc/lysn/argc/gluc/aspn/chymotrypsin/none，默认 trypsin）。
# 注入 ISO_ENZYME 环境变量，子进程（阶段 A/B、主表导出、pFind 后处理）经
# iso_enzyme() 读取，用于漏切计数与 N 端非特异分类；第二步搜库由 user_step2
# 读 base_cfg 改写各引擎模板。未知值回退 trypsin 并警告。
# Enzyme type the config_step1 key enzyme (trypsin/
# trypsin_p/lysarginase/lysc/lysn/argc/gluc/aspn/chymotrypsin/none, default
# trypsin). Injected as ISO_ENZYME for the child processes (stages A/B, master
# table export, pFind postprocess) to read via iso_enzyme() for the missed-
# cleavage count and the N-terminal non-specific classification; the step-2
# search rewrites each engine template from base_cfg. Unknown values fall back
# to trypsin with a warning.
en_cfg <- cfg1$enzyme
en_src <- "默认"
if (!is.null(en_cfg) && length(en_cfg) > 0 && nzchar(trimws(en_cfg[1]))) {
  # B-2未知酶值报错而非静默回退 trypsin（suppressWarnings 会
  # 消化警告、使报错分支不可达）
  w_enr <- character(0)
  en <- withCallingHandlers(enzyme_resolve(en_cfg[1]),
                            warning = function(ww) {
                              w_enr <<- c(w_enr, conditionMessage(ww)); invokeRestart("muffleWarning")
                            })
  if (length(w_enr) > 0)
    stop("未知酶切类型：", trimws(en_cfg[1]), "\nunknown enzyme type: ", trimws(en_cfg[1]))
  en_src <- "config_step1 键 enzyme"
} else {
  en <- "trypsin"
  # v1.0键缺失 = 模板即真相——第二步不再改写引擎
  # 模板；此处提示内部口径按 trypsin 计数，非胰蛋白酶酶切样本须显式填键
  message("提示：config_step1 未配置 enzyme 键——第二步搜库保留各引擎模板酶定义；",
          "管线内部漏切/N 端非特异判据按 trypsin 口径。若样本非胰蛋白酶酶切，请显式填写 enzyme 键。\n",
          "note: no enzyme key in config_step1 - step 2 keeps each engine template's enzyme; ",
          "the in-pipeline missed-cleavage/N-terminal rules use trypsin. Fill the enzyme key ",
          "explicitly for non-tryptic samples.")
}
Sys.setenv(ISO_ENZYME = en)
cat("酶切类型 enzyme:", en, "（来源:", en_src, "）\n")

# ---- 结果文件夹（pGlycoQuant 风格 result_YYYYMMDDHHMMSS，用户 规格）----
# 全部流程产出统一写在该文件夹内（清洗件 Competitive_search_result/、merge/、
# 1.unlabeled_result/、2.labeled_result/）；已存在 result_* 时取最新一个（断点
# 续跑就地复用），删除/移走全部 result_* 后重跑即生成新的时间戳文件夹
# ---- result folder (pGlycoQuant-style result_YYYYMMDDHHMMSS, user 
# spec) - every pipeline output lives inside it (cleaned copies under
# Competitive_search_result/, merge/, 1.unlabeled_result/, 2.labeled_result/);
# an existing result_* folder means the newest one is reused (resume in place),
# and removing all result_* folders before a rerun yields a fresh timestamped one
# RESULT_DIR 环境变量优先同项目多套口径并行重跑时
# 指定各自新文件夹，避免复用旧 merge 造成"假重跑"（get_result_dir 优先
# RESULT_DIR，未设时才按最新 result_* 定位）
# the RESULT_DIR environment variable takes priority 
# parallel re-runs with different settings get their own fresh folders instead
# of resuming the old merge (get_result_dir prefers RESULT_DIR, then falls
# back to the newest result_* folder)
result_dir <- get_result_dir(cfg_dir)
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
cat("结果文件夹：", result_dir, "\n")
cat("result folder:", result_dir, "\n")

# ---- 运行遥测接线TELEMETRY=1 时在每步边界写 run_telemetry 标记 ----
# 采样器由 run_telemetry.R start 独立启动（5 秒间隔）；此处只在步骤边界打 mark，
# 遥测汇总（telemetry_summary.txt）即得每步 CPU/内存峰值均值与外存占用。
# TELEMETRY=1 wires run_telemetry marks at every step boundary; the sampler is a
# separate `run_telemetry.R start` process (5 s interval); marks give the per-step
# CPU/memory peaks and disk usage in telemetry_summary.txt.
tele_path <- file.path(script_dir, "run_telemetry.R")
tele_mark <- function(lbl) {
  # B4 修复cmd 的 "set TELEMETRY=1 &" 值带尾空格（"1 "），
  # identical("1 ","1") 恒假曾致打点静默失效；此处 trimws 比较。
  # 目录固定用 cfg_dir（采样器约定 start 指向 config 目录；RESULT_DIR 指向
  # 结果子目录，marks 写到那里会与采样器文件分离、attach 聚不到）。
  # B4 fix cmd's "set TELEMETRY=1 &" carries a trailing space
  # ("1 "), so identical("1 ","1") was always FALSE and marks silently died;
  # compare after trimws. The directory is always cfg_dir (the sampler
  # convention is `start <config dir>`; RESULT_DIR points at the result
  # subfolder, where marks would split from the sampler file and never
  # aggregate on attach).
  if (!identical(trimws(Sys.getenv("TELEMETRY")), "1")) return(invisible(NULL))
  td <- cfg_dir
  st <- suppressWarnings(system2(rscript, c(shQuote(tele_path), "mark",
                                            shQuote(lbl), shQuote(td)),
                                 stdout = FALSE, stderr = FALSE))
  if (!isTRUE(st == 0)) {
    cat("[遥测] mark 写入失败 (exit ", st, ")：", lbl, "\n", sep = "")
    cat("[telemetry] mark failed (exit ", st, "): ", lbl, "\n", sep = "")
  }
  invisible(NULL)
}
if (identical(trimws(Sys.getenv("TELEMETRY")), "1")) {
  cat("[遥测] TELEMETRY=1 步骤边界打点已启用（采样器目录：", cfg_dir, "）\n", sep = "")
  cat("[telemetry] TELEMETRY=1 step-boundary marks enabled (sampler dir: ",
      cfg_dir, ")\n", sep = "")
}

# ---- B5 代码指纹防新旧混合态----
# 断点续跑复用旧结果文件夹；若代码已更新而旧结果仍在，下游会按新口径消费
# 旧输出形成混合态。启动时对比核心代码指纹：不一致时醒目告警（由用户决定
# 删除全部 result_* 重跑，还是接受混合态继续）；首次运行记录当前指纹。
# ---- B5 code fingerprint against stale-resume mixing ----
# resume reuses the existing result folder; after a code update the stale
# outputs are consumed under the new rules and form a mixed state. Compare
# the core-code fingerprint at startup: a mismatch prints a loud warning (the
# user decides whether to delete all result_* folders and rerun, or accept
# the mixed state); the first run records the current fingerprint.
fp_file <- file.path(result_dir, "code_fingerprint.txt")
fp_now <- code_fingerprint()
fp_old <- if (file.exists(fp_file)) readLines(fp_file, warn = FALSE) else character(0)
if (length(fp_old) > 0 && nzchar(fp_old[1]) && identical(fp_old[1], fp_now)) {
  cat("[代码指纹] 与上次运行一致，断点续跑安全。\n")
  cat("[code fingerprint] matches the last run; resume is safe.\n")
} else if (length(fp_old) > 0 && nzchar(fp_old[1])) {
  cat("############################################################\n")
  cat("# [警告] 核心代码自本结果文件夹生成后已更新！\n")
  cat("#        断点续跑会把旧结果与新代码混用（新旧口径混合态）。\n")
  cat("#        建议：删除全部 result_* 文件夹后重跑；\n")
  cat("#        如确认继续，本次运行结果可能口径不一致。\n")
  cat("# [WARNING] the core code changed since this result folder\n")
  cat("#          was produced! Resuming mixes stale outputs with\n")
  cat("#          new code. Consider deleting all result_* folders\n")
  cat("#          and rerunning; continuing risks mixed semantics.\n")
  cat("############################################################\n")
} else {
  writeLines(fp_now, fp_file)
  cat("[代码指纹] 已记录（首次运行）：", substr(fp_now, 1, 16), "...\n")
  cat("[code fingerprint] recorded (first run):", substr(fp_now, 1, 16), "...\n")
}

# PDF_ON 状态留档（网页报告"谱图 PDF"卡片显示用；read_pdf_on_flag 读取）
# record the PDF_ON state (shown by the web report's spectrum-PDF card)
writeLines(if (identical(pdf_on, "1")) "PDF_ON=1" else "PDF_ON=0",
           file.path(result_dir, "pdf_on_flag.txt"))

# 单组样本流程：目录不体现组名
# single-sample workflow: plain folders, no group names anywhere
step1_done <- file.path(cfg_dir, "unlabeled", "step1_done.txt")
peptide_txt <- file.path(result_dir, "merge", "peptide.txt")
cfg2_path   <- file.path(cfg_dir, "config_step2.txt")
cfg3_path   <- file.path(cfg_dir, "config_step3.txt")

# 以子进程方式运行各步驱动并注入环境变量；子进程继承终端输出，进度实时可见
# run one step driver as a subprocess with the given environment variables; the subprocess
# inherits stdout/stderr so its bilingual progress is live
# 函数：run_step —— 以子进程方式运行某一步驱动脚本。
# 作用：注入给定的环境变量、实时转发子进程输出，并检查退出码。
# 输入：tag = 步骤显示名；script = 待运行的 R 脚本路径；envs = 环境变量命名向量。
# 输出：无显式返回值；成功时打印完成标记，失败时以 stop 报错。
# Function: run_step -- run one step driver script as a subprocess.
# Purpose: inject the given environment variables, forward the child output live, and check the exit code.
# Inputs: tag = step display name; script = path of the R script to run; envs = named vector of environment variables.
# Output: no explicit return value; prints a completion marker on success and stops with an error on failure.
run_step <- function(tag, tag_en, script, envs) {
  # 模块开始横幅：附带墙钟时间戳（网页壳子据此解析运行时间线）；
  # 中文/英文标签分两行独立打印，便于解析与展示
  # step start banner with a wall-clock timestamp (parsed by the web shell into
  # the run timeline); the Chinese and English labels print on separate lines
  cat("\n==== [总程序 / Master] ", tag, " 开始 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
  cat("==== [Master] ", tag_en, " start @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
  tele_mark(paste0(tag, "_START"))
  # 保存即将被覆盖的环境变量的旧值，以便退出时恢复
  # Save the old values of the environment variables about to be overwritten, for later restore
  old <- Sys.getenv(names(envs), unset = NA)
  names(old) <- names(envs)
  # 无论函数如何退出，都把环境变量恢复为原值
  # Restore the environment variables to their original values no matter how this function exits
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else {
        v <- old[[nm]]
        names(v) <- nm
        do.call(Sys.setenv, as.list(v))
      }
    }
  }, add = TRUE)
  # 注入本步所需的环境变量后启动子进程并检查退出码
  # Inject the environment variables needed by this step, launch the subprocess, and check its exit code
  do.call(Sys.setenv, as.list(envs))
  st <- system2(rscript, shQuote(script))
  # 退出码 7＝停点约定（子进程等待用户输入后正常停止）：总程序打印停点并正常收尾，不再继续后续步骤
  # exit status 7 = the stop-point convention (the subprocess stopped for user input):
  # the master prints the stop point and ends cleanly, without running the later steps
  if (st == 7) {
    tele_mark(paste0(tag, "_STOP"))
    cat("==== [总程序 / Master] ", tag, " 停点（等待用户输入）@ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
    cat("==== [Master] ", tag_en, " stop point (awaiting user input) @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
    quit(status = 0)
  }
  # 非零退出码表示该步失败，立即报错中止
  # A nonzero exit code means the step failed: report it and abort
  if (st != 0) {
    tele_mark(paste0(tag, "_FAIL"))
    stop("[总程序 / Master] ", tag, " 失败 (exit ", st, ")\n[Master] ", tag_en, " failed (exit ", st, ")")
  }
  tele_mark(paste0(tag, "_END"))
  cat("==== [总程序 / Master] ", tag, " 完成 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
  cat("==== [Master] ", tag_en, " done @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
}

cat("==== [总程序 / Master] 无标记 isoPSV 三步流程总入口 ====\n")
cat("==== [Master] three-step pipeline master entry ====\n")
# P3-1/N-12启动打印 VERSION，日志可自证版本
vf <- file.path(script_dir, "..", "VERSION")
if (file.exists(vf)) {
  ver_line <- readLines(vf, warn = FALSE, n = 1L)
  cat("[版本] ", ver_line, "\n[version] ", ver_line, "\n", sep = "")
}
cat("CONFIG:", cfg_path, "\n")
cat("参数文件目录：", cfg_dir, "\n")
cat("config dir: ", cfg_dir, "\n")

# ---- 第一步：完成标记（unlabeled/step1_done.txt）已存在则跳过 ----------
# ---- step 1: skip when the completion marker (unlabeled/step1_done.txt) exists
# 纯蛋白质基因组学/已注释矫正模式的搜库前早退输入
# psm_file 若无目标侧行（pg=1 无新肽候选 / pg=2 无已注释谱），任何搜库都不
# 必做，直接写说明文件并正常结束（退出码 0）。早退判断需要 Proteins 列；
# 无该列时回退到第三步的门控路径。
pg_mode_u <- "0"
if (!is.null(cfg1$proteogenomics_mode) && length(cfg1$proteogenomics_mode) > 0 &&
    tail(cfg1$proteogenomics_mode, 1) %in% c("0", "1", "2")) {
  pg_mode_u <- tail(cfg1$proteogenomics_mode, 1)
}
if (pg_mode_u %in% c("1", "2")) {
  psm_u <- resolve_cfg_paths(cfg1, "psm_file")[1]
  if (nzchar(psm_u) && file.exists(psm_u)) {
    hdr_u <- strsplit(readLines(psm_u, n = 1), "\t")[[1]]
    if ("Proteins" %in% hdr_u) {
      cc_u <- rep("NULL", length(hdr_u))
      cc_u[which(hdr_u == "Proteins")] <- "character"
      psmu <- read.delim(psm_u, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                         colClasses = cc_u, na.strings = character(0))
      fa_u <- character(0)
      if (!is.null(cfg1$fasta) && length(cfg1$fasta) > 0) {
        fa_f_u <- resolve_cfg_paths(cfg1, "fasta")[1]
        if (nzchar(fa_f_u) && file.exists(fa_f_u)) fa_u <- unique(names(read_fasta(fa_f_u)))
      }
      seq_u <- if ("Sequence" %in% hdr_u) as.character(psmu$Sequence)
               else rep("", nrow(psmu))
      in_fa_u <- fa_annotated_proteins(as.character(psmu$Proteins), seq_u,
                                       if (nzchar(fa_f_u)) fa_f_u else character(0))
      n_target_u <- if (pg_mode_u == "1") sum(!in_fa_u) else sum(in_fa_u)
      if (n_target_u == 0L) {
        dir.create(file.path(result_dir, "1.unlabeled_result"),
                   showWarnings = FALSE, recursive = TRUE)
        sumf_u <- file.path(result_dir, "1.unlabeled_result", "proteogenomics_summary.txt")
        msg_u <- if (pg_mode_u == "1") {
          sprintf("纯蛋白质基因组学模式：输入 psm_file 的 %d 行全部为已注释肽段、无新肽段候选，任何搜库与竞争验证都不必执行，流程直接结束。\npure proteogenomics mode: all %d rows of the psm_file are annotated - no novel candidates, nothing to search or validate; the run ends here.",
                  nrow(psmu), nrow(psmu))
        } else {
          sprintf("已注释肽段矫正模式：输入 psm_file 的 %d 行全部为新肽段候选、无已注释肽段，无需矫正，流程直接结束。\nannotated-correction mode: all %d rows of the psm_file are novel candidates - nothing annotated to correct; the run ends here.",
                  nrow(psmu), nrow(psmu))
        }
        writeLines(msg_u, sumf_u, useBytes = TRUE)
        cat(msg_u, "\n")
        quit(save = "no", status = 0)
      }
      cat(sprintf("[蛋白质基因组学] 搜库前预检：目标侧行 %d / 总行 %d\n",
                  n_target_u, nrow(psmu)))
    }
  }
}
if (file.exists(step1_done)) {
  cat("[跳过] 第一步：完成标记已存在 -> ", step1_done, "\n")
  cat("[Skip] step 1: step-1 marker exists -> ", step1_done, "\n")
} else {
  run_step("[用户第一步] 无标记谱图提取",
           "[User Step 1] unlabeled-spectrum extraction",
           file.path(script_dir, "unlabeled_extract.R"),
           c(CONFIG = cfg_path))
}

# 第二步需要用户填写的 search_file 条目
# step 2 needs the user-filled search_file entries
n_sf <- 0L
if (!file.exists(cfg2_path)) {
  # 第一步被跳过而第二步参数文件缺失时，自动重新生成（保留自动发现的 pFind 模板键）
  # when step 1 is skipped but the step-2 parameter file is missing, regenerate
  # it (the auto-discovered pFind template keys come along)
  cat("重新生成第二步参数文件：", cfg2_path, "\n")
  cat("regenerating the step-2 parameter file: ", cfg2_path, "\n")
  make_step2_config(cfg_path, cfg1, cfg2_path)
}
# 统计用户已在第二步参数文件中填写的 search_file 条目数
# Count the search_file entries the user has already filled into the step-2 parameter file
if (file.exists(cfg2_path)) {
  cfg2 <- read_param_file(cfg2_path)
  sf <- cfg2$search_file
  if (!is.null(sf)) n_sf <- sum(!is.na(sf) & sf != "")
}
if (n_sf == 0L) {
  # 搜库方式取自用户第一步参数文件的 auto_search（已链式写入第二步文件）：1=自动、0=手动
  # search mode comes from the user's step-1 parameter file (chained into the
  # step-2 file): auto_search=1 runs pFind automatically, 0 is manual
  auto_on <- FALSE
  for (src in list(cfg1, if (exists("cfg2", inherits = FALSE)) cfg2 else NULL)) {
    if (!is.null(src) && !is.null(src$auto_search) && length(src$auto_search) > 0 &&
        tail(src$auto_search, 1) == "1") auto_on <- TRUE
  }
  # pfind_exe：用户填 pFind 安装路径（版本目录自动补 bin），这里规范化并验证 Searcher.exe 存在
  # pfind_exe: the user fills the pFind install path (a version folder gets bin appended automatically); it is normalized here and Searcher.exe is verified to exist
  pfind_ok <- FALSE
  for (src in list(if (exists("cfg2", inherits = FALSE)) cfg2 else NULL, cfg1)) {
    if (!is.null(src) && !is.null(src$pfind_exe) && length(src$pfind_exe) > 0 &&
        tail(src$pfind_exe, 1) != "") {
      # v100a：pfind_exe 自愈——resolve_software_path 返回的目录已含
      # Searcher.exe（levels_up=0 → 签名文件所在目录），直接拼文件名，
      # 不再重复套 "bin"（早期实现 file.path(resolved,"bin","Searcher.exe")
      # 在 resolved=…/bin 时拼出 …/bin/bin/Searcher.exe 导致失败）
      pf.raw <- tail(src$pfind_exe, 1)
      pf.resolved <- tryCatch({
        resolve_software_path(pf.raw, "Searcher.exe", levels_up = 0L)
      }, error = function(e) "")
      exe <- if (nzchar(pf.resolved)) {
        file.path(pf.resolved, "Searcher.exe")
      } else {
        try(normalize_pfind_exe(pf.raw), silent = TRUE)
      }
      if (!inherits(exe, "try-error") && nzchar(exe)) pfind_ok <- TRUE
      break
    }
  }
  # has_auto：自动搜库模式可用（auto_search=1 且 pFind 可执行文件可用；
  # pfind_cfg_open/pfind_cfg_restrict 缺键时由第二步从内置模板自动生成，无需停点）
  # has_auto: automatic search is usable (auto_search=1 with the pFind executable
  # available); a missing pfind_cfg_open/pfind_cfg_restrict is auto-generated by
  # step 2 from the built-in template, so there is no stop point
  has_auto <- auto_on && pfind_ok
  if (!has_auto) {
    cat("\n[等待输入] 请填写 ", cfg2_path,
        if (auto_on) {
          " 的 pfind_exe（自动搜库模式，缺 pFind 版本目录——填含 bin/Searcher.exe 的那一级目录，如 C:/pFindStudio/pFind3 或 D:/Bioinformatics/pFind3.1）"
        } else {
          " 的 search_file（手动搜库模式：四个补充搜库结果的绝对路径）；如需软件自动搜库，请在 config_step1.txt 把 auto_search 设为 1 并填写 pfind_exe"
        },
        "，然后重新运行本程序。\n")
    cat("[Awaiting input] fill in ", cfg2_path,
        if (auto_on) {
          " pfind_exe (auto-search mode, pFind executable missing)"
        } else {
          " search_file (manual mode: absolute paths of the four supplementary search results); for auto search set auto_search = 1 in config_step1.txt and fill pfind_exe"
        },
        ", then rerun this program.\n")
    quit(status = 0)
  }
}

# ---- 第二步：pDeep2 输入 peptide.txt 已存在则跳过 -----------------
# ---- step 2: skip when the pDeep2 peptide.txt exists -----------------------
# 功能排队①：跳过判据升级为"策略集合一致"——peptide.txt 存在但
# 当前 config_step2 的策略集合（search_file+类型）与侧车不符（新增/删改策略、
# 或早期实现本产物无侧车）时重跑步骤 2；引擎搜索结果经参数指纹复用，仅重算合并。
# sidecar = merge/peptide.txt.strategies，由 user_step2 写、与
# strategy_identity_lines 同一实现比对。
if (file.exists(peptide_txt)) {
  sidecar <- paste0(peptide_txt, ".strategies")
  strategies_ok <- file.exists(sidecar) &&
    identical(trimws(readLines(sidecar, warn = FALSE)),
              strategy_identity_lines(cfg2))
  if (!strategies_ok) {
    cat("[检测] 策略集合变化（或策略侧车缺失），重跑第二步（引擎结果按参数指纹复用）\n")
    cat("[Detect] strategy set changed (or sidecar missing), rerunning step 2 (engine results reused via fingerprints)\n")
    run_step("[用户第二步] revise 修正 + 多引擎合并",
             "[User Step 2] revise + multi-engine merge",
             file.path(script_dir, "user_step2.R"),
             c(CONFIG2 = cfg2_path, RESULT_DIR = result_dir))
  } else {
    cat("[跳过] 第二步：pDeep2 输入 peptide.txt 已存在且策略集合未变 -> ", peptide_txt, "\n")
    cat("[Skip] step 2: pDeep2 input peptide.txt exists and the strategy set is unchanged -> ", peptide_txt, "\n")
  }
} else {
  run_step("[用户第二步] revise 修正 + 多引擎合并",
           "[User Step 2] revise + multi-engine merge",
           file.path(script_dir, "user_step2.R"),
           c(CONFIG2 = cfg2_path, RESULT_DIR = result_dir))
}

# ---- 第三步需要用户填写的 predict_file 条目 -----------------------
# ---- step 3 needs the user-filled predict_file entry -----------------------
if (!file.exists(cfg3_path)) {
  # 第二步被跳过且第三步参数文件缺失时：重新生成
  # （predict_file 占位），停下让用户填写
  # step 2 was skipped and the step-3 parameter file is missing: regenerate
  # it (predict_file placeholder) and stop for the user to fill it
  cfg2 <- read_param_file(cfg2_path)
  make_step3_config(cfg2_path, cfg2, cfg3_path)
} else {
  # 参数文件存在但 merge_result/predict_file 指向已不存在的文件时（删除
  # result_* 文件夹重跑的典型场景）：按当前结果文件夹重建 cfg3，并保留已填的
  # weights_A/threshold_A；失效的 predict_file 行先清掉
  # when the parameter file exists but its merge_result/predict_file points at
  # a missing file (the typical case after removing the result_* folder and
  # rerunning): rebuild cfg3 against the current result folder, keeping any
  # filled weights_A/threshold_A; a stale predict_file line is dropped first
  # (fix: self-healing, no manual step)
  cfg3_chk <- read_param_file(cfg3_path)
  mr <- resolve_cfg_paths(cfg3_chk, "merge_result")
  mr <- mr[!is.na(mr) & mr != ""]
  pf <- resolve_cfg_paths(cfg3_chk, "predict_file")
  pf <- pf[!is.na(pf) & pf != ""]
  if ((length(mr) == 0 || !file.exists(tail(mr, 1))) ||
      (length(pf) > 0 && !file.exists(tail(pf, 1)))) {
    cat("第三步参数文件的 merge_result/predict_file 已失效，按当前结果文件夹重建：", cfg3_path, "\n")
    cat("merge_result/predict_file in the step-3 file is stale; rebuilding against the current result folder: ", cfg3_path, "\n")
    if (length(pf) > 0 && !file.exists(tail(pf, 1))) {
      cfg3_lines <- readLines(cfg3_path, warn = FALSE)
      writeLines(cfg3_lines[!grepl("^predict_file\t", cfg3_lines)],
                 cfg3_path, useBytes = TRUE)
    }
    cfg2 <- read_param_file(cfg2_path)
    make_step3_config(cfg2_path, cfg2, cfg3_path)
  }
}
# 读取第三步参数文件；use_pdeep2=1（自动，predict_file 已由上方自动执行写入）、
# =2（手动，用户填 predict_file）时确认 predict_file 已填；use_pdeep2=0 时跳过
# 该停点——预测指标按缺省处理（缺失不填 0、predict_missing 标注、投票降级）
# Read the step-3 parameter file; when use_pdeep2=1 (auto, predict_file already
# written by the auto-run above) or =2 (manual, the user fills predict_file)
# verify predict_file is filled; when use_pdeep2=0 this stop point is skipped -
# the predict metric is treated as missing (never 0, flagged by predict_missing,
# and dropped from the votes)
cfg3 <- read_param_file(cfg3_path)
use_pd2 <- if (!is.null(cfg1$use_pdeep2) && length(cfg1$use_pdeep2) > 0) tail(cfg1$use_pdeep2, 1) else "1"
pd2_exe <- cfg1$pdeep2_exe
pd2_exe <- if (!is.null(pd2_exe)) tail(pd2_exe[!is.na(pd2_exe) & pd2_exe != ""], 1) else ""
# v100：pdeep2_exe 自愈解析（用户填任意层级都能找到正确的 pDeep2 根）
if (nzchar(pd2_exe)) {
  pd2_resolved <- resolve_pdeep2_exe(pd2_exe)
  if (nzchar(pd2_resolved) && !identical(pd2_resolved, pd2_exe)) {
    cat("[pDeep2] 路径自愈：", pd2_exe, " -> ", pd2_resolved, "
", sep = "")
    cat("[pDeep2] path healed: ", pd2_exe, " -> ", pd2_resolved, "
", sep = "")
    pd2_exe <- pd2_resolved
  }
}
pd2_inst <- if (!is.null(cfg1$instrument)) tail(cfg1$instrument,1)[1] else "Lumos"
pd2_inst <- if (!is.null(pd2_inst) && length(pd2_inst) > 0 && nzchar(tail(pd2_inst, 1))) tail(pd2_inst, 1) else "Lumos"
# N-1/P0pep_txt 此前只在自动 pDeep2 块内定义，手动模式/
# use_pdeep2=0/断点续跑时下游预测器循环（length(pids)>0 分支）引用即
# object 'pep_txt' not found。提升到块外无条件定义。
# pep_txt: step-2 合并肽表路径（自动 pDeep2 的输入、各预测器适配器的输入）
pep_txt <- file.path(result_dir, "merge", "peptide.txt")
# ---- pDeep2 自动执行（use_pdeep2=1 且 config1 填了 pdeep2_exe）----
# 对第二步的 peptide.txt 运行 pDeep2 生成 predict.txt，并自动写入第三步参数文件——
# 第三步不再停点；predict.txt 已存在则复用（断点续跑）。安装目录布局按安装文档：
# <pdeep2_exe>\py37\python.exe ＋ <pdeep2_exe>\pDeep\pDeep2\predict.py
# ---- auto-run pDeep2 (use_pdeep2=1 and pdeep2_exe filled in config1) ---------
# runs pDeep2 on the step-2 peptide.txt, writes predict.txt and patches it into the
# step-3 parameter file - step 3 no longer stops; an existing predict.txt is reused
# (resume). Install layout per the install record: <pdeep2_exe>\py37\python.exe plus
# <pdeep2_exe>\pDeep\pDeep2\predict.py
  if (identical(use_pd2, "1") && nzchar(pd2_exe)) {
  # 预测器分文件夹制——pDeep2 产物迁入
  # merge/predictors/pdeep2/（peptides.txt＝本轮提供的输入列表＋predict.txt＋
  # 指纹侧车）。旧布局 merge/predict.txt 指纹匹配时仍复用（断点续跑不重算），
  # 新预测一律写新位置。
  #: per-predictor folders - pDeep2 output moves to
  # merge/predictors/pdeep2/ (peptides.txt = the list fed this run, predict.txt,
  # fingerprint sidecar). A legacy merge/predict.txt with a matching fingerprint
  # is still reused (resume); fresh predictions go to the new location.
  pd2_dir <- file.path(result_dir, "merge", "predictors", "pdeep2")
  predict_out_legacy <- file.path(result_dir, "merge", "predict.txt")
  predict_out <- file.path(pd2_dir, "predict.txt")
  # P1-2复用判据从"文件存在"升级为
  # 肽表指纹——新增策略/换库后 peptide.txt 变化但 predict.txt 残留时，旧预测
  # 会静默缺列（阶段 A 按缺失处理）。侧车记 peptide.txt 的 mtime/size。
  # predict_sidecar_file：predict.txt 的指纹侧车；peptide_file_stat：肽表文件信息；
  # predict_fingerprint_lines：本次应匹配的指纹行；predict_reusable：可复用判定
  fingerprint_matches <- function(out_path) {
    sc <- paste0(out_path, ".run_params")
    file.exists(out_path) && file.exists(sc) &&
      identical(trimws(readLines(sc, warn = FALSE)), predict_fingerprint_lines)
  }
  peptide_file_stat <- file.info(pep_txt)
  predict_fingerprint_lines <- c(paste0("pep=", format(peptide_file_stat$mtime, "%Y-%m-%d %H:%M:%OS3"), "|", peptide_file_stat$size),
              paste0("inst=", pd2_inst))
  reusable_new <- fingerprint_matches(predict_out)
  reusable_legacy <- fingerprint_matches(predict_out_legacy)
  predict_reusable <- reusable_new || reusable_legacy
  if (predict_reusable) predict_out <- if (reusable_new) predict_out else predict_out_legacy
  if (!predict_reusable) {
    if (file.exists(predict_out))
      cat("[pDeep2] 检测到肽表/安装变化，将重新预测\n[pDeep2] peptide-table change detected, re-predicting\n")
    # python 路径三级解析：环境变量 PDEEP2_PYTHON（不同电脑可覆盖）→ 安装文档
    # 标准布局 <pdeep2_exe>\py37\python.exe；P1-7再补 py38/py39/
    # py310 变体与 pDeep2/predict.py 备选位置（Python 3.7 已 EOL）
    # python resolution: PDEEP2_PYTHON env override -> standard layout
    # <pdeep2_exe>\py37\python.exe -> py38/py39/py310 variants 
    py_exe <- Sys.getenv("PDEEP2_PYTHON")
    if (!nzchar(py_exe) || !file.exists(py_exe)) {
      # candidate_env_dir：候选 python 环境子目录；
      # candidate_python：候选中命中的解释器
      for (candidate_env_dir in c("py37", "py38", "py39", "py310")) {
        candidate_python <- file.path(pd2_exe, candidate_env_dir, "python.exe")
        if (file.exists(candidate_python)) { py_exe <- candidate_python; break }
      }
      if (!file.exists(py_exe)) py_exe <- file.path(pd2_exe, "py37", "python.exe")
    }
    # v99h：附件包布局——python37 是 pDeep2 的兄弟目录（predictor_env/python37）
    # pd2_exe=…/predictor_env/pDeep2 → python37 在 …/predictor_env/python37
    if (!file.exists(py_exe)) {
      py_exe <- file.path(dirname(pd2_exe), "python37", "python.exe")
    }
    pd_script <- file.path(pd2_exe, "pDeep", "pDeep2", "predict.py")
    if (!file.exists(pd_script)) pd_script <- file.path(pd2_exe, "predict.py")
    if (!file.exists(py_exe) || !file.exists(pd_script)) {
      stop("[总程序 / Master] pDeep2 安装目录不完整：", pd2_exe,
           "（python 3.7：", py_exe, "；predict.py：", pd_script,
           "）。请按安装文档的标准布局安装，或设置环境变量 PDEEP2_PYTHON 指向 python 3.7。\n",
           "[Master] pDeep2 install incomplete under: ", pd2_exe,
           " (python 3.7: ", py_exe, "; predict.py: ", pd_script,
           "). Follow the standard install layout, or set the PDEEP2_PYTHON environment variable to the python 3.7 executable.")
    }
    cat("==== [自动执行 pDeep2] 对以下文件运行预测：", pep_txt, " ====\n")
    cat("==== [Auto pDeep2] predict on: ", pep_txt, " ====\n")
    cat("python:", py_exe, "\n")
    #：本轮提供给 pDeep2 的输入列表落其独立文件夹（peptides.txt）
    #: the list fed to pDeep2 this run lands in its own folder
    dir.create(pd2_dir, recursive = TRUE, showWarnings = FALSE)
    pd2_list <- file.path(pd2_dir, "peptides.txt")
    if (!file.copy(pep_txt, pd2_list, overwrite = TRUE))
      pd2_list <- pep_txt   # 复制失败（如只读盘）退回主表，不阻断
    old_wd2 <- getwd()
    setwd(dirname(pd_script))
    st2 <- system2(py_exe, c(shQuote(pd_script), "-e", "0.3", "-i", pd2_inst,
                             "-in", shQuote(pd2_list), "-out", shQuote(predict_out)))
    setwd(old_wd2)
    if (st2 != 0 || !file.exists(predict_out)) stop("pDeep2 运行失败\npDeep2 run failed")
    writeLines(predict_fingerprint_lines, paste0(predict_out, ".run_params"),
               useBytes = TRUE)   # P1-2 肽表指纹侧车（ 修复：predict_sidecar_file 已内联）
  } else {
    cat("[跳过] pDeep2 predict.txt 已存在且肽表未变 -> ", predict_out, "\n")
    cat("[Skip] pDeep2 predict.txt exists and peptide table unchanged -> ", predict_out, "\n")
  }
  # 把 predict_file 自动写入第三步参数文件（已有键替换、缺键追加）
  # patch predict_file into the step-3 parameter file (replace or append)
  cfg3_lines <- readLines(cfg3_path, warn = FALSE)
  i_pf <- which(grepl("^predict_file\t", cfg3_lines))
  if (length(i_pf) > 0) {
    cfg3_lines[i_pf[1]] <- paste0("predict_file\t", predict_out)
  } else {
    cfg3_lines <- c(cfg3_lines, paste0("predict_file\t", predict_out))
  }
  writeLines(cfg3_lines, cfg3_path, useBytes = TRUE)
  cfg3 <- read_param_file(cfg3_path)
  cat("[自动执行 pDeep2] predict_file 已写入：", cfg3_path, "\n")
  cat("[Auto pDeep2] predict_file written to: ", cfg3_path, "\n")
}
# v1.0多理论谱图预测器——每个启用的额外预测器（config_step1
# predictorN_* 键）跑其 L1 适配器（code/predictors/pred_<type>.R），产物与
# pDeep2 predict.txt 同构； 起各入独立文件夹 merge/predictors/<type>/；
# predictor_primary 决定 stage A 的 predict 列（predict_file 回写）；失败只记
# manifest、不打断整轮。
# v1.0 multi-predictor loop - each enabled extra predictor
# (config_step1 predictorN_* keys) runs its L1 adapter; outputs are pDeep2-
# shaped predict files; predictor_primary selects the stage-A predict column;
# failures only mark the manifest and the column is treated as missing.
#：STEP1 质谱仪键（instrument，默认 Lumos）作为全体预测器的默认
# 仪器（适配器 params 留空即用它；pDeep3 仅 QE/Lumos，其余映射 Lumos）
Sys.setenv(ISOPSV_INSTRUMENT = toupper(
  if (!is.null(cfg1$instrument)) tail(cfg1$instrument, 1)[1] else "Lumos"))
pn <- names(cfg1)
pids <- unique(sub("^predictor([0-9]+)_.*$", "\\1", pn[grepl("^predictor[0-9]+_", pn)]))
pids <- pids[order(suppressWarnings(as.integer(pids)))]
if (length(pids) > 0 && file.exists(pep_txt)) {
  manifest <- character(0)
  for (pid in pids) {
    getp <- function(k) {
      v <- cfg1[[paste0("predictor", pid, "_", k)]]
      if (is.null(v) || length(v) == 0) "" else trimws(tail(v, 1))
    }
    # N-9typ 拼进 pred_<typ>.R 路径前先消毒（先 tolower
    # 再剔非 [a-z0-9_]——先 gsub 会删掉大写字母）；含 ../ 的恶意键无法再
    # 指到任意脚本。N-10：mode 键缺失时 getp 返回 ""，此前会默认执行——改为
    # 只在 mode 明确为 "1" 时启用（模板语义：1=启用，0=停用）
    typ <- gsub("[^a-z0-9_]", "", tolower(getp("type"))); mode <- getp("mode")
    if (!nzchar(typ) || !identical(mode, "1") || identical(typ, "pdeep2")) next
    py <- getp("python"); params <- getp("params")
    # （设计约定）：预测器路径一律"软件安装文件夹"级，不到 python.exe
    # 级——自愈链扩展两级：pDeep3 安装夹的同级 pDeep2/py37（本地
    # D:\Bioinformatics 布局）与附件包嵌套 python37。
    #: predictor paths are at the install-folder level, never
    # python.exe - the heal chain gains two candidates: the pDeep2/py37
    # sibling of a pDeep3 install folder (D:\Bioinformatics layout) and the
    # nested addon python37.
    py_raw <- py
    if (nzchar(py) && !grepl("[.](exe|bat)$", tolower(py))) {
      # 依次尝试：<path>/python.exe、<path>/python37、<path>/py37、
      # dirname(<path>)/python37、dirname(<path>)/pDeep2/py37
      py_cands <- c(file.path(py, "python.exe"),
                    file.path(py, "python37", "python.exe"),
                    file.path(py, "py37", "python.exe"),
                    file.path(dirname(py), "python37", "python.exe"),
                    file.path(dirname(py), "pDeep2", "py37", "python.exe"))
      py_hit <- py_cands[file.exists(py_cands)]
      if (length(py_hit) > 0) {
        cat("[预测器] python 路径自愈：", py, " -> ", py_hit[1], "
", sep = "")
        cat("[predictor] python path healed: ", py, " -> ", py_hit[1], "
", sep = "")
        py <- py_hit[1]
      }
    }
    #：用户填的是 pDeep3 安装文件夹时，直接把 PDEEP3_ROOT 指过去——
    # driver 不再靠 python 兄弟目录猜根（附件包/自定义布局两可）。
    #: when the user gave the pDeep3 install folder, point
    # PDEEP3_ROOT straight at it - the driver no longer has to guess the root
    # from the python siblings (works for both addon and custom layouts).
    if (typ == "pdeep3" && nzchar(py_raw) && dir.exists(py_raw) &&
        (file.exists(file.path(py_raw, "tune_and_predict.py")) ||
         dir.exists(file.path(py_raw, "pDeep")))) {
      Sys.setenv(PDEEP3_ROOT = normalizePath(py_raw, winslash = "/"))
      cat("[预测器] pDeep3 安装文件夹识别：PDEEP3_ROOT=", py_raw, "\n", sep = "")
      cat("[predictor] pDeep3 install folder detected: PDEEP3_ROOT=", py_raw, "\n", sep = "")
    }
    adapter <- file.path(getwd(), "predictors", paste0("pred_", typ, ".R"))
    if (!file.exists(adapter)) {
      manifest <- c(manifest, sprintf("%s\tFAIL\t适配器缺失: %s", typ, adapter))
      cat("[预测器] ", typ, " 适配器缺失，跳过：", adapter, "\n", sep = "")
      next
    }
    # 预测器分文件夹制——每个理论谱图预测器独立目录
    # merge/predictors/<type>/：peptides.txt＝本轮提供给该预测器的输入列表
    # （各预测器格式不同，由适配器在此列表基础上再派生自己的中间格式）、
    # predict.txt＝产物、run.log＝运行日志。旧布局 merge/predict_<type>.txt
    # 不再读写（残留无害）。
    #: per-predictor folders - each predictor gets its own
    # merge/predictors/<type>/ folder: peptides.txt = the input list fed to
    # THIS predictor this run (adapters derive their own formats from it),
    # predict.txt = output, run.log = the run log. The legacy flat
    # merge/predict_<type>.txt layout is no longer read or written.
    pred_dir <- file.path(result_dir, "merge", "predictors", typ)
    dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
    pred_list <- file.path(pred_dir, "peptides.txt")
    if (!file.copy(pep_txt, pred_list, overwrite = TRUE)) pred_list <- pep_txt
    out_pred <- file.path(pred_dir, "predict.txt")
    # N-3iso_system2 内部对每个参数再 shQuote，这里
    # 再包一层 shQuote 会让含空格路径（如 C:\Program Files\...）双重引号坏掉。
    # 一律传原样路径，交给 iso_system2 统一加引号。
    a <- c(adapter, "--in", pred_list, "--out", out_pred)
    if (nzchar(py)) a <- c(a, "--python", py)
    if (nzchar(params)) a <- c(a, "--params", params)
    cat("[预测器] 运行 ", typ, "：", adapter, "\n", sep = "")
    # P0-Cstdout=TRUE 会把输出向量当返回值，多行时
    # if(st==0) 抛 length>1——改输出到日志文件、只取退出码（ 入预测器夹）
    predictor_log <- file.path(pred_dir, "run.log")
    predictor_status <- iso_system2(rscript, a, stdout = predictor_log, stderr = predictor_log)
    if (predictor_status == 0 && file.exists(out_pred)) {
      manifest <- c(manifest, sprintf("%s\tOK\t%s", typ, out_pred))
    } else {
      manifest <- c(manifest, sprintf("%s\tFAIL\texit %s", typ, predictor_status))
      cat("[预测器] ", typ, " 失败（exit ", predictor_status, "），该预测列按缺失处理\n", sep = "")
    }
  }
  writeLines(manifest, file.path(result_dir, "merge", "predictor_run_manifest.txt"),
             useBytes = TRUE)
  # 主预测器：非 pDeep2 时用其产物回写 predict_file（stage A 的 predict 列）
  # N-8回写只在自动模式（use_pdeep2=1）生效——手动模式
  # （=2）predict_file 是用户手填值，无条件替换会覆盖它；=0 时 predict 按缺失
  # 处理，照写 config_step3 与之口径不一。非自动模式只提示不写入。
  pp <- cfg1$predictor_primary
  if (!is.null(pp) && length(pp) > 0) {
    # P3-2pp 与 typ 同源消毒——既拼进
    # predict_<pp>.txt 路径又作 grepl 正则（^ 锚定）；gsub 后只剩 [a-z0-9_]，
    # 路径恒落在 merge 下、正则无元字符
    pp <- gsub("[^a-z0-9_]", "", tolower(trimws(tail(pp, 1))))
    if (nzchar(pp) && pp != "pdeep2") {
      #：主预测器产物在 merge/predictors/<pp>/predict.txt
      out_pp <- file.path(result_dir, "merge", "predictors", pp, "predict.txt")
      # P0-2仅查文件存在会消费上一轮残留的陈旧产物
      # （本轮适配器失败时旧 predict_<type>.txt 仍在磁盘上）——必须以本轮
      # manifest 的 OK 行为准；pp 未配置/本轮失败时保持 pDeep2
      pp_ok <- any(grepl(paste0("^", pp, "\tOK\t"), manifest))
      if (!identical(use_pd2, "1")) {
        cat("[预测器] use_pdeep2 != 1（手动/关闭），主预测器 ", pp,
            " 不回写 predict_file（保持现状）\n", sep = "")
      } else if (pp_ok && file.exists(out_pp)) {
        cfg3_lines <- readLines(cfg3_path, warn = FALSE)
        i_pf <- which(grepl("^predict_file\t", cfg3_lines))
        if (length(i_pf) > 0) cfg3_lines[i_pf[1]] <- paste0("predict_file\t", out_pp)
        else cfg3_lines <- c(cfg3_lines, paste0("predict_file\t", out_pp))
        writeLines(cfg3_lines, cfg3_path, useBytes = TRUE)
        cfg3 <- read_param_file(cfg3_path)
        cat("[预测器] 主预测器 ", pp, " 的 predict_file 已写入：", out_pp, "\n", sep = "")
      } else {
        cat("[预测器] 主预测器 ", pp, " 本轮未成功（未配置/运行失败），predict_file 保持 pDeep2\n", sep = "")
      }
    }
  }
}
if (!identical(use_pd2, "0")) {
  pf <- cfg3$predict_file
  if (is.null(pf)) pf <- character(0)
  pf <- pf[!is.na(pf) & pf != ""]
  if (length(pf) == 0L) {
    if (identical(use_pd2, "1")) {
      # use_pdeep2=1（自动模式）但未填 pdeep2_exe：提示两种解法
      # use_pdeep2=1 (auto mode) without pdeep2_exe: suggest both fixes
      cat("\n[等待输入] 请填写 ", cfg_path,
          " 的 pdeep2_exe（use_pdeep2=1 自动模式需要 pDeep2 安装目录）；或把 use_pdeep2 改为 2，再填写 ", cfg3_path,
          " 的 predict_file（手动运行 pDeep2 生成的 predict.txt 绝对路径），然后重新运行本程序。\n")
      cat("[Awaiting input] fill in ", cfg_path,
          " pdeep2_exe (the auto mode use_pdeep2=1 needs the pDeep2 install folder); or set use_pdeep2 = 2 and fill in ", cfg3_path,
          " predict_file (absolute path of the predict.txt generated by running pDeep2 manually), then rerun this program.\n")
    } else {
      # use_pdeep2=2（手动模式）：等待用户对 merge\peptide.txt 运行 pDeep2 后填写 predict_file
      # use_pdeep2=2 (manual mode): await predict_file after the user runs pDeep2 on merge\peptide.txt
      cat("\n[等待输入] 请填写 ", cfg3_path,
          " 的 predict_file（对 merge\\peptide.txt 手动运行 pDeep2 生成的 predict.txt 绝对路径），然后重新运行本程序。\n")
      cat("[Awaiting input] fill in ", cfg3_path,
          " predict_file (absolute path of the predict.txt generated by running pDeep2 on merge\\peptide.txt manually), then rerun this program.\n")
    }
    quit(status = 0)
  }
} else {
  cat("[跳过] pDeep2 理论预测谱未启用（use_pdeep2=0）：预测指标按缺省处理\n")
  cat("[Skip] pDeep2 predicted spectra not used (use_pdeep2=0): predict metric treated as missing\n")
}

# ---- 第三步（始终运行）＋后处理扩充 ---------------------------------
# ---- step 3 (always) + extend ----------------------------------------------
run_step("[用户第三步] 阶段 A 评估",
         "[User Step 3] stage A evaluation",
         file.path(script_dir, "user_step3.R"),
         c(CONFIG3 = cfg3_path, RESULT_DIR = result_dir))

# extend 自动定位 SAP 清单（merge/SAP_spectra_index.txt，第二步生成）
# extend locates the SAP index automatically (merge/SAP_spectra_index.txt,
# generated by step 2) - no user-provided SAP file
run_step("[后处理] report 参考格式扩充",
         "[Post-process] report reference-format extension",
         file.path(script_dir, "extend_report_A.R"),
         c(CONFIG3 = cfg3_path, RESULT_DIR = result_dir))

# ---- 阶段 A 结果详细分类（后处理之后执行，纯表格运算、无 MGF 加载）----
# 对 report_A.tsv 的全部谱图逐谱给出结论状态：无竞争匹配（策略一致）/
# 存在竞争但原始结果正确/原始结果错误（附错误原因与修饰差异）/
# 混合谱/谱图质量差难区分/子离子完全一致难区分
# ---- stage A detailed classification (runs after the post-process; pure
# table computation, no MGF loading) - labels every report_A.tsv spectrum:
# no competition / original correct / original wrong (with error reason and
# modification detail) / chimeric / poor-quality ambiguous / identical-ions
# ambiguous
run_step("[阶段A分类] 阶段 A 结果分类",
         "[Stage A classification] stage-A result classification",
         file.path(script_dir, "stageA_classify.R"),
         c(EXT_REPORT = file.path(result_dir, "1.unlabeled_result", "report_A.tsv"),
           OUT = file.path(result_dir, "1.unlabeled_result", "classification_A.tsv"),
           UNLABELED_DIR = file.path(cfg_dir, "unlabeled")))

# 纯蛋白质基因组学模式
# 第三步已把已注释谱写入 pg_skipped.tsv——此处按清单补记到分类表，结论
# annotated_skipped、指标列留空（已注释谱不做错误校正）。
# pure proteogenomics mode: step 3 wrote the annotated spectra into
# pg_skipped.tsv - append them here to the classification table with the
# annotated_skipped conclusion and empty metric columns.
pg_skipped_f <- file.path(result_dir, "1.unlabeled_result", "pg_skipped.tsv")
if (file.exists(pg_skipped_f)) {
  sk <- readLines(pg_skipped_f, warn = FALSE)
  sk <- unique(sk[nzchar(sk)])
  if (length(sk) > 0) {
    clf <- file.path(result_dir, "1.unlabeled_result", "classification_A.tsv")
    if (file.exists(clf)) {
      cdat <- read.delim(clf, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE, na.strings = character(0))
      scol <- grep("conclusion|status", colnames(cdat), value = TRUE)[1]
      if (length(scol) > 0 && !is.na(scol) && "File_Name" %in% colnames(cdat)) {
        add <- cdat[rep(1, length(sk)),, drop = FALSE]
        for (cc in colnames(add)) add[[cc]] <- rep(NA, length(sk))
        add[["File_Name"]] <- sk
        # 补记结论按工作流程模式：纯蛋白质基因组学=已注释跳过；
        # 已注释矫正=新肽跳过
        skip_lab <- if (identical(pg_mode_u, "2")) "novel_skipped" else "annotated_skipped"
        add[[scol]] <- skip_lab
        write.table(rbind(cdat, add), clf, sep = "\t", row.names = FALSE,
                    col.names = TRUE, quote = FALSE)
        cat(sprintf("[流程门控] 分类表补记跳过行 %d（%s）\n", length(sk), skip_lab))
      }
    }
  }
}

# ---- 阶段 A 谱图 PDF（无标记侧判定：判错＝candidate_misidentification、
# 存疑＝打分差近阈值/为 0；输出到 1.unlabeled_result）----
# 默认关闭；
# 设 PDF_ON=1 显式开启。
# ---- stage-A spectrum PDF (unlabeled-side judgment: error = the
# candidate_misidentification decision, ambiguous = the by-score delta near
# the threshold or zero; written into 1.unlabeled_result)
# OFF by default (user ruling: PDF export not needed for now,
# pending a redesign); set PDF_ON=1 to enable.
label_type <- if (!is.null(cfg1$label_type) && length(cfg1$label_type) > 0) {
  tolower(tail(cfg1$label_type, 1))
} else "unlabeled"
# ---- 无标记样本两分法报告与标记样本同名同构产出
# 1.annotation_correction_report.tsv（已注释矫正）与
# 2.novel_peptide_competition_report.tsv（新肽段竞争）＋ two_part_report_summary
# ——落 1.unlabeled_result/（v1.0 口径：report 路径），网页两分法面板自动
# 呈现；可信新肽段三列（novel_doubt/novel_status/credible_novel）按设计约定
# 修订与标记与否无关，同六重存疑规则计算（数据源＝阶段 A
# 逐谱候选＋第一步无标记谱；见 unlabeled_two_part_report.R 头注释）。
# P2本步必须排在 PDF-A 之前——第三类 PDF（可信
# 新肽段）依赖本报告的 credible_novel 列；顺序颠倒时 PDF 读到上一轮旧报告、
# 渲染 0 页且文件不落地。
# ---- unlabeled two-part report same names and
# structure as the labeled samples - annotation-correction + novel-competition
# TSVs plus the summary, written into 1.unlabeled_result/ (v1.0 ruling:
# report path) so the web two-part panel renders them unchanged; the credible-
# novel columns are computed with the SAME six doubt rules (ruling:
# independent of labeling; stage-A per-spectrum candidates + step-1 unlabeled
# spectra - see the unlabeled_two_part_report.R header).
# P2 this step MUST run before PDF-A -
# the third PDF (credible novel peptides) depends on this report's
# credible_novel column; reversed, the PDF reads the previous round's report,
# renders 0 pages and writes no file.
# （无标记通道=0）：全标记样本同样走无标记侧两分法报告
ul_ch_ur <- if (!is.null(cfg1$unlabeled_channel) && length(cfg1$unlabeled_channel) > 0)
  tail(cfg1$unlabeled_channel, 1) else "1"
if (identical(label_type, "unlabeled") || identical(ul_ch_ur, "0")) {
  run_step("[两分法报告] 无标记样本两分法报告（矫正/新肽段）",
           "[two-part report] unlabeled two-part report (correction/novel)",
           file.path(script_dir, "unlabeled_two_part_report.R"),
           c(CONFIG3 = cfg3_path, RESULT_DIR = result_dir))
}
if (identical(pdf_on, "1")) run_step("[谱图 PDF-A] 阶段 A 错误匹配与存疑谱图导出",
         "[spectrum PDF-A] stage-A error-match and ambiguous-match export",
         file.path(script_dir, "spectrum_pdf_export_A.R"),
         c(CONFIG3 = cfg3_path, RESULT_DIR = result_dir))

# ---- 阶段 B：成对标记谱验证（标记样本）----------------------
# ---- stage B: paired labeled-spectrum verification (labeled samples) -------
# 阶段 B 只对标记样本有意义：成对标记谱匹配＋五指标多数投票＋成对余弦判定；
# 报告已存在时跳过；标记样本缺 all_mgf 时停点让用户填写（与 predict 停点一致）。
# stage B only applies to labeled samples: paired labeled-spectrum matching +
# five-metric majority vote + paired-cosine decision; skipped when the report
# already exists; for labeled samples a missing all_mgf stops for user input
# (same await-input pattern as predict_file).
# （无标记通道=0）：全标记样本没有无标记侧可配对——阶段 B 与冲突
# 报告跳过（与无标记样本同口径），下游走无标记侧两分法报告（ul_ch_ur 已在上方定义）。
if (!identical(label_type, "unlabeled") && !identical(ul_ch_ur, "0")) {
  stb_out <- file.path(result_dir, "2.labeled_result", paste0("paired_report_", label_type, ".tsv"))
  all_mgf <- cfg1$all_mgf
  if (is.null(all_mgf)) all_mgf <- character(0)
  all_mgf <- all_mgf[!is.na(all_mgf) & all_mgf != ""]
  if (file.exists(stb_out)) {
    cat("[跳过] 阶段 B：成对匹配报告已存在 -> ", stb_out, "\n")
    cat("[Skip] stage B: paired report exists -> ", stb_out, "\n")
  } else if (length(all_mgf) == 0L) {
    cat("\n[等待输入] 请填写 ", cfg_path,
        " 的 all_mgf（原始 MGF 文件夹，含无标记＋标记谱），然后重新运行本程序。\n")
    cat("[Awaiting input] fill in ", cfg_path,
        " all_mgf (original MGF folder with unlabeled + labeled spectra), then rerun this program.\n")
    quit(status = 0)
  } else {
    run_step("[阶段 B] 成对标记谱匹配",
             "[Stage B] paired labeled-spectrum matching",
             file.path(script_dir, "stageB_paired_match.R"),
             c(CONFIG3 = cfg3_path, ALL_MGF = tail(all_mgf, 1),
               RESULT_DIR = result_dir))
  }
  # ---- 冲突谱结论标注（阶段 B 完成后执行，纯表格运算、无 MGF 加载）----
  # 在 report_A.tsv 冲突行的基础上合并阶段 B 的新增列，逐谱给出结论状态
  # （混合谱/子离子完全一致难区分/投票并列难区分/谱图质量差难区分/原始结果
  # 错误/存在竞争但原始结果正确），错误行附错误原因与修饰差异细节
  # ---- conflict-spectrum conclusion annotation (runs after stage B; pure
  # table computation, no MGF loading) - merges the stage-B columns onto the
  # report_A conflict rows and labels each spectrum with a conclusion status
  # (chimeric / identical-ions ambiguous / vote-tie ambiguous / poor-quality
  # ambiguous / original result wrong / original result correct); wrong rows
  # carry the error reason and modification-difference detail
  run_step("[冲突报告] 冲突谱结论标注",
           "[Conflict report] conflict-spectrum conclusion annotation",
           file.path(script_dir, "stageB_conflict_report.R"),
           c(EXT_REPORT = file.path(result_dir, "1.unlabeled_result", "report_A.tsv"),
             PAIRED = stb_out,
             OUT = file.path(result_dir, "2.labeled_result", paste0("conflict_report_", label_type, ".tsv")),
             UNLABELED_DIR = file.path(cfg_dir, "unlabeled")))
} else {
  cat("[跳过] 阶段 B：无标记样本\n")
  cat("[Skip] stage B: unlabeled sample\n")
  # v1.0 修复：无标记主表接线引用 all_mgf，但该变量只在标记分支（上方
  # if 块内）赋值——无标记跑到主表步即 object not found 崩溃。此处提前按同口径赋值（可为空向量）。
  # v1.0 fix: the unlabeled master-table call references all_mgf, which
  # was assigned only inside the labeled branch above - unlabeled runs
  # crashed with object-not-found at the master step . Assign it here with the same rules (may be empty).
  all_mgf <- cfg1$all_mgf
  if (is.null(all_mgf)) all_mgf <- character(0)
  all_mgf <- all_mgf[!is.na(all_mgf) & all_mgf != ""]
}

# ---- 全谱主表导出（仅标记样本；无标记样本跳过）----
# 生成与 8.衣藻全谱主表.txt 同格式的 master_table.tsv：候选 1–5 全部参数、
# 判定结论、蛋白定位、差异类型。无标记样本暂不支持（阶段 A 报告为 candidate_*
# 列名体系，与主表所依赖的成对报告列名 unlabel_ions_match.* 之间尚无列映射；
# 未完成功能，起按设计约定跳过）——无标记样本的科学需求已由
# stageA/report_A.tsv 与 stageA/classification_A.tsv 两张表覆盖。
# ---- master-table export (labeled samples only; unlabeled samples are
# skipped) - writes master_table.tsv in the same format as 8.衣藻全谱主表.txt:
# all candidate-1..5 parameters, the judgment conclusion, protein mapping and
# the difference type. Unlabeled samples are not supported yet (the stage-A
# report uses the candidate_* column system while the master table expects the
# paired-report columns such as unlabel_ions_match.*; the column mapping is
# missing - unfinished feature, skipped since ) - the unlabeled needs
# are fully covered by stageA/report_A.tsv and stageA/classification_A.tsv.
# （无标记通道=0）：全标记样本无成对报告，主表（依赖成对列）一并
# 跳过——科学需求由 report_A/classification_A 与两分法报告覆盖。
if (identical(ul_ch_ur, "0")) {
  cat("[跳过] 主表导出：无标记通道的全标记样本暂不支持（依赖成对报告列）。\n")
  cat("[Skip] master-table export: not available for fully labeled samples",
      "without an unlabeled channel (depends on paired-report columns).\n")
} else if (label_type == "unlabeled") {
  # v1.0无标记流程也出主表——由阶段 A 的
  # scan_details_A（逐谱×逐候选）透视 58 列宽表（cand1..cand5＋判定列，
  # 标记侧列留空），落在 1.unlabeled_result/master_table.tsv。
  # (user ruling): unlabeled runs also get a master table - pivots
  # the stage-A scan_details_A into a 58-column wide table (cand1..cand5 +
  # conclusion columns, labeled columns blank) at
  # 1.unlabeled_result/master_table.tsv.
  run_step("[无标记主表] 阶段 A 全谱主表（58 列）",
           "[unlabeled master] stage-A full-spectrum master table (58 cols)",
           file.path(script_dir, "unlabeled_master_export.R"),
           c(CONFIG3 = cfg3_path, ALL_MGF = tail(all_mgf, 1),
             RESULT_DIR = result_dir))
} else {
  run_step("[主表导出] 全谱主表导出",
           "[master table] full-spectrum master-table export",
           file.path(script_dir, "master_table_export.R"),
           c(CONFIG3 = cfg3_path, ALL_MGF = tail(all_mgf, 1),
             RESULT_DIR = result_dir))
}

# ---- 谱图 PDF 导出（标记样本：阶段 B 与冲突报告完成后自动执行）----
# 用 Fig.7 定稿判定模型对成对报告逐谱判定，自动渲染两类 PDF：
#   error_match_spectra.pdf ＝ 判需校正谱；ambiguous_matches_spectra.pdf ＝ 难区分带谱
# 默认关闭；
# 设 PDF_ON=1 显式开启（PDF_LIMIT 仍可限制张数，自检用）
# ---- spectrum PDF export (labeled samples: runs automatically after stage B
# and the conflict report) - applies the finalized Fig.7 judgment model to the
# paired report and renders two PDFs: error_match_spectra.pdf (need-correction
# spectra) and ambiguous_matches_spectra.pdf (hard-to-judge band)
# OFF by default (user ruling: PDF export not needed for now,
# pending a redesign); set PDF_ON=1 to enable (PDF_LIMIT still caps pages)
if (exists("stb_out", inherits = FALSE) && label_type != "unlabeled" &&
    file.exists(stb_out) && identical(pdf_on, "1")) {
  run_step("[谱图 PDF] 错误匹配与存疑谱图导出",
           "[spectrum PDF] error-match and ambiguous-match spectrum export",
           file.path(script_dir, "spectrum_pdf_export.R"),
           c(CONFIG3 = cfg3_path, ALL_MGF = tail(all_mgf, 1),
             RESULT_DIR = result_dir))
}

# ---- SAP 索引留档到最终结果文件夹（用户 规格：SAP 单独记录、
# 不影响理论谱预测、结果可还原）----
# 第二步已在 merge/ 与 Competitive_search_result/ 写出 SAP_spectra_index.txt
#（原序列＋"pos,XXX->YYY[site]" 修饰＋revised_peptide/ref_protein 映射）；
# 这里再复制一份到本项目的最终结果文件夹，使任何下游 SAP 替换序列都能在
# 结果旁直接还原为原始形式
# ---- SAP index archived next to the final results (user spec:
# SAP recorded separately, prediction unaffected, results restorable) ----
# step 2 already wrote SAP_spectra_index.txt to merge/ and
# Competitive_search_result/ (original sequence + "pos,XXX->YYY[site]"
# modification + revised_peptide/ref_protein mapping); copy it once more into
# this project's final-result folders so every downstream substituted SAP
# sequence can be restored to its original form right next to the results
sap_idx <- file.path(result_dir, "merge", "SAP_spectra_index.txt")
if (file.exists(sap_idx)) {
  for (rf in c("1.unlabeled_result", "2.labeled_result")) {
    rfd <- file.path(result_dir, rf)
    if (dir.exists(rfd)) {
      file.copy(sap_idx, file.path(rfd, "SAP_spectra_index.txt"), overwrite = TRUE)
      cat("[SAP 索引留档]", rfd, "\n")
      cat("[SAP index archived]", rfd, "\n")
    }
  }
}

# ---- v1.0 多预测器余弦叠加------------------
# 各预测器 predict.txt → 逐谱算余弦 → master_table 加 cand1_predict_<type> 列
# 列名直接含软件名（cand1_predict_pdeep2 / cand1_predict_pdeep3）
tryCatch({
  ov_script <- file.path(getwd(), "predict_cosine_overlay.R")
  if (file.exists(ov_script)) {
    cat("\n==== [多预测器余弦] 开始 ====\n")
    system2(file.path(R.home("bin"), "Rscript"),
            c(shQuote(ov_script), paste0("--result_dir=", shQuote(result_dir))),
            stdout = "", stderr = "")
    cat("==== [多预测器余弦] 完成 ====\n")
  }
}, error = function(e) {
  cat("[多预测器余弦] 跳过：", conditionMessage(e), "\n")
})

cat("\n==== [总程序 / Master] 全部完成 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n")
cat("==== [Master] all done @ ", format(Sys.time(), "%H:%M:%S"), " ====\n")
