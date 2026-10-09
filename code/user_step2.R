# =============================================================================
# user_step2.R - 用户第二步驱动程序：读取第二步参数文件（四个
# user_step2.R - user step 2 driver: read the step-2 parameter file (four
# 补充搜库结果）与基础配置（策略 1），
# supplementary search results) together with the base config (strategy 1),
# 随后按以下已确认的流程执行：
# then per the confirmed flow:
#   1. 校正四个原始搜库结果（生产版修饰校正，revise_search_result——
#   1. REVISE the four original search results (production modification
#      1.modification_exchange.R 的移植版），
#      corrections, revise_search_result - port of 1.modification_exchange.R)
#      并把 _revised.spectra 副本写到 merge/revised/；
#      and write the _revised.spectra copies to merge/revised/;
#   2. 将校正后的结果与策略 1 合并为
#   2. MERGE the revised results with strategy 1 into
#      1.combine_result.spectra；
#      1.combine_result.spectra;
#   3. 生成 pDeep2 的 peptide.txt。
#   3. generate peptide.txt for pDeep2.
#
# 用法：
# Usage:
#   运行方式：Rscript user_step2.R
#   Rscript user_step2.R
# 第二步参数文件路径放在环境变量 CONFIG2 中。
# with the step-2 parameter file path in the environment variable CONFIG2.
# =============================================================================

# 加载通用工具函数（配置解析、路径解析、合并/校正/SAP 提取等）
# Load the shared utility functions (config parsing, path resolution, merge/revise/SAP extraction, etc.)
# 核心库统一入口（N-1，2026-09-25）：按 --file= 定位代码目录、绝对路径加载
# （boot.R），替代裸相对 source——此前仅当启动 cwd 恰为 code/ 时成立。
# core-library entry (N-1, 2026-09-25): locate the code dir via --file= and
# load by absolute path (boot.R), replacing the bare relative source().
sd0 <- dirname(tryCatch(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"),
    error = function(e) ""))
if (is.na(sd0) || !dir.exists(sd0)) sd0 <- getwd()
source(file.path(sd0, "boot.R"))

# 读取第二步参数文件（路径来自环境变量 CONFIG2），并解析基础配置（策略 1）
# Read the step-2 parameter file (path from env var CONFIG2) and the base config (strategy 1)
cfg2_path <- Sys.getenv("CONFIG2")
if (!nzchar(cfg2_path)) stop("set CONFIG2 to the step-2 parameter file path")
cfg2 <- read_param_file(cfg2_path)
# 记下本脚本所在目录（随客户端分发、含 pUniFind_compat_mods.txt 等资源文件）；
# run_pfind_auto 内部会临时切到 pFind bin 目录，故须提前保存
# remember this script's folder (shipped with the client, holds resources such as
# pUniFind_compat_mods.txt); run_pfind_auto temporarily changes to the pFind bin
# folder, so it must be captured up front
code_dir <- getwd()
base_path <- resolve_cfg_paths(cfg2, "base_config")[1]
base_cfg <- read_param_file(base_path)
# candidate_scope（2026-09-27 用户口径）：filtered=现状；top5=忽略引擎 1% FDR
# 过滤、每引擎每谱取引擎自身排序前 5 候选进竞争池。top5 下经 ISOPSV_CAND_ALL
# 环境变量让适配器子进程输出全部候选（合并步封顶 5，由 merge_search_results
# 的 candidate_scope 参数接收）。须在引擎自动运行循环之前设置——适配器经
# system2 子进程调用，PSOCK 集群 worker 继承本进程环境。
cand_scope <- "filtered"
for (src in list(base_cfg, cfg2)) {
  if (!is.null(src$candidate_scope) && length(src$candidate_scope) > 0 &&
      tail(src$candidate_scope, 1) %in% c("filtered", "top5")) {
    cand_scope <- tail(src$candidate_scope, 1)
  }
}
if (cand_scope == "top5") {
  Sys.setenv(ISOPSV_CAND_ALL = "1")
  cat("[候选口径] candidate_scope=top5：忽略引擎 FDR，每引擎每谱前 5 候选进竞争池\n")
  cat("[candidate scope] candidate_scope=top5: engine FDR ignored, top-5 candidates per engine per spectrum enter the pool\n")
}

# 解析 PSM 文件、补充搜库结果列表（仅 .spectra）与输出目录（merge/revised/）
# 输出目录位于结果文件夹 result_YYYYMMDDHHMMSS 内（用户 2026-09-05 规格：
# 全部流程产出统一入结果文件夹；RESULT_DIR 由总程序注入，独立运行按时间戳定位）
# Resolve the PSM file, the supplementary search-result list (.spectra only) and the output dir (merge/revised/)
# the output dir lives inside the result folder result_YYYYMMDDHHMMSS (user
# 2026-09-05 spec: all pipeline outputs go into the result folder; RESULT_DIR
# is injected by the master, standalone runs locate by timestamp)
psm_file <- resolve_cfg_paths(base_cfg, "psm_file")[1]
# search_file 先保留原始条目（文件或文件夹均可）；按 search_file_type 的类型
# 相关文件发现（文件夹 → 递归找目标文件、取最新）推迟到 search_types 确定之后。
# search_file keeps its raw entries first (a file OR a folder both work); the
# per-type file discovery (folder -> recursive target-file search, newest wins)
# runs after search_types is settled.
search_files <- resolve_cfg_paths(cfg2, "search_file")
# 多引擎版：并行键 search_file_type 指定每个补充搜库结果的策略类型
# （pfind_open/pfind_restrict/punifind_rescore/punifind_denovo/maxquant/
# msfragger_open/msfragger/msfragger_boost/sage/comet/xtandem/metamorpheus/peaks/mascot/pnovo/
# other）；未填时按文件名/目录规则推断（restrict 目录→pfind_restrict、
# denovo→punifind_denovo、文件名含 _punifind→punifind_rescore、其余→pfind_open），
# 旧配置零改动可用。槽位字母由策略目录表 strategy_slots 硬编码分配。
# multi-engine edition: the parallel key search_file_type names the strategy
# type of each supplementary result (pfind_open/pfind_restrict/punifind_rescore/
# punifind_denovo/maxquant/msfragger_open/msfragger/msfragger_boost/sage/comet/xtandem/
# metamorpheus/peaks/mascot/pnovo/other); when unfilled the file/folder naming
# rules infer it (restrict folder -> pfind_restrict, denovo -> punifind_denovo,
# _punifind in the name -> punifind_rescore, else pfind_open), so legacy configs
# work unchanged. Slot letters come from the hardcoded strategy catalog
# strategy_slots.
sf_type_cfg <- cfg2$search_file_type
sf_type_cfg <- if (!is.null(sf_type_cfg)) sf_type_cfg[!is.na(sf_type_cfg) & sf_type_cfg != ""] else character(0)
if (length(sf_type_cfg) > 0 && length(sf_type_cfg) != length(search_files)) {
  stop("search_file_type 行数须与 search_file 一致\nsearch_file_type must match search_file")
}
infer_strategy_type <- function(sf) {
  sf_dir <- tolower(basename(dirname(dirname(sf))))
  sf_name <- tolower(basename(sf))
  if (grepl("denovo", sf_dir, fixed = TRUE)) return("punifind_denovo")
  if (grepl("boost", sf_dir, fixed = TRUE)) return("msfragger_boost")
  if (grepl("pnovo", sf_dir, fixed = TRUE)) return("pnovo")
  if (grepl("restrict", sf_dir, fixed = TRUE)) return("pfind_restrict")
  if (grepl("_punifind", sf_name, fixed = TRUE)) return("punifind_rescore")
  "pfind_open"
}
search_types <- if (length(sf_type_cfg) > 0) sf_type_cfg else {
  if (length(search_files) > 0) vapply(search_files, infer_strategy_type, character(1))
  else character(0)
}
# 按策略类型的文件夹解析（用户 2026-09-11）：search_file 填文件夹时，按类型
# 对应的目标文件模式递归（≤4 层）发现结果文件，多个候选取最新，找不到报错；
# 填文件路径时原样通过。使 CLI 与网页壳都可用"文件或文件夹"配置。
# per-type folder resolution (user 2026-09-11): a directory in search_file is
# resolved by the type-specific file pattern (recursive, depth <= 4, newest
# wins, clear error when nothing matches); a file path passes through. Both
# the CLI and the web shell accept a file OR a folder.
resolve_per_type <- function(sf, ty) {
  # v1.0：auto 标记直通——auto-search 段会替换为真实产物路径
  # v1.0: auto markers pass through; the auto-search section replaces
  # them with real result paths
  if (startsWith(sf, "(auto:") || startsWith(basename(sf), "(auto:")) return(sf)
  if (dir.exists(sf)) {
    pat <- switch(ty,
      maxquant = "msms\\.txt$",
      msfragger = , msfragger_open = , msfragger_boost = "(\\.tsv$|\\.pepXML$|\\.pep\\.xml$)",
      sage = "results\\.sage\\.tsv$",
      comet = "(pep\\.tsv$|\\.txt$)",
      metamorpheus = "AllPeptides\\.psmtsv$",
      xtandem = "\\.xml$",
      peaks = "(\\.csv$|\\.tsv$)",
      mascot = "(\\.csv$|\\.dat$)",
      pnovo = "\\.csv$",
      "\\.spectra$")
    fs <- unique(list.files(sf, pattern = pat, recursive = TRUE, full.names = TRUE))
    fs <- fs[!dir.exists(fs)]
    if (length(fs) == 0)
      stop("search_file 文件夹内未找到 ", ty, " 结果文件（模式 ", pat, "）：", sf,
           "\nno ", ty, " result file found in folder (pattern ", pat, "): ", sf)
    if (length(fs) > 1) {
      fs <- fs[which.max(file.mtime(fs))]
      cat("search_file 文件夹内多个候选，取最新 newest wins:", fs, "\n")
    }
    return(fs)
  }
  if (file.exists(sf)) return(sf)
  stop("search_file 不存在 missing:", sf)
}
if (length(search_files) > 0) {
  search_files <- mapply(resolve_per_type, search_files, search_types, USE.NAMES = FALSE)
}
result_dir <- get_result_dir(cfg2$config_dir)
out_dir <- file.path(result_dir, "merge")
rev_dir <- file.path(out_dir, "revised")
dir.create(rev_dir, recursive = TRUE, showWarnings = FALSE)
# denovo CSV 转换（2026-09-23 网页版全流程冒烟修复）：punifind_denovo 类型的
# search_file 若填 GUI 产出的 .csv（无表头 9 列），先经 denovo_csv_to_spectra
# 转成 .spectra 再合并——此前 CSV 直接进合并，read 解析错列致 names 长度 0 崩溃
# (2026-09-23 web-shell smoke fix): a punifind_denovo search_file pointing at the
# GUI-produced .csv must first be converted by denovo_csv_to_spectra; the raw CSV
# went straight into the merge before, whose column parse broke with names[3] vs 0
# MiniMax P2-1（2026-09-26）：转换失败降级为跳过该 CSV 并打印详情，不整轮崩溃。
# 注意转换失败必须把该条目从 search_files/search_types 两个平行向量里一起移除
# （循环内删会错位，先记 drop_indices 循环后再删）——否则下游合并仍会读到坏 CSV。
# converted_spectra：转换产物；NULL = 转换失败（降级跳过该 CSV）
drop_indices <- integer(0)
search_files_orig <- search_files   # P3-5：原始清单留底（报错列被移除的源文件用）
# v1.0：从合并清单中移除 (auto:...) 标记——auto-search 段独立运行引擎并
# 产出策略文件，标记行不参与合并（防止 merge 读不到文件报错）
# v1.0: drop (auto:...) markers from the merge list; the auto-search
# section runs engines independently and produces strategy files
auto_marker_idx <- which(startsWith(basename(as.character(search_files)), "(auto:") |
                          startsWith(as.character(search_files), "(auto:"))
if (length(auto_marker_idx) > 0) {
  cat("[auto-search] 检测到 ", length(auto_marker_idx),
      " 个自动引擎标记，将从合并清单移除（auto-search 段独立处理）\n", sep = "")
  for (mi in auto_marker_idx) {
    cat("  marker: ", search_files[mi], " -> type: ", search_types[mi], "\n", sep = "")
  }
  search_files <- search_files[-auto_marker_idx]
  search_types <- search_types[-auto_marker_idx]
}
for (k in seq_along(search_files)) {
  if (identical(as.character(search_types[[k]]), "punifind_denovo") &&
      grepl("\\.csv$", search_files[[k]], ignore.case = TRUE)) {
    converted_spectra <- tryCatch(denovo_csv_to_spectra(search_files[[k]], rev_dir),
                                  error = function(e) NULL)
    if (is.null(converted_spectra)) {
      cat("punifind_denovo csv 转换失败，跳过：", search_files[[k]], "\n")
      drop_indices <- c(drop_indices, k)
    } else {
      search_files[[k]] <- converted_spectra
      cat("punifind_denovo csv -> spectra:", search_files[[k]], "\n")
    }
  }
}
if (length(drop_indices)) {
  search_files <- search_files[-drop_indices]
  search_types <- search_types[-drop_indices]
  cat("punifind_denovo 转换失败条目已从合并清单移除，剩余策略文件：",
      length(search_files), " 个\n")
  # P3-5（2026-09-26，deepseek v1.0 审查）：全部条目都被移除时立即报错并
  # 列出被移除文件——否则错误推迟到下游合并处，信息更远更难读
  if (length(search_files) == 0) {
    stop("punifind_denovo CSV 转换全部失败，合并清单为空；被移除文件：\n",
         paste(search_files_orig[drop_indices], collapse = "\n"))
  }
}
# Mascot / PEAKS 原始导出自动转换（2026-09-24）：这两个引擎只有 GUI 导出件
# （Mascot "Export Search Results" CSV、PEAKS DB search 肽段 CSV），合并前由
# 适配器转成交换格式。Mascot .dat 无法解析，给导出指引后停点。
# Mascot / PEAKS raw-export conversion (2026-09-24): these two engines only
# ship GUI exports (Mascot "Export Search Results" CSV, PEAKS DB-search peptide
# CSV); the adapters convert them to the interchange format before the merge.
# Mascot .dat cannot be parsed - stop with export guidance.
rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
for (k in seq_along(search_files)) {
  ty <- as.character(search_types[[k]])
  sf <- search_files[[k]]
  if (identical(ty, "mascot")) {
    if (grepl("\\.csv$", sf, ignore.case = TRUE)) {
      out_i <- file.path(rev_dir, "mascot_auto.spectra")
      st <- system2(rscript_bin, c(shQuote(file.path(code_dir, "search_adapter_mascot.R")),
                                   shQuote(sf), shQuote(out_i)))
      if (st != 0) stop("Mascot CSV 转换失败（exit ", st, "）：", sf,
                        "\nMascot CSV conversion failed (exit ", st, "): ", sf)
      search_files[[k]] <- out_i
      cat("mascot csv -> spectra:", out_i, "\n")
    } else {
      stop("Mascot 结果请导出为 CSV 后填入：在 Mascot 中 Export Search Results 选择 CSV，" ,
           "勾选 Query title、Peptide sequence 与 Variable modifications（.dat 文件无法直接解析）\n",
           "export the Mascot result as CSV (Export Search Results, tick Query title, ",
           "Peptide sequence and Variable modifications); .dat files cannot be parsed directly")
    }
  }
  if (identical(ty, "peaks") && grepl("\\.(csv|tsv)$", sf, ignore.case = TRUE)) {
    out_i <- file.path(rev_dir, "peaks_auto.spectra")
    st <- system2(rscript_bin, c(shQuote(file.path(code_dir, "search_adapter_peaks.R")),
                                 shQuote(sf), shQuote(out_i)))
    if (st != 0) stop("PEAKS CSV 转换失败（exit ", st, "）：", sf,
                      "\nPEAKS CSV conversion failed (exit ", st, "): ", sf)
    search_files[[k]] <- out_i
    cat("peaks csv -> spectra:", out_i, "\n")
  }
  if (identical(ty, "xtandem") && grepl("\\.xml$", sf, ignore.case = TRUE)) {
    out_i <- file.path(rev_dir, "xtandem_auto.spectra")
    # P1（2026-09-26，deepseek 实弹）：新版 X!Tandem 的 group label 是蛋白描述
    # 而非谱图标题，无 map 时按 label 回退会丢失谱级信息——同目录
    # <同名>_renum_map.tsv（new_scan<TAB>original_title）存在时一并传入，
    # 适配器按 <group id> 解析标题
    map_xt <- sub("\\.xml$", "_renum_map.tsv", sf, ignore.case = TRUE)
    xt_args <- c(shQuote(file.path(code_dir, "search_adapter_xtandem.R")),
                 shQuote(sf), shQuote(out_i))
    if (file.exists(map_xt)) xt_args <- c(xt_args, shQuote(map_xt))
    st <- system2(rscript_bin, xt_args)
    if (st != 0) stop("X!Tandem XML 转换失败（exit ", st, "）：", sf,
                      "\nX!Tandem XML conversion failed (exit ", st, "): ", sf)
    search_files[[k]] <- out_i
    cat("xtandem xml -> spectra:", out_i, "\n")
  }
}
# 打印本步骤运行信息：PSM 文件路径与搜库结果文件个数
# Print the run info for this step: the PSM file path and the number of search files
cat("==== [用户第二步] revise 修正 + 多策略合并 -> peptide.txt ====\n")
cat("==== [User Step 2] revise + multi-strategy merge -> peptide.txt ====\n")
cat("psm_file:", psm_file, "\n")
cat("search files:", length(search_files), "\n")

# 第二步各阶段的计时标记（双语输出）
# bilingual timing marks for the step-2 phases
pt0 <- proc.time()
# pt_mark：打印某个阶段的耗时标记（双语），并重置计时起点
# 作用：记录并打印从上次标记到当前的墙钟/CPU 耗时，便于跟踪各阶段运行时间
# 输入：tag——阶段名称（字符串，输出时作为标记打印）
# 输出：无返回值；副作用为打印耗时行并更新外层计时变量 pt0
# pt_mark: print an elapsed-time mark (bilingual) for a phase and reset the timer
# Purpose: record and print wall/CPU time since the previous mark, for tracking per-phase runtime
# Input: tag - the phase name (character, printed as the mark label)
# Output: no return value; side effects are the printed timing line and the updated outer pt0
pt_mark <- function(tag, tag_en = tag) {
  # 计算与上一标记之间的墙钟/CPU 耗时，打印后把计时起点推进到当前
  # Compute the wall/CPU time since the previous mark, print it, then advance the timer origin
  pt1 <- proc.time()
  d <- pt1 - pt0
  cat(sprintf("[计时] %s：已用 %.1fs / CPU %.1fs\n", tag, d[3], d[1]))
  cat(sprintf("[timing] %s: elapsed %.1fs / CPU %.1fs\n", tag_en, d[3], d[1]))
  pt0 <<- pt1
}

# ---- 1. 校正四个原始搜库结果（SAP 提取与校正共用同一次读入：每个文件只读一遍，
# 而不是两遍） -------------------------------------------------------------------
# ---- 1. revise the four original search results (SAP extraction shares the
# same read: one read per file instead of two) -------------------------------
fasta_path <- resolve_cfg_paths(base_cfg, "fasta")[1]

# ---- 自动搜库：用户未提供搜库结果时，软件按模板 cfg 自动跑 pFind（开放式＋限定式）。
# ---- 模板仅替换 msmspath1 / fastapath / outputpath / outputname，其余搜索参数保持用户保存值。
# Auto search: when the user has not supplied search results, run pFind (open + restricted
# multi-mod) from the template cfgs automatically. The template cfg only needs its msmspath1 /
# fastapath / outputpath / outputname replaced - everything else (search parameters) stays as
# the user saved it.
# run_pfind_auto：按模板 cfg 自动运行 pFind（开放式 / 限定式）
# 作用：改写模板 cfg（msmspath 逐 raw、fastapath、outputpath、outputname），
#       在 pFind 的 bin 目录运行 Searcher.exe，并返回结果 .spectra 文件路径；
#       结果已存在时直接复用（删除 result/ 可强制重搜）
# 输入：cfg_tpl——pFind 模板 cfg 文件路径（字符串）
# 输出：搜库结果 .spectra 文件路径（优先 Filtered；失败时抛错中止）
# run_pfind_auto: run pFind automatically from a template cfg (open / restrict)
# Purpose: rewrite the template cfg (per-raw msmspath, fastapath, outputpath, outputname),
#          run Searcher.exe from pFind's bin folder and return the result .spectra path;
#          an existing result is reused (delete result/ to force a re-search)
# Input: cfg_tpl - the pFind template cfg file path (character)
# Output: the search-result .spectra path (Filtered preferred; errors abort on failure)
run_pfind_auto <- function(cfg_tpl) {
  # outputname 取模板所在项目目录名（pFind_open / pFind_restrict），与参数文件路径保持一致；
  # 结果输出到该目录的 result（与参考 pFind 输出一致），自动生成的运行参数写在同套文件夹的 param 内
  # outputname follows the template's project folder name (pFind_open / pFind_restrict); results go
  # into its result/ subfolder (same as the reference pFind run), the auto-run cfg stays in param/
  outputname <- basename(dirname(dirname(cfg_tpl)))
  out_dir_auto <- file.path(dirname(dirname(cfg_tpl)), "result")
  # 搜库已完成则直接复用（删除 result/ 可强制重新搜库）——用户 GUI 完成 pUniFind 后重跑第二步时，
  # pFind 搜库不再重复跑（open 约 10 分钟、restrict 约 3 小时）
  # reuse a completed search (delete result/ to force a re-search) - when the user
  # reruns step 2 after finishing pUniFind in the GUI, the pFind searches are not
  # repeated (open ~10 min, restrict ~3 h)
  # 结果目录中已有 .spectra 时直接复用（优先 Filtered）
  # Reuse an existing .spectra in the result folder when present (Filtered preferred)
  cand0 <- list.files(out_dir_auto, pattern = "\\.spectra$", full.names = TRUE)
  # N-18（2026-09-26，GLM）：复用判据加 MGF 新鲜度——同指纹但 MGF 比产物新
  # （换了 raw 数据）时同样强制重搜；此前只查"存在+指纹"，静默复用旧结果。
  # 逐 raw 纯 MGF 清单（排除合并文件与 pFind *_check.mgf），下方搜索段复用。
  mgf_dir <- file.path(cfg2$config_dir, "unlabeled")
  mgfs <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
  mgfs <- mgfs[basename(mgfs) != "unlabeled.mgf" &
               !grepl("_check\\.mgf$", basename(mgfs))]
  # P2-2（2026-09-26）：复用前比对参数指纹——config 的 enzyme 键 + 模板文件
  # mtime/size + fasta 路径；任一变化都强制重跑（旧行为只看产物存在与否，
  # 改酶/改库后静默复用旧结果）
  rp_file <- file.path(out_dir_auto, "pfind.run_params")
  en_raw <- { v <- base_cfg$enzyme; if (!is.null(v) && length(v) > 0) trimws(v[1]) else "" }
  # B-11（2026-09-26）：指纹存解析后键（"Trypsin/P" 与 "trypsin" 同义，存原始串
  # 会无谓重跑 restrict ~3h）
  en_fp <- if (nzchar(en_raw)) suppressWarnings(tryCatch(enzyme_resolve(en_raw), error = function(e) "unknown")) else "<missing>"
  tpl_stat <- file.info(cfg_tpl)
  fp_lines <- c(paste0("enzyme=", en_fp),
                paste0("fasta=", fasta_path),
                paste0("tpl=", format(tpl_stat$mtime, "%Y-%m-%d %H:%M:%OS3"), "|", tpl_stat$size))
  rp_ok <- file.exists(rp_file) &&
    identical(trimws(readLines(rp_file, warn = FALSE)), fp_lines)
  # N-18（2026-09-26）：MGF 新鲜度（容差 1s）——MGF 比产物新则不复用
  mgf_newest <- if (length(mgfs)) max(file.info(mgfs)$mtime) else NA
  cand_mtime <- if (length(cand0)) max(file.info(cand0)$mtime) else NA
  mgf_fresh_ok <- is.na(mgf_newest) || is.na(cand_mtime) ||
    cand_mtime >= mgf_newest - 1
  if (length(cand0) > 0 && rp_ok && mgf_fresh_ok) {
    flt0 <- cand0[grepl("Filtered", basename(cand0))]
    cat("pFind 结果已存在，直接复用：", out_dir_auto, "\n")
    cat("reusing existing pFind result: ", out_dir_auto, "\n")
    return(if (length(flt0) > 0) tail(flt0, 1) else tail(cand0, 1))
  }
  if (length(cand0) > 0 && !rp_ok) {
    cat("pFind 检测到参数变化（enzyme/fasta/模板），将重新搜库\n")
    cat("pFind parameter change detected (enzyme/fasta/template), re-searching\n")
  } else if (length(cand0) > 0 && !mgf_fresh_ok) {
    cat("pFind 检测到 MGF 比产物新（原始谱数据变化），将重新搜库\n")
    cat("pFind MGF newer than the product (raw data changed), re-searching\n")
  }
  if (length(cand0) > 0 && (!rp_ok || !mgf_fresh_ok)) {
    # B-6（2026-09-26）：指纹触发的重搜须清除陈旧兼容过滤标记，否则新产物
    # 跳过 pUniFind 兼容修饰过滤、非兼容名流入 rescore
    unlink(file.path(out_dir_auto, ".compat_filter_done"))
  }
  # pFind.cfg 为 GBK 编码：先按字节读入（行结构正确）再 iconv 转 UTF-8 处理，写回时以 GBK 连接
  # 输出，保证全文件单一编码（混合编码会让 pFind 解析路径失败）
  # pFind.cfg is GBK-encoded: read as bytes (line structure intact), iconv to UTF-8 for the
  # replacements, write back through a GBK connection so the file stays single-encoded
  lines <- readLines(cfg_tpl, warn = FALSE)
  lines <- iconv(lines, from = "GBK", to = "UTF-8", sub = "")   # 非法字节容错（sub 防 NA 行）
  # 每 raw 一个 mgf（排除合并文件 unlabeled.mgf）：pUniFind 的 open 结果校正（rescore）
  # 只支持逐 raw 搜库格式——res 文件数量必须与 raw 数一致，合并文件格式 res 数不够
  # one mgf per raw (merged unlabeled.mgf excluded): pUniFind rescore only supports
  # the per-raw search layout - the res-file count must match the raw count, which
  # the single merged file cannot provide
  # （逐 raw 纯 MGF 清单已在上方复用判据处枚举——N-18 前置；此处沿用 mgfs）
  if (length(mgfs) == 0) stop("unlabeled 逐 raw MGF 缺失：", mgf_dir, "\nper-raw MGFs missing in: ", mgf_dir)
  # 改写模板 cfg：msmspath（逐 raw）、msmsnum、fastapath、outputpath、outputname
  # Rewrite the template cfg: msmspath (per raw), msmsnum, fastapath, outputpath, outputname
  i_ms <- which(grepl("^msmspath1=", lines))
  if (length(i_ms) == 1) {
    # i_ms 恰为最后一行时 (i_ms+1):length(lines) 是降序序列（如 6:5），须单独处理
    # when i_ms is the last line, (i_ms+1):length(lines) is a DESCENDING
    # sequence (e.g. 6:5) and must be handled separately
    tail_part <- if (i_ms < length(lines)) lines[(i_ms + 1):length(lines)] else character(0)
    lines <- c(lines[seq_len(i_ms - 1)],
               paste0("msmspath", seq_along(mgfs), "=", mgfs),
               tail_part)
  } else {
    lines <- c(lines, paste0("msmspath", seq_along(mgfs), "=", mgfs))
  }
  i_num <- which(grepl("^msmsnum=", lines))
  if (length(i_num) == 1) lines[i_num] <- paste0("msmsnum=", length(mgfs))
  i_fa <- which(grepl("^fastapath=", lines))
  if (length(i_fa) == 1) lines[i_fa] <- paste0("fastapath=", fasta_path)
  i_out <- which(grepl("^outputpath=", lines))
  if (length(i_out) == 1) lines[i_out] <- paste0("outputpath=", out_dir_auto, "\\")
  i_name <- which(grepl("^outputname=", lines))
  if (length(i_name) == 1) lines[i_name] <- paste0("outputname=", outputname)
  # 线程数写死为总核数-1（覆盖模板 cfg 里的 thread 值）：留一核给系统与交互
  # （GUI/网页壳子/其他程序），pFind 搜库期间机器保持可用、不会满负荷卡死
  # threads are pinned to total cores - 1 (overriding the template cfg's thread
  # value): one core stays free for the OS and interactivity (GUI / web shell /
  # other programs), so the machine stays usable during the pFind search
  thread_n <- max(parallel::detectCores() - 1L, 1L)
  i_thr <- which(grepl("^thread=", lines))
  if (length(i_thr) == 1) {
    lines[i_thr] <- paste0("thread=", thread_n)
  } else {
    lines <- c(lines, paste0("thread=", thread_n))
  }
  # 修饰表一律改写为随软件分发的项目内文件（消除模板中对安装目录的硬编码路径）：
  # 开放式搜库用 pUniFind 兼容修饰库（生产 open 结果提取的 302 名修饰构建）——
  # 从搜库侧保证结果只含 pUniFind rescore 认识的修饰名（如 Methyl[D] 类问题在源头
  # 消除），搜库后的兼容过滤保留为保险丝（pFind 仍可能动态生成 Delta_H(...) 等库外
  # 名）；限定式搜库（selectmod 逐名指定）用项目内 pFind_modification.ini（随软件分发
  # 的完整修饰表副本，Carbofuran[S] 已禁用）。
  # The modpath is always rewritten to a project-local file shipped with the client
  # (removing the template's hardcoded installation path): the open search uses the
  # pUniFind-compatible library (built from the 302 production-verified names), so
  # the search itself only reports names pUniFind rescore knows (the Methyl[D]-class
  # problem is removed at the source) with the post-search compat filter kept as a
  # safety net (pFind can still synthesize dynamic Delta_H(...) names); the restricted
  # search (mods named in selectmod) uses the shipped pFind_modification.ini - a copy
  # of the full pFind modification table with Carbofuran[S] disabled.
  mod_ini_name <- if (tolower(outputname) == "pfind_open")
    "pUniFind_mods.ini" else "pFind_modification.ini"
  compat_ini <- file.path(code_dir, mod_ini_name)
  i_mod <- which(grepl("^modpath=", lines))
  if (length(i_mod) == 1 && file.exists(compat_ini)) {
    lines[i_mod] <- paste0("modpath=", compat_ini)
    cat("modpath → 项目内修饰库 ", mod_ini_name, "：", compat_ini, "\n")
    cat("modpath -> project-local mods ", mod_ini_name, ": ", compat_ini, "\n")
  } else if (length(i_mod) == 1) {
    stop("项目内修饰库缺失：", compat_ini)
  }
  # 酶切类型（v1.0，2026-09-26）：按 config_step1 键 enzyme 改写 pFind cfg 的
  # enzyme 行（格式 enzyme=<名> <切点> <禁点> <C/N>，无切点/禁点用 _）；none
  # 保持模板行不变（pFind 无酶切语法随版本而异，留模板）。
  # Enzyme type (v1.0, 2026-09-26): rewrite the pFind cfg enzyme line from
  # the config_step1 key enzyme (format enzyme=<name> <cut> <nocut> <C/N>,
  # "_" for empty cut/nocut); none keeps the template line (pFind's no-enzyme
  # syntax varies across versions).
  en_v <- base_cfg$enzyme
  # v1.0（2026-09-26，审查 P1-3）：键缺失 = 模板即真相——不再无条件按
  # trypsin 改写，否则会覆盖用户手写进模板的酶行（工作站 R29 曾预置
  # LysargiNase 模板行，键缺失时被改写成 Trypsin）。键缺失且模板酶行非
  # Trypsin 时直接报错要求显式填键（模板搜 LysargiNase 而内部按 trypsin
  # 计数会造成口径错位）。
  en_explicit <- !is.null(en_v) && length(en_v) > 0 && nzchar(trimws(en_v[1]))
  i_en <- which(grepl("^enzyme=", lines))
  if (en_explicit) {
    # B-2（2026-09-26，GLM 复核）：未知酶值必须报错而非静默回退 trypsin——
    # 旧写法 suppressWarnings 消化了 enzyme_resolve 的警告、把已解析的
    # "trypsin" 传下去，报错分支不可达
    w_en <- character(0)
    en <- withCallingHandlers(enzyme_resolve(en_v[1]),
                              warning = function(ww) {
                                w_en <<- c(w_en, conditionMessage(ww)); invokeRestart("muffleWarning")
                              })
    if (length(w_en) > 0)
      stop("未知酶切类型：", trimws(en_v[1]), "
unknown enzyme type: ", trimws(en_v[1]))
    if (length(i_en) == 1) {
      lines[i_en] <- pfind_enzyme_line(en)
      cat("pFind enzyme → ", lines[i_en], "（来源 config_step1 键 enzyme）\n")
    }
  } else if (length(i_en) == 1) {
    tpl_name <- trimws(sub("^enzyme=\\s*([^ ]+).*", "\\1", lines[i_en]))
    tpl_key <- gsub("[^a-z0-9]", "", tolower(tpl_name))
    if (!tpl_key %in% c("trypsin", "trypsinp")) {
      stop("config_step1 缺少 enzyme 键，且 pFind 模板酶定义非 Trypsin（", tpl_name,
           "）。请显式填写 enzyme 键（如 enzyme\tLysargiNase）。\n",
           "config_step1 has no enzyme key, and the pFind template enzyme is not Trypsin (",
           tpl_name, "); fill the enzyme key explicitly.")
    }
    cat("pFind enzyme → 保持模板（config_step1 无 enzyme 键，模板为 Trypsin）\n")
  }
  # 以 GBK 连接写出自动生成的运行 cfg，保证文件单一编码（pFind 可读）
  # Write the generated run cfg through a GBK connection, keeping a single file encoding (pFind-readable)
  dir.create(out_dir_auto, recursive = TRUE, showWarnings = FALSE)
  cfg_run <- file.path(dirname(cfg_tpl), paste0("auto_", outputname, ".cfg"))
  con <- file(cfg_run, open = "w", encoding = "GBK")
  writeLines(lines, con)
  close(con)
  # ---- P4（2026-09-27 用户口径：全部搜库流程改逐个 raw/mgf 单独搜库、并行搜库）----
  # pfind_workers（config_step1 键，默认 1=旧行为）>1 时：按 raw 分组，每组一份
  # cfg（msmspath 只含该组 mgf、outputname=<name>_w<k>、线程按 worker 数摊薄），
  # PSOCK 并行跑 Searcher，完成后把各组产物按 raw 顺序合并回 <outputname>.spectra
  # 与 <outputname>-Filtered.spectra（与旧单 cfg 产物等价，下游合并不变）。
  pw <- suppressWarnings(as.integer(tail(c(base_cfg$pfind_workers, "1"), 1)))
  if (is.na(pw) || pw < 1L) pw <- 1L
  if (pw > 1L && length(mgfs) > 1L) {
    groups <- split(seq_along(mgfs), rep(seq_len(pw), length.out = length(mgfs)))
    thread_w <- max(floor(thread_n / pw), 1L)
    cfg_files <- character(length(groups))
    i_ms0 <- which(grepl("^msmspath", lines))
    i_num0 <- which(grepl("^msmsnum=", lines))
    i_name0 <- which(grepl("^outputname=", lines))
    i_thr0 <- which(grepl("^thread=", lines))
    for (w in seq_along(groups)) {
      g <- groups[[w]]
      lw <- lines[!grepl("^msmspath", lines)]
      lw <- append(lw, paste0("msmspath", seq_along(g), "=", mgfs[g]),
                   after = min(i_ms0) - 1L)
      if (length(i_num0) == 1) {
        lw[i_num0] <- paste0("msmsnum=", length(g))
      } else {
        lw <- c(lw, paste0("msmsnum=", length(g)))
      }
      if (length(i_name0) == 1) lw[i_name0] <- paste0("outputname=", outputname, "_w", w)
      if (length(i_thr0) == 1) lw[i_thr0] <- paste0("thread=", thread_w)
      cfg_w <- file.path(dirname(cfg_tpl), sprintf("auto_%s_w%d.cfg", outputname, w))
      con_w <- file(cfg_w, open = "w", encoding = "GBK")
      writeLines(lw, con_w)
      close(con_w)
      cfg_files[w] <- cfg_w
    }
    pfind_bin <- dirname(pfind_exe)
    log_base <- file.path(out_dir_auto, outputname)
    cl <- parallel::makeCluster(length(cfg_files))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(cl, c("pfind_bin", "pfind_exe", "log_base"),
                            envir = environment())
    sts <- parallel::parLapply(cl, cfg_files, function(cf) {
      old <- setwd(pfind_bin)
      on.exit(setwd(old), add = TRUE)
      system2(pfind_exe, shQuote(cf),
              stdout = paste0(log_base, "_", basename(cf), ".log"),
              stderr = paste0(log_base, "_", basename(cf), ".log"))
    })
    cat("运行 pFind（并行 ", length(groups), " worker，每组线程 ", thread_w,
        "）：", pfind_exe, " auto_", outputname, "_w*.cfg\n", sep = "")
    cat("running pFind (parallel ", length(groups), " workers, ", thread_w,
        " threads each): ", pfind_exe, " auto_", outputname, "_w*.cfg\n", sep = "")
    if (any(unlist(sts) != 0)) stop("pFind auto search failed (parallel): ", outputname)
    # 按组序合并各组 Filtered 与 .spectra（首件带头、其余去头追加；组序=raw 序）
    flt_files <- unlist(lapply(seq_along(groups), function(w) {
      ff <- file.path(out_dir_auto, paste0(outputname, "_w", w, "-Filtered.spectra"))
      if (file.exists(ff)) ff else character(0)
    }))
    if (length(flt_files) == 0)
      stop("pFind auto search failed (parallel): 无 Filtered 产物 ", outputname)
    heads <- lapply(flt_files, readLines, n = 1, warn = FALSE)
    bodies <- lapply(flt_files, function(p) readLines(p, warn = FALSE)[-1])
    writeLines(c(heads[[1]], unlist(bodies, use.names = FALSE)),
               file.path(out_dir_auto, paste0(outputname, "-Filtered.spectra")),
               useBytes = TRUE)
    sp_files <- unlist(lapply(seq_along(groups), function(w) {
      ff <- file.path(out_dir_auto, paste0(outputname, "_w", w, ".spectra"))
      if (file.exists(ff)) ff else character(0)
    }))
    if (length(sp_files) > 0) {
      sh <- readLines(sp_files[1], n = 1, warn = FALSE)
      sb <- unlist(lapply(sp_files, function(p) readLines(p, warn = FALSE)[-1]),
                   use.names = FALSE)
      writeLines(c(sh, sb), file.path(out_dir_auto, paste0(outputname, ".spectra")),
                 useBytes = TRUE)
    }
  } else {
    # 调用 pFind：在其 bin 目录运行（保证 ini 相对路径解析正确），失败则中止
    # Invoke pFind: run from its bin folder (correct ini relative-path resolution), aborting on failure
    cat("运行 pFind：", pfind_exe, cfg_run, "（逐 raw 搜库：", length(mgfs), " 个 mgf 文件，线程 ", thread_n, "）\n")
    cat("running pFind: ", pfind_exe, cfg_run, " (per-raw search: ", length(mgfs), " mgf files, ", thread_n, " threads)\n")
    # v1.0：pFind 孤儿实例防撞车——上轮运行被外杀（如排障 taskkill 连带
    # R 主控）后，Searcher 作为孤儿继续搜库；此时重跑再启第二个实例会被
    # pFind 单实例语义顶掉（秒退非 0）误报"搜库失败"（Chlamy_allengines_test
    # 实弹：11:08 孤儿跑到 F6，13:05 重跑秒败）。改为检测到活实例就等待＋复用：
    # 结束后有 .spectra 直接用；无产物才真正启动（pFind 按已完成 .qry.res 分块
    # 续跑）。探活走 shell 整条命令（v1.0 tasklist 拆参坑同口径）。
    # v1.0: orphan-instance guard - when a previous master was killed
    # externally (e.g. a debugging taskkill), the Searcher it spawned keeps
    # running as an orphan; a rerun would launch a second instance that pFind's
    # single-instance semantics rejects (instant nonzero exit), misreported as
    # "search failed" (live fire: an 11:08 orphan reached chunk F6 while the
    # 13:05 rerun died at once). Now wait for the live instance and reuse:
    # if .spectra appears after it exits, use it; launch only when there is
    # nothing to reuse (pFind resumes finished chunks from .qry.res files).
    # Liveness probe via a whole-command shell call (same pattern as the
    # v1.0 tasklist argument-splitting fix).
    pfind_alive <- function() {
      o <- suppressWarnings(shell('tasklist /NH | findstr /I /C:"earcher.exe"',
                                  intern = TRUE, mustWork = FALSE))
      length(o) > 0 && any(grepl("earcher", o))
    }
    if (pfind_alive()) {
      cat("[pFind] 检测到既有 Searcher/PSSearcher 实例运行中，等待其完成后复用（每分钟心跳，最长 12h）\n")
      cat("[pFind] a live Searcher/PSSearcher instance is running; waiting to reuse its output (heartbeat per minute, max 12h)\n")
      waited <- 0L
      while (pfind_alive()) {
        Sys.sleep(60); waited <- waited + 1L
        if (waited %% 10L == 0L)
          cat("[pFind] 仍在等待既有实例（", waited, " 分钟）\n", sep = "")
        if (waited > 720L)
          stop("pFind 既有实例超过 12h 未结束，放弃等待——请手动核查 Searcher/PSSearcher 进程\n",
               "the live pFind instance did not finish within 12h; inspect Searcher/PSSearcher manually")
      }
    }
    cand_pre <- list.files(out_dir_auto, pattern = "\\.spectra$")
    if (length(cand_pre) > 0) {
      cat("[pFind] 结果已存在（含既有实例刚完成），跳过搜库：", out_dir_auto, "\n")
      cat("[pFind] result already present (incl. just finished by the live instance), skipping the search: ", out_dir_auto, "\n")
    } else {
      # pFind 相对自身工作目录解析 aa.ini / modification.ini，故在其 bin 目录运行（同用户习惯：cd bin && Searcher.exe cfg）
      # pFind resolves aa.ini / modification.ini relative to its CWD, so run it from
      # its own bin folder (as the user would: cd bin && Searcher.exe cfg)
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(dirname(pfind_exe))
      st <- system2(pfind_exe, shQuote(cfg_run),
                    stdout = file.path(out_dir_auto, paste0(outputname, ".log")),
                    stderr = file.path(out_dir_auto, paste0(outputname, ".log")))
      setwd(old_wd)
      if (st != 0)
        stop("pFind auto search failed (exit ", st, "): ", outputname,
             "\n日志 log: ", file.path(out_dir_auto, paste0(outputname, ".log")),
             "\n若系统里另有 Searcher/PSSearcher 在跑，pFind 单实例会立即退出——等它结束再重跑",
             "\nif another Searcher/PSSearcher is running, pFind exits at once; wait for it and rerun")
    }
  }
  # 开放式搜库以全量 modification.ini 为修饰表时会解析出 pUniFind rescore 不认识的修饰名
  # （如 Methyl[D]，导致 KeyError）；对 open 项目按兼容集合（生产 open 结果提取）过滤，
  # 删除含集合外修饰名的行——qry.res 的修饰名经查均为兼容集合内，无需处理。
  # The open search over the full modification.ini reports modification names that
  # pUniFind rescore does not know (e.g. Methyl[D], a KeyError); for the open project
  # the results are filtered against the compatible set (built from the production
  # open results), dropping rows with outside names - the qry.res names were checked
  # to lie inside the set, so they need no treatment.
  if (tolower(outputname) == "pfind_open") {
    compat_file <- file.path(code_dir, "pUniFind_compat_mods.txt")
    # 复用检测（2026-09-19 优化）：过滤结果随搜库结果一起复用，已完成过滤的
    # 项目打标记，重跑时跳过重过滤（原实现每次重跑都重新扫整份 pFind.spectra）
    # reuse detection (opt 2026-09-19): the filtered output is reused together
    # with the search result; a marker is written after filtering and later runs
    # skip re-filtering (the old code re-scanned the whole pFind.spectra on every
    # rerun)
    compat_done <- file.path(out_dir_auto, ".compat_filter_done")
    if (file.exists(compat_done)) {
      cat("pUniFind 兼容修饰过滤：已完成，直接复用 ", out_dir_auto, "\n")
      cat("pUniFind compat-mod filter: already done, reused ", out_dir_auto, "\n")
    } else if (file.exists(compat_file)) {
      for (sf in c("pFind.spectra", "pFind-Filtered.spectra")) {
        sp <- file.path(out_dir_auto, sf)
        if (file.exists(sp)) {
          r <- filter_mods_compatible(sp, compat_file)
          cat("pUniFind 兼容修饰过滤：", sf,
              " 行 ", r["rows"], "，保留 ", r["kept"], "，剔除 ", r["dropped"], "\n")
          cat("pUniFind compat-mod filter: ", sf,
              " rows ", r["rows"], " kept ", r["kept"], " dropped ", r["dropped"], "\n")
        }
      }
      writeLines("done", compat_done)
    } else {
      cat("注意：未找到 pUniFind 兼容修饰集合 ", compat_file, "\n")
      cat("note: pUniFind compat mods set missing ", compat_file, "\n")
    }
  }
  # 返回最终 .spectra 结果（优先 Filtered）；pFind 退出码 0 但无产出时报错，
  # 避免该策略被静默丢弃
  # Return the final .spectra result (Filtered preferred); error out when pFind
  # exited 0 but produced nothing, so the strategy is never silently dropped
  cand <- list.files(out_dir_auto, pattern = "\\.spectra$", full.names = TRUE)
  flt <- cand[grepl("Filtered", basename(cand))]
  if (length(cand) > 0) writeLines(fp_lines, rp_file, useBytes = TRUE)   # P2-2 参数指纹留档
  if (length(flt) > 0) return(tail(flt, 1))
  if (length(cand) > 0) return(tail(cand, 1))
  stop("pFind 搜库完成但无 .spectra 结果：", out_dir_auto,
       "\npFind finished but no .spectra in: ", out_dir_auto)
}

# 从第二步/基础配置定位 pFind 可执行文件（pfind_exe 键），取最后一个非空值并规范化
# Locate the pFind executable from the step-2/base configs (pfind_exe key), keep the last non-empty value and normalize it
pfind_dir <- NULL
for (src in list(cfg2, base_cfg)) {
  if (!is.null(src$pfind_exe) && length(src$pfind_exe) > 0 &&
      tail(src$pfind_exe, 1) != "") {
    pfind_dir <- tail(src$pfind_exe, 1)
    break
  }
}
pfind_exe <- if (!is.null(pfind_dir)) {
  try(normalize_pfind_exe(pfind_dir), silent = TRUE)
} else ""
if (inherits(pfind_exe, "try-error")) pfind_exe <- ""
# 判断是否启用自动搜库：任一配置的 auto_search 取值为 "1" 即启用
# Decide whether auto search is on: enabled when any config's auto_search equals "1"
auto_on <- FALSE
for (src in list(cfg2, base_cfg)) {
  if (!is.null(src$auto_search) && length(src$auto_search) > 0 &&
      tail(src$auto_search, 1) == "1") auto_on <- TRUE
}
  # 第二步起不再输出任何 mgf：逐 raw 纯文件已由第一步在 unlabeled/ 下产出（按标题 raw 前缀重切），
  # 本步只读不写；pUniFind GUI 的 mgf 目录直接选用 unlabeled/
  # From step 2 on no mgf is ever written: the pure per-raw files were produced by
  # step 1 under unlabeled/ (re-split by title raw prefix); this step only reads
  # them, and the pUniFind GUI uses unlabeled/ directly as its mgf folder
  # 自动搜库分支：用户未提供搜库结果、开启自动搜库且 pFind 可用时执行
  # （2026-09-13 修正：自动搜库的 b/c 与用户提供的外部引擎条目是并列策略，
  # 用户提供了 search_file 时同样要跑——此前互斥导致 b/c 缺位）
  # Auto-search branch: runs when auto search is on and pFind is available
  # (2026-09-13 fix: the auto b/c are parallel strategies alongside the
  # user-provided engine files and must run even when search_file is given -
  # the old mutual exclusion silently dropped b/c)
  auto_files <- character(0)
  auto_types <- character(0)
  if (auto_on && nzchar(pfind_exe)) {
  mgf_dir_chk <- file.path(cfg2$config_dir, "unlabeled")
  if (length(list.files(mgf_dir_chk, pattern = "\\.mgf$")) == 0) {
    stop("unlabeled 逐 raw MGF 缺失：", mgf_dir_chk, "\nper-raw MGFs missing in: ", mgf_dir_chk)
  }
  # 依次处理开放式（open）与限定式（restrict）两个模板，收集自动搜库产生的 .spectra；
  # 约定目录模板缺失或 config_step2 缺键时从内置模板自动生成（无需停点手填）
  # Loop over the open/restrict templates and collect the .spectra produced by the auto
  # search; a missing template (convention folder or config key) is auto-generated from
  # the built-in copy (no stop point for manual filling)
  # v1.0（工作站反馈#19）：qry.res 读取段默认顺序（实测比 6 路 PSOCK 快
  # 约 2.5 倍且省内存）；config_step1 的 pfind_reader_workers >=1 时按该值并行。
  rdr_w <- cfg_get_preferred(base_cfg, cfg2, "pfind_reader_workers", "")$value
  rdr_w <- if (length(rdr_w) > 0 && nzchar(trimws(as.character(rdr_w[1]))))
    suppressWarnings(as.integer(rdr_w[1])) else NA_integer_
  if (is.na(rdr_w) || rdr_w < 1L) { rdr_par <- FALSE; rdr_n <- NULL }
  else { rdr_par <- TRUE; rdr_n <- rdr_w }

  # ---- v1.0 两队列调度器（2026-10-02，工作站设计稿＋开发侧评审四修复）----
  # 实现＝两条物理串行道：master 前台跑 CPU 道（open→open-pp→restrict→
  # restrict-pp→引擎），GPU 管家子进程独占 GPU 道（denovo→[等 open-pp 标记]
  # →rescore）。道间并行等价于提案的补位轨迹：denovo 完即转 rescore（GPU 无
  # 空窗）、rescore ∥ restrict/引擎；亲和与跨队列唤醒由物理道结构消解。
  # 评审修复：③ rescore 前置含 open 后处理（管家等 _tq_open_pp_done.txt）；
  # ④ 引擎槽无假依赖（master 原生顺序即 FIFO，无 pFind 谓词）。
  # two_queue_scheduler=1 启用；默认 0＝与 v1.0 逐行为一致。
  # ---- v1.0 two-queue scheduler (workstation design + dev review fixes) --
  # Realized as two physical serial lanes: the master runs the CPU lane in the
  # foreground (open -> open-pp -> restrict -> restrict-pp -> engines) while a
  # GPU steward child owns the GPU lane (denovo -> [wait for the open-pp marker]
  # -> rescore). Inter-lane overlap equals the proposal's refill trajectory:
  # denovo hands straight to rescore (no GPU idle window) and rescore overlaps
  # restrict/engines. Affinity and cross-queue wake-up dissolve into the
  # two-lane structure. Review fixes: (3) rescore is gated on open post-
  # processing (the steward waits for _tq_open_pp_done.txt); (4) engine slots
  # carry no false pFind dependency (the master's native order is the FIFO).
  # Enable with two_queue_scheduler=1 in config_step1; default 0 = unchanged.
  tq_on <- FALSE
  for (src_tq in list(base_cfg, cfg2)) {
    if (!is.null(src_tq$two_queue_scheduler) && length(src_tq$two_queue_scheduler) > 0 &&
        identical(trimws(as.character(src_tq$two_queue_scheduler[1])), "1")) tq_on <- TRUE
  }
  # 无条件加载：旧 DAG join 也需要 tq_pid_alive（v1.0 修复其 tasklist 拆参 bug）
  # loaded unconditionally: the legacy DAG join also needs tq_pid_alive (v1.0
  # fix for its argument-splitting tasklist bug)
  source(file.path(code_dir, "two_queue_steward.R"), encoding = "UTF-8")
  if (tq_on) {
    cat("\n==== [两队列调度器] GPU 管家启动（CPU 道随后在前台展开）====\n")
    cat("==== [two-queue scheduler] GPU steward launching (CPU lane follows in the foreground) ====\n")
    tq_launch_steward(cfg2, cfg2_path, code_dir)
  }

  # ---- v1.0 DAG 双通道流水线（2026-10-01，用户裁定）------------------
  # 队列 1（CPU，前台）：open-pFind → restrict-pFind → 引擎槽
  # 队列 2（GPU，后台）：denovo（T=0 立即启动，零依赖）→ rescore（等 open 完
  #   成后且 denovo 释放 GPU 后执行）
  # dag_pipeline=1 时启用；默认 0=现行顺序（零行为变化）。
  # ---- v1.0 DAG dual-lane pipeline (user decision 2026-10-01) ---------
  # Lane 1 (CPU, foreground): open-pFind -> restrict-pFind -> engine slots
  # Lane 2 (GPU, background):  denovo (starts at T=0, zero dependency)
  #   -> rescore (waits for open results + denovo GPU release)
  # Enable with dag_pipeline=1 in config_step1; default 0 = sequential.
  dag_on <- FALSE
  for (src_dag in list(base_cfg, cfg2)) {
    if (!is.null(src_dag$dag_pipeline) && length(src_dag$dag_pipeline) > 0 &&
        identical(trimws(as.character(src_dag$dag_pipeline[1])), "1")) dag_on <- TRUE
  }
  denovo_bg_pid <- NULL
  denovo_bg_log <- NULL

  # DAG 前置：从 config 直接读取 pUniFind 设置（不依赖后续定义的变量）
  # DAG pre-read: pull pUniFind settings directly from config (no dependency on later-defined vars)
  if (dag_on) {
    puni_exe_dag <- ""
    puni_py_dag <- ""
    denovo_on_dag <- FALSE
    for (src_dag in list(base_cfg, cfg2)) {
      if (!is.null(src_dag$punifind_exe) && length(src_dag$punifind_exe) > 0 &&
          nzchar(trimws(as.character(src_dag$punifind_exe[1])))) {
        puni_exe_dag <- normalizePath(trimws(as.character(src_dag$punifind_exe[1])), winslash="/")
      }
      if (!is.null(src_dag$punifind_python) && length(src_dag$punifind_python) > 0 &&
          nzchar(trimws(as.character(src_dag$punifind_python[1])))) {
        puni_py_dag <- normalizePath(trimws(as.character(src_dag$punifind_python[1])), winslash="/")
      }
      if (!is.null(src_dag$punifind_denovo) && length(src_dag$punifind_denovo) > 0 &&
          identical(trimws(as.character(src_dag$punifind_denovo[1])), "1")) denovo_on_dag <- TRUE
    }
    # 读取并行旋钮（用默认值兜底）
    pf_cp_dag <- 1L
    for (src_dag in list(base_cfg, cfg2)) {
      v_raw <- if (!is.null(src_dag$punifind_chunk_parallel) &&
                     length(src_dag$punifind_chunk_parallel) > 0)
        as.character(src_dag$punifind_chunk_parallel[1]) else ""
      v <- suppressWarnings(as.integer(trimws(v_raw)))
      if (!is.na(v) && v >= 1) pf_cp_dag <- v
    }
  } else {
    denovo_on_dag <- FALSE; puni_exe_dag <- ""
  }

  if (!tq_on && dag_on && denovo_on_dag && nzchar(puni_exe_dag)) {
    # 检查 denovo 是否已有缓存产物
    # （v1.0：two_queue_scheduler=1 时由 GPU 管家接管，本块跳过）
    # (v1.0: superseded by the GPU steward when two_queue_scheduler=1)
    dproj_dag <- file.path(cfg2$config_dir, "pUniFind_denovo")
    csvs_dag <- list.files(file.path(dproj_dag, "pUniFind_result"),
                           pattern = "\\.csv$", full.names = TRUE)
    if (length(csvs_dag) == 0) {
      cat("\n==== [DAG 流水线] GPU 通道启动：pUniFind denovo 后台运行 ====\n")
      cat("==== [DAG pipeline] GPU lane starting: pUniFind denovo in background ====\n")

      # 清理 stray MGF（与现有 denovo 执行前清理同口径）
      unlab_dag <- file.path(cfg2$config_dir, "unlabeled")
      stray_dag <- list.files(unlab_dag, pattern = "_check\\.mgf$", full.names = TRUE)
      if (length(stray_dag) > 0) {
        stray_dir_dag <- file.path(cfg2$config_dir, "_stray_mgf")
        dir.create(stray_dir_dag, recursive = TRUE, showWarnings = FALSE)
        file.rename(stray_dag, file.path(stray_dir_dag, basename(stray_dag)))
        cat("[DAG denovo] 移出 ", length(stray_dag), " 个 *_check.mgf\n", sep = "")
      }

      # 生成后台启动脚本。v1.0 fix 5＋fix 6（实弹四轮定位）：
      # (1) 脚本内不能出现中文字符量——code_dir 含"正式流程"，无 BOM 的
      #     UTF-8 字节被子进程 R 按本地码页误读 → setwd 失败退出。改为脚本
      #     零 setwd：工作目录由启动器 -WorkingDirectory 指定，getwd() 取。
      # (2) fix 6（2026-10-02 DAG 实弹）：**不能加 BOM**——本机 Rscript 不
      #     跳过脚本首的 UTF-8 BOM（unexpected input in "\ufeff"），加了
      #     必崩。非 ASCII 路径字面量的兜底改为 \uXXXX/\UXXXXXXXX 转义
      #     （R 解析期还原、与文件编码无关），未来项目路径含中文也安全。
      # (3) 启动前清 LC_CTYPE/LC_ALL/LANG——Git Bash 会话的 LC_CTYPE=
      #     C.UTF-8 会被继承，Windows R 设置失败回落 C locale，中文注释
      #     的 source() 直接 invalid input。
      # generate the background launch script. v1.0 fix 5 + fix 6 (four
      # rounds of live debugging): (1) no non-ASCII literals inside the
      # script — the code_dir contains Chinese and a BOM-less UTF-8 literal
      # mojibakes in the child R. Zero setwd instead: the working dir comes
      # from the launcher's -WorkingDirectory and getwd() recovers it.
      # (2) fix 6 (live DAG test 2026-10-02): do NOT add a BOM — Rscript here
      # does not skip a leading UTF-8 BOM ("unexpected input in ..."), the
      # BOM itself crashes the child. The fallback for non-ASCII path
      # literals is \uXXXX/\UXXXXXXXX escaping (restored at R parse time,
      # independent of file encoding), so future Chinese project paths stay
      # safe. (3) LC_CTYPE/LC_ALL/LANG are cleared before launching — an
      # inherited LC_CTYPE=C.UTF-8 from a Git Bash session makes Windows R
      # fall back to the C locale, where source()-ing files with Chinese
      # comments fails with invalid input.
      esc_nonascii <- function(s) {
        s <- enc2utf8(as.character(s))
        paste(vapply(utf8ToInt(s), function(cp)
          if (cp < 128L) intToUtf8(cp)
          else if (cp <= 0xFFFFL) sprintf("\\u%04X", cp)
          else sprintf("\\U%08X", cp),
          character(1)), collapse = "")
      }
      dag_script <- file.path(cfg2$config_dir, "_dag_denovo_bg.R")
      writeLines(c(
        "code_dir <- getwd()",
        "source(file.path(code_dir, 'punifind_gpu_compat.R'), encoding='UTF-8')",
        sprintf("pf_run('denovo', '%s', '%s', '%s', '%s',",
                esc_nonascii(puni_exe_dag),
                esc_nonascii(file.path(dproj_dag, "result")),
                esc_nonascii(unlab_dag), esc_nonascii(dproj_dag)),
        sprintf("  python_exe = '%s', num_workers = 0, num_proc = 8,",
                esc_nonascii(puni_py_dag)),
        sprintf("  batch_size = 64, chunk_parallel = %d, chunk_cleanup = FALSE)",
                pf_cp_dag)
      ), dag_script, useBytes = TRUE)

      # v1.0 fix 4：启动前必须先建目录——重定向目标目录不存在时
      # cmd/PowerShell 都会立即失败（实弹：PID=0 + "The system cannot
      # find the path specified"，denovo 从未启动）。pf_run 在子进程内
      # 才建目录，来不及。
      # v1.0 fix 4: create the target dir BEFORE launching — a missing
      # redirect dir makes cmd/PowerShell fail instantly (live evidence:
      # PID=0 + "The system cannot find the path specified"). pf_run only
      # creates dirs inside the child, which is too late.
      dir.create(dproj_dag, recursive = TRUE, showWarnings = FALSE)
      # v1.0 fix 3：路径含中文字符（isoPSV正式流程），cmd.exe .bat 无法处理。
      # 改用 PowerShell Start-Process（原生 Unicode 路径支持），取真实 PID。
      # v1.0 fix 3: path contains Chinese chars (isoPSV正式流程) which
      # cmd.exe .bat cannot handle. Use PowerShell Start-Process (native
      # Unicode path support) and capture the real child PID.
      rscript_full <- file.path(R.home("bin"), "Rscript")
      # Start-Process 不允许 stdout/stderr 指向同一文件，须分离
      # Start-Process forbids the same file for stdout and stderr; keep them apart
      denovo_bg_out <- file.path(dproj_dag, "denovo_bg.out.log")
      denovo_bg_err <- file.path(dproj_dag, "denovo_bg.err.log")
      # -WorkingDirectory 让子进程 getwd() 落在 code_dir（路径含中文也无
      # 妨——命令行经 CreateProcessW 传 Unicode，不经文件解析）
      # -WorkingDirectory makes getwd() land in code_dir inside the child
      # (Chinese in the path is fine — the command line travels as Unicode
      # via CreateProcessW, never parsed from a file)
      ps_cmd <- sprintf(
        paste0('$env:LC_CTYPE=$null; $env:LC_ALL=$null; $env:LANG=$null; ',
               'Start-Process -FilePath "%s" -ArgumentList \'"%s"\' ',
               '-WorkingDirectory "%s" -WindowStyle Hidden ',
               '-RedirectStandardOutput "%s" -RedirectStandardError "%s" ',
               '-PassThru | Select-Object -ExpandProperty Id'),
        rscript_full, dag_script, code_dir, denovo_bg_out, denovo_bg_err)
      denovo_bg_pid <- suppressWarnings(as.integer(
        system2("powershell", c("-NoProfile", "-Command", shQuote(ps_cmd)),
                stdout = TRUE, stderr = FALSE, invisible = TRUE)))
      if (is.na(denovo_bg_pid)) denovo_bg_pid <- -1L
      cat("[DAG 流水线] denovo 后台 PID: ", denovo_bg_pid,
          "  日志: ", denovo_bg_out, "\n", sep = "")
      cat("[DAG pipeline] denovo background PID: ", denovo_bg_pid, "\n", sep = "")
    } else {
      cat("[DAG 流水线] denovo 已有缓存产物（", length(csvs_dag), " 个 csv），跳过后台启动\n", sep = "")
      cat("[DAG pipeline] denovo cache found (", length(csvs_dag), " csvs), skipping background launch\n", sep = "")
    }
  }

  for (nm in c("pfind_cfg_open", "pfind_cfg_restrict")) {
    tpl <- cfg2[[nm]]
    if (is.null(tpl) || length(tpl) == 0 || tpl[1] == "") {
      kind <- if (nm == "pfind_cfg_open") "open" else "restrict"
      tpl <- ensure_pfind_template(kind, fasta_path, cfg2$config_dir)
    } else {
      tpl <- tail(resolve_cfg_paths(cfg2, nm), 1)
    }
    tag <- if (nm == "pfind_cfg_open") "open" else "restrict"
    cat("==== [自动搜库] pFind ", tag, " ====\n")
    cat("==== [Auto search] pFind ", tag, " ====\n")
    sp_auto <- run_pfind_auto(tpl)
    # pFind 结果后处理（用户口径 2026-09-16）：未过滤 pFind.spectra 全量替换
    # Filtered 参与竞争；其中 decoy 行用 qry.res 中最高分非 decoy 候选救援；
    # 处理结果写回 pFind-Filtered.spectra（官方原版备份 .origin），后续合并、
    # 理论谱图生成与竞争打分自动消费新表
    # pFind post-processing (user ruling 2026-09-16): the unfiltered
    # pFind.spectra replaces Filtered for the competition; decoy rows are
    # rescued with the best non-decoy candidate from qry.res; the result is
    # written back as pFind-Filtered.spectra (official copy backed up as
    # .origin), which the merge / spectrum prediction / competition scoring
    # consume automatically
    pp_dir <- dirname(sp_auto)
    pp_mark_file <- file.path(pp_dir, ".pfind_postprocess_done")
    if (pfind_postprocess_fresh(pp_dir, fasta_path, pp_mark_file)) {
      cat("[pFind 后处理] 指纹一致，复用既有后处理产物（跳过约 3h 重算）：", pp_dir, "\n")
      cat("[pFind postprocess] fingerprint match, reusing output (skips ~3h): ", pp_dir, "\n")
    } else {
      # v1.0（工作站反馈#18/#21）：复用判据＋失败醒目标记。失败不再静默
      # "跳过"——写 .pfind_postprocess_failed 并醒目提示口径降级（合并继续用
      # 未救援的 pFind-Filtered.spectra）。
      pp_ok <- tryCatch({
        pfind_res_postprocess(pp_dir, fasta_path,
                              reader_parallel = rdr_par, reader_workers = rdr_n)
        TRUE
      }, error = function(e) {
        cat("!! pFind 后处理失败：", conditionMessage(e), "\n",
            "!! 合并将继续使用未救援的 pFind-Filtered.spectra（口径降级，\n",
            "!! 详见标记文件 .pfind_postprocess_failed）\n", sep = "")
        cat("!! pFind postprocess FAILED: ", conditionMessage(e), "\n",
            "!! the merge will proceed with the un-rescued pFind-Filtered (degraded;\n",
            "!! see marker .pfind_postprocess_failed)\n", sep = "")
        writeLines(c(format(Sys.time()), conditionMessage(e)),
                   file.path(pp_dir, ".pfind_postprocess_failed"), useBytes = TRUE)
        FALSE
      })
      if (pp_ok) pfind_postprocess_mark(pp_dir, fasta_path, pp_mark_file)
    }
    auto_files <- c(auto_files, sp_auto)
    auto_types <- c(auto_types, if (nm == "pfind_cfg_open") "pfind_open"
                    else "pfind_restrict")
    # v1.0 两队列：open 搜库＋后处理完成后写标记，GPU 管家据此解锁 rescore
    # （评审修复③：rescore 前置含后处理；标记幂等，off 模式不产生）
    # v1.0 two-queue: after the open search + post-processing complete, drop
    # the marker that unlocks the steward's rescore phase (review fix 3; the
    # marker is idempotent and never written when the scheduler is off)
    if (tq_on && identical(tag, "open"))
      writeLines(format(Sys.time()),
                 file.path(cfg2$config_dir, "_tq_open_pp_done.txt"), useBytes = TRUE)
    pt_mark(paste0("pFind 自动搜库 (", tag, ")"),
            paste0("pFind auto search (", tag, ")"))
  }

  # ---- v1.0 两队列汇合门：等 GPU 管家终态（denovo＋rescore 都完成）----
  # 放在 pFind 前台链（open/restrict）之后：rescore 与 restrict 的并行由管家
  # 自行完成（denovo 完即转 rescore），master 到这里时通常只剩 rescore 尾段。
  # 管家活着则持续等待（无自动放行，防合并/rescore 写盘竞态）；死而无终态按
  # v1.0 失败容忍口径放行，下方发现段自会处理缺失产物。
  # ---- v1.0 two-queue join gate: wait for the GPU steward terminal state --
  # Placed after the foreground pFind chain: the rescore/restrict overlap is
  # the steward's own doing (denovo hands straight to rescore), so at most the
  # rescore tail remains when the master arrives. A live steward is waited on
  # with no auto-bail (no merge-vs-rescore write races); a dead steward without
  # a terminal state proceeds per the v1.0 failure-tolerance semantics and
  # the discovery sections below tolerate the missing outputs.
  if (tq_on) {
    cat("\n==== [两队列调度器] CPU 前台链完成，等 GPU 管家汇合 ====\n")
    cat("==== [two-queue scheduler] CPU foreground chain done, joining the GPU steward ====\n")
    tq_join_gate(cfg2)
  }

  # ---- v1.0 DAG 汇合：等后台 denovo 完成 ----
  # pFind open/restrict（CPU 通道）已在前台跑完；现在等 GPU 通道的 denovo。
  # 等完后 rescore 在下方现有代码中执行（此时 GPU 已释放，open 结果已就绪）。
  # ---- v1.0 DAG join: wait for background denovo to complete ----
  # The CPU lane (open + restrict) has finished in the foreground; now wait
  # for the GPU lane (denovo). After this join, rescore runs via the existing
  # code below (GPU is free, open results are ready).
  if (!tq_on && !is.null(denovo_bg_pid) && denovo_bg_pid > 0) {
    cat("\n==== [DAG 流水线] CPU 通道完成，等待 GPU 通道（denovo PID=", denovo_bg_pid, "）====\n", sep = "")
    cat("==== [DAG pipeline] CPU lane done, waiting for GPU lane (denovo PID=", denovo_bg_pid, ") ====\n", sep = "")
    dag_elapsed <- 0
    repeat {
      # v1.0：system2("tasklist", c("/FI","PID eq N")) 把 "PID eq N" 按空格
      # 拆成多参 → tasklist 报 Invalid argument 'eq' 恒判死（此前被 denovo 早死
      # 掩盖，two_queue 实弹暴露）。改用 tq_pid_alive 整条命令探测。
      # v1.0: system2("tasklist", c("/FI", "PID eq N")) splits the filter at
      # spaces and tasklist always errors out (masked until now by denovo dying
      # early; exposed by the two-queue live tests). Probe via tq_pid_alive.
      if (!tq_pid_alive(denovo_bg_pid)) break
      Sys.sleep(30)
      dag_elapsed <- dag_elapsed + 30
      if (dag_elapsed %% 300 == 0) {
        cat("[DAG 流水线] denovo 后台运行中... 已等 ", dag_elapsed / 60, " min\n", sep = "")
        cat("[DAG pipeline] denovo background running... waited ", dag_elapsed / 60, " min\n", sep = "")
      }
      if (dag_elapsed > 36000) { # 10h 超时保护
        warning("DAG: denovo background process did not exit within 10h")
        break
      }
    }
    cat("[DAG 流水线] GPU 通道（denovo）完成，耗时 ", dag_elapsed / 60, " min\n", sep = "")
    cat("[DAG pipeline] GPU lane (denovo) done after ", dag_elapsed / 60, " min\n", sep = "")
    pt_mark("pUniFind 自动 denovo (DAG background)", "pUniFind auto denovo (DAG background)")
  }

  # ---- pUniFind 两项：策略 4 rescore（对 pFind 开放式结果重打分）、策略 5 denovo（从头测序） ----
  # 仅在用户 config1 选择 use_punifind=1 时启用；0＝全 pFind 流程（忽略 pUniFind 输入）。
  # pUniFind 由用户在 GUI（piUniMSGUI.exe）中手工完成；软件负责项目准备与产出自动发现：
  #   rescore 项目＝pFind_open（.pac 复制到项目根；mgf 目录在 GUI 中选 unlabeled/——第一步产出的
  #   逐 raw 纯文件；rescore 只支持逐 raw 搜库格式，res 文件数须与 raw 数一致）；
  #   denovo 项目＝pUniFind_denovo（mgf 目录同样选 unlabeled/）。第二步起软件不再输出任何 mgf。
  #   GUI 产出下次运行自动登记：rescore <项目名>fdr0.01_pUniFind.spectra→策略 4；
  #   denovo _merged.csv→denovo_csv_to_spectra→策略 5。未产出时打印双语 GUI 待办。
  # pUniFind two auto items: strategy 4 rescore (pFind open result rescoring), strategy 5
  # denovo. Enabled only when the user chose use_punifind=1 in config1; 0 = pure-pFind
  # workflow (pUniFind inputs ignored). pUniFind itself is run by the user in the Windows
  # GUI (piUniMSGUI.exe); the software prepares the projects and auto-discovers the outputs:
  #   rescore project = pFind_open (.pac copied to the project root; select unlabeled/ as
  #   the mgf folder in the GUI - the pure per-raw files from step 1; rescore only supports
  #   the per-raw layout); denovo project = pUniFind_denovo (mgf folder = unlabeled/ too).
  #   From step 2 on the software writes no mgf files at all.
  #   GUI outputs are auto-registered next run: rescore <project>fdr0.01_pUniFind.spectra ->
  #   strategy 4; denovo _merged.csv -> denovo_csv_to_spectra -> strategy 5. When they are
  #   missing a bilingual GUI to-do is printed.
  # 判断是否启用 pUniFind（策略 4 rescore＋策略 5 denovo）：任一配置 use_punifind=1 即启用
  # Decide whether pUniFind is enabled (strategy 4 rescore + strategy 5 denovo): on when any config has use_punifind=1
  use_puni <- FALSE
  for (src in list(cfg2, base_cfg)) {
    if (!is.null(src$use_punifind) && length(src$use_punifind) > 0 &&
        tail(src$use_punifind, 1) == "1") use_puni <- TRUE
  }
  # 读取 pUniFind 各项运行模式（任一配置提供即可）：punifind_exe 为安装的 punifind
  # 目录；punifind_rescore / punifind_denovo 为三值：0=不运行、1=自动运行、
  # 2=手动运行；自动模式（=1）必须提供 punifind_exe，否则报错
  # Read the per-item pUniFind run modes (any config may provide them): punifind_exe
  # is the installed punifind folder; punifind_rescore / punifind_denovo take three
  # values: 0=off, 1=auto, 2=manual; auto mode (=1) requires punifind_exe
  # 键优先级（2026-09-26 功能排队②，deepseek P3-2）：config_step1（base_cfg）
  # 优先于 config_step2（cfg2）——此前 exe/python"先命中即 break"（cfg2 赢）、
  # 模式键"后读覆盖"（step1 赢），两配置共存时来源不对称；统一 step1→step2
  # 并打印生效值与来源。
  # key precedence (2026-09-26, queue item 2, deepseek P3-2): config_step1 wins
  # over config_step2 - previously the exe/python keys took the first hit (cfg2
  # won) while the mode keys were overwritten by the later read (step1 won), an
  # asymmetric source when both configs are present; unified to step1->step2
  # with the effective value and its source printed.
  cfg_get <- function(key, default = NULL, validate = NULL) {
    r <- cfg_get_preferred(base_cfg, cfg2, key, default, validate)
    if (nzchar(r$source)) {
      cat("[pUniFind 配置] ", key, " = ", as.character(r$value),
          "（来源：", r$source, "）\n", sep = "")
    }
    r$value
  }
  puni_exe_raw <- cfg_get("punifind_exe", default = "",
                          validate = function(m) nzchar(trimws(as.character(m))))
  puni_exe <- if (nzchar(puni_exe_raw)) normalizePath(puni_exe_raw, winslash = "/") else ""
  puni_rescore_mode <- cfg_get("punifind_rescore", default = "0",
                               validate = function(m) as.character(m) %in% c("0", "1", "2"))
  puni_denovo_mode  <- cfg_get("punifind_denovo", default = "0",
                               validate = function(m) as.character(m) %in% c("0", "1", "2"))
  auto_rescore   <- puni_rescore_mode == "1"
  manual_rescore <- puni_rescore_mode == "2"
  auto_denovo    <- puni_denovo_mode == "1"
  manual_denovo  <- puni_denovo_mode == "2"
  if ((auto_rescore || auto_denovo) && !nzchar(puni_exe)) {
    stop("punifind_rescore / punifind_denovo 设为 1（自动运行）时必须填写 punifind_exe（pUniFind 安装的 punifind 目录）\n",
         "auto mode (=1) for punifind_rescore / punifind_denovo requires punifind_exe (the installed pUniFind punifind folder)")
  }
  # rescore 打分 pkl→fdr0.01 .spectra 的转换器（2026-09-06）：官方 CLI 只产出
  # 打分 pkl，转换脚本随客户端分发在 code 目录；python 由 punifind_python 配置
  # （需装有 torch/numpy，如 D:/Bioinformatics/torch_cpu_venv/Scripts/python.exe），
  # 未配置时退回 GUI 导出待办口径
  # the rescore pkl -> fdr0.01 .spectra converter (2026-09-06): the official CLI
  # writes only score pkls; the shipped converter script lives in the code dir;
  # python is named by punifind_python (torch/numpy required), and when unset the
  # GUI-export to-do fallback is kept
  puni_py <- ""
  puni_py_raw <- cfg_get("punifind_python", default = "",
                         validate = function(m) nzchar(trimws(as.character(m))))
  if (nzchar(puni_py_raw)) puni_py <- normalizePath(puni_py_raw, winslash = "/")
  converter_py <- file.path(code_dir, "punifind_rescore_convert.py")
  # v1.0（2026-09-26，deepseek GPU 适配模块配套）：加载离线适配模块并解析
  # 数据加载/CPU 进程数（config_step1 键 punifind_num_workers / punifind_num_proc
  # 优先，环境变量 PUNIFIND_NUM_WORKERS / PUNIFIND_NUM_PROC 兜底，默认 8）
  source(file.path(code_dir, "punifind_gpu_compat.R"), encoding = "UTF-8")
  pf_nw_env <- trimws(Sys.getenv("PUNIFIND_NUM_WORKERS", unset = ""))
  pf_nw_cfg <- cfg_get("punifind_num_workers", default = pf_nw_env,
                       validate = function(m) nzchar(trimws(as.character(m))))
  pf_nw <- suppressWarnings(as.integer(pf_nw_cfg))
  if (is.na(pf_nw) || pf_nw < 0) {
    # v1.0（2026-10-01，用户需求：全自动平台判定）：未显式配置时，
    # Windows 默认 0（dataloader 多进程在 Windows spawn 模式下 pickle 报错
    # cannot pickle 'Environment' object），非 Windows 默认 8。
    if (.Platform$OS.type == "windows") {
      pf_nw <- 0L
      cat("[pUniFind 适配] Windows 平台，punifind_num_workers 默认 0（避免多进程 pickle 错误）；如需启用请在 config_step1 设置 punifind_num_workers 8\n")
      cat("[pUniFind compat] Windows detected, defaulting punifind_num_workers to 0 (avoids multiprocessing pickle errors); set punifind_num_workers 8 in config_step1 to override\n")
    } else {
      pf_nw <- 8L
    }
  }
  pf_np_env <- trimws(Sys.getenv("PUNIFIND_NUM_PROC", unset = ""))
  pf_np_cfg <- cfg_get("punifind_num_proc", default = pf_np_env,
                       validate = function(m) nzchar(trimws(as.character(m))))
  pf_np <- suppressWarnings(as.integer(pf_np_cfg))
  if (is.na(pf_np) || pf_np < 1) pf_np <- 8L
  # v1.0（2026-09-30，工作站 A①③/D）：GPU 供给瓶颈三旋钮——推理批大小
  # （punifind_batch_size，默认 64，可试 128–256）、片级并行路数
  # （punifind_chunk_parallel，默认 1＝顺序；16 GB 卡建议 2–3）、合并后可选
  # 清理片内数据（punifind_chunk_cleanup，1＝清理，保留 done 标记与合并产物）
  pf_bs <- suppressWarnings(as.integer(cfg_get("punifind_batch_size", default = "64",
      validate = function(m) nzchar(trimws(as.character(m))))))
  if (is.na(pf_bs) || pf_bs < 1) pf_bs <- 64L
  pf_cp <- suppressWarnings(as.integer(cfg_get("punifind_chunk_parallel", default = "1",
      validate = function(m) nzchar(trimws(as.character(m))))))
  if (is.na(pf_cp) || pf_cp < 1) pf_cp <- 1L
  pf_cu <- identical(trimws(cfg_get("punifind_chunk_cleanup", default = "0",
      validate = function(m) nzchar(trimws(as.character(m))))), "1")
  # pUniFind 项目准备分支：仅启用 pUniFind 时执行（denovo 与 rescore 两个项目）
  # pUniFind project preparation branch: runs only when pUniFind is enabled (denovo +
  # rescore); no mgf is copied here - step 1 already wrote the pure per-raw files under
  # unlabeled/, which is the mgf folder to select in the GUI for both items
  if (use_puni) {
    # rescore 项目准备：.pac 在项目根（官方项目结构 .pac、param/、result/；mgf 目录选 unlabeled/）
    # prepare the rescore project: .pac at the project root (official layout .pac,
    # param/, result/; select unlabeled/ as the mgf folder)
    pfo_tpl <- cfg2$pfind_cfg_open
    pfo_tpl <- if (!is.null(pfo_tpl) && length(pfo_tpl) > 0 && pfo_tpl[1] != "") {
      tail(resolve_cfg_paths(cfg2, "pfind_cfg_open"), 1)
    } else {
      # 缺键时从内置模板自动生成（与自动搜库循环同一生成逻辑）
      # auto-generate from the built-in template when the key is missing
      ensure_pfind_template("open", fasta_path, cfg2$config_dir)
    }
    proj <- if (nzchar(pfo_tpl)) dirname(dirname(pfo_tpl)) else ""
    got4 <- FALSE
    if (nzchar(proj)) {
      # pFind 生成的 td.pac 为追加式命名（<fasta>_td.pac，如 xxx.fasta_td.pac）
      # pFind names the td.pac by APPENDING _td.pac to the fasta name (xxx.fasta_td.pac)
      pac_src <- paste0(fasta_path, "_td.pac")
      if (file.exists(pac_src) && !file.exists(file.path(proj, basename(pac_src)))) {
        file.copy(pac_src, file.path(proj, basename(pac_src)))
      }
      # 产物自动发现（先查后跑）：rescore 的 <项目名>fdr0.01_pUniFind.spectra
      # 已存在时直接登记为策略 4，跳过重复运行（与 pFind 搜库复用同口径）
      # auto-discover the rescore output <project>fdr0.01_pUniFind.spectra first:
      # an existing file is registered as strategy 4 and the run is skipped
      # (the same reuse rule as the pFind search)
      rr <- list.files(file.path(proj, "result"), pattern = "pUniFind\\.spectra$",
                       full.names = TRUE)
      if (length(rr) > 0) {
        auto_files <- c(auto_files, tail(rr, 1)); auto_types <- c(auto_types, "punifind_rescore")
        got4 <- TRUE
      }
      # 自动 rescore：启用开关且有安装路径且产物不存在时由软件调官方 CLI 完成
      # （产物写回 result）
      # auto rescore: when enabled with an install path and no existing output, the
      # official CLI is run by the software (the output lands back in result)
      if (!got4 && auto_rescore && nzchar(puni_exe) && !tq_on) {
        # v1.0：two_queue_scheduler=1 时 rescore 归 GPU 管家，master 不重试
        # v1.0: with the two-queue scheduler on, rescore belongs to the GPU
        # steward; the master never retries it in the foreground
        # 时间线模块横幅（网页时间线按此解析）：pUniFind 单独成行，附墙钟开始/结束
        # timeline module banners (parsed by the web timeline): pUniFind gets its
        # own row with wall-clock start/end; failure falls through to the stop point
        cat("\n==== [总程序 / Master] pUniFind 自动 rescore 开始 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
        cat("==== [Master] pUniFind auto rescore start @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
        ok4 <- tryCatch({
          # v1.0（2026-09-26）：改走 GPU 适配模块（数据目录前置 + flag 面探测
          # + python 直调路线 + 官方开关环境变量 + qry.res 序号压实）；旧
          # run_punifind_cli 保留为回退实现
          pf_run("rescore", puni_exe, file.path(proj, "result"),
                 file.path(cfg2$config_dir, "unlabeled"), proj,
                 converter_py, puni_py,
                 num_workers = pf_nw, num_proc = pf_np,
                 batch_size = pf_bs, chunk_parallel = pf_cp, chunk_cleanup = pf_cu)
          TRUE
        }, error = function(e) {
          msg <- conditionMessage(e)
          mm <- regmatches(msg, regexpr("exit [0-9]+", msg))
          code <- if (length(mm) > 0) sub("exit ", "", mm) else "1"
          cat("[总程序 / Master] pUniFind 自动 rescore 失败 (exit ", code, ")\n", sep = "")
          cat("[Master] pUniFind auto rescore failed (exit ", code, ")\n", sep = "")
          cat("pUniFind 自动 rescore 错误信息：", msg, "\n")
          cat("pUniFind auto rescore error: ", msg, "\n")
          FALSE
        })
        if (isTRUE(ok4)) {
          cat("==== [总程序 / Master] pUniFind 自动 rescore 完成 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
          cat("==== [Master] pUniFind auto rescore done @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
          pt_mark("pUniFind 自动 rescore", "pUniFind auto rescore")
        } else {
          # v1.0（2026-10-01，用户需求：全自动版本判定）：rescore 失败时不再
          # 停点等待手动干预。检测到老构建特征（AssertionError / 参数不兼容 /
          # Windows 多进程 pickle 错误）后自动跳过 rescore 并继续执行 denovo，
          # 以 3 引擎（pFind open + restrict + pUniFind denovo）完成合并。
          # 下次运行时如 rescore 产物仍不存在会再次尝试。
          cat("[pUniFind 适配] rescore 自动运行失败，已自动跳过（不影响 denovo 与后续合并）\n")
          cat("[pUniFind compat] rescore auto-run failed, auto-skipping (denovo and merge continue)\n")
          cat("[pUniFind 适配] 当前引擎组合：pFind open + pFind restrict + pUniFind denovo（3 策略）\n")
          cat("[pUniFind compat] current engine set: pFind open + pFind restrict + pUniFind denovo (3 strategies)\n")
        }
      }
      # 自动运行后再发现一次（自动模式产物）
      # re-discover after the run (the auto-mode output)
      rr <- list.files(file.path(proj, "result"), pattern = "pUniFind\\.spectra$",
                       full.names = TRUE)
      if (length(rr) > 0 && !got4) {
        auto_files <- c(auto_files, tail(rr, 1)); auto_types <- c(auto_types, "punifind_rescore")
        got4 <- TRUE
      }
    }
    # GUI 产出自动发现：denovo 的 pUniFind_result/*.csv（GUI 无表头 9 列格式，取最新一个）→ 策略 5
    # auto-discover the GUI denovo output pUniFind_result/*.csv (the GUI's headerless
    # 9-column layout; the newest file wins) -> strategy 5
    dproj <- file.path(cfg2$config_dir, "pUniFind_denovo")
    # 预先创建 denovo 结果接收文件夹（自动模式产物与手动拷贝共用）
    # create the denovo result inbox folder up front (shared by the auto run and manual copy-in)
    dir.create(file.path(dproj, "pUniFind_result"), recursive = TRUE, showWarnings = FALSE)
    # 自动 denovo：启用开关且有安装路径时由软件调官方 CLI 完成
    # auto denovo: when enabled with an install path, the official CLI is run by the software
    got5 <- FALSE
    # 产物自动发现（先查后跑）：pUniFind_result 里已有 csv 时直接登记为策略 5，
    # 跳过重复运行（与 pFind 搜库复用同口径）
    # auto-discover the output first: an existing csv under pUniFind_result is
    # registered as strategy 5 and the run is skipped (same reuse rule as pFind)
    csvs <- list.files(file.path(dproj, "pUniFind_result"), pattern = "\\.csv$",
                       full.names = TRUE)
    if (length(csvs) > 0) {
      # P2-1（2026-09-26）：优先登记收尾产物（_5_filtered/_5_merged），仅当不存在
      # 时退回流式 *_5_0.csv 并双语提示
      # final_csvs：收尾产物（推理+过滤结束才出现）；chosen_csv：本次登记的文件。
      # P2-1（2026-09-26）：确定性优先 _5_merged（pUniFind 工作流的规范输出，
      # 其 e2e golden 即比对 merged；filtered 是 percolator/FDR 过滤件，与
      # isoPSV 自身的 FDR 流程叠加语义不清）——merged 缺失才退回 filtered，
      # 再缺失才退回流式 *_5_0.csv。
      final_csvs <- csvs[grepl("_5_(filtered|merged)\\.csv$", basename(csvs))]
      merged_csvs <- final_csvs[grepl("_5_merged\\.csv$", basename(final_csvs))]
      filtered_csvs <- final_csvs[grepl("_5_filtered\\.csv$", basename(final_csvs))]
      chosen_csv <- if (length(merged_csvs) > 0) tail(merged_csvs, 1)
        else if (length(filtered_csvs) > 0) tail(filtered_csvs, 1)
        else {
          cat("[pUniFind denovo] 未发现收尾产物（_5_filtered/_5_merged），退回流式 *_5_0.csv；若 denovo 仍在运行，请等收尾后重跑本步\n")
          tail(csvs, 1)
        }
      cat("[pUniFind denovo] 登记产物：", basename(chosen_csv), "（",
          file.info(chosen_csv)$size, " 字节）\n", sep = "")
      auto_files <- c(auto_files, denovo_csv_to_spectra(chosen_csv,
                                                       file.path(dproj, "pUniFind_result")))
      auto_types <- c(auto_types, "punifind_denovo")
      got5 <- TRUE
    }
    if (!got5 && auto_denovo && nzchar(puni_exe) && !tq_on) {
      # v1.0：two_queue_scheduler=1 时 denovo 归 GPU 管家，master 不重试
      # v1.0: with the two-queue scheduler on, denovo belongs to the GPU
      # steward; the master never retries it in the foreground
      # 时间线模块横幅（网页时间线按此解析）：pUniFind 单独成行，附墙钟开始/结束
      # timeline module banners (parsed by the web timeline): pUniFind gets its
      # own row with wall-clock start/end; failure falls through to the stop point
      cat("\n==== [总程序 / Master] pUniFind 自动 denovo 开始 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
      cat("==== [Master] pUniFind auto denovo start @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
      ok5 <- tryCatch({
        # v1.0（2026-09-26）：同 rescore，改走 GPU 适配模块（python 直调路线
        # 需 puni_py，否则回退 exe 启动器——新构建下 exe 路线会重组 argv 失败）
        # 全量检测 FAIL 修复（2026-09-26，deepseek v1.0 全量检测报告 C5）：
        # unlabeled 目录在步骤 2 期间会混入 pFind 写的 *_check.mgf（无
        # RTINSECONDS），pUniFind denovo 逐个读目录内每个 MGF → KeyError:
        # 'rtinseconds'。调用前移入 _stray 子目录（同卷移动瞬时；下游各
        # 枚举点与 load_unlabeled_spectra 均已排除该模式）
        unlab_dir <- file.path(cfg2$config_dir, "unlabeled")
        stray_mgfs <- list.files(unlab_dir, pattern = "_check\\.mgf$", full.names = TRUE)
        if (length(stray_mgfs) > 0) {
          # 实弹修复（2026-09-26，deepseek v1.0 报告）：_stray 建在 unlabeled/
          # 之内会被 denovo 加载器当文件打开（PermissionError）——移到
          # config_dir 下（unlabeled 之外）
          stray_dir <- file.path(cfg2$config_dir, "_stray_mgf")
          dir.create(stray_dir, recursive = TRUE, showWarnings = FALSE)
          file.rename(stray_mgfs, file.path(stray_dir, basename(stray_mgfs)))
          cat("[pUniFind denovo] 移出 ", length(stray_mgfs),
              " 个 *_check.mgf 残留（无 RT，denovo 无法读取）→ _stray_mgf/\n", sep = "")
        }
        pf_run("denovo", puni_exe, file.path(dproj, "result"),
               unlab_dir, dproj,
               python_exe = puni_py,
               num_workers = pf_nw, num_proc = pf_np,
               batch_size = pf_bs, chunk_parallel = pf_cp, chunk_cleanup = pf_cu)
        TRUE
      }, error = function(e) {
        msg <- conditionMessage(e)
        mm <- regmatches(msg, regexpr("exit [0-9]+", msg))
        code <- if (length(mm) > 0) sub("exit ", "", mm) else "1"
        cat("[总程序 / Master] pUniFind 自动 denovo 失败 (exit ", code, ")\n", sep = "")
        cat("[Master] pUniFind auto denovo failed (exit ", code, ")\n", sep = "")
        cat("pUniFind 自动 denovo 错误信息：", msg, "\n")
        cat("pUniFind auto denovo error: ", msg, "\n")
        FALSE
      })
      if (isTRUE(ok5)) {
        cat("==== [总程序 / Master] pUniFind 自动 denovo 完成 @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
        cat("==== [Master] pUniFind auto denovo done @ ", format(Sys.time(), "%H:%M:%S"), " ====\n", sep = "")
        pt_mark("pUniFind 自动 denovo", "pUniFind auto denovo")
      }
    }
    # 自动运行后再发现一次（自动模式产物）
    # re-discover after the run (the auto-mode output)
    csvs <- list.files(file.path(dproj, "pUniFind_result"), pattern = "\\.csv$",
                       full.names = TRUE)
    if (length(csvs) > 0 && !got5) {
      # P2-1（2026-09-26）：优先登记收尾产物（_5_filtered/_5_merged），仅当不存在
      # 时退回流式 *_5_0.csv 并双语提示
      # final_csvs：收尾产物（推理+过滤结束才出现）；chosen_csv：本次登记的文件。
      # P2-1（2026-09-26）：确定性优先 _5_merged（pUniFind 工作流的规范输出，
      # 其 e2e golden 即比对 merged；filtered 是 percolator/FDR 过滤件，与
      # isoPSV 自身的 FDR 流程叠加语义不清）——merged 缺失才退回 filtered，
      # 再缺失才退回流式 *_5_0.csv。
      final_csvs <- csvs[grepl("_5_(filtered|merged)\\.csv$", basename(csvs))]
      merged_csvs <- final_csvs[grepl("_5_merged\\.csv$", basename(final_csvs))]
      filtered_csvs <- final_csvs[grepl("_5_filtered\\.csv$", basename(final_csvs))]
      chosen_csv <- if (length(merged_csvs) > 0) tail(merged_csvs, 1)
        else if (length(filtered_csvs) > 0) tail(filtered_csvs, 1)
        else {
          cat("[pUniFind denovo] 未发现收尾产物（_5_filtered/_5_merged），退回流式 *_5_0.csv；若 denovo 仍在运行，请等收尾后重跑本步\n")
          tail(csvs, 1)
        }
      cat("[pUniFind denovo] 登记产物：", basename(chosen_csv), "（",
          file.info(chosen_csv)$size, " 字节）\n", sep = "")
      auto_files <- c(auto_files, denovo_csv_to_spectra(chosen_csv,
                                                       file.path(dproj, "pUniFind_result")))
      auto_types <- c(auto_types, "punifind_denovo")
      got5 <- TRUE
    }
    # 有缺项：按各项模式（0=不运行不提示；1=自动已尝试；2=手动待办）打印双语指引；
    # 未齐时本步不合并、不生成 peptide.txt，直接重跑总程序续接
    # anything missing: print the bilingual guidance per item mode (0=off, no prompt;
    # 1=auto already attempted; 2=manual to-do); while items are missing this step
    # does not merge and writes no peptide.txt - just rerun the master program
    need4 <- puni_rescore_mode != "0" && !got4
    need5 <- puni_denovo_mode != "0" && !got5
    if (need4 || need5) {
      cat("\n[待办] pUniFind 以下项目未齐：\n")
      cat("[To-do] the following pUniFind items are missing:\n")
      if (need4) {
        if (manual_rescore) {
          cat("  1) rescore（手动）：GUI 运行目录的 result\\<项目名>fdr0.01_pUniFind.spectra（如 testfdr0.01_pUniFind.spectra）\n")
          cat("  1) rescore (manual): from the GUI run folder, result\\<project>fdr0.01_pUniFind.spectra (e.g. testfdr0.01_pUniFind.spectra)\n")
          cat("     → 拷贝到：", file.path(proj, "result"), "\\\n")
          cat("     -> copy to: ", file.path(proj, "result"), "\\\n")
        } else {
          cat("  1) rescore（自动）：未产出结果，请检查上一步 rescore 日志后重跑本程序。\n")
          cat("  1) rescore (auto): no output produced; check the rescore log above and rerun.\n")
        }
      }
      if (need5) {
        if (manual_denovo) {
          cat("  2) denovo（手动）：GUI 运行目录的 pUniFind_result\\test_5_0.csv（取最新一次运行）\n")
          cat("  2) denovo (manual): from the GUI run folder, pUniFind_result\\test_5_0.csv (the newest run)\n")
          cat("     → 拷贝到：", file.path(cfg2$config_dir, "pUniFind_denovo", "pUniFind_result"), "\\\n")
          cat("     -> copy to: ", file.path(cfg2$config_dir, "pUniFind_denovo", "pUniFind_result"), "\\\n")
        } else {
          cat("  2) denovo（自动）：未产出结果，请检查上一步 denovo 日志后重跑本程序。\n")
          cat("  2) denovo (auto): no output produced; check the denovo log above and rerun.\n")
        }
      }
      cat("  完成后直接重跑总程序：软件自动登记为策略 4/5 并五策略合并（pFind 搜库自动复用）。\n")
      cat("  Then simply rerun the master program: the files are auto-registered as strategies 4/5 and the five-strategy merge runs (pFind results are reused).\n")
    }
    # 启用项未齐：本步停在此处，不合并、不生成 peptide.txt——
    # 完成后直接重跑总程序续接（pFind 搜库自动复用、pUniFind 产出自动登记为策略 4/5）
    # When enabled items are not all in: stop here WITHOUT merging or generating
    # peptide.txt - after completing them simply rerun the master program (pFind
    # results are reused, the pUniFind outputs are auto-registered as strategies 4/5)
    if (use_puni && (need4 || need5)) {
      cat("\n==== [等待输入] pUniFind 两项 GUI 完成后直接重跑本程序（无需删除文件）====\n")
      cat("==== [Awaiting input] after finishing both pUniFind GUI items, simply rerun this program (no file needs deleting) ====\n")
      # 退出码 7＝停点约定：总程序识别后正常收尾（不再继续后续步骤）
      # exit status 7 = the stop-point convention: the master recognizes it and ends cleanly
      quit(status = 7)
    }
  }
  # 用户未提供搜库结果时：以自动搜库＋pUniFind GUI 产出作为后续校正/合并的输入
  # When no user search results exist: use the auto-search and pUniFind GUI outputs as the revise/merge input
  if (length(search_files) == 0) {
    search_files <- auto_files
    search_types <- auto_types
  }
}
# 自动搜库产物（b/c）与用户提供的外部引擎条目是并列策略：两者都必须进合并
# （2026-09-13 修复：此前用户提供 search_file 时自动搜库的 b/c 被丢弃，
# 竞争判定少了 pFind 开放/限定两路证据）
# the auto-searched outputs (b/c) and the user-provided engine files are
# parallel strategies: BOTH must enter the merge (2026-09-13 fix: previously
# the auto b/c were dropped whenever the user supplied any search_file, so the
# competition lost the pFind open/restrict evidence)

# ---- 外部引擎自动运行（用户口径 2026-09-17）：isoPSV 调用引擎软件自动执行 ----
# ---- external-engine auto-run (user ruling 2026-09-17): isoPSV invokes the
# ---- engine software itself -------------------------------------------------
# 键：<engine>_exe（可执行文件/jar 路径）、<engine>_params（参数模板路径）、
# <engine>_mode（0 手动 / 1 自动）。自动模式：逐 unlabeled MGF 重编号→运行
# 引擎→适配器转换→拼接为 <engine>_auto.spectra 进 auto_files（产物存在且新于
# MGF 时复用）。msfragger/comet 直接吃 MGF；sage 经 msconvert 转 mzML；
# msfragger_boost（FragPipe headless 工作流）与 maxquant（mqpar.xml）在提供
# params 时运行，否则打印双语待办。
# Keys: <engine>_exe (executable/jar), <engine>_params (parameter template),
# <engine>_mode (0 manual / 1 auto). Auto mode: per unlabeled MGF renumber ->
# run engine -> adapter convert -> concatenate to <engine>_auto.spectra into
# auto_files (reused when the output is newer than the MGFs). msfragger/comet
# take MGFs directly; sage goes through msconvert to mzML; msfragger_boost
# (FragPipe headless workflow) and maxquant (mqpar.xml) run when params are
# provided, otherwise a bilingual to-do is printed.
cfg_get <- function(key, default = NULL) {
  v <- cfg2[[key]]
  if (is.null(v) || length(v) == 0) default else tail(v, 1)
}
# Python 解释器探测（P1-4，2026-09-25）：python → py -3 启动器 → python3 回退；
# 返回 list(cmd, pre) 或 NULL。F-6（2026-09-26）公共化至 compat.R，与
# sysinfo_probe.R 共用同一实现。
engine_auto_run <- function(engine, exe, params_tpl, mode, mgfs, work_dir,
                            fasta_path, mgf_dir, code_dir, mgf_name = NULL,
                            enzyme = "trypsin", en_explicit = TRUE,
                            label_mass = NULL, cand_scope = "filtered") {
  mode <- if (length(mode) == 0) "0" else mode
  # v1.0（2026-09-28）：无标记通道模式（unlabeled_channel=0）下把 NeuCode
  # 类标记修饰注入引擎参数；解析失败（质量 NA）时不注入（模板即真相）
  if (!is.null(label_mass) && (nrow(label_mass) == 0 || anyNA(label_mass$mass)))
    label_mass <- NULL
  # v1.0（2026-09-26，审查 P1-1/P1-2/P1-3）：en_explicit=FALSE 时"模板即真相"
  # ——不改写任何引擎模板，en 仅用于日志/侧车；en_explicit=TRUE 时未知酶报错，
  # 不支持的酶×引擎组合显式跳过该槽位（不再静默按模板胰蛋白酶搜库、结果也不
  # 进入竞争合并）。
  if (en_explicit) {
    w <- character(0)
    en <- withCallingHandlers(enzyme_resolve(enzyme),
                              warning = function(ww) {
                                w <<- c(w, conditionMessage(ww)); invokeRestart("muffleWarning")
                              })
    if (length(w) > 0) stop("未知酶切类型：", enzyme, "\nunknown enzyme type: ", enzyme)
  } else en <- "template"
  if (mode != "1") return(character(0))
  if (length(exe) == 0 || !nzchar(exe) || !file.exists(exe))
    stop(engine, " 自动模式（mode=1）需有效 ", engine, "_exe\n",
         engine, " auto mode requires a valid ", engine, "_exe")
  if (en_explicit) {
    if (en == "none" && engine %in% c("sage", "maxquant")) {
      cat("[", engine, " 自动] enzyme=none 该引擎不支持（无酶切语法随版本而异），跳过该槽位、不参与合并\n", sep = "")
      cat("[", engine, " auto] enzyme=none unsupported by this engine, slot skipped (excluded from the merge)\n")
      return(character(0))
    }
    if (en == "lysarginase" && engine == "metamorpheus") {
      cat("[metamorpheus 自动] enzyme=lysarginase 不在 MetaMorpheus 蛋白酶列表内，跳过该槽位、不参与合并\n")
      cat("[metamorpheus auto] enzyme=lysarginase not in MetaMorpheus proteases, slot skipped (excluded from the merge)\n")
      return(character(0))
    }
  }
  out_spec <- file.path(work_dir, paste0(engine, "_auto.spectra"))
  newest_mgf <- if (length(mgfs) > 0) max(file.mtime(mgfs)) else 0
  # P2-2（2026-09-26）：复用前比对参数侧车（enzyme 变化则强制重跑）
  rp_file <- file.path(work_dir, paste0(engine, "_auto.run_params"))
  # B-3/G-2（2026-09-26）：侧车加 fasta 与模板指纹——换库或"模板即真相"模式下
  # 手改模板酶行都强制重跑（旧侧车只比 enzyme，换 FASTA 静默复用旧结果）
  tpl_stat2 <- file.info(params_tpl)
  rp_want <- c(paste0("enzyme=", en),
               paste0("fasta=", fasta_path),
               paste0("tpl=", format(tpl_stat2$mtime, "%Y-%m-%d %H:%M:%OS3"), "|", tpl_stat2$size),
               # v1.0（工作站反馈#14）：候选口径入侧车——从 filtered 切 top5
               # 必须重搜，复用旧口径产物＝静默吞掉全候选开关
               paste0("cand_scope=", cand_scope))
  # v1.0：标记注入变化须强制重跑（旧侧车只比 enzyme/fasta/模板）
  if (!is.null(label_mass))
    rp_want <- c(rp_want, paste0("label=", paste(sprintf("%s|%s|%.6f",
                label_mass$name, label_mass$residue, label_mass$mass),
                collapse = ",")))
  rp_ok <- file.exists(rp_file) &&
    identical(trimws(readLines(rp_file, warn = FALSE)), rp_want)
  if (file.exists(out_spec) && file.mtime(out_spec) >= newest_mgf && rp_ok) {
    cat("[", engine, " 自动] 产物存在且较新，复用：", out_spec, "\n")
    cat("[", engine, " auto] output exists and is fresh, reuse: ", out_spec, "\n")
    return(out_spec)
  }
  if (file.exists(out_spec) && file.mtime(out_spec) >= newest_mgf && !rp_ok) {
    cat("[", engine, " 自动] 检测到参数变化（enzyme 侧车不匹配），将重新运行\n")
    cat("[", engine, " auto] parameter change detected (enzyme sidecar mismatch), rerunning\n")
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  # MaxQuant 单次全 raw 运行（一次读全部 raw，不走逐 MGF 循环）
  # MaxQuant single all-raw run (reads every raw at once, not the per-MGF loop)
  if (engine == "maxquant") {
    # v1.0：兼容网页壳子的 search_file_raw 键（旧键 raw_folder 仍优先）
    # v1.0: also accept the web-shell's search_file_raw key
    raw_dir <- cfg_get("raw_folder", "")
    raws <- if (nzchar(raw_dir)) sort(list.files(raw_dir, pattern = "\\.raw$",
                                                 full.names = TRUE)) else character(0)
    if (length(raws) == 0) {
      # 回退到 search_file_raw（网页写入；可为文件夹或分号/换行分隔的文件列表）
      sfr <- cfg_get("search_file_raw", "")
      if (nzchar(sfr)) {
        sfr_items <- unique(trimws(strsplit(gsub(";", "\n", sfr), "\n")[[1]]))
        sfr_items <- sfr_items[nzchar(sfr_items)]
        dir_items <- sfr_items[dir.exists(sfr_items)]
        file_items <- sfr_items[file.exists(sfr_items) & grepl("\\.raw$", sfr_items, ignore.case = TRUE)]
        for (dd in dir_items)
          file_items <- c(file_items, sort(list.files(dd, pattern = "\\.raw$",
                                                       full.names = TRUE)))
        raws <- sort(unique(file_items))
      }
    }
    if (length(raws) == 0) {
      cat("[maxquant 自动] raw_folder/search_file_raw 未配置或无 .raw 文件，跳过\n")
      cat("[maxquant auto] raw_folder/search_file_raw unset or no .raw files, skipped\n")
      return(character(0))
    }
    newest_raw <- max(file.mtime(raws))
    if (file.exists(out_spec) && file.mtime(out_spec) >= newest_raw && rp_ok) {
      cat("[maxquant 自动] 产物存在且较新，复用：", out_spec, "\n")
      cat("[maxquant auto] output exists and is fresh, reuse: ", out_spec, "\n")
      return(out_spec)
    }
    if (file.exists(out_spec) && file.mtime(out_spec) >= newest_raw && !rp_ok) {
      cat("[maxquant 自动] 检测到参数变化（enzyme 侧车不匹配），将重新运行\n")
    }
    logf <- file.path(work_dir, "maxquant.log")
    tpl <- readLines(params_tpl, warn = FALSE)
    # MaxQuant 的 mqpar 内路径必须是反斜杠 Windows 风格（正斜杠会在其
    # Substring 派生基目录时崩溃）；fasta/raw/输出目录统一转换
    fa_bs <- gsub("/", "\\\\", fasta_path)
    raw_bs <- gsub("/", "\\\\", raws)
    out_bs <- gsub("/", "\\\\", work_dir)
    # 占位替换用 fixed=TRUE（替换串含反斜杠，正则替换会吃掉 \2 等序列）
    # placeholder substitution with fixed=TRUE (the replacements carry
    # backslashes, which the regex engine would consume as \2 etc.)
    for (k in seq_along(tpl)) {
      tpl[k] <- sub("{FASTA}", fa_bs, tpl[k], fixed = TRUE)
      tpl[k] <- sub("{OUT}", out_bs, tpl[k], fixed = TRUE)
    }
    raw_lines <- paste0("      <string>", raw_bs, "</string>")
    exp_lines <- paste0("      <string>", sub("\\.raw$", "", basename(raws)), "</string>")
    frac_lines <- rep("      <short>32767</short>", length(raws))
    rep_blk <- function(tag, lines) {
      i <- grep(tag, tpl)
      if (length(i) == 1) {
        tpl <<- c(tpl[seq_len(i - 1L)], lines,
                  if (i < length(tpl)) tpl[(i + 1L):length(tpl)])
      }
    }
    rep_blk("\\{RAWS\\}", raw_lines)
    rep_blk("\\{EXP\\}", exp_lines)
    rep_blk("\\{FRAC\\}", frac_lines)
    # v1.0（2026-09-26）：酶切类型改写 mqpar 的 <enzymes> 块第一个 <string>
    # （模板为 Trypsin/P；none/MaxQuant 不支持的酶保持模板并提示）
    # v1.0（2026-09-26）：酶切类型改写 mqpar 的 <enzymes> 块第一个 <string>。
    # v1.0：锚定 <enzymes> 块（E-3 消除未锚定全局 grep）；仅 en_explicit 时
    # 改写（键缺失 = 模板即真相），none 已在入口跳过。
    if (en_explicit) {
      mq_name <- switch(en, trypsin = "Trypsin",
                        lysc = "LysC", lysn = "LysN", argc = "ArgC", aspn = "AspN",
                        gluc = "GluC", chymotrypsin = "Chymotrypsin",
                        lysarginase = "LysargiNase", "")
      if (nzchar(mq_name)) {
        i_start <- which(trimws(tpl) == "<enzymes>")
        i_end <- which(trimws(tpl) == "</enzymes>")
        if (length(i_start) == 1 && length(i_end) == 1 && i_end[1] > i_start[1]) {
          js <- which(grepl("<string>", tpl[(i_start[1] + 1):(i_end[1] - 1)], fixed = TRUE))
          if (length(js) >= 1) {
            j <- i_start[1] + js[1]
            tpl[j] <- sub("<string>[^<]*</string>", paste0("<string>", mq_name, "</string>"), tpl[j])
            cat("[maxquant 自动] enzyme → ", mq_name, "\n", sep = "")
          }
        }
      } else {
        cat("[maxquant 自动] enzyme ", en, " 不在 MaxQuant 支持列表内，跳过该槽位、不参与合并\n", sep = "")
        return(character(0))
      }
    }
    mqpar <- file.path(work_dir, "mqpar.xml")
    writeLines(tpl, mqpar, useBytes = TRUE)
    st <- system2(exe, shQuote(mqpar), stdout = logf, stderr = logf)
    # MaxQuantCmd 忽略 fixedCombinedFolder、固定把 combined 写到第一个 raw 所在
    # 目录（实测 2.8.1.0），msms.txt 按该位置取
    msms <- file.path(dirname(raws[1]), "combined", "txt", "msms.txt")
    if (!file.exists(msms)) {
      cat("[maxquant 自动] 无 msms.txt（exit ", st, "）\n", sep = "")
      return(character(0))
    }
    out_i <- file.path(work_dir, "maxquant_conv.spectra")
    rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
    st2 <- iso_system2(rscript_bin,
                       c(file.path(code_dir, "search_adapter_maxquant.R"), msms, out_i),
                       stdout = logf, stderr = logf)
    if (st2 == 0 && file.exists(out_i)) {
      cat("[maxquant 自动] 完成：", out_i, "\n")
      cat("[maxquant auto] done: ", out_i, "\n")
      writeLines(rp_want, rp_file, useBytes = TRUE)   # P2-2 参数侧车（enzyme+fasta+模板）
      return(out_i)
    }
    cat("[maxquant 自动] 适配器转换失败（exit ", st2, "）\n", sep = "")
    return(character(0))
  }
  parts <- character(0)
  # N-2（2026-09-25）：-Xmx 只在需要时探测一次（旧实现放在逐 MGF 循环内，
  # 每个文件触发一次慢速 systeminfo）；提出循环并缓存于本函数作用域。
  # N-2 (2026-09-25): probe -Xmx once, lazily (the old code sat inside the
  # per-MGF loop, triggering the slow systeminfo once per file); hoisted here.
  xmx_cached <- NULL
  msfragger_xmx <- function() {
    if (!is.null(xmx_cached)) return(xmx_cached)
    xmx <- cfg_get("msfragger_xmx", "")
    if (!nzchar(xmx)) {
      tm <- system_total_mb()
      gb <- if (!is.na(tm) && tm > 0) min(0.6 * tm / 1024, 15) else 4
      xmx <- paste0(round(gb, 1), "G")
    }
    xmx_cached <<- xmx
    xmx
  }
  # F-7（2026-09-26）：Python 探测同样提出循环（原在逐 MGF 循环内，每文件
  # Sys.which ×3）；函数定义已公共化至 compat.R（F-6），此处只调用一次。
  # v1.0：先看 config_step1 的 python_exe 显式键（用户 Python 不在 PATH
  # 时填完整路径），再回退环境目录探测。
  # v1.0（工作站反馈#16）：此处原经 cfg2-only 的 cfg_get 读取，与"step1
  # 优先"口径不符；改 cfg_get_preferred，模板占位符（含 <>）视为空。
  py_hint <- cfg_get_preferred(base_cfg, cfg2, "python_exe", "")$value
  py_hint <- as.character(py_hint)[1]
  if (is.na(py_hint) || grepl("[<>]", py_hint)) py_hint <- ""
  pyenv <- py_detect(py_hint)
  # ---- P4（2026-09-27 用户口径）：逐 raw/mgf 单独搜库 + 并行搜库 ----
  # 循环体抽为闭包；<engine>_workers>1 时经 PSOCK 并行执行，产物按 raw 序收集；
  # 并行时各引擎的线程与内存按 worker 数摊薄（num_threads 与 -Xmx 在分支内缩放）。
  ew0 <- suppressWarnings(as.integer(tail(c(cfg_get(paste0(engine, "_workers"), "1"), "1"), 1)))
  if (is.na(ew0) || ew0 < 1L) ew0 <- 1L
  ncores_avail <- max(parallel::detectCores() - 1L, 1L)
  thr_per_w <- max(floor(ncores_avail / ew0), 1L)
  run_one_raw <- function(mgf) {
    stem <- sub("\\.mgf$", "", basename(mgf))
    rmgf <- file.path(work_dir, paste0(stem, "_renum.mgf"))
    mapf <- file.path(work_dir, paste0(stem, "_renum_map.tsv"))
    # P1-4（2026-09-25）：Python 探测 + 返回码强制检查（旧写法硬用 "python"
    # 且忽略失败，下游拿到的是不存在的重编号 MGF）
    if (is.null(pyenv)) {
      stop("renumber_mgf.py 需要 Python（python / py -3 / python3 均未检出）；",
           "若 Python 未加入环境目录（PATH），请在 config_step1.txt 设 ",
           "python_exe 为解释器完整路径，或安装 Python 3.8+ 后重跑。\n",
           "renumber_mgf.py needs Python (none of python / py -3 / python3 found); ",
           "if Python is not on PATH, set python_exe in config_step1.txt to the ",
           "full interpreter path, or install Python 3.8+ and rerun.")
    }
    stp <- iso_system2(pyenv$cmd,
                       c(pyenv$pre, file.path(code_dir, "renumber_mgf.py"),
                         mgf, rmgf, mapf),
                       stdout = file.path(work_dir, paste0(stem, "_renum.log")),
                       stderr = file.path(work_dir, paste0(stem, "_renum.log")))
    if (stp != 0 || !file.exists(rmgf)) {
      stop("重编号 MGF 生成失败（exit ", stp, "）：", stem, "\n",
           "renumber_mgf.py failed (exit ", stp, "): ", stem)
    }
    logf <- file.path(work_dir, paste0(stem, "_", engine, ".log"))
    if (engine == "msfragger") {
      prm <- file.path(work_dir, paste0(stem, "_msfragger.params"))
      pl <- readLines(params_tpl, warn = FALSE)
      i_db <- grep("^database_name", pl)
      if (length(i_db) == 1) pl[i_db] <- paste0("database_name = ", fasta_path)
      # v1.0（2026-09-28 用户口径）：全候选模式（ISOPSV_CAND_ALL=1）下把
      # 引擎输出深度改写为每谱前 5（output_report_topN），否则引擎自身只吐
      # 每谱 1 个候选、top5 合并池无从膨胀
      if (identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
        i_top <- grep("^output_report_topN", pl)
        if (length(i_top) == 1) pl[i_top] <- "output_report_topN = 5"
      }
      # v1.0：标记修饰注入（编号续接模板已有 variable_mod_NN，避免撞车）
      # v1.0：engine="msfragger"——端部标记折算 [^ / ]$ / n<c>前缀语法
      if (!is.null(label_mass))
        pl <- label_mod_lines(pl, label_mass, "^variable_mod_[0-9]",
                              "variable_mod_%02d = %.6f %s 3", engine = "msfragger")
      # v1.0（2026-09-26）：酶切类型改写 MSFragger 酶 1 四行（none → nonspecific）；
      # v1.0：仅 en_explicit 时改写（键缺失 = 模板即真相）
      if (en_explicit) {
        d <- enzyme_def(en)
        set_line <- function(key, val) {
          i <- grep(paste0("^", key), pl)
          if (length(i) == 1) pl[i] <<- paste0(key, " = ", val)
        }
        set_line("search_enzyme_name_1", if (en == "none") "nonspecific" else tolower(d$label))
        set_line("search_enzyme_cut_1", d$cut)
        set_line("search_enzyme_nocut_1", d$nocut)
        set_line("search_enzyme_sense_1", d$sense)
      }
      # P4（2026-09-27）：并行时 num_threads 按 worker 数摊薄
      if (ew0 > 1L)
        pl <- gsub("num_threads = [0-9]+",
                   paste0("num_threads = ", thr_per_w), pl)
      writeLines(pl, prm, useBytes = TRUE)
      # P1-5（2026-09-25）：-Xmx 可配置（msfragger_xmx），留空按可用物理内存
      # 60% 自适应、封顶 15G；N-2：探测结果缓存，只跑一次
      xmx <- msfragger_xmx()
      # P4（2026-09-27）：并行时 -Xmx 按 worker 数摊薄（最低 2G）
      if (ew0 > 1L) {
        xm <- suppressWarnings(as.numeric(sub("[gG]$", "", as.character(xmx))))
        if (is.finite(xm)) xmx <- paste0(max(2, floor(xm / ew0)), "G")
      }
      # 本地化+口径修正（2026-09-26）：MSFragger 只读 mzML，此前把 MGF 直接
      # 喂给 java -jar 属错误口径（该分支从未实弹验证）——先转 mzML 再调用；
      # msconvert 配置了就用它，否则用内置纯 R 转换器 mgf_to_mzml
      mzml_msf <- file.path(work_dir, paste0(stem, "_renum_pw.mzML"))
      msconv <- cfg_get("msconvert_exe", "")
      if (nzchar(msconv) && file.exists(msconv)) {
        iso_system2(msconv, c(rmgf, "-o", work_dir, "--outfile", basename(mzml_msf)),
                    stdout = logf, stderr = logf)
      } else {
        # v1.0：共享缓存（0.19 MB/s 的纯 R 转换，跨引擎不重复转）
        mgf_to_mzml_cached(rmgf, mzml_msf)
        cat("[msfragger 自动] 内置 R 转换 MGF→mzML（未配置 msconvert，走共享缓存）\n")
      }
      st <- iso_system2("java", c(paste0("-Xmx", xmx), "-jar", exe, prm, mzml_msf),
                        stdout = logf, stderr = logf)
      cat("[msfragger 自动] -Xmx", xmx, "（物理内存自适应，可经 msfragger_xmx 覆盖）\n", sep = "")
      # P1-2（2026-09-26，工作站 v1.0 复测）：MSFragger 按输入文件名命名
      # 产物（输入 <stem>_renum_pw.mzML ⇒ <stem>_renum_pw.tsv / .pepXML）；
      # 旧回退由 rmgf 推出同名等于没有。取 <stem>_renum_pw.tsv，缺失时在
      # work_dir 下取最新 *_renum_pw.tsv / *.pepXML
      res <- file.path(work_dir, paste0(stem, "_renum_pw.tsv"))
      if (!file.exists(res)) {
        # P4（2026-09-27）：并行模式下 work_dir 共享，回退 glob 必须锚定 stem
        # （否则会捡到其他 raw 并发写出的同名产物）
        cand <- list.files(work_dir,
                           pattern = paste0("^", stem, "_renum_pw\\.(tsv|pepXML|pep\\.xml)$"),
                           full.names = TRUE)
        res <- if (length(cand) > 0) cand[which.max(file.mtime(cand))] else ""
      }
    } else if (engine == "comet") {
      # 实弹修复（2026-09-26，deepseek v1.0 报告）：模板 database_name 为
      # 占位符 __FASTAPATH__，自动分支必须注入项目 FASTA（此前不注入，读到
      # 占位符/残留路径导致 18/18 失败）
      pl0 <- readLines(params_tpl, warn = FALSE)
      i_db0 <- grep("^database_name", pl0)
      if (length(i_db0) == 1) pl0[i_db0] <- paste0("database_name = ", fasta_path)
      # v1.0（2026-09-28 用户口径）：全候选模式下每谱输出前 5
      if (identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
        i_ol <- grep("^num_output_lines", pl0)
        if (length(i_ol) == 1) pl0[i_ol] <- "num_output_lines = 5"
      }
      # v1.0：标记修饰注入（编号续接模板已有 variable_modNN）
      # v1.0：engine="comet"——端部标记折算裸 n/c（任意残基端）
      if (!is.null(label_mass))
        pl0 <- label_mod_lines(pl0, label_mass, "^variable_mod[0-9]",
                               "variable_mod%02d = %.6f %s 0 3 -1 0 0 0.0",
                               engine = "comet")
      prm0 <- file.path(work_dir, paste0(stem, "_comet.params"))
      writeLines(pl0, prm0, useBytes = TRUE)
      params_use <- prm0
      # v1.0/v1.0（2026-09-26）：酶切类型按 Comet 参数文件自身的
      # COMET_ENZYME_INFO 表按名取编号（E-2 硬编码消除：用户改酶表/换 Comet
      # 版本都不会错配）；表内没有的酶（LysargiNase/GluC）按同格式追加编号。
      # none → 0（Cut_everywhere）。仅 en_explicit 时改写（键缺失 = 模板即真相）。
      if (en_explicit) {
        pl <- readLines(prm0, warn = FALSE)
        num <- comet_enzyme_number(en, pl)
        if (is.na(num)) {
          pl <- c(pl, comet_enzyme_append_line(en, pl))
          num <- as.integer(sub("^([0-9]+)\\..*", "\\1", tail(pl, 1)))
        }
        i_en <- grep("^search_enzyme_number", pl)
        if (length(i_en) == 1) pl[i_en] <- paste0("search_enzyme_number = ", num)
        if (ew0 > 1L)
          pl <- gsub("num_threads = [0-9]+",
                     paste0("num_threads = ", thr_per_w), pl)
        writeLines(pl, prm0, useBytes = TRUE)
        cat("[comet 自动] enzyme → ", en, "（编号 ", num, "）\n", sep = "")
      }
      st <- iso_system2(exe, c(paste0("-P", params_use), rmgf),
                        stdout = logf, stderr = logf)
      res <- list.files(work_dir, pattern = paste0("^", stem, ".*\\.txt$"),
                        full.names = TRUE)
      res <- res[!grepl("renum_map", res)]
      res <- if (length(res) > 0) res[order(file.mtime(res), decreasing = TRUE)[1]] else ""
    } else if (engine == "sage") {
      mzml <- file.path(work_dir, paste0(stem, "_renum_pw.mzML"))
      msconv <- cfg_get("msconvert_exe", "")
      # 本地化（2026-09-26）：msconvert 配置了就用它；否则用内置纯 R
      # MGF→mzML 转换器（管线输入本就是 MGF，不再硬依赖 ProteoWizard）
      if (nzchar(msconv) && file.exists(msconv)) {
        iso_system2(msconv, c(rmgf, "-o", work_dir, "--outfile", basename(mzml)),
                    stdout = logf, stderr = logf)
      } else {
        # v1.0：共享缓存（与 msfragger 槽共用同一转换产物）
        mgf_to_mzml_cached(rmgf, mzml)
        cat("[sage 自动] 内置 R 转换 MGF→mzML（未配置 msconvert，走共享缓存）\n")
      }
      js <- file.path(work_dir, paste0(stem, "_sage.json"))
      jt <- readLines(params_tpl, warn = FALSE)
      # sage 0.15 的 object_store 不识别 Windows 盘符路径（小写/大写均报
      # "Unable to recognise URL"，实测 file:/// URL 形式可正确解析为
      # Windows 路径）；json 内路径与 json 文件本身都转 file:/// 形式
      file_url <- function(p) {
        p <- gsub("\\\\", "/", p)
        paste0("file:///", p)
      }
      rmgf <- file_url(rmgf); mzml <- file_url(mzml)
      fa_url <- file_url(fasta_path)
      out_dir <- file.path(work_dir, paste0("sage_", stem))
      jt <- gsub("\\{[A-Za-z_]*MGF\\}", rmgf, jt)
      jt <- gsub("\\{[A-Za-z_]*MZML\\}", mzml, jt)
      jt <- gsub("\\{[A-Za-z_]*FASTA\\}", fa_url, jt)
      jt <- gsub("\\{[A-Za-z_]*OUT\\}", file_url(out_dir), jt)
      # v1.0：标记修饰按残基聚合并入 variable_mods JSON 对象
      if (!is.null(label_mass)) {
        if (any(grepl('"variable_mods"', jt, fixed = TRUE))) {
          jt <- sub('"variable_mods"\\s*:\\s*\\{',
                    paste0('"variable_mods": {\n',
                           sage_label_mod_entries(label_mass), ','), jt)
        } else {
          cat("[sage 自动] 模板无 variable_mods 键，标记修饰注入跳过\n")
        }
      }
      # v1.0（2026-09-26）：酶切类型改写 Sage json 的 cleave_at/c_terminal/
      # restrict（restrict 即 P 位阻断集）；v1.0：仅 en_explicit 时改写
      # （键缺失 = 模板即真相），none 已在入口跳过
      if (en_explicit) {
        d <- enzyme_def(en)
        jt <- sub('"cleave_at": "[^"]*"', sprintf('"cleave_at": "%s"', d$cut), jt)
        jt <- sub('"c_terminal": (true|false)',
                  sprintf('"c_terminal": %s', if (d$sense == "C") "true" else "false"), jt)
        jt <- sub('"restrict": "[^"]*"', sprintf('"restrict": "%s"', d$nocut), jt)
        cat("[sage 自动] enzyme → ", en, "\n", sep = "")
      }
      writeLines(jt, js, useBytes = TRUE)
      # v1.0（2026-10-01，Macaca 报告#1）：sage 大库 OOM 防护。
      # Macaca 实测 72K 蛋白库单笔 28.78 GiB 分配失败；小库 1.6 s 通过。
      # sage_max_db_entries 旋钮：>0 时检查 fasta 条目数，超限则跳过 sage 槽
      # 并打印警告（不终止管线）；0＝不检查（默认）。
      sage_max_db <- suppressWarnings(as.integer(cfg_get("sage_max_db_entries", default = "0")))
      if (!is.na(sage_max_db) && sage_max_db > 0) {
        fa_n <- length(readLines(fasta_path, warn = FALSE)[grep("^>", readLines(fasta_path, warn = FALSE))])
        if (fa_n > sage_max_db) {
          cat("[sage 自动] 数据库条目数 ", fa_n, " 超过 sage_max_db_entries=", sage_max_db,
              "，跳过 sage 槽（大库 OOM 防护）\n", sep = "")
          cat("[sage auto] database has ", fa_n, " entries > sage_max_db_entries=", sage_max_db,
              ", skipping sage slot (large-db OOM guard)\n", sep = "")
          res <- ""
        } else {
          st <- system2(exe, shQuote(file_url(js)), stdout = logf, stderr = logf)
        }
      } else {
        st <- system2(exe, shQuote(file_url(js)), stdout = logf, stderr = logf)
      }
      res <- file.path(work_dir, paste0(stem, "_sage.tsv"))
      if (!file.exists(res)) {
        # P4（2026-09-27）：并行模式下回退 glob 锚定本 raw 的 sage_<stem> 子目录
        # （否则会捡到其他 raw 并发写出的 results.sage.tsv）
        rr <- list.files(file.path(work_dir, paste0("sage_", stem)),
                         pattern = "results\\.sage\\.tsv$",
                         full.names = TRUE, recursive = TRUE)
        res <- if (length(rr) > 0) rr[order(file.mtime(rr), decreasing = TRUE)[1]] else ""
      }
    } else if (engine == "metamorpheus") {
      # mzML 转档（与 sage 同口径，2026-09-26 本地化：msconvert 优先、
      # 内置纯 R 转换器回退，不再硬停）；数据库与谱图经 CLI 选项传入（-d/-s），
      # 模板为默认 SearchTask.toml（-g 生成口径），无需占位替换
      mzml <- file.path(work_dir, paste0(stem, "_renum_pw.mzML"))
      msconv <- cfg_get("msconvert_exe", "")
      if (nzchar(msconv) && file.exists(msconv)) {
        iso_system2(msconv, c(rmgf, "-o", work_dir, "--outfile", basename(mzml)),
                    stdout = logf, stderr = logf)
      } else {
        mgf_to_mzml(rmgf, mzml)
        cat("[metamorpheus 自动] 内置 R 转换 MGF→mzML（未配置 msconvert）\n")
      }
      mm_out <- file.path(work_dir, paste0("mm_", stem))
      # v1.0（2026-09-26）：酶切类型改写 SearchTask.toml 的 Protease /
      # SpecificProtease；v1.0：仅 en_explicit 时改写、空格容差 \\s*=\\s*
      # （E-4 的 Trypsin/P 问题随 trypsin_p 变体删除已闭环；lysarginase 已在
      # 入口跳过）
      params_use <- params_tpl
      if (en_explicit) {
        mm_name <- switch(en, trypsin = "trypsin", lysc = "lys-c", lysn = "lys-n",
                          argc = "arg-c", aspn = "asp-n", gluc = "glu-c",
                          chymotrypsin = "chymotrypsin", none = "non-specific", "")
        if (nzchar(mm_name)) {
          tpl2 <- file.path(work_dir, paste0(stem, "_searchtask.toml"))
          tl <- readLines(params_tpl, warn = FALSE)
          # (^|\\s) 锚定，避免匹配 SpecificProtease 里的 "Protease" 子串
          tl <- gsub('(^|\\s)Protease\\s*=\\s*"[^"]*"', paste0('\\1Protease = "', mm_name, '"'), tl)
          tl <- gsub('(^|\\s)SpecificProtease\\s*=\\s*"[^"]*"', paste0('\\1SpecificProtease = "', mm_name, '"'), tl)
          writeLines(tl, tpl2, useBytes = TRUE)
          params_use <- tpl2
          cat("[metamorpheus 自动] enzyme → ", mm_name, "\n", sep = "")
        } else {
          cat("[metamorpheus 自动] enzyme ", en, " 不在 MetaMorpheus 支持列表内，跳过该槽位、不参与合并\n", sep = "")
          return("")
        }
      }
      st <- system2(exe, c("-t", shQuote(params_use),
                           "-d", shQuote(fasta_path),
                           "-s", shQuote(mzml),
                           "-o", shQuote(mm_out)),
                    stdout = logf, stderr = logf)
      res <- list.files(mm_out, pattern = "AllPSMs\\.psmtsv$",
                        full.names = TRUE, recursive = TRUE)
      res <- if (length(res) > 0) res[order(file.mtime(res), decreasing = TRUE)[1]] else ""
    } else stop("未知引擎自动类型：", engine)
    if (!nzchar(res) || !file.exists(res)) {
      # P1-2（2026-09-26）：日志区分"引擎失败"与"产物未发现"（此前 18/18
      # exit 0 但找错文件名也打成"无结果"，误导排查方向）
      if (isTRUE(st != 0)) {
        cat("[", engine, " 自动] ", stem, " 引擎失败（exit ", st, "）\n", sep = "")
      } else {
        cat("[", engine, " 自动] ", stem, " 引擎 exit 0 但产物未发现（结果文件不在预期位置，日志：", logf, "）\n", sep = "")
      }
      return("")
    }
    adapter <- switch(engine, msfragger = "search_adapter_msfragger.R",
                      comet = "search_adapter_comet.R",
                      sage = "search_adapter_sage.R",
                      metamorpheus = "search_adapter_metamorpheus.R")
    out_i <- file.path(work_dir, paste0(stem, "_conv.spectra"))
    rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
    st2 <- if (engine == "metamorpheus") {
      # 适配器签名：<psmtsv or folder> <out.spectra> [qval] [map.tsv]
      iso_system2(rscript_bin,
                  c(file.path(code_dir, adapter), res, out_i, "0.01", mapf),
                  stdout = logf, stderr = logf)
    } else {
      iso_system2(rscript_bin,
                  c(file.path(code_dir, adapter), res, mapf, out_i),
                  stdout = logf, stderr = logf)
    }
    if (st2 == 0 && file.exists(out_i)) {
      # P3（2026-09-26，工作站反馈）：适配器可能只写表头（results 无有效行）
      # ——仅表头的产物不并入合并，日志指明跳过
      nl <- length(readLines(out_i, n = 2, warn = FALSE))
      if (nl < 2L) {
        cat("[", engine, " 自动] ", stem, " 适配产物仅表头 0 行，跳过该 raw\n", sep = "")
        return("")
      } else return(out_i)
    }
  }
  # ---- P4（2026-09-27）：<engine>_workers（config_step2 键，默认 1=旧行为）>1 时
  # 逐 raw 并行执行；线程与内存由各引擎分支按 worker 数摊薄（xmx/thread 已在
  # 闭包外按 ew 预缩放？——此处只负责并行调度，各分支的线程改写后续版本再摊薄）
  ew <- suppressWarnings(as.integer(tail(c(cfg_get(paste0(engine, "_workers"), "1"), "1"), 1)))
  if (is.na(ew) || ew < 1L) ew <- 1L
  if (ew > 1L && length(mgfs) > 1L) {
    cl <- parallel::makeCluster(min(ew, length(mgfs)))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    ispsv_bootstrap_cluster(cl)
    parallel::clusterExport(cl, c("run_one_raw", "cfg_get"), envir = environment())
    resl <- parallel::parLapply(cl, mgfs, run_one_raw)
  } else {
    resl <- lapply(mgfs, run_one_raw)
  }
  parts <- unlist(resl, use.names = FALSE)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    cat("[", engine, " 自动] 无任何转换产物，跳过\n")
    cat("[", engine, " auto] no converted outputs, skipped\n")
    return(character(0))
  }
  # 拼接（首件带头，其余去头追加）
  heads <- lapply(parts, readLines, n = 1, warn = FALSE)
  bodies <- lapply(parts, function(p) readLines(p, warn = FALSE)[-1])
  lines <- c(heads[[1]], unlist(bodies, use.names = FALSE))
  writeLines(lines, out_spec, useBytes = TRUE)
  writeLines(rp_want, rp_file, useBytes = TRUE)   # P2-2 参数侧车（enzyme+fasta+模板）
  cat("[", engine, " 自动] 完成：", out_spec, "\n")
  cat("[", engine, " auto] done: ", out_spec, "\n")
  out_spec
}
mgf_auto_dir <- file.path(cfg2$config_dir, "unlabeled")
mgfs_auto <- sort(list.files(mgf_auto_dir, pattern = "\\.mgf$", full.names = TRUE))
mgfs_auto <- mgfs_auto[basename(mgfs_auto) != "unlabeled.mgf" &
                    !grepl("_check\\.mgf$", basename(mgfs_auto))]   # P1-1（2026-09-26）
# v1.0（2026-09-26）：酶切类型（config_step1 键 enzyme）；v1.0：键缺失时
# 不改写任何引擎模板（模板即真相），en_explicit=FALSE 传入 engine_auto_run
en_v <- base_cfg$enzyme
en_auto_explicit <- !is.null(en_v) && length(en_v) > 0 && nzchar(trimws(en_v[1]))
if (en_auto_explicit) {
  w_en2 <- character(0)
  en_auto <- withCallingHandlers(enzyme_resolve(en_v[1]),
                                 warning = function(ww) {
                                   w_en2 <<- c(w_en2, conditionMessage(ww)); invokeRestart("muffleWarning")
                                 })
  if (length(w_en2) > 0)
    stop("未知酶切类型：", trimws(en_v[1]), "
unknown enzyme type: ", trimws(en_v[1]))
} else en_auto <- "trypsin"
# v1.0（2026-09-28）：无标记通道模式（unlabeled_channel=0）下，把 NeuCode
# 类标记修饰解析成 (name,residue,mass) 注入外部引擎自动搜库参数
ul_ch_auto <- if (!is.null(base_cfg$unlabeled_channel) && length(base_cfg$unlabeled_channel) > 0)
  tail(base_cfg$unlabeled_channel, 1) else "1"
label_mass_auto <- if (!identical(ul_ch_auto, "1"))
  label_mod_masses(parse_label_spec(base_cfg$label_type)$mods) else NULL
for (eng in c("msfragger", "comet", "sage", "metamorpheus", "msfragger_boost", "maxquant")) {
  exe_v <- cfg_get(paste0(eng, "_exe"), "")
  prm_v <- cfg_get(paste0(eng, "_params"), "")
  mod_v <- cfg_get(paste0(eng, "_mode"), "0")
  if (!nzchar(exe_v)) next
  # 网页版只填 <engine>_exe 时按内置模板回退（2026-09-24 网页自动填充配套）
  # fall back to the bundled template when the web only filled <engine>_exe
  if (!nzchar(prm_v)) {
    tpl_map <- c(msfragger = "engine_templates/msfragger_con.params",
                 comet = "engine_templates/comet_con.params",
                 sage = "engine_templates/sage_con.json",
                 metamorpheus = "engine_templates/metamorpheus_con.toml",
                 maxquant = "engine_templates/maxquant_con.mqpar.xml")
    if (eng %in% names(tpl_map)) {
      cand <- file.path(code_dir, tpl_map[[eng]])
      if (file.exists(cand)) prm_v <- cand
    }
  }
  wd <- file.path(cfg2$config_dir, paste0(eng, "_auto"))
  if (eng %in% c("msfragger", "comet", "sage", "metamorpheus", "maxquant")) {
    sp <- engine_auto_run(eng, exe_v, prm_v, mod_v, mgfs_auto, wd,
                          fasta_path, mgf_auto_dir, code_dir,
                          enzyme = en_auto, en_explicit = en_auto_explicit,
                          label_mass = label_mass_auto, cand_scope = cand_scope)
    if (length(sp) > 0) {
      auto_files <- c(auto_files, sp)
      auto_types <- c(auto_types, eng)
    }
  } else {
    if (mod_v != "1") next
    if (!nzchar(prm_v) || !file.exists(prm_v)) {
      cat("[待办] ", eng, " 自动模式需 ", eng, "_params（FragPipe headless 工作流文件）\n", sep = "")
      next
    }
    out_spec <- file.path(wd, paste0(eng, "_auto.spectra"))
    if (file.exists(out_spec)) {
      auto_files <- c(auto_files, out_spec)
      auto_types <- c(auto_types, eng)
      next
    }
    cat("[", eng, " 自动] 参数已提供但自动运行未实现——请先用 FragPipe（--headless --workflow）跑完后把结果放 search_file 手动登记。\n", sep = "")
  }
}

# ---- 0. 诱饵库去除（用户 2026-09-05 规格）：先清洗，再校正/合并 ----------------
# ---- 0. decoy removal (user 2026-09-05 spec): clean first, then revise/merge --
# P0-1（2026-09-26，GLM 审查）：合并块原在外部引擎自动运行循环之前，引擎循环
# 末尾才把 <engine>_auto.spectra 追加进 auto_files——此后无任何合并语句，六个
# 外部引擎的自动搜库产物被完整计算后静默丢弃（v1.0 的多引擎自动运行端到端
# 实际无效）。合并块移到引擎循环之后统一收口。
# P0-1 (2026-09-26, GLM review): the merge block used to sit BEFORE the
# external-engine auto-run loop, whose tail appends <engine>_auto.spectra to
# auto_files - with no merge afterwards, the six external engines' auto-search
# outputs were computed and then silently dropped (the v1.0 multi-engine
# auto-run was effectively dead end-to-end). The merge now closes the loop
# after the engine loop.
if (length(auto_files) > 0) {
  comb_f <- c(search_files, auto_files)
  comb_t <- c(search_types, auto_types)
  keep <- !duplicated(comb_f)
  search_files <- comb_f[keep]
  search_types <- comb_t[keep]
}
# 对策略 1（原始 pFind 结果）与全部补充搜库结果逐行清洗：
#   Proteins 全为 REV_ 诱饵的行整行删除；正库与诱饵库并存的行仅保留正库蛋白
#   （Positions 按相同下标同步剔除）；pUniFind 转换件无 Proteins 列（转换时已按
#   pscore_is_target 过滤），原样复制留档。
# 清洗件统一写入 cfg_dir/Competitive_search_result/，并作为后续 revise 校正、
# SAP 提取与五策略合并的唯一输入——诱饵不再进入任何下游步骤，理论谱预测
# （pDeep2）与阶段 A/B 判定自然不受其影响。
# Row-level cleaning for strategy 1 (the original pFind result) and every
# supplementary search result: rows whose Proteins are all REV_ decoy are
# dropped; rows holding both target and decoy proteins keep only the target
# proteins (Positions stripped at the same subscripts); pUniFind converted
# files carry no Proteins column (already filtered by pscore_is_target at
# conversion) and are copied as-is for the record.
# The cleaned copies are written to cfg_dir/Competitive_search_result/ and are
# the sole input of the revise corrections, SAP extraction and the
# five-strategy merge - decoys never enter any downstream step, so the
# theoretical-spectrum prediction (pDeep2) and the stage A/B judgments are
# untouched by them.
csr_dir <- file.path(result_dir, "Competitive_search_result")
dir.create(csr_dir, showWarnings = FALSE, recursive = TRUE)
clean_jobs <- list(list(src = psm_file, tag = "1_original_pFind-Filtered"))
for (sf in search_files) {
  sf_dir <- tolower(basename(dirname(dirname(sf))))
  tag <- sub("\\.[^.]*$", "", paste0(sf_dir, "_", basename(sf)))
  clean_jobs[[length(clean_jobs) + 1L]] <- list(src = sf, tag = tag)
}
# v1.0（性能第二轮）：fasta 此前每个 punifind 清洗件读一遍、SAP 预处理再读
# 一遍——循环前读一次，两份变换分别复用（fa_chk = 去诱饵+toupper；
# sap_fasta = 原始 + I→L/N→GG，均与旧口径逐字节一致）
fasta_raw <- read_fasta(fasta_path)
fasta_fwd_toupper <- toupper(fasta_raw[!grepl("^(REV_|DECOY_|reversed_)", names(fasta_raw))])
cleaned <- vector("list", length(clean_jobs))
for (j in seq_along(clean_jobs)) {
  src <- clean_jobs[[j]]$src
  tag <- clean_jobs[[j]]$tag
  # pUniFind rescore 转换件（无 Proteins 列，pscore_is_target 过滤不严）：
  # 用同目录未过滤 qry.res（pFind.spectra）的蛋白注释按 File_Name+Sequence
  # 回查判定诱饵；注释缺失的行用正库蛋白组兜底（序列不在正库则删除）。
  # denovo 不适用该回查（用户规格不含 denovo，其库外序列由合并步骤回帖处理）。
  # The pUniFind rescore converted file (no Proteins column; the
  # pscore_is_target filter is loose): judge decoys by the protein annotation
  # of the unfiltered qry.res (pFind.spectra) in the same folder, matched on
  # File_Name+Sequence; annot-missing rows fall back to the forward proteome
  # (dropped when the sequence is absent from it). denovo is exempt (the
  # user's spec excludes denovo; its library-out sequences are handled by the
  # look-back at merge).
  annot <- NULL
  fa_chk <- NULL
  if (!grepl("denovo", tag, fixed = TRUE) &&
      grepl("punifind", tolower(basename(src)), fixed = TRUE)) {
    qpath <- file.path(dirname(src), "pFind.spectra")
    if (file.exists(qpath)) {
      qd <- read.delim(qpath, sep = "\t", stringsAsFactors = FALSE,
                       colClasses = "character", na.strings = character(0))
      need <- c("File_Name", "Sequence", "Proteins")
      if (all(need %in% colnames(qd))) annot <- qd[, need, drop = FALSE]
    }
    fa_chk <- fasta_fwd_toupper   # v1.0：复用循环前读入（不再每文件读一遍）
  }
  r <- clean_decoy_spectra(src, csr_dir, tag, annot = annot, fa = fa_chk)
  cleaned[[j]] <- r
  cat(sprintf("[诱饵去除] %s: %d -> %d 行（删 %d 行、剥 %d 行、涉及谱图 %d 张）\n",
              tag, r$n_in, r$n_out, r$n_dropped, r$n_stripped,
              length(r$dropped_spectra)))
  cat(sprintf("[decoy removal] %s: %d -> %d rows (dropped %d, stripped %d, %d spectra)\n",
              tag, r$n_in, r$n_out, r$n_dropped, r$n_stripped,
              length(r$dropped_spectra)))
}
pt_mark("诱饵库去除", "decoy removal")
psm_file <- cleaned[[1]]$file
clean_search_files <- vapply(cleaned[-1], function(z) z$file, character(1))
rm(cleaned, clean_jobs)

# ---- 读入 fasta 并做 SAP 预处理：I→L、N→GG（与生产版规则一致） ----
# ---- Read the fasta and apply the SAP pre-processing: I->L, N->GG (production rule) ----
sap_fasta <- gsub("N", "GG", gsub("I", "L", fasta_raw, fixed = TRUE), fixed = TRUE)
# 依次处理每个搜库结果（诱饵去除后的清洗件）：同一次读入先提取 SAP 行，
# 再执行 revise 校正并写出 _revised.spectra；模式判定与副本命名沿用原始来源路径
# For each search result (the decoy-cleaned copy): extract SAP rows from the
# single read, then run the revise corrections and write _revised.spectra; the
# mode choice and copy naming still use the original source path
sap_parts <- list()
revised_files <- character(0)
for (j in seq_along(search_files)) {
  sf <- search_files[j]
  sf_clean <- clean_search_files[j]
  cat("revising", basename(sf), "(cleaned)\n")
  # 全字符列读入：更快（不做类型转换），写出时逐字节一致
  # all-character read: faster (no type conversion) and byte-faithful on write
  dat <- read.delim(sf_clean, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = "character", na.strings = character(0))
  n0 <- nrow(dat)
  sap_r <- extract_sap_rows_df(dat, sap_fasta)
  if (!is.null(sap_r) && nrow(sap_r) > 0) {
    sap_parts[[length(sap_parts) + 1]] <- sap_r
  }
  # 执行 revise 校正，随后写出 _revised.spectra 副本
  # Apply the revise corrections, then write the _revised.spectra copy
  # 多引擎版：后处理模式由策略目录表 STRATEGY_REVISE 按 search_file_type 决定
  # （pfind_open=full、pfind_restrict/punifind_denovo/其余引擎=none、
  # punifind_rescore=terminal）；revised 副本命名沿用原始来源路径。
  # multi-engine edition: the post-processing mode comes from the strategy
  # catalog STRATEGY_REVISE keyed by search_file_type (pfind_open = full,
  # pfind_restrict/punifind_denovo/other engines = none, punifind_rescore =
  # terminal); the revised-copy naming still uses the original source path.
  sf_name <- tolower(basename(sf))
  sf_dir  <- tolower(basename(dirname(dirname(sf))))
  revise_mode <- if (search_types[j] %in% names(STRATEGY_REVISE)) {
    STRATEGY_REVISE[[search_types[j]]]
  } else "none"
  dat <- revise_search_result(dat, revise_mode)
  cat("  type:", search_types[j], "| revise mode:", revise_mode, "\n")
  rev_name <- sub("\\.spectra$", "_revised.spectra", paste0(sf_dir, "_", basename(sf)))
  write.table(dat, file.path(rev_dir, rev_name), sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE)
  cat("  rows:", n0, "->", nrow(dat), "\n")
  revised_files <- c(revised_files, file.path(rev_dir, rev_name))
}
pt_mark("revise 四个补充搜库结果", "revise the four search results")

# ---- 1b. 用同一次读入构建 SAP 索引（生产版 1.SAP_extract.R 的移植） --------------------
# ---- 1b. SAP index from the same reads (production 1.SAP_extract.R port) ---
# 合并各文件的 SAP 行并写出 SAP 索引 SAP_spectra_index.txt（UTF-8）
# Combine the SAP rows across files and write the SAP index SAP_spectra_index.txt (UTF-8)
sap_df <- bind_sap_rows(sap_parts)
sap_path <- file.path(out_dir, "SAP_spectra_index.txt")
write.table(sap_df, sap_path, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)
cat("SAP index rows:", nrow(sap_df), "->", sap_path, "\n")
# SAP 单独记录（用户 2026-09-05 规格）：索引同时留档一份到
# Competitive_search_result/——行内保留原始序列＋"pos,XXX->YYY[site]" 修饰与
# revised_peptide/ref_protein 映射，可把任何下游 SAP 替换序列还原回原始形式；
# pDeep2 的 peptide.txt 只用替换后序列，理论谱预测不受 SAP 记录影响
# Separate SAP record (user 2026-09-05 spec): the index is also archived in
# Competitive_search_result/ - each row keeps the original sequence + the
# "pos,XXX->YYY[site]" modification together with the revised_peptide /
# ref_protein mapping, so any downstream substituted SAP sequence can be
# restored to its original form; pDeep2's peptide.txt only carries the
# substituted sequences, so the SAP record never affects prediction
file.copy(sap_path, file.path(csr_dir, "SAP_spectra_index.txt"),
          overwrite = TRUE)
pt_mark("自动提取 SAP 清单", "extract SAP index")

# ---- 2-3. 合并多策略结果并写出 peptide.txt --------------------------------
# ---- 2-3. merge and write peptide.txt --------------------------------------
# 多引擎版：槽位字母由策略目录表硬编码分配（a=原始结果，b..z=补充策略），
# 合并按字母序自然排序；denovo 槽（punifind_denovo）在合并步骤执行数据库回帖。
# multi-engine edition: slot letters come from the hardcoded strategy catalog
# (a = original result, b..z = supplementary); the merge orders columns by
# letter naturally; denovo slots (punifind_denovo) get the database look-back
# at merge time.
slot_letters <- strategy_slots(search_types)
denovo_slots <- slot_letters[search_types == "punifind_denovo"]
# denovo 回帖模式（用户口径 2026-09-16）：keep＝仅保留能回帖的（默认）；
# all＝放开全部（no_match 的 denovo 序列带标记参与竞争，SAP 类变异可进仲裁）
# denovo look-back mode (user ruling 2026-09-16): keep = only matching survive
# (default); all = release everything (no-match denovo sequences compete with
# a DB_match=no_match flag, so SAP-type variants can enter the arbitration)
denovo_lb_mode <- "keep"
for (src in list(cfg2, base_cfg)) {
  if (!is.null(src$denovo_lookback) && length(src$denovo_lookback) > 0 &&
      tail(src$denovo_lookback, 1) %in% c("keep", "all")) {
    denovo_lb_mode <- tail(src$denovo_lookback, 1)
  }
}
# v1.0（无标记通道=0）：peptide.txt 剥除标记修饰后再给 pDeep2
ul_ch2 <- if (!is.null(base_cfg$unlabeled_channel) && length(base_cfg$unlabeled_channel) > 0)
  tail(base_cfg$unlabeled_channel, 1) else "1"
strip_labels2 <- if (!identical(ul_ch2, "1")) parse_label_spec(base_cfg$label_type)$mods else NULL
res <- merge_search_results(psm_file, revised_files, out_dir,
                            fasta_file = fasta_path,
                            slot_letters = slot_letters,
                            denovo_slots = denovo_slots,
                            denovo_lookback_mode = denovo_lb_mode,
                            candidate_scope = cand_scope,
                            strip_labels = strip_labels2)
pt_mark("多策略合并与 peptide.txt", "multi-strategy merge + peptide.txt")
# 策略清单指纹侧车（2026-09-26 功能排队①）：步骤 2 跳过判据从"peptide.txt
# 存在"升级为"策略集合一致"——新增/删改策略后重跑步骤 2（引擎结果经参数
# 指纹复用，仅重算合并）。user_run 用同一 strategy_identity_lines 比对。
writeLines(strategy_identity_lines(cfg2),
           file.path(out_dir, "peptide.txt.strategies"), useBytes = TRUE)
# 自动生成用户第三步参数文件（predict_file 为占位符；
# 用户已填写的 predict_file/weights 会被保留）
# auto-generate the user-step-3 parameter file (predict_file placeholder;
# existing user-filled predict_file/weights are preserved)
make_step3_config(cfg2_path, cfg2,
                  file.path(cfg2$config_dir, "config_step3.txt"))
# 打印本步骤完成信息与输出目录
# Print the completion message and output directory for this step
cat("==== [用户第二步] 完成，输出目录：", out_dir, " ====\n")
cat("==== [User Step 2] done, output: ", out_dir, " ====\n")
cat("DONE\n")
