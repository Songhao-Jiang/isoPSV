# predict_cosine_overlay.R — 多预测器余弦叠加到 master_table（v1.0）
# 在 unlabeled_master_export.R 之后运行：
#   1. 扫描 merge/predictors/<type>/predict.txt 的所有成功预测器
#   2. 读 MGF → spectrum_params → 各预测器 predict.txt → 余弦
#   3. 为 cand1 追加 cand1_predict_<type> 列（列名即软件名）
# 用法：Rscript predict_cosine_overlay.R <result_dir> <mgf_folder>
# 输出：master_table.tsv 新增列，原列不动
CODE <- NULL
for (.p in c(
  file.path(getwd(), "code"),
  file.path(getwd(), "..", "code"),
  getwd())) {
  if (file.exists(file.path(.p, "ReadMGF_script.R"))) { CODE <- .p; break }
}
if (is.null(CODE)) stop("ReadMGF_script.R not found from getwd: ", getwd())
# （工作站 170 报告新发现）：独立 CLI 调用时的 LC_CTYPE 兜底——
# 继承 LC_CTYPE=C.UTF-8 的环境里 source(..., encoding="UTF-8") 对中文注释
# 文件报 invalid input、ReadMGF 只部分载入、后续 "could not find function"
# 炸。探针对
# "locale 名义 utf8 但 CRT 实际坏" 的场景不触发，故加哨兵验证：源后若
# 关键函数缺席，强制切换候选 locale 重源一遍；仍失败则以可操作信息终止。
# master 启动器场景本就受 user_run.R 探针保护，本块只服务独立调用。
#  (workstation 170-round finding): LC_CTYPE fallback for standalone
# CLI invocation - under an inherited LC_CTYPE=C.UTF-8, source(...,
# encoding="UTF-8") fails with "invalid input" on Chinese-commented files,
# ReadMGF loads partially and the run later dies with "could not find
# function" (hit by both the workstation CLI test method and the dev-side
# local E2E; clearing LC_* heals it). A name-only probe cannot catch a locale
# that claims utf8 while the CRT is broken, so a sentinel check follows the
# sources: if key functions are missing, force-switch through candidate
# locales and re-source; abort with an actionable message if still broken.
# The master launcher already probes (user_run.R); this block only serves
# standalone calls.
.isopsv_utf8_locale <- function()
  grepl("utf-?8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)
if (!.isopsv_utf8_locale()) {
  for (.lc in c("English_United States.UTF-8", "C.UTF-8")) {
    try(Sys.setlocale("LC_CTYPE", .lc), silent = TRUE)
    if (.isopsv_utf8_locale()) break
  }
}
.isopsv_sources <- function() {
  source(file.path(CODE, "compat.R"), encoding = "UTF-8")
  source(file.path(CODE, "ReadMGF_script.R"), encoding = "UTF-8")
  source(file.path(CODE, "paired_params.R"), encoding = "UTF-8")
}
try(.isopsv_sources(), silent = TRUE)
if (!exists("parse_lines_vectorized", inherits = FALSE) ||
    !exists("read_predict_txt", inherits = FALSE)) {
  for (.lc in c("English_United States.UTF-8", "C.UTF-8",
                "Chinese (Simplified)_China.936", "")) {
    try(Sys.setlocale("LC_CTYPE", .lc), silent = TRUE)
    try(.isopsv_sources(), silent = TRUE)
    if (exists("parse_lines_vectorized", inherits = FALSE) &&
        exists("read_predict_txt", inherits = FALSE)) break
  }
  if (!exists("parse_lines_vectorized", inherits = FALSE) ||
      !exists("read_predict_txt", inherits = FALSE)) {
    stop("UTF-8 sources failed to load under LC_CTYPE=",
         Sys.getlocale("LC_CTYPE"),
         " - clear LC_* before standalone invocation (or run under the master)")
  }
}

args <- commandArgs(FALSE)
arg_of <- function(flag) {
  hits <- args[startsWith(args, paste0("--", flag, "="))]
  if (length(hits) > 0) sub(paste0("^--", flag, "="), "", hits[1]) else ""
}
rd <- arg_of("result_dir")
if (!nzchar(rd)) stop("usage: Rscript predict_cosine_overlay.R --result_dir=<dir>")

aa <- read.delim(file.path(CODE, "AA_mass_dataframe.txt"), sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- data.frame(
    Full_Name = c("Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
                  "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
                  "Leucine","Lysine","Methionine","Phenylyalanine","Proline",
                  "Serine","Threonine","Tryptophan","Tyrosine","Valine","Ile_or_Leu"),
    Three_Letter = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
                     "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp",
                     "Tyr","Val","Xle"),
    One_Letter = c("A","R","N","D","C","Q","E","G","H","I","L","K","M",
                   "F","P","S","T","W","Y","V","X"))

mt_path <- file.path(rd, "1.unlabeled_result", "master_table.tsv")
if (!file.exists(mt_path)) stop("master_table.tsv not found: ", mt_path)
mt <- read.delim(mt_path, sep = "\t", check.names = FALSE,
                 stringsAsFactors = FALSE, na.strings = character(0),
                 colClasses = "character")
cat("[predict_cosine_overlay] master_table: ", nrow(mt), " 行\n", sep = "")

# 扫描预测器
pred_dirs <- list.dirs(file.path(rd, "merge", "predictors"), recursive = FALSE)
pred_types <- character(0)
for (pd in pred_dirs) {
  ptxt <- file.path(pd, "predict.txt")
  if (file.exists(ptxt)) pred_types <- c(pred_types, basename(pd))
}
# pDeep2 产物在 merge/predict.txt（旧布局）或 predictors/pdeep2/（新布局）
pd2_new <- file.path(rd, "merge", "predictors", "pdeep2", "predict.txt")
pd2_old <- file.path(rd, "merge", "predict.txt")
if (file.exists(pd2_old) && !("pdeep2" %in% pred_types)) pred_types <- c("pdeep2", pred_types)
cat("[predict_cosine_overlay] 预测器: ", paste(pred_types, collapse = ", "), "\n", sep = "")
if (length(pred_types) == 0) { cat("无预测器产物，跳过\n"); quit(status = 0) }

# 读各预测器 predict.txt
preds <- list()
for (pt in pred_types) {
  pf <- if (pt == "pdeep2" && file.exists(pd2_old)) pd2_old else
        file.path(rd, "merge", "predictors", pt, "predict.txt")
  pr <- try(read_predict_txt(pf), silent = TRUE)
  if (!inherits(pr, "try-error") && !is.null(pr$idx)) {
    preds[[pt]] <- pr
    cat("  ", pt, ": ", length(pr$int), " entries\n", sep = "")
  }
}

# MGF 位置：rd 可能是 result_* 时间戳子目录——先试 rd/unlabeled，再试
# dirname(rd)/unlabeled（项目根），最后读 config_step1 的 mgf 键
mgf_dir <- file.path(rd, "unlabeled")
mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
if (length(mgf_files) == 0) {
  mgf_dir <- file.path(dirname(rd), "unlabeled")
  mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
}
if (length(mgf_files) == 0) {
  # 从 config_step1 读 mgf 键
  cfg1_path <- file.path(dirname(rd), "config_step1.txt")
  if (file.exists(cfg1_path)) {
    cfg1 <- read_param_file(cfg1_path)
    mgf_key <- if (!is.null(cfg1$mgf)) tail(cfg1$mgf, 1) else ""
    if (nzchar(mgf_key) && dir.exists(mgf_key))
      mgf_files <- sort(list.files(mgf_key, pattern = "\\.mgf$", full.names = TRUE))
  }
}
cat("[predict_cosine_overlay] MGF: ", length(mgf_files), " 个\n", sep = "")
if (length(mgf_files) == 0) {
  cat("[predict_cosine_overlay] 未找到 MGF，跳过余弦计算（master_table 已含可用性列）\n")
  cat("[predict_cosine_overlay] no MGF found, skipping (master_table already has availability columns)\n")
  quit(status = 0)
}

# 建 File_Name → MGF spectrum 映射
spec_map <- new.env(parent = emptyenv())
for (mf in mgf_files) {
  mgfs <- parse_lines_vectorized(readLines(mf, warn = FALSE))
  for (s in mgfs) assign(s@title, s, envir = spec_map)
  rm(mgfs); gc()
}
cat("[predict_cosine_overlay] 谱图索引: ", ls(spec_map), " 条\n", sep = "")

# ----：行主序重构＋心跳（工作站 run8 实测单列 30-65 分钟）----------
# 旧列主序在每行×每列重复调 spectrum_params——同一行的实验 b/y 离子强度
# 对全部预测器列是同一个向量，却被算了一遍又一遍（4 列＝4 遍；65 分钟
# 零日志单核 grind 的病根）。行主序：每行提取一次离子强度向量，内层
# 只做各预测器的键查＋余弦（O(1) 查表）；总 spectrum_params 次数＝行数，
# 与列数无关。每 2,000 行打心跳，列开始前打印列清单。
# ----: row-major restructure + heartbeat (workstation run8
# measured 30-65 min PER COLUMN). The old column-major order called
# spectrum_params once per row per column - the same row's experimental
# b/y ion intensities are one vector shared by every predictor column but
# were recomputed for each (4 columns = 4 passes; the root of the 65-minute
# zero-log single-core grind). Row-major: extract the intensity vector once
# per row, and the inner loop only does per-predictor key lookups and
# cosines (O(1)); total spectrum_params calls = row count, independent of
# the column count. Heartbeat every 2,000 rows; the column list prints up
# front.
new_cols <- paste0("cand1_predict_", names(preds))
new_cols <- setdiff(new_cols, names(mt))
if (length(new_cols) == 0) {
  cat("[predict_cosine_overlay] 无新增列\n")
  quit(status = 0)
}
cat("[predict_cosine_overlay] 计算列: ", paste(new_cols, collapse = ", "),
    "\n", sep = "")
vals <- setNames(rep(list(rep("", nrow(mt))), length(new_cols)), new_cols)
n_hit <- setNames(integer(length(new_cols)), new_cols)
pr_list <- preds[sub("^cand1_predict_", "", new_cols)]
t0 <- proc.time()[3]
n_spec <- 0L
for (i in seq_len(nrow(mt))) {
  seq_p <- mt$cand1_sequence[i]
  if (is.na(seq_p) || !nzchar(seq_p)) next
  mod_p <- sub(";+$", "", as.character(mt$cand1_modification[i]))
  if (is.na(mod_p)) mod_p <- ""
  # 查谱图＋实验离子强度——每行一次，全部预测器列共用
  # spectrum lookup + experimental ion intensities - once per row, shared
  spc <- get0(mt$File_Name[i], envir = spec_map, inherits = FALSE,
              ifnotfound = NULL)
  if (is.null(spc)) next
  sp <- try(spectrum_params(seq_p, mod_p, "unlabeled", aa, aa_tab,
                            spc@spectrum, 20, "ppm", labels = NULL),
            silent = TRUE)
  if (is.null(sp) || inherits(sp, "try-error")) next
  v <- sp$b_y_ions_intensity
  if (length(v) == 0) next
  n_spec <- n_spec + 1L
  key1 <- paste0(seq_p, "|", mod_p, "|", mt$Charge[i])
  key2 <- paste0(seq_p, "||", mt$Charge[i])
  for (j in seq_along(new_cols)) {
    pr <- pr_list[[j]]
    pi <- unname(pr$idx[key1])
    if (is.na(pi)) pi <- unname(pr$idx[key2])
    if (is.na(pi)) next
    pv <- pr$int[[pi]]
    nmin <- min(length(pv), length(v))
    if (nmin <= 0) next
    vals[[j]][i] <- sprintf("%.4f", cosine(v[seq_len(nmin)], pv[seq_len(nmin)]))
    n_hit[j] <- n_hit[j] + 1L
  }
  if (i %% 2000L == 0L) {
    cat(sprintf("  [overlay] 行 %d/%d（spectrum_params %d 次，已用 %.0fs）\n",
                i, nrow(mt), n_spec, proc.time()[3] - t0))
  }
}
for (j in seq_along(new_cols)) mt[[new_cols[j]]] <- vals[[j]]
n_new_cols <- length(new_cols)
for (j in seq_along(new_cols)) {
  cat("  ", new_cols[j], ": ", n_hit[j], "/", nrow(mt), " 覆盖\n", sep = "")
}
cat(sprintf("[predict_cosine_overlay] 行主序总耗时 %.0fs（spectrum_params %d 次＝行级一次，与列数无关）\n",
            proc.time()[3] - t0, n_spec))

# 写回
if (n_new_cols > 0) {
  write.table(mt, mt_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[predict_cosine_overlay] 已写回 ", n_new_cols, " 列 -> ", mt_path, "\n", sep = "")
} else {
  cat("[predict_cosine_overlay] 无新增列\n")
}
