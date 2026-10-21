# =============================================================================
# pure_prediction.R — 纯预测模式
# pure_prediction.R - pure prediction mode 
#
# 触发：竞争性搜库槽位一个都不选（config_step2 的 pure_prediction=1）。
# 本脚本两段式运行（由 user_run.R 的 run_step 分别调用）：
#   stage=peptides  基础鉴定表（config_step1 的 psm_file）→ merge/peptide.txt
#                   （3 列 peptide/modification/charge，去重；min_pep_len 与
#                   无标记通道的标记修饰剥除口径与 user_step2 一致），并写
#                   peptide.txt.strategies 侧车（pure_prediction）与
#                   merge/pure_input.path（导出段回指输入表）。
#   stage=export    输入表原样 + 每个预测器一列 predict_cosine_<type>
#                   （观测谱 vs 理论谱余弦，键 seq|mod|z → seq||z 回退，
#                   与 predict_cosine_overlay.R 同一实现口径）→
#                   <result_dir>/pure_prediction.tsv。
#
# Trigger: no competitive-search slot selected (pure_prediction=1 in
# config_step2). Two stages, each invoked by user_run.R via run_step:
#   stage=peptides  base identification table (psm_file of config_step1) ->
#                   merge/peptide.txt (3 columns peptide/modification/charge,
#                   deduplicated; min_pep_len and the label-mod stripping for
#                   the unlabeled channel follow user_step2), plus the
#                   peptide.txt.strategies sidecar (pure_prediction) and
#                   merge/pure_input.path pointing back at the input table.
#   stage=export    the input table as-is + one predict_cosine_<type> column
#                   per predictor (observed-vs-predicted cosine; key
#                   seq|mod|z with a seq||z fallback, the same implementation
#                   convention as predict_cosine_overlay.R) ->
#                   <result_dir>/pure_prediction.tsv.
#
# 用法 / Usage:
#   Rscript pure_prediction.R   （环境变量 STAGE / RESULT_DIR / CFG1）
#   env vars: STAGE = peptides | export, RESULT_DIR = result folder,
#             CFG1 = config_step1.txt path
# =============================================================================
args0 <- commandArgs(FALSE)
farg <- sub("^--file=", "", args0[grepl("^--file=", args0)][1])
if (is.na(farg)) farg <- "pure_prediction.R"
CODE <- normalizePath(dirname(farg), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(CODE, "ReadMGF_script.R"))) {
  d2 <- dirname(CODE)
  if (file.exists(file.path(d2, "code", "ReadMGF_script.R"))) CODE <- file.path(d2, "code")
}
suppressMessages(source(file.path(CODE, "compat.R"), encoding = "UTF-8"))
suppressMessages(source(file.path(CODE, "ReadMGF_script.R"), encoding = "UTF-8"))
suppressMessages(source(file.path(CODE, "paired_params.R"), encoding = "UTF-8"))

STAGE <- toupper(Sys.getenv("STAGE", unset = ""))
RESULT_DIR <- Sys.getenv("RESULT_DIR", unset = "")
CFG1 <- Sys.getenv("CFG1", unset = "")
if (!nzchar(STAGE) || !nzchar(RESULT_DIR)) {
  stop("pure_prediction.R: env STAGE (peptides|export) and RESULT_DIR are required")
}
if (!nzchar(CFG1)) CFG1 <- file.path(dirname(RESULT_DIR), "config_step1.txt")
if (!file.exists(CFG1)) stop("config_step1 not found: ", CFG1)
base_cfg <- read_param_file(CFG1)
aa <- read.delim(file.path(CODE, "AA_mass_dataframe.txt"), sep = "\t",
                  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# 输入表：与 user_step2 同源（cfg1 的 psm_file 首个有效值）
psm_raw <- resolve_cfg_paths(base_cfg, "psm_file")
psm_raw <- psm_raw[!is.na(psm_raw) & psm_raw != ""]
if (length(psm_raw) == 0) stop("psm_file not set in config_step1")
psm_file <- psm_raw[1]
if (!file.exists(psm_file)) stop("psm_file does not exist: ", psm_file)

read_psm <- function(path) {
  # 列类型守卫：与 merge_search_results 同款（Modification 全空列防 logical 化）
  hdr0 <- strsplit(readLines(path, n = 1, warn = FALSE), "\t")[[1]]
  cc0 <- rep(NA_character_, length(hdr0))
  for (col0 in c("File_Name", "Sequence", "Modification")) {
    i0 <- which(hdr0 == col0)
    if (length(i0) == 1) cc0[i0] <- "character"
  }
  d <- read.delim(path, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, colClasses = cc0)
  if ("Modification" %in% names(d)) {
    d$Modification[is.na(d$Modification)] <- ""
    d$Modification[d$Modification == "NA"] <- ""
  }
  d
}

if (STAGE == "PEPTIDES") {
  cat("==== [纯预测/peptides] 输入表：", psm_file, " ====\n", sep = "")
  cat("==== [pure-prediction/peptides] input table: ", psm_file, " ====\n", sep = "")
  dat <- read_psm(psm_file)
  cat("[纯预测] 输入行数：", nrow(dat), "\n", sep = "")
  cat("[pure-prediction] input rows: ", nrow(dat), "\n", sep = "")
  seqv <- as.character(dat$Sequence)
  modv <- if ("Modification" %in% names(dat)) as.character(dat$Modification) else rep("", nrow(dat))
  chv <- suppressWarnings(as.integer(dat$Charge))
  # 长度阈值（config_step1 键 min_pep_len，缺省 7）：短肽预测软件做不了
  # length gate (config_step1 key min_pep_len, default 7): predictors cannot
  # handle very short peptides
  mpl <- 7L
  if (!is.null(base_cfg$min_pep_len) && length(base_cfg$min_pep_len) > 0) {
    v <- suppressWarnings(as.integer(tail(base_cfg$min_pep_len, 1)))
    if (!is.na(v) && v >= 4L) mpl <- v
  }
  keep <- nzchar(seqv) & nchar(seqv) >= mpl & !is.na(chv) & chv >= 1
  cat("[纯预测] 长度阈值 min_pep_len=", mpl, "，进入预测的行 ", sum(keep),
      "/", nrow(dat), "\n", sep = "")
  cat("[pure-prediction] length gate min_pep_len=", mpl, ", rows entering prediction ",
      sum(keep), "/", nrow(dat), "\n", sep = "")
  seqv <- seqv[keep]; modv <- modv[keep]; chv <- chv[keep]
  # 无标记通道≠1 时剥除标记修饰（与 user_step2 的 peptide.txt 口径一致）
  ul_ch <- if (!is.null(base_cfg$unlabeled_channel) && length(base_cfg$unlabeled_channel) > 0)
    tail(base_cfg$unlabeled_channel, 1) else "1"
  if (!identical(ul_ch, "1")) {
    strip_mods <- parse_label_spec(base_cfg$label_type)$mods
    if (length(strip_mods) > 0) {
      modv <- strip_label_mods(modv, strip_mods)
      cat("[纯预测] 无标记通道=", ul_ch, "，已剥除标记修饰\n", sep = "")
      cat("[pure-prediction] unlabeled channel=", ul_ch, ", label mods stripped\n", sep = "")
    }
  }
  # 去重成三元组（peptide/modification/charge）
  trip <- unique(data.frame(peptide = seqv, modification = modv, charge = chv,
                            stringsAsFactors = FALSE))
  out_dir <- file.path(RESULT_DIR, "merge")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_pep <- file.path(out_dir, "peptide.txt")
  write.table(trip, out_pep, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  writeLines("pure_prediction", file.path(out_dir, "peptide.txt.strategies"), useBytes = TRUE)
  writeLines(psm_file, file.path(out_dir, "pure_input.path"), useBytes = TRUE)
  cat("[纯预测] peptide.txt：", nrow(trip), " 个肽段形式 -> ", out_pep, "\n", sep = "")
  cat("[pure-prediction] peptide.txt: ", nrow(trip), " peptide forms -> ", out_pep, "\n", sep = "")
  quit(status = 0)
}

if (STAGE == "EXPORT") {
  inp_path <- file.path(RESULT_DIR, "merge", "pure_input.path")
  in_file <- if (file.exists(inp_path)) readLines(inp_path, warn = FALSE)[1] else psm_file
  if (!file.exists(in_file)) in_file <- psm_file
  dat <- read_psm(in_file)
  cat("[纯预测/export] 输入表：", in_file, "（", nrow(dat), " 行）\n", sep = "")
  cat("[pure-prediction/export] input table: ", in_file, " (", nrow(dat), " rows)\n", sep = "")
  # 预测器产物发现（与 predict_cosine_overlay.R 同一布局规则）
  pred_types <- character(0)
  pd_root <- file.path(RESULT_DIR, "merge", "predictors")
  if (dir.exists(pd_root)) {
    for (pd in list.dirs(pd_root, recursive = FALSE)) {
      if (file.exists(file.path(pd, "predict.txt")))
        pred_types <- c(pred_types, basename(pd))
    }
  }
  pd2_old <- file.path(RESULT_DIR, "merge", "predict.txt")
  if (file.exists(pd2_old) && !("pdeep2" %in% pred_types)) pred_types <- c("pdeep2", pred_types)
  cat("[纯预测/export] 预测器：", paste(pred_types, collapse = ", "), "\n", sep = "")
  cat("[pure-prediction/export] predictors: ", paste(pred_types, collapse = ", "), "\n", sep = "")
  if (length(pred_types) == 0) {
    cat("[纯预测/export] 无任何预测器产物——输出仅原样输入表\n")
    cat("[pure-prediction/export] no predictor output - writing the input table unchanged\n")
  }
  preds <- list()
  for (pt in pred_types) {
    pf <- if (pt == "pdeep2" && file.exists(pd2_old)) pd2_old else
      file.path(pd_root, pt, "predict.txt")
    pr <- try(read_predict_txt(pf), silent = TRUE)
    if (!inherits(pr, "try-error") && !is.null(pr$idx)) {
      preds[[pt]] <- pr
      cat("  ", pt, ": ", length(pr$int), " entries\n", sep = "")
    }
  }
  # MGF 谱图索引
  mgf_dir <- file.path(RESULT_DIR, "unlabeled")
  mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
  if (length(mgf_files) == 0) {
    mgf_dir <- file.path(dirname(RESULT_DIR), "unlabeled")
    mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
  }
  if (length(mgf_files) == 0) {
    mgf_key <- if (!is.null(base_cfg$mgf)) tail(base_cfg$mgf, 1) else ""
    if (nzchar(mgf_key) && dir.exists(mgf_key))
      mgf_files <- sort(list.files(mgf_key, pattern = "\\.mgf$", full.names = TRUE))
  }
  cat("[纯预测/export] MGF：", length(mgf_files), " 个\n", sep = "")
  cat("[pure-prediction/export] MGF files: ", length(mgf_files), "\n", sep = "")
  spec_map <- new.env(parent = emptyenv())
  if (length(mgf_files) > 0) {
    for (mf in mgf_files) {
      mgfs <- parse_lines_vectorized(readLines(mf, warn = FALSE))
      for (s in mgfs) assign(s@title, s, envir = spec_map)
      rm(mgfs); gc()
    }
  }
  seqv <- as.character(dat$Sequence)
  modv <- if ("Modification" %in% names(dat)) as.character(dat$Modification) else rep("", nrow(dat))
  modv[is.na(modv)] <- ""
  chv <- suppressWarnings(as.integer(dat$Charge))
  fnv <- as.character(dat$File_Name)
  for (pt in names(preds)) {
    col_name <- paste0("predict_cosine_", pt)
    vals <- rep("", nrow(dat))
    pr <- preds[[pt]]
    n_hit <- 0L
    if (length(mgf_files) > 0) {
      for (i in seq_len(nrow(dat))) {
        sq <- seqv[i]
        if (!nzchar(sq) || is.na(sq) || is.na(chv[i])) next
        md <- sub(";+$", "", modv[i])
        if (is.na(md)) md <- ""
        pi <- unname(pr$idx[paste0(sq, "|", md, "|", chv[i])])
        if (is.na(pi)) pi <- unname(pr$idx[paste0(sq, "||", chv[i])])
        #：pDeep2 的 predict.py 按 peptide.txt 原样回写 pepinfo，
        # 修饰串带尾分号（"…8,Oxidation[M];|2"）——补带尾分号的第三键再查
        #: pDeep2's predict.py echoes pepinfo verbatim, keeping the
        # trailing semicolon - try that third key before giving up
        if (is.na(pi) && nzchar(modv[i]))
          pi <- unname(pr$idx[paste0(sq, "|", modv[i], "|", chv[i])])
        if (is.na(pi)) next
        #：get() 未命中会中止整段——改 tryCatch 计数后继续
        # （真实数据里 File_Name 与 MGF title 错位时不再全表崩掉）
        #: a get() miss aborted the whole stage - catch, count, continue
        # (a File_Name / MGF-title mismatch on real data no longer crashes out)
        spc <- tryCatch(get(fnv[i], envir = spec_map, inherits = FALSE),
                        error = function(e) NULL)
        if (is.null(spc)) next
        sp <- try(spectrum_params(sq, md, "unlabeled", aa, amino_acid_table,
                                  spc@spectrum, 20, "ppm", labels = NULL),
                  silent = TRUE)
        if (is.null(sp) || inherits(sp, "try-error")) next
        v <- sp$b_y_ions_intensity
        pv <- pr$int[[pi]]
        nmin <- min(length(pv), length(v))
        if (nmin <= 0) next
        vals[i] <- sprintf("%.4f", cosine(v[seq_len(nmin)], pv[seq_len(nmin)]))
        n_hit <- n_hit + 1L
      }
    }
    dat[[col_name]] <- vals
    cat("  ", col_name, ": ", n_hit, "/", nrow(dat), " 覆盖\n", sep = "")
    cat("  ", col_name, ": ", n_hit, "/", nrow(dat), " covered\n", sep = "")
  }
  out_tsv <- file.path(RESULT_DIR, "pure_prediction.tsv")
  write.table(dat, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("==== [纯预测/export] 输出：", out_tsv, " ====\n", sep = "")
  cat("==== [pure-prediction/export] output: ", out_tsv, " ====\n", sep = "")
  quit(status = 0)
}

stop("unknown STAGE: ", STAGE)
