# =============================================================================
# unlabeled_master_export.R —— 无标记流程全谱主表（v1.0，2026-10-05）
#
# 背景：94/96 列 master_table.tsv 历来仅标记样本（依赖成对报告 unlabel_
# ions_match.* 列）；无标记样本按 2026-09-06 口径跳过。用户口径 2026-10-05：
# 无标记也要主表。本脚本从阶段 A 的 scan_details_A.tsv（逐谱×逐候选长表，
# 14 列）透视出与标记主表同风格的宽表：
#   - 每谱一行；cand1..cand5 装按 is_best/顺序排定的候选（sequence/
#     modification/coverage/PSM_score_by/PSM_score_all/predict_cosine/
#     specific_ratio/ACS_A/decision_A）；
#   - 标记侧列（lb_*、pcos、paired_state）一律留空；
#   - 判定列并入 classification_A.tsv 的 conclusion_status/vote_says_wrong/
#     error_reason/ambiguous_residue/chimeric（有则并入）。
# 输出：<result>/1.unlabeled_result/master_table.tsv
# 用法（由 user_run.R 无标记分支调用；也可独立）：
#   CONFIG3=<config_step3> ALL_MGF=<mgf 或空> Rscript unlabeled_master_export.R
# =============================================================================

args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
if (!dir.exists(script_dir)) {
  script_dir <- dirname(if (grepl("^[A-Za-z]:", script_arg)) script_arg else file.path(getwd(), script_arg))
  script_dir <- gsub("\\\\", "/", script_dir)
}
old_wd <- setwd(script_dir)
source("compat.R")
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x[1]

read_kv <- function(p) {
  if (!file.exists(p)) return(list())
  kv <- readLines(p, warn = FALSE)
  kv <- kv[!startsWith(kv, "#") & nzchar(kv)]
  keys <- sub("\t.*$", "", kv); vals <- sub("^\t?", "", sub("^[^\t]*\t", "", kv))
  vals <- sub("^\t", "", vals)
  as.list(stats::setNames(vals, keys))
}
resolve_first <- function(cfg, base_dir, key) {
  v <- cfg[[key]][1]; if (is.null(v) || !nzchar(v)) return("")
  v <- gsub("\\\\", "/", v)
  if (grepl("^[A-Za-z]:", v)) v else file.path(base_dir, v)
}

cfg3_path <- Sys.getenv("CONFIG3")
if (!nzchar(cfg3_path)) stop("CONFIG3 not set")
cfg3 <- read_kv(cfg3_path)
base_dir <- dirname(cfg3_path)
bc <- resolve_first(cfg3, base_dir, "base_config")
cfg1 <- if (nzchar(bc)) read_kv(bc) else list()
merge_result <- resolve_first(cfg3, base_dir, "merge_result")
# v100d：结果目录解析——四级回退（pipeline 会在 result_dir 下建
# result_时间戳/ 子目录，config_step1 的 result_dir 只是根，不是结果本体）
# ① merge_result 上溯 2 级（最精确：…/result_YYY/merge/… → …/result_YYY/）
# ② cfg1$result_dir（可能是根而非结果本体，需验 1.unlabeled_result 存在）
# ③ cfg1$result_dir 下扫最新的 result_* 含 1.unlabeled_result 的目录
# ④ config_step3 自身 result_dir 键
find_unl_dir <- function(base) {
  d <- file.path(base, "1.unlabeled_result")
  if (dir.exists(d)) return(d)
  NULL
}
unl_dir <- NULL
# ① merge_result 上溯
if (nzchar(merge_result)) {
  cand <- dirname(dirname(merge_result))
  unl_dir <- find_unl_dir(cand)
  if (!is.null(unl_dir)) cat("[无标记主表] 结果目录（merge_result 推导）：", cand, "\n", sep = "")
}
# ② cfg1$result_dir 直下
if (is.null(unl_dir) && !is.null(cfg1$result_dir)) {
  rd <- gsub("\\\\", "/", trimws(cfg1$result_dir))
  unl_dir <- find_unl_dir(rd)
}
# ③ result_dir 下最新的 result_* 含 1.unlabeled_result
if (is.null(unl_dir) && !is.null(cfg1$result_dir)) {
  rd <- gsub("\\\\", "/", trimws(cfg1$result_dir))
  if (dir.exists(rd)) {
    subs <- list.dirs(rd, recursive = FALSE)
    subs <- subs[grepl("^result_", basename(subs))]
    if (length(subs) > 0) {
      subs <- subs[order(file.info(subs)$mtime, decreasing = TRUE)]
      for (s in subs) {
        unl_dir <- find_unl_dir(s)
        if (!is.null(unl_dir)) {
          cat("[无标记主表] 结果目录（result_* 时间戳扫描）：", s, "\n", sep = "")
          break
        }
      }
    }
  }
}
# ④ config_step3 自身
if (is.null(unl_dir)) {
  rd3 <- resolve_first(cfg3, base_dir, "result_dir")
  if (nzchar(rd3)) unl_dir <- find_unl_dir(rd3)
}
if (is.null(unl_dir))
  stop("scan_details_A.tsv not found: no 1.unlabeled_result under result_dir or result_* subfolders (CONFIG3: ",
       cfg3_path, ")")
sd_path <- file.path(unl_dir, "scan_details_A.tsv")
if (!file.exists(sd_path)) stop("scan_details_A.tsv not found: ", sd_path)

cat("[无标记主表] 读取 ", sd_path, "\n", sep = "")
sd <- read.delim(sd_path, sep = "\t", check.names = FALSE,
                 stringsAsFactors = FALSE, na.strings = character(0),
                 colClasses = "character")
sd$File_Name <- as.character(sd$File_Name)
sd$is_best <- as.character(sd$is_best)

# classification_A 并入（有则取判定列）
# v1.0：report_A 的 candidate_predict_cosine 覆盖 cand1 的 predict 列
# （patch_predict_cols.R 只重算 report_A/classification_A，不动 scan_details）
ra_path <- file.path(unl_dir, "report_A.tsv")
ra_pred <- NULL
if (file.exists(ra_path)) {
  ra <- read.delim(ra_path, sep = "	", check.names = FALSE,
                   stringsAsFactors = FALSE, na.strings = character(0),
                   colClasses = "character")
  if ("candidate_predict_cosine" %in% names(ra))
    ra_pred <- stats::setNames(as.character(ra$candidate_predict_cosine),
                               as.character(ra$File_Name))
}
cls_path <- file.path(unl_dir, "classification_A.tsv")
cls <- NULL
if (file.exists(cls_path)) {
  cls <- read.delim(cls_path, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = character(0),
                    colClasses = "character")
  cls$File_Name <- as.character(cls$File_Name)
}

want_cols <- c("conclusion_status", "vote_says_wrong", "error_reason",
                "ambiguous_residue", "chimeric", "best_sequence", "decision_A",
                "best_group")
CAND_COLS <- c("sequence", "modification", "coverage", "PSM_score_by",
               "PSM_score_all", "predict_cosine", "specific_ratio",
               "ACS_A", "decision_A")
gen_cand_cols <- function(k) paste0("cand", k, "_",
  c("sequence", "modification", "ul_coverage", "ul_by", "ul_all",
    "predict", "specific_ratio", "ACS_A", "decision_A"))
out_cols <- c("File_Name", "Scan_No", "Charge",
              unlist(lapply(1:5, gen_cand_cols)),
              "n_candidates", "best_cand_idx",
              paste0("cls_", want_cols))

split_spec <- function(v) strsplit(as.character(v), "\\|")[[1]]
# 行序＝长表首现顺序（候选已按打分降序产出；is_best=1 的行排 cand1 保障位）
by_spec <- split(sd, sd$File_Name, drop = TRUE)
rows <- vector("list", length(by_spec))
cls_idx <- if (!is.null(cls)) stats::setNames(seq_len(nrow(cls)),
                    make.unique(as.character(cls$File_Name))) else NULL
  cr <- NULL
for (i in seq_along(by_spec)) {
  g <- by_spec[[i]]
  fn <- names(by_spec)[i]
  scan_no <- as.character(g$Scan_No[1]) %||% ""
  charge  <- as.character(g$Charge[1]) %||% ""
  # 候选排序：is_best==1 优先，其余按原序
  ord <- order(!(g$is_best == "1"), seq_len(nrow(g)))
  g <- g[ord, , drop = FALSE]
  row <- c(fn, scan_no, charge)
  for (k in 1:5) {
    if (k <= nrow(g)) {
      for (j in seq_along(CAND_COLS))
        row <- c(row, as.character(g[[CAND_COLS[j]]][k]) %||% "")
    } else row <- c(row, rep("", length(CAND_COLS)))
  }
  row <- c(row, nrow(g), if (any(g$is_best == "1")) which(g$is_best == "1")[1] else 1)
  # cand1_predict 用 report_A 的重算值覆盖（scan_details 的该列未被 patch）
  if (!is.null(ra_pred) && fn %in% names(ra_pred) &&
      nzchar(ra_pred[[fn]]) && !is.na(ra_pred[[fn]])) {
    idx <- 3 + 6                       # cand1_predict 在行向量中的位置
    row[[idx]] <- ra_pred[[fn]]
  }
  cr <- NULL
  if (!is.null(cls_idx) && fn %in% names(cls_idx) && !is.na(cls_idx[[fn]])) {
    cr <- cls[cls_idx[[fn]], , drop = FALSE]
  }
  for (wc in want_cols)
    row <- c(row, if (!is.null(cr) && wc %in% names(cr) &&
                     !is.na(cr[[wc]][1])) as.character(cr[[wc]][1]) else "")
  rows[[i]] <- row
}
out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(out) <- out_cols

# ---- v1.0（李治反馈）：多预测器分列——从 merge/predictors/<type>/ 读取
# 各预测器产物，为 cand1 加 predict 可用性列（后续版本可扩展为各自余弦值）
# ⚠ v1.0 修复：不能用 result_dir（此脚本未定义该变量，同类 bug 第 5 次）；
# unl_dir = <结果根>/1.unlabeled_result，merge 在结果根下
pred_dirs <- list.dirs(file.path(dirname(dirname(unl_dir)), "merge", "predictors"),
                       recursive = FALSE)
# 回退旧布局（merge/predictors 在 unl_dir 同级）
if (length(pred_dirs) == 0) {
  alt <- file.path(dirname(unl_dir), "merge", "predictors")
  if (dir.exists(alt)) pred_dirs <- list.dirs(alt, recursive = FALSE)
}
if (length(pred_dirs) > 0) {
  for (pd in pred_dirs) {
    ptype <- basename(pd)
    ptxt <- file.path(pd, "predict.txt")
    if (!file.exists(ptxt)) next
    # 读 TITLE 行建键集（seq|mod|charge）
    pkeys <- character(0)
    con <- file(ptxt, "r", encoding = "UTF-8")
    repeat {
      ln <- readLines(con, n = 1)
      if (length(ln) == 0) break
      if (startsWith(ln, "TITLE=")) pkeys <- c(pkeys, substr(ln, 7, nchar(ln)))
    }
    close(con)
    if (length(pkeys) == 0) next
    pset <- unique(sub("\\|[^|]*$", "", pkeys))  # seq|mod（去 charge 便于匹配）
    # 为每个 cand1 加可用性
    col_nm <- paste0("cand1_has_predict_", ptype)
    out[[col_nm]] <- vapply(out$cand1_sequence, function(s) {
      if (!nzchar(as.character(s))) return("")
      hit <- any(pset == as.character(s) | pset == paste0(as.character(s), "|"))
      if (hit) "1" else "0"
    }, character(1))
    cat("[无标记主表] 多预测器列：", ptype, " ", sum(out[[col_nm]] == "1"),
        "/", nrow(out), " 覆盖\n", sep = "")
  }
}

out_path <- file.path(unl_dir, "master_table.tsv")
write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[无标记主表] 完成：", nrow(out), " 行 × ", ncol(out), " 列 -> ", out_path, "\n", sep = "")
setwd(old_wd)
