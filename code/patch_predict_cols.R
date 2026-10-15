# 补丁脚本 v2predict 谱图相似性列全量重算。
# Patch script v2 recompute ALL predict similarity columns.
#
# v2 变更（设计约定"576 两文件夹所有行的理论谱图相似性都重新计算"）：
#   不再只处理 predict_missing=1 的行，而是对 report_A/classification_A/
#   conflict 的全部行重算 candidates 列中各候选的 pDeep2 余弦；对
#   master 与两分法报告的 cand1..cand5_predict 按各自候选序列与修饰直接
#   重算（包括不在 stage A candidates 列中的 a 槽候选）。键回退链：
#   seq|mod|z → seq||z → SAP 替换/Xle 展开变体 → I/L 互换变体，取最高余弦。
#   投票与判定列一律不动（predict 不参与投票）。
# v2 changes (user ruling "recompute the theoretical-spectrum similarity for
#   ALL rows of the two 576 folders"): instead of only the predict_missing=1
#   rows, all rows of report_A / classification_A / conflict get their
#   candidates-column cosines recomputed, and master / the two-part reports
#   get cand1..cand5_predict recomputed directly from each candidate's own
#   sequence and modification (including a-slot candidates absent from the
#   stage-A candidates column). Key fallback chain: seq|mod|z -> seq||z ->
#   SAP-substituted / Xle-expanded variants -> I/L-swapped variants, taking
#   the highest cosine. Vote and judgment columns are never touched (predict
#   does not vote).
#
# 用法: Rscript patch_predict_cols.R <result_dir> <all_mgf_dir> [predict_path]
# Usage: Rscript patch_predict_cols.R <result_dir> <all_mgf_dir> [predict_path]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript patch_predict_cols.R <result_dir> <all_mgf_dir> [predict_path]")
rd <- args[1]
mgf_dir <- args[2]
# 代码目录取自脚本自身位置（不硬编码中文路径：Windows 控制台源文件编码
# 不一定是 UTF-8，硬编码中文路径会被读成乱码）
# code dir derived from the script's own location (no hard-coded Chinese path:
# the Windows console source-file encoding is not guaranteed UTF-8, a
# hard-coded Chinese path can be read as mojibake)
argv0 <- commandArgs(FALSE)
farg <- sub("^--file=", "", argv0[grepl("^--file=", argv0)][1])
code_dir <- normalizePath(dirname(farg), winslash = "/", mustWork = TRUE)
setwd_try <- function(d) tryCatch({ setwd(d); TRUE }, error = function(e) FALSE)
ok <- setwd_try(code_dir)
if (!ok) {
  d2 <- dirname(farg)
  if (nzchar(d2) && d2 != ".") ok <- setwd_try(d2)
}
if (!ok && !file.exists(file.path(getwd(), "compat.R"))) {
  stop("cannot locate the code directory; run this script from the code folder or by its absolute path")
}
suppressMessages(source("compat.R"))
suppressMessages(source("acs_core.R"))
suppressMessages(source("paired_params.R"))
suppressMessages(source("ReadMGF_script.R"))
suppressMessages(source("chimeric_check.R"))
aa_path <- if (file.exists("AA_mass_dataframe.txt")) {
  "AA_mass_dataframe.txt"
} else {
  file.path(code_dir, "AA_mass_dataframe.txt")
}
aa <- read.delim(aa_path, sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# predict 文件：第三参数显式指定（576 流 predict_file 在项目外），缺省按
#  新布局 merge/predictors/pdeep2/predict.txt → 旧布局 merge/predict.txt
# 两级回退（旧项目仍可用）
predict_path <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  cands <- c(file.path(rd, "merge", "predictors", "pdeep2", "predict.txt"),
             file.path(rd, "merge", "predict.txt"))
  hit <- cands[file.exists(cands)]
  if (length(hit) > 0) hit[1] else cands[2]
}
pred <- read_predict_txt(predict_path)
pred_idx <- pred$idx
pred_int <- pred$intensities
cat("predict entries:", length(pred_int), "\n")

# ---- 谱图一次性加载（本项目涉及的全部 raw）------------------------------
# ---- load every spectrum once (all raws touched by this project) ---------
fns <- character(0)
mt_first <- list()   #：首遍读入的三大报告暂存，第二遍重算复用
ra_path <- file.path(rd, "1.unlabeled_result", "report_A.tsv")
if (file.exists(ra_path)) {
  ra0 <- read.delim(ra_path, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "",
                    na.strings = character(0))
  fns <- c(fns, as.character(ra0$File_Name))
}
for (f in c("master_table.tsv", "1.annotation_correction_report.tsv",
            "2.novel_peptide_competition_report.tsv")) {
  p <- file.path(rd, "2.labeled_result", f)
  if (file.exists(p)) {
    mt0 <- read.delim(p, sep = "\t", check.names = FALSE,
                      stringsAsFactors = FALSE, quote = "", comment.char = "",
                      na.strings = character(0))
    fns <- c(fns, as.character(mt0$File_Name))
    mt_first[[f]] <- mt0   #：第二遍全量重算直接复用，不再二次读
  }
}
fns <- unique(fns[!is.na(fns) & nzchar(fns)])
raws <- unique(unname(vapply(fns, function(t) unname(title_raw_scan(t)["raw"]),
                             character(1))))
raws <- raws[!is.na(raws)]
lu <- load_unlabeled_spectra(mgf_dir, only = raws)
cat("spectra loaded for raws:", length(raws), "\n")

# ---- 单候选余弦：键回退链 + 最高余弦 -------------------------------------
# ---- per-candidate cosine: key fallback chain, take the maximum ----------
compute_one <- function(fn, sq, md) {
  sq <- as.character(sq)
  md <- as.character(md)
  if (is.na(sq) || sq == "") return(NA_real_)
  if (is.na(md)) md <- ""
  z <- sub(".*\\.(\\d)\\.\\d\\.dta$", "\\1", fn)
  md_ <- sub(";+$", "", md)
  raw <- unname(title_raw_scan(fn)["raw"])
  sp <- NULL
  objs <- lu$objs[[raw]]
  if (!is.null(objs)) {
    li <- match(fn, lu$titles[[raw]])
    if (!is.na(li)) sp <- objs[[li]]@spectrum
  }
  if (is.null(sp)) return(NA_real_)
  r <- tryCatch(suppressWarnings(invisible(capture.output(
    sp2 <- spectrum_params(sq, md, "unlabeled", aa, aa_tab, sp, 20, "ppm")))),
    error = function(e) NULL)
  sp2 <- if (is.null(r)) NULL else sp2
  if (is.null(sp2)) return(NA_real_)
  v <- sp2$b_y_ions_intensity
  if (length(v) == 0) return(NA_real_)
  keys <- c(paste0(sq, "|", md_, "|", z), paste0(sq, "||", z))
  vars <- tryCatch(apply_sap_variants(sq, md), error = function(e) NULL)
  if (!is.null(vars) && length(vars) >= 1) {
    for (vr in vars) {
      keys <- c(keys, paste0(vr$seq, "|", sub(";+$", "", vr$mod), "|", z),
                paste0(vr$seq, "||", z))
    }
  }
  il1 <- gsub("I", "L", sq)
  il2 <- gsub("L", "I", sq)
  if (il1 != sq) keys <- c(keys, paste0(il1, "|", md_, "|", z), paste0(il1, "||", z))
  if (il2 != sq) keys <- c(keys, paste0(il2, "|", md_, "|", z), paste0(il2, "||", z))
  keys <- unique(keys)
  best <- NA_real_
  for (k in keys) {
    pi <- unname(pred_idx[k])
    if (is.na(pi)) next
    pv <- pred_int[[pi]]
    nmin <- min(length(pv), length(v))
    if (nmin > 0) {
      c <- cosine(v[seq_len(nmin)], pv[seq_len(nmin)])
      if (is.na(best) || c > best) best <- c
    }
  }
  best
}

fmt <- function(x) ifelse(is.na(x), "", sprintf("%.4f", x))

# ---- report_A / classification_A / conflict：全部行重算 candidates 列 ----
# ---- report_A / classification_A / conflict: recompute every row --------
if (exists("ra0")) {
  pm1 <- seq_len(nrow(ra0))
  fixed <- 0L
  pc_new <- rep("", nrow(ra0))
  for (i in pm1) {
    fn <- as.character(ra0$File_Name[i])
    cands <- strsplit(as.character(ra0$candidates[i]), "|", fixed = TRUE)[[1]]
    vals <- character(length(cands))
    for (k in seq_along(cands)) {
      c <- cands[k]
      if (c == "") next
      sq <- sub("\\(.*$", "", c)
      md <- if (grepl("(", c, fixed = TRUE)) {
        sc <- sub("^.*\\(", "(", c)
        m <- regmatches(sc, regexec("^\\((.*)\\)$", sc))[[1]]
        if (length(m) > 1) m[2] else ""
      } else ""
      v <- compute_one(fn, sq, md)
      if (!is.na(v)) {
        vals[k] <- sprintf("%.4f", v)
        fixed <- fixed + 1L
      }
    }
    pc_new[i] <- paste(vals, collapse = "|")
  }
  ra0$candidate_predict_cosine <- pc_new
  ra0$predict_missing <- vapply(strsplit(pc_new, "|", fixed = TRUE),
                                function(x) any(x == ""), integer(1))
  bk <- paste0(ra_path, ".pre_predict_v2.bak")
  if (!file.exists(bk)) file.copy(ra_path, bk)
  write.table(ra0, ra_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("report_A 重算候选余弦:", fixed, "个（行数", nrow(ra0), "）\n")

  ca_path <- file.path(rd, "1.unlabeled_result", "classification_A.tsv")
  if (file.exists(ca_path)) {
    ca <- read.delim(ca_path, sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE, quote = "", comment.char = "",
                     na.strings = character(0))
    ii <- match(as.character(ra0$File_Name), as.character(ca$File_Name))
    ok <- !is.na(ii)
    if (any(ok)) {
      ca$candidate_predict_cosine[ii[ok]] <- pc_new[ok]
      ca$predict_missing[ii[ok]] <- ra0$predict_missing[ok]
      bk2 <- paste0(ca_path, ".pre_predict_v2.bak")
      if (!file.exists(bk2)) file.copy(ca_path, bk2)
      write.table(ca, ca_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
      cat("classification_A 已同步:", sum(ok), "行\n")
    }
  }

  cr_path <- file.path(rd, "2.labeled_result", "conflict_report_n15.tsv")
  if (file.exists(cr_path)) {
    cr <- read.delim(cr_path, sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE, quote = "", comment.char = "",
                     na.strings = character(0))
    ii <- match(as.character(ra0$File_Name), as.character(cr$File_Name))
    ok <- !is.na(ii)
    if (any(ok)) {
      cr$candidate_predict_cosine[ii[ok]] <- pc_new[ok]
      if ("predict_cosines" %in% names(cr))
        cr$predict_cosines[ii[ok]] <- pc_new[ok]
      cr$predict_missing[ii[ok]] <- ra0$predict_missing[ok]
      bk3 <- paste0(cr_path, ".pre_predict_v2.bak")
      if (!file.exists(bk3)) file.copy(cr_path, bk3)
      write.table(cr, cr_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
      cat("conflict_report_n15 已同步:", sum(ok), "行\n")
    }
  }
}

# ---- master + 两分法报告：cand1..cand5_predict 直接全量重算 -------------
# ---- master + two-part reports: cand1..cand5_predict recomputed ----------
for (f in c("master_table.tsv", "1.annotation_correction_report.tsv",
            "2.novel_peptide_competition_report.tsv")) {
  p <- file.path(rd, "2.labeled_result", f)
  if (!file.exists(p)) next
  mt <- mt_first[[f]]   #：复用首遍读入（文件在本段之前未被改动）
  if (is.null(mt)) next
  filled <- 0L
  for (k in 1:5) {
    scol <- paste0("cand", k, "_sequence")
    mcol <- paste0("cand", k, "_modification")
    pcol <- paste0("cand", k, "_predict")
    if (!(scol %in% names(mt)) || !(pcol %in% names(mt))) next
    for (i in seq_len(nrow(mt))) {
      sq <- as.character(mt[[scol]][i])
      if (is.na(sq) || sq == "") next
      md <- if (mcol %in% names(mt)) as.character(mt[[mcol]][i]) else ""
      if (is.na(md)) md <- ""
      v <- compute_one(as.character(mt$File_Name[i]), sq, md)
      if (!is.na(v)) {
        mt[[pcol]][i] <- sprintf("%.4f", v)
        filled <- filled + 1L
      } else {
        mt[[pcol]][i] <- ""
      }
    }
  }
  bk <- paste0(p, ".pre_predict_v2.bak")
  if (!file.exists(bk)) file.copy(p, bk)
  write.table(mt, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("已更新:", f, "（填充候选余弦", filled, "个）\n")
}
cat("PATCH V2 DONE:", rd, "\n")
