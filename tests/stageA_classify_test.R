# stageA_classify_test.R — 阶段 A 分类向量化对拍（v1.0 优化）
# 子进程运行真实 stageA_classify.R（UNLABELED_DIR 置空走报告列口径），
# 与测试内旧逐行实现逐列比对 7 个结论列（conclusion_status / vote_says_wrong /
# error_reason / modification_diff / chimeric / identical_ions / quality_low）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
suppressMessages(source(file.path(cd, "boot.R")))
setwd(cd)   # AA_mass_dataframe.txt / 子进程 cwd 都在 code/

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 旧逐行实现（v1.0 版本原样照抄，仅去掉文件读写）
old_classify <- function(ext, aa, aa_tab, lu = NULL) {
  parse_mods <- function(s) {
    s <- as.character(s)
    if (!nzchar(s) || is.na(s)) return(character(0))
    parts <- strsplit(s, ";", fixed = TRUE)[[1]]
    parts <- trimws(parts); parts[nzchar(parts)]
  }
  mod_diff_detail <- function(a, b) {
    pa <- parse_mods(a); pb <- parse_mods(b)
    da <- setdiff(pa, pb); db <- setdiff(pb, pa)
    paste(c(if (length(da) > 0) paste0("original extra: ", paste(da, collapse = ";")),
            if (length(db) > 0) paste0("winner extra: ", paste(db, collapse = ";"))),
          collapse = " | ")
  }
  identical_ion_sets <- function(s1, m1, s2, m2) {
    if (s1 == "" || s2 == "") return(FALSE)
    t1 <- try(suppressWarnings(invisible(capture.output(
      r1 <- by_ions_position(s1, m1, "Unlabel", aa, aa_tab)))), silent = TRUE)
    t2 <- try(suppressWarnings(invisible(capture.output(
      r2 <- by_ions_position(s2, m2, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(t1, "try-error") || inherits(t2, "try-error")) return(FALSE)
    if (length(r1) == 1 || length(r2) == 1) return(FALSE)
    a <- as.matrix(r1$ions_table[, -1, drop = FALSE])
    b <- as.matrix(r2$ions_table[, -1, drop = FALSE])
    identical(dim(a), dim(b)) &&
      isTRUE(all.equal(a, b, tolerance = 1e-6, check.attributes = FALSE))
  }
  same_composition <- function(s1, s2) {
    comp <- function(s) sort(strsplit(gsub("I", "L", as.character(s)), "", fixed = TRUE)[[1]])
    identical(comp(s1), comp(s2))
  }
  n <- nrow(ext)
  res <- vector("list", n)
  slot_seq_cols <- grep("^[a-z]\\.Sequence$", colnames(ext), value = TRUE)
  slot_mod_cols <- grep("^[a-z]\\.Modification$", colnames(ext), value = TRUE)
  for (i in seq_len(n)) {
    is_conf <- (ext$Peptide_state[i] != "Same") | (ext$Modification_state[i] != "Same")
    sqs <- vapply(slot_seq_cols, function(cn) as.character(ext[i, cn]), character(1))
    mds <- vapply(slot_mod_cols, function(cn) as.character(ext[i, cn]), character(1))
    cl <- if ("candidates" %in% names(ext)) {
      strsplit(as.character(ext$candidates[i]), "|", fixed = TRUE)[[1]]
    } else character(0)
    cand_sqs <- sub("\\(.*$", "", cl)
    has_mod <- grepl("\\(", cl)
    cand_mds <- ifelse(has_mod, sub("^.*\\((.*)\\)$", "\\1", cl), "")
    cand_mds <- rep_len(cand_mds, length(cand_sqs))
    x1s <- if (length(sqs) >= 1 && any(slot_seq_cols == "a.Sequence")) sqs[which(slot_seq_cols == "a.Sequence")] else
      (if (length(cand_sqs) >= 1) cand_sqs[1] else "")
    x1m <- if (length(mds) >= 1 && any(slot_mod_cols == "a.Modification")) mds[which(slot_mod_cols == "a.Modification")] else
      (if (length(cand_mds) >= 1) cand_mds[1] else "")
    status <- "no competition (strategies agree)"
    vote_wrong <- 0L; err_reason <- ""; mod_diff <- ""
    chimeric <- 0L; ident_ions <- 0L; qual_low <- 0L
    if (is_conf && length(cand_sqs) >= 1) {
      win <- suppressWarnings(as.integer(ext$best_group[i]))
      if (is.na(win) || win < 1 || win > length(cand_sqs)) win <- 1L
      ws <- cand_sqs[win]; wm <- cand_mds[win]
      vote_wrong <- as.integer(!(ws == x1s && wm == x1m))
      if (!is.null(lu)) {
        chimeric <- chimeric_check(cand_sqs, cand_mds, lu,
                                   as.character(ext$File_Name[i]), aa, aa_tab)
      } else {
        sr <- if ("candidate_specific_ratio" %in% names(ext)) {
          suppressWarnings(as.numeric(strsplit(
            as.character(ext$candidate_specific_ratio[i]), "|", fixed = TRUE)[[1]]))
        } else numeric(0)
        sr <- sort(sr[!is.na(sr)], decreasing = TRUE)
        if (length(sr) >= 2 && sr[1] > 0.1 && sr[2] > 0.1) chimeric <- 1L
      }
      cov_all <- suppressWarnings(as.numeric(strsplit(
        as.character(ext$candidate_coverage[i]), "|", fixed = TRUE)[[1]]))
      cov <- if (length(cov_all) >= win && !is.na(cov_all[win])) cov_all[win] else NA_real_
      if (!is.na(cov) && cov < 20) qual_low <- 1L
      byv <- suppressWarnings(as.numeric(strsplit(
        as.character(ext$candidate_PSM_score_by[i]), "|", fixed = TRUE)[[1]]))
      others <- setdiff(which(!is.na(byv) & seq_along(byv) <= length(cand_sqs)), win)
      if (length(others) > 0) {
        runner <- others[which.max(byv[others])]
        rs <- cand_sqs[runner]; rm_ <- cand_mds[runner]
        if (ws != "" && rs != "" && ws != rs) {
          ident_ions <- as.integer(identical_ion_sets(ws, wm, rs, rm_))
        }
      }
      if (chimeric == 1L) {
        status <- "chimeric spectrum"
      } else if (ident_ions == 1L) {
        status <- "ambiguous: identical fragment ions"
      } else if (qual_low == 1L) {
        status <- "ambiguous: poor spectrum quality"
      } else if (vote_wrong == 1L) {
        status <- "original result wrong"
      } else {
        status <- "competition but original correct"
      }
      if (vote_wrong == 1L) {
        if (nzchar(as.character(ext$Error_reason[i]))) {
          err_reason <- as.character(ext$Error_reason[i])
        }
      if (ws != "" && identical(unname(x1s), ws)) {   # P2-1 修复后口径：只比值
        if (x1m != wm) {
          if (!nzchar(err_reason)) err_reason <- "modification_change"
          mod_diff <- mod_diff_detail(x1m, wm)
        }
      } else if (ws != "" && x1s != "") {
          if (!nzchar(err_reason)) {
            err_reason <- if (same_composition(x1s, ws)) {
              "Same-weight AA replacement or permutation"
            } else {
              "sequence_change"
            }
          }
        }
      }
    }
    res[[i]] <- data.frame(
      ext[i, , drop = FALSE],
      conclusion_status = status, vote_says_wrong = vote_wrong,
      error_reason = err_reason, modification_diff = mod_diff,
      chimeric = chimeric, identical_ions = ident_ions, quality_low = qual_low,
      stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  }
  do.call(rbind, res)
}

aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

td <- tempfile("stga_"); dir.create(td)
ext_f <- file.path(td, "report_A.tsv")
out_f <- file.path(td, "classification_A.tsv")

# 合成报告表：覆盖 无竞争/修饰改变/序列改变/同重置换/混合谱/质量差/子离子
# 一致/best_group 越界/空候选 等路径
ext <- data.frame(
  File_Name = c("raw1.100.100.2.dta", "raw1.200.200.2.dta", "raw1.300.300.2.dta",
                "raw1.400.400.2.dta", "raw1.500.500.2.dta", "raw1.600.600.2.dta",
                "raw1.700.700.2.dta", "raw1.800.800.2.dta", "raw1.900.900.2.dta"),
  Peptide_state = c("Same", "Diff", "Diff", "Diff", "Diff", "Same", "Diff", "Diff", "Diff"),
  Modification_state = c("Same", "Same", "Diff", "Same", "Same", "Same", "Diff", "Same", "Same"),
  a.Sequence = c("PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "PEPTIDEK",
                 "PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "PEPTIDEK"),
  a.Modification = c("", "4,Oxidation[M];", "", "", "", "", "", "", ""),
  b.Sequence = c("", "PEPTIDER", "PEPTIDEK", "PEPTIDEK", "PEPTIDEK", "",
                 "PEPTIDEK", "PEPTIDEL", "PEPTIDEK"),
  b.Modification = c("", "", "4,Oxidation[M];", "", "", "", "", "", ""),
  candidates = c("PEPTIDEK()",
                 "PEPTIDEK(4,Oxidation[M])|PEPTIDER()",
                 "PEPTIDEK()|PEPTIDEK(4,Oxidation[M])",
                 "PEPTIDEK()|PEPTIDER()",
                 "PEPTIDEK()|PEPTIDER()",
                 "PEPTIDEK()",
                 "PEPTIDEK()|PEPTIDEK(4,Oxidation[M])",
                 "PEPTIDEK()|PEPTIDEL()",
                 ""),
  best_group = c("1", "2", "2", "2", "2", "1", "2", "99", "1"),
  candidate_specific_ratio = c("1", "0.2|0.05", "0.5|0.4", "0.9|0.01", "0.9|0.01",
                               "1", "0.2|0.05", "0.2|0.05", ""),
  candidate_coverage = c("80", "50|30", "60|40", "50|10", "60|30",
                         "80", "60|40", "50|30", ""),
  candidate_PSM_score_by = c("0.9", "0.5|0.8", "0.6|0.7", "0.5|0.8", "0.5|0.8",
                             "0.9", "0.6|0.7", "0.5|0.8", ""),
  Error_reason = c("", "", "", "", "", "", "", "", ""),
  stringsAsFactors = FALSE, check.names = FALSE)
write.table(ext, ext_f, sep = "\t", row.names = FALSE, quote = FALSE)

# 子进程运行真实（新）实现
rscript <- file.path(R.home("bin"), "Rscript.exe")
Sys.setenv(EXT_REPORT = ext_f, OUT = out_f, UNLABELED_DIR = "")
logf <- file.path(td, "run.log")
st <- system2(rscript, c(shQuote(file.path(cd, "stageA_classify.R"))),
              stdout = logf, stderr = logf)
ck("子进程运行成功", st == 0 && file.exists(out_f))
new_out <- read.delim(out_f, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, na.strings = character(0))
old_out <- old_classify(ext, aa, aa_tab)

concl_cols <- c("conclusion_status", "vote_says_wrong", "error_reason",
                "modification_diff", "chimeric", "identical_ions", "quality_low")
# 读回伪影校正：本夹具下 modification_diff 可能全空，read.delim 读成 logical
# NA——与生产脚本入口相同的 logical 恢复后再比对
new_out[] <- lapply(new_out, function(x) if (is.logical(x)) rep("", length(x)) else x)
for (cn in concl_cols) {
  ck(paste0("结论列逐行一致：", cn),
     identical(as.character(new_out[[cn]]), as.character(old_out[[cn]])))
}
ck("行数一致", nrow(new_out) == nrow(old_out))
ck("状态分布一致",
   identical(table(new_out$conclusion_status), table(old_out$conclusion_status)))
# P2-1 修复后口径（2026-09-26，用户拍板）：同序列仅修饰不同的行（夹具 3/7
# 行）必须标 modification_change 且 modification_diff 出值，不再误标
# Same-weight AA replacement（该分支是此前 identical 比属性的死代码路径）
ck("同序列仅修饰变化标 modification_change（P2-1）",
   as.character(new_out$error_reason[3]) == "modification_change" &&
   as.character(new_out$error_reason[7]) == "modification_change")
ck("modification_diff 列出值（P2-1）",
   grepl("winner extra: 4,Oxidation[M]",
         as.character(new_out$modification_diff[3]), fixed = TRUE) &&
   grepl("winner extra: 4,Oxidation[M]",
         as.character(new_out$modification_diff[7]), fixed = TRUE))

unlink(td, recursive = TRUE)
cat(sprintf("=== stageA_classify_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
