# sanitize_mods_equiv_test.R — 合并路径修饰越界过滤向量化对拍（v1.0 优化）
# 新展平判定与旧 sanitize_pdeep2_mods 逐 token 等价（空 token/不可识别位置
# 保留、越界整数位点丢弃、分号结构原样保留）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
suppressMessages(source(file.path(cd, "compat.R")))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 新实现（与 ReadMGF 内联逻辑一致）
new_sanitize_all <- function(tmp) {
  nonempty <- which(nzchar(tmp$modification))
  if (length(nonempty) > 0) {
    mods_n <- as.character(tmp$modification)[nonempty]
    lens_n <- nchar(tmp$peptide)[nonempty]
    toks <- strsplit(mods_n, ";", fixed = TRUE)
    flat_t <- unlist(toks, use.names = FALSE)
    g <- rep(seq_along(nonempty), lengths(toks))
    pos <- suppressWarnings(as.integer(sub(",.*$", "", flat_t)))
    bad <- !is.na(pos) & (pos < 0L | pos > lens_n[g] + 1L)
    if (any(bad)) {
      bad_rows <- unique(g[bad])
      keep_tok <- !bad
      pos_of <- cumsum(lengths(toks)) - lengths(toks) + 1L
      rebuilt <- vapply(bad_rows, function(ri) {
        sel <- pos_of[ri] + seq_len(lengths(toks)[ri]) - 1L
        paste(flat_t[sel][keep_tok[sel]], collapse = ";")
      }, "")
      mods_n[bad_rows] <- rebuilt
      clean_rows <- setdiff(seq_along(nonempty), bad_rows)
      mods_n[clean_rows] <- sub(";$", "", mods_n[clean_rows])
    } else {
      mods_n <- sub(";$", "", mods_n)
    }
    tmp$modification[nonempty] <- mods_n
  }
  tmp$modification
}

set.seed(7)
n <- 4000
mods <- c("4,Oxidation[M];", "0,Acetyl[AnyN-term];", "64,Oxidation[M];",
          "-2,Phospho[S];3,Oxidation[M];", "", "9,Unknown[X];;;", "33,Deamidated[N];")
tmp <- data.frame(
  peptide = sprintf("PEP%03d", sample(5:30, n, replace = TRUE)),
  modification = sample(mods, n, replace = TRUE),
  stringsAsFactors = FALSE)
new_v <- new_sanitize_all(tmp)
old_v <- vapply(seq_len(n), function(r)
  sanitize_pdeep2_mods(tmp$modification[r], nchar(tmp$peptide[r])), "")
ck("向量化过滤与旧逐格过滤逐行一致", identical(as.character(new_v), as.character(old_v)))
ck("越界 64 位点被丢弃", any(grepl("64,Oxidation", tmp$modification)) &&
   !any(grepl("64,Oxidation", new_v)))
# 空表/全空修饰不报错
tmp0 <- data.frame(peptide = "PEPTIDEK", modification = "", stringsAsFactors = FALSE)
ck("空修饰表不报错", identical(new_sanitize_all(tmp0), tmp0$modification))

cat(sprintf("=== sanitize_mods_equiv_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
