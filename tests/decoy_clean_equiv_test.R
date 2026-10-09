# decoy_clean_equiv_test.R — clean_decoy_spectra 向量化对拍（v1.0 优化）
# 用合成 .spectra 夹具断言：新向量化实现与旧逐行循环参考实现逐字节同义
# （含无 Proteins 列 + 注释分支、fa 正库兜底、全诱饵/混合/空条目/Positions
# 不等长等边界）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 旧逐行循环参考实现（两种分支原样照抄）
old_clean_no_proteins <- function(dat, annot, fa) {
  n_in <- nrow(dat)
  m <- match(paste(dat$File_Name, dat$Sequence),
             paste(annot$File_Name, annot$Sequence))
  ap <- rep("", n_in); ap[!is.na(m)] <- annot$Proteins[m[!is.na(m)]]
  parts <- strsplit(ap, "/", fixed = TRUE)
  keep <- rep(TRUE, n_in); n_stripped <- 0L
  for (i in seq_len(n_in)) {
    ps <- trimws(parts[[i]]); ps <- ps[nzchar(ps)]
    if (length(ps) == 0) {
      if (!is.null(fa) && length(fa) > 0 &&
          !any(grepl(dat$Sequence[i], fa, fixed = TRUE))) keep[i] <- FALSE
      next
    }
    isrev <- grepl("^REV_", ps)
    if (all(isrev)) keep[i] <- FALSE
    else if (any(isrev)) n_stripped <- n_stripped + 1L
  }
  list(keep = keep, n_stripped = n_stripped)
}
old_clean_with_proteins <- function(dat) {
  n_in <- nrow(dat)
  pr <- dat$Proteins
  parts <- strsplit(pr, "/", fixed = TRUE)
  keep <- rep(TRUE, n_in); n_stripped <- 0L
  new_pr <- pr
  new_po <- if ("Positions" %in% colnames(dat)) dat$Positions else NULL
  po_parts <- if (!is.null(new_po)) strsplit(new_po, "/", fixed = TRUE) else NULL
  for (i in seq_len(n_in)) {
    ps <- trimws(parts[[i]]); ps <- ps[nzchar(ps)]
    if (length(ps) == 0) next
    isrev <- grepl("^REV_", ps)
    if (all(isrev)) { keep[i] <- FALSE; next }
    if (any(isrev)) {
      n_stripped <- n_stripped + 1L
      ps2 <- ps[!isrev]
      new_pr[i] <- paste0(paste(ps2, collapse = "/"), "/")
      if (!is.null(po_parts)) {
        posparts <- trimws(po_parts[[i]]); posparts <- posparts[nzchar(posparts)]
        if (length(posparts) == length(ps)) {
          new_po[i] <- paste0(paste(posparts[!isrev], collapse = "/"), "/")
        }
      }
    }
  }
  list(keep = keep, n_stripped = n_stripped, new_pr = new_pr, new_po = new_po)
}

td <- tempfile("dcd_"); dir.create(td)
cln <- file.path(td, "cleaned")
dir.create(cln, recursive = TRUE, showWarnings = FALSE)

# 夹具 1：Proteins 列分支（混合/全诱饵/空条目/Positions 不等长）
dat1 <- data.frame(
  File_Name = c("raw1.100.100.2.dta", "raw1.200.200.2.dta", "raw1.300.300.2.dta",
                "raw1.400.400.2.dta", "raw1.500.500.2.dta"),
  Sequence = c("PEPTIDEK", "PEPTIDER", "PEPTIDES", "PEPTIDET", "PEPTIDEV"),
  Proteins = c("sp|A|P1/REV_sp|B|P2/", "REV_sp|C|P3/REV_sp|D|P4/",
               "sp|E|P5/", "", "sp|F|P6//"),
  Positions = c("3/5/", "1/2/", "4/", "9/", "7/8/"),
  stringsAsFactors = FALSE, check.names = FALSE)
write.table(dat1, file.path(td, "pfind-Filtered.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
new1 <- clean_decoy_spectra(file.path(td, "pfind-Filtered.spectra"), cln, "t1")
old1 <- old_clean_with_proteins(dat1)
new1_rows <- read.delim(new1$file, sep = "\t", colClasses = "character",
                        na.strings = character(0))
ck("Proteins 分支 keep 掩码一致", identical(new1_rows$Sequence, dat1$Sequence[old1$keep]))
ck("Proteins 分支 n_stripped 一致", new1$n_stripped == old1$n_stripped)
ck("Proteins 分支重建 Proteins 一致", identical(new1_rows$Proteins, old1$new_pr[old1$keep]))
ck("Proteins 分支重建 Positions 一致", identical(new1_rows$Positions, old1$new_po[old1$keep]))

# 夹具 2：无 Proteins 列 + 注释分支（含 fa 正库兜底）
dat2 <- data.frame(
  File_Name = c("raw2.100.100.2.dta", "raw2.200.200.2.dta", "raw2.300.300.2.dta",
                "raw2.400.400.2.dta"),
  Sequence = c("AAAAAA", "BBBBBB", "CCCCCC", "DDDDDD"),
  stringsAsFactors = FALSE, check.names = FALSE)
write.table(dat2, file.path(td, "nop.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
annot2 <- data.frame(
  File_Name = c("raw2.100.100.2.dta", "raw2.200.200.2.dta", "raw2.300.300.2.dta"),
  Sequence = c("AAAAAA", "BBBBBB", "CCCCCC"),
  Proteins = c("sp|A|P1/REV_sp|B|P2/", "REV_sp|C|P3/", ""),
  stringsAsFactors = FALSE)
fa2 <- c(">sp|A|P1", "AAAAAACCCC", ">sp|E|P9", "DDDDDD")
new2 <- clean_decoy_spectra(file.path(td, "nop.spectra"), cln, "t2",
                            annot = annot2, fa = fa2)
old2 <- old_clean_no_proteins(dat2, annot2, fa2)
new2_rows <- read.delim(new2$file, sep = "\t", colClasses = "character",
                        na.strings = character(0))
ck("无 Proteins 分支 keep 掩码一致", identical(new2_rows$Sequence, dat2$Sequence[old2$keep]))
ck("无 Proteins 分支 n_stripped 一致", new2$n_stripped == old2$n_stripped)
ck("无 Proteins 分支 fa 兜底一致（CCCCCC 删、DDDDDD 留）",
   all(c("CCCCCC" %in% dat2$Sequence[!old2$keep], "DDDDDD" %in% new2_rows$Sequence)))

# 夹具 3：空表（n_in = 0）不报错
dat3 <- dat1[0, , drop = FALSE]
write.table(dat3, file.path(td, "empty.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
new3 <- clean_decoy_spectra(file.path(td, "empty.spectra"), cln, "t3")
ck("空表分支 n_out=0", new3$n_out == 0L)

unlink(td, recursive = TRUE)
cat(sprintf("=== decoy_clean_equiv_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
