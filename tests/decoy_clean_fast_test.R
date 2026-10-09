# decoy_clean_fast_test.R — v1.0 诱饵清洗行级过滤快速路径回归（2026-09-29）
# 覆盖：纯目标透传（字节保真）、纯诱饵整行删、混合行剥诱饵（Positions 条目数
# 同步/失配两分支）、空 Proteins 行、"REV_" 出现在他列的假阳性（判定后原样
# 保留）、行序保持、fast 与 _base 的计数与解析后表相等、无 Proteins 列回落
# _base、annot 场景不走快速路径。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "boot.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

td <- tempfile("dcf_"); dir.create(td)
dir.create(file.path(td, "f")); dir.create(file.path(td, "b"))
hdr <- "File_Name\tScan\tSequence\tModification\tProteins\tPositions\tOthers"
mk <- function(fn, sc, seq, mod, pr, po, ot) {
  paste(fn, sc, seq, mod, pr, po, ot, sep = "\t")
}
rows <- c(
  mk("r1.1.1", 1, "AAAAK", "", "sp|P1|A_HUMAN/", "1,K,I/", "10"),             # 纯目标
  mk("r1.2.2", 2, "BBBBK", "", "REV_sp|P2|B_HUMAN/", "2,K,I/", "20"),         # 纯诱饵 → 删
  mk("r1.3.3", 3, "CCCCK", "", "sp|P3|C_HUMAN/REV_sp|P4|D_HUMAN/", "3,K,I/3,K,I/", "30"),  # 混合 → 剥（同步）
  mk("r1.4.4", 4, "DDDDK", "", "sp|P5|E_HUMAN/REV_sp|P6|F_HUMAN/", "4,K,I/", "40"),        # 混合但 Positions 条目数失配 → 只剥 Proteins
  mk("r1.5.5", 5, "EEEE K", "", "", "", "50"),                                 # 空 Proteins → 保留
  mk("r1.6.6", 6, "FFFFK", "", "sp|P7|REVK_HUMAN/", "6,K,I/", "60"),          # 蛋白名含 REV 但非前缀 → 保留
  mk("r1.7.7", 7, "GGGGK", "1,Ox[REV_M];", "sp|P8|G_HUMAN/", "7,K,I/", "70"),  # 他列假阳性 → 原样保留
  paste("r1.8.8", "8", "HHHHK", "1,K,I", "sp|P9|I_HUMAN/", "8,K,I", sep = "\t"),  # 短行（5 字段）→ 补齐保留
  paste("r1.9.9", "9", "IIIIK", "1,K,I", "REV_sp|P10|J_HUMAN", sep = "\t")        # 短行含 REV（无 Proteins 字段）→ 保留补齐
)
src <- file.path(td, "in.spectra")
writeLines(c(hdr, rows), src, useBytes = TRUE)

rf <- clean_decoy_spectra(src, file.path(td, "f"), "t")
rb <- clean_decoy_spectra_base(src, file.path(td, "b"), "t")
ck("计数一致（2 删 2 剥 → 7 出；短行第 5 字段落 Proteins 列，全诱饵即删）",
   identical(rf$n_in, rb$n_in) &&
   rf$n_out == 7L && identical(rf$n_out, rb$n_out) &&
   identical(rf$n_dropped, rb$n_dropped) && rf$n_dropped == 2L &&
   identical(rf$n_stripped, rb$n_stripped) && rf$n_stripped == 2L)
ck("dropped_spectra 一致", setequal(rf$dropped_spectra, rb$dropped_spectra) &&
   setequal(rf$dropped_spectra, c("r1.2.2", "r1.9.9")))
lf <- readLines(rf$file, warn = FALSE); lb <- readLines(rb$file, warn = FALSE)
ck("fast 表头保留原始", identical(lf[1L], hdr) && length(lf) == length(lb))
ck("纯目标行字节透传", identical(lf[2], rows[1]))
ck("行序保持：剥后 r1.3.3（Positions 同步剔除）",
   grepl("^r1\\.3\\.3\t3\tCCCCK\t\tsp\\|P3\\|C_HUMAN/\t3,K,I/\t30$", lf[3]))
ck("r1.4.4 仅剥 Proteins（Positions 失配保留）",
   grepl("^r1\\.4\\.4\t4\tDDDDK\t\tsp\\|P5\\|E_HUMAN/\t4,K,I/\t40$", lf[4]))
ck("空 Proteins 行字节透传", identical(lf[5], rows[5]))
ck("蛋白名含 REV 非前缀：原样保留", identical(lf[6], rows[6]))
ck("他列假阳性：原样保留", identical(lf[7], rows[7]))
ck("短行补齐到表头宽（r1.8.8，6 字段补 1 空）", grepl("\t8,K,I\t$", lf[8]) &&
   lengths(regmatches(lf[8], gregexpr("\t", lf[8], fixed = TRUE))) == 6L)
ck("短行含 REV 且第 5 字段为全诱饵：整行删（r1.9.9 不在输出）",
   !any(grepl("^r1\\.9\\.9\t", lf)))
df_f <- read.delim(rf$file, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = "character", na.strings = character(0))
df_b <- read.delim(rb$file, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = "character", na.strings = character(0))
ck("解析后表完全相等（fast vs base）", identical(as.matrix(df_f), as.matrix(df_b)))

# 无 Proteins 列：回落 _base（原样复制）
src2 <- file.path(td, "np.spectra")
writeLines(c("File_Name\tSequence", "r2.1\tAAAK"), src2, useBytes = TRUE)
r2 <- clean_decoy_spectra(src2, file.path(td, "f"), "np")
r2b <- clean_decoy_spectra_base(src2, file.path(td, "f"), "np2")
ck("无 Proteins 列回落 _base 原样复制", r2$n_out == 1L &&
   identical(readLines(r2$file, warn = FALSE), readLines(r2b$file, warn = FALSE)))
# annot 场景不走快速路径（交给 _base 的回查逻辑）
ck("annot 传入时走 _base", !is.null(clean_decoy_spectra(src, file.path(td, "f"), "t3",
   annot = data.frame(File_Name = "x", Sequence = "y", Proteins = "z"))$file))

cat("==== decoy_clean_fast_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
