# mzml_cache_test.R — v1.0 内置 MGF→mzML 转换共享缓存回归（2026-09-28）
# 覆盖：首转换落缓存、跨引擎目录第二槽缓存命中（复制且内容一致）、
# 产物新于源 MGF 的 0 成本跳过（#18 口径）、源字节数变化换缓存键重转。
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

# 项目级布局：proj/{msfragger_auto,sage_auto}/<stem>_renum.mgf（模拟两引擎
# 各自的重编号产物——同名同内容、目录不同、mtime 不同）
proj <- tempfile("mzc_"); dir.create(file.path(proj, "msfragger_auto"), recursive = TRUE)
dir.create(file.path(proj, "sage_auto"), recursive = TRUE)
mgf_lines <- c(
  "BEGIN IONS", "TITLE=raw9.1.1.2.0.dta", "PEPMASS=456.789 12345", "CHARGE=2+",
  "RTINSECONDS=12.34", "100.1 1000", "200.2 2000", "300.3 3000", "END IONS",
  "BEGIN IONS", "TITLE=raw9.2.2.3.0.dta", "PEPMASS=789.123 678", "CHARGE=3+",
  "RTINSECONDS=56.78", "400.4 4000", "500.5 5000", "END IONS")
rmgf1 <- file.path(proj, "msfragger_auto", "raw9_renum.mgf")
rmgf2 <- file.path(proj, "sage_auto", "raw9_renum.mgf")
writeLines(mgf_lines, rmgf1, useBytes = TRUE)
writeLines(mgf_lines, rmgf2, useBytes = TRUE)
out1 <- file.path(proj, "msfragger_auto", "raw9_renum_pw.mzML")
out2 <- file.path(proj, "sage_auto", "raw9_renum_pw.mzML")

# 1) 首槽：真转换 + 落缓存
msg1 <- capture.output(mgf_to_mzml_cached(rmgf1, out1))
cache_dir <- file.path(proj, "mzml_cache")
cf <- list.files(cache_dir, pattern = "^raw9_renum_.*\\.mzML$", full.names = TRUE)
ck("首转换产出 mzML（含 spectrum）", file.exists(out1) &&
   grepl("<spectrum", paste(readLines(out1, warn = FALSE), collapse = " "), fixed = TRUE))
ck("缓存文件已落盘（键含字节数）", length(cf) == 1L &&
   grepl(paste0(format(file.size(rmgf1), scientific = FALSE, trim = TRUE), ".mzML$"), cf))

# 2) 第二槽：缓存命中、复制且逐字节一致
msg2 <- capture.output(mgf_to_mzml_cached(rmgf2, out2))
ck("第二槽命中共享缓存（提示语）", any(grepl("命中共享缓存", msg2)))
ck("两槽产物逐字节一致", identical(readLines(out1, warn = FALSE),
                                    readLines(out2, warn = FALSE)))

# 3) 产物新于源 MGF → 0 成本跳过（#18 口径）
msg3 <- capture.output(mgf_to_mzml_cached(rmgf1, out1))
ck("产物已新于源 → 跳过", any(grepl("跳过", msg3)))

# 4) 源字节数变化 → 新键重转（缓存条目 +1）
cat("EXTRA 0 0\n", file = rmgf2, append = TRUE)
msg4 <- capture.output(mgf_to_mzml_cached(rmgf2, out2))
cf2 <- list.files(cache_dir, pattern = "^raw9_renum_.*\\.mzML$", full.names = TRUE)
ck("源大小变化换键重转（缓存 2 条）", length(cf2) == 2L &&
   !any(grepl("命中共享缓存", msg4)))

unlink(proj, recursive = TRUE)
cat("==== mzml_cache_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
