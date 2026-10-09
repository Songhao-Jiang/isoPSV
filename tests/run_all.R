# tests/run_all.R — 统一测试入口（P3-3，2026-09-25；F-1 locale 卫生，2026-09-26）
# 遍历 tests/ 下全部 *_test.R 与 lu_cluster_smoke.R，子进程逐个执行，
# PASS/FAIL 汇总写入 last_run.log（含时间戳、R 版本与失败套件的 FAIL 明细行）。
# Run every *_test.R and lu_cluster_smoke.R under tests/ as subprocesses,
# write a PASS/FAIL summary (timestamp, R version and per-suite FAIL detail
# lines) to last_run.log.
options(stringsAsFactors = FALSE)
# F-1（2026-09-26）：与 run_all.bat 的 locale 清理保持一致——测试脚本以
# 中文字面量 grep 源码，父进程继承的 LC_*（Git Bash / MSYS2 / CI 常见）会让
# readLines(encoding=) 与脚本字面量的编码标记不一致，产生假 FAIL。这里先清掉。
# Match run_all.bat's locale hygiene: inherited LC_* (common under Git Bash /
# MSYS2 / CI) can make encoding marks disagree and produce false FAILs.
Sys.unsetenv(c("LC_CTYPE", "LC_ALL", "LANG", "LC_COLLATE", "LC_MONETARY", "LC_TIME"))
# F-1 补充（2026-09-26）：父进程自身若在 LC_*=C.UTF-8 等 Windows 无法设置的
# locale 下启动（Git Bash 常见），R 回退为 C（ASCII）locale，此时 --file= 里的
# 中文路径会被 normalizePath 弄乱、last_run.log 写不出去。恢复系统默认 locale；
# 子进程则由上面的 unsetenv 保证干净。
# F-1 addendum (2026-09-26): if the harness itself starts under an unsettable
# locale (e.g. LC_*=C.UTF-8 in Git Bash), R falls back to the ASCII "C" locale
# and mangles Chinese characters in --file=, breaking the last_run.log write.
# Restore the system default locale; children are kept clean by the unsetenv above.
tryCatch({
  Sys.setlocale("LC_CTYPE", ""); Sys.setlocale("LC_COLLATE", "")
  Sys.setlocale("LC_MONETARY", ""); Sys.setlocale("LC_TIME", "")
}, error = function(e) NULL)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
tests <- sort(list.files(sd, pattern = "(_test\\.R|lu_cluster_smoke\\.R)$",
                         full.names = TRUE))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
res <- character(0); npass <- 0L; nfail <- 0L
detail <- character(0)
for (t in tests) {
  cat("==== ", basename(t), " ====\n", sep = "")
  out <- suppressWarnings(system2(rscript, shQuote(t),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  # F-1：崩溃/空输出不得判 PASS——即使退出码为 0，输出也必须非空
  # （套件结尾行各不一致：多数打"汇总：n PASS / m FAIL"，个别打
  # "全部通过/SMOKE PASS"，故以"有输出"为完成标记而非特定词）。
  # A suite that crashed or produced no output must not be marked PASS
  # (endings vary: most print "汇总: n PASS / m FAIL", a few print
  # "全部通过"/"SMOKE PASS", so non-empty output is the completion marker).
  ran <- length(out) > 0L
  ok <- identical(as.integer(status), 0L) && ran
  res <- c(res, sprintf("%s\t%s", basename(t), if (ok) "PASS" else "FAIL"))
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", basename(t)))
  if (!ok) {
    nfail <- nfail + 1L
    fl <- grep("FAIL", out, value = TRUE, fixed = TRUE)
    if (!length(fl)) fl <- tail(out, 6)
    detail <- c(detail, sprintf("# ---- %s ----", basename(t)), fl)
    cat(paste(tail(out, 8), collapse = "\n"), "\n")
  } else npass <- npass + 1L
}
log_lines <- c(
  paste("# run_all", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R.version.string),
  paste0("summary\t", npass, " PASS / ", nfail, " FAIL"),
  res)
if (length(detail)) log_lines <- c(log_lines, detail)
writeLines(log_lines, file.path(sd, "last_run.log"), useBytes = TRUE)
cat(sprintf("=== run_all 汇总：%d PASS / %d FAIL（详见 last_run.log）===\n",
            npass, nfail))
if (nfail > 0L) quit(status = 1L)
