# =============================================================================
# run_telemetry.R - 运行遥测记录脚本（用户口径 2026-09-21：软件质控之二）
#
# 记录每个运行步骤：起止时间、耗时、CPU 峰值、内存峰值、步骤内外存占用峰值、
# 结束后外存增加值。
#
# 用法（三态，独立进程调用；与 user_run.R 的接线在批次空闲后统一完成）：
#   Rscript run_telemetry.R start  <结果目录>          # 启动后台采样器（5 秒/次）
#   Rscript run_telemetry.R mark   <步骤标签> [结果目录] # 追加步骤边界标记（目录缺省取 RESULT_DIR）
#   Rscript run_telemetry.R stop                      # 汇总：逐步骤聚合表落盘
#   Rscript run_telemetry.R attach                     # 已有采样文件时直接重汇总
#
# 采样内容（每次）：整机 CPU%、R/引擎进程树 RSS 合计与最大单进程、
#   结果盘剩余空间；mark 时另记结果目录当前磁盘占用。
# 输出：<结果目录>\telemetry_samples.tsv（原始采样）与 telemetry_summary.txt
#   （逐步骤：elapsed / CPU 峰值 / 内存峰值 / 外存占用峰值 / 步骤后外存增量）。
# 仅基础包 + PowerShell；不依赖 Python。
# =============================================================================
args <- commandArgs(TRUE)
mode <- if (length(args) >= 1) args[1] else ""
# 参数位次按模式区分：start <目录>；mark <步骤标签> [目录]；stop [目录]；attach [目录]
# mark 的目录缺省取 RESULT_DIR 环境变量（与管线接线时的用法一致）
# argument positions differ per mode: start <dir>; mark <step label> [dir];
# stop [dir]; attach [dir]; mark falls back to the RESULT_DIR env var
rdir <- if (mode == "mark") {
  if (length(args) >= 3) args[3] else Sys.getenv("RESULT_DIR")
} else if (length(args) >= 2) args[2] else Sys.getenv("RESULT_DIR")
if (!nzchar(rdir) || is.na(rdir)) rdir <- getwd()
dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
smpl_f <- file.path(rdir, "telemetry_samples.tsv")

ps1 <- function(cmd) {
  r <- suppressWarnings(system2("powershell", c("-NoProfile","-Command", cmd),
                                stdout = TRUE, stderr = TRUE))
  r <- r[nzchar(trimws(r))]
  if (length(r)) trimws(r[1]) else NA_character_
}
cpu_pct <- function() suppressWarnings(as.numeric(ps1(
  "(Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average")))
free_mem_gb <- function() suppressWarnings(as.numeric(ps1(
  "[math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB,2)")))
tot_mem_gb <- function() suppressWarnings(as.numeric(ps1(
  "[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,2)")))
# 管线相关进程（R/引擎）树内存：GB，两值＝合计/最大单进程
proc_mem <- function() ps1(paste0(
  "$p = Get-Process Rscript,R,python,Searcher,sage,comet.win64,java,MaxQuantCmd,pNovo3_Search,Morpheus -ErrorAction SilentlyContinue; ",
  "if ($p) { $s = ($p | Measure-Object WorkingSet64 -Sum).Sum/1GB; ",
  "$m = ($p | Measure-Object WorkingSet64 -Maximum).Maximum/1GB; ",
  "'{0} {1}' -f [math]::Round($s,3), [math]::Round($m,3) } else { '0 0' }"))
disk_free_gb <- function(path) suppressWarnings(as.numeric(ps1(sprintf(
  "[math]::Round((Get-PSDrive '%s').Free/1GB,2)",
  substr(normalizePath(path, winslash="/"), 1, 1)))))
dir_size_gb <- function(path) suppressWarnings(as.numeric(ps1(sprintf(
  "[math]::Round((Get-ChildItem -LiteralPath '%s' -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum/1GB,3)",
  gsub("'", "''", normalizePath(path))))))

if (mode == "start") {
  tm <- suppressWarnings(as.numeric(Sys.getenv("TELEMETRY_INTERVAL", "5")))
  if (is.na(tm) || tm <= 0) tm <- 5   # 空/非法值兜底 5 秒（防 Sys.sleep(NA) 首拍即死）
  writeLines(c("ts\tevent\tstep\tcpu_pct\tfree_mem_gb\ttot_mem_gb\tproc_rss_sum\tproc_rss_max\tdisk_free_gb\tdir_size_gb"),
             smpl_f, useBytes = TRUE)
  tot <- tot_mem_gb()
  while (TRUE) {
    pm <- proc_mem(); pmv <- strsplit(pm, " ")[[1]]
    line <- sprintf("%s\tsample\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getenv("TELEMETRY_STEP", ""),
      cpu_pct(), free_mem_gb(), tot,
      if (length(pmv) >= 1) pmv[1] else NA,
      if (length(pmv) >= 2) pmv[2] else NA,
      disk_free_gb(rdir), dir_size_gb(rdir))
    cat(line, "
", file = smpl_f, append = TRUE, sep = "")
    Sys.sleep(tm)
  }
} else if (mode == "mark") {
  step <- if (length(args) >= 2) args[2] else "step"
  pm <- proc_mem(); pmv <- strsplit(pm, " ")[[1]]
  line <- sprintf("%s\tmark\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step,
    cpu_pct(), free_mem_gb(), tot_mem_gb(),
    if (length(pmv) >= 1) pmv[1] else NA,
    if (length(pmv) >= 2) pmv[2] else NA,
    disk_free_gb(rdir), dir_size_gb(rdir))
  cat(line, "
", file = smpl_f, append = TRUE, sep = "")
  cat("marked:", step, "\n")
} else if (mode %in% c("stop", "attach")) {
  if (!file.exists(smpl_f)) stop("telemetry_samples.tsv 不存在 / missing (先 start)")
  d <- read.delim(smpl_f, sep = "\t", stringsAsFactors = FALSE, na.strings = c("NA",""),
                  check.names = FALSE, quote = "")
  mk <- which(d$event == "mark")
  steps <- data.frame()
  for (i in seq_along(mk)) {
    a <- mk[i]
    b <- if (i < length(mk)) mk[i + 1] - 1 else nrow(d)
    seg <- if (b >= a + 1) d[(a + 1):b, , drop = FALSE] else d[0, , drop = FALSE]
    if (nrow(seg) == 0) seg <- d[a, , drop = FALSE]   # 相邻 mark（无采样间隔）不并入下段
    steps <- rbind(steps, data.frame(
      step = d$step[a],
      start = d$ts[a],
      elapsed_s = max(0, as.numeric(difftime(as.POSIXct(d$ts[min(b, nrow(d))]),
                                             as.POSIXct(d$ts[a]), units = "secs"))),
      cpu_peak_pct = if (nrow(seg)) suppressWarnings(max(seg$cpu_pct, na.rm = TRUE)) else NA_real_,
      cpu_avg_pct = if (nrow(seg)) suppressWarnings(round(mean(seg$cpu_pct, na.rm = TRUE), 1)) else NA_real_,
      proc_rss_sum_peak_gb = if (nrow(seg)) suppressWarnings(max(seg$proc_rss_sum, na.rm = TRUE)) else NA_real_,
      proc_rss_max_peak_gb = if (nrow(seg)) suppressWarnings(max(seg$proc_rss_max, na.rm = TRUE)) else NA_real_,
      dir_size_peak_gb = if (nrow(seg)) suppressWarnings(max(seg$dir_size_gb, na.rm = TRUE)) else NA_real_,
      dir_size_after_gb = suppressWarnings(min(tail(seg$dir_size_gb, 1), d$dir_size_gb[min(b, nrow(d))], na.rm = TRUE)),
      disk_free_min_gb = if (nrow(seg)) suppressWarnings(min(seg$disk_free_gb, na.rm = TRUE)) else NA_real_))
  }
  outf <- file.path(rdir, "telemetry_summary.txt")
  con <- file(outf, open = "wt", encoding = "UTF-8")
  writeLines(c("# 运行遥测汇总 / run telemetry summary",
               sprintf("# 结果目录 / result dir: %s", rdir),
               sprintf("# 采样行数 / samples: %d；步骤标记 / marks: %d", nrow(d), length(mk)), ""),
    con, useBytes = TRUE)
  if (nrow(steps)) write.table(steps, con, sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(c("", "# 步骤外存占用峰值＝dir_size_peak；结束后外存增加值＝末步 dir_size_after − 首步前基线（见首行采样）"), con, useBytes = TRUE)
  close(con)
  if (nrow(steps)) print(steps, row.names = FALSE, digits = 3)
  cat("written:", outf, "\n")
} else {
  cat("用法 / usage:\n  Rscript run_telemetry.R start <result_dir>\n",
      "  Rscript run_telemetry.R mark <step_label>\n",
      "  Rscript run_telemetry.R stop\n  Rscript run_telemetry.R attach\n")
}
