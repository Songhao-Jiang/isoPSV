# =============================================================================
# unlabeled_extract.R - unlabeled isoPSV CLI driver for steps 1-2:
#   step 1: resource assessment (cores / memory / MGF sizes -> max_cores,
#           printed to the console and used for the reads; not written back
#           into the parameter file)
#   step 2: per-input-MGF unlabeled extraction (strict subset by the regular
#           search result's Label==1 titles; same raw+scan patch blocks with
#           TITLE rewritten; one unlabeled MGF per input MGF; no peptide.txt)
#
# Usage:
#   Rscript unlabeled_extract.R
# with the parameter file path in the environment variable CONFIG (Chinese
# paths must not go through the command line on Windows).
# =============================================================================
#
# unlabeled_extract.R —— 无标记 isoPSV 的第一步/第二步命令行驱动：
#   第一步：资源评估（核数 / 内存 / MGF 体积 -> max_cores，
#           打印到控制台并用于读取；不会写回参数文件）
#   第二步：对每个输入 MGF 做无标记提取（按常规搜库结果中 Label==1
#           的谱题名做严格子集提取；同一 raw+scan 分块仅改写 TITLE；
#           每个输入 MGF 输出一个无标记 MGF；不产出 peptide.txt）
#
# 用法：
#   Rscript unlabeled_extract.R
# 参数文件路径放在环境变量 CONFIG 中（Windows 上中文路径不能走命令行）。

# 载入库文件（兼容层 + MGF 读取模块）
# Load the library files (compatibility layer + MGF reader)
# 核心库统一入口按 --file= 定位代码目录、绝对路径加载
# （boot.R），替代裸相对 source——此前仅当启动 cwd 恰为 code/ 时成立。
# core-library entry locate the code dir via --file= and
# load by absolute path (boot.R), replacing the bare relative source().
sd0 <- dirname(tryCatch(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"),
    error = function(e) ""))
if (is.na(sd0) || !dir.exists(sd0)) sd0 <- getwd()
source(file.path(sd0, "boot.R"))

# 从环境变量 CONFIG 读取参数文件路径并解析该文件
# Read the parameter-file path from the CONFIG environment variable and parse that file
cfg_path <- Sys.getenv("CONFIG")
if (!nzchar(cfg_path)) stop("set CONFIG to the parameter file path")
cfg <- read_param_file(cfg_path)
# candidate_scope=top5第一步 psm_file 若为需转换的
# 非表格格式（X!Tandem XML / Mascot CSV / PEAKS CSV），转换适配器经
# ISOPSV_CAND_ALL 环境变量输出全部候选（第二步合并步封顶 5）
cs1 <- if (!is.null(cfg$candidate_scope) && length(cfg$candidate_scope) > 0) tail(cfg$candidate_scope, 1) else "filtered"
if (identical(cs1, "top5")) {
  Sys.setenv(ISOPSV_CAND_ALL = "1")
}

# 解析 psm_file 与 MGF 列表；常规搜库结果缺失时立即报错
# Resolve the psm_file and the MGF list; stop immediately if the regular search result is missing
psm_file <- resolve_cfg_paths(cfg, "psm_file")[1]
# 全场景谱图入口检测（detect_spectra.R）
# 文件夹/文件列表 × raw/mgf/混合 × MaxQuant raw 警告＋自动传递
source(file.path(getwd(), "detect_spectra.R"))
mgf_entries <- resolve_cfg_paths(cfg, "mgf")
if (length(mgf_entries) && any(nzchar(trimws(mgf_entries)))) {
  ds <- detect_spectra_entries(mgf_entries, cfg)
  mgf_files <- ds$mgf_files
} else {
  d0 <- dirname(psm_file)
  mgs <- gsub("\\\\", "/", list.files(d0, pattern = "\\.mgf$",
                ignore.case = TRUE, recursive = TRUE, full.names = TRUE))
  cat("[第一步] mgf 键为空：按 psm_file 所在文件夹自动查找 ", d0,
      "（mgf ", length(mgs), " 个）
", sep = "")
  if (length(mgs)) {
    mgf_files <- sort(mgs)
  } else {
    raws0 <- gsub("\\\\", "/", list.files(d0, pattern = "\\.raw$",
                  ignore.case = TRUE, recursive = TRUE, full.names = TRUE))
    mgf_files <- if (length(raws0)) detect_spectra_entries(raws0, cfg)$mgf_files else character(0)
  }
  if (!length(mgf_files))
    stop("[第一步] mgf 键为空，且 psm_file 所在文件夹内未找到 .mgf 或 .raw：", d0,
         "
请在 mgf 键填写质谱 raw/MGF 文件或文件夹。
")
}
if (!file.exists(psm_file)) stop("psm_file not found: ", psm_file)
cat("==== [用户第一步] 无标记谱图提取：常规搜库结果 -> unlabeled 谱图 ====\n")
cat("==== [User Step 1] unlabeled-spectrum extraction: regular search result -> unlabeled spectra ====\n")
cat("psm_file:", psm_file, "\n")
cat("mgf files:", length(mgf_files), "\n")

# max_cores 在运行时自动评估（打印到控制台）并用于读取，
# 但不会写回参数文件；参数文件中已生效的 max_cores 行
# 仍可作为手动覆盖值。
# max_cores is assessed at run time (printed to the console) and used for the
# reads, but NOT written back into the parameter file; an active max_cores
# line in the parameter file still works as a manual override.
max_cores <- if (!is.null(cfg$max_cores) && length(cfg$max_cores) > 0) {
  mc <- suppressWarnings(as.integer(tail(cfg$max_cores, 1)))
  if (!is.na(mc) && mc >= 1L) mc else assess_resources(mgf_files)
} else {
  assess_resources(mgf_files)
}

out_dir <- file.path(cfg$config_dir, "unlabeled")
cat("max_cores:", max_cores, "| out:", out_dir, "\n")

# 第一步各阶段双语计时
# bilingual timing marks for the step-1 phases
pt0 <- proc.time()
# 函数：pt_mark —— 打印某一阶段的双语计时标记。
# 作用：计算自上一标记以来经过的时间，并更新全局计时起点 pt0。
# 输入：tag = 阶段名称字符串。
# 输出：无返回值；仅在控制台打印一行双语耗时信息。
# Function: pt_mark -- print the bilingual timing mark of one phase.
# Purpose: compute the elapsed time since the previous mark and update the global timing start pt0.
# Inputs: tag = the phase name string.
# Output: no return value; only prints one bilingual elapsed-time line to the console.
pt_mark <- function(tag, tag_en = tag) {
  pt1 <- proc.time()
  # 计算自上一标记以来的墙钟与 CPU 耗时
  # Compute the wall-clock and CPU time elapsed since the previous mark
  d <- pt1 - pt0
  cat(sprintf("[计时] %s：已用 %.1fs / CPU %.1fs\n", tag, d[3], d[1]))
  cat(sprintf("[timing] %s: elapsed %.1fs / CPU %.1fs\n", tag_en, d[3], d[1]))
  # 打印双语耗时并重置全局计时起点
  # Print the bilingual elapsed time and reset the global timing start point
  pt0 <<- pt1
}

# 对每个输入 MGF 做无标记谱图提取，先输出到临时目录（逐输入文件名与重切目标
# <raw>_HCDFT.mgf 可能同名，直接写 unlabeled 会边读边写同一文件导致死循环）
# Extract unlabeled spectra for every input MGF into a TEMP directory first (the
# per-input file names may equal the re-split target <raw>_HCDFT.mgf; writing both
# into unlabeled/ would read and write the same file, an endless loop)
tmp_dir <- file.path(cfg$config_dir, "unlabeled_tmp")
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
#：目录创建失败不再静默——G 盘不存在等场景此前 showWarnings=FALSE
# 吞掉真因，写文件时才报 cannot open the connection（用户无法定位）
if (!dir.exists(tmp_dir))
  stop("无法创建结果目录：", tmp_dir,
       "（请检查 result_dir 路径是否存在/可写）
",
       "cannot create result directory: ", tmp_dir,
       " (check that result_dir exists and is writable)", call. = FALSE)
# 非表格格式（X!Tandem XML 等）自动经对应适配器转换为交换格式后再提取
# 
# P1-3转换产物落在本临时目录内——此调用
# 必须放在 unlink+dir.create 之后、读取之前（此前放在清目录之前，刚转好的
# 文件被随后 unlink 删掉，xtandem/mascot/peaks 三条可转路径全废，报错
# cannot open the connection 且不指向真因）
psm_file <- psm_prepare_input(psm_file, cfg, sd0)
# （无标记通道=0）：不按 Label 过滤，输入结果文件的全部鉴定谱都提取
ul_ch <- if (!is.null(cfg$unlabeled_channel) && length(cfg$unlabeled_channel) > 0)
  tail(cfg$unlabeled_channel, 1) else "1"
summ <- extract_unlabeled_per_mgf(psm_file, mgf_files, tmp_dir, max_cores,
                                  label_type = if (!is.null(cfg$label_type) &&
                                      length(cfg$label_type) > 0)
                                    tolower(tail(cfg$label_type, 1)) else "unlabeled",
                                  psm_all_unlabeled = if (!is.null(cfg$psm_all_unlabeled) &&
                                      length(cfg$psm_all_unlabeled) > 0)
                                    identical(tail(cfg$psm_all_unlabeled, 1), "1") else FALSE,
                                  keep_all = identical(ul_ch, "0"))
print(summ, row.names = FALSE)
pt_mark("提取无标记谱图", "extract unlabeled spectra")

# mgf 输出只发生在本步：把临时目录中的逐输入文件按标题 raw 前缀重切成逐 raw 纯文件
# （unlabeled\<raw>_HCDFT.mgf），随后删除临时目录——第二步起只读、不再输出任何 mgf。
# The mgf outputs happen ONLY in this step: the per-input intermediates in the temp
# dir are re-split by title raw prefix into pure per-raw files (unlabeled\<raw>_HCDFT.mgf),
# then the temp dir is removed - from step 2 on the mgfs are only read, never written.
out_files <- sort(list.files(tmp_dir, pattern = "\\.mgf$", full.names = TRUE))
if (length(out_files) == 0) {
  unlink(tmp_dir, recursive = TRUE)
  stop("第一步提取未产出任何无标记 MGF（psm_file 中无 Label 首字符为 1 的谱图）\nstep 1 produced no unlabeled MGFs (no spectra with Label first char == 1 in the psm_file)")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
per_raw_files <- write_per_raw_mgfs(out_files, out_dir)
unlink(tmp_dir, recursive = TRUE)
cat("逐 raw 纯无标记 MGF：", length(per_raw_files), " 个文件\n")
cat("per-raw pure unlabeled MGFs: ", length(per_raw_files), " files\n")
pt_mark("逐 raw 重切", "per-raw re-split")

# 写第一步完成标记（总程序以它作为跳过检查：存在即跳过第一步）
# write the step-1 completion marker (the master program skips step 1 when it exists)
writeLines(c("isoPSV step 1 done",
             paste0("files: ", length(per_raw_files)),
             paste0("at: ", Sys.time())),
           file.path(cfg$config_dir, "unlabeled", "step1_done.txt"))
cat("完成标记：unlabeled/step1_done.txt\n")
cat("completion marker: unlabeled/step1_done.txt\n")

# 自动生成第二步参数文件（search_file 占位；保留用户已填条目）
# auto-generate the user-step-2 parameter file (search_file placeholders;
# existing user-filled search_file entries are preserved)
make_step2_config(cfg_path, cfg,
                  file.path(cfg$config_dir, "config_step2.txt"))
cat("==== [用户第一步] 完成，输出目录：", out_dir, " ====\n")
cat("==== [User Step 1] done, output: ", out_dir, " ====\n")
cat("DONE\n")
