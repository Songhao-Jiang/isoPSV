# =============================================================================
# stageB_paired_match.R - 阶段 B 驱动程序：成对标记谱匹配、
# stageB_paired_match.R - stage B driver: paired labeled-spectrum matching,
# 打分、筛选与错误归因，针对五策略候选在序列或修饰上
# scoring, screening and error attribution for the spectra whose five-strategy
# 不一致的谱（isoPSV 阶段 B）。
# candidates disagree in sequence or modification (isoPSV stage B).
#
# 用法：
# Usage:
#   CONFIG3=<第三步参数文件> ALL_MGF=<原始 MGF 目录> \
#   CONFIG3=<step-3 parameter file> ALL_MGF=<original-MGF folder> \
#     Rscript stageB_paired_match.R
#     Rscript stageB_paired_match.R
# 标记质量模型取自第一步参数文件的 label_type
# The label mass model comes from the step-1 parameter file's label_type
#；输出写到 <config dir>/stageB/。
# (N15 / C13 / SILAC); outputs go to <config dir>/stageB/.
# =============================================================================

# ---- 首先定位本脚本所在文件夹 ---------------------------------------------
# ---- locate this script's folder first -------------------------------------
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
if (is.na(script_arg)) script_arg <- args0[grepl("\\.R$", args0)][1]
script_arg <- sub("^--file=", "", script_arg)
setwd(dirname(normalizePath(script_arg, winslash = "/")))

# 加载通用工具与核心计算模块（ACS 核心、成对参数、MGF 读取等）
# Load the shared utilities and the core computation modules (ACS core, paired params, MGF reader, etc.)
source("compat.R")
source("ReadMGF_script.R")

pt0 <- proc.time()
# pt_mark：打印某个阶段的耗时标记（双语），并重置计时起点
# 作用：记录并打印从上次标记到当前的墙钟/CPU 耗时，便于跟踪各阶段运行时间
# 输入：tag——阶段名称（字符串，输出时作为标记打印）
# 输出：无返回值；副作用为打印耗时行并更新外层计时变量 pt0
# pt_mark: print an elapsed-time mark (bilingual) for a phase and reset the timer
# Purpose: record and print wall/CPU time since the previous mark, for tracking per-phase runtime
# Input: tag - the phase name (character, printed as the mark label)
# Output: no return value; side effects are the printed timing line and the updated outer pt0
pt_mark <- function(tag, tag_en = tag) {
  # 计算与上一标记之间的耗时，打印后把计时起点推进到当前
  # Compute the elapsed time since the previous mark, print it, then advance the timer origin
  pt1 <- proc.time()
  d <- pt1 - pt0
  cat(sprintf("[计时] %s：已用 %.1fs / CPU %.1fs\n", tag, d[3], d[1]))
  cat(sprintf("[timing] %s: elapsed %.1fs / CPU %.1fs\n", tag_en, d[3], d[1]))
  pt0 <<- pt1
}

# 读取第三步参数文件与基础配置，确定标记类型与原始 MGF 目录
# Read the step-3 parameter file and the base config; determine the label type (default N15) and the raw MGF folder
cfg3 <- read_param_file(if (nzchar(Sys.getenv("CONFIG3"))) Sys.getenv("CONFIG3") else
                           file.path(getwd(), "config_step3.txt"))
cfg_dir <- cfg3$config_dir
if (is.null(cfg_dir) || length(cfg_dir) == 0 || !nzchar(cfg_dir[1]))   # 审查 P0-3/P3-6：config_dir 缺失/为空时给可读报错
  stop("config_step3 的 config_dir 缺失或为空，请检查配置（或删除 result_* 后重跑）
config_dir missing/empty in config_step3")
# 结果文件夹（RESULT_DIR 由总程序注入；独立运行按时间戳定位/新建）
# result folder (RESULT_DIR injected by the master; standalone runs locate/create by timestamp)
result_dir <- get_result_dir(cfg_dir)
base_cfg <- read_param_file(resolve_cfg_paths(cfg3, "base_config")[1])
# v1.0label_type 为多行键（每行一个标记）。旧单值
# 口径走 legacy 路径、行为逐字节不变；新多行口径：
# 类 1 全标记作理论离子标签，类 2 修饰标记清单经 labels 追加进
# 修饰串（paired_match_candidate 内按适用残基/N 端展开）。
label_spec <- parse_label_spec(base_cfg$label_type)
if (!is.null(label_spec$legacy)) {
  label_type <- label_spec$legacy
  labels_b <- NULL
} else {
  label_type <- if (!is.null(label_spec$class1)) label_spec$class1 else "Unlabel"
  labels_b <- label_spec$mods
}
label_name <- if (!is.null(labels_b) && length(labels_b) > 0) {
  paste(tolower(c(label_spec$class1, labels_b)), collapse = "_")
} else tolower(label_type)
all_mgf_dir <- Sys.getenv("ALL_MGF")
if (!nzchar(all_mgf_dir)) stop("set ALL_MGF to the original-MGF folder")
# 建立阶段 B 输出目录 stageB/
# Create the stage B output directory stageB/
out_dir <- file.path(result_dir, "2.labeled_result")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 打印本阶段运行信息：标记类型、原始 MGF 目录与输出目录
# Print the run info for this stage: label type, raw MGF folder and output directory
cat("==== [阶段 B] 成对标记谱图匹配 ====\n")
cat("==== [Stage B] paired labeled-spectrum matching ====\n")
cat("label_type:", label_type, "| ALL_MGF:", all_mgf_dir, "| out:", out_dir, "\n")

# 读入五策略合并结果、pDeep2 预测（索引＋强度）与氨基酸质量表
# Read the five-strategy merge result, the pDeep2 predictions (index + intensities) and the amino-acid mass table
merged <- read.delim(resolve_cfg_paths(cfg3, "merge_result")[1], sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE,
                     na.strings = character(0), colClasses = "character")
# colClasses="character"某策略列全空时 read.table 会把它定成
# logical 型、空值变 NA——含该策略候选的行在 by_ions_position 里
# if(modification != "") 崩溃（576 运行 f 槽候选全崩、整块连坐的根因）
# use_pdeep2=0 时构造空预测表：预测余弦全 NA（缺省处理：不填 0、标注、投票降级）
# when use_pdeep2=0 an empty prediction table is built: all predict cosines are
# NA (treated as missing: never 0, flagged, dropped from the votes)
use_pd2 <- if (!is.null(cfg3$use_pdeep2) && length(cfg3$use_pdeep2) > 0) tail(cfg3$use_pdeep2, 1) else "1"
if (identical(use_pd2, "0")) {
  pred_idx <- setNames(integer(0), character(0))
  pred_int <- list()
  cat("use_pdeep2=0：预测指标按缺省处理\n")
  cat("use_pdeep2=0: predict metric treated as missing\n")
} else {
  pred <- read_predict_txt(resolve_cfg_paths(cfg3, "predict_file")[1])
  pred_idx <- pred$idx
  pred_int <- pred$intensities
}
aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# 目标谱集：五策略候选在序列或修饰上不一致的行
# target spectra: five-strategy candidates disagree in sequence or modification
# 按 Peptide_state / Modification_state 不一致筛选目标谱行，并打印数量
# Select the target rows where Peptide_state / Modification_state disagree, and print the counts
is_disagree <- merged$Peptide_state != "Same" | merged$Modification_state != "Same"
cat("五策略不一致谱行：", sum(is_disagree), "/", nrow(merged), "\n")
cat("disagreeing spectra: ", sum(is_disagree), "/", nrow(merged), "\n")
#：阶段 B 对全部谱做原始鉴定配对——一致谱仅 a 候选（combine 对
# Same 状态返回 "a"），使主表的成对状态覆盖每一行。设计约定：筛选级联第 2 步
# 的配对证据是"原始 pFind 结果自身配对（cos>0.7）"，与竞争仲裁无关；此前仅
# 不一致谱入配对报告，主表其余行成对状态为占位符，与逐谱配对文件不一致。
#: stage B pair-checks EVERY spectrum with its original
# identification (consistent rows carry only the "a" candidate, which is what
# combine() returns for "Same" states), so the master table's paired state
# covers every row. The cascade step-2 evidence is the ORIGINAL pFind result's
# own pairing (cos > 0.7), independent of the competitive arbitration;
# previously only disagreeing spectra entered the paired report and the master
# left placeholders for the rest - inconsistent with the per-spectrum pairing.
target <- merged

# 定位目标表中各策略的 Sequence / Modification 列位置（后续按策略号取列用）
# Locate the Sequence / Modification column positions per strategy in the target table (used later for per-strategy access)
seq_cols <- which(grepl("\\.Sequence$", colnames(target)))
mod_cols <- which(grepl("\\.Modification$", colnames(target)))
# 多引擎版：槽位字母（列名前缀）随 worker 导出，供状态串字母→槽列定位
# multi-engine edition: the slot letters (column prefixes) are exported to the
# workers for the state-letter -> slot-column lookup
slot_syms <- sub("\\.Sequence$", "", colnames(target)[seq_cols])

# 目标谱对应的 raw 名与原始 MGF 文件的对应关系
# raw stems of the target spectra and the corresponding original MGFs
# 由谱标题解析 raw 名，枚举原始 MGF 目录中的 mgf，建立 raw→mgf 映射
# Parse the raw names from the spectrum titles, list the mgfs under the raw MGF folder and build the raw->mgf mapping
raws <- unique(vapply(target$File_Name, function(t) title_raw_scan(t)["raw"],
                      character(1)))
mgf_files <- sort(list.files(all_mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
if (length(mgf_files) == 0) {
  cat("警告：原始 MGF 目录为空或不存在：", all_mgf_dir, "——全部谱图将不参与成对核查\n")
  cat("warning: raw MGF folder empty or missing: ", all_mgf_dir, " - no spectra will be pair-checked\n")
}
mgf_stems <- vapply(mgf_files, function(f) sub("_HCDFT$", "",
                    sub("\\.mgf$", "", basename(f))), character(1))
need_idx <- match(raws, mgf_stems)
# 检查缺漏：缺少对应 MGF 的 raw 打印提示，其目标谱不参与成对核查
# （下游冲突报告将这些谱标为 not pair-checked，不会静默丢失）
# Check for gaps: print the raws whose MGF is missing; their target spectra
# skip the paired check (the conflict report marks them not pair-checked, so
# no row is silently lost)
if (any(is.na(need_idx))) {
  cat("缺少原始 MGF 的 raw：", paste(raws[is.na(need_idx)], collapse = ", "), "\n")
  cat("raws missing their MGFs: ", paste(raws[is.na(need_idx)], collapse = ", "), "\n")
}
need_idx <- need_idx[!is.na(need_idx)]
cat("需要加载的原始 MGF：", length(need_idx), "\n")
cat("raw MGFs to load: ", length(need_idx), "\n")

# 确定工作进程数（检测核数-1，上限 16），创建 PSOCK 集群并在退出时关闭
# Determine the worker count (detected cores - 1, capped at 16); create a PSOCK cluster and stop it on exit
nc <- parallel::detectCores()
if (is.na(nc) || nc < 2L) nc <- 2L
use_cores <- max(1L, min(nc - 1L, 16L))
cl <- parallel::makeCluster(use_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
# 统一 bootstrap（绝对路径 source + 健康检查），
# worker 错误处理器单独设置
ispsv_bootstrap_cluster(cl)
parallel::clusterEvalQ(cl, {
  options(error = function() {
    #：worker stdout 会被丢弃，堆栈改写入 worker_error_log.txt
    # 
    tb <- paste(vapply(sys.calls(), function(cc) {
      h <- deparse(cc[[1]], nlines = 1L)
      if (length(h) == 0) "?" else substr(h[1], 1, 80)
    }, character(1)), collapse = " -> ")
    msg <- paste0("WORKER-FATAL: ", conditionMessage(geterrmessage()))
    try(write(c(msg, tb), file.path(get("out_dir", envir = .GlobalEnv),
                                    "worker_error_log.txt"), append = TRUE),
        silent = TRUE)
    message("WORKER-ERROR-TRACEBACK-START")
    message(tb)
    message("WORKER-ERROR-TRACEBACK-END")
    q("no", 1, FALSE)
  })
  invisible(NULL)
})
# 原始 MGF 谱缓存：缓存存在且 raw 集合＋文件 mtime 一致时复用否则重新加载并保存
# Raw MGF cache: reuse when the raw set AND file mtimes match (B2 optimization
#: file-freshness validation added so a modified MGF invalidates the
# cache); otherwise reload and save
loaded <- NULL
cache_file <- file.path(out_dir, paste0("raw_loaded_", label_name, ".rds"))
cache_sig_file <- file.path(out_dir, paste0("raw_loaded_", label_name, ".sig"))
if (file.exists(cache_file) && file.exists(cache_sig_file)) {
  tmp <- try(readRDS(cache_file), silent = TRUE)
  sig_old <- try(readLines(cache_sig_file, warn = FALSE), silent = TRUE)
  sig_new <- paste(mgf_files[need_idx],
                   format(file.mtime(mgf_files[need_idx]), "%Y%m%d%H%M%S"),
                   sep = ":", collapse = "|")
  if (!inherits(tmp, "try-error") && !inherits(sig_old, "try-error") &&
      identical(sig_old, sig_new) && setequal(names(tmp), mgf_stems[need_idx])) {
    loaded <- tmp
    cat("raw cache loaded (freshness verified):", length(loaded), "files\n")
  } else {
    cat("raw cache stale (MGF changed) - reloading\n")
  }
}
if (is.null(loaded)) {
  loaded <- load_mgf_shared(mgf_files[need_idx], cl)
  names(loaded) <- mgf_stems[need_idx]
  saveRDS(loaded, cache_file)
  writeLines(paste(mgf_files[need_idx],
                   format(file.mtime(mgf_files[need_idx]), "%Y%m%d%H%M%S"),
                   sep = ":", collapse = "|"), cache_sig_file, useBytes = TRUE)
}
pt_mark("加载原始 MGF", "load raw MGFs")
# 构建按 raw 分组的工作负载：目标表行＋该 raw 的谱对象与谱参数
# build per-raw payloads: rows of the target table + raw spectrum objects
# 将目标行按 raw 分组，逐个 raw 打包（行、谱对象、谱参数）供工作进程使用
# Split the target rows by raw and pack per-raw payloads (rows, spectrum objects, spectrum parameters) for the workers
raw_of <- vapply(target$File_Name, function(t) title_raw_scan(t)["raw"],
                 character(1))
raw_idx <- split(seq_len(nrow(target)), factor(raw_of, levels = unique(raw_of)))
payloads <- lapply(names(raw_idx), function(rw) {
  rows <- raw_idx[[rw]]
  objs <- loaded[[rw]]
  if (is.null(objs)) return(NULL)   # 缺失 raw：目标谱跳过成对核查（已打印警告）
  pars <- extractParameters_MGF(objs)
  list(raw = rw, rows = rows, datsub = target[rows,, drop = FALSE],
       pars = pars, objs = objs)
})
n_miss_raws <- sum(vapply(payloads, is.null, logical(1)))
payloads <- payloads[!vapply(payloads, is.null, logical(1))]
if (n_miss_raws > 0) {
  cat("跳过成对核查的 raw：", n_miss_raws, " 个（其目标谱由下游标为 not pair-checked）\n")
  cat("raws skipped in the paired check: ", n_miss_raws, " (their spectra are marked not pair-checked downstream)\n")
}
rm(loaded, target)
invisible(gc())

# P1 优化raw 内再分块（每块 ≤2000 行），充分利用全部 worker。
# 此前 7 个 raw = 7 个并行任务，16 个 worker 中 9 个空闲；分块后 ~35 个任务。
# P1 optimization split each raw's rows into chunks of <=2000
# to utilise all workers. Previously 7 raws = 7 parallel tasks with 9 of 16
# workers idle; chunking yields ~35 tasks.
CHUNK_N <- 2000L
chunked <- list()
for (p in payloads) {
  nr <- nrow(p$datsub)
  if (nr <= CHUNK_N) {
    chunked[[length(chunked) + 1]] <- p
  } else {
    nc <- ceiling(nr / CHUNK_N)
    for (ci in seq_len(nc)) {
      idx <- ((ci - 1L) * CHUNK_N + 1L):min(ci * CHUNK_N, nr)
      chunked[[length(chunked) + 1]] <- list(
        raw = p$raw, datsub = p$datsub[idx,, drop = FALSE],
        pars = p$pars, objs = p$objs)
    }
  }
}
payloads <- chunked
cat("[P1] 分块后 payload 数：", length(payloads), "\n")

# scan_worker：阶段 B 并行工作函数，处理一个 raw 的成对匹配任务
# 作用：对该 raw 下的每条目标谱，枚举五策略候选；逐候选做无标记离子匹配与成对标记谱匹配，
#       再按五指标多数投票选出最佳组，输出参考格式的一行结果
# 输入：pp——按 raw 打包的工作负载列表（raw、rows、datsub、pars、objs）
# 输出：该 raw 全部目标谱的结果数据框（无结果时为空数据框）
# scan_worker: the stage B parallel worker that processes one raw's paired-matching tasks
# Purpose: for every target spectrum of this raw, enumerate the five-strategy candidates; per candidate
#          run the unlabeled ion match and the paired labeled-spectrum match, then pick the best group
#          by the five-metric majority vote and emit one row in the reference format
# Input: pp - the per-raw payload list (raw, rows, datsub, pars, objs)
# Output: a data.frame of all target spectra of this raw (empty when there are no results)
scan_worker <- function(pp) {
  # P2 优化标题 O(1) 命名向量索引取代线性 which() 扫描
  # P2 optimization O(1) named-vector title index replaces
  # the linear which() scan
  title_ix <- setNames(seq_along(pp$pars$Title), pp$pars$Title)
  # v1.0RT 排序每 raw 分块算一次，随
  # 参数传入 paired_match_candidate（函数内挂 attr 因写时复制从不跨调用生效，
  # 每候选重排整 raw 参数表 O(R log R)）
  pp$rt_ord <- order(pp$pars$RTinSeconds, na.last = NA)
  # P3 优化worker 级理论离子缓存——同一 (seq, mod, label) 在
  # 预测余弦与 Xle 变体分支反复调用 by_ions_position，哈希缓存避免重算。
  # （性能第二轮）：缓存本体提升到库级进程缓存 cached_by_ions_position
  # （跨分块复用，此前每 chunk 重建即失效）；本闭包只保留 NULL 语义包装。
  # P3 optimization worker-level theoretical-ion cache - the
  # predicted-cosine and Xle-variant branches call by_ions_position repeatedly
  # for the same (seq, mod, label); a hash cache avoids recomputation.
  # (perf round 2): the cache itself moved to the library-level
  # cached_by_ions_position (reused across chunks - the old per-chunk env was
  # rebuilt every chunk); this closure keeps only the NULL-semantics wrapper.
  theo_cached <- function(seq, mod, lab) {
    v <- cached_by_ions_position(seq, mod, lab, aa, aa_tab)
    if (identical(v, "Not a peptide.") || is.null(v) || length(v) <= 1) return(NULL)
    v
  }
  # 取出该 raw 的目标行子表，逐行处理
  # Take the target-row subset of this raw and process row by row
  dat <- pp$datsub
  n <- nrow(dat)
  rep <- vector("list", n)
  #：行体提取为闭包，循环内 tryCatch 按行隔离（单谱崩溃不再
  # 连坐整个分块）；行体内 next 改为 return(NULL)
  row_fun <- function(j) {
    fn <- dat$File_Name[j]
    li <- unname(title_ix[fn])
    if (is.na(li)) {
      key <- title_scan_key(fn)
      li <- which(startsWith(pp$pars$Title, paste0(key, ".")))
      if (length(li) == 0) return(NULL)
      li <- li[1]
    }
    light_spec <- pp$objs[[li]]@spectrum
    light_rt <- pp$pars$RTinSeconds[li]
    z <- as.integer(dat$Charge[j])
    # 生产版候选分组：合并肽/修饰状态组，
    # 每组取一个代表性（序列、修饰）候选
    # production candidate grouping: combine the peptide/modification state
    # groups, one representative (sequence, modification) per group
    pms <- combine_peptide_modification_state(dat$Peptide_state[j],
                                              dat$Modification_state[j])
    state_split <- strsplit(pms, ";", fixed = TRUE)[[1]]
    # 逐状态组解析策略号并取出对应序列/修饰，非空的加入候选列表
    # Parse the strategy letter per state group and fetch the corresponding sequence/modification; add non-empty ones to the candidate list
    cand <- list()
    for (st in state_split) {
      if (!nzchar(st)) return(NULL)
      k1 <- substr(st, 1, 1)
      ci <- seq_cols[match(k1, slot_syms)]
      sq <- as.character(dat[j, ci])
      md <- as.character(dat[j, ci + 1])
      if (is.na(md)) md <- ""            # 空修饰列被读成 logical NA 的防御
      if (is.na(sq) || sq == "") return(NULL)
      cand[[length(cand) + 1]] <- list(idx = strsplit(st, "", fixed = TRUE)[[1]],
                                       seq = sq, mod = md)
    }
    if (length(cand) == 0) return(NULL)
    # P4 优化结果向量预分配＋下标赋值（取代逐次 c() 的 O(n²) 增长）
    # P4 optimization preallocated result vectors with indexed
    # assignment (replaces the O(n^2) growth of repeated c() appends)
    nc <- length(cand)
    states <- character(nc); ul_cov <- character(nc); ul_by <- character(nc)
    ul_all <- character(nc); lb_fn <- character(nc); lb_cov <- character(nc)
    lb_by <- character(nc); lb_all <- character(nc)
    dcos <- numeric(nc); prd <- numeric(nc)
    # 逐候选计算：pDeep2 预测余弦＋无标记离子匹配＋成对标记谱匹配
    # Per candidate: the pDeep2 predicted cosine + unlabeled ion matching + paired labeled-spectrum matching
    for (k in seq_along(cand)) {
      sq <- cand[[k]]$seq; md <- cand[[k]]$mod
      # 由序列|修饰|电荷查 pDeep2 预测索引；查不到时退化为不带修饰的键
      # Look up the pDeep2 prediction index by seq|mod|charge; fall back to the no-modification key when missing
      pi <- unname(pred_idx[paste0(sq, "|", md, "|", z)])
      if (is.na(pi)) pi <- unname(pred_idx[paste0(sq, "||", z)])
      if (is.na(pi)) {                                            # SAP 代入键
        fxk <- apply_sap_to_sequence(sq, md)
        pi <- unname(pred_idx[paste0(fxk$seq, "|", fxk$mod, "|", z)])
        if (is.na(pi)) pi <- unname(pred_idx[paste0(fxk$seq, "||", z)])
      }
      # 预测谱缺失时不填 0，记 NA（结果文件中以空值显式标注、并在 predict_missing 列标记）
      # a missing predicted spectrum is NA, never 0 (annotated empty in the outputs and
      # flagged by the predict_missing column)
      pc <- NA_real_
      pred_cos_one <- function(pj, sq_key, md_key) {              # 单预测条目余弦
        if (is.na(pj)) return(NA_real_)
        pv <- pred_int[[pj]]
        im0 <- try(ions_match(theo_cached(sq_key, md_key, "Unlabel"),
                              light_spec, 20, "ppm"), silent = TRUE)
        if (inherits(im0, "try-error")) return(NA_real_)
        v <- im0[["by_ions_long.match_intensity"]]
        nmin <- min(length(pv), length(v))
        if (nmin > 0) cosine(v[seq_len(nmin)], pv[seq_len(nmin)]) else NA_real_
      }
      if (!is.na(pi)) {
        pc <- pred_cos_one(pi, sq, md)
      }
      if (is.na(pc) && grepl("Xle->", md, fixed = TRUE)) {       # Xle 双条目取最高
        best <- NA_real_
        for (v_k in apply_sap_variants(sq, md)) {
          pj <- unname(pred_idx[paste0(v_k$seq, "|", v_k$mod, "|", z)])
          if (is.na(pj)) pj <- unname(pred_idx[paste0(v_k$seq, "||", z)])
          c_k <- pred_cos_one(pj, v_k$seq, v_k$mod)
          if (!is.na(c_k) && (is.na(best) || c_k > best)) best <- c_k
        }
        pc <- best
      }
      # 成对匹配核心：在无标记谱附近寻找标记谱并按标记质量模型打分
      # Paired-match core: find the labeled spectrum near the unlabeled one and score it with the label mass model
      pm <- paired_match_candidate(sq, md, z, light_spec, light_rt,
                                   pp$pars, pp$objs, label_type,
                                   aa = aa, aa_tab = aa_tab,
                                   rt_ord = pp$rt_ord,
                                   labels = labels_b)
      # 无成对结果时按 Single 记录空值（缺失不填 0：标记侧指标与成对余弦记 NA，
      # 输出为空值并由 pair_missing 列标记）；有成对结果时记录双谱打分与余弦
      # When the paired match is missing, record empty Single values (missing stays
      # missing: the label-side metrics and the paired cosine are NA, shown empty and
      # flagged by pair_missing); otherwise record the dual-spectrum scores and cosines
      if (is.null(pm)) {
        states[k] <- "single"; ul_cov[k] <- ""              # 二值小写
        ul_by[k] <- ""; ul_all[k] <- ""
        lb_fn[k] <- ""; lb_cov[k] <- ""; lb_by[k] <- ""
        lb_all[k] <- ""; dcos[k] <- NA_real_; prd[k] <- round(pc, 4)
      } else {
        states[k] <- if (pm$paired) "paired" else "single"  # 二值小写
        ul_cov[k] <- pm$light$coverage
        ul_by[k] <- pm$light$PSM_score_by
        ul_all[k] <- pm$light$PSM_score_all
        lb_fn[k] <- pm$label_title
        lb_cov[k] <- if (pm$paired) pm$label$coverage else ""
        lb_by[k] <- if (pm$paired) pm$label$PSM_score_by else ""
        lb_all[k] <- if (pm$paired) pm$label$PSM_score_all else ""
        dcos[k] <- pm$paired_cosine
        prd[k] <- round(pc, 4)
      }
    }
    # 五指标多数投票（按谱降级）：任一候选缺某指标（无配对→标记侧指标缺；无预测→预测
    # 余弦缺），则该指标维度本谱弃用；缺失以空值显式标注、不填 0；并列全部保留
    # five-metric majority vote with per-spectrum degradation: a metric missing for
    # ANY candidate (no pair -> label-side metrics; no prediction -> predict cosine)
    # is dropped from this spectrum's vote; missing stays empty, never 0; ties kept
    vm <- list(ul_by = as.numeric(ul_by), lb_by = as.numeric(lb_by),
               ul_all = as.numeric(ul_all), lb_all = as.numeric(lb_all),
               predict = prd)
    use <- !vapply(vm, function(v) any(is.na(v)), logical(1))
    winners <- integer(0)
    for (v in vm[use]) {
      winners <- c(winners, which(v == max(v, na.rm = TRUE)))
    }
    if (length(winners) == 0) winners <- 1L
    tb <- table(winners)
    bg <- as.integer(names(tb)[tb == max(tb)])
    bg_s <- paste(bg, collapse = ";")
    bg1 <- bg[1]
    # 汇总该谱的参考格式结果行（候选、匹配状态、各指标、最佳组、缺失标记与投票指标等）
    # Assemble the reference-format result row for this spectrum (candidates, match
    # state, per-metric values, best group, missing flags and the voting metrics)
    rep[[j]] <- data.frame(
      File_Name = fn,
      candidates = paste(vapply(cand, function(x) {
        fx <- normalize_terminal_loss(x$seq, x$mod)  # 输出格式规范化（内部匹配不变）
        paste(fx$s, if (nzchar(fx$m)) paste0("(", fx$m, ")") else "", sep = "")
      }, character(1)), collapse = "|"),
      match_state = paste(states, collapse = ";"),
      unlabel_ions_match.by_ions_coverage = paste(ul_cov, collapse = ";"),
      unlabel_ions_match.PSM_score_by = paste(ul_by, collapse = ";"),
      unlabel_ions_match.PSM_score_all = paste(ul_all, collapse = ";"),
      label_File_names = paste(lb_fn, collapse = ";"),
      label_ions_match.by_ions_coverage = paste(lb_cov, collapse = ";"),
      label_ions_match.PSM_score_by = paste(lb_by, collapse = ";"),
      label_ions_match.PSM_score_all = paste(lb_all, collapse = ";"),
      paired_label_cosine = paste(ifelse(is.na(dcos), "", round(dcos, 4)), collapse = ";"),
      predict_cosines = paste(ifelse(is.na(prd), "", round(prd, 4)), collapse = ";"),
      choosen_predict_cosine = if (is.na(prd[bg1])) "" else prd[bg1],
      choosen_paired_label_cosine = if (is.na(dcos[bg1])) "" else round(dcos[bg1], 4),
      best_paired_label_cosine = if (all(is.na(dcos))) "" else max(dcos, na.rm = TRUE),
      paired_state = if (isTRUE(dcos[bg1] >= 0.7)) "paired" else "single",  # 二值小写
      best_group = bg_s,
      predict_missing = if (any(is.na(prd))) 1L else 0L,
      pair_missing = if (any(is.na(dcos))) 1L else 0L,
      vote_metrics = paste(names(vm)[use], collapse = ";"),
      stringsAsFactors = FALSE)
  }
  for (j in seq_len(n)) {
    rep[[j]] <- tryCatch(row_fun(j), error = function(e) {
      #：按行隔离——单谱崩溃不再连坐整个分块（576 运行 73 行
      # 被整块丢弃的根因排查发现）；行号＋出错函数名堆栈写入日志
      tb <- paste(vapply(sys.calls(), function(cc) {
        h <- deparse(cc[[1]], nlines = 1L)
        if (length(h) == 0) "?" else substr(h[1], 1, 80)
      }, character(1)), collapse = " -> ")
      try(write(c(sprintf("ROW-ERROR raw=%s fn=%s: %s", pp$raw,
                          dat$File_Name[j], conditionMessage(e)), tb),
                file.path(out_dir, "worker_error_log.txt"), append = TRUE),
          silent = TRUE)
      NULL
    })
  }
  # 合并该 raw 全部谱的结果行并返回
  # Combine all result rows of this raw and return
  do.call(rbind, rep[!vapply(rep, is.null, logical(1))])
}
# 将工作函数挂到全局环境并显式导出所需变量，避免闭包序列化问题
# Bind the worker to the global environment and explicitly export the required variables to avoid closure serialization issues
environment(scan_worker) <- globalenv()
# 小样本空配对报告根因：B1 的 tryCatch 匿名闭包把 scan_worker 变成
# 自由符号，工作进程全局环境里没有该函数（旧代码直接传 scan_worker 作为 FUN，
# 按值序列化所以可用）。out_dir 同样被错误处理分支引用，必须一并导出，
# 否则 worker_error_log.txt 写入会静默失败（try(silent=TRUE) 吞掉）。
parallel::clusterExport(cl, c("scan_worker", "out_dir", "seq_cols", "mod_cols", "slot_syms",
                              "aa", "aa_tab", "pred_idx", "pred_int", "label_type", "labels_b"),
                        envir = environment())

# 打印成对匹配任务规模，随后并行执行所有分块工作负载
# Print the paired-matching task size, then run all chunked payloads in parallel
cat("[stage B] 成对匹配：", length(payloads), " 个分块，使用 ", use_cores, " 个进程\n")
cat("[stage B] paired matching: ", length(payloads), " chunks with ", use_cores, " workers\n")
# B1 优化tryCatch 防单块崩溃杀死全管线；parLapplyLB 动态分配
# B1 optimization tryCatch prevents a single chunk crash from
# killing the entire pipeline; parLapplyLB distributes chunks dynamically
res <- parallel::parLapplyLB(cl, payloads, function(p) {
  tryCatch(scan_worker(p), error = function(e) {
    # 错误除打印外同时写入 worker_error_log.txt（parLapply 工作进程 stdout
    # 会被丢弃，仅打印会导致错误静默——小样本空配对报告追查发现）；
    # 附函数名堆栈（snow 协议会在 worker 层吞掉内层帧，
    # 行级 tryCatch 才记录真实出错函数，这里只记分块级）
    msg <- paste0("WORKER-ERROR raw=", p$raw, ": ", conditionMessage(e))
    tb <- paste(vapply(sys.calls(), function(cc) {
      h <- deparse(cc[[1]], nlines = 1L)
      if (length(h) == 0) "?" else substr(h[1], 1, 80)
    }, character(1)), collapse = " -> ")
    cat(msg, "
")
    try(write(c(msg, tb), file.path(out_dir, "worker_error_log.txt"), append = TRUE),
        silent = TRUE)
    NULL
  })
})
res <- Filter(Negate(is.null), res)
pt_mark("成对匹配打分", "paired matching + scoring")

# 合并各 raw 的并行结果并写出配对报告 paired_report_<label_type>.tsv
# Combine the per-raw parallel results and write the paired report paired_report_<label_type>.tsv
out <- do.call(rbind, res)
if (is.null(out)) out <- data.frame()
write.table(out, file.path(out_dir, paste0("paired_report_", tolower(label_type), ".tsv")), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("paired report rows:", nrow(out), "\n")
# 有结果时打印成对状态分布与余弦统计
# When results exist, print the paired-state distribution and cosine summary
if (nrow(out) > 0) {
  cat("成对状态分布：\n")
  cat("paired_state distribution:\n")
  print(table(out$paired_state))
  cat("成对余弦分布（choosen）：\n")
  cat("paired cosine distribution (choosen):\n")
  print(summary(out$choosen_paired_label_cosine))
}
# 打印阶段 B 完成标记
# Print the stage B completion marker
cat("==== [阶段 B] 完成 ====\n")
cat("==== [Stage B] done ====\n")
