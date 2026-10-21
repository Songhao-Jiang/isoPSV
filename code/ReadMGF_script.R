library(parallel)

# =============================================================================
# ReadMGF_script.R - ReadMGF 核心库（R 包化候选：逐函数 Roxygen 双语文档，含标题/描述/参数/返回/示例）
#   ReadMGF core library (package-ready: per-function bilingual Roxygen docs with
#   title / description / params / return / examples)
# 谱图读写与解析（read_large_mgf / readMGFFile / writeMGFFile / extractParameters_MGF）
#   spectrum reading, writing and parsing
# 理论离子与打分（by_ions_position / all_ions_position / ions_match / ions_match_selected / cosine）
#   theoretical ions and scoring
# 无标记 isoPSV 流程 CLI 函数（参数文件、第一步提取、第二步 revise/合并/SAP/自动搜库、阶段 A/阶段 B）
#   unlabeled isoPSV workflow CLI helpers (parameter files, step-1 extraction,
#   step-2 revise/merge/SAP/auto-search, stage A / stage B)
# stringr / reshape2 调用由 compat.R 的 base R 实现提供（保留原函数名）；source 本文件前先 source compat.R。
#   stringr / reshape2 calls are provided by the base-R implementations in
#   compat.R (same original function names); source compat.R before this file.
# =============================================================================

#' MGF S4 类定义。
#' Definition of the S4 class MGF.
#'
#' 定义用于存储 MGF（Mascot Generic Format）文件数据的 S4 类；含验证函数，
#' 确保 charge 为正整数、rtInSeconds 为非负数。
#' Defines the S4 class holding MGF (Mascot Generic Format) file data; a
#' validity function ensures the charge is a positive integer and rtInSeconds
#' is non-negative.
#'
#' @slot title
#'   character 类型，谱图标题。
#'   Character, the spectrum title.
#' @slot charge
#'   integer 类型，电荷状态。
#'   Integer, the charge state.
#' @slot rtInSeconds
#'   numeric 类型，以秒为单位的保留时间。
#'   Numeric, the retention time in seconds.
#' @slot pepMass
#'   numeric 类型，肽段质量。
#'   Numeric, the peptide mass.
#' @slot spectrum
#'   data.frame 类型，包含 m/z-强度对的数据框。
#'   Data frame, the m/z-intensity pairs.
#' @examples
#' # 创建 MGF 对象实例
#' # create an MGF object instance
#' mgf_obj <- new("MGF",
#'                title = "Sample Spectrum",
#'                charge = 2L,
#'                rtInSeconds = 120.5,
#'                pepMass = 500.25,
#'                spectrum = data.frame(mz = c(100.1, 101.2), intensity = c(1000, 2000)))
#' @export
setClass("MGF", 
         slots = list(title = "character", charge = "integer", 
                      rtInSeconds = "numeric", pepMass = "numeric",
                      spectrum = "data.frame"), 
  prototype = list(title = character(0), charge = integer(0),
                   rtInSeconds = numeric(0), pepMass = numeric(0),
                   spectrum = data.frame()),
  validity = function(object) {
    # 检查 charge 是否为正整数
    # check that the charge is a positive integer
    if (length(object@charge) > 0 && any(object@charge <= 0)) {
      return("Charge must be a positive integer")
      }
    # 检查 rtInSeconds 是否为非负数
    # check that rtInSeconds is non-negative
    if (length(object@rtInSeconds) > 0 && any(object@rtInSeconds < 0)) {
      return("rtInSeconds must be a non-negative number")
      }
    return(TRUE)
    })
#' 块级向量化解析：一次调用解析整块 MGF 行，返回 S4 MGF 对象列表。
#' Vectorized chunk-level parsing: parse a whole block of MGF lines in one
#' call, returning a list of S4 MGF objects.
#'
#' 用 cumsum 分组、tabulate 计数、一次 strsplit 解析全部峰行，
#' 按谱区间切片组装；块须按谱边界对齐（首行为 BEGIN IONS）。
#' 与 parse_spectrum 逐谱逐行解析的输出等价，无峰行的空谱跳过。
#' Group with cumsum, count with tabulate, split ALL peak lines in one
#' strsplit call, then slice per-spectrum ranges; the chunk must be
#' spectrum-aligned (starting with BEGIN IONS). Equivalent to parsing each
#' spectrum line-by-line with parse_spectrum; empty spectra are skipped.
#'
#' @param lines
#'   字符向量，MGF 行（首行为 BEGIN IONS，块尾为完整谱块）。
#'   Character vector of MGF lines (starting at BEGIN IONS, ending on a
#'   spectrum boundary).
#' @return
#'   返回 S4 MGF 对象列表。
#'   A list of S4 MGF objects.
#' @keywords internal
parse_lines_vectorized <- function(lines) {
  if (length(lines) == 0) return(list())
  is_peak <- grepl("^[0-9.]", lines)
  is_beg  <- startsWith(lines, "BEGIN IONS")
  spec_id <- cumsum(is_beg)
  in_spec <- spec_id > 0L
  max_spectra <- sum(is_beg)

  pk    <- lines[is_peak & in_spec]
  pk_id <- spec_id[is_peak & in_spec]
  n_peaks <- length(pk)
  cnt <- tabulate(pk_id, nbins = max_spectra)
  peak_ends <- as.integer(cumsum(cnt))

  hdr_val <- function(prefix) {
    i <- which(in_spec & startsWith(lines, prefix))
    i <- i[!duplicated(spec_id[i], fromLast = TRUE)]
    v <- rep(NA_character_, max_spectra)
    v[spec_id[i]] <- sub(prefix, "", lines[i], fixed = TRUE)
    v
  }
  titles  <- hdr_val("TITLE="); titles[is.na(titles)] <- ""
  ch_raw  <- hdr_val("CHARGE=")
  charges <- rep(NA_integer_, max_spectra)
  ch_ok   <- grepl("^\\d+\\+$", ch_raw)
  charges[ch_ok] <- as.integer(sub("\\+$", "", ch_raw[ch_ok]))
  rt_raw  <- hdr_val("RTINSECONDS=")
  rts     <- ifelse(is.na(rt_raw), 0, suppressWarnings(as.numeric(rt_raw)))
  pms     <- suppressWarnings(as.numeric(sub("[ \t]+.*$", "", hdr_val("PEPMASS="))))

  if (n_peaks > 0L) {
    pairs <- strsplit(pk, "[ \t]+")
    nums <- as.numeric(unlist(pairs, use.names = FALSE))
    if (length(nums) == 2L * n_peaks) {
      mz_all  <- nums[seq.int(1L, length(nums), 2L)]
      int_all <- nums[seq.int(2L, length(nums), 2L)]
      rm(pairs, nums)
    } else {
      # 畸形峰行（非两段）：退回逐行安全解析，避免奇偶位错位。
      # Malformed peak lines (not two tokens): fall back to a safe per-line
      # parse so odd/even positions cannot shift.
      mz_all  <- vapply(pairs, function(p) as.numeric(p[1]), 0.0)
      int_all <- vapply(pairs, function(p) as.numeric(p[2]), 0.0)
    }
  } else {
    mz_all <- numeric(0)
    int_all <- numeric(0)
  }

  keep <- which(cnt > 0L)
  lo <- c(0L, peak_ends[-length(peak_ends)])[keep] + 1L
  hi <- peak_ends[keep]
  out <- vector("list", length(keep))
  for (q in seq_along(keep)) {
    seg <- seq.int(lo[q], hi[q])
    out[[q]] <- new("MGF",
                    title = titles[keep[q]],
                    charge = charges[keep[q]],
                    rtInSeconds = rts[keep[q]],
                    pepMass = pms[keep[q]],
                    spectrum = structure(
                      list(m_z = mz_all[seg], Intensity = int_all[seg]),
                      class = "data.frame", row.names = .set_row_names(length(seg))))
  }
  out
}
#' 把大 MGF 按行分块。
#' Split a big MGF into line chunks.
#'
#' 按 chunk_size 个谱图把文件的 BEGIN/END IONS 行区间切成分块；
#' 每块只携带自己的行与该块内谱图数，作为 parLapply 的 payload 发送。
#' Cuts the BEGIN/END IONS line ranges into chunks of chunk_size spectra;
#' each chunk carries only its own lines plus its spectrum count, sent as a
#' parLapply payload.
#'
#' @param filename
#'   字符串，MGF 文件的路径。
#'   Character string, the path to the MGF file.
#' @param chunk_size
#'   整数，每个块的谱图数量，默认为 500。
#'   Integer, spectra per chunk; defaults to 500.
#' @return
#'   返回 list(lines=字符向量, n=谱图数) 的列表。
#'   A list of list(lines = character vector, n = spectrum count).
#' @examples
#' if (file.exists("spectra.mgf")) {
#'   ch <- mgf_line_chunks("spectra.mgf", chunk_size = 1000)
#'   ch[[1]]$n
#' }
#' @export
mgf_line_chunks <- function(filename, chunk_size = 500) {
  all_lines <- readLines(filename, warn = FALSE)
  b <- which(all_lines == "BEGIN IONS")
  e <- which(all_lines == "END IONS")
  if (length(b) == 0) return(list())
  # 用户工作站实测尾部 BEGIN 无 END（截断/损坏块）时 e 短于
  # b，e[rng[2]] 取 NA → "NA/NaN argument" 崩第一步。截齐：丢弃无 END 的
  # 尾部谱块并双语提示。
  # live fire: a truncated tail (BEGIN without END) makes e shorter than b and
  # e[rng[2]] NA, crashing step 1. Trim unmatched trailing blocks with a note.
  if (length(e) < length(b)) {
    cat("[MGF] ", basename(filename), "：尾部 ", length(b) - length(e),
        " 个谱块缺 END IONS（截断/损坏），已跳过\n", sep = "")
    cat("[MGF] ", basename(filename), ": ", length(b) - length(e),
        " trailing block(s) lack END IONS (truncated/corrupt), skipped\n", sep = "")
    b <- b[seq_len(length(e))]
  }
  n <- length(b)
  grp <- ceiling(seq_len(n) / chunk_size)
  lapply(split(seq_len(n), grp), function(ix) {
    rng <- range(ix)
    list(lines = all_lines[b[rng[1]]:e[rng[2]]], n = length(ix))
  })
}
#' worker 端分块解析 MGF。
#' Parse an MGF chunk on the worker side.
#'
#' 把 payload 的行拆成谱图并逐个用 parse_spectrum 解析；
#' MGF S4 类在 worker 端定义，不依赖主进程环境。
#' Splits the payload lines into spectra and parses each with parse_spectrum;
#' the MGF S4 class is defined on the worker, independent of the master.
#'
#' @param payload
#'   list(lines=字符向量, n=谱图数)，mgf_line_chunks 的输出。
#'   A list(lines = character vector, n = spectrum count) from mgf_line_chunks.
#' @return
#'   返回 MGF 对象列表。
#'   A list of MGF objects.
#' @examples
#' # 先分块再解析
#' # chunk first, then parse
#' if (file.exists("spectra.mgf")) {
#'   ch <- mgf_line_chunks("spectra.mgf")
#'   objs <- parse_mgf_chunk(ch[[1]])
#' }
#' @export
parse_mgf_chunk <- function(payload) {
  if (length(payload$lines) == 0) return(list())
  parse_lines_vectorized(payload$lines)
}
#' 在共享 PSOCK 集群上加载全部无标记 MGF。
#' Load all unlabeled MGF files on a shared PSOCK cluster.
#'
#' 经 load_mgf_shared 三档策略读取（小文件每文件一个 worker，大文件分块），
#' 再展平成标题/峰表/电荷三个向量供阶段 A 哈希索引使用。
#' Reads via the three-tier load_mgf_shared (one worker per small file, chunks
#' for big files) and flattens the result into title/data/charge vectors for
#' the stage-A hash indexing.
#'
#' @param files
#'   字符向量，无标记 MGF 文件路径。
#'   Character vector of unlabeled MGF file paths.
#' @param cl
#'   PSOCK 集群对象。
#'   A PSOCK cluster object.
#' @param chunk_size
#'   整数，每个块的谱图数量，默认为 500。
#'   Integer, spectra per chunk; defaults to 500.
#' @return
#'   返回 list(titles=标题向量, data=逐谱峰表列表, charge=整数向量)。
#'   A list(titles = vector, data = per-spectrum peak data.frames, charge = integer vector).
#' @keywords internal
load_unlabeled_mgfs <- function(files, cl, chunk_size = 500) {
  loaded <- load_mgf_shared(files, cl, chunk_size)
  # v1.0 优化三向量先按块谱数预分配、分块填入，避免
  # c() 每次整向量复制（11.4s 载入里的 GC 损耗之一）
  n_spec <- sum(vapply(loaded, length, integer(1)))
  titles <- character(n_spec)
  data <- vector("list", n_spec)
  charge <- integer(n_spec)
  pos <- 1L
  for (m in loaded) {
    nm <- length(m)
    if (nm == 0L) next
    idx <- pos:(pos + nm - 1L)
    titles[idx] <- vapply(m, function(x) x@title, character(1))
    data[idx] <- lapply(m, function(x) x@spectrum)
    charge[idx] <- vapply(m, function(x) as.integer(x@charge), integer(1))
    pos <- pos + nm
  }
  list(titles = titles, data = data, charge = charge)
}
#' worker 端整读一个 MGF 文件并取其中一份。
#' Worker-side whole-file read taking one equal part.
#'
#' 整读一个 MGF 文件、取 n 个等分谱块中的第 k 块并解析；
#' 每个 worker 重读文件行，仅用于内存仍安全的中等文件（15-200MB）。
#' Reads ONE whole MGF file, takes the k-th of n equal spectral chunks and
#' parses those spectra; each worker re-reads the file lines, used only for
#' medium files (15-200MB) where this stays memory-safe.
#'
#' @param f
#'   字符串，MGF 文件路径。
#'   Character string, the MGF file path.
#' @param part
#'   整数，取第几份（1 起）。
#'   Integer, which part to take (1-based).
#' @param parts
#'   整数，总份数。
#'   Integer, total number of parts.
#' @return
#'   返回 MGF 对象列表（该份内的谱图）。
#'   A list of MGF objects (the spectra of that part).
#' @keywords internal
load_mgf_part_worker <- function(f, part, parts) {
  all_lines <- readLines(f, warn = FALSE)
  b <- which(all_lines == "BEGIN IONS")
  e <- which(all_lines == "END IONS")
  n <- length(b)
  lo <- floor((part - 1) * n / parts) + 1
  hi <- floor(part * n / parts)
  if (lo > hi) return(list())
  parse_lines_vectorized(all_lines[b[lo]:e[hi]])
}
#' 进程级缓存的 MGF 读取集群（一次建群＋一次库加载，后续调用复用）。
#' Process-level cached cluster for MGF reading (built and sourced once).
#'
#' load_unlabeled_spectra 等多处会反复建群并让每个 worker 重新 source 整个
#' 库（260KB，实测 40-60 s/次）；同一进程内第一次真正建群后全部复用，
#' 后续调用只付 parLapply 本身的成本。worker 内存代价：8 个空闲 worker
#' 常驻（约 1-1.6GB），换来每次调用省 40-60 s。
#' Several callers (load_unlabeled_spectra etc.) rebuild the cluster and
#' re-source the whole library on every worker each time (260KB, measured
#' 40-60 s per call); this helper builds once per process and every later
#' call reuses it, paying only the parLapply cost. Memory trade-off: up to
#' 8 idle workers stay resident (~1-1.6GB) in exchange for 40-60 s per call.
#'
#' @return
#'   返回 PSOCK 集群对象（同一进程内恒定）。
#'   A PSOCK cluster object (constant within one process).
#' @keywords internal
#' 探测 compat.R 所在代码目录（worker bootstrap 用）。
#' Probe the code directory holding compat.R (for worker bootstrap).
#'
#' 起用于 mgf_read_cluster；提为顶层函数，
#' 供全部集群 bootstrap 共用——worker 里 source 相对路径只当启动目录恰为
#' code/ 时成立，从其他目录启动（tests/、上级目录 Rscript code/xxx.R）
#' 会全体报 "cannot open the connection"。
#' Used by mgf_read_cluster since; promoted to top level on
#' (P0-2) for every cluster bootstrap - relative source() inside
#' workers only works when the process starts from code/; any other cwd
#' (tests/, or Rscript code/xxx.R from above) fails every worker with
#' "cannot open the connection".
probe_code_dir <- function() {
  d <- getwd()
  if (file.exists(file.path(d, "compat.R"))) return(d)
  a0 <- commandArgs(FALSE)
  f <- a0[grepl("^--file=", a0)][1]
  if (length(f) > 0 && !is.na(f)) {
    dd <- dirname(sub("^--file=", "", f))
    cands <- c(dd, file.path(d, dd),
               file.path(d, "..", "code"), file.path(dd, "..", "code"))
    for (cc in cands) {
      if (nzchar(cc) && file.exists(file.path(cc, "compat.R"))) return(cc)
    }
  }
  NULL
}
#' 统一的工作进程 bootstrap：绝对路径 source 全部核心库 + 健康检查。
#' Unified worker bootstrap: source all core libraries by absolute path plus
#' a health check.
#'
#' 旧实现三处集群裸用相对路径 source（依赖 cwd 恰为
#' code/），全部替换为本函数；末尾断言 parse_lines_vectorized 已载入，
#' 失败立即给出可读提示而非 "cannot open the connection"。
#' (P0-2): three clusters used bare relative-path source() (only
#' correct when cwd happens to be code/); all replaced by this function. A
#' final assertion checks parse_lines_vectorized is loaded, failing with a
#' readable message instead of "cannot open the connection".
ispsv_bootstrap_cluster <- function(cl) {
  cdir <- probe_code_dir()
  if (is.null(cdir)) {
    parallel::clusterEvalQ(cl, {
      source("compat.R")
      suppressMessages(source("acs_core.R"))
      suppressMessages(source("paired_params.R"))
      suppressMessages(source("ReadMGF_script.R"))
      invisible(NULL)
    })
  } else {
    ispsv_code_dir <- cdir
    parallel::clusterExport(cl, "ispsv_code_dir", envir = environment())
    parallel::clusterEvalQ(cl, {
      source(file.path(ispsv_code_dir, "compat.R"))
      suppressMessages(source(file.path(ispsv_code_dir, "acs_core.R")))
      suppressMessages(source(file.path(ispsv_code_dir, "paired_params.R")))
      suppressMessages(source(file.path(ispsv_code_dir, "ReadMGF_script.R")))
      invisible(NULL)
    })
  }
  ok <- parallel::clusterEvalQ(cl, exists("parse_lines_vectorized", inherits = TRUE))
  if (!all(vapply(ok, isTRUE, logical(1)))) {
    stop("worker bootstrap 失败：parse_lines_vectorized 未载入。请从 code/ 目录启动，或检查代码目录探测（probe_code_dir）。\n",
         "worker bootstrap failed: parse_lines_vectorized not loaded. Start from code/ or check probe_code_dir.")
  }
  invisible(cl)
}
mgf_read_cluster <- local({
  cl <- NULL
  function() {
    if (!is.null(cl)) return(cl)
    # 工作进程 bootstrap 用的代码目录（加固；提为顶层
    # probe_code_dir，见上）。在创建集群前探测 compat.R 所在目录并改用绝对路径。
    # code dir for the worker bootstrap (hardened; promoted to the
    # top-level probe_code_dir above on ). Probe for the directory
    # holding compat.R and source by absolute path instead.
    ncores <- parallel::detectCores()
    if (is.na(ncores) || ncores < 2L) ncores <- 2L
    cl <<- parallel::makeCluster(max(1L, min(8L, ncores - 1L)))
    # P3-9与其余三处集群统一走 ispsv_bootstrap_cluster
    # （绝对路径 + 健康检查），消除最后一份自有 bootstrap 分支
    ispsv_bootstrap_cluster(cl)
    cl
  }
})
#' 核心代码指纹对决定结果口径的核心脚本取 md5 串接。
#' Core code fingerprint concatenated md5 of the scripts
#' that determine the result schema.
#'
#' 断点续跑会复用旧结果文件夹；代码更新后旧结果仍被下游按新口径消费就会
#' 形成新旧混合态。总程序启动时对比本指纹（result_*/code_fingerprint.txt），
#' 不一致时醒目告警，让用户决定删除 result_* 重跑还是接受混合态。
#' Resume reuses the existing result folder; after a code update the stale
#' outputs are consumed by the new code and produce a mixed state. The master
#' program compares this fingerprint (result_*/code_fingerprint.txt) at
#' startup and warns loudly on mismatch, leaving the delete-and-rerun decision
#' to the user.
#'
#' @param code_dir
#'   字符串，代码目录，默认为当前工作目录。
#'   Character string, the code folder; defaults to the working directory.
#' @return
#'   返回指纹字符串（核心文件缺失时返回空串）。
#'   The fingerprint string ("" when the core files are missing).
#' @keywords internal
code_fingerprint <- function(code_dir = getwd()) {
  # F-5纳入 boot.R（6 个驱动脚本的唯一加载入口，决定加载顺序）、
  # master_table_export.R（结果列口径）与 stageB_paired_match.R（配对裁决口径）。
  core <- c("compat.R", "ReadMGF_script.R", "user_run.R", "user_step2.R",
            "boot.R", "master_table_export.R", "stageB_paired_match.R")
  # G-3适配器与引擎/模板目录同样决定结果口径（酶口径、修饰
  # 位点语义都在其中），纳入指纹。
  # 附带发现修复预测器适配器
  # （code/predictors/*.R，如 pred_ms2pip.R）直接决定 stage A 的 predict 列
  # 口径——此前不在指纹内，换包续跑时预测器改动不会触发混合态告警
  extra <- c(sort(list.files(file.path(code_dir, "search_adapter_*.R"), full.names = TRUE)),
             sort(list.files(file.path(code_dir, "engine_templates"), full.names = TRUE)),
             sort(list.files(file.path(code_dir, "pFind_templates"), full.names = TRUE)),
             sort(list.files(file.path(code_dir, "predictors"), pattern = "\\.R$",
                             full.names = TRUE)))
  paths <- c(file.path(code_dir, core), extra)
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return("")
  paste0(tools::md5sum(paths), collapse = "")
}
#' 在共享 PSOCK 集群上加载 MGF 文件（三档并行策略）。
#' Load MGF files on a shared PSOCK cluster (three-tier parallel strategy).
#'
#' 保留完整 S4 对象（title/charge/rtInSeconds/pepMass/spectrum），
#' 返回每输入文件一个列表。三档策略：≤15MB 每文件一个 worker；
#' 15-200MB 文件按 ceil(大小/15MB) 份拆分（上限 4），每份一个 worker 任务；
#' >200MB 主进程 readLines＋行块发送给 worker（超大文件内存安全）。
#' Keeps the FULL S4 objects (title/charge/rtInSeconds/pepMass/spectrum) and
#' returns one list per input file. Tiers: <=15MB one worker per file;
#' 15-200MB split into ceil(size/15MB) parts (cap 4), one worker task each;
#' >200MB master-side readLines with line chunks sent to the workers
#' (memory-safe for very large files).
#'
#' @param files
#'   字符向量，MGF 文件路径。
#'   Character vector of MGF file paths.
#' @param cl
#'   PSOCK 集群对象。
#'   A PSOCK cluster object.
#' @param chunk_size
#'   整数，每个块的谱图数量，默认为 500。
#'   Integer, spectra per chunk; defaults to 500.
#' @return
#'   返回列表（每输入文件一个元素，各为 MGF S4 对象列表）。
#'   A list with one element per input file, each a list of MGF S4 objects.
#' @examples
#' cl <- parallel::makeCluster(2)
#' res <- load_mgf_shared(c("a.mgf", "b.mgf"), cl)
#' parallel::stopCluster(cl)
#' @export
load_mgf_shared <- function(files, cl, chunk_size = 500) {
  sizes <- file.size(files) / (1024^2)
  parts_all <- pmin(4L, pmax(1L, ceiling(sizes / 15)))
  parts_all[sizes >= 200] <- 0L
  big <- files[sizes >= 200]
  plan <- list()
  for (k in seq_along(files)) {
    if (parts_all[k] == 0L) next
    if (parts_all[k] == 1L) {
      plan[[length(plan) + 1]] <- list(f = files[k], part = 1L, parts = 1L)
    } else {
      for (p in seq_len(parts_all[k])) {
        plan[[length(plan) + 1]] <- list(f = files[k], part = as.integer(p),
                                         parts = parts_all[k])
      }
    }
  }
  res <- vector("list", length(files))
  if (length(plan) > 0) {
    cat("### 并行读取 ", length(plan), " 个分块任务（", length(files), " 个文件）\n")
    cat("### parallel read of ", length(plan), " chunk tasks (", length(files), " files)\n")
    r <- parallel::parLapply(cl, plan, function(task)
      load_mgf_part_worker(task$f, task$part, task$parts))
    # idx 与 plan 一一对应：parts_all 为 0 的文件不产生任何 plan 条目，
    # rep() 已天然跳过——不能再按 [parts_all > 0L] 过滤（长度不等的逻辑
    # 向量在索引时会循环复用，导致其他文件的份错位丢弃）
    # idx is 1:1 with plan: files with parts_all 0 contribute no plan entries
    # and rep() already skips them - do not filter by [parts_all > 0L]
    # (a shorter logical recycles when indexing and drops the wrong parts)
    idx <- rep(seq_along(files), parts_all)
    for (k in seq_along(files)) {
      if (parts_all[k] > 0L) res[[k]] <- unlist(r[idx == k], recursive = FALSE)
    }
  }
  if (length(big) > 0) {
    # 主进程先串行读取全部大文件的行（紧凑），再对所有跨文件块做一次并行处理
    # master-side line reading for ALL big files first (serial but compact),
    # then ONE parallel pass over all chunks across files
    big_chunks <- list()
    big_file_id <- integer(0)
    for (f in big) {
      cat("### 大文件分块读取: ", f, "\n")
      cat("### big file chunked read: ", f, "\n")
      ch <- mgf_line_chunks(f, chunk_size)
      big_file_id <- c(big_file_id, rep(match(f, files), length(ch)))
      big_chunks <- c(big_chunks, ch)
    }
    r <- parallel::parLapply(cl, big_chunks, parse_mgf_chunk)
    for (k in which(parts_all == 0L)) {
      res[[k]] <- unlist(r[big_file_id == k], recursive = FALSE)
    }
    rm(r, big_chunks)
    invisible(gc())
  }
  res
}
#' 把 MGF 对象写入一个 MGF 文件。
#' Write MGF objects to an MGF file.
#'
#' 向量化写出：每个谱图一次 writeLines，字节级直通（不改变峰行内容）。
#' Vectorized output: one writeLines call per spectrum; byte-faithful pass-through.
#'
#' @param MGF_objects
#'   MGF 对象列表。
#'   A list of MGF objects.
#' @param filename
#'   字符串，输出 MGF 文件的路径。
#'   Character string, the path to the output MGF file.
#' @return
#'   无返回值（写文件）。
#'   Invisible NULL (writes the file).
#' @examples
#' mgfs <- list(parse_spectrum(c("BEGIN IONS", "TITLE=t", "PEPMASS=100", "END IONS")))
#' writeMGFFile(mgfs, tempfile(fileext = ".mgf"))
#' @export
writeMGFFile<-function (MGF_objects, filename) 
{
  # 向量化：每谱一次 writeLines，字节级直通
  # vectorized: one writeLines call per spectrum, byte-faithful pass-through
  con <- file(filename, open = "wb")
  for (MGF in MGF_objects) {
    blk <- c("BEGIN IONS",
             paste0("TITLE=", MGF@title),
             paste0("CHARGE=", MGF@charge, "+"),
             paste0("RTINSECONDS=", MGF@rtInSeconds),
             paste0("PEPMASS=", MGF@pepMass),
             paste0(MGF@spectrum$m_z, " ", MGF@spectrum$Intensity),
             "END IONS")
    writeLines(blk, con, useBytes = TRUE)
  }
  close(con)
}
#' 从 MGF 对象提取参数。
#' Extract parameters from MGF objects.
#'
#' 向量化提取各对象的标题/电荷/保留时间/肽段质量/峰数，
#' 返回可直接索引的参数表。
#' Vectorized extraction of title/charge/retention time/peptide mass/peak count;
#' returns an indexable parameter table.
#'
#' @param MGF_objects
#'   MGF 对象列表。
#'   A list of MGF objects.
#' @return
#'   返回含提取参数的数据框。
#'   A data frame containing the extracted parameters.
#' @examples
#' mgfs <- list(parse_spectrum(c("BEGIN IONS", "TITLE=t", "CHARGE=2+",
#'                              "PEPMASS=447.9", "150.04 1468.3", "END IONS")))
#' extractParameters_MGF(mgfs)
#' @export
extractParameters_MGF<-function (MGF_objects) 
{
  # 向量化提取
  # vectorized extraction
  data.frame(Title = vapply(MGF_objects, function(x) x@title, character(1)),
             Charge = vapply(MGF_objects, function(x) x@charge, integer(1)),
             RTinSeconds = vapply(MGF_objects, function(x) x@rtInSeconds, numeric(1)),
             PepMass = vapply(MGF_objects, function(x) x@pepMass, numeric(1)),
             stringsAsFactors = FALSE)
}
#' pFind qry.res 候选行的列名对应表。
#' Named field map of a pFind qry.res candidate row.
#'
#' 返回 res 候选行固定列字段名表：前 7 列固定（rank / 序列 / Calc.MH+ /
#' Raw_Score / E-value / 局部 FDR / 修饰数）；随后是 修饰数对（位点,名称）、
#' 4 个辅助列（第一列＝Specificity、末列＝蛋白命中数）、若干蛋白列
#' （"蛋白id;诱饵旗标;位置;前一残基;后一残基"，旗标 1=decoy）、其余打分明细列。
#' Returns the field-name map of the fixed leading columns: the first 7 are
#' fixed (rank / sequence / Calc.MH+ / Raw_Score / E-value / local FDR /
#' modification count); then one (site, name) pair per modification, 4
#' auxiliary columns (first = Specificity, last = protein-hit count), one
#' column per protein hit ("protein_id;decoy_flag;position;prev_AA;next_AA",
#' flag 1 = decoy) and trailing score columns.
#'
#' @return
#'   命名字符向量：列下标 → 字段名。
#'   Named character vector: column index -> field name.
#' @export
res_field_map <- function() {
  c(rank = 1L, sequence = 2L, calc_mh = 3L, raw_score = 4L,
    e_value = 5L, fdr_local = 6L, n_mods = 7L)
}
#' R 会话已用内存（MB），用 gc() 统计口径估算。
#' Estimate the R session's used memory (MB) from gc() accounting.
#'
#' R >= 4.1 的 gc() 直接给出 "used (Mb)" 列；早期实现本给出 cell 数
#' （Ncells/Vcells，按 8 字节/单元估算）。只作预算守卫用，不求精确。
#' gc() of R >= 4.1 reports a "used (Mb)" column directly; older R reports
#' cell counts (Ncells/Vcells, estimated at 8 bytes per cell). Approximate -
#' only used as a budget guard.
#'
#' @return
#'   数字，已用内存 MB 数。
#'   Numeric, used memory in MB.
#' @export
iso_mem_mb <- function() {
  g <- gc(full = FALSE)
  cn <- colnames(g)
  if ("used (Mb)" %in% cn) return(sum(g[, "used (Mb)"]))
  if ("used" %in% cn) return(sum(g[, "used"]) * 8 / 1024^2)
  NA_real_
}
#' 内存预算守卫：R 堆超过阈值即中止并打印上下文行号。
#' Memory budget guard: abort with the context line when the R heap
#' exceeds the threshold.
#'
#' 工作站 1.4M 谱全量轮曾以 75GB+ 内存跑 14 小时零输出（合并段逐行
#' rbind/cbind 拷贝反模式），本守卫把失控提前变成可定位的中止。
#' 阈值可经环境变量 ISOMEM_LIMIT_MB 覆盖；超过 70%
#' 阈值先打一次预警并触发 gc()，超过上限才中止。
#' The 1.4M-spectrum workstation run once burned 75GB+ for 14h with zero
#' output (per-row rbind/cbind copy anti-pattern); this guard turns runaway
#' memory into a locatable abort. (P1-8): the limit is overridable
#' via the ISOMEM_LIMIT_MB environment variable; crossing 70% of the limit
#' warns once and triggers gc(); crossing the limit aborts.
#'
#' @param limit_mb
#'   数字，内存上限 MB。默认按物理内存 80% 自适应（解析失败退 40000）；
#'   环境变量 ISOMEM_LIMIT_MB 覆盖。显式传参（如 40000）恒以传参为准。
#'   Numeric, memory ceiling in MB. Defaults to 80% of physical memory
#'   (fallback 40000 when probing fails); the ISOMEM_LIMIT_MB environment
#'   variable overrides. An explicit argument (e.g. 40000) always wins.
#' @param ctx
#'   字符串，中止时打印的上下文（如当前行号/阶段）。
#'   Character string, context printed on abort (e.g. current row/stage).
#' @export
iso_mem_guard <- function(limit_mb = NULL, ctx = "") {
  # N-9默认上限按物理内存 80% 自适应（小内存机器同样有保护）；
  # 显式传参（如 40000）仍按物理内存 80% 封顶——工作站 128GB 下 40GB 预算不变、
  # 16GB 机器下自动降到 ~13GB；ISOMEM_LIMIT_MB 为用户显式覆盖，不封顶。
  # 预警按 ctx 重置（换阶段可再提示一次）；预警用 message 不用 warning
  # （避免与调用方 options(warn=2)/tryCatch(warning=) 冲突）。
  env_l <- suppressWarnings(as.numeric(trimws(Sys.getenv("ISOMEM_LIMIT_MB", unset = ""))))
  if (!is.na(env_l) && env_l > 0) {
    limit_mb <- env_l
    # P3-2显式覆盖超物理内存时守卫形同虚设——探测到物理总量
    # 且超限则只警告不拦截（用户显式覆盖，后果自负），但必须让用户看见。
    tm_env <- system_total_mb()
    if (!is.na(tm_env) && tm_env > 0 && limit_mb > tm_env) {
      message(sprintf("注意：ISOMEM_LIMIT_MB=%.0f 超过物理内存 %.0f MB，内存守卫将不会生效\nnote: ISOMEM_LIMIT_MB=%.0f exceeds the %.0f MB physical RAM; the memory guard is effectively disabled",
                      limit_mb, tm_env, limit_mb, tm_env))
    }
  } else {
    if (is.null(limit_mb) || !is.finite(limit_mb) || limit_mb <= 0) limit_mb <- 40000
    tm <- system_total_mb()
    if (!is.na(tm) && tm > 0) limit_mb <- min(limit_mb, 0.8 * tm)
  }
  mb <- iso_mem_mb()
  if (is.na(mb)) return(invisible(NULL))
  if (mb > 0.7 * limit_mb) {
    last_ctx <- getOption("iso_mem_warned", "")
    if (!identical(last_ctx, ctx)) {
      options(iso_mem_warned = ctx)
      message(sprintf("内存已用 %.0f MB，达上限的 %.0f%%（%s），触发 gc() 并继续监视\nmemory at %.0f MB (%.0f%% of the limit, %s); gc() triggered, monitoring continues",
                      mb, 100 * mb / limit_mb, ctx))
      gc()
    }
  }
  if (mb > limit_mb) {
    stop(sprintf(paste0("内存预算守卫：R 堆 %.0f MB 超过上限 %.0f MB（%s）。\n",
                        "memory budget guard: R heap %.0f MB over the %.0f MB limit (%s)."),
                 mb, limit_mb, ctx))
  }
  invisible(NULL)
}
#' 酶切类型定义表。
#' Enzyme definition table .
#'
#' cut 为酶切残基集合（""＝无酶切）；nocut 为酶切残基后禁止位点（如 Trypsin 的 P，
#' 仅用于引擎配置行，不影响管线内部漏切计数——保持与历史口径一致）；sense 为
#' 酶切端：C=在 cut 残基 C 端切（肽段以 cut 残基结尾），N=在 cut 残基 N 端切
#' （肽段以 cut 残基开头）。label 为各引擎使用的标准名（pFind/Comet 等以实际
#' 模板/参数行兼容为准，见 user_step2.R 的引擎映射）。
#' cut = the cleavage-residue set ("" = no enzyme); nocut = the blocking residues
#' after a cut residue (e.g. P for trypsin) - engine config lines only, NOT used
#' by the in-pipeline missed-cleavage count (kept identical to the historical
#' K/R behaviour); sense = the side cut: C = C-terminal to the cut residue
#' (peptides end with a cut residue), N = N-terminal (peptides start with one).
#' label = the canonical name engines use.
enzyme_defs <- list(
  trypsin      = list(label = "Trypsin",       cut = "KR",    nocut = "P", sense = "C"),
  # 设计约定：Trypsin 统一为 P 切不开（不再保留无 P 例外的 Trypsin/P
  # 变体；旧配置里的 Trypsin/P 写法经别名归一为 trypsin）。LysargiNase 同 Trypsin，
  # K/R 后遇 P 不切。LysC/ArgC 同族（K/R C 端切）也带 P 阻断（与 Comet 内置表
  # 一致）。GluC 的缓冲液依赖（碳酸氢铵切 DE/磷酸缓冲只切 E）与 AspN 的 N 切点
  # （部分定义含 N：Comet 表为 DN）、Chymotrypsin 是否含 M（Comet 表为 FWYL）为
  # 口径待确认项，暂按 isoPSV 表执行（见交接文件）。
  lysarginase  = list(label = "LysargiNase",   cut = "KR",    nocut = "P", sense = "N"),
  lysc         = list(label = "LysC",          cut = "K",     nocut = "P", sense = "C"),
  lysn         = list(label = "LysN",          cut = "K",     nocut = "",  sense = "N"),
  argc         = list(label = "ArgC",          cut = "R",     nocut = "P", sense = "C"),
  gluc         = list(label = "GluC",          cut = "DE",    nocut = "P", sense = "C"),
  aspn         = list(label = "AspN",          cut = "D",     nocut = "",  sense = "N"),
  chymotrypsin = list(label = "Chymotrypsin",  cut = "FWYML", nocut = "P", sense = "C"),
  none         = list(label = "No enzyme",     cut = "",      nocut = "",  sense = "C")
)
#' 酶名归一：大小写/连字符/下划线/中文别名归一为内部键；未知值回退 trypsin 并警告。
#' Normalize an enzyme name: case, hyphens, underscores and the Chinese alias
#' normalize to the internal key; unknown values fall back to trypsin with a
#' warning.
enzyme_resolve <- function(x, default = "trypsin") {
  if (is.null(x) || length(x) == 0) return(default)
  x <- trimws(as.character(x)[1])
  if (!nzchar(x)) return(default)
  # 中文别名（v1.0 扩充）：糜蛋白酶/胰凝乳蛋白酶、胰蛋白酶、无酶切
  if (x %in% c("糜蛋白酶", "胰凝乳蛋白酶")) return("chymotrypsin")
  if (x == "胰蛋白酶") return("trypsin")
  if (x == "无酶切") return("none")
  # 注意顺序：先 tolower 再去非字母数字（否则大写字母被当作非字母数字删掉）
  key <- gsub("[^a-z0-9]", "", tolower(x))
  # G-1Trypsin/P 写法已统一为带 P 阻断的 Trypsin归一为 trypsin 时给一条 message 提示语义变化
  if (key %in% c("trypsinp", "trypsin_p")) {
    message("Trypsin/P 已统一为 Trypsin（K/R 后遇 P 不切）
Trypsin/P is normalized to Trypsin (no cleavage after K/R-P)")
  }
  map <- c(trypsin = "trypsin", trypsinp = "trypsin", trypsin_p = "trypsin",
           lysarginase = "lysarginase",
           lysargnase = "lysarginase", lysc = "lysc", lysn = "lysn",
           argc = "argc", gluc = "gluc", aspn = "aspn",
           chymotrypsin = "chymotrypsin", chymo = "chymotrypsin",
           none = "none", noenzyme = "none", unspecific = "none",
           nonspecific = "none", nospecific = "none")
  v <- unname(map[key])
  if (is.na(v) || length(v) == 0) {
    warning("未知酶切类型，按 trypsin 处理：", x,
            "\nunknown enzyme type, treating as trypsin: ", x)
    return(default)
  }
  v
}
#' 取酶定义列表（按酶名）。
#' Return the enzyme definition list for a name.
enzyme_def <- function(x = "trypsin") enzyme_defs[[enzyme_resolve(x)]]
#' 当前流程的酶切类型（ISO_ENZYME 环境变量，user_run 启动时注入；独立运行默认
#' trypsin）。按环境变量值记忆化——值不变时零解析开销，值变化时重新解析
#' （也便于测试内切换酶）。
#' The enzyme of the current run (the ISO_ENZYME environment variable injected
#' by user_run at startup; standalone runs default to trypsin). Memoized on the
#' raw env value: no re-parse cost while unchanged, re-resolved when it changes
#' (which also lets tests switch enzymes in-process).
iso_enzyme <- local({
  cached <- NULL
  function() {
    key <- Sys.getenv("ISO_ENZYME", unset = "trypsin")
    if (is.null(cached) || !identical(cached$raw, key)) {
      cached <<- list(raw = key, enz = enzyme_resolve(key))
    }
    cached$enz
  }
})
#' N 端酶切特异性判定（五类分类/Error_reason 的 N 端非特异判据参数化）。
#' prev_aa 为肽段前一残基（"_"/"-" 表示蛋白 N 端或未知），first_aa 为肽段
#' 首残基。C 端酶：前一残基必须是 cut 残基（Trypsin=K/R，与历史口径一致）；
#' N 端酶：肽段首残基必须是 cut 残基；无酶切：恒特异。蛋白 N 端（prev 为
#' "_"/"-"）视为特异起点。
#' N-terminal cleavage specificity (parameterizes the N-terminal non-specific
#' rule of the five-class classification / Error_reason). prev_aa = the residue
#' before the peptide ("_"/"-" = protein N-terminus or unknown), first_aa = the
#' peptide's first residue. C-terminal enzymes: prev must be a cut residue
#' (Trypsin = K/R, identical to the historical rule); N-terminal enzymes: the
#' first residue must be a cut residue; no enzyme: always specific. The protein
#' N-terminus (prev in "_"/"-") counts as a specific start.
enzyme_nterm_ok <- function(x = iso_enzyme(), prev_aa = "", first_aa = "") {
  d <- enzyme_def(x)
  prot_start <- c("_", "-", "\u2212")
  if (identical(d$sense, "N")) {
    return(first_aa %in% strsplit(d$cut, "", fixed = TRUE)[[1]] ||
             prev_aa %in% prot_start)
  }
  if (!nzchar(d$cut)) return(TRUE)   # 无酶切：无特异性规则
  prev_aa %in% c(strsplit(d$cut, "", fixed = TRUE)[[1]], prot_start)
}
#' 漏切位点计数（与历史 [KR] 口径同构的参数化）：C 端酶统计出现在非末位的 cut
#' 残基（Trypsin 即"非 C 端的 K/R"，与早期实现逐字节一致）；N 端酶统计出现在非首
#' 位的 cut 残基；无酶切为 0。nocut 不参与（保持历史口径）。
#' Missed-cleavage count (the parameterization of the historical [KR] rule):
#' C-terminal enzymes count cut residues NOT at the C-terminus (trypsin = K/R
#' not at the last position, byte-identical to the old version); N-terminal
#' enzymes count cut residues NOT at the first position; no enzyme = 0. The
#' nocut set is ignored (historical parity).
enzyme_missed_cleavages <- function(x = iso_enzyme(), seq = "") {
  d <- enzyme_def(x)
  if (!nzchar(d$cut) || !nzchar(seq)) return(0L)
  # G-6非正则实现——避免将来扩酶表时 cut 含特殊字符
  # （-、]、^）破坏字符类
  pos <- which(strsplit(seq, "", fixed = TRUE)[[1]] %in%
                 strsplit(d$cut, "", fixed = TRUE)[[1]])
  if (length(pos) == 0) return(0L)
  if (identical(d$sense, "N")) sum(pos != 1L) else sum(pos != nchar(seq))
}
#' 生成 pFind cfg 的 enzyme 行。
#' Build the pFind cfg enzyme line .
#'
#' pFind 行格式 enzyme=<名> <切点> <禁点> <C/N>；none 写 pFind 自带的无酶切条目
#' NoEnzyme。
#' The pFind line format is enzyme=<name> <cut> <nocut> <C/N>; none writes
#' pFind's built-in no-enzyme entry NoEnzyme (pFind3 enzyme.ini EN16, verified
#' rc=0 by the deepseek review).
pfind_enzyme_line <- function(x = "trypsin") {
  en <- enzyme_resolve(x)
  if (en == "none") return("enzyme=NoEnzyme U _ C")
  d <- enzyme_def(en)
  paste0("enzyme=", d$label, " ",
         if (nzchar(d$cut)) d$cut else "_", " ",
         if (nzchar(d$nocut)) d$nocut else "_", " ", d$sense)
}
#' 从 Comet 参数文件的 [COMET_ENZYME_INFO] 酶表取酶编号。
#' Look up a Comet enzyme number from the [COMET_ENZYME_INFO] table of a
#' params file .
#'
#' 按名称匹配（大小写不敏感）返回表中编号；未命中返回 NA（调用方自行追加）。
#' 相比早期实现硬编码编号（trypsin_p=2/lysc=3/…），本实现不依赖随包模板的编号布局，
#' 用户自改酶表或换 Comet 版本都不会错配。
#' Matches by name (case-insensitive) and returns the table's number; NA on
#' miss (the caller appends). Unlike the old hardcoded numbers, this does not
#' depend on the bundled template's numbering, so a user-edited enzyme table
#' or another Comet version cannot silently mis-map.
#' @param en
#'   内部酶键。
#'   The internal enzyme key.
#' @param enzyme_table_lines
#'   字符向量，COMET_ENZYME_INFO 块的全部行（含编号行）。
#'   A character vector of the whole COMET_ENZYME_INFO block lines.
#' @return 整数编号或 NA。
#'   The integer number or NA.
comet_enzyme_number <- function(en = "trypsin", enzyme_table_lines = character(0)) {
  label <- switch(en, trypsin = "Trypsin", lysc = "Lys_C", lysn = "Lys_N",
                  argc = "Arg_C", aspn = "Asp_N", gluc = "GluC",
                  chymotrypsin = "Chymotrypsin", none = "Cut_everywhere", "")
  if (!nzchar(label)) return(NA_integer_)
  rows <- enzyme_table_lines[grepl("^[0-9]+\\.", enzyme_table_lines)]
  for (r in rows) {
    num <- as.integer(sub("^([0-9]+)\\..*", "\\1", r))
    nm <- sub("^[0-9]+\\.\\s*([^ ]+).*", "\\1", r)
    if (!is.na(num) && tolower(nm) == tolower(label)) return(num)
  }
  NA_integer_
}
#' 追加一条 Comet 酶表行。
#' Append one Comet enzyme-table entry .
#' @return 追加行文本。
#'   The appended line text.
comet_enzyme_append_line <- function(en = "trypsin", enzyme_table_lines = character(0)) {
  d <- enzyme_def(en)
  rows <- enzyme_table_lines[grepl("^[0-9]+\\.", enzyme_table_lines)]
  nums <- suppressWarnings(as.integer(sub("^([0-9]+)\\..*", "\\1", rows)))
  num <- max(c(nums, 0L), na.rm = TRUE) + 1L
  sprintf("%d.  %-22s %d      %-10s %s", num, d$label,
          if (d$sense == "C") 1L else 0L, d$cut,
          if (nzchar(d$nocut)) d$nocut else "-")
}
#' 前体 m/zpep_mass 为质子化单同位素质量
#' （M+H），前体 m/z = (M + z·H)/z = (pep_mass + (z-1)·H)/z。旧公式
#' (pep_mass + H)/z 只在 z=2 时正确（z=1 偏高 1.0078 Da、z=3 偏高 0.336 Da）。
#' Precursor m/z pep_mass is the protonated
#' monoisotopic mass (M+H), so the precursor m/z = (M + z*H)/z =
#' (pep_mass + (z-1)*H)/z. The old formula (pep_mass + H)/z was only right at
#' z=2 (+1.0078 Da at z=1, +0.336 Da at z=3).
precursor_mz <- function(pep_mass, charge, H = 1.007825) {
  z <- as.integer(charge)
  if (is.na(z) || z < 1L) stop("precursor_mz：非法电荷 ", charge)   # B-12
  (pep_mass + (z - 1L) * H) / z
}
#' qry.res 消毒pFind 开放式搜库对修饰表
#' 外的 delta 质量写出 ";name#0" 占位修饰名，pUniFind 预处理 multiprocessing
#' 断言崩溃。处理：含 name#0 的候选整行删除；无幸存候选的谱块整块删除；
#' 幸存候选的 rank 逐块重排为 1..M 连续（rank 级 FDR 语义不变）。幂等。
#' 由 legacy run_punifind_cli 与 punifind_gpu_compat.R 的 pf_run 共用。
#' Sanitize qry.res files in place 
#' pFind's open search writes a ";name#0" placeholder for delta masses outside
#' the mod table, which trips pUniFind's preprocessing assertion. Candidate rows
#' carrying name#0 are dropped; blocks with no surviving candidate are dropped
#' whole; surviving candidates are renumbered 1..M per block (rank-level FDR
#' semantics unchanged). Idempotent. Shared by the legacy run_punifind_cli and
#' punifind_gpu_compat.R's pf_run.
pf_sanitize_qryres <- function(qry_res_dir, verbose = TRUE) {
  qry_res_files <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)
  if (!length(qry_res_files)) return(invisible(NULL))
  for (res_file in qry_res_files) {
    file_lines <- readLines(res_file, warn = FALSE)
    bad_rows <- grepl("name#0", file_lines, fixed = TRUE)
    is_block_head <- grepl("^S\t", file_lines)
    keep_rows <- rep(FALSE, length(file_lines))
    changed <- FALSE
    # block_rows：当前谱块的行下标（按 S 头切块）；surviving：块内可留的候选行
    for (block_rows in split(seq_along(file_lines), cumsum(is_block_head))) {
      if (!any(is_block_head[block_rows])) { keep_rows[block_rows] <- TRUE; next }  # 文件头杂行
      surviving_rows <- block_rows[!bad_rows[block_rows]]
      surviving_candidates <- surviving_rows[grepl("^[0-9]+\t", file_lines[surviving_rows])]
      if (length(surviving_candidates) == 0) { changed <- TRUE; next }   # 整块删除
      keep_rows[surviving_rows] <- TRUE
      old_ranks <- sub("\t.*$", "", file_lines[surviving_candidates])
      if (!identical(old_ranks, as.character(seq_along(surviving_candidates)))) {
        file_lines[surviving_candidates] <- sub("^[0-9]+\t", "", file_lines[surviving_candidates])
        file_lines[surviving_candidates] <- paste0(seq_along(surviving_candidates), "\t",
                                                   file_lines[surviving_candidates])
        changed <- TRUE
      }
    }
    if (changed) {
      writeLines(file_lines[keep_rows], res_file, useBytes = TRUE)
      if (verbose)
        cat("sanitized qry.res (name#0 rows + empty blocks dropped, ranks renumbered):",
            basename(res_file), "removed", sum(!keep_rows), "rows\n")
    }
  }
  invisible(NULL)
}

#' 快速读取 pFind qry.res 文件为结构化列表。
#' Fast reader of pFind qry.res files into a structured list.
#'
#' 每个谱块解析为一个列表元素：title（谱标题）、precursor_mass、charge 与
#' candidates 数据框（rank/sequence/calc_mh/raw_score/e_value/fdr_local/
#' n_mods/specificity/modification（"位点,名称;" 拼接）/proteins（蛋白命中
#' 表：pid/is_decoy/pos/prev_aa/next_aa 的嵌套列表列）/is_decoy（任一蛋白
#' 旗标为 1 即 TRUE））。文件头杂行自动跳过；支持一次传入多个 res 文件。
#' Each spectrum block becomes one list element: title, precursor_mass, charge
#' and a candidates data frame (rank/sequence/calc_mh/raw_score/e_value/
#' fdr_local/n_mods/specificity/modification ("site,name;" collapsed)/
#' proteins (a list column of per-hit pid/is_decoy/pos/prev_aa/next_aa)/
#' is_decoy (TRUE when any protein flag is 1)). Leading junk lines are
#' skipped; several res files may be passed at once.
#'
#' @param res_files
#'   字符向量，qry.res 文件路径（一个或多个）。
#'   Character vector of qry.res file paths (one or more).
#' @param keep_spectra
#'   可选的标题集合：仅保留这些谱的块。
#'   Optional title set: only the blocks of these spectra are kept
#'   (memory control for large-table rescue; R2 75GB fix).
#' @param drop_decoy_candidates
#'   逻辑，TRUE 时不构造诱饵候选列表项（救援只需要非诱饵候选）。
#'   Logical; when TRUE decoy candidates are not built into the list
#'   (rescue only needs non-decoy candidates).
#' @return
#'   列表：每谱一个元素（title 键）；附属性 names/charges 便捷向量由调用方
#'   自 vapply 提取。
#'   A list with one element per spectrum (keyed by title).
#' @export
pFind_res_read_one_file <- function(rf, keep_vec, drop, fmap) {
  # 单文件解析核心
  # per-file parsing core (shared by the serial and PSOCK paths; split out
  # for the R3 P0 parallel reader)
  kset <- NULL
  if (!is.null(keep_vec)) {
    kset <- new.env(parent = emptyenv(), hash = TRUE)
    for (t0 in keep_vec) kset[[t0]] <- TRUE
  }
  out <- list()
  L <- readLines(rf, warn = FALSE)
  L <- enc2utf8(L)
  n <- length(L)
  i <- 1L
    while (i <= n) {
      if (!grepl("^S\t", L[i], useBytes = TRUE)) { i <- i + 1L; next }
      sh <- strsplit(L[i], "\t", fixed = TRUE)[[1]]
      title <- if (i + 1L <= n) L[i + 1L] else ""
      # 标题不在保留集合：整块跳过（只扫下一 "S\t" 行，不解析候选）
      # title not in the keep set: skip the whole block (scan for the next
      # "S\t" line only, no candidate parsing)
      if (!is.null(kset) &&
          !exists(title, envir = kset, inherits = FALSE)) {
        j <- i + 2L
        while (j <= n && !grepl("^S\t", L[j], useBytes = TRUE)) j <- j + 1L
        i <- j
        next
      }
      # 候选行：从标题行后连续的 "数字\t" 行
      j <- i + 2L
      cand <- list()
      while (j <= n) {
        p <- strsplit(L[j], "\t", fixed = TRUE)[[1]]
        if (length(p) < 8L || !grepl("^[0-9]+$", p[1])) break
        nmods <- suppressWarnings(as.integer(p[fmap[["n_mods"]]]))
        if (is.na(nmods) || nmods < 0L) nmods <- 0L
        mod_seg <- character(0)
        invalid_delta_site <- FALSE
        off <- 8L                               # 修饰对起始（1-based：col7＝位点）
        if (nmods > 0L) {
          # 修饰对是（位点,名称）顺序
          # the modification pairs are (site, name) - verified spectrum by
          # spectrum against the official Filtered table 
          sts <- p[seq.int(off, off + 2L * nmods - 1L, 2L)]
          raw_nms <- p[seq.int(off + 1L, off + 2L * nmods - 1L, 2L)]
          delta0 <- grepl("#[0-9]+$", raw_nms)
          nms <- sub("#[0-9]+$", "", raw_nms)
          mod_seg <- paste0(sts, ",", nms, ";")
          # 设计约定：开放搜库 delta 命名（#0）把 AnyN-term /
          # AnyC-term 名分配到内部位点的候选剔除。pFind 位点惯例0＝N 端槽、len+1＝C 端槽——C 端名在 len+1、
          # N 端名在 0/1 均合法，不剔除
          # user ruling: a #0 delta-named terminal mod localized
          # at an internal site is invalid. pFind site convention (verified on
          # the full restrict F1): 0 = N-terminal slot, len+1 = C-terminal
          # slot - C-term names at len+1 and N-term names at 0/1 are legal
          if (any(delta0)) {
            seqlen <- nchar(p[2])
            sts_num <- suppressWarnings(as.integer(sts))
            for (mi in which(delta0)) {
              s <- sts_num[mi]; nm <- nms[mi]
              if (is.na(s)) { invalid_delta_site <- TRUE; next }
              if (s > seqlen + 1L) invalid_delta_site <- TRUE
              else if (grepl("AnyN-term", nm, fixed = TRUE) && s > 1L)
                invalid_delta_site <- TRUE
              else if (grepl("AnyC-term", nm, fixed = TRUE) &&
                       s >= 1L && s <= seqlen)
                invalid_delta_site <- TRUE
            }
          }
        }
        # 4 辅助列：第一列＝Specificity、末列＝蛋白命中数；蛋白列随后
        # 4 auxiliary columns: first = Specificity, last = protein-hit count;
        # protein columns follow
        specificity <- suppressWarnings(as.integer(p[off + 2L * nmods]))
        pstart <- off + 2L * nmods + 4L        # 蛋白列起始
        prots <- list()
        for (k in seq.int(pstart, length(p))) {
          seg <- strsplit(p[k], ";", fixed = TRUE)[[1]]
          if (length(seg) >= 5L && grepl("^[0-9]+$", seg[1]) &&
              seg[2] %in% c("0", "1") && grepl("^[0-9]+$", seg[3])) {
            prots[[length(prots) + 1L]] <- list(pid = seg[1],
              is_decoy = seg[2] == "1", pos = seg[3], prev_aa = seg[4],
              next_aa = seg[5])
          } else break                          # 非蛋白列 = 打分明细段开始
        }
        this_decoy <- any(vapply(prots, function(x) x$is_decoy, logical(1)))
        # drop：诱饵候选不构造（救援只用非诱饵候选）
        if (drop && this_decoy) { j <- j + 1L; next }
        cand[[length(cand) + 1L]] <- list(
          rank = as.integer(p[1]), sequence = p[2],
          calc_mh = suppressWarnings(as.numeric(p[3])),
          raw_score = suppressWarnings(as.numeric(p[4])),
          e_value = suppressWarnings(as.numeric(p[5])),
          fdr_local = suppressWarnings(as.numeric(p[6])),
          specificity = specificity,
          modification = paste(mod_seg, collapse = ""),
          proteins = prots,
          invalid_delta_site = invalid_delta_site,
          is_decoy = this_decoy)
        j <- j + 1L
      }
      out[[title]] <- list(title = title,
        precursor_mass = suppressWarnings(as.numeric(sh[2])),
        charge = suppressWarnings(as.integer(sh[3])),
        candidates = cand)
      i <- j
    }
  out
}
#' PSOCK 工作函数
#' B7worker 结果落盘 RDS、只回传文件路径——百万级
#' 嵌套列表经 PSOCK 回传时主进程单核反序列化 30+ 分钟；RDS 由 worker 本地写盘，
#' 主进程按序 readRDS，反序列化与网络回传解耦。
#' PSOCK worker function (top-level: closures are not passed around, avoiding
#' whole-frame serialization). B7 the worker saves its
#' result to an RDS file and returns only the path - million-element nested
#' lists took 30+ min of single-core deserialization on the main process when
#' returned through PSOCK; the RDS is written locally by the worker and the
#' main process reads the files back in order.
pFind_res_worker <- function(a) {
  # PSOCK 分支补逐文件进度（此前 2.5h 只打首尾两行）
  cat(sprintf("[pFind_res_reader] worker 开始：%s\n", basename(a$rf)))
  t_w <- proc.time()
  one <- pFind_res_read_one_file(a$rf, a$keep, a$drop, res_field_map())
  cat(sprintf("[pFind_res_reader] worker 完成：%s（%.1fs，%d 谱）\n",
              basename(a$rf), (proc.time() - t_w)[3], length(one)))
  saveRDS(one, a$rds, version = 2)
  a$rds
}
#' 快速读取 pFind qry.res 文件为结构化列表（串行/按文件 PSOCK 并行双路径）。
#' Fast reader of pFind qry.res files into a structured list (serial or
#' per-file PSOCK parallel).
#'
#' R3 P0：多文件时 parallel=TRUE 按文件分派 PSOCK worker（文件间
#' 互不依赖），主进程按文件序合并——后文件同标题覆盖，与串行语义逐字节一致。
#' R3 P0: with several files parallel=TRUE dispatches one PSOCK
#' worker per file (files are mutually independent); the main process merges
#' in file order, so later files overwrite duplicate titles exactly as the
#' serial path does.
pFind_res_reader <- function(res_files, keep_spectra = NULL,
                             drop_decoy_candidates = FALSE, parallel = FALSE,
                             n_workers = NULL) {
  fmap <- res_field_map()
  n_files <- length(res_files)
  if (!parallel || n_files <= 1L) {
    out <- list()
    for (rf_idx in seq_along(res_files)) {
      rf <- res_files[rf_idx]
      t_rf <- proc.time()
      cat(sprintf("[pFind_res_reader] 读取 %d/%d：%s\n",
                  rf_idx, n_files, basename(rf)))
      one <- pFind_res_read_one_file(rf, keep_spectra, drop_decoy_candidates, fmap)
      for (nm in names(one)) out[[nm]] <- one[[nm]]
      cat(sprintf("  done in %.1fs：本文件保留 %d 谱，累计 %d 谱\n",
                  (proc.time() - t_rf)[3], length(one), length(out)))
    }
    return(out)
  }
  # PSOCK 按文件并行：6 个 qry.res 串行 60-90min → 受最慢单文件约束（约 15min）
  nw <- if (is.null(n_workers)) min(n_files, max(1L, parallel::detectCores() - 1L))
        else as.integer(n_workers)
  nw <- max(1L, min(nw, n_files))
  cat(sprintf("[pFind_res_reader] PSOCK 并行读取 %d 文件（%d worker）\n", n_files, nw))
  t_par <- proc.time()
  rds_dir <- tempfile("res_rds_"); dir.create(rds_dir)
  on.exit(unlink(rds_dir, recursive = TRUE), add = TRUE)   # F-4：RDS 用后即清
  args_list <- lapply(seq_along(res_files), function(k)
    list(rf = res_files[[k]], keep = keep_spectra, drop = drop_decoy_candidates,
         rds = file.path(rds_dir, sprintf("%03d.rds", k))))
  # v1.0outfile="" 把 worker stdout 接回主进程，
  # 否则 makePSOCKcluster 默认丢弃 worker 输出，逐文件进度（pFind_res_worker 的
  # cat）在日志里不可见——与 #19 要解决的"长阻塞无进度"同形
  cl <- parallel::makePSOCKcluster(nw, outfile = "")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # 顶层函数从 .GlobalEnv 导出；args_list 是本地帧变量，单独导出
  # top-level functions come from .GlobalEnv; args_list is a local frame
  # variable and is exported separately
  parallel::clusterExport(cl, c("pFind_res_read_one_file", "pFind_res_worker",
                                "res_field_map"), envir = .GlobalEnv)
  parallel::clusterExport(cl, "args_list", envir = environment())
  rds_paths <- parallel::parLapply(cl, args_list, pFind_res_worker)
  out <- list()
  for (rp in rds_paths) {
    # deepseek §2.1worker 结果缺失必须显式报错，不得静默少并
    if (!file.exists(rp)) stop("worker 结果缺失（RDS 未落盘）：", rp,
                               "
worker result missing (RDS not written): ", rp)
    one <- readRDS(rp)
    for (nm in names(one)) out[[nm]] <- one[[nm]]
  }
  unlink(rds_dir, recursive = TRUE)
  cat(sprintf("[pFind_res_reader] 并行读取完成：%.0fs，保留 %d 谱\n",
              (proc.time() - t_par)[3], length(out)))
  out
}
#' pFind 搜库结果后处理：未过滤表替换 Filtered ＋ 诱饵行用 res 最佳非诱饵救援。
#' Post-process pFind results: the unfiltered table replaces Filtered and decoy
#' rows are rescued with the best non-decoy candidate from the res files.
#'
#' 步骤①读取 result 目录全部 *.qry.res；②以
#' pFind.spectra（未过滤全量 PSM）为基表；③Target/Decoy＝decoy 的行，从该谱
#' res 候选中取 Raw_Score 最高的非 decoy 候选替换（序列/修饰/特异性/蛋白/
#' 位置/分数；蛋白名与位置经目标库序列回贴解析，不依赖 res 的蛋白 id 映射；
#' res 无非 decoy 候选的行原样保留，交由下游诱饵去除步骤处理）；④结果写为
#' pFind-Filtered.spectra（原官方 Filtered 备份为 pFind-Filtered.origin.spectra）
#' 供结果合并、理论谱图生成与竞争打分使用。幂等：origin 备份存在时不重复备份。
#' Steps (1) read every *.qry.res in the result
#' folder; (2) start from the unfiltered pFind.spectra table; (3) for rows with
#' Target/Decoy = decoy, replace with the highest Raw_Score non-decoy candidate
#' of that spectrum (sequence/modification/specificity/proteins/position/
#' scores; protein names and positions are resolved by mapping the peptide back
#' onto the target FASTA - the res protein ids are not trusted; rows without a
#' non-decoy candidate are kept as-is for the downstream decoy-removal step);
#' (4) write the table as pFind-Filtered.spectra (the official Filtered is
#' backed up as pFind-Filtered.origin.spectra) for the merge, spectrum
#' prediction and competition scoring. Idempotent: the origin backup is only
#' written once.
#'
#' @param result_dir
#'   字符串，pFind result 目录（含 pFind.spectra 与 *.qry.res）。
#'   Character, the pFind result folder (pFind.spectra and *.qry.res inside).
#' @param fasta
#'   字符串，目标蛋白库 FASTA（非 td 版；用于救援候选的蛋白名与位置解析）。
#'   Character, the target FASTA (non-td; used to resolve protein names and
#'   positions of the rescued candidates).
#' @return
#'   命名向量 c(rows=全量行, decoy=诱饵行, rescued=已救援, kept=无救援候选保留,
#'   written=输出路径)。
#'   Named vector c(rows, decoy, rescued, kept, written).
#' @export
# v1.0pFind 后处理复用判据。指纹＝
# pFind.spectra＋全部 *.qry.res＋FASTA 的 mtime|size；与标记文件逐行一致即可
# 跳过（单槽约省 3–3.4h）。pFind.spectra 与 FASTA 只读、pf_sanitize_qryres 只在
# 内容变化时回写 qry.res——指纹自第二轮起稳定（首轮消毒后的 mtime 变化会多付
# 一次，可接受）。
pfind_postprocess_fingerprint <- function(result_dir, fasta) {
  qr <- sort(list.files(result_dir, pattern = "\\.qry\\.res$", full.names = TRUE))
  fp <- file.path(result_dir, "pFind.spectra")
  fa <- if (!is.null(fasta) && length(fasta) > 0 && nzchar(fasta[1]) &&
            file.exists(fasta[1])) fasta[1] else character(0)
  fp <- c(fp, qr, fa)
  sprintf("%s|%.0f|%.0f", normalizePath(fp, winslash = "/"),
          file.mtime(fp), file.size(fp))
}
pfind_postprocess_fresh <- function(result_dir, fasta, marker) {
  if (!file.exists(marker)) return(FALSE)
  identical(trimws(readLines(marker, warn = FALSE)),
            pfind_postprocess_fingerprint(result_dir, fasta))
}
pfind_postprocess_mark <- function(result_dir, fasta, marker) {
  writeLines(pfind_postprocess_fingerprint(result_dir, fasta), marker, useBytes = TRUE)
  invisible(TRUE)
}

pfind_res_postprocess <- function(result_dir, fasta,
                                  reader_parallel = FALSE, reader_workers = NULL) {
  sp_path <- file.path(result_dir, "pFind.spectra")
  if (!file.exists(sp_path)) stop("pFind.spectra 不存在：", sp_path, "\npFind.spectra missing: ", sp_path)
  res_files <- sort(list.files(result_dir, pattern = "\\.qry\\.res$", full.names = TRUE))
  if (length(res_files) == 0) stop("qry.res 缺失：", result_dir, "\nqry.res missing in: ", result_dir)
  sp <- read.delim(sp_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "", comment.char = "",
                   na.strings = character(0))
  # 只加载诱饵行对应谱的块、并丢弃诱饵候选
  # load only the decoy spectra's blocks and drop decoy candidates
  # (memory fix: the full ~14M-candidate list
  # once pushed RSS past 75GB)
  td_col <- grep("^Target/Decoy$", colnames(sp), value = TRUE)
  if (length(td_col) == 0) td_col <- "Target.Decoy"
  is_dec <- sp[[td_col[1]]] == "decoy"
  cat(sprintf("[pFind后处理] 全量 %d 行，诱饵 %d 行；读取对应 qry.res 块 ...\n",
              nrow(sp), sum(is_dec)))
  res <- pFind_res_reader(res_files,
                          keep_spectra = unique(sp$File_Name[is_dec]),
                          drop_decoy_candidates = TRUE,
                          parallel = reader_parallel,
                          n_workers = reader_workers)
  cat("[pFind后处理] qry.res 块读取完成\n")
  # 目标库蛋白序列（救援候选的蛋白名/位置按序列回贴解析）
  fa_names <- character(0); fa_seqs <- character(0)
  if (!is.null(fasta) && nzchar(fasta) && file.exists(fasta)) {
    ln <- readLines(fasta, warn = FALSE)
    ln <- enc2utf8(ln)
    hid <- grepl("^>", ln)
    fa_names <- sub("^>", "", ln[hid])
    fa_names <- vapply(strsplit(fa_names, " "), function(x) x[1], "")
    grp <- cumsum(hid)
    fa_seqs <- tapply(ln[!hid], grp[!hid], function(x) paste(x, collapse = ""))
    fa_seqs <- as.character(fa_seqs)
  }
  # 蛋白组拼接成单串＋各蛋白起始偏移：一条序列的全部命中由一次 C 级
  # gregexpr 得到，替代早期实现"逐救援行 vapply 扫描整蛋白组"
  # the proteome joined into one string plus per-protein start offsets: all
  # hits of one sequence come from a single C-level gregexpr, replacing the
  # old per-rescued-row whole-proteome vapply scan (the R2 14h single-core
  # stall: ~300K decoy rows x 20-30K proteins = billions of grepl calls)
  big <- ""
  cum_start <- integer(0)
  if (length(fa_seqs) > 0) {
    big <- paste(fa_seqs, collapse = "\n")
    cum_start <- cumsum(c(1L, nchar(fa_seqs) + 1L))[seq_along(fa_seqs)]
  }
  hit_env <- new.env(parent = emptyenv(), hash = TRUE)
  prot_hits <- function(seq0) {
    # 每唯一序列回贴一次（缓存命中直接复用），语义与原逐蛋白扫描一致：
    # 含该序列的蛋白按出现序取前 5，每蛋白首现位；返回 list(prots, poss)
    # 或 NULL（无命中）。空序列不回贴（pFind 行序列非空，仅防御）。
    # each unique sequence is mapped once (cache); semantics match the old
    # per-protein scan: first 5 containing proteins in order, first
    # occurrence per protein; returns list(prots, poss) or NULL.
    if (exists(seq0, envir = hit_env, inherits = FALSE))
      return(get(seq0, envir = hit_env, inherits = FALSE))
    val <- NULL
    if (length(fa_seqs) > 0 && nzchar(seq0)) {
      offs <- gregexpr(seq0, big, fixed = TRUE)[[1]]
      offs <- offs[offs > 0]
      if (length(offs) > 0) {
        pi <- findInterval(offs, cum_start)
        u <- !duplicated(pi)
        offs <- offs[u]; pi <- pi[u]
        if (length(offs) > 5) { offs <- offs[1:5]; pi <- pi[1:5] }
        pn <- character(0); pp <- character(0)
        for (h in seq_along(offs)) {
          st0 <- offs[h] - cum_start[pi[h]] + 1L
          s_len <- nchar(fa_seqs[pi[h]])
          prev_aa <- if (st0 > 1) substr(fa_seqs[pi[h]], st0 - 1L, st0 - 1L) else "-"
          next_aa <- if (st0 + nchar(seq0) <= s_len)
            substr(fa_seqs[pi[h]], st0 + nchar(seq0), st0 + nchar(seq0)) else "-"
          pn <- c(pn, fa_names[pi[h]])
          pp <- c(pp, sprintf("%d,%s,%s", st0, prev_aa, next_aa))
        }
        val <- list(prots = paste0(paste(pn, collapse = "/"), "/"),
                    poss  = paste0(paste(pp, collapse = "/"), "/"))
      }
    }
    assign(seq0, val, envir = hit_env)
    val
  }
  dec_idx <- which(is_dec)
  n_dec <- length(dec_idx)
  # 批量写入缓冲：救援值写入预分配向量——
  # 早期实现 13 个累加向量逐行 c() 追加，每次整份复制、单向量 O(n²)（实测 5 万行
  # 7.7s / 10 万行 30s，外推 45 万诱饵 ≈2.2h）；预分配后按 n_rescue 下标填写、
  # 末尾截断。res 查找改为循环外一次性 match()（逐行 match 每次重建哈希，
  # 实测比线性扫描还慢 1.8×，同样不可用）。
  # batched write buffers into PREALLOCATED vectors the old 13 accumulating vectors grew with per-row
  # c(), copying the whole vector each time - O(n²) per vector (measured 7.7s
  # at 50k and 30s at 100k rows, ~2.2h extrapolated for 452k decoys);
  # preallocation writes by a rescue counter and trims at the end. The res
  # lookup becomes ONE batch match() outside the loop (per-row match() rebuilds
  # the hash every call - measured 1.8x SLOWER than the linear scan).
  r_idx <- integer(n_dec); r_seq <- character(n_dec); r_calc <- numeric(n_dec)
  r_raw <- numeric(n_dec); r_eval <- numeric(n_dec); r_fdr <- numeric(n_dec)
  r_mod <- character(n_dec); r_spec <- character(n_dec); r_shift <- numeric(n_dec)
  r_ph_idx <- integer(n_dec); r_prots <- character(n_dec); r_poss <- character(n_dec)
  r_miss <- character(n_dec)
  n_rescue <- 0L; n_keep <- 0L
  t_loop <- proc.time()
  # B6res 是百万级命名列表，按名查找线性扫描；
  # 名称向量与诱饵标题一次取出、循环外一次性 match（哈希，860×），循环内按下标取
  res_names <- names(res)
  dec_fns <- as.character(sp$File_Name[dec_idx])
  dec_j <- match(dec_fns, res_names)
  for (pos in seq_along(dec_idx)) {
    i <- dec_idx[pos]
    j <- dec_j[pos]
    blk <- if (is.na(j)) NULL else res[[j]]
    if (is.null(blk)) { n_keep <- n_keep + 1L; next }
    best <- NULL
    for (cd in blk$candidates) {
      if (isTRUE(cd$is_decoy)) next
      # 末端限定名出现在非末端位点的 #0 delta 候选剔除
      # skip #0-delta candidates whose terminal-named mods sit at internal sites
      # 
      if (isTRUE(cd$invalid_delta_site)) next
      if (is.null(best) || (cd$raw_score > best$raw_score)) best <- cd
    }
    if (is.null(best)) { n_keep <- n_keep + 1L; next }
    # 蛋白名/位置：目标库序列回贴（每唯一序列一次；首现位置；prev/next 取上下残基）
    ph <- prot_hits(best$sequence)
    n_rescue <- n_rescue + 1L
    r_idx[n_rescue] <- i
    r_seq[n_rescue] <- best$sequence
    r_calc[n_rescue] <- best$calc_mh
    r_raw[n_rescue] <- best$raw_score
    r_eval[n_rescue] <- best$e_value
    r_fdr[n_rescue] <- if (is.finite(best$fdr_local) && best$fdr_local > 0)
      best$fdr_local else best$e_value
    r_mod[n_rescue] <- best$modification
    r_spec[n_rescue] <- as.character(best$specificity)
    r_shift[n_rescue] <- suppressWarnings(as.numeric(sp[["Exp.MH+"]][i]))
    if (length(fa_seqs) > 0) {
      # 旧口径：fasta 存在时救援行一律写回蛋白注释——无命中也写 "/"
      # old behavior: with a fasta present the rescued row always gets the
      # protein annotation rewritten - no-hit rows get "/" too (an artifact of
      # the old nzchar("/") check; kept for byte-exact equivalence)
      r_ph_idx[n_rescue] <- i
      r_prots[n_rescue] <- if (!is.null(ph)) ph$prots else "/"
      r_poss[n_rescue] <- if (!is.null(ph)) ph$poss else "/"
    }
    r_miss[n_rescue] <- as.character(enzyme_missed_cleavages(iso_enzyme(), best$sequence))
    if (pos %% 10000L == 0L) {
      cat(sprintf("[pFind后处理] 诱饵扫描 %d/%d（救援 %d，保留 %d），已用 %.0fs，R 堆 %.0f MB\n",
                  pos, n_dec, n_rescue, n_keep,
                  (proc.time() - t_loop)[3], iso_mem_mb()))
      iso_mem_guard(40000, sprintf("pFind后处理 诱饵扫描第 %d/%d 行", pos, n_dec))
    }
  }
  # 截断到实际救援数（预分配尾部置空）
  if (n_rescue < n_dec) {
    r_idx <- r_idx[seq_len(n_rescue)]; r_seq <- r_seq[seq_len(n_rescue)]
    r_calc <- r_calc[seq_len(n_rescue)]; r_raw <- r_raw[seq_len(n_rescue)]
    r_eval <- r_eval[seq_len(n_rescue)]; r_fdr <- r_fdr[seq_len(n_rescue)]
    r_mod <- r_mod[seq_len(n_rescue)]; r_spec <- r_spec[seq_len(n_rescue)]
    r_shift <- r_shift[seq_len(n_rescue)]; r_ph_idx <- r_ph_idx[seq_len(n_rescue)]
    r_prots <- r_prots[seq_len(n_rescue)]; r_poss <- r_poss[seq_len(n_rescue)]
    r_miss <- r_miss[seq_len(n_rescue)]
  }
  # 一次性批量写回（每列一次拷贝，替代逐行整列拷贝）
  # batched write-back (one copy per column instead of one per row)
  if (n_rescue > 0L) {
    sp$Sequence[r_idx] <- r_seq
    sp[["Calc.MH+"]][r_idx] <- r_calc
    if ("Mass_Shift(Exp.-Calc.)" %in% colnames(sp)) {
      ok <- is.finite(r_shift)
      if (any(ok)) sp[["Mass_Shift(Exp.-Calc.)"]][r_idx[ok]] <-
        r_shift[ok] - r_calc[ok]
    }
    sp$Raw_Score[r_idx] <- r_raw
    if ("Final_Score" %in% colnames(sp)) sp$Final_Score[r_idx] <- r_eval
    q_col <- grep("^Q[-.]value$", colnames(sp), value = TRUE)[1]
    if (!is.na(q_col)) sp[[q_col]][r_idx] <- r_fdr
    sp$Modification[r_idx] <- r_mod
    sp$Specificity[r_idx] <- r_spec
    if (length(r_ph_idx) > 0) {
      sp$Proteins[r_ph_idx] <- r_prots
      sp$Positions[r_ph_idx] <- r_poss
    }
    sp[[td_col[1]]][r_idx] <- "target"
    mc_col <- grep("^Miss[.]Clv[.]Sites$", colnames(sp), value = TRUE)[1]
    if (!is.na(mc_col)) sp[[mc_col]][r_idx] <- r_miss
  }
  rm(res, hit_env, big, cum_start)
  invisible(gc())
  cat(sprintf("[pFind后处理] 救援扫描完成：%.0fs，R 堆 %.0f MB\n",
              (proc.time() - t_loop)[3], iso_mem_mb()))
  out_path <- file.path(result_dir, "pFind-Filtered.spectra")
  bak <- file.path(result_dir, "pFind-Filtered.origin.spectra")
  if (file.exists(out_path) && !file.exists(bak)) file.rename(out_path, bak)
  write.table(sp, out_path, sep = "\t", row.names = FALSE, quote = FALSE,
              fileEncoding = "UTF-8", na = "")
  cat(sprintf("[pFind后处理] %s：全量 %d 行（诱饵 %d，救援 %d，保留 %d）-> %s\n",
              basename(result_dir), nrow(sp), sum(is_dec), n_rescue, n_keep, out_path))
  cat(sprintf("[pFind postprocess] %s: %d rows (decoy %d, rescued %d, kept %d) -> %s\n",
              basename(result_dir), nrow(sp), sum(is_dec), n_rescue, n_keep, out_path))
  c(rows = nrow(sp), decoy = sum(is_dec), rescued = n_rescue,
    kept = n_keep, written = out_path)
}
#' 计算肽段 b、y 离子的理论位置与质量。
#' Calculate theoretical positions and masses of b and y ions for a peptide.
#'
#' 给定肽段序列（考虑修饰）计算 b、y 离子的理论位置与质量；允许对特定
#' 氨基酸或末端加标记质量。
#' Computes the theoretical positions and masses of the b and y ions for a
#' given peptide sequence with modifications; supports residue/terminal labels
#' (Unlabel / N15 / C13 / custom LetterMass forms).
#'
#' @param peptide
#'   字符串，肽段氨基酸序列。
#'   Character string, the amino-acid sequence of the peptide.
#' @param modification
#'   字符串，修饰；以分号分隔，格式 "position,modification_name"，
#'   如 "1,Oxidation[M];2,Carbamidomethyl[C];"。
#'   Character string of modifications; semicolon-separated, formatted as
#'   "position,modification_name", e.g. "1,Oxidation[M];2,Carbamidomethyl[C];".
#' @param label
#'   字符串，标记类型（见详情），默认为 "Unlabel"。
#'   Character string, the label type (see Details); defaults to "Unlabel".
#' @param AA_mass_dataframe
#'   氨基酸与修饰质量数据框（AA_mass_dataframe.txt）。
#'   The amino-acid/modification mass data frame (AA_mass_dataframe.txt).
#' @param amino_acid_table
#'   氨基酸对照表（三字母/单字母/全名）。
#'   The amino-acid lookup table (three-letter/one-letter/full name).
#' @return
#'   返回 list(ions_table=数据框, pep_mass=理论单同位素质量)。
#'   A list(ions_table = data frame, pep_mass = theoretical monoisotopic mass).
#' @examples
#' # 计算一条肽段的无标记 b/y 离子
#' # b/y ions of one peptide, unlabeled
#' res <- by_ions_position("PEPTIDE", "", "Unlabel")
#' res$pep_mass
#' @export
by_ions_position<-function (peptide = "", modification = "", label = "Unlabel",
                            AA_mass_dataframe, amino_acid_table ) {
  if (nchar(peptide) == 0) {
    return("Not a peptide.")
  }
  terminal_label <- 0
  if (label != "Unlabel") {
    if (label == "N15") {
      AA_mass_dataframe$Mono_weight <- AA_mass_dataframe$N15_weight
    }
    else {
      if (label == "C13") {
        AA_mass_dataframe$Mono_weight <- AA_mass_dataframe$C13_weight
      }
      else {
        matches <- str_match_all(label, "([A-Za-z_]+)([0-9]+(?:\\.[0-9]+)?)")[[1]]
        label_df <- data.frame(Letter = matches[, 2], 
                               Number = as.numeric(matches[, 3]), stringsAsFactors = FALSE)
        for (row_index in seq_len(nrow(label_df))) {
          w <- which(AA_mass_dataframe$AA == label_df$Letter[row_index])
          if (length(w) > 0) {
            AA_mass_dataframe$Mono_weight[w] <- AA_mass_dataframe$Mono_weight[w] +
              label_df$Number[row_index]
          }
        }
        if (any(str_detect(label_df$Letter, "terminal"))) {
          terminal_label <- 1
        }
      }
    }
  }
  #：自定义修饰行（label_type "|" 行）补进质量表——修饰串里
  # Name[STY]/Name[AnyN-term] 等带标签名经此可查到组成式算出的 Δ质量
  AA_mass_dataframe <- augment_aa_mass_table(AA_mass_dataframe)
  H_atom_weight <- 1.007825
  O_atom_weight <- 15.994914
  N_terminal_weight <- H_atom_weight
  C_terminal_weight <- H_atom_weight + O_atom_weight
  if (terminal_label == 1) {
    if ("N_terminal" %in% label_df$Letter) {
      N_terminal_weight <- N_terminal_weight + label_df$Number[which(label_df$Letter == 
                                                                       "N_terminal")]
    }
    if ("C_terminal" %in% label_df$Letter) {
      C_terminal_weight <- C_terminal_weight + label_df$Number[which(label_df$Letter == 
                                                                       "C_terminal")]
    }
  }
  AA_weight <- AA_mass_dataframe$Mono_weight[match(strsplit(peptide, "", fixed = TRUE)[[1]],
                                                    AA_mass_dataframe$AA)]
  if (modification != "") {
    modification_split <- str_split(modification, ";")[[1]]
    modification_split <- modification_split[which(modification_split != 
                                                     "")]
    for (modification_index in 1:length(modification_split)) {
      temp_modification_split <- str_split(modification_split[modification_index], 
                                           ",")[[1]]
      if(!temp_modification_split[2]%in%AA_mass_dataframe$AA &&
         !str_detect(temp_modification_split[2],"->")){
        # 残基变体修饰名兜底如 Carbamidomethyl[E]
        # 表中只有 Carbamidomethyl[C] 等同族行——按基础名前缀取非末端变体行
        # 的同族质量（避免映射到 AnyN-term 行被误当末端修饰）
        # residue-variant mod-name fallback e.g.
        # Carbamidomethyl[E] has only Carbamidomethyl[C]-style family rows in
        # the table - take a NON-terminal variant of the same base name (so a
        # site-localised mod never maps to an AnyN-term row by accident)
        base_nm <- strsplit(temp_modification_split[2], "[", fixed = TRUE)[[1]][1]
        fam_i <- which(startsWith(AA_mass_dataframe$AA, base_nm) &
                       nzchar(AA_mass_dataframe$AA))
        nt_i <- fam_i[!grepl("term", AA_mass_dataframe$AA[fam_i])]
        if (length(nt_i) > 0) fam_i <- nt_i
        if (length(fam_i) > 0) temp_modification_split[2] <- AA_mass_dataframe$AA[fam_i[1]]
      }
      if(!temp_modification_split[2]%in%AA_mass_dataframe$AA){
        if(str_detect(temp_modification_split[2],"->")){
          AA_to_AA<-str_split(temp_modification_split[2],"->")[[1]]
          AA_changed<-str_split(AA_to_AA[2],"\\[")[[1]][1]
          AA_changed<-amino_acid_table$One_Letter[which(amino_acid_table$Three_Letter == AA_changed)]
          if (length(AA_changed) != 1) {
            # 不规则替换目标（如 Trp->Kynurenin、Cys->Dha）
            # irregular substitution targets (e.g. Trp->Kynurenin, Cys->Dha)
            cat("Not found modification",temp_modification_split[2],sep = "\n")
            return(NA)
          }
          mp <- suppressWarnings(as.integer(temp_modification_split[1]))
          if (is.na(mp) || mp < 1L || mp > nchar(peptide)) {   # P0-4位点越界守卫，防 NA 静默污染
            cat("修饰位点越界：", temp_modification_split[1], " / 肽长 ", nchar(peptide), "（", peptide, "）
", sep = "")
            return(NA)
          }
          AA_weight[as.numeric(temp_modification_split[1])]<-AA_mass_dataframe$Mono_weight[which(AA_mass_dataframe$AA == AA_changed)]
        }else{
          cat("Not found modification",temp_modification_split[2],sep = "\n")
          return(NA)
        }
      }else{
        mass_shift<-AA_mass_dataframe$Mono_weight[which(AA_mass_dataframe$AA==temp_modification_split[2])]
        if(str_detect(temp_modification_split[2],"N-term")){
          N_terminal_weight <- N_terminal_weight + mass_shift
        }else{
          if(str_detect(temp_modification_split[2],"C-term")){
            C_terminal_weight <- C_terminal_weight + mass_shift
          }else{
          mp <- suppressWarnings(as.integer(temp_modification_split[1]))
          if (is.na(mp) || mp < 1L || mp > nchar(peptide)) {   # P0-4位点越界守卫，防 NA 静默污染
            cat("修饰位点越界：", temp_modification_split[1], " / 肽长 ", nchar(peptide), "（", peptide, "）
", sep = "")
            return(NA)
          }
            AA_weight[as.numeric(temp_modification_split[1])]<- 
              AA_weight[as.numeric(temp_modification_split[1])] + mass_shift
          }
        }
      }
    }
  }
  ions_table <- as.data.frame(matrix(NA, nrow = max(0L, nchar(peptide) - 
                                                       1L), ncol = 5))
  colnames(ions_table) <- c("Number", "b_1_charge", "b_2_charge", 
                            "y_1_charge", "y_2_charge")
  ions_table$Number <- seq_len(max(0L, nchar(peptide) - 1L))
  ions_table$b_1_charge <- cumsum(AA_weight)[seq_len(nchar(peptide) - 1L)] +
    N_terminal_weight
  # y[x] = 倒数 x 个残基之和：对倒序残基做累加即可（无需二次倒序）
  # y[x] = the sum of the last x residues: cumsum over the reversed residues
  ions_table$y_1_charge <- cumsum(rev(AA_weight))[seq_len(nchar(peptide) - 1L)] +
    C_terminal_weight
  ions_table$y_1_charge <- ions_table$y_1_charge + H_atom_weight * 2
  ions_table$b_2_charge <- ions_table$b_1_charge/2 + H_atom_weight/2
  ions_table$y_2_charge <- ions_table$y_1_charge/2 + H_atom_weight/2
  pep_mass <- sum(AA_weight, N_terminal_weight, C_terminal_weight)
  pep_mass <- pep_mass + H_atom_weight
  return(list(ions_table = ions_table, pep_mass = pep_mass))
}
#' 计算肽段 [M]、a、b、y、-NH3、-H2O 与内离子的理论位置与质量。
#' Calculate theoretical positions and masses of [M], a, b, y, -NH3, -H2O and internal ions.
#'
#' 在 b/y 之外补充 a 离子、丢失离子与内离子表（internal_aa_table），
#' 供 ions_match_selected 做全离子打分。
#' Adds the a ions, loss ions and the internal-ion table (internal_aa_table) on
#' top of b/y, for the all-ion scoring in ions_match_selected.
#'
#' @param peptide
#'   字符串，肽段氨基酸序列。
#'   Character string, the amino-acid sequence of the peptide.
#' @param modification
#'   字符串，修饰（格式同 by_ions_position）。
#'   Character string of modifications (same format as by_ions_position).
#' @param label
#'   字符串，标记类型，默认为 "Unlabel"。
#'   Character string, the label type; defaults to "Unlabel".
#' @param AA_mass_dataframe
#'   氨基酸与修饰质量数据框。
#'   The amino-acid/modification mass data frame.
#' @param amino_acid_table
#'   氨基酸对照表。
#'   The amino-acid lookup table.
#' @return
#'   返回 list(ions_table, internal_aa_table, pep_mass)。
#'   A list(ions_table, internal_aa_table, pep_mass).
#' @examples
#' res <- all_ions_position("PEPTIDE", "", "Unlabel")
#' head(res$ions_table)
#' @export
# =============================================================================
# 候选命中的观测峰集合
# 特异子离子比例必须按"观测峰 m/z"判断，不能按肽段坐标位点判断：两条肽段的
# 同一碎片（如 L/I 等重或序列平移）落在不同位点，但谱图里命中的是同一个峰。
# 本函数返回该候选的理论 b/y 离子（含 +1/+2 电荷）在观测谱中命中的峰下标
# 集合（20 ppm 窗口）；修饰不被质量表识别时按无修饰重试一次。
# The set of observed peaks a candidate matches (m/z basis; user revision of
# ). Specific-ion ratios must be judged by the observed peak m/z,
# not by the position coordinate on the peptide: the same fragment of two
# peptides (e.g., an L/I swap or a sequence shift) sits at different positions
# yet hits the SAME peak in the spectrum. Returns the indices of the observed
# peaks hit by the candidate's theoretical b/y ions (+1/+2 charges) within a
# 20 ppm window; retries without the modification when the mass table rejects
# it.
matched_peak_set <- function(seq, mod, sp, aa, aa_tab, tol_ppm = 20,
                             mz_sorted = NULL, mz_ord = NULL) {
  for (mod_variant in unique(c(mod, ""))) {
    ion_positions <- NULL
    position_call <- try(suppressWarnings(invisible(capture.output(
      ion_positions <- cached_by_ions_position(seq, mod_variant, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(position_call, "try-error") || is.null(ion_positions) ||
        length(ion_positions) == 1 || is.null(ion_positions$ions_table)) next
    ion_table <- ion_positions$ions_table
    observed_mz <- as.numeric(sp[[1]])
    theoretical_mz <- c(ion_table$b_1_charge, ion_table$b_2_charge,
                        ion_table$y_1_charge, ion_table$y_2_charge)
    if (length(observed_mz) == 0 || length(theoretical_mz) == 0) return(integer(0))
    # O(n×m) 逐离子线性扫描改为 findInterval
    # 二分范围查询——观测峰排序后每个理论窗口 O(log n)，整体 O((n+m) log n)；
    # 返回值为观测峰在原始 sp[[1]] 中的下标（与原实现语义一致）。
    # （性能第二轮）：同谱多候选调用时排序可复用——调用方一次
    # order() 后经 mz_sorted/mz_ord 传入（映射回原始下标，语义不变）。
    # observed_mz_order：排序位置 → 原始下标映射；window_lower/window_upper：
    # 各理论离子的匹配窗口边界；left_index/right_index：窗口在排序峰表上的
    # [左闭, 右闭] 区间；matched_indices：命中的原始峰下标集合。
    if (!is.null(mz_sorted) && !is.null(mz_ord)) {
      observed_mz_order <- mz_ord
      observed_mz_sorted <- mz_sorted
    } else {
      observed_mz_order <- order(observed_mz, na.last = NA)
      observed_mz_sorted <- observed_mz[observed_mz_order]
    }
    window_lower <- theoretical_mz * (1 - tol_ppm * 1e-6)
    window_upper <- theoretical_mz * (1 + tol_ppm * 1e-6)
    # GLM v1.0 P1findInterval 对等于下界的峰返回该峰下标，
    # +1 后正好把它跳过——旧实现 which(mz >= lower) 含等值峰。left.open=TRUE
    # 返回严格小于 lower 的计数，+1 即第一个 >= lower 的位置（含等值峰）。
    left_index <- findInterval(window_lower, observed_mz_sorted, left.open = TRUE) + 1L
    right_index <- findInterval(window_upper, observed_mz_sorted)
    # P0-D理论离子 m/z 为 NA 时 findInterval 返回
    # NA，if(NA<=NA) 抛 missing value——旧实现 which(NA 比较) 静默丢弃。这里
    # 逐窗口跳过 NA 边界，与旧行为逐字节同义。
    matched_indices <- unlist(lapply(seq_along(theoretical_mz), function(k) {
      if (!is.na(left_index[k]) && !is.na(right_index[k]) &&
          left_index[k] <= right_index[k])
        observed_mz_order[left_index[k]:right_index[k]]
      else integer(0)
    }))
    return(sort(unique(matched_indices)))
  }
  integer(0)
}
all_ions_position<-function (peptide = "", modification = "", label = "Unlabel",
                            AA_mass_dataframe, amino_acid_table ) {
  if (nchar(peptide) == 0) {
    return("Not a peptide.")
  }
  terminal_label <- 0
  if (label != "Unlabel") {
    if (label == "N15") {
      AA_mass_dataframe$Mono_weight <- AA_mass_dataframe$N15_weight
    }
    else {
      if (label == "C13") {
        AA_mass_dataframe$Mono_weight <- AA_mass_dataframe$C13_weight
      }
      else {
        matches <- str_match_all(label, "([A-Za-z_]+)([0-9]+(?:\\.[0-9]+)?)")[[1]]
        label_df <- data.frame(Letter = matches[, 2], 
                               Number = as.numeric(matches[, 3]), stringsAsFactors = FALSE)
        for (row_index in seq_len(nrow(label_df))) {
          w <- which(AA_mass_dataframe$AA == label_df$Letter[row_index])
          if (length(w) > 0) {
            AA_mass_dataframe$Mono_weight[w] <- AA_mass_dataframe$Mono_weight[w] +
              label_df$Number[row_index]
          }
        }
        if (any(str_detect(label_df$Letter, "terminal"))) {
          terminal_label <- 1
        }
      }
    }
  }
  #：自定义修饰行（label_type "|" 行）补进质量表——修饰串里
  # Name[STY]/Name[AnyN-term] 等带标签名经此可查到组成式算出的 Δ质量
  AA_mass_dataframe <- augment_aa_mass_table(AA_mass_dataframe)
  H_atom_weight <- 1.007825
  O_atom_weight <- 15.994914
  # 计算理论质量损失
  # calculate theoretical mass loss
  H2O_loss <- 18.010565  # 水分子的质量
  NH3_loss <- 17.026549  # 氨分子的质量
  CO_loss <- 27.994915   # CO分子的质量（用于计算a离子）
  if(label == "N15"){NH3_loss<-NH3_loss+15.000109-14.003073}
  if(label == "C13"){CO_loss<-CO_loss+13.003355-12.000000}
  ion_length<-(nchar(peptide) - 1)
  N_terminal_weight <- H_atom_weight
  C_terminal_weight <- H_atom_weight + O_atom_weight
  if (terminal_label == 1) {
    if ("N_terminal" %in% label_df$Letter) {
      N_terminal_weight <- N_terminal_weight + label_df$Number[which(label_df$Letter == 
                                                                       "N_terminal")]
    }
    if ("C_terminal" %in% label_df$Letter) {
      C_terminal_weight <- C_terminal_weight + label_df$Number[which(label_df$Letter == 
                                                                       "C_terminal")]
    }
  }
  AA_weight <- AA_mass_dataframe$Mono_weight[match(strsplit(peptide, "", fixed = TRUE)[[1]],
                                                    AA_mass_dataframe$AA)]
  if (modification != "") {
    modification_split <- str_split(modification, ";")[[1]]
    modification_split <- modification_split[which(modification_split != "")]
    for (modification_index in 1:length(modification_split)) {
      temp_modification_split <- str_split(modification_split[modification_index],",")[[1]]
      if(!temp_modification_split[2]%in%AA_mass_dataframe$AA &&
         !str_detect(temp_modification_split[2],"->")){
        # 残基变体修饰名兜底如 Carbamidomethyl[E]
        # 表中只有 Carbamidomethyl[C] 等同族行——按基础名前缀取非末端变体行
        # 的同族质量（避免映射到 AnyN-term 行被误当末端修饰）
        # residue-variant mod-name fallback e.g.
        # Carbamidomethyl[E] has only Carbamidomethyl[C]-style family rows in
        # the table - take a NON-terminal variant of the same base name (so a
        # site-localised mod never maps to an AnyN-term row by accident)
        base_nm <- strsplit(temp_modification_split[2], "[", fixed = TRUE)[[1]][1]
        fam_i <- which(startsWith(AA_mass_dataframe$AA, base_nm) &
                       nzchar(AA_mass_dataframe$AA))
        nt_i <- fam_i[!grepl("term", AA_mass_dataframe$AA[fam_i])]
        if (length(nt_i) > 0) fam_i <- nt_i
        if (length(fam_i) > 0) temp_modification_split[2] <- AA_mass_dataframe$AA[fam_i[1]]
      }
      if(!temp_modification_split[2]%in%AA_mass_dataframe$AA){
        if(str_detect(temp_modification_split[2],"->")){
          AA_to_AA<-str_split(temp_modification_split[2],"->")[[1]]
          AA_changed<-str_split(AA_to_AA[2],"\\[")[[1]][1]
          AA_changed<-amino_acid_table$One_Letter[which(amino_acid_table$Three_Letter == AA_changed)]
          if (length(AA_changed) != 1) {
            # 不规则替换目标（如 Trp->Kynurenin、Cys->Dha）
            # irregular substitution targets (e.g. Trp->Kynurenin, Cys->Dha)
            cat("Not found modification",temp_modification_split[2],sep = "\n")
            return(NA)
          }
          mp <- suppressWarnings(as.integer(temp_modification_split[1]))
          if (is.na(mp) || mp < 1L || mp > nchar(peptide)) {   # P0-4位点越界守卫，防 NA 静默污染
            cat("修饰位点越界：", temp_modification_split[1], " / 肽长 ", nchar(peptide), "（", peptide, "）
", sep = "")
            return(NA)
          }
          AA_weight[as.numeric(temp_modification_split[1])]<-AA_mass_dataframe$Mono_weight[which(AA_mass_dataframe$AA == AA_changed)]
        }else{
          cat("Not found modification",temp_modification_split[2],sep = "\n")
          return(NA)
        }
      }else{
        mass_shift<-AA_mass_dataframe$Mono_weight[which(AA_mass_dataframe$AA==temp_modification_split[2])]
        if(str_detect(temp_modification_split[2],"N-term")){
          N_terminal_weight <- N_terminal_weight + mass_shift
        }else{
          if(str_detect(temp_modification_split[2],"C-term")){
            C_terminal_weight <- C_terminal_weight + mass_shift
          }else{
          mp <- suppressWarnings(as.integer(temp_modification_split[1]))
          if (is.na(mp) || mp < 1L || mp > nchar(peptide)) {   # P0-4位点越界守卫，防 NA 静默污染
            cat("修饰位点越界：", temp_modification_split[1], " / 肽长 ", nchar(peptide), "（", peptide, "）
", sep = "")
            return(NA)
          }
            AA_weight[as.numeric(temp_modification_split[1])]<- 
              AA_weight[as.numeric(temp_modification_split[1])] + mass_shift
          }
        }
      }
    }
  }
  
  # b 离子
  #ions b
  b_1_charge <- cumsum(AA_weight)[seq_len(ion_length)] + N_terminal_weight
  b_2_charge <- b_1_charge/2 + H_atom_weight/2
  b_1_charge_NH3 <- b_1_charge - NH3_loss
  b_2_charge_NH3 <- b_2_charge - NH3_loss/2
  b_1_charge_H2O <- b_1_charge - H2O_loss
  b_2_charge_H2O <- b_2_charge - H2O_loss/2
  
  # a 离子
  #ions a
  a_1_charge <- b_1_charge - CO_loss
  a_2_charge <- b_2_charge - CO_loss/2
  a_1_charge_NH3 <- a_1_charge - NH3_loss
  a_2_charge_NH3 <- a_2_charge - NH3_loss/2
  a_1_charge_H2O <- a_1_charge - H2O_loss
  a_2_charge_H2O <- a_2_charge - H2O_loss/2
  
  # y 离子
  #ions y
  y_1_charge <- cumsum(rev(AA_weight))[seq_len(ion_length)] +
    C_terminal_weight + H_atom_weight * 2
  y_2_charge <- y_1_charge/2 + H_atom_weight/2
  y_1_charge_NH3 <- y_1_charge - NH3_loss
  y_2_charge_NH3 <- y_2_charge - NH3_loss/2
  y_1_charge_H2O <- y_1_charge - H2O_loss
  y_2_charge_H2O <- y_2_charge - H2O_loss/2
  
  pep_mass <- sum(AA_weight, N_terminal_weight, C_terminal_weight) + H_atom_weight
  M_2_charge <- pep_mass/2 + H_atom_weight/2
  M_2_charge_NH3 <- M_2_charge - NH3_loss/2
  M_2_charge_H2O <- M_2_charge - H2O_loss/2
  
  # 离子表
  ##ions_table
  ions_table <- as.data.frame(matrix(NA,nrow = ion_length*18+3,ncol = 3))
  colnames(ions_table )<- c("Class","Charge","m_z")
  ions_table$Class<-c(rep("b",2*ion_length),rep("b_NH3",2*ion_length),rep("b_H2O",2*ion_length),
                      rep("y",2*ion_length),rep("y_NH3",2*ion_length),rep("y_H2O",2*ion_length),
                      rep("a",2*ion_length),rep("a_NH3",2*ion_length),rep("a_H2O",2*ion_length),
                      "M","M_NH3","M_H2O")
  ions_table$Charge <- c(rep(c(rep(1,ion_length),rep(2,ion_length)),9),rep(2,3))
  ions_table$m_z<-c(b_1_charge,b_2_charge,b_1_charge_NH3,b_2_charge_NH3,b_1_charge_H2O,b_2_charge_H2O,
                    y_1_charge,y_2_charge,y_1_charge_NH3,y_2_charge_NH3,y_1_charge_H2O,y_2_charge_H2O,
                    a_1_charge,a_2_charge,a_1_charge_NH3,a_2_charge_NH3,a_1_charge_H2O,a_2_charge_H2O,
                    M_2_charge,M_2_charge_NH3,M_2_charge_H2O)
  
  # 内离子氨基酸表
  ##internal aa
  #：短肽（≤4 残基）守卫——两个坑：① 1:(n-2) 在 n<2 时变成倒序
  # 区间（如 2:1 取 2 个元素），paste 结果赋进 0 行 data.frame 报
  # "replacement has N rows, data has 0"；② R 4.x 对 0 行 data.frame 做
  # $col <- 标量 赋值同样报该错（denovo 短肽候选触发，整簇 worker 崩溃）。
  # 统一走 mk_int_tab：0 行直接返回带列空表，不经过任何赋值。
  #: short-peptide (<=4 residues) guards - two traps: (a) 1:(n-2)
  # becomes a REVERSED index run when n < 2 (2:1 yields 2 elements) whose
  # paste result is assigned into a 0-row frame; (b) R 4.x errors the same
  # way on $col <- scalar into a 0-row frame. Triggered by short de-novo
  # candidates, killing the whole worker cluster. mk_int_tab returns a
  # fully-typed empty frame for 0 rows without any assignment.
  mk_int_tab <- function(cls, mz) {
    if (length(cls) == 0L)
      return(data.frame(Class = character(0), Charge = integer(0),
                        m_z = numeric(0), stringsAsFactors = FALSE))
    d <- data.frame(Class = as.character(cls),
                    Charge = rep(1L, length(cls)),
                    m_z = as.numeric(mz), stringsAsFactors = FALSE)
    d_NH3 <- d; d_NH3$Class <- paste0(d_NH3$Class, "_NH3")
    d_NH3$m_z <- d_NH3$m_z - NH3_loss
    d_H2O <- d; d_H2O$Class <- paste0(d_H2O$Class, "_H2O")
    d_H2O$m_z <- d_H2O$m_z - H2O_loss
    rbind(d, d_NH3, d_H2O)
  }
  internal_aa <- if (ion_length >= 2L)
    str_split(str_sub(peptide, 2, ion_length), "")[[1]] else character(0)
  AA_weight_int <- AA_weight[if (ion_length >= 2L) 2:ion_length else integer(0)]
  names(AA_weight_int) <- internal_aa

  k2 <- max(0L, ion_length - 2L)
  i2 <- seq_len(k2)
  k3 <- max(0L, ion_length - 3L)
  i3 <- seq_len(k3)
  internal_aa_table <- rbind(
    mk_int_tab(internal_aa, AA_weight_int + H_atom_weight),
    mk_int_tab(internal_aa[i2],
               AA_weight_int[i2] + AA_weight_int[i2 + 1L] + H_atom_weight),
    mk_int_tab(paste(internal_aa[i3], internal_aa[i3 + 1L],
                     internal_aa[i3 + 2L], sep = "_"),
               AA_weight_int[i3] + AA_weight_int[i3 + 1L] +
                 AA_weight_int[i3 + 2L] + H_atom_weight))

  return(list(ions_table = ions_table, internal_aa_table = internal_aa_table, pep_mass = pep_mass))
}
#' 计算两向量间余弦相似度。
#' Calculate the cosine similarity between two vectors.
#'
#' 即两非零向量夹角余弦＝内积除以模长之积；
#' 用于预测谱与实际谱、无标记谱与标记谱的相似性度量。
#' The cosine of the angle between two non-zero vectors, i.e. their inner
#' product divided by the product of their magnitudes; used as the similarity
#' between predicted and observed spectra, and between unlabeled and labeled spectra.
#'
#' @param x
#'   数值向量。
#'   A numeric vector.
#' @param y
#'   另一个数值向量；长度应与 x 相同。
#'   Another numeric vector; its length should equal that of x.
#' @return
#'   返回 x、y 余弦相似度的数值。
#'   A numeric value, the cosine similarity between x and y.
#' @examples
#' cosine(c(1, 2, 3), c(4, 5, 6))
#' @export
cosine<-function (x = c(1, 2, 3), y = c(4, 5, 6)) 
{
  if (length(x) != length(y)) {
    stop("Vectors x and y must be of the same length.")
  }
  return(sum(x * y)/(sqrt(sum(x * x)) * sqrt(sum(y * y))))
}
#' 闭区间窗口 [low, high] 内的强度求和（向量化）。
#' Windowed intensity sums over closed intervals [low, high] (vectorized).
#'
#' 对实验谱 m/z 排序一次、强度做前缀和，再用 findInterval 二分定位每个
#' 窗口，取 [lo, hi] 闭区间内的强度之和；多个窗口一次调用完成。
#' Sorts the observed m/z once, builds the intensity prefix sums and locates
#' every window with findInterval, summing the intensities inside each closed
#' [lo, hi] window; all windows are handled in one call.
#'
#' @param lo
#'   数值向量，窗口下界。
#'   Numeric vector, the window lower bounds.
#' @param hi
#'   数值向量，窗口上界（与 lo 等长）。
#'   Numeric vector, the window upper bounds (same length as lo).
#' @param mz
#'   数值向量，实验谱 m/z。
#'   Numeric vector, the observed m/z values.
#' @param int
#'   数值向量，实验谱强度（与 mz 等长）。
#'   Numeric vector, the observed intensities (same length as mz).
#' @return
#'   返回与 lo 等长的强度向量。
#'   A numeric vector of the summed intensities, one per window.
#' @keywords internal
win_intensity <- function(lo, hi, mz, int) {
  o <- order(mz)
  mz_s <- mz[o]
  ps <- c(0, cumsum(int[o]))
  n_lt <- findInterval(lo - 1e-9, mz_s)
  n_le <- findInterval(hi, mz_s)
  out <- ps[n_le + 1L] - ps[n_lt + 1L]
  out[is.na(lo) | is.na(hi)] <- 0
  out
}
#' 每个 m/z 点被多少窗口覆盖（向量化）。
#' Count how many windows cover each m/z point (vectorized).
#'
#' 覆盖数 = 满足 low <= m/z <= high 的窗口个数；窗口下界/上界各排序一次后
#' 用 findInterval 定位，逐点两次二分即可。
#' The cover count of a point is the number of windows with low <= m/z <= high;
#' the window bounds are each sorted once and every point is located with two
#' findInterval calls.
#'
#' @param px
#'   数值向量，待查询的 m/z 点。
#'   Numeric vector, the query m/z points.
#' @param lo
#'   数值向量，窗口下界。
#'   Numeric vector, the window lower bounds.
#' @param hi
#'   数值向量，窗口上界（与 lo 等长）。
#'   Numeric vector, the window upper bounds (same length as lo).
#' @return
#'   返回与 px 等长的整数覆盖计数。
#'   An integer vector of cover counts, one per query point.
#' @keywords internal
cover_count <- function(px, lo, hi) {
  ok <- !is.na(lo) & !is.na(hi)   # NA 窗口（lo/hi 配对）整体剔除
  lo <- lo[ok]
  hi <- hi[ok]
  findInterval(px, sort(lo)) - findInterval(px - 1e-9, sort(hi))
}

# 氨基酸字母
#AA letters
{
  # 20种氨基酸三字母缩写
  # 20 amino-acid three-letter abbreviations
  amino_acids_three_letter <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val","Xle")
  # 20种氨基酸单字母缩写
  # 20 amino-acid one-letter abbreviations
  amino_acids_one_letter <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","I")
  # 氨基酸全名
  # amino-acid full names
  amino_acid_names <- c("Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
                        "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
                        "Leucine","Lysine","Methionine","Phenylalanine","Proline",
                        "Serine","Threonine","Tryptophan","Tyrosine","Valine","Isoleucine_or_Leucine")
  # 创建对照表
  # create the lookup table
  amino_acid_table <- data.frame(  Full_Name = amino_acid_names,  Three_Letter = amino_acids_three_letter,
                                   One_Letter = amino_acids_one_letter)
  rm(amino_acids_three_letter)
  rm(amino_acids_one_letter)
  rm(amino_acid_names)
  
}
#' b/y 离子匹配与打分。
#' b/y ion matching and scoring.
#'
#' 把 by_ions_position 的理论 b/y 离子与实验谱匹配（ppm 或 Da 容差），
#' 计算覆盖度、PSM_score_by（b/y 离子强度占比）与 PSM_score_all（全匹配离子强度占比）；
#' 支持 N15/C13 标记态的理论质量重算。
#' Matches the theoretical b/y ions from by_ions_position against the observed
#' spectrum (ppm or Da tolerance) and computes the coverage, PSM_score_by
#' (b/y ion-intensity share) and PSM_score_all (all matched-ion share); supports
#' recomputing theoretical masses under the N15/C13 label states.
#'
#' @param therotical_by_ions
#'   by_ions_position 的输出。
#'   The output of by_ions_position.
#' @param mgf_ions
#'   数据框，实验谱的 m/z-强度表。
#'   A data frame, the observed m/z-intensity table.
#' @param mass_shift_tolerence
#'   数值，质量容差，默认为 20。
#'   Numeric, the mass tolerance; defaults to 20.
#' @param ppm_or_Da
#'   字符向量 c("ppm","Da")，容差单位，默认取 "ppm"。
#'   Character vector c("ppm","Da"), the tolerance unit; "ppm" by default.
#' @param N15_label_state
#'   逻辑值，是否按 N15 标记态重算理论质量，默认为 FALSE。
#'   Logical, whether to recompute masses under the N15 label state; FALSE.
#' @param C13_label_state
#'   逻辑值，是否按 C13 标记态重算理论质量，默认为 FALSE。
#'   Logical, whether to recompute masses under the C13 label state; FALSE.
#' @return
#'   返回含 by_ions_coverage、PSM_score_by、PSM_score_all 与
#'   by_ions_long.match_intensity 等成分的列表。
#'   A list with by_ions_coverage, PSM_score_by, PSM_score_all and
#'   by_ions_long.match_intensity among its components.
#' @examples
#' ions <- by_ions_position("PEPTIDE", "", "Unlabel")
#' spec <- data.frame(m_z = c(226.11, 227.11, 370.17, 371.17), Intensity = c(1000, 900, 800, 700))
#' ions_match(ions, spec, mass_shift_tolerence = 20, ppm_or_Da = c("ppm", "Da"))
#' @export

# ---- 进程级理论离子缓存------------------
# by_ions_position / all_ions_position 都是 (seq|mod|label) 纯函数（内部的
# AA_mass_dataframe 改写是函数内局部复制，不影响调用方），阶段 A/B 的 57k 谱
# 里 (seq|mod) 候选大量重复（SJH 实测 19k 唯一肽/57k 行）。哈希缓存跨候选、
# 跨分块复用——每个 worker 进程各持一份（库 source 时创建）。R 写时复制保证
# 调用方对返回对象后续修改（如加 match 列）不会污染缓存。容量上限防 worker
# 内存膨胀：超限整体清空（纯函数，清空只是丢加速）。
# Process-level theoretical-ion cache both
# functions are pure in (seq|mod|label) (their internal AA-mass rewrites are
# function-local copies); (seq|mod) candidates repeat heavily across the 57k
# spectra (SJH: 19k unique peptides / 57k rows). The hash cache reuses results
# across candidates and chunks - one cache per worker process, created when the
# library is sourced. R copy-on-modify keeps later caller mutations (e.g. the
# match column) from polluting the cache. A size cap keeps worker memory
# bounded: when full the cache is cleared wholesale (pure functions, so
# clearing merely drops the speedup).
theo_cache_env <- new.env(hash = TRUE, parent = emptyenv())
theo_cache_max <- 5000L
cached_by_ions_position <- function(peptide, modification, label,
                                    AA_mass_dataframe, amino_acid_table) {
  # P0-1键加 "by\r" 命名空间前缀——与
  # cached_all_ions_position 共用环境但返回结构不同（宽表 vs 长表），键格式
  # 相同会互相覆盖：all 先入缓存时 by 取到 all 对象 → 读列全空 →
  # candidate_specific_ratio/spec_ratio_* 整列塌空、混合谱类归零
  key <- paste0("by\r", peptide, "\r", modification, "\r", label)
  if (exists(key, envir = theo_cache_env, inherits = FALSE))
    return(get(key, envir = theo_cache_env, inherits = FALSE))
  val <- by_ions_position(peptide, modification, label,
                          AA_mass_dataframe, amino_acid_table)
  if (length(ls(theo_cache_env, all.names = TRUE)) >= theo_cache_max)
    rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)
  assign(key, val, envir = theo_cache_env)
  val
}
cached_all_ions_position <- function(peptide, modification, label,
                                     AA_mass_dataframe, amino_acid_table) {
  # P0-1键加 "all\r" 命名空间前缀
  key <- paste0("all\r", peptide, "\r", modification, "\r", label)
  if (exists(key, envir = theo_cache_env, inherits = FALSE))
    return(get(key, envir = theo_cache_env, inherits = FALSE))
  val <- all_ions_position(peptide, modification, label,
                           AA_mass_dataframe, amino_acid_table)
  if (length(ls(theo_cache_env, all.names = TRUE)) >= theo_cache_max)
    rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)
  assign(key, val, envir = theo_cache_env)
  val
}
ions_match<-function(therotical_by_ions = by_ions_position(...), mgf_ions = NULL,
                     mass_shift_tolerence = 20, ppm_or_Da = c("ppm","Da"), 
                     N15_label_state = FALSE, C13_label_state = FALSE){
  if (is.null(mgf_ions)) stop("mgf_ions required (the legacy read_large_mgf_auto default was removed )")
  ions_table<-therotical_by_ions$ions_table
  by_ions_long<-suppressMessages(melt(ions_table[,2:ncol(ions_table)]))
  colnames(by_ions_long)<-c("Type","m_z")
  by_ions_long$Type<-as.character(by_ions_long$Type)
  by_ions_long$m_z<-as.numeric(by_ions_long$m_z)
  by_ions_long<-cbind(by_ions_long,Charge = 0)
  by_ions_long$Charge[which(by_ions_long$Type == "b_1_charge")]<-1
  by_ions_long$Charge[which(by_ions_long$Type == "y_1_charge")]<-1
  by_ions_long$Charge[which(by_ions_long$Type == "b_2_charge")]<-2
  by_ions_long$Charge[which(by_ions_long$Type == "y_2_charge")]<-2
  by_ions_long$Type[which(by_ions_long$Type == "b_1_charge")]<-"b"
  by_ions_long$Type[which(by_ions_long$Type == "y_1_charge")]<-"y"
  by_ions_long$Type[which(by_ions_long$Type == "b_2_charge")]<-"b"
  by_ions_long$Type[which(by_ions_long$Type == "y_2_charge")]<-"y"
  
  if(C13_label_state == TRUE){
    by_ions_long_1<-by_ions_long
    by_ions_long_1$m_z<-by_ions_long_1$m_z-1/by_ions_long_1$Charge
    by_ions_long_2<-by_ions_long
    by_ions_long_2$m_z<-by_ions_long_2$m_z-2/by_ions_long_2$Charge
  }
  
  
  # 计算理论质量损失
  # calculate theoretical mass loss
  H2O_loss <- 18.010565  # 水分子的质量
  NH3_loss <- 17.026549  # 氨分子的质量
  CO_loss <- 27.994915   # CO分子的质量（用于计算a离子）
  if(N15_label_state == TRUE){NH3_loss<-NH3_loss+1}
  if(C13_label_state == TRUE){CO_loss<-CO_loss+1}
  
  ions_H2O_loss<-by_ions_long
  ions_H2O_loss$m_z<-ions_H2O_loss$m_z-H2O_loss/ions_H2O_loss$Charge
  
  ions_NH3_loss<-by_ions_long
  ions_NH3_loss$m_z<-ions_NH3_loss$m_z-NH3_loss/ions_NH3_loss$Charge
  
  a_ions<-by_ions_long[which(by_ions_long$Type == "b"),]
  a_ions$m_z<-a_ions$m_z-CO_loss/a_ions$Charge
  
  combined_ions<-rbind(by_ions_long,ions_H2O_loss,ions_NH3_loss,a_ions)
  if(C13_label_state == TRUE){
    combined_ions_1<-combined_ions
    combined_ions_1$m_z<-combined_ions_1$m_z-1/combined_ions_1$Charge
    combined_ions_2<-combined_ions
    combined_ions_2$m_z<-combined_ions_2$m_z-2/combined_ions_2$Charge
  }
  
  if(!tolower(ppm_or_Da)%in%c("ppm","da")){cat("Error mass shift tolerence!");return(NA)}
  if(!is.numeric(mass_shift_tolerence)){cat("Error mass shift tolerence!");return(NA)}
  
  if(tolower(ppm_or_Da) == "ppm"){
    by_ions_long<-cbind(by_ions_long, low_m_z = by_ions_long$m_z*(1-mass_shift_tolerence*10^-6))
    by_ions_long<-cbind(by_ions_long, high_m_z = by_ions_long$m_z*(1+mass_shift_tolerence*10^-6))
    combined_ions<-cbind(combined_ions, low_m_z = combined_ions$m_z*(1-mass_shift_tolerence*10^-6))
    combined_ions<-cbind(combined_ions, high_m_z = combined_ions$m_z*(1+mass_shift_tolerence*10^-6))
  }else{
    by_ions_long<-cbind(by_ions_long, low_m_z = by_ions_long$m_z-mass_shift_tolerence/by_ions_long$Charge)
    by_ions_long<-cbind(by_ions_long, high_m_z = by_ions_long$m_z+mass_shift_tolerence/by_ions_long$Charge)
    combined_ions<-cbind(combined_ions, low_m_z = combined_ions$m_z-mass_shift_tolerence/combined_ions$Charge)
    combined_ions<-cbind(combined_ions, high_m_z = combined_ions$m_z+mass_shift_tolerence/combined_ions$Charge)
  }
  
  if(C13_label_state == TRUE){
    if(tolower(ppm_or_Da) == "ppm"){
      by_ions_long_1<-cbind(by_ions_long_1, low_m_z = by_ions_long_1$m_z*(1-mass_shift_tolerence*10^-6))
      by_ions_long_1<-cbind(by_ions_long_1, high_m_z = by_ions_long_1$m_z*(1+mass_shift_tolerence*10^-6))
      combined_ions_1<-cbind(combined_ions_1, low_m_z = combined_ions_1$m_z*(1-mass_shift_tolerence*10^-6))
      combined_ions_1<-cbind(combined_ions_1, high_m_z = combined_ions_1$m_z*(1+mass_shift_tolerence*10^-6))
      by_ions_long_2<-cbind(by_ions_long_2, low_m_z = by_ions_long_2$m_z*(1-mass_shift_tolerence*10^-6))
      by_ions_long_2<-cbind(by_ions_long_2, high_m_z = by_ions_long_2$m_z*(1+mass_shift_tolerence*10^-6))
      combined_ions_2<-cbind(combined_ions_2, low_m_z = combined_ions_2$m_z*(1-mass_shift_tolerence*10^-6))
      combined_ions_2<-cbind(combined_ions_2, high_m_z = combined_ions_2$m_z*(1+mass_shift_tolerence*10^-6))
    }else{
      by_ions_long_1<-cbind(by_ions_long_1, low_m_z = by_ions_long_1$m_z-mass_shift_tolerence/by_ions_long_1$Charge)
      by_ions_long_1<-cbind(by_ions_long_1, high_m_z = by_ions_long_1$m_z+mass_shift_tolerence/by_ions_long_1$Charge)
      combined_ions_1<-cbind(combined_ions_1, low_m_z = combined_ions_1$m_z-mass_shift_tolerence/combined_ions_1$Charge)
      combined_ions_1<-cbind(combined_ions_1, high_m_z = combined_ions_1$m_z+mass_shift_tolerence/combined_ions_1$Charge)
      by_ions_long_2<-cbind(by_ions_long_2, low_m_z = by_ions_long_2$m_z-mass_shift_tolerence/by_ions_long_2$Charge)
      by_ions_long_2<-cbind(by_ions_long_2, high_m_z = by_ions_long_2$m_z+mass_shift_tolerence/by_ions_long_2$Charge)
      combined_ions_2<-cbind(combined_ions_2, low_m_z = combined_ions_2$m_z-mass_shift_tolerence/combined_ions_2$Charge)
      combined_ions_2<-cbind(combined_ions_2, high_m_z = combined_ions_2$m_z+mass_shift_tolerence/combined_ions_2$Charge)
    }
  }
  
  
  
  by_ions_long$match_intensity<-win_intensity(by_ions_long$low_m_z, by_ions_long$high_m_z,
                                              mgf_ions$m_z, mgf_ions$Intensity)
  
  if(C13_label_state == TRUE){
    by_ions_long$match_intensity<-by_ions_long$match_intensity+
      win_intensity(by_ions_long_1$low_m_z, by_ions_long_1$high_m_z,
                    mgf_ions$m_z, mgf_ions$Intensity)+
      win_intensity(by_ions_long_2$low_m_z, by_ions_long_2$high_m_z,
                    mgf_ions$m_z, mgf_ions$Intensity)
  }
  
  
  mgf_ions$Match<-cover_count(mgf_ions$m_z, by_ions_long$low_m_z, by_ions_long$high_m_z)
  
  mgf_ions$All_Match<-cover_count(mgf_ions$m_z, combined_ions$low_m_z, combined_ions$high_m_z)
  
  if(C13_label_state == TRUE){
    mgf_ions$Match<-cover_count(mgf_ions$m_z,
                                c(by_ions_long$low_m_z, by_ions_long_1$low_m_z,
                                  by_ions_long_2$low_m_z),
                                c(by_ions_long$high_m_z, by_ions_long_1$high_m_z,
                                  by_ions_long_2$high_m_z))
    mgf_ions$All_Match<-cover_count(mgf_ions$m_z,
                                    c(combined_ions$low_m_z, combined_ions_1$low_m_z,
                                      combined_ions_2$low_m_z),
                                    c(combined_ions$high_m_z, combined_ions_1$high_m_z,
                                      combined_ions_2$high_m_z))
  }
  
  by_ions_coverage<-round(length(which(by_ions_long$match_intensity>0))/nrow(by_ions_long)*100,2)
  PSM_score_by<-round(sum(mgf_ions$Intensity[which(mgf_ions$Match>0)])/sum(mgf_ions$Intensity)*100,2)
  PSM_score_all<-round(sum(mgf_ions$Intensity[which(mgf_ions$All_Match>0)])/sum(mgf_ions$Intensity)*100,2)
  return(c(by_ions_long = by_ions_long,by_ions_coverage = by_ions_coverage,PSM_score_by = PSM_score_by,PSM_score_all = PSM_score_all))
}
#' 选定理论离子的匹配。
#' Match a selected set of theoretical ions.
#'
#' 对给定的一列理论 m/z（通常为 all_ions_position 的 b/y/a/丢失离子）做容差匹配，
#' 返回逐离子匹配强度向量与打标后的实验谱表（Match/All_Match 列）。
#' Matches a given vector of theoretical m/z values (typically the b/y/a/loss
#' ions from all_ions_position) within tolerance and returns the per-ion match
#' intensities plus the annotated observed table (Match/All_Match columns).
#'
#' @param therotical_ions
#'   数值向量，理论离子 m/z。
#'   Numeric vector, the theoretical ion m/z values.
#' @param mgf_ions
#'   数据框，实验谱的 m/z-强度表。
#'   A data frame, the observed m/z-intensity table.
#' @param mass_shift_tolerence
#'   数值，质量容差，默认为 20。
#'   Numeric, the mass tolerance; defaults to 20.
#' @param ppm_or_Da
#'   字符向量 c("ppm","Da")，容差单位。
#'   Character vector c("ppm","Da"), the tolerance unit.
#' @return
#'   返回 list(match_result=逐离子强度向量, mgf_ions=带 Match 列的实验谱表)。
#'   A list(match_result = per-ion intensity vector, mgf_ions = observed table
#'   with the Match column).
#' @examples
#' spec <- data.frame(m_z = c(226.11, 227.11), Intensity = c(1000, 900))
#' ions_match_selected(c(226.11, 370.17), spec, 20, c("ppm", "Da"))
#' @export
ions_match_selected<-function(therotical_ions = all_ions_position(...)$ions_table$m_z,
                              mgf_ions = NULL,
                     mass_shift_tolerence = 20, ppm_or_Da = c("ppm","Da")){
  if (is.null(mgf_ions)) stop("mgf_ions required (the legacy read_large_mgf_auto default was removed )")
  if(!tolower(ppm_or_Da)%in%c("ppm","da")){cat("Error mass shift tolerence!");return(NA)}
  if(!is.numeric(mass_shift_tolerence)){cat("Error mass shift tolerence!");return(NA)}
  
  if(tolower(ppm_or_Da) == "ppm"){
    low_m_z <- therotical_ions*(1-mass_shift_tolerence*10^-6)
    high_m_z <- therotical_ions*(1+mass_shift_tolerence*10^-6)
  }else{
    low_m_z <- therotical_ions - mass_shift_tolerence
    high_m_z <- therotical_ions + mass_shift_tolerence
    }
  if(!"Match"%in%colnames(mgf_ions)){
    mgf_ions<-cbind(mgf_ions,Match = 0)
  }
  match_result <- win_intensity(low_m_z, high_m_z, mgf_ions$m_z, mgf_ions$Intensity)
  covered <- cover_count(mgf_ions$m_z, low_m_z, high_m_z) > 0
  mgf_ions$Match[covered] <- 1
  
  return(list(match_result = match_result, mgf_ions = mgf_ions))
}
#' pFind .dta 标题 → raw + scan。
#' pFind .dta title -> raw + scan.
#'
#' 标题格式：'RAW.SCAN.SCAN.CHARGE.RANK.dta'（RAW 名不含点），
#' 返回 c(raw=..., scan=...) 命名向量，供按 raw 分组与标记谱查找。
#' Title format: 'RAW.SCAN.SCAN.CHARGE.RANK.dta' (the RAW stem contains no
#' dots); returns a c(raw=..., scan=...) named vector for raw grouping and
#' labeled-spectrum lookup.
#'
#' @param title
#'   字符串，pFind .dta 风格谱图标题。
#'   Character string, a pFind .dta-style spectrum title.
#' @return
#'   返回命名向量 c(raw=..., scan=...)。
#'   A named vector c(raw=..., scan=...).
#' @examples
#' title_raw_scan("RAW1.1234.1234.2.0.dta")
#' @export
title_raw_scan <- function(title) {
  p <- strsplit(title, ".", fixed = TRUE)[[1]]
  c(raw = p[1], scan = p[2])
}
#' pFind .dta 标题 → 'RAW.SCAN' 键。
#' pFind .dta title -> 'RAW.SCAN' key.
#'
#' 取标题前两段组成 'RAW.SCAN' 键，用于同 raw＋scan 的谱块查找。
#' Takes the first two fields as a 'RAW.SCAN' key for same-raw-plus-scan
#' block lookup.
#'
#' @param title
#'   字符串，pFind .dta 风格谱图标题。
#'   Character string, a pFind .dta-style spectrum title.
#' @return
#'   返回 'RAW.SCAN' 字符串。
#'   A 'RAW.SCAN' string.
#' @examples
#' title_scan_key("RAW1.1234.1234.2.0.dta")
#' @export
# 同 scan 多 MS2（pParse 共洗脱/多电荷拆分）
# 的回退选择——前体离子质量差 ppm 最小者优先；ppm 并列再取价态一致者；
# 无质量信息时退回"价态一致优先，否则首个"（旧口径）。psm_mh 为 Exp.MH.
# （单电荷 M+H+），按 psm 电荷折算 m/z 后比较。
# Fallback selection among same-scan MS2 blocks: minimal precursor-ppm first,
# tie-broken by matching charge; without mass info, same-charge-then-first.
pick_best_src <- function(cands, chg, psm_mh = NA_real_) {
  if (length(cands) == 1) return(cands[[1]])
  chg_n <- suppressWarnings(as.numeric(chg))
  if (is.na(chg_n)) chg_n <- 2
  mz <- vapply(cands, function(x) suppressWarnings(as.numeric(x[["mz"]])), 0)
  if (!is.na(psm_mh) && psm_mh > 0 && any(!is.na(mz) & mz > 0)) {
    psm_mz <- (psm_mh + (chg_n - 1) * 1.007276) / chg_n
    z <- vapply(cands, function(x) suppressWarnings(as.numeric(x[["charge"]])), 0)
    z <- ifelse(is.na(z) | z < 1, 1, z)
    ppm <- abs(mz - psm_mz) / psm_mz * 1e6
    ppm[is.na(ppm)] <- Inf
    same_z <- abs(z - chg_n) < 1e-9
    # 1 ppm 内视为并列（设计约定"多个候选再比价态"）；并列集内价态一致优先、
    # 再取更小 ppm
    tie <- ppm <= min(ppm) + 1
    idx <- which(tie)
    ord2 <- idx[order(!same_z[idx], ppm[idx])]
    return(cands[[ord2[1]]])
  }
  if (!is.na(chg_n))
    for (x in cands) if (suppressWarnings(as.numeric(x[["charge"]])) == chg_n) return(x)
  cands[[1]]
}
title_scan_key <- function(title) {
  p <- strsplit(title, ".", fixed = TRUE)[[1]]
  paste0(p[1], ".", p[2])
}
#' 把参数文件条目展开为排序文件列表。
#' Expand parameter-file entries into a sorted file list.
#'
#' 每条目可为文件路径（直接收）或文件夹路径（展开为其中全部匹配 pattern 的文件并排序）；
#' 缺失路径或无匹配文件的文件夹报错并列出。
#' Each entry is either a file path (taken as-is) or a folder path (expanded to
#' all files inside matching the pattern, sorted); missing paths or folders
#' without matches raise an error listing them.
#'
#' @param entries
#'   字符向量，参数文件某键的条目（每行一个）。
#'   Character vector, the entries of a parameter-file key (one per line).
#' @param pattern
#'   字符串，文件名正则（如 "\\.mgf$"）。
#'   Character string, the file-name pattern (e.g. "\\.mgf$").
#' @return
#'   返回排序后的文件绝对路径字符向量。
#'   A sorted character vector of absolute file paths.
#' @keywords internal
resolve_file_list <- function(entries, pattern) {
  out <- character(0)
  bad <- character(0)
  for (e in entries) {
    if (dir.exists(e)) {
      fs <- sort(list.files(e, pattern = pattern, full.names = TRUE))
      if (length(fs) == 0) bad <- c(bad, e)
      out <- c(out, fs)
    } else if (file.exists(e)) {
      out <- c(out, e)
    } else {
      bad <- c(bad, e)
    }
  }
  if (length(bad) > 0) {
    stop("entries missing or without matching files: ", paste(bad, collapse = "; "))
  }
  unique(out)
}
#' 把参数文件 mgf 条目展开为排序的 MGF 文件列表。
#' Expand parameter-file mgf entries into a sorted MGF file list.
#'
#' resolve_file_list 的 MGF 特化（pattern 固定为 "\\.mgf$"）。
#' The MGF specialization of resolve_file_list (pattern fixed to "\\.mgf$").
#'
#' @param entries
#'   字符向量，参数文件的 mgf 条目。
#'   Character vector, the mgf entries of a parameter file.
#' @return
#'   返回排序后的 MGF 文件绝对路径。
#'   A sorted character vector of absolute MGF file paths.
#' @keywords internal
resolve_mgf_list <- function(entries) resolve_file_list(entries, "\\.mgf$")

# gf 键可直接填 raw 文件/文件夹——第一步
# 用 pParse（随 pFind 安装；pfind_exe 是本软件必装依赖，零额外安装）自动转
# MGF；无 pFind 时回退 msconvert（msconvert_exe 键）；两者皆无则明确报错。
# pParse 输出标题（raw.scan.scan.charge.fraction.dta）＝本管线原生 File_Name
# 约定，与手工经 pParse 转 MGF 的既有口径完全一致（无新语义）。
# 调用契约（实测）：pParse.exe <cfg路径>（位置参数；cfg 为 GBK＋CRLF；进程
# 对失败也返回 0——以产物 MGF 是否落盘为成功判据）。
#: the mgf key accepts raw files/folders directly - step 1 converts
# them with pParse (ships with pFind, a required dependency; nothing extra to
# install); falls back to msconvert (msconvert_exe key); a clear error when
# neither exists. pParse titles follow this pipeline's native File_Name
# convention, identical to the established manual-pParse workflow.
# Contract (verified live): pParse.exe <cfg-path> (positional; GBK + CRLF;
# the process returns 0 even on failure - success is judged by output MGfs).
convert_raw_entries_to_mgf <- function(mgf_entries, cfg) {
  ent <- mgf_entries[trimws(mgf_entries) != ""]
  direct_raws <- ent[grepl("\\.raw$", ent, ignore.case = TRUE) & file.exists(ent)]
  dirs <- ent[dir.exists(ent)]
  dir_raws <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = "\\.raw$", ignore.case = TRUE,
               recursive = TRUE, full.names = TRUE)), use.names = FALSE)
  raws <- unique(gsub("\\\\", "/", normalizePath(c(direct_raws, dir_raws), winslash = "/")))
  if (!length(raws)) return(resolve_mgf_list(ent))
  out_dir <- file.path(cfg$config_dir, "pparse_mgf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("[第一步] 检测到 ", length(raws), " 个 raw 文件，自动转换 MGF（pParse）\n", sep = "")
  cat("[step 1] ", length(raws), " raw file(s) detected; converting to MGF (pParse)\n", sep = "")
  # 定位 pParse：<pfind_exe>/bin/pParse.exe（pfind_exe 约定为安装根；兼容直指 bin）
  pf_root <- if (!is.null(cfg$pfind_exe) && length(cfg$pfind_exe) > 0)
    trimws(tail(cfg$pfind_exe, 1)) else ""
  pparse <- ""
  if (nzchar(pf_root)) {
    cand <- file.path(gsub("\\\\", "/", pf_root), "bin", "pParse.exe")
    cand2 <- if (grepl("/bin$", cand)) file.path(dirname(cand), "pParse.exe") else cand
    pparse <- if (file.exists(cand)) cand else if (file.exists(cand2)) cand2 else ""
  }
  logf <- file.path(out_dir, "pparse_run.log")
  if (nzchar(pparse)) {
    cfg_lines <- c("# auto-generated by isoPSV (step-1 raw->MGF)",
                   "[Basic Options]",
                   paste0("datanum = ", length(raws)),
                   paste0("datapath", seq_along(raws), " = ", raws),
                   "[Advanced Options]",
                   "input_format = raw",
                   "output_mgf = 1",
                   "output_pf = 0",
                   "output_msb = 0",
                   "delete_msn = 0",
                   "co-elute = 1",
                   paste0("outputpath = ", gsub("\\\\", "/", out_dir)),
                   "rewrite_files = 0")
    pcfg <- file.path(out_dir, "pparse_auto.cfg")
    writeLines(cfg_lines, pcfg, useBytes = TRUE, sep = "\r\n")
    # cfg 需 GBK（pParse 中文版读取约定）：UTF-8 键均为 ASCII，仅路径可能含
    # 非 ASCII——统一转 GBK 写
    iconv_lines <- iconv(cfg_lines, from = "UTF-8", to = "GBK")
    if (any(is.na(iconv_lines))) iconv_lines <- cfg_lines
    con <- file(pcfg, open = "wb"); writeLines(iconv_lines, con, sep = "\r\n"); close(con)
    old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
    setwd(dirname(pparse))
    st <- suppressWarnings(system2(pparse, shQuote(pcfg), stdout = logf, stderr = logf))
    setwd(old_wd)
    mgs <- list.files(out_dir, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE)
    if (!length(mgs))
      stop("[第一步] pParse 转换失败（无 MGF 产出；exit ", st, "）：", basename(raws[1]),
           "（日志：", logf, "）\n",
           "[step 1] pParse conversion produced no MGF (exit ", st, "): ",
           basename(raws[1]), " (log: ", logf, ")")
    cat("[第一步] pParse 转换完成：", length(mgs), " 个 MGF -> ", out_dir, "\n", sep = "")
    cat("[step 1] pParse conversion done: ", length(mgs), " MGF(s) -> ", out_dir, "\n", sep = "")
    # 已有手工 MGF 条目与转换产物合并去重
    manual <- resolve_mgf_list(ent[!grepl("\\.raw$", ent, ignore.case = TRUE) & !dir.exists(ent)])
    return(sort(unique(c(manual, gsub("\\\\", "/", mgs)))))
  }
  # 回退：msconvert（检测键存在时逐个转换）
  msc <- if (!is.null(cfg$msconvert_exe) && length(cfg$msconvert_exe) > 0)
    trimws(tail(cfg$msconvert_exe, 1)) else ""
  if (nzchar(msc) && file.exists(msc)) {
    for (r in raws) {
      suppressWarnings(system2(msc, c(shQuote(r), "--mgf", "-o", shQuote(out_dir)),
                               stdout = logf, stderr = logf))
    }
    mgs <- list.files(out_dir, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE)
    if (length(mgs)) {
      cat("[第一步] msconvert 转换完成：", length(mgs), " 个 MGF\n", sep = "")
      manual <- resolve_mgf_list(ent[!grepl("\\.raw$", ent, ignore.case = TRUE) & !dir.exists(ent)])
      return(sort(unique(c(manual, gsub("\\\\", "/", mgs)))))
    }
  }
  stop("[第一步] mgf 键填了 raw 文件但无法转换：未找到 pParse（pfind_exe/bin/pParse.exe）",
       " 与 msconvert（msconvert_exe）\n",
       "请安装 pFind Studio（http://www.pfind.net，软件本身需要它）或 ProteoWizard",
       "（https://proteowizard.sourceforge.io），或先把 raw 转成 MGF 再填 mgf 键。\n",
       "[step 1] raw files given in the mgf key but no converter found: neither",
       " pParse (pfind_exe/bin/pParse.exe) nor msconvert (msconvert_exe)\n",
       "Install pFind Studio (http://www.pfind.net - required by this software anyway)",
       " or ProteoWizard (https://proteowizard.sourceforge.io), or convert the raw",
       " files to MGF before filling the mgf key.")
}
#' 读 pFind.cfg 风格参数文件。
#' Read a pFind.cfg-style parameter file.
#'
#' 制表符分隔的键值行＋'#' 注释；重复键的值按序收集；
#' 配置所在目录以 config_dir 返回（相对路径以其为基准解析）。
#' Tab-separated key/value lines with '#' comments; values of repeated keys are
#' collected in order; the config folder is returned as config_dir (relative
#' paths resolve against it).
#'
#' @param path
#'   字符串，参数文件路径。
#'   Character string, the parameter-file path.
#' @return
#'   返回列表（每个键一个字符向量）＋config_dir 元素。
#'   A list (one character vector per key) plus the config_dir element.
#' @keywords internal
read_param_file <- function(path) {
  # N-17先给可读的文件不存在提示（早期实现 readBin(n=NA) 报含糊错误）
  if (!file.exists(path)) stop("参数文件不存在\nparameter file not found: ", path)
  # P3-5编码安全读入——UTF-8（含 BOM）按 UTF-8 解码并剥 BOM，
  # 非 UTF-8 经 iconv(GBK→UTF-8) 转码，失败才退回本地编码 readLines；二者都会
  # 正确保留中文。（同步：注释此前仍写"非 UTF-8 退回本地编码"，
  # 与实际 iconv 分支不符。）
  # P3-5 encoding-safe read - UTF-8 (with BOM) is decoded as
  # UTF-8 with the BOM stripped; non-UTF-8 goes through iconv(GBK->UTF-8),
  # falling back to native-encoding readLines only on failure; both preserve
  # Chinese correctly. (sync: the comment used to say non-UTF-8
  # falls back to native readLines, which did not match the iconv branch.)
  # N-4rawToChar 结果无编码标记，必须 Encoding<-"UTF-8" 再
  # enc2native，否则 GBK 代码页下 UTF-8 中文路径仍被按本地字节串解释。
  # N-4 rawToChar results carry no encoding mark; mark UTF-8
  # then enc2native, or Chinese paths in UTF-8 configs are still interpreted
  # as native bytes under a GBK codepage.
  raw <- readBin(path, "raw", n = file.info(path)$size)
  has_bom <- length(raw) >= 3 && identical(raw[1:3], as.raw(c(0xEF, 0xBB, 0xBF)))
  # 注意：无 BOM 时必须保留全量字节——raw[-seq_len(0)] 在 R 里是"取空集"
  # 而非"取全部"（曾把无 BOM 文件整体清空、全部误入回退分支）
  # careful: without a BOM the full bytes must be kept - raw[-seq_len(0)] is
  # an EMPTY selection in R, not "everything" (it once emptied every non-BOM
  # file into the fallback branch)
  body <- if (has_bom) raw[-seq_len(3L)] else raw
  txt <- tryCatch(suppressWarnings(rawToChar(body)), error = function(e) "")
  if (has_bom || (nzchar(txt) && isTRUE(all(validUTF8(txt))))) {
    Encoding(txt) <- "UTF-8"
    txt <- enc2native(txt)
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  } else {
    # P0-BUTF-8 原生构建的 R（如工作站
    # R 4.6.1）里 readLines 读 GBK 文件返回非法 UTF-8 串，随后的 trimws 直接
    # 报 "input string 1 is invalid UTF-8"；改为对已读原始字节 iconv(GBK→UTF-8)，
    # 失败才退回 readLines 并警告。
    txt_gbk <- tryCatch(suppressWarnings(iconv(txt, from = "GBK", to = "UTF-8")),
                        error = function(e) NA_character_)
    if (!is.na(txt_gbk) && nzchar(txt_gbk)) {
      lines <- strsplit(txt_gbk, "\n", fixed = TRUE)[[1]]
    } else {
      warning("参数文件无法按 UTF-8/GBK 解码，按本地编码读取（可能出现乱码）：", path,
              call. = FALSE)
      lines <- readLines(path, warn = FALSE)
    }
  }
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  if (length(lines) == 0) stop("no active entries in ", path)
  #：缩进续行并入前置键（重写为重复键行）——多行键（如 4 标记
  # label_type）的 \t/空格缩进续行写法此前被解析成空键垃圾行、标记静默丢失
  #（fixmod 只剩第一行、Stage A 理论标记不全）。重写后走既有重复键向量机制，
  # 下游 label_type 消费方零改动。
  #: indented continuation lines fold into the preceding key (rewritten
  # as repeated key rows) - multi-row keys (e.g. a 4-label label_type) written
  # with tab/space-indented continuation lines used to become empty-key junk
  # lines, silently dropping labels (fixmod kept only the first row, stage-A
  # theoretical labels incomplete). Rewritten rows ride the existing
  # duplicate-key vector mechanism; downstream consumers need no change.
  cont <- grepl("^[ \t\r]", lines)
  if (any(cont)) {
    if (cont[1])
      warning("参数文件首行即续行（无前置键），已忽略 / leading continuation line dropped: ",
              trimws(lines[1]), call. = FALSE)
    drop <- logical(length(lines))
    for (i in which(cont)) {
      j <- i - 1L
      while (j >= 1L && drop[j]) j <- j - 1L
      if (j >= 1L) {
        key_j <- strsplit(lines[j], "\t", fixed = TRUE)[[1]][1]
        lines[i] <- paste0(key_j, "\t", trimws(lines[i]))
        cont[i] <- FALSE
      } else drop[i] <- TRUE
    }
    lines <- lines[!drop]
  }
  kv <- strsplit(lines, "\t", fixed = TRUE)
  keys <- vapply(kv, function(x) trimws(x[1]), character(1))
  # P3-7值保留第 3 段及之后内容（含制表符的值不再被截断）；
  # 无制表符的行给出空格推断提示并按连续空格切分；重复键 warning 说明取值规则。
  # P3-7 keep segment 3 onwards (tab-containing values are no
  # longer truncated); tab-less lines are split on 2+ spaces with a warning;
  # duplicate keys warn about the take-all-values rule.
  vals <- vapply(kv, function(x) {
    if (length(x) >= 2) trimws(paste(x[-1], collapse = "\t")) else ""
  }, character(1))
  # P1值被切出 3 段以上时提示——第 3 段起并入值，
  # 尾列注释/多余列会被并进路径导致"文件不存在"；提示用户值内勿用制表符。
  if (any(lengths(kv) > 2)) {
    warning("参数文件存在值含制表符的行（第 3 段起并入值；若为尾列注释请删除，值内勿用制表符）：",
            paste(which(lengths(kv) > 2), collapse = ", "), call. = FALSE)
  }
  if (any(!grepl("\t", lines, fixed = TRUE))) {
    notab <- which(!grepl("\t", lines, fixed = TRUE))
    sp <- strsplit(trimws(lines[notab]), "[ \t]{2,}")
    has2 <- lengths(sp) >= 2
    if (any(has2)) {
      warning("参数文件 ", basename(path), " 有 ", sum(has2),
              " 行无制表符，已按连续空格推断键值分隔（请改用制表符）",
              call. = FALSE)
      ii <- notab[has2]
      keys[ii] <- vapply(sp[has2], function(x) trimws(x[1]), character(1))
      vals[ii] <- vapply(sp[has2], function(x) trimws(paste(x[-1], collapse = " ")),
                         character(1))
    }
  }
  dupk <- unique(keys[duplicated(keys)])
  if (length(dupk) > 0)
    # N-8文案与实现对齐——实现保留全部同名值为向量，各消费方
    # 自行取 [1]/tail()（如 user_step2 的 cfg_get 用 tail 取最后一个非空）
    warning("参数文件存在重复键（同名值全部保留为向量，各消费方自行取第一个或最后一个非空）：",
            paste(dupk, collapse = ", "), call. = FALSE)
  uk <- unique(keys)
  out <- setNames(lapply(uk, function(k) vals[keys == k]), uk)
  out$config_dir <- dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
  out
}
#' 把参数文件某键的条目按配置目录解析为绝对路径。
#' Resolve the entries of a parameter-file key against the config folder.
#'
#' 相对路径以参数文件所在目录（config_dir）为基准转绝对路径。
#' Relative paths are made absolute against the parameter file's folder
#' (config_dir).
#'
#' @param cfg
#'   列表，read_param_file 的输出。
#'   A list, the output of read_param_file.
#' @param key
#'   字符串，要解析的键名。
#'   Character string, the key to resolve.
#' @return
#'   返回绝对路径字符向量。
#'   A character vector of absolute paths.
#' @keywords internal
resolve_cfg_paths <- function(cfg, key) {
  v <- cfg[[key]]
  if (is.null(v) || length(v) == 0) return(character(0))
  vapply(v, function(p) {
    if (grepl("^[A-Za-z]:[/\\\\]", p) || startsWith(p, "/")) p
    else file.path(cfg$config_dir, p)
  }, character(1), USE.NAMES = FALSE)
}
#' 步骤 2 策略集合指纹同一配置的"策略集合"身份——
#' search_file 原始条目（经路径解析）+ 类型（显式或命名规则推断），排序去重。
#' user_run（跳过判据）与 user_step2（写侧车 merge/peptide.txt.strategies）
#' 共用同一实现，保证两侧逐字节一致；策略集合变化时步骤 2 重跑（引擎结果
#' 经参数指纹复用，仅重算合并）。
#' Step-2 strategy-set fingerprint the identity of
#' the configured strategy set - resolved search_file entries + types (explicit
#' or inferred by the naming rules), sorted and deduplicated. user_run (the
#' skip criterion) and user_step2 (the sidecar writer) share this one
#' implementation so the two sides agree byte-for-byte; a changed strategy set
#' reruns step 2 (engine results are reused via their parameter fingerprints,
#' only the merge is recomputed).
strategy_identity_lines <- function(cfg2) {
  search_files <- resolve_cfg_paths(cfg2, "search_file")
  sf_type_cfg <- cfg2$search_file_type
  sf_type_cfg <- if (!is.null(sf_type_cfg)) sf_type_cfg[!is.na(sf_type_cfg) & sf_type_cfg != ""] else character(0)
  if (length(sf_type_cfg) > 0 && length(sf_type_cfg) != length(search_files)) {
    stop("search_file_type 行数须与 search_file 一致\nsearch_file_type must match search_file")
  }
  infer_strategy_type <- function(sf) {
    sf_dir <- tolower(basename(dirname(dirname(sf))))
    sf_name <- tolower(basename(sf))
    if (grepl("denovo", sf_dir, fixed = TRUE)) return("punifind_denovo")
    if (grepl("boost", sf_dir, fixed = TRUE)) return("msfragger_boost")
    if (grepl("pnovo", sf_dir, fixed = TRUE)) return("pnovo")
    if (grepl("restrict", sf_dir, fixed = TRUE)) return("pfind_restrict")
    if (grepl("standard", sf_dir, fixed = TRUE)) return("pfind_standard")
    if (grepl("_punifind", sf_name, fixed = TRUE)) return("punifind_rescore")
    "pfind_open"
  }
  search_types <- if (length(sf_type_cfg) > 0) sf_type_cfg else {
    if (length(search_files) > 0) vapply(search_files, infer_strategy_type, character(1))
    else character(0)
  }
  out <- sort(unique(paste0(as.character(search_types), "\t",
                            gsub("\\\\", "/", as.character(search_files)))))
  #：候选长度阈值入清单——改 min_pep_len 须触发步骤 2 重跑合并
  # （搜库结果按参数指纹复用，仅重算合并表与 peptide.txt）
  #: the candidate length gate joins the manifest - a min_pep_len
  # change must re-run the step-2 merge (search results are reused via the
  # parameter fingerprints; only the merged table and peptide.txt recompute)
  mpl <- 7L
  try({
    bc <- resolve_cfg_paths(cfg2, "base_config")[1]
    if (nzchar(bc) && file.exists(bc)) {
      b <- read_param_file(bc)
      if (!is.null(b$min_pep_len) && length(b$min_pep_len) > 0) {
        v <- suppressWarnings(as.integer(tail(b$min_pep_len, 1)))
        if (!is.na(v) && v >= 4L) mpl <- v
      }
    }
  }, silent = TRUE)
  #：pFind 臂开关入清单——换臂组合（只跑 open / 只跑 restrict /
  # 加跑 standard）须触发步骤 2 重跑合并；读取口径与主流程一致
  # （config_step1 优先，值须为 0/1，默认 open=1/restrict=1/standard=0）
  #: the pFind arm switches join the manifest - changing the arm
  # combo must re-run the step-2 merge; same precedence as the main flow
  # (step1 first, 0/1 values only, defaults 1/1/0).
  arm_lines <- character(0)
  try({
    if (!exists("b", inherits = FALSE)) {
      bc2 <- resolve_cfg_paths(cfg2, "base_config")[1]
      b <- if (nzchar(bc2) && file.exists(bc2)) read_param_file(bc2) else list()
    }
    arm_get <- function(key, default) {
      for (s in list(b, cfg2)) {
        if (!is.null(s[[key]]) && length(s[[key]]) > 0 &&
            trimws(as.character(tail(s[[key]], 1))) %in% c("0", "1"))
          return(trimws(as.character(tail(s[[key]], 1))))
      }
      default
    }
    arm_lines <- paste0("pfind_arms=open=", arm_get("pfind_open", "1"),
                        "|restrict=", arm_get("pfind_restrict", "1"),
                        "|standard=", arm_get("pfind_standard", "0"))
  }, silent = TRUE)
  c(out, paste0("min_pep_len=", mpl), arm_lines)
}
#' pUniFind 配置键优先级读取config_step1（step1_cfg）
#' 优先于 config_step2（step2_cfg）；validate 过滤无效值（如越界模式码）；两处
#' 都无有效值返回 default。返回 list(value=, source=)（source 为
#' "config_step1"/"config_step2"/""），user_step2 据此打印生效值与来源。
#' pUniFind config-key precedence read 
#' config_step1 wins over config_step2; validate drops invalid values; the
#' default is returned when neither side holds a valid value. Returns
#' list(value=, source=) with source in "config_step1"/"config_step2"/"" so the
#' caller can print the effective value and where it came from.
# v1.0参考库匹配三口径——accession 命中、
# 序列级匹配（大写＋I/L 归一子串；六框翻译库 orf| accession 与参考库
# WP_/NP_/tr| 不同体系时仅比 accession 会整侧漏判）、污染库排除（CON_ 前缀
# accession 一律视为非新肽段）。fa_file 缺失时只做 accession 与 CON_ 判定。
fa_annotated_proteins <- function(proteins, seqs, fa_file) {
  fa_acc <- character(0); fa_all <- ""
  if (!is.null(fa_file) && length(fa_file) > 0 && nzchar(fa_file[1]) &&
      file.exists(fa_file[1])) {
    fa_obj <- tryCatch(read_fasta(fa_file[1]), error = function(e) character(0))
    fa_acc <- unique(names(fa_obj))
    fa_all <- paste(toupper(unname(as.character(fa_obj))), collapse = "\n")
  }
  n <- max(length(proteins), length(seqs))
  proteins <- rep_len(as.character(proteins), n)
  seqs <- rep_len(as.character(seqs), n)
  vapply(seq_len(n), function(i) {
    pr <- proteins[i]
    if (is.na(pr) || !nzchar(pr) || pr %in% c("-", "NA")) return(FALSE)
    p <- trimws(strsplit(pr, "/", fixed = TRUE)[[1]])
    p <- p[nzchar(p)]
    if (length(p) == 0) return(FALSE)
    if (any(startsWith(p, "CON_"))) return(TRUE)
    if (length(fa_acc) > 0 && any(p %in% fa_acc)) return(TRUE)
    s <- toupper(seqs[i])
    if (nzchar(fa_all) && !is.na(s) && nzchar(s))
      if (grepl(s, fa_all, fixed = TRUE) ||
          grepl(gsub("I", "L", s), fa_all, fixed = TRUE)) return(TRUE)
    FALSE
  }, logical(1))
}

cfg_get_preferred <- function(step1_cfg, step2_cfg, key, default = NULL,
                              validate = NULL) {
  # v1.0step1 与 step2 同键不同值时显式告警——
  # 旧行为静默取 step1，Macaca 把正确的新路径写在 config_step2 却被
  # config_step1 的旧值覆盖、GPU 段跑旧 CLI 崩掉白等一夜。告警不改取值口径
  # （仍 step1 优先），只把覆盖变成可见诊断。
  v1 <- step1_cfg[[key]]
  v2 <- step2_cfg[[key]]
  v1v <- if (!is.null(v1) && length(v1) > 0) as.character(tail(v1, 1)) else NA_character_
  v2v <- if (!is.null(v2) && length(v2) > 0) as.character(tail(v2, 1)) else NA_character_
  if (!is.na(v1v) && !is.na(v2v) && !identical(v1v, v2v)) {
    cat("[配置告警] ", key, " 在 config_step1 与 config_step2 的值不同：\n",
        "  step1 = ", v1v, "\n  step2 = ", v2v,
        "\n  按口径取 step1（config_step1 优先）。若 step2 的值为新路径，请同时更正 config_step1。\n",
        sep = "")
    cat("[config warning] ", key, " differs between config_step1 and config_step2: ",
        "step1=", v1v, " step2=", v2v,
        "; step1 wins by rule - if step2 holds the new value, fix config_step1 too.\n", sep = "")
  }
  for (nm in c("config_step1", "config_step2")) {
    src <- if (nm == "config_step1") step1_cfg else step2_cfg
    v <- src[[key]]
    if (!is.null(v) && length(v) > 0) {
      val <- tail(v, 1)
      if (is.null(validate) || isTRUE(validate(val))) {
        return(list(value = val, source = nm))
      }
    }
  }
  list(value = default, source = "")
}
# ---- MGF → mzML 本地化转换------------------------
# 背景：Sage/MetaMorpheus/MSFragger 只读 mzML，而管线输入本身是 MGF；此前依赖
# 外部 msconvert（ProteoWizard，未装机即硬停），代码目录的
# mgf_to_mzml_builtin.py 从未被调用。本函数在 R 内完成 MGF→mzML(1.1.0)：
# 每谱 MS1 谱、m/z 与 intensity 双 64 位小端 base64 数组、scan start time/
# selected ion m/z/charge/title 保留，并写规范 indexList 字节偏移索引。
# 注：vendor RAW→MGF 属管线上游（R 无法解析专有二进制），仍由用户侧 msconvert
# 完成；本函数覆盖的是管线内部的 MGF→mzML 这一唯一必需转换。
# ---------------------------------------------------------------------------
b64_alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1]]
#' base64 编码（raw 向量 → 字符串；纯 base R，向量化 6 位分组）
base64_encode_raw <- function(bytes) {
  if (!length(bytes)) return("")
  m8 <- matrix(as.integer(rawToBits(bytes)), nrow = 8L)   # 每字节一列、位 LSB 在前
  m8 <- m8[8:1,, drop = FALSE]                            # 翻转为 MSB 在前
  bits <- as.logical(as.vector(m8))                        # 按字节序展平
  npad <- (6L - length(bits) %% 6L) %% 6L
  bits <- c(bits, rep(FALSE, npad))
  m6 <- matrix(bits, nrow = 6L)
  idx <- colSums(m6 * (2^(5:0)))                           # 6 位组 → 0..63
  out <- paste0(b64_alphabet[idx + 1L], collapse = "")
  paste0(out, strrep("=", (4L - nchar(out) %% 4L) %% 4L))
}
#' mgf_to_mzml：MGF 文件 → 规范最小 mzML（indexedmzML + indexList 字节偏移）
# v1.0内置 MGF→mzML 转换共享缓存。
# 工作站实测纯 R 转换吞吐 0.19 MB/s，msfragger/sage 两槽对同一重编号 MGF
# 各转一遍，Macaca 规模（143 万谱）重复转换是 4–9 h 级浪费——转换产物缓存
# 到项目级 mzml_cache/ 目录，第二槽起直接复制；产物已存在且新于源 MGF 时
# 0 成本跳过（同 #18 标记跳过口径）。缓存键＝源 MGF 文件名+字节数（重编号
# 对同一输入确定性生成，跨引擎目录同名同内容）；并发双写以 临时文件+rename
# 原子落盘，后写胜出、内容一致无碍正确性。msconvert 外部路径不走本缓存。
# (workstation receipt priority): shared cache for the built-in
# MGF->mzML converter. The pure-R converter runs at 0.19 MB/s (measured);
# msfragger and sage each converting the same renumbered MGF wastes 4-9 h at
# Macaca scale. The converted mzML is cached at project level (mzml_cache/);
# later slots copy instead of re-converting, and an output newer than its
# source skips at zero cost (#18-style marker skip). Key = source MGF file
# name + byte size (renumbering is deterministic for the same input, so the
# per-engine copies share content); concurrent writers go through
# temp-file + rename. The external msconvert path bypasses this cache.
mgf_to_mzml_cached <- function(rmgf, out_mzml) {
  if (file.exists(out_mzml) && file.mtime(out_mzml) > file.mtime(rmgf)) {
    cat("[mzML 转换] 产物已存在且新于源 MGF，跳过（", basename(out_mzml), "）\n",
        "[mzML convert] output exists and is newer than the source MGF; skipping (",
        basename(out_mzml), ")\n", sep = "")
    return(invisible(TRUE))
  }
  cache_dir <- file.path(dirname(dirname(out_mzml)), "mzml_cache")
  cf <- file.path(cache_dir,
                  paste0(sub("\\.mgf$", "", basename(rmgf)), "_",
                         format(file.size(rmgf), scientific = FALSE, trim = TRUE),
                         ".mzML"))
  if (file.exists(cf)) {
    dir.create(dirname(out_mzml), showWarnings = FALSE, recursive = TRUE)
    if (file.copy(cf, out_mzml, overwrite = TRUE)) {
      cat("[mzML 转换] 命中共享缓存，复制（", basename(cf), " → ",
          basename(out_mzml), "）\n",
          "[mzML convert] shared cache hit, copying\n", sep = "")
      return(invisible(TRUE))
    }
  }
  mgf_to_mzml(rmgf, out_mzml)
  try({
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    tmp <- file.path(cache_dir, paste0(".tmp_", Sys.getpid(), "_", basename(cf)))
    if (file.copy(out_mzml, tmp, overwrite = TRUE)) file.rename(tmp, cf)
  }, silent = TRUE)
  invisible(TRUE)
}
mgf_to_mzml_base <- function(mgf_file, out_mzml) {
  lines <- readLines(mgf_file, warn = FALSE)
  esc_xml <- function(s) {
    gsub('"', "&quot;", gsub(">", "&gt;", gsub("<", "&lt;",
         gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
  }
  specs <- list(); cur <- NULL
  for (ln in lines) {
    ls <- trimws(ln)
    if (ls == "BEGIN IONS") {
      cur <- list(mz = numeric(0), int = numeric(0), charge = "2",
                  rt = "0", pep = "0", title = "")
    } else if (ls == "END IONS") {
      if (!is.null(cur)) specs[[length(specs) + 1L]] <- cur
      cur <- NULL
    } else if (!is.null(cur)) {
      if (startsWith(ls, "TITLE=")) {
        cur$title <- sub("^TITLE=", "", ls)
      } else if (startsWith(ls, "CHARGE=")) {
        cur$charge <- sub("\\+.*$", "", sub("^CHARGE=", "", ls))
      } else if (startsWith(ls, "RTINSECONDS=")) {
        cur$rt <- sub("^RTINSECONDS=", "", ls)
      } else if (startsWith(ls, "PEPMASS=")) {
        cur$pep <- strsplit(sub("^PEPMASS=", "", ls), "[ \t]+")[[1]][1]
      } else if (grepl("^[0-9]", ls)) {
        p <- strsplit(ls, "[ \t]+")[[1]]
        v1 <- suppressWarnings(as.numeric(p[1]))
        if (length(p) >= 2 && !is.na(v1)) {
          cur$mz <- c(cur$mz, v1)
          cur$int <- c(cur$int, suppressWarnings(as.numeric(p[2])))
        }
      }
    }
  }
  n <- length(specs)
  head1 <- paste0(
    '<?xml version="1.0" encoding="utf-8"?>\n',
    '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" ',
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
    'xsi:schemaLocation="http://psi.hupo.org/ms/mzml ',
    'http://psi.hupo.org/schema/mzml/1.1.0/mzML1.1.0.xsd">\n',
    '<mzML version="1.1.0">\n',
    '<cvList count="2">',
    '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" ',
    'URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/',
    'controlledVocabulary/psi-ms.obo"/>',
    '<cv id="UO" fullName="Unit Ontology" ',
    'URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>',
    '</cvList>\n',
    '<fileDescription><fileContent/></fileDescription>\n',
    # 实测验收补齐batmass-io 等
    # reader 要求标准头部链——softwareList / instrumentConfigurationList（run
    # 挂 defaultInstrumentConfigurationRef）/ dataProcessingList；cvParam 均带
    # value 属性（空值也写 value=""，否则 reader 取 getValue() 时 NPE）
    '<softwareList count="1"><software id="isoPSV" version="1.0">',
    '<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="isoPSV"/>',
    '</software></softwareList>\n',
    '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">',
    '<cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value="unknown"/>',
    '<componentList count="3">',
    '<source order="1"><cvParam cvRef="MS" accession="MS:1000844" name="unknown ionization type" value=""/></source>',
    '<analyzer order="2"><cvParam cvRef="MS" accession="MS:1000443" name="mass analyzer type" value=""/>',
    '<cvParam cvRef="MS" accession="MS:1000594" name="low energy collision-induced dissociation" value=""/></analyzer>',
    '<detector order="3"><cvParam cvRef="MS" accession="MS:1000026" name="detector type" value=""/></detector>',
    '</componentList></instrumentConfiguration></instrumentConfigurationList>\n',
    '<dataProcessingList count="1"><dataProcessing id="DP1">',
    '<processingMethod order="0" softwareRef="isoPSV"/>',
    '</dataProcessing></dataProcessingList>\n',
    '<run id="run1" defaultInstrumentConfigurationRef="IC1">\n',
    '<spectrumList count="', n, '">\n')
  blocks <- character(n); offs <- integer(n)
  pos <- nchar(head1, type = "bytes")
  for (i in seq_len(n)) {
    sp <- specs[[i]]; np <- length(sp$mz)
    mzb <- base64_encode_raw(writeBin(sp$mz, raw(), size = 8L, endian = "little"))
    inb <- base64_encode_raw(writeBin(sp$int, raw(), size = 8L, endian = "little"))
    blk <- paste0(
      '<spectrum index="', i - 1L, '" id="scan=', i, '" defaultArrayLength="', np, '">\n',
      # 实测修复MGF 谱是二级谱——写 ms level=2 + MSn spectrum，
      # 而非 MS1（此前把全部谱标成 MS1，Sage/MSFragger 都找不到 MS2 谱）
      '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>\n',
      '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n',
      # P1-1MGF 峰本就 centroided——补逐谱
      # centroid spectrum 标记。缺它 Sage panic（contains profile data，exit
      # 1033）、MetaMorpheus 计 All MS2 Scans: 0；补上后 Sage 实测出 1,457 PSM
      '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>\n',
      '<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="',
      if (np > 0) format(sum(sp$int), scientific = FALSE, trim = TRUE) else "0", '"/>\n',
      '<scanList count="1"><scan>\n',
      '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="',
      esc_xml(sp$rt), '" unitName="second" unitAccession="UO:0000010"/>\n',
      '<userParam name="title" value="', esc_xml(sp$title), '" type="xsd:string"/>\n',
      '</scan></scanList>\n',
      '<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>\n',
      '<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="',
      esc_xml(sp$pep), '"/>\n',
      '<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="',
      esc_xml(sp$charge), '"/>\n',
      '</selectedIon></selectedIonList></precursor></precursorList>\n',
      '<binaryDataArrayList count="2">\n',
      '<binaryDataArray encodedLength="', nchar(mzb, type = "bytes"), '">',
      '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
      '<binary>', mzb, '</binary></binaryDataArray>\n',
      '<binaryDataArray encodedLength="', nchar(inb, type = "bytes"), '">',
      '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
      '<binary>', inb, '</binary></binaryDataArray>\n',
      '</binaryDataArrayList>\n</spectrum>\n')
    blocks[i] <- blk
    offs[i] <- pos
    pos <- pos + nchar(blk, type = "bytes")
  }
  tail1 <- '</spectrumList>\n</run>\n</mzML>\n'
  idx_block <- paste0('<indexList count="1"><index name="spectrum">',
                      paste(sprintf('<offset idRef="scan=%d" spotID="" scanTime="0">%d</offset>',
                                    seq_len(n), offs), collapse = ""),
                      '</index></indexList>\n</indexedmzML>\n')
  con <- file(out_mzml, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(c(head1, blocks, tail1, idx_block), con, sep = "")
  n
}
# v1.0base64 快速编码器——3 字节组整数查表替代
# 逐位矩阵（旧 base64_encode_raw 每谱 ~1.6 KB 字节要构造 1.28 万个 logical
# 再 colSums，是纯 R 转换 0.19 MB/s 的主根因之一）。输出与旧实现逐字节
# 一致（含 1/2 字节尾组与 "=" 补位；mzml_convert_vec_test 对拍随机 raw）。
# (perf round 3): fast base64 via 3-byte-group integer table lookup,
# replacing the per-bit matrix of base64_encode_raw (which materialized ~12.8k
# logicals per average spectrum). Byte-identical output incl. 1/2-byte tails
# and "=" padding; cross-checked against the old encoder on random inputs.
base64_encode_fast <- function(bytes) {
  nb <- length(bytes)
  if (nb == 0L) return("")
  v <- as.integer(bytes)
  ng <- nb %/% 3L
  if (ng > 0L) {
    i3 <- seq_len(ng) * 3L
    b1 <- v[i3 - 2L]; b2 <- v[i3 - 1L]; b3 <- v[i3]
    # rbind＋列主序展平＝逐组交错（c1₁c2₁c3₁c4₁c1₂…）——此前按段拼接错序
    out <- b64_alphabet[as.vector(rbind((b1 %/% 4L) + 1L,
                                        ((b1 %% 4L) * 16L + b2 %/% 16L) + 1L,
                                        ((b2 %% 16L) * 4L + b3 %/% 64L) + 1L,
                                        (b3 %% 64L) + 1L))]
  } else out <- character(0)
  r <- nb %% 3L
  if (r == 1L) {
    b1 <- v[nb]
    out <- c(out, b64_alphabet[(b1 %/% 4L) + 1L],
             b64_alphabet[((b1 %% 4L) * 16L) + 1L], "=", "=")
  } else if (r == 2L) {
    b1 <- v[nb - 1L]; b2 <- v[nb]
    out <- c(out, b64_alphabet[(b1 %/% 4L) + 1L],
             b64_alphabet[((b1 %% 4L) * 16L + b2 %/% 16L) + 1L],
             b64_alphabet[((b2 %% 16L) * 4L) + 1L], "=")
  }
  paste0(out, collapse = "")
}
# v1.0mgf_to_mzml 向量化版。解析段：整文件
# trimws 一次、BEGIN/END 向量定位（良构输入＝严格交替；否则回落 _base 旧
# 实现），参数行（TITLE=/CHARGE=/RTINSECONDS=/PEPMASS=）向量子提取、每块
# 取最后一次出现值（与旧循环覆盖语义一致）；峰行批量 strsplit＋首两 token
# 取 mz/int（旧循环同判：token<2 或 mz 非数值的行丢弃，多余 token 忽略），
# split 按谱分组。编码段：base64_encode_fast。输出与 _base 逐字节一致
# （mgf_to_mzml_test 结构断言＋真件 diff 双重对拍）。
# (perf round 3): vectorized mgf_to_mzml. Parsing: one trimws pass,
# vectorized BEGIN/END locating (well-formed input = strictly alternating;
# malformed input falls back to _base), vectorized parameter extraction with
# last-occurrence-wins per block (same as the old loop's overwrite semantics);
# peak lines are bulk-strsplit with first-two-token mz/int (same discards as
# the old loop), grouped per spectrum via split. Encoding via
# base64_encode_fast. Output is byte-identical to _base (structural
# assertions in mgf_to_mzml_test plus real-file diff).
mgf_to_mzml <- function(mgf_file, out_mzml) {
  lines <- readLines(mgf_file, warn = FALSE)
  ls <- trimws(lines)
  begins <- which(ls == "BEGIN IONS")
  ends <- which(ls == "END IONS")
  n <- length(begins)
  # 良构性：BEGIN/END 数相等且严格交替（b1<e1<b2<e2<…）；否则走旧实现
  well_formed <- length(ends) == n && n > 0L &&
    all(begins < ends) && (n == 1L || all(begins[-1L] > ends[-n]))
  if (!well_formed) return(mgf_to_mzml_base(mgf_file, out_mzml))
  esc_xml <- function(s) {
    gsub('"', "&quot;", gsub(">", "&gt;", gsub("<", "&lt;",
         gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
  }
  nl <- length(ls)
  posv <- seq_len(nl)
  blk_of <- findInterval(posv, begins)
  in_blk <- blk_of > 0L & posv < ends[pmax(blk_of, 1L)]
  # 参数行（块内、匹配前缀）；每块取最后一次出现（旧循环覆盖语义）
  pick_last <- function(flag, subpat, default) {
    hit <- which(flag & in_blk)
    if (length(hit) == 0L) return(rep(default, n))
    b <- blk_of[hit]
    keep <- !duplicated(b, fromLast = TRUE)   # 每块保留最后一次
    vals <- sub(subpat, "", ls[hit[keep]])
    out <- rep(default, n)
    out[b[keep]] <- vals
    out
  }
  title_v <- pick_last(startsWith(ls, "TITLE="), "^TITLE=", "")
  charge_v <- pick_last(startsWith(ls, "CHARGE="), "^CHARGE=", "2")
  charge_v <- sub("\\+.*$", "", charge_v)
  rt_v <- pick_last(startsWith(ls, "RTINSECONDS="), "^RTINSECONDS=", "0")
  pep_v <- pick_last(startsWith(ls, "PEPMASS="), "^PEPMASS=", "0")
  pep_v <- vapply(strsplit(pep_v, "[ \t]+"), `[`, character(1), 1L)
  # 峰行：块内、首字符为数字（旧循环仅对块内行做 grepl("^[0-9]")）
  is_peak <- grepl("^[0-9]", ls) & in_blk
  pk <- which(is_peak)
  npks <- integer(n)
  mz_list <- int_list <- vector("list", n)
  if (length(pk) > 0L) {
    toks <- strsplit(ls[pk], "[ \t]+")
    tl <- lengths(toks)
    st <- cumsum(tl) - tl + 1L
    flat <- unlist(toks, use.names = FALSE)
    mz_t <- suppressWarnings(as.numeric(flat[st]))
    int_t <- suppressWarnings(as.numeric(flat[st + 1L]))
    ok <- !is.na(mz_t) & tl >= 2L
    if (any(ok)) {
      g <- factor(blk_of[pk[ok]], levels = seq_len(n))
      mz_list <- split(mz_t[ok], g)
      int_list <- split(int_t[ok], g)
      npks <- tabulate(blk_of[pk[ok]], nbins = n)
    }
  }
  tic_v <- vapply(seq_len(n), function(i) {
    x <- int_list[[i]]
    if (length(x)) sum(x) else 0
  }, numeric(1))
  # 逐值格式化＝与旧标量 format 调用完全同义（向量 format 的对齐规则不引入）
  tic_fmt <- vapply(tic_v, function(x) format(x, scientific = FALSE, trim = TRUE), "")
  title_e <- esc_xml(title_v); rt_e <- esc_xml(rt_v)
  pep_e <- esc_xml(pep_v); chg_e <- esc_xml(charge_v)
  head1 <- paste0(
    '<?xml version="1.0" encoding="utf-8"?>\n',
    '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" ',
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
    'xsi:schemaLocation="http://psi.hupo.org/ms/mzml ',
    'http://psi.hupo.org/schema/mzml/1.1.0/mzML1.1.0.xsd">\n',
    '<mzML version="1.1.0">\n',
    '<cvList count="2">',
    '<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" ',
    'URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/',
    'controlledVocabulary/psi-ms.obo"/>',
    '<cv id="UO" fullName="Unit Ontology" ',
    'URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>',
    '</cvList>\n',
    '<fileDescription><fileContent/></fileDescription>\n',
    '<softwareList count="1"><software id="isoPSV" version="1.0">',
    '<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="isoPSV"/>',
    '</software></softwareList>\n',
    '<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">',
    '<cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value="unknown"/>',
    '<componentList count="3">',
    '<source order="1"><cvParam cvRef="MS" accession="MS:1000844" name="unknown ionization type" value=""/></source>',
    '<analyzer order="2"><cvParam cvRef="MS" accession="MS:1000443" name="mass analyzer type" value=""/>',
    '<cvParam cvRef="MS" accession="MS:1000594" name="low energy collision-induced dissociation" value=""/></analyzer>',
    '<detector order="3"><cvParam cvRef="MS" accession="MS:1000026" name="detector type" value=""/></detector>',
    '</componentList></instrumentConfiguration></instrumentConfigurationList>\n',
    '<dataProcessingList count="1"><dataProcessing id="DP1">',
    '<processingMethod order="0" softwareRef="isoPSV"/>',
    '</dataProcessing></dataProcessingList>\n',
    '<run id="run1" defaultInstrumentConfigurationRef="IC1">\n',
    '<spectrumList count="', n, '">\n')
  blocks <- character(n); offs <- integer(n)
  pos <- nchar(head1, type = "bytes")
  for (i in seq_len(n)) {
    mz <- mz_list[[i]]; inten <- int_list[[i]]
    mzb <- base64_encode_fast(writeBin(mz, raw(), size = 8L, endian = "little"))
    inb <- base64_encode_fast(writeBin(inten, raw(), size = 8L, endian = "little"))
    blk <- paste0(
      '<spectrum index="', i - 1L, '" id="scan=', i, '" defaultArrayLength="', npks[i], '">\n',
      '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>\n',
      '<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n',
      '<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>\n',
      '<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="',
      tic_fmt[i], '"/>\n',
      '<scanList count="1"><scan>\n',
      '<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="',
      rt_e[i], '" unitName="second" unitAccession="UO:0000010"/>\n',
      '<userParam name="title" value="', title_e[i], '" type="xsd:string"/>\n',
      '</scan></scanList>\n',
      '<precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>\n',
      '<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="',
      pep_e[i], '"/>\n',
      '<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="',
      chg_e[i], '"/>\n',
      '</selectedIon></selectedIonList></precursor></precursorList>\n',
      '<binaryDataArrayList count="2">\n',
      '<binaryDataArray encodedLength="', nchar(mzb, type = "bytes"), '">',
      '<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
      '<binary>', mzb, '</binary></binaryDataArray>\n',
      '<binaryDataArray encodedLength="', nchar(inb, type = "bytes"), '">',
      '<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>',
      '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>',
      '<binary>', inb, '</binary></binaryDataArray>\n',
      '</binaryDataArrayList>\n</spectrum>\n')
    blocks[i] <- blk
    offs[i] <- pos
    pos <- pos + nchar(blk, type = "bytes")
  }
  tail1 <- '</spectrumList>\n</run>\n</mzML>\n'
  idx_block <- paste0('<indexList count="1"><index name="spectrum">',
                      paste(sprintf('<offset idRef="scan=%d" spotID="" scanTime="0">%d</offset>',
                                    seq_len(n), offs), collapse = ""),
                      '</index></indexList>\n</indexedmzML>\n')
  con <- file(out_mzml, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(c(head1, blocks, tail1, idx_block), con, sep = "")
  n
}
#' 经 systeminfo 获取可用物理内存（MB）。
#' Available physical memory in MB via systeminfo.
#'
#' R ≥4.3 不再支持 memory.size()/memory.limit()；systeminfo 输出可能为系统
#' ANSI 代码页（此处 GBK），先转 UTF-8 再解析；失败退回 PowerShell 数值查询，
#' 再失败返回 NA_real_。
#' memory.size()/memory.limit() are no longer supported on R >= 4.3; the
#' systeminfo output may arrive in the OS ANSI code page (GBK here) and is
#' converted to UTF-8 before parsing; falls back to a PowerShell numeric query,
#' then to NA_real_ on failure.
#'
#' @return
#'   返回可用内存 MB 数值（失败为 NA_real_）。
#'   The available memory in MB (NA_real_ on failure).
#' @keywords internal
system_avail_mb <- function() {
  if (.Platform$OS.type != "windows") return(NA_real_)
  out <- tryCatch(
    system("systeminfo", intern = TRUE, show.output.on.console = FALSE),
    error = function(e) character(0))
  mem_lines <- out[grepl("物理内存|Physical Memory", out)]
  if (length(mem_lines) == 0) {
    out <- tryCatch(iconv(out, from = "GBK", to = "UTF-8", sub = ""),
                    error = function(e) out)
    mem_lines <- out[grepl("物理内存|Physical Memory", out)]
  }
  if (length(mem_lines) > 0) cat(mem_lines, sep = "\n")
  line <- mem_lines[grepl("可用|Available", mem_lines)][1]
  if (!is.na(line) && nzchar(line)) {
    num <- regmatches(line, regexpr("[0-9,]+", line))
    if (length(num) > 0) return(as.numeric(gsub(",", "", num)))
  }
  kb <- tryCatch(
    system(paste0("powershell -NoProfile -Command \"",
                  "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory\""),
           intern = TRUE),
    error = function(e) character(0))
  v <- suppressWarnings(as.numeric(trimws(kb)))
  if (length(v) > 0 && !is.na(v[1]) && v[1] > 0) return(v[1] / 1024)
  NA_real_
}
#' 物理内存总量（MB；PowerShell CIM 优先、systeminfo 兜底）。
#' Total physical memory (MB; PowerShell CIM first, systeminfo fallback).
#'
#' MSFragger -Xmx 自适应需要总量。systeminfo 在
#' Windows 上是公认慢命令（10-30 秒），因此优先 PowerShell CIM（秒级）。
#' (P1-5/N-2): the adaptive MSFragger -Xmx needs the total.
#' systeminfo is notoriously slow on Windows (10-30 s), so PowerShell CIM
#' runs first (sub-second).
system_total_mb <- local({
  cached <- NULL
  fail_at <- NULL
  function() {
    # 进程级记忆化iso_mem_guard 逐万行调用，探测只做一次。
    # 失败值不永久缓存——首探若因瞬时原因失败（权限/瞬断），
    # 30 分钟后允许重试；期间返回 NA 由调用方走默认值。
    # Cache only successes: a transient first-probe failure (permissions, hiccup)
    # would otherwise be frozen for the whole run. Retry after 30 min; in between
    # return NA and let callers use their defaults.
    if (!is.null(cached)) return(cached)
    if (!is.null(fail_at) && as.numeric(Sys.time() - fail_at, units = "mins") < 30)
      return(NA_real_)
    if (.Platform$OS.type != "windows") return(NA_real_)
    kb <- tryCatch(
      system(paste0("powershell -NoProfile -Command \"",
                    "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory\""),
             intern = TRUE),
      error = function(e) character(0))
    v <- suppressWarnings(as.numeric(trimws(kb)))
    if (length(v) > 0 && !is.na(v[1]) && v[1] > 0) {
      cached <<- v[1] / 1024^2
      return(cached)
    }
    out <- tryCatch(
      system("systeminfo", intern = TRUE, show.output.on.console = FALSE),
      error = function(e) character(0))
    mem_lines <- out[grepl("物理内存|Physical Memory", out)]
    if (length(mem_lines) == 0) {
      out <- tryCatch(iconv(out, from = "GBK", to = "UTF-8", sub = ""),
                      error = function(e) out)
      mem_lines <- out[grepl("物理内存|Physical Memory", out)]
    }
    line <- mem_lines[grepl("总量|Total", mem_lines)][1]
    if (!is.na(line) && nzchar(line)) {
      num <- regmatches(line, regexpr("[0-9,]+", line))
      if (length(num) > 0) {
        cached <<- as.numeric(gsub(",", "", num))
        return(cached)
      }
    }
    fail_at <<- Sys.time()
    NA_real_
  }
})
#' system2 安全封装：参数统一 shQuote。
#' Safe system2 wrapper: every argument is shQuote'd .
#'
#' 旧写法大量未加引号的路径参数在中文/含空格路径下被 shell 拆分而静默失败；
#' 已用 shQuote 的调用点不要经本封装（会二次引用）。标志类参数（如 -Xmx15G）
#' 无特殊字符，shQuote 原样返回，安全。
#' Old call sites passed many path arguments unquoted, which the shell splits
#' on Chinese/space-containing paths and fails silently; call sites that
#' already shQuote must NOT go through this wrapper (double quoting). Flag
#' arguments such as -Xmx15G contain no special characters and pass through
#' shQuote unchanged.
iso_system2 <- function(command, args = character(), ...) {
  system2(command, args = vapply(args, shQuote, character(1), USE.NAMES = FALSE), ...)
}
#' 步骤 1 资源评估：核数与内存 → 读取核数。
#' Step-1 resource assessment: cores and memory -> reader cores.
#'
#' 硬编码简单算术：单次并行读取的峰值内存约为最大 MGF 文件的 6 倍；
#' 上限为 CPU 并行数减一（留一核给系统），内存上限为可用内存可容纳的
#' 峰值数（向下取整），不小于 1。
#' Simple hardcoded arithmetic: the peak memory of one parallel read is about
#' 6x the largest MGF file; the upper bound is the CPU parallel count minus one
#' (leave one core for the system) and the memory bound is how many such peaks
#' fit in the available memory (floor); never less than 1.
#'
#' @param mgf_files
#'   字符向量，输入 MGF 文件路径。
#'   Character vector, the input MGF file paths.
#' @return
#'   返回 list(max_cores=整数, avail_mb=数值, ...)。
#'   A list(max_cores = integer, avail_mb = numeric, ...).
#' @keywords internal
assess_resources <- function(mgf_files) {
  cores_cpu <- max(parallel::detectCores() - 1L, 1L)
  avail_mb <- system_avail_mb()
  if (is.na(avail_mb) || avail_mb <= 0) {
    avail_mb <- 8 * 1024
    cat("available memory unknown, fallback", avail_mb, "MB\n")
  }
  biggest_mb <- max(file.size(mgf_files)) / 1024^2
  est_peak_mb <- 6 * biggest_mb
  # 内存余量 0.7：并行峰值吃满可用内存会导致 worker 被杀（socket 连接错误），
  # 因此最多只允许占用约 70% 的可用内存
  # 0.7 headroom: letting the parallel peak fill all available memory kills the
  # workers (socket connection errors), so use at most ~70% of the free memory
  cores_mem <- max(1L, as.integer(floor(avail_mb * 0.7 / est_peak_mb)))
  max_cores <- max(1L, min(cores_cpu, cores_mem))
  cat("resource assessment: cpu cores =", cores_cpu, "| avail mem (MB) ~=",
      round(avail_mb), "| biggest MGF (MB) =", round(biggest_mb, 1),
      "| est peak per read (MB) ~=", round(est_peak_mb),
      "-> max_cores =", max_cores, "\n")
  max_cores
}
#' 每输入 MGF 一个无标记 MGF（手工步骤 2 的库实现）。
#' One unlabeled MGF per input MGF (library form of manual step 2).
#'
#' 严格子集：保留标题在常规搜库结果中且 Label 首字符=="1" 的谱块；缺失标题
#' 以同 raw＋scan 谱块补齐（优先同价态），TITLE 改写为 pFind 报告标题、MS2 峰
#' 原样（isoPSV 对 MS1 无严格要求）；同标题重复块只留首次；空产出文件的旧结果删除。
#' Strict subset: keep the blocks whose title is in the regular search result
#' with Label first char == "1"; missing titles are patched with the same
#' raw + scan block (preferring the same charge), TITLE rewritten to the pFind
#' report title, MS2 peaks untouched; duplicate-titled blocks keep only the
#' first; stale outputs of files that kept nothing are removed.
#'
#' @param psm_file
#'   字符串，常规搜库结果（.spectra）。
#'   Character string, the regular search result (.spectra).
#' @param mgf_files
#'   字符向量，输入 MGF 文件。
#'   Character vector, the input MGF files.
#' @param out_dir
#'   字符串，输出目录。
#'   Character string, the output directory.
#' @param max_cores
#'   整数，并行核数，默认为 1。
#'   Integer, worker count; defaults to 1.
#' @param cl
#'   PSOCK 集群对象或 NULL（NULL 时自建）。
#'   A PSOCK cluster object or NULL (built when NULL).
#' @return
#'   返回写出的无标记 MGF 文件路径。
#'   The written unlabeled MGF file paths.
#' @keywords internal
#' 第一步 psm_file 格式校验接受任何带 File_Name
#' 列的制表符结果表——pFind 家族表（File_Name+Label）、适配器交换格式
#' （.spectra）与其他搜库软件的表格输出均可用；非表格格式（XML 等）须先经
#' 对应适配器转交换格式。无标记/标记谱区分：有 Label 列按 pFind 语义（首字符
#' 1＝无标记）；无 Label 列时无标记样本全部视为无标记谱，标记样本需显式设
#' psm_all_unlabeled=1 确认全为无标记谱，否则报错。返回 list(i_fn, i_lab,
#' all_unlabeled)。
#' Step-1 psm_file format check any
#' tab-separated result table with a File_Name column is accepted - pFind
#' family tables (File_Name + Label), adapter exchange .spectra and other
#' engines' tabular outputs; non-tabular formats (XML etc.) must first pass
#' their adapter. Label semantics: with a Label column the pFind rule applies
#' (first char 1 = unlabeled); without it unlabeled samples treat all spectra
#' as unlabeled, labeled samples must confirm via psm_all_unlabeled=1.
check_psm_header <- function(psm_file, label_type = "unlabeled",
                             psm_all_unlabeled = FALSE) {
  first <- readLines(psm_file, n = 1, warn = FALSE)
  is_xml <- length(first) > 0 && grepl("^\\s*<", first)
  if (is_xml) {
    stop("[第一步] psm_file 不能是非表格格式（XML 等）：", psm_file,
         "\n请先用对应引擎适配器转换为交换格式（或直接填该引擎的表格输出）。\n",
         "[step 1] psm_file must not be a non-tabular format (XML etc.): ", psm_file,
         "\nconvert it with the corresponding adapter first (or fill that engine's tabular output).")
  }
  hdr <- strsplit(first, "\t")[[1]]
  i_fn <- which(hdr == "File_Name")
  i_lab <- which(hdr == "Label")
  if (length(i_fn) != 1) {
    stop("[第一步] psm_file 缺 File_Name 列：", psm_file,
         "\n第一步接受任何带 File_Name 列的制表符结果表（pFind 家族表、交换格式、",
         "各引擎表格输出）；非表格格式（XML 等）请先经对应适配器转换。\n",
         "[step 1] psm_file lacks the File_Name column: ", psm_file,
         "\nstep 1 accepts any tab-separated result table with a File_Name column",
         " (pFind tables, the exchange format, engines' tabular outputs); non-tabular",
         " formats must first go through their adapter.")
  }
  if (length(i_lab) != 1) {
    # 无 Label 列：无标记样本全收；标记样本须显式确认全为无标记谱
    if (identical(label_type, "unlabeled") || isTRUE(psm_all_unlabeled)) {
      return(list(i_fn = i_fn, i_lab = NA_integer_, all_unlabeled = TRUE))
    }
    stop("[第一步] psm_file 无 Label 列，无法区分无标记/标记谱：", psm_file,
         "\n标记样本请使用 pFind 家族表（File_Name+Label）；如该结果表对应的谱全部",
         "为无标记，请在 config_step1 设 psm_all_unlabeled\t1 显式确认。\n",
         "[step 1] psm_file has no Label column, cannot distinguish unlabeled/labeled: ",
         psm_file, "\nlabeled samples need a pFind-family table (File_Name + Label);",
         " if all spectra behind this table are unlabeled, set psm_all_unlabeled 1 in",
         " config_step1 to confirm.")
  }
  list(i_fn = i_fn, i_lab = i_lab, all_unlabeled = FALSE)
}
#' 第一步 psm_file 格式嗅探判断是否需要先经适配器转换。
#' 返回 ""＝可直接当表读（非 XML 且含制表符）；"xtandem"（<bioml）；"msfragger"
#' （<msms_pipeline_analysis）；其余未知返回 "unknown"。
psm_sniff_engine <- function(psm_file, first_line = NULL) {
  if (is.null(first_line)) {
    first_line <- readLines(psm_file, n = 3, warn = FALSE)
  }
  if (!length(first_line)) return("unknown")
  # P2真实 XML 首行常是 <?xml …?>
  # 声明（<bioml> 等在下一行）——只看首行恒嗅探失败；取前 3 行、剥声明、
  # 再看首个标签行
  lines <- sub("^\\s*<\\?xml[^>]*\\?>", "", as.character(first_line))
  lines <- lines[grepl("^\\s*<", lines)]
  if (length(lines)) {
    l1 <- lines[1]
    if (grepl("^\\s*<bioml", l1)) return("xtandem")
    if (grepl("^\\s*<msms_pipeline_analysis", l1)) return("msfragger")
    if (grepl("^\\s*<mzIdentML", l1)) return("mzidentml")
    return("unknown")
  }
  if (any(grepl("\t", as.character(first_line)))) return("")
  "unknown"
}
#' 第一步 psm_file 自动转换非表格格式（XML/CSV）可以
#' 转换——bioml 自动走 X!Tandem 适配器；CSV/其他格式经 config_step1 的
#' psm_type 键指定适配器（xtandem/mascot/peaks），转换到交换格式后作为第一步
#' 输入；msfragger/pnovo 因需要重编号映射/参数文件（第一步不产出），给出定向
#' 指引（走其表格输出或 config_step2 search_file 槽位）。返回可用的表格路径。
psm_prepare_input <- function(psm_file, cfg, code_dir) {
  first_line <- readLines(psm_file, n = 3, warn = FALSE)
  if (!length(first_line)) return(psm_file)
  # 用户检测反馈：网页 STEP1 的 search_engine 键此前 R 侧零消费——
  # MaxQuant msms.txt 这类表格嗅探不出引擎，直通到 File_Name 校验才报错。
  # 现按用户所选引擎指认：maxquant 走两参适配器自动转换（无需映射）；
  # sage/comet/msfragger 的结果需扫描号映射（第一步不产出），给出双路指引。
  # The STEP1 search_engine key was never consumed R-side: tabular results
  # like MaxQuant msms.txt sniff to no engine and only died at the
  # File_Name check. Now the selected engine names it: maxquant converts
  # via its two-arg adapter (no map needed); sage/comet/msfragger results
  # need the scan map (not produced at step 1) and get a two-path guidance.
  se_key <- if (!is.null(cfg$search_engine) && length(cfg$search_engine) > 0)
    tolower(trimws(tail(cfg$search_engine, 1))) else ""
  engine <- psm_sniff_engine(psm_file, first_line)
  if (!nzchar(engine)) {
    has_fn <- any(grepl("File_Name", first_line, fixed = TRUE))
    if (!has_fn && se_key == "maxquant") {
      engine <- "maxquant"
    } else if (!has_fn && se_key %in% c("sage", "comet", "msfragger")) {
      stop("[第一步] search_engine=", se_key, " 的结果文件需要扫描号映射才能转换（第一步不产出映射）：",
           psm_file,
           "\n两条路径任选其一：",
           "\n  (1) 把该引擎结果放入 config_step2 的 search_file 槽位接入（自动链路含映射与适配）；",
           "\n  (2) 改填任何带 File_Name 列的制表符结果表（pFind 家族表或交换格式）。",
           "\n[step 1] search_engine=", se_key, " results need the scan map to convert, ",
           "which step 1 does not produce: ", psm_file,
           "\nEither (1) feed the engine result through a config_step2 search_file slot ",
           "(the automated path provides the map and adapter), or (2) fill any ",
           "tab-separated table with a File_Name column (pFind-family or interchange).")
    }
  }
  if (!nzchar(engine)) return(psm_file)              # 表格：直接可用（File_Name 校验在后）
  psm_type <- if (!is.null(cfg$psm_type) && length(cfg$psm_type) > 0)
    tolower(trimws(tail(cfg$psm_type, 1))) else ""
  if (engine == "unknown") {
    if (psm_type %in% c("xtandem", "mascot", "peaks")) engine <- psm_type
    else if (psm_type %in% c("msfragger", "pnovo", "mzidentml")) engine <- psm_type
    else {
      stop("[第一步] 无法自动识别 psm_file 格式：", psm_file,
           "\n可自动转换：X!Tandem XML（<bioml>，自动识别）与 Mascot/PEAKS CSV",
           "（请在 config_step1 设 psm_type 为 mascot/peaks）；",
           "其余格式请先经对应适配器转换或走 config_step2 search_file 槽位。\n",
           "[step 1] cannot auto-detect psm_file format: ", psm_file,
           "\nauto-convertible: X!Tandem XML (<bioml>, auto) and Mascot/PEAKS CSV",
           " (set psm_type mascot/peaks in config_step1); other formats must first",
           " pass their adapter or go into config_step2 search_file.")
    }
  }
  if (engine %in% c("msfragger", "pnovo", "mzidentml")) {
    stop("[第一步] ", engine, " 结果的自动转换需要重编号映射/参数文件（第一步不产出）：",
         psm_file, "\n请填该引擎的表格输出，或经 config_step2 的 search_file 槽位接入。\n",
         "[step 1] auto-converting ", engine, " output needs the renumber map/params file",
         " (not produced at step 1): ", psm_file, "\nfill that engine's tabular output,",
         " or put it into config_step2 search_file.")
  }
  adapter <- switch(engine, xtandem = "search_adapter_xtandem.R",
                    mascot = "search_adapter_mascot.R",
                    peaks = "search_adapter_peaks.R",
                    maxquant = "search_adapter_maxquant.R")
  out_ex <- file.path(cfg$config_dir, "unlabeled_tmp", "psm_exchange.spectra")
  dir.create(dirname(out_ex), recursive = TRUE, showWarnings = FALSE)
  logf <- paste0(out_ex, ".log")
  # 空格路径修复裸 system2 参数不自动加引号，
  # 代码目录/输入/输出路径含空格时被 cmd 拆断（适配器 exit 5）——三个参数
  # 全部 shQuote（裸 system2 的既有安全模式）
  # P2xtandem 转换与 step2 同口径——同目录
  # <同名>_renum_map.tsv（new_scan<TAB>original_title）存在时一并传入，适配器按
  # <group id> 解析谱标题（否则 label=蛋白描述，谱级信息退化）
  map_xt <- if (engine == "xtandem")
    sub("\\.xml$", "_renum_map.tsv", psm_file, ignore.case = TRUE) else ""
  args <- c(shQuote(file.path(code_dir, adapter)), shQuote(psm_file),
            shQuote(out_ex))
  if (nzchar(map_xt) && file.exists(map_xt)) args <- c(args, shQuote(map_xt))
  st <- system2(file.path(R.home("bin"), "Rscript.exe"), args,
                stdout = logf, stderr = logf)
  if (st != 0 || !file.exists(out_ex)) {
    stop("[第一步] psm_file 自动转换失败（", engine, " 适配器 exit ", st, "）：",
         psm_file, "（见 ", logf, "）\n",
         "[step 1] psm_file auto-conversion failed (", engine, " adapter exit ", st,
         "): ", psm_file, " (see ", logf, ")")
  }
  cat("[第一步] psm_file 已自动转换（", engine, " → 交换格式）：", out_ex, "\n",
      sep = "")
  out_ex
}
extract_unlabeled_per_mgf <- function(psm_file, mgf_files, out_dir, max_cores = 1, cl = NULL,
                                      label_type = "unlabeled",
                                      psm_all_unlabeled = FALSE,
                                      keep_all = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  own_cl <- is.null(cl)
  if (own_cl) {
    cl <- parallel::makeCluster(max_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ispsv_bootstrap_cluster(cl)
  }
  idx <- check_psm_header(psm_file, label_type, psm_all_unlabeled)
  hdr <- strsplit(readLines(psm_file, n = 1), "\t")[[1]]
  i_fn <- idx$i_fn; i_lab <- idx$i_lab
  cc <- rep("NULL", length(hdr))
  cc[i_fn] <- "character"
  if (!is.na(i_lab)) cc[i_lab] <- "character"
  #：读 Exp.MH./Charge（有则读）——回退选择的 ppm 比对用
  i_mh <- match("Exp.MH.", hdr); i_z <- match("Charge", hdr)
  if (!is.na(i_mh)) cc[i_mh] <- "character"
  if (!is.na(i_z)) cc[i_z] <- "character"
  ps <- read.delim(psm_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   colClasses = cc, na.strings = character(0))
  if (isTRUE(keep_all)) {
    wanted <- unique(as.character(ps$File_Name))
    cat("unlabeled_channel=0：输入结果文件的全部鉴定谱都提取：", length(wanted), " 张\n")
  } else if (isTRUE(idx$all_unlabeled)) {
    wanted <- unique(as.character(ps$File_Name))
    cat("psm_file 无 Label 列，全部谱按无标记处理：", length(wanted), " 张\n")
  } else {
    wanted <- unique(ps$File_Name[substr(ps$Label, 1, 1) == "1"])
  }
  cat("wanted titles:", length(wanted), "\n")

  wanted_stems <- unique(vapply(wanted, function(t) title_raw_scan(t)["raw"],
                                character(1)))
  per_file_kept <- list()          # basename -> list of MGF objects
  per_file_input <- list()         # basename -> input block count
  #：两个 scan 键改登记"全部候选"（file/title/charge/mz），
  # 回退时按 ppm 最小＋价态一致选优（pick_best_src）
  scan_c_map <- new.env(hash = TRUE)  # 'RAW.SCAN.CHARGE' -> list of entries
  scan_map <- new.env(hash = TRUE)    # 'RAW.SCAN' -> list of entries
  mass_map <- new.env(hash = TRUE)    # title -> Exp.MH.（ppm 比对用）
  if (!is.na(i_mh)) {
    mh_vals <- suppressWarnings(as.numeric(ps[["Exp.MH."]]))
    fn_vals <- as.character(ps$File_Name)
    for (j in seq_len(nrow(ps)))
      if (!is.na(mh_vals[j]) && is.null(mass_map[[fn_vals[j]]]))
        mass_map[[fn_vals[j]]] <- mh_vals[j]
  }
  seen_env <- new.env(hash = TRUE)    # kept titles (dedupe duplicate blocks)
  missing <- wanted

  stems <- vapply(mgf_files, function(f) sub("_HCDFT$", "",
                    sub("\\.mgf$", "", basename(f))), character(1))
  need_idx <- which(stems %in% wanted_stems)
  # 一次并行读取全部所需文件（每文件一个 worker；墙钟时间＝最慢文件）——后续循环只从内存选择/补齐
  # load ALL needed files in one parallel pass (one worker per file; wall
  # time = slowest file) - the loop below only selects/patches from memory
  # （性能第二轮）：补丁分支的 src 文件必在 need_idx 内——loaded 按
  # basename 命名保留到补丁阶段结束，二次整读 MGF 去掉（原实现同文件读两遍）
  loaded <- load_mgf_shared(mgf_files[need_idx], cl)
  names(loaded) <- basename(mgf_files[need_idx])
  for (k0 in seq_along(need_idx)) {
    f <- mgf_files[need_idx[k0]]
    m <- loaded[[k0]]
    pars <- extractParameters_MGF(m)
    sel <- match(pars$Title, wanted, nomatch = 0L) > 0L
    if (any(sel)) {
      tt <- pars$Title[sel]
      first <- vapply(tt, function(t) {
        if (exists(t, envir = seen_env, inherits = FALSE)) return(FALSE)
        assign(t, TRUE, envir = seen_env)
        TRUE
      }, logical(1))
      sel[which(sel)] <- first
      if (any(sel)) {
        per_file_kept[[basename(f)]] <- m[sel]
        missing <- setdiff(missing, pars$Title[sel])
      }
    }
    keys <- vapply(pars$Title, title_scan_key, character(1))
    kc_vec <- paste0(keys, ".", pars$Charge)
    pm_vec <- if (!is.null(pars$pepMass)) pars$pepMass else rep(NA_real_, length(pars$Title))
    for (j in seq_along(pars$Title)) {
      ent <- c(file = basename(f), title = pars$Title[j],
               charge = as.character(pars$Charge[j]),
               mz = as.character(pm_vec[j]))
      kc <- kc_vec[j]
      prev <- scan_c_map[[kc]]
      scan_c_map[[kc]] <- if (is.null(prev)) list(ent) else c(prev, list(ent))
      kk <- keys[j]
      prev2 <- scan_map[[kk]]
      scan_map[[kk]] <- if (is.null(prev2)) list(ent) else c(prev2, list(ent))
    }
    per_file_input[[basename(f)]] <- length(m)
    cat(sprintf("  %-55s in=%d kept=%d\n", basename(f), length(m), sum(sel)))
    rm(m, pars)
    invisible(gc())
  }
  cat("missing after exact match:", length(missing), "\n")

  patched <- character(0)
  if (length(missing) > 0) {
    need <- list()   # file -> named vector c(missing_title = source_title)
    for (t in missing) {
      k <- title_scan_key(t)
      chg <- strsplit(t, ".", fixed = TRUE)[[1]][4]
      cands <- scan_c_map[[paste0(k, ".", chg)]]
      if (is.null(cands)) cands <- scan_map[[k]]
      if (!is.null(cands)) {
        src <- pick_best_src(cands, chg, mass_map[[t]])
        f0 <- unname(src["file"])
        if (is.null(need[[f0]])) need[[f0]] <- character(0)
        need[[f0]][t] <- unname(src["title"])
      }
    }
    # （性能第二轮）：src_files 必在 need_idx 内——直接按 basename 从
    # 已加载的 loaded 取（不再二次整读同一批 MGF）
    src_load <- loaded[names(need)]
    names(src_load) <- names(need)
    for (f0 in names(need)) {
      m <- src_load[[f0]]
      pars <- extractParameters_MGF(m)
      for (t in names(need[[f0]])) {
        idx <- match(need[[f0]][t], pars$Title)
        if (is.na(idx)) next
        obj <- m[[idx]]
        obj@title <- t
        per_file_kept[[f0]] <- c(per_file_kept[[f0]], list(obj))
        patched <- c(patched, t)
      }
      rm(m, pars)
      invisible(gc())
    }
  }
  cat("patched (same raw+scan):", length(patched), "\n")
  rm(loaded)
  invisible(gc())

  final_missing <- setdiff(missing, patched)
  for (f0 in names(per_file_kept)) {
    writeMGFFile(per_file_kept[[f0]], file.path(out_dir, f0))
  }
  # 删除本次运行中什么都没保留的文件的旧产物
  # remove stale outputs of files that kept nothing in this run
  for (f0 in setdiff(basename(mgf_files), names(per_file_kept))) {
    fp <- file.path(out_dir, f0)
    if (file.exists(fp)) unlink(fp)
  }
  if (length(final_missing) > 0) {
    writeLines(sort(final_missing), file.path(out_dir, "missing_titles.txt"))
  }
  summ <- data.frame(
    source = names(per_file_input),
    input_spectra = unlist(per_file_input, use.names = FALSE),
    stringsAsFactors = FALSE)
  kept_all <- setNames(rep(0L, nrow(summ)), summ$source)
  for (f0 in names(per_file_kept)) kept_all[f0] <- length(per_file_kept[[f0]])
  summ$kept <- unname(kept_all)
  cat("total wanted:", length(wanted), "| kept:", sum(summ$kept),
      "| patched:", length(patched), "| missing:", length(final_missing), "\n")
  summ
}
#' 合并前对搜库结果表应用生产修饰校正（1.modification_exchange.R 的忠实移植）。
#' Apply the production modification corrections to a search-result table before merging (faithful port of 1.modification_exchange.R).
#'
#' 四种校正：① Modification 含 "_13C" 的行删除；② "0,Met-" 行：删首字母、
#' 修饰位点减 1、加 "0,Acetyl[AnyN-term];"（位点减到 0 的行删除）；③ Modification
#' 恰为 "0,Arg[AnyN-term];" 的行：序列前补 R、清修饰；④ "N,XXX->YYY[.]" 单残基替换：
#' 该位字母换成目标残基单字母码、替换条目移除。含 pFind_* 重复列的文件（pUniFind_rescore）
#' 跳过校正（生产行为）。
#' Four corrections: 1. rows whose Modification contains "_13C" are removed;
#' 2. "0,Met-" rows drop the first letter, decrement positions by 1 and add
#' "0,Acetyl[AnyN-term];" (rows ending at position 0 are removed); 3. rows whose
#' Modification is exactly "0,Arg[AnyN-term];" get R prepended and the mod
#' cleared; 4. "N,XXX->YYY[.]" substitutions replace the sequence letter and
#' drop the entry. Files with extra pFind_* columns (pUniFind_rescore) are
#' skipped (production behavior).
#'
#' @param dat
#'   数据框，原始搜库结果表。
#'   A data frame, the original search-result table.
#' @return
#'   返回校正后的数据框。
#'   The corrected data frame.
#' @keywords internal
revise_search_result <- function(dat, mode = "full") {
  # 后处理模式"none"＝无后处理
  # （限定式搜库、denovo 由回帖步骤处理）；"full"＝13C 删行＋末端丢失归一化
  # （开放式搜库）；"terminal"＝同 full 的两步（pUniFind rescore）。
  # SAP 替换已不再在任何模式施进序列（设计约定：X->Y 全程保持
  # SAP 修饰形式，仅在 pDeep2 理论谱图生成时代入——见下方第 4 步注释）。
  # Post-processing modes 
  # "none" = no post-processing (restricted search; denovo is handled by the
  # database look-back); "full" = 13C drop + terminal-loss normalization
  # (open search); "terminal" = the same two steps (pUniFind rescore).
  # SAP substitutions are no longer applied into sequences in ANY mode
  # (user ruling: X->Y stays in the SAP-modification form
  # throughout; substitution only when building pDeep2 spectra - see the
  # step-4 note below).
  if (identical(mode, "none")) return(dat)
  # 精确列名（pUniFind_rescore 表另含 pFind_Sequence / pFind_Modification 列；校正与合并都只针对普通 Sequence / Modification 列）
  # exact column names (pUniFind_rescore tables also carry pFind_Sequence /
  # pFind_Modification columns; the corrections and the merge both target
  # the plain Sequence / Modification columns only)
  col_mod <- which(colnames(dat) == "Modification")
  col_seq <- which(colnames(dat) == "Sequence")
  if (length(col_mod) != 1 || length(col_seq) != 1) {
    stop("expected exactly one Sequence and one Modification column")
  }
  # 1) 删除 _13C 行
  # 1) drop _13C rows
  c13 <- which(str_detect(dat[, col_mod], "_13C"))
  if (length(c13) > 0) dat <- dat[-c13, ]
  # 2) N 端 Met 丢失 → N 端乙酰化
  # 2) N-terminal Met loss -> N-terminal acetylation
  met <- which(str_detect(dat[, col_mod], "0,Met-"))
  if (length(met) > 0) {
    error_index <- c()
    # v1.0 优化原实现逐行 dat[i, col] <- 每行整表复制
    # （2.1s 的主因）——改为本地字符向量改写，循环后一次性写回两列
    seq_v <- as.character(dat[, col_seq]); mod_v <- as.character(dat[, col_mod])
    for (i in met) {
      seq0 <- seq_v[i]
      mod0 <- mod_v[i]
      seq_v[i] <- str_sub(seq0, 2, nchar(seq0))
      sp <- str_split(mod0, ";")[[1]]
      sp <- sp[-1]
      sp <- sp[sp != ""]
      if (length(sp) == 0) {
        mod_v[i] <- "0,Acetyl[AnyN-term];"
        next
      }
      sp_df <- t(as.data.frame(sapply(sp, function(x) str_split(x, ",")[[1]])))
      sp_df[, 1] <- as.numeric(sp_df[, 1]) - 1
      if ("0" %in% sp_df[, 1]) {
        error_index <- c(error_index, i)
        next
      }
      mod_v[i] <- paste(c("0,Acetyl[AnyN-term]",
                          sapply(seq_len(nrow(sp_df)), function(x)
                            paste(sp_df[x, ], collapse = ",")),
                          ""), collapse = ";")
    }
    if (length(error_index) > 0) {
      # 原口径：报错行整表删除（不只删这两列）——本地向量与全表同步删行
      seq_v <- seq_v[-error_index]; mod_v <- mod_v[-error_index]
      dat <- dat[-error_index,, drop = FALSE]
    }
    dat[[col_seq]] <- seq_v
    dat[[col_mod]] <- mod_v
  }
  # 3) N 端 Arg → 前补 R
  # 3) N-terminal Arg -> prepend R
  r_add <- which(str_detect(dat[, col_mod], "0,Arg"))
  if (length(r_add) > 0) {
    # v1.0 优化：同上，本地向量改写
    seq_v <- as.character(dat[, col_seq]); mod_v <- as.character(dat[, col_mod])
    for (i in r_add) {
      if (mod_v[i] == "0,Arg[AnyN-term];") {
        seq_v[i] <- paste0("R", seq_v[i])
        mod_v[i] <- ""
      } else {
        cat(i, "0,Arg[AnyN-term]; error", sep = "\t", "\n")
      }
    }
    dat[[col_seq]] <- seq_v
    dat[[col_mod]] <- mod_v
  }
  # 4) XXX->YYY 替换 → 序列字母替换：已删除——
  #    X->Y 全程保持 SAP 修饰形式（参考序列＋替换记号），参考蛋白/位置不丢失，
  #    分类按真 SAP 判据走；仅在生成 pDeep2 理论谱图时才代入序列变常规形式
  #    （pDeep2 输入构建处的 X->Y 代入与 stageB 的 apply_sap_to_sequence 已覆盖；
  #    打分的 by_ions_position 原生支持 "->" 记号，按目标残基质量计）。
  #    此前 full 模式在此处把替换施进序列，产生参考库中不存在的变体肽段，
  #    导致其蛋白归属为空（corrected_protein 空行的根因）。
  # 4) XXX->YYY substitution -> sequence letter replacement: REMOVED - X->Y tokens now stay in the SAP-modification form
  #    (reference sequence + substitution token) throughout the pipeline, so
  #    protein/position attribution is preserved and classification follows
  #    the real-SAP rule; the substitution into the plain sequence happens
  #    ONLY when building pDeep2 theoretical spectra (covered by the X->Y
  #    substitution at the pDeep2 input construction and by
  #    apply_sap_to_sequence in stage B; scoring via by_ions_position reads
  #    "->" tokens natively, using the target-residue mass).
  #    The former full-mode substitution here created variant peptides absent
  #    from the reference proteome - the root cause of empty corrected_protein.
  dat
}
#' denovo 数据库回帖denovo 序列必须回帖到搜库库
#' （I/L 等价＋N=GG 双向替换，组合上限 kmax），不能匹配的记"没有匹配"——
#' 该谱的 denovo 槽序列与修饰清空（不参与合并竞争），并加 "<slot>.DB_match"
#' 状态列（match:direct / match:GG->N / match:N->GG / no_match）。
#' Database look-back for de-novo sequences every denovo
#' sequence must map back into the search database under I/L equivalence with
#' N=GG substitution in both directions (combo cap kmax); a non-matching entry
#' is recorded as "no match" - its slot sequence and modification are blanked
#' (out of the merge competition) and a "<slot>.DB_match" status column is
#' added (match:direct / match:GG->N / match:N->GG / no_match).
#'
#' @param datalist data.frame, the merged table with "<slot>.Sequence" and
#'   "<slot>.Modification" columns.
#'   data.frame, 含 "<slot>.Sequence"/"<slot>.Modification" 列的合并表。
#' @param fasta_file character, the search FASTA (target-decoy tolerated;
#'   REV_/DECOY_/reversed_ prefixed accessions are dropped).
#'   字符串, 搜库 FASTA 路径（兼容 target-decoy，反库前缀自动剔除）。
#' @param slot integer, the denovo strategy slot number (default 5).
#'   整数, denovo 策略槽号（默认 5, 即第五策略）。
#' @param kmax integer, max number of isobaric-swap sites to expand (default 5).
#'   整数, N=GG 替换位点组合上限（默认 5）。
#' @return The table with the status column added and unmatched slots blanked.
#'   加状态列后的表；未匹配槽已清空。
#' @keywords internal
denovo_db_lookback <- function(datalist, fasta_file, slot = "e",
                               mode = c("keep", "all")) {
  mode <- match.arg(mode[1], c("keep", "all"))
  # denovo 数据库回帖肽段做 I→L、N→GG 替换，
  # 蛋白组数据库同样生成替换后数据库；保留在替换后库内有匹配的替换肽段序列，
  # 并标注对应蛋白与前后氨基酸（与 pFind-Filtered.spectra 的 Proteins/Positions
  # 列同格式：Proteins＝"ACC/"、Positions＝"pos,前,后/"）；无匹配记 no_match、
  # 该槽序列与修饰清空。位置映射：N→GG 的两个 G 均映射到同一原始残基位，
  # 故标注的前后氨基酸与位置均为原始蛋白坐标（pFind 0 起口径：前导数）。
  # Database look-back for denovo (user spec item 4): the peptide gets I->L
  # and N->GG substitutions and the proteome gets the same transformed
  # database; substituted peptides that match inside it are kept, annotated
  # with protein and flanking residues in the pFind-Filtered.spectra format
  # (Proteins="ACC/", Positions="pos,pre,post/"); non-matching entries are
  # marked no_match and their slot is blanked. Position mapping: both G's of
  # an N->GG map to the same original residue, so the annotation uses
  # original-protein coordinates (pFind 0-based: number of preceding residues).
  cn_s <- paste0(slot, ".Sequence"); cn_m <- paste0(slot, ".Modification")
  cn_p <- paste0(slot, ".Proteins");   cn_po <- paste0(slot, ".Positions")
  if (!cn_s %in% colnames(datalist))
    stop("denovo slot column missing: ", cn_s)
  # denovo=0 空槽跳过该槽序列/修饰全空说明 denovo 策略
  # 没有任何结果，直接返回原表，不做无意义的库构建与匹配
  # skip an all-empty slot when the denovo strategy
  # produced no results the slot is entirely blank - return the table as-is
  # instead of building the transformed database for nothing
  if (all(is.na(datalist[[cn_s]]) | !nzchar(as.character(datalist[[cn_s]]))) &&
      (!cn_m %in% colnames(datalist) ||
       all(is.na(datalist[[cn_m]]) | !nzchar(as.character(datalist[[cn_m]])))))
    return(datalist)
  fa <- read_fasta(fasta_file)
  fa <- fa[!grepl("^(REV_|DECOY_|reversed_)", names(fa))]
  fa <- toupper(fa)
  fa_orig <- fa                                            # 原始拼写（规格④第 3 点替换用）
  fa <- gsub("I", "L", fa, fixed = TRUE)                   # I→L（蛋白端）
  n_prot <- length(fa)
  # 替换后数据库与原始位置映射：逐蛋白 N→GG，两个 G 记同一原始位
  # transformed database + original-position mapping: N->GG per protein,
  # both G's of a pair map to the same original index
  tr_seq <- character(n_prot)
  tr_map <- vector("list", n_prot)
  for (k in seq_len(n_prot)) {
    s0 <- strsplit(fa[k], "", fixed = TRUE)[[1]]
    segs <- lapply(s0, function(ch) if (ch == "N") "GG" else ch)
    tr_seq[k] <- paste(unlist(segs), collapse = "")
    tr_map[[k]] <- rep(seq_along(s0), vapply(segs, nchar, integer(1)))  # N→GG 两位映射
  }
  big <- paste(tr_seq, collapse = "\n")
  cum_start <- cumsum(c(1L, nchar(tr_seq) + 1L))
  if (!cn_p %in% colnames(datalist)) {                      # 预建标注列（pFind 格式）
    datalist[[cn_p]] <- ""
    datalist[[cn_po]] <- ""
  }
  status <- rep("", nrow(datalist))
  # P8按变换后序列缓存匹配结果。实测 48,203 行 denovo 仅有
  # 22,581 个唯一序列（最高重复 91 次），命中蛋白/位置只查一次整库；no_match
  # 同样缓存。I→L 折叠与 N→GG 均为同向替换，同变换串的原始蛋白坐标一致。
  # P8 cache the match by the transformed sequence. 48,203
  # denovo rows hold only 22,581 unique sequences (top repeat 91x), so each
  # proteome-wide scan runs once per unique sequence; no_match is cached too.
  # Both I->L and N->GG are one-way, so equal transformed strings share the
  # same original-protein coordinates.
  sq_cache <- new.env(parent = emptyenv(), hash = TRUE)
  # v1.0 优化原实现逐行 datalist[i, col] <- 每行整表复制
  # （48k 行 × 4 列 O(n²) 拷贝）——改为本地向量改写，循环后一次写回。
  # cn_m 缺失且 mode≠all 且有 no_match 行时，原实现会顺带创建全空修饰列
  # （no_match 行写 ""）——按同样条件创建，保持列存在性口径一致。
  seq_v <- as.character(datalist[[cn_s]])
  mod_v <- if (cn_m %in% colnames(datalist)) as.character(datalist[[cn_m]])
           else rep("", nrow(datalist))
  pr_v <- as.character(datalist[[cn_p]])
  po_v <- as.character(datalist[[cn_po]])
  cn_m_was_present <- cn_m %in% colnames(datalist)
  for (i in seq_len(nrow(datalist))) {
    sq <- toupper(trimws(seq_v[i]))
    if (is.na(sq) || !nzchar(sq)) next
    sq_t <- gsub("N", "GG", gsub("I", "L", sq, fixed = TRUE), fixed = TRUE)
    hit <- if (exists(sq_t, envir = sq_cache, inherits = FALSE)) {
      get(sq_t, envir = sq_cache, inherits = FALSE)
    } else {
      m <- regexpr(sq_t, big, fixed = TRUE)                # 替换后库子串匹配
      v <- if (m[1] < 1L) NULL else {
        pi <- findInterval(m[1], cum_start)                # 命中蛋白
        toff <- m[1] - cum_start[pi] + 1L                  # 替换库内偏移
        tlen <- nchar(sq_t)
        oi <- tr_map[[pi]][toff:(toff + tlen - 1L)]        # 原始位映射
        list(pi = pi, o_start = oi[1], o_end = oi[length(oi)])
      }
      assign(sq_t, v, envir = sq_cache)
      v
    }
    if (is.null(hit)) {
      status[i] <- "no_match"
      # mode="all"no_match 的 denovo 序列不清空，
      # 带 DB_match=no_match 标记保留进竞争（SAP 类库外变异肽段由此参与仲裁）
      # mode="all" a no-match denovo sequence stays
      # in the competition flagged DB_match=no_match (SAP-type out-of-database
      # variants thereby enter the arbitration)
      if (!identical(mode, "all")) {
        seq_v[i] <- ""; mod_v[i] <- ""
        pr_v[i] <- ""; po_v[i] <- ""
      }
      next
    }
    pi <- hit$pi; o_start <- hit$o_start; o_end <- hit$o_end
    o_seq <- fa[pi]
    pre <- if (o_start > 1L) substr(fa_orig[pi], o_start - 1L, o_start - 1L) else "-"
    post <- if (o_end < nchar(fa_orig[pi])) substr(fa_orig[pi], o_end + 1L, o_end + 1L) else "-"
    status[i] <- "match:transformed"
    # 规格④第 3 点：denovo 序列替换为原始蛋白质数据库中的对应序列
    # spec item 4.3: replace the denovo sequence with the ORIGINAL database
    # fragment (untransformed form, as the protein has it)
    seq_v[i] <- substr(fa_orig[pi], o_start, o_end)
    pr_v[i] <- paste0(names(fa)[pi], "/")                  # Proteins 格式
    po_v[i] <- sprintf("%d,%s,%s/", o_start - 1L, pre, post)
  }
  datalist[[cn_s]] <- seq_v
  if (cn_m_was_present || (!identical(mode, "all") && any(status == "no_match"))) {
    datalist[[cn_m]] <- mod_v
  }
  datalist[[cn_p]] <- pr_v
  datalist[[cn_po]] <- po_v
  datalist[[paste0(slot, ".DB_match")]] <- status
  cat(sprintf("denovo database look-back / denovo 数据库回帖: %d entries, %d matched, %d no-match blanked\n",
              sum(!is.na(status) & nzchar(status)), sum(status == "match:transformed"),
              sum(status == "no_match")))
  datalist
}

#' 用户步骤 2：多策略搜库结果合并（多引擎版）。
#' User step 2: merge the search strategies (multi-engine edition).
#'
#' 策略目录表每种策略类型硬编码为单字母槽号
#' （a=原始结果，b..z=补充策略），合并按字母序自然排序；字母开头即合法 R 列名，
#' 无需 X 前缀。p..z 为通用槽（type="other"，可重复挂任意引擎）。
#' The strategy catalog every strategy type maps to a
#' fixed single-letter slot (a = original result, b..z = supplementary);
#' the merge orders columns by letter naturally, and letter-leading column
#' names are valid R names without an X prefix. p..z are generic slots
#' (type = "other", reusable for any engine).
#'
#' @keywords internal
STRATEGY_SLOT <- c(original = "a", pfind_open = "b", pfind_restrict = "c",
                   punifind_rescore = "d", punifind_denovo = "e",
                   #  R5pFind 系列槽位不动；
                   # 其余引擎按软件名字典序重排 f–p（Comet→X!Tandem）
                   comet = "f", mascot = "g", maxquant = "h",
                   metamorpheus = "i", msfragger = "j", msfragger_boost = "k",
                   msfragger_open = "l", peaks = "m", pnovo = "n",
                   sage = "o", xtandem = "p", pfind_standard = "q")
STRATEGY_REVISE <- c(pfind_open = "full", pfind_restrict = "none",
                     pfind_standard = "none",
                     punifind_rescore = "terminal", punifind_denovo = "none",
                     maxquant = "none", msfragger_open = "none",
                     msfragger = "none", sage = "none", comet = "none",
                     xtandem = "none", metamorpheus = "none", peaks = "none",
                     mascot = "none", pnovo = "none", msfragger_boost = "none",
                   other = "none")
STRATEGY_LOOKBACK <- c(punifind_denovo = TRUE, pnovo = TRUE)
#' 由策略类型向量求槽位字母：固定类型按目录表取字母（重复报错）；
#' other 类型依次占用 q..z 通用槽。original 只属于策略 1（psm_file），
#' 不允许出现在补充策略类型中。
#' Resolve slot letters from a strategy-type vector: fixed types take their
#' catalog letters (duplicates are an error); "other" types fill the q..z
#' generic slots in order. "original" belongs to strategy 1 (psm_file) only
#' and is not allowed among the supplementary types.
#' @keywords internal
strategy_slots <- function(types) {
  types <- as.character(types)
  bad <- setdiff(types, c(names(STRATEGY_SLOT), "other"))
  if (length(bad) > 0) stop("未知策略类型：", paste(bad, collapse = ", "))
  if (any(types == "original")) stop("original 只属于策略 1（psm_file）")
  out <- character(length(types))
  used <- character(0)
  free_other <- letters[17:26]
  for (i in seq_along(types)) {
    ty <- types[i]
    if (ty == "other") {
      fr <- free_other[!free_other %in% used]
      if (length(fr) == 0) stop("other 通用槽已用完（最多 26 个策略）")
      out[i] <- fr[1]
    } else {
      out[i] <- STRATEGY_SLOT[[ty]]
      if (out[i] %in% used) stop("策略类型重复：", ty)
    }
    used <- c(used, out[i])
  }
  out
}

#' 用户步骤 2：多策略搜库结果合并。
#' User step 2: merge the search strategies.
#'
#' 把常规搜库结果（策略 a）与各补充搜库结果按谱键控合并；调用方应传 REVISED 表；
#' 序列存原始（I→L / N→GG 归一属于 SAP 与蛋白匹配步骤）；输出 1.combine_result.spectra
#' 与 peptide.txt（pDeep2 输入，仅在本次合并后生成——生产约定）。
#' Merges the regular search result (strategy a) with the supplementary
#' results into one spectrum-keyed table; callers pass the REVISED tables;
#' sequences are stored raw; outputs 1.combine_result.spectra and peptide.txt
#' (pDeep2 input, generated only after this merge - production convention).
#'
#' @param psm_file
#'   字符串，策略 a 常规搜库结果。
#'   Character string, the strategy-a regular search result.
#' @param search_files
#'   字符向量，补充搜库结果（每项一个策略，槽位由 slot_letters 指定）。
#'   Character vector, the supplementary search results (one strategy each; the
#'   slot letters are given by slot_letters).
#' @param out_dir
#'   字符串，输出目录（merge\）。
#'   Character string, the output directory (merge\).
#' @param slot_letters
#'   字符向量，与 search_files 等长的槽位字母（缺省 b、c、…）。
#'   Character vector, the slot letters matching search_files (default b, c, ...).
#' @param denovo_slots
#'   字符向量，需要数据库回帖的 denovo 槽位字母（如 "e"）。
#'   Character vector, the denovo slot letters needing database look-back (e.g. "e").
#' @param denovo_lookback_mode
#'   字符串，denovo 回帖模式"keep"＝仅保留能回帖到
#'   数据库的肽段（默认，历史行为）；"all"＝放开全部结果——no_match 的
#'   denovo 序列不清空、带 DB_match=no_match 标记参与竞争（SAP 类库外变异
#'   肽段由此进入仲裁）。
#'   Character, the denovo look-back mode "keep" =
#'   only peptides mapping back into the database survive (default,
#'   historical behaviour); "all" = release everything - no-match denovo
#'   sequences stay in the competition flagged DB_match=no_match (so SAP-type
#'   out-of-database variants enter the arbitration).
#' @return
#'   返回合并表的数据框。
#'   The merged table as a data frame.
#' @keywords internal
sanitize_pdeep2_mods <- function(mod, len) {
  if (is.na(mod) || !nzchar(mod)) return(mod)
  toks <- strsplit(mod, ";", fixed = TRUE)[[1]]
  keep <- vapply(toks, function(tk) {
    if (!nzchar(tk)) return(TRUE)
    pos <- suppressWarnings(as.integer(sub(",.*$", "", tk)))
    is.na(pos) || (pos >= 0 && pos <= len + 1L)
  }, logical(1))
  paste(toks[keep], collapse = ";")
}
#' 槽位状态列（Peptide_state / Modification_state）的向量化计算。
#' Vectorized computation of the slot-state columns (Peptide_state /
#' Modification_state).
#'
#' 输入已归一化的 n×k 字符矩阵（每槽一列），输出每行分组标签：全部非空
#' 值相同或全空 → "Same"；否则按首现顺序把同值槽的字母并串、组间用 ";"
#' 连接（如 "ab;c"）。逐行 sapply 版本在大表上耗时小时级；本版把每行编码为整数键，标签闭包只对唯一键执行
#' （唯一键数远小于行数），再按键索引还原。
#' Input: a normalized n-by-k character matrix (one column per slot). Output:
#' per-row group labels - "Same" when all non-empty values agree or the row is
#' empty, else the slot letters of equal values concatenated in
#' first-occurrence order, groups joined by ";" (e.g. "ab;c"). The per-row
#' sapply version cost hours on large tables; this one encodes each row as an integer key, runs the label closure
#' once per unique key (far fewer than rows) and maps back by key.
#'
#' @param m
#'   字符矩阵，n×k，已按调用方规则归一化。
#'   Character matrix n x k, already normalized by the caller's rules.
#' @param sym
#'   字符向量，槽字母（按列顺序）。
#'   Character vector, slot letters in column order.
#' @return
#'   字符向量，长度 n 的状态标签。
#'   Character vector of length n with the state labels.
#' @export
slot_state_matrix <- function(m, sym) {
  n <- nrow(m); k <- ncol(m)
  if (n == 0L) return(character(0))
  f <- factor(m)                              # 全矩阵唯一值 → 整数编码
  codes <- as.integer(f)
  dim(codes) <- dim(m)
  codes[is.na(codes)] <- 0L                   # NA 视同空值
  if ("" %in% levels(f)) codes[!is.na(m) & m == ""] <- 0L  # 空串视同空值
  if (k == 1L) {
    keys <- as.character(codes)
  } else {
    keymat <- matrix(as.character(codes), n, k)
    keys <- do.call(paste, c(lapply(seq_len(k), function(j) keymat[, j]),
                             sep = ","))
  }
  uk <- unique(keys)
  lab <- vapply(uk, function(ky) {
    cc <- as.integer(strsplit(ky, ",", fixed = TRUE)[[1]])
    u <- unique(cc[cc != 0L])
    if (length(u) <= 1L) return("Same")
    paste(vapply(u, function(z) paste(sym[cc == z], collapse = ""), ""),
          collapse = ";")
  }, "")
  unname(setNames(lab, uk)[keys])
}
# 全候选模式合并表装配忽略各
# 引擎 1% FDR 过滤，每个策略槽每张谱取引擎自身排序（文件序）前 5 个肽段+
# 修饰组合进竞争池。输出表与常规合并完全同构（逐谱 1..K 行、K 为各槽候选数
# 最大者且封顶 5，槽位字母列），下游 stageA/peptide.txt/报告通用消费。
# top-5 candidate table assembly (candidate_scope=top5): each strategy slot
# contributes its first 5 (file-order) peptide+modification candidates per
# spectrum, ignoring the engines' 1% FDR filtering. The output is shaped
# exactly like the regular merge (1..K rows per spectrum, K = max candidate
# count across slots capped at 5) so stageA / peptide.txt and the reports
# consume it unchanged.
build_top5_datalist <- function(datalist, search_files, slot_letters) {
  cap5 <- function(x) x[seq_len(min(5L, length(x)))]
  fn0 <- as.character(datalist$File_Name)
  spec_key <- unique(fn0)
  n_spec <- length(spec_key)
  first0 <- match(spec_key, fn0)
  keep_annot <- setdiff(names(datalist),
                        c("File_Name",
                          grep("\\.Sequence$", names(datalist), value = TRUE),
                          grep("\\.Modification$", names(datalist), value = TRUE)))
  slot_data <- vector("list", length(search_files))
  K <- rep(1L, n_spec)
  for (index in seq_along(search_files)) {
    sf <- search_files[index]
    sl <- slot_letters[index]
    cat("merging (top5) [", sl, "]", basename(sf), "\n")
    hdr <- strsplit(readLines(sf, n = 1), "\t")[[1]]
    cc <- rep("NULL", length(hdr))
    for (col in c("File_Name", "Sequence", "Modification")) {
      i <- which(hdr == col)
      if (length(i) == 1) cc[i] <- "character"
    }
    temp <- read.delim(sf, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       colClasses = cc, na.strings = character(0))
    if (nrow(temp) == 0) next
    if (!"Modification" %in% names(temp)) temp$Modification <- ""
    if (!"Sequence" %in% names(temp)) temp$Sequence <- ""
    fn <- as.character(temp$File_Name)
    g <- split(seq_len(nrow(temp)), match(fn, spec_key))
    g <- g[as.character(seq_len(n_spec))]
    g[lengths(g) == 0L] <- list(integer(0))
    g <- lapply(g, cap5)
    K <- pmax(K, lengths(g))
    slot_data[[index]] <- list(g = g, seq = as.character(temp$Sequence),
                               mod = as.character(temp$Modification))
  }
  # 槽 a：psm 自身多行（每谱候选按文件序），同样封顶 5
  ga <- split(seq_len(nrow(datalist)), match(fn0, spec_key))
  ga <- ga[as.character(seq_len(n_spec))]
  ga[lengths(ga) == 0L] <- list(integer(0))
  ga <- lapply(ga, cap5)
  K <- pmin(pmax(K, lengths(ga), 1L), 5L)
  a_seq <- as.character(datalist[["a.Sequence"]])
  a_mod <- as.character(datalist[["a.Modification"]])
  spec_rep <- rep(seq_len(n_spec), K)
  k_rep <- unlist(lapply(K, seq_len), use.names = FALSE)
  base <- datalist[first0, keep_annot, drop = FALSE]
  out <- data.frame(File_Name = spec_key[spec_rep], stringsAsFactors = FALSE)
  out[keep_annot] <- base[spec_rep,, drop = FALSE]
  pick <- function(g, v) {
    r <- rep("", length(spec_rep))
    for (j in seq_len(n_spec)) {
      idx <- g[[j]]
      if (length(idx)) {
        kk <- k_rep[spec_rep == j]
        ok <- kk <= length(idx)
        r[spec_rep == j][ok] <- v[idx[kk[ok]]]
      }
    }
    r
  }
  out[["a.Sequence"]] <- pick(ga, a_seq)
  out[["a.Modification"]] <- pick(ga, a_mod)
  for (index in seq_along(search_files)) {
    sln <- paste0(slot_letters[index], ".Sequence")
    ml <- paste0(slot_letters[index], ".Modification")
    if (is.null(slot_data[[index]])) {
      out[[sln]] <- rep("", nrow(out))
      out[[ml]] <- rep("", nrow(out))
      next
    }
    out[[sln]] <- pick(slot_data[[index]]$g, slot_data[[index]]$seq)
    out[[ml]] <- pick(slot_data[[index]]$g, slot_data[[index]]$mod)
  }
  out
}

merge_search_results <- function(psm_file, search_files, out_dir,
                             fasta_file = NULL, slot_letters = NULL,
                             denovo_slots = NULL,
                             denovo_lookback_mode = c("keep", "all"),
                             candidate_scope = "filtered",
                             strip_labels = NULL,
                             min_pep_len = 7L) {
  # 多引擎版：槽位用单字母（a=原始结果，b..z=补充策略），字母开头即合法 R 列名
  # multi-engine edition: slots are single letters (a = original, b..z =
  # supplementary); a letter-leading column name is a valid R name
  if (length(search_files) < 2) stop("need at least two supplementary search results")
  if (is.null(slot_letters)) slot_letters <- letters[seq_along(search_files) + 1L]
  if (length(slot_letters) != length(search_files) ||
      any(!slot_letters %in% letters) || anyDuplicated(slot_letters)) {
    stop("slot_letters must be unique letters matching search_files")
  }
  if (any(slot_letters == "a")) stop("slot 'a' is reserved for the original result")
  if (is.null(denovo_slots)) denovo_slots <- character(0)
  if (any(!denovo_slots %in% slot_letters)) {
    stop("denovo_slots must be a subset of slot_letters")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # 列类型守卫psm 的 Modification 列全空时（如
  # pFind 3.2.3 双通道搜库的某些输出）R 会把该列推断为 logical 并整体读成
  # NA，随后 write.table 默认 na="NA" 写成字面 "NA"，令下游 pDeep2 崩溃。
  # 显式把 File_Name/Sequence/Modification 读为 character。
  # column-type guard when the psm Modification column
  # is entirely empty (e.g. certain pFind 3.2.3 two-channel outputs) R infers
  # it as logical and reads it all-NA; write.table's default na="NA" then
  # emits literal "NA" strings that crash pDeep2 downstream. Force the three
  # columns to character explicitly.
  hdr0 <- strsplit(readLines(psm_file, n = 1), "\t")[[1]]
  cc0 <- rep(NA_character_, length(hdr0))
  for (col0 in c("File_Name", "Sequence", "Modification")) {
    i0 <- which(hdr0 == col0)
    if (length(i0) == 1) cc0[i0] <- "character"
  }
  datalist <- read.delim(psm_file, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = cc0)
  # 修饰列规范化：NA / 字面 "NA" 一律按空修饰处理
  # normalize modification columns: NA / literal "NA" become empty
  for (mc in c("a.Modification", "Modification")) {
    if (mc %in% colnames(datalist)) {
      datalist[[mc]][is.na(datalist[[mc]])] <- ""
      datalist[[mc]][datalist[[mc]] == "NA"] <- ""
    }
  }
  spectra_list <- datalist$File_Name
  # 原位改名（不复制整表）：Sequence/Modification -> a.Sequence/a.Modification；
  # 槽列追加与列序重排都只做一次（R2 内存修复：早期实现每引擎一次
  # cbind 整表拷贝，十来个引擎累计 70GB+ 级拷贝流量，RSS 被 GC 节律推高）
  # rename in place (no whole-table copy): Sequence/Modification ->
  # a.Sequence/a.Modification; slot appends and the final column reorder each
  # happen once (R2 memory fix: the old per-engine cbind copied the
  # whole growing table, ~70GB+ cumulative copy traffic across a dozen engines)
  names(datalist)[names(datalist) %in% c("Sequence", "Modification")] <-
    c("a.Sequence", "a.Modification")

  # 全候选模式跳过常规逐槽 top-1 映射，改由
  # build_top5_datalist 装配逐谱 1..K 行、每槽第 k 个候选的合并表
  if (is.null(candidate_scope) || !length(candidate_scope)) candidate_scope <- "filtered"
  scope_top5 <- identical(candidate_scope[1], "top5")
  if (scope_top5) {
    datalist <- build_top5_datalist(datalist, search_files, slot_letters)
  }
  if (!scope_top5)
    for (index in seq_along(search_files)) {
    sf <- search_files[index]
    sl <- slot_letters[index]
    cat("merging [", sl, "]", basename(sf), "\n")
    hdr <- strsplit(readLines(sf, n = 1), "\t")[[1]]
    cc <- rep("NULL", length(hdr))
    for (col in c("File_Name", "Sequence", "Modification")) {
      i <- which(hdr == col)
      if (length(i) == 1) cc[i] <- "character"
    }
    temp <- read.delim(sf, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       colClasses = cc, na.strings = character(0))
    # 列齐性守卫部分适配器产物（如 MetaMorpheus
    # 修订文件）只有 File_Name/Sequence 而无 Modification 列，旧逻辑直接按三列
    # 改名致 names[3] vs 2 崩溃——缺失列补空
    # column-completeness guard some adapter
    # outputs (e.g. MetaMorpheus revised) lack the Modification column and the
    # old 3-name assignment crashed with names[3] vs 2 - fill the missing empty
    # 空结果容错引擎结果文件仅表头（如
    # MetaMorpheus 自动搜库 0 产出）时该槽整列置空，不中断合并——须先于列守卫
    # （0 行 data.frame 上 `$<-` 赋长度 1 会崩）
    # empty-result tolerance: a header-only engine file blanks the whole slot;
    # it MUST come before the column guards ($<- with length 1 on a 0-row frame
    # crashes)
    if (nrow(temp) == 0) {
      temp <- data.frame(File_Name = spectra_list,
                         Sequence = rep("", length(spectra_list)),
                         Modification = rep("", length(spectra_list)),
                         stringsAsFactors = FALSE)
      cat("  (empty strategy file: slot filled blank)\n")
      # 空结果分支：行序与 spectra_list 一致，排序恒等（含重复 File_Name 亦同）
      ord <- seq_len(nrow(temp))
    } else {
      # 列齐性守卫：部分适配器产物（如 MetaMorpheus 修订文件）只有
      # File_Name/Sequence 而无 Modification 列，缺失列补空
      # column-completeness guard: some adapter outputs lack the Modification
      # column - fill the missing one empty
      if (!"Modification" %in% names(temp)) temp$Modification <- ""
      if (!"Sequence" %in% names(temp)) temp$Sequence <- ""
      # v1.0 优化原实现对 57k 键向量哈希三次
      # （%in% + setdiff + order(match)），改为一次 match 复用于三处：
      # mt2 = temp 行在 spectra_list 中的位置；hit_pos = 被命中的位置；
      # 未命中 = unique(spectra_list) 中未出现的位置（与 setdiff 口径一致，
      # 含重复 File_Name 场景）；最终排序直接 order 位置向量（不再二次哈希）
      mt2 <- match(temp$File_Name, spectra_list)
      v <- !is.na(mt2)
      temp <- temp[v,, drop = FALSE]
      mt2 <- mt2[v]
      # P2-3unmatched 必须用 setdiff
      # 值域口径——此前拿"唯一序下标"去比"位置"跨索引空间比较，重复
      # File_Name 时由旧实现的响亮报错变静默错位；setdiff 恢复旧语义
      # （%in% 本就是哈希实现，无性能损失）
      unmatched <- setdiff(spectra_list, temp$File_Name)
      if (length(unmatched) > 0) {
        u <- setNames(data.frame(File_Name = unmatched, X = "", Y = "",
                                 stringsAsFactors = FALSE), colnames(temp))
        temp <- rbind(temp, u)
      }
      # 排序位置向量 = 已命中行位置 + 未命中行的 match 首现位置（与旧
      # order(match(...)) 逐字节同义）
      ord_positions <- c(mt2, if (length(unmatched) > 0)
        match(unmatched, spectra_list) else integer(0))
      ord <- order(ord_positions)
    }
    colnames(temp) <- c("File_Name", paste0(sl, ".Sequence"),
                        paste0(sl, ".Modification"))
    # 槽列按名追加：只移动新列指针，不复制整表（早期实现 cbind 每引擎全帧拷贝）
    # append by name: only the two new column vectors move, no whole-table
    # copy (the old cbind duplicated the entire frame per engine)
    datalist[[colnames(temp)[2]]] <- temp[[2]][ord]
    datalist[[colnames(temp)[3]]] <- temp[[3]][ord]
  }
  # 列序一次性对齐：槽对（Sequence 与 Modification 相邻）移到末尾——仅此一次
  # 整表列指针拷贝，替代早期实现逐引擎 cbind 整表数据拷贝
  # one-shot column order: slot pairs (Sequence and Modification adjacent)
  # move to the end - a single pointer-level copy replaces the old per-engine
  # whole-table data copies
  seq_nm <- grep("\\.Sequence$", names(datalist), value = TRUE)
  mod_nm <- grep("\\.Modification$", names(datalist), value = TRUE)
  keep_nm <- setdiff(names(datalist), c(seq_nm, mod_nm))
  datalist <- datalist[c(keep_nm, as.vector(rbind(seq_nm, mod_nm)))]
  iso_mem_guard(40000, "合并：主表与槽列装配完成")

  # 候选长度阈值过滤——所有槽位序列长度
  # 低于 min_pep_len（默认 7，config_step1 键）的候选置空（退出竞争、不进
  # peptide.txt 预测输入）。pFind 搜库侧 len_lower 只约束搜库自身，denovo/
  # psm_file 候选（如 4 肽 TRHW）不受控——短肽预测软件做不了（Prosit 硬限
  # 7-30 aa）、理论离子信息量不足；须在状态分组前过滤，状态反映过滤后集合。
  #  candidate length gate - every slot
  # sequence shorter than min_pep_len (default 7, config_step1 key) is
  # blanked (out of the competition and of the peptide.txt prediction input).
  # pFind's len_lower only constrains the search itself; denovo/psm_file
  # candidates (e.g. the 4-mer TRHW) were uncontrolled - predictors cannot
  # handle them (Prosit hard-limits 7-30 aa) and their theoretical ions carry
  # too little information. Filtering precedes the state grouping so the
  # states reflect the filtered set.
  {
    mpl <- suppressWarnings(as.integer(min_pep_len[1]))
    if (!is.na(mpl) && mpl >= 4L) {
      n_short <- 0L; n_rows_hit <- 0L
      for (ci in grep("\\.Sequence$", colnames(datalist))) {
        sq <- as.character(datalist[[ci]])
        bad <- nzchar(sq) & nchar(sq) < mpl
        if (any(bad)) {
          mod_cn <- sub("\\.Sequence$", ".Modification", colnames(datalist)[ci])
          datalist[[ci]][bad] <- ""
          if (mod_cn %in% colnames(datalist))
            datalist[[mod_cn]][bad] <- ""
          n_short <- n_short + sum(bad)
          n_rows_hit <- n_rows_hit + 1L
        }
      }
      if (n_short > 0L) {
        cat(sprintf("[合并] 长度阈值过滤：min_pep_len=%d，置空短肽候选 %d 个（涉及槽位 %d 列；不进竞争与预测输入）\n",
                    mpl, n_short, n_rows_hit))
        cat(sprintf("[merge] length gate: min_pep_len=%d, blanked %d short candidates across %d slot column(s)\n",
                    mpl, n_short, n_rows_hit))
      }
    }
  }

  # denovo 数据库回帖无库内匹配的 denovo 槽清空并记状态；
  # 需在 Peptide_state/Modification_state 计算前执行（状态反映过滤后的集合）。
  # 多引擎版：denovo 槽为字母向量（如 "e"），可多个
  # denovo database look-back: blank non-matching denovo slots before the
  # Peptide_state/Modification_state grouping (states reflect the filtered set).
  # multi-engine edition: denovo slots are a letter vector (e.g. "e"), several
  # are allowed
  if (!is.null(fasta_file)) {
    for (ds in denovo_slots)
      datalist <- denovo_db_lookback(datalist, fasta_file, ds,
                                      mode = denovo_lookback_mode[1])
  }
  datalist$Peptide_state <- ""
  datalist$Modification_state <- ""
  col_names <- colnames(datalist)
  seq_idx <- which(grepl("\\.Sequence$", col_names))
  seq_sym <- sapply(col_names[seq_idx], function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
  mod_idx <- which(grepl("\\.Modification$", col_names))

  cat(sprintf("[合并] 向量化计算 Peptide_state / Modification_state（%d 行 × %d 槽）...\n",
              nrow(datalist), length(seq_idx)))
  t_st <- proc.time()
  # 逐行 sapply 曾在大表上跑出小时级耗时；此处把每行 k 个
  # 槽值编码为整数键，标签闭包只对唯一键执行（唯一键数远小于行数），
  # 再按键索引还原
  # the per-row sapply once cost hours on large tables (same family as the R2
  # stall); each row's k slot values are encoded as an integer key, the label
  # closure runs once per unique key (far fewer than rows) and labels are
  # mapped back by key
  # 生产在 I→L 归一化后分组状态（2.combine_result_spectra.R）
  # 规格④：denovo 槽存替换后序列（N→GG），分组同步做 N→GG 归一
  # production groups the states after I->L normalization
  # (2.combine_result_spectra.R); spec item 4: denovo slots hold substituted
  # sequences (N->GG); the state grouping therefore also normalizes N->GG so
  # equivalent forms unify
  s_mat <- as.matrix(datalist[, seq_idx, drop = FALSE])
  s_mat[] <- gsub("N", "GG", gsub("I", "L", s_mat, fixed = TRUE), fixed = TRUE)
  datalist$Peptide_state <- slot_state_matrix(s_mat, seq_sym)
  m_mat <- as.matrix(datalist[, mod_idx, drop = FALSE])
  m_mat[] <- gsub("ProteinN-term", "AnyN-term",
                  gsub("I", "L", m_mat, fixed = TRUE), fixed = TRUE)
  datalist$Modification_state <- slot_state_matrix(m_mat, seq_sym)
  cat(sprintf("[合并] 状态计算完成：%.0fs，R 堆 %.0f MB\n",
              (proc.time() - t_st)[3], iso_mem_mb()))

  write.table(datalist, file.path(out_dir, "1.combine_result.spectra"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)

  cat("[合并] 构建 peptide.txt（SAP 展开 + 越界修饰清理）...\n")
  t_pd <- proc.time()
  pd_parts <- list()
  for (pc in seq_idx) {
    tmp <- datalist[, c(pc, pc + 1, which(colnames(datalist) == "Charge"))]
    colnames(tmp) <- c("peptide", "modification", "charge")
    # （无标记通道=0）：pDeep2 对去标记肽段预测——肽表修饰剥除标记修饰
    if (!is.null(strip_labels) && length(strip_labels) > 0)
      tmp$modification <- strip_label_mods(tmp$modification, strip_labels)
    # SAP 型替换修饰（X->Y）代入序列后再给 pDeep2（如 DASAVILVSDAR+1,Asp->Asn[D]
    # -> NASAVILVSDAR），避免 pDeep2 无法预测返回空白
    # apply SAP-type substitutions (X->Y) into the sequence before handing to
    # pDeep2 so predictions are not blank 
    # R2 内存修复：逐行循环只对含 "->" 的行执行（apply_sap_to_sequence
    # 对无 "->" 的行原样返回，全表逐行调用曾是大表小时级开销的一部分）
    # R2 memory fix: the per-row loop only runs on rows carrying
    # "->" (apply_sap_to_sequence returns rows without "->" unchanged; the
    # full-table per-row calls were part of the large-table hour-scale cost)
    sap <- grepl("->", tmp$modification, fixed = TRUE)
    extras <- NULL
    if (any(sap)) {
      #：Xle 歧义行列表收集 + 循环后单次拼接（原逐行 rbind 累积拷贝，
      # SAP 行数大时二次方退化——合并段规模化审查）
      # collect Xle-ambiguity rows in a list and bind once after the loop (the
      # per-row rbind accumulated copies and went quadratic on SAP-heavy runs)
      xle_rows <- vector("list", sum(sap)); xi <- 0L
      for (r in which(sap)) {
        vars <- apply_sap_variants(as.character(tmp$peptide[r]),
                                   as.character(tmp$modification[r]))
        if (length(vars) == 1) {
          tmp$peptide[r] <- vars[[1]]$seq
          tmp$modification[r] <- vars[[1]]$mod
        } else {
          # Xle 歧义：展开为 I/L 两行（pDeep2 分别预测，stageB 取最高余弦）
          # Xle ambiguity: expand to I/L rows (pDeep2 predicts both; stageB takes max)
          extra <- tmp[r, ]
          tmp$peptide[r] <- vars[[1]]$seq; tmp$modification[r] <- vars[[1]]$mod
          extra$peptide <- vars[[2]]$seq; extra$modification <- vars[[2]]$mod
          xi <- xi + 1L
          xle_rows[[xi]] <- extra
        }
      }
      if (xi > 0L) extras <- do.call(rbind, xle_rows[seq_len(xi)])
    }
    # 防御过滤位点超出肽长的修饰——个别引擎会输出错位修饰
    # （如 33 肽上的 64,Oxidation[M]），pDeep2 按位点对列表赋值会越界崩溃；
    # 越界修饰是无意义产物，直接丢弃。：只对非空修饰行执行
    # defence drop mods whose position exceeds the peptide -
    # some engines emit mispositioned mods (e.g. 64,Oxidation[M] on a 33-mer)
    # which crash pDeep2's position-indexed assignment; such mods are
    # meaningless and are discarded.: non-empty mods only
    nonempty <- which(nzchar(tmp$modification))
    if (length(nonempty) > 0) {
      # v1.0 优化原 vapply 逐格调用 sanitize_pdeep2_mods
      # （23 万次 strsplit+vapply 级）改为展平一次判定。注意旧函数的隐藏
      # 归一化：R 的 strsplit 恰好剥掉一个尾空串（"x;;;"→"x","",""），故旧
      # 结果 = 剥一个尾分号 + 丢越界 token + 按 ";" 重组——无越界 token 的行
      # 用 sub(";$","") 与之逐字节等价；仅含越界 token 的行（极少）逐行重组。
      # 保留判据逐 token 等价（空 token、位置不可识别保留，仅越界整数位点丢弃）
      mods_n <- as.character(tmp$modification)[nonempty]
      lens_n <- nchar(tmp$peptide)[nonempty]
      toks <- strsplit(mods_n, ";", fixed = TRUE)
      flat_t <- unlist(toks, use.names = FALSE)
      g <- rep(seq_along(nonempty), lengths(toks))
      pos <- suppressWarnings(as.integer(sub(",.*$", "", flat_t)))
      bad <- !is.na(pos) & (pos < 0L | pos > lens_n[g] + 1L)
      if (any(bad)) {
        # 越界行按旧 paste 重组（结果可能以 ";" 结尾，不再剥）；其余行只做
        # 单尾分号剥离（strsplit 恰剥一个尾空串的等价物）
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
    if (!is.null(extras)) tmp <- rbind(tmp, extras)
    pd_parts[[length(pd_parts) + 1L]] <- unique(tmp)
  }
  # 列表收集 + 一次拼接（早期实现逐槽 rbind 累积拷贝）
  # collect in a list and bind once (the old per-slot rbind accumulated copies)
  peptide_df <- do.call(rbind, pd_parts)
  peptide_df <- unique(peptide_df)
  peptide_df <- peptide_df[peptide_df$peptide != "", ]
  write.table(peptide_df, file.path(out_dir, "peptide.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat(sprintf("[合并] peptide.txt 完成：%.0fs（%d 行），R 堆 %.0f MB\n",
              (proc.time() - t_pd)[3], nrow(peptide_df), iso_mem_mb()))
  cat("written: 1.combine_result.spectra (", nrow(datalist),
      "rows) and peptide.txt (", nrow(peptide_df), "peptides)\n")
  list(merged_rows = nrow(datalist), peptides = nrow(peptide_df))
}
#' 读 pDeep2 predict.txt（pDeep2 唯一输出）。
#' Read a pDeep2 predict.txt file (the ONLY pDeep2 output).
#'
#' 每块一个条目，按 TITLE（"peptide|modification|charge"）键控；强度取
#' "m/z intensity label" 行的第二列。返回 list(keys, intensities, idx)。
#' One entry per block, keyed by the TITLE ("peptide|modification|charge");
#' intensities = the second column of the "m/z intensity label" lines.
#' Returns list(keys, intensities, idx).
#'
#' @param path
#'   字符串，predict.txt 的路径。
#'   Character string, the path to predict.txt.
#' @return
#'   返回 list(keys=标题向量, intensities=强度向量列表, idx=键→下标命名向量)。
#'   A list(keys = titles, intensities = list of intensity vectors,
#'   idx = named key-to-index vector).
#' @keywords internal
read_predict_txt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begins <- which(lines == "BEGIN IONS")
  ends <- which(lines == "END IONS")
  # P1-1BEGIN/END 计数不等时下方
  # rep() 长度不匹配报 invalid 'times' argument，用户无法定位——对齐长度并
  # 给出可行动错误（截断的 predict.txt 删掉重跑 pDeep2）
  if (length(begins) != length(ends)) {
    stop("predict.txt 的 BEGIN IONS / END IONS 不配对（", length(begins),
         " vs ", length(ends), "）：文件可能被截断，请删除后重跑 pDeep2：", path)
  }
  n <- length(begins)
  keys <- rep(NA_character_, n)
  intensities <- vector("list", n)
  if (n > 0) {
    # v1.0 优化原实现逐块拷贝 lines[begins[i]:ends[i]] +
    # 块内 grepl/startsWith 扫描（19k 块、predict.txt 解析数十秒级）——改为
    # 一次性展平分组：blk_idx = 每行所属块号；blk_lines = 全块行（按序拼接）
    # block_lens 钳制非负（病态块 ends<begins 不崩溃，该块键留 NA）
    block_lens <- pmax(ends - begins + 1L, 0L)
    blk_idx <- rep(seq_len(n), block_lens)
    blk_lines <- lines[rep(begins, block_lens) + sequence(block_lens) - 1L]
    # TITLE：每块第一个 ^TITLE= 行（块内行序即组内序，duplicated 取首现）
    t <- which(startsWith(blk_lines, "TITLE="))
    if (length(t) > 0) {
      tgrp <- blk_idx[t]
      keep_t <- !duplicated(tgrp)
      keys[tgrp[keep_t]] <- sub("^TITLE=", "", blk_lines[t[keep_t]])
    }
    # 每类强度行：全展平向量上 startsWith 一次命中全部，再按块取首行拆分
    # per_prefix_values 返回长度 n 的列表（缺该类型的块为 NULL → 组装时 numeric(0)）
    per_prefix_values <- function(prefix) {
      out <- vector("list", n)
      hit_lines <- which(startsWith(blk_lines, paste0(prefix, "=")))
      if (length(hit_lines) > 0) {
        grp <- blk_idx[hit_lines]
        keep <- !duplicated(grp)
        hit_lines <- hit_lines[keep]; grp <- grp[keep]
        vals <- lapply(strsplit(sub(paste0(prefix, "="), "",
                                    blk_lines[hit_lines], fixed = TRUE),
                                ",", fixed = TRUE), as.numeric)
        out[grp] <- vals
      }
      out
    }
    b1_all <- per_prefix_values("b+1"); b2_all <- per_prefix_values("b+2")
    y1_all <- per_prefix_values("y+1"); y2_all <- per_prefix_values("y+2")
    # 生产顺序：b+1 / b+2 / y+1 / y+2 的 TYPE 行，各为逐位强度逗号列表（4×肽长个值，与 by_ions_long 同形）。"m/z intensity label" 峰行不用——数字相同但按 m/z 排序，会改变向量顺序、从而改变与 by_ions_long.match_intensity 的余弦。
    # production order: the b+1 / b+2 / y+1 / y+2 TYPE lines, each a comma
    # list of per-position intensities (4 x peptide-length values, same shape
    # as by_ions_long). The "m/z intensity label" peak lines are NOT used -
    # they carry the same numbers in m/z order, which changes the vector
    # order and therefore the cosine against by_ions_long.match_intensity.
    for (i in seq_len(n)) {
      b1v <- b1_all[[i]]; if (is.null(b1v)) b1v <- numeric(0)
      b2v <- b2_all[[i]]; if (is.null(b2v)) b2v <- numeric(0)
      y1v <- y1_all[[i]]; if (is.null(y1v)) y1v <- numeric(0)
      y2v <- y2_all[[i]]; if (is.null(y2v)) y2v <- numeric(0)
      # y 向量 reversepDeep2 的 y+1/y+2 向量按
      # "位点顺序"存储（index 1 = 最长 y 离子 y_{n-1}，index n-1 = y1），与
      # by_ions_long / b_y_ions_intensity 的 y1→y_{n-1} 正序相反——此前直接
      # 拼接未 reverse，导致 predict 余弦的 y 离子强度逐位错位（理论注入谱上
      # 实测定位：MFRPAQR y4 实测 0.992 被分配到 y3 位置）。
      # y-vector reverse pDeep2 stores the
      # y+1/y+2 vectors in position order (index 1 = longest y ion y_{n-1},
      # index n-1 = y1), the opposite of the y1->y_{n-1} forward order used by
      # by_ions_long / b_y_ions_intensity - the raw concatenation below used
      # to leave every y-ion intensity at the wrong position in the predict
      # cosine (measured on the synthetic spectra: MFRPAQR's y4 intensity
      # 0.992 was assigned to the y3 position).
      intensities[[i]] <- c(b1v, b2v, rev(y1v), rev(y2v))
    }
  }
  list(keys = keys, intensities = intensities, idx = setNames(seq_len(n), keys))
}
#' 写入生成参数文件的绝对路径。
#' Absolute path of 'p' for writing into a generated parameter file.
#'
#' 软件生成的文件一律绝对路径（与用户填写规则相同）；相对 "../" 引用在
#' 带点的中文目录名下不可靠。
#' The software-generated files always carry ABSOLUTE paths (same rule as the
#' user-written ones); relative "../" references proved unreliable under the
#' dotted Chinese folder names.
#'
#' @param p
#'   字符串，任意路径。
#'   Character string, any path.
#' @param out_path
#'   字符串，生成文件的路径（仅用于语义对应，不影响结果）。
#'   Character string, the generated file's path (semantic only).
#' @return
#'   返回规范化的绝对路径。
#'   The normalized absolute path.
#' @keywords internal
rel_to <- function(p, out_path) {
  normalizePath(p, winslash = "/", mustWork = TRUE)
}
#' 按标题 raw 前缀把输入无标记 mgf 重切成逐 raw 纯文件。
#' Re-split the input unlabeled mgfs into pure per-raw files by title prefix.
#'
#' 写出 unlabeled\<raw>_HCDFT.mgf（覆盖旧的按输入文件版本）。输入 MGF
#' （pParse 共洗脱拆分）与旧的逐输入文件输出可能 raw 混杂；pUniFind rescore 按
#' qry.res 内的 raw 名加载 <raw>_HCDFT.mgf 并以其中 TITLE 建字典，要求 mgf 与 raw
#' 严格一一对应（否则 KeyError）。流式处理：只缓存当前块，逐 raw 追加写。
#' Writes unlabeled\<raw>_HCDFT.mgf (overwriting the old per-input-file
#' versions). The input MGFs (pParse co-elution split) and the old per-file
#' outputs can mix raws; pUniFind rescore loads <raw>_HCDFT.mgf by the raw name
#' inside the qry.res and indexes its TITLEs, so the mgf must be raw-pure
#' (otherwise KeyError). Streaming: only the current block is buffered, each
#' raw is appended to its own file.
#'
#' @param mgf_files
#'   字符向量，输入无标记 mgf 文件（逐输入提取的中间产物，可多个）。
#'   A character vector of the input unlabeled mgfs (per-input intermediates).
#' @param out_dir
#'   字符串，输出目录（unlabeled\）。
#'   Character string, the output directory (unlabeled\).
#' @return
#'   返回写出的文件绝对路径字符向量。
#'   A character vector of the written absolute file paths.
#' @keywords internal
write_per_raw_mgfs <- function(mgf_files, out_dir) {
  cons <- list()
  written <- character(0)
  written_raws <- character(0)
  block <- character(1024)   # 预分配块缓冲（超出一谱行数时自动扩展）
  bidx <- 1L
  cur_raw <- ""
  flush_block <- function(raw) {
    if (bidx <= 1L || !nzchar(raw)) return(invisible(NULL))
    if (is.null(cons[[raw]])) {
      #：raw 数超过 R 连接上限（128）时批量关闭再按 append 重开——
      # 加标实验 511 个 raw 触发 "all 128 connections are in use" 崩溃的修复
      if (length(cons) >= 100L) {
        for (cc in names(cons)) try(close(cons[[cc]]), silent = TRUE)
        cons <<- list()
      }
      out_path <- file.path(out_dir, paste0(raw, "_HCDFT.mgf"))
      if (file.exists(out_path) && !(raw %in% written_raws)) unlink(out_path)
      cons[[raw]] <<- file(out_path, open = "ab")
      if (!(raw %in% written_raws)) {
        written <<- c(written, out_path)
        written_raws <<- c(written_raws, raw)
      }
    }
    writeLines(block[seq_len(bidx - 1L)], cons[[raw]])
    invisible(NULL)
  }
  for (f in mgf_files) {
    con_in <- file(f, open = "r")
    repeat {
      chunk <- readLines(con_in, n = 20000, warn = FALSE)
      if (length(chunk) == 0) break
      for (l in chunk) {
        if (l == "BEGIN IONS") {
          block[1L] <- "BEGIN IONS"
          bidx <- 2L
          cur_raw <- ""
        } else if (l == "END IONS") {
          block[bidx] <- "END IONS"
          bidx <- bidx + 1L
          flush_block(cur_raw)
          bidx <- 1L
        } else {
          block[bidx] <- l
          bidx <- bidx + 1L
          if (startsWith(l, "TITLE=")) {
            cur_raw <- strsplit(sub("^TITLE=", "", l), ".", fixed = TRUE)[[1]][1]
          }
        }
      }
    }
    close(con_in)
  }
  for (con in cons) try(close(con), silent = TRUE)
  written
}
#' 定位/创建本次运行的结果文件夹（pGlycoQuant 风格 result_YYYYMMDDHHMMSS）。
#' Locate/create the result folder for this run (pGlycoQuant-style
#' result_YYYYMMDDHHMMSS).
#'
#' 复用语义：项目下已存在 result_* 文件夹时取最新一个（断点续跑就地复用，
#' 各步跳过检查保持不变）；不存在时以当前时间戳新建。删除/移走全部
#' result_* 后重跑即得到新的时间戳文件夹，实现"结果与时间唯一性"——
#' 全部流程产出（清洗件、merge、1.unlabeled_result、2.labeled_result）都
#' 写在该文件夹内；搜索原始产物（pFind_open 等）与输入（unlabeled、config）
#' 保留在项目根，供断点复用。
#' Reuse semantics: when result_* folders already exist under the project the
#' newest one is used (resume reuses it in place; the per-step skip checks are
#' unchanged); otherwise a new one is created with the current timestamp.
#' Removing all result_* folders and rerunning yields a fresh timestamped
#' folder, i.e. result-time uniqueness - every pipeline output (cleaned copies,
#' merge, 1.unlabeled_result, 2.labeled_result) lives inside it, while the raw
#' search outputs (pFind_open etc.) and the inputs (unlabeled, config) stay at
#' the project root for resume reuse.
#'
#' @param cfg_dir
#'   字符标量，配置目录。
#'   A character scalar, the config folder.
#' @return
#'   字符标量，结果文件夹路径。
#'   A character scalar, the result folder path.
#' @keywords internal
resolve_result_dir <- function(cfg_dir) {
  if (!dir.exists(cfg_dir)) {
    dir.create(cfg_dir, recursive = TRUE, showWarnings = FALSE)
  }
  cands <- list.dirs(cfg_dir, recursive = FALSE, full.names = TRUE)
  cands <- cands[grepl("^result_[0-9]{14}$", basename(cands))]
  if (length(cands) > 0) {
    return(cands[order(basename(cands), decreasing = TRUE)[1]])
  }
  nd <- file.path(cfg_dir, paste0("result_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(nd, recursive = TRUE, showWarnings = FALSE)
  nd
}
#' 从环境变量 RESULT_DIR 或 resolve_result_dir 取结果文件夹。
#' Get the result folder from the RESULT_DIR environment variable or
#' resolve_result_dir.
#'
#' 总程序把 RESULT_DIR 注入各步子进程，保证同一文件夹贯穿全流程；
#' 脚本独立运行时退化为按时间戳定位/新建。
#' The master injects RESULT_DIR into every step subprocess so one folder spans
#' the whole run; standalone script runs fall back to locate/create by
#' timestamp.
#'
#' @param cfg_dir
#'   字符标量，配置目录。
#'   A character scalar, the config folder.
#' @return
#'   字符标量，结果文件夹路径。
#'   A character scalar, the result folder path.
#' @keywords internal
get_result_dir <- function(cfg_dir) {
  rd <- Sys.getenv("RESULT_DIR")
  if (nzchar(rd)) rd else resolve_result_dir(cfg_dir)
}
#' 对一张搜库结果表执行诱饵库去除。
#' Decoy removal on ONE search-result table .
#'
#' 逐行处理 Proteins 列（pFind 格式：条目以 "/" 结尾、直接相邻拼接，如
#' "XP_.../REV_.../"）：列表全部为 REV_ 诱饵的行整行删除；正库与诱饵库并存
#' 的行仅保留正库蛋白，Positions（"pos,pre,post/" 拼接）按相同下标同步剔除。
#' 无 Proteins 列的转换件（pUniFind rescore/denovo，转换时已按
#' pscore_is_target 过滤）在提供 annot（File_Name+Sequence→Proteins 注释，
#' 通常来自同一项目的未过滤 pFind.spectra/qry.res）时按同一规则删除全诱饵行
#'（保留行维持原列布局，不新增蛋白列）；注释缺失的行在提供 fa（正库蛋白组
#' 序列向量）时做兜底：序列不在正库则删除。清洗件以 "<tag>.spectra" 写入
#' cleaned_dir，返回 list(file=, n_in=, n_out=, n_dropped=, n_stripped=,
#' dropped_spectra=)。
#' Processes the Proteins column row by row (pFind format: entries are
#' "/"-terminated and directly concatenated, e.g. "XP_.../REV_.../"): rows
#' whose protein list is all REV_ decoy are dropped entirely; rows holding
#' both target and decoy proteins keep only the target proteins, with
#' Positions ("pos,pre,post/" concatenated) stripped at the same subscripts.
#' Converted files without a Proteins column (pUniFind rescore/denovo,
#' already filtered by pscore_is_target at conversion) drop their all-decoy
#' rows by the same rule when annot (File_Name+Sequence -> Proteins, usually
#' from the same project's unfiltered pFind.spectra/qry.res) is supplied
#' (kept rows keep the original column layout, no protein column is added);
#' annot-missing rows fall back to fa (forward-proteome sequence vector) when
#' provided: dropped if the sequence is absent from the forward proteome.
#' The cleaned copy is written to cleaned_dir as "<tag>.spectra"; returns
#' list(file=, n_in=, n_out=, n_dropped=, n_stripped=, dropped_spectra=).
#'
#' @param x
#'   字符标量，.spectra 文件路径。
#'   A character scalar, the .spectra file path.
#' @param cleaned_dir
#'   字符标量，清洗件输出目录。
#'   A character scalar, the output folder for cleaned copies.
#' @param tag
#'   字符标量，输出文件前缀（不含扩展名）。
#'   A character scalar, the output file prefix (without extension).
#' @param annot
#'   可选数据框，含 File_Name/Sequence/Proteins 三列；无 Proteins 列的输入
#'   用它判定诱饵。
#'   Optional data frame with File_Name/Sequence/Proteins; used to judge
#'   decoys when the input lacks a Proteins column.
#' @param fa
#'   可选字符向量，正库蛋白组序列；注释缺失行的兜底判定。
#'   Optional character vector of forward-proteome sequences; fallback for
#'   rows without annotation.
#' @return
#'   返回清洗统计列表。
#'   The cleaning statistics list.
#' @keywords internal
clean_decoy_spectra_base <- function(x, cleaned_dir, tag, annot = NULL, fa = NULL) {
  dat <- read.delim(x, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = "character", na.strings = character(0))
  n_in <- nrow(dat)
  out <- file.path(cleaned_dir, paste0(tag, ".spectra"))
  # 无 Proteins 列：有注释则按注释做同一套全诱饵删除（保留行维持原列布局）；
  # 无注释则原样复制（诱饵过滤已在转换时完成，仅留档）
  # no Proteins column: with annotation apply the same all-decoy drop (kept
  # rows keep the original column layout); without annotation copy as-is
  # (decoy filtering already done at conversion, just for the record)
  if (!"Proteins" %in% colnames(dat)) {
    if (!is.null(annot) && nrow(annot) > 0 &&
        all(c("File_Name", "Sequence", "Proteins") %in% colnames(annot))) {
      m <- match(paste(dat$File_Name, dat$Sequence),
                 paste(annot$File_Name, annot$Sequence))
      ap <- rep("", n_in)
      ap[!is.na(m)] <- annot$Proteins[m[!is.na(m)]]
      parts <- strsplit(ap, "/", fixed = TRUE)
      keep <- rep(TRUE, n_in)
      n_stripped <- 0L
      # v1.0 优化原逐行 trimws+grepl 循环（57k×5 文件约
      # 20.9s）改为展平一次向量化——flat：全部蛋白片段；nz：非空标记；g：
      # 片段所属行号；nrev/ntot：每行诱饵数/非空条目数。注意 strsplit 对
      # 空串返回 character(0)，空行会缺席 rowsum 分组——按行名回填全长向量
      if (n_in > 0) {
        flat <- trimws(unlist(parts, use.names = FALSE))
        nz <- nzchar(flat)
        isrev <- nz & grepl("^REV_", flat)
        nrev <- integer(n_in); ntot <- integer(n_in)
        if (length(flat) > 0) {
          g <- rep(seq_len(n_in), lengths(parts))
          grp_ids <- as.integer(names(rowsum(rep(1L, length(g)), g)[, 1L]))
          nrev[grp_ids] <- rowsum(as.integer(isrev), g)[, 1L]
          ntot[grp_ids] <- rowsum(as.integer(nz), g)[, 1L]
        }
        keep <- !(ntot > 0L & nrev == ntot)
        n_stripped <- sum(nrev > 0L & nrev < ntot)
        # 注释缺失（无蛋白条目）：正库兜底——序列不在正库蛋白组则删除；
        # 序列先去重再逐个判（同序列结果一致），避免逐行扫全蛋白组
        empty_rows <- which(ntot == 0L)
        if (length(empty_rows) > 0 && !is.null(fa) && length(fa) > 0) {
          uniq_seq <- unique(dat$Sequence[empty_rows])
          in_fa <- vapply(uniq_seq, function(s) any(grepl(s, fa, fixed = TRUE)),
                          logical(1))
          keep[empty_rows] <- in_fa[match(dat$Sequence[empty_rows], uniq_seq)]
        }
      }
      dropped <- unique(dat$File_Name[!keep])
      dat <- dat[keep,, drop = FALSE]
      write.table(dat, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      return(list(file = out, n_in = n_in, n_out = nrow(dat),
                  n_dropped = n_in - nrow(dat), n_stripped = n_stripped,
                  dropped_spectra = dropped))
    }
    write.table(dat, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    return(list(file = out, n_in = n_in, n_out = n_in, n_dropped = 0L,
                n_stripped = 0L, dropped_spectra = character(0)))
  }
  # pFind 的 Proteins/Positions 均为 "/" 结尾条目直接拼接（实测无 ";" 分隔符），
  # 因此按 "/" 切分并丢弃空串；重建时条目加回结尾 "/"
  # pFind's Proteins/Positions are "/"-terminated entries concatenated
  # directly (no ";" separator in practice), so split on "/" and drop empty
  # strings; the trailing "/" is restored on rebuild
  pr <- dat$Proteins
  parts <- strsplit(pr, "/", fixed = TRUE)
  keep <- rep(TRUE, n_in)
  n_stripped <- 0L
  new_pr <- pr
  new_po <- if ("Positions" %in% colnames(dat)) dat$Positions else NULL
  po_parts <- if (!is.null(new_po)) strsplit(new_po, "/", fixed = TRUE) else NULL
  # v1.0 优化同上向量化——只对"正库+诱饵并存"的剥离行
  # 用 vapply 重建 Proteins/Positions（保持与旧循环逐字节同义：片段先去空、
  # 再按 isrev 掩码剔除、结尾补 "/"）。空行（strsplit 返回 character(0)）
  # 按行名回填全长 nrev/ntot，避免 keep 向量长度错位被下标回收
  if (n_in > 0) {
    flat <- trimws(unlist(parts, use.names = FALSE))
    nz <- nzchar(flat)
    isrev <- nz & grepl("^REV_", flat)
    nrev <- integer(n_in); ntot <- integer(n_in)
    if (length(flat) > 0) {
      g <- rep(seq_len(n_in), lengths(parts))
      grp_ids <- as.integer(names(rowsum(rep(1L, length(g)), g)[, 1L]))
      nrev[grp_ids] <- rowsum(as.integer(isrev), g)[, 1L]
      ntot[grp_ids] <- rowsum(as.integer(nz), g)[, 1L]
    }
    keep <- !(ntot > 0L & nrev == ntot)
    strip_rows <- which(nrev > 0L & nrev < ntot)
    n_stripped <- length(strip_rows)
    if (length(strip_rows) > 0) {
      pos_of <- cumsum(lengths(parts)) - lengths(parts) + 1L
      #：长度向量预算一次——闭包内 lengths(parts)[i] 每次迭代对全表
      # list 重算 lengths，剥离行上万时二次方退化（open 件实测 +35 min）
      parts_len <- lengths(parts)
      new_pr[strip_rows] <- vapply(strip_rows, function(i) {
        sel <- pos_of[i] + seq_len(parts_len[i]) - 1L
        pnz <- nz[sel]
        paste0(paste(flat[sel][pnz][!isrev[sel][pnz]], collapse = "/"), "/")
      }, character(1))
      if (!is.null(po_parts)) {
        po_flat <- trimws(unlist(po_parts, use.names = FALSE))
        po_nz <- nzchar(po_flat)
        po_start <- cumsum(lengths(po_parts)) - lengths(po_parts) + 1L
        po_len <- lengths(po_parts)
        new_po[strip_rows] <- vapply(strip_rows, function(i) {
          sel <- pos_of[i] + seq_len(parts_len[i]) - 1L
          pnz <- nz[sel]
          po_sel <- po_start[i] + seq_len(po_len[i]) - 1L
          qnz <- po_nz[po_sel]
          if (sum(pnz) == sum(qnz)) {
            # positions 非空条目数与 proteins 一致才同步剔除（旧口径）
            paste0(paste(po_flat[po_sel][qnz][!isrev[sel][pnz]],
                         collapse = "/"), "/")
          } else new_po[i]
        }, character(1))
      }
    }
  }
  dropped <- unique(dat$File_Name[!keep])
  dat <- dat[keep,, drop = FALSE]
  dat$Proteins <- new_pr[keep]
  if (!is.null(new_po)) dat$Positions <- new_po[keep]
  write.table(dat, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  list(file = out, n_in = n_in, n_out = nrow(dat),
       n_dropped = n_in - nrow(dat), n_stripped = n_stripped,
       dropped_spectra = dropped)
}
# v1.0诱饵清洗行级过滤快速路径。
# 真实瓶颈不是诱饵判定（v1.0 已向量化），而是 read.delim 对 515 MB/
# 1.8M 行×18 列字符表的逐字段解析与 write.table 全量回写（实测整段
# 91–97 min ≈ 3 ms/行）。新路径：readLines 整读（~200 MB/s）→ 一次向量化
# grepl("REV_", fixed=TRUE) 预筛（诱饵判定只看 Proteins 片段前缀 REV_，
# 片段含 REV_ 则行内必现——预筛无漏检；他列偶含 "REV_" 只会多进慢枝，
# 判定后原样保留）→ 未命中行原样透写（保输入字节），仅命中行进入列解析
# 与既有 nrev/ntot 向量化判定与 Proteins/Positions 重建（口径与 _base 逐条
# 同义）。输出与旧路径唯一差别＝表头：旧 read.delim 会把 "Exp.MH+" 规范化
# 成 "Exp.MH."，新路径保留原始表头——下游读表时两者同样被 check.names
# 规范化，语义等价（decoy_clean_equiv_test 断言解析后表相等）。列数异常
# （防御）或无 Proteins 列（punifind 回查/兜底）回落 _base 旧实现。
# (perf round 1): line-level filtering for decoy cleaning. The real
# cost was read.delim parsing + write.table re-serializing every field of a
# 515 MB / 1.8M-row character table (~3 ms/row); the decoy logic itself has
# been vectorized since v1.0. New path: read the file as lines, one
# vectorized fixed-string prefilter grepl("REV_") (the verdict only looks at
# REV_ prefixes of Proteins fragments, so a hit line must contain "REV_" -
# no false negatives; occasional hits in other columns just take the parse
# branch and pass through unchanged), verbatim passthrough for non-hit lines,
# and the existing vectorized verdict/rebuild for hit lines only. The only
# output difference vs _base is the header (original "Exp.MH+" kept instead
# of read.delim's mangled "Exp.MH."; both mangle identically under
# check.names downstream). Malformed rows or missing Proteins column fall
# back to _base.
clean_decoy_spectra_fast <- function(x, cleaned_dir, tag) {
  out <- file.path(cleaned_dir, paste0(tag, ".spectra"))
  hdr <- readLines(x, n = 1L, warn = FALSE)
  hdr_f <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
  ipr <- match("Proteins", hdr_f)
  ipo <- match("Positions", hdr_f)
  ifn <- match("File_Name", hdr_f)
  lines <- readLines(x, warn = FALSE)
  body <- lines[-1L]
  n_in <- length(body)
  ncols <- length(hdr_f)
  # 参差行（真实 open 件实测 ~15% 短行：表头 19 列、6.5 万行 6 列、21.4 万行
  # 18 列）：短行按旧 read.delim(fill=TRUE) 补 NA→write.table na="" 的口径
  # 在写侧补尾空字段（字节等价）；超长行回落 _base（旧 read.table 遇超长
  # 即错，正常产物不会出现）。
  # ragged rows are real (~15% of an open file): short lines are padded with
  # trailing empty fields at write time, byte-equal to the old
  # read.delim(fill=TRUE) + write.table(na="") behaviour; over-long lines
  # fall back to _base.
  # 计数用定长字面 gsub（single-char fixed 扫描）——勿用 "[^\t]" 负类正则
  # （每字符一次匹配，1.8M 行分钟级）；亦勿依赖 gregexpr（构造 1.8M 匹配表）。
  # count tabs with a fixed single-char gsub scan; a negated character class
  # costs one regex match PER CHARACTER (minutes on 1.8M lines).
  n_tab <- nchar(body) - nchar(gsub("\t", "", body, fixed = TRUE))
  ragged <- any(n_tab != ncols - 1L)
  if (any(n_tab > ncols - 1L))
    return(list(fallback = TRUE))
  pre <- grepl("REV_", body, fixed = TRUE)
  hit <- which(pre)
  n_stripped <- 0L
  dropped <- character(0)
  res <- body
  keep_all <- rep(TRUE, n_in)
  if (length(hit) > 0) {
    flds <- strsplit(body[hit], "\t", fixed = TRUE)
    # 短行 Proteins/Positions 字段缺失按空处理（旧口径 NA→ntot=0→保留）
    pr <- vapply(flds, function(f) if (length(f) >= ipr) f[ipr] else "", character(1))
    parts <- strsplit(pr, "/", fixed = TRUE)
    flat <- trimws(unlist(parts, use.names = FALSE))
    nz <- nzchar(flat)
    isrev <- nz & grepl("^REV_", flat)
    nh <- length(hit)
    nrev <- integer(nh); ntot <- integer(nh)
    if (length(flat) > 0) {
      g <- rep(seq_len(nh), lengths(parts))
      grp_ids <- as.integer(names(rowsum(rep(1L, length(g)), g)[, 1L]))
      nrev[grp_ids] <- rowsum(as.integer(isrev), g)[, 1L]
      ntot[grp_ids] <- rowsum(as.integer(nz), g)[, 1L]
    }
    keep <- !(ntot > 0L & nrev == ntot)
    strip_rows <- which(nrev > 0L & nrev < ntot)
    n_stripped <- length(strip_rows)
    if (length(strip_rows) > 0) {
      pos_of <- cumsum(lengths(parts)) - lengths(parts) + 1L
      #：长度向量预算一次——闭包内 lengths(parts)[i] 每次迭代对全表
      # list 重算 lengths，剥离行上万时二次方退化（open 件实测 +35 min）
      parts_len <- lengths(parts)
      new_pr <- vapply(strip_rows, function(i) {
        sel <- pos_of[i] + seq_len(parts_len[i]) - 1L
        pnz <- nz[sel]
        paste0(paste(flat[sel][pnz][!isrev[sel][pnz]], collapse = "/"), "/")
      }, character(1))
      for (j in seq_along(strip_rows)) flds[[strip_rows[j]]][ipr] <- new_pr[j]
      if (!is.na(ipo)) {
        po <- vapply(flds, function(f) if (length(f) >= ipo) f[ipo] else "", character(1))
        po_parts <- strsplit(po, "/", fixed = TRUE)
        po_flat <- trimws(unlist(po_parts, use.names = FALSE))
        po_nz <- nzchar(po_flat)
        po_start <- cumsum(lengths(po_parts)) - lengths(po_parts) + 1L
        po_len <- lengths(po_parts)
        for (i in strip_rows) {
          sel <- pos_of[i] + seq_len(parts_len[i]) - 1L
          pnz <- nz[sel]
          po_sel <- po_start[i] + seq_len(po_len[i]) - 1L
          qnz <- po_nz[po_sel]
          if (sum(pnz) == sum(qnz)) {
            # positions 非空条目数与 proteins 一致才同步剔除（与 _base 同口径）
            flds[[i]][ipo] <- paste0(paste(po_flat[po_sel][qnz][!isrev[sel][pnz]],
                                           collapse = "/"), "/")
          }
        }
      }
    }
    if (!is.na(ifn) && any(!keep))
      dropped <- unique(vapply(flds[!keep], function(f) f[ifn], character(1)))
    res[hit[keep]] <- vapply(flds[keep], paste, character(1), collapse = "\t")
    keep_all[hit] <- keep
  }
  out_body <- res[keep_all]
  if (ragged) {
    oi <- which(keep_all)
    pad <- (ncols - 1L) - n_tab[oi]
    sp <- which(pad > 0L)
    if (length(sp) > 0L)
      out_body[sp] <- paste0(out_body[sp], strrep("\t", pad[sp]))
  }
  con <- file(out, open = "wb")
  # CRLF 行尾对齐旧 write.table 文本模式输出（字节级对拍前提）
  writeLines(c(hdr, out_body), con, sep = "\r\n", useBytes = TRUE)
  close(con)
  list(fallback = FALSE, file = out, n_in = n_in, n_out = length(out_body),
       n_dropped = n_in - length(out_body), n_stripped = n_stripped,
       dropped_spectra = dropped)
}

clean_decoy_spectra <- function(x, cleaned_dir, tag, annot = NULL, fa = NULL) {
  # v1.0 调度：含 Proteins 列且非注释回查场景走行级过滤快速路径，
  # 否则回落 _base（口径同义，测试 decoy_clean_equiv_test/fast 双断言）。
  if (is.null(annot)) {
    hdr1 <- tryCatch(readLines(x, n = 1L, warn = FALSE), error = function(e) "")
    if (length(hdr1) == 1L && !is.na(hdr1) &&
        "Proteins" %in% strsplit(hdr1, "\t", fixed = TRUE)[[1]]) {
      r <- clean_decoy_spectra_fast(x, cleaned_dir, tag)
      if (!isTRUE(r$fallback)) return(r)
    }
  }
  clean_decoy_spectra_base(x, cleaned_dir, tag, annot = annot, fa = fa)
}
#' 对一张已读入搜库结果表执行 SAP 提取的核心。
#' Core of the SAP extraction on ONE already-read search-result table.
#'
#' 不修改调用方的 dat；fasta 须已做 I→L / N→GG 归一；返回带 revised_* 列的
#' SAP 行；表内无 "->" 替换行时返回 NULL。
#' Does not modify the caller's dat; fasta must already be I->L / N->GG
#' normalized; returns the SAP rows with the revised_* columns added, or NULL
#' when the table has no "->" substitution rows.
#'
#' @param dat
#'   数据框，搜库结果表。
#'   A data frame, a search-result table.
#' @param fasta
#'   命名字符向量，已归一化的蛋白组（read_fasta 输出）。
#'   A named character vector, the normalized proteome (read_fasta output).
#' @return
#'   返回 SAP 行数据框或 NULL。
#'   The SAP-row data frame or NULL.
#' @keywords internal
extract_sap_rows_df <- function(dat, fasta) {
  aa3 <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
           "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
           "Tyr", "Val", "Xle")
  aa1 <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
           "F", "P", "S", "T", "W", "Y", "V", "I")
  pat <- "^[0-9]+,[A-Za-z]{3}->[A-Za-z]{3}\\[.[A-Za-z]?\\]$"
  if (nrow(dat) == 0) return(NULL)
  si <- which(grepl("Sequence", colnames(dat)))[1]
  mi <- which(grepl("Modification", colnames(dat)))[1]
  mod <- as.character(dat[[mi]])
  keep <- grepl("->", mod, fixed = TRUE)
  if (!any(keep)) return(NULL)
  dat <- dat[keep,, drop = FALSE]
  mod <- as.character(dat[[mi]])
  seq <- as.character(dat[[si]])

  # 生产 revise 校正（1.SAP_extract.R 使用的子集）
  # production revise corrections (the subset used by 1.SAP_extract.R)
  drop <- grepl("_13C", mod, fixed = TRUE)
  met <- which(grepl("0,Met-", mod, fixed = TRUE))
  for (j in met) {
    seq[j] <- substr(seq[j], 2, nchar(seq[j]))
    parts <- strsplit(mod[j], ";", fixed = TRUE)[[1]]
    parts <- parts[-1]
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) {
      mod[j] <- "0,Acetyl[AnyN-term];"
    } else {
      pos <- as.integer(sub(",.*", "", parts))
      rest <- sub("^[0-9]+,", "", parts)
      pos <- pos - 1L
      if (any(pos == 0L)) {
        drop[j] <- TRUE
      } else {
        mod[j] <- paste0("0,Acetyl[AnyN-term];",
                         paste(paste0(pos, ",", rest), collapse = ";"),
                         ";")
      }
    }
  }
  rarg <- which(mod == "0,Arg[AnyN-term];")
  for (j in rarg) {
    seq[j] <- paste0("R", seq[j])
    mod[j] <- ""
  }
  if (any(drop)) {
    dat <- dat[!drop,, drop = FALSE]
    mod <- mod[!drop]
    seq <- seq[!drop]
  }
  # 替换条目：每行恰一条 "pos,XXX->YYY[site]"
  # replacement entries: exactly one "pos,XXX->YYY[site]" per row
  n <- length(seq)
  rp <- rep("", n)
  rm2 <- rep("", n)
  refp <- rep("", n)
  for (j in seq_len(n)) {
    parts <- strsplit(mod[j], ";", fixed = TRUE)[[1]]
    m <- grepl(pat, parts)
    if (sum(m) != 1) next
    sp <- strsplit(parts[which(m)], ",", fixed = TRUE)[[1]]
    target <- sub(".*->", "", sp[2])          # "YYY[site]"
    one <- aa1[match(substr(target, 1, 3), aa3)]
    pos <- suppressWarnings(as.integer(sp[1]))
    if (length(one) != 1 || is.na(one) || is.na(pos) ||
        pos < 1 || pos > nchar(seq[j])) next
    substr(seq[j], pos, pos) <- one
    rp[j] <- seq[j]
    rm2[j] <- paste(parts[-which(m)], collapse = ";")
    # 把 I→L / N→GG 归一后的肽段与归一后的蛋白组匹配（肽段：先 I→L 再 N→GG；生产顺序）
    # match the I->L / N->GG normalized peptide against the normalized
    # proteome (peptide: I->L then N->GG; production order)
    ns <- gsub("N", "GG", gsub("I", "L", seq[j], fixed = TRUE), fixed = TRUE)
    hits <- grepl(ns, fasta, fixed = TRUE)
    if (any(hits)) refp[j] <- paste(sort(names(fasta)[hits]), collapse = "/")
  }
  keep2 <- rp != ""
  if (!any(keep2)) return(NULL)
  dat <- dat[keep2,, drop = FALSE]
  dat$revised_peptide <- rp[keep2]
  dat$revised_modification <- rm2[keep2]
  dat$ref_protein <- refp[keep2]
  dat$six_protein <- ""
  dat
}
#' 把逐文件 SAP 表统一列布局并行拼接。
#' Unify a list of per-file SAP tables to one column layout and row-bind.
#'
#' pFind 文件带全 19 列；3 列 pUniFind 文件保留 File_Name/Sequence/Modification、
#' 其余列填空。
#' pFind files carry the full 19 columns; 3-column pUniFind files keep
#' File_Name/Sequence/Modification and get empty columns for the rest.
#'
#' @param parts
#'   数据框列表，逐文件 SAP 表。
#'   A list of data frames, the per-file SAP tables.
#' @return
#'   返回统一列布局的行拼接数据框。
#'   The row-bound data frame with the unified column layout.
#' @keywords internal
bind_sap_rows <- function(parts) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) {
    return(data.frame(File_Name = character(0), stringsAsFactors = FALSE))
  }
  base_cols <- colnames(parts[[1]])
  for (k in seq_along(parts)) {
    miss <- setdiff(base_cols, colnames(parts[[k]]))
    if (length(miss) > 0) {
      for (mc in miss) parts[[k]][[mc]] <- ""
    }
    parts[[k]] <- parts[[k]][, base_cols, drop = FALSE]
  }
  extra <- setdiff(unique(unlist(lapply(parts, colnames))), base_cols)
  if (length(extra) > 0) {
    for (k in seq_along(parts)) {
      e2 <- setdiff(extra, colnames(parts[[k]]))
      for (ec in e2) parts[[k]][[ec]] <- ""
    }
    parts <- lapply(parts, function(d) d[, c(base_cols, extra), drop = FALSE])
  }
  do.call(rbind, parts)
}
#' pFind 安装目录规范化 → Searcher.exe 绝对路径。
#' Normalize the pFind install directory to the Searcher.exe absolute path.
#'
#' 以 bin 结尾则直接使用；若填的是版本目录（如 pFind3.1）则自动补 bin；
#' 不存在时报错。
#' Use it as-is when it ends with bin; append bin when a version folder (e.g.
#' pFind3.1) is given; error when Searcher.exe is not found.
#'
#' @param p
#'   字符串，pFind 安装路径或 Searcher.exe 路径。
#'   Character string, the pFind install path or the Searcher.exe path.
#' @return
#' 通用软件路径解析器：接受任意层级的用户输入，DFS 找到签名文件后
#' 返回正确的约定层级目录。用户填多深或多浅都能自愈。
#' Generic software path resolver: accepts any-level user input, DFS-finds
#' the signature file, and returns the directory at the correct convention
#' level. Heals both too-deep and too-shallow user inputs.
#'
#' @param user_path 用户填写的路径（可以是 exe 本身、正确目录、过深/过浅目录）
#' @param signature_file 要找的签名文件名（如 "Searcher.exe"、"predict.py"）
#' @param levels_up 从签名文件所在目录向上取几级作为约定目录（0=签名文件
#'   所在目录本身就是约定目录；正数=向上 N 级）
#' @param max_depth DFS 最大深度（默认 5 层）
#' @param max_nodes DFS 最大节点数（默认 3000，防超时）
#' @return 约定层级的目录路径；找不到时返回 ""（调用方决定后续）
#' @keywords internal
resolve_software_path <- function(user_path, signature_file,
                                  levels_up = 0L,
                                  max_depth = 5L, max_nodes = 3000L) {
  if (!nzchar(user_path)) return("")
  p <- gsub("\\\\", "/", trimws(user_path))
  # ① 用户给的就是签名文件本身
  if (tolower(basename(p)) == tolower(signature_file)) {
    hit_dir <- dirname(p)
  } else if (file.exists(file.path(p, signature_file))) {
    # ② 签名文件直接在路径下（正确层级或零层偏移）
    hit_dir <- p
  } else if (dir.exists(p)) {
    # ③ DFS 向下找（list 装路径＋深度，避免 c() 把深度变字符）
    hit_dir <- ""
    stack <- list(list(p, 0L)); n <- 0L
    while (length(stack) > 0 && n < max_nodes) {
      n <- n + 1L
      cur <- stack[[1]]; stack[[1]] <- NULL
      d <- cur[[2]]
      cand <- file.path(cur[[1]], signature_file)
      if (file.exists(cand)) { hit_dir <- cur[[1]]; break }
      if (d < max_depth) {
        subs <- list.dirs(cur[[1]], recursive = FALSE)
        for (s in rev(subs)) stack <- c(list(list(s, d + 1L)), stack)
      }
    }
    if (!nzchar(hit_dir)) return("")
  } else {
    return("")
  }
  # 从签名文件所在目录向上取 levels_up 级
  out <- hit_dir
  for (i in seq_len(levels_up)) out <- dirname(out)
  normalizePath(out, winslash = "/")
}

#' pDeep2 专用解析：不管用户填哪一级，自愈到含 pDeep/pDeep2/predict.py
#' 的根目录（predict.py 上 3 级）。同时验证 python 可达。
#' pDeep2-specific resolver: regardless of input level, heals to the root
#' directory containing pDeep/pDeep2/predict.py (3 levels up from predict.py).
#' Also verifies python is reachable.
#'
#' @keywords internal
resolve_pdeep2_exe <- function(user_path) {
  resolved <- resolve_software_path(user_path, "predict.py", levels_up = 3L)
  if (!nzchar(resolved)) return(user_path)  # 找不到时原样返回（让调用方报具体错）
  # 验证 predict.py 确实在 <root>/pDeep/pDeep2/predict.py
  if (file.exists(file.path(resolved, "pDeep", "pDeep2", "predict.py"))) {
    return(resolved)
  }
  # 附件包布局（predict.py 在 <addon>/pDeep2/pDeep/pDeep2/predict.py，
  # resolve 返回的 levels_up=3 可能指向 addon 根而非 pDeep2 子目录）
  cand <- file.path(resolved, "pDeep2")
  if (file.exists(file.path(cand, "pDeep", "pDeep2", "predict.py"))) {
    return(normalizePath(cand, winslash = "/"))
  }
  # 反向：如果 levels_up=3 给的是 pDeep2 根的父目录，再试一层下
  cand2 <- file.path(dirname(resolved), "pDeep2")
  if (file.exists(file.path(cand2, "pDeep", "pDeep2", "predict.py"))) {
    return(normalizePath(cand2, winslash = "/"))
  }
  resolved
}
normalize_pfind_exe <- function(p) {
  # v1.0约定层级＝pFind 版本目录（含
  # bin/Searcher.exe 的那一级，如 C:/pFindStudio/pFind3）。宽容化：bin
  # 目录、Searcher.exe/pFind.exe 全路径、只给到 Studio 根（自动扫一层
  # 版本目录）全部可接；既有输入的结果保持不变。
  p <- sub("[\\\\/]+$", "", p)
  if (grepl("\\.exe$", tolower(p))) {
    if (tolower(basename(p)) == "searcher.exe" && file.exists(p)) return(p)
    p <- dirname(p)
    if (tolower(basename(p)) != "bin") p <- file.path(p, "bin")
  }
  cands <- c(file.path(p, "bin", "Searcher.exe"),
             file.path(p, "Searcher.exe"))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) {
    for (d in list.dirs(p, recursive = FALSE)) {
      cand <- file.path(d, "bin", "Searcher.exe")
      if (file.exists(cand)) { hit <- cand; break }
    }
  }
  if (length(hit) == 0)
    stop("pFind Searcher.exe not found under: ", p,
         " (expected the pFind version folder, e.g. C:/pFindStudio/pFind3, containing bin/Searcher.exe)")
  hit[[1]]
}

#' pUniFind denovo 的 _merged.csv 转为三列 .spectra。
#' Convert a pUniFind denovo _merged.csv into a three-column .spectra table.
#'
#' 输出 File_Name / Sequence / Modification 三列；列名按常见标题自动识别
#' （title/spectrum、peptide/sequence、modif）。
#' Outputs File_Name / Sequence / Modification; columns are matched by common
#' names (title/spectrum, peptide/sequence, modif).
#'
#' @param csv_path
#'   字符串，_merged.csv 的路径。
#'   Character string, the _merged.csv path.
#' @param out_dir
#'   字符串，输出目录。
#'   Character string, the output directory.
#' @return
#'   返回写出的 .spectra 文件路径。
#'   The written .spectra file path.
#' @keywords internal
denovo_csv_to_spectra <- function(csv_path, out_dir) {
  # 兼容两种 pUniFind denovo CSV：①带表头的常见格式（title/peptide/modif 列）；②GUI 产出的
  # 无表头 9 列格式：列1 标题、列8 修订肽段（可能带 "_位点_修饰名" 后缀）、列9 修饰串（空则无）。
  # Accepts both pUniFind denovo CSV layouts: 1. a headed table with title/peptide/modif
  # columns; 2. the GUI's headerless 9-column format: col1 title, col8 revised peptide
  # (possibly with a "_pos_ModName" suffix), col9 the modification string (empty = none).
  first <- readLines(csv_path, n = 1, warn = FALSE)
  if (grepl("title|spectrum|peptide|modif", tolower(first))) {
    dat <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
    cn <- tolower(colnames(dat))
    i_fn <- which(grepl("title|spectrum", cn))[1]
    i_sq <- which(grepl("peptide|sequence", cn))[1]
    i_md <- which(grepl("modif", cn))[1]
    if (is.na(i_fn) || is.na(i_sq)) {
      stop("denovo csv columns not recognized: ", paste(colnames(dat), collapse = ","))
    }
    fn <- as.character(dat[[i_fn]])
    sq <- as.character(dat[[i_sq]])
    md <- if (!is.na(i_md)) as.character(dat[[i_md]]) else ""
  } else {
    dat <- read.delim(csv_path, sep = ",", header = FALSE, stringsAsFactors = FALSE,
                      na.strings = character(0), quote = "\"", check.names = FALSE)
    fn <- as.character(dat[[1]])
    sq <- sub("_\\d.*$", "", as.character(dat[[8]]))
    md <- if (ncol(dat) >= 9) as.character(dat[[9]]) else rep("", length(fn))
    md[is.na(md)] <- ""
  }
  out <- data.frame(File_Name = fn, Sequence = sq, Modification = md,
                    stringsAsFactors = FALSE)
  # 去重denovo merged csv 对同一谱可有多行（并列候选），
  # 合并步骤按标题对齐会因重复行报行数不一致——每谱保留首个（top）候选
  # dedup the denovo merged csv may hold several rows for one
  # spectrum (tied candidates), which breaks the title-aligned merge with a
  # row-count mismatch - keep the first (top) candidate per spectrum
  out <- out[!duplicated(out$File_Name),, drop = FALSE]
  out_path <- file.path(out_dir, "pUniFind_denovo.spectra")
  write.table(out, out_path, sep = "	", row.names = FALSE, quote = FALSE)
  out_path
}

#' 按 pUniFind 兼容修饰集合过滤开放式搜库结果表。
#' Filter an open-search result table against the pUniFind-compatible modification set.
#'
#' pFind 开放式搜库若以全量 modification.ini 为修饰表，会解析出大量 pUniFind rescore
#' 不认识的修饰名（如 Methyl[D]）并导致 KeyError；生产 open 搜库以受限修饰表（1.mod）
#' 运行、其结果经 pUniFind 验证兼容。本函数逐单元格提取修饰名并与兼容集合（由生产
#' open 结果提取、随客户端分发）比对，删除含集合外修饰名的行；原文件就地覆盖、表头保留。
#' When the pFind open search uses the full modification.ini it reports many
#' modification names that pUniFind rescore does not know (e.g. Methyl[D]), raising
#' a KeyError; the production open search ran with a restricted table (1.mod) whose
#' results are pUniFind-verified. This function extracts the modification names cell
#' by cell, compares them with the compatible set (built from the production open
#' results and shipped with the client) and drops the rows carrying outside names;
#' the file is overwritten in place and the header is preserved.
#'
#' @param path
#'   字符串，.spectra 表路径（制表符分隔、首行表头）。
#'   Character string, the .spectra table path (tab-separated, header row).
#' @param compat_file
#'   字符串，兼容修饰集合文件（每行一个修饰名）。
#'   Character string, the compatible-modification set file (one name per line).
#' @return
#'   返回 c(rows=原数据行数, kept=保留行数, dropped=删除行数) 命名向量。
#'   A named vector c(rows, kept, dropped) with the original, kept and dropped counts.
#' @keywords internal
filter_mods_compatible <- function(path, compat_file) {
  compat <- readLines(compat_file, warn = FALSE)
  compat <- compat[nzchar(trimws(compat))]
  dat <- readLines(path, warn = FALSE)
  header <- dat[1]
  body <- dat[-1]
  # （性能第二轮）：原逐行×逐单元格×逐分片嵌套循环改展平一次判定
  # （rows×cells×pieces 三维展成一次子集运算，再映射回行）
  cells <- strsplit(body, "\t", fixed = TRUE)
  flat <- unlist(cells, use.names = FALSE)          # 全部单元格
  pieces <- strsplit(flat, ";", fixed = TRUE)
  flat2 <- unlist(pieces, use.names = FALSE)        # 全部修饰分片
  has_bracket <- grepl("[", flat2, fixed = TRUE)
  nm <- sub("^[^,]*?,", "", flat2)
  bad_piece <- has_bracket & nzchar(nm) & !(nm %in% compat)
  # piece → cell → row 两级映射（strsplit 的分组顺序即首现顺序）
  cell_row <- rep(seq_along(cells), lengths(cells))
  piece_row <- cell_row[rep(seq_along(flat), lengths(pieces))]
  bad_rows <- unique(piece_row[bad_piece])
  row_ok <- !(seq_along(body) %in% bad_rows)
  writeLines(c(header, body[row_ok]), path)
  c(rows = length(body), kept = sum(row_ok), dropped = sum(!row_ok))
}
#' 确保 pFind 搜库模板存在：约定目录缺模板时从内置模板自动生成。
#' Ensure the pFind search template exists: auto-generate it from the built-in
#' template when the convention folder is missing.
#'
#' 模板约定目录为 <project_dir>/pFind_open/param/pFind.cfg（open）与
#' <project_dir>/pFind_restrict/param/pFind.cfg（restrict）；已存在则直接复用
#' （用户可自行定制），不存在时从客户端分发的内置模板（code/pFind_templates/）
#' 生成：替换 fasta 路径、输出目录与输出名，以 GBK 编码写出（pFind 按 GBK 解析）。
#' 路径一律正斜杠，规避 pFind 对反斜杠转义序列（如 \1、\r）的解析。
#' The convention folders are <project_dir>/pFind_open/param/pFind.cfg (open) and
#' <project_dir>/pFind_restrict/param/pFind.cfg (restrict); an existing template is
#' reused as-is (user-customizable), otherwise it is generated from the shipped
#' built-in template (code/pFind_templates/) with the fasta path, output folder and
#' output name substituted, written as GBK. All paths use forward slashes to avoid
#' pFind's backslash-escape parsing (e.g. \1, \r).
#'
#' @param kind
#'   字符串，"open" 或 "restrict"。
#'   Character string, "open" or "restrict".
#' @param fasta_path
#'   字符串，参考蛋白库 FASTA 的绝对路径。
#'   Character string, the absolute path of the reference FASTA.
#' @param project_dir
#'   字符串，项目目录（通常为参数文件的 config_dir）。
#'   Character string, the project folder (usually the config folder).
#' @return
#'   返回模板 pFind.cfg 的路径。
#'   The path of the template pFind.cfg.
#' @keywords internal
ensure_pfind_template <- function(kind, fasta_path, project_dir) {
  #：三臂模板目录——open（开放式）/ restrict（半酶切多参数限定式）/
  # standard（常规参数限定式：全特异酶切＋闭合质量＋常规修饰组合）
  #: three arm templates - open / restrict (semi-specific multi-mod)
  # / standard (conventional closed search, fully-specific digestion)
  km <- c(open = "pFind_open", restrict = "pFind_restrict",
          standard = "pFind_standard")
  tm <- c(open = "open_pFind.cfg", restrict = "restrict_pFind.cfg",
          standard = "standard_pFind.cfg")
  if (!(kind %in% names(km)))
    stop("未知 pFind 模板类型：", kind, "\nunknown pFind template kind: ", kind)
  proj_name <- km[[kind]]
  dst <- file.path(project_dir, proj_name, "param", "pFind.cfg")
  if (file.exists(dst)) return(dst)
  tpl_rel <- file.path("pFind_templates", tm[[kind]])
  # P3模板路径不再依赖 cwd——先按全局 code_dir/
  # sd0 解析（脚本经 --file= 定位的代码目录），都不可用才回退 cwd 相对（旧行为）
  tpl <- character(0)
  for (base in c(if (exists("code_dir", envir = .GlobalEnv, inherits = FALSE))
                   get("code_dir", envir = .GlobalEnv) else "",
                 if (exists("sd0", envir = .GlobalEnv, inherits = FALSE))
                   get("sd0", envir = .GlobalEnv) else "")) {
    if (nzchar(base)) {
      cand <- file.path(base, tpl_rel)
      if (file.exists(cand)) { tpl <- cand; break }
    }
  }
  if (!length(tpl)) tpl <- tpl_rel
  if (!file.exists(tpl)) {
    stop("内置 pFind 模板缺失：", tpl,
         "\nbuilt-in pFind template missing: ", tpl)
  }
  out_dir <- file.path(project_dir, proj_name, "result")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lines <- readLines(tpl, warn = FALSE)
  lines <- gsub("__FASTAPATH__", fasta_path, lines, fixed = TRUE)
  lines <- gsub("__OUTPUTPATH__", paste0(out_dir, "/"), lines, fixed = TRUE)
  lines <- gsub("__OUTPUTNAME__", proj_name, lines, fixed = TRUE)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  con <- file(dst, open = "w", encoding = "GBK")
  writeLines(lines, con)
  close(con)
  cat("已自动生成 pFind 模板：", dst, "\n")
  cat("auto-generated pFind template: ", dst, "\n")
  dst
}

#' 用 pUniFind 官方命令行完成 rescore 或 denovo（参数与 checkpoint 对齐）。
#' Run the pUniFind official CLI for rescore or denovo (args aligned with the checkpoint).
#'
#' 在 pUniFind 安装的 punifind 目录运行官方 exe，参数集与 checkpoint_rank.pt 的
#' 训练参数对齐（use-clip=1、num-layers=9 等）；rescore 产物写回 qry_res_dir
#' （<项目名>fdr0.01_pUniFind.spectra），denovo 产物写到 out_dir/pUniFind_result
#' （test_5_0.csv）。命令行直接执行失败时，若路径含非 ASCII 字符，则在
#' C:\\temp 下新建英文临时目录、把 qry.res/mgf 拷入（大文件优先符号链接）后
#' 重试一次，结束后删除临时目录。
#' Runs the official exe from the pUniFind punifind folder with an argument set
#' aligned to the checkpoint_rank.pt training args (use-clip=1, num-layers=9, ...);
#' the rescore output is written back into qry_res_dir (<project>fdr0.01_pUniFind
#' .spectra) and the denovo output into out_dir/pUniFind_result (test_5_0.csv).
#' When the direct run fails AND the paths contain non-ASCII characters, a fresh
#' English-only temp folder is created under C:\\temp, the qry.res/mgf files are
#' copied in (large files via symlink first) and the run is retried once; the
#' temp folder is removed afterwards.
#'
#' @param kind
#'   字符串，"rescore" 或 "denovo"。
#'   Character string, "rescore" or "denovo".
#' @param punifind_dir
#'   字符串，pUniFind 安装的 punifind 目录（含 official_*.exe、ckpts、consts）。
#'   Character string, the pUniFind punifind folder (official_*.exe, ckpts, consts).
#' @param qry_res_dir
#'   字符串，qry.res 所在目录（rescore 输入；rescore 产物也写回此处）。
#'   Character string, the qry.res folder (rescore input; its output goes back here).
#' @param mgf_dir
#'   字符串，逐 raw MGF 目录（unlabeled）。
#'   Character string, the per-raw MGF folder (unlabeled).
#' @param out_dir
#'   字符串，输出目录（denovo 的 test_5_0.csv 写到其 pUniFind_result 子目录）。
#'   Character string, the output folder (denovo writes test_5_0.csv into its pUniFind_result subfolder).
#' @return
#'   返回 c(rescore_spectra=路径) 或 c(denovo_csv=路径)。
#'   A named vector: c(rescore_spectra=path) or c(denovo_csv=path).
#' @keywords internal
run_punifind_cli <- function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir,
                             converter_py = NULL, python_exe = NULL) {
  # N-16回退实现复用 pf_build_args——适配模块随包分发，
  # 模块文件存在时直接委托 pf_run（参数面单一真源：数据目录前置、--flag=value
  # 连写、flag 面探测、消毒/压实/新鲜度/缓存失效全链路），不再维护本函数内
  # 漂移的独立参数面；仅模块文件缺失时走下方 legacy 路径（仅修 worker flag
  # 漂移，不再扩展新参数）。
  code_dir_here <- if (exists("code_dir", envir = .GlobalEnv, inherits = FALSE))
    get("code_dir", envir = .GlobalEnv) else ""
  mod_file <- if (nzchar(code_dir_here))
    file.path(code_dir_here, "punifind_gpu_compat.R") else ""
  if (nzchar(mod_file) && file.exists(mod_file)) {
    mod_env <- new.env(parent = baseenv())
    if (exists("pf_run", mode = "function", inherits = TRUE)) {
      mod_env$pf_run <- pf_run
    } else {
      source(mod_file, local = mod_env, encoding = "UTF-8")
    }
    return(mod_env$pf_run(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir,
                          converter_py = converter_py, python_exe = python_exe))
  }
  # ---- legacy 回退（仅模块文件缺失时可达） ----
  exe <- file.path(punifind_dir, if (kind == "rescore") "official_score_workflow.exe"
                   else "official_denovo_workflow.exe")
  ckpt <- file.path(punifind_dir, "ckpts", "checkpoint_rank.pt")
  tok  <- file.path(punifind_dir, "consts", "tokenize1.pkl")
  meta <- file.path(punifind_dir, "consts", "modification_meta_dict.pkl")
  for (f in c(exe, ckpt, tok, meta)) {
    if (!file.exists(f)) stop("pUniFind 组件缺失：", f, "\npUniFind component missing: ", f)
  }
  task <- if (kind == "rescore") "score_inference" else "denovo_inference"
  data_tmp <- file.path(out_dir, "pUniFind_data_tmp")
  res_dir  <- file.path(out_dir, "pUniFind_result")
  dir.create(data_tmp, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  # pUniFind 推理的 CPU 并行度（设计约定：GUI Parameters 的
  # CPU Processors＝8，与 GPU Batch Size 64 同面板；此前 CLI 写死 --num-proc 1
  # 导致推理近单核、rescore/denovo 各约 3 h）。默认 8，可用环境变量
  # PUNIFIND_NUM_PROC 覆盖。
  # CPU parallelism of pUniFind inference (user ruling: the GUI
  # Parameters panel sets CPU Processors = 8 alongside GPU Batch Size 64; the
  # CLI had hardcoded --num-proc 1, keeping inference near single-core at
  # ~3 h per stage). Default 8; override with the PUNIFIND_NUM_PROC env var.
  num_proc <- Sys.getenv("PUNIFIND_NUM_PROC")
  if (!nzchar(num_proc) || suppressWarnings(is.na(as.integer(num_proc))) ||
      as.integer(num_proc) < 1) num_proc <- "8"
    # P2-1--num-workers 0 让 GPU 挨饿（8 分钟只
    # 推进 100/822 batch），pUniFind 自身默认 16。PUNIFIND_NUM_WORKERS 环境变量
    # 可配，默认 8；留空则省略该 flag、交给 pUniFind 默认值。
    num_workers <- trimws(Sys.getenv("PUNIFIND_NUM_WORKERS", unset = ""))

  mk_args <- function(q, m, o) {
    # system2 的 args 每个元素一个参数（含空格的值由 system2 自行加引号）
    # each element of args is ONE argument (values with spaces are quoted by system2)
    # P1-1新 pUniFind 构建不再定义 --load-contra，
    # "flag 与值分离"写法下其值 0 泄漏成位置参数 data → 0\test.parquet 找不到、
    # rescore/denovo 100% 失败。改法：① 数据目录放第一条（官方 CLI 同款）；
    # ② 删除 --load-contra；③ 数值型参数一律 --flag=value 连写（防同类泄漏）。
    # N-160 不再传出（旧口径 >=0 会把 0 传给 pUniFind），
    # 且改 --flag=value 连写（与 pf_build_args 同口径，防值泄漏成位置参数）
    nw <- if (nzchar(num_workers) && !suppressWarnings(is.na(as.integer(num_workers))) &&
              as.integer(num_workers) > 0) paste0("--num-workers=", num_workers) else character(0)
    c(file.path(o, "pUniFind_data_tmp"),
      "--user-dir", file.path(o, "pUniFind_data_tmp"),
      "--valid-subset", "test",
      "--results-path", file.path(o, "pUniFind_result"),
      "--tmp-data-path", file.path(o, "pUniFind_data_tmp", "test.lmdb"),
      "--project-name", "test",
      "--weight-path", ckpt,
      "--qry-res-path", q,
      "--mgf-path", m,
      paste0("--num-proc=", num_proc), "--data-buffer-size=256", nw,
      "--ddp-backend=c10d", "--batch-size=64",
      "--task", task, "--loss", "triplet", "--arch", "pep_ms",
      "--node-dim=512", "--num-layers=9", "--num-head=8",
      "--use-nce=1", "--use-ins=1", "--use-clip=1", "--use-hard-clip=0",
      "--denovo-pred=1", "--spectrum-pred=1", "--denoise-pred=0",
      "--intensity-pred=0", "--dataset-name", "_", "--max-charges=4",
      "--use-rope=0", "--noise-peak-pred=0", "--half-spectrum-pred=0",
      "--from-scratch=1", "--norm-first=0", "--modification-pred=1",
      "--modification-cls-pred=0", "--multi-class-noise-peak-pred=1",
      "--res-num-pred=0", "--res-type-pred=0", "--res-seq-len-pred=1",
      "--cutoff-spectra=300", "--comp-res-num-pred=0", "--train-pfind=1",
      "--hard-neg=1", "--shift=1", "--head-clip=1",
      "--fdr-thread=0.1",
      "--pkl-path", file.path(o, "pUniFind_data_tmp"),
      "--joint-encoding-layers=4", "--num-pept-layers=9",
      "--inference=1", "--check=0", "--joint-encoding=1",
      "--range-pred=5", "--tof=0",
      "--tokenize1-pkl-path", tok,
      "--modification-meta-dict-path", meta,
      "--no-progress-bar",
      # 小输入集（如只跑 SPK 命中谱的定向 denovo）不满足默认的批大小整除
      # 校验时，官方迭代器会除零崩溃——跳过该校验（fairseq 标准开关）。
      # small inputs (e.g. the targeted denovo over just the SPK-hit spectra)
      # trip a divide-by-zero in the official iterator when they fail the
      # default batch-size divisibility check - skip that check (standard
      # fairseq switch).
      "--skip-invalid-size-inputs-valid-test")
  }

  run_once <- function(q, m, o) {
    args <- mk_args(q, m, o)
    # N-16统一经 iso_system2 全量 shQuote（旧写法只对"含空格且
    # 不以 -- 开头"的条件加引号，含 &/(/)/^/% 的路径仍会被 cmd 截断；shQuote 对
    # 无特殊字符的参数原样返回，flag 不受影响）
    logf <- file.path(o, if (kind == "rescore") "rescore.log" else "denovo.log")
    st <- iso_system2(exe, args, stdout = logf, stderr = logf,
                      wait = TRUE)
    logtxt <- if (file.exists(logf)) {
      paste(readLines(logf, warn = FALSE, n = -1L), collapse = "\n")
    } else ""
    list(ok = (st == 0), status = st, log = logtxt)
  }

  # 消毒 qry.respFind 开放式搜库对修饰表外的 delta 质量写出
  # ";name#0" 占位修饰名（新加入的边缘谱所致），pUniFind 预处理的
  # multiprocessing 断言崩溃。这些 PSM 的唯一修饰即未命名 delta，无法被
  # pUniFind 评分——整行删除（token 清除会留下空修饰字段继续触发断言）。
  # sanitize the qry.res in place pFind's open search writes a
  # ";name#0" placeholder for delta masses outside the mod table (caused by the
  # marginal spectra added to the new unlabeled set), which trips pUniFind's
  # preprocessing assertion. Such PSMs carry nothing but the unnamed delta and
  # cannot be scored by pUniFind - the whole ROW is dropped (stripping just the
  # token leaves an empty modification field that keeps tripping the assertion).
  pf_sanitize_qryres(qry_res_dir)   # P0-A消毒提为公共函数，与 pf_run 共用

  r <- run_once(qry_res_dir, mgf_dir, out_dir)
  if (!r$ok && any(grepl("[^\x01-\x7F]", c(qry_res_dir, mgf_dir, out_dir)))) {
    # 路径含非 ASCII 且直接运行失败：英文临时目录兜底重试一次
    # non-ASCII paths failed directly: retry once from an English-only temp folder
    cat("pUniFind 直接运行失败，改用英文临时目录重试（C:\\temp）\n")
    cat("pUniFind direct run failed, retrying from an English temp folder (C:\\temp)\n")
    base_t <- if (dir.exists("C:/temp") || dir.create("C:/temp", showWarnings = FALSE))
      "C:/temp" else tempdir()                               # C:/temp 优先（纯 ASCII），
    tdir <- file.path(base_t, paste0("punifind_",            # 不可用退 R tempdir()
      format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(tdir, recursive = TRUE), add = TRUE)
    tq <- file.path(tdir, "qry"); dir.create(tq, recursive = TRUE)
    tm <- file.path(tdir, "mgf"); dir.create(tm, recursive = TRUE)
    if (kind == "rescore") {
      qrs <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)
      if (length(qrs) == 0) stop("rescore 输入 qry.res 缺失：", qry_res_dir)
      file.copy(qrs, tq)
    }
    mgfs <- list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE)
    for (f in mgfs) {
      dst <- file.path(tm, basename(f))
      if (!file.exists(dst)) {
        ok <- suppressWarnings(try(file.symlink(f, dst), silent = TRUE))
        if (inherits(ok, "try-error")) file.copy(f, tm)
      }
    }
    to <- file.path(tdir, "out"); dir.create(to, recursive = TRUE)
    # CLI 参数引用的数据/结果子目录必须在 out 下预先存在（官方代码按路径直读，
    # 缺失报 FileNotFoundError）
    # the data/result subfolders referenced by the CLI args must exist under out
    # up front (the official code reads them by path; missing -> FileNotFoundError)
    dir.create(file.path(to, "pUniFind_data_tmp"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(to, "pUniFind_result"), recursive = TRUE, showWarnings = FALSE)
    r <- run_once(tq, tm, to)
    if (r$ok) {
      # 把产物搬回原输出目录
      # move the outputs back into the original output folder
      if (kind == "rescore") {
        sp <- list.files(tq, pattern = "pUniFind\\.spectra$", full.names = TRUE)
        if (length(sp) > 0) file.copy(sp, qry_res_dir, overwrite = TRUE)
      } else {
        cs <- list.files(file.path(to, "pUniFind_result"), pattern = "\\.csv$",
                         full.names = TRUE)
        if (length(cs) > 0) file.copy(cs, res_dir, overwrite = TRUE)
      }
    }
  }
  # CLI 主进程退出时其 multiprocessing worker 可能仍在收尾写 .spectra
  # （实测：主进程 exit 0 后产物稍后才落盘，随即检查会误报缺失）：
  # 退出码为 0 时轮询等待产物落盘（最多约 60 s）再判定。
  # the CLI's multiprocessing workers may still be flushing the .spectra when
  # the main process exits (observed: the output lands moments after
  # the exit, so an immediate check falsely reports it missing): when the exit
  # code is 0, poll for the output (up to ~60 s) before deciding.
  if (r$ok && kind == "rescore") {
    for (w in 1:12) {
      if (length(list.files(qry_res_dir, pattern = "pUniFind\\.spectra$")) > 0) break
      Sys.sleep(5)
    }
  }
  if (!r$ok) {
    stop("pUniFind ", kind, " 运行失败（exit ", r$status, "）\npUniFind ",
         kind, " failed (exit ", r$status, "):\n",
         substr(r$log, max(1, nchar(r$log) - 2000), nchar(r$log)))
  }
  if (kind == "rescore") {
    sp <- list.files(qry_res_dir, pattern = "pUniFind\\.spectra$", full.names = TRUE)
    #：官方 CLI 只产出打分 pkl、不产出 .spectra。用随客户端分发的
    # punifind_rescore_convert.py 把谱级打分（pscore_res_*.pkl）转成 fdr0.01
    # .spectra；python 由 punifind_python 配置（需装有 torch/numpy），未配置或
    # 缺失时退回原有报错（等价于 GUI 导出待办）。
    # the official CLI writes only score pkls, not .spectra the
    # shipped punifind_rescore_convert.py turns the spectrum-level scores
    # (pscore_res_*.pkl) into an fdr0.01 .spectra; the python named by
    # punifind_python (torch/numpy required) is used, and when missing the
    # original error is kept (equivalent to the GUI-export to-do).
    if (length(sp) == 0 && !is.null(converter_py) && nzchar(converter_py) &&
        !is.null(python_exe) && nzchar(python_exe) && file.exists(converter_py)) {
      # 官方 CLI 写打分 pkl 存在延迟（实测竞态：文件在 CLI 返回后几十毫秒才
      # 落盘），轮询等待最多 120 s 再转换
      # the official CLI writes the score pkls with a delay (a measured race:
      # the file lands tens of milliseconds after the CLI returns); poll up
      # to 120 s before converting
      for (tries in seq_len(25)) {
        pks <- list.files(file.path(out_dir, "pUniFind_data_tmp"),
                          pattern = "^pscore_res_.*\\.pkl$", full.names = TRUE)
        if (length(pks) > 0) break
        if (tries == 25) break
        Sys.sleep(5)
      }
      if (length(pks) > 0) {
        pk <- pks[which.max(file.info(pks)$mtime)]
        out_sp <- file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra")
        # system2 取 stdout=TRUE 时返回值是捕获的输出行、退出码在 "status"
        # 属性里（实测竞态：转换实际成功却因 st==0 恒假误报失败，
        # spectra 在报错后 1 s 落盘）——先取属性，再判转换产物
        # with stdout=TRUE system2 returns the captured lines and the exit
        # code rides the "status" attribute (measured race: a
        # successful conversion was reported as failure because st==0 was
        # always FALSE; the spectra landed 1 s after the error) - read the
        # attribute first, then check the converted output
        conv_out <- suppressWarnings(system2(python_exe,
                          c(shQuote(converter_py), shQuote(pk), shQuote(out_sp)),
                          stdout = TRUE, stderr = TRUE))
        conv_st <- attr(conv_out, "status")
        if (is.null(conv_st)) conv_st <- 0
        if (conv_st == 0 && file.exists(out_sp)) {
          sp <- out_sp
          cat("pUniFind rescore pkl 已转换为 fdr0.01 .spectra：", out_sp, "\n")
          cat("pUniFind rescore pkl converted to the fdr0.01 .spectra: ", out_sp, "\n")
        }
      }
    }
    if (length(sp) == 0) stop("rescore 完成但未产出 pUniFind.spectra：", qry_res_dir)
    cat("pUniFind rescore 完成：", tail(sp, 1), "\n")
    return(c(rescore_spectra = tail(sp, 1)))
  }
  cs <- list.files(res_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(cs) == 0) stop("denovo 完成但未产出 csv：", res_dir)
  cat("pUniFind denovo 完成：", tail(cs, 1), "\n")
  c(denovo_csv = tail(cs, 1))
}

#' 第一步完成后写出用户步骤 2 参数文件。
#' Write the user-step-2 parameter file after step 1 completes.
#'
#' 四行 search_file 为待填占位；文件已存在时保留用户已填的 search_file（仅缺的
#' 条目加占位）；pFind 模板自动发现约定（pFind_open\param\pFind.cfg 与
#' pFind_restrict*\param\pFind.cfg）发现即写 active 键。单样本流程：文件为
#' 步骤 1 参数文件旁的 config_step2.txt。
#' The four search_file lines are placeholders for the user to fill; when the
#' file already exists the user-filled search_file entries are PRESERVED; the
#' pFind templates are auto-discovered (pFind_open\param\pFind.cfg and
#' pFind_restrict*\param\pFind.cfg) and written as active keys when found.
#' Single-sample workflow: the file is config_step2.txt next to the step-1 file.
#'
#' @param cfg1_path
#'   字符串，步骤 1 参数文件路径。
#'   Character string, the step-1 parameter-file path.
#' @param cfg1
#'   列表，read_param_file(cfg1_path) 的输出。
#'   A list, the output of read_param_file(cfg1_path).
#' @param out_path
#'   字符串，输出 config_step2.txt 的路径。
#'   Character string, the output config_step2.txt path.
#' @return
#'   返回写出的文件路径。
#'   The written file path.
#' @keywords internal
make_step2_config <- function(cfg1_path, cfg1, out_path) {
  keep <- character(0)
  keep_ty <- character(0)
  #：纯预测模式键随既有配置保留（网页零策略行写入；第一步重生成
  # 参数文件时不得冲掉——否则纯模式被静默降级为手动搜库停点）
  #: the pure-prediction flag survives regeneration (written by the
  # web shell on zero strategy rows; step 1 must not clobber it - otherwise
  # pure mode silently degrades into the manual-search stop point)
  keep_pure <- FALSE
  if (file.exists(out_path)) {
    old <- read_param_file(out_path)
    keep <- old$search_file
    keep <- keep[!is.na(keep) & keep != ""]
    kp <- if (!is.null(old$pure_prediction)) old$pure_prediction else ""
    keep_pure <- length(keep) == 0 && length(kp) > 0 && !is.na(kp[1]) &&
      identical(trimws(tail(kp, 1)), "1")
    # 保留 search_file_type此前只保留路径、类型键被
    # 丢弃，step1 重跑后类型全部按文件名规则推成 pfind_open -> "策略类型重复"报错
    keep_ty <- old$search_file_type
    if (is.null(keep_ty)) keep_ty <- rep("", length(keep))
    keep_ty <- keep_ty[!is.na(keep_ty)]
    if (length(keep_ty) < length(keep)) keep_ty <- c(keep_ty, rep("", length(keep) - length(keep_ty)))
    keep_ty <- keep_ty[seq_along(keep)]
    keep_ty[keep_ty == ""] <- "pfind_open"
  }
  # 自动发现 pFind 模板参数文件（约定目录：pFind_open/param/pFind.cfg 与 pFind_restrict*/param/pFind.cfg），
  # 发现即写入 active 键；约定目录缺失时从内置模板自动生成（fasta 取自第一步参数）——
  # 初始输入即可全自动搜库、无停点
  # auto-discover the pFind template cfgs (convention folders: pFind_open/param/pFind.cfg
  # and pFind_restrict*/param/pFind.cfg); when found they become active keys; when the
  # convention folder is missing the template is auto-generated from the built-in copy
  # (fasta taken from the step-1 parameters) - the initial inputs then run the whole
  # search chain with no stop point
  has_fasta <- !is.null(cfg1$fasta) && length(cfg1$fasta) > 0 && cfg1$fasta[1] != ""
  pfo <- file.path(cfg1$config_dir, "pFind_open", "param", "pFind.cfg")
  if (!file.exists(pfo) && has_fasta) {
    pfo <- ensure_pfind_template("open",
                                 resolve_cfg_paths(cfg1, "fasta")[1],
                                 cfg1$config_dir)
  }
  if (!file.exists(pfo)) pfo <- ""
  pfr <- ""
  for (d in list.files(cfg1$config_dir, pattern = "^pFind_restrict",
                       full.names = TRUE)) {
    f <- file.path(d, "param", "pFind.cfg")
    if (file.exists(f)) { pfr <- f; break }
  }
  if (!nzchar(pfr) && has_fasta) {
    pfr <- ensure_pfind_template("restrict",
                                 resolve_cfg_paths(cfg1, "fasta")[1],
                                 cfg1$config_dir)
  }
  lines <- c(
    paste0("# user step 2 parameter file (auto-generated after step 1)"),
    paste0("# read together with the step-1 parameter file"),
    paste0("base_config\t", rel_to(cfg1_path, out_path)),
    paste0("# auto_search: 1=软件自动调 pFind 搜库；0=用户手动搜库后填 search_file"),
    paste0("#              1=automatic pFind search; 0=manual (fill search_file)"),
    paste0("auto_search\t", if (!is.null(cfg1$auto_search) && length(cfg1$auto_search) > 0) tail(cfg1$auto_search, 1) else "1"),
    paste0("# pfind_exe: pFind 安装路径（以 bin 结尾；填版本目录如 pFind3.1 时软件自动补 bin）"),
    paste0("#            pFind install path (must end with bin; a version folder like pFind3.1 gets bin appended automatically)"),
    if (!is.null(cfg1$pfind_exe) && length(cfg1$pfind_exe) > 0 && cfg1$pfind_exe[1] != "") {
      paste0("pfind_exe\t", tail(cfg1$pfind_exe, 1))
    } else {
      c(paste0("# pfind_exe\t<pFind 安装路径>"),
        paste0("# pfind_exe\t<pFind install path>"))
    },
    paste0("# use_punifind: 1=用 pUniFind rescore/denovo；0=不用（全 pFind 流程）"),
    paste0("#               1=use pUniFind rescore/denovo; 0=no (pure-pFind workflow)"),
    paste0("# 1 时软件备好两个项目（逐 raw mgf＋.pac）并打印 GUI 待办，用户 GUI 完成后删除 merge\\peptide.txt 重跑即续接"),
    paste0("# =1 prepares the two projects (per-raw mgfs + .pac) and prints the GUI to-do; after the GUI delete merge\\peptide.txt and rerun"),
    paste0("use_punifind\t", if (!is.null(cfg1$use_punifind) && length(cfg1$use_punifind) > 0) tail(cfg1$use_punifind, 1) else "1"),
    paste0("# punifind_exe: pUniFind 安装的 punifind 目录（含 official_*.exe、ckpts、consts）；留空即手动模式"),
    paste0("#              the installed pUniFind punifind folder (official_*.exe, ckpts, consts); empty = manual mode"),
    if (!is.null(cfg1$punifind_exe) && length(cfg1$punifind_exe) > 0 && cfg1$punifind_exe[1] != "") {
      paste0("punifind_exe\t", tail(cfg1$punifind_exe, 1))
    } else {
      c(paste0("# punifind_exe\t<pUniFind 安装的 punifind 目录>"),
        paste0("# punifind_exe\t<the installed pUniFind punifind folder>"))
    },
    paste0("# punifind_rescore: 开放式搜库结果重新打分校正的运行方式（0=不运行；1=自动；2=手动）"),
    paste0("#                   how to run the open-search re-scoring (0=off; 1=auto; 2=manual)"),
    paste0("punifind_rescore\t", if (!is.null(cfg1$punifind_rescore) && length(cfg1$punifind_rescore) > 0) tail(cfg1$punifind_rescore, 1) else "0"),
    paste0("# punifind_denovo: 从头测序的运行方式（0=不运行；1=自动；2=手动）"),
    paste0("#                 how to run de novo sequencing (0=off; 1=auto; 2=manual)"),
    paste0("punifind_denovo\t", if (!is.null(cfg1$punifind_denovo) && length(cfg1$punifind_denovo) > 0) tail(cfg1$punifind_denovo, 1) else "0"),
    paste0("# use_pdeep2: 是否使用 pDeep2 理论预测谱（1=用，按 pdeep2_exe 自动运行；2=用，手填 predict_file；0=不用，预测指标按缺省处理）"),
    paste0("#             use pDeep2 predicted spectra (1=yes, auto-run via pdeep2_exe; 2=yes, fill predict_file manually; 0=no, the predict metric is treated as missing)"),
    paste0("use_pdeep2\t", if (!is.null(cfg1$use_pdeep2) && length(cfg1$use_pdeep2) > 0) tail(cfg1$use_pdeep2, 1) else "1"),
    paste0("# search_file: 每行一个；文件路径直接收，文件夹路径展开为其中全部匹配文件"),
    paste0("#              one entry per line; a file path is taken as-is, a folder path is expanded"),
    paste0("# 软件对原始结果自行执行 I→L、N→GG 替换与四种修饰校正"),
    paste0("# the software applies the I->L / N->GG replacements and the four modification corrections to the ORIGINAL results itself"),
    if (length(keep) > 0)
      unlist(Map(function(a, b) c(paste0("search_file\t", a),
                                 paste0("search_file_type\t", b)),
                 keep, keep_ty)) else c(
      paste0("# search_file\t<strategy 2: pFind open search>"),
      paste0("# search_file\t<strategy 3: pFind restricted multi-mod search>"),
      paste0("# search_file\t<strategy 4: pUniFind open search correction>"),
      paste0("# search_file\t<strategy 5: pUniFind de novo>")),
    paste0("# ---- 自动搜库（可选）：不填策略 2/3 的 search_file、改填以下三键即可由软件自动跑 pFind"),
    paste0("# ---- auto pFind search (optional): fill these instead of the strategy 2/3 search_file entries"),
    if (nzchar(pfo)) paste0("pfind_cfg_open\t", pfo) else c(
      paste0("# pfind_cfg_open\t<pFind 开放式搜库模板 pFind.cfg 绝对路径>"),
      paste0("# pfind_cfg_open\t<absolute path of the pFind open-search template pFind.cfg>")),
    if (nzchar(pfr)) paste0("pfind_cfg_restrict\t", pfr) else c(
      paste0("# pfind_cfg_restrict\t<pFind 限定式搜库模板 pFind.cfg 绝对路径>"),
      paste0("# pfind_cfg_restrict\t<absolute path of the pFind restricted-search template pFind.cfg>")),
    #：纯预测模式键回写（零策略行时保持；用户一旦填了 search_file
    # 即自动清除——回到竞争搜库流程）
    #: re-emit the pure-prediction flag (kept on zero strategy rows;
    # automatically dropped once the user fills any search_file - back to the
    # competitive-search flow)
    if (keep_pure) paste0("pure_prediction\t1") else NULL)
  writeLines(lines, out_path, useBytes = TRUE)
  cat("[用户第一步] 已生成第二步参数文件：", out_path,
      if (length(keep) > 0) paste0("（保留已填 search_file ", length(keep), " 条）") else "",
      if (keep_pure) "（纯预测模式保留）" else "",
      "\n")
  cat("[User Step 1] step-2 parameter file generated: ", out_path,
      if (length(keep) > 0) paste0(" (kept ", length(keep), " filled search_file entries)") else "",
      "\n")
  invisible(out_path)
}
#' 第二步完成后写出用户步骤 3 参数文件。
#' Write the user-step-3 parameter file after step 2 completes.
#'
#' predict_file 为待填占位；文件已存在时保留用户已填的 predict_file /
#' weights_A / threshold_A。base_config 指向步骤 1 参数文件（fasta 键在那里）。
#' predict_file is a placeholder for the user to fill; when the file already
#' exists the user-filled predict_file / weights_A / threshold_A entries are
#' PRESERVED. base_config points at the STEP-1 parameter file (the fasta key).
#'
#' @param cfg2_path
#'   字符串，步骤 2 参数文件路径。
#'   Character string, the step-2 parameter-file path.
#' @param cfg2
#'   列表，read_param_file(cfg2_path) 的输出。
#'   A list, the output of read_param_file(cfg2_path).
#' @param out_path
#'   字符串，输出 config_step3.txt 的路径。
#'   Character string, the output config_step3.txt path.
#' @return
#'   返回写出的文件路径。
#'   The written file path.
#' @keywords internal
make_step3_config <- function(cfg2_path, cfg2, out_path) {
  # base_config 指向步骤 1 参数文件（fasta 键在那里；步骤 2 文件没有 fasta 键而 extend_report_A 需要）
  # base_config points at the STEP-1 parameter file (the fasta key lives
  # there; the step-2 file has no fasta key and extend_report_A needs it)
  cfg1_path <- resolve_cfg_paths(cfg2, "base_config")[1]
  # merge_result 位于结果文件夹 result_YYYYMMDDHHMMSS/merge 内（与第二步输出一致）
  # merge_result lives inside the result folder result_YYYYMMDDHHMMSS/merge
  # (consistent with the step-2 output)
  merge_abs <- file.path(get_result_dir(cfg2$config_dir), "merge",
                         "1.combine_result.spectra")
  old_pf <- old_wa <- old_ta <- old_wb <- old_tb <- old_pg <- NULL
  if (file.exists(out_path)) {
    old <- read_param_file(out_path)
    old_pf <- old$predict_file
    old_pf <- old_pf[!is.na(old_pf) & old_pf != ""]
    old_wa <- old$weights_A
    old_wa <- old_wa[!is.na(old_wa) & old_wa != ""]
    old_ta <- old$threshold_A
    old_ta <- old_ta[!is.na(old_ta) & old_ta != ""]
    old_wb <- old$weights_B
    old_wb <- old_wb[!is.na(old_wb) & old_wb != ""]
    old_tb <- old$threshold_B
    old_tb <- old_tb[!is.na(old_tb) & old_tb != ""]
    # v1.0proteogenomics_mode 等用户键此前会被
    # 本函数重建时冲掉（只保留 5 键）——与 predict_file 同口径保留。
    old_pg <- old$proteogenomics_mode
    old_pg <- old_pg[!is.na(old_pg) & old_pg != ""]
  }
  lines <- c(
    paste0("# user step 3 parameter file (auto-generated after step 2)"),
    paste0("base_config\t", rel_to(cfg1_path, out_path)),
    paste0("merge_result\t", rel_to(merge_abs, out_path)),
    paste0("# stage A weights/threshold: four unlabeled-side metrics ul_by,ul_all,ul_cov,predict"),
    paste0("# (unified across label types, gold-calibrated; threshold_A = ACS_A gap gate)"),
    paste0("weights_A\t", if (length(old_wa) > 0) tail(old_wa, 1) else "0.1306,0.4418,0.2381,0.1895"),
    paste0("threshold_A\t", if (length(old_ta) > 0) tail(old_ta, 1) else "0.0204"),
    paste0("# stage B per-label weights (optional; six metrics ul_by,lb_by,ul_all,"),
    paste0("# lb_all,predict,pcos) and per-label tau; absent = unified W* defaults"),
    if (length(old_wb) > 0) paste0("weights_B\t", tail(old_wb, 1)) else NULL,
    if (length(old_tb) > 0) paste0("threshold_B\t", tail(old_tb, 1)) else NULL,
    paste0("# predict_file: the pDeep2 predict.txt path (filled by the user;"),
    paste0("# pDeep2 only outputs predict.txt, the software parses it directly)"),
    if (length(old_pf) > 0) paste0("predict_file\t", old_pf) else
      paste0("# predict_file\t<pDeep2 predict.txt path>"),
    paste0("# use_pdeep2: 是否使用 pDeep2 理论预测谱（1=用，按 pdeep2_exe 自动运行；2=用，手填 predict_file；0=不用，预测指标按缺省处理）"),
    paste0("#             use pDeep2 predicted spectra (1=yes, auto-run via pdeep2_exe; 2=yes, fill predict_file manually; 0=no, the predict metric is treated as missing)"),
    paste0("use_pdeep2\t", if (!is.null(cfg2$use_pdeep2) && length(cfg2$use_pdeep2) > 0) tail(cfg2$use_pdeep2, 1) else "1"),
    if (length(old_pg) > 0) paste0("# proteogenomics_mode 保留（v1.0 起不被重建冲掉；生效读取见 user_step3 口径）"),
    if (length(old_pg) > 0) paste0("proteogenomics_mode\t", tail(old_pg, 1)) else NULL)
  writeLines(lines, out_path, useBytes = TRUE)
  cat("[用户第二步] 已生成第三步参数文件：", out_path,
      if (length(old_pf) > 0) "（保留已填 predict_file）" else "",
      "\n")
  cat("[User Step 2] step-3 parameter file generated: ", out_path,
      if (length(old_pf) > 0) " (kept the filled predict_file)" else "",
      "\n")
  invisible(out_path)
}
#' 把 Peptide_state / Modification_state 字符串展开为策略组。
#' Expand a Peptide_state / Modification_state string into strategy groups.
#'
#' "Same"/"" → 仅策略 1；"13;25" → list(c(1,3), c(2,5))。
#' "Same"/"" -> strategy 1 alone; "13;25" -> list(c(1,3), c(2,5)).
#'
#' @param state_str
#'   字符串，状态分组字符串。
#'   Character string, the state-group string.
#' @return
#'   返回整数向量的列表。
#'   A list of integer vectors.
#' @keywords internal
expand_groups <- function(state_str) {
  if (is.na(state_str) || state_str == "" || state_str == "Same") return(list("a"))
  lapply(strsplit(state_str, ";")[[1]], function(z) strsplit(z, "")[[1]])
}
#' 用户步骤 3：无标记侧 isoPSV 评估（阶段 A）。
#' User step 3: unlabeled-side isoPSV evaluation (stage A).
#'
#' 按 scan 并行（PSOCK worker）。已确认口径：候选＝（序列＋修饰）组合；逐候选指标列
#' 按候选顺序以 "|" 连接；candidate_specific_ratio 首个（策略 1）候选为空、其余候选
#' 相对它（仅嵌合证据）；ACS 列留空（权重待成对阶段后固定）；判定 consistent /
#' candidate_misidentification（无 low_confidence）。输出 report_A.tsv、summary_A.txt、
#' scan_details_A.tsv（每 scan 逐候选一行，阈值重标定免重读谱图）。
#' Parallelized per scan (PSOCK workers). Confirmed conventions: candidates are
#' (sequence + modification) combinations; per-candidate metric columns are
#' "|"-joined in candidate order; candidate_specific_ratio leaves the first
#' (strategy-1) candidate empty and compares the others against it (chimeric
#' evidence only); ACS columns stay EMPTY (weights fixed after the paired
#' stage); decision is consistent / candidate_misidentification (no
#' low_confidence). Outputs report_A.tsv, summary_A.txt and scan_details_A.tsv
#' (one row per candidate per scan, so threshold re-calibration never re-reads).
#'
#' @param merged_file
#'   字符串，第二步合并表 1.combine_result.spectra。
#'   Character string, the step-2 merged table 1.combine_result.spectra.
#' @param unlabeled_mgf_files
#'   字符向量，无标记 MGF 文件。
#'   Character vector, the unlabeled MGF files.
#' @param predict_file
#'   字符串，pDeep2 predict.txt。
#'   Character string, the pDeep2 predict.txt.
#' @param weights_A
#'   数值向量，阶段 A 三权重（ul_by、ul_all、predict；设计约定 /21，
#'   统一跨标记类型；金标定权见 stageA权重_20260921）。
#'   Numeric vector, the three stage-A weights (ul_by, ul_all, predict; user
#'   ruling /21, unified across label types; gold-calibrated).
#' @param threshold_A
#'   数值，阶段 A 判定阈值（ΔACS_A：最强挑战者与 a 槽候选的加权分差门限）。
#'   Numeric, the stage-A decision threshold (the ACS_A gap between the top
#'   challenger and the a-slot candidate required to overturn the original).
#' @param out_dir
#'   字符串，输出目录（stageA\）。
#'   Character string, the output directory (stageA\).
#' @param max_cores
#'   整数，资源评估核数，默认为 1。
#'   Integer, the assessed core count; defaults to 1.
#' @param worker_cores
#'   整数或 NULL，扫描并行核数；NULL 时自动取核数−1。
#'   Integer or NULL, the scan worker count; NULL auto-takes cores-1.
#' @return
#'   返回 list(report, summary, details)——三份输出的数据。
#'   A list(report, summary, details) with the three outputs.
#' @keywords internal
stageA_evaluate <- function(merged_file, unlabeled_mgf_files,
                            predict_file,
                            weights_A, threshold_A, out_dir, max_cores = 1,
                            worker_cores = NULL,
                            label_tag = "unlabeled", label_mods = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # 双语计时标记：各阶段墙钟（elapsed）与 CPU 秒
  # bilingual timing marks: wall (elapsed) and CPU seconds per phase
  pt0 <- proc.time()
  mark <- function(tag, tag_en = tag) {
    pt1 <- proc.time()
    d <- pt1 - pt0
    cat(sprintf("[计时] %s：已用 %.1fs / CPU %.1fs\n", tag, d[3], d[1]))
    cat(sprintf("[timing] %s: elapsed %.1fs / CPU %.1fs\n", tag_en, d[3], d[1]))
    pt0 <<- pt1
  }
  aa <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"), sep = "\t", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
  aa_tab <- amino_acid_table
  dat <- read.delim(merged_file, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = character(0))
  seq_cols <- which(grepl("\\.Sequence$", colnames(dat)))
  mod_cols <- which(grepl("\\.Modification$", colnames(dat)))
  # 多引擎版：槽位字母来自合并表列前缀（a=原始结果，其余为补充策略）
  # multi-engine edition: slot letters come from the merged column prefixes
  # (a = original result, the rest are supplementary strategies)
  slot_letters_all <- sub("\\.Sequence$", "", colnames(dat)[seq_cols])
  supp_slots <- setdiff(slot_letters_all, "a")
  mark("读取合并表", "read merged table")

  # ---- pDeep2 predict.txt（pDeep2 唯一输出；直接解析）----------------------
  # ---- pDeep2 predict.txt (the ONLY pDeep2 output; parsed directly) ---------
  # use_pdeep2=0（predict_file 为空）时构造空预测表：全部候选预测余弦为 NA，
  # 下游按缺省处理（predict_missing 标注、预测指标不参与投票、不填 0）
  # when use_pdeep2=0 (empty predict_file) an empty prediction table is built:
  # every candidate's predict cosine is NA and downstream treats the metric as
  # missing (flagged by predict_missing, dropped from the votes, never 0)
  if (nzchar(predict_file) && file.exists(predict_file)) {
    pred <- read_predict_txt(predict_file)
    pred_idx <- pred$idx
    pred_int <- pred$intensities
    cat("predict entries loaded:", length(pred_int), "\n")
  } else {
    pred_idx <- setNames(integer(0), character(0))
    pred_int <- list()
    cat("use_pdeep2=0：预测指标按缺省处理\n")
    cat("use_pdeep2=0: predict metric treated as missing\n")
  }
  mark("解析 predict.txt", "parse predict.txt")

  # ---- 一个共享 PSOCK 集群同时用于 MGF 加载与扫描 ------------------------
  # ---- one shared PSOCK cluster for both the MGF loading and the scan ------
  use_cores <- if (is.null(worker_cores)) {
    max(1L, min(max_cores, parallel::detectCores() - 1L))
  } else {
    max(1L, min(worker_cores, parallel::detectCores() - 1L))
  }
  cl <- parallel::makeCluster(use_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ispsv_bootstrap_cluster(cl)
  cat("parallel workers:", use_cores, "\n")

  # ---- 加载无标记谱（共享集群；逐 worker 分块 payload）-------------------
  # ---- load unlabeled spectra (shared cluster; chunk payloads per worker) --
  ul <- load_unlabeled_mgfs(unlabeled_mgf_files, cl)
  ul_data <- ul$data
  ul_charge <- ul$charge
  ul_title_idx <- setNames(seq_along(ul_data), ul$titles)   # named vector = O(1) lookup
  cat("unlabeled spectra loaded:", length(ul_data), "\n")
  mark("加载无标记谱图", "load unlabeled spectra")

  # ---- 构建逐块 payload（块行＋该块谱图/predict 条目；以 parLapply payload 发送而非全部导出给每个 worker）----
  # ---- build per-chunk payloads (chunk rows + that chunk's spectra/predict
  # entries only; sent as parLapply payloads instead of exporting everything
  # to every worker) ----------------------------------------------------------
  chunks <- split(seq_len(nrow(dat)),
                   droplevels(cut(seq_len(nrow(dat)), 100)))  # 空组剔除：<99 行小样本不再产生空块
  mk_payload <- function(ch) {
    datsub <- dat[ch,, drop = FALSE]
    fn <- dat$File_Name[ch]
    li <- unname(ul_title_idx[fn])
    uniq <- sort(unique(li[!is.na(li)]))
    loc <- setNames(seq_along(uniq), uniq)
    # predict 子集：本块内所有非空候选条目的键
    # predict subset: keys of all non-empty candidate entries in this chunk
    keys <- character(0)
    for (k in seq_along(seq_cols)) {
      sq <- as.character(datsub[[seq_cols[k]]])
      md <- sub(";+$", "", as.character(datsub[[mod_cols[k]]]))  # 键修饰去尾分号
      nz <- !is.na(sq) & sq != ""                           # NA 守卫（空候选防 NA 传播）
      if (any(nz)) {
        # （无标记通道=0）：pDeep2 对去标记肽段预测，键修饰同步去标记
        if (!is.null(label_mods) && length(label_mods) > 0)
          md <- strip_label_mods(md, label_mods)
        z <- datsub$Charge[nz]
        keys <- c(keys, paste0(sq[nz], "|", md[nz], "|", z),
                  paste0(sq[nz], "||", z))
      }
    }
    keys <- unique(keys)
    pii <- unname(pred_idx[keys])
    ok <- !is.na(pii)
    uniq_pi <- unique(pii[ok])
    list(datsub = datsub, li = li, loc = loc,
         specs = ul_data[uniq], charges = ul_charge[uniq],
         pred_idx = setNames(match(pii[ok], uniq_pi), keys[ok]),
         pred_int = pred_int[uniq_pi])
  }
  payloads <- lapply(chunks, mk_payload)
  rm(ul_data, ul_charge, ul_title_idx, pred_idx, pred_int)
  invisible(gc())
  mark("准备逐块扫描载荷", "build scan payloads")

  scan_chunk <- function(pp) {
    mod_sites <- build_mod_sites(aa)
    dat <- pp$datsub
    #：观测前体（Exp.MH. 带 pFind 尾点；兼容无点写法）——传给
    # spectrum_params 做按质量推断的标记追加（多标记互斥口径）
    #: observed precursor (pFind writes "Exp.MH."; dotless tolerated)
    # - fed to spectrum_params for mass-inferred label appending
    mh_v <- if ("Exp.MH." %in% colnames(dat)) suppressWarnings(as.numeric(dat[["Exp.MH."]]))
            else if ("Exp.MH" %in% colnames(dat)) suppressWarnings(as.numeric(dat[["Exp.MH"]]))
            else rep(NA_real_, nrow(dat))
    if (length(mh_v) == 0) mh_v <- rep(NA_real_, nrow(dat))
    rep_rows <- vector("list", nrow(dat))
    det_rows <- vector("list", nrow(dat))
    for (j in seq_len(nrow(dat))) {
      li <- pp$li[j]
      if (is.na(li)) next
      pos <- unname(pp$loc[[as.character(li)]])
      if (length(pos) == 0) next
      spec <- pp$specs[[pos]]
      z <- pp$charges[pos]
      # （性能第二轮）：同谱多候选的 matched_peak_set 复用一次排序
      # （spec_ord 映射回原始峰下标，与函数内自排序语义一致）
      spec_mz <- as.numeric(spec[[1]])
      spec_ord <- order(spec_mz, na.last = NA)
      spec_sorted <- spec_mz[spec_ord]
      p_groups <- expand_groups(dat$Peptide_state[j])
      m_groups <- expand_groups(dat$Modification_state[j])
      # "Same" 肽段＝全部槽一致：肽段组应覆盖全部槽，否则"仅修饰冲突/仅序列冲突"
      # 会被塌缩成单 a 槽候选而失去仲裁
      # a "Same" peptide means every slot agrees: the peptide group must span ALL
      # slots, otherwise mod-only/seq-only conflicts collapse to the lone a-slot
      if (identical(as.character(dat$Peptide_state[j]), "Same"))
        p_groups <- list(unique(unlist(m_groups)))
      if (identical(as.character(dat$Modification_state[j]), "Same"))
        m_groups <- list(unique(unlist(p_groups)))
      all_slots <- c("a", supp_slots)              # 槽位字母（a 恒为原始结果）
      cand <- list()
      for (pg in p_groups) for (mg in m_groups) {
        idx <- intersect(pg, mg)                   # 字母交集＝共同策略槽
        if (length(idx) == 0) next
        ci <- seq_cols[match(idx[1], all_slots)]
        sq <- as.character(dat[j, ci])
        md <- as.character(dat[j, mod_cols[match(idx[1], all_slots)]])
        if (sq == "") next
        cand[[length(cand) + 1]] <- list(idx = idx, seq = sq, mod = md,
                                        sp = NULL, predict_cosine = NA_real_, sspec = "")  # 缺失记 NA 不填 0（sp=NULL 候选保持 NA）
      }
      if (length(cand) == 0) next
      # 等价归一后相同的（序列, 修饰）候选合并（如 0/1,Carbamidomethyl[AnyN-term] 与 1,Carbamidomethyl[H] 及 I/L 序列变体是同一个候选）
      # merge candidates that are the same (sequence, modification) after
      # equivalence normalization (e.g. 0/1,Carbamidomethyl[AnyN-term] and
      # 1,Carbamidomethyl[H] plus I/L sequence variants are ONE candidate)
      if (length(cand) > 1) {
        key <- vapply(cand, function(x) paste0(gsub("I", "L", x$seq), "@",
                        normalize_mod_str(x$mod, x$seq, mod_sites)), character(1))
        keep <- !duplicated(key)
        if (any(!keep)) {
          for (k in which(!keep)) {
            kk <- match(key[k], key)
            cand[[kk]]$idx <- sort(unique(c(cand[[kk]]$idx, cand[[k]]$idx)))
          }
          cand <- cand[keep]
        }
      }
      for (k in seq_along(cand)) {
        sq <- cand[[k]]$seq; md <- cand[[k]]$mod
        sp <- spectrum_params(sq, md, label_tag, aa, aa_tab, spec, 20, "ppm",
                              labels = label_mods, precursor_mh = mh_v[j])
        if (is.null(sp)) next
        # 预测谱缺失时不填 0，记 NA（结果文件中以空值显式标注、并在 predict_missing 列标记）
        # a missing predicted spectrum is recorded as NA, never 0 (annotated as empty in
        # the outputs and flagged by the predict_missing column)
        pr <- NA_real_
        md_ <- sub(";+$", "", md)
        if (!is.null(label_mods) && length(label_mods) > 0)
          md_ <- sub(";+$", "", strip_label_mods(md_, label_mods))
        pi <- unname(pp$pred_idx[paste0(sq, "|", md_, "|", z)])
        if (is.na(pi)) pi <- unname(pp$pred_idx[paste0(sq, "||", z)])
        if (is.na(pi)) pi <- unname(pp$pred_idx[paste0(sq, "|", md_, "|", dat$Charge[j])])  # 键电荷回退：表内 Charge
        if (!is.na(pi)) {
          pv <- pp$pred_int[[pi]]
          v <- sp$b_y_ions_intensity
          nmin <- min(length(pv), length(v))
          if (nmin > 0) {
            # （无标记通道=0）：按子离子顺序对齐，缺失位点 0 补齐
            if (!is.null(label_mods) && length(label_mods) > 0) {
              vv <- v[seq_len(nmin)]; vv[is.na(vv)] <- 0
              pvv <- pv[seq_len(nmin)]; pvv[is.na(pvv)] <- 0
              pr <- cosine(vv, pvv)
            } else {
              pr <- cosine(v[seq_len(nmin)], pv[seq_len(nmin)])
            }
          }
        }
        cand[[k]]$sp <- sp
        cand[[k]]$predict_cosine <- pr
      }

      # 相对策略 a 候选的特异离子比（嵌合证据；m/z 口径：按观测峰
      # 集合，特异峰数 ÷ 两候选命中峰并集，取 min——同一峰不论位点只计一次）
      # specific-ion ratios vs the strategy-a candidate (chimeric evidence;
      # m/z basis since: observed-peak sets, specific peaks over
      # the union of both candidates' matched peaks, taking the minimum - a
      # shared peak counts once regardless of its position in either peptide)
      spec_ratios <- setNames(rep("", length(supp_slots)),
                              paste0("spec_ratio_", supp_slots))
      c1 <- which(vapply(cand, function(x) "a" %in% x$idx, logical(1)))[1]
      if (!is.na(c1)) {
        p1 <- matched_peak_set(cand[[c1]]$seq, cand[[c1]]$mod, spec, aa, aa_tab,
                               mz_sorted = spec_sorted, mz_ord = spec_ord)
        for (k in seq_along(cand)) {
          if (k == c1) next
          pk <- matched_peak_set(cand[[k]]$seq, cand[[k]]$mod, spec, aa, aa_tab,
                                 mz_sorted = spec_sorted, mz_ord = spec_ord)
          if (length(p1) == 0 || length(pk) == 0) next
          U <- union(p1, pk)
          spec1 <- length(setdiff(p1, pk))
          speck <- length(setdiff(pk, p1))
          r <- if (length(U) > 0) as.character(round(min(spec1, speck) / length(U), 4)) else "0"
          cand[[k]]$sspec <- r
          for (sym in cand[[k]]$idx) {
            if (sym %in% supp_slots) spec_ratios[[paste0("spec_ratio_", sym)]] <- r
          }
        }
      }

      # 三指标多数投票（按谱降级）：任一候选缺某指标则本谱弃用该指标维度——缺预测时
      # 预测余弦不参与投票；缺失以空值在候选列显式标注（不填 0），并在 predict_missing 列标记
      # three-metric majority vote with per-spectrum degradation: a metric whose value
      # is missing for ANY candidate is dropped from this spectrum's vote - when the
      # prediction is missing the predict cosine does not vote; missing values are shown
      # as empty in the candidate columns (never 0) and flagged by predict_missing
      n <- length(cand)
      u_by <- sapply(cand, function(x) if (is.null(x$sp)) NA_real_ else x$sp$PSM_score_by)
      u_all <- sapply(cand, function(x) if (is.null(x$sp)) NA_real_ else x$sp$PSM_score_all)
      # 竞争匹配不引入 predict任何报告结果都不如
      # pDeep2 原生生成的谱图相似性高（理论注入谱即 pDeep2 产物、余弦理论
      # 1.00），predict 指标会系统性偏向 pDeep2 生成谱而失实——投票仅用
      # by/all 两指标；predict 仍计算并写入候选列供报告，不参与投票
      # the competition vote excludes predict no
      # reported result reaches the similarity of a pDeep2-native spectrum
      # (the theoretical injections ARE pDeep2 outputs, cosine ~1.00 by
      # construction), so the predict metric systematically favours generated
      # spectra - the vote uses by/all only; predict is still computed and
      # written into the candidate columns for reporting, just not voted
      u_pred <- sapply(cand, function(x) x$predict_cosine)
      # ---- 加权 ACS_A 判定----
      # 四指标 = ul_by/100、ul_all/100、ul_cov/100、predict（无标记侧可得，八参数
      # 版配套口径）；任一候选缺某指标则本谱弃用该维度并重归一权重（与阶段 B 的
      # 可用性降级同构）。胜者 = ACS_A 最高候选；仅当胜者非 a 槽且
      # ΔACS_A ≥ threshold_A 时改判 candidate_misidentification，否则保持原始。
      # ---- weighted ACS_A decision (user ruling: stage A carries
      # its own UNIFIED four-metric weights): metrics are ul_by/100, ul_all/100,
      # ul_cov/100 and predict; a metric missing for ANY candidate is dropped
      # for this spectrum with the remaining weights renormalized. The winner
      # is the highest-ACS_A candidate; the verdict flips only when the winner
      # is not the a-slot AND the gap clears threshold_A.
      u_cov <- sapply(cand, function(x) if (is.null(x$sp)) NA_real_ else x$sp$b_y_coverage)
      mm <- cbind(u_by / 100, u_all / 100, u_cov / 100, u_pred)
      avail <- apply(mm, 2, function(cv) all(!is.na(cv)))
      a_pos <- which(vapply(cand, function(x) "a" %in% x$idx, logical(1)))
      best_group <- if (length(a_pos) > 0) a_pos[1] else 1L
      acs <- rep(NA_real_, length(cand)); acs_gap <- NA_real_
      ww <- weights_A[avail]
      if (any(avail) && sum(ww) > 0) {           # 权重和 0 的退化配置不判（保持 a 槽）
        ww <- ww / sum(ww)
        acs <- as.numeric(mm[, avail, drop=FALSE] %*% ww)
        w_i <- which.max(acs)
        acs_gap <- acs[w_i] - acs[best_group]
        if (w_i != best_group && acs_gap >= threshold_A) best_group <- w_i
      }
      best_seq <- cand[[best_group]]$seq
      flag <- if (!("a" %in% cand[[best_group]]$idx)) "candidate_misidentification" else "consistent"

      # 逐候选列以 "|" 连接（候选＝序列＋修饰）
      # "|"-joined per-candidate columns (candidate = sequence + modification)
      cnd_show <- vapply(cand, function(x) {
        if (is.null(x$mod) || is.na(x$mod) || x$mod == "") x$seq else paste0(x$seq, "(", x$mod, ")")
      }, character(1))
      cnd_cov <- vapply(cand, function(x) if (is.null(x$sp)) "" else as.character(x$sp$b_y_coverage), character(1))
      cnd_by <- vapply(cand, function(x) if (is.null(x$sp)) "" else as.character(x$sp$PSM_score_by), character(1))
      cnd_all <- vapply(cand, function(x) if (is.null(x$sp)) "" else as.character(x$sp$PSM_score_all), character(1))
      cnd_pred <- vapply(cand, function(x) if (is.na(x$predict_cosine)) "" else as.character(round(x$predict_cosine, 4)), character(1))
      cnd_spec <- vapply(cand, function(x) x$sspec, character(1))

      rep_row <- c(
        list(
          File_Name = dat$File_Name[j],
          candidates = paste(cnd_show, collapse = "|"),
          candidate_coverage = paste(cnd_cov, collapse = "|"),
          candidate_PSM_score_by = paste(cnd_by, collapse = "|"),
          candidate_PSM_score_all = paste(cnd_all, collapse = "|"),
          candidate_predict_cosine = paste(cnd_pred, collapse = "|"),
          candidate_specific_ratio = paste(cnd_spec, collapse = "|"),
          candidate_ACS_A = paste(ifelse(is.na(acs), "",
            as.character(round(acs, 4))), collapse = "|"),
          best_group = best_group,
          best_sequence = best_seq,
          best_ACS_A = if (is.na(acs[best_group])) ""
            else as.character(round(acs[best_group], 4))),
        as.list(spec_ratios),
        list(
          predict_missing = if (any(is.na(u_pred))) 1L else 0L,
          vote_metrics = paste(c("ul_by", "ul_all", "ul_cov", "predict")[avail], collapse = ";"),
          ACS_A_gap = if (is.na(acs_gap)) "" else as.character(round(acs_gap, 4)),
          decision_A = flag))
      rep_rows[[j]] <- as.data.frame(rep_row, stringsAsFactors = FALSE)

      dts <- vector("list", n)
      for (k in seq_along(cand)) {
        spk <- cand[[k]]$sp
        dts[[k]] <- data.frame(
          File_Name = dat$File_Name[j],
          Scan_No = dat$Scan_No[j],
          Charge = dat$Charge[j],
          candidates = paste(cand[[k]]$idx, collapse = ""),
          sequence = cand[[k]]$seq,
          modification = cand[[k]]$mod,
          coverage = if (!is.null(spk)) spk$b_y_coverage else NA,
          PSM_score_by = if (!is.null(spk)) spk$PSM_score_by else NA,
          PSM_score_all = if (!is.null(spk)) spk$PSM_score_all else NA,
          predict_cosine = round(cand[[k]]$predict_cosine, 4),
          specific_ratio = cand[[k]]$sspec,
          ACS_A = if (is.na(acs[k])) "" else as.character(round(acs[k], 4)),
          is_best = k == best_group,
          decision_A = flag,
          stringsAsFactors = FALSE)
      }
      det_rows[[j]] <- do.call(rbind, dts)
    }
    list(report = if (length(rep_rows) > 0) do.call(rbind, rep_rows[!vapply(rep_rows, is.null, logical(1))]) else NULL,
         detail = if (length(det_rows) > 0) do.call(rbind, det_rows[!vapply(det_rows, is.null, logical(1))]) else NULL)
  }
  # worker 在 workers 的全局环境中评估，小共享对象（aa、aa_tab、seq_cols、mod_cols、
  # supp_slots）已导出
  # the worker evaluates in the workers' global environment, where the small
  # shared objects (aa, aa_tab, seq_cols, mod_cols, supp_slots) were exported
  environment(scan_chunk) <- globalenv()
  parallel::clusterExport(cl, c("aa", "aa_tab", "seq_cols", "mod_cols",
                                "supp_slots", "weights_A", "threshold_A",
                                "label_tag", "label_mods",
                                "scan_chunk"), envir = environment())
  mark("并行集群建立与数据分发", "cluster setup + data export")

  # 带实时进度条的批量 parLapply（已完成谱数 / 总数）
  # batched parLapply with a live progress bar (spectra done / total)
  cat("[stage A] 扫描 ", nrow(dat), " 张谱图，使用 ", use_cores, " 个进程\n")
  cat("[stage A] scanning ", nrow(dat), " spectra with ", use_cores, " workers\n")
  # 预分配结果列表（审计 P4：原 res <- c(res, ...) 每批浅拷贝整个增长列表，
  # 批次数 = 谱数/核数（如 142k 谱 ≈ 7,500 批）时累计为 O(n²) 指针拷贝）
  # pre-allocate the result list (audit P4: the old res <- c(res, ...) shallow-
  # copied the whole growing list per batch, which accumulates to O(n^2)
  # pointer copies at batch count = spectra / cores, e.g. ~7,500 for 142k)
  res <- vector("list", length(payloads))
  batch_size <- use_cores
  last_pct <- -1L
  for (b0 in seq(1, length(payloads), by = batch_size)) {
    b1 <- min(b0 + batch_size - 1, length(payloads))
    res[b0:b1] <- parallel::parLapply(cl, payloads[b0:b1], scan_chunk)
    last_scan <- max(chunks[[b1]])
    # 每推进 1% 打印一次（中英各一行）；不再用 \r 覆盖式刷屏，网页日志里不会挤成一行
    # print once per 1% gain (one Chinese + one English line); the \r-overwriting
    # progress is gone so the web log no longer squashes everything into one line
    pct <- as.integer(100 * last_scan / nrow(dat))
    if (pct != last_pct) {
      cat(sprintf("  进度：%.1f%% (%d / %d 谱图)\n",
                  100 * last_scan / nrow(dat), last_scan, nrow(dat)))
      cat(sprintf("  progress: %.1f%% (%d / %d spectra)\n",
                  100 * last_scan / nrow(dat), last_scan, nrow(dat)))
      flush.console()
      last_pct <- pct
    }
  }
  mark("逐谱扫描评估", "per-scan evaluation")

  report <- do.call(rbind, lapply(res, function(x) x$report))
  detail <- do.call(rbind, lapply(res, function(x) x$detail))
  mark("汇总块结果", "assemble chunk results")
  if (!is.null(report)) {
    write.table(report, file.path(out_dir, "report_A.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE)
    write.table(detail, file.path(out_dir, "scan_details_A.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE)
    #: spectrum-level summary is the primary count - a top5
    # multi-engine merged table holds up to K candidate rows per spectrum
    # and row-level counting inflates "errors" with rows 2..K flips (live
    # case: 42% row-level vs 2% spectrum-level). summary_A.txt counts
    # first-row verdicts per spectrum; row-level counts move to the sidecar.
    sm <- as.data.frame(table(report$decision_A[!duplicated(report$File_Name)]),
                        stringsAsFactors = FALSE)
    write.table(sm, file.path(out_dir, "summary_A.txt"), sep = "	",
                row.names = FALSE, quote = FALSE)
    sm_row <- as.data.frame(table(report$decision_A), stringsAsFactors = FALSE)
    write.table(sm_row, file.path(out_dir, "summary_A.rowlevel.txt"), sep = "	",
                row.names = FALSE, quote = FALSE)
    cat("stage A reported rows:", nrow(report),
        "| spectra:", sum(sm$Freq), "| detail rows:", nrow(detail), "
")
    print(sm)
    mark("写出报告", "write outputs")
  } else {
    cat("no reportable spectra\n")
  }
  invisible(list(report = report, details = detail))
}
#' 把 FASTA 文件读为命名字符向量。
#' Read a FASTA file into a named character vector.
#'
#' base R 实现（无需 Biostrings）；名称取 "^>" 行的第一段。
#' Base R only (no Biostrings); names come from the first field of the ">" lines.
#'
#' @param fasta_path
#'   字符串，FASTA 文件路径。
#'   Character string, the FASTA file path.
#' @return
#'   返回命名字符向量（蛋白名→序列）。
#'   A named character vector (protein name -> sequence).
#' @keywords internal
read_fasta <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  hdr <- which(startsWith(lines, ">"))
  n <- length(hdr)
  out <- character(n)
  starts <- hdr + 1L
  ends <- c(hdr[-1] - 1L, length(lines))
  names(out) <- sub("^>", "", sub(" .*$", "", lines[hdr]))
  for (i in seq_len(n)) {
    if (starts[i] <= ends[i]) out[i] <- paste(lines[starts[i]:ends[i]], collapse = "")
  }
  out
}
#' 把一条肽段匹配到蛋白序列。
#' Match one peptide against the protein sequences.
#'
#' 调用方已做 I→L 归一；返回 list(proteins, positions, contexts, ref_string)。
#' The caller has already I->L normalized the peptide; returns
#' list(proteins, positions, contexts, ref_string).
#'
#' @param peptide
#'   字符串，肽段序列。
#'   Character string, the peptide sequence.
#' @param proteins
#'   命名字符向量，蛋白组序列。
#'   A named character vector of protein sequences.
#' @return
#'   返回 list(proteins, positions, contexts, ref_string)；无命中时各字段为空字符串。
#'   A list(proteins, positions, contexts, ref_string); empty strings when no hit.
#' @keywords internal
match_peptide_protein <- function(peptide, proteins) {
  empty <- list(proteins = "", positions = "", contexts = "", ref_string = "")
  if (!nzchar(peptide)) return(empty)
  hits <- grepl(peptide, proteins, fixed = TRUE)
  if (!any(hits)) return(empty)
  pn <- sort(names(proteins)[hits])
  pos <- character(0)
  ctx <- character(0)
  for (p in pn) {
    seq <- unname(proteins[p])
    s <- as.integer(regexpr(peptide, seq, fixed = TRUE)[1])
    e <- s + nchar(peptide) - 1L
    # 前导撇号让 "start-end" 在 Excel 打开时保持文本
    # leading apostrophe keeps "start-end" as text when opened in Excel
    pos <- c(pos, paste0("'", s, "-", e))
    ctx <- c(ctx, paste0("(", substr(seq, s - 1, s - 1), ")", peptide,
                         "(", substr(seq, e + 1, e + 1), ")"))
  }
  list(proteins = paste(pn, collapse = "/"),
       positions = paste(pos, collapse = "/"),
       contexts = paste(ctx, collapse = "/"),
       ref_string = paste(paste(pn, collapse = "/"), paste(ctx, collapse = "/"),
                          sep = ";"))
}
#' 从 AA_mass_dataframe 构建修饰位点知识库。
#' Build the modification-site knowledge base from AA_mass_dataframe.
#'
#' 取带位点标记的修饰行（如 Carbamidomethyl[C]），记录哪些名称具有
#' 末端形式与残基形式，供修饰等价判断（canon_mod_entry）使用。
#' Collects the site-carrying modification rows (e.g. Carbamidomethyl[C]) and
#' records which names have terminal and residue forms, feeding the
#' modification-equivalence logic in canon_mod_entry.
#'
#' @param aa_df
#'   数据框，AA_mass_dataframe.txt 的内容。
#'   A data frame, the contents of AA_mass_dataframe.txt.
#' @return
#'   返回位点知识库（列表）。
#'   The site knowledge base (a list).
#' @keywords internal
build_mod_sites <- function(aa_df) {
  aa <- aa_df$AA
  # perl=TRUE：R 4.5.1 UCRT 的 TRE 引擎对否定类匹配失败
  has_br <- grepl("\\[[^\\[]*\\]$", aa, perl = TRUE)
  aa <- aa[has_br]
  nm <- character(length(aa))
  site <- character(length(aa))
  for (i in seq_along(aa)) {
    open <- gregexpr("\\[", aa[i])[[1]]
    last_open <- open[length(open)]
    nm[i] <- substr(aa[i], 1, last_open - 1)
    site[i] <- substr(aa[i], last_open + 1, nchar(aa[i]) - 1)
  }
  site <- gsub("[^A-Za-z-]", "", site)   # drop trailing artifacts like "])"
  sites <- tapply(site, nm, unique, simplify = FALSE)
  sites
}
#' 单条 "p,Name[Site]" 的等价比较规范形。
#' Canonical form of one "p,Name[Site]" entry for equivalence comparison.
#'
#' 末端条目对位置不敏感：写在 0 位或 1 位的 N 端修饰相同；通用 AnyN-term/AnyC-term
#' 形式在该残基能携带修饰时解析为残基字母。
#' Terminal entries are position-insensitive: N-terminal mods written at
#' position 0 or 1 are the same, and the generic AnyN-term/AnyC-term form is
#' resolved to the residue letter when that residue can carry the modification.
#'
#' @param entry
#'   字符串，一条修饰条目。
#'   Character string, one modification entry.
#' @param seq
#'   字符串，肽段序列。
#'   Character string, the peptide sequence.
#' @param sites
#'   位点知识库（build_mod_sites 输出）。
#'   The site knowledge base (build_mod_sites output).
#' @return
#'   返回规范形字符串。
#'   The canonical-form string.
#' @keywords internal
canon_mod_entry <- function(entry, seq, sites) {
  parts <- strsplit(entry, ",", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(entry)
  p <- suppressWarnings(as.integer(parts[1]))
  open <- gregexpr("\\[", parts[2])[[1]]
  if (length(open) == 0 || open[1] < 1) return(entry)
  last_open <- open[length(open)]
  nm <- substr(parts[2], 1, last_open - 1)
  site <- substr(parts[2], last_open + 1, nchar(parts[2]) - 1)
  site <- gsub("[^A-Za-z-]", "", site)
  if (site == "ProteinN-term") site <- "AnyN-term"
  letter <- ""
  if (grepl("^AnyN-term", site)) {
    p <- 1L                              # N-terminal mods: position 0 == 1
    letter <- sub("^AnyN-term", "", site)
    if (letter == "") {
      forms <- sites[[nm]]
      if (!is.null(forms)) {
        res_forms <- setdiff(forms, c("AnyN-term", "ProteinN-term", "AnyC-term"))
        aa_p <- substr(seq, 1, 1)
        if (aa_p %in% res_forms) letter <- aa_p
      }
    }
    site <- if (letter != "") letter else "AnyN-term"
  } else if (grepl("^AnyC-term", site)) {
    p <- nchar(seq)                      # C-terminal mods: last position
    letter <- sub("^AnyC-term", "", site)
    if (letter == "") {
      forms <- sites[[nm]]
      if (!is.null(forms)) {
        res_forms <- setdiff(forms, c("AnyN-term", "ProteinN-term", "AnyC-term"))
        aa_p <- substr(seq, nchar(seq), nchar(seq))
        if (aa_p %in% res_forms) letter <- aa_p
      }
    }
    site <- if (letter != "") letter else "AnyC-term"
  }
  paste0(p, ",", nm, "[", site, "]")
}
#' 整条修饰串的规范形（等价比较用）。
#' Canonical form of a full modification string (for equivalence comparison).
#'
#' 按分号拆分后逐条 canon_mod_entry 再按原序拼回。
#' Splits on semicolons, canonicalizes each entry with canon_mod_entry and
#' rejoins in the original order.
#'
#' @param mod_str
#'   字符串，修饰串。
#'   Character string, the modification string.
#' @param seq
#'   字符串，肽段序列。
#'   Character string, the peptide sequence.
#' @param sites
#'   位点知识库。
#'   The site knowledge base.
#' @return
#'   返回规范形修饰串。
#'   The canonical modification string.
#' @keywords internal
normalize_mod_str <- function(mod_str, seq, sites) {
  if (is.na(mod_str) || mod_str == "") return("")
  entries <- strsplit(mod_str, ";", fixed = TRUE)[[1]]
  entries <- entries[nzchar(entries)]
  if (length(entries) == 0) return("")
  paste(vapply(entries, canon_mod_entry, character(1), seq = seq, sites = sites),
        collapse = ";")
}
#' revise_form 分类（与 SI 对齐）。
#' revise_form classification (SI-aligned).
#'
#' 修饰比较使用等价归一后的字符串；分类：Modification location changed /
#' Different modification 等（SI 1.8 口径）。
#' Modification comparison uses the equivalence-normalized strings; classes
#' include Modification location changed / Different modification (SI 1.8).
#'
#' @param orig_seq
#'   字符串，原始序列。
#'   Character string, the original sequence.
#' @param orig_mod
#'   字符串，原始修饰串。
#'   Character string, the original modification string.
#' @param rev_seq
#'   字符串，修订序列。
#'   Character string, the revised sequence.
#' @param rev_mod
#'   字符串，修订修饰串。
#'   Character string, the revised modification string.
#' @param sites
#'   位点知识库。
#'   The site knowledge base.
#' @return
#'   返回分类字符串。
#'   The classification string.
#' @keywords internal
revise_form_classify <- function(orig_seq, orig_mod, rev_seq, rev_mod, sites) {
  om <- normalize_mod_str(orig_mod, orig_seq, sites)
  rm <- normalize_mod_str(rev_mod, rev_seq, sites)
  if (rev_seq == "" || rev_seq == orig_seq) return("modification_location_changed")
  if (rm != om) {
    if (rev_seq != orig_seq) return("different_modification")
    return("modification_location_changed")
  }
  if (rev_seq != orig_seq) return("same_weight_aa_pairs")
  "Others"
}
#' 差异字母对（3.match_count.R 的 revise_AA_pairs 移植）。
#' Differing-letter pairs (port of 3.match_count.R revise_AA_pairs).
#'
#' I/L 视作相同后比较两序列，返回差异位置与字母对信息。
#' Compares the two sequences with I and L treated as equal and returns the
#' differing positions with the letter pairs.
#'
#' @param orig_seq
#'   字符串，原始序列。
#'   Character string, the original sequence.
#' @param rev_seq
#'   字符串，修订序列。
#'   Character string, the revised sequence.
#' @return
#'   返回差异字母对字符串（或空串）。
#'   The differing-letter-pair string (or empty).
#' @keywords internal
revise_aa_pairs <- function(orig_seq, rev_seq) {
  rs <- strsplit(gsub("I", "L", rev_seq), "")[[1]]
  os <- strsplit(gsub("I", "L", orig_seq), "")[[1]]
  for (one in rs) {
    k <- match(one, os)
    if (!is.na(k)) {
      os <- os[-k]
      rs <- rs[-match(one, rs)]
    }
  }
  if (length(os) == 0 || length(rs) == 0) {
    return(paste0(paste0(sort(os), collapse = ""), "=",
                  paste0(sort(rs), collapse = "")))
  }
  if (sort(os)[1] < sort(rs)[1]) {
    return(paste0(paste0(sort(os), collapse = ""), "=",
                  paste0(sort(rs), collapse = "")))
  }
  paste0(paste0(sort(rs), collapse = ""), "=", paste0(sort(os), collapse = ""))
}
#' Error_reason - SI 1.8(3)-(5) 错误类型分类。
#' Error_reason - SI 1.8(3)-(5) error-type classification.
#'
#' 逐级判定，以 ";" 堆叠，"same_weight_aa_pairs_or_permutation" 除外（其他原因
#' 伴随的结果，仅单独报告）：① SAP（严格：File_Name 在 SAP 索引中且 revise_peptide
#' == 索引 revised_peptide）——先堆叠；② N 端非特异酶切（ref 匹配、左翼非 K/R）；
#' ③ 修饰校正（修饰类型与位置差异）；④ 等重氨基酸对或互换（序列不同、修饰相同）——
#' 仅在无其他原因时。
#' Judged level by level, stackable with ";" EXCEPT "Same weight AA pairs or
#' permutation" (an attendant result, reported only alone): 1. SAP (strict);
#' 2. N-terminal non specific cleavage (ref match, left flank not K/R);
#' 3. Modification correction (type and location differences); 4. Same weight
#' AA pairs or permutation (sequences differ, modification identical) - only
#' when no other reason was assigned.
#'
#' @param rev_form
#'   字符串，revise_form 分类。
#'   Character string, the revise_form class.
#' @param rev_pep
#'   字符串，修订肽段。
#'   Character string, the revised peptide.
#' @param ref_protein
#'   字符串，ref_Protein 字段。
#'   Character string, the ref_Protein field.
#' @param sap_hit
#'   字符串，SAP 命中标记（空串为无）。
#'   Character string, the SAP hit marker (empty = none).
#' @return
#'   返回 Error_reason 字符串。
#'   The Error_reason string.
#' @keywords internal
error_reason_compute <- function(rev_form, rev_pep, ref_protein, sap_hit) {
  parts <- character(0)
  if (nzchar(sap_hit)) parts <- c(parts, "SAP")
  if (nzchar(ref_protein)) {
    ctx_part <- strsplit(ref_protein, ";")[[1]]
    if (length(ctx_part) >= 2) {
      flanks <- vapply(strsplit(ctx_part[2], "/")[[1]], function(c) substr(c, 2, 2),
                       character(1))
      # v1.0N 端非特异判据按酶切类型参数化（早期实现写死 K/R）；
      # C 端酶看前一残基，N 端酶看肽段首残基，无酶切恒特异
      first_aa <- substr(rev_pep, 1, 1)
      if (!any(vapply(flanks, function(fl) enzyme_nterm_ok(iso_enzyme(), fl, first_aa),
                      logical(1))))
        parts <- c(parts, "n_terminal_non_specific_cleavage")
    }
  }
  if (rev_form %in% c("modification_location_changed", "different_modification")) {
    parts <- c(parts, "modification_correction")
  }
  if (length(parts) == 0 && rev_form == "same_weight_aa_pairs") {
    parts <- "same_weight_aa_pairs_or_permutation"
  }
  paste(parts, collapse = ";")
}
#' 把阶段 A 报告扩充参考格式列（向量化）。
#' Extend the stage A report with the reference-format columns (vectorized).
#'
#' 追加 X1-X5.Sequence/Modification、Peptide_state/Modification_state、revise_*、
#' ref_Protein、ref_Protein_position（'XX-XX 转义）、Error_reason、SAP 等列，
#' 使 report_A 达到参考格式 38 列；纯 pFind 流程（缺策略列）列填空、结构不变。
#' Appends X1-X5.Sequence/Modification, Peptide_state/Modification_state, the
#' revise_* columns, ref_Protein, ref_Protein_position (with the leading
#' apostrophe), Error_reason, SAP and so on, reaching the reference 38 columns;
#' for the pure-pFind flow the missing strategy columns are filled empty with
#' the structure unchanged.
#'
#' @param report
#'   数据框，阶段 A 的 report_A。
#'   A data frame, the stage-A report_A.
#' @param details
#'   数据框，scan_details_A。
#'   A data frame, the scan_details_A.
#' @param merged
#'   数据框，第二步合并表。
#'   A data frame, the step-2 merged table.
#' @param fasta
#'   命名字符向量，参考蛋白组。
#'   A named character vector, the reference proteome.
#' @param sap
#'   数据框或 NULL，SAP 索引。
#'   A data frame or NULL, the SAP index.
#' @return
#'   返回扩充后的 report 数据框。
#'   The extended report data frame.
#' @keywords internal
extend_report_A <- function(report, details, merged, fasta, sap = NULL) {
  # 多引擎版：合并表列名已是字母开头（a.Sequence 等），合法 R 列名，无需 X 前缀；
  # 按合并表实际槽列动态透传
  # multi-engine edition: merged columns are already letter-leading (a.Sequence
  # etc.) and valid R names - no X prefix; the actual slot columns are copied
  # through dynamically
  slot_cols <- colnames(merged)[grepl("^[a-z]\\.(Sequence|Modification)$",
                                      colnames(merged))]
  aa_df <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"), sep = "\t", header = TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE)
  sites <- build_mod_sites(aa_df)
  n <- nrow(report)
  new_cols <- c(slot_cols, "Peptide_state",
                "Modification_state", "revise_state", "revise_peptide",
                "revise_modification", "revise_form", "revise_AA_pairs",
                "Replaced_peptide", "ref_Protein", "ref_Protein_position",
                "Error_reason", "SAP")
  for (cn in new_cols) report[[cn]] <- ""

  idx_m <- match(report$File_Name, merged$File_Name)
  if (any(!is.na(idx_m))) {
    # 只启用的策略槽才有列：缺失槽不占列
    # only the enabled strategy slots carry columns: absent slots take no column
    for (sc in slot_cols) {
      report[[sc]] <- if (sc %in% colnames(merged)) {
        as.character(merged[idx_m, sc])
      } else rep("", nrow(report))
    }
    report$Peptide_state <- as.character(merged[idx_m, "Peptide_state"])
    report$Modification_state <- as.character(merged[idx_m, "Modification_state"])
  }

  wb <- which(as.logical(details$is_best))
  det_best <- details[wb,, drop = FALSE]
  idx_d <- match(report$File_Name, det_best$File_Name)
  best_syms <- rep("a", n)
  best_mod <- rep("", n)
  has_d <- !is.na(idx_d)
  if (any(has_d)) {
    best_syms[has_d] <- as.character(det_best$candidates[idx_d[has_d]])
    best_mod[has_d] <- as.character(det_best$modification[idx_d[has_d]])
  }
  x1s <- report[["a.Sequence"]]
  x1m <- report[["a.Modification"]]
  best_seq <- as.character(report$best_sequence)

  is_orig <- grepl("a", best_syms)
  revise_state <- ifelse(is_orig, "Original", "Changed")
  revise_pep <- ifelse(is_orig, "", best_seq)
  revise_mod <- ifelse(is_orig, "", best_mod)
  revert <- !is_orig & (grepl("loss", revise_mod) | grepl("1,Acetyl", revise_mod))
  revise_pep[revert] <- x1s[revert]
  revise_mod[revert] <- x1m[revert]
  mod_eq <- rep(FALSE, n)
  wh0 <- which(!is_orig)
  if (length(wh0) > 0) {
    # v1.0>2000 行时并行化 normalize_mod_str 逐行比较。
    # Macaca 890K 行的 extend_report_A 主要耗时段之一。
    # v1.0 parallelize normalize_mod_str comparison for >2000 rows.
    if (length(wh0) > 2000 && requireNamespace("parallel", quietly = TRUE)) {
      nw_eq <- min(4L, parallel::detectCores() - 1L)
      cat("[extend_report_A] mod_eq comparison: ", length(wh0), " rows, ", nw_eq, " workers\n", sep = "")
      cl_eq <- parallel::makeCluster(nw_eq)
      on.exit(parallel::stopCluster(cl_eq), add = TRUE)
      #：clusterExport 必须传 envir = environment()（局部变量在
      # .GlobalEnv 找不到）；normalize_mod_str 依赖顶层 canon_mod_entry，
      # worker 侧无闭包链须一并导出
      parallel::clusterExport(cl_eq, c("normalize_mod_str", "canon_mod_entry",
                                       "revise_mod", "revise_pep",
                                       "x1m", "x1s", "sites"),
                              envir = environment())
      eq_chunks <- split(wh0, ceiling(seq_along(wh0) / (ceiling(length(wh0) / (nw_eq * 4)))))
      eq_res <- parallel::parLapplyLB(cl_eq, eq_chunks, function(idx) {
        vapply(idx, function(i)
          normalize_mod_str(revise_mod[i], revise_pep[i], sites) ==
          normalize_mod_str(x1m[i], x1s[i], sites), logical(1))
      })
      mod_eq[wh0] <- unlist(eq_res, use.names = FALSE)
    } else {
      mod_eq[wh0] <- vapply(wh0, function(i)
        normalize_mod_str(revise_mod[i], revise_pep[i], sites) ==
        normalize_mod_str(x1m[i], x1s[i], sites), logical(1))
    }
  }
  back_orig <- !is_orig & (gsub("I", "L", revise_pep) == gsub("I", "L", x1s) & mod_eq)
  revise_state[back_orig] <- "Original"
  revise_pep[back_orig] <- ""
  revise_mod[back_orig] <- ""
  changed <- revise_state == "Changed"

  rev_form <- rep("", n)
  rev_pairs <- rep("", n)
  replaced <- rep("", n)
  if (any(changed)) {
    wh <- which(changed)
    # v1.0>2000 行时并行化 revise_form_classify + revise_aa_pairs。
    # v1.0 parallelize form classification and AA-pair analysis for >2000 rows.
    if (length(wh) > 2000 && requireNamespace("parallel", quietly = TRUE)) {
      nw_cls <- min(4L, parallel::detectCores() - 1L)
      cat("[extend_report_A] form classify + AA pairs: ", length(wh), " rows, ", nw_cls, " workers\n", sep = "")
      cl_cls <- parallel::makeCluster(nw_cls)
      on.exit(parallel::stopCluster(cl_cls), add = TRUE)
      #：同上——revise_form_classify 间接经 normalize_mod_str 需要
      # canon_mod_entry，两者均须导出给 worker
      parallel::clusterExport(cl_cls, c("revise_form_classify", "revise_aa_pairs",
                                         "normalize_mod_str", "canon_mod_entry",
                                         "x1s", "x1m", "revise_pep", "revise_mod", "sites"),
                              envir = environment())
      cls_chunks <- split(wh, ceiling(seq_along(wh) / (ceiling(length(wh) / (nw_cls * 4)))))
      cls_res <- parallel::parLapplyLB(cl_cls, cls_chunks, function(idx) {
        fm <- vapply(idx, function(i)
          revise_form_classify(x1s[i], x1m[i], revise_pep[i], revise_mod[i], sites),
          character(1))
        pr <- vapply(idx, function(i)
          revise_aa_pairs(x1s[i], revise_pep[i]), character(1))
        list(fm = fm, pr = pr)
      })
      rev_form[wh] <- unlist(lapply(cls_res, `[[`, "fm"), use.names = FALSE)
      rev_pairs[wh] <- unlist(lapply(cls_res, `[[`, "pr"), use.names = FALSE)
    } else {
      rev_form[wh] <- vapply(wh, function(i)
        revise_form_classify(x1s[i], x1m[i], revise_pep[i], revise_mod[i], sites),
        character(1))
      rev_pairs[wh] <- vapply(wh, function(i)
        revise_aa_pairs(x1s[i], revise_pep[i]), character(1))
    }
    replaced[wh] <- gsub("I", "L", revise_pep[wh])
  }

  norm_peps <- unique(gsub("I", "L", revise_pep[nzchar(revise_pep)]))
  ref_str <- setNames(rep("", length(norm_peps)), norm_peps)
  ref_pos <- setNames(rep("", length(norm_peps)), norm_peps)
  if (length(norm_peps) > 0) {
    # v1.0肽→蛋白匹配并行化。
    # Macaca 报告此处为单核 2.6h（50K+ unique peptides × 72K proteins）。
    # 每条肽段的蛋白匹配完全独立（embarrassingly parallel），>200 条时
    # 用 parLapplyLB 分发到多 worker。低于 200 条时走原顺序循环（序列化
    # 开销不划算）。match_peptide_protein 保持原逻辑不变。
    # v1.0 parallelize peptide→protein matching.
    # Macaca reports this as single-core 2.6h (50K+ peptides × 72K proteins).
    # Each peptide's protein scan is fully independent; >200 peptides triggers
    # parLapplyLB with multiple workers. Below 200 runs the original sequential
    # loop (serialization overhead not worthwhile). match_peptide_protein unchanged.
    if (length(norm_peps) > 200 && requireNamespace("parallel", quietly = TRUE)) {
      nw <- min(5L, parallel::detectCores() - 1L)
      cat("[extend_report_A] peptide-protein matching: ", length(norm_peps),
          " peptides, ", nw, " workers (parallel)\n", sep = "")
      cat("[extend_report_A] peptide-protein matching: ", length(norm_peps),
          " peptides, ", nw, " workers (parallel)\n", sep = "")
      cl <- parallel::makeCluster(nw)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      #：fasta 为函数局部变量，须显式传 environment()
      parallel::clusterExport(cl, c("fasta", "match_peptide_protein"),
                              envir = environment())
      hits_list <- parallel::parLapplyLB(cl, norm_peps, function(p) {
        match_peptide_protein(p, fasta)
      })
      for (pi in seq_along(norm_peps)) {
        ref_str[[pi]] <- hits_list[[pi]]$ref_string
        ref_pos[[pi]] <- hits_list[[pi]]$positions
      }
    } else {
      for (p in norm_peps) {
        m <- match_peptide_protein(p, fasta)
        ref_str[[p]] <- m$ref_string
        ref_pos[[p]] <- m$positions
      }
    }
  }
  key <- gsub("I", "L", revise_pep)
  ref_protein <- ifelse(nzchar(revise_pep), unname(ref_str[key]), "")
  ref_positions <- ifelse(nzchar(revise_pep), unname(ref_pos[key]), "")

  sap_hit <- rep("", n)
  if (!is.null(sap) && nrow(sap) > 0) {
    idx_s <- match(report$File_Name, sap$File_Name)
    hit <- !is.na(idx_s) & revise_pep == as.character(sap$revised_peptide[idx_s])
    hit[is.na(hit)] <- FALSE
    if (any(hit)) {
      sap_hit[hit] <- paste(as.character(sap$Proteins[idx_s[hit]]),
                            as.character(sap$Modification[idx_s[hit]]), sep = ";")
    }
  }

  err <- vapply(seq_len(n), function(i)
    error_reason_compute(rev_form[i], revise_pep[i], ref_protein[i], sap_hit[i]),
    character(1))

  report$revise_state <- revise_state
  report$revise_peptide <- revise_pep
  report$revise_modification <- revise_mod
  report$revise_form <- rev_form
  report$revise_AA_pairs <- rev_pairs
  report$Replaced_peptide <- replaced
  report$ref_Protein <- ref_protein
  report$ref_Protein_position <- ref_positions
  report$Error_reason <- err
  report$SAP <- ifelse(sap_hit != "" & !nzchar(ref_protein), sap_hit, "")
  report[] <- lapply(report, function(x) {
    v <- as.character(x)
    v[is.na(v)] <- ""
    v[v == "NA"] <- ""
    v
  })
  report
}
#' 把 Peptide_state / Modification_state 策略组合并为候选组字符串。
#' Combine the Peptide_state / Modification_state strategy groups into one candidate-group string.
#'
#' 按（肽组, 修饰组）组合，每组一个代表（生产 4.heavy_spectra_match.R 移植）。
#' One representative per (peptide group, modification group) combination
#' (production 4.heavy_spectra_match.R port).
#'
#' @param peptide_state
#'   字符串，Peptide_state。
#'   Character string, the Peptide_state.
#' @param modification_state
#'   字符串，Modification_state。
#'   Character string, the Modification_state.
#' @return
#'   返回候选组字符串（如 "1;23"）。
#'   The candidate-group string (e.g. "1;23").
#' @keywords internal
combine_peptide_modification_state <- function(peptide_state, modification_state) {
  # 多引擎版：状态串符号为单字母 a-z；"Same" 是占位词，其字母（a/m/e）不能
  # 计入策略号——只拆分非 Same 的状态串
  # multi-engine edition: state-string symbols are single letters a-z; "Same"
  # is a placeholder word whose letters (a/m/e) must NOT count as strategy
  # symbols - split only the non-Same state strings
  ps_chars <- if (is.na(peptide_state) || peptide_state %in% c("", "Same"))
    character(0) else strsplit(peptide_state, "", fixed = TRUE)[[1]]
  ms_chars <- if (is.na(modification_state) || modification_state %in% c("", "Same"))
    character(0) else strsplit(modification_state, "", fixed = TRUE)[[1]]
  all_nums <- unique(c(ps_chars, ms_chars))
  all_nums <- all_nums[all_nums %in% letters]
  all_nums <- all_nums[all_nums != ""]
  if (length(all_nums) == 0) return("a")
  if (peptide_state == "Same") {
    num_peptide_group <- rep(1, length(all_nums))
  } else {
    peptide_groups <- strsplit(peptide_state, ";", fixed = TRUE)[[1]]
    num_peptide_group <- sapply(all_nums, function(num) {
      idx <- which(sapply(peptide_groups, function(g) num %in% strsplit(g, "", fixed = TRUE)[[1]]))
      if (length(idx) == 0) return(NA)
      idx
    })
  }
  if (modification_state == "Same") {
    num_mod_group <- rep(1, length(all_nums))
  } else {
    mod_groups <- strsplit(modification_state, ";", fixed = TRUE)[[1]]
    num_mod_group <- sapply(all_nums, function(num) {
      idx <- which(sapply(mod_groups, function(g) num %in% strsplit(g, "", fixed = TRUE)[[1]]))
      if (length(idx) == 0) return(0)
      idx
    })
  }
  combo_id <- paste(num_peptide_group, num_mod_group, sep = "_")
  unique_combos <- unique(combo_id)
  result_groups <- sapply(unique_combos, function(combo) {
    nums_in_group <- all_nums[combo_id == combo]
    paste(sort(nums_in_group), collapse = "")
  })
  paste(result_groups, collapse = ";")
}
#' SILAC 基础模型：在 K/R 位追加 K6/R10 标记修饰。
#' SILAC basic model: append the K6/R10 labels as modifications at the K/R positions.
#'
#' 生产 zeroP 模型（Pro 转换暂未建模）。
#' The production zeroP model (Pro conversion not modeled yet).
#'
#' @param seq
#'   字符串，肽段序列。
#'   Character string, the peptide sequence.
#' @param mod
#'   字符串，原始修饰串。
#'   Character string, the original modification string.
#' @return
#'   返回追加 K6/R10 后的修饰串。
#'   The modification string with the K6/R10 labels appended.
#' @keywords internal
silac_labeled_mod <- function(seq, mod) {
  letters <- strsplit(seq, "", fixed = TRUE)[[1]]
  pos <- which(letters %in% c("K", "R"))
  if (length(pos) == 0) return(mod)
  extra <- vapply(pos, function(p) {
    if (letters[p] == "K") paste0(p, ",Label_13C(6)[K];") else paste0(p, ",Label_13C(6)15N(4)[R];")
  }, character(1))
  paste0(mod, paste(extra, collapse = ""))
}
#' 为一张无标记谱的单个候选定位并打分标记谱对（阶段 B）。
#' Locate and score the labeled spectrum pair for ONE candidate of ONE unlabeled spectrum (stage B).
#'
#' raw_params/raw_objs 属于同一 raw 的原始 MGF（含无标记与标记谱）。候选不可用时返回 NULL
#' （未知修饰或无标记谱子离子不足），否则返回含无标记谱匹配指标、选定标记谱标题与其指标、
#' 成对余弦（无标记 vs 标记 by_ions_long 匹配强度向量）的列表。注：正文级 ">=4 子离子"
#' 规则属于下游筛选阶段；生产 R 配对脚本没有逐候选子离子过滤，这里同样不设。
#' raw_params/raw_objs belong to the SAME raw's original MGF (unlabeled + labeled).
#' Returns NULL when the candidate is unusable, otherwise a list with the light
#' match metrics, the chosen labeled spectrum title and its metrics, and the
#' paired cosine (unlabeled vs labeled by_ions_long matched-intensity vectors).
#' NOTE: the manuscript-level ">= 4 fragments" rule applies at the downstream
#' screening stage; the production R pairing script has NO per-candidate
#' fragment filter, so none is applied here either.
#'
#' @param seq
#'   字符串，候选序列。
#'   Character string, the candidate sequence.
#' @param mod
#'   字符串，候选修饰串。
#'   Character string, the candidate modification string.
#' @param charge
#'   整数，母离子电荷。
#'   Integer, the precursor charge.
#' @param light_spec
#'   数据框，无标记谱峰表。
#'   A data frame, the unlabeled spectrum peaks.
#' @param light_rt
#'   数值，无标记谱保留时间（秒）。
#'   Numeric, the unlabeled retention time (seconds).
#' @param raw_params
#'   数据框，同一 raw 原始 MGF 的参数表。
#'   A data frame, the parameter table of the same raw's original MGF.
#' @param raw_objs
#'   MGF 对象列表，同一 raw 的原始谱图（无标记＋标记）。
#'   A list of MGF objects, the same raw's original spectra (unlabeled + labeled).
#' @param label_type
#'   字符串，标记类型默认为 "N15"。
#'   Character string, the label type (N15 / C13 / SILAC); defaults to "N15".
#' @param aa
#'   数据框或 NULL，AA_mass_dataframe 内容（NULL 时自动读）。
#'   A data frame or NULL, the AA_mass_dataframe contents (auto-read when NULL).
#' @param aa_tab
#'   氨基酸对照表或 NULL。
#'   The amino-acid table or NULL.
#' @param rt_tol
#'   数值，RT 容差（秒），默认为 20。
#'   Numeric, the RT tolerance (seconds); defaults to 20.
#' @param mass_tol
#'   数值，标记母离子中性质量窗口（Da），默认为 10。
#'   Numeric, the labeled precursor neutral-mass window (Da); defaults to 10.
#' @param ppm
#'   数值，标记母离子 m/z 窗口（ppm），默认为 20。
#'   Numeric, the labeled precursor m/z window (ppm); defaults to 20.
#' @param min_frags
#'   整数，无标记谱最少匹配子离子（保留参数，当前不设过滤），默认为 4。
#'   Integer, the unlabeled-spectrum minimum fragments (kept as a parameter; no
#'   filter is applied); defaults to 4.
#' @return
#'   返回列表（paired/light/label/paired_cosine 等）或 NULL。
#'   A list (paired/light/label/paired_cosine etc.) or NULL.
#' @keywords internal
# RT 窗口纯函数findInterval 区间 + 下标升序排序，与旧
# which(rt >= lo & rt <= hi) 逐位等价（集合同、枚举顺序同；NA RT 丢弃；
# rt_tol/light_rt 为 NA 返回空）。rt_ord 由调用方每 raw 预排序传入。
rt_window_indices <- function(rt, light_rt, rt_tol, rt_ord = NULL) {
  if (is.na(light_rt) || is.na(rt_tol) || length(rt) == 0) return(integer(0))
  if (is.null(rt_ord)) rt_ord <- order(rt, na.last = NA)
  rt_sorted <- rt[rt_ord]
  i_lo <- findInterval(light_rt - rt_tol, rt_sorted, left.open = TRUE) + 1L
  i_hi <- findInterval(light_rt + rt_tol, rt_sorted)
  if (i_lo <= i_hi) sort(rt_ord[i_lo:i_hi]) else integer(0)
}
# ---- 多标记行式口径------------------------
# label_type 为多行键：每行一个标记。取值两类——第一类全原子稳定同位素替换：
# "N15"/"C13"；第二类修饰形式标记（格式与 modification.ini 一致，名称=原子
# 类型与数量）：如 "Dimethyl[AnyN-term]"、"Label_13C(6)[K]"、
# "Label_13C(6)15N(4)[R]"。空行/无键＝无标记。兼容旧单值口径
#。
parse_label_spec <- function(label_type) {
  lt <- trimws(as.character(label_type))
  lt <- lt[!is.na(lt) & nzchar(lt)]
  #：自定义 "|" 行注册到 option——理论离子构建器经
  # augment_aa_mass_table() 把 名称→计算质量 补进 AA_mass_dataframe，
  # 修饰串按名查质量的全部既有逻辑对自定义名无感生效。
  options(isopsv_custom_labels = if (any(grepl("|", lt, fixed = TRUE))) {
    rows <- lapply(lt[grepl("|", lt, fixed = TRUE)],
                   function(x) tryCatch(parse_label_row(x), error = function(e) NULL))
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) do.call(rbind, rows) else NULL
  } else NULL)
  if (length(lt) == 0) return(list(legacy = "unlabeled", class1 = NULL, mods = character(0)))
  if (length(lt) == 1 && toupper(lt[1]) %in% c("UNLABELED", "N15", "C13", "SILAC")) {
    c1 <- if (toupper(lt[1]) %in% c("N15", "C13")) toupper(lt[1]) else NULL
    return(list(legacy = toupper(lt[1]), class1 = c1, mods = character(0)))
  }
  up <- toupper(lt)
  c1 <- if (any(up %in% c("N15", "C13"))) up[up %in% c("N15", "C13")][1] else NULL
  mods <- lt[!up %in% c("N15", "C13", "UNLABELED")]
  list(legacy = NULL, class1 = c1, mods = unique(mods))
}
# ---- v1.0第二类修饰形式标记——用户自定义修饰行 ----
# 行语法（label_type 的值，单行四字段、竖线分隔）：<名称>|<位置>|<氨基酸集合>|<元素组成>
#   位置：ANY 任意位置 / PEP_N 肽段N端 / PEP_C 肽段C端 / PRO_N 蛋白N端 / PRO_C 蛋白C端
#         （中文短语与 pFind 风格 NORMAL/AnyN-term/AnyC-term/ProteinN-term/ProteinC-term 同收）
#   氨基酸集合：* = 全部氨基酸；否则单字母集合（如 STY、KR）
#   元素组成：modification.ini 风格 H(2)C(2)O(1)，计数可负、支持同位素元素
#         （2H/D、13C、15N、18O、34S）；标记质量 = Σ 元素单价质量 × 计数
# 例：Phospho|ANY|STY|H(1)O(3)P(1)；NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)
# 旧预设名（无竖线，如 Label_13C(6)[K]）仍按 AA_mass_dataframe 查表，完全兼容。
# (user ruling): class-2 modification-form labels accept user-defined
# rows name|position|aa-set|composition; mass is computed from element counts.
# Legacy preset names (no "|") keep the AA_mass_dataframe lookup.
ELEMENT_MONO_MASSES <- c(
  H = 1.00782503223, `2H` = 2.01410177812, D = 2.01410177812,
  C = 12.0, `13C` = 13.00335483507, N = 14.00307400443,
  `15N` = 15.00010889888, O = 15.99491461957, `18O` = 17.99915961286,
  S = 31.9720711744, `34S` = 33.967867004, P = 30.97376163,
  Se = 79.9165218, Na = 22.989769282, K = 38.9637064864,
  Ca = 39.96259086, Fe = 55.93493633, Zn = 63.92914201,
  Cl = 34.968852682, Br = 78.9183371, I = 126.9044719)

# 元素组成串 → 命名计数向量（H(2)C(2)O(1)、2H(4)、C(-6)13C(6)、裸式 H2O 均可；
# 非法语法/未知元素/空串一律 stop（双语））。计数为 0 的元素剔除；同名元素累加。
parse_comp_str <- function(s) {
  s <- gsub("[[:space:]]", "", as.character(s))
  if (!nzchar(s))
    stop("元素组成为空 / empty composition string", call. = FALSE)
  patt <- "^([0-9]{0,2}[A-Z][a-z]?)(?:\\((-?[0-9]+)\\)|([0-9]*))(.*)$"
  els <- character(0); cnts <- integer(0); guard <- 0L
  while (nzchar(s)) {
    guard <- guard + 1L
    if (guard > 300L)
      stop("元素组成解析超长（疑似非法串）/ composition too long: ", s, call. = FALSE)
    m <- regmatches(s, regexec(patt, s))[[1]]
    if (length(m) != 5L)
      stop("元素组成语法非法 / invalid composition syntax: ", s, call. = FALSE)
    el <- m[2]
    cnt <- if (nzchar(m[3])) as.integer(m[3])
           else if (nzchar(m[4])) as.integer(m[4]) else 1L
    if (!(el %in% names(ELEMENT_MONO_MASSES)))
      stop("未知元素 / unknown element: ", el, "（支持 / supported: ",
           paste(names(ELEMENT_MONO_MASSES), collapse = " "), "）", call. = FALSE)
    i <- match(el, els)
    if (is.na(i)) { els <- c(els, el); cnts <- c(cnts, cnt) }
    else cnts[i] <- cnts[i] + cnt
    s <- m[5]
  }
  keep <- cnts != 0L
  setNames(cnts[keep], els[keep])
}

comp_mass <- function(s) {
  cnt <- parse_comp_str(s)
  if (length(cnt) == 0L) return(0)
  sum(ELEMENT_MONO_MASSES[names(cnt)] * cnt)
}

label_pos_normalize <- function(p) {
  p1 <- trimws(as.character(p)); up <- toupper(p1)
  if (up %in% c("ANY", "ANYWHERE", "NORMAL") || p1 %in% c("任意位置", "任意", "全部位置")) return("ANY")
  if (up %in% c("PEP_N", "PEPN", "ANYN-TERM", "ANYNTERM", "NTERM", "N-TERM") ||
      p1 %in% c("肽段N端", "肽段N末端")) return("PEP_N")
  if (up %in% c("PEP_C", "PEPC", "ANYC-TERM", "ANYCTERM", "CTERM", "C-TERM") ||
      p1 %in% c("肽段C端", "肽段C末端")) return("PEP_C")
  if (up %in% c("PRO_N", "PRON", "PROTEINN-TERM", "PROTEIN_N") ||
      p1 %in% c("蛋白N端", "蛋白质N端", "蛋白N末端")) return("PRO_N")
  if (up %in% c("PRO_C", "PROC", "PROTEINC-TERM", "PROTEIN_C") ||
      p1 %in% c("蛋白C端", "蛋白质C端", "蛋白C末端")) return("PRO_C")
  stop("未知修饰位置 / unknown modification position: ", p1,
       "（ANY/PEP_N/PEP_C/PRO_N/PRO_C 或 任意位置/肽段N端/肽段C端/蛋白N端/蛋白C端）",
       call. = FALSE)
}

is_custom_label_row <- function(x) grepl("|", as.character(x), fixed = TRUE)

# 自定义修饰行 → (name, residue, mass, pos)。名称不得含 [ ];, 与制表符
# （修饰串格式 "pos,name;" 的安全前提）；氨基酸集合去重保序、限标准字母。
parse_label_row <- function(x) {
  f <- trimws(strsplit(as.character(x), "|", fixed = TRUE)[[1]])
  if (length(f) != 4L)
    stop("自定义修饰行应为 4 个竖线分隔字段（名称|位置|氨基酸|组成）/ a custom label row needs 4 pipe-separated fields (name|position|aa|composition): ",
         x, call. = FALSE)
  name <- f[1]
  if (!nzchar(name) || grepl("[][;,\t]", name))
    stop("修饰名称非法（不得为空或含 [ ];, 制表符）/ invalid modification name: ",
         f[1], call. = FALSE)
  pos <- label_pos_normalize(f[2])
  aas <- toupper(f[3])
  if (!nzchar(aas) || aas %in% c("*", "ALL", "ANY", "全部", "全部氨基酸")) {
    aas <- "*"
  } else {
    letters_v <- strsplit(aas, "")[[1]]
    allowed <- strsplit("ACDEFGHIKLMNPQRSTVWYUO", "")[[1]]
    if (anyNA(match(letters_v, allowed)))
      stop("氨基酸集合非法（仅标准单字母）/ invalid amino-acid set: ", f[3],
           call. = FALSE)
    aas <- paste(unique(letters_v), collapse = "")
  }
  data.frame(name = name, residue = aas, mass = comp_mass(f[4]), pos = pos,
             stringsAsFactors = FALSE)
}

#：自定义行进修饰串用的带标签名（pFind 命名惯例）——理论离子构建器
# 靠名字里的 "N-term"/"C-term" 路由端基质量、按名查表取质量：
#   ANY+集合 → Name[STY]；ANY+* → Name（裸名）；PEP_N/PRO_N → Name[AnyN-term]<AA?>
#   PEP_C/PRO_C → Name[AnyC-term]<AA?>（蛋白端理论侧按肽段端近似，与引擎侧一致）
label_mod_tag <- function(r) {
  if (r$pos %in% c("PEP_N", "PRO_N"))
    return(paste0("[AnyN-term", if (r$residue == "*") "" else r$residue, "]"))
  if (r$pos %in% c("PEP_C", "PRO_C"))
    return(paste0("[AnyC-term", if (r$residue == "*") "" else r$residue, "]"))
  if (r$residue == "*") return("")
  paste0("[", r$residue, "]")
}
# 该自定义行可能出现在表/修饰串里的全部名字（裸名 + 带标签名）；
# r 为 parse_label_row 的输出行（不再重解析组成式）。
label_table_names <- function(r) unique(c(r$name, paste0(r$name, label_mod_tag(r))))
# strip/匹配用别名集（原始行值 + 裸名 + 带标签名）
label_mod_aliases <- function(lb) {
  if (!is_custom_label_row(lb)) return(lb)
  unique(c(lb, label_table_names(parse_label_row(lb))))
}
# 把 option 里的自定义行（parse_label_spec 注册）补进 AA_mass_dataframe：
# 名字 = 裸名与带标签名，Mono/N15/C13 三个质量列同值（组成式算出的 Δ）。
# 表已含同名行时不覆盖（用户表优先）。
augment_aa_mass_table <- function(aa_df) {
  cl <- getOption("isopsv_custom_labels")
  if (is.null(cl) || nrow(cl) == 0L) return(aa_df)
  key <- paste(cl$name, cl$pos, cl$residue, cl$mass, collapse = ";")
  if (identical(attr(aa_df, "isopsv_augmented"), key)) return(aa_df)
  add_names <- unique(unlist(lapply(seq_len(nrow(cl)),
                                    function(k) label_table_names(cl[k,, drop = FALSE]))))
  new <- add_names[!add_names %in% aa_df$AA]
  if (length(new) == 0L) {
    attr(aa_df, "isopsv_augmented") <- key
    return(aa_df)
  }
  bare <- sub("\\[[^]]*\\]$", "", new)
  extra <- do.call(rbind, lapply(seq_along(new), function(i) {
    k <- match(bare[i], cl$name); if (is.na(k)) k <- 1L
    row <- as.list(rep(NA, ncol(aa_df))); names(row) <- names(aa_df)
    row$AA <- new[i]
    row$Mono_weight <- cl$mass[k]; row$N15_weight <- cl$mass[k]; row$C13_weight <- cl$mass[k]
    dfi <- as.data.frame(row, stringsAsFactors = FALSE)
    names(dfi) <- names(aa_df)  # 空列名/重复列名会被 as.data.frame 改写——强制还原
    dfi
  }))
  out <- rbind(aa_df, extra)
  attr(out, "isopsv_augmented") <- key
  out
}

# 标记修饰的质量与残基解析（引擎自动搜库注入 NeuCode 标记用；
#：自定义 "|" 行走组成式计算，旧预设名仍查 AA_mass_dataframe）。
# 旧名 [AA] 括号给残基；[AnyN-termQ] 等端部括号折算为 pos 码 + 纯残基。
label_mod_masses <- function(labels) {
  if (is.null(labels) || length(labels) == 0) return(NULL)
  aa <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"), sep = "\t",
                   header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  out <- lapply(labels, function(lab) {
    if (is_custom_label_row(lab)) return(parse_label_row(lab))
    res <- if (grepl("\\[[^]]+\\]$", lab)) sub("^.*\\[([^]]+)\\]$", "\\1", lab) else ""
    i <- match(lab, aa$AA)
    mass <- if (!is.na(i)) as.numeric(aa$Mono_weight[i]) else NA_real_
    pos <- "ANY"; r2 <- res
    if (grepl("^AnyN-term", res))       { pos <- "PEP_N"; r2 <- sub("^AnyN-term", "", res) }
    else if (grepl("^AnyC-term", res))  { pos <- "PEP_C"; r2 <- sub("^AnyC-term", "", res) }
    else if (grepl("^ProteinN-term", res)) { pos <- "PRO_N"; r2 <- sub("^ProteinN-term", "", res) }
    else if (grepl("^ProteinC-term", res)) { pos <- "PRO_C"; r2 <- sub("^ProteinC-term", "", res) }
    data.frame(name = lab, residue = r2, mass = mass, pos = pos,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# v1.0标记修饰行生成。comet/msfragger 共用编号续接逻辑
# （编号接在模板已有 variable_mod(NN) 之后避免撞车）；sage 走 JSON 条目。
#：残基字段先经 engine_mod_token 按引擎语法折算端部（comet 裸 n/c =
# 任意残基端；msfragger [^ / ]$ = 任意残基端、n/c 前缀 = 特定残基端）；
# 蛋白端（PRO_N/PRO_C）comet/msfragger 无原生支持，按肽段端近似并在注入时提示。
engine_mod_token <- function(residue, pos, engine = c("comet", "msfragger")) {
  engine <- match.arg(engine)
  all20 <- "ACDEFGHIKLMNPQRSTVWY"
  aa <- if (is.na(residue) || !nzchar(residue) || residue == "*") all20 else residue
  nterm <- pos %in% c("PEP_N", "PRO_N"); cterm <- pos %in% c("PEP_C", "PRO_C")
  if (engine == "comet") {
    if (nterm) return("n")
    if (cterm) return("c")
    return(aa)
  }
  if (nterm) return(if (identical(aa, all20)) "[^" else paste0("n", aa))
  if (cterm) return(if (identical(aa, all20)) "]$" else paste0("c", aa))
  aa
}
label_mod_lines <- function(pl, label_mass, pattern, fmt, engine = "comet") {
  i_vm <- grep(pattern, pl)
  n0 <- if (length(i_vm) > 0)
    max(as.integer(sub("^variable_mod_?([0-9]+).*", "\\1", pl[i_vm]))) else 0L
  pos_v <- if (!is.null(label_mass$pos)) label_mass$pos else rep("ANY", nrow(label_mass))
  approx <- pos_v %in% c("PRO_N", "PRO_C")
  if (any(approx))
    cat("[标记注入] 蛋白端标记按肽段端近似注入 ", engine,
        "（该引擎无蛋白端可变修饰原生语法）：",
        paste(label_mass$name[approx], collapse = "、"), "\n",
        "[label inject] protein-terminal labels approximated as peptide-terminal for ",
        engine, ": ", paste(label_mass$name[approx], collapse = ", "), "\n", sep = "")
  for (k in seq_len(nrow(label_mass))) {
    n0 <- n0 + 1L
    pl <- c(pl, sprintf(fmt, n0, label_mass$mass[k],
                         engine_mod_token(label_mass$residue[k], pos_v[k], engine)))
  }
  pl
}
sage_label_mod_entries <- function(label_mass) {
  # sage 模板 variable_mods 为残基键对象，端部标记无对应语法——跳过并提示
  #（其余引擎正常携带；合并口径不受影响）。
  pos_v <- if (!is.null(label_mass$pos)) label_mass$pos else rep("ANY", nrow(label_mass))
  term <- pos_v %in% c("PEP_N", "PEP_C", "PRO_N", "PRO_C")
  if (any(term))
    cat("[sage 标记注入] sage variable_mods 仅支持残基位点，端部标记不注入 sage：",
        paste(label_mass$name[term], collapse = "、"), "\n",
        "[sage label inject] sage variable_mods is residue-keyed only; terminal labels not injected into sage: ",
        paste(label_mass$name[term], collapse = ", "), "\n", sep = "")
  lm2 <- label_mass[!term,, drop = FALSE]
  grp <- split(round(lm2$mass, 6), lm2$residue)
  paste0(sprintf('      "%s": [%s]', names(grp),
                 vapply(grp, function(m) paste(sprintf("%.6f", m), collapse = ", "),
                        character(1))),
         collapse = ",\n")
}
label_tau_of <- function(spec) {
  if (!is.null(spec$legacy)) {
    if (spec$legacy == "SILAC") return(0.0183)
    return(0.0083)
  }
  if (is.null(spec$class1) && length(spec$mods) > 0) return(0.0183)
  0.0083
}
# 把修饰形式的标记按适用残基/N 端追加进修饰串（内部格式 pos,name;...）
#：自定义 "|" 行按 位置码+氨基酸集合 定位——ANY 集合内所有残基位、
# PEP_N/PRO_N → 0、PEP_C/PRO_C → 肽长（蛋白端在此无蛋白上下文，按肽段端
# 近似，理论侧口径与引擎侧一致）；修饰串里写带标签名（Name[STY] /
# Name[AnyN-term]，与 pFind 命名惯例一致——理论端基质量路由靠名字里的
# N-term/C-term，质量查表经 augment_aa_mass_table 补行）。
append_label_mods <- function(seq, mod, labels) {
  for (lab in labels) {
    if (is_custom_label_row(lab)) {
      r <- parse_label_row(lab)
      nm <- paste0(r$name, label_mod_tag(r))
      locs <- if (r$pos %in% c("PEP_N", "PRO_N")) {
        0L
      } else if (r$pos %in% c("PEP_C", "PRO_C")) {
        nchar(seq)
      } else if (r$residue == "*") {
        seq_len(nchar(seq))
      } else {
        v <- as.integer(gregexpr(paste0("[", r$residue, "]"), seq)[[1]])
        v[v > 0]
      }
      for (pos in locs) {
        seg <- paste0(pos, ",", nm, ";")
        mod <- if (nzchar(mod)) paste0(mod, seg) else seg
      }
      next
    }
    m <- regmatches(lab, regexec("\\[([^]]+)\\]$", lab))[[1]]
    aa <- if (length(m) == 2) m[2] else ""
    locs <- if (aa == "AnyN-term") {
      0L
    } else if (nzchar(aa) && nchar(aa) == 1) {
      v <- as.integer(unlist(gregexpr(aa, seq, fixed = TRUE)))
      v[v > 0]
    } else integer(0)
    for (pos in locs) {
      seg <- paste0(pos, ",", lab, ";")
      mod <- if (nzchar(mod)) paste0(mod, seg) else seg
    }
  }
  mod
}

# 从修饰串中去除指定标记修饰（v1.0 无标记通道口径：NeuCode 等全标记样本
# 做 pDeep2 理论预测前，把肽段修饰列表里的标记修饰去掉——如
# "0,Acetyl[AnyN-term];2,K602[K];3,Oxidation[M]" 去掉 K602[K] 后保留
# "0,Acetyl[AnyN-term];3,Oxidation[M]"）。
strip_label_mods <- function(mod, labels) {
  if (is.null(labels) || length(labels) == 0 || is.null(mod)) return(mod)
  if (length(mod) > 1L)
    return(vapply(mod, strip_label_mods, character(1), labels = labels,
                  USE.NAMES = FALSE))
  if (is.na(mod) || !nzchar(mod)) return(mod)
  toks <- strsplit(mod, ";", fixed = TRUE)[[1]]
  toks <- toks[nzchar(toks)]
  drop <- vapply(toks, function(tk) {
    nm <- sub("^[^,]*,\\s*", "", tk)
    #：自定义行按别名集匹配（原始行值/裸名/带标签名 Name[STY] 等）
    any(vapply(labels, function(lb) nm %in% label_mod_aliases(lb), logical(1)))
  }, logical(1))
  out <- toks[!drop]
  if (length(out) == 0) "" else paste0(paste(out, collapse = ";"), ";")
}

#：按观测前体质量推断每肽的标记组合（全标记样品多标记互斥口径）。
# 背景：多行 label_type此前逐行
# 循环把每个适用位点上叠加全部同族标记（K 位 K602+K080 双计质量）。真实
# 样品同一位点只带一种标记，且同肽各 K 位通常同通道同标记。
# 口径：
#   ① 修饰串里已有的同族标记保留（denovo/pFind-Filtered 结果自带 K602[K]
#      等修饰名时按原位点精确计质，不重复追加——先 strip 防双计）；
#   ② 无标记的适用位点按观测前体（Exp.MH，M+H）与无标记理论 pep_mass 的
#      质量差，在"各族选无/选族内一行 × 全部适用位点"的组合网格里取
#      |Δ| 最小且过 ppm 容差（默认 20 ppm，K602 vs K080 差 36 mDa 在
#      1000 Da 下 20 ppm 可分）的组合追加；
#   ③ 推断失败（无前体、理论质量不可得或差值不过容差）：保留原修饰串
#      （不追加），忠实于已有标记名。
#: infer the per-peptide label combination from the observed
# precursor mass (mutually-exclusive multi-label, fully-labeled samples).
# The old append loop stacked EVERY same-family label row onto every
# applicable site (K got K602+K080 double-counted). Real samples carry one
# label per site, and one peptide's sites are usually all one channel.
# Rules: (1) label names already present in the mod string are kept verbatim
# (strip-then-extend prevents double counting); (2) unlabelled applicable
# sites get the grid combination ("none / one family row x all sites" per
# residue-family) closest to observed-minus-theoretical within a ppm
# tolerance; (3) on inference failure the mod string is returned unchanged.
append_label_mods_by_mass <- function(seq, mod, labels, obs_mh,
                                      aa = NULL, aa_tab = NULL,
                                      tol_ppm = 20) {
  if (is.null(labels) || length(labels) == 0) return(mod)
  if (is.null(mod) || is.na(mod)) mod <- ""
  if (is.na(obs_mh) || obs_mh <= 0)
    return(append_label_mods(seq, mod, labels))
  lm <- label_mod_masses(labels)
  if (is.null(lm) || nrow(lm) == 0 || any(is.na(lm$mass)))
    return(append_label_mods(seq, mod, labels))
  #：token 生成与 append_label_mods 同口径——custom 行
  # （名称|位置|残基|组成）给带标签名 Name[K]/Name[AnyN-term]，legacy
  # 预设名保留原文；label_mod_masses 的 name 列是裸名（"K602"），直接
  # 写进修饰串会与质量表 Name[K] 行失配，理论离子按 0 质量错算。
  #: tokens mirror append_label_mods - custom rows get the
  # tagged name Name[K]/Name[AnyN-term], legacy preset names stay verbatim;
  # label_mod_masses' bare name column ("K602") would miss the mass-table
  # Name[K] row and the theoretical ions would carry a wrong mass.
  lm <- do.call(rbind, lapply(labels, function(lab) {
    if (is_custom_label_row(lab)) {
      r <- parse_label_row(lab)
      data.frame(token = paste0(r$name, label_mod_tag(r)),
                 residue = r$residue, pos = r$pos, mass = r$mass,
                 stringsAsFactors = FALSE)
    } else {
      r <- label_mod_masses(lab)
      data.frame(token = lab, residue = r$residue, pos = r$pos, mass = r$mass,
                 stringsAsFactors = FALSE)
    }
  }))
  mod0 <- strip_label_mods(mod, labels)
  if (is.null(aa)) aa <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"),
                                    sep = "\t", header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)
  if (is.null(aa_tab)) aa_tab <- amino_acid_table
  bp <- cached_by_ions_position(seq, mod0, "Unlabel", aa, aa_tab)
  if (length(bp) == 1 || is.null(bp$pep_mass) || is.na(bp$pep_mass))
    return(mod)
  diff <- obs_mh - bp$pep_mass
  tol <- tol_ppm * 1e-6 * obs_mh
  if (abs(diff) > 50) return(mod)   # 前体差超过 50 Da：非标记问题，不追加
  # 族 = 残基集 + 位置码（同族行互斥）；逐族给"无 / 某行×全部位点"选项
  fam_key <- paste(lm$residue, lm$pos, sep = "@")
  fams <- split(lm, fam_key, drop = TRUE)
  choices <- lapply(fams, function(f) {
    first <- f[1,, drop = FALSE]
    locs <- if (first$pos %in% c("PEP_N", "PRO_N")) 0L
            else if (first$pos %in% c("PEP_C", "PRO_C")) nchar(seq)
            else if (first$residue == "*") seq_len(nchar(seq))
            else {
              v <- as.integer(gregexpr(paste0("[", first$residue, "]"), seq)[[1]])
              v[v > 0]
            }
    if (length(locs) == 0) return(list(list(rows = integer(0), mass = 0)))
    c(list(list(rows = integer(0), mass = 0)),
      lapply(seq_len(nrow(f)), function(k)
        list(rows = k, mass = f$mass[k] * length(locs))))
  })
  # 笛卡尔组合（族数 × 行数通常 ≤ 3×3 = 9）取 |Δ| 最小者
  grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE)
  best_j <- NA_integer_; best_d <- Inf
  for (g in seq_len(nrow(grid))) {
    m <- 0
    for (fi in seq_along(choices)) m <- m + choices[[fi]][[grid[g, fi]]]$mass
    d <- abs(diff - m)
    if (d < best_d) { best_d <- d; best_j <- g }
  }
  if (is.na(best_j) || best_d > tol) return(mod)
  # 追加各族选中行（带标签名，与 append_label_mods 同名口径）
  out <- mod0
  for (fi in seq_along(choices)) {
    k <- choices[[fi]][[grid[best_j, fi]]]$rows
    if (length(k) == 0) next
    f <- fams[[fi]][k,, drop = FALSE]
    nm <- f$token
    locs <- if (f$pos %in% c("PEP_N", "PRO_N")) 0L
            else if (f$pos %in% c("PEP_C", "PRO_C")) nchar(seq)
            else if (f$residue == "*") seq_len(nchar(seq))
            else {
              v <- as.integer(gregexpr(paste0("[", f$residue, "]"), seq)[[1]])
              v[v > 0]
            }
    for (pos in locs) out <- paste0(out, pos, ",", nm, ";")
  }
  out
}

paired_match_candidate <- function(seq, mod, charge, light_spec, light_rt,
                                   raw_params, raw_objs, label_type = "N15",
                                   aa = NULL, aa_tab = NULL,
                                   rt_tol = 20, mass_tol = 10, ppm = 20,
                                   min_frags = 4, rt_ord = NULL,
                                   labels = NULL) {
  if (is.null(aa)) aa <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"), sep = "\t",
                                    header = TRUE, check.names = FALSE,
                                    stringsAsFactors = FALSE)
  if (is.null(aa_tab)) aa_tab <- amino_acid_table
  lab <- label_type
  # 理论离子：无标记 vs 标记质量模型。SILAC 并行跑生产 zeroP与（含 Pro 肽段的）fullP模型；每张候选标记谱按各自 m/z 窗口对应的模型打分。
  # theoretical ions: unlabeled vs labeled mass model. SILAC runs the
  # production zeroP (K6/R10 only) and, for Pro-containing peptides, the
  # fullP (all prolines labeled, "P_C5N1") model in parallel; each candidate
  # labeled spectrum is scored under the model of its own m/z window.
  # （性能第二轮）：理论离子走进程级哈希缓存（跨候选/跨分块复用）
  ti_ul <- cached_by_ions_position(seq, mod, "Unlabel", aa, aa_tab)
  ti_lb_p <- NULL
  # v1.0多标记行式口径——labels 为修饰形式的标记
  # 清单（名称=原子类型与数量，如 "Dimethyl[AnyN-term]"、"Label_13C(6)[K]"），
  # 全部标记按适用残基/N 端追加进修饰串后生成标记态理论离子；lab 为类 1
  # 全标记标签或 "Unlabel"（仅修饰标记时）。
  if (!is.null(labels) && length(labels) > 0) {
    mod_lb <- append_label_mods(seq, mod, labels)
    ti_lb <- cached_by_ions_position(seq, mod_lb, lab, aa, aa_tab)
  } else if (lab == "SILAC") {
    mod_lb <- silac_labeled_mod(seq, mod)
    ti_lb <- cached_by_ions_position(seq, mod_lb, "Unlabel", aa, aa_tab)
    if (grepl("P", seq, fixed = TRUE)) {
      p_pos <- which(strsplit(seq, "", fixed = TRUE)[[1]] == "P")
      mod_fullp <- paste0(mod_lb, paste0(";", p_pos, ",P_C5N1;", collapse = ""))
      ti_lb_p <- cached_by_ions_position(seq, mod_fullp, "Unlabel", aa, aa_tab)
    }
  } else {
    ti_lb <- cached_by_ions_position(seq, mod, lab, aa, aa_tab)
  }
  if (length(ti_lb) == 1) return(NULL)
  im_ul <- ions_match(ti_ul, light_spec, mass_shift_tolerence = ppm, ppm_or_Da = "ppm")
  n_ul_frag <- sum(im_ul[["by_ions_long.match_intensity"]] > 0)
  # 注：正文级 ">=4 子离子" 规则属于下游筛选阶段；生产 R 配对脚本没有逐候选子离子过滤，这里同样不设
  # 标记母离子 m/z 窗口：标记 m/z 的 ±20 ppm（生产规则）；C13 下界做 -2 Da 质量亏损补偿
  # NOTE: the manuscript-level ">= 4 fragments" rule applies at the downstream
  # screening stage; the production R pairing script has NO per-candidate
  # fragment filter, so none is applied here either
  # labeled precursor m/z window(s): +/-20 ppm around the labeled m/z
  # (production rule); C13 lower bound gets the -2 Da mass-defect compensation
  H <- 1.007825
  # P0-3pep_mass 为质子化单同位素质量（M+H），正确前体
  # m/z = (M + z·H)/z = (pep_mass + (z-1)·H)/z。旧公式 (pep_mass + H)/z 只在
  # z=2 时正确——z=1 偏高 1.0078 Da、z=3 偏高 0.336 Da，阶段 B ±20 ppm 配对
  # 窗口下 1+/3+ 谱永远配不上。
  mz_of <- function(ti) precursor_mz(ti$pep_mass, charge)
  # RT 窗口排序由调用方每
  # raw 算一次经 rt_ord 传入（函数内挂 attr 因 R 写时复制从不跨调用生效）；
  # 窗口取 findInterval 区间 O(log R + 区间长)；left.open=TRUE 保证下界含
  # 等值峰（与 which(>=) 同口径）；NA RT 行按 which 口径丢弃；**窗口下标
  # sort 恢复"原始下标升序"**——旧 which() 的枚举顺序，下游决胜对候选顺序
  # 敏感（并列取首个，近平手案例真实存在）。纯函数 rt_window_indices 承载
  # 全部窗口逻辑（对拍测试直接覆盖）
  rt_idx <- rt_window_indices(raw_params$RTinSeconds, light_rt, rt_tol, rt_ord)
  win <- function(ti, is_p) {
    mz <- mz_of(ti)
    hi <- mz * (1 + ppm * 1e-6)
    lo <- mz * (1 - ppm * 1e-6) - if (lab == "C13") 2 / charge else 0
    ci <- rt_idx[raw_params$PepMass[rt_idx] >= lo &
                   raw_params$PepMass[rt_idx] <= hi]
    data.frame(ci = ci, pstate = rep(if (is_p) 1L else 0L, length(ci)))
  }
  win0 <- win(ti_lb, FALSE)
  if (!is.null(ti_lb_p) && length(ti_lb_p) > 1) {
    cand <- rbind(win0, win(ti_lb_p, TRUE))
  } else {
    cand <- win0
  }
  light_metrics <- list(coverage = im_ul[["by_ions_coverage"]],
                        PSM_score_by = im_ul[["PSM_score_by"]],
                        PSM_score_all = im_ul[["PSM_score_all"]],
                        n_frag = n_ul_frag,
                        match_int = im_ul[["by_ions_long.match_intensity"]])
  if (length(cand$ci) == 0) {
    return(list(paired = FALSE, light = light_metrics, label_title = "",
                label = NULL, paired_cosine = 0, n_candidates = 0L))
  }
  # 给每张候选标记谱打分；生产取 PSM_score_all 最高（并列取首个）——每候选的成对余弦仍报告，最佳配对余弦入 best_paired_label_cosine
  # score every candidate labeled spectrum; production keeps the highest
  # PSM_score_all (first on ties) - the paired cosine of each candidate is
  # still reported, the best pair cosine goes to best_paired_label_cosine
  best <- NULL
  for (k in seq_len(nrow(cand))) {
    ci <- cand$ci[k]
    ti_use <- if (cand$pstate[k] == 1L) ti_lb_p else ti_lb
    im_lb <- ions_match(ti_use, raw_objs[[ci]]@spectrum,
                        mass_shift_tolerence = ppm, ppm_or_Da = "ppm",
                        N15_label_state = (lab == "N15"),
                        C13_label_state = (lab == "C13"))
    dc <- cosine(im_ul[["by_ions_long.match_intensity"]],
                 im_lb[["by_ions_long.match_intensity"]])
    if (is.nan(dc)) dc <- 0
    # 近平手决胜标记侧 PSM 分差 <1% 视为不可区分，改按成对余弦定胜负，
    # 避免稀疏谱以毫厘分差压过成对证据更好的候选（案例 pH6.14907 61.10 vs 14932 61.07，
    # 余弦 0.65 vs 0.90）；分差 >=1% 时维持 PSM_score_all 最高规则，与历史结果保持一致。
    # near-tie tiebreak when labeled-side PSM_score_all values are
    # within 1% (indistinguishable), the higher paired cosine decides; otherwise
    # the higher PSM_score_all wins as before, keeping prior runs reproducible.
    if (is.null(best)) {
      best <- list(ci = ci, dc = dc, im = im_lb, pstate = cand$pstate[k])
    } else {
      s_new <- im_lb[["PSM_score_all"]]; s_old <- best$im[["PSM_score_all"]]
      hi <- max(abs(s_new), abs(s_old))
      near_tie <- is.finite(hi) && hi > 0 && abs(s_new - s_old) / hi < 0.01
      if ((s_new > s_old && (!near_tie || dc >= best$dc)) ||
          (near_tie && dc > best$dc)) {
        best <- list(ci = ci, dc = dc, im = im_lb, pstate = cand$pstate[k])
      }
    }
  }
  # 生产怪癖：共享 N15/C13 脚本仅在候选谱多于 1 张时才用硬编码 N15_label_state=TRUE 重算最终标记匹配（重算位于 length>1 分支内）；恰一张候选时保留循环自身的标记校正匹配。为对账复现；SILAC 脚本改用标记校正状态重算。
  # production quirk: the shared N15/C13 script recomputes the FINAL label
  # match with N15_label_state=TRUE hardcoded ONLY when more than one
  # candidate spectrum exists (the recompute sits inside the length>1
  # branch); with exactly one candidate the loop's own labeling-correct
  # match is kept. Reproduced for reconciliation; the SILAC script
  # recomputes with the labeling-correct states instead.
  if (lab %in% c("N15", "C13") && nrow(cand) > 1) {
    im_f <- ions_match(if (best$pstate == 1L) ti_lb_p else ti_lb,
                       raw_objs[[best$ci]]@spectrum,
                       mass_shift_tolerence = ppm, ppm_or_Da = "ppm",
                       N15_label_state = TRUE, C13_label_state = FALSE)
    dc_f <- cosine(im_ul[["by_ions_long.match_intensity"]],
                   im_f[["by_ions_long.match_intensity"]])
    if (is.nan(dc_f)) dc_f <- 0
    best$im <- im_f
    best$dc <- dc_f
  }
  list(paired = TRUE, light = light_metrics,
       label_title = raw_params$Title[best$ci],
       label = list(coverage = best$im[["by_ions_coverage"]],
                    PSM_score_by = best$im[["PSM_score_by"]],
                    PSM_score_all = best$im[["PSM_score_all"]],
                    match_int = best$im[["by_ions_long.match_intensity"]]),
       paired_cosine = round(best$dc, 4),
       pro_state = if (best$pstate == 1L) "fullP" else "zeroP",
       n_candidates = nrow(cand))
}


# ---- 终端丢失规范化----------------
# Terminal-loss normalization (added; output format aligned with the
# panel_f master table):
#   1) 位置 1 的末端类修饰（[AnyN-term] / [ProteinN-term）归一为位置 0；
#   2) Met-loss[ProteinN-termM]（位置 0/1 且序列以 M 开头）＝起始甲硫氨酸丢失：
#      序列去首位 M；+Acetyl 版新增 0,Acetyl[AnyN-term]，纯 Met-loss 无新增；
#      其余位置 >=2 的修饰前移 1 位；旧位置 1 的残基类修饰（锚定在被去除的 M）丢弃；
#   3) Arg-loss[AnyC-termR]（位置＝序列长且序列以 R 结尾）＝C 端精氨酸丢失：
#      序列去尾 R、无修饰残留，其余修饰位置不动（C 端去除不影响 N 端锚定位置）。
normalize_terminal_loss <- function(seqv, modv) {
  term1 <- function(tok) grepl("[AnyN-term]", tok, fixed = TRUE) ||
                          grepl("[ProteinN-term", tok, fixed = TRUE)
  toks <- if (is.na(modv) || !nzchar(modv)) character(0) else {
    t <- trimws(strsplit(modv, ";", fixed = TRUE)[[1]]); t[nzchar(t)]
  }
  if (length(toks) == 0) return(list(s = seqv, m = ""))
  pos_of <- function(tok) suppressWarnings(as.integer(sub(",.*$", "", tok)))
  typ_of <- function(tok) sub("^[0-9]+,", "", tok)
  toks <- vapply(toks, function(tt) {
    pp <- pos_of(tt)
    if (!is.na(pp) && pp == 1L && term1(tt)) sprintf("0,%s", typ_of(tt)) else tt
  }, character(1), USE.NAMES = FALSE)
  if (!any(grepl("Met-loss", toks, fixed = TRUE)) &&
      !any(grepl("Arg-loss[AnyC-termR]", toks, fixed = TRUE))) {
    mm <- paste(toks, collapse = ";")
    if (nzchar(mm)) mm <- paste0(mm, ";")                        # 源格式尾分号
    return(list(s = seqv, m = mm))
  }
  # 两遍法：先判 M 丢失是否触发（决定其他修饰是否前移；R 丢失不移位）
  mfire <- grepl("^M", seqv) && any(vapply(toks, function(tt) {
    pp <- pos_of(tt)
    (grepl("Met-loss+Acetyl[ProteinN-termM]", tt, fixed = TRUE) ||
     grepl("Met-loss[ProteinN-termM]", tt, fixed = TRUE)) && pp %in% c(0L, 1L)
  }, logical(1)))
  out_t <- character(0); add_t <- character(0)
  conv_r <- FALSE
  for (tt in toks) {
    pp <- pos_of(tt)
    if (grepl("Met-loss+Acetyl[ProteinN-termM]", tt, fixed = TRUE) &&
        pp %in% c(0L, 1L) && grepl("^M", seqv)) {
      add_t <- c(add_t, "0,Acetyl[AnyN-term]")
    } else if (grepl("Met-loss[ProteinN-termM]", tt, fixed = TRUE) &&
               pp %in% c(0L, 1L) && grepl("^M", seqv)) {
      next
    } else if (grepl("Arg-loss[AnyC-termR]", tt, fixed = TRUE) &&
               pp == nchar(seqv) && grepl("R$", seqv)) {
      conv_r <- TRUE
    } else if (mfire && !is.na(pp) && pp == 1L && !term1(tt)) {
      next
    } else if (mfire && !is.na(pp) && pp >= 2L) {
      out_t <- c(out_t, sprintf("%d,%s", pp - 1L, typ_of(tt)))
    } else {
      out_t <- c(out_t, tt)
    }
  }
  if (mfire) seqv <- substr(seqv, 2, nchar(seqv))
  if (conv_r) seqv <- substr(seqv, 1, nchar(seqv) - 1)
  mm <- paste(c(out_t, add_t), collapse = ";")
  if (nzchar(mm)) mm <- paste0(mm, ";")                           # 保留源格式尾分号
  list(s = seqv, m = mm)
}


# ---- SAP 修饰代入序列（增；修 pDeep2 对 X->Y 型修饰无法预测的问题）----
# Apply SAP-type substitutions (X->Y notation) into the sequence so pDeep2 sees
# plain peptides: e.g. DASAVILVSDAR + "1,Asp->Asn[D]" -> NASAVILVSDAR with the
# SAP token removed from the modification. Ambiguous Xle tokens are kept as-is.
# 返回 list(seq, mod)：代入后的序列与去掉 SAP token 的修饰串（无法代入时原样返回）。
# Returns list(seq, mod): substituted sequence and mod string without the SAP
# tokens (returns inputs unchanged when no unambiguous substitution applies).
apply_sap_to_sequence <- function(sq, md) {
  AA3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q",
           Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K",
           Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
           Tyr = "Y", Val = "V")
  if (is.na(md) || !nzchar(md) || !grepl("->", md, fixed = TRUE))
    return(list(seq = sq, mod = md))
  toks <- trimws(strsplit(md, ";", fixed = TRUE)[[1]])
  toks <- toks[nzchar(toks)]
  keep <- rep(TRUE, length(toks)); changed <- FALSE
  for (ti in seq_along(toks)) {
    tt <- toks[ti]
    m <- regmatches(tt, regexpr("^([0-9]+),([A-Za-z]+)->([A-Za-z]+)[[]", tt))
    if (length(m) == 0) next
    parts <- strsplit(sub("[[]$", "", m[[1]]), ",", fixed = TRUE)[[1]]
    pos <- suppressWarnings(as.integer(parts[1]))
    ft <- strsplit(parts[2], "->", fixed = TRUE)[[1]]  # From->To 在 parts[2]
    if (length(ft) < 2) next
    fr <- ft[1]; to <- ft[2]
    if (is.na(pos) || pos < 1 || pos > nchar(sq)) next
    if (!(fr %in% names(AA3)) || !(to %in% names(AA3))) next
    if (fr == "Xle") next                                     # I/L 歧义不可代入
    if (substr(sq, pos, pos) != AA3[[fr]]) next               # 原残基不符则不动
    sq <- paste0(substr(sq, 1, pos - 1L), AA3[[to]], substr(sq, pos + 1L, nchar(sq)))
    keep[ti] <- FALSE; changed <- TRUE
  }
  if (!changed) return(list(seq = sq, mod = md))
  list(seq = sq, mod = paste(toks[keep], collapse = ";"))
}


# ---- SAP 变体枚举Xle->Y 歧义时分别生成 I 和 L 两条 -----------
# Enumerate SAP variants: an Xle->Y substitution is ambiguous (I or L), so
# generate BOTH entries for prediction; the caller keeps the higher cosine.
# 返回 list 每元素 list(seq, mod)；无 Xle 歧义时只含 1 个条目。
# Returns a list of list(seq, mod); a single entry when no Xle ambiguity exists.
apply_sap_variants <- function(sq, md) {
  base <- apply_sap_to_sequence(sq, md)
  if (!grepl("Xle->", md, fixed = TRUE)) return(list(base))
  AA3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q",
           Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K",
           Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
           Tyr = "Y", Val = "V")
  toks <- trimws(strsplit(md, ";", fixed = TRUE)[[1]])
  toks <- toks[nzchar(toks)]
  out <- list()
  for (sub in c("I", "L")) {                                     # I/L 两个替换
    sq2 <- sq; keep <- rep(TRUE, length(toks))
    for (ti in seq_along(toks)) {
      tt <- toks[ti]
      if (!grepl("Xle->", tt, fixed = TRUE)) next
      ft <- strsplit(sub("^[0-9]+,", "", tt), "->", fixed = TRUE)[[1]]
      to_aa <- sub("[[].*$", "", ft[2])                          # 去 [X] 后缀
      pos <- suppressWarnings(as.integer(sub(",.*$", "", tt)))
      if (is.na(pos) || pos < 1 || pos > nchar(sq2)) next
      if (!(to_aa %in% names(AA3))) next
      sq2 <- paste0(substr(sq2, 1, pos - 1L), AA3[[to_aa]],
                    substr(sq2, pos + 1L, nchar(sq2)))
      keep[ti] <- FALSE
    }
    out[[length(out) + 1]] <- list(seq = sq2, mod = paste(toks[keep], collapse = ";"))
  }
  out
}
