# =============================================================================
# chimeric_check.R - 混合谱判定（用户口径，2026-09-06 m/z 口径修订）
# "两个序列的各自特异性子离子比例的最小值 > 10%" 判为混合谱：
# 对每个竞争序列，其特异子离子比例 = 只被该序列命中的观测峰数 ÷
# 两个序列命中观测峰的并集数（按观测峰 m/z 判断，不看肽段坐标位点——
# 同一碎片在两序列中位点不同但命中同一峰时按同一离子计）；两个序列各自的
# 该比例都 > 10%（即最小值 > 10%）才算混合谱。
#
# Chimeric detection (the user's definition, revised to the m/z basis on
# 2026-09-06): a spectrum is chimeric when the MINIMUM of the two competing
# sequences' own specific-ion ratios is > 10%. Each sequence's specific-ion
# ratio = the observed peaks matched ONLY by it divided by the UNION of the
# observed peaks matched by the two sequences (judged by observed peak m/z,
# not by position coordinate on the peptide - the same fragment that sits at
# different positions in two sequences counts as ONE ion when it hits the
# same peak); both ratios must be > 10%.
#
# 依赖：compat.R / acs_core.R / paired_params.R / ReadMGF_script.R 已 source，
# aa、aa_tab 已就绪。用法见 load_unlabeled_spectra / chimeric_top2。
# Requires compat.R / acs_core.R / paired_params.R / ReadMGF_script.R to be
# sourced and aa / aa_tab to be ready. See load_unlabeled_spectra / chimeric_top2.
# =============================================================================

# 载入逐 raw 无标记 MGF（文件夹内平铺的 per-raw mgf，排除合并文件），
# 返回 list(objs, titles)；目录为空返回 NULL。
# Load the per-raw unlabeled MGFs (flat files in one folder, merged file
# excluded); returns list(objs, titles), or NULL when the folder is empty.
load_unlabeled_spectra <- function(unlabeled_dir, only = NULL) {
  # only：可选 raw 茎向量——只加载需要的 raw（主表导出仅需嵌合谱所在 raw；
  # 全文件夹 7 组分全量加载实测 4-5 分钟、内存 6GB+，2026-09-11 提效）
  # only: optional raw-stem vector - load just the needed raws (the master
  # export needs only the raws of chimeric spectra; loading the whole folder
  # measured 4-5 min and 6 GB+ RSS, optimized 2026-09-11)
  files <- sort(list.files(unlabeled_dir, pattern = "\\.mgf$", full.names = TRUE))
  files <- files[basename(files) != "unlabeled.mgf"]
  # 统一排除 pFind 残留 *_check.mgf（无 RTINSECONDS；全量检测报告附带发现，
  # 2026-09-26——与 pFind 模板改写路径、denovo 调用前移出同口径）
  files <- files[!grepl("_check\\.mgf$", basename(files))]
  if (!is.null(only) && length(files) > 0) {
    stems <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(files)))
    files <- files[stems %in% only]
  }
  if (length(files) == 0) return(NULL)
  # 复用进程级缓存集群：不重复建群、不重复 source 260KB 库（40-60 s/次省掉；
  # 集群对象随进程结束释放，不再每调用 stopCluster）
  # reuse the process-level cached cluster: no per-call rebuild and no per-call
  # re-source of the 260KB library (40-60 s saved per call; the cluster lives
  # until the process ends instead of being stopped on every call)
  cl <- mgf_read_cluster()
  loaded <- load_mgf_shared(files, cl)
  parallel::clusterCall(cl, gc)   # worker 垃圾回收，长驻集群不越积越多
                                   # worker gc keeps the resident cluster lean
  # 键用去 _HCDFT 后缀的文件茎（与标题解析出的 raw 名一致：raw.scan.scan.价态.拆分.dta）
  # key by the file stem with the _HCDFT suffix stripped (matching the raw name
  # parsed from the title: raw.scan.scan.charge.split.dta)
  names(loaded) <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(files)))
  titles <- lapply(loaded, function(objs)
    vapply(objs, function(o) o@title, character(1)))
  # 标题索引（v1.0 优化，2026-09-26）：chimeric_check 逐冲突行
  # match(fn, titles[[raw]]) 是 O(谱数) 线性扫——每 raw 预建命名向量索引，
  # 查谱 O(1)；ix 与 titles 同构（list[raw] = named int vector）。
  ix <- lapply(titles, function(t) setNames(seq_along(t), t))
  list(objs = loaded, titles = titles, ix = ix)
}

# 单个候选的匹配槽位集合（20 ppm 窗口、排序＋前缀和，与 ions_match 同口径；
# 修饰不被质量表识别时按无修饰重试一次）
# the matched b/y slot set of one candidate (20 ppm windows via sorting +
# prefix sums, same definition as ions_match; retries without the modification
# when the mass table rejects it)
matched_positions_of <- function(seq, mod, sp, aa, aa_tab, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- cached_by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) ||
        length(bip) == 1 || is.null(bip$ions_table)) next
    t <- bip$ions_table
    mz <- as.numeric(sp[[1]]); ints <- as.numeric(sp[[2]])
    ord <- order(mz, na.last = NA)
    mz <- mz[ord]; ints <- ints[ord]
    cum2 <- c(0, cumsum(ints))
    mz_all <- c(t$b_1_charge, t$b_2_charge, t$y_1_charge, t$y_2_charge)
    # 键位口径（2026-09-15 修复）：b_j = 键 j；y_j = 键 n-j（n=肽长）。
    # 旧版对 y 系列复用离子序号 j，把 y2/y3/y4 误标为键 2/3/4、与 b2/b3/b4
    # 混为一谈——C 端连续缺口（如 y1 与 y5 之间的 4 残基）被 N 端 b 离子
    # 错误补位，consecutive_missing 规则系统性漏判。
    # bond-wise positions (fixed 2026-09-15): b_j is bond j; y_j is bond n-j
    # (n = peptide length). The old code reused the ion number j for the y
    # series, mislabeling y2/y3/y4 as bonds 2/3/4 and confusing them with
    # b2/b3/b4 - C-terminal consecutive gaps (e.g. the four residues between
    # y1 and y5) were wrongly covered by N-terminal b ions, so the
    # consecutive_missing rule systematically missed them.
    pos_all <- c(t$Number, t$Number,
                 nchar(seq) - t$Number, nchar(seq) - t$Number)
    lo <- mz_all * (1 - tol_ppm * 1e-6)
    hi <- mz_all * (1 + tol_ppm * 1e-6)
    i_lo <- findInterval(lo, mz)
    i_hi <- findInterval(hi, mz)
    # 先清 NA（未知修饰质量 → findInterval 得 NA，下标赋值不允许 NA）
    # clear NA first (an unknown modification mass gives NA from findInterval;
    # NA is not allowed in subscripted assignment)
    na_flag <- is.na(i_lo) | is.na(i_hi)
    i_lo[na_flag] <- 0L
    i_hi[na_flag] <- 0L
    mi <- cum2[i_hi + 1] - cum2[i_lo + 1]
    mi[i_hi <= i_lo] <- 0
    return(unique(pos_all[mi > 0]))
  }
  integer(0)
}

# 对一条谱，对全部不同序列两两判定"各自特异子离子比例的最小值 ≥ 10%"；
# 任一序列对满足即判为混合谱（每条不同序列只匹配一次）。
# For one spectrum, check every pair of DISTINCT sequences by "min of the two
# specific-ion ratios > 10%" (m/z basis, union denominator); any pair
# satisfying it marks the spectrum chimeric.
chimeric_check <- function(seqs, mods, lu, fn, aa, aa_tab) {
  if (is.null(lu) || length(seqs) < 2) return(0L)
  raw <- unname(title_raw_scan(fn)["raw"])
  objs <- lu$objs[[raw]]
  if (is.null(objs)) return(0L)
  # v1.0 优化：优先走 O(1) 命名向量索引（load_unlabeled_spectra 新返回
  # 的 lu$ix），旧 match 线性扫保留为回退（无 ix 时语义不变）
  li <- if (!is.null(lu$ix) && !is.null(lu$ix[[raw]])) {
    unname(lu$ix[[raw]][fn])
  } else {
    match(fn, lu$titles[[raw]])
  }
  if (is.na(li)) {
    # raw+scan 兜底（价态/拆分序号可能不同）
    key <- title_scan_key(fn)
    li <- which(startsWith(lu$titles[[raw]], paste0(key, ".")))
    if (length(li) == 0) return(0L)
    li <- li[1]
  }
  sp <- objs[[li]]@spectrum
  uniq <- which(!duplicated(seqs) & seqs != "")
  if (length(uniq) < 2) return(0L)
  # 各不同序列命中的观测峰集合（m/z 口径；同一峰不论位点只计一次）
  # the observed-peak sets matched by each distinct sequence (m/z basis; a
  # shared peak counts once regardless of its position in either peptide)
  peaks <- lapply(seq_along(uniq), function(u) {
    k <- uniq[u]
    matched_peak_set(seqs[k], mods[k], sp, aa, aa_tab)
  })
  for (a in seq_along(peaks)) {
    for (b in seq_len(a - 1)) {
      A <- peaks[[a]]; B <- peaks[[b]]
      if (length(A) == 0 || length(B) == 0) next
      U <- union(A, B)
      rA <- length(setdiff(A, B)) / length(U)
      rB <- length(setdiff(B, A)) / length(U)
      if (min(rA, rB) > 0.1) return(1L)
    }
  }
  0L
}
