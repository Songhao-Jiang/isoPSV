# =============================================================================
# stageA_classify.R - stage A detailed result classification
#
# 对扩展 report_A.tsv 的每一行（全部谱图，非仅冲突谱）给出逐谱结论状态：
#   无竞争匹配（各策略一致） | 存在竞争但原始结果正确 | 原始结果错误
#   （附错误原因与修饰差异细节） | 混合谱 | 谱图质量差难区分 |
#   子离子完全一致难区分
# 纯表格运算（无 MGF 加载）；错误行按图 4 口径分类：修饰改变（给出差异的
# 修饰种类与位置）、等重氨基酸替换或顺序互换、序列改变等。
#
# Labels EVERY row of the extended report_A.tsv (all spectra, not just the
# conflicts) with a per-spectrum conclusion status:
#   no competition (all strategies agree) | competition but original correct |
#   original result wrong (with error reason and modification-difference
#   detail) | chimeric | poor-quality ambiguous | identical-fragment-ions
#   ambiguous
# Pure table computation (no MGF loading). Wrong rows are classified in the
# Fig. 4 style: modification change (with the differing modification names and
# positions), same-weight AA replacement/permutation, sequence change, etc.
#
# Environment variables: EXT_REPORT, OUT.
# =============================================================================
# 核心库统一入口按 --file= 定位代码目录、绝对路径加载
# （boot.R），替代裸相对 source——此前仅当启动 cwd 恰为 code/ 时成立。
# core-library entry locate the code dir via --file= and
# load by absolute path (boot.R), replacing the bare relative source().
sd0 <- dirname(tryCatch(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"),
    error = function(e) ""))
if (is.na(sd0) || !dir.exists(sd0)) sd0 <- getwd()
source(file.path(sd0, "boot.R"))

ext_f <- Sys.getenv("EXT_REPORT")
out_f <- Sys.getenv("OUT")

aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# 混合谱判定需要逐 raw 无标记 MGF（UNLABELED_DIR 提供；为空则退回报告列口径）
# chimeric detection needs the per-raw unlabeled MGFs (via UNLABELED_DIR;
# when empty, fall back to the report-column definition)
ul_dir <- Sys.getenv("UNLABELED_DIR")
lu <- if (nzchar(ul_dir)) load_unlabeled_spectra(ul_dir) else NULL
if (!is.null(lu)) cat("unlabeled spectra loaded for chimeric check\n")

# 解析修饰串（"2,Oxidation[M];" -> "2,Oxidation[M]" 向量）
parse_mods <- function(s) {
  s <- as.character(s)
  if (!nzchar(s) || is.na(s)) return(character(0))
  parts <- strsplit(s, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

# 两修饰串的对称差异（英文：original extra X / winner extra Y）
mod_diff_detail <- function(a, b) {
  pa <- parse_mods(a)
  pb <- parse_mods(b)
  da <- setdiff(pa, pb)
  db <- setdiff(pb, pa)
  paste(c(if (length(da) > 0) paste0("original extra: ", paste(da, collapse = ";")),
          if (length(db) > 0) paste0("winner extra: ", paste(db, collapse = ";"))),
        collapse = " | ")
}

# 两个候选的理论 b/y 离子集合是否完全一致
identical_ion_sets <- function(s1, m1, s2, m2) {
  if (s1 == "" || s2 == "") return(FALSE)
  t1 <- try(suppressWarnings(invisible(capture.output(
    r1 <- by_ions_position_memo(s1, m1, "Unlabel", aa, aa_tab)))), silent = TRUE)
  t2 <- try(suppressWarnings(invisible(capture.output(
    r2 <- by_ions_position_memo(s2, m2, "Unlabel", aa, aa_tab)))), silent = TRUE)
  if (inherits(t1, "try-error") || inherits(t2, "try-error")) return(FALSE)
  if (length(r1) == 1 || length(r2) == 1) return(FALSE)
  a <- as.matrix(r1$ions_table[, -1, drop = FALSE])
  b <- as.matrix(r2$ions_table[, -1, drop = FALSE])
  identical(dim(a), dim(b)) &&
    isTRUE(all.equal(a, b, tolerance = 1e-6, check.attributes = FALSE))
}

# 序列的氨基酸组成是否一致（I 归一为 L，与生产口径一致）
# 直接比较排序后的字母向量：与出现顺序无关，顺序互换视为同组成
same_composition <- function(s1, s2) {
  comp <- function(s) sort(strsplit(gsub("I", "L", as.character(s)), "", fixed = TRUE)[[1]])
  identical(comp(s1), comp(s2))
}

# 理论离子缓存by_ions_position 是纯函数，同
# (seq|mod|label) 的候选在 5.7 万谱中大量重复（identical_ion_sets 每次两调）。
# 进程内哈希缓存命中直接复用返回值——语义不变（缓存只存成功结果，错误仍走
# try 路径）。ion_cache：哈希环境；key：seq|mod|label 拼接键。
ion_cache <- new.env(hash = TRUE, parent = emptyenv())
by_ions_position_memo <- function(seq, mod, label, aa, aa_tab) {
  key <- paste(seq, mod, label, sep = "|")
  if (exists(key, envir = ion_cache, inherits = FALSE))
    return(get(key, envir = ion_cache))
  val <- by_ions_position(seq, mod, label, aa, aa_tab)
  assign(key, val, envir = ion_cache)
  val
}

ext <- read.delim(ext_f, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0))
# 全空字符列会被 read.delim 读成逻辑 NA（空字段规则），恢复为空字符串，
# 否则 NA == "" 的比较产生 NA 并导致判定崩溃
# fully-empty character columns come back as logical NA (read.delim blank-field
# rule); restore them to empty strings so NA == "" cannot poison the judgment
ext[] <- lapply(ext, function(x) if (is.logical(x)) rep("", length(x)) else x)
cat("report_A rows:", nrow(ext), "\n")

n <- nrow(ext)
# 多引擎版：槽列名 a.Sequence 等（字母开头，无需 X 前缀），NSLOT 由报告列动态推断
# multi-engine edition: slot columns are a.Sequence etc. (letter-leading, no X
# prefix); NSLOT is inferred dynamically from the report columns
slot_seq_cols <- grep("^[a-z]\\.Sequence$", colnames(ext), value = TRUE)
slot_mod_cols <- grep("^[a-z]\\.Modification$", colnames(ext), value = TRUE)

# ---- v1.0 优化原实现每行 5 次 strsplit + 每行 data.frame
# 再 5.7 万个 rbind（实测 65–109s 的主因）。改为：
# ① 循环外一次性向量化 strsplit（cand/sr/cov/by 四个列表）；② 槽位两矩阵
#    sq_mat/md_mat 一次性 as.character；③ 7 个结论列预分配，末次一次 cbind。
# 逐行剩余逻辑与原来逐语句等价（仅把 ext[i, cn] 换成预取矩阵/列表下标）。
cand_lists <- if ("candidates" %in% names(ext))
  strsplit(as.character(ext$candidates), "|", fixed = TRUE) else
  replicate(n, character(0), simplify = FALSE)
num_lists <- function(col) {
  if (col %in% names(ext)) {
    lapply(strsplit(as.character(ext[[col]]), "|", fixed = TRUE),
           function(v) suppressWarnings(as.numeric(v)))
  } else replicate(n, numeric(0), simplify = FALSE)
}
sr_lists <- num_lists("candidate_specific_ratio")
cov_lists <- num_lists("candidate_coverage")
by_lists <- num_lists("candidate_PSM_score_by")
wins <- suppressWarnings(as.integer(ext$best_group))
# 注意：sqs 继承槽列名、x1s 因此带名（v1.0 对拍保留的旧语义）——P2-1
# 修复后，判定处用 identical(unname(x1s),
# ws) 只比值，同序列仅修饰不同的行正确进入 modification_change 分支
# （modification_diff 列自此生效）。矩阵仍 setNames 保留行名（供槽位取值）。
sq_mat <- if (length(slot_seq_cols)) do.call(rbind, setNames(
  lapply(slot_seq_cols, function(cn) as.character(ext[[cn]])), slot_seq_cols)) else NULL
md_mat <- if (length(slot_mod_cols)) do.call(rbind, setNames(
  lapply(slot_mod_cols, function(cn) as.character(ext[[cn]])), slot_mod_cols)) else NULL
i_a_seq <- match("a.Sequence", slot_seq_cols)
i_a_mod <- match("a.Modification", slot_mod_cols)
err_reasons <- if ("Error_reason" %in% names(ext)) as.character(ext$Error_reason) else NULL

conclusion_status <- rep("no competition (strategies agree)", n)
vote_says_wrong <- integer(n); error_reason <- character(n)
modification_diff <- character(n); chimeric <- integer(n)
identical_ions <- integer(n); quality_low <- integer(n)

# ---- v1.0分类循环并行化 ＋ 重建 ----
# 事故：v1.0 重写版把数据前导（env 读取 / read.delim / n <- nrow(ext) /
# cand_lists 等构建）整体丢失，独立运行时自动判定行立即引用未定义的 n，
# Stage A 秒崩（v1.0 全量测试 06:47 实测，Error: object 'n' not found）；
# 其"顺序回退"还 source 了不存在的 stageA_classify_sequential.R。
# 本次重建＝v1.0 原版前导（含 boot.R 加载与全部辅助函数）＋ v1.0
# 并行分支 ＋ v1.0 原版顺序循环 ＋ v1.0 的 X/U 标记输出尾。
# 并行 worker 体与原版循环体经规范化 diff 逐语句等价（仅 cl→cl_i 重命名
# 避免遮蔽集群对象）。
# ---- v1.0 parallel classify loop; rebuilt ----
# The v1.0 rewrite accidentally dropped the whole data-prep prologue, so
# the standalone run crashed instantly on an undefined `n` (live at 06:47 in
# the v1.0 full test) and its sequential fallback sourced a nonexistent
# stageA_classify_sequential.R. This rebuild = the v1.0 original prologue
# (boot.R load + all helpers) + the v1.0 parallel branch + the original
# sequential loop + the v1.0 X/U-marking output tail. The parallel worker
# body is statement-for-statement equivalent to the original loop (normalized
# diff; only cl -> cl_i renamed to avoid shadowing the cluster object).
stageA_workers_val <- 1L
wc <- Sys.getenv("STAGEA_WORKERS", unset = "")
if (nzchar(wc)) {
  stageA_workers_val <- suppressWarnings(as.integer(wc))
  if (is.na(stageA_workers_val) || stageA_workers_val < 1) stageA_workers_val <- 1L
}
# 自动判定：大表（>10K 行）且多核时启用 5 worker
if (stageA_workers_val == 1L && n > 10000) {
  stageA_workers_val <- min(5L, parallel::detectCores() - 1L)
  cat("[Stage A 并行] 大表（", n, " 行），自动启用 ", stageA_workers_val, " worker\n", sep = "")
  cat("[Stage A parallel] large table (", n, " rows), auto-enabling ", stageA_workers_val, " workers\n", sep = "")
}
# results 初始化必须在并行 if 之前——小表（≤10K 行）不进并行分支，
# 若无此初始化，下方 is.null(results) 直接 object not found（
# run_all 实测：12 行夹具秒崩；大表因并行分支赋值而掩盖）
# initialize results BEFORE the parallel branch — a small table (≤10K rows)
# never enters it, and without the init the is.null(results) gate below dies
# with object not found (live in run_all: a 12-row fixture crashed
# instantly; large tables masked it because the parallel branch assigns)
results <- NULL
if (stageA_workers_val > 1 && n > 10000 && requireNamespace("parallel", quietly = TRUE)) {
  # ---- 并行路径 ----
  chunk_size <- ceiling(n / (stageA_workers_val * 4))  # 4× 过采样给 LB 均衡
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  cat("[Stage A 并行] ", length(chunks), " chunks, ", stageA_workers_val, " workers\n", sep = "")

  cl <- parallel::makeCluster(stageA_workers_val)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # 导出共享只读数据＋全部辅助函数（一次性序列化）。
  # v1.0 重建修正原显式清单漏传递依赖——worker 内
  # chimeric_check → title_raw_scan / title_scan_key（boot.R 经 ReadMGF_script.R
  # 载入）→ matched_peak_set → cached_by_ions_position → by_ions_position 一整条
  # 链不在清单里，parLapplyLB 10/20 chunk 报 "could not find function
  # "title_raw_scan"" 回退顺序。逐个枚举不可维护（改一处库函数就再漏），
  # 改为导出全环境（仅排除集群对象自身）。附带好处：每个 worker 得到独立
  # 的 ion_cache 环境，免共享竞争。
  # export the shared read-only data + ALL helper functions in one shot.
  # v1.0 rebuild fix the explicit list missed transitive
  # deps — inside the worker chimeric_check needs title_raw_scan /
  # title_scan_key (loaded via boot.R from ReadMGF_script.R), then
  # matched_peak_set -> cached_by_ions_position -> by_ions_position; 10/20
  # chunks died with "could not find function "title_raw_scan"" and fell back
  # to sequential. Enumerating names is unmaintainable (any library refactor
  # re-breaks it), so export the whole environment instead (minus the cluster
  # object itself). Bonus: every worker gets its own ion_cache environment,
  # so there is no shared-cache contention.
  parallel::clusterExport(cl, setdiff(ls(), "cl"))

  results <- tryCatch({
    parallel::parLapplyLB(cl, chunks, function(idx) {
      n_local <- length(idx)
      cs <- character(n_local); vw <- integer(n_local)
      er <- character(n_local); md <- character(n_local)
      ch <- integer(n_local); ii <- integer(n_local); ql <- integer(n_local)

      for (j in seq_along(idx)) {
        i <- idx[j]
        is_conf <- (ext$Peptide_state[i] != "Same") | (ext$Modification_state[i] != "Same")
        sqs <- if (length(slot_seq_cols)) sq_mat[, i] else character(0)
        mds <- if (length(slot_mod_cols)) md_mat[, i] else character(0)
        cl_i <- cand_lists[[i]]
        cand_sqs <- sub("\\(.*$", "", cl_i)
        has_mod <- grepl("\\(", cl_i)
        cand_mds <- ifelse(has_mod, sub("^.*\\((.*)\\)$", "\\1", cl_i), "")
        cand_mds <- rep_len(cand_mds, length(cand_sqs))
        x1s <- if (length(sqs) >= 1 && !is.na(i_a_seq)) sqs[i_a_seq] else
          (if (length(cand_sqs) >= 1) cand_sqs[1] else "")
        x1m <- if (length(mds) >= 1 && !is.na(i_a_mod)) mds[i_a_mod] else
          (if (length(cand_mds) >= 1) cand_mds[1] else "")

        status <- "no competition (strategies agree)"
        vote_wrong <- 0L; err_reason <- ""; mod_diff <- ""
        chimeric_i <- 0L; ident_ions <- 0L; qual_low <- 0L

        if (is_conf && length(cand_sqs) >= 1) {
          win <- wins[i]
          if (is.na(win) || win < 1 || win > length(cand_sqs)) win <- 1L
          ws <- cand_sqs[win]
          wm <- cand_mds[win]
          vote_wrong <- as.integer(!(ws == x1s && wm == x1m))
          if (!is.null(lu)) {
            chimeric_i <- chimeric_check(cand_sqs, cand_mds, lu,
                                         as.character(ext$File_Name[i]), aa, aa_tab)
          } else {
            sr <- sort(sr_lists[[i]][!is.na(sr_lists[[i]])], decreasing = TRUE)
            if (length(sr) >= 2 && sr[1] > 0.1 && sr[2] > 0.1) chimeric_i <- 1L
          }
          cov_all <- cov_lists[[i]]
          cov <- if (length(cov_all) >= win && !is.na(cov_all[win])) cov_all[win] else NA_real_
          if (!is.na(cov) && cov < 20) qual_low <- 1L
          byv <- by_lists[[i]]
          others <- setdiff(which(!is.na(byv) & seq_along(byv) <= length(cand_sqs)), win)
          if (length(others) > 0) {
            runner <- others[which.max(byv[others])]
            rs <- cand_sqs[runner]
            rm_ <- cand_mds[runner]
            if (ws != "" && rs != "" && ws != rs) {
              ident_ions <- as.integer(identical_ion_sets(ws, wm, rs, rm_))
            }
          }
          if (chimeric_i == 1L) {
            status <- "chimeric spectrum"
          } else if (ident_ions == 1L) {
            status <- "ambiguous: identical fragment ions"
          } else if (qual_low == 1L) {
            status <- "ambiguous: poor spectrum quality"
          } else if (vote_wrong == 1L) {
            status <- "original result wrong"
          } else {
            status <- "competition but original correct"
          }
          if (vote_wrong == 1L) {
            if (!is.null(err_reasons) && nzchar(err_reasons[i])) {
              err_reason <- err_reasons[i]
            }
            if (ws != "" && identical(unname(x1s), ws)) {
              if (x1m != wm) {
                if (!nzchar(err_reason)) err_reason <- "modification_change"
                mod_diff <- mod_diff_detail(x1m, wm)
              }
            } else if (ws != "" && x1s != "") {
              if (!nzchar(err_reason)) {
                err_reason <- if (same_composition(x1s, ws)) {
                  "Same-weight AA replacement or permutation"
                } else {
                  "sequence_change"
                }
              }
            }
          }
        }

        cs[j] <- status; vw[j] <- vote_wrong; er[j] <- err_reason
        md[j] <- mod_diff; ch[j] <- chimeric_i; ii[j] <- ident_ions; ql[j] <- qual_low
      }

      list(cs = cs, vw = vw, er = er, md = md, ch = ch, ii = ii, ql = ql)
    })
  }, error = function(e) {
    cat("[Stage A 并行] 失败，回退顺序模式：", conditionMessage(e), "\n")
    cat("[Stage A parallel] failed, falling back to sequential: ", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(results)) {
    conclusion_status <- unlist(lapply(results, `[[`, "cs"), use.names = FALSE)
    vote_says_wrong   <- unlist(lapply(results, `[[`, "vw"), use.names = FALSE)
    error_reason      <- unlist(lapply(results, `[[`, "er"), use.names = FALSE)
    modification_diff <- unlist(lapply(results, `[[`, "md"), use.names = FALSE)
    chimeric          <- unlist(lapply(results, `[[`, "ch"), use.names = FALSE)
    identical_ions    <- unlist(lapply(results, `[[`, "ii"), use.names = FALSE)
    quality_low       <- unlist(lapply(results, `[[`, "ql"), use.names = FALSE)
    cat("[Stage A 并行] 完成，合并 ", length(conclusion_status), " 行\n", sep = "")
  } else {
    results <- NULL
  }
}

if (is.null(results)) {
  # ---- 顺序路径：v1.0 原版循环，零行为变化（并行失败/小表/单 worker）----
  # sequential path: the original v1.0 loop, zero behavior change (taken
  # when parallel failed / small table / single worker)
for (i in seq_len(n)) {
  is_conf <- (ext$Peptide_state[i] != "Same") | (ext$Modification_state[i] != "Same")
  # 多策略候选（候选列表＝非空策略槽按字母序排列；候选 1＝策略 a 的原始结果）
  sqs <- if (length(slot_seq_cols)) sq_mat[, i] else character(0)
  mds <- if (length(slot_mod_cols)) md_mat[, i] else character(0)
  # 审查修正：best_group 是"去重候选列表"的下标（含 I/L 归一与同
  # 内容槽合并），槽位列表与它不同序不同长——候选一律取 candidates 列（与
  # best_group 同源同序）；a 槽内容仅作原始结果参照
  # best_group indexes the DEDUPED candidate list; slots are a different list.
  # Candidates therefore come from the candidates column (same list, same order)
  cl <- cand_lists[[i]]
  cand_sqs <- sub("\\(.*$", "", cl)
  has_mod <- grepl("\\(", cl)
  cand_mds <- ifelse(has_mod, sub("^.*\\((.*)\\)$", "\\1", cl), "")
  cand_mds <- rep_len(cand_mds, length(cand_sqs))
  x1s <- if (length(sqs) >= 1 && !is.na(i_a_seq)) sqs[i_a_seq] else
    (if (length(cand_sqs) >= 1) cand_sqs[1] else "")
  x1m <- if (length(mds) >= 1 && !is.na(i_a_mod)) mds[i_a_mod] else
    (if (length(cand_mds) >= 1) cand_mds[1] else "")

  status <- "no competition (strategies agree)"
  vote_wrong <- 0L; err_reason <- ""; mod_diff <- ""
  chimeric_i <- 0L; ident_ions <- 0L; qual_low <- 0L

  if (is_conf && length(cand_sqs) >= 1) {
    # 阶段 A 多数投票最佳组（单一值，候选列表索引）
    win <- wins[i]
    if (is.na(win) || win < 1 || win > length(cand_sqs)) win <- 1L
    ws <- cand_sqs[win]
    wm <- cand_mds[win]
    # 内容口径：胜选候选与原始结果序列＋修饰一致即视为正确
    vote_wrong <- as.integer(!(ws == x1s && wm == x1m))
    # 混合谱：不同序列两两"各自特异子离子比例的最小值 > 10%"（设计约定，
    # m/z 口径：特异峰数 ÷ 两候选命中峰并集；修订）
    # chimeric: any pair of distinct sequences with min of the two
    # specific-ion ratios > 10% (the user's definition; m/z basis: specific
    # peaks over the union of both candidates' matched peaks; revised )
    if (!is.null(lu)) {
      chimeric_i <- chimeric_check(cand_sqs, cand_mds, lu,
                                 as.character(ext$File_Name[i]), aa, aa_tab)
    } else {
      sr <- sort(sr_lists[[i]][!is.na(sr_lists[[i]])], decreasing = TRUE)
      if (length(sr) >= 2 && sr[1] > 0.1 && sr[2] > 0.1) chimeric_i <- 1L
    }
    # 谱图质量差：选中候选覆盖度 < 20%（面板 e 校准的判别阈值）
    cov_all <- cov_lists[[i]]
    cov <- if (length(cov_all) >= win && !is.na(cov_all[win])) cov_all[win] else NA_real_
    if (!is.na(cov) && cov < 20) qual_low <- 1L
    # 子离子完全一致：胜选候选 vs 其余候选中最优者（b/y 打分最高）
    byv <- by_lists[[i]]
    others <- setdiff(which(!is.na(byv) & seq_along(byv) <= length(cand_sqs)), win)
    if (length(others) > 0) {
      runner <- others[which.max(byv[others])]
      rs <- cand_sqs[runner]
      rm_ <- cand_mds[runner]
      if (ws != "" && rs != "" && ws != rs) {
        ident_ions <- as.integer(identical_ion_sets(ws, wm, rs, rm_))
      }
    }
    # 结论状态（优先级：混合谱 > 子离子一致 > 质量差 > 错误 > 正确；英文输出）
    if (chimeric_i == 1L) {
      status <- "chimeric spectrum"
    } else if (ident_ions == 1L) {
      status <- "ambiguous: identical fragment ions"
    } else if (qual_low == 1L) {
      status <- "ambiguous: poor spectrum quality"
    } else if (vote_wrong == 1L) {
      status <- "original result wrong"
    } else {
      status <- "competition but original correct"
    }
    # 错误行：按图 4 口径给错误原因与修饰差异细节
    if (vote_wrong == 1L) {
      if (!is.null(err_reasons) && nzchar(err_reasons[i])) {
        err_reason <- err_reasons[i]
      }
      # P2-1 修复x1s 取自
      # 槽位矩阵带列名属性，identical() 默认比属性 → 恒 FALSE，同序列仅修饰
      # 不同的行永远掉进 same_composition 分支被误标"Same-weight AA
      # replacement"、modification_diff 死列。unname 后只比值——同序列命中
      # 第 1 层：标 modification_change 并填充修饰差异详情
      if (ws != "" && identical(unname(x1s), ws)) {
        if (x1m != wm) {
          if (!nzchar(err_reason)) err_reason <- "modification_change"
          mod_diff <- mod_diff_detail(x1m, wm)
        }
      } else if (ws != "" && x1s != "") {
        if (!nzchar(err_reason)) {
          err_reason <- if (same_composition(x1s, ws)) {
            "Same-weight AA replacement or permutation"
          } else {
            "sequence_change"
          }
        }
      }
    }
  }

  conclusion_status[i] <- status
  vote_says_wrong[i] <- vote_wrong
  error_reason[i] <- err_reason
  modification_diff[i] <- mod_diff
  chimeric[i] <- chimeric_i
  identical_ions[i] <- ident_ions
  quality_low[i] <- qual_low
}

}
# ---- 输出组装 ----
out <- ext
out$conclusion_status <- conclusion_status
out$vote_says_wrong <- vote_says_wrong
out$error_reason <- error_reason
out$modification_diff <- modification_diff
out$chimeric <- chimeric
out$identical_ions <- identical_ions
out$quality_low <- quality_low

# v1.0零质量/歧义残基标记列
# 含 X（Xle，I/L 歧义）或 U（Sec，此前 pFind 赋零质量）的肽段行加标记
# mark rows containing X (Xle) or U (Sec) residues
has_XU <- grepl("X", ext$candidates %||% "") | grepl("U", ext$candidates %||% "")
out$ambiguous_residue <- ifelse(has_XU, "X_or_U", "")
if (any(has_XU)) {
  cat("[Stage A] 含 X/U 残基行: ", sum(has_XU), "（详见 ambiguous_residue 列）\n", sep = "")
  cat("[Stage A] rows with X/U residues: ", sum(has_XU), " (see ambiguous_residue column)\n", sep = "")
}

write.table(out, out_f, sep = "\t", row.names = FALSE, quote = FALSE)
cat("stage A classification written:", out_f, "(", nrow(out), "rows )\n")
cat("=== conclusion status distribution 结论状态分布 ===\n")
print(table(out$conclusion_status))
cat("=== wrong rows by error reason 错误原因分布（原始结果错误判定行） ===\n")
print(table(out$error_reason[out$vote_says_wrong == 1]))
