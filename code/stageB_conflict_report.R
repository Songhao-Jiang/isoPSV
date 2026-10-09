# =============================================================================
# stageB_conflict_report.R - conflict-spectrum conclusion annotation
#
# 在扩展 report_A.tsv 的基础上提取"冲突谱"行（Peptide_state 或
# Modification_state != Same），逐行合并阶段 B 配对报告的新增列，并给出
# 逐谱结论状态：
#   混合谱 | 子离子完全一致难区分 | 投票并列难区分 | 谱图质量差难区分 |
#   原始结果错误（附错误原因与修饰差异细节） | 存在竞争但原始结果正确 |
#   未进行成对核查（无谱/单候选被阶段 B 剔除）
# 纯表格运算（无 MGF 加载）；错误行按图 4 口径分类：修饰改变（给出差异的
# 修饰种类与位置）、等重氨基酸替换或顺序互换、序列改变、N 端非特异酶切等。
#
# Extracts the conflict rows (Peptide_state or Modification_state != Same)
# from the extended report_A.tsv, merges the stage-B paired-report columns
# onto each row, and labels every spectrum with a conclusion status:
#   chimeric | identical-fragment-ions ambiguous | vote-tie ambiguous |
#   poor-quality ambiguous | original result wrong (with the error reason
#   and the modification-difference detail) | original result correct |
#   not pair-checked (dropped by stage B: no spectrum / single candidate)
# Pure table computation (no MGF loading). Wrong rows are classified in the
# Fig. 4 style: modification change (with the differing modification names and
# positions), same-weight AA replacement/permutation, sequence change,
# N-terminal non-specific cleavage, etc.
#
# Environment variables: EXT_REPORT, PAIRED, OUT.
# =============================================================================
# 核心库统一入口（N-1，2026-09-25）：按 --file= 定位代码目录、绝对路径加载
# （boot.R），替代裸相对 source——此前仅当启动 cwd 恰为 code/ 时成立。
# core-library entry (N-1, 2026-09-25): locate the code dir via --file= and
# load by absolute path (boot.R), replacing the bare relative source().
sd0 <- dirname(tryCatch(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"),
    error = function(e) ""))
if (is.na(sd0) || !dir.exists(sd0)) sd0 <- getwd()
source(file.path(sd0, "boot.R"))

ext_f <- Sys.getenv("EXT_REPORT")
paired_f <- Sys.getenv("PAIRED")
out_f <- Sys.getenv("OUT")

aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# 混合谱判定需要逐 raw 无标记 MGF（UNLABELED_DIR 提供；为空则退回报告列口径）
ul_dir <- Sys.getenv("UNLABELED_DIR")
lu <- if (nzchar(ul_dir)) load_unlabeled_spectra(ul_dir) else NULL
if (!is.null(lu)) cat("unlabeled spectra loaded for chimeric check\n")

# 解析修饰串（"2,Oxidation[M];" -> "2,Oxidation[M]" 向量）
# parse a modification string ("2,Oxidation[M];" -> vector of "pos,name")
parse_mods <- function(s) {
  s <- as.character(s)
  if (!nzchar(s) || is.na(s)) return(character(0))
  parts <- strsplit(s, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

# 两修饰串的对称差异（原结果多 X / 胜选候选多 Y）
# symmetric difference of two modification strings
mod_diff_detail <- function(a, b) {
  pa <- parse_mods(a)
  pb <- parse_mods(b)
  da <- setdiff(pa, pb)
  db <- setdiff(pb, pa)
  paste(c(if (length(da) > 0) paste0("original extra: ", paste(da, collapse = ";")),
          if (length(db) > 0) paste0("winner extra: ", paste(db, collapse = ";"))),
        collapse = " | ")
}

# 两个候选的理论 b/y 离子集合是否完全一致（等重替换如 I/L 会产生完全相同的碎片）
# whether two candidates have identical theoretical b/y ion sets (same-weight
# replacements such as I/L produce exactly the same fragments)
identical_ion_sets <- function(s1, m1, s2, m2) {
  if (s1 == "" || s2 == "") return(FALSE)
  t1 <- try(suppressWarnings(invisible(capture.output(
    r1 <- by_ions_position(s1, m1, "Unlabel", aa, aa_tab)))), silent = TRUE)
  t2 <- try(suppressWarnings(invisible(capture.output(
    r2 <- by_ions_position(s2, m2, "Unlabel", aa, aa_tab)))), silent = TRUE)
  if (inherits(t1, "try-error") || inherits(t2, "try-error")) return(FALSE)
  if (length(r1) == 1 || length(r2) == 1) return(FALSE)
  a <- as.matrix(r1$ions_table[, -1, drop = FALSE])
  b <- as.matrix(r2$ions_table[, -1, drop = FALSE])
  identical(dim(a), dim(b)) &&
    isTRUE(all.equal(a, b, tolerance = 1e-6, check.attributes = FALSE))
}

# 序列的氨基酸组成是否一致（I 归一为 L，与生产口径一致）
# 直接比较排序后的字母向量：与出现顺序无关，顺序互换视为同组成
# identical amino-acid composition (I folded to L, per the production rule);
# compare the sorted letter vectors directly: order-independent
same_composition <- function(s1, s2) {
  comp <- function(s) sort(strsplit(gsub("I", "L", as.character(s)), "", fixed = TRUE)[[1]])
  identical(comp(s1), comp(s2))
}

ext <- read.delim(ext_f, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0))
# 全空字符列会被 read.delim 读成逻辑 NA（空字段规则），恢复为空字符串
# fully-empty character columns come back as logical NA (read.delim blank-field
# rule); restore them to empty strings
ext[] <- lapply(ext, function(x) if (is.logical(x)) rep("", length(x)) else x)
if (!file.exists(paired_f)) {
  stop("paired report not found: ", paired_f,
       "\n（请先运行阶段 B 生成 paired_report_<label>.tsv）\n(run stage B first)")
}
pr <- read.delim(paired_f, sep = "\t", stringsAsFactors = FALSE,
                 check.names = FALSE, na.strings = character(0))
pr[] <- lapply(pr, function(x) if (is.logical(x)) rep("", length(x)) else x)

conf <- ext[ext$Peptide_state != "Same" | ext$Modification_state != "Same", ,
            drop = FALSE]
# 终端丢失规范化（全部槽列 a.Sequence 等）：起始 M 丢失/C 端 R 丢失转为序列变化，
# 终端类位置 1 归 0——内容判定与输出列统一用规范化形式（2026-09-02 增）
# Terminal-loss normalization across the slot columns a.Sequence etc. so both
# the content-based calls and the output columns use the canonical form
# (added 2026-09-02)
slot_seq_cols <- grep("^[a-z]\\.Sequence$", names(conf), value = TRUE)
slot_mod_cols <- grep("^[a-z]\\.Modification$", names(conf), value = TRUE)
for (j in seq_along(slot_seq_cols)) {
  scn <- slot_seq_cols[j]; mcn <- slot_mod_cols[j]
  if (scn %in% names(conf)) {
    sq <- as.character(conf[[scn]])
    md <- if (mcn %in% names(conf)) as.character(conf[[mcn]]) else rep("", length(sq))
    md[is.na(md)] <- ""
    nz <- nzchar(sq)
    fx <- vapply(which(nz), function(k) {
      r <- normalize_terminal_loss(sq[k], md[k])
      paste(r$s, if (nzchar(r$m)) paste0("	", r$m) else "", sep = "")
    }, character(1))
    sq[nz] <- sub("	.*$", "", fx)
    md[nz] <- ifelse(grepl("	", fx), sub("^[^	]*	", "", fx), "")
    conf[[scn]] <- sq
    if (mcn %in% names(conf)) conf[[mcn]] <- md
  }
}
ip <- match(conf$File_Name, pr$File_Name)
cat("report_A rows:", nrow(ext), "| conflict rows:", nrow(conf),
    "| with paired result:", sum(!is.na(ip)), "\n")

# v1.0（2026-10-01）：冲突报告并行化（逐冲突谱独立）
# 冲突谱数 >10K 且多核时自动启用 parLapplyLB；否则走原顺序循环
conflict_workers_val <- 1L
cw <- Sys.getenv("CONFLICT_WORKERS", unset = "")
if (nzchar(cw)) {
  conflict_workers_val <- suppressWarnings(as.integer(cw))
  if (is.na(conflict_workers_val) || conflict_workers_val < 1) conflict_workers_val <- 1L
}
if (conflict_workers_val == 1L && nrow(conf) > 10000) {
  conflict_workers_val <- min(5L, parallel::detectCores() - 1L)
  cat("[冲突报告并行] 大表（", nrow(conf), " 行），自动启用 ", conflict_workers_val, " worker\n", sep = "")
  cat("[Conflict parallel] large table (", nrow(conf), " rows), auto-enabling ", conflict_workers_val, " workers\n", sep = "")
}

if (conflict_workers_val > 1 && nrow(conf) > 10000 && requireNamespace("parallel", quietly = TRUE)) {
  # ---- 并行路径 ----
  n_conf <- nrow(conf)
  chunk_size <- ceiling(n_conf / (conflict_workers_val * 4))
  chunks <- split(seq_len(n_conf), ceiling(seq_len(n_conf) / chunk_size))
  cat("[冲突报告并行] ", length(chunks), " chunks, ", conflict_workers_val, " workers\n", sep = "")

  cl <- parallel::makeCluster(conflict_workers_val)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # v1.0（2026-10-02，随 stageA 重建一并修）：显式清单漏 chimeric_check
  # 的传递依赖（title_raw_scan / title_scan_key / matched_peak_set /
  # cached_by_ions_position / by_ions_position 一整条链，均经 boot.R 载入）。
  # 大表触发并行时会整体失败并静默回退顺序（Macaca 规模＝白丢并行收益）。
  # 改导出全环境（排除集群对象），与 stageA 同口径。
  # v1.0 (2026-10-02, fixed alongside the stageA rebuild): the explicit
  # list missed chimeric_check's transitive deps (title_raw_scan /
  # title_scan_key / matched_peak_set / cached_by_ions_position /
  # by_ions_position, all loaded via boot.R). On a large table the parallel
  # lane would fail wholesale and silently fall back to sequential (on Macaca
  # scale that forfeits the entire parallel gain). Export the whole
  # environment instead (minus the cluster object), same as stageA.
  parallel::clusterExport(cl, setdiff(ls(), "cl"))

  res <- tryCatch({
    parallel::parLapplyLB(cl, chunks, function(idx) {
      lapply(idx, function(i) {
        fn <- conf$File_Name[i]
        k <- ip[i]
        sb <- data.frame(match_state = "", ul_cov_B = "", ul_by_B = "", ul_all_B = "",
                         lb_fn_B = "", lb_cov_B = "", lb_by_B = "", lb_all_B = "",
                         paired_label_cosine_B = "", choosen_paired_label_cosine_B = "",
                         paired_state_B = "", best_group_B = "",
                         predict_missing_B = "", pair_missing_B = "",
                         vote_metrics_B = "", stringsAsFactors = FALSE)
        status <- "competition but original correct"
        vote_wrong <- 0L; err_reason <- ""; mod_diff <- ""
        chimeric <- 0L; ident_ions <- 0L; tie_flag <- 0L; qual_low <- 0L
        ws <- ""; wm <- ""

        if (!is.na(k)) {
          sb <- data.frame(
            match_state = pr$match_state[k],
            ul_cov_B = pr$unlabel_ions_match.by_ions_coverage[k],
            ul_by_B = pr$unlabel_ions_match.PSM_score_by[k],
            ul_all_B = pr$unlabel_ions_match.PSM_score_all[k],
            lb_fn_B = pr$label_File_names[k],
            lb_cov_B = pr$label_ions_match.by_ions_coverage[k],
            lb_by_B = pr$label_ions_match.PSM_score_by[k],
            lb_all_B = pr$label_ions_match.PSM_score_all[k],
            paired_label_cosine_B = pr$paired_label_cosine[k],
            choosen_paired_label_cosine_B = pr$choosen_paired_label_cosine[k],
            paired_state_B = pr$paired_state[k],
            best_group_B = pr$best_group[k],
            predict_missing_B = pr$predict_missing[k],
            pair_missing_B = pr$pair_missing[k],
            vote_metrics_B = pr$vote_metrics[k],
            stringsAsFactors = FALSE)
          bg <- suppressWarnings(as.integer(strsplit(as.character(pr$best_group[k]),
                                                     ";", fixed = TRUE)[[1]]))
          bg <- bg[!is.na(bg)]
          if (length(bg) == 0) bg <- 1L
          if (length(bg) > 1) tie_flag <- 1L
          sqs <- vapply(slot_seq_cols, function(cn) as.character(conf[i, cn]), character(1))
          mds <- vapply(slot_mod_cols, function(cn) as.character(conf[i, cn]), character(1))
          strat_present <- which(sqs != "")
          cand_sqs <- sqs[strat_present]
          cand_mds <- mds[strat_present]
          win <- bg[1]
          if (win > length(cand_sqs)) win <- 1L
          ws <- cand_sqs[win]
          wm <- cand_mds[win]
          x1s <- if (length(cand_sqs) >= 1) cand_sqs[1] else ""
          x1m <- if (length(cand_mds) >= 1) cand_mds[1] else ""
          vote_wrong <- as.integer(!(ws == x1s && wm == x1m))
          if (!is.null(lu)) {
            chimeric <- chimeric_check(cand_sqs, cand_mds, lu,
                                       as.character(fn), aa, aa_tab)
          } else {
            sr <- if ("candidate_specific_ratio" %in% names(conf)) {
              suppressWarnings(as.numeric(strsplit(
                as.character(conf$candidate_specific_ratio[i]), "|", fixed = TRUE)[[1]]))
            } else numeric(0)
            sr <- sort(sr[!is.na(sr)], decreasing = TRUE)
            if (length(sr) >= 2 && sr[1] > 0.1 && sr[2] > 0.1) chimeric <- 1L
          }
          ul_cov <- suppressWarnings(as.numeric(strsplit(
            as.character(pr$unlabel_ions_match.by_ions_coverage[k]), ";",
            fixed = TRUE)[[1]]))
          if (win <= length(ul_cov) && !is.na(ul_cov[win]) && ul_cov[win] < 20) {
            qual_low <- 1L
          }
          ul_by <- suppressWarnings(as.numeric(strsplit(
            as.character(pr$unlabel_ions_match.PSM_score_by[k]), ";",
            fixed = TRUE)[[1]]))
          others <- setdiff(which(!is.na(ul_by) & seq_along(ul_by) <= length(cand_sqs)), win)
          if (length(others) > 0) {
            runner <- others[which.max(ul_by[others])]
            rs <- cand_sqs[runner]
            rm_ <- cand_mds[runner]
            if (ws != "" && rs != "" && ws != rs) {
              ident_ions <- as.integer(identical_ion_sets(ws, wm, rs, rm_))
            }
          }
          if (chimeric == 1L) {
            status <- "chimeric spectrum"
          } else if (ident_ions == 1L) {
            status <- "ambiguous: identical fragment ions"
          } else if (tie_flag == 1L) {
            status <- "ambiguous: vote tie"
          } else if (qual_low == 1L) {
            status <- "ambiguous: poor spectrum quality"
          } else if (vote_wrong == 1L) {
            status <- "original result wrong"
          }
          if (vote_wrong == 1L) {
            if (nzchar(as.character(conf$Error_reason[i]))) {
              err_reason <- as.character(conf$Error_reason[i])
            }
            if (ws != "" && identical(x1s, ws)) {
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
        } else {
          status <- "not pair-checked (dropped by stage B)"
        }

        data.frame(
          File_Name = fn,
          conf[i, setdiff(colnames(conf), "File_Name"), drop = FALSE],
          sb,
          conclusion_status = status,
          winner_sequence = ws, winner_modification = wm,
          vote_says_wrong = vote_wrong, error_reason = err_reason,
          modification_diff = mod_diff, chimeric = chimeric,
          identical_ions = ident_ions, vote_tie = tie_flag,
          quality_low = qual_low,
          stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
      })
    })
  }, error = function(e) {
    cat("[冲突报告并行] 失败，回退顺序：", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(res)) {
    # 回退：走原顺序循环
    source(sub("stageB_conflict_report_parallel\\.R$", "stageB_conflict_report_sequential\\.R",
               basename(sys.frames()[[1]]$ofile %||% "stageB_conflict_report_sequential.R")),
           local = TRUE, encoding = "UTF-8")
  } else {
    res <- unlist(res, recursive = FALSE)
    cat("[冲突报告并行] 完成，合并 ", length(res), " 行\n", sep = "")
  }
} else {
  # ---- 顺序路径（原代码） ----
  res <- vector("list", nrow(conf))
  for (i in seq_len(nrow(conf))) {
  fn <- conf$File_Name[i]
  k <- ip[i]
  # 阶段 B 新增列（无成对结果时留空）
  sb <- data.frame(match_state = "", ul_cov_B = "", ul_by_B = "", ul_all_B = "",
                   lb_fn_B = "", lb_cov_B = "", lb_by_B = "", lb_all_B = "",
                   paired_label_cosine_B = "", choosen_paired_label_cosine_B = "",
                   paired_state_B = "", best_group_B = "",
                   predict_missing_B = "", pair_missing_B = "",
                   vote_metrics_B = "", stringsAsFactors = FALSE)
  status <- "competition but original correct"
  vote_wrong <- 0L; err_reason <- ""; mod_diff <- ""
  chimeric <- 0L; ident_ions <- 0L; tie_flag <- 0L; qual_low <- 0L
  ws <- ""; wm <- ""

  if (!is.na(k)) {
    sb <- data.frame(
      match_state = pr$match_state[k],
      ul_cov_B = pr$unlabel_ions_match.by_ions_coverage[k],
      ul_by_B = pr$unlabel_ions_match.PSM_score_by[k],
      ul_all_B = pr$unlabel_ions_match.PSM_score_all[k],
      lb_fn_B = pr$label_File_names[k],
      lb_cov_B = pr$label_ions_match.by_ions_coverage[k],
      lb_by_B = pr$label_ions_match.PSM_score_by[k],
      lb_all_B = pr$label_ions_match.PSM_score_all[k],
      paired_label_cosine_B = pr$paired_label_cosine[k],
      choosen_paired_label_cosine_B = pr$choosen_paired_label_cosine[k],
      paired_state_B = pr$paired_state[k],
      best_group_B = pr$best_group[k],
      predict_missing_B = pr$predict_missing[k],
      pair_missing_B = pr$pair_missing[k],
      vote_metrics_B = pr$vote_metrics[k],
      stringsAsFactors = FALSE)
    # 多数投票最佳组（best_group 是候选列表索引：按策略 1..5 顺序跳过空序列；
    # 候选 1＝策略 1 的原始结果，与阶段 B 的候选顺序一致）
    # the majority-vote best group (best_group indexes the candidate list built
    # in strategy order 1..5 skipping empty ones; candidate 1 = the strategy-1
    # original result, matching stage B's candidate order)
    bg <- suppressWarnings(as.integer(strsplit(as.character(pr$best_group[k]),
                                               ";", fixed = TRUE)[[1]]))
    bg <- bg[!is.na(bg)]
    if (length(bg) == 0) bg <- 1L
    if (length(bg) > 1) tie_flag <- 1L
    # 多策略候选（槽列 a.Sequence 等与阶段 B 候选同源；候选列表＝非空策略按序排列）
    sqs <- vapply(slot_seq_cols, function(cn) as.character(conf[i, cn]), character(1))
    mds <- vapply(slot_mod_cols, function(cn) as.character(conf[i, cn]), character(1))
    strat_present <- which(sqs != "")
    cand_sqs <- sqs[strat_present]
    cand_mds <- mds[strat_present]
    win <- bg[1]
    if (win > length(cand_sqs)) win <- 1L
    ws <- cand_sqs[win]
    wm <- cand_mds[win]
    x1s <- if (length(cand_sqs) >= 1) cand_sqs[1] else ""
    x1m <- if (length(cand_mds) >= 1) cand_mds[1] else ""
    # 内容口径判定：胜选候选与原始结果（候选 1）序列＋修饰一致即视为原始结果正确
    #（不同策略组可提出完全相同的候选内容，仅索引不同不算错误）
    # content-based call: when the winner's sequence AND modification equal the
    # original (candidate 1), the original result counts as correct (different
    # strategy groups may propose the same content - a mere index difference is
    # not an error)
    vote_wrong <- as.integer(!(ws == x1s && wm == x1m))
    # 混合谱：不同序列两两"各自特异子离子比例的最小值 > 10%"（用户口径，
    # m/z 口径：特异峰数 ÷ 两候选命中峰并集；2026-09-06 修订）
    # chimeric: any pair of distinct sequences with min of the two
    # specific-ion ratios > 10% (the user's definition; m/z basis: specific
    # peaks over the union of both candidates' matched peaks; revised 2026-09-06)
    if (!is.null(lu)) {
      chimeric <- chimeric_check(cand_sqs, cand_mds, lu,
                                 as.character(fn), aa, aa_tab)
    } else {
      sr <- if ("candidate_specific_ratio" %in% names(conf)) {
        suppressWarnings(as.numeric(strsplit(
          as.character(conf$candidate_specific_ratio[i]), "|", fixed = TRUE)[[1]]))
      } else numeric(0)
      sr <- sort(sr[!is.na(sr)], decreasing = TRUE)
      if (length(sr) >= 2 && sr[1] > 0.1 && sr[2] > 0.1) chimeric <- 1L
    }
    # 谱图质量差：胜选候选无标记侧覆盖度 < 20%（面板 e 校准的判别阈值）
    ul_cov <- suppressWarnings(as.numeric(strsplit(
      as.character(pr$unlabel_ions_match.by_ions_coverage[k]), ";",
      fixed = TRUE)[[1]]))
    if (win <= length(ul_cov) && !is.na(ul_cov[win]) && ul_cov[win] < 20) {
      qual_low <- 1L
    }
    # 子离子完全一致：胜选候选 vs 其余候选中最优者（无标记侧 b/y 打分最高）
    ul_by <- suppressWarnings(as.numeric(strsplit(
      as.character(pr$unlabel_ions_match.PSM_score_by[k]), ";",
      fixed = TRUE)[[1]]))
    others <- setdiff(which(!is.na(ul_by) & seq_along(ul_by) <= length(cand_sqs)), win)
    if (length(others) > 0) {
      runner <- others[which.max(ul_by[others])]
      rs <- cand_sqs[runner]
      rm_ <- cand_mds[runner]
      if (ws != "" && rs != "" && ws != rs) {
        ident_ions <- as.integer(identical_ion_sets(ws, wm, rs, rm_))
      }
    }
    # 结论状态（优先级：混合谱 > 子离子一致 > 投票并列 > 质量差 > 错误 > 正确；英文输出）
    if (chimeric == 1L) {
      status <- "chimeric spectrum"
    } else if (ident_ions == 1L) {
      status <- "ambiguous: identical fragment ions"
    } else if (tie_flag == 1L) {
      status <- "ambiguous: vote tie"
    } else if (qual_low == 1L) {
      status <- "ambiguous: poor spectrum quality"
    } else if (vote_wrong == 1L) {
      status <- "original result wrong"
    }
    # 错误行：按图 4 口径给错误原因与修饰差异细节（原始结果＝候选 1）
    if (vote_wrong == 1L) {
      if (nzchar(as.character(conf$Error_reason[i]))) {
        err_reason <- as.character(conf$Error_reason[i])
      }
      if (ws != "" && identical(x1s, ws)) {
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
  } else {
    status <- "not pair-checked (dropped by stage B)"
  }

  res[[i]] <- data.frame(
    File_Name = fn,
    conf[i, setdiff(colnames(conf), "File_Name"), drop = FALSE],
    sb,
    conclusion_status = status,
    winner_sequence = ws, winner_modification = wm,
    vote_says_wrong = vote_wrong, error_reason = err_reason,
    modification_diff = mod_diff, chimeric = chimeric,
    identical_ions = ident_ions, vote_tie = tie_flag,
    quality_low = qual_low,
    stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  }  # close for loop
}  # close else branch (sequential path)
if (nrow(conf) == 0) {
  # 零冲突行：直接写带结论列的空表，避免 rbind 空列表得到 NULL
  # no conflict rows: write the empty table with the conclusion columns directly
  out <- conf
  for (cn in c("match_state", "ul_cov_B", "ul_by_B", "ul_all_B", "lb_fn_B",
               "lb_cov_B", "lb_by_B", "lb_all_B", "paired_label_cosine_B",
               "choosen_paired_label_cosine_B", "paired_state_B", "best_group_B",
               "predict_missing_B", "pair_missing_B", "vote_metrics_B",
               "conclusion_status", "winner_sequence", "winner_modification",
               "vote_says_wrong", "error_reason",
               "modification_diff", "chimeric", "identical_ions", "vote_tie",
               "quality_low")) {
    out[[cn]] <- character(0)
  }
  out$vote_says_wrong <- integer(0)
  out$chimeric <- integer(0)
  out$identical_ions <- integer(0)
  out$vote_tie <- integer(0)
  out$quality_low <- integer(0)
} else {
  out <- do.call(rbind, res)
}
write.table(out, out_f, sep = "\t", row.names = FALSE, quote = FALSE)
cat("conclusion report written:", out_f, "(", nrow(out), "rows )\n")
cat("=== conclusion status distribution 结论状态分布 ===\n")
print(table(out$conclusion_status))
cat("=== wrong rows by error reason 错误原因分布（原始结果错误判定行） ===\n")
print(table(out$error_reason[out$vote_says_wrong == 1]))

# ---- 分组与肽段层面统计（Single / Paired 分开 ＋ 肽段三 panel）--------------
# ---- grouped and peptide-level statistics (Single vs Paired + 3 panels) -----
# 谱图级按 paired_state_B 分组（Single＝无有效成对标记证据，Paired＝成对余弦
# ≥0.7）；肽段级取每谱最终肽段（投票胜选序列，无阶段 B 结果时退回原始
# X1.Sequence）；三个 panel：Single 谱贡献的肽段、Paired 谱贡献的肽段、合并
# （全部冲突谱的去重肽段）。同一肽段可出现在多个 panel（其谱图分属两组），
# 也可落在多个结论类目（其谱图结论不同），均按各自谱图计入对应类目。
# Spectrum level groups by paired_state_B (Single = no valid paired evidence,
# Paired = paired cosine >= 0.7); peptide level takes each spectrum's final
# peptide (the vote-winner sequence, falling back to X1.Sequence when stage B
# dropped the spectrum); three panels: peptides from Single spectra, from
# Paired spectra, and combined (unique peptides over ALL conflict spectra).
# The same peptide may appear in several panels (its spectra split across
# groups) and in several conclusion categories (its spectra concluded
# differently) - each spectrum counts toward its own category.
if (nrow(out) > 0) {
  pep <- ifelse(nzchar(as.character(out$winner_sequence)),
                as.character(out$winner_sequence),
                as.character(out[["a.Sequence"]]))
  grp <- as.character(out$paired_state_B)
  con <- as.character(out$conclusion_status)
  er  <- as.character(out$error_reason)
  wr  <- suppressWarnings(as.integer(out$vote_says_wrong))

  panel_stats <- function(rows, tag) {
    lines <- c(
      sprintf("\n==== [%s] ====", tag),
      sprintf("spectra rows 谱图行数: %d", sum(rows)),
      sprintf("unique peptides 去重肽段数: %d",
              length(unique(pep[rows]))),
      "conclusion status (spectra) 结论状态分布（谱图级）:",
      paste(capture.output(print(table(con[rows]))), collapse = "\n"),
      "unique peptides per conclusion 每类结论的去重肽段数:",
      paste(capture.output(print(sort(tapply(pep[rows], con[rows],
                                       function(x) length(unique(x))),
                                decreasing = TRUE))), collapse = "\n"),
      "error reasons of wrong rows 错误原因分布（原始结果错误行）:",
      paste(capture.output(print(table(er[rows & wr == 1]))), collapse = "\n")
    )
    lines
  }
  lines <- c(
    "================ 冲突谱分组与肽段层面统计 ================",
    "================ conflict spectra: grouped + peptide-level stats ================",
    panel_stats(grp == "single", "Single（无有效成对证据 no valid paired evidence）"),
    panel_stats(grp == "paired", "Paired（有成对标记证据 paired-label evidence）"),
    panel_stats(rep(TRUE, nrow(out)), "Combined 合并（全部冲突谱 all conflict spectra）"),
    sprintf("\nnot pair-checked rows 未核查谱数: %d", sum(!nzchar(grp)))
  )
  cat(paste(lines, collapse = "\n"), "\n")
  sum_f <- sub("\\.tsv$", "_summary.txt", out_f)
  writeLines(enc2utf8(paste(lines, collapse = "\n")), sum_f, useBytes = TRUE)
  cat("summary written:", sum_f, "\n")
}
