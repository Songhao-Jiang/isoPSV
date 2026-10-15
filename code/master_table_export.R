# =============================================================================
# master_table_export.R - 全谱主表导出（与 8.衣藻全谱主表.txt 同格式，全谱行）
# master_table_export.R - full-spectrum master-table export (same format as the
# 8.衣藻全谱主表.txt master table, one row per spectrum)
#
# 用途：
# 生成与 "8.衣藻全谱主表.txt"（94 列）同格式的全谱主表，**每一张已鉴定谱一行**：
#   - 冲突谱（多引擎不一致，进入阶段 B 的行）：候选 1–5 全部参数、判定结论、
#     蛋白定位、差异类型全部填值；成对标记相关列按阶段 B 结果填写；
#   - 无差异谱（多引擎一致，未进入阶段 B 的行）：候选按阶段 A 报告列出，
#     结论列留空（无需校正）、corrected_* 复制候选 1；标记相关列留空；
#   - 无标记样本：全表同格式输出，标记相关列全部留空。
#
# Purpose:
# Generates the same 94-column format as "8.衣藻全谱主表.txt", **one row for
# every identified spectrum**:
#   - conflict spectra (five strategies disagree, i.e. the stage-B rows): all
#     candidate-1..5 parameters, the judgment conclusion, the protein mapping
#     and the difference type are filled; the paired-label columns come from
#     stage B;
#   - non-conflict spectra (five strategies agree, never enter stage B): the
#     candidates are listed from the stage-A report, the conclusion columns
#     stay blank (no correction needed) and corrected_* copies candidate 1;
#     the labeled columns stay blank;
#   - unlabeled samples: the same format with every labeled column blank.
#
# 输入（环境变量，由 user_run.R 设置）：
#   CONFIG3 ＝ 第三步参数文件路径（定位结果目录、标记类型与搜库 FASTA）
#   ALL_MGF ＝ 原始 MGF 文件夹（仅嵌合特异比例需要；可为空）
#
# Inputs (environment variables, set by user_run.R):
#   CONFIG3 = path of the step-3 parameter file (locates the result folder,
#             the label type and the search FASTA)
#   ALL_MGF = the original MGF folder (only for the chimeric specific-ratio
#             columns; may be empty)
#
# 输出：
#   <结果目录>/master_table.tsv（94 列；全谱一行；无差异谱结论列留空）
#
# Outputs:
#   <result folder>/master_table.tsv (94 columns; one row per spectrum; the
#   conclusion columns stay blank for non-conflict spectra)
# =============================================================================
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
# 同 user_run.R：normalizePath 在中文路径＋R 4.5 C.UTF-8 下会产出腐蚀路径，
# 校验不通过时按字符串解析回退
# same fallback as user_run.R: normalizePath corrupts CJK paths under R 4.5 +
# C.UTF-8; resolve by string when the candidate dir does not exist
if (!dir.exists(script_dir)) {
  script_dir <- dirname(if (grepl("^[A-Za-z]:", script_arg)) script_arg else file.path(getwd(), script_arg))
  script_dir <- gsub("\\\\", "/", script_dir)
}
# ---- 库环境 ----
# ---- library environment
old_wd <- setwd(script_dir)
source("compat.R")
source("ReadMGF_script.R")
suppressMessages(source("acs_core.R")); suppressMessages(source("paired_params.R"))
source("chimeric_check.R")
setwd(old_wd)
aa <- read.delim(file.path(script_dir, "AA_mass_dataframe.txt"), sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# ---- 环境与参数 ----
# ---- environment and parameters
cfg3_path <- Sys.getenv("CONFIG3")
if (!nzchar(cfg3_path)) stop("CONFIG3 not set (set by user_run.R)")
cfg3 <- read_param_file(cfg3_path)
cfg_dir <- cfg3$config_dir
if (is.null(cfg_dir) || length(cfg_dir) == 0 || !nzchar(cfg_dir[1]))   # 审查 P0-3/P3-6：config_dir 缺失/为空时给可读报错
  stop("config_step3 的 config_dir 缺失或为空，请检查配置（或删除 result_* 后重跑）
config_dir missing/empty in config_step3")
# 结果文件夹（RESULT_DIR 由总程序注入；独立运行按时间戳定位/新建）
# result folder (RESULT_DIR injected by the master; standalone runs locate/create by timestamp)
result_dir <- get_result_dir(cfg_dir)
base_cfg <- read_param_file(resolve_cfg_paths(cfg3, "base_config")[1])
label_type <- if (!is.null(base_cfg$label_type) && length(base_cfg$label_type) > 0)
  tolower(tail(base_cfg$label_type, 1)) else "unlabeled"
all_mgf_dir <- Sys.getenv("ALL_MGF")
fasta_file <- if (!is.null(base_cfg$fasta) && length(base_cfg$fasta) > 0)
  tail(resolve_cfg_paths(base_cfg, "fasta"), 1) else ""

# P0-3无标记样本检查前置——旧位置在阶段 A 报告读入之后，用户
# 跑完大半解析才失败；改为解析完配置立即判断（成对报告缺失＝无标记样本）。
# 无标记样本暂不支持全谱主表本表依赖成对报告的
# unlabel_ions_match.* 等列名，而阶段 A 报告是 candidate_* 列名体系，列映射缺失；
# 无标记样本的科学需求由 stageA/report_A.tsv 与 classification_A.tsv 覆盖。
# 直接调用本脚本且无成对报告时，以明确提示退出而不是产出缺列的表。
# P0-3 the unlabeled check moved to the front - it used to sit
# after the stage-A report read, so the user failed only after most parsing
# was done; now it runs right after config resolution (no paired report =
# unlabeled sample).
stb_dir <- file.path(result_dir, "2.labeled_result")
pr_path <- file.path(stb_dir, paste0("paired_report_", label_type, ".tsv"))
labeled_run <- file.exists(pr_path)
if (!labeled_run) {
  stop("[主表导出] 无标记样本暂不支持全谱主表（未完成功能）；",
       "请使用 stageA/report_A.tsv 与 stageA/classification_A.tsv。",
       "\n[master table] the full-spectrum master table is not yet available for ",
       "unlabeled samples (unfinished feature); use stageA/report_A.tsv and ",
       "stageA/classification_A.tsv instead.")
}

# ---- 定稿判定模型参数 ----
# ---- finalized judgment-model parameters
W6 <- c(ul_by = 0.6200, lb_by = 0.0340, ul_all = 0.1186, lb_all = 0.1700,
        predict = 0, pcos = 0.0573)
TAU_MAP <- c(n15 = 0.0083, c13 = 0.0083, silac = 0.0183, unlabeled = 0.0083)
THETA <- 0.30
tau_v <- label_tau_of(parse_label_spec(label_type))
if (is.null(tau_v)) tau_v <- 0.0083
# ---- 分组权重覆盖--------------------------
# config_step3 可选键 weights_B（六指标逗号分隔，顺序 ul_by,lb_by,ul_all,
# lb_all,predict,pcos）与 threshold_B（本标签类型的 τ）优先于上方统一默认值；
# 旧配置无这两键时行为与原版完全一致（统一 W\*）
# ---- per-label weight override -----------
# optional config keys weights_B (six comma-separated metrics in the order
# above) and threshold_B (the τ of THIS label type) take precedence; configs
# without them reproduce the unified W\* behavior exactly
wb_raw <- if (!is.null(cfg3$weights_B) && length(cfg3$weights_B) > 0)
  tail(cfg3$weights_B, 1) else ""
if (nzchar(wb_raw)) {
  wb <- suppressWarnings(as.numeric(strsplit(wb_raw, ",")[[1]]))
  if (!length(wb) %in% c(6L, 8L) || any(is.na(wb)))
    stop(sprintf("weights_B 应为 6 或 8 个逗号分隔数值（六指标或六指标＋ul_cov,lb_cov），当前：%s
weights_B must be 6 or 8 comma-separated numbers, got: %s",
        wb_raw, wb_raw))
  if (length(wb) == 6L) {
    W6 <- setNames(wb, names(W6))
  } else {
    W6 <- setNames(wb, c("ul_by","lb_by","ul_all","lb_all",
                         "predict","pcos","ul_cov","lb_cov"))
  }
  cat("[分组权重] weights_B 覆盖：", paste(sprintf("%.4f", W6), collapse = ","), "\n")
  cat("[per-label] weights_B override: ", paste(sprintf("%.4f", W6), collapse = ","), "\n")
}
tb_raw <- if (!is.null(cfg3$threshold_B) && length(cfg3$threshold_B) > 0)
  tail(cfg3$threshold_B, 1) else ""
if (nzchar(tb_raw)) {
  tau_v <- as.numeric(tb_raw)
  cat("[分组权重] threshold_B 覆盖：", tau_v, "\n")
  cat("[per-label] threshold_B override: ", tau_v, "\n")
}

# ---- 谱图全集：阶段 A 报告（每张已鉴定谱一行，含候选与无标记指标）----
# ---- the spectrum universe: the stage-A report (one row per identified
# spectrum with the candidates and the unlabeled metrics)
p <- read.delim(file.path(result_dir, "1.unlabeled_result", "report_A.tsv"), sep = "\t",
                check.names = FALSE, na.strings = character(0), quote = "")
if (!"candidates" %in% names(p) && any(grepl("^[a-z]\\.Sequence$", names(p)))) {
  # 多引擎版：以槽列 a.Sequence 等重建候选串（阶段 A 报告无 candidates 列时的兜底）
  # multi-engine edition: rebuild the candidates string from the slot columns
  # a.Sequence etc. (fallback when the stage-A report has no candidates column)
  slot_seq_cols <- grep("^[a-z]\\.Sequence$", names(p), value = TRUE)
  slot_mod_cols <- grep("^[a-z]\\.Modification$", names(p), value = TRUE)
  p$candidates <- vapply(seq_len(nrow(p)), function(i) {
    parts <- character(0)
    for (j in seq_along(slot_seq_cols)) {
      sq <- as.character(p[[slot_seq_cols[j]]][i])
      md <- as.character(p[[slot_mod_cols[j]]][i])
      if (!nzchar(sq)) next
      parts <- c(parts, if (nzchar(md)) paste0(sq, "(", md, ")") else sq)
    }
    paste(parts, collapse = "|")
  }, character(1))
}
# 多引擎版：候选槽数（NSLOT）由报告槽列动态推断，主表候选列相应扩展
# multi-engine edition: the candidate slot count (NSLOT) is inferred from the
# report slot columns and the master-table candidate columns extend accordingly
NSLOT <- length(grep("^[a-z]\\.Sequence$", names(p), value = TRUE))
if (NSLOT == 0) NSLOT <- 1L
n <- nrow(p)

# ---- 冲突行定位：成对报告（阶段 B）与冲突报告按 File_Name 对齐 ----
# ---- conflict-row mapping: align the paired report (stage B) and the conflict
# report by File_Name onto the stage-A rows
pr_i <- rep(NA_integer_, n)                                    # 阶段A行 → 成对报告行
pr <- NULL; conf <- NULL
if (labeled_run) {
  pr <- read.delim(pr_path, sep = "\t", check.names = FALSE,
                   na.strings = character(0), quote = "")
  pr_i <- match(as.character(p$File_Name), as.character(pr$File_Name))
  cr_path <- file.path(stb_dir, paste0("conflict_report_", label_type, ".tsv"))
  if (file.exists(cr_path)) {
    conf <- read.delim(cr_path, sep = "\t", check.names = FALSE,
                       na.strings = character(0), quote = "")
    conf <- conf[match(as.character(p$File_Name), as.character(conf$File_Name)), ]
  }
}
conflict_row <- !is.na(pr_i)                                   # 冲突行（进入阶段 B）
cval <- function(df, cn, dflt = "-") {                         # 缺失列/行安全取值
  if (is.null(df) || !cn %in% names(df)) return(rep(dflt, n))
  x <- as.character(df[[cn]])
  x[is.na(x)] <- dflt
  x
}

# ---- 指标解析（无标记指标来自阶段 A 报告；标记指标来自成对报告）----
# ---- metric parsing (unlabeled metrics from the stage-A report; labeled
# metrics from the paired report)
parse_col <- function(x) lapply(strsplit(as.character(x), "[|;]"), function(v)
  suppressWarnings(as.numeric(v)))
ncand <- lengths(strsplit(as.character(p$candidates), "|", fixed = TRUE))
# 冲突行的候选数以成对报告 match_state 为准（阶段 A 的 candidates 列多数行只有
# 1 条，直接用它会让主表丢失竞争结构、漏判需校正——修复，与谱图
# PDF 判定同口径）
# the candidate count of conflict rows follows the paired report's match_state
# (the stage-A candidates column holds a single entry for most rows; using it
# directly drops the competition structure and under-judges corrections -
# fixed the same basis as the spectrum-PDF judgment)
if (!is.null(pr) && "match_state" %in% names(pr)) {
  pr_nc <- lengths(strsplit(as.character(pr$match_state), ";", fixed = TRUE))
  for (i in which(conflict_row)) ncand[i] <- pr_nc[pr_i[i]]
}
# 无标记指标：冲突行优先用成对报告（与 panel_f 主表同源，阶段 B 在原始 MGF 上
# 重算）；无差异行用阶段 A 报告（分隔符 | 或; 均兼容）
# unlabeled metrics: conflict rows prefer the paired report (the same source
# as the panel-f master table; stage B recomputes them on the raw MGF), while
# non-conflict rows use the stage-A report (both | and; separators accepted)
ul_col <- list(ul_by = "candidate_PSM_score_by", ul_all = "candidate_PSM_score_all",
               predict = "candidate_predict_cosine")
pr_ul <- list(ul_by = "unlabel_ions_match.PSM_score_by",
              ul_all = "unlabel_ions_match.PSM_score_all",
              predict = "predict_cosines")
LUL <- lapply(names(ul_col), function(mn) {
  cn <- ul_col[[mn]]
  if (!cn %in% names(p)) return(lapply(seq_len(n), function(i) rep(NA_real_, ncand[i])))
  L <- parse_col(p[[cn]])
  out <- Map(function(v, k) {
    if (length(v) >= k) v[seq_len(k)] else c(v, rep(NA_real_, k - length(v)))
  }, L, ncand)
  if (!is.null(pr) && pr_ul[[mn]] %in% names(pr)) {           # 冲突行换用成对报告
    LP <- parse_col(pr[[pr_ul[[mn]]]])
    for (i in which(conflict_row)) {
      v <- LP[[pr_i[i]]]
      out[[i]] <- if (length(v) >= ncand[i]) v[seq_len(ncand[i])] else
        c(v, rep(NA_real_, ncand[i] - length(v)))
    }
  }
  out
})
names(LUL) <- names(ul_col)
# 覆盖度：冲突行用成对报告，无差异行用阶段 A 报告。
# 注意：整列只解析一次（逐行循环内重复 parse_col 整列是 75k×75k 量级，极慢）
# coverage: the paired report for conflict rows, the stage-A report otherwise.
# Note: each whole column is parsed only once (re-parsing a whole column inside
# the per-row loop is 75k x 75k scale and extremely slow)
pr_cov_parsed <- if (!is.null(pr) &&
                     "unlabel_ions_match.by_ions_coverage" %in% names(pr))
  parse_col(pr[["unlabel_ions_match.by_ions_coverage"]]) else NULL
p_cov_parsed <- if ("candidate_coverage" %in% names(p))
  parse_col(p[["candidate_coverage"]]) else NULL
covL <- lapply(seq_len(n), function(i) {
  v <- if (conflict_row[i] && !is.null(pr_cov_parsed) && length(pr_cov_parsed) >= pr_i[i] &&
           !is.na(pr_i[i])) {
    pr_cov_parsed[[pr_i[i]]]
  } else if (!is.null(p_cov_parsed) && length(p_cov_parsed) >= i) {
    p_cov_parsed[[i]]
  } else rep(NA_real_, ncand[i])   # 缺列/越界一律 NA（原写法 else 分支返回值被丢弃）
  if (length(v) >= ncand[i]) v[seq_len(ncand[i])] else
    c(v, rep(NA_real_, ncand[i] - length(v)))
})
# 标记指标：仅冲突行有值（成对报告），其余行整谱 NA（权重自动降级）
# labeled metrics: only conflict rows carry values (from the paired report);
# other rows are all-NA so the weights degrade automatically
lb_col <- list(lb_by = "label_ions_match.PSM_score_by",
               lb_all = "label_ions_match.PSM_score_all",
               pcos = "paired_label_cosine",
               lb_cov = "label_ions_match.by_ions_coverage")
LLB <- lapply(lb_col, function(cn) {
  out <- lapply(seq_len(n), function(i) rep(NA_real_, ncand[i]))
  if (is.null(pr) || !cn %in% names(pr)) return(out)
  L <- parse_col(pr[[cn]])
  for (i in which(conflict_row)) {
    v <- L[[pr_i[i]]]
    out[[i]] <- if (length(v) >= ncand[i]) v[seq_len(ncand[i])] else
      c(v, rep(NA_real_, ncand[i] - length(v)))
  }
  out
})
# 指标总表（按谱）——weights_B 给 8 值时纳入双覆盖率；covL 即无标记 by 覆盖率（上文已解析），lb_cov 由 LLB 新增
# the combined metric table (per spectrum) - with an 8-value weights_B the two
# coverage metrics join covL is the unlabeled by
# coverage parsed above; lb_cov comes from the extended LLB
LL <- list(ul_by = LUL[["ul_by"]], lb_by = LLB[["lb_by"]],
           ul_all = LUL[["ul_all"]], lb_all = LLB[["lb_all"]],
           predict = LUL[["predict"]], pcos = LLB[["pcos"]])
if ("ul_cov" %in% names(W6)) LL$ul_cov <- covL
if ("lb_cov" %in% names(W6)) LL$lb_cov <- LLB[["lb_cov"]]
DIV <- c(ul_by = 100, lb_by = 100, ul_all = 100, lb_all = 100, predict = 1,
         pcos = 1, ul_cov = 100, lb_cov = 100)
# 配对状态与标记谱号：冲突行取成对报告，其余行留空
# paired state and labeled scan: conflict rows from the paired report, blank
# elsewhere
#：阶段 B 现覆盖全部谱，成对状态一律取成对报告（"matched" 占位
# 仅保留给报告里确实没有的行，如 raw 缺失）
mstL <- lapply(seq_len(n), function(i) {
  if (!is.na(pr_i[i]))
    strsplit(as.character(pr$match_state[pr_i[i]]), ";", fixed = TRUE)[[1]] else
      rep("matched", ncand[i])
})
labL <- lapply(seq_len(n), function(i) {
  if (!is.na(pr_i[i]))
    strsplit(as.character(pr$label_File_names[pr_i[i]]), ";", fixed = TRUE)[[1]] else
      character(0)
})

# ---- 逐谱打分与判定（仅冲突行；无差异行结论留空、保持原始）----
# ---- per-spectrum scoring and judgment (conflict rows only; non-conflict rows
# keep the original with blank conclusion columns)
score_of <- function(i) {
  M <- sapply(names(W6), function(mn) LL[[mn]][[i]] / DIV[[mn]])
  # 单候选谱 sapply 返回向量，需强制为矩阵再按列判断可用性
  # single-candidate rows collapse sapply to a vector; force a matrix before
  # the per-column availability check
  if (is.null(dim(M))) M <- matrix(M, nrow = 1L)
  avail <- apply(M, 2, function(colv) all(!is.na(colv)))
  ww <- W6[avail]
  if (sum(ww) <= 0) ww <- rep(1, sum(avail))
  ww <- ww / sum(ww)
  rowSums(sweep(M[, avail, drop = FALSE], 2, ww, "*"))
}
scores <- vector("list", n)
rev_idx <- integer(n)
delta <- rep(NA_real_, n)
for (i in which(conflict_row)) {
  s <- score_of(i)
  scores[[i]] <- s
  if (length(s) >= 2) {
    k2 <- which.max(s[-1]) + 1L
    rev_idx[i] <- k2
    delta[i] <- s[k2] - s[1]
  }
}
dpred <- vapply(seq_len(n), function(i) {
  prm <- LL[["predict"]][[i]]
  if (conflict_row[i] && length(prm) >= 2 && rev_idx[i] >= 2 &&
      !is.na(prm[1]) && !is.na(prm[rev_idx[i]]))
    prm[rev_idx[i]] - prm[1] else NA_real_
}, numeric(1))
judged <- rep(FALSE, n)
rescue <- rep(FALSE, n)
for (i in which(conflict_row)) {
  if (!is.na(delta[i]) && delta[i] >= tau_v) judged[i] <- TRUE
  if (!is.na(delta[i]) && delta[i] < tau_v && !is.na(dpred[i]) && dpred[i] >= THETA) {
    judged[i] <- TRUE; rescue[i] <- TRUE
  }
}
ions_identical <- !is.na(delta) & delta == 0
ions_small_diff <- !is.na(delta) & delta != 0 & abs(tau_v - delta) <= 0.015
cov_below20 <- vapply(seq_len(n), function(i) {
  v <- covL[[i]]
  !is.na(v[1]) && v[1] < 20
}, logical(1))


# ---- 蛋白与前后氨基酸定位（I/L 等价去重＋搜库 FASTA 直接定位＋互换前继承）----
# ---- protein and flanking-residue mapping (I/L-equivalence dedup + direct
# search in the search FASTA + pre-swap inheritance)
fa_env <- new.env(parent = emptyenv())
if (nzchar(fasta_file) && file.exists(fasta_file)) {
  ln <- readLines(fasta_file, warn = FALSE)
  hd <- grep("^>", ln)
  fa_acc <- trimws(gsub("\r", "", sub("^>([^ \t]+).*", "\\1", ln[hd])))
  ends <- c(hd[-1] - 1L, length(ln))
  fa_seq <- gsub("\r", "", vapply(seq_along(hd), function(k)
    toupper(paste(ln[(hd[k] + 1L):ends[k]], collapse = "")), character(1)))
  keep_t <- !startsWith(fa_acc, "REV_")
  fa_acc <- fa_acc[keep_t]; fa_seq <- fa_seq[keep_t]
  big_il <- chartr("LI", "II", paste(fa_seq, collapse = "\n"))
  cum_start <- cumsum(c(1L, nchar(fa_seq) + 1L))
  fmap <- function(sq) {                                       # c(蛋白, 定位) 或 NULL
    q <- chartr("LI", "II", toupper(sq))
    if (!nzchar(q)) return(NULL)
    m <- regexpr(q, big_il, fixed = TRUE)
    if (m[1] < 1L) return(NULL)
    pi <- findInterval(m[1], cum_start)
    st <- m[1] - cum_start[pi] + 1L
    nq <- nchar(q); ns <- nchar(fa_seq[pi])
    pre <- if (st > 1L) substr(fa_seq[pi], st - 1L, st - 1L) else "_"
    post <- if (st + nq <= ns) substr(fa_seq[pi], st + nq, st + nq) else "_"
    # pFind 原生 Positions 格式"位点,前一残基,后一残基/"
    # 位点＝0 起前置残基数，缺相邻残基用 "_"
    c(fa_acc[pi], sprintf("%d,%s,%s/", st - 1L, pre, post))
  }
  assign("fmap", fmap, envir = fa_env)
}
loc_cache <- new.env(parent = emptyenv())                    # 肽段级缓存（I/L 键）
# B-7环境不能存 NULL——assign(key, NULL) 在 R 里是"删除绑定"，
# 漏命中缓存形同虚设（v1.0 的修复无效）。改用 NA_character_ 哨兵：get 后
# is.na 判定即视为未命中。
loc_of <- function(sq) {                                       # 直接＋GG→N 回退
  if (!nzchar(sq) || !exists("fmap", envir = fa_env)) return(NULL)
  key <- chartr("LI", "II", sq)                               # 唯一肽段只定位一次
  if (exists(key, envir = loc_cache, inherits = FALSE)) {
    hit <- get(key, envir = loc_cache, inherits = FALSE)
    return(if (length(hit) == 1L && is.na(hit)) NULL else hit)
  }
  fmap <- get("fmap", envir = fa_env)
  v <- fmap(sq)
  if (!is.null(v)) { assign(key, v, envir = loc_cache); return(v) }
  gg <- gregexpr("GG", toupper(sq), fixed = TRUE)[[1]]; gg <- gg[gg > 0]
  if (length(gg) > 0 && length(gg) <= 5) {
    for (mask in seq_len(2^length(gg)) - 1L) {
      b <- as.integer(intToBits(mask)[seq_along(gg)])
      vv <- sq; off <- 0L
      for (j in seq_along(gg)) if (b[j] == 1L) {
        pp <- gg[j] + off
        vv <- paste0(substr(vv, 1L, pp - 1L), "N", substr(vv, pp + 2L, nchar(vv)))
        off <- off - 1L
      }
      w <- fmap(vv)
      if (!is.null(w)) { assign(key, w, envir = loc_cache); return(w) }
    }
  }
  # 直接与回退都未命中时也缓存 NULL，避免重复扫描
  # （审查 P2-3：原 assign 误落在 if 块内，gg=0 或 >5 的肽段每条
  # 漏命中都重扫全库 5.55ms；移到 if 外并缓存全路径）
  # cache the NULL result too so the fallback is not re-scanned
  # (review P2-3: the assign used to sit INSIDE the if block, so
  # peptides with 0 or >5 GG re-scanned the whole DB at 5.55 ms per miss;
  # moved out of the if so every path is cached)
  assign(key, NA_character_, envir = loc_cache)   # 哨兵＝未命中
  NULL
}
# 互换前继承：校正序列直接未命中但与 cand1 为衍生变体时报告 cand1 位点
# pre-swap inheritance: a corrected sequence that misses directly but is a
# derived variant of cand1 takes cand1's location
inherit_loc <- function(a, b, loc1) {
  if (!nzchar(a) || !nzchar(b) || !nzchar(loc1)) return(NULL)
  if (adist(chartr("LI", "II", a), chartr("LI", "II", b))[1, 1] > nchar(b) %/% 3) return(NULL)
  # pFind 原生 Positions 格式 "位点,前一残基,后一残基/"
  f <- regmatches(loc1, regexec("^([0-9]+),([A-Z_]),([A-Z_])/?$", loc1))[[1]]
  if (length(f) == 0L) return(NULL)
  pos0 <- as.integer(f[2]); pre0 <- f[3]; post0 <- f[4]
  if (b == substr(a, 2L, nchar(a)))
    return(sprintf("%d,%s,%s/", pos0 + 1L, substr(a, 1L, 1L), post0))
  if (b == substr(a, 1L, nchar(a) - 1L))
    return(sprintf("%d,%s,%s/", pos0, pre0, substr(a, nchar(a), nchar(a))))
  loc1
}

# ---- diff_type 分类树（与 5.全谱主表v2.R 同口径）----
# ---- the diff_type classification tree (same definition as the panel-f
# master-table script)
AA_MASS <- c(G = 57.02146, A = 71.03711, S = 87.03203, P = 97.05276,
             V = 99.06841, T = 101.04768, C = 103.00919, L = 113.08406,
             I = 113.08406, N = 114.04293, D = 115.02694, Q = 128.05858,
             K = 128.09496, E = 129.04259, M = 131.04049, H = 137.05891,
             F = 147.06841, R = 156.10111, Y = 163.06333, W = 186.07931)
AA3_TAB <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q",
             Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K",
             Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
             Tyr = "Y", Val = "V", Xle = "L")
pep_mass <- function(s) {
  a1 <- unlist(strsplit(toupper(s), ""))
  if (length(a1) == 0 || any(!a1 %in% names(AA_MASS))) NA_real_ else sum(AA_MASS[a1])
}
# 含修饰的完整分子 [M+H]+（残基和＋修饰＋H2O＋H；未知修饰返回 NA）
# [M+H]+ of the intact molecule with modifications (residue sum + mods + H2O +
# H; unknown modifications return NA)
pep_mh_mod <- function(s, mstr) {
  mass <- pep_mass(s)
  if (is.na(mass)) return(NA_real_)
  a1 <- unlist(strsplit(toupper(s), ""))
  for (tok in mod_tokens(mstr)) {
    parts <- strsplit(tok, ",", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    nm <- parts[2]; pos <- suppressWarnings(as.integer(parts[1]))
    if (grepl("->", nm, fixed = TRUE)) {                  # X->Y 替换
      ab <- strsplit(nm, "->", fixed = TRUE)[[1]]
      f1 <- sub("\\[.*$", "", ab[1]); t1r <- sub("\\[.*$", "", ab[2])
      if (t1r == "pyro-Glu") {
        # N 端谷氨酰胺环化：残基质量替换为 pyro-Glu（Q - NH3 = 111.03203）
        if (is.na(pos) || pos == 0) {
          if (length(a1) >= 1) mass <- mass + 111.03203 - AA_MASS[[a1[1]]]
        } else if (pos >= 1 && pos <= length(a1)) {
          mass <- mass + 111.03203 - AA_MASS[[a1[pos]]]
        } else return(NA_real_)
      } else {
        f1x <- if (f1 %in% names(AA3_TAB)) AA3_TAB[[f1]] else f1
        t1 <- if (t1r %in% names(AA3_TAB)) AA3_TAB[[t1r]] else NA_character_
        if (is.na(t1) || !f1x %in% names(AA_MASS)) return(NA_real_)
        if (is.na(pos) || pos == 0) {
          if (length(a1) >= 1) mass <- mass + AA_MASS[[t1]] - AA_MASS[[a1[1]]]
        } else if (pos >= 1 && pos <= length(a1)) {
          mass <- mass + AA_MASS[[t1]] - AA_MASS[[a1[pos]]]
        } else return(NA_real_)
      }
    } else if (grepl("Oxidation", nm, fixed = TRUE)) mass <- mass + 15.994915
    else if (grepl("Acetyl", nm, fixed = TRUE)) mass <- mass + 42.010565
    else if (grepl("Carbamidomethyl", nm, fixed = TRUE)) mass <- mass + 57.021464
    else if (grepl("Deamidat", nm, fixed = TRUE)) mass <- mass + 0.984016
    else if (grepl("pyro-Glu", nm, fixed = TRUE)) mass <- mass - 17.026549
    else return(NA_real_)                                 # 未知修饰：放弃判定
  }
  mass + 18.010565 + 1.007276
}
mod_tokens <- function(mstr) {
  if (is.na(mstr) || mstr == "") return(character(0))
  t <- trimws(strsplit(mstr, ";", fixed = TRUE)[[1]])
  t[t != ""]
}
mod_type <- function(tok) sub("^[0-9]+,", "", tok)
sub0 <- function(tb, nm) as.integer(if (nm %in% names(tb)) tb[[nm]] else 0)
mod_cancel <- function(t1, t2) {
  tab1 <- table(mod_type(t1)); tab2 <- table(mod_type(t2))
  keep1 <- rep(names(tab1), pmax(0, tab1 - pmin(tab1, vapply(names(tab1),
             function(nm) sub0(tab2, nm), integer(1)))))
  keep2 <- rep(names(tab2), pmax(0, tab2 - pmin(tab2, vapply(names(tab2),
             function(nm) sub0(tab1, nm), integer(1)))))
  list(l = sort(unname(keep1)), r = sort(unname(keep2)))
}
aa_cancel <- function(s1, s2) {
  a1 <- unlist(strsplit(toupper(s1), "")); a2 <- unlist(strsplit(toupper(s2), ""))
  if (length(a1) == 0 || length(a2) == 0) return(list(l = character(0), r = character(0)))
  t1 <- table(a1); t2 <- table(a2)
  k1 <- rep(names(t1), pmax(0, t1 - pmin(t1, vapply(names(t1),
           function(nm) sub0(t2, nm), integer(1)))))
  k2 <- rep(names(t2), pmax(0, t2 - pmin(t2, vapply(names(t2),
           function(nm) sub0(t1, nm), integer(1)))))
  list(l = sort(unname(k1)), r = sort(unname(k2)))
}

# ---- 嵌合特异比例（ALL_MGF 可用且行标记为 chimeric 时计算）----
# ---- chimeric specific-ratio (computed when ALL_MGF is available and the row
# is flagged chimeric)
chim_flag <- cval(conf, "chimeric") == "1"
lu <- NULL
# ---- 异价嵌合判据----
# 被替换肽按物理质量（含修饰）在隔离窗（≤1.0 m/z）内与谱图母离子吻合的价态 z*
# 若不同于谱图标题价态 → "chimera_different_charge"；否则维持原嵌合标签。
# Different-charge chimera test the charge z* whose
# m/z (physical mass with mods) falls within 1.0 m/z of the spectrum precursor
# m/z differs from the spectrum title's charge -> "chimera_different_charge";
# otherwise the row keeps the plain chimeric label.
mer_mh <- NULL
mg_path <- file.path(result_dir, "merge", "1.combine_result.spectra")
if (file.exists(mg_path)) {
  mer <- read.delim(mg_path, sep = "\t", check.names = FALSE,
                    na.strings = character(0), quote = "")
  if ("File_Name" %in% names(mer) && "Exp.MH." %in% names(mer)) {
    mer_mh <- setNames(as.numeric(as.character(mer[["Exp.MH."]])),
                       as.character(mer[["File_Name"]]))
  }
}
diff_charge <- function(fn, seq, mstr) {
  if (is.null(mer_mh)) return(FALSE)
  tt <- strsplit(fn, ".", fixed = TRUE)[[1]]
  if (length(tt) < 4) return(FALSE)
  z_s <- suppressWarnings(as.integer(tt[4]))
  if (is.na(z_s)) return(FALSE)
  em <- unname(mer_mh[fn])
  if (is.na(em)) return(FALSE)
  mz_obs <- (em + (z_s - 1) * 1.007276) / z_s
  mh <- pep_mh_mod(seq, mstr)
  if (is.na(mh)) return(FALSE)
  mzs <- (mh + (0:5) * 1.007276) / (1:6)
  d <- abs(mzs - mz_obs)
  which.min(d) != z_s && min(d) <= 1.0
}
if (nzchar(all_mgf_dir) && any(chim_flag)) {
  # 仅加载嵌合谱所在 raw
  # load only the raws holding chimeric spectra (P7,: loading the
  # whole folder measured 4-5 minutes)
  needed_raws <- unique(unname(vapply(as.character(p$File_Name[chim_flag]),
                                      function(t) title_raw_scan(t)["raw"], character(1))))
  old3 <- setwd(script_dir)
  lu <- tryCatch(load_unlabeled_spectra(all_mgf_dir, only = needed_raws), error = function(e) {
    # 加载失败不再静默：特异比例列将整体留空，必须让用户知情
    # a failed load is no longer silent: the specific-ratio columns stay blank,
    # which the user must be told about 
    cat("[主表导出] 警告：无标记 MGF 加载失败，特异比例列留空：", conditionMessage(e), "
")
    cat("[master table] WARNING: unlabeled MGF load failed, specific-ratio columns left blank: ",
        conditionMessage(e), "
")
    NULL
  })
  setwd(old3)
}
# matched_positions_of 结果缓存主表逐行循环中同一
# (谱图, 序列, 修饰) 会被反复匹配（cand1 与多个候选、跨行重复候选），
# 缓存按键命中后直接返回整数向量，把特异比例计算的重复谱图匹配消除。
# matched_positions_of cache the row loop re-matches the
# same (spectrum, sequence, modification) over and over (cand1 against each
# candidate, repeated candidates across rows); memoize by key and return the
# cached integer vector, removing the duplicate spectrum matches behind the
# specific-ratio columns.
mp_cache <- new.env(parent = emptyenv())
mp_cached <- function(fn, sq, md, sp) {
  key <- paste(fn, sq, md, sep = "\r")
  hit <- get0(key, mp_cache, inherits = FALSE)
  if (!is.null(hit)) return(hit)
  v <- tryCatch(matched_positions_of(sq, md, sp, aa, aa_tab),
                error = function(e) integer(0))
  assign(key, v, mp_cache)
  v
}
spec_ratio <- function(i, k) {                                 # 与 cand1 相比的特异比例
  if (is.null(lu) || !chim_flag[i] || length(cand_sm[[i]]) < k) return("")
  fn <- as.character(p$File_Name[i])
  raw <- unname(title_raw_scan(fn)["raw"])
  objs <- lu$objs[[raw]]
  if (is.null(objs)) return("")
  # （性能第二轮）：O(1) 命名向量索引（lu$ix 由 load_unlabeled_spectra
  # 提供），旧 match 线性扫保留为回退
  li <- if (!is.null(lu$ix) && !is.null(lu$ix[[raw]])) {
    unname(lu$ix[[raw]][fn])
  } else {
    match(fn, lu$titles[[raw]])
  }
  if (is.na(li)) return("")
  sp <- objs[[li]]@spectrum
  s1 <- cand_sm[[i]][[1]][1]; m1 <- cand_sm[[i]][[1]][2]
  sk <- cand_sm[[i]][[k]][1]; mk <- cand_sm[[i]][[k]][2]
  A1 <- mp_cached(fn, s1, m1, sp)
  if (length(A1) == 0) return("")
  Ak <- mp_cached(fn, sk, mk, sp)
  if (length(Ak) == 0) return("")
  rk <- length(setdiff(Ak, A1)) / length(Ak)
  r1 <- length(setdiff(A1, Ak)) / length(A1)
  sprintf("%.4f", min(rk, r1))
}

# ---- 逐行组装 94 列主表 ----
# ---- assemble the 94-column master table row by row
cand_raw <- strsplit(as.character(p$candidates), "|", fixed = TRUE)
cand_sm <- lapply(seq_len(n), function(i) {                    # [[seq,mod],[seq,mod],...]
  cc <- cand_raw[[i]]
  out <- list()
  for (s in cc) {
    s <- trimws(s)
    m <- regmatches(s, regexpr("[(].*[)]", s))
    out[[length(out) + 1L]] <- c(sub("[(].*$", "", s),
                                 if (length(m) == 1L) substr(m, 2, nchar(m) - 1L) else "")
  }
  out
})
# ---- 预计算全部候选序列的蛋白定位----
# 每个唯一序列只查一次并回填 loc_cache（GG→N 回退一并完成），主循环全部命中缓存。
# 此前用 18 worker 集群：仅 4499 条序列，集群启动与整库串×18 次序列化开销反而
# 占大头；独立计时实测 4411 条定位仅 3.5s（全部命中，regexpr 单线程毫秒级），
# 直接串行即可，worker 闭包复制（fmap_w/loc_w）随之删除。
# ---- precompute the protein location of every candidate sequence each unique sequence is resolved once (GG->N fallback included)
# and cached, so the main loop always hits the cache. The former 18-worker
# cluster was removed - for ~4.5k sequences its startup plus serializing the
# concatenated proteome to every worker dominated the runtime; a standalone
# measurement shows all 4,411 sequences locate in 3.5 s single-threaded
# (every one hits, so one regexpr scan each), hence the direct serial loop
# and the worker-closure duplicates (fmap_w/loc_w) are gone.
if (nzchar(fasta_file) && file.exists(fasta_file)) {
  all_sq <- unique(unlist(lapply(cand_sm, function(cs)
    vapply(cs, `[`, character(1), 1L))))
  all_sq <- all_sq[nzchar(all_sq)]
  if (length(all_sq) > 0) {
    invisible(lapply(all_sq, loc_of))   # loc_of 自行回填 loc_cache
    cat(sprintf("[主表导出] 蛋白定位预计算 %d 条唯一序列完成\n", length(all_sq)))
  }
}
# ---- 混合谱与非混合谱分道竞争--------------------
# 把每个冲突行的候选按"与 cand1 的特异离子比例"分成两道：比例 ≥0.1 的候选
# 为混合谱竞争者（D/E 道），其余为非混合谱竞争者（B/C 道）。两道各自按同一
# 判定（ΔScore≥τ 或理论谱图找回）独立决出胜者，并给出整体结论列，从而区分
# "混合谱且原始错误"（两道都改判）与"混合谱但原始正确"（非混合道保留原始）。
# ---- split-lane competition ----------------------
# Each conflict row's candidates are split by their specific-ion ratio against
# cand1: candidates at >=0.1 form the chimeric lane (D/E) and the rest form the
# non-chimeric lane (B/C). Each lane runs the same judgment (ΔScore>=tau or
# predicted-spectrum rescue) independently, and overall conclusion columns
# separate "chimeric AND original wrong" (both lanes correct the original) from
# "chimeric BUT original consistent" (the non-chimeric lane keeps the original).
spec_r <- vector("list", n)
for (i in which(conflict_row & chim_flag)) {
  v <- rep(0, max(1L, ncand[i]))
  if (ncand[i] >= 2) {
    for (k in 2:ncand[i]) {
      rk <- suppressWarnings(as.numeric(spec_ratio(i, k)))
      v[k] <- if (!is.na(rk)) rk else 0
    }
  }
  spec_r[[i]] <- v
}
nonchim_win <- rep(0L, n)
chim_win <- rep(0L, n)
for (i in which(conflict_row)) {
  s <- scores[[i]]
  if (length(s) < 2) next
  sr <- spec_r[[i]]
  if (length(sr) < length(s)) sr <- c(sr, rep(0, length(s) - length(sr)))
  idx_nc <- c(1L, which(sr[-1] < 0.1) + 1L)
  idx_c <- c(1L, which(sr[-1] >= 0.1) + 1L)
  for (lane in c("nonchim", "chim")) {
    idx <- if (lane == "nonchim") idx_nc else idx_c
    if (length(idx) < 2) next
    s2 <- s[idx]
    k2 <- which.max(s2[-1]) + 1L
    win_idx <- idx[k2]
    d <- s2[k2] - s2[1]
    prm <- LL[["predict"]][[i]]
    dpr <- NA_real_
    if (length(prm) >= max(1L, win_idx) && !is.na(prm[1]) && !is.na(prm[win_idx]))
      dpr <- prm[win_idx] - prm[1]
    j <- (!is.na(d) && d >= tau_v) ||
         (!is.na(d) && d < tau_v && !is.na(dpr) && dpr >= THETA)
    if (lane == "nonchim") nonchim_win[i] <- if (j) win_idx else 0L
    else chim_win[i] <- if (j) win_idx else 0L
  }
}

rows <- vector("list", n)
for (i in seq_len(n)) {
  r <- list(File_Name = as.character(p$File_Name[i]))
  # 冲突行的候选以成对报告 candidates 列为准（紧凑列表，与打分/指标 token 顺序
  # 一一对应；冲突报告 X1..X5 按策略列会有空槽，直接取会把候选与打分错位——
  # 二次修复，与谱图 PDF 判定同源）
  # the conflict-row candidates come from the paired report's candidates column
  # (the compacted list, aligned one-to-one with the score/metric tokens; the
  # conflict report's per-strategy X1..X5 columns have gaps that misalign the
  # candidates with the scores - fixed again same source as the
  # spectrum-PDF judgment)
  xs <- NULL; xm <- NULL
  if (conflict_row[i] && !is.null(pr)) {
    cc <- strsplit(as.character(pr$candidates[pr_i[i]]), "|", fixed = TRUE)[[1]]
    cc <- trimws(cc)
    xs <- sub("[(].*$", "", cc)
    xm <- vapply(cc, function(s) {
      m <- regmatches(s, regexpr("[(].*[)]", s))
      if (length(m) == 1L) substr(m, 2, nchar(m) - 1L) else ""
    }, character(1))
  }
  n_avail <- if (!is.null(xs)) length(xs) else length(cand_sm[[i]])
  for (k in seq_len(NSLOT)) {
    if (k > n_avail) {                 # 无该候选：整组置空，不复制最后一条
      r[[sprintf("cand%d_sequence", k)]] <- ""
      r[[sprintf("cand%d_modification", k)]] <- ""
      for (mn in c("ul_by", "lb_by", "ul_all", "lb_all", "pcos", "predict"))
        r[[sprintf("cand%d_%s", k, mn)]] <- ""
      r[[sprintf("cand%d_paired_state", k)]] <- ""
      r[[sprintf("cand%d_unlabeled_coverage", k)]] <- ""
      r[[sprintf("cand%d_labeled_scan", k)]] <- ""
      r[[sprintf("cand%d_score", k)]] <- ""
      if (k >= 2) {
        r[[sprintf("cand%d_delta_score", k)]] <- ""
        r[[sprintf("cand%d_pred_rescue", k)]] <- ""
        r[[sprintf("cand%d_specific_ratio_vs_cand1", k)]] <- ""
      }
      next
    }
    kk <- k
    if (!is.null(xs)) {
      sm <- c(if (!is.na(xs[kk])) xs[kk] else "",
              if (!is.na(xm[kk])) xm[kk] else "")
    } else {
      sm <- cand_sm[[i]][[kk]]
    }
    loc1 <- if (k == 1L) {
      v <- loc_of(sm[1])
      if (is.null(v)) c("NA", "-") else v
    } else c("NA", "-")
    r[[sprintf("cand%d_sequence", k)]] <- sm[1]
    r[[sprintf("cand%d_modification", k)]] <- sm[2]
    for (mn in c("ul_by", "lb_by", "ul_all", "lb_all", "pcos", "predict")) {
      v <- LL[[mn]][[i]]
      r[[sprintf("cand%d_%s", k, mn)]] <-
        if (length(v) >= kk && !is.na(v[kk])) sprintf("%.4f", v[kk]) else ""
    }
    r[[sprintf("cand%d_paired_state", k)]] <-
      if (length(mstL[[i]]) >= kk) tolower(trimws(mstL[[i]][kk])) else ""
    v <- covL[[i]]
    r[[sprintf("cand%d_unlabeled_coverage", k)]] <-
      if (length(v) >= kk && !is.na(v[kk])) sprintf("%.4f", v[kk]) else ""
    r[[sprintf("cand%d_labeled_scan", k)]] <-
      if (length(labL[[i]]) >= kk) trimws(labL[[i]][kk]) else ""
    if (k == 1L) {
      r[["cand1_protein"]] <- loc1[1]
      r[["cand1_location"]] <- loc1[2]
    }
    if (length(scores[[i]]) >= kk) r[[sprintf("cand%d_score", k)]] <-
      sprintf("%.6f", scores[[i]][kk]) else r[[sprintf("cand%d_score", k)]] <- ""
    if (k >= 2) {
      r[[sprintf("cand%d_delta_score", k)]] <-
        if (length(scores[[i]]) >= k && length(scores[[i]]) >= 1)
          sprintf("%.6f", scores[[i]][k] - scores[[i]][1]) else ""
      prm <- LL[["predict"]][[i]]
      r[[sprintf("cand%d_pred_rescue", k)]] <-
        if (length(prm) >= k && !is.na(prm[1]) && !is.na(prm[k]))
          ifelse(prm[k] - prm[1] >= THETA, "TRUE", "") else ""
      r[[sprintf("cand%d_specific_ratio_vs_cand1", k)]] <- spec_ratio(i, k)
    }
  }
  # 结论 16 列：仅冲突行判定；无差异行结论留空、corrected_* 复制候选 1
  # the 16 conclusion columns: judged for conflict rows only; non-conflict rows
  # keep blank conclusion columns with corrected_* copying candidate 1
  k2 <- min(rev_idx[i], ncand[i])
  if (!is.null(xs)) {
    sm1 <- c(if (!is.na(xs[1])) xs[1] else "", if (!is.na(xm[1])) xm[1] else "")
    smr <- if (k2 >= 2 && !is.na(xs[k2])) c(xs[k2], xm[k2]) else sm1
  } else {
    sm1 <- if (length(cand_sm[[i]]) >= 1) cand_sm[[i]][[1]] else c("", "")
    smr <- if (k2 >= 2) cand_sm[[i]][[k2]] else sm1
  }
  # 结论列：冲突行**全部**填写（含判保留/存疑行，便于筛选）；无差异行留空
  # conclusion columns: filled for **every** conflict row (including keep and
  # ambiguous-match rows so they can be filtered); blank for non-conflict rows
  r$rev_idx <- if (conflict_row[i] && k2 >= 2) k2 else ""
  r$max_delta_score <- if (conflict_row[i] && !is.na(delta[i])) sprintf("%.6f", delta[i]) else ""
  r$tau <- if (conflict_row[i]) sprintf("%.4f", tau_v) else ""
  r$gap <- if (conflict_row[i] && !is.na(delta[i])) sprintf("%.6f", tau_v - delta[i]) else ""
  r$call_base_on_score <- if (conflict_row[i] && !is.na(delta[i]) && delta[i] >= tau_v) "TRUE" else ""
  r$delta_pred <- if (conflict_row[i] && !is.na(dpred[i])) sprintf("%.6f", dpred[i]) else ""
  r$pred_rescue <- if (rescue[i]) "TRUE" else ""
  r[["call_base_on_score&pred"]] <- if (judged[i]) "TRUE" else ""
  r$ions_identical <- if (conflict_row[i] && ions_identical[i]) "TRUE" else ""
  r$ions_small_diff <- if (conflict_row[i] && ions_small_diff[i]) "TRUE" else ""
  r$coverage_below_20 <- if (cov_below20[i]) "TRUE" else ""
  locc <- if (judged[i]) {
    v <- loc_of(smr[1])
    if (is.null(v)) {
      v <- inherit_loc(sm1[1], smr[1], r[["cand1_location"]])
      if (is.null(v)) c("NA", "-") else c(r[["cand1_protein"]], v)
    } else v
  } else c(r[["cand1_protein"]], r[["cand1_location"]])
  r$corrected_protein <- locc[1]
  r$corrected_location <- locc[2]
  r$corrected_sequence <- if (judged[i]) smr[1] else sm1[1]
  r$corrected_modification <- if (judged[i]) smr[2] else sm1[2]
  r$corrected_labeled_scan <- if (judged[i] && k2 >= 2 && length(labL[[i]]) >= k2)
    trimws(labL[[i]][k2]) else ""
  # diff_type / mod_changed / aa_changed（仅判需校正行；与主表口径一致）
  # diff_type / mod_changed / aa_changed (judged rows only; master-table defs)
  r$mod_type_changed <- ""                                   # 修饰类型列（全行）
  if (!judged[i]) {
    r$diff_type <- ""; r$mod_changed <- ""; r$aa_changed <- ""
  } else {
    t1 <- mod_tokens(sm1[2]); t2 <- mod_tokens(smr[2])
    types_same <- setequal(mod_type(t1), mod_type(t2)) &&
      all(table(mod_type(t1))[names(table(mod_type(t2)))] == table(mod_type(t2)))
    swap <- nchar(sm1[1]) == nchar(smr[1]) && !identical(sm1[1], smr[1]) &&
      identical(sort(unlist(strsplit(sm1[1], ""))), sort(unlist(strsplit(smr[1], ""))))
    # 顺序互换序列长度相同、氨基酸组成与各自数量相同、
    # 仅顺序不同（sort 后多集相等即满足）
    # AA order change same length, same amino-acid
    # composition with the same counts, differing only in order (equal sorted
    # multisets)
    # N 端非特异性酶切corrected 肽段前一位氨基酸非 K/R
    # 且不在蛋白 N 端（corrected_location 首字符 ∉ {K,R,_,-}）。
    # v1.0判据按酶切类型参数化——C 端酶看前一残基是否 cut 残基
    # （Trypsin 即 K/R，与历史口径一致）；N 端酶（LysN/AspN/LysargiNase）看肽段
    # 首残基；无酶切恒特异。
    # N-terminal non-specific cleavage (user's definition): the residue before
    # the corrected peptide is neither K nor R and not at the protein N-terminus
    # (first char of corrected_location not in {K,R,_,-}). v1.0 
    # parameterized by enzyme - C-terminal enzymes check the previous residue
    # against the cut set (trypsin = K/R, historical parity); N-terminal enzymes
    # (LysN/AspN/LysargiNase) check the peptide's first residue; no enzyme is
    # always specific.
    loc_pre <- tryCatch(strsplit(as.character(r$corrected_location), ",", fixed = TRUE)[[1]][2],
                       error = function(e) "")
    nterm <- nzchar(loc_pre) &&
      !enzyme_nterm_ok(iso_enzyme(), loc_pre, substr(smr[1], 1, 1))
    # 真 SAP 判据看 "->" 目标是否为标准
    # 氨基酸名（单字母或三字母，如 Pro->Xle、Ala->Ser——其 [P]/[A] 是残基标记，非修饰
    # 证据）；目标是修饰名（如 pyro-Glu）则是修饰转换（N 端谷氨酰胺环化），归修饰类。
    # 注意不能用"有无方括号"区分：真 SAP 记法 7,Pro->Xle[P] 同样带方括号。
    # Real-SAP test the target
    # after "->" must be a standard amino-acid name (one- or three-letter, e.g.
    # Pro->Xle, Ala->Ser - their [P]/[A] tags mark the residue, not a
    # modification); a modification name as target (e.g. pyro-Glu) is an N-
    # terminal cyclisation conversion and belongs to the modification class.
    # Do NOT test for brackets: real-SAP tokens like 7,Pro->Xle[P] carry them too.
    tk_all <- c(t1, t2)
    tgt <- sub("[[].*$", "", sub("^.*->", "", tk_all))
    AA_names <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P",
                  "S","T","W","Y","V",
                  "Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
                  "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val",
                  "Xle","Sec","Pyl")
    sap <- any(grepl("->", tk_all, fixed = TRUE) & tgt %in% AA_names)
    dt <- if (chim_flag[i]) {
      if (diff_charge(as.character(p$File_Name[i]), smr[1], smr[2]))
        "chimera_different_charge" else "chimera_same_charge"
    }
      else if (sap) "sap"
      # 优先级修饰差异高于 N 端非特异性酶切——
      # 既带修饰差异又 N 端非特的行归修饰类，不挪入 N 端类
      # Priority modification difference outranks
      # N-terminal non-specificity - rows with both stay in the modification
      # class and are NOT moved to the N-terminal class
      else if (!types_same) "different_modification"
      else if (nterm) "n_terminal_non_specific_cleavage"
      else if (nchar(smr[1]) < 7) "peptide_length_less_than_7aa"
      # 修饰位置改变序列相同、修饰类型一致（前支
      # !types_same 已拦截类型差异），仅修饰位置存在差异
      # Modification location changed identical
      # sequence with the same modification types (type differences are caught
      # by the !types_same branch above), differing only in the positions
      else if (identical(sm1[1], smr[1]) && !identical(sm1[2], smr[2])) "modification_location_changed"
      else if (swap) "aa_order_change"
      else "same_weight_aa_pairs"      # 序列不同的替换
    r$diff_type <- dt
    mc <- mod_cancel(t1, t2); ac <- aa_cancel(sm1[1], smr[1])
    mcs <- function(v) if (length(v) == 0) "" else paste0(paste0(v, ";"), collapse = "")
    if (dt == "modification_location_changed") {
      r$mod_changed <- sprintf("%s -> %s", paste(t1, collapse = ""), paste(t2, collapse = ""))
      r$aa_changed <- ""
    } else {
      r$mod_changed <- if (length(mc$l) == 0 && length(mc$r) == 0) "" else
        sprintf("%s%s=%s%s", mcs(mc$l), paste(ac$l, collapse = ""),
                mcs(mc$r), paste(ac$r, collapse = ""))
      r$aa_changed <- if (length(mc$l) == 0 && length(mc$r) == 0)
        sprintf("%s=%s", paste(ac$l, collapse = ""), paste(ac$r, collapse = "")) else ""
    }
    # 修饰类型：mod_changed 按 [=;] 切分、保留含 [ 的 token（去氨基酸与等号）
    # modification types: split mod_changed on [=;] and keep the tokens that
    # contain "[" (amino-acid letters and the equals sign are dropped)
    toks <- trimws(strsplit(r$mod_changed, "[=;]")[[1]])
    toks <- toks[nzchar(toks)]
    r$mod_type_changed <- paste(toks[grepl("\\[", toks)], collapse = ";")
  }
  # 分道竞争列仅混合谱行填结论，其余留空
  # split-lane columns filled only for chimeric rows
  if (chim_flag[i] && conflict_row[i]) {
    # 冲突行候选以 pr$candidates（xs/xm）为准（与主装配同口径）；阶段 A 的
    # cand_sm 可能比冲突候选短（越界来源），仅在无 xs 时回退
    # conflict-row candidates follow pr$candidates (xs/xm, same as the main
    # assembly); stage-A cand_sm can be shorter than the conflict list (the
    # out-of-bounds source), so fall back only when xs is absent
    sm <- if (!is.null(xs) && length(xs) >= 2) {
      lapply(seq_along(xs), function(k) c(xs[k], xm[k]))
    } else cand_sm[[i]]
    pick <- function(w) {
      if (is.na(w) || w < 2 || w > length(sm)) return(c("", ""))
      sm[[w]]
    }
    pc_nc <- pick(nonchim_win[i])
    pc_c <- pick(chim_win[i])
    r$nonchim_lane_win_idx <- if (nonchim_win[i] >= 2 && nonchim_win[i] <= length(sm)) as.character(nonchim_win[i]) else ""
    r$nonchim_lane_corrected_sequence <- pc_nc[1]
    r$nonchim_lane_corrected_modification <- pc_nc[2]
    r$chim_lane_win_idx <- if (chim_win[i] >= 2 && chim_win[i] <= length(sm)) as.character(chim_win[i]) else ""
    r$chim_lane_partner_sequence <- pc_c[1]
    r$chim_lane_partner_modification <- pc_c[2]
    r$chim_and_original_wrong <- if (nonchim_win[i] >= 2 && chim_win[i] >= 2) "TRUE" else ""
    r$chim_but_original_consistent <- if (nonchim_win[i] == 0L) "TRUE" else ""
  } else {
    r$nonchim_lane_win_idx <- ""; r$nonchim_lane_corrected_sequence <- ""
    r$nonchim_lane_corrected_modification <- ""; r$chim_lane_win_idx <- ""
    r$chim_lane_partner_sequence <- ""; r$chim_lane_partner_modification <- ""
    r$chim_and_original_wrong <- ""; r$chim_but_original_consistent <- ""
  }
  rows[[i]] <- as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE)
}
mt <- do.call(rbind, rows)

# ---- 最小特异子离子比例列（cand2..5 与 cand1 相比的最小值；无竞争为空）----
# ---- the minimum specific-ion-ratio column (minimum of the cand2..5 vs cand1
# specific ratios; blank when there is no competitor)
spec_cols <- intersect(paste0("cand", seq_len(max(2L, NSLOT))[-1], "_specific_ratio_vs_cand1"), names(mt))
if (length(spec_cols) > 0) mt$min_specific_ratio <- apply(mt[, spec_cols, drop = FALSE], 1, function(x) {
  x <- x[nzchar(x)]
  if (length(x) == 0) "" else sprintf("%.4f", min(as.numeric(x)))
})

# ---- 列序定稿（与 8.衣藻全谱主表.txt 一致）----
# ---- finalize the column order (identical to 8.衣藻全谱主表.txt)
cand_block <- function(k) {
  v <- sprintf("cand%d_%s", k, c("sequence", "modification", "ul_by", "lb_by",
    "ul_all", "lb_all", "pcos", "predict", "paired_state",
    "unlabeled_coverage", "labeled_scan", "score"))
  if (k >= 2) v <- c(v, sprintf("cand%d_delta_score", k),
                     sprintf("cand%d_pred_rescue", k),
                     sprintf("cand%d_specific_ratio_vs_cand1", k))
  v
}
concl <- c("rev_idx", "max_delta_score", "tau", "gap", "call_base_on_score",
           "delta_pred", "pred_rescue", "call_base_on_score&pred",
           "ions_identical", "ions_small_diff", "coverage_below_20",
           "min_specific_ratio",
           "corrected_protein", "corrected_location",
           "corrected_sequence", "corrected_modification",
           "corrected_labeled_scan", "diff_type", "mod_changed",
           "mod_type_changed", "aa_changed",
           "nonchim_lane_win_idx", "nonchim_lane_corrected_sequence",
           "nonchim_lane_corrected_modification", "chim_lane_win_idx",
           "chim_lane_partner_sequence", "chim_lane_partner_modification",
           "chim_and_original_wrong", "chim_but_original_consistent")
ord <- c("File_Name", "cand1_protein", "cand1_location",
         unlist(lapply(seq_len(NSLOT), cand_block)), concl)
ord <- ord[ord %in% names(mt)]
mt <- mt[, ord]
# ---- Excel 防公式误判：定位列占位横线改成下划线 ----
# 定位列已是 pFind 原生格式 "pos,prev,next/"；此步把
# 历史遗留的 "-"/U+2212 占位统一归并为 "_" 作保险丝（新产出不再含 "-"，除非
# 整列缺失占位 "-"），并保持 Excel 公式安全。只作用于 _location 结尾的列。
# ---- Excel-safe export: in the location columns (format pre[pos]post) the
# placeholder hyphens for a missing N/C-terminal flanking residue become
# underscores (plain ASCII): Excel never parses "_" as a formula and text
# search works directly; any legacy Unicode minus U+2212 is merged in too.
# Only columns ending in _location are touched; numeric columns are not.
loc_cols <- grepl("_location$", names(mt))
for (j in which(loc_cols)) {
  mt[[j]] <- gsub("−", "_", mt[[j]], fixed = TRUE)
  mt[[j]] <- gsub("-", "_", mt[[j]], fixed = TRUE)
}

# ---- 新肽段五列证据/状态----
# 对原始鉴定蛋白不在用户 FASTA 内的行（新肽段部分）计算并随主表写出：
#   n_frag＝新肽段匹配子离子数（b/y 位点 1..n-1 中有命中的位数）；
#   closest_annotated_competitor＝ΔScore 最大（最强）的已注释竞争
#   （"序列|修饰|Δ"，Δ＝竞争者分数−新肽分数，正＝竞争者分数更高）；
#   novel_doubt / novel_status＝存疑详情与状态分类——判错行写具体校正
#   类型（diff_type，如 "chimera_same_charge"），其余照旧（chimeric /
#   ions_identical / close_annotated / consecutive_missing）；
#   credible_novel＝无存疑即 TRUE。已注释行这五列留空。
# ---- five novel-peptide evidence/status columns in the master itself
# ----
# computed for rows whose original-identification protein is NOT in the user
# FASTA and written WITH the master:
#   n_frag = matched fragment-ion count (b/y positions 1..n-1 with any hit);
#   closest_annotated_competitor = the strongest annotated competitor
#   ("seq|mod|Δ"; Δ = competitor score minus the novel score, positive =
#   competitor higher);
#   novel_doubt / novel_status = doubt detail / status class - judged-wrong
#   rows carry their specific correction type (diff_type, e.g.
#   "chimera_same_charge"); chimeric / ions_identical / close_annotated /
#   consecutive_missing as before;
#   credible_novel = TRUE when no doubt. Annotated rows stay blank.
fa_path2 <- NULL
if (!is.null(base_cfg$fasta) && length(base_cfg$fasta) > 0 && base_cfg$fasta[1] != "") {
  fa_path2 <- tryCatch(resolve_cfg_paths(base_cfg, "fasta")[1], error = function(e) NULL)
}
# out_f 在文件后段才定义（master_table.tsv 写出点），此处提前按 result_dir
# 构造等价目录（原 dirname(out_f) 会在 out_f 赋值前引用——两流程均在此崩溃）
out_dir2 <- file.path(result_dir, "2.labeled_result")
if (!is.null(fa_path2) && file.exists(fa_path2)) {
  # N-5FASTA 真·只读一次（旧实现此处与下方再各读一次；
  # 大库解析两遍纯浪费且两次之间文件改动会产生不一致）
  # N-5 the FASTA is genuinely read once (the old code read it
  # here AND again below; parsing a large library twice is pure waste and a
  # file change in between creates inconsistency)
  fa_obj <- tryCatch(read_fasta(fa_path2), error = function(e) character(0))
  fa_acc <- if (length(fa_obj) > 0) names(fa_obj) else character(0)
  cp2 <- as.character(mt$cand1_protein)
  # 设计约定：两分法以"输入原始鉴定"（a 槽蛋白）定侧——Part2
  # 必须包含全部新肽段输入行（如理论注入的 novel-kind 全 2,000 张），投票把
  # 个别行带向已注释候选不改变侧别；Part2 是在 Part1 行集上补充成对标记
  # 数值的报告，不是另做筛选
  # user ruling: the two-part split follows the ORIGINAL input
  # identification (a-slot protein) - Part2 must contain every novel input
  # row (e.g. all 2,000 novel-kind injections); a vote pointing some row at an
  # annotated candidate does not change its side. Part2 supplements the Part1
  # row set with the paired-label values rather than re-filtering it
  if (exists("mer", inherits = FALSE) && !is.null(mer) &&
      "File_Name" %in% names(mer)) {
    # 原始蛋白列：多引擎合并表用 "a.Proteins"，纯 pFind 合并表为 "Proteins"
    ap_col <- if ("a.Proteins" %in% names(mer)) "a.Proteins"
              else if ("Proteins" %in% names(mer)) "Proteins" else NULL
    if (!is.null(ap_col)) {
      ap <- as.character(mer[[ap_col]][match(as.character(mt$File_Name),
                                             as.character(mer[["File_Name"]]))])
      # 污染蛋白（CON_ 前缀，cRAP 类）属已注释搜索的校正域，不进新肽段竞争
      con_flag <- grepl("^CON_", ap)
      ap <- sub("^(ANN_|CON_)", "", ap)   # 加标库 accession 前缀与参考库归一
    } else {
      ap <- rep(NA_character_, nrow(mt))
      con_flag <- rep(FALSE, nrow(mt))
    }
    ap <- sub("/.*$", "", ap)
    cp2o <- ifelse(is.na(ap) | ap == "", cp2, ap)
  } else {
    cp2o <- cp2
    con_flag <- rep(FALSE, nrow(mt))
  }
  in_fasta <- (!is.na(cp2o) & cp2o != "NA" & cp2o != "-" & cp2o != "" &
                 cp2o %in% fa_acc) | con_flag
  nfrag <- character(nrow(mt))
  closest <- character(nrow(mt))
  closest_d <- rep(-Inf, nrow(mt))
  doubt <- character(nrow(mt))
  maxrun <- integer(nrow(mt))
  # 已注释序列空间＝用户 FASTA 全部序列拼接（子串匹配判定；I/L 归一）
  # P2-2/N-7分隔符 "|"→"\n"（蛋白序列不含换行，杜绝跨分隔符
  # 巧合命中）；fa_obj 已于块首读入（N-5）；annotated_seq 加环境缓存（同一肽段
  # 在数万张谱中反复出现，去重后判定次数下降一个量级）。
  # P2-2/N-7 separator "|" -> "\n" (protein sequences never
  # contain newlines, removing coincidental cross-boundary hits); fa_obj was
  # already read at the top of this block (N-5); annotated_seq gains an
  # environment cache (the same peptide recurs across tens of thousands of
  # spectra, so dedup cuts the judgment count by an order of magnitude).
  fa_seqs <- unname(as.character(fa_obj))
  fa_all <- paste(toupper(fa_seqs), collapse = "\n")
  ann_cache <- new.env(parent = emptyenv(), hash = TRUE)
  annotated_seq <- function(s) {
    s <- toupper(s)
    if (exists(s, envir = ann_cache, inherits = FALSE))
      return(get(s, envir = ann_cache, inherits = FALSE))
    v <- grepl(s, fa_all, fixed = TRUE) || grepl(gsub("I", "L", s), fa_all, fixed = TRUE)
    assign(s, v, envir = ann_cache)
    v
  }
  seq_cols <- grep("^cand[0-9]+_sequence$", colnames(mt), value = TRUE)
  seq_cols <- seq_cols[seq_cols != "cand1_sequence"]
  nov_rows <- which(!in_fasta)
  nov_raws <- unique(vapply(as.character(mt$File_Name[nov_rows]),
                            function(t) unname(title_raw_scan(t)["raw"]), character(1)))
  lu2 <- NULL
  if (length(nov_raws) > 0) {
    lu2 <- tryCatch(load_unlabeled_spectra(all_mgf_dir, only = nov_raws),
                    error = function(e) NULL)
  }
  for (k in nov_rows) {
    fn <- as.character(mt$File_Name[k])
    s1 <- as.character(mt$cand1_sequence[k])
    # 谱图定位（判错/嵌合行同样补证据列）
    sp <- NULL
    if (!is.null(lu2) && nzchar(s1) && !is.na(s1)) {
      raw <- unname(title_raw_scan(fn)["raw"])
      objs <- lu2$objs[[raw]]
      if (!is.null(objs)) {
        # （性能第二轮）：O(1) 命名向量索引，旧 match 线性扫保留为回退
        li <- if (!is.null(lu2$ix) && !is.null(lu2$ix[[raw]])) {
          unname(lu2$ix[[raw]][fn])
        } else {
          match(fn, lu2$titles[[raw]])
        }
        if (is.na(li)) {
          key <- title_scan_key(fn)
          ci <- which(startsWith(lu2$titles[[raw]], paste0(key, ".")))
          if (length(ci) > 0) li <- ci[1]
        }
        if (!is.na(li)) sp <- objs[[li]]@spectrum
      }
    }
    # 证据列①：新肽段匹配子离子数
      md1 <- as.character(mt$cand1_modification[k])
    pos_hit <- if (!is.null(sp)) {
      mp_cached(fn, s1, md1, sp)
    } else integer(0)
    if (!is.null(sp)) nfrag[k] <- as.character(length(unique(pos_hit)))
    # 证据列②：打分最相近的已注释竞争（ΔScore 最大者＝最强竞争者）
    for (ci in seq_along(seq_cols)) {
      cs <- as.character(mt[k, seq_cols[ci]])
      if (is.na(cs) || !nzchar(cs) || cs == s1) next
      if (!annotated_seq(cs)) next
      num <- sub("_sequence$", "", seq_cols[ci])
      dv <- suppressWarnings(as.numeric(as.character(mt[k, paste0(num, "_delta_score")])))
      if (length(dv) != 1 || is.na(dv)) next
      if (dv > closest_d[k]) {
        closest_d[k] <- dv
        closest[k] <- sprintf("%s|%s|%.6g", cs,
                              as.character(mt[k, paste0(num, "_modification")]), dv)
      }
    }
    # ③ 状态判定（优先级：判错 > 嵌合 > 峰集全同 > 已注释有力竞争 >
    #    连续缺失≥2 > 低分卡值）
    if (judged[k]) {
      dt <- as.character(mt$diff_type[k])
      if (is.na(dt) || !nzchar(dt)) dt <- "judged_wrong"
      doubt[k] <- dt
      next
    }
    if (chim_flag[k]) { doubt[k] <- "chimeric"; next }
    # 规则 1：竞争候选（cand2..N）的子离子匹配峰集与 cand1 完全一致
    if (!is.null(sp)) {
      A <- tryCatch(matched_peak_set(s1, md1, sp, aa, aa_tab), error = function(e) integer(0))
      comps <- unique(as.character(unlist(mt[k, seq_cols], use.names = FALSE)))
      comps <- comps[!is.na(comps) & nzchar(comps) & comps != s1]
      if (length(A) > 0) {
        for (s2 in comps) {
          B <- tryCatch(matched_peak_set(s2, "", sp, aa, aa_tab), error = function(e) integer(0))
          if (length(B) > 0 && setequal(A, B)) { doubt[k] <- paste0("ions_identical:", s2); break }
        }
      }
    }
    if (nzchar(doubt[k])) next
    # 规则 2：打分最相近的已注释竞争 |ΔScore| ≤ τ → 无法排除已注释解释
    tau_k <- suppressWarnings(as.numeric(as.character(mt$tau[k])))
    if (nzchar(closest[k]) && !is.na(tau_k) && abs(closest_d[k]) <= tau_k) {
      doubt[k] <- paste0("close_annotated:", closest[k])
    }
    if (nzchar(doubt[k])) next
    # 规则 3：连续 ≥2 个位点无子离子命中（≈ 连续 ≥3 个残基无子离子证据；
    # 设计约定收紧：GDGELDEELAK 类中间缺段的谱也应排除）
    npos <- nchar(s1) - 1L
    if (npos >= 2 && length(pos_hit) < npos) {
      missing <- setdiff(seq_len(npos), pos_hit)
      ml <- 0L; cur <- 0L; prev <- -99L
      for (p in missing) {
        cur <- if (p == prev + 1L) cur + 1L else 1L
        if (cur > ml) ml <- cur
        prev <- p
      }
      maxrun[k] <- ml
      if (ml >= 2L) doubt[k] <- sprintf("consecutive_missing:%d", ml)
    }
    if (nzchar(doubt[k])) next
    # 规则 4：低质量谱卡值
    # cand1_score < 0.2 → low_score（谱图级，与该肽段鉴定到几张谱无关）
    sc_k <- suppressWarnings(as.numeric(as.character(mt$cand1_score[k])))
    if (!is.na(sc_k) && sc_k < 0.2) doubt[k] <- sprintf("low_score:%.4f", sc_k)
  }
  mt$n_frag <- nfrag
  mt$closest_annotated_competitor <- closest
  mt$novel_doubt <- doubt
  mt$novel_status <- vapply(doubt, function(d) if (!nzchar(d)) "" else sub(":.*$", "", d),
                            character(1))
  mt$credible_novel <- ifelse(doubt == "" & !in_fasta, "TRUE", "")
  # 补充参数列
  #   paired_evidence＝成对证据一键筛选（paired 且余弦>0.7 → "1"）；
  #   max_missing_run＝最大连续缺失位点数（数值列，卡值 ≥2）；
  #   rank_in_peptide / n_spectra_total＝同肽段内按分数降序、子离子数降序
  #   的排名与该肽段总谱数（取子离子缺失最少/分数最高的报告）。
  # extra parameter columns (user ruling: the report file is
  # self-sufficient for step-by-step Excel filtering):
  #   paired_evidence = one-click paired-evidence filter (paired AND cos>0.7 -> "1");
  #   max_missing_run = the longest run of consecutive positions without any
  #     matched ion (numeric; the cutoff is >= 2);
  #   rank_in_peptide / n_spectra_total = the row's rank within its peptide
  #     (by score desc, then fragment count desc) and that peptide's total
  #     spectrum count (report the best-spectrum row first).
  mt$paired_evidence <- ifelse(as.character(mt$cand1_paired_state) == "paired" &
    suppressWarnings(as.numeric(as.character(mt$cand1_pcos))) > 0.7, "1", "")
  mt$max_missing_run <- maxrun
  seq_f <- as.character(mt$cand1_sequence)
  sc_all <- suppressWarnings(as.numeric(as.character(mt$cand1_score)))
  nf_all <- suppressWarnings(as.numeric(mt$n_frag))
  ord_all <- order(seq_f, -ifelse(is.na(sc_all), -Inf, sc_all),
                   -ifelse(is.na(nf_all), -Inf, nf_all), seq_along(seq_f))
  rk <- ave(seq_along(ord_all), seq_f[ord_all], FUN = seq_along)[order(ord_all)]
  nt <- ave(seq_along(ord_all), seq_f[ord_all], FUN = length)[order(ord_all)]
  mt$rank_in_peptide <- as.integer(rk)
  mt$n_spectra_total <- as.integer(nt)
} else {
  in_fasta <- rep(TRUE, nrow(mt))
}

out_f <- file.path(result_dir, "2.labeled_result", "master_table.tsv")
write.table(mt, out_f, sep = "\t", row.names = FALSE, quote = FALSE,
            fileEncoding = "UTF-8", na = "")
cat(sprintf("[主表导出] master_table.tsv：%d 行 %d 列（冲突行 %d；无差异行 %d；%s）\n",
            nrow(mt), ncol(mt), sum(conflict_row), sum(!conflict_row),
            ifelse(labeled_run, "标记样本：成对标记列已填", "无标记样本：标记列留空")))
cat(sprintf("[master table] master_table.tsv: %d rows x %d columns (conflict %d; non-conflict %d; %s)\n",
            nrow(mt), ncol(mt), sum(conflict_row), sum(!conflict_row),
            ifelse(labeled_run, "labeled: paired columns filled", "unlabeled: labeled columns blank")))

# ---- 判校正行的成对标记证据统计----
# 校正后保留序列＝cand1；统计其谱图有成对标记对应体（paired_state）及余弦>0.7 的数量，
# 按合并五类（混合谱/序列替换/修饰改变/SAP/N 端非特异）分列，写入同名 summary 文件。
# ---- Paired-label evidence among judged-corrected rows (user ruling
#: the master-table export must always report this alongside):
# the post-correction retained sequence is cand1; count spectra with a paired
# labeled counterpart (paired_state) and with paired cosine > 0.7, broken down
# by the merged five error categories; written next to the master table.
if (labeled_run && "diff_type" %in% names(mt)) {
  cd <- mt[!is.na(mt$diff_type) & nzchar(as.character(mt$diff_type)),, drop = FALSE]
  if (nrow(cd) > 0) {
    st <- as.character(cd$cand1_paired_state)
    pc <- suppressWarnings(as.numeric(cd$cand1_pcos))
    cat5 <- ifelse(cd$diff_type == "chimera_same_charge", "chimera_same_charge",
      ifelse(cd$diff_type == "chimera_different_charge", "chimera_different_charge",
      ifelse(cd$diff_type == "n_terminal_non_specific_cleavage", "n_term",
      ifelse(cd$diff_type %in% c("different_modification", "modification_location_changed"), "modification_change",
      ifelse(cd$diff_type %in% c("same_weight_aa_pairs", "peptide_length_less_than_7aa",
                                 "aa_order_change"),
             "sequence_substitution",
      ifelse(cd$diff_type == "sap", "sap", "other"))))))  # 未知类型显式归 other，不再误入 sap
    sum_f <- file.path(result_dir, "2.labeled_result", "correction_paired_summary.txt")
    con <- file(sum_f, open = "wt", encoding = "UTF-8")
    hdr <- sprintf("corrections_total\tpaired_labeled\tpaired_cos_gt_0.7\n%d\t%d\t%d\n",
                   nrow(cd), sum(st == "paired", na.rm = TRUE),
                   sum(st == "paired" & pc > 0.7, na.rm = TRUE))
    writeLines(hdr, con)
    writeLines("category\tn\tpaired_labeled\tpaired_cos_gt_0.7", con)
    for (g in c("chimera_same_charge", "chimera_different_charge", "sequence_substitution", "modification_change", "sap", "n_term")) {
      i <- which(cat5 == g)
      writeLines(sprintf("%s\t%d\t%d\t%d", g, length(i),
                         sum(st[i] == "paired", na.rm = TRUE),
                         sum(st[i] == "paired" & pc[i] > 0.7, na.rm = TRUE)), con)
    }
    close(con)
    cat(sprintf("[主表导出] 判校正 %d 行：成对标记 %d（余弦>0.7 %d）-> %s\n",
                nrow(cd), sum(st == "paired", na.rm = TRUE),
                sum(st == "paired" & pc > 0.7, na.rm = TRUE), sum_f))
    cat(sprintf("[master table] corrected rows %d: paired-labeled %d (cos>0.7 %d) -> %s\n",
                nrow(cd), sum(st == "paired", na.rm = TRUE),
                sum(st == "paired" & pc > 0.7, na.rm = TRUE), sum_f))
  }
}

# ---- 两分法结果呈现----
# 第 1 部分 1.annotation_correction_report.tsv ＝ 注释校正报告：原始鉴定蛋白在
#   用户 FASTA 内的行——判校正行（diff_type 非空）即错误匹配清单，判保留/存疑行
#   附全部证据列；
# 第 2 部分 2.novel_peptide_competition_report.tsv ＝ 新肽段竞争报告：原始鉴定
#   蛋白不在用户 FASTA 内的行（含无法定位蛋白的 "NA" 行——肽段未落在用户 FASTA
#   即新肽候选），证据列 n_frag / closest_annotated_competitor / novel_doubt /
#   novel_status / credible_novel 已随主表计算；
# 两分法摘要 two_part_report_summary.txt 供网页报告呈现两个面板。
# ---- two-part result presentation ----
# Part 1 1.annotation_correction_report.tsv = the annotation-correction report:
#   rows whose original-identification protein is IN the user FASTA -
#   judged-corrected rows (non-empty diff_type) form the misidentification list;
#   keep/ambiguous rows carry all the evidence columns;
# Part 2 2.novel_peptide_competition_report.tsv = the novel-peptide competition
#   report: rows whose original-identification protein is NOT in the user FASTA
#   (including "NA"-protein rows - a peptide that did not localise in the user
#   FASTA is a novel candidate); the evidence columns n_frag /
#   closest_annotated_competitor / novel_doubt / novel_status / credible_novel
#   were computed with the master;
# two_part_report_summary.txt feeds the two web-report panels.
if (!is.null(fa_path2) && file.exists(fa_path2)) {
  # 第 1 部分：注释校正报告
  p1 <- mt[in_fasta,, drop = FALSE]
  out1 <- file.path(out_dir2, "1.annotation_correction_report.tsv")
  write.table(p1, out1, sep = "	", row.names = FALSE, quote = FALSE,
              fileEncoding = "UTF-8", na = "")
  n1_err <- sum(nzchar(as.character(p1$diff_type)), na.rm = TRUE)
  # 第 2 部分：新肽段竞争报告
  p2 <- mt[!in_fasta,, drop = FALSE]
  # 设计约定：嵌合谱图判定为新肽段谱图中的错误谱图，加标注列
  # （结果不修改——行归属与可信判定不变，仅显式标注错误谱图性质）
  # user ruling: chimeric spectra are judged WRONG spectra among
  # the novel-peptide spectra - add an annotation column (the results are
  # unchanged: row assignment and credibility stay as-is, the wrong-spectrum
  # nature is just labelled explicitly)
  p2$spectrum_verdict <- ifelse(p2$novel_status %in% c("chimeric", "chimera_same_charge"),
                                "wrong_spectrum_mixture", "")
  out2 <- file.path(out_dir2, "2.novel_peptide_competition_report.tsv")
  write.table(p2, out2, sep = "	", row.names = FALSE, quote = FALSE,
              fileEncoding = "UTF-8", na = "")
  n2_cred <- sum(p2$credible_novel == "TRUE")
  n2_wrong_sp <- sum(p2$spectrum_verdict != "")
  n2_err <- sum(nzchar(as.character(p2$diff_type)), na.rm = TRUE)
  n2_doubt <- sum(nzchar(as.character(p2$novel_doubt)), na.rm = TRUE)
  n2_close <- sum(p2$novel_status == "close_annotated")
  cat(sprintf("[两分法报告] Part2 嵌合错误谱图标注：%d 行（混合谱判定，结果不变）\n",
              n2_wrong_sp))
  cat(sprintf("[two-part report] Part2 chimeric wrong-spectrum annotation: %d rows (mixture ruling, results unchanged)\n",
              n2_wrong_sp))
  # ---- pBuild 谱图清单可信新肽段的无标记行（pFind
  # 格式字段取自原始搜库结果）＋成对标记行（Label=2、标题/电荷/前体质量换成
  # 标记对应体），供 pBuild 画谱。无标记行按 psm_file 的 File_Name 回查原行。----
  # pBuild spectrum list the credible novel peptides'
  # unlabeled rows (pFind-format fields from the original search result) plus the
  # paired labeled rows (Label=2, title/charge/precursor swapped to the labeled
  # counterpart) for pBuild spectrum rendering. Unlabeled rows are looked up in
  # psm_file by File_Name.
  p2_cred <- p2[p2$credible_novel == "TRUE",, drop = FALSE]
  psm_src <- if (!is.null(base_cfg$psm_file) && length(base_cfg$psm_file) > 0 &&
                 base_cfg$psm_file[1] != "") {
    tryCatch(resolve_cfg_paths(base_cfg, "psm_file")[1], error = function(e) "")
  } else ""
  if (nrow(p2_cred) > 0 && nzchar(psm_src) && file.exists(psm_src)) {
    orig_psm <- tryCatch(read.delim(psm_src, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE, check.names = FALSE,
                                    quote = "", comment.char = ""), error = function(e) NULL)
    if (!is.null(orig_psm) && "File_Name" %in% colnames(orig_psm)) {
      om <- orig_psm[match(as.character(p2_cred$File_Name),
                           as.character(orig_psm$File_Name)),, drop = FALSE]
      lab_t <- character(nrow(p2_cred))
      if (exists("pr", inherits = FALSE) && !is.null(pr) &&
          "label_File_names" %in% colnames(pr) && exists("pr_i", inherits = FALSE)) {
        ii <- pr_i[!in_fasta]
        ii <- ii[p2$credible_novel == "TRUE"]
        lt <- tryCatch(strsplit(as.character(pr$label_File_names[ii]), ";", fixed = TRUE),
                       error = function(e) NULL)
        if (!is.null(lt)) {
          # R 4.5 实测 nzchar(NA) 返回 TRUE：pr_i 无对应成对报告的 NA 行会漏过
          # 下方 nzchar 守卫并以 "" 降级；长度不足补 "" 防越界 NA
          # nzchar(NA) is TRUE in R 4.5: NA entries (no paired-report row via
          # pr_i) must be demoted to "" here, and any shortfall padded, so the
          # loop below never sees NA
          lab_t <- vapply(lt, function(x)
            if (length(x) > 0 && !is.na(x[1])) x[1] else "", character(1))
          if (length(lab_t) < nrow(p2_cred))
            lab_t <- c(lab_t, rep("", nrow(p2_cred) - length(lab_t)))
        }
      }
      rows_out <- om
      n_lab <- 0L
      # （性能第二轮）：标记行预分配列表、循环后一次 rbind（原逐行
      # rbind 增长 O(n²)）
      lab_rows <- vector("list", nrow(p2_cred))
      for (j in seq_len(nrow(p2_cred))) {
        if (!isTRUE(nzchar(lab_t[j]))) next
        lr <- om[j,, drop = FALSE]
        lr$File_Name <- lab_t[j]
        if ("Label" %in% colnames(lr)) lr$Label <- sub("^1", "2", as.character(lr$Label))
        ch <- sub("^.*\\.(\\d)\\.\\d\\.dta$", "\\1", lab_t[j])
        if (!is.na(ch) && ch != lab_t[j] && "Charge" %in% colnames(lr)) lr$Charge <- ch
        if (!is.null(lu2)) {
          raw2 <- unname(title_raw_scan(lab_t[j])["raw"])
          objs2 <- lu2$objs[[raw2]]
          if (!is.null(objs2)) {
            # （性能第二轮）：O(1) 命名向量索引
            li2 <- if (!is.null(lu2$ix) && !is.null(lu2$ix[[raw2]])) {
              unname(lu2$ix[[raw2]][lab_t[j]])
            } else {
              match(lab_t[j], lu2$titles[[raw2]])
            }
            if (!is.na(li2)) {
              pm2 <- tryCatch(objs2[[li2]]@precursorMass, error = function(e) NA_real_)
              if (is.finite(pm2) && "Exp.MH." %in% colnames(lr)) lr$Exp.MH. <- pm2
            }
          }
        }
        lab_rows[[n_lab + 1L]] <- lr
        n_lab <- n_lab + 1L
      }
      if (n_lab > 0) rows_out <- rbind(rows_out, do.call(rbind, lab_rows[seq_len(n_lab)]))
      pb_path <- file.path(out_dir2, "novel_spectra_pBuild.spectra")
      write.table(rows_out, pb_path, sep = "\t", row.names = FALSE, quote = FALSE,
                  fileEncoding = "UTF-8", na = "")
      cat(sprintf("[两分法报告] pBuild 谱图清单：%d 行（无标记 %d＋标记 %d）-> %s\n",
                  nrow(rows_out), nrow(p2_cred), n_lab, pb_path))
      cat(sprintf("[two-part report] pBuild spectrum list: %d rows (unlabeled %d + labeled %d) -> %s\n",
                  nrow(rows_out), nrow(p2_cred), n_lab, pb_path))
    }
  }
  # 两分法摘要（纯 ASCII，网页报告按 part 名渲染中英标题）
  sum2 <- file.path(out_dir2, "two_part_report_summary.txt")
  con2 <- file(sum2, open = "wt", encoding = "UTF-8")
  writeLines(c("part	rows	headline",
               sprintf("1.annotation_correction_report	%d	annotation corrections (misidentifications): %d rows",
                       nrow(p1), n1_err),
               sprintf("2.novel_peptide_competition_report	%d	credible novel peptides: %d rows (doubt %d; judged wrong %d; close annotated competition %d)",
                       nrow(p2), n2_cred, n2_doubt, n2_err, n2_close)), con2)
  close(con2)
  cat(sprintf("[两分法报告] Part1 注释校正：%d 行（判校正 %d）-> %s
",
              nrow(p1), n1_err, out1))
  cat(sprintf("[two-part report] part1 annotation correction: %d rows (corrected %d) -> %s
",
              nrow(p1), n1_err, out1))
  cat(sprintf("[两分法报告] Part2 新肽段竞争：%d 行（可信 %d；判错 %d；已注释有力竞争 %d）-> %s
",
              nrow(p2), n2_cred, n2_err, n2_close, out2))
  cat(sprintf("[two-part report] part2 novel peptide competition: %d rows (credible %d; judged wrong %d; close annotated competition %d) -> %s
",
              nrow(p2), n2_cred, n2_err, n2_close, out2))
} else {
  cat("[两分法报告] 未找到用户 FASTA（config3 fasta 缺失或路径不存在），跳过两分法报告
")
  cat("[two-part report] user FASTA missing (config3 fasta empty or absent); two-part report skipped
")
  # 自查整改跳过时清掉上一轮遗留的两分法三件套，防止后续
  # 谱图 PDF 模块读到旧一轮的 credible_novel 集
  for (stale in c(file.path(out_dir2, "1.annotation_correction_report.tsv"),
                  file.path(out_dir2, "2.novel_peptide_competition_report.tsv"),
                  file.path(out_dir2, "two_part_report_summary.txt"))) {
    if (file.exists(stale)) file.remove(stale)
  }
}
