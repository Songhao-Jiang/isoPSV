# =============================================================================
# unlabeled_two_part_report.R - 无标记样本两分法报告
#
# 与标记样本同名同构地产出两份报告（文件名、目录、摘要文件、网页呈现一致）：
#   第 1 部分 1.annotation_correction_report.tsv ＝ 已注释矫正报告：
#     原始鉴定蛋白在用户 FASTA 内的行——判错行（vote_says_wrong=1）即错误
#     匹配清单（diff_type = error_reason），其余行附全部证据列；
#   第 2 部分 2.novel_peptide_competition_report.tsv ＝ 新肽段竞争报告：
#     原始鉴定蛋白不在用户 FASTA 内的行（含空蛋白/无法定位行），附 denovo
#     库外标记（<slot>.DB_match=no_match）与谱图判定列 spectrum_verdict。
#
# 与标记版的差异（如实声明，不伪造证据）：
#   - 判定来源：classification_A.tsv（阶段 A 逐谱结论）＋ merge 表槽列蛋白；
#   - novel_doubt / novel_status / credible_novel 与标记版同标准计算
#     （设计约定修订：可信新肽段判定与标记与否无关——六重存疑
#     规则的数据源换为 classification_A 逐谱候选＋第一步无标记谱；
#     paired_evidence 等成对专属列不适用、不产出）；
#   - spectrum_verdict 与标记版同口径（混合谱 = wrong_spectrum_mixture）。
#
# 输入环境变量：CONFIG3（第三步参数文件）、RESULT_DIR（结果目录）。
# 输出（v1.0 设计约定：无标记版落在 report 路径，即 1.unlabeled_result/
# 报告夹；网页两分法面板已加目录回退）：
#   <result>/1.unlabeled_result/{1.annotation_correction_report.tsv,
#   2.novel_peptide_competition_report.tsv, two_part_report_summary.txt}
# =============================================================================
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
if (!dir.exists(script_dir)) {
  script_dir <- dirname(if (grepl("^[A-Za-z]:", script_arg)) script_arg
                        else file.path(getwd(), script_arg))
  script_dir <- gsub("\\\\", "/", script_dir)
}
old_wd <- setwd(script_dir)
source("compat.R")
source("ReadMGF_script.R")
setwd(old_wd)

cfg3_path <- Sys.getenv("CONFIG3")
if (!nzchar(cfg3_path)) stop("CONFIG3 not set (set by user_run.R)")
cfg3 <- read_param_file(cfg3_path)
cfg_dir <- cfg3$config_dir
if (is.null(cfg_dir) || length(cfg_dir) == 0 || !nzchar(cfg_dir[1]))
  stop("config_step3 的 config_dir 缺失或为空，请检查配置（或删除 result_* 后重跑）\nconfig_dir missing/empty in config_step3")
result_dir <- get_result_dir(cfg_dir)
base_cfg <- read_param_file(resolve_cfg_paths(cfg3, "base_config")[1])
fasta_file <- if (!is.null(base_cfg$fasta) && length(base_cfg$fasta) > 0)
  tail(resolve_cfg_paths(base_cfg, "fasta"), 1) else ""

cls_f <- file.path(result_dir, "1.unlabeled_result", "classification_A.tsv")
mrg_f <- file.path(result_dir, "merge", "1.combine_result.spectra")
if (!file.exists(cls_f))
  stop("classification_A.tsv 缺失，无法生成无标记两分法报告：", cls_f,
       "\nclassification_A.tsv missing, cannot build the unlabeled two-part report: ", cls_f)
cls <- read.delim(cls_f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "", na.strings = character(0))
cls[] <- lapply(cls, function(x) if (is.logical(x)) rep("", length(x)) else x)
cat("classification_A rows:", nrow(cls), "\n")

# 原始鉴定蛋白归属：merge 表槽 a 的 Proteins（"ACC/..."）
proteins <- rep("", nrow(cls))
a_seqs <- rep("", nrow(cls))
denovo_seq <- rep("", nrow(cls))
denovo_db  <- rep("", nrow(cls))
if (file.exists(mrg_f)) {
  mrg <- read.delim(mrg_f, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "", na.strings = character(0))
  mrg[] <- lapply(mrg, function(x) if (is.logical(x)) rep("", length(x)) else x)
  if (!"File_Name" %in% colnames(mrg))
    stop("merge 表缺 File_Name 列：", mrg_f, "\nmerge table lacks File_Name: ", mrg_f)
  mm <- match(as.character(cls$File_Name), as.character(mrg$File_Name))
  # 合并表槽 a 蛋白列就是基表 Proteins（无
  # a.Proteins 列），旧条件恒假致报告该列整列恒空；按 master_table_export
  # 同口径回退 Proteins
  if ("a.Proteins" %in% colnames(mrg)) proteins <- as.character(mrg$a.Proteins)[mm]
  else if ("Proteins" %in% colnames(mrg)) proteins <- as.character(mrg$Proteins)[mm]
  # 槽 a 肽序列（序列级参考库匹配用）
  if ("a.Sequence" %in% colnames(mrg)) a_seqs <- as.character(mrg$a.Sequence)[mm]
  # denovo 槽：任一 <slot>.DB_match 列（取第一个）
  db_cols <- grep("^[a-z]\\.DB_match$", colnames(mrg), value = TRUE)
  if (length(db_cols) > 0) {
    slot <- sub("\\.DB_match$", "", db_cols[1])
    s_col <- paste0(slot, ".Sequence")
    denovo_db <- as.character(mrg[[db_cols[1]]])[mm]
    if (s_col %in% colnames(mrg)) denovo_seq <- as.character(mrg[[s_col]])[mm]
  }
} else {
  cat("[提示] merge 表缺失（", mrg_f, "），蛋白归属列留空、按新肽候选侧处理\n", sep = "")
}
# FASTA 登录名集合（与标记版同口径：原始蛋白登录名命中用户 FASTA 即"已注释"）
fa_acc <- character(0)
if (nzchar(fasta_file) && file.exists(fasta_file)) {
  fa_acc <- unique(names(read_fasta(fasta_file)))
}

# 已注释序列空间＝用户 FASTA 全部序列拼接（I/L 归一子串匹配，与主表同口径）
fa_obj <- if (nzchar(fasta_file) && file.exists(fasta_file))
  tryCatch(read_fasta(fasta_file), error = function(e) character(0)) else character(0)
fa_all <- paste(toupper(unname(as.character(fa_obj))), collapse = "\n")
ann_cache <- new.env(parent = emptyenv(), hash = TRUE)
annotated_seq_ul <- function(s) {
  s <- toupper(s)
  if (!nzchar(s)) return(FALSE)   # 空串不参与子串匹配（grepl("") 恒真）
  if (exists(s, envir = ann_cache, inherits = FALSE))
    return(get(s, envir = ann_cache, inherits = FALSE))
  v <- grepl(s, fa_all, fixed = TRUE) || grepl(gsub("I", "L", s), fa_all, fixed = TRUE)
  assign(s, v, envir = ann_cache)
  v
}
acc_of <- function(pr) {
  if (is.na(pr) || !nzchar(pr) || pr %in% c("-", "NA")) return(character(0))
  p <- trimws(strsplit(pr, "/", fixed = TRUE)[[1]])
  p[nzchar(p)]
}
# accession 之外补序列级匹配（六框库 orf| 与
# 参考库 WP_/NP_/tr| 不同体系）与污染库排除（CON_ 前缀 accession＝非新肽段）
in_fasta <- vapply(proteins, function(pr) {
  acc <- acc_of(pr)
  length(acc) > 0 && any(acc %in% fa_acc)
}, logical(1))
con_acc <- vapply(proteins, function(pr) {
  acc <- acc_of(pr)
  length(acc) > 0 && any(startsWith(acc, "CON_"))
}, logical(1))
seq_in <- vapply(a_seqs, function(s) {
  nzchar(s) && !is.na(s) && annotated_seq_ul(s)
}, logical(1))
in_fasta <- in_fasta | con_acc | seq_in

# 谱图判定列（与标记版同口径：混合谱 = 错误谱图性质）
spec_verdict <- ifelse(as.integer(cls$chimeric) == 1L, "wrong_spectrum_mixture", "")
# diff_type（判错行 = 错误匹配清单；与标记版"diff_type 非空即判校正"同语义）
diff_type <- ifelse(as.integer(cls$vote_says_wrong) == 1L,
                    as.character(cls$error_reason), "")

out_dir2 <- file.path(result_dir, "1.unlabeled_result")
dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)

p1 <- cls[in_fasta,, drop = FALSE]
p1$a.Proteins <- proteins[in_fasta]
p1$diff_type <- diff_type[in_fasta]
write.table(p1, file.path(out_dir2, "1.annotation_correction_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8", na = "")

p2 <- cls[!in_fasta,, drop = FALSE]
p2$a.Proteins <- proteins[!in_fasta]
p2$denovo_sequence <- denovo_seq[!in_fasta]
p2$denovo_DB_match <- denovo_db[!in_fasta]
p2$spectrum_verdict <- spec_verdict[!in_fasta]

# ---- 可信新肽段六重存疑规则----
# Stage A 侧与 Stage B 侧同标准判定，全部通过者 credible_novel=TRUE；判定
# 优先级：判错 > 嵌合 > 峰集全同 > 已注释有力竞争 > 连续缺失≥2 > 低分卡值。
# 数据源＝classification_A 的逐谱候选（candidates 紧凑 "seq(mod)|…" 列）＋
# 第一步无标记谱（cfg_dir/unlabeled）。τ 与 PDF-A/主表同口径
# （unlabeled=0.0083，可经 config_step3 threshold_B 覆盖）。
aa <- read.delim(file.path(script_dir, "AA_mass_dataframe.txt"), sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table
TAU_UL <- c(n15 = 0.0083, c13 = 0.0083, silac = 0.0183, unlabeled = 0.0083)
tau_ul <- TAU_UL[["unlabeled"]]
tb_raw <- if (!is.null(cfg3$threshold_B) && length(cfg3$threshold_B) > 0)
  tail(cfg3$threshold_B, 1) else ""
if (nzchar(tb_raw)) {
  tv <- suppressWarnings(as.numeric(tb_raw))
  if (!is.na(tv)) tau_ul <- tv
}
# v1.0新肽段侧可能 0 行（序列级分侧生效后，参考库内
# 鉴定肽整侧归已注释侧）——列初始化必须按行数扩空，`<- ""` 在 0 行数据框上报
# "replacement has 1 row, data has 0"
p2$cand1_sequence <- rep("", nrow(p2)); p2$cand1_modification <- rep("", nrow(p2))
p2$n_frag <- rep("", nrow(p2)); p2$max_missing_run <- rep("", nrow(p2))
p2$closest_annotated_competitor <- rep("", nrow(p2))
p2$novel_doubt <- rep("", nrow(p2)); p2$novel_status <- rep("", nrow(p2))
p2$credible_novel <- rep("", nrow(p2))
raw_sm_ul <- function(s) {
  m <- regmatches(s, regexpr("[(].*[)]", s))
  c(sub("[(].*$", "", s), if (length(m) == 1L) substr(m, 2, nchar(m) - 1L) else "")
}
nov_idx <- which(!in_fasta)
cand_l <- strsplit(as.character(cls$candidates), "|", fixed = TRUE)
sc_l <- if ("candidate_PSM_score_by" %in% colnames(cls)) {
  lapply(strsplit(as.character(cls$candidate_PSM_score_by), "[|;]"),
         function(v) suppressWarnings(as.numeric(v)))
} else vector("list", nrow(cls))
# 只加载新肽候选行涉及的 raw 谱（only 参数，与主表同口径避免全量重读）
nov_raws <- unique(unname(vapply(as.character(cls$File_Name[nov_idx]),
                                 function(t) unname(title_raw_scan(t)["raw"]),
                                 character(1))))
nov_raws <- nov_raws[!is.na(nov_raws) & nzchar(nov_raws)]
lu_ul <- NULL
unlab_dir <- file.path(cfg_dir, "unlabeled")
if (length(nov_raws) > 0 && dir.exists(unlab_dir)) {
  lu_ul <- tryCatch(load_unlabeled_spectra(unlab_dir, only = nov_raws),
                    error = function(e) NULL)
}
for (k in nov_idx) {
  p2i <- match(k, nov_idx)
  fn <- as.character(cls$File_Name[k])
  cc <- cand_l[[k]]
  sm1 <- if (length(cc) > 0 && !is.na(cc[1])) raw_sm_ul(trimws(cc[1])) else c("", "")
  s1 <- sm1[1]; md1 <- sm1[2]
  p2$cand1_sequence[p2i] <- s1; p2$cand1_modification[p2i] <- md1
  # 谱图定位（用于离子证据规则）
  sp <- NULL
  if (!is.null(lu_ul) && nzchar(s1)) {
    raw <- unname(title_raw_scan(fn)["raw"])
    objs <- lu_ul$objs[[raw]]
    if (!is.null(objs)) {
      li <- match(fn, lu_ul$titles[[raw]])
      if (is.na(li)) {
        key <- title_scan_key(fn)
        ci <- which(startsWith(lu_ul$titles[[raw]], paste0(key, ".")))
        if (length(ci) > 0) li <- ci[1]
      }
      if (!is.na(li)) sp <- objs[[li]]@spectrum
    }
  }
  pos_hit <- integer(0)
  if (!is.null(sp) && nzchar(s1)) {
    ph <- tryCatch(matched_peak_set(s1, md1, sp, aa, aa_tab),
                   error = function(e) integer(0))
    pos_hit <- sort(unique(ph))
    if (length(pos_hit) > 0) p2$n_frag[p2i] <- as.character(length(pos_hit))
  }
  # 状态判定（优先级同标记版）
  doubt_k <- ""
  if (identical(as.integer(cls$vote_says_wrong[k]), 1L)) {
    er <- as.character(cls$error_reason[k])
    doubt_k <- if (is.na(er) || !nzchar(er)) "judged_wrong" else er
  } else if (identical(as.integer(cls$chimeric[k]), 1L)) {
    doubt_k <- "chimeric"
  } else if (nzchar(s1) && annotated_seq_ul(s1)) {
    # cand1 自身在参考库内——竞争结论即"有力已
    # 注释解释"，不能判可信新肽段；旧口径只在竞争者上做序列判定、cand1 从不
    # 被检，参考库肽段直通 credible_novel（实测约占 97%）
    doubt_k <- paste0("annotated_cand1:", s1)
  } else {
    # 规则 1：竞争候选与 cand1 的匹配峰集完全一致
    comps <- if (length(cc) > 1) unique(trimws(cc[-1])) else character(0)
    comps <- comps[!is.na(comps) & nzchar(comps)]
    if (!is.null(sp) && length(pos_hit) > 0 && length(comps) > 0) {
      for (s2 in comps) {
        sm2 <- raw_sm_ul(s2)
        B <- tryCatch(matched_peak_set(sm2[1], sm2[2], sp, aa, aa_tab),
                      error = function(e) integer(0))
        if (length(B) > 0 && setequal(pos_hit, B)) {
          doubt_k <- paste0("ions_identical:", sm2[1]); break
        }
      }
    }
    # 规则 2：打分最相近的已注释竞争 |ΔScore| ≤ τ → 无法排除已注释解释
    closest_k <- ""; closest_dk <- -Inf
    if (!nzchar(doubt_k) && length(cc) > 1 && length(sc_l[[k]]) >= length(cc)) {
      s1sc <- sc_l[[k]][1] / 100
      for (j in seq(2, length(cc))) {
        cs <- trimws(cc[j]); sm2 <- raw_sm_ul(cs)
        if (!nzchar(sm2[1]) || sm2[1] == s1) next
        if (!annotated_seq_ul(sm2[1])) next
        dv <- sc_l[[k]][j] / 100 - s1sc
        if (!is.finite(dv)) next
        if (dv > closest_dk) {
          closest_dk <- dv
          closest_k <- sprintf("%s|%s|%.6g", sm2[1], sm2[2], dv)
        }
      }
      if (nzchar(closest_k) && abs(closest_dk) <= tau_ul)
        doubt_k <- paste0("close_annotated:", closest_k)
    }
    if (nzchar(closest_k)) p2$closest_annotated_competitor[p2i] <- closest_k
    # 规则 3：连续 ≥2 个位点无子离子命中
    if (!nzchar(doubt_k) && !is.null(sp) && nzchar(s1)) {
      npos <- nchar(s1) - 1L
      if (npos >= 2 && length(pos_hit) < npos) {
        missing <- setdiff(seq_len(npos), pos_hit)
        ml <- 0L; cur <- 0L; prev <- -99L
        for (p in missing) {
          cur <- if (p == prev + 1L) cur + 1L else 1L
          if (cur > ml) ml <- cur
          prev <- p
        }
        p2$max_missing_run[p2i] <- as.character(ml)
        if (ml >= 2L) doubt_k <- sprintf("consecutive_missing:%d", ml)
      }
    }
    # 规则 4：低质量谱卡值（cand1 score < 0.2）
    if (!nzchar(doubt_k) && length(sc_l[[k]]) > 0 && !is.na(sc_l[[k]][1])) {
      sc_k <- sc_l[[k]][1] / 100
      if (sc_k < 0.2) doubt_k <- sprintf("low_score:%.4f", sc_k)
    }
  }
  if (nzchar(doubt_k)) {
    p2$novel_doubt[p2i] <- doubt_k
    p2$novel_status[p2i] <- sub(":.*$", "", doubt_k)
  } else {
    p2$credible_novel[p2i] <- "TRUE"
  }
}
write.table(p2, file.path(out_dir2, "2.novel_peptide_competition_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8", na = "")

n1_err <- sum(nzchar(p1$diff_type))
n2_cred <- sum(p2$credible_novel == "TRUE")
n2_doubt <- sum(nzchar(p2$novel_doubt))
sum2 <- file.path(out_dir2, "two_part_report_summary.txt")
con2 <- file(sum2, open = "wt", encoding = "UTF-8")
writeLines(c("part\trows\theadline",
             sprintf("1.annotation_correction_report\t%d\tannotation corrections (misidentifications): %d rows",
                     nrow(p1), n1_err),
             sprintf("2.novel_peptide_competition_report\t%d\tnovel competition rows: %d (credible novel peptides: %d; doubted %d)",
                     nrow(p2), nrow(p2), n2_cred, n2_doubt)), con2)
close(con2)

cat(sprintf("[两分法报告-无标记] Part1 注释校正：%d 行（判错 %d）-> %s\n",
            nrow(p1), n1_err, file.path(out_dir2, "1.annotation_correction_report.tsv")))
cat(sprintf("[two-part report (unlabeled)] Part1 annotation correction: %d rows (judged wrong %d) -> %s\n",
            nrow(p1), n1_err, file.path(out_dir2, "1.annotation_correction_report.tsv")))
cat(sprintf("[两分法报告-无标记] Part2 新肽段竞争：%d 行（可信新肽段 %d、存疑 %d）-> %s\n",
            nrow(p2), n2_cred, n2_doubt,
            file.path(out_dir2, "2.novel_peptide_competition_report.tsv")))
cat(sprintf("[two-part report (unlabeled)] Part2 novel competition: %d rows (credible %d, doubted %d) -> %s\n",
            nrow(p2), n2_cred, n2_doubt,
            file.path(out_dir2, "2.novel_peptide_competition_report.tsv")))
