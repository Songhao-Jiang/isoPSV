# =============================================================================
# spectrum_pdf_export.R - 错误匹配与存疑谱图 PDF 导出（流程结尾自动运行）
# spectrum_pdf_export.R - error-match and ambiguous-match spectrum PDF export (runs
# automatically at the end of the pipeline)
#
# 用途：
# 在 user_run.R 三步流程（阶段 A＋阶段 B＋冲突报告）完成后自动执行本模块，
# 用 Fig.7 定稿判定模型（五指标加权得分＋阈值 τ＋预测救援 θ）对成对报告
# 逐谱判定，并把两类谱渲染成 PDF：
#   ① error_match_spectra.pdf ＝ 判需校正（call_v2）的全部谱图；
#   ② ambiguous_matches_spectra.pdf ＝ 需要人工判断的难区分带谱（ions_small_diff，
#      |τ−ΔScore|≤0.015 且 ΔScore≠0）。
# 每谱一页四面板（候选 1 谱＋其 ppm 图、胜选候选谱＋其 ppm 图），含信息行
# （谱图号/打分/结论）、离子梯（价态前缀、黑色 1 磅切口竖条、b 下标在下/
# y 下标在上）、Base peak 与 Mass error 右上角注释。
#
# Purpose:
# This module runs automatically after the three-step pipeline (stage A +
# stage B + conflict report) in user_run.R. It applies the finalized Fig.7
# judgment model (five-metric weighted score + threshold tau + prediction
# rescue theta) to every spectrum of the paired report and renders two PDFs:
#   ① error_match_spectra.pdf = all spectra judged "need correction" (call_v2);
#   ② ambiguous_matches_spectra.pdf = spectra in the hard-to-judge band (ions_small_diff,
#     |tau - DeltaScore| <= 0.015 and DeltaScore != 0).
# Each spectrum is one page with four panels (candidate-1 spectrum + its ppm
# panel, winner spectrum + its ppm panel), the info lines (spectrum id / scores
# / conclusion), the ion ladder (charge prefix, black 1-pt cleavage ticks, b
# subscripts below / y subscripts above), and the Base peak + Mass error
# annotations at the top right.
#
# 输入（环境变量，由 user_run.R 设置）：
#   CONFIG3  ＝ 第三步参数文件路径（由此定位结果目录与标记类型）
#   ALL_MGF  ＝ 原始 MGF 文件夹（无标记＋标记谱）
#   PDF_LIMIT＝ 可选，最多渲染张数（默认全部；用于快速自检）
#
# Inputs (environment variables, set by user_run.R):
#   CONFIG3  = path of the step-3 parameter file (locates the result folder
#              and the label type)
#   ALL_MGF  = the original MGF folder (unlabeled + labeled spectra)
#   PDF_LIMIT = optional, maximum pages to render (default: all; for smoke tests)
#
# 输出：
#   <结果目录>/1.unlabeled_result/error_match_spectra.pdf
#   <结果目录>/1.unlabeled_result/ambiguous_matches_spectra.pdf
#   <结果目录>/1.unlabeled_result/novel_peptide_spectra.pdf（可信新肽段，
#     2026-09-26 起——与标记与否无关，两分法报告 credible_novel=TRUE 行；
#     上=新肽段本身、下=最强已注释竞争者）
#
# Outputs:
#   <result folder>/1.unlabeled_result/error_match_spectra.pdf
#   <result folder>/1.unlabeled_result/ambiguous_matches_spectra.pdf
#   <result folder>/1.unlabeled_result/novel_peptide_spectra.pdf (credible
#     novel peptides, from 2026-09-26 - independent of labeling; top = the
#     novel peptide, bottom = its strongest annotated competitor)
# =============================================================================
args0 <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "",
  args0[grepl("--file=", args0)][1]), winslash = "/"))
# ---- 库环境：本模块与全流程代码同目录，直接加载同目录库函数 ----
# ---- library environment: this module lives in the pipeline code folder and
# loads the sibling library files directly
old_wd <- setwd(script_dir)
source("compat.R")
source("ReadMGF_script.R")
suppressMessages(source("acs_core.R")); suppressMessages(source("paired_params.R"))
source("chimeric_check.R")
setwd(old_wd)
aa <- read.delim(file.path(script_dir, "AA_mass_dataframe.txt"), sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# ---- 读取环境与参数文件 ----
# ---- read the environment and the parameter files
cfg3_path <- Sys.getenv("CONFIG3")
if (!nzchar(cfg3_path)) stop("CONFIG3 not set (set by user_run.R)")
cfg3 <- read_param_file(cfg3_path)
cfg_dir <- cfg3$config_dir
if (is.null(cfg_dir) || length(cfg_dir) == 0 || !nzchar(cfg_dir[1]))   # 2026-09-26 审查 P0-3/P3-6：config_dir 缺失/为空时给可读报错
  stop("config_step3 的 config_dir 缺失或为空，请检查配置（或删除 result_* 后重跑）
config_dir missing/empty in config_step3")
# 结果文件夹（RESULT_DIR 由总程序注入；独立运行按时间戳定位/新建）
# result folder (RESULT_DIR injected by the master; standalone runs locate/create by timestamp)
result_dir <- get_result_dir(cfg_dir)
base_cfg <- read_param_file(resolve_cfg_paths(cfg3, "base_config")[1])
label_type <- if (!is.null(base_cfg$label_type) && length(base_cfg$label_type) > 0)
  tolower(tail(base_cfg$label_type, 1)) else "n15"
# 无标记谱用第一步提取的逐 raw 文件（cfg_dir/unlabeled；两样本类型均有）
# the unlabeled spectra come from the step-1 per-raw files (cfg_dir/unlabeled,
# present for both sample types)
all_mgf_dir <- file.path(cfg_dir, "unlabeled")
pdf_limit <- suppressWarnings(as.numeric(Sys.getenv("PDF_LIMIT", "Inf")))
if (is.na(pdf_limit)) pdf_limit <- Inf
# PDF 拆分（2026-09-26 用户口径）：新肽段/错误谱 PDF 页数太大、单文件加载过重——
# 按 raw 分组输出，组内超过 pdf_pages_per_file（默认 1000 页）再切 part；
# 谱图 PDF 索引 spectrum_pdf_index.tsv 记录每谱所属类别/文件/页码
ppf_raw <- if (!is.null(base_cfg$pdf_pages_per_file) &&
               length(base_cfg$pdf_pages_per_file) > 0)
  suppressWarnings(as.numeric(tail(base_cfg$pdf_pages_per_file, 1))) else NA_real_
ppf <- if (is.finite(ppf_raw) && ppf_raw >= 1) as.integer(ppf_raw) else 1000L
IDX <- list()

# ---- 定稿判定模型参数（Fig.7 b–e；不同标记组各自阈值，权重统一）----
# ---- finalized judgment-model parameters (Fig.7 b-e; per-label thresholds,
# shared weights)
W6 <- c(ul_by = 0.6200, lb_by = 0.0340, ul_all = 0.1186, lb_all = 0.1700,
        predict = 0, pcos = 0.0573)
# P1（2026-09-26，deepseek R29 实弹）：TAU_MAP 缺 unlabeled——[[ 对原子向量
# 取不存在的名字是报错而非 NULL，旧兜底恒不生效，无标记样本 PDF 导出必崩；
# 与 master_table_export.R 同口径补 unlabeled，并改防御式查找
TAU_MAP <- c(n15 = 0.0083, c13 = 0.0083, silac = 0.0183, unlabeled = 0.0083)
THETA <- 0.30
tau_v <- label_tau_of(parse_label_spec(label_type))
# 分组权重覆盖（用户口径 2026-09-21）：config_step3 可选键 weights_B（六指标
# 逗号分隔，顺序 ul_by,lb_by,ul_all,lb_all,predict,pcos）与 threshold_B 优先；
# 旧配置无这两键时行为与原版一致（与 master_table_export.R 同一约定）
# per-label override (user ruling 2026-09-21): optional config keys weights_B /
# threshold_B take precedence; absent keys reproduce the unified behavior
wb_raw <- if (!is.null(cfg3$weights_B) && length(cfg3$weights_B) > 0)
  tail(cfg3$weights_B, 1) else ""
if (nzchar(wb_raw)) {
  wb <- as.numeric(strsplit(wb_raw, ",")[[1]])
  if (length(wb) == 8L)
    stop("weights_B 为八指标（含 ul_cov/lb_cov），谱图 PDF 导出尚未支持八指标打分；请设 PDF_ON=0 或改用六指标权重
weights_B carries eight metrics; the spectrum-PDF export does not support the 8-metric scoring yet (set PDF_ON=0)")
  stopifnot(length(wb) == 6L, !any(is.na(wb)))
  W6 <- setNames(wb, names(W6))
  cat("[分组权重] weights_B 覆盖：", paste(sprintf("%.4f", W6), collapse = ","), "
")
}
tb_raw <- if (!is.null(cfg3$threshold_B) && length(cfg3$threshold_B) > 0)
  tail(cfg3$threshold_B, 1) else ""
if (nzchar(tb_raw)) {
  tau_v <- as.numeric(tb_raw)
  cat("[分组权重] threshold_B 覆盖：", tau_v, "
")
}

# ---- 读入阶段 A 报告（无标记侧评估；判定用 decision_A）----
# ---- read the stage-A report (unlabeled-side evaluation; the judgment uses
# decision_A)
sta_dir <- file.path(result_dir, "1.unlabeled_result")
repA <- read.delim(file.path(sta_dir, "report_A.tsv"), sep = "\t",
                   check.names = FALSE, na.strings = character(0), quote = "")
cval <- function(v, dflt = "-") {                              # 缺失行安全取值
  x <- as.character(v)
  x[is.na(x)] <- dflt
  x
}
# ---- 候选与打分解析：candidates 为 "seq(mod)|seq(mod)" 紧凑列表；打分取
# candidate_PSM_score_by（PSM 类百分数，除以 100 与 Fig.7 口径一致）----
parse_col <- function(x) lapply(strsplit(as.character(x), "[|;]"), function(v)
  suppressWarnings(as.numeric(v)))
cand_raw <- strsplit(as.character(repA$candidates), "|", fixed = TRUE)
n <- nrow(repA)
scores <- parse_col(repA$candidate_PSM_score_by)
scores <- lapply(scores, function(v) v / 100)
rev_idx <- integer(n)
delta <- rep(NA_real_, n)
for (i in seq_len(n)) {
  s <- scores[[i]]
  if (length(s) >= 2 && !all(is.na(s)) && !all(is.na(s[-1]))) {
    k2 <- which.max(s[-1]) + 1L
    rev_idx[i] <- k2
    delta[i] <- s[k2] - s[1]
  }
}
predL <- parse_col(repA$candidate_predict_cosine)
dpred <- vapply(seq_len(n), function(i) {
  p <- predL[[i]]
  if (length(p) >= 2 && rev_idx[i] >= 2 &&
      !is.na(p[1]) && !is.na(p[rev_idx[i]])) p[rev_idx[i]] - p[1] else NA_real_
}, numeric(1))
# 判错＝阶段 A 判定 candidate_misidentification；存疑＝打分差近阈值或恰为 0
# error = the stage-A decision "candidate_misidentification"; ambiguous =
# the by-score delta near the threshold or exactly zero
judged <- !is.na(repA$decision_A) & repA$decision_A == "candidate_misidentification"
rescue <- rep(FALSE, n)
ions_identical <- !is.na(delta) & delta == 0
ions_small_diff <- !is.na(delta) & delta != 0 & abs(tau_v - delta) <= 0.015

# ---- 两类谱集合 ----
# ---- the two spectrum sets
err_idx <- which(judged)
sus_idx <- which(ions_small_diff)
if (is.finite(pdf_limit)) {
  err_idx <- head(err_idx, pdf_limit)
  sus_idx <- head(sus_idx, pdf_limit)
}
cat(sprintf("[谱图 PDF] 错误匹配 %d 张；存疑 %d 张（τ=%.4f）\n", length(err_idx),
            length(sus_idx), tau_v))
cat(sprintf("[spectrum PDF] error-match %d; ambiguous-match %d (tau=%.4f)\n", length(err_idx),
            length(sus_idx), tau_v))

# ---- 谱图读取（库内三档并行读取器）----
# ---- spectrum loading (the library three-tier parallel reader)
old2 <- setwd(script_dir)
lu <- load_unlabeled_spectra(all_mgf_dir)
setwd(old2)
if (is.null(lu)) stop("no MGF spectra loaded from ALL_MGF")

# ---- 绘制引擎（与 Fig.7 谱图定稿版式一致）----
# ---- the drawing engine (identical to the finalized Fig.7 spectrum layout)
ION_COL <- c(a = "#004000", b = "#008000", c = "#CA1D52", x = "#650E29",
             y = "#CA1D52", z = "#943AA4", M = "#808080", I = "#FF8C00", O = "#9ACD32")
ION_PRI <- c(b = 1, y = 2, a = 3, c = 4, x = 5, z = 6, M = 7, I = 8, O = 9)
CO_M <- 27.994915; NH3_M <- 17.026549; H2O_M <- 18.010565; H_M <- 1.007276

# 全部理论离子（a/b/c/x/y/z 1+/2+ 与 M 系）＋20 ppm 窗最近峰匹配
# full theoretical ion ladder (a/b/c/x/y/z 1+/2+ and the M series) matched to
# the nearest peak inside a 20 ppm window
all_ions <- function(seq, mod, mz, ints, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) || length(bip) == 1 ||
        is.null(bip$ions_table)) next
    t <- bip$ions_table; num <- t$Number
    theo <- data.frame(
      mz = c(t$b_1_charge, t$b_1_charge - CO_M, t$b_1_charge + NH3_M,
             t$b_2_charge, t$b_2_charge - CO_M / 2, t$b_2_charge + NH3_M / 2,
             t$y_1_charge, t$y_1_charge + CO_M, t$y_1_charge - NH3_M,
             t$y_2_charge, t$y_2_charge + CO_M / 2, t$y_2_charge - NH3_M / 2),
      lab = c(paste0("b", num, "+"), paste0("a", num, "+"), paste0("c", num, "+"),
              paste0("b", num, "++"), paste0("a", num, "++"), paste0("c", num, "++"),
              paste0("y", num, "+"), paste0("x", num, "+"), paste0("z", num, "+"),
              paste0("y", num, "++"), paste0("x", num, "++"), paste0("z", num, "++")),
      ion = c(rep(c("b", "a", "c"), each = length(num)),
              rep(c("b", "a", "c"), each = length(num)),
              rep(c("y", "x", "z"), each = length(num)),
              rep(c("y", "x", "z"), each = length(num))),
      num = c(rep(num, 3), rep(num, 3), rep(num, 3), rep(num, 3)),
      chg = c(rep(1L, 3 * length(num)), rep(2L, 3 * length(num)),
              rep(1L, 3 * length(num)), rep(2L, 3 * length(num))),
      stringsAsFactors = FALSE)
    if (is.list(bip) && !is.null(bip$pep_mass) && is.finite(bip$pep_mass)) {
      pm <- bip$pep_mass
      theo <- rbind(theo, data.frame(
        mz = c(pm + H_M, pm + H_M - H2O_M, pm + H_M - NH3_M,
               (pm + 2 * H_M) / 2, (pm + 2 * H_M - H2O_M) / 2, (pm + 2 * H_M - NH3_M) / 2),
        lab = c("M+", "M-H2O+", "M-NH3+", "M++", "M-H2O++", "M-NH3++"),
        ion = "M", num = NA_integer_, chg = c(1L, 1L, 1L, 2L, 2L, 2L),
        stringsAsFactors = FALSE))
    }
    theo <- theo[is.finite(theo$mz), ]
    lo <- theo$mz * (1 - tol_ppm * 1e-6); hi <- theo$mz * (1 + tol_ppm * 1e-6)
    i_lo <- findInterval(lo, mz); i_hi <- findInterval(hi, mz)
    peak_i <- rep(NA_integer_, nrow(theo))
    for (j in seq_len(nrow(theo))) {
      if (i_hi[j] <= i_lo[j]) next
      cand <- (i_lo[j] + 1L):i_hi[j]
      if (length(cand)) peak_i[j] <- cand[which.max(ints[cand])]
    }
    theo$peak <- peak_i
    return(theo)
  }
  NULL
}

# 序列离子梯：黑色 1 磅切口竖条（仅限有 b/y 命中的切口），b 下标在下、
# y 下标在上，价态前缀在梯行左端，修饰残基橙色加粗
# sequence ion ladder: black 1-pt cleavage ticks (only at b/y-matched sites),
# b subscripts below and y subscripts above, charge prefix at the ladder's
# left end, modified residues bold orange
draw_ladder <- function(sq, ions, mod = "", chg_v = 0L) {
  nch <- nchar(sq)
  if (nch < 2L) return(invisible())
  u <- par("usr")
  x0 <- u[1] + 0.015 * (u[2] - u[1]); x1 <- u[2] - 0.015 * (u[2] - u[1])
  u_per_in <- abs(diff(grconvertX(c(0, 1), from = "inches", to = "user")))
  slot_u <- min(0.16 * u_per_in, (x1 - x0) / nch)
  xc0 <- (x0 + x1) / 2 - nch * slot_u / 2
  xs <- xc0 + (seq_len(nch) - 0.5) * slot_u
  gap <- slot_u; ymid <- 1.30; v_h <- 0.045
  site_by <- b_sites <- y_sites <- integer(0)
  if (!is.null(ions)) {
    hit <- ions[!is.na(ions$peak), ]
    if (nrow(hit)) {
      hb <- hit[hit$ion == "b", ]
      if (nrow(hb)) b_sites <- sort(unique(as.integer(hb$num)))
      hy <- hit[hit$ion == "y", ]
      if (nrow(hy)) y_sites <- sort(unique(nch - as.integer(hy$num)))
      site_by <- sort(unique(c(b_sites, y_sites)))
    }
  }
  op <- par(xpd = TRUE); on.exit(par(op), add = TRUE)
  slot_in <- slot_u / u_per_in
  cexs <- min(8 / 12, max(6 / 12, 0.9 * slot_in / strwidth("M", cex = 1, units = "inches")))
  cexn <- min(6 / 12, max(5 / 12, 0.9 * slot_in / strwidth("b88", cex = 1, units = "inches")))
  if (length(site_by)) {
    vx <- (xs[site_by] + xs[site_by + 1L]) / 2
    segments(vx, ymid - v_h, vx, ymid + v_h, col = "black", lwd = 1.33)
  }
  if (length(b_sites)) {
    vxb <- (xs[b_sites] + xs[b_sites + 1L]) / 2
    segments(vxb, ymid - v_h, vxb - gap * 0.5, ymid - v_h, col = "black", lwd = 1.33)
    text(xs[b_sites], ymid - 0.085,
         as.expression(lapply(b_sites, function(k) bquote(b[.(k)]))),
         cex = cexn, col = ION_COL[["b"]], font = 2)
  }
  if (length(y_sites)) {
    vxy <- (xs[y_sites] + xs[y_sites + 1L]) / 2
    segments(vxy, ymid + v_h, vxy + gap * 0.5, ymid + v_h, col = "black", lwd = 1.33)
    text(xs[y_sites + 1L], ymid + 0.085,
         as.expression(lapply(y_sites, function(k) bquote(y[.(nch - k)]))),
         cex = cexn, col = ION_COL[["y"]], font = 2)
  }
  if (is.finite(chg_v) && chg_v > 0L)
    text(xc0 - gap * 1.5, ymid, sprintf("%d+ ", chg_v), cex = cexs, font = 2,
         adj = 1, col = "#CA1D52")
  lab_col <- rep("black", nch); lab_font <- rep(1, nch)
  if (!is.null(mod) && nzchar(mod)) {
    pm <- suppressWarnings(as.integer(sub(",.*$", "", strsplit(mod, ";", fixed = TRUE)[[1]])))
    pm <- pm[!is.na(pm) & pm >= 1L & pm <= nch]
    lab_col[pm] <- ION_COL[["I"]]; lab_font[pm] <- 2
  }
  text(xs, ymid, strsplit(sq, "")[[1]], cex = cexs, font = lab_font, col = lab_col)
}

# 单谱一页：两候选谱＋各自 ppm 面板＋信息行（v7 定稿版式）
# one page per spectrum: two candidate spectra + their ppm panels + info lines
# (the finalized v7 layout)
draw_page <- function(sp_obj, top, bot, lines5, chg_v, ann) {
  layout(matrix(c(1, 2, 3, 4), 4, 1), heights = c(1.0, 0.35, 1.0, 0.35))
  par(oma = c(4, 3, 9, 3))
  par(cex = 1)
  mz <- as.numeric(sp_obj[[1]]); ints <- as.numeric(sp_obj[[2]])
  o <- order(mz, na.last = NA); mz <- mz[o]; ints <- ints[o]
  ints_n <- ints / max(ints[is.finite(ints) & ints > 0])
  xmax <- max(mz, na.rm = TRUE) * 1.03
  for (pi in 1:2) {
    pn <- list(top, bot)[[pi]]
    par(mar = c(0.5, 4, 4, 1))
    plot(0, 0, type = "n", xlim = c(0, xmax), ylim = c(0, 1.15),
         xlab = "", ylab = "", cex.lab = 0.83, cex.axis = 0.83,
         yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
    title(ylab = "Relative intensity (%)", cex.lab = 0.83, line = 2.6)
    par(mgp = c(2.0, 0.8, 0))
    axis(2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 100, 20)),
         cex.axis = 0.83, las = 1)
    text(xmax * 0.995, 1.10, ann[[1]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
    text(xmax * 0.995, 1.00, ann[[2]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
    segments(mz, 0, mz, ints_n, lwd = 0.3, col = "grey60")
    th <- pn$ions
    if (!is.null(th)) {
      hit <- th[!is.na(th$peak), ]
      if (nrow(hit)) {
        o <- order(ION_PRI[hit$ion], hit$num, hit$chg, na.last = TRUE)
        hit <- hit[o, ]
        spl <- split(hit, hit$peak)
        pk_i <- as.integer(names(spl))
        labs <- vapply(spl, function(d) paste(d$lab, collapse = ""), "")
        cols <- vapply(spl, function(d) ION_COL[[d$ion[1]]], "")
        pmz <- mz[pk_i]; pin <- ints_n[pk_i]
        segments(pmz, 0, pmz, pin, lwd = 0.8, col = cols)
        text(pmz, pin + 0.03, labs, cex = 8 / 12, srt = 90, col = cols, adj = 0)
      }
    }
    draw_ladder(pn$seq, th, pn$mod, chg_v)
    mtext(pn$head, side = 3, line = 3, cex = 0.83, adj = 0)
    # ppm 面板：横轴与谱面板一致，纵轴严格 ±20（匹配窗），命中离子色点
    # the ppm panel: same x axis as the spectrum, strict +/-20 ppm (the match
    # window), one colored point per matched ion
    par(mar = c(2.6, 4, 0.3, 1), xpd = FALSE)
    plot(NA, xlim = c(0, xmax), ylim = c(-20, 20),
         xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i", yaxs = "i")
    title(xlab = expression(italic("m/z")), cex.lab = 0.83, line = 1.8)
    title(ylab = "ppm", cex.lab = 0.83, line = 2.6)
    abline(h = 0, col = "#999999", lwd = 0.6)
    abline(h = c(-20, 20), col = "#BBBBBB", lty = 3, lwd = 0.6)
    if (!is.null(th)) {
      hm <- th[!is.na(th$peak), ]
      if (nrow(hm)) {
        px <- mz[hm$peak]
        py <- (px - hm$mz) / hm$mz * 1e6
        points(px, py, col = ION_COL[hm$ion], pch = 16, cex = 0.5)
      }
    }
    atk <- pretty(c(0, xmax), n = 10)
    atk <- atk[atk >= 0 & atk <= xmax]
    par(mgp = c(2.0, 0.8, 0))
    axis(1, at = atk, cex.axis = 0.83)
    axis(2, at = c(-20, 0, 20), las = 1, cex.axis = 0.83)
    box()
  }
  title(lines5[1], outer = TRUE, line = 5.2, cex.main = 1.5)
  title(lines5[2], outer = TRUE, line = 3.4, cex.main = 0.83)
  title(lines5[3], outer = TRUE, line = 2.4, cex.main = 0.83)
  title(lines5[4], outer = TRUE, line = 1.4, cex.main = 0.83)
  title(lines5[5], outer = TRUE, line = 0.4, cex.main = 0.83)
}

# ---- 候选拆解与蛋白侧翼 ----
# ---- candidate parsing and protein flanking
raw_sm <- function(i, k) {
  cc <- cand_raw[[i]]
  if (length(cc) < k) return(c("", ""))
  s <- trimws(cc[k])
  m <- regmatches(s, regexpr("[(].*[)]", s))
  c(sub("[(].*$", "", s), if (length(m) == 1L) substr(m, 2, nchar(m) - 1L) else "")
}
# 蛋白与侧翼的权威来源＝主表（cand1_*/corrected_* 覆盖 100%；阶段 A 报告本身
# 无蛋白列）。主表定位格式 "K[280]C"（前残基[0 起位置]后残基；"_" 为 Excel
# 防公式占位，还原为 "-"）；无标记样本没有主表时退回 "-"
# the authoritative protein/flanking source is the master table (cand1_*/
# corrected_* cover 100%; the stage-A report itself has no protein column).
# The master location format is "K[280]C" (flank[0-based pos]flank; "_" is the
# Excel-safe placeholder restored to "-"); unlabeled samples have no master
# table and fall back to "-"
mt <- NULL
mt_path <- file.path(result_dir, "2.labeled_result", "master_table.tsv")
if (file.exists(mt_path)) {
  mt <- read.delim(mt_path, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = character(0), quote = "")
}
mcol <- function(col) {
  if (is.null(mt) || !col %in% colnames(mt)) return(rep("-", n))
  v <- as.character(mt[[col]][match(as.character(repA$File_Name),
                                    as.character(mt$File_Name))])
  v[is.na(v) | v %in% c("", "NA")] <- "-"
  v
}
p1_prot <- mcol("cand1_protein");  p1_loc <- mcol("cand1_location")
p2_prot <- mcol("corrected_protein"); p2_loc <- mcol("corrected_location")
fmt_head <- function(sq, md, acc, pos) {
  pre <- "-"; post <- "-"
  if (nzchar(pos) && pos != "-") {
    p2 <- gsub("_", "-", pos, fixed = TRUE)
    if (grepl("\\[[0-9]+\\]", p2)) {
      pre <- sub("\\[.*$", "", p2)
      post <- sub("^.*\\]", "", p2)
      if (!nzchar(pre) || pre == "-") pre <- "-"
      if (!nzchar(post) || post == "-") post <- "-"
    }
  }
  sprintf("| %s | (%s)%s(%s) | %s |",
          if (nzchar(acc) && acc != "-") acc else "-", pre, sq, post,
          if (nzchar(md)) md else "-")
}
fmt_num <- function(v, dig) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)
}

# ---- 渲染批次：rows_idx 索引的谱写入指定 PDF；返回实际画出的谱名（顺序）----
# ---- render batch: write the spectra of rows_idx into the given PDF; return
# the spectra actually drawn (in order)
render_batch <- function(rows_idx, file_tag, pdf_name = NULL) {
  if (is.null(pdf_name)) pdf_name <- paste0(file_tag, ".pdf")
  pdf_f <- file.path(sta_dir, pdf_name)
  pdf(pdf_f, width = 8.27, height = 11.69)
  drawn <- character(0)
  for (r in rows_idx) {
    fn <- as.character(repA$File_Name[r])
    raw <- unname(title_raw_scan(fn)["raw"])
    objs <- lu$objs[[raw]]
    if (is.null(objs)) next
    li <- match(fn, lu$titles[[raw]])
    if (is.na(li)) {
      key <- title_scan_key(fn)
      ci <- which(startsWith(lu$titles[[raw]], paste0(key, ".")))
      if (length(ci) == 0) next
      li <- ci[1]
    }
    sp_obj <- objs[[li]]@spectrum
    mz0 <- as.numeric(sp_obj[[1]]); ints0 <- as.numeric(sp_obj[[2]])
    # 价态与注释：Base peak＝谱峰最大强度；Mass error＝Exp.MH+（PEPMASS×z 还原）
    # − Calc.MH+（by_ions_position 的 pep_mass 已含质子），Da 与 ppm 双口径
    # charge and annotations: Base peak = the largest intensity; Mass error =
    # Exp.MH+ (PEPMASS scaled back by z) - Calc.MH+ (pep_mass already includes
    # the proton), reported in Da and ppm
    chg_v <- tryCatch(as.integer(objs[[li]]@charge), error = function(e) NA_integer_)
    if (is.na(chg_v) || chg_v < 1L) chg_v <- 0L
    exp_mh <- tryCatch(as.numeric(objs[[li]]@pepMass), error = function(e) NA_real_)
    if (is.finite(exp_mh) && chg_v > 1L)
      exp_mh <- exp_mh * chg_v - (chg_v - 1L) * 1.007276
    s1 <- raw_sm(r, 1); k2 <- rev_idx[r]
    s2 <- if (k2 >= 2) raw_sm(r, k2) else s1
    calc_mh <- NA_real_
    btmp <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(s1[1], s1[2], "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (!inherits(btmp, "try-error") && is.list(bip) &&
        !is.null(bip$pep_mass) && is.finite(bip$pep_mass))
      calc_mh <- bip$pep_mass
    base_pk <- max(ints0[is.finite(ints0) & ints0 > 0])
    bp_lab <- sub("E\\+0?", "E", sprintf("%.2E", base_pk))
    ann <- c(sprintf("Base peak = %s", bp_lab),
             if (is.finite(exp_mh) && is.finite(calc_mh) && calc_mh > 0)
               sprintf("Mass error = %+.4f Da / %+.2f ppm",
                       exp_mh - calc_mh, (exp_mh - calc_mh) / calc_mh * 1e6)
             else "Mass error = -")
    top <- list(head = fmt_head(s1[1], s1[2], p1_prot[r], p1_loc[r]),
                seq = s1[1], mod = s1[2],
                ions = all_ions(s1[1], s1[2], mz0, ints0))
    bot <- list(head = fmt_head(s2[1], s2[2], p2_prot[r], p2_loc[r]),
                seq = s2[1], mod = s2[2],
                ions = all_ions(s2[1], s2[2], mz0, ints0))
    # 信息行五行：Type（本页判定口径的标签）、谱图号、打分、修饰差异、判定
    # five info lines: Type (labeled by THIS page's judgment), spectrum id,
    # scores, modification difference, judgment
    type_txt <- if (judged[r]) "candidate misidentification"
      else "ambiguous: score near threshold"
    lines5 <- c(
      sprintf("Type: %s", type_txt),
      sprintf("Spectrum: %s", basename(fn)),
      sprintf("Scores: candidate 1 = %.4f | candidate %d score = %s | Delta score = %s | Tau = %.4f",
              scores[[r]][1], k2,
              if (k2 >= 2 && length(scores[[r]]) >= k2) fmt_num(scores[[r]][k2], 4L) else "-",
              fmt_num(delta[r], 6L), tau_v),
      "Modification change: -",
      sprintf("Amino acid change: - | Judgment: %s%s",
              if (judged[r]) "need correction" else "keep original",
              if (rescue[r]) " + prediction rescue" else ""))
    # 自查整改（2026-09-26）：单谱渲染失败只跳过该谱并记录，不中断整个批次
    ok <- tryCatch({
      draw_page(sp_obj, top, bot, lines5, chg_v, ann)
      TRUE
    }, error = function(e) {
      cat("[谱图 PDF] 跳过渲染失败谱：", fn, "（", conditionMessage(e), "）\n")
      FALSE
    })
    if (ok) drawn <- c(drawn, fn)
  }
  dev.off()
  # 自查整改（2026-09-26）：批次内无谱可画时删掉刚建的空文件，不留孤儿 PDF
  if (length(drawn) == 0 && file.exists(pdf_f)) file.remove(pdf_f)
  drawn
}
# 自查整改（2026-09-26）：raw 名原样进文件名不安全（Thermo 标题可含盘符/冒号
# 等非法字符），清洗限长并去重
raw_token <- function(rw) {
  s <- substr(rw, 1, 80)
  s <- gsub("[^A-Za-z0-9_.-]", "_", s)
  if (!nzchar(s)) s <- "_"
  s
}
# ---- 按 raw 分组 + 组内超阈值切片的驱动（2026-09-26 拆分口径）----
render_grouped <- function(rows_idx, file_tag, get_fn) {
  if (length(rows_idx) == 0) return(invisible(character(0)))
  fns <- get_fn(rows_idx)
  raws <- unname(vapply(fns, function(t) unname(title_raw_scan(t)["raw"]),
                        character(1)))
  used <- character(0)
  for (rw in unique(raws)) {
    sel <- rows_idx[raws == rw]
    n <- length(sel)
    nparts <- ceiling(n / ppf)
    tok <- raw_token(rw)
    while (tok %in% used) tok <- paste0(tok, "_2")
    used <- c(used, tok)
    for (prt in seq_len(nparts)) {
      chunk <- sel[seq.int(1 + (prt - 1) * ppf, min(n, prt * ppf))]
      pdf_name <- if (nparts == 1L) paste0(file_tag, "__", tok, ".pdf")
                  else sprintf("%s__%s__part%d.pdf", file_tag, tok, prt)
      drawn <- render_batch(chunk, file_tag, pdf_name)
      for (j in seq_along(drawn))
        IDX[[length(IDX) + 1L]] <<- c(drawn[j], file_tag, pdf_name,
                                      as.character(j))
    }
  }
  invisible(character(0))
}

# ---- 第三类：可信新肽段谱图（novel_peptide_spectra.pdf，2026-09-26 用户口径：
# 可信新肽段判定与标记与否无关——无标记两分法报告同标准计算 credible_novel）----
# 每页复用四面板管线：上=新肽段本身（cand1_sequence 口径），下=最强已注释
# 竞争者（closest_annotated_competitor）；竞争者缺失时下层面板留空。
nov <- NULL
nov_path <- file.path(sta_dir, "2.novel_peptide_competition_report.tsv")
if (file.exists(nov_path)) {
  nov <- read.delim(nov_path, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = character(0), quote = "")
}
nov_idx <- integer(0)
if (!is.null(nov) && nrow(nov) > 0 && "credible_novel" %in% colnames(nov)) {
  nov_idx <- which(tolower(as.character(nov$credible_novel)) %in% c("true", "1"))
  if (is.finite(pdf_limit)) nov_idx <- head(nov_idx, pdf_limit)
}
render_novel_batch <- function(rows_idx, pdf_name = NULL) {
  gcol <- function(col) {
    if (col %in% colnames(nov)) as.character(nov[[col]]) else rep("", nrow(nov))
  }
  fns <- gcol("File_Name"); seqs <- gcol("cand1_sequence")
  mods <- gcol("cand1_modification"); comps <- gcol("closest_annotated_competitor")
  prots <- gcol("a.Proteins")
  if (is.null(pdf_name)) pdf_name <- "novel_peptide_spectra.pdf"
  pdf_f <- file.path(sta_dir, pdf_name)
  pdf(pdf_f, width = 8.27, height = 11.69)
  drawn <- character(0)
  for (k in rows_idx) {
    fn <- fns[k]
    if (!nzchar(fn)) next
    raw <- unname(title_raw_scan(fn)["raw"])
    objs <- lu$objs[[raw]]
    if (is.null(objs)) next
    li <- match(fn, lu$titles[[raw]])
    if (is.na(li)) {
      key <- title_scan_key(fn)
      ci <- which(startsWith(lu$titles[[raw]], paste0(key, ".")))
      if (length(ci) == 0) next
      li <- ci[1]
    }
    sp_obj <- objs[[li]]@spectrum
    mz0 <- as.numeric(sp_obj[[1]]); ints0 <- as.numeric(sp_obj[[2]])
    chg_v <- tryCatch(as.integer(objs[[li]]@charge), error = function(e) NA_integer_)
    if (is.na(chg_v) || chg_v < 1L) chg_v <- 0L
    exp_mh <- tryCatch(as.numeric(objs[[li]]@pepMass), error = function(e) NA_real_)
    if (is.finite(exp_mh) && chg_v > 1L)
      exp_mh <- exp_mh * chg_v - (chg_v - 1L) * 1.007276
    seq1 <- seqs[k]; mod1 <- mods[k]
    # 竞争者列 "seq|mod|delta" → 取序列
    comp_raw <- comps[k]
    seq2 <- if (nzchar(comp_raw))
      strsplit(comp_raw, "|", fixed = TRUE)[[1]][1] else ""
    calc_mh <- NA_real_
    btmp <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq1, mod1, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (!inherits(btmp, "try-error") && is.list(bip) &&
        !is.null(bip$pep_mass) && is.finite(bip$pep_mass))
      calc_mh <- bip$pep_mass
    base_pk <- max(ints0[is.finite(ints0) & ints0 > 0])
    bp_lab <- sub("E\\+0?", "E", sprintf("%.2E", base_pk))
    ann <- c(sprintf("Base peak = %s", bp_lab),
             if (is.finite(exp_mh) && is.finite(calc_mh) && calc_mh > 0)
               sprintf("Mass error = %+.4f Da / %+.2f ppm",
                       exp_mh - calc_mh, (exp_mh - calc_mh) / calc_mh * 1e6)
             else "Mass error = -")
    prot1 <- prots[k]
    acc1 <- if (is.na(prot1) || !nzchar(prot1) || prot1 %in% c("-", "NA")) "-"
            else sub("/.*$", "", prot1)
    top <- list(head = fmt_head(seq1, mod1, acc1, "-"),
                seq = seq1, mod = mod1,
                ions = all_ions(seq1, mod1, mz0, ints0))
    bot <- if (nzchar(seq2) && seq2 != "-") {
      list(head = fmt_head(seq2, "", "-", "-"),
           seq = seq2, mod = "",
           ions = all_ions(seq2, "", mz0, ints0))
    } else {
      list(head = "-", seq = "", mod = "", ions = NULL)
    }
    lines5 <- c(
      "Type: credible novel peptide",
      sprintf("Spectrum: %s", basename(fn)),
      sprintf("Novel peptide: %s | Modification: %s",
              seq1, if (nzchar(mod1) && mod1 != "NA") mod1 else "-"),
      sprintf("Closest annotated competitor: %s",
              if (nzchar(seq2) && seq2 != "-") seq2 else "-"),
      "Judgment: credible novel peptide (passed all six doubt rules)")
    # 自查整改（2026-09-26）：单谱渲染失败只跳过该谱，不中断批次
    ok <- tryCatch({
      draw_page(sp_obj, top, bot, lines5, chg_v, ann)
      TRUE
    }, error = function(e) {
      cat("[谱图 PDF] 跳过渲染失败谱：", fn, "（", conditionMessage(e), "）\n")
      FALSE
    })
    if (ok) drawn <- c(drawn, fn)
  }
  dev.off()
  if (length(drawn) == 0 && file.exists(pdf_f)) file.remove(pdf_f)
  drawn
}
# ---- 新肽段按 raw 分组 + 超阈值切片的驱动（2026-09-26 拆分口径）----
render_novel_grouped <- function(rows_idx) {
  if (length(rows_idx) == 0) return(invisible(character(0)))
  fns <- vapply(rows_idx, function(k) {
    v <- if ("File_Name" %in% colnames(nov)) as.character(nov$File_Name[k]) else ""
    if (is.na(v)) "" else v
  }, character(1))
  raws <- unname(vapply(fns, function(t) unname(title_raw_scan(t)["raw"]),
                        character(1)))
  tag <- "novel_peptide_spectra"
  used <- character(0)
  for (rw in unique(raws)) {
    sel <- rows_idx[raws == rw]
    n <- length(sel)
    nparts <- ceiling(n / ppf)
    tok <- raw_token(rw)
    while (tok %in% used) tok <- paste0(tok, "_2")
    used <- c(used, tok)
    for (prt in seq_len(nparts)) {
      chunk <- sel[seq.int(1 + (prt - 1) * ppf, min(n, prt * ppf))]
      pdf_name <- if (nparts == 1L) paste0(tag, "__", tok, ".pdf")
                  else sprintf("%s__%s__part%d.pdf", tag, tok, prt)
      drawn <- render_novel_batch(chunk, pdf_name)
      for (j in seq_along(drawn))
        IDX[[length(IDX) + 1L]] <<- c(drawn[j], tag, pdf_name, as.character(j))
    }
  }
  invisible(character(0))
}

# ---- 三类 PDF 输出（2026-09-26 拆分口径：按 raw 分组＋阈值切片）----
render_grouped(err_idx, "error_match_spectra",
               function(i) as.character(repA$File_Name[i]))
render_grouped(sus_idx, "ambiguous_matches_spectra",
               function(i) as.character(repA$File_Name[i]))
if (length(nov_idx)) render_novel_grouped(nov_idx)
if (length(IDX) > 0) {
  idx_df <- do.call(rbind, IDX)
  colnames(idx_df) <- c("File_Name", "category", "pdf_file", "page")
  idx_df <- idx_df[order(idx_df[, 1], idx_df[, 2], idx_df[, 4]), , drop = FALSE]
  write.table(idx_df, file.path(sta_dir, "spectrum_pdf_index.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("[谱图 PDF] 索引 spectrum_pdf_index.tsv（%d 行）\n", nrow(idx_df)))
  cat(sprintf("[spectrum PDF] index spectrum_pdf_index.tsv (%d rows)\n", nrow(idx_df)))
}
cat("[谱图 PDF] 完成\n")
cat("[spectrum PDF] done\n")
