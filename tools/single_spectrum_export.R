# =============================================================================
# 单谱图输出.R - pFind-Filtered.spectra 单谱图输出软件（与谱图可视化同格式）
# 单谱图输出.R - standalone single-spectrum figure exporter (same figure format)
#
# 输入：① pFind-Filtered.spectra（列：File_Name、Scan_No、Exp.MH+、Charge、
#   Sequence、Modification、Raw_Score、Final_Score、Proteins、Positions）；
#   ② MGF 文件夹（按 File_Name 标题定位谱图，文件夹内全部 .mgf）。
# 输出：每谱一页 PDF＋300dpi PNG，文件名＝scan<Scan_No>（谱图号）；
#   页面内容＝谱图号/参数/价态等信息行＋离子梯＋谱图（Base peak/Mass error 注释、
#   0–100 相对强度、同峰合并与窜开标注）＋子离子 ppm 图。
# 用法（环境变量）：
#   PFFILE=<pFind-Filtered.spectra 路径>  MGFDIR=<MGF 文件夹>
#   SELECT=<行号或 Scan_No，逗号分隔；缺省全部>  LIMIT=<最多张数>
# 例：SELECT=33516,8094 只画这两个 scan；不设 LIMIT 则全部输出。
# 依赖：isoPSV 客户端库（compat.R/ReadMGF_script.R/chimeric_check.R，仅 base R）。
# =============================================================================
args0 <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "",
  args0[grepl("--file=", args0)][1]), winslash = "/"))
OUT <- script_dir
CODE <- script_dir                                            # 库版与全流程同目录，直接使用自身所在 code 目录
PF_FILE <- Sys.getenv("PFFILE", "")
MGF_DIR <- Sys.getenv("MGFDIR", "")
if (!nzchar(PF_FILE) || !nzchar(MGF_DIR))
  stop("请设置环境变量 PFFILE 与 MGFDIR（见脚本头部用法）")
SEL_ENV <- Sys.getenv("SELECT", "")
LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))

old <- setwd(CODE)
source("compat.R"); source("ReadMGF_script.R")
source("chimeric_check.R")
setwd(old)
aa <- read.delim(file.path(CODE, "AA_mass_dataframe.txt"), sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# ---- 输入解析 -------------------------------------------------------------------
pf <- read.delim(PF_FILE, sep = "\t", check.names = FALSE,
                 na.strings = character(0), quote = "", fill = TRUE)
need_cols <- c("File_Name", "Scan_No", "Exp.MH+", "Charge", "Sequence", "Modification")
miss_c <- setdiff(need_cols, names(pf))
if (length(miss_c)) stop("pFind-Filtered.spectra 缺列：", paste(miss_c, collapse = ","))
for (cn in need_cols) pf[[cn]] <- as.character(pf[[cn]])
pf$Scan_No <- suppressWarnings(as.numeric(pf$Scan_No))
pf$Charge <- suppressWarnings(as.integer(pf$Charge))
pf$ExpMH <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pf[["Exp.MH+"]])))
if (!"Raw_Score" %in% names(pf)) pf$Raw_Score <- "-"
if (!"Final_Score" %in% names(pf)) pf$Final_Score <- "-"
if (!"Proteins" %in% names(pf)) pf$Proteins <- ""
if (!"Positions" %in% names(pf)) pf$Positions <- ""

if (nzchar(SEL_ENV)) {                                       # 选择：行号或 Scan_No
  want <- suppressWarnings(as.numeric(trimws(strsplit(SEL_ENV, "[,;]")[[1]])))
  sel <- which(pf$Scan_No %in% want | seq_len(nrow(pf)) %in% want)
} else sel <- seq_len(nrow(pf))
if (is.finite(LIMIT)) sel <- head(sel, LIMIT)
cat(sprintf("pFind-Filtered %d 行；选择 %d 行\n", nrow(pf), length(sel)))

old2 <- setwd(CODE)
lu <- load_unlabeled_spectra(MGF_DIR)
setwd(old2)
if (is.null(lu)) stop("MGF 文件夹无 .mgf 文件")
find_spec <- function(fn) {
  raw <- unname(title_raw_scan(fn)["raw"])
  objs <- lu$objs[[raw]]
  if (is.null(objs)) return(NULL)
  li <- match(fn, lu$titles[[raw]])
  if (is.na(li)) {
    key <- title_scan_key(fn)
    ci <- which(startsWith(lu$titles[[raw]], paste0(key, ".")))
    if (length(ci) == 0) return(NULL)
    li <- ci[1]
  }
  objs[[li]]
}

# ---- 与谱图可视化一致的绘制引擎（单候选版）---------------------------------------
ION_COL <- c(a = "#004000", b = "#008000", c = "#CA1D52", x = "#650E29",
             y = "#CA1D52", z = "#943AA4", M = "#808080", I = "#FF8C00", O = "#9ACD32")
ION_PRI <- c(b = 1, y = 2, a = 3, c = 4, x = 5, z = 6, M = 7, I = 8, O = 9)
CO_M <- 27.994915; NH3_M <- 17.026549; H2O_M <- 18.010565; H_M <- 1.007276

all_ions <- function(seq, mod, mz, ints, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) ||
        length(bip) == 1 || is.null(bip$ions_table)) next
    t <- bip$ions_table; n <- t$Number
    theo <- data.frame(
      mz  = c(t$b_1_charge, t$b_1_charge - CO_M, t$b_1_charge + NH3_M,
              t$b_2_charge, t$b_2_charge - CO_M/2, t$b_2_charge + NH3_M/2,
              t$y_1_charge, t$y_1_charge + CO_M, t$y_1_charge - NH3_M,
              t$y_2_charge, t$y_2_charge + CO_M/2, t$y_2_charge - NH3_M/2),
      lab = c(paste0("b", n, "+"), paste0("a", n, "+"), paste0("c", n, "+"),
              paste0("b", n, "++"), paste0("a", n, "++"), paste0("c", n, "++"),
              paste0("y", n, "+"), paste0("x", n, "+"), paste0("z", n, "+"),
              paste0("y", n, "++"), paste0("x", n, "++"), paste0("z", n, "++")),
      ion = c(rep(c("b", "a", "c"), each = length(n)),
              rep(c("b", "a", "c"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n))),
      num = c(rep(n, 3), rep(n, 3), rep(n, 3), rep(n, 3)),
      chg = c(rep(1L, 3 * length(n)), rep(2L, 3 * length(n)),
              rep(1L, 3 * length(n)), rep(2L, 3 * length(n))),
      stringsAsFactors = FALSE)
    if (is.finite(bip$pep_mass)) {
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

draw_ladder <- function(sq, ions, mod = "", chg_v = 0L) {
  n <- nchar(sq)
  if (n < 2L) return(invisible())
  u <- par("usr")
  x0 <- u[1] + 0.015 * (u[2] - u[1]); x1 <- u[2] - 0.015 * (u[2] - u[1])
  u_per_in <- abs(diff(grconvertX(c(0, 1), from = "inches", to = "user")))
  slot_u <- min(0.16 * u_per_in, (x1 - x0) / n)
  xc0 <- (x0 + x1) / 2 - n * slot_u / 2
  xs <- xc0 + (seq_len(n) - 0.5) * slot_u
  gap <- slot_u; ymid <- 1.30; v_h <- 0.045
  site_by <- b_sites <- y_sites <- integer(0)
  if (!is.null(ions)) {
    hit <- ions[!is.na(ions$peak), ]
    if (nrow(hit)) {
      hb <- hit[hit$ion == "b", ]
      if (nrow(hb)) b_sites <- sort(unique(as.integer(hb$num)))
      hy <- hit[hit$ion == "y", ]
      if (nrow(hy)) y_sites <- sort(unique(n - as.integer(hy$num)))
      site_by <- sort(unique(c(b_sites, y_sites)))
    }
  }
  op <- par(xpd = TRUE); on.exit(par(op), add = TRUE)
  slot_in <- slot_u / u_per_in
  cexs <- min(8 / 12, max(6 / 12,
           0.9 * slot_in / strwidth("M", cex = 1, units = "inches")))
  cexn <- min(6 / 12, max(5 / 12,
           0.9 * slot_in / strwidth("b88", cex = 1, units = "inches")))
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
         as.expression(lapply(y_sites, function(k) bquote(y[.(n - k)]))),
         cex = cexn, col = ION_COL[["y"]], font = 2)
  }
  if (is.finite(chg_v) && chg_v > 0L)
    text(xc0 - gap * 1.5, ymid, sprintf("%d+ ", chg_v),
         cex = cexs, font = 2, adj = 1, col = "#CA1D52")
  lab_col <- rep("black", n); lab_font <- rep(1, n)
  if (!is.null(mod) && nzchar(mod)) {
    pm <- suppressWarnings(as.integer(sub(",.*$", "",
      strsplit(mod, ";", fixed = TRUE)[[1]])))
    pm <- pm[!is.na(pm) & pm >= 1L & pm <= n]
    lab_col[pm] <- ION_COL[["I"]]; lab_font[pm] <- 2
  }
  text(xs, ymid, strsplit(sq, "")[[1]], cex = cexs, font = lab_font, col = lab_col)
}

# 单谱一页：信息行＋离子梯＋谱图＋ppm
draw_single <- function(spec_obj, row, lines3, chg_v, ann, head_txt) {
  layout(matrix(c(1, 2), 2, 1), heights = c(1.0, 0.35))
  par(oma = c(3, 3, 5.2, 3))
  par(cex = 1)                                    # layout 自动 cex 缩放复位
  mz <- as.numeric(spec_obj@spectrum[[1]]); ints <- as.numeric(spec_obj@spectrum[[2]])
  o <- order(mz, na.last = NA); mz <- mz[o]; ints <- ints[o]
  ints_n <- ints / max(ints[is.finite(ints) & ints > 0])
  xmax <- max(mz, na.rm = TRUE) * 1.03
  # —— 谱图面板 ——
  par(mar = c(0.5, 4, 3.9, 1))
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
  th <- all_ions(as.character(row$Sequence), as.character(row$Modification), mz, ints)
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
      y0 <- pin + 0.03
      y_per_in <- abs(diff(grconvertY(c(0, 1), from = "inches", to = "user")))
      len <- strwidth(labs, units = "inches", cex = 8/12) * y_per_in   # 8 号
      cl <- cumsum(c(0, diff(pmz) > 50))
      for (g in unique(cl)) {
        k <- which(cl == g)
        if (length(k) < 2) next
        k <- k[order(pin[k])]
        for (j in seq_along(k)[-1L]) {
          i0 <- k[j - 1L]; i1 <- k[j]
          floor_y <- y0[i0] + len[i0] + 0.008
          if (y0[i1] < floor_y) y0[i1] <- min(floor_y, 1.10)
        }
      }
      text(pmz, y0, labs, cex = 8/12, srt = 90, col = cols, adj = 0)  # 8 号
    }
  }
  draw_ladder(as.character(row$Sequence), th, as.character(row$Modification), chg_v)
  mtext(head_txt, side = 3, line = 4.8, cex = 0.83, adj = 0)
  # —— ppm 面板 ——
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
  # —— 信息行三行全 10 号（信息↔序列行零间距；序列↔梯、梯↔谱正常间距）——
  title(lines3[1], outer = TRUE, line = 4.2, cex.main = 0.83)
  title(lines3[2], outer = TRUE, line = 3.1, cex.main = 0.83)
  title(lines3[3], outer = TRUE, line = 2.0, cex.main = 0.83)
}

# ---- 主循环：逐行输出 PDF＋300dpi PNG --------------------------------------------
fmt_num <- function(v, dig) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)
}
n_ok <- 0L
for (i in sel) {
  row <- pf[i, ]
  fn <- as.character(row$File_Name)
  spec_obj <- find_spec(fn)
  if (is.null(spec_obj)) { cat("跳过（MGF 无该谱）:", fn, "\n"); next }
  mz0 <- as.numeric(spec_obj@spectrum[[1]]); ints0 <- as.numeric(spec_obj@spectrum[[2]])
  chg_v <- tryCatch(as.integer(spec_obj@charge), error = function(e) NA_integer_)
  if (is.na(chg_v)) chg_v <- row$Charge
  if (is.na(chg_v) || chg_v < 1L) chg_v <- 0L
  exp_mh <- tryCatch(as.numeric(spec_obj@pepMass), error = function(e) NA_real_)
  calc_mh <- NA_real_
  btmp <- try(suppressWarnings(invisible(capture.output(
    bip <- by_ions_position(as.character(row$Sequence),
                            as.character(row$Modification), "Unlabel", aa, aa_tab)))),
    silent = TRUE)
  if (!inherits(btmp, "try-error") && !is.null(bip) && is.finite(bip$pep_mass))
    calc_mh <- bip$pep_mass
  if (is.finite(exp_mh)) exp_mh <- exp_mh * max(1L, chg_v) - (max(1L, chg_v) - 1L) * 1.007276
  base_pk <- max(ints0[is.finite(ints0) & ints0 > 0])
  bp_lab <- sub("E\\+0?", "E", sprintf("%.2E", base_pk))
  ann <- c(sprintf("Base peak = %s", bp_lab),
           if (is.finite(exp_mh) && is.finite(calc_mh) && calc_mh > 0)
             sprintf("Mass error = %+.4f Da / %+.2f ppm",
                     exp_mh - calc_mh, (exp_mh - calc_mh) / calc_mh * 1e6)
           else "Mass error = -")
  # 蛋白与侧翼（Positions＝pos,前,后/ 与 Proteins＝ACC/）
  acc <- "-"; pos_info <- "-"
  if (!is.na(row$Proteins) && nzchar(as.character(row$Proteins)))
    acc <- strsplit(as.character(row$Proteins), "/", fixed = TRUE)[[1]][1]
  if (!is.na(row$Positions) && nzchar(as.character(row$Positions))) {
    ft <- strsplit(strsplit(as.character(row$Positions), "/", fixed = TRUE)[[1]][1],
                   ",", fixed = TRUE)[[1]]
    if (length(ft) >= 3) pos_info <- sprintf("%s[%s]%s", ft[2], ft[1], ft[3])
  }
  pre <- sub("\\[.*$", "", pos_info); post <- sub("^.*\\]", "", pos_info)
  head_txt <- sprintf("| %s | (%s)%s(%s) | %s |", acc, pre,
                      as.character(row$Sequence), post,
                      if (!is.na(row$Modification) && nzchar(as.character(row$Modification)))
                        as.character(row$Modification) else "-")
  psm_sc <- "-"                                           # PSM score（b/y 离子强度占比）
  btmp2 <- try(suppressWarnings(invisible(capture.output(
    b2 <- ions_match(bip, spec_obj@spectrum, 20, "ppm")))), silent = TRUE)
  if (!inherits(btmp2, "try-error") && !is.null(b2) && is.finite(b2$PSM_score_by))
    psm_sc <- sprintf("%.2f", as.numeric(b2$PSM_score_by))
  lines3 <- c(
    sprintf("Spectrum: %s", fn),                           # 谱图号（仅 Spectrum）
    sprintf("Scores: PSM score = %s | Raw score = %s | Final score = %s",
            psm_sc, fmt_num(row$Raw_Score, 2L), fmt_num(row$Final_Score, 2L)),
    sprintf("Protein: %s | %s | %s",                       # 蛋白质|序列|修饰
            acc, as.character(row$Sequence),
            if (!is.na(row$Modification) && nzchar(as.character(row$Modification)))
              as.character(row$Modification) else "-"))
  tag <- if (is.na(row$Scan_No)) i else row$Scan_No
  pdf_f <- file.path(OUT, sprintf("scan%s.pdf", tag))
  png_f <- file.path(OUT, sprintf("scan%s.png", tag))
  pdf(pdf_f, width = 12.8, height = 7.2)        # 16:9 横向画布
  draw_single(spec_obj, row, lines3, chg_v, ann, head_txt)
  dev.off()
  png(png_f, width = 12.8, height = 7.2, units = "in", res = 300)  # 300dpi
  draw_single(spec_obj, row, lines3, chg_v, ann, head_txt)
  dev.off()
  n_ok <- n_ok + 1L
  if (n_ok %% 20 == 0L) cat(sprintf("进度 %d/%d\n", n_ok, length(sel)))
}
cat(sprintf("完成：输出 %d 张（PDF＋300dpi PNG，文件名 scan<谱图号>）\n", n_ok))
