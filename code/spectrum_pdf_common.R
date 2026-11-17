# spectrum_pdf_common.R - shared rendering machinery for the two
# spectrum-PDF exporters (v2.0.181). The nine functions below were
# byte-identical copies in spectrum_pdf_export.R (labeled side) and
# spectrum_pdf_export_A.R (unlabeled side); they now live here once.
# The judgment/driver logic and the five functions that genuinely
# differ (parse_col separators, fmt_head signature, mcol, render_batch,
# render_novel_batch) stay in each script.
# 本文件＝两个谱图 PDF 导出器的共享渲染机械（v2.0.181）：以下九个函数
# 原为两脚本中的逐字拷贝，现收敛于此；判定/驱动逻辑与五个真正有语义
# 差异的函数（parse_col 分隔符、fmt_head 签名、mcol、render_batch、
# render_novel_batch）仍留在各自脚本。

cval <- function(v, dflt = "-") {                              # 缺失行安全取值
  x <- as.character(v)
  x[is.na(x)] <- dflt
  x
}

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

raw_sm <- function(i, k) {
  cc <- cand_raw[[i]]
  if (length(cc) < k) return(c("", ""))
  s <- trimws(cc[k])
  m <- regmatches(s, regexpr("[(].*[)]", s))
  c(sub("[(].*$", "", s), if (length(m) == 1L) substr(m, 2, nchar(m) - 1L) else "")
}

fmt_num <- function(v, dig) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)
}

raw_token <- function(rw) {
  s <- substr(rw, 1, 80)
  s <- gsub("[^A-Za-z0-9_.-]", "_", s)
  if (!nzchar(s)) s <- "_"
  s
}

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
