# =============================================================================
# search_adapter_peaks.R - PEAKS DB search 肽段 CSV -> isoPSV 交换格式（m 槽候选）
# search_adapter_peaks.R - PEAKS DB-search peptide CSV -> isoPSV interchange (m-slot)
#
# PEAKS Studio 为 GUI 软件（Bioinformatics Solutions Inc.）。本适配器消费
# PEAKS DB search 的肽段导出 CSV（peptide.csv 或 protein-peptides.csv，列随
# 版本略有差异）。必需列：Peptide（修饰以 "+Δ" 括号内联在残基后，如
# M(+15.99)）、Scan（扫描号）、Source File（原始数据文件名）、z（电荷）。
# 可选列：-10lgP（置信分，每谱取最高一行）、Accession（蛋白号）、PTM（修饰
# 概要，无内联信息时仅作留存列）。
# PEAKS Studio is GUI software (Bioinformatics Solutions Inc.). This adapter
# consumes the PEAKS DB-search peptide export CSV (peptide.csv or
# protein-peptides.csv; column set varies slightly across versions). Required
# columns: Peptide (modifications inline as "(+delta)" right after the residue,
# e.g. M(+15.99)), Scan, Source File, z. Optional: -10lgP (best row kept per
# spectrum), Accession (protein IDs), PTM (kept as a trace column when the
# inline information is absent).
#
# 谱图标题按 raw.scan.scan.charge.0.dta 口径由 Source File + Scan + z 构建
# （与 MaxQuant 适配器一致；用户 MGF 标题须同口径，阶段 B 另有 raw+scan 兜底）。
# The spectrum title is built as raw.scan.scan.charge.0.dta from Source File +
# Scan + z (same convention as the MaxQuant adapter; the user's MGF titles must
# follow it, stage B also falls back to a raw+scan prefix).
#
# 用法：Rscript search_adapter_peaks.R <CSV 或文件夹> <输出.spectra>
# usage: Rscript search_adapter_peaks.R <CSV or folder> <out.spectra>
# =============================================================================
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- getwd()
if (length(script_arg) > 0 && !is.na(script_arg)) {
  cand <- tryCatch(dirname(normalizePath(sub("^--file=", "", script_arg),
                                         winslash = "/")), error = function(e) "")
  if (nzchar(cand) && dir.exists(cand)) script_dir <- cand
}
compat_src <- NULL
probe_dirs <- unique(c(script_dir, getwd(), dirname(script_dir),
                       file.path(dirname(script_dir), "code")))
for (d in probe_dirs) {
  p <- file.path(d, "compat.R")
  if (file.exists(p)) { compat_src <- p; break }
}
if (is.null(compat_src)) stop("compat.R not found next to the adapter")
source(compat_src)


# PEAKS 内联修饰位移 -> isoPSV 修饰名（与 Sage 适配器同一张表 + 标记修饰）
# PEAKS inline mod delta -> isoPSV mod name (same table as the Sage adapter
# plus label modifications)
PEAKS_MOD_NAME <- c(
  "15.9949" = "Oxidation", "15.99" = "Oxidation",
  "57.0215" = "Carbamidomethyl", "57.02" = "Carbamidomethyl",
  "0.984" = "Deamidated", "0.98" = "Deamidated",
  "42.0106" = "Acetyl", "42.01" = "Acetyl",
  "79.9663" = "Phospho", "79.97" = "Phospho",
  "-17.0265" = "Gln->pyro-Glu", "-17.03" = "Gln->pyro-Glu",
  "-18.0106" = "Glu->pyro-Glu", "-18.01" = "Glu->pyro-Glu",
  "14.0157" = "Methyl", "28.0313" = "Dimethyl",
  "10.0083" = "Label_13C(6)15N(4)", "6.0201" = "Label_13C(6)",
  "-1.0034" = "1_13C", "-2.0067" = "2_13C"
)

peaks_mod_string <- function(pep) {
  # "A(+42.01)DVLDLHE(+15.99)AGG..." -> "1,Acetyl[A];8,Oxidation[H];"
  out <- character(0)
  clean <- gsub("\\([^)]*\\)", "", pep)
  m <- gregexpr("\\([+-]?[0-9.]+?\\)", pep)[[1]]
  if (m[1] != -1) {
    for (i in seq_along(m)) {
      tok <- substr(pep, m[i], m[i] + attr(m, "match.length")[i] - 1L)
      dlt <- as.numeric(gsub("[()]", "", tok))
      hit <- which(abs(dlt - as.numeric(names(PEAKS_MOD_NAME))) < 0.05)
      if (length(hit) == 0) {
        cat("[PEAKS 适配] 未知修饰位移（忽略）：", dlt, "\n")
        next
      }
      nm <- PEAKS_MOD_NAME[[hit[1]]]
      # 修饰残基＝括号前一位残基（N 端修饰无前位，按位置 0/AnyN-term）
      before <- substr(pep, 1, m[i] - 1L)
      pos <- nchar(gsub("\\([^)]*\\)", "", before))   # 纯净序列中的位置
      if (pos == 0) {
        # 括号在序列首＝N 端修饰（无前位残基）
        out <- c(out, paste0("0,", nm, "[AnyN-term];"))
      } else {
        res <- if (pos >= 1L && pos <= nchar(clean)) substr(clean, pos, pos) else "X"
        out <- c(out, paste0(pos, ",", nm, "[", res, "];"))
      }
    }
  }
  paste(sort(out), collapse = "")
}

convert_peaks_csv <- function(csv_file, out_file) {
  csv_file <- resolve_adapter_input(csv_file, "\\.csv$")
  # P3-5编码容错读取（UTF-8 BOM/UTF-8/本地编码）
  d <- read_csv_tolerant(csv_file, stringsAsFactors = FALSE, check.names = FALSE,
                         quote = "\"", comment.char = "", fill = TRUE,
                         na.strings = character(0), colClasses = "character")
  need <- c("Peptide", "Scan", "Source File", "z")
  if (!all(need %in% names(d))) {
    stop("PEAKS CSV 缺列：", paste(setdiff(need, names(d)), collapse = ", "),
         "\n请从 PEAKS DB search 结果导出肽段 CSV（含 Peptide/Scan/Source File/z）\n",
         "missing columns; export the PEAKS DB-search peptide CSV including Peptide/Scan/Source File/z")
  }
  d <- d[nzchar(d[["Peptide"]]),, drop = FALSE]
  fn <- paste0(d[["Source File"]], ".", d[["Scan"]], ".", d[["Scan"]],
               ".", d[["z"]], ".0.dta")
  # 每谱取 -10lgP 最高一行（全候选模式 ISOPSV_CAND_ALL=1 时保留全部候选）
  sc <- if ("-10lgP" %in% names(d)) {
    s <- suppressWarnings(as.numeric(d[["-10lgP"]])); s[is.na(s)] <- -Inf; s
  } else rep(0, nrow(d))
  d <- d[order(sc, decreasing = TRUE),, drop = FALSE]
  fn <- fn[order(sc, decreasing = TRUE)]
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    keep <- !duplicated(fn)
    d <- d[keep,, drop = FALSE]; fn <- fn[keep]
  }
  seqs <- gsub("\\([^)]*\\)", "", d[["Peptide"]])
  mods <- vapply(d[["Peptide"]], peaks_mod_string, character(1))
  prot <- if ("Accession" %in% names(d)) {
    vapply(strsplit(d[["Accession"]], ";", fixed = TRUE), function(ps) {
      ps <- trimws(ps); ps <- ps[nzchar(ps)]
      if (length(ps) == 0) "" else paste0(paste(ps, collapse = "/"), "/")
    }, character(1))
  } else rep("", nrow(d))
  out <- data.frame(File_Name = fn, Sequence = seqs, Modification = mods,
                    Proteins = prot, Charge = d[["z"]],
                    stringsAsFactors = FALSE, check.names = FALSE)
  # PTM 概要留存列（内联修饰缺失时便于溯源；合并步骤只读固定列，多余列无害）
  if ("PTM" %in% names(d)) out$PEAKS_PTM <- d[["PTM"]]
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("PEAKS CSV -> %s：%d 行（-10lgP 最佳每谱去重）\n", out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_peaks_csv(a[1], a[2])
}
