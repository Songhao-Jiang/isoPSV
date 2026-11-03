# =============================================================================
# search_adapter_denovo.R —— Casanovo / PrimeNovo de novo 结果转换件
#
# 两个外部 de novo 引擎（用户自添加槽，默认无行）：
#   casanovo：casanovo.exe sequence <mgf> -m <ckpt> -d <dir> -o <stem> -f
#               → mzTab（PSM 段）；spectra_ref=ms_run[1]:index=N 需按 MGF
#               谱序回映射标题；modifications 列为 mzTab 风格
#               "14-Deamidated (N):UNIMOD:7"（null＝无修饰）
#   primenovo：python -m PrimeNovo.PrimeNovo --mode=denovo --peak_path=<mgf>
#               --model=<ckpt>（repo 根目录运行）→ denovo.tsv
#               列＝label（＝MGF TITLE，直接可用）/prediction（ProForma）/
#               charge/score
#
# 本件把两类输出统一转成与 denovo_csv_to_spectra 相同的三列 .spectra
# （File_Name / Sequence / Modification），供合并层按 punifind_denovo 同
# 口径消费（denovo 槽、数据库回帖）。未知名/质量位移修饰剥除并计数（与
#  AlphaPeptDeep 映射层同纪律）。
#
# 用法：
#   Rscript search_adapter_denovo.R --convert casanovo  --in <mztab> --mgf <mgf> --out <spectra>
#   Rscript search_adapter_denovo.R --convert primenovo --in <tsv>             --out <spectra>
#   Rscript search_adapter_denovo.R --self-test
# 退出码：0 成功；2 参数/输入错误；3 输出为空（仍写表头）。
# 依赖：仅 base R。
# =============================================================================

parse_cli <- function(args = commandArgs(TRUE)) {
  parsed <- list(); i <- 1L
  while (i <= length(args)) {
    tok <- args[i]
    if (startsWith(tok, "--")) {
      k <- sub("^--", "", tok)
      if (grepl("=", k, fixed = TRUE)) {
        kv <- strsplit(k, "=", fixed = TRUE)[[1]]
        parsed[[kv[1]]] <- paste(kv[-1], collapse = "=")
      } else {
        parsed[[k]] <- if (i < length(args) && !startsWith(args[i + 1], "--")) {
          v <- args[i + 1]; i <- i + 1L; v
        } else TRUE
      }
    }
    i <- i + 1L
  }
  parsed
}

# ---- ProForma / mzTab 修饰 → pFind 族 "pos,Name[res]" ------------------------
# 未知名/未映射质量位移：剥除并计数（返回属性 unknown 供日志）
DN_NAME_ALIAS <- c(
  "Deamidated" = "Deamidated", "Oxidation" = "Oxidation",
  "Carbamidomethyl" = "Carbamidomethyl", "Acetyl" = "Acetyl",
  "Carbamyl" = "Carbamyl", "Phospho" = "Phospho",
  "Gln->pyro-Glu" = "Gln->pyro-Glu", "Glu->pyro-Glu" = "Glu->pyro-Glu"
)
# 常见质量位移 → pFind 名（残基在转换时按所在残基写入）
DN_MASS_ALIAS <- c(
  "57.021" = "Carbamidomethyl", "57.022" = "Carbamidomethyl",
  "15.995" = "Oxidation", "15.994" = "Oxidation",
  "0.984" = "Deamidated", "0.981" = "Deamidated",
  "42.011" = "Acetyl", "42.012" = "Acetyl",
  "43.006" = "Carbamyl", "43.007" = "Carbamyl",
  "79.966" = "Phospho", "79.963" = "Phospho",
  "-17.027" = "Gln->pyro-Glu", "-18.011" = "Glu->pyro-Glu",
  # 质量位移查表按绝对值（N 端前缀 [-17.027]PEP… 的 pyro-Glu 类）
  "17.027" = "Gln->pyro-Glu", "18.011" = "Glu->pyro-Glu"
)

# ProForma 肽段（"EDL…EN[Deamidated]GR" / "PEP[+57.021]TIDE" / N 端
# "Acetyl@…"-前缀暂按无修饰处理并计数）→ list(seq, mod, unknown)
proforma_to_pf <- function(p) {
  unknown <- character(0)
  if (is.na(p) || !nzchar(p)) return(list(seq = "", mod = "", unknown = unknown))
  p <- trimws(as.character(p))
  mods <- character(0); out <- character(0)
  # 逐字符扫描，命中 [..] 即解析一个修饰
  i <- 1L; pos0 <- 0L
  while (i <= nchar(p)) {
    ch <- substr(p, i, i)
    if (identical(ch, "[")) {
      j <- regexpr("]", substring(p, i), fixed = TRUE)
      if (j < 0) break
      inner <- substring(p, i + 1L, i + j - 2L)
      aa <- if (pos0 >= 1L && pos0 <= length(out)) out[pos0] else ""
      nm <- ""
      if (grepl("^[+-]?[0-9.]+$", inner)) {  # 质量位移风格
        key <- substring(sprintf("%.3f", abs(as.numeric(inner))), 1, 6)
        nm <- unname(DN_MASS_ALIAS[key] %||% "")
      } else {
        nm0 <- sub("\\(.*\\)$", "", inner)
        nm <- unname(DN_NAME_ALIAS[nm0] %||% "")
      }
      if (!nzchar(nm)) {
        unknown <- c(unknown, inner)
      } else {
        site <- if (nzchar(aa)) aa else "AnyN-term"
        mods <- c(mods, sprintf("%d,%s[%s]", pos0, nm, site))
      }
      i <- i + j                            # 跳过 "]"
    } else if (ch %in% LETTERS) {
      out <- c(out, ch); pos0 <- pos0 + 0L; i <- i + 1L
      # 注：残基位置 pos0 在遇到 [..] 时先加一，这里不加——统一在修饰分支推进
      pos0 <- length(out)
    } else {
      i <- i + 1L                           # 其他字符（如 N 端前缀分隔）跳过
    }
  }
  list(seq = paste(out, collapse = ""), mod = paste(mods, collapse = ";"),
       unknown = unknown)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a

# mzTab modifications 列（"14-Deamidated (N):UNIMOD:7;0-Acetyl:UNIMOD:1"）
# → list(mod="14,Deamidated[N];0,Acetyl[AnyN-term];", unknown=…)
# 手工解析（pos-名(残基):UNIMOD:id 四段式），不用正则——shell/heredoc 的
# 反斜杠转义坑与 TRE/PCRE 差异一并绕开。
mztab_mods_to_pf <- function(s) {
  if (is.na(s) || !nzchar(s) || identical(tolower(s), "null"))
    return(list(mod = "", unknown = character(0)))
  out <- character(0); unknown <- character(0)
  for (tok in strsplit(s, ";", fixed = TRUE)[[1]]) {
    tok <- trimws(tok)
    if (!nzchar(tok)) next
    parts <- strsplit(tok, "-", fixed = TRUE)[[1]]
    pos <- parts[1]; rest <- paste(parts[-1], collapse = "-")
    ok <- grepl("^[0-9]+$", pos)
    # 剥 ":UNIMOD:…" 尾巴
    if (grepl(":", rest, fixed = TRUE)) rest <- strsplit(rest, ":", fixed = TRUE)[[1]][1]
    # 残基段 " (N)" → 名 + 残基
    res <- "AnyN-term"
    op <- regexpr(" (", rest, fixed = TRUE)
    if (op > 0L && endsWith(rest, ")")) {
      res <- substring(rest, op + 2L, nchar(rest) - 1L)
      rest <- trimws(substring(rest, 1L, op - 1L))
    }
    if (!ok || !nzchar(rest)) { unknown <- c(unknown, tok); next }
    nm <- unname(DN_NAME_ALIAS[trimws(rest)] %||% "")
    if (!nzchar(nm)) { unknown <- c(unknown, tok); next }
    out <- c(out, sprintf("%s,%s[%s]", pos, nm, res))
  }
  list(mod = paste(out, collapse = ";"), unknown = unknown)
}

# ---- MGF 谱序 → 标题（Casanovo spectra_ref 索引回映射） ----------------------
mgf_titles <- function(mgf) {
  t <- character(0)
  for (ln in readLines(mgf, warn = FALSE)) {
    if (startsWith(ln, "TITLE=")) t <- c(t, substring(ln, 7))
  }
  t
}

# ---- Casanovo mzTab → 三列 .spectra ----------------------------------------
casanovo_mztab_to_spectra <- function(mztab, mgf, out_path) {
  lines <- readLines(mztab, warn = FALSE, encoding = "UTF-8")
  psm <- lines[startsWith(lines, "PSM\t")]
  hdr <- strsplit(lines[startsWith(lines, "PSH\t")][1], "\t", fixed = TRUE)[[1]]
  if (length(hdr) == 0L || length(psm) == 0L)
    stop("mzTab 无 PSM 行 / no PSM rows: ", mztab)
  titles <- mgf_titles(mgf)
  # 注：PSH/PSM 前缀列保留在位（表头与数据行含前缀逐列对齐，mzTab 数据行
  # 可能省略空尾列，剥离前缀反而引入错位）；按列名 match 取列。
  rows <- strsplit(psm, "\t", fixed = TRUE)
  ci <- function(nm) match(nm, hdr)
  i_sq <- ci("sequence"); i_md <- ci("modifications"); i_ref <- ci("spectra_ref")
  get_idx <- function(ref) {
    m <- regmatches(ref, regexec("index=([0-9]+)", ref))[[1]]
    if (length(m) == 0L) NA_integer_ else as.integer(m[2])
  }
  idx0 <- vapply(rows, function(r) get_idx(r[i_ref]), integer(1))
  hit0 <- sum(!is.na(idx0) & idx0 + 1L <= length(titles), na.rm = TRUE)
  hit1 <- sum(!is.na(idx0) & idx0 <= length(titles), na.rm = TRUE)
  base_off <- if (hit0 >= hit1) 1L else 0L
  out <- data.frame(File_Name = character(0), Sequence = character(0),
                    Modification = character(0), stringsAsFactors = FALSE)
  n_unk <- 0L; unk_names <- character(0)
  for (r in rows) {
    k <- idx0 <- get_idx(r[i_ref])
    if (is.na(k)) next
    ti <- k + base_off
    if (ti < 1L || ti > length(titles)) next
    sq <- r[i_sq]
    if (!grepl("^[A-Z]+$", sq)) next          # 含非常规残基的预测弃用
    md_raw <- r[i_md]
    mdl <- mztab_mods_to_pf(md_raw)
    md <- mdl$mod
    u <- mdl$unknown
    n_unk <- n_unk + length(u); unk_names <- c(unk_names, u)
    out <- rbind(out, data.frame(File_Name = titles[ti], Sequence = sq,
                                 Modification = md, stringsAsFactors = FALSE))
  }
  out <- out[!duplicated(out$File_Name),, drop = FALSE]
  write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[casanovo] 转换 ", nrow(out), " 谱（mzTab ", length(psm), " 行；索引基 +",
      base_off - 1L, "；修饰未知剥除 ", n_unk, "）-> ", out_path, "\n", sep = "")
  if (length(unk_names))
    cat("[casanovo] 未知修饰（前 10）：", paste(head(unique(unk_names), 10), collapse = " | "), "\n", sep = "")
  nrow(out)
}

# ---- PrimeNovo denovo.tsv → 三列 .spectra ----------------------------------
primenovo_tsv_to_spectra <- function(tsv, out_path) {
  raw <- read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE,
                    na.strings = character(0))
  cn <- tolower(colnames(raw))
  i_lb <- grep("^label$|^title", cn)[1]
  i_pr <- grep("^prediction$|^peptide", cn)[1]
  if (is.na(i_lb) || is.na(i_pr))
    stop("primenovo tsv columns not recognized: ", paste(colnames(raw), collapse = ","))
  out <- data.frame(File_Name = character(0), Sequence = character(0),
                    Modification = character(0), stringsAsFactors = FALSE)
  n_unk <- 0L; unk_names <- character(0)
  for (i in seq_len(nrow(raw))) {
    pf <- proforma_to_pf(raw[[i_pr]][i])
    n_unk <- n_unk + length(pf$unknown)
    unk_names <- c(unk_names, pf$unknown)
    if (!nzchar(pf$seq) || !grepl("^[A-Z]+$", pf$seq)) next
    out <- rbind(out, data.frame(File_Name = as.character(raw[[i_lb]][i]),
                                 Sequence = pf$seq, Modification = pf$mod,
                                 stringsAsFactors = FALSE))
  }
  out <- out[!duplicated(out$File_Name),, drop = FALSE]
  write.table(out, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[primenovo] 转换 ", nrow(out), " 谱（修饰未知剥除 ", n_unk, "）-> ", out_path, "\n", sep = "")
  if (length(unk_names))
    cat("[primenovo] 未知修饰（前 10）：", paste(head(unique(unk_names), 10), collapse = " | "), "\n", sep = "")
  nrow(out)
}

self_test <- function() {
  cat("[self-test] search_adapter_denovo\n")
  pf <- proforma_to_pf("EDLHHEGFGYVEHN[Deamidated]GR")
  stopifnot(pf$seq == "EDLHHEGFGYVEHNGR", pf$mod == "14,Deamidated[N]")
  cat("  [ok] ProForma 具名修饰\n")
  pf2 <- proforma_to_pf("PEPT[+57.021]IDEK")
  stopifnot(pf2$seq == "PEPTIDEK", pf2$mod == "4,Carbamidomethyl[T]")
  cat("  [ok] ProForma 质量位移\n")
  pf3 <- proforma_to_pf("PEP[ZzzMod]TIDE")
  stopifnot(pf3$seq == "PEPTIDE", pf3$mod == "", length(pf3$unknown) == 1L)
  cat("  [ok] 未知名剥除计数\n")
  mdl <- mztab_mods_to_pf("14-Deamidated (N):UNIMOD:7;0-Acetyl:UNIMOD:1")
  stopifnot(mdl$mod == "14,Deamidated[N];0,Acetyl[AnyN-term]", length(mdl$unknown) == 0L)
  cat("  [ok] mzTab 修饰列\n")
  tf <- tempfile(fileext = ".mgf"); tf2 <- tempfile(fileext = ".mztab")
  writeLines(c("BEGIN IONS", "TITLE=specA.2.0.dta", "END IONS",
               "BEGIN IONS", "TITLE=specB.2.0.dta", "END IONS"), tf)
  writeLines(c("PSH\tsequence\tPSM_ID\tmodifications\tspectra_ref",
               "PSM\tPEPTIDEK\t1\tnull\tms_run[1]:index=0",
               "PSM\tEDITIDEK\t2\t3-Oxidation (M):UNIMOD:35\tms_run[1]:index=1"), tf2)
  o <- tempfile(fileext = ".spectra")
  n <- casanovo_mztab_to_spectra(tf2, tf, o)
  d <- read.delim(o, stringsAsFactors = FALSE)
  stopifnot(n == 2L, d$File_Name[1] == "specA.2.0.dta", d$Sequence[2] == "EDITIDEK")
  cat("  [ok] mzTab index 回映射（0 基）\n")
  cat("[self-test] 全部通过\n")
  TRUE
}

main <- function() {
  a <- parse_cli()
  if (isTRUE(a$`self-test`)) return(if (self_test()) 0L else 1L)
  a$infile <- if (is.null(a$infile)) a[["in"]] else a$infile
  if (is.null(a$convert) || is.null(a$infile) || is.null(a$out)) {
    cat("用法/usage: search_adapter_denovo.R --convert casanovo|primenovo --in <out> [--mgf <mgf>] --out <spectra>\n")
    return(2L)
  }
  if (identical(a$convert, "casanovo")) {
    if (is.null(a$mgf)) { cat("--convert casanovo 需要 --mgf\n"); return(2L) }
    n <- tryCatch(casanovo_mztab_to_spectra(a$infile, a$mgf, a$out),
                  error = function(e) { cat("[casanovo] 转换失败：", conditionMessage(e), "\n"); NULL })
  } else if (identical(a$convert, "primenovo")) {
    n <- tryCatch(primenovo_tsv_to_spectra(a$infile, a$out),
                  error = function(e) { cat("[primenovo] 转换失败：", conditionMessage(e), "\n"); NULL })
  } else { cat("未知 --convert：", a$convert, "\n"); return(2L) }
  if (is.null(n)) return(2L)
  if (n == 0L) return(3L)
  0L
}

if (sys.nframe() == 0L) quit(status = main())
