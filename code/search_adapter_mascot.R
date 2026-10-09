# =============================================================================
# search_adapter_mascot.R - Mascot CSV 导出 -> isoPSV 交换格式（n 槽候选）
# search_adapter_mascot.R - Mascot CSV export -> isoPSV interchange (n-slot)
#
# Mascot 为 GUI/服务器软件（Matrix Science），无命令行自动搜库。本适配器消费
# Mascot 的 "Export Search Results" CSV 导出件（列规范见 Matrix Science
# csv_headers.html）。必需列：pep_scan_title（谱图标题，取自峰列表，直通
# File_Name）、pep_seq、pep_rank（每谱保留 top-1）、pep_var_mod（修饰名列表）
# 与 pep_var_mod_pos（等长数字串，非零数字＝1 基索引指向修饰名列表，多位点
# 同修饰时 pep_summed_mod_pos 给第二修饰）。prot_acc/pep_exp_z 可选补充。
# Mascot is GUI/server software (Matrix Science) with no CLI auto-search. This
# adapter consumes the "Export Search Results" CSV (column spec at Matrix
# Science csv_headers.html). Required columns: pep_scan_title (the peak-list
# scan title, passed through to File_Name), pep_seq, pep_rank (top-1 kept),
# pep_var_mod (modification name list) and pep_var_mod_pos (equal-length digit
# string, non-zero digit = 1-based index into the name list; pep_summed_mod_pos
# holds the second modification at multi-mod sites). prot_acc / pep_exp_z are
# optional extras.
#
# 固定修饰注入：Mascot 导出只含可变修饰；按 isoPSV 常规搜库口径为序列中每个
# C 注入 Carbamidomethyl[C]（第三个参数可改，如 "none" 关闭）。
# Fixed-mod injection: the export lists variable mods only; per the isoPSV
# standard search config Carbamidomethyl[C] is injected for every C (3rd
# argument overrides, e.g. "none" to disable).
#
# 用法：Rscript search_adapter_mascot.R <CSV 或文件夹> <输出.spectra> [固定修饰名，默认 Carbamidomethyl (C)]
# usage: Rscript search_adapter_mascot.R <CSV or folder> <out.spectra> [fixed mod name, default "Carbamidomethyl (C)"]
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


# Mascot 修饰名 → isoPSV 修饰名（基础名；括号内残基说明弃用，位点残基由序列给出）
# Mascot mod names -> isoPSV mod names (base name only; the parenthetical
# residue spec is dropped, the site residue comes from the sequence)
MASCOT_MOD_NAME <- c(
  "Oxidation" = "Oxidation", "Carbamidomethyl" = "Carbamidomethyl",
  "Deamidated" = "Deamidated", "Acetyl" = "Acetyl", "Phospho" = "Phospho",
  "Pyro-carbamidomethyl" = "Pyro-carbamidomethyl",
  "Gln->pyro-Glu" = "Gln->pyro-Glu", "Glu->pyro-Glu" = "Glu->pyro-Glu",
  "Propionamide" = "Propionamide", "Methyl" = "Methyl",
  "Dimethyl" = "Dimethyl", "Trimethyl" = "Trimethyl",
  "Methylthio" = "Methylthio", "Oxidation (M)" = "Oxidation",
  "Label:13C(6)" = "Label_13C(6)",
  "Label:13C(6)15N(4)" = "Label_13C(6)15N(4)",
  "Label:2H(4)" = "Label_2H(4)", "Label:18O(1)" = "Label_18O(1)"
)

mascot_mod_string <- function(seq0, var_mod, pos_str, summed_str = NULL, fixed_name = NULL) {
  out <- character(0)
  n <- nchar(seq0)
  names_l <- trimws(strsplit(var_mod, "[;,]", fixed = FALSE)[[1]])
  names_l <- names_l[nzchar(names_l)]
  add_mod <- function(pos, name) {
    base <- trimws(sub("\\(.*$", "", trimws(name)))   # 去括号说明后 trim
    if (base %in% names(MASCOT_MOD_NAME)) base <- MASCOT_MOD_NAME[[base]]
    if (grepl("N-term", name, fixed = TRUE)) {
      out <<- c(out, paste0("0,", base, "[AnyN-term];"))
    } else if (grepl("C-term", name, fixed = TRUE)) {
      out <<- c(out, paste0(n + 1L, ",", base, "[AnyC-term];"))
    } else {
      site <- if (pos >= 1L && pos <= n) substr(seq0, pos, pos) else "X"
      out <<- c(out, paste0(pos, ",", base, "[", site, "];"))
    }
  }
  dig <- strsplit(if (nzchar(pos_str)) pos_str else "0", "", fixed = TRUE)[[1]]
  dig2 <- if (!is.null(summed_str) && nzchar(summed_str))
    strsplit(summed_str, "", fixed = TRUE)[[1]] else rep("0", max(1L, n))
  if (length(dig) != n || length(dig2) != n) {
    cat("[Mascot 适配] 修饰位置串长度与序列不符（", length(dig), "/", length(dig2),
        " vs ", n, "），该行按无修饰处理\n", sep = "")
    dig <- rep("0", n); dig2 <- rep("0", n)
  }
  for (i in seq_len(n)) {
    d <- suppressWarnings(as.integer(dig[i]))
    if (!is.na(d) && d >= 1L && d <= length(names_l))
      add_mod(i, names_l[d])
    d2 <- suppressWarnings(as.integer(dig2[i]))
    if (!is.na(d2) && d2 >= 1L && d2 <= length(names_l))
      add_mod(i, names_l[d2])
  }
  # 固定修饰注入（默认 Carbamidomethyl (C)）
  if (!is.null(fixed_name) && nzchar(fixed_name) && tolower(fixed_name) != "none") {
    fx <- trimws(sub("\\(.*$", "", trimws(fixed_name)))
    res <- regmatches(fixed_name, regexec("\\(([A-Za-z]+)\\)", fixed_name))[[1]]
    res <- if (length(res) == 2) res[2] else "C"
    for (pos in which(strsplit(seq0, "", fixed = TRUE)[[1]] %in%
                      strsplit(res, "", fixed = TRUE)[[1]])) {
      out <- c(out, paste0(pos, ",", fx, "[", substr(seq0, pos, pos), "];"))
    }
  }
  paste(sort(unique(out)), collapse = "")
}

convert_mascot_csv <- function(csv_file, out_file, fixed_name = "Carbamidomethyl (C)") {
  csv_file <- resolve_adapter_input(csv_file, "\\.csv$")
  # P3-5（2026-09-25）：编码容错读取（UTF-8 BOM/UTF-8/本地编码）
  d <- read_csv_tolerant(csv_file, stringsAsFactors = FALSE, check.names = FALSE,
                         quote = "\"", comment.char = "", fill = TRUE,
                         na.strings = character(0), colClasses = "character")
  need <- c("pep_scan_title", "pep_seq")
  if (!all(need %in% names(d))) {
    stop("Mascot CSV 缺列：", paste(setdiff(need, names(d)), collapse = ", "),
         "\n请在 Mascot 导出时勾选 Query title 与 Peptide sequence（以及 Variable modifications / Modifications）\n",
         "missing columns; tick Query title, Peptide sequence and Variable modifications in the Mascot export")
  }
  # top-1：pep_rank 存在时取 1；缺省时按每谱标题取 pep_score 最高一行
  # 全候选模式（ISOPSV_CAND_ALL=1）：跳过 rank 过滤（保留全部候选，文件序）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    if ("pep_rank" %in% names(d)) {
      rk <- suppressWarnings(as.integer(d[["pep_rank"]]))
      keep <- !is.na(rk) & rk == 1L
      if (any(keep)) d <- d[keep, , drop = FALSE]
    }
  }
  d <- d[nzchar(d[["pep_seq"]]) & nzchar(d[["pep_scan_title"]]), , drop = FALSE]
  if ("pep_score" %in% names(d)) {
    sc <- suppressWarnings(as.numeric(d[["pep_score"]]))
    sc[is.na(sc)] <- -Inf
    d <- d[order(sc, decreasing = TRUE), , drop = FALSE]
  }
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[!duplicated(d[["pep_scan_title"]]), , drop = FALSE]
  }
  seqs <- d[["pep_seq"]]
  vmod <- if ("pep_var_mod" %in% names(d)) d[["pep_var_mod"]] else rep("", nrow(d))
  vpos <- if ("pep_var_mod_pos" %in% names(d)) d[["pep_var_mod_pos"]] else rep("", nrow(d))
  spos <- if ("pep_summed_mod_pos" %in% names(d)) d[["pep_summed_mod_pos"]] else rep("", nrow(d))
  mods <- vapply(seq_len(nrow(d)), function(i)
    mascot_mod_string(seqs[i], vmod[i], vpos[i], spos[i], fixed_name),
    character(1))
  prot <- if ("prot_acc" %in% names(d)) {
    vapply(strsplit(d[["prot_acc"]], ";", fixed = TRUE), function(ps) {
      ps <- trimws(ps); ps <- ps[nzchar(ps)]
      if (length(ps) == 0) "" else paste0(paste(ps, collapse = "/"), "/")
    }, character(1))
  } else rep("", nrow(d))
  charge <- if ("pep_exp_z" %in% names(d)) d[["pep_exp_z"]] else rep("", nrow(d))
  out <- data.frame(File_Name = d[["pep_scan_title"]], Sequence = seqs,
                    Modification = mods, Proteins = prot, Charge = charge,
                    stringsAsFactors = FALSE, check.names = FALSE)
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Mascot CSV -> %s：%d 行（top-1 每谱去重）\n", out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  fx <- if (length(a) >= 3) a[3] else "Carbamidomethyl (C)"
  convert_mascot_csv(a[1], a[2], fx)
}
