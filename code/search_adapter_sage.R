# =============================================================================
# search_adapter_sage.R - Sage results.sage.tsv + scan 映射表 -> isoPSV 交换格式
# search_adapter_sage.R - Sage TSV + scan map -> isoPSV interchange format
#
# 多引擎版（2026-09-09）：把 Sage（0.15.x）的 results.sage.tsv 转换为 isoPSV
# 补充策略文件（File_Name/Sequence/Modification/Proteins/Charge）。
#   - Sage 只读 mzML：先用 msconvert（或内置转换器）把重编号 MGF 转 mzML；
#     scannr 列为 "index=N"（Sage 0.15.0-beta.2 实测为 "scan=N"，裸整数亦可），
#     N 对应重编号 MGF 的顺序 scan，按 map.tsv 回查原始标题；
#   - peptide 列内联修饰 "M[+15.9949]"（质量位移）；静态修饰（CAM C）不标注——
#     适配器为序列中每个 C 注入 Carbamidomethyl[C]；
#   - label 列 -1 为诱饵（Sage 内产 REV_ 蛋白），默认过滤。
# multi-engine edition (2026-09-09): converts the Sage (0.15.x) results.sage.tsv
# into an isoPSV supplementary-strategy file.
#   - Sage reads mzML only: convert the renumbered MGF with msconvert (or the
#     built-in converter) first; the scannr column is "index=N" where N is the
#     sequential scan of the renumbered MGF, joined back via map.tsv;
#   - the peptide column carries inline mods "M[+15.9949]" (mass shifts);
#     static mods (CAM C) are NOT annotated - the adapter injects
#     Carbamidomethyl[C] for every C;
#   - the label column is -1 for decoys (Sage-generated REV_ proteins), filtered
#     by default.
#
# 用法：Rscript search_adapter_sage.R <results.sage.tsv> <map.tsv> <输出.spectra>
# usage: Rscript search_adapter_sage.R <results.sage.tsv> <map.tsv> <out.spectra>
# =============================================================================
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- getwd()
if (length(script_arg) > 0 && !is.na(script_arg)) {
  cand <- dirname(sub("^--file=", "", script_arg))
  if (dir.exists(cand)) script_dir <- cand
}
# compat.R 与适配器同目录；逐候选探测，不依赖易腐蚀的 normalizePath（Bug 6 同源）
compat_src <- NULL
for (d in unique(c(script_dir, getwd()))) {
  p <- file.path(d, "compat.R")
  if (file.exists(p)) { compat_src <- p; break }
}
if (is.null(compat_src)) stop("compat.R not found next to the adapter")
source(compat_src)

SAGE_MOD_NAME <- c(
  "15.9949"  = "Oxidation",
  "0.984"    = "Deamidated",
  "0.984016" = "Deamidated",
  "57.02146" = "Carbamidomethyl",
  "57.021464" = "Carbamidomethyl",
  "42.0106"  = "Acetyl",
  "42.010565" = "Acetyl",
  "8.050220" = "K080",
  "8.014204" = "K602",
  "3.988142" = "R004",
  "4.025110" = "R040",
  "-17.0265" = "Gln->pyro-Glu",
  "-17.026549" = "Gln->pyro-Glu",
  "-18.0106" = "Glu->pyro-Glu",
  "-18.010565" = "Glu->pyro-Glu"
)

sage_mod_string <- function(pep) {
  # "YNKDDPAM[+15.9949]EEEESGWK" -> "8,Oxidation[M];"
  # 另为序列中每个 C 注入静态 Carbamidomethyl
  out <- character(0)
  clean <- gsub("\\[[^]]*\\]", "", pep)     # 去修饰标注得纯序列
  # 内联修饰
  m <- gregexpr("[A-Z]\\[-?[+0-9.]+\\]", pep)[[1]]
  if (m[1] != -1) {
    for (i in seq_along(m)) {
      tok <- substr(pep, m[i], m[i] + attr(m, "match.length")[i] - 1)
      rm <- regmatches(tok, regexec("^([A-Z])\\[([+-]?)([0-9.]+)\\]$", tok))[[1]]
      if (length(rm) == 0) next
      res <- rm[2]
      shift <- paste0(rm[3], rm[4]) %||% rm[4]
      shift <- as.numeric(sub("^\\+", "", shift))
      hit <- which(abs(shift - as.numeric(names(SAGE_MOD_NAME))) < 1e-3)
      if (length(hit) == 0) {
        cat("[Sage 适配] 未知修饰位移（忽略）：", shift, "\n")
        next
      }
      nm <- SAGE_MOD_NAME[[hit[1]]]
      pos <- nchar(gsub("\\[[^]]*\\]", "", substr(pep, 1, m[i])))  # 修饰残基位置
      out <- c(out, paste0(pos, ",", nm, "[", res, "];"))
    }
  }
  # 静态 CAM
  if (grepl("C", clean, fixed = TRUE)) {
    for (p in gregexpr("C", clean, fixed = TRUE)[[1]]) {
      out <- c(out, paste0(p, ",Carbamidomethyl[C];"))
    }
  }
  paste(sort(out), collapse = "")
}


# 文件夹解析（用户 2026-09-11）：输入填文件夹时按目标模式递归（<=4 层）发现
# 结果文件、多个候选取最新、找不到报错；文件路径原样通过。
# folder resolution (user 2026-09-11): a directory input is resolved by the
# target pattern (recursive depth <= 4, newest wins, clear error otherwise);
# a file path passes through unchanged.

convert_sage_tsv <- function(tsv_file, map_file, out_file, keep_decoy = FALSE) {
  tsv_file <- resolve_adapter_input(tsv_file, "results\\.sage\\.tsv$")
  d <- read.delim(tsv_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0), quote = "",
                  colClasses = "character")
  # P3（2026-09-26，工作站反馈：sage 适配器产物 0 行无诊断）——缺列时列出
  # 实际列名，Sage 版本差异（列改名）一眼可辨
  need_cols <- c("peptide", "proteins", "scannr", "charge", "label", "rank")
  miss <- setdiff(need_cols, colnames(d))
  if (length(miss) > 0)
    stop("results.sage.tsv 缺列：", paste(miss, collapse = ","),
         "（实际列：", paste(head(colnames(d), 20), collapse = ","), "）",
         "\nresults.sage.tsv missing columns: ", paste(miss, collapse = ","),
         " (actual: ", paste(head(colnames(d), 20), collapse = ","), ")")
  mp <- read.delim(map_file, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = character(0), quote = "",
                   colClasses = "character")
  stopifnot(all(c("new_scan", "original_title") %in% colnames(mp)))
  n_raw <- nrow(d)
  # 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（引擎自身排序，文件序）；
  # 常规模式：每谱 top-1（rank==1 + 首行去重）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[d$rank == "1", , drop = FALSE]                       # 每谱 top-1
  }
  n_rank <- nrow(d)
  if (!keep_decoy) d <- d[d$label != "-1", , drop = FALSE]   # 过滤诱饵
  n_dec <- nrow(d)
  # v1.0（工作站根因定位）：Sage 0.15.0-beta.2 的 scannr 写 "scan=N"，
  # 老版本写 "index=N"——只剥 index= 会导致 map 命中恒 0。剥任意 ^[^=]*=
  # 前缀再取整数（裸数字不受影响）；v1.0 自查整改：剥后 trimws 与
  # "scan = 123" 这类等号带空格的形式对齐，map 键同样 trimws 归一
  d$scan_num <- trimws(sub("^[^=]*=\\s*", "", d$scannr))
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[!duplicated(d$scan_num), , drop = FALSE]
  }
  fn <- mp$original_title[match(d$scan_num, mp$new_scan)]
  n_map <- sum(!is.na(fn))
  seqs <- gsub("\\[[^]]*\\]", "", d$peptide)
  keep <- !is.na(fn) & nzchar(seqs)
  d <- d[keep, , drop = FALSE]; fn <- fn[keep]; seqs <- seqs[keep]
  mods <- vapply(d$peptide, sage_mod_string, character(1), USE.NAMES = FALSE)
  prot <- vapply(strsplit(d$proteins, ";", fixed = TRUE), function(ps) {
    acc <- sub("\\s.*$", "", trimws(ps))
    acc <- acc[nzchar(acc)]
    if (length(acc) == 0) "" else paste0(paste(acc, collapse = "/"), "/")
  }, character(1))
  out <- data.frame(File_Name = fn, Sequence = seqs,
                    Modification = mods, Proteins = prot,
                    Charge = d$charge,
                    stringsAsFactors = FALSE, check.names = FALSE)
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Sage TSV -> %s：%d 行（原始 %d | top-1 %d | 去诱饵 %d | map 命中 %d）\n",
              out_file, nrow(out), n_raw, n_rank, n_dec, n_map))
  if (nrow(out) == 0 && n_raw > 0)
    cat("[Sage 适配] 警告：结果 0 行但输入有 ", n_raw,
        " 行——核对 scannr 前缀（Sage 0.15 实测为 scan=N，老版本为 index=N，裸数字亦可）与 map 键，或检查结果是否全为诱饵/非 rank=1\n",
        sep = "")
  invisible(out)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 3) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_sage_tsv(a[1], a[2], a[3])
}
