
# 文件夹解析（用户 2026-09-11）：输入填文件夹时按目标模式递归（<=4 层）发现
# 结果文件、多个候选取最新、找不到报错；文件路径原样通过。
# folder resolution (user 2026-09-11): a directory input is resolved by the
# target pattern (recursive depth <= 4, newest wins, clear error otherwise);
# a file path passes through unchanged.

# search_adapter_metamorpheus.R - MetaMorpheus AllPeptides.psmtsv -> isoPSV interchange
# 用法：Rscript search_adapter_metamorpheus.R <psmtsv 或文件夹> <输出.spectra> [qval] [map.tsv]
# usage: Rscript search_adapter_metamorpheus.R <psmtsv or folder> <out.spectra> [qval] [map.tsv]
# map.tsv（可选）：重编号映射（new_scan→original_title），MetaMorpheus 搜的是重编号
# mzML，其 Scan Number 为顺序编号——有映射时 File_Name 还原为 pFind 原始标题
# （与 MSFragger/Sage/Comet 适配器同口径）。
convert_metamorpheus_psm <- function(input_file, out_file, qval_cutoff = 0.01,
                                     map_file = NULL) {
  input_file <- resolve_adapter_input(input_file, "AllPeptides\\.psmtsv$")
  if (!file.exists(input_file)) stop("file not found: ", input_file)
  mp <- NULL
  if (!is.null(map_file) && file.exists(map_file)) {
    mp <- read.delim(map_file, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, na.strings = character(0), quote = "",
                     colClasses = "character")
  }

  d <- read.delim(input_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0), quote = "")
  cat("read:", nrow(d), "rows\n")
  d <- d[d[["Decoy/Contaminant/Target"]] == "T" &
         suppressWarnings(as.numeric(d[["QValue"]])) <= qval_cutoff, ]
  cat("filtered:", nrow(d), "rows\n")

# Sequence: use Base Sequence column (already mod-stripped)
seqs <- d[["Base Sequence"]]

# Modification: parse Full Sequence for mod positions relative to Base.
# 位点＝修饰 token 前已消费的 Base 残基数（1 起，前缀计数法）——旧实现的
# 累计逻辑在连续 token 时会错位（33 肽写出 29/64 位点致 pDeep2 越界崩溃，
# 2026-09-13 修复）。
extract_mods <- function(full, base) {
  if (is.na(full) || !nzchar(full)) return("")
  if (!grepl("\\[", full)) return("")
  # perl=TRUE：R 4.5.1 UCRT 的 TRE 引擎对 "[^\]]" 否定转义类匹配失败
  # （实测 FALSE），PCRE 正常——此处必须用 perl（2026-09-13 修复）
  tok_pos <- gregexpr("\\[[^\\]]+\\]", full, perl = TRUE)[[1]]
  tokens <- regmatches(full, gregexpr("\\[[^\\]]+\\]", full, perl = TRUE))[[1]]
  if (length(tokens) == 0) return("")
  mods <- c()
  for (k in seq_along(tokens)) {
    prefix <- substr(full, 1, tok_pos[k] - 1)
    # perl=TRUE 必须：R 4.5.1 UCRT 的 TRE 引擎对 "[^\]]" 否定类匹配失败
    aa_before <- nchar(gsub("\\[[^\\]]*\\]", "", prefix, perl = TRUE))
    pos <- aa_before + 1L
    clean <- gsub("^\\[|\\]$", "", tokens[k])
    clean <- sub("^Common Variable:", "", clean)
    clean <- sub("^Common Fixed:", "", clean)
    pf <- NULL
    if (grepl("Oxidation on M", clean, fixed = TRUE)) {
      pf <- paste0(pos, ",Oxidation[M];")
    } else if (grepl("Acetyl on Protein N", clean, fixed = TRUE)) {
      pf <- "0,Acetyl[AnyN-term];"
    } else if (grepl("Deamidated on N", clean, fixed = TRUE)) {
      pf <- paste0(pos, ",Deamidated[N];")
    } else if (grepl("Deamidated on Q", clean, fixed = TRUE)) {
      pf <- paste0(pos, ",Deamidated[Q];")
    }
    # Carbamidomethyl = fixed mod, skip
    if (!is.null(pf)) mods <- c(mods, pf)
  }
  paste(mods, collapse = "")
}
mods <- mapply(extract_mods, d[["Full Sequence"]], seqs, USE.NAMES = FALSE)

# title：有映射时还原为 pFind 原始标题；无映射时用 mzML 文件名＋顺序扫描号拼接
# title: restore the original pFind title via the map when given; otherwise join
# the mzML file name with the sequential scan number
if (!is.null(mp)) {
  sc <- as.character(d[["Scan Number"]])
  fn0 <- mp$original_title[match(sc, mp$new_scan)]
  fn0[is.na(fn0) | fn0 == ""] <- paste0(d[["File Name"]], ".", sc, ".", sc,
      ".", as.integer(as.numeric(d[["Precursor Charge"]])), ".0.dta")[is.na(fn0) | fn0 == ""]
  fn <- fn0
} else {
  fn <- paste0(d[["File Name"]], ".", d[["Scan Number"]], ".", d[["Scan Number"]],
               ".", as.integer(as.numeric(d[["Precursor Charge"]])), ".0.dta")
}
# protein accession only
prots <- sub("\\s.*$", "/", d[["Accession"]])
prots[is.na(prots) | prots == "/"] <- "unknown/"

out <- data.frame(File_Name = fn, Sequence = seqs, Modification = mods,
                  Proteins = prots, Charge = as.integer(as.numeric(d[["Precursor Charge"]])),
                  stringsAsFactors = FALSE, check.names = FALSE)
# 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（文件序）
if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
  out <- out[!duplicated(out$File_Name), ]
}
write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("MetaMorpheus ->", out_file, ":", nrow(out), "rows\n")
invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  qv <- if (length(a) >= 3) suppressWarnings(as.numeric(a[3])) else 0.01
  if (is.na(qv)) qv <- 0.01
  mp <- if (length(a) >= 4) a[4] else NULL
  convert_metamorpheus_psm(a[1], a[2], qv, mp)
}
