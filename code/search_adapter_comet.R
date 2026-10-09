# =============================================================================
# search_adapter_comet.R - Comet pepXML/txt + scan 映射表 -> isoPSV 交换格式
# search_adapter_comet.R - Comet txt + scan map -> isoPSV interchange format
#
# 多引擎版（2026-09-09）：把 Comet（2026.02 rev.2）的 tab 分隔 txt 输出转换为
# isoPSV 补充策略文件（File_Name/Sequence/Modification/Proteins/Charge）。
#   - Comet 原生支持 MGF，但输出只有 scan 号；pParse 标题同 scan 多价态/拆分，
#     故与 MSFragger 相同先用 _renumber_mgf.py 重编号并留 map.tsv，本适配器按
#     scan 回查原始标题；
#   - txt 结构：第 1 行版本元信息、第 2 行列表头、第 3 行起数据；peptide 列带
#     侧翼 "R.PEPTIDE.K"，可变修饰内联 "[15.9949]"、N 端记 "n[42.0106]"；
#     静态修饰（CAM C）不标注——适配器为序列中每个 C 注入 Carbamidomethyl[C]；
#   - 质量位移→pFind 名与 MSFragger 适配器同一张表。
# multi-engine edition (2026-09-09): converts the Comet (2026.02 rev.2)
# tab-delimited txt output into an isoPSV supplementary-strategy file.
#   - Comet reads MGF natively but reports only scan numbers; pParse titles carry
#     duplicate scans across charges/splits, so renumber the MGF first with
#     _renumber_mgf.py (unique sequential scans + map.tsv) exactly as for
#     MSFragger; this adapter joins the scan back to the original title;
#   - txt layout: row 1 version meta, row 2 column header, row 3+ data; the
#     peptide column is flanked "R.PEPTIDE.K" with inline variable mods
#     "[15.9949]" and the N-term as "n[42.0106]"; static mods (CAM C) are NOT
#     annotated - the adapter injects Carbamidomethyl[C] for every C;
#   - the mass-shift -> pFind name table is shared with the MSFragger adapter.
#
# 用法：Rscript search_adapter_comet.R <comet.txt> <map.tsv> <输出.spectra>
# usage: Rscript search_adapter_comet.R <comet.txt> <map.tsv> <out.spectra>
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

CMT_MOD_NAME <- c(
  "15.9949"  = "Oxidation",
  "0.9840"   = "Deamidated",
  "0.984016" = "Deamidated",
  "57.02146" = "Carbamidomethyl",
  "42.0106"  = "Acetyl",
  "8.050220" = "K080",
  "8.014204" = "K602",
  "3.988142" = "R004",
  "4.025110" = "R040",
  "-17.0265" = "Gln->pyro-Glu",
  "-18.0106" = "Glu->pyro-Glu",
  "21.9819"  = "Cation_Na"
)
CMT_TERM_SITE <- c("42.0106" = "AnyN-term",
                   "-17.0265" = "AnyN-termQ",
                   "-18.0106" = "AnyN-termE")

comet_mod_string <- function(flanged, seq) {
  # "R.n[42.0106]AE[-18.0106]ADERTK.A" -> "0,Acetyl[AnyN-term];2,Glu->pyro-Glu[E];"
  # 另为序列中每个 C 注入静态 Carbamidomethyl（Comet 不标注静态修饰）
  # plus inject static Carbamidomethyl at every C (Comet omits static mods)
  core <- substr(flanged, 3, nchar(flanged) - 2)     # 去两侧翼 "R." ".A"
  out <- character(0)
  pos <- 0L
  i <- 1L
  ch <- strsplit(core, "", fixed = TRUE)[[1]]
  while (i <= length(ch)) {
    cc <- ch[i]
    if (grepl("[A-Z]", cc)) {
      pos <- pos + 1L
      i <- i + 1L
      if (i <= length(ch) && ch[i] == "[") {
        j <- i + 1L
        buf <- ""
        while (j <= length(ch) && ch[j] != "]") { buf <- paste0(buf, ch[j]); j <- j + 1L }
        hit <- which(abs(as.numeric(buf) - as.numeric(names(CMT_MOD_NAME))) < 1e-3)
        if (length(hit) > 0) {
          nm <- CMT_MOD_NAME[[hit[1]]]
          key <- names(CMT_MOD_NAME)[hit[1]]
          if (pos == 0L || key %in% names(CMT_TERM_SITE)) {
            site <- CMT_TERM_SITE[[key]]
            out <- c(out, paste0("0,", nm, "[", site, "];"))
          } else {
            site <- substr(seq, pos, pos)
            out <- c(out, paste0(pos, ",", nm, "[", site, "];"))
          }
        }
        i <- j + 1L
      }
    } else if (cc == "n" && i + 1L <= length(ch) && ch[i + 1L] == "[") {
      # 肽 N 端修饰 "n[...]"（尚未计入残基，位置 0）
      j <- i + 2L
      buf <- ""
      while (j <= length(ch) && ch[j] != "]") { buf <- paste0(buf, ch[j]); j <- j + 1L }
      hit <- which(abs(as.numeric(buf) - as.numeric(names(CMT_MOD_NAME))) < 1e-3)
      if (length(hit) > 0) {
        nm <- CMT_MOD_NAME[[hit[1]]]
        key <- names(CMT_MOD_NAME)[hit[1]]
        site <- if (key %in% names(CMT_TERM_SITE)) CMT_TERM_SITE[[key]] else "AnyN-term"
        out <- c(out, paste0("0,", nm, "[", site, "];"))
      }
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  # 静态 CAM：序列中每个 C
  if (grepl("C", seq, fixed = TRUE)) {
    for (p in gregexpr("C", seq, fixed = TRUE)[[1]]) {
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

convert_comet_txt <- function(txt_file, map_file, out_file) {
  txt_file <- resolve_adapter_input(txt_file, "(pep\\.tsv$|\\.txt$)")
  d <- read.delim(txt_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0), quote = "",
                  colClasses = "character", header = FALSE, skip = 2)
  stopifnot(ncol(d) >= 16)
  # V1=scan V3=charge V12=纯序列 V13=带侧翼肽（修饰标注） V16=蛋白
  # V1=scan V3=charge V12=plain sequence V13=flanked peptide (mod tags) V16=proteins
  names(d)[c(1, 3, 12, 13, 16)] <- c("scan", "charge", "seq", "peptide", "proteins")
  mp <- read.delim(map_file, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = character(0), quote = "",
                   colClasses = "character")
  stopifnot(all(c("new_scan", "original_title") %in% colnames(mp)))
  # 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（文件序＝引擎 rank 序）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[!duplicated(d$scan), , drop = FALSE]                 # 每谱 top-1
  }
  seqs <- d$seq
  fn <- mp$original_title[match(d$scan, mp$new_scan)]
  keep <- !is.na(fn) & nzchar(seqs)
  d <- d[keep, , drop = FALSE]; fn <- fn[keep]; seqs <- seqs[keep]
  mods <- mapply(comet_mod_string, d$peptide, seqs, USE.NAMES = FALSE)
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
  cat(sprintf("Comet txt -> %s：%d 行\n", out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 3) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_comet_txt(a[1], a[2], a[3])
}
