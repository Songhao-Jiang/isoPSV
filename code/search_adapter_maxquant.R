# =============================================================================
# search_adapter_maxquant.R - MaxQuant msms.txt -> isoPSV 交换格式转换器
# search_adapter_maxquant.R - MaxQuant msms.txt -> isoPSV interchange-format converter
#
# 多引擎版把 MaxQuant 的 msms.txt 转换为 isoPSV 补充策略文件
# （最少三列 File_Name/Sequence/Modification，另带 Proteins/Charge）。
#   - 过滤：Reverse（诱饵）与 Potential contaminant 排除、PEP <= 阈值；
#   - 修饰映射：MaxQuant Modified sequence 的 (ox)/(ac)/(de)/(Gln->pyro-Glu) 等
#     标注转换为 pFind 修饰串（位置 1 起、N 端修饰记 0、[AA]/[AnyN-term] 形式）；
#   - 谱图号映射：Raw file + Scan number 按 pFind 标题约定拼接
#     "<raw>.<scan>.<scan>.<charge>.0.dta"（阶段 B 有 raw+scan 前缀兜底匹配，
#     价态/拆分序号差异不影响定位）。
# multi-engine edition converts MaxQuant msms.txt into an isoPSV
# supplementary-strategy file (at least File_Name/Sequence/Modification, plus
# Proteins/Charge).
#   - filters: exclude Reverse (decoy) and Potential contaminant, PEP <= cutoff;
#   - modification mapping: the (ox)/(ac)/(de)/(Gln->pyro-Glu) tags of the
#     MaxQuant modified sequence become pFind modification strings (1-based
#     positions, N-terminal mods at 0, [AA]/[AnyN-term] forms);
#   - spectrum-id mapping: Raw file + Scan number join to the pFind title
#     convention "<raw>.<scan>.<scan>.<charge>.0.dta" (stage B falls back to a
#     raw+scan prefix match, so charge/fraction differences do not matter).
#
# 用法：Rscript search_adapter_maxquant.R <msms.txt> <输出.spectra> [PEP阈值]
# usage: Rscript search_adapter_maxquant.R <msms.txt> <out.spectra> [PEP cutoff]
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

# MaxQuant Modified sequence 用修饰全名（如 "(Oxidation (M))"）；映射到 pFind 名
# the MaxQuant modified sequence uses full mod names (e.g. "(Oxidation (M))");
# map them to the pFind names
# P0-2补 MaxQuant 最常见的内联修饰全名（旧表仅 5 种，
# 未知名如 S(Phospho (ST)) 会进入逐残基循环、括号内字母被误计为残基致位点错位）
MQ_FULL <- c("Oxidation (M)" = "Oxidation",
             "Deamidation (NQ)" = "Deamidated",
             "Acetyl (Protein N-term)" = "Acetyl",
             "Gln->pyro-Glu (N-term Q)" = "Gln->pyro-Glu",
             "Glu->pyro-Glu (N-term E)" = "Glu->pyro-Glu",
             "Phospho (STY)" = "Phospho",
             "Methyl (KR)" = "Methyl",
             "Dimethyl (KR)" = "Dimethyl",
             "Trimethyl (K)" = "Trimethyl",
             "Carbamidomethyl (C)" = "Carbamidomethyl",
             # MiniMax P0-3MaxQuant 短标签（ac=乙酰化 N 端、
             # ox=氧化、ph=磷酸化、de=脱酰胺）——旧表只认全名，短标签被静默丢弃
             "ac" = "Acetyl",
             "ox" = "Oxidation",
             "ph" = "Phospho",
             "de" = "Deamidated")
MQ_TERM_SITE <- c("Acetyl (Protein N-term)" = "AnyN-term",
                  "Gln->pyro-Glu (N-term Q)" = "AnyN-termQ",
                  "Glu->pyro-Glu (N-term E)" = "AnyN-termE",
                  "ac" = "AnyN-term")

parse_modified_sequence <- function(mseq, seq) {
  # 返回 pFind 风格修饰串："pos,Name[site];..."；无法解析时返回空
  # Returns a pFind-style modification string "pos,Name[site];..."; empty on failure
  if (!nzchar(mseq) || mseq == seq) return("")
  mods <- character(0)
  us <- gregexpr("_", mseq, fixed = TRUE)[[1]]
  us <- us[us > 0]
  mid <- mseq
  if (length(us) > 0) {
    # 去掉首尾 "_" 外框（N 端修饰标签位于第一个 "_" 之后、C 端修饰在最后一个 "_" 之后）
    # strip the surrounding "_" frame (N-terminal tags follow the first "_",
    # C-terminal mods follow the last "_")
    mid <- substr(mseq, if (us[1] >= 1) us[1] + 1 else 1,
                  if (us[length(us)] > us[1]) us[length(us)] - 1 else nchar(mseq))
  }
  # 先吞掉开头的 N 端修饰标签（"(Acetyl (Protein N-term))" 等，含大写字母，
  # 不能按残基解析）——按已知全名锚定匹配
  # first consume the leading N-terminal tags ("(Acetyl (Protein N-term))"
  # etc., they contain uppercase letters and must not be parsed as residues) -
  # anchor on the known full names
  repeat {
    if (!nzchar(mid) || substr(mid, 1, 1) != "(") break
    hit <- FALSE
    for (nm in names(MQ_FULL)) {
      if (substr(mid, 2, 1 + nchar(nm)) == nm) {
        site <- if (nm %in% names(MQ_TERM_SITE)) MQ_TERM_SITE[[nm]] else ""
        mods <- c(mods, paste0("0,", MQ_FULL[[nm]],
                               if (nzchar(site)) paste0("[", site, "]") else "",
                               ";"))
        mid <- substr(mid, nchar(nm) + 3L, nchar(mid))
        hit <- TRUE
        break
      }
    }
    if (!hit) {
      # B-4前导未知 N 端标签（不在 MQ_FULL 内）整体跳过配对括号，
      # 否则 "Protein N-term" 里的 P/N 等大写字母会被逐字符循环误计为残基。
      # N-19这里须 continue 而非 break——连续两个未知前导
      # 标签时 break 只吞掉第一个，第二个会漏进逐残基解析导致位点漂移。
      dep2 <- 1L; k <- 2L
      while (k <= nchar(mid) && dep2 > 0) {
        cc2 <- substr(mid, k, k)
        if (cc2 == "(") dep2 <- dep2 + 1L else if (cc2 == ")") dep2 <- dep2 - 1L
        k <- k + 1L
      }
      mid <- substr(mid, k, nchar(mid))
      next
    }
  }
  # 逐残基解析内联修饰标签
  # parse the inline tags residue by residue
  pos <- 0L
  ch <- strsplit(mid, "", fixed = TRUE)[[1]]
  i <- 1L
  while (i <= length(ch)) {
    cc <- ch[i]
    if (grepl("[A-Z]", cc)) {
      pos <- pos + 1L
      i <- i + 1L
      if (i <= length(ch) && ch[i] == "(") {
        hit <- FALSE
        for (nm in names(MQ_FULL)) {
          if (substr(mid, i + 1L, i + nchar(nm)) == nm) {
            if (nm %in% c("Gln->pyro-Glu (N-term Q)", "Glu->pyro-Glu (N-term E)") && pos == 1L) {
              mods <- c(mods, paste0("0,", MQ_FULL[[nm]], "[",
                                     MQ_TERM_SITE[[nm]], "];"))
            } else {
              site <- substr(seq, pos, pos)
              mods <- c(mods, paste0(pos, ",", MQ_FULL[[nm]], "[", site, "];"))
            }
            i <- i + nchar(nm) + 2L
            hit <- TRUE
            break
          }
        }
        if (!hit) {
          # 未知名内联修饰：整体跳过配对括号内容（旧实现只跳过 "(" 一个字符，
          # 括号内大写字母被误计为残基，其后所有修饰位点系统性错位）
          dep <- 1L; j <- i + 1L
          while (j <= length(ch) && dep > 0) {
            if (ch[j] == "(") dep <- dep + 1L
            else if (ch[j] == ")") dep <- dep - 1L
            j <- j + 1L
          }
          cat("未知 MaxQuant 内联修饰，跳过括号内容：", substr(mid, i, min(j - 1L, nchar(mid))), "
")
          i <- j
        }
      }
    } else {
      i <- i + 1L
    }
  }
  paste(mods, collapse = "")
}


# 文件夹解析输入填文件夹时按目标模式递归（<=4 层）发现
# 结果文件、多个候选取最新、找不到报错；文件路径原样通过。
# folder resolution a directory input is resolved by the
# target pattern (recursive depth <= 4, newest wins, clear error otherwise);
# a file path passes through unchanged.

convert_msms <- function(msms_file, out_file, pep_cutoff = 0.01) {
  msms_file <- resolve_adapter_input(msms_file, "msms\\.txt$")
  d <- read.delim(msms_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0), quote = "")
  names(d)[1] <- sub("^\ufeff", "", names(d)[1])
  # MaxQuant 2.8 把诱饵列改名为 Decoy（早期实现为 Reverse），兼容两种表头
  # MaxQuant 2.8 renamed the decoy column to Decoy (older versions used
  # Reverse); accept both headers
  if ("Reverse" %in% names(d)) names(d)[names(d) == "Reverse"] <- "Decoy"
  need <- c("Raw file", "Scan number", "Sequence", "Modifications",
            "Modified sequence", "Proteins", "Charge", "PEP", "Decoy")
  if (!all(need %in% names(d))) {
    stop("msms.txt 缺列：", paste(setdiff(need, names(d)), collapse = ", "))
  }
  # 过滤：非诱饵、非污染、PEP <= 阈值
  # filters: keep non-decoy, non-contaminant rows with PEP <= cutoff
  rev_ok <- !(d[["Decoy"]] %in% c("+"))
  cont_ok <- if ("Potential contaminant" %in% names(d))
    !(d[["Potential contaminant"]] %in% c("+")) else rep(TRUE, nrow(d))
  pep_ok <- suppressWarnings(as.numeric(d[["PEP"]])) <= pep_cutoff
  keep <- rev_ok & cont_ok & pep_ok
  keep[is.na(keep)] <- FALSE
  d <- d[keep,, drop = FALSE]
  # 同一谱图多条 PSM（如不同价态同扫描号重复）时保留 PEP 最佳一行——交换格式
  # 每谱每策略一行
  # several PSMs per spectrum (e.g. duplicated scan numbers across charges)
  # keep the best-PEP row - the interchange format holds one row per spectrum
  # per strategy
  d$pep_num <- suppressWarnings(as.numeric(d[["PEP"]]))
  d <- d[order(d$pep_num, na.last = TRUE),, drop = FALSE]
  # 谱图号：raw.scan.scan.charge.0.dta（阶段 B 有 raw+scan 前缀兜底）
  # spectrum id: raw.scan.scan.charge.0.dta (stage B falls back to raw+scan)
  # 修复fn 必须在排序之后构建。此前 fn 在排序前构建、排序后未随之
  # 重排，duplicated 掩码把 File_Name 与肽段行错位配对（VEVFEY 真实在 scan 23928
  # 却被写到 26377 名下）。先排序、后建 fn、再按 fn 去重，行对齐得以保持。
  # FIX fn must be built AFTER the sort. Previously fn was built
  # before the PEP sort and not reordered with it, so the duplicated() mask
  # paired File_Names with the wrong peptide rows (VEVFEY truly at scan 23928
  # was written under 26377). Sort first, then build fn, then dedup on fn -
  # keeping rows aligned.
  fn <- paste0(d[["Raw file"]], ".", d[["Scan number"]], ".", d[["Scan number"]],
               ".", d[["Charge"]], ".0.dta")
  # 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（PEP 升序＝引擎 rank 序）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[!duplicated(fn),, drop = FALSE]
    fn <- fn[!duplicated(fn)]
  }
  mods <- vapply(seq_len(nrow(d)), function(i)
    parse_modified_sequence(d[["Modified sequence"]][i], d[["Sequence"]][i]),
    character(1))
  # Proteins：MaxQuant 以 ";" 分隔的 ID → pFind "ACC1/ACC2/" 拼接
  # Proteins: MaxQuant ";"-separated IDs -> pFind "ACC1/ACC2/" concatenation
  prot <- vapply(strsplit(d[["Proteins"]], ";", fixed = TRUE), function(ps) {
    ps <- trimws(ps); ps <- ps[nzchar(ps)]
    if (length(ps) == 0) "" else paste0(paste(ps, collapse = "/"), "/")
  }, character(1))
  out <- data.frame(File_Name = fn, Sequence = d[["Sequence"]],
                    Modification = mods, Proteins = prot,
                    Charge = d[["Charge"]],
                    stringsAsFactors = FALSE, check.names = FALSE)
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("MaxQuant msms.txt -> %s：%d 行（过滤后）\n", out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  pep <- if (length(a) >= 3) as.numeric(a[3]) else 0.01
  convert_msms(a[1], a[2], pep)
}
