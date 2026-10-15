# =============================================================================
# search_adapter_msfragger.R - MSFragger TSV + scan 映射表 -> isoPSV 交换格式
# search_adapter_msfragger.R - MSFragger TSV + scan map -> isoPSV interchange format
#
# 多引擎版把 MSFragger（4.4.1）的 TSV 输出转换为 isoPSV 补充
# 策略文件（File_Name/Sequence/Modification/Proteins/Charge）。
#   - MSFragger 可直接搜 MGF，但 pParse 标题（scan.scan.charge.fraction.dta）
#     的同 scan 多价态/拆分触发 "duplicated scan numbers"——先用 _renumber_mgf.py
#     把 MGF 重编号为唯一顺序 scan 并产出 map.tsv（new_scan→原始标题），本适配器
#     按 scannum 列回查原始标题（pFind File_Name 口径，阶段 B 有 raw+scan 兜底）；
#   - modification_info 编码 "20M(15.9949), 6M(15.9949)"（pos+残基+质量位移，
#     固定修饰也在内），按质量位移映射表转 pFind 修饰串（pos,Name[site]）；
#   - 蛋白列 ";" 分隔且带描述，取首个空格前的 accession，"ACC1/ACC2/" 拼接。
# multi-engine edition converts the MSFragger (4.4.1) TSV output
# into an isoPSV supplementary-strategy file.
#   - MSFragger searches MGF directly, but pParse titles (scan.scan.charge.
#     fraction.dta) carry duplicate scan numbers across charges/splits and trip
#     the "duplicated scan numbers" check - renumber the MGF first with
#     _renumber_mgf.py (unique sequential scans + map.tsv of new_scan ->
#     original title); this adapter joins the TSV scannum back to the original
#     title (the pFind File_Name convention; stage B falls back to raw+scan);
#   - modification_info "20M(15.9949), 6M(15.9949)" (position + residue + mass
#     shift, fixed mods included) becomes a pFind string via the shift map
#     (pos,Name[site];);
#   - the proteins column is ";"/description-bearing - take the accession
#     before the first space and join as "ACC1/ACC2/".
#
# 用法：Rscript search_adapter_msfragger.R <msfragger.tsv> <map.tsv> <输出.spectra>
# usage: Rscript search_adapter_msfragger.R <msfragger.tsv> <map.tsv> <out.spectra>
# =============================================================================
args0 <- commandArgs(FALSE)
script_arg <- args0[grepl("--file=", args0)][1]
script_dir <- getwd()
if (length(script_arg) > 0 && !is.na(script_arg)) {
  cand <- dirname(sub("^--file=", "", script_arg))
  # normalizePath 在中文路径＋R 4.5 C.UTF-8 下会产出腐蚀路径（Bug 6 同源），
  # 只接受真实存在的候选目录
  # normalizePath corrupts CJK paths under R 4.5 + C.UTF-8 (same root as
  # Bug 6): accept the candidate only when it really exists
  if (dir.exists(cand)) script_dir <- cand
}
# compat.R 与适配器同目录；逐候选（脚本目录、工作目录）探测，不依赖易腐蚀的路径转换
# compat.R lives next to the adapter; probe the candidates instead of trusting
# the corruption-prone path conversion
compat_src <- NULL
for (d in unique(c(script_dir, getwd()))) {
  p <- file.path(d, "compat.R")
  if (file.exists(p)) { compat_src <- p; break }
}
if (is.null(compat_src)) stop("compat.R not found next to the adapter")
source(compat_src)

# 质量位移 → pFind 修饰名（与 closed/open isoPSV 参数一致）
# mass shift -> pFind modification name (matching the isoPSV param sets)
MSF_MOD_NAME <- c(
  "15.9949"  = "Oxidation",
  "0.984016" = "Deamidated",
  "57.02146" = "Carbamidomethyl",
  "42.0106"  = "Acetyl",
  "-17.0265" = "Gln->pyro-Glu",
  "-18.0106" = "Glu->pyro-Glu",
  "21.9819"  = "Cation_Na",
  "8.050220" = "K080",
  "8.014204" = "K602",
  "3.988142" = "R004",
  "4.025110" = "R040"
)
# N 端类修饰的位点形式（位移命中且 token 无残基/为 n 时使用）
# site forms for the N-terminal-like mods (used when the token has no residue)
MSF_TERM_SITE <- c("42.0106" = "AnyN-term",
                   "-17.0265" = "AnyN-termQ",
                   "-18.0106" = "AnyN-termE")

# 独立 AA 质量/三字母表（适配器独立运行，不依赖主表导出模块）。
# Standalone AA mass / three-letter tables (the adapter runs standalone).
MSF_AA_MASS <- c(G = 57.02146, A = 71.03711, S = 87.03203, P = 97.05276,
                 V = 99.06841, T = 101.04768, C = 103.00919, L = 113.08406,
                 I = 113.08406, N = 114.04293, D = 115.02694, Q = 128.05858,
                 K = 128.09496, E = 129.04259, M = 131.04049, H = 137.05891,
                 F = 147.06841, R = 156.10111, Y = 163.06333, W = 186.07931)
MSF_AA3 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q",
             Glu = "E", Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K",
             Met = "M", Phe = "F", Pro = "P", Ser = "S", Thr = "T", Trp = "W",
             Tyr = "Y", Val = "V")

# 未知位移 → 残基替换 X->Y 解释
# 位移 delta 在位置残基 X 上若满足 mass(X) + delta ∈ 某标准氨基酸质量 ±0.05 Da，
# 视为序列替换并输出三字母 "pos,X->Y"（I/L 同质量不可区分时按 pFind SAP 惯例
# 记 Xle）——与 isoPSV 的 pep_mh_mod 质量计算和 SAP 判定（目标 ∈ AA_names）
# 完全兼容；不满足时返回 ""（调用方按"未知修饰"丢弃，诚实放弃判定）。
# unknown shift -> residue substitution X->Y (g-slot open-search delta
# convention, settled ): when mass(X) + delta hits a standard AA
# mass within 0.05 Da the token becomes "pos,X->Y" in three-letter names
# (I/L both 113.08406 -> "Xle" per the pFind SAP convention), compatible with
# isoPSV's pep_mh_mod mass math and the SAP test (target in AA_names);
# otherwise "" is returned (the caller drops it - an honest non-judgement).
msf_delta_sub <- function(res, delta) {
  if (!res %in% names(MSF_AA_MASS)) return("")
  tm <- MSF_AA_MASS[[res]] + delta
  hits <- names(MSF_AA_MASS)[abs(MSF_AA_MASS - tm) <= 0.05]
  if (length(hits) == 0) return("")
  if (length(hits) > 1 && all(hits %in% c("I", "L"))) hits <- "Xle"
  y1 <- hits[1]
  if (y1 == "Xle") {
    y3 <- "Xle"
  } else {
    y3 <- names(MSF_AA3)[MSF_AA3 == y1][1]
    if (is.na(y3)) return("")
  }
  x3 <- names(MSF_AA3)[MSF_AA3 == res][1]
  if (is.na(x3)) return("")
  paste0(x3, "->", y3)
}

msf_mod_string <- function(minfo, seq) {
  # token 形式（4.4.1 实测）："20M(15.9949)"（残基修饰）、"N-term(42.0106)"
  # （肽/蛋白 N 端修饰），多 token 以 ", " 分隔 -> "20,Oxidation[M];..."
  # token forms (4.4.1): "20M(15.9949)" (residue) and "N-term(42.0106)"
  # (N-terminal), ", "-separated -> "20,Oxidation[M];..."
  if (is.na(minfo) || !nzchar(minfo)) return("")
  toks <- trimws(strsplit(minfo, ",", fixed = TRUE)[[1]])
  toks <- toks[nzchar(toks)]
  out <- character(0)
  for (tk in toks) {
    m <- regmatches(tk, regexec("^([0-9]+)([A-Z]+)\\((-?[0-9.]+)\\)$", tk))[[1]]
    is_nterm <- FALSE
    if (length(m) == 0) {
      m <- regmatches(tk, regexec("^N-term\\((-?[0-9.]+)\\)$", tk))[[1]]
      if (length(m) == 0) next
      is_nterm <- TRUE
    }
    shift <- if (is_nterm) m[2] else m[4]
    nm_shifts <- as.numeric(names(MSF_MOD_NAME))
    hit <- which(abs(as.numeric(shift) - nm_shifts) < 1e-4)
    if (length(hit) == 0) {
      # 已知修饰表外位移：先按残基替换 X->Y 解释（g 槽 delta 口径），
      # 无解才丢弃——替换命中时按目标位置残基出 "pos,X->Y" 记号。
      # shifts outside the known-mod table: try the residue-substitution
      # interpretation first (the g-slot delta convention); only drop when
      # nothing matches - a hit emits "pos,X->Y" at the token position.
      if (!is_nterm) {
        res_tok <- m[3]
        sub_tok <- msf_delta_sub(res_tok, as.numeric(shift))
        if (nzchar(sub_tok)) {
          pos <- as.integer(m[2])
          out <- c(out, paste0(pos, ",", sub_tok, ";"))
          next
        }
      }
      cat("[MSFragger 适配] 未知修饰位移（忽略该 token）：", shift, "\n")
      next
    }
    nm <- MSF_MOD_NAME[[hit[1]]]
    if (is_nterm || (as.integer(m[2]) == 0L && names(MSF_MOD_NAME)[hit[1]] %in%
                     names(MSF_TERM_SITE))) {
      site <- MSF_TERM_SITE[[names(MSF_MOD_NAME)[hit[1]]]]
      out <- c(out, paste0("0,", nm, "[", site, "];"))
    } else {
      pos <- as.integer(m[2])
      res <- if (!is_nterm) m[3] else ""
      site <- if (nzchar(res) && grepl("[A-Z]", res)) res else substr(seq, pos, pos)
      out <- c(out, paste0(pos, ",", nm, "[", site, "];"))
    }
  }
  paste(out, collapse = "")
}


# 文件夹解析输入填文件夹时按目标模式递归（<=4 层）发现
# 结果文件、多个候选取最新、找不到报错；文件路径原样通过。
# folder resolution a directory input is resolved by the
# target pattern (recursive depth <= 4, newest wins, clear error otherwise);
# a file path passes through unchanged.

convert_msfragger_tsv <- function(tsv_file, map_file, out_file) {
  # MSFragger TSV 的表头比数据行少一列（4.4.1 实测：header 29 名 / data 30 列），
  # 不能按列名读——跳过表头按实测硬列位取数：
  #   V1=scannum V4=charge V7=hit_rank V8=peptide V11=proteins V18=modification_info
  # the MSFragger TSV header carries one name fewer than the data rows (4.4.1:
  # 29 header names / 30 data columns), so read headerless and index the
  # empirically verified positions
  tsv_file <- resolve_adapter_input(tsv_file, "(_renum\\.tsv$|\\.tsv$|\\.pepXML$|\\.pep\\.xml$)")
  d <- read.delim(tsv_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0), quote = "",
                  colClasses = "character", header = FALSE, skip = 1)
  names(d)[c(1, 4, 7, 8, 11, 18)] <- c("scannum", "charge", "hit_rank",
                                       "peptide", "proteins", "modification_info")
  mp <- read.delim(map_file, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = character(0), quote = "",
                   colClasses = "character")
  stopifnot(ncol(d) >= 18,
            all(c("new_scan", "original_title") %in% colnames(mp)))
  # 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（文件序＝引擎 rank 序）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[is.na(d$hit_rank) | d$hit_rank == "1",, drop = FALSE]   # 每谱 top-1
    d <- d[!duplicated(d$scannum),, drop = FALSE]
  }
  fn <- mp$original_title[match(d$scannum, mp$new_scan)]
  keep <- !is.na(fn) & nzchar(d$peptide)
  d <- d[keep,, drop = FALSE]; fn <- fn[keep]
  mods <- mapply(msf_mod_string, d$modification_info, d$peptide,
                 USE.NAMES = FALSE)
  prot <- vapply(strsplit(d$proteins, ";", fixed = TRUE), function(ps) {
    acc <- sub("\\s.*$", "", trimws(ps))
    acc <- acc[nzchar(acc)]
    if (length(acc) == 0) "" else paste0(paste(acc, collapse = "/"), "/")
  }, character(1))
  out <- data.frame(File_Name = fn, Sequence = d$peptide,
                    Modification = mods, Proteins = prot,
                    Charge = d$charge,
                    stringsAsFactors = FALSE, check.names = FALSE)
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("MSFragger TSV -> %s：%d 行\n", out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 3) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_msfragger_tsv(a[1], a[2], a[3])
}
