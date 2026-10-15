# =============================================================================
# search_adapter_xtandem.R - X!Tandem XML -> isoPSV 交换格式（k 槽候选）
# search_adapter_xtandem.R - X!Tandem XML -> isoPSV interchange (k-slot)
#
# X!Tandem（thegpm.org/TANDEM）输出 XML（tandem.xml）。每个 <group> 块含
# label（谱图标题，取自输入峰列表，直通 File_Name）、mh/z/expect 与若干
# <protein>/<peptide>；<peptide> 内 <domain> 的 seq 属性为肽段序列，可变修饰
# 以 "[delta]" 内联在残基后（如 M[15.9949]）；N 端修饰写为前导 "n[delta]"、
# C 端修饰写为尾随 "c[delta]"。每谱取第一个 <peptide> 的第一个 <domain>
# （top-1），按 label 去重（expect 最小者优先）。
# X!Tandem (thegpm.org/TANDEM) outputs XML (tandem.xml). Each <group> block
# carries label (the spectrum title taken from the input peak list, passed
# through to File_Name), mh/z/expect and <protein>/<peptide> entries; the seq
# attribute of the <domain> inside <peptide> is the sequence, variable mods
# inline as "[delta]" after the residue (e.g. M[15.9949]); N-terminal mods are
# a leading "n[delta]" token and C-terminal mods a trailing "c[delta]". The
# first <domain> of the first <peptide> is kept (top-1); dedup by label keeps
# the lowest expect.
#
# 用法：Rscript search_adapter_xtandem.R <tandem.xml 或文件夹> <输出.spectra>
# usage: Rscript search_adapter_xtandem.R <tandem.xml or folder> <out.spectra>
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


# X!Tandem 内联修饰位移 -> isoPSV 修饰名（与 Sage/PEAKS 适配器同一张表）
XT_MOD_NAME <- c(
  "15.9949" = "Oxidation", "15.994915" = "Oxidation",
  "57.021464" = "Carbamidomethyl", "57.02146" = "Carbamidomethyl",
  "0.984016" = "Deamidated", "0.984" = "Deamidated",
  "42.010565" = "Acetyl", "42.0106" = "Acetyl",
  "79.966331" = "Phospho", "79.9663" = "Phospho",
  "-17.026549" = "Gln->pyro-Glu", "-17.0265" = "Gln->pyro-Glu",
  "-18.010565" = "Glu->pyro-Glu", "-18.0106" = "Glu->pyro-Glu",
  "14.01565" = "Methyl", "28.0313" = "Dimethyl",
  "71.037114" = "Propionamide"
)

xtandem_mod_string <- function(seq0) {
  out <- character(0)
  pick <- function(dlt) {                        # 容差内取第一个命中（表内近义键并存）
    hit <- which(abs(dlt - as.numeric(names(XT_MOD_NAME))) < 0.05)
    if (length(hit) == 0) 0L else hit[1]
  }
  # 前导 n[delta]＝N 端修饰、尾随 c[delta]＝C 端修饰
  # （regmatches 取捕获字符串；regexec 裸返回的是位置向量，勿直接 as.numeric）
  m_n <- regmatches(seq0, regexec("^n\\[([+-]?[0-9.]+)\\]", seq0))[[1]]
  if (length(m_n) == 2) {
    hit <- pick(as.numeric(m_n[2]))
    if (hit > 0L) {
      out <- c(out, paste0("0,", XT_MOD_NAME[[hit]], "[AnyN-term];"))
      seq0 <- sub("^n\\[[+-]?[0-9.]+\\]", "", seq0)
    }
  }
  m_c <- regmatches(seq0, regexec("c\\[([+-]?[0-9.]+)\\]$", seq0))[[1]]
  if (length(m_c) == 2) {
    hit <- pick(as.numeric(m_c[2]))
    if (hit > 0L) {
      seq0 <- sub("c\\[[+-]?[0-9.]+\\]$", "", seq0)
      out <- c(out, paste0(nchar(seq0) + 1L, ",", XT_MOD_NAME[[hit]],
                           "[AnyC-term];"))
    }
  }
  clean <- gsub("\\[[^]]*\\]", "", seq0)
  m <- gregexpr("[A-Z]\\[[+-]?[0-9.]+\\]", seq0)[[1]]
  if (m[1] != -1) {
    for (i in seq_along(m)) {
      tok <- substr(seq0, m[i], m[i] + attr(m, "match.length")[i] - 1L)
      rm <- regmatches(tok, regexec("^([A-Z])\\[([+-]?)([0-9.]+)\\]$", tok))[[1]]
      if (length(rm) == 0) next
      dlt <- as.numeric(paste0(rm[3], rm[4]))
      hit <- pick(dlt)
      if (hit == 0L) {
        cat("[X!Tandem 适配] 未知修饰位移（忽略）：", dlt, "\n")
        next
      }
      res <- rm[2]
      pos <- nchar(gsub("\\[[^]]*\\]", "", substr(seq0, 1, m[i])))
      out <- c(out, paste0(pos, ",", XT_MOD_NAME[[hit]], "[", res, "];"))
    }
  }
  list(seq = clean, mod = paste(sort(unique(out)), collapse = ""))
}

convert_xtandem_xml <- function(xml_file, out_file, map_file = NULL) {
  xml_file <- resolve_adapter_input(xml_file, "\\.xml$")
  L <- readLines(xml_file, warn = FALSE)
  L <- enc2utf8(L)
  txt <- paste(L, collapse = "\n")
  # P0-5TRE 正则的 . 不匹配换行，真实 tandem.xml 是
  # 多行格式（测试夹具单行掩盖了问题）——perl=TRUE + (?s) 使 . 跨行匹配；同时
  # 用 [[:space:]] 前瞻替代 \\b（\\b 在 C locale 下行为不可靠，实测 R 4.5.1
  # LC_CTYPE=C 时恒不匹配）。
  groups <- regmatches(txt,
                       gregexpr("(?s)<group(?=[[:space:]])[^>]*>.*?</group>",
                                txt, perl = TRUE))[[1]]
  if (length(groups) == 0) stop("no <group> blocks found in ", xml_file)
  n <- length(groups)
  fn <- character(n); seqs <- character(n); mods <- character(n)
  expect <- numeric(n); charge <- character(n)
  # P1新版 X!Tandem 的 <group label> 是蛋白描述
  # 而非谱图标题（按蛋白名去重会把 46,610 组塌缩成 1,488 行）——提供可选 map
  # 参数（new_scan<TAB>original_title）时按 <group id>（= 输入谱 1-based 序）
  # 解析标题；无 map 或 id 缺失/越界时保留 label 回退并计数上报
  mp <- NULL
  if (!is.null(map_file) && nzchar(map_file) && file.exists(map_file)) {
    mp <- read.delim(map_file, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, na.strings = character(0),
                     quote = "", colClasses = "character")
    if (!all(c("new_scan", "original_title") %in% colnames(mp))) mp <- NULL
  }
  n_map_hit <- 0L; n_map_miss <- 0L
  # 全候选模式（ISOPSV_CAND_ALL=1）：每组取前 5 个 <domain>（引擎自身排序，
  # expect 升序）；常规模式取第一个 <domain>（top-1，现状）
  all_mode <- identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")
  max_dom <- if (all_mode) 5L else 1L
  out_rows <- list()
  for (i in seq_len(n)) {
    g <- groups[i]
    head_txt <- sub(">.*$", ">", regmatches(g, regexpr("<group(?=[[:space:]])[^>]*>", g, perl = TRUE)))
    lab <- regmatches(head_txt, regexec('label="([^"]*)"', head_txt))[[1]]
    fn[i] <- if (length(lab) == 2) lab[2] else paste0("group_", i)
    # map 优先：id = 输入峰列表 1-based 序
    if (!is.null(mp)) {
      gid <- regmatches(head_txt, regexec('id="([0-9]+)"', head_txt))[[1]]
      if (length(gid) == 2) {
        mi <- match(gid[2], mp$new_scan)
        if (!is.na(mi)) {
          fn[i] <- mp$original_title[mi]
          n_map_hit <- n_map_hit + 1L
        } else {
          n_map_miss <- n_map_miss + 1L
        }
      } else {
        n_map_miss <- n_map_miss + 1L
      }
    }
    zz <- regmatches(head_txt, regexec('z="([0-9]+)"', head_txt))[[1]]
    charge[i] <- if (length(zz) == 2) zz[2] else ""
    ee <- regmatches(head_txt, regexec('expect="([0-9.eE+-]+)"', head_txt))[[1]]
    expect[i] <- if (length(ee) == 2) suppressWarnings(as.numeric(ee[2])) else 1
    # <peptide> 候选域：常规取第一个；全候选模式取前 max_dom 个
    doms <- regmatches(g, gregexpr("<domain(?=[[:space:]])[^>]*>", g, perl = TRUE))[[1]]
    if (length(doms) == 0) next
    doms <- doms[seq_len(min(max_dom, length(doms)))]
    for (dom in doms) {
      sq <- regmatches(dom, regexec('seq="([^"]*)"', dom))[[1]]
      if (length(sq) != 2) next
      conv <- xtandem_mod_string(sq[2])
      out_rows[[length(out_rows) + 1L]] <-
        c(fn[i], charge[i], as.character(expect[i]), conv$seq, conv$mod)
    }
  }
  if (!is.null(mp)) {
    cat(sprintf("[X!Tandem 适配] map 命中 %d / 未命中回退 label %d\n",
                n_map_hit, n_map_miss))
  }
  if (length(out_rows) == 0L) {
    d <- data.frame(File_Name = character(0), Sequence = character(0),
                    Modification = character(0), Proteins = character(0),
                    Charge = character(0), expect = numeric(0),
                    stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    m <- do.call(rbind, out_rows)
    d <- data.frame(File_Name = m[, 1], Sequence = m[, 4],
                    Modification = m[, 5], Proteins = "",
                    Charge = m[, 2], expect = as.numeric(m[, 3]),
                    stringsAsFactors = FALSE, check.names = FALSE)
  }
  d <- d[order(d$expect, na.last = TRUE),, drop = FALSE]
  if (!all_mode) d <- d[!duplicated(d$File_Name),, drop = FALSE]
  d$expect <- NULL
  d$Proteins[d$Proteins == ""] <- ""
  write.table(d, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("X!Tandem XML -> %s：%d 行（%s）\n", out_file, nrow(d),
              if (all_mode) "全候选（前 5/组）" else "top-1 每谱去重"))
  invisible(d)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_xtandem_xml(a[1], a[2], if (length(a) >= 3) a[3] else NULL)
}
