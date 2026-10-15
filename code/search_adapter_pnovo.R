# =============================================================================
# search_adapter_pnovo.R - pNovo3 results.res -> isoPSV 交换格式（o 槽候选）
# search_adapter_pnovo.R - pNovo3 results.res -> isoPSV interchange (o-slot)
#
# pNovo+ 3.1.3（pFind 家族）：pNovo3_Search.exe 从安装目录启动（无参数，
# 读 temp/pNovo.tsk + temp/param/pNovo.param，结果写 temp/result/）。
# results.res 每行一个候选（tab 分隔）：标题 <tab> 序列 <tab> 分值...
# 每谱多行按名次排列——取每标题首行（top-1）。
# 序列字母口径：大写＝标准残基（FixMod 不改变字母）；小写 a-z＝param 里
# VarMod 定义的可变修饰残基（"a=147.035405 Oxidation[M]"→a 即 M）；
# 序列前导数字 0-9＝N 端修饰标记。本适配器按 param 的 VarMod 行把可还原的
# 小写字母还原成残基（[X] 后缀），无法还原的字母丢弃（诚实降级：回帖口径
# 与 e 槽 denovo 一致，denovo_db_lookback 再作 I→L、N→GG 替换回帖）。
# pNovo+ 3.1.3 (pFind family): pNovo3_Search.exe runs from the install dir
# with no args (reads temp/pNovo.tsk + temp/param/pNovo.param, writes
# temp/result/). results.res holds one ranked candidate per line
# (tab-separated): title <tab> sequence <tab> scores... - take the first line
# of each title (top-1). Sequence alphabet: uppercase = standard residues
# (FixMods keep their letter); lowercase a-z = VarMod residues defined in the
# param ("a=147.035405 Oxidation[M]" -> a is M); leading digits 0-9 mark
# N-terminal mods. This adapter restores lowercase letters to residues via
# the param VarMod lines ([X] suffix); letters it cannot resolve are dropped
# (honest degradation: the look-back matches the e-slot denovo semantics,
# denovo_db_lookback then applies I->L / N->GG and matches the database).
#
# 用法：Rscript search_adapter_pnovo.R <results.res 或结果文件夹> <pNovo.param> <输出.spectra>
# usage: Rscript search_adapter_pnovo.R <results.res or folder> <pNovo.param> <out.spectra>
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

# 文件夹解析（与其余适配器同口径）：目录时按模式递归（<=4 层）取最新文件

# 从 pNovo.param 解析小写字母 → 残基映射（VarMod 行 "a=147.035405 Oxidation[M]"；
# 名字以 [X] 结尾时 X 即残基；FixMod 行 "C=160.030654 Carbamidomethyl[C]" 的
# 大写字母不需转换）
# letter -> residue map from pNovo.param (VarMod lines "a=147.035405
# Oxidation[M]"; the [X] suffix names the residue; FixMod lines such as
# "C=160.030654 Carbamidomethyl[C]" need no conversion)
pnovo_var_residue <- function(param_lines) {
  map <- character(0)
  for (l in param_lines) {
    l <- trimws(l)
    if (!grepl("^[a-z][0-9]*=", l)) next
    letter <- sub("=.*$", "", l)
    rest <- sub("^[^=]*=", "", l)
    m <- regmatches(rest, regexec("\\[([A-Z])\\]\\s*$", rest))[[1]]
    if (length(m) == 2) map[letter] <- m[2]
  }
  map
}

# 序列还原：前导数字（N 端修饰标记）去除，小写字母按映射还原，无法还原丢弃
# restore the sequence: strip leading digits (N-terminal mod marks), restore
# lowercase letters via the map, drop letters the map cannot resolve
pnovo_restore_seq <- function(sq, vmap) {
  sq <- sub("^[0-9]+", "", sq)
  ch <- strsplit(sq, "", fixed = TRUE)[[1]]
  out <- character(0)
  for (x in ch) {
    if (x %in% LETTERS) {
      out <- c(out, x)
    } else if (x %in% names(vmap)) {
      out <- c(out, vmap[[x]])
    }   # 其余字符（无法还原的小写字母）丢弃
  }
  paste(out, collapse = "")
}

convert_pnovo_res <- function(res_file, param_file, out_file) {
  res_file <- resolve_adapter_input(res_file, "\\.res$|\\.txt$")
  pl <- readLines(param_file, warn = FALSE, encoding = "UTF-8")
  vmap <- pnovo_var_residue(pl)
  cat("pNovo VarMod letter map:", paste(names(vmap), vmap, sep = "->", collapse = ", "), "\n")
  d <- read.delim(res_file, sep = "\t", stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = character(0),
                  quote = "", colClasses = "character", header = FALSE)
  keep <- which(!is.na(d[[1]]) & nzchar(d[[1]]) & !is.na(d[[2]]) & nzchar(d[[2]]))
  d <- d[keep,, drop = FALSE]
  # 全候选模式（ISOPSV_CAND_ALL=1）：保留全部候选（每标题按名次排列，文件序）
  if (!identical(Sys.getenv("ISOPSV_CAND_ALL"), "1")) {
    d <- d[!duplicated(d[[1]]),, drop = FALSE]      # 每谱首行＝top-1
  }
  fn <- d[[1]]
  seq <- vapply(d[[2]], pnovo_restore_seq, character(1), vmap = vmap,
                USE.NAMES = FALSE)
  charge <- vapply(strsplit(fn, ".", fixed = TRUE), function(f) {
    if (length(f) >= 5) f[length(f) - 2] else ""
  }, character(1), USE.NAMES = FALSE)
  out <- data.frame(File_Name = fn, Sequence = seq, Modification = "",
                    Proteins = "", Charge = charge,
                    stringsAsFactors = FALSE, check.names = FALSE)
  out <- out[!is.na(out$Sequence) & nzchar(out$Sequence),, drop = FALSE]
  write.table(out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("pNovo results.res -> %s：%d 行（top-1，还原后非空）\n",
              out_file, nrow(out)))
  invisible(out)
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 3) {
  a <- commandArgs(trailingOnly = TRUE)
  convert_pnovo_res(a[1], a[2], a[3])
}
