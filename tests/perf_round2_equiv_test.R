# perf_round2_equiv_test.R — 性能第二轮对拍（v1.0）
# 覆盖：matched_peak_set 预排序参数等价、理论离子缓存 CoW 不污染、
# filter_mods_compatible 展平与旧嵌套循环逐字节等价。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))
setwd(cd)

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

# 1) 预排序参数与自排序逐字节等价（同一谱多候选场景）
set.seed(23)
for (rep_i in 1:8) {
  mz <- runif(300, 100, 1500)
  ord <- order(mz, na.last = NA); mz_s <- mz[ord]
  r_base <- matched_peak_set("PEPTIDEK", "", list(mz), aa, aa_tab, 20)
  r_pre  <- matched_peak_set("PEPTIDEK", "", list(mz), aa, aa_tab, 20,
                             mz_sorted = mz_s, mz_ord = ord)
  if (!identical(as.integer(r_base), as.integer(r_pre))) break
}
ck("预排序参数与自排序逐字节等价（8 组随机谱）", rep_i == 8)

# 2) 理论离子缓存：同一 (seq|mod|label) 返回同一对象，写时复制不污染缓存
rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)
v1 <- cached_by_ions_position("PEPTIDEKM", "7,Oxidation[M];", "Unlabel", aa, aa_tab)
v2 <- cached_by_ions_position("PEPTIDEKM", "7,Oxidation[M];", "Unlabel", aa, aa_tab)
ck("缓存命中返回同一结果（identical）", identical(v1, v2))
v2$ions_table$b_1_charge[1] <- 0    # 调用方修改
v3 <- cached_by_ions_position("PEPTIDEKM", "7,Oxidation[M];", "Unlabel", aa, aa_tab)
ck("调用方修改不污染缓存（CoW）", identical(v1, v3) &&
   v1$ions_table$b_1_charge[1] != 0)
ck("缓存容量上限后仍正确（清空重填）", {
  rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)
  v4 <- cached_by_ions_position("PEPTIDEKM", "7,Oxidation[M];", "Unlabel", aa, aa_tab)
  identical(v4, v1)
})

# 3) filter_mods_compatible 展平 vs 旧嵌套循环（合成文件逐字节对拍）
old_filter <- function(path, compat_file) {
  compat <- readLines(compat_file, warn = FALSE)
  compat <- compat[nzchar(trimws(compat))]
  dat <- readLines(path, warn = FALSE)
  header <- dat[1]; body <- dat[-1]
  row_ok <- vapply(body, function(line) {
    for (cell in strsplit(line, "\t", fixed = TRUE)[[1]]) {
      for (piece in strsplit(cell, ";", fixed = TRUE)[[1]]) {
        if (!grepl("[", piece, fixed = TRUE)) next
        nm <- sub("^[^,]*?,", "", piece)
        if (nzchar(nm) && !(nm %in% compat)) return(FALSE)
      }
    }
    TRUE
  }, logical(1))
  writeLines(c(header, body[row_ok]), path)
  c(rows = length(body), kept = sum(row_ok), dropped = sum(!row_ok))
}
td <- tempfile("fmc_"); dir.create(td)
compat_f <- file.path(td, "compat.txt")
writeLines(c("Oxidation[M]", "Acetyl[AnyN-term]"), compat_f)
mk_lines <- function(n) {
  c("H1\tH2\tH3", vapply(seq_len(n), function(i) {
    paste(sample(c("3,Oxidation[M];", "0,Acetyl[AnyN-term];", "5,Phospho[S];",
                   "9,Unknown[X];", "", "1,Carbamidomethyl[C];"), 3, replace = TRUE),
          collapse = "\t")
  }, character(1)))
}
set.seed(9)
lines <- mk_lines(600)
res <- list()
for (ver in c("new", "old")) {
  f <- file.path(td, paste0(ver, ".spectra"))
  writeLines(lines, f)
  r <- if (ver == "new") filter_mods_compatible(f, compat_f)
       else old_filter(f, compat_f)
  res[[ver]] <- list(r = r, out = readLines(f, warn = FALSE))
}
ck("filter_mods_compatible 展平与旧循环计数一致",
   identical(res$new$r, res$old$r))
ck("filter_mods_compatible 输出文件逐行一致",
   identical(res$new$out, res$old$out))
unlink(td, recursive = TRUE)

# 4) RT 窗口（v1.0，P1-2 修复）：与旧 which() 逐位等价（集合同 + 枚举
#    顺序同），NA RT / NA tol / NA light_rt 边界，预排序参数路径等价
set.seed(31)
rt_ok <- TRUE
rt_last <- c(runif(60, 0, 100), NA_real_)
for (rep_j in 1:10) {
  rt <- c(runif(60, 0, 100), NA_real_)
  lr <- runif(1, 20, 80); tol <- runif(1, 2, 15)
  old_w <- which(rt >= lr - tol & rt <= lr + tol)
  new_w <- rt_window_indices(rt, lr, tol)
  if (!identical(as.integer(old_w), as.integer(new_w))) { rt_ok <- FALSE; break }
  rt_last <- rt
}
ck("RT 窗口与旧 which() 逐位等价（10 组随机含 NA，P1-2）", rt_ok)
ck("RT 窗口 NA 输入返回空（NA light_rt/NA tol）",
   length(rt_window_indices(rt_last, NA_real_, 20)) == 0L &&
   length(rt_window_indices(rt_last, 50, NA_real_)) == 0L)
ort <- order(rt_last, na.last = NA)
ck("RT 窗口预排序传入与自排序等价",
   identical(rt_window_indices(rt_last, 50, 10),
             rt_window_indices(rt_last, 50, 10, ort)))

cat(sprintf("=== perf_round2_equiv_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
