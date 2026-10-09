# matched_peak_set_test.R — findInterval 化后的新旧等价性回归（2026-09-26）
# 用同一肽段与合成谱图，断言新实现与旧线性扫描实现逐字节同义（含多修饰变体、
# 无匹配窗口、峰重复等边界）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R"))) file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))
setwd(cd)   # AA_mass_dataframe.txt 在 code/ 下
on.exit(setwd(sd), add = TRUE)
aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 旧实现（线性扫描）本地副本，用于逐字节对比
matched_peak_set_old <- function(seq, mod, sp, aa, aa_tab, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) ||
        length(bip) == 1 || is.null(bip$ions_table)) next
    t <- bip$ions_table
    mz <- as.numeric(sp[[1]])
    mz_all <- c(t$b_1_charge, t$b_2_charge, t$y_1_charge, t$y_2_charge)
    lo <- mz_all * (1 - tol_ppm * 1e-6)
    hi <- mz_all * (1 + tol_ppm * 1e-6)
    hits <- lapply(seq_along(mz_all), function(k)
      which(mz >= lo[k] & mz <= hi[k]))
    return(sort(unique(unlist(hits))))
  }
  integer(0)
}

set.seed(11)
ck_equiv <- function(seq, mod, sp_mz, tol_ppm = 20) {
  new_res <- matched_peak_set(seq, mod, list(sp_mz), aa, aa_tab, tol_ppm)
  old_res <- matched_peak_set_old(seq, mod, list(sp_mz), aa, aa_tab, tol_ppm)
  identical(as.integer(new_res), as.integer(old_res))
}
ck("无修饰肽随机谱等价", ck_equiv("PEPTIDEK", "", runif(400, 100, 1500)))
ck("修饰肽随机谱等价", ck_equiv("PEPTIDEKM", "7,Oxidation[M];", runif(400, 100, 1500)))
ck("双修饰变体随机谱等价", ck_equiv("PEPTIDEKM", "7,Oxidation[M];7,Phospho[M];", runif(400, 100, 1500)))
ck("无匹配窗口（谱远离理论离子）等价", ck_equiv("PEPTIDEK", "", rep(9999, 50)))
ck("窄容差等价", ck_equiv("PEPTIDEK", "", runif(300, 100, 1500), tol_ppm = 5))
ck("宽容差等价", ck_equiv("PEPTIDEK", "", runif(300, 100, 1500), tol_ppm = 50))
ck("空谱等价", ck_equiv("PEPTIDEK", "", numeric(0)))
# 含 NA 的观测谱（新实现 na.last=NA 剔除，旧实现 which(NA 比较) 会报错——新实现更稳）
new_na <- tryCatch(matched_peak_set("PEPTIDEK", "", list(c(100, NA, 200)), aa, aa_tab, 20),
                   error = function(e) "error")
ck("含 NA 观测谱不报错", is.integer(new_na))

# P0-D（2026-09-26）：理论离子表含 NA 的情形——用桩 by_ions_position 注入
# NA m/z 离子，新旧实现都不报错且结果一致（旧实现 which 静默丢弃 NA 窗口）。
# 注意：桩必须返回长度 >1 的 list，否则 matched_peak_set 的
# length(ion_positions)==1 守卫会直接整函数退出（返回 integer(0)），
# 测试形同虚设——extra 元素就是为此存在。
by_ions_position_stub <- function(seq, mod, label, aa, aa_tab) {
  list(ions_table = data.frame(
    b_1_charge = c(100, 200, NA, 400),
    b_2_charge = c(50, 100, 150, NA),
    y_1_charge = c(80, 160, 240, 320),
    y_2_charge = c(40, 80, 120, 160),
    stringsAsFactors = FALSE),
    extra = NULL)
}
by_ions_position_saved <- by_ions_position
assign("by_ions_position", by_ions_position_stub, envir = globalenv())
on.exit(assign("by_ions_position", by_ions_position_saved, envir = globalenv()), add = TRUE)
# v1.0：matched_peak_set 走进程级缓存——装桩后清一次缓存，否则前面随机
# 等价测试存入的真实 by_ions_position 结果会让桩失效
if (exists("theo_cache_env", inherits = FALSE)) {
  rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)
}
new_na_theo <- tryCatch(matched_peak_set("PEPTIDEK", "", list(c(100, 200, 50, 80, 160, 320, 450)),
                                         aa, aa_tab, 20),
                        error = function(e) "error")
old_na_theo <- tryCatch(matched_peak_set_old("PEPTIDEK", "", list(c(100, 200, 50, 80, 160, 320, 450)),
                                             aa, aa_tab, 20),
                        error = function(e) "error")
ck("理论离子含 NA 不报错（P0-D）", is.integer(new_na_theo))
ck("理论离子含 NA 新旧等价（P0-D）", identical(as.integer(new_na_theo), as.integer(old_na_theo)))

# GLM v1.0 P1（2026-09-26）：观测峰恰等于窗口下界（mz == theoretical*(1-ppm)）
# ——旧实现 which(mz >= lower) 含等值峰；v1.0 的 findInterval+1 会把它漏掉。
# 构造 4 个峰恰落在 y_1 离子（80/160/240/320）的下界上，新旧实现都必须命中 4 个。
boundary_mz <- c(80, 160, 240, 320) * (1 - 20 * 1e-6)
new_boundary <- tryCatch(matched_peak_set("PEPTIDEK", "", list(boundary_mz), aa, aa_tab, 20),
                         error = function(e) "error")
old_boundary <- tryCatch(matched_peak_set_old("PEPTIDEK", "", list(boundary_mz), aa, aa_tab, 20),
                         error = function(e) "error")
ck("恰在下界的峰被匹配（GLM v1.0 P1）",
   is.integer(new_boundary) && length(new_boundary) == 4L)
ck("恰在下界的峰新旧等价（GLM v1.0 P1）",
   identical(as.integer(new_boundary), as.integer(old_boundary)))

cat(sprintf("=== matched_peak_set_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
