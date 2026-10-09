# theo_cache_isolation_test.R — 理论离子缓存键隔离回归（deepseek-flash 20260926 审查提出）
# 背景：cached_by_ions_position 与 cached_all_ions_position 共用 theo_cache_env 且
# 键格式相同（"<pep>\r<mod>\r<label>"），但返回结构不同——by=宽表(b_1_charge…)、
# all=长表(Class/Charge/m_z)+internal_aa_table。stageA_evaluate/scan_chunk 先调
# spectrum_params（all）再调 matched_peak_set（by），后者取到 all 对象 →
# theoretical_mz 为空 → 直接 return(integer(0))。
# 期望：修复后本套件全 PASS（缓存键须含函数身份，或两函数各用独立环境）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "paired_params.R")))
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
reset_cache <- function() rm(list = ls(theo_cache_env, all.names = TRUE), envir = theo_cache_env)

# 1) 结构隔离：all 先入缓存后，by 必须仍返回宽表（含 b_1_charge）
reset_cache()
a1 <- cached_all_ions_position("PEPTIDEK", "", "Unlabel", aa, aa_tab)
b1 <- cached_by_ions_position("PEPTIDEK", "", "Unlabel", aa, aa_tab)
ck("all 入缓存后 by 仍返回宽表（含 b_1_charge）", !is.null(b1$ions_table$b_1_charge))
ck("by 不与 all 共用键（返回结构不同）", !identical(names(b1), names(a1)))

# 2) 反向：by 先入缓存后，all 必须仍返回长表（含 internal_aa_table）
reset_cache()
b2 <- cached_by_ions_position("PEPTIDEK", "", "Unlabel", aa, aa_tab)
a2 <- cached_all_ions_position("PEPTIDEK", "", "Unlabel", aa, aa_tab)
ck("by 入缓存后 all 仍含 internal_aa_table", !is.null(a2$internal_aa_table))
ck("all 的 ions_table 含 Class/Charge/m_z",
   all(c("Class", "Charge", "m_z") %in% colnames(a2$ions_table)))

# 3) 端到端：复刻 scan_chunk 顺序（spectrum_params 先，matched_peak_set 后）
reset_cache()
all_ions <- all_ions_position("PEPTIDEK", "", "Unlabel", aa, aa_tab)
spec <- data.frame(m_z = c(all_ions$ions_table$m_z[1:6], 1e4), Intensity = rep(1000, 7))
spec_mz <- as.numeric(spec[[1]]); spec_ord <- order(spec_mz, na.last = NA)
sp_par <- spectrum_params("PEPTIDEK", "", "unlabeled", aa, aa_tab, spec, 20, "ppm")
pk <- matched_peak_set("PEPTIDEK", "", spec, aa, aa_tab, tol_ppm = 20,
                       mz_sorted = spec_mz[spec_ord], mz_ord = spec_ord)
ck("spectrum_params 后 matched_peak_set 仍能命中理论离子（期望 6）",
   length(pk) == 6L)

# 4) 反向顺序不得报错（spectrum_params 不得取到 by 对象）
reset_cache()
pk2 <- matched_peak_set("PEPTIDEK", "", spec, aa, aa_tab, tol_ppm = 20,
                        mz_sorted = spec_mz[spec_ord], mz_ord = spec_ord)
sp2 <- tryCatch(spectrum_params("PEPTIDEK", "", "unlabeled", aa, aa_tab, spec, 20, "ppm"),
                error = function(e) NULL)
ck("matched_peak_set 后 spectrum_params 不报错且含 internal_aa_table",
   !is.null(sp2) && !is.null(sp2$internal_aa_table))

cat(sprintf("=== theo_cache_isolation_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
