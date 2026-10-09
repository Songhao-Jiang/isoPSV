# enzyme_test.R — 酶切类型选项单元测试（v1.0，2026-09-26）
# 覆盖：酶名归一（含中文别名/未知回退）、漏切计数（各酶 + 与历史 [KR] 口径逐字节
# 一致）、N 端特异性判据（C 端酶/N 端酶/无酶切）、error_reason_compute 的酶依赖
# 分类（经 ISO_ENZYME 环境变量切换）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R"))) file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# ---- 1. 酶名归一 ----
ck("Trypsin → trypsin", identical(enzyme_resolve("Trypsin"), "trypsin"))
ck("Trypsin/P → trypsin（统一 P 切不开）", identical(enzyme_resolve("Trypsin/P"), "trypsin"))
ck("Lys-C → lysc（连字符归一）", identical(enzyme_resolve("Lys-C"), "lysc"))
ck("Lys_N → lysn（下划线归一）", identical(enzyme_resolve("Lys_N"), "lysn"))
ck("GluC → gluc", identical(enzyme_resolve("GluC"), "gluc"))
ck("糜蛋白酶 → chymotrypsin", identical(enzyme_resolve("糜蛋白酶"), "chymotrypsin"))
ck("LysargiNase → lysarginase", identical(enzyme_resolve("LysargiNase"), "lysarginase"))
# 中文别名（v1.0 扩充）
ck("胰蛋白酶 → trypsin", identical(enzyme_resolve("胰蛋白酶"), "trypsin"))
ck("胰凝乳蛋白酶 → chymotrypsin", identical(enzyme_resolve("胰凝乳蛋白酶"), "chymotrypsin"))
ck("无酶切 → none", identical(enzyme_resolve("无酶切"), "none"))
# P 位阻断（用户口径 2026-09-26：LysargiNase 同 Trypsin，K/R 后遇 P 不切；
# Trypsin 统一为 P 切不开，旧写法 Trypsin/P 归一为 trypsin）
ck("trypsin nocut = P", identical(enzyme_def("trypsin")$nocut, "P"))
ck("lysarginase nocut = P（同 Trypsin）", identical(enzyme_def("lysarginase")$nocut, "P"))
ck("Trypsin/P 写法归一为 trypsin（P 切不开）",
   identical(enzyme_resolve("Trypsin/P"), "trypsin") &&
     identical(enzyme_def("Trypsin/P")$nocut, "P"))
ck("gluc nocut = P", identical(enzyme_def("gluc")$nocut, "P"))
ck("chymotrypsin nocut = P", identical(enzyme_def("chymotrypsin")$nocut, "P"))
# nocut 不影响管线内部计数（历史口径），仅进引擎配置行
ck("nocut 不进漏切计数：trypsin KPKA → 2（历史口径，KP 不豁免）",
   enzyme_missed_cleavages("trypsin", "KPKA") == 2)
ck("lysarginase KPKA → 1（K@3 非首位；P 位不影响计数）",
   enzyme_missed_cleavages("lysarginase", "KPKA") == 1)
ck("空串 → trypsin", identical(enzyme_resolve(""), "trypsin"))
ck("NULL → trypsin", identical(enzyme_resolve(NULL), "trypsin"))
wk <- tryCatch({ enzyme_resolve("bogus") }, warning = function(w) conditionMessage(w))
ck("未知酶警告并回退 trypsin", identical(enzyme_resolve("bogus"), "trypsin") && length(wk) == 1)

# ---- 2. 漏切计数 ----
ck("trypsin AAKAAK → 1（K@3 非末位）", enzyme_missed_cleavages("trypsin", "AAKAAK") == 1)
ck("trypsin KAAAK → 1（K@1 非末位；与历史 [KR] 口径一致）", enzyme_missed_cleavages("trypsin", "KAAAK") == 1)
ck("trypsin KAKK → 2（K@1,K@3 非末位）", enzyme_missed_cleavages("trypsin", "KAKK") == 2)
ck("lysc AARK → 0（仅 K 为切点，K 在末位）", enzyme_missed_cleavages("lysc", "AARK") == 0)
ck("lysc KAAR → 1（K@1 非末位）", enzyme_missed_cleavages("lysc", "KAAR") == 1)
ck("lysn KAAK → 1（K@4 非首位）", enzyme_missed_cleavages("lysn", "KAAK") == 1)
ck("lysn KAKA → 1（K@1 是切点，K@3 非首位）", enzyme_missed_cleavages("lysn", "KAKA") == 1)
ck("gluc AADEAE → 2（D@3/E@4 非末位，E@6 末位）", enzyme_missed_cleavages("gluc", "AADEAE") == 2)
ck("aspn DAAD → 1（D@1 是切点，D@4 非首位）", enzyme_missed_cleavages("aspn", "DAAD") == 1)
ck("chymotrypsin AAFWY → 2（F/W 非末位，Y 末位）", enzyme_missed_cleavages("chymotrypsin", "AAFWY") == 2)
ck("none → 0", enzyme_missed_cleavages("none", "AARKK") == 0)
# 与历史 [KR] 口径逐字节一致（trypsin 默认）：随机序列 2000 条
set.seed(7)
aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
old_miss <- function(s) { k <- gregexpr("[KR]", s, fixed = FALSE)[[1]]; sum(k > 0 & k != nchar(s)) }
seqs <- replicate(2000, paste0(sample(aa, sample(5:25, 1), replace = TRUE), collapse = ""))
ck("trypsin 漏切与历史 [KR] 口径 2000 条随机序列一致",
   all(vapply(seqs, function(s) enzyme_missed_cleavages("trypsin", s) == old_miss(s), logical(1))))

# ---- 3. N 端特异性判据 ----
ck("trypsin prev=K 特异", enzyme_nterm_ok("trypsin", "K", "A"))
ck("trypsin prev=R 特异", enzyme_nterm_ok("trypsin", "R", "A"))
ck("trypsin prev=A 非特异", !enzyme_nterm_ok("trypsin", "A", "A"))
ck("trypsin prev=_（蛋白 N 端）特异", enzyme_nterm_ok("trypsin", "_", "M"))
ck("lysc prev=K 特异", enzyme_nterm_ok("lysc", "K", "A"))
ck("lysc prev=R 非特异", !enzyme_nterm_ok("lysc", "R", "A"))
ck("gluc prev=D 特异", enzyme_nterm_ok("gluc", "D", "A"))
ck("gluc prev=E 特异", enzyme_nterm_ok("gluc", "E", "A"))
ck("gluc prev=K 非特异", !enzyme_nterm_ok("gluc", "K", "A"))
ck("lysn first=K 特异（前残基任意）", enzyme_nterm_ok("lysn", "A", "K"))
ck("lysn first=A 非特异", !enzyme_nterm_ok("lysn", "A", "A"))
ck("lysn 蛋白 N 端（prev=_）特异（首残基任意）", enzyme_nterm_ok("lysn", "_", "M"))
ck("lysarginase first=R 特异", enzyme_nterm_ok("lysarginase", "A", "R"))
ck("lysarginase first=A 非特异", !enzyme_nterm_ok("lysarginase", "A", "A"))
ck("none 恒特异", enzyme_nterm_ok("none", "A", "A") && enzyme_nterm_ok("none", "_", "M"))

# ---- 4b. pFind/Comet 引擎行 helper（v1.0）----
ck("pfind_enzyme_line trypsin → Trypsin KR P C",
   identical(pfind_enzyme_line("trypsin"), "enzyme=Trypsin KR P C"))
ck("pfind_enzyme_line lysarginase → LysargiNase KR P N",
   identical(pfind_enzyme_line("lysarginase"), "enzyme=LysargiNase KR P N"))
ck("pfind_enzyme_line none → NoEnzyme U _ C",
   identical(pfind_enzyme_line("none"), "enzyme=NoEnzyme U _ C"))
# 真实 comet_con.params 的 COMET_ENZYME_INFO 表按名取号（E-2 硬编码消除）
cd2 <- if (file.exists(file.path(cd, "engine_templates", "comet_con.params")))
  file.path(cd, "engine_templates", "comet_con.params") else ""
if (nzchar(cd2)) {
  cpl <- readLines(cd2, warn = FALSE)
  ck("comet 查表 trypsin → 1", identical(comet_enzyme_number("trypsin", cpl), 1L))
  ck("comet 查表 lysc → 3", identical(comet_enzyme_number("lysc", cpl), 3L))
  ck("comet 查表 lysn → 4", identical(comet_enzyme_number("lysn", cpl), 4L))
  ck("comet 查表 aspn → 6", identical(comet_enzyme_number("aspn", cpl), 6L))
  ck("comet 查表 chymotrypsin → 10", identical(comet_enzyme_number("chymotrypsin", cpl), 10L))
  ck("comet 查表 none → 0", identical(comet_enzyme_number("none", cpl), 0L))
  ck("comet 查表 lysarginase 未命中 → NA", is.na(comet_enzyme_number("lysarginase", cpl)))
  ap <- comet_enzyme_append_line("lysarginase", cpl)
  ck("comet 追加 lysarginase 编号 = max+1 且含 KR/P",
     grepl("^12\\.\\s+LysargiNase\\s+0\\s+KR\\s+P$", ap))
  ap2 <- comet_enzyme_append_line("gluc", cpl)
  ck("comet 追加 gluc 编号 = max+1 且含 DE/P",
     grepl("^12\\.\\s+GluC\\s+1\\s+DE\\s+P$", ap2))
}

# ---- 4. error_reason_compute 酶依赖分类（经 ISO_ENZYME 切换）----
old_env <- Sys.getenv("ISO_ENZYME", unset = "")
on.exit(if (nzchar(old_env)) Sys.setenv(ISO_ENZYME = old_env) else Sys.unsetenv("ISO_ENZYME"), add = TRUE)
# trypsin（历史口径）：左翼 K/R 任一 → 不报 N 端非特异
Sys.setenv(ISO_ENZYME = "trypsin")
e1 <- error_reason_compute("different_modification", "PEPTIDE", "prot1/prot2;(K)PEPTIDE(R)/(A)PEPTIDE(K)", "")
ck("trypsin：任一左翼 K/R → 无 N 端非特异（历史口径）",
   !grepl("n_terminal_non_specific_cleavage", e1) && grepl("modification_correction", e1))
e2 <- error_reason_compute("different_modification", "PEPTIDE", "prot1;(A)PEPTIDE(A)", "")
ck("trypsin：左翼非 K/R → N 端非特异",
   grepl("n_terminal_non_specific_cleavage", e2))
# lysc：左翼 K 特异、R 不特异
Sys.setenv(ISO_ENZYME = "lysc")
e3 <- error_reason_compute("different_modification", "PEPTIDE", "prot1;(K)PEPTIDE(K)", "")
ck("lysc：左翼 K → 特异", !grepl("n_terminal_non_specific_cleavage", e3))
e4 <- error_reason_compute("different_modification", "PEPTIDE", "prot1;(R)PEPTIDE(K)", "")
ck("lysc：左翼 R → N 端非特异", grepl("n_terminal_non_specific_cleavage", e4))
# lysn：看首残基——rev_pep 首残基 K → 特异
Sys.setenv(ISO_ENZYME = "lysn")
e5 <- error_reason_compute("different_modification", "KPEPTIDE", "prot1;(A)KPEPTIDE(K)", "")
ck("lysn：首残基 K → 特异（左翼任意）", !grepl("n_terminal_non_specific_cleavage", e5))
e6 <- error_reason_compute("different_modification", "APEPTIDE", "prot1;(K)APEPTIDE(K)", "")
ck("lysn：首残基 A → N 端非特异（即使左翼 K）", grepl("n_terminal_non_specific_cleavage", e6))
# none：恒特异
Sys.setenv(ISO_ENZYME = "none")
e7 <- error_reason_compute("different_modification", "APEPTIDE", "prot1;(A)APEPTIDE(A)", "")
ck("none：恒特异（无 N 端非特异）", !grepl("n_terminal_non_specific_cleavage", e7))

cat(sprintf("=== enzyme_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
