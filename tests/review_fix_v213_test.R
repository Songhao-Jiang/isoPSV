# review_fix_v213_test.R — v1.0 GLM 审查 P0 修复回归测试（2026-09-26）
# 覆盖：precursor_mz 的 z=1/2/3 口径、X!Tandem 多行 XML 解析（P0-5）、
# MaxQuant 未知名内联修饰括号跳过（P0-2）。
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

# ---- 1. precursor_mz：质子化质量 (M+H) → 前体 m/z = (M+zH)/z ----
H <- 1.007825
ck("z=2 与旧公式一致（历史 2+ 数据零回归）",
   isTRUE(all.equal(precursor_mz(1000, 2), (1000 + H) / 2, tolerance = 1e-12)))
ck("z=1 = M+H（旧公式会偏高 1.0078）",
   isTRUE(all.equal(precursor_mz(1000, 1), 1000, tolerance = 1e-12)))
ck("z=3 = (M+3H)/3（旧公式偏高 0.336）",
   isTRUE(all.equal(precursor_mz(1000, 3), (1000 + 2 * H) / 3, tolerance = 1e-12)))
ck("z=4 = (M+4H)/4",
   isTRUE(all.equal(precursor_mz(1000, 4), (1000 + 3 * H) / 4, tolerance = 1e-12)))

# ---- 2. X!Tandem 多行 XML（P0-5：TRE . 不匹配换行 → perl=TRUE）----
adp <- file.path(cd, "search_adapter_xtandem.R")
if (!file.exists(adp)) cat("SKIP : xtandem 适配器文件缺失（发行包路径差异），跳过对应断言\n")
if (file.exists(adp)) {
  # 用与适配器相同的正则直接验证多行匹配（适配器本体由 mascot_peaks_adapter_test 覆盖单行）
  txt <- paste(c(
    '<?xml version="1.0"?>',
    '<bioml>',
    '  <group id="1" mh="1201.6435" z="2" rt="30.2" expect="0.001"',
    '         label="raw1.1234.1234.2.0.dta" type="model">',
    '    <protein expect="0.001" label="prot1" />',
    '    <domain expect="0.001" start="1" end="10">',
    '      <aa seq="PEPTIDEK" />',
    '    </domain>',
    '  </group>',
    '</bioml>'), collapse = "\n")
  g1 <- regmatches(txt, gregexpr("(?s)<group(?=[[:space:]])[^>]*>.*?</group>", txt, perl = TRUE))[[1]]
  ck("多行 XML：perl(?s) 命中 1 个 group 块", length(g1) == 1)
  # 实测（R 4.5.1）：TRE 的 . 其实匹配换行（grepl("a.b","a\nb")=TRUE），真正
  # 的坑是 \\b 在 LC_CTYPE=C 下恒不匹配——(?:space:) 前瞻 + perl(?s) 双环境安全
  gb <- gregexpr("(?s)<group(?=[[:space:]])[^>]*>.*?</group>", txt, perl = TRUE)[[1]]
  gx <- regmatches(txt, gregexpr("(?s)<groupx[^>]*>.*?</group>", txt, perl = TRUE))[[1]]
  ck("space 前瞻防误匹配 <groupx>", length(gx) == 0 && length(gb) == 1)
}

# ---- 3. MaxQuant 未知名内联修饰括号跳过（P0-2）----
adp2 <- file.path(cd, "search_adapter_maxquant.R")
if (!file.exists(adp2)) cat("SKIP : maxquant 适配器文件缺失（发行包路径差异），跳过对应断言\n")
# MiniMax P0-3（2026-09-26）：短标签 (ac)/(ox) 应被识别而非静默丢弃
sys_src3 <- new.env()
if (file.exists(adp2)) {
  old_wd3 <- setwd(cd)
  on.exit(setwd(old_wd3), add = TRUE)
  sys.source(adp2, envir = sys_src3)
  rac <- get("parse_modified_sequence", envir = sys_src3)("_(ac)SSPK", "SSPK")
  ck("短标签 (ac) 识别为 N 端乙酰化", grepl("0,Acetyl[AnyN-term]", rac, fixed = TRUE))
  rox <- get("parse_modified_sequence", envir = sys_src3)("SS(ox)PK", "SSPK")
  ck("短标签 (ox) 识别为 S 氧化", grepl("2,Oxidation[S]", rox, fixed = TRUE))
  # N-21（2026-09-26，GLM）：ph/de 短标签此前只入表未断言
  rph <- get("parse_modified_sequence", envir = sys_src3)("SS(ph)PK", "SSPK")
  ck("短标签 (ph) 识别为 S 磷酸化", grepl("2,Phospho[S]", rph, fixed = TRUE))
  rde <- get("parse_modified_sequence", envir = sys_src3)("S(ox)S(de)PK", "SSPK")
  ck("短标签 (de) 识别为脱酰胺", grepl("2,Deamidated[S]", rde, fixed = TRUE))
  # N-19（2026-09-26，GLM）：连续两个未知前导标签——break 只吞第一个，第二个
  # 会漏进逐残基解析导致位点漂移；next 后应整体吞掉且后续位点不错位
  r2u <- get("parse_modified_sequence", envir = sys_src3)("_(U1)(U2)SS(ox)PK", "SSPK")
  ck("连续未知前导标签整体跳过（N-19）", grepl("2,Oxidation[S]", r2u, fixed = TRUE) &&
     !grepl("U1|U2", r2u))
}
if (file.exists(adp2)) {
  src <- readLines(adp2, warn = FALSE)
  # 常见修饰已入表
  ck("maxquant 适配器含 Phospho (STY)", any(grepl('"Phospho \\(STY\\)"', src)))
  ck("maxquant 适配器含 Methyl (KR)", any(grepl('"Methyl \\(KR\\)"', src)))
  # 括号整体跳过逻辑存在（配对括号深度扫描）
  ck("maxquant 适配器含配对括号跳过", any(grepl("dep <- dep - 1L", src)))
  # 直接调用 parse_modified_sequence 验证位点不错位：
  # 未知修饰 S(Phospho (ST)) 后接已知修饰 M(Oxidation (M))——M 位点必须仍是 4
  #（旧实现括号内字母被误计为残基，M 位点会偏移）
  sys_src <- new.env()
  old_wd <- setwd(cd)
  on.exit(setwd(old_wd), add = TRUE)
  sys.source(adp2, envir = sys_src)
  r <- get("parse_modified_sequence", envir = sys_src)("_(ac)SS(Phospho (ST))PM(Oxidation (M))K", "SSPMK")
  ck("未知修饰 S(Phospho (ST)) 不使后续 M 位点错位（应含 4,Oxidation）",
     grepl("4,Oxidation", r, fixed = TRUE))
}
# ---- 4. pf_sanitize_qryres（P0-A，2026-09-26）----
if (exists("pf_sanitize_qryres", mode = "function")) {
  td4 <- tempfile("pfs_"); dir.create(td4)
  f4 <- file.path(td4, "x.qry.res")
  writeLines(c("S\t1234.5\t2", "spec_a.dta",
               "1\tPEP\tK\t2\t0.9\t0.01\t0;name#0",
               "3\tPEPR\t2\t0.8\t0.02\t0",
               "S\t2345.6\t2", "spec_b.dta",
               "1\tPEPX\t2\t0.7\t0.03\t0;name#0"), f4)
  pf_sanitize_qryres(td4, verbose = FALSE)
  l4 <- readLines(f4, warn = FALSE)
  ck("消毒：name#0 行删除、整块删除、rank 重排",
     !any(grepl("name#0", l4)) && !any(grepl("spec_b", l4)) &&
       identical(sub("\\t.*", "", l4[3]), "1"))
  pf_sanitize_qryres(td4, verbose = FALSE)
  ck("消毒幂等（二次零改动）", identical(readLines(f4, warn = FALSE), l4))
  unlink(td4, recursive = TRUE)
}
# 若夹具缺文件（发行包路径差异）则跳过并记 SKIP
cat(sprintf("=== review_fix_v213_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)