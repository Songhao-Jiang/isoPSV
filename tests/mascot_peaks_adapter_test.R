# mascot_peaks_adapter_test.R — Mascot/PEAKS/X!Tandem 适配器规格测试（2026-09-24）
# 覆盖：Mascot CSV 导出（Matrix Science csv_headers.html 列规范）、PEAKS DB
# search 肽段 CSV（PEAKS 手册列规范）与 X!Tandem XML（tandem.xml <group> 块）
# → 交换格式的转换正确性。
# 自包含：夹具按厂商文档格式合成，不依赖外部文件；--file= 相对定位。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- file.path(dirname(sd), "code")
if (!dir.exists(cd)) cd <- file.path(sd, "..", "code")
source(file.path(cd, "compat.R"))
source(file.path(cd, "search_adapter_mascot.R"))
source(file.path(cd, "search_adapter_peaks.R"))
source(file.path(cd, "search_adapter_xtandem.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, ok) {
  cat(ifelse(ok, "PASS", "FAIL"), ":", label, "\n")
  if (ok) npass <<- npass + 1L else nfail <<- nfail + 1L
}

td <- tempfile("mpa_"); dir.create(td)

# ============================================================================
# 1. Mascot CSV 导出夹具（含逗号字段、rank 过滤、双修饰位点、N 端修饰、固定修饰注入）
# ============================================================================
mascot_lines <- c(
  "prot_hit_num,prot_acc,prot_desc,prot_score,prot_matches,prot_matches_sig,prot_sequences,prot_sequences_sig,prot_cover,prot_len,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_mr,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_miss,pep_score,pep_homol,pep_ident,pep_expect,pep_res_before,pep_seq,pep_res_after,pep_var_mod,pep_var_mod_pos,pep_summed_mod_pos,pep_num_match,pep_scan_title",
  "1,sp|P12345|GENE1,\"GENE1, protein description with comma\",58,5,4,3,2,12,200,5,1,1,1,601.8324,1201.6502,2,1201.6435,0.0067,10,21,0,58.1,45,45,1.2e-05,K,AMCLPEPTIDEK,R,Oxidation (M),010000000000,,10,raw1.1234.1234.2.0.dta",
  "1,sp|P12345|GENE1,\"GENE1, protein description with comma\",58,5,4,3,2,12,200,5,2,0,0,601.8331,1201.6516,2,1201.6435,0.0081,10,21,0,31.2,45,45,0.002,K,AMKLDEPTIDEK,R,,000000000000,,5,raw1.1234.1234.2.0.dta",
  "2,sp|Q67890|GENE2,GENE2 protein,44,3,3,2,2,9,160,9,1,1,1,700.4000,1458.7842,2,1458.7776,0.0066,2,16,0,44.0,40,40,0.001,R,ACETYLAMEPTIDEK,S,\"Acetyl (Protein N-term);Oxidation (M);Phospho (ST)\",100000020000000,000000030000000,8,raw2.2222.2222.2.0.dta",
  "3,sp|Q11111|GENE3,GENE3 protein,20,1,1,1,1,8,100,12,1,1,0,500.2500,1038.4769,2,1038.4703,0.0066,1,10,0,20.0,30,30,0.5,R,ACDEPTIDEK,R,Carbamidomethyl (C),0100000000,,4,raw3.3333.3333.3.0.dta"
)
mcsv <- file.path(td, "mascot_export.csv")
writeLines(mascot_lines, mcsv, useBytes = TRUE)
mout <- file.path(td, "mascot.spectra")
convert_mascot_csv(mcsv, mout, "Carbamidomethyl (C)")
m <- read.delim(mout, sep = "\t", check.names = FALSE, na.strings = character(0))

ck("Mascot rank 过滤：每谱仅 top-1 行", nrow(m) == 3L)
ck("Mascot 标题直通 File_Name",
   identical(m$File_Name, c("raw3.3333.3333.3.0.dta",
                            "raw1.1234.1234.2.0.dta",
                            "raw2.2222.2222.2.0.dta")) ||
   setequal(m$File_Name, c("raw1.1234.1234.2.0.dta", "raw2.2222.2222.2.0.dta",
                           "raw3.3333.3333.3.0.dta")))
m1 <- m[m$File_Name == "raw1.1234.1234.2.0.dta", ]
ck("Mascot Oxidation(M) 位点 + 固定 Carbamidomethyl 注入",
   m1$Modification == "2,Oxidation[M];3,Carbamidomethyl[C];")
m2 <- m[m$File_Name == "raw2.2222.2222.2.0.dta", ]
ck("Mascot N 端 Acetyl -> 0,Acetyl[AnyN-term] + 同一位点双修饰（summed_mod_pos）",
   m2$Modification == "0,Acetyl[AnyN-term];2,Carbamidomethyl[C];8,Oxidation[M];8,Phospho[M];")
m3 <- m[m$File_Name == "raw3.3333.3333.3.0.dta", ]
ck("Mascot 可变+固定同修饰去重合并", m3$Modification == "2,Carbamidomethyl[C];")
ck("Mascot 蛋白号 pFind 拼接", m1$Proteins == "sp|P12345|GENE1/")
ck("Mascot 电荷透传", identical(as.character(m1$Charge), "2"))

# ============================================================================
# 2. PEAKS DB search 肽段 CSV 夹具（内联修饰、-10lgP 去重、标题构建、Accession）
# ============================================================================
peaks_lines <- c(
  'Peptide,-10lgP,Mass,Length,ppm,m/z,z,RT,Area,Fraction,Scan,Source File,Accession,PTM,AScore,Found By,Id',
  'AM(+15.99)CLPEPTIDEK,45.2,1201.65,11,3.1,601.83,2,30.2,1e6,3,1234,raw1.raw,sp|P12345|GENE1,M(+15.99):1,20,DB,1',
  'AMCLPEPTIDEK,25.1,1201.64,11,3.0,601.82,2,30.2,1e6,3,1234,raw1.raw,sp|P12345|GENE1,,,DB,2',
  'ACETYLATEPTIDEK,40.0,1398.78,13,2.5,700.40,2,31.0,8e5,4,2222,raw2.raw,sp|Q67890|GENE2,Acetyl (Protein N-term):1,15,DB,3',
  'AD(+42.01)VLPEPTIDEK,38.0,1243.66,11,3.2,622.84,2,29.5,7e5,5,3333,raw3.raw,sp|Q11111|GENE3,Acetylation:1,12,DB,4'
)
pcsv <- file.path(td, "peaks_peptide.csv")
writeLines(peaks_lines, pcsv, useBytes = TRUE)
pout <- file.path(td, "peaks.spectra")
convert_peaks_csv(pcsv, pout)
p <- read.delim(pout, sep = "\t", check.names = FALSE, na.strings = character(0))

ck("PEAKS -10lgP 最佳每谱去重", nrow(p) == 3L)
p1 <- p[p$Sequence == "AMCLPEPTIDEK", ]
ck("PEAKS 内联 M(+15.99) -> 2,Oxidation[M];",
   p1$Modification == "2,Oxidation[M];")
ck("PEAKS 标题 raw.scan.scan.z.0.dta",
   p1$File_Name == "raw1.raw.1234.1234.2.0.dta")
p2 <- p[p$Sequence == "ACETYLATEPTIDEK", ]
ck("PEAKS 无内联修饰行保持空修饰", p2$Modification == "")
p3 <- p[p$Sequence == "ADVLPEPTIDEK", ]
ck("PEAKS 内联 D(+42.01) 乙酰化位点", p3$Modification == "2,Acetyl[D];")
ck("PEAKS Accession -> Proteins", p1$Proteins == "sp|P12345|GENE1/")
ck("PEAKS 电荷透传", identical(as.character(p1$Charge), "2"))

# ============================================================================
# 3. X!Tandem XML 夹具（<group> 块、内联修饰、N 端 n[delta]、label 去重）
# ============================================================================
xt_lines <- c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  '<bioml>',
  '<note label="output" type="model"/>',
  '<group id="1" mh="1201.6435" z="2" rt="30.2" expect="0.001" label="raw1.1234.1234.2.0.dta" type="model">',
  '  <protein expect="0.01" label="sp|P12345|GENE1" sumI="4.5">',
  '    <peptide start="10" end="21">',
  '      <domain start="10" end="21" expect="0.001" seq="AM[15.9949]CLPEPTIDEK" mh="1201.6435"/>',
  '    </peptide>',
  '  </protein>',
  '</group>',
  '<group id="2" mh="998.4788" z="2" rt="31.0" expect="0.05" label="raw2.2222.2222.2.0.dta" type="model">',
  '  <protein expect="0.06" label="sp|Q67890|GENE2">',
  '    <peptide start="1" end="10">',
  '      <domain start="1" end="10" expect="0.05" seq="n[42.010565]ACDEPTIDEK"/>',
  '    </peptide>',
  '  </protein>',
  '</group>',
  '<group id="3" mh="1300.0" z="3" rt="32.0" expect="0.5" label="raw2.2222.2222.2.0.dta" type="model">',
  '  <protein expect="0.6" label="sp|Q11111|GENE3">',
  '    <peptide start="1" end="11">',
  '      <domain start="1" end="11" expect="0.5" seq="AAAADEPTIDEK"/>',
  '    </peptide>',
  '  </protein>',
  '</group>',
  '</bioml>'
)
xout <- file.path(td, "xtandem.spectra")
writeLines(xt_lines, file.path(td, "tandem.xml"), useBytes = TRUE)
convert_xtandem_xml(file.path(td, "tandem.xml"), xout)
x <- read.delim(xout, sep = "\t", check.names = FALSE, na.strings = character(0))
ck("X!Tandem 每谱 top-1 且 label 去重（expect 最小优先）",
   nrow(x) == 2L && setequal(x$File_Name, c("raw1.1234.1234.2.0.dta",
                                             "raw2.2222.2222.2.0.dta")))
x1 <- x[x$File_Name == "raw1.1234.1234.2.0.dta", ]
ck("X!Tandem 内联 M[15.9949] -> 2,Oxidation[M];",
   x1$Sequence == "AMCLPEPTIDEK" && x1$Modification == "2,Oxidation[M];")
x2 <- x[x$File_Name == "raw2.2222.2222.2.0.dta", ]
ck("X!Tandem 前导 n[42.010565] -> 0,Acetyl[AnyN-term];",
   x2$Sequence == "ACDEPTIDEK" && x2$Modification == "0,Acetyl[AnyN-term];")
ck("X!Tandem 电荷透传", identical(as.character(x1$Charge), "2"))

# ============================================================================
# 4. 交换格式并入合并口径（三列最小契约：File_Name/Sequence/Modification）
# ============================================================================
ck("Mascot 产物三列契约", all(c("File_Name", "Sequence", "Modification") %in% names(m)))
ck("PEAKS 产物三列契约", all(c("File_Name", "Sequence", "Modification") %in% names(p)))

unlink(td, recursive = TRUE)
cat(sprintf("=== mascot_peaks_adapter_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
