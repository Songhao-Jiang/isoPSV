# predictor_adapter_test.R — 多理论谱图预测器适配层回归（v1.0，2026-09-26；
# v1.0 按 GLM N-2/N-11/N-12 与 MiniMax P0-3 更新口径）
# 覆盖 pred_ms2pip.R 的纯函数（无需 ms2pip 安装）：修饰名映射、peprec 构造、
# 预测块布局（b 正序 / y 位点顺序）、合成 CSV 组装、消费端键契约、电荷 3+ 计数。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
ad <- file.path(sd, "..", "code", "predictors", "pred_ms2pip.R")
if (!file.exists(ad)) ad <- file.path(sd, "..", "..", "code", "predictors", "pred_ms2pip.R")
if (!file.exists(ad)) {
  cat("SKIP : pred_ms2pip.R 缺失（发行包路径差异），跳过对应断言\n")
  quit(status = 0L)
}
source(ad)

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 1) pFind 修饰 → peprec 修饰串
ck("pFind 修饰转 peprec",
   identical(pfind_mods_to_peprec("3,Oxidation[M];0,Acetyl[AnyN-term];"),
             "3|Oxidation;0|Acetyl"))
ck("未知修饰名丢弃",
   identical(pfind_mods_to_peprec("9,Unknown[X];", warn = FALSE), "-"))
ck("未修饰肽写 -（N-11）", identical(pfind_mods_to_peprec(""), "-"))

# 2) peprec 构造（含 N-2 的 orig_mod 列：修饰原串去尾分号、保留 [M]）
df <- data.frame(peptide = "PEPTIDEK", modification = "4,Oxidation[M];", charge = "2",
                 stringsAsFactors = FALSE)
pr <- build_peprec(df)
ck("peprec 构造", nrow(pr) == 1 && pr$peptide[1] == "PEPTIDEK" && pr$charge[1] == "2")
ck("orig_mod 保留修饰原串（去尾分号）", identical(pr$orig_mod[1], "4,Oxidation[M]"))

# 3) 块布局：b 正序、y 写入时反转（read_predict_txt 读取时再 reverse）
bl <- block_lines("PEPTIDEK||2", 1:7, 11:17, 21:27, 31:37)
ck("块 7 行且 TITLE 正确", length(bl) == 7 && identical(bl[2], "TITLE=PEPTIDEK||2"))
b1 <- as.numeric(strsplit(sub("b+1=", "", bl[startsWith(bl, "b+1=")][1], fixed = TRUE),
                          ",", fixed = TRUE)[[1]])
ck("b 行正序", identical(b1, as.numeric(1:7)))
y1 <- as.numeric(strsplit(sub("y+1=", "", bl[startsWith(bl, "y+1=")][1], fixed = TRUE),
                          ",", fixed = TRUE)[[1]])
ck("y 行按位点顺序（写入时反转）", identical(y1, as.numeric(27:21)))

# 4) 合成 CSV 组装（模拟 ms2pip 输出；P0-3 夹具含电荷 3+ 行）
td <- tempfile("pfm_"); dir.create(td)
csv <- file.path(td, "pred.csv")
writeLines(c("spec_id,charge,ion,b_ion,b_ion_mz,b_int,y_ion,y_ion_mz,y_int",
             "p000001,1,1,b1,100,10,y1,700,20",
             "p000001,1,2,b2,200,30,y2,800,40",
             "p000001,1,3,b3,300,50,y3,900,60",
             "p000001,2,1,b1,100,11,y1,700,21",
             "p000001,2,2,b2,200,31,y2,800,41",
             "p000001,2,3,b3,300,51,y3,900,61",
             "p000001,3,1,b1,100,12,y1,700,22",
             "p000001,3,2,b2,200,32,y2,800,42",
             "p000001,3,3,b3,300,52,y3,900,62"), csv)
parsed <- parse_ms2pip_csv(csv)
ck("解析返回 by_id 命名列表（N-12 split 路径）",
   identical(names(parsed$by_id), "p000001"))
ck("电荷 3+ 行计数上报（P0-3）", parsed$charge3plus_ids == 1L)
peprec <- data.frame(spec_id = "p000001", modifications = "", peptide = "PEPT",
                     charge = "3", orig_mod = "", stringsAsFactors = FALSE)
built <- build_predict_txt(peprec, parsed$by_id)
lines <- built$lines
ck("合成 CSV → 预测块 7 行（电荷3+ 谱不整谱跳过）", length(lines) == 7)
b2 <- as.numeric(strsplit(sub("b+2=", "", lines[4], fixed = TRUE), ",", fixed = TRUE)[[1]])
ck("b+2 取电荷 2 行", identical(b2, as.numeric(c(11, 31, 51))))
y2 <- as.numeric(strsplit(sub("y+2=", "", lines[6], fixed = TRUE), ",", fixed = TRUE)[[1]])
ck("y+2 反转后为位点顺序", identical(y2, as.numeric(c(61, 41, 21))))
unlink(td, recursive = TRUE)

# 5) N-2 键契约：适配器 TITLE 键 = <seq>|<修饰原串>|<charge>，与消费端
#    ReadMGF_script.R read_predict_txt 的 sub(";+$","",...) 口径逐字节一致
peprec_mod <- data.frame(spec_id = "p000002", modifications = "4|Oxidation",
                         peptide = "PEPT", charge = "2",
                         orig_mod = "4,Oxidation[M]", stringsAsFactors = FALSE)
parsed2 <- parsed$by_id["p000001"]; names(parsed2) <- "p000002"
lines2 <- build_predict_txt(peprec_mod, parsed2)$lines
ck("修饰肽键与消费端逐字节一致（保留 [M]、无尾分号）",
   length(lines2) == 7 && identical(lines2[2], "TITLE=PEPT|4,Oxidation[M]|2"))

# 6) 跳过计数：无产物 / 离子数不符均上报（不再静默）
bad <- data.frame(spec_id = "p000009", modifications = "-", peptide = "PEPT",
                  charge = "2", orig_mod = "", stringsAsFactors = FALSE)
pbad <- list(p000009 = list(b1 = 1, b2 = 1, y1 = 1, y2 = 1))
ck("离子数不符计数上报", build_predict_txt(bad, pbad)$skip_len == 1L)

cat(sprintf("=== predictor_adapter_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
