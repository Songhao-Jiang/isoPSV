# y_vector_fix_test.R — 验证 y 向量 reverse 修复后 predict 解析正确（自包含版 2026-09-23）
# 自包含：合成 pDeep2 存储序夹具（y+1 按倒序存储），不依赖任何外部数据文件；
# 路径按 --file= 相对定位（可随发行包任意布置，无开发机硬编码）。
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
code <- file.path(dirname(sd), "code")
if (!dir.exists(code)) code <- file.path(sd, "..", "code")
source(file.path(code, "compat.R"))
source(file.path(code, "ReadMGF_script.R"))

# ---- 合成夹具：pDeep2 存储序（y+1 段 = y6..y1 倒序存放；修复前解析错位）----
fx <- file.path(tempdir(), "predict_fixture.txt")
writeLines(c(
  "BEGIN IONS",
  "TITLE=MFRPAQR||2",
  "CHARGE=2",
  "pepinfo=MFRPAQR||2",
  "b+1=0.00000,0.10006,0.30881,0.00008,0.12148,0.32354",
  "b+2=0.00000,0.00000,0.00000,0.00000,0.00000,0.00162",
  "y+1=0.00000,0.00161,0.99248,0.02581,0.02255,0.42056",
  "y+2=0.06541,0.08913,0.03674,0.00000,0.00000,0.00000",
  "END IONS"), fx, useBytes = TRUE)

pred <- read_predict_txt(fx)
key <- "MFRPAQR||2"
if (!key %in% pred$keys) stop("MFRPAQR not found in predict keys")
pv <- pred$intensities[[pred$idx[key]]]

cat("总向量长度:", length(pv), "\n")
# MFRPAQR 7 残基 → 6 b/y 离子；修复后顺序: b+1(6) + b+2(6) + rev(y+1)(6) + rev(y+2)(6)
y1_fixed <- pv[13:18]  # 修复后 y+1 段 = [y1, y2, y3, y4, y5, y6]
cat("y+1 fixed:", sprintf("%.5f", y1_fixed), "\n")
cat("expected:  y1=0.421 y2=0.023 y3=0.026 y4=0.992 y5=0.002 y6=0.000\n")

pass <- TRUE
if (abs(y1_fixed[1] - 0.421) < 0.01) {
  cat("PASS: y1 = 0.421 (was 0.000 before fix)\n")
} else { cat("FAIL: y1 =", y1_fixed[1], "\n"); pass <- FALSE }
if (abs(y1_fixed[4] - 0.992) < 0.01) {
  cat("PASS: y4 = 0.992 (was at wrong position before fix)\n")
} else { cat("FAIL: y4 =", y1_fixed[4], "\n"); pass <- FALSE }
if (length(pv) == 24L) {
  cat("PASS: total vector length 24 (4 charge-type segments x 6)\n")
} else { cat("FAIL: length", length(pv), "\n"); pass <- FALSE }
if (abs(pv[1] - 0.00000) < 1e-4 && abs(pv[7] - 0.00000) < 1e-4) {
  cat("PASS: b+1 / b+2 leading zero preserved\n")
} else { cat("FAIL: b leading zeros\n"); pass <- FALSE }
if (abs(y1_fixed[6] - 0.000) < 0.01 && abs(y1_fixed[5] - 0.002) < 0.01) {
  cat("PASS: y+1 tail (y5,y6) at expected positions\n")
} else { cat("FAIL: y+1 tail\n"); pass <- FALSE }

file.remove(fx)
cat("=== y 向量 reverse 修复验证", if (pass) "全部 PASS" else "FAIL", "===\n")
if (!pass) quit(status = 1)
