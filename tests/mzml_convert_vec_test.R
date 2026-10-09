# mzml_convert_vec_test.R — v1.0 ③mzML 转换向量化回归（2026-09-29）
# 覆盖：base64_encode_fast 与 base64_encode_raw 在随机 raw（含 0/1/2 字节尾组、
# 8/24 字节对齐）上逐字节一致；合成 MGF（含边角谱）新旧转换产物逐字节一致；
# 非良构 MGF（BEGIN 无 END）回落 _base 且产物一致。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
suppressMessages(source(file.path(cd, "compat.R")))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}
set.seed(268)

# 1) 编码器对拍：尺寸覆盖 0、1、2、3、4、5、7、8k、24k、随机 0-600
sizes <- c(0L, 1L, 2L, 3L, 4L, 5L, 7L, 8L, 16L, 24L, 48L, 1600L,
           sample(0:600, 40L))
same <- vapply(sizes, function(k) {
  r <- as.raw(sample(0:255, k, replace = TRUE))
  identical(base64_encode_fast(r), base64_encode_raw(r))
}, logical(1))
ck("base64 快编码器与旧实现逐字节一致（53 组尺寸）", all(same))

# 2) 合成 MGF → 新旧转换产物逐字节一致（含边角谱）
td <- tempfile("mzv_"); dir.create(td)
mgf <- file.path(td, "t.mgf")
writeLines(c(
  "BEGIN IONS",
  "TITLE=raw1.1.1.2.0.dta",
  "PEPMASS=456.789 12345",
  "CHARGE=2+",
  "RTINSECONDS=12.34",
  "100.1 1000",
  "200.2 2000 999",          # 三 token：第三忽略
  "300.3 3000",
  "END IONS",
  "7 7",                     # 块外数字行：忽略
  "BEGIN IONS",
  "TITLE=a<b>&\"q\".2.2.3.0.dta",
  "PEPMASS=789.123 678",
  "CHARGE=3+",
  "RTINSECONDS=56.78",
  "END IONS",                # 无峰谱
  "BEGIN IONS",
  "TITLE=raw1.3.3.2.0.dta",  # 缺 PEPMASS/CHARGE/RT（默认 0/2/0）
  "400.4 4000",
  "abc 100",                 # mz 非数值：丢弃
  "500.5",                   # 单 token：丢弃
  "END IONS",
  "BEGIN IONS",
  "TITLE=raw1.4.4.2.0.dta",  # 参数行重复：最后一次生效
  "PEPMASS=111.1 1",
  "PEPMASS=222.2 2",
  "CHARGE=2+",
  "CHARGE=4+",
  "600.6 600",
  "END IONS"), mgf, useBytes = TRUE)
o1 <- file.path(td, "new.mzML"); o2 <- file.path(td, "old.mzML")
n1 <- mgf_to_mzml(mgf, o1); n2 <- mgf_to_mzml_base(mgf, o2)
ck("合成 MGF 谱数一致（4）", identical(n1, 4L) && identical(n1, n2))
ck("合成 MGF 产物逐字节一致", identical(readLines(o1, warn = FALSE),
                                        readLines(o2, warn = FALSE)))
x <- paste(readLines(o1, warn = FALSE), collapse = "\n")
ck("重复参数行取最后值（PEPMASS=222.2 / charge=4）",
   grepl('selected ion m/z" value="222.2"', x, fixed = TRUE) &&
   grepl('charge state" value="4"', x, fixed = TRUE))
ck("无峰谱 TIC=0 且空 binary",
   grepl('total ion current" value="0"', x, fixed = TRUE) &&
   grepl('defaultArrayLength="0"', x, fixed = TRUE))

# 3) 非良构（BEGIN 无 END / 双 END）→ 回落 _base，产物一致
mgf2 <- file.path(td, "bad.mgf")
writeLines(c("BEGIN IONS", "TITLE=b.1.1.2.0.dta", "100.1 5",
             "BEGIN IONS", "TITLE=b.2.2.2.0.dta", "200.2 6", "END IONS"),
           mgf2, useBytes = TRUE)
o3 <- file.path(td, "bad_new.mzML"); o4 <- file.path(td, "bad_old.mzML")
n3 <- mgf_to_mzml(mgf2, o3); n4 <- mgf_to_mzml_base(mgf2, o4)
ck("非良构回落 _base 且产物一致", identical(n3, n4) &&
   identical(readLines(o3, warn = FALSE), readLines(o4, warn = FALSE)))

cat("==== mzml_convert_vec_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
