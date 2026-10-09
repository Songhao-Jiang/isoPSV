# mgf_to_zml_test.R — MGF→mzML 本地化转换回归（v1.0，纯 R）
# 覆盖：base64 标准向量、谱数与结构、m/z 与 intensity 往返解码逐值一致、
# charge/rt/title/pepmass 保留、indexList 字节偏移正确、CV 词齐备。
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

# 1) base64 标准向量（RFC 4648 测试向量：0x00..0x0F）
ck("base64 标准向量 0x00..0x0F",
   identical(base64_encode_raw(as.raw(0:15)), "AAECAwQFBgcICQoLDA0ODw=="))

# 测试侧解码器（与编码器互逆，用于往返对拍）
b64_decode_raw <- function(s) {
  s <- gsub("=", "", s)
  vals <- match(strsplit(s, "")[[1]], b64_alphabet) - 1L
  m6 <- sapply(vals, function(v) as.integer(intToBits(v)[1:6]))  # LSB 在前
  m6 <- m6[6:1, , drop = FALSE]                                   # 翻 MSB 在前
  bits <- as.logical(as.vector(m6))
  nb <- (length(bits) %/% 8L) * 8L
  bits <- bits[seq_len(nb)]
  m8 <- matrix(bits, nrow = 8L)
  m8 <- m8[8:1, , drop = FALSE]
  packBits(as.integer(as.vector(m8)), "raw")
}

# 2) 合成 MGF → 转换 → 结构断言
td <- tempfile("mzm_"); dir.create(td)
mgf_f <- file.path(td, "t.mgf"); mzml_f <- file.path(td, "t.mzML")
mz1 <- c(100.1, 200.2, 300.3); in1 <- c(1000, 2000, 3000)
mz2 <- c(400.4, 500.5); in2 <- c(4000, 5000)
writeLines(c(
  "BEGIN IONS", "TITLE=raw1.1.1.2.0.dta", "PEPMASS=456.789 12345", "CHARGE=2+",
  "RTINSECONDS=12.34", "100.1 1000", "200.2 2000", "300.3 3000", "END IONS",
  "BEGIN IONS", "TITLE=raw1.2.2.3.0.dta", "PEPMASS=789.123 678", "CHARGE=3+",
  "RTINSECONDS=56.78", "400.4 4000", "500.5 5000", "END IONS"), mgf_f)
n_out <- mgf_to_mzml(mgf_f, mzml_f)
ck("转换返回谱数 2", n_out == 2L)
txt <- paste(readLines(mzml_f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
ck("mzML 谱表计数", grepl('<spectrumList count="2">', txt, fixed = TRUE))
ck("谱 id 与 index 正确",
   grepl('id="scan=1" defaultArrayLength="3"', txt, fixed = TRUE) &&
   grepl('id="scan=2" defaultArrayLength="2"', txt, fixed = TRUE))
ck("CV 词齐备（m/z/intensity/64bit/no compression）",
   all(vapply(c("MS:1000514", "MS:1000515", "MS:1000523", "MS:1000576",
                "MS:1000016", "MS:1000744", "MS:1000041"),
              function(t) grepl(t, txt, fixed = TRUE), logical(1))))
# 实弹验收（2026-09-26，deepseek v1.0 报告）：谱级 ms level 必须为 2
# （MSn），不得标 MS1；头部链四件齐备
ck("谱级 ms level=2 + MSn spectrum（不再误标 MS1）",
   grepl('accession="MS:1000511" name="ms level" value="2"', txt, fixed = TRUE) &&
   grepl('accession="MS:1000580" name="MSn spectrum" value=""', txt, fixed = TRUE) &&
   !grepl("MS1 spectrum", txt, fixed = TRUE))
ck("总离子流 TIC 写入（MS:1000285）",
   grepl('accession="MS:1000285" name="total ion current"', txt, fixed = TRUE))
# P1-1（2026-09-26，工作站 v1.0 复测）：每谱必须带 centroid spectrum
# （MS:1000127）——缺它 Sage panic（contains profile data）、MetaMorpheus 计
# 0 二级谱；补上后 Sage 实测 1,457 PSM
ck("每谱 centroid spectrum（MS:1000127 计数＝谱数 2）",
   length(regmatches(txt,
     gregexpr('accession="MS:1000127" name="centroid spectrum" value=""', txt,
              fixed = TRUE))[[1]]) == 2L)
ck("头部链齐备（software/instrumentConfiguration/defaultInstrument/dataProcessing）",
   all(vapply(c("softwareList", "instrumentConfigurationList",
                'defaultInstrumentConfigurationRef="IC1"', "dataProcessingList"),
              function(t) grepl(t, txt, fixed = TRUE), logical(1))))
ck("charge/rt/title/pepmass 保留",
   grepl('value="2"', txt, fixed = TRUE) && grepl('name="scan start time" value="12.34"', txt, fixed = TRUE) &&
   grepl('value="raw1.1.1.2.0.dta"', txt, fixed = TRUE) &&
   grepl('value="456.789"', txt, fixed = TRUE))

# 3) 往返解码：谱 1 的 m/z 与 intensity 逐值一致
bins <- regmatches(txt, gregexpr("<binary>[^<]+</binary>", txt))[[1]]
b64_1_mz <- sub("<binary>|</binary>", "", bins[1])
b64_1_in <- sub("<binary>|</binary>", "", bins[2])
dec_mz <- readBin(b64_decode_raw(b64_1_mz), "double", n = 3L, size = 8L, endian = "little")
dec_in <- readBin(b64_decode_raw(b64_1_in), "double", n = 3L, size = 8L, endian = "little")
ck("谱 1 m/z 往返逐值一致", identical(dec_mz, mz1))
ck("谱 1 intensity 往返逐值一致", identical(dec_in, in1))
b64_2_mz <- sub("<binary>|</binary>", "", bins[3])
dec2 <- readBin(b64_decode_raw(b64_2_mz), "double", n = 2L, size = 8L, endian = "little")
ck("谱 2 m/z 往返逐值一致", identical(dec2, mz2))

# 4) indexList 偏移：offset 文本值 == 文件中 '<spectrum index=' 的字节偏移
spec_pos <- gregexpr("<spectrum index=", txt, fixed = TRUE)[[1]]
off_m <- regmatches(txt, gregexpr('<offset idRef="scan=[0-9]+"[^>]*>[0-9]+</offset>', txt))[[1]]
offs <- as.integer(sub("^.*>([0-9]+)</offset>$", "\\1", off_m))
exp_offs <- vapply(seq_along(spec_pos), function(k)
  nchar(substr(txt, 1, spec_pos[k] - 1), type = "bytes"), numeric(1))
ck("indexList 字节偏移与谱元素位置一致", identical(offs, as.integer(exp_offs)))
unlink(td, recursive = TRUE)

cat(sprintf("=== mgf_to_mzml_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
