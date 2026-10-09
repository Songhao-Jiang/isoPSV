# workstation_v262_test.R — v1.0 工作站反馈修复回归（2026-09-28）
# 覆盖：#14 引擎侧车候选口径、#16 python_exe step1 优先、#18 后处理复用判据
# （指纹三函数单元测试）、#19 读段默认顺序+旋钮+worker 进度、#21 失败醒目标记、
# 模板新键与版本号。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R"))) file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "boot.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# ---- 1. #18 指纹三函数单元测试 ----
tmp <- tempfile("pp_"); dir.create(tmp)
spf <- file.path(tmp, "pFind.spectra")
qr1 <- file.path(tmp, "a.qry.res"); qr2 <- file.path(tmp, "b.qry.res")
fa  <- file.path(tmp, "db.fasta")
writeLines("x", spf); writeLines("y", qr1); writeLines("z", qr2); writeLines(">a", fa)
mrk <- file.path(tmp, ".pfind_postprocess_done")
ck("无标记 → fresh=FALSE（需要跑）", !pfind_postprocess_fresh(tmp, fa, mrk))
pfind_postprocess_mark(tmp, fa, mrk)
ck("标记后 → fresh=TRUE（跳过）", pfind_postprocess_fresh(tmp, fa, mrk))
Sys.sleep(1.2)
writeLines("xx", spf)   # pFind.spectra 变化（重搜）
ck("产物变化 → fresh=FALSE", !pfind_postprocess_fresh(tmp, fa, mrk))
pfind_postprocess_mark(tmp, fa, mrk)
writeLines("yy", qr2)   # qry.res 变化
ck("qry.res 变化 → fresh=FALSE", !pfind_postprocess_fresh(tmp, fa, mrk))
pfind_postprocess_mark(tmp, fa, mrk)
unlink(fa)              # FASTA 变化（换库）
ck("FASTA 变化 → fresh=FALSE", !pfind_postprocess_fresh(tmp, fa, mrk))
fp_lines <- trimws(readLines(mrk, warn = FALSE))
ck("标记内容为 逐行 路径|mtime|size 格式", length(fp_lines) == 4L &&
   all(grepl("\\|[0-9]+\\|[0-9]+$", fp_lines)))
unlink(tmp, recursive = TRUE)

# ---- 2. 静态断言 ----
s2 <- paste(readLines(file.path(cd, "user_step2.R"), warn = FALSE), collapse = "\n")
rd <- paste(readLines(file.path(cd, "ReadMGF_script.R"), warn = FALSE), collapse = "\n")

# #14
ck("#14 签名带 cand_scope", grepl("label_mass = NULL, cand_scope = \"filtered\") {", s2, fixed = TRUE))
ck("#14 侧车含 cand_scope=", grepl('paste0("cand_scope=", cand_scope)', s2, fixed = TRUE))
ck("#14 调用点传 cand_scope = cand_scope", grepl("label_mass = label_mass_auto, cand_scope = cand_scope", s2, fixed = TRUE))

# #16
ck("#16 python_exe 走 cfg_get_preferred",
   grepl('cfg_get_preferred(base_cfg, cfg2, "python_exe", "")$value', s2, fixed = TRUE) &&
   grepl("grepl(\"[<>]\", py_hint)", s2, fixed = TRUE))

# #18/#21
ck("#18 调用 pfind_postprocess_fresh + 标记文件",
   grepl(".pfind_postprocess_done", s2, fixed = TRUE) &&
   grepl("pfind_postprocess_fresh(pp_dir, fasta_path, pp_mark_file)", s2, fixed = TRUE) &&
   grepl("pfind_postprocess_mark(pp_dir, fasta_path, pp_mark_file)", s2, fixed = TRUE))
ck("#21 失败写 .pfind_postprocess_failed 并醒目标记",
   grepl(".pfind_postprocess_failed", s2, fixed = TRUE) &&
   grepl("口径降级", s2, fixed = TRUE))

# #19
ck("#19 读段旋钮：pfind_reader_workers 解析",
   grepl('cfg_get_preferred(base_cfg, cfg2, "pfind_reader_workers", "")$value', s2, fixed = TRUE))
ck("#19 后处理传 reader_parallel/reader_workers",
   grepl("reader_parallel = rdr_par, reader_workers = rdr_n", s2, fixed = TRUE))
ck("#19 reader 不再硬编码 parallel=TRUE",
   grepl("parallel = reader_parallel", rd, fixed = TRUE) &&
   !grepl("parallel = TRUE)", rd, fixed = TRUE))
ck("#19 PSOCK worker 逐文件进度",
   grepl("[pFind_res_reader] worker 开始", rd, fixed = TRUE) &&
   grepl("[pFind_res_reader] worker 完成", rd, fixed = TRUE))
ck("#19(v1.0) PSOCK outfile=\"\" 使 worker 输出可见（工作站验收 P2）",
   grepl("makePSOCKcluster(nw, outfile = \"\")", rd, fixed = TRUE))

# 模板与版本
tpl <- paste(readLines(file.path(cd, "..", "config_step1_template.txt"), warn = FALSE),
             collapse = "\n")
ck("模板含 pfind_reader_workers 键", grepl("pfind_reader_workers\t0", tpl, fixed = TRUE))
ver <- readLines(file.path(cd, "..", "VERSION"), warn = FALSE)
ck("VERSION = v1.0.NN", grepl("v2\\.0\\.[0-9]{2}", ver))

cat("==== workstation_v262_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
