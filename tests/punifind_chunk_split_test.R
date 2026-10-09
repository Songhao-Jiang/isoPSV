# punifind_chunk_split_test.R — v1.0/70 pUniFind 单 MGF 分片回归（2026-09-29）
# 覆盖（mock runner 注入，不跑真实推理）：多 MGF 逐片调用（每片恰 1 个 mgf＋
# 全量 qry.res 副本）、denovo CSV 无表头拼接、rescore 谱表表头去重合并、
# 片级完成标记与失败续跑（失败片重跑、成功片跳过）、单 MGF 直通原路径；
# v1.0：qry.res 片内改写不穿透真件（P2）、输入指纹失效与旧标记重跑、
# rescore 表头不一致即停（P3）、denovo 三参数门控三态（P1）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "punifind_gpu_compat.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

mkproj <- function(n = 3L) {
  td <- tempfile("pfk_"); dir.create(td)
  mgfs <- file.path(td, "mgf"); qry <- file.path(td, "qry")
  dir.create(mgfs); dir.create(qry)
  for (i in seq_len(n)) writeLines(c("BEGIN IONS", paste0("TITLE=r", i, ".1.1.2.0.dta"),
                                     "100.1 5", "END IONS"),
                                   file.path(mgfs, paste0("raw", i, "_HCDFT.mgf")))
  # v1.0（2026-10-01 P1）：片内 qry.res 按本片 raw 内容级过滤——夹具须带
  # per-raw 块（S 行＋raw 前缀 dta 行）才能体现"拷贝后过滤"语义：p1 含 raw1+raw2
  # 块、p2 含 raw3 块 → 每片恰 1 件（自身块的文件），全量拷贝时代是 2 件。
  # v1.0 (P1 2026-10-01): chunk qry.res is content-filtered per raw — the
  # fixture needs per-raw blocks ("S" line + raw-prefixed dta line): p1 holds
  # raw1+raw2, p2 holds raw3, so each chunk keeps exactly the 1 file carrying
  # its own block (the old full-copy behavior saw 2 files per chunk).
  # 块头须 "S	"（制表符）——过滤器按 startsWith(ln, "S	") 识别块首
  # block headers must be "S	" (tab) — the filter detects blocks via
  # startsWith(ln, "S	")
  q1 <- c("S	1	2", "raw1.1.1.2.0.dta	1	2", "S	3	4", "raw2.1.1.2.0.dta	3	4")
  q2 <- c("S	5	6", "raw3.1.1.2.0.dta	5	6")
  writeLines(q1, file.path(qry, "p1.qry.res"))
  writeLines(q2, file.path(qry, "p2.qry.res"))
  list(td = td, mgfs = mgfs, qry = qry)
}
mock <- function(log) function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
  cat(paste0(kind, "|", length(mgfs_in(mgf_dir)), "|", basename(out_dir), "|",
             length(list.files(qry_res_dir, pattern = "\\.qry\\.res$")), "\n"),
      file = log, append = TRUE)
  if (identical(kind, "denovo")) {
    dir.create(file.path(out_dir, "pUniFind_result"), recursive = TRUE)
    writeLines(paste0(basename(mgf_dir), ",row"),
               file.path(out_dir, "pUniFind_result", "test_5_0.csv"))
  } else {
    writeLines(c("File_Name\tScan_No", paste0(basename(mgf_dir), "\t1")),
               file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra"))
  }
  NULL
}
mgfs_in <- function(d) list.files(d, pattern = "\\.mgf$")

# 1) denovo：三片调用＋无表头拼接
p <- mkproj(3L)
logf <- tempfile(); invisible(file.create(logf))
r <- pf_run("denovo", "PF", p$qry, p$mgfs, p$td, .runner = mock(logf))
m <- readLines(r[["denovo_csv"]], warn = FALSE)
fld <- do.call(rbind, strsplit(readLines(logf), "|", fixed = TRUE))
# v1.0 起断言改为片内过滤口径：每片 qry 目录恰 1 件（带本片块的文件；
# p1 在 raw1/raw2 片各保留 1 件，p2 仅 raw3 片保留）
# since v1.0 assert the filtered convention: each chunk sees exactly 1
# qry.res (the file carrying its own block; p1 survives in raw1/raw2 chunks,
# p2 only in the raw3 chunk)
ck("denovo 三片各调用一次（每片恰 1 mgf、qry 片内过滤后各 1 件）",
   nrow(fld) == 3L && all(fld[, 1] == "denovo") &&
   all(fld[, 2] == "1") && all(fld[, 4] == "1") &&
   setequal(fld[, 3], c("raw1", "raw2", "raw3")))
ck("denovo 合并＝3 行无表头拼接", length(m) == 3L &&
   identical(m, c("mgf,row", "mgf,row", "mgf,row")))
ck("合并产物落回原接口位置", basename(r[["denovo_csv"]]) == "test_5_merged.csv" &&
   grepl("pUniFind_result$", dirname(r[["denovo_csv"]])))

# 2) rescore：表头只留一次
logf2 <- tempfile(); invisible(file.create(logf2))
r2 <- pf_run("rescore", "PF", p$qry, p$mgfs, p$td, .runner = mock(logf2))
m2 <- readLines(r2[["rescore_spectra"]], warn = FALSE)
ck("rescore 合并去重表头（1 表头＋3 行）", sum(grepl("^File_Name\t", m2)) == 1L &&
   length(m2) == 4L)
ck("rescore 合并落回真实 qry 目录",
   normalizePath(dirname(r2[["rescore_spectra"]]), winslash = "/") ==
   normalizePath(p$qry, winslash = "/"))

# 3) 完成标记：重跑全跳过，合并产物不变
logf3 <- tempfile(); invisible(file.create(logf3))
r3 <- pf_run("denovo", "PF", p$qry, p$mgfs, p$td, .runner = mock(logf3))
ck("重跑零调用（标记跳过）", length(readLines(logf3)) == 0L &&
   identical(readLines(r3[["denovo_csv"]], warn = FALSE), m))

# 4) 失败续跑：chunk2 失败 → 报 chunk 身份；重跑只补失败片
p4 <- mkproj(3L)
logf4 <- tempfile(); invisible(file.create(logf4))
fail2 <- function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
  if (identical(basename(out_dir), "raw2")) stop("boom")
  mock(logf4)(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...)
}
err <- tryCatch(pf_run("denovo", "PF", p4$qry, p4$mgfs, p4$td, .runner = fail2),
                error = function(e) e)
ck("失败报 chunk 身份且带续跑指引", inherits(err, "error") &&
   grepl("raw2", conditionMessage(err)) && grepl("续跑", conditionMessage(err)))
ck("失败后成功片有标记（1/3）",
   length(list.files(file.path(p4$td, "chunks"), pattern = "^done_denovo",
                     recursive = TRUE)) == 1L)
logf5 <- tempfile(); invisible(file.create(logf5))
r5 <- pf_run("denovo", "PF", p4$qry, p4$mgfs, p4$td, .runner = mock(logf5))
ck("续跑只补失败片（2 次新调用）", length(readLines(logf5)) == 2L &&
   length(readLines(r5[["denovo_csv"]], warn = FALSE)) == 3L)

# 5) 单 MGF：直通原路径（runner 收到原目录本体）
p6 <- mkproj(1L)
logf6 <- tempfile(); invisible(file.create(logf6))
pf_run("denovo", "PF", p6$qry, p6$mgfs, p6$td, .runner = mock(logf6))
ck("单 MGF 直通（原 mgf/qry 目录、不建 chunks）",
   length(readLines(logf6)) == 1L &&
   !dir.exists(file.path(p6$td, "chunks")))

# 6) v1.0 P2：片内就地改写 qry.res 不得穿透项目真件（copy 而非硬链接）
p7 <- mkproj(2L)
qry1 <- file.path(p7$qry, "p1.qry.res")
before <- readLines(qry1, warn = FALSE)
rewrite_mock <- function(log) function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
  q <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)[1]
  writeLines(c(readLines(q, warn = FALSE), "S 9 9"), q)   # 模拟消毒/压实就地覆写
  mock(log)(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...)
}
logf7 <- tempfile(); invisible(file.create(logf7))
pf_run("denovo", "PF", p7$qry, p7$mgfs, p7$td, .runner = rewrite_mock(logf7))
ck("片内改写 qry.res 不穿透真件（内容不变）",
   identical(readLines(qry1, warn = FALSE), before))

# 7) v1.0：done 标记带输入指纹——输入变化即作废重跑、未变即跳过
logf8 <- tempfile(); invisible(file.create(logf8))
pf_run("denovo", "PF", p7$qry, p7$mgfs, p7$td, .runner = mock(logf8))
ck("输入未变重跑零调用", length(readLines(logf8)) == 0L)
writeLines(c(readLines(qry1, warn = FALSE), "S 8 8"), qry1)  # 换库重搜模拟
logf9 <- tempfile(); invisible(file.create(logf9))
r9 <- pf_run("denovo", "PF", p7$qry, p7$mgfs, p7$td, .runner = mock(logf9))
ck("输入变化后全片重跑（指纹作废）", length(readLines(logf9)) == 2L &&
   length(readLines(r9[["denovo_csv"]], warn = FALSE)) == 2L)
mk1 <- file.path(p7$td, "chunks", "raw1", "done_denovo.txt")
writeLines(c("chunk: raw1", paste0("at: ", Sys.time())), mk1)  # v1.0 无指纹旧标记
logf10 <- tempfile(); invisible(file.create(logf10))
pf_run("denovo", "PF", p7$qry, p7$mgfs, p7$td, .runner = mock(logf10))
ck("旧版无指纹标记 → 该片重跑、其余片照常跳过",
   length(readLines(logf10)) == 1L &&
   grepl("^fingerprint: ", readLines(mk1, warn = FALSE)[3]))

# 8) v1.0 P3：rescore 各片表头不一致 → 合并即停
p8 <- mkproj(2L)
logf11 <- tempfile(); invisible(file.create(logf11))
hdr_mock <- function(log) function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
  if (identical(basename(out_dir), "raw2"))
    writeLines(c("File_Name\tScan_No\tEXTRA", paste0(basename(mgf_dir), "\t1\tX")),
               file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra"))
  else
    writeLines(c("File_Name\tScan_No", paste0(basename(mgf_dir), "\t1")),
               file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra"))
  cat(paste0(kind, "|", basename(out_dir), "\n"), file = log, append = TRUE)
  NULL
}
err8 <- tryCatch(pf_run("rescore", "PF", p8$qry, p8$mgfs, p8$td, .runner = hdr_mock(logf11)),
                 error = function(e) e)
ck("rescore 表头不一致即停（报片名）", inherits(err8, "error") &&
   grepl("表头不一致", conditionMessage(err8)))

# 9) v1.0 P1：denovo 三参数门控三态（flags=NULL 时的定向探针结果）
d1 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL, denovo_len_ok = TRUE)
ck("flags=NULL＋探针TRUE → 三参带上",
   all(c("--min-length=6", "--max-length=40") %in% d1))
d2 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL, denovo_len_ok = FALSE)
ck("flags=NULL＋探针FALSE → 三参不带", !any(grepl("min-length|max-length", d2)))
d3 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL)
ck("flags=NULL＋不可判 → 按官方口径带上",
   all(c("--min-length=6", "--max-length=40") %in% d3))
d4 <- pf_build_args("denovo", "P", "Q", "M", "O",
                    flags = c("--user-dir", "--task"))
ck("flags 探测成功且无三参 → 不带", !any(grepl("min-length|max-length", d4)))

# 10) denovo 合并行数守恒（merged ＝ Σ 各片 csv 行数；p7 为 denovo 全跑项目）
csv_sum <- 0L
for (st in c("raw1", "raw2")) {
  fs <- list.files(file.path(p7$td, "chunks", st, "pUniFind_result"),
                   pattern = "\\.csv$", full.names = TRUE)
  for (f in fs) csv_sum <- csv_sum + length(readLines(f, warn = FALSE))
}
mg <- file.path(p7$td, "pUniFind_result", "test_5_merged.csv")
ck("denovo 合并行数守恒", file.exists(mg) &&
   length(readLines(mg, warn = FALSE)) == csv_sum)

# 11) 指纹纯函数：mtime|size 变化可检
tdf <- tempfile("pff_"); dir.create(tdf)
fq <- file.path(tdf, "a.qry.res"); fm <- file.path(tdf, "m.mgf")
writeLines("S 1 2", fq); writeLines("BEGIN IONS", fm)
f1 <- pf_chunk_fingerprint(fq, fm)
Sys.sleep(1.1); writeLines(c("S 1 2", "S 3 4"), fq)
f2 <- pf_chunk_fingerprint(fq, fm)
ck("指纹随输入 mtime/size 变化", !identical(f1, f2) &&
   identical(pf_chunk_fingerprint(fq, fm), f2))
unlink(tdf, recursive = TRUE)


# 12) v1.0 A①：batch-size 参数化
b1 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL, denovo_len_ok = TRUE,
                    batch_size = 256)
ck("batch_size=256 → --batch-size=256", any(grepl("^--batch-size=256$", b1)))
b2 <- pf_build_args("rescore", "P", "Q", "M", "O", flags = NULL)
ck("默认 batch_size=64 保持", any(grepl("^--batch-size=64$", b2)))

# 13) v1.0 A③：片级并行（chunk_parallel=2，mock 带 0.2s 睡眠）
p13 <- mkproj(4L)
logf13 <- tempfile(); invisible(file.create(logf13))
self_mock <- function(log, sleep = 0) {
  force(log); force(sleep)
  function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
    # 自包含：并行 worker 上无法引用测试脚本的其他全局函数；force 后形参
    # 随闭包序列化（clusterExport 不落地未强制 promise，2026-09-30 实测）。
    # 计数写片内 ran.txt——多 worker 并发 append 同一日志会互相覆盖（实测）
    if (sleep > 0) Sys.sleep(sleep)
    writeLines("ran", file.path(out_dir, "ran.txt"))
    if (identical(kind, "denovo")) {
      dir.create(file.path(out_dir, "pUniFind_result"), recursive = TRUE)
      writeLines("mgf,row", file.path(out_dir, "pUniFind_result", "test_5_0.csv"))
    } else {
      writeLines(c("File_Name\tScan_No", paste0(basename(mgf_dir), "\t1")),
                 file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra"))
    }
    NULL
  }
}
sleep_mock <- function(log) self_mock(log, sleep = 0.2)
r13 <- pf_run("denovo", "PF", p13$qry, p13$mgfs, p13$td, .runner = sleep_mock(logf13),
              chunk_parallel = 2L)
ck("并行 4 片全跑＋合并守恒",
   length(list.files(file.path(p13$td, "chunks"), pattern = "^ran[.]txt$",
                     recursive = TRUE)) == 4L &&
   length(readLines(r13[["denovo_csv"]], warn = FALSE)) == 4L &&
   length(list.files(file.path(p13$td, "chunks"), pattern = "^done_denovo",
                     recursive = TRUE)) == 4L)

# 14) 并行失败汇总：raw2 失败 → 报全部失败片、成功片标记保留
p14 <- mkproj(3L)
logf14 <- tempfile(); invisible(file.create(logf14))
fail_raw2 <- local({
  log <- logf14; force(log)
  function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...) {
    if (identical(basename(out_dir), "raw2")) {
      writeLines(c("Traceback REAL_ROOT_CAUSE_LINE in child", "OOM killer"),
                 file.path(out_dir, "denovo.log"))
      stop("boom")
    }
    # 自包含产出（并行 worker 上无法引用测试脚本全局；log 经 local+force 烘焙）
    cat(paste0(basename(out_dir), "\n"), file = log, append = TRUE)
    if (identical(kind, "denovo")) {
      dir.create(file.path(out_dir, "pUniFind_result"), recursive = TRUE)
      writeLines("mgf,row", file.path(out_dir, "pUniFind_result", "test_5_0.csv"))
    } else {
      writeLines(c("File_Name\tScan_No", paste0(basename(mgf_dir), "\t1")),
                 file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra"))
    }
    NULL
  }
})
err14 <- tryCatch(pf_run("denovo", "PF", p14$qry, p14$mgfs, p14$td, .runner = fail_raw2,
                         chunk_parallel = 2L), error = function(e) e)
ck("并行失败：报片名＋子日志尾透出", inherits(err14, "error") &&
   grepl("raw2", conditionMessage(err14)) &&
   grepl("REAL_ROOT_CAUSE_LINE", conditionMessage(err14)))
ck("并行失败：成功片标记保留（2/3）",
   length(list.files(file.path(p14$td, "chunks"), pattern = "^done_denovo",
                     recursive = TRUE)) == 2L)

# 15) 顺序路径同样附子日志尾（B 项）
err15 <- tryCatch(pf_run("denovo", "PF", p14$qry, p14$mgfs, p14$td, .runner = fail_raw2),
                  error = function(e) e)
ck("顺序失败也带子日志尾", inherits(err15, "error") &&
   grepl("子日志尾", conditionMessage(err15)) &&
   grepl("REAL_ROOT_CAUSE_LINE", conditionMessage(err15)))

# 16) v1.0 D：chunk_cleanup=1 → 片内数据清、标记与合并产物留
logf16 <- tempfile(); invisible(file.create(logf16))
r16 <- pf_run("denovo", "PF", p14$qry, p14$mgfs, p14$td, .runner = mock(logf16),
              chunk_cleanup = TRUE)
left16 <- list.files(file.path(p14$td, "chunks"), recursive = TRUE)
data_dirs16 <- list.files(file.path(p14$td, "chunks"),
                          pattern = "^(mgf|qry|pUniFind_data_tmp|pUniFind_result)$",
                          recursive = FALSE)
ck("清理后数据目录全清、done 标记与合并产物保留",
   length(data_dirs16) == 0L &&
   sum(grepl("^done_denovo[.]txt$", basename(left16))) >= 1L &&
   file.exists(r16[["denovo_csv"]]))

# 17) v1.0 C：探针日志双编码读取（GBK 字节行不再炸）
td17 <- tempfile("pfenc_"); dir.create(td17)
f17 <- file.path(td17, "probe.log")
# "注释" 的 GBK 字节 D7 A2 CA CD ＋ 一行 ASCII flag
writeBin(as.raw(c(0xD7, 0xA2, 0xCA, 0xCD, 0x0A, 0x2D, 0x2D, 0x6D, 0x69,
                 0x6E, 0x2D, 0x6C, 0x65, 0x6E, 0x67, 0x74, 0x68, 0x0A)), f17)
ln17 <- pf_read_log_lines(f17)
ck("GBK 日志双编码读取成功", !is.null(ln17) && length(ln17) == 2L &&
   grepl("min-length", ln17[2]))
unlink(td17, recursive = TRUE)

cat("==== punifind_chunk_split_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
