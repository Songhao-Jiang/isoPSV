# punifind_compat_test.R — pUniFind GPU 适配模块回归测试（v1.0，2026-09-26）
# 适配层合成断言（无需 pUniFind 安装）：环境开关、argv 规则（数据目录第一、
# 无 --load-contra、num-workers 生效、flag 过滤）、qry.res 序号压实幂等。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "punifind_gpu_compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "punifind_gpu_compat.R"), encoding = "UTF-8")

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# 1) 环境变量：两个 torch 官方开关 + num_workers
ev <- pf_env_vars(8)
evkv <- paste0(names(ev), "=", ev)
ck("子进程环境含 USE_LIBUV=0", "USE_LIBUV=0" %in% evkv)
ck("子进程环境含 TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1",
   "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1" %in% evkv)
ck("子进程环境含 PUNIFIND_NUM_WORKERS=8", "PUNIFIND_NUM_WORKERS=8" %in% evkv)
ev2 <- pf_env_vars(distributed = TRUE)
ck("python 直调路线含单进程分布式变量",
   all(c("RANK=0", "WORLD_SIZE=1", "MASTER_ADDR=127.0.0.1") %in% paste0(names(ev2), "=", ev2)))

# 2) argv 规则
a <- pf_build_args("rescore", "P", "Q", "M", "O", flags = NULL, num_workers = 8)
ck("argv[1] 是数据目录（位置参数前置）", substr(a[1], 1, 1) != "-")
ck("不含已消失的 --load-contra", !any(grepl("load-contra", a)))
ck("num-workers 已注入 --num-workers=8", any(a == "--num-workers=8"))
ck("数值参数 --flag=value 连写（--batch-size=64）", any(a == "--batch-size=64"))
a2 <- pf_build_args("rescore", "P", "Q", "M", "O",
                    flags = c("--user-dir", "--task", "--num-proc"), num_workers = 4)
ck("只保留探测到的 flag（过滤生效）",
   any(grepl("^--task=", a2)) && !any(grepl("^--batch-size", a2)) && any(a2 == "--num-proc=8"))
a3 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL)
ck("denovo task=denovo_inference", any(a3 == "--task=denovo_inference"))

# 3) qry.res 序号压实：合成夹具 + 幂等 + 内容不变式
td <- tempfile("pfst_"); dir.create(td)
f <- file.path(td, "x.qry.res")
writeLines(c("S 1234.5 2", "spec_a.dta", "1 PEP K 2 0.9 0.01 0",
             "3 PEPR 2 0.8 0.02 0", "S 2345.6 2", "spec_b.dta", "1 PEPTIDE 2 0.7 0.03 0"), f)
r1 <- pf_qryres_renumber(td, verbose = FALSE)
r2 <- pf_qryres_renumber(td, verbose = FALSE)
ln <- readLines(f, warn = FALSE)
ck("序号压实生效且幂等（首次 1 行 / 二次 0 行）", r1$fixed[1] == 1 && r2$fixed[1] == 0)
ck("候选序号已连续 1..n", identical(sub("^([0-9]+) .*", "\\1", ln[3:4]), c("1", "2")))
ck("压实不改候选内容（除行首序号外字节原样）",
   identical(sub("^[0-9]+", "", ln[4]), sub("^[0-9]+", "", "3 PEPR 2 0.8 0.02 0")))
unlink(td, recursive = TRUE)

# 4) N-21（2026-09-26，GLM）：pf_run 产物发现此前零覆盖——抽出的三个纯函数
#    （N-4 回查正则 / P0-B 新鲜度 / P0-E 缓存失效）用临时目录夹具直接测
td2 <- tempfile("pfrun_"); dir.create(td2)
qry_dir <- file.path(td2, "qry"); res_dir2 <- file.path(td2, "res"); dtm <- file.path(td2, "dt")
dir.create(qry_dir); dir.create(res_dir2); dir.create(dtm)
# N-4：fdr0.01_（新构建自写）与 fdr.01_（旧写）两写都认
writeLines("x", file.path(qry_dir, "testfdr0.01_pUniFind.spectra"))
writeLines("x", file.path(qry_dir, "oldfdr.01_pUniFind.spectra"))
ck("N-4 回查两写都认（fdr0.01_/fdr.01_）", length(pf_match_final_spectra(qry_dir)) == 2L)
# P0-B：陈旧 .spectra 被新鲜度门槛拒绝
writeLines("x", file.path(res_dir2, "testfdr0.01_pUniFind.spectra"))
Sys.setFileTime(file.path(res_dir2, "testfdr0.01_pUniFind.spectra"), Sys.time() - 300)
ck("P0-B 陈旧产物被拒（早于运行开始）",
   length(pf_fresh_own_spectra(res_dir2, Sys.time())) == 0L)
Sys.setFileTime(file.path(res_dir2, "testfdr0.01_pUniFind.spectra"), Sys.time())
ck("P0-B 新产物被接受", length(pf_fresh_own_spectra(res_dir2, Sys.time())) == 1L)
# P0-E：qry.res 新于缓存 → 删除 parquet/lmdb/pkl；输入不新 → 保留
writeLines("x", file.path(qry_dir, "a.qry.res"))
writeLines("x", file.path(dtm, "test.parquet")); writeLines("x", file.path(dtm, "test.lmdb"))
writeLines("x", file.path(res_dir2, "pscore_res_1.pkl"))
# P0-2（2026-09-26，deepseek v1.0 审查）：推理检查点 pkl 真实命名为
# {subset}_{subset}_{rank}.pkl（isoPSV 侧 subset=test → test_test_0.pkl），
# 此前按 pscore_res_* 收窄会漏清使 P0-E 复现——本夹具把"推理检查点必被清"
# 锁成门禁；keys pkl 同族失效一并清（W1）
writeLines("x", file.path(res_dir2, "test_test_0.pkl"))
writeLines("x", file.path(dtm, "test_FDR0.1_keys.pkl"))
Sys.setFileTime(file.path(dtm, "test.parquet"), Sys.time() - 300)
Sys.setFileTime(file.path(dtm, "test.lmdb"), Sys.time() - 300)
Sys.setFileTime(file.path(res_dir2, "pscore_res_1.pkl"), Sys.time() - 300)
Sys.setFileTime(file.path(res_dir2, "test_test_0.pkl"), Sys.time() - 300)
Sys.setFileTime(file.path(dtm, "test_FDR0.1_keys.pkl"), Sys.time() - 300)
deleted_files <- pf_clear_stale_cache(qry_dir, dtm, res_dir2, verbose = FALSE)
ck("P0-E 输入更新清空三类缓存（parquet/lmdb/pkl）", length(deleted_files) == 5L)
ck("P0-2 推理检查点 test_test_0.pkl 必被清（W1）",
   file.path(res_dir2, "test_test_0.pkl") %in% deleted_files)
ck("P0-2 keys pkl 同族失效一并清",
   file.path(dtm, "test_FDR0.1_keys.pkl") %in% deleted_files)
ck("P0-E 缓存确实被删除",
   !any(file.exists(c(file.path(dtm, "test.parquet"), file.path(dtm, "test.lmdb"),
                      file.path(res_dir2, "pscore_res_1.pkl"),
                      file.path(res_dir2, "test_test_0.pkl"),
                      file.path(dtm, "test_FDR0.1_keys.pkl")))))
# 输入未变（缓存新于 qry.res）：不删
writeLines("x", file.path(dtm, "test.parquet")); writeLines("x", file.path(res_dir2, "pscore_res_1.pkl"))
writeLines("x", file.path(res_dir2, "test_test_0.pkl"))
Sys.setFileTime(file.path(qry_dir, "a.qry.res"), Sys.time() - 600)
kept_files <- pf_clear_stale_cache(qry_dir, dtm, res_dir2, verbose = FALSE)
ck("P0-E 输入未变缓存保留", length(kept_files) == 0L &&
   file.exists(file.path(dtm, "test.parquet")) && file.exists(file.path(res_dir2, "pscore_res_1.pkl")) &&
   file.exists(file.path(res_dir2, "test_test_0.pkl")))
# .spectra 产物不算输入（qry.res 目录里放更新的 .spectra 不触发清缓存）
writeLines("x", file.path(qry_dir, "b.spectra"))
extra_kept <- pf_clear_stale_cache(qry_dir, dtm, res_dir2, verbose = FALSE)
ck("P0-E .spectra 产物不计入输入新鲜度", length(extra_kept) == 0L)
unlink(td2, recursive = TRUE)

# 5) N-15（2026-09-26）：pf_help_flags timeout 生效 + 正常探测 + py 发行版探测
rscript <- file.path(R.home("bin"), "Rscript.exe")
# v1.0（2026-10-01）自动版本判定后的探测契约（原单断言"正常探测出 flag
# 集合"按旧契约编写，已过时）：
# ① 探测面缺关键 flag（--check/--fdr-thread/--cutoff-spectra，老构建特征）
#    → 自动降级返回 NULL（走全参数盲传）
# ② env=0 强制启用过滤 → 返回探测到的 flag 集（旧契约行为，显式请求才生效）
# ③ 探测面含全部关键 flag（新构建特征）→ 不降级，返回含自定义 flag 的全集
# v1.0 auto version detection changed the probe contract (the old single
# assertion predates it): ① help surface missing critical flags (old-build
# signature) auto-degrades to NULL (blind full-arg mode); ② env=0 forces the
# flag filter on and returns the probed set; ③ surface containing all
# critical flags (new build) returns the full set un-degraded.
fl_old <- pf_help_flags(rscript, c("--vanilla", "-e",
                                   shQuote("cat('--foo --bar-baz')")), timeout = 30)
ck("pf_help_flags 缺关键 flag 自动降级 NULL（v1.0）", is.null(fl_old))
Sys.setenv(PUNIFIND_NO_FLAGS_PROBE = "0")
fl_force <- pf_help_flags(rscript, c("--vanilla", "-e",
                                     shQuote("cat('--foo --bar-baz')")), timeout = 30)
Sys.unsetenv("PUNIFIND_NO_FLAGS_PROBE")
ck("pf_help_flags 强制启用（env=0）探测出 flag 集合",
   !is.null(fl_force) && all(c("--foo", "--bar-baz") %in% fl_force))
fl_new <- pf_help_flags(rscript, c("--vanilla", "-e",
      shQuote("cat('--foo --bar-baz --check --fdr-thread --cutoff-spectra')")),
      timeout = 30)
ck("pf_help_flags 含全部关键 flag 不降级（新构建路径）",
   !is.null(fl_new) && all(c("--foo", "--bar-baz", "--check",
                             "--fdr-thread", "--cutoff-spectra") %in% fl_new))
# 挂起超时：Rscript 睡 120s，timeout=5 → 5 秒左右返回 NULL（按 PID 杀进程树）
t_hang <- system.time(fl_hang <- pf_help_flags(rscript, c("--vanilla", "-e",
                                    shQuote("Sys.sleep(120)")), timeout = 5))
ck("pf_help_flags 挂起限时返回 NULL（N-15 timeout 生效）",
   is.null(fl_hang) && t_hang["elapsed"] < 30)
# py 发行版探测：无 exe、有 workflow.py + consts/ckpts → ok 且 py_route_ok
td5 <- tempfile("pfp_"); dir.create(td5)
dir.create(file.path(td5, "ckpts")); dir.create(file.path(td5, "consts"))
writeLines("x", file.path(td5, "ckpts", "checkpoint_rank.pt"))
writeLines("x", file.path(td5, "consts", "tokenize1.pkl"))
writeLines("x", file.path(td5, "consts", "modification_meta_dict.pkl"))
writeLines("print('wf')", file.path(td5, "official_score_workflow.py"))
writeLines("print('wf')", file.path(td5, "official_denovo_workflow.py"))
pr5 <- pf_probe(td5, python_exe = rscript)
ck("pf_probe 纯 .py 发行版（无 exe）通过且 py_route_ok（N-15）",
   isTRUE(pr5$ok) && isTRUE(pr5$py_route_ok))
pr5b <- pf_probe(td5, python_exe = NULL)
ck("pf_probe 无 python 时回退要求 exe（报缺）", !isTRUE(pr5b$ok))
unlink(td5, recursive = TRUE)
# 端口参数传递（N-15 MASTER_PORT 可换的基础）
ev5 <- pf_env_vars(8, distributed = TRUE, master_port = 41234)
ck("pf_env_vars 自定义 master_port 生效",
   "MASTER_PORT=41234" %in% paste0(names(ev5), "=", ev5))

cat(sprintf("=== punifind_compat_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
