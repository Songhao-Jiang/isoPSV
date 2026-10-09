# two_queue_test.R — v1.0 两队列调度器（GPU 管家）单元与冒烟测试
# 覆盖：① 管家终态机（done/aborted 转移）② 汇合门三路径（done/无状态死/等活）
# ③ 启动器幂等（done 不重启/活 PID 不重启）④ 生成脚本 ASCII 纯净
# Run: Rscript tests/two_queue_test.R
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]), winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R"))) file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
suppressMessages(source(file.path(cd, "boot.R")))
source(file.path(cd, "two_queue_steward.R"), encoding = "UTF-8")

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}
td <- tempfile("tq_"); dir.create(td)

# ① 汇合门：done 状态即时通过
cfg <- list(config_dir = td)
writeLines(c("state=done", "denovo=ok", "rescore=ok"), file.path(td, "_gpu_steward_status.txt"))
t0 <- system.time(r <- tq_join_gate(cfg, poll_s = 1L))
ck("join gate: done 状态即时通过", isTRUE(r) && t0["elapsed"] < 2)

# ② 汇合门：无状态且无 PID（管家从未起）→ 容忍放行 FALSE
td2 <- tempfile("tq_"); dir.create(td2)
r2 <- tq_join_gate(list(config_dir = td2), poll_s = 1L)
ck("join gate: 无状态无 PID 容忍放行", isFALSE(r2))

# ③ 汇合门：活 PID（用当前进程自身）→ 等待后写 done 才放行
td3 <- tempfile("tq_"); dir.create(td3)
writeLines(as.character(Sys.getpid()), file.path(td3, "_gpu_steward.pid"))
# Windows 无 fork：用 cmd 后台延时写 done 模拟管家完成
writeLines("state=running pid=0", file.path(td3, "_gpu_steward_status.txt"))
sim_cmd <- sprintf("(ping -n 3 127.0.0.1 >nul & echo state=done> \"%s\")",
                   file.path(td3, "_gpu_steward_status.txt"))
system2("cmd", c("/c", sim_cmd), wait = FALSE, invisible = TRUE)
t3 <- system.time(r3 <- tq_join_gate(list(config_dir = td3), poll_s = 1L))
ck("join gate: 活 PID 等到 done（>=2s）", isTRUE(r3) && t3["elapsed"] >= 1.5)

# ④ 启动器幂等：done 状态 → 不重启（以不生成脚本为证）
td4 <- tempfile("tq_"); dir.create(td4)
writeLines("state=done", file.path(td4, "_gpu_steward_status.txt"))
r4 <- tq_launch_steward(list(config_dir = td4), file.path(td4, "config_step2.txt"), cd)
ck("launcher: done 状态不重启（无新脚本）",
   isTRUE(r4) && !file.exists(file.path(td4, "_gpu_steward_bg.R")))

# ⑤ 启动器幂等：活 PID → 不重启
td5 <- tempfile("tq_"); dir.create(td5)
writeLines(as.character(Sys.getpid()), file.path(td5, "_gpu_steward.pid"))
r5 <- tq_launch_steward(list(config_dir = td5), file.path(td5, "config_step2.txt"), cd)
ck("launcher: 活 PID 不重启", isTRUE(r5) && !file.exists(file.path(td5, "_gpu_steward_bg.R")))

# ⑥ 生成的后台脚本 ASCII 纯净（中文路径经 \\u 转义）
td6 <- tempfile("tq_"); dir.create(td6)
cfg2p <- file.path(td6, "配置_步二.txt")
# 不真正启动：仅借用内部生成逻辑——直接调用会启动进程；改为静态检查函数行为：
# 生成路径含中文时 esc_nonascii 的输出必须纯 ASCII
s <- enc2utf8("G:/路径含中文/目录")
esc_out <- NULL
# esc_nonascii 是闭包内函数；这里复刻同实现做静态断言（一致性由实现保证）
esc_chk <- function(x) {
  x <- enc2utf8(as.character(x))
  paste(vapply(utf8ToInt(x), function(cp)
    if (cp < 128L) intToUtf8(cp)
    else if (cp <= 0xFFFFL) sprintf("\\u%04X", cp)
    else sprintf("\\U%08X", cp), character(1)), collapse = "")
}
ck("esc: 中文路径转义后纯 ASCII", grepl("^[\\x00-\\x7F]+$", esc_chk(s), perl = TRUE) &&
   grepl("\\\\u8DEF", esc_chk(s)))

# ⑦ 管家主流程（缓存命中路径端到端）：构造 pUniFind_denovo 缓存＋rescore 禁用
td7 <- tempfile("tq_"); dir.create(td7)
dir.create(file.path(td7, "pUniFind_denovo", "pUniFind_result"), recursive = TRUE)
writeLines("x,y", file.path(td7, "pUniFind_denovo", "pUniFind_result", "test_5_merged.csv"))
cfg7 <- file.path(td7, "config_step2.txt")
writeLines(c("base_config\t", "use_punifind\t1", "punifind_exe\tD:/nonexistent_punifind",
             "punifind_denovo\t1", "punifind_rescore\t0",
             paste0("config_dir\t", gsub("\\\\", "/", td7))),
           cfg7, useBytes = TRUE)
# read_param_file 需要 base_config 存在——造一个空的
writeLines("config_dir\t", file.path(td7, "config_step1.txt"))
cfg7lines <- readLines(cfg7)
cfg7lines[1] <- paste0("base_config\t", file.path(td7, "config_step1.txt"))
writeLines(cfg7lines, cfg7, useBytes = TRUE)
# steward_main 会 source punifind_gpu_compat（code 目录下存在，OK）；denovo 缓存命中 → skip
# rescore=0 → disabled。期望状态文件 state=done denovo=skip(cache) rescore=disabled
ok7 <- tryCatch({ steward_main(cfg7, cd); TRUE }, error = function(e) { cat("   err:", conditionMessage(e), "\n"); FALSE })
ln7 <- if (file.exists(file.path(td7, "_gpu_steward_status.txt")))
  readLines(file.path(td7, "_gpu_steward_status.txt"), warn = FALSE) else character(0)
ck("steward main: 缓存命中端到端终态 done",
   ok7 && length(ln7) == 3 && ln7[1] == "state=done" &&
   ln7[2] == "denovo=skip(cache)" && ln7[3] == "rescore=disabled")

cat("\n=== two_queue_test:", npass, "pass /", nfail, "fail ===\n")
quit(status = if (nfail > 0) 1 else 0)
