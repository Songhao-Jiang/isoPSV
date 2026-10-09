# python_exe_test.R — python_exe 显式键与 py_detect 回退单元测试
# （v1.0，2026-09-27：用户 Python 未写入环境目录 PATH 时，config_step1
# 的 python_exe 指向解释器完整路径即优先使用；路径不存在报消息并回退 PATH）
# Unit tests for the python_exe explicit key and py_detect fallback
# (v1.0, 2026-09-27: when Python is NOT on PATH, the config_step1
# python_exe key points at the full interpreter path and wins; a missing
# path is reported and PATH probing is the fallback).
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R"))) file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

# ---- 1. hint 指向存在的文件：优先返回该路径，pre 为空 ----
tmp_exe <- tempfile(fileext = ".exe")
file.create(tmp_exe)
p1 <- py_detect(tmp_exe)
ck("python_exe 存在 → cmd 用显式路径", is.list(p1) && identical(p1$cmd, tmp_exe))
ck("python_exe 存在 → pre 为空", is.list(p1) && identical(p1$pre, character(0)))
unlink(tmp_exe)

# ---- 2. hint 指向不存在的文件：报消息并回退 PATH 探测（不回退成坏路径）----
bad <- file.path(tempdir(), "no_such_python_9x7.exe")
msg <- capture.output(p2 <- py_detect(bad), type = "message")
ck("python_exe 不存在 → 不得返回坏路径",
   is.null(p2) || (!identical(p2$cmd, bad)))
ck("python_exe 不存在 → 有回退提示消息", length(msg) > 0L &&
   any(grepl("回退", msg)) || (is.null(p2) && length(msg) == 0L && TRUE))

# ---- 3. 无 hint 与空串 hint：与旧行为一致（PATH 探测或 NULL）----
p3 <- py_detect()
p4 <- py_detect("")
if (!is.null(p3)) {
  ck("无 hint → PATH 探测结果合法", is.list(p3) &&
     p3$cmd %in% c("python", "py", "python3"))
} else {
  ck("无 hint → 本机无 Python（NULL 合法）", is.null(p3))
}
ck("空串 hint → 等同无 hint", identical(p3, p4))

# ---- 4. NA / NULL hint 不报错 ----
p5 <- suppressWarnings(py_detect(NA_character_))
ck("NA hint → 不报错且结构合法", is.null(p5) || is.list(p5))
p6 <- py_detect(NULL)
ck("NULL hint → 等同无 hint", identical(p3, p6))

cat("汇总：", npass, " PASS / ", nfail, " FAIL\n", sep = "")
if (nfail > 0L) quit(status = 1L)
