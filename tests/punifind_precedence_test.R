# punifind_precedence_test.R — pUniFind 配置键优先级回归（v1.0 功能排队②）
# 覆盖：step1 优先、step2 兜底、无效值跳过（validate）、双空取 default、
# 来源标记正确、workers 键环境变量兜底语义。
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

s1 <- list(punifind_exe = "C:/pf/step1", punifind_rescore = "1",
           punifind_num_workers = "6")
s2 <- list(punifind_exe = "C:/pf/step2", punifind_rescore = "2",
           punifind_denovo = "1", punifind_python = "C:/py")
# 1) step1 优先
r1 <- cfg_get_preferred(s1, s2, "punifind_exe", "")
ck("exe 键 step1 优先", r1$value == "C:/pf/step1" && r1$source == "config_step1")
ck("rescore 模式 step1 优先", cfg_get_preferred(s1, s2, "punifind_rescore", "0")$value == "1")
# 2) step2 兜底
r2 <- cfg_get_preferred(s1, s2, "punifind_denovo", "0")
ck("step1 缺键时 step2 兜底", r2$value == "1" && r2$source == "config_step2")
r3 <- cfg_get_preferred(s1, s2, "punifind_python", "")
ck("python 键 step2 兜底（此前先命中即 break 会取 step2，现仍可取但来源明确）",
   r3$value == "C:/py" && r3$source == "config_step2")
# 3) 无效值跳过（与生产 user_step2 同款模式校验）
s3 <- list(punifind_rescore = "9")   # 非法模式码
mode_valid <- function(m) as.character(m) %in% c("0", "1", "2")
ck("step1 无效值被跳过、step2 有效值生效（rescore）",
   identical(cfg_get_preferred(s3, s2, "punifind_rescore", "0", mode_valid),
             list(value = "2", source = "config_step2")))
ck("step1 无效值时 step2 有效值生效（denovo）",
   cfg_get_preferred(s3, s2, "punifind_denovo", "0", mode_valid)$value == "1")
ck("双端都无效时落 default",
   identical(cfg_get_preferred(s3, list(punifind_rescore = "7"),
                               "punifind_rescore", "0", mode_valid),
             list(value = "0", source = "")))
# 4) 双空取 default
ck("双空取 default（值+来源空）",
   identical(cfg_get_preferred(list(), list(), "punifind_exe", "DEF"),
             list(value = "DEF", source = "")))
# 5) workers 键环境变量兜底语义（default=环境值 → 来源空、值=环境值）
ck("workers 键 default（环境变量）语义",
   identical(cfg_get_preferred(list(), list(), "punifind_num_workers", "8"),
             list(value = "8", source = "")))
ck("workers 键 step1 优先于环境 default",
   cfg_get_preferred(s1, list(), "punifind_num_workers", "8")$value == "6")

cat(sprintf("=== punifind_precedence_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
