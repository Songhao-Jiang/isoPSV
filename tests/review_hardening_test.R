# review_hardening_test.R — WorkBuddy 审核加固回归（2026-09-25）
# 覆盖：strategy_slots 全部 15 种引擎类型槽位映射（P3-4；original 为保留槽不计入）、read_param_file
# 边界行为（值含制表符不截断、重复键警告、无制表符行推断、BOM 剥离、空文件报错）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- file.path(dirname(sd), "code")
if (!dir.exists(cd)) cd <- file.path(sd, "..", "code")
source(file.path(cd, "compat.R"))
source(file.path(cd, "ReadMGF_script.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, ok) {
  cat(ifelse(ok, "PASS", "FAIL"), ":", label, "\n")
  if (ok) npass <<- npass + 1L else nfail <<- nfail + 1L
}

# ---- 1. strategy_slots 全映射（16 类型） ----
expect_map <- c(pfind_open = "b", pfind_restrict = "c", punifind_rescore = "d",
                punifind_denovo = "e", maxquant = "f", comet = "g",
                msfragger = "h", sage = "i", msfragger_boost = "j",
                xtandem = "k", metamorpheus = "l", peaks = "m", mascot = "n",
                pnovo = "o", msfragger_open = "p")
got <- strategy_slots(names(expect_map))
ck("strategy_slots 15 类型逐个映射",
   identical(unname(got), unname(expect_map)))
ck("strategy_slots other 通用槽 q 起",
   identical(strategy_slots(c("other", "other", "other")),
             c("q", "r", "s")))
ck("strategy_slots 未知类型报错",
   inherits(try(strategy_slots("no_such_type"), silent = TRUE), "try-error"))
ck("strategy_slots 重复类型报错",
   inherits(try(strategy_slots(c("comet", "comet")), silent = TRUE), "try-error"))
ck("strategy_slots original 禁止",
   inherits(try(strategy_slots("original"), silent = TRUE), "try-error"))
ck("STRATEGY_REVISE 全类型有定义（mascot/peaks/xtandem=none）",
   all(c("mascot", "peaks", "xtandem") %in% names(STRATEGY_REVISE)) &&
     all(STRATEGY_REVISE[c("mascot", "peaks", "xtandem")] == "none"))

# ---- 2. read_param_file 边界 ----
td <- tempfile("rh_"); dir.create(td)
p1 <- file.path(td, "cfg_tabvalue.txt")
writeLines(c("key1\tval_a\tval_b\tval_c", "key2\tplain"), p1, useBytes = TRUE)
c1 <- read_param_file(p1)
ck("值含制表符不截断（第 3 段起保留）", identical(c1$key1, "val_a\tval_b\tval_c"))

p2 <- file.path(td, "cfg_dup.txt")
writeLines(c("key1\tfirst", "key1\tsecond"), p2, useBytes = TRUE)
w2 <- character(0)
c2 <- withCallingHandlers(read_param_file(p2),
                          warning = function(w) { w2 <<- c(w2, conditionMessage(w)); invokeRestart("muffleWarning") })
ck("重复键全部保留 + 警告", identical(c2$key1, c("first", "second")) && length(w2) > 0)

p3 <- file.path(td, "cfg_notab.txt")
writeLines(c("key1  val by spaces", "key2\tok"), p3, useBytes = TRUE)
w3 <- character(0)
c3 <- withCallingHandlers(read_param_file(p3),
                          warning = function(w) { w3 <<- c(w3, conditionMessage(w)); invokeRestart("muffleWarning") })
ck("无制表符行按连续空格推断 + 警告",
   identical(c3$key1, "val by spaces") && length(w3) > 0)

p4 <- file.path(td, "cfg_bom.txt")
writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)),
           charToRaw("key1\tval1\n")), p4)
c4 <- read_param_file(p4)
ck("UTF-8 BOM 剥离", "key1" %in% names(c4) && identical(c4$key1, "val1"))

p5 <- file.path(td, "cfg_empty.txt")
writeLines(character(0), p5, useBytes = TRUE)
ck("空文件明确报错", inherits(try(read_param_file(p5), silent = TRUE), "try-error"))

p6 <- file.path(td, "cfg_comment.txt")
writeLines(c("# 注释行", "key1\tval1", "  # 缩进注释"), p6, useBytes = TRUE)
c6 <- read_param_file(p6)
ck("注释行跳过", identical(names(c6), c("key1", "config_dir")))

# N-4 配套（2026-09-25）：UTF-8 无 BOM 配置里的中文路径，读回后转回 UTF-8
# 必须与原文一致（Encoding 标记缺失时 GBK 代码页下会按本地字节串解释）
p7 <- file.path(td, "cfg_cn.txt")
cn_val <- "G:/数据/蛋白库.fasta"
Encoding(cn_val) <- "UTF-8"   # F-1（2026-09-26）：locale 非 UTF-8 时字面量按 native 标记，显式标 UTF-8 保证与读回值同标记可比
writeBin(charToRaw(paste0("fasta\t", cn_val, "\n")), p7)
c7 <- read_param_file(p7)
ck("UTF-8 中文路径读回一致", identical(enc2utf8(c7$fasta), cn_val))

# N-17（2026-09-25）：文件不存在的可读提示
err_missing <- tryCatch({ read_param_file(file.path(td, "no_such.txt")); "" },
                        error = function(e) conditionMessage(e))
ck("文件不存在给可读提示", grepl("parameter file not found", err_missing))

# P0-B（2026-09-25 代码审查实测）：GBK 编码配置（中文值）经 iconv(GBK→UTF-8)
# 分支读取，转回 UTF-8 与原文一致（UTF-8 原生构建的 R 下 readLines 直读会
# 报 "input string 1 is invalid UTF-8"）
p8 <- file.path(td, "cfg_gbk.txt")
cn_gbk <- "G:/数据/蛋白库.fasta"
Encoding(cn_gbk) <- "UTF-8"   # F-1（2026-09-26）：同上
gbk_raw <- iconv(paste0("fasta\t", cn_gbk, "\n"), from = "UTF-8", to = "GBK")
writeBin(charToRaw(gbk_raw), p8)
w8 <- character(0)
c8 <- withCallingHandlers(read_param_file(p8),
                          warning = function(w) { w8 <<- c(w8, conditionMessage(w)); invokeRestart("muffleWarning") })
ck("GBK 中文路径读回一致（iconv 分支）", identical(enc2utf8(c8$fasta), cn_gbk))

unlink(td, recursive = TRUE)
cat(sprintf("=== review_hardening_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
