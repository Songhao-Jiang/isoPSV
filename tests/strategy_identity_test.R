# strategy_identity_test.R — 步骤 2 策略集合指纹回归（v1.0 功能排队①）
# 覆盖：同配置指纹稳定、增删策略改变指纹、类型推断（restrict/denovo/
# _punifind/boost/pnovo/默认 pfind_open）、显式类型覆盖、排序去重、侧车往返一致。
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

mk_cfg <- function(config_dir, search_file, search_file_type = NULL) {
  c(list(config_dir = config_dir, search_file = search_file),
    if (!is.null(search_file_type)) list(search_file_type = search_file_type) else list())
}

cfg <- mk_cfg("G:/proj",
              c("G:/proj/pfind_open/sub/x.spectra",
                "G:/proj/pfind_restrict/sub/y.spectra",
                "G:/proj/pfind_open/sub/z.spectra"))
ck("同配置指纹稳定", identical(strategy_identity_lines(cfg),
                               strategy_identity_lines(cfg)))
ck("类型推断：restrict 文件夹 → pfind_restrict",
   any(grepl("^pfind_restrict\t", strategy_identity_lines(cfg))))
ck("类型推断：默认 → pfind_open",
   any(grepl("^pfind_open\tG:/proj/pfind_open/sub/x.spectra$",
             strategy_identity_lines(cfg))))
# 增策略 → 指纹变化
cfg2 <- mk_cfg("G:/proj",
               c("G:/proj/pfind_open/sub/x.spectra",
                 "G:/proj/pfind_restrict/sub/y.spectra",
                 "G:/proj/pfind_open/sub/z.spectra",
                 "G:/proj/pUniFind_denovo/sub/d.csv"))
ck("新增策略改变指纹", !identical(strategy_identity_lines(cfg),
                                  strategy_identity_lines(cfg2)))
ck("denovo 文件夹推断 → punifind_denovo",
   any(grepl("^punifind_denovo\t", strategy_identity_lines(cfg2))))
# 显式类型覆盖推断
cfg3 <- mk_cfg("G:/proj", c("G:/proj/anything/x.spectra"), search_file_type = "other")
ck("显式类型覆盖推断", any(grepl("^other\tG:/proj/anything/x.spectra$",
                                 strategy_identity_lines(cfg3))))
# 排序去重
cfg4 <- mk_cfg("G:/proj", c("G:/proj/b.spectra", "G:/proj/a.spectra", "G:/proj/a.spectra"))
l4 <- strategy_identity_lines(cfg4)
ck("排序去重（重复条目只留一条）", length(l4) == 2L &&
   identical(l4, sort(l4)))
# 类型行数不匹配报错
cfg5 <- mk_cfg("G:/proj", c("G:/proj/a.spectra", "G:/proj/b.spectra"),
               search_file_type = "other")
ck("类型行数不匹配显式报错",
   inherits(try(strategy_identity_lines(cfg5), silent = TRUE), "try-error"))
# 侧车往返：写读一致（user_run 跳过判据的比对契约）
td <- tempfile("sti_"); dir.create(td)
writeLines(strategy_identity_lines(cfg2), file.path(td, "peptide.txt.strategies"),
           useBytes = TRUE)
ck("侧车往返逐字节一致（user_run 比对契约）",
   identical(trimws(readLines(file.path(td, "peptide.txt.strategies"), warn = FALSE)),
             strategy_identity_lines(cfg2)))
unlink(td, recursive = TRUE)

cat(sprintf("=== strategy_identity_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
