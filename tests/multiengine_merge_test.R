# 合成测试：6 策略（b-e pFind 系 + f maxquant）合并的字母列序/状态串/回帖槽
# 以及策略目录表校验（重复类型报错、other 通用槽、original 禁止）
options(stringsAsFactors = FALSE)
# 代码目录取自脚本自身位置（不硬编码中文路径：Windows 控制台源文件编码不一定是
# UTF-8，硬编码中文路径会被读成乱码致 setwd 失败）
args0 <- commandArgs(FALSE)
sd0 <- dirname(sub("^--file=", "", args0[grepl("--file=", args0)][1]))
code <- file.path("..", "code")
if (!dir.exists(code)) code <- file.path(sd0, "..", "code")
source(file.path(code, "compat.R"))
source(file.path(code, "ReadMGF_script.R"))

td <- tempfile("me_test_"); dir.create(td)
psm <- data.frame(File_Name = c("s1.dta", "s2.dta", "s3.dta"),
                  Charge = c(2, 2, 2),
                  Sequence = c("PEPTIDEK", "ANOTHERR", "COMMONPK"),
                  Modification = c("", "1,Oxidation[M];", ""),
                  stringsAsFactors = FALSE, check.names = FALSE)
write.table(psm, file.path(td, "psm.spectra"), sep = "\t", row.names = FALSE, quote = FALSE)
mk <- function(tag, seqs, mods = NULL) {
  d <- data.frame(File_Name = c("s1.dta", "s2.dta", "s3.dta"),
                  Sequence = seqs,
                  Modification = if (is.null(mods)) rep("", 3) else mods,
                  stringsAsFactors = FALSE, check.names = FALSE)
  write.table(d, file.path(td, paste0(tag, ".spectra")), sep = "\t",
              row.names = FALSE, quote = FALSE)
}
mk("open",     c("PEPTIDEK", "ANOTHERR", "COMMONPK"))          # b：与 a 一致
mk("restrict", c("PEPTIDKK", "ANOTHERR", ""))                  # c：s1 差异
mk("rescore",  c("PEPTIDEK", "ZZZOTHERR", "COMMONPK"))         # d：s2 差异
mk("denovo",   c("", "", "COMMONPK"))                          # e：s3 与 a 一致（回帖槽，无 fasta 时跳过）
mk("maxquant", c("MAXQUANTP", "", "COMMONPK"))                 # f：s1 另一差异

slots <- strategy_slots(c("pfind_open", "pfind_restrict", "punifind_rescore",
                          "punifind_denovo", "maxquant"))
stopifnot(identical(slots, c("b", "c", "d", "e", "f")))
# 目录表校验（2026-09-17 槽位定稿：g=comet、j=msfragger_boost、p=msfragger_open、
# other 通用槽 q..z）
stopifnot(strategy_slots("other") == "q")
stopifnot(identical(strategy_slots(c("other", "other")), c("q", "r")))
stopifnot(identical(strategy_slots("msfragger_open"), "p"))
stopifnot(identical(strategy_slots("comet"), "g"))
stopifnot(identical(strategy_slots("msfragger_boost"), "j"))
err1 <- try(strategy_slots(c("pfind_open", "pfind_open")), silent = TRUE)
stopifnot(inherits(err1, "try-error"))
err2 <- try(strategy_slots("original"), silent = TRUE)
stopifnot(inherits(err2, "try-error"))

res <- merge_search_results(file.path(td, "psm.spectra"),
                            file.path(td, paste0(c("open", "restrict", "rescore",
                                                   "denovo", "maxquant"),
                                                 ".spectra")),
                            file.path(td, "merge"),
                            fasta_file = NULL,
                            slot_letters = slots,
                            denovo_slots = "e")
m <- read.delim(file.path(td, "merge", "1.combine_result.spectra"), sep = "\t",
                check.names = FALSE, na.strings = character(0))
cn <- colnames(m)
stopifnot(all(c("a.Sequence", "b.Sequence", "c.Sequence", "d.Sequence",
                "e.Sequence", "f.Sequence") %in% cn))
seq_cols <- which(grepl("\\.Sequence$", cn))
stopifnot(identical(sub("\\.Sequence$", "", cn[seq_cols]), letters[1:6]))
# 状态串：s1 有 a/b/d/f(PEPTIDEK)+c(PEPTIDKK)+f(MAXQUANTP) -> 各字母
cat("s1 Peptide_state:", m$Peptide_state[1], "\n")
stopifnot(grepl("PEPTIDEK", m$a.Sequence[1]))
stopifnot(m$f.Sequence[1] == "MAXQUANTP")
stopifnot(identical(sort(strsplit(m$Peptide_state[1], ";")[[1]]),
                    c("abd", "c", "f")))
# s3：a/b/d/e 一致 -> Same
stopifnot(m$Peptide_state[3] == "Same")
# peptide.txt 覆盖 6 槽候选
pep <- read.delim(file.path(td, "merge", "peptide.txt"), sep = "\t",
                  check.names = FALSE, na.strings = character(0))
stopifnot("MAXQUANTP" %in% pep$peptide)
cat("合成 6 策略合并单测全部通过\n")
unlink(td, recursive = TRUE)
