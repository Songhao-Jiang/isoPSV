# unlabeled_two_part_test.R — 无标记样本两分法报告回归（v1.0 新功能；
# v1.0 口径修订：可信新肽段与标记与否无关，六重存疑规则同标准计算）
# 子进程运行真实 unlabeled_two_part_report.R，断言：
# Part1/Part2 按 in_fasta 拆分、diff_type=error_reason（判错行）、
# spectrum_verdict 混合谱口径、denovo 库外列、可信度三列按六重规则计算
# （判错/嵌合存疑、干净行 credible_novel=TRUE）、摘要文件格式。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

td <- tempfile("utp_"); dir.create(td)
proj <- file.path(td, "proj"); dir.create(proj)
res <- file.path(td, "res")
dir.create(file.path(res, "1.unlabeled_result"), recursive = TRUE)
dir.create(file.path(res, "merge"), recursive = TRUE)
fasta_f <- file.path(proj, "test.fasta")
writeLines(c(">sp|A|P1", "MAAAAAAAAAA", ">sp|B|P2", "MBBBBBBBBBB"), fasta_f)
writeLines(c(paste0("fasta\t", fasta_f), "label_type\tunlabeled"),
           file.path(proj, "config_step1.txt"))
writeLines(c(paste0("config_dir\t", proj),
             paste0("base_config\t", file.path(proj, "config_step1.txt"))),
           file.path(proj, "config_step3.txt"))

cls <- data.frame(
  File_Name = c("raw.1.1.2.dta", "raw.2.2.2.dta", "raw.3.3.2.dta", "raw.4.4.2.dta",
                "raw.5.5.2.dta", "raw.6.6.2.dta", "raw.7.7.2.dta"),
  a.Sequence = c("PEPTIDEK", "PEPTIDER", "PEPTIDES", "PEPTIDET", "PEPTIDEV",
                 "PEPTIDEW", "PEPTIDEX"),
  candidates = c("PEPTIDEK()", "PEPTIDER()", "PEPTIDES()", "PEPTIDET()",
                 "PEPTIDEV()", "PEPTIDEW()", "PEPTIDEX()"),
  candidate_PSM_score_by = c("95", "90", "85", "85", "70", "92", "80"),
  best_group = rep("1", 7),
  vote_says_wrong = c(1L, 0L, 1L, 0L, 1L, 0L, 0L),
  error_reason = c("modification_correction", "", "sequence_change", "",
                   "n_terminal_non_specific_cleavage", "", ""),
  modification_diff = c("winner extra: 4,Oxidation[M]", "", "", "", "", "", ""),
  chimeric = c(0L, 1L, 0L, 0L, 1L, 0L, 1L),
  conclusion_status = c("original result wrong", "chimeric spectrum",
                        "original result wrong", "no competition (strategies agree)",
                        "chimeric spectrum", "no competition (strategies agree)",
                        "chimeric spectrum"),
  stringsAsFactors = FALSE, check.names = FALSE)
write.table(cls, file.path(res, "1.unlabeled_result", "classification_A.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
mrg <- data.frame(
  File_Name = cls$File_Name,
  a.Proteins = c("sp|A|P1/", "sp|A|P1/", "sp|C|P3/", "", "NA", "sp|B|P2/", ""),
  e.Sequence = c("", "", "PEPTIDES", "", "PEPTIDEV", "", "PEPTIDEX"),
  e.DB_match = c("", "", "no_match", "", "no_match", "", "no_match"),
  stringsAsFactors = FALSE, check.names = FALSE)
write.table(mrg, file.path(res, "merge", "1.combine_result.spectra"),
            sep = "\t", row.names = FALSE, quote = FALSE)

rscript <- file.path(R.home("bin"), "Rscript.exe")
Sys.setenv(CONFIG3 = file.path(proj, "config_step3.txt"),
           RESULT_DIR = res)
st <- system2(rscript, c(shQuote(file.path(cd, "unlabeled_two_part_report.R"))),
              stdout = file.path(td, "run.log"), stderr = file.path(td, "run.log"))
ck("子进程运行成功", st == 0)
out2 <- file.path(res, "1.unlabeled_result")
p1 <- read.delim(file.path(out2, "1.annotation_correction_report.tsv"),
                 sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                 quote = "", comment.char = "", na.strings = character(0))
p2 <- read.delim(file.path(out2, "2.novel_peptide_competition_report.tsv"),
                 sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                 quote = "", comment.char = "", na.strings = character(0))
ck("两份报告同名产出且行数互补", nrow(p1) == 3L && nrow(p2) == 4L &&
   nrow(p1) + nrow(p2) == nrow(cls))
ck("Part1 = 原始蛋白在 FASTA 内的行", all(c("raw.1.1.2.dta", "raw.2.2.2.dta",
   "raw.6.6.2.dta") %in% p1$File_Name))
ck("Part2 = 不在 FASTA（含空蛋白/NA 行）", all(c("raw.3.3.2.dta", "raw.4.4.2.dta",
   "raw.5.5.2.dta", "raw.7.7.2.dta") %in% p2$File_Name))
ck("Part1 diff_type = 判错行 error_reason、其余为空",
   identical(as.character(p1$diff_type), c("modification_correction", "",
                                           "")) ||
   all(p1$diff_type[as.integer(p1$vote_says_wrong) == 1L] != ""))
ck("Part2 denovo 库外列透出（no_match）",
   any(p2$denovo_DB_match == "no_match"))
ck("Part2 spectrum_verdict 混合谱口径（chimeric→wrong_spectrum_mixture）",
   identical(as.character(p2$spectrum_verdict),
             ifelse(as.integer(p2$chimeric) == 1L, "wrong_spectrum_mixture", "")))
# v1.0 口径修订：可信度三列与标记与否无关，六重存疑规则同标准计算
ck("Part2 判错行 novel_doubt=error_reason 且非可信",
   p2$novel_doubt[p2$File_Name == "raw.3.3.2.dta"] == "sequence_change" &&
   !isTRUE(p2$credible_novel[p2$File_Name == "raw.3.3.2.dta"] == "TRUE"))
ck("Part2 判错优先于嵌合（raw.5 双标记取 error_reason）",
   p2$novel_doubt[p2$File_Name == "raw.5.5.2.dta"] ==
     "n_terminal_non_specific_cleavage")
ck("Part2 纯嵌合行（vote=0 chimeric=1）novel_doubt=chimeric",
   p2$novel_doubt[p2$File_Name == "raw.7.7.2.dta"] == "chimeric")
ck("Part2 干净行（无判错/嵌合/竞争/低分）credible_novel=TRUE",
   isTRUE(p2$credible_novel[p2$File_Name == "raw.4.4.2.dta"] == "TRUE") &&
   !isTRUE(p2$novel_doubt[p2$File_Name == "raw.4.4.2.dta"] != ""))
ck("Part2 证据/竞争者列齐备（cand1_sequence/n_frag/max_missing_run/closest）",
   all(c("cand1_sequence", "cand1_modification", "n_frag", "max_missing_run",
         "closest_annotated_competitor") %in% colnames(p2)) &&
   p2$cand1_sequence[p2$File_Name == "raw.4.4.2.dta"] == "PEPTIDET")
sum2 <- readLines(file.path(out2, "two_part_report_summary.txt"), warn = FALSE)
ck("摘要文件三列格式且两行齐全",
   length(sum2) == 3L && startsWith(sum2[1], "part\trows\theadline") &&
   all(vapply(sum2[2:3], function(l) length(strsplit(l, "\t", fixed = TRUE)[[1]]) == 3L,
              logical(1))))
unlink(td, recursive = TRUE)
cat(sprintf("=== unlabeled_two_part_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
