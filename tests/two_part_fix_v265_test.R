# two_part_fix_v265_test.R — v1.0 两分法报告/新肽段判定四缺陷修复回归（2026-09-28）
# 覆盖：fa_annotated_proteins 三口径（accession/序列级 I-L 归一/CON_ 污染库）、
# a.Proteins 回退、分侧序列级、cand1 自身参考库检查、step3/user_run PG 分侧统一。
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

# ---- 1. fa_annotated_proteins 单元（三口径） ----
tmp <- tempfile("fa_"); dir.create(tmp)
fa <- file.path(tmp, "ref.fasta")
writeLines(c(">WP_001234 protein one", "AAAKKKDE",
             ">NP_002345 protein two", "GGLLLRPP"), fa)
h <- fa_annotated_proteins
r1 <- h(c("WP_001234/", "orf|1|-|1048161-1049147/", "orf|1|-|1-90/",
          "CON_Trypsin/", "orf|9|-|5-40/", ""),
        c("", "AAAKKKDE", "GGILIRPP", "XYZZ", "NOTPRESENT", "AAAKKKDE"),
        fa)
ck("accession 命中 → TRUE", isTRUE(r1[1]))
ck("六框库 orf| + 序列在参考库 → TRUE", isTRUE(r1[2]))
ck("序列级 I/L 归一命中 → TRUE", isTRUE(r1[3]))
ck("CON_ 污染库 → TRUE（非新肽段）", isTRUE(r1[4]))
ck("accession 与序列都未命中 → FALSE", identical(r1[5], FALSE))
ck("蛋白列为空 → FALSE（保守归新肽段侧）", identical(r1[6], FALSE))
ck("fa_file 缺失时 CON_ 仍判 TRUE", identical(h(c("CON_Trypsin/"), "X", "")[1], TRUE))
ck("fa_file 缺失时 accession 判 FALSE", identical(h(c("WP_001234/"), "X", "")[1], FALSE))
unlink(tmp, recursive = TRUE)

# ---- 2. 静态断言 ----
rd <- paste(readLines(file.path(cd, "ReadMGF_script.R"), warn = FALSE), collapse = "\n")
ck("共享帮助器含三口径分支",
   grepl("startsWith(p, \"CON_\")", rd, fixed = TRUE) &&
   grepl('gsub("I", "L", s)', rd, fixed = TRUE) &&
   grepl("p %in% fa_acc", rd, fixed = TRUE))

rp <- paste(readLines(file.path(cd, "unlabeled_two_part_report.R"), warn = FALSE),
            collapse = "\n")
ck("缺陷1：a.Proteins 回退 Proteins",
   grepl('else if ("Proteins" %in% colnames(mrg)) proteins <- as.character(mrg$Proteins)[mm]',
         rp, fixed = TRUE))
ck("缺陷2：分侧补 a.Sequence 序列级匹配",
   grepl("a_seqs <- as.character(mrg$a.Sequence)[mm]", rp, fixed = TRUE) &&
   grepl("in_fasta <- in_fasta | con_acc | seq_in", rp, fixed = TRUE))
ck("缺陷4：CON_ 前缀判非新肽段",
   grepl('any(startsWith(acc, "CON_"))', rp, fixed = TRUE))
ck("缺陷3：cand1 自身参考库检查（annotated_cand1）",
   grepl('doubt_k <- paste0("annotated_cand1:", s1)', rp, fixed = TRUE) &&
   grepl("nzchar(s1) && annotated_seq_ul(s1)", rp, fixed = TRUE))
ck("annotated_seq_ul 空串守卫", grepl("if (!nzchar(s)) return(FALSE)", rp, fixed = TRUE))
ck("新肽段侧 0 行时列初始化按行数扩空（实弹暴露）",
   grepl("p2$cand1_sequence <- rep(\"\", nrow(p2))", rp, fixed = TRUE) &&
   grepl("p2$credible_novel <- rep(\"\", nrow(p2))", rp, fixed = TRUE))

s3 <- paste(readLines(file.path(cd, "user_step3.R"), warn = FALSE), collapse = "\n")
ck("step3 PG 分侧用共享帮助器",
   grepl("fa_annotated_proteins(as.character(prot_col), seq_col", s3, fixed = TRUE))
ur <- paste(readLines(file.path(cd, "user_run.R"), warn = FALSE), collapse = "\n")
ck("user_run PG 早退用共享帮助器",
   grepl("fa_annotated_proteins(as.character(psmu$Proteins), seq_u", ur, fixed = TRUE))

cat("==== two_part_fix_v265_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
