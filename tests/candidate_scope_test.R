# candidate_scope_test.R — 全候选模式（candidate_scope=top5）与纯蛋白质基因组学
# 门控单元测试（v1.0，2026-09-27 用户口径）
# 覆盖：build_top5_datalist 逐谱 1..K 行装配（第 k 候选对齐、空位补齐、注释复制、
# 封顶 5）；merge_search_results top5 分支端到端；适配器全候选模式（Mascot CSV
# 子进程，ISOPSV_CAND_ALL=1 时保留全部候选）；蛋白质基因组学定侧口径
# （Proteins 全命中 FASTA＝已注释；空/未命中＝新肽候选侧）；默认模式不受影响。
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
td <- tempfile("cs5_"); dir.create(td)

# ---- 1. build_top5_datalist 装配（真实口径：psm 逐谱一行，槽文件可多候选）----
psm <- data.frame(
  File_Name = c("S1", "S2", "S3"),
  Sequence = c("AAA", "BBB", "CCC"),
  Modification = c("m1", "m2", ""),
  Charge = c(2L, 2L, 3L),
  Proteins = c("P1/", "P2/", "P3/"),
  stringsAsFactors = FALSE)
write.table(psm, file.path(td, "psm.spectra"), sep = "	",
            row.names = FALSE, quote = FALSE)
slotb <- data.frame(File_Name = c("S1", "S1", "S1", "S2", "S2", "S3"),
                    Sequence = c("b1", "b2", "b3", "c1", "c2", "e1"),
                    Modification = rep("", 6),
                    stringsAsFactors = FALSE)
write.table(slotb, file.path(td, "slotb.spectra"), sep = "	",
            row.names = FALSE, quote = FALSE)
slotc <- data.frame(File_Name = c("S1", "S3"), Sequence = c("d1", "f1"),
                    Modification = c("", ""),
                    stringsAsFactors = FALSE)
write.table(slotc, file.path(td, "slotc.spectra"), sep = "	",
            row.names = FALSE, quote = FALSE)

names(psm)[names(psm) %in% c("Sequence", "Modification")] <-
  c("a.Sequence", "a.Modification")
top5 <- build_top5_datalist(psm, c(file.path(td, "slotb.spectra"),
                                   file.path(td, "slotc.spectra")), c("b", "c"))
ck("top5 行数 = 3+2+1（K=max 各槽候选数）", nrow(top5) == 6L)
ck("S1 三行：k2 取 b 槽第 2 候选 b2、a 槽空位", top5$b.Sequence[2] == "b2" && top5$a.Sequence[2] == "")
ck("S1 k3：b 槽第 3 候选 b3", top5$b.Sequence[3] == "b3")
ck("S1 k1 的 c 槽 = d1、k2 空位", top5$c.Sequence[1] == "d1" && top5$c.Sequence[2] == "")
ck("S2 两行：k2 取 b 槽第 2 候选 c2", top5$b.Sequence[5] == "c2" && top5$a.Sequence[5] == "")
ck("S3 一行：b=e1、c=f1", top5$b.Sequence[6] == "e1" && top5$c.Sequence[6] == "f1")
ck("注释复制：Charge 全谱一致", all(top5$Charge == c(2L, 2L, 2L, 2L, 2L, 3L)))
ck("槽列齐备", all(c("a.Sequence", "a.Modification", "b.Sequence",
                     "b.Modification", "c.Sequence", "c.Modification") %in% names(top5)))

# ---- 2. merge_search_results top5 分支端到端 ----
outd5 <- file.path(td, "m5"); outd1 <- file.path(td, "m1")
res5 <- merge_search_results(file.path(td, "psm.spectra"),
                             c(file.path(td, "slotb.spectra"),
                               file.path(td, "slotc.spectra")),
                             outd5, slot_letters = c("b", "c"),
                             candidate_scope = "top5")
ck("top5 合并产物存在", file.exists(file.path(outd5, "1.combine_result.spectra")))
m5 <- read.delim(file.path(outd5, "1.combine_result.spectra"), sep = "	",
                 stringsAsFactors = FALSE, check.names = FALSE)
ck("top5 合并行数 = 6", nrow(m5) == 6L)
ck("top5 合并 Peptide_state 已计算", all(c("Peptide_state", "Modification_state") %in% names(m5)))
# 默认模式要求槽文件每谱一行（top-1），用单行件验证现状不受影响
slotb1 <- data.frame(File_Name = c("S1", "S2", "S3"),
                     Sequence = c("b1", "c1", "e1"),
                     Modification = rep("", 3),
                     stringsAsFactors = FALSE)
write.table(slotb1, file.path(td, "slotb_top1.spectra"), sep = "	",
            row.names = FALSE, quote = FALSE)
res1 <- merge_search_results(file.path(td, "psm.spectra"),
                             c(file.path(td, "slotb_top1.spectra"),
                               file.path(td, "slotc.spectra")),
                             outd1, slot_letters = c("b", "c"),
                             candidate_scope = "filtered")
m1 <- read.delim(file.path(outd1, "1.combine_result.spectra"), sep = "	",
                 stringsAsFactors = FALSE, check.names = FALSE)
ck("默认模式合并行数 = psm 行数（3，逐槽 top-1）", nrow(m1) == 3L)
ck("默认模式 b 槽首行 = b1", m1$b.Sequence[1] == "b1")
ck("两模式列集合一致", setequal(names(m5), names(m1)))

# ---- 3. Mascot 适配器全候选模式（子进程 + 环境变量）----
mascot_csv <- file.path(td, "mascot_export.csv")
writeLines(c("pep_query,pep_rank,pep_score,pep_scan_title,pep_seq,pep_var_mod",
             "1,1,50,T1,AAA,",
             "1,2,40,T1,BBB,",
             "2,1,60,T2,CCC,",
             "2,2,30,T2,DDD,"),
           mascot_csv, useBytes = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
adapter <- file.path(cd, "search_adapter_mascot.R")
out_t1 <- file.path(td, "mascot_top1.spectra")
out_all <- file.path(td, "mascot_all.spectra")
Sys.setenv(ISOPSV_CAND_ALL = "0")
system2(rscript, c(shQuote(adapter), shQuote(mascot_csv), shQuote(out_t1)),
        stdout = TRUE, stderr = TRUE)
Sys.setenv(ISOPSV_CAND_ALL = "1")
system2(rscript, c(shQuote(adapter), shQuote(mascot_csv), shQuote(out_all)),
        stdout = TRUE, stderr = TRUE)
Sys.unsetenv("ISOPSV_CAND_ALL")
d_t1 <- read.delim(out_t1, sep = "\t", stringsAsFactors = FALSE)
d_all <- read.delim(out_all, sep = "\t", stringsAsFactors = FALSE)
ck("默认模式：Mascot 每谱 top-1（2 行）", nrow(d_t1) == 2L)
ck("全候选模式：Mascot 保留全部候选（4 行）", nrow(d_all) == 4L)
ck("全候选模式：每谱内高分在前（AAA 先于 BBB、CCC 先于 DDD）",
   match("AAA", d_all$Sequence) < match("BBB", d_all$Sequence) &&
   match("CCC", d_all$Sequence) < match("DDD", d_all$Sequence))

# ---- 4. 蛋白质基因组学定侧口径 ----
fa <- c(">P1", "AAAAAA", ">P3", "CCCCCC")
fa_f <- file.path(td, "mini.fasta")
writeLines(fa, fa_f, useBytes = TRUE)
fa_acc <- unique(names(read_fasta(fa_f)))
acc_of <- function(pr) {
  if (is.na(pr) || !nzchar(pr) || pr %in% c("-", "NA")) return(character(0))
  p <- trimws(strsplit(pr, "/", fixed = TRUE)[[1]])
  p[nzchar(p)]
}
in_fasta <- vapply(c("P1/", "P3/", "P9/", "", "P1/P9/"), function(pr) {
  acc <- acc_of(pr)
  length(acc) > 0 && any(acc %in% fa_acc)
}, logical(1))
novel <- !in_fasta
ck("全命中 FASTA＝已注释（P1、P3）", isTRUE(in_fasta[1]) && isTRUE(in_fasta[2]))
ck("全部未命中＝新肽侧（P9）", isTRUE(novel[3]))
ck("空蛋白＝新肽候选侧", isTRUE(novel[4]))
ck("混合（P1/P9）按两分法口径＝已注释", isTRUE(in_fasta[5]))

# ---- 5. 第三步门控配置来源（v1.0 勘误回归：config_step3 的 base_config
# 直接指向 config_step1，故 base_cfg 即第一步；不得再对 base_cfg 找 base_config）----
s3 <- paste(readLines(file.path(cd, "user_step3.R"), warn = FALSE), collapse = "
")
ck("user_step3 门控循环含 base_cfg（=config_step1）与 cfg3", grepl("list\\(base_cfg, cfg3\\)", s3))
ck("user_step3 不再对 base_cfg 二次找 base_config", !grepl("cfg1_pg", s3))

# ---- 6. 引擎输出深度改写（v1.0：CAND_ALL 下 comet/msfragger 输出前 5）----
s2 <- paste(readLines(file.path(cd, "user_step2.R"), warn = FALSE), collapse = "
")
ck("comet num_output_lines 改写", grepl("num_output_lines = 5", s2))
ck("msfragger output_report_topN 改写", grepl("output_report_topN = 5", s2))

# ---- 7. step1 优先覆盖告警 + config_step3 重建保留（v1.0 工作站反馈）----
w <- capture.output(
  r <- cfg_get_preferred(list(punifind_exe = "OLD"), list(punifind_exe = "NEW"), "punifind_exe", ""),
  type = "output")
ck("step1/step2 同键不同值 → 显式告警", any(grepl("配置告警", w)) || any(grepl("config warning", w)))
ck("告警不改口径（仍取 step1）", identical(r$value, "OLD"))
ck("同值不告警", !any(grepl("告警|warning", capture.output(cfg_get_preferred(list(punifind_exe="X"), list(punifind_exe="X"), "punifind_exe", ""), type="output"))))
# make_step3_config 保留 proteogenomics_mode
td2 <- tempfile("cfg3_")
dir.create(td2)
cfg2f <- file.path(td2, "config_step2.txt")
writeLines(c(paste0("base_config	", td2, "/config_step1.txt"), paste0("config_dir	", td2)), cfg2f, useBytes = TRUE)
writeLines("proteogenomics_mode	1", file.path(td2, "config_step1.txt"), useBytes = TRUE)
c2 <- read_param_file(cfg2f)
out3 <- file.path(td2, "config_step3.txt")
writeLines(c("base_config	x", "proteogenomics_mode	1"), out3, useBytes = TRUE)
suppressWarnings(try(make_step3_config(cfg2f, c2, out3), silent = TRUE))
ck("config_step3 重建保留 proteogenomics_mode", any(grepl("^proteogenomics_mode", readLines(out3, warn = FALSE))))

cat("汇总：", npass, " PASS / ", nfail, " FAIL\n", sep = "")
if (nfail > 0L) quit(status = 1L)
