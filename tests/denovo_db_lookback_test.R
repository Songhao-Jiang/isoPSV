# =============================================================================
# denovo_db_lookback_test.R - denovo 数据库回帖单元测试（用户 2026-09-03 规格④）
# denovo_db_lookback_test.R - unit tests for the denovo database look-back
# 运行：在本目录 Rscript denovo_db_lookback_test.R（预期全部 PASS）
# =============================================================================
# 代码目录定位（2026-09-18 加固）：约定在本目录运行；若从其他目录启动，
# 按 --file= 的目录定位。不用 normalizePath（受限控制台下中文路径会被读成
# 乱码，normalizePath/setwd 均失败）。
args0 <- commandArgs(FALSE)
sd0 <- dirname(sub("^--file=", "", args0[grepl("--file=", args0)][1]))
code <- file.path("..", "code")
if (!dir.exists(code)) code <- file.path(sd0, "..", "code")
source(file.path(code, "compat.R"))
source(file.path(code, "ReadMGF_script.R"))

# 临时 FASTA：
#   T1 = MAGGLPEPTIDEK（含真 GG）
#   T2 = MPEPTIDEKLAR（L）
#   T3 = MANLPEPTIDEK（N 位，供 N→GG 等价）
#   T4 = REV_ 反库蛋白（应被剔除，不作为命中源）
fa <- c(
  ">T1",
  "KMAGGLPEPTIDEKPK",
  ">T2",
  "KMPEPTIDEKLARP",
  ">T3",
  "MANLPEPTIDEK",
  ">REV_T4",
  "MGGLPEPTIDEK",
  ">T5",
  "KSANLPEPTIDEKP"
)
fa_f <- tempfile(fileext = ".fasta")
writeLines(fa, fa_f)

dl <- data.frame(
  File_Name = c("a.dta", "b.dta", "c.dta", "d.dta", "e.dta", "f.dta"),
  "e.Sequence" = c("MPEPTIDEKIAR",      # I→L 直接
                   "MANLPEPTIDEK",      # N→GG 后命中 T1（替换库）
                   "QQQPEPTIDEK",       # 无匹配 -> no_match 清空
                   "",                  # 空槽
                   "MGGLPEPTIDEK",      # 与 REV_T4 同串：反库剔除后应 no_match
                   "KSANLPEPTIDEK"),    # N 位蛋白（T5）：N→GG 后命中、标注跨 N 位
  "e.Modification" = rep("", 6),
  check.names = FALSE, stringsAsFactors = FALSE)

out <- denovo_db_lookback(dl, fa_f)
st <- out[["e.DB_match"]]
pass1 <- st[1] == "match:transformed" && out[1, "e.Sequence"] == "MPEPTIDEKLAR" &&
         out[1, "e.Proteins"] == "T2/" &&
         grepl("^1,K,P/", out[1, "e.Positions"])           # 起始位 2（pFind 0 起=1）、前 K、后 P
pass2 <- st[2] == "match:transformed" && out[2, "e.Sequence"] == "MAGGLPEPTIDEK" &&
         out[2, "e.Proteins"] == "T1/" &&
         grepl("^1,K,P/", out[2, "e.Positions"])           # 起始位 2（pFind 0 起=1）、前 K、后 P
pass3 <- st[3] == "no_match" && out[3, "e.Sequence"] == "" &&
         out[3, "e.Modification"] == "" &&
         out[3, "e.Proteins"] == "" && out[3, "e.Positions"] == ""
pass4 <- st[4] == ""
pass5 <- st[5] == "no_match"                                # 反库蛋白不作为命中源
pass6 <- st[6] == "match:transformed" && out[6, "e.Sequence"] == "KSANLPEPTIDEK" &&
         out[6, "e.Proteins"] == "T5/" &&
         grepl("^0,-,P/", out[6, "e.Positions"])           # 原始库序列（N 回原）：起于蛋白首（0 起 0）

# revise 模式测试（规格①②③）：none 不动 / terminal 跳过 SAP 交换 / full 全做
rv <- data.frame(File_Name = "x.dta",
                 Sequence = c("AMSTGK", "PEPK", "ADPVK"),
                 Modification = c("0,Met-;", "", "2,Asp->Asn[A];"),
                 stringsAsFactors = FALSE)
r0 <- revise_search_result(rv, "none")
p0 <- identical(r0$Sequence, rv$Sequence) && identical(r0$Modification, rv$Modification)
rt <- revise_search_result(rv, "terminal")
p1 <- rt$Sequence[1] == "MSTGK" && rt$Modification[1] == "0,Acetyl[AnyN-term];" &&
      rt$Sequence[3] == "ADPVK" && rt$Modification[3] == "2,Asp->Asn[A];"  # 交换保留
rf <- revise_search_result(rv, "full")
# 2026-09-10 用户口径：X->Y 全程保持 SAP 修饰形式，不再施进序列
# user ruling 2026-09-10: X->Y stays in the SAP-modification form, never
# substituted into the sequence in any revise mode
p2 <- rf$Sequence[3] == "ADPVK" && grepl("Asp->Asn", rf$Modification[3], fixed = TRUE)

cat(sprintf("[%s] I→L 直接＋Proteins/Positions 标注  seq=%s prot=%s pos=%s\n",
            if (pass1) "PASS" else "FAIL", out[1, "e.Sequence"], out[1, "e.Proteins"], out[1, "e.Positions"]))
cat(sprintf("[%s] N→GG 替换库命中＋原始坐标标注   seq=%s prot=%s pos=%s\n",
            if (pass2) "PASS" else "FAIL", out[2, "e.Sequence"], out[2, "e.Proteins"], out[2, "e.Positions"]))
cat(sprintf("[%s] 无匹配清空四字段                 status=%s\n", if (pass3) "PASS" else "FAIL", st[3]))
cat(sprintf("[%s] 空槽状态空                       status=[%s]\n", if (pass4) "PASS" else "FAIL", st[4]))
cat(sprintf("[%s] REV_ 反库剔除（不作为命中源）    status=%s\n", if (pass5) "PASS" else "FAIL", st[5]))
cat(sprintf("[%s] revise none 模式（限定式无后处理）\n", if (p0) "PASS" else "FAIL"))
cat(sprintf("[%s] revise terminal 模式（rescore 仅末端丢失）\n", if (p1) "PASS" else "FAIL"))
cat(sprintf("[%s] revise full 模式（open 含 SAP 交换）\n", if (p2) "PASS" else "FAIL"))
unlink(fa_f)
cat(if (all(pass1, pass2, pass3, pass4, pass5, p0, p1, p2)) "全部 PASS\n" else "存在 FAIL\n")
if (!all(pass1, pass2, pass3, pass4, pass5, p0, p1, p2)) quit(status = 1L)   # P0-6（2026-09-26）：失败必须非零退出，否则 run_all 假 PASS
