# sage_adapter_test.R — Sage 适配器 scannr 前缀回归（v1.0）
# 工作站根因（2026-09-26）：Sage 0.15.0-beta.2 的 scannr 写 "scan=N"，老版本写
# "index=N"——适配器须剥任意 ^[^=]*= 前缀，三种写法（scan=/index=/裸整数）都
# 能与 _renum_map.tsv 的 new_scan 对上。
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

td <- tempfile("sga_"); dir.create(td)
tsv_f <- file.path(td, "results.sage.tsv")
writeLines(c(
  "peptide\tproteins\tscannr\tcharge\tlabel\trank",
  "PEPTIDEK\tsp|P1|\tscan=2\t2\t1\t1",
  "PEPTIDER\tsp|P1|\tindex=3\t2\t1\t1",
  "PEPTIDES\tsp|P2|\t4\t2\t1\t1",
  "PEPTIDET\tsp|P2|\tscan=5\t2\t1\t1",
  "PEPTIDEV\tsp|P2|\tscan = 6\t2\t1\t1"),
  tsv_f, useBytes = TRUE)
map_f <- file.path(td, "t_renum_map.tsv")
writeLines(c("new_scan\toriginal_title",
             "2\traw1.1234.1234.2.0.dta",
             "3\traw1.1235.1235.3.0.dta",
             "4\traw1.1236.1236.4.0.dta",
             "5\traw1.1237.1237.5.0.dta",
             "6\traw1.1238.1238.6.0.dta"),
           map_f, useBytes = TRUE)
out_f <- file.path(td, "sage_out.spectra")
rscript <- file.path(R.home("bin"), "Rscript.exe")
st <- system2(rscript, c(shQuote(file.path(cd, "search_adapter_sage.R")),
                         shQuote(tsv_f), shQuote(map_f), shQuote(out_f)),
              stdout = file.path(td, "a.log"), stderr = file.path(td, "a.log"))
d <- if (st == 0 && file.exists(out_f)) read.delim(out_f, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
ck("四种 scannr 前缀（含等号带空格 scan = 6）全部经 map 解析出原始标题（v1.0）",
   st == 0 && nrow(d) == 5L &&
   all(c("raw1.1234.1234.2.0.dta", "raw1.1235.1235.3.0.dta",
         "raw1.1236.1236.4.0.dta", "raw1.1237.1237.5.0.dta",
         "raw1.1238.1238.6.0.dta") %in% d$File_Name))
ck("分级计数出现在日志（原始 5 | top-1 5 | 去诱饵 5 | map 命中 5）",
   grepl("原始 5 | top-1 5 | 去诱饵 5 | map 命中 5",
         paste(readLines(file.path(td, "a.log"), warn = FALSE), collapse = "\n")))
unlink(td, recursive = TRUE)
cat(sprintf("=== sage_adapter_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
