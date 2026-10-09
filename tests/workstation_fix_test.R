# workstation_fix_test.R — 工作站三份报告整改回归（v1.0）
# 覆盖：TAU_MAP 无标记口径（源断言）、Comet FASTA 注入（源断言）、_stray 移出
# unlabeled（源断言）、模板占位符与自洽（文件断言）、XTandem map 参数实转。
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

# 1) TAU_MAP：两个 PDF 模块均含 unlabeled 且防御式查找（源断言）
pdfA <- readLines(file.path(cd, "spectrum_pdf_export_A.R"), warn = FALSE)
pdfB <- readLines(file.path(cd, "spectrum_pdf_export.R"), warn = FALSE)
ck("PDF-A TAU_MAP 含 unlabeled（P1）",
   any(grepl("unlabeled = 0.0083", pdfA, fixed = TRUE)))
ck("PDF-A 多标记口径（label_tau_of(parse_label_spec)）",
   any(grepl("label_tau_of(parse_label_spec(label_type))", pdfA, fixed = TRUE)))
ck("PDF-B TAU_MAP 含 unlabeled（P1）",
   any(grepl("unlabeled = 0.0083", pdfB, fixed = TRUE)))
ck("PDF-B 默认 label_type 与主表统一为 unlabeled",
   any(grepl('else "unlabeled"', pdfB, fixed = TRUE)))

# 2) Comet 自动分支注入 FASTA + 模板占位符 + msfragger 模板自洽
st2 <- readLines(file.path(cd, "user_step2.R"), warn = FALSE)
ck("comet 自动分支注入 database_name = fasta_path（实弹修复）",
   any(grepl('paste0("database_name = ", fasta_path)', st2, fixed = TRUE)))
comet_t <- readLines(file.path(cd, "engine_templates", "comet_con.params"), warn = FALSE)
msf_t <- readLines(file.path(cd, "engine_templates", "msfragger_con.params"), warn = FALSE)
ck("comet 模板 database_name 为占位符（无残留绝对路径）",
   any(grepl("database_name = __FASTAPATH__", comet_t, fixed = TRUE)) &&
   !any(grepl("G:/JKY", comet_t, fixed = TRUE)))
ck("msfragger 模板 database_name 为占位符",
   any(grepl("database_name = __FASTAPATH__", msf_t, fixed = TRUE)) &&
   !any(grepl("G:/JKY", msf_t, fixed = TRUE)))
ck("msfragger 模板 precursor_mass_mode 自洽（selected）",
   any(grepl("precursor_mass_mode = selected", msf_t, fixed = TRUE)))

# 3) _stray 移出 unlabeled（denovo 加载器 EACCES 修复）
ck("_stray 目录位于 config_dir（unlabeled 之外）",
   any(grepl('file.path(cfg2$config_dir, "_stray_mgf")', st2, fixed = TRUE)) &&
   !any(grepl('file.path(unlab_dir, "_stray")', st2, fixed = TRUE)))

# 4) ensure_pfind_template 按 code_dir/sd0 解析（cwd 依赖修复）
rmg <- readLines(file.path(cd, "ReadMGF_script.R"), warn = FALSE)
ck("ensure_pfind_template 不再依赖 cwd（code_dir/sd0 解析）",
   any(grepl("tpl_rel <- file.path", rmg, fixed = TRUE)) &&
   any(grepl('exists("code_dir", envir = .GlobalEnv', rmg, fixed = TRUE)))

# 5) XTandem 适配器 map 参数实转：id→标题 vs label 回退
adp <- file.path(cd, "search_adapter_xtandem.R")
td <- tempfile("wsf_"); dir.create(td)
xml_f <- file.path(td, "t.xml")
writeLines(c(
  "<bioml>",
  '  <group id="2" mh="1201.6435" z="2" rt="30.2" expect="0.001"',
  '         label="sp|P1|PROT1 Putative" type="model">',
  '    <domain expect="0.001" seq="PEPTIDEK" start="1" end="10" />',
  "  </group>",
  '  <group id="3" mh="1000.0" z="3" rt="40.2" expect="0.002"',
  '         label="sp|P2|PROT2 Putative" type="model">',
  '    <domain expect="0.001" seq="PEPTIDER" start="1" end="9" />',
  "  </group>",
  "</bioml>"), xml_f, useBytes = TRUE)
map_f <- file.path(td, "t_renum_map.tsv")
writeLines(c("new_scan\toriginal_title",
             "2\traw1.1234.1234.2.0.dta",
             "3\traw1.1235.1235.3.0.dta"), map_f, useBytes = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
out_map <- file.path(td, "with_map.spectra")
st_a <- system2(rscript, c(shQuote(adp), shQuote(xml_f), shQuote(out_map),
                           shQuote(map_f)), stdout = file.path(td, "a.log"),
                stderr = file.path(td, "a.log"))
d_map <- if (st_a == 0 && file.exists(out_map)) read.delim(out_map, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
ck("XTandem 带 map 转换成功且 File_Name 取自 map（id→标题，P1）",
   st_a == 0 && nrow(d_map) == 2L &&
   all(c("raw1.1234.1234.2.0.dta", "raw1.1235.1235.3.0.dta") %in% d_map$File_Name))
out_nomap <- file.path(td, "no_map.spectra")
st_b <- system2(rscript, c(shQuote(adp), shQuote(xml_f), shQuote(out_nomap)),
                stdout = file.path(td, "b.log"), stderr = file.path(td, "b.log"))
d_nomap <- if (st_b == 0 && file.exists(out_nomap)) read.delim(out_nomap, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
ck("XTandem 无 map 保留 label 回退",
   st_b == 0 && nrow(d_nomap) == 2L &&
   all(grepl("PROT", d_nomap$File_Name)))
unlink(td, recursive = TRUE)

cat(sprintf("=== workstation_fix_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
