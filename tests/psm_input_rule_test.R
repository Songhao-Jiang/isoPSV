# psm_input_rule_test.R — 第一步 psm_file 输入口径回归（v1.0，口径修订）
# 第一步接受任何带 File_Name 列的制表符结果表（pFind 家族表、交换格式、
# 其他引擎表格输出）；无 Label 列时无标记样本全收、标记样本须 psm_all_unlabeled=1；
# 非表格格式（XML）拒收。
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

td <- tempfile("psr_"); dir.create(td)
mk <- function(name, lines) {
  f <- file.path(td, name)
  writeLines(lines, f, useBytes = TRUE)
  f
}

# 1) pFind 家族两列表 → 通过（Label 语义）
f_ok <- mk("pfind.spectra", c("File_Name\tLabel\tScore", "raw.1.1.2.dta\t1\t0.9"))
r_ok <- check_psm_header(f_ok)
ck("pFind 家族两列表通过（Label 语义）",
   identical(unname(r_ok$i_fn), 1L) && identical(unname(r_ok$i_lab), 2L) &&
   !isTRUE(r_ok$all_unlabeled))

# 2) 交换格式（File_Name/Sequence/Modification，无 Label）→ 无标记样本通过（全收）
f_ada <- mk("adapter.spectra", c("File_Name\tSequence\tModification",
                                 "raw.1.1.2.dta\tPEPTIDEK\t"))
r_ada <- check_psm_header(f_ada, label_type = "unlabeled")
ck("交换格式在无标记样本下通过且全收（all_unlabeled=TRUE）",
   identical(unname(r_ada$i_fn), 1L) && isTRUE(r_ada$all_unlabeled) &&
   is.na(r_ada$i_lab))

# 3) 无 Label 列 + 标记样本 → 报错要求 psm_all_unlabeled
r_lab <- tryCatch(check_psm_header(f_ada, label_type = "n15"),
                  error = function(e) conditionMessage(e))
ck("无 Label 列 + 标记样本报错并给 psm_all_unlabeled 指引",
   grepl("psm_all_unlabeled", r_lab))

# 4) 无 Label 列 + 标记样本 + psm_all_unlabeled=1 → 通过
r_conf <- check_psm_header(f_ada, label_type = "n15", psm_all_unlabeled = TRUE)
ck("psm_all_unlabeled=1 显式确认后通过", isTRUE(r_conf$all_unlabeled))

# 5) 其他引擎表格输出（如 Comet/Sage 的 tsv，含 File_Name）→ 无标记样本通过
f_eng <- mk("engine.tsv", c("File_Name\tscore\tcharge", "raw.1.1.2.dta\t0.9\t2"))
r_eng <- check_psm_header(f_eng, label_type = "unlabeled")
ck("其他引擎表格输出通过（File_Name 列即可）", isTRUE(r_eng$all_unlabeled))

# 6) 非表格格式（XML）→ 拒收并指引走适配器
f_xml <- mk("raw.xml", c("<?xml version=\"1.0\"?>", "<bioml>"))
r_xml <- tryCatch(check_psm_header(f_xml), error = function(e) conditionMessage(e))
ck("非表格格式（XML）拒收且指引适配器转换",
   grepl("XML", r_xml) && grepl("适配器", r_xml))

# 7) 无 File_Name 列 → 拒收并说明新口径
f_fn <- mk("fnonly.txt", c("Title\tScore", "raw.1.1.2.dta\t0.9"))
r_fn <- tryCatch(check_psm_header(f_fn), error = function(e) conditionMessage(e))
ck("缺 File_Name 列拒收且报错含新口径说明", grepl("File_Name", r_fn))

# 8) 格式嗅探（v1.0 修订二：可转换的不拒收，先经适配器）
ck("嗅探：表格直接可用", identical(psm_sniff_engine(f_ok), ""))
ck("嗅探：bioml → xtandem",
   identical(psm_sniff_engine(mk("x.xml", c("<bioml><note/></bioml>"))), "xtandem"))
ck("嗅探：msms_pipeline_analysis → msfragger",
   identical(psm_sniff_engine(mk("p.xml", c("<msms_pipeline_analysis/>"))), "msfragger"))
ck("嗅探：mzIdentML → mzidentml",
   identical(psm_sniff_engine(mk("m.xml", c("<mzIdentML/>"))), "mzidentml"))
ck("嗅探：未知 XML → unknown",
   identical(psm_sniff_engine(mk("u.xml", c("<something/>"))), "unknown"))
# 8b) P2（2026-09-26，工作站 v1.0 复测）：真实 X!Tandem XML 首行是
# <?xml …?> 声明（<bioml> 在下一行）——嗅探须剥声明、看前 3 行
ck("嗅探：<?xml 声明首行 + 第二行 <bioml → xtandem",
   identical(psm_sniff_engine(
     mk("x2.xml", c("<?xml version=\"1.0\"?>", "<bioml><note/></bioml>"))),
     "xtandem"))
ck("嗅探：<?xml 声明 + 空行 + <bioml → xtandem",
   identical(psm_sniff_engine(
     mk("x3.xml", c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", "", "<bioml/>"))),
     "xtandem"))
ck("嗅探：<?xml 声明 + <msms_pipeline_analysis → msfragger",
   identical(psm_sniff_engine(
     mk("p2.xml", c("<?xml version=\"1.0\"?>", "<msms_pipeline_analysis/>"))),
     "msfragger"))

# 9) 表格 pass-through
cfg0 <- list(config_dir = td)
ck("psm_prepare_input 表格原样通过",
   identical(psm_prepare_input(f_ok, cfg0, cd), f_ok))

# 10) Mascot CSV + psm_type=mascot → 自动转换出交换格式（含 File_Name 列）
mascot_lines <- c(
  "prot_hit_num,prot_acc,prot_desc,prot_score,prot_matches,prot_matches_sig,prot_sequences,prot_sequences_sig,prot_cover,prot_len,pep_query,pep_rank,pep_isbold,pep_isunique,pep_exp_mz,pep_exp_mr,pep_exp_z,pep_calc_mr,pep_delta,pep_start,pep_end,pep_miss,pep_score,pep_homol,pep_ident,pep_expect,pep_res_before,pep_seq,pep_res_after,pep_var_mod,pep_var_mod_pos,pep_summed_mod_pos,pep_num_match,pep_scan_title",
  "1,sp|P12345|GENE1,GENE1 protein,58,5,4,3,2,12,200,5,1,1,1,601.8324,1201.6502,2,1201.6435,0.0067,10,21,0,58.1,45,45,1.2e-05,K,AMCLPEPTIDEK,R,Oxidation (M),010000000000,,10,raw1.1234.1234.2.0.dta")
f_mas <- mk("mascot_export.csv", mascot_lines)
cfg_m <- list(config_dir = td, psm_type = "mascot")
conv <- psm_prepare_input(f_mas, cfg_m, cd)
conv_first <- strsplit(readLines(conv, n = 1, warn = FALSE), "\t")[[1]]
conv_ok <- file.exists(conv) && "File_Name" %in% conv_first
ck("Mascot CSV 自动转换出交换格式（File_Name 列）", conv_ok)

# 11) msfragger（需重编号映射）→ 定向报错指引走 search_file 槽位
f_pxml <- mk("msf.pepXML", c("<msms_pipeline_analysis/>"))
r_msf <- tryCatch(psm_prepare_input(f_pxml, cfg0, cd),
                  error = function(e) conditionMessage(e))
ck("msfragger XML 定向报错并指引 search_file", grepl("search_file", r_msf))

# 12) 未知格式且无 psm_type → 报错给出可自动转换清单
f_unk <- mk("unk.dat", c("BINARYDATA"))
r_unk <- tryCatch(psm_prepare_input(f_unk, cfg0, cd),
                  error = function(e) conditionMessage(e))
ck("未知格式无 psm_type 报错并给出可转换清单", grepl("psm_type", r_unk))

# 13) P2（2026-09-26，工作站 v1.0 复测）：<?xml 声明首行的真实 X!Tandem
# XML——自动嗅探识别后走转换，且转换分支自动传同目录 <同名>_renum_map.tsv
# （File_Name 按 map 的 id→标题解析，而非 label 蛋白描述）
xml_xt <- mk("tx.xml", c(
  '<?xml version="1.0" encoding="UTF-8"?>',
  "<bioml>",
  '  <group id="2" mh="1201.6435" z="2" rt="30.2" expect="0.001"',
  '         label="sp|P1|PROT1 Putative" type="model">',
  '    <domain expect="0.001" seq="PEPTIDEK" start="1" end="10" />',
  "  </group>",
  "</bioml>"))
map_tx <- file.path(td, "tx_renum_map.tsv")
writeLines(c("new_scan\toriginal_title", "2\traw1.9999.9999.2.0.dta"),
           map_tx, useBytes = TRUE)
cfg_tx <- list(config_dir = td)
conv_tx <- psm_prepare_input(xml_xt, cfg_tx, cd)
d_tx <- if (file.exists(conv_tx)) read.delim(conv_tx, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
ck("声明首行 XML：嗅探识别 + 自动转换带 map 出真实标题（P2）",
   nrow(d_tx) == 1L && "raw1.9999.9999.2.0.dta" %in% d_tx$File_Name &&
   !grepl("PROT", d_tx$File_Name[1]))
unlink(td, recursive = TRUE)

cat(sprintf("=== psm_input_rule_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
