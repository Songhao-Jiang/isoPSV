# label_custom_rows_test.R — v1.0 第二类自定义修饰行回归（2026-09-28）
# 覆盖：parse_comp_str/comp_mass（元素组成解析与质量计算，含同位素与负计数）、
# label_pos_normalize（五位置码 + 中文短语）、parse_label_row（4 字段与非法输入）、
# label_mod_masses 混合（自定义行 vs 旧预设名质量等价）、engine_mod_token/
# label_mod_lines（comet/msfragger 端部语法折算与编号续接）、sage 端部跳过、
# append_label_mods/strip_label_mods（自定义行的理论侧定位与剥除）、
# server.py label_type 多行修复与 index.html 表单的静态断言。
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
err <- function(expr) {
  inherits(tryCatch(force(expr), error = function(e) e), "error")
}

# ---- 1. parse_comp_str / comp_mass ----
c1 <- parse_comp_str("H(2)C(2)O(1)")
ck("H(2)C(2)O(1) 计数", identical(unname(c1["H"]), 2L) &&
   identical(unname(c1["C"]), 2L) && identical(unname(c1["O"]), 1L))
ck("C(-6)13C(6) 质量 6.020129",
   identical(sprintf("%.6f", comp_mass("C(-6)13C(6)")), "6.020129"))
ck("H(-3)N(-1) 质量 -17.026549",
   identical(sprintf("%.6f", comp_mass("H(-3)N(-1)")), "-17.026549"))
ck("2H(4)C(2) 质量 32.056407",
   identical(sprintf("%.6f", comp_mass("2H(4)C(2)")), "32.056407"))
ck("H(-2)2H(6)C(2) 质量 34.068960",
   sprintf("%.5f", comp_mass("H(-2)2H(6)C(2)")) == "34.06896")
ck("裸式 H2O 质量 18.010565",
   identical(sprintf("%.6f", comp_mass("H2O")), "18.010565"))
ck("同名元素累加 C(2)C(1) = C(3)",
   identical(unname(parse_comp_str("C(2)C(1)")["C"]), 3L))
ck("零计数剔除 C(0) → 空", length(parse_comp_str("C(0)")) == 0L &&
   comp_mass("C(0)") == 0)
ck("未知元素报错", err(parse_comp_str("H(2)Xy(1)")))
ck("语法非法报错", err(parse_comp_str("H(2)%%%")))
ck("空组成报错", err(comp_mass("  ")))
ck("Phospho H(1)O(3)P(1) ≈ 79.966331（容差 1e-4）",
   abs(comp_mass("H(1)O(3)P(1)") - 79.966331) < 1e-4)

# ---- 2. label_pos_normalize ----
ck("中文位置短语", identical(label_pos_normalize("任意位置"), "ANY") &&
   identical(label_pos_normalize("肽段N端"), "PEP_N") &&
   identical(label_pos_normalize("肽段C端"), "PEP_C") &&
   identical(label_pos_normalize("蛋白N端"), "PRO_N") &&
   identical(label_pos_normalize("蛋白C端"), "PRO_C"))
ck("pFind 风格别名", identical(label_pos_normalize("anyN-term"), "PEP_N") &&
   identical(label_pos_normalize("ProteinC-term"), "PRO_C") &&
   identical(label_pos_normalize("NORMAL"), "ANY"))
ck("未知位置报错", err(label_pos_normalize("随便哪里")))

# ---- 3. parse_label_row ----
r1 <- parse_label_row("NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)")
ck("NeuCode8 行四字段", identical(r1$name, "NeuCode8") &&
   identical(r1$residue, "K") && identical(r1$pos, "ANY") &&
   identical(sprintf("%.6f", r1$mass), "8.014199"))
r2 <- parse_label_row("MeN|肽段N端|*|H(4)C(2)")
ck("中文位置 + 全部氨基酸", identical(r2$pos, "PEP_N") &&
   identical(r2$residue, "*") &&
   identical(sprintf("%.6f", r2$mass), "28.031300"))
ck("氨基酸集合非法报错", err(parse_label_row("Bad|ANY|BZ|H(1)")))
ck("字段数不足报错", err(parse_label_row("Bad|ANY|K")))
ck("名称含分号报错", err(parse_label_row("Ba;d|ANY|K|H(1)")))
ck("is_custom_label_row 判别", is_custom_label_row("a|ANY|K|H(1)") &&
   !is_custom_label_row("Label_13C(6)[K]"))

# ---- 4. label_mod_masses：自定义行与旧预设等价/兼容 ----
lm <- label_mod_masses(c("Label_13C(6)15N(2)[K]", "NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)"))
ck("预设与自定义质量等价（8.014199）",
   nrow(lm) == 2L && identical(sprintf("%.6f", lm$mass[1]), sprintf("%.6f", lm$mass[2])) &&
   identical(sprintf("%.6f", lm$mass[1]), "8.014199"))
ck("旧预设残基/位置不变（K/ANY）", identical(lm$residue[1], "K") &&
   identical(lm$pos[1], "ANY") && identical(lm$name[1], "Label_13C(6)15N(2)[K]"))
lmt <- label_mod_masses("Dimethyl_2H(4)[AnyN-term]")
ck("旧 [AnyN-term] 折算 PEP_N + 空残基", identical(lmt$pos, "PEP_N") &&
   identical(lmt$residue, "") && identical(sprintf("%.6f", lmt$mass), "32.056407"))
lm0 <- label_mod_masses(c("K080[K]", "NoSuch[Z]"))
ck("旧名查表行为保持（NA 守卫依赖）", !is.na(lm0$mass[1]) && is.na(lm0$mass[2]))

# ---- 5. engine_mod_token / label_mod_lines 端部折算 ----
ck("comet 端部裸 n/c", identical(engine_mod_token("K", "PEP_N", "comet"), "n") &&
   identical(engine_mod_token("*", "PEP_C", "comet"), "c"))
ck("comet 任意位置 = 残基集合", identical(engine_mod_token("STY", "ANY", "comet"), "STY") &&
   identical(engine_mod_token("*", "ANY", "comet"), "ACDEFGHIKLMNPQRSTVWY"))
ck("msfragger 任意残基端 [^ / ]$、特定残基端 n 前缀",
   identical(engine_mod_token("*", "PEP_N", "msfragger"), "[^") &&
   identical(engine_mod_token("*", "PEP_C", "msfragger"), "]$") &&
   identical(engine_mod_token("Q", "PEP_N", "msfragger"), "nQ"))
pl_msf <- c("variable_mod_05 = 15.9949 M 3")
lmm <- label_mod_masses(c("MeN|PEP_N|*|H(4)C(2)", "NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)"))
out_msf <- label_mod_lines(pl_msf, lmm, "^variable_mod_[0-9]",
                           "variable_mod_%02d = %.6f %s 3", engine = "msfragger")
ck("msfragger 注入行：[^ 任意残基端 + K 残基",
   identical(out_msf[2], "variable_mod_06 = 28.031300 [^ 3") &&
   identical(out_msf[3], "variable_mod_07 = 8.014199 K 3"))
pl_com <- c("variable_mod05 = 15.9949 M 0 3 -1 0 0 0.0")
out_com <- label_mod_lines(pl_com, lmm, "^variable_mod[0-9]",
                           "variable_mod%02d = %.6f %s 0 3 -1 0 0 0.0", engine = "comet")
ck("comet 注入行：裸 n 端 + K 残基",
   identical(out_com[2], "variable_mod06 = 28.031300 n 0 3 -1 0 0 0.0") &&
   identical(out_com[3], "variable_mod07 = 8.014199 K 0 3 -1 0 0 0.0"))
lmp <- label_mod_masses("MeP|PRO_N|K|C(2)H(4)")
out_p <- label_mod_lines(c("variable_mod_05 = 15.9949 M 3"), lmp, "^variable_mod_[0-9]",
                         "variable_mod_%02d = %.6f %s 3", engine = "msfragger")
ck("蛋白端 msfragger 按肽段端近似注入", identical(out_p[2], "variable_mod_06 = 28.031300 nK 3"))
out_l <- label_mod_lines(c("a = b"), label_mod_masses("K602[K]"),
                         "^variable_mod_[0-9]", "variable_mod_%02d = %.6f %s 3")
ck("无既有编号从 01 起 + 纯残基 token 与旧行为一致",
   identical(out_l[2], "variable_mod_01 = 8.014204 K 3"))

# ---- 6. sage：端部跳过、残基行照常 ----
lmx <- label_mod_masses(c("NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)", "MeN|PEP_N|*|H(4)C(2)"))
ent <- sage_label_mod_entries(lmx)
ck("sage 条目只含残基行", grepl('"K": \\[8\\.014199\\]', ent) &&
   !grepl("MeN", ent) && !grepl("28.031300", ent))

# ---- 7. append_label_mods / strip_label_mods 理论侧 ----
a1 <- append_label_mods("KASTYK", "", c("Phos|ANY|STY|H(1)O(3)P(1)"))
ck("ANY+STY 定位 3/4/5（带标签名）", identical(a1, "3,Phos[STY];4,Phos[STY];5,Phos[STY];"))
a2 <- append_label_mods("KASTYK", "", c("MeN|PEP_N|*|H(4)C(2)"))
ck("PEP_N 定位 0（AnyN-term 标签）", identical(a2, "0,MeN[AnyN-term];"))
a3 <- append_label_mods("KASTYK", "", c("MeC|PEP_C|*|H(4)C(2)"))
ck("PEP_C 定位肽长（AnyC-term 标签）", identical(a3, "6,MeC[AnyC-term];"))
a4 <- append_label_mods("KASTYK", "", c("K602[K]"))
ck("旧预设定位行为不变", identical(a4, "1,K602[K];6,K602[K];"))
a5 <- append_label_mods("AST", "2,Oxidation[M];", c("Phos|ANY|STY|H(1)O(3)P(1)"))
ck("已有修饰串续接", identical(a5, "2,Oxidation[M];2,Phos[STY];3,Phos[STY];"))
s1 <- strip_label_mods("0,MeN[AnyN-term];3,Phos[STY];2,Oxidation[M];",
                       c("MeN|PEP_N|*|H(4)C(2)", "Phos|ANY|STY|H(1)O(3)P(1)"))
ck("自定义行按别名集剥除", identical(s1, "2,Oxidation[M];"))

# ---- 7b. 质量表增广 + 离子构建器实调 ----
invisible(parse_label_spec(c("NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)",
                             "MeN|PEP_N|*|H(4)C(2)")))
cl_opt <- getOption("isopsv_custom_labels")
ck("parse_label_spec 注册自定义行 option", !is.null(cl_opt) && nrow(cl_opt) == 2L)
aa_df <- read.delim(file.path(probe_code_dir(), "AA_mass_dataframe.txt"), sep = "\t",
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aug <- augment_aa_mass_table(aa_df)
ck("增广表含裸名与带标签名", all(c("NeuCode8", "NeuCode8[K]", "MeN[AnyN-term]") %in% aug$AA))
ck("增广质量 = 组成式 Δ（8.014199/28.031300）",
   identical(sprintf("%.6f", aug$Mono_weight[aug$AA == "NeuCode8[K]"]), "8.014199") &&
   identical(sprintf("%.6f", aug$Mono_weight[aug$AA == "MeN[AnyN-term]"]), "28.031300"))
ck("增广幂等（重复调用不重复追加）",
   identical(nrow(augment_aa_mass_table(aug)), nrow(aug)))
res0 <- by_ions_position("PEPTIDEK", "", "Unlabel", aa_df, amino_acid_table)
res1 <- by_ions_position("PEPTIDEK", "8,NeuCode8[K];", "Unlabel", aa_df, amino_acid_table)
ck("离子构建器吃下自定义修饰（Δ=8.014199）",
   abs((res1$pep_mass - res0$pep_mass) - 8.014199) < 1e-4)
res2 <- by_ions_position("PEPTIDEK", "0,MeN[AnyN-term];", "Unlabel", aa_df, amino_acid_table)
ck("端部自定义修饰走 N 端质量路由（Δ=28.031300）",
   abs((res2$pep_mass - res0$pep_mass) - 28.031300) < 1e-4)
s2 <- strip_label_mods("2,K602[K];3,Oxidation[M];", "K602[K]")
ck("旧预设剥除行为不变", identical(s2, "3,Oxidation[M];"))

# ---- 8. parse_label_spec：行值原样进入 mods（去重）----
sp <- parse_label_spec(c("C13", "NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)",
                         "NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)"))
ck("自定义行进 mods 并去重", identical(sp$class1, "C13") &&
   identical(sp$mods, "NeuCode8|ANY|K|C(-6)13C(6)N(-2)15N(2)"))

# ---- 9. 静态断言：server.py 多行修复 + index.html 表单 ----
srv <- paste(readLines(file.path(dirname(cd), "web", "server.py"), warn = FALSE,
                       encoding = "UTF-8"), collapse = "\n")
ck("server.py：label_type 已移出 single 名单",
   !grepl('"fasta", "label_type"', srv, fixed = TRUE) &&
   grepl('"fasta", "enzyme"', srv, fixed = TRUE))
ck("server.py：逐行净化 lt_rows", grepl("lt_rows", srv, fixed = TRUE) &&
   grepl('replace("\\t", " ")', srv, fixed = TRUE))
idx <- paste(readLines(file.path(dirname(cd), "web", "index.html"), warn = FALSE,
                      encoding = "UTF-8"), collapse = "\n")
ck("index.html：4 格行控件（name/pos/aa/comp）",
   grepl('class="label-name"', idx, fixed = TRUE) &&
   grepl('class="label-pos"', idx, fixed = TRUE) &&
   grepl('class="label-aa-mode"', idx, fixed = TRUE) &&
   grepl('class="label-comp"', idx, fixed = TRUE))
ck("index.html：五位置选项与预设快填",
   grepl('["PRO_C","蛋白C端 Protein C-term"]', idx) &&
   grepl('["ANY","任意位置 Anywhere"]', idx) &&
   grepl("LABEL_PRESETS", idx, fixed = TRUE))
ck("index.html：序列化含竖线语法",
   grepl('r.name+"|"+r.pos+"|"+r.aa+"|"+r.comp', idx, fixed = TRUE))
ck("user_step2 注入调用带 engine 参数",
   grepl('engine = "msfragger"', paste(readLines(file.path(cd, "user_step2.R"),
       warn = FALSE), collapse = "\n"), fixed = TRUE))

cat("==== label_custom_rows_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
