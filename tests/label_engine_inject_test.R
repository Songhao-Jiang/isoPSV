# label_engine_inject_test.R — v1.0 外部引擎自动搜库标记注入回归（2026-09-28）
# 覆盖：label_mod_masses（名称→残基/质量解析）、label_mod_lines（comet/msfragger
# 编号续接）、sage_label_mod_entries（JSON 条目聚合）、三个适配器质量→名称映射、
# user_step2 线程传递与侧车指纹的静态断言。
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

# ---- 1. label_mod_masses：名称 → (name, residue, mass) ----
lm0 <- label_mod_masses(character(0))
ck("空标签 → NULL", is.null(lm0))
lm1 <- label_mod_masses(c("K080[K]", "K602[K]", "R004[R]", "R040[R]"))
ck("4 标签解析出 4 行", is.data.frame(lm1) && nrow(lm1) == 4L)
ck("残基解析 K/K/R/R",
   identical(lm1$residue, c("K", "K", "R", "R")))
ck("质量解析 8.050220/8.014204/3.988142/4.025110",
   identical(sprintf("%.6f", lm1$mass), c("8.050220", "8.014204", "3.988142", "4.025110")))
lm2 <- label_mod_masses("R003[R]")
ck("R003 解析 2.991107", identical(sprintf("%.6f", lm2$mass), "2.991107") &&
   identical(lm2$residue, "R"))
lm3 <- label_mod_masses("NoSuch[Z]")
ck("未知标签 → 质量 NA（守卫会禁用注入）", isTRUE(is.na(lm3$mass)) &&
   identical(lm3$residue, "Z"))

# ---- 2. label_mod_lines：编号续接（comet/msfragger 共用） ----
pl_msf <- c("search_enzyme_name_1 = trypsin", "variable_mod_05 = 15.9949 M 3")
out_msf <- label_mod_lines(pl_msf, lm1, "^variable_mod_[0-9]",
                           "variable_mod_%02d = %.6f %s 3")
ck("msfragger：续接 06-09 共 4 行", length(out_msf) == 6L &&
   identical(out_msf[3:6],
             c("variable_mod_06 = 8.050220 K 3",
               "variable_mod_07 = 8.014204 K 3",
               "variable_mod_08 = 3.988142 R 3",
               "variable_mod_09 = 4.025110 R 3")))
pl_com <- c("num_output_lines = 5", "variable_mod05 = 15.9949 M 0 3 -1 0 0 0.0")
out_com <- label_mod_lines(pl_com, lm1, "^variable_mod[0-9]",
                           "variable_mod%02d = %.6f %s 0 3 -1 0 0 0.0")
ck("comet：续接 06 起且格式为 0 3 -1 0 0 0.0",
   identical(out_com[3], "variable_mod06 = 8.050220 K 0 3 -1 0 0 0.0") &&
   identical(out_com[6], "variable_mod09 = 4.025110 R 0 3 -1 0 0 0.0"))
out0 <- label_mod_lines(c("a = b"), lm2, "^variable_mod_[0-9]", "variable_mod_%02d = %.6f %s 3")
ck("无既有编号 → 从 01 开始", identical(out0[2], "variable_mod_01 = 2.991107 R 3"))
# 真实模板实弹：包内 comet/msfragger 模板编号续接不撞车
tpl_msf <- readLines(file.path(cd, "engine_templates", "msfragger_con.params"), warn = FALSE)
tpl_com <- readLines(file.path(cd, "engine_templates", "comet_con.params"), warn = FALSE)
t2_msf <- label_mod_lines(tpl_msf, lm1, "^variable_mod_[0-9]", "variable_mod_%02d = %.6f %s 3")
t2_com <- label_mod_lines(tpl_com, lm1, "^variable_mod[0-9]", "variable_mod%02d = %.6f %s 0 3 -1 0 0 0.0")
n_old_msf <- grep("^variable_mod_[0-9]", tpl_msf)
n_old_com <- grep("^variable_mod[0-9]", tpl_com)
ck("真实 msfragger 模板：编号从 05 后续接",
   max(as.integer(sub("^variable_mod_([0-9]+).*", "\\1", tpl_msf[n_old_msf]))) == 5L &&
   grepl("^variable_mod_06 = 8.050220 K 3$", tail(t2_msf, 4)[1]))
ck("真实 comet 模板：编号从 05 后续接",
   max(as.integer(sub("^variable_mod([0-9]+).*", "\\1", tpl_com[n_old_com]))) == 5L &&
   grepl("^variable_mod06 = 8.050220 K 0 3 -1 0 0 0.0$", tail(t2_com, 4)[1]))

# ---- 3. sage_label_mod_entries：按残基聚合 JSON 条目 ----
sg <- sage_label_mod_entries(lm1)
ck("sage 条目：K/R 两组聚合",
   identical(sg, '      "K": [8.050220, 8.014204],\n      "R": [3.988142, 4.025110]'))

# ---- 4. 静态断言：适配器质量→名称映射 ----
for (f in c("search_adapter_comet.R", "search_adapter_msfragger.R", "search_adapter_sage.R")) {
  txt <- paste(readLines(file.path(cd, f), warn = FALSE), collapse = "\n")
  ck(paste(f, "含 4 个标签质量映射"),
     all(vapply(c('"8.050220" = "K080"', '"8.014204" = "K602"',
                  '"3.988142" = "R004"', '"4.025110" = "R040"'),
                function(pat) grepl(pat, txt, fixed = TRUE), logical(1))))
}

# ---- 5. 静态断言：user_step2 线程传递与侧车指纹 ----
s2 <- paste(readLines(file.path(cd, "user_step2.R"), warn = FALSE), collapse = "\n")
ck("engine_auto_run 带 label_mass 形参",
   grepl("label_mass = NULL, cand_scope = \"filtered\") {", s2, fixed = TRUE))
ck("侧车含 label= 指纹", grepl('paste0\\("label="', s2))
ck("调用点传入 label_mass = label_mass_auto", grepl("label_mass = label_mass_auto", s2, fixed = TRUE))
ck("线程计算：unlabeled_channel=0 时解析 label_mod_masses",
   grepl("label_mass_auto <- if (!identical(ul_ch_auto, \"1\"))", s2, fixed = TRUE) &&
   grepl("parse_label_spec(base_cfg$label_type)$mods", s2, fixed = TRUE))
ck("NA 质量守卫禁用注入", grepl("anyNA(label_mass$mass)", s2, fixed = TRUE))

cat("==== label_engine_inject_test:", npass, "pass /", nfail, "fail ====\n")
if (nfail > 0) quit(status = 1)
