# labels_rows_test.R — 多标记行式口径单元测试（v1.0，2026-09-28 用户口径）
# 覆盖：parse_label_spec（空/旧单值 legacy/新多行向量）、label_tau_of、
# append_label_mods（N 端与残基位点展开）、网页与 server 的静态断言。
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

# ---- 1. parse_label_spec ----
sp0 <- parse_label_spec(character(0))
ck("空 label_type → legacy unlabeled", identical(sp0$legacy, "unlabeled") && is.null(sp0$class1))
sp1 <- parse_label_spec("N15")
ck("单值 N15 → legacy N15", identical(sp1$legacy, "N15") && identical(sp1$class1, "N15"))
sp2 <- parse_label_spec("silac")
ck("单值 silac → legacy SILAC", identical(sp2$legacy, "SILAC") && is.null(sp2$class1))
sp3 <- parse_label_spec(c("N15", "Dimethyl[AnyN-term]"))
ck("多行 N15+二甲基 → class1=N15、mods=1", identical(sp3$class1, "N15") &&
   identical(sp3$mods, "Dimethyl[AnyN-term]"))
sp4 <- parse_label_spec(c("Label_13C(6)[K]", "Label_13C(6)15N(4)[R]"))
ck("多行双修饰 → class1=NULL、mods=2", is.null(sp4$class1) && length(sp4$mods) == 2L)

# ---- 2. label_tau_of ----
ck("legacy SILAC → 0.0183", identical(label_tau_of(sp2), 0.0183))
ck("legacy N15 → 0.0083", identical(label_tau_of(sp1), 0.0083))
ck("纯修饰标记 → 0.0183", identical(label_tau_of(sp4), 0.0183))
ck("类1+修饰 → 0.0083", identical(label_tau_of(sp3), 0.0083))

# ---- 3. append_label_mods ----
m1 <- append_label_mods("AAAKK", "", c("Dimethyl[AnyN-term]"))
ck("N 端二甲基 → 0,Dimethyl[AnyN-term];", identical(m1, "0,Dimethyl[AnyN-term];"))
m2 <- append_label_mods("AKAK", "", c("Label_13C(6)[K]"))
ck("K 位点展开 → 2 与 4 位", identical(m2, "2,Label_13C(6)[K];4,Label_13C(6)[K];"))
m3 <- append_label_mods("AR", "3,Oxidation[M];", c("Label_13C(6)15N(4)[R]"))
ck("叠加既有修饰", identical(m3, "3,Oxidation[M];2,Label_13C(6)15N(4)[R];"))

# ---- 4. 静态断言：网页与 server ----
html <- paste(readLines(file.path(cd, "..", "web", "index.html"), warn = FALSE), collapse = "\n")
ck("网页类1下拉存在", grepl('id="label_class1"', html))
ck("网页标记行模块存在", grepl("LABEL_PRESETS", html) && grepl("function addLabelRow", html))
ck("网页 label_type 合成多行", grepl("labelTypeValue()", html))
srv <- paste(readLines(file.path(cd, "..", "web", "server.py"), warn = FALSE), collapse = "\n")
ck("server label_type 多行写入", grepl('"label_type\\\\t" \\+ "\\\\n".join', srv))

# ---- 5. v1.0 无标记通道：剥标记 + 0 补齐余弦 ----
sm <- strip_label_mods("0,Acetyl[AnyN-term];2,Label_13C(6)[K];3,Oxidation[M];", c("Label_13C(6)[K]"))
ck("剥标记：去除 Label_13C(6)[K]", identical(sm, "0,Acetyl[AnyN-term];3,Oxidation[M];"))
sm2 <- strip_label_mods("0,Acetyl[AnyN-term];2,K602[K];3,Oxidation[M];", c("K602[K]"))
ck("剥标记：去除 K602[K]", identical(sm2, "0,Acetyl[AnyN-term];3,Oxidation[M];"))
sm3 <- strip_label_mods("2,Label_13C(6)[K];", c("Label_13C(6)[K]"))
ck("剥标记：全部为标记 → 空串", identical(sm3, ""))
sm4 <- strip_label_mods("3,Oxidation[M];", c("Label_13C(6)[K]"))
ck("无标记不动", identical(sm4, "3,Oxidation[M];"))
v_pad <- c(1, NA, 2); pv_pad <- c(1, 2, 0)
ck("0 补齐余弦（缺失位点补 0）", is.finite(cosine(v_pad[!is.na(v_pad)], pv_pad[!is.na(v_pad)])) &&
   abs(cosine(c(1,0,2), c(1,2,0)) - 0.2) < 1e-6)

cat("汇总：", npass, " PASS / ", nfail, " FAIL\n", sep = "")
if (nfail > 0L) quit(status = 1L)
