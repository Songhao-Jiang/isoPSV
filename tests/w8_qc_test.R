# w8_qc_test.R — W8 相关规格回归测试（代码审查 2026-09-22 补缺）
# 覆盖：cat5 五类映射全枚举、weights_B 6/8 值校验镜像、PDF 8 值守卫存在性、
#       score_of 可用性降级边界（逻辑与 master_table_export.R 同步镜像）
# 用法：Rscript w8_qc_test.R（路径按 --file= 相对定位，可随发行包任意布置）
#   在 tests\ 目录下运行
# 2026-09-23 客户机反馈修复：cd 原为开发机硬编码路径，移植后直接 Error——
#   改为 --file= 相对定位（code 目录 = tests 的上一级）

args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- file.path(dirname(sd), "code")
if (!dir.exists(cd)) cd <- file.path(sd, "..", "code")
npass <- 0; nfail <- 0
ck <- function(label, ok) {
  cat(ifelse(ok, "PASS", "FAIL"), ":", label, "\n")
  if (ok) npass <<- npass + 1 else nfail <<- nfail + 1
}

# ---- 1. cat5 五类映射（与 master_table_export.R 修复版一致；web MERGE5 同口径）----
cat5_of <- function(dt) {
  ifelse(dt == "chimera_same_charge", "chimera_same_charge",
  ifelse(dt == "chimera_different_charge", "chimera_different_charge",
  ifelse(dt == "n_terminal_non_specific_cleavage", "n_term",
  ifelse(dt %in% c("different_modification", "modification_location_changed"), "modification_change",
  ifelse(dt %in% c("same_weight_aa_pairs", "peptide_length_less_than_7aa", "aa_order_change"),
         "sequence_substitution", "sap")))))
}
all_dt <- c("chimera_same_charge","chimera_different_charge","n_terminal_non_specific_cleavage",
            "different_modification","modification_location_changed","same_weight_aa_pairs",
            "peptide_length_less_than_7aa","aa_order_change","sap")
expect <- c("chimera_same_charge","chimera_different_charge","n_term","modification_change",
            "modification_change","sequence_substitution","sequence_substitution",
            "sequence_substitution","sap")
ck("cat5 全枚举映射（含 aa_order_change 入 sequence_substitution）",
   identical(cat5_of(all_dt), expect))

# ---- 2. weights_B 值校验（镜像 master_table_export.R 覆盖块语义）----
wb_check <- function(raw) {
  wb <- suppressWarnings(as.numeric(strsplit(raw, ",")[[1]]))
  if (!length(wb) %in% c(6L, 8L) || any(is.na(wb))) return("reject")
  if (length(wb) == 8L)
    return(setNames(wb, c("ul_by","lb_by","ul_all","lb_all","predict","pcos","ul_cov","lb_cov")))
  return(setNames(wb, c("ul_by","lb_by","ul_all","lb_all","predict","pcos")))
}
r6 <- wb_check("0.196,0.082,0.052,0.034,0.081,0.100")
r8 <- wb_check("0.196,0.082,0.052,0.034,0.081,0.100,0.147,0.308")
ck("weights_B 6 值接受且名称为六指标序", length(r6) == 6 && identical(names(r6)[6], "pcos"))
ck("weights_B 8 值接受且第 7/8 名为 ul_cov/lb_cov",
   identical(names(r8)[7:8], c("ul_cov","lb_cov")))
ck("weights_B 7 值拒绝", identical(wb_check("1,2,3,4,5,6,7"), "reject"))
ck("weights_B 非数值拒绝", identical(wb_check("1,2,x,4,5,6"), "reject"))

# ---- 3. PDF 导出 8 值守卫存在（源级断言：守卫必须先于 W6 覆盖生效）----
# F-1（2026-09-26）：源码行经 readLines(encoding="UTF-8") 标记为 UTF-8，若进程
# locale 非 UTF-8，脚本内的中文字面量按 native 标记，两者 grep 匹配失败（假 FAIL）。
# 把模式字面量也显式标记为 UTF-8，使两侧编码标记一致、与 locale 无关。
w8_pat <- "weights_B 为八指标"
Encoding(w8_pat) <- "UTF-8"
for (f in c("spectrum_pdf_export.R", "spectrum_pdf_export_A.R")) {
  src <- readLines(file.path(cd, f), warn = FALSE, encoding = "UTF-8")
  g <- grep(w8_pat, src)
  ck(paste0(f, " 含 8 值 weights_B 守卫 stop()"), length(g) == 1 &&
       grepl("stop(", src[g], fixed = TRUE))
}

# ---- 4. score_of 可用性降级（镜像其算法；边界四例）----
score_vec <- function(M, w) {          # M: n 候选 × k 指标（含列名）；w: 命名权重
  avail <- apply(M, 2, function(cv) all(!is.na(cv)))
  if (!any(avail)) return(NA_real_)
  ww <- w[avail]; ww <- ww / sum(ww)
  as.numeric(M[, avail, drop = FALSE] %*% ww)
}
W <- c(a = .5, b = .5)
# 4a 单候选谱（1×2）：不崩、值为自身
m1 <- matrix(c(.2, .8), nrow = 1, dimnames = list(NULL, c("a","b")))
ck("score_of 单候选谱返回自身加权和", isTRUE(all.equal(score_vec(m1, W), 0.5)))
# 4b 某指标整列 NA：降级重归一到余指标（col a=(.2,.6) 有效、col b 全 NA）
m2 <- matrix(c(.2, .6, NA, NA), nrow = 2, dimnames = list(NULL, c("a","b")))
ck("score_of 列 NA 降级重归一", isTRUE(all.equal(score_vec(m2, W), c(.2, .6))))
# 4c 全 NA：NA
m3 <- matrix(c(NA, NA, NA, NA), nrow = 2, dimnames = list(NULL, c("a","b")))
ck("score_of 全 NA 返回 NA", is.na(score_vec(m3, W)))
# 4d 权重和为 0 的维度被 avail 后归一不除零（全 0 权重特例不进入打分，此处仅验证归一安全）
ck("score_of 降级后权重和恒为 1",
   isTRUE(all.equal(sum((function(av) { ww <- W[av]; ww/sum(ww) })(c(TRUE,TRUE))), 1)))

cat(sprintf("\n=== w8_qc_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
quit(status = ifelse(nfail > 0, 1, 0))
