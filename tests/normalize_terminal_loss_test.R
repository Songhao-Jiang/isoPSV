# =============================================================================
# normalize_terminal_loss_test.R - 终端丢失规范化单元测试（四类典型用例）
# normalize_terminal_loss_test.R - unit tests for terminal-loss normalization
# 运行：在本目录 Rscript normalize_terminal_loss_test.R（预期全部 PASS）
# =============================================================================
args0 <- commandArgs(FALSE)
sd0 <- dirname(sub("^--file=", "", args0[grepl("--file=", args0)][1]))
code <- file.path("..", "code")
if (!dir.exists(code)) code <- file.path(sd0, "..", "code")
source(file.path(code, "compat.R"))
source(file.path(code, "ReadMGF_script.R"))

chk <- function(name, got, exp_s, exp_m) {
  ok <- identical(got$s, exp_s) && identical(got$m, exp_m)
  cat(sprintf("[%s] %s  seq=%s mod=%s\n", ifelse(ok, "PASS", "FAIL"), name, got$s, got$m))
  if (!ok) cat(sprintf("  期望 seq=%s mod=%s\n", exp_s, exp_m))
  ok
}
ok <- TRUE
# 1) Met-loss+Acetyl：序列去首位 M、新增 0,Acetyl、其余位置前移、旧位置1残基修饰丢弃
ok <- chk("Met-loss+Acetyl+位移+丢弃",
          normalize_terminal_loss("MAMLSTTTIGFPR", "0,Met-loss+Acetyl[ProteinN-termM];4,Oxidation[M];2,Carbamidomethyl[C];"),
          "AMLSTTTIGFPR", "3,Oxidation[M];1,Carbamidomethyl[C];0,Acetyl[AnyN-term];") && ok  # 源串顺序保留
# 2) 纯 Met-loss：无新增修饰
ok <- chk("纯 Met-loss",
          normalize_terminal_loss("MDPELQQLATVFQ", "1,Met-loss[ProteinN-termM];"),
          "DPELQQLATVFQ", "") && ok
# 3) Arg-loss：序列去尾 R、无修饰残留、其余位置不动
ok <- chk("Arg-loss",
          normalize_terminal_loss("VLLAPEDEEGGDSGDYPLKR", "20,Arg-loss[AnyC-termR];12,Oxidation[M];"),
          "VLLAPEDEEGGDSGDYPLK", "12,Oxidation[M];") && ok
# 4) 位置 1 末端类归 0（无任何丢失，仅归一）
ok <- chk("末端类 1->0",
          normalize_terminal_loss("AAEPEAVATEADDFTF", "1,Acetyl[AnyN-term];3,Oxidation[M];"),
          "AAEPEAVATEADDFTF", "0,Acetyl[AnyN-term];3,Oxidation[M];") && ok
# 5) 无修饰/无丢失：原样
ok <- chk("无修饰", normalize_terminal_loss("PEPTIDEK", ""), "PEPTIDEK", "") && ok
ok <- chk("普通修饰不动",
          normalize_terminal_loss("PEPTIDEKM", "8,Oxidation[M];"),
          "PEPTIDEKM", "8,Oxidation[M];") && ok
cat(ifelse(ok, "全部 PASS\n", "存在 FAIL\n"))
if (!ok) quit(status = 1L)   # P0-6（2026-09-26）：失败必须非零退出，否则 run_all 假 PASS
