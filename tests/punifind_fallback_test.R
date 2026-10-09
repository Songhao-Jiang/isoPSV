# punifind_fallback_test.R — run_punifind_cli 回退复用回归（v1.0 排队⑤-4）
# 覆盖：模块文件存在时委托 pf_run（报错形态 = pf_probe 的多组件清单，证明走
# 了模块参数面）；模块文件缺失时走 legacy（报错形态 = 首组件 exe 缺失）。
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- if (file.exists(file.path(sd, "..", "code", "compat.R")))
  file.path(sd, "..", "code") else file.path(sd, "..", "..", "code")
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))
suppressMessages(source(file.path(cd, "punifind_gpu_compat.R"), encoding = "UTF-8"))

npass <- 0L; nfail <- 0L
ck <- function(label, cond) {
  if (isTRUE(cond)) { npass <<- npass + 1L; cat("PASS :", label, "\n") }
  else { nfail <<- nfail + 1L; cat("FAIL :", label, "\n") }
}

td <- tempfile("pfb_"); dir.create(td)
fake_dir <- file.path(td, "punifind_fake"); dir.create(fake_dir)
qry <- file.path(td, "qry"); mgf <- file.path(td, "mgf"); out <- file.path(td, "out")
dir.create(qry); dir.create(mgf)

# 1) 模块存在（随包分发）→ 委托 pf_run：pf_probe 报多组件缺失清单
assign("code_dir", cd, envir = .GlobalEnv)
on.exit(rm("code_dir", envir = .GlobalEnv), add = TRUE)
r1 <- tryCatch(run_punifind_cli("rescore", fake_dir, qry, mgf, out),
               error = function(e) conditionMessage(e))
ck("模块存在时委托 pf_run（报错含多组件清单而非 legacy 首组件）",
   grepl("ckpts/checkpoint_rank.pt", r1) && grepl("official_denovo_workflow.exe", r1))

# 2) 模块缺失 → legacy 路径（首组件 exe 缺失报错形态）
assign("code_dir", file.path(td, "no_module"), envir = .GlobalEnv)
r2 <- tryCatch(run_punifind_cli("rescore", fake_dir, qry, mgf, out),
               error = function(e) conditionMessage(e))
ck("模块缺失时走 legacy（首组件 exe 报错形态）",
   grepl("pUniFind 组件缺失：", r2) && grepl("official_score_workflow.exe", r2) &&
   !grepl("ckpts", r2))
unlink(td, recursive = TRUE)

cat(sprintf("=== punifind_fallback_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
if (nfail > 0L) quit(status = 1L)
