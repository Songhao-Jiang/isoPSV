# =============================================================================
# boot.R — 核心库唯一加载入口
# boot.R - single entry point loading the core libraries 
#
# 六个驱动脚本（user_step2/stageA_classify/stageB_conflict_report/user_step3/
# extend_report_A/unlabeled_extract）此前在文件顶部裸用相对路径 source 核心库，
# 仅当启动 cwd 恰为 code/ 时才成立（此前依赖 user_run.R setwd 的副作用）。
# 本文件按 --file= 定位脚本目录、以绝对路径加载 compat/acs_core/paired_params/
# ReadMGF_script，替代全部裸相对 source。各脚本顶部：
#   sd0 <- dirname(normalizePath(sub("^--file=", "", grep("^--file=",
#         commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
#   if (!dir.exists(sd0)) sd0 <- getwd()
#   source(file.path(sd0, "boot.R"))
# Six driver scripts used to source the core libraries by bare relative path,
# which only works when the cwd happens to be code/ (a side effect of
# user_run.R's setwd). This file locates the script folder via --file= and
# loads compat/acs_core/paired_params/ReadMGF_script by absolute path.
# =============================================================================
sd <- dirname(tryCatch(
  normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
                winslash = "/"), error = function(e) ""))
if (is.na(sd) || !dir.exists(sd)) {
  f1 <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  sd <- if (!is.na(f1)) dirname(sub("^--file=", "", f1)) else getwd()
  if (!dir.exists(sd)) sd <- getwd()
}
cd <- if (file.exists(file.path(sd, "compat.R"))) sd else file.path(sd, "..", "code")
if (!file.exists(file.path(cd, "compat.R"))) {
  stop("boot.R 无法定位代码目录（compat.R 未找到）：", cd,
       "\nboot.R cannot locate the code directory (compat.R not found): ", cd)
}
source(file.path(cd, "compat.R"))
suppressMessages(source(file.path(cd, "acs_core.R")))
suppressMessages(source(file.path(cd, "paired_params.R")))
suppressMessages(source(file.path(cd, "ReadMGF_script.R")))
suppressMessages(source(file.path(cd, "chimeric_check.R")))
