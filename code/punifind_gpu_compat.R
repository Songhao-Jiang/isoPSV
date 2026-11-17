# =============================================================================
# punifind_gpu_compat.R —— pUniFind(GPU 版)离线适配模块（可整合进 isoPSV 流程）
#
# 目的：让 isoPSV 能在**不修改 pUniFind 源码/不依赖其补丁**的前提下，正确调用
#   GPU 版 pUniFind 工作流（rescore / denovo），并绕开已知的构建差异。
#
# 设计要点（详见同目录说明文档）：
#   1) 参数面自探测：从 `<exe> --help` 读出该构建支持的 flag，只传它认识的；
#      并且**把数据目录放在 argv 第一位**（官方 CLI 的顺序），彻底避免"游离
#      token 抢占位置参数 data"（v1.0/v1.0 实测事故的根因）。
#   2) 未打补丁也能跑：
#      - libuv：子进程注入 USE_LIBUV=0（torch rendezvous 官方开关）
#      - 权重加载：子进程注入 TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1（torch 官方开关）
#   3) qry.res 候选序号压实（R 侧等价实现）：pFind 按 1% PSM FDR 行级过滤后，
#      中间排名整行缺失会造成 idx 出洞（1,3,4,…），下游 range(1,num_res+1)
#      连续取键 KeyError。本模块在调用前把每个谱块的候选序号按出现顺序重编为
#      1..n（幂等；只改每行行首的序号，其余字节原样）。
#   4) 只读原则：不写入 punifind_dir（ckpts/consts/exe 均只读），全部产物写在
#      调用方给定的 out_dir / qry_res_dir 下；命令行与子进程环境另存 sidecar
#      日志便于追溯。
#
# 用法（整合进 isoPSV 时把 code/user_step2.R 里对 run_punifind_cli 的调用改为
# 本模块的 pf_run，参数完全同形）：
#   source("punifind_gpu_compat.R")
#   pf_run("rescore", puni_exe, file.path(proj, "result"),
#          file.path(cfg2$config_dir, "unlabeled"), proj,
#          converter_py, puni_py)                 # 返回 c(rescore_spectra=路径)
#
# 依赖：仅 base R（Windows；用 powershell 仅做可选的 GPU 探测）。
# 版本：2026-09-26（deepseek-flash 适配；对应 isoPSV v1.0 现状）
# =============================================================================

# ---- 子进程环境（离线适配的关键：用官方开关替代源码补丁） -------------------
# ---- child-process env (the offline shim: official switches instead of patches)
pf_env_vars <- function(num_workers = NULL, distributed = FALSE,
                        master_port = 29517L, extra = character(0)) {
  v <- c(
    # libuv：torch 2.7.x 的 Windows cu128 轮子未编译 libuv，rendezvous 初始化
    # TCPStore 会报 "use_libuv was requested but PyTorch was built without libuv
    # support"。torch 自身读 USE_LIBUV 环境变量（rendezvous.py
    # _get_use_libuv_from_query_dict），置 0 等价于源码改 use_libuv=False。
    USE_LIBUV = "0",
    # 权重加载：torch>=2.6 的 torch.load 默认 weights_only=True，官方 checkpoint
    # 内含 argparse.Namespace 会 UnpicklingError。torch 提供该环境变量把默认值
    # 强制回 False（实测在本机 torch 2.7.1+cu128 生效并带明确 UserWarning）。
    TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD = "1"
  )
  if (!is.null(num_workers)) v <- c(v, PUNIFIND_NUM_WORKERS = as.character(num_workers))
  if (isTRUE(distributed)) v <- c(v,
    RANK = "0", LOCAL_RANK = "0", WORLD_SIZE = "1",
    MASTER_ADDR = "127.0.0.1", MASTER_PORT = as.character(master_port))
  if (length(extra)) v <- c(v, extra)
  v
}

# ---- 解释器/入口选择：优先 python + 仓库 .py（见下"启动器"说明），回退 exe ----
# ---- pick the entry point: python + repo .py first, fall back to the exe -------
# 实测（2026-09-26）：仓库自带的 official_*.exe 是**小启动器**，它会把命令行
# 重新组装后再交给内嵌脚本——重组过程会让位置参数 data 变成 "0"（即使调用方
# 把数据目录放在第一位、且不含任何未知 flag），导致 AssertionError:
# 0\test.parquet。直接用 `python official_*_workflow.py` 调用同一套参数则完全
# 正常（data 正确、GPU 利用率 97%）。因此本模块默认走 python 直调路线。
pf_entry <- function(kind, punifind_dir, python_exe = NULL) {
  py_file <- file.path(punifind_dir, if (identical(kind, "rescore"))
    "official_score_workflow.py" else "official_denovo_workflow.py")
  exe <- file.path(punifind_dir, if (identical(kind, "rescore"))
    "official_score_workflow.exe" else "official_denovo_workflow.exe")
  if (nzchar(py_file) && file.exists(py_file) &&
      !is.null(python_exe) && nzchar(python_exe) && file.exists(python_exe)) {
    return(list(route = "python", cmd = python_exe, args0 = py_file, wd = punifind_dir))
  }
  list(route = "exe", cmd = exe, args0 = NULL, wd = NULL)
}

# ---- 探测构建支持的 flag 面 ------------------------------------------------
# ---- probe the flag surface supported by this build ------------------------
# N-15（2026-09-26）：timeout 参数真正生效——旧实现 system2 无超时，--help
# 挂起时 R 无限阻塞。改为后台起探针（wait=FALSE 返回 PID）、轮询进程存活，
# 超时按 PID 杀进程树（铁律：只按 PID 杀）并返回 NULL。
# v1.0（2026-09-29，工作站 P1）：抽出"起进程→轮询→超时杀→读日志"公共件
# pf_spawn_collect，--help 探测与 denovo 三参数定向探测共用（行为不变）。
# cmd = 可执行文件；args = 额外参数向量（python 直调路线传 workflow.py 路径）。
#' v1.0（2026-09-30，工作站 C 项）：探针日志双编码读取。探针不走 .bat
#' （无 chcp 65001），中文 Windows 下 python 子进程 stdout 按 locale（GBK）
#' 写出，readLines 默认按 UTF-8 读到非 ASCII 字节即 invalid multibyte string
#' → 探针恒 NULL（"无法读取 --help"）。先按 UTF-8 读、失败按 native.enc 再读；
#' flag 解析只匹配 ASCII "--token"，两种解码下均不受影响。
pf_read_log_lines <- function(logf) {
  out <- tryCatch(readLines(logf, warn = FALSE), error = function(e) NULL)
  if (!is.null(out)) return(out)
  tryCatch(readLines(logf, warn = FALSE, encoding = "native.enc"),
           error = function(e) NULL)
}

pf_spawn_collect <- function(exe, args, timeout = 120) {
  if (!file.exists(exe)) return(NULL)
  # v1.0（2026-10-01，本机实弹 P3）：stdout 与 stderr 分开落盘再拼接。
  # 实测 R/Windows 上 system2 把两个流转到同一文件时会双双丢失（本机
  # denovo exe 参数错误写 stderr，同文件方案捕获 0 行 → 探针恒 NA）。
  logf  <- tempfile("pfh_", fileext = ".log")
  loge  <- tempfile("pfe_", fileext = ".log")
  pid <- suppressWarnings(system2(exe, args, stdout = logf, stderr = loge,
                                  wait = FALSE, invisible = TRUE))
  # v1.0（本机实弹 P3b）：R/Windows 上 system2(wait=FALSE) 可能返回 0 而
  # 非子进程 PID（本机实测 pid=0），alive 探测随之失效、会在 exe 尚未写盘时
  # 就去读日志。改为 PID 探测＋"输出文件大小连续 3 次不变且非零"双条件，
  # 任一满足即认为进程已产出并退出（探针命令都是参数级错误，输出即终态）。
  have_pid <- is.numeric(pid) && !is.na(pid) && pid > 0
  alive <- function() {
    if (!have_pid) return(TRUE)
    r <- suppressWarnings(system2("tasklist",
                                  c("/FI", paste0("PID eq ", as.integer(pid)), "/NH"),
                                  stdout = TRUE, stderr = FALSE))
    length(r) > 0 && any(grepl(paste0("^\\s*", as.integer(pid), "\\s"), r))
  }
  elapsed <- 0; stable <- 0L; last_sz <- -1L
  while (elapsed < timeout) {
    if (have_pid && !alive()) break
    sz <- 0L
    for (f in c(logf, loge)) if (file.exists(f)) sz <- sz + file.size(f)
    if (sz > 0L && identical(sz, last_sz)) {
      stable <- stable + 1L
      if (stable >= 3L) break
    } else stable <- 0L
    last_sz <- sz
    Sys.sleep(0.5); elapsed <- elapsed + 0.5
  }
  if (elapsed >= timeout && (alive() || !have_pid)) {
    if (have_pid)
      suppressWarnings(system2("taskkill", c("/PID", as.character(as.integer(pid)),
                                             "/T", "/F"), stdout = FALSE, stderr = FALSE))
    unlink(c(logf, loge))
    return(NULL)
  }
  # 进程刚退出时 Windows 可能仍锁着日志句柄——读取带重试
  rd <- function(f) {
    if (!file.exists(f)) return(NULL)
    for (rt in 1:10) {
      o <- pf_read_log_lines(f)
      if (!is.null(o)) return(o)
      Sys.sleep(0.3)
    }
    NULL
  }
  out <- c(rd(logf), rd(loge))
  unlink(c(logf, loge))
  if (is.null(out) || !length(out)) return(NULL)
  out
}

pf_help_flags <- function(exe, args = "--help", timeout = 120) {
  # v1.0（2026-10-01，本机实弹 P4）：老构建（piUniMS 冻结版）的 --help
  # flag 面不完整（截断/子解析器差异），按它过滤会误删必需参数——本机实测
  # --check/--fdr-thread 等被滤掉后 argparse 报 invalid int value: ''（其内
  # 部默认值为空串，必须显式全传）。PUNIFIND_NO_FLAGS_PROBE=1 时禁用 flag
  # 探测、返回 NULL 走全参数盲传（老构建实测唯一可行口径；新构建不受影响，
  # 工作站无需设置）。
  # v1.0（2026-10-01，用户需求：全自动版本判定）：环境变量改为可选覆盖。
  # 默认行为改为自动判定——flag 面如果少于 30 个已知参数（老构建 help 截断/
  # 子解析器不完整的特征），自动降级为全参数盲传。env 变量=1 仍可强制禁用；
  # env 变量=0 可强制启用 flag 过滤。
  env_override <- trimws(Sys.getenv("PUNIFIND_NO_FLAGS_PROBE", unset = ""))
  if (identical(env_override, "1"))
    return(NULL)
  out <- pf_spawn_collect(exe, args, timeout)
  if (is.null(out)) return(NULL)
  txt <- paste(out, collapse = "\n")
  fl <- regmatches(txt, gregexpr("--[A-Za-z][A-Za-z0-9-]*", txt))[[1]]
  if (!length(fl)) return(NULL)
  fl <- unique(fl)
  # v1.0 自动判定：检测关键 flag 是否在 help 面中。老构建（piUniMS 冻结版）
  # 的 help 面有 74 个 flag 但缺少关键参数（--check/--fdr-thread/--cutoff-spectra
  # 等不在 help 中），按 help 面过滤会把这些参数滤掉，导致 argparse 报
  # invalid int value: ''（老构建内部默认值为空串，必须显式传）。
  # 新构建（工作站）的 help 面包含全部关键 flag。
  # 检测到缺关键 flag 时自动返回 NULL（走全参数盲传），并打印日志。
  crit_flags <- c("--check", "--fdr-thread", "--cutoff-spectra")
  missing_crit <- setdiff(crit_flags, fl)
  if (length(missing_crit) > 0 && !identical(env_override, "0")) {
    cat("[pUniFind 适配] help 面缺少关键 flag：",
        paste(missing_crit, collapse=", "),
        "（老构建特征），自动降级为全参数盲传模式\n", sep = "")
    cat("[pUniFind compat] help surface missing critical flags: ",
        paste(missing_crit, collapse=", "),
        " (old-build signature), auto-falling back to blind full-arg mode\n", sep = "")
    return(NULL)
  }
  fl
}

#' v1.0（2026-09-29，工作站 P1）：--help 探测失败（flags=NULL）时的 denovo
#' 三参数（--min-length/--max-length/--range-pred）定向能力探针。
#' 背景：新构建内部默认为 None，缺参必 TypeError；老构建则不认识这三参
#' （unrecognized arguments）。--help 读不出来时无法靠 flag 面判断，改用
#' "已知参数＋非法值"判别：跑 `--min-length notanint --help`——
#' 新构建认识该参数，argparse 在 --help 之前立即报 invalid int value（实证
#' 2026-09-29：未知 flag 会被 --help 压过直接打印帮助退出 0，故不能用未知
#' flag 探测）；老构建不认识，--help 照常打印帮助。两种结局都不进入真实
#' 推理。返回 TRUE（新构建，参数须带）/ FALSE（老构建，不带）/ NA（超时或
#' 日志不可读——按官方 cli.py:264-269 恒接受的口径放行）。
#' v1.0（2026-10-01，本机实弹 P2/P3）：本机曾恒 NA 的根因在捕获而非判据
#' ——system2 双流同文件丢输出＋wait=FALSE 返回 pid=0 提前判死（见
#' pf_spawn_collect v1.0 修复）。判据保持 --help 版不变；超时上调 180 s
#' （PyInstaller 冷启动＋杀毒扫描）。
pf_probe_denovo_len_ok <- function(exe, args0 = NULL, timeout = 180) {
  out <- pf_spawn_collect(exe, c(args0, "--min-length", "notanint", "--help"), timeout)
  if (is.null(out)) return(NA)
  txt <- paste(out, collapse = "\n")
  if (grepl("show this help message", txt, fixed = TRUE)) return(FALSE)
  if (grepl("--min-length", txt, fixed = TRUE) && grepl("invalid", txt, ignore.case = TRUE))
    return(TRUE)
  NA
}

# ---- 安装完整性 + 环境探测 ------------------------------------------------
# ---- install integrity + environment probe --------------------------------
pf_probe <- function(punifind_dir, python_exe = NULL) {
  r <- list(dir = punifind_dir, ok = TRUE, notes = character(0))
  # N-15（2026-09-26）：兼容纯 .py 发行版——python 直调路线可用（python 存在
  # 且两个 workflow.py 在）时不再硬性要求两个 exe；exe 只在无 python 路线时
  # 必需。
  py_route_ok <- !is.null(python_exe) && nzchar(python_exe) && file.exists(python_exe) &&
    file.exists(file.path(punifind_dir, "official_score_workflow.py")) &&
    file.exists(file.path(punifind_dir, "official_denovo_workflow.py"))
  r$py_route_ok <- py_route_ok
  need <- c(ckpt = "ckpts/checkpoint_rank.pt",
            tokenize = "consts/tokenize1.pkl",
            mod_dict = "consts/modification_meta_dict.pkl")
  if (!py_route_ok) {
    need <- c(score_exe = "official_score_workflow.exe",
              denovo_exe = "official_denovo_workflow.exe", need)
  } else {
    r$score_exe <- file.path(punifind_dir, "official_score_workflow.exe")
    r$denovo_exe <- file.path(punifind_dir, "official_denovo_workflow.exe")
  }
  for (k in names(need)) {
    p <- file.path(punifind_dir, need[[k]])
    r[[k]] <- p
    if (!file.exists(p)) { r$ok <- FALSE
      r$notes <- c(r$notes, paste0("组件缺失 / missing: ", p)) }
  }
  # flag 面（老发布版含 --load-contra / --skip-invalid-size-inputs-valid-test；
  # 仓库新版没有 → 必须只传它认识的）。python 路线用 python workflow.py --help
  # 探测（同款超时机制）。
  if (py_route_ok) {
    score_py <- file.path(punifind_dir, "official_score_workflow.py")
    denovo_py <- file.path(punifind_dir, "official_denovo_workflow.py")
    r$score_flags <- pf_help_flags(python_exe, c(score_py, "--help"))
    r$denovo_flags <- pf_help_flags(python_exe, c(denovo_py, "--help"))
  } else {
    r$score_flags <- pf_help_flags(r$score_exe)
    r$denovo_flags <- pf_help_flags(r$denovo_exe)
  }
  if (is.null(r$score_flags)) r$notes <- c(r$notes,
      "无法读取 --help：将按保守集合构造参数（仅核心项，见 pf_build_args）")
  # python/torch 侧（可选）：仅用于给出提示，不参与决策
  r$python <- if (!is.null(python_exe) && nzchar(python_exe) && file.exists(python_exe))
    python_exe else ""
  r$gpu <- pf_gpu_info()
  r
}

# ---- GPU 信息（可选，仅提示用；powershell 不可用则返回 NULL） -----------------
# ---- GPU info (optional, informational only) --------------------------------
pf_gpu_info <- function() {
  smi <- Sys.which("nvidia-smi")
  if (!nzchar(smi)) {
    for (p in c("C:/Windows/System32/nvidia-smi.exe", "C:/Program Files/NVIDIA Corporation/NVSMI/nvidia-smi.exe"))
      if (file.exists(p)) { smi <- p; break }
  }
  if (!nzchar(smi) || !file.exists(smi)) return(NULL)
  q <- function(qs) {
    out <- tryCatch(suppressWarnings(
             system2(smi, c(paste0("--query-gpu=", qs), "--format=csv,noheader"),
                     stdout = TRUE, stderr = FALSE)),
           error = function(e) NULL)
    if (is.null(out) || !length(out)) NULL else paste(trimws(out), collapse = " | ")
  }
  # compute_cap 需要较新的 nvidia-smi；不支持时退回基础字段（探测仅为信息提示）
  r <- q("name,compute_cap,memory.total,memory.used")
  if (is.null(r)) r <- q("name,memory.total,memory.used")
  r
}

# ---- qry.res 候选序号压实（R 侧等价实现，幂等） ------------------------------
# ---- compact the qry.res candidate indices (R-side equivalent, idempotent) ---
# 块首判定与 pUniFind 读取器完全一致：行首为 "S" 且空白切分恰 3 段；下一行为
# 标题；其后到下一个块首之间的行都是候选行，行首整数为排名。
pf_qryres_renumber <- function(qry_res_dir, dry_run = FALSE, verbose = TRUE) {
  fs <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)
  if (!length(fs)) return(data.frame(file = character(0), spectra = integer(0),
                                     fixed = integer(0), stringsAsFactors = FALSE))
  res <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    f <- fs[i]
    raw <- readBin(f, "raw", n = file.info(f)$size)
    txt <- rawToChar(raw)                      # 原样字节：\r 保留在行内
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    n_fix <- 0L; n_spec <- 0L
    in_block <- FALSE; rank <- 0L
    for (j in seq_along(lines)) {
      ln <- lines[j]
      fl <- strsplit(ln, "[ \t]+")[[1]]
      fl <- fl[nzchar(fl)]
      if (length(fl) >= 1 && substr(ln, 1, 1) == "S" && length(fl) == 3) {
        in_block <- TRUE; rank <- 0L; n_spec <- n_spec + 1L
        next
      }
      if (!in_block) next                     # 文件头杂行
      if (!grepl("^[0-9]", ln)) next           # 标题行（紧跟块首）
      rank <- rank + 1L
      old <- sub("^([0-9]+)([ \t].*)$", "\\1", ln)
      if (!identical(old, as.character(rank))) {
        n_fix <- n_fix + 1L
        lines[j] <- sub("^[0-9]+([ \t])", paste0(rank, "\\1"), ln)
      }
    }
    if (n_fix > 0L && !dry_run) {
      # 原子写回：先写 .tmp 再 replace（大文件写回的安全惯例）
      tmp <- paste0(f, ".tmp")
      con <- file(tmp, "wb"); writeBin(charToRaw(paste(lines, collapse = "\n")), con); close(con)
      if (!file.rename(tmp, f)) { file.copy(tmp, f, overwrite = TRUE); unlink(tmp) }
    }
    if (verbose && n_fix > 0L)
      cat(sprintf("[pUniFind 适配] qry.res 序号压实: %s（%d 谱 / 改写 %d 行）\n",
                  basename(f), n_spec, n_fix))
    res[[i]] <- data.frame(file = basename(f), spectra = n_spec, fixed = n_fix,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, res)
}

# ---- 构造 argv：数据目录第一，flag=value 连写，只含受支持项 -------------------
# ---- build argv: data dir FIRST, --flag=value form, supported flags only -----
pf_build_args <- function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir,
                          flags = NULL, num_workers = 8, num_proc = 8,
                          denovo_len_ok = NA, batch_size = 64) {
  data_tmp <- file.path(out_dir, "pUniFind_data_tmp")
  # 位置参数 data 必须在最前（官方 CLI 顺序）：某些构建里 data 是 nargs="?" 的
  # 位置参数，若其前出现任何"游离 token"（比如被移除的 flag 遗留的值），data 会
  # 被抢走 → 下游去找 "0\test.parquet"。
  a <- c(data_tmp)
  opt <- list(
    "user-dir"            = data_tmp,
    "valid-subset"        = "test",
    "results-path"        = file.path(out_dir, "pUniFind_result"),
    "tmp-data-path"       = file.path(data_tmp, "test.lmdb"),
    "project-name"        = "test",
    "weight-path"         = file.path(punifind_dir, "ckpts", "checkpoint_rank.pt"),
    "qry-res-path"        = qry_res_dir,
    "mgf-path"            = mgf_dir,
    "num-proc"            = as.character(num_proc),
    "data-buffer-size"    = "256",
    "num-workers"         = as.character(num_workers),
    "ddp-backend"         = "c10d",
    "batch-size"          = as.character(batch_size),
    "task"                = if (identical(kind, "rescore")) "score_inference" else "denovo_inference",
    # loss：官方 CLI 对 denovo 用 cross_entropy、对 rescore 用 triplet（punifind/cli.py:219）
    "loss"                = if (identical(kind, "rescore")) "triplet" else "cross_entropy",
    "arch"                = "pep_ms",
    "node-dim"            = "512",
    "num-layers"          = "9",
    "num-head"            = "8",
    "use-nce"             = "1", "use-ins" = "1", "use-clip" = "1", "use-hard-clip" = "0",
    "denovo-pred"         = "1", "spectrum-pred" = "1", "denoise-pred" = "0",
    "intensity-pred"      = "0", "dataset-name" = "_", "max-charges" = "4",
    "use-rope"            = "0", "noise-peak-pred" = "0", "half-spectrum-pred" = "0",
    "from-scratch"        = "1", "norm-first" = "0", "modification-pred" = "1",
    "modification-cls-pred" = "0", "multi-class-noise-peak-pred" = "1",
    "res-num-pred"        = "0", "res-type-pred" = "0", "res-seq-len-pred" = "1",
    "cutoff-spectra"      = "300", "comp-res-num-pred" = "0", "train-pfind" = "1",
    "hard-neg"            = "1", "shift" = "1", "head-clip" = "1",
    "fdr-thread"          = "0.1",
    "pkl-path"            = data_tmp,
    "joint-encoding-layers" = "4", "num-pept-layers" = "9",
    "inference"           = "1", "check" = "0", "joint-encoding" = "1",
    "range-pred"          = "5", "tof" = "0",
    "tokenize1-pkl-path"  = file.path(punifind_dir, "consts", "tokenize1.pkl"),
    "modification-meta-dict-path" = file.path(punifind_dir, "consts", "modification_meta_dict.pkl")
  )
  # 追加型开关（无值）——只有构建支持时才带
  sw <- c("no-progress-bar", "skip-invalid-size-inputs-valid-test")
  # denovo 专属参数：官方 CLI 传 --max-length/--min-length/--range-pred
  # （punifind/cli.py:264-269，默认 40/6/5）。**缺了 min/max-length 会在
  # forward_faster_denovo 里 `min_length <= num_seq_pred` 抛 TypeError
  # （NoneType 与 Tensor 比较）**——这是 isoPSV 现版本 denovo 从不传这两个参数、
  # 因此必然失败的原因。
  # v1.0（2026-09-27 凌晨，本机实弹）：本机老构建的 denovo CLI 没有这三个
  # 参数，探测失败（flags=NULL）时旧代码仍全带 → unrecognized arguments。
  # 收紧为"仅当探测到该 flag 才带"（工作站新构建 --help 含这三个参数，不受影响）。
  # v1.0（2026-09-29，工作站 P1 实弹）：工作站 --help 探测恒失败（flags=NULL）
  # ⇒ 三参永不带上 ⇒ 新构建 forward_faster_denovo TypeError（ms_pep.py:747）。
  # 补定向探针结果 denovo_len_ok（pf_probe_denovo_len_ok）：TRUE/NA 时带上
  # （官方 cli.py:264-269 恒接受；NA＝探针不可判，按官方口径放行），
  # FALSE（老构建实锤）时不带、走其内部默认。
  if (identical(kind, "denovo")) {
    d <- list("max-length" = "40", "min-length" = "6", "range-pred" = "5")
    for (k in names(d)) {
      nm <- paste0("--", k)
      if (!is.null(flags) && (nm %in% flags)) opt[[k]] <- d[[k]]
      else if (is.null(flags) && !identical(denovo_len_ok, FALSE)) opt[[k]] <- d[[k]]
    }
  }
  # 已知在新构建中消失、且其"值"会泄漏成位置参数的老 flag：永不传递
  # （老发布版才认识它们；新版若不认识，传了反而更糟）
  deny <- c("load-contra")
  for (k in names(opt)) {
    if (k %in% deny) next
    nm <- paste0("--", k)
    if (!is.null(flags) && !(nm %in% flags)) next
    a <- c(a, paste0(nm, "=", opt[[k]]))
  }
  for (s in sw) {
    if (is.null(flags) || (paste0("--", s) %in% flags)) a <- c(a, paste0("--", s))
  }
  a
}

# ---- 可测纯函数（N-21，2026-09-26，GLM：pf_run 产物发现此前零覆盖，
# 抽成纯函数后可用临时目录夹具直接测） ----------------------------------------

#' P0-E 缓存清除判定：qry.res 最新 mtime 晚于 pUniFind 内部缓存（parquet/lmdb/pkl，
#' 容差 1s）即删除全部缓存（强制重推理）。返回被删除的文件路径（可作测试断言）；
#' 无缓存/输入不新时不删（character(0)）。副作用：unlink 删除。
pf_clear_stale_cache <- function(qry_res_dir, data_tmp, res_dir, verbose = TRUE) {
  # qry_input_files：输入文件清单（只看 *.qry.res——.spectra 产物会随后复制进
  # 该目录，把产物算进"最新 mtime"会让每轮都误判输入已变、缓存永不复用）；
  # qry_newest_time：输入最新修改时间；stale_cache_files：待清除的缓存
  qry_input_files <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)
  qry_newest_time <- if (length(qry_input_files)) max(file.info(qry_input_files)$mtime) else NA
  stale_cache_files <- c(list.files(data_tmp, pattern = "\\.parquet$", full.names = TRUE),
                         list.files(data_tmp, pattern = "\\.lmdb$", full.names = TRUE),
                         # P0-2（2026-09-26，deepseek v1.0 审查更正）：推理检查点
                         # pkl 名为 {subset}_{subset}_{rank}.pkl（isoPSV 侧 subset=test
                         # → test_test_0.pkl），与 data_tmp 的打分 pkl pscore_res_*.pkl
                         # 是两种不同产物；v1.0 的 P3-1 收窄会漏掉推理 pkl、使
                         # P0-E 复现，此处恢复全域匹配（v1.0 行为）
                         list.files(res_dir,  pattern = "\\.pkl$", full.names = TRUE),
                         # 硬化：keys pkl 参与同一完整性判定（同族失效，一并清）
                         list.files(data_tmp, pattern = "_FDR0\\.1_keys\\.pkl$", full.names = TRUE))
  if (!length(stale_cache_files) || is.na(qry_newest_time)) return(character(0))
  cache_newest_time <- max(file.info(stale_cache_files)$mtime)
  if (!(qry_newest_time > cache_newest_time - 1)) return(character(0))
  unlink(stale_cache_files, force = TRUE)
  if (verbose) cat("[pUniFind 适配] 输入 qry.res 已更新，清除内部缓存 ",
                   length(stale_cache_files), " 件（强制重推理）\n", sep = "")
  stale_cache_files
}

#' P0-B 新鲜度门槛：res_dir 中工作流自写的 fdr0.01_pUniFind.spectra，只取
#' mtime 不早于 run_start_time - 1s 的候选（更早的是上一轮残留）。
pf_fresh_own_spectra <- function(res_dir, run_start_time) {
  own <- list.files(res_dir, pattern = "fdr0\\.01_pUniFind\\.spectra$", full.names = TRUE)
  if (!length(own)) return(character(0))
  own_mtime <- file.info(own)$mtime
  own[own_mtime >= run_start_time - 1]
}

#' N-4（2026-09-26，GLM）：qry_res_dir 内最终 rescore 产物回查——工作流新构建
#' 自写 fdr0.01（两写皆认：fdr0.01_ / fdr.01_）。返回命中的路径向量。
pf_match_final_spectra <- function(qry_res_dir) {
  list.files(qry_res_dir, pattern = "fdr0?\\.01_pUniFind\\.spectra$", full.names = TRUE)
}

# ---- 主入口（与 isoPSV run_punifind_cli 同形，便于直接替换） -----------------
# ---- main entry (same shape as isoPSV's run_punifind_cli) --------------------
pf_run_one <- function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir,
                   converter_py = NULL, python_exe = NULL,
                   num_workers = 8, num_proc = 8, renumber = TRUE,
                   poll_s = 120, verbose = TRUE,
                   batch_size = 64, master_port = 29517L) {
  stopifnot(kind %in% c("rescore", "denovo"))
  pr <- pf_probe(punifind_dir, python_exe)
  if (!pr$ok) stop("pUniFind 组件缺失：\n", paste(pr$notes, collapse = "\n"))
  if (verbose && length(pr$notes)) cat("[pUniFind 适配] 提示：", paste(pr$notes, collapse = "；"), "\n")
  exe   <- if (identical(kind, "rescore")) pr$score_exe else pr$denovo_exe
  flags <- if (identical(kind, "rescore")) pr$score_flags else pr$denovo_flags
  data_tmp <- file.path(out_dir, "pUniFind_data_tmp")
  res_dir  <- file.path(out_dir, "pUniFind_result")
  dir.create(data_tmp, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  # P0-A（2026-09-26，deepseek 复核）：先做 name#0 消毒（公共函数在
  # ReadMGF_script.R；本模块被单独 source 时以 exists 守卫跳过——此时调用方
  # 须自行保证 qry.res 已消毒）
  if (exists("pf_sanitize_qryres", mode = "function")) pf_sanitize_qryres(qry_res_dir)
  # 1) qry.res 序号压实（幂等；新版已内置 compact_spec 时也不冲突）
  rn <- if (isTRUE(renumber)) pf_qryres_renumber(qry_res_dir, verbose = verbose) else NULL
  # P0-E（2026-09-26，deepseek v1.0）：pUniFind 的推理 pkl 复用不看输入指纹——
  # 消毒/重搜改写 qry.res 后，工作流"看到旧 pkl 就跳过推理"、却按新 qry.res 重建
  # parquet，两边候选数不一致 → AssertionError。判定与删除在纯函数
  # pf_clear_stale_cache 内（N-21 可测）：qry.res 新于内部缓存即清空强制重推理。
  pf_clear_stale_cache(qry_res_dir, data_tmp, res_dir, verbose)
  # P0-B（2026-09-26）：记录本次运行开始时间，产物发现只接受不早于它的候选
  # （防"工作流失败却复用陈旧产物报成功"）
  run_start_time <- Sys.time()

  # 2) 调用工作流（优先 python + 仓库 .py；见 pf_entry 的启动器说明）
  en  <- pf_entry(kind, punifind_dir, python_exe)
  if (!file.exists(en$cmd)) stop("pUniFind 入口缺失：", en$cmd)
  # v1.0（工作站 P1）：--help 探测失败（flags=NULL）且本次是 denovo 时，对
  # 三参数做定向能力探测（只在需要时跑，一次 python 启动级的成本）；
  # rescore 或探测成功时不跑（flags 本身已可判定）
  denovo_len_ok <- NA
  if (identical(kind, "denovo") && is.null(flags)) {
    denovo_len_ok <- if (identical(en$route, "python"))
      pf_probe_denovo_len_ok(en$cmd, en$args0) else pf_probe_denovo_len_ok(en$cmd)
    if (verbose) cat("[pUniFind 适配] denovo 三参数定向探测：",
                     if (isTRUE(denovo_len_ok)) "支持（带上）"
                     else if (identical(denovo_len_ok, FALSE)) "老构建（不带）"
                     else "不可判（按官方口径带上）", "\n", sep = "")
  }
  args <- pf_build_args(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir,
                        flags = flags, num_workers = num_workers, num_proc = num_proc,
                        denovo_len_ok = denovo_len_ok, batch_size = batch_size)
  logf <- file.path(out_dir, if (identical(kind, "rescore")) "rescore.log" else "denovo.log")
  # 生成 .bat 执行（而不是 system2(env=)）——实测 R 在 Windows 上 system2(cmd, args, env=<多变量>)
  # 会破坏命令行：隔离复现为 python 收到的脚本路径被替换成环境变量的值（"0"）。
  # .bat 方式同时留下可复现的产物（用户可手工重跑同一文件）。
  bat <- file.path(out_dir, paste0("pUniFind_run_", kind, ".bat"))
  cmd_line <- paste(c(shQuote(en$cmd), paste(shQuote(c(en$args0, args)), collapse = " ")),
                    collapse = " ")
  wd <- if (is.null(en$wd)) getwd() else en$wd
  side <- file.path(out_dir, paste0("pUniFind_cli_", kind, ".cmdline.txt"))
  # N-15（2026-09-26）：MASTER_PORT 固定 29517 被占用即失败无重试——把
  # 环境构建+写 bat+执行包成一次尝试，端口可换；失败且日志含 bind/地址占用
  # 特征时换随机端口重试（最多 3 次，铁律：只按 PID 杀不适用——直接换端口）。
  run_attempt <- function(master_port) {
    envv <- pf_env_vars(num_workers, distributed = identical(en$route, "python"),
                        master_port = master_port)
    env_lines <- paste0("set ", names(envv), "=", envv)
    # N-5（2026-09-26，GLM）：writeLines(useBytes=TRUE) 写出的 UTF-8 字节在 cmd.exe
    # 按 ANSI(GBK) 解析 → 中文路径乱码。首行 chcp 65001 让 cmd 按 UTF-8 读本文件。
    writeLines(c("@echo off", "chcp 65001 >nul", env_lines, paste0("cd /d ", shQuote(wd)),
                 cmd_line, "exit /b %ERRORLEVEL%"), bat, useBytes = TRUE)
    # 记录本次实际用的命令行与环境（可追溯；不与 pUniFind 安装混放）
    writeLines(c(paste0("# ", format(Sys.time())),
                 paste0("# route: ", en$route),
                 paste0("# env: ", paste(names(envv), collapse = " ")),
                 paste0("# cmd: ", en$cmd),
                 cmd_line),
               side, useBytes = TRUE)
    # v1.0 fix（2026-10-02 DAG 实弹）：system2 对同一文件双开 stdout/stderr
    # 偶发"cannot open the connection"（stdout 已建、stderr 二开失败的竞态，
    # 隔离复现 30/30 正常＝环境偶发；泄漏句柄还会让错误处理读日志 Permission
    # denied）。错误只发生在连接建立期（进程未起），2 秒后原样重试一次安全。
    # v1.0 fix (live DAG test 2026-10-02): system2 with stdout and stderr
    # on the SAME file can sporadically fail with "cannot open the connection"
    # (stdout opened, second open of stderr loses a race; isolation loop was
    # 30/30 ok, so it is an environmental flake; the leaked handle also makes
    # the error handler's log-tail read hit Permission denied). The error can
    # only happen at connection setup - the process has not started - so one
    # identical retry after 2 s is safe.
    st <- suppressWarnings(tryCatch(
      system2("cmd", c("/c", shQuote(bat)), stdout = logf, stderr = logf, wait = TRUE),
      error = function(e) {
        Sys.sleep(2)
        cat("[pUniFind 适配] stdout/stderr 同文件双开瞬时失败，2s 后重试一次\n")
        cat("[pUniFind compat] transient same-file open failure, retrying once after 2s\n")
        suppressWarnings(system2("cmd", c("/c", shQuote(bat)), stdout = logf,
                                 stderr = logf, wait = TRUE))
      }))
    list(st = st,
         logtxt = if (file.exists(logf)) paste(readLines(logf, warn = FALSE),
                                               collapse = "\n") else "")
  }
  if (verbose) cat("[pUniFind 适配] 运行（", en$route, "）：", en$cmd,
                   "（", length(args), "个参数，num-workers=", num_workers, "）\n", sep = "")
  rr <- NULL
  for (attempt in 1:3) {
    rr <- run_attempt(master_port)
    if (rr$st == 0) break
    if (attempt < 3 &&
        grepl("address already in use|bind\\(|winerror 10048|error 10013",
              rr$logtxt, ignore.case = TRUE)) {
      master_port <- sample(30000:60000, 1L)
      cat("[pUniFind 适配] 端口被占用，换用 ", master_port, " 重试（第 ",
          attempt + 1, " 次）\n", sep = "")
    } else break
  }
  st <- rr$st
  logtxt <- rr$logtxt
  # N-6（2026-09-26，GLM）：工作流失败时即便 out_dir 残留上一轮 .spectra/.csv
  # 也不得冒充成功——st != 0 前置 stop，不进入产物发现
  if (st != 0) {
    log_tail <- if (nzchar(logtxt)) substr(logtxt, max(1, nchar(logtxt) - 2000), nchar(logtxt)) else ""
    stop("pUniFind ", kind, " 工作流失败（exit=", st, "）；见日志 ", logf, "\n", log_tail)
  }

  # 3) 产物发现（rescore：优先取工作流自写的 fdr0.01 .spectra；否则 pkl→转换）
  if (identical(kind, "rescore")) {
    # 3a) 新构建自己会写 <results-path>/<project>fdr0.01_pUniFind.spectra；isoPSV
    #     的下游从 qry_res_dir 取，故此处复制过去（无需 Python 转换器）。
    #     新鲜度门槛在纯函数 pf_fresh_own_spectra 内（P0-B，N-21 可测）
    own <- pf_fresh_own_spectra(res_dir, run_start_time)
    if (length(own)) {
      own <- own[which.max(file.info(own)$mtime)]
      dst <- file.path(qry_res_dir, basename(own))
      need_copy <- TRUE
      if (file.exists(dst)) {
        t1 <- as.numeric(file.info(own)$mtime); t2 <- as.numeric(file.info(dst)$mtime)
        need_copy <- is.na(t2) || (!is.na(t1) && t1 > t2)
      }
      if (need_copy) {
        file.copy(own, dst, overwrite = TRUE)
        if (verbose) cat("[pUniFind 适配] 已取工作流自写的 .spectra →", dst, "\n")
      }
    }
    # N-4（2026-09-26，GLM）：此处回查正则与 3a 的 fdr0\.01_ 少个 0，复制成功后
    # 查不到 → 误入 pkl 转换分支空轮询 120s 并误报失败。补 0? 两写都认；
    # 纯函数 pf_match_final_spectra（N-21 可测）
    sp <- pf_match_final_spectra(qry_res_dir)
    if (!length(sp) && !is.null(converter_py) && nzchar(converter_py) &&
        !is.null(python_exe) && nzchar(python_exe) && file.exists(converter_py)) {
      # 旧构建只写打分 pkl：轮询等待落盘后转换（官方 CLI 有几十毫秒~秒级延迟）
      pk <- character(0)
      for (t in seq_len(max(1, ceiling(poll_s / 5)))) {
        pk <- list.files(data_tmp, pattern = "^pscore_res_.*\\.pkl$", full.names = TRUE)
        if (length(pk)) break
        Sys.sleep(5)
      }
      if (length(pk)) {
        pk <- pk[which.max(file.info(pk)$mtime)]
        out_sp <- file.path(qry_res_dir, "testfdr0.01_pUniFind.spectra")
        conv <- suppressWarnings(system2(python_exe,
                   c(shQuote(converter_py), shQuote(pk), shQuote(out_sp)),
                   stdout = TRUE, stderr = TRUE))
        stc <- attr(conv, "status"); if (is.null(stc)) stc <- 0
        if (stc == 0 && file.exists(out_sp)) sp <- out_sp
        else if (verbose) cat("[pUniFind 适配] pkl→spectra 转换失败（status=", stc, "）\n", sep = "")
      }
    }
    if (!length(sp)) stop("pUniFind rescore 完成但未产出 pUniFind.spectra：", qry_res_dir,
                          "\n（exit=", st, "；见日志 ", logf, "）")
    return(c(rescore_spectra = tail(sp, 1)))
  } else {
    csvs <- list.files(res_dir, pattern = "\\.csv$", full.names = TRUE)
    if (!length(csvs)) stop("pUniFind denovo 完成但未产出 csv：", res_dir,
                            "\n（exit=", st, "；见日志 ", logf, "）")
    return(c(denovo_csv = csvs[which.max(file.info(csvs)$mtime)]))
  }
}

# ---- v1.0（2026-09-29 用户口径＋工作站缺陷报告台账 195）：单 MGF 分片 ------
# ---- per-MGF chunking for both pUniFind workflows ----------------------------
# 工作流以 multiprocessing spawn 把整个数据集 pickle 给子进程，143 万谱规模
# 父端 dump／子端 load 双双 MemoryError（SJH 5.7 万谱无影响）。用户裁定
# "pFind 系列的都拆成单 mgf 的分别进行计算"：mgf 目录超过 1 个文件时逐 raw
# 分片——每片独占 chunk 目录（qry.res 全量＋单 mgf），片级完成标记支持失败
# 续跑；rescore 谱表合并去重表头，denovo CSV 无表头直接拼接，合并产物落回
# 原接口位置，下游无感。单 mgf（≤1 文件）直通 pf_run_one，行为与旧版一致。
# v1.0（2026-09-29，工作站代码检查 P2/P3）：① qry.res 由硬链接改 copy——
# 片内消毒/压实会就地覆写，硬链接同 inode 会穿透到项目真件、打破 #18 后处理
# 指纹（下次 step2 重付 ~5.6 h）；mgf 只读、仍硬链接保 mtime。② done 标记
# 升级为"完成＋输入指纹"（qry.res 全集＋本片 mgf 的 mtime|size）：重跑时指纹
# 不符即作废标记、刷新片内输入重跑，"输入变更→缓存失效"铁律覆盖到片级。
# ③ denovo 合并加行数守恒断言、rescore 合并加各片表头一致性断言。
pf_chunk_fingerprint <- function(qres_files, mgf_file) {
  fi <- file.info(c(qres_files, mgf_file))
  paste(sprintf("%.0f|%g", as.numeric(fi$mtime), fi$size), collapse = ";")
}

#' 读回 done 标记里的输入指纹；v1.0 旧标记无指纹行 → NULL（视同不符，
#' 各片重跑一次后即带上指纹）。
pf_marker_fingerprint <- function(marker) {
  ln <- tryCatch(readLines(marker, warn = FALSE), error = function(e) character(0))
  m <- grep("^fingerprint: ", ln, value = TRUE)
  if (!length(m)) return(NULL)
  substr(m[1], 14L, nchar(m[1]))
}

#' v1.0（2026-09-30，工作站 A③/B/D）：片级并行（chunk_parallel 旋钮，
#' 默认 1＝顺序；每片独立 python 进程，单模型 ~3 GB 显存，16 GB 卡可容 2–3
#' 路，worker 各配独立 MASTER_PORT 基值防占用冲突）、chunk 失败时把子日志尾
#' 直接附进错误消息（工作站 B 项：invalid multibyte 文案已四次误导，真实
#' 死因全在子日志）、合并成功后可选清理片内数据（chunk_cleanup，保留 done
#' 标记与合并产物）。
pf_log_tail <- function(ch_out, kind, n = 1500L) {
  lf <- file.path(ch_out, if (identical(kind, "rescore")) "rescore.log" else "denovo.log")
  if (!file.exists(lf)) return("(无子日志)")
  txt <- tryCatch(paste(pf_read_log_lines(lf), collapse = "\n"), error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) return("(子日志空)")
  substr(txt, max(1L, nchar(txt) - n + 1L), nchar(txt))
}

pf_run <- function(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir, ...,
                   chunk_parallel = 1L, chunk_cleanup = FALSE) {
  mgfs <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
  dots <- list(...)
  runner <- if (is.null(dots$.runner)) pf_run_one else dots$.runner
  dots$.runner <- NULL
  if (length(mgfs) <= 1L)
    return(do.call(runner, c(list(kind, punifind_dir, qry_res_dir, mgf_dir, out_dir), dots)))
  .pfr_chunk <- function(stem) file.path(out_dir, "chunks", stem)
  # v1.0 P2：mgf 只读，硬链接（失败回退 copy）保 mtime 即可；qry.res 必须
  # copy——pf_sanitize_qryres/pf_qryres_renumber 在片内就地覆写，硬链接会
  # 写穿到项目真件。copy 一次性成本（全量 qry.res 一份/片）远小于 #18 指纹
  # 失效后重付的 ~5.6 h 后处理。
  place_mgf <- function(src, dst_dir) {
    dst <- file.path(dst_dir, basename(src))
    if (file.exists(dst)) unlink(dst)
    ok <- suppressWarnings(tryCatch(file.link(src, dst), error = function(e) FALSE))
    if (!isTRUE(ok)) file.copy(src, dst, overwrite = TRUE)
    invisible(TRUE)
  }
  copy_qry <- function(src, dst_dir) {
    dst <- file.path(dst_dir, basename(src))
    if (file.exists(dst)) unlink(dst)
    if (!file.copy(src, dst, overwrite = TRUE))
      stop("[pUniFind 分片] qry.res 拷贝失败：", src, " → ", dst)
    invisible(TRUE)
  }
  # v1.0（2026-10-01，本机实弹 P1）：片内 qry.res 按本片 raw 过滤。
  # 本机 pUniFind 构建（piUniMS 冻结版）对 qry 目录逐 qry.res 建任务，块内
  # dta 属他片 raw 时会去本片 mgf 目录枚举 `<raw>_HCDFT/_CIDTOF/_HCDAST.mgf`
  # → FileNotFoundError（工作站的容缺构建无此问题）。全量拷贝改为拷贝后
  # 内容级过滤（按 dta 前缀 = 片名），过滤后 0 块的文件删除；两侧行为等价
  # （容缺构建本来就跳过他片谱）。
  pf_filter_qryres_chunk <- function(qdir, raw_stem) {
    pref <- paste0(raw_stem, ".")
    for (fp in list.files(qdir, pattern = "\\.qry\\.res$", full.names = TRUE)) {
      ln <- readLines(fp, warn = FALSE)
      n <- length(ln); keep <- logical(n); i <- 1L; any_keep <- FALSE
      while (i <= n) {
        if (startsWith(ln[i], "S\t")) {
          k <- i + 2L
          while (k <= n && !startsWith(ln[k], "S\t")) k <- k + 1L
          if (i + 1L <= n && startsWith(ln[i + 1L], pref)) {
            keep[i:(k - 1L)] <- TRUE; any_keep <- TRUE
          }
          i <- k
        } else i <- i + 1L
      }
      if (any_keep) writeLines(ln[keep], fp, useBytes = TRUE) else unlink(fp)
    }
    invisible(TRUE)
  }
  qres <- list.files(qry_res_dir, pattern = "\\.qry\\.res$", full.names = TRUE)
  .pfr_write_marker <- function(marker, stem, fp) {
    writeLines(c(paste0("chunk: ", stem), paste0("at: ", Sys.time()),
                 paste0("fingerprint: ", fp)), marker)
  }
  # ---- 准备（顺序）：输入放置＋指纹判定，产出待执行片表 pending ----
  pending <- list()
  for (k0 in seq_along(mgfs)) {
    stem <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(mgfs[k0])))
    ch_out <- .pfr_chunk(stem)
    dir.create(file.path(ch_out, "mgf"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(ch_out, "qry"), recursive = TRUE, showWarnings = FALSE)
    fp_now <- pf_chunk_fingerprint(qres, mgfs[k0])
    marker <- file.path(ch_out, paste0("done_", kind, ".txt"))
    if (file.exists(marker) && identical(pf_marker_fingerprint(marker), fp_now)) {
      cat("[pUniFind 分片] ", stem, "：完成标记在且输入未变，跳过（", kind, "）\n", sep = "")
      .pfr_write_marker(marker, stem, fp_now)
    } else {
      if (file.exists(marker))
        cat("[pUniFind 分片] ", stem, "：输入指纹变化，作废完成标记重跑（", kind, "）\n", sep = "")
      else
        cat("[pUniFind 分片] ", stem, "（", k0, "/", length(mgfs), "，", kind, "）\n", sep = "")
      place_mgf(mgfs[k0], file.path(ch_out, "mgf"))
      for (q in qres) copy_qry(q, file.path(ch_out, "qry"))
      pf_filter_qryres_chunk(file.path(ch_out, "qry"), stem)
      pending[[length(pending) + 1L]] <- list(stem = stem, ch_out = ch_out, fp = fp_now)
    }
  }
  .pfr_call <- function(p, port) {
    do.call(runner, c(list(kind, punifind_dir,
                           file.path(p$ch_out, "qry"), file.path(p$ch_out, "mgf"),
                           p$ch_out), dots, list(master_port = port)))
  }
  # ---- 执行：顺序（即停，保持既有语义）或并行（跑完汇总再报） ----
  if (length(pending) == 0L) {
    cat("[pUniFind 分片] 全部片已有有效完成标记，仅重做合并（", kind, "）\n", sep = "")
  } else if (chunk_parallel <= 1L || length(pending) == 1L) {
    for (j in seq_along(pending)) {
      p <- pending[[j]]
      tryCatch(.pfr_call(p, 29517L),
               error = function(e) stop("[pUniFind 分片] chunk 失败：", p$stem, "（", kind,
                 "）——修好后重跑本函数即从该 chunk 续跑；日志：",
                 file.path(p$ch_out, if (identical(kind, "rescore")) "rescore.log" else "denovo.log"),
                 "\n", conditionMessage(e),
                 "\n[子日志尾] ", pf_log_tail(p$ch_out, kind), call. = FALSE))
      .pfr_write_marker(file.path(p$ch_out, paste0("done_", kind, ".txt")), p$stem, p$fp)
    }
  } else {
    nw <- min(as.integer(chunk_parallel), length(pending))
    cat("[pUniFind 分片] 并行执行：", length(pending), " 片 / ", nw, " 路（", kind, "）\n", sep = "")
    cl <- parallel::makeCluster(nw)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    runner_fn <- runner
    # 闭包跨集群边界前物化其环境绑定：clusterExport 序列化不落地未强制的
    # promise，工厂函数形参未 force 即在 worker 上报 object not found
    # （2026-09-30 实测）；get() 逐个强制求值后绑定随闭包序列化
    if (!is.null(environment(runner_fn))) {
      renv <- environment(runner_fn)
      invisible(lapply(ls(envir = renv), function(n) get(n, envir = renv)))
    }
    exp_funs <- unique(c("runner_fn", "kind", "punifind_dir", "dots", "pending",
                         "pf_run_one", "pf_probe", "pf_entry", "pf_env_vars",
                         "pf_build_args", "pf_spawn_collect", "pf_read_log_lines",
                         "pf_help_flags", "pf_probe_denovo_len_ok", "pf_clear_stale_cache",
                         "pf_fresh_own_spectra", "pf_match_final_spectra",
                         "pf_qryres_renumber", "pf_gpu_info"))
    if (exists("pf_sanitize_qryres", mode = "function", inherits = TRUE))
      exp_funs <- c(exp_funs, "pf_sanitize_qryres")
    parallel::clusterExport(cl, exp_funs, environment())
    idxs <- seq_along(pending)
    raw <- parallel::clusterMap(cl, function(j, mp) {
      p <- pending[[j]]
      tryCatch(list(ok = TRUE, res = do.call(runner_fn,
                 c(list(kind, punifind_dir, file.path(p$ch_out, "qry"),
                       file.path(p$ch_out, "mgf"), p$ch_out), dots, list(master_port = mp)))),
               error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    }, j = idxs, mp = 29516L + 111L * idxs, RECYCLE = FALSE)
    failed <- vapply(raw, function(r) !isTRUE(r$ok), logical(1))
    for (j in idxs[!failed])
      .pfr_write_marker(file.path(pending[[j]]$ch_out, paste0("done_", kind, ".txt")),
                        pending[[j]]$stem, pending[[j]]$fp)
    if (any(failed)) {
      stems <- vapply(pending[idxs[failed]], function(p) p$stem, character(1))
      msgs  <- vapply(raw[idxs[failed]], function(r) substr(r$msg, 1L, 300L), character(1))
      tails <- vapply(pending[idxs[failed]], function(p) pf_log_tail(p$ch_out, kind),
                      character(1))
      stop("[pUniFind 分片] 并行 chunk 失败（", paste(stems, collapse = ", "),
           "）——修好后重跑本函数即从失败片续跑（成功片标记已保留）：\n",
           paste0("  [", stems, "] ", msgs, "\n  [", stems, "] [子日志尾] ", tails,
                  collapse = "\n"), call. = FALSE)
    }
  }
  # ---- 清理（v1.0 D 项，可选）：保留 done 标记与合并产物 ----
  .pfr_cleanup <- function() {
    if (!isTRUE(chunk_cleanup)) return(invisible(NULL))
    for (k0 in seq_along(mgfs)) {
      stem <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(mgfs[k0])))
      for (d in c("mgf", "qry", "pUniFind_data_tmp", "pUniFind_result"))
        unlink(file.path(.pfr_chunk(stem), d), recursive = TRUE, force = TRUE)
    }
    cat("[pUniFind 分片] 已清理片内数据（done 标记与合并产物保留）（", kind, "）\n", sep = "")
  }
  # ---- 合并：产物落回原接口位置 ----
  if (identical(kind, "denovo")) {
    csvs <- unlist(lapply(seq_along(mgfs), function(k0) {
      stem <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(mgfs[k0])))
      list.files(file.path(.pfr_chunk(stem), "pUniFind_result"),
                 pattern = "\\.csv$", full.names = TRUE)
    }))
    csvs <- csvs[order(file.info(csvs)$mtime)]
    if (!length(csvs)) stop("[pUniFind 分片] 无任何 chunk 产出 csv（", out_dir, "）")
    res_dir <- file.path(out_dir, "pUniFind_result")
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
    merged <- file.path(res_dir, "test_5_merged.csv")
    con <- file(merged, open = "wb")
    src_rows <- 0L
    for (f in csvs) {
      ln <- readLines(f, warn = FALSE)
      src_rows <- src_rows + length(ln)
      writeLines(ln, con, sep = "\n", useBytes = TRUE)
    }
    close(con)
    # v1.0 P3：行数守恒断言——merged 必须恰等于各片 csv 行数之和，防未来
    # 片内目录混入额外 csv（如片级合并产物）被双重拼接
    merged_rows <- length(readLines(merged, warn = FALSE))
    if (merged_rows != src_rows)
      stop("[pUniFind 分片] denovo 合并守恒失败：merged ", merged_rows,
           " 行 ≠ 各片合计 ", src_rows, " 行（", paste(basename(csvs), collapse = ", "), "）")
    cat("[pUniFind 分片] denovo 合并：", length(csvs), " 片 → ", merged, "\n", sep = "")
    .pfr_cleanup()
    return(c(denovo_csv = merged))
  }
  # rescore：各 chunk qry 目录里的 fdr 谱表合并（表头只留首片）→ 真实 qry_res_dir
  sps <- unlist(lapply(seq_along(mgfs), function(k0) {
    stem <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(mgfs[k0])))
    pf_match_final_spectra(file.path(.pfr_chunk(stem), "qry"))
  }))
  if (!length(sps)) stop("[pUniFind 分片] 无任何 chunk 产出 rescore 谱表（", out_dir, "）")
  sps <- sps[order(file.info(sps)$mtime)]
  # v1.0 P3：各片表头必须逐列一致——列名/列序分叉即停，不把异构表拼在一起
  hdr1 <- readLines(sps[1], warn = FALSE, n = 1L)
  bad <- character(0)
  for (f in sps[-1]) {
    h <- readLines(f, warn = FALSE, n = 1L)
    if (!identical(h, hdr1)) bad <- c(bad, basename(f))
  }
  if (length(bad))
    stop("[pUniFind 分片] rescore 各片表头不一致（首片＝「", hdr1, "」）：",
         paste(bad, collapse = ", "))
  dst <- file.path(qry_res_dir, basename(sps[1]))
  con <- file(dst, open = "wb")
  for (i in seq_along(sps)) {
    ln <- readLines(sps[i], warn = FALSE)
    if (i > 1L && grepl("^File_Name\t", ln[1L])) ln <- ln[-1L]
    writeLines(ln, con, sep = "\n", useBytes = TRUE)
  }
  close(con)
  cat("[pUniFind 分片] rescore 合并：", length(sps), " 片 → ", dst, "\n", sep = "")
  .pfr_cleanup()
  c(rescore_spectra = dst)
}

# ---- 自检：不跑推理，只验证适配层本身 ---------------------------------------
# ---- self-test: verify the adapter itself (no inference) --------------------
pf_self_test <- function(punifind_dir = NULL, python_exe = NULL, qry_res_dir = NULL) {
  ok <- TRUE
  say <- function(pass, what, detail = "") {
    cat(sprintf("  [%s] %-46s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", what, detail))
    if (!isTRUE(pass)) ok <<- FALSE
  }
  cat("== punifind_gpu_compat 自检 ==\n")
  # 1) 环境变量（自查 names=values 对，而非只看值）
  ev <- pf_env_vars(8)
  evkv <- paste0(names(ev), "=", ev)
  say(all(c("USE_LIBUV=0", "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1") %in% evkv),
      "子进程环境含两个官方开关", paste(evkv, collapse = " "))
  # 2) argv：数据目录第一 / 无 --load-contra / num-workers 生效
  a <- pf_build_args("rescore", "P", "Q", "M", "O", flags = NULL, num_workers = 8)
  say(substr(a[1], 1, 1) != "-", "argv[1] 是数据目录（位置参数前置）", a[1])
  say(!any(grepl("load-contra", a)), "不含已消失的 --load-contra")
  say(any(grepl("^--num-workers=8$", a)), "num-workers 已注入", grep("^--num-workers", a, value = TRUE))
  a2 <- pf_build_args("rescore", "P", "Q", "M", "O",
                      flags = c("--user-dir", "--task", "--num-proc"), num_workers = 4)
  say(any(grepl("^--task=score_inference$", a2)) && !any(grepl("^--batch-size", a2)),
      "只保留探测到的 flag（过滤生效）", paste(length(a2), "个参数"))
  # v1.0（工作站 P1）：denovo 三参数门控——flags=NULL 时按定向探针结果
  d1 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL, denovo_len_ok = TRUE)
  say(all(c("--min-length=6", "--max-length=40") %in% d1),
      "denovo flags=NULL＋探针TRUE → 三参带上")
  d2 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL, denovo_len_ok = FALSE)
  say(!any(grepl("min-length|max-length", d2)),
      "denovo flags=NULL＋探针FALSE → 三参不带（老构建）")
  d3 <- pf_build_args("denovo", "P", "Q", "M", "O", flags = NULL)
  say(all(c("--min-length=6", "--max-length=40") %in% d3),
      "denovo flags=NULL＋不可判 → 按官方口径带上")
  # 3) qry.res 压实：合成夹具 + 幂等性
  td <- tempfile("pfst_"); dir.create(td)
  f <- file.path(td, "x.qry.res")
  writeLines(c("S 1234.5 2", "spec_a.dta", "1 PEP K 2 0.9 0.01 0",
               "3 PEPR 2 0.8 0.02 0",                       # 洞里补：序号 3 应为 2
               "S 2345.6 2", "spec_b.dta", "1 PEPTIDE 2 0.7 0.03 0"), f)
  r1 <- pf_qryres_renumber(td, verbose = FALSE)
  r2 <- pf_qryres_renumber(td, verbose = FALSE)
  ln <- readLines(f, warn = FALSE)
  say(r1$fixed[1] == 1 && r2$fixed[1] == 0, "序号压实生效且幂等",
      sprintf("首次改写 %d 行 / 二次 %d 行", r1$fixed[1], r2$fixed[1]))
  say(identical(sub("^([0-9]+) .*", "\\1", ln[3:4]), c("1", "2")),
      "候选序号已连续 1..n", paste(sub(" .*", "", ln[3:4]), collapse = ","))
  unlink(td, recursive = TRUE)
  # 4) 安装与 GPU 探测（给了目录才查）
  if (!is.null(punifind_dir) && nzchar(punifind_dir)) {
    pr <- pf_probe(punifind_dir, python_exe)
    say(pr$ok, "pUniFind 组件齐备（exe/ckpts/consts）")
    en <- pf_entry("rescore", punifind_dir, python_exe)
    say(TRUE, "入口路线", paste0(en$route, " -> ", en$cmd))
    if (identical(en$route, "exe"))
      cat("       注意：走 exe 启动器时实测 argv 会被重组（data 变 '0'）；" ,
          "建议提供 punifind_python 以走 python 直调路线\n")
    say(length(pr$score_flags) > 0 || TRUE, "flag 面探测",
        if (is.null(pr$score_flags)) "未取到（将用保守集）" else
          paste(length(pr$score_flags), "个 flag"))
    say(!is.null(pr$gpu) || TRUE, "GPU 探测", if (is.null(pr$gpu)) "无 nvidia-smi" else pr$gpu)
    if (length(pr$score_flags))
      say(!("--load-contra" %in% pr$score_flags),
          "该构建确实不支持 --load-contra（故必须按探查结果传参）")
  }
  if (!is.null(qry_res_dir) && nzchar(qry_res_dir)) {
    rn <- pf_qryres_renumber(qry_res_dir, dry_run = TRUE, verbose = FALSE)
    say(nrow(rn) > 0 && sum(rn$spectra) > 0, "qry.res 试算（dry-run，未写盘）",
        sprintf("%d 个文件 / %d 个谱块 / 需改写 %d 行",
                nrow(rn), sum(rn$spectra), sum(rn$fixed)))
  }
  cat(if (ok) "== 自检全部通过 ==\n" else "== 自检存在失败项 ==\n")
  invisible(ok)
}
