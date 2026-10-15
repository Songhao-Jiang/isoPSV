# =============================================================================
# two_queue_steward.R — v1.0 GPU 通道管家（两队列调度器的 GPU 道实现）
# two_queue_steward.R — v1.0 GPU-lane steward (GPU lane of the two-queue scheduler)
#
# 设计依据：WorkBuddy《isoPSV_DAG双通道两队列调度器提案_工作站设计稿_20261001》
# ＋开发侧评审四修复（亲和补位/跨队列唤醒由"两条物理串行道"结构性消解；
# rescore 前置含 open 后处理；引擎槽无假依赖——由 master 原生顺序保证）。
#
# 形态：master（CPU 道，前台串行：open→open-pp→restrict→restrict-pp→引擎→汇合）
#       ＋ steward（GPU 道，本模块，后台串行：denovo→[等 open-pp 标记]→rescore）。
# 两道各自内部串行、道间并行——等价于提案的两队列调度在该任务集上的执行轨迹：
#   T=0   master 起 open（前台）∥ steward 起 denovo（后台）
#   denovo 完 ∧ open-pp 完 → steward 立即转 rescore（GPU 无空窗）
#   rescore ∥ master 的 restrict/引擎段
#   master 到汇合门等 steward 状态 done 后走原有发现/合并逻辑（零合并口径变化）
#
# 失败语义（与 v1.0/76 现行一致）：denovo 失败→记录并继续（合并缺失策略 5，
# 下次运行重试）；rescore 失败→自动跳过（3 策略合并）；两者都不阻断 master。
#
# 产物（均在 config_dir）：
#   _tq_open_pp_done.txt     master 在 open 后处理完成后写入（管家等待的门牌）
#   _gpu_steward.pid         管家 PID（master 复跑时判活，防双管家）
#   _gpu_steward_status.txt  终态（state=done; denovo=…; rescore=…）
#   _gpu_steward.out/err.log 管家 stdout/stderr
# =============================================================================
options(stringsAsFactors = FALSE)

# Windows 下 source(file, encoding="UTF-8") 的
# connection 层偶发 "invalid input found on input connection"（文件本体合法
# UTF-8、iconv 逐行无坏字符仍复现），致 punifind_gpu_compat.R 的 pf_run 未定
# 义、自动 rescore 假失败。改为 readLines+parse+eval，语义等同 source。
#  source(file, encoding="UTF-8") sporadically fails at
# the connection layer on Windows ("invalid input found on input connection"
# even for valid UTF-8 with clean per-line iconv), leaving pf_run undefined and
# faking an auto-rescore failure. Replaced with readLines+parse+eval.
source_utf8 <- function(path) {
  l <- readLines(path, warn = FALSE, encoding = "UTF-8")
  eval(parse(text = paste(l, collapse = "
")), envir = globalenv())
}

# %||% 兜底：R 4.4+ 原生；4.1–4.3 目标环境下给同语义定义
# %||% fallback: native since R 4.4; provide identical semantics for 4.1-4.3
`%||%` <- function(a, b) if (is.null(a)) b else a

# PID 判活（跨平台）。v1.0 join 的 system2("tasklist", c("/FI","PID eq N")) 有
# 实测缺陷：system2 把 "PID eq N" 按空格拆成多参 → tasklist 报 Invalid argument
# 'eq'、恒判死（此前被"denovo 早死"掩盖）。改为整条命令经 shell/ps 执行。
# PID liveness probe (cross-platform). The v1.0 join pattern
# system2("tasklist", c("/FI", "PID eq N")) splits "PID eq N" at spaces into
# separate arguments, tasklist errors with Invalid argument 'eq' and the probe
# always reports dead (previously masked because denovo had died early).
# Pass the whole command through shell/ps instead.
tq_pid_alive <- function(pid) {
  if (length(pid) == 0 || is.na(pid) || pid <= 0) return(FALSE)
  out <- if (.Platform$OS.type == "windows") {
    suppressWarnings(shell(sprintf('tasklist /FI "PID eq %d" /NH', as.integer(pid)),
                           intern = TRUE))
  } else {
    suppressWarnings(system2("ps", c("-p", as.character(as.integer(pid))),
                             stdout = TRUE, stderr = FALSE))
  }
  any(grepl(as.character(as.integer(pid)), out, fixed = TRUE))
}

# ---- 自含配置解析（复刻 user_step2 口径；不 source 驱动脚本，避免副作用）----
# ---- self-contained config resolution mirroring user_step2 semantics ----
steward_read_cfg <- function(cfg2_path) {
  cfg2 <- read_param_file(cfg2_path)
  base_path <- resolve_cfg_paths(cfg2, "base_config")[1]
  base_cfg <- if (!is.na(base_path) && nzchar(base_path) && file.exists(base_path))
    read_param_file(base_path) else list()
  list(cfg2 = cfg2, base_cfg = base_cfg)
}

steward_kv <- function(src, key) {
  v <- src[[key]]
  if (is.null(v) || length(v) == 0) "" else trimws(as.character(tail(v, 1)))
}

steward_main <- function(cfg2_path, code_dir) {
  cfgs <- steward_read_cfg(cfg2_path)
  cfg2 <- cfgs$cfg2; base_cfg <- cfgs$base_cfg
  cdir <- cfg2$config_dir
  status_f <- file.path(cdir, "_gpu_steward_status.txt")
  writeLines(c(paste0("state=running pid=", Sys.getpid()), paste0("at=", format(Sys.time()))),
             status_f, useBytes = TRUE)
  on.exit({
    if (!grepl("^state=done", readLines(status_f, warn = FALSE)[1] %||% "")) {
      writeLines(c("state=aborted", paste0("at=", format(Sys.time()))), status_f, useBytes = TRUE)
    }
  }, add = TRUE)

  # ---- pUniFind 开关与旋钮（口径同 master）----
  use_puni <- any(vapply(list(base_cfg, cfg2), function(s)
    steward_kv(s, "use_punifind") == "1", logical(1)))
  puni_exe <- ""; puni_py <- ""
  for (s in list(base_cfg, cfg2)) {
    v <- steward_kv(s, "punifind_exe");  if (nzchar(v)) puni_exe <- v
    v <- steward_kv(s, "punifind_python"); if (nzchar(v)) puni_py <- v
  }
  puni_exe <- normalizePath(puni_exe, winslash = "/", mustWork = FALSE)
  if (nzchar(puni_py)) puni_py <- normalizePath(puni_py, winslash = "/", mustWork = FALSE)
  auto_rescore <- any(vapply(list(base_cfg, cfg2), function(s)
    steward_kv(s, "punifind_rescore") == "1", logical(1)))
  auto_denovo <- any(vapply(list(base_cfg, cfg2), function(s)
    steward_kv(s, "punifind_denovo") == "1", logical(1)))
  knob <- function(key, default) {
    for (s in list(base_cfg, cfg2)) {
      v <- suppressWarnings(as.integer(steward_kv(s, key)))
      if (!is.na(v)) return(v)
    }
    default
  }
  pf_nw <- if (.Platform$OS.type == "windows") knob("punifind_num_workers", 0L)
           else knob("punifind_num_workers", 8L)
  pf_np <- knob("punifind_num_proc", 8L)
  pf_bs <- knob("punifind_batch_size", 64L)
  pf_cp <- knob("punifind_chunk_parallel", 1L)
  pf_cu <- identical(steward_kv(cfg2, "punifind_chunk_cleanup"), "1") ||
           identical(steward_kv(base_cfg, "punifind_chunk_cleanup"), "1")

  source_utf8(file.path(code_dir, "punifind_gpu_compat.R"))
  unlab_dir <- file.path(cdir, "unlabeled")
  fasta_path <- resolve_cfg_paths(base_cfg, "fasta")[1]

  denovo_state <- "disabled"; rescore_state <- "disabled"
  if (!use_puni || !nzchar(puni_exe)) {
    cat("[GPU 管家] pUniFind 未启用或无安装路径，GPU 道空跑\n")
    cat("[GPU steward] pUniFind disabled or no install path; GPU lane no-op\n")
    writeLines(c("state=done", "denovo=disabled", "rescore=disabled"),
               status_f, useBytes = TRUE)
    return(invisible(TRUE))
  }

  # ---- 阶段 1：denovo（与 master 前台/dag 版完全同参；缓存命中即跳过）----
  dproj <- file.path(cdir, "pUniFind_denovo")
  if (auto_denovo) {
    stray <- list.files(unlab_dir, pattern = "_check\\.mgf$", full.names = TRUE)
    if (length(stray) > 0) {
      stray_dir <- file.path(cdir, "_stray_mgf")
      dir.create(stray_dir, recursive = TRUE, showWarnings = FALSE)
      file.rename(stray, file.path(stray_dir, basename(stray)))
    }
    dir.create(file.path(dproj, "pUniFind_result"), recursive = TRUE, showWarnings = FALSE)
    csvs <- list.files(file.path(dproj, "pUniFind_result"), pattern = "\\.csv$",
                       full.names = TRUE)
    if (length(csvs) > 0) {
      denovo_state <- "skip(cache)"
      cat("[GPU 管家] denovo 缓存命中（", length(csvs), " csv），跳过\n", sep = "")
      cat("[GPU steward] denovo cache hit (", length(csvs), " csv), skipped\n", sep = "")
    } else {
      cat("[GPU 管家] denovo 开始 @ ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
      cat("[GPU steward] denovo start @ ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
      ok <- tryCatch({
        pf_run("denovo", puni_exe, file.path(dproj, "result"), unlab_dir, dproj,
               python_exe = puni_py, num_workers = pf_nw, num_proc = pf_np,
               batch_size = pf_bs, chunk_parallel = pf_cp, chunk_cleanup = pf_cu)
        TRUE
      }, error = function(e) {
        cat("[GPU 管家] denovo 失败（记录后继续；下次运行重试）：", conditionMessage(e), "\n")
        cat("[GPU steward] denovo failed (recorded; retried next run): ",
            conditionMessage(e), "\n")
        FALSE
      })
      denovo_state <- if (ok) "ok" else "fail"
    }
  }

  # ---- 阶段门：等 master 写 open-pp 完成标记（评审修复③：rescore 前置含后处理）----
  marker <- file.path(cdir, "_tq_open_pp_done.txt")
  if (auto_rescore) {
    waited <- 0L
    while (!file.exists(marker) && waited < 12L * 3600L) {
      Sys.sleep(30L); waited <- waited + 30L
      if (waited %% 300L == 0L)
        cat("[GPU 管家] 等 open 后处理标记... 已等 ", waited / 60, " min\n", sep = "")
    }
    if (!file.exists(marker)) {
      cat("[GPU 管家] 12h 未见 open-pp 标记（master 异常？），rescore 跳过\n")
      cat("[GPU steward] open-pp marker absent after 12h; rescore skipped\n")
      writeLines(c("state=done", paste0("denovo=", denovo_state), "rescore=fail(no-marker)"),
                 status_f, useBytes = TRUE)
      return(invisible(TRUE))
    }
  }

  # ---- 阶段 2：rescore（项目解析与调用逐参复刻 master；缓存命中即跳过）----
  if (auto_rescore) {
    pfo_tpl <- cfg2$pfind_cfg_open
    pfo_tpl <- if (!is.null(pfo_tpl) && length(pfo_tpl) > 0 && pfo_tpl[1] != "") {
      tail(resolve_cfg_paths(cfg2, "pfind_cfg_open"), 1)
    } else {
      ensure_pfind_template("open", fasta_path, cdir)
    }
    proj <- if (nzchar(pfo_tpl)) dirname(dirname(pfo_tpl)) else ""
    if (nzchar(proj)) {
      pac_src <- paste0(fasta_path, "_td.pac")
      if (file.exists(pac_src) && !file.exists(file.path(proj, basename(pac_src))))
        file.copy(pac_src, file.path(proj, basename(pac_src)))
      rr <- list.files(file.path(proj, "result"), pattern = "pUniFind\\.spectra$",
                       full.names = TRUE)
      if (length(rr) > 0) {
        rescore_state <- "skip(cache)"
        cat("[GPU 管家] rescore 缓存命中，跳过：", tail(rr, 1), "\n")
        cat("[GPU steward] rescore cache hit, skipped: ", tail(rr, 1), "\n")
      } else {
        cat("[GPU 管家] rescore 开始 @ ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
        cat("[GPU steward] rescore start @ ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
        ok <- tryCatch({
          pf_run("rescore", puni_exe, file.path(proj, "result"), unlab_dir, proj,
                 file.path(code_dir, "punifind_rescore_convert.py"), puni_py,
                 num_workers = pf_nw, num_proc = pf_np,
                 batch_size = pf_bs, chunk_parallel = pf_cp, chunk_cleanup = pf_cu)
          TRUE
        }, error = function(e) {
          cat("[GPU 管家] rescore 失败，自动跳过（不影响合并）：", conditionMessage(e), "\n")
          cat("[GPU steward] rescore failed, auto-skipped (merge unaffected): ",
              conditionMessage(e), "\n")
          FALSE
        })
        rescore_state <- if (ok) "ok" else "fail"
      }
    } else rescore_state <- "fail(no-project)"
  }

  writeLines(c("state=done", paste0("denovo=", denovo_state),
               paste0("rescore=", rescore_state)), status_f, useBytes = TRUE)
  cat("[GPU 管家] GPU 道完成：denovo=", denovo_state, " rescore=", rescore_state,
      " @ ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
  cat("[GPU steward] GPU lane done: denovo=", denovo_state, " rescore=", rescore_state, "\n")
  invisible(TRUE)
}

# master 侧：管家启动（生成 ASCII 后台脚本＋PowerShell 启动器，全部沿用
# v1.0 修订-6 的 Windows 实测口径：真 PID/无 BOM/\\u 转义/locale 清理/
# out-err 分离/目录预建/WorkingDirectory 定 code 目录）
# master side: launch the steward (ASCII bg script + PowerShell launcher reusing
# all six v1.0 Windows launch fixes verbatim)
tq_launch_steward <- function(cfg2, cfg2_path, code_dir) {
  cdir <- cfg2$config_dir
  status_f <- file.path(cdir, "_gpu_steward_status.txt")
  pid_f <- file.path(cdir, "_gpu_steward.pid")
  # 已完成或仍在跑：不重复起
  if (file.exists(status_f) &&
      grepl("^state=done", readLines(status_f, warn = FALSE)[1] %||% "")) {
    cat("[两队列] GPU 管家已完成（状态文件在），不重启\n")
    cat("[two-queue] GPU steward already done (status file present), not relaunched\n")
    return(invisible(TRUE))
  }
  if (file.exists(pid_f)) {
    pid_old <- suppressWarnings(as.integer(readLines(pid_f, warn = FALSE)[1]))
    if (!is.na(pid_old)) {
      if (tq_pid_alive(pid_old)) {
        cat("[两队列] GPU 管家仍在运行（PID ", pid_old, "），不重启，汇合门将等待\n", sep = "")
        cat("[two-queue] GPU steward still running (PID ", pid_old,
            "), not relaunched; join gate will wait\n", sep = "")
        return(invisible(TRUE))
      }
    }
  }
  esc_nonascii <- function(s) {
    s <- enc2utf8(as.character(s))
    paste(vapply(utf8ToInt(s), function(cp)
      if (cp < 128L) intToUtf8(cp)
      else if (cp <= 0xFFFFL) sprintf("\\u%04X", cp)
      else sprintf("\\U%08X", cp), character(1)), collapse = "")
  }
  script <- file.path(cdir, "_gpu_steward_bg.R")
  # 不能 source boot.R：boot 按 --file= 定位代码目录，而本入口脚本位于项目目
  # 录（非 code/），boot 会错找 <项目>/../code（冒烟实测）。改为按 code_dir
  # 绝对路径镜像 boot 的五个加载（清单与 boot.R 保持一致）。
  # boot.R must NOT be sourced here: it locates the code dir from the --file=
  # path, and this entry script lives in the PROJECT dir (not code/), so boot
  # probes <project>/../code and fails (live in the smoke run). Mirror boot's
  # five loads by absolute code_dir instead (keep in sync with boot.R).
  writeLines(c(
    "code_dir <- getwd()",
    "source(file.path(code_dir, 'compat.R'))",
    "suppressMessages(source(file.path(code_dir, 'acs_core.R')))",
    "suppressMessages(source(file.path(code_dir, 'paired_params.R')))",
    "suppressMessages(source(file.path(code_dir, 'ReadMGF_script.R')))",
    "suppressMessages(source(file.path(code_dir, 'chimeric_check.R')))",
    "source(file.path(code_dir, 'two_queue_steward.R'), encoding='UTF-8')",
    sprintf("writeLines(as.character(Sys.getpid()), '%s', useBytes=TRUE)",
            esc_nonascii(pid_f)),
    sprintf("steward_main(cfg2_path='%s', code_dir=getwd())", esc_nonascii(cfg2_path))
  ), script, useBytes = TRUE)

  rscript_full <- file.path(R.home("bin"), "Rscript")
  out_log <- file.path(cdir, "_gpu_steward.out.log")
  err_log <- file.path(cdir, "_gpu_steward.err.log")
  ps_cmd <- sprintf(
    paste0('$env:LC_CTYPE=$null; $env:LC_ALL=$null; $env:LANG=$null; ',
           'Start-Process -FilePath "%s" -ArgumentList \'"%s"\' ',
           '-WorkingDirectory "%s" -WindowStyle Hidden ',
           '-RedirectStandardOutput "%s" -RedirectStandardError "%s" ',
           '-PassThru | Select-Object -ExpandProperty Id'),
    rscript_full, script, code_dir, out_log, err_log)
  pid <- suppressWarnings(as.integer(
    system2("powershell", c("-NoProfile", "-Command", shQuote(ps_cmd)),
            stdout = TRUE, stderr = FALSE, invisible = TRUE)))
  if (is.na(pid)) pid <- -1L
  if (pid > 0) writeLines(as.character(pid), pid_f, useBytes = TRUE)
  cat("[两队列] GPU 管家已启动 PID=", pid, "  日志=", out_log, "\n", sep = "")
  cat("[two-queue] GPU steward launched, PID=", pid, "\n", sep = "")
  invisible(pid > 0)
}

# master 侧：汇合门——等管家终态；管家进程死而无终态→告警放行（发现段自会
# 处理缺失产物，保持 v1.0 失败容忍语义）。管家活着则持续等待（无自动放行，
# 防止合并与 rescore 并发写盘竞态）。
# master side: join gate - wait for the steward terminal state; a dead steward
# without a terminal state warns and proceeds (the discovery sections tolerate
# missing outputs per v1.0 semantics). A live steward is waited on with no
# auto-bail (prevents merge-vs-rescore write races).
tq_join_gate <- function(cfg2, poll_s = 30L) {
  cdir <- cfg2$config_dir
  status_f <- file.path(cdir, "_gpu_steward_status.txt")
  pid_f <- file.path(cdir, "_gpu_steward.pid")
  waited <- 0L
  repeat {
    if (file.exists(status_f)) {
      ln <- readLines(status_f, warn = FALSE)
      if (length(ln) > 0 && grepl("^state=done", ln[1])) {
        cat("[两队列] 汇合：GPU 道完成（", paste(ln[-1], collapse = "; "), "）\n", sep = "")
        cat("[two-queue] join: GPU lane done (", paste(ln[-1], collapse = "; "), ")\n", sep = "")
        return(invisible(TRUE))
      }
      if (length(ln) > 0 && grepl("^state=aborted", ln[1])) {
        warning("GPU 管家中止（见 _gpu_steward.err.log）；按缓存缺失口径继续")
        return(invisible(FALSE))
      }
    }
    alive <- FALSE
    if (file.exists(pid_f)) {
      pid <- suppressWarnings(as.integer(readLines(pid_f, warn = FALSE)[1]))
      if (!is.na(pid)) alive <- tq_pid_alive(pid)
    }
    if (!alive) {
      if (!file.exists(status_f)) {
        cat("[两队列] GPU 管家未留状态即退出（异常）；按缓存缺失口径继续\n")
        cat("[two-queue] GPU steward exited without a status (abnormal); continue on-cache\n")
      } else {
        cat("[两队列] GPU 管家进程已消失而状态未 done；按缓存缺失口径继续\n")
        cat("[two-queue] steward gone before done; continue on-cache\n")
      }
      return(invisible(FALSE))
    }
    Sys.sleep(poll_s); waited <- waited + poll_s
    if (waited %% 300L == 0L)
      cat("[两队列] 等 GPU 道... 已等 ", waited / 60, " min\n", sep = "")
  }
}
