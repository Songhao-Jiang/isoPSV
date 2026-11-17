# =============================================================================
# pred_alphapeptdeep.R —— AlphaPeptDeep(peptdeep) 理论谱图预测器 L1 适配器
# （工作站自写 2026-10-05；结构同 pred_pdeep3.R，契约同 pred_ms2pip.R）
#
# 把 isoPSV 的 merge/peptide.txt（peptide / modification / charge 三列）转成
# driver 输入，经 alphapeptdeep_driver.py 调用 peptdeep（MannLabs，
# ModelManager.predict_ms2），再把输出转回与 pDeep2 predict.txt 同构的
# isoPSV 预测表（BEGIN IONS / TITLE=<seq>|<mod>|<charge> / b+1 / b+2 /
# y+1 / y+2——b 正序、y 位点序＝rev(y 正序)，与 read_predict_txt 对齐）。
#
# 用法：
#   Rscript pred_alphapeptdeep.R --in merge/peptide.txt \
#          --out merge/predict_alphapeptdeep.txt \
#          --python <python-with-peptdeep> --params <instrument,ce，如 Lumos,30>
#   Rscript pred_alphapeptdeep.R --self-test   # 无需 peptdeep，只测适配层
# 退出码：0 成功；非 0 失败（user_run 据此把该预测列标缺失并继续）。
# 依赖：仅 base R；peptdeep 运行于用户提供的 python 环境（本机
# C:\software\miniconda3\envs\alphapeptdeep）。
# 行序口径（2026-10-05 实测定标）：peptdeep 碎片表 b 行＝正序 b1..b_{n-1}，
# y 行＝位点序 y_{n-1}..y1——driver 已翻转为 y 正序，本适配器组块时再翻回。
# =============================================================================

# v2.0.179：parse_cli／指纹侧车三件套（input_fp/sidecar_matches/
# write_run_params）已上收到 predictors/_common.R（三适配器原为逐字三胞胎）；
# --file= 同目录定位，任何 cwd 均可。
# v2.0.179: parse_cli and the fingerprint-sidecar trio moved up to
# predictors/_common.R (they were verbatim triplets); located via the
# --file= sibling dir, any cwd works.
.sc_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", tail(f, 1))) else "predictors"
})
source(file.path(.sc_dir, "_common.R"))
# （parse_cli 自 v2.0.179 起来自 _common.R——上方 source 引导）

pfind_mods_strip <- function(mod) {
  if (is.null(mod) || is.na(mod) || !nzchar(mod)) return("")
  sub(";+$", "", trimws(as.character(mod)))
}

build_table <- function(pep_df) {
  n <- nrow(pep_df)
  data.frame(
    spec_id = sprintf("p%06d", seq_len(n)),
    peptide = as.character(pep_df$peptide),
    modinfo = vapply(seq_len(n), function(i)
      pfind_mods_strip(pep_df$modification[i]), character(1)),
    charge = as.character(pep_df$charge),
    orig_mod = sub(";+$", "", as.character(pep_df$modification)),
    stringsAsFactors = FALSE)
}

block_lines <- function(title, b1, b2, y1, y2) {
  # 块格式参考实现（主/兜底路径均为向量化直写；此函数保留作文档与单测锚）
  # reference block layout (both live paths write vectorized; kept for docs)
  c("BEGIN IONS",
    paste0("TITLE=", title),
    paste0("b+1=", paste(b1, collapse = ",")),
    paste0("b+2=", paste(b2, collapse = ",")),
    paste0("y+1=", paste(rev(y1), collapse = ",")),
    paste0("y+2=", paste(rev(y2), collapse = ",")),
    "END IONS")
}

# v2.0.166（工作站 11-05 T2 计时反例）：driver 产物 → predict.txt 的向量化
# 重建。165 线性化只修掉了 c() 累积（输入装填 ≤1 分钟已达标），但逐行
# strsplit×4＋as.numeric×4＋format×4 仍在——66.8 万行重入实测 79 分钟单核。
# 本版两处换掉：①主路径＝driver 直写 predict.txt（--predict-out，python 侧
# 写字符串快一个量级，适配器退化为透传）；②fallback（driver 产物在、
# predict 缺失时的重建）用本函数——数值串原样透传（driver %.5f 与旧
# format(as.numeric) 回读后同一双精度，仅尾零字节差），无逐行数值解析。
# 心跳打到 stdout（165 打 stderr 在 master 捕获链下不可见——工作站实测
# 全程零心跳行）。
# v2.0.166 (workstation 11-05 T2 counterexample): vectorized rebuild from the
# driver output. The 165 linearization only fixed the c() accumulation (input
# loading <=1 min, on target), but the per-row strsplit/as.numeric/format
# remained - 79 single-core minutes on a 668k-row re-entry. This version
# replaces both: (1) main path = the driver writes predict.txt directly
# (--predict-out, python string writing is an order of magnitude faster, the
# adapter becomes a passthrough); (2) the fallback rebuild (driver output
# present, predict missing) uses this function - value strings pass through
# verbatim (driver %.5f reparses to the same doubles as the old
# format(as.numeric) output, trailing-zero byte differences only). The
# heartbeat goes to stdout (165's stderr was invisible under the master's
# capture chain - zero heartbeat lines in the workstation run).
build_predict_txt_fast <- function(input_table, drv_df) {
  n <- nrow(input_table)
  m <- match(input_table$spec_id, drv_df$spec_id)
  have <- !is.na(m)
  skip_no_rows <- sum(!have)
  residue_count <- nchar(input_table$peptide)
  comma_n <- function(x) nchar(x) - nchar(gsub(",", "", x, fixed = TRUE)) + 1L
  b1s <- rep("", n); b2s <- rep("", n); y1s <- rep("", n); y2s <- rep("", n)
  hh <- which(have)
  if (length(hh)) {
    b1s[hh] <- drv_df$b1[m[hh]]; b2s[hh] <- drv_df$b2[m[hh]]
    y1s[hh] <- drv_df$y1[m[hh]]; y2s[hh] <- drv_df$y2[m[hh]]
  }
  len_ok <- rep(FALSE, n)
  if (length(hh)) {
    len_ok[hh] <- nzchar(b1s[hh]) & nzchar(y1s[hh]) &
      comma_n(b1s[hh]) == residue_count[hh] - 1L &
      comma_n(y1s[hh]) == residue_count[hh] - 1L
  }
  keep <- which(have & len_ok)
  skip_len <- sum(have & !len_ok)
  nk <- length(keep)
  # y 行反转（正序 → 位点序）：唯一残余的逐行循环，每 10 万行打心跳
  # y reversal (positive -> site order): the only per-row loop left,
  # heartbeat every 100k rows
  y1r <- character(nk); y2r <- character(nk)
  if (nk) {
    sp1 <- strsplit(y1s[keep], ",", fixed = TRUE)
    sp2 <- strsplit(y2s[keep], ",", fixed = TRUE)
    t0 <- proc.time()[3]
    for (k in seq_len(nk)) {
      y1r[k] <- paste(rev(sp1[[k]]), collapse = ",")
      y2r[k] <- paste(rev(sp2[[k]]), collapse = ",")
      if (k %% 100000L == 0L) {
        cat(sprintf("[pred_alphapeptdeep] 重建进度 %d/%d（已用 %.0fs）\n",
                    k, nk, proc.time()[3] - t0))
      }
    }
  }
  out <- character(7L * nk)
  if (nk) {
    i1 <- seq.int(1L, by = 7L, length.out = nk)
    out[i1] <- "BEGIN IONS"
    out[i1 + 1L] <- paste0("TITLE=", input_table$peptide[keep], "|",
                           input_table$orig_mod[keep], "|",
                           input_table$charge[keep])
    out[i1 + 2L] <- paste0("b+1=", b1s[keep])
    out[i1 + 3L] <- paste0("b+2=", b2s[keep])
    out[i1 + 4L] <- paste0("y+1=", y1r)
    out[i1 + 5L] <- paste0("y+2=", y2r)
    out[i1 + 6L] <- "END IONS"
  }
  list(lines = out, skip_no_rows = skip_no_rows, skip_len = skip_len)
}

# ---- v2.0.168：产物指纹复用（升级 166 的两层存在性复用；工作站 run7 实锤：----
# ---- "predict＋driver 均在"复用了 11-04 23:59 的旧 predict 而 peptide.txt ----
# ---- 12:45 刚重建，陈旧预测流入 step3——策略增删/换库场景必错）。指纹＝     ----
# ---- 输入 md5＋size（master 的 file.copy 每轮刷新 peptides.txt mtime，      ----
# ---- mtime 序恒失效），侧车 <out>.run_params。旧版产物无侧车＝不可证新鲜→   ----
# ---- 重跑（一次性代价，升版后首轮全量）。
# ---- v2.0.168: fingerprint-based artifact reuse (upgrading the 166
# ---- existence-only two-tier check; workstation run7 proved it reused a
# ---- predict from 11-04 23:59 while peptide.txt had just been rebuilt at
# ---- 12:45, letting stale predictions flow into step3 - guaranteed wrong
# ---- under strategy/library changes). Fingerprint = input md5 + size in an
# ---- <out>.run_params sidecar (the master's file.copy refreshes peptides.txt
# ---- mtime every run, so mtime ordering never holds). Legacy artifacts
# ---- without a sidecar are unverifiable and rerun once.
# （指纹侧车三件套自 v2.0.179 起来自 _common.R）

self_test <- function() {
  cat("[self-test] pred_alphapeptdeep 适配层\n")
  pep_df <- data.frame(peptide = c("PEPTIDEK", "ALAMAR"),
                       modification = c("2,Oxidation[M];", ""),
                       charge = c(2L, 3L), stringsAsFactors = FALSE)
  tab <- build_table(pep_df)
  stopifnot(tab$modinfo[1] == "2,Oxidation[M]", tab$modinfo[2] == "")
  cat("  [ok] 修饰串转换（剥尾分号/空串）\n")
  # v2.0.166：driver 产物直连（spec_id 对齐）＋数值串原样透传（%.5f 不重排）
  # v2.0.166: direct driver-table join on spec_id + verbatim value strings
  # (%.5f passes through, not reformatted)
  drv <- data.frame(spec_id = tab$spec_id,
                    b1 = c("1.00000,2.00000,3.00000,4.00000,5.00000,6.00000,7.00000",
                           "1,2,3,4,5"),
                    b2 = c("0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000",
                           "0.1,0.2,0.3,0.4,0.5"),
                    y1 = c("7,6,5,4,3,2,1", "5,4,3,2,1"),
                    y2 = c("0.7,0.6,0.5,0.4,0.3,0.2,0.1", "0.5,0.4,0.3,0.2,0.1"),
                    stringsAsFactors = FALSE)
  built <- build_predict_txt_fast(tab, drv)
  lines <- built$lines
  stopifnot(built$skip_no_rows == 0L, built$skip_len == 0L,
            length(lines) == 2L * 7L,
            sum(lines == "BEGIN IONS") == 2L, sum(lines == "END IONS") == 2L,
            lines[2] == "TITLE=PEPTIDEK|2,Oxidation[M]|2",
            lines[3] == "b+1=1.00000,2.00000,3.00000,4.00000,5.00000,6.00000,7.00000",
            lines[4] == "b+2=0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000",
            lines[5] == "y+1=1,2,3,4,5,6,7")
  cat("  [ok] 块结构/TITLE 键/y 位点序反转/数值串透传\n")
  key2 <- paste(tab$peptide[2], sub(";+$", "", pep_df$modification[2]),
                pep_df$charge[2], sep = "|")
  stopifnot(any(lines == paste0("TITLE=", key2)))
  cat("  [ok] TITLE 键与消费端口径一致\n")
  # 跳过计数：无产物（spec_id 缺席）与长度不符各一例
  tab2 <- rbind(tab, data.frame(spec_id = "p000009", peptide = "PEPTIDEK",
                                modinfo = "", orig_mod = "", charge = "2",
                                stringsAsFactors = FALSE))
  drv2 <- rbind(drv, data.frame(spec_id = "p000009",
                                b1 = "1,2", b2 = "1,2", y1 = "1,2", y2 = "1,2",
                                stringsAsFactors = FALSE))
  b2t <- build_predict_txt_fast(tab2, drv2)
  stopifnot(b2t$skip_no_rows == 0L, b2t$skip_len == 1L,
            length(b2t$lines) == 2L * 7L)
  tab3 <- rbind(tab, data.frame(spec_id = "p000010", peptide = "ALAMAR",
                                modinfo = "", orig_mod = "", charge = "3",
                                stringsAsFactors = FALSE))
  b3t <- build_predict_txt_fast(tab3, drv)
  stopifnot(b3t$skip_no_rows == 1L, b3t$skip_len == 0L)
  cat("  [ok] 跳过计数（无产物/长度不符）\n")
  cat("[self-test] 全部通过\n")
  TRUE
}

main <- function() {
  cli_args <- parse_cli()
  if (isTRUE(cli_args$`self-test`)) return(if (self_test()) 0L else 1L)
  required_args <- c("in", "out", "python")
  missing_args <- setdiff(required_args, names(cli_args))
  if (length(missing_args)) stop("缺参数: ", paste(paste0("--", missing_args), collapse = " "))

  pep_df <- read.delim(cli_args$`in`, sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = character(0))
  pep_df$modification[is.na(pep_df$modification)] <- ""
  missing_cols <- setdiff(c("peptide", "modification", "charge"), colnames(pep_df))
  if (length(missing_cols)) stop("peptide.txt 缺列: ", paste(missing_cols, collapse = ","))
  # v2.0.170（工作站 169 报告建议，ms2pip 管道符防线推广）：pFind 二选一
  # 异构体记法（肽序列含 "|"，全表实测 918 条）对 peptdeep 是非法残基——
  # 此前靠 164 的二分重试逐批丢弃（代价高且无计数）；入侧统一跳过＋计数。
  # v2.0.170 (workstation 169-round request, generalizing the ms2pip pipe
  # guard): pFind's either-or isoform notation ("|" inside the peptide, 918
  # rows on the full table) is an illegal residue for peptdeep - previously
  # absorbed by the v164 bisection retry (expensive and uncounted); skip at
  # the input side with a count.
  pep_df <- guard_pipe_peptides(pep_df, "pred_alphapeptdeep")

  input_table <- build_table(pep_df)
  in_tsv <- paste0(cli_args$`out`, ".apd.tsv")
  # v2.0.166：orig_mod 列随输入 TSV 下发（driver 直写 TITLE 用，与消费端
  # read_predict_txt 的 sub(";+$","",...) 同口径；旧 driver 忽略该列）
  # v2.0.166: orig_mod rides along in the input TSV for the driver's direct
  # TITLE writing (same sub(";+$","",...) key as read_predict_txt; old
  # drivers ignore the extra column)
  write.table(input_table[, c("spec_id", "peptide", "modinfo", "charge", "orig_mod")],
              in_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  if (isTRUE(cli_args$`dry-run`)) {
    cat("dry-run: 输入 TSV 已写", in_tsv, "（", nrow(input_table), "行），未调用 peptdeep\n")
    return(0L)
  }

  params_raw <- if (is.null(cli_args$params)) "" else as.character(cli_args$params)
  params_parts <- strsplit(params_raw, ",", fixed = TRUE)[[1]]
  instrument <- if (length(params_parts) >= 1L && nzchar(trimws(params_parts[1])))
    trimws(params_parts[1]) else "Lumos"
  ce <- if (length(params_parts) >= 2L && nzchar(trimws(params_parts[2])))
    trimws(params_parts[2]) else "30"

  drv_out <- paste0(cli_args$`out`, ".driver.tsv")
  drv_log <- paste0(cli_args$`out`, ".driver.log")
  # v2.0.168：指纹制复用（升级 166 的存在性两层）。
  # ① predict 在且指纹匹配 → 复用（删除 driver.tsv 可强制重跑）；
  # ② predict 缺/删而 driver 产物在且指纹匹配 → R 侧向量化重建（分钟级）；
  # ③ 指纹失配（肽表变化）或无侧车（旧版产物）→ 全量重跑（旧 predict/driver
  #    删除防误读；一次性代价，升版后首轮全量）。
  # v2.0.168: fingerprint-based reuse (upgrading the 166 existence tiers).
  # (1) predict present with a matching fingerprint -> reuse (delete
  #     driver.tsv to force a rerun);
  # (2) predict absent but the driver output present with a matching
  #     fingerprint -> vectorized R rebuild (minutes);
  # (3) fingerprint mismatch (peptide table changed) or no sidecar (legacy
  #     artifact) -> full rerun with the stale predict/driver removed so the
  #     downstream cannot read them (one-time cost on the first run after
  #     upgrading).
  sc_ok <- sidecar_matches(cli_args$`out`, input_fp(cli_args$`in`))
  if (file.exists(cli_args$`out`) && sc_ok) {
    cat("[pred_alphapeptdeep] 复用已有产物（指纹匹配；删除 predict.txt.driver.tsv 可强制重跑）：",
        cli_args$`out`, "\n", sep = "")
    return(0L)
  }
  if (file.exists(cli_args$`out`)) {
    cat("[pred_alphapeptdeep] 肽表变化或旧产物无指纹——旧 predict 失效，重新预测\n")
    cat("[pred_alphapeptdeep] stale predict invalidated (peptide table changed or no fingerprint); re-predicting\n")
    unlink(cli_args$`out`)
  }
  need_rebuild <- TRUE
  if (file.exists(drv_out) && sc_ok) {
    cat("[pred_alphapeptdeep] predict 缺失、driver 产物指纹仍匹配——R 侧向量化重建：",
        drv_out, "\n", sep = "")
  } else {
    if (file.exists(drv_out)) {
      cat("[pred_alphapeptdeep] driver 产物同样不可证新鲜——删除重跑\n")
      unlink(drv_out)
    }
    driver <- if (!is.null(cli_args$driver)) cli_args$driver else
      file.path(getwd(), "predictors", "alphapeptdeep_driver.py")
    if (!file.exists(driver)) stop("alphapeptdeep_driver.py 未找到: ", driver)

    a <- c(shQuote(normalizePath(driver)), "--in", shQuote(in_tsv),
           "--out", shQuote(drv_out), "--instrument", instrument, "--ce", ce,
           "--predict-out", shQuote(normalizePath(cli_args$`out`)))
    cat("[pred_alphapeptdeep] 调用 driver：", paste(a, collapse = " "), "\n", sep = "")
    st <- system2(cli_args$python, a, stdout = drv_log, stderr = drv_log)
    if (st != 0L || !file.exists(drv_out)) {
      cat("[pred_alphapeptdeep] driver 失败（exit ", st, "），日志：", drv_log, "\n", sep = "")
      return(1L)
    }
    if (file.exists(cli_args$`out`)) {
      write_run_params(cli_args$`out`, cli_args$`in`)   # v2.0.168：指纹侧车
      cat("[pred_alphapeptdeep] driver 直写完成（输入 ", nrow(input_table),
          " 行）-> ", cli_args$`out`, "（详记 ", drv_log, "）\n", sep = "")
      return(0L)
    }
    cat("[pred_alphapeptdeep] driver 未直写 predict（旧版 driver？），R 侧重建\n")
  }

  if (need_rebuild) {
    drv_df <- read.delim(drv_out, sep = "\t", check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = character(0),
                         colClasses = "character")
    missing_dcols <- setdiff(c("spec_id", "b1", "b2", "y1", "y2"), colnames(drv_df))
    if (length(missing_dcols)) stop("alphapeptdeep driver 输出缺列: ",
                                    paste(missing_dcols, collapse = ","))
    built <- build_predict_txt_fast(input_table, drv_df)
    writeLines(built$lines, cli_args$`out`, useBytes = TRUE)
    write_run_params(cli_args$`out`, cli_args$`in`)   # v2.0.168：指纹侧车
    cat("[pred_alphapeptdeep] 重建完成：", length(built$lines), " 行（跳过：无产物 ",
        built$skip_no_rows, "，长度不符 ", built$skip_len, "）-> ", cli_args$`out`, "\n", sep = "")
  }
  0L
}

if (sys.nframe() == 0L) quit(status = main())
