# =============================================================================
# pred_pdeep3.R —— pDeep3 理论谱图预测器 L1 适配器
#
# 把 isoPSV 的 merge/peptide.txt（peptide / modification / charge 三列）转成
# pDeep3 输入，经 pdeep3_driver.py 调用 pDeep3（pFindStudio/pDeep3，
# python -m pDeep.cmd.tune_and_predict 的 API 形式），再把输出转回与 pDeep2
# predict.txt 同构的 isoPSV 预测表：
#   块格式 BEGIN IONS / TITLE=<seq>|<mod>|<charge> /
#   b+1= / b+2=（正序 b1..b_{n-1}）/ y+1= / y+2=（位点顺序，即写入 rev(y 正序)）
#   / END IONS——与 read_predict_txt 的 y-reverse 语义逐字对齐
# （同 pred_ms2pip.R 的契约；m/z 峰行不需要，消费端只读四类 TYPE 行）。
#
# 用法（由 user_run.R 多预测器循环自动拼装）：
#   Rscript pred_pdeep3.R --in merge/peptide.txt --out merge/predict_pdeep3.txt \
#          --python <python-with-pDeep3> --params <instrument,ce，如 Lumos,30>
#   Rscript pred_pdeep3.R --self-test   # 无需 pDeep3 安装，只测适配层
# 退出码：0 成功；非 0 失败（user_run 据此把该预测列标缺失并继续）。
# 依赖：仅 base R；pDeep3 运行于用户提供的 python 环境。pDeep3 为用户自装独立软件，本适配器不随包分发其本体。
# 状态：适配层自测通过；pDeep3 侧返回结构的强度排序假设（b+1,b+2,y+1,y+2）
# 待 pDeep3 实测环境核验（driver 内有两代返回形态的防御分支与明确报错）。
# =============================================================================

# ---- 参数解析（--flag=value 与 --flag 开关，纯 base R；与 pred_ms2pip 同款）----
parse_cli <- function(args = commandArgs(TRUE)) {
  parsed_args <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (startsWith(token, "--")) {
      flag_key <- sub("^--", "", token)
      if (grepl("=", flag_key, fixed = TRUE)) {
        flag_value_pair <- strsplit(flag_key, "=", fixed = TRUE)[[1]]
        parsed_args[[flag_value_pair[1]]] <- paste(flag_value_pair[-1], collapse = "=")
      } else {
        parsed_args[[flag_key]] <- if (i < length(args) && !startsWith(args[i + 1], "--")) {
          flag_value <- args[i + 1]; i <- i + 1L; flag_value
        } else TRUE
      }
    }
    i <- i + 1L
  }
  parsed_args
}

#' pFind/pDeep 修饰串（"3,Oxidation[M];0,Acetyl[AnyN-term];"）→ pDeep3 modinfo
#' （同族格式，剥尾分号；空串返回 ""——driver 侧按无修饰处理）。
pfind_mods_to_pdeep3 <- function(mod) {
  # 实测教训read.delim 缺省 na.strings 会把空修饰转 NA，
  # TITLE 键随之变 "PEPTIDEK|NA|2"——与消费端 "PEPTIDEK||2" 失配、预测列全缺。
  # 两道保险：读表 na.strings=character(0)＋此处 NA 兜底归 ""。
  if (is.null(mod) || is.na(mod) || !nzchar(mod)) return("")
  sub(";+$", "", trimws(as.character(mod)))
}

#' peptide.txt → pDeep3 输入 TSV 表（spec_id / peptide / modinfo / charge）。
#' orig_mod 保留原串（剥一层尾分号）用于预测块 TITLE 键——与消费端
#' read_predict_txt 的 sub(";+$","",...) 键口径逐字节一致（N-2 同款教训）。
build_pdeep3_table <- function(pep_df) {
  n <- nrow(pep_df)
  data.frame(
    spec_id = sprintf("p%06d", seq_len(n)),
    peptide = as.character(pep_df$peptide),
    modinfo = vapply(seq_len(n), function(i)
      pfind_mods_to_pdeep3(pep_df$modification[i]), character(1)),
    charge = as.character(pep_df$charge),
    orig_mod = sub(";+$", "", as.character(pep_df$modification)),
    stringsAsFactors = FALSE)
}

#' 解析 driver 输出 TSV（spec_id / b1 / b2 / y1 / y2；强度为逗号列表，
#' b、y 均为正序 y1..y_{n-1}）。返回 spec_id → list(b1,b2,y1,y2)。
parse_driver_tsv <- function(tsv_file) {
  raw <- read.delim(tsv_file, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = character(0))
  missing_cols <- setdiff(c("spec_id", "b1", "b2", "y1", "y2"), colnames(raw))
  if (length(missing_cols)) stop("pdeep3 driver 输出缺列: ",
                                 paste(missing_cols, collapse = ","))
  to_num <- function(x) {
    v <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
    v[is.na(v)] <- 0
    v
  }
  out <- vector("list", nrow(raw))
  for (i in seq_len(nrow(raw))) {
    out[[i]] <- list(b1 = to_num(raw$b1[i]), b2 = to_num(raw$b2[i]),
                     y1 = to_num(raw$y1[i]), y2 = to_num(raw$y2[i]))
  }
  setNames(out, raw$spec_id)
}

#' 四条离子向量 → isoPSV predict.txt 一个谱块行组（与 read_predict_txt 对齐：
#' b 正序；y 按位点顺序写入＝rev(y 正序)；无需 m/z 行）。
block_lines <- function(title, b1, b2, y1, y2) {
  c("BEGIN IONS",
    paste0("TITLE=", title),
    paste0("b+1=", paste(format(b1, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("b+2=", paste(format(b2, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("y+1=", paste(rev(format(y1, scientific = FALSE, trim = TRUE)), collapse = ",")),
    paste0("y+2=", paste(rev(format(y2, scientific = FALSE, trim = TRUE)), collapse = ",")),
    "END IONS")
}

#' 输入表 + driver 结果 → predict.txt 行组。离子向量长度须为 n-1；不符或无
#' 产物的行跳过并计数（下游按缺失处理，不阻断）。
#'：组块累积改 list 预分配＋一次 unlist——早期实现逐行 c() 追加是二次方
#' 拷贝，79k 谱（55 万行）实测要拖数分钟；O(n) 后秒级。
#': blocks accumulate into a pre-allocated list + one unlist - the old
#' per-line c() append is quadratic copying (minutes at 79k spectra / 554k
#' lines); O(n) assembly takes seconds.
build_predict_txt <- function(input_table, by_id) {
  blocks <- vector("list", nrow(input_table))
  n_ok <- 0L; skip_no_rows <- 0L; skip_len <- 0L
  for (i in seq_len(nrow(input_table))) {
    ions <- by_id[[input_table$spec_id[i]]]
    if (is.null(ions)) { skip_no_rows <- skip_no_rows + 1L; next }
    residue_count <- nchar(input_table$peptide[i])
    if (length(ions$b1) != residue_count - 1L ||
        length(ions$y1) != residue_count - 1L) { skip_len <- skip_len + 1L; next }
    predict_key <- paste(input_table$peptide[i],
                         input_table$orig_mod[i],
                         input_table$charge[i], sep = "|")
    n_ok <- n_ok + 1L
    blocks[[n_ok]] <- block_lines(predict_key, ions$b1, ions$b2, ions$y1, ions$y2)
  }
  list(lines = if (n_ok > 0L) unlist(blocks[seq_len(n_ok)], use.names = FALSE)
              else character(0),
       skip_no_rows = skip_no_rows, skip_len = skip_len)
}

#' 适配层自测（无需 pDeep3）：合成肽表→构造→假 driver 输出→组块→断言结构。
self_test <- function() {
  cat("[self-test] pred_pdeep3 适配层\n")
  ok <- TRUE
  pep_df <- data.frame(peptide = c("PEPTIDEK", "ALAMAR"),
                       modification = c("2,Oxidation[M];", ""),
                       charge = c(2L, 3L), stringsAsFactors = FALSE)
  tab <- build_pdeep3_table(pep_df)
  stopifnot(tab$modinfo[1] == "2,Oxidation[M]", tab$modinfo[2] == "")
  cat("  [ok] 修饰串转换（剥尾分号/空串）\n")
  # 假 driver 输出：b/y 各 n-1 个正序强度
  drv <- data.frame(spec_id = tab$spec_id,
                    b1 = c("1,2,3,4,5,6,7", "1,2,3,4,5"),
                    b2 = c("0.1,0.2,0.3,0.4,0.5,0.6,0.7", "0.1,0.2,0.3,0.4,0.5"),
                    y1 = c("7,6,5,4,3,2,1", "5,4,3,2,1"),
                    y2 = c("0.7,0.6,0.5,0.4,0.3,0.2,0.1", "0.5,0.4,0.3,0.2,0.1"),
                    stringsAsFactors = FALSE)
  tf1 <- tempfile(fileext = ".tsv"); tf2 <- tempfile(fileext = ".tsv")
  write.table(drv, tf1, sep = "\t", row.names = FALSE, quote = FALSE)
  by_id <- parse_driver_tsv(tf1)
  built <- build_predict_txt(tab, by_id)
  lines <- built$lines
  stopifnot(built$skip_no_rows == 0L, built$skip_len == 0L,
            length(lines) == 2L * 7L,
            sum(lines == "BEGIN IONS") == 2L, sum(lines == "END IONS") == 2L,
            lines[2] == "TITLE=PEPTIDEK|2,Oxidation[M]|2",
            lines[3] == "b+1=1,2,3,4,5,6,7",
            lines[5] == "y+1=1,2,3,4,5,6,7")   # y 正序 7..1 → 位点序写出 1..7
  cat("  [ok] 块结构/TITLE 键/y 位点序反转（与 read_predict_txt 对齐）\n")
  # 键口径回验：消费端 sub(\";+$\",\"\") 后应逐字节等于 TITLE
  key2 <- paste(tab$peptide[2], sub(";+$", "", pep_df$modification[2]),
                pep_df$charge[2], sep = "|")
  stopifnot(grepl(paste0("TITLE=", key2, "$"), lines, fixed = TRUE)[1] ||
              any(lines == paste0("TITLE=", key2)))
  cat("  [ok] TITLE 键与消费端口径一致\n")
  cat("[self-test] 全部通过\n")
  TRUE
}

# ---- 主入口 ----
main <- function() {
  cli_args <- parse_cli()
  if (isTRUE(cli_args$`self-test`)) return(if (self_test()) 0L else 1L)
  required_args <- c("in", "out", "python")
  missing_args <- setdiff(required_args, names(cli_args))
  if (length(missing_args)) stop("缺参数: ", paste(paste0("--", missing_args), collapse = " "))

  pep_df <- read.delim(cli_args$`in`, sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE, na.strings = character(0))
  # 空修饰列会被类型推断成 logical NA（na.strings 拦不住）——统一归空串，
  # 保 TITLE 键 "PEPTIDEK||2" 与消费端口径一致
  pep_df$modification[is.na(pep_df$modification)] <- ""
  missing_cols <- setdiff(c("peptide", "modification", "charge"), colnames(pep_df))
  if (length(missing_cols)) stop("peptide.txt 缺列: ", paste(missing_cols, collapse = ","))

  input_table <- build_pdeep3_table(pep_df)
  #：中间件入预测器独立夹（与 --out 同目录），peptides.tsv＝pDeep3 专属
  # 输入列表（spec_id/peptide/modinfo/charge 格式，与 pDeep2 的 peptides.txt
  # 同源不同形）——目录不存在时（旧布局调用）先建。
  #: intermediates live beside --out in the predictor folder;
  # peptides.tsv is the pDeep3-specific input list (same source, different
  # format than pDeep2's peptides.txt). Folder is created when missing.
  dir.create(dirname(cli_args$`out`), recursive = TRUE, showWarnings = FALSE)
  in_tsv <- file.path(dirname(cli_args$`out`), "peptides.tsv")
  write.table(input_table[, c("spec_id", "peptide", "modinfo", "charge")],
              in_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  if (isTRUE(cli_args$`dry-run`)) {
    cat("dry-run: 输入 TSV 已写", in_tsv, "（", nrow(input_table), "行），未调用 pDeep3\n")
    return(0L)
  }

  # --params：instrument,ce（缺省 Lumos,30）。pDeep3 的 EThcD 模型须
  # instrument=Lumos＋指定 ce；HCD 模型按实际仪器填（QE/Lumos）。
  params_raw <- if (is.null(cli_args$params)) "" else as.character(cli_args$params)
  params_parts <- strsplit(params_raw, ",", fixed = TRUE)[[1]]
  instrument <- if (length(params_parts) >= 1L && nzchar(trimws(params_parts[1])))
    trimws(params_parts[1]) else "Lumos"
  ce <- if (length(params_parts) >= 2L && nzchar(trimws(params_parts[2])))
    trimws(params_parts[2]) else "30"

  drv_out <- paste0(cli_args$`out`, ".driver.tsv")
  drv_log <- paste0(cli_args$`out`, ".driver.log")
  # 适配器与 driver 同目录（code/predictors/）。user_run 在 code/ 下以
  # Rscript predictors/pred_<type>.R 调用（getwd()＝code/），故默认相对
  # 解析；非该布局时由调用方显式传 --driver。
  driver <- if (!is.null(cli_args$driver)) cli_args$driver else
    file.path(getwd(), "predictors", "pdeep3_driver.py")
  if (!file.exists(driver)) stop("pdeep3_driver.py 未找到: ", driver)

  a <- c(shQuote(normalizePath(driver)), "--in", shQuote(in_tsv),
         "--out", shQuote(drv_out), "--instrument", instrument, "--ce", ce)
  cat("[pred_pdeep3] 调用 driver：", paste(a, collapse = " "), "\n", sep = "")
  st <- system2(cli_args$python, a, stdout = drv_log, stderr = drv_log)
  if (st != 0L || !file.exists(drv_out)) {
    cat("[pred_pdeep3] driver 失败（exit ", st, "），日志：", drv_log, "\n", sep = "")
    return(1L)
  }

  by_id <- parse_driver_tsv(drv_out)
  built <- build_predict_txt(input_table, by_id)
  writeLines(built$lines, cli_args$`out`, useBytes = TRUE)
  cat("[pred_pdeep3] 完成：", length(built$lines), " 行（跳过：无产物 ",
      built$skip_no_rows, "，长度不符 ", built$skip_len, "）-> ", cli_args$`out`, "\n", sep = "")
  0L
}

if (sys.nframe() == 0L) quit(status = main())
