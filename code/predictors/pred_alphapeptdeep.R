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

parse_driver_tsv <- function(tsv_file) {
  raw <- read.delim(tsv_file, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = character(0))
  missing_cols <- setdiff(c("spec_id", "b1", "b2", "y1", "y2"), colnames(raw))
  if (length(missing_cols)) stop("alphapeptdeep driver 输出缺列: ",
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

block_lines <- function(title, b1, b2, y1, y2) {
  c("BEGIN IONS",
    paste0("TITLE=", title),
    paste0("b+1=", paste(format(b1, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("b+2=", paste(format(b2, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("y+1=", paste(rev(format(y1, scientific = FALSE, trim = TRUE)), collapse = ",")),
    paste0("y+2=", paste(rev(format(y2, scientific = FALSE, trim = TRUE)), collapse = ",")),
    "END IONS")
}

build_predict_txt <- function(input_table, by_id) {
  predict_lines <- character(0)
  skip_no_rows <- 0L; skip_len <- 0L
  for (i in seq_len(nrow(input_table))) {
    ions <- by_id[[input_table$spec_id[i]]]
    if (is.null(ions)) { skip_no_rows <- skip_no_rows + 1L; next }
    residue_count <- nchar(input_table$peptide[i])
    if (length(ions$b1) != residue_count - 1L ||
        length(ions$y1) != residue_count - 1L) { skip_len <- skip_len + 1L; next }
    predict_key <- paste(input_table$peptide[i],
                         input_table$orig_mod[i],
                         input_table$charge[i], sep = "|")
    predict_lines <- c(predict_lines,
                       block_lines(predict_key, ions$b1, ions$b2, ions$y1, ions$y2))
  }
  list(lines = predict_lines, skip_no_rows = skip_no_rows, skip_len = skip_len)
}

self_test <- function() {
  cat("[self-test] pred_alphapeptdeep 适配层\n")
  pep_df <- data.frame(peptide = c("PEPTIDEK", "ALAMAR"),
                       modification = c("2,Oxidation[M];", ""),
                       charge = c(2L, 3L), stringsAsFactors = FALSE)
  tab <- build_table(pep_df)
  stopifnot(tab$modinfo[1] == "2,Oxidation[M]", tab$modinfo[2] == "")
  cat("  [ok] 修饰串转换（剥尾分号/空串）\n")
  drv <- data.frame(spec_id = tab$spec_id,
                    b1 = c("1,2,3,4,5,6,7", "1,2,3,4,5"),
                    b2 = c("0.1,0.2,0.3,0.4,0.5,0.6,0.7", "0.1,0.2,0.3,0.4,0.5"),
                    y1 = c("7,6,5,4,3,2,1", "5,4,3,2,1"),
                    y2 = c("0.7,0.6,0.5,0.4,0.3,0.2,0.1", "0.5,0.4,0.3,0.2,0.1"),
                    stringsAsFactors = FALSE)
  tf1 <- tempfile(fileext = ".tsv")
  write.table(drv, tf1, sep = "\t", row.names = FALSE, quote = FALSE)
  by_id <- parse_driver_tsv(tf1)
  built <- build_predict_txt(tab, by_id)
  lines <- built$lines
  stopifnot(built$skip_no_rows == 0L, built$skip_len == 0L,
            length(lines) == 2L * 7L,
            sum(lines == "BEGIN IONS") == 2L, sum(lines == "END IONS") == 2L,
            lines[2] == "TITLE=PEPTIDEK|2,Oxidation[M]|2",
            lines[3] == "b+1=1,2,3,4,5,6,7",
            lines[5] == "y+1=1,2,3,4,5,6,7")
  cat("  [ok] 块结构/TITLE 键/y 位点序反转\n")
  key2 <- paste(tab$peptide[2], sub(";+$", "", pep_df$modification[2]),
                pep_df$charge[2], sep = "|")
  stopifnot(any(lines == paste0("TITLE=", key2)))
  cat("  [ok] TITLE 键与消费端口径一致\n")
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

  input_table <- build_table(pep_df)
  in_tsv <- paste0(cli_args$`out`, ".apd.tsv")
  write.table(input_table[, c("spec_id", "peptide", "modinfo", "charge")],
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
  driver <- if (!is.null(cli_args$driver)) cli_args$driver else
    file.path(getwd(), "predictors", "alphapeptdeep_driver.py")
  if (!file.exists(driver)) stop("alphapeptdeep_driver.py 未找到: ", driver)

  a <- c(shQuote(normalizePath(driver)), "--in", shQuote(in_tsv),
         "--out", shQuote(drv_out), "--instrument", instrument, "--ce", ce)
  cat("[pred_alphapeptdeep] 调用 driver：", paste(a, collapse = " "), "\n", sep = "")
  st <- system2(cli_args$python, a, stdout = drv_log, stderr = drv_log)
  if (st != 0L || !file.exists(drv_out)) {
    cat("[pred_alphapeptdeep] driver 失败（exit ", st, "），日志：", drv_log, "\n", sep = "")
    return(1L)
  }

  by_id <- parse_driver_tsv(drv_out)
  built <- build_predict_txt(input_table, by_id)
  writeLines(built$lines, cli_args$`out`, useBytes = TRUE)
  cat("[pred_alphapeptdeep] 完成：", length(built$lines), " 行（跳过：无产物 ",
      built$skip_no_rows, "，长度不符 ", built$skip_len, "）-> ", cli_args$`out`, "\n", sep = "")
  0L
}

if (sys.nframe() == 0L) quit(status = main())
