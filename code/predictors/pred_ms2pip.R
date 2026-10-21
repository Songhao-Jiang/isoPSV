# =============================================================================
# pred_ms2pip.R —— MS2PIP 理论谱图预测器 L1 适配器
#
# 多理论谱图预测器接入（deepseek 总体设计）的第一个外部适配器。把 isoPSV 的
# merge/peptide.txt（peptide / modification / charge 三列）转成 MS2PIP peprec，
# 调用 MS2PIP 预测，再把输出转回与 pDeep2 predict.txt 同构的 isoPSV 预测表：
#   块格式 BEGIN IONS / TITLE=<seq>|<mod>|<charge> / b+1= / b+2= / y+1= / y+2= / END IONS
#   b 行正序（b1..b_{n-1}）、y 行按"位点顺序"存（y_{n-1}..y1）——与
#   read_predict_txt 的 reverse 语义一致（y 向量写入时需反转）。
#
# 用法：
#   Rscript pred_ms2pip.R --in merge/peptide.txt --out merge/predict_ms2pip.txt \
#          --python <python-with-ms2pip> --params <ms2pip 模型目录> [--threads 8]
#   Rscript pred_ms2pip.R --self-test   # 无需 ms2pip 安装，只测适配层
# 退出码：0 成功；非 0 失败（user_run 据此把该预测列标缺失并继续）。
# 依赖：仅 base R；MS2PIP 运行于用户提供的 python 环境（ms2pip>=3，含模型）。
# =============================================================================

# ---- 参数解析（--flag=value 与 --flag 开关，纯 base R）----
parse_cli <- function(args = commandArgs(TRUE)) {
  # parsed_args：解析结果（flag 名 → 值/TRUE）；token：当前命令行片段；
  # flag_key：去掉 -- 的 flag 名；flag_value_pair：--k=v 的键值对
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

# ---- pFind 修饰名 → MS2PIP/Unimod 名 ----
mod_name_map <- c(Oxidation = "Oxidation", Deamidated = "Deamidated",
                  Carbamidomethyl = "Carbamidomethyl", Acetyl = "Acetyl",
                  Phospho = "Phospho", Methyl = "Methyl", Dimethyl = "Dimethyl",
                  Trimethyl = "Trimethyl", `Gln->pyro-Glu` = "Gln->pyro-Glu",
                  `Glu->pyro-Glu` = "Glu->pyro-Glu")
#' pFind 修饰串（"3,Oxidation[M];0,Acetyl[AnyN-term];"）→ peprec modifications 串
#' （"3|Oxidation;0|Acetyl"）。未知名丢弃并警告（与 pDeep2 路径的 sanitize 同理）。
#' 无修饰返回 "-"（N-11：MS2PIP peprec 对未修饰肽要求写 "-"，空串会解析失败）。
pfind_mods_to_peprec <- function(mod, warn = TRUE) {
  if (is.null(mod) || !nzchar(mod)) return("-")
  mod_tokens <- strsplit(mod, ";", fixed = TRUE)[[1]]
  mod_tokens <- mod_tokens[nzchar(mod_tokens)]
  mapped_mods <- character(0)
  for (mod_token in mod_tokens) {
    site_position <- sub("^([0-9]+),.*$", "\\1", mod_token)
    full_name <- sub("^[0-9]+,(.*)$", "\\1", mod_token)
    mod_name <- sub("\\[[^]]*\\]$", "", full_name)   # 去 [site]/[AnyN-term]
    if (mod_name %in% names(mod_name_map)) {
      mapped_mods <- c(mapped_mods, paste0(site_position, "|", unname(mod_name_map[mod_name])))
    } else if (warn) {
      cat("未知修饰名（跳过）:", mod_token, "\n")
    }
  }
  #：多修饰以 "|" 扁平连接（ms2pip 4.x/psm_utils 的 PEPREC 语法＝
  # 位置|名称|位置|名称…；分号连接会被当成单个标签名解析失败）。位置 0＝N 端
  # （psm_utils peprec_to_proforma 的 peptide_list[0] 即 N-term）。
  #: join multiple mods with "|" (ms2pip 4.x / psm_utils PEPREC
  # grammar = flat position|name|position|name...; the previous "; " join was
  # parsed as one giant label and failed). Position 0 = N-terminus.
  if (length(mapped_mods)) paste(mapped_mods, collapse = "|") else "-"
}

#' peptide.txt → MS2PIP peprec 数据框（spec_id / modifications / peptide / charge /
#' orig_mod）。N-2orig_mod 保留 peptide.txt 修饰原串（尾分号
#' 剥掉一层，与消费端 read_predict_txt 的 sub(";+$","",...) 完全同口径），预测块
#' TITLE 键用它拼接——保证与下游键逐字节一致（写 peprec 文件时不含此列）。
build_peprec <- function(pep_df) {
  # pep_df：peptide.txt 三列表（peptide/modification/charge）；返回 peprec 表
  n <- nrow(pep_df)
  data.frame(
    spec_id = sprintf("p%06d", seq_len(n)),
    modifications = vapply(seq_len(n), function(i)
      pfind_mods_to_peprec(as.character(pep_df$modification[i])), character(1)),
    peptide = as.character(pep_df$peptide),
    charge = as.character(pep_df$charge),
    orig_mod = sub(";+$", "", as.character(pep_df$modification)),
    stringsAsFactors = FALSE)
}

#' 解析 MS2PIP 预测 CSV（列：spec_id, charge, ion, b_ion, b_ion_mz, b_int,
#' y_ion, y_ion_mz, y_int；每个 ion 行按电荷 1..z 各一行，z=母离子电荷）。
#' 返回 list(by_id=按 spec_id 分组的命名 list(seq, b1, b2, y1, y2),
#' charge3plus_ids=含电荷 3+ 行的 spec_id 数)——y1/y2 为 y1..y_{n-1} 正序。
#' N-12原实现对每个 spec_id 全表扫描 O(ids×rows)，
#' 改为一次 split()；分组顺序与 unique() 的首现顺序一致。
#' P0-3isoPSV 预测块只承载 b+1/b+2/y+1/y+2（与
#' pDeep2 同口径），电荷 3+ 行的强度不进入预测表——但必须计数上报，不得静默。
parse_ms2pip_csv <- function(csv_file) {
  raw_csv <- read.delim(csv_file, sep = ",", check.names = FALSE,
                        stringsAsFactors = FALSE, na.strings = character(0))
  required_cols <- c("spec_id", "charge", "ion", "b_int", "y_int")
  missing_cols <- setdiff(required_cols, colnames(raw_csv))
  if (length(missing_cols)) stop("MS2PIP 输出缺列: ", paste(missing_cols, collapse = ","))
  raw_csv$charge <- as.integer(as.character(raw_csv$charge))
  raw_csv$ion <- as.integer(as.character(raw_csv$ion))
  raw_csv$b_int <- suppressWarnings(as.numeric(as.character(raw_csv$b_int)))
  raw_csv$y_int <- suppressWarnings(as.numeric(as.character(raw_csv$y_int)))
  raw_csv$b_int[is.na(raw_csv$b_int)] <- 0
  raw_csv$y_int[is.na(raw_csv$y_int)] <- 0
  # spec_groups：spec_id → 行子集（split 保持首现顺序，与 unique 同序）
  spec_groups <- split(raw_csv, raw_csv$spec_id)
  # charge3plus_ids：含电荷 3+ 行的 spec_id（上报用，不阻断；P3-4 加
  # na.rm——畸形 CSV 的 charge 空值 as.integer 为 NA，any(NA)=NA 会把
  # sum 传染成 NA，下游 if(NA>0) 崩）
  charge3plus_ids <- sum(vapply(spec_groups, function(s)
    any(s$charge >= 3L, na.rm = TRUE), logical(1)))
  by_id <- lapply(spec_groups, function(spec_rows) {
    # spec_rows：该 spec_id 的行；row_order：按电荷+离子号排序；
    # charge1_rows/charge2_rows：单/双电荷子集；charge1_sorted/charge2_sorted：
    # 按 ion 号排序后的子集（b_int/y_int 即 b1.. / y1.. 正序）
    row_order <- order(spec_rows$charge, spec_rows$ion)
    spec_rows <- spec_rows[row_order, ]
    charge1_rows <- spec_rows[spec_rows$charge == 1L, ]
    charge2_rows <- spec_rows[spec_rows$charge == 2L, ]
    charge1_sorted <- charge1_rows[order(charge1_rows$ion), ]
    charge2_sorted <- charge2_rows[order(charge2_rows$ion), ]
    list(b1 = charge1_sorted$b_int, b2 = charge2_sorted$b_int,
         y1 = charge1_sorted$y_int, y2 = charge2_sorted$y_int)
  })
  list(by_id = by_id, charge3plus_ids = charge3plus_ids)
}

#' 把四条离子向量写成 isoPSV predict.txt 的一个谱块行组。
#' b 行正序；y 行按 read_predict_txt 的期望存"位点顺序"（rev(y 正序)）。
block_lines <- function(title, b1, b2, y1, y2) {
  c("BEGIN IONS",
    paste0("TITLE=", title),
    paste0("b+1=", paste(format(b1, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("b+2=", paste(format(b2, scientific = FALSE, trim = TRUE), collapse = ",")),
    paste0("y+1=", paste(rev(format(y1, scientific = FALSE, trim = TRUE)), collapse = ",")),
    paste0("y+2=", paste(rev(format(y2, scientific = FALSE, trim = TRUE)), collapse = ",")),
    "END IONS")
}

#' peprec + 解析结果 → isoPSV predict.txt 内容。离子向量长度须为 n-1；
#' 长度不符或无产物的行跳过并按计数上报（下游按缺失处理）。
#' 返回 list(lines=块行组, skip_no_rows=无 CSV 产物数, skip_len=离子数不符数)。
build_predict_txt <- function(peprec_table, parsed_by_id) {
  # peprec_table：peprec 表；parsed_by_id：spec_id → 离子强度列表；
  # parsed_ions：当前肽的离子强度；predict_key：isoPSV 预测表 TITLE 键
  # （N-2：用 orig_mod 原串拼接，与消费端逐字节一致）；residue_count：肽长
  # （各离子向量须为 n-1，不符的行跳过按缺失处理）
  predict_lines <- character(0)
  skip_no_rows <- 0L; skip_len <- 0L
  for (i in seq_len(nrow(peprec_table))) {
    spec_id0 <- peprec_table$spec_id[i]
    parsed_ions <- parsed_by_id[[spec_id0]]
    if (is.null(parsed_ions)) { skip_no_rows <- skip_no_rows + 1L; next }
    residue_count <- nchar(peprec_table$peptide[i])
    if (length(parsed_ions$b1) != residue_count - 1L ||
        length(parsed_ions$y1) != residue_count - 1L) { skip_len <- skip_len + 1L; next }
    predict_key <- paste(peprec_table$peptide[i],
                         peprec_table$orig_mod[i],
                         peprec_table$charge[i], sep = "|")
    predict_lines <- c(predict_lines,
                       block_lines(predict_key, parsed_ions$b1, parsed_ions$b2,
                                   parsed_ions$y1, parsed_ions$y2))
  }
  list(lines = predict_lines, skip_no_rows = skip_no_rows, skip_len = skip_len)
}

#' MS2PIP v4 现代契约的长表解析。
#' v4.0+ 正式版 predict-batch -f tsv 输出：psm_index/ion_type/ion_number/
#' mz/predicted（log2 强度、单电荷 b/y）。b2/y2 以 0 填充（v4 批量输出无
#' 电荷维度）；psm_index 0 起按首现顺序对齐 peprec 的 p000001 起 spec_id。
parse_ms2pip_long <- function(tsv_file) {
  raw <- read.delim(tsv_file, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = character(0))
  missing_cols <- setdiff(c("psm_index", "ion_type", "ion_number", "predicted"),
                          colnames(raw))
  if (length(missing_cols)) stop("MS2PIP v4 长表缺列: ",
                                 paste(missing_cols, collapse = ","))
  raw$predicted <- suppressWarnings(as.numeric(as.character(raw$predicted)))
  raw$predicted[is.na(raw$predicted)] <- 0
  raw$intensity <- 2^raw$predicted
  grp_order <- as.integer(unique(raw$psm_index))
  grp_order <- sort(grp_order)
  by_id <- lapply(grp_order, function(idx) {
    g <- raw[raw$psm_index == idx,, drop = FALSE]
    n_frag <- max(g$ion_number)
    mk <- function(ion) {
      v <- numeric(n_frag)
      sub_g <- g[g$ion_type == ion,, drop = FALSE]
      v[sub_g$ion_number] <- sub_g$intensity
      v
    }
    list(b1 = mk("b"), b2 = numeric(n_frag),
         y1 = mk("y"), y2 = numeric(n_frag))
  })
  stats::setNames(by_id, sprintf("p%06d", seq_along(grp_order)))
}

# ---- 主入口 ----
main <- function() {
  cli_args <- parse_cli()
  if (isTRUE(cli_args$`self-test`)) return(if (self_test()) 0L else 1L)
  #：--params 改为可选——网页的 MS2PIP 行按设计留空（默认取适配器
  # 内置 HCD 模型），此前把它列入必填，导致"网页填参数→跑"这条路径上 ms2pip
  # 恒失败（exit 1，缺 --params），该列整列丢失。空/缺省一律按 HCD 处理。
  #: --params is optional - the web form intentionally leaves the
  # MS2PIP row empty (default = the adapter's built-in HCD model). Requiring it
  # made ms2pip fail on every form-driven run (exit 1, missing --params) and
  # dropped the whole column. Missing/empty now means HCD.
  required_args <- c("in", "out", "python")
  missing_args <- setdiff(required_args, names(cli_args))
  if (length(missing_args)) stop("缺参数: ", paste(paste0("--", missing_args), collapse = " "))
  if (is.null(cli_args$params)) cli_args$params <- "HCD"
  if (!nzchar(cli_args$params)) cli_args$params <- "HCD"

  pep_df <- read.delim(cli_args$`in`, sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE)
  required_cols2 <- c("peptide", "modification", "charge")
  missing_cols2 <- setdiff(required_cols2, colnames(pep_df))
  if (length(missing_cols2)) stop("peptide.txt 缺列: ", paste(missing_cols2, collapse = ","))
  #：modification 列全空时 read.delim 推断为 logical 全 NA（与
  # merge_search_results 的列类型守卫同源坑）——orig_mod 拼出 "seq|NA|z" 键，
  # 下游 seq||z 查键全部落空、ms2pip 列整列丢失。此处统一 NA/字面 "NA" → ""。
  #: an all-empty modification column is inferred as logical NA by
  # read.delim (the same trap as the column-type guard in merge_search_results)
  # - orig_mod then builds "seq|NA|z" keys and every downstream seq||z lookup
  # misses, dropping the whole ms2pip column. Normalize NA / literal "NA" here.
  pep_df$modification[is.na(pep_df$modification)] <- ""
  pep_df$modification[pep_df$modification == "NA"] <- ""
  threads <- if (!is.null(cli_args$threads)) as.integer(cli_args$threads) else 8L
  if (is.na(threads) || threads < 1) threads <- 8L

  peprec_table <- build_peprec(pep_df)
  peprec_file <- paste0(cli_args$`out`, ".peprec.tsv")
  # orig_mod 只用于本适配器的键拼接，不进 peprec 文件（MS2PIP 只认四列）
  write.table(peprec_table[, c("spec_id", "modifications", "peptide", "charge")],
              peprec_file, sep = "\t", row.names = FALSE, quote = FALSE)
  ms2pip_csv_file <- paste0(cli_args$`out`, ".ms2pip.csv")

  if (isTRUE(cli_args$`dry-run`)) {
    cat("dry-run: peprec 已写", peprec_file, "（", nrow(peprec_table), "行），未调用 ms2pip\n")
    return(0L)
  }
  # N-11CLI 契约按官方两代——4.x 用 predict-batch +
  # --model-dir/-o；3.x 无 predict 子命令（-c 配置 + -n 并行）。--threads 两代
  # 都不存在。先试 4.x，失败自动回退 3.x；两代都不行则报错并留两份日志。
  # P2-2两代日志分开——此前共用一条
  # 路径，3.x 运行覆盖 4.x 失败证据，恰好丢失最需要诊断的那半。
  # ms2pip_log_4x/3x：两代独立日志路径；cli_contract_used：实际生效的契约
  ms2pip_log_4x <- paste0(cli_args$`out`, ".ms2pip.log")
  ms2pip_log_3x <- paste0(cli_args$`out`, ".ms2pip.3x.log")
  run_ms2pip <- function(arg_list, log_path) {
    suppressWarnings(system2(cli_args$python, arg_list,
      stdout = log_path, stderr = log_path))
  }
  # v4 现代契约位置参数 PSMS ＋ -t peprec ＋
  # --model <名>（params 非目录时作模型名，默认 HCD）＋ -f tsv；-o 为输出名
  # 前缀（实际写 <prefix>.tsv）。官方 4.0-4.2 无 --peprec 旗标（原 4.x 契约
  # 系预发布口径，保留于下方作回退）。
  # v4 -o semantics = Path.with_suffix(".tsv") (replaces the last suffix)
  ms2pip_tsv_file <- sub("\\.[^.]+$", ".tsv", cli_args$`out`)
  if (ms2pip_tsv_file == cli_args$`out`) ms2pip_tsv_file <- paste0(cli_args$`out`, ".tsv")
  ms2pip_log_v4 <- paste0(cli_args$`out`, ".ms2pip.v4.log")
  model_name <- if (!is.null(cli_args$params) && nzchar(cli_args$params) &&
                    !dir.exists(cli_args$params)) cli_args$params else "HCD"
  unlink(c(ms2pip_tsv_file, ms2pip_log_v4))
  stv4 <- run_ms2pip(c("-m", "ms2pip", "predict-batch", "-t", "peprec",
                       "--model", model_name, "-f", "tsv",
                       "-o", cli_args$`out`,
                       peprec_file), ms2pip_log_v4)
  if (stv4 == 0 && file.exists(ms2pip_tsv_file)) {
    long_by_id <- parse_ms2pip_long(ms2pip_tsv_file)
    predict_result <- build_predict_txt(peprec_table, long_by_id)
    if (!length(predict_result$lines))
      stop("MS2PIP v4 结果与 peptide.txt 无匹配（肽长/离子数不符）")
    writeLines(predict_result$lines, cli_args$`out`, useBytes = TRUE)
    cat("[MS2PIP] v4 现代契约完成（b+2/y+2 单电荷无维度以 0 填充）；跳过：无产物 ",
        predict_result$skip_no_rows, "，长度不符 ", predict_result$skip_len,
        " -> ", cli_args$`out`, "\n", sep = "")
    return(0L)
  }
  cat("[MS2PIP] v4 现代契约失败（exit=", stv4, "），回退旧契约；日志 ",
      ms2pip_log_v4, "\n", sep = "")
  # 4.x：python -m ms2pip predict-batch --model-dir <模型目录> --peprec <tsv> -o <csv>
  unlink(c(ms2pip_csv_file, ms2pip_log_3x))   # 先清上一轮残留，防失败后误判"产物已存在"
  st4 <- run_ms2pip(c("-m", "ms2pip", "predict-batch",
                      "--model-dir", cli_args$params, "--peprec", peprec_file,
                      "-o", ms2pip_csv_file), ms2pip_log_4x)
  cli_contract_used <- "4x"
  if (st4 != 0 || !file.exists(ms2pip_csv_file)) {
    cat("[MS2PIP] 4.x 契约失败（exit=", st4, "），回退 3.x（-c 配置 + -n 并行）\n", sep = "")
    # 3.x：config 优先取 params 目录下的 config.txt；params 本身是文件则直接用
    config3 <- if (dir.exists(cli_args$params) && file.exists(file.path(cli_args$params, "config.txt"))) {
      file.path(cli_args$params, "config.txt")
    } else cli_args$params
    unlink(ms2pip_csv_file)
    st3 <- run_ms2pip(c("-m", "ms2pip", "-c", config3, peprec_file, "-n", as.character(threads)),
                      ms2pip_log_3x)
    cli_contract_used <- "3x"
    if (st3 != 0 || !file.exists(ms2pip_csv_file))
      stop("ms2pip 两代契约均未产出预测 CSV（4.x exit=", st4, "，见 ",
           ms2pip_log_4x, "；3.x exit=", st3, "，见 ", ms2pip_log_3x, "）")
  }
  parsed_result <- parse_ms2pip_csv(ms2pip_csv_file)
  predict_result <- build_predict_txt(peprec_table, parsed_result$by_id)
  predict_lines <- predict_result$lines
  if (!length(predict_lines)) stop("MS2PIP 结果与 peptide.txt 无匹配（肽长/离子数不符）")
  writeLines(predict_lines, cli_args$`out`, useBytes = TRUE)
  # P0-3电荷 3+ 行的强度不进入预测表（表格式与 pDeep2
  # 同口径只含 +1/+2）——计数上报，绝不静默；无产物/长度不符的跳过数同样上报
  if (parsed_result$charge3plus_ids > 0L || predict_result$skip_no_rows > 0L ||
      predict_result$skip_len > 0L) {
    cat("[MS2PIP] 跳过计数：电荷3+（仅取+1/+2 行）", parsed_result$charge3plus_ids,
        " 条；无 CSV 产物 ", predict_result$skip_no_rows,
        " 条；离子数不符 ", predict_result$skip_len, " 条\n", sep = "")
  }
  writeLines(c(
    paste0("tool=ms2pip"),
    paste0("cli_contract=", cli_contract_used),
    paste0("python=", cli_args$python),
    paste0("params=", cli_args$params),
    paste0("peptides=", nrow(peprec_table)),
    paste0("written=", length(predict_lines) / 7L),
    paste0("charge3plus=", parsed_result$charge3plus_ids),
    paste0("time=", format(Sys.time()))),
    paste0(cli_args$`out`, ".meta"), useBytes = TRUE)
  cat("MS2PIP 预测完成（契约 ", cli_contract_used, "）：", length(predict_lines) / 7L,
      " 条 →", cli_args$`out`, "\n", sep = "")
  0L
}

# ---- 自检（无需 ms2pip 安装；只测适配层纯函数）----
self_test <- function() {
  ok <- TRUE
  say <- function(pass, what, detail = "") {
    cat(sprintf("  [%s] %-46s %s\n", if (isTRUE(pass)) "PASS" else "FAIL", what, detail))
    if (!isTRUE(pass)) ok <<- FALSE
  }
  cat("== pred_ms2pip 自检 ==\n")
  # 1) pFind 修饰 → peprec
  say(identical(pfind_mods_to_peprec("3,Oxidation[M];0,Acetyl[AnyN-term];"), "3|Oxidation;0|Acetyl"),
      "pFind 修饰转 peprec", pfind_mods_to_peprec("3,Oxidation[M];0,Acetyl[AnyN-term];"))
  # 2) peprec 构造（含未知修饰跳过）
  df <- data.frame(peptide = "PEPTIDEK", modification = "4,Oxidation[M];9,Unknown[X];",
                   charge = "2", stringsAsFactors = FALSE)
  pr <- build_peprec(df)
  say(nrow(pr) == 1 && pr$peptide[1] == "PEPTIDEK" && pr$charge[1] == "2",
      "peprec 构造", paste(pr[1, ], collapse = " | "))
  # 3) 块布局：b 正序、y 反转（与 read_predict_txt 的 reverse 语义一致）
  bl <- block_lines("PEPTIDEK||2", 1:7, 11:17, 21:27, 31:37)
  y1 <- bl[which(startsWith(bl, "y+1="))]
  y1v <- as.numeric(strsplit(sub("y+1=", "", y1, fixed = TRUE), ",", fixed = TRUE)[[1]])
  say(identical(y1v, as.numeric(27:21)), "y 行按位点顺序（写入时反转）", paste(range(y1v), collapse = ".."))
  b1 <- as.numeric(strsplit(sub("b+1=", "", bl[startsWith(bl, "b+1=")][1], fixed = TRUE),
                             ",", fixed = TRUE)[[1]])
  say(identical(b1, as.numeric(1:7)), "b 行正序")
  # 4) 合成 CSV 解析 + 组装（模拟 ms2pip 输出；含电荷 3+ 行的夹具）
  td <- tempfile("pfm_"); dir.create(td)
  csv <- file.path(td, "pred.csv")
  writeLines(c("spec_id,charge,ion,b_ion,b_ion_mz,b_int,y_ion,y_ion_mz,y_int",
               "p000001,1,1,b1,100,10,y1,700,20",
               "p000001,1,2,b2,200,30,y2,800,40",
               "p000001,1,3,b3,300,50,y3,900,60",
               "p000001,2,1,b1,100,11,y1,700,21",
               "p000001,2,2,b2,200,31,y2,800,41",
               "p000001,2,3,b3,300,51,y3,900,61",
               "p000001,3,1,b1,100,12,y1,700,22",
               "p000001,3,2,b2,200,32,y2,800,42",
               "p000001,3,3,b3,300,52,y3,900,62"), csv)
  parsed <- parse_ms2pip_csv(csv)
  say(identical(names(parsed$by_id), "p000001") && parsed$charge3plus_ids == 1L,
      "解析含电荷3+行并计数上报", paste("charge3plus=", parsed$charge3plus_ids))
  peprec <- data.frame(spec_id = "p000001", modifications = "", peptide = "PEPT",
                       charge = "3", orig_mod = "", stringsAsFactors = FALSE)
  built <- build_predict_txt(peprec, parsed$by_id)
  lines <- built$lines
  say(length(lines) == 7 && startsWith(lines[2], "TITLE=PEPT||3"),
      "合成 CSV → isoPSV 预测块", paste(length(lines), "行"))
  b2 <- as.numeric(strsplit(sub("b+2=", "", lines[4], fixed = TRUE), ",", fixed = TRUE)[[1]])
  say(identical(b2, as.numeric(c(11, 31, 51))), "b+2 取电荷 2 行（电荷3+行不混入）")
  y2 <- as.numeric(strsplit(sub("y+2=", "", lines[6], fixed = TRUE), ",", fixed = TRUE)[[1]])
  say(identical(y2, as.numeric(c(61, 41, 21))), "y+2 反转后为位点顺序")
  # 5) N-2 键口径：TITLE 用肽表修饰原串（保留 [M]、去尾分号），与消费端
  #    read_predict_txt 的 sub(";+$","",...) 逐字节一致
  peprec_mod <- data.frame(spec_id = "p000002", modifications = "4|Oxidation",
                           peptide = "PEPT", charge = "2",
                           orig_mod = "4,Oxidation[M]", stringsAsFactors = FALSE)
  parsed2 <- parsed$by_id["p000001"]; names(parsed2) <- "p000002"
  lines2 <- build_predict_txt(peprec_mod, parsed2)$lines
  say(length(lines2) == 7 && identical(lines2[2], "TITLE=PEPT|4,Oxidation[M]|2"),
      "修饰肽键与消费端逐字节一致", lines2[2])
  # 6) N-11：未修饰肽 peprec modifications 写 "-" 而非空串
  say(identical(pfind_mods_to_peprec(""), "-"), "未修饰肽写 -")
  # 7) 跳过计数：无产物 / 离子数不符均上报
  skip_built <- build_predict_txt(peprec_mod[0, ], parsed2)
  say(skip_built$skip_no_rows == 0L, "空表不误报跳过")
  bad <- data.frame(spec_id = "p000009", modifications = "-", peptide = "PEPT",
                    charge = "2", orig_mod = "", stringsAsFactors = FALSE)
  pbad <- list(b1 = 1, b2 = 1, y1 = 1, y2 = 1)
  names(pbad) <- NULL; pbad <- list(p000009 = pbad)
  say(build_predict_txt(bad, pbad)$skip_len == 1L, "离子数不符计数上报")
  unlink(td, recursive = TRUE)
  cat(if (ok) "== pred_ms2pip 自检全部通过 ==\n" else "== pred_ms2pip 自检存在失败项 ==\n")
  ok
}

# 仅当本文件作为脚本直接运行（--file= 指向自身）时执行 main；被 source 时只定义函数
cli_detect <- commandArgs(FALSE)
is_main_script <- length(cli_detect) > 0 &&
  any(grepl("pred_ms2pip\\.R$", cli_detect[grepl("^--file=", cli_detect)]))
if (is_main_script) {
  st <- tryCatch(main(), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); 1L })
  quit(status = as.integer(st))
}
