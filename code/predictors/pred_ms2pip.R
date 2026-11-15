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
#  #7：词表扩张——pFind 名与 Unimod 标题基本同名（psm_utils 直接认），
# 常见修饰直通补齐（现场曾点名 Carbamyl/Formyl/ICPL_13C(6) 等"未知名跳过"）；
# 引擎槽特有名（非 Unimod 标题形态）仍走跳过＋警告。
#  #7: vocabulary expansion - pFind names are mostly Unimod titles
# (accepted by psm_utils directly); common mods pass through. Engine-slot-only
# names (non-Unimod-looking) keep the skip-and-warn path.
mod_name_map <- c(Oxidation = "Oxidation", Deamidated = "Deamidated",
                  Carbamidomethyl = "Carbamidomethyl", Acetyl = "Acetyl",
                  Phospho = "Phospho", Methyl = "Methyl", Dimethyl = "Dimethyl",
                  Trimethyl = "Trimethyl", `Gln->pyro-Glu` = "Gln->pyro-Glu",
                  `Glu->pyro-Glu` = "Glu->pyro-Glu",
                  Carbamyl = "Carbamyl", Formyl = "Formyl",
                  Crotonyl = "Crotonyl", Succinyl = "Succinyl",
                  Malonyl = "Malonyl", Butyryl = "Butyryl",
                  Propionyl = "Propionyl", Propionamide = "Propionamide",
                  `Pyro-carbamidomethyl` = "Pyro-carbamidomethyl",
                  GG = "GG", Amidated = "Amidated", Biotin = "Biotin",
                  `Sodium` = "Sodium", `Potassium` = "Potassium",
                  `Cysteine->Selenocysteine` = "Cysteine->Selenocysteine")
#：同位点双修饰括号记法拆分。工作站 11-05 测试报告一例：
# "11,Carbamidomethyl[K](Gly[K])" —— pFind 对同一位点的第二个修饰以括号
# 内嵌；整条 token 在旧解析里既匹配不上映射表（尾部只剥到 [K]，剩下
# "Carbamidomethyl[K](Gly"），连本可映射的 Carbamidomethyl 也被连坐丢弃，
# 且该肽以无修饰态送去预测（预测口径失真）。拆成两个独立 token 后各自
# 走既有映射/未知名跳过路径：主修饰正常映射，次修饰（如 Gly）按未知名
# 跳过并警告——与分号分隔的多修饰同纪律。
#: split the same-site double-modification paren notation seen in
# the workstation 11-05 report ("11,Carbamidomethyl[K](Gly[K])"). The whole
# token previously failed to map (the trailing-bracket strip left
# "Carbamidomethyl[K](Gly"), dragging the mappable Carbamidomethyl down with
# it and sending the peptide to prediction as unmodified. After the split
# each token walks the existing map/skip path: the primary maps normally,
# the secondary (e.g. Gly) skips with a warning - same discipline as
# semicolon-separated multi-mods.
split_paren_mod <- function(mod_token) {
  m <- regmatches(mod_token,
                  regexec("^([0-9]+),([^()]*)\\((.*)\\)$", mod_token))[[1]]
  if (length(m) == 0L) return(mod_token)
  inner <- sub(";+$", "", m[4])
  # 内嵌修饰自带位置（"(2,Acetyl[K])"）则尊重之；无位置（"(Gly[K])"）＝
  # 同位点，继承外层位置
  # an inner mod with its own position ("(2,Acetyl[K])") keeps it; a bare
  # "(Gly[K])" is same-site and inherits the outer position
  inner_tok <- if (grepl("^[0-9]+,", inner)) inner else paste0(m[2], ",", inner)
  c(paste0(m[2], ",", m[3]), inner_tok)
}

#' pFind 修饰串（"3,Oxidation[M];0,Acetyl[AnyN-term];"）→ peprec modifications 串
#' （"3|Oxidation|0|Acetyl"， 起多修饰以 "|" 扁平连接）。未知名丢弃并
#' 警告（与 pDeep2 路径的 sanitize 同理）。
#' 无修饰返回 "-"（N-11：MS2PIP peprec 对未修饰肽要求写 "-"，空串会解析失败）。
pfind_mods_to_peprec <- function(mod, warn = TRUE) {
  if (is.null(mod) || !nzchar(mod)) return("-")
  mod_tokens <- strsplit(mod, ";", fixed = TRUE)[[1]]
  mod_tokens <- mod_tokens[nzchar(mod_tokens)]
  mod_tokens <- unlist(lapply(mod_tokens, split_paren_mod), use.names = FALSE)
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
  # 同位点双修饰防线。psm_utils 把同一
  # peprec 位置的两个修饰序列化成堆叠标签（如 N 端 "[Acetyl][Gln->pyro-Glu]-"），
  # pyteomics 的 ProForma 状态机解析不了（post_tag_before unexpected [），
  # ms2pip v4 读 PSM 阶段即整批崩——实测 peprec 667,889 行中 73,276 行
  # （11.0%）含同位点双修饰、12,709 行为 N 端双修饰，其中大量系 open-search
  # 把同一修饰写两遍的伪迹。两级处理：①完全重复的 pos|Name 条目折叠（纯
  # 伪迹，零语义变化）；②折叠后仍有同位点异名修饰→该肽按未知名同款纪律
  # 整条跳过（返回 NA，build_peprec 剔行＋计数）——牺牲 11% 换 89% 全量产出。
  #  (workstation 11-05 R4 root cause): same-position double-mod
  # guard. psm_utils serializes two mods at one peprec position as stacked
  # tags (e.g. "[Acetyl][Gln->pyro-Glu]-" at the N-terminus), which the
  # pyteomics ProForma state machine cannot parse (post_tag_before
  # unexpected [), and ms2pip v4 dies batch-wide while reading PSMs - of
  # 667,889 peprec rows, 73,276 (11.0%) carry same-position doubles, 12,709
  # of them N-terminal, many being open-search artifacts writing the same
  # mod twice. Two-tier handling: (1) exact duplicate pos|Name entries are
  # collapsed (pure artifact, zero semantic change); (2) peptides still
  # holding two different mods at one position after collapsing are skipped
  # whole like the unknown-name discipline (NA return; build_peprec drops
  # and counts the row) - trading 11% for the other 89% of the batch.
  if (length(mapped_mods)) {
    mapped_mods <- mapped_mods[!duplicated(mapped_mods)]
    pos_v <- sub("\\|.*$", "", mapped_mods)
    if (anyDuplicated(pos_v) > 0L) {
      if (warn) cat("同位点双修饰（整肽跳过）:", mod, "\n")
      return(NA_character_)
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
  # pep_df：peptide.txt 三列表（peptide/modification/charge）；返回 peprec 表。
  #：同位点双修饰肽返回 NA → 剔行＋spec_id 重编号（v4 长表的
  # psm_index 按首现顺序对齐 p000001..pK 连续编号，剔行后必须重编）＋
  # attr(same_site_dropped) 计数（warn=FALSE 关闭逐行告警防刷屏）。
  #: same-position double-mod peptides map to NA and are dropped
  # with spec_id renumbered (the v4 long table aligns psm_index to contiguous
  # p000001..pK by first appearance); attr(same_site_dropped) carries the
  # count (warn=FALSE silences the per-row flood).
  n <- nrow(pep_df)
  mods_v <- vapply(seq_len(n), function(i)
    pfind_mods_to_peprec(as.character(pep_df$modification[i]), warn = FALSE),
    character(1))
  # （工作站 run8 实锤，全表崩因之二）：肽序列含 "|"＝pFind 二选一
  # 异构体记法（如 QDPAAAQ…EEGEK|QDPAAAQ…GEGEK，实测 918 条、最多 5+ 变体
  # 连写）——pyteomics 不认，ms2pip v4 在 Reading PSMs 阶段整批崩。167 的
  # 修饰防线通过后（跳 8,687 条）仍死于此。取首变体会发明搜库未做的指认，
  # 按未知名同款纪律整肽跳过＋计数。
  #  (workstation run8, batch-killer #2): a "|" inside the peptide is
  # pFind's either-or isoform notation (e.g. QDPAAAQ...EEGEK|QDPAAAQ...GEGEK,
  # 918 rows observed, up to 5+ variants chained) - pyteomics cannot parse it
  # and ms2pip v4 dies batch-wide at Reading PSMs even after the 167 mod guard
  # (8,687 skipped) passed. Taking the first variant would invent an
  # assignment the search never made; skip the peptide whole like the
  # unknown-name discipline, with a count.
  pipe_hit <- grepl("|", as.character(pep_df$peptide), fixed = TRUE)
  # U/X 哨兵残基第四道
  # 卫兵。T@174 全表实测"长度不符 584,043（89%）"的根因链＝v4 把 39 条
  # U/X PSM 从输出丢弃但保留原 psm_index（留洞）→ 本适配器"排序后顺序
  # 命名 p%06d"在第一个洞（≈28,339 位）之后全部错位 → 长度门大量巧合性
  # 失配（幸存 74,202＝洞前前缀＋邻位长度巧合，Σp²≈11% 吻合）；"无产物
  # 39"＝末尾 39 个 spec_id 无组而非 U/X 计数。4,200 样本零 U/X→零洞→
  # 97.5% 过，同证。适配器侧先跳过 U/X → v4 永无洞、对齐恒成立。
  #  (workstation 174-round P1 verdict, chained locally): the
  # fourth guard for U/X sentinel residues. The full-table "length mismatch
  # 584,043 (89%)" root chain: v4 drops the 39 U/X PSMs from its output
  # while KEEPING their original psm_index values (gaps) -> this adapter's
  # sorted-index sequential naming misaligns everything after the first gap
  # (~position 28,339) -> the length gate fails coincidentally en masse
  # (survivors 74,202 = pre-gap prefix + adjacent-length coincidences,
  # matching the ~11% sum-of-squares collision rate); the "no-output 39"
  # was the trailing 39 unmatched spec_ids, not a U/X count. The 4,200-row
  # sample had zero U/X -> zero gaps -> 97.5% pass, confirming the chain.
  # Skipping U/X adapter-side means v4 never has gaps and the alignment
  # always holds.
  ux_hit <- grepl("[UX]", as.character(pep_df$peptide))
  keep <- !is.na(mods_v) & !pipe_hit & !ux_hit
  n_drop <- sum(!keep)
  out <- data.frame(
    spec_id = sprintf("p%06d", seq_len(sum(keep))),
    modifications = mods_v[keep],
    peptide = as.character(pep_df$peptide)[keep],
    charge = as.character(pep_df$charge)[keep],
    orig_mod = sub(";+$", "", as.character(pep_df$modification))[keep],
    stringsAsFactors = FALSE)
  attr(out, "same_site_dropped") <- sum(is.na(mods_v))
  attr(out, "pipe_dropped") <- sum(pipe_hit & !is.na(mods_v))
  attr(out, "ux_dropped") <- sum(ux_hit & !pipe_hit & !is.na(mods_v))
  out
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
  # （工作站 173 报告发现②，T@171/run11 实测）：原逐行 c() 累积
  # ＝平方级（65.8 万块逐块全量重扫），与 parse_ms2pip_long 的逐组全表扫描
  # 同为 166 在 APD 侧修过的病根在 ms2pip 侧重演。本版向量化：match 对齐＋
  # 预分配交错写入＋分片 10 万块心跳。语义不变（TITLE 键/跳过计数/y 位点序
  # 反转/format 口径逐字保留）。
  #  (workstation 173-round find 2, T@171/run11): the per-row c()
  # accumulation was quadratic (658k blocks), the same disease 166 fixed on
  # the APD side replaying here. Vectorized: match alignment + pre-allocated
  # interleaved writes + a heartbeat every 100k blocks. Semantics unchanged
  # (TITLE key, skip counts, y site-order reversal, format() wording).
  n <- nrow(peprec_table)
  m <- match(peprec_table$spec_id, names(parsed_by_id))
  have <- !is.na(m)
  skip_no_rows <- sum(!have)
  residue_count <- nchar(peprec_table$peptide)
  b1_len <- integer(n)
  hh <- which(have)
  if (length(hh)) {
    b1_len[hh] <- vapply(hh, function(i)
      length(parsed_by_id[[m[i]]]$b1), integer(1))
  }
  keep <- which(have & b1_len == residue_count - 1L)
  skip_len <- sum(have & b1_len != residue_count - 1L)
  nk <- length(keep)
  out <- character(7L * nk)
  if (nk) {
    i1 <- seq.int(1L, by = 7L, length.out = nk)
    out[i1] <- "BEGIN IONS"
    out[i1 + 1L] <- paste0("TITLE=", peprec_table$peptide[keep], "|",
                           peprec_table$orig_mod[keep], "|",
                           peprec_table$charge[keep])
    out[i1 + 6L] <- "END IONS"
    t0 <- proc.time()[3]
    idx_keep <- m[keep]
    CH <- 100000L
    n_ch <- ceiling(nk / CH)
    for (ci in seq_len(n_ch)) {
      s0 <- (ci - 1L) * CH + 1L
      s1 <- min(ci * CH, nk)
      quad <- vapply(idx_keep[s0:s1], function(j) {
        x <- parsed_by_id[[j]]
        c(paste0("b+1=", paste(format(x$b1, scientific = FALSE, trim = TRUE),
                               collapse = ",")),
          paste0("b+2=", paste(format(x$b2, scientific = FALSE, trim = TRUE),
                               collapse = ",")),
          paste0("y+1=", paste(rev(format(x$y1, scientific = FALSE, trim = TRUE)),
                               collapse = ",")),
          paste0("y+2=", paste(rev(format(x$y2, scientific = FALSE, trim = TRUE)),
                               collapse = ",")))
      }, character(4))
      out[i1[s0:s1] + 2L] <- quad[1L, ]
      out[i1[s0:s1] + 3L] <- quad[2L, ]
      out[i1[s0:s1] + 4L] <- quad[3L, ]
      out[i1[s0:s1] + 5L] <- quad[4L, ]
      cat(sprintf("[MS2PIP] 组装进度 %d/%d 块（已用 %.0fs）\n",
                  s1, nk, proc.time()[3] - t0))
      flush.console()
    }
  }
  list(lines = out, skip_no_rows = skip_no_rows, skip_len = skip_len)
}

#' MS2PIP v4 现代契约的长表解析。
#' v4.0+ 正式版 predict-batch -f tsv 输出：psm_index/ion_type/ion_number/
#' mz/predicted（log2 强度、单电荷 b/y）。b2/y2 以 0 填充（v4 批量输出无
#' 电荷维度）；psm_index 0 起按首现顺序对齐 peprec 的 p000001 起 spec_id。
parse_ms2pip_long <- function(tsv_file) {
  #：split-once O(n) 重写——原 raw[raw$psm_index==idx,] 逐组全表
  # 扫描＝组数×行数（65.8 万组×1,850 万行），843MB 长表实测 65min+ 仍在磨
  # ＝实质挂起（T@171 与 run11 的 2h+ 零 python 空转同源——v4 本体 2.5 分钟
  # 已完成、R 卡在此处）。列先抽取向量、data.frame 早期释放；每 10 万组
  # 心跳。
  #: split-once O(n) rewrite - the original per-group full-table
  # scan cost groups x rows (658k x 18.5M); the 843 MB long table ground for
  # 65+ minutes without finishing (an effective hang; T@171 and run11's 2h+
  # zero-python spins share this root - the v4 prediction itself finishes in
  # 2.5 minutes, the R side stalled here). Columns are pulled into vectors
  # with the data.frame released early; heartbeat every 100k groups.
  raw <- read.delim(tsv_file, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = character(0))
  missing_cols <- setdiff(c("psm_index", "ion_type", "ion_number", "predicted"),
                          colnames(raw))
  if (length(missing_cols)) stop("MS2PIP v4 长表缺列: ",
                                 paste(missing_cols, collapse = ","))
  psm_index <- suppressWarnings(as.integer(as.character(raw$psm_index)))
  ion_type <- as.character(raw$ion_type)
  ion_number <- suppressWarnings(as.integer(as.character(raw$ion_number)))
  predicted <- suppressWarnings(as.numeric(as.character(raw$predicted)))
  predicted[is.na(predicted)] <- 0
  intensity <- 2^predicted
  rm(raw, predicted)
  invisible(gc())
  grp <- split(seq_along(psm_index), psm_index)
  rm(psm_index)
  t0 <- proc.time()[3]
  by_id <- lapply(seq_along(grp), function(k) {
    ii <- grp[[k]]
    n_frag <- max(ion_number[ii])
    mk <- function(ion) {
      v <- numeric(n_frag)
      sel <- ion_type[ii] == ion
      v[ion_number[ii][sel]] <- intensity[ii][sel]
      v
    }
    if (k %% 100000L == 0L) {
      cat(sprintf("[MS2PIP] 长表解析进度 %d/%d 组（已用 %.0fs）\n",
                  k, length(grp), proc.time()[3] - t0))
      flush.console()
    }
    list(b1 = mk("b"), b2 = numeric(n_frag),
         y1 = mk("y"), y2 = numeric(n_frag))
  })
  stats::setNames(by_id, sprintf("p%06d", seq_along(by_id)))
}

# ----：产物指纹复用（与 pred_pdeep3/pred_alphapeptdeep 同款；      ----
# ---- 工作站 167 轮建议：四个适配器统一"产物在即复用"语义）。指纹＝输入       ----
# ---- md5＋size（master 的 file.copy 每轮刷新 peptides.txt mtime，mtime 序    ----
# ---- 恒失效），侧车 <out>.run_params。早期实现产物无侧车＝不可证新鲜→重跑。
# ----: artifact fingerprint reuse (same pattern as the pdeep3 and
# ---- alphapeptdeep adapters; the workstation's 167-round request to unify
# ---- "artifact present means reusable" across all four adapters). The
# ---- fingerprint is input md5 + size (the master's file.copy refreshes
# ---- peptides.txt mtime every run, so mtime ordering never holds) in an
# ---- <out>.run_params sidecar. Legacy artifacts without a sidecar rerun once.
input_fp <- function(in_path) {
  c(paste0("in_md5=", unname(tools::md5sum(in_path))),
    paste0("in_size=", file.info(in_path)$size))
}
sidecar_matches <- function(out_path, fp) {
  sc <- paste0(out_path, ".run_params")
  file.exists(sc) && identical(trimws(readLines(sc, warn = FALSE)), fp)
}
write_run_params <- function(out_path, in_path) {
  try(writeLines(input_fp(in_path), paste0(out_path, ".run_params"),
                 useBytes = TRUE), silent = TRUE)
}

# ---- 主入口 ----
main <- function() {
  cli_args <- parse_cli()
  if (isTRUE(cli_args$`self-test`)) return(if (self_test()) 0L else 1L)
  required_args <- c("in", "out", "python")
  missing_args <- setdiff(required_args, names(cli_args))
  if (length(missing_args)) stop("缺参数: ", paste(paste0("--", missing_args), collapse = " "))

  #：产物指纹复用——输入未变则整跑跳过（ms2pip 本体小时级 CPU）
  #: fingerprint reuse - an unchanged input skips the whole run
  # (the ms2pip prediction itself costs CPU-hours on full tables)
  if (file.exists(cli_args$`out`) && sidecar_matches(cli_args$`out`, input_fp(cli_args$`in`))) {
    cat("[MS2PIP] 复用已有产物（指纹匹配；删除 predict.txt 可强制重跑）：",
        cli_args$`out`, "\n", sep = "")
    return(0L)
  }
  if (file.exists(cli_args$`out`)) {
    cat("[MS2PIP] 肽表变化或旧产物无指纹——旧产物失效，重新预测\n")
    cat("[MS2PIP] stale predict invalidated (peptide table changed or no fingerprint); re-predicting\n")
  }
  #：--params 改为可选——网页的 MS2PIP 行按设计留空（默认取适配器
  # 内置 HCD 模型），此前把它列入必填，导致"网页填参数→跑"这条路径上 ms2pip
  # 恒失败（exit 1，缺 --params），该列整列丢失。空/缺省一律按 HCD 处理。
  #: --params is optional - the web form intentionally leaves the
  # MS2PIP row empty (default = the adapter's built-in HCD model). Requiring it
  # made ms2pip fail on every form-driven run (exit 1, missing --params) and
  # dropped the whole column. Missing/empty now means HCD.
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
  #：同位点双修饰整条跳过计数上报（ms2pip v4 对堆叠标签整批崩，
  # 跳过是换取其余肽全量产出的必要降级）
  #: report the same-position double-mod skips (ms2pip v4 dies
  # batch-wide on stacked tags; skipping is the necessary degrade)
  ss_dropped <- attr(peprec_table, "same_site_dropped")
  if (!is.null(ss_dropped) && ss_dropped > 0L) {
    cat("[MS2PIP] 同位点双修饰肽整条跳过 ", ss_dropped,
        " 条（ProForma 堆叠标签会使 ms2pip v4 整批崩溃，详见上游逐条告警）\n", sep = "")
    cat("[MS2PIP] skipped ", ss_dropped,
        " same-position double-mod peptides (stacked ProForma tags crash ms2pip v4 batch-wide)\n", sep = "")
  }
  #：管道符异构体肽计数上报（全表崩因之二）
  pipe_dropped <- attr(peprec_table, "pipe_dropped")
  if (!is.null(pipe_dropped) && pipe_dropped > 0L) {
    cat("[MS2PIP] 管道符异构体肽（pFind 二选一记法，pyteomics 不认）整条跳过 ",
        pipe_dropped, " 条\n", sep = "")
    cat("[MS2PIP] skipped ", pipe_dropped,
        " pipe-notation isoform peptides (pFind either-or notation, unparseable by pyteomics)\n", sep = "")
  }
  #：U/X 哨兵残基计数上报（对齐防线——见 build_peprec 注）
  #: U/X sentinel-residue count (the alignment guard - see build_peprec)
  ux_dropped <- attr(peprec_table, "ux_dropped")
  if (!is.null(ux_dropped) && ux_dropped > 0L) {
    cat("[MS2PIP] U/X 哨兵残基肽（v4 会静默丢弃并使 psm_index 留洞、错位全部对齐）整条跳过 ",
        ux_dropped, " 条\n", sep = "")
    cat("[MS2PIP] skipped ", ux_dropped,
        " U/X sentinel-residue peptides (v4 silently drops them and gaps psm_index, misaligning everything)\n", sep = "")
  }
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
    #：组数一致性告警——v4 若静默丢弃任何 PSM（如未卫兵的非法残基），
    # psm_index 留洞会使顺序命名错位、长度门大量巧合失配（174 轮 89% 实锤）。
    # 卫兵齐备时组数应恰等于 peprec 行数；不等＝仍有被丢的 PSM，大声告警。
    #: group-count consistency warning - if v4 silently drops any
    # PSM (e.g. an unguarded illegal residue), the gapped psm_index breaks
    # the sequential naming and the length gate fails coincidentally en
    # masse (the 174-round 89%). With the guards in place the group count
    # must equal the peprec rows; anything else means dropped PSMs remain.
    if (length(long_by_id) != nrow(peprec_table)) {
      cat("[MS2PIP][警告] v4 输出组数 ", length(long_by_id), " ≠ peprec 行数 ",
          nrow(peprec_table), "（差 ", nrow(peprec_table) - length(long_by_id),
          "）——仍有 PSM 被 v4 静默丢弃，psm_index 对齐已破坏，长度门将大量失配；",
          "请把被丢残基并入卫兵\n", sep = "")
      cat("[MS2PIP][WARN] v4 output groups ", length(long_by_id), " != peprec rows ",
          nrow(peprec_table), " (", nrow(peprec_table) - length(long_by_id),
          " dropped) - v4 still silently dropped PSMs, psm_index alignment is",
          " broken and the length gate will fail en masse; add the dropped",
          " residues to the guards\n", sep = "")
    }
    predict_result <- build_predict_txt(peprec_table, long_by_id)
    if (!length(predict_result$lines))
      stop("MS2PIP v4 结果与 peptide.txt 无匹配（肽长/离子数不符）")
    writeLines(predict_result$lines, cli_args$`out`, useBytes = TRUE)
    write_run_params(cli_args$`out`, cli_args$`in`)   #：指纹侧车
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
  write_run_params(cli_args$`out`, cli_args$`in`)   #：指纹侧车
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
  # 1) pFind 修饰 → peprec（期望串随  的 "|" 扁平连接更新——旧断言
  #    "3|Oxidation;0|Acetyl" 系 154c 改连接符前遗留，此前未随版本重跑）
  say(identical(pfind_mods_to_peprec("3,Oxidation[M];0,Acetyl[AnyN-term];"), "3|Oxidation|0|Acetyl"),
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
  # 4b)：v4 长表 split-once 解析（此前无自检覆盖）＋向量化组装
  td4 <- tempfile("ms2pip_long_"); dir.create(td4)
  lt <- file.path(td4, "pred.tsv")
  writeLines(c("psm_index\tion_type\tion_number\tmz\tpredicted",
               "0\tb\t1\t100.0\t1.0", "0\tb\t2\t200.0\t2.0",
               "0\ty\t1\t900.0\t3.0", "0\ty\t2\t800.0\t0.0",
               "1\tb\t1\t110.0\t0.5", "1\ty\t1\t910.0\t1.5"), lt)
  lg <- parse_ms2pip_long(lt)
  say(identical(names(lg), c("p000001", "p000002")) &&
        identical(lg$p000001$b1, c(2, 4)) &&          # 2^1, 2^2
        identical(lg$p000001$y1, c(8, 1)) &&          # ion1→2^3, ion2→2^0
        identical(lg$p000002$b1, 2^0.5) &&
        identical(lg$p000001$b2, c(0, 0)),
      "v4 长表 split-once 解析（b/y 装配＋2^predicted＋p%06d 命名＋b2 零填）",
      paste(sapply(lg, function(z) paste(z$b1, collapse = ",")), collapse = " / "))
  unlink(td4, recursive = TRUE)
  # 6b)：同位点双修饰括号记法——主修饰保住、次修饰按未知名跳过
  say(identical(pfind_mods_to_peprec("11,Carbamidomethyl[K](Gly[K])"),
                "11|Carbamidomethyl"),
      "同位点双修饰括号拆分（Gly 跳过、Carbamidomethyl 保留）",
      pfind_mods_to_peprec("11,Carbamidomethyl[K](Gly[K])"))
  say(identical(pfind_mods_to_peprec("3,Oxidation[M](2,Acetyl[K])"),
                "3|Oxidation|2|Acetyl"),
      "双修饰括号拆分后扁平连接",
      pfind_mods_to_peprec("3,Oxidation[M](2,Acetyl[K])"))
  # 6c)：同位点双修饰两级防线——完全重复折叠 / 真双修饰整肽 NA
  say(identical(pfind_mods_to_peprec("0,Gln->pyro-Glu[AnyN-term];0,Gln->pyro-Glu[AnyN-term];"),
                "0|Gln->pyro-Glu"),
      "同位点同一修饰写两遍＝伪迹折叠",
      pfind_mods_to_peprec("0,Gln->pyro-Glu[AnyN-term];0,Gln->pyro-Glu[AnyN-term];"))
  say(is.na(pfind_mods_to_peprec("0,Acetyl[AnyN-term];0,Gln->pyro-Glu[Q];")),
      "同位点异名双修饰→整肽 NA（ProForma 堆叠标签防线）")
  # 6d)：build_peprec 剔行＋spec_id 重编号＋计数
  df_ss <- data.frame(peptide = c("PEPT", "PEPT", "PEPT"),
                      modification = c("", "0,Acetyl[AnyN-term];0,Gln->pyro-Glu[Q];", ""),
                      charge = "2", stringsAsFactors = FALSE)
  pr_ss <- build_peprec(df_ss)
  say(nrow(pr_ss) == 2 &&
        identical(pr_ss$spec_id, c("p000001", "p000002")) &&
        identical(attr(pr_ss, "same_site_dropped"), 1L),
      "同位点双修饰剔行＋spec_id 重编号＋attr 计数",
      paste0(nrow(pr_ss), " 行 / dropped=", attr(pr_ss, "same_site_dropped")))
  # 6e)：管道符异构体肽（pFind 二选一记法）整肽跳过＋计数
  df_pipe <- data.frame(peptide = c("PEPT", "QDPAAAQEKG|QDPAAAQEGK", "ALAMAR"),
                        modification = c("", "", ""),
                        charge = "2", stringsAsFactors = FALSE)
  pr_pipe <- build_peprec(df_pipe)
  say(nrow(pr_pipe) == 2 &&
        !any(grepl("|", pr_pipe$peptide, fixed = TRUE)) &&
        identical(attr(pr_pipe, "pipe_dropped"), 1L) &&
        identical(attr(pr_pipe, "same_site_dropped"), 0L),
      "管道符异构体肽整肽跳过＋计数＋重编号",
      paste0(nrow(pr_pipe), " 行 / pipe=", attr(pr_pipe, "pipe_dropped")))
  # 6f)：U/X 哨兵残基第四道卫兵（对齐防线——v4 静默丢弃留洞）
  df_ux <- data.frame(peptide = c("PEPT", "PEPTIDXK", "ALAMAR"),
                      modification = c("", "", ""),
                      charge = "2", stringsAsFactors = FALSE)
  pr_ux <- build_peprec(df_ux)
  say(nrow(pr_ux) == 2 &&
        !any(grepl("[UX]", pr_ux$peptide)) &&
        identical(attr(pr_ux, "ux_dropped"), 1L) &&
        identical(attr(pr_ux, "pipe_dropped"), 0L),
      "U/X 哨兵残基肽整肽跳过＋计数＋重编号",
      paste0(nrow(pr_ux), " 行 / ux=", attr(pr_ux, "ux_dropped")))
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
