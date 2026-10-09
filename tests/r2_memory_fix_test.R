# r2_memory_fix_test.R — R2 工作站内存爆炸修复的冒烟与等价性测试（2026-09-24）
# 覆盖：
#   1. slot_state_matrix 向量化版 vs 逐行朴素参照（随机矩阵，含空串/NA/重复）
#   2. pfind_res_postprocess 新版（reader 过滤＋big-string 回帖缓存＋批量写回）
#      vs 旧算法逐行参照——小夹具全表 identical 比对
#   3. pFind_res_reader keep_spectra/drop_decoy_candidates 过滤语义等价
#   4. merge_search_results 新版（无逐引擎 cbind、状态向量化、peptide 列表收集）
#      vs 旧算法参照——小夹具输出 identical 比对
#   5. 放大应激：15 万行合并 / 6 万行×2 万诱饵后处理——线性完成、堆内存 < 2GB、
#      后处理另做 300 行子集逐行参照抽查
# 用法：在 tests\ 目录下 Rscript r2_memory_fix_test.R（--file= 相对定位）
options(stringsAsFactors = FALSE)
args0 <- commandArgs(FALSE)
sd <- dirname(normalizePath(sub("^--file=", "", args0[grepl("--file=", args0)][1]),
                            winslash = "/"))
cd <- file.path(dirname(sd), "code")
if (!dir.exists(cd)) cd <- file.path(sd, "..", "code")
source(file.path(cd, "compat.R"))
source(file.path(cd, "ReadMGF_script.R"))

npass <- 0L; nfail <- 0L
ck <- function(label, ok) {
  cat(ifelse(ok, "PASS", "FAIL"), ":", label, "\n")
  if (ok) npass <<- npass + 1L else nfail <<- nfail + 1L
}

# ============================================================================
# 1. slot_state_matrix 向量化版 vs 逐行朴素参照
# ============================================================================
naive_state <- function(m, sym) {
  # 旧 merge 里的逐行 sapply 语义（NA 视同空值：生产读入用 na.strings=character(0)
  # 不会产生 NA，向量化版以空值口径对齐，参照同样先 NA→""）
  apply(m, 1L, function(s) {
    s[is.na(s)] <- ""
    u <- unique(s)
    u <- u[u != ""]
    if (length(u) <= 1L) return("Same")
    paste(vapply(u, function(z) paste(sym[which(s == z)], collapse = ""), ""),
          collapse = ";")
  })
}
set.seed(42)
for (k in 1:6) {
  for (rep in 1:3) {
    n <- 80
    pool <- c("", "A", "B", "C", "D")
    m <- matrix(sample(pool, n * k, replace = TRUE, prob = c(.25, rep(.1875, 4))),
                n, k)
    m[sample(seq_len(n * k), 6)] <- NA
    got <- slot_state_matrix(m, letters[seq_len(k)])
    exp <- naive_state(m, letters[seq_len(k)])
    ck(sprintf("slot_state 等价 k=%d rep=%d（含空串/NA/重复）", k, rep),
       identical(got, exp))
  }
}
ck("slot_state 全空行 = Same",
   identical(slot_state_matrix(matrix("", 4, 3), letters[1:3]),
             rep("Same", 4)))
ck("slot_state 单列 k=1",
   identical(slot_state_matrix(matrix(c("A", "B", "A"), 3, 1), "a"),
             c("Same", "Same", "Same")))
ck("slot_state 0 行返回空向量",
   identical(slot_state_matrix(matrix(character(0), 0, 3), letters[1:3]),
             character(0)))

# ============================================================================
# 2. pfind_res_postprocess 小夹具：新版 vs 旧算法逐行参照（全表 identical）
# ============================================================================
td <- tempfile("r2fix_"); dir.create(td)
# --- fasta：8 蛋白（含共享肽、重复蛋白内命中） ---
fa_lines <- c(
  ">p1", "MKPEPTIDEKAAAFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
  ">p2", "MRANOTHERRLLLLPEPTIDEKLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL",
  ">p3", "MFLAGK", ">p4", "MVVVVVVV", ">p5", "MPEPTIDEKPEPTIDEKGGGG",
  ">p6", "MZZZOTHERRKK", ">p7", "MCCCCCCCC", ">p8", "MDDDDDDDD")
fa_file <- file.path(td, "test.fasta"); writeLines(fa_lines, fa_file, useBytes = TRUE)

# --- 诱饵候选的 qry.res（12 谱，6 诱饵谱各含诱饵+非诱饵候选） ---
titles <- sprintf("sp%02d.dta", 1:12)
decoy_seq <- c("PEPTIDEK", "ANOTHERR", "FLAGK", "ZZZOTHERR", "CCCCCCCC", "DDDDDDDD")
mk_block <- function(title, seqs, raws, mods_n = NULL, decoy_flags = NULL,
                     invalid_flags = NULL) {
  # 每个候选：rank seq calc raw eval fdr nmods [mod 对] 4aux 蛋白列 打分明细列
  out <- c("S\t1234.567\t2", title)
  np <- length(seqs)
  for (j in seq_len(np)) {
    dec <- if (!is.null(decoy_flags)) decoy_flags[j] else FALSE
    inv <- if (!is.null(invalid_flags)) invalid_flags[j] else FALSE
    nmods <- if (!is.null(mods_n)) mods_n[j] else 0L
    body <- sprintf("%d\t%s\t%.2f\t%.2f\t%.2g\t%.3f\t%d",
                    j, seqs[j], 1000 + j * 0.1, raws[j], 1e-9 - j * 1e-10,
                    0.001 + j * 0.001, nmods)
    if (nmods > 0L) body <- paste0(body, "\t3\tAnyN-term#0") # 内部位点 #0 末端名 → invalid
    body <- paste0(body, "\t", if (inv) "1" else "0", "\t0\t0\t0")
    # 蛋白列：诱饵候选挂 1;1 旗标，非诱饵挂 0;0（蛋白序号按 j 取 1..8）
    body <- paste0(body, sprintf("\t%d;%d;10;K;P", (j %% 8L) + 1L, if (dec) 1L else 0L))
    body <- paste0(body, "\t88.5\t0.99")
    out <- c(out, body)
  }
  out
}
writeLines(c(
  mk_block(titles[1], c(decoy_seq[1], "PEPTIDEK"), c(60, 90),
           decoy_flags = c(TRUE, FALSE)),
  mk_block(titles[2], c(decoy_seq[2], "ANOTHERR"), c(70, 95),
           decoy_flags = c(TRUE, FALSE)),
  mk_block(titles[3], c(decoy_seq[3], "FLAGK"), c(50, 85),
           decoy_flags = c(TRUE, FALSE)),
  mk_block(titles[4], c(decoy_seq[4], "ZZZOTHERR"), c(40, 80),
           decoy_flags = c(TRUE, FALSE)),
  mk_block(titles[5], c(decoy_seq[5], "CCCCCCCC"), c(30, 75),
           decoy_flags = c(TRUE, FALSE)),
  # sp06：诱饵谱但只有诱饵候选 → n_keep；sp07：非诱饵谱不入 keep 集
  mk_block(titles[6], c("QEFEFEFEF"), c(99), decoy_flags = c(TRUE)),
  mk_block(titles[7], c("WYWYWYWY"), c(88)),
  # sp08：含 invalid delta 的候选（应跳过，取下一个有效候选）
  mk_block(titles[8], c(decoy_seq[6], "DDDDDDDD", "EEEEEEEE"),
           c(20, 77, 95), mods_n = c(0, 0, 1), decoy_flags = c(TRUE, FALSE, FALSE),
           invalid_flags = c(FALSE, FALSE, TRUE)),
  # sp09-12：非诱饵谱（不参与救援，仅占位）
  mk_block(titles[9], c("GGGGGGGG"), c(50)),
  mk_block(titles[10], c("HHHHHHHH"), c(50)),
  mk_block(titles[11], c("IIIIIIII"), c(50)),
  mk_block(titles[12], c("JJJJJJJJ"), c(50))
), file.path(td, "test.qry.res"), useBytes = TRUE)

# --- pFind.spectra：12 行，sp01-08 诱饵、sp09-12 目标 ---
sp0 <- data.frame(
  File_Name = titles,
  Sequence = c(decoy_seq, "WYWYWYWY", decoy_seq[6], "GGGGGGGG", "HHHHHHHH",
               "IIIIIIII", "JJJJJJJJ"),
  Modification = rep("", 12),
  `Calc.MH+` = rep(999.0, 12),
  `Exp.MH+` = rep(1000.0, 12),
  `Mass_Shift(Exp.-Calc.)` = rep(1.0, 12),
  Raw_Score = rep(10, 12),
  Final_Score = rep(0.5, 12),
  `Q.value` = rep(0.01, 12),
  Specificity = rep(0L, 12),
  Proteins = rep("1;0;10;K;P/", 12),
  Positions = rep("10,K,P/", 12),
  `Target/Decoy` = c(rep("decoy", 8), rep("target", 4)),
  `Miss.Clv.Sites` = rep(0L, 12),
  check.names = FALSE, stringsAsFactors = FALSE)
sp_file <- file.path(td, "pFind.spectra")
write.table(sp0, sp_file, sep = "\t", row.names = FALSE, quote = FALSE)

# --- 旧算法逐行参照（旧版函数体原样移植，作等价性基准） ---
ref_postprocess <- function(sp, res, fa_names, fa_seqs, rows = NULL) {
  td_col <- grep("^Target/Decoy$", colnames(sp), value = TRUE)
  if (length(td_col) == 0) td_col <- "Target.Decoy"
  is_dec <- sp[[td_col[1]]] == "decoy"
  idx <- if (is.null(rows)) which(is_dec) else intersect(which(is_dec), rows)
  n_rescue <- 0L; n_keep <- 0L
  for (i in idx) {
    fn <- as.character(sp$File_Name[i])
    blk <- res[[fn]]
    if (is.null(blk)) { n_keep <- n_keep + 1L; next }
    best <- NULL
    for (cd in blk$candidates) {
      if (cd$is_decoy) next
      if (isTRUE(cd$invalid_delta_site)) next
      if (is.null(best) || (cd$raw_score > best$raw_score)) best <- cd
    }
    if (is.null(best)) { n_keep <- n_keep + 1L; next }
    prots <- ""; poss <- ""
    if (length(fa_seqs) > 0) {
      hit_i <- which(vapply(fa_seqs, function(s) grepl(best$sequence, s, fixed = TRUE), logical(1)))
      if (length(hit_i) > 5) hit_i <- hit_i[1:5]
      pn <- character(0); pp <- character(0)
      for (hi in hit_i) {
        p0 <- regexpr(best$sequence, fa_seqs[hi], fixed = TRUE)
        st <- as.integer(p0)
        prev_aa <- if (st > 1) substr(fa_seqs[hi], st - 1, st - 1) else "-"
        next_aa <- if (st + nchar(best$sequence) <= nchar(fa_seqs[hi]))
          substr(fa_seqs[hi], st + nchar(best$sequence), st + nchar(best$sequence)) else "-"
        pn <- c(pn, fa_names[hi]); pp <- c(pp, sprintf("%d,%s,%s", st, prev_aa, next_aa))
      }
      prots <- paste0(paste(pn, collapse = "/"), "/")
      poss <- paste0(paste(pp, collapse = "/"), "/")
    }
    exp_mh <- suppressWarnings(as.numeric(sp[["Exp.MH+"]][i]))
    sp$Sequence[i] <- best$sequence
    sp[["Calc.MH+"]][i] <- best$calc_mh
    if ("Mass_Shift(Exp.-Calc.)" %in% colnames(sp) && is.finite(exp_mh))
      sp[["Mass_Shift(Exp.-Calc.)"]][i] <- exp_mh - best$calc_mh
    sp$Raw_Score[i] <- best$raw_score
    if ("Final_Score" %in% colnames(sp)) sp$Final_Score[i] <- best$e_value
    q_col <- grep("^Q[-.]value$", colnames(sp), value = TRUE)[1]
    if (!is.na(q_col))
      sp[[q_col]][i] <- if (is.finite(best$fdr_local) && best$fdr_local > 0)
        best$fdr_local else best$e_value
    sp$Modification[i] <- best$modification
    sp$Specificity[i] <- as.character(best$specificity)
    if (!is.na(prots) && nzchar(prots)) { sp$Proteins[i] <- prots; sp$Positions[i] <- poss }
    sp[[td_col[1]]][i] <- "target"
    mc_col <- grep("^Miss[.]Clv[.]Sites$", colnames(sp), value = TRUE)[1]
    if (!is.na(mc_col)) {
      kr <- gregexpr("[KR]", best$sequence, fixed = FALSE)[[1]]
      sp[[mc_col]][i] <- as.character(sum(kr > 0 & kr != nchar(best$sequence)))
    }
    n_rescue <- n_rescue + 1L
  }
  sp
}

# 参照读入（未过滤 reader + 旧循环）
res_full <- pFind_res_reader(file.path(td, "test.qry.res"))
sp_ref <- read.delim(sp_file, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, quote = "", comment.char = "",
                     na.strings = character(0))
fa_ln <- enc2utf8(readLines(fa_file, warn = FALSE))
hid <- grepl("^>", fa_ln)
fa_names <- vapply(strsplit(sub("^>", "", fa_ln[hid]), " "), function(x) x[1], "")
grp <- cumsum(hid)
fa_seqs <- as.character(tapply(fa_ln[!hid], grp[!hid], function(x) paste(x, collapse = "")))
sp_ref2 <- ref_postprocess(sp_ref, res_full, fa_names, fa_seqs)

# 新版：直接跑 pfind_res_postprocess 再读回输出
r <- pfind_res_postprocess(td, fa_file)
sp_new <- read.delim(file.path(td, "pFind-Filtered.spectra"), sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "", comment.char = "", na.strings = character(0))
# 参照结果同样经 write/read 往返（列类型推断一致后再 identical）
ref_file <- file.path(td, "ref_out.spectra")
write.table(sp_ref2, ref_file, sep = "\t", row.names = FALSE, quote = FALSE,
            fileEncoding = "UTF-8", na = "")
sp_ref3 <- read.delim(ref_file, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, quote = "", comment.char = "",
                      na.strings = character(0))
ck("postprocess 小夹具 输出与旧算法逐行参照 identical", identical(sp_new, sp_ref3))
ck("postprocess 小夹具 救援计数一致", r[["rescued"]] == 7L &&
     r[["kept"]] == 1L && r[["decoy"]] == 8L)

# --- 3. reader 过滤语义等价 ---
res_filt <- pFind_res_reader(file.path(td, "test.qry.res"),
                             keep_spectra = titles[1:8],
                             drop_decoy_candidates = TRUE)
ck("reader keep_spectra 只保留集合内谱",
   identical(sort(names(res_filt)), sort(titles[1:8])))
ok_filt <- TRUE
for (t0 in titles[1:8]) {
  full_cand <- res_full[[t0]]$candidates
  keep_cand <- full_cand[!vapply(full_cand, function(x) x$is_decoy, logical(1))]
  got <- res_filt[[t0]]$candidates
  if (length(got) != length(keep_cand)) { ok_filt <- FALSE; break }
  for (j in seq_along(got))
    if (!identical(got[[j]]$sequence, keep_cand[[j]]$sequence) ||
        !identical(got[[j]]$raw_score, keep_cand[[j]]$raw_score)) {
      ok_filt <- FALSE; break
    }
}
ck("reader drop_decoy_candidates 过滤后与非诱饵子集一致", ok_filt)
rm(res_full, res_filt, sp_ref, sp_ref2, sp_ref3, sp_new)

# ============================================================================
# 4. merge_search_results 小夹具：新版 vs 旧算法参照（输出 identical）
# ============================================================================
ref_merge <- function(psm_file, search_files, out_dir, slot_letters, denovo_slots,
                      fasta_file, mode = "keep") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  hdr0 <- strsplit(readLines(psm_file, n = 1), "\t")[[1]]
  cc0 <- rep(NA_character_, length(hdr0))
  for (col0 in c("File_Name", "Sequence", "Modification")) {
    i0 <- which(hdr0 == col0)
    if (length(i0) == 1) cc0[i0] <- "character"
  }
  datalist <- read.delim(psm_file, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = cc0)
  for (mc in c("a.Modification", "Modification")) {
    if (mc %in% colnames(datalist)) {
      datalist[[mc]][is.na(datalist[[mc]])] <- ""
      datalist[[mc]][datalist[[mc]] == "NA"] <- ""
    }
  }
  spectra_list <- datalist$File_Name
  seqmod <- which(colnames(datalist) %in% c("Sequence", "Modification"))
  datalist <- cbind(datalist[, -seqmod], datalist[, c("Sequence", "Modification")])
  colnames(datalist)[which(colnames(datalist) %in% c("Sequence", "Modification"))] <-
    c("a.Sequence", "a.Modification")
  for (index in seq_along(search_files)) {
    sf <- search_files[index]; sl <- slot_letters[index]
    hdr <- strsplit(readLines(sf, n = 1), "\t")[[1]]
    cc <- rep("NULL", length(hdr))
    for (col in c("File_Name", "Sequence", "Modification")) {
      i <- which(hdr == col)
      if (length(i) == 1) cc[i] <- "character"
    }
    temp <- read.delim(sf, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                       colClasses = cc, na.strings = character(0))
    if (nrow(temp) == 0) {
      temp <- data.frame(File_Name = spectra_list,
                         Sequence = rep("", length(spectra_list)),
                         Modification = rep("", length(spectra_list)),
                         stringsAsFactors = FALSE)
    } else {
      if (!"Modification" %in% names(temp)) temp$Modification <- ""
      if (!"Sequence" %in% names(temp)) temp$Sequence <- ""
      temp <- temp[temp$File_Name %in% spectra_list, ]
      unmatched <- setdiff(spectra_list, temp$File_Name)
      if (length(unmatched) > 0) {
        u <- setNames(data.frame(File_Name = unmatched, X = "", Y = "",
                                 stringsAsFactors = FALSE), colnames(temp))
        temp <- rbind(temp, u)
      }
    }
    colnames(temp) <- c("File_Name", paste0(sl, ".Sequence"),
                        paste0(sl, ".Modification"))
    ord <- order(match(temp$File_Name, datalist$File_Name))
    datalist <- cbind(datalist, temp[ord, 2:3])
  }
  if (!is.null(fasta_file)) {
    for (ds in denovo_slots)
      datalist <- denovo_db_lookback(datalist, fasta_file, ds, mode = mode)
  }
  datalist <- cbind(datalist, Peptide_state = "", Modification_state = "")
  col_names <- colnames(datalist)
  seq_idx <- which(grepl("\\.Sequence$", col_names))
  seq_sym <- sapply(col_names[seq_idx], function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
  mod_idx <- which(grepl("\\.Modification$", col_names))
  datalist$Peptide_state <- sapply(seq_len(nrow(datalist)), function(x) {
    s <- as.character(datalist[x, seq_idx])
    s <- gsub("N", "GG", gsub("I", "L", s))
    u <- unique(s); u <- u[u != ""]
    if (length(u) <= 1) return("Same")
    paste(sapply(u, function(z) paste(seq_sym[which(s == z)], collapse = "")), collapse = ";")
  })
  datalist$Modification_state <- sapply(seq_len(nrow(datalist)), function(x) {
    s <- as.character(datalist[x, mod_idx])
    s <- gsub("I", "L", s)
    s <- gsub("ProteinN-term", "AnyN-term", s)
    u <- unique(s); u <- u[u != ""]
    if (length(u) <= 1) return("Same")
    paste(sapply(u, function(z) paste(seq_sym[which(s == z)], collapse = "")), collapse = ";")
  })
  write.table(datalist, file.path(out_dir, "1.combine_result.spectra"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  peptide_df <- c()
  for (pc in seq_idx) {
    tmp <- datalist[, c(pc, pc + 1, which(colnames(datalist) == "Charge"))]
    colnames(tmp) <- c("peptide", "modification", "charge")
    for (r in seq_len(nrow(tmp))) {
      vars <- apply_sap_variants(as.character(tmp$peptide[r]),
                                 as.character(tmp$modification[r]))
      if (length(vars) == 1) {
        tmp$peptide[r] <- vars[[1]]$seq
        tmp$modification[r] <- vars[[1]]$mod
      } else {
        extra <- tmp[r, ]
        tmp$peptide[r] <- vars[[1]]$seq; tmp$modification[r] <- vars[[1]]$mod
        extra$peptide <- vars[[2]]$seq; extra$modification <- vars[[2]]$mod
        tmp <- rbind(tmp, extra)
      }
      tmp$modification[r] <- sanitize_pdeep2_mods(tmp$modification[r],
                                                  nchar(tmp$peptide[r]))
    }
    peptide_df <- rbind(peptide_df, unique(tmp))
  }
  peptide_df <- unique(peptide_df)
  peptide_df <- peptide_df[peptide_df$peptide != "", ]
  write.table(peptide_df, file.path(out_dir, "peptide.txt"), sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  list(merged_rows = nrow(datalist), peptides = nrow(peptide_df))
}

set.seed(7)
mt <- tempfile("r2merge_"); dir.create(mt)
nrow_s <- 60L
titles_s <- sprintf("m%03d.dta", seq_len(nrow_s))
seqpool <- c("PEPTIDEK", "ANOTHERR", "MFLAGK", "MVVVVVV", "GGGNNGGG", "XLEPEPTIDE")
psm_s <- data.frame(File_Name = titles_s, Charge = rep(2L, nrow_s),
                    Sequence = sample(seqpool, nrow_s, replace = TRUE),
                    Modification = sample(c("", "1,Oxidation[M];", "2,Xle->I[X];"),
                                          nrow_s, replace = TRUE),
                    stringsAsFactors = FALSE, check.names = FALSE)
psm_s$Modification[psm_s$Sequence != "XLEPEPTIDE" & grepl("Xle", psm_s$Modification)] <- ""
write.table(psm_s, file.path(mt, "psm.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
mk_s <- function(tag, keep = TRUE, header_only = FALSE) {
  if (header_only) {
    writeLines("File_Name\tSequence\tModification",
               file.path(mt, paste0(tag, ".spectra")), useBytes = TRUE)
    return(invisible())
  }
  d <- data.frame(File_Name = titles_s,
                  Sequence = sample(seqpool, nrow_s, replace = TRUE),
                  Modification = sample(c("", "3,Oxidation[M];", "64,Oxidation[M];"),
                                        nrow_s, replace = TRUE),
                  stringsAsFactors = FALSE, check.names = FALSE)
  if (!keep) d <- d[sample(nrow_s, nrow_s - 7L), , drop = FALSE]   # 缺 7 谱 → 空槽补齐
  write.table(d, file.path(mt, paste0(tag, ".spectra")), sep = "\t",
              row.names = FALSE, quote = FALSE)
}
mk_s("open"); mk_s("restrict", keep = FALSE); mk_s("rescore")
mk_s("denovo", header_only = TRUE); mk_s("maxquant", keep = FALSE)
slots_s <- c("b", "c", "d", "e", "f")
# 修正 denovo 槽：用非空 denovo 文件（回帖场景）
d2 <- data.frame(File_Name = titles_s[1:40],
                 Sequence = c("MAGGLPEPTIDEK", sample(seqpool, 39, replace = TRUE)),
                 Modification = rep("", 40),
                 stringsAsFactors = FALSE, check.names = FALSE)
write.table(d2, file.path(mt, "denovo.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)

args_s <- list(file.path(mt, "psm.spectra"),
               file.path(mt, paste0(c("open", "restrict", "rescore", "denovo",
                                      "maxquant"), ".spectra")),
               file.path(mt, "merge_new"), slots_s, "e", fa_file, "keep")
res_new <- do.call(merge_search_results, c(args_s[1:3],
                                           list(fasta_file = args_s[[6]],
                                                slot_letters = slots_s,
                                                denovo_slots = "e",
                                                denovo_lookback_mode = "keep")))
res_ref <- ref_merge(args_s[[1]], args_s[[2]], file.path(mt, "merge_ref"),
                     slots_s, "e", fa_file, "keep")
mn <- read.delim(file.path(mt, "merge_new", "1.combine_result.spectra"), sep = "\t",
                 check.names = FALSE, na.strings = character(0))
mr <- read.delim(file.path(mt, "merge_ref", "1.combine_result.spectra"), sep = "\t",
                 check.names = FALSE, na.strings = character(0))
ck("merge 小夹具 合并表与旧算法参照 identical", identical(mn, mr))
pn <- read.delim(file.path(mt, "merge_new", "peptide.txt"), sep = "\t",
                 check.names = FALSE, na.strings = character(0))
pr <- read.delim(file.path(mt, "merge_ref", "peptide.txt"), sep = "\t",
                 check.names = FALSE, na.strings = character(0))
ck("merge 小夹具 peptide.txt 与旧算法参照 identical", identical(pn, pr))
ck("merge 小夹具 计数一致", identical(res_new, res_ref))
rm(mn, mr, pn, pr)

# ============================================================================
# 5a. merge 放大应激：15 万行 × 5 槽，线性完成 + 堆内存 < 2GB
# ============================================================================
set.seed(11)
mt2 <- tempfile("r2mergeb_"); dir.create(mt2)
nb <- 150000L
tb <- sprintf("b%06d.dta", seq_len(nb))
psm_b <- data.frame(File_Name = tb, Charge = rep(2L, nb),
                    Sequence = sample(seqpool, nb, replace = TRUE),
                    Modification = sample(c("", "1,Oxidation[M];"), nb, replace = TRUE),
                    stringsAsFactors = FALSE, check.names = FALSE)
write.table(psm_b, file.path(mt2, "psm.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
for (tg in c("e1", "e2", "e3", "e4")) {
  d <- data.frame(File_Name = tb,
                  Sequence = sample(seqpool, nb, replace = TRUE),
                  Modification = rep("", nb),
                  stringsAsFactors = FALSE, check.names = FALSE)
  d$Sequence[sample(nb, 3000L)] <- ""
  write.table(d, file.path(mt2, paste0(tg, ".spectra")), sep = "\t",
              row.names = FALSE, quote = FALSE)
}
t_b <- system.time(res_b <- merge_search_results(
  file.path(mt2, "psm.spectra"),
  file.path(mt2, paste0(c("e1", "e2", "e3", "e4"), ".spectra")),
  file.path(mt2, "merge"),
  slot_letters = c("b", "c", "d", "e"), denovo_slots = character(0)))
cat(sprintf("[stress] merge 150K 行：%.0fs，堆 %.0f MB\n",
            t_b[3], iso_mem_mb()))
ck("merge 应激 15 万行完成且行数正确", res_b$merged_rows == nb)
ck("merge 应激 15 万行 < 8 分钟", t_b[3] < 480)
ck("merge 应激 堆内存 < 2GB", iso_mem_mb() < 2000)
rm(psm_b); invisible(gc())

# ============================================================================
# 5b. postprocess 放大应激：6 万行/2 万诱饵 + 1200 蛋白，300 行子集参照抽查
# ============================================================================
set.seed(23)
mt3 <- tempfile("r2post_"); dir.create(mt3)
np_tot <- 60000L; np_dec <- 20000L
tt3 <- sprintf("q%06d.dta", seq_len(np_tot))
# fasta：1200 蛋白 × ~400 aa（候选序列以蛋白序列片段保证回帖命中）
nprot <- 1200L; protlen <- 400L
aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
prot_seqs <- vapply(seq_len(nprot), function(k) {
  paste(sample(aa, protlen, replace = TRUE), collapse = "")
}, "")
fa3 <- paste0(">prot", seq_len(nprot), "\n", prot_seqs)
fa3_file <- file.path(mt3, "big.fasta")
writeLines(fa3, fa3_file, useBytes = TRUE)
# 诱饵谱的救援候选序列：取各蛋白第 10..17 位 8 肽（必然命中）
rescue_seqs <- substr(prot_seqs[((seq_len(np_dec) - 1L) %% nprot) + 1L], 10, 17)
# qry.res：只写诱饵谱块（每块 4 诱饵候选 + 4 非诱饵候选）
blocks <- character(0)
for (j in seq_len(np_dec)) {
  blk <- c("S\t1234.567\t2", tt3[j])
  for (cand in 1:8) {
    dec <- cand <= 4
    sq <- if (dec) paste0("D", sprintf("%07d", j), sample(aa, 7, replace = TRUE), collapse = "")
          else rescue_seqs[j]
    blk <- c(blk, sprintf("%d\t%s\t%.2f\t%.2f\t%.2g\t%.3f\t0\t0\t0\t0\t0\t%d;%d;10;K;P\t88.5\t0.99",
                          cand, sq, 1000 + cand, if (dec) 50 + cand else 90 - cand,
                          1e-9 - cand * 1e-10, 0.001 + cand * 0.001,
                          (j %% nprot) + 1L, if (dec) 1L else 0L))
  }
  blocks <- c(blocks, blk)
}
qry3 <- file.path(mt3, "big.qry.res")
writeLines(blocks, qry3, useBytes = TRUE)
# pFind.spectra：np_tot 行（前 np_dec 诱饵，其余 target）
sp3 <- data.frame(File_Name = tt3,
                  Sequence = c(rescue_seqs, rep("TTTTTTTT", np_tot - np_dec)),
                  Modification = rep("", np_tot),
                  `Calc.MH+` = rep(999.0, np_tot),
                  `Exp.MH+` = rep(1000.0, np_tot),
                  `Mass_Shift(Exp.-Calc.)` = rep(1.0, np_tot),
                  Raw_Score = rep(10, np_tot), Final_Score = rep(0.5, np_tot),
                  `Q.value` = rep(0.01, np_tot), Specificity = rep(0L, np_tot),
                  Proteins = rep("1;0;10;K;P/", np_tot),
                  Positions = rep("10,K,P/", np_tot),
                  `Target/Decoy` = c(rep("decoy", np_dec), rep("target", np_tot - np_dec)),
                  `Miss.Clv.Sites` = rep(0L, np_tot),
                  check.names = FALSE, stringsAsFactors = FALSE)
write.table(sp3, file.path(mt3, "pFind.spectra"), sep = "\t",
            row.names = FALSE, quote = FALSE)
t_p <- system.time(r_p <- pfind_res_postprocess(mt3, fa3_file))
cat(sprintf("[stress] postprocess 60K 行/20K 诱饵：%.0fs，堆 %.0f MB，救援 %s\n",
            t_p[3], iso_mem_mb(), r_p[["rescued"]]))
ck("postprocess 应激 20K 诱饵全救援", r_p[["rescued"]] == np_dec)
ck("postprocess 应激 < 5 分钟", t_p[3] < 300)
ck("postprocess 应激 堆内存 < 2GB", iso_mem_mb() < 2000)
# 300 行子集参照抽查
set.seed(99)
sub <- sort(sample(seq_len(np_dec), 300L))
res3 <- pFind_res_reader(qry3)
sp3_ref <- read.delim(file.path(mt3, "pFind.spectra"), sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      quote = "", comment.char = "", na.strings = character(0))
fa3_ln <- enc2utf8(readLines(fa3_file, warn = FALSE))
hid3 <- grepl("^>", fa3_ln)
fa3_names <- vapply(strsplit(sub("^>", "", fa3_ln[hid3]), " "), function(x) x[1], "")
grp3 <- cumsum(hid3)
fa3_seqs <- as.character(tapply(fa3_ln[!hid3], grp3[!hid3],
                                function(x) paste(x, collapse = "")))
sp3_ref2 <- ref_postprocess(sp3_ref, res3, fa3_names, fa3_seqs, rows = sub)
ref3_file <- file.path(mt3, "ref3_out.spectra")
write.table(sp3_ref2, ref3_file, sep = "\t", row.names = FALSE, quote = FALSE,
            fileEncoding = "UTF-8", na = "")
sp3_ref3 <- read.delim(ref3_file, sep = "\t", stringsAsFactors = FALSE,
                       check.names = FALSE, quote = "", comment.char = "",
                       na.strings = character(0))
sp3_new <- read.delim(file.path(mt3, "pFind-Filtered.spectra"), sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      quote = "", comment.char = "", na.strings = character(0))
ck("postprocess 应激 300 行子集与逐行参照 identical",
   identical(sp3_new[sub, , drop = FALSE], sp3_ref3[sub, , drop = FALSE]))

# ============================================================================
# 6. reader PSOCK 并行（2026-09-24 R3 P0）：并行 vs 串行逐元素 identical，
#    含跨文件同标题覆盖顺序语义（后文件胜出）
# ============================================================================
mt4 <- tempfile("r2par_"); dir.create(mt4)
mk_pfile <- function(tag, blocks) {
  writeLines(unlist(blocks), file.path(mt4, tag), useBytes = TRUE)
}
# 复用前文 mk_block（单块＝S 行+标题+候选行）
p1 <- mk_block("pa001.dta", c("PEPTIDEK", "PEPTIDEK"), c(60, 90),
               decoy_flags = c(TRUE, FALSE))
p2 <- mk_block("pa002.dta", c("ANOTHERR", "ANOTHERR"), c(70, 95),
               decoy_flags = c(TRUE, FALSE))
p3 <- mk_block("pa003.dta", c("MFLAGK", "MFLAGK"), c(50, 85),
               decoy_flags = c(TRUE, FALSE))
p4 <- mk_block("pa004.dta", c("ZZZOTHERR", "ZZZOTHERR"), c(40, 80),
               decoy_flags = c(TRUE, FALSE))
p5 <- mk_block("pa005.dta", c("CCCCCCCC", "CCCCCCCC"), c(30, 75),
               decoy_flags = c(TRUE, FALSE))
p6 <- mk_block("pa006.dta", c("WYWYWYWY"), c(88))
# 文件 3 里重复 pa001 标题、候选不同——合并后必须以文件 3 版本为准（覆盖顺序）
p7 <- mk_block("pa001.dta", c("QQQQQQQQ"), c(99), decoy_flags = c(FALSE))
p8 <- mk_block("pa007.dta", c("DDDDDDDD", "DDDDDDDD"), c(20, 77),
               decoy_flags = c(TRUE, FALSE))
mk_pfile("f1.qry.res", list(p1, p2, p3))
mk_pfile("f2.qry.res", list(p4, p5, p6))
mk_pfile("f3.qry.res", list(p7, p8))
files4 <- file.path(mt4, c("f1.qry.res", "f2.qry.res", "f3.qry.res"))
keep4 <- c("pa001.dta", "pa002.dta", "pa005.dta", "pa007.dta")
rs <- pFind_res_reader(files4, keep_spectra = keep4,
                       drop_decoy_candidates = TRUE)
rp <- pFind_res_reader(files4, keep_spectra = keep4,
                       drop_decoy_candidates = TRUE, parallel = TRUE,
                       n_workers = 3L)
ck("reader 并行 vs 串行逐元素 identical", identical(rp, rs))
ck("reader 并行 keep 集合过滤一致", identical(sort(names(rp)), sort(keep4)))
ck("reader 跨文件同标题覆盖顺序语义（后文件胜出）",
   identical(rp[["pa001.dta"]]$candidates[[1]]$sequence, "QQQQQQQQ"))
rm(rs, rp)

cat(sprintf("=== r2_memory_fix_test 汇总：%d PASS / %d FAIL ===\n", npass, nfail))
unlink(td, recursive = TRUE); unlink(mt, recursive = TRUE)
unlink(mt2, recursive = TRUE); unlink(mt3, recursive = TRUE)
unlink(mt4, recursive = TRUE)
if (nfail > 0L) quit(status = 1L)
