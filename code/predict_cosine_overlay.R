# predict_cosine_overlay.R — 多预测器余弦叠加到 master_table（v1.0）
# 在 unlabeled_master_export.R 之后运行：
#   1. 扫描 merge/predictors/<type>/predict.txt 的所有成功预测器
#   2. 读 MGF → spectrum_params → 各预测器 predict.txt → 余弦
#   3. 为 cand1 追加 cand1_predict_<type> 列（列名即软件名）
# 用法：Rscript predict_cosine_overlay.R <result_dir> <mgf_folder>
# 输出：master_table.tsv 新增列，原列不动
CODE <- NULL
for (.p in c(
  file.path(getwd(), "code"),
  file.path(getwd(), "..", "code"),
  getwd())) {
  if (file.exists(file.path(.p, "ReadMGF_script.R"))) { CODE <- .p; break }
}
if (is.null(CODE)) stop("ReadMGF_script.R not found from getwd: ", getwd())
source(file.path(CODE, "compat.R"), encoding = "UTF-8")
source(file.path(CODE, "ReadMGF_script.R"), encoding = "UTF-8")
source(file.path(CODE, "paired_params.R"), encoding = "UTF-8")

args <- commandArgs(FALSE)
arg_of <- function(flag) {
  hits <- args[startsWith(args, paste0("--", flag, "="))]
  if (length(hits) > 0) sub(paste0("^--", flag, "="), "", hits[1]) else ""
}
rd <- arg_of("result_dir")
if (!nzchar(rd)) stop("usage: Rscript predict_cosine_overlay.R --result_dir=<dir>")

aa <- read.delim(file.path(CODE, "AA_mass_dataframe.txt"), sep = "\t",
                 header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- data.frame(
    Full_Name = c("Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
                  "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
                  "Leucine","Lysine","Methionine","Phenylyalanine","Proline",
                  "Serine","Threonine","Tryptophan","Tyrosine","Valine","Ile_or_Leu"),
    Three_Letter = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
                     "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp",
                     "Tyr","Val","Xle"),
    One_Letter = c("A","R","N","D","C","Q","E","G","H","I","L","K","M",
                   "F","P","S","T","W","Y","V","X"))

mt_path <- file.path(rd, "1.unlabeled_result", "master_table.tsv")
if (!file.exists(mt_path)) stop("master_table.tsv not found: ", mt_path)
mt <- read.delim(mt_path, sep = "\t", check.names = FALSE,
                 stringsAsFactors = FALSE, na.strings = character(0),
                 colClasses = "character")
cat("[predict_cosine_overlay] master_table: ", nrow(mt), " 行\n", sep = "")

# 扫描预测器
pred_dirs <- list.dirs(file.path(rd, "merge", "predictors"), recursive = FALSE)
pred_types <- character(0)
for (pd in pred_dirs) {
  ptxt <- file.path(pd, "predict.txt")
  if (file.exists(ptxt)) pred_types <- c(pred_types, basename(pd))
}
# pDeep2 产物在 merge/predict.txt（旧布局）或 predictors/pdeep2/（新布局）
pd2_new <- file.path(rd, "merge", "predictors", "pdeep2", "predict.txt")
pd2_old <- file.path(rd, "merge", "predict.txt")
if (file.exists(pd2_old) && !("pdeep2" %in% pred_types)) pred_types <- c("pdeep2", pred_types)
cat("[predict_cosine_overlay] 预测器: ", paste(pred_types, collapse = ", "), "\n", sep = "")
if (length(pred_types) == 0) { cat("无预测器产物，跳过\n"); quit(status = 0) }

# 读各预测器 predict.txt
preds <- list()
for (pt in pred_types) {
  pf <- if (pt == "pdeep2" && file.exists(pd2_old)) pd2_old else
        file.path(rd, "merge", "predictors", pt, "predict.txt")
  pr <- try(read_predict_txt(pf), silent = TRUE)
  if (!inherits(pr, "try-error") && !is.null(pr$idx)) {
    preds[[pt]] <- pr
    cat("  ", pt, ": ", length(pr$int), " entries\n", sep = "")
  }
}

# MGF 位置：rd 可能是 result_* 时间戳子目录——先试 rd/unlabeled，再试
# dirname(rd)/unlabeled（项目根），最后读 config_step1 的 mgf 键
mgf_dir <- file.path(rd, "unlabeled")
mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
if (length(mgf_files) == 0) {
  mgf_dir <- file.path(dirname(rd), "unlabeled")
  mgf_files <- sort(list.files(mgf_dir, pattern = "\\.mgf$", full.names = TRUE))
}
if (length(mgf_files) == 0) {
  # 从 config_step1 读 mgf 键
  cfg1_path <- file.path(dirname(rd), "config_step1.txt")
  if (file.exists(cfg1_path)) {
    cfg1 <- read_param_file(cfg1_path)
    mgf_key <- if (!is.null(cfg1$mgf)) tail(cfg1$mgf, 1) else ""
    if (nzchar(mgf_key) && dir.exists(mgf_key))
      mgf_files <- sort(list.files(mgf_key, pattern = "\\.mgf$", full.names = TRUE))
  }
}
cat("[predict_cosine_overlay] MGF: ", length(mgf_files), " 个\n", sep = "")
if (length(mgf_files) == 0) {
  cat("[predict_cosine_overlay] 未找到 MGF，跳过余弦计算（master_table 已含可用性列）\n")
  cat("[predict_cosine_overlay] no MGF found, skipping (master_table already has availability columns)\n")
  quit(status = 0)
}

# 建 File_Name → MGF spectrum 映射
spec_map <- new.env(parent = emptyenv())
for (mf in mgf_files) {
  mgfs <- parse_lines_vectorized(readLines(mf, warn = FALSE))
  for (s in mgfs) assign(s@title, s, envir = spec_map)
  rm(mgfs); gc()
}
cat("[predict_cosine_overlay] 谱图索引: ", ls(spec_map), " 条\n", sep = "")

# 逐行计算
n_new_cols <- 0
for (pt in names(preds)) {
  col_name <- paste0("cand1_predict_", pt)
  if (col_name %in% names(mt)) next
  vals <- rep("", nrow(mt))
  pr <- preds[[pt]]
  n_hit <- 0
  for (i in seq_len(nrow(mt))) {
    fn <- mt$File_Name[i]
    seq_p <- mt$cand1_sequence[i]
    mod_p <- mt$cand1_modification[i]
    if (!nzchar(seq_p) || is.na(seq_p)) next
    mod_p <- sub(";+$", "", as.character(mod_p))
    if (is.na(mod_p)) mod_p <- ""
    # 查预测器
    key1 <- paste0(seq_p, "|", mod_p, "|", mt$Charge[i])
    pi <- unname(pr$idx[key1])
    if (is.na(pi)) pi <- unname(pr$idx[paste0(seq_p, "||", mt$Charge[i])])
    if (is.na(pi)) next
    # 查谱图
    spc <- get(fn, envir = spec_map, inherits = FALSE)
    if (is.null(spc)) next
    # spectrum_params
    sp <- try(spectrum_params(seq_p, mod_p, "unlabeled", aa, aa_tab,
                              spc@spectrum, 20, "ppm", labels = NULL),
              silent = TRUE)
    if (is.null(sp) || inherits(sp, "try-error")) next
    v <- sp$b_y_ions_intensity
    pv <- pr$int[[pi]]
    nmin <- min(length(pv), length(v))
    if (nmin <= 0) next
    vals[i] <- sprintf("%.4f", cosine(v[seq_len(nmin)], pv[seq_len(nmin)]))
    n_hit <- n_hit + 1
  }
  mt[[col_name]] <- vals
  n_new_cols <- n_new_cols + 1
  cat("  ", col_name, ": ", n_hit, "/", nrow(mt), " 覆盖\n", sep = "")
}

# 写回
if (n_new_cols > 0) {
  write.table(mt, mt_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[predict_cosine_overlay] 已写回 ", n_new_cols, " 列 -> ", mt_path, "\n", sep = "")
} else {
  cat("[predict_cosine_overlay] 无新增列\n")
}
