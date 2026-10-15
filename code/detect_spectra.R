#  全场景谱图入口检测
# 按优先级依次检测，覆盖文件夹/文件列表 × raw/mgf/混合 × MaxQuant raw 警告
# 全场景检测
# 输入为文件夹：查 mgf+raw→两者有则保存双值；仅 raw→pParse 转；仅 mgf→MaxQuant 警告；空→报错
# 输入为文件列表：raw→先查同文件夹 mgf→无则 pParse 转；mgf→MaxQuant 警告
# mgf 键空：从 psm_file 旁边找（ 已有）
detect_spectra_entries <- function(mgf_entries, cfg) {
  ent <- mgf_entries[trimws(mgf_entries) != ""]

  # ── 分类输入条目 ──
  direct_raws  <- ent[grepl("\\.raw$",  ent, ignore.case = TRUE) & file.exists(ent)]
  direct_mgfs  <- ent[grepl("\\.mgf$",  ent, ignore.case = TRUE) & file.exists(ent)]
  other_files  <- ent[file.exists(ent) & !grepl("\\.(raw|mgf)$", ent, ignore.case = TRUE)]
  dirs         <- ent[dir.exists(ent)]

  # 文件夹展开：同时找 mgf 和 raw
  dir_mgfs <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = "\\.mgf$", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)),
    use.names = FALSE)
  dir_raws <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = "\\.raw$", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)),
    use.names = FALSE)

  all_mgfs <- unique(gsub("\\\\", "/", normalizePath(c(direct_mgfs, dir_mgfs), winslash = "/")))
  all_raws <- unique(gsub("\\\\", "/", normalizePath(c(direct_raws, dir_raws), winslash = "/")))

  # ── 检测流程是否涉及 MaxQuant ──
  mq_needed <- {
    se <- if (!is.null(cfg$search_engine)) tolower(trimws(tail(cfg$search_engine, 1))) else ""
    # 检查 config_step2 是否已存在且含 maxquant
    cs2 <- file.path(cfg$config_dir, "config_step2.txt")
    cs2_mq <- FALSE
    if (file.exists(cs2)) {
      cs2_lines <- readLines(cs2, warn = FALSE)
      cs2_mq <- any(grepl("maxquant", cs2_lines, ignore.case = TRUE))
    }
    identical(se, "maxquant") || cs2_mq
  }

  # ── 检查 pparse_mgf 已有转换产物（指纹复用） ──
  pp_dir <- file.path(cfg$config_dir, "pparse_mgf")
  pp_mgfs <- if (dir.exists(pp_dir))
    list.files(pp_dir, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE) else character(0)

  # ── 优先级 1：MGF 可用（直接用 + 保留 raw 给 MaxQuant） ──
  mgf_result <- character(0)
  raws_for_mq <- character(0)

  if (length(all_mgfs) > 0) {
    mgf_result <- all_mgfs
    cat("[第一步] 检测到 ", length(all_mgfs), " 个 MGF 文件\n", sep = "")
  }

  # ── 优先级 2：raw 存在 → 先查同文件夹已有 MGF，再补转换 ──
  if (length(all_raws) > 0) {
    raws_for_mq <- all_raws
    cat("[第一步] 检测到 ", length(all_raws), " 个 raw 文件（供 MaxQuant 使用）\n", sep = "")

    # 检查哪些 raw 已有对应 MGF（同文件夹或 pparse_mgf 产物）
    raw_stems  <- tools::file_path_sans_ext(basename(all_raws))
    mgf_stems  <- tools::file_path_sans_ext(basename(c(all_mgfs, pp_mgfs)))
    need_conv  <- all_raws[!raw_stems %in% mgf_stems]

    if (length(need_conv) > 0) {
      cat("[第一步] ", length(need_conv), " 个 raw 无对应 MGF，启动 pParse 转换\n", sep = "")
      # 调用已有转换链（ 的 pParse/msconvert 回退）
      converted <- pp_convert(need_conv, cfg)
      mgf_result <- c(mgf_result, converted)
    } else {
      cat("[第一步] 全部 raw 已有对应 MGF，跳过转换\n", sep = "")
    }
  }

  # ── 优先级 3：仅 MGF、无 raw，且 MaxQuant 需要 raw → 警告 ──
  if (mq_needed && length(all_raws) == 0) {
    cat("⚠ [第一步] MaxQuant 自动搜库需要 raw 文件，但 mgf 键中未检测到 raw。\n",
        "  请在 mgf 键补充 raw 文件路径或文件夹（与 MGF 可混合填写），\n",
        "  或在 config_step2 的 search_file_raw 键填写 raw 路径。\n",
        "  MaxQuant 无法读取 MGF——无 raw 时该引擎将跳过。\n",
        "  [step 1] MaxQuant auto-search needs raw files, but none detected.\n",
        "  Add raw file paths/folders to the mgf key (can mix with MGF),\n",
        "  or fill search_file_raw in config_step2.\n", sep = "")
  }

  # ── 优先级 4：什么都没有 → 报错 ──
  if (length(mgf_result) == 0 && length(all_raws) == 0 && length(other_files) == 0) {
    stop("[第一步] mgf 键填写了条目但未找到任何 .mgf 或 .raw 文件：\n  ",
         paste(ent, collapse = "\n  "), "\n",
         "请确认路径正确且文件夹内含有谱图文件。\n",
         "[step 1] mgf key has entries but no .mgf or .raw files found.\n",
         "Check the paths and ensure the folders contain spectrum files.")
  }

  # ── 自动传递 raw 路径给 MaxQuant（写 config_step2 的 search_file_raw） ──
  if (mq_needed && length(raws_for_mq) > 0) {
    cs2 <- file.path(cfg$config_dir, "config_step2.txt")
    raw_line <- paste(raws_for_mq, collapse = ";")
    if (file.exists(cs2)) {
      lines <- readLines(cs2, warn = FALSE)
      i_sfr <- grep("^search_file_raw\\t", lines)
      if (length(i_sfr) == 0) {
        lines <- c(lines, paste0("search_file_raw\t", raw_line))
      } else {
        lines[i_sfr] <- paste0("search_file_raw\t", raw_line)
      }
      writeLines(lines, cs2)
      cat("[第一步] 已自动写入 search_file_raw（", length(raws_for_mq), " 个 raw）→ MaxQuant\n", sep = "")
    } else {
      cat("[第一步] raw 路径已保存（", length(raws_for_mq),
          " 个），config_step2 生成时将自动传递给 MaxQuant\n", sep = "")
    }
  }

  list(mgf_files = sort(unique(mgf_result)),
       raw_files = raws_for_mq,
       mq_needed = mq_needed)
}

# pParse 转换辅助（从 convert_raw_entries_to_mgf 提取核心，供 detect_spectra_entries 调用）
pp_convert <- function(raws, cfg) {
  out_dir <- file.path(cfg$config_dir, "pparse_mgf")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pf_root <- if (!is.null(cfg$pfind_exe) && length(cfg$pfind_exe) > 0)
    trimws(tail(cfg$pfind_exe, 1)) else ""
  pparse <- ""
  if (nzchar(pf_root)) {
    cand <- file.path(gsub("\\\\", "/", pf_root), "bin", "pParse.exe")
    cand2 <- if (grepl("/bin$", cand)) file.path(dirname(cand), "pParse.exe") else cand
    pparse <- if (file.exists(cand)) cand else if (file.exists(cand2)) cand2 else ""
  }
  logf <- file.path(out_dir, "pparse_run.log")
  if (nzchar(pparse)) {
    cfg_lines <- c("# auto-generated by isoPSV (step-1 raw->MGF)",
                   "[Basic Options]",
                   paste0("datanum = ", length(raws)),
                   paste0("datapath", seq_along(raws), " = ", raws),
                   "[Advanced Options]",
                   "input_format = raw", "output_mgf = 1", "output_pf = 0",
                   "output_msb = 0", "delete_msn = 0", "co-elute = 1",
                   paste0("outputpath = ", gsub("\\\\", "/", out_dir)),
                   "rewrite_files = 0")
    pcfg <- file.path(out_dir, "pparse_auto.cfg")
    writeLines(cfg_lines, pcfg, useBytes = TRUE, sep = "\r\n")
    iconv_lines <- iconv(cfg_lines, from = "UTF-8", to = "GBK")
    if (any(is.na(iconv_lines))) iconv_lines <- cfg_lines
    con <- file(pcfg, open = "wb"); writeLines(iconv_lines, con, sep = "\r\n"); close(con)
    old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
    setwd(dirname(pparse))
    st <- suppressWarnings(system2(pparse, shQuote(pcfg), stdout = logf, stderr = logf))
    setwd(old_wd)
    mgs <- list.files(out_dir, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE)
    if (!length(mgs))
      stop("[第一步] pParse 转换失败：", basename(raws[1]), "（日志：", logf, "）\n")
    cat("[第一步] pParse 转换完成：", length(mgs), " 个 MGF -> ", out_dir, "\n", sep = "")
    return(gsub("\\\\", "/", mgs))
  }
  # msconvert 回退
  msc <- if (!is.null(cfg$msconvert_exe) && length(cfg$msconvert_exe) > 0)
    trimws(tail(cfg$msconvert_exe, 1)) else ""
  if (nzchar(msc) && file.exists(msc)) {
    for (r in raws)
      suppressWarnings(system2(msc, c(shQuote(r), "--mgf", "-o", shQuote(out_dir)),
                               stdout = logf, stderr = logf))
    mgs <- list.files(out_dir, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE)
    if (length(mgs)) {
      cat("[第一步] msconvert 转换完成：", length(mgs), " 个 MGF\n", sep = "")
      return(gsub("\\\\", "/", mgs))
    }
  }
  stop("[第一步] raw 转换失败：未找到 pParse 或 msconvert\n",
       "[step 1] raw conversion failed: neither pParse nor msconvert found")
}
