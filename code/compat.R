# =============================================================================
# compat.R - base-R implementations of the stringr / reshape2 functions used
# by the client code, keeping the ORIGINAL function names so the calling
# code does not need to change. Source this file before the pipeline
# scripts; the client then runs without the stringr / reshape2 external
# packages (parallel is a base package and stays available for the PSOCK
# workers).
# =============================================================================
#
# 本文件用 base R 重新实现客户端用到的 stringr / reshape2 函数，并保留原始
# 函数名，使调用代码无需改动；在流水线脚本之前 source 本文件后，客户端即可
# 在不安装 stringr / reshape2 外部包的情况下运行（parallel 为 base 包，
# PSOCK 并行工作进程仍可正常使用）。
# This file reimplements in base R the stringr / reshape2 functions used by
# the client, keeping the original names so the calling code needs no change;
# source it before the pipeline scripts and the client runs without the
# stringr / reshape2 external packages (parallel is a base package and stays
# available for the PSOCK workers).

# 检测字符串是否匹配正则模式（stringr::str_detect 的 base R 替代）
# Test whether strings match a regular expression (base-R replacement for stringr::str_detect)
# 输入：
#   string: 待检测字符串向量；pattern: 正则模式；negate: 是否取反（逻辑）
# 输出：
#   逻辑向量：匹配（或按 negate 取反后）为 TRUE
# Input:
#   string: character vector to test; pattern: regular expression; negate: negate the result (logical)
# Output:
#   logical vector: TRUE where matched (or where negated)
str_detect <- function(string, pattern, negate = FALSE) {
  # 用 base R 的 grepl 完成匹配
  # Use base R grepl() for matching
  out <- grepl(pattern, string)
  # 需要取反时整体翻转结果
  # Flip the result when negation is requested
  if (negate) !out else out
}

# 按位置提取子串（stringr::str_sub 的 base R 替代，支持负索引）
# Extract substrings by position (base-R replacement for stringr::str_sub, negative indices supported)
# 输入：
#   string: 字符串；start / end: 起止位置（1 基；负值表示从末尾倒数）
# 输出：
#   提取出的子串向量
# Input:
#   string: character string; start / end: positions (1-based; negatives count from the end)
# Output:
#   vector of extracted substrings
str_sub <- function(string, start = 1L, end = -1L) {
  # 负索引换算为正索引（负数从末尾倒数）
  # Convert negative indices to positive ones (negatives count from the end)
  n <- nchar(string)
  s <- ifelse(start >= 0, start, n + start + 1)
  e <- ifelse(end >= 0, end, n + end + 1)
  # 调用 substr 完成截取
  # Call substr() to do the extraction
  substr(string, s, e)
}

# str_sub 的替换形式 str_sub(x, s, e) <- value（供 revise_search_result 使用）
# replacement form str_sub(x, s, e) <- value (used by revise_search_result)
# str_sub 的替换版本：把字符串指定区段替换为新值（保留赋值接口）
# Replacement form of str_sub: assign a new value into a substring region (keeps the assignment interface)
# 输入：
#   string: 原始字符串；start / end: 待替换区段位置；value: 新值
# 输出：
#   替换后的字符串
# Input:
#   string: original string; start / end: region to replace; value: new value
# Output:
#   the string after replacement
`str_sub<-` <- function(string, start = 1L, end = -1L, value) {
  # 与 str_sub 相同的负索引换算
  # Same negative-index conversion as str_sub
  n <- nchar(string)
  s <- ifelse(start >= 0, start, n + start + 1)
  e <- ifelse(end >= 0, end, n + end + 1)
  # 用 substr<- 就地替换区段后返回整串
  # Replace the region in place with substr<- and return the whole string
  substr(string, s, e) <- value
  string
}

# 按模式切分字符串（stringr::str_split 的 base R 替代）
# Split strings by a pattern (base-R replacement for stringr::str_split)
# 输入：
#   string: 待切分字符串；pattern: 切分模式；n: 切分块数上限；simplify: 是否简化为矩阵
# 输出：
#   simplify = FALSE 时为字符向量列表，TRUE 时为字符矩阵
# Input:
#   string: strings to split; pattern: splitting pattern; n: max pieces; simplify: return a matrix
# Output:
#   list of character vectors when simplify = FALSE, character matrix when TRUE
str_split <- function(string, pattern, n = Inf, simplify = FALSE) {
  # 模式不含正则元字符时按固定字符串处理，避免特殊字符被误当正则
  # Treat the pattern as fixed text when it has no regex metacharacters
  fixed <- !grepl("[][\\^$.*+?(){}|]", pattern)
  # simplify 开启时按行合并成矩阵
  # Combine into a matrix by rows when simplify is on
  if (simplify) {
    return(do.call(rbind, strsplit(as.character(string), pattern, fixed = fixed)))
  }
  # 默认返回按字符串切分出的向量列表
  # Return a list of split character vectors by default
  strsplit(as.character(string), pattern, fixed = fixed)
}

# 用替换串替换所有匹配（stringr::str_replace_all 的 base R 替代）
# Replace all matches with a replacement (base-R replacement for stringr::str_replace_all)
# 输入：
#   string: 字符串；pattern: 正则模式；replacement: 替换串
# 输出：
#   替换后的字符串
# Input:
#   string: character string; pattern: regular expression; replacement: replacement string
# Output:
#   the string after replacement
str_replace_all <- function(string, pattern, replacement) {
  # 使用 perl 正则的 gsub 做全局替换
  # Use gsub() with perl regex for the global replacement
  gsub(pattern, replacement, string, perl = TRUE)
}

# 定位首个匹配的起止位置（stringr::str_locate 的 base R 替代）
# Locate the start/end of the first match (base-R replacement for stringr::str_locate)
# 输入：
#   string: 字符串；pattern: 正则模式
# 输出：
#   两列矩阵（start, end）；无匹配时 start 为 -1
# Input:
#   string: character string; pattern: regular expression
# Output:
#   two-column matrix (start, end); start is -1 when there is no match
str_locate <- function(string, pattern) {
  # regexpr 找到首个匹配的位置与长度
  # regexpr() finds the position and length of the first match
  m <- regexpr(pattern, string, perl = TRUE)
  # 起止位置拼成两列矩阵
  # Combine start/end into a two-column matrix
  cbind(start = as.integer(m),
        end = as.integer(m) + attr(m, "match.length") - 1L)
}

# 定位全部匹配的起止位置（stringr::str_locate_all 的 base R 替代）
# Locate start/end of all matches (base-R replacement for stringr::str_locate_all)
# 输入：
#   string: 字符串向量；pattern: 正则模式
# 输出：
#   每个元素一个两列矩阵（start, end）；无匹配时为零行矩阵
# Input:
#   string: character vector; pattern: regular expression
# Output:
#   one two-column matrix (start, end) per element; zero-row matrix when there is no match
str_locate_all <- function(string, pattern) {
  # gregexpr 返回每个元素全部匹配的位置
  # gregexpr() returns all match positions for each element
  m <- gregexpr(pattern, string, perl = TRUE)
  lapply(seq_along(m), function(i) {
    # 无匹配（-1）时返回空的两列矩阵
    # Return an empty two-column matrix when there is no match (-1)
    if (m[[i]][1] == -1) return(matrix(integer(0), ncol = 2))
    # 组合全部匹配的起止位置
    # Combine start/end positions of all matches
    cbind(start = as.integer(m[[i]]),
          end = as.integer(m[[i]]) + attr(m[[i]], "match.length") - 1L)
  })
}

# 提取首个匹配及其捕获组（stringr::str_match 的 base R 替代）
# Extract the first match and its capture groups (base-R replacement for stringr::str_match)
# 输入：
#   string: 字符串；pattern: 含捕获组的正则模式
# 输出：
#   两列矩阵（完整匹配, 第一个捕获组）；无匹配时为 NA
# Input:
#   string: character string; pattern: regex with capture groups
# Output:
#   two-column matrix (full match, first capture group); NA when there is no match
str_match <- function(string, pattern) {
  # regexpr 给出匹配位置与捕获组起止
  # regexpr() gives match positions and capture-group spans
  m <- regexpr(pattern, string, perl = TRUE)
  s <- attr(m, "capture.start")
  l <- attr(m, "capture.length")
  # 分别截取完整匹配与第一个捕获组
  # Extract the full match and the first capture group separately
  cbind(ifelse(m > 0, substr(string, m, m + attr(m, "match.length") - 1L), NA),
        ifelse(s > 0, substr(string, s, s + l - 1L), NA))
}

# 提取全部匹配及其捕获组（stringr::str_match_all 的 base R 替代）
# Extract all matches and their capture groups (base-R replacement for stringr::str_match_all)
# 输入：
#   string: 字符串向量；pattern: 含捕获组的正则模式
# 输出：
#   每个元素一个矩阵（首列为完整匹配，其余列为捕获组）
# Input:
#   string: character vector; pattern: regex with capture groups
# Output:
#   one matrix per element (first column full match, remaining columns capture groups)
str_match_all <- function(string, pattern) {
  # gregexpr 返回全部匹配位置及捕获组信息
  # gregexpr() returns all match positions and capture-group info
  m <- gregexpr(pattern, string, perl = TRUE)
  lapply(m, function(mm) {
    # 无匹配时返回零列矩阵
    # Return a zero-column matrix when there is no match
    if (mm[1] == -1) return(matrix(character(0), ncol = 0))
    s <- attr(mm, "capture.start")
    l <- attr(mm, "capture.length")
    # 首列为完整匹配，其余列为各捕获组
    # First column full match, remaining columns the capture groups
    cbind(substr(string, mm, mm + attr(mm, "match.length") - 1L),
          matrix(ifelse(s > 0, substr(string, s, s + l - 1L), NA), ncol = ncol(s)))
  })
}

# reshape2::melt 的无 id.vars 数据框版本：按列（列主序）堆叠成
# reshape2::melt for a data.frame without id.vars: stack the columns
# 两列数据框（variable, value）
# (column-major) into a two-column data.frame (variable, value)
# 数据框长表化（reshape2::melt 的 base R 替代，不支持 id.vars）
# Melt a data.frame to long form (base-R replacement for reshape2::melt, no id.vars support)
# 输入：
#   data: 数据框；na.rm: 是否丢弃 NA（仅保留接口兼容）；value.name: 数值列名
# 输出：
#   两列数据框（variable, value）
# Input:
#   data: data.frame; na.rm: drop NA (kept only for interface compatibility); value.name: value column name
# Output:
#   two-column data.frame (variable, value)
melt <- function(data, ..., na.rm = FALSE, value.name = "value") {
  # 数据框分支：堆叠为长表
  # data.frame branch: stack into long form
  if (is.data.frame(data)) {
    # 构造 variable 因子并按列主序拼接数值
    # Build the variable factor and stack values in column-major order
    vars <- names(data)
    n <- nrow(data)
    out <- data.frame(
      variable = factor(rep(vars, each = n), levels = vars),
      value = unlist(lapply(data, as.vector), use.names = FALSE),
      stringsAsFactors = FALSE)
    # 数值列重命名为 value.name
    # Rename the value column to value.name
    names(out)[2] <- value.name
    return(out)
  }
  # 非数据框输入直接报错
  # Stop with an error for non-data.frame input
  stop("melt: unsupported input (data.frame expected)")
}
# 编码容错 CSV 读取UTF-8 BOM / UTF-8 / 本地编码（GBK 等）
# 依次判定，供搜索引擎适配器读取外部导出件（Mascot/PEAKS CSV 常带 BOM 或
# ANSI/GBK）。其余参数原样传给 read.csv。
# Encoding-tolerant CSV read tries UTF-8 BOM, UTF-8, then
# the native encoding (GBK etc.), for the engine adapters reading external
# exports (Mascot/PEAKS CSVs often carry a BOM or ANSI/GBK). Extra arguments
# pass through to read.csv unchanged.
read_csv_tolerant <- function(path, ...) {
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) stop("empty input: ", path)
  # N-3编码判定只需文件头——只读前 64 KB
  fh <- file(path, "rb"); on.exit(close(fh), add = TRUE)
  raw <- readBin(fh, "raw", n = min(sz, 65536L))
  has_bom <- length(raw) >= 3 && identical(raw[1:3], as.raw(c(0xEF, 0xBB, 0xBF)))
  # N-3 配套64 KB 截断点可能切在多字节字符中间，此时
  # validUTF8 误判为 FALSE 而错误落入本地编码分支（GBK 代码页 R ≤ 4.1 上
  # UTF-8 导出件会整篇乱码）。校验前先退到最后一个完整字符边界再判定。
  # The 64 KB cutoff may split a multi-byte character; back off to the last
  # complete character before the validUTF8 check.
  txt <- ""
  if (length(raw)) {
    for (k in 0:3) {
      h <- if (k) raw[seq_len(length(raw) - k)] else raw
      t0 <- tryCatch(suppressWarnings(rawToChar(h)), error = function(e) "")
      if (nzchar(t0) && isTRUE(all(validUTF8(t0)))) { txt <- t0; break }
    }
    if (!nzchar(txt)) txt <- tryCatch(suppressWarnings(rawToChar(raw)), error = function(e) "")
  }
  if (has_bom) return(read.csv(path, fileEncoding = "UTF-8-BOM", ...))
  if (nzchar(txt) && isTRUE(all(validUTF8(txt))))
    return(read.csv(path, fileEncoding = "UTF-8", ...))
  read.csv(path, ...)   # 本地编码（如 GBK 控制台导出件）
}

# Python 解释器探测
# python → py -3 启动器 → python3 回退；返回 list(cmd, pre) 或 NULL。
# 旧写法硬用 "python"，官方安装器/商店版常只有 py 启动器在 PATH 上，
# renumber_mgf.py 静默失败且返回码被忽略。
# Python interpreter probe (P1-4,; shared here per F-6, ):
# python -> py -3 launcher -> python3 fallback; returns list(cmd, pre) or NULL.
# The old code hardcoded "python", but installer/store Python often exposes only
# the py launcher on PATH, so renumber_mgf.py silently failed with its exit code
# ignored. Shared by user_step2.R and sysinfo_probe.R so the environment probe
# and the pipeline always agree on the interpreter.
py_detect <- function(hint = NULL) {
  # 自带 runtime 优先锁定——与 start_isoPSV.bat 的
  # "bundled portable runtimes take precedence" 同规则（网页层一直如此，
  # 管线此前唯独依赖 PATH，PATH 首个 python 损坏时报 initfsencoding/
  # No module named 'encodings' 而流程中断）。自带 python 为嵌入式 3.12
  # 纯标准库全隔离——本函数消费方（renumber_mgf.py 等）均 stdlib-only；
  # 需三方包的场景走各自专属键（如 punifind_python 的 torch 环境），
  # 不经本函数。runtime 缺失（开发树/裸副本）时按下方 hint → PATH 链回退。
  #  the bundled runtime wins (same rule as
  # start_isoPSV.bat's "bundled portable runtimes take precedence" - the web
  # layer always did this; the pipeline alone used PATH and aborted when the
  # first PATH python was broken). The bundled python is an embedded,
  # stdlib-only, fully-isolated 3.12 - every consumer of this function is
  # stdlib-only; third-party needs use their dedicated keys (e.g.
  # punifind_python's torch env) and never come here. When runtime is absent
  # (dev tree / bare copies), fall back to the hint -> PATH chain below.
  # 探测锚点＝代码目录（调用方 user_step2.R / sysinfo_probe.R 均以 code 为
  # 工作目录）的上一级包根。
  bund <- file.path(dirname(getwd()), "runtime", "python", "python.exe")
  if (file.exists(bund)) return(list(cmd = bund, pre = character(0)))
  # v1.0python_exe 显式键（runtime 缺失布局时的覆盖）。
  # 用户 Python 未写入环境目录（PATH）时，config_step1.txt 的 python_exe
  # 指向解释器完整路径即可；hint 指向的文件不存在时报消息并回退 PATH。
  # v1.0 the explicit python_exe key (override for layouts
  # without the bundled runtime). When the user's Python is NOT on PATH,
  # config_step1.txt's python_exe points at the full interpreter path; a hint
  # that does not exist is reported and PATH probing is used as fallback.
  hp <- if (is.null(hint)) "" else as.character(hint)[1]
  if (is.na(hp)) hp <- ""
  hp <- trimws(hp)
  if (nzchar(hp)) {
    if (file.exists(hp)) return(list(cmd = hp, pre = character(0)))
    message("python_exe 指向的 Python 不存在（", hp, "），回退环境目录探测\n",
            "python_exe points to a missing Python (", hp,
            "), falling back to PATH probing")
  }
  if (nzchar(Sys.which("python"))) return(list(cmd = "python", pre = character(0)))
  if (nzchar(Sys.which("py"))) return(list(cmd = "py", pre = "-3"))
  if (nzchar(Sys.which("python3"))) return(list(cmd = "python3", pre = character(0)))
  NULL
}

# 失败侧日志回贴引擎/步骤失败在主日志里只有一行
# exit code 时，真正的错因藏在各自的侧日志里（Mch_TL 实测：sage 自动段逐文件
# "引擎失败（exit 1）"，而 sage_exe 指错这一事实只在侧日志里可见）。把侧日志
# 尾部回贴主日志——与  的 renumber 回贴同格式、同诊断目的。日志可能是
# GBK 输出，按字节读取直通（不做转码，避免多字节标点编码错配）。
# Echo the tail of a step/engine side log into the main log when a failure is
# otherwise reported as a bare exit code (same format and purpose as the
#  renumber echo). Logs may be GBK; read as bytes and pass through
# (no transcoding - avoids multi-byte punctuation encodes mismatches).
log_fail_tail <- function(logf, tag = "", n = 8L) {
  # file.exists 前置守卫：文件缺失时 readLines 会经 warning 通道漏一条
  # "cannot open file" 到主日志（tryCatch 只收 error）——先判存在再读。
  lines <- if (file.exists(logf))
    tryCatch(tail(readLines(logf, warn = FALSE, encoding = "bytes"), n),
             error = function(e) character(0))
  else character(0)
  tryCatch({
    if (length(lines) > 0) {
      cat("[", tag, "] 失败日志尾部 / failure log tail (", logf, "):\n", sep = "")
      cat(paste0("  ", lines, collapse = "\n"), "\n", sep = "")
    } else {
      cat("[", tag, "] 日志为空或不可读 / log empty or unreadable: ", logf, "\n", sep = "")
    }
  }, error = function(e) NULL)
  invisible(NULL)
}

# 引擎适配器共用输入解析。
# 输入路径可为文件或文件夹：文件夹按模式递归找目标文件、多候选取最新 mtime、
# 无文件报错；文件路径直接返回。
# Shared input resolution for the engine adapters (: nine duplicated
# copies consolidated; variable names made self-documenting). The input path may
# be a file or a folder: a folder is searched recursively for files matching the
# pattern (newest mtime wins among multiple hits; error when none); a file path
# is returned as-is.
resolve_adapter_input <- function(input_path, file_pattern) {
  if (dir.exists(input_path)) {
    candidate_files <- unique(list.files(input_path, pattern = file_pattern,
                                         recursive = TRUE, full.names = TRUE))
    candidate_files <- candidate_files[!dir.exists(candidate_files)]
    if (length(candidate_files) == 0)
      stop("no file matching ", file_pattern, " found in folder: ", input_path)
    if (length(candidate_files) > 1)
      candidate_files <- candidate_files[which.max(file.mtime(candidate_files))]
    return(candidate_files)
  }
  if (file.exists(input_path)) return(input_path)
  stop("input missing: ", input_path)
}
