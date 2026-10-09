# =============================================================================
# sysinfo_probe.R - 独立运行环境探测脚本（用户口径 2026-09-21：软件质控之一）
#
# 检测：CPU（型号/核数/线程）、内存（总量/可用）、磁盘（各盘剩余）、GPU、
#       R 与 Python 是否安装及版本、管线依赖的外部引擎路径是否存在。
# 用法：Rscript sysinfo_probe.R [输出目录（默认当前目录）]
# 输出：sysinfo_report.txt（双语人读）+ sysinfo_report.tsv（机器可读）
# 独立脚本：不 source 管线任何文件，不写配置，可单独分发。
# =============================================================================
OUT_DIR <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else getwd()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ps <- function(cmd) {
  r <- suppressWarnings(system2("powershell",
        c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = TRUE))
  paste(trimws(r[r != ""]), collapse = " | ")
}

lines <- c()
add <- function(zh, en, val) {
  lines <<- c(lines, sprintf("%s / %s: %s", zh, en, val))
}

# ---- CPU ----
cpu_name <- ps("(Get-CimInstance Win32_Processor | Select-Object -First 1).Name")  # 双路 CPU 取首路
cpu_cores <- ps("(Get-CimInstance Win32_Processor).NumberOfCores")
cpu_threads <- ps("(Get-CimInstance Win32_Processor).NumberOfLogicalProcessors")
add("CPU 型号", "CPU model", cpu_name)
add("物理核", "physical cores", cpu_cores)
add("逻辑线程", "logical threads", cpu_threads)

# ---- 内存 ----
mem_gb <- ps("[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)")
free_mem_gb <- ps("[math]::Round((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1MB,1)")
add("内存总量", "total RAM (GB)", mem_gb)
add("可用内存", "free RAM (GB)", free_mem_gb)

# ---- 磁盘 ----
disks <- strsplit(ps("(Get-CimInstance Win32_LogicalDisk -Filter 'DriveType=3' | ForEach-Object { $_.DeviceID + '=' + [math]::Round($_.FreeSpace/1GB,1) + '/' + [math]::Round($_.Size/1GB,1) + 'GB' })"), " \\| ")[[1]]
add("磁盘剩余/总量", "disk free/total", if (length(disks)) paste(disks, collapse="  ") else "n/a")

# ---- GPU ----
gpu <- ps("(Get-CimInstance Win32_VideoController | Select-Object -First 1).Name")
add("GPU", "GPU", gpu)

# ---- R / Python ----
rv <- paste(R.version$major, R.version$minor, sep = ".")
rv <- paste(R.version$major, R.version$minor, sep = ".")
rv_parts <- as.numeric(strsplit(rv, ".", fixed = TRUE)[[1]])
add("R 版本", "R version", rv)
add("R 版本判定", "R version verdict",
    if (rv_parts[1] > 4 || (rv_parts[1] == 4 && rv_parts[2] >= 1))
      "满足（≥4.1；4.3/4.4/4.5 均可）/ satisfied (>=4.1)" else
      "低于 4.1：部分基础函数可能缺失，建议升级 / below 4.1: upgrade recommended")
py  <- tryCatch(suppressWarnings(system2("python", c("--version"), stdout = TRUE, stderr = TRUE)),
                error = function(e) character(0))   # 未安装时 system2 抛 error 而非 warning
py3 <- tryCatch(suppressWarnings(system2("py", c("-3", "--version"), stdout = TRUE, stderr = TRUE)),
                error = function(e) character(0))
py4 <- tryCatch(suppressWarnings(system2("python3", c("--version"), stdout = TRUE, stderr = TRUE)),
                error = function(e) character(0))   # F-6（2026-09-26）：与 py_detect 三段回退对齐
add("Python (python)", "Python via python", if (length(py)) trimws(py[1]) else "未安装 / not installed")
add("Python (py -3)", "Python via py launcher", if (length(py3)) trimws(py3[1]) else "未安装 / not installed")
add("Python (python3)", "Python via python3", if (length(py4)) trimws(py4[1]) else "未安装 / not installed")

# ---- 外部引擎存在性（常见安装位置探测；未配置路径时报未检出而非错误）----
# 2026-09-24 R3 P3：固定路径之外，按安装根做有界深度模式搜索（覆盖工作站
# C:/software 下的新装布局 MaxQuant/MSFragger/Comet/Sage/MetaMorpheus/MSBooster/
# X!Tandem）；固定路径命中优先（不触发深搜）。
# 2026-09-24 R3 P3: beyond the fixed paths, a bounded-depth pattern search
# under the install roots (covers the workstation C:/software layouts of
# MaxQuant/MSFragger/Comet/Sage/MetaMorpheus/MSBooster/X!Tandem); fixed-path
# hits win (no deep walk).
ENG <- c(
  pFind        = "D:/Bioinformatics/pFindStudio",
  pFind_alt    = "C:/software/pFind3",
  pDeep2       = "D:/Bioinformatics/pDeep2",
  pDeep2_alt   = "C:/software/pDeep2",
  pUniFind     = "D:/Bioinformatics/piUniMS/punifind",
  pUniFind_alt = "C:/software/pUniFind-main",
  MSFragger    = "D:/Bioinformatics/MSFragger",
  Comet        = "D:/Bioinformatics/Comet",
  Sage         = "D:/Bioinformatics/Sage",
  MaxQuant     = "D:/Bioinformatics/MaxQuant",
  MetaMorpheus = "D:/Bioinformatics/MetaMorpheus",
  ProteoWizard = "C:/Program Files/ProteoWizard")
ENG_ROOTS <- c("D:/Bioinformatics", "C:/software", "C:/Program Files",
               "C:/Program Files (x86)")
ENG_PAT <- c(
  pFind        = "pFind\\.exe$",
  pDeep2       = "predict\\.py$",
  pUniFind     = "official_score_workflow\\.exe$",
  MSFragger    = "(?i)MSFragger[-.][0-9].*\\.jar$",
  Comet        = "comet\\.win64\\.exe$",
  Sage         = "(^|[/\\\\])sage\\.exe$",   # B8：锚定段首，排除 usage.exe 类假阳性
  MaxQuant     = "MaxQuantCmd\\.exe$",
  MetaMorpheus = "(?i)MetaMorpheus.*\\.exe$",
  MSBooster    = "(?i)MSBooster.*\\.jar$",
  XTandem      = "tandem\\.exe$",
  ProteoWizard = "msconvert\\.exe$")
bounded_find <- function(rt, rx, maxd = 6L) {
  walk <- function(d, depth) {
    if (depth > maxd) return(character(0))
    fs <- list.files(d, full.names = TRUE, no.. = TRUE, include.dirs = TRUE)
    hit <- character(0)
    for (f in fs) {
      if (dir.exists(f)) {
        sub <- walk(f, depth + 1L)
        if (length(sub) > 0) return(sub[1])       # 首个命中即返回
      } else if (grepl(rx, basename(f), perl = TRUE)) {
        return(f)
      }
    }
    hit
  }
  tryCatch(walk(rt, 1L), error = function(e) character(0))
}
eng_status <- character(0)
all_eng <- unique(c(names(ENG), names(ENG_PAT)))
for (nm in all_eng) {
  ok <- FALSE; path <- ""
  if (nm %in% names(ENG)) {
    ok <- dir.exists(ENG[[nm]])
    path <- ENG[[nm]]
  }
  if (!ok && nm %in% names(ENG_PAT)) {
    for (rt in ENG_ROOTS) {
      f <- bounded_find(rt, ENG_PAT[[nm]])
      if (length(f) > 0 && nzchar(f)) { ok <- TRUE; path <- f; break }
    }
  }
  eng_status[[nm]] <- ok
  add(paste0("引擎 ", nm), paste0("engine ", nm),
      if (ok) paste0("检出 / found: ", path) else "未检出（不影响核心流程，仅该引擎不可用）/ not found")
}

# ---- 结论 ----
nc <- suppressWarnings(as.integer(cpu_threads))
nr <- suppressWarnings(as.numeric(gsub(",", ".", mem_gb)))
core_ok <- all(c(eng_status[["pFind"]], eng_status[["pDeep2"]],
                 eng_status[["pUniFind"]]))
verdict <- if (!is.na(nc) && !is.na(nr)) {
  if (nr >= 32 && nc >= 16 && core_ok) "满足全要素运行（核心引擎齐备）/ full run capable"
  else if (nr >= 16 && nc >= 8 && core_ok) "满足核心流程（外部引擎按需）/ core pipeline capable"
  else if (nr >= 16 && nc >= 8) "资源满足但核心引擎未齐（pFind/pDeep2/pUniFind 缺）——先补装引擎 / resources OK but core engines missing"
  else "低于最低建议（16GB/8线程），仅小数据集可用 / below minimum, small datasets only"
} else "无法判定 / undetermined"
add("结论", "verdict", verdict)
add("缺 R 时", "if R missing", "管线主体为 R：请安装 R ≥ 4.1（推荐 4.5.x；installer /VERYSILENT 静默装）或使用便携包（见部署说明 deploy_notes.md）")
add("缺 Python 时", "if Python missing", "核心 R 流程不受影响；pUniFind pkl 转换器与网页界面需要 Python 3.8+，可装或用便携嵌入包")

writeLines(lines, file.path(OUT_DIR, "sysinfo_report.txt"), useBytes = TRUE)
kv <- sub(": .*$", "", lines)   # 键＝首个 ": " 之前（原写法把值当成了键）
writeLines(paste0(kv, "\t", sub("^[^:]+: ", "", lines)),
           file.path(OUT_DIR, "sysinfo_report.tsv"), useBytes = TRUE)
cat(lines, sep = "\n")
cat("\nsysinfo_report.txt / .tsv written to: ", normalizePath(OUT_DIR), "\n")
