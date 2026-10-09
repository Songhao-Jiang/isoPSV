
# 构建 SILAC 标记修饰字符串（与 4.heavy_spectra_match_SILAC.R 一致）
#   Build SILAC labeling modification string (mirrors 4.heavy_spectra_match_SILAC.R);
# p_model = FALSE → zeroP（无 P 校正），TRUE → fullP（加 P_C5N1）
#   p_model = FALSE -> zeroP (no P correction), TRUE -> fullP (P_C5N1 added)
# 输入：
#   seq: 肽段序列（字符）；mod: 初始修饰字符串（字符，可为空）
#   p_model: 是否追加脯氨酸 P_C5N1 修正（逻辑，默认 TRUE）
# 输出：
#   完整修饰字符串（字符，多个修饰片段以分号分隔拼接）
# Input:
#   seq: peptide sequence (character); mod: initial modification string (may be empty)
#   p_model: whether to add the P_C5N1 proline correction (logical, default TRUE)
# Output:
#   full modification string (character, segments joined with semicolons)
build_silac_mod <- function(seq, mod, p_model = TRUE) {
  # 初始化输出修饰串，并定义按位置追加修饰片段的小工具函数
  # Initialize the output string and define a helper that appends a (pos, name) segment
  mod_silac <- mod
  add_mod <- function(m, pos, name) {
    seg <- paste0(pos, ",", name, ";")
    if (nchar(m) == 0) seg else paste0(m, seg)
  }
  # 对序列中的每个 K 位置追加 13C(6) 标记修饰
  # Append the 13C(6) label at every K position in the sequence
  if (str_detect(seq, "K")) {
    locs <- as.numeric(str_locate_all(seq, "K")[[1]][, 1])
    for (pos in locs) mod_silac <- add_mod(mod_silac, pos, "Label_13C(6)[K]")
  }
  # 对序列中的每个 R 位置追加 13C(6)15N(4) 标记修饰
  # Append the 13C(6)15N(4) label at every R position in the sequence
  if (str_detect(seq, "R")) {
    locs <- as.numeric(str_locate_all(seq, "R")[[1]][, 1])
    for (pos in locs) mod_silac <- add_mod(mod_silac, pos, "Label_13C(6)15N(4)[R]")
  }
  # p_model 开启时对每个 P 位置追加 P_C5N1 修正
  # Add the P_C5N1 correction at every P position when p_model is enabled
  if (p_model && str_detect(seq, "P")) {
    locs <- as.numeric(str_locate_all(seq, "P")[[1]][, 1])
    for (pos in locs) mod_silac <- add_mod(mod_silac, pos, "P_C5N1")
  }
  mod_silac
}
