# =============================================================================
# paired_params.R - 成对谱参数系列与中间数据提取的合并函数
#   merged function for paired-spectrum parameter series and
# 基于 Fig.3 辅助代码 7.x.split_paired_search_PSM_*.R（用 all_ions_position +
#   intermediate data extraction (based on the Figure-3 helper code
# ions_match_selected）
#   7.x.split_paired_search_PSM_*.R, which uses all_ions_position +
#   ions_match_selected).
#
# 模式：unlabeled / SILAC / N15 / C13 / C13_comp
#   Modes: unlabeled / SILAC / N15 / C13 / C13_comp
#   C13_comp = C13 标记侧带质量补偿，等价于追加 -1 Da / -2 Da 修饰
#     C13_comp = C13 labeled side with mass compensation equivalent to adding
#     （1_13C / 2_13C），即额外的匹配窗口整体偏移
#                -1 Da / -2 Da modifications (1_13C / 2_13C), i.e. extra matching
#     -1.003355 与 -2.006711 Da（中性质量）
#                windows shifted by -1.003355 and -2.006711 Da (neutral).
#
# 依赖：ReadMGF 函数（all_ions_position、ions_match_selected）与 acs_core.R 的 build_silac_mod
#   Requires: ReadMGF functions loaded (all_ions_position, ions_match_selected)
#             and build_silac_mod (from acs_core.R).
# =============================================================================

# 计算单谱的理论离子匹配参数（b/y 覆盖度、PSM 得分、补偿字段等）
# Compute theoretical-ion match parameters for one spectrum (b/y coverage, PSM scores, compensation fields)
# 输入：
#   peptide / modification: 肽段序列与修饰
#   mode: 标记模式
#   AA_mass_dataframe / amino_acid_table: 氨基酸质量数据
#   spectrum: 目标谱（MGF 谱对象）
#   mass_shift_tolerence / ppm_or_da: 质量容差及单位
# 输出：
#   参数列表（ions_table、internal_aa_table、mgf_ions、b/y 强度与覆盖度、
#   PSM 得分、pep_mass、mode）；理论离子为空时返回 NULL
# Input:
#   peptide / modification: peptide sequence and modification
#   mode: labeling mode (unlabeled / SILAC / N15 / C13 / C13_comp)
#   AA_mass_dataframe / amino_acid_table: amino-acid mass data
#   spectrum: target spectrum (MGF spectrum object)
#   mass_shift_tolerence / ppm_or_da: mass tolerance and its unit
# Output:
#   parameter list (ions_table, internal_aa_table, mgf_ions, b/y intensities and coverage,
#   PSM scores, pep_mass, mode); NULL when the theoretical ions are empty
spectrum_params <- function(peptide, modification = "", mode = "unlabeled",
                            AA_mass_dataframe, amino_acid_table, spectrum,
                            mass_shift_tolerence = 20, ppm_or_da = "ppm",
                            labels = NULL, precursor_mh = NA_real_) {
  # 修饰为空/缺失时统一为空串
  # Normalize empty/NA modifications to ""
  if (is.null(modification) || is.na(modification)) modification <- ""
  # SILAC 标记作为 K/R/P 残基修饰处理
  # SILAC labeling handled as K/R/P residue modifications
  # 调 acs_core.R 的 build_silac_mod 生成 SILAC 修饰串
  # Call build_silac_mod from acs_core.R to assemble the SILAC modification
  if (mode == "SILAC") modification <- build_silac_mod(peptide, modification)
  # （无标记通道=0 口径）：标记作为额外固定修饰追加进修饰串
  #：带观测前体（precursor_mh，Exp.MH 即 M+H）时改走按质量推断的
  # 追加——多行标记同位互斥，旧逐行追加会双计同族质量（K 位 K602+K080）
  # (unlabeled_channel=0): labels append as extra fixed modifications
  #: with an observed precursor the append goes through the
  # mass-inference path - multi-row labels are mutually exclusive per site
  # and the old row loop double-counted same-family masses (K602+K080 on K)
  if (!is.null(labels) && length(labels) > 0) {
    if (!is.na(precursor_mh) && precursor_mh > 0)
      modification <- append_label_mods_by_mass(peptide, modification, labels,
                                                precursor_mh,
                                                aa = AA_mass_dataframe,
                                                aa_tab = amino_acid_table)
    else
      modification <- append_label_mods(peptide, modification, labels)
  }
  # 各标记模式映射到 ReadMGF 的理论离子标签
  # Map each labeling mode to the ReadMGF theoretical-ion label
  label <- switch(mode,
                  unlabeled = "Unlabel",
                  SILAC     = "Unlabel",
                  N15       = "N15",
                  C13       = "C13",
                  C13_comp  = "C13",
                  stop("unknown mode: ", mode))

  # 生成全部理论离子（含内离子表）；结果为空则直接返回 NULL
  # Generate all theoretical ions (including the internal-ion table); NULL if empty
  # （性能第二轮）：走进程级哈希缓存（(peptide|mod|label) 纯函数，
  # 57k 谱候选大量重复；R 写时复制保证后续 match 列修改不污染缓存）
  pos <- cached_all_ions_position(peptide, modification, label, AA_mass_dataframe, amino_acid_table)
  if (length(pos) == 1) return(NULL)

  # 先匹配全部离子（b/y/a/丢失），再匹配内离子，同 Fig.3 代码
  # match all ions (b/y/a/loss) then internal ions, as in the Figure-3 code
  # 对全部离子（b/y/a/丢失）与谱图做匹配
  # Match all ions (b/y/a/loss) against the spectrum
  match_res <- ions_match_selected(pos$ions_table$m_z, spectrum,
                                   mass_shift_tolerence, ppm_or_da)
  # 记录主离子匹配结果，并对剩余谱峰继续匹配内离子
  # Record the main-ion match result, then match internal ions against the remaining peaks
  pos$ions_table$match <- match_res$match_result
  int_match <- ions_match_selected(pos$internal_aa_table$m_z, match_res$mgf_ions,
                                   mass_shift_tolerence, ppm_or_da)
  pos$internal_aa_table$match <- int_match$match_result
  mgf_ions <- int_match$mgf_ions

  # 取 b/y 子表，计算 b/y 强度、覆盖度与三类 PSM 得分
  # Take the b/y sub-table and compute b/y intensities, coverage and the three PSM scores
  b_y <- pos$ions_table[pos$ions_table$Class %in% c("b", "y"),, drop = FALSE]
  b_y_ions_intensity <- b_y$match
  b_y_coverage <- if (length(b_y_ions_intensity) > 0) {
    round(length(which(b_y_ions_intensity > 0)) / length(b_y_ions_intensity) * 100, 2)
  } else 0
  tot_int <- sum(mgf_ions$Intensity)
  PSM_score_by <- if (tot_int > 0) round(sum(b_y$match) / tot_int * 100, 2) else 0
  PSM_score_all <- if (tot_int > 0) round(sum(mgf_ions$Intensity[mgf_ions$Match == 1]) / tot_int * 100, 2) else 0
  PSM_score_all_internal <- if (tot_int > 0) round(sum(mgf_ions$Intensity[mgf_ions$All_Match == 1]) / tot_int * 100, 2) else 0

  # C13 质量补偿：额外匹配窗口偏移 -1/-2 Da（中性质量）
  # C13 mass compensation: extra windows shifted by -1/-2 Da (neutral mass)
  # 初始化质量补偿偏移向量
  # Initialize the mass-compensation shift vector
  shifts <- 0
  if (mode == "C13_comp") shifts <- c(0, -1.003355, -2.006711)
  # 补偿结果初始为 NULL，仅 C13_comp 模式会填充
  # Compensation results start as NULL and are filled only in C13_comp mode
  b_y_ions_intensity_comp <- NULL
  b_y_coverage_comp <- NULL
  # 存在补偿时：对每个偏移独立匹配后累加
  # When compensation is present: match each shift separately and accumulate
  if (length(shifts) > 1) {
    comb <- NULL
    for (s in shifts) {
      tb <- pos$ions_table
      # 偏移量按电荷换算后叠加到理论 m/z 上
      # Convert the shift by charge and add it to the theoretical m/z
      tb$m_z <- tb$m_z + s / tb$Charge
      sm <- ions_match_selected(tb$m_z, spectrum, mass_shift_tolerence, ppm_or_da)
      comb <- if (is.null(comb)) sm$match_result else comb + sm$match_result
    }
    # 汇总补偿匹配结果，并计算补偿后的 b/y 覆盖度
    # Combine the compensated match results and compute compensated b/y coverage
    pos$ions_table$match_comp <- comb
    b_y_c <- pos$ions_table[pos$ions_table$Class %in% c("b", "y"),, drop = FALSE]
    b_y_ions_intensity_comp <- b_y_c$match_comp
    b_y_coverage_comp <- round(length(which(b_y_ions_intensity_comp > 0)) / length(b_y_ions_intensity_comp) * 100, 2)
  }

  # 汇总全部参数并返回
  # Collect all parameters and return
  list(ions_table = pos$ions_table,
       internal_aa_table = pos$internal_aa_table,
       mgf_ions = mgf_ions,
       b_y_ions_intensity = b_y_ions_intensity,
       b_y_ions_intensity_comp = b_y_ions_intensity_comp,
       b_y_coverage = b_y_coverage,
       b_y_coverage_comp = b_y_coverage_comp,
       PSM_score_by = PSM_score_by,
       PSM_score_all = PSM_score_all,
       PSM_score_all_internal = PSM_score_all_internal,
       pep_mass = pos$pep_mass,
       mode = mode)
}
