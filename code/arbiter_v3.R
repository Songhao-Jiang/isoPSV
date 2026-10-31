# arbiter_v3.R - isoPSV  judgment layer (R port of isoPSV_arbiter_v3.py)
#
# Codifies the per-spectrum conflict-adjudication doctrine fixed on
# /29 through the benchmark (ii)/(iii) per-spectrum reviews
# (54 individual + 180 batch rulings); Python reference =
# Agent/isoPSV_arbiter_v3.py (98.0% agreement with the user-ruled final
# classification on benchmark (iii)).
#
# DECISION FLOW (target = strategy-a candidate; competitor = best rival):
#   0. target low quality   unmodified-ladder missing-run >= 3 or bond
#                           coverage < 20%  -> doubtful flag
#   1. no competitor        -> non-chimeric, original correct
#   2. same molecule        competitor ladder reproduces the target ladder
#                           (substitution mods / isobaric terminal adducts)
#                           -> notation equivalent, original correct
#   3. competition gate     shared quantitative fraction > 50% ENTERS the
#                           competition evaluation (entry, not a verdict)
#   4. structural verdict:
#        shared < 50%       independent evidence dominant:
#                             competitor passes quality -> chimeric (both
#                             components confirmed, original correct)
#                             competitor low quality     -> original correct
#        shared >= 50%      containment -> original correct;
#                           competitor low quality -> original correct;
#                           one-sided anchors (target exclusive < 10% while
#                           competitor exclusive >= 20%) -> competition;
#                           few specific ions both sides (<= 5) -> competition;
#                           dual complete ladders, ratio <= 2 -> chimeric;
#                           ratio >= 2 -> competition; default keep original
#   5. target low quality overrides any correct verdict to doubtful.
# Thresholds: gate 50; lowq miss>=3 / cov<20; one-sided 10/20; balance 2.0;
# few-specific 5 peaks.

arb3_gate <- 50.0
arb3_lowq_miss <- 3L
arb3_lowq_cov <- 20.0
arb3_onesided_target <- 10.0
arb3_onesided_comp <- 20.0
arb3_balance_ratio <- 2.0
arb3_few_specific <- 5L

# quality of one candidate AS CLAIMED (mod-aware): (max missing run, bond cov %)
arb3_quality <- function(seq, mod, sp, aa, aa_tab, mz_sorted = NULL, mz_ord = NULL) {
  for (mod_variant in unique(c(mod, ""))) {
    ip <- NULL
    call <- try(suppressWarnings(invisible(capture.output(
      ip <- cached_by_ions_position(seq, mod_variant, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(call, "try-error") || is.null(ip) || length(ip) == 1 ||
        is.null(ip$ions_table)) next
    it <- ip$ions_table
    observed <- as.numeric(sp[[1]])
    if (is.null(mz_sorted)) {
      o <- order(observed, na.last = NA); obs_s <- observed[o]
    } else { o <- mz_ord; obs_s <- mz_sorted }
    hit <- function(m) {
      if (is.na(m)) return(FALSE)
      li <- findInterval(m * (1 - 20e-6), obs_s, left.open = TRUE) + 1L
      ri <- findInterval(m * (1 + 20e-6), obs_s)
      !is.na(li) && !is.na(ri) && li <= ri
    }
    nb <- length(it$b_1_charge)
    bond <- logical(max(nb, 1L))
    for (i in seq_len(nb)) {
      j <- nb + 1L - i  # complementary y index for bond i
      bond[i] <- hit(it$b_1_charge[i]) || hit(it$b_2_charge[i]) ||
        (j >= 1L && j <= length(it$y_1_charge) &&
           (hit(it$y_1_charge[j]) || hit(it$y_2_charge[j])))
    }
    miss <- 0L; cur <- 0L
    for (f in bond) { cur <- if (f) 0L else cur + 1L; if (cur > miss) miss <- cur }
    return(list(miss = miss, cov = 100 * sum(bond) / length(bond)))
  }
  NULL
}

# same-molecule notation: competitor mod-aware ladder equals the target's
# (every b/y entry within 0.01 Da) - covers substitution mods and isobaric
# terminal adducts alike, since equality is asserted on the ladder itself
arb3_same_molecule <- function(seq1, mod1, seq2, mod2, aa, aa_tab) {
  if (identical(seq2, seq1) && identical(mod2, mod1)) return(TRUE)
  if (is.null(mod2) || is.na(mod2) || !nzchar(mod2)) return(FALSE)
  i1 <- i2 <- NULL
  c1 <- try(suppressWarnings(invisible(capture.output(
    i1 <- cached_by_ions_position(seq1, "", "Unlabel", aa, aa_tab)))), silent = TRUE)
  c2 <- try(suppressWarnings(invisible(capture.output(
    i2 <- cached_by_ions_position(seq2, mod2, "Unlabel", aa, aa_tab)))), silent = TRUE)
  if (inherits(c1, "try-error") || inherits(c2, "try-error") ||
      is.null(i1) || is.null(i2) ||
      length(i1) == 1 || length(i2) == 1 ||
      is.null(i1$ions_table) || is.null(i2$ions_table)) return(FALSE)
  t1 <- i1$ions_table; t2 <- i2$ions_table
  if (length(t1$b_1_charge) != length(t2$b_1_charge)) return(FALSE)
  all(abs(c(t1$b_1_charge, t1$b_2_charge, t1$y_1_charge, t1$y_2_charge) -
          c(t2$b_1_charge, t2$b_2_charge, t2$y_1_charge, t2$y_2_charge)) < 0.01)
}

# full adjudication. cand = candidate list (seq/mod/sp), a_idx = index of the
# strategy-a (original) candidate. Returns list(class, rule, keep_original,
# competitor_wins, gate_shared, metrics).
arbitrate_v3 <- function(spec, spec_sorted, spec_ord, cand, a_idx, aa, aa_tab) {
  if (length(cand) < 2 || is.na(a_idx)) {
    return(list(class = "non-chimeric, original correct", rule = "no-competitor",
                keep_original = TRUE, competitor_wins = FALSE,
                gate_shared = NA_real_))
  }
  # best competitor: highest ACS outside the a slot is chosen by the caller;
  # here take the strongest rival by PSM_score_by
  riv <- setdiff(seq_along(cand), a_idx)
  if (length(riv) == 0) {
    return(list(class = "non-chimeric, original correct", rule = "no-competitor",
                keep_original = TRUE, competitor_wins = FALSE,
                gate_shared = NA_real_))
  }
  b_idx <- riv[which.max(vapply(cand[riv], function(x)
    if (is.null(x$sp)) NA_real_ else x$sp$PSM_score_by, numeric(1)))]
  ints <- as.numeric(spec[[2]])
  A <- matched_peak_set(cand[[a_idx]]$seq, cand[[a_idx]]$mod, spec, aa, aa_tab,
                        mz_sorted = spec_sorted, mz_ord = spec_ord)
  B <- matched_peak_set(cand[[b_idx]]$seq, cand[[b_idx]]$mod, spec, aa, aa_tab,
                        mz_sorted = spec_sorted, mz_ord = spec_ord)
  if (length(A) == 0 || length(B) == 0) {
    return(list(class = "non-chimeric, original correct", rule = "no-evidence",
                keep_original = TRUE, competitor_wins = FALSE,
                gate_shared = NA_real_))
  }
  S <- intersect(A, B)
  totA <- sum(ints[A]); totS <- if (length(S)) sum(ints[S]) else 0
  shared_pct <- if (totA > 0) 100 * totS / totA else 0
  ce <- setdiff(B, A); te <- setdiff(A, B)
  base <- if (length(ints)) max(ints) else 1
  ce_i <- if (length(ce)) 100 * sum(ints[ce]) / base else 0
  te_i <- if (length(te)) 100 * sum(ints[te]) / base else 0

  # 0. target low quality (unmodified ladder)
  q1 <- arb3_quality(cand[[a_idx]]$seq, "", spec, aa, aa_tab, spec_sorted, spec_ord)
  target_lowq <- !is.null(q1) && (q1$miss >= arb3_lowq_miss || q1$cov < arb3_lowq_cov)
  # competitor quality as claimed
  q2 <- arb3_quality(cand[[b_idx]]$seq, cand[[b_idx]]$mod, spec, aa, aa_tab,
                     spec_sorted, spec_ord)
  comp_lowq <- is.null(q2) || q2$miss >= arb3_lowq_miss || q2$cov < arb3_lowq_cov
  # 2. same molecule
  if (arb3_same_molecule(cand[[a_idx]]$seq, cand[[a_idx]]$mod,
                         cand[[b_idx]]$seq, cand[[b_idx]]$mod, aa, aa_tab)) {
    cls <- if (target_lowq) "doubtful: target low quality" else "non-chimeric, original correct"
    return(list(class = cls, rule = "same-molecule", keep_original = TRUE,
                competitor_wins = FALSE, gate_shared = round(shared_pct, 1)))
  }
  # 4. structural verdict
  ratio <- if (te_i > 0) ce_i / te_i else 99
  if (shared_pct < arb3_gate) {
    if (comp_lowq) {
      verdict <- "target correct"; rule <- "competitor-low-quality"
    } else {
      verdict <- "chimeric"; rule <- "independent-ladders"
    }
  } else if (length(ce) == 0 || (length(ce) <= 1 && te_i >= arb3_onesided_target)) {
    verdict <- "target correct"; rule <- "containment"
  } else if (comp_lowq) {
    verdict <- "target correct"; rule <- "competitor-low-quality"
  } else if (te_i < arb3_onesided_target && ce_i >= arb3_onesided_comp) {
    verdict <- "competition"; rule <- "one-sided-anchors"
  } else if (max(length(ce), length(te)) <= arb3_few_specific) {
    verdict <- "competition"; rule <- "high-shared-few-specific"
  } else if (length(ce) >= arb3_few_specific && length(te) >= arb3_few_specific &&
             ratio <= arb3_balance_ratio) {
    verdict <- "chimeric"; rule <- "dual-complete-ladders"
  } else if (ratio >= arb3_balance_ratio) {
    verdict <- "competition"; rule <- "dominant-contrast"
  } else {
    verdict <- "target correct"; rule <- "default-keep"
  }
  keep <- verdict != "competition"
  if (target_lowq && verdict %in% c("target correct", "chimeric")) {
    cls <- "doubtful: target low quality"
  } else if (verdict == "target correct") {
    cls <- "non-chimeric, original correct"
  } else if (verdict == "chimeric") {
    cls <- "chimeric, both components confirmed, original correct"
  } else {
    cls <- "non-chimeric, original doubtful (annotated competition)"
  }
  list(class = cls, rule = rule, keep_original = keep,
       competitor_wins = !keep, gate_shared = round(shared_pct, 1),
       winner_seq = cand[[b_idx]]$seq, winner_mod = cand[[b_idx]]$mod,
       ce_n = length(ce), ce_i = round(ce_i, 1), te_n = length(te),
       te_i = round(te_i, 1), target_lowq = target_lowq,
       comp_miss = if (is.null(q2)) NA else q2$miss,
       comp_cov = if (is.null(q2)) NA else round(q2$cov, 1))
}
