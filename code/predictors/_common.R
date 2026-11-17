# _common.R - shared layer for the three predictor adapters (v2.0.179).
# parse_cli / the fingerprint-sidecar trio (input_fp, sidecar_matches,
# write_run_params) / the pipe-notation peptide guard were verbatim copies in
# pred_pdeep3.R, pred_alphapeptdeep.R and pred_ms2pip.R; they live here now.
# This file is ASCII-ONLY by design: the bilingual guard strings are written
# as \uXXXX escapes, restored at R parse time independently of file encoding
# (same pattern as the step-2 background-script generator).
# Adapters source it via the --file= bootstrap at their head.
#
# Sidecar contract (unchanged): <out>.run_params holds in_md5=/in_size= lines;
# the master's file.copy refreshes peptides.txt mtime every run, so mtime
# ordering never holds - reuse is decided by CONTENT fingerprint only, and
# legacy artifacts without a sidecar rerun once.

parse_cli <- function(args = commandArgs(TRUE)) {
  # parsed_args: flag name -> value/TRUE; token: current argv piece;
  # flag_key: flag name without --; flag_value_pair: the --k=v pair halves
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

# Pipe-notation isoform peptides (pFind either-or notation, e.g. PEPT|DEK)
# are skipped whole with a count: the drivers cannot express the either-or
# form and a first-variant pick would silently invent a sequence (v2.0.169;
# counts must stay visible in the run log for the master's reconciliation
# line). Returns the filtered data.frame.
guard_pipe_peptides <- function(pep_df, tag) {
  hit <- grepl("|", pep_df$peptide, fixed = TRUE)
  if (any(hit)) {
    cat(sprintf(paste0("[%s] \u7ba1\u9053\u7b26\u5f02\u6784\u4f53\u80bd\uff08pFind ",
                       "\u4e8c\u9009\u4e00\u8bb0\u6cd5\uff09\u6574\u6761\u8df3\u8fc7 %d \u6761\n"),
                tag, sum(hit)))
    cat(sprintf("[%s] skipped %d pipe-notation isoform peptides (pFind either-or notation)\n", tag, sum(hit)))
    pep_df <- pep_df[!hit, , drop = FALSE]
  }
  pep_df
}
