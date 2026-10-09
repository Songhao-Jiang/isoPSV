# unit test: MSFragger adapter modification encoding incl. g-slot delta
# run from the tests dir; the adapter probes compat.R next to itself, so switch
# into the code dir (relative path avoids CJK-path corruption) before sourcing
args0 <- commandArgs(FALSE)
sd0 <- dirname(sub("^--file=", "", args0[grepl("--file=", args0)][1]))
code <- file.path("..", "code")
if (!dir.exists(code)) code <- file.path(sd0, "..", "code")
old <- setwd(code)
source("search_adapter_msfragger.R")
setwd(old)

cases <- list(
  list(tok = "20M(15.9949)",            want = "20,Oxidation[M];"),
  list(tok = "3V(14.0157)",             want = "3,Val->Xle;"),
  list(tok = "7P(40.0062)",             want = "7,Pro->His;"),
  list(tok = "2A(999.0)",               want = ""),
  list(tok = "N-term(42.0106)",         want = "0,Acetyl[AnyN-term];"),
  list(tok = "10T(-18.0106)",           want = "10,Glu->pyro-Glu[T];"),
  list(tok = "6G(30.0106)",             want = "6,Gly->Ser;"),
  list(tok = "4M(15.9949), 8Q(0.9840)", want = "4,Oxidation[M];8,Deamidated[Q];")
)
seq8 <- "MGVPTQAM"
npass <- 0L
for (k in seq_along(cases)) {
  got <- msf_mod_string(cases[[k]]$tok, seq8)
  ok <- identical(got, cases[[k]]$want)
  cat(sprintf("case %d %s | tok=%s | got='%s' want='%s'\n",
              k, if (ok) "PASS" else "FAIL",
              cases[[k]]$tok, got, cases[[k]]$want))
  npass <- npass + ok
}
cat(sprintf("msfragger delta test: %d/%d passed\n", npass, length(cases)))
if (npass < length(cases)) quit(status = 1)
