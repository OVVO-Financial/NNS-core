#!/usr/bin/env Rscript
# tests/fidelity/generate_truth.R
#
# Generates JSON "truth" by calling the installed (live) R NNS package. The sync
# workflow installs the *incoming* R checkout with `R CMD INSTALL` first, so the
# numbers here come from the exact R source being synced. scripts/run_fidelity.py
# then drives the portable C++ core (tests/fidelity/nns_eval) on the same inputs
# and asserts the outputs match within tolerance.
#
# Usage:
#   Rscript generate_truth.R <out.json> [r_commit] [r_version]
#
# Dependency-free JSON emission (no jsonlite required).

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else "truth.json"
r_commit <- if (length(args) >= 2) args[[2]] else "unknown"
r_version <- if (length(args) >= 3) args[[3]] else as.character(utils::packageVersion("NNS"))

suppressMessages(library(NNS))

TOL <- 1e-9

num <- function(v) sprintf("%.17g", as.numeric(v))
arr <- function(v) paste0("[", paste(num(v), collapse = ","), "]")

cases <- list()
add_case <- function(name, func, expected, scalars = list(), vectors = list()) {
  parts <- c(
    sprintf('"name":"%s"', name),
    sprintf('"func":"%s"', func),
    sprintf('"expected":%s', num(expected))
  )
  for (k in names(scalars)) parts <- c(parts, sprintf('"%s":%s', k, num(scalars[[k]])))
  for (k in names(vectors)) parts <- c(parts, sprintf('"%s":%s', k, arr(vectors[[k]])))
  cases[[length(cases) + 1]] <<- paste0("    {", paste(parts, collapse = ","), "}")
}

# ---- fixed, reproducible inputs -------------------------------------------
set.seed(123)
x1 <- c(1, 2, 3, 4, 5)
x2 <- c(-2.5, 0.0, 1.25, 3.75, 8.0, 8.0, 11.0)
xr <- round(rnorm(40), 6)
yr <- round(rnorm(40), 6)

# ---- univariate LPM / UPM --------------------------------------------------
for (deg in c(0, 1, 2)) {
  for (tgt in c(3, mean(x1))) {
    add_case(sprintf("lpm_d%g_t%g", deg, tgt), "lpm",
             LPM(deg, tgt, x1),
             scalars = list(degree = deg, target = tgt),
             vectors = list(x = x1))
    add_case(sprintf("upm_d%g_t%g", deg, tgt), "upm",
             UPM(deg, tgt, x1),
             scalars = list(degree = deg, target = tgt),
             vectors = list(x = x1))
  }
}
for (deg in c(1, 2)) {
  add_case(sprintf("lpm_x2_d%g", deg), "lpm", LPM(deg, mean(x2), x2),
           scalars = list(degree = deg, target = mean(x2)), vectors = list(x = x2))
  add_case(sprintf("upm_x2_d%g", deg), "upm", UPM(deg, mean(x2), x2),
           scalars = list(degree = deg, target = mean(x2)), vectors = list(x = x2))
}

# ---- bivariate Co / D partial moments -------------------------------------
tx <- mean(xr); ty <- mean(yr)
for (deg in c(0, 1)) {
  add_case(sprintf("co_lpm_d%g", deg), "co_lpm", Co.LPM(deg, xr, yr, tx, ty),
           scalars = list(degree = deg, target_x = tx, target_y = ty),
           vectors = list(x = xr, y = yr))
  add_case(sprintf("co_upm_d%g", deg), "co_upm", Co.UPM(deg, xr, yr, tx, ty),
           scalars = list(degree = deg, target_x = tx, target_y = ty),
           vectors = list(x = xr, y = yr))
}
for (dl in c(0, 1)) {
  for (du in c(0, 1)) {
    add_case(sprintf("d_lpm_l%g_u%g", dl, du), "d_lpm", D.LPM(dl, du, xr, yr, tx, ty),
             scalars = list(degree_lpm = dl, degree_upm = du, target_x = tx, target_y = ty),
             vectors = list(x = xr, y = yr))
    add_case(sprintf("d_upm_l%g_u%g", dl, du), "d_upm", D.UPM(dl, du, xr, yr, tx, ty),
             scalars = list(degree_lpm = dl, degree_upm = du, target_x = tx, target_y = ty),
             vectors = list(x = xr, y = yr))
  }
}

body <- paste(unlist(cases), collapse = ",\n")
json <- paste0(
  "{\n",
  sprintf('  "r_version": "%s",\n', r_version),
  sprintf('  "r_commit": "%s",\n', r_commit),
  sprintf('  "tolerance": %s,\n', num(TOL)),
  '  "cases": [\n',
  body, "\n",
  "  ]\n",
  "}\n"
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
writeLines(json, out_path)
cat(sprintf("wrote %d truth cases to %s (NNS %s @ %s)\n",
            length(cases), out_path, r_version, r_commit))
