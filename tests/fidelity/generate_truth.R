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
# Every case records the R call that produced `expected` (a scalar or a
# vector; matrices are flattened column-major). Cases tagged with "r_fix" cover
# an R bug that NNS-core fixes ahead of R: they are expected to fail against an
# R build without that fix (run_fidelity.py --xfail-r-fixes reports them as
# XFAIL/XPASS instead of failing).
#
# Dependency-free JSON emission (no jsonlite required).

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[[1]] else "truth.json"
r_commit <- if (length(args) >= 2) args[[2]] else "unknown"
r_version <- if (length(args) >= 3) args[[3]] else as.character(utils::packageVersion("NNS"))

suppressMessages(library(NNS))

TOL <- 1e-9

num <- function(v) {
  v <- as.numeric(v)
  s <- sprintf("%.17g", v)
  s[is.nan(v) | is.na(v)] <- "NaN"
  s[!is.na(v) & v == Inf] <- "Infinity"
  s[!is.na(v) & v == -Inf] <- "-Infinity"
  s
}
arr <- function(v) paste0("[", paste(num(v), collapse = ","), "]")
str_json <- function(s) sprintf('"%s"', gsub('"', '\\\\"', s))

cases <- list()
# scalars: named numbers/logicals (logical -> 0/1); vectors: named numeric
# vectors; matrices: named matrices (emitted column-major plus <name>_nrow);
# strings: named character scalars. `expected` may be a scalar or a vector.
add_case <- function(name, func, expected, scalars = list(), vectors = list(),
                     matrices = list(), strings = list(), tol = NULL, r_fix = NULL) {
  parts <- c(
    sprintf('"name":"%s"', name),
    sprintf('"func":"%s"', func),
    sprintf('"expected":%s', arr(expected))
  )
  if (!is.null(tol)) parts <- c(parts, sprintf('"tol":%s', if (length(tol) > 1) arr(tol) else num(tol)))
  if (!is.null(r_fix)) parts <- c(parts, sprintf('"r_fix":%s', str_json(r_fix)))
  for (k in names(scalars)) parts <- c(parts, sprintf('"%s":%s', k, num(as.numeric(scalars[[k]]))))
  for (k in names(vectors)) parts <- c(parts, sprintf('"%s":%s', k, arr(vectors[[k]])))
  for (k in names(matrices)) {
    m <- as.matrix(matrices[[k]])
    parts <- c(parts, sprintf('"%s":%s', k, arr(as.numeric(m))),
               sprintf('"%s_nrow":%d', k, nrow(m)))
  }
  for (k in names(strings)) parts <- c(parts, sprintf('"%s":%s', k, str_json(strings[[k]])))
  cases[[length(cases) + 1]] <<- paste0("    {", paste(parts, collapse = ","), "}")
}

# ---- fixed, reproducible inputs -------------------------------------------
set.seed(123)
x1 <- c(1, 2, 3, 4, 5)
x2 <- c(-2.5, 0.0, 1.25, 3.75, 8.0, 8.0, 11.0)
xr <- round(rnorm(40), 6)
yr <- round(rnorm(40), 6)
ret <- round(rt(600, df = 4) * 0.01, 8)                      # heavy-tailed returns
price <- round(400 * cumprod(1 + rnorm(1500, 4e-4, 0.012)), 4) # price levels (offset)
disc <- as.numeric(sample(-3:3, 300, TRUE))                   # discrete with ties

# ---- univariate LPM / UPM (scalar direct path) -----------------------------
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
for (deg in c(0.5, 1.5, 3)) {
  add_case(sprintf("lpm_ret_d%g", deg), "lpm", LPM(deg, 0.001, ret),
           scalars = list(degree = deg, target = 0.001), vectors = list(x = ret))
  add_case(sprintf("upm_ret_d%g", deg), "upm", UPM(deg, 0.001, ret),
           scalars = list(degree = deg, target = 0.001), vectors = list(x = ret))
}

# ---- vectorised LPM / UPM (> 32 targets: prefix backend) -------------------
# Targets: a grid across (and beyond) the support plus order statistics.
grid_targets <- function(v) c(seq(min(v) - 0.1 * diff(range(v)), max(v) + 0.1 * diff(range(v)),
                                  length.out = 60), sort(v)[seq(1, length(v), by = 15)])
for (data_name in c("ret", "price", "disc")) {
  v <- get(data_name)
  tv <- grid_targets(v)
  for (deg in c(0, 1, 2, 3, 4, 8, 12, 0.5, 2.5)) {
    # Degrees >= 3 on these data hit the NNS 13.0-13.3 binomial-expansion
    # cancellation (CORE-01); fractional degrees use the direct path.
    fix <- if (deg >= 3 && deg == round(deg)) "CORE-01" else NULL
    # Degrees 1-2 agree to rounding with the old backend except in far tails.
    tol <- if (deg %in% c(1, 2)) 1e-7 else NULL
    add_case(sprintf("lpm_v_%s_d%g", data_name, deg), "lpm_v", LPM(deg, tv, v),
             scalars = list(degree = deg), vectors = list(target = tv, x = v), tol = tol, r_fix = fix)
    add_case(sprintf("upm_v_%s_d%g", data_name, deg), "upm_v", UPM(deg, tv, v),
             scalars = list(degree = deg), vectors = list(target = tv, x = v), tol = tol, r_fix = fix)
  }
  for (deg in c(0, 1, 2, 4)) {
    fix <- if (deg >= 3) "CORE-01" else NULL
    tol <- if (deg %in% c(1, 2)) 1e-7 else NULL
    add_case(sprintf("lpm_ratio_v_%s_d%g", data_name, deg), "lpm_ratio_v", LPM.ratio(deg, tv, v),
             scalars = list(degree = deg), vectors = list(target = tv, x = v), tol = tol, r_fix = fix)
    add_case(sprintf("upm_ratio_v_%s_d%g", data_name, deg), "upm_ratio_v", UPM.ratio(deg, tv, v),
             scalars = list(degree = deg), vectors = list(target = tv, x = v), tol = tol, r_fix = fix)
  }
  # <= 32 targets stay on the direct path at every degree.
  t20 <- tv[seq(1, length(tv), length.out = 20)]
  for (deg in c(3, 8)) {
    add_case(sprintf("upm_v20_%s_d%g", data_name, deg), "upm_v", UPM(deg, t20, v),
             scalars = list(degree = deg), vectors = list(target = t20, x = v))
  }
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
# unequal, fractional and higher degrees; ties
xd <- disc[1:120]; yd <- disc[121:240]
for (dg in list(c(2, 0.5), c(0, 3), c(1.5, 1.5), c(3, 2))) {
  add_case(sprintf("co_lpm_%g_%g", dg[1], dg[2]), "co_lpm", Co.LPM(dg[1], xr, yr, tx, ty, dg[2]),
           scalars = list(degree = dg[1], degree_y = dg[2], target_x = tx, target_y = ty),
           vectors = list(x = xr, y = yr))
  add_case(sprintf("co_upm_ties_%g_%g", dg[1], dg[2]), "co_upm", Co.UPM(dg[1], xd, yd, 0, 0, dg[2]),
           scalars = list(degree = dg[1], degree_y = dg[2], target_x = 0, target_y = 0),
           vectors = list(x = xd, y = yd))
  add_case(sprintf("d_lpm_%g_%g", dg[1], dg[2]), "d_lpm", D.LPM(dg[1], dg[2], xd, yd, 0, 0),
           scalars = list(degree_lpm = dg[1], degree_upm = dg[2], target_x = 0, target_y = 0),
           vectors = list(x = xd, y = yd))
  add_case(sprintf("d_upm_%g_%g", dg[1], dg[2]), "d_upm", D.UPM(dg[1], dg[2], xr, yr, tx, ty),
           scalars = list(degree_lpm = dg[1], degree_upm = dg[2], target_x = tx, target_y = ty),
           vectors = list(x = xr, y = yr))
}
# vectorised bivariate (target recycling: 40 x targets, 1 y target)
txv <- sort(xr); tyv <- ty
for (dg in list(c(0, 0), c(1, 2), c(0.5, 1))) {
  add_case(sprintf("co_lpm_v_%g_%g", dg[1], dg[2]), "co_lpm_v",
           Co.LPM(dg[1], xr, yr, txv, tyv, dg[2]),
           scalars = list(degree_1 = dg[1], degree_2 = dg[2]),
           vectors = list(x = xr, y = yr, target_x = txv, target_y = tyv))
  add_case(sprintf("co_upm_v_%g_%g", dg[1], dg[2]), "co_upm_v",
           Co.UPM(dg[1], xr, yr, txv, tyv, dg[2]),
           scalars = list(degree_1 = dg[1], degree_2 = dg[2]),
           vectors = list(x = xr, y = yr, target_x = txv, target_y = tyv))
  add_case(sprintf("d_lpm_v_%g_%g", dg[1], dg[2]), "d_lpm_v",
           D.LPM(dg[1], dg[2], xr, yr, txv, tyv),
           scalars = list(degree_1 = dg[1], degree_2 = dg[2]),
           vectors = list(x = xr, y = yr, target_x = txv, target_y = tyv))
  add_case(sprintf("d_upm_v_%g_%g", dg[1], dg[2]), "d_upm_v",
           D.UPM(dg[1], dg[2], xr, yr, txv, tyv),
           scalars = list(degree_1 = dg[1], degree_2 = dg[2]),
           vectors = list(x = xr, y = yr, target_x = txv, target_y = tyv))
}

# ---- n-dimensional co-partial moments, batch, PM.matrix --------------------
X3 <- cbind(ret[1:200], ret[201:400], ret[401:600])
X3t <- cbind(disc[1:100], disc[101:200], disc[201:300])
for (nm in c("X3", "X3t")) {
  X <- get(nm); tgt <- colMeans(X)
  for (deg in c(0, 0.5, 1, 2)) {
    for (nrm in c(TRUE, FALSE)) {
      sc <- list(degree = deg, norm = nrm)
      add_case(sprintf("clpm_nd_%s_d%g_n%d", nm, deg, nrm), "clpm_nd", Co.LPM_nD(X, tgt, deg, nrm),
               scalars = sc, vectors = list(target = tgt), matrices = list(data = X))
      add_case(sprintf("cupm_nd_%s_d%g_n%d", nm, deg, nrm), "cupm_nd", Co.UPM_nD(X, tgt, deg, nrm),
               scalars = sc, vectors = list(target = tgt), matrices = list(data = X))
      add_case(sprintf("dpm_nd_%s_d%g_n%d", nm, deg, nrm), "dpm_nd", DPM_nD(X, tgt, deg, nrm),
               scalars = sc, vectors = list(target = tgt), matrices = list(data = X))
    }
  }
  targets <- X[seq(1, nrow(X), by = 7), , drop = FALSE]
  for (deg in c(0, 1)) {
    add_case(sprintf("clpm_nd_batch_%s_d%g", nm, deg), "clpm_nd_batch",
             NNS:::Co.LPM_nD.batch(X, targets, deg, TRUE),
             scalars = list(degree = deg, norm = TRUE), matrices = list(data = X, targets = targets))
  }
  for (cfg in list(c(1, 1, 1, 0), c(0, 0, 0, 0), c(2, 1, 1, 1), c(0.5, 1.5, 0, 1))) {
    pm <- PM.matrix(cfg[1], cfg[2], tgt, X, pop_adj = as.logical(cfg[3]), norm = as.logical(cfg[4]))
    add_case(sprintf("pm_matrix_%s_%g_%g_%g_%g", nm, cfg[1], cfg[2], cfg[3], cfg[4]), "pm_matrix",
             c(pm$cupm, pm$dupm, pm$dlpm, pm$clpm, pm$cov.matrix),
             scalars = list(degree_lpm = cfg[1], degree_upm = cfg[2], pop_adj = cfg[3], norm = cfg[4]),
             vectors = list(target = tgt), matrices = list(variable = X))
  }
}

# Degree-0 counts with targets ON data values (ties): CLPM + CUPM + DPM must
# partition the sample. R 13.3 counted ties as both lower and discordant/upper;
# fixed in the CORE-09 R patch, so these cases are tagged r_fix.
Xtie <- cbind(disc[1:100], disc[101:200], disc[201:300])
ttie <- apply(Xtie, 2, median)
for (nrm in c(TRUE, FALSE)) {
  sc <- list(degree = 0, norm = nrm)
  add_case(sprintf("clpm_nd_tie_n%d", nrm), "clpm_nd", Co.LPM_nD(Xtie, ttie, 0, nrm),
           scalars = sc, vectors = list(target = ttie), matrices = list(data = Xtie))
  add_case(sprintf("cupm_nd_tie_n%d", nrm), "cupm_nd", Co.UPM_nD(Xtie, ttie, 0, nrm),
           scalars = sc, vectors = list(target = ttie), matrices = list(data = Xtie), r_fix = "CORE-09")
  add_case(sprintf("dpm_nd_tie_n%d", nrm), "dpm_nd", DPM_nD(Xtie, ttie, 0, nrm),
           scalars = sc, vectors = list(target = ttie), matrices = list(data = Xtie), r_fix = "CORE-09")
}

# ---- gravity / mode --------------------------------------------------------
cvec <- list(
  normal = round(rnorm(500), 6), ret = ret, disc = disc, small4 = c(0, 0, 1, 1),
  small5 = c(1, 1, 1, 2, 9), bimodal = round(c(rnorm(300, -2), rnorm(300, 2)), 5),
  cauchy = round(rcauchy(400), 5), tight = c(5 + 1e-5 * round(rnorm(300), 6), 6, -30)
)
for (nm in names(cvec)) {
  v <- cvec[[nm]]
  for (dsc in c(FALSE, TRUE)) {
    add_case(sprintf("gravity_%s_d%d", nm, dsc), "gravity", NNS.gravity(v, discrete = dsc),
             scalars = list(discrete = dsc), vectors = list(x = v))
    for (mlt in c(FALSE, TRUE)) {
      add_case(sprintf("mode_%s_d%d_m%d", nm, dsc, mlt), "mode", NNS.mode(v, discrete = dsc, multi = mlt),
               scalars = list(discrete = dsc, multi = mlt), vectors = list(x = v))
    }
  }
}
big <- c(3e9, 3e9, 3e9, 1, 2)
add_case("mode_discrete_3e9", "mode", NNS.mode(big, discrete = TRUE, multi = FALSE),
         scalars = list(discrete = TRUE, multi = FALSE), vectors = list(x = big), r_fix = "CORE-07")
# Equally high peaks: the single mode averages them instead of taking one by sort order.
eqp <- rep(0:19, each = 400)
add_case("mode_equal_peaks", "mode", NNS.mode(eqp, discrete = FALSE, multi = FALSE),
         scalars = list(discrete = FALSE, multi = FALSE), vectors = list(x = eqp), r_fix = "MODE-TIES")

# ---- partition (NNS_part_cpp) ----------------------------------------------
px <- ret[1:300]; py <- ret[301:600]
for (noise in c("off", "mean", "median", "mode")) {
  for (xo in c(FALSE, TRUE)) {
    p <- NNS:::NNS_part_cpp(px, py, if (xo) "XONLY" else NULL, NULL, 8L, FALSE, noise)
    q <- as.character(p$dt$quadrant)
    qc <- ifelse(nchar(q) > 1, as.numeric(substring(q, 2)), 0)
    add_case(sprintf("partition_%s_x%d", noise, xo), "partition",
             c(p$order, qc, p$regression.points$x, p$regression.points$y),
             scalars = list(xonly = xo, obs_req = 8), vectors = list(x = px, y = py),
             strings = list(noise_reduction = noise))
  }
}
ix <- c(1:20, Inf, -Inf); iy <- c(20:1, 5, 7)
p <- NNS:::NNS_part_cpp(ix, iy, NULL, 2L, 8L, FALSE, "off")
q <- as.character(p$dt$quadrant)
add_case("partition_inf", "partition",
         c(p$order, ifelse(nchar(q) > 1, as.numeric(substring(q, 2)), 0),
           p$regression.points$x, p$regression.points$y),
         scalars = list(order = 2, obs_req = 8), vectors = list(x = ix, y = iy),
         strings = list(noise_reduction = "off"))

# ---- dependence (NNS_dep_matrix_cpp) ---------------------------------------
# NNS_dep_matrix_cpp is not called by any R function and was removed from R as
# dead code in a later patch; only emit these cases when the R build has it.
# The random draws happen either way so the cases after this block keep their
# inputs.
Xd <- cbind(ret[1:300], ret[301:600], 0.5 * ret[1:300] + ret[301:600])
Xi <- matrix(round(rnorm(900), 6), ncol = 3)
if (exists("NNS_dep_matrix_cpp", envir = asNamespace("NNS"), inherits = FALSE)) {
  dm <- NNS:::NNS_dep_matrix_cpp(Xd)
  add_case("dep_matrix_ret", "dep_matrix", c(dm$Correlation, dm$Dependence), matrices = list(X = Xd))
  dm <- NNS:::NNS_dep_matrix_cpp(Xi)
  add_case("dep_matrix_indep", "dep_matrix", c(dm$Correlation, dm$Dependence), matrices = list(X = Xi))
}

# ---- dependence (NNS_dep_pair_cpp, which NNS.dep calls) ---------------------
# Quadrant labels go to R as strings and to the core as integer cell ids.
dep_quads <- function(x, y) {
  l <- length(x)
  q <- NNS:::NNS_part_cpp(x, y, "XONLY", as.integer(max(ceiling(log(l, 2)), 1L)),
                          as.integer(max(8L, as.integer(l / 8L))), FALSE, "off", TRUE)$quadrant
  match(q, unique(q))
}
add_dep_pair <- function(name, x, y, qa, qb, r_fix = NULL) {
  for (asym in c(FALSE, TRUE)) {
    d <- NNS:::NNS_dep_pair_cpp(as.numeric(x), as.numeric(y), as.character(qa), as.character(qb), asym)
    add_case(sprintf("dep_pair_%s_a%d", name, asym), "dep_pair", c(d$Correlation, d$Dependence),
             scalars = list(asym = asym),
             vectors = list(x = x, y = y, quad_xy = qa, quad_yx = qb), r_fix = r_fix)
  }
}
add_dep_pair("ret", ret[1:300], ret[301:600], dep_quads(ret[1:300], ret[301:600]),
             dep_quads(ret[301:600], ret[1:300]))
# Every point discordant: R 13.3 reported Correlation +1 (DEP-01).
add_dep_pair("negative", ret[1:300], -ret[1:300], dep_quads(ret[1:300], -ret[1:300]),
             dep_quads(-ret[1:300], ret[1:300]), r_fix = "DEP-01")
# 49 * (1/49) != 1, so R 13.3 gave dep(b, b) = 0.51 (DEP-08).
bdep <- c(rep(0, 51), rep(1, 49))
add_dep_pair("binary_self", bdep, bdep, bdep + 1, bdep + 1, r_fix = "DEP-08")
# Points on the cell mean were counted as co-lower and discordant (CORE-09).
xtie <- rep(c(-1, 0, 1), 33); ytie <- rep(c(0, 1, -1, 0, -1, 1, 1, 0, -1), 11)
add_dep_pair("ties", xtie, ytie, rep(1, 99), rep(1, 99), r_fix = "CORE-09")

# ---- distance (NNS_distance_path_cpp / bulk) --------------------------------
RPM <- matrix(round(runif(120 * 3), 6), ncol = 3); yh <- round(rnorm(120), 6)
Xt <- matrix(round(runif(15 * 3), 6), ncol = 3)
add_case("distance_path", "distance_path", NNS:::NNS_distance_path_cpp(RPM, yh, Xt, 8L, FALSE),
         scalars = list(k = 8, is_class = FALSE), vectors = list(yhat = yh),
         matrices = list(RPM = RPM, Xtest = Xt))
add_case("distance_bulk", "distance_bulk", NNS:::NNS_distance_bulk_cpp(RPM, yh, Xt, 5L, FALSE),
         scalars = list(k = 5, is_class = FALSE), vectors = list(yhat = yh),
         matrices = list(RPM = RPM, Xtest = Xt))

# ---- fast_lm -------------------------------------------------------------------
i500 <- 1:500
fx <- 1e9 + i500 * 1e-3; fy <- ((i500 * 7919) %% 1000) / 1000
z <- NNS:::fast_lm(fx, fy)
# A fitted value cancels terms of size |b0| + |b1 x|, so its rounding (FMA or not) scales with that.
fit_tol <- 4 * .Machine$double.eps * (abs(z$coef[1]) + abs(z$coef[2] * fx))
add_case("fast_lm_offset", "fast_lm", c(z$coef, z$fitted.values), vectors = list(x = fx, y = fy),
         tol = c(1e-12, 1e-12, fit_tol))
z <- NNS:::fast_lm(xr, yr)
add_case("fast_lm_small", "fast_lm", c(z$coef, z$fitted.values), vectors = list(x = xr, y = yr))

# ---- M.reg / stack kernels (R 13.1-13.3) -----------------------------------
Xm <- cbind(a = round(runif(400), 3), b = round(rnorm(400), 3), c = round(rnorm(400), 2))
Xm <- cbind(Xm, a2 = Xm[, "a"])                     # duplicated column
ym <- round(2 * Xm[, 1] - Xm[, 2] + rnorm(400, 0, 0.3), 6)
add_case("duplicate_column_map", "duplicate_column_map",
         NNS:::NNS_duplicate_column_map_cpp(Xm), matrices = list(X = Xm))
bounds <- lapply(seq_len(ncol(Xm)), function(j) unique(quantile(Xm[, j], seq(0, 1, 0.125), names = FALSE)))
bvec <- unlist(bounds); blen <- vapply(bounds, length, 1)
for (red in 0:3) {
  for (cls in c(FALSE, TRUE)) {
    yy <- if (cls) round(ym) else ym
    s <- NNS:::NNS_mreg_setup_cpp(Xm, yy, bounds, red, cls)
    add_case(sprintf("mreg_setup_r%d_c%d", red, cls), "mreg_setup", s$RPM,
             scalars = list(reducer = red, is_class = cls),
             vectors = list(y = yy, boundaries = bvec, boundary_lengths = blen),
             matrices = list(X = Xm),
             # reducer 2 ("mode"): R 13.3's native setup took the most frequent
             # exact value; fixed to NNS.mode in the CORE-03 R patch.
             r_fix = if (red == 2L) "CORE-03" else NULL)
  }
}
s <- NNS:::NNS_mreg_setup_cpp(Xm[, 1:3], ym, bounds[1:3], 3L, FALSE)
rx <- s$RPM[, 1:3]; ry <- s$RPM[, 4]
mins <- apply(Xm[, 1:3], 2, min); maxs <- apply(Xm[, 1:3], 2, max)
Xq <- cbind(round(runif(40), 3), round(rnorm(40), 3), round(rnorm(40), 2))
Xq[1:5, ] <- rx[1:5, ]                               # exact hits (distance 0 ties)
for (dc in 0:3) {
  for (k in c(1, 3, 7)) {
    for (cls in c(FALSE, TRUE)) {
      yv <- if (cls) round(ry) else ry
      Xuse <- if (dc == 3) round(Xq) else Xq
      Ruse <- if (dc == 3) round(rx) else rx
      add_case(sprintf("mreg_predict_d%d_k%d_c%d", dc, k, cls), "mreg_predict",
               NNS:::NNS_mreg_predict_v2_cpp(Ruse, yv, Xuse, as.integer(k), dc, mins, maxs, cls, 1L),
               scalars = list(k = k, dist_code = dc, is_class = cls),
               vectors = list(yhat = yv, mins = mins, maxs = maxs),
               matrices = list(rpm_x = Ruse, Xtest = Xuse))
    }
  }
  add_case(sprintf("mreg_predict_path_d%d", dc), "mreg_predict_path",
           NNS:::NNS_mreg_predict_path_v2_cpp(rx, ry, Xq, 12L, dc, mins, maxs, FALSE, 1L),
           scalars = list(k = 12, dist_code = dc, is_class = FALSE),
           vectors = list(yhat = ry, mins = mins, maxs = maxs),
           matrices = list(rpm_x = rx, Xtest = Xq))
}
# 1-D fast path: duplicated coordinate (cbind(X*, X*)), with unequal ranges.
x1d <- sort(round(runif(150), 4)); r1d <- cbind(x1d, x1d); y1d <- round(sin(6 * x1d), 6)
q1d <- cbind(round(runif(30), 4)); q1d <- cbind(q1d, q1d)
for (dc in 0:2) {
  add_case(sprintf("mreg_predict_path_1d_d%d", dc), "mreg_predict_path",
           NNS:::NNS_mreg_predict_path_v2_cpp(r1d, y1d, q1d, 9L, dc, c(0, 0), c(1, 2), FALSE, 1L),
           scalars = list(k = 9, dist_code = dc, is_class = FALSE),
           vectors = list(yhat = y1d, mins = c(0, 0), maxs = c(1, 2)),
           matrices = list(rpm_x = r1d, Xtest = q1d))
  add_case(sprintf("mreg_predict_1d_d%d_k4", dc), "mreg_predict",
           NNS:::NNS_mreg_predict_v2_cpp(r1d, y1d, q1d, 4L, dc, c(0, 0), c(1, 2), FALSE, 1L),
           scalars = list(k = 4, dist_code = dc, is_class = FALSE),
           vectors = list(yhat = y1d, mins = c(0, 0), maxs = c(1, 2)),
           matrices = list(rpm_x = r1d, Xtest = q1d))
}
# Range scaling makes (.07, .07) equidistant from .01 and .13; the raw differences are not.
xt1 <- cbind(c(.01, .13, 1), c(.01, .13, 1))
add_case("mreg_predict_1d_tie", "mreg_predict",
         NNS:::NNS_mreg_predict_v2_cpp(xt1, c(0, 100, 999), matrix(.07, 1, 2), 1L, 0L,
                                       c(.01, .01), c(1, 1), FALSE, 1L),
         scalars = list(k = 1, dist_code = 0, is_class = FALSE),
         vectors = list(yhat = c(0, 100, 999), mins = c(.01, .01), maxs = c(1, 1)),
         matrices = list(rpm_x = xt1, Xtest = matrix(.07, 1, 2)), r_fix = "MREG-1D")
# Dividing by the range puts .03 strictly nearer (.26, .26) than .49; 1/range ties them.
xd <- cbind(c(.03, .49, 1), c(.03, .49, 1)); qd <- rbind(c(.26, .26), c(.2, .3))
add_case("mreg_predict_range_divide", "mreg_predict",
         NNS:::NNS_mreg_predict_v2_cpp(xd, c(0, 100, 999), qd, 1L, 0L, c(.03, .03), c(1, 1), FALSE, 1L),
         scalars = list(k = 1, dist_code = 0, is_class = FALSE),
         vectors = list(yhat = c(0, 100, 999), mins = c(.03, .03), maxs = c(1, 1)),
         matrices = list(rpm_x = xd, Xtest = qd), r_fix = "MREG-DIV")
# A query about 58 from every regression point: the Student-t weights underflow.
xf <- matrix(0, 1000, 60); yf <- c(numeric(999), 100); qf <- matrix(0, 1, 60)
for (i in 0:999) {
  xf[i + 1, (i + 0:28) %% 60 + 1] <- 1
  xf[i + 1, (i + 29) %% 60 + 1] <- (sqrt(1 + 4 * (59.57 - i / 999 - 58)) - 1) / 2
}
add_case("mreg_predict_far_underflow", "mreg_predict",
         NNS:::NNS_mreg_predict_v2_cpp(xf, yf, qf, 1000L, 0L, rep(0, 60), rep(1, 60), FALSE, 1L),
         scalars = list(k = 1000, dist_code = 0, is_class = FALSE),
         vectors = list(yhat = yf, mins = rep(0, 60), maxs = rep(1, 60)),
         matrices = list(rpm_x = xf, Xtest = qf))

# X* projection path
trd <- matrix(round(rnorm(60 * 4), 5), ncol = 4); ted <- matrix(round(rnorm(10 * 4), 5), ncol = 4)
for (cf in list(c(0.5, -1.2, 0, 2), c(0, 0, 0, 0), c(1, 1, 1, 1))) {
  ord <- order(abs(cf), decreasing = TRUE, method = "radix")
  xs <- NNS:::NNS_xstar_path_cpp(trd, ted, cf, as.integer(ord), 1L)
  add_case(sprintf("xstar_path_%s", paste(cf, collapse = "_")), "xstar_path",
           c(xs$train, xs$test, xs$denominator),
           vectors = list(coef = cf, column_order = ord), matrices = list(train = trd, test = ted))
}

# ---- stochastic dominance (NNS_FSD/SSD/TSD_uni_cpp) --------------------------
# x is never longer than y: R 13.3 walked y with x's length, which reads past
# the end of a shorter y.
add_sd <- function(name, degree, x, y, type = "discrete", r_fix = NULL) {
  v <- switch(degree, NNS.FSD.uni(x, y, type), NNS.SSD.uni(x, y), NNS.TSD.uni(x, y))
  add_case(sprintf("sd_uni_%s", name), "sd_uni", v,
           scalars = list(degree = degree, discrete = type == "discrete"),
           vectors = list(x = x, y = y), r_fix = r_fix)
}
sa <- round(rnorm(50), 4); sb <- round(sa - abs(rnorm(50)), 4)
for (deg in 1:3) {
  add_sd(sprintf("shift_d%d", deg), deg, sa, sb)
  add_sd(sprintf("shift_rev_d%d", deg), deg, sb, sa)
}
add_sd("shift_cont", 1, sa, sb, "continuous")
su <- round(runif(30, 1, 3), 3); sv <- round(runif(45, 0, 2), 3)
add_sd("unequal_d1", 1, su, sv, r_fix = "SDA-08")
add_sd("unequal_short_d1", 1, 4, c(3, 0), r_fix = "SDA-08")
add_sd("unequal_short_d2", 2, c(1, 4), c(3, 1, 3, 0, 2), r_fix = "SDA-08")
add_sd("unequal_short_d3", 3, c(4, 1), c(0, 3, 3, 4), r_fix = "SDA-08")
# Equal curves, and ratios equal at every knot, are not dominance.
add_sd("coincide_cont", 1, c(1, 2), c(1, 1, 2, 2), "continuous", r_fix = "SDA-07")
add_sd("coincide_d2", 2, c(1, 2), c(1, 1, 2, 2), r_fix = "SDA-07")
add_sd("coincide_d3", 3, c(1, 2), c(1, 1, 2, 2), r_fix = "SDA-07")
add_sd("knots_cross_xy", 1, c(0, 0, 3, 3), c(1, 1, 3, 0), "continuous", r_fix = "SDA-07")
add_sd("knots_cross_yx", 1, c(1, 1, 3, 0), c(0, 0, 3, 3), "continuous", r_fix = "SDA-07")
add_sd("midpoint_xy", 1, c(0, 1, 1), c(0, 0, 1), "continuous", r_fix = "SDA-07")
add_sd("midpoint_yx", 1, c(0, 0, 1), c(0, 1, 1), "continuous", r_fix = "SDA-07")
add_sd("midpoint_mass_xy", 1, c(1, 1), c(0, 1), "continuous", r_fix = "SDA-07")
# Continuous ratio 0/0 at a point mass is the CDF value.
add_sd("ratio_00_mass", 1, rep(4, 6), c(5, 5, 0, 4, 0, 2), "continuous", r_fix = "SDA-13")
# The curves cross between sample points.
add_sd("interior_d3", 3, c(-0.85, -1.52, 0.68), c(-0.21, -0.33, -1.67), r_fix = "SDA-05")
add_sd("interior_cont", 1, c(0.15, 1.55, -0.86, -0.8, 0.25, 0.33),
       c(-0.15, 1.17, 0.26, -0.04, -2.44, 1.19), "continuous", r_fix = "SDA-05")
# Rounding-level differences are ties; the large sample needs compensated
# prefix sums to stay within the tolerance.
add_sd("tol_sum", 2, c(0.1, 0.7), c(0, 0.8), r_fix = "SDA-09")
zq <- qnorm(ppoints(100))
add_sd("tol_spread", 2, zq + 0.37, 2 * zq + 0.37, r_fix = "SDA-09")
zr <- round(qnorm(ppoints(1e4)), 1)
add_sd("tol_spread_1e4_d2", 2, 0.1 * zr + 0.3, 0.2 * zr + 0.3, r_fix = "SDA-09")
add_sd("tol_spread_1e4_d3", 3, 0.1 * zr + 0.3, 0.2 * zr + 0.3, r_fix = "SDA-09")

# ---- seasonality (NNS_seas_cpp) ---------------------------------------------
# R 13.3 scored three reverse-stride subsamples per lag; R now scores all m
# component series, weighted by length, with < and a near-zero-mean sd rule.
add_seas <- function(name, x, modulo = NULL, mod_only = TRUE, r_fix = "TS-04") {
  s <- NNS:::NNS_seas_cpp(as.numeric(x), if (is.null(modulo)) NULL else as.integer(modulo), mod_only)
  ap <- s$all.periods
  vec <- list(x = x)
  if (!is.null(modulo)) vec$modulo <- modulo
  add_case(sprintf("seas_%s", name), "seasonality",
           c(s$best.period, ap$Period, ap$Coefficient.of.Variation,
             ap$Variable.Coefficient.of.Variation),
           scalars = list(mod_only = mod_only), vectors = vec, r_fix = r_fix)
}
ap_x <- as.numeric(AirPassengers)
add_seas("airpassengers", ap_x)
add_seas("airpassengers_mod12", ap_x, 12)
add_seas("airpassengers_mod12_all", ap_x, 12, mod_only = FALSE)
add_seas("lynx", as.numeric(lynx))
add_seas("toy", c(10, 18, 11, 21, 12, 20, 13, 23))
add_seas("odd", ap_x[1:143])
sz <- sin(2 * pi * (1:96) / 12)
add_seas("zero_mean", sz - mean(sz))
add_seas("constant", rep(5, 20))
add_seas("short", c(1, 2, 3, 4), r_fix = NULL)

# ---- ARMA seasonal weights (ARMA.seas.weighting) -----------------------------
# NNS.ARMA falls back to equal weights in its R wrapper when
# (VarCV - CV) / sqrt(lag) sums to 0; R 13.3 has no wrapper and returns zeros.
add_arma_w <- function(name, M, r_fix = NULL) {
  ns <- asNamespace("NNS")
  fn <- if (exists(".arma_lag_weights", envir = ns, inherits = FALSE)) ns$.arma_lag_weights
        else ns$ARMA.seas.weighting
  w <- fn(FALSE, M)
  add_case(sprintf("arma_seas_weighting_%s", name), "arma_seas_weighting", c(w$lag, w$Weights),
           vectors = list(periods = M$Period, covar = M$Coefficient.of.Variation,
                          varcovar = M$Variable.Coefficient.of.Variation), r_fix = r_fix)
}
add_arma_w("airpassengers", NNS:::NNS_seas_cpp(ap_x, NULL, TRUE)$all.periods)
add_arma_w("sentinel", data.frame(Period = c(7, 14, 21), Coefficient.of.Variation = 0.3,
                                  Variable.Coefficient.of.Variation = 0.3), r_fix = "TS-11")
add_arma_w("constant", data.frame(Period = c(12, 24), Coefficient.of.Variation = 0,
                                  Variable.Coefficient.of.Variation = 0), r_fix = "TS-11")

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
