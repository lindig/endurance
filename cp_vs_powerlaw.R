#!/usr/bin/env Rscript
# cp_vs_powerlaw.R
#
# Power law (PL), Paul's Law (PP) or critical power (CP)?  Compares the one-parameter
# power-law model validated in power_law_analysis.R (t = k * d^1.064) with Paul's Law
# (pace rises by a constant number of seconds per doubling of distance) and with the
# two-parameter critical-power (CP) model, on the SAME athlete-seasons and events.
#
# Design: leave-one-event-out.  For every athlete-season each event is held out in
# turn; every model is fitted to that athlete's remaining events and asked to
# predict the time for the distance of the held-out result.
#
# Models (number of free parameters per athlete in brackets):
#   PL1     power law, exponent fixed at BETA0 = 1.064             [1: level]
#   PL1cal  power law, exponent fixed per sex x age group          [1: level]
#           (exponent from other athletes only, 10-fold cross-fitting)
#   PL1_2k  power law from the athlete's 2 km result alone         [1: level]
#   PL2     power law, free exponent                               [2]
#   PP1     Paul's Law, constant fixed at 5 s per 500 m per doubling [1: level]
#   PPcal   Paul's Law, constant fixed per sex x age group         [1: level]
#           (constant from other athletes only, 10-fold cross-fitting)
#   PP_2k   Paul's Law from the athlete's 2 km result alone        [1: level]
#   PP2     Paul's Law, free constant per athlete                  [2]
#   CP      critical power  P = CP + W'/t,  P = 2.8 (d/t)^3,
#           fitted by the classical regression of work on time     [2]
#   CPlog   same CP model, fitted by least squares on log time     [2]
#           (same loss as the power laws)
#   CPr     CP, classical fit, using only efforts of 2-20 min      [2]
#   CS      critical speed  d = D' + CS * t                        [2]
#
# Input : data/rankings_long.csv  (same file as power_law_analysis.R)
# Output: results/cp_*.csv, results/cp_key_numbers.json, results/fig_cp_*.png
# Run   : Rscript cp_vs_powerlaw.R        (needs power_law_analysis.R alongside)

options(powerlaw.no_main = TRUE)
source("power_law_analysis.R")   # data preparation, eligibility rule, BETA0, ...

## ---------------------------------------------------------------- config ----
CP_CONST   <- 2.8            # Concept2: watts = 2.8 * (metres per second)^3
CP_WINDOW  <- c(120, 1200)   # seconds; usual validity range of the CP model
N_BOOT_CP  <- 1000           # cluster-bootstrap resamples (athletes are clusters)
BETA_FOLDS <- 10             # cross-fitting folds for the calibrated exponent
PP_C500    <- 5              # Paul's Law: seconds per 500 m added per doubling of distance
PP_C       <- PP_C500 / 500  # the same in seconds per metre
PACE_FLOOR <- 0.01           # s/m; a predicted pace below this counts as a (huge) error
N_PACE_BINS <- 8             # pace bins for the "does the constant depend on pace?" test
MODELS     <- c("PL1", "PL1cal", "PL1_2k", "PL2", "PP1", "PPcal", "PP_2k", "PP2",
                "CP", "CPlog", "CPr", "CS")
MAIN       <- c("PL1", "PL1cal", "PL2", "PP1", "PPcal", "PP2", "CP", "CPlog", "CS")
BEST       <- c(MAIN, "CPr")                                    # adds CP on 2-20 min efforts
SINGLE     <- c("PL1", "PL1_2k", "PL1cal", "PP1", "PP_2k", "PPcal")  # incl. 2 km-only models
PAIRS <- data.frame(
  a   = c("PP1", "PPcal", "PP2", "PP2", "CPlog", "PP_2k"),
  b   = c("PL1", "PL1cal", "PL2", "CPlog", "PL2", "PL1_2k"),
  set = c("main", "main", "main", "main", "main", "single"), stringsAsFactors = FALSE)
MODEL_LABELS <- c(PL1 = "PL1 (1.064)", PL1cal = "PL1 (calibrated exponent)",
                  PL1_2k = "PL1 (from 2 km only)", PL2 = "PL2 (free exponent)",
                  PP1 = "PP1 (5 s per 500 m)", PPcal = "PP (calibrated constant)",
                  PP_2k = "PP (from 2 km only)", PP2 = "PP2 (free constant)",
                  CP = "CP (classical fit)", CPlog = "CP (log-time fit)",
                  CPr = "CP (2-20 min efforts)", CS = "Critical speed")
MODEL_COLS <- c(PL1 = "#1f4e9c", PL1cal = "#6fa8dc", PL1_2k = "#999999", PL2 = "#2e8b57",
                PP1 = "#b03a8f", PPcal = "#e08fc0", PP_2k = "#cccccc", PP2 = "#6c2c5a",
                CP = "#c0392b", CPlog = "#e67e22", CPr = "#d4ac0d", CS = "#7f6000")

## ------------------------------------------------------- model fitting -------
# y = s * x + c with s, c >= 0, least squares
nnls2 <- function(x, y) {
  mx <- mean(x); my <- mean(y); sxx <- sum((x - mx)^2)
  if (sxx > 0) {
    s <- sum((x - mx) * (y - my)) / sxx; c0 <- my - s * mx
    if (s >= 0 && c0 >= 0) return(c(slope = s, intercept = c0))
  }
  s1 <- sum(x * y) / sum(x^2); c2 <- my           # boundary candidates: c = 0 or s = 0
  if (sum((y - s1 * x)^2) <= sum((y - c2)^2)) c(slope = s1, intercept = 0)
  else c(slope = 0, intercept = c2)
}

# time for distance d under CP t^3 + W' t^2 = K d^3  (positive root; Newton from an
# upper bound, which converges monotonically because the function is convex)
solve_t <- function(d, cp, wp) {
  rhs <- CP_CONST * d^3
  ta <- if (cp > 0) (rhs / cp)^(1 / 3) else rep(Inf, length(d))
  tb <- if (wp > 0) sqrt(rhs / wp)      else rep(Inf, length(d))
  t <- pmin(ta, tb)
  if (any(!is.finite(t))) return(rep(NA_real_, length(d)))
  for (i in 1:80) {
    f  <- cp * t^3 + wp * t^2 - rhs
    fp <- 3 * cp * t^2 + 2 * wp * t
    tn <- t - f / fp
    if (all(abs(tn - t) <= 1e-12 * t)) { t <- tn; break }
    t <- tn
  }
  t
}

fit_cp <- function(d, t) {                       # classical: work = CP * t + W'
  f <- nnls2(t, CP_CONST * d^3 / t^2)
  c(cp = f[["slope"]], wp = f[["intercept"]])
}
fit_cp_log <- function(d, t, start) {            # same model, least squares on log time
  loss <- function(p) {
    tp <- solve_t(d, exp(p[1]), exp(p[2]))
    if (anyNA(tp)) return(1e10)
    sum((log(tp) - log(t))^2)
  }
  p0 <- unname(log(pmax(start, c(10, 100))))
  o <- optim(p0, loss, method = "Nelder-Mead", control = list(maxit = 300, reltol = 1e-10))
  c(cp = exp(unname(o$par[1])), wp = exp(unname(o$par[2])))
}
fit_cs <- function(d, t) {                       # d = CS * t + D'
  f <- nnls2(t, d)
  c(cs = f[["slope"]], dp = f[["intercept"]])
}
pred_cs <- function(par, dn) {
  if (par[["cs"]] <= 0) return(rep(NA_real_, length(dn)))
  tt <- (dn - par[["dp"]]) / par[["cs"]]
  ifelse(tt > 0, tt, NA_real_)
}
fit_pl2 <- function(d, t) {
  x <- log(d / 2000); y <- log(t)
  b <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
  c(a = mean(y) - b * mean(x), b = b)
}
pred_pl <- function(par, dn) exp(par[["a"]] + par[["b"]] * log(dn / 2000))
fit_pl1 <- function(d, t, beta) c(a = mean(log(t) - beta * log(d / 2000)), b = beta)

## Paul's Law: pace(d) = pace_2km + c * log2(d / 2000)   (pace in s/m, c in s/m per doubling)
pp_level <- function(d, t, cc) {            # pace at 2 km for a fixed c, least squares on log time
  L <- log2(d / 2000); p <- t / d
  lb <- max(1e-4, max(-cc * L) + 1e-6)
  ub <- max(3 * max(p - cc * L), 2 * lb)
  optimize(function(q) sum((log(q + cc * L) - log(p))^2), c(lb, ub), tol = 1e-10)$minimum
}
pp_pred <- function(pref, cc, dn) dn * pmax(pref + cc * log2(dn / 2000), PACE_FLOOR)
fit_pp2 <- function(d, t) {                 # level and constant free, least squares on log time
  L <- log2(d / 2000); p <- t / d; X <- cbind(1, L)
  b <- tryCatch(as.numeric(solve(crossprod(X), crossprod(X, p))), error = function(e) c(mean(p), 0))
  pref0 <- b[1]; crel0 <- b[2] / pref0
  if (!is.finite(pref0) || pref0 <= 0 || any(1 + crel0 * L <= 0)) { pref0 <- mean(p); crel0 <- 0 }
  loss <- function(q) { f <- 1 + q[2] * L; if (any(f <= 0)) return(1e10); sum((q[1] + log(f) - log(p))^2) }
  o <- optim(c(log(pref0), crel0), loss, method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-12))
  c(pref = exp(o$par[1]), crel = o$par[2])
}
pp2_pred <- function(par, dn) dn * pmax(par[["pref"]] * (1 + par[["crel"]] * log2(dn / 2000)), PACE_FLOOR)

# slope per sex x age group estimated on other athletes only (cross-fitting). The slope of a
# group is  sum(num) / sum(den)  over athlete-seasons: exponent (xy/xx) or Paul constant (slp/sll).
make_calib_fun <- function(AS, den, num, fold) {
  AS$fold <- fold[AS$athlete_id]; AS$sg <- paste(AS$sex, AS$age_grp)
  tabs <- lapply(seq_len(BETA_FOLDS), function(f) {
    tr <- AS[AS$fold != f, ]
    lvl <- function(key) {
      r <- rowsum(cbind(den = tr[[den]], num = tr[[num]], one = 1), key)
      list(slope = setNames(r[, "num"] / r[, "den"], rownames(r)), n = setNames(r[, "one"], rownames(r)))
    }
    list(sg = lvl(tr$sg), sex = lvl(tr$sex), all = lvl(rep("all", nrow(tr))))
  })
  function(athlete_id, sex, age_grp) {
    tb <- tabs[[fold[[athlete_id]]]]; k <- paste(sex, age_grp)
    if (!is.na(tb$sg$n[k])  && tb$sg$n[k]  >= MIN_GROUP) return(tb$sg$slope[[k]])
    if (!is.na(tb$sex$n[sex]) && tb$sex$n[sex] >= MIN_GROUP) return(tb$sex$slope[[sex]])
    tb$all$slope[["all"]]
  }
}

## ---------------------------------------------- leave-one-event-out ----------
predict_one <- function(tr, te, beta_cal, c_cal) {
  d <- tr$dist; t <- tr$t; dn <- te$dist; lt <- log(te$t)
  e <- setNames(rep(NA_real_, length(MODELS)), MODELS)
  e["PL1"]    <- log(pred_pl(fit_pl1(d, t, BETA0), dn)) - lt
  e["PL1cal"] <- log(pred_pl(fit_pl1(d, t, beta_cal), dn)) - lt
  i <- which(tr$event == ANCHOR_EVENT)
  if (length(i) == 1) e["PL1_2k"] <- log(t[i] * (dn / d[i])^BETA0) - lt
  pl2 <- fit_pl2(d, t)
  e["PL2"] <- log(pred_pl(pl2, dn)) - lt
  cp <- fit_cp(d, t)
  e["CP"] <- log(solve_t(dn, cp[["cp"]], cp[["wp"]])) - lt
  cpl <- fit_cp_log(d, t, cp)
  e["CPlog"] <- log(solve_t(dn, cpl[["cp"]], cpl[["wp"]])) - lt
  w <- t >= CP_WINDOW[1] & t <= CP_WINDOW[2]
  if (sum(w) >= 2 && length(unique(round(t[w], 6))) >= 2) {
    cpr <- fit_cp(d[w], t[w]); e["CPr"] <- log(solve_t(dn, cpr[["cp"]], cpr[["wp"]])) - lt
  }
  cs <- fit_cs(d, t); e["CS"] <- log(pred_cs(cs, dn)) - lt
  # Paul's Law
  e["PP1"]   <- log(pp_pred(pp_level(d, t, PP_C), PP_C, dn)) - lt
  cc <- max(c_cal, 1e-5)
  e["PPcal"] <- log(pp_pred(pp_level(d, t, cc), cc, dn)) - lt
  if (length(i) == 1) e["PP_2k"] <- log(dn * max(t[i] / d[i] + PP_C * log2(dn / d[i]), PACE_FLOOR)) - lt
  pp2 <- fit_pp2(d, t); e["PP2"] <- log(pp2_pred(pp2, dn)) - lt
  pos <- if (te$t < min(t)) 2 else if (te$t > max(t)) 3 else 1
  c(n_train = length(t), train_ratio = max(t) / min(t), pos = pos, e,
    pl2_b = pl2[["b"]], cp = cp[["cp"]], wp = cp[["wp"]],
    cp_boundary = as.numeric(cp[["cp"]] == 0 || cp[["wp"]] == 0),
    pp2_c500 = 500 * pp2[["pref"]] * pp2[["crel"]])
}

run_loeo <- function(el, beta_fun, c_fun) {
  idx <- split(seq_len(nrow(el)), el$id)
  idx <- idx[vapply(idx, length, 1L) >= 3]
  total <- sum(vapply(idx, length, 1L))
  out <- vector("list", total); rows <- integer(total); k <- 0L
  pb <- 0
  for (g in idx) {
    sub <- el[g, ]
    bc  <- beta_fun(sub$athlete_id[1], sub$sex[1], sub$age_grp[1])
    cc  <- c_fun(sub$athlete_id[1], sub$sex[1], sub$age_grp[1])
    for (j in seq_len(nrow(sub))) {
      k <- k + 1L
      out[[k]] <- predict_one(sub[-j, ], sub[j, ], bc, cc)
      rows[k] <- g[j]
    }
    pb <- pb + 1
    if (pb %% 500 == 0) message(sprintf("   %d / %d athlete-seasons", pb, length(idx)))
  }
  P <- as.data.frame(do.call(rbind, out))
  names(P) <- ifelse(names(P) %in% MODELS, paste0("e_", names(P)), names(P))
  cbind(data.frame(row = rows), P)
}

## ------------------------------------------------------------ scoring -------
pct_err <- function(E) (exp(E) - 1) * 100

score_models <- function(P, models, B = N_BOOT_CP, ref = "PL1") {
  E <- as.matrix(P[, paste0("e_", models), drop = FALSE])
  A <- abs(pct_err(E)); S <- pct_err(E)
  r <- match(ref, models)
  est <- data.frame(model = models, n = colSums(!is.na(A)),
                    mdape = apply(A, 2, median), p90 = apply(A, 2, quantile, 0.9, names = FALSE),
                    bias = apply(S, 2, median),
                    within1 = 100 * colMeans(A <= 1), within2 = 100 * colMeans(A <= 2),
                    row.names = NULL)
  est$win_vs_ref <- vapply(seq_along(models), function(m) {
    w <- sum(A[, m] < A[, r]); l <- sum(A[, m] > A[, r]); if (w + l == 0) NA_real_ else 100 * w / (w + l)
  }, 0)
  idx <- split(seq_len(nrow(P)), P$athlete_idx)
  bm <- bd <- matrix(NA_real_, B, length(models))
  for (b in seq_len(B)) {
    i <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
    m <- apply(A[i, , drop = FALSE], 2, median)
    bm[b, ] <- m; bd[b, ] <- m - m[r]
  }
  q <- function(M, p) apply(M, 2, quantile, p, na.rm = TRUE, names = FALSE)
  est$mdape_lo <- q(bm, .025); est$mdape_hi <- q(bm, .975)
  est$diff <- est$mdape - est$mdape[r]; est$diff_lo <- q(bd, .025); est$diff_hi <- q(bd, .975)
  est
}

# paired comparison of named model pairs (median absolute error of a minus that of b)
pair_diffs <- function(Pmain, Psingle, pairs, B = N_BOOT_CP) {
  do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
    P <- if (pairs$set[k] == "main") Pmain else Psingle
    ea <- abs(pct_err(P[[paste0("e_", pairs$a[k])]])); eb <- abs(pct_err(P[[paste0("e_", pairs$b[k])]]))
    idx <- split(seq_len(nrow(P)), P$athlete_idx)
    bd <- replicate(B, {
      i <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
      median(ea[i]) - median(eb[i])
    })
    w <- sum(ea < eb); l <- sum(ea > eb)
    data.frame(a = pairs$a[k], b = pairs$b[k], set = pairs$set[k], n = nrow(P),
               mdape_a = median(ea), mdape_b = median(eb), diff = median(ea) - median(eb),
               lo = quantile(bd, .025, names = FALSE), hi = quantile(bd, .975, names = FALSE),
               win_a = if (w + l == 0) NA_real_ else 100 * w / (w + l))
  }))
}

breakdown <- function(P, models, dims) {
  E <- as.matrix(P[, paste0("e_", models)]); A <- abs(pct_err(E)); S <- pct_err(E)
  do.call(rbind, lapply(names(dims), function(dn) {
    lv <- dims[[dn]]
    do.call(rbind, lapply(sort(unique(lv[!is.na(lv)])), function(l) {
      s <- which(lv == l)
      data.frame(dimension = dn, level = l, model = models, n = length(s),
                 mdape = apply(A[s, , drop = FALSE], 2, median),
                 p90 = apply(A[s, , drop = FALSE], 2, quantile, 0.9, names = FALSE),
                 bias = apply(S[s, , drop = FALSE], 2, median), row.names = NULL)
    }))
  }))
}

## ---------------------------------------- in-sample fit and parameters --------
insample <- function(el) {
  idx <- split(seq_len(nrow(el)), el$id); idx <- idx[vapply(idx, length, 1L) >= 3]
  res <- lapply(idx, function(g) {
    d <- el$dist[g]; t <- el$t[g]; n <- length(g); lt <- log(t)
    sd_of <- function(tp, p) if (anyNA(tp)) NA_real_ else sqrt(sum((log(tp) - lt)^2) / (n - p))
    p1 <- fit_pl1(d, t, BETA0); p2 <- fit_pl2(d, t)
    cp <- fit_cp(d, t); cpl <- fit_cp_log(d, t, cp); cs <- fit_cs(d, t)
    pr1 <- pp_level(d, t, PP_C); pq <- fit_pp2(d, t)
    c(n = n, sd_PL1 = sd_of(pred_pl(p1, d), 1), sd_PL2 = sd_of(pred_pl(p2, d), 2),
      sd_PP1 = sd_of(pp_pred(pr1, PP_C, d), 1), sd_PP2 = sd_of(pp2_pred(pq, d), 2),
      sd_CP = sd_of(solve_t(d, cp[["cp"]], cp[["wp"]]), 2),
      sd_CPlog = sd_of(solve_t(d, cpl[["cp"]], cpl[["wp"]]), 2),
      sd_CS = sd_of(pred_cs(cs, d), 2),
      pl2_b = p2[["b"]], pp2_c500 = 500 * pq[["pref"]] * pq[["crel"]],
      cp = cpl[["cp"]], wp_kj = cpl[["wp"]] / 1000,
      boundary = as.numeric(cp[["cp"]] == 0 || cp[["wp"]] == 0))
  })
  r <- as.data.frame(do.call(rbind, res)); r$id <- names(res); r
}

# single-effort transfer: source event -> target event with the fixed-constant laws
transfer <- function(el) {
  a <- el[, c("id", "event", "dist", "t")]
  m <- merge(a, a, by = "id", suffixes = c("_s", "_e"))
  m <- m[m$event_s != m$event_e, ]
  m$err_pl <- (m$t_s * (m$dist_e / m$dist_s)^BETA0 / m$t_e - 1) * 100
  m$err_pp <- (m$dist_e * pmax(m$t_s / m$dist_s + PP_C * log2(m$dist_e / m$dist_s), PACE_FLOOR) / m$t_e - 1) * 100
  agg <- do.call(rbind, lapply(split(m, list(m$event_s, m$event_e), drop = TRUE), function(g)
    data.frame(source = g$event_s[1], target = g$event_e[1], n = nrow(g),
               mdape = median(abs(g$err_pl)), bias = median(g$err_pl),
               mdape_pp = median(abs(g$err_pp)), bias_pp = median(g$err_pp))))
  rownames(agg) <- NULL
  agg$source_lbl <- EVENT_LABELS[agg$source]; agg$target_lbl <- EVENT_LABELS[agg$target]
  agg
}

## ----------------------- does Paul's constant depend on pace? -----------------
# Within each athlete-season regress pace (s/m) on log2(distance); pool across athletes.
# Paul's Law: the slope (added seconds per doubling) is the same for everybody.
# Power law: the slope is proportional to the athlete's pace.
paul_stats <- function(el) {
  L <- log2(el$dist / 2000); p <- el$t / el$dist; id <- el$id
  n  <- ag(rep(1, nrow(el)), id)
  Lc <- L - ag(L, id) / n; pc <- p - ag(p, id) / n
  S  <- rowsum(cbind(sll = Lc * Lc, slp = Lc * pc), id)
  info <- el[!duplicated(el$id), c("id", "athlete_id", "sex", "age_grp")]
  A <- cbind(info[match(rownames(S), info$id), ], as.data.frame(S), row.names = NULL)
  lev <- tapply(log(el$t) - BETA0 * log(el$dist / 2000), id, mean)   # level under the fixed law
  A$pace500 <- exp(lev[A$id]) / 2000 * 500                           # 2 km-equivalent pace, s/500 m
  A
}
paul_const <- function(A, label) {
  n_ath <- length(unique(A$athlete_id)); if (n_ath < MIN_GROUP) return(NULL)
  cl <- rowsum(as.matrix(A[, c("sll", "slp")]), A$athlete_id)
  r <- cluster_boot(cl, function(s) c(c500 = 500 * s[["slp"]] / s[["sll"]]))
  data.frame(group = label, n_athletes = n_ath, c500 = r$est, lo = r$lo, hi = r$hi,
             median_pace500 = median(A$pace500))
}
paul_by_pace <- function(A) {
  br <- unique(quantile(A$pace500, seq(0, 1, length.out = N_PACE_BINS + 1), names = FALSE))
  bin <- cut(A$pace500, br, include.lowest = TRUE)
  out <- do.call(rbind, lapply(levels(bin), function(b) {
    r <- paul_const(A[bin == b, ], b); if (is.null(r)) return(NULL)
    r$bin <- b; r
  }))
  out$implied_exponent <- 1 + (out$c500 / out$median_pace500) / log(2)
  out$implied_lo <- 1 + (out$lo / out$median_pace500) / log(2)
  out$implied_hi <- 1 + (out$hi / out$median_pace500) / log(2)
  out
}
# constant against proportional shape, weighted by the bootstrap precision of each bin
paul_shape_test <- function(pb) {
  se <- (pb$hi - pb$lo) / (2 * 1.96); w <- 1 / se^2; x <- pb$median_pace500; y <- pb$c500
  cst <- sum(w * y) / sum(w); k <- sum(w * y * x) / sum(w * x^2)
  list(const_c500 = cst, sse_const = sum(w * (y - cst)^2),
       prop_k = k, beta_prop = 1 + k / log(2), sse_prop = sum(w * (y - k * x)^2), n_bins = nrow(pb))
}

## ------------------------------------------------------------- plots --------
plot_summary <- function(S, file) {
  S <- S[order(S$mdape), ]; yy <- rev(seq_len(nrow(S)))
  png(file.path(OUT_DIR, file), width = 10, height = 4.8, units = "in", res = 200)
  op <- par(mfrow = c(1, 2)); on.exit({ par(op); dev.off() })
  par(mar = c(4, 12, 3, 1))
  plot(S$mdape, yy, xlim = range(c(S$mdape_lo, S$mdape_hi)), yaxt = "n", ylab = "", pch = 19,
       col = MODEL_COLS[S$model], xlab = "median absolute error (%)",
       main = "Accuracy", cex.main = 0.95)
  segments(S$mdape_lo, yy, S$mdape_hi, yy, col = MODEL_COLS[S$model])
  axis(2, at = yy, labels = MODEL_LABELS[S$model], las = 1, cex.axis = 0.8)
  par(mar = c(4, 1, 3, 1))
  plot(S$diff, yy, xlim = range(c(S$diff_lo, S$diff_hi, 0)), yaxt = "n", ylab = "", pch = 19,
       col = MODEL_COLS[S$model], xlab = "difference from PL1 (percentage points)",
       main = "Paired difference vs PL1", cex.main = 0.95)
  segments(S$diff_lo, yy, S$diff_hi, yy, col = MODEL_COLS[S$model]); abline(v = 0, lty = 2)
}

plot_event <- function(B, file) {
  b <- B[B$dimension == "target event", ]
  lv <- unique(b$level); ord <- order(match(lv, ev_order)); lv <- lv[ord]
  png(file.path(OUT_DIR, file), width = 9, height = 5, units = "in", res = 200)
  on.exit(dev.off())
  ymax <- min(max(b$mdape, na.rm = TRUE), quantile(b$mdape, 0.98, na.rm = TRUE) * 1.2)
  plot(NA, xlim = c(1, length(lv)), ylim = c(0, ymax), xaxt = "n", xlab = "held-out event (shortest to longest)",
       ylab = "median absolute error (%)", main = "Prediction error by held-out event")
  axis(1, at = seq_along(lv), labels = lv, cex.axis = 0.8)
  for (m in MAIN) {
    g <- b[b$model == m, ]; y <- g$mdape[match(lv, g$level)]
    lines(seq_along(lv), y, col = MODEL_COLS[[m]], lwd = 2); points(seq_along(lv), y, col = MODEL_COLS[[m]], pch = 19)
  }
  legend("topright", MODEL_LABELS[MAIN], col = MODEL_COLS[MAIN], lwd = 2, pch = 19, bty = "n", cex = 0.7)
}

plot_position <- function(B, file) {
  b <- B[B$dimension == "position", ]; lv <- c("interpolation", "extrapolation to shorter", "extrapolation to longer")
  lv <- lv[lv %in% b$level]
  M <- sapply(MAIN, function(m) sapply(lv, function(l) b$mdape[b$model == m & b$level == l]))
  M <- matrix(M, nrow = length(lv), dimnames = list(lv, MAIN))
  png(file.path(OUT_DIR, file), width = 10, height = 4.8, units = "in", res = 200)
  on.exit(dev.off())
  par(mar = c(4, 4, 3, 1))
  barplot(M, beside = TRUE, col = c("grey85", "grey60", "grey35")[seq_along(lv)], names.arg = MAIN,
          cex.names = 0.75, ylab = "median absolute error (%)", main = "Interpolation versus extrapolation")
  legend("topleft", lv, fill = c("grey85", "grey60", "grey35")[seq_along(lv)], bty = "n", cex = 0.8)
}

heat <- function(M, ev, cols, main, file, fmt = "%.1f", zlim = NULL) {
  png(file.path(OUT_DIR, file), width = 8, height = 6.5, units = "in", res = 200)
  on.exit(dev.off())
  par(mar = c(5, 6, 3, 1))
  if (is.null(zlim)) zlim <- c(0, quantile(M, 0.98, na.rm = TRUE))
  image(seq_along(ev), seq_along(ev), pmin(pmax(t(M), zlim[1]), zlim[2]), col = cols, zlim = zlim,
        axes = FALSE, xlab = "target event (predicted)", ylab = "source event (known)", main = main)
  axis(1, seq_along(ev), ev, las = 2, cex.axis = 0.8); axis(2, seq_along(ev), ev, las = 1, cex.axis = 0.8)
  for (i in seq_along(ev)) for (j in seq_along(ev)) if (!is.na(M[i, j])) text(j, i, sprintf(fmt, M[i, j]), cex = 0.7)
}
plot_transfer <- function(Tm, file_pl, file_diff) {
  ev <- ev_order[ev_order %in% unique(c(Tm$source_lbl, Tm$target_lbl))]
  mk <- function(col) { M <- matrix(NA_real_, length(ev), length(ev), dimnames = list(ev, ev))
    for (i in seq_len(nrow(Tm))) M[Tm$source_lbl[i], Tm$target_lbl[i]] <- Tm[[col]][i]; M }
  Mpl <- mk("mdape"); Mpp <- mk("mdape_pp")
  heat(Mpl, ev, hcl.colors(50, "YlOrRd", rev = TRUE), "PL1 from a single effort: median absolute error (%)", file_pl)
  D <- Mpp - Mpl; z <- max(abs(D), na.rm = TRUE)
  heat(D, ev, hcl.colors(51, "Blue-Red 3"), "Paul's Law minus PL1 (percentage points; blue: Paul better)",
       file_diff, fmt = "%.1f", zlim = c(-z, z))
}

plot_paul <- function(pb, file) {
  png(file.path(OUT_DIR, file), width = 8, height = 5.2, units = "in", res = 200)
  on.exit(dev.off())
  par(mar = c(4.5, 4.5, 3, 1))
  plot(pb$median_pace500, pb$c500, ylim = range(c(pb$lo, pb$hi, 0)), xlim = range(pb$median_pace500) + c(-4, 4),
       pch = 19, xlab = "2 km-equivalent pace (s per 500 m)", ylab = "added seconds per 500 m per doubling of distance",
       main = "Does the constant in Paul's Law depend on pace?")
  segments(pb$median_pace500, pb$lo, pb$median_pace500, pb$hi)
  abline(h = PP_C500, col = "#b03a8f", lwd = 2, lty = 2)
  xs <- seq(par("usr")[1], par("usr")[2], length.out = 50)
  lines(xs, (BETA0 - 1) * log(2) * xs, col = "#1f4e9c", lwd = 2)
  legend("topleft", c("Paul's Law (constant 5 s)", sprintf("power law (exponent %.3f)", BETA0)),
         col = c("#b03a8f", "#1f4e9c"), lwd = 2, lty = c(2, 1), bty = "n", cex = 0.85)
}

plot_examples <- function(el, file) {
  all_idx <- split(seq_len(nrow(el)), el$id)
  idx <- all_idx[vapply(all_idx, length, 1L) >= 5]
  if (length(idx) < 4) idx <- all_idx[vapply(all_idx, length, 1L) >= 4]
  if (length(idx) == 0) return(invisible(NULL))
  pick <- idx[sample(length(idx), min(4, length(idx)))]
  png(file.path(OUT_DIR, file), width = 9, height = 7, units = "in", res = 200)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1)); on.exit({ par(op); dev.off() })
  for (g in pick) {
    d <- el$dist[g]; t <- el$t[g]; P <- CP_CONST * (d / t)^3
    cp <- fit_cp(d, t); cpl <- fit_cp_log(d, t, cp); p2 <- fit_pl2(d, t); p1 <- fit_pl1(d, t, BETA0)
    pr <- pp_level(d, t, PP_C)
    tg <- exp(seq(log(min(t) * 0.7), log(max(t) * 1.4), length.out = 200))
    pl_power <- function(par) { dd <- 2000 * (tg / exp(par[["a"]]))^(1 / par[["b"]]); CP_CONST * (dd / tg)^3 }
    pp_power <- function() {            # Paul: t = d * (pref + c log2(d/2000)); invert numerically
      dgrid <- exp(seq(log(max(min(d) * 0.4, 50)), log(max(d) * 2.5), length.out = 400))
      tt <- pp_pred(pr, PP_C, dgrid); approx(tt, CP_CONST * (dgrid / tt)^3, tg, rule = 2)$y }
    pwp <- pp_power()
    plot(t, P, log = "x", pch = 19, xlab = "time (s, log scale)", ylab = "average power (W)",
         ylim = range(c(P, cpl[["cp"]] + cpl[["wp"]] / tg, pl_power(p2), pl_power(p1), pwp), na.rm = TRUE),
         main = sprintf("%s, %d events", el$sex[g[1]], length(g)))
    lines(tg, pl_power(p1), col = MODEL_COLS[["PL1"]], lwd = 2)
    lines(tg, pwp, col = MODEL_COLS[["PP1"]], lwd = 2)
    lines(tg, pl_power(p2), col = MODEL_COLS[["PL2"]], lwd = 2)
    lines(tg, cpl[["cp"]] + cpl[["wp"]] / tg, col = MODEL_COLS[["CPlog"]], lwd = 2)
  }
  legend("topright", c("PL1 (1.064)", "PP1 (5 s)", "PL2 (free exponent)", "CP (log-time fit)"),
         col = MODEL_COLS[c("PL1", "PP1", "PL2", "CPlog")], lwd = 2, bty = "n", cex = 0.7)
}

## ---------------------------------------------------------------- main ------
main_cp <- function() {
  raw <- read_data(); dd <- prepare(raw)
  el <- eligible_rows(dd, MIN_EVENTS, MIN_RATIO)
  if (nrow(el) == 0) stop("No eligible athlete-seasons.", call. = FALSE)
  st <- athlete_stats(el); AS <- st$AS
  message(sprintf("%d athlete-seasons, %d results", nrow(AS), nrow(el)))

  # sufficient statistics for Paul's constant, and cross-fitted calibration of both laws
  PS <- paul_stats(el)
  AS$sll <- PS$sll[match(AS$id, PS$id)]; AS$slp <- PS$slp[match(AS$id, PS$id)]
  ath  <- unique(AS$athlete_id)
  fold <- setNames(sample(rep_len(seq_len(BETA_FOLDS), length(ath))), ath)
  beta_fun <- make_calib_fun(AS, "xx", "xy", fold)       # exponent (log-log slope)
  c_fun    <- make_calib_fun(AS, "sll", "slp", fold)     # Paul constant (s/m per doubling)

  message("Leave-one-event-out predictions ...")
  P <- run_loeo(el, beta_fun, c_fun)
  P$athlete_idx <- match(el$athlete_id[P$row], unique(el$athlete_id))
  P$event <- el$event[P$row]; P$sex <- el$sex[P$row]; P$age_grp <- el$age_grp[P$row]
  P$season <- el$season[P$row]; P$d <- el$dist[P$row]; P$t <- el$t[P$row]
  P$event_lbl <- EVENT_LABELS[P$event]
  P$position <- c("interpolation", "extrapolation to shorter", "extrapolation to longer")[P$pos]
  P$spread <- as.character(cut(P$train_ratio, c(0, 2, 4, 10, Inf), c("< 2", "2-4", "4-10", ">= 10"), right = FALSE))
  P$n_train_grp <- ifelse(P$n_train >= 4, "4+", as.character(P$n_train))
  med_t <- tapply(P$t, P$event_lbl, median)
  ev_order <<- names(sort(med_t))
  write.csv(data.frame(event = ev_order, median_time_s = as.numeric(sort(med_t))),
            file.path(OUT_DIR, "cp_event_order.csv"), row.names = FALSE)

  # predictions file: no athlete identifiers, only an anonymous index
  write.csv(P[, c("athlete_idx", "season", "sex", "age_grp", "event", "d", "t", "n_train", "train_ratio",
                  "position", paste0("e_", MODELS), "pl2_b", "cp", "wp", "cp_boundary", "pp2_c500")],
            file.path(OUT_DIR, "cp_predictions.csv"), row.names = FALSE)

  # common prediction sets
  Pm <- P[complete.cases(P[, paste0("e_", MAIN)]), ]
  Pb <- P[complete.cases(P[, paste0("e_", BEST)]), ]
  Ps <- P[complete.cases(P[, paste0("e_", SINGLE)]), ]
  message(sprintf("%d predictions; %d common to the main models; %d with CP on 2-20 min efforts; %d with a 2 km effort",
                  nrow(P), nrow(Pm), nrow(Pb), nrow(Ps)))

  message("Scoring ...")
  S_main   <- score_models(Pm, MAIN);   write.csv(S_main,   file.path(OUT_DIR, "cp_summary_main.csv"),   row.names = FALSE)
  S_best   <- score_models(Pb, BEST);   write.csv(S_best,   file.path(OUT_DIR, "cp_summary_best.csv"),   row.names = FALSE)
  S_single <- score_models(Ps, SINGLE); write.csv(S_single, file.path(OUT_DIR, "cp_summary_single.csv"), row.names = FALSE)
  PD <- pair_diffs(Pm, Ps, PAIRS);      write.csv(PD,       file.path(OUT_DIR, "cp_pairs.csv"),          row.names = FALSE)
  dims <- list("target event" = Pm$event_lbl, position = Pm$position, "training spread" = Pm$spread,
               "training events" = Pm$n_train_grp, sex = Pm$sex, "age group" = Pm$age_grp)
  B <- breakdown(Pm, MAIN, dims); write.csv(B, file.path(OUT_DIR, "cp_breakdowns.csv"), row.names = FALSE)

  message("In-sample fit and parameters ...")
  IS <- insample(el)
  sdcols <- c(PL1 = "sd_PL1", PL2 = "sd_PL2", PP1 = "sd_PP1", PP2 = "sd_PP2", CP = "sd_CP", CPlog = "sd_CPlog", CS = "sd_CS")
  ins <- data.frame(model = names(sdcols), params = c(1, 2, 1, 2, 2, 2, 2),
                    median_resid_sd_pct = 100 * vapply(sdcols, function(cn) median(IS[[cn]], na.rm = TRUE), 0))
  two <- IS[, c("sd_PL2", "sd_PP2", "sd_CPlog", "sd_CS")]
  best2 <- apply(two, 1, function(r) if (anyNA(r)) NA_character_ else c("PL2", "PP2", "CPlog", "CS")[which.min(r)])
  share <- prop.table(table(best2)) * 100
  ins$best_two_param_share_pct <- as.numeric(share[ins$model])
  write.csv(ins, file.path(OUT_DIR, "cp_insample.csv"), row.names = FALSE)
  par_tab <- data.frame(
    median_pl2_exponent = median(IS$pl2_b), pl2_exponent_q25 = quantile(IS$pl2_b, .25, names = FALSE),
    pl2_exponent_q75 = quantile(IS$pl2_b, .75, names = FALSE),
    median_pp2_c500 = median(IS$pp2_c500), pp2_c500_q25 = quantile(IS$pp2_c500, .25, names = FALSE),
    pp2_c500_q75 = quantile(IS$pp2_c500, .75, names = FALSE),
    median_cp_w = median(IS$cp), median_wprime_kj = median(IS$wp_kj),
    wprime_q25_kj = quantile(IS$wp_kj, .25, names = FALSE), wprime_q75_kj = quantile(IS$wp_kj, .75, names = FALSE),
    boundary_share_pct = 100 * mean(IS$boundary), n_athlete_seasons = nrow(IS))
  write.csv(par_tab, file.path(OUT_DIR, "cp_params.csv"), row.names = FALSE)

  message("Paul's constant ...")
  pgr <- split_groups(PS); pgr <- pgr[vapply(pgr, function(a) length(unique(a$athlete_id)), 1L) >= MIN_GROUP]
  PG <- do.call(rbind, Map(paul_const, pgr, names(pgr)))
  write.csv(PG, file.path(OUT_DIR, "pp_constant_groups.csv"), row.names = FALSE)
  PB <- paul_by_pace(PS); write.csv(PB, file.path(OUT_DIR, "pp_constant_pace.csv"), row.names = FALSE)
  shape <- paul_shape_test(PB)

  Tm <- transfer(el); write.csv(Tm, file.path(OUT_DIR, "cp_transfer.csv"), row.names = FALSE)

  plot_summary(S_main, "fig_cp_summary.png")
  plot_event(B, "fig_cp_event.png")
  plot_position(B, "fig_cp_position.png")
  plot_transfer(Tm, "fig_cp_transfer.png", "fig_pp_transfer.png")
  plot_paul(PB, "fig_pp_constant.png")
  plot_examples(el, "fig_cp_examples.png")

  key <- list(
    config = list(beta0 = BETA0, pp_c500 = PP_C500, cp_const = CP_CONST, cp_window_min = CP_WINDOW[1] / 60,
                  cp_window_max = CP_WINDOW[2] / 60, n_boot = N_BOOT_CP, beta_folds = BETA_FOLDS,
                  min_events = MIN_EVENTS, min_ratio = MIN_RATIO, pace_bins = N_PACE_BINS,
                  seasons = paste(sort(unique(el$season)), collapse = ", ")),
    data = list(athlete_seasons = nrow(AS), athletes = length(unique(AS$athlete_id)), results = nrow(el),
                predictions = nrow(P), predictions_main = nrow(Pm), predictions_best = nrow(Pb),
                predictions_single = nrow(Ps),
                share_interpolation = mean(Pm$pos == 1), share_extrap_short = mean(Pm$pos == 2),
                share_extrap_long = mean(Pm$pos == 3), cp_boundary_share = mean(P$cp_boundary)),
    paul_shape = shape
  )
  writeLines(to_json(key), file.path(OUT_DIR, "cp_key_numbers.json"))
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "cp_session_info.txt"))

  cat("\n== data ==\n"); cat(sprintf("%d athlete-seasons; %d predictions (%d common to the %d main models; %d with 2 km; %d with CP 2-20 min)\n",
                                     nrow(AS), nrow(P), nrow(Pm), length(MAIN), nrow(Ps), nrow(Pb)))
  cat("\n== accuracy, main models (median abs. error of predicted time, %) ==\n")
  print(S_main[order(S_main$mdape), c("model", "n", "mdape", "mdape_lo", "mdape_hi", "p90", "bias", "within1", "diff", "diff_lo", "diff_hi", "win_vs_ref")],
        digits = 3, row.names = FALSE)
  cat("\n== paired comparisons (a minus b; negative: a more accurate) ==\n")
  print(PD[, c("a", "b", "set", "n", "mdape_a", "mdape_b", "diff", "lo", "hi", "win_a")], digits = 3, row.names = FALSE)
  cat("\n== from the 2 km result alone ==\n")
  print(S_single[order(S_single$mdape), c("model", "n", "mdape", "bias", "diff", "diff_lo", "diff_hi")], digits = 3, row.names = FALSE)
  cat("\n== incl. CP on 2-20 min efforts ==\n")
  print(S_best[order(S_best$mdape), c("model", "n", "mdape", "bias", "diff", "diff_lo", "diff_hi")], digits = 3, row.names = FALSE)
  cat("\n== by position ==\n"); print(B[B$dimension == "position", c("level", "model", "n", "mdape", "bias")], digits = 3, row.names = FALSE)
  cat("\n== Paul's constant (s per 500 m per doubling) by group ==\n"); print(PG, digits = 3, row.names = FALSE)
  cat("\n== Paul's constant by 2 km-equivalent pace ==\n")
  print(PB[, c("bin", "n_athletes", "median_pace500", "c500", "lo", "hi", "implied_exponent")], digits = 3, row.names = FALSE)
  cat(sprintf("weighted fit: constant %.2f s (SSE %.1f) vs proportional exponent %.3f (SSE %.1f)\n",
              shape$const_c500, shape$sse_const, shape$beta_prop, shape$sse_prop))
  cat("\n== in-sample fit ==\n"); print(ins, digits = 3, row.names = FALSE)
  cat("\n== parameters ==\n"); print(par_tab, digits = 3, row.names = FALSE)
  cat(sprintf("\nResults written to %s/\n", OUT_DIR))
  invisible(NULL)
}

if (!interactive()) main_cp()
