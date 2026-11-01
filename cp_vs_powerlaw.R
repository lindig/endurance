#!/usr/bin/env Rscript
# cp_vs_powerlaw.R
#
# One parameter or two?  Compares the one-parameter power-law model validated in
# power_law_analysis.R (t = k * d^1.064, one free parameter k per athlete) with
# two-parameter models, above all the critical-power (CP) model, on the SAME
# athlete-seasons and the SAME events.
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
MODELS     <- c("PL1", "PL1cal", "PL1_2k", "PL2", "CP", "CPlog", "CPr", "CS")
MAIN       <- c("PL1", "PL1cal", "PL2", "CP", "CPlog", "CS")   # compared on common set
BEST       <- c("PL1", "PL1_2k", "PL1cal", "PL2", "CP", "CPlog", "CPr", "CS")
MODEL_LABELS <- c(PL1 = "PL1 (1.064)", PL1cal = "PL1 (calibrated exponent)",
                  PL1_2k = "PL1 (from 2 km only)", PL2 = "PL2 (free exponent)",
                  CP = "CP (classical fit)", CPlog = "CP (log-time fit)",
                  CPr = "CP (2-20 min efforts)", CS = "Critical speed")
MODEL_COLS <- c(PL1 = "#1f4e9c", PL1cal = "#6fa8dc", PL1_2k = "#999999", PL2 = "#2e8b57",
                CP = "#c0392b", CPlog = "#e67e22", CPr = "#8e44ad", CS = "#7f6000")

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

# exponent per sex x age group, estimated on other athletes only (cross-fitting)
make_beta_fun <- function(AS) {
  ath  <- unique(AS$athlete_id)
  fold <- setNames(sample(rep_len(seq_len(BETA_FOLDS), length(ath))), ath)
  AS$fold <- fold[AS$athlete_id]; AS$sg <- paste(AS$sex, AS$age_grp)
  tabs <- lapply(seq_len(BETA_FOLDS), function(f) {
    tr <- AS[AS$fold != f, ]
    lvl <- function(key) {
      r <- rowsum(cbind(xx = tr$xx, xy = tr$xy, one = 1), key)
      list(slope = setNames(r[, "xy"] / r[, "xx"], rownames(r)), n = setNames(r[, "one"], rownames(r)))
    }
    list(sg = lvl(tr$sg), sex = lvl(tr$sex), all = lvl(rep("all", nrow(tr))))
  })
  list(fold = fold, get = function(athlete_id, sex, age_grp) {
    tb <- tabs[[fold[[athlete_id]]]]; k <- paste(sex, age_grp)
    if (!is.na(tb$sg$n[k])  && tb$sg$n[k]  >= MIN_GROUP) return(tb$sg$slope[[k]])
    if (!is.na(tb$sex$n[sex]) && tb$sex$n[sex] >= MIN_GROUP) return(tb$sex$slope[[sex]])
    tb$all$slope[["all"]]
  })
}

## ---------------------------------------------- leave-one-event-out ----------
predict_one <- function(tr, te, beta_cal) {
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
  pos <- if (te$t < min(t)) 2 else if (te$t > max(t)) 3 else 1
  c(n_train = length(t), train_ratio = max(t) / min(t), pos = pos, e,
    pl2_b = pl2[["b"]], cp = cp[["cp"]], wp = cp[["wp"]],
    cp_boundary = as.numeric(cp[["cp"]] == 0 || cp[["wp"]] == 0))
}

run_loeo <- function(el, beta) {
  idx <- split(seq_len(nrow(el)), el$id)
  idx <- idx[vapply(idx, length, 1L) >= 3]
  total <- sum(vapply(idx, length, 1L))
  out <- vector("list", total); rows <- integer(total); k <- 0L
  pb <- 0
  for (g in idx) {
    sub <- el[g, ]
    bc  <- beta$get(sub$athlete_id[1], sub$sex[1], sub$age_grp[1])
    for (j in seq_len(nrow(sub))) {
      k <- k + 1L
      out[[k]] <- predict_one(sub[-j, ], sub[j, ], bc)
      rows[k] <- g[j]
    }
    pb <- pb + 1
    if (pb %% 500 == 0) message(sprintf("   %d / %d athlete-seasons", pb, length(idx)))
  }
  P <- as.data.frame(do.call(rbind, out))
  names(P) <- sub("^(PL1|PL1cal|PL1_2k|PL2|CP|CPlog|CPr|CS)$", "e_\\1", names(P))
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
    c(n = n, sd_PL1 = sd_of(pred_pl(p1, d), 1), sd_PL2 = sd_of(pred_pl(p2, d), 2),
      sd_CP = sd_of(solve_t(d, cp[["cp"]], cp[["wp"]]), 2),
      sd_CPlog = sd_of(solve_t(d, cpl[["cp"]], cpl[["wp"]]), 2),
      sd_CS = sd_of(pred_cs(cs, d), 2),
      pl2_b = p2[["b"]], cp = cpl[["cp"]], wp_kj = cpl[["wp"]] / 1000,
      boundary = as.numeric(cp[["cp"]] == 0 || cp[["wp"]] == 0),
      a_pl1 = p1[["a"]], a_pl2 = p2[["a"]], cp_start = cp[["cp"]], wp_start = cp[["wp"]])
  })
  r <- as.data.frame(do.call(rbind, res)); r$id <- names(res); r
}

# single-effort transfer with the fixed exponent: source event -> target event
transfer <- function(el) {
  a <- el[, c("id", "event", "dist", "t")]
  m <- merge(a, a, by = "id", suffixes = c("_s", "_e"))
  m <- m[m$event_s != m$event_e, ]
  m$err <- (m$t_s * (m$dist_e / m$dist_s)^BETA0 / m$t_e - 1) * 100
  agg <- do.call(rbind, lapply(split(m, list(m$event_s, m$event_e), drop = TRUE), function(g)
    data.frame(source = g$event_s[1], target = g$event_e[1], n = nrow(g),
               mdape = median(abs(g$err)), bias = median(g$err))))
  rownames(agg) <- NULL
  agg$source_lbl <- EVENT_LABELS[agg$source]; agg$target_lbl <- EVENT_LABELS[agg$target]
  agg
}

## ------------------------------------------------------------- plots --------
plot_summary <- function(S, file) {
  S <- S[order(S$mdape), ]; yy <- rev(seq_len(nrow(S)))
  png(file.path(OUT_DIR, file), width = 10, height = 4.2, units = "in", res = 200)
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
  legend("topright", MODEL_LABELS[MAIN], col = MODEL_COLS[MAIN], lwd = 2, pch = 19, bty = "n", cex = 0.8)
}

plot_position <- function(B, file) {
  b <- B[B$dimension == "position", ]; lv <- c("interpolation", "extrapolation to shorter", "extrapolation to longer")
  lv <- lv[lv %in% b$level]
  M <- sapply(MAIN, function(m) sapply(lv, function(l) b$mdape[b$model == m & b$level == l]))
  M <- matrix(M, nrow = length(lv), dimnames = list(lv, MAIN))
  png(file.path(OUT_DIR, file), width = 9, height = 4.5, units = "in", res = 200)
  on.exit(dev.off())
  par(mar = c(4, 4, 3, 1))
  bp <- barplot(M, beside = TRUE, col = c("grey85", "grey60", "grey35")[seq_along(lv)], names.arg = MODEL_LABELS[MAIN],
                cex.names = 0.6, ylab = "median absolute error (%)", main = "Interpolation versus extrapolation")
  legend("topleft", lv, fill = c("grey85", "grey60", "grey35")[seq_along(lv)], bty = "n", cex = 0.8)
}

plot_transfer <- function(Tm, file) {
  ev <- ev_order[ev_order %in% unique(c(Tm$source_lbl, Tm$target_lbl))]
  M <- matrix(NA_real_, length(ev), length(ev), dimnames = list(ev, ev))
  for (i in seq_len(nrow(Tm))) M[Tm$source_lbl[i], Tm$target_lbl[i]] <- Tm$mdape[i]
  png(file.path(OUT_DIR, file), width = 8, height = 6.5, units = "in", res = 200)
  on.exit(dev.off())
  par(mar = c(5, 6, 3, 1))
  zmax <- quantile(M, 0.98, na.rm = TRUE)
  image(seq_along(ev), seq_along(ev), pmin(t(M), zmax), col = hcl.colors(50, "YlOrRd", rev = TRUE),
        axes = FALSE, xlab = "target event (predicted)", ylab = "source event (known)",
        main = "PL1 from a single effort: median absolute error (%)")
  axis(1, seq_along(ev), ev, las = 2, cex.axis = 0.8); axis(2, seq_along(ev), ev, las = 1, cex.axis = 0.8)
  for (i in seq_along(ev)) for (j in seq_along(ev)) if (!is.na(M[i, j])) text(j, i, sprintf("%.1f", M[i, j]), cex = 0.7)
}

plot_examples <- function(el, file) {
  idx <- split(seq_len(nrow(el)), el$id); idx <- idx[vapply(idx, length, 1L) >= 5]
  if (length(idx) < 4) idx <- split(seq_len(nrow(el)), el$id)[vapply(split(seq_len(nrow(el)), el$id), length, 1L) >= 4]
  if (length(idx) == 0) return(invisible(NULL))
  pick <- idx[sample(length(idx), min(4, length(idx)))]
  png(file.path(OUT_DIR, file), width = 9, height = 7, units = "in", res = 200)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1)); on.exit({ par(op); dev.off() })
  for (g in pick) {
    d <- el$dist[g]; t <- el$t[g]; P <- CP_CONST * (d / t)^3
    cp <- fit_cp(d, t); cpl <- fit_cp_log(d, t, cp); p2 <- fit_pl2(d, t); p1 <- fit_pl1(d, t, BETA0)
    tg <- exp(seq(log(min(t) * 0.7), log(max(t) * 1.4), length.out = 200))
    pl_power <- function(par) { dd <- 2000 * (tg / exp(par[["a"]]))^(1 / par[["b"]]); CP_CONST * (dd / tg)^3 }
    plot(t, P, log = "x", pch = 19, xlab = "time (s, log scale)", ylab = "average power (W)",
         ylim = range(c(P, cpl[["cp"]] + cpl[["wp"]] / tg, pl_power(p2), pl_power(p1)), na.rm = TRUE),
         main = sprintf("%s, %d events", el$sex[g[1]], length(g)))
    lines(tg, pl_power(p1), col = MODEL_COLS[["PL1"]], lwd = 2)
    lines(tg, pl_power(p2), col = MODEL_COLS[["PL2"]], lwd = 2)
    lines(tg, cpl[["cp"]] + cpl[["wp"]] / tg, col = MODEL_COLS[["CPlog"]], lwd = 2)
  }
  legend("topright", c("PL1 (1.064)", "PL2 (free exponent)", "CP (log-time fit)"),
         col = MODEL_COLS[c("PL1", "PL2", "CPlog")], lwd = 2, bty = "n", cex = 0.7)
}

## ---------------------------------------------------------------- main ------
main_cp <- function() {
  raw <- read_data(); dd <- prepare(raw)
  el <- eligible_rows(dd, MIN_EVENTS, MIN_RATIO)
  if (nrow(el) == 0) stop("No eligible athlete-seasons.", call. = FALSE)
  st <- athlete_stats(el); AS <- st$AS
  message(sprintf("%d athlete-seasons, %d results", nrow(AS), nrow(el)))

  beta <- make_beta_fun(AS)
  message("Leave-one-event-out predictions ...")
  P <- run_loeo(el, beta)
  P$athlete_idx <- match(el$athlete_id[P$row], unique(el$athlete_id))
  P$event <- el$event[P$row]; P$sex <- el$sex[P$row]; P$age_grp <- el$age_grp[P$row]
  P$season <- el$season[P$row]; P$d <- el$dist[P$row]; P$t <- el$t[P$row]
  P$event_lbl <- EVENT_LABELS[P$event]
  P$position <- c("interpolation", "extrapolation to shorter", "extrapolation to longer")[P$pos]
  P$spread <- cut(P$train_ratio, c(0, 2, 4, 10, Inf), c("< 2", "2-4", "4-10", ">= 10"), right = FALSE)
  P$spread <- as.character(P$spread)
  P$n_train_grp <- ifelse(P$n_train >= 4, "4+", as.character(P$n_train))
  med_t <- tapply(P$t, P$event_lbl, median)
  ev_order <<- names(sort(med_t))
  write.csv(data.frame(event = ev_order, median_time_s = as.numeric(sort(med_t))),
            file.path(OUT_DIR, "cp_event_order.csv"), row.names = FALSE)

  # predictions file: no athlete identifiers, only an anonymous index
  write.csv(P[, c("athlete_idx", "season", "sex", "age_grp", "event", "d", "t", "n_train", "train_ratio",
                  "position", paste0("e_", MODELS), "pl2_b", "cp", "wp", "cp_boundary")],
            file.path(OUT_DIR, "cp_predictions.csv"), row.names = FALSE)

  # common prediction sets
  ok_main <- complete.cases(P[, paste0("e_", MAIN)])
  ok_best <- complete.cases(P[, paste0("e_", BEST)])
  Pm <- P[ok_main, ]; Pb <- P[ok_best, ]
  message(sprintf("%d predictions; %d with all six main models; %d with all eight", nrow(P), nrow(Pm), nrow(Pb)))

  message("Scoring ...")
  S_main <- score_models(Pm, MAIN);  write.csv(S_main, file.path(OUT_DIR, "cp_summary_main.csv"), row.names = FALSE)
  S_best <- score_models(Pb, BEST);  write.csv(S_best, file.path(OUT_DIR, "cp_summary_best.csv"), row.names = FALSE)
  dims <- list("target event" = Pm$event_lbl, position = Pm$position, "training spread" = Pm$spread,
               "training events" = Pm$n_train_grp, sex = Pm$sex, "age group" = Pm$age_grp)
  B <- breakdown(Pm, MAIN, dims); write.csv(B, file.path(OUT_DIR, "cp_breakdowns.csv"), row.names = FALSE)

  message("In-sample fit and parameters ...")
  IS <- insample(el)
  ins <- data.frame(model = c("PL1", "PL2", "CP", "CPlog", "CS"),
                    params = c(1, 2, 2, 2, 2),
                    median_resid_sd_pct = 100 * c(median(IS$sd_PL1, na.rm = TRUE), median(IS$sd_PL2, na.rm = TRUE),
                                                  median(IS$sd_CP, na.rm = TRUE), median(IS$sd_CPlog, na.rm = TRUE),
                                                  median(IS$sd_CS, na.rm = TRUE)))
  two <- IS[, c("sd_PL2", "sd_CPlog", "sd_CS")]
  best2 <- apply(two, 1, function(r) if (anyNA(r)) NA_character_ else c("PL2", "CPlog", "CS")[which.min(r)])
  share <- prop.table(table(best2)) * 100
  ins$best_two_param_share_pct <- as.numeric(share[ins$model]); ins$best_two_param_share_pct[is.na(ins$best_two_param_share_pct)] <- NA
  write.csv(ins, file.path(OUT_DIR, "cp_insample.csv"), row.names = FALSE)
  par_tab <- data.frame(
    median_pl2_exponent = median(IS$pl2_b), pl2_exponent_q25 = quantile(IS$pl2_b, .25, names = FALSE),
    pl2_exponent_q75 = quantile(IS$pl2_b, .75, names = FALSE),
    median_cp_w = median(IS$cp), median_wprime_kj = median(IS$wp_kj),
    wprime_q25_kj = quantile(IS$wp_kj, .25, names = FALSE), wprime_q75_kj = quantile(IS$wp_kj, .75, names = FALSE),
    boundary_share_pct = 100 * mean(IS$boundary), n_athlete_seasons = nrow(IS))
  write.csv(par_tab, file.path(OUT_DIR, "cp_params.csv"), row.names = FALSE)

  Tm <- transfer(el); write.csv(Tm, file.path(OUT_DIR, "cp_transfer.csv"), row.names = FALSE)

  plot_summary(S_main, "fig_cp_summary.png")
  plot_event(B, "fig_cp_event.png")
  plot_position(B, "fig_cp_position.png")
  plot_transfer(Tm, "fig_cp_transfer.png")
  plot_examples(el, "fig_cp_examples.png")

  ref <- S_main[S_main$model == "PL1", ]
  key <- list(
    config = list(beta0 = BETA0, cp_const = CP_CONST, cp_window_min = CP_WINDOW[1] / 60, cp_window_max = CP_WINDOW[2] / 60,
                  n_boot = N_BOOT_CP, beta_folds = BETA_FOLDS, min_events = MIN_EVENTS, min_ratio = MIN_RATIO,
                  seasons = paste(sort(unique(el$season)), collapse = ", ")),
    data = list(athlete_seasons = nrow(AS), athletes = length(unique(AS$athlete_id)), results = nrow(el),
                predictions = nrow(P), predictions_main = nrow(Pm), predictions_best = nrow(Pb),
                share_interpolation = mean(Pm$pos == 1), share_extrap_short = mean(Pm$pos == 2),
                share_extrap_long = mean(Pm$pos == 3), cp_boundary_share = mean(P$cp_boundary))
  )
  writeLines(to_json(key), file.path(OUT_DIR, "cp_key_numbers.json"))
  writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "cp_session_info.txt"))

  fm <- function(x) format(x, digits = 3)
  cat("\n== data ==\n"); cat(sprintf("%d athlete-seasons; %d predictions (%d common to the six main models, %d to all eight)\n",
                                     nrow(AS), nrow(P), nrow(Pm), nrow(Pb)))
  cat("\n== accuracy, six main models (median abs. error of predicted time, %) ==\n")
  print(S_main[order(S_main$mdape), c("model", "n", "mdape", "mdape_lo", "mdape_hi", "p90", "bias", "within1", "within2", "diff", "diff_lo", "diff_hi", "win_vs_ref")],
        digits = 3, row.names = FALSE)
  cat("\n== accuracy incl. best-practice variants ==\n")
  print(S_best[order(S_best$mdape), c("model", "n", "mdape", "mdape_lo", "mdape_hi", "p90", "bias", "diff", "diff_lo", "diff_hi", "win_vs_ref")],
        digits = 3, row.names = FALSE)
  cat("\n== by position ==\n"); print(B[B$dimension == "position", c("level", "model", "n", "mdape", "bias")], digits = 3, row.names = FALSE)
  cat("\n== by target event ==\n"); print(B[B$dimension == "target event", c("level", "model", "n", "mdape", "bias")], digits = 3, row.names = FALSE)
  cat("\n== in-sample fit ==\n"); print(ins, digits = 3, row.names = FALSE)
  cat("\n== parameters ==\n"); print(par_tab, digits = 3, row.names = FALSE)
  cat(sprintf("\nResults written to %s/\n", OUT_DIR))
  invisible(NULL)
}

if (!interactive()) main_cp()
