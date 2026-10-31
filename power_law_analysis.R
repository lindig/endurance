#!/usr/bin/env Rscript
# power_law_analysis.R
#
# Tests whether, for one athlete in one season, time follows  t ~ d^BETA0
# (BETA0 = 1.064) across the distances/durations of the Concept2 RowErg rankings.
#
# Input : data/rankings_long.csv   (from concept2_rankings.R)
# Output: results/*.csv, results/key_numbers.json, results/fig_*.png,
#         results/session_info.txt   (all read by the Typst report)
#
# Run:  Rscript power_law_analysis.R
# Base R only; the optional mixed model needs lme4 (skipped if not installed).
#
# Model, per athlete-season:  log t = a + b * log(d / 2000 m)
#   H0: b = BETA0.   Timed events (1, 4, 30, 60 min) enter as (t = fixed time,
#   d = distance rowed); distance events as (d = nominal distance, t = time).

## ---------------------------------------------------------------- config ----
IN_CSV           <- "data/rankings_long.csv"
OUT_DIR          <- "results"
BETA0            <- 1.064
MIN_EVENTS       <- 3      # events needed per athlete-season
MIN_RATIO        <- 3      # longest/shortest distance needed per athlete-season
WINDOW_MIN_RATIO <- 2      # same, for the distance-window analysis
ANCHOR_EVENT     <- "2000" # anchor for the "implied exponent" analysis
N_BOOT           <- 2000   # cluster-bootstrap resamples (athletes are clusters)
N_BOOT_ROBUST    <- 200    # resamples for the (slower) robust-regression variant
MIN_GROUP        <- 30     # smallest stratum/window (athletes) that gets a CI
RUN_MIXED        <- TRUE   # lme4 random-slope model (can take minutes)
SEED             <- 1
TIMED_EVENTS     <- c("1", "4", "30", "60")            # minutes
EVENT_LABELS     <- c("1" = "1 min", "4" = "4 min", "30" = "30 min", "60" = "60 min",
                      "100" = "100 m", "500" = "500 m", "1000" = "1000 m", "2000" = "2000 m",
                      "5000" = "5000 m", "6000" = "6000 m", "10000" = "10,000 m",
                      "21097" = "Half", "42195" = "Marathon", "100000" = "100,000 m")
AGE_BREAKS <- c(0, 18, 29, 39, 49, 59, Inf)
AGE_LABELS <- c("<=18", "19-29", "30-39", "40-49", "50-59", "60+")
WINDOWS <- list(                 # windows on the distance actually rowed (m)
  "all events"    = c(0, Inf),
  "500 m - 10 km" = c(500, 10000),
  "1 km - 10 km"  = c(1000, 10000),
  "2 km - 21 km"  = c(2000, 22000),
  ">= 2 km"       = c(2000, Inf),
  "<= 6 km"       = c(0, 6000)
)

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## --------------------------------------------------------------- helpers ----
FUNNEL <- NULL
step <- function(label, d) {
  ok <- !is.na(d$athlete_id)
  n_as <- length(unique(paste(d$athlete_id[ok], d$season[ok])))
  FUNNEL <<- rbind(FUNNEL, data.frame(step = label, rows = nrow(d), athlete_seasons = n_as))
}

# broadcast of per-group sums back to rows (fast: rowsum + match)
ag <- function(v, id) { s <- rowsum(v, id); s[match(id, rownames(s)), 1] }

# cluster bootstrap: S = matrix of per-cluster sufficient statistics, f(colSums) -> named numbers
cluster_boot <- function(S, f, B = N_BOOT) {
  n <- nrow(S)
  est <- f(colSums(S))
  bs <- replicate(B, f(colSums(S[sample.int(n, n, replace = TRUE), , drop = FALSE])))
  if (is.null(dim(bs))) bs <- matrix(bs, nrow = 1)
  data.frame(stat = names(est) %||% "est", est = unname(est),
             lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE, names = FALSE),
             hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE, names = FALSE),
             se = apply(bs, 1, sd, na.rm = TRUE), row.names = NULL)
}

to_json <- function(x, indent = 0) {
  pad <- strrep("  ", indent)
  if (is.list(x)) {
    if (!is.null(names(x))) {
      items <- vapply(seq_along(x), function(i)
        sprintf('%s  "%s": %s', pad, names(x)[i], to_json(x[[i]], indent + 1)), "")
      paste0("{\n", paste(items, collapse = ",\n"), "\n", pad, "}")
    } else paste0("[", paste(vapply(x, to_json, "", indent + 1), collapse = ", "), "]")
  } else if (length(x) != 1) {
    paste0("[", paste(vapply(x, function(v) to_json(v, indent), ""), collapse = ", "), "]")
  } else if (is.na(x)) "null"
  else if (is.logical(x)) tolower(as.character(x))
  else if (is.numeric(x)) sprintf("%.10g", x)
  else sprintf('"%s"', gsub('"', '\\\\"', as.character(x)))
}

## ----------------------------------------------------------- data prep ------
read_data <- function() {
  hdr <- names(read.csv(IN_CSV, nrows = 1))
  cc <- c(sex = "character", event = "character", athlete_id = "character",
          result_id = "character", type = "character", verified = "character")
  d <- read.csv(IN_CSV, stringsAsFactors = FALSE, colClasses = cc[names(cc) %in% hdr])
  need <- c("season", "event", "sex", "athlete_id", "time_s", "distance_m")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("CSV lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  d
}

prepare <- function(d) {
  step("rows in CSV", d)
  d <- d[!is.na(d$athlete_id) & nzchar(d$athlete_id), ]
  step("rows with an athlete id", d)

  if ("verified" %in% names(d)) {
    d <- d[d$verified %in% "Yes", ]
  } else if (!("filter_status" %in% names(d)) || !all(d$filter_status %in% "verified")) {
    warning("No 'verified' column and no verified-only download filter: unverified rows may be included.")
  }
  step("verified (individual) results", d)
  if ("type" %in% names(d)) {
    d <- d[d$type %in% "R", ]
  } else if (!("filter_rower" %in% names(d)) || !all(d$filter_rower %in% "rower")) {
    warning("No 'type' column and no RowErg-only download filter: dynamic/slides rows may be included.")
  }
  step("RowErg (static) only", d)

  timed <- d$event %in% TIMED_EVENTS
  d$t    <- ifelse(timed, as.numeric(d$event) * 60, d$time_s)
  d$dist <- ifelse(timed, d$distance_m, suppressWarnings(as.numeric(d$event)))
  d <- d[is.finite(d$t) & is.finite(d$dist) & d$t > 0 & d$dist > 0, ]
  step("valid time and distance", d)

  # one result per athlete-season-event: keep the fastest (highest average speed)
  d <- d[order(-(d$dist / d$t)), ]
  d <- d[!duplicated(paste(d$athlete_id, d$season, d$event)), ]
  step("one result per athlete-season-event", d)

  d$id  <- paste(d$athlete_id, d$season, sep = "|")
  d$age <- suppressWarnings(as.numeric(d$age))
  d$age_grp <- as.character(cut(d$age, AGE_BREAKS, AGE_LABELS, right = TRUE))
  d
}

eligible_rows <- function(d, min_events, min_ratio) {
  n  <- ag(rep(1, nrow(d)), d$id)
  mx <- tapply(d$dist, d$id, max)[d$id]
  mn <- tapply(d$dist, d$id, min)[d$id]
  d[n >= min_events & (mx / mn) >= min_ratio, ]
}

# per athlete-season sufficient statistics and slope; also row-level centred values
athlete_stats <- function(d) {
  x <- log(d$dist / 2000); y <- log(d$t); q <- x^2
  n  <- ag(rep(1, nrow(d)), d$id)
  xc <- x - ag(x, d$id) / n
  yc <- y - ag(y, d$id) / n
  qc <- q - ag(q, d$id) / n
  S <- rowsum(cbind(xx = xc * xc, xq = xc * qc, qq = qc * qc,
                    xy = xc * yc, qy = qc * yc, yy = yc * yc, n = 1), d$id)
  info <- d[!duplicated(d$id), c("id", "athlete_id", "season", "sex", "age", "age_grp")]
  AS <- cbind(info[match(rownames(S), info$id), ], as.data.frame(S), row.names = NULL)
  AS$slope <- AS$xy / AS$xx
  rss <- pmax(AS$yy - AS$xy^2 / AS$xx, 0)
  AS$se <- ifelse(AS$n > 2, sqrt(rss / (AS$n - 2) / AS$xx), NA_real_)
  rows <- cbind(d, x = x, y = y, xc = xc, yc = yc, qc = qc)
  list(AS = AS, rows = rows)
}

cluster_table <- function(AS) {
  rowsum(cbind(as.matrix(AS[, c("xx", "xq", "qq", "xy", "qy", "slope")]), one = 1), AS$athlete_id)
}

summarise_fits <- function(AS, label) {
  n_ath <- length(unique(AS$athlete_id))
  out <- data.frame(group = label, n_athlete_seasons = nrow(AS), n_athletes = n_ath,
                    median_slope = median(AS$slope), sd_slope = sd(AS$slope),
                    share_below_beta0 = mean(AS$slope < BETA0),
                    mean_slope = NA_real_, mean_lo = NA_real_, mean_hi = NA_real_,
                    pooled_slope = NA_real_, pooled_lo = NA_real_, pooled_hi = NA_real_)
  if (n_ath >= MIN_GROUP) {
    cl <- cluster_table(AS)
    m <- cluster_boot(cl, function(s) c(mean = s[["slope"]] / s[["one"]]))
    p <- cluster_boot(cl, function(s) c(pooled = s[["xy"]] / s[["xx"]]))
    out[, c("mean_slope", "mean_lo", "mean_hi")]       <- m[, c("est", "lo", "hi")]
    out[, c("pooled_slope", "pooled_lo", "pooled_hi")] <- p[, c("est", "lo", "hi")]
  }
  out$beta0_in_pooled_ci <- with(out, pooled_lo <= BETA0 & BETA0 <= pooled_hi)
  out
}

## ----------------------------------------------------------- analyses -------
fit_mixed <- function(rows) {
  if (!RUN_MIXED || !requireNamespace("lme4", quietly = TRUE)) return(NULL)
  rows$y_adj <- rows$y - BETA0 * rows$x       # so the slope on x estimates (b - BETA0)
  ctl <- lme4::lmerControl(check.conv.singular = "ignore")
  f <- tryCatch(lme4::lmer(y_adj ~ x + (1 + x | id), data = rows, control = ctl), error = function(e) NULL)
  type <- "random intercept + random slope"
  if (is.null(f)) {
    f <- tryCatch(lme4::lmer(y_adj ~ x + (1 | id), data = rows, control = ctl), error = function(e) NULL)
    type <- "random intercept only"
  }
  if (is.null(f)) return(NULL)
  est <- unname(lme4::fixef(f)["x"]); se <- sqrt(diag(as.matrix(vcov(f))))[["x"]]
  vc <- as.data.frame(lme4::VarCorr(f))
  slope_sd <- vc$sdcor[vc$grp == "id" & vc$var1 == "x" & is.na(vc$var2)]
  data.frame(model = type, beta = BETA0 + est, se = se,
             lo = BETA0 + est - 1.96 * se, hi = BETA0 + est + 1.96 * se,
             slope_sd = if (length(slope_sd)) slope_sd else NA_real_,
             singular = lme4::isSingular(f))
}

# implied local exponent between the anchor event and every other event
anchored <- function(d) {
  a <- d[d$event == ANCHOR_EVENT, c("id", "t")]; names(a)[2] <- "t0"
  m <- merge(d[d$event != ANCHOR_EVENT, ], a, by = "id")
  a_dist <- as.numeric(ANCHOR_EVENT)
  m$lr <- log(m$dist / a_dist)
  m <- m[abs(m$lr) > 0.1, ]
  m$implied <- log(m$t / m$t0) / m$lr
  m$pct_dev <- (m$t / (m$t0 * (m$dist / a_dist)^BETA0) - 1) * 100
  res <- lapply(split(m, m$event), function(g) {
    idx <- split(seq_len(nrow(g)), g$athlete_id)
    est <- c(implied = median(g$implied), pct_dev = median(g$pct_dev))
    bs <- replicate(N_BOOT, {
      i <- unlist(idx[sample.int(length(idx), length(idx), replace = TRUE)], use.names = FALSE)
      c(median(g$implied[i]), median(g$pct_dev[i]))
    })
    data.frame(event = g$event[1], n = nrow(g), median_dist = median(g$dist),
               implied_exponent = est[["implied"]],
               implied_lo = quantile(bs[1, ], 0.025, names = FALSE),
               implied_hi = quantile(bs[1, ], 0.975, names = FALSE),
               median_pct_dev = est[["pct_dev"]],
               dev_lo = quantile(bs[2, ], 0.025, names = FALSE),
               dev_hi = quantile(bs[2, ], 0.975, names = FALSE))
  })
  out <- do.call(rbind, res)
  out <- out[order(out$median_dist), ]
  out$label <- EVENT_LABELS[out$event]
  rownames(out) <- NULL
  out
}

## ---------------------------------------------------------------- plots -----
open_png <- function(name, w = 9, h = 6) {
  png(file.path(OUT_DIR, name), width = w, height = h, units = "in", res = 200)
}

forest_plot <- function(df, file, title, xlab = "pooled within-athlete exponent") {
  df <- df[!is.na(df$pooled_slope), ]
  if (nrow(df) == 0) return(invisible(NULL))
  open_png(file, 8, max(3, 0.32 * nrow(df) + 1.5))
  op <- par(mar = c(4.2, 9, 2.5, 1)); on.exit({ par(op); dev.off() })
  yy <- rev(seq_len(nrow(df)))
  plot(df$pooled_slope, yy, xlim = range(c(df$pooled_lo, df$pooled_hi, BETA0), na.rm = TRUE),
       yaxt = "n", ylab = "", xlab = xlab, main = title, pch = 19)
  segments(df$pooled_lo, yy, df$pooled_hi, yy)
  abline(v = BETA0, col = "red", lty = 2)
  axis(2, at = yy, labels = sprintf("%s (n=%d)", df$group, df$n_athletes), las = 1, cex.axis = 0.8)
}


f_quad <- function(s) {
  A <- matrix(c(s[["xx"]], s[["xq"]], s[["xq"]], s[["qq"]]), 2)
  b <- tryCatch(solve(A, c(s[["xy"]], s[["qy"]])), error = function(e) c(NA_real_, NA_real_))
  c(exponent_at_2km = b[1], curvature = b[2])
}

# local exponent at 2 km and curvature (change of exponent per unit log distance, x2)
quad_summary <- function(AS, label) {
  n_ath <- length(unique(AS$athlete_id))
  if (n_ath < MIN_GROUP) return(NULL)
  q <- cluster_boot(cluster_table(AS), f_quad)
  data.frame(group = label, n_athletes = n_ath,
             exp2k = q$est[1], exp2k_lo = q$lo[1], exp2k_hi = q$hi[1],
             curv = q$est[2], curv_lo = q$lo[2], curv_hi = q$hi[2])
}

# all / by sex / by age group / by sex x age group
split_groups <- function(AS) {
  g <- list(all = AS)
  sx <- sort(unique(AS$sex)); ag_ <- AGE_LABELS[AGE_LABELS %in% AS$age_grp]
  for (s in sx) g[[s]] <- AS[AS$sex == s, ]
  for (a in ag_) g[[a]] <- AS[AS$age_grp %in% a, ]
  for (s in sx) for (a in ag_) g[[paste(s, a)]] <- AS[AS$sex == s & AS$age_grp %in% a, ]
  g[vapply(g, nrow, 1L) > 0]
}

age_plot <- function(strata) {
  st <- strata[grepl("^[MF] ", strata$group) & !is.na(strata$pooled_slope), ]
  if (nrow(st) == 0) return(invisible(NULL))
  st$age <- sub("^[MF] ", "", st$group); st$sx <- substr(st$group, 1, 1)
  open_png("fig_age.png", 8, 5)
  plot(NA, xlim = c(0.7, length(AGE_LABELS) + 0.3),
       ylim = range(c(st$pooled_lo, st$pooled_hi, BETA0)),
       xaxt = "n", xlab = "age group", ylab = "pooled within-athlete exponent",
       main = "Exponent by age group and sex")
  axis(1, at = seq_along(AGE_LABELS), labels = AGE_LABELS)
  abline(h = BETA0, col = "red", lty = 2)
  cols <- c(F = "darkorange", M = "steelblue"); off <- c(F = -0.08, M = 0.08)
  for (sx in c("F", "M")) {
    g <- st[st$sx == sx, ]; if (nrow(g) == 0) next
    x <- match(g$age, AGE_LABELS) + off[[sx]]
    segments(x, g$pooled_lo, x, g$pooled_hi, col = cols[[sx]])
    lines(x, g$pooled_slope, col = cols[[sx]]); points(x, g$pooled_slope, pch = 19, col = cols[[sx]])
  }
  legend("topright", c("women", "men"), col = cols[c("F", "M")], pch = 19, bty = "n")
  dev.off()
}


## ------------------------------------------- effort / robustness analysis ----
# Not every logged result is an all-out effort. Sub-maximal results add noise and
# a bias whose sign depends on which events are affected, so we check how far the
# estimates move when such results are down-weighted, trimmed or excluded.

# Huber M-estimator via IRLS (no extra packages)
huber_fit <- function(X, y, k = 1.345, iter = 50) {
  b <- solve(crossprod(X), crossprod(X, y))
  for (i in seq_len(iter)) {
    r <- as.vector(y - X %*% b); sc <- mad(r)
    if (!is.finite(sc) || sc == 0) break
    w <- pmin(1, k * sc / abs(r))
    bn <- solve(crossprod(X, w * X), crossprod(X, w * y))
    done <- max(abs(bn - b)) < 1e-9; b <- bn
    if (done) break
  }
  as.vector(b)
}

robust_fit <- function(M) {   # M columns: xc, qc, yc (centred within athlete-season)
  b2 <- huber_fit(M[, 1:2, drop = FALSE], M[, 3])
  b1 <- huber_fit(M[, 1, drop = FALSE], M[, 3])
  c(exp2k = b2[1], curv = b2[2], linear = b1[1])
}

robust_row <- function(label, rows, B = N_BOOT_ROBUST) {
  M <- cbind(rows$xc, rows$qc, rows$yc)
  idx <- split(seq_len(nrow(M)), rows$athlete_id)
  est <- robust_fit(M)
  bs <- replicate(B, robust_fit(M[unlist(idx[sample.int(length(idx), length(idx), TRUE)], use.names = FALSE), , drop = FALSE]))
  q <- function(i, p) quantile(bs[i, ], p, na.rm = TRUE, names = FALSE)
  data.frame(variant = label, n_athletes = length(idx),
             exp2k = est[["exp2k"]], exp2k_lo = q(1, .025), exp2k_hi = q(1, .975),
             curv = est[["curv"]], curv_lo = q(2, .025), curv_hi = q(2, .975),
             linear = est[["linear"]], linear_lo = q(3, .025), linear_hi = q(3, .975))
}

ols_row <- function(label, AS) {
  cl <- cluster_table(AS)
  q <- cluster_boot(cl, f_quad)
  l <- cluster_boot(cl, function(s) c(linear = s[["xy"]] / s[["xx"]]))
  data.frame(variant = label, n_athletes = length(unique(AS$athlete_id)),
             exp2k = q$est[1], exp2k_lo = q$lo[1], exp2k_hi = q$hi[1],
             curv = q$est[2], curv_lo = q$lo[2], curv_hi = q$hi[2],
             linear = l$est, linear_lo = l$lo, linear_hi = l$hi)
}

skew <- function(v) { v <- v - mean(v); mean(v^3) / mean(v^2)^1.5 }

effort_sensitivity <- function(d, AS, rows, quad) {
  curv <- quad$est[quad$stat == "curvature"]
  b2   <- quad$est[quad$stat == "exponent_at_2km"]
  # residual scatter of each athlete-season after removing the common curvature
  # and the athlete's own level and slope: large values suggest an off effort
  AS$resid_sd <- with(AS, {
    rss <- pmax(yy - 2 * curv * qy + curv^2 * qq - (xy - curv * xq)^2 / xx, 0)
    ifelse(n > 2, sqrt(rss / (n - 2)), NA_real_)
  })
  lim <- function(v, ...) { dw <- d[v, ]; if (!nrow(dw)) return(NULL); eligible_rows(dw, MIN_EVENTS, WINDOW_MIN_RATIO) }
  timed <- d$event %in% TIMED_EVENTS
  out <- list(
    ols_row("baseline (least squares)", AS),
    robust_row("Huber-robust regression", rows),
    ols_row("drop worst-fitting 10% of athletes", AS[AS$resid_sd <= quantile(AS$resid_sd, 0.90, na.rm = TRUE), ]),
    ols_row("drop worst-fitting 25% of athletes", AS[AS$resid_sd <= quantile(AS$resid_sd, 0.75, na.rm = TRUE), ])
  )
  e6 <- lim(d$dist <= 6000)
  if (!is.null(e6) && nrow(e6)) out[[length(out) + 1]] <- ols_row("efforts <= 6 km only", athlete_stats(e6)$AS)
  ed <- lim(!timed)
  if (!is.null(ed) && nrow(ed)) out[[length(out) + 1]] <- ols_row("distance events only (no timed events)", athlete_stats(ed)$AS)
  es <- do.call(rbind, out)

  r <- rows$yc - b2 * rows$xc - curv * rows$qc
  diag <- list(resid_sd_median = median(AS$resid_sd, na.rm = TRUE),
               resid_sd_p90 = quantile(AS$resid_sd, 0.90, na.rm = TRUE, names = FALSE),
               resid_sd_p99 = quantile(AS$resid_sd, 0.99, na.rm = TRUE, names = FALSE),
               residual_skewness = skew(r))
  list(table = es, diag = diag)
}

effort_plot <- function(es) {
  open_png("fig_effort.png", 9, 4.2)
  op <- par(mfrow = c(1, 2)); on.exit({ par(op); dev.off() })
  yy <- rev(seq_len(nrow(es)))
  for (m in c("linear", "exp2k")) {
    est <- es[[m]]; lo <- es[[paste0(m, "_lo")]]; hi <- es[[paste0(m, "_hi")]]
    par(mar = c(4, if (m == "linear") 15 else 1, 3, 1))
    plot(est, yy, xlim = range(c(lo, hi, BETA0), na.rm = TRUE), yaxt = "n", ylab = "", pch = 19,
         xlab = "exponent",
         main = if (m == "linear") "Linear exponent" else "Exponent at 2 km", cex.main = 0.95)
    segments(lo, yy, hi, yy); abline(v = BETA0, col = "red", lty = 2)
    if (m == "linear") axis(2, at = yy, labels = es$variant, las = 1, cex.axis = 0.65)
  }
}

make_plots <- function(AS, rows, anch, strata, windows, overall) {
  # 1. distribution of athlete slopes
  open_png("fig_slopes.png", 8, 5)
  s <- pmin(pmax(AS$slope, 0.7), 1.5)
  hist(s, breaks = seq(0.7, 1.5, by = 0.02), col = "grey85", border = "white",
       main = "Within-athlete exponent (log time vs log distance)",
       xlab = "exponent (values clipped to 0.7-1.5)")
  abline(v = BETA0, col = "red", lwd = 2); abline(v = median(AS$slope), col = "blue", lwd = 2, lty = 2)
  legend("topright", c(sprintf("H0: %.3f", BETA0), sprintf("median: %.3f", median(AS$slope))),
         col = c("red", "blue"), lty = c(1, 2), lwd = 2, bty = "n")
  dev.off()

  # 2. implied local exponent vs the 2 km anchor, by event
  if (nrow(anch) > 0) {
    open_png("fig_implied_exponent.png", 8, 5)
    plot(anch$median_dist, anch$implied_exponent, log = "x", pch = 19,
         ylim = { r <- range(c(anch$implied_lo, anch$implied_hi, BETA0)); r + c(0, 0.1 * diff(r)) },
         xlab = "median distance rowed (m, log scale)",
         ylab = sprintf("implied exponent vs %s anchor", EVENT_LABELS[ANCHOR_EVENT]),
         main = "Where does the power law hold?")
    segments(anch$median_dist, anch$implied_lo, anch$median_dist, anch$implied_hi)
    abline(h = BETA0, col = "red", lty = 2)
    text(anch$median_dist, anch$implied_hi, anch$label, pos = 3, cex = 0.7)
    dev.off()
  }

  # 3. spaghetti of centred log-log curves
  ids <- unique(rows$id); ids <- sample(ids, min(300, length(ids)))
  open_png("fig_loglog.png", 7, 6)
  plot(NA, xlim = range(rows$xc), ylim = quantile(rows$yc, c(0.002, 0.998)),
       xlab = "log distance, centred per athlete-season", ylab = "log time, centred per athlete-season",
       main = "Sample of athlete-seasons vs the hypothesised slope")
  sub <- rows[rows$id %in% ids, ]; sub <- sub[order(sub$id, sub$xc), ]
  for (g in split(sub, sub$id)) lines(g$xc, g$yc, col = rgb(0, 0, 0, 0.12))
  abline(0, BETA0, col = "red", lwd = 2)
  dev.off()

  forest_plot(strata, "fig_strata.png", "Pooled exponent by sex and age group")
  age_plot(strata)
  forest_plot(windows, "fig_windows.png", "Pooled exponent by distance window")
}

## ---------------------------------------------------------------- main ------
main <- function() {
  raw <- read_data()
  d <- prepare(raw)

  # who is in, and are they typical?
  all2k <- d[d$event == ANCHOR_EVENT, ]
  el <- eligible_rows(d, MIN_EVENTS, MIN_RATIO)
  step(sprintf(">= %d events and distance ratio >= %g", MIN_EVENTS, MIN_RATIO), el)
  if (nrow(el) == 0) stop("No eligible athlete-seasons; check the funnel.", call. = FALSE)
  st <- athlete_stats(el); AS <- st$AS; rows <- st$rows

  sel <- rbind(
    data.frame(group = "all verified athletes with a 2 km result", all2k[, c("sex", "t")]),
    data.frame(group = "athletes in the analysis (with a 2 km result)",
               el[el$event == ANCHOR_EVENT, c("sex", "t")]))
  selection <- aggregate(t ~ group + sex, sel, function(v) c(n = length(v), median_s = median(v)))
  selection <- do.call(data.frame, selection)
  names(selection) <- c("group", "sex", "n", "median_2k_time_s")
  write.csv(selection, file.path(OUT_DIR, "selection_check.csv"), row.names = FALSE)
  write.csv(FUNNEL, file.path(OUT_DIR, "funnel.csv"), row.names = FALSE)

  # per-athlete fits
  write.csv(AS[, c("athlete_id", "season", "sex", "age", "n", "slope", "se")],
            file.path(OUT_DIR, "athlete_fits.csv"), row.names = FALSE)

  # overall, by sex, age group, sex x age group
  groups  <- split_groups(AS)
  strata  <- do.call(rbind, Map(summarise_fits, groups, names(groups)))
  overall <- strata[strata$group == "all", ]
  write.csv(strata, file.path(OUT_DIR, "slopes_by_stratum.csv"), row.names = FALSE)

  # curvature: local exponent at 2 km and its change with log distance
  cl <- cluster_table(AS)
  quad <- cluster_boot(cl, f_quad)
  write.csv(quad, file.path(OUT_DIR, "curvature.csv"), row.names = FALSE)
  quad_strata <- do.call(rbind, Map(quad_summary, groups, names(groups)))
  write.csv(quad_strata, file.path(OUT_DIR, "curvature_by_stratum.csv"), row.names = FALSE)

  # which events are actually in the analysis
  ev_all <- aggregate(id ~ event, d,  function(v) length(unique(v)))
  ev_el  <- aggregate(id ~ event, el, function(v) length(unique(v)))
  events_used <- merge(setNames(ev_all, c("event", "athlete_seasons_verified")),
                       setNames(ev_el,  c("event", "athlete_seasons_in_analysis")), all.x = TRUE)
  events_used$label <- EVENT_LABELS[events_used$event]
  events_used <- events_used[order(suppressWarnings(as.numeric(events_used$event))), ]
  write.csv(events_used, file.path(OUT_DIR, "events_used.csv"), row.names = FALSE)

  # distance windows
  windows <- do.call(rbind, lapply(names(WINDOWS), function(w) {
    lim <- WINDOWS[[w]]
    dw <- d[d$dist >= lim[1] & d$dist <= lim[2], ]
    ew <- if (nrow(dw)) eligible_rows(dw, MIN_EVENTS, WINDOW_MIN_RATIO) else dw
    if (nrow(ew) == 0) return(NULL)
    summarise_fits(athlete_stats(ew)$AS, w)
  }))
  write.csv(windows, file.path(OUT_DIR, "slopes_by_window.csv"), row.names = FALSE)

  # is the age/sex pattern just a different event mix? repeat it inside fixed distance windows
  age_window <- do.call(rbind, lapply(c("1 km - 10 km", "2 km - 21 km"), function(w) {
    lim <- WINDOWS[[w]]
    dw <- d[d$dist >= lim[1] & d$dist <= lim[2], ]
    ew <- if (nrow(dw)) eligible_rows(dw, MIN_EVENTS, WINDOW_MIN_RATIO) else dw
    if (nrow(ew) == 0) return(NULL)
    g <- split_groups(athlete_stats(ew)$AS)
    r <- do.call(rbind, Map(summarise_fits, g, names(g))); r$window <- w; r
  }))
  write.csv(age_window, file.path(OUT_DIR, "slopes_by_age_window.csv"), row.names = FALSE)

  # sensitivity to the number of events required
  sens <- do.call(rbind, lapply(c(2, 3, 4, 5), function(k) {
    ek <- eligible_rows(d, k, MIN_RATIO)
    if (nrow(ek) == 0) return(NULL)
    r <- summarise_fits(athlete_stats(ek)$AS, sprintf(">= %d events", k)); r
  }))
  write.csv(sens, file.path(OUT_DIR, "sensitivity_min_events.csv"), row.names = FALSE)

  # anchored analysis (does not need >= 3 events)
  anch <- anchored(d)
  write.csv(anch, file.path(OUT_DIR, "implied_exponent_by_event.csv"), row.names = FALSE)

  eff <- effort_sensitivity(d, AS, rows, quad)
  write.csv(eff$table, file.path(OUT_DIR, "effort_sensitivity.csv"), row.names = FALSE)
  effort_plot(eff$table)

  mixed <- fit_mixed(rows)
  if (!is.null(mixed)) write.csv(mixed, file.path(OUT_DIR, "mixed_model.csv"), row.names = FALSE)

  make_plots(AS, rows, anch, strata[strata$group != "all", ], windows, overall)

  key <- list(
    config = list(beta0 = BETA0, min_events = MIN_EVENTS, min_ratio = MIN_RATIO,
                  anchor_event = ANCHOR_EVENT, n_boot = N_BOOT,
                  seasons = paste(sort(unique(d$season)), collapse = ", ")),
    data = list(rows_in_csv = FUNNEL$rows[1], rows_used = nrow(el),
                athlete_seasons = nrow(AS), athletes = length(unique(AS$athlete_id)),
                events_per_athlete_median = median(AS$n)),
    overall = list(median_slope = overall$median_slope, sd_slope = overall$sd_slope,
                   mean = list(est = overall$mean_slope, lo = overall$mean_lo, hi = overall$mean_hi),
                   pooled = list(est = overall$pooled_slope, lo = overall$pooled_lo, hi = overall$pooled_hi),
                   beta0_in_pooled_ci = overall$beta0_in_pooled_ci,
                   share_below_beta0 = overall$share_below_beta0),
    curvature = list(exponent_at_2km = quad$est[1], exponent_at_2km_lo = quad$lo[1], exponent_at_2km_hi = quad$hi[1],
                     curvature = quad$est[2], curvature_lo = quad$lo[2], curvature_hi = quad$hi[2])
  )
  key$effort <- eff$diag
  if (!is.null(mixed)) key$mixed <- as.list(mixed)
  writeLines(to_json(key), file.path(OUT_DIR, "key_numbers.json"))
  writeLines(c(capture.output(sessionInfo()), "", capture.output(str(key))),
             file.path(OUT_DIR, "session_info.txt"))

  cat("\n== data funnel ==\n"); print(FUNNEL, row.names = FALSE)
  cat("\n== selection check (median 2 km time, s) ==\n"); print(selection, row.names = FALSE)
  cat("\n== overall / strata ==\n")
  print(format(strata[, c("group", "n_athletes", "median_slope", "mean_slope", "pooled_slope", "pooled_lo", "pooled_hi")],
               digits = 4), row.names = FALSE)
  cat("\n== events in the analysis ==\n"); print(events_used, row.names = FALSE)
  cat("\n== curvature (all) ==\n"); print(quad, digits = 4)
  cat("\n== exponent at 2 km and curvature by group ==\n")
  print(format(quad_strata, digits = 4), row.names = FALSE)
  cat("\n== pooled exponent by age group inside fixed windows ==\n")
  print(format(age_window[age_window$group %in% c("all", "M", "F", AGE_LABELS),
                          c("window", "group", "n_athletes", "pooled_slope", "pooled_lo", "pooled_hi")], digits = 4),
        row.names = FALSE)
  cat("\n== windows ==\n")
  print(format(windows[, c("group", "n_athletes", "pooled_slope", "pooled_lo", "pooled_hi")], digits = 4), row.names = FALSE)
  cat("\n== effort sensitivity ==\n")
  print(format(eff$table[, c("variant", "n_athletes", "linear", "linear_lo", "linear_hi", "exp2k", "exp2k_lo", "exp2k_hi", "curv")], digits = 4), row.names = FALSE)
  cat(sprintf("residual sd per athlete-season: median %.4f, p90 %.4f, p99 %.4f; pooled residual skewness %.2f\n",
              eff$diag$resid_sd_median, eff$diag$resid_sd_p90, eff$diag$resid_sd_p99, eff$diag$residual_skewness))
  cat("\n== implied exponent vs anchor ==\n")
  print(format(anch[, c("label", "n", "implied_exponent", "implied_lo", "implied_hi", "median_pct_dev")], digits = 4), row.names = FALSE)
  if (!is.null(mixed)) { cat("\n== mixed model ==\n"); print(mixed, digits = 4) }
  cat(sprintf("\nResults written to %s/\n", OUT_DIR))
  invisible(NULL)
}

if (!interactive()) main()
