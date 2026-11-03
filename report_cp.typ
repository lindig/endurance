// report_cp.typ — power law, Paul's Law or critical power?
// Compares the one-parameter power law (t = k d^1.064), Paul's Law (a fixed number of
// seconds per doubling of distance) and two-parameter models (critical power, critical
// speed, free-exponent power law, free-constant Paul's Law).
// Reads everything from results/ (written by power_law_analysis.R and cp_vs_powerlaw.R).
// Compile with:  typst compile report_cp.typ
// Numbers, tables and figures come from the results files; the interpretive sentences
// are conditional on those numbers but should still be re-read.

#set page(paper: "a4", margin: (x: 2.2cm, y: 2.3cm), numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#show heading: set block(above: 1.3em, below: 0.7em)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)

// ------------------------------------------------------------------ data ----
#let key    = json("results/cp_key_numbers.json")
#let sm     = csv("results/cp_summary_main.csv", row-type: dictionary)
#let sb     = csv("results/cp_summary_best.csv", row-type: dictionary)
#let ss     = csv("results/cp_summary_single.csv", row-type: dictionary)
#let pairs  = csv("results/cp_pairs.csv", row-type: dictionary)
#let bd     = csv("results/cp_breakdowns.csv", row-type: dictionary)
#let ins    = csv("results/cp_insample.csv", row-type: dictionary)
#let par    = csv("results/cp_params.csv", row-type: dictionary).first()
#let trn    = csv("results/cp_transfer.csv", row-type: dictionary)
#let evo    = csv("results/cp_event_order.csv", row-type: dictionary)
#let pg     = csv("results/pp_constant_groups.csv", row-type: dictionary)
#let pb     = csv("results/pp_constant_pace.csv", row-type: dictionary)
#let sel    = csv("results/selection_check.csv", row-type: dictionary)
#let strata = csv("results/slopes_by_stratum.csv", row-type: dictionary)
#let shape  = key.paul_shape

#let b0 = key.config.beta0
#let c500 = key.config.pp_c500

// --------------------------------------------------------------- helpers ----
#let num(x) = if x == none or x == "NA" or x == "" { none } else if type(x) == str { float(x) } else { x }
#let f(x, d: 3) = {
  let v = num(x)
  if v == none { return "–" }
  let s = str(calc.round(v, digits: d))
  let parts = s.split(".")
  let dec = if parts.len() > 1 { parts.at(1) } else { "" }
  while dec.len() < d { dec += "0" }
  let r = parts.at(0) + if d > 0 { "." + dec } else { "" }
  r.replace("-", "−")
}
#let n(x) = {
  let v = num(x)
  if v == none { return "–" }
  let s = str(int(v))
  let k = s.len()
  let out = ""
  for (i, c) in s.clusters().enumerate() {
    if i > 0 and calc.rem(k - i, 3) == 0 { out += "," }
    out += c
  }
  out
}
#let pct(x, d: 1) = f(x, d: d) + "%"
#let ci(lo, hi, d: 2) = "[" + f(lo, d: d) + ", " + f(hi, d: d) + "]"
#let mmss(x, d: 0) = {
  let m = calc.floor(x / 60)
  let s = x - 60 * m
  str(m) + ":" + (if s < 10 { "0" } else { "" }) + f(s, d: d)
}
#let lv(g) = {
  let s = g.replace("<=", "≤").replace(">=", "≥")
  s = s.replace(regex("(\\d)-(\\d)"), m => m.captures.at(0) + "–" + m.captures.at(1))
  if s == "F" { "Women" } else if s == "M" { "Men" } else if s == "all" { "All athletes" }
  else if s.starts-with("F ") { "Women " + s.slice(2) } else if s.starts-with("M ") { "Men " + s.slice(2) }
  else { let cs = s.clusters(); upper(cs.first()) + cs.slice(1).join() }
}
#let tbl(cols, head, body, align: auto, size: 9pt) = text(size: size, table(
  columns: cols, align: align, stroke: none, inset: (x: 4pt, y: 3pt),
  table.hline(), table.header(..head.map(h => [*#h*])), table.hline(stroke: 0.5pt),
  ..body.flatten(), table.hline()))

#let mlab = ("PL1": "PL1 (exponent 1.064)", "PL1cal": "PL1 (calibrated exponent)", "PL1_2k": "PL1 (from 2 km only)",
  "PL2": "PL2 (free exponent)", "PP1": "PP1 (Paul, 5 s per 500 m)", "PPcal": "PP (calibrated constant)",
  "PP_2k": "PP (from 2 km only)", "PP2": "PP2 (free constant)", "CP": "CP (classical fit)",
  "CPlog": "CP (log-time fit)", "CPr": "CP (2–20 min efforts only)", "CS": "Critical speed")
#let mshort = ("PL1": "PL1", "PL1cal": "PL1cal", "PL1_2k": "PL1 2k", "PL2": "PL2", "PP1": "PP1", "PPcal": "PPcal",
  "PP_2k": "PP 2k", "PP2": "PP2", "CP": "CP", "CPlog": "CP log", "CPr": "CP 2–20", "CS": "CS")
#let main-models = ("PL1", "PL1cal", "PL2", "PP1", "PPcal", "PP2", "CP", "CPlog", "CS")

#let rm(model) = sm.find(r => r.model == model)
#let rb(model) = sb.find(r => r.model == model)
#let rs(model) = ss.find(r => r.model == model)
#let pr(a, b) = pairs.find(r => r.a == a and r.b == b)
// verdict of a model against PL1 (paired bootstrap interval of the difference)
#let vs(r) = if num(r.diff_hi) < 0 [more accurate than] else if num(r.diff_lo) > 0 [less accurate than] else [not clearly different from]
#let vsp(r) = if num(r.hi) < 0 [more accurate than] else if num(r.lo) > 0 [less accurate than] else [not clearly different from]
#let dtext(r, base) = [#f(r.mdape, d: 2)% against #f(base.mdape, d: 2)%; difference #f(r.diff, d: 2) percentage points #ci(r.diff_lo, r.diff_hi)]
#let ptext(r) = [#f(r.mdape_a, d: 2)% against #f(r.mdape_b, d: 2)%; difference #f(r.diff, d: 2) percentage points #ci(r.lo, r.hi)]

#let ranked = sm.sorted(key: r => num(r.mdape))
#let best = ranked.first()
#let rank-of(model) = ranked.position(r => r.model == model) + 1
#let pp-wins-fixed = num(pr("PP1", "PL1").hi) < 0
#let pp-loses-fixed = num(pr("PP1", "PL1").lo) > 0
#let ppcal-wins = num(pr("PPcal", "PL1cal").hi) < 0
#let ppcal-loses = num(pr("PPcal", "PL1cal").lo) > 0
#let pp2-wins = num(pr("PP2", "PL2").hi) < 0
#let pp2-loses = num(pr("PP2", "PL2").lo) > 0
#let cp-wins = (num(rm("CP").diff_hi) < 0) or (num(rm("CPlog").diff_hi) < 0)
#let cp-loses = (num(rm("CP").diff_lo) > 0) and (num(rm("CPlog").diff_lo) > 0)
#let prop-fits = num(shape.sse_prop) < num(shape.sse_const)
#let all-row = pg.find(r => r.group == "all")

#let bdv(dim, level, model) = {
  let r = bd.find(r => r.dimension == dim and r.level == level and r.model == model)
  if r == none { none } else { num(r.mdape) }
}
#let best-bd(dim, level) = bd.filter(r => r.dimension == dim and r.level == level).sorted(key: r => num(r.mdape)).first()
#let present(dim, levels) = levels.filter(l => bd.any(r => r.dimension == dim and r.level == l))
#let bd-rows(dim, levels) = levels.map(l => {
  let vals = main-models.map(m => bdv(dim, l, m))
  let mn = calc.min(..vals.filter(v => v != none))
  let cnt = bd.find(r => r.dimension == dim and r.level == l).n
  (lv(l), n(cnt)) + vals.map(v => if v == none [–] else if v == mn [*#f(v, d: 2)*] else [#f(v, d: 2)])
})
#let bd-table(dim, levels, first-col) = tbl(
  (auto, auto, ..main-models.map(_ => auto)),
  (first-col, "n", ..main-models.map(m => mshort.at(m))),
  bd-rows(dim, present(dim, levels)), align: (left, right, ..main-models.map(_ => right)), size: 8.5pt)

#let ev-levels = evo.map(r => r.event)
#let positions = ("interpolation", "extrapolation to shorter", "extrapolation to longer")
#let age-levels = ("<=18", "19-29", "30-39", "40-49", "50-59", "60+")

// Paul's Law expressed as a power law
#let local-exp(p500) = 1 + c500 / (calc.ln(2) * p500)
#let eff-exp(p500, r) = 1 + calc.ln(1 + c500 * calc.log(r, base: 2) / p500) / calc.ln(r)

// ----------------------------------------------------------------- title ----
#align(center)[
  #text(17pt, weight: "bold")[Power law, Paul's Law or critical power?]
  #v(0.2em)
  #text(12pt)[Predicting erg performance with one- and two-parameter models]
  #v(0.2em)
  #text(10pt)[Verified Concept2 RowErg rankings, season #key.config.seasons]
  #v(0.2em)
  #text(9pt, fill: gray)[Generated #datetime.today().display("[year]-[month]-[day]")]
]
#v(0.8em)

= Summary

The power law validated in the companion analysis, $t = k thin d^(#f(b0))$, has one free parameter per athlete (the level $k$). Paul's Law is its best-known rival: pace rises by a fixed number of seconds, #c500 s per 500 m, for every doubling of distance. We compare both with each other and with the two-parameter critical-power (CP) model, using the same #n(key.data.athlete_seasons) athlete-seasons as the companion analysis. Each model is fitted to an athlete's other efforts and asked to predict the time for the distance of a held-out effort (leave-one-event-out); #n(key.data.predictions_main) predictions are common to the nine main models.

- *Overall accuracy.* The lowest median absolute error of the predicted time is #f(best.mdape, d: 2)%, for #mlab.at(best.model). PL1 has #f(rm("PL1").mdape, d: 2)% (rank #rank-of("PL1") of #sm.len()), Paul's Law with its published constant PP1 has #f(rm("PP1").mdape, d: 2)% (rank #rank-of("PP1")), and the best CP variant #f(calc.min(num(rm("CP").mdape), num(rm("CPlog").mdape)), d: 2)%.
- *Paul's Law against the power law.* With the published constant, Paul's Law is #vsp(pr("PP1", "PL1")) the power law with exponent #f(b0) (#ptext(pr("PP1", "PL1"))). With the constant (or exponent) calibrated on other athletes by sex and age group, Paul's Law is #vsp(pr("PPcal", "PL1cal")) the power law (#ptext(pr("PPcal", "PL1cal"))). With a free constant or exponent per athlete (two parameters), Paul's Law is #vsp(pr("PP2", "PL2")) the power law (#ptext(pr("PP2", "PL2"))).
- *The constant.* The constant that best fits these athletes is #f(all-row.c500, d: 2) s per 500 m per doubling #ci(all-row.lo, all-row.hi) against the published #c500 s. It #if prop-fits [grows with pace, as a power law implies, rather than staying fixed] else [stays close to fixed across paces, as Paul's Law assumes]: a constant fits the pace bins with weighted SSE #f(shape.sse_const, d: 0), a constant proportional to pace (a power law, exponent #f(shape.beta_prop)) with #f(shape.sse_prop, d: 0).
- *From the 2 km result alone,* Paul's Law is #vsp(pr("PP_2k", "PL1_2k")) the power law (#ptext(pr("PP_2k", "PL1_2k"))).
- *Critical power.* The classical CP fit is #vs(rm("CP")) PL1 (#dtext(rm("CP"), rm("PL1"))); fitted on log time like the other models, CP is #vs(rm("CPlog")) PL1 (#dtext(rm("CPlog"), rm("PL1"))).

= Paul's Law as a power law

Paul's Law in its usual form predicts the time $t_1$ over distance $d_1$ from a known effort $t_0$ over $d_0$:

$ t_1 = d_1 (t_0 / d_0 + 1/100 thin log_2 (d_1 / d_0)) . $

Written in pace $p = t / d$ (seconds per metre) it adds $1 slash 100$ s per metre for each doubling of distance, i.e. #c500 s per 500 m. The power law multiplies the pace instead: $p_1 = p_0 (d_1 slash d_0)^(beta - 1)$. For a distance ratio $r = d_1 slash d_0$ close to 1, $r^(beta - 1) approx 1 + (beta - 1) ln r$, so the two laws agree when

$ beta = 1 + c / (p_0 ln 2) , $

with $c$ Paul's constant (#c500 s per 500 m) and $p_0$ the athlete's pace. Paul's Law is therefore a power law whose exponent depends on the athlete's pace: the slower the athlete, the smaller the exponent. An exponent of #f(b0) corresponds to a pace of #f(c500 / (calc.ln(2) * (b0 - 1)), d: 1) s per 500 m, which is a 2 km time of #mmss(4 * c500 / (calc.ln(2) * (b0 - 1)), d: 1).

The equivalence is local. Paul's pace increase is linear in $ln d$, the power law's is exponential, so over a finite range the equivalent exponent is a little smaller than the local one. @equiv gives both.

#figure(
  tbl((auto, auto, auto, auto, auto, auto), ("2 km time", "Pace (s/500 m)", "Local exponent", "2 km → 6 km", "2 km → 10 km", "2 km → marathon"),
    (360, 390, 420, 450, 480, 510, 540, 570, 600).map(t => (mmss(t), f(t / 4, d: 1), f(local-exp(t / 4)),
      f(eff-exp(t / 4, 3)), f(eff-exp(t / 4, 5)), f(eff-exp(t / 4, 42195 / 2000)))),
    align: (right, right, right, right, right, right)),
  caption: [Power-law exponent equivalent to Paul's Law (#c500 s per 500 m per doubling) for an athlete with a given 2 km time: local exponent and exponent that gives the same prediction over a given distance ratio.],
) <equiv>

Two consequences can be checked against the data. First, Paul's Law predicts that exponents are smaller for slower athletes. @sexexp compares the exponent Paul's Law implies for the median 2 km time of the analysed men and women with the pooled exponents estimated in the companion analysis.

#figure(
  tbl((auto, auto, auto, auto, auto, auto), ("Sex", "Median 2 km time (s)", "Pace (s/500 m)", "Paul-implied exponent", "Observed pooled exponent", "Observed − implied"),
    ("F", "M").map(sx => {
      let t2 = num(sel.find(r => r.group.starts-with("athletes in the analysis") and r.sex == sx).median_2k_time_s)
      let ob = num(strata.find(r => r.group == sx).pooled_slope)
      (lv(sx), f(t2, d: 1), f(t2 / 4, d: 1), f(local-exp(t2 / 4)), f(ob), f(ob - local-exp(t2 / 4)))
    }),
    align: (left, right, right, right, right, right)),
  caption: [Exponent implied by Paul's Law for the median analysed athlete of each sex, against the observed pooled exponent.],
) <sexexp>

Second, Paul's Law makes the exponent fall with distance: the local exponent changes by $-(beta - 1)^2$ per unit of $ln d$ even without any extra parameter, a built-in concavity of the log-log curve that a constant-exponent power law lacks. The companion analysis found such a concavity in the data.

= Models

All models describe one athlete in one season and are fitted separately for each athlete. Average power on the Concept2 RowErg is $P = 2.8 thin (d \/ t)^3$; the constant 2.8 does not affect any prediction because the CP parameters rescale with it.

#figure(
  tbl((auto, auto, auto, auto), ("Model", "Free parameters", "Prediction", "Fitting"), (
    ("PL1", "1 (level)", [$t = k thin d^(1.064)$], "geometric mean of the training efforts"),
    ("PL1 calibrated", "1 (level)", [$t = k thin d^(beta_g)$], [$beta_g$: pooled exponent of the athlete's sex and age group from other athletes only (10-fold cross-fitting)]),
    ("PL1 from 2 km", "1 (2 km time)", [$t = t_(2000) (d \/ 2000)^(1.064)$], "the 2 km result alone"),
    ("PL2", "2", [$ln t = a + b ln(d \/ 2000)$], "least squares on log time"),
    ("PP1", "1 (level)", [$p = p_(2000) + c log_2 (d \/ 2000)$, $c$ = 5 s per 500 m], "level by least squares on log time"),
    ("PP calibrated", "1 (level)", [as PP1 with $c_g$], [$c_g$: pooled within-athlete regression of pace on $log_2 d$ for the athlete's sex and age group, from other athletes only]),
    ("PP from 2 km", "1 (2 km time)", [$t_1 = d_1 (t_0 \/ d_0 + 1/100 log_2 (d_1 \/ d_0))$], "the 2 km result alone"),
    ("PP2", "2 (level, $c$)", "as PP1", "least squares on log time"),
    ("CP, classical", "2 (CP, W′)", [$P = "CP" + W' \/ t$, i.e. $"CP" t^3 + W' t^2 = 2.8 d^3$], [regression of work $W = P t$ on $t$, CP and W′ ≥ 0]),
    ("CP, log-time fit", "2 (CP, W′)", "as above", "least squares on log time"),
    ("CP, 2–20 min", "2 (CP, W′)", "as above", "classical fit on training efforts of 2–20 minutes (at least two needed)"),
    ("Critical speed", "2 (CS, D′)", [$d = D' + "CS" thin t$], [regression of $d$ on $t$, D′ ≥ 0]),
  ), align: (left, left, left, left), size: 8.5pt),
  caption: [Models compared.],
)

The CP cubic is solved numerically for the time at the held-out distance, so every model answers the same question: _how long would this athlete take for this distance?_ A predicted pace below 0.01 s per metre (an impossible 5 s per 500 m) is floored, so that a failed Paul fit counts as a large error rather than being dropped. The calibrated variants use constants estimated on other athletes only, so they cost the athlete's own held-out result nothing. Where a one-parameter power law and a one-parameter Paul's Law are compared, they are fitted the same way: level by least squares on log time, the shape constant fixed.

= Evaluation design

*Data.* The same athlete-seasons as the companion analysis: verified individual RowErg results from season #key.config.seasons with at least #key.config.min_events events spanning a distance ratio of at least #key.config.min_ratio (#n(key.data.athlete_seasons) athlete-seasons, #n(key.data.results) results).

*Leave-one-event-out.* Each event is held out in turn, every model is fitted to the athlete's remaining events and predicts the held-out result, giving #n(key.data.predictions) predictions. Two-parameter models therefore have at least two training efforts. For a timed event (4, 30 or 60 minutes) the held-out result is the distance rowed in a fixed time; the model predicts the time for that distance and is compared with the fixed time.

*Error and uncertainty.* The error is the predicted time divided by the observed time minus one, in percent. We report the median absolute error, its 90th percentile, the median signed error (bias) and the share of predictions within 1% and 2%. Intervals come from a cluster bootstrap (#n(key.config.n_boot) resamples) that resamples athletes. Model comparisons are paired: the difference of median absolute errors is bootstrapped on the same predictions.

*Interpolation or extrapolation.* A held-out result is an interpolation if its duration lies between the shortest and longest training efforts. #pct(100 * key.data.share_interpolation, d: 0) of predictions are interpolations, #pct(100 * key.data.share_extrap_short, d: 0) extrapolate to shorter and #pct(100 * key.data.share_extrap_long, d: 0) to longer efforts.

*Prediction sets.* The nine main models are compared on the #n(key.data.predictions_main) predictions for which all are defined. The models that need a 2 km result are compared on the #n(key.data.predictions_single) predictions where the athlete has a 2 km result among the training efforts, and CP restricted to 2–20 minute efforts on the #n(key.data.predictions_best) predictions where it exists as well.

= Results

== Overall accuracy

#figure(
  tbl((auto, auto, auto, auto, auto, auto, auto, auto), ("Model", "Median error", "90th pct", "Bias", "Within 1%", "Within 2%", "Difference vs PL1", "Wins vs PL1"),
    ranked.map(r => (mlab.at(r.model), pct(r.mdape, d: 2), pct(r.p90, d: 1), pct(r.bias, d: 2),
      pct(r.within1, d: 0), pct(r.within2, d: 0),
      if r.model == "PL1" [–] else [#f(r.diff, d: 2) #ci(r.diff_lo, r.diff_hi)],
      if r.model == "PL1" [–] else [#pct(r.win_vs_ref, d: 0)])),
    align: (left, right, right, right, right, right, left, right), size: 8.5pt),
  caption: [Out-of-sample accuracy on #n(key.data.predictions_main) predictions (median absolute error of the predicted time; difference in percentage points with 95% bootstrap interval; "wins" is the share of predictions, excluding ties, where the model is closer than PL1).],
) <acc>

#figure(image("results/fig_cp_summary.png", width: 100%),
  caption: [Median absolute error with 95% intervals (left) and paired difference from PL1 (right).]) <fig-summary>

@fig-examples shows what the fitted curves look like for four athlete-seasons.

#figure(image("results/fig_cp_examples.png", width: 88%),
  caption: [Four randomly chosen athlete-seasons: measured average power against duration with the PL1, PP1, PL2 and CP (log-time fit) curves fitted to all of the athlete's events.]) <fig-examples>

== Like-for-like comparisons

Models should be compared at equal complexity. @pairs pairs the Paul's Law version of each model with its power-law counterpart, and the best two-parameter models with each other.

#let plabel = ("PP1-PL1": "Paul, fixed constant vs power law, fixed exponent (one parameter each)",
  "PPcal-PL1cal": "Paul vs power law, both calibrated by sex and age group (one parameter each)",
  "PP2-PL2": "Paul vs power law, both free per athlete (two parameters each)",
  "PP2-CPlog": "Paul, free constant vs critical power (two parameters each)",
  "CPlog-PL2": "Critical power vs power law, free exponent (two parameters each)",
  "PP_2k-PL1_2k": "Paul vs power law from the 2 km result alone (one parameter each)")
#figure(
  tbl((auto, auto, auto, auto, auto), ("Comparison (a vs b)", "n", "Error a vs b", "Difference a − b", "a closer"),
    pairs.map(r => (plabel.at(r.a + "-" + r.b), n(r.n), [#pct(r.mdape_a, d: 2) vs #pct(r.mdape_b, d: 2)],
      [#f(r.diff, d: 2) #ci(r.lo, r.hi)], pct(r.win_a, d: 0))),
    align: (left, right, left, left, right), size: 8.5pt),
  caption: [Paired comparisons of median absolute error (percentage points; 95% bootstrap interval; negative: the first model is more accurate; "a closer": share of predictions where a beats b).],
) <pairs>

== The constant in Paul's Law

The published constant of #c500 s per 500 m need not be optimal. We estimate it directly, as the pooled within-athlete slope of pace on $log_2$ distance, by group (@pgroups).

#figure(
  tbl((auto, auto, auto, auto), ("Group", "Athletes", "Constant (s per 500 m per doubling)", "95% CI"),
    pg.map(r => (lv(r.group), n(r.n_athletes), f(r.c500, d: 2), ci(r.lo, r.hi))),
    align: (left, right, right, left), size: 8.5pt),
  caption: [Best-fitting Paul constant by group, estimated on the data (published value: #c500).],
) <pgroups>

A decisive question is whether the constant is the same for everyone, as Paul's Law assumes, or grows with pace, as a power law implies. Athletes were binned by their 2 km-equivalent pace and the constant estimated within each bin (@pbins, @fig-pace).

#figure(image("results/fig_pp_constant.png", width: 80%),
  caption: [Added seconds per 500 m per doubling of distance, by the athletes' 2 km-equivalent pace, with 95% bootstrap intervals. Dashed: Paul's Law (constant). Solid: power law with exponent #f(b0) (proportional to pace).]) <fig-pace>

#figure(
  tbl((auto, auto, auto, auto, auto), ("Pace bin (s/500 m)", "Athletes", "Median pace", "Constant (95% CI)", "Equivalent exponent"),
    pb.map(r => (r.bin.replace(",", ", "), n(r.n_athletes), f(r.median_pace500, d: 1), [#f(r.c500, d: 2) #ci(r.lo, r.hi)], f(r.implied_exponent))),
    align: (left, right, right, left, right), size: 8.5pt),
  caption: [Paul constant by pace bin, and the power-law exponent that corresponds to it at the bin's median pace.],
) <pbins>

Fitting the #n(shape.n_bins) bins with a single constant gives #f(shape.const_c500, d: 2) s and a weighted sum of squares of #f(shape.sse_const, d: 0); fitting a constant proportional to pace (a power law with exponent #f(shape.beta_prop)) gives #f(shape.sse_prop, d: 0). #if prop-fits [The proportional shape fits better, so the added time per doubling grows with pace and a power law describes these data better than a single Paul constant.] else [The constant shape fits better, so the added time per doubling does not grow with pace and Paul's Law describes these data better than a power law.] The bins are few and the weights come from bootstrap intervals, so this is a descriptive comparison, not a formal test.

== From the 2 km result alone

Paul's Law was designed to predict from one known effort. @single compares the one-effort versions with the versions that use all the athlete's other efforts.

#figure(
  tbl((auto, auto, auto, auto, auto), ("Model", "Median error", "Bias", "Difference vs PL1", "Wins vs PL1"),
    ss.sorted(key: r => num(r.mdape)).map(r => (mlab.at(r.model), pct(r.mdape, d: 2), pct(r.bias, d: 2),
      if r.model == "PL1" [–] else [#f(r.diff, d: 2) #ci(r.diff_lo, r.diff_hi)],
      if r.model == "PL1" [–] else [#pct(r.win_vs_ref, d: 0)])),
    align: (left, right, right, left, right), size: 8.5pt),
  caption: [Predictions where the athlete has a 2 km result among the training efforts (#n(key.data.predictions_single) predictions). "From 2 km only" uses just that result; the other models use all other efforts.],
) <single>

#let from2k = trn.filter(r => r.source_lbl == "2000 m").sorted(key: r => ev-levels.position(e => e == r.target_lbl))
#if from2k.len() > 0 [
  Using only the 2 km result as the source, @from2k lists the error for each target event with the published constants of both laws.

  #figure(
    tbl((auto, auto, auto, auto, auto, auto), ("Target event", "n", "PL1 error", "PL1 bias", "Paul error", "Paul bias"),
      from2k.map(r => (r.target_lbl, n(r.n), pct(r.mdape, d: 2), pct(r.bias, d: 2), pct(r.mdape_pp, d: 2), pct(r.bias_pp, d: 2))),
      align: (left, right, right, right, right, right), size: 8.5pt),
    caption: [Median absolute error and median signed error when the 2 km result alone predicts each other event (published constants; positive bias: predicted too slow).],
  ) <from2k>
]

#figure(image("results/fig_cp_transfer.png", width: 70%),
  caption: [PL1 (exponent 1.064) predicting one event from another single event: median absolute error (%).]) <fig-transfer>

#figure(image("results/fig_pp_transfer.png", width: 70%),
  caption: [Paul's Law minus PL1, same sources and targets (percentage points of median absolute error; blue: Paul's Law more accurate).]) <fig-transfer-diff>

== Interpolation and extrapolation

#figure(image("results/fig_cp_position.png", width: 90%),
  caption: [Median absolute error by whether the held-out effort lies inside or outside the range of the training efforts.]) <fig-position>

#figure(
  bd-table("position", positions, "Held-out effort"),
  caption: [Median absolute error (%) by position of the held-out effort relative to the training efforts. Bold: lowest in the row.],
) <pos>

#figure(
  bd-table("training spread", ("< 2", "2-4", "4-10", ">= 10"), "Longest / shortest training effort"),
  caption: [Median absolute error (%) by spread of the training efforts' durations. Bold: lowest in the row.],
) <spread>

#figure(
  bd-table("training events", ("2", "3", "4+"), "Training efforts"),
  caption: [Median absolute error (%) by number of training efforts. Bold: lowest in the row.],
) <ntrain>

== By held-out event

#figure(image("results/fig_cp_event.png", width: 88%),
  caption: [Median absolute error for each held-out event, ordered from shortest to longest typical duration.]) <fig-event>

#figure(
  bd-table("target event", ev-levels, "Held-out event"),
  caption: [Median absolute error (%) by held-out event. Bold: lowest in the row.],
) <event>

The lowest-error model differs by event: #ev-levels.map(e => [#e: #mshort.at(best-bd("target event", e).model)]).join(", ").

== Critical power on its usual range

The CP model is usually applied to efforts of about 2 to 20 minutes. @best adds the variant fitted only to such efforts.

#figure(
  tbl((auto, auto, auto, auto, auto), ("Model", "Median error", "Bias", "Difference vs PL1", "Wins vs PL1"),
    sb.sorted(key: r => num(r.mdape)).map(r => (mlab.at(r.model), pct(r.mdape, d: 2), pct(r.bias, d: 2),
      if r.model == "PL1" [–] else [#f(r.diff, d: 2) #ci(r.diff_lo, r.diff_hi)],
      if r.model == "PL1" [–] else [#pct(r.win_vs_ref, d: 0)])),
    align: (left, right, right, left, right), size: 8.5pt),
  caption: [The main models and CP restricted to 2–20 minute training efforts, on the #n(key.data.predictions_best) predictions where all are defined.],
) <best>

On this set CP restricted to 2–20 minute efforts is #vs(rb("CPr")) PL1 (#dtext(rb("CPr"), rb("PL1"))).

== Fit in the athlete's own data and fitted parameters

@insample asks how well each model describes an athlete's efforts in the data it was fitted to, adjusted for the number of parameters (residual standard deviation of log time with $n - p$ degrees of freedom).

#figure(
  tbl((auto, auto, auto, auto), ("Model", "Parameters", "Median residual SD", "Best two-parameter fit"),
    ins.map(r => (mlab.at(r.model), r.params, pct(r.median_resid_sd_pct, d: 2),
      if num(r.best_two_param_share_pct) == none { "–" } else { pct(r.best_two_param_share_pct, d: 0) })),
    align: (left, right, right, right), size: 8.5pt),
  caption: [In-sample residual scatter and, among the two-parameter models, the share of athlete-seasons in which each fits best.],
) <insample>

Fitting each athlete-season to all of its events, the median exponent of the free power law is #f(par.median_pl2_exponent) (interquartile range #f(par.pl2_exponent_q25)–#f(par.pl2_exponent_q75)); the median free Paul constant is #f(par.median_pp2_c500, d: 2) s per 500 m per doubling (interquartile range #f(par.pp2_c500_q25, d: 2)–#f(par.pp2_c500_q75, d: 2)). The median critical power is #f(par.median_cp_w, d: 0) W and the median W′ is #f(par.median_wprime_kj, d: 1) kJ (interquartile range #f(par.wprime_q25_kj, d: 1)–#f(par.wprime_q75_kj, d: 1) kJ); #pct(par.boundary_share_pct, d: 1) of classical CP fits had CP or W′ at the zero boundary.

== Sex and age

#figure(
  bd-table("sex", ("F", "M"), "Sex"),
  caption: [Median absolute error (%) by sex. Bold: lowest in the row.],
) <sex>

#figure(
  bd-table("age group", age-levels, "Age group"),
  caption: [Median absolute error (%) by age group. Bold: lowest in the row.],
) <age>

= Discussion

*Paul's Law and the power law.* #if pp-wins-fixed [With the published constants Paul's Law predicts held-out efforts better than the power law with exponent #f(b0).] else if pp-loses-fixed [With the published constants Paul's Law predicts held-out efforts worse than the power law with exponent #f(b0).] else [With the published constants Paul's Law and the power law with exponent #f(b0) are not clearly different.] #if ppcal-wins [After calibrating both on other athletes, Paul's Law remains the more accurate.] else if ppcal-loses [After calibrating both on other athletes, the power law is the more accurate.] else [After calibrating both on other athletes, the two are not clearly different.] #if pp2-wins [When each athlete gets a free exponent or constant, Paul's Law is again more accurate.] else if pp2-loses [When each athlete gets a free exponent or constant, the power law is more accurate.] else [When each athlete gets a free exponent or constant, the two are not clearly different.] The laws are close over a limited range of distances and athlete paces (@equiv), so small differences are to be expected; they differ most when the athlete's pace is far from the pace at which Paul's constant happens to match the power law's exponent.

*Is the constant constant?* The data on the constant (@pbins) address the one assumption on which the two laws differ: #if prop-fits [added time per doubling grows with pace, so a fixed #c500 s over-corrects fast athletes' paces and under-corrects slow athletes', which is what a power law predicts.] else [added time per doubling stays roughly fixed across paces, as Paul's Law assumes.] #if prop-fits [The age and sex pattern of exponents in the companion analysis is consistent with this: slower athletes (older athletes and women on average) have smaller power-law exponents under a pace-dependent law, but a pace-proportional constant cannot by itself be shown to explain it.]

*Calibration.* The constants of both laws were treated as tunable; calibrated by sex and age group on other athletes they are free to adapt to the population. The comparison of the calibrated and uncalibrated versions in @acc shows how much of the error of each published law is due to its constant rather than its form.

*What a second parameter buys.* #if cp-wins [The critical-power model is more accurate than the one-parameter power law.] else if cp-loses [The critical-power model is less accurate than the one-parameter power law. CP implies a power that levels off at a critical value for long efforts and rises without bound for short ones, while the power laws decline smoothly over all durations.] else [The critical-power model is not clearly better or worse than the one-parameter power law.] The CP model is normally applied to efforts of about 2 to 20 minutes, whereas the events here run from about 1.5 to 60 minutes; the restricted fit in @best tests that.

*Interpolation and extrapolation.* A two-parameter model has more freedom and, with two or three training efforts, little data to pin it down. @pos and @spread show where that matters. Predictions that extrapolate beyond the training range are the usual practical need (for example a 10 km time from a 2 km result) and are where the choice of functional form matters most.

*Effort.* All models assume maximal efforts. A sub-maximal result can distort a fit to two or three efforts, either by being absorbed into a distorted curve or by being amplified in extrapolation. The analysis cannot separate this from true model error; the errors reported include both.

*Limitations.* The athletes are those with at least three verified events across a threefold range of distance, mostly male masters rowers, from one season. The comparison is on the typical error of one held-out result, not the best case for any athlete. The calibrated constants are cross-fitted by athlete but come from the same data source as the evaluation. The pace-bin analysis uses eight bins and bootstrap weights, and pace is the athlete's 2 km-equivalent pace under a fixed exponent. Intervals reflect sampling variation only.

= Conclusion

On #n(key.data.predictions_main) held-out efforts the model with the lowest median error is #mlab.at(best.model) (#f(best.mdape, d: 2)%); the one-parameter power law has #f(rm("PL1").mdape, d: 2)% and Paul's Law with the published constant #f(rm("PP1").mdape, d: 2)%. Paul's Law is, to first order, a power law whose exponent depends on pace (@equiv). In these data the added time per doubling #if prop-fits [grows with pace, which favours the power law's proportional form over Paul's fixed constant.] else [stays close to fixed across paces, which favours Paul's constant over the power law's proportional form.] #if cp-wins [A critical-power model is more accurate than the one-parameter laws on these data.] else if cp-loses [A critical-power model is less accurate than the one-parameter laws on these data.] else [A critical-power model is not clearly better or worse than the one-parameter laws on these data.] The one-parameter models need only a single effort, which counts in their favour when a second well-spaced maximal effort is not available. The comparison does not say that any model is physiologically correct, only how well each predicts the next effort.

= Appendix: reproducibility

Run `power_law_analysis.R` first (it prepares the data, defines the athlete-season eligibility rule and writes the results that @sexexp uses), then `cp_vs_powerlaw.R`, which sources it and writes `results/cp_*.csv`, `results/pp_*.csv`, `results/cp_key_numbers.json` and `results/fig_cp_*.png`, `results/fig_pp_*.png`. Compile this report with `typst compile report_cp.typ`. Settings: #n(key.config.n_boot) bootstrap resamples, #n(key.config.beta_folds)-fold cross-fitting of the calibrated constants, #n(key.config.pace_bins) pace bins, Paul's constant #key.config.pp_c500 s per 500 m, power constant #f(key.config.cp_const, d: 1), CP efforts restricted to #key.config.cp_window_min–#key.config.cp_window_max minutes in the restricted variant. The predictions file contains an anonymous athlete index rather than athlete identifiers.
