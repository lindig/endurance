// report_cp.typ — one parameter or two?
// Compares the one-parameter power-law model (t = k d^1.064) with two-parameter
// models (critical power, critical speed, free-exponent power law).
// Reads everything from results/ (written by cp_vs_powerlaw.R).
// Compile with:  typst compile report_cp.typ
// Numbers, tables and figures come from the results files; the interpretive
// sentences are conditional on those numbers but should still be re-read.

#set page(paper: "a4", margin: (x: 2.2cm, y: 2.3cm), numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#show heading: set block(above: 1.3em, below: 0.7em)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)

// ------------------------------------------------------------------ data ----
#let key = json("results/cp_key_numbers.json")
#let sm   = csv("results/cp_summary_main.csv", row-type: dictionary)
#let sb   = csv("results/cp_summary_best.csv", row-type: dictionary)
#let bd   = csv("results/cp_breakdowns.csv", row-type: dictionary)
#let ins  = csv("results/cp_insample.csv", row-type: dictionary)
#let par  = csv("results/cp_params.csv", row-type: dictionary).first()
#let trn  = csv("results/cp_transfer.csv", row-type: dictionary)
#let evo  = csv("results/cp_event_order.csv", row-type: dictionary)

#let b0 = key.config.beta0

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
#let lv(g) = {
  let s = g.replace("<=", "≤").replace(">=", "≥")
  s = s.replace(regex("(\\d)-(\\d)"), m => m.captures.at(0) + "–" + m.captures.at(1))
  if s == "F" { "Women" } else if s == "M" { "Men" }
  else { let cs = s.clusters(); upper(cs.first()) + cs.slice(1).join() }
}
#let tbl(cols, head, body, align: auto) = text(size: 9pt, table(
  columns: cols, align: align, stroke: none, inset: (x: 5pt, y: 3pt),
  table.hline(), table.header(..head.map(h => [*#h*])), table.hline(stroke: 0.5pt),
  ..body.flatten(), table.hline()))

#let mlab = ("PL1": "PL1 (exponent 1.064)", "PL1cal": "PL1 (calibrated exponent)", "PL1_2k": "PL1 (from 2 km only)",
  "PL2": "PL2 (free exponent)", "CP": "CP (classical fit)", "CPlog": "CP (log-time fit)",
  "CPr": "CP (2–20 min efforts only)", "CS": "Critical speed")
#let mshort = ("PL1": "PL1", "PL1cal": "PL1cal", "PL1_2k": "PL1 2k", "PL2": "PL2", "CP": "CP", "CPlog": "CP log", "CPr": "CP 2–20", "CS": "CS")
#let main-models = ("PL1", "PL1cal", "PL2", "CP", "CPlog", "CS")

#let rm(model) = sm.find(r => r.model == model)
#let rb(model) = sb.find(r => r.model == model)
// verdict of a model against PL1 (paired bootstrap interval of the difference)
#let vs(r) = if num(r.diff_hi) < 0 [more accurate than] else if num(r.diff_lo) > 0 [less accurate than] else [not clearly different from]
#let diff-text(r) = [#f(r.mdape, d: 2)% against #f(rm("PL1").mdape, d: 2)%; difference #f(r.diff, d: 2) percentage points #ci(r.diff_lo, r.diff_hi)]
#let diff-text-b(r) = [#f(r.mdape, d: 2)% against #f(rb("PL1").mdape, d: 2)%; difference #f(r.diff, d: 2) percentage points #ci(r.diff_lo, r.diff_hi)]

#let ranked = sm.sorted(key: r => num(r.mdape))
#let best = ranked.first()
#let pl1-rank = ranked.position(r => r.model == "PL1") + 1
#let two-param = ("PL2", "CP", "CPlog", "CS")
#let best2 = sm.filter(r => r.model in two-param).sorted(key: r => num(r.mdape)).first()
#let cp-wins = (num(rm("CP").diff_hi) < 0) or (num(rm("CPlog").diff_hi) < 0)
#let cp-loses = (num(rm("CP").diff_lo) > 0) and (num(rm("CPlog").diff_lo) > 0)
#let pl2-wins = num(rm("PL2").diff_hi) < 0
#let pl2-loses = num(rm("PL2").diff_lo) > 0

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
  bd-rows(dim, present(dim, levels)), align: (left, right, ..main-models.map(_ => right)))

#let ev-levels = evo.map(r => r.event)
#let positions = ("interpolation", "extrapolation to shorter", "extrapolation to longer")
#let age-levels = ("<=18", "19-29", "30-39", "40-49", "50-59", "60+")

// ----------------------------------------------------------------- title ----
#align(center)[
  #text(17pt, weight: "bold")[One parameter or two?]
  #v(0.2em)
  #text(12pt)[Predicting erg performance with a power law versus a critical-power model]
  #v(0.2em)
  #text(10pt)[Verified Concept2 RowErg rankings, season #key.config.seasons]
  #v(0.2em)
  #text(9pt, fill: gray)[Generated #datetime.today().display("[year]-[month]-[day]")]
]
#v(0.8em)

= Summary

The power law validated in the companion analysis, $t = k thin d^(#f(b0))$, has one free parameter per athlete (the level $k$). Given one effort it predicts the time for another distance. We ask whether a model with a second free parameter, above all the critical-power (CP) model, predicts a held-out effort better. Every model is fitted to an athlete's other efforts and asked to predict the time for the distance of the held-out one (leave-one-event-out), on the same #n(key.data.athlete_seasons) athlete-seasons as the companion analysis. There are #n(key.data.predictions_main) predictions common to the six main models.

- *Overall accuracy.* The lowest median absolute error of the predicted time is #f(best.mdape, d: 2)%, for #mlab.at(best.model). The one-parameter power law (PL1) has #f(rm("PL1").mdape, d: 2)% and ranks #pl1-rank of #sm.len().
- *Critical power against the power law.* The classical CP fit is #vs(rm("CP")) PL1 (#diff-text(rm("CP"))). Fitted by least squares on log time, the same loss as the power laws, CP is #vs(rm("CPlog")) PL1 (#diff-text(rm("CPlog"))).
- *A free exponent.* Letting the exponent vary by athlete (PL2) is #vs(rm("PL2")) PL1 (#diff-text(rm("PL2"))).
- *Where the models differ.* For interpolation the lowest error is #mlab.at(best-bd("position", "interpolation").model) (#f(best-bd("position", "interpolation").mdape, d: 2)%); when extrapolating to longer efforts it is #mlab.at(best-bd("position", "extrapolation to longer").model) (#f(best-bd("position", "extrapolation to longer").mdape, d: 2)%) and to shorter efforts #mlab.at(best-bd("position", "extrapolation to shorter").model) (#f(best-bd("position", "extrapolation to shorter").mdape, d: 2)%).
- *Best practice.* Restricting CP to 2–20 minute efforts gives #f(rb("CPr").mdape, d: 2)% and predicting from the 2 km result alone with PL1 gives #f(rb("PL1_2k").mdape, d: 2)%, against #f(rb("PL1").mdape, d: 2)% for PL1 on all other efforts (on the #n(key.data.predictions_best) predictions where all eight models exist).

= Models

All models describe one athlete in one season and are fitted separately for each athlete. Average power on the Concept2 RowErg is $P = 2.8 thin (d \/ t)^3$ (watts, metres, seconds); the constant 2.8 does not affect any prediction because the CP parameters rescale with it.

#figure(
  tbl((auto, auto, auto, auto), ("Model", "Free parameters", "Prediction", "Fitting"), (
    ("PL1", "1 (level)", [$t = k thin d^(1.064)$], "geometric mean of the training efforts"),
    ("PL1 calibrated", "1 (level)", [$t = k thin d^(beta_g)$], [$beta_g$ = pooled exponent of the athlete's sex and age group, estimated on other athletes only (10-fold cross-fitting)]),
    ("PL1 from 2 km", "1 (2 km time)", [$t = t_(2000) (d \/ 2000)^(1.064)$], "the athlete's 2 km result alone"),
    ("PL2", "2", [$ln t = a + b ln(d \/ 2000)$], "least squares on log time"),
    ("CP, classical", "2 (CP, W′)", [$P = "CP" + W' \/ t$, i.e. $"CP" t^3 + W' t^2 = 2.8 d^3$], [regression of work $W = P t$ on $t$, with CP and W′ ≥ 0]),
    ("CP, log-time fit", "2 (CP, W′)", "as above", "least squares on log time"),
    ("CP, 2–20 min", "2 (CP, W′)", "as above", "classical fit using only training efforts of 2–20 minutes (at least two needed)"),
    ("Critical speed", "2 (CS, D′)", [$d = D' + "CS" thin t$], [regression of $d$ on $t$, with D′ ≥ 0]),
  ), align: (left, left, left, left)),
  caption: [Models compared.],
)

The cubic $"CP" t^3 + W' t^2 = 2.8 d^3$ is solved numerically for the time at the held-out distance, so every model answers the same question: _how long would this athlete take for this distance?_ The power laws are one- or two-parameter descriptions of the same curve, so PL2 is the like-for-like two-parameter rival of CP. The log-time CP fit removes a possible disadvantage of the classical regression, which weights long efforts more heavily than the power-law fit does.

= Evaluation design

*Data.* The same athlete-seasons as the companion analysis: verified individual RowErg results from season #key.config.seasons with at least #key.config.min_events events spanning a distance ratio of at least #key.config.min_ratio (#n(key.data.athlete_seasons) athlete-seasons, #n(key.data.results) results).

*Leave-one-event-out.* For each athlete-season each event is held out in turn, every model is fitted to the remaining events and predicts the held-out result, giving #n(key.data.predictions) predictions. Two-parameter models therefore have at least two training efforts. For a timed event (4, 30 or 60 minutes) the held-out result is the distance rowed in a fixed time; the model predicts the time for that distance and is compared with the fixed time.

*Error.* The error is the predicted time divided by the observed time minus one, in percent. We report the median absolute error, its 90th percentile, the median signed error (bias) and the share of predictions within 1% and 2%. Intervals are from a cluster bootstrap (#n(key.config.n_boot) resamples) that resamples athletes. The comparison with PL1 is paired: the difference of median absolute errors is bootstrapped on the same predictions.

*Interpolation or extrapolation.* A held-out result is an interpolation if its duration lies between the shortest and longest training efforts, otherwise an extrapolation to shorter or longer efforts. #pct(100 * key.data.share_interpolation, d: 0) of predictions are interpolations, #pct(100 * key.data.share_extrap_short, d: 0) extrapolate to shorter and #pct(100 * key.data.share_extrap_long, d: 0) to longer efforts.

*Common prediction sets.* Models are compared on the same predictions. The six main models are compared on every prediction for which all are defined. The best-practice variants (CP restricted to 2–20 minutes, PL1 from 2 km only) need extra data and are compared on the smaller set where all eight models exist.

= Results

== Overall accuracy

#figure(
  tbl((auto, auto, auto, auto, auto, auto, auto, auto), ("Model", "Median error", "90th pct", "Bias", "Within 1%", "Within 2%", "Difference vs PL1", "Wins vs PL1"),
    ranked.map(r => (mlab.at(r.model), pct(r.mdape, d: 2), pct(r.p90, d: 1), pct(r.bias, d: 2),
      pct(r.within1, d: 0), pct(r.within2, d: 0),
      if r.model == "PL1" [–] else [#f(r.diff, d: 2) #ci(r.diff_lo, r.diff_hi)],
      if r.model == "PL1" [–] else [#pct(r.win_vs_ref, d: 0)])),
    align: (left, right, right, right, right, right, left, right)),
  caption: [Out-of-sample accuracy on #n(key.data.predictions_main) predictions (median absolute error of the predicted time; difference in percentage points with 95% bootstrap interval; "wins" is the share of predictions, excluding ties, where the model is closer than PL1).],
) <acc>

#figure(image("results/fig_cp_summary.png", width: 100%),
  caption: [Median absolute error with 95% intervals (left) and paired difference from PL1 (right).]) <fig-summary>

#if cp-wins [The critical-power model predicts held-out efforts better than the fixed-exponent power law (@acc).] else if cp-loses [The critical-power model predicts held-out efforts worse than the fixed-exponent power law (@acc).] else [The critical-power model is not clearly better or worse than the fixed-exponent power law (@acc).] #if pl2-wins [A power law with a free exponent is more accurate than the fixed 1.064.] else if pl2-loses [A power law with a free exponent is less accurate than the fixed 1.064, so the extra parameter costs more in noise than it gains.] else [A power law with a free exponent is not clearly different from the fixed 1.064.] The two-parameter model with the lowest error is #mlab.at(best2.model) (#f(best2.mdape, d: 2)%). @fig-examples shows what the fitted curves look like for four athlete-seasons.

#figure(image("results/fig_cp_examples.png", width: 88%),
  caption: [Four randomly chosen athlete-seasons: measured average power against duration with the PL1, PL2 and CP (log-time fit) curves fitted to all of the athlete's events.]) <fig-examples>

== Interpolation and extrapolation

#figure(image("results/fig_cp_position.png", width: 85%),
  caption: [Median absolute error by whether the held-out effort lies inside or outside the range of the training efforts.]) <fig-position>

#figure(
  bd-table("position", positions, "Held-out effort"),
  caption: [Median absolute error (%) by position of the held-out effort relative to the training efforts. Bold: lowest in the row.],
) <pos>

A second parameter can only help if the extra flexibility is identified by the training data. @spread splits the predictions by how widely the training efforts are spread in duration, and by how many there are.

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

== Best-practice variants

#figure(
  tbl((auto, auto, auto, auto, auto), ("Model", "Median error", "Bias", "Difference vs PL1", "Wins vs PL1"),
    sb.sorted(key: r => num(r.mdape)).map(r => (mlab.at(r.model), pct(r.mdape, d: 2), pct(r.bias, d: 2),
      if r.model == "PL1" [–] else [#f(r.diff, d: 2) #ci(r.diff_lo, r.diff_hi)],
      if r.model == "PL1" [–] else [#pct(r.win_vs_ref, d: 0)])),
    align: (left, right, right, left, right)),
  caption: [All eight models on the #n(key.data.predictions_best) predictions where all are defined (the athlete has a 2 km result among the training efforts and at least two training efforts of 2–20 minutes).],
) <best>

On this set CP restricted to 2–20 minute efforts is #vs(rb("CPr")) PL1 (#diff-text-b(rb("CPr"))), and predicting from the 2 km result alone is #vs(rb("PL1_2k")) PL1 using all other efforts (#diff-text-b(rb("PL1_2k"))).

== Fit in the athlete's own data and fitted parameters

The out-of-sample test asks about prediction. @insample asks how well each model describes an athlete's efforts in the data it was fitted to, adjusted for the number of parameters (residual standard deviation of log time with degrees of freedom $n - p$).

#figure(
  tbl((auto, auto, auto, auto), ("Model", "Parameters", "Median residual SD", "Best two-parameter fit"),
    ins.map(r => (mlab.at(r.model), r.params, pct(r.median_resid_sd_pct, d: 2),
      if num(r.best_two_param_share_pct) == none { "–" } else { pct(r.best_two_param_share_pct, d: 0) })),
    align: (left, right, right, right)),
  caption: [In-sample residual scatter and, among the two-parameter models, the share of athlete-seasons in which each fits best.],
) <insample>

Fitting each athlete-season to all of its events, the median exponent of the free power law is #f(par.median_pl2_exponent, d: 3) (interquartile range #f(par.pl2_exponent_q25, d: 3)–#f(par.pl2_exponent_q75, d: 3)). The median critical power is #f(par.median_cp_w, d: 0) W and the median W′ is #f(par.median_wprime_kj, d: 1) kJ (interquartile range #f(par.wprime_q25_kj, d: 1)–#f(par.wprime_q75_kj, d: 1) kJ). #pct(par.boundary_share_pct, d: 1) of classical CP fits had CP or W′ at the zero boundary.

== One effort is enough for the power law

The power law needs only one effort. @fig-transfer shows how well PL1 transfers a single known effort (rows) to every other event (columns) for athletes who have both.

#figure(image("results/fig_cp_transfer.png", width: 78%),
  caption: [Median absolute error (%) when PL1 (exponent 1.064) predicts the target event from a single source event.]) <fig-transfer>

#let from2k = trn.filter(r => r.source_lbl == "2000 m")
#if from2k.len() > 0 [
  From the 2 km result alone the median absolute error is between #f(calc.min(..from2k.map(r => num(r.mdape))), d: 1)% and #f(calc.max(..from2k.map(r => num(r.mdape))), d: 1)% depending on the target event, and the median signed error runs from #f(calc.min(..from2k.map(r => num(r.bias))), d: 1)% to #f(calc.max(..from2k.map(r => num(r.bias))), d: 1)% (positive: predicted too slow).
]

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

*What the second parameter buys.* #if cp-wins [The critical-power model beats the fixed-exponent power law, so for these athletes and events a second parameter describing a limited anaerobic reserve on top of a sustainable power predicts better than a single curve with a fixed exponent.] else if cp-loses [The critical-power model predicts worse than the fixed-exponent power law. The two models differ in shape as well as parameter count: CP implies a power that levels off at a critical value for long efforts and rises without bound for short ones, while the power law implies a smooth decline across all durations.] else [Neither the critical-power model nor the fixed-exponent power law is clearly better.] #if pl2-wins [Letting the power-law exponent vary by athlete also helps, which is consistent with exponents that differ between athletes.] else if pl2-loses [Letting the exponent vary by athlete hurts: with only a few training efforts the estimated exponent is too noisy to beat the fixed value.] else [Letting the exponent vary by athlete makes no clear difference.]

*Interpolation and extrapolation.* A two-parameter model has more freedom and, with two or three training efforts, little data to pin it down. @pos and @spread show where that matters: whether the advantage or disadvantage of the two-parameter models sits in interpolation or in extrapolation, and whether it depends on how widely the training efforts are spread. Predictions that extrapolate beyond the training range are the usual practical need (a 2 km time from short tests, or a 10 km time from a 2 km result) and are where the choice of functional form matters most.

*Validity range of CP.* The CP model is normally applied to efforts of roughly 2 to 20 minutes, where the severe-intensity domain is assumed; outside that range it is known to misrepresent both very short efforts and very long ones. The events here run from about 1.5 minutes (500 m) to 60 minutes, so the unrestricted CP fits are used well outside their usual range. @best tests the restricted variant, and @event shows how the error depends on the event.

*Effort.* Both families assume maximal efforts. A single sub-maximal result can distort a two-parameter fit to two or three efforts, either by being absorbed into a distorted curve or by being amplified in extrapolation. The analysis cannot separate this from true model error; the errors reported here include both.

*Limitations.* The athletes are those with at least three verified events across a threefold range of distance, mostly male masters rowers, from one season. The comparison is on the typical error of one held-out result, not on the best case for any single athlete. The calibrated exponent is cross-fitted by athlete, but it uses the same data source as the evaluation. Performance of CP with more training efforts, with a purpose-designed testing protocol, or with efforts chosen inside its validity range may be better than reported here. Intervals reflect sampling variation only.

= Conclusion

On #n(key.data.predictions_main) held-out efforts the model with the lowest median error is #mlab.at(best.model) (#f(best.mdape, d: 2)%), while the one-parameter power law with exponent #f(b0) has #f(rm("PL1").mdape, d: 2)%. #if cp-wins [A critical-power model is more accurate than the one-parameter power law on these data] else if cp-loses [A critical-power model is less accurate than the one-parameter power law on these data] else [A critical-power model is not clearly better or worse than the one-parameter power law on these data]#if pl2-wins [, and a free exponent also improves on the fixed 1.064]#if pl2-loses [, and a free exponent does not improve on the fixed 1.064]. The one-parameter model needs only a single effort, which counts in its favour when a second well-spaced maximal effort is not available. The comparison does not say that either model is physiologically correct, only how well each predicts the next effort.

= Appendix: reproducibility

Run `power_law_analysis.R` first (it prepares the data and defines the athlete-season eligibility rule), then `cp_vs_powerlaw.R`, which sources it and writes `results/cp_*.csv`, `results/cp_key_numbers.json` and `results/fig_cp_*.png`. Compile this report with `typst compile report_cp.typ`. Settings: #n(key.config.n_boot) bootstrap resamples, #n(key.config.beta_folds)-fold cross-fitting of the calibrated exponent, CP efforts restricted to #key.config.cp_window_min–#key.config.cp_window_max minutes in the restricted variant, power constant #f(key.config.cp_const, d: 1). The predictions file contains an anonymous athlete index rather than athlete identifiers.
