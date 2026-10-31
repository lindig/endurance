// report.typ — does erg performance follow a 1.064 power law?
// Reads everything from results/ (written by power_law_analysis.R).
// Compile with:  typst compile report.typ
// All numbers, tables and figures come from the results files. The interpretive
// wording was written for the season-2026 run; re-read it if you change the data.

#set page(paper: "a4", margin: (x: 2.2cm, y: 2.3cm), numbering: "1")
#set text(size: 10.5pt)
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#show heading: set block(above: 1.3em, below: 0.7em)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)

// ------------------------------------------------------------------ data ----
#let key      = json("results/key_numbers.json")
#let funnel   = csv("results/funnel.csv", row-type: dictionary)
#let events   = csv("results/events_used.csv", row-type: dictionary)
#let selection = csv("results/selection_check.csv", row-type: dictionary)
#let strata   = csv("results/slopes_by_stratum.csv", row-type: dictionary)
#let cstrata  = csv("results/curvature_by_stratum.csv", row-type: dictionary)
#let awin     = csv("results/slopes_by_age_window.csv", row-type: dictionary)
#let windows  = csv("results/slopes_by_window.csv", row-type: dictionary)
#let anch     = csv("results/implied_exponent_by_event.csv", row-type: dictionary)
#let effort   = csv("results/effort_sensitivity.csv", row-type: dictionary)
#let sens     = csv("results/sensitivity_min_events.csv", row-type: dictionary)

#let b0 = key.config.beta0
#let P  = key.overall.pooled
#let Cv = key.curvature
#let E  = key.effort

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
#let ci(lo, hi, d: 3) = "[" + f(lo, d: d) + ", " + f(hi, d: d) + "]"
#let lbl(g) = {
  let s = g.replace("<=", "≤").replace(">=", "≥").replace(" - ", " – ")
  s = s.replace(regex("(\\d)-(\\d)"), m => m.captures.at(0) + "–" + m.captures.at(1))
  if s == "all" { "All athletes" }
  else if s == "F" { "Women" }
  else if s == "M" { "Men" }
  else if s.starts-with("F ") { "Women " + s.slice(2) }
  else if s.starts-with("M ") { "Men " + s.slice(2) }
  else { s }
}
#let row-of(rows, col, val) = rows.find(r => r.at(col) == val)
#let get(rows, col, val, field) = {
  let r = row-of(rows, col, val)
  if r == none { none } else { num(r.at(field)) }
}
#let tbl(cols, head, body, align: auto) = text(size: 9pt, table(
  columns: cols, align: align, stroke: none, inset: (x: 5pt, y: 3pt),
  table.hline(), table.header(..head.map(h => [*#h*])), table.hline(stroke: 0.5pt),
  ..body.flatten(), table.hline()))
#let aw(w, g) = num(awin.find(r => r.window == w and r.group == g).pooled_slope)
#let ratio(b, d) = calc.pow(d / 2000, b)
#let curved(d) = calc.exp(Cv.exponent_at_2km * calc.ln(d / 2000) + Cv.curvature * calc.pow(calc.ln(d / 2000), 2))

// derived numbers used in the text
#let age-groups = ("<=18", "19-29", "30-39", "40-49", "50-59", "60+")
#let p-young = get(strata, "group", "19-29", "pooled_slope")
#let p-old   = get(strata, "group", "60+", "pooled_slope")
#let e-young = get(cstrata, "group", "19-29", "exp2k")
#let e-old   = get(cstrata, "group", "60+", "exp2k")
#let dF = (get(strata, "group", "F", "pooled_slope"), get(strata, "group", "M", "pooled_slope"))
#let base = effort.at(0)
#let max-dev = calc.max(..effort.map(r => calc.abs(num(r.linear) - num(base.linear))))
#let age-grad     = p-young != none and p-old != none and p-young - p-old > 0.01
#let age-grad-win = aw("2 km - 21 km", "<=18") - aw("2 km - 21 km", "60+") > 0.01
#let women-lower  = dF.at(0) < dF.at(1) - 0.003
#let a4    = anch.find(r => r.event == "4")
#let a1000 = anch.find(r => r.event == "1000")
#let low4  = a4 != none and a1000 != none and num(a1000.implied_exponent) - num(a4.implied_exponent) > 0.01
#let win-range = calc.max(..windows.map(r => num(r.pooled_slope))) - calc.min(..windows.map(r => num(r.pooled_slope)))
#let sens-range = calc.max(..sens.map(r => num(r.pooled_slope))) - calc.min(..sens.map(r => num(r.pooled_slope)))
#let anch-sorted = anch.sorted(key: r => num(r.implied_exponent))
#let share-men  = num(row-of(strata, "group", "M").n_athletes) / num(row-of(strata, "group", "all").n_athletes)
#let share-over40 = (num(row-of(strata, "group", "40-49").n_athletes) + num(row-of(strata, "group", "50-59").n_athletes) + num(row-of(strata, "group", "60+").n_athletes)) / num(row-of(strata, "group", "all").n_athletes)
#let share-elig = key.data.athlete_seasons / num(funnel.at(funnel.len() - 2).athlete_seasons)
#let traits = {
  let t = ()
  if age-grad { t.push("falls with age") }
  if women-lower { t.push("is lower for women") }
  if Cv.curvature < -0.002 { t.push("decreases with distance") } else if Cv.curvature > 0.002 { t.push("increases with distance") }
  t
}
#let verdict = if P.lo <= b0 and b0 <= P.hi { "consistent with" } else if calc.abs(P.est - b0) > 0.01 { if P.est > b0 { "above" } else { "below" } } else { if P.est > b0 { "slightly above" } else { "slightly below" } }

// ----------------------------------------------------------------- title ----
#align(center)[
  #text(17pt, weight: "bold")[Does erg performance follow a 1.064 power law?]
  #v(0.2em)
  #text(12pt)[Evidence from verified Concept2 RowErg rankings, season #key.config.seasons]
  #v(0.2em)
  #text(9pt, fill: gray)[Generated #datetime.today().display("[year]-[month]-[day]")]
]
#v(0.8em)

= Summary

We tested whether, for the same athlete in the same season, race time grows with distance as $t prop d^(1.064)$. The analysis uses #n(key.data.athlete_seasons) athlete-seasons (#n(key.data.rows_used) verified results) from the Concept2 RowErg rankings, covering both sexes and all age groups.

- *On average the exponent is #verdict 1.064.* The pooled within-athlete exponent is #f(P.est) #ci(P.lo, P.hi). Using 1.064 instead of #f(P.est) changes a predicted 6 km time from a 2 km time by #pct(100 * (ratio(P.est, 6000) / ratio(b0, 6000) - 1), d: 2) and a predicted 10 km time by #pct(100 * (ratio(P.est, 10000) / ratio(b0, 10000) - 1), d: 2).
- *Does it vary between athletes?* The pooled exponent is #f(p-young) for ages 19–29 and #f(p-old) for 60+#if age-grad-win [, and this age gradient remains inside fixed distance windows and in the exponent at 2 km (#f(e-young) vs #f(e-old))]#if not age-grad-win [, with no clear gradient once distance windows are held fixed]. Women minus men is #f(dF.at(0) - dF.at(1)) over all ages.
- *It depends on distance.* The local exponent at 2 km is #f(Cv.exponent_at_2km) #ci(Cv.exponent_at_2km_lo, Cv.exponent_at_2km_hi) and it #if Cv.curvature < 0 [falls] else [rises] by #f(2 * Cv.curvature * calc.ln(2), d: 4) per doubling of distance (curvature #f(Cv.curvature, d: 4)). A single power law therefore fits a limited range of distances well and the pooled slope depends on which events are included.
- *Sub-maximal efforts are a limitation.* Across five alternative treatments of unusual results the pooled exponent moves by at most #f(max-dev, d: 3) from the baseline (@sec-effort), but the direction of any bias is unknown.

= Hypothesis and model

For one athlete, the hypothesis is that time depends on distance as a power law, $t = k thin d^beta$ with $beta = #f(b0)$. Equivalent forms are average pace $prop d^(beta - 1) = d^(#f(b0 - 1))$ and average power $prop d^(-3(beta - 1)) = d^(#f(-3 * (b0 - 1)))$. For each athlete-season $i$ with results $j$ we fit

$ ln t_(i j) = a_i + b thin ln (d_(i j) / "2000 m") + epsilon_(i j) , $

and test $b = beta$. Timed events (4, 30 and 60 minutes) enter as the pair (fixed time, distance rowed); distance events as (nominal distance, time). To look for departures from a straight line we add a quadratic term, $ln t = a_i + b thin x + c thin x^2$ with $x = ln(d "/" 2000)$. Then $b$ is the local exponent at 2 km and $2c$ is how fast the local exponent changes with $ln d$.

= Data

*Source.* Concept2 logbook rankings (log.concept2.com/rankings), RowErg only, season #key.config.seasons (the Concept2 season, not the calendar year). We used the "verified only" lists and kept individual results that carry an athlete identifier. Race-verified rows mostly lack an identifier and could not be linked reliably across events. Each ranking lists a person's best verified result for an event and season, so athletes appear once per event.

*Events.* The events at the extremes (1 minute, 100 m, half marathon, marathon, 100 km) were left out on purpose: they are dominated by anaerobic capacity, pacing and endurance limits that no single power law is expected to cover. The remaining events span roughly 500 m to 10 km (4 to 60 minutes), as in @events.

#figure(
  tbl((auto, auto, auto), ("Event", "Athlete-seasons (verified)", "In analysis"),
    events.map(r => (r.label, n(r.athlete_seasons_verified), n(r.athlete_seasons_in_analysis))),
    align: (left, right, right)),
  caption: [Events used. "In analysis" counts athlete-seasons that also meet the eligibility rule.],
) <events>

*Eligibility.* An athlete-season is analysed if it has at least #key.config.min_events events whose longest and shortest distances differ by a factor of at least #key.config.min_ratio, so that a slope is identifiable. @funnel shows how the data narrow down.

#figure(
  tbl((auto, auto, auto), ("Step", "Rows", "Athlete-seasons"),
    funnel.map(r => (r.step, n(r.rows), n(r.athlete_seasons))), align: (left, right, right)),
  caption: [Data funnel.],
) <funnel>

*Who is in the sample.* Only a minority of athlete-seasons qualify, and the sample is dominated by men and by masters rowers (@strata). To see whether qualifying athletes are unusual, @selection compares the median 2 km time of qualifying athletes with all verified athletes who have a 2 km result.

#figure(
  tbl((auto, auto, auto, auto), ("Group", "Sex", "n", "Median 2 km (s)"),
    selection.map(r => (r.group, r.sex, n(r.n), f(r.median_2k_time_s, d: 1))), align: (left, left, right, right)),
  caption: [Median 2 km time of qualifying athletes against all verified athletes with a 2 km result.],
) <selection>

= Methods

*Per-athlete exponents.* Each athlete-season gets its own least-squares slope. We report their median and mean, and the share below $beta$.

*Pooled within-athlete exponent.* The pooled slope is $sum S_(x y) "/" sum S_(x x)$ over athlete-seasons, where $S$ are within-athlete sums of squares and cross-products. It removes each athlete's overall level and weights athletes by the spread of their distances. Confidence intervals come from a cluster bootstrap (#n(key.config.n_boot) resamples) that resamples athletes, so an athlete with several seasons is not counted as independent.

*Distance dependence.* We estimate the quadratic model above, pooled and by group. Separately, for every event we compute the exponent implied between an athlete's #key.config.anchor_event m result (the anchor) and that event, $ln(t_e "/" t_0) "/" ln(d_e "/" d_0)$, together with the percentage by which the actual time differs from the 1.064 prediction. These two analyses do not need the same athletes to have three events.

*Event mix.* Because curvature makes the slope depend on which events an athlete rowed, we repeat the age and sex comparison within fixed distance windows (1–10 km and 2–21 km), using the same eligibility rule with a distance ratio of at least 2.

*Effort sensitivity.* Not every logged result is an all-out effort. We re-estimate with (i) a Huber-robust regression, (ii) dropping the 10% or 25% of athlete-seasons that fit worst after removing the common curvature, (iii) efforts up to 6 km only, and (iv) distance events only. These variants show how much the answer depends on unusual results; they cannot prove that all remaining results are maximal.

= Results

== Overall exponent

#let ov = row-of(strata, "group", "all")
@strata gives the pooled exponent for all athletes and by sex and age group. The pooled exponent for all athletes is #f(P.est) #ci(P.lo, P.hi); the median athlete slope is #f(key.overall.median_slope) with a standard deviation of #f(key.overall.sd_slope) across athletes, and #pct(100 * key.overall.share_below_beta0, d: 0) of athletes have a slope below #f(b0). @fig-slopes shows the distribution.

#figure(image("results/fig_slopes.png", width: 78%),
  caption: [Within-athlete exponents. Red: hypothesised 1.064; blue: median.]) <fig-slopes>

== Age and sex

#figure(
  tbl((auto, auto, auto, auto, auto, auto), ("Group", "Athletes", "Median", "Pooled", "95% CI", "Exponent at 2 km"),
    strata.map(r => (lbl(r.group), n(r.n_athletes), f(r.median_slope), f(r.pooled_slope),
      ci(r.pooled_lo, r.pooled_hi), f(get(cstrata, "group", r.group, "exp2k")))),
    align: (left, right, right, right, left, right)),
  caption: [Within-athlete exponent by sex and age group. Groups with fewer than 30 athletes have no confidence interval.],
) <strata>

#if age-grad [The exponent falls with age] else [There is no clear age gradient in the exponent]#if women-lower [ and is lower for women than for men] (@strata, @fig-age). Because older athletes tend to row longer events and curvature makes slopes depend on the event mix, we checked whether any difference could be an artefact. @agewin repeats the comparison inside fixed distance windows, and the last column of @strata gives the exponent at 2 km from the common-shape quadratic model.

#figure(image("results/fig_age.png", width: 80%),
  caption: [Pooled exponent by age group and sex (95% cluster-bootstrap intervals).]) <fig-age>

#figure(
  tbl((auto, auto, auto, auto, auto), ("Group", "All events", "1–10 km window", "2–21 km window", "Exponent at 2 km"),
    ("all", "F", "M", ..age-groups).map(g => (lbl(g),
      f(get(strata, "group", g, "pooled_slope")),
      f(aw("1 km - 10 km", g)),
      f(aw("2 km - 21 km", g)),
      f(get(cstrata, "group", g, "exp2k")))),
    align: (left, right, right, right, right)),
  caption: [Pooled exponent for the same athletes grouped by age, in fixed distance windows.],
) <agewin>

#if age-grad-win [Inside both windows the age gradient persists] else [Inside the windows the age gradient is much weaker or absent]: in the 2–21 km window the pooled exponent goes from #f(aw("2 km - 21 km", "<=18")) (under 19) to #f(aw("2 km - 21 km", "60+")) (60 and over). Women minus men is #f(aw("1 km - 10 km", "F") - aw("1 km - 10 km", "M")) in the 1–10 km window and #f(aw("2 km - 21 km", "F") - aw("2 km - 21 km", "M")) in the 2–21 km window. #if age-grad-win [So the age difference is not explained by event choice alone. We have not tested why it arises; training history, event specialisation and body mass (not available in the rankings) are all candidates.] else [So the apparent age difference in the overall estimate may largely reflect which events athletes row.]

== Distance dependence

The quadratic model gives a local exponent of #f(Cv.exponent_at_2km) at 2 km and a curvature of #f(Cv.curvature, d: 4) #ci(Cv.curvature_lo, Cv.curvature_hi, d: 4). Using that curved model instead of a constant 1.064 changes predicted times relative to the 2 km time as in @pred.

#figure(
  tbl((auto, auto, auto, auto), ("Distance", "Constant 1.064", "Curved model", "Difference"),
    (500, 1000, 5000, 6000, 10000).map(d => (str(d) + " m", f(ratio(b0, d), d: 3), f(curved(d), d: 3),
      pct(100 * (curved(d) / ratio(b0, d) - 1), d: 1))),
    align: (left, right, right, right)),
  caption: [Predicted time as a multiple of the 2 km time.],
) <pred>

@fig-implied and @anch show the exponent implied between the 2 km anchor and each other event. A value near 1.064 means that event is consistent with the hypothesis relative to the 2 km result.

#figure(image("results/fig_implied_exponent.png", width: 80%),
  caption: [Exponent implied between the 2 km result and each other event (median with 95% bootstrap interval); red line: 1.064.]) <fig-implied>

#figure(
  tbl((auto, auto, auto, auto, auto), ("Event", "n", "Implied exponent", "95% CI", "Median time vs 1.064 prediction"),
    anch.map(r => (r.label, n(r.n), f(r.implied_exponent), ci(r.implied_lo, r.implied_hi), pct(num(r.median_pct_dev), d: 1))),
    align: (left, right, right, left, right)),
  caption: [Implied exponent relative to the 2 km result. Positive percentages mean slower than predicted.],
) <anch>

The implied exponent ranges from #f(anch-sorted.first().implied_exponent) (#anch-sorted.first().label) to #f(anch-sorted.last().implied_exponent) (#anch-sorted.last().label), so how well the hypothesis holds depends on which event is compared with the 2 km result.

== Distance windows

#figure(
  tbl((auto, auto, auto, auto), ("Window (distance rowed)", "Athletes", "Pooled exponent", "95% CI"),
    windows.map(r => (lbl(r.group), n(r.n_athletes), f(r.pooled_slope), ci(r.pooled_lo, r.pooled_hi))),
    align: (left, right, right, left)),
  caption: [Pooled exponent in different distance windows.],
) <win>

The pooled exponent varies with the window by about #f(win-range) (@win)#if win-range > calc.abs(P.est - b0) [, which is more than the gap between the overall estimate and 1.064] else [, which is less than the gap between the overall estimate and 1.064]. The hypothesis is therefore best stated for a defined range of distances rather than as a universal constant.

== Effort sensitivity and other checks <sec-effort>

#figure(image("results/fig_effort.png", width: 100%),
  caption: [Linear exponent (left) and exponent at 2 km from the quadratic fit (right) under different treatments of unusual results. Red: 1.064.]) <fig-effort>

#figure(
  tbl((auto, auto, auto, auto, auto, auto), ("Variant", "Athletes", "Linear exponent", "95% CI", "Exponent at 2 km", "Curvature"),
    effort.map(r => (lbl(r.variant), n(r.n_athletes), f(r.linear), ci(r.linear_lo, r.linear_hi), f(r.exp2k), f(r.curv, d: 4))),
    align: (left, right, right, left, right, right)),
  caption: [Sensitivity of the estimates to unusual or sub-maximal results.],
) <effort>

The baseline and the alternatives differ by at most #f(max-dev) in the linear exponent (@effort). The residual scatter of an athlete-season around its fitted curve has a median of #pct(100 * E.resid_sd_median, d: 2) per result, a 90th percentile of #pct(100 * E.resid_sd_p90, d: 2) and a 99th percentile of #pct(100 * E.resid_sd_p99, d: 2). The skewness of the pooled residuals is #f(E.residual_skewness, d: 2)#if E.residual_skewness > 0.5 [, which points to a tail of unusually slow results] else [, so there is no strong tail of slow results].

Requiring more events per athlete-season (@sens) moves the pooled exponent by at most #f(sens-range)#if sens-range < 0.005 [, which suggests the main result is not driven by how athletes qualify] else [, so the result does depend on how athletes qualify].

#figure(
  tbl((auto, auto, auto, auto), ("Minimum events", "Athletes", "Pooled exponent", "95% CI"),
    sens.map(r => (lbl(r.group), n(r.n_athletes), f(r.pooled_slope), ci(r.pooled_lo, r.pooled_hi))),
    align: (left, right, right, left)),
  caption: [Pooled exponent when more events per athlete-season are required.],
) <sens>

#if "mixed" in key [
  A random-slope mixed model (#key.mixed.model) gives an overall exponent of #f(key.mixed.beta) #ci(key.mixed.lo, key.mixed.hi), in agreement with the pooled estimate. The estimated spread of true athlete slopes is #f(key.mixed.slope_sd).
]

= Discussion

*Sub-maximal efforts.* The hypothesis describes maximal efforts, but a logged verified result can be an easy piece. Because each ranking keeps an athlete's best verified result for an event and season, repeated easy sessions are mostly filtered out, but an athlete whose only verified 10 km is a steady row still contributes it. The sign of the resulting bias depends on where it occurs: sub-maximal results in long events raise the implied exponent; in short events they lower it; and a sub-maximal 2 km anchor shifts every implied exponent. These effects add noise (part of the athlete-to-athlete spread of slopes is probably effort rather than physiology) and could shift the pooled values in either direction. #if max-dev < 0.01 [The sensitivity analyses suggest the headline estimates are not driven by obviously unusual results, but they cannot identify results that are consistently easy.] else [The sensitivity analyses show the estimates do move when unusual results are down-weighted or removed, so the headline values should be treated cautiously.] The estimates should be read as exponents of typical verified efforts rather than of strictly maximal ones.

*Selection.* Athletes choose their events. Only about #pct(100 * share-elig, d: 0) of verified athlete-seasons have three or more events across a threefold range of distances, and the analysed sample is #pct(100 * share-men, d: 0) male and #pct(100 * share-over40, d: 0) aged 40 or over. Qualifying athletes can be compared with all verified athletes by their 2 km time (@selection), but they may still differ in motivation and training specificity. Results for women (#n(row-of(strata, "group", "F").n_athletes) athletes) and for athletes under 30 rest on much smaller groups.

#if low4 [*The 4-minute test.* The implied exponent for 4 minutes (#f(a4.implied_exponent)) is clearly lower than for the 1000 m event (#f(a1000.implied_exponent)), although the two efforts last nearly the same time. This suggests the people who record verified 4-minute results are a different group, or that timed and distance results are not equivalent; it deserves a separate look before being interpreted physiologically.]

*Range.* We excluded the extremes by design, so conclusions apply to roughly 500 m to 10 km (about 1.5 to 60 minutes), with the 60-minute event reaching about 14 km. The finding that the local exponent changes with distance means extrapolation to the marathon or to 100 m would not be safe.

*Other limitations.* One season is analysed. Body mass, drag factor and sleeping/training state are not available. Race-verified results were excluded because they cannot be linked to athletes. Age is as recorded in the ranking. Confidence intervals reflect sampling variation only, not model error.

= Conclusion

For verified RowErg results from #key.config.seasons the within-athlete exponent is #f(P.est) on average (#ci(P.lo, P.hi)), #if calc.abs(P.est - b0) < 0.01 [close to] else [clearly different from] the 1.064 used so far. Using 1.064 instead changes a predicted 6 km time by #pct(100 * (ratio(P.est, 6000) / ratio(b0, 6000) - 1), d: 1) relative to the 2 km time. #if traits.len() > 0 [But the exponent is not a universal constant: it #traits.join(", ", last: " and "). A prediction rule that uses the 2 km exponent (#f(Cv.exponent_at_2km)) with the curvature term, or one specific to age and sex, should do better than a single value; we did not test prediction accuracy directly.] Whether efforts are all maximal remains the main uncontrolled factor.

= Appendix: reproducibility

The report is produced by two scripts. `concept2_rankings.R` downloads and caches the ranking lists (with `status=verified` and `rower=rower`) into `data/rankings_long.csv`; `power_law_analysis.R` reads that file and writes `results/` (tables, `key_numbers.json`, figures). Settings: minimum #key.config.min_events events, distance ratio at least #key.config.min_ratio, anchor event #key.config.anchor_event m, #n(key.config.n_boot) bootstrap resamples, exponent under test #f(b0). Compile this report with `typst compile report.typ`.

#figure(image("results/fig_loglog.png", width: 70%),
  caption: [A random sample of athlete-seasons: log time against log distance, centred per athlete-season. The red line has the hypothesised slope.])
