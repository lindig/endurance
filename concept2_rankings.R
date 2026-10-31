#!/usr/bin/env Rscript
# concept2_rankings.R
#
# Download Concept2 logbook ranking lists (RowErg) into one long table, one row
# per ranking entry, with the athlete id taken from the row's detail link.
#
#   Rscript concept2_rankings.R inspect   # list the search form's field names/values
#   Rscript concept2_rankings.R           # download + parse + pilot checks
#
# Interactive: source("concept2_rankings.R"); main("inspect"); main(character())
#
# Needs R >= 4.1 and the packages httr2, rvest, xml2:
#   install.packages(c("httr2", "rvest", "xml2"))
#
# Raw HTML is cached in CACHE_DIR, so re-running resumes where it stopped and
# never re-downloads a page. Delete the cache (or set REFRESH <- TRUE) to refresh.
# Check the site's terms and robots.txt before running at scale.

suppressPackageStartupMessages({
  library(httr2)
  library(rvest)
  library(xml2)
})

## ---------------------------------------------------------------- config ----
SEASONS     <- 2026          # Concept2 "season" (not calendar year)
SEXES       <- "F"           # "F", "M" or c("F", "M")
EVENTS      <- c("1000", "2000")  # path codes; timed events: find codes with inspect
MACHINE     <- "rower"       # path segment for RowErg
MAX_PAGES   <- 2             # smoke test; set to Inf for the full pilot
EXTRA_QUERY <- list()        # server-side filters, e.g. list(verified = "1")
                             # (find the real names/values with `inspect`; if left
                             # empty everything is downloaded and the verified /
                             # type columns are used to filter afterwards)

SEASONS     <- 2026
SEXES       <- "F"
EVENTS      <- c("2000", "1", "4", "30", "60")
MAX_PAGES   <- 1
EXTRA_QUERY <- list(status = "verified", rower = "rower")

SEXES     <- c("F", "M")
EVENTS    <- c("4", "30", "60", "500", "1000", "2000", "5000", "6000", "10000")
MAX_PAGES <- 100

DELAY_S     <- 2             # pause before every real request ...
JITTER_S    <- 1             # ... plus random 0..JITTER_S seconds
REFRESH     <- FALSE         # TRUE = ignore the cache
KEEP_NAMES  <- FALSE         # FALSE = names are not written to the output CSV
CACHE_DIR   <- "cache/concept2"
OUT_DIR     <- "data"
USER_AGENT  <- "erg-power-law-research-script (personal research; contact: you@example.com)"
BASE        <- "https://log.concept2.com"

## ---------------------------------------------------------------- helpers ---
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a[1])) b else a

cache_path <- function(url) {
  key <- gsub("[^A-Za-z0-9]+", "_", sub("^https?://", "", url))
  file.path(CACHE_DIR, paste0(substr(key, 1, 200), ".html"))
}

# Returns the path of the cached HTML for `url`, downloading it if needed.
fetch_to_cache <- function(url) {
  path <- cache_path(url)
  if (!REFRESH && file.exists(path) && file.size(path) > 0) return(path)

  Sys.sleep(DELAY_S + runif(1, 0, JITTER_S))
  resp <- request(url) |>
    req_user_agent(USER_AGENT) |>
    req_retry(
      max_tries = 4,
      is_transient = \(r) resp_status(r) %in% c(429, 500, 502, 503, 504),
      backoff = \(i) 5 * 2^(i - 1)
    ) |>
    req_error(is_error = \(r) FALSE) |>
    req_perform()

  status <- resp_status(resp)
  if (status != 200) {
    stop(sprintf("HTTP %d for %s\nStopping; re-run later to resume from the cache.",
                 status, url), call. = FALSE)
  }
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  writeBin(resp_body_raw(resp), path)
  path
}

read_doc <- function(url) read_html(fetch_to_cache(url), encoding = "UTF-8")

list_url <- function(season, event, sex, page = 1) {
  query <- c(list(gender = sex), EXTRA_QUERY)
  if (page > 1) query$page <- page   # page 1 has no page parameter on the site
  req <- request(sprintf("%s/rankings/%s/%s/%s", BASE, season, MACHINE, event))
  req <- do.call(req_url_query, c(list(req), query))
  req$url
}

## ------------------------------------------------------------- inspect ------
# Prints every form control on the rankings search page, so you can see the
# parameter names/values behind "Verified", "Slides only", the timed events, etc.
inspect_form <- function() {
  doc <- read_doc(paste0(BASE, "/rankings"))
  forms <- xml_find_all(doc, "//form")
  for (i in seq_along(forms)) {
    f <- forms[[i]]
    cat(sprintf("\n== form %d  action=%s  method=%s\n", i,
                xml_attr(f, "action") %||% "(none)", xml_attr(f, "method") %||% "get"))
    for (s in xml_find_all(f, ".//select")) {
      opts <- xml_find_all(s, ".//option")
      cat(sprintf("  <select name=%s>  (%d options)\n", xml_attr(s, "name") %||% "?", length(opts)))
      show <- head(opts, 25)
      vals <- xml_attr(show, "value"); vals[is.na(vals)] <- ""
      cat(sprintf("      %-14s %s\n", vals, trimws(xml_text(show))), sep = "")
      if (length(opts) > 25) cat(sprintf("      ... %d more\n", length(opts) - 25))
    }
    for (inp in xml_find_all(f, ".//input")) {
      cat(sprintf("  <input name=%s type=%s value=%s>\n", xml_attr(inp, "name") %||% "?",
                  xml_attr(inp, "type") %||% "text", xml_attr(inp, "value") %||% ""))
    }
  }
  invisible(NULL)
}

## -------------------------------------------------------------- parsing -----
KNOWN <- c("Pos." = "pos", "Name" = "name", "Age" = "age", "Location" = "location",
           "Country" = "country", "Club/Affiliation" = "club",
           "Type" = "type", "Verified" = "verified")

n_pages <- function(doc) {
  hrefs <- html_attr(html_elements(doc, "a[href*='page=']"), "href")
  p <- suppressWarnings(as.integer(sub(".*[?&]page=([0-9]+).*", "\\1", hrefs)))
  if (length(p) == 0 || all(is.na(p))) 1L else max(p, na.rm = TRUE)
}

total_people <- function(doc) {
  txt <- html_text2(html_element(doc, "body"))
  m <- regmatches(txt, regexpr("[0-9,]+(?=\\s*Total people)", txt, perl = TRUE))
  if (length(m) == 0) NA_integer_ else as.integer(gsub(",", "", m))
}

parse_page <- function(doc, season, event, sex, page) {
  tables <- xml_find_all(doc, "//table")
  tbl <- NULL
  for (t in tables) {
    h <- trimws(xml_text(xml_find_all(t, ".//th")))
    if (all(c("Pos.", "Name") %in% h)) { tbl <- t; break }
  }
  if (is.null(tbl)) stop("No ranking table found (page layout changed?).", call. = FALSE)

  headers <- trimws(xml_text(xml_find_all(tbl, ".//th")))
  other <- !(headers %in% names(KNOWN))
  col_names <- ifelse(other, "result_raw", unname(KNOWN[headers]))
  col_names <- make.unique(col_names, sep = "_")
  result_header <- paste(headers[other], collapse = "|")

  rows <- xml_find_all(tbl, ".//tr[td]")
  out <- lapply(rows, function(r) {
    cells <- xml_find_all(r, "./td")
    if (length(cells) != length(headers)) return(NULL)
    txt <- trimws(xml_text(cells))
    names(txt) <- col_names
    href <- xml_attr(xml_find_first(cells[[match("name", col_names)]], ".//a"), "href")
    m <- regmatches(href, regexec("/rankings/detail/([0-9]+)/(individual|race)(?:/([0-9]+))?", href))[[1]]
    data.frame(
      as.list(txt),
      result_id   = if (length(m)) m[2] else NA_character_,
      result_kind = if (length(m)) m[3] else NA_character_,
      athlete_id  = if (length(m) && nzchar(m[4])) m[4] else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  skipped <- sum(vapply(out, is.null, logical(1)))
  if (skipped > 0) warning(sprintf("%d row(s) skipped on %s %s %s page %d (cell count mismatch)",
                                   skipped, season, event, sex, page))
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out) || nrow(out) == 0) return(NULL)
  cbind(data.frame(season = season, event = event, sex = sex, page = page,
                   result_header = result_header,
                   table_headers = paste(headers, collapse = "|"),
                   stringsAsFactors = FALSE), out)
}

# "3:09.4" / "1:23:45.6" -> seconds; "4,321" (or "4321 m") -> metres.
parse_result <- function(x) {
  x <- gsub("[,[:space:]]|m$", "", x)
  time_s <- rep(NA_real_, length(x)); dist_m <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    v <- x[i]
    if (is.na(v) || !nzchar(v)) next
    if (grepl(":", v)) {
      p <- suppressWarnings(as.numeric(strsplit(v, ":", fixed = TRUE)[[1]]))
      if (!anyNA(p)) time_s[i] <- sum(p * 60^((length(p) - 1):0))
    } else if (grepl("^[0-9]+$", v)) {
      dist_m[i] <- as.numeric(v)
    } else if (grepl("^[0-9]+\\.[0-9]+$", v)) {
      time_s[i] <- as.numeric(v)
    }
  }
  data.frame(time_s = time_s, distance_m = dist_m)
}

## ------------------------------------------------------------ download ------
download_list <- function(season, event, sex) {
  doc1 <- read_doc(list_url(season, event, sex, 1))
  np <- n_pages(doc1)
  expected <- total_people(doc1)
  last <- min(np, MAX_PAGES)
  message(sprintf("[%s %s %s] %d page(s), %s people; fetching %d",
                  season, event, sex, np, format(expected), last))
  pages <- vector("list", last)
  for (p in seq_len(last)) {
    doc <- if (p == 1) doc1 else read_doc(list_url(season, event, sex, p))
    pages[[p]] <- parse_page(doc, season, event, sex, p)
    if (p %% 10 == 0) message(sprintf("   page %d/%d", p, last))
  }
  res <- do.call(rbind, pages)
  res$filter_status <- EXTRA_QUERY$status %||% "all"
  res$filter_rower  <- EXTRA_QUERY$rower  %||% "all"
  if (last == np && !is.na(expected) && nrow(res) != expected) {
    warning(sprintf("[%s %s %s] parsed %d rows but page says %d people",
                    season, event, sex, nrow(res), expected))
  }
  res
}

## ------------------------------------------------------------- checks -------
pilot_report <- function(df) {
  cat("\n---- rows per event / sex ----\n"); print(table(df$event, df$sex))
  cat("\n---- table headers and result column (format check) ----\n")
  for (ev in unique(df$event)) {
    d <- df[df$event == ev, ]
    cat(sprintf("event %s: headers = %s\n   result header = '%s'; examples: %s\n", ev,
                d$table_headers[1], d$result_header[1],
                paste(head(d$result_raw, 3), collapse = ", ")))
  }
  has_vt <- all(c("verified", "type") %in% names(df))
  if (has_vt) {
    cat("\n---- verified x type ----\n"); print(table(df$verified, df$type, useNA = "ifany"))
    cat("\n---- share of rows with an athlete id, by verified ----\n")
    print(round(tapply(!is.na(df$athlete_id), df$verified, mean), 3))
  } else {
    cat("\n(no 'verified' and/or 'type' column in these tables; the server-side filters",
        "\n seem to drop them. Filters requested: status =", unique(df$filter_status),
        ", rower =", unique(df$filter_rower), ")\n")
    cat(sprintf("share of rows with an athlete id: %.3f\n", mean(!is.na(df$athlete_id))))
  }

  ok <- df[!is.na(df$athlete_id), ]
  if (has_vt) ok <- ok[ok$verified == "Yes", ]
  if (nrow(ok) > 0) {
    if ("name" %in% names(ok)) {
      nm <- tapply(ok$name, ok$athlete_id, function(z) length(unique(z)))
      cat(sprintf("\nathlete ids with more than one distinct name: %d of %d\n", sum(nm > 1), length(nm)))
    }
    dup <- aggregate(pos ~ athlete_id + event + season + sex, ok, length)
    cat(sprintf("athlete/event/season combinations listed more than once: %d\n", sum(dup$pos > 1)))
    ev <- unique(ok$event)
    if (length(ev) >= 2) {
      sets <- lapply(ev, \(e) unique(ok$athlete_id[ok$event == e]))
      cat(sprintf("athletes present in all of events %s: %d\n",
                  paste(ev, collapse = ", "), length(Reduce(intersect, sets))))
      n_ev <- table(ok$athlete_id[!duplicated(paste(ok$athlete_id, ok$event))])
      cat("number of events per athlete:\n"); print(table(n_ev))
    }
  }
  invisible(NULL)
}

## --------------------------------------------------------------- main -------
main <- function(args = character()) {
  if ("inspect" %in% args) { inspect_form(); return(invisible(NULL)) }
  if (length(EXTRA_QUERY) == 0) {
    message("Note: EXTRA_QUERY is empty, so unverified results are downloaded too; ",
            "filter on `verified` / `type` afterwards.")
  }
  parts <- list()
  for (season in SEASONS) for (sex in SEXES) for (event in EVENTS) {
    parts[[length(parts) + 1]] <- download_list(season, event, sex)
  }
  df <- do.call(rbind, parts)
  df <- cbind(df, parse_result(df$result_raw))

  pilot_report(df)

  if (!KEEP_NAMES) df$name <- NULL
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(OUT_DIR, "rankings_long.csv")
  write.csv(df, out, row.names = FALSE, fileEncoding = "UTF-8")
  message(sprintf("\nWrote %d rows to %s", nrow(df), out))
  invisible(df)
}

if (!interactive()) main(commandArgs(trailingOnly = TRUE))
