# ============================================================================
# Checks for the app's wrappers. Base R only; run from the app root:
#   Rscript tests/test_matching.R
# ============================================================================

for (f in sort(list.files("R", full.names = TRUE, pattern = "\\.R$"))) source(f)

n_fail <- 0
check <- function(ok, what) {
  cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", what))
  if (!isTRUE(ok)) n_fail <<- n_fail + 1
}
same_pairs <- function(a, b) {
  key <- function(m) sort(paste(m$id_A, m$id_B))
  identical(key(a), key(b))
}

det <- simulate_example()
check(all(c("A", "B") %in% det$observer) && !anyDuplicated(det$id), "example data is well formed")

# ---- 1. Wrapper reproduces the plain function calls ------------------------
set.seed(1)
pairs <- build_candidate_pairs(det, 15, 300)
fit   <- fit_em(pairs, sum(det$observer == "A"), sum(det$observer == "B"), verbose = FALSE)
direct <- assign_matches(fit, threshold = 0.5)            # one global assignment
res <- run_matching(det, "fellegi_sunter", "hungarian", 15, 300, seed = 1)
check(isTRUE(all.equal(fit$params, res$fit$params)), "run_matching: same EM estimates as fit_em()")
check(same_pairs(direct, res$draws[[1]]),
      "per-transect Hungarian assignment equals the global assignment")

set.seed(1)
md <- match_duplicates(det, max_time_gap_s = 15, max_space_gap_m = 300)
check(same_pairs(md$matches, res$draws[[1]]), "run_matching agrees with match_duplicates()")

ev <- evaluate_truth(det, res$draws[[1]])
cat(sprintf("    FS + Hungarian: %d declared, precision %.3f, recall %.3f\n",
            ev$n_declared, ev$precision, ev$recall))
check(ev$precision > 0.95 && ev$recall > 0.95, "example data: precision and recall above 0.95")

# ---- 2. Every scoring x assignment combination -----------------------------
margs <- list(fellegi_sunter = list(),
              discrepancy = list(time_scale = 5, time_cap = 15, centre_time = TRUE,
                                 space_scale = 50, space_field = "d_m"),
              fixed = list(time_tol = 4, space_tol = 40, space_field = "d_m",
                           time_offset = 1.5))
for (m in names(margs)) for (a in c("hungarian", "stochastic_cluster", "mh_sample")) {
  r <- tryCatch(run_matching(det, m, a, 15, 300, method_args = margs[[m]],
                             K = if (a == "hungarian") 1 else 5,
                             n_burn = 500, n_thin = 20, seed = 2),
                error = function(e) conditionMessage(e))
  ok <- is.list(r)
  if (ok) {
    h  <- build_history(det, r$draws[[length(r$draws)]])
    mt <- match_table(r, 1)
    e  <- evaluate_truth(det, r$draws[[1]])
    n_det <- sum(lengths(strsplit(c(na.omit(h$id_A), na.omit(h$id_B)), "+", fixed = TRUE)))
    ok <- all(h$seen %in% c("11", "10", "01")) && n_det == nrow(det) &&
          nrow(mt) == nrow(r$draws[[1]]) && length(r$draws) == r$settings$K
    cat(sprintf("    %-15s %-19s n=%3d precision %.3f recall %.3f\n", m, a,
                e$n_declared, e$precision, e$recall))
    if (a != "hungarian") ok <- ok && nrow(pair_frequency(r)) > 0 && nrow(all_draws_long(r)) > 0
  } else cat("    error:", r, "\n")
  check(ok, sprintf("%s + %s runs and every detection appears once in the history", m, a))
}

# ---- 3. Column mapping ------------------------------------------------------
hms <- function(s) sprintf("%02d:%02d:%05.2f", s %/% 3600, (s %% 3600) %/% 60, s %% 60)
raw <- data.frame(Seat = ifelse(det$observer == "A", "front", "rear"),
                  Clock = hms(det$time + 9 * 3600), Longitude = det$lon,
                  Latitude = det$lat, Count = det$group_size, Line = det$transect,
                  Spp = det$species, PerpDist = det$distance, truth_id = det$truth_id,
                  stringsAsFactors = FALSE)
raw <- rbind(raw, transform(raw[1, ], Seat = "pilot"))        # a row to be dropped
check(guess_column("observer", names(raw)) == "Seat" &&
      guess_column("lon", names(raw)) == "Longitude" &&
      guess_column("transect", names(raw)) == "Line" &&
      guess_column("id", names(raw)) == NONE, "column guessing")
map <- list(observer = "Seat", time = "Clock", lon = "Longitude", lat = "Latitude",
            group_size = "Count", id = NONE, transect = "Line", side = NONE,
            species = "Spp", distance = "PerpDist", angle = NONE)
s <- standardise_detections(raw, map, "front", "rear")
check(nrow(s$data) == nrow(det) && any(grepl("1 rows dropped", s$notes)),
      "rows for other observers are dropped and reported")
check(isTRUE(all.equal(s$data$time - 9 * 3600, det$time, tolerance = 1e-6)),
      "clock times are parsed to seconds")
r2 <- run_matching(s$data, "fellegi_sunter", "hungarian", 15, 300, seed = 1)
truth_of <- function(d, m) sort(paste(d$truth_id[match(m$id_A, d$id)], d$truth_id[match(m$id_B, d$id)]))
check(identical(truth_of(s$data, r2$draws[[1]]), truth_of(det, res$draws[[1]])),
      "mapped upload gives the same matches as the original table")

check(isTRUE(all.equal(as.numeric(parse_time(c("2025-06-01 09:15:32.5", "2025-06-01T09:15:34"))),
                       c(0, 1.5) + as.numeric(as.POSIXct("2025-06-01 09:15:32.5", tz = "UTC")))),
      "date-times are parsed, with fractional seconds")
check(identical(attr(parse_time(c(1.5, 2)), "format"), "seconds"), "numeric time is taken as seconds")

# Missing optional fields and missing values
map2 <- modifyList(map, list(species = NONE, transect = NONE, distance = NONE))
raw3 <- raw[raw$Line == "T01" & raw$Seat != "pilot", ]
raw3$Count[3] <- NA; raw3$Longitude[5] <- NA
s3 <- standardise_detections(raw3, map2, "front", "rear")
r3 <- tryCatch(run_matching(s3$data, "fellegi_sunter", "hungarian", 15, 300, seed = 1),
               error = function(e) conditionMessage(e))
check(is.list(r3) && nrow(build_history(s3$data, r3$draws[[1]])) > 0 &&
        all(is.finite(r3$fit$pairs$p_match)),
      "runs without species or transect columns and with missing values")

# Guards
err <- function(expr) tryCatch({ expr; "" }, error = function(e) conditionMessage(e))
check(grepl("unique", err(standardise_detections(raw, modifyList(map, list(id = "Line", transect = NONE)),
                                                "front", "rear"))), "duplicate IDs are rejected")
check(grepl("gate", err(run_matching(det, max_time_gap_s = 0.001, max_space_gap_m = 0.001))),
      "an empty gate gives a readable error")
check(grepl("angle", err(run_matching(det, "fixed", method_args = list(space_field = "angle")))),
      "asking for a missing spatial field gives a readable error")
big <- det[rep(1:nrow(det), 12), ]; big$transect <- NULL
check(grepl("million", err(run_matching(big))), "oversized jobs are refused")

cat(sprintf("\n%s\n", if (n_fail) paste(n_fail, "check(s) FAILED") else "All checks passed"))
if (n_fail) quit(status = 1)
