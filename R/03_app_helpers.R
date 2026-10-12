# ============================================================================
# App helpers: example data, column mapping, time parsing, and a wrapper that
# runs the matcher. The statistical work is all done by the functions in
# 01_mrds_duplicate_match.R and 02_hamilton_matching.R.
# ============================================================================

# Refuse jobs whose within-transect A x B table would exhaust browser memory.
MAX_CARTESIAN_PAIRS <- 5e6

# ---------------------------------------------------------------------------
# Example data
# ---------------------------------------------------------------------------
# A simulated tandem-observer survey with known truth, adapted from
# simulate_survey() in scripts/test_hamilton_options.R: one side watched,
# front (A) and rear (B) observers with independent half-normal detection,
# the rear observer calling about 1.5 s later.

simulate_example <- function(n_tr = 6, L_m = 30000, n_groups = 120, w = 300,
                             speed = 50, mu_lag = 1.5, sd_lag = 0.8,
                             g0 = c(A = 0.9, B = 0.8), sigma = 110,
                             species_pool = c("EKR", "WGK", "EUR", "WAL"),
                             species_prob = c(0.5, 0.3, 0.15, 0.05),
                             species_correct = 0.98, seed = 2026) {
  if (!is.null(seed)) set.seed(seed)
  m_per_deg_lon <- 111320 * cos(-30 * pi / 180)
  truth <- bind_rows(lapply(seq_len(n_tr), function(k) tibble(
    transect = sprintf("T%02d", k),
    group_id = sprintf("T%02d_G%03d", k, seq_len(n_groups)),
    along = runif(n_groups, 0, L_m), perp = runif(n_groups, 0, w),
    species = sample(species_pool, n_groups, replace = TRUE, prob = species_prob),
    size  = 1 + rpois(n_groups, 3),
    # transects are 0.02 degrees of latitude (~2.2 km) apart
    lat0  = -30 + (k - 1) * 0.02)))
  rec <- function(obs) {
    p <- g0[[obs]] * exp(-truth$perp^2 / (2 * sigma^2))
    g <- truth[runif(nrow(truth)) < p, ]
    n <- nrow(g)
    lag  <- if (obs == "A") rnorm(n, 0, 0.5) else rnorm(n, mu_lag, sd_lag)
    perp <- pmax(0, g$perp + rnorm(n, 0, 4))
    sz   <- pmax(1L, as.integer(round(g$size * exp(rnorm(n, 0, 0.15)))))
    wrong <- runif(n) > species_correct
    sp <- g$species
    sp[wrong] <- vapply(sp[wrong], function(s) sample(setdiff(species_pool, s), 1), "")
    tibble(id = paste0(obs, "_", g$group_id), observer = obs,
           transect = g$transect,
           time = round(g$along / speed + lag, 2),
           lon = round(145 + (g$along + rnorm(n, 0, 5)) / m_per_deg_lon, 6),
           lat = round(g$lat0 + (perp + rnorm(n, 0, 5)) / 111320, 6),
           distance = round(perp, 1), group_size = sz, species = sp,
           truth_id = g$group_id)
  }
  bind_rows(rec("A"), rec("B")) |> arrange(transect, time)
}

# ---------------------------------------------------------------------------
# Column mapping
# ---------------------------------------------------------------------------

# role -> (label, required, name patterns used to guess the column)
COLUMN_ROLES <- list(
  observer   = list(label = "Observer", required = TRUE,
                    pat = "^obs(erver)?(_?(id|code|name|pos(ition)?))?$|^seat$|^position$"),
  time       = list(label = "Time of sighting", required = TRUE,
                    pat = "^(date_?)?time$|^secs?$|^seconds$|^timestamp$|^t$|time|clock|^utc$"),
  lon        = list(label = "Longitude (decimal degrees)", required = TRUE,
                    pat = "^lon(g|gitude)?$|^x$"),
  lat        = list(label = "Latitude (decimal degrees)", required = TRUE,
                    pat = "^lat(itude)?$|^y$"),
  group_size = list(label = "Group size", required = TRUE,
                    pat = "^group_?size$|^size$|^count$|^n(umber)?$|^cluster_?size$|size"),
  id         = list(label = "Detection ID", required = FALSE,
                    pat = "^id$|^det(ection)?_?id$|^record(_?id)?$|^object$"),
  transect   = list(label = "Transect", required = FALSE,
                    pat = "^transect(_?(id|label|no|number))?$|^line$|^sample\\.?label$"),
  side       = list(label = "Side of aircraft", required = FALSE,
                    pat = "^side$"),
  species    = list(label = "Species", required = FALSE,
                    pat = "^species$|^spp?$|^taxon$"),
  distance   = list(label = "Perpendicular distance", required = FALSE,
                    pat = "^(perp(endicular)?_?)?dist(ance)?$"),
  angle      = list(label = "Sighting angle (degrees)", required = FALSE,
                    pat = "^angle$|^declination$")
)

NONE <- "(none)"

guess_column <- function(role, cols) {
  hit <- cols[grepl(COLUMN_ROLES[[role]]$pat, cols, ignore.case = TRUE)]
  exact <- hit[tolower(hit) == role]
  if (length(exact)) exact[1] else if (length(hit)) hit[1] else NONE
}

# ---------------------------------------------------------------------------
# Time parsing: everything becomes numeric seconds
# ---------------------------------------------------------------------------

parse_time <- function(x, format = c("auto", "seconds", "clock", "datetime")) {
  format <- match.arg(format)
  clock_re <- "^\\s*\\d{1,2}:\\d{2}(:\\d{2}(\\.\\d+)?)?\\s*$"
  if (format == "auto") {
    format <- if (is.numeric(x)) "seconds"
              else if (all(grepl(clock_re, x[!is.na(x)]))) "clock"
              else "datetime"
  }
  out <- switch(format,
    seconds = suppressWarnings(as.numeric(x)),
    clock = vapply(strsplit(trimws(as.character(x)), ":", fixed = TRUE), function(p) {
      p <- suppressWarnings(as.numeric(p))
      if (length(p) < 2 || anyNA(p)) return(NA_real_)
      sum(p * c(3600, 60, 1)[seq_along(p)])
    }, 0),
    datetime = {
      # One format is chosen for the whole column (the one that reads most
      # values), so day-first and year-first dates are never mixed.
      fmts <- c("%Y-%m-%d %H:%M:%OS", "%Y/%m/%d %H:%M:%OS",
                "%d/%m/%Y %H:%M:%OS", "%d-%m-%Y %H:%M:%OS")
      xc <- sub("Z$", "", trimws(as.character(x)))
      xc <- sub("^(\\d{4}-\\d{2}-\\d{2})T", "\\1 ", xc)
      best <- rep(NA_real_, length(xc))
      for (f in fmts) {
        v <- as.numeric(as.POSIXct(strptime(xc, f, tz = "UTC")))
        if (sum(!is.na(v)) > sum(!is.na(best))) best <- v
      }
      best
    })
  attr(out, "format") <- format
  out
}

# ---------------------------------------------------------------------------
# Standardise an uploaded table to the columns the matcher expects
# ---------------------------------------------------------------------------
# `map` is a named list role -> column name (or NONE). Returns
# list(data = tibble, notes = character()). Stops with a readable message on
# anything the matcher could not handle.

standardise_detections <- function(raw, map, obs_A, obs_B, time_format = "auto") {
  notes <- character()
  pick <- function(role) {
    v <- map[[role]]
    if (is.null(v) || identical(v, NONE) || !nzchar(v)) NULL else v
  }
  for (role in names(COLUMN_ROLES)) if (COLUMN_ROLES[[role]]$required && is.null(pick(role)))
    stop("Choose a column for: ", COLUMN_ROLES[[role]]$label, call. = FALSE)
  used <- unlist(lapply(names(COLUMN_ROLES), pick))
  if (anyDuplicated(used))
    stop("The same column is mapped to two fields: ",
         paste(unique(used[duplicated(used)]), collapse = ", "), call. = FALSE)
  if (is.null(obs_A) || is.null(obs_B) || !nzchar(obs_A) || !nzchar(obs_B) || obs_A == obs_B)
    stop("Choose two different observer values for A (front) and B (rear).", call. = FALSE)

  obs_raw <- as.character(raw[[pick("observer")]])
  keep <- !is.na(obs_raw) & obs_raw %in% c(obs_A, obs_B)
  if (sum(!keep))
    notes <- c(notes, sprintf("%d rows dropped: observer is not '%s' or '%s'.",
                              sum(!keep), obs_A, obs_B))
  raw <- raw[keep, , drop = FALSE]
  obs_raw <- obs_raw[keep]
  if (!nrow(raw)) stop("No rows left after selecting the two observers.", call. = FALSE)

  num <- function(role) suppressWarnings(as.numeric(raw[[pick(role)]]))
  tm <- parse_time(raw[[pick("time")]], time_format)
  if (mean(is.na(tm)) > 0.5)
    stop("Could not read the time column as ", attr(tm, "format"),
         ". Check the time format setting.", call. = FALSE)

  out <- tibble(
    observer   = ifelse(obs_raw == obs_A, "A", "B"),
    time       = as.numeric(tm),
    lon        = num("lon"),
    lat        = num("lat"),
    group_size = num("group_size"))
  if (any(abs(out$lon) > 180, na.rm = TRUE) || any(abs(out$lat) > 90, na.rm = TRUE))
    stop("Longitude/latitude are outside the range of decimal degrees. ",
         "Projected coordinates (eastings/northings) are not supported.", call. = FALSE)

  if (!is.null(pick("transect"))) out$transect <- as.character(raw[[pick("transect")]])
  if (!is.null(pick("side")))     out$side     <- as.character(raw[[pick("side")]])
  # A constant species is exactly neutral: every pair agrees under both the
  # match and the non-match model.
  out$species <- if (!is.null(pick("species"))) as.character(raw[[pick("species")]])
                 else "unknown"
  if (is.null(pick("species")))
    notes <- c(notes, "No species column: species is not used to separate matches.")
  if (!is.null(pick("distance"))) out$distance <- num("distance")
  if (!is.null(pick("angle")))    out$angle    <- num("angle")

  if (!is.null(pick("id"))) {
    id <- as.character(raw[[pick("id")]])
    if (anyNA(id) || anyDuplicated(id))
      stop("Detection IDs must be unique and non-missing. Set Detection ID to ",
           NONE, " to have them generated.", call. = FALSE)
    if (any(grepl("+", id, fixed = TRUE)))
      stop("Detection IDs must not contain '+', which marks merged sightings.", call. = FALSE)
  } else {
    id <- character(nrow(out))
    for (o in c("A", "B")) {
      i <- which(out$observer == o)
      id[i] <- sprintf("%s%05d", o, seq_along(i))
    }
    notes <- c(notes, "Detection IDs generated as A00001..., B00001... in file order.")
  }
  out <- bind_cols(tibble(id = id), out)

  # Carry the user's other columns through, renaming any that would clash
  # with columns the matcher creates.
  reserved <- c(names(out), "dt", "d_m", "log_size_ratio", "species_agree",
                "p_match", "log_LR", "discrep", "seen", "group_id", "size",
                "id_A", "id_B", "min", "best", "max")
  extra <- raw[, setdiff(names(raw), used), drop = FALSE]
  clash <- names(extra) %in% reserved
  names(extra)[clash] <- paste0(names(extra)[clash], "_orig")
  names(extra) <- make.unique(names(extra))
  out <- bind_cols(out, as_tibble(extra))

  bad <- is.na(out$time) | is.na(out$lon) | is.na(out$lat)
  if (any(bad))
    notes <- c(notes, sprintf(paste("%d detections have a missing time or position.",
               "They cannot be matched and are kept as single-observer records."), sum(bad)))
  bad_g <- is.na(out$group_size) | out$group_size <= 0
  if (any(bad_g))
    notes <- c(notes, sprintf(paste("%d detections have a missing or non-positive group size.",
               "Group size is treated as uninformative for those."), sum(bad_g)))
  if (!all(c("A", "B") %in% out$observer))
    stop("Both observers need at least one detection.", call. = FALSE)
  if (is.null(pick("transect")))
    notes <- c(notes, paste("No transect column: every A detection is compared with every",
                            "B detection, which is slow for large files."))
  list(data = out, notes = notes)
}

# ---------------------------------------------------------------------------
# Run the matcher
# ---------------------------------------------------------------------------
# Same steps as match_duplicates(), with two differences that do not change
# the model: (1) assignment is solved transect by transect when a transect
# column exists, which is equivalent because candidate pairs never cross
# transects, and keeps the Hungarian cost matrix small enough for a browser;
# (2) the burn-in and thinning of the MH sampler are exposed.

n_cartesian <- function(det) {
  if ("transect" %in% names(det)) {
    tb <- table(det$transect, det$observer)
    if (!all(c("A", "B") %in% colnames(tb))) return(0)
    sum(as.numeric(tb[, "A"]) * as.numeric(tb[, "B"]))
  } else as.numeric(sum(det$observer == "A")) * sum(det$observer == "B")
}

run_matching <- function(det,
                         method = "fellegi_sunter", assignment = "hungarian",
                         max_time_gap_s = 30, max_space_gap_m = 500,
                         method_args = list(), known_nonmatch = "none",
                         threshold = 0.5, K = 1, n_burn = 200, n_thin = 10,
                         seed = NULL, progress = function(msg) NULL) {
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  nc <- n_cartesian(det)
  if (nc > MAX_CARTESIAN_PAIRS)
    stop(sprintf(paste("%.1f million A x B comparisons would be needed, more than a",
                       "browser can hold. Map a transect column, or split the file."),
                 nc / 1e6), call. = FALSE)

  sf <- method_args$space_field
  if (!is.null(sf) && sf %in% c("angle", "distance") && !sf %in% names(det))
    stop("Spatial field '", sf, "' needs a mapped ", sf, " column.", call. = FALSE)
  if (known_nonmatch == "opposite_side" && !"side" %in% names(det))
    stop("Opposite-side anchors need a mapped side column.", call. = FALSE)

  progress("Building candidate pairs")
  pairs <- build_candidate_pairs(det, max_time_gap_s, max_space_gap_m)
  if (!nrow(pairs))
    stop("No candidate pairs fall inside the gate. Widen the time or distance gate, ",
         "and check the time format and coordinates.", call. = FALSE)
  n_A <- sum(det$observer == "A"); n_B <- sum(det$observer == "B")

  progress("Scoring pairs")
  fit <- switch(method,
    fellegi_sunter = do.call(fit_em, c(list(pairs = pairs, n_A = n_A, n_B = n_B,
                                            verbose = FALSE), method_args)),
    discrepancy = {
      kn <- if (known_nonmatch == "none") NULL
            else build_known_nonmatch_pairs(det, max_time_gap_s, max_space_gap_m,
                                            known_nonmatch)
      do.call(fit_discrepancy_mixture, c(list(pairs = pairs, known_nonmatch = kn),
                                         method_args))
    },
    fixed = do.call(score_fixed_tolerance, c(list(pairs = pairs), method_args)),
    stop("Unknown scoring method: ", method))
  fit$method <- method

  progress("Assigning duplicates")
  by_tr <- "transect" %in% names(fit$pairs)
  draws <- switch(assignment,
    hungarian = list(assign_matches(fit, threshold = threshold,
                                    transects = if (by_tr) TRUE else NULL)),
    stochastic_cluster = replicate(K, resolve_clusters_stochastic(fit), simplify = FALSE),
    mh_sample = {
      if (by_tr) {
        per_tr <- lapply(split(fit$pairs, fit$pairs$transect), function(p)
          sample_matchings(list(pairs = p), K = K, n_burn = n_burn, n_thin = n_thin,
                           init_threshold = threshold))
        lapply(seq_len(K), function(k) bind_rows(lapply(per_tr, `[[`, k)))
      } else sample_matchings(fit, K = K, n_burn = n_burn, n_thin = n_thin,
                              init_threshold = threshold)
    },
    stop("Unknown assignment method: ", assignment))

  list(detections = det, fit = fit, draws = draws, method = method,
       assignment = assignment,
       settings = list(max_time_gap_s = max_time_gap_s, max_space_gap_m = max_space_gap_m,
                       threshold = threshold, K = length(draws), seed = seed))
}

# One draw as a display table: pair fields joined on for simple pairs; merged
# sightings ("+"-joined ids) keep only their ids and probability.
match_table <- function(res, k = 1) {
  m <- res$draws[[k]]
  cols <- intersect(c("id_A", "id_B", "transect", "dt", "d_m", "group_size_A",
                      "group_size_B", "species_A", "species_B", "p_match"),
                    names(res$fit$pairs))
  out <- tibble(id_A = as.character(m$id_A), id_B = as.character(m$id_B)) |>
    left_join(res$fit$pairs[, cols], by = c("id_A", "id_B"))
  if ("resolution" %in% names(m)) {
    out$resolution <- m$resolution
    out$p_match <- ifelse(is.na(out$p_match), m$p_match, out$p_match)
  }
  out
}

# How often each pair (or merged set) is declared a duplicate across draws.
pair_frequency <- function(res) {
  K <- length(res$draws)
  all <- bind_rows(lapply(res$draws, function(m)
    tibble(id_A = as.character(m$id_A), id_B = as.character(m$id_B))))
  all |> count(id_A, id_B, name = "n_draws") |>
    mutate(proportion = n_draws / K) |>
    left_join(res$fit$pairs[, c("id_A", "id_B", "p_match")], by = c("id_A", "id_B")) |>
    arrange(desc(proportion), id_A)
}

all_draws_long <- function(res) {
  bind_rows(lapply(seq_along(res$draws), function(k) {
    m <- res$draws[[k]]
    tibble(draw = k, id_A = as.character(m$id_A), id_B = as.character(m$id_B))
  }))
}

# Precision and recall against a truth_id column (example data only).
evaluate_truth <- function(det, matches) {
  if (!"truth_id" %in% names(det)) return(NULL)
  truth <- setNames(as.character(det$truth_id), det$id)
  n_true <- sum(table(det$truth_id[det$observer %in% c("A", "B")]) == 2)
  ids_a <- strsplit(as.character(matches$id_A), "+", fixed = TRUE)
  ids_b <- strsplit(as.character(matches$id_B), "+", fixed = TRUE)
  tp <- mapply(function(a, b) any(truth[a] %in% truth[b]), ids_a, ids_b)
  tp <- as.logical(unlist(tp))
  list(n_true = n_true, n_declared = nrow(matches), true_pos = sum(tp),
       precision = if (length(tp)) mean(tp) else NA_real_,
       recall = if (n_true) sum(tp) / n_true else NA_real_)
}

PARAM_LABELS <- c(
  mu_t = "Mean time offset, B minus A (s)", sigma_t = "SD of time offset (s)",
  sigma_s = "Position scale, Rayleigh (m)", sigma_g = "SD of log group-size ratio",
  eps = "Species disagreement rate", pi_match = "Share of candidate pairs that match",
  mu1 = "Near (duplicate) component mean", sig1 = "Near component SD",
  mu2 = "Far (non-duplicate) component mean", sig2 = "Far component SD",
  dupP = "Near component weight",
  time_tol = "Time tolerance (s)", space_tol = "Spatial tolerance",
  space_field = "Spatial field", time_offset = "Time offset (s)")

param_table <- function(fit) {
  p <- fit$params
  tibble(Parameter = unname(ifelse(names(p) %in% names(PARAM_LABELS),
                                   PARAM_LABELS[names(p)], names(p))),
         Estimate = vapply(p, function(v) if (is.numeric(v)) format(signif(v, 4)) else
                             as.character(v), ""))
}
