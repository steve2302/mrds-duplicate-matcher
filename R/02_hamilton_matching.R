# ============================================================================
# Hamilton et al. (2018) duplicate-matching options: app subset
# ----------------------------------------------------------------------------
# Sections 0-3 and 7 of scripts/hamilton_options.R in the analysis project,
# copied unchanged: fixed-tolerance scoring, the discrepancy mixture model,
# stochastic cluster resolution and the match_duplicates() dispatcher.
# Group-size imputation, MRDS fitting, bootstrap and MI are left out because
# the app only does matching (and leaving them out keeps the mrds package out
# of the WebAssembly bundle).
#
# Hamilton, O. N. P., Kincaid, S. E., Constantine, R., Kozmian-Ledward, L.,
# Walker, C. G. & Fewster, R. M. (2018) Accounting for uncertainty in
# duplicate identification and group size judgments in mark-recapture
# distance sampling. Methods in Ecology and Evolution 9: 354-362.
#
# Needs 01_mrds_duplicate_match.R to be sourced first (Shiny sources R/ in
# alphabetical order).
# ============================================================================

# ---------------------------------------------------------------------------
# 0. Helpers
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Spatial disagreement between the two records of each pair. Hamilton et al.
# used the difference in downward sighting angle (degrees). Fallbacks are the
# difference in perpendicular distance and then the GPS separation d_m.
.space_diff <- function(pairs, space_field = c("auto", "angle", "distance", "d_m")) {
  space_field <- match.arg(space_field)
  if (space_field == "auto") {
    space_field <- if (all(c("angle_A", "angle_B") %in% names(pairs))) "angle"
                   else if (all(c("distance_A", "distance_B") %in% names(pairs))) "distance"
                   else "d_m"
  }
  v <- switch(space_field,
              angle    = abs(pairs$angle_A - pairs$angle_B),
              distance = abs(pairs$distance_A - pairs$distance_B),
              d_m      = pairs$d_m)
  attr(v, "field") <- space_field
  v
}

# Same-observer pairs within the gate. Assuming an observer reports each
# group once, these pairs cannot be duplicates. They serve as "known
# non-matches" anchoring the far component of the discrepancy mixture, which
# is the role opposite-side pairs play in Hamilton et al.
build_known_nonmatch_pairs <- function(detections,
                                       max_time_gap_s  = 30,
                                       max_space_gap_m = 500,
                                       type = c("same_observer", "opposite_side")) {
  type <- match.arg(type)
  if (type == "opposite_side") {
    if (!"side" %in% names(detections)) stop("opposite_side needs a `side` column.")
    A <- detections |> filter(observer == "A")
    B <- detections |> filter(observer == "B")
    by <- intersect("transect", names(detections))
    pairs <- inner_join(A, B, by = by, suffix = c("_A", "_B"),
                        relationship = "many-to-many") |>
      filter(side_A != side_B)
  } else {
    d <- detections
    by <- intersect(c("transect", "observer"), names(d))
    pairs <- inner_join(d, d, by = by, suffix = c("_A", "_B"),
                        relationship = "many-to-many") |>
      filter(id_A < id_B)
  }
  dt <- if (inherits(pairs$time_A, "POSIXct"))
          as.numeric(difftime(pairs$time_B, pairs$time_A, units = "secs"))
        else pairs$time_B - pairs$time_A
  pairs |>
    mutate(dt = dt,
           d_m = geosphere::distHaversine(cbind(lon_A, lat_A), cbind(lon_B, lat_B))) |>
    filter(abs(dt) <= max_time_gap_s, d_m <= max_space_gap_m)
}

# ---------------------------------------------------------------------------
# 1. Fixed-tolerance scoring (Hamilton "fixedTol = TRUE")
# ---------------------------------------------------------------------------
# A pair is a duplicate if |dt - time_offset| <= time_tol and the spatial
# disagreement <= space_tol. Hamilton et al. used T = 5-10 s and A = 10-15
# degrees of sighting angle. space_tol is in the units of `space_field`.

score_fixed_tolerance <- function(pairs, time_tol = 10, space_tol = 15,
                                  space_field = "auto", time_offset = 0) {
  sd <- .space_diff(pairs, space_field)
  pairs$p_match <- as.numeric(abs(pairs$dt - time_offset) <= time_tol &
                              sd <= space_tol)
  list(pairs  = pairs,
       params = list(time_tol = time_tol, space_tol = space_tol,
                     space_field = attr(sd, "field"), time_offset = time_offset),
       method = "fixed")
}

# ---------------------------------------------------------------------------
# 2. Discrepancy mixture model (Hamilton discrepancy.func)
# ---------------------------------------------------------------------------
# Discrepancy  r = min(|dt|, time_cap) / time_scale + space_diff / space_scale
# (Hamilton: time_cap = time_scale = 25 s, space_scale = 65 degrees).
# r is modelled as a two-component mixture of Normals truncated to
# [0, discrep_max]: a "near" (duplicate) and a "far" (non-duplicate)
# component. Known non-matches contribute only to the far component.
# P(duplicate | r) is the near component's posterior share. Pairs with
# r > discrep_max get P = 0, as in Hamilton.
#
# space_scale defaults to 65 for angles and otherwise to the spatial gate, so
# the spatial term spans roughly [0, 1] like the time term.
#
# The scales set how sharply duplicates separate from non-duplicates and
# should be tuned to the survey. Hamilton's 25 s / 65 degrees suit their
# dolphin survey. On the tandem simulation in test_hamilton_options.R,
# time_scale = 5, time_cap = 15, centre_time = TRUE, space_scale = 50 (m)
# raised recall from 0.75 to 0.92 under stochastic_cluster.

fit_discrepancy_mixture <- function(pairs,
                                    time_cap    = 25,
                                    time_scale  = 25,
                                    space_scale = NULL,
                                    space_field = "auto",
                                    discrep_max = 1,
                                    centre_time = FALSE,
                                    known_nonmatch = NULL,
                                    start_cut   = 0.4,
                                    verbose     = FALSE) {

  discrep_of <- function(p, mu_t = 0) {
    sd <- .space_diff(p, space_field)
    sc <- space_scale %||%
      if (attr(sd, "field") == "angle") 65 else (attr(pairs, "max_space_gap_m") %||% max(sd))
    pmin(abs(p$dt - mu_t), time_cap) / time_scale + sd / sc
  }
  # Hamilton uses |dt|; centre_time = TRUE first removes the median lag of
  # close pairs, which helps when the rear observer calls systematically late.
  mu_t <- if (centre_time) median(pairs$dt[pairs$d_m <= quantile(pairs$d_m, 0.25)]) else 0

  r_all <- discrep_of(pairs, mu_t)
  keep  <- r_all <= discrep_max
  r     <- r_all[keep]
  r_kn  <- if (!is.null(known_nonmatch) && nrow(known_nonmatch)) {
             x <- discrep_of(known_nonmatch, mu_t); x[x <= discrep_max]
           } else numeric(0)
  if (length(r) < 5) stop("Too few candidate pairs within discrep_max to fit the mixture.")

  ldtn <- function(x, m, s)   # log truncated-normal density on [0, discrep_max]
    dnorm(x, m, s, log = TRUE) -
      log(pmax(pnorm(discrep_max, m, s) - pnorm(0, m, s), 1e-300))

  nll <- function(th) {
    mN <- th[1]; sN <- th[2]; mF <- th[3]; sF <- th[4]; pD <- th[5]
    a <- log(pD)     + ldtn(r, mN, sN)
    b <- log(1 - pD) + ldtn(r, mF, sF)
    mx <- pmax(a, b)
    -(sum(mx + log(exp(a - mx) + exp(b - mx))) + sum(ldtn(r_kn, mF, sF)))
  }

  # Start values from a fixed split, as Hamilton did with T = 10, A = 15
  near <- r[r <= start_cut]; far <- c(r[r > start_cut], r_kn)
  if (length(near) < 3 || length(far) < 3) {
    near <- r[r <= quantile(r, 0.3)]; far <- r[r > quantile(r, 0.3)]
  }
  st <- c(mean(near), max(sd(near), 0.02), mean(far), max(sd(far), 0.05),
          min(max(length(near) / length(r), 0.01), 0.99))
  fit <- optim(st, nll, method = "L-BFGS-B",
               lower = c(0, 1e-3, 0, 1e-3, 1e-4),
               upper = c(5 * discrep_max, 5 * discrep_max,
                         5 * discrep_max, 5 * discrep_max, 1 - 1e-4),
               control = list(maxit = 2000))
  if (fit$convergence != 0) {   # restart once from the stopping point
    fit2 <- optim(fit$par, nll, method = "L-BFGS-B", lower = c(0, 1e-3, 0, 1e-3, 1e-4),
                  upper = c(5 * discrep_max, 5 * discrep_max,
                            5 * discrep_max, 5 * discrep_max, 1 - 1e-4),
                  control = list(maxit = 2000))
    if (fit2$value <= fit$value) fit <- fit2
  }
  if (verbose && fit$convergence != 0) message("optim: ", fit$message)
  th <- fit$par
  if (th[1] > th[3]) th <- c(th[3], th[4], th[1], th[2], 1 - th[5])   # label switch
  names(th) <- c("mu1", "sig1", "mu2", "sig2", "dupP")
  if (verbose) print(round(th, 4))

  a <- log(th[5])     + ldtn(r_all, th[1], th[2])
  b <- log(1 - th[5]) + ldtn(r_all, th[3], th[4])
  p <- 1 / (1 + exp(b - a))
  p[!keep] <- 0
  pairs$discrep <- r_all
  pairs$p_match <- p

  list(pairs = pairs, params = as.list(th), loglik = -fit$value,
       converged = fit$convergence == 0, optim_code = fit$convergence,
       optim_message = fit$message, method = "discrepancy",
       settings = list(time_cap = time_cap, time_scale = time_scale,
                       space_scale = space_scale, discrep_max = discrep_max,
                       centre_time = centre_time, mu_t = mu_t,
                       n_known_nonmatch = length(r_kn)))
}

# ---------------------------------------------------------------------------
# 3. Stochastic selection and cluster resolution (Hamilton
#    compile.data.func + build.dup.data.func + dup.cluster.func)
# ---------------------------------------------------------------------------
# 1. Each candidate pair is declared a duplicate with probability p_match.
# 2. Selected pairs form connected clusters. While a cluster has more than
#    one pair: draw an anchor sighting with probability proportional to the
#    summed p_match of the pairs it is in, merge every sighting linked to
#    the anchor into one record, and drop all pairs touching the anchor or
#    the merged sightings.
# Merging lets one observer's single sighting match two sightings by the
# other observer (a group recorded as split), which one-to-one assignment
# cannot express. Each call returns one random realisation.

resolve_clusters_stochastic <- function(fit, select = TRUE) {
  pairs <- fit$pairs |> filter(p_match > 0)
  if (select) pairs <- pairs[runif(nrow(pairs)) < pairs$p_match, , drop = FALSE]
  sel <- tibble(a = paste0("A:", pairs$id_A), b = paste0("B:", pairs$id_B),
                p = pairs$p_match)

  out <- list(); status <- list()
  while (nrow(sel)) {
    cl <- .components(sel$a, sel$b)
    sizes <- table(cl)
    multi <- names(sizes)[sizes > 1]
    if (!length(multi)) break
    rows   <- which(cl == multi[1])
    inds   <- c(sel$a[rows], sel$b[rows])
    anchor <- inds[sample.int(length(inds), 1, prob = rep(sel$p[rows], 2))]
    linked <- rows[sel$a[rows] == anchor | sel$b[rows] == anchor]
    merged <- setdiff(c(sel$a[linked], sel$b[linked]), anchor)
    out[[length(out) + 1]] <- list(anchor = anchor, merged = merged,
                                   p = max(sel$p[linked]))
    gone <- c(anchor, merged)
    sel <- sel[!(sel$a %in% gone | sel$b %in% gone), , drop = FALSE]
  }

  strip <- function(x) sub("^[AB]:", "", x)
  side  <- function(x) substr(x, 1, 1)
  res_merged <- bind_rows(lapply(out, function(o) {
    ids <- c(o$anchor, o$merged)
    tibble(id_A = paste(strip(ids[side(ids) == "A"]), collapse = "+"),
           id_B = paste(strip(ids[side(ids) == "B"]), collapse = "+"),
           p_match = o$p, resolution = "merged")
  }))
  res_simple <- tibble(id_A = strip(sel$a), id_B = strip(sel$b),
                       p_match = sel$p, resolution = rep("pair", nrow(sel)))
  bind_rows(res_simple, res_merged)
}

# Connected components of an undirected graph given as an edge list, by
# minimum-label propagation. Returns the component label of each edge.
.components <- function(from, to) {
  nodes <- unique(c(from, to))
  f <- match(from, nodes); t <- match(to, nodes)
  lab <- seq_along(nodes)
  repeat {
    m <- pmin(lab[f], lab[t])
    best <- tapply(c(m, m), c(f, t), min)          # smallest label per node
    idx <- as.integer(names(best))
    new <- lab
    new[idx] <- pmin(lab[idx], best)
    if (identical(new, lab)) break
    lab <- new
  }
  lab[f]
}

# ---------------------------------------------------------------------------
# 7. Dispatcher
# ---------------------------------------------------------------------------

match_duplicates <- function(detections,
                             method     = c("fellegi_sunter", "discrepancy", "fixed"),
                             assignment = c("hungarian", "stochastic_cluster", "mh_sample"),
                             max_time_gap_s  = 30,
                             max_space_gap_m = 500,
                             method_args = list(),
                             known_nonmatch = c("none", "same_observer", "opposite_side"),
                             threshold = 0.5, K = 1, verbose = FALSE) {
  method <- match.arg(method); assignment <- match.arg(assignment)
  known_nonmatch <- match.arg(known_nonmatch)
  pairs <- build_candidate_pairs(detections, max_time_gap_s, max_space_gap_m)
  n_A <- sum(detections$observer == "A"); n_B <- sum(detections$observer == "B")

  fit <- switch(method,
    fellegi_sunter = do.call(fit_em, c(list(pairs = pairs, n_A = n_A, n_B = n_B,
                                            verbose = verbose), method_args)),
    discrepancy = {
      kn <- if (known_nonmatch == "none") NULL
            else build_known_nonmatch_pairs(detections, max_time_gap_s,
                                            max_space_gap_m, known_nonmatch)
      do.call(fit_discrepancy_mixture, c(list(pairs = pairs, known_nonmatch = kn,
                                              verbose = verbose), method_args))
    },
    fixed = do.call(score_fixed_tolerance, c(list(pairs = pairs), method_args)))

  draws <- switch(assignment,
    hungarian          = list(assign_matches(fit, threshold = threshold)),
    stochastic_cluster = replicate(K, resolve_clusters_stochastic(fit), simplify = FALSE),
    mh_sample          = sample_matchings(fit, K = K))
  list(fit = fit, matches = draws[[1]], draws = draws,
       method = method, assignment = assignment)
}

