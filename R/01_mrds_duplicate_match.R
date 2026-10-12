# ============================================================================
# Probabilistic duplicate matching for tandem-observer MRDS
# ----------------------------------------------------------------------------
# Implements a Fellegi-Sunter record linkage model with EM estimation of
# match/non-match field distributions, followed by one-to-one assignment
# via the Hungarian algorithm. Designed for aerial surveys with two
# independent observers on the same side.
#
# Required matching fields per detection:
#   id, observer (A/B), transect, time (POSIXct or numeric seconds),
#   lon, lat (WGS84 decimal degrees), group_size (integer), species (factor/char)
#
# Author: drafted for Steven McLeod, VPRU/DPIRD
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(geosphere)   # haversine for lon/lat -> metres
  library(clue)        # solve_LSAP for Hungarian matching
})

# ---------------------------------------------------------------------------
# 1. Candidate pair generation
# ---------------------------------------------------------------------------
# Only build pairs that are physically plausible. A hard time/space gate
# avoids O(n_A * n_B) blow-up and stops EM being dominated by obviously
# non-matching pairs. The gate must be GENEROUS — set to several times the
# expected match-class SD so very few true matches are excluded.
#
# Returns a long tibble of candidate (a, b) pairs with comparison fields.

build_candidate_pairs <- function(detections,
                                  max_time_gap_s   = 30,
                                  max_space_gap_m  = 500,
                                  by_transect      = TRUE) {

  stopifnot(c("id", "observer", "time", "lon", "lat",
              "group_size", "species") %in% names(detections))

  A <- detections |> filter(observer == "A")
  B <- detections |> filter(observer == "B")

  # Cartesian join within transect (or globally if by_transect = FALSE).
  # If a `side` column is present, only same-side pairs can be duplicates.
  join_by_cols <- intersect(c(if (by_transect) "transect", "side"),
                            names(detections))
  if (length(join_by_cols)) {
    pairs <- inner_join(A, B, by = join_by_cols,
                        suffix = c("_A", "_B"),
                        relationship = "many-to-many")
  } else {
    pairs <- A |>
      mutate(.k = 1) |>
      inner_join(B |> mutate(.k = 1), by = ".k",
                 suffix = c("_A", "_B"),
                 relationship = "many-to-many") |>
      select(-.k)
  }

  # Convert time to numeric seconds if POSIXct
  if (inherits(pairs$time_A, "POSIXct")) {
    pairs <- pairs |>
      mutate(dt = as.numeric(difftime(time_B, time_A, units = "secs")))
  } else {
    pairs <- pairs |> mutate(dt = time_B - time_A)
  }

  # Great-circle distance in metres.
  # Missing or non-positive group sizes give a neutral size ratio (0) and
  # missing species count as agreement, so neither field can veto a match.
  pairs <- pairs |>
    mutate(
      d_m = distHaversine(cbind(lon_A, lat_A),
                          cbind(lon_B, lat_B)),
      log_size_ratio = ifelse(is.finite(group_size_A) & group_size_A > 0 &
                              is.finite(group_size_B) & group_size_B > 0,
                              log(group_size_B / group_size_A), 0),
      species_agree  = as.integer(is.na(species_A) | is.na(species_B) |
                                  species_A == species_B)
    ) |>
    filter(abs(dt) <= max_time_gap_s,
           d_m   <= max_space_gap_m)

  # Record the gate so fit_em() can normalise the non-match densities over
  # the same region the pairs were drawn from.
  attr(pairs, "max_time_gap_s")  <- max_time_gap_s
  attr(pairs, "max_space_gap_m") <- max_space_gap_m
  pairs
}

# ---------------------------------------------------------------------------
# 2. EM for the Fellegi-Sunter model
# ---------------------------------------------------------------------------
# Match-class densities:
#   dt              ~ N(mu_t, sigma_t^2)
#   d_m             ~ Rayleigh-like via 2D Gaussian on (dx, dy);
#                     we use d_m directly with chi(2) scaling for tractability
#   log_size_ratio  ~ N(0, sigma_g^2)   (mean fixed at 0 by symmetry)
#   species_agree   ~ Bernoulli(1 - eps)
#
# Non-match densities are uniform/empirical (computed once from the data).
# We work on the LOG scale throughout for numerical stability.

fit_em <- function(pairs,
                   n_A, n_B,
                   # initial values
                   init = list(mu_t = 0, sigma_t = 3,
                               sigma_s = 50, sigma_g = 0.3,
                               eps = 0.02, pi_match = NULL),
                   # search window for non-match densities; defaults to the
                   # gate used in build_candidate_pairs() so u is normalised
                   # over the region the pairs were actually drawn from
                   time_window_s   = attr(pairs, "max_time_gap_s"),
                   space_window_m  = attr(pairs, "max_space_gap_m"),
                   max_iter = 200,
                   tol = 1e-6,
                   verbose = TRUE) {

  if (is.null(time_window_s) || is.null(space_window_m))
    stop("time_window_s / space_window_m not supplied and pairs carry no gate attributes.")

  # --- non-match (u) densities, fixed -------------------------------------
  # u_t: uniform on [-time_window, +time_window]
  log_u_t <- -log(2 * time_window_s)
  # u_s: distance under uniform 2D over a disc of radius space_window
  #      f(d) = 2d / R^2, log_u_s = log(2*d/R^2)
  R <- space_window_m
  # u_g: null distribution of log_size_ratio, from pairing A's sizes with a
  #      random permutation of B's sizes. Uses the caller's RNG stream.
  ok_g <- is.finite(pairs$group_size_A) & pairs$group_size_A > 0 &
          is.finite(pairs$group_size_B) & pairs$group_size_B > 0
  shuffled <- log(sample(pairs$group_size_B[ok_g]) / pairs$group_size_A[ok_g])
  if (length(shuffled) >= 2) {
    u_g_mu <- mean(shuffled); u_g_sd <- max(sd(shuffled), 0.05)
  } else {
    u_g_mu <- 0; u_g_sd <- 1
  }
  # u_sp: probability random pair shares species
  sp_freq <- table(c(pairs$species_A, pairs$species_B))
  sp_p    <- sp_freq / sum(sp_freq)
  u_sp_agree <- sum(sp_p^2)

  # --- initialise match-class parameters ----------------------------------
  mu_t    <- init$mu_t
  sigma_t <- init$sigma_t
  sigma_s <- init$sigma_s
  sigma_g <- init$sigma_g
  eps     <- init$eps
  if (is.null(init$pi_match)) {
    # rough guess: 30% of the smaller list are duplicates
    n_pairs <- nrow(pairs)
    pi_match <- min(n_A, n_B) * 0.3 / n_pairs
  } else {
    pi_match <- init$pi_match
  }
  pi_match <- max(min(pi_match, 0.5), 1e-4)

  ll_prev <- -Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # --- E-step: match-class log densities -------------------------------
    log_m_t  <- dnorm(pairs$dt,  mu_t,  sigma_t, log = TRUE)
    # Spatial: 2D isotropic Gaussian, marginal density of distance is
    #   f(d) = d/sigma_s^2 * exp(-d^2/(2 sigma_s^2))
    log_m_s  <- log(pairs$d_m + 1e-9) - 2*log(sigma_s) -
                pairs$d_m^2 / (2 * sigma_s^2)
    log_m_g  <- dnorm(pairs$log_size_ratio, 0, sigma_g, log = TRUE)
    log_m_sp <- ifelse(pairs$species_agree == 1, log(1 - eps), log(eps + 1e-9))

    # Non-match: spatial Rayleigh on uniform disc
    log_u_s  <- log(2 * pairs$d_m / R^2 + 1e-12)
    log_u_g  <- dnorm(pairs$log_size_ratio, u_g_mu, u_g_sd, log = TRUE)
    log_u_sp <- ifelse(pairs$species_agree == 1, log(u_sp_agree),
                                                 log(1 - u_sp_agree))

    log_lik_M <- log_m_t + log_m_s + log_m_g + log_m_sp
    log_lik_U <- log_u_t + log_u_s + log_u_g + log_u_sp

    # Cap extreme log-LR to avoid degeneracy. Capping the ratio (rather than
    # M and U separately) keeps the ordering of pairs intact.
    log_LR <- pmin(pmax(log_lik_M - log_lik_U, -50), 50)

    # Posterior match probability
    log_pi   <- log(pi_match)
    log_1mpi <- log(1 - pi_match)
    a <- log_pi   + log_lik_U + log_LR
    b <- log_1mpi + log_lik_U
    mx <- pmax(a, b)
    p  <- exp(a - mx) / (exp(a - mx) + exp(b - mx))

    # observed-data log-likelihood (for monitoring)
    ll <- sum(mx + log(exp(a - mx) + exp(b - mx)))

    # --- M-step ----------------------------------------------------------
    Sp     <- sum(p)
    if (Sp < 1) {
      if (verbose) message("EM: expected matches < 1 — stopping early.")
      break
    }
    mu_t    <- sum(p * pairs$dt) / Sp
    sigma_t <- sqrt(max(sum(p * (pairs$dt - mu_t)^2) / Sp, 0.25))   # floor 0.5 s
    sigma_s <- sqrt(max(sum(p * pairs$d_m^2) / (2 * Sp), 25))       # floor 5 m
    sigma_g <- sqrt(max(sum(p * pairs$log_size_ratio^2) / Sp, 0.01))
    eps     <- 1 - sum(p * pairs$species_agree) / Sp
    eps     <- min(max(eps, 1e-4), 0.5)
    pi_match <- Sp / nrow(pairs)
    pi_match <- min(max(pi_match, 1e-4), 0.5)

    if (verbose && (iter %% 10 == 1 || iter <= 5))
      cat(sprintf("iter %3d  ll=%.3f  pi=%.4f  mu_t=%.2f  sig_t=%.2f  sig_s=%.1f  sig_g=%.2f  eps=%.3f\n",
                  iter, ll, pi_match, mu_t, sigma_t, sigma_s, sigma_g, eps))

    if (abs(ll - ll_prev) < tol * max(1, abs(ll))) {
      converged <- TRUE
      break
    }
    ll_prev <- ll
  }

  pairs$p_match <- p
  pairs$log_LR  <- log_LR

  # --- diagnostic warnings ---
  warnings_emitted <- character()
  if (pi_match >= 0.49)
    warnings_emitted <- c(warnings_emitted,
      "pi_match hit upper cap (0.5): match/non-match classes may not be separable. Check gate width and field informativeness.")
  if (pi_match <= 2e-4)
    warnings_emitted <- c(warnings_emitted,
      "pi_match collapsed to near zero: EM found no match class. Try better initial values or check for systematic offset (mu_t).")
  if (eps >= 0.3)
    warnings_emitted <- c(warnings_emitted,
      "Species disagreement among inferred matches is high (eps > 0.3): species field is uninformative or matches contain many false positives.")
  if (sigma_t >= time_window_s / 3)
    warnings_emitted <- c(warnings_emitted,
      "sigma_t is large relative to time gate: time field provides little discrimination.")
  if (length(warnings_emitted) && verbose) {
    cat("\n--- EM diagnostics ---\n")
    for (w in warnings_emitted) cat("WARNING: ", w, "\n", sep = "")
  }

  list(pairs   = pairs,
       params  = list(mu_t = mu_t, sigma_t = sigma_t,
                      sigma_s = sigma_s, sigma_g = sigma_g,
                      eps = eps, pi_match = pi_match),
       loglik  = ll,
       n_iter  = iter,
       converged = converged,
       warnings  = warnings_emitted)
}

# ---------------------------------------------------------------------------
# 3. One-to-one assignment (Hungarian)
# ---------------------------------------------------------------------------
# Each detection in A matches at most one detection in B. We maximise the
# total log-posterior-odds across matched pairs, subject to a minimum
# per-pair probability (threshold).

assign_matches <- function(em_fit,
                           threshold = 0.5,
                           transects = NULL) {

  pairs <- em_fit$pairs
  if (!is.null(transects) && "transect" %in% names(pairs)) {
    out_list <- lapply(split(pairs, pairs$transect),
                       assign_matches_one,
                       threshold = threshold)
    return(bind_rows(out_list))
  }
  assign_matches_one(pairs, threshold)
}

assign_matches_one <- function(pairs, threshold = 0.5) {
  if (nrow(pairs) == 0) return(pairs[0, ])
  ids_A <- unique(pairs$id_A)
  ids_B <- unique(pairs$id_B)
  nA <- length(ids_A); nB <- length(ids_B)

  # Augmented (nA + nB) square cost matrix so every detection always has its
  # own "unmatched" slot, however dense the candidate pairs are:
  #
  #             B_1..B_nB            dummy_A_1..dummy_A_nA
  #   A_i   [ -logodds / BIG   |   -thr_logit on diagonal, BIG off ]
  #   dB_j  [ -thr_logit diag  |   0                               ]
  #
  # solve_LSAP minimises. Leaving A_i and B_j both unmatched costs
  # 2 * (-thr_logit / 2) = -thr_logit, so the pair (i, j) is preferred only
  # when its log-odds exceed thr_logit, i.e. p_match > threshold.
  thr_logit <- log(threshold / (1 - threshold))
  BIG <- 1e6
  N <- nA + nB
  cost <- matrix(BIG, N, N)
  for (k in seq_len(nrow(pairs))) {
    i <- match(pairs$id_A[k], ids_A)
    j <- match(pairs$id_B[k], ids_B)
    p <- min(max(pairs$p_match[k], 1e-9), 1 - 1e-9)
    cost[i, j] <- -log(p / (1 - p))   # negative log-odds (to minimise)
  }
  half <- -thr_logit / 2
  cost[cbind(seq_len(nA), nB + seq_len(nA))] <- half         # A_i unmatched
  cost[cbind(nA + seq_len(nB), seq_len(nB))] <- half         # B_j unmatched
  cost[nA + seq_len(nB), nB + seq_len(nA)]   <- 0            # dummy-dummy
  sol <- as.integer(clue::solve_LSAP(cost - min(cost)))

  col_A <- sol[seq_len(nA)]
  real  <- col_A <= nB
  matched <- tibble(id_A = ids_A[real], id_B = ids_B[col_A[real]]) |>
    left_join(pairs, by = c("id_A", "id_B")) |>
    filter(!is.na(p_match), p_match >= threshold)
  matched
}


# ---------------------------------------------------------------------------
# 4. Build final detection history
# ---------------------------------------------------------------------------
# Classifies every detection as: 11 (seen by both), 10 (A only), 01 (B only)
# and collapses to one row per physical group.
#
# `matches` may contain merged sightings from resolve_clusters_stochastic(),
# where id_A or id_B is several ids joined by "+" (one observer split a group
# that the other saw as one). Merged sides follow Hamilton et al.: group-size
# fields are summed (NA only if all NA) and distance is a size-weighted mean.
#
# Output keeps the primary detection's columns (A's record for a duplicate)
# plus:
#   id_A, id_B            contributing ids ("+"-joined if merged)
#   group_size_A/_B       per-observer size (summed if merged)
#   size                  consensus size: rounded mean of the two observers
#   distance_A/_B, distance   per-observer and consensus (mean) distance,
#                             if detections carry a `distance` column
#   min_A..max_B          per-observer min/best/max, if present, for
#                         impute_group_size()

build_history <- function(detections, matches) {
  split_ids <- function(x) strsplit(as.character(x), "+", fixed = TRUE)
  matched_A <- unlist(split_ids(matches$id_A))
  matched_B <- unlist(split_ids(matches$id_B))
  out <- detections |>
    mutate(seen = case_when(
      observer == "A" & id %in% matched_A ~ "11",
      observer == "B" & id %in% matched_B ~ "11",
      observer == "A"                      ~ "10",
      observer == "B"                      ~ "01"
    ))

  size_fields <- intersect(c("group_size", "min", "best", "max"), names(detections))
  has_dist    <- "distance" %in% names(detections)

  # Summarise one side (A or B) of each matched group, merging if needed.
  side_summary <- function(ids_col, suffix) {
    long <- tibble(.row = seq_along(ids_col), id = split_ids(ids_col)) |>
      tidyr::unnest(id) |>
      left_join(detections, by = "id")
    long$.wt <- .size_weight(best = long[["best"]], min = long[["min"]],
                             max = long[["max"]], gs = long[["group_size"]],
                             n = nrow(long))
    sum_na <- function(v) if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
    agg <- long |>
      group_by(.row) |>
      summarise(across(all_of(size_fields), sum_na),
                distance = if (has_dist) sum(distance * .wt) / sum(.wt) else NA_real_,
                .groups = "drop")
    if (!has_dist) agg$distance <- NULL
    names(agg)[-1] <- paste0(names(agg)[-1], "_", suffix)
    agg
  }

  pair_rows <- matches |>
    select(id_A, id_B) |>
    mutate(.row = row_number(), group_id = paste0("g", .row))
  if (nrow(pair_rows)) {
    pair_rows <- pair_rows |>
      left_join(side_summary(pair_rows$id_A, "A"), by = ".row") |>
      left_join(side_summary(pair_rows$id_B, "B"), by = ".row")
  }
  # Primary record = first A id of each matched group
  first_A <- vapply(split_ids(pair_rows$id_A), `[`, "", 1)
  out_pairs <- out |>
    filter(seen == "11", observer == "A", id %in% first_A) |>
    select(-any_of(c(size_fields, if (has_dist) "distance"))) |>
    left_join(pair_rows |> mutate(.first = first_A) |> select(-.row),
              by = c("id" = ".first"))

  # Singletons: their own fields become the _A or _B columns
  out_single <- out |>
    filter(seen != "11") |>
    mutate(group_id = paste0("s", row_number()),
           id_A = if_else(observer == "A", as.character(id), NA_character_),
           id_B = if_else(observer == "B", as.character(id), NA_character_))
  for (f in c(size_fields, if (has_dist) "distance")) {
    out_single[[paste0(f, "_A")]] <- ifelse(out_single$observer == "A", out_single[[f]], NA)
    out_single[[paste0(f, "_B")]] <- ifelse(out_single$observer == "B", out_single[[f]], NA)
  }

  res <- bind_rows(out_pairs, out_single) |> arrange(group_id)

  # Consensus fields
  mean2 <- function(a, b) ifelse(is.na(a), b, ifelse(is.na(b), a, (a + b) / 2))
  if ("group_size" %in% size_fields) {
    res$group_size <- ifelse(is.na(res$group_size_A), res$group_size_B, res$group_size_A)
    res$size <- round(mean2(res$group_size_A, res$group_size_B))
  }
  if (has_dist) res$distance <- mean2(res$distance_A, res$distance_B)
  res
}

# Group-size weights for merging distances (Hamilton et al.): best if
# present, else mean(min, max), else group_size, else 1.
.size_weight <- function(best = NULL, min = NULL, max = NULL, gs = NULL, n) {
  w <- rep(NA_real_, n)
  if (!is.null(best)) w <- best
  if (!is.null(min) || !is.null(max)) {
    mm <- rowMeans(cbind(if (is.null(min)) NA else min,
                         if (is.null(max)) NA else max), na.rm = TRUE)
    w <- ifelse(is.na(w), mm, w)
  }
  if (!is.null(gs)) w <- ifelse(is.na(w), gs, w)
  w[is.na(w) | w <= 0] <- 1
  w
}

# ---------------------------------------------------------------------------
# 5. Diagnostics
# ---------------------------------------------------------------------------

diag_residual_dependence <- function(em_fit) {
  # Residual within-class correlations: should be ~0 if conditional
  # independence holds. Computed as weighted correlations using p_match
  # as weights for the match class and (1-p_match) for non-match.
  p <- em_fit$pairs$p_match
  vars <- with(em_fit$pairs,
               cbind(dt = dt, d_m = d_m,
                     lsr = log_size_ratio,
                     sp  = species_agree))
  wcor <- function(w) {
    wm <- colSums(vars * w) / sum(w)
    centred <- sweep(vars, 2, wm)
    v <- (t(centred) %*% (centred * w)) / sum(w)
    sds <- sqrt(diag(v))
    r <- v / (sds %o% sds)
    r[is.nan(r)] <- NA   # zero-variance field (e.g. species always agrees)
    r
  }
  list(match_class     = wcor(p),
       nonmatch_class  = wcor(1 - p))
}

diag_classification_table <- function(em_fit, breaks = c(0, 0.1, 0.5, 0.9, 1)) {
  p <- em_fit$pairs$p_match
  table(cut(p, breaks = breaks, include.lowest = TRUE))
}

# ---------------------------------------------------------------------------
# 6. Multiple-imputation sampling of matchings
# ---------------------------------------------------------------------------
# Draw K matchings from the posterior over partial bipartite assignments.
# Uses a simple Metropolis sampler over the matching: propose a swap of two
# match assignments (or add/drop a match), accept with the usual MH ratio.
# Output: list of K data frames each giving one matching.
#
# Use these in MRDS as a multiple-imputation: fit your detection function
# once per matching and combine using Rubin's rules. This propagates
# matching uncertainty into your abundance estimates.

sample_matchings <- function(em_fit, K = 50, n_burn = 200, n_thin = 10,
                             init_threshold = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pairs <- em_fit$pairs
  ids_A <- unique(pairs$id_A); ids_B <- unique(pairs$id_B)
  nA <- length(ids_A); nB <- length(ids_B)

  # Build log-odds matrix
  logodds <- matrix(-Inf, nA, nB, dimnames = list(ids_A, ids_B))
  for (k in seq_len(nrow(pairs))) {
    i <- match(pairs$id_A[k], ids_A)
    j <- match(pairs$id_B[k], ids_B)
    p <- min(max(pairs$p_match[k], 1e-9), 1 - 1e-9)
    logodds[i, j] <- log(p / (1 - p))
  }

  # Initialise from Hungarian solution
  init <- assign_matches(em_fit, threshold = init_threshold)
  # state: integer vector of length nA, value = column index in B (0 = unmatched)
  state <- integer(nA)
  if (nrow(init) > 0) {
    for (k in seq_len(nrow(init))) {
      state[match(init$id_A[k], ids_A)] <- match(init$id_B[k], ids_B)
    }
  }

  # Log-odds contribution of A-row r assigned to B-column j (0 = unmatched)
  lo <- function(r, j) if (j == 0L) 0 else logodds[r, j]

  samples <- vector("list", K)
  s_idx <- 0L
  total_iter <- n_burn + K * n_thin

  for (it in seq_len(total_iter)) {
    # Propose: pick a random A-index, propose to (un)assign or swap to a
    # random B-index (0 means unmatched). To preserve one-to-one, if proposed
    # B is already assigned to another A, swap them. The proposal is
    # symmetric, so the MH ratio is the posterior ratio.
    i <- sample.int(nA, 1)
    new_j <- sample.int(nB + 1, 1) - 1L   # 0..nB
    old_j <- state[i]

    if (new_j != old_j) {
      conflict <- if (new_j == 0L) integer(0) else which(state == new_j)
      # Only rows i and conflict change, so the log-posterior difference
      # involves just those terms.
      delta <- lo(i, new_j) - lo(i, old_j)
      if (length(conflict) > 0)
        delta <- delta + lo(conflict, old_j) - lo(conflict, new_j)

      if (is.finite(delta) && log(runif(1)) < delta) {
        if (length(conflict) > 0) state[conflict] <- old_j
        state[i] <- new_j
      }
    }

    # Always reached, so exactly K samples are recorded.
    if (it > n_burn && ((it - n_burn) %% n_thin) == 0) {
      s_idx <- s_idx + 1L
      matched_i <- which(state > 0)
      samples[[s_idx]] <- tibble(
        id_A = ids_A[matched_i],
        id_B = ids_B[state[matched_i]]
      )
    }
  }
  samples
}

