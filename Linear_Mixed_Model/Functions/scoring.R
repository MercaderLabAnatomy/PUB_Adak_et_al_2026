# =============================================================================
# scoring.R
#
# The evidence-based scoring engine.
#
# Every sub-score is pre-scaled to [0, 1] with "higher = better" before
# weighting, so the total evidence score is a simple weighted sum with no
# cross-candidate min-max normalization step:
#
#     total_evidence = sum(weight_i * subscore_i)
#
# Missing sub-scores (e.g. Shapiro not computed for n > 5000, or a DHARMa
# test unavailable) have their weight share redistributed among the
# remaining available criteria for that candidate, and a warning is
# recorded.
#
# Hard exclusions (singularity, non-convergence, structural infeasibility)
# are handled entirely upstream in model_candidates.R / is_hard_excluded();
# this file only scores candidates that survived that filter.
# =============================================================================

# -----------------------------------------------------------------------
# Hard-coded scoring weights.
# -----------------------------------------------------------------------
model_selection_weights <- list(
  positivity           = 2,
  skewness             = 2,
  heavy_tails          = 1,
  empirical_fit        = 1,
  marginal_r2          = 2,
  conditional_r2       = 1,
  icc                  = 1,
  residual_diagnostics = 3,
  stability            = 2,
  singularity_penalty  = Inf   # not part of the weighted sum - hard exclusion only
)

# -----------------------------------------------------------------------
# Sub-score builders. Each returns a value in [0,1] or NA (NA = "not
# available for this candidate", handled by weight-redistribution below).
# -----------------------------------------------------------------------

.subscore_positivity <- function(family) {
  if (family == "gaussian") 0 else 1
}

.subscore_skewness <- function(family, skewness) {
  if (family == "gaussian") return(0)
  if (is.na(skewness)) return(NA_real_)
  # Only right-skewness counts as evidence for Gamma/Log-normal; caps at
  # skewness = 2 for full evidence.
  pmin(pmax(skewness, 0) / 2, 1)
}

.subscore_heavy_tails <- function(family, kurtosis_excess) {
  if (family == "gaussian") return(0)
  if (is.na(kurtosis_excess)) return(NA_real_)
  pmin(pmax(kurtosis_excess, 0) / 3, 1)
}

.subscore_empirical_fit <- function(family, akaike_weights) {
  if (is.null(akaike_weights)) return(NA_real_)
  key <- switch(family, gaussian = "normal", gamma = "gamma", lognormal = "lognormal")
  if (!key %in% names(akaike_weights)) return(NA_real_)
  unname(akaike_weights[key])
}

# marginal/conditional R^2 and ICC from performance:: are already
# proportions of variance in [0,1] - used directly, no rescaling needed.
.subscore_identity <- function(x) if (is.na(x)) NA_real_ else pmin(pmax(x, 0), 1)

# Residual diagnostics: p-values are already "higher = more consistent with
# the assumption" and already bounded in [0,1]. For Gamma/Log-normal, this
# averages whichever of KS/dispersion/outlier are available.
.subscore_residual_diagnostics <- function(family, shapiro_p = NA, ks_p = NA, dispersion_p = NA, outlier_p = NA) {
  if (family == "gaussian") {
    return(if (is.na(shapiro_p)) NA_real_ else shapiro_p)
  }
  vals <- c(ks_p, dispersion_p, outlier_p)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

# Stability: penalizes convergence/optimizer warnings continuously rather
# than a pass/fail cutoff. Singular fits never reach here (hard-excluded).
.subscore_stability <- function(n_warnings) {
  if (is.na(n_warnings)) return(NA_real_)
  max(0, 1 - 0.25 * n_warnings)
}

# -----------------------------------------------------------------------
# Scores one evaluated (non-excluded) candidate against
# distribution_metrics (raw-data evidence) and its own evaluate_candidate()
# output. Returns the total evidence score plus a breakdown data.frame row,
# and any redistribution warning.
# -----------------------------------------------------------------------
score_candidate <- function(evald, distribution_metrics, weights = model_selection_weights) {

  raw_scores <- c(
    positivity           = .subscore_positivity(evald$family),
    skewness             = .subscore_skewness(evald$family, distribution_metrics$skewness),
    heavy_tails          = .subscore_heavy_tails(evald$family, distribution_metrics$kurtosis_excess),
    empirical_fit        = .subscore_empirical_fit(evald$family, distribution_metrics$empirical_akaike_weights),
    marginal_r2          = .subscore_identity(evald$marginal_r2),
    conditional_r2       = .subscore_identity(evald$conditional_r2),
    icc                  = .subscore_identity(evald$icc),
    residual_diagnostics = .subscore_residual_diagnostics(evald$family, evald$shapiro_p, evald$ks_p, evald$dispersion_p, evald$outlier_p),
    stability            = .subscore_stability(evald$n_convergence_warnings)
  )

  w <- unlist(weights[names(raw_scores)])
  available <- !is.na(raw_scores)
  warning_msg <- NA_character_

  if (!all(available)) {
    missing_criteria <- names(raw_scores)[!available]
    warning_msg <- paste0(
      evald$family_label, ": criteria unavailable (", paste(missing_criteria, collapse = ", "),
      ") - their weight was redistributed proportionally among the remaining available criteria."
    )
  }

  w_available <- w[available]
  contrib <- setNames(rep(NA_real_, length(raw_scores)), names(raw_scores))
  if (sum(w_available) == 0) {
    total_evidence <- 0
  } else {
    w_rescaled <- w_available / sum(w_available) * sum(w)  # redistribute, keep total weight mass constant
    contrib[available] <- raw_scores[available] * w_rescaled
    total_evidence <- sum(contrib, na.rm = TRUE)
  }

  # Weighted point contribution per criterion (raw subscore x its rescaled
  # weight). This determines Total_evidence, and lets a reader see exactly
  # which criterion a candidate lost ground on relative to another, not
  # just its unweighted 0-1 subscore.
  contrib_pts <- contrib
  names(contrib_pts) <- paste0(names(contrib_pts), "_pts")

  breakdown <- data.frame(
    Family = evald$family_label,
    t(raw_scores),
    t(round(contrib_pts, 4)),
    Total_evidence = round(total_evidence, 4),
    stringsAsFactors = FALSE
  )

  list(total_evidence = total_evidence, breakdown = breakdown, warning = warning_msg)
}

# -----------------------------------------------------------------------
# Builds the comparison table and selects the best candidate across all
# surviving candidates for one measurement (fixed-slope candidates plus
# any independently-fitted random-slope rescues).
# -----------------------------------------------------------------------
compare_candidates <- function(evaluated_candidates, distribution_metrics, weights = model_selection_weights) {

  scored <- lapply(evaluated_candidates, function(evald) {
    if (isTRUE(evald$excluded)) {
      list(total_evidence = NA_real_,
           breakdown = data.frame(Family = evald$family_label, Total_evidence = NA_real_,
                                  Excluded = TRUE, Reason = evald$exclusion_reason, stringsAsFactors = FALSE),
           warning = NA_character_)
    } else {
      s <- score_candidate(evald, distribution_metrics, weights)
      s$breakdown$Excluded <- FALSE
      s$breakdown$Reason <- NA_character_
      s
    }
  })

  comparison_table <- dplyr::bind_rows(lapply(scored, `[[`, "breakdown"))
  warnings_out <- vapply(scored, function(s) if (is.null(s$warning)) NA_character_ else s$warning, character(1))
  warnings_out <- warnings_out[!is.na(warnings_out)]

  eligible_idx <- which(!vapply(evaluated_candidates, function(x) isTRUE(x$excluded), logical(1)))
  if (length(eligible_idx) == 0) {
    return(list(comparison_table = comparison_table, best_index = NA_integer_,
               best_candidate = NULL, warnings = warnings_out))
  }

  totals <- vapply(scored[eligible_idx], `[[`, numeric(1), "total_evidence")
  best_within <- which.max(totals)
  best_index <- eligible_idx[best_within]

  list(
    comparison_table = comparison_table,
    best_index = best_index,
    best_candidate = evaluated_candidates[[best_index]],
    warnings = warnings_out
  )
}
