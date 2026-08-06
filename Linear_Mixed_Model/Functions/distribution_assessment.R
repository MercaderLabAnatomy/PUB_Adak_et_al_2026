# =============================================================================
# distribution_assessment.R
#
# Raw-data assessment, computed once per measurement before any mixed model
# is fit. Depends only on the response vector - not on comparison_parameter
# or confounder_vars.
#
# calculate_distribution_metrics() feeds two consumers downstream:
#   - the report's "Raw Data Assessment" table
#   - scoring.R, which turns skewness/kurtosis/positivity/empirical-fit into
#     evidence contributions for the Gamma/Log-normal candidates
#
# The optional fitdistrplus fits are not the inferential model - they ignore
# comparison_parameter/confounder_vars entirely and just ask "which simple
# distribution shape (Normal/Gamma/Log-normal) best matches the raw
# numbers, with no experimental design at all". Their only role is to
# contribute one row of "empirical fit" evidence to the scoring matrix.
# =============================================================================

# Fisher-Pearson skewness (bias-corrected) and excess kurtosis, computed by
# hand so this file has no hard dependency on {moments}/{e1071}.
.skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  g1 <- sum((x - m)^3) / n / s^3
  sqrt(n * (n - 1)) / (n - 2) * g1
}

.kurtosis_excess <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x); s <- sd(x)
  if (s == 0) return(NA_real_)
  g2 <- sum((x - m)^4) / n / s^4
  ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * g2 - 3 * (n - 1)) + 3 - 3
}

calculate_distribution_metrics <- function(y, fit_empirical = TRUE) {

  y_complete <- y[!is.na(y)]
  n <- length(y_complete)

  metrics <- list(
    n                  = n,
    min                = if (n > 0) min(y_complete) else NA_real_,
    max                = if (n > 0) max(y_complete) else NA_real_,
    mean               = if (n > 0) mean(y_complete) else NA_real_,
    median             = if (n > 0) stats::median(y_complete) else NA_real_,
    variance           = if (n > 1) stats::var(y_complete) else NA_real_,
    sd                 = if (n > 1) stats::sd(y_complete) else NA_real_,
    cv                 = NA_real_,
    skewness           = .skewness(y_complete),
    kurtosis_excess    = .kurtosis_excess(y_complete),
    pct_zero           = if (n > 0) mean(y_complete == 0) * 100 else NA_real_,
    pct_negative       = if (n > 0) mean(y_complete < 0) * 100 else NA_real_,
    strictly_positive  = if (n > 0) all(y_complete > 0) else NA,
    right_skewed          = NA,
    strongly_right_skewed = NA,
    heavy_tailed          = NA,
    empirical_fits = NULL,
    empirical_akaike_weights = NULL
  )

  if (n > 0 && !is.na(metrics$mean) && metrics$mean != 0) {
    metrics$cv <- metrics$sd / metrics$mean
  }

  if (!is.na(metrics$skewness)) {
    metrics$right_skewed          <- metrics$skewness > 0.5
    metrics$strongly_right_skewed <- metrics$skewness > 1
  }
  if (!is.na(metrics$kurtosis_excess)) {
    metrics$heavy_tailed <- metrics$kurtosis_excess > 1
  }

  # ---------------------------------------------------------------------
  # Optional empirical distribution fits (fitdistrplus). Structure-free:
  # no random effects, no treatment groups - just the raw response.
  # Only attempted if fit_empirical = TRUE and the package is available;
  # Gamma/Log-normal fits additionally require strictly positive data.
  # ---------------------------------------------------------------------
  if (fit_empirical && .has_pkg("fitdistrplus") && n >= 10) {

    fit_one <- function(distname, data) {
      tryCatch({
        f <- fitdistrplus::fitdist(data, distname)
        list(aic = f$aic, bic = f$bic, loglik = f$loglik, fit = f, error = NULL)
      }, error = function(e) list(aic = NA_real_, bic = NA_real_, loglik = NA_real_, fit = NULL, error = conditionMessage(e)))
    }

    fits <- list(normal = fit_one("norm", y_complete))
    if (isTRUE(metrics$strictly_positive)) {
      fits$gamma    <- fit_one("gamma", y_complete)
      fits$lognormal <- fit_one("lnorm", y_complete)
    }
    metrics$empirical_fits <- fits

    # Akaike weights across whichever fits succeeded, normalized to sum to
    # 1 so they can be used directly as an "empirical fit" evidence
    # contribution per family (see scoring.R).
    aics <- vapply(fits, function(f) f$aic, numeric(1))
    aics <- aics[!is.na(aics)]
    if (length(aics) > 0) {
      delta <- aics - min(aics)
      w <- exp(-0.5 * delta)
      w <- w / sum(w)
      metrics$empirical_akaike_weights <- w
    }
  }

  metrics
}
