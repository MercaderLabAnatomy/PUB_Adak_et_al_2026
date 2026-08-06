# =============================================================================
# diagnostics.R
#
# evaluate_candidate() extracts model-specific diagnostics into one
# standardized object per candidate, regardless of family. Evaluation logic
# only - model fitting lives in model_candidates.R.
#
# Uses {performance} for r2_nakagawa() / icc() so marginal/conditional R^2
# and ICC are computed on the correct (link/latent) scale for the GLMM
# candidates.
#
# Log-normal AIC/BIC receive a Jacobian correction so they are comparable
# to Gaussian/Gamma on the response scale:
#   AIC_corrected = AIC(model on log(y)) + 2 * sum(log(y))
#   BIC_corrected = BIC(model on log(y)) + 2 * sum(log(y))
# (adding 2*sum(log(y)) converts the log-scale log-likelihood back to the
# response-scale log-likelihood before the AIC/BIC penalty is reapplied).
# =============================================================================

.shapiro_safe <- function(resid_vals, alpha = 0.05) {
  n <- length(resid_vals)
  if (n < 3 || n > 5000) return(list(W = NA_real_, p = NA_real_, note = "Shapiro-Wilk not computed (n < 3 or n > 5000)"))
  sh <- tryCatch(shapiro.test(resid_vals), error = function(e) NULL)
  if (is.null(sh)) return(list(W = NA_real_, p = NA_real_, note = "Shapiro-Wilk failed"))
  list(W = round(unname(sh$statistic), 4), p = sh$p.value, note = NA_character_)
}

.dharma_safe <- function(model, n_sim = 1000) {
  out <- list(ks_p = NA_real_, dispersion_p = NA_real_, outlier_p = NA_real_, note = NA_character_, sim = NULL, seed = NA_integer_)
  if (!.has_pkg("DHARMa")) { out$note <- "DHARMa not installed"; return(out) }
  # DHARMa's simulation draws from R's current RNG stream, so results vary
  # run to run unless seeded. No seed is forced here; instead, one is drawn
  # from whatever stream is already running and recorded, so the exact draw
  # can be reproduced later via set.seed(<that value>) if needed.
  seed_used <- sample.int(.Machine$integer.max, 1)
  set.seed(seed_used)
  sim <- tryCatch(DHARMa::simulateResiduals(model, n = n_sim, plot = FALSE), error = function(e) NULL)
  if (is.null(sim)) { out$note <- "DHARMa simulation failed"; return(out) }
  out$sim <- sim
  out$seed <- seed_used
  ks   <- tryCatch(DHARMa::testUniformity(sim, plot = FALSE), error = function(e) NULL)
  disp <- tryCatch(DHARMa::testDispersion(sim, plot = FALSE), error = function(e) NULL)
  outl <- tryCatch(DHARMa::testOutliers(sim, plot = FALSE), error = function(e) NULL)
  if (!is.null(ks))   out$ks_p         <- ks$p.value
  if (!is.null(disp)) out$dispersion_p <- disp$p.value
  if (!is.null(outl)) out$outlier_p    <- outl$p.value
  out
}

.r2_icc_safe <- function(model) {
  out <- list(marginal_r2 = NA_real_, conditional_r2 = NA_real_, icc = NA_real_)
  if (!.has_pkg("performance")) return(out)
  r2 <- tryCatch(performance::r2_nakagawa(model), error = function(e) NULL)
  if (!is.null(r2)) {
    out$marginal_r2    <- unname(r2$R2_marginal)
    out$conditional_r2 <- unname(r2$R2_conditional)
  }
  icc <- tryCatch(performance::icc(model), error = function(e) NULL)
  if (!is.null(icc)) out$icc <- unname(icc$ICC_adjusted)
  out
}

# -----------------------------------------------------------------------
# evaluate_candidate(): dispatcher producing the standardized object
# consumed by scoring.R and the report.
# -----------------------------------------------------------------------
evaluate_candidate <- function(fit, alpha = 0.05) {

  label <- paste0(
    if (fit$family == "gaussian") "Gaussian LMM" else if (fit$family == "gamma") "Gamma GLMM" else "Log-normal LMM",
    " (", describe_re_structure(fit$random_slopes, fit$correlated_slopes), ")"
  )

  excluded <- is_hard_excluded(fit)
  exclusion_reason <- if (excluded) {
    if (!is.null(fit$error)) paste0("Fitting error: ", fit$error)
    else if (!is.na(fit$singular) && fit$singular) "Singular fit"
    else if (!is.na(fit$adequacy_note)) fit$adequacy_note
    else if (is.na(fit$converged) || !fit$converged) "Failed to converge"
    else "Excluded (unspecified structural failure)"
  } else NA_character_

  out <- list(
    family = fit$family, family_label = label,
    random_slopes = fit$random_slopes, formula = fit$formula,
    model = fit$model, excluded = excluded, exclusion_reason = exclusion_reason,
    aic = NA_real_, bic = NA_real_, marginal_r2 = NA_real_, conditional_r2 = NA_real_, icc = NA_real_,
    shapiro_W = NA_real_, shapiro_p = NA_real_,
    ks_p = NA_real_, dispersion_p = NA_real_, outlier_p = NA_real_, dharma_seed = NA_integer_,
    n_convergence_warnings = if (is.na(fit$fit_warnings)) 0L else length(strsplit(fit$fit_warnings, " \\| ")[[1]]),
    convergence_warnings = fit$fit_warnings,
    converged = fit$converged, singular = fit$singular, sim = NULL
  )

  if (excluded || is.null(fit$model)) return(out)

  model <- fit$model
  r2icc <- .r2_icc_safe(model)
  out$marginal_r2    <- r2icc$marginal_r2
  out$conditional_r2 <- r2icc$conditional_r2
  out$icc            <- r2icc$icc

  if (fit$family == "gaussian") {
    diag <- .shapiro_safe(residuals(model), alpha = alpha)
    out$shapiro_W <- diag$W; out$shapiro_p <- diag$p
    out$aic <- tryCatch(stats::AIC(model), error = function(e) NA_real_)
    out$bic <- tryCatch(stats::BIC(model), error = function(e) NA_real_)
  } else {
    dh <- .dharma_safe(model)
    out$ks_p <- dh$ks_p; out$dispersion_p <- dh$dispersion_p; out$outlier_p <- dh$outlier_p
    out$sim <- dh$sim
    out$dharma_seed <- dh$seed
    raw_aic <- tryCatch(stats::AIC(model), error = function(e) NA_real_)
    raw_bic <- tryCatch(stats::BIC(model), error = function(e) NA_real_)
    if (fit$family == "lognormal") {
      # Jacobian correction back to the response scale. Response column
      # name is recovered from the formula's "log(`meas`)" term and looked
      # up in the (already positive-filtered) fitting data.
      resp_col <- sub("^log\\(`(.*)`\\)\\s*~.*$", "\\1", fit$formula)
      y_response <- fit$data_used[[resp_col]]
      jac <- 2 * sum(log(y_response), na.rm = TRUE)
      raw_aic <- raw_aic + jac
      raw_bic <- raw_bic + jac
    }
    out$aic <- raw_aic
    out$bic <- raw_bic
  }

  out
}
