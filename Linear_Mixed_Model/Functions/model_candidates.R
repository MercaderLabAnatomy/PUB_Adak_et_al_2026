# =============================================================================
# model_candidates.R
#
# Candidate determination (structural feasibility only, no fitting yet) and
# the independent fitting functions for each candidate family.
#
# Each fit_*_candidate() returns the fitted model plus convergence/
# singularity/warning/error information, and does not evaluate anything -
# that happens in diagnostics.R. This separation means a random-slopes
# refit (see scoring.R) is just another call to the same fit_*_candidate()
# function with random_slopes = TRUE, and can be added as an independent
# extra candidate without touching this file.
# =============================================================================

# -----------------------------------------------------------------------
# Candidate list based purely on structural feasibility: Gaussian is always
# included; Gamma and Log-normal are added only if the response is
# strictly positive.
# -----------------------------------------------------------------------
determine_candidate_families <- function(distribution_metrics) {
  families <- "gaussian"
  if (isTRUE(distribution_metrics$strictly_positive)) {
    families <- c(families, "gamma", "lognormal")
  }
  families
}

# -----------------------------------------------------------------------
# Independent fitting functions, family in {"gaussian","gamma","lognormal"}.
# Gaussian and Log-normal both use lmerTest::lmer() (Log-normal on
# log(y)); Gamma uses lme4::glmer(family = Gamma(link = "log")).
# -----------------------------------------------------------------------
.fit_candidate <- function(family, data, measurement, comparison_parameter, confounder_vars,
                           random_slopes = FALSE, correlated_slopes = TRUE) {
  needs_positive <- family %in% c("gamma", "lognormal")

  df <- data
  n_dropped <- 0L
  if (needs_positive) {
    n_before <- nrow(df)
    df <- df[!is.na(df[[measurement]]) & df[[measurement]] > 0, , drop = FALSE]
    n_dropped <- n_before - nrow(df)
  }

  base_result <- list(
    family = family, random_slopes = random_slopes, correlated_slopes = correlated_slopes,
    model = NULL, formula = NA_character_, data_used = NULL, n_dropped = n_dropped,
    converged = NA, singular = NA, adequacy_note = NA_character_,
    fit_warnings = NA_character_, error = NULL
  )

  if (nrow(df) == 0) {
    base_result$error <- "No rows remain after filtering (NA / non-positive values)"
    return(base_result)
  }

  df[[comparison_parameter]] <- factor(df[[comparison_parameter]])
  for (rv in confounder_vars) df[[rv]] <- factor(df[[rv]])

  formula_str <- build_model_formula(measurement, comparison_parameter, confounder_vars,
                                     random_slopes = random_slopes, correlated_slopes = correlated_slopes,
                                     log_transform = (family == "lognormal"))
  base_result$formula <- formula_str
  base_result$data_used <- df

  fit_warnings <- character(0)
  model <- withCallingHandlers(
    tryCatch({
      if (family == "gamma") {
        lme4::glmer(as.formula(formula_str), data = df, family = Gamma(link = "log"),
                    control = lme4::glmerControl(optimizer = "bobyqa"))
      } else {
        lmerTest::lmer(as.formula(formula_str), data = df, REML = TRUE)
      }
    }, error = function(e) e),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  base_result$fit_warnings <- if (length(fit_warnings) == 0) NA_character_ else paste(unique(fit_warnings), collapse = " | ")

  if (inherits(model, "error")) {
    base_result$error <- conditionMessage(model)
    base_result$converged <- FALSE
    return(base_result)
  }

  base_result$model <- model
  base_result$converged <- TRUE
  base_result$singular <- tryCatch(lme4::isSingular(model), error = function(e) NA)
  if (random_slopes) {
    base_result$adequacy_note <- check_re_structure_adequacy(model, confounder_vars)
  }

  base_result
}

fit_gaussian_candidate <- function(data, measurement, comparison_parameter, confounder_vars,
                                   random_slopes = FALSE, correlated_slopes = TRUE) {
  .fit_candidate("gaussian", data, measurement, comparison_parameter, confounder_vars,
                random_slopes, correlated_slopes)
}

fit_gamma_candidate <- function(data, measurement, comparison_parameter, confounder_vars,
                                random_slopes = FALSE, correlated_slopes = TRUE) {
  .fit_candidate("gamma", data, measurement, comparison_parameter, confounder_vars,
                random_slopes, correlated_slopes)
}

fit_lognormal_candidate <- function(data, measurement, comparison_parameter, confounder_vars,
                                    random_slopes = FALSE, correlated_slopes = TRUE) {
  .fit_candidate("lognormal", data, measurement, comparison_parameter, confounder_vars,
                random_slopes, correlated_slopes)
}

# -----------------------------------------------------------------------
# Hard exclusion rules. Only structural failures disqualify a candidate -
# residual diagnostics are never used here; they are evidence, handled
# entirely in scoring.R.
# -----------------------------------------------------------------------
is_hard_excluded <- function(fit) {
  if (!is.null(fit$error)) return(TRUE)
  if (is.na(fit$converged) || !fit$converged) return(TRUE)
  if (!is.na(fit$singular) && fit$singular) return(TRUE)
  if (!is.na(fit$adequacy_note)) return(TRUE)  # underpowered random-slope structure
  FALSE
}
