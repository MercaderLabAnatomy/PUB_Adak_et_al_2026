# =============================================================================
# run_adaptive_model_selection_v2.R
#
# Main entry point for the evidence-based model selection engine. Implements:
#   - input validation
#   - the random-slopes rescue logic: fixed-slope candidates are fitted and
#     evaluated first; any fixed-slope candidate whose primary residual
#     diagnostic fails (p < alpha) triggers an independent random-slopes
#     refit of the same family. If that refit is not hard-excluded
#     (converges, non-singular, adequate group counts) it is added to the
#     candidate set as an independent extra candidate and scored alongside
#     everything else; if it is excluded, it is dropped and logged, never
#     scored.
#   - downstream statistics for the winning candidate
#   - a scoring_log audit trail recording every fit attempt, rescue
#     attempt, and scoring-time warning
#
# Usage: source helper_functions.R, distribution_assessment.R,
# model_candidates.R, diagnostics.R, scoring.R, reporting.R, plotting.R and
# THEN this file, all from one place in your execution script/notebook:
#   result <- run_adaptive_model_selection(
#     datamatrix = my_data, measurement = "meas1",
#     comparison_parameter = "treatment", confounder_vars = c("plateID", "batchNo")
#   )
#   batch <- run_adaptive_model_selection_batch(
#     datamatrix = my_data, measurement_variables = c("meas1", "meas2"),
#     comparison_parameter = "treatment", confounder_vars = c("plateID", "batchNo")
#   )
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(car)
  library(emmeans)
  library(dplyr)
})

# NOTE: this file does not auto-source its sibling modules
# (helper_functions.R, distribution_assessment.R, model_candidates.R,
# diagnostics.R, scoring.R, reporting.R, plotting.R). Source all seven
# files yourself, in that order, from a single place in your execution
# script or notebook, before sourcing this file.

# -----------------------------------------------------------------------
# Support for restricted/defined pairwise comparisons instead of all-pairs.
# `level_pairs` is a list of length-2 character vectors, e.g.:
#   list(c("control_inner", "laser_inner"), c("control_middle", "laser_middle"))
# Each becomes exactly one emmeans contrast (first level minus second
# level), named "A - B" the same way emmeans::pairs() would. Any pair
# referencing a level not present in this measurement's emmGrid (e.g.
# dropped because it had no data at this timepoint) is skipped and
# reported back, rather than erroring out or silently vanishing.
# -----------------------------------------------------------------------
build_defined_contrasts <- function(emm, level_pairs, comparison_parameter) {
  grid_df <- as.data.frame(emm)
  lvls <- as.character(grid_df[[comparison_parameter]])

  contrast_list <- list()
  skipped <- character(0)

  for (pair in level_pairs) {
    a <- as.character(pair[1]); b <- as.character(pair[2])
    if (!(a %in% lvls) || !(b %in% lvls)) {
      skipped <- c(skipped, paste0(a, " - ", b))
      next
    }
    coefs <- rep(0, length(lvls))
    coefs[lvls == a] <- 1
    coefs[lvls == b] <- -1
    contrast_list[[paste0(a, " - ", b)]] <- coefs
  }

  list(contrasts = contrast_list, skipped = skipped)
}

# -----------------------------------------------------------------------
# Downstream statistics for whichever candidate wins scoring: ANOVA
# (Satterthwaite F for Gaussian/Log-normal, Wald Chisq for Gamma), emmeans
# pairwise contrasts with response-scale back-transform for Gamma/
# Log-normal, and ICC per random effect.
#
# `defined_contrasts`: optional list of length-2 character vectors
# specifying exactly which comparisons to run, e.g.
# list(c("control_inner","laser_inner"), c("control_middle","laser_middle"),
# ...). Leave NULL (the default) to run all pairwise comparisons.
# -----------------------------------------------------------------------
extract_downstream_stats <- function(model, family, comparison_parameter, p_adjust_method = "BH",
                                     defined_contrasts = NULL) {
  anova_tbl <- tryCatch({
    if (family == "gamma") {
      tbl <- as.data.frame(car::Anova(model, type = 3, test.statistic = "Chisq"))
      tbl$Effect <- rownames(tbl); rownames(tbl) <- NULL
      tbl <- tbl[!tbl$Effect %in% c("(Intercept)"), , drop = FALSE]
      tbl %>% rename(p_value = `Pr(>Chisq)`) %>%
        mutate(pval_rounded = sapply(p_value, format_pval), signif = sapply(p_value, signif_stars))
    } else {
      tbl <- as.data.frame(anova(model, type = 3))
      tbl$Effect <- rownames(tbl); rownames(tbl) <- NULL
      tbl %>% rename(NumDf = NumDF, DenDf = DenDF, F_value = `F value`, p_value = `Pr(>F)`) %>%
        mutate(pval_rounded = sapply(p_value, format_pval), signif = sapply(p_value, signif_stars))
    }
  }, error = function(e) NULL)

  skipped_contrasts <- character(0)

  pw_tbl <- tryCatch({
    emm <- if (family == "gaussian") {
      emmeans::emmeans(model, specs = as.formula(paste0("~ ", .bt(comparison_parameter))))
    } else {
      emmeans::emmeans(model, specs = as.formula(paste0("~ ", .bt(comparison_parameter))), type = "response")
    }

    if (is.null(defined_contrasts)) {
      tmp <- as.data.frame(pairs(emm, adjust = p_adjust_method))
    } else {
      built <- build_defined_contrasts(emm, defined_contrasts, comparison_parameter)
      skipped_contrasts <<- built$skipped
      if (length(built$contrasts) == 0) {
        tmp <- NULL
      } else {
        tmp <- as.data.frame(emmeans::contrast(emm, method = built$contrasts, adjust = p_adjust_method))
      }
    }

    if (is.null(tmp)) NULL else tmp %>% mutate(pval_rounded = sapply(p.value, format_pval), signif = sapply(p.value, signif_stars))
  }, error = function(e) NULL)

  random_eff_tbl <- tryCatch({
    vc <- as.data.frame(VarCorr(model))
    resid_var <- vc$vcov[vc$grp == "Residual"]
    if (length(resid_var) == 0) resid_var <- NA_real_
    total_var <- sum(vc$vcov, na.rm = TRUE)
    rows <- lapply(unique(vc$grp[vc$grp != "Residual"]), function(rv) {
      re_var <- vc$vcov[vc$grp == rv][1]
      data.frame(Random_effect = rv, Variance = round(re_var, 5),
                 Residual_var = round(resid_var, 5), ICC = round(re_var / total_var, 4),
                 stringsAsFactors = FALSE)
    })
    bind_rows(rows)
  }, error = function(e) NULL)

  list(anova_table = anova_tbl, pairwise_table = pw_tbl, random_eff_table = random_eff_tbl,
       skipped_contrasts = skipped_contrasts)
}

# -----------------------------------------------------------------------
# Primary residual p-value used to decide whether a random-slopes rescue is
# worth attempting for a given fixed-slope candidate.
# -----------------------------------------------------------------------
.primary_residual_p <- function(evald) {
  if (evald$family == "gaussian") evald$shapiro_p else evald$ks_p
}

# -----------------------------------------------------------------------
# Core orchestrator - runs the full evidence-based selection for one
# measurement.
# -----------------------------------------------------------------------
run_adaptive_model_selection <- function(datamatrix,
                                         measurement,
                                         comparison_parameter,
                                         confounder_vars,
                                         correlated_slopes = TRUE,
                                         alpha = 0.05,
                                         fit_empirical = TRUE,
                                         weights = model_selection_weights,
                                         p_adjust_method = "BH",
                                         defined_contrasts = NULL) {

  log_rows <- list()
  add_log <- function(step, description, decision = NA_character_, notes = NA_character_) {
    log_rows[[length(log_rows) + 1]] <<- data.frame(
      Measurement = measurement, Step = step, Description = description,
      Decision = decision, Notes = notes, stringsAsFactors = FALSE
    )
  }

  empty_result <- function(status_note) {
    list(measurement = measurement, scoring_log = bind_rows(log_rows),
         distribution_metrics = NULL, raw_data_table = NULL, sample_size_table = NULL,
         evaluated_candidates = NULL, candidate_table = NULL, comparison_table = NULL,
         performance_table = NULL, residual_table = NULL, score_gap_table = NULL,
         final_model = NULL, final_family = NA_character_, final_family_label = NA_character_,
         final_formula = NA_character_, decision_summary = NA_character_,
         diagnostics_plot = NULL, anova_table = NULL, pairwise_table = NULL, random_eff_table = NULL,
         status = status_note)
  }

  # --- Input validation ---------------------------------------------------
  needed_cols <- c(measurement, comparison_parameter, confounder_vars)
  missing_cols <- setdiff(needed_cols, names(datamatrix))
  if (length(missing_cols) > 0) {
    add_log(1, "Input validation", decision = "ABORTED",
            notes = paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    return(empty_result("Aborted: missing columns"))
  }

  df <- clean_dataset(datamatrix, factor_cols = c(comparison_parameter, confounder_vars))
  df <- df %>% filter(if_all(all_of(needed_cols), ~ !is.na(.)))
  if (nrow(df) == 0) {
    add_log(1, "Input validation", decision = "ABORTED", notes = "No complete-case rows after removing NAs")
    return(empty_result("Aborted: no complete-case rows"))
  }
  add_log(1, "Input validation", decision = "Passed", notes = paste0("n = ", nrow(df), " complete-case rows"))

  sample_size_table <- build_sample_size_table(measurement, df, comparison_parameter)
  add_log(1.5, "Sample sizes", decision = "Computed",
          notes = paste(paste0(sample_size_table$Group, "=", sample_size_table$N), collapse = ", "))

  # --- Raw data assessment (computed once, cached) -------------------------
  distribution_metrics <- calculate_distribution_metrics(df[[measurement]], fit_empirical = fit_empirical)
  raw_data_table <- build_raw_data_table(measurement, distribution_metrics)
  add_log(2, "Raw data assessment", decision = "Computed",
          notes = paste0("Strictly positive = ", distribution_metrics$strictly_positive,
                        "; skewness = ", round(distribution_metrics$skewness, 3),
                        "; kurtosis(excess) = ", round(distribution_metrics$kurtosis_excess, 3)))

  # --- Candidate families based on structural feasibility -------------------
  families <- determine_candidate_families(distribution_metrics)
  add_log(3, "Determine candidate families", decision = paste(families, collapse = ", "))

  # --- Fit + evaluate all fixed-slope candidates ----------------------------
  evaluated_candidates <- list()
  for (fam in families) {
    fit_fun <- switch(fam, gaussian = fit_gaussian_candidate, gamma = fit_gamma_candidate, lognormal = fit_lognormal_candidate)
    fit <- fit_fun(df, measurement, comparison_parameter, confounder_vars,
                  random_slopes = FALSE, correlated_slopes = correlated_slopes)
    evald <- evaluate_candidate(fit, alpha = alpha)
    evaluated_candidates[[length(evaluated_candidates) + 1]] <- evald
    add_log(4, paste0("Fit + evaluate fixed-slope candidate: ", evald$family_label),
            decision = if (evald$excluded) paste0("Excluded (", evald$exclusion_reason, ")") else "Evaluated",
            notes = paste0("AIC = ", round(evald$aic, 2), "; residual diagnostic p = ", format_pval(.primary_residual_p(evald))))
  }

  # --- Random-slopes rescue -------------------------------------------------
  # For each fixed-slope candidate that is not excluded but whose primary
  # residual test fails (p < alpha), attempt a random-slopes refit of the
  # same family. Add as an independent candidate if viable; drop + log if not.
  fixed_candidates <- evaluated_candidates
  for (evald in fixed_candidates) {
    if (evald$excluded) next
    p <- .primary_residual_p(evald)
    if (is.na(p) || p >= alpha) next  # fixed-slope candidate already adequate - no rescue needed

    fit_fun <- switch(evald$family, gaussian = fit_gaussian_candidate, gamma = fit_gamma_candidate, lognormal = fit_lognormal_candidate)
    fit_rs <- fit_fun(df, measurement, comparison_parameter, confounder_vars,
                      random_slopes = TRUE, correlated_slopes = correlated_slopes)
    evald_rs <- evaluate_candidate(fit_rs, alpha = alpha)

    if (evald_rs$excluded) {
      add_log(4.5, paste0("Random-slopes rescue attempt: ", evald$family_label),
              decision = paste0("Dropped (", evald_rs$exclusion_reason, ")"),
              notes = paste0("Triggered by residual diagnostic p = ", format_pval(p), " < alpha = ", alpha))
    } else {
      evaluated_candidates[[length(evaluated_candidates) + 1]] <- evald_rs
      add_log(4.5, paste0("Random-slopes rescue attempt: ", evald$family_label),
              decision = "Added as independent candidate",
              notes = paste0("Triggered by residual diagnostic p = ", format_pval(p), " < alpha = ", alpha))
    }
  }

  # --- Evidence scoring + comparison ----------------------------------------
  comparison_result <- compare_candidates(evaluated_candidates, distribution_metrics, weights = weights)
  if (length(comparison_result$warnings) > 0) {
    for (w in comparison_result$warnings) {
      add_log(6, "Scoring warning: missing criteria", decision = "Weight redistributed", notes = w)
    }
  }
  candidate_table <- build_candidate_table(measurement, evaluated_candidates)
  performance_table <- build_model_performance_table(measurement, evaluated_candidates)
  residual_table <- build_residual_diagnostics_table(measurement, evaluated_candidates)
  score_gap_table <- build_score_gap_table(comparison_result$comparison_table)

  if (is.na(comparison_result$best_index)) {
    add_log(8, "Compare candidates", decision = "ABORTED - no viable candidate")
    result <- empty_result("Aborted: no candidate survived hard-exclusion checks")
    result$distribution_metrics <- distribution_metrics
    result$raw_data_table <- raw_data_table
    result$sample_size_table <- sample_size_table
    result$evaluated_candidates <- evaluated_candidates
    result$candidate_table <- candidate_table
    result$comparison_table <- comparison_result$comparison_table
    result$performance_table <- performance_table
    result$residual_table <- residual_table
    result$score_gap_table <- score_gap_table
    result$decision_summary <- generate_decision_summary(measurement, distribution_metrics, comparison_result, evaluated_candidates)
    return(result)
  }

  best <- comparison_result$best_candidate
  add_log(8, "Compare candidates and select best", decision = paste0("Selected: ", best$family_label))

  decision_summary <- generate_decision_summary(measurement, distribution_metrics, comparison_result, evaluated_candidates)
  add_log(9, "Generate automatic decision summary", decision = "Done")

  # --- Downstream statistics (optional defined_contrasts restricting -------
  # pairwise comparisons) ----------------------------------------------------
  stats_out <- extract_downstream_stats(best$model, best$family, comparison_parameter, p_adjust_method,
                                        defined_contrasts = defined_contrasts)
  add_log(10, "Downstream statistics",
          decision = if (is.null(defined_contrasts)) "All-pairwise comparisons" else paste0(length(defined_contrasts) - length(stats_out$skipped_contrasts), " of ", length(defined_contrasts), " defined comparison(s) computed"),
          notes = if (length(stats_out$skipped_contrasts) > 0) paste0("Skipped (level not present for this measurement): ", paste(stats_out$skipped_contrasts, collapse = "; ")) else NA_character_)

  diag_plot <- if (best$family == "gaussian") {
    plot_gaussian_diagnostics(best$model, measurement)
  } else {
    plot_dharma_diagnostics(best$model, measurement, best$sim, family = best$family)
  }

  list(
    measurement = measurement, scoring_log = bind_rows(log_rows),
    distribution_metrics = distribution_metrics, raw_data_table = raw_data_table,
    sample_size_table = sample_size_table,
    evaluated_candidates = evaluated_candidates, candidate_table = candidate_table,
    comparison_table = comparison_result$comparison_table,
    performance_table = performance_table, residual_table = residual_table,
    score_gap_table = score_gap_table,
    final_model = best$model, final_family = best$family, final_family_label = best$family_label,
    final_formula = best$formula, decision_summary = decision_summary,
    diagnostics_plot = diag_plot,
    anova_table = stats_out$anova_table, pairwise_table = stats_out$pairwise_table,
    random_eff_table = stats_out$random_eff_table,
    status = "OK"
  )
}

# -----------------------------------------------------------------------
# Batch wrapper - runs the evidence-based selection across several
# measurement variables and aggregates everything the report needs.
# -----------------------------------------------------------------------
run_adaptive_model_selection_batch <- function(datamatrix,
                                               measurement_variables,
                                               comparison_parameter,
                                               confounder_vars,
                                               correlated_slopes = TRUE,
                                               alpha = 0.05,
                                               fit_empirical = TRUE,
                                               weights = model_selection_weights,
                                               p_adjust_method = "BH",
                                               defined_contrasts = NULL) {
  results <- list()
  scoring_log_rows <- list()
  summary_rows <- list()
  raw_data_rows <- list()
  sample_size_rows <- list()
  candidate_rows <- list()
  comparison_rows <- list()
  performance_rows <- list()
  residual_rows <- list()
  score_gap_rows <- list()
  decision_summary_rows <- list()
  anova_rows <- list()
  pairwise_rows <- list()
  random_eff_rows <- list()
  diagnostics_plots <- list()

  for (meas in measurement_variables) {
    res <- run_adaptive_model_selection(
      datamatrix = datamatrix, measurement = meas,
      comparison_parameter = comparison_parameter, confounder_vars = confounder_vars,
      correlated_slopes = correlated_slopes, alpha = alpha, fit_empirical = fit_empirical,
      weights = weights, p_adjust_method = p_adjust_method, defined_contrasts = defined_contrasts
    )
    results[[meas]] <- res
    scoring_log_rows[[length(scoring_log_rows) + 1]] <- res$scoring_log

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Measurement = meas,
      Final_family = if (!is.null(res$final_family_label)) res$final_family_label else NA_character_,
      Formula = if (!is.null(res$final_formula)) res$final_formula else NA_character_,
      Status = res$status, stringsAsFactors = FALSE
    )

    decision_summary_rows[[length(decision_summary_rows) + 1]] <- data.frame(
      Measurement = meas,
      Decision_summary = if (!is.null(res$decision_summary)) res$decision_summary else NA_character_,
      stringsAsFactors = FALSE
    )

    if (!is.null(res$raw_data_table))   raw_data_rows[[length(raw_data_rows) + 1]] <- res$raw_data_table
    if (!is.null(res$sample_size_table)) sample_size_rows[[length(sample_size_rows) + 1]] <- res$sample_size_table
    if (!is.null(res$candidate_table))  candidate_rows[[length(candidate_rows) + 1]] <- res$candidate_table
    if (!is.null(res$comparison_table)) comparison_rows[[length(comparison_rows) + 1]] <- res$comparison_table %>% mutate(Measurement = meas, .before = 1)
    if (!is.null(res$performance_table)) performance_rows[[length(performance_rows) + 1]] <- res$performance_table
    if (!is.null(res$residual_table))    residual_rows[[length(residual_rows) + 1]] <- res$residual_table
    if (!is.null(res$score_gap_table))   score_gap_rows[[length(score_gap_rows) + 1]] <- res$score_gap_table %>% mutate(Measurement = meas, .before = 1)
    if (!is.null(res$anova_table))      anova_rows[[length(anova_rows) + 1]] <- res$anova_table %>% mutate(Measurement = meas, .before = 1)
    if (!is.null(res$pairwise_table))   pairwise_rows[[length(pairwise_rows) + 1]] <- res$pairwise_table %>% mutate(Measurement = meas, .before = 1)
    if (!is.null(res$random_eff_table)) random_eff_rows[[length(random_eff_rows) + 1]] <- res$random_eff_table %>% mutate(Measurement = meas, .before = 1)
    if (!is.null(res$diagnostics_plot)) diagnostics_plots[[meas]] <- res$diagnostics_plot
  }

  agg_final_summary <- bind_rows(summary_rows)
  agg_sample_size_table <- bind_rows(sample_size_rows)
  agg_candidate_table <- bind_rows(candidate_rows)
  agg_performance_table <- bind_rows(performance_rows)
  agg_residual_table <- bind_rows(residual_rows)
  agg_anova_table <- bind_rows(anova_rows)
  agg_pairwise_table <- bind_rows(pairwise_rows)

  supplementary_table <- build_supplementary_table(
    final_summary = agg_final_summary, sample_size_table = agg_sample_size_table,
    performance_table = agg_performance_table, residual_table = agg_residual_table,
    candidate_table = agg_candidate_table, anova_table = agg_anova_table,
    pairwise_table = agg_pairwise_table
  )

  list(
    results = results,
    scoring_log = bind_rows(scoring_log_rows),
    final_summary = agg_final_summary,
    decision_summary_table = bind_rows(decision_summary_rows),
    raw_data_table = bind_rows(raw_data_rows),
    sample_size_table = agg_sample_size_table,
    candidate_table = agg_candidate_table,
    comparison_table = bind_rows(comparison_rows),
    performance_table = agg_performance_table,
    residual_table = agg_residual_table,
    score_gap_table = bind_rows(score_gap_rows),
    diagnostics_plots = diagnostics_plots,
    anova_table = agg_anova_table,
    pairwise_table = agg_pairwise_table,
    random_eff_table = bind_rows(random_eff_rows),
    supplementary_table = supplementary_table
  )
}
