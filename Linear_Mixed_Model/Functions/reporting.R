# =============================================================================
# reporting.R
#
# Generates the automatic decision-summary text and the table builders that
# feed the report template. The .Rmd itself only renders these; all data
# preparation happens here.
# =============================================================================

.format_pct <- function(x, digits = 1) if (is.na(x)) "NA" else paste0(round(x, digits), "%")

# -----------------------------------------------------------------------
# Automatic decision summary, generated from distribution_metrics, the
# comparison table, and the winning candidate. Wording is templated but
# driven entirely by the actual numbers for that measurement.
# -----------------------------------------------------------------------
generate_decision_summary <- function(measurement, distribution_metrics, comparison_result, evaluated_candidates) {

  if (is.na(comparison_result$best_index)) {
    return(paste0(
      "For '", measurement, "', no candidate model survived the hard-exclusion checks ",
      "(convergence, singularity, or structural feasibility). Manual review is required - ",
      "see the comparison table for exclusion reasons per candidate."
    ))
  }

  best <- evaluated_candidates[[comparison_result$best_index]]

  positivity_txt <- if (isTRUE(distribution_metrics$strictly_positive)) {
    "strictly positive, making Gamma and Log-normal plausible candidate distributions alongside the Gaussian baseline"
  } else {
    "not strictly positive, so only the Gaussian LMM was structurally feasible"
  }

  skew_txt <- if (!is.na(distribution_metrics$skewness)) {
    if (isTRUE(distribution_metrics$strongly_right_skewed)) "strongly right-skewed"
    else if (isTRUE(distribution_metrics$right_skewed)) "moderately right-skewed"
    else "not markedly skewed"
  } else "of undetermined skewness (insufficient data)"

  n_candidates <- sum(!vapply(evaluated_candidates, function(x) isTRUE(x$excluded), logical(1)))
  n_excluded   <- sum(vapply(evaluated_candidates, function(x) isTRUE(x$excluded), logical(1)))

  excluded_txt <- if (n_excluded > 0) {
    reasons <- vapply(evaluated_candidates, function(x) {
      if (isTRUE(x$excluded)) paste0(x$family_label, " was excluded (", x$exclusion_reason, ").") else NA_character_
    }, character(1))
    paste(stats::na.omit(reasons), collapse = " ")
  } else ""

  paste0(
    "The response variable for '", measurement, "' was ", positivity_txt, " and ", skew_txt, ". ",
    length(evaluated_candidates), " candidate model(s) were fitted, of which ", n_candidates,
    " survived hard-exclusion checks. ", excluded_txt,
    " The ", best$family_label, " achieved the strongest overall evidence score (",
    round(comparison_result$comparison_table$Total_evidence[comparison_result$comparison_table$Family == best$family_label][1], 3),
    ") and was therefore selected for downstream statistical inference."
  )
}

# -----------------------------------------------------------------------
# Raw data assessment table.
# -----------------------------------------------------------------------
build_raw_data_table <- function(measurement, distribution_metrics) {
  data.frame(
    Measurement = measurement,
    N = distribution_metrics$n,
    Mean = round(distribution_metrics$mean, 4),
    SD = round(distribution_metrics$sd, 4),
    CV = round(distribution_metrics$cv, 4),
    Skewness = round(distribution_metrics$skewness, 4),
    Kurtosis_excess = round(distribution_metrics$kurtosis_excess, 4),
    Strictly_positive = distribution_metrics$strictly_positive,
    Right_skewed = distribution_metrics$right_skewed,
    Heavy_tailed = distribution_metrics$heavy_tailed,
    Best_empirical_fit = if (!is.null(distribution_metrics$empirical_akaike_weights)) {
      names(distribution_metrics$empirical_akaike_weights)[which.max(distribution_metrics$empirical_akaike_weights)]
    } else NA_character_,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------
# Candidate models table: which candidates were fitted, converged,
# singular, and - if excluded - why.
# -----------------------------------------------------------------------
build_candidate_table <- function(measurement, evaluated_candidates) {
  dplyr::bind_rows(lapply(evaluated_candidates, function(x) {
    data.frame(
      Measurement = measurement, Candidate = x$family_label,
      Converged = isTRUE(x$converged), Singular = x$singular,
      Excluded = isTRUE(x$excluded), Reason = x$exclusion_reason,
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------------------------------------------------
# Model performance table (marginal R2, conditional R2, ICC, AIC, BIC),
# pulled directly from each evaluated candidate rather than from the
# evidence comparison table, since AIC/BIC are not part of the weighted
# evidence score.
# -----------------------------------------------------------------------
build_model_performance_table <- function(measurement, evaluated_candidates) {
  dplyr::bind_rows(lapply(evaluated_candidates, function(x) {
    data.frame(
      Measurement = measurement, Candidate = x$family_label, Excluded = isTRUE(x$excluded),
      AIC = round(x$aic, 2), BIC = round(x$bic, 2),
      Marginal_R2 = round(x$marginal_r2, 4), Conditional_R2 = round(x$conditional_r2, 4),
      ICC = round(x$icc, 4),
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------------------------------------------------
# Residual diagnostics table (Shapiro-Wilk for Gaussian; DHARMa KS/
# dispersion/outlier for Gamma & Log-normal). DHARMa_Seed records the
# random seed drawn for that candidate's simulation, so its exact
# DHARMa-based diagnostics can be reproduced later via set.seed().
# -----------------------------------------------------------------------
build_residual_diagnostics_table <- function(measurement, evaluated_candidates) {
  dplyr::bind_rows(lapply(evaluated_candidates, function(x) {
    data.frame(
      Measurement = measurement, Candidate = x$family_label, Excluded = isTRUE(x$excluded),
      Shapiro_W = round(x$shapiro_W, 4), Shapiro_p = x$shapiro_p,
      DHARMa_KS_p = x$ks_p, DHARMa_Dispersion_p = x$dispersion_p, DHARMa_Outlier_p = x$outlier_p,
      DHARMa_Seed = x$dharma_seed,
      stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------------------------------------------------
# Score gap / ranking table: takes the per-measurement comparison_table
# (which carries both the raw 0-1 subscores and their weighted point
# contributions in "<criterion>_pts" columns) and turns it into a ranked
# leaderboard - which candidate came 1st/2nd/3rd, how far behind the
# winner each one is, and how much closer the runner-up is than the
# next-lowest candidate. Comparing the "<criterion>_pts" columns row by
# row shows exactly where one candidate lost ground relative to another.
# -----------------------------------------------------------------------
build_score_gap_table <- function(comparison_table) {
  if (is.null(comparison_table) || nrow(comparison_table) == 0) return(comparison_table)

  pts_cols <- grep("_pts$", names(comparison_table), value = TRUE)
  keep_cols <- c("Family", pts_cols, "Total_evidence", "Excluded", "Reason")
  keep_cols <- intersect(keep_cols, names(comparison_table))
  tbl <- comparison_table[, keep_cols, drop = FALSE]

  eligible <- tbl[!isTRUE_vec(tbl$Excluded) & !is.na(tbl$Total_evidence), , drop = FALSE]
  excluded <- tbl[isTRUE_vec(tbl$Excluded) | is.na(tbl$Total_evidence), , drop = FALSE]

  if (nrow(eligible) == 0) {
    excluded$Rank <- NA_integer_; excluded$Trail_behind_best <- NA_real_; excluded$Lead_over_next <- NA_real_
    return(excluded)
  }

  ord <- order(-eligible$Total_evidence)
  eligible <- eligible[ord, , drop = FALSE]
  eligible$Rank <- seq_len(nrow(eligible))
  best_score <- eligible$Total_evidence[1]
  eligible$Trail_behind_best <- round(best_score - eligible$Total_evidence, 4)
  eligible$Lead_over_next <- c(round(-diff(eligible$Total_evidence), 4), NA_real_)

  if (nrow(excluded) > 0) {
    excluded$Rank <- NA_integer_; excluded$Trail_behind_best <- NA_real_; excluded$Lead_over_next <- NA_real_
  }

  result <- dplyr::bind_rows(eligible, excluded)
  # Ranking columns go up front, right after Family, for readability.
  front <- c("Family", "Rank", "Total_evidence", "Trail_behind_best", "Lead_over_next")
  result[, c(intersect(front, names(result)), setdiff(names(result), front))]
}

isTRUE_vec <- function(x) vapply(x, isTRUE, logical(1))

# -----------------------------------------------------------------------
# Per-group and total sample sizes for one measurement, computed from the
# complete-case dataset used for fitting (post NA-filtering, before any
# further per-family positivity filtering). Long format (one row per
# group, plus a "Total" row) so it stays flexible across measurements with
# different numbers of comparison_parameter levels.
# -----------------------------------------------------------------------
build_sample_size_table <- function(measurement, df, comparison_parameter) {
  grp_counts <- table(df[[comparison_parameter]])
  dplyr::bind_rows(
    data.frame(Measurement = measurement, Group = names(grp_counts),
               N = as.integer(grp_counts), stringsAsFactors = FALSE),
    data.frame(Measurement = measurement, Group = "Total",
               N = nrow(df), stringsAsFactors = FALSE)
  )
}

# -----------------------------------------------------------------------
# R version and package versions actually used, for reproducibility.
# -----------------------------------------------------------------------
build_session_info_table <- function() {
  pkgs <- c("lme4", "lmerTest", "car", "emmeans", "dplyr", "performance",
           "DHARMa", "fitdistrplus", "ggplot2", "patchwork", "kableExtra",
           "knitr", "rmarkdown", "openxlsx")
  pkg_versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))

  data.frame(
    Field = c("R_version", "Platform", "Run_date", pkgs),
    Value = c(paste(R.version$major, R.version$minor, sep = "."),
             R.version$platform,
             as.character(Sys.Date()),
             unname(pkg_versions)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------
# Adds any columns in `cols` that don't already exist in df, filled with
# NA. Used to safely reference columns (like "ratio" vs "estimate", or
# "F_value" vs "Chisq") that only exist for some measurements after
# bind_rows-ing heterogeneous per-family tables across a batch run.
# -----------------------------------------------------------------------
.ensure_cols <- function(df, cols) {
  for (co in setdiff(cols, names(df))) df[[co]] <- NA
  df
}

# -----------------------------------------------------------------------
# One consolidated, publication-ready table: one row per pairwise
# contrast, with the measurement's selected model, total N, omnibus
# (Type III) test, R^2/ICC, the contrast's effect size (ratio for Gamma/
# Log-normal, difference for Gaussian - unified into one column), and a
# Caveat column flagging measurements where model selection was a
# fallback (only one structurally feasible candidate) or where the
# selected model's residual diagnostics indicated a poor fit.
#
# Built once per batch run (needs the already-aggregated, cross-
# measurement tables) - called from run_adaptive_model_selection_batch()
# after all the other tables exist.
# -----------------------------------------------------------------------
build_supplementary_table <- function(final_summary, sample_size_table, performance_table,
                                      residual_table, candidate_table, anova_table, pairwise_table) {
  if (is.null(pairwise_table) || nrow(pairwise_table) == 0) return(NULL)

  model_info <- final_summary %>% select(Measurement, Model_selected = Final_family)

  n_info <- if (!is.null(sample_size_table) && nrow(sample_size_table) > 0) {
    sample_size_table %>% filter(Group == "Total") %>% select(Measurement, N_total = N)
  } else data.frame(Measurement = character(0), N_total = integer(0))

  perf_info <- if (!is.null(performance_table) && nrow(performance_table) > 0) {
    performance_table %>%
      left_join(model_info, by = "Measurement") %>%
      filter(Candidate == Model_selected) %>%
      select(Measurement, Marginal_R2, Conditional_R2, ICC)
  } else data.frame(Measurement = character(0))

  anova_norm <- if (!is.null(anova_table) && nrow(anova_table) > 0) {
    .ensure_cols(anova_table, c("F_value", "Chisq", "NumDf", "DenDf", "Df", "p_value", "pval_rounded")) %>%
      mutate(
        Omnibus_stat_type = ifelse(!is.na(F_value), "F", "Chisq"),
        Omnibus_stat = ifelse(!is.na(F_value), round(F_value, 3), round(Chisq, 3)),
        Omnibus_df = ifelse(!is.na(NumDf), paste0(round(NumDf, 1), ", ", round(DenDf, 1)), as.character(Df)),
        # Uses the already-computed, format_pval()-safe string rather than
        # the raw numeric p_value: a very small p-value rounds to exactly 0
        # as a number (displaying as literal "0" in Excel/Word), whereas
        # format_pval() renders it as "<0.001" instead.
        Omnibus_p = pval_rounded
      ) %>%
      select(Measurement, Omnibus_stat_type, Omnibus_stat, Omnibus_df, Omnibus_p)
  } else data.frame(Measurement = character(0))

  # Caveats: (a) selection was a fallback because only one candidate was
  # structurally feasible/survived, (b) the selected model's own residual
  # diagnostics indicate a poor fit despite winning the comparison.
  fallback_meas <- character(0)
  if (!is.null(candidate_table) && nrow(candidate_table) > 0) {
    fallback_meas <- candidate_table %>% filter(!Excluded) %>%
      count(Measurement, name = "n_candidates") %>%
      filter(n_candidates <= 1) %>% pull(Measurement)
  }
  poor_fit_meas <- character(0)
  if (!is.null(residual_table) && nrow(residual_table) > 0) {
    poor_fit_meas <- residual_table %>% left_join(model_info, by = "Measurement") %>%
      filter(Candidate == Model_selected) %>%
      .ensure_cols(c("Shapiro_p", "DHARMa_KS_p")) %>%
      filter((!is.na(Shapiro_p) & Shapiro_p < 0.05) | (!is.na(DHARMa_KS_p) & DHARMa_KS_p < 0.05)) %>%
      pull(Measurement)
  }
  caveat_text <- function(meas) {
    notes <- c()
    if (meas %in% fallback_meas) notes <- c(notes, "Only one structurally feasible candidate (fallback selection)")
    if (meas %in% poor_fit_meas) notes <- c(notes, "Selected model's residual diagnostics indicate imperfect fit")
    if (length(notes) == 0) NA_character_ else paste(notes, collapse = "; ")
  }

  pw <- .ensure_cols(pairwise_table, c("ratio", "estimate", "SE", "df", "t.ratio", "z.ratio",
                                       "p.value", "pval_rounded", "signif")) %>%
    mutate(
      Effect_size = ifelse(!is.na(ratio), ratio, estimate),
      Effect_size_type = ifelse(!is.na(ratio), "ratio (response scale)", "difference (response scale)"),
      # GLMM contrasts (Gamma) are asymptotic Wald z-tests, not t-tests -
      # emmeans reports these with df = Inf and a "z.ratio" column instead
      # of "t.ratio" (a t-distribution with df = Inf is the normal
      # distribution, so this is a real, correct statistical distinction,
      # not missing data). Writing a literal Inf to Excel produces #NUM!,
      # so the test statistic and its type are unified into one column
      # pair (like Effect_size/Effect_size_type above), and df is left
      # blank specifically for z-tests rather than written as Inf.
      Contrast_stat = ifelse(!is.na(t.ratio), t.ratio, z.ratio),
      Contrast_stat_type = ifelse(!is.na(t.ratio), "t", "z"),
      Contrast_df = ifelse(is.infinite(df), NA_real_, df)
    )

  pw %>%
    left_join(model_info, by = "Measurement") %>%
    left_join(n_info, by = "Measurement") %>%
    left_join(perf_info, by = "Measurement") %>%
    left_join(anova_norm, by = "Measurement") %>%
    rowwise() %>%
    mutate(Caveat = caveat_text(Measurement)) %>%
    ungroup() %>%
    transmute(
      Measurement, Model_selected, N_total,
      Omnibus_stat_type, Omnibus_stat, Omnibus_df,
      Omnibus_p,  # already a format_pval()-safe string (e.g. "<0.001"), not rounded further
      Marginal_R2 = round(Marginal_R2, 4), Conditional_R2 = round(Conditional_R2, 4), ICC = round(ICC, 4),
      Contrast = contrast, Effect_size = round(Effect_size, 4), Effect_size_type,
      SE = round(SE, 4),
      Contrast_stat_type, Contrast_stat = round(Contrast_stat, 3),
      Contrast_df = round(Contrast_df, 2),
      # Same fix as Omnibus_p: uses the format_pval()-safe string instead
      # of a bare round(p.value, 4), which would otherwise display as
      # literal "0" for very small p-values. This also makes the separate
      # pval_rounded column (identical to this) redundant, so it's
      # dropped here.
      Contrast_p_adj = pval_rounded, Signif = signif,
      Caveat
    )
}

# -----------------------------------------------------------------------
# Compresses the wide supplementary_table (20 columns, meant for full
# machine-readable transparency) into a lean, manuscript-ready version:
# - Omnibus test and contrast test statistic each collapse from 3 columns
#   into one conventionally-formatted string ("F(7, 238.6) = 2.63",
#   "z = -2.25"), matching standard reporting style.
# - Effect_size and SE collapse into one "estimate (SE)" column.
# - pval_rounded is dropped (redundant with Contrast_p_adj).
# - Caveat text is replaced with a superscript-style footnote marker
#   (a, b, c, ...); a second small lookup table maps each marker back to
#   its full text, meant to be placed as footnotes below the table rather
#   than repeating the same sentence on every row for a given measurement.
# Returns list(table = <lean data.frame>, footnotes = <marker/text lookup>).
# -----------------------------------------------------------------------
build_manuscript_table <- function(supplementary_table) {
  if (is.null(supplementary_table) || nrow(supplementary_table) == 0) {
    return(list(table = supplementary_table, footnotes = NULL))
  }

  st <- supplementary_table

  omnibus_str <- ifelse(
    is.na(st$Omnibus_stat_type), NA_character_,
    ifelse(st$Omnibus_stat_type == "F",
          paste0("F(", st$Omnibus_df, ") = ", st$Omnibus_stat),
          paste0("\u03c7\u00b2(", st$Omnibus_df, ") = ", st$Omnibus_stat))
  )

  contrast_str <- ifelse(
    is.na(st$Contrast_stat_type), NA_character_,
    ifelse(st$Contrast_stat_type == "t" & !is.na(st$Contrast_df),
          paste0("t(", st$Contrast_df, ") = ", st$Contrast_stat),
          paste0(st$Contrast_stat_type, " = ", st$Contrast_stat))
  )

  effect_str <- ifelse(
    is.na(st$Effect_size), NA_character_,
    paste0(st$Effect_size, " (", st$SE, ")")
  )

  # Footnote markers: one letter per unique non-missing Caveat, assigned in
  # first-encountered order so they read a, b, c... down the table.
  unique_caveats <- unique(st$Caveat[!is.na(st$Caveat)])
  markers <- setNames(letters[seq_along(unique_caveats)], unique_caveats)
  note_col <- ifelse(is.na(st$Caveat), "", unname(markers[st$Caveat]))

  lean <- data.frame(
    target = if ("target" %in% names(st)) st$target else NA,
    Measurement = st$Measurement,
    Contrast = st$Contrast,
    Model = st$Model_selected,
    N = st$N_total,
    Omnibus = omnibus_str,
    Omnibus_p = st$Omnibus_p,
    R2m = st$Marginal_R2, R2c = st$Conditional_R2, ICC = st$ICC,
    `Effect (SE)` = effect_str,
    Test_stat = contrast_str,
    p_adj = st$Contrast_p_adj,
    Sig = st$Signif,
    Note = note_col,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!"target" %in% names(st)) lean$target <- NULL

  footnotes <- if (length(unique_caveats) == 0) {
    NULL
  } else {
    data.frame(Marker = unname(markers), Text = names(markers), stringsAsFactors = FALSE)
  }

  list(table = lean, footnotes = footnotes)
}
