# =============================================================================
# power_analysis_functions.R
#
# Shared functions for post-hoc power analysis: required sample size and
# achieved power for two-sample comparisons, based on Cohen's d effect size
# (unequal-variance correction) and a two-sample t-test power calculation.
# Used identically by the 3d and 7d execution notebooks.
# =============================================================================

# -----------------------------------------------------------------------
# Runs a two-sample power analysis for group1 vs group2 on measurement_col,
# optionally repeated across every combination of grouping_vars (e.g. one
# row per timepoint). Returns one data.frame row per combination with N,
# mean, SD, CI, Cohen's d, required n, achieved power, and an Adequate flag.
# -----------------------------------------------------------------------
perform_power_analysis <- function(
    data,
    measurement_col,
    treatment_col,
    group1,
    group2,
    grouping_vars = NULL,
    power = 0.9,
    alpha = 0.01,
    save_excel = FALSE,
    excel_path = NULL,
    analysis_name = "Power Analysis"
) {
  # Setup grouping levels
  if (is.null(grouping_vars) || length(grouping_vars) == 0) {
    grouping_levels <- list(dummy = "all")
  } else {
    grouping_levels <- lapply(grouping_vars[grouping_vars %in% colnames(data)],
                              function(v) unique(data[[v]]))
    names(grouping_levels) <- grouping_vars[grouping_vars %in% colnames(data)]
    if (length(grouping_levels) == 0) grouping_levels <- list(dummy = "all")
  }

  # Generate all combinations
  combinations <- expand.grid(rev(grouping_levels), stringsAsFactors = FALSE)
  combinations <- combinations[, rev(names(combinations)), drop = FALSE]

  results <- data.frame()

  for (i in 1:nrow(combinations)) {
    df_subset <- data

    # Apply filters for this combination
    for (var in names(combinations)) {
      if (var != "dummy") df_subset <- df_subset[df_subset[[var]] == combinations[i, var], ]
    }

    # Extract and clean groups
    g1 <- na.omit(df_subset[[measurement_col]][df_subset[[treatment_col]] == group1])
    g2 <- na.omit(df_subset[[measurement_col]][df_subset[[treatment_col]] == group2])

    if (length(g1) < 2 || length(g2) < 2) next

    n1 <- length(g1); n2 <- length(g2)
    mean1 <- mean(g1); mean2 <- mean(g2)
    sd1 <- sd(g1); sd2 <- sd(g2)

    if (sd1 == 0 || sd2 == 0) next

    # Confidence intervals
    sem1 <- sd1 / sqrt(n1); sem2 <- sd2 / sqrt(n2)
    t_crit1 <- qt(1 - alpha / 2, n1 - 1); t_crit2 <- qt(1 - alpha / 2, n2 - 1)

    # Effect size and power
    es <- cohensD(g1, g2, method = "unequal")
    if (is.na(es)) next

    # Effect sizes near zero are kept in the output but use |d| = 0.01 for
    # the power/sample-size calculations, which are undefined at d = 0.
    small_es_flag <- FALSE
    if (abs(es) < 0.01) {
      small_es_flag <- TRUE
      es_used <- sign(es) * 0.01
    } else {
      es_used <- es
    }

    pwr_req <- tryCatch({
      pwr.t.test(d = es_used, sig.level = alpha, power = power, type = "two.sample")
    }, error = function(e) list(n = NA))

    pwr_ach <- tryCatch({
      pwr.t.test(n = min(n1, n2), d = es_used, sig.level = alpha, type = "two.sample")
    }, error = function(e) list(power = NA))

    # Build result row
    result_row <- combinations[i, , drop = FALSE]
    result_row$Group1_Name <- group1
    result_row$Group2_Name <- group2

    result_row[[paste0(group1, "_N")]] <- n1
    result_row[[paste0(group2, "_N")]] <- n2
    result_row[[paste0(group1, "_Mean")]] <- round(mean1, 3)
    result_row[[paste0(group2, "_Mean")]] <- round(mean2, 3)
    result_row[[paste0(group1, "_SD")]] <- round(sd1, 3)
    result_row[[paste0(group2, "_SD")]] <- round(sd2, 3)
    result_row[[paste0(group1, "_CI_Low")]] <- round(mean1 - t_crit1 * sem1, 3)
    result_row[[paste0(group1, "_CI_High")]] <- round(mean1 + t_crit1 * sem1, 3)
    result_row[[paste0(group2, "_CI_Low")]] <- round(mean2 - t_crit2 * sem2, 3)
    result_row[[paste0(group2, "_CI_High")]] <- round(mean2 + t_crit2 * sem2, 3)

    # Cohen's d stored as text for Excel compatibility
    cohens_text <- ifelse(small_es_flag || is.na(es), "NA", format(round(es, 3), nsmall = 3))
    result_row$Cohens_d <- as.character(cohens_text)
    result_row$Required_n <- ifelse(is.na(pwr_req$n), NA, ceiling(pwr_req$n))
    result_row$Achieved_power <- ifelse(is.na(pwr_ach$power), NA, round(pwr_ach$power, 3))
    result_row$Adequate <- ifelse(small_es_flag, "Effect too small",
                                   ifelse(is.na(pwr_req$n), "Calculation error",
                                          ifelse(min(n1, n2) >= ceiling(pwr_req$n), "Yes", "No")))

    results <- rbind(results, result_row)
  }

  if ("dummy" %in% colnames(results)) results$dummy <- NULL

  # Optional direct-to-Excel export. Unused by the execution notebooks,
  # which export separately via save_power_results() after combining
  # several perform_power_analysis() calls; kept here for standalone use.
  if (save_excel && !is.null(excel_path) && nrow(results) > 0) {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      install.packages("writexl", repos = "https://cloud.r-project.org/")
    }

    get_pkg_version <- function(pkg) {
      tryCatch(as.character(packageVersion(pkg)), error = function(e) "not installed")
    }

    params_df <- data.frame(
      Parameter = c(
        "Analysis Name", "Measurement Column", "Group Column", "Group 1", "Group 2",
        "Power", "Alpha", "Grouping Variables", "",
        "R Version", "Package: pwr", "Package: readxl", "Package: writexl",
        "Package: lsr", "Package: dplyr", "Package: stringr", "",
        "Methods Description"
      ),
      Value = c(
        analysis_name, measurement_col, treatment_col, group1, group2, power, alpha,
        paste(grouping_vars, collapse = ", "), "",
        R.version.string, get_pkg_version("pwr"), get_pkg_version("readxl"),
        get_pkg_version("writexl"), get_pkg_version("lsr"), get_pkg_version("dplyr"),
        get_pkg_version("stringr"), "",
        "Power analysis for two-sample comparisons using Cohen's d effect size (unequal variance correction) and two-sample t-test power calculations. Required sample sizes calculated for specified power level. Confidence intervals based on t-distribution."
      )
    )

    writexl::write_xlsx(
      list(Results = results, Parameters = params_df),
      gsub("\\.csv$", ".xlsx", excel_path)
    )
  }

  return(results)
}

# -----------------------------------------------------------------------
# Converts one perform_power_analysis() result into the common long-format
# schema shared across comparisons: one Treatment/Comparison/Group1.../
# Group2... row per combination, so results from different calls (e.g. one
# per timepoint pair) can be row-bound together. treatment_value can be a
# single value (recycled) or a vector matching nrow(result).
# -----------------------------------------------------------------------
standardize_result <- function(result, treatment_value, comparison_label, group1, group2) {
  if (nrow(result) == 0) return(result)
  data.frame(
    Treatment = treatment_value,
    Comparison = comparison_label,
    Group1_Name = group1,
    Group2_Name = group2,
    N_Group1 = result[[paste0(group1, "_N")]],
    N_Group2 = result[[paste0(group2, "_N")]],
    Mean_Group1 = result[[paste0(group1, "_Mean")]],
    Mean_Group2 = result[[paste0(group2, "_Mean")]],
    SD_Group1 = result[[paste0(group1, "_SD")]],
    SD_Group2 = result[[paste0(group2, "_SD")]],
    CI_Low_Group1 = result[[paste0(group1, "_CI_Low")]],
    CI_High_Group1 = result[[paste0(group1, "_CI_High")]],
    CI_Low_Group2 = result[[paste0(group2, "_CI_Low")]],
    CI_High_Group2 = result[[paste0(group2, "_CI_High")]],
    Cohens_d = result$Cohens_d,
    Required_n = result$Required_n,
    Achieved_power = result$Achieved_power,
    Adequate = result$Adequate,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------
# Renames Label_GroupN columns (e.g. "Mean_Group1") to Group_Label form
# (e.g. "Group1_Mean") for the final combined results table.
# -----------------------------------------------------------------------
normalize_group_cols <- function(df) {
  cn <- colnames(df)
  labels <- c("N", "Mean", "SD", "CI_Low", "CI_High")
  for (lab in labels) {
    pref_regex <- paste0("^", lab, "_(.+)$")
    matches <- grepl(pref_regex, cn)
    if (any(matches)) {
      groups <- sub(pref_regex, "\\1", cn[matches])
      cn[matches] <- paste0(groups, "_", lab)
    }
  }
  colnames(df) <- cn
  df
}

# -----------------------------------------------------------------------
# Runs the within-treatment timepoint-vs-timepoint comparisons: for every
# combination of group_levels (e.g. treatments) x timepoint_pairs, subsets
# the data to that group and calls perform_power_analysis() comparing the
# two timepoints. Prints a per-pair summary and returns one combined,
# standardized results table.
# -----------------------------------------------------------------------
run_timepoint_comparisons <- function(data, measurement_col, group_levels,
                                      timepoint_pairs, power, alpha,
                                      group_col = "treatment", compare_col = "timepoint",
                                      label_prefix = "") {
  results_list <- list()

  for (grp in group_levels) {
    cat("\n--- Treatment:", grp, "---\n\n")
    df_group <- data[data[[group_col]] == grp, ]

    for (pair in timepoint_pairs) {
      tp1 <- pair[1]; tp2 <- pair[2]

      result <- perform_power_analysis(
        data = df_group, measurement_col = measurement_col, treatment_col = compare_col,
        group1 = tp1, group2 = tp2, grouping_vars = NULL,
        power = power, alpha = alpha, save_excel = FALSE,
        analysis_name = paste0(label_prefix, grp, ": ", tp1, " vs ", tp2)
      )

      if (nrow(result) > 0) {
        result_std <- standardize_result(result, treatment_value = grp,
                                         comparison_label = paste0(tp1, " vs ", tp2),
                                         group1 = tp1, group2 = tp2)
        results_list[[length(results_list) + 1]] <- result_std

        cat("  ", tp1, "vs", tp2, ":\n")
        cat("    Cohen's d:", result$Cohens_d, "\n")
        cat("    Required n:", result$Required_n, "\n")
        cat("    Achieved power:", result$Achieved_power, "\n")
        cat("    Adequate:", result$Adequate, "\n\n")
      }
    }
  }

  do.call(rbind, results_list)
}

# -----------------------------------------------------------------------
# Combines the treatment-comparison and timepoint-comparison result tables
# into one long-format table, writes it to Excel alongside a Parameters
# sheet documenting the analysis settings and package versions, and
# returns the combined table invisibly.
# -----------------------------------------------------------------------
save_power_results <- function(treatment_result, treatment_group1, treatment_group2,
                               timepoint_result, output_path, measurement_col,
                               data_source, desired_power, alpha_level,
                               n_total_comparisons, alpha_corrected,
                               comparison1_desc, comparison2_desc) {

  treatment_std <- if (nrow(treatment_result) > 0) {
    standardize_result(
      treatment_result, treatment_value = treatment_result$timepoint,
      comparison_label = paste0(treatment_group2, " vs ", treatment_group1),
      group1 = treatment_group1, group2 = treatment_group2
    )
  } else treatment_result

  results_combined <- normalize_group_cols(rbind(treatment_std, timepoint_result))

  params_df <- data.frame(
    Parameter = c(
      "Analysis Date", "Data Source", "Measurement Column", "Desired Power",
      "Original Alpha Level", "Multiple Testing Correction", "Total Number of Comparisons",
      "Bonferroni-Corrected Alpha", "",
      "Comparison 1", "Comparison 2", "",
      "R Version", "Package: pwr", "Package: lsr"
    ),
    Value = c(
      as.character(Sys.time()), data_source, measurement_col, desired_power, alpha_level,
      "Bonferroni", n_total_comparisons, round(alpha_corrected, 5), "",
      comparison1_desc, comparison2_desc, "",
      R.version.string, as.character(packageVersion("pwr")), as.character(packageVersion("lsr"))
    )
  )

  writexl::write_xlsx(list(Results = results_combined, Parameters = params_df), output_path)
  cat("\nResults saved to:", output_path, "\n")

  invisible(results_combined)
}
