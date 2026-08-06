# =============================================================================
# helper_functions.R
#
# Shared low-level utilities used across the evidence-based model selection
# engine.
# =============================================================================

.has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.bt <- function(x) paste0("`", x, "`")

format_pval <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  as.character(signif(p, 3))
}

signif_stars <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  if (p < 0.1) return(".")
  ""
}

clean_dataset <- function(datamatrix, factor_cols = NULL) {
  names(datamatrix) <- trimws(names(datamatrix))
  df <- datamatrix
  if (is.null(factor_cols)) factor_cols <- character(0)
  for (col in factor_cols) {
    if (col %in% names(df)) df[[col]] <- factor(trimws(as.character(df[[col]])))
  }
  df
}

# Human-readable label for a random-effects structure.
describe_re_structure <- function(random_slopes, correlated_slopes = TRUE) {
  rs <- isTRUE(random_slopes); cs <- isTRUE(correlated_slopes)
  if (!rs) {
    "Fixed slopes (random intercept only)"
  } else if (cs) {
    "Random slopes, correlated with intercept (fully random)"
  } else {
    "Random slopes, decorrelated from intercept (slope levels within the block remain mutually correlated)"
  }
}

# -----------------------------------------------------------------------
# Formula building. The random intercept is always included.
# -----------------------------------------------------------------------
.build_re_formula <- function(confounder_vars, comparison_parameter,
                              random_slopes = FALSE, correlated_slopes = TRUE) {
  if (length(confounder_vars) == 0) {
    stop(".build_re_formula(): confounder_vars is empty - a mixed model needs at least one grouping variable.")
  }
  if (!random_slopes) {
    return(paste(paste0("(1 | ", .bt(confounder_vars), ")"), collapse = " + "))
  }
  if (correlated_slopes) {
    re_term <- paste0("1 + ", .bt(comparison_parameter))
    return(paste(paste0("(", re_term, " | ", .bt(confounder_vars), ")"), collapse = " + "))
  }
  paste(
    paste0("(1 | ", .bt(confounder_vars), ")"),
    paste0("(0 + ", .bt(comparison_parameter), " | ", .bt(confounder_vars), ")"),
    collapse = " + ", sep = " + "
  )
}

build_model_formula <- function(measurement, comparison_parameter, confounder_vars,
                                random_slopes = FALSE, correlated_slopes = TRUE,
                                log_transform = FALSE) {
  response_term <- if (log_transform) paste0("log(", .bt(measurement), ")") else .bt(measurement)
  re_part <- .build_re_formula(confounder_vars, comparison_parameter, random_slopes, correlated_slopes)
  paste0(response_term, " ~ ", .bt(comparison_parameter), " + ", re_part)
}

# -----------------------------------------------------------------------
# Post-fit check: are there enough grouping levels to support the number of
# variance/covariance parameters a random-slope structure estimates?
# Read-only - contributes to hard-exclusion of random-slope candidates,
# never refits anything.
# -----------------------------------------------------------------------
check_re_structure_adequacy <- function(model, confounder_vars) {
  if (is.null(model) || !inherits(model, "merMod")) return(NA_character_)
  n_theta <- tryCatch(length(lme4::getME(model, "theta")), error = function(e) NA_integer_)
  if (is.na(n_theta) || n_theta <= length(confounder_vars)) return(NA_character_)

  grp_counts <- tryCatch({
    ng <- lme4::ngrps(model)
    ng[confounder_vars[confounder_vars %in% names(ng)]]
  }, error = function(e) NULL)
  if (is.null(grp_counts) || length(grp_counts) == 0) return(NA_character_)

  min_grp <- min(grp_counts, na.rm = TRUE)
  min_grp_name <- names(grp_counts)[which.min(grp_counts)]

  if (min_grp < n_theta + 2) {
    paste0(
      "Random-slope structure estimates ", n_theta, " variance/covariance parameter(s) ",
      "from only ", min_grp, " level(s) of '", min_grp_name, "'. Likely underpowered/unstable ",
      "(rule of thumb: need noticeably more groups than variance/covariance parameters)."
    )
  } else {
    NA_character_
  }
}
