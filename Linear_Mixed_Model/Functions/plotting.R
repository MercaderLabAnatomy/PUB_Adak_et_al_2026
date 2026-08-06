# =============================================================================
# plotting.R
#
# Diagnostic plots for whichever candidate is ultimately selected.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

plot_gaussian_diagnostics <- function(model, measurement) {
  resid_vals <- residuals(model)
  fitted_vals <- fitted(model)
  sh <- .shapiro_safe(resid_vals)

  p_qq <- ggplot(data.frame(resid = resid_vals), aes(sample = resid)) +
    stat_qq(alpha = 0.6) +
    stat_qq_line(color = "firebrick") +
    labs(title = paste0(measurement, " - Residual QQ-plot"),
         subtitle = if (!is.na(sh$p)) paste0("Shapiro-Wilk p = ", format_pval(sh$p)) else "Shapiro-Wilk not computed",
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal(base_size = 10)

  p_rvf <- ggplot(data.frame(fitted = fitted_vals, resid = resid_vals), aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, color = "firebrick", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "steelblue", linewidth = 0.6) +
    labs(title = paste0(measurement, " - Residuals vs Fitted"), x = "Fitted values", y = "Residuals") +
    theme_minimal(base_size = 10)

  wrap_plots(p_qq, p_rvf, nrow = 1) +
    plot_annotation(title = paste0("Measurement: ", measurement, " - Gaussian LMM diagnostics"),
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))
}

plot_dharma_diagnostics <- function(model, measurement, sim = NULL, family = NULL) {
  if (is.null(sim)) sim <- tryCatch(DHARMa::simulateResiduals(model, n = 1000, plot = FALSE), error = function(e) NULL)
  if (is.null(sim)) {
    return(wrap_plots(ggplot() + annotate("text", x = 1, y = 1, label = "DHARMa simulation failed") + theme_void()))
  }
  resid_scaled <- sim$scaledResiduals

  # fitted() on a Gamma GLMM (glmer, log link) already back-transforms to
  # the response scale by default. fitted() on a Log-normal LMM returns
  # values on the log scale, because the model itself was fit to log(y) -
  # exponentiate here so the x-axis is actually on the response scale.
  fitted_vals <- fitted(model)
  if (!is.null(family) && family == "lognormal") fitted_vals <- exp(fitted_vals)

  ks_p   <- tryCatch(DHARMa::testUniformity(sim, plot = FALSE)$p.value, error = function(e) NA_real_)
  disp_p <- tryCatch(DHARMa::testDispersion(sim, plot = FALSE)$p.value, error = function(e) NA_real_)

  n <- length(resid_scaled)
  p_qq <- ggplot(data.frame(resid = sort(resid_scaled), theoretical = ppoints(n)), aes(x = theoretical, y = resid)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
    labs(title = paste0(measurement, " - DHARMa scaled-residual QQ-plot"),
         subtitle = if (!is.na(ks_p)) paste0("KS test (uniformity) p = ", format_pval(ks_p)) else "KS test not computed",
         x = "Expected (Uniform quantiles)", y = "Observed (scaled residuals)") +
    theme_minimal(base_size = 10)

  p_rvf <- ggplot(data.frame(fitted = fitted_vals, resid = resid_scaled), aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), color = "firebrick", linetype = "dashed", linewidth = 0.4) +
    geom_smooth(method = "loess", se = FALSE, color = "steelblue", linewidth = 0.6) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste0(measurement, " - DHARMa residuals vs Fitted"),
         subtitle = if (!is.na(disp_p)) paste0("Dispersion test p = ", format_pval(disp_p)) else "Dispersion test not computed",
         x = "Fitted values (response scale)", y = "DHARMa scaled residual") +
    theme_minimal(base_size = 10)

  wrap_plots(p_qq, p_rvf, nrow = 1) +
    plot_annotation(title = paste0("Measurement: ", measurement, " - model diagnostics (DHARMa)"),
                    theme = theme(plot.title = element_text(size = 13, face = "bold")))
}
