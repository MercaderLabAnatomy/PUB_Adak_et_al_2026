# =============================================================================
# Mercader Lab — Visualization Functions
# viz_functions.R
#
# Functions:
#   - PCA()          : PCA plot from a count/expression matrix
#   - draw_volcano() : Volcano plot from DESeq2-style results
#
# Usage: source("R/viz_functions.R")
# =============================================================================

library(ggplot2)
library(ggrepel)

options(ggrepel.max.overlaps = 20)

# -----------------------------------------------------------------------------
# PCA()
#
# Args:
#   mat          : numeric matrix, samples as columns
#   color_pca    : vector mapping samples to color groups
#   shape_pca    : vector mapping samples to shape groups
#   label_pca    : vector of sample labels for repel text
#   save_plot    : "yes" saves a PDF; "no" (default) prints and returns the plot
#   name_of_plot : base filename (no extension) when save_plot = "yes"
#   comp1, comp2 : which PCs to plot (default 1, 2)
#   pdf_height, pdf_width : dimensions in inches when saving
#
# Returns: ggplot object
# -----------------------------------------------------------------------------
PCA <- function(mat, color_pca = "", shape_pca = "", label_pca = "",
                save_plot = "no", name_of_plot = "PCA",
                comp1 = 1, comp2 = 2, pdf_height = 12, pdf_width = 12) {
  dt <- mat
  pca_dt <- prcomp(t(dt))
  cat("PCA running...\n")
  percentVar_dt <- pca_dt$sdev^2 / sum(pca_dt$sdev^2)
  cat("Percents calculated...\n")
  dt_f <- data.frame(
    PC1       = pca_dt$x[, comp1],
    PC2       = pca_dt$x[, comp2],
    color_pca = color_pca,
    shape_pca = shape_pca,
    label_pca = label_pca
  )
  cat("Data frame built...\n")
  cat("Plotting...\n")

  pca_p <- ggplot(data = dt_f,
                  aes_string(x = "PC1", y = "PC2",
                             color = "color_pca", shape = "shape_pca",
                             label = "label_pca")) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 50,
                    box.padding = 1.5, point.padding = 0.5, force = 50) +
    xlab(paste0("PC", comp1, ": ", round(percentVar_dt[comp1] * 100), "% variance")) +
    ylab(paste0("PC", comp2, ": ", round(percentVar_dt[comp2] * 100), "% variance")) +
    theme_bw()

  if (save_plot == "yes") {
    pdf(paste0(name_of_plot, ".pdf"), width = pdf_width, height = pdf_height)
    cat("Saving plot as:", paste0(name_of_plot, ".pdf\n"))
    print(pca_p)
    dev.off()
  }

  cat("Done\n")
  print(pca_p)
}


# -----------------------------------------------------------------------------
# draw_volcano()
#
# Args:
#   fileinput : data frame with columns log2FoldChange and padj
#   title     : string in the form "GroupA vs GroupB" — group names are parsed
#               automatically to annotate axes (red = higher in GroupA,
#               blue = higher in GroupB)
#   lfcthres  : absolute log2FC threshold for coloring points (default 1)
#
# Returns: ggplot object (does NOT save — wrap in ggsave() or pdf() at call site)
# -----------------------------------------------------------------------------
draw_volcano <- function(fileinput, title, lfcthres = 1) {
  # Clean data
  fileinput <- fileinput[is.finite(fileinput$log2FoldChange), ]
  fileinput$padj[is.na(fileinput$padj) | !is.finite(fileinput$padj) | fileinput$padj < 0] <- 1
  fileinput$padj[fileinput$padj < 1e-100] <- 1e-110
  fileinput$log2FoldChange[fileinput$log2FoldChange >  7] <-  7
  fileinput$log2FoldChange[fileinput$log2FoldChange < -7] <- -7

  xlim <- c(round(min(fileinput$log2FoldChange, na.rm = TRUE) - 1),
            round(max(fileinput$log2FoldChange, na.rm = TRUE) + 1))
  ylim <- c(-1, round(-log10(min(fileinput$padj, na.rm = TRUE)) + 1))

  p <- ggplot(data = fileinput, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", col = "darkgoldenrod") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.5, size = 2, color = "grey51") +
    geom_point(data = fileinput[fileinput$padj < 0.05 & fileinput$log2FoldChange < -lfcthres, ],
               aes(x = log2FoldChange, y = -log10(padj)),
               shape = 21, color = "#3A5FCD", fill = "#3A5FCD",
               alpha = 0.5, size = 2, inherit.aes = FALSE) +
    geom_point(data = fileinput[fileinput$padj < 0.05 & fileinput$log2FoldChange >  lfcthres, ],
               aes(x = log2FoldChange, y = -log10(padj)),
               shape = 21, color = "#CD0000", fill = "#CD0000",
               alpha = 0.5, size = 2, inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(xlim[1], xlim[2], by = 1), limits = xlim) +
    scale_y_continuous(breaks = seq(0, ylim[2], by = 4), limits = ylim) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5),
          axis.text    = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))

  # Auto-annotate axes from "GroupA vs GroupB" title
  match <- regexec("^(.*)\\s+[Vv][Ss]\\s+(.*)$", title)
  parts <- regmatches(title, match)[[1]]
  if (length(parts) == 3) {
    group1  <- trimws(parts[2])
    group2  <- trimws(parts[3])
    y_annot <- -0.75
    p <- p +
      annotate("text", x = xlim[2] - 0.5, y = y_annot,
               label = paste("Higher in", group1),
               hjust = 1, size = 3.5, fontface = "italic", color = "red3") +
      annotate("text", x = xlim[1] + 0.5, y = y_annot,
               label = paste("Higher in", group2),
               hjust = 0, size = 3.5, fontface = "italic", color = "royalblue3")
  }
  return(p)
}

cat("✓ viz_functions.R loaded: PCA(), draw_volcano()\n")
