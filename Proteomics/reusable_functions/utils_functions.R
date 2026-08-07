# =============================================================================
# Mercader Lab — Data Export & Set Analysis Functions
# utils_functions.R
#
# Functions:
#   - make_gct()           : Create and write GCT files from a DESeqDataSet
#   - process_venn_sets()  : Venn diagram + intersection CSV export
#
# Usage: source("R/utils_functions.R")
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggVennDiagram)
library(cmapR)   # for GCT, write_gct()
library(DESeq2)  # for vst(), assay()


# -----------------------------------------------------------------------------
# make_gct()
#
# Creates a GCT object from a DESeqDataSet and writes output files.
# Call twice (once with use_vst = FALSE, once with TRUE) to get both
# log-counts and VST-normalised versions — do not duplicate the code blocks.
#
# Args:
#   dds        : DESeqDataSet (batch-corrected if needed, e.g. dds_combat)
#   metadata   : data frame containing sample metadata
#   sample_col : column in metadata used as rownames for cdesc (default "Sample_Name")
#   output_dir : directory for all output files (created if absent)
#   prefix     : base name for all output files e.g. "MyProject_RNAseq"
#   use_vst    : FALSE → log(counts + 1); TRUE → VST-normalised matrix
#
# Output files written to output_dir/:
#   <prefix>_log_counts.gct / <prefix>_vst_counts.gct
#   <prefix>_log_counts.csv / <prefix>_vst_counts.csv
#   <prefix>_metadata.csv
#   <prefix>_dds.rds
#
# Returns: GCT object (invisibly)
#
# Example:
#   make_gct(dds_combat, metadata, output_dir = "MyProject", prefix = "MyProject")
#   make_gct(dds_combat, metadata, output_dir = "MyProject", prefix = "MyProject", use_vst = TRUE)
# -----------------------------------------------------------------------------
make_gct <- function(dds, metadata, sample_col = "Sample_Name",
                     output_dir = ".", prefix = "project",
                     use_vst = FALSE) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  mat    <- if (use_vst) assay(vst(dds)) else log(assay(dds) + 1)
  suffix <- if (use_vst) "vst_counts" else "log_counts"

  cdesc <- metadata
  rownames(cdesc) <- cdesc[[sample_col]]
  cdesc <- cdesc %>% mutate(across(-all_of(sample_col), as.character))

  gct_obj <- new("GCT", mat = mat, cdesc = cdesc)

  gct_path  <- file.path(output_dir, paste0(prefix, "_", suffix, ".gct"))
  csv_path  <- file.path(output_dir, paste0(prefix, "_", suffix, ".csv"))
  meta_path <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  rds_path  <- file.path(output_dir, paste0(prefix, "_dds.rds"))

  write_gct(gct_obj, gct_path)
  write.csv(mat,     csv_path,  row.names = TRUE)
  write.csv(cdesc,   meta_path, row.names = FALSE)
  saveRDS(dds,       rds_path)

  cat("✓ Written:", gct_path, "\n")
  cat("✓ Written:", csv_path, "\n")
  cat("✓ Written:", meta_path, "\n")
  cat("✓ Saved dds:", rds_path, "\n")

  invisible(gct_obj)
}


# -----------------------------------------------------------------------------
# process_venn_sets()
#
# Computes all pairwise and higher-order intersections for a list of gene sets,
# exports three CSV formats, and optionally saves a Venn diagram.
#
# Args:
#   sets_list    : named list of character vectors (gene sets to compare)
#   friendly_name: string appended to all output filenames
#   output_dir   : directory for output files (default "Venn_diagrams")
#   make_plot    : whether to save PDF + PNG Venn diagram (default TRUE)
#
# Output files written to output_dir/:
#   venn_regions_summary_<friendly_name>.csv  : one row per intersection combo
#   venn_regions_long_<friendly_name>.csv     : one row per gene per combo
#   venn_regions_wide_<friendly_name>.csv     : genes as columns per combo
#   venn_plot_<friendly_name>.pdf / .png      : Venn diagram (if make_plot = TRUE)
#
# Returns: invisible list with intersections_summary, intersections_long,
#          intersections_wide
# -----------------------------------------------------------------------------
process_venn_sets <- function(sets_list, friendly_name,
                               output_dir = "Venn_diagrams", make_plot = TRUE) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(names(sets_list))) names(sets_list) <- paste0("Set", seq_along(sets_list))

  set_names <- names(sets_list)

  combos <- unlist(lapply(seq_along(set_names), function(k) {
    combn(set_names, k, simplify = FALSE)
  }), recursive = FALSE)

  intersections <- purrr::map_dfr(combos, function(x) {

    # intersection of selected sets
    inter <- Reduce(intersect, sets_list[x])

    # remove genes present in any non-selected sets
    other_sets <- setdiff(set_names, x)

    if (length(other_sets) > 0) {
      others <- unique(unlist(sets_list[other_sets]))
      inter  <- setdiff(inter, others)
    }

    tibble::tibble(
      combo = paste(x, collapse = "&"),
      genes = list(inter),
      n = length(inter)
    )
  })

  intersections_long <- intersections %>% tidyr::unnest(genes)

  intersections_csv <- intersections %>%
    dplyr::mutate(genes = sapply(genes, function(x) paste(x, collapse = ", ")))

  max_len <- max(intersections$n)

  pad_to_max <- function(x, n) { length(x) <- n; x }

  intersections_wide <- intersections %>%
    dplyr::mutate(genes = purrr::map(genes, ~ pad_to_max(.x, max_len))) %>%
    dplyr::select(combo, genes) %>%
    tidyr::unnest_wider(genes, names_sep = "_")

  intersections_wide <- as.data.frame(intersections_wide)
  rownames(intersections_wide) <- intersections_wide$combo

  intersections_wide <- as.data.frame(t(intersections_wide[, -1]))
  colnames(intersections_wide) <- intersections$combo
  rownames(intersections_wide) <- NULL

  summary_path <- file.path(output_dir, paste0("venn_regions_summary_", friendly_name, ".csv"))
  long_path    <- file.path(output_dir, paste0("venn_regions_long_",    friendly_name, ".csv"))
  wide_path    <- file.path(output_dir, paste0("venn_regions_wide_",    friendly_name, ".csv"))

  write.csv(intersections_csv,  summary_path, row.names = FALSE)
  write.csv(intersections_long, long_path,    row.names = FALSE)
  write.csv(intersections_wide, wide_path,    row.names = FALSE)

  cat("✓ Written:", summary_path, "\n")
  cat("✓ Written:", long_path,    "\n")
  cat("✓ Written:", wide_path,    "\n")

  if (make_plot) {
    pdf_path <- file.path(output_dir, paste0("venn_plot_", friendly_name, ".pdf"))
    png_path <- file.path(output_dir, paste0("venn_plot_", friendly_name, ".png"))

    p <- ggVennDiagram::ggVennDiagram(sets_list) +
      ggplot2::scale_fill_gradient(low = "white", high = "dodgerblue") +
      ggplot2::ggtitle(paste("Venn Diagram -", friendly_name))

    ggplot2::ggsave(pdf_path, p, width = 6, height = 6)
    ggplot2::ggsave(png_path, p, width = 6, height = 6, dpi = 300)

    cat("✓ Written:", pdf_path, "\n")
    cat("✓ Written:", png_path, "\n")
  }

  cat("✓ Processed Venn sets for:", friendly_name, "\n")
  
  invisible(list(
    intersections_summary = intersections_csv,
    intersections_long    = intersections_long,
    intersections_wide    = intersections_wide
  ))

}

cat("✓ Loaded utils_functions.R\n")
