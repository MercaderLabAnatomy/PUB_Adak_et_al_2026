# =============================================================================
# Mercader Lab — Source All Reusable Functions
# source_all.R
#
# Sources all shared function files.
# Call at the top of every analysis script / notebook cell.
#
# Usage: source("reusable_functions/source_all.R")
# =============================================================================

base_dir <- "reusable_functions"

source(file.path(base_dir, "viz_functions.R"))
source(file.path(base_dir, "enrichment_functions.R"))
source(file.path(base_dir, "utils_functions.R"))
source(file.path(base_dir, "orthology_functions.R"))
source(file.path(base_dir, "sc_functions.R"))

cat("✓ All Mercader Lab reusable functions loaded.\n")
cat("  Bulk:        PCA(), draw_volcano(), make_gct(), process_venn_sets()\n")
cat("  Enrichment:  ORA_function(), GSEA_function(), run_ORA_loop(), run_GSEA_loop()\n")
cat("  Orthology:   convert_zf_to_mouse_alliance(), convert_zf_to_human_alliance()\n")
cat("  Single-cell: run_liana_zebrafish(), run_liana_loop(), plot_liana_dotplot(),\n")
cat("               pseudobulk_de(), run_ORA_loop_sc(), run_GSEA_loop_sc()\n")
