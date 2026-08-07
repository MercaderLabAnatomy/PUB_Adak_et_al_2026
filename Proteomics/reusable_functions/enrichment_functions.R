# =============================================================================
# Mercader Lab — Pathway Enrichment Functions
# enrichment_functions.R
#
# Functions:
#   - ORA_function()     : Single ORA (GO / Reactome) enrichment call
#   - GSEA_function()    : Single GSEA call (GO, Reactome, or MSigDB)
#   - run_ORA_loop()     : Loop ORA_function() over all contrasts × thresholds
#   - run_GSEA_loop()    : Loop GSEA_function() over all contrasts × thresholds
#
# Output folder structure:
#   pathway_analysis/
#     ORA/
#       lfc2/
#         TreatmentVsControl/
#           GO/
#             BP/  CC/  MF/
#           Reactome/
#     GSEA/
#       lfc2/
#         TreatmentVsControl/
#           GO/
#             BP/  CC/  MF/
#           Reactome/
#
# Usage: source("reusable_functions/enrichment_functions.R")
# =============================================================================

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


# -----------------------------------------------------------------------------
# .make_enrichment_dir()  [internal helper]
#
# Builds and creates the output directory path for a given analysis type.
#
# Args:
#   analysis_type : "ORA" or "GSEA"
#   lfc_label     : threshold key e.g. "lfc2", "lfc1", "lfc0.5"
#   contrast      : contrast name e.g. "TreatmentVsControl"
#   onto          : "BP", "MF", "CC", "reactome", or "msigdb"
#
# Returns: directory path string (also creates the directory)
# -----------------------------------------------------------------------------
.make_enrichment_dir <- function(analysis_type, lfc_label, contrast, onto) {
  onto_folder <- switch(onto,
    "BP"       = "GO/BP",
    "MF"       = "GO/MF",
    "CC"       = "GO/CC",
    "reactome" = "Reactome",
    "msigdb"   = "MSigDB",
    onto
  )
  path <- file.path("pathway_analysis", analysis_type, lfc_label, contrast, onto_folder)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  return(path)
}


# -----------------------------------------------------------------------------
# ORA_function()
#
# Runs GO (BP/MF/CC) or Reactome overrepresentation analysis on a gene list.
# Auto-detects ENSEMBL vs SYMBOL input. Saves CSV + dotplot PDF.
#
# Args:
#   gene_list    : character vector of gene identifiers
#   pval         : adjusted p-value cutoff (default 0.05)
#   onto         : "BP", "MF", "CC", or "reactome"
#   contrast     : contrast name e.g. "TreatmentVsControl" — used in folder path
#   lfc_label    : threshold key e.g. "lfc2" — used in folder path
#   gene_set     : gene set direction e.g. "up", "down", "all" — used in filename
#   org          : "human" (default) or "mouse"
#   num_pathways : max pathways to plot (default 15)
#
# Output files (in ORA/<lfc_label>/<contrast>/GO/<onto>/):
#   <contrast>_<lfc_label>_<gene_set>_ORA_<onto>.csv
#   <contrast>_<lfc_label>_<gene_set>_ORA_<onto>.pdf
#
# Returns: enrichResult object or NULL
# -----------------------------------------------------------------------------
ORA_function <- function(gene_list, pval = 0.05, onto = "MF",
                         contrast = "contrast", lfc_label = "lfc1",
                         gene_set = "all", org = "human", num_pathways = 15) {

  if (length(gene_list) == 0) {
    cat("⚠ Warning: Empty gene list provided to ORA_function\n"); return(NULL)
  }

  orgdb      <- if (org == "human") org.Hs.eg.db else org.Mm.eg.db
  is_ensembl <- grepl("^ENSG\\d+", gene_list)
  key_type   <- if (sum(is_ensembl) > length(gene_list) * 0.5) "ENSEMBL" else "SYMBOL"
  compGO     <- NULL

  if (onto %in% c("MF", "CC", "BP")) {
    compGO <- tryCatch(
      enrichGO(gene = gene_list, pvalueCutoff = pval, keyType = key_type,
               pAdjustMethod = "BH", OrgDb = orgdb, ont = onto),
      error = function(e) { cat("✗ Error in ORA enrichment:", e$message, "\n"); NULL })

  } else if (onto == "reactome") {
    gene_list_df <- tryCatch(
      bitr(gene_list, fromType = key_type, toType = "ENTREZID", OrgDb = orgdb),
      error = function(e) { cat("✗ Error in gene ID conversion:", e$message, "\n"); NULL })

    if (!is.null(gene_list_df) && nrow(gene_list_df) > 0) {
      compGO <- tryCatch(
        enrichPathway(gene = gene_list_df$ENTREZID, pvalueCutoff = pval,
                      organism = if (org == "human") "human" else "mouse",
                      readable = TRUE),
        error = function(e) { cat("✗ Error in Reactome ORA:", e$message, "\n"); NULL })
    }
  }

  if (!is.null(compGO) && nrow(as.data.frame(compGO)) > 0) {
    compGO_df <- as.data.frame(compGO)
    compGO_df$GeneRatio_decimal <- sapply(compGO_df$GeneRatio,
                                          function(x) eval(parse(text = x)))

    out_dir  <- .make_enrichment_dir("ORA", lfc_label, contrast, onto)
    prefix   <- paste0(contrast, "_", lfc_label, "_", gene_set)

    write.csv(compGO_df,
              file.path(out_dir, paste0(prefix, "_ORA_", onto, ".csv")),
              row.names = FALSE)

    pathways_selected <- unique(compGO_df$Description)[1:min(num_pathways, nrow(compGO_df))]
    compGO_df_subset  <- compGO_df[compGO_df$Description %in% pathways_selected, ]

    p <- ggplot(compGO_df_subset,
                aes(x = GeneRatio_decimal,
                    y = reorder(Description, GeneRatio_decimal),
                    color = p.adjust, size = Count)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue",
                           guide = guide_colourbar(reverse = TRUE)) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9)) +
      labs(title    = paste0("ORA: ", onto),
           subtitle = paste(contrast, lfc_label, gene_set),
           y = "Pathway", x = "Gene Ratio",
           color = "Adj. P-value", size = "Count")

    ggsave(file.path(out_dir, paste0(prefix, "_ORA_", onto, ".pdf")),
           p, width = 12, height = max(6, num_pathways * 0.4))

    cat("✓ ORA completed:", onto, "-", nrow(compGO_df), "pathways,",
        nrow(compGO_df_subset), "plotted →", out_dir, "\n")
    return(compGO)
  } else {
    cat("✗ No ORA enrichment found for", onto, "\n"); return(NULL)
  }
}


# -----------------------------------------------------------------------------
# GSEA_function()
#
# Runs Gene Set Enrichment Analysis using a ranked gene list (by log2FC).
# Supports GO (BP/MF/CC), Reactome, and MSigDB collections.
# Saves CSV, dotplot PDF, and ranked enrichment PDFs for top pathways.
#
# Args:
#   df                : DESeq2 results data frame; requires columns:
#                       padj, log2FoldChange, gene
#   pval_deg          : padj cutoff for building the ranked list (default 0.05)
#   pval_enrich       : enrichment p-value cutoff (default 0.2)
#   onto              : "BP", "MF", "CC", "reactome", or "msigdb"
#   contrast          : contrast name e.g. "TreatmentVsControl" — used in folder path
#   lfc_label         : threshold key e.g. "lfc2" — used in folder path
#   organism          : "human" (default) or "mouse"
#   msigdb_collection : MSigDB category string e.g. "H" (Hallmarks);
#                       only used when onto = "msigdb"
#   num_pathways      : max pathways in dotplot (default 15)
#   top_n_ranked      : number of individual enrichment plots to save (default 3)
#
# Output files (in GSEA/<lfc_label>/<contrast>/GO/<onto>/):
#   <contrast>_<lfc_label>_GSEA_<onto>.csv
#   <contrast>_<lfc_label>_GSEA_<onto>_dotplot.pdf
#   <contrast>_<lfc_label>_GSEA_<onto>_ranked_<term>.pdf
#
# Returns: gseaResult object or NULL
# -----------------------------------------------------------------------------

# Internal helper: replace ENTREZID with SYMBOL in core_enrichment column
.entrezid_to_symbol <- function(result_df, OrgDb) {
  if (!"core_enrichment" %in% names(result_df)) return(result_df)
  all_ids <- unique(unlist(strsplit(as.character(result_df$core_enrichment), "/")))
  id_map  <- tryCatch(
    bitr(all_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb),
    error = function(e) NULL)
  if (is.null(id_map) || nrow(id_map) == 0) return(result_df)
  result_df$core_enrichment <- sapply(result_df$core_enrichment, function(ids_str) {
    ids  <- strsplit(ids_str, "/")[[1]]
    syms <- id_map$SYMBOL[match(ids, id_map$ENTREZID)]
    syms[is.na(syms)] <- ids[is.na(syms)]
    paste(syms, collapse = "/")
  })
  return(result_df)
}

GSEA_function <- function(df, pval_deg = 0.05, pval_enrich = 0.2, onto = "MF",
                          contrast = "contrast", lfc_label = "lfc1",
                          lfc_threshold = 0,
                          organism = "human", msigdb_collection = NULL,
                          num_pathways = 15, top_n_ranked = 3) {

  gene_list_df <- df %>%
    filter(padj <= pval_deg, !is.na(gene)) %>%
    filter(abs(log2FoldChange) >= lfc_threshold) %>%
    arrange(desc(log2FoldChange)) %>%
    filter(!duplicated(gene))
  gene_list_df %>% summarise(total_genes = n(), upregulated = sum(log2FoldChange > 0),
                 downregulated = sum(log2FoldChange < 0)) %>%
    print()
  if (nrow(gene_list_df) == 0) {
    cat("⚠ No significant genes found for GSEA\n"); return(NULL)
  }
  gene_list <- sort(setNames(gene_list_df$log2FoldChange, gene_list_df$gene),
                   decreasing = TRUE)
  cat("  Ranked gene list:", length(gene_list), "genes\n")

  # Auto-detect input ID type (same logic as ORA_function)
  is_ensembl <- grepl("^ENSG\\d+", names(gene_list))
  key_type   <- if (sum(is_ensembl) > length(gene_list) * 0.5) "ENSEMBL" else "SYMBOL"
  cat("  Detected ID type:", key_type, "\n")

  OrgDb  <- if (organism == "human") org.Hs.eg.db else org.Mm.eg.db
  compGO <- NULL

  if (onto %in% c("MF", "CC", "BP")) {
    # Convert to ENTREZID — more robust than SYMBOL keyType in gseGO
    gene_entrez_go <- tryCatch(
      bitr(names(gene_list), fromType = key_type, toType = "ENTREZID", OrgDb = OrgDb),
      error = function(e) { cat("  ✗ Error in", key_type, "->ENTREZID:", e$message, "\n"); NULL })

    if (is.null(gene_entrez_go) || nrow(gene_entrez_go) == 0) {
      cat("  ✗ No genes mapped to ENTREZID for", onto, "\n"); return(NULL)
    }
    gene_list_ez <- gene_list[names(gene_list) %in% gene_entrez_go[[key_type]]]
    names(gene_list_ez) <- gene_entrez_go$ENTREZID[
      match(names(gene_list_ez), gene_entrez_go[[key_type]])]
    gene_list_ez <- sort(gene_list_ez, decreasing = TRUE)
    gene_list_ez <- gene_list_ez[!duplicated(names(gene_list_ez))]
    cat("  Mapped:", length(gene_list_ez), "genes to ENTREZID\n")

    compGO <- tryCatch(
      gseGO(geneList = gene_list_ez, pvalueCutoff = pval_enrich, keyType = "ENTREZID",
            pAdjustMethod = "BH", OrgDb = OrgDb, ont = onto),
      error = function(e) { cat("  ✗ Error in gseGO:", e$message, "\n"); NULL })

  } else if (onto == "reactome") {
    gene_entrez <- tryCatch(
      bitr(names(gene_list), fromType = key_type, toType = "ENTREZID", OrgDb = OrgDb),
      error = function(e) { cat("  ✗ Error in gene ID conversion:", e$message, "\n"); NULL })

    if (!is.null(gene_entrez) && nrow(gene_entrez) > 0) {
      gene_list_entrez <- gene_list[names(gene_list) %in% gene_entrez[[key_type]]]
      names(gene_list_entrez) <- gene_entrez$ENTREZID[
        match(names(gene_list_entrez), gene_entrez[[key_type]])]
      gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
      gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
      compGO <- tryCatch(
        gsePathway(geneList = gene_list_entrez, pvalueCutoff = pval_enrich,
                   organism = organism),
        error = function(e) { cat("  ✗ Error in gsePathway:", e$message, "\n"); NULL })
    }

  } else if (onto == "msigdb" && !is.null(msigdb_collection)) {
    # MSigDB TERM2GENE uses gene symbols — convert ENSEMBL to SYMBOL if needed
    gene_list_sym <- gene_list
    if (key_type == "ENSEMBL") {
      sym_map <- tryCatch(
        bitr(names(gene_list), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = OrgDb),
        error = function(e) { cat("  ✗ Error in ENSEMBL->SYMBOL:", e$message, "\n"); NULL })
      if (is.null(sym_map) || nrow(sym_map) == 0) {
        cat("  ✗ Failed to convert ENSEMBL to SYMBOL for MSigDB\n"); return(NULL)
      }
      gene_list_sym <- gene_list[names(gene_list) %in% sym_map$ENSEMBL]
      names(gene_list_sym) <- sym_map$SYMBOL[match(names(gene_list_sym), sym_map$ENSEMBL)]
      gene_list_sym <- sort(gene_list_sym, decreasing = TRUE)
      gene_list_sym <- gene_list_sym[!duplicated(names(gene_list_sym))]
    }
    msigdbr_df       <- msigdbr(species = if (organism == "human") "Homo sapiens" else "Mus musculus")
    cat_col          <- if ("gs_cat" %in% names(msigdbr_df)) "gs_cat" else "gs_collection"
    msig_df_filtered <- msigdbr_df %>% filter(.data[[cat_col]] == msigdb_collection)
    term2gene        <- msig_df_filtered %>% dplyr::select(gs_name, gene_symbol)
    compGO <- tryCatch(
      GSEA(geneList = gene_list_sym, TERM2GENE = term2gene,
           pvalueCutoff = pval_enrich, pAdjustMethod = "BH"),
      error = function(e) { cat("  ✗ Error in GSEA:", e$message, "\n"); NULL })
  }

  if (!is.null(compGO) && nrow(as.data.frame(compGO)) > 0) {
    compGO@result$.sign <- ifelse(compGO@result$NES > 0, "Activated", "Suppressed")

    out_dir <- .make_enrichment_dir("GSEA", lfc_label, contrast, onto)
    prefix  <- paste0(contrast, "_", lfc_label)

    # Convert ENTREZID in core_enrichment back to gene symbols for readability
    if (onto %in% c("BP", "MF", "CC", "reactome")) {
      compGO@result <- .entrezid_to_symbol(compGO@result, OrgDb)
    }

    write.csv(compGO@result,
              file.path(out_dir, paste0(prefix, "_GSEA_", onto, ".csv")),
              row.names = FALSE)

    compGO_df_subset <- compGO@result %>%
      slice_max(order_by = abs(NES), n = num_pathways, with_ties = FALSE)

    if (nrow(compGO_df_subset) > 0) {
      p <- ggplot(compGO_df_subset,
                  aes(x = NES, y = reorder(Description, NES),
                      color = p.adjust, size = setSize)) +
        geom_point(alpha = 0.7) +
        scale_color_gradient(low = "red", high = "blue") +
        facet_grid(. ~ .sign, scales = "free_y", space = "free_y") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 8)) +
        labs(title    = paste0("GSEA: ", onto),
             subtitle = paste(contrast, lfc_label),
             y = "Pathway", x = "NES",
             color = "Adj. P-value", size = "Set Size")
      ggsave(file.path(out_dir, paste0(prefix, "_GSEA_", onto, "_dotplot.pdf")),
             p, width = 12, height = max(6, num_pathways * 0.35))
    }

    for (term in head(compGO@result$ID, top_n_ranked)) {
      term_name <- gsub("[^A-Za-z0-9_-]", "_", term)
      pdf_file  <- file.path(out_dir, paste0(prefix, "_GSEA_", onto,
                                             "_ranked_", term_name, ".pdf"))
      tryCatch({
        pdf(pdf_file, width = 10, height = 6)
        print(enrichplot::gseaplot2(compGO, geneSetID = term, title = term))
        dev.off()
      }, error = function(e) cat("  ⚠ Could not create enrichment plot for", term, "\n"))
    }

    cat("  ✓ GSEA completed:", onto, "-", nrow(compGO@result), "pathways,",
        nrow(compGO_df_subset), "plotted →", out_dir, "\n")
    return(compGO)
  } else {
    cat("  ✗ No GSEA enrichment found for", onto, "\n"); return(NULL)
  }
}


# -----------------------------------------------------------------------------
# run_ORA_loop()
#
# Iterates ORA_function() over all contrasts × thresholds × ontologies × gene sets.
#
# Args:
#   filtered_results : named list [[contrast]][[thresh_name]], each a DEG data
#                      frame with columns gene, log2FoldChange
#   contrasts        : character vector of contrast names
#   thresholds       : named list keyed by LFC label e.g. list(lfc0.5=0.5, lfc1=1, lfc2=2)
#   org              : "human" (default) or "mouse"
#   num_pathways     : passed through to ORA_function() (default 15)
# -----------------------------------------------------------------------------
run_ORA_loop <- function(filtered_results, contrasts, thresholds,
                         org = "human", num_pathways = 15) {
  ontologies <- c("BP", "MF", "CC", "reactome")
  gene_sets  <- c("all", "up", "down")
  total <- length(contrasts) * length(ontologies) * length(gene_sets) * length(thresholds)

  cat("Starting ORA enrichment analysis...\n")
  cat("Total analyses to run:", length(contrasts), "×", length(ontologies), "×",
      length(gene_sets), "×", length(thresholds), "=", total, "\n\n")

  for (contrast in contrasts) {
    cat("\n=== Processing contrast:", contrast, "===\n")

    for (thresh_name in names(thresholds)) {
      cat("\n--- Threshold:", thresh_name, "---\n")
      deg_data <- filtered_results[[contrast]][[thresh_name]]

      if (is.null(deg_data) || nrow(deg_data) == 0) {
        cat("⚠ No DEGs found for", contrast, thresh_name, "- skipping\n"); next
      }

      gene_sets_list <- list(
        all  = deg_data$gene,
        up   = deg_data$gene[deg_data$log2FoldChange >  0],
        down = deg_data$gene[deg_data$log2FoldChange <  0]
      )

      for (onto in ontologies) {
        for (set_name in gene_sets) {
          genes <- gene_sets_list[[set_name]]
          if (length(genes) < 3) {
            cat("  ⊘", onto, "-", set_name, ": Too few genes (", length(genes), ")\n"); next
          }
          cat("  →", onto, "-", set_name, ":", length(genes), "genes... ")
          ORA_function(gene_list = genes, pval = 0.05, onto = onto,
                       contrast = contrast, lfc_label = thresh_name,
                       gene_set = set_name, org = org, num_pathways = num_pathways)
        }
      }
    }
    cat("\n✓ Completed all ORA analyses for", contrast, "\n")
  }
  cat("\n✓✓✓ All ORA enrichment analyses completed ✓✓✓\n")
}


# -----------------------------------------------------------------------------
# run_GSEA_loop()
#
# Iterates GSEA_function() over all contrasts × thresholds × analysis types.
# Analysis types: BP, MF, CC, Reactome, MSigDB Hallmarks.
#
# Args:
#   results_list  : named list [[contrast]] — each element the full (unfiltered)
#                   DESeq2 results data frame with columns padj, log2FoldChange, gene
#   contrasts     : character vector of contrast names
#   thresholds    : named list keyed by LFC label — same convention as run_ORA_loop()
#   organism      : "human" (default) or "mouse"
#   num_pathways  : passed through to GSEA_function() (default 15)
#   top_n_ranked  : number of individual enrichment plots per run (default 3)
# -----------------------------------------------------------------------------
run_GSEA_loop <- function(results_list, contrasts, thresholds,
                          organism = "human", num_pathways = 15, top_n_ranked = 3,
                          pval_deg = 1.0) {
  gsea_analyses <- list(
    list(onto = "BP",       msigdb = NULL),
    list(onto = "MF",       msigdb = NULL),
    list(onto = "CC",       msigdb = NULL),
    list(onto = "reactome", msigdb = NULL),
    list(onto = "msigdb",   msigdb = "H")
  )

  total <- length(contrasts) * length(gsea_analyses) * length(thresholds)
  cat("Starting GSEA analysis...\n")
  cat("Total analyses:", length(contrasts), "×", length(gsea_analyses), "×",
      length(thresholds), "=", total, "\n\n")

  for (contrast in contrasts) {
    cat("\n=== GSEA for contrast:", contrast, "===\n")
    df_full <- results_list[[contrast]]

    for (thresh_name in names(thresholds)) {
      cat("\n--- Threshold:", thresh_name, "---\n")

      lfc_thresh_val <- thresholds[[thresh_name]]
      for (analysis in gsea_analyses) {
        cat("  →", analysis$onto, "... ")
        GSEA_function(df                = df_full,
                      pval_deg          = pval_deg,
                      pval_enrich       = 0.2,
                      onto              = analysis$onto,
                      contrast          = contrast,
                      lfc_label         = thresh_name,
                      lfc_threshold     = lfc_thresh_val,
                      organism          = organism,
                      msigdb_collection = analysis$msigdb,
                      num_pathways      = num_pathways,
                      top_n_ranked      = top_n_ranked)
      }
    }
    cat("\n✓ Completed all GSEA analyses for", contrast, "\n")
  }
  cat("\n✓✓✓ All GSEA analyses completed ✓✓✓\n")
}

cat("✓ enrichment_functions.R loaded: ORA_function(), GSEA_function(), run_ORA_loop(), run_GSEA_loop()\n")
