# =============================================================================
# Mercader Lab — Single-Cell Functions
# sc_functions.R
#
# Functions:
#   ── Clustering quality ──────────────────────────────────────────────────────
#   compute_clustering_qc()      : Silhouette + normalized Shannon entropy
#   plot_clustering_qc()         : Two-panel QC plot (silhouette + entropy)
#
#   ── Cell type annotation ────────────────────────────────────────────────────
#   annotate_celltypes_enrichr() : Automated cell-type labeling via Enrichr
#                                  (zf→human, fuzzy consensus, XLSX output)
#
#   ── DE method selection ─────────────────────────────────────────────────────
#   compute_data_characteristics(): Sparsity & depth per cell type
#   recommend_de_method()        : Nguyen 2023 lookup → recommended DE methods
#
#   ── Nguyen 2023 DE implementations ──────────────────────────────────────────
#   run_pseudobulk_limma()       : Pseudobulk_limma  — pseudobulk + limma-trend
#   run_limmatrend_cov()         : limmatrend / limmatrend_Cov
#   run_limmavoom_cov()          : limmavoom / limmavoom_Cov
#   run_edgeR_cov()              : Pseudobulk_edgeR / ZW_edgeR_Cov
#   run_mast_cov()               : MAST / MAST_Cov
#   run_wilcox_raw()             : Raw_Wilcox
#   run_de_per_celltype()        : Dispatcher — selects & runs recommended method
#
#   ── LIANA ───────────────────────────────────────────────────────────────────
#   run_liana_zebrafish()        : Build zebrafish LIANA ortholog resource
#   run_liana_loop()             : LIANA per condition → CSVs + dotplots
#   plot_liana_dotplot()         : Save publication LIANA dotplot PDF
#
#   ── SC enrichment wrappers ──────────────────────────────────────────────────
#   run_ORA_loop_sc()            : ORA with <CellType> subfolder
#   run_GSEA_loop_sc()           : GSEA with <CellType> subfolder
#
#   ── Legacy ──────────────────────────────────────────────────────────────────
#   pseudobulk_de()              : Alias for run_pseudobulk_limma()
#
# Usage: source("reusable_functions/sc_functions.R")
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(ggplot2); library(patchwork)
})


# =============================================================================
# 1. CLUSTERING QUALITY METRICS
# =============================================================================

# -----------------------------------------------------------------------------
# compute_clustering_qc()
#
# Silhouette width (cluster compactness, PCA/Harmony space, subsampled) and
# normalized Shannon entropy (sample mixing per cluster = integration quality).
#
# Args:
#   seurat_obj    : Seurat object with Harmony reduction
#   resolutions   : numeric vector — must already be run via FindClusters()
#   reduction     : reduction embedding used for silhouette (default "harmony")
#   prefix        : cluster column prefix (default "SCT_snn_res.")
#   subsample_n   : max cells for silhouette distance calc (default 5000)
#   sample_col    : metadata column for sample ID used in entropy (default "sample")
#
# Returns: data frame — resolution, n_clusters, mean_silhouette, norm_entropy
# -----------------------------------------------------------------------------
compute_clustering_qc <- function(seurat_obj,
                                   resolutions,
                                   reduction   = "harmony",
                                   prefix      = "SCT_snn_res.",
                                   subsample_n = 5000,
                                   sample_col  = "sample") {
  if (!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")
  library(cluster)

  emb  <- Embeddings(seurat_obj, reduction)
  meta <- seurat_obj@meta.data

  results <- lapply(resolutions, function(res) {
    col <- paste0(prefix, res)
    if (!col %in% colnames(meta)) {
      cat("  ⚠ Column not found:", col, "— skipping\n"); return(NULL)
    }
    clust   <- as.integer(factor(meta[[col]]))
    n_clust <- length(unique(clust))

    # ── Silhouette (subsampled) ───────────────────────────────────────────────
    if (nrow(emb) > subsample_n) {
      idx  <- sample(nrow(emb), subsample_n)
      sil  <- silhouette(clust[idx], dist(emb[idx, ]))
    } else {
      sil  <- silhouette(clust, dist(emb))
    }
    mean_sil <- mean(sil[, 3], na.rm = TRUE)

    # ── Normalized Shannon entropy (sample mixing) ────────────────────────────
    if (!sample_col %in% colnames(meta)) {
      norm_ent <- NA_real_
    } else {
      cluster_ent <- tapply(meta[[sample_col]], clust, function(s) {
        tbl <- table(s); p <- tbl / sum(tbl); p <- p[p > 0]
        ent <- -sum(p * log2(p))
        n_s <- length(unique(s))
        if (n_s < 2) return(0)
        ent / log2(n_s)   # normalize to [0, 1]
      })
      norm_ent <- mean(unlist(cluster_ent), na.rm = TRUE)
    }

    data.frame(resolution = res, n_clusters = n_clust,
               mean_silhouette = mean_sil, norm_entropy = norm_ent)
  })
  do.call(rbind, Filter(Negate(is.null), results))
}


# -----------------------------------------------------------------------------
# plot_clustering_qc()
#
# Two-panel publication plot: silhouette width (left) + normalized entropy
# (right) across resolutions. Saves PDF + PNG.
#
# Args:
#   qc_df   : output of compute_clustering_qc()
#   out_dir : output directory (default "QC")
# -----------------------------------------------------------------------------
plot_clustering_qc <- function(qc_df, out_dir = "QC") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  p_sil <- ggplot(qc_df, aes(x = factor(resolution), y = mean_silhouette)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_point(size = 3, color = "black") +
    labs(title = "Silhouette Width", x = "Resolution", y = "Mean Silhouette") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_ent <- ggplot(qc_df, aes(x = factor(resolution), y = norm_entropy)) +
    geom_col(fill = "coral3", alpha = 0.8) +
    geom_point(size = 3, color = "black") +
    labs(title = "Normalized Sample Entropy\n(Integration Quality)",
         x = "Resolution", y = "Norm. Shannon Entropy") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p <- p_sil | p_ent
  ggsave(file.path(out_dir, "clustering_qc_silhouette_entropy.pdf"),  p, width = 12, height = 5)
  ggsave(file.path(out_dir, "clustering_qc_silhouette_entropy.png"),  p, width = 12, height = 5, dpi = 300)
  cat("✓ Clustering QC plot saved →", out_dir, "\n")
  print(qc_df)
  invisible(p)
}


# =============================================================================
# 2. AUTOMATED CELL TYPE ANNOTATION VIA ENRICHR
# =============================================================================

# -----------------------------------------------------------------------------
# annotate_celltypes_enrichr()
#
# Per-cluster automated cell-type labeling:
#   1. Top markers (avg_log2FC-ranked, FDR < p_adj_cut) selected per cluster
#   2. Zebrafish genes → human orthologs (convert_zf_to_human_alliance)
#   3. 15 Enrichr databases queried
#   4. Fuzzy consensus (Jaro-Winkler >= sim_thresh, spanning >= min_libs)
#   5. Per-cluster XLSX saved to cluster_enrichr/ (color-coded by DB)
#   6. Summary saved to marker_results/enrichr_predictions.xlsx
#
# Args:
#   all_markers   : FindAllMarkers() output — needs cluster, gene, avg_log2FC, p_val_adj
#   alliance_file : path to ORTHOLOGY-ALLIANCE_COMBINED_*.tsv
#   top_n         : max markers sent per cluster (default 250)
#   n_top_terms   : top Enrichr hits collected per DB (default 10)
#   p_adj_cut     : FDR threshold for marker selection (default 0.05)
#   min_genes     : min human orthologs required (default 5)
#   sim_thresh    : Jaro-Winkler threshold for fuzzy term grouping (default 0.85)
#   min_libs      : min libraries required for confident consensus (default 2)
#   out_dir       : directory for per-cluster XLSX (default "cluster_enrichr")
#   summary_dir   : directory for summary XLSX (default "marker_results")
#
# Returns: summary data frame (cluster, n_markers_sent, cell_type, n_libs_support)
# -----------------------------------------------------------------------------
annotate_celltypes_enrichr <- function(all_markers,
                                        alliance_file,
                                        top_n       = 250,
                                        n_top_terms = 10,
                                        p_adj_cut   = 0.05,
                                        min_genes   = 5,
                                        sim_thresh  = 0.85,
                                        min_libs    = 2,
                                        out_dir     = "cluster_enrichr",
                                        summary_dir = "marker_results") {

  for (pkg in c("enrichR", "stringdist", "openxlsx")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  suppressPackageStartupMessages({
    library(enrichR); library(stringdist); library(openxlsx)
  })

  dir.create(out_dir,     showWarnings = FALSE, recursive = TRUE)
  dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

  # ── Ortholog table ──────────────────────────────────────────────────────────
  if (!file.exists(alliance_file))
    stop("Alliance orthology file not found: ", alliance_file,
         "\nDownload from: https://www.alliancegenome.org/downloads")

  ortholog_table <- convert_zf_to_human_alliance(unique(all_markers$gene), alliance_file)

  zf_to_human <- function(genes) {
    hits <- ortholog_table[ortholog_table$zebrafish_gene_id %in% genes, "human_gene_symbol"]
    unique(hits[!is.na(hits) & hits != ""])
  }

  # ── Enrichr setup ───────────────────────────────────────────────────────────
  setEnrichrSite("Enrichr")
  target_dbs <- c(
    "CellMarker_2024", "Descartes_Cell_Types_and_Tissue_2021",
    "CCLE_Proteomics_2020", "Tabula_Sapiens",
    "HuBMAP_ASCTplusB_augmented_2022", "CellMarker_Augmented_2021",
    "PanglaoDB_Augmented_2021", "Tabula_Muris", "Human_Gene_Atlas",
    "Azimuth_Cell_Types_2021", "ProteomicsDB_2020", "Mouse_Gene_Atlas",
    "ARCHS4_Tissues", "GTEx_Tissues_V8_2023", "Azimuth_2023"
  )
  available_dbs <- tryCatch(listEnrichrDbs()$libraryName, error = function(e) character(0))
  dbs_use       <- intersect(target_dbs, available_dbs)
  cat("Querying", length(dbs_use), "Enrichr databases\n")

  # ── XLSX helper (single sheet, color-coded by DB) ───────────────────────────
  db_palette <- c("FFE4E1","E1F0FF","E1FFE8","FFF8E1","F0E1FF",
                  "E1FFFF","FFE8D6","E8FFE1","FFE1F0","E1E8FF",
                  "FFFDE1","F5E1FF","E1FFF5","FFE1E1","E1FFEB")

  .write_xlsx <- function(per_db_results, filepath) {
    combined    <- do.call(rbind, per_db_results)
    dbs_present <- unique(combined$database)
    wb          <- createWorkbook(); addWorksheet(wb, "Enrichr_results")
    hdr_style   <- createStyle(fontColour="#000000", fgFill="#D9D9D9", halign="left",
                               textDecoration="Bold", border="Bottom", borderColour="#888888")
    writeData(wb, sheet=1, combined, headerStyle=hdr_style)
    color_map <- setNames(db_palette[seq_along(dbs_present)], dbs_present)
    for (db in dbs_present) {
      rows  <- which(combined$database == db) + 1L
      style <- createStyle(fgFill=paste0("#", color_map[[db]]))
      addStyle(wb, sheet=1, style=style, rows=rows, cols=seq_len(ncol(combined)),
               gridExpand=TRUE, stack=TRUE)
    }
    freezePane(wb, sheet=1, firstRow=TRUE)
    setColWidths(wb, sheet=1, cols=seq_len(ncol(combined)), widths="auto")
    saveWorkbook(wb, filepath, overwrite=TRUE)
  }

  # ── Fuzzy consensus ─────────────────────────────────────────────────────────
  .fuzzy_consensus <- function(top_terms, sim_thr=sim_thresh, min_l=min_libs) {
    all_terms <- unlist(top_terms)
    clean     <- trimws(gsub("\\s+(Human|Mouse|FACS|10x|Smart-seq2)$", "", all_terms, ignore.case=TRUE))
    if (length(clean) <= 1) return(NULL)
    sim_mat  <- as.matrix(1 - stringdistmatrix(tolower(clean), method="jw", p=0.1))
    groups   <- integer(length(clean)); gid <- 0L
    for (i in seq_along(clean)) {
      if (groups[i] != 0L) next
      gid <- gid + 1L; groups[i] <- gid
      for (j in seq_along(clean))
        if (j != i && groups[j] == 0L && sim_mat[i,j] >= sim_thr) groups[j] <- gid
    }
    db_per_term <- rep(names(top_terms), lengths(top_terms))
    group_libs  <- tapply(db_per_term, groups, function(x) length(unique(x)))
    best_gid    <- as.integer(names(which.max(group_libs)))
    if (group_libs[as.character(best_gid)] < min_l) return(NULL)
    members <- clean[groups == best_gid]; orig <- all_terms[groups == best_gid]
    list(label=members[which.min(nchar(members))],
         n_libs=as.integer(group_libs[as.character(best_gid)]),
         supporting=unique(orig))
  }

  # ── Main cluster loop ────────────────────────────────────────────────────────
  clusters     <- sort(unique(all_markers$cluster))
  summary_rows <- list()

  for (cl in clusters) {
    cat("\n-- Cluster", cl, "--\n")
    blank_row <- data.frame(cluster=cl, n_markers_sent=0, cell_type=NA_character_,
                             n_libs_support=0L, supporting_terms=NA_character_,
                             stringsAsFactors=FALSE)
    cl_df <- all_markers[all_markers$cluster == cl & !is.na(all_markers$p_val_adj) &
                         all_markers$p_val_adj < p_adj_cut & all_markers$avg_log2FC > 0, ]

    if (nrow(cl_df) < min_genes) {
      cat("  Too few markers — skipping\n"); summary_rows[[as.character(cl)]] <- blank_row; next
    }
    cl_df       <- cl_df[order(-cl_df$avg_log2FC), ]
    zf_genes    <- head(cl_df$gene, top_n)
    human_genes <- zf_to_human(zf_genes)

    if (length(human_genes) < min_genes) {
      cat("  Too few orthologs — skipping\n")
      blank_row$n_markers_sent <- length(human_genes)
      summary_rows[[as.character(cl)]] <- blank_row; next
    }
    cat("  Sending", length(human_genes), "human genes to Enrichr\n")

    enriched <- tryCatch(enrichr(human_genes, dbs_use),
                         error=function(e) { cat("  Enrichr error:", e$message, "\n"); NULL })
    if (is.null(enriched)) {
      blank_row$n_markers_sent <- length(human_genes)
      summary_rows[[as.character(cl)]] <- blank_row; next
    }

    per_db <- lapply(dbs_use, function(db) {
      res <- enriched[[db]]
      if (is.null(res) || nrow(res) == 0) return(NULL)
      res$database <- db; res$zf_genes <- paste(zf_genes, collapse=";")
      res$human_genes <- paste(human_genes, collapse=";")
      head(res[order(res$Adjusted.P.value), ], n_top_terms)
    })
    names(per_db) <- dbs_use; per_db <- Filter(Negate(is.null), per_db)

    if (length(per_db) > 0)
      .write_xlsx(per_db, file.path(out_dir, paste0("cluster_", cl, "_enrichr.xlsx")))

    consensus <- tryCatch(.fuzzy_consensus(lapply(per_db, function(d) d$Term)),
                          error=function(e) NULL)

    if (!is.null(consensus)) {
      cat("  ->", consensus$label, "(", consensus$n_libs, "libraries)\n")
      summary_rows[[as.character(cl)]] <- data.frame(
        cluster=cl, n_markers_sent=length(human_genes),
        cell_type=consensus$label, n_libs_support=consensus$n_libs,
        supporting_terms=paste(consensus$supporting, collapse=" | "),
        stringsAsFactors=FALSE)
    } else {
      cat("  -> No confident consensus (< 2 libraries) — NA assigned\n")
      blank_row$n_markers_sent <- length(human_genes)
      summary_rows[[as.character(cl)]] <- blank_row
    }
  }

  # ── Summary XLSX ─────────────────────────────────────────────────────────────
  enrichr_summary <- do.call(rbind, summary_rows)
  cat("\n=== Enrichr predictions ===\n")
  print(enrichr_summary[, c("cluster","cell_type","n_libs_support")])

  wb_s <- createWorkbook(); addWorksheet(wb_s, "predictions")
  hdr  <- createStyle(fontColour="#000000", fgFill="#D9D9D9", halign="left",
                      textDecoration="Bold", border="Bottom", borderColour="#888888")
  writeData(wb_s, sheet=1, enrichr_summary, headerStyle=hdr)
  freezePane(wb_s, sheet=1, firstRow=TRUE)
  setColWidths(wb_s, sheet=1, cols=seq_len(ncol(enrichr_summary)), widths="auto")
  saveWorkbook(wb_s, file.path(summary_dir, "enrichr_predictions.xlsx"), overwrite=TRUE)
  cat("Summary saved:", file.path(summary_dir, "enrichr_predictions.xlsx"), "\n")
  return(invisible(enrichr_summary))
}


# =============================================================================
# 3. DATA CHARACTERISTICS & NGUYEN 2023 DE METHOD SELECTION
# =============================================================================

# -----------------------------------------------------------------------------
# compute_data_characteristics()
#
# Compute sparsity (zero rate) and median UMI depth per cell type.
# Call before run_de_per_celltype() to inform method selection.
#
# Returns: data frame — cell_type, sparsity, median_depth, n_cells
# -----------------------------------------------------------------------------
compute_data_characteristics <- function(seurat_obj, cell_type_col = "cell_type") {
  DefaultAssay(seurat_obj) <- "RNA"
  ct_list <- levels(factor(seurat_obj@meta.data[[cell_type_col]]))
  results <- lapply(ct_list, function(ct) {
    cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == ct]
    if (length(cells) < 10) return(NULL)
    mat  <- GetAssayData(seurat_obj[, cells], assay="RNA", slot="counts")
    data.frame(cell_type    = ct,
               sparsity     = round(mean(mat == 0), 3),
               median_depth = round(median(Matrix::colSums(mat)), 1),
               n_cells      = length(cells))
  })
  out <- do.call(rbind, Filter(Negate(is.null), results))
  cat("Data characteristics:\n"); print(out); invisible(out)
}


# -----------------------------------------------------------------------------
# recommend_de_method()
#
# Nguyen 2023 (doi:10.1093/bib/bbad218) lookup table.
# Categorizes sparsity (>= 0.6 = high), depth (>= 40 = high, 8-39 = medium,
# < 8 = low) and batch status to recommend DE workflows.
#
# Returns: list(primary, recommended, profile)
# -----------------------------------------------------------------------------
recommend_de_method <- function(sparsity, depth, has_batch = FALSE) {

  nguyen_table <- list(
    # ── High sparsity (≈80%) ──────────────────────────────────────────────────
    high_high_substantial = list(
      label   = "High sparsity / High depth / Substantial batch",
      methods = c("MAST_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","scVI_limmatrend","DESeq2_FEM","limmatrend_Cov"),
      primary = "limmatrend_Cov"),
    high_high_small = list(
      label   = "High sparsity / High depth / Small batch",
      methods = c("MAST","ZW_edgeR","ZW_DESeq2","Pseudobulk_limma","DESeq2_FEM","limmatrend_Cov"),
      primary = "Pseudobulk_limma"),
    high_medium_substantial = list(
      label   = "High sparsity / Medium depth / Substantial batch",
      methods = c("DESeq2_Cov","limmatrend_Cov","DESeq2_wFisher","LogN_FEM","MAST_Cov","Raw_Wilcox"),
      primary = "limmatrend_Cov"),
    high_medium_small = list(
      label   = "High sparsity / Medium depth / Small batch",
      methods = c("DESeq2","limmatrend_Cov","LogN_FEM","Pseudobulk_edgeR","Pseudobulk_limma","Raw_Wilcox"),
      primary = "Pseudobulk_limma"),
    high_low_substantial = list(
      label   = "High sparsity / Low depth / Substantial batch",
      methods = c("LogN_FEM","limmatrend","Raw_Wilcox","RISC_QP"),
      primary = "limmatrend"),
    high_low_small = list(
      label   = "High sparsity / Low depth / Small batch",
      methods = c("LogN_FEM","limmatrend","Raw_Wilcox","RISC_QP"),
      primary = "limmatrend"),

    # ── Low sparsity (≈40%) ───────────────────────────────────────────────────
    low_high_substantial = list(
      label   = "Low sparsity / High depth / Substantial batch",
      methods = c("MAST_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","limma_BEC_Wilcox","Scanorama_limmatrend","LogN_FEM"),
      primary = "limmatrend_Cov"),
    low_high_small = list(
      label   = "Low sparsity / High depth / Small batch",
      methods = c("Pseudobulk_limma","Raw_Wilcox","Pseudobulk_edgeR","DESeq2","LogN_FEM"),
      primary = "Pseudobulk_limma"),
    low_low_substantial = list(
      label   = "Low sparsity / Low depth / Substantial batch",
      methods = c("limmavoom_Cov","limmatrend_Cov","ZW_edgeR_Cov","ZW_DESeq2_Cov","LogN_FEM","limma_BEC_Wilcox"),
      primary = "limmatrend_Cov"),
    low_low_small = list(
      label   = "Low sparsity / Low depth / Small batch",
      methods = c("Pseudobulk_limma","limmavoom","limmatrend","Pseudobulk_edgeR","Raw_Wilcox","edgeR"),
      primary = "Pseudobulk_limma")
  )

  sp_cat <- if (sparsity >= 0.6) "high" else "low"
  dp_cat <- if (depth >= 40) "high" else if (depth >= 8) "medium" else "low"
  if (sp_cat == "low" && dp_cat == "medium") dp_cat <- "low"   # no medium row for low sparsity
  bt_cat <- if (has_batch) "substantial" else "small"
  key    <- paste(sp_cat, dp_cat, bt_cat, sep="_")

  if (!key %in% names(nguyen_table)) {
    cat("⚠ No exact match for key:", key, "— falling back to high_medium_small\n")
    key <- "high_medium_small"
  }
  row <- nguyen_table[[key]]

  cat("─────────────────────────────────────────────────────\n")
  cat("DE method recommendation (Nguyen 2023)\n")
  cat("  Profile  :", row$label, "\n")
  cat("  Sparsity :", round(sparsity*100), "% →", sp_cat, "\n")
  cat("  Depth    :", depth, "UMIs →", dp_cat, "\n")
  cat("  Batch    :", bt_cat, "\n")
  cat("  PRIMARY  :", row$primary, "\n")
  cat("  All      :", paste(row$methods, collapse=", "), "\n")
  cat("─────────────────────────────────────────────────────\n")

  list(primary=row$primary, recommended=row$methods,
       profile=list(sparsity_cat=sp_cat, depth_cat=dp_cat, batch_cat=bt_cat))
}


# =============================================================================
# 4. NGUYEN 2023 DE METHOD IMPLEMENTATIONS
# =============================================================================
# All single-cell methods expect a Seurat object with RNA assay (raw counts).
# Standard output columns: gene, log2FoldChange, pvalue, padj, method

# -----------------------------------------------------------------------------
# run_pseudobulk_limma()   ← Pseudobulk_limma
# Aggregates per sample → limma-trend. Best all-around primary method.
# -----------------------------------------------------------------------------
run_pseudobulk_limma <- function(seurat_obj, cell_type,
                                  group_col="group", sample_col="sample",
                                  cell_type_col="cell_type", contrast_str=NULL,
                                  min_cells=5) {
  library(limma); library(edgeR)

  cells  <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < min_cells*2) { cat("  ⚠ Too few cells for", cell_type, "\n"); return(NULL) }
  obj    <- seurat_obj[, cells]; DefaultAssay(obj) <- "RNA"
  meta   <- obj@meta.data; counts <- GetAssayData(obj, assay="RNA", slot="counts")

  pb_list <- lapply(unique(meta[[sample_col]]), function(s) {
    sc <- rownames(meta)[meta[[sample_col]] == s]
    if (length(sc) < min_cells) return(NULL)
    Matrix::rowSums(counts[, sc, drop=FALSE])
  })
  names(pb_list) <- unique(meta[[sample_col]]); pb_list <- pb_list[!sapply(pb_list, is.null)]
  if (length(pb_list) < 3) { cat("  ⚠ Too few samples for", cell_type, "\n"); return(NULL) }

  pb_mat  <- do.call(cbind, pb_list)
  pb_meta <- meta[!duplicated(meta[[sample_col]]), ]
  rownames(pb_meta) <- pb_meta[[sample_col]]; pb_meta <- pb_meta[colnames(pb_mat),, drop=FALSE]

  y      <- edgeR::DGEList(counts=pb_mat); keep <- edgeR::filterByExpr(y, group=pb_meta[[group_col]])
  y      <- y[keep, ]
  design <- model.matrix(~ 0 + factor(pb_meta[[group_col]]))
  colnames(design) <- levels(factor(pb_meta[[group_col]]))

  if (is.null(contrast_str)) {
    gs <- levels(factor(pb_meta[[group_col]]))
    if (length(gs) != 2) stop("contrast_str required when >2 groups")
    contrast_str <- paste(rev(gs), collapse=" - ")
  }
  cont_mat <- limma::makeContrasts(contrasts=contrast_str, levels=design)
  logcpm   <- edgeR::cpm(y, log=TRUE, prior.count=3)
  fit      <- limma::lmFit(logcpm, design) %>% limma::contrasts.fit(cont_mat) %>%
              limma::eBayes(trend=TRUE, robust=TRUE)
  res <- limma::topTable(fit, number=Inf, sort.by="P") %>%
    dplyr::rename(log2FoldChange=logFC, pvalue=P.Value, padj=adj.P.Val) %>%
    tibble::rownames_to_column("gene") %>% dplyr::mutate(method="Pseudobulk_limma")

  cat("  ✓ Pseudobulk_limma:", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_limmatrend_cov()   ← limmatrend / limmatrend_Cov
# Single-cell limma-trend + optional batch covariate (TMM-normalized logCPM).
# -----------------------------------------------------------------------------
run_limmatrend_cov <- function(seurat_obj, cell_type,
                                group_col="group", sample_col="sample",
                                cell_type_col="cell_type", has_batch=TRUE) {
  library(limma); library(edgeR)

  cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj   <- seurat_obj[, cells]; DefaultAssay(obj) <- "RNA"
  meta  <- obj@meta.data
  cellinfo <- data.frame(Group=factor(meta[[group_col]]), Batch=factor(meta[[sample_col]]))

  y      <- edgeR::DGEList(GetAssayData(obj, assay="RNA", slot="counts"))
  y      <- edgeR::calcNormFactors(y)
  logCPM <- edgeR::cpm(y, log=TRUE, prior.count=3)
  design <- if (has_batch) model.matrix(~ Group + Batch, data=cellinfo) else
                           model.matrix(~ Group, data=cellinfo)
  fit    <- limma::lmFit(logCPM, design) %>% limma::eBayes(trend=TRUE, robust=TRUE)
  res    <- limma::topTable(fit, n=Inf, adjust.method="BH", coef=2) %>%
    dplyr::rename(log2FoldChange=logFC, pvalue=P.Value, padj=adj.P.Val) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=ifelse(has_batch, "limmatrend_Cov", "limmatrend"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_limmavoom_cov()   ← limmavoom / limmavoom_Cov
# limma-voom with TMM normalization + optional batch covariate.
# -----------------------------------------------------------------------------
run_limmavoom_cov <- function(seurat_obj, cell_type,
                               group_col="group", sample_col="sample",
                               cell_type_col="cell_type", has_batch=TRUE) {
  library(limma); library(edgeR)

  cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj   <- seurat_obj[, cells]; DefaultAssay(obj) <- "RNA"
  meta  <- obj@meta.data
  cellinfo  <- data.frame(Group=factor(meta[[group_col]]), Batch=factor(meta[[sample_col]]))
  count_mat <- as.matrix(GetAssayData(obj, assay="RNA", slot="counts"))
  nf        <- edgeR::calcNormFactors(count_mat)
  design    <- if (has_batch) model.matrix(~ Group + Batch, data=cellinfo) else
                              model.matrix(~ Group, data=cellinfo)
  voom_d    <- limma::voom(count_mat, design=design, lib.size=colSums(count_mat)*nf)
  fit       <- limma::lmFit(voom_d, design) %>% limma::eBayes()
  res       <- limma::topTable(fit, n=Inf, adjust.method="BH", coef=2) %>%
    dplyr::rename(log2FoldChange=logFC, pvalue=P.Value, padj=adj.P.Val) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=ifelse(has_batch, "limmavoom_Cov", "limmavoom"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_edgeR_cov()   ← Pseudobulk_edgeR / ZW_edgeR_Cov
# glmQLFit edgeR. pseudobulk=TRUE → Pseudobulk_edgeR; FALSE → ZW_edgeR_Cov.
# -----------------------------------------------------------------------------
run_edgeR_cov <- function(seurat_obj, cell_type,
                           group_col="group", sample_col="sample",
                           cell_type_col="cell_type", has_batch=TRUE, pseudobulk=FALSE) {
  library(edgeR)

  cells  <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj    <- seurat_obj[, cells]; DefaultAssay(obj) <- "RNA"
  meta   <- obj@meta.data; counts <- GetAssayData(obj, assay="RNA", slot="counts")

  if (pseudobulk) {
    pb_list <- lapply(unique(meta[[sample_col]]), function(s) {
      sc <- rownames(meta)[meta[[sample_col]] == s]
      if (length(sc) < 3) return(NULL); Matrix::rowSums(counts[, sc, drop=FALSE])
    })
    names(pb_list) <- unique(meta[[sample_col]]); pb_list <- pb_list[!sapply(pb_list, is.null)]
    if (length(pb_list) < 3) { cat("  ⚠ Too few samples\n"); return(NULL) }
    count_mat <- do.call(cbind, pb_list)
    pb_meta   <- meta[!duplicated(meta[[sample_col]]), ]; rownames(pb_meta) <- pb_meta[[sample_col]]
    pb_meta   <- pb_meta[colnames(count_mat),, drop=FALSE]
    cellGroup  <- factor(pb_meta[[group_col]]); cellBatch <- factor(pb_meta[[sample_col]])
  } else {
    count_mat <- as.matrix(counts)
    cellGroup  <- factor(meta[[group_col]]); cellBatch <- factor(meta[[sample_col]])
  }

  y      <- edgeR::DGEList(counts=count_mat, group=cellGroup) %>% edgeR::calcNormFactors()
  design <- if (has_batch) model.matrix(~ cellGroup + cellBatch) else model.matrix(~ cellGroup)
  rownames(design) <- colnames(y)
  y      <- edgeR::estimateDisp(y, design, robust=TRUE)
  qlf    <- edgeR::glmQLFit(y, design, robust=TRUE, prior.df=0) %>% edgeR::glmQLFTest(coef=2)
  res    <- qlf$table %>%
    dplyr::mutate(padj=p.adjust(PValue, method="BH")) %>%
    dplyr::rename(log2FoldChange=logFC, pvalue=PValue) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=if(pseudobulk) "Pseudobulk_edgeR" else
                  ifelse(has_batch, "ZW_edgeR_Cov", "edgeR"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_mast_cov()   ← MAST / MAST_Cov
# MAST via Seurat FindMarkers. has_batch=TRUE adds sample as latent variable.
# -----------------------------------------------------------------------------
run_mast_cov <- function(seurat_obj, cell_type,
                          group_col="group", sample_col="sample",
                          cell_type_col="cell_type", has_batch=TRUE) {
  if (!requireNamespace("MAST", quietly=TRUE)) stop("MAST package required")
  library(MAST)

  cells  <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj    <- seurat_obj[, cells]
  Idents(obj) <- obj@meta.data[[group_col]]
  groups <- sort(levels(factor(obj@meta.data[[group_col]])))

  markers <- tryCatch({
    if (has_batch)
      FindMarkers(obj, ident.1=groups[2], ident.2=groups[1], test.use="MAST",
                  logfc.threshold=0, latent.vars=sample_col,
                  min.cells.feature=0, min.cells.group=0, min.pct=0)
    else
      FindMarkers(obj, ident.1=groups[2], ident.2=groups[1], test.use="MAST",
                  logfc.threshold=0, min.cells.feature=0, min.cells.group=0, min.pct=0)
  }, error=function(e) { cat("  ✗ MAST error:", e$message, "\n"); NULL })

  if (is.null(markers)) return(NULL)
  res <- markers %>%
    dplyr::rename(log2FoldChange=avg_log2FC, pvalue=p_val, padj=p_val_adj) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(method=ifelse(has_batch, "MAST_Cov", "MAST"))

  cat("  ✓", unique(res$method), ":", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_wilcox_raw()   ← Raw_Wilcox
# Wilcoxon on log-normalized counts. Fallback for very low depth data only.
# WARNING: does not correct for batch effects.
# -----------------------------------------------------------------------------
run_wilcox_raw <- function(seurat_obj, cell_type,
                            group_col="group", cell_type_col="cell_type") {
  cells <- colnames(seurat_obj)[seurat_obj@meta.data[[cell_type_col]] == cell_type]
  if (length(cells) < 10) { cat("  ⚠ Too few cells\n"); return(NULL) }
  obj   <- seurat_obj[, cells]; DefaultAssay(obj) <- "RNA"
  obj   <- NormalizeData(obj, verbose=FALSE)
  Idents(obj) <- obj@meta.data[[group_col]]
  groups <- sort(levels(factor(obj@meta.data[[group_col]])))

  markers <- tryCatch(
    FindMarkers(obj, ident.1=groups[2], ident.2=groups[1], test.use="wilcox",
                logfc.threshold=0, min.cells.feature=0, min.cells.group=0, min.pct=0),
    error=function(e) { cat("  ✗ Wilcox error:", e$message, "\n"); NULL })

  if (is.null(markers)) return(NULL)
  res <- markers %>%
    dplyr::rename(log2FoldChange=avg_log2FC, pvalue=p_val, padj=p_val_adj) %>%
    tibble::rownames_to_column("gene") %>% dplyr::mutate(method="Raw_Wilcox")

  cat("  ✓ Raw_Wilcox:", nrow(res), "genes,", sum(res$padj < 0.05, na.rm=TRUE), "sig\n")
  return(res)
}


# -----------------------------------------------------------------------------
# run_de_per_celltype()
#
# Master dispatcher:
#   1. Runs compute_data_characteristics() per cell type
#   2. Calls recommend_de_method() for each
#   3. Dispatches to the recommended primary method
#   4. Falls back to Pseudobulk_limma if method not implemented
#
# Args:
#   seurat_obj    : Seurat object
#   cell_types    : character vector of cell types (NULL = all in cell_type_col)
#   group_col     : metadata group column
#   sample_col    : metadata sample/batch column
#   cell_type_col : metadata cell type column
#   contrast_str  : contrast string for pseudobulk methods
#   has_batch     : logical override (NULL = auto-detect from sample/group counts)
#
# Returns: named list [[cell_type]] of DE result data frames
# -----------------------------------------------------------------------------
run_de_per_celltype <- function(seurat_obj,
                                 cell_types    = NULL,
                                 group_col     = "group",
                                 sample_col    = "sample",
                                 cell_type_col = "cell_type",
                                 contrast_str  = NULL,
                                 has_batch     = NULL) {
  if (is.null(cell_types))
    cell_types <- levels(factor(seurat_obj@meta.data[[cell_type_col]]))

  if (is.null(has_batch)) {
    n_samples  <- length(unique(seurat_obj@meta.data[[sample_col]]))
    n_groups   <- length(unique(seurat_obj@meta.data[[group_col]]))
    has_batch  <- (n_samples > n_groups)
    cat("Auto-detected batch:", has_batch, "(", n_samples, "samples,", n_groups, "groups)\n")
  }

  char_df <- compute_data_characteristics(seurat_obj, cell_type_col)
  results <- list()

  for (ct in cell_types) {
    cat("\n=== DE for:", ct, "===\n")
    ct_char <- char_df[char_df$cell_type == ct, ]
    rec <- if (nrow(ct_char) == 0) {
      cat("  ⚠ No characteristics — defaulting to Pseudobulk_limma\n")
      list(primary="Pseudobulk_limma")
    } else {
      recommend_de_method(ct_char$sparsity, ct_char$median_depth, has_batch)
    }

    results[[ct]] <- switch(rec$primary,
      "Pseudobulk_limma"  = run_pseudobulk_limma(seurat_obj, ct, group_col, sample_col, cell_type_col, contrast_str),
      "limmatrend_Cov"    = run_limmatrend_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=TRUE),
      "limmatrend"        = run_limmatrend_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=FALSE),
      "limmavoom_Cov"     = run_limmavoom_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=TRUE),
      "limmavoom"         = run_limmavoom_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=FALSE),
      "Pseudobulk_edgeR"  = run_edgeR_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=FALSE, pseudobulk=TRUE),
      "ZW_edgeR_Cov"      = run_edgeR_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=TRUE,  pseudobulk=FALSE),
      "MAST_Cov"          = run_mast_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=TRUE),
      "MAST"              = run_mast_cov(seurat_obj, ct, group_col, sample_col, cell_type_col, has_batch=FALSE),
      "Raw_Wilcox"        = run_wilcox_raw(seurat_obj, ct, group_col, cell_type_col),
      {
        cat("  ⚠ Method", rec$primary, "not implemented — falling back to Pseudobulk_limma\n")
        run_pseudobulk_limma(seurat_obj, ct, group_col, sample_col, cell_type_col, contrast_str)
      }
    )
  }
  cat("\n✓✓✓ DE complete for all cell types ✓✓✓\n")
  return(results)
}

# Legacy alias
pseudobulk_de <- function(...) run_pseudobulk_limma(...)


# =============================================================================
# 5. LIANA — LIGAND-RECEPTOR ANALYSIS
# =============================================================================

run_liana_zebrafish <- function(target_organism=7955, resource_name="Consensus") {
  library(liana)
  cat("Generating zebrafish LIANA ortholog resource...\n")
  or <- liana::generate_homologs(op_resource=liana::select_resource(resource_name)[[1]],
                                  target_organism=target_organism)
  cat("✓ Resource ready:", nrow(or), "interactions\n"); return(or)
}

run_liana_loop <- function(seurat_obj, conditions, ortholog_resource,
                            source_groups, target_groups,
                            condition_col="condition", cell_type_col="cell_type",
                            agg_rank_cutoff=0.05, ntop=20, out_dir="LR_analysis") {
  library(liana)
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  if (is.character(source_groups)) source_groups <- list(all=source_groups)
  if (is.character(target_groups)) target_groups <- list(all=target_groups)
  if (length(source_groups) != length(target_groups))
    stop("source_groups and target_groups must have same length")

  results_all <- list()
  for (cond in conditions) {
    cat("\n=== LIANA: condition =", cond, "===\n")
    cells_cond <- colnames(seurat_obj)[seurat_obj@meta.data[[condition_col]] == cond]
    if (length(cells_cond) < 10) { cat("  ⚠ Too few cells\n"); next }
    seurat_sub <- seurat_obj[, cells_cond]
    Idents(seurat_sub) <- seurat_sub@meta.data[[cell_type_col]]
    DefaultAssay(seurat_sub) <- "SCT"

    liana_res <- tryCatch(
      liana::liana_wrap(seurat_sub, resource="custom", external_resource=ortholog_resource) %>%
        liana::liana_aggregate(),
      error=function(e) { cat("  ✗ LIANA error:", e$message, "\n"); NULL })
    if (is.null(liana_res) || nrow(liana_res) == 0) next

    results_all[[cond]] <- liana_res
    liana_sig <- liana_res %>% dplyr::filter(aggregate_rank < agg_rank_cutoff)
    write.csv(liana_res, file.path(out_dir, paste0("liana_", cond, ".csv")), row.names=FALSE)
    write.csv(liana_sig, file.path(out_dir, paste0("liana_", cond, "_sig_agg", agg_rank_cutoff, ".csv")), row.names=FALSE)
    cat("  ✓", nrow(liana_res), "total,", nrow(liana_sig), "significant\n")

    for (pair_name in names(source_groups))
      plot_liana_dotplot(liana_sig, source_groups[[pair_name]], target_groups[[pair_name]],
                         cond, pair_name, ntop, out_dir)
  }
  cat("\n✓✓✓ LIANA complete ✓✓✓\n"); return(invisible(results_all))
}

plot_liana_dotplot <- function(liana_result, source_groups, target_groups,
                                condition="cond", pair_name="pair", ntop=20,
                                out_dir="LR_analysis", width=12, height=14) {
  library(liana)
  src_str <- paste(source_groups, collapse="_"); tgt_str <- paste(target_groups, collapse="_")
  fname   <- file.path(out_dir, paste0("liana_", condition, "_", pair_name,
                                        "_", src_str, "__TO__", tgt_str, "_agg05.pdf"))
  plot_data <- liana_result %>% dplyr::filter(source %in% source_groups | target %in% target_groups)
  if (nrow(plot_data) == 0) { cat("  ⚠ No interactions to plot\n"); return(invisible(NULL)) }
  tryCatch({
    p <- liana::liana_dotplot(plot_data, source_groups=source_groups, target_groups=target_groups,
           specificity="natmi.edge_specificity", magnitude="aggregate_rank",
           colour.label="Aggregate Rank\n(p-val)", ntop=ntop) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    ggsave(fname, p, width=width, height=height, limitsize=FALSE)
    cat("  ✓ Dotplot →", fname, "\n")
  }, error=function(e) cat("  ✗ Dotplot error:", e$message, "\n"))
}


# =============================================================================
# 6. SC ENRICHMENT WRAPPERS (adds <CellType> subfolder)
# =============================================================================

run_ORA_loop_sc <- function(de_results_sc, contrast, thresholds, org="mouse", num_pathways=15) {
  ontologies <- c("BP","MF","CC","reactome"); gene_sets <- c("all","up","down")
  cat("Starting single-cell ORA for contrast:", contrast, "\n")
  for (cell_type in names(de_results_sc)) {
    cat("\n  ---", cell_type, "---\n")
    res_df <- de_results_sc[[cell_type]]
    if (is.null(res_df) || nrow(res_df) == 0) { cat("  ⚠ No results\n"); next }
    ct_safe <- gsub("[^A-Za-z0-9_.-]","_", cell_type)
    for (thresh_name in names(thresholds)) {
      thresh  <- thresholds[[thresh_name]]
      gs_list <- list(
        all  = res_df %>% dplyr::filter(padj < 0.05) %>% dplyr::pull(gene),
        up   = res_df %>% dplyr::filter(padj < 0.05, log2FoldChange >  thresh) %>% dplyr::pull(gene),
        down = res_df %>% dplyr::filter(padj < 0.05, log2FoldChange < -thresh) %>% dplyr::pull(gene)
      )
      for (onto in ontologies)
        for (set_name in gene_sets) {
          genes <- gs_list[[set_name]]
          if (length(genes) < 3) next
          ORA_function(gene_list=genes, pval=0.05, onto=onto,
                       contrast=file.path(contrast, ct_safe), lfc_label=thresh_name,
                       gene_set=set_name, org=org, num_pathways=num_pathways)
        }
    }
  }
  cat("\n✓✓✓ Single-cell ORA complete ✓✓✓\n")
}

run_GSEA_loop_sc <- function(de_results_sc, contrast, thresholds,
                              organism="mouse", num_pathways=15, top_n_ranked=3) {
  gsea_analyses <- list(list(onto="BP",msigdb=NULL), list(onto="MF",msigdb=NULL),
                        list(onto="CC",msigdb=NULL), list(onto="reactome",msigdb=NULL),
                        list(onto="msigdb",msigdb="H"))
  cat("Starting single-cell GSEA for contrast:", contrast, "\n")
  for (cell_type in names(de_results_sc)) {
    cat("\n  ---", cell_type, "---\n")
    res_df <- de_results_sc[[cell_type]]
    if (is.null(res_df) || nrow(res_df) == 0) { cat("  ⚠ No results\n"); next }
    ct_safe <- gsub("[^A-Za-z0-9_.-]","_", cell_type)
    for (thresh_name in names(thresholds)) {
      lfc_thresh_val <- thresholds[[thresh_name]]
      for (analysis in gsea_analyses)
        GSEA_function(df=res_df, pval_deg=1.0, pval_enrich=0.2, onto=analysis$onto,
                      contrast=file.path(contrast, ct_safe), lfc_label=thresh_name,
                      lfc_threshold=lfc_thresh_val, organism=organism,
                      msigdb_collection=analysis$msigdb,
                      num_pathways=num_pathways, top_n_ranked=top_n_ranked)
    }
  }
  cat("\n✓✓✓ Single-cell GSEA complete ✓✓✓\n")
}


cat("✓ sc_functions.R loaded\n")
cat("  Clustering QC  : compute_clustering_qc(), plot_clustering_qc()\n")
cat("  Annotation     : annotate_celltypes_enrichr()\n")
cat("  DE selection   : compute_data_characteristics(), recommend_de_method()\n")
cat("  DE dispatcher  : run_de_per_celltype()\n")
cat("  DE methods     : run_pseudobulk_limma(), run_limmatrend_cov(),\n")
cat("                   run_limmavoom_cov(), run_edgeR_cov(),\n")
cat("                   run_mast_cov(), run_wilcox_raw()\n")
cat("  LIANA          : run_liana_zebrafish(), run_liana_loop(), plot_liana_dotplot()\n")
cat("  SC enrichment  : run_ORA_loop_sc(), run_GSEA_loop_sc()\n")
