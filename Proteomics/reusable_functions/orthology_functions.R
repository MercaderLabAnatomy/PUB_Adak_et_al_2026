# =============================================================================
# Mercader Lab — Orthology Conversion Functions
# orthology_functions.R
#
# Functions for converting gene identifiers across species using orthology
# databases (Alliance Genome, Ensembl biomaRt, etc.)
# =============================================================================

#' Convert Zebrafish Genes to Mouse Orthologs Using Alliance Database
#'
#' Converts zebrafish gene symbols to mouse orthologs using high-quality
#' Alliance Genome orthology data with bidirectional best score filtering.
#'
#' @param zf_genes Character vector of zebrafish gene symbols
#' @param alliance_file Path to Alliance orthology TSV file
#'        Default: "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv"
#'
#' @return Data frame with columns:
#'   - zebrafish_gene_id: Input zebrafish gene symbol
#'   - mouse_gene_symbol: Corresponding mouse gene symbol
#'   Returns NULL if file not found or parsing fails.
#'
#' @details
#' Uses Alliance Genome orthology data filtered for high confidence:
#'   - IsBestScore = "Yes" (best zebrafish→mouse match)
#'   - IsBestRevScore = "Yes" (best mouse→zebrafish match)
#'
#' Bidirectional mapping ensures orthologs are reciprocal best hits.
#'
#' Download Alliance file from:
#' https://www.alliancegenome.org/downloads
#'
#' @examples
#' zf_genes <- c("tp53", "myca", "nkx2.5")
#' orthologs <- convert_zf_to_mouse_alliance(zf_genes)
#'
#' @export
convert_zf_to_mouse_alliance <- function(zf_genes, 
                                         alliance_file = "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv") {
    
    # Check file exists
    if (!file.exists(alliance_file)) {
        cat("ERROR: Alliance orthology file not found:", alliance_file, "\n")
        cat("Please download from Alliance Genome (https://www.alliancegenome.org/downloads)\n")
        return(NULL)
    }
    
    cat("Reading Alliance orthology file...\n")
    
    tryCatch({
        # Read Alliance TSV file
        alliance_data <- read.delim(alliance_file, 
                                   header = TRUE, 
                                   sep = "\t", 
                                   stringsAsFactors = FALSE, 
                                   quote = "", 
                                   comment.char = "#")
        
        # Handle single-column parsing edge case
        if (ncol(alliance_data) == 1) {
            alliance_data <- read.table(alliance_file, 
                                       header = TRUE, 
                                       sep = "\t", 
                                       stringsAsFactors = FALSE, 
                                       quote = "", 
                                       fill = TRUE, 
                                       comment.char = "")
        }
        
        # Validate expected columns
        expected_cols <- c("Gene1Symbol", "Gene2Symbol", "Gene1SpeciesName", 
                          "Gene2SpeciesName", "IsBestScore", "IsBestRevScore")
        if (any(!expected_cols %in% colnames(alliance_data))) {
            cat("ERROR: Missing required columns in Alliance file\n")
            cat("Expected:", paste(expected_cols, collapse = ", "), "\n")
            cat("Found:", paste(colnames(alliance_data), collapse = ", "), "\n")
            return(NULL)
        }
        
        # Filter for high-quality orthologs (bidirectional best scores)
        high_quality <- alliance_data[
            alliance_data$IsBestScore == "Yes" & alliance_data$IsBestRevScore == "Yes", 
        ]
        
        cat("High-quality ortholog pairs (bidirectional best scores):", nrow(high_quality), "\n")
        
        # Zebrafish → mouse direction
        zf_mouse <- high_quality[
            high_quality$Gene1SpeciesName == "Danio rerio" &
            high_quality$Gene2SpeciesName == "Mus musculus", 
        ]
        
        # Mouse → zebrafish direction (reverse)
        mouse_zf <- high_quality[
            high_quality$Gene1SpeciesName == "Mus musculus" &
            high_quality$Gene2SpeciesName == "Danio rerio", 
        ]
        
        # Combine both directions for comprehensive mapping
        mapping <- rbind(
            data.frame(zebrafish_gene_id = zf_mouse$Gene1Symbol,
                      mouse_gene_symbol = zf_mouse$Gene2Symbol,
                      stringsAsFactors = FALSE),
            data.frame(zebrafish_gene_id = mouse_zf$Gene2Symbol,
                      mouse_gene_symbol = mouse_zf$Gene1Symbol,
                      stringsAsFactors = FALSE)
        )
        
        # Clean and deduplicate
        mapping <- mapping[!is.na(mapping$zebrafish_gene_id) & 
                          !is.na(mapping$mouse_gene_symbol) &
                          mapping$zebrafish_gene_id != "" & 
                          mapping$mouse_gene_symbol != "", ]
        mapping <- mapping[!duplicated(mapping), ]
        
        cat("Total unique ortholog pairs:", nrow(mapping), "\n")
        
        # Filter for input genes
        result <- mapping[mapping$zebrafish_gene_id %in% zf_genes, ]
        cat("Mapped", length(unique(result$zebrafish_gene_id)), "of", 
            length(zf_genes), "input genes to mouse orthologs\n")
        
        return(result)
        
    }, error = function(e) {
        cat("ERROR reading Alliance file:", e$message, "\n")
        return(NULL)
    })
}

#' Convert Zebrafish Genes to Human Orthologs Using Alliance Database
#'
#' Converts zebrafish gene symbols to human orthologs using high-quality
#' Alliance Genome orthology data with bidirectional best score filtering.
#' Mirror of convert_zf_to_mouse_alliance() for the Danio rerio <-> Homo sapiens
#' direction; used by annotate_celltypes_enrichr() for Enrichr cell-type calls.
#'
#' @param zf_genes Character vector of zebrafish gene symbols
#' @param alliance_file Path to Alliance orthology TSV file
#'        Default: "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv"
#'
#' @return Data frame with columns:
#'   - zebrafish_gene_id: Input zebrafish gene symbol
#'   - human_gene_symbol: Corresponding human gene symbol
#'   Returns NULL if file not found or parsing fails.
#'
#' @details
#' Uses Alliance Genome orthology data filtered for high confidence:
#'   - IsBestScore = "Yes" (best zebrafish->human match)
#'   - IsBestRevScore = "Yes" (best human->zebrafish match)
#'
#' Bidirectional mapping ensures orthologs are reciprocal best hits.
#'
#' Download Alliance file from:
#' https://www.alliancegenome.org/downloads
#'
#' @examples
#' zf_genes <- c("tp53", "myca", "nkx2.5")
#' orthologs <- convert_zf_to_human_alliance(zf_genes)
#'
#' @export
convert_zf_to_human_alliance <- function(zf_genes,
                                         alliance_file = "ORTHOLOGY-ALLIANCE_COMBINED_7_23052025.tsv") {

    # Check file exists
    if (!file.exists(alliance_file)) {
        cat("ERROR: Alliance orthology file not found:", alliance_file, "\n")
        cat("Please download from Alliance Genome (https://www.alliancegenome.org/downloads)\n")
        return(NULL)
    }

    cat("Reading Alliance orthology file...\n")

    tryCatch({
        # Read Alliance TSV file
        alliance_data <- read.delim(alliance_file,
                                   header = TRUE,
                                   sep = "\t",
                                   stringsAsFactors = FALSE,
                                   quote = "",
                                   comment.char = "#")

        # Handle single-column parsing edge case
        if (ncol(alliance_data) == 1) {
            alliance_data <- read.table(alliance_file,
                                       header = TRUE,
                                       sep = "\t",
                                       stringsAsFactors = FALSE,
                                       quote = "",
                                       fill = TRUE,
                                       comment.char = "")
        }

        # Validate expected columns
        expected_cols <- c("Gene1Symbol", "Gene2Symbol", "Gene1SpeciesName",
                          "Gene2SpeciesName", "IsBestScore", "IsBestRevScore")
        if (any(!expected_cols %in% colnames(alliance_data))) {
            cat("ERROR: Missing required columns in Alliance file\n")
            cat("Expected:", paste(expected_cols, collapse = ", "), "\n")
            cat("Found:", paste(colnames(alliance_data), collapse = ", "), "\n")
            return(NULL)
        }

        # Filter for high-quality orthologs (bidirectional best scores)
        high_quality <- alliance_data[
            alliance_data$IsBestScore == "Yes" & alliance_data$IsBestRevScore == "Yes",
        ]

        cat("High-quality ortholog pairs (bidirectional best scores):", nrow(high_quality), "\n")

        # Zebrafish → human direction
        zf_human <- high_quality[
            high_quality$Gene1SpeciesName == "Danio rerio" &
            high_quality$Gene2SpeciesName == "Homo sapiens",
        ]

        # Human → zebrafish direction (reverse)
        human_zf <- high_quality[
            high_quality$Gene1SpeciesName == "Homo sapiens" &
            high_quality$Gene2SpeciesName == "Danio rerio",
        ]

        # Combine both directions for comprehensive mapping
        mapping <- rbind(
            data.frame(zebrafish_gene_id = zf_human$Gene1Symbol,
                      human_gene_symbol = zf_human$Gene2Symbol,
                      stringsAsFactors = FALSE),
            data.frame(zebrafish_gene_id = human_zf$Gene2Symbol,
                      human_gene_symbol = human_zf$Gene1Symbol,
                      stringsAsFactors = FALSE)
        )

        # Clean and deduplicate
        mapping <- mapping[!is.na(mapping$zebrafish_gene_id) &
                          !is.na(mapping$human_gene_symbol) &
                          mapping$zebrafish_gene_id != "" &
                          mapping$human_gene_symbol != "", ]
        mapping <- mapping[!duplicated(mapping), ]

        cat("Total unique ortholog pairs:", nrow(mapping), "\n")

        # Filter for input genes
        result <- mapping[mapping$zebrafish_gene_id %in% zf_genes, ]
        cat("Mapped", length(unique(result$zebrafish_gene_id)), "of",
            length(zf_genes), "input genes to human orthologs\n")

        return(result)

    }, error = function(e) {
        cat("ERROR reading Alliance file:", e$message, "\n")
        return(NULL)
    })
}
