# ==============================================================================
# LLM Prompt Builder Functions for FragPipe-Analyst
#
# These functions extract analysis results from the app's reactive objects and
# assemble them into structured, context-rich prompts ready for LLM input.
# ==============================================================================

# ------------------------------------------------------------------------------
# Gene name resolution
# ------------------------------------------------------------------------------

#' Resolve the best available gene name column from rowData
#'
#' Handles the different column naming conventions across LFQ/TMT/DIA and
#' protein/peptide/site levels. Returns a character vector of gene labels
#' aligned to the rows of rd.
#'
#' @param rd data.frame from as.data.frame(rowData(dep))
#' @param exp  character: metadata(dep)$exp
#' @param level character: metadata(dep)$level
resolve_gene_labels <- function(rd, exp, level) {
  # TMT gene-level: gene is the row ID
  if (exp == "TMT" && level == "gene") {
    return(as.character(rd$ID))
  }
  # TMT protein/peptide/site: dedicated Gene column
  if (exp == "TMT" && "Gene" %in% colnames(rd)) {
    return(as.character(rd$Gene))
  }
  # DIA peptide/site with Gene column
  if ("Gene" %in% colnames(rd)) {
    return(as.character(rd$Gene))
  }
  # DIA protein: sometimes Genes (plural)
  if ("Genes" %in% colnames(rd)) {
    return(as.character(rd$Genes))
  }
  # LFQ / DIA protein: use 'name', strip isoform suffixes (e.g. "EGFR.1" -> "EGFR")
  if ("name" %in% colnames(rd)) {
    return(gsub("[.].*", "", as.character(rd$name)))
  }
  # Last resort
  as.character(rd$ID)
}


# ------------------------------------------------------------------------------
# DE gene extraction
# ------------------------------------------------------------------------------

#' Extract significant DE genes for a specific contrast from a SummarizedExperiment
#'
#' @param dep   SummarizedExperiment with rowData containing DE columns
#' @param contrast  Character like "Drug_A_vs_DMSO"
#' @param lfc_thresh  Minimum absolute log2 fold change (used for up/down split)
#' @param alpha  Adjusted p-value cutoff
#' @return Named list: $up and $down, each a data.frame with columns
#'         gene, lfc, padj, (description)
extract_de_genes <- function(dep, contrast, lfc_thresh = 1, alpha = 0.05) {
  rd <- as.data.frame(rowData(dep, use.names = FALSE))

  diff_col <- paste0(contrast, "_diff")
  padj_col <- paste0(contrast, "_p.adj")
  sig_col  <- paste0(contrast, "_significant")

  if (!all(c(diff_col, padj_col) %in% colnames(rd))) {
    return(list(up = data.frame(), down = data.frame(), n_up = 0, n_down = 0))
  }

  exp   <- metadata(dep)$exp
  level <- metadata(dep)$level

  rd$gene  <- resolve_gene_labels(rd, exp, level)
  rd$lfc   <- as.numeric(rd[[diff_col]])
  rd$padj  <- as.numeric(rd[[padj_col]])

  # Use pre-computed significance flag when available; fall back to threshold filtering
  if (sig_col %in% colnames(rd)) {
    rd$sig <- !is.na(rd[[sig_col]]) & rd[[sig_col]]
  } else {
    rd$sig <- !is.na(rd$padj) & !is.na(rd$lfc) &
              rd$padj < alpha & abs(rd$lfc) > lfc_thresh
  }

  # Pull optional protein description (helps LLM with obscure gene names)
  has_desc <- "Description" %in% colnames(rd)

  sig_rows <- rd[rd$sig & !is.na(rd$lfc) & !is.na(rd$padj), , drop = FALSE]

  up   <- sig_rows[sig_rows$lfc >  0, ]
  down <- sig_rows[sig_rows$lfc <  0, ]

  # Sort by magnitude
  up   <- up[order(-up$lfc), ]
  down <- down[order(down$lfc), ]

  # Keep only needed columns
  cols <- c("gene", "lfc", "padj")
  if (has_desc) cols <- c(cols, "Description")

  list(
    up     = up[, cols, drop = FALSE],
    down   = down[, cols, drop = FALSE],
    n_up   = nrow(up),
    n_down = nrow(down)
  )
}


# ------------------------------------------------------------------------------
# Formatting helpers
# ------------------------------------------------------------------------------

#' Format a ranked gene data.frame as human-readable text
#'
#' @param df        data.frame with columns gene, lfc, padj, (Description)
#' @param top_n     Maximum rows to include
#' @param show_desc Logical: include protein description snippet when available
format_gene_list <- function(df, top_n = 30, show_desc = TRUE) {
  if (is.null(df) || nrow(df) == 0) return("  (none)")
  df <- head(df, top_n)

  has_desc <- show_desc && "Description" %in% colnames(df) &&
              !all(is.na(df$Description)) && !all(df$Description == "")

  lines <- mapply(function(gene, lfc, padj, desc) {
    stat_part <- sprintf("%-14s  log2FC=%+.2f  adj.p=%s",
                         gene, lfc, formatC(padj, format = "e", digits = 1))
    if (has_desc && !is.na(desc) && nchar(trimws(desc)) > 0) {
      # Abbreviate long descriptions
      desc_short <- substr(trimws(desc), 1, 60)
      if (nchar(trimws(desc)) > 60) desc_short <- paste0(desc_short, "...")
      paste0("  ", stat_part, "\n    (", desc_short, ")")
    } else {
      paste0("  ", stat_part)
    }
  },
  df$gene, df$lfc, df$padj,
  if (has_desc) df$Description else rep("", nrow(df)),
  SIMPLIFY = TRUE)

  paste(lines, collapse = "\n")
}


#' Format enrichment results for a given contrast as human-readable text
#'
#' Handles both clusterProfiler output (columns: Adjusted.P.value, IN, geneID)
#' and Enrichr output (columns: Adjusted.P.value, IN, Genes).
#'
#' @param enrich_df  data.frame from go_results() / pathway_results()
#' @param contrast   character: contrast name to filter on
#' @param top_n      max rows to include
#' @param alpha      significance cutoff for adj.p (or P.value if adjust=FALSE)
#' @param adjust     logical: use adjusted p-values (Adjusted.P.value)
format_enrichment_table <- function(enrich_df, contrast,
                                    top_n = 20, alpha = 0.05, adjust = TRUE) {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    return("  (no enrichment results available)")
  }

  # Column name convention from test_ora_mod():
  #   clusterProfiler: Adjusted.P.value, P.value, IN (count), geneID ("/"-sep)
  #   enrichr:         Adjusted.P.value, P.value, IN (count), Genes  (";"-sep)
  p_col <- if (adjust && "Adjusted.P.value" %in% colnames(enrich_df)) {
    "Adjusted.P.value"
  } else if ("P.value" %in% colnames(enrich_df)) {
    "P.value"
  } else {
    # Unexpected column naming — try common alternatives
    intersect(c("p.adjust_hyper", "p.adjust", "pvalue"), colnames(enrich_df))[1]
  }
  if (is.na(p_col)) return("  (unrecognised enrichment result format)")

  # Detect gene column and its separator
  gene_col <- if ("geneID" %in% colnames(enrich_df)) "geneID"  # clusterProfiler "/"
              else if ("Genes" %in% colnames(enrich_df)) "Genes" # enrichr ";"
              else NULL
  gene_sep <- if (!is.null(gene_col) && gene_col == "Genes") ";" else "/"

  # Detect count column
  count_col <- if ("IN" %in% colnames(enrich_df)) "IN"
               else if ("Count" %in% colnames(enrich_df)) "Count"
               else NULL

  # Match contrast: strip trailing _significant if present
  contrast_clean <- gsub("_significant$", "", contrast)
  df <- enrich_df[enrich_df$contrast == contrast_clean, , drop = FALSE]

  if (nrow(df) == 0) {
    return(sprintf(
      "  (no enrichment results for contrast '%s';\n   run Pathway Enrichment or GO analysis first)",
      contrast_clean
    ))
  }

  # Filter by significance
  df <- df[!is.na(df[[p_col]]) & df[[p_col]] < alpha, , drop = FALSE]
  if (nrow(df) == 0) {
    return(sprintf("  (no terms pass adj.p < %.2f threshold)", alpha))
  }

  df <- df[order(df[[p_col]]), ]
  df <- head(df, top_n)

  # Determine database column name
  db_col <- if ("var" %in% colnames(df)) "var" else
            if ("database" %in% colnames(df)) "database" else NULL

  lines <- mapply(function(i) {
    row    <- df[i, , drop = FALSE]
    term   <- as.character(row$Term)
    db_tag <- if (!is.null(db_col)) paste0("[", row[[db_col]], "] ") else ""
    pv     <- formatC(as.numeric(row[[p_col]]), format = "e", digits = 1)

    # Gene ratio / overlap string
    gr <- if ("GeneRatio" %in% colnames(row)) as.character(row$GeneRatio)
          else if ("Overlap"   %in% colnames(row)) as.character(row$Overlap)
          else "N/A"

    # Top overlapping genes
    gene_str <- ""
    if (!is.null(gene_col)) {
      raw <- as.character(row[[gene_col]])
      if (!is.na(raw) && nchar(raw) > 0) {
        genes_vec <- strsplit(raw, gene_sep)[[1]]
        top6      <- head(genes_vec, 6)
        n_total   <- if (!is.null(count_col)) as.integer(row[[count_col]])
                     else length(genes_vec)
        gene_str  <- paste0("\n    Genes: ", paste(top6, collapse = ", "),
                            if (n_total > 6) paste0(" ... (", n_total, " total)") else "")
      }
    }

    sprintf("  %s%s\n    Overlap=%s  adj.p=%s%s", db_tag, term, gr, pv, gene_str)
  }, seq_len(nrow(df)), SIMPLIFY = TRUE)

  paste(lines, collapse = "\n\n")
}


# ------------------------------------------------------------------------------
# Main prompt builder
# ------------------------------------------------------------------------------

#' Build a structured LLM-ready prompt from FragPipe-Analyst results
#'
#' Assembles experiment metadata, ranked DE gene lists, and enrichment results
#' into a single, well-formatted prompt text suitable for any LLM.
#'
#' @param dep           SummarizedExperiment with DE results in rowData
#' @param contrast      Character: which contrast to interpret, e.g. "Drug_A_vs_DMSO"
#' @param enrich_df     Optional: combined enrichment data.frame (GO + pathway results).
#'                      Pass NULL to omit enrichment section.
#' @param include_de    Logical: include DE gene lists
#' @param include_enrich Logical: include enrichment section
#' @param top_n_genes   Max DE genes per direction
#' @param top_n_paths   Max enriched terms
#' @param lfc_thresh    log2FC threshold (used for up/down split annotation)
#' @param alpha         Adjusted p-value cutoff
#' @param user_context  Optional free-text: additional experiment context
#' @param prompt_style  One of "summary", "hypothesis", "methods", "review"
#'
#' @return Character string: the complete prompt text
build_llm_prompt <- function(dep,
                              contrast,
                              enrich_df      = NULL,
                              include_de     = TRUE,
                              include_enrich = TRUE,
                              top_n_genes    = 30,
                              top_n_paths    = 20,
                              lfc_thresh     = 1,
                              alpha          = 0.05,
                              user_context   = "",
                              prompt_style   = "summary") {

  # ---- Experiment metadata ------------------------------------------------
  exp_type   <- metadata(dep)$exp
  data_level <- metadata(dep)$level
  conditions <- paste(sort(unique(as.character(colData(dep)$condition))), collapse = ", ")
  n_features <- nrow(dep)
  n_samples  <- ncol(dep)

  # Parse the contrast string into readable condition names
  parts  <- strsplit(contrast, "_vs_")[[1]]
  cond_A <- parts[1]
  cond_B <- if (length(parts) > 1) paste(parts[-1], collapse = "_vs_") else "reference"

  # Friendly data type label
  dtype_label <- switch(paste(exp_type, data_level, sep = "-"),
    "LFQ-protein"    = "Label-free (LFQ) protein-level",
    "TMT-protein"    = "TMT isobaric protein-level",
    "TMT-gene"       = "TMT isobaric gene-level",
    "TMT-peptide"    = "TMT isobaric peptide-level",
    "TMT-site"       = "TMT isobaric phosphosite-level",
    "DIA-protein"    = "Data-independent acquisition (DIA) protein-level",
    "DIA-peptide"    = "DIA peptide-level",
    "DIA-site"       = "DIA phosphosite-level",
    "LFQ-peptide"    = "Label-free (LFQ) peptide-level",
    paste(exp_type, data_level, sep = " ")  # fallback
  )

  context_block <- if (nchar(trimws(user_context)) > 0) {
    paste0("\nUser-provided context: ", trimws(user_context))
  } else ""

  # ---- Header -------------------------------------------------------------
  sep <- paste(rep("=", 66), collapse = "")

  header <- paste0(
    "You are an expert proteomics data analyst with deep knowledge of\n",
    "cell biology, signaling pathways, and disease mechanisms. Please\n",
    "analyze the following quantitative proteomics results.\n\n",
    sep, "\n",
    "EXPERIMENT OVERVIEW\n",
    sep, "\n",
    "Data type   : ", dtype_label, "\n",
    "Conditions  : ", conditions, "\n",
    "Comparison  : ", cond_A, " vs. ", cond_B, "\n",
    "Features    : ", format(n_features, big.mark = ","),
    " quantified across ", n_samples, " samples\n",
    "Thresholds  : adj.p < ", alpha, ",  |log2FC| > ", lfc_thresh,
    context_block
  )

  # ---- DE section ---------------------------------------------------------
  de_section <- ""
  if (include_de) {
    genes <- extract_de_genes(dep, contrast, lfc_thresh, alpha)
    n_up   <- genes$n_up
    n_down <- genes$n_down

    up_shown   <- min(n_up, top_n_genes)
    down_shown <- min(n_down, top_n_genes)

    up_text   <- format_gene_list(genes$up,   top_n = top_n_genes)
    down_text <- format_gene_list(genes$down, top_n = top_n_genes)

    de_section <- paste0(
      "\n\n", sep, "\n",
      "DIFFERENTIAL EXPRESSION  (", cond_A, " vs. ", cond_B, ")\n",
      sep, "\n",
      "Total significant: ", n_up + n_down,
      "  (", n_up, " up-regulated,  ", n_down, " down-regulated)\n\n",
      "UP-REGULATED in '", cond_A, "'",
      if (n_up > top_n_genes) paste0("  [top ", up_shown, " of ", n_up, " shown]"),
      ":\n",
      up_text,
      "\n\n",
      "DOWN-REGULATED in '", cond_A, "'",
      if (n_down > top_n_genes) paste0("  [top ", down_shown, " of ", n_down, " shown]"),
      ":\n",
      down_text
    )
  }

  # ---- Enrichment section -------------------------------------------------
  enrich_section <- ""
  if (include_enrich && !is.null(enrich_df) && nrow(enrich_df) > 0) {
    # Which databases contributed results for this contrast?
    contrast_clean <- gsub("_significant$", "", contrast)
    available_dbs  <- unique(
      if ("var" %in% colnames(enrich_df)) enrich_df$var[enrich_df$contrast == contrast_clean]
      else if ("database" %in% colnames(enrich_df)) enrich_df$database[enrich_df$contrast == contrast_clean]
      else character(0)
    )
    db_label <- if (length(available_dbs) > 0) paste(available_dbs, collapse = ", ") else "multiple databases"

    enrich_text <- format_enrichment_table(
      enrich_df, contrast,
      top_n  = top_n_paths,
      alpha  = alpha,
      adjust = TRUE
    )

    enrich_section <- paste0(
      "\n\n", sep, "\n",
      "ENRICHMENT ANALYSIS  (", db_label, ")\n",
      sep, "\n",
      enrich_text
    )
  }

  # ---- Task section -------------------------------------------------------
  task_section <- switch(prompt_style,

    "summary" = paste0(
      "\n\n", sep, "\n",
      "YOUR TASK\n",
      sep, "\n",
      "Based on the data above, please provide:\n\n",
      "1. BIOLOGICAL THEMES\n",
      "   Identify the main processes, pathways, or functional modules\n",
      "   represented by the regulated proteins.\n\n",
      "2. KEY PROTEINS\n",
      "   Highlight 3-5 proteins of particular biological or clinical\n",
      "   significance and explain their relevance.\n\n",
      "3. BIOLOGICAL INTERPRETATION\n",
      "   Provide a concise 2-3 paragraph narrative interpreting what\n",
      "   the proteome changes tell us about the biology of '", cond_A,
      "' vs. '", cond_B, "'."
    ),

    "hypothesis" = paste0(
      "\n\n", sep, "\n",
      "YOUR TASK\n",
      sep, "\n",
      "Based on the data above, propose 3 mechanistic hypotheses:\n\n",
      "For each hypothesis:\n",
      "  - State it as a specific, testable claim\n",
      "  - Cite 2-3 supporting proteins or pathways from the data\n",
      "  - Suggest one feasible follow-up experiment to test it\n",
      "  - Note any contradictory findings that need to be reconciled\n\n",
      "Also: identify the most unexpected finding and explain its significance."
    ),

    "methods" = paste0(
      "\n\n", sep, "\n",
      "YOUR TASK\n",
      sep, "\n",
      "Write a Results section paragraph (150-250 words) suitable for\n",
      "a proteomics manuscript describing the key findings from the\n",
      "'", cond_A, "' vs. '", cond_B, "' comparison.\n\n",
      "Requirements:\n",
      "  - Formal scientific language\n",
      "  - Report key findings with exact statistics (fold change, p-value)\n",
      "  - Interpret the biological significance\n",
      "  - Reference enriched pathways when available\n",
      "  - End with a sentence connecting the findings to broader biology"
    ),

    "review" = paste0(
      "\n\n", sep, "\n",
      "YOUR TASK\n",
      sep, "\n",
      "Act as a critical peer reviewer examining these proteomics results:\n\n",
      "1. COHERENCE CHECK\n",
      "   Are the DE findings and enriched pathways internally consistent?\n",
      "   Flag any contradictions.\n\n",
      "2. UNEXPECTED FINDINGS\n",
      "   Identify results that are surprising or hard to interpret.\n\n",
      "3. ALTERNATIVE EXPLANATIONS\n",
      "   Propose at least one alternative biological explanation for\n",
      "   the main finding.\n\n",
      "4. RECOMMENDED FOLLOW-UPS\n",
      "   List 3 specific orthogonal experiments that would validate or\n",
      "   challenge the key conclusions."
    ),

    # Default fallback
    paste0(
      "\n\n", sep, "\n",
      "YOUR TASK\n",
      sep, "\n",
      "Please interpret these proteomics results and provide biological insights."
    )
  )

  # ---- Assemble & return --------------------------------------------------
  paste0(header, de_section, enrich_section, task_section)
}


# ------------------------------------------------------------------------------
# Utility: combine enrichment data frames safely
# ------------------------------------------------------------------------------

#' Safely combine go_results() and pathway_results() into one data frame
#'
#' Both reactives are eventReactive — they throw errors if accessed before
#' their trigger button has been clicked. This function wraps each in
#' tryCatch and returns NULL if neither has been computed yet.
#'
#' @param go_res      Result of go_results() or NULL
#' @param pathway_res Result of pathway_results() or NULL
#' @return Combined data.frame or NULL
combine_enrich_results <- function(go_res, pathway_res) {
  dfs <- list()

  if (!is.null(go_res) && is.data.frame(go_res) && nrow(go_res) > 0) {
    dfs[[length(dfs) + 1]] <- go_res
  }
  if (!is.null(pathway_res) && is.data.frame(pathway_res) && nrow(pathway_res) > 0) {
    dfs[[length(dfs) + 1]] <- pathway_res
  }

  if (length(dfs) == 0) return(NULL)

  # rbind only columns that exist in both — different databases may produce
  # slightly different column sets
  common_cols <- Reduce(intersect, lapply(dfs, colnames))
  do.call(rbind, lapply(dfs, function(d) d[, common_cols, drop = FALSE]))
}


# ------------------------------------------------------------------------------
# Utility: character / token count summary
# ------------------------------------------------------------------------------

#' Return a readable summary of prompt length and rough token estimate
#'
#' @param prompt_text Character string
#' @return Character string like "2,341 characters · ~585 tokens"
prompt_length_summary <- function(prompt_text) {
  if (is.null(prompt_text) || nchar(prompt_text) == 0) return("")
  nch        <- nchar(prompt_text)
  token_est  <- round(nch / 4)
  sprintf("%s characters  ·  ~%s tokens (estimated)",
          format(nch, big.mark = ","),
          format(token_est, big.mark = ","))
}
