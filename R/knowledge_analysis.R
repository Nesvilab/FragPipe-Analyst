###############################################################################
## Knowledge Based Analysis Functions
## Provides: read_gmt, plot_gsva_heatmap, plot_ppi_network,
##           plot_ppi_network_multi, plot_kinase_substrate_network
###############################################################################

# --------------------------------------------------------------------------- #
#  GMT reader
# --------------------------------------------------------------------------- #

#' Read a GMT (Gene Matrix Transposed) file
#'
#' @param file_path Path to a .gmt file
#' @return A named list of character vectors (name → gene symbols)
read_gmt <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  result <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3L) return(NULL)
    genes <- parts[seq(3L, length(parts))]
    genes <- genes[nchar(trimws(genes)) > 0L]
    list(name = parts[1L], genes = genes)
  })
  result <- Filter(Negate(is.null), result)
  setNames(
    lapply(result, `[[`, "genes"),
    vapply(result, `[[`, "name", FUN.VALUE = character(1L))
  )
}

# --------------------------------------------------------------------------- #
#  Internal helper: retrieve msigdbr gene sets as a named list of symbols
# --------------------------------------------------------------------------- #
.get_msig_genesets <- function(database) {
  species <- if (endsWith(database, "(Mouse)")) "Mus musculus" else "Homo sapiens"

  db_map <- list(
    "KEGG"                 = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
    "Reactome"             = list(collection = "C2", subcollection = "CP:REACTOME"),
    "WikiPathways"         = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
    "Hallmark"             = list(collection = "H",  subcollection = NULL),
    "MF"                   = list(collection = "C5", subcollection = "GO:MF"),
    "BP"                   = list(collection = "C5", subcollection = "GO:BP"),
    "CC"                   = list(collection = "C5", subcollection = "GO:CC"),
    "KEGG (Mouse)"         = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
    "WikiPathways (Mouse)" = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS")
  )

  info <- db_map[[database]]
  if (is.null(info)) stop("Unknown database: ", database)

  df <- if (is.null(info$subcollection)) {
    msigdbr::msigdbr(species = species, collection = info$collection)
  } else {
    msigdbr::msigdbr(species = species, collection = info$collection,
                     subcollection = info$subcollection)
  }

  split(df$gene_symbol, df$gs_name)
}


# --------------------------------------------------------------------------- #
#  ORA Heatmap  (rows = terms, cols = contrasts)
# --------------------------------------------------------------------------- #

#' ORA Pathway Heatmap across all contrasts
#'
#' Takes the data frame returned by \code{test_ora_mod()} and draws a
#' ComplexHeatmap where rows = pathway terms and columns = contrasts.
#'
#' @param ora_results  data.frame from \code{test_ora_mod()}
#' @param top_n        Maximum number of terms to show (selected by smallest
#'                     p_hyper across any contrast)
#' @param value_type   One of \code{"log2OR"} (log2 odds ratio),
#'                     \code{"pvalue"} (-log10 hypergeometric p),
#'                     or \code{"size"} (hit count)
#' @param alpha        p-value cutoff — only terms significant in at least one
#'                     contrast are shown
#' @param use_adjp     If TRUE, use adjusted p-value for significance filter
#' @return ComplexHeatmap object (call \code{draw()} inside renderPlot)
plot_ora_heatmap <- function(ora_results, value_type = "log2OR",
                             alpha = 0.05, use_adjp = TRUE) {

  required_cols <- c("Term", "contrast", "log_odds", "p_hyper", "IN")
  missing <- setdiff(required_cols, colnames(ora_results))
  if (length(missing) > 0)
    stop("ORA results missing columns: ", paste(missing, collapse = ", "))

  # ---- 1. Select terms significant in at least one contrast --------------- #

  # If direction column exists (UP/DOWN), make Term unique per direction
  has_direction <- "direction" %in% colnames(ora_results)
  if (has_direction) {
    ora_results <- ora_results %>%
      dplyr::mutate(Term = paste0(Term, " (", direction, ")"))
  }

  p_col <- if (use_adjp && "p.adjust_hyper" %in% colnames(ora_results))
    "p.adjust_hyper" else "p_hyper"

  sig_terms <- ora_results %>%
    dplyr::filter(.data[[p_col]] < alpha) %>%
    dplyr::pull(Term) %>%
    unique()

  if (length(sig_terms) == 0)
    stop("No ORA terms with ", if (use_adjp) "adj. " else "",
         "p < ", alpha, " in any contrast.",
         "\nTry raising the cutoff or choosing a different database.")

  df_sub <- ora_results %>% dplyr::filter(Term %in% sig_terms)

  # ---- 2. Pivot to wide matrix -------------------------------------------- #
  if (value_type == "log2OR") {
    fill_label <- "log2 OR"
    fill_na    <- 0
    df_wide    <- df_sub %>%
      dplyr::select(Term, contrast, val = log_odds) %>%
      tidyr::pivot_wider(names_from  = contrast,
                         values_from = val,
                         values_fill = fill_na)
  } else if (value_type == "pvalue") {
    fill_label <- "-log10(p)"
    fill_na    <- 0
    df_wide    <- df_sub %>%
      dplyr::mutate(val = -log10(pmax(p_hyper, 1e-300))) %>%
      dplyr::select(Term, contrast, val) %>%
      tidyr::pivot_wider(names_from  = contrast,
                         values_from = val,
                         values_fill = fill_na)
  } else {  # "size"
    fill_label <- "Hits (size)"
    fill_na    <- 0L
    df_wide    <- df_sub %>%
      dplyr::select(Term, contrast, val = IN) %>%
      tidyr::pivot_wider(names_from  = contrast,
                         values_from = val,
                         values_fill = fill_na)
  }

  mat          <- as.matrix(df_wide[, -1, drop = FALSE])
  rownames(mat) <- df_wide$Term

  # Shorten long term names for readability
  rownames(mat) <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOMF_|^GOCC_",
    "", rownames(mat)
  )
  rownames(mat) <- gsub("_", " ", rownames(mat))

  # ---- 3. Color scale ----------------------------------------------------- #
  if (value_type == "log2OR") {
    max_val <- max(abs(mat), na.rm = TRUE)
    if (!is.finite(max_val) || max_val == 0) max_val <- 1
    col_fun <- circlize::colorRamp2(
      c(-max_val, 0, max_val),
      c("#3498db", "white", "#e74c3c")
    )
  } else {
    max_val <- max(mat, na.rm = TRUE)
    if (!is.finite(max_val) || max_val == 0) max_val <- 1
    col_fun <- circlize::colorRamp2(
      c(0, max_val),
      c("white", "#e74c3c")
    )
  }

  # ---- 4. Draw heatmap ---------------------------------------------------- #
  ComplexHeatmap::Heatmap(
    mat,
    name              = fill_label,
    col               = col_fun,
    cluster_rows      = TRUE,
    cluster_columns   = FALSE,        # preserve contrast order
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    row_names_max_width = grid::unit(12, "cm"),
    na_col            = "grey90",
    rect_gp           = grid::gpar(col = "white", lwd = 0.5),
    column_title      = paste0("ORA Heatmap \u2013 ", fill_label),
    column_title_gp   = grid::gpar(fontsize = 11, fontface = "bold")
  )
}


# --------------------------------------------------------------------------- #
#  Helpers for knowledge-based ORA across all contrasts (GMT-aware)
# --------------------------------------------------------------------------- #

# Extract gene symbols from dep for background (contrast=NULL) or a contrast.
.kb_extract_genes <- function(dep, contrast = NULL, direction = "UP",
                               log2_threshold = 0, alpha = 0.05,
                               adjust_alpha = FALSE) {
  rd  <- as.data.frame(SummarizedExperiment::rowData(dep))
  exp <- SummarizedExperiment::metadata(dep)$exp
  lvl <- SummarizedExperiment::metadata(dep)$level

  get_sym <- function(df) {
    if (!is.null(exp) && !is.null(lvl)) {
      if (exp == "TMT" && lvl == "gene")    return(unique(df$ID))
      if (exp == "TMT" && lvl == "protein") return(unique(gsub("[.].*", "", df$Gene)))
    }
    if ("Gene" %in% colnames(df)) return(unique(df$Gene))
    unique(gsub("[.].*", "", df$name))
  }

  background <- get_sym(rd)
  background <- background[!is.na(background) & nchar(trimws(background)) > 0]

  if (is.null(contrast)) return(list(background = background))

  diff_col <- paste0(contrast, "_diff")
  p_col    <- if (adjust_alpha) paste0(contrast, "_p.adj") else paste0(contrast, "_p.val")
  if (!diff_col %in% colnames(rd) || !p_col %in% colnames(rd))
    return(list(background = background, sig_genes = character(0)))

  p_ok <- !is.na(rd[[p_col]]) & rd[[p_col]] < alpha
  dir_ok <- if (direction == "UP") {
    !is.na(rd[[diff_col]]) & rd[[diff_col]] > log2_threshold
  } else if (direction == "DOWN") {
    !is.na(rd[[diff_col]]) & rd[[diff_col]] < -log2_threshold
  } else {
    !is.na(rd[[diff_col]]) & abs(rd[[diff_col]]) > log2_threshold
  }

  sig_genes <- get_sym(rd[p_ok & dir_ok, , drop = FALSE])
  sig_genes <- sig_genes[!is.na(sig_genes) & nchar(trimws(sig_genes)) > 0]
  list(background = background, sig_genes = sig_genes)
}

# Run ORA using a custom GMT collection across every contrast in dep.
.run_ora_gmt_all <- function(dep, gmt_list, gmt_name, alpha = 0.05,
                              selected_contrasts = NULL) {
  term2gene <- data.frame(
    term = rep(names(gmt_list), lengths(gmt_list)),
    gene = unlist(gmt_list, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  rd            <- as.data.frame(SummarizedExperiment::rowData(dep))
  contrast_cols <- grep("_significant$", colnames(rd), value = TRUE)
  contrasts     <- gsub("_significant$", "", contrast_cols)

  if (!is.null(selected_contrasts) && length(selected_contrasts) > 0)
    contrasts <- intersect(contrasts, selected_contrasts)

  bg_info <- .kb_extract_genes(dep)
  bg      <- bg_info$background
  bg_n    <- length(bg)

  bg_enrich <- tryCatch(
    suppressMessages(clusterProfiler::enricher(
      gene = bg, TERM2GENE = term2gene,
      pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3)),
    error = function(e) NULL)

  bg_df <- if (!is.null(bg_enrich) && nrow(as.data.frame(bg_enrich)) > 0) {
    as.data.frame(bg_enrich) %>%
      dplyr::mutate(bg_IN = Count, bg_OUT = bg_n - Count) %>%
      dplyr::select(ID, bg_IN, bg_OUT)
  } else {
    data.frame(ID = character(0), bg_IN = integer(0), bg_OUT = integer(0))
  }

  results <- lapply(contrasts, function(ct) {
    info    <- .kb_extract_genes(dep, contrast = ct, direction = "UP",
                                 log2_threshold = 0, alpha = alpha,
                                 adjust_alpha   = FALSE)
    genes   <- info$sig_genes
    n_genes <- length(genes)
    if (n_genes == 0) return(NULL)

    res <- tryCatch(
      suppressMessages(clusterProfiler::enricher(
        gene = genes, TERM2GENE = term2gene,
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 3)),
      error = function(e) NULL)
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)

    df %>%
      dplyr::rename(Term = Description, IN = Count) %>%
      dplyr::mutate(contrast = ct, OUT = n_genes - IN) %>%
      dplyr::left_join(bg_df, by = "ID") %>%
      dplyr::mutate(
        bg_IN          = dplyr::if_else(is.na(bg_IN),  1L, as.integer(bg_IN)),
        bg_OUT         = dplyr::if_else(is.na(bg_OUT), bg_n - 1L, as.integer(bg_OUT)),
        log_odds       = log2((IN * bg_OUT) / (OUT * bg_IN)),
        p_hyper        = phyper(q = IN - 1, m = bg_IN, n = bg_OUT,
                                k = IN + OUT, lower.tail = FALSE),
        p.adjust_hyper = p.adjust(p_hyper, method = "BH"),
        var            = gmt_name
      )
  })
  dplyr::bind_rows(Filter(Negate(is.null), results))
}

#' Run ORA for a single contrast + single database (Gene Set Explorer)
#'
#' @param dep            SummarizedExperiment after test_diff()
#' @param database       Database name (built-in or GMT collection name)
#' @param contrast       Single contrast name (e.g. "A_vs_B")
#' @param direction      "UP", "DOWN", or "both"
#' @param alpha          p-value cutoff for selecting DE genes
#' @param log2_threshold log2 fold-change cutoff for DE gene selection
#' @param adjust_alpha   Use adjusted p-values for DE selection?
#' @param backend        "clusterProfiler" (default) or "enrichr"
#' @return data.frame with ORA results; has a 'direction' column when direction="both"
run_ora_kb_single <- function(dep, database, contrast,
                               direction = "both", alpha = 0.05,
                               log2_threshold = 1, adjust_alpha = TRUE,
                               backend = "clusterProfiler") {
  loaded_gmt <- tryCatch(
    get("LOADED_GMT_FILES", envir = .GlobalEnv, inherits = FALSE),
    error = function(e) list()
  )
  is_gmt <- database %in% names(loaded_gmt)

  run_one_dir <- function(dir) {
    if (is_gmt) {
      # Custom GMT
      term2gene <- data.frame(
        term = rep(names(loaded_gmt[[database]]), lengths(loaded_gmt[[database]])),
        gene = unlist(loaded_gmt[[database]], use.names = FALSE),
        stringsAsFactors = FALSE
      )
      info  <- .kb_extract_genes(dep, contrast = contrast, direction = dir,
                                  log2_threshold = log2_threshold, alpha = alpha,
                                  adjust_alpha = adjust_alpha)
      genes <- info$sig_genes
      bg    <- info$background
      if (length(genes) == 0) return(NULL)
      res <- tryCatch(
        suppressWarnings(suppressMessages(
          clusterProfiler::enricher(gene = genes, universe = bg,
                                     TERM2GENE = term2gene,
                                     pvalueCutoff = 1, qvalueCutoff = 1,
                                     minGSSize = 3)
        )), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      bg_IN  <- as.numeric(gsub("/.*", "", df$BgRatio))
      bg_OUT <- as.numeric(gsub(".*/", "", df$BgRatio)) - bg_IN
      df %>%
        dplyr::rename(Term = Description, IN = Count) %>%
        dplyr::mutate(
          contrast       = contrast,
          OUT            = length(genes) - IN,
          Odds.Ratio     = (IN * bg_OUT) / (OUT * bg_IN),
          log_odds       = log2(Odds.Ratio),
          p_hyper        = pvalue,
          p.adjust_hyper = p.adjust,
          var            = database
        )
    } else {
      # Built-in database via test_ora_mod (all contrasts at once)
      res <- test_ora_mod(dep, databases = database, contrasts = TRUE,
                          direction = dir, log2_threshold = log2_threshold,
                          alpha = alpha, adjust_alpha = adjust_alpha,
                          backend = backend)
      if (is.null(res) || nrow(res) == 0) return(NULL)
      res
    }
  }

  if (direction == "both") {
    up   <- run_one_dir("UP")
    down <- run_one_dir("DOWN")
    if (!is.null(up))   up$direction   <- "UP"
    if (!is.null(down))  down$direction <- "DOWN"
    dplyr::bind_rows(up, down)
  } else {
    res <- run_one_dir(direction)
    if (!is.null(res)) res$direction <- direction
    res
  }
}

#' Run ORA across all contrasts — built-in databases or custom GMT files
#'
#' @param dep               SummarizedExperiment after test_diff()
#' @param database          Database name (built-in or GMT collection name)
#' @param alpha             p-value threshold for the foreground gene list
#' @param selected_contrasts Optional character vector of contrast names to keep
#' @param backend           "clusterProfiler" (default) or "enrichr"
#' @return data.frame compatible with \code{plot_ora_heatmap()}
run_ora_kb <- function(dep, database, alpha = 0.05, selected_contrasts = NULL,
                       backend = "clusterProfiler") {
  loaded_gmt <- tryCatch(
    get("LOADED_GMT_FILES", envir = .GlobalEnv, inherits = FALSE),
    error = function(e) list()
  )
  if (database %in% names(loaded_gmt)) {
    .run_ora_gmt_all(dep, gmt_list = loaded_gmt[[database]],
                     gmt_name = database, alpha = alpha,
                     selected_contrasts = selected_contrasts)
  } else {
    res <- test_ora_mod(dep, databases = database, contrasts = TRUE,
                        direction = "UP", log2_threshold = 0,
                        alpha = alpha, adjust_alpha = FALSE,
                        backend = backend)
    if (!is.null(selected_contrasts) && length(selected_contrasts) > 0 &&
        "contrast" %in% colnames(res))
      res <- dplyr::filter(res, contrast %in% selected_contrasts)
    res
  }
}

#' Pre-compute ORA across every built-in database and loaded GMT file
#'
#' Called when the Knowledge Based Analysis tab is first entered (or on
#' manual refresh).  Returns a named list where each element is the ORA
#' result data.frame (or NULL on failure) for that database.
#'
#' @param dep      SummarizedExperiment after test_diff()
#' @param alpha    p-value threshold for foreground gene selection
#' @param backend  "clusterProfiler" or "enrichr"
#' @param progress Optional shiny Progress object; incremented after each db
#' @return Named list: database_name -> data.frame | NULL
run_ora_kb_all <- function(dep, alpha = 0.05, backend = "clusterProfiler",
                           progress = NULL) {
  builtin_dbs <- c("KEGG", "Reactome", "WikiPathways", "Hallmark",
                   "MF", "BP", "CC")
  loaded_gmt  <- tryCatch(
    get("LOADED_GMT_FILES", envir = .GlobalEnv, inherits = FALSE),
    error = function(e) list()
  )
  all_dbs <- c(builtin_dbs, names(loaded_gmt))
  n       <- length(all_dbs)

  results <- setNames(vector("list", n), all_dbs)
  t_total <- proc.time()[["elapsed"]]
  for (i in seq_along(all_dbs)) {
    db <- all_dbs[[i]]
    if (!is.null(progress))
      progress$inc(1 / n, detail = paste("Running", db, "\u2026"))
    t0 <- proc.time()[["elapsed"]]
    results[[db]] <- tryCatch(
      suppressWarnings(suppressMessages(
        run_ora_kb(dep, database = db, alpha = alpha, backend = backend)
      )),
      error = function(e) { message("ORA skipped [", db, "]: ", e$message); NULL }
    )
    t1 <- proc.time()[["elapsed"]]
    nrows <- if (!is.null(results[[db]])) nrow(results[[db]]) else 0L
    message(sprintf("[ORA] %-20s %6.1fs  (%d results)", db, t1 - t0, nrows))
  }
  message(sprintf("[ORA] TOTAL %43.1fs", proc.time()[["elapsed"]] - t_total))
  results
}


# --------------------------------------------------------------------------- #
#  PTM-SEA (PTM Signature Enrichment Analysis) via ssGSEA2
# --------------------------------------------------------------------------- #

#' Get path to the appropriate PTMsigDB GMT file
#'
#' @param species One of "human", "mouse", "rat"
#' @return File path to the GMT file
get_ptmsigdb_path <- function(species = "human") {
  species <- tolower(species)
  if (!species %in% c("human", "mouse", "rat")) species <- "human"
  file.path("data", "gene_sets", "ptmsigdb",
            paste0("ptm.sig.db.all.flanking.", species, ".v2.0.0.gmt"))
}

#' Build a named stat vector for PTM-SEA from DE results
#'
#' Extracts sign(logFC) * -log10(p) per site, names by flanking sequence,
#' deduplicates, and returns a clean named numeric vector.
#'
#' @param rd       Row data data.frame from SummarizedExperiment
#' @param contrast Contrast name
#' @param mod_type Modification type suffix (e.g. "p", "ac")
#' @return Named numeric vector (flanking_seq-mod_type → stat)
.build_ptmsea_stat <- function(rd, contrast, mod_type) {
  diff_col <- paste0(contrast, "_diff")
  p_col    <- paste0(contrast, "_p.val")
  if (!diff_col %in% colnames(rd) || !p_col %in% colnames(rd))
    stop("Contrast columns not found for: ", contrast)

  lfc  <- rd[[diff_col]]
  pval <- rd[[p_col]]
  seqw <- rd$SequenceWindow

  pval[pval < 1e-300] <- 1e-300
  stat <- sign(lfc) * (-log10(pval))

  mod_suffix <- paste0("-", mod_type)
  names(stat) <- paste0(toupper(seqw), mod_suffix)

  keep <- !is.na(stat) & !is.na(seqw) & nchar(trimws(seqw)) > 0
  stat <- stat[keep]

  if (anyDuplicated(names(stat))) {
    df_dedup <- data.frame(nm = names(stat), val = stat, absval = abs(stat),
                           stringsAsFactors = FALSE)
    df_dedup <- df_dedup[order(-df_dedup$absval), ]
    df_dedup <- df_dedup[!duplicated(df_dedup$nm), ]
    stat <- setNames(df_dedup$val, df_dedup$nm)
  }
  stat
}

#' Write a stat vector to a GCT v1.3 file for ssGSEA2
#'
#' @param stat     Named numeric vector (site IDs → values)
#' @param gct_path Output file path
#' @return The file path (invisibly)
.write_stat_gct <- function(stat, gct_path) {
  mat <- matrix(stat, ncol = 1, dimnames = list(names(stat), "sample"))
  gct <- new("GCT", mat = mat)
  cmapR::write_gct(gct, gct_path, appenddim = FALSE)
  invisible(gct_path)
}

#' Parse ssGSEA2 output GCT files into a tidy data.frame
#'
#' Reads the scores, p-values, and FDR GCT files produced by run_ssGSEA2()
#' and combines them into a long-format data.frame. Supports multi-column
#' (multi-contrast) output.
#'
#' @param output_prefix The prefix used for run_ssGSEA2 output
#' @param output_dir    The directory containing output files
#' @return data.frame with columns: set, set_size, ES, NES, p_value,
#'         adj_p_value, contrast
.parse_ssgsea2_output <- function(output_prefix, output_dir) {
  scores_path   <- file.path(output_dir, paste0(output_prefix, "-scores.gct"))
  pvals_path    <- file.path(output_dir, paste0(output_prefix, "-pvalues.gct"))
  fdr_path      <- file.path(output_dir, paste0(output_prefix, "-fdr-pvalues.gct"))
  combined_path <- file.path(output_dir, paste0(output_prefix, "-combined.gct"))

  if (!file.exists(scores_path))
    stop("ssGSEA2 scores output not found: ", scores_path)

  scores_gct <- cmapR::parse_gctx(scores_path)
  scores_mat <- scores_gct@mat  # rows = gene sets, cols = contrasts

  pvals_mat <- if (file.exists(pvals_path)) {
    cmapR::parse_gctx(pvals_path)@mat
  } else {
    matrix(NA_real_, nrow = nrow(scores_mat), ncol = ncol(scores_mat),
           dimnames = dimnames(scores_mat))
  }

  fdr_mat <- if (file.exists(fdr_path)) {
    cmapR::parse_gctx(fdr_path)@mat
  } else {
    matrix(NA_real_, nrow = nrow(scores_mat), ncol = ncol(scores_mat),
           dimnames = dimnames(scores_mat))
  }

  # Try to extract set_size from combined GCT row annotations
  set_sizes <- if (file.exists(combined_path)) {
    combined_gct <- cmapR::parse_gctx(combined_path)
    rdesc <- combined_gct@rdesc
    if ("n.geneset.genes" %in% colnames(rdesc)) {
      as.integer(rdesc[["n.geneset.genes"]])
    } else if ("Geneset.size" %in% colnames(rdesc)) {
      as.integer(rdesc[["Geneset.size"]])
    } else {
      rep(NA_integer_, nrow(scores_mat))
    }
  } else {
    rep(NA_integer_, nrow(scores_mat))
  }

  # Build long-format data.frame across all columns (contrasts)
  contrast_names <- colnames(scores_mat)
  all_res <- lapply(seq_along(contrast_names), function(j) {
    data.frame(
      set         = rownames(scores_mat),
      set_size    = set_sizes,
      ES          = scores_mat[, j],
      NES         = scores_mat[, j],
      p_value     = pvals_mat[, j],
      adj_p_value = fdr_mat[, j],
      contrast    = contrast_names[j],
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(all_res)
}

#' Run PTM-SEA across all contrasts using ssGSEA2
#'
#' Builds stat vectors for all contrasts, writes them as columns in a single
#' GCT file, and runs ssGSEA2 once with parallel gene-set processing.
#'
#' @param dep      SummarizedExperiment after test_diff() (site-level)
#' @param species  "human", "mouse", or "rat"
#' @param mod_type PTM modification type suffix for PTMsigDB matching.
#' @param nperm    Number of permutations for p-value estimation
#' @param min_size Minimum gene set size (min.overlap in ssGSEA2)
#' @param subsets  Character vector of PTMsigDB subset prefixes to include.
#'   Valid values: "PERT", "PATH", "DISEASE", "KINASE". If NULL, all are used.
#' @return data.frame with columns: contrast, set, set_size, ES, NES,
#'         p_value, adj_p_value
run_ptmsea_all <- function(dep, species = "human", mod_type = "p",
                           nperm = 1000L, min_size = 5L, subsets = NULL) {
  rd <- as.data.frame(SummarizedExperiment::rowData(dep))

  if (!"SequenceWindow" %in% colnames(rd))
    stop("SequenceWindow column missing — PTM-SEA requires site-level data.")

  all_contrasts <- gsub("_significant$", "",
                        grep("_significant$", colnames(rd), value = TRUE))
  if (length(all_contrasts) == 0)
    stop("No contrast columns found in data.")

  gmt_path <- get_ptmsigdb_path(species)
  if (!file.exists(gmt_path))
    stop("PTMsigDB file not found: ", gmt_path)

  # Filter GMT to selected subsets by writing a filtered GMT file
  # PTMsigDB set names are prefixed: KINASE-*, PATH-*, PERT-*, DISEASE-*
  if (!is.null(subsets) && length(subsets) > 0 &&
      length(subsets) < 4) {
    all_lines <- readLines(gmt_path, warn = FALSE)
    prefix_pattern <- paste0("^(", paste(subsets, collapse = "|"), ")-")
    keep <- grepl(prefix_pattern, all_lines)
    if (sum(keep) == 0)
      stop("No gene sets match the selected subsets: ",
           paste(subsets, collapse = ", "))
    filtered_gmt <- tempfile(fileext = ".gmt")
    writeLines(all_lines[keep], filtered_gmt)
    gmt_path <- filtered_gmt
    on.exit(unlink(filtered_gmt), add = TRUE)
    message(sprintf("[PTM-SEA] Filtered to %d gene sets (subsets: %s)",
                    sum(keep), paste(subsets, collapse = ", ")))
  }

  # Build stat vectors for all contrasts and collect into a matrix
  stat_list <- list()
  for (contrast in all_contrasts) {
    diff_col <- paste0(contrast, "_diff")
    p_col    <- paste0(contrast, "_p.val")
    if (!diff_col %in% colnames(rd) || !p_col %in% colnames(rd)) next

    stat <- tryCatch(.build_ptmsea_stat(rd, contrast, mod_type),
                     error = function(e) { message("[PTM-SEA] ", e$message); NULL })
    if (is.null(stat)) next

    if (length(stat) < 10) {
      message("[PTM-SEA] Skipping ", contrast, ": too few valid sites (",
              length(stat), ")")
      next
    }
    stat_list[[contrast]] <- stat
  }

  if (length(stat_list) == 0)
    stop("No contrasts had enough valid sites for PTM-SEA.")

  # Merge stat vectors into a single matrix (union of all site IDs)
  all_sites <- unique(unlist(lapply(stat_list, names)))
  mat <- matrix(NA_real_, nrow = length(all_sites), ncol = length(stat_list),
                dimnames = list(all_sites, names(stat_list)))
  for (cn in names(stat_list)) {
    s <- stat_list[[cn]]
    mat[names(s), cn] <- s
  }
  # Replace NA with 0 (sites not in a contrast get neutral stat)
  mat[is.na(mat)] <- 0

  # Write multi-column GCT
  tmp_dir <- tempfile("ptmsea_all_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  gct_path <- file.path(tmp_dir, "input.gct")
  t0 <- proc.time()[["elapsed"]]
  gct <- new("GCT", mat = mat)
  cmapR::write_gct(gct, gct_path, appenddim = FALSE)
  t_gct <- proc.time()[["elapsed"]] - t0
  message(sprintf("[PTM-SEA] Writing GCT (%d sites x %d contrasts): %.1fs",
                  nrow(mat), ncol(mat), t_gct))

  # Run ssGSEA2 once for all contrasts, with parallel gene-set processing
  t0 <- proc.time()[["elapsed"]]
  tryCatch(
    run_ssGSEA2(
      input.ds           = gct_path,
      output.prefix      = "ptmsea",
      gene.set.databases = gmt_path,
      output.directory   = tmp_dir,
      sample.norm.type   = "none",
      weight             = 0.75,
      statistic          = "area.under.RES",
      output.score.type  = "NES",
      nperm              = as.integer(nperm),
      min.overlap        = as.integer(min_size),
      correl.type        = "rank",
      global.fdr         = FALSE,
      par                = TRUE,
      spare.cores        = 1,
      export.signat.gct  = FALSE,
      param.file         = FALSE,
      log.file           = file.path(tmp_dir, "run.log")
    ),
    error = function(e) stop("ssGSEA2 failed: ", e$message)
  )
  t_sea <- proc.time()[["elapsed"]] - t0
  message(sprintf("[PTM-SEA] ssGSEA2 computation: %.1fs", t_sea))

  .parse_ssgsea2_output("ptmsea", tmp_dir)
}

#' Parse a PTMsigDB gene set into expected-up and expected-down site lists
#'
#' @param gmt_list  Named list from read_gmt()
#' @param set_name  Name of the gene set to parse
#' @param mod_type  Modification type to filter for (e.g. "p", "ac", "ub").
#'                  If NULL, all modification types are included.
#' @return List with $up and $down character vectors of bare flanking sequences
parse_ptmsea_set_sites <- function(gmt_list, set_name, mod_type = NULL) {
  members <- gmt_list[[set_name]]
  if (is.null(members)) return(list(up = character(0), down = character(0)))

  # Filter by modification type if specified (e.g. keep only "-p;u"/"-p;d")
  if (!is.null(mod_type)) {
    mod_pattern <- paste0("-", mod_type, ";[ud]$")
    members <- members[grepl(mod_pattern, members)]
  }

  is_up   <- grepl(";u$", members)
  is_down <- grepl(";d$", members)
  # Strip the direction and modification suffixes (e.g., "-p;u" -> bare flanking)
  strip <- function(x) toupper(gsub("-[a-z]+;[ud]$", "", x))

  list(
    up   = strip(members[is_up]),
    down = strip(members[is_down])
  )
}

#' GSVA Pathway Heatmap
#'
#' Runs GSVA on the normalised expression matrix of \code{dep} and displays
#' enrichment scores as a heatmap.  Rows = gene sets, columns = samples.
#' A condition annotation bar is drawn at the top of the heatmap.
#'
#' @param dep             SummarizedExperiment after test_diff()
#' @param database        One of "KEGG", "Reactome", "WikiPathways", "Hallmark",
#'                        "MF", "BP", "CC", "KEGG (Mouse)", "WikiPathways (Mouse)"
#' @param top_n           Number of most-variable gene sets to display
#' @param order_by_condition If TRUE, sort columns by condition; otherwise
#'                        cluster columns hierarchically
#' @return ComplexHeatmap object (call draw() or print() in renderPlot)
plot_gsva_heatmap <- function(dep, database = "Hallmark", gene_sets = NULL,
                              top_n = 25, order_by_condition = TRUE) {

  # ---- 1. Expression matrix (genes x samples) ----------------------------- #
  expr_mat <- SummarizedExperiment::assay(dep)
  col_data <- as.data.frame(SummarizedExperiment::colData(dep))
  row_data <- as.data.frame(SummarizedExperiment::rowData(dep))

  # Map rows to gene symbols
  if ("Gene" %in% colnames(row_data)) {
    gene_syms <- row_data$Gene
  } else {
    gene_syms <- gsub("[.].*", "", row_data$name)
  }

  # Keep only rows with a valid symbol
  valid <- !is.na(gene_syms) & nchar(gene_syms) > 0
  expr_mat  <- expr_mat[valid, , drop = FALSE]
  gene_syms <- gene_syms[valid]
  rownames(expr_mat) <- gene_syms

  # For duplicate symbols, retain the one with highest inter-sample variance
  gene_var  <- apply(expr_mat, 1, var, na.rm = TRUE)
  expr_mat  <- expr_mat[order(-gene_var), , drop = FALSE]
  expr_mat  <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

  # Impute remaining NAs per row with the row median
  expr_mat <- t(apply(expr_mat, 1, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  }))

  # ---- 2. Gene sets -------------------------------------------------------- #
  if (is.null(gene_sets)) gene_sets <- .get_msig_genesets(database)

  # Restrict to gene sets with >= 5 expressed genes
  expressed <- rownames(expr_mat)
  gene_sets <- gene_sets[
    vapply(gene_sets, function(gs) length(intersect(gs, expressed)) >= 5L, logical(1))
  ]
  if (length(gene_sets) == 0)
    stop("No gene sets with >= 5 overlapping genes found for the selected database.")

  # ---- 3. GSVA ------------------------------------------------------------ #
  param     <- GSVA::gsvaParam(exprData = expr_mat, geneSets = gene_sets,
                               kcdf = "Gaussian", minSize = 5L)
  gsva_mat  <- GSVA::gsva(param, verbose = FALSE)

  # ---- 4. Top-N most variable gene sets ----------------------------------- #
  row_vars  <- apply(gsva_mat, 1, var, na.rm = TRUE)
  top_idx   <- order(-row_vars)[seq_len(min(top_n, nrow(gsva_mat)))]
  gsva_mat  <- gsva_mat[top_idx, , drop = FALSE]

  # Clean pathway name prefixes for readability
  rownames(gsva_mat) <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOMF_|^GOCC_",
    "", rownames(gsva_mat)
  )
  rownames(gsva_mat) <- gsub("_", " ", rownames(gsva_mat))

  # ---- 5. Column ordering ------------------------------------------------- #
  condition <- col_data$condition
  if (is.null(condition)) condition <- rep("Sample", ncol(gsva_mat))

  if (order_by_condition) {
    col_order <- order(condition)
    gsva_mat  <- gsva_mat[, col_order, drop = FALSE]
    condition <- condition[col_order]
  }

  # ---- 6. Condition annotation bar ---------------------------------------- #
  cond_levels <- unique(condition)
  n_cond      <- length(cond_levels)
  palette     <- if (n_cond <= 8) {
    RColorBrewer::brewer.pal(max(3L, n_cond), "Set2")[seq_len(n_cond)]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_cond)
  }
  cond_colors <- stats::setNames(palette, cond_levels)

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition = condition,
    col       = list(Condition = cond_colors),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = 9)
  )

  # ---- 7. Color scale ----------------------------------------------------- #
  max_val <- max(abs(gsva_mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  col_fun <- circlize::colorRamp2(
    c(-max_val, 0, max_val),
    c("#3498db", "white", "#e74c3c")
  )

  # ---- 8. Draw heatmap ---------------------------------------------------- #
  ComplexHeatmap::Heatmap(
    gsva_mat,
    name              = "GSVA score",
    col               = col_fun,
    top_annotation    = top_anno,
    cluster_rows      = TRUE,
    cluster_columns   = TRUE,
    show_column_names = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 9),
    column_names_rot  = 45,
    row_names_max_width = grid::unit(12, "cm"),
    na_col            = "grey90",
    rect_gp           = grid::gpar(col = "white", lwd = 0.5),
    column_title      = paste0("GSVA \u2013 ", database),
    column_title_gp   = grid::gpar(fontsize = 11, fontface = "bold")
  )
}


#' Query STRING REST API for protein interactions
#'
#' @param gene_symbols Character vector of gene symbols
#' @param species_id   NCBI taxonomy ID (9606 = human, 10090 = mouse)
#' @param string_score Minimum combined STRING score (0–1, e.g. 0.4)
#' @return list with $edges data.frame (from, to, score, evidence columns)
#'         and $api_error character or NULL
.query_string_api <- function(gene_symbols, species_id = 9606,
                              string_score = 0.4) {
  identifiers <- paste(gene_symbols, collapse = "%0d")
  api_score <- as.integer(string_score * 1000)  # API expects 0-1000
  api_url <- paste0(
    "https://string-db.org/api/json/network",
    "?identifiers=", identifiers,
    "&species=", species_id,
    "&required_score=", api_score,
    "&caller_identity=FragPipe-Analyst"
  )

  resp <- tryCatch(
    httr::GET(api_url, httr::timeout(30)),
    error = function(e) list(error = e$message)
  )

  if (is.list(resp) && !is.null(resp$error))
    return(list(edges = data.frame(), api_error = resp$error))
  if (httr::status_code(resp) != 200)
    return(list(edges = data.frame(),
                api_error = paste0("HTTP ", httr::status_code(resp))))

  net_json <- rjson::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
  if (length(net_json) == 0)
    return(list(edges = data.frame(), api_error = NULL))

  # Evidence score labels
  ev_labels <- c(nscore = "neighborhood", fscore = "gene fusion",
                 pscore = "phylogenetic", ascore = "coexpression",
                 escore = "experimental", dscore = "database",
                 tscore = "textmining")

  .safe_num <- function(x) if (is.null(x)) 0 else as.numeric(x)
  edges <- data.frame(
    from  = vapply(net_json, function(x) x$preferredName_A, ""),
    to    = vapply(net_json, function(x) x$preferredName_B, ""),
    score = vapply(net_json, function(x) .safe_num(x$score), 0),
    escore = vapply(net_json, function(x) .safe_num(x$escore), 0),
    dscore = vapply(net_json, function(x) .safe_num(x$dscore), 0),
    tscore = vapply(net_json, function(x) .safe_num(x$tscore), 0),
    ascore = vapply(net_json, function(x) .safe_num(x$ascore), 0),
    nscore = vapply(net_json, function(x) .safe_num(x$nscore), 0),
    fscore = vapply(net_json, function(x) .safe_num(x$fscore), 0),
    pscore = vapply(net_json, function(x) .safe_num(x$pscore), 0),
    stringsAsFactors = FALSE
  )

  # Build human-readable evidence tooltip
  edges$evidence <- vapply(seq_len(nrow(edges)), function(i) {
    scores <- c(escore = edges$escore[i], dscore = edges$dscore[i],
                tscore = edges$tscore[i], ascore = edges$ascore[i],
                nscore = edges$nscore[i], fscore = edges$fscore[i],
                pscore = edges$pscore[i])
    present <- scores[scores > 0]
    if (length(present) == 0) return("")
    paste(paste0(ev_labels[names(present)], ": ",
                 round(present, 3)), collapse = "<br>")
  }, "")

  edges$title <- paste0("<b>Score: ", round(edges$score, 3), "</b>",
                        ifelse(edges$evidence != "",
                               paste0("<br><br>", edges$evidence), ""))
  edges$width <- edges$score * 5

  list(edges = edges, api_error = NULL)
}

#' Extract significant DE genes for a contrast
#' @return data.frame with gene_symbol, diff, padj columns
.get_de_genes <- function(dep, contrast, lfc_threshold = 1, alpha = 0.05,
                          use_adjp = TRUE) {
  row_data <- as.data.frame(SummarizedExperiment::rowData(dep))
  diff_col <- paste0(contrast, "_diff")
  padj_col <- paste0(contrast, if (use_adjp) "_p.adj" else "_p.val")

  if ("Gene" %in% colnames(row_data)) {
    row_data$gene_symbol <- row_data$Gene
  } else if ("Genes" %in% colnames(row_data)) {
    row_data$gene_symbol <- row_data$Genes
  } else {
    row_data$gene_symbol <- gsub("[.].*", "", row_data$name)
  }

  row_data %>%
    dplyr::filter(!is.na(.data[[diff_col]]) & !is.na(.data[[padj_col]])) %>%
    dplyr::filter(abs(.data[[diff_col]]) > lfc_threshold,
                  .data[[padj_col]] < alpha) %>%
    dplyr::select(gene_symbol,
                  diff  = dplyr::all_of(diff_col),
                  padj  = dplyr::all_of(padj_col)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "")
}

# --------------------------------------------------------------------------- #
#  Visual helpers for PPI network plots
# --------------------------------------------------------------------------- #

#' Map fold-change to a red-white-blue color
#' @param fc     Numeric vector of log2 fold-changes
#' @param max_fc Fold-change at which color fully saturates (default 3)
#' @return Character vector of hex colors
.fc_to_color <- function(fc, max_fc = 3) {
  max_fc <- max(max_fc, 0.1)  # avoid division by zero
  vapply(fc, function(v) {
    if (is.na(v)) return("#DCDCDC")  # grey for NA
    v <- max(min(v, max_fc), -max_fc)
    if (v > 0) {
      t <- v / max_fc
      r <- as.integer(220 - t * (220 - 211))
      g <- as.integer(220 - t * (220 - 47))
      b <- as.integer(220 - t * (220 - 47))
    } else {
      t <- abs(v) / max_fc
      r <- as.integer(220 - t * (220 - 41))
      g <- as.integer(220 - t * (220 - 128))
      b <- as.integer(220 - t * (220 - 185))
    }
    sprintf("#%02x%02x%02x", r, g, b)
  }, "")
}

#' Map STRING combined score to edge color (light grey → dark teal)
#'
#' Normalizes scores to the range [score_min, 1] so the full color gradient
#' is used regardless of the cutoff.
#' @param score     Numeric vector of STRING scores
#' @param score_min Minimum score (user cutoff), used as the low end of the range
#' @return Character vector of hex colors
.score_to_edge_color <- function(score, score_min = 0) {
  range_width <- max(1 - score_min, 0.01)  # avoid division by zero
  vapply(score, function(s) {
    s <- max(min(s, 1), score_min)
    t <- (s - score_min) / range_width  # normalize to 0–1 within the range
    # Light grey (#d5d8dc) → dark teal (#2c3e50)
    r <- as.integer(213 - t * (213 - 44))
    g <- as.integer(216 - t * (216 - 62))
    b <- as.integer(220 - t * (220 - 80))
    sprintf("#%02x%02x%02x", r, g, b)
  }, "")
}

#' Build CSS gradient string from a vector of hex colors
.css_gradient <- function(colors) {
  n <- length(colors)
  stops <- paste(colors, paste0(round(seq(0, 100, length.out = n)), "%"),
                 collapse = ", ")
  paste0("linear-gradient(to right, ", stops, ")")
}

#' Add custom select/deselect events to a visNetwork object
#'
#' Replaces the default highlightNearest behaviour:
#' - Click node → highlight only the node, its direct edges, and neighbor nodes
#' - Click edge → highlight the edge and its two endpoint nodes
#' - Click background → restore everything
#' @param vis A visNetwork object
#' @return The same visNetwork with events attached
.vis_select_events <- function(vis) {
  # JS: dim all nodes/edges to faded, keep selected ones at original color
  dim_color <- "rgba(200,200,200,0.15)"
  dim_font  <- "rgba(160,160,160,0.3)"

  # Helper JS functions shared across events (injected as IIFE)
  js_helpers <- paste0(
    "var _dim = '", dim_color, "';",
    "var _dimFont = '", dim_font, "';",
    # Save originals once
    "function _saveOrig(ctx) {",
    "  if (ctx._origSaved) return;",
    "  ctx._origSaved = true;",
    "  ctx._origNodes = {};",
    "  ctx._origEdges = {};",
    "  var nds = ctx.body.data.nodes.get();",
    "  nds.forEach(function(n) {",
    "    ctx._origNodes[n.id] = {",
    "      color: n.color, font: n.font || {}",
    "    };",
    "  });",
    "  var eds = ctx.body.data.edges.get();",
    "  eds.forEach(function(e) {",
    "    ctx._origEdges[e.id] = {",
    "      color: e.color, width: e.width",
    "    };",
    "  });",
    "}",
    # Dim everything
    "function _dimAll(ctx) {",
    "  _saveOrig(ctx);",
    "  var nUpd = ctx.body.data.nodes.get().map(function(n) {",
    "    return {id: n.id, color: _dim,",
    "      font: {color: _dimFont, strokeColor: _dimFont, strokeWidth: 0}};",
    "  });",
    "  var eUpd = ctx.body.data.edges.get().map(function(e) {",
    "    return {id: e.id, color: {color: _dim}, width: 0.5};",
    "  });",
    "  ctx.body.data.nodes.update(nUpd);",
    "  ctx.body.data.edges.update(eUpd);",
    "}",
    # Restore selected nodes & edges to originals
    "function _highlight(ctx, nodeIds, edgeIds) {",
    "  var nSet = {};",
    "  nodeIds.forEach(function(id){ nSet[id] = true; });",
    "  var eSet = {};",
    "  edgeIds.forEach(function(id){ eSet[id] = true; });",
    "  var nUpd = nodeIds.map(function(id) {",
    "    var o = ctx._origNodes[id] || {};",
    "    return {id: id, color: o.color,",
    "      font: {color: o.font.color || '#222222',",
    "        strokeColor: o.font.strokeColor || '#ffffff',",
    "        strokeWidth: o.font.strokeWidth != null ? o.font.strokeWidth : 3}",
    "    };",
    "  });",
    "  var eUpd = edgeIds.map(function(id) {",
    "    var o = ctx._origEdges[id] || {};",
    "    return {id: id, color: o.color, width: o.width};",
    "  });",
    "  ctx.body.data.nodes.update(nUpd);",
    "  ctx.body.data.edges.update(eUpd);",
    "}",
    # Restore all
    "function _restoreAll(ctx) {",
    "  if (!ctx._origSaved) return;",
    "  var nUpd = Object.keys(ctx._origNodes).map(function(id) {",
    "    var o = ctx._origNodes[id];",
    "    return {id: id, color: o.color,",
    "      font: {color: o.font.color || '#222222',",
    "        strokeColor: o.font.strokeColor || '#ffffff',",
    "        strokeWidth: o.font.strokeWidth != null ? o.font.strokeWidth : 3}",
    "    };",
    "  });",
    "  var eUpd = Object.keys(ctx._origEdges).map(function(id) {",
    "    var o = ctx._origEdges[id];",
    "    return {id: id, color: o.color, width: o.width};",
    "  });",
    "  ctx.body.data.nodes.update(nUpd);",
    "  ctx.body.data.edges.update(eUpd);",
    "}"
  )

  vis %>%
    visNetwork::visInteraction(
      hover = TRUE, tooltipDelay = 100,
      selectConnectedEdges = FALSE
    ) %>%
    visNetwork::visEvents(
      selectNode = paste0(
        "function(params) {",
        js_helpers,
        "  var nid = params.nodes[0];",
        "  var allEdges = this.body.data.edges.get();",
        "  var connEdges = allEdges.filter(function(e){",
        "    return e.from === nid || e.to === nid;",
        "  });",
        "  var eids = connEdges.map(function(e){ return e.id; });",
        "  var nids = connEdges.map(function(e){",
        "    return e.from === nid ? e.to : e.from;",
        "  });",
        "  nids.push(nid);",
        "  _dimAll(this);",
        "  _highlight(this, nids, eids);",
        "}"
      ),
      selectEdge = paste0(
        "function(params) {",
        js_helpers,
        "  if (params.nodes.length > 0) return;",
        "  var eid = params.edges[0];",
        "  var edge = this.body.data.edges.get(eid);",
        "  if (!edge) return;",
        "  _dimAll(this);",
        "  _highlight(this, [edge.from, edge.to], [eid]);",
        "}"
      ),
      deselectNode = paste0(
        "function(params) {",
        js_helpers,
        "  _restoreAll(this);",
        "}"
      ),
      deselectEdge = paste0(
        "function(params) {",
        js_helpers,
        "  _restoreAll(this);",
        "}"
      )
    )
}

#' Generate HTML legend for single-bait PPI network
#'
#' @param max_fc     Maximum |FC| for color scale (derived from data)
#' @param score_min  Minimum STRING score shown (0–1, same as user cutoff)
#' @param p_min      Minimum -log10(p.adj) in data (smallest dot)
#' @param p_max      Maximum -log10(p.adj) in data (largest dot)
#' @return shiny tagList
ppi_legend_html <- function(max_fc = 3, score_min = 0.4,
                            p_min = 1, p_max = 10) {
  max_fc <- max(round(max_fc, 1), 0.5)

  fc_vals <- seq(-max_fc, max_fc, length.out = 11)
  fc_colors <- .fc_to_color(fc_vals, max_fc = max_fc)
  fc_grad <- .css_gradient(fc_colors)

  sc_vals <- seq(score_min, 1, length.out = 11)
  sc_colors <- .score_to_edge_color(sc_vals, score_min = score_min)
  sc_grad <- .css_gradient(sc_colors)

  bar_style   <- "height:12px; border-radius:3px; border:1px solid #ccc; "
  label_style <- "font-size:13px; font-weight:600; color:#444; margin-bottom:2px;"
  tick_style  <- "font-size:11px; color:#666;"

  p_min_label <- round(p_min, 1)
  p_max_label <- round(p_max, 1)

  htmltools::tagList(
    htmltools::tags$div(
      style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-start; padding:6px 0;",

      # -- Node color (FC) --
      htmltools::tags$div(
        style = "min-width:160px; flex:1;",
        htmltools::tags$div(style = label_style, "Node color: log2 FC"),
        htmltools::tags$div(style = paste0(bar_style, "background:", fc_grad, ";")),
        htmltools::tags$div(
          style = "display:flex; justify-content:space-between;",
          htmltools::tags$span(style = tick_style, paste0("-", max_fc)),
          htmltools::tags$span(style = tick_style, "0"),
          htmltools::tags$span(style = tick_style, paste0("+", max_fc))
        )
      ),

      # -- Edge color/width (STRING score) --
      htmltools::tags$div(
        style = "min-width:160px; flex:1;",
        htmltools::tags$div(style = label_style, "Edge: STRING score"),
        htmltools::tags$div(style = paste0(bar_style, "background:", sc_grad, ";")),
        htmltools::tags$div(
          style = "display:flex; justify-content:space-between;",
          htmltools::tags$span(style = tick_style, score_min),
          htmltools::tags$span(style = tick_style, "1.0")
        )
      ),

      # -- Node size with labeled p-values --
      htmltools::tags$div(
        style = "min-width:160px; flex:0 0 auto;",
        htmltools::tags$div(style = label_style, "Node size: -log10(p.adj)"),
        htmltools::tags$div(
          style = "display:flex; align-items:center; gap:6px; margin-top:2px;",
          htmltools::tags$div(style = "text-align:center;",
            htmltools::tags$div(
              style = "width:10px; height:10px; border-radius:50%; background:#aab7b8; border:1px solid #888; margin:0 auto;"
            ),
            htmltools::tags$div(style = tick_style, p_min_label)
          ),
          htmltools::tags$span(style = tick_style, "\u2014"),
          htmltools::tags$div(style = "text-align:center;",
            htmltools::tags$div(
              style = "width:24px; height:24px; border-radius:50%; background:#aab7b8; border:1px solid #888; margin:0 auto;"
            ),
            htmltools::tags$div(style = tick_style, p_max_label)
          )
        )
      ),

      # -- Bait shape --
      htmltools::tags$div(
        style = "min-width:80px; flex:0 0 auto; text-align:center;",
        htmltools::tags$div(style = label_style, "Bait"),
        htmltools::tags$div(
          style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-top:4px;",
          htmltools::tags$div(
            style = paste0("width:18px; height:18px; background:#aab7b8; ",
                           "transform:rotate(45deg); border:2px solid #666;")
          ),
          htmltools::tags$span(style = tick_style, "diamond")
        )
      )
    )
  )
}

#' Generate HTML legend for multi-bait PPI network
#'
#' @param bait_names  Character vector of bait names
#' @param bait_colors Character vector of hex colors for each bait
#' @param score_min   Minimum STRING score shown
#' @param p_min       Minimum -log10(p.adj) in data
#' @param p_max       Maximum -log10(p.adj) in data
#' @return shiny tagList
ppi_multi_legend_html <- function(bait_names, bait_colors,
                                  score_min = 0.4, p_min = 1, p_max = 10) {
  sc_vals <- seq(score_min, 1, length.out = 11)
  sc_colors <- .score_to_edge_color(sc_vals, score_min = score_min)
  sc_grad <- .css_gradient(sc_colors)

  bar_style   <- "height:12px; border-radius:3px; border:1px solid #ccc; "
  label_style <- "font-size:13px; font-weight:600; color:#444; margin-bottom:2px;"
  tick_style  <- "font-size:11px; color:#666;"

  p_min_label <- round(p_min, 1)
  p_max_label <- round(p_max, 1)

  # Bait color swatches
  bait_items <- lapply(seq_along(bait_names), function(i) {
    htmltools::tags$div(
      style = "display:inline-flex; align-items:center; gap:4px; margin-right:12px;",
      htmltools::tags$div(
        style = paste0("width:12px; height:12px; border-radius:50%; ",
                       "background:", bait_colors[i], "; border:1px solid #888;")
      ),
      htmltools::tags$span(style = tick_style, bait_names[i])
    )
  })

  htmltools::tagList(
    htmltools::tags$div(
      style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-start; padding:6px 0;",

      # -- Bait colors --
      htmltools::tags$div(
        style = "min-width:180px; flex:1;",
        htmltools::tags$div(style = label_style, "Node color: bait specificity"),
        htmltools::tags$div(
          style = "display:flex; flex-wrap:wrap; align-items:center; gap:2px;",
          bait_items,
          htmltools::tags$div(
            style = "display:inline-flex; align-items:center; gap:4px; margin-right:12px;",
            htmltools::tags$div(
              style = "width:12px; height:12px; background:#f1c40f; border:1px solid #c9a800; border-radius:2px;"
            ),
            htmltools::tags$span(style = tick_style, "Shared")
          ),
          htmltools::tags$div(
            style = "display:inline-flex; align-items:center; gap:4px;",
            htmltools::tags$div(
              style = "width:12px; height:12px; border-radius:50%; background:#d5d8dc; border:1px solid #aaa;"
            ),
            htmltools::tags$span(style = tick_style, "STRING only")
          )
        )
      ),

      # -- Edge score --
      htmltools::tags$div(
        style = "min-width:150px; flex:0.7;",
        htmltools::tags$div(style = label_style, "Edge: STRING score"),
        htmltools::tags$div(style = paste0(bar_style, "background:", sc_grad, ";")),
        htmltools::tags$div(
          style = "display:flex; justify-content:space-between;",
          htmltools::tags$span(style = tick_style, score_min),
          htmltools::tags$span(style = tick_style, "1.0")
        )
      ),

      # -- Node size with labeled p-values --
      htmltools::tags$div(
        style = "min-width:160px; flex:0 0 auto;",
        htmltools::tags$div(style = label_style, "Node size: -log10(p.adj)"),
        htmltools::tags$div(
          style = "display:flex; align-items:center; gap:6px; margin-top:2px;",
          htmltools::tags$div(style = "text-align:center;",
            htmltools::tags$div(
              style = "width:10px; height:10px; border-radius:50%; background:#aab7b8; border:1px solid #888; margin:0 auto;"
            ),
            htmltools::tags$div(style = tick_style, p_min_label)
          ),
          htmltools::tags$span(style = tick_style, "\u2014"),
          htmltools::tags$div(style = "text-align:center;",
            htmltools::tags$div(
              style = "width:24px; height:24px; border-radius:50%; background:#aab7b8; border:1px solid #888; margin:0 auto;"
            ),
            htmltools::tags$div(style = tick_style, p_max_label)
          )
        )
      ),

      # -- Shapes --
      htmltools::tags$div(
        style = "min-width:100px; flex:0 0 auto; text-align:center;",
        htmltools::tags$div(style = label_style, "Shape"),
        htmltools::tags$div(
          style = "display:flex; align-items:center; gap:5px; margin-top:4px;",
          htmltools::tags$div(
            style = "width:14px; height:14px; background:#aab7b8; transform:rotate(45deg); border:2px solid #666;"
          ),
          htmltools::tags$span(style = tick_style, "Bait"),
          htmltools::tags$div(
            style = "width:12px; height:12px; background:#f1c40f; border:1px solid #c9a800; border-radius:2px; margin-left:6px;"
          ),
          htmltools::tags$span(style = tick_style, "Shared")
        )
      )
    )
  )
}

#' PPI Network via STRING REST API (bait-centric, single contrast)
#'
#' Fetches protein-protein interactions for significant DE genes from the
#' STRING database REST API and renders an interactive bait-centric visNetwork
#' graph.  The bait protein is placed at the center; prey are arranged around
#' it sized by interaction confidence (-log10 p.adj).
#'
#' @param dep            SummarizedExperiment object after test_diff()
#' @param contrast       Name of the contrast (e.g. "BaitA_vs_Control")
#' @param bait           Bait gene symbol (placed at center; optional — inferred
#'                       from contrast name if NULL)
#' @param lfc_threshold  Absolute log2 fold-change cutoff
#' @param alpha          Adjusted p-value cutoff
#' @param string_score   Minimum combined STRING interaction score (0–1)
#' @param species_id     NCBI taxonomy ID (9606 = human, 10090 = mouse)
#' @return list with $network (visNetwork object) and $meta (list of
#'         max_fc, score_min, p_min, p_max for legend rendering)
plot_ppi_network <- function(dep, contrast, bait = NULL,
                             lfc_threshold = 1, alpha = 0.05,
                             string_score = 0.4, species_id = 9606,
                             use_adjp = TRUE,
                             bait_connected_only = FALSE) {

  de_genes <- .get_de_genes(dep, contrast, lfc_threshold, alpha, use_adjp)

  if (nrow(de_genes) == 0) {
    nodes <- data.frame(id = 1, label = "No significant DE genes found",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(list(
      network = visNetwork::visNetwork(nodes, data.frame()),
      meta    = list(max_fc = 3, score_min = string_score,
                     p_min = 1, p_max = 10)
    ))
  }

  # Infer bait from contrast name — try both sides of _vs_
  if (is.null(bait) || bait == "") {
    left_side  <- gsub("_vs_.*", "", contrast)
    right_side <- gsub(".*_vs_", "", contrast)
    bait <- left_side
  } else {
    left_side  <- bait
    right_side <- NULL
  }

  # Try to match bait to a gene symbol (exact, case-insensitive, or partial)
  # Check both sides of the contrast name
  .match_bait <- function(candidate, gene_syms) {
    if (candidate %in% gene_syms) return(candidate)
    ci <- gene_syms[toupper(gene_syms) == toupper(candidate)]
    if (length(ci) > 0) return(ci[1])
    uc <- toupper(candidate)
    ugs <- toupper(gene_syms)
    # candidate is substring of gene, or gene is substring of candidate
    partial <- gene_syms[grepl(uc, ugs, fixed = TRUE) |
                         vapply(ugs, function(g) grepl(g, uc, fixed = TRUE), FALSE)]
    if (length(partial) > 0) return(partial[1])
    NULL
  }

  gene_syms <- de_genes$gene_symbol
  bait_gene <- .match_bait(left_side, gene_syms)
  if (is.null(bait_gene) && !is.null(right_side))
    bait_gene <- .match_bait(right_side, gene_syms)
  if (is.null(bait_gene)) {
    message("[PPI] Could not match bait from contrast '", contrast, "' to any DE gene")
    bait_gene <- left_side  # fallback — no diamond will be shown
  }

  # --- STRING REST API ---
  string_res <- .query_string_api(de_genes$gene_symbol, species_id, string_score)
  edges <- string_res$edges

  # --- Compute data-driven ranges ---
  max_fc    <- ceiling(max(abs(de_genes$diff), na.rm = TRUE))
  score_min <- string_score
  conf      <- pmin(-log10(pmax(de_genes$padj, 1e-300)), 20)
  p_min_val <- min(conf, na.rm = TRUE)
  p_max_val <- max(conf, na.rm = TRUE)

  # --- Build nodes ---
  network_genes <- if (nrow(edges) > 0) unique(c(edges$from, edges$to))
                   else character(0)

  # Confidence-based sizing
  conf_norm   <- conf / p_max_val
  conf_scaled <- 8 + conf_norm * 25  # range 8–33

  is_bait <- de_genes$gene_symbol == bait_gene

  # All nodes colored by FC (including bait — bait distinguished by shape)
  node_colors <- .fc_to_color(de_genes$diff, max_fc = max_fc)

  nodes <- data.frame(
    id    = de_genes$gene_symbol,
    label = de_genes$gene_symbol,
    title = paste0("<b>", de_genes$gene_symbol, "</b>",
                   ifelse(is_bait, " (bait)", ""),
                   "<br>log2FC: ", round(de_genes$diff, 3),
                   "<br>p.adj: ",  signif(de_genes$padj, 3)),
    color = node_colors,
    value = ifelse(is_bait, max(conf_scaled) * 1.3, conf_scaled),
    shape = ifelse(is_bait, "diamond", "dot"),
    borderWidth = ifelse(is_bait, 3,
                         ifelse(de_genes$gene_symbol %in% network_genes, 2, 1)),
    font.size = ifelse(is_bait, 18, 14),
    font.color = "#222222",
    font.strokeWidth = 3,
    font.strokeColor = "#ffffff",
    fixed = ifelse(is_bait, TRUE, FALSE),
    x     = ifelse(is_bait, 0, NA_real_),
    y     = ifelse(is_bait, 0, NA_real_),
    stringsAsFactors = FALSE
  )

  # Add STRING-only neighbors not in our DE list
  extra <- setdiff(network_genes, nodes$id)
  if (length(extra) > 0) {
    extra_nodes <- data.frame(
      id = extra, label = extra, title = extra,
      color = "#d5d8dc", value = 8, shape = "dot",
      borderWidth = 1, font.size = 12,
      font.color = "#222222",
      font.strokeWidth = 3,
      font.strokeColor = "#ffffff",
      fixed = FALSE, x = NA_real_, y = NA_real_,
      stringsAsFactors = FALSE
    )
    nodes <- dplyr::bind_rows(nodes, extra_nodes)
  }

  # Style edges: color by score (light grey → dark blue-grey)
  if (nrow(edges) > 0) {
    edges$color <- .score_to_edge_color(edges$score, score_min = score_min)
    range_w <- max(1 - score_min, 0.01)
    edges$width <- 1 + ((edges$score - score_min) / range_w) * 5   # 1–6
  }

  # Filter to bait-connected subnetwork if requested
  if (bait_connected_only && nrow(edges) > 0) {
    bait_id <- if (bait_gene %in% nodes$id) bait_gene else nodes$id[1]
    adj <- list()
    for (nid in nodes$id) adj[[nid]] <- character(0)
    for (r in seq_len(nrow(edges))) {
      f <- edges$from[r]; t <- edges$to[r]
      if (f %in% names(adj)) adj[[f]] <- c(adj[[f]], t)
      if (t %in% names(adj)) adj[[t]] <- c(adj[[t]], f)
    }
    queue <- bait_id; visited <- character(0)
    while (length(queue) > 0) {
      cur <- queue[1]; queue <- queue[-1]
      if (cur %in% visited) next
      visited <- c(visited, cur)
      queue <- c(queue, setdiff(adj[[cur]], visited))
    }
    nodes <- nodes[nodes$id %in% visited, ]
    edges <- edges[edges$from %in% visited & edges$to %in% visited, ]
  }

  vis <- visNetwork::visNetwork(
    nodes, edges,
    main = list(text = paste0("PPI Network: ", contrast),
                style = "font-size:14px; font-weight:bold;")
  ) %>%
    visNetwork::visOptions(nodesIdSelection = FALSE) %>%
    .vis_select_events() %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           stabilization = list(iterations = 200),
                           forceAtlas2Based = list(
                             gravitationalConstant = -30,
                             centralGravity = 0.01,
                             springLength = 100
                           )) %>%
    visNetwork::visEdges(smooth = list(type = "continuous"))

  list(
    network = vis,
    meta    = list(max_fc = max_fc, score_min = score_min,
                   p_min = round(p_min_val, 1), p_max = round(p_max_val, 1))
  )
}


#' Multi-bait PPI Network
#'
#' Fetches STRING interactions for significant DE genes across multiple
#' contrasts (baits) and renders a combined network.  Each bait is shown
#' as a diamond node; prey are colored by which bait(s) they associate with.
#' Shared prey (found in 2+ baits) are highlighted.
#'
#' @param dep            SummarizedExperiment object after test_diff()
#' @param contrasts      Character vector of contrast names
#' @param lfc_threshold  Absolute log2 fold-change cutoff
#' @param alpha          Adjusted p-value cutoff
#' @param string_score   Minimum combined STRING interaction score (0–1)
#' @param species_id     NCBI taxonomy ID
#' @return list with $network (visNetwork object) and $meta
plot_ppi_network_multi <- function(dep, contrasts, lfc_threshold = 1,
                                   alpha = 0.05, string_score = 0.4,
                                   species_id = 9606, use_adjp = TRUE,
                                   bait_connected_only = FALSE) {

  score_min <- string_score
  default_meta <- list(score_min = score_min, p_min = 1, p_max = 10)

  if (length(contrasts) == 0) {
    nodes <- data.frame(id = 1, label = "No contrasts selected",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(list(network = visNetwork::visNetwork(nodes, data.frame()),
                meta = default_meta))
  }

  # Bait palette (up to 8 distinct baits)
  bait_palette <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12",
                    "#9b59b6", "#1abc9c", "#e67e22", "#34495e")

  # Collect DE genes per contrast
  bait_names <- gsub("_vs_.*", "", contrasts)
  all_de <- list()
  for (i in seq_along(contrasts)) {
    de <- .get_de_genes(dep, contrasts[i], lfc_threshold, alpha, use_adjp)
    if (nrow(de) > 0) {
      de$bait     <- bait_names[i]
      de$contrast <- contrasts[i]
      all_de[[i]] <- de
    }
  }

  if (length(all_de) == 0) {
    nodes <- data.frame(id = 1, label = "No significant DE genes found",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(list(network = visNetwork::visNetwork(nodes, data.frame()),
                meta = default_meta))
  }

  combined_de <- dplyr::bind_rows(all_de)
  unique_genes <- unique(combined_de$gene_symbol)

  # Query STRING for all genes at once
  string_res <- .query_string_api(unique_genes, species_id, string_score)
  edges <- string_res$edges

  # Build per-gene summary: which baits, best FC/p across contrasts
  gene_summary <- combined_de %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(
      baits    = paste(unique(bait), collapse = ", "),
      n_baits  = dplyr::n_distinct(bait),
      best_lfc = diff[which.min(padj)],
      best_p   = min(padj),
      .groups  = "drop"
    )

  bait_color_map <- setNames(
    bait_palette[seq_along(bait_names)], bait_names)

  # Assign node color: single-bait → bait color, multi-bait → gold highlight
  gene_summary$color <- vapply(seq_len(nrow(gene_summary)), function(i) {
    if (gene_summary$n_baits[i] > 1) return("#f1c40f")  # shared = gold
    bait_color_map[strsplit(gene_summary$baits[i], ", ")[[1]][1]]
  }, "")

  # Confidence sizing
  conf <- pmin(-log10(pmax(gene_summary$best_p, 1e-300)), 20)
  p_min_val <- min(conf, na.rm = TRUE)
  p_max_val <- max(conf, na.rm = TRUE)
  conf_norm <- conf / p_max_val
  conf_scaled <- 8 + conf_norm * 25

  is_bait <- gene_summary$gene_symbol %in% bait_names

  nodes <- data.frame(
    id    = gene_summary$gene_symbol,
    label = gene_summary$gene_symbol,
    title = paste0("<b>", gene_summary$gene_symbol, "</b>",
                   ifelse(is_bait, " (bait)", ""),
                   "<br>Baits: ", gene_summary$baits,
                   "<br>Best log2FC: ", round(gene_summary$best_lfc, 3),
                   "<br>Best p.adj: ", signif(gene_summary$best_p, 3)),
    color = ifelse(is_bait, "#f39c12", gene_summary$color),
    value = ifelse(is_bait, max(conf_scaled) * 1.3, conf_scaled),
    shape = ifelse(is_bait, "diamond",
                   ifelse(gene_summary$n_baits > 1, "square", "dot")),
    borderWidth = ifelse(is_bait, 3,
                         ifelse(gene_summary$n_baits > 1, 3, 1)),
    font.size = ifelse(is_bait, 18, 14),
    font.color = "#222222",
    font.strokeWidth = 3,
    font.strokeColor = "#ffffff",
    stringsAsFactors = FALSE
  )

  # Add STRING-only neighbors
  network_genes <- if (nrow(edges) > 0) unique(c(edges$from, edges$to))
                   else character(0)
  extra <- setdiff(network_genes, nodes$id)
  if (length(extra) > 0) {
    extra_nodes <- data.frame(
      id = extra, label = extra, title = extra,
      color = "#d5d8dc", value = 8, shape = "dot",
      borderWidth = 1, font.size = 12,
      font.color = "#222222",
      font.strokeWidth = 3,
      font.strokeColor = "#ffffff",
      stringsAsFactors = FALSE
    )
    nodes <- dplyr::bind_rows(nodes, extra_nodes)
  }

  # Style edges
  if (nrow(edges) > 0) {
    edges$color <- .score_to_edge_color(edges$score, score_min = score_min)
    range_w <- max(1 - score_min, 0.01)
    edges$width <- 1 + ((edges$score - score_min) / range_w) * 5
  }

  # Filter to bait-connected subnetwork if requested
  if (bait_connected_only && nrow(edges) > 0) {
    # Build adjacency list and find connected components via BFS
    all_node_ids <- nodes$id
    adj <- list()
    for (nid in all_node_ids) adj[[nid]] <- character(0)
    for (r in seq_len(nrow(edges))) {
      f <- edges$from[r]; t <- edges$to[r]
      if (f %in% names(adj)) adj[[f]] <- c(adj[[f]], t)
      if (t %in% names(adj)) adj[[t]] <- c(adj[[t]], f)
    }
    # BFS from each bait to find all reachable nodes
    bait_ids <- intersect(bait_names, all_node_ids)
    # Also match baits by the same logic as single-bait
    if (length(bait_ids) == 0) bait_ids <- all_node_ids[1:min(2, length(all_node_ids))]
    reachable <- character(0)
    for (bid in bait_ids) {
      queue <- bid; visited <- character(0)
      while (length(queue) > 0) {
        cur <- queue[1]; queue <- queue[-1]
        if (cur %in% visited) next
        visited <- c(visited, cur)
        neighbors <- adj[[cur]]
        queue <- c(queue, setdiff(neighbors, visited))
      }
      reachable <- base::union(reachable, visited)
    }
    # Keep only reachable nodes and their edges
    nodes <- nodes[nodes$id %in% reachable, ]
    edges <- edges[edges$from %in% reachable & edges$to %in% reachable, ]
  }

  vis <- visNetwork::visNetwork(
    nodes, edges,
    main = list(text = paste0("Multi-bait PPI (",
                              paste(bait_names, collapse = ", "), ")"),
                style = "font-size:14px; font-weight:bold;")
  ) %>%
    visNetwork::visOptions(nodesIdSelection = FALSE) %>%
    .vis_select_events() %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           stabilization = list(iterations = 250),
                           forceAtlas2Based = list(
                             gravitationalConstant = -40,
                             centralGravity = 0.005,
                             springLength = 120
                           )) %>%
    visNetwork::visEdges(smooth = list(type = "continuous"))

  list(
    network = vis,
    meta    = list(score_min = score_min,
                   p_min = round(p_min_val, 1), p_max = round(p_max_val, 1))
  )
}


# --------------------------------------------------------------------------- #
#  Kinase-Substrate Network (from PTM-SEA results + PTMsigDB)
# --------------------------------------------------------------------------- #

#' Build a kinase–substrate network from PTM-SEA results
#'
#' Takes significant kinase sets from PTM-SEA, retrieves their substrate
#' members from PTMsigDB, matches substrates to the user's DE site data,
#' and renders a visNetwork with kinase hub nodes and substrate leaf nodes.
#'
#' @param ptmsea_results  Data.frame returned by run_ptmsea() (columns: set,
#'                        NES, p_value, adj_p_value, set_size, ES)
#' @param dep             SummarizedExperiment (site-level, after test_diff())
#' @param contrast        Contrast name
#' @param species         "human", "mouse", or "rat"
#' @param mod_type        Modification type ("p", "ac", "ub", etc.)
#' @param nes_cutoff      Minimum |NES| to include a kinase (default 1.5)
#' @param p_cutoff        Adjusted p-value cutoff for kinases (default 0.05)
#' @param max_kinases     Maximum number of kinases to display (default 15)
#' @return list with $network (visNetwork) and $meta (list for legend)
plot_kinase_substrate_network <- function(ptmsea_results, dep, contrast,
                                          species = "human", mod_type = "p",
                                          nes_cutoff = 1.5, p_cutoff = 0.05,
                                          use_adjp = TRUE,
                                          site_lfc_cutoff = 0,
                                          site_p_cutoff = 1,
                                          site_use_adjp = TRUE,
                                          hide_empty_kinases = TRUE) {

  # --- 1. Filter to current contrast and significant kinase sets ---
  if ("contrast" %in% colnames(ptmsea_results))
    ptmsea_results <- ptmsea_results[ptmsea_results$contrast == contrast, ]
  kinase_res <- ptmsea_results[grepl("^KINASE-", ptmsea_results$set), ]

  if (nrow(kinase_res) == 0) {
    nodes <- data.frame(id = 1, label = "No KINASE sets in PTM-SEA results",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(list(
      network = visNetwork::visNetwork(nodes, data.frame()),
      meta    = list(max_nes = 3, has_substrates = FALSE)
    ))
  }

  ks_p_col <- if (use_adjp) "adj_p_value" else "p_value"
  message(sprintf("[KS-Net] %d KINASE sets found. Filtering: |NES| >= %.2f, %s <= %.4f",
                  nrow(kinase_res), nes_cutoff, ks_p_col, p_cutoff))

  # Log kinases that fail cutoffs for debugging
  failed <- kinase_res[is.na(kinase_res[[ks_p_col]]) |
                        kinase_res[[ks_p_col]] > p_cutoff |
                        abs(kinase_res$NES) < nes_cutoff, ]
  if (nrow(failed) > 0) {
    top_failed <- head(failed[order(-abs(failed$NES)), ], 5)
    for (r in seq_len(nrow(top_failed))) {
      message(sprintf("[KS-Net]   Excluded: %s (NES=%.2f, %s=%.4f)",
                      top_failed$set[r], top_failed$NES[r],
                      ks_p_col, top_failed[[ks_p_col]][r]))
    }
  }

  kinase_res <- kinase_res[!is.na(kinase_res[[ks_p_col]]) &
                           kinase_res[[ks_p_col]] <= p_cutoff &
                           abs(kinase_res$NES) >= nes_cutoff, ]

  message(sprintf("[KS-Net] %d kinases pass cutoffs", nrow(kinase_res)))

  if (nrow(kinase_res) == 0) {
    nodes <- data.frame(id = 1,
                        label = "No kinases pass significance cutoffs",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(list(
      network = visNetwork::visNetwork(nodes, data.frame()),
      meta    = list(max_nes = 3, has_substrates = FALSE)
    ))
  }

  # Deduplicate PSP vs iKiP kinase sets: when the same kinase appears in both
  # sources, keep the PSP (literature-curated) set; use iKiP only when PSP
  # is absent.  Extract canonical gene symbol from the set name to group them.
  kinase_res$source <- ifelse(grepl("^KINASE-PSP_", kinase_res$set), "PSP", "iKiP")
  kinase_res$canon <- toupper(gsub("^KINASE-(PSP|iKiP)_", "", kinase_res$set))
  # Normalize PSP aliases: "Akt1/AKT1" → "AKT1" (take the part after "/")
  kinase_res$canon <- ifelse(
    grepl("/", kinase_res$canon),
    toupper(sub(".*/", "", kinase_res$canon)),
    kinase_res$canon
  )
  # Remove isoform suffixes like "_ISO2"
  kinase_res$canon <- gsub("_ISO\\d+$", "", kinase_res$canon)
  # Also normalize iKiP complex names: "CDK1-CCNB1" → "CDK1"
  kinase_res$canon <- gsub("-.*$", "", kinase_res$canon)
  # Dot separator in iKiP (e.g., "MAPK1.ERK2") → take first part
  kinase_res$canon <- gsub("\\..*$", "", kinase_res$canon)

  n_before_dedup <- nrow(kinase_res)
  # For each canonical kinase, keep PSP if available, otherwise iKiP
  kinase_res <- kinase_res[order(kinase_res$canon,
                                  ifelse(kinase_res$source == "PSP", 0, 1)), ]
  kinase_res <- kinase_res[!duplicated(kinase_res$canon), ]
  message(sprintf("[KS-Net] After PSP/iKiP dedup: %d → %d kinases",
                  n_before_dedup, nrow(kinase_res)))

  # Sort by |NES|
  kinase_res <- kinase_res[order(-abs(kinase_res$NES)), ]

  # --- 2. Load PTMsigDB and prepare DE data ---
  gmt_path <- get_ptmsigdb_path(species)

  # Parse the GMT to extract flanking sequence members per set
  gmt_list <- read_gmt(gmt_path)

  rd <- as.data.frame(SummarizedExperiment::rowData(dep))
  diff_col <- paste0(contrast, "_diff")
  padj_col <- paste0(contrast, if (site_use_adjp) "_p.adj" else "_p.val")

  # Build flanking key for matching (same format as PTM-SEA input)
  mod_suffix <- paste0("-", mod_type)
  rd$flank_key <- paste0(toupper(rd$SequenceWindow), mod_suffix)

  # Build site label for display (e.g., "AKT1_S473")
  gene_sym <- if ("Gene" %in% colnames(rd)) rd$Gene else gsub("[.].*", "", rd$name)
  site_pos <- regmatches(rd$Index, regexpr("[STY][0-9]+", rd$Index))
  rd$site_label <- ifelse(nchar(site_pos) > 0 & !is.na(gene_sym),
                          paste0(gene_sym, "_", site_pos),
                          rd$Index)

  # Build lookup: flank_key → row index (for fast matching)
  flank_lookup <- setNames(seq_len(nrow(rd)), rd$flank_key)
  # Keep first occurrence for duplicates
  flank_lookup <- flank_lookup[!duplicated(names(flank_lookup))]

  message(sprintf("[KS-Net] DE sites with flanking key: %d / %d",
                  sum(!is.na(rd$SequenceWindow) & nchar(trimws(rd$SequenceWindow)) > 0),
                  nrow(rd)))

  # --- 3. Build kinase and substrate nodes + edges ---
  max_nes <- max(abs(kinase_res$NES), na.rm = TRUE)

  kinase_label <- gsub("^KINASE-(PSP|iKiP)_", "", kinase_res$set)
  kinase_nodes <- data.frame(
    id    = kinase_res$set,
    label = kinase_label,
    title = paste0("<b>", kinase_label, "</b>",
                   "<br>Set: ", kinase_res$set,
                   "<br>NES: ", round(kinase_res$NES, 3),
                   "<br>adj.p: ", signif(kinase_res$adj_p_value, 3),
                   "<br>Set size: ", kinase_res$set_size),
    color = .nes_to_color(kinase_res$NES, max_nes = max_nes),
    value = 30 + abs(kinase_res$NES) / max_nes * 20,
    shape = "diamond",
    borderWidth = 3,
    font.size = 16,
    font.color = "#222222",
    font.strokeWidth = 3,
    font.strokeColor = "#ffffff",
    type  = "kinase",
    stringsAsFactors = FALSE
  )

  all_substrate_nodes <- list()
  all_edges <- list()

  for (i in seq_len(nrow(kinase_res))) {
    set_name <- kinase_res$set[i]
    members <- gmt_list[[set_name]]
    n_gmt_members <- length(members)
    if (is.null(members) || n_gmt_members == 0) {
      message(sprintf("[KS-Net]   %s: 0 members in GMT", set_name))
      next
    }

    # Filter by mod type
    mod_pattern <- paste0("-", mod_type, "(;[ud])?$")
    members <- members[grepl(mod_pattern, members)]
    if (length(members) == 0) {
      message(sprintf("[KS-Net]   %s: %d GMT members, 0 match mod type '%s'",
                      set_name, n_gmt_members, mod_type))
      next
    }

    n_matched <- 0L; n_na <- 0L; n_lfc_fail <- 0L; n_p_fail <- 0L; n_pass <- 0L

    for (j in seq_along(members)) {
      member <- members[j]

      # Match by flanking sequence: strip ";u"/";d" suffix, uppercase the
      # flanking portion only (keep mod-type suffix like "-p" lowercase)
      flank_key <- gsub(";[ud]$", "", member)
      dash_pos <- regexpr("-", flank_key)
      if (dash_pos > 0) {
        flank_upper <- paste0(toupper(substr(flank_key, 1, dash_pos - 1)),
                              substr(flank_key, dash_pos, nchar(flank_key)))
      } else {
        flank_upper <- toupper(flank_key)
      }
      rd_idx <- flank_lookup[flank_upper]

      if (is.na(rd_idx)) next
      n_matched <- n_matched + 1L
      rd_idx <- as.integer(rd_idx)

      # Build substrate node
      sid     <- rd$flank_key[rd_idx]
      fc_val  <- rd[[diff_col]][rd_idx]
      p_val   <- rd[[padj_col]][rd_idx]
      slabel  <- rd$site_label[rd_idx]

      # Skip substrates with NA values or not passing cutoffs
      if (is.na(fc_val) || is.na(p_val)) { n_na <- n_na + 1L; next }
      if (abs(fc_val) < site_lfc_cutoff) { n_lfc_fail <- n_lfc_fail + 1L; next }
      if (p_val > site_p_cutoff) { n_p_fail <- n_p_fail + 1L; next }
      n_pass <- n_pass + 1L

      if (!sid %in% names(all_substrate_nodes)) {
        all_substrate_nodes[[sid]] <- data.frame(
          id    = sid,
          label = slabel,
          title = paste0("<b>", slabel, "</b>",
                         "<br>log2FC: ", round(fc_val, 3),
                         "<br>p.adj: ", signif(p_val, 3)),
          color = .fc_to_color(fc_val, max_fc = 3),
          value = 10,
          shape = "dot",
          borderWidth = 1,
          font.size = 11,
          font.color = "#222222",
          font.strokeWidth = 2,
          font.strokeColor = "#ffffff",
          type  = "substrate",
          stringsAsFactors = FALSE
        )
      }

      all_edges[[length(all_edges) + 1]] <- data.frame(
        from   = set_name,
        to     = sid,
        arrows = "to",
        color  = "#888888",
        title  = paste0(slabel, " \u2014 kinase substrate"),
        width  = 1.5,
        dashes = FALSE,
        stringsAsFactors = FALSE
      )
    }
    message(sprintf("[KS-Net]   %s: %d members, %d matched DE, %d NA, %d fail LFC, %d fail p, %d pass",
                    set_name, length(members), n_matched, n_na, n_lfc_fail, n_p_fail, n_pass))
  }

  substrate_df <- if (length(all_substrate_nodes) > 0)
    do.call(rbind, all_substrate_nodes) else data.frame()
  edges_df <- if (length(all_edges) > 0)
    do.call(rbind, all_edges) else data.frame()

  # Hide kinases that have no substrates passing the filters
  n_kinases_before <- nrow(kinase_nodes)
  if (hide_empty_kinases && nrow(edges_df) > 0) {
    kinases_with_subs <- unique(edges_df$from)
    kinase_nodes <- kinase_nodes[kinase_nodes$id %in% kinases_with_subs, ]
  } else if (hide_empty_kinases && nrow(edges_df) == 0) {
    kinase_nodes <- kinase_nodes[FALSE, ]
  }
  message(sprintf("[KS-Net] Kinases: %d passed enrichment cutoff, %d have substrates in data",
                  n_kinases_before, nrow(kinase_nodes)))

  if (nrow(substrate_df) == 0) {
    nodes <- dplyr::bind_rows(
      kinase_nodes,
      data.frame(id = "no_sub", label = "No substrate matches in DE data",
                 color = "#95a5a6", value = 10, shape = "dot",
                 borderWidth = 1, font.size = 12,
                 font.color = "#222222", font.strokeWidth = 2,
                 font.strokeColor = "#ffffff", type = "info",
                 stringsAsFactors = FALSE)
    )
    return(list(
      network = visNetwork::visNetwork(nodes, data.frame()) %>%
        .vis_select_events() %>%
        visNetwork::visLayout(randomSeed = 42),
      meta = list(max_nes = round(max_nes, 1), has_substrates = FALSE)
    ))
  }

  # Recalculate max_nes from remaining kinases (for legend)
  if (nrow(kinase_nodes) > 0) {
    remaining_nes <- kinase_res$NES[kinase_res$set %in% kinase_nodes$id]
    if (length(remaining_nes) > 0) max_nes <- max(abs(remaining_nes), na.rm = TRUE)
  }

  nodes <- dplyr::bind_rows(kinase_nodes, substrate_df)

  # --- 4. Render visNetwork ---
  vis <- visNetwork::visNetwork(
    nodes, edges_df,
    main = list(text = paste0("Kinase-Substrate Network: ", contrast),
                style = "font-size:14px; font-weight:bold;")
  ) %>%
    visNetwork::visOptions(nodesIdSelection = FALSE) %>%
    .vis_select_events() %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           stabilization = list(iterations = 200),
                           forceAtlas2Based = list(
                             gravitationalConstant = -40,
                             centralGravity = 0.005,
                             springLength = 120
                           )) %>%
    visNetwork::visEdges(smooth = list(type = "continuous"))

  list(
    network = vis,
    meta    = list(max_nes = round(max_nes, 1), has_substrates = TRUE)
  )
}


#' Map NES to color: purple (negative) → grey → orange (positive)
#' Distinct from .fc_to_color (blue → grey → red) to avoid confusion.
#' @param nes     Numeric vector of NES values
#' @param max_nes Maximum |NES| for saturation
#' @return Character vector of hex colors
.nes_to_color <- function(nes, max_nes = 3) {
  max_nes <- max(max_nes, 0.1)
  vapply(nes, function(v) {
    if (is.na(v)) return("#DCDCDC")
    v <- max(min(v, max_nes), -max_nes)
    if (v > 0) {
      t <- v / max_nes
      # grey (#DCDCDC) → orange (#E67E22)
      r <- as.integer(220 + t * (230 - 220))
      g <- as.integer(220 - t * (220 - 126))
      b <- as.integer(220 - t * (220 - 34))
    } else {
      t <- abs(v) / max_nes
      # grey (#DCDCDC) → purple (#8E44AD)
      r <- as.integer(220 - t * (220 - 142))
      g <- as.integer(220 - t * (220 - 68))
      b <- as.integer(220 - t * (220 - 173))
    }
    sprintf("#%02x%02x%02x", r, g, b)
  }, "")
}


#' Generate HTML legend for kinase-substrate network
#'
#' @param max_nes Maximum |NES| shown in color scale
#' @return shiny tagList
kinase_substrate_legend_html <- function(max_nes = 3) {
  max_nes <- max(round(max_nes, 1), 0.5)

  # Kinase NES color gradient (purple → grey → orange)
  nes_vals   <- seq(-max_nes, max_nes, length.out = 11)
  nes_colors <- .nes_to_color(nes_vals, max_nes = max_nes)
  nes_grad   <- .css_gradient(nes_colors)

  # Substrate FC color gradient (blue → grey → red)
  fc_vals   <- seq(-3, 3, length.out = 11)
  fc_colors <- .fc_to_color(fc_vals, max_fc = 3)
  fc_grad   <- .css_gradient(fc_colors)

  bar_style   <- "height:12px; border-radius:3px; border:1px solid #ccc; "
  label_style <- "font-size:13px; font-weight:600; color:#444; margin-bottom:2px;"
  tick_style  <- "font-size:11px; color:#666;"

  htmltools::tagList(
    htmltools::tags$div(
      style = "padding:8px 12px; font-family:sans-serif; max-width:320px;",

      # Kinase NES color bar
      htmltools::tags$div(style = label_style,
        htmltools::tags$span(
          style = "display:inline-block; width:12px; height:12px; background:#d5d8dc; border:2px solid #666; transform:rotate(45deg); margin-right:6px; vertical-align:middle;"
        ),
        "Kinase node: NES"
      ),
      htmltools::tags$div(style = paste0(bar_style, "background:", nes_grad, ";")),
      htmltools::tags$div(
        style = "display:flex; justify-content:space-between; margin-bottom:10px;",
        htmltools::tags$span(style = tick_style, paste0("-", max_nes)),
        htmltools::tags$span(style = tick_style, "0"),
        htmltools::tags$span(style = tick_style, paste0("+", max_nes))
      ),

      # Substrate FC color bar
      htmltools::tags$div(style = label_style,
        htmltools::tags$span(
          style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#d5d8dc; border:1px solid #999; margin-right:6px; vertical-align:middle;"
        ),
        "Substrate node: log2 FC"
      ),
      htmltools::tags$div(style = paste0(bar_style, "background:", fc_grad, ";")),
      htmltools::tags$div(
        style = "display:flex; justify-content:space-between; margin-bottom:10px;",
        htmltools::tags$span(style = tick_style, "-3"),
        htmltools::tags$span(style = tick_style, "0"),
        htmltools::tags$span(style = tick_style, "+3")
      ),

      # Edge legend
      htmltools::tags$div(
        style = "font-size:12px; color:#555;",
        htmltools::tags$span(
          style = "display:inline-block; width:20px; height:2px; background:#888888; vertical-align:middle; margin-right:4px;"
        ),
        htmltools::tags$span("\u2192 kinase \u2192 substrate")
      )
    )
  )
}
