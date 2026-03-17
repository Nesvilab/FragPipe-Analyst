###############################################################################
## Knowledge Based Analysis Functions
## Provides: read_gmt, plot_gsva_heatmap, plot_ppi_network,
##           plot_kegg_pathway_network
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
      # Built-in database via test_ora_mod (single contrast)
      res <- test_ora_mod(dep, databases = database, contrasts = TRUE,
                          direction = dir, log2_threshold = log2_threshold,
                          alpha = alpha, adjust_alpha = adjust_alpha,
                          backend = backend)
      if (is.null(res) || nrow(res) == 0) return(NULL)
      # Filter to the requested contrast
      if ("contrast" %in% colnames(res))
        res <- dplyr::filter(res, contrast == !!contrast)
      if (nrow(res) == 0) return(NULL)
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
#  PTM-SEA (PTM Signature Enrichment Analysis) via fast.ssgsea
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

#' Run PTM-SEA on site-level DE results
#'
#' Computes a signed statistic (sign(log2FC) * -log10(p)) per site, maps sites
#' to PTMsigDB flanking format, and runs fast.ssgsea.
#'
#' @param dep       SummarizedExperiment after test_diff() (site-level)
#' @param contrast  Contrast name (e.g., "A_vs_B")
#' @param species   "human", "mouse", or "rat"
#' @param mod_type  PTM modification type suffix for PTMsigDB matching.
#'                  One of "p" (phospho), "ac" (acetyl), "ub" (ubiquitin),
#'                  "me" (methyl), "sm" (sumoyl), etc.
#' @param nperm     Number of permutations for p-value estimation
#' @param min_size  Minimum gene set size
#' @return data.frame from fast_ssgsea (set, set_size, ES, NES, p_value, adj_p_value, etc.)
run_ptmsea <- function(dep, contrast, species = "human", mod_type = "p",
                       nperm = 1000L, min_size = 5L) {
  rd <- as.data.frame(SummarizedExperiment::rowData(dep))

  # Validate SequenceWindow column exists
  if (!"SequenceWindow" %in% colnames(rd))
    stop("SequenceWindow column missing — PTM-SEA requires site-level data ",
         "with flanking sequence information.")

  diff_col <- paste0(contrast, "_diff")
  p_col    <- paste0(contrast, "_p.val")
  if (!diff_col %in% colnames(rd) || !p_col %in% colnames(rd))
    stop("Contrast columns not found for: ", contrast)

  lfc   <- rd[[diff_col]]
  pval  <- rd[[p_col]]
  seqw  <- rd$SequenceWindow

  # Build signed statistic: sign(logFC) * -log10(p)
  # Use a floor for p-values to avoid Inf
  pval[pval < 1e-300] <- 1e-300
  stat <- sign(lfc) * (-log10(pval))

  # Name by SequenceWindow + modification type marker
  # PTMsigDB members are e.g. "ARQSRRSTQGVTLTD-p;d" for phospho,
  # "ARQSRRSTQGVTLTD-ac;u" for acetyl — the suffix must match.
  mod_suffix <- paste0("-", mod_type)
  names(stat) <- paste0(toupper(seqw), mod_suffix)

  # Filter out invalid entries
  keep <- !is.na(stat) & !is.na(seqw) & nchar(trimws(seqw)) > 0
  stat <- stat[keep]

  # Handle duplicate flanking sequences: keep largest |stat|
  if (anyDuplicated(names(stat))) {
    df_dedup <- data.frame(nm = names(stat), val = stat, absval = abs(stat),
                           stringsAsFactors = FALSE)
    df_dedup <- df_dedup[order(-df_dedup$absval), ]
    df_dedup <- df_dedup[!duplicated(df_dedup$nm), ]
    stat <- setNames(df_dedup$val, df_dedup$nm)
  }

  if (length(stat) < 10)
    stop("Too few valid sites (", length(stat), ") for PTM-SEA.")

  # Load PTMsigDB
  gmt_path <- get_ptmsigdb_path(species)
  if (!file.exists(gmt_path))
    stop("PTMsigDB file not found: ", gmt_path)
  gene_sets <- fast.ssgsea::read_gmt(gmt_path)

  # Run fast_ssgsea
  fast.ssgsea::fast_ssgsea(
    stats     = stat,
    gene_sets = gene_sets,
    alpha     = 1,
    nperm     = as.integer(nperm),
    min_size  = as.integer(min_size),
    seed      = 42L
  )
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


#' PPI Network via STRING REST API
#'
#' Fetches protein-protein interactions for significant DE genes from the
#' STRING database REST API and renders an interactive visNetwork graph.
#'
#' @param dep            SummarizedExperiment object after test_diff()
#' @param contrast       Name of the contrast (e.g. "CondA_vs_CondB")
#' @param lfc_threshold  Absolute log2 fold-change cutoff
#' @param alpha          Adjusted p-value cutoff
#' @param string_score   Minimum combined STRING interaction score (0–1000)
#' @return visNetwork object
plot_ppi_network <- function(dep, contrast, lfc_threshold = 1, alpha = 0.05,
                             string_score = 400) {
  row_data <- as.data.frame(SummarizedExperiment::rowData(dep))
  diff_col <- paste0(contrast, "_diff")
  padj_col <- paste0(contrast, "_p.adj")

  # Resolve gene symbols
  if ("Gene" %in% colnames(row_data)) {
    row_data$gene_symbol <- row_data$Gene
  } else {
    row_data$gene_symbol <- gsub("[.].*", "", row_data$name)
  }

  de_genes <- row_data %>%
    dplyr::filter(!is.na(.data[[diff_col]]) & !is.na(.data[[padj_col]])) %>%
    dplyr::filter(abs(.data[[diff_col]]) > lfc_threshold,
                  .data[[padj_col]] < alpha) %>%
    dplyr::select(gene_symbol,
                  diff  = dplyr::all_of(diff_col),
                  padj  = dplyr::all_of(padj_col)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "")

  if (nrow(de_genes) == 0) {
    nodes <- data.frame(id = 1, label = "No significant DE genes found",
                        color = "#95a5a6", stringsAsFactors = FALSE)
    return(visNetwork::visNetwork(nodes, data.frame()))
  }

  # --- STRING REST API ---
  identifiers <- paste(de_genes$gene_symbol, collapse = "%0d")
  api_url <- paste0(
    "https://string-db.org/api/json/network",
    "?identifiers=", identifiers,
    "&species=9606",
    "&required_score=", string_score,
    "&caller_identity=FragPipe-Analyst"
  )

  resp <- tryCatch(
    httr::GET(api_url, httr::timeout(30)),
    error = function(e) NULL
  )

  edges <- data.frame()
  if (!is.null(resp) && httr::status_code(resp) == 200) {
    # STRING returns Content-Type: text/json which httr has no automatic
    # parser for — read as text and parse explicitly.
    net_json <- rjson::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
    if (length(net_json) > 0) {
      edges <- data.frame(
        from  = sapply(net_json, function(x) x$preferredName_A),
        to    = sapply(net_json, function(x) x$preferredName_B),
        width = sapply(net_json, function(x) as.numeric(x$score) * 5),
        title = sapply(net_json, function(x) paste0("score: ", round(as.numeric(x$score), 3))),
        stringsAsFactors = FALSE
      )
    }
  }

  # Build nodes (all DE genes; highlight those in the network)
  network_genes <- unique(c(edges$from, edges$to))
  nodes <- de_genes %>%
    dplyr::mutate(
      id    = gene_symbol,
      label = gene_symbol,
      title = paste0(gene_symbol,
                     "<br>log2FC: ", round(diff, 3),
                     "<br>p.adj: ",  signif(padj, 3)),
      color = ifelse(diff > 0, "#e74c3c", "#3498db"),
      value = pmax(abs(diff) * 3, 3),
      borderWidth = ifelse(gene_symbol %in% network_genes, 2, 1)
    ) %>%
    dplyr::select(id, label, title, color, value, borderWidth)

  # Add any extra genes returned by STRING that are not in our DE list
  extra <- setdiff(network_genes, nodes$id)
  if (length(extra) > 0) {
    extra_nodes <- data.frame(
      id = extra, label = extra, title = extra,
      color = "#95a5a6", value = 5, borderWidth = 1,
      stringsAsFactors = FALSE
    )
    nodes <- dplyr::bind_rows(nodes, extra_nodes)
  }

  visNetwork::visNetwork(nodes, edges,
                         main = list(text = paste0("PPI: ", contrast),
                                     style = "font-size:13px; font-weight:bold;")) %>%
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = FALSE) %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           stabilization = list(iterations = 150)) %>%
    visNetwork::visEdges(smooth = list(type = "continuous")) %>%
    visNetwork::visLegend(
      addNodes = data.frame(
        label = c("Up-regulated", "Down-regulated", "Not in network"),
        color = c("#e74c3c", "#3498db", "#95a5a6"),
        stringsAsFactors = FALSE
      ),
      useGroups = FALSE
    )
}


#' KEGG Pathway Internal Network
#'
#' Fetches the KGML for a given KEGG pathway and renders an interactive
#' visNetwork graph where nodes are genes and edges are biological relations
#' (activation, inhibition, phosphorylation, etc.).  DE fold-change from the
#' current contrast is overlaid as node colour and size.
#'
#' @param pathway_id  KEGG pathway ID, e.g. "hsa04010" (from enrichment ID col)
#' @param dep         SummarizedExperiment object after test_diff()
#' @param contrast    Name of the contrast
#' @return visNetwork object
plot_kegg_pathway_network <- function(pathway_id, dep, contrast) {
  # --- Fetch KGML ---
  kgml_raw <- tryCatch(
    KEGGREST::keggGet(pathway_id, "kgml"),
    error = function(e) stop("Failed to fetch KEGG KGML for ", pathway_id, ": ", e$message)
  )
  kgml_str <- if (is.list(kgml_raw)) kgml_raw[[1]] else kgml_raw

  doc <- tryCatch(
    xml2::read_xml(kgml_str),
    error = function(e) stop("Failed to parse KEGG KGML: ", e$message)
  )

  # --- Parse entries (nodes) ---
  entries <- xml2::xml_find_all(doc, ".//entry")
  if (length(entries) == 0) {
    nodes <- data.frame(id = 1, label = "No pathway entries found", color = "#95a5a6")
    return(visNetwork::visNetwork(nodes, data.frame()))
  }

  entry_rows <- lapply(entries, function(e) {
    graphics   <- xml2::xml_find_first(e, ".//graphics")
    gene_names <- if (!is.na(graphics)) xml2::xml_attr(graphics, "name") else NA_character_
    data.frame(
      id         = xml2::xml_attr(e, "id"),
      kegg_name  = xml2::xml_attr(e, "name"),
      type       = xml2::xml_attr(e, "type"),
      label_raw  = gene_names,
      stringsAsFactors = FALSE
    )
  })
  entry_df <- do.call(rbind, entry_rows)

  # Keep only gene entries
  gene_df <- entry_df[!is.na(entry_df$type) & entry_df$type == "gene", ]
  if (nrow(gene_df) == 0) {
    nodes <- data.frame(id = 1, label = "No gene entries in pathway", color = "#95a5a6")
    return(visNetwork::visNetwork(nodes, data.frame()))
  }

  # Extract first gene symbol from "GENE1, GENE2, ..." label
  gene_df$label <- sapply(gene_df$label_raw, function(x) {
    if (is.na(x) || x == "") return("?")
    trimws(strsplit(x, ",")[[1]][1])
  })

  # --- Overlay DE data ---
  row_data <- as.data.frame(SummarizedExperiment::rowData(dep))
  diff_col <- paste0(contrast, "_diff")
  padj_col <- paste0(contrast, "_p.adj")

  if ("Gene" %in% colnames(row_data)) {
    row_data$gene_symbol <- row_data$Gene
  } else {
    row_data$gene_symbol <- gsub("[.].*", "", row_data$name)
  }

  de_map <- row_data %>%
    dplyr::select(gene_symbol,
                  diff_val = dplyr::all_of(diff_col),
                  padj_val = dplyr::all_of(padj_col)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  gene_df <- gene_df %>%
    dplyr::left_join(de_map, by = c("label" = "gene_symbol"))

  # Build visNetwork nodes
  nodes <- gene_df %>%
    dplyr::mutate(
      title = paste0(
        label,
        dplyr::if_else(
          !is.na(diff_val),
          paste0("<br>log2FC: ", round(diff_val, 3),
                 "<br>p.adj: ",  signif(padj_val, 3)),
          "<br>(not DE)"
        )
      ),
      color = dplyr::case_when(
        !is.na(diff_val) & diff_val > 0  ~ "#e74c3c",
        !is.na(diff_val) & diff_val < 0  ~ "#3498db",
        TRUE                              ~ "#95a5a6"
      ),
      value = dplyr::if_else(!is.na(diff_val), pmax(abs(diff_val) * 4, 3), 5)
    ) %>%
    dplyr::select(id, label, title, color, value)

  # --- Parse relations (edges) ---
  relations <- xml2::xml_find_all(doc, ".//relation")
  edge_color_map <- c(
    "activation"         = "#27ae60",
    "inhibition"         = "#c0392b",
    "phosphorylation"    = "#8e44ad",
    "dephosphorylation"  = "#d35400",
    "ubiquitination"     = "#16a085",
    "binding/association"= "#f39c12",
    "indirect effect"    = "#bdc3c7",
    "state change"       = "#1abc9c",
    "compound"           = "#e67e22",
    "hidden compound"    = "#e67e22",
    "missing interaction"= "#ecf0f1",
    "unknown"            = "#95a5a6"
  )

  if (length(relations) > 0) {
    rel_rows <- lapply(relations, function(r) {
      subtypes   <- xml2::xml_find_all(r, ".//subtype")
      subtype_nm <- if (length(subtypes) > 0) xml2::xml_attr(subtypes[[1]], "name") else "unknown"
      data.frame(
        from    = xml2::xml_attr(r, "entry1"),
        to      = xml2::xml_attr(r, "entry2"),
        rel_type = xml2::xml_attr(r, "type"),
        subtype  = subtype_nm,
        stringsAsFactors = FALSE
      )
    })
    edges_raw <- do.call(rbind, rel_rows)
    # Keep only gene-to-gene edges
    edges_raw <- edges_raw[edges_raw$from %in% gene_df$id &
                           edges_raw$to   %in% gene_df$id, ]

    edges <- edges_raw %>%
      dplyr::mutate(
        color = dplyr::if_else(
          subtype %in% names(edge_color_map),
          edge_color_map[subtype], "#95a5a6"
        ),
        title = subtype,
        arrows = dplyr::if_else(
          subtype %in% c("activation", "inhibition",
                         "phosphorylation", "dephosphorylation", "indirect effect"),
          "to", "")
      ) %>%
      dplyr::select(from, to, color, title, arrows)
  } else {
    edges <- data.frame(from = character(), to = character(),
                        color = character(), title = character(),
                        arrows = character())
  }

  visNetwork::visNetwork(
    nodes, edges,
    main = list(text  = pathway_id,
                style = "font-size:13px; font-weight:bold;")
  ) %>%
    visNetwork::visOptions(highlightNearest = TRUE) %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           stabilization = list(iterations = 100)) %>%
    visNetwork::visLegend(
      addNodes = data.frame(
        label = c("Up-regulated", "Down-regulated", "Not significant"),
        color = c("#e74c3c", "#3498db", "#95a5a6"),
        stringsAsFactors = FALSE
      ),
      useGroups = FALSE
    )
}
