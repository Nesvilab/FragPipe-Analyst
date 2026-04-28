# update_gene_sets.R
#
# Maintains local gene set GMT files for ORA enrichment.
# - If files are missing: downloads/builds them.
# - If files exist: checks versions and updates only outdated ones.
#
# Usage:
#   source("update_gene_sets.R")                 # runs on source
#   check_genesets_updates()                     # re-run manually anytime

# ── Configuration ─────────────────────────────────────────────────────────────

DEST_DIR <- "data/gene_sets"

# ── Packages ──────────────────────────────────────────────────────────────────

required_pkgs <- c("KEGGREST", "org.Hs.eg.db", "org.Mm.eg.db", "AnnotationDbi", "dplyr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: BiocManager::install(c('",
       paste(missing_pkgs, collapse = "', '"), "'))")

suppressPackageStartupMessages({
  library(KEGGREST)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
})

# ── Helpers ───────────────────────────────────────────────────────────────────

write_gmt <- function(term2gene, file) {
  pathways <- split(term2gene$gene_symbol, term2gene$gs_name)
  pathways <- pathways[sapply(pathways, length) > 0]
  lines    <- sapply(names(pathways), function(nm)
    paste(c(nm, "NA", pathways[[nm]]), collapse = "\t"))
  writeLines(lines, file)
  message(sprintf("    Written %d pathways -> %s", length(lines), basename(file)))
}

build_go_term2gene <- function(orgdb, ontology) {
  suppressMessages(
    AnnotationDbi::select(orgdb,
      keys    = AnnotationDbi::keys(orgdb, "SYMBOL"),
      columns = c("SYMBOL", "GO", "ONTOLOGY"),
      keytype = "SYMBOL")
  ) %>%
    filter(ONTOLOGY == ontology, !is.na(GO)) %>%
    select(gs_name = GO, gene_symbol = SYMBOL) %>%
    distinct()
}

build_wikipathways_term2gene <- function(gmt_file, orgdb) {
  lines  <- readLines(gmt_file, warn = FALSE)
  parts  <- strsplit(lines, "\t", fixed = TRUE)
  # Pathway names arrive as "name%WikiPathways_DATE%WPxxx%species" — keep only the readable name
  terms  <- vapply(parts, function(p) gsub("%.*$", "", p[[1]]), character(1))
  eids   <- lapply(parts, function(p) p[-(1:2)])  # skip name and URL columns

  all_entrez <- unique(unlist(eids))
  message(sprintf("    Converting %d ENTREZID -> SYMBOL...", length(all_entrez)))
  mapping <- suppressMessages(
    AnnotationDbi::select(orgdb,
      keys    = all_entrez,
      columns = "SYMBOL",
      keytype = "ENTREZID")
  )

  data.frame(
    gs_name  = rep(terms, lengths(eids)),
    ENTREZID = unlist(eids, use.names = FALSE),
    stringsAsFactors = FALSE
  ) %>%
    left_join(mapping, by = "ENTREZID") %>%
    filter(!is.na(SYMBOL)) %>%
    select(gs_name, gene_symbol = SYMBOL) %>%
    distinct()
}

build_kegg_term2gene <- function(organism_kegg, orgdb) {
  message(sprintf("    Fetching gene-pathway links for '%s'...", organism_kegg))
  gene_pathway <- keggLink("pathway", organism_kegg)
  entrezids    <- gsub(paste0(organism_kegg, ":"), "", names(gene_pathway))
  pathway_ids  <- gsub("path:", "", as.character(gene_pathway))

  message("    Fetching pathway names...")
  pathway_list <- keggList("pathway", organism_kegg)
  pathway_df   <- data.frame(
    id   = gsub("path:", "", names(pathway_list)),
    name = gsub(" - .*$", "", as.character(pathway_list)),
    stringsAsFactors = FALSE)

  message("    Converting ENTREZID -> SYMBOL...")
  mapping <- suppressMessages(
    AnnotationDbi::select(orgdb,
      keys    = unique(entrezids),
      columns = "SYMBOL",
      keytype = "ENTREZID"))

  data.frame(id = pathway_ids, ENTREZID = entrezids, stringsAsFactors = FALSE) %>%
    left_join(mapping,    by = "ENTREZID") %>%
    left_join(pathway_df, by = "id") %>%
    filter(!is.na(SYMBOL), !is.na(name)) %>%
    select(gs_name = name, gene_symbol = SYMBOL) %>%
    distinct()
}

# ── Main ──────────────────────────────────────────────────────────────────────

check_genesets_updates <- function() {
  dir.create(DEST_DIR, showWarnings = FALSE, recursive = TRUE)
  versions_file <- file.path(DEST_DIR, "versions.rds")
  versions      <- if (file.exists(versions_file)) readRDS(versions_file) else list()

  # ── 1. WikiPathways ──────────────────────────────────────────────────────────
  message("\n[1/4] WikiPathways")
  for (cfg in list(
    list(label = "Homo_sapiens", key = "WikiPathways_human"),
    list(label = "Mus_musculus",  key = "WikiPathways_mouse")
  )) {
    tryCatch({
      base_url    <- "https://data.wikipathways.org/current/gmt/"
      listing     <- paste(readLines(base_url, warn = FALSE), collapse = "\n")
      filename    <- regmatches(listing,
                      regexpr(sprintf("wikipathways-[0-9]{8}-gmt-%s\\.gmt", cfg$label),
                              listing))
      if (!length(filename) || !nchar(filename))
        stop("filename not found in directory listing")
      current_ver <- gsub("wikipathways-([0-9]{8})-.*", "\\1", filename)
      stored      <- versions[[cfg$key]]
      file_ok     <- !is.null(stored) && file.exists(file.path(DEST_DIR, stored$file))

      if (file_ok && identical(stored$version, current_ver)) {
        message(sprintf("  %s: up to date (%s)", cfg$label, current_ver))
      } else {
        msg <- if (!file_ok) "not found" else
          sprintf("%s -> %s", stored$version, current_ver)
        message(sprintf("  %s: %s, downloading...", cfg$label, msg))
        if (!is.null(stored$file) && stored$file != filename)
          suppressWarnings(file.remove(file.path(DEST_DIR, stored$file)))
        tmp <- tempfile(fileext = ".gmt")
        download.file(paste0(base_url, filename), tmp, quiet = TRUE, mode = "wb")
        orgdb <- if (cfg$label == "Homo_sapiens") org.Hs.eg.db else org.Mm.eg.db
        t2g <- build_wikipathways_term2gene(tmp, orgdb)
        write_gmt(t2g, file.path(DEST_DIR, filename))
        unlink(tmp)
        versions[[cfg$key]] <- list(version    = current_ver,
                                     file       = filename,
                                     downloaded = format(Sys.time()))
        message(sprintf("    Done -> %s", filename))
      }
    }, error = function(e)
      message(sprintf("  ERROR (%s): %s", cfg$label, e$message)))
  }

  # ── 2. Reactome ──────────────────────────────────────────────────────────────
  message("\n[2/4] Reactome")
  tryCatch({
    current_ver <- trimws(readLines(
      "https://reactome.org/ContentService/data/database/version", warn = FALSE)[1])
    stored      <- versions$Reactome_human
    gmt_file    <- "Reactome_human.gmt"
    file_ok     <- file.exists(file.path(DEST_DIR, gmt_file))

    if (file_ok && identical(stored$version, current_ver)) {
      message(sprintf("  up to date (v%s)", current_ver))
    } else {
      msg <- if (!file_ok) "not found" else
        sprintf("v%s -> v%s", stored$version, current_ver)
      message(sprintf("  %s, downloading...", msg))
      zip_path <- file.path(DEST_DIR, "ReactomePathways.gmt.zip")
      download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",
                    zip_path, quiet = TRUE, mode = "wb")
      unzip(zip_path, files = "ReactomePathways.gmt", exdir = DEST_DIR)
      file.rename(file.path(DEST_DIR, "ReactomePathways.gmt"),
                  file.path(DEST_DIR, gmt_file))
      file.remove(zip_path)
      versions$Reactome_human <- list(version    = current_ver,
                                       file       = gmt_file,
                                       downloaded = format(Sys.time()))
      message(sprintf("    Done -> %s", gmt_file))
    }
  }, error = function(e) message(sprintf("  ERROR: %s", e$message)))

  # ── 3. KEGG ──────────────────────────────────────────────────────────────────
  message("\n[3/4] KEGG")
  tryCatch({
    kegg_info   <- keggInfo("pathway")
    m           <- regexpr("Release [0-9]+\\.[0-9]+[^\\n]*", kegg_info)
    current_ver <- trimws(regmatches(kegg_info, m))
    if (!nchar(current_ver)) stop("could not parse KEGG version")

    for (cfg in list(
      list(organism = "hsa", orgdb = org.Hs.eg.db,
           key = "KEGG_human", file = "KEGG_human.gmt"),
      list(organism = "mmu", orgdb = org.Mm.eg.db,
           key = "KEGG_mouse", file = "KEGG_mouse.gmt")
    )) {
      tryCatch({
        stored  <- versions[[cfg$key]]
        file_ok <- file.exists(file.path(DEST_DIR, cfg$file))

        if (file_ok && identical(stored$version, current_ver)) {
          message(sprintf("  %s: up to date (%s)", cfg$organism, current_ver))
        } else {
          msg <- if (!file_ok) "not found" else
            sprintf("%s -> %s", stored$version, current_ver)
          message(sprintf("  %s: %s, building...", cfg$organism, msg))
          t2g <- build_kegg_term2gene(cfg$organism, cfg$orgdb)
          write_gmt(t2g, file.path(DEST_DIR, cfg$file))
          versions[[cfg$key]] <- list(version    = current_ver,
                                       file       = cfg$file,
                                       downloaded = format(Sys.time()))
        }
      }, error = function(e)
        message(sprintf("  ERROR (%s): %s", cfg$organism, e$message)))
    }
  }, error = function(e) message(sprintf("  ERROR: %s", e$message)))

  # ── 4. GO ────────────────────────────────────────────────────────────────────
  message("\n[4/4] GO")
  go_meta   <- AnnotationDbi::metadata(org.Hs.eg.db)
  go_ver    <- go_meta[go_meta$name == "GOSOURCEDATE", "value"]
  orgdb_h_v <- as.character(packageVersion("org.Hs.eg.db"))
  orgdb_m_v <- as.character(packageVersion("org.Mm.eg.db"))
  message(sprintf("  GO source: %s | org.Hs.eg.db %s | org.Mm.eg.db %s",
                  go_ver, orgdb_h_v, orgdb_m_v))

  for (ont in c("BP", "MF", "CC")) {
    for (cfg in list(
      list(orgdb = org.Hs.eg.db, label = "human", orgdb_ver = orgdb_h_v),
      list(orgdb = org.Mm.eg.db, label = "mouse", orgdb_ver = orgdb_m_v)
    )) {
      key      <- sprintf("GO_%s_%s", ont, cfg$label)
      gmt_file <- sprintf("GO_%s_%s.gmt", ont, cfg$label)
      tryCatch({
        stored  <- versions[[key]]
        file_ok <- file.exists(file.path(DEST_DIR, gmt_file))

        if (file_ok && identical(stored$orgdb_version, cfg$orgdb_ver)) {
          message(sprintf("  GO %s %s: up to date (%s)", ont, cfg$label, cfg$orgdb_ver))
        } else {
          msg <- if (!file_ok) "not found" else
            sprintf("%s -> %s", stored$orgdb_version, cfg$orgdb_ver)
          message(sprintf("  GO %s %s: %s, building...", ont, cfg$label, msg))
          t2g <- build_go_term2gene(cfg$orgdb, ont)
          write_gmt(t2g, file.path(DEST_DIR, gmt_file))
          versions[[key]] <- list(version       = go_ver,
                                   orgdb_version = cfg$orgdb_ver,
                                   file          = gmt_file,
                                   downloaded    = format(Sys.time()))
        }
      }, error = function(e)
        message(sprintf("  ERROR (GO %s %s): %s", ont, cfg$label, e$message)))
    }
  }

  # ── Save versions ─────────────────────────────────────────────────────────
  saveRDS(versions, versions_file)
  message("\nDone. Gene set versions:")
  for (nm in names(versions)) {
    v <- versions[[nm]]
    message(sprintf("  %-30s %s", nm, v$version))
  }
  invisible(versions)
}

check_genesets_updates()
