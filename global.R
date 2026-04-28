library("SummarizedExperiment")
library("tidyverse")
library("testthat")
library("shiny")
library("shinydashboard")
library("shinyjs")
library("shinyalert")
library("ComplexHeatmap")
library("dplyr")
library("limma")
library("DT")
library("ggrepel")
library("httr")
library("rjson")
library("svglite")
library("ensembldb")
library(EnsDb.Hsapiens.v86)
library("conflicted")
library("plotly")
library("shinyWidgets")# new added
library("ggVennDiagram") # new added
library("rhandsontable") # new added
library("shinyBS") # new added
library("shinycssloaders") # new added
library("shiny.info")
library("fastcluster")
library("factoextra")
library("UpSetR")
library("MSnbase")
library("clusterProfiler")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("ReactomePA")
library(pcaMethods) # for bpca imputation
library(missForest) # for RF imputation
library(vegan)
library(assertthat)
library(RColorBrewer)
library(data.table)

conflict_prefer("box", "shinydashboard")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")
conflict_prefer("setdiff", "base")
conflict_prefer("unite", "tidyr")
conflict_prefer("intersect", "base")
conflict_prefer("colMedians", "matrixStats")
conflicts_prefer(base::as.factor)

source("R/functions.R")
source("R/volcano_function.R")
source("R/customized.R")
source("R/tests.R")
source("R/demo_functions.R")
source("R/enrichment_functions.R")
source("R/filter.R")

# Pre-load local GMT gene sets for ORA enrichment
.read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  terms <- vapply(parts, `[[`, character(1), 1)
  genes <- lapply(parts, function(p) p[seq_along(p)[-(1:2)]])
  data.frame(
    gs_name     = rep(terms, lengths(genes)),
    gene_symbol = unlist(genes, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

LOADED_GMT_FILES <- local({
  gmt_dir <- "data/gene_sets"
  files   <- list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
  if (length(files) == 0) {
    message("No GMT files in data/gene_sets/. Run update_gene_sets.R to build them.")
    return(list())
  }
  message("Loading ", length(files), " GMT gene set files...")
  setNames(lapply(files, function(f) {
    df <- .read_gmt(f)
    df$gene_symbol <- toupper(df$gene_symbol)
    df
  }), sub("\\.gmt$", "", basename(files)))
})

GO_TERM2NAME <- tryCatch(
  suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys    = AnnotationDbi::keys(GO.db::GO.db, "GOID"),
    columns = "TERM", keytype = "GOID")),
  error = function(e) { message("GO.db unavailable; GO terms will show IDs."); NULL })
