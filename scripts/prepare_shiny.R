
library("data.table")
library("tidyverse")
library("DBI")
library("RSQLite")
library("Seurat")


############################
## Prepare Data for Shiny ##
############################

## Load seurat object.

seurat_obj <- readRDS(file.path("results/r_objects/seurat_integrated.RDS"))

## Connect to SQLite Server.

con <- dbConnect(SQLite(), "/N/slate/rpolicas/kevin_scRNAseq_shiny/scRNAseqShiny/data/murach.sqlite")

## Prepare count data.

iwalk(seurat_obj, function(x, y) {
  counts <- x[["SCT"]]@counts
  counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")

  counts <- melt(
    counts, id.vars = "cell_id", variable.name = "gene",
    value.name = "exp"
  )

  copy_to(
    con, counts, str_c(y, "_counts"), temporary = FALSE, overwrite = TRUE,
    indexes = list("gene")
  )
})

## Prepare meta-data.

iwalk(seurat_obj, function(x, y) {
  meta_data <- as.data.table(x[[]], keep.rownames = "cell_id")

  copy_to(
    con, meta_data, str_c(y, "_metadata"), temporary = FALSE, overwrite = TRUE,
    indexes = list("seurat_clusters", "orig.ident")
  )
})

## Prepare UMAP.

iwalk(seurat_obj, function(x, y) {
  reductions <- as.data.table(
    Embeddings(x, "umap"), keep.rownames = "cell_id",
    key = "cell_id"
  )

  copy_to(
    con, reductions, str_c(y, "_reductions"), temporary = FALSE, overwrite = TRUE
  )
})

## Markers.

markers <- fread(file.path("results", "markers", "marker_table.tsv"))
markers <- split(markers, by = "experiment", keep.by = FALSE)

iwalk(markers, function(x, y) {
  markers <- x[,
    .(cluster, gene, pct.1, pct.2, p_val,
    p_val_adj, avg_logFC, avg_log2FC)
  ]

  copy_to(
    con, markers, str_c(y, "_markers"), temporary = FALSE, overwrite = TRUE,
    indexes = list("cluster", "gene")
  )
})

## Prepare enrichment data.

enriched <- list(
	GO = "/N/project/sc_sequencing/kevin_scRNAseq/results/enrichment/go_enrichment.tsv",
	Reactome = "/N/project/sc_sequencing/kevin_scRNAseq/results/enrichment/reactome_enrichment.tsv"
)

enriched <- map(enriched, fread, sep = "\t")
enriched <- rbindlist(enriched, idcol = "database")
enriched <- split(enriched, by = "experiment", keep.by = FALSE)

iwalk(enriched, function(x, y) {
  copy_to(
    con, x, str_c(y, "_enriched"), temporary = FALSE, overwrite = TRUE,
    indexes = list("database", "group")
  )
})

## Add general data to database.

# Samples.
samples <- imap(seurat_obj, ~data.table(
  experiment = .y, samples = unique(.x$orig.ident)
))
samples <- rbindlist(samples)

copy_to(con, samples, "samples", temporary = FALSE, overwrite = TRUE)

## Turn of db connection.

dbDisconnect(con)
