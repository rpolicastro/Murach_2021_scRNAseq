
library("data.table")
library("tidyverse")
library("DBI")
library("RSQLite")
library("Seurat")


############################
## Prepare Data for Shiny ##
############################

## Load Seurat Objects
## ----------

## Load seurat object.

seurat_obj <- readRDS(file.path("results/r_objects/seurat_expanded.RDS"))
seurat_velocity <- readRDS(file.path("results/r_objects/seurat_integrated_spliced.RDS"))

## Annotate clusters.

seurat_obj$tdT_Expression$seurat_clusters <- seurat_obj$tdT_Expression$integrated_snn_res.0.5
seurat_obj$tdT_4Day$seurat_clusters <- seurat_obj$tdT_4Day$integrated_snn_res.0.4 %>% {case_when(
    . == 1 ~ "M2 Macrophages",
    . == 2 ~ "Neutrophils",
    . == 3 ~ "Neutrophils",
    . == 4 ~ "Myeloid",
    . == 5 ~ "Tenocytes",
    . == 6 ~ "Myeloid",
    . == 7 ~ "T Cells",
    . == 8 ~ "Dead/Dying",
    . == 9 ~ "Endothelial",
    . == 10 ~ "Myonuclei",
    . == 11 ~ "FAPs",
    . == 12 ~ "NK Cells",
    . == 13 ~ "Monocytes",
    . == 14 ~ "Myeloid",
    . == 15 ~ "Fibrogenic",
    . == 16 ~ "Pericytes"
)}

seurat_velocity$tdT_Expression$seurat_clusters <- seurat_velocity$tdT_Expression$integrated_snn_res.0.5
seurat_velocity$tdT_4Day$seurat_clusters <- seurat_velocity$tdT_4Day$integrated_snn_res.0.4 %>% {case_when(
    . == 1 ~ "M2 Macrophages",
    . == 2 ~ "Myeloid",
    . == 3 ~ "FAPs",
    . == 4 ~ "Myeloid",
    . == 5 ~ "Neutrophils",
    . == 6 ~ "Neutrophils",
    . == 7 ~ "Myeloid",
    . == 8 ~ "T Cells",
    . == 9 ~ "Endothelial",
    . == 10 ~ "Myonuclei",
    . == 11 ~ "NK Cells",
    . == 12 ~ "Monocytes",
    . == 13 ~ "Myeloid",
    . == 14 ~ "Fibrogenic",
    . == 15 ~ "Pericytes",
    . == 16 ~ "Tenocytes"
)}

## Make list of samples.

seu <- list(
  tdT_14Day_Normal = seurat_obj[["tdT_Expression"]],
  tdT_4Day_Normal = seurat_obj[["tdT_4Day"]],
  tdT_14Day_Spliced = seurat_velocity[["tdT_Expression"]],
  tdT_4Day_Spliced = seurat_velocity[["tdT_4Day"]]
)

## Save seurat object.

saveRDS(seu, file.path("results", "r_objects", "seurat_complete.RDS"))

## Connect to SQLite Server
## ----------

con <- dbConnect(SQLite(), "shiny_app/murach.sqlite")

## Prepare Gene Counts
## ----------

## Extract counts.

counts <- imap(seu, function(x, y) {
  counts <- x[["SCT"]]@counts
  counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")

  counts <- melt(
    counts, id.vars = "cell_id", variable.name = "gene",
    value.name = "exp"
  )

  return(counts)
})

## Save to database.

iwalk(counts, function(x, y) {
  copy_to(
    con, x, str_c(y, "_counts"), temporary = FALSE, overwrite = TRUE,
    indexes = list("gene")
  )
})

## Prepare Meta-Data.
## ----------

metadata <- imap(seu, function(x, y) {
  meta_data <- as.data.table(x[[]], keep.rownames = "cell_id")
  return(meta_data)
})

## Save to database.

iwalk(metadata, function(x, y) {
  copy_to(
    con, x, str_c(y, "_metadata"), temporary = FALSE, overwrite = TRUE,
    indexes = list("seurat_clusters", "orig.ident")
  )
})

## Prepare UMAP
## ----------

reductions <- map(seu, function(x) {
  reductions <- as.data.table(
    Embeddings(x, "umap"), keep.rownames = "cell_id",
    key = "cell_id"
  )

  return(reductions)
})

## Save to database.

iwalk(reductions, function(x, y) {
  copy_to(con, x, str_c(y, "_reductions"), temporary = FALSE, overwrite = TRUE)
})

## Prepare Markers.
## ----------

library("future")

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

## Find markers.

markers <- imap(seu, function(x, y) {

  Idents(x) <- "seurat_clusters"
  markers <- FindAllMarkers(
    x, assay = "SCT", slot = "data", logfc.threshold = log(1.5),
    min.pct = 0.25, return.thresh = 0.05
  )

  setDT(markers)
  markers[, avg_log2FC := log2(exp(avg_logFC))]
  markers <- markers[p_val_adj < 0.05]
  markers <- markers[order(cluster, p_val_adj)]
  markers <- markers[,
    .(cluster, gene, pct.1, pct.2, p_val,
    p_val_adj, avg_logFC, avg_log2FC)
  ]

  return(markers)
})

## Write out file.

if (!dir.exists("shiny_app")) dir.create("shiny_app")

fwrite(
  rbindlist(markers, idcol = "experiment"),
  file.path("shiny_app", "markers.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Save to database.

iwalk(markers, function(x, y) {
  copy_to(
    con, x, str_c(y, "_markers"), temporary = FALSE, overwrite = TRUE,
    indexes = list("cluster", "gene")
  )
})

## Prepare Enrichment
## ----------

library("clusterProfiler")
library("ReactomePA")
library("org.Mm.eg.db")

## prepare markers for analysis.

enr_markers <- map(markers, function(x) {
  x[, change := fifelse(avg_log2FC > 0, "up", "down")]
  x[,
        group := str_c("cluster", cluster, change, sep = "_"),
        by = seq_len(nrow(x))
  ]

  markers <- split(x, x[["group"]])
  return(markers)
})

## GO analysis.

go_enrichment <- map(enr_markers, function(x) {
  go_enrichment <- map(x, function(y) {
    enriched <- enrichGO(
      gene = y[["gene"]], OrgDb = "org.Mm.eg.db",
      keyType = "SYMBOL", ont = "BP"
    )
    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  go_enrichment <- rbindlist(go_enrichment, idcol = "group")
  return(go_enrichment)
})

go_enrichment <- rbindlist(go_enrichment, idcol = "experiment")

## Reactome analysis.

reactome_enrichment <- map(enr_markers, function(x) {
  enriched <- map(x, function(y) {
    entrez <- bitr(
      y[["gene"]], fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db"
    )
    entrez <- entrez[["ENTREZID"]]

    enriched <- enrichPathway(gene = entrez, organism = "mouse", readable = TRUE)

    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  enriched <- rbindlist(enriched, idcol = "group")
  return(enriched)
})

reactome_enrichment <- rbindlist(reactome_enrichment, idcol = "experiment")

## Save file of enrichment.

fwrite(
  rbindlist(list(go_enrichment, reactome_enrichment), idcol = "database"),
  file.path("shiny_app", "enriched.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Save to database.

enriched <- rbindlist(list(go_enrichment, reactome_enrichment), idcol = "database")

iwalk(enriched, function(x, y) {
  copy_to(
    con, x, str_c(y, "_enriched"), temporary = FALSE, overwrite = TRUE,
    indexes = list("database", "group")
  )
})

## Sample Sheet
## ----------

## Samples.

samples <- imap(seu, ~data.table(
  experiment = .y, samples = unique(.x$orig.ident)
))
samples <- rbindlist(samples)

copy_to(con, samples, "samples", temporary = FALSE, overwrite = TRUE)

## Turn of db connection.

dbDisconnect(con)

