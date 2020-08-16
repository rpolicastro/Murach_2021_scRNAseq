
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

# 14 Day Unspliced.
seurat_obj$tdT_Expression$seurat_clusters <- seurat_obj$tdT_Expression$integrated_snn_res.0.5 %>% {case_when(
    . == 1 ~ "FAPs",
    . == 2 ~ "Fibroblasts",
    . == 3 ~ "Myonuclei",
    . == 4 ~ "FAPs",
    . == 5 ~ "Myeloid (H2-Ab1+ H2-Eb1+)",
    . == 6 ~ "Macrophages (Mrc1 High)",
    . == 7 ~ "Endothelial (Pecam1+)",
    . == 8 ~ "Satellite Cells",
    . == 9 ~ "Myeloid",
    . == 10 ~ "Tenocytes",
    . == 11 ~ "Myeloid",
    . == 12 ~ "Unknown",
    . == 13 ~ "Cycling Cells",
    . == 14 ~ "Smooth Muscle Cells",
    . == 15 ~ "Satellite Cells (MyoG+)",
    . == 16 ~ "Cycling Endothelial Cells"
)}
seurat_obj$tdT_Expression$seurat_clusters_merged <- str_replace(
  seurat_obj$tdT_Expression$seurat_clusters,
  "\\s\\(.+\\)$", ""
)

# 14 Day Spliced.
seurat_velocity$tdT_Expression$seurat_clusters <- seurat_velocity$tdT_Expression$integrated_snn_res.0.5 %>%
{case_when(
    . == 1 ~ "FAPs",
    . == 2 ~ "FAPs (Tnxb High)",
    . == 3 ~ "Myeloid (H2-Ab1+ H2-Eb1+)",
    . == 4 ~ "Fibroblasts",
    . == 5 ~ "Macrophages (Mrc1 High)",
    . == 6 ~ "Endothelial Cell (Pecam1+)",
    . == 7 ~ "Myonuclei",
    . == 8 ~ "Satellite Cells",
    . == 9 ~ "Tenocytes",
    . == 10 ~ "Myeloid",
    . == 11 ~ "Myeloid",
    . == 12 ~ "Unknown",
    . == 13 ~ "Cycling Cells",
    . == 14 ~ "Satellite Cells (MyoG+)",
    . == 15 ~ "Smooth Muscle Cells",
    . == 16 ~ "Cycling Endothelial Cells"
)}
seurat_velocity$tdT_Expression$seurat_clusters_merged <- str_replace(
  seurat_velocity$tdT_Expression$seurat_clusters,
  "\\s\\(.+\\)$", ""
)

# 4 Day Spliced.
seurat_velocity$tdT_4Day$seurat_clusters <- seurat_velocity$tdT_4Day$integrated_snn_res.0.4 %>% {case_when(
    . == 1 ~ "Macrophages",
    . == 2 ~ "Myeloid (H2-Eb1- H2-Ab1-)",
    . == 3 ~ "FAPs (Pdgfra and Lbp High)",
    . == 4 ~ "Myeloid (H2-Eb1+ H2-Ab1+)",
    . == 5 ~ "Neutrophils (s100a8/9 High)",
    . == 6 ~ "Neutrophils (s100a8/9 Low)",
    . == 7 ~ "Myeloid (H2-Eb1+ H2-Eb1+)",
    . == 8 ~ "B Cells (Cd74+ Ccr7+)",
    . == 9 ~ "Endothelial (Pecam1+)",
    . == 10 ~ "Myonuclei (Myh4+/Myh1+)",
    . == 11 ~ "Tregs (Ms4a4b+)",
    . == 12 ~ "Cycling Monocytes (Pclaf High)",
    . == 13 ~ "Myeloid (Cst3 High)",
    . == 14 ~ "Fibrogenic (Pdgfra-)",
    . == 15 ~ "Pericytes (Rgs5+)",
    . == 16 ~ "Tenocytes"
)}
seurat_velocity$tdT_4Day$seurat_clusters_merged <- str_replace(
  seurat_velocity$tdT_4Day$seurat_clusters,
  "\\s\\(.+\\)$", ""
)

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
enriched <- split(enriched, by="experiment", keep.by=FALSE)

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

