
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

counts <- Assays(seurat_obj, "SCT")@counts
counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")

counts <- melt(
  counts, id.vars = "cell_id", variable.name = "gene",
  value.name = "exp"
)

fwrite(
	counts, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/counts.csv",
	sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, counts, "counts", temporary = FALSE, overwrite = TRUE,
	indexes = list("gene"),
)

## Prepare meta-data.

meta_data <- as.data.table(seurat_obj[[]], keep.rownames = "cell_id")
meta_data[, seurat_clusters := integrated_snn_res.0.5]

fwrite(
        meta_data, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/metadata.csv",
        sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, meta_data, "metadata", temporary = FALSE, overwrite = TRUE,
	indexes = list("seurat_clusters", "orig.ident")
)

## Prepare UMAP.

reductions <- as.data.table(
	Embeddings(seurat_obj, "umap"), keep.rownames = "cell_id",
	key = "cell_id"
)

fwrite(
        reductions, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/umap.csv",
        sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(con, reductions, "reductions", temporary = FALSE, overwrite = TRUE)

## Markers.

markers <- fread(file.path("results", "markers", "marker_table.tsv"))

markers <- markers[,
	.(cluster, gene, pct.1, pct.2, p_val,
	p_val_adj, avg_logFC, avg_log2FC)
]

fwrite(
	markers, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/markers.csv",
	sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, markers, "markers", temporary = FALSE, overwrite = TRUE,
	indexes = list("cluster", "gene")
)

## Prepare enrichment data.

enriched <- list(
	GO = "/N/project/sc_sequencing/kevin_scRNAseq/results/enrichment/go_enrichment.tsv",
	Reactome = "/N/project/sc_sequencing/kevin_scRNAseq/results/enrichment/reactome_enrichment.tsv"
)

enriched <- map(enriched, ~fread(., sep = "\t"))
enriched <- rbindlist(enriched, idcol = "database")

fwrite(
	enriched, "/N/slate/rpolicas/kevin_scRNAseq_shiny/scRNAseqShiny/data/enriched.csv",
	sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, enriched, "enriched", temporary = FALSE, overwrite = TRUE,
	indexes = list("database", "group")
)

## Turn of db connection.

dbDisconnect(con)
