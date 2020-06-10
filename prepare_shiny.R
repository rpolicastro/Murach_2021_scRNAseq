
library("data.table")
library("tidyverse")
library("DBI")
library("odbc")
library("Seurat")

############################
## Prepare Data for Shiny ##
############################

## Load seurat object.

seurat_obj <- readRDS(file.path("results/r_objects/seurat_integrated.RDS"))

## Connect to MySQL Server.

con <- dbConnect(
        odbc(),
        Driver = "/usr/lib64/libmyodbc5.so",
        database = "rpolicas_murach",
        Server = "108.167.189.70",
        UID = "",
        PWD = "",
        Port = 3306
)

## Prepare count data.

counts <- Assays(seurat_obj, "SCT")@counts
counts <- as.data.table(t(as.matrix(counts)), keep.rownames = "cell_id")

counts <- melt(
  counts, id.vars = "cell_id", variable.name = "gene",
  value.name = "exp"
)

gene_ids <- data.table(gene = as.character(unique(counts[["gene"]])), key = "gene")
gene_ids[, gene_unid := seq_len(.N)]

cell_ids <- data.table(cell_id = as.character(unique(counts[["cell_id"]])), key = "cell_id")
cell_ids[, cell_unid := seq_len(.N)]

ids <- list(genes = gene_ids, cells = cell_ids)

setkey(counts, gene, cell_id)
counts <- merge(counts, ids[["cells"]], by = "cell_id")
counts <- merge(counts, ids[["genes"]], by = "gene")

fwrite(
	counts, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/counts.csv",
	sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, counts, "counts", temporary = FALSE, overwrite = TRUE,
	indexes = list("gene_unid", "cell_unid"),
)

## Prepare meta-data.

meta_data <- as.data.table(seurat_obj[[]], keep.rownames = "cell_id")

ident <- data.table(orig.ident = as.character(unique(meta_data[["orig.ident"]])), key = "orig.ident")
ident[, ident_unid := seq_len(.N)]

ids <- c(ids, list(ident = ident))

setkey(meta_data, cell_id, orig.ident)
meta_data <- merge(meta_data, ids[["cells"]], by = "cell_id")
meta_data <- merge(meta_data, ids[["ident"]], by = "orig.ident")

fwrite(
        meta_data, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/metadata.csv",
        sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, meta_data, "metadata", temporary = FALSE,
	overwrite = TRUE, indexes = list("cell_unid", "ident_unid")
)

## Prepare UMAP.

reductions <- as.data.table(
	Embeddings(seurat_obj, "umap"), keep.rownames = "cell_id",
	key = "cell_id"
)

reductions <- merge(reductions, ids[["cells"]], by = "cell_id")

fwrite(
        reductions, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/umap.csv",
        sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, reductions, "reductions", temporary = FALSE, overwrite = TRUE,
	indexes = list("cell_unid")
)

## Markers.

markers <- fread(file.path("results", "markers", "marker_table.tsv"))

markers <- markers[,
	.(cluster, gene, pct.1, pct.2, p_val,
	p_val_adj, avg_logFC, avg_log2FC)
]

markers <- merge(markers, ids[["genes"]], by = "gene")

fwrite(
	markers, "/N/slate/rpolicas/kevin_scRNAseq_shiny/data/markers.csv",
	sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
)

copy_to(
	con, markers, "markers", temporary = FALSE, overwrite = TRUE,
	indexes = list("cluster", "gene_unid")
)
