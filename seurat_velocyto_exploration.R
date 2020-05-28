
## Load singularity container.
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh && conda activate seurat && R

library("Seurat")
library("tidyverse")
library("data.table")
library("future")
library("wesanderson")
library("readxl")
library("scProportionTest")
library("velocyto.R")
library("SeuratWrappers")

###################################
## Exploration of scRNA-seq Data ##
###################################

## Load Integrated Data
## ----------

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

##################
## Marker Genes ##
##################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 6)

## Find Markers
## ----------

Idents(seurat_integrated) <- "integrated_snn_res.0.5"

markers <- FindAllMarkers(
	seurat_integrated, assay = "SCT", slot = "data", min.pct = 0.25,
	return.thresh = 0.05, logfc.threshold = log(1.5)
)

saveRDS(markers, file.path("results", "r_objects", "markers_spliced.RDS"))

## Load and prepare marker list.

setDT(markers)
markers <- markers[p_val_adj < 0.05 & abs(avg_logFC) > log(1.5)]
markers[, c("avg_log2FC", "cluster") := list(log2(exp(avg_logFC)), str_c("cluster_", cluster))]
markers <- markers[order(cluster, -avg_log2FC)]

## Split the list based on cluster and save results.

if (!dir.exists(file.path("results", "markers"))) {
        dir.create(file.path("results", "markers"), recursive = TRUE)
}

fwrite(
	markers, file.path("results", "markers", "marker_table_spliced.tsv"),
	col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t"
)

##################
## RNA Velocity ##
##################

if (!dir.exists(file.path("results", "rna_velocity"))) {
        dir.create(file.path("results", "rna_velocity"), recursive = TRUE)
}

## Run RNA velocity.

Idents(seurat_integrated) <- "integrated_snn_res.0.5"

walk(unique(seurat_integrated[["orig.ident"]][[1]]), function(sc_sample) {

	test_subset <- subset(seurat_integrated, subset = orig.ident == sc_sample, downsample = 250)

	seurat_velocity <- RunVelocity(
		test_subset, ambiguous = "ambiguous", ncores = 6,
		deltaT = 1, kCells = 25, fit.quantile = 0.02,
		group.by = "integrated_snn_res.0.5"
	)

	ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat_velocity)))
	names(x = ident.colors) <- levels(x = seurat_velocity)
	cell.colors <- ident.colors[Idents(object = seurat_velocity)]
	names(x = cell.colors) <- colnames(x = seurat_velocity)

	pdf(file.path("results", "rna_velocity", str_c(sc_sample, "_velocity.pdf")), height = 8, width = 8)
	show.velocity.on.embedding.cor(
		emb = Embeddings(object = seurat_velocity, reduction = "umap"),
		vel = Tool(object = seurat_velocity, slot = "RunVelocity"),
		n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
		cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE,
		min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE,
		cell.border.alpha = 0.1, n.cores = 6
	)
	dev.off()
})
