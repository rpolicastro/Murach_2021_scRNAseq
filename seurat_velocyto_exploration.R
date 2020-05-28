
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

###################################
## Exploration of scRNA-seq Data ##
###################################

## Load Integrated Data
## ----------

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

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
