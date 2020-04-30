
## Prepare Signularity Container
## singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:seurat_velocytor_0.3
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library("cerebroApp")

#####################
## Kevin scRNA-seq ##
#####################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

## Prepare Counts
## ----------

## Load raw counts.

counts_10X_dir <- file.path("aligned", "new", "KY_mononuclear", "outs", "filtered_feature_bc_matrix")
counts_10X <- Read10X(counts_10X_dir)

## Create seurat object.

seurat_obj <- CreateSeuratObject(
	counts = counts_10X, project = "kevin_scRNAseq",
	min.cells = 10, min.features = 250
)

## Cell Quality Control
## ----------

## Add mitochondrial percentage.

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, "^mt-")

## Plot mitochondrial percentage versus feature counts.

mt_data <- as_tibble(seurat_obj@meta.data, .name_repair = "unique")

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_RNA)) +
	geom_point(size = 0.1)

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"), recursive = TRUE)
}

pdf(file.path("results", "preprocessing", "mt_content.pdf"), height = 4, width = 4)
p
dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- subset(seurat_obj, subset = percent.mt <= 25 & nFeature_RNA >= 1000)

## Save the seurat object after cell quality control.

if (!dir.exists(file.path("results", "r_objects"))) {
	dir.create(file.path("results", "r_objects"))
}

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))

## Normamlize and Cluster Data
## ----------

## SCTransform to normalize data.

seurat_obj <- SCTransform(seurat_obj)

## PCA dimension reduction for clustering.

seurat_obj <- RunPCA(seurat_obj, npcs = 100)

## Elbow plot to explore PCA dimensions.

p <- ElbowPlot(seurat_obj, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
	dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot.pdf"), height = 5, width = 5)
p
dev.off()

## Clustering the data.

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

seurat_obj <- FindClusters(
	seurat_obj, resolution = seq(0.2, 1.2, 0.1),
	algorithm = 4
)

## Plotting a cluster tree.

p <- clustree(seurat_obj, prefix = "SCT_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree.pdf"), height = 16, width = 12)
p
dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_obj) <- "SCT_snn_res.0.7"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

## Plot Clusters.

p <- DimPlot(seurat_obj, group.by = "ident", pt.size = 0.1)

pdf(file.path("results", "clustering", "clusters.pdf"), height = 5, width = 7.5)
p
dev.off()

## Save seurat object with the dimension reduction and clusters.

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))
