
## Prepare Signularity Container
## singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:seurat_velocytor_0.3
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
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
plan("multiprocess", workers = 2)

## Prepare Counts
## ----------

sample_files <- list(
	"KY_mononuclear" = file.path("aligned", "new", "KY_mononuclear", "outs", "filtered_feature_bc_matrix"),
	"Pax7_DTA_4Day" = file.path("aligned", "new", "pax7_dta_4day", "outs", "filtered_feature_bc_matrix"),
	"Pax7_tdT_4Day" = file.path("aligned", "new", "pax7_tdt_4day", "outs", "filtered_feature_bc_matrix"),
	"tdT_Parental" = file.path("aligned", "new", "tdt_parental", "outs", "filtered_feature_bc_matrix")
)

## Load raw counts.

counts_10X <- map(sample_files, Read10X)

## Create seurat object.

seurat_obj <- imap(counts_10X, function(x, y) {
	x <- CreateSeuratObject(counts = x, project = y, min.cells = 10, min.features = 250)
	return(x)
})

## Cell Quality Control
## ----------

## Add mitochondrial percentage.

seurat_obj <- map(seurat_obj, function(x) {
	x[["percent.mt"]] <- PercentageFeatureSet(x, "^mt-")
	return(x)
})

## Plot mitochondrial percentage versus feature counts.

mt_data <- map(seurat_obj, function(x) {
	x <- as_tibble(x@meta.data, .name_repair = "unique")
	return(x)
}) %>%
	bind_rows(.id = "sample")

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_RNA)) +
	geom_point(size = 0.1) +
	facet_wrap(~ sample, ncol = 2)

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"), recursive = TRUE)
}

pdf(file.path("results", "preprocessing", "mt_content.pdf"), height = 4, width = 4)
p
dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- map(seurat_obj, function(x) {
	x <- subset(x, subset = percent.mt <= 25 & nFeature_RNA >= 750)
	return(x)
})

## Save the seurat object after cell quality control.

if (!dir.exists(file.path("results", "r_objects"))) {
	dir.create(file.path("results", "r_objects"))
}

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj.RDS"))

## Normamlize and Integrate Data
## ----------

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, SCTransform)

## Prepare for integration.

integration_features <- SelectIntegrationFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- PrepSCTIntegration(seurat_obj, anchor.features = integration_features)

## Integrate the reference dataset.

reference_datasets <- which(names(seurat_obj) == "tdT_Parental")
integration_anchors <- FindIntegrationAnchors(
	seurat_obj, normalization.method = "SCT", anchor.features = integration_features,
	reference = reference_datasets
)

seurat_integrated <- IntegrateData(integration_anchors, normalization.method = "SCT")

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))

## Dimensionality Reduction and Clustering
## ----------

## PCA dimension reduction for clustering.

seurat_integrated <- RunPCA(seurat_integrated, npcs = 100)

## Elbow plot to explore PCA dimensions.

p <- ElbowPlot(seurat_integrated, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
	dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot.pdf"), height = 5, width = 5)
p
dev.off()

## Clustering the data.

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:35)

seurat_integrated <- FindClusters(
	seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
	algorithm = 4, method = "igraph", weights = TRUE
)

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree.pdf"), height = 16, width = 12)
p
dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.5"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:35)

## Plot Clusters.

p <- DimPlot(seurat_integrated, group.by = "ident", pt.size = 0.1)

pdf(file.path("results", "clustering", "clusters.pdf"), height = 5, width = 7.5)
p
dev.off()

## Save seurat object with the dimension reduction and clusters.

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))
