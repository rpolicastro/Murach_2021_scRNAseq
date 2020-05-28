
## Load singularity container.
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh && conda activate seurat && R

library("Seurat")
library("tidyverse")
library("clustree")
library("future")
library("unixtools")
library("cerebroApp")
library("loomR")
library("velocyto.R")
library("SeuratWrappers")

###########################
## Seurat Analysis Kevin ##
###########################

## Loading Data
## ----------

## Samples to analyze.

velocyto_samples <- list(
	"tdT_Parental" = "tdt_parental",
	"KY_Mononuclear" = "KY_mononuclear",
	"Pax7_tdT_4day" = "pax7_tdt_4day",
	"Pax7_DTA_4day" = "pax7_dta_4day"
)

velocyto_samples <- map(velocyto_samples, function(x) {
	x <- file.path("aligned", "new", x, "velocyto", str_c(x, ".loom"))
	return(x)
})

## Read in data.

velocyto_counts <- map(velocyto_samples, ReadVelocity)

## Convert to seurat object.

seurat_obj <- map(velocyto_counts, as.Seurat)

seurat_obj <- imap(seurat_obj, function(x, y) {
	x@meta.data$orig.ident <- y
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
	meta_data <- x@meta.data
	meta_data <- as_tibble(meta_data, .name_repair = "unique")
	return(meta_data)
})

mt_data <- bind_rows(mt_data)

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_spliced)) +
	geom_point(size = 0.1) +
	facet_wrap(. ~ orig.ident, ncol = 3, scales = "free")

if (!dir.exists(file.path("results", "preprocessing"))) {
	dir.create(file.path("results", "preprocessing"))
}

png(
	file.path("results", "preprocessing", "mt_content_spliced.png"),
	height = 10, width = 10, res = 300, units = "in"
)
p; dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- map(seurat_obj, function(x) {
	x <- subset(x, subset = percent.mt <= 25 & nFeature_spliced >= 750)
	return(x)
})

## Save the seurat object after cell quality control.

if (!dir.exists(file.path("results", "r_objects"))) {
	dir.create(file.path("results", "r_objects"))
}

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))

## Integrating Data
## ----------

## SCTransform to normalize data.

seurat_obj <- map(seurat_obj, ~ SCTransform(., assay = "spliced"))

saveRDS(seurat_obj, file.path("results", "r_objects", "seurat_obj_spliced.RDS"))

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

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Dimensionality Reduction and Clustering
## ----------

## PCA dimension reduction for clustering.

seurat_integrated <- RunPCA(seurat_integrated, npcs = 100)

## Elbow plot of PCA dimensions.

p <- ElbowPlot(seurat_integrated, ndims = 100)

if (!dir.exists(file.path("results", "clustering"))) {
	dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot_spliced.pdf"), height = 5, width = 5)
p; dev.off()

## Clustering the data.

seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:35)
seurat_integrated <- FindClusters(
	seurat_integrated, resolution = seq(0.2, 1.2, 0.1),
	method = "igraph", algorithm = 4, weights = TRUE
)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plotting a cluster tree.

p <- clustree(seurat_integrated, prefix = "integrated_snn_res.")

pdf(file.path("results", "clustering", "cluster_tree_spliced.pdf"), height = 16, width = 12)
p ;dev.off()

## Switch identity to a presumptive good clustering resolution.

Idents(seurat_integrated) <- "integrated_snn_res.0.5"

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:35)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

## Plot Clusters.

p <- DimPlot(
	seurat_integrated, group.by = "ident", split.by = "orig.ident",
	ncol = 2, label = TRUE
)

pdf(file.path("results", "clustering", "clusters_spliced.pdf"), height = 12, width = 10)
p; dev.off()
