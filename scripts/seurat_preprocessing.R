
## Prepare Singularity Container
## singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:seurat_velocytor_0.3
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("data.table")
library("clustree")
library("future")
library("unixtools")

#####################
## Kevin scRNA-seq ##
#####################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

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
  CreateSeuratObject(counts = x, project = y, min.cells = 10, min.features = 250)
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
  as.data.table(x[[]], keep.rownames = "cell_id")
})
mt_data <- rbindlist(mt_data, idcol = "sample")

p <- ggplot(mt_data, aes(x = percent.mt, y = nFeature_RNA)) +
  geom_point(size = 0.1) +
  facet_wrap(~ sample, ncol = 2)

if (!dir.exists(file.path("results", "preprocessing"))) {
  dir.create(file.path("results", "preprocessing"), recursive = TRUE)
}

pdf(file.path("results", "preprocessing", "mt_content.pdf"), height = 4, width = 4)
p; dev.off()

## Filter the data based on number of features and mitochondrial content.

seurat_obj <- map(
  seurat_obj, subset,
  subset = percent.mt <= 25 & nFeature_RNA >= 750
)

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

seurat_obj <- list(
  "tdT_Expression" = list(
    "tdT_Parental" = seurat_obj[["tdT_Parental"]],
    "tdT_Expressed" = seurat_obj[["KY_mononuclear"]]
  ),
  "tdT_4Day" = list(
    "Pax7_tdT_4Day" = seurat_obj[["Pax7_tdT_4Day"]],
    "Pax7_DTA_4Day" = seurat_obj[["Pax7_DTA_4Day"]]
 )  
)

integration_features <- map(seurat_obj, SelectIntegrationFeatures, nfeatures = 3000)
seurat_obj <- map2(seurat_obj, integration_features, function(x, y) {
  PrepSCTIntegration(x, anchor.features = y)
})

## Integrate the reference dataset.

integration_anchors <- map2(seurat_obj, integration_features, function(x, y) {
  anchors <- FindIntegrationAnchors(
    x, normalization.method = "SCT", anchor.features = y,
    reference = which(names(x) %in% c("tdT_Parental", "Pax7_tdT_4Day"))
  )
  return(anchors)
})

seurat_integrated <- map(integration_anchors, IntegrateData, normalization.method = "SCT")

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))

## Dimensionality Reduction and Clustering
## ----------

## PCA dimension reduction for clustering.

seurat_integrated <- map(seurat_integrated, RunPCA, npcs = 100)

## Elbow plot to explore PCA dimensions.

if (!dir.exists(file.path("results", "clustering"))) {
  dir.create(file.path("results", "clustering"))
}

pdf(file.path("results", "clustering", "pca_elbow_plot.pdf"), height = 5, width = 5)
 imap(seurat_integrated, function(x, y) {
    ElbowPlot(x, ndims = 100) + ggtitle(y)
  })
dev.off()

## Clustering the data.

seurat_integrated <- map(seurat_integrated, FindNeighbors, dims = 1:30)

seurat_integrated <- map(
  seurat_integrated, FindClusters,
  resolution = seq(0.2, 1.0, 0.1),
  algorithm = 4, method = "igraph", weights = TRUE
)

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))

## Plotting a cluster tree.

pdf(file.path("results", "clustering", "cluster_tree.pdf"), height = 16, width = 12)
  imap(seurat_integrated, function(x, y) {
    clustree(x, prefix = "integrated_snn_res.") + ggtitle(y)
  })
dev.off()

## Switch identity to a presumptive good clustering resolution.

seurat_integrated <- imap(seurat_integrated, function(x, y) {
  if (y == "tdT_Expression") {
    Idents(x) <- "integrated_snn_res.0.5"
    x$seurat_clusters <- x$integrated_snn_res.0.5
  } else {
    Idents(x) <- "integrated_snn_res.0.4"
    x$seurat_clusters <- x$integrated_snn_res.0.4
  }
  return(x)
})

## UMAP dimension reduction for visualization.

if (!dir.exists("tempdir")) dir.create("tempdir")
set.tempdir("tempdir")

seurat_integrated <- map(seurat_integrated, RunUMAP, dims = 1:30)

## Plot Clusters.

pdf(file.path("results", "clustering", "clusters.pdf"), height = 5, width = 7.5)
  imap(seurat_integrated, function(x, y) {
    DimPlot(x, group.by = "ident", pt.size = 0.1, label = TRUE) + ggtitle(y)
  })
dev.off()

pdf(file.path("results", "clustering", "clusters_per_sample.pdf"), height = 6, width = 12)
  imap(seurat_integrated, function(x, y) {
    DimPlot(
	x, group.by = "ident", pt.size = 0.1,
	split.by = "orig.ident", label = TRUE, ncol = 2
    ) +
    ggtitle(y)
  })
dev.off()

## Save seurat object with the dimension reduction and clusters.

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_integrated.RDS"))
