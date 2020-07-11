
library("tidyverse")
library("Seurat")
library("slingshot")
library("tradeSeq")
library("patchwork")
library("data.table")

#########################
## Trajectory Analysis ##
#########################

setwd("..")

## Loading the data.

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

## Slingshot
## ----------

## Trajectory analysis.

trajectory <- map(seurat_obj, function(x) {
  traj <- slingshot(
    Embeddings(x, "umap"), clusterLabels = x$seurat_clusters
  )
  return(traj)
})

saveRDS(
  trajectory,
  file.path("results", "r_objects", "slingshot_trajectory.RDS")
)

## Calculate pseudotime over the various lineages.

pseudotime <- map(trajectory, slingPseudotime)

## Calculate the sling curves.

curves <- map(trajectory, slingCurves)

crv <- as.data.table(slingCurves(trajectory[[1]])[[1]]$s[slingCurves(trajectory[[1]])[[1]]$ord, ])

p <- ggplot() +
  geom_line(data = as.data.frame(slingCurves(trajectory[[1]])[[1]]$s[slingCurves(trajectory[[1]])[[1]]$ord, ]), aes(x=UMAP_1, y=UMAP_2))

#########################
## ggplot trajectories ##
#########################

## Preparing Data
## ----------

## Find centers and connectivity.

trj_data <- map(trajectory, function(x) {

  X <- reducedDim(x)
  clusterLabels <- slingClusterLabels(x)

  connectivity <- slingAdjacency(x)
  clusters <- rownames(connectivity)

  centers <- t(vapply(clusters,function(clID){
    w <- clusterLabels[,clID]
    return(apply(X, 2, weighted.mean, w = w))
    },
    rep(0,ncol(X))
  ))

  rownames(centers) <- clusters
  X <- X[rowSums(clusterLabels) > 0, , drop = FALSE]
  clusterLabels <- clusterLabels[rowSums(clusterLabels) > 0, , drop = FALSE]

  trj_data <- list(
    connectivity = connectivity,
    centers = as.data.table(centers, keep.rownames = "clust"),
    umap = as.data.table(X, keep.rownames = "cell_id"),
    clusters = clusters
  )

  return(trj_data)

})

## Add meta-data to umap coords.

metadata <- map(seurat_obj, function(x) {
  x <- as.data.table(x[[]], keep.rownames = "cell_id")[,
    .(cell_id, orig.ident, seurat_clusters)
  ]
  return(x)
})

trj_data <- map2(trj_data, metadata, function(x, y) {
  x[["umap"]] <- merge(y, x[["umap"]], by = "cell_id")
  return(x)
})

## Add curves.

trj_data <- imap(trj_data, function(x, y) {

  curves <- slingCurves(trajectory[[y]])
  curves <- map(curves, function(j) {
    curve <- as.data.table(j$s[j$ord, ], keep.rownames = "cell_id")
  })
  curves <- rbindlist(curves, idcol = "curve")

  x[["curves"]] <- curves

  return(x)
 
})

## Prepare line segment data.table.

trj_data <- map(trj_data, function(x) {
  combos <- gtools::combinations(length(x[["clusters"]]), 2, x[["clusters"]])

  segments <- apply(combos, 1, function(y) {
    if (x[["connectivity"]][y[1], y[2]] == 1) {
      DT <- data.table(
        xstart = x[["centers"]][clust == y[1]]$UMAP_1,
        ystart = x[["centers"]][clust == y[1]]$UMAP_2,
        xend = x[["centers"]][clust == y[2]]$UMAP_1,
        yend = x[["centers"]][clust == y[2]]$UMAP_2
      )
      return(DT)
    }
  })

  segments <- purrr::compact(segments)
  segments <- rbindlist(segments)
  x$segments <- segments

  return(x)

})

## Plot Overall
## ----------

## Make base plot.

overall_plots <- imap(trj_data, function(x, y) {

  # Base plot.
  p <- ggplot(x[["umap"]], aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = seurat_clusters), size = 0.2) +
    theme_classic() +
    ggtitle(y)

  # Center points.
  p <- p + geom_point(data = x[["centers"]], aes(x = UMAP_1, y = UMAP_2), size = 3)

  # Lines.
  p <- p + geom_segment(data = x[["segments"]], aes(x = xstart, y = ystart, xend = xend, yend = yend))

  # Return plot.
  return(p)

})

p <- wrap_plots(overall_plots, ncol = 1)

ggsave(
  file.path("results", "trajectory", "slingshot_overall.pdf"),
  plot = p, device = cairo_pdf, height = 10, width = 8
)

## Make base curves plot.

curve_plots <- imap(trj_data, function(x, y) {

  # Base plot.
  p <- ggplot(x[["umap"]], aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = seurat_clusters), size = 0.2) +
    theme_classic() +
    ggtitle(y)

  # Lines.
  p <- p + geom_path(data = x[["curves"]], aes(x = UMAP_1, y = UMAP_2, group = curve))

  # Return plot.
  return(p)

})

curve_plots <- wrap_plots(curve_plots, ncol = 1)

ggsave(
  file.path("results", "trajectory", "slingshot_overall_curves.pdf"),
  plot = curve_plots, device = cairo_pdf, height = 10, width = 8
)

## Pseudotime Plots
## ----------

## Merge pseudotime data into umap coords and meta-data.

pseudotime_dt <- map(pseudotime, as.data.table, keep.rownames = "cell_id")

trj_data <- map2(trj_data, pseudotime_dt, function(x, y) {
  x[["umap"]] <- merge(x[["umap"]], y, by = "cell_id")
  return(x)
})

## Plot lineages.

lineage_plots <- imap(trj_data, function(x, y) {

  curves <- keep(colnames(x[["umap"]]), ~str_starts(., "curve"))

  lineages <- map(curves, function(z) {
    p <- ggplot(x[["umap"]], aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes_string(color = z), size = 0.2) +
      theme_classic() +
      scale_color_viridis_c(na.value = "grey") +
      ggtitle(str_c(z, y, sep = "_"))

    p <- p + geom_point(data = x[["centers"]], aes(x = UMAP_1, y = UMAP_2), size = 3)
    p <- p + geom_segment(data = x[["segments"]], aes(x = xstart, y = ystart, xend = xend, yend = yend))

    return(p)
  })

  lineages <- wrap_plots(lineages, ncol = 3)

  return(lineages)

})

lineage_plots <- wrap_plots(lineage_plots, ncol = 1)

ggsave(
  file.path("results", "trajectory", "slingshot_pseudotime.png"),
  plot = lineage_plots, device = "png", type = "cairo",
  dpi = 300, height = 18, width = 12
)

## Plot lineages with curves.

lineage_plots_curve <- imap(trj_data, function(x, y) {

  curves <- keep(colnames(x[["umap"]]), ~str_starts(., "curve"))
  curve_lines <- split(x[["curves"]], by = "curve")

  lineages <- map(curves, function(z) {
    p <- ggplot(x[["umap"]], aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes_string(color = z), size = 0.2) +
      theme_classic() +
      scale_color_viridis_c(na.value = "grey") +
      ggtitle(str_c(z, y, sep = "_"))

    p <- p + geom_path(data = curve_lines[[z]], aes(x = UMAP_1, y = UMAP_2))

    return(p)
  })

  lineages <- wrap_plots(lineages, ncol = 3)

  return(lineages)

})

lineage_plots_curve <- wrap_plots(lineage_plots_curve, ncol = 1)

ggsave(
  file.path("results", "trajectory", "slingshot_pseudotime_curves.png"),
  plot = lineage_plots_curve, device = "png", type = "cairo",
  dpi = 300, height = 18, width = 12
)

##############
## tradeSeq ##
##############

## Find the ideal number of knots.

pdf(file.path("results", "trajectory", "knot_plots.pdf"), height = 6, width = 14)
  knots <- map2(seurat_obj, trajectory, function(x, y) {
    knts <- evaluateK(
      counts = as.matrix(x[["RNA"]]@counts),
      sds = y, plot = TRUE, k = seq(3, 10),
      conditions = factor(x$orig.ident)
    )
    return(knts)
  })
dev.off()

saveRDS(knots, file.path("results", "r_objects", "knots.RDS"))

## Fit the negative binomial model.

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 6

nb_model <- map2(seurat_obj, trajectory, function(x, y) {
  nbm <- fitGAM(
    counts = as.matrix(x[["RNA"]]@counts), sds = y,
    conditions = factor(x$orig.ident), nknots = 6,
    parallel = TRUE, BPPARAM = BPPARAM, sce = TRUE
  )
  return(nbm)
})

nb_model <- map2(seurat_obj[1], trajectory[1], function(x, y) {
  nbm <- fitGAM(
    counts = as.matrix(x[["RNA"]]@counts), sds = y,
    conditions = factor(x$orig.ident), nknots = 6,
    genes = seq_len(10), sce = TRUE,
    parallel = TRUE, BPPARAM = BPPARAM
  )
  return(nbm)
})

