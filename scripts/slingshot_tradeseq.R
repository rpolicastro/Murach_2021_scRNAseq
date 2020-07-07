
library("tidyverse")
library("Seurat")
library("slingshot")
library("tradeSeq")
library("tidymodels")
library("ranger")
library("janitor")
library("patchwork")
library("gam")
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

###################
## Random Forest ##
###################

## Get Counts.

gene_counts <- map(seurat_obj, function(x) {
  x <- GetAssayData(x, assay = "SCT", slot = "data")
  x <- as.data.table(t(as.matrix(x)), keep.rownames = "cell_id")
  return(x)
})

## Merge with pseudotimes.

lineages <- map(pseudotime, as.data.table, keep.rownames = "cell_id")

gene_counts <- map2(gene_counts, lineages, merge, by = "cell_id")
gene_counts <- map(gene_counts, clean_names)

## Split the data for curve1.

split_data <- map(gene_counts, function(x) {
  x <- x[!is.na(curve1)]
  x <- initial_split(x)
  x <- list(
    training = training(x),
    testing = testing(x)
  )
  return(x)
})

## Make the model.

models <- map(split_data, function(x) {

  x <- x[["training"]]

  remove_cols <- keep(colnames(x), ~str_detect(.x, "^curve|cell_id"))

  rf_formula <- as.formula(str_c(
    "curve1 ~ . - ", str_c(remove_cols, collapse = " - "),
    sep = " "
  ))

  model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
    set_engine("ranger", importance = "impurity", num.threads = 8) %>%
    fit(rf_formula, data = x)

  return(model)

})

##############
## GAM Gene ##
##############

## Gene counts.

gene_counts <- map(seurat_obj, function(x) {
  x <- GetAssayData(x, assay = "SCT", slot = "data")
  x <- as.data.table(t(as.matrix(x)), keep.rownames = "cell_id")
  return(x)
})

og_names <- map(gene_counts, function(x) {
  x <- data.table(original = colnames(x))
  return(x)
})

gene_counts <- map(gene_counts, clean_names)

og_names <- map2(og_names, gene_counts, function(x, y) {
  x[["gene"]] <- colnames(y)
  return(x)
})

## Pseudotimes.

lineages <- map(pseudotime, as.data.table, keep.rownames = "cell_id")

## Combine required data.

model_data <- imap(gene_counts, function(x, y) {
  genes <- keep(colnames(x), ~!str_detect(., "cell_id"))

  x <- merge(x, lineages[[y]], by = "cell_id")
  x <- merge(x, metadata[[y]], by = "cell_id")

  data_list <- list(
    counts = x, genes = genes,
    lineages = keep(colnames(x), ~str_detect(., "curve\\d+$"))
  )
  return(data_list)
})

## Run the models.

model_results <- map(model_data, function(x) {

  count_data <- x[["counts"]][!is.na(curve1)]

  pvals <- map(x[["genes"]], function(y) {
    model_formula <- as.formula(str_c(y, "~ lo(curve1)", sep = " "))
    model <- gam(model_formula, data = count_data)
    pval <- summary(model)$anova$`Pr(F)`[2]
    return(data.table(pvalue = pval, gene = y))
  })

  pvals <- rbindlist(pvals)
  return(pvals)

})

model_results <- map2(model_results, og_names, function(x, y) {
  x <- merge(x, y, by = "gene")
  x[, gene := NULL]
  setnames(x, old = "original", new = "gene")
  return(x)
})

## Correct for multiple comparisons.

model_results <- map(model_results, function(x) {
  x <- x[, .(gene, pvalue, FDR = p.adjust(pvalue, "fdr"))]
  return(x)
})

################################

Y <- log1p(assays(sim)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
})

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

## Calculate cell weights.

cell_weights <- map(trajectory, slingCurveWeights)

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
