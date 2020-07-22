
## Prepare Signularity Container
## singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:seurat_velocytor_0.3
##
## singularity shell -eCB `pwd` -H `pwd` scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh && conda activate seurat && R

library("Seurat")
library("tidyverse")
library("data.table")
library("scProportionTest")
library("future")

##################################
## Exploring the scRNA-seq Data ##
##################################

## Loading the Data
## ----------

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

##########################
## tdTomato Exploration ##
##########################

## Exploring tdTomato
## ----------

if (!dir.exists(file.path("results", "tdTomato"))) {
	dir.create(file.path("results", "tdTomato"))
}

## Plotting the expression of tdTomato and tdTomato 5' UTR.

DefaultAssay(seurat_integrated) <- "SCT"

p <- FeaturePlot(
	seurat_integrated, feature = c("tdTomato", "tdTomatoStop", "Pax7"),
	pt.size = 0.1, ncol = 2, split.by = "orig.ident"
)

pdf(file.path("results", "tdTomato", "tdTomato_dimplot.pdf"), height = 12, width = 16)
p
dev.off()

## Modeling tdTomato Expression
## ----------

## preparing counts.

raw_counts <- map(seurat_integrated, function(x) {
  x <- x[["SCT"]]@counts
  counts <- as.data.table(t(as.matrix(x)), keep.rownames = "cell_id")
  counts <- counts[, .(cell_id, tdTomato, tdTomatoStop)]
  counts[, tdT_Log2_Ratio := log2(tdTomato + 1) - log2(tdTomatoStop + 1)]
  return(counts)
})

## Grabbing cell meta data.

meta_data <- map(seurat_integrated, function(x) {
  metadata <- as.data.table(x[[]], keep.rownames = "cell_id")
  metadata <- metadata[, .(cell_id, orig.ident)]
  return(metadata)
})

## Merge the meta.data into the cell data.

merged <- map2(raw_counts, meta_data, function(x, y) {
  merged <- merge(y, x, by = "cell_id")
  return(merged)
})

merged <- rbindlist(merged, idcol = "experiment")

## PU Learning
## Code courtesy of Professor Daniel McDonald
## Based on Jaskie et al., 2018 (DOI: 10.1109/IEEECONF44664.2019.9048765)
## ----------

## PU learning functions.

weird_logit <- function(theta, s, x){
  p = 1 / (1 + theta[1]^2 + exp(-(theta[2] + theta[3]*x)))
  -mean(dbinom(s, 1, p, TRUE))
}

nontrad_class <- function(dat){
  y = dat$tdT_Log2_Ratio
  n = length(y)
  s = as.numeric(dat$orig.ident =="tdT_Parental")
  out = optim(c(1,0,0), weird_logit, method="BFGS", s=s, x=y)
  theta = out$par
  chat = 1 / (1+theta[1]^2)
  ps0 = 1 / (1 + theta[1]^2 + exp(-(theta[2] + theta[3]*y)))
  py0 = ps0 / chat
  py0
}

## Get probability values.

merged[, PU_learning_pvalue :=  nontrad_class(merged)]
merged[,
  PU_learning_FDR := p.adjust(PU_learning_pvalue, "fdr"),
  by = orig.ident
]


## Export the results table.

fwrite(
	merged, file.path("results", "tdTomato", "tdT_Ratio_PU_Learning.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Prediction Interval Test
## ----------

## Preparing the negative control set.

neg_control <- merged[orig.ident == "tdT_Parental"][, .(orig.ident, cell_id, tdT_Log2_Ratio)]

## prediction interval test.

exact_results <- map(unique(merged[["orig.ident"]]), function(x) {
	pred_interval <- pmap_df(merged[orig.ident == x], function(...) {
		row_args <- data.table(...)
		pval <- sum(neg_control[["tdT_Log2_Ratio"]] >= row_args[["tdT_Log2_Ratio"]]) / nrow(neg_control)
		row_args[["Pred_interval_pvalue"]] <- pval
		return(row_args)
	})
	pred_interval[, Pred_interval_FDR := p.adjust(Pred_interval_pvalue, "fdr")]
	return(pred_interval)
})

merged <- rbindlist(exact_results)

## Export results table.

fwrite(
	exact_results, file.path("results", "tdTomato", "tdT_Ratio.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Add results back to seurat object.
## ----------

## Prepare data.

merged[,
  c("PU_learning_sig", "Pred_interval_sig") := list(
    PU_learning_FDR < 0.05, Pred_interval_FDR < 0.05
  )
]

merged <- split(merged, by = "experiment", keep.by = FALSE)

row_order <- map2(merged, seurat_integrated, function(x, y) {
  x <- x[,
    .(cell_id, tdT_Log2_Ratio, PU_learning_pvalue, PU_learning_FDR,
    PU_learning_sig, Pred_interval_pvalue, Pred_interval_FDR,
    Pred_interval_sig
    )
  ][order(match(cell_id, rownames(y@meta.data)))]
  x[, cell_id := NULL]
  return(x)
})

seurat_integrated[[1]][[]] <- cbind(
  seurat_integrated[[1]][[]],
  row_order[[1]]
)

## Add test info back into seurat object.

seurat_integrated <- map2(seurat_integrated, row_order, function(x, y) {
  x@meta.data <- cbind(x[[]], y)
  return(x)
})

## Save expanded seurat object.

saveRDS(seurat_integrated, file.path("results", "r_objects", "seurat_expanded.RDS"))

##########################
## Cell Proportion Test ##
##########################

if (!dir.exists(file.path("results", "cell_counts"))) {
	dir.create(file.path("results", "cell_counts"))
}

## Compare cell counts for each cluster between samples.

sc_utils_obj <- sc_utils(seurat_integrated)

comparisons <- list(
        c("tdT_Parental", "KY_mononuclear"),
        c("Pax7_tdT_4Day", "Pax7_DTA_4Day")
)

walk(comparisons, function(x) {
        sc_utils_obj <- permutation_test(
                sc_utils_obj, cluster_identity = "integrated_snn_res.0.5",
                sample_1 = x[1], sample_2 = x[2]
        )

        p <- permutation_plot(sc_utils_obj, log2FD_threshold = 1)

        file_name <- str_c(x[1], "_vs_", x[2], ".pdf")
        pdf(file.path("results", "cell_counts", file_name), height = 3, width = 8)
        print(p); dev.off()
})

#####################
## Marker Analysis ##
#####################

options(future.globals.maxSize = 10000 * 1024 ^2)
plan("multiprocess", workers = 4)

## Find Markers
## ----------

## Find markers for all clusters.

markers <- map(seurat_integrated, function(x) {
  markers <- FindAllMarkers(
    x, assay = "SCT", slot = "data",
    logfc.threshold = log(1.5), min.pct = 0.25, return.thresh = 0.05
  )
  return(markers)
})

## Process markers.

markers <- map(markers, function(x) {
  setDT(x)
  x[, avg_log2FC := log2(exp(avg_logFC))]
  x <- x[p_val_adj < 0.05]
  x <- x[order(cluster, p_val_adj)]
  return(x)
})

## Save markers.

if (!dir.exists(file.path("results", "markers"))) {
        dir.create(file.path("results", "markers"))
}

markers <- rbindlist(markers, idcol = "experiment")
fwrite(
	markers, file.path("results", "markers", "marker_table.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Marker Plots
## ----------

## Get top markers.

top_markers <- markers[avg_log2FC >= log2(1.5)]
top_markers <- top_markers[order(cluster, -avg_log2FC)]
top_markers <- top_markers[, head(.SD, 5), by = cluster]
top_markers <- top_markers[["gene"]]

## Marker heatmap.

p <- DoHeatmap(
	seurat_integrated, features = top_markers, group.by = "integrated_snn_res.0.5",
	assay = "SCT", slot = "data"
) +
	scale_fill_viridis_c()

pdf(file.path("results", "markers", "marker_heatmap.pdf"), width = 16, height = 10)
p; dev.off()

#####################
## Term Enrihcment ##
#####################

library("clusterProfiler")
library("ReactomePA")
library("org.Mm.eg.db")

if (!dir.exists(file.path("results", "enrichment"))) {
	dir.create(file.path("results", "enrichment"))
}

## Load up the markers and prepare for analysis.

markers <- fread(file.path("results", "markers", "marker_table.tsv"), sep = "\t")
markers <- split(markers, by = "experiment", keep.by = FALSE)

markers <- map(markers, function(x) {
  x[, change := fifelse(avg_log2FC > 0, "up", "down")]
  x[,
	group := str_c("cluster", cluster, change, sep = "_"),
	by = seq_len(nrow(x))
  ]

  markers <- split(x, x[["group"]])
  return(markers)
})

## GO analysis.

go_enrichment <- map(markers, function(x) {
  go_enrichment <- map(x, function(y) {
    enriched <- enrichGO(
      gene = y[["gene"]], OrgDb = "org.Mm.eg.db",
      keyType = "SYMBOL", ont = "BP"
    )
    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  go_enrichment <- rbindlist(go_enrichment, idcol = "group")
  return(go_enrichment)
})

go_enrichment <- rbindlist(go_enrichment, idcol = "experiment")

fwrite(
  go_enrichment, file.path("results", "enrichment", "go_enrichment.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Reactome analysis.

reactome_enrichment <- map(markers, function(x) {
  enriched <- map(x, function(y) {
    entrez <- bitr(
      y[["gene"]], fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db"
    )
    entrez <- entrez[["ENTREZID"]]

    enriched <- enrichPathway(gene = entrez, organism = "mouse", readable = TRUE)

    enriched <- as.data.table(enriched)
    enriched <- enriched[p.adjust < 0.05]
    return(enriched)
  })

  enriched <- rbindlist(enriched, idcol = "group")
  return(enriched)
})

reactome_enrichment <- rbindlist(reactome_enrichment, idcol = "experiment")

fwrite(
        reactome_enrichment, file.path("results", "enrichment", "reactome_enrichment.tsv"),
        sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
