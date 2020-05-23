
## Prepare Signularity Container
## singularity pull --arch amd64 library://rpolicastro/default/scrnaseq_software:seurat_velocytor_0.3
##
## singularity shell -eCB "$(pwd)" -H "$(pwd)" scrnaseq_software_seurat_velocytor_0.3.sif
##
## . /opt/conda/etc/profile.d/conda.sh
## conda activate seurat; R

library("Seurat")
library("tidyverse")
library("data.table")
library("scProportionTest")

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

DefaultAssay(seurat_integrated) <- "RNA"

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

raw_counts <- seurat_integrated@assays$RNA[
	rownames(seurat_integrated@assays$RNA) %in% c("tdTomato", "tdTomatoStop", "Pax7")
]

raw_counts <- raw_counts %>%
	as.data.frame %>%
	rownames_to_column("cell_id") %>%
	as.data.table %>%
	transpose(keep.names = "cell_id", make.names = 1)

raw_counts[, tdT_Log2_Ratio := log2(tdTomato + 1) - log2(tdTomatoStop + 1)]

## Grabbing cell meta data.

meta_data <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[,
	.(orig.ident, cell_id)
]

## Merge the meta.data into the cell data.

merged <- merge(meta_data, raw_counts, by = "cell_id")

## Preparing the negative control set.

neg_control <- merged[orig.ident == "tdT_Parental"][, .(orig.ident, cell_id, tdT_Log2_Ratio)]

## prediction interval test.

exact_results <- map(unique(merged[["orig.ident"]]), function(x) {
	pred_interval <- pmap_df(merged[orig.ident == x], function(...) {
		row_args <- data.table(...)
		pval <- sum(neg_control[["tdT_Log2_Ratio"]] >= row_args[["tdT_Log2_Ratio"]]) / nrow(neg_control)
		row_args[["pval"]] <- pval
		return(row_args)
	})
	pred_interval[, tdT_Ratio_FDR := p.adjust(pval, "fdr")]
	pred_interval[, tdT_Ratio_Sig := ifelse(tdT_Ratio_FDR < 0.05, TRUE, FALSE)]
	return(pred_interval)
})

exact_results <- rbindlist(exact_results)

## Export results table.

fwrite(
	exact_results, file.path("results", "tdTomato", "tdT_Ratio.tsv"),
	sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

## Add results back to seurat object.

plot_results <- exact_results[,
	.(cell_id, tdT_Log2_Ratio, tdT_Ratio_FDR, tdT_Ratio_Sig)
][
	order(match(cell_id, rownames(seurat_integrated@meta.data)))
]

seurat_integrated[["tdT_Log2_Ratio"]] <- plot_results[["tdT_Log2_Ratio"]]
seurat_integrated[["tdT_Ratio_FDR"]] <- plot_results[["tdT_Ratio_FDR"]]
seurat_integrated[["tdT_Ratio_Sig"]] <- plot_results[["tdT_Ratio_Sig"]]

## Plot Results.

p <- FeaturePlot(
	seurat_integrated, split.by = "orig.ident", pt.size = 0.01,
	features = c("tdT_Log2_Ratio", "tdT_Ratio_FDR", "tdT_Ratio_Sig")
)

pdf(file.path("results", "tdTomato", "tdT_Ratio.pdf"), height = 12, width = 16)
p
dev.off()

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
