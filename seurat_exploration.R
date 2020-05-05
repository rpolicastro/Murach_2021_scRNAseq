
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

##################################
## Exploring the scRNA-seq Data ##
##################################

## Loading the Data
## ----------

seurat_integrated <- readRDS(file.path("results", "r_objects", "seurat_integrated.RDS"))

## Exploring tdTomato
## ----------

if (!dir.exists(file.path("results", "tdTomato"))) {
	dir.create(file.path("results", "tdTomato"))
}

## Plotting the expression of tdTomato and tdTomato 5' UTR.

p <- FeaturePlot(
	seurat_integrated, feature = c("tdTomato", "tdTomatoStop", "Pax7"),
	pt.size = 0.1, ncol = 2, split.by = "ident"
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

raw_counts[, tdT_ratio := (tdTomato + 1) / (tdTomatoStop + 1)]

## Grabbing cell meta data.

meta_data <- as.data.table(seurat_integrated@meta.data, keep.rownames = "cell_id")[,
	.(orig.ident, cell_id)
]

## Merge the meta.data into the cell data.

merged <- merge(meta_data, raw_counts, by = "cell_id")

## Preparing the negative control and test set.

neg_control <- merged[orig.ident == "tdT_Parental"][, .(orig.ident, cell_id, tdT_ratio)]
experimental <- merged[orig.ident == "KY_mononuclear"][, .(orig.ident, cell_id, tdT_ratio)]

## Permutation test.

permute_rank <- function(combined_data) {
	permuted_rank <- sample_frac(combined_data, 1L)[1, "rank"]
	return(permuted_rank)
}

permutation_results <- map_df(experimental[["cell_id"]], function(x) {

	combined <- rbind(neg_control, row_args)
	combined[, rank := dense_rank(-tdT_ratio)]

	actual_rank <- combined[cell_id == row_args[1, "cell_id"]][["rank"]]
	perm_result <- replicate(10, permute_rank(combined))
	perm_result <- sum(perm_result <= actual_rank)

	row_args[["pval"]] <- (perm_result + 1) / (10 + 1)
	return(row_args)
})
