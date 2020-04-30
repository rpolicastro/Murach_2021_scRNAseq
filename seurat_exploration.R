
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

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_obj.RDS"))

## Exploring tdTomato
## ----------

if (!dir.exists(file.path("results", "tdTomato"))) {
	dir.create(file.path("results", "tdTomato"))
}

## Plotting the expression of tdTomato and tdTomato 5' UTR.

p <- FeaturePlot(
	seurat_obj, feature = c("tdTomato", "tdTomato-fiveprime", "Pax7"),
	pt.size = 0.1, ncol = 2
)

pdf(file.path("results", "tdTomato", "tdTomato_dimplot.pdf"), height = 6, width = 8)
p
dev.off()

## Creating a ratio of tdTomato 5' UTR to CDS.

raw_counts <- seurat_obj@assays$RNA[
	rownames(seurat_obj@assays$RNA) %in% c("tdTomato", "tdTomato-fiveprime", "Pax7")
]

seurat_obj$tdTCDS_tdT5UTR_ratio <- apply(raw_counts, 2, function(x) {
	ratio <- log2(x["tdTomato"] + 1) - log2(x["tdTomato-fiveprime"] + 1)
	return(ratio)
})

p <- FeaturePlot(seurat_obj, feature = "tdTCDS_tdT5UTR_ratio", pt.size = 0.1) +
	scale_color_viridis_c()

pdf(file.path("results", "tdTomato", "tdTomato_ratio.pdf"), height = 5, width = 7.5)
p
dev.off()
