
library("tidyverse")
library("Seurat")
library("SeuratWrappers")
library("SeuratDisk")

##############################
## Export Seurat for scanpy ##
##############################

## Load seurat object.

seurat_obj <- readRDS(file.path("results", "r_objects", "seurat_integrated_spliced.RDS"))

seurat_obj <- list(
  tdT_Parental = subset(seurat_obj[[1]], subset = orig.ident == "tdT_Parental"),
  KY_Mononuclear = subset(seurat_obj[[1]], subset = orig.ident == "KY_Mononuclear"),
  Pax7_tdT_4day = subset(seurat_obj[[2]], subset = orig.ident == "Pax7_tdT_4day"),
  Pax7_DTA_4day = subset(seurat_obj[[2]], subset = orig.ident == "Pax7_DTA_4day")
)

## Save as H5Seurat.

if (!dir.exists(file.path("results", "py_objects"))) {
  dir.create(file.path("results", "py_objects"))
}

iwalk(seurat_obj, function(x, y) {
  x[["RNA"]] <- x[["spliced"]]
  DefaultAssay(x) <- "RNA"
  SaveH5Seurat(x, file.path("results", "py_objects", str_c(y, ".h5seurat")))
})

## Convert to h5ad.

walk(names(seurat_obj), function(x) {
  file_name <- file.path("results", "py_objects", str_c(x, ".h5seurat"))
  Convert(file_name, dest = "h5ad")
})
