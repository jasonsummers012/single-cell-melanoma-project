library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)

#Read in the file
data = Read10X_h5(filename = "data/5k_human_melanoma_tumor_CNIK_5pv2_filtered_feature_bc_matrix.h5")
counts = data

#Create the Seurat object
melanoma_seurat = CreateSeuratObject(
  counts = counts,
  project = "melanoma",
  min.cells = 3,
  min.features = 200
)

#Create a column for mitochondrial RNA percentage
melanoma_seurat[["percent.mt"]] = PercentageFeatureSet(
  melanoma_seurat,
  pattern = "^MT-"
)

#Filter for more than 200 genes, less than 500 genes, and less than 5% mitochondrial RNA
melanoma_seurat = subset(
  melanoma_seurat,
  subset = nFeature_RNA > 200 &
  nFeature_RNA < 2500 &
  percent.mt < 5
)

#Plot the resulting object
VlnPlot(
  melanoma_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

#Save the results
saveRDS(melanoma_seurat, file = "results/melanoma_seruat_filtered.rds")