library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)

data = Read10X_h5(filename = "data/5k_human_melanoma_tumor_CNIK_5pv2_filtered_feature_bc_matrix.h5")
counts = data

melanoma_seurat = CreateSeuratObject(
  counts = counts,
  project = "melanoma",
  min.cells = 3,
  min.features = 200
)

melanoma_seurat[["percent.mt"]] = PercentageFeatureSet(
  melanoma_seurat,
  pattern = "^MT-"
)

VlnPlot(
  melanoma_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)