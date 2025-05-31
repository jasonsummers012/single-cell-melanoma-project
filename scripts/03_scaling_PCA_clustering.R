library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(clustree)

#Load the Seurat object
melanoma_seurat = readRDS("results/melanoma_seurat_normalized.rds")

#Scale the data (only variable features by default)
melanoma_seurat = ScaleData(melanoma_seurat)

#Perform dimensionality reduction - PCA
melanoma_seurat = RunPCA(
  melanoma_seurat,
  features = VariableFeatures(
    object = melanoma_seurat
  )
)

#Save the top contributing genes to PCA
sink("results/pca_top_genes.txt")
print(melanoma_seurat[["pca"]],
      dims = 1:5,
      nfeatures = 5
)
sink()

#Visualize PCA results
DimHeatmap(
  melanoma_seurat,
  dims = 1,
  cells = 500,
  balanced = TRUE
)

ElbowPlot(melanoma_seurat)

#Perform clustering
melanoma_seurat = FindNeighbors(
  melanoma_seurat,
  dims = 1:15
)

melanoma_seurat = FindClusters(
  melanoma_seurat,
  resolution = c(0.1, 0.3, 0.5, 0.7, 1)
)

#Compares different resolutions to pick the best one
clustree(melanoma_seurat)
ggsave("results/clustree_plot.png")

DimPlot(
  melanoma_seurat,
  group.by = "RNA_snn_res.0.5",
  label = TRUE
)

Idents(melanoma_seurat) = "RNA_snn_res.0.5"
Idents(melanoma_seurat)

#Non-linear dimensionality reduction
melanoma_seurat = RunUMAP(
  melanoma_seurat,
  dims = 1:15
)

DimPlot(melanoma_seurat,
        reduction = "umap"
)

#Identify marker genes for each cluster
markers = FindAllMarkers(
  melanoma_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers = markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

#Save the results of top markers
sink("results/top_markers.txt")
print(top_markers)
sink()