library(Seurat)
library(Matrix)
library(tidyverse)
library(ggplot2)

#Load the Seurat object
melanoma_seurat = readRDS("results/melanoma_seurat_filtered.rds")

#Normalize data
melanoma_seurat = NormalizeData(
  melanoma_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

#Identify the 2000 most variable features
melanoma_seurat = FindVariableFeatures(
  melanoma_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

#Identify the top 10 most variable features
top10 = head(VariableFeatures(melanoma_seurat), 10)

#Plot the most variable features and label the top 10
feature_plot = VariableFeaturePlot(melanoma_seurat)
LabelPoints(
  plot = feature_plot,
  points = top10,
  repel = TRUE,
  max.overlaps = 20
)

