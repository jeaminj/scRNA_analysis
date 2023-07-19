# Script performs normalization of filtered scRNAseq data
# input: filtered and merged seurat object, output: normalized seurat obj (and integrated if needed), ready for clustering 

# scRNAseq Analysis: 5k A549 Lung Carcinoma Cells, Treated (resveratrol) & untreated - Initial Quality Control
# data subsetted from source: https://www.10xgenomics.com/resources/datasets/30-k-a-549-lung-carcinoma-cells-treatments-transduced-with-a-crispr-pool-multiplexed-6-cm-os-3-1-standard-6-0-0

library(Seurat)
library(tidyverse)

# Normalization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filtered_seurat <- NormalizeData(object = filtered_seurat)

# Creates a split filtered_seurat object in case integration in necessary:
filtered_seurat_split <- SplitObject(filtered_seurat, split.by = "sample")

# Continuing from previous step (No integration):
filtered_seurat <- FindVariableFeatures(filtered_seurat) # identifies highly variable features (default nfeatures=2000)
filtered_seurat <- ScaleData(filtered_seurat) # Scaling
filtered_seurat <- RunPCA(filtered_seurat) # linear dimensionaility reduction by principal component analysis

# Visualizing PCA results:
DimPlot(filtered_seurat, reduction = 'pca', group.by = 'sample') # condition specific clustering is seen -> need to integrate

# Integration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# using the filtered_seurat_split object created earlier:
for (i in 1:length(filtered_seurat_split)) {
    filtered_seurat_split[[i]] <- SCTransform(filtered_seurat_split[[i]], vst.flavor = "v2") # SCTransform to normalize and stabilize variance
}

# Selecting for most variable features to use for integration (default nfeatures=2000)
integ_features <- SelectIntegrationFeatures(object.list = filtered_seurat_split)

filtered_seurat_split <- PrepSCTIntegration(object.list = filtered_seurat_split,
                                           anchor.features = integ_features)
# finds integration anchors 
integ_anchors <- FindIntegrationAnchors(object.list = filtered_seurat_split,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
# Integrating across conditions
integrated_seurat <- IntegrateData(anchorset = integ_anchors, 
                                  normalization.method = "SCT")

# ~~~ quick check to see if integration worked ~~~
integrated_seurat <- RunPCA(integrated_seurat)
DimPlot(integrated_seurat, reduction = 'pca', group.by = 'sample') # condition specific clustering is no longer seen


