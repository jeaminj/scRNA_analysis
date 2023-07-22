# Script performs clustering of scRNAseq normalized/integrated data
# input: normalized/integrated seurat object, output: 

# scRNAseq Analysis: 5k A549 Lung Carcinoma Cells, Treated (resveratrol) & untreated - Initial Quality Control
# data subsetted from source: https://www.10xgenomics.com/resources/datasets/30-k-a-549-lung-carcinoma-cells-treatments-transduced-with-a-crispr-pool-multiplexed-6-cm-os-3-1-standard-6-0-0

library(Seurat)
library(SingleR)
library(celldex)
library(pheatmap)

# Clustering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determines dimensionality of the data (number of clusters that should be chosen)
ElbowPlot(integrated_seurat) 

# 
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))

# check cluster resolutions
View(integrated_seurat@meta.data)
DimPlot(integrated_seurat, group.by = "integrated_snn_res.0.2", label = TRUE)
DimPlot(integrated_seurat, group.by = "integrated_snn_res.0.4", label = TRUE)
DimPlot(integrated_seurat, group.by = "integrated_snn_res.0.6", label = TRUE)
DimPlot(integrated_seurat, group.by = "integrated_snn_res.0.8", label = TRUE)
DimPlot(integrated_seurat, group.by = "integrated_snn_res.1", label = TRUE)

# assigning identity of clusters
Idents(integrated_seurat) <- "seurat_clusters"

# non linear dimensionality reduction
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:10)

DimPlot(integrated_seurat, reduction = 'umap', label = TRUE)
DimPlot(integrated_seurat, reduction = 'umap', group.by = 'sample')

# Auto-annotation with SingleR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# reference data 
ref <- celldex::HumanPrimaryCellAtlasData()

# get assay data
integrated_seurat_data <- GetAssayData(integrated_seurat, slot = 'data')

pred <- SingleR(test = integrated_seurat_data,
                ref = ref,
                labels = ref$label.main)

# Create and add singleR annotations column to the integrated_seurat metadata from predicted annotations 
integrated_seurat$singleR_annotation <- pred$labels[match(rownames(integrated_seurat@meta.data), rownames(pred))]

# Visualize annotated plot
DimPlot(integrated_seurat, reduction = 'umap', group.by = 'singleR_annotation', label =TRUE) + NoLegend()


# Annotation Quality Check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotScoreHeatmap(pred)
plotDeltaDistribution(pred)

# Comparing to unsupervised clustering
tab <- table(Assigned=pred$labels, Clusters=integrated_seurat$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','brown'))(10))

