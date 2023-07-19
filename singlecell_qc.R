# Script performs initial QC of raw data
# input: raw 10x matrix data, output: filtered and merged seurat object, ready for normalization/integration 

# scRNAseq Analysis: 5k A549 Lung Carcinoma Cells, Treated (resveratrol) & untreated - Initial Quality Control
# data subsetted from source: https://www.10xgenomics.com/resources/datasets/30-k-a-549-lung-carcinoma-cells-treatments-transduced-with-a-crispr-pool-multiplexed-6-cm-os-3-1-standard-6-0-0

# loading libraries
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)

# Pre QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Reads 10x data and creates Seurat object for each sample group (control & treated)
for (file in c("ctrl_feature_bc_matrix", "resveratrol_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, #Extracting just the Gene Expression datatype
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj)
}

# Merging the two Seurat objects into one
merged_seurat <- merge(x = ctrl_feature_bc_matrix,
                       y = resveratrol_feature_bc_matrix,
                       add.cell.id = c("control", "treated"))

# QC metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. log10GenesPerUMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# 2. MT reads (ratio)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio/100 #converts the percentage value to a ratio

# Adding additional column info to metadata:

metadata <- merged_seurat@meta.data #Creating a separate metadata dataframe to work with so as not to risk affecting the data in the merged_seurat object

# Creating sample ID column
metadata$sample <- NA
metadata$sample[which(str_detect(rownames(metadata), "^control_"))] <- "control"
metadata$sample[which(str_detect(rownames(metadata), "^treated_"))] <- "treated"

# Renaming columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Adding the metadata dataframe back to merged_seurat object
merged_seurat@meta.data <- metadata

# Visualizing QC Metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Number of UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Distribution of genes detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10000) +
  geom_hline(yintercept = 3000) +
  facet_wrap(~sample)

# Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subsetting raw data based on thresholds from QC
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 10000) & 
                                  (nGene >= 3000) & 
                                  (log10GenesPerUMI > 0.80) & 
                                  (mitoRatio < 0.20))

# Visualizing filtered data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
metadata_clean <- filtered_seurat@meta.data

# Number of cell counts per sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Number of UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 10000)

# Distribution of genes detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 3000)

# Overall complexity of gene expression by visualizing the genes detected per UMI (novelty score)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Distribution of mitochondrial gene expression detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 10000) +
  geom_hline(yintercept = 3000) +
  facet_wrap(~sample)

# 
