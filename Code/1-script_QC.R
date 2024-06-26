# -------------------------------------------------------
#     
#             ZSCAN 4 - ANALYSIS (I)
# -------------------------------------------------------
## ------------------------------------------------------
#                   First script of the Analysis

# Libraries
library(Seurat)
library(dplyr)
library(patchwork) # design and combine visualizations 
library(Matrix)
library(ggplot2)
library(umap)
library(Rtsne)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(Biostrings)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(data.table)
library(SeuratObject)
library(psychTools)
library(data.table)
library(RColorBrewer)

# -----------------------------
#   1. Create Seurat Object 
# -----------------------------
sce.data <- Read10X(data.dir = "YourPath") # feature_bc_matrix
colnames(sce.data) <- gsub("-1", "", colnames(sce.data))


sce <- CreateSeuratObject(counts = sce.data, project = "Single-cells", min.cells = 3, min.features = 200)
dim(sce)
#  23338 18259
# Inspect the data
# head(sce.data[1:3, 1:3])


# -----------------------------
#   2. QUALITY CONTROL
# -----------------------------
# delete the column of the origin

sce@meta.data$orig.ident <- NULL

## Genes
genes <- rownames(sce@assays$RNA)
genes_df <- data.frame(Genes = genes)
View(genes_df)

# 1. mitochondrial genes
# ........................................................................................
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")

# 2. Filter the samples depending on the nFeatures and percentage of mit genes
##........................................................................................
View(sce@meta.data)
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, col = "mediumorchid1")

counts <- as.numeric(colSums(GetAssayData(sce)))
hist(counts, breaks = 50, col = "mediumorchid1", xlab = "Counts", ylab = "Frequency", main = "Histogram of counts before filtering")
# plots for the filteringz
p1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method='lm')
p1+p2


sce_filtered <- subset(sce, subset = nFeature_RNA > 1500 & percent.mt < 10 & nCount_RNA > 3000)

VlnPlot(sce_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, col = "mediumorchid1")

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "skyblue3")
plot1 + geom_hline(yintercept = max(sce_filtered@meta.data$percent.mt), linetype = "dashed", color = "red") +
  geom_vline(xintercept = max(sce_filtered@meta.data$nFeature_RNA), linetype = "dashed", color = "red") 


counts_filtered <- as.numeric(colSums(GetAssayData(sce_filtered)))
hist(counts_filtered, breaks = 50, col = "mediumorchid1", xlab = "Counts", ylab = "Frequency", main = "Histogram of counts after filtering")

n_cells <- ncol(sce) # Number of cells before filtering
n_filt_cells <- ncol(sce_filtered) # Number of resulted cells after filtering
perc_filtered <- (n_filt_cells/n_cells)*100
cat("\n", "number of cells before filtering =", n_cells, "\n", 
    "number of cells after filtering = ", n_filt_cells, "\n",
    "Percentage of filtered cells:", perc_filtered, "%\n",
    "number of genes:", dim(sce_filtered))

# number of cells before filtering = 18259 
# number of cells after filtering =  11662 
# Percentage of filtered cells: 63.86987 %

# number of cells before filtering = 18259 
# number of cells after filtering =  11413 
# Percentage of filtered cells: 62.50616 %


# colnames(sce@assays$RNA$counts[, "GFP"] == 0)
# View(sce@assays$RNA$counts[["GFP"]])ç


# MAYBE WE CAN REMOVED THE MITOCONDRIAL GENES AS IN THE PAPER¿?
# Obtener los índices de las características que no comienzan con "mt-"
sce <- subset(sce, features = !grepl("^mt-"))


# 3. Normalization and Scaling
##........................................................................................
sce_norm <- sce_filtered
sce_norm <- NormalizeData(sce_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
# the normalized data of the counts for each gene is saved in sce_norm@assays$RNA$data
sce_norm <- FindVariableFeatures(sce_norm, selection.method = "vst", nfeatures = 2000)
# head(sce_norm@assays$RNA$data)
top10 <- head(VariableFeatures(sce_norm), 10)

sce_norm <- ScaleData(sce_norm, features = genes)
sce_norm <- RunPCA(sce_norm, features = VariableFeatures(object = sce_norm))


# 4. Clustering
##........................................................................................
ElbowPlot(sce_norm) # ndims = 50)
sce_norm <- FindNeighbors(sce_norm, dims= 1:15)
sce_norm <- FindClusters(sce_norm, resolution = c(0.05), algorithm = 2)

DimPlot(sce_norm, group.by = "RNA_snn_res.0.05", label = TRUE,  cols = "Set2")

Idents(sce_norm) <- "RNA_snn_res.0.05"

sce_norm <- RunUMAP(sce_norm, dims = 1:15)
sce_norm <- RunTSNE(sce_norm, dims = 1:15)
DimPlot(sce_norm, reduction = "umap")
DimPlot(sce_norm, reduction = "tsne") #,  cols = "Accent")

sce_norm@meta.data$TSNE_cluster <- sce_norm$RNA_snn_res.0.05
# View(sce_norm@meta.data)

# cells in cluster 0 and in cluster 1
cat(" Cluster 0:", sum(sce_norm@meta.data$TSNE_cluster == 0), "\n",
    "Cluster 1:", sum(sce_norm@meta.data$TSNE_cluster == 1), "\n",
    "Cluster 2:", sum(sce_norm@meta.data$TSNE_cluster == 2), "\n")

# Cluster 0: 9841 
# Cluster 1: 1525 
# Cluster 2: 47 


# ------------------------
# SAVE THE SEURAT OBJECT
# ------------------------
saveRDS(sce_norm, file = "YourPath") # to save


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     FUNCTION PARA EXTRAER EL FEATURE PLOT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_metadata_feature <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "tsne"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  
  max_feature <- max(feature_value)
  
  ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = so@meta.data[[feature_name]])) +
    geom_point() +
    scale_color_gradient(low = "gray", high = "red", limits = c(0, max_feature)) +  # Ajustar los límites y la paleta de colores según lo deseado
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +  # Etiquetas de los ejes y la leyenda
    theme_minimal()
  
}
# FOR UMAP
plot_metadata_feature_umap <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "umap"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  
  max_feature <- max(feature_value)
  
  ggplot(plot_data, aes(x = umap_1, y = umap_2, color = so@meta.data[[feature_name]])) +
    geom_point() +
    scale_color_gradient(low = "gray", high = "red", limits = c(0, max_feature)) +  # Ajustar los límites y la paleta de colores según lo deseado
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +  # Etiquetas de los ejes y la leyenda
    theme_minimal()
  
}

plot_metadata_feature(sce_norm, "nFeature_RNA")


# let's do this but with GFP values
genes_gfp <- genes[grepl("^GFP", genes)]
genes_values <- GetAssayData(sce_norm, slot = "data")[genes_gfp,]

sce_norm@meta.data$gfp <- genes_values
plot_metadata_feature(sce_norm, "gfp")


# ...............................................................................................
#       SUBSET OF THE CELLS CLUSTER 1 GFP POSITIVES AND CLUSTER 0 AND CLUSTER 2 GFP NEGATIVES
# ...............................................................................................
sce_filt_1 <- sce_norm
sce_filt_1 <- subset(sce_norm, (TSNE_cluster == 1 & gfp == 0) | (TSNE_cluster == 0 & gfp > 0) | (TSNE_cluster == 2 & gfp > 0))

dim(sce_norm) # 23338 11413
dim(sce_filt_1) #  23338  4956

View(sce_filt_1@meta.data)
plot_metadata_feature(sce_filt_1, "gfp") #            with TSNE
plot_metadata_feature_umap(sce_filt_1, "gfp") #       WITH UMAP


cat(" Cluster 0:", sum(sce_filt_1@meta.data$TSNE_cluster == 0), "\n",
    "Cluster 1:", sum(sce_filt_1@meta.data$TSNE_cluster == 1), "\n",
    "Cluster 2:", sum(sce_filt_1@meta.data$TSNE_cluster == 2), "\n")

# Cluster 0: 3506 
# Cluster 1: 1410 
# Cluster 2: 40 


###### WE RUN AGAIN PCA TO SEE:

sce_filt_1 <- FindVariableFeatures(sce_filt_1, selection.method = "vst", nfeatures = 2000)
sce_filt_1 <- RunPCA(sce_filt_1, features = VariableFeatures(object = sce_filt_1))

ElbowPlot(sce_filt_1) # ndims = 50)
sce_filt_1 <- FindNeighbors(sce_filt_1, dims= 1:15)
sce_filt_1 <- FindClusters(sce_filt_1, resolution = 0.03, algorithm = 2)

DimPlot(sce_filt_1, group.by = "RNA_snn_res.0.03", label = TRUE,  cols = "Set2")

Idents(sce_filt_1) <- "RNA_snn_res.0.03"

sce_filt_1 <- RunUMAP(sce_filt_1, dims = 1:15)
sce_filt_1 <- RunTSNE(sce_filt_1, dims = 1:15)
DimPlot(sce_filt_1, reduction = "umap")
DimPlot(sce_filt_1, reduction = "tsne")


library(RColorBrewer)
palete <-brewer.pal(n = 3, name = "Set2")

DimPlot(sce_filt_1, reduction = "tsne", cols = palete) #,  cols = "Accent")

cat(" Cluster 0:", sum(sce_filt_1@meta.data$TSNE_cluster == 0), "\n",
    "Cluster 1:", sum(sce_filt_1@meta.data$TSNE_cluster == 1), "\n",
    "Cluster 2:", sum(sce_filt_1@meta.data$TSNE_cluster == 2), "\n")

# Cluster 0: 3506 
# Cluster 1: 1410 
# Cluster 2: 40 

# 2) SAVE THE OBJECT
#  .....................
saveRDS(sce_filt_1, file = "YourPath")


