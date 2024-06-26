
## -------------------------------------------------------
#     
#             ZSCAN 4 - ANALYSIS  (III)
# -------------------------------------------------------
## ------------------------------------------------------
######################################################################################################

#         Analysis of the SeuratObject with Cluster control (1), pluri (0) and toti (2)

#######################################################################################################


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
library(ComplexHeatmap)

# Load the seurat object
sce_zscan4_reads <- readRDS("YourPath")

# =====================
#     ZSCAN4 GENES
# =====================
all.genes<-rownames(sce_zscan4_reads)
genes_df <- data.frame(Gene = all.genes)

genes_zscan <- all.genes[grepl("^Zscan4", all.genes)]
genes_zscan <- genes_zscan[!grepl("-ps[123]", genes_zscan)]


FeaturePlot(sce_zscan4_reads, features = genes_zscan, reduction = "tsne")

zscan4_expression <- GetAssayData(object = sce_zscan4_reads, layer = "data")[genes_zscan, ]  
valores_no_cero <- zscan4_expression[zscan4_expression != 0]


frame <- as.data.frame(t(zscan4_expression))
frame$mean <- rowMeans(frame)
# add it to the metadata
sce_zscan4_reads@meta.data$zscan4_expr <- frame$mean 
zscan4_expression <- sce_zscan4_reads@meta.data$zscan4_expr
sce_zscan4_reads@meta.data$zscan4_expr <- zscan4_expression 


# SPLIT THE SCE WITHOUT CLUSTER 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sce_cl02 <- subset(sce_zscan4_reads, TSNE_cluster != 1)
sce_cl02 <- sce_zscan4_reads


zd <- all.genes[grepl("^Zscan4d", all.genes)]
zscan4d_expression <- GetAssayData(object = sce_cl02, slot = "data")[zd, ] 
sce_cl02@meta.data$zscan4d<- zscan4d_expression 

zc <-all.genes[grepl("^Zscan4c", all.genes)]
zscan4c_expression <- GetAssayData(object = sce_cl02, slot = "data")[zc, ] 
sce_cl02@meta.data$zscan4c<- zscan4c_expression 

zf <-all.genes[grepl("^Zscan4f", all.genes)]
zscan4f_expression <- GetAssayData(object = sce_cl02, slot = "data")[zf, ] 
sce_cl02@meta.data$zscan4f<- zscan4f_expression 

genes_zscan <- c("Zscan4d", "Zscan4c", "Zscan4f", "Zscan4b", "Zscan4e")
for (gene_name in genes_zscan) {
  gene_ids <- grep(paste0("^", gene_name), all.genes, value = TRUE)
  gene_expression <- GetAssayData(object = sce_cl02, slot = "data")[gene_ids, ]
  sce_cl02@meta.data[[tolower(gene_name)]] <- gene_expression
}

# stat3
Stat3 <- all.genes[grepl("^Stat3", all.genes)]
Stat3_expression <- GetAssayData(object = sce_zscan4_reads, slot = "data")[Stat3, ] 
sce_zscan4_reads@meta.data$Stat3 <- Stat3_expression 
plot_metadata_feature_tsne(sce_zscan4_reads, "Zdbf2")



sce_cl00 <- subset(sce_cl02, TSNE_cluster != 2)
sce_z4c <- subset(sce_cl02, zscan4c != 0)
dim(sce_z4c)
correlation_result <- cor.test(sce_z4c@meta.data$zscan4c, sce_z4c@meta.data$Diversity, method = "spearman")
print(paste("Spearman's rho:", correlation_result$estimate))
print(paste("P-value:", correlation_result$p.value)) 


# =======================
#       1)  Plots --> TSNE
# =======================

# we delete the cells in cluster 1 with zscan4 expression
sce_zscan4_reads <- subset(sce_zscan4_reads, !(sce_zscan4_reads@meta.data$TSNE_cluster == 1 & zscan4_expr > 0))

DimPlot(sce_zscan4_reads, reduction = "umap", pt.size = 2, cols = "Accent")

#   1.1) Plot functions to represent a feature of the metadata in tsne or umap reduction
library(viridis)
plot_metadata_feature_tsne <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "tsne"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  max_feature <- max(feature_value)
  ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = so@meta.data[[feature_name]])) +
    geom_point(size = 3) +
    scale_color_gradient(low = "gray", high = "#CD3333", limits = c(0, max_feature)) +  
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 19),  # Ajusta el tamaño del título del eje X
      axis.title.y = element_text(size = 19),  # Ajusta el tamaño del título del eje Y
      legend.text = element_text(size = 13)) #+
  #guides(color = guide_legend(title = feature_name, title.theme = element_text(size = 19)))
}


plot_metadata_feature_tsne <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "tsne"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  max_feature <- max(feature_value)
  
  ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = so@meta.data[[feature_name]])) +
    geom_point(size = 3) +
    scale_color_viridis(option = "B", direction = -1, limits = c(0, max_feature)) +  
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 19),  # Ajusta el tamaño del título del eje X
      axis.title.y = element_text(size = 19),  # Ajusta el tamaño del título del eje Y
      legend.text = element_text(size = 13)) #+
  #guides(color = guide_legend(title = feature_name, title.theme = element_text(size = 19)))
}

# SAVE
##########
ggsave("YourPath", plot = p, device = "pdf", width = 8, height = 7)
#################
plot_metadata_feature_tsne(sce_zscan4_reads, "Diversity")
p<-plot_metadata_feature_tsne(sce_zscan4_reads, "Diversity")
stat <- plot_metadata_feature_tsne(sce_zscan4_reads, "Stat3")
plot_metadata_feature_tsne(sce_zscan4_reads, "Coverage")
plot_metadata_feature_tsne(sce_zscan4_reads, "Uncuts")



plot_gfp <- plot_metadata_feature_tsne(sce_zscan4_reads, "Cas9")
plot_metadata_feature_tsne(sce_DEA, "zscan4_expr")
plot_zscan <- plot_metadata_feature_tsne(sce_zscan4_reads, "zscan4_expr")
stat <- plot_metadata_feature_tsne(sce_zscan4_reads, "Stat3")



plot_metadata_feature_umap <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "umap"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  max_feature <- max(feature_value)
  ggplot(plot_data, aes(x = umap_1, y = umap_2, color = so@meta.data[[feature_name]])) +
    geom_point() +
    scale_color_gradient(low = "gray", high = "red", limits = c(0, max_feature)) + 
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +  
    theme_minimal()
}

plot_metadata_feature_umap(sce_zscan4_reads, "Diversity")
plot_metadata_feature_umap(sce_zscan4_reads, "Coverage")
plot_metadata_feature_umap(sce_zscan4_reads, "Uncuts")
plot_metadata_feature_umap(sce_zscan4_reads, "gfp")
plot_metadata_feature_umap(sce_zscan4_reads, "zscan4_expr")

plot_metadata_feature_umap(sce_cl02, "Diversity")

sce_cl02_2 <- subset(sce_cl02, !(sce_cl02@meta.data$TSNE_cluster == 2 & zscan4_expr == 0))
plot_metadata_feature_umap(sce_cl02_2, "Diversity")
plot_metadata_feature_umap(sce_cl02_2, "zscan4_expr")


#   1.2) Boxplots of the features
# --------------------------------------

# To change the order of levels in the metadata
sce_zscan4_reads@meta.data[["TSNE_cluster"]] <- factor(sce_zscan4_reads@meta.data[["TSNE_cluster"]], levels = c("1", "2", "0"))# c("1", "2", "0")
sce_zscan4_reads@meta.data[["TSNE_cluster"]] # Levels: 1 2 0

# WITHOUT OUTLIERS    

boxplot_metadata_feature <- function(seurat_obj, feature_name) {
  metadata <- seurat_obj@meta.data
  
  plot_data <- data.frame(
    Feature = metadata[[feature_name]],
    TSNE_cluster = factor(metadata$TSNE_cluster)
  )
  levels(plot_data$TSNE_cluster) <- c("Control", "Totipotents\nZscan4+", "Pluripotents\nZscan4-")
  
  # Calcular los cuantiles y el rango intercuartílico para cada cluster
  for (cluster in unique(plot_data$TSNE_cluster)) {
    subset_data <- plot_data[plot_data$TSNE_cluster == cluster, ]
    Q1 <- quantile(subset_data$Feature, 0.25, na.rm = TRUE)
    Q3 <- quantile(subset_data$Feature, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    
    # Calcular los límites para identificar los valores atípicos
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    # Filtrar los datos para excluir los valores atípicos
    plot_data <- plot_data[!(plot_data$TSNE_cluster == cluster & 
                               (plot_data$Feature < lower_bound | plot_data$Feature > upper_bound)), ]
  }
  
  # Access specific colors from Accent palette
  accent_colors <- brewer.pal(8, "Accent")
  colors <- c(accent_colors[2], accent_colors[3], accent_colors[1])
  
  gg <- ggplot(plot_data, aes(x = TSNE_cluster, y = Feature, fill = TSNE_cluster)) +
    geom_boxplot(alpha = 0.5, width = 0.5, aes(color = TSNE_cluster), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.25, aes(color = TSNE_cluster), size = 1) + # Más grueso
    xlab("") +
    ylab("Zscan4 expression") + # "%Original sequence" "Diversity per cell" Mean length (bp)
    # labs(title = paste0(feature_name, " distribution in each cluster")) +
    theme_minimal() +
    geom_jitter(aes(color = TSNE_cluster), width = 0.2, alpha = 1) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    #ggtitle(feature_name) +
    theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"),
          axis.title = element_text(size = 19, color = "black"),
          axis.text = element_text(size = 19, color = "black"),
          panel.grid.major = element_line(color = "gray85", size = 0.5), 
          panel.grid.minor = element_line(color = "gray85", size = 0.25),
          legend.position = "none")
  
  return(gg)
}

b1 <- boxplot_metadata_feature(sce_zscan4_reads, "zscan4_expr")
b2<- boxplot_metadata_feature(sce_zscan4_reads, "Uncuts")
b3<-boxplot_metadata_feature(sce_zscan4_reads, "Diversity")
b4<-boxplot_metadata_feature(sce_zscan4_reads, "Mean_length")



# ===========================================
#     2) DIFFERENTIAL EXPRESSION CL0-CL2
# ============================================

# Differential expression analysis between the cluster 0 and cluster 3 (the whole cluster 0 and whole cluster 3)
# GOAL: Run FindMarkers() function between the whole cluster 0 and 3 (extracting the 12 zscan(+) cells)

sce_proces_2 <- sce_cl02
Idents(sce_proces_2)
markers_02 <- FindMarkers(sce_proces_2, ident.1 = "0", ident.2 = "2")

markers_02$p_val_adj = p.adjust(markers_02$p_val, method='fdr')
markers_02$gene <- rownames(markers_02) 

markers02<- subset(markers_02, p_val_adj < 0.05)

# ORDER 
markers02_ordered <- markers02 %>% 
  arrange(desc(avg_log2FC))
head(markers02_ordered)

#         AFTER
#                       p_val avg_log2FC pct.1 pct.2    p_val_adj          gene
# Gm11627       5.448506e-05   2.660516 0.431 0.000 1.001308e-03       Gm11627
# Filip1l       5.100793e-04   2.400613 0.350 0.000 7.021740e-03       Filip1l
# B230118H07Rik 6.854134e-06   2.191469 0.596 0.125 1.493936e-04 B230118H07Rik
# Ntn1          8.385522e-05   2.136834 0.558 0.125 1.479792e-03          Ntn1
# Zdbf2         6.307286e-09   2.131956 0.860 0.333 2.321176e-07         Zdbf2
# Gm28960       1.438176e-04   2.115025 0.507 0.125 2.289671e-03       Gm28960

#SAVE DATAFRAME ORDERED
write.csv(markers02_ordered, file = "YourPath")
markers02_ordered <- read.csv("YourPath")

#     2.2.) see what genes are in pluripotency
# ------------------------------------------> SUBSET DOWN OR UP (totigenes)
all.genes<-rownames(sce_proces_2)
genes_df <- data.frame(Gene = all.genes)
pluri_genes <- markers02_ordered$gene[markers02_ordered$avg_log2FC>1]
toti_genes <- markers02_ordered$gene[markers02_ordered$avg_log2FC < -2]
all <- markers02_ordered$gene[markers02_ordered$avg_log2FC < -2 | markers02_ordered$avg_log2FC>1]
markers_02_pluri <- subset(markers02_ordered, markers02_ordered$X  %in% pluri_genes)
markers_02_toti <- subset(markers02_ordered, markers02_ordered$X %in% toti_genes)

#toti
sorted_markers <- markers_02_toti[order(markers_02_toti$avg_log2FC), ]
top_200_genes <- head(sorted_markers$gene, 200)

#pluri
sorted_markers_plu <- markers_02_pluri[order(markers_02_pluri$avg_log2FC), ]
top_200_genes_plu <- head(sorted_markers_plu$gene, 200)

markers_02_all <- subset(markers02_ordered, markers02_ordered$X %in% all)

#     2.3.) PLOTS DEA
# -----------------------------------------
# - Barplot with number up/down DEG
# ______________________________________
# Step 1: Count the number of upregulated and downregulated genes
num_upregulated <- length(toti_genes)
num_downregulated <- length(pluri_genes)

# Step 2: Create a data frame with the counts
deg_counts <- data.frame(
  Category = c("Upregulated", "Downregulated"),
  Count = c(num_upregulated, num_downregulated)
)

# Step 3: Generate the barplot
barplot_DE_genes <- ggplot(deg_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"), # Título más oscuro y grande
        axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
        axis.text = element_text(size = 19, color = "black"),
        panel.grid.major = element_line(color = "gray85", size = 0.5),  # Adjusts major grid lines
        panel.grid.minor = element_line(color = "gray85", size = 0.25),# Texto de ejes más oscuro
        legend.position = "none")+
  labs(title = "",
       x = "",
       y = "Number of Differentially Expressed Genes",
       main =  "DEA Cluster 2 vs Cluster")+
  scale_fill_manual(values = c("Upregulated" = "#fa9fb5", "Downregulated" = "#bcbddc"))+
  geom_text(aes(label = Count), vjust = -0.2, size = 7)




# - Select top up-regulated and down-regulated DEGs and show a boxplot of their expression in cluster 0 and cluster 2.
# _______________________________________________________________________________________________________________________
top_100_toti <- head(sorted_markers$gene, 100)
top_100_plu <- head(sorted_markers_plu$gene, 100)

selected_genes <- unique(c(top_100_toti, top_100_plu))
library(tidyr)

# Extract expression data for the selected genes in clusters 0 and 2
cluster_0_cells <- WhichCells(sce_proces_2, ident = "0")
cluster_2_cells <- WhichCells(sce_proces_2, ident = "2")



expression_data_toti_cluster_0 <- FetchData(sce_proces_2, vars = top_100_toti, cells = cluster_0_cells)
expression_data_toti_cluster_2 <- FetchData(sce_proces_2, vars = top_100_toti, cells = cluster_2_cells)

expression_data_pluri_cluster_0 <- FetchData(sce_proces_2, vars = top_100_plu, cells = cluster_0_cells)
expression_data_pluri_cluster_2 <- FetchData(sce_proces_2, vars = top_100_plu, cells = cluster_2_cells)

expression_data_toti_cluster_0$mean_expression <- rowMeans(expression_data_toti_cluster_0[, top_100_toti])
expression_data_toti_cluster_2$mean_expression <- rowMeans(expression_data_toti_cluster_2[, top_100_toti])

expression_data_pluri_cluster_0$mean_expression <- rowMeans(expression_data_pluri_cluster_0[, top_100_plu])
expression_data_pluri_cluster_2$mean_expression <- rowMeans(expression_data_pluri_cluster_2[, top_100_plu])

expression_data_toti_cluster_0$Cell <- rownames(expression_data_toti_cluster_0)
expression_data_toti_cluster_2$Cell <- rownames(expression_data_toti_cluster_2)

expression_data_pluri_cluster_0$Cell <- rownames(expression_data_pluri_cluster_0)
expression_data_pluri_cluster_2$Cell <- rownames(expression_data_pluri_cluster_2)

# Combine the expression data into a single data frame for plotting
expression_data_toti_cluster_0$Cluster <- "Cluster 0"
expression_data_toti_cluster_2$Cluster <- "Cluster 2"

expression_data_pluri_cluster_0$Cluster <- "Cluster 0"
expression_data_pluri_cluster_2$Cluster <- "Cluster 2"

expression_data_toti <- rbind(expression_data_toti_cluster_0[, c("Cell", "mean_expression", "Cluster")], 
                              expression_data_toti_cluster_2[, c("Cell", "mean_expression", "Cluster")])

expression_data_pluri <- rbind(expression_data_pluri_cluster_0[, c("Cell", "mean_expression", "Cluster")], 
                               expression_data_pluri_cluster_2[, c("Cell", "mean_expression", "Cluster")])

# Reorder the factor levels for Cluster
expression_data_toti$Cluster <- factor(expression_data_toti$Cluster, levels = c("Cluster 2", "Cluster 0"))
expression_data_pluri$Cluster <- factor(expression_data_pluri$Cluster, levels = c("Cluster 2", "Cluster 0"))



# Create the boxplots
accent_colors <- brewer.pal(8, "Accent")
colors <- c( accent_colors[3], accent_colors[1])

boxplot_toti <- ggplot(expression_data_toti, aes(x = Cluster, y = mean_expression, fill = Cluster)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Cluster), outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Expression of Top 100 Up-regulated DEGs",
       x = "",
       y = "Mean Expression") +
  geom_jitter(aes(color = Cluster), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"), 
        axis.title = element_text(size = 19, color = "black"),
        axis.text = element_text(size = 19, color = "black"),
        panel.grid.major = element_line(color = "gray85", size = 0.5),  
        panel.grid.minor = element_line(color = "gray85", size = 0.25),
        legend.position = "none")

boxplot_toti <- ggplot(expression_data_toti, aes(x = Cluster, y = mean_expression, fill = Cluster)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Cluster), outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Expression of Top 100 Up-regulated DEGs",
       x = "",
       y = "Mean Expression") +
  geom_jitter(aes(color = Cluster), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 19, color = "black"), 
        axis.text = element_text(size = 19, color = "black"),
        panel.grid.major = element_line(color = "gray85", size = 0.5),  
        panel.grid.minor = element_line(color = "gray85", size = 0.25),
        legend.position = "none")

boxplot_toti

boxplot_pluri <- ggplot(expression_data_pluri, aes(x = Cluster, y = mean_expression, fill = Cluster)) +
  geom_boxplot(alpha = 0.5, width = 0.5, aes(color = Cluster), outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Expression of Top 100 Down-regulated DEGs",
       x = "",
       y = "Mean Expression") +
  geom_jitter(aes(color = Cluster), width = 0.2, alpha = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme(plot.title = element_text(hjust = 0.6, size = 16, color = "black", face = "bold"), # Título más oscuro y grande
        axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
        axis.text = element_text(size = 19, color = "black"),
        panel.grid.major = element_line(color = "gray85", size = 0.5),  # Adjusts major grid lines
        panel.grid.minor = element_line(color = "gray85", size = 0.25),# Texto de ejes más oscuro
        legend.position = "none")


# - Plot gene ontology enrichments and/or provide a list. Moreover, add a section in methods to know how it has been performed.
# ______________________________________________________________________________________________________________________________________



#     2.4.) GENE ONTOLOGY
# ------------------------------------------> GENE ONTOLOGY
library(GOstats)
library(org.Mm.eg.db)
library(clusterProfiler)

# PLURIGENES
pluri_genes
entrez_ids <- mapIds(org.Mm.eg.db, keys = pluri_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# GENES UNIVERSE
length(intersect(all.genes, all))
length(all)
all.genes_2 <- setdiff(all.genes, all )
length(all.genes)
length(all.genes_2)

entrez_ids <- mapIds(org.Mm.eg.db, keys = all, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")


#         Apply some filters in genes
#         -----------------------------------
sce_proces_2
counts_genes <- GetAssayData(object = sce_proces_2, layer = "counts")
length(rownames(counts_genes))
filtered_genes <- rownames(counts_genes)[rowSums(counts_genes >= 1)]
length(filtered_genes)

filtered_genes <- rownames(counts_genes)
filtered_genes <- sample(all.genes, 22000)

# subset count matrix and see genes present in more than 3 cells 
counts_filtered <- counts_genes[filtered_genes, ]
filtered_genes_cells <- rownames(counts_filtered)[rowSums(counts_filtered > 0) >= 3]
length(filtered_genes_cells)

all.genes_2 <- intersect(all.genes, rownames(counts_genes)[rowSums(counts_genes >= 2)])
length(all.genes_2)
#         UNIVERSE
#     -------------------
background_entrez_ids <- mapIds(org.Mm.eg.db, keys = all.genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
background_entrez_ids <- background_entrez_ids[!is.na(background_entrez_ids)]
background_entrez_ids <- unique(background_entrez_ids)
background_entrez_ids <- background_entrez_ids[background_entrez_ids != ""]

ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  universe = background_entrez_ids,
  readable = TRUE
)

length(intersect(background_entrez_ids, entrez_ids))
summary(as.data.frame(ego))
go_df_2 <- as.data.frame(summary(ego))
# Visualize the results
dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment Analysis")
fit <- plot(barplot(ego, showCategory = 20))

# TOTI GENES
toti_genes
entrez_ids_toti <- mapIds(org.Mm.eg.db, keys = toti_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
background_entrez_ids <- mapIds(org.Mm.eg.db, keys = all.genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")


ego_toti <- enrichGO(gene = entrez_ids_toti,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     #universe = background_entrez_ids,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)

summary(as.data.frame(ego_toti))
fit <- plot(barplot(ego_toti, showCategory = 20))








#--------------------------------------------------------------------------------------------------------------

# ======================================
#            GROUPS OF CELLS 
# ======================================

#### ---->>>> GROUPS IN SCE_CL02
toti <- rownames(sce_cl02@meta.data[sce_cl02@meta.data$TSNE_cluster == 2, ]) 
transient <- rownames(sce_cl02@meta.data[sce_cl02@meta.data$zscan4_expr > 0 & sce_cl02@meta.data$TSNE_cluster == 0, ]) 
pluri<- rownames(sce_cl02@meta.data[sce_cl02@meta.data$TSNE_cluster == 0 & sce_cl02@meta.data$zscan4_expr == 0, ]) 

cat(" number of TOTIPOTENT cells -->", length(toti), "\n",
    "number of TRANSIENT cells --> ", length(transient), "\n",
    "number of PLURIPOTENT cells --> ", length(pluri), "\n")

cells_to_keep <- c(toti, pluri, transient)
length(cells_to_keep) # 1584

sce_exp1<- sce_cl02[, cells_to_keep]
dim(sce_exp1) # 23338  2967
sce_exp1@meta.data$grupo <- ifelse(colnames(sce_exp1) %in% toti, "totipotents", # number for each group
                                   ifelse(colnames(sce_exp1) %in% transient, "transient", 
                                          ifelse(colnames(sce_exp1) %in% pluri,"pluripotents", NA)))

# number of TOTIPOTENT cells --> 24 
# number of TRANSIENT cells -->  15 
# number of PLURIPOTENT cells -->  1545 


my_colors <- brewer.pal(12, "Set3")[c(8, 3, 1,6)]  # Por ejemplo, colores 1, 4 y 7
# my_colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

DimPlot(sce_exp1, reduction = "tsne", group.by = "grupo", pt.size = 3, cols = my_colors)+
  ggtitle("Groups of cells")



# ===================================================================================
#     4) CL 0: differences between the expression of toti and pluri genes in cells
# ===================================================================================

#   3.1.) compute the expression of the toti and pluri genes in each cell inside sce_exp1
totipotent_expr <- FetchData(object = sce_exp1, vars = toti_genes)
pluripotent_expr <- FetchData(object = sce_exp1, vars = pluri_genes)

#   3.2.) Add it to the metadata

gene_expr_data <- data.frame(
  totipotent = rowMeans(totipotent_expr),
  pluripotent = rowMeans(pluripotent_expr),
  diversity = sce_exp1@meta.data[["Diversity"]],
  coverage = sce_exp1@meta.data[["Coverage"]],
  ncount = sce_exp1@meta.data[["nCount_RNA"]],
  tsne_cluster = sce_exp1@meta.data[["TSNE_cluster"]],
  zscan4 = sce_exp1@meta.data[["zscan4_expr"]],
  mit = sce_exp1@meta.data$percent.mt,
  z4d = sce_exp1@meta.data$zscan4d,
  z4c = sce_exp1@meta.data$zscan4c,
  z4f = sce_exp1@meta.data$zscan4f,
  uncuts = sce_exp1@meta.data[["Uncuts"]],
  group = sce_exp1@meta.data[["grupo"]]
)
dim(gene_expr_data) #  1584   13

# 3.3.) extract the cells with zscan4 expresion

rows_with_expression <- rownames(sce_exp1@meta.data[sce_exp1@meta.data[["zscan4_expr"]] > 0, ])
gene_expr_data_zscan <- subset(gene_expr_data, rownames(gene_expr_data) %in% rows_with_expression)
# 32 cells with zscan4 expression
# palete <-brewer.pal(n = 3, name = "Set2")
# palete <-brewer.pal(n = 3, name = "Reds")
palete <- brewer.pal(n=6, "Purples")[2:6]

shapes <- c(21, 24)

# TOTI vs PLURI
png("Documents/zscan4/plots/toti_vs_pluri_diversity.png", width = 800, height = 700, res =150)  # El tamaño está en pulgadas



# PLURI VS TOI
palete1 <- brewer.pal(n=6, "Purples")[2:6]
p1<-ggplot(gene_expr_data_zscan, aes(x = pluripotent, y = totipotent, color = diversity, fill = diversity, shape = tsne_cluster)) +
  geom_point(size = 4, aes(shape = tsne_cluster)) +
  scale_color_gradientn(colours = palete1, name = "Diversity") +
  scale_fill_gradientn(colours = palete1, name = "Diversity") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes")+
  theme_classic()+
  #ggtitle("Zscan4+ mESC")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 19, color = "black", face = "bold"), # Título más oscuro y grande
    axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
    axis.text = element_text(size = 19, color = "black"),
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto en la leyenda
    legend.title = element_text(size = 16)  # Aumenta el tamaño del título de la leyenda
  ) +
  scale_shape_manual(values = shapes)

p1


# COLORED BY ZSCAN4 EXPRESSION
palete1 <- brewer.pal(n=6, "PuRd")[2:6]
p11<-ggplot(gene_expr_data_zscan, aes(x = pluripotent, y = totipotent, color = z4c, fill = z4c, shape = tsne_cluster)) +
  geom_point(size = 4, aes(shape = tsne_cluster)) +
  scale_color_gradientn(colours = palete1, name = "zscan4") +
  scale_fill_gradientn(colours = palete1, name = "zscan4") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes")+
  theme_classic()+
  ggtitle("Zscan4+ mESC")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 19, color = "black", face = "bold"), # Título más oscuro y grande
    axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
    axis.text = element_text(size = 19, color = "black"),
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto en la leyenda
    legend.title = element_text(size = 16)  # Aumenta el tamaño del título de la leyenda
  ) +
  scale_shape_manual(values = shapes)

p11
p1+p11
# 
gene_z4c <- subset(gene_expr_data, !z4c == 0)

ggplot(gene_z4c, aes(x = diversity, y = z4c, color = uncuts, fill = uncuts)) +
  geom_point(size = 4, aes(shape = tsne_cluster)) +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete, name = "% Uncuts") +
  labs(x = "Diversity", y = "Log Zscan4c expression") +
  theme_classic()
ggtitle("Zscan4c+ cells")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 14, color = "black"), # Título más oscuro y grande
    axis.title = element_text(size = 12, color = "black", face = "bold"), # Títulos de ejes más oscuros y grandes
    axis.text = element_text(size = 12, color = "black")+
      scale_shape_manual(values = shapes)  
  )



palete <- brewer.pal(n=9, "RdPu")[3:9]
p2<-ggplot(gene_expr_data, aes(x = diversity, y = z4c, color = uncuts, fill = uncuts)) +
  geom_point(size = 4, aes(shape = tsne_cluster)) +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete, name = "% Uncuts") +
  labs(x = "Diversity", y = "Log Zscan4c expression") +
  theme_classic()+
  ggtitle("Zscan4+ mESC")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 19, color = "black", face = "bold"), # Título más oscuro y grande
    axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
    axis.text = element_text(size = 19, color = "black")+
      scale_shape_manual(values = shapes)  
  )
ggsave("Documents/zscan4/sup_data/zscan_vs_diversity_plot.pdf", plot = p2, device = "pdf", width = 8, height = 7)


# dev.off()
p2
p2/p1

p2 <- ggplot(gene_expr_data, aes(x = diversity, y = z4c, color = uncuts, fill = uncuts)) +
  geom_point(size = 4, aes(shape = tsne_cluster)) +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete, name = "% Uncuts") +
  labs(x = "Diversity", y = "Log Zscan4c expression") +
  theme_classic()+
  #ggtitle("Zscan4+ mESC")+
  theme(
    plot.title = element_text(hjust = 0.6, size = 19, color = "black", face = "bold"), # Título más oscuro y grande
    axis.title = element_text(size = 19, color = "black"), # Títulos de ejes más oscuros y grandes
    axis.text = element_text(size = 19, color = "black"),
    legend.text = element_text(size = 14),  # Aumenta el tamaño del texto en la leyenda
    legend.title = element_text(size = 16)  # Aumenta el tamaño del título de la leyenda
  ) +
  scale_shape_manual(values = shapes) 

p2
#CORRELATION:
correlation_r <- cor.test(gene_expr_data$totipotent, gene_expr_data$pluripotent, method = "spearman")
print(paste("Spearman's rho:", correlation_r$estimate)) # -0.88
print(paste("P-value:", correlation_r$p.value)) # 1.46818579502331e-07

# ======================================
#            BOXPLOTS
# ======================================


my_colors <- brewer.pal(12, "Set3")[c(8, 3, 1,6)]  # Por ejemplo, colores 1, 4 y 7
# my_colors <- c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3")

DimPlot(sce_exp1, reduction = "tsne", group.by = "grupo", pt.size = 3, cols = my_colors)+
  ggtitle("Groups of cells")


boxplot_metadata_feature_grupo <- function(seurat_obj, feature_name) {
  metadata <- seurat_obj@meta.data
  
  plot_data <- data.frame(
    Feature = metadata[[feature_name]],
    grupo = factor(metadata$grupo)
  )
  levels(plot_data$grupo) <- c("Pluripotents", "Totipotents", "Transients")
  
  gg <- ggplot(plot_data, aes(x = grupo, y = Feature, fill = grupo)) +
    geom_boxplot(alpha = 0.5, width = 0.5, outliers = FALSE, aes(color = grupo)) +
    xlab("") +
    ylab(feature_name) +
    labs(title = paste0(feature_name, "distribution in each cluster")) +
    theme_minimal()
  
  # Definir los colores específicos del paquete Set3
  my_colors <- c("#FCCDE5","#BEBADA", "#8DD3C7", "#FDB462")
  
  gg + geom_jitter(data = plot_data, aes(color = grupo), width = 0.2, alpha = 1) +
    scale_fill_manual(values = my_colors) +  # Establecer los colores para los boxplots
    scale_color_manual(values = my_colors) +  # Establecer los colores para los puntos de datos
    ggtitle(feature_name) +
    theme(plot.title = element_text(hjust = 0.6))
}



boxplot_metadata_feature_grupo(sce_exp1, "Uncuts")
boxplot_metadata_feature_grupo(sce_exp1, "Diversity")
boxplot_metadata_feature_grupo(sce_exp1, "Mean_length")
boxplot_metadata_feature_grupo(sce_exp1, "gfp")



# ==============================================================
#       CORRELATION GENES EXPRESSION WITH DIVERSITY 
# ==============================================================

sce_cor <- sce_exp1[toti_genes,]
sce_cor <- sce_exp1

# 1) Extract vector of expression of each gene (37 cells)
# ----------------------------------------------
matrix_expresion <- as.matrix(sce_cor@assays$RNA$data)
# matrix_expresion <- as.matrix(log(sce_cor@assays$RNA$data+1))

order(colnames(matrix_expresion)) == order(colnames(sce_cor)) # verify the order is the same 
order(colnames(matrix_expresion)) == order(rownames(sce_cor@meta.data)) # verify the order is the same 

genes <- rownames(matrix_expresion)

vectores_expresion <- list() # create a list
for (gen in genes) {
  vector_gen <- matrix_expresion[gen, , drop = TRUE]  # drop = TRUE to mantain the vector structure
  vectores_expresion[[gen]] <- vector_gen
}

#example:
vectores_expresion$Sox21

# 2) vector of diversity of the 37 cells:
# -------------------------------------------------
diversidad <- sce_cor@meta.data$Diversity
dim(sce_cor)
length(diversidad)

# 3) Perform the correlation and save each result in a data frame
# -----------------------------------------------------------------
correlations <- sapply(vectores_expresion, function(gen_expresion) {
  cor.test(gen_expresion, diversidad, method = "spearman")
})


correlation_result <- cor.test(sce_z4c@meta.data$zscan4c, sce_z4c@meta.data$Diversity, method = "spearman")
print(paste("Spearman's rho:", correlation_result$estimate)) # -0.70
print(paste("P-value:", correlation_result$p.value)) # 0.00038


View(t(correlations))
# Warning messages:
#   1: In cor(gen_expresion, diversidad, method = "spearman") :
#   the standard deviation is zero
# 2: In cor(gen_expresion, diversidad, method = "spearman") :
#   the standard deviation is zero

df_correlations <- data.frame(
  Gen = names(correlations),
  Spearman_Correlation = correlations,
  stringsAsFactors = FALSE
)


genes <- colnames(matrix_expresion)  
diversidad <- sce_cor@meta.data$Diversity
vectores_expresion <- lapply(genes, function(gen) matrix_expresion[gen, , drop = TRUE])

resultados_correlacion <- lapply(names(vectores_expresion), function(gen) {
  gen_expresion <- vectores_expresion[[gen]]
  test_correlacion <- cor.test(gen_expresion, diversidad, method = "spearman")
  
  list(
    gen = gen,
    coeficiente_correlacion = test_correlacion$estimate,
    p_value = test_correlacion$p.value
  )
})

df_resultados <- do.call(rbind, lapply(resultados_correlacion, function(x) data.frame(x, row.names = NULL)))

View(df_resultados)



# >>>>>>>>>>>>>>>>>>>>>>>>
#         HEATMAP
# >>>>>>>>>>>>>>>>>>>>>>>>

# we come back to sce_exp1 where we have all the groups of cells OR, if we dont have that sce_exp1, we can do it from sce_proces as following:


# 1) Add to the metdata the number of group to compare between them 
# ---------------------------------------------------
sce_heatmap <- subset(sce_exp1, grupo == "transient" | grupo == "totipotents")


# 2) Heatmap
# ---------------------------------------------------
# ALL
sce_zscan4_heatmap<- sce_heatmap[all,]
expression_matrix<- GetAssayData(sce_zscan4_heatmap, layer = "data")
expression_matrix_ordered <- expression_matrix[match(all, rownames(expression_matrix)), ]

# _____________Diversity ________________________________________________________________
sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "Diversity")
barcode.diversity <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","Diversity")])
diversity <- barcode.diversity$Diversity

barcode.grupo <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","TSNE_cluster")])
group <- barcode.grupo$TSNE_cluster

color.df <- data.frame(
  COLOR_VALUE = diversity,
  barcode = barcode.diversity$Barcode
)
# color.df$Rangos <- cut(color.df$COLOR_VALUE, breaks = c(rangos, Inf), labels = FALSE, include.lowest = TRUE)
my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))
# color.df<-data.frame(COLOR_VALUE=diversity, color.name=my.colors(length(diversity)))or.df$Rangos]
color.df$color.name <- my.colors(length(diversity))

expression <- expression_matrix_ordered[ ,color.df$barcode]
expression <- expression_matrix_ordered[, order(color.df$COLOR_VALUE)]


order(colnames(expression)) == order(color.df$barcode)
expression_binary <- expression

expression_binary@x <- ifelse(expression_binary@x != 0, 1, 0)

pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df$color.name, col=pal, 
        main = "Pluripotent genes - z(+) cells ordered by Diversity", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")

# ___________________________________________________________________________________________


sce_zscan4_heatmap<- sce_heatmap[pluri_genes,]
dim(sce_zscan4_heatmap) #  1305   37
expression_matrix<- GetAssayData(sce_zscan4_heatmap, layer = "data")
expression_matrix_ordered <- expression_matrix[match(pluri_genes, rownames(expression_matrix)), ]

# _____________Diversity ________________________________________________________________
sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "Diversity")
barcode.diversity <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","Diversity")])
diversity <- barcode.diversity$Diversity

barcode.grupo <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","TSNE_cluster")])
group <- barcode.grupo$TSNE_cluster

color.df <- data.frame(
  COLOR_VALUE = diversity,
  barcode = barcode.diversity$Barcode
)
# color.df$Rangos <- cut(color.df$COLOR_VALUE, breaks = c(rangos, Inf), labels = FALSE, include.lowest = TRUE)
my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))
# color.df<-data.frame(COLOR_VALUE=diversity, color.name=my.colors(length(diversity)))or.df$Rangos]
color.df$color.name <- my.colors(length(diversity))

expression <- expression_matrix_ordered[ ,color.df$barcode]
# expression <- expression_matrix_ordered[, order(color.df$COLOR_VALUE)]

order(colnames(expression)) == order(color.df$barcode)
expression_binary <- expression

# expression_binary@x <- ifelse(expression_binary@x != 0, 1, 0)

pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
pal <- brewer.pal(n=3, "RdYlBu")
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df$color.name, col=pal, 
        main = "Pluripotent genes - z(+) cells ordered by Diversity", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")





####################################################
library(ComplexHeatmap)


# Crear datos de ejemplo: matriz de expresión y nombres de genes
set.seed(123)
expression_matrix <- matrix(rnorm(100 * 10), nrow = 100)
dim(expression_matrix)
rownames(expression_matrix) <- paste0("Gene", 1:100)
colnames(expression_matrix) <- paste0("Sample", 1:10)

pluri <- paste0("Gene", sample(1:50, 10, replace = FALSE))
toti <- paste0("Gene", sample(51:100, 10, replace = FALSE))

# Crear un vector lógico para los genes pluripotentes y totipotentes
genes_seleccionados <- rep(FALSE, nrow(expression_matrix))
genes_seleccionados[rownames(expression_matrix) %in% c(pluri, toti)] <- TRUE

# Seleccionar las filas de expression_matrix usando el vector lógico
exp <- expression_matrix[genes_seleccionados, ]

dim(exp)

# # Crear texto para anotaciones
# pluri_text <- lapply(unique(row_split), function(x) random_text(10))
# toti_text <- lapply(unique(row_split), function(x) random_text(10))

pluri_text <- pluri
toti_text <- toti


row_split <- c(rep("Pluripotent", length(pluri)), rep("Totipotent", length(toti)))
# row_split <- grupos


# names(pluri_text) <- rep("Pluripotent", length(pluri_text))
# names(toti_text) <- rep("Totipotent", length(toti_text))

grupos <- list(Pluripotent = pluri_text, Totipotent = toti_text)

text <- grupos

Heatmap(exp, name = "Expression", cluster_rows = FALSE,
        row_split = row_split,
        show_row_names = FALSE,
        right_annotation = rowAnnotation(textbox = anno_textbox(row_split, text)))

pdf("Documents/heatmap_PRUEBA_2.pdf")

# Crear el heatmap y guardarlo en el PDF
Heatmap(exp, name = "Expression", cluster_rows = FALSE,
        row_split = row_split,
        show_row_names = FALSE,
        right_annotation = rowAnnotation(textbox = anno_textbox(row_split, text)))

# Cerrar el dispositivo gráfico PDF
dev.off()

####################################################
# ___________ GROUP _________________________________________________________________________
sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "grupo")
barcode.grupo <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","grupo")])
group <- barcode.grupo$grupo

# NO SE SI ESTAS 3 LINEAS DE ABAJO CONTARIAN
color.df<-data.frame(COLOR_VALUE=group,  barcode = barcode.grupo$Barcode)

my.colors <- colorRampPalette(c("#FA8072", "#6C7B8B"))

unique_colors <- my.colors(length(unique(group)))
color.df$color.name <- unique_colors[match(color.df$COLOR_VALUE, unique(group))]

color.df_filtered <- color.df[color.df$COLOR_VALUE %in% c(1, 2), ]

expression <- expression_matrix[ ,color.df_filtered$barcode]
order(colnames(expression)) == order(color.df_filtered$barcode)
expression_binary <- expression

expression_binary@x <- ifelse(expression_binary@x != 0, 1, 0)

pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df_filtered$color.name, col=pal, 
        main = "Zscan (+) D=0 (13) vs D>0 (24)(ordered by group)", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")

# ___________________________________________________________________________________________

# ___________ CLUSTER _______________________________________________________________________
sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "TSNE_cluster")
barcode.grupo <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","TSNE_cluster")])
group <- barcode.grupo$TSNE_cluster

color.df<-data.frame(COLOR_VALUE=group,  barcode = barcode.grupo$Barcode)

my.colors <- colorRampPalette(c("#EEA2AD", "darkseagreen2"))

unique_colors <- my.colors(length(unique(group)))
color.df$color.name <- unique_colors[match(color.df$COLOR_VALUE, unique(group))]

color.df_filtered <- color.df[color.df$COLOR_VALUE %in% c(0, 3), ]

expression <- expression_matrix[ ,color.df_filtered$barcode]
order(colnames(expression)) == order(color.df_filtered$barcode)
expression_binary <- expression

pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df_filtered$color.name, col=pal, 
        main = "Zscan (+) D=0 vs D>0", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")

heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df_filtered$color.name, col=pal, 
        main = "cluster 0 (16 cells) and 3 (21 cells)", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4, margins = c(5, 5)) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")

# ___________________________________________________________________________________________


#                      DIVERSITY
# ________________________________________________

sce_zscan4_heatmap<- sce_heatmap[pluri_genes,]
dim(sce_zscan4_heatmap) #  1305   37
expression_matrix<- GetAssayData(sce_zscan4_heatmap, layer = "data")
expression_matrix_ordered <- expression_matrix[match(pluri_genes, rownames(expression_matrix)), ]

# _____________Diversity ________________________________________________________________
sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "Diversity")
barcode.diversity <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","Diversity")])
diversity <- barcode.diversity$Diversity

color.df <- data.frame(
  COLOR_VALUE = diversity,
  barcode = barcode.diversity$Barcode
)

my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))
color.df$color.name <- my.colors(length(diversity))

expression <- expression_matrix_ordered[ ,color.df$barcode]
order(colnames(expression)) == order(color.df$barcode)
expression_binary <- expression

# pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
pal <- brewer.pal(n=11, "RdYlBu")[8:11]
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df$color.name, col=pal, 
        main = "Pluripotent genes - z(+) cells ordered by Diversity", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")



# ________________________________________________

library(Seurat)
library(data.table)
library(RColorBrewer)
library(gplots)

pluri_genes
genestoti
top_200_genes
top_200_genes_plu

# Cargamos y preparamos los datos
sce_zscan4_heatmap <- sce_heatmap[top_200_genes_plu,]
expression_matrix <- GetAssayData(sce_zscan4_heatmap, layer = "data")
expression_matrix_ordered <- expression_matrix[match(top_200_genes_plu, rownames(expression_matrix)), ]

sce_zscan4_heatmap@meta.data <- sce_zscan4_heatmap@meta.data[order(sce_zscan4_heatmap@meta.data$Diversity), ]
barcode.diversity <- as.data.table(sce_zscan4_heatmap@meta.data[,c("Barcode","Diversity")])
color.df <- data.frame(
  COLOR_VALUE = barcode.diversity$Diversity,
  barcode = barcode.diversity$Barcode
)
my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))(length(diversity))
color.df$color.name <- my.colors#(length(diversity))

expression <- expression_matrix_ordered[,match(color.df$barcode, colnames(expression_matrix_ordered))]
pal <- brewer.pal(n=9, "RdYlBu")

# Creación del heatmap
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row",
        ColSideColors=color.df$color.name, col=pal,
        main = "Pluripotent genes - z(+) cells ordered by Diversity", 
        cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4)

# Añadir leyenda para colores de expresión
legend("topright", legend=c("Upregulation", "Downregulation"), 
       fill=c(pal[6], pal[3]), title="Expression", cex=0.8)  # Usando colores específicos de la paleta

col_fun <- colorRamp2(c(0, 0.5555), c("blue","white", "red")) 

pal <- rev(brewer.pal(n=11, "RdYlBu")[2:11])

Diversity_values <- sce_zscan4_heatmap@meta.data$Diversity


# Crear una función de colores que mapee los valores numéricos a un gradiente de colores
col_fun <- colorRamp2(c(min(Diversity_values), max(Diversity_values)),  c("#79CDCD", "#2F4F4F"))
full_pal <- brewer.pal(n=11, "RdYlBu")

# Selecciona un subconjunto específico asegurando que el rojo oscuro y el naranja están incluidos en los índices superiores
selected_pal <- full_pal[c(6, 1)]

# Ahora revierte la selección para que el naranja oscuro esté en los valores altos y el azul claro en los bajos
pal <- selected_pal
ht_toti <- Heatmap(as.matrix(expression),
                   name = "Relative expression",
                   col = pal,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   top_annotation = HeatmapAnnotation(df = data.frame(
                     Diversity = Diversity_values),
                     col = list(Diversity = col_fun),
                     show_legend = TRUE),
                   #col = list(Group = c("1" = "blue", "2" = "red"))),#HeatmapAnnotation(barcode = anno_block(gp = gpar(fill = color.df$color.name))),
                   
                   heatmap_legend_param = list(title = "Expression", at = c(-2, 0, 2), labels = c("Low", "Medium", "High")),
                   column_title = "Zscan4+ mESC ordered by Diversity",
                   row_title = "Top 200 Totipotent Genes",
                   column_title_gp = gpar(fontsize = 19, fontface = "bold"),  # Aumenta el tamaño de fuente y pone en negrita
                   row_title_gp = gpar(fontsize = 19, fontface = "bold")  # Aumenta el tamaño de fuente y pone en negrita
)
ht_toti

ht_pluri <- Heatmap(as.matrix(expression),
                    name = "Relative expression",
                    col = pal,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    top_annotation = HeatmapAnnotation(df = data.frame(
                      Diversity = Diversity_values),
                      col = list(Diversity = col_fun),
                      show_legend = TRUE),
                    #col = list(Group = c("1" = "blue", "2" = "red"))),#HeatmapAnnotation(barcode = anno_block(gp = gpar(fill = color.df$color.name))),
                    
                    heatmap_legend_param = list(title = "Expression", at = c(-2, 0, 2), labels = c("Low", "Medium", "High")),
                    column_title = "Zscan4+ mESC ordered by Diversity",
                    row_title = "Top 200 Pluripotent Genes",
                    column_title_gp = gpar(fontsize = 19, fontface = "bold"),  # Aumenta el tamaño de fuente y pone en negrita
                    row_title_gp = gpar(fontsize = 19, fontface = "bold")  # Aumenta el tamaño de fuente y pone en negrita
)

ht_pluri


ht_combined <- ht_toti + ht_pluri

# Ahora, para renderizar correctamente los dos heatmaps combinados
draw(ht_combined, heatmap_legend_side = "bot", annotation_legend_side = "bot")
ht_toti + ht_pluri

pu_rd_colors
Heatmap(as.matrix(expression), 
        name = "Relative expression", 
        col = color_palette, 
        show_row_names = FALSE,
        show_column_names = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        # clustering_method = "complete",
        top_annotation = HeatmapAnnotation(df = data.frame(
          Diversity = factor(sce_zscan4_heatmap@meta.data$Diversity, levels = c("1", "2"))
        ), col = list(Group = c("1" = "blue", "2" = "red"))),
        # bottom_annotation = rowAnnotation(link = anno_link(1:ncol(expression))),
        column_title = "Pluripotent genes - z(+) cells ordered by Diversity",
        row_title = "Genes",
        heatmap_legend_param = list(title = "Expression", at = c(-2, 0, 2), labels = c("Low", "Medium", "High"))
)

library(Seurat)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# Asumimos que ya tienes el objeto Seurat cargado y preprocesado
expression_matrix <- GetAssayData(sce_zscan4_heatmap, layer = "data")
expression_matrix_ordered <- expression_matrix[match(pluri_genes, rownames(expression_matrix)), ]

# Ordenando los datos basados en la metainformación
sce_zscan4_heatmap@meta.data <- sce_zscan4_heatmap@meta.data[order(sce_zscan4_heatmap@meta.data$Diversity), ]

# Ajuste de la matriz de expresión según el orden
expression <- expression_matrix_ordered[,match(sce_zscan4_heatmap@meta.data$Barcode, colnames(expression_matrix_ordered))]

# Definición de la paleta de colores para el gradiente
color_palette <- colorRampPalette(rev(brewer.pal(n = 11, "RdYlBu")))(100)

# Heatmap usando ComplexHeatmap
Heatmap(as.matrix(expression), 
        name = "Relative expression", 
        col = color_palette, 
        show_row_names = FALSE,
        show_column_names = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        # clustering_method = "complete",
        top_annotation = HeatmapAnnotation(df = data.frame(
          Diversity = factor(sce_zscan4_heatmap@meta.data$Diversity, levels = c("1", "2"))
        ), col = list(Group = c("1" = "blue", "2" = "red"))),
        # bottom_annotation = rowAnnotation(link = anno_link(1:ncol(expression))),
        column_title = "Pluripotent genes - z(+) cells ordered by Diversity",
        row_title = "Genes",
        heatmap_legend_param = list(title = "Expression", at = c(-2, 0, 2), labels = c("Low", "Medium", "High"))
)

#  ___________________________________________________________


