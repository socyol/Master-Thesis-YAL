##
## -------------------------------------------------------
#     
#             ZSCAN 4 - ANALYSIS  (II)
# -------------------------------------------------------
## ------------------------------------------------------

#                   Second script of the Analysis (TFM) 08.06.24

#               #   Script to analyze the cells and reads form the analysis of UMI (umi file extracted from cluster and compared to the clles that suprass QC)

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
library(ComplexHeatmap)

# load Seurat object
sce <- readRDS("YourPath")

# load reads table
umi_reads <- read.csv("YourPath")

# put barcode column in sce@metadata
sce@meta.data$Barcode <- rownames(sce@meta.data)

# intersect barcodes metadata and reads table:
colnames_to_keep <- intersect(colnames(sce), umi_reads$Barcode)

umi_reads_2 <- subset(umi_reads, Barcode %in% colnames_to_keep)
foo <- umi_reads_2[1:100,]
subset_foo <- subset(foo, select = c("Barcode", "UMI", "Match", "repetition_count", "total_count", "percentage"))

subset_umi <- subset(umi_reads_2, select = c("Barcode", "UMI", "Match", "repetition_count", "total_count", "percentage"))

# function reads table with the diversity and alignment:
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
read_table_function_modified <- function(row) {
  
  # Extract the BARCODES CASETE IN OTHER COLUMN TO PERFORM THE ALIGNMENT
  match_column <- row[["Match"]]
  row[["Length"]] <- nchar(match_column)
  
  # GLOBAL ALIGNMENT
  spacer <- DNAString("GGTATGCGGATGCAATCTCCGGGG")  
  seq_to_align <- DNAStringSet(as.character(match_column))  
  
  if (!exists("seq_to_align", envir = environment())) {
    return(NULL)  # Retorna NULL para eliminar la row
  } else {
    
    # sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
    alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global-local", substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
    
    if (!exists("alignment", envir = environment())) {
      return(NULL)
    } else {
      row$Alineamiento_Score <- alignment@score  
      row$Original <- ifelse(as.character(alignment@subject) == as.character(spacer), "YES", "NO")   
    }
  }
  
  return(row)
}

result_list <- lapply(1:nrow(subset_umi), function(i) read_table_function_modified(subset_umi[i, ]))
read_table <- do.call(rbind, result_list)
View(read_table)

# SAVE:
write.csv(read_table,"YourPath")
read_table <- read.csv("YourPath")

# Filters:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     1 - PROCESSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) delete those with length < 4 and >35
data <- subset(read_table, Length > 4 & Length < 35)
data <- na.omit(data)

# Function for indels:
spacer <- DNAString("GGTATGCGGATGCAATCTCCGGGG")  
indels_function <- function(row) {
  
  # Extract the BARCODES CASETE IN OTHER COLUMN TO PERFORM THE ALIGNMENT
  match_column <- row[["Match"]]
  
  # GLOBAL ALIGNMENT
  seq_to_align <- DNAStringSet(as.character(match_column))  
  
  if (!exists("seq_to_align", envir = environment())) {
    return(NULL)  # Return NULL to remove the row
  } else {
    if (row[["Length"]] == 24) {
      row$insertions <- 0
      row$deletions <- 0
      return(row)
    }
    
    # sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
    alignment <- pairwiseAlignment(spacer, seq_to_align, type = "global", substitutionMatrix = sigma, gapOpening = 5, gapExtension = 5)
    
    if (!exists("alignment", envir = environment())) {
      return(NULL)
    } else {
      subject <- toString(alignedSubject(alignment))
      pattern <- toString(alignedPattern(alignment))
      
      row$insertions <- sum(strsplit(pattern, "")[[1]] == "-")
      row$deletions <- sum(strsplit(subject, "")[[1]] == "-")
    }
  }
  return(row)
}

indel_result <- lapply(1:nrow(data), function(i) indels_function(data[i, ]))
indel_table <- do.call(rbind, indel_result)


write.csv(indel_table, "YourPath")

data <- read.csv("YourPath") # same as the one you saved before
#- ------------------------------------------------------------
colnames_to_keep <- intersect(colnames(sce), data$Barcode)
data <- subset(data, Barcode %in% colnames_to_keep)


# 2) Add a column "Diversity" with different scale of "Alignment_score" where the higher number, higher diversity
data$Diversity <- -data$Alineamiento_Score
data$Diversity <- (data$Diversity - min(data$Diversity)) / 
  (max(data$Diversity) - min(data$Diversity))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     2 - TABLE by BARCODE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subset_resultado <- data
subset_resultado %>% distinct(Barcode) %>% nrow() # 10976

cells_output <- subset_resultado %>%
  group_by(Barcode) %>%
  summarise(
    Coverage = n(),
    Diversity = mean(Diversity),
    Mean_length = mean(Length),
    Uncuts = sum(Original == "YES")/ n() *100,
    Insertions = sum(insertions),
    Deletions = sum(deletions))

cells_output <- cells_output[cells_output$Coverage > 8,]
write.csv(cells_output,"YourPath")
cells_output <- read.csv("YourPath")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     3 - MERGE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colnames_to_keep <- intersect(colnames(sce), cells_output$Barcode)
sce_zscan4_reads <- sce[, colnames_to_keep]

# sce_zscan4_reads@meta.data$RNA_snn_res.0.02 <- NULL
sce_zscan4_reads@meta.data$Barcode <- rownames(sce_zscan4_reads@meta.data)
sce_zscan4_reads@meta.data <- merge(sce_zscan4_reads@meta.data, cells_output, by = "Barcode")


saveRDS(sce_zscan4_reads, file ="YourPath")



###############################################################################################################################
#     FROM here is all the analysis of the sce_zscan4_UMIS-wreads.UNCUTS_abril
###############################################################################################################################

# ==============================
#             PLOTS
# ==============================

library(RColorBrewer)
palete <-brewer.pal(n = 3, name = "Set2")

FeaturePlot(sce_zscan4_reads, features = "Diversity" , reduction = "tsne", cols=  c("lightblue", "red"), pt.size = 3)

FeaturePlot(sce_zscan4_reads, features = "Coverage" , reduction = "tsne", cols=  c("lightblue", "red"), pt.size = 3)

DimPlot(sce_zscan4_reads, reduction = "tsne", pt.size = 3, cols = "Accent")
DimPlot(sce_zscan4_reads, reduction = "tsne", pt.size = 2, cols = "Accent")



library(Seurat)
library(ggplot2)

# Asumiendo que sce_zscan4_reads es tu objeto Seurat
# Genera el DimPlot para ver la estructura de los datos
p <- DimPlot(sce_zscan4_reads, reduction = "tsne", pt.size = 2, cols = "Accent", group.by = "seurat_clusters")
tsne_data <- as.data.frame(Embeddings(sce_zscan4_reads, reduction = "tsne"))
tsne_data$seurat_clusters <- sce_zscan4_reads@meta.data$seurat_clusters
centroids <- aggregate(cbind(tSNE_1 = tSNE_1, tSNE_2 = tSNE_2) ~ seurat_clusters, data = tsne_data, FUN = mean)
centroids$seurat_clusters <- paste("Cluster", centroids$seurat_clusters)
p + geom_text(data = centroids, aes(x = tSNE_1, y = tSNE_2, label = seurat_clusters), vjust = -0.2, hjust = -0.1,fontface = "bold", color = "black")


# 2) TSNE FUNCTION (UMAP IN zscan4_plots_abri.R)
# .......................................................
plot_metadata_feature <- function(so, feature_name) {
  tsne_coords <- as.data.frame(Embeddings(object = so, reduction = "tsne"))
  feature_value <- as.data.frame(so@meta.data[[feature_name]])
  plot_data <- cbind(tsne_coords, feature_value)
  
  max_feature <- max(feature_value)
  
  ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = so@meta.data[[feature_name]])) +
    geom_point(size = 2) +
    scale_color_gradient(low = "gray", high = "red", limits = c(0, max_feature)) +  # Ajustar los límites y la paleta de colores según lo deseado
    labs(x = "t-SNE 1", y = "t-SNE 2", color = feature_name) +  # Etiquetas de los ejes y la leyenda
    theme_classic()
  
}


plot_metadata_feature(sce_zscan4_reads, "Diversity")
plot_metadata_feature(sce_zscan4_reads, "Coverage")
plot_metadata_feature(sce_zscan4_reads, "Uncuts")
plot_metadata_feature(sce_zscan4_reads, "gfp")
plot_metadata_feature(sce_zscan4_reads, "zscan")

DimPlot(sce_zscan4_reads, reduction = "umap", pt.size = 3, cols = "Accent")



cat(" Cluster 0:", sum(sce_zscan4_reads@meta.data$TSNE_cluster == 0), "\n",
    "Cluster 1:", sum(sce_zscan4_reads@meta.data$TSNE_cluster == 1), "\n",
    "Cluster 2:", sum(sce_zscan4_reads@meta.data$TSNE_cluster == 2), "\n")


# Cluster 0: 1560 
# Cluster 1: 1391 
# Cluster 2: 24 


# 
# boxplot_data <- data.frame(
#   Uncuts = sce_zscan4_reads@meta.data$Uncuts,
#   TSNE_cluster = factor(sce_zscan4_reads@meta.data$TSNE_cluster),
#   diversity
# )
# levels(boxplot_data$TSNE_cluster) <- c("Cluster 0", "Cluster 1", "Cluster 2")
# # Create the ggplot with boxplot
# gg <- ggplot(boxplot_data, aes(x = TSNE_cluster, y = Uncuts, fill = TSNE_cluster)) +
#   geom_boxplot(alpha = 0.5, width = 0.5, outliers = FALSE, aes(color = TSNE_cluster)) +
#   stat_boxplot(geom = "errorbar", width = 0.25)+
#   xlab("") +
#   ylab("Uncuts") +
#   labs(title = "% Uncuts distribution in each cluster") +
#   theme_minimal()
# 
# # Add data points for each level of TSNE_cluster
# gg + geom_jitter(data = boxplot_data, aes(color = TSNE_cluster), width = 0.2, alpha = 0.6)+
#   scale_fill_brewer(palette = "Accent") +  # Establecer la paleta de colores para los boxplots
#   scale_color_brewer(palette = "Accent")+
#   ggtitle("% Uncuts distribution in each cluster") + 
#   theme(plot.title = element_text(hjust = 0.6))  

# PLOTS ---> ZSCAN4_PLOTS_ABRIL.R


# for (cluster in unique(boxplot_data$TSNE_cluster)) {
#   subset_data <- boxplot_data[boxplot_data$TSNE_cluster == cluster, ]
#   Q1 <- quantile(subset_data$Uncuts, 0.25)
#   Q3 <- quantile(subset_data$Uncuts, 0.75)
#   IQR <- Q3 - Q1
#   
#   # Calcular los límites para identificar los valores atípicos
#   lower_bound <- Q1 - 1.5 * IQR
#   upper_bound <- Q3 + 1.5 * IQR
#   
#   # Identificar los valores atípicos para este cluster
#   outliers <- subset_data$Uncuts[subset_data$Uncuts < lower_bound | subset_data$Uncuts > upper_bound]
#   
#   # Imprimir los valores atípicos para este cluster
#   cat("Cluster:", cluster, "\n")
#   cat("Outliers:", outliers, "\n\n")
# }


# 2) BOXPLOT function WITHOUT outliers:
# ................................................
# FUNCTION BOXPLOT WITHOUT OUTLIERS
boxplot_metadata_feature <- function(seurat_obj, feature_name) {
  metadata <- seurat_obj@meta.data
  
  plot_data <- data.frame(
    Feature = metadata[[feature_name]],
    TSNE_cluster = factor(metadata$TSNE_cluster)
  )
  levels(plot_data$TSNE_cluster) <- c("Cluster 0", "Cluster 1", "Cluster 2")
  
  # Calcular los cuantiles y el rango intercuartílico para cada cluster
  for (cluster in unique(plot_data$TSNE_cluster)) {
    subset_data <- plot_data[plot_data$TSNE_cluster == cluster, ]
    Q1 <- quantile(subset_data$Feature, 0.25)
    Q3 <- quantile(subset_data$Feature, 0.75)
    IQR <- Q3 - Q1
    
    # Calcular los límites para identificar los valores atípicos
    lower_bound <- Q1 - 1.5 * IQR
    upper_bound <- Q3 + 1.5 * IQR
    
    # Filtrar los datos para excluir los valores atípicos
    plot_data <- plot_data[!(plot_data$TSNE_cluster == cluster & 
                               (plot_data$Feature < lower_bound | plot_data$Feature > upper_bound)), ]
  }
  
  gg <- ggplot(plot_data, aes(x = TSNE_cluster, y = Feature, fill = TSNE_cluster)) +
    geom_boxplot(alpha = 0.5, width = 0.5, aes(color = TSNE_cluster), outliers = FALSE) +
    stat_boxplot(geom = "errorbar", width = 0.25)+
    xlab("") +
    ylab(feature_name) +
    labs(title = paste0(feature_name, " distribution in each cluster")) +
    theme_minimal() +
    geom_jitter(aes(color = TSNE_cluster), width = 0.2, alpha = 1) +
    scale_fill_brewer(palette = "Accent") +
    scale_color_brewer(palette = "Accent") +
    ggtitle(feature_name) +
    theme(plot.title = element_text(hjust = 0.6))
  
  return(gg)
}


boxplot_metadata_feature(sce_zscan4_reads, "Uncuts")
boxplot_metadata_feature(sce_zscan4_reads, "Diversity")
boxplot_metadata_feature(sce_zscan4_reads, "Mean_length")







# ==============================
#        ZSCAN4 GENES
# ==============================
all.genes<-rownames(sce_zscan4_reads)
genes_df <- data.frame(Gene = all.genes)

# we only take the zscan4c(+) cells (we try later with the DEA genes
genes_zscan <- all.genes[grepl("^Zscan4", all.genes)]
FeaturePlot(sce_zscan4_reads, features = genes_zscan, reduction = "tsne")

zscan4_expression <- GetAssayData(object = sce_zscan4_reads, slot = "data")[genes_zscan, ]  # Reemplaza 'gene_id' por el ID de Zscan4f
valores_no_cero <- zscan4_expression[zscan4_expression != 0]

max(sce_zscan4_reads@assays$RNA$counts["gRNA",])

length(valores_no_cero) #113

# Create a dataframeto sum the values of each gene and then add it to the metadata
frame <- as.data.frame(t(zscan4_expression))
View(frame)
frame$suma <- rowSums(frame)
sum(frame$suma >0) # 44 cells with zscan4 expression
head(frame)

# add it to the metadata
sce_zscan4_reads@meta.data$zscan4_expr <- frame$suma # si es la suma de todos las las celulas con la expresion de cualquier zscan4 gene
zscan4_expression <- sce_zscan4_reads@meta.data$zscan4_expr
sce_zscan4_reads@meta.data$zscan4_expr <- zscan4_expression # si es para la expresion de un zscan4 gene concreto

# Logarithm ---------------------------------
sce_zscan4_reads@meta.data$zscan4_expr_log <- log2(sce_zscan4_reads@meta.data$zscan4_expr + 1)  # Se suma 1 para evitar log(0)
zscan4_expression_log <- sce_zscan4_reads@meta.data$zscan4_expr_log 
# -------------------------------------------


# foo is after deleting cluster 3 (totipotent cluster) to see in the feature plot the distribution of zscan4 expression better
foo <- subset(sce_zscan4_reads, TSNE_cluster != 2)

FeaturePlot(sce_zscan4_reads, features = "zscan4_expr_log" , reduction = "tsne", cols=  c("lightblue", "red"), pt.size = 3)

FeaturePlot(foo, features = "zscan4_expr_log" , reduction = "tsne", cols=  c("lightblue", "red"), pt.size = 3)
View(foo@meta.data) #
#       -------------> WE CAN SEE 5 CELLS FROM CLUSTER 1 WITH ZSCAN4 EXPRESSION
# 

#     I DECIDED TO DELETE THESE CELLS:

plot_metadata_feature(sce_zscan4_reads_filtered, "gfp")
plot_metadata_feature(sce_zscan4_reads_filtered, "Uncuts")
plot_metadata_feature(sce_zscan4_reads_filtered, "Diversity")

# ==========================================
#      DELETE CELLS CLUSTER 1 ZSCAN4(+)
# ==========================================
sce_zscan4_reads_filtered <- subset(sce_zscan4_reads, TSNE_cluster != 1 | zscan4_expr <=0)
View(sce_zscan4_reads_filtered@meta.data)

dim(sce_zscan4_reads_filtered@meta.data) # 2967
dim(sce_zscan4_reads@meta.data) # 2975

# 1) repeat the plots:
# ....................................
DimPlot(sce_zscan4_reads_filtered, reduction = "tsne", pt.size = 3, cols = "Accent")

cat(" Cluster 0:", sum(sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 0), "\n",
    "Cluster 1:", sum(sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 1), "\n",
    "Cluster 2:", sum(sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 2), "\n")

boxplot_metadata_feature(sce_zscan4_reads_filtered, "Uncuts")
boxplot_metadata_feature(sce_zscan4_reads_filtered, "Diversity")
boxplot_metadata_feature(sce_zscan4_reads_filtered, "Mean_length")
boxplot_metadata_feature(sce_zscan4_reads_filtered, "gfp")

# 2) CHANGE THE NAMES ¿?

# 2)  Remove cluster 1
sce_cl02 <- subset(sce_zscan4_reads_filtered, TSNE_cluster != 1)




# ALL THE ZSCAN4 SEPARADOS
zd <- all.genes[grepl("^Zscan4d", all.genes)]
zscan4d_expression <- GetAssayData(object = sce_cl02, slot = "data")[zd, ] 
sce_cl02@meta.data$zscan4d<- zscan4d_expression 

zc <-all.genes[grepl("^Zscan4c", all.genes)]
zscan4c_expression <- GetAssayData(object = sce_cl02, slot = "data")[zc, ] 
sce_cl02@meta.data$zscan4c<- zscan4c_expression 

zf <-all.genes[grepl("^Zscan4f", all.genes)]
zscan4f_expression <- GetAssayData(object = sce_cl02, slot = "data")[zf, ] 
sce_cl02@meta.data$zscan4f<- zscan4f_expression 




#----------------------------------------------------------------
# Plots to see diversity and coverage with ZSCAN4 EXPRESSION:
# ----------------------------------------------------------------
uncuts_percentage <- sce_cl02@meta.data$Uncuts
diversity <- sce_cl02@meta.data$Diversity
coverage <- sce_cl02@meta.data$Coverage
cluster <- sce_cl02@meta.data$TSNE_cluster
zscan4_expression_log <- sce_cl02@meta.data$zscan4_expr 
mit <- sce_cl02@meta.data$percent.mt
z4d <- sce_cl02@meta.data$zscan4d
z4c <- sce_cl02@meta.data$zscan4c
z4f <- sce_cl02@meta.data$zscan4f

sce_cl02@meta.data$B <- rownames(sce_cl02@meta.data)

# View(sce_zscan4_reads_filtered@meta.data)

# dispersion plot
plot(uncuts_percentage, zscan4_expression_log, 
     xlab = "% Uncuts", ylab = "Log Zscan4 expression")

# View(sce_zscan4_reads@meta.data)

data <- data.frame(uncuts_percentage, zscan4_expression_log, coverage, diversity, cluster, z4d, z4f, z4c)

# color_palette <- colorRampPalette(c("gray", "blue"))(100)
# color_palette <- colorRampPalette(c("#C1CDC1","#9F79EE", "red"))(100)

palete2 <- brewer.pal(n=4, "Purples")

ggplot(data, aes(x = diversity, y = zscan4_expression_log, color = uncuts_percentage, fill = uncuts_percentage)) +
  geom_point(size = 4, shape = 21, color = "gray") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()

ggplot(data, aes(x = coverage, y = zscan4_expression_log, color = uncuts_percentage, fill = uncuts_percentage)) +
  geom_point(size = 4, shape = 21, color = "gray") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()


# colored by cluster
ggplot(data, aes(x = diversity, y = zscan4_expression_log, color = cluster)) +
  geom_point(size = 5) +
  scale_color_manual(values = palete, name = "cluster") +
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()


# Zscan4D
ggplot(data, aes(x = diversity, y = z4d, color = uncuts_percentage, fill = uncuts_percentage)) +
  geom_point(size = 4, shape = 21, color = "gray") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "Zscan4d expression") +
  theme_minimal()
# Zscan4C
ggplot(data, aes(x = diversity, y = z4c, color = uncuts_percentage, fill = uncuts_percentage)) +
  geom_point(size = 4, shape = 21, color = "gray") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "Zscan4c expression") +
  theme_minimal()
# Zscan4F
ggplot(data, aes(x = diversity, y = z4f, color = uncuts_percentage, fill = uncuts_percentage)) +
  geom_point(size = 4, shape = 21, color = "gray") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "Zscan4f expression") +
  theme_minimal()


shapes <- c(21, 24)
# CLUSTER AND UNCUTS
ggplot(data, aes(x = diversity, y = z4c, color = uncuts_percentage, fill = uncuts_percentage, shape = TSNE-cluster)) +
  geom_point(size = 4, aes(shape = cluster), color ="gray59") +  # Contorno negro por defecto
  scale_color_gradientn(colours = palete2, name = "% Uncuts") +
  scale_fill_gradientn(colours = palete2, name = "% Uncuts") +
  labs(x = "Diversity", y = "Zscan4c expression", title = "Zscan4 c") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = shapes)  


# ------------------------------------------------------------------------------------------



# ..................................................................................................................................................
# .................................................................................................................................................
# ..................................................................................................................................................
# ..................................................................................................................................................



#             ============================================
#                        DIFFERENTIAL EXPRESSION
#            ============================================

# Differential expression analysis between the cluster 0 and cluster 3 (the whole cluster 0 and whole cluster 3)
# GOAL: Run FindMarkers() function between the whole cluster 0 and 3 (extracting the 12 zscan(+) cells)

# sce_zscan4_reads <- readRDS("Documents/sce_zscan_values_DEF")


# 1) Perform find markers between the cluster 0 and 2
# -------------------------------------------------------------------------------------------------
sce_proces_2 <- sce_cl02
Idents(sce_proces_2)
markers_02 <- FindMarkers(sce_proces_2, ident.1 = "0", ident.2 = "2")

markers_02$p_val_adj = p.adjust(markers_02$p_val, method='fdr')
markers_02$gene <- rownames(markers_02) 

markers02<- subset(markers_02, p_val_adj < 0.05)
markers02 # MARKERS BETWEEN CLUSTER 0 AND 2
dim(markers02) # 1796    6

View(markers02) 

# ordenar de mayor logFC a menor:
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
write.csv(markers02_ordered, file = "Documents/zscan4/markers02_abril.csv")

markers02_ordered <- read.csv("Documents/markers02.csv")

# 2) see what genes are in pluripotency
# ------------------------------------------> SUBSET DOWN OR UP (totigenes)
all.genes<-rownames(sce_proces_2)
genes_df <- data.frame(Gene = all.genes)
pluri_genes <- markers02_ordered$gene[markers02_ordered$avg_log2FC>1]
toti_genes <- markers02_ordered$gene[markers02_ordered$avg_log2FC < -1]
all <- markers02_ordered$gene[markers02_ordered$avg_log2FC < -1 | markers02_ordered$avg_log2FC>1]


cat(" Genes pluri:", length(pluri_genes), "\n",
    "Genes toti:", length(toti_genes))

# NOW:
# Genes pluri: 230 
# Genes toti: 1075

#       BEFORE:
# Genes pluri: 226
# Genes toti: 1065

# markers_02_pluri <- subset(markers02, rownames(markers02) %in% pluri)
# markers_02_toti <- subset(markers02, rownames(markers02) %in% toti)
# 
# markers_02_all <- subset(markers02, rownames(markers02) %in% markers02$gene)

# View(markers_03_pluri)
# View(markers_03_toti)







# >>>>>>>>>>>>>>>>>>>>>>>>
#     GROUP OF CELLS
# >>>>>>>>>>>>>>>>>>>>>>>>

library(readxl)
library(openxlsx)
library(data.table)
library(SeuratObject)
library(psychTools)
library(data.table)
library(RColorBrewer)



# Exp 27-02-24:

#       Heatmap con las celulas divididas en 4 grupos diferentes (luego se tienen que ordenar los barcodes como lo de diversity)

# --------------------------------
#  1. SELECT CELLS
# --------------------------------

# 1) SELECT CELLS WITH ZSCAN GENES -------> Group 1 and 2
# '''''''''''''''''''''''''''''''''''''''''''''''''''''
# sce_zscan4_reads <- readRDS("Documents/sce_zscan4_wreads-UNCUTS.rds")

all.genes<-rownames(sce_zscan4_reads_filtered)
genes_df <- data.frame(Gene = all.genes)

# save the sce with zscan values
saveRDS(sce_zscan4_reads, "Documents/zscan4/sce_zscan4_reads_filtered_DEFabril.rds") # -------> this object has both reads, umis and zscan4 expression levels

toti <- rownames(sce_zscan4_reads_filtered@meta.data[sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 2, ]) 
transient <- rownames(sce_zscan4_reads_filtered@meta.data[sce_zscan4_reads_filtered@meta.data$zscan4_expr > 0 & sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 0, ]) 
control<- rownames(sce_zscan4_reads_filtered@meta.data[sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 1, ]) 
pluri<- rownames(sce_zscan4_reads_filtered@meta.data[sce_zscan4_reads_filtered@meta.data$TSNE_cluster == 0 & sce_zscan4_reads_filtered@meta.data$zscan4_expr == 0, ]) 

cat(" number of TOTIPOTENT cells -->", length(toti), "\n",
    "number of TRANSIENT cells --> ", length(transient), "\n",
    "number of GFP(-) cells --> ", length(control), "\n",
    "number of PLURIPOTENT cells --> ", length(pluri), "\n")

# number of TOTIPOTENT cells --> 24 
# number of TRANSIENT cells -->  15 
# number of GFP(-) cells -->  1383 
# number of PLURIPOTENT cells -->  1545 

# --------------------------------
# 2. Create metadata only with those cells
# --------------------------------

cells_to_keep <- c(toti, pluri, control,transient)
length(cells_to_keep) # 2967

sce_exp1<- sce_zscan4_reads_filtered[, cells_to_keep]
dim(sce_exp1) # 23338  2967
sce_exp1@meta.data$grupo <- ifelse(colnames(sce_exp1) %in% toti, "totipotents", # number for each group
                                   ifelse(colnames(sce_exp1) %in% transient, "transient", 
                                          ifelse(colnames(sce_exp1) %in% pluri,"pluripotents", "GFP(-)")))

View(sce_exp1@meta.data)
cat(" Number of cells in each Group", "\n",
    "..... Totipotents =", sum(sce_exp1@meta.data$grupo == "totipotents"), "\n",
    "..... Transient =", sum(sce_exp1@meta.data$grupo == "transient"), "\n",
    "..... Pluripotents =", sum(sce_exp1@meta.data$grupo == "pluripotents"), "\n",
    "..... control GFP(-) =", sum(sce_exp1@meta.data$grupo == "GFP(-)"), "\n")

# Number of cells in each Group 
# ..... Totipotents = 24 
# ..... Transient = 15 
# ..... Pluripotents = 1545 
# ..... control GFP(-) = 1383 

# AFTER UMIS:
# Number of cells in each Group 
# ..... D00_zpos (1) = 13 
# ..... D_zpos (2) = 24 
# ..... D00_zneg (3) = 267 
# ..... control (4) = 127 


# --------------------------------
#         3. PLOTS
# --------------------------------
# tSNE of the sce colored by grupo
DimPlot(sce_exp1, group.by = "grupo", reduction = "tsne", pt.size = 3)+
  ggtitle("Groups of cells")

sum(sce_exp1@meta.data$grupo == 1 & sce_exp1@meta.data$TSNE_cluster == 2 )




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
  levels(plot_data$grupo) <- c("GFP(-)", "Pluripotents", "Totipotents", "Transients")
  
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





# >>>>>>>>>>>>>>>>>>>>>>>>
#         HEATMAP
# >>>>>>>>>>>>>>>>>>>>>>>>

# we come back to sce_exp1 where we have all the groups of cells OR, if we dont have that sce_exp1, we can do it from sce_proces as following:


# 1) Add to the metdata the number of group to compare between them 
# ---------------------------------------------------
View(sce_exp1@meta.data)
sce_heatmap <- subset(sce_exp1, grupo == "transient" | grupo == "totipotents")
dim(sce_heatmap) #  23338    39


# 2) Heatmap
# ---------------------------------------------------

# ALL
sce_zscan4_heatmap<- sce_heatmap[all,]
dim(sce_zscan4_heatmap) #  1305   37
expression_matrix<- GetAssayData(sce_zscan4_heatmap, slot = "data")
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

#               FOR A HEATMAP WITH GENE ANNOTATION:
#               ........................................................................
color_values <- colorRamp(my_colors)(100)
color_vector <- color_values[findInterval(color.df$COLOR_VALUE, seq(min(color.df$COLOR_VALUE), max(color.df$COLOR_VALUE), length.out = 100))]

ha_diversity <- HeatmapAnnotation(Diversity = diversity,
                                  Cluster = group)

expression <- expression[,order(color.df$COLOR_VALUE)]
colnames(expression)

ha_diversity <- HeatmapAnnotation(Diversity = color.df$COLOR_VALUE)
col_names <- color.df$color.name


names(col_names) <- color.df$COLOR_VALUE
draw(ha_diversity)
genes_totipotent <- intersect(rownames(expression), toti)
genes_pluripotent <- intersect(rownames(expression), pluri)

row_split <- c(rep("Pluripotent", length(genes_pluripotent)), rep("Totipotent", length(genes_totipotent)))

grupos <- list(Pluripotent = genes_pluripotent, Totipotent = genes_totipotent)

text <- grupos
class(expression)
my_colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))(100)


pdf("Documents/heatmap_170424.pdf")
Heatmap(as.matrix(expression), name = "Expression", cluster_rows = FALSE, cluster_columns = FALSE,
        row_split = row_split,
        show_row_names = FALSE,
        color_space = "LAB",
        column_names_gp = gpar(fontsize = 6),
        top_annotation = ha_diversity)

# right_annotation = rowAnnotation(textbox = anno_textbox(row_split, text)))

dev.off()

#             .............................................................................

order(colnames(expression)) == order(color.df$barcode)
expression_binary <- expression

expression_binary@x <- ifelse(expression_binary@x != 0, 1, 0)

pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale="row", ColSideColors=color.df$color.name, col=pal, 
        main = "Pluripotent genes - z(+) cells ordered by Diversity", cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4) #, main.title = "Diversity Heatmap", main.subtitle = "Subtitle Text")

# ___________________________________________________________________________________________





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







#
#             DIFFERENTIAL EXPRESSION BETWEEN ZSCAN4 CELLS
#           --------------------------------------------------

# Gonna do the DEA between the cells zscan4 with high diversity and low diversity, and then between those in cluster 0 and cluster 3:
View(sce_zscan4_groups@meta.data)

cat(" Cells with D == 0 -->", sum(sce_zscan4_groups@meta.data$grupo == 1), "\n",
    "Cells with D > 0 -->", sum(sce_zscan4_groups@meta.data$grupo == 2), "\n")
# Cells with D == 0 --> 13 
# Cells with D > 0 --> 24 

cat(" Cells in toti -->", sum(sce_zscan4_groups@meta.data$TSNE_cluster == 2), "\n",
    "Cells in pluri -->", sum(sce_zscan4_groups@meta.data$TSNE_cluster == 0), "\n")

# Cells in toti --> 21 
# Cells in pluri --> 16 

cat(" Cells in toti with D==0 -->", sum(sce_zscan4_groups@meta.data$TSNE_cluster == 2 & sce_zscan4_groups@meta.data$grupo == 1), "\n",
    "Cells in toti WITH D>0 -->", sum(sce_zscan4_groups@meta.data$TSNE_cluster == 2 & sce_zscan4_groups@meta.data$grupo == 2), "\n")

# Cells in toti with D==0 --> 11 
# Cells in toti WITH D>0 --> 10 

# -------> PERFORM THE DEA INSIDE CLUSTER 3:

sce_3 <- subset(sce_zscan4_groups, TSNE_cluster == 2)
boxplot(sce_3@meta.data$Diversity ~ sce_3@meta.data$TSNE_cluster, 
        xlab = "TSNE cluster", ylab = "Diversity")

diversity001<- rownames(sce_3@meta.data[sce_3@meta.data$Diversity > 0.01, ])
diversity0<- rownames(sce_3@meta.data[sce_3@meta.data$Diversity < 0.01, ])

cells <- c(diversity0, diversity001)
length(cells) # 21

sce_3_1<- sce_3[, cells]
dim(sce_3_1) #  23338    21
sce_3_1@meta.data$exp_div <- ifelse(colnames(sce_3_1) %in% diversity001, 1, 0) 

Idents(sce_3_1) <- "exp_div"
markers_0_1 <- FindMarkers(sce_3_1, ident.1 = "0", ident.2 = "1")
markers_0_1$p_val_adj = p.adjust(markers_0_1$p_val, method='fdr')
markers_03$gene <- rownames(markers_03) 

markers01<- subset(markers_0_1, p_val_adj < 0.05)
# 0
head(markers_0_1) # -----------------------> CONCLUSION: THERE IS NOT ANY DE GENE

# p_val avg_log2FC pct.1 pct.2 p_val_adj
# Mlana  0.0002358873  -3.534817 0.000  0.75         1
# Csrnp1 0.0007763394  -3.140727 0.176  1.00         1
# Scrn2  0.0011363448  -2.508342 0.176  1.00         1
# Polr2a 0.0013366750   1.426936 1.000  1.00         1
# Itm2c  0.0013731517  -3.455344 0.059  0.75         1
# Ggh    0.0013731517  -2.958925 0.059  0.75         1


# -------> PERFORM THE DEA with ZSCAN4 CELLS (group 1 and 2):
#                       we take the sce from the one done with the heatmap:



sce_12 <- sce_zscan4_groups
boxplot(sce_12@meta.data$Diversity ~ sce_12@meta.data$grupo, 
        xlab = "Group", ylab = "Diversity")


Idents(sce_3_1) <- "grupo"
markers_0_1 <- FindMarkers(sce_3_1, ident.1 = "1", ident.2 = "2")
markers_0_1$p_val_adj = p.adjust(markers_0_1$p_val, method='fdr')
markers_2<- subset(markers_0_1, p_val_adj < 0.05)
# -----------------------> CONCLUSION: THERE IS ONLY ONE GENE BUT WITH A LFC OF -1.8
#              p_val avg_log2FC pct.1 pct.2  p_val_adj
# Lsm12 3.507701e-06  -1.819663 0.538     1 0.04565624


### ---------------------------> NOW--> gonna do a plot to see what cells inside this category have more diversity


# ________________________________________________________________
# Plots to see diversity and coverage with ZSCAN4 EXPRESSION:
# ----------------------------------------------------------------
zscan4_expression_log <- sce_12@meta.data$zscan4_expr_log 
uncuts_percentage <- sce_12@meta.data$Uncuts
group <- sce_12@meta.data$grupo
cluster <- sce_12@meta.data$TSNE_cluster
diversity <- sce_12@meta.data$Diversity
coverage <- sce_12@meta.data$Coverage
head(diversity)

# sce_12@meta.data$B <- rownames(sce_12@meta.data)

View(sce_12@meta.data)

# dispersion plot
plot(uncuts_percentage, zscan4_expression_log, 
     xlab = "% Uncuts", ylab = "Zscan4f expression")

# View(sce_zscan4_reads@meta.data)

data <- data.frame(uncuts_percentage, zscan4_expression_log, coverage)

color_palette <- colorRampPalette(c("gray", "blue"))(100)

color_palette <- colorRampPalette(c("#C1CDC1","#9F79EE", "red"))(100)

# ----------->>>> DIVERSITY colored by UNCUTS
data <- subset(sce_12@meta.data[, c("Uncuts", "zscan4_expr_log", "Diversity")])

ggplot(data, aes(x = Diversity, y = zscan4_expr_log, color = Uncuts)) +
  geom_point(size = 3) +
  scale_color_gradientn(colours = color_palette, name = "% Uncuts") +
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Diversity")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

# ------------>>>> COVERAGE colored by UNCUTS
data <- subset(sce_12@meta.data[, c("Uncuts", "zscan4_expr_log", "Coverage")])

ggplot(data, aes(x = Coverage, y = zscan4_expr_log, color = Uncuts)) +
  geom_point(size = 3) +
  scale_color_gradientn(colours = color_palette, name = "% Uncuts") +
  labs(x = "Coverage", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Coverage")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )


# COLORED BY GROUP
# .....................
data <- subset(sce_12@meta.data[, c("grupo", "zscan4_expr_log", "Diversity", "Coverage")])
data$grupo <- as.factor(data$grupo)
color_palette <- c("#FA8072", "#6C7B8B")

ggplot(data, aes(x = Diversity, y = zscan4_expr_log, color = grupo)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, name = "Group") +  # Escala de color discreto
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Diversity")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(data, aes(x = Coverage, y = zscan4_expr_log, color = grupo)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, name = "Group") +  # Escala de color discreto
  labs(x = "Coverage", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Coverage")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )



# COLORED BY CLUSTER
# ........................
data <- subset(sce_12@meta.data[, c("TSNE_cluster", "zscan4_expr_log", "Diversity", "Coverage")])

data$TSNE_cluster <- as.factor(data$TSNE_cluster) # convert that column into factor

color_palette <-c("#EEA2AD", "darkseagreen3")

ggplot(data, aes(x = Diversity, y = zscan4_expr_log, color = TSNE_cluster)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, name = "cluster") +  # Escala de color discreto
  labs(x = "Diversity", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Diversity")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(data, aes(x = Coverage, y = zscan4_expr_log, color = TSNE_cluster)) +
  geom_point(size = 5) +
  scale_color_manual(values = color_palette, name = "cluster") +  # Escala de color discreto
  labs(x = "Coverage", y = "log Zscan4 expression") +
  theme_minimal()+
  ggtitle("Zscan4 expression vs Coverage")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

# ---------------> MAYBE WE CAN STUDY THE DIFFERENCES BETWEEN THIS CELLS THAT HAVE AROUND 50% OF UNCUTS
sce_3 <- sce_zscan4_groups
boxplot(sce_3@meta.data$Diversity ~ sce_3@meta.data$TSNE_cluster, 
        xlab = "TSNE cluster", ylab = "Diversity")


# 1) Select groups to separate the cells depending on the diversity they have 
# # .............................................................................
diversity002<- rownames(sce_3@meta.data[sce_3@meta.data$Diversity > 0.25, ])
diversity0<- rownames(sce_3@meta.data[sce_3@meta.data$Diversity < 0.25, ])

cells <- c(diversity0, diversity002)
length(cells) # 21

sce_3_1<- sce_3[, cells]
dim(sce_3_1) #  23338    21
sce_3_1@meta.data$exp_div <- ifelse(colnames(sce_3_1) %in% diversity002, 1, 0) 


# 2) perform the DEA
# ..........................
Idents(sce_3_1) <- "exp_div"
markers_0_1 <- FindMarkers(sce_3_1, ident.1 = "0", ident.2 = "1")

dim(markers_0_1) # 14154     5

markers_0_1$p_val_adj = p.adjust(markers_0_1$p_val, method='fdr')
markers_0_1$gene <- rownames(markers_0_1) 
markers_2<- subset(markers_0_1, p_val_adj < 0.05)
markers_2

#                      p_val avg_log2FC pct.1 pct.2   p_val_adj          gene
# Atxn7         5.571366e-07  -4.134822 0.029 1.000 0.001999908         Atxn7
# Hhat          1.978149e-06  -4.230408 0.000 0.667 0.001999908          Hhat
# Scn3a         1.978149e-06  -4.214468 0.000 0.667 0.001999908         Scn3a
# Pde11a        1.978149e-06  -4.230408 0.000 0.667 0.001999908        Pde11a
# Srrm4         1.978149e-06  -4.214468 0.000 0.667 0.001999908         Srrm4
# Gm12784       1.978149e-06  -4.834995 0.000 0.667 0.001999908       Gm12784
# Fes           1.978149e-06  -4.284498 0.000 0.667 0.001999908           Fes
# Lmo1          1.978149e-06  -4.284498 0.000 0.667 0.001999908          Lmo1


View(markers_2)

# these genes corresponds to PLURIPOTENT genes, since they are DOWNREGULATED to the group with low diversity


# ---------> DUDA: ALOMEJOR TENOG QUE NORMALIZAR EL COVERAGE????





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FUNCTION TO DO THE HEATMAP ORDERED BY DIVERSITY BUT COLORED BY TSNE CLUSTER OR GROUP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#           Change the c(1,2) in case we want clusters instead of groups.
create_zscan_heatmap <- function(seurat_object, cells_subset, cluster_column) {
  sce_zscan4_heatmap <- seurat_object[cells_subset, ]
  expression_matrix <- GetAssayData(seurat_object, layer = "data")
  
  # Order metadata by cluster
  sce_zscan4_heatmap@meta.data <- dfOrder(sce_zscan4_heatmap@meta.data, "Diversity") #cluster_column)
  barcode_cluster <- as.data.table(sce_zscan4_heatmap@meta.data[, c("Barcode", cluster_column)])
  group <- barcode_cluster[[cluster_column]]
  
  # Order the expression matrix
  expression <- expression_matrix[, rownames(sce_zscan4_heatmap@meta.data)]
  
  # Set up color palette
  # my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))
  # for group:
  # my.colors <- colorRampPalette(c("#FA8072", "#6C7B8B"))
  # for TSNE:
  my.colors <- colorRampPalette(c("#EEA2AD", "darkseagreen2"))
  
  
  color.df <- data.frame(COLOR_VALUE = group, barcode = barcode_cluster$Barcode)
  
  unique_colors <- my.colors(length(unique(group)))
  color.df$color.name <- unique_colors[match(color.df$COLOR_VALUE, unique(group))]
  
  # Filter color data frame for specific clusters and order the expression matrix
  color.df_filtered <- color.df[color.df$COLOR_VALUE %in% c(0,3), ]
  expression_filtered <- expression_matrix[, color.df_filtered$barcode]
  
  # Create binary expression matrix
  expression_binary <- expression_filtered
  expression_binary@x <- ifelse(expression_binary@x != 0, 1, 0)
  # expression_binary <- ifelse(expression_binary != 0, 1, 0)
  
  
  # Set up color palette for the heatmap
  pal <- colorRampPalette(c("cornsilk", "purple4"))(100)
  heatmap(as.matrix(expression), Colv = NA, Rowv = NA, scale = "row", 
          ColSideColors = color.df_filtered$color.name, col = pal, 
          main = paste("Cluster", paste(unique(color.df_filtered$COLOR_VALUE), collapse = " and "), "cells"),
          cexRow = 0.8, keep.dendro = FALSE, cexCol = 0.4, margins = c(5, 5))
}

# Example usage:
# Replace 'your_seurat_object', 'your_cells_subset', and 'YourClusterColumnName' with actual values
create_zscan_heatmap(sce_zscan4_groups, pluri, "TSNE_cluster")

create_zscan_heatmap(sce_zscan4_groups, pluri, "grupo")

create_zscan_heatmap(sce_zscan4_heatmap, pluri_genes, "TSNE_cluster")


sce_zscan4_heatmap <- sce_heatmap[pluri_genes,]



# ___________________________________________________________
#
#   Plot genes toti vs genes pluti
# ------------------------------------------------------------
# genes: toti and pluri

# Calcula las expresiones de los genes totipotentes y pluripotentes para cada célula
totipotent_expr <- FetchData(object = sce_12, vars = toti)
pluripotent_expr <- log(FetchData(object = sce_12, vars = pluri)+1)

# Crea un nuevo conjunto de datos con las expresiones de los genes totipotentes y pluripotentes
gene_expr_data <- data.frame(
  totipotent = rowMeans(totipotent_expr),
  pluripotent = rowMeans(pluripotent_expr),
  diversity = sce_12@meta.data[["Diversity"]],
  coverage = sce_12@meta.data[["Coverage"]],
  ncount = sce_12@meta.data[["nCount_RNA"]]
)

color_palette <- colorRampPalette(c("#eaac8b","#e56b6f", "#b56576", "#6d597a", "#355070"))(100)
color_palette <- colorRampPalette(c("#ffbf69","purple", "#355070"))(100)
palete <-brewer.pal(n = 3, name = "Set2")
palete <-brewer.pal(n = 3, name = "Reds")



# Crea un gráfico de dispersión con ggplot2
ggplot(gene_expr_data, aes(x = pluripotent, y = totipotent, color = diversity)) +
  geom_point(size = 6) +
  scale_color_gradientn(colours = palete, name = "Diversity") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes",
       title = "Expression of Totipotent vs. Pluripotent Genes")+
  theme_minimal()+
  ggtitle("Expression of Totipotent vs. Pluripotent Genes")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(gene_expr_data, aes(y = pluripotent, x = diversity, color = totipotent)) +
  geom_point(size=6) +
  scale_color_gradientn(colours = color_palette, name = "Expression of totipotent genes") +
  labs(x = "Diversity",
       y = "Pluripotent genes",
       title = "Diversity vs Pluripotency")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")

ggplot(gene_expr_data, aes(y = totipotent, x = diversity, color = pluripotent)) +
  geom_point(size=6) +
  scale_color_gradientn(colours = color_palette, name = "Expression of pluripotent genes") +
  labs(x = "Diversity",
       y = "Totipotent genes",
       title = "Diversity vs Totipotency")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")

ggplot(gene_expr_data, aes(x = diversity, y = coverage, color = totipotent)) +
  geom_point(size = 5) +
  scale_color_gradientn(colours = color_palette, name = "Expression of totipotent genes") +
  labs(x = "Diversity",
       y = "Coverage",
       title = "Coverage vs Diversity")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(gene_expr_data, aes(y = ncount, x = coverage, color = totipotent)) +
  geom_point() +
  scale_color_gradientn(colours = color_palette, name = "totipotent genes") +
  labs(y = "nCount_RNA",
       x = "Coverage (UMI)",
       title = "Coverage vs Pluripotent genes")+
  theme_minimal()

# SOLO CON LAS QUE ESTAN EN EL CLUSTER 0
# ----------------------------------------------------------------------------------------------------
sce_cluster0 <- subset(sce_12, subset = TSNE_cluster == 0)
# Calcula las expresiones de los genes totipotentes y pluripotentes para cada célula
totipotent_expr <- FetchData(object = sce_cluster0, vars = toti)
pluripotent_expr <- FetchData(object = sce_cluster0, vars = pluri)

# Crea un nuevo conjunto de datos con las expresiones de los genes totipotentes y pluripotentes
gene_expr_data <- data.frame(
  totipotent = rowMeans(totipotent_expr),
  pluripotent = rowMeans(pluripotent_expr),
  diversity = sce_cluster0@meta.data[["Diversity"]],
  coverage = sce_cluster0@meta.data[["Coverage"]]
)

palete <-brewer.pal(n = 3, name = "Set2")

color_palette <- colorRampPalette(c("#eaac8b","#e56b6f", "#b56576", "#6d597a", "#355070"))(100)
color_palette <- colorRampPalette(c("#ffbf69","purple", "#355070"))(100)


# Crea un gráfico de dispersión con ggplot2
ggplot(gene_expr_data, aes(x = pluripotent, y = totipotent, color = diversity)) +
  geom_point() +
  scale_color_gradientn(colours = palete, name = "Diversity") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes",
       title = "Expression of Totipotent vs. Pluripotent Genes")+
  theme_minimal()

ggplot(gene_expr_data, aes(y = diversity, x = pluripotent, color = totipotent)) +
  geom_point() +
  scale_color_gradientn(colours = color_palette, name = "Expression of totipotent genes") +
  labs(x = "Pluripotent Genes",
       y = "Diversity",
       title = "Diversity vs Pluripotency")+
  theme_minimal()

ggplot(gene_expr_data, aes(x = diversity, y = coverage, color = totipotent)) +
  geom_point() +
  scale_color_gradientn(colours = color_palette, name = "Expression of totipotent genes") +
  labs(x = "Diversity",
       y = "Coverage",
       title = "Coverage vs Diversity")+
  theme_minimal()

# ----------------------------------------------------
#              CON TODAS LAS CELULAS:
# ----------------------------------------------------
totipotent_expr <- FetchData(object = sce_zscan4_reads, vars = toti)
pluripotent_expr <- FetchData(object = sce_zscan4_reads, vars = pluri)

# Crea un nuevo conjunto de datos con las expresiones de los genes totipotentes y pluripotentes
gene_expr_data <- data.frame(
  totipotent = rowMeans(totipotent_expr),
  pluripotent = rowMeans(pluripotent_expr),
  diversity = sce_zscan4_reads@meta.data[["Diversity"]],
  zscan = sce_zscan4_reads@meta.data[["zscan4_expr_log"]],
  coverage = sce_zscan4_reads@meta.data[["Coverage"]],
  ncount = sce_zscan4_reads@meta.data[["nCount_RNA"]]
)

color_palette <- colorRampPalette(c("#eaac8b","#e56b6f", "#b56576", "#6d597a", "#355070"))(100)
color_palette <- colorRampPalette(c("#ffbf69","red", "#355070"))(100)
palete <-brewer.pal(n = 3, name = "Set2")
palete <-brewer.pal(n = 3, name = "Reds")

# Crea un gráfico de dispersión con ggplot2
ggplot(gene_expr_data, aes(x = pluripotent, y = totipotent, color = zscan)) +
  geom_point(size = 5) +
  scale_color_gradientn(colours = palete, name = "zscan4 levels") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes",
       title = "Expression of Totipotent vs. Pluripotent Genes")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(gene_expr_data, aes(x = pluripotent, y = totipotent, color = diversity)) +
  geom_point(size = 5) +
  scale_color_gradientn(colours = palete, name = "Diverstiy") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Total Expression of Totipotent Genes",
       title = "Expression of Totipotent vs. Pluripotent Genes colored by Diversity")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(gene_expr_data, aes(x = diversity, y = coverage, color = pluripotent)) +
  geom_point(size = 5) +
  scale_color_gradientn(colours = palete, name = "Expression of Pluri genes") +
  labs(x = "Diversity",
       y = "Coverage",
       title = "Coverage vs Diversity")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )


ggplot(gene_expr_data, aes(x = pluripotent, y = coverage, color = diversity)) +
  geom_point(size = 5) +
  scale_color_gradientn(colours = color_palette, name = "Diversity (spacer)") +
  labs(x = "Total Expression of Pluripotent Genes",
       y = "Coverage",
       title = "Coverage vs Pluripotent genes")+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )


ggplot(gene_expr_data, aes(y = diversity, x = coverage, color = totipotent)) +
  geom_point() +
  scale_color_gradientn(colours = color_palette, name = "totipotent") +
  labs(y = "diversity",
       x = "Coverage (UMI)",
       title = "Coverage vs Pluripotent genes")+
  theme_minimal()
# .....................................................................................................
# .....................................................................................................
# .....................................................................................................


# REVISARMOS EL COVERGA DE TODOS LOS GENES DE LAS CELULAS EN EL SCE ANTES DE NORMALIZAR !!

# sce.data <- Read10X(data.dir = "yolanda/Documentos/TFM/zscan4/filtered_feature_bc_matrix/")
# colnames(sce.data) <- gsub("-1", "", colnames(sce.data))
# sce <- CreateSeuratObject(counts = sce.data, project = "Zscan4", min.cells = 3, min.features = 200)

# DESPUES DE FILTROS Y TAL ---> sce_filtered


rownames(sce_zscan4_reads@meta.data)
colnames_to_keep <- intersect(colnames(sce_zscan4_reads), colnames(sce_filtered))
length(colnames_to_keep) # 1696

sce_cells <- sce_filtered[, colnames_to_keep]
dim(sce_cells)

View(as.matrix(sce_cells@assays$RNA$counts))

# Acceder a la matriz de expresión
expression_matrix <- as.matrix(sce_cells@assays$RNA$counts)

# Calcular el recuento total de genes en cada célula sumando todas las filas por columna
total_genes_per_cell <- colSums(expression_matrix)

# Imprimir los resultados
print(total_genes_per_cell)


sce_zscan4_reads@meta.data$nCounts <- total_genes_per_cell

View(sce_zscan4_reads@meta.data)












# ----------------------------------------------
#     COEFICIENTE DE CORRELACION DIVERSITY --- GEN                                
# ----------------------------------------------
# ----------------------------------------------

# Suponiendo que tienes tus vectores de expresión de genes y diversidad
expresion_genes <- list(
  Gen1 = c(10, 7, 2, 4, 3, 1),
  Gen2 = c(12, 18, 24, 30, 36, 42)
  # Agrega los vectores para otros genes aquí
)
# Vector de diversidad
diversidad <- c(0.2, 0.3, 0.5, 0.4, 0.6, 0.7)
sce_12@meta.data$Diversity

sce_cor <- sce_12[all,]
dim(sce_cor) #  1291   37


# 1) Extract vector of expression of each gene (37 cells)
# ----------------------------------------------
matrix_expresion <- as.matrix(log(sce_cor@assays$RNA$data+1))
order(colnames(matrix_expresion)) == order(colnames(sce_cor)) # verify the order is the same 
order(colnames(matrix_expresion)) == order(rownames(sce_cor@meta.data)) # verify the order is the same 

genes <- rownames(matrix_expresion)

vectores_expresion <- list() # create a list
for (gen in genes) {
  vector_gen <- matrix_expresion[gen, , drop = TRUE]  # drop = TRUE to mantain the vector structure
  vectores_expresion[[gen]] <- vector_gen
}

#example:
vectores_expresion$Zscan4d

# 2) vector of diversity of the 37 cells:
# -------------------------------------------------
diversidad <- sce_cor@meta.data$Diversity

# 3) Perform the correlation and save each result in a data frame
# -----------------------------------------------------------------
correlations <- sapply(vectores_expresion, function(gen_expresion) {
  cor(gen_expresion, diversidad, method = "spearman")
})
View(correlations)
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

View(df_correlaciones)
head(df_correlations)

# 4) see if the genes above 0.5 are pluri and bellow -0.5 are toti
# ------------------------------------------------------------------
neg <- subset(df_correlations, Spearman_Correlation < -0.6)$Gen
neg %in% toti

# "Gtf2b"   "Zscan4d" "Gm8300"  "Gm5662"  "Haus4"   "Ankrd22" "GFP"    

pos <- subset(df_correlations, Spearman_Correlation > 0.6)$Gen
pos %in% pluri
# "Sbf2"    "Cyb5b"   "Arid5b"  "Pitpnc1" "Arhgdia" "Trerf1" 

cat(" Genes totipotency with corr < -0.6 -->", neg, "\n",
    "Genes pluripotency with corr > 0.6 -->", pos, "\n")

pos_neg <- subset(df_correlations, Spearman_Correlation > 0.6 | Spearman_Correlation < -0.6)$Gen

# 5) PLOT SOME OF THESE GENES
# --------------------------------
all.genes<-rownames(sce_cor)
genes_df <- data.frame(Gene = all.genes)

zscan4d <- all.genes[grepl("^Zscan4d", all.genes)]

zscan4d_expression <- log(GetAssayData(object = sce_cor, slot = "data")[zscan4d, ] +1)  # Reemplaza 'gene_id' por el ID de Zscan4f


# add it to the metadata
sce_cor@meta.data$zscan4d <- zscan4d_expression # si es la suma de todos las las celulas con la expresion de cualquier zscan4 gene

#       ---->>>>>>> colored by the gene 

totipotent_expr <- log(FetchData(object = sce_cor, vars = toti)+1)
pluripotent_expr <- log(FetchData(object = sce_cor, vars = pluri)+1)

# Crea un nuevo conjunto de datos con las expresiones de los genes totipotentes y pluripotentes
gene_expr_data_cor <- data.frame(
  totipotent = rowMeans(totipotent_expr),
  pluripotent = rowMeans(pluripotent_expr),
  diversity = sce_cor@meta.data[["Diversity"]],
  coverage = sce_cor@meta.data[["Coverage"]],
  ncount = sce_cor@meta.data[["nCount_RNA"]],
  zscan4d = sce_cor@meta.data[["zscan4d"]]
)

color_palette <- colorRampPalette(c("#eaac8b","#e56b6f", "#b56576", "#6d597a", "#355070"))(100)
color_palette <- colorRampPalette(c("#ffbf69","purple", "#355070"))(100)
palete <-brewer.pal(n = 3, name = "Set2")
palete <-brewer.pal(n = 3, name = "Reds")



# Crea un gráfico de dispersión con ggplot2
ggplot(gene_expr_data_cor, aes(x = diversity, y = totipotent, color = zscan4d)) +
  geom_point(size = 6) +
  scale_color_gradientn(colours = palete, name = "Diversity") +
  labs(x = "Diversity",
       y = "log Zscan4d expression",
       title = "Zscan4d expression vs diversity")+
  theme_minimal()+
  ggtitle("Expression of Totipotent vs. Pluripotent Genes")+
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )

ggplot(gene_expr_data_cor, aes(x = diversity, y = zscan4d)) +
  geom_point(size = 6, color = "#b56576") +  # Establecer un color específico para todos los puntos
  scale_color_gradientn(colours = palete, name = "Diversity") +  # Si deseas mantener la leyenda, puedes dejar esta línea
  labs(x = "Diversity",
       y = "log Zscan4d expression",
       title = "Zscan4d expression vs diversity") +
  theme_minimal() +
  ggtitle("Zscan4d expression vs Diversity") +
  theme(
    plot.title = element_text(hjust = 0.5)  # Esto centra el título horizontalmente
  )+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")



library(ggplot2)
library(viridis)

# Función para visualizar la expresión de un gen en relación con la diversidad
plot_gene_diversity <- function(gen_name, sce_data, color_palette) {
  # Obtener la expresión del gen especificado
  gen_expr <- log(GetAssayData(object = sce_data, slot = "data")[gen_name, ] + 1)
  
  gene_expr_data <- data.frame(
    expression = gen_expr,
    diversity = sce_data@meta.data[["Diversity"]]
  )
  
  # ggplot(gene_expr_data, aes(x = diversity, y = expression)) +
  #   geom_point(size = 6, aes(color = expression)) +
  #   scale_color_viridis(option = "mako", name = "Gene Expression") +
  #   labs(x = "Diversity",
  #        y = paste("log", gen_name, "expression"),
  #        title = paste(gen_name, "Expression vs. Diversity")) +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5)
  #   )
  
  ggplot(gene_expr_data, aes(x = diversity, y = expression)) +
    geom_point(size = 6, color = "#b56576") +
    labs(x = "Diversity",
         y = paste("log", gen_name, "expression"),
         title = paste(gen_name, "Expression vs. Diversity")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)+
        geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed")
    )
}

# Uso de la función para visualizar la expresión del gen Zscan4d en relación con la diversidad
plot_gene_diversity("Trerf1", sce_cor, )



##########################################
#               HEATMAP 
##########################################


sce_cor_heatmap<- sce_cor[pos_neg,]
dim(sce_cor_heatmap) #  13   37
expression_matrix<- log(GetAssayData(sce_cor_heatmap, slot = "data")+1)
expression_matrix_ordered <- expression_matrix[match(pos_neg, rownames(expression_matrix)), ]

# _____________Diversity ________________________________________________________________
sce_cor_heatmap@meta.data <- dfOrder(sce_cor_heatmap@meta.data, "Diversity")
barcode.diversity <- as.data.table(sce_cor_heatmap@meta.data[,c("Barcode","Diversity")])
diversity <- barcode.diversity$Diversity

barcode.grupo <- as.data.table(sce_cor_heatmap@meta.data[,c("Barcode","TSNE_cluster")])
group <- sce_cor_heatmap$TSNE_cluster

color.df <- data.frame(
  COLOR_VALUE = diversity,
  barcode = barcode.diversity$Barcode
)
# color.df$Rangos <- cut(color.df$COLOR_VALUE, breaks = c(rangos, Inf), labels = FALSE, include.lowest = TRUE)
my.colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))
# color.df<-data.frame(COLOR_VALUE=diversity, color.name=my.colors(length(diversity)))or.df$Rangos]
color.df$color.name <- my.colors(length(diversity))

expression <- expression_matrix_ordered[ ,color.df$barcode]

#               FOR A HEATMAP WITH GENE ANNOTATION:
#               ........................................................................
color_values <- colorRamp(my_colors)(100)
color_vector <- color_values[findInterval(color.df$COLOR_VALUE, seq(min(color.df$COLOR_VALUE), max(color.df$COLOR_VALUE), length.out = 100))]

ha_diversity <- HeatmapAnnotation(Diversity = diversity,
                                  cluster = group)

expression <- expression[,order(color.df$COLOR_VALUE)]
colnames(expression)

ha_diversity <- HeatmapAnnotation(Diversity = color.df$COLOR_VALUE)
col_names <- color.df$color.name


names(col_names) <- color.df$COLOR_VALUE
draw(ha_diversity)
genes_totipotent <- intersect(rownames(expression), neg)
genes_pluripotent <- intersect(rownames(expression), pos)

row_split <- c(rep("Pluripotent", length(genes_pluripotent)), rep("Totipotent", length(genes_totipotent)))

grupos <- list(Pluripotent = genes_pluripotent, Totipotent = genes_totipotent)

text <- grupos
class(expression)
my_colors <- colorRampPalette(c("lightskyblue1", "hotpink4"))(100)



pdf("Documents/heatmap_PRUEBA_4.pdf")
Heatmap(as.matrix(expression), name = "Expression", cluster_rows = FALSE, cluster_columns = FALSE,
        row_split = row_split,
        show_row_names = TRUE,
        color_space = "LAB",
        column_names_gp = gpar(fontsize = 6),
        top_annotation = ha_diversity)


# right_annotation = rowAnnotation(textbox = anno_textbox(row_split, text)))

dev.off()






#                                      ---------------------------------------
#                                               PEARSON CORREALTION 
#                                             --------------------------

sce_cor <- sce_12[all,]
dim(sce_cor) #  1291   37


# 1) Extract vector of expression of each gene (37 cells)
# ----------------------------------------------
matrix_expresion <- as.matrix(log(sce_cor@assays$RNA$data+1))
order(colnames(matrix_expresion)) == order(colnames(sce_cor)) # verify the order is the same 
order(colnames(matrix_expresion)) == order(rownames(sce_cor@meta.data)) # verify the order is the same 

genes <- rownames(matrix_expresion)

vectores_expresion <- list() # create a list
for (gen in genes) {
  vector_gen <- matrix_expresion[gen, , drop = TRUE]  # drop = TRUE to mantain the vector structure
  vectores_expresion[[gen]] <- vector_gen
}

#example:
vectores_expresion$Zscan4d

# 2) vector of diversity of the 37 cells:
# -------------------------------------------------
diversidad <- sce_cor@meta.data$Diversity

# 3) Perform the correlation and save each result in a data frame
# -----------------------------------------------------------------
correlations <- sapply(vectores_expresion, function(gen_expresion) {
  cor(gen_expresion, diversidad, method = "pearson")
})
# Warning messages:
#   1: In cor(gen_expresion, diversidad, method = "spearman") :
#   the standard deviation is zero
# 2: In cor(gen_expresion, diversidad, method = "spearman") :
#   the standard deviation is zero

df_correlations_2 <- data.frame(
  Gen = names(correlations),
  Pearson_correlation = correlations,
  stringsAsFactors = FALSE
)

View(df_correlations_2)
head(df_correlations_2)

# 4) see if the genes above 0.5 are pluri and bellow -0.5 are toti
# ------------------------------------------------------------------
neg <- subset(df_correlations, Spearman_Correlation < -0.6)$Gen
neg %in% toti

# "Gtf2b"   "Zscan4d" "Gm8300"  "Gm5662"  "Haus4"   "Ankrd22" "GFP"    

pos <- subset(df_correlations, Spearman_Correlation > 0.6)$Gen
pos %in% pluri
# "Sbf2"    "Cyb5b"   "Arid5b"  "Pitpnc1" "Arhgdia" "Trerf1" 

cat(" Genes totipotency with corr < -0.6 -->", neg, "\n",
    "Genes pluripotency with corr > 0.6 -->", pos, "\n")

pos_neg <- subset(df_correlations, Spearman_Correlation > 0.6 | Spearman_Correlation < -0.6)$Gen

# 5) PLOT SOME OF THESE GENES
# --------------------------------





# NO SE SI ESTO DE ABAJO SE BORRARIA O QUE 
# Obtener el orden de las células (suponiendo que tienes el mismo orden en tus vectores de expresión de genes)
orden_celulas <- seq_along(diversidad)

# Ordenar los vectores de expresión de genes y diversidad según el orden de las células
expresion_genes_ordenado <- lapply(expresion_genes, function(gen_exp) gen_exp[orden_celulas])
diversidad_ordenada <- diversidad[orden_celulas]

# Calcular la correlación de Spearman entre cada gen y la diversidad
resultados_correlacion <- sapply(expresion_genes_ordenado, function(gen_exp) cor(gen_exp, diversidad_ordenada, method = "spearman"))

# Crear un data frame con los resultados
resultados_df <- data.frame(
  Gen = names(resultados_correlacion),
  Correlacion_Spearman = resultados_correlacion,
  row.names = NULL
)

# Imprimir el data frame con los resultados
print(resultados_df)




















