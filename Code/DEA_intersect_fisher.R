##################################--> TFM- Yolanda Andrés

#           -------------------------------------------
#             OVERLAPPING DEA GENES AND FU ET AL 2020 !
#           -------------------------------------------

t1_subset <- read.csv("table_s3")
all_genes <- read.csv("table_s2")
# 1. intersect
# ---------------
genes_to_keep<- intersect(t1_subset$Transcript_Name, genes_df$Gene) 
t1_subset_1 <- t1_subset[t1_subset$Transcript_Name %in% genes_to_keep, ]


# TOTI
# ------------
t1_subset_up <- subset(t1_subset_1, isDE == "Up") 
t1_subset_up <- subset(t1_subset_up, logFC >= 2) 
genestoti <- intersect(t1_subset_up$Transcript_Name, toti_genes)


# PLURI
# ----------
t1_subset_down <- subset(t1_subset_1, isDE == "Down")
t1_subset_down <- subset(t1_subset_down, logFC < -1.1)
genespluri <- intersect(t1_subset_down$Transcript_Name, pluri_genes) 


# EULER DIAGRAM
#..............................
overlap_toti <-length(t)
data_toti <- c(
  `Cluster 2 markers` = length(toti_genes),
  `2C embryo markers \n Fu et al. (2020)` = length(t1_subset_up$Transcript_Name),
  `Cluster 2 markers&2C embryo markers \n Fu et al. (2020)` = overlap_toti)

overlap_pluri <- length(genespluri)
data_pluri <- c(
  `Cluster 0 \n markers` = length(pluri_genes),
  `Pluripotent markers \n Fu et al. (2020)` = length(t1_subset_down$Transcript_Name),
  `Cluster 0 \n markers&Pluripotent markers \n Fu et al. (2020)` = overlap_pluri)

fit_toti <- euler(data_toti)
fit_pluri <- euler(data_pluri)


plot(fit_toti, fills = c("#fa9fb5", "gray89"), edges = list(col = c("#fa9fb9", "gray3")), 
     labels = list(font = 2, fontsize= 25), quantities = list(fontsize = 30))
grid.text("Overlap between Toti Genes and Up-regulated Genes", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

plot(fit_pluri, fills = c("#bcbddc", "gray89"), edges = list(col = c("#bcbddc", "gray3")), 
     labels = list(font = 2, fontsize= 30), quantities = list(fontsize = 35))
grid.text("Overlap between Toti Genes and Up-regulated Genes", x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))




#   
#           -------------------------------------------
#                         PERMUTATION TEST
#           -------------------------------------------

#### TEST------------> TOTI
# ......................
all_genes <- genes_df$Gene

set.seed(123)
n_perm <- 10000
overlap_counts <- numeric(n_perm)

# permutation test:
for (i in 1:n_perm) {
  sampled_genes <- sample(all_genes, length(toti_genes))
  overlap <- length(intersect(sampled_genes, genes_2c_permu))
  overlap_counts[i] <- overlap
}

mean_overlap <- mean(overlap_counts)
hist(overlap_counts, main = "Permutation Test of Overlap Cluster 2 markers \n with Totipotent markers", xlab = "Number of Overlapping Genes", breaks = 20)
abline(v = mean_overlap, col = "red", lwd = 2)
obs_overlap <- length(intersect(toti_genes, t1_subset_1$Transcript_Name))  # 53 genes comparten
abline(v = obs_overlap, col = "blue", lwd = 2)

p_value <- mean(overlap_counts >= obs_overlap)





#### TEST ------------> PLURI
# ......................

set.seed(123)
genes_2c_permu <- t1_subset_1$Transcript_Name
n_perm <- 10000
overlap_counts <- numeric(n_perm)


for (i in 1:n_perm) {
  sampled_genes <- sample(all, length(pluri_genes))
  overlap <- length(intersect(sampled_genes, genes_2c_permu))
  overlap_counts[i] <- overlap
}

mean_overlap <- mean(overlap_counts)
hist(overlap_counts, main = "Permutation Test of Overlap Cluster 2 markers \n with Totipotent markers", xlab = "Number of Overlapping Genes", breaks = 20)
abline(v = mean_overlap, col = "red", lwd = 2)
obs_overlap <- length(intersect(pluri_genes, t1_subset_1$Transcript_Name))  # 53 genes comparten
abline(v = obs_overlap, col = "blue", lwd = 2)

p_value <- mean(overlap_counts >= obs_overlap)


#   
#           -------------------------------------------
#                            FISHER TEST 
#           -------------------------------------------


genest <- toti_genes
subset_genes <- t1_subset_up$Transcript_Name

genes_total <- unique(c(toti_genes, subset_genes))

genes_total <- all_genes
# genest <- pluri_genes
in_toti_and_subset <- length(intersect(genest, subset_genes))
in_toti_not_subset <- length(setdiff(genest, subset_genes))
not_in_toti_in_subset <- length(setdiff(subset_genes, genest))
not_in_toti_not_subset <- length(genes_total) - in_toti_and_subset - in_toti_not_subset - not_in_toti_in_subset

contingency_table <- matrix(c(in_toti_and_subset, not_in_toti_in_subset,
                              in_toti_not_subset, not_in_toti_not_subset), 
                            nrow = 2, 
                            dimnames = list('Toti' = c('In', 'Not In'),
                                            'Subset' = c('In', 'Not In')))

fisher_result <- fisher.test(contingency_table)
print(fisher_result)

#CHI TEST
chi_squared_result <- chisq.test(contingency_table)
print(chi_squared_result)
