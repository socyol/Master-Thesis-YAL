## -------------------------------------------------------
#     
#             FILTERING UMIS - (I)
# -------------------------------------------------------
## ------------------------------------------------------


# Libraries

library(dplyr)

library(ggplot2)


# FUNCTIONS
#~~~~~~~~~~~~~~

calcular_porcentaje_perdida <- function(t1, t2, barcode_col, umi_col) {
  t1_grouped <- t1 %>% group_by(!!sym(barcode_col))
  t2_grouped <- t2 %>% group_by(!!sym(barcode_col))
  unique_UMI_t1 <- t1_grouped %>% summarize(unique_UMI = n_distinct(!!sym(umi_col)))
  unique_UMI_t2 <- t2_grouped %>% summarize(unique_UMI = n_distinct(!!sym(umi_col)))
  
  result <- left_join(unique_UMI_t1, unique_UMI_t2, by = barcode_col) %>%
    mutate(percentage_lost = ifelse(is.na(unique_UMI.y), 100, (1 - unique_UMI.y / unique_UMI.x) * 100))
  

}


print_porcentaje_perdida <- function(t1, t2, barcode_col, umi_col) {

  result <- calcular_porcentaje_perdida(t1, t2, barcode_col, umi_col)
  
}

# 1) load
#-----------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
nombre_archivo <- args[1]


output <- read.table(yourpath, sep = "\t")  

colnames(output)[colnames(output) == "V1"] <- "Barcode"
colnames(output)[colnames(output) == "V2"] <- "UMI"
colnames(output)[colnames(output) == "V3"] <- "Cassette"
colnames(output)[colnames(output) == "V4"] <- "Quality"


# 2. Filters
#-----------------------------------------------------------------------------------
# 2.1 DELETE those without the fln regions and with N in the sequence
# ..................................................................
# Original sequence spacer: GGTATGCGGATGCAATCTCCG
sfl1 <- "ACACC" # flanquing regions
sfl2 <- "TTAGA"
output_filtered <- output[grepl(paste0(sfl1, ".*", sfl2), output$Cassette), ] # to filter those without the flanquing reagions
output_filtered <- output_filtered[!grepl("N", output_filtered$Cassette), ] # filter those with "N" in the sequencing column
# 336200      4

# Statistics whole dataset
n_reads <- nrow(output) # Number of cells before filtering
n_filt_reads <- nrow(output_filtered) # Number of resulted cells after filtering
perc_filtered <- (n_filt_reads/n_reads)*100
length(unique(output$Barcode))

UMI_before <- length(unique(output$UMI))
UMI_after <- length(unique(output_filtered$UMI))
perc_UMI <- (1-UMI_after/UMI_before)*100
cat("\n", "AFTER FILTERING FLN REGIONS:", "\n",
    "________________________________________________", "\n",
    "-------------------------------------------------", "\n",
    "number of reads before filtering =", n_reads, "\n", 
    "number of reads after filtering = ", n_filt_reads, "\n",
    "Percentage of final reads:", perc_filtered, "%\n",
    "......................................................", "\n",
    "......................................................", "\n",
    "number of UMIs before filtering =", length(unique(output$UMI)), "\n", 
    "number of UMIs after filtering = ", length(unique(output_filtered$UMI)), "\n",
    "Percentage LOST of UMIS:", perc_UMI, "%\n")



# ---------------------------
#   create a TABLE TO SAVE THE RESULTS 
# --------------------------------------------------
library(data.table)
mi_tabla <- data.table(
  STEP = character(),
  Reads = numeric(),
  UMI = numeric(),
  perc_lost_UMI=numeric(),
  Cells = numeric()
)

# NO FILTERS
result1 <- data.table(
  STEP = "NO filters",
  Reads = n_reads,
  UMI = UMI_before,
  Cells = length(unique(output$Barcode)),
  perc_lost_UMI = 0
)

mi_tabla <- rbind(mi_tabla, result1)

# FILTER FLN
result2 <- data.table(
  STEP = "FLN Filter",
  Reads = n_filt_reads,
  UMI = UMI_after,
  Cells = length(unique(output_filtered$Barcode)),
  perc_lost_UMI = perc_UMI
)

mi_tabla <- rbind(mi_tabla, result2)

# ---------------------------------------------------------



# 2.2. Delete the UMI with less than 3 reads
# ..................................................................
# First we do a new column with Barcode and UMI pasted by a "-":
# Suponiendo que tu tabla filtrada se llama output_filtered
output_filtered <- output_filtered %>%
  mutate(Barcode_UMI = paste(Barcode, UMI, sep = "-"))

output_filtered2 <- output_filtered %>%
  group_by(Barcode_UMI) %>%
  filter(n() >= 3) %>%
  ungroup()

st <- print_porcentaje_perdida(output_filtered, output_filtered2, "Barcode", "UMI")


# reads statistics
n_reads_2 <- nrow(output_filtered) # Number of cells before filtering
n_filt_reads_2 <- nrow(output_filtered2) # Number of resulted cells after filtering
perc_filtered_2 <- (n_filt_reads_2/n_reads_2)*100

UMI_before_2 <- length(unique(output_filtered$Barcode_UMI))
UMI_after_2 <- length(unique(output_filtered2$Barcode_UMI))
perc_UMI_2 <- (1 - UMI_after_2/UMI_before_2)*100
cat("\n", "________________________________________________", "\n", 
    "FILTERING OF >3 READS PER UMI:", "\n",
    "________________________________________________", "\n",
    "number of UMIs before filtering =", UMI_before_2, "\n", 
    "number of UMIs after filtering = ", UMI_after_2, "\n",
    "Percentage LOST of UMIS:", perc_UMI_2, "%\n",
    "......................................................", "\n",
    "MEAN of LOST of UMIs per cell =", mean(st$percentage_lost), "\n")

result3 <- data.table(
  STEP = ">3 reads Filter",
  Reads = n_filt_reads_2,
  UMI = UMI_after_2,
  Cells = length(unique(output_filtered2$Barcode)),
  perc_lost_UMI = perc_UMI_2
)

mi_tabla <- rbind(mi_tabla, result3)



# 3. Apply the function OF READS
#-----------------------------------------------------------------------------------

source("yourpath/function-readtable.R")

result_list <- lapply(1:nrow(output_filtered2), function(i) read_table_function(output_filtered2[i, ]))
read_table <- do.call(rbind, result_list)

args <- commandArgs(trailingOnly = TRUE)
nombre_archivo <- args[1]

numero_archivo <- sub("subset_", "", sub("_output.txt", "", basename(nombre_archivo)))

ruta_salida <- "yourpath" 

nombre_salida <- paste0(ruta_salida, "OUTPUT_", numero_archivo, ".csv")
write.csv(read_table, file = nombre_salida)

# -------------------------------------------------------------------------------------

# PROCESSING reads table:

data <- read_table
data_filtered <- data %>%
  mutate(Barcode_UMI_read = paste(Barcode, UMI, Match, sep = "-"))


output_summary <- data_filtered %>%
  group_by(Barcode_UMI_read, Barcode_UMI) %>%
  summarise(repetition_count = n()) %>%
  ungroup()

total_counts <- output_summary %>%
  group_by(Barcode_UMI) %>%
  summarise(total_count = sum(repetition_count))

output_summary_with_total <- left_join(output_summary, total_counts, by = "Barcode_UMI")

output_summary_with_percentage <- output_summary_with_total %>%
  mutate(percentage = repetition_count / total_count * 100)


# 5.Merge and keep just 1 UMI-read
#-----------------------------------------------------------------------------------
data_filtered2 <- merge(data_filtered, output_summary_with_percentage, by="Barcode_UMI_read")

data_filtered2 <- na.omit(data_filtered2)


# Keep UNIQUES
data_1 <- data_filtered2 %>%
  distinct(Barcode_UMI_read, .keep_all = TRUE)


# Filtrar las filas que tienen menos de un 50% en percentage
data_2 <- subset(data_1, percentage > 50)

st2<- print_porcentaje_perdida(data_1, data_2, "Barcode", "UMI")


# AFTER KEEPING UNIQUES: ................................
# reads statistics
n_reads_3 <- nrow(data_1) # Number of cells before filtering
n_filt_reads_3 <- nrow(data_2) # Number of resulted cells after filtering
perc_filtered_3 <- (n_filt_reads_3/n_reads_3)*100

UMI_before_3 <- length(unique(data_1$Barcode_UMI.x))
UMI_after_3 <- length(unique(data_2$Barcode_UMI.x))
perc_UMI_3 <- (1 - UMI_after_3/UMI_before_3)*100
cat("\n", "________________________________________________", "\n", 
    "ºFILTER OF CONSENSUS", "\n",
    "________________________________________________", "\n",
    "number of UMIs before filtering =", UMI_before_3, "\n", 
    "number of UMIs after filtering = ", UMI_after_3, "\n",
    "Percentage LOST of UMIS:", perc_UMI_3, "%\n",
    "......................................................", "\n",
    "MEAN of LOST of UMIs per cell =", mean(st2$percentage_lost), "\n")

result4 <- data.table(
  STEP = "Consensus and unique UMI Filter",
  Reads = n_filt_reads_3,
  UMI = UMI_after_3,
  Cells = length(unique(data_2$Barcode)),
  perc_lost_UMI = perc_UMI_3
)

mi_tabla <- rbind(mi_tabla, result4)

# SAVE THE TABLE WITH RESULTS
outputpath <- "yourpath"
nombre_salida <- paste0(outputpath, "TABLE_", numero_archivo, ".csv")
write.csv(mi_tabla, file = nombre_salida)

# SAVE THE OUTPUT
outputpath_2 <- "yourpath"
out_name <- paste0(outputpath_2, "FINAL_OUTPUT_", numero_archivo, ".csv")
write.csv(data_2, file = out_name)

