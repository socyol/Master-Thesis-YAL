########################

# Function to extract information about CRISPR-Barcode features
# .............................................................

read_table_function <- function(row) {
  
  sfl1.1 <- "ACACC"
  sfl1.2 <- "AACACC"
  sfl1.3 <- "AAACACC"
  
  sfl2.1 <- "TTAGA"
  sfl2.2 <- "TTAGAG"
  sfl2.3 <- "TTAGAGC"
  
  cassette <- row[["Cassette"]]
  
  # 1- FIND POSITIONS OF THE FLANKING REGIONS
  positions_sfl1 <- gregexpr(sfl1.1, cassette)[[1]]
  count_sfl1 <- length(positions_sfl1)
  
  if (count_sfl1 != 1) {
    positions_sfl1_2 <- gregexpr(sfl1.2, cassette)[[1]]
    count_sfl1_2 <- length(positions_sfl1_2)
    
    if (count_sfl1_2 == 1) {
      row[["POSICION_SFL1"]] <- positions_sfl1_2
      row[["SFL1"]] <- 2
      sfl1 <- sfl1.2
      
    } else {
      
      positions_sfl1_3 <- gregexpr(sfl1.3, cassette)[[1]]
      count_sfl1_3 <- length(positions_sfl1_3)
      
      if (count_sfl1_3 == 1) {
        row[["POSICION_SFL1"]] <- positions_sfl1_3
        row[["SFL1"]] <- 3
        sfl1 <- sfl1.3
      }
    }
  } else {
    row[["POSICION_SFL1"]] <- positions_sfl1
    row[["SFL1"]] <- 1
    sfl1 <- sfl1.1
  }
  positions_sfl2 <- gregexpr(sfl2.1, cassette)[[1]]
  count_sfl2 <- length(positions_sfl2)
  
  if (count_sfl2 != 1) {
    positions_sfl2_2 <- gregexpr(sfl2.2, cassette)[[1]]
    count_sfl2_2 <- length(positions_sfl2_2)
    
    if (count_sfl2_2 == 1) {
      row[["POSICION_SFL2"]] <- positions_sfl2_2
      row[["SFL2"]] <- 2
      sfl2 <- sfl2.2
      
    } else {
      
      positions_sfl2_3 <- gregexpr(sfl2.3, cassette)[[1]]
      count_sfl2_3 <- length(positions_sfl2_3)
      
      if (count_sfl2_3 == 1) {
        row[["POSICION_SFL2"]] <- positions_sfl2_3
        row[["SFL2"]] <- 3
        sfl2 <- sfl2.3
      }
    }
  } else {
    row[["POSICION_SFL2"]] <- positions_sfl2
    row[["SFL2"]] <- 1
    sfl2 <- sfl2.1
  }
  
  if (!("POSICION_SFL1" %in% names(row)) || !("POSICION_SFL2" %in% names(row)) ||
      is.null(row$POSICION_SFL1) || is.null(row$POSICION_SFL2) ||
      length(row$POSICION_SFL1) == 0 || length(row$POSICION_SFL2) == 0) {
    return(NULL)  # Retorna NULL para eliminar la row
  }
  
  row[["Match"]] <- gsub(paste0(".*(?:", sfl1, ")(.*?)(?:", sfl2, ").*"), "\\1", cassette)
  
  return(row)
}

