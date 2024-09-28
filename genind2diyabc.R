genind2diyabc <- function(genind_object, remove_intrapop_mono = FALSE, remove_intrapop_allmissing = FALSE, file_name = NULL) {
  remove_intrapop_monomorphic_loci <- function(data_tab, populations) {
    loci_to_remove <- logical(ncol(data_tab))
    for (pop in populations) {
      subset_data <- data_tab[substr(rownames(data_tab), 1, 2) == as.character(pop), , drop = FALSE]
      loci_to_remove <- loci_to_remove | sapply(subset_data, function(col) length(unique(col[col != 9])) < 2)
    }
    return(data_tab[, !loci_to_remove, drop = FALSE])
  }

  remove_intrapop_missing_loci <- function(data_tab, populations) {
    loci_to_remove <- logical(ncol(data_tab))
    for (pop in populations) {
      subset_data <- data_tab[substr(rownames(data_tab), 1, 2) == as.character(pop), , drop = FALSE]
      loci_to_remove <- loci_to_remove | sapply(subset_data, function(col) length(unique(col)) == 1)
    }
    return(data_tab[, !loci_to_remove, drop = FALSE])
  }

  prepare_diyabc_dataframe <- function(genind_object, data_tab) {
    ind_data <- data.frame(IND = rownames(data_tab), SEX = rep("9", nrow(data_tab)), POP = genind_object@pop)
    output_df <- cbind(ind_data, data_tab)
    colnames(output_df)[-(1:3)] <- rep("A", ncol(output_df) - 3)
    return(output_df)
  }

  write_diyabc_output <- function(df, populations, file_name) {
    if (is.null(file_name)) {
      file_name <- paste(gsub(", ", "_", toString(populations)), "_DIYABC.snp", sep = "")  
    }
    write("<NM=1.0NF> <MAF=hudson>", file_name, append = FALSE)
    write.table(df, file_name, sep = "\t", quote = FALSE, na = "9", row.names = FALSE, append = TRUE)
    print(sprintf("Output written to %s", file_name))
  }

  populations <- unique(genind_object@pop)
  data_tab <- as.data.frame(genind_object[genind_object@pop %in% populations]@tab)
  data_tab[is.na(data_tab)] <- 9
  num_loci_start <- ncol(data_tab)

  if (remove_intrapop_mono) {
    data_tab <- remove_intrapop_monomorphic_loci(data_tab, populations)
  }

  if (remove_intrapop_allmissing) {
    data_tab <- remove_intrapop_missing_loci(data_tab, populations)
  }

  num_loci_end <- ncol(data_tab)
  loci_removed <- num_loci_start - num_loci_end
  cat(sprintf("%d loci removed\n", loci_removed))

  df <- prepare_diyabc_dataframe(genind_object, data_tab)
  write_diyabc_output(df, populations, file_name)
}
