#' imprints_phoQP_ksea
#'
#' Function to put the analysis tab in a summary format suitable to use the \link[https://casecpb.shinyapps.io/ksea/]{KSEA app}
#'
#' @param data The analysis tab output from the imprints_phoQP_hit_peptide
#'   (either a data frame or a path to the file)
#' @param save Logical to tell if you want to save the ksea file in the csv format for each treatment
#'
#' @return A data frame in the KSEA format
#'
#' @export

imprints_phoQP_ksea <- function(data, save = FALSE){
  if(inherits(data, "character") & length(data) == 1){
    if(grepl("\\.txt$", data)){
      data <- readr::read_tsv(data)
    }
    else if(grepl("\\.csv$", data)){
      readr::read_csv(data)
    }
    else if(grepl("\\.xlsx$", data)){
      data <- openxlsx::read.xlsx(data)
    }
  }
  else if(!inherits(data, "data.frame")){
    stop("Input data can only be a data.frame or a path to an existing file")
  }

  cn_needed <- c("Gene", "description", "Master.Protein.Accessions",
                 "Positions.in.Master.Proteins", "Sequence", "Annotated.Sequence", "Modifications")
  if(!all(cn_needed %in% colnames(data))){
    cn_missing <- cn_needed[!(cn_needed %in% colnames(data))]
    stop(paste("The", ifelse(length(cn_missing) > 1, "columns", "column"),
               paste(cn_missing, collapse = ", "),
               ifelse(length(cn_missing) > 1, "are", "is"),
               "missing !")
         )
  }

  check_columns <- colnames(data)[!(colnames(data) %in% cn_needed)]
  if(sum(sapply(strsplit(check_columns, "_"), length) == 2) < 2){ # at least one treatment giving p-value + fold-change
    stop("Your data doesn't follow the format from the 'analysis_tab' output from imprints_phoQP_hit_peptide function")
  }
  check_columns <- unique(sapply(strsplit(check_columns, "_"), "[[", 1))
  if(!("pval" %in% check_columns) | !any(grepl("^\\d{2}C", check_columns))){
    stop("Your data doesn't follow the format from the 'analysis_tab' output from imprints_phoQP_hit_peptide function")
  }
  FC_column <- grep("^\\d{2}C", check_columns, value = TRUE)

  data <- data[,c(cn_needed,
                  grep(paste0(FC_column, "_"), colnames(data), value = TRUE),
                  grep("pval_", colnames(data), value = TRUE))]

  # reading
  ksea <- data %>%
    tidyr::gather("key", "value",
                  -Gene, -description, -Master.Protein.Accessions, -Positions.in.Master.Proteins,
                  -Sequence, -Annotated.Sequence, -Modifications) %>%
    tidyr::separate(key, into = c("key", "treatment"), sep = "_") %>%
    tidyr::spread(key, value)

  # extract phosphorylation site (only keeping site with position)
  ksea$Residue.Both <- gsub(".*Phospho \\[|\\]$|\\(.{,4})", "", ksea$Modifications)
  ksea$Residue.Both <- gsub("[A-Z]\\D{1,}|[A-Z]$", "", ksea$Residue.Both)
  ksea$Residue.Both <- gsub(";$|; $| ", "", ksea$Residue.Both)

  # selecting needed columns and renaming
  ksea <- ksea[,c("Master.Protein.Accessions", "Gene", "Sequence", "Residue.Both", FC_column, "pval", "treatment")]
  colnames(ksea)[c(1,3,5,6)] <- c("Protein","Peptide", "FC", "p")

  # taking non log2 tranformed fold change
  ksea$FC <- 2**ksea$FC
  # if any NA, row is discared anyway
  ksea <- ksea[!apply(ksea, 1, function(x) any(is.na(x))),]

  if(save){
    for(i in unique(ksea$treatment)){
      x <- ksea %>%
        dplyr::filter(treatment == i) %>%
        select(-treatment)

      write.table(x, paste0(format(Sys.time(), "%y%m%d_%H%M_"), "ksea_", i, ".csv"),
                  row.names = FALSE,
                  quote = FALSE, sep = ",", qmethod = "double")

    }
  }

  return(ksea)
}


