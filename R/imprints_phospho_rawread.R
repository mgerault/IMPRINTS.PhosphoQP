#' imprints_phospho_rawread
#'
#' Function to read the phosphoproteomics data from Proteome Discoverer
#'
#' @param phospho_file The path to the PeptidesGroups file from PD.
#' @param treatment The treatment corresponding to each fo your TMT channel. It should be
#'   the same length as your number of TMT channel and should follow the format 'Rep_treatment'
#'   like 'B1_Vehicle' for example.
#' @param prefixcontaminant Character corresponding to the prefix used to identify contaminants.
#' @param countthreshold The minimal threshold number of associated abundance count of peptides, default is 1.
#' @param dataset_name The name of your dataset
#'
#' @return Your data filtered and in the right format to use the IMPRINTS.CETSA package.
#'
#' @export

imprints_phospho_rawread <- function(phospho_file, treatment, prefixcontaminant = "",
                                     countthreshold = 1, dataset_name = "phospho"){
  if(countthreshold < 0){
    countthreshold <- 0
    messsage("Parameter countthreshold can only be positive; value set to 0, hence no filtering will be applied.")
  }

  # data reading
  if(inherits(phospho_file, "character")){
    stopifnot(file.exists(phospho_file))

    if(grepl("\\.txt$", phospho_file)){
      phospho <- readr::read_tsv(phospho_file, show_col_types = FALSE)
    }
    else if(grepl("\\.xlsx$", phospho_file)){
      phospho <- openxlsx::read.xlsx(phospho_file, check.names = FALSE, sep.names = " ")
    }
    else if(grepl("\\.csv$", phospho_file)){
      phospho <- readr::read_csv(phospho_file, check.names = FALSE, show_col_types = FALSE)
    }
  }
  else{
    stop("Your input can only be a file name !")
  }

  # data formatting
  abundance_columns <- grep("^Abundance: |^Abundances: ", colnames(phospho), value = TRUE)
  if(!length(abundance_columns)){
    abundance_columns <- grep("^Abundances \\(Grouped\\):", colnames(phospho), value = TRUE)
    if(!length(abundance_columns)){
      message(paste("Error: The file",  phospho_file,
                    "doesn't contain any abundance information ! It should contain columns starting with 'Abundance:' or 'Abundances (Grouped):'"))
      return()
    }
  }
  info_col <- c("Master Protein Accessions", "Master Protein Descriptions",
                "Positions in Master Proteins",
                "Annotated Sequence", "Modifications")
  miss_col <- info_col[!(info_col %in% colnames(phospho))]
  if(length(miss_col)){
    if("Master Protein Descriptions" %in% miss_col){
      message(paste("Warning: 'Master Protein Descriptions' wasn't found in", phospho_file))
      info_col <- info_col[-2]
      miss_col <- info_col[!(info_col %in% colnames(phospho))]
    }
    if(length(miss_col)){
      message(paste("Error:", paste(miss_col, collapse = ", "),
                    ifelse(length(miss_col) > 1, "are", "is"),
                    "missing in", phospho_file, "!"))
      return()
    }
  }

  # adding Modification all sit prob if there
  if("Modifications (all possible sites)" %in% colnames(phospho))
    info_col <- c(info_col, "Modifications (all possible sites)")
  else
    message("Warning: 'Modifications (all possible sites)' column was missing, you might miss some PTM probabilities!")

  # extracting peptide abundance count and PSM count
  countNum <- grep(pattern = "^Abundances Count: [A-z0-9,. -]+: 126[A-z0-9,. -]+$",
                   colnames(phospho), value = TRUE)
  psm <- grep(pattern = " PSM", colnames(phospho), value = TRUE)
  if(length(countNum) == 1){
    if(length(psm) == 1){
      phospho <- phospho[,c(info_col, abundance_columns, countNum, psm)]
      colnames(phospho)[ncol(phospho)-1] <- "countNum"
      colnames(phospho)[ncol(phospho)] <- "sumPSMs"
    }
    else{
      phospho <- phospho[,c(info_col, abundance_columns, countNum)]
      colnames(phospho)[ncol(phospho)] <- "countNum"
      message("Warning: There is no PSM count information in the original data!")
      phospho$sumPSMs <- NA
    }
  }
  else{
    if(length(psm) == 1){
      phospho <- phospho[,c(info_col, abundance_columns, psm)]
      colnames(phospho)[ncol(phospho)] <- "sumPSMs"
      message("Warning: There is no peptide abundance count information in the original data!\nNo filter will be applied.")
      phospho$countNum <- NA
      phospho <- phospho[,c(info_col, abundance_columns, "countNum", "sumPSMs")]
    }
    else{
      phospho <- phospho[,c(info_col, abundance_columns)]
      message("Warning: There are neither peptide abundance count and PSM count information in the original data!")
      phospho$countNum <- NA
      phospho$sumPSMs <- NA
    }
  }

  # removing contaminants
  if(nchar(prefixcontaminant)) {
    ncont <- nrow(phospho)
    pattern <- grep(paste0("(^|(;|; ))", prefixcontaminant),
                    phospho$`Master Protein Accessions`)
    if (length(pattern) > 0) {
      phospho <- phospho[-pattern, ]
    }
    if (nrow(phospho) == 0) {
      stop("After removing Contaminant proteins, the dataset was empty.")
    }
    pattern <- grep("(B|b)os taurus", phospho$`Master Protein Descriptions`)
    if (length(pattern) > 0) {
      phospho <- phospho[-pattern, ]
    }
    if (nrow(phospho) == 0) {
      stop("After removing Bos taurus proteins, the dataset was empty.")
    }
    message(paste0(ncont - nrow(phospho), " contaminants were removed from the data"))
  }

  # removing mix and/or empty channels if any
  colnames(phospho)[(length(info_col)+1):(ncol(phospho) - 2)] <- paste0("25C_", treatment) # adding fake temperature to use other functions
  if(length(grep("_Mix|_Empty", colnames(phospho)))){
    message(paste(paste(sub(".*_", "",
                            grep("_Mix|_Empty", colnames(phospho), value = TRUE)
                            ),
                        collapse = ", "),
                  "channel(s) has(ve) been removed."))
    phospho <- phospho[, -grep("_Mix|_Empty", colnames(phospho))]
  }

  # order data according protein id
  phospho <- phospho[order(phospho$`Master Protein Accessions`),]

  ## renaming column description if there
  descr <- grep("Master Protein Descriptions", colnames(phospho))
  if(length(descr)){
    colnames(phospho)[descr] <- "description"
  }
  else if(!("description" %in% colnames(phospho))){
    cn <- colnames(phospho)
    phospho$description <- NA
    phospho <- phospho[,c(cn[1], "description", cn[-1])]
  }

  ## filtering
  # filtering on abundance count
  if(!all(is.na(phospho$countNum))){
    n_abc <- nrow(phospho)
    phospho <- phospho[which(phospho$countNum >= countthreshold),]
    message(paste0(n_abc - nrow(phospho), " peptides didn't pass the cutoff of abundance count number ", countthreshold))
  }

  # only keep phospho peptides
  phospho <- phospho[grep("Phospho", phospho$Modifications),]

  # remove TMT modification information from modifications
  phospho$Modifications <- unlist(
    lapply(strsplit(phospho$Modifications, "\\]; "),
           function(x){
             x <- x[-grep("TMT", x)]
             if(length(grep("\\]$", x, invert = TRUE))){
               x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
             }
             x <- paste(x, collapse = "; ")
             x
           })
    )

  # remove peptides with only NAs
  n <- nrow(phospho)
  phospho <- phospho[which(apply(phospho[,6:(ncol(phospho)-2)], 1,
                                 function(x) sum(is.na(x)) < ncol(phospho) - 7)),]
  n <- n - nrow(phospho)
  message(paste(n, "peptides without quantitative information were removed"))

  message("Formatting Modifications...")
  ## Add probability to sites that eventually don't have any (PD remove this info when sveral with same prob)
  if("Modifications (all possible sites)" %in% colnames(phospho)){
    # remove TMT modification information from modifications
    phospho$`Modifications (all possible sites)` <- unlist(
      lapply(strsplit(phospho$`Modifications (all possible sites)`, "\\]; "),
             function(x){
               x <- x[-grep("TMT", x)]
               if(length(grep("\\]$", x, invert = TRUE))){
                 x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
               }
               x <- paste(x, collapse = "; ")
               x
             })
    )

    ambiguous_prob <- grep("(S|T|Y)(?!\\d{1,})", phospho$Modifications, perl = TRUE)
    if(length(ambiguous_prob)){
      for(ambp in ambiguous_prob){
        current_modif <- phospho$Modifications[ambp]
        n_phospho <- as.numeric(sub(".* ", "", sub("xPhospho.*", "", current_modif)))

        phospho_allprob <- gsub(".* \\[|\\]", "",
                                grep("Phospho", strsplit(phospho$`Modifications (all possible sites)`[ambp],
                                                         "\\]; ")[[1]],
                                     value = TRUE)
                                )
        phospho_allprob <- strsplit(phospho_allprob, "; ")[[1]]
        phospho_allprob <- sapply(phospho_allprob, function(p) as.numeric(gsub(".*\\(|\\)", "", p)))

        phospho_top <- sort(unique(phospho_allprob), decreasing = TRUE)[1:n_phospho]
        phospho_top <- phospho_top[which(phospho_top >= 10)]
        phospho_top <- phospho_allprob[which(!is.na(match(phospho_allprob, phospho_top)))]
        phospho_top <- paste(names(phospho_top), collapse = "; ")

        update_modif <- sub("(?<=Phospho \\[).*(?=\\])", phospho_top, current_modif, perl = TRUE)
        phospho$Modifications[ambp] <- update_modif
      }
    }
    phospho$`Modifications (all possible sites)` <- NULL
  }

  ## Add correct site in Modifications
  phospho$Modifications <- apply(phospho[,c("Positions in Master Proteins", "Modifications")],
        1, function(x){
          start_pos <- as.numeric(gsub(".* \\[|-.*", "",
                                       sub(";.*", "", x[["Positions in Master Proteins"]]))
                                  )

          modif_pos <- strsplit(x[["Modifications"]], "(?<=(\\[| )[A-Z])|(?=(\\(|\\]))", perl = TRUE)[[1]]
          modif_pos[grep("^\\d{1,3}$", modif_pos)] <- as.character(as.numeric(modif_pos[grep("^\\d{1,3}$", modif_pos)]) +
                                                                     start_pos - 1);
          paste(modif_pos, collapse = "")
        })

  message("Saving files")
  colnames(phospho)[1:5] <- gsub(" ", ".", colnames(phospho)[1:5])
  readr::write_tsv(phospho,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "peptides_", dataset_name, ".txt")
                   )
  return(phospho)
}
