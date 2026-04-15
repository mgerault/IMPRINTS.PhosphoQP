#' imprints_QPpep_rawread
#'
#' Function to read the quantitative proteomics data on the peptide level from Proteome Discoverer
#'
#' @param QP_file The path to the PeptidesGroups file from PD.
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

imprints_QPpep_rawread <- function(QP_file, treatment, prefixcontaminant = "",
                                   countthreshold = 1, dataset_name = "QPpeptide"){
  if(countthreshold < 0){
    countthreshold <- 0
    messsage("Parameter countthreshold can only be positive; value set to 0, hence no filtering will be applied.")
  }

  # data reading
  if(inherits(QP_file, "character")){
    stopifnot(file.exists(QP_file))

    if(grepl("\\.txt$", QP_file)){
      QP <- readr::read_tsv(QP_file, show_col_types = FALSE)
    }
    else if(grepl("\\.xlsx$", QP_file)){
      QP <- openxlsx::read.xlsx(QP_file, check.names = FALSE, sep.names = " ")
    }
    else if(grepl("\\.csv$", QP_file)){
      QP <- readr::read_csv(QP_file, check.names = FALSE, show_col_types = FALSE)
    }
  }
  else{
    stop("Your input can only be a file name !")
  }

  # data formatting
  abundance_columns <- grep("^Abundance: |^Abundances: ", colnames(QP), value = TRUE)
  if(!length(abundance_columns)){
    abundance_columns <- grep("^Abundances \\(Grouped\\):", colnames(QP), value = TRUE)
    if(!length(abundance_columns)){
      message(paste("Error: The file",  QP_file,
                    "doesn't contain any abundance information ! It should contain columns starting with 'Abundance:' or 'Abundances (Grouped):'"))
      return()
    }
  }
  info_col <- c("Master Protein Accessions", "Master Protein Descriptions",
                "Positions in Master Proteins",
                "Annotated Sequence", "Modifications")
  miss_col <- info_col[!(info_col %in% colnames(QP))]
  if(length(miss_col)){
    if("Master Protein Descriptions" %in% miss_col){
      message(paste("Warning: 'Master Protein Descriptions' wasn't found in", QP_file))
      info_col <- info_col[-2]
      miss_col <- info_col[!(info_col %in% colnames(QP))]
    }
    if(length(miss_col)){
      message(paste("Error:", paste(miss_col, collapse = ", "),
                    ifelse(length(miss_col) > 1, "are", "is"),
                    "missing in", QP_file, "!"))
      return()
    }
  }

  # extracting peptide abundance count and PSM count
  countNum <- grep(pattern = "^Abundances Count: [A-z0-9,. -]+: 126[A-z0-9,. -]+$",
                   colnames(QP), value = TRUE)
  psm <- grep(pattern = " PSM", colnames(QP), value = TRUE)
  if(length(countNum) == 1){
    if(length(psm) == 1){
      QP <- QP[,c(info_col, abundance_columns, countNum, psm)]
      colnames(QP)[ncol(QP)-1] <- "countNum"
      colnames(QP)[ncol(QP)] <- "sumPSMs"
    }
    else{
      QP <- QP[,c(info_col, abundance_columns, countNum)]
      colnames(QP)[ncol(QP)] <- "countNum"
      message("Warning: There is no PSM count information in the original data!")
      QP$sumPSMs <- NA
    }
  }
  else{
    if(length(psm) == 1){
      QP <- QP[,c(info_col, abundance_columns, psm)]
      colnames(QP)[ncol(QP)] <- "sumPSMs"
      message("Warning: There is no peptide abundance count information in the original data!\nNo filter will be applied.")
      QP$countNum <- NA
      QP <- QP[,c(info_col, abundance_columns, "countNum", "sumPSMs")]
    }
    else{
      QP <- QP[,c(info_col, abundance_columns)]
      message("Warning: There are neither peptide abundance count and PSM count information in the original data!")
      QP$countNum <- NA
      QP$sumPSMs <- NA
    }
  }

  # removing contaminants
  if (nchar(prefixcontaminant)) {
    ncont <- nrow(QP)
    pattern <- grep(paste0("(^|(;|; ))", prefixcontaminant),
                    QP$`Master Protein Accessions`)
    if (length(pattern) > 0) {
      QP <- QP[-pattern, ]
    }
    if (nrow(QP) == 0) {
      stop("After removing Contaminant proteins, the dataset was empty.")
    }
    pattern <- grep("(B|b)os taurus", QP$`Master Protein Descriptions`)
    if (length(pattern) > 0) {
      QP <- QP[-pattern, ]
    }
    if (nrow(QP) == 0) {
      stop("After removing Bos taurus proteins, the dataset was empty.")
    }
    message(paste0(ncont - nrow(QP), " contaminants were removed from the data"))
  }

  # removing mix and/or empty channels if any
  colnames(QP)[(length(info_col)+1):(ncol(QP) - 2)] <- paste0("26C_", treatment) # adding fake temperature to use other functions
  if(length(grep("_Mix|_Empty", colnames(QP)))){
    message(paste(paste(sub(".*_", "",
                            grep("_Mix|_Empty", colnames(QP), value = TRUE)
    ),
    collapse = ", "),
    "channel(s) has(ve) been removed."))
    QP <- QP[, -grep("_Mix|_Empty", colnames(QP))]
  }

  # order data according protein id
  QP <- QP[order(QP$`Master Protein Accessions`),]

  ## renaming column description if there
  descr <- grep("Master Protein Descriptions", colnames(QP))
  if(length(descr)){
    colnames(QP)[descr] <- "description"
  }
  else if(!("description" %in% colnames(QP))){
    cn <- colnames(QP)
    QP$description <- NA
    QP <- QP[,c(cn[1], "description", cn[-1])]
  }

  ## filtering
  # filtering on abundance count
  if(!all(is.na(QP$countNum))){
    n_abc <- nrow(QP)
    QP <- QP[which(QP$countNum >= countthreshold),]
    message(paste0(n_abc - nrow(QP), " peptides didn't pass the cutoff of abundance count number ", countthreshold))
  }

  # remove TMT modification information from modifications
  QP$Modifications <- unlist(
    lapply(strsplit(QP$Modifications, "\\]; "),
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
  n <- nrow(QP)
  QP <- QP[which(apply(QP[,6:(ncol(QP)-2)], 1,
                                 function(x) sum(is.na(x)) < ncol(QP) - 7)),]
  n <- n - nrow(QP)
  message(paste(n, "peptides without quantitative information were removed"))

  message("Add correct site in Modifications")
  QP$Modifications <- apply(QP[,c("Positions in Master Proteins", "Modifications")],
                            1, function(x){
                              if(nchar(x[["Modifications"]])){
                                start_pos <- as.numeric(gsub(".* \\[|-.*", "",
                                                             sub(";.*", "", x[["Positions in Master Proteins"]]))
                                )
                                modif_pos <- strsplit(x[["Modifications"]], "(?<=(\\[| )[A-Z])|(?=(\\(|\\]))", perl = TRUE)[[1]]
                                modif_pos[grep("^\\d{1,3}$", modif_pos)] <- as.character(as.numeric(modif_pos[grep("^\\d{1,3}$", modif_pos)]) +
                                                                                           start_pos - 1);
                                paste(modif_pos, collapse = "")
                              }
                              else
                                ""
                              })

  message("Saving files")
  colnames(QP)[1:5] <- gsub(" ", ".", colnames(QP)[1:5])
  readr::write_tsv(QP,
                   file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                                 "peptides_", dataset_name, ".txt")
  )
  return(QP)
}
