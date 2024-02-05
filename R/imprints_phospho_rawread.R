#' imprints_phospho_rawread
#' 
#' Function to read the phosphoproteomics data from Proteome Discoverer
#' 
#' @param data Your input data. It can the path to the output file from PD or a 
#'    data.frame corresponding to that output.
#' @param treatment The treatment corresponding to each fo your TMT channel. It should be 
#'   the same length as your number of TMT channel and should follow the format 'Rep_treatment'
#'   like 'B1_Vehicle' for example.
#' @param contanimant Character to tell which suffix indicating a protein is a contaminant.
#' 
#' @return Your data filtered and in the right format to use the IMPRINTS.CETSA package.
#'
#' @export

imprints_phospho_rawread <- function(data, treatment, contaminant = "Cont"){
  stopifnot(!any(duplicated(treatment)))
  
  if(inherits(data, "character")){
    stopifnot(file.exists(data))
    
    if(grepl("\\.txt$", data)){
      data <- readr::read_tsv(data)
    }
    else if(grepl("\\.xlsx$", data)){
      data <- openxlsx::read.xlsx(data)
    }
    else if(grepl("\\.csv$", data)){
      data <- readr::read_csv(data)
    }
  }
  else if(!inherits(data, "data.frame")){
    stop("Your input can only be a file or a dataframe !")
  }
  
  cn_needed <- c("^Sequence$", "^Modifications$", "^Modifications in Master Proteins$",
                 "^Positions in Master Proteins$", "^Master Protein Descriptions$",
                 "PSMs$")
  cn <- colnames(data)
  check_cn_needed <- sapply(cn_needed,
                            function(x) length(grep(x, cn)) == 1)
  
  if(!(all(check_cn_needed))){
    cn_missing <- cn_needed[!check_cn_needed]
    stop(paste("The", ifelse(length(cn_missing) > 1, "columns", "column"),
               paste(cn_missing, collapse = ", "),
               ifelse(length(cn_missing) > 1, "are", "is"),
               "missing in your data !"))
  }
  else if(length(grep("Abundance:|Abundance F\\d{1,}", cn)) != length(treatment)){
    stop("Your number of raw abundance columns doesn't match your number of treatments !")
  }
  
  # remove potential contaminants
  contaminants_torm <- grep("Bos taurus", data[["Master Protein Descriptions"]])
  contaminants_torm <- c(contaminants_torm, grep(paste0("^", contaminant), data[["Positions in Master Proteins"]]))
  if("Contaminant" %in% cn){
    contaminants_torm <- c(contaminants_torm, which(data$Contaminant))
  }
  if(length(contaminants_torm)){
    contaminants_torm <- unique(contaminants_torm)
    message(paste(length(contaminants_torm), "contaminants have been removed"))
    data <- data[-contaminants_torm,]
  }
  
  # select columns
  data <- data[,grep("^Sequence$|^Modifications$|^(Modifications|Positions) in Master Proteins$|^Master Protein Descriptions$|PSMs$|(Abundance:|Abundance F\\d{1,})",
                     colnames(data)
  )
  ]
  
  ## renaming columns
  colnames(data)[grep("PSMs$", colnames(data))] <- "sumPSMs"
  
  colnames(data)[grep("^Abundance:|Abundance F\\d{1,}", colnames(data))] <- treatment
  
  colnames(data) <- gsub(" ", ".", colnames(data))
  
  data$Mix <- NULL
  
  ## filtering
  # only keep pep with quantitative info
  data <- data[apply(data[7:ncol(data)], 1, function(x) sum(is.na(x)) < length(treatment) - 1),] 
  
  # only keep phospho peptides
  data <- data[grep("Phospho", data$Modifications.in.Master.Proteins),]
  
  
  # remove TMT modification information from modifications
  data$Modifications <- unlist(
    lapply(strsplit(data$Modifications, "\\]; "), 
           function(x){
             x <- x[-grep("TMTpro", x)]
             if(length(grep("\\]$", x, invert = TRUE))){
               x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
             }
             x <- paste(x, collapse = "; ")
             x
           })
  )
  
  data$Modifications.in.Master.Proteins <- NULL
  
  
  # put in right format in order to use existing package IMPRINTS.CETSA
  colnames(data)[c(1,2,4,5)] <- c("id", "sumUniPeps", "description", "countNum")
  data <- data[,c(1,4, 6:ncol(data), 2,3,5)]
  colnames(data)[3:(ncol(data)-3)] <- paste0("25C_", colnames(data)[3:(ncol(data)-3)])
  
  # make sure every id is different
  data$id <- paste0(data$id, "_", 1:nrow(data))
  # need to order according id for imprints_normalization function (expression used)
  data <- data[order(data$id),] 
  
  return(data)
}
