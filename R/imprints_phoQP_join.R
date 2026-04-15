#' imprints_phoQP_join
#'
#' Function to join phosphoproteomics data, quantitative data on protein and peptide level
#'
#' @param phospho Phosphoproteomics data processed after \code{imprints_phospho_rawread}
#' @param QPpep Quantitative proteomics data on the peptide level processed after \code{imprints_QPpep_rawread}
#' @param QP The quantitative proteomics data on the protein level; it can also contains the CETSA data
#' @param method Either full or inner join
#'
#' @return The joined dataset
#'
#' @export

imprints_phoQP_join <- function(phospho = NULL, QPpep = NULL, QP = NULL,
                                method = c("full", "inner")){
  method <- match.arg(method)

  data_list <- list("phospho" =  phospho,
                    "QPpep" = QPpep,
                    "QP" = QP)
  check_data <- unlist(lapply(data_list, is.null))

  if(all(check_data)){ # no data input
    return(NULL)
  }
  else{
    data_type <- names(data_list)[which(!check_data)]
    if(length(data_type) == 1){ # no merge
      data <- data_list[[data_type]]
      if(data_type == "phospho"){ # separating Phospho from the other PTMs
        data$Phosphorylation <- sapply(strsplit(data$Modifications, "\\]; "),
                                       function(x){
                                         x <- x[grep("Phospho", x)]
                                         if(length(grep("\\]$", x, invert = TRUE))){
                                           x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                         }
                                         x <- paste(x, collapse = "; ");
                                         x
                                       })

        data$Modifications <- sapply(strsplit(data$Modifications, "\\]; "),
                                     function(x){
                                       x <- x[-grep("Phospho", x)]
                                       if(length(grep("\\]$", x, invert = TRUE))){
                                         x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                       }
                                       x <- paste(x, collapse = "; ");
                                       x
                                     })
      }

      return(data)
    }
    else if(all(c("phospho", "QPpep") %in% data_type)){
      ## merge phospho and QP peptide
      # prepare phospho
      phospho <- data_list$phospho
      colnames(phospho)[c(ncol(phospho) - 1, ncol(phospho))] <- paste0(colnames(phospho)[c(ncol(phospho) - 1, ncol(phospho))],
                                                                       ".phospho")

      phospho$Phosphorylation <- sapply(strsplit(phospho$Modifications, "\\]; "),
                                        function(x){
                                          x <- x[grep("Phospho", x)]
                                          if(length(grep("\\]$", x, invert = TRUE))){
                                            x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                          }
                                          x <- paste(x, collapse = "; ");
                                          x
                                        })

      phospho$Modifications <- sapply(strsplit(phospho$Modifications, "\\]; "),
                                      function(x){
                                        x <- x[-grep("Phospho", x)]
                                        if(length(grep("\\]$", x, invert = TRUE))){
                                          x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                        }
                                        x <- paste(x, collapse = "; ");
                                        x
                                      })

      # prepare QPpep
      QPpep <- data_list$QPpep
      colnames(QPpep)[c(ncol(QPpep) - 1, ncol(QPpep))] <- paste0(colnames(QPpep)[c(ncol(QPpep) - 1, ncol(QPpep))],
                                                                       ".QPpep")

      # join phospho and QPpep
      if(method == "full"){
        data <- dplyr::full_join(phospho, QPpep,
                                 by = c("Master.Protein.Accessions", "description",
                                        "Positions.in.Master.Proteins", "Annotated.Sequence", "Modifications"))

        nb_common <- dplyr::inner_join(phospho, QPpep,
                                       by = c("Master.Protein.Accessions", "description",
                                              "Positions.in.Master.Proteins", "Annotated.Sequence", "Modifications"))
        nb_common <- nrow(nb_common)
      }
      else if(method == "inner"){
        data <- dplyr::inner_join(phospho, QPpep,
                                  by = c("Master.Protein.Accessions", "description",
                                         "Positions.in.Master.Proteins", "Annotated.Sequence", "Modifications"))
        nb_common <- nrow(data)
      }

      message(paste0(nb_common, " peptides are in common with the phospho data and the QP data. \nHence, ",
                     round(100*nb_common/nrow(phospho), 2), "% of the phospho peptides don't have a corresponding peptide in the QP data."))
    }
    else{
      peptide_type <- c("phospho", "QPpep")[c("phospho", "QPpep") %in% data_type]
      data <- data_list[[peptide_type]]

      if(peptide_type == "phospho"){
        colnames(data)[c(ncol(data) - 1, ncol(data))] <- paste0(colnames(data)[c(ncol(data) - 1, ncol(data))],
                                                                ".phospho")

        data$Phosphorylation <- sapply(strsplit(data$Modifications, "\\]; "),
                                    function(x){
                                      x <- x[grep("Phospho", x)]
                                      if(length(grep("\\]$", x, invert = TRUE))){
                                        x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                      }
                                      x <- paste(x, collapse = "; ");
                                      x
                                    })

        data$Modifications <- sapply(strsplit(data$Modifications, "\\]; "),
                                     function(x){
                                       x <- x[-grep("Phospho", x)]
                                       if(length(grep("\\]$", x, invert = TRUE))){
                                         x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                       }
                                       x <- paste(x, collapse = "; ");
                                       x
                                      })
      }
      else if(peptide_type == "QPpep"){
        colnames(data)[c(ncol(data) - 1, ncol(data))] <- paste0(colnames(data)[c(ncol(data) - 1, ncol(data))],
                                                                ".QPpep")
      }
    }

    if("QP" %in% data_type){
      # merge QP with peptide level data
      QP <- data_list$QP
      colnames(QP)[grep("sumUniPeps|sumPSMs|countNum", colnames(QP))] <- paste0(grep("sumUniPeps|sumPSMs|countNum", colnames(QP), value = TRUE),
                                                                                ".QP")
      colnames(QP) <- sub("^id$", "Master.Protein.Accessions", colnames(QP))

      if(method == "full"){
        data <- dplyr::full_join(data, QP, by = c("Master.Protein.Accessions", "description"))
      }
      else if(method == "inner"){
        data <- dplyr::inner_join(data, QP, by = c("Master.Protein.Accessions", "description"))
      }
    }
  }

  data <- data[order(data$Master.Protein.Accessions),
               c("Master.Protein.Accessions", "description",
                 "Positions.in.Master.Proteins",
                 "Annotated.Sequence", "Modifications",
                 grep("Phosphorylation", colnames(data), value = TRUE),
                 grep("^\\d{2}C_", colnames(data), value = TRUE),
                 grep("sumUniPeps", colnames(data), value = TRUE),
                 grep("sumPSMs", colnames(data), value = TRUE),
                 grep("countNum", colnames(data), value = TRUE))]

  return(data)
}

