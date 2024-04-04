#' imprints_phoQP_join
#'
#' Function to join phosphoproteomics data, quantitative data on protein and peptide level
#'
#' @param phospho Phosphoproteomics data processed  after \code{imprints_phospho_rawread}
#' @param QPpep Quantitative proteomics data on the peptide level processed  after \code{imprints_QPpep_rawread}
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
      ### could put in same format as joined --> then use imprints_barplotting phospho !
      data <- data_list[[data_type]]
      if(data_type == "phospho"){
        data$phospho.site <- unlist(lapply(strsplit(data$sumUniPeps, "\\]; "),
                                           function(x){
                                             x <- x[grep("Phospho", x)]
                                             if(length(grep("\\]$", x, invert = TRUE))){
                                               x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                             }
                                             x <- paste(x, collapse = "; ");
                                             x
                                           }
        )
        )
        data$sumUniPeps <- unlist(lapply(strsplit(data$sumUniPeps, "\\]; "),
                                         function(x){
                                           x <- x[-grep("Phospho", x)]
                                           if(length(grep("\\]$", x, invert = TRUE))){
                                             x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                           }
                                           x <- paste(x, collapse = "; ");
                                           x
                                         }
        )
        )
        colnames(data) <- sub("sumPSMs", "sumPSMs.phospho", colnames(data))
      }
      else if(data_type == "QPpep"){
        colnames(data) <- sub("sumPSMs", "sumPSMs.QPpep", colnames(data))
      }

      if(data_type != "QP"){
        data$Position.in.Master.Protein <- gsub("[A-Z].{5,9} ", "", data$description)
        data$Master.Protein <- gsub(" \\[.{,9}\\]", "", data$description)
        data$Master.Protein <- sub(";$", "", data$Master.Protein)
        data$description <- NULL

        colnames(data) <- sub("^id$", "Sequence", colnames(data))
        colnames(data) <- sub("^countNum$", "description", colnames(data))
        colnames(data) <- sub("^sumUniPeps$", "Modifications", colnames(data))

        data <- data[,c("Sequence", grep("phospho.site", colnames(data), value = TRUE),
                        "Modifications", "Master.Protein",
                        "Position.in.Master.Protein", "description",
                        grep("^\\d{2}C_", colnames(data), value = TRUE),
                        grep("^sumPSMs", colnames(data), value = TRUE)
        )
        ]
      }

      return(data)
    }
    else if(all(c("phospho", "QPpep") %in% data_type)){
      ## merge phospho and QP peptide

      # prepare phospho
      phospho <- data_list$phospho
      phospho$id <- sub("_.*", "", phospho$id)
      phospho$phospho.site <- unlist(lapply(strsplit(phospho$sumUniPeps, "\\]; "),
                                            function(x){
                                              x <- x[grep("Phospho", x)]
                                              if(length(grep("\\]$", x, invert = TRUE))){
                                                x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                              }
                                              x <- paste(x, collapse = "; ");
                                              x
                                            }
      )
      )
      phospho$sumUniPeps <- unlist(lapply(strsplit(phospho$sumUniPeps, "\\]; "),
                                          function(x){
                                            x <- x[-grep("Phospho", x)]
                                            if(length(grep("\\]$", x, invert = TRUE))){
                                              x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                            }
                                            x <- paste(x, collapse = "; ");
                                            x
                                          }
      )
      )
      colnames(phospho) <- sub("sumPSMs", "sumPSMs.phospho", colnames(phospho))

      # prepare QPpep
      QPpep <- data_list$QPpep
      QPpep$id <- sub("_.*", "", QPpep$id)
      colnames(QPpep) <- sub("sumPSMs", "sumPSMs.QPpep", colnames(QPpep))

      # join phospho and QPpep
      if(method == "full"){
        data <- dplyr::full_join(phospho, QPpep,
                                 by = c("id", "description", "sumUniPeps", "countNum"))

        nb_common <- dplyr::inner_join(phospho, QPpep,
                                       by = c("id", "description", "sumUniPeps", "countNum"))
        nb_common <- nrow(nb_common)
      }
      else if(method == "inner"){
        data <- dplyr::inner_join(phospho, QPpep,
                                  by = c("id", "description", "sumUniPeps", "countNum"))
        nb_common <- nrow(data)
      }

      message(paste0(nb_common, " peptides are in common with the phospho data and the QP data. \nHence, ", round(100*nb_common/nrow(phospho), 2), "% of the phospho peptides don't have a corresponding peptide in the QP data."))
    }
    else{
      peptide_type <- c("phospho", "QPpep")[c("phospho", "QPpep") %in% data_type]
      data <- data_list[[peptide_type]]
      data$id <- sub("_.*", "", data$id)

      if(peptide_type == "phospho"){
        data$phospho.site <- unlist(lapply(strsplit(data$sumUniPeps, "\\]; "),
                                           function(x){
                                             x <- x[grep("Phospho", x)]
                                             if(length(grep("\\]$", x, invert = TRUE))){
                                               x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                             }
                                             x <- paste(x, collapse = "; ");
                                             x
                                           }
        )
        )
        data$sumUniPeps <- unlist(lapply(strsplit(data$sumUniPeps, "\\]; "),
                                         function(x){
                                           x <- x[-grep("Phospho", x)]
                                           if(length(grep("\\]$", x, invert = TRUE))){
                                             x[grep("\\]$", x, invert = TRUE)] <- paste0(x[grep("\\]$", x, invert = TRUE)], "]")
                                           }
                                           x <- paste(x, collapse = "; ");
                                           x
                                         }
        )
        )
        colnames(data) <- sub("sumPSMs", "sumPSMs.phospho", colnames(data))
      }
      else if(peptide_type == "QPpep"){
        colnames(data) <- sub("sumPSMs", "sumPSMs.QPpep", colnames(data))
      }
    }

    data$Position.in.Master.Protein <- gsub("[A-Z].{5,9} ", "", data$description)
    data$Master.Protein <- gsub(" \\[.{,9}\\]", "", data$description)
    data$Master.Protein <- sub(";$", "", data$Master.Protein)
    data$description <- NULL

    colnames(data) <- sub("^id$", "Sequence", colnames(data))
    colnames(data) <- sub("^countNum$", "description", colnames(data))
    colnames(data) <- sub("^sumUniPeps$", "Modifications", colnames(data))

    data <- data[,c("Sequence", grep("phospho.site", colnames(data), value = TRUE),
                    "Modifications", "Master.Protein",
                    "Position.in.Master.Protein", "description",
                    grep("^\\d{2}C_", colnames(data), value = TRUE),
                    grep("^sumPSMs", colnames(data), value = TRUE)
    )
    ]

    if("QP" %in% data_type){
      # merge QP with peptide level data

      QP <- data_list$QP
      colnames(QP)[grep("sumUniPeps|sumPSMs|countNum", colnames(QP))] <- paste0(grep("sumUniPeps|sumPSMs|countNum", colnames(QP), value = TRUE),
                                                                                ".QP")
      colnames(QP) <- sub("^id$", "Master.Protein", colnames(QP))

      if(method == "full"){
        data <- dplyr::full_join(data, QP, by = c("Master.Protein", "description"))
      }
      else if(method == "inner"){
        data <- dplyr::inner_join(data, QP, by = c("Master.Protein", "description"))
      }

      data <- data[,c("Sequence", grep("phospho.site", colnames(data), value = TRUE),
                      "Modifications", "Master.Protein",
                      "Position.in.Master.Protein", "description",
                      grep("^\\d{2}C_", colnames(data), value = TRUE),
                      "sumUniPeps.QP",
                      grep("^sumPSMs", colnames(data), value = TRUE),
                      "countNum.QP"
                      )

      ]
    }
  }

  data <- data[order(unlist(sapply(data$description, IMPRINTS.CETSA:::getGeneName))),]

  if(!is.null(data$Modifications)){
    modif_pos <- mapply(function(org, pos){
      if(!is.na(pos)){
        if(nchar(pos)){
          org <- as.numeric(gsub("\\[|-.*", "", org)) - 1

          pos_char <- paste(regmatches(pos,  gregexpr("\\[.*?\\]", pos))[[1]], collapse = " ")
          pos_char <-  regmatches(pos_char,  gregexpr("\\D{1}\\d{1,2}?(\\]|; )", pos_char))[[1]]
          if(length(pos_char)){
            pos_num <- as.numeric(gsub("\\D", "", pos_char))
            pos_num_new <- pos_num + org

            pos_modif <- mapply(function(old, new, oldchr){
              sub(old, new, oldchr)
            },
            as.character(pos_num), as.character(pos_num_new), pos_char,
            SIMPLIFY = FALSE)

            pos_char <- gsub("\\]", "\\\\]", pos_char)
            for(i in 1:length(pos_char)){
              pos <- sub(pos_char[i], pos_modif[i], pos)
            }
          }
        }
      };
      pos
    },
    data$Position.in.Master.Protein, data$Modifications,
    SIMPLIFY = FALSE)
    data$Modifications <- unlist(modif_pos)
  }

  if(!is.null(data$phospho.site)){
    pho_pos <- mapply(function(org, pos){
      if(!is.na(pos)){
        if(nchar(pos)){
          org <- as.numeric(gsub("\\[|-.*", "", org)) - 1

          pos_char <- regmatches(pos,  gregexpr("\\[.*?\\]", pos))[[1]]
          pos_char <-  regmatches(pos_char,  gregexpr("\\D{1}\\d{1,2}?\\(", pos_char))[[1]]
          if(length(pos_char)){
            pos_num <- as.numeric(gsub("\\D", "", pos_char))
            pos_num_new <- pos_num + org

            pos_modif <- mapply(function(old, new, oldchr){
              sub(old, new, oldchr)
            },
            as.character(pos_num), as.character(pos_num_new), pos_char,
            SIMPLIFY = FALSE)

            pos_char <- gsub("\\(", "\\\\(", pos_char)
            for(i in 1:length(pos_char)){
              pos <- sub(pos_char[i], pos_modif[i], pos)
            }
          }
        }
      };
      pos
    },
    data$Position.in.Master.Protein, data$phospho.site,
    SIMPLIFY = FALSE)
    data$phospho.site <- unlist(pho_pos)
  }

  return(data)
}

