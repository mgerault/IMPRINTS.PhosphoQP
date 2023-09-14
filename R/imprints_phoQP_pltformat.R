#' imprints_phoQP_pltformat
#'
#' Function to put the processed data after imprints_QPpep_rawread or imprints_phospho_rawread
#' in the right format to use \code{\link[IMPRINTS.CETSA.app:imprints_barplotting_app]{imprints_barplotting_app}}
#'
#' @param data Processed data after \code{imprints_QPpep_rawread} or \code{imprints_phospho_rawread}
#' @param hits A vector of ids that you would like to keep 
#' @param reftreatment A character corresponding to the reference treatment that you 
#'   wish to remove (can be NULL)
#'
#' @return A filter data frame ready to use in \code{\link[IMPRINTS.CETSA.app:imprints_barplotting_app]{imprints_barplotting_app}}
#' 
#' @export


imprints_phoQP_pltformat <- function(data, hits = NULL, reftreatment = NULL){
  if(length(hits)){
    data <- data[which(!is.na(match(data$id, unique(hits)))),]
  }
  
  # order according genes
  data <- data[order(unname(unlist(sapply(data$countNum, IMPRINTS.CETSA:::getGeneName)))),]
  
  data$id <- sub("_.*", "", data$id)
  data$id <- paste(data$id, data$sumUniPeps, data$description, 
                   sep = "\n"
  )
  data$id <- sub("\n\n", "\n", data$id)
  colnames(data)[c(2,ncol(data))] <- c("countNum", "description")
  
  if(length(grep(paste0("_", reftreatment, "$"), colnames(data)))){
    data <- data[,-grep(paste0("_", reftreatment, "$"), colnames(data))]
  }
  
  return(data)
}