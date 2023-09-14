#' imprints_phoQP_barplotting
#' 
#' Function to print or save barplots from peptide data and protein data like 
#' phosphoproteomics and quantitative proteomics and CETSA.
#' 
#' @param data Dataframe obtained after \code{imprints_phoQP_join} or any output from 
#'   \code{\link[IMPRINTS.CETSA::imprints_caldiff]{imprints_caldiff}}
#' @param treatmentlevel A vector of treatment labels, such as c("DMSO","TNFa","AT26533")
#'                       the order determines the arrangement, so in this case DMSO
#'                       group would be the first group
#' @param printBothName A logical to tell if you want to print the both protein names on the plot
#' @param printGeneName A logical to tell if you want to print the gene names on the plot
#' @param pfdatabase A logical for using pdf database or not
#' @param witherrorbar A logical to print or not the error bar on the plot
#' @param layout A vector indicating the panel layout for multi-panel plots per page,
#'               default value is c(2,3) for set containing data, otherwise c(4,3), use when save_pdf = TRUE
#' @param colorpanel A vector of customizable color scheme provided by default with the function PaletteWithoutGrey
#' @param ratio Aspect ratio of the plot, default set to 0.6
#' @param ret_plot Logical to tell if you want to return the last plot
#' @param save_pdf A logical to tell if you want to save plots in a pdf file
#' @param toplabel Textual label at the top part of the page
#' @param leftlabel Textual label at the left side of the page
#' @param bottomlabel Textual label at the bottom part of the page
#' @param pdfname Textual label of the pdf file
#' @param pdfheight A number indicate the height of pdf file, default value 12
#' @param pdfwidth A number indicate the width of pdf file, default value 12
#' 
#' @return The imprints barplot
#' 
#' @export

imprints_phoQP_barplotting <- function (data, treatmentlevel = get_treat_level(data),
                                        printBothName = TRUE, printGeneName = FALSE,
                                        pfdatabase = FALSE, witherrorbar = TRUE, layout = NULL,
                                        colorpanel = PaletteWithoutGrey(treatmentlevel),
                                        ratio = 0.6, ret_plot = TRUE,
                                        save_pdf = FALSE, toplabel = "IMPRINTS-CETSA bar plotting",
                                        leftlabel = "", bottomlabel = "", pdfname = "barplot",
                                        pdfheight = 12, pdfwidth = 12){
  # barplotting function
  barplotting <- function(d1) {
    if(any(is.na(d1$mean))){
      d1$mean[which(is.na(d1$mean))] <- NA
    }
    minreading = -0.5
    maxreading = 0.5
    
    legendscale = c(min(min(d1$mean, na.rm = T) -
                          0.1, minreading), max(max(d1$mean, na.rm = T) +
                                                  0.1, maxreading))
    
    d1$QP <- FALSE
    if("26C" %in% d1$temperature){
      d1$QP[which(d1$temperature == "26C")] <- TRUE
      lvl_tokeep <- levels(d1$condition)
      lvl_tokeep <- gsub("26C", "QPpep", lvl_tokeep)
      d1$condition <- as.character(d1$condition)
      d1$condition[which(d1$temperature == "26C")] <- gsub("26C", "QPpep",
                                                           d1$condition[which(d1$temperature == "26C")]
                                                           )
      d1$condition <- factor(d1$condition, levels = lvl_tokeep)
      
      d1$axis1 <- d1$mean
      d1$axis1[grep("QPpep", d1$condition)] <- NA
      d1$axis2 <- d1$mean
      d1$axis2[-grep("QPpep", d1$condition)] <- NA
    }
    if("36C" %in% d1$temperature){
      d1$QP[which(d1$temperature == "36C")] <- TRUE
      lvl_tokeep <- levels(d1$condition)
      lvl_tokeep <- gsub("36C", "QP", lvl_tokeep)
      d1$condition <- as.character(d1$condition)
      d1$condition[which(d1$temperature == "36C")] <- gsub("36C", "QP",
                                                           d1$condition[which(d1$temperature == "36C")]
                                                           )
      d1$condition <- factor(d1$condition, levels = lvl_tokeep)
    }
    if("25C" %in% d1$temperature){
      lvl_tokeep <- levels(d1$condition)
      lvl_tokeep <- gsub("25C", "Phospho", lvl_tokeep)
      d1$condition <- as.character(d1$condition)
      d1$condition[which(d1$temperature == "25C")] <- gsub("25C", "Phospho", 
                                                           d1$condition[which(d1$temperature == "25C")]
                                                           )
      d1$condition <- factor(d1$condition, levels = lvl_tokeep)
      
      d1$axis1 <- d1$mean
      d1$axis1[grep("Phospho", d1$condition)] <- NA
      d1$axis2 <- d1$mean
      d1$axis2[-grep("Phospho", d1$condition)] <- NA
    }
    
    if(any(c("25C", "26C") %in% d1$temperature)){
      d1$level <- NA
      d1$level[grep("Phospho|QPpep", d1$condition)] <- "peptide"
      d1$level[-grep("Phospho|QPpep", d1$condition)] <- "protein"
      d1$level <- factor(d1$level, levels = c("protein", "peptide"))
      
      d1_label <- function(string){
        lb <- strsplit(as.character(d1$id[1]), "\n")[[1]]
        lb_pep <- paste(lb[c(1,2,3,5,6)], collapse = "\n")
        lb_prot <- paste(lb[c(4,5,6)], collapse = "\n")
        
        string <- gsub("^protein$", lb_prot, string)
        string <- gsub("^peptide$", lb_pep, string)
      }
      
      legendscale1 = c(min(min(d1$axis1, na.rm = T) - 0.1, minreading),
                       max(max(d1$axis1, na.rm = T) + 0.1, maxreading))
      legendscale2 = c(min(min(d1$axis2, na.rm = T) - 0.1, minreading),
                       max(max(d1$axis2, na.rm = T) + 0.1, maxreading))
      d1$ymin <- NA
      d1$ymax <- NA
      d1$ymin[d1$level == "protein"] <- legendscale1[1]
      d1$ymax[d1$level == "protein"] <- legendscale1[2]
      d1$ymin[d1$level == "peptide"] <- legendscale2[1]
      d1$ymax[d1$level == "peptide"] <- legendscale2[2]
    }
    
    q <- ggplot(d1, aes(x = condition, y = mean, fill = treatment)) +
      geom_bar(stat = "identity", aes(color = QP), size = rel(0.85)) +
      scale_fill_manual(drop = FALSE, values = colorpanel) +
      scale_color_manual(values = c("TRUE" = "#656565", "FALSE" = "#FFFFFF00")) +
      guides(color = "none") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9))
    
    
    
    if(any(c("25C", "26C") %in% d1$temperature)){
      q <- q + facet_wrap(~level, scales = "free", labeller = labeller(level = d1_label)) +
        scale_x_discrete(labels = gsub("_.{1,}", "", levels(d1$condition)),
                         breaks = levels(d1$condition))  +
        geom_blank(aes(y = ymin)) +
        geom_blank(aes(y = ymax))
    }
    else{
      q <- q + coord_cartesian(ylim = legendscale) +
        ggtitle(as.character(unique(d1$id))) +
        scale_x_discrete(labels = gsub("_.{1,}", "", levels(d1$condition)))
    }
    q <- q + labs(y = "fold change(log2)",
                  subtitle = subt[as.character(unique(d1$id)), "category"]) +
      cowplot::theme_cowplot() + 
      theme(text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "plain"),
            strip.text = element_text(size = rel(0.8)),
            legend.background = element_rect(fill = NULL),
            legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.15,"cm"),
            legend.title = element_text(face = "bold"),
            legend.text = element_text(size = rel(0.7)),
            legend.justification = "center", panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), strip.background = element_blank(),
            axis.line.x = element_line(), axis.line.y = element_line(),
            axis.text.x = element_text(angle = 45, hjust = 1,
                                       size = rel(0.7)),
            axis.title.x = element_blank(),
            aspect.ratio = 0.6)
    
    return(q)
  }

  if(save_pdf){
    dataname <- deparse(substitute(data))
  }

  nrowdata <- nrow(data)
  if (nrowdata == 0) {
    message("Make sure there are more than one experimental condition in dataset.")
    stop("Otherwise specify remsinglecondprot==FALSE !")
  }
  
  # checking columns and adding id
  if(!("id" %in% colnames(data))){
    cn_needed <- c("Sequence", "Master.Protein", "Position.in.Master.Protein")
    if(all(cn_needed %in% colnames(data))){
      if(NA %in% data$phospho.site){
        data$phospho.site[which(is.na(data$phospho.site))] <- ""
      }
      ids <- paste(paste(data$Sequence, data$Modifications),
                   data$phospho.site,
                   paste(data$Master.Protein, data$Position.in.Master.Protein),
                   data$Master.Protein,
                   sep = "\n")
      
      data$id <- ids
    }
    else{
      cn_missing <- cn_needed[!(cn_needed %in% colnames(data))]
      stop(paste0("If you don't have the column 'id' in your data, you need at least the columns
                  Sequence, Master.Protein and Position.in.Master.Protein 
                  The ", ifelse(length(cn_missing) > 1, "columns", "column"), paste(cn_missing, collapse = ", "),
                  ifelse(length(cn_missing) > 1, "are", "is"), "missing !"))
    }
  }
  
  if (printBothName & !pfdatabase) {
    data <- data %>% dplyr::rowwise() %>% 
      dplyr::mutate(description1 = getProteinName(description,pfdatabase)) %>%
      dplyr::mutate(description2 = getGeneName(description)) %>%
      dplyr::mutate(id = paste(id, description1, description2, sep = "\n"))
      data$description1 <- NULL
      data$description2 <- NULL
  }
  else if (printGeneName & !pfdatabase) {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getGeneName(description)) %>%
      dplyr::mutate(id = paste(id, description, sep = "\n"))
  }
  else {
    data <- data %>% dplyr::rowwise() %>%
      dplyr::mutate(description = getProteinName(description, pfdatabase)) %>%
      dplyr::mutate(id = paste(id, description, sep = "\n"))
  }

  if(length(grep("^category", names(data)))){
    if(length(grep("^score", names(data)))){
      subt <- data[, c(1, grep("^category", names(data)), grep("^score", names(data)))]
      subt <- as.data.frame(subt)
      colnames(subt) <- c("id", "category", "score")
      subt$category <- paste("Category :", subt$category, ", Score :", round(subt$score,4))
      subt$score <- NULL
      rownames(subt) <- subt$id
      data <- data[,-grep("^score", names(data))]
      data[grep("^category", names(data))] <- subt$category
    }
    else{
      subt <- data[, c(1, grep("^category", names(data)))]
      subt <- as.data.frame(subt)
      colnames(subt) <- c("id", "category")
      subt$category <- paste("Category :", subt$category)
      rownames(subt) <- subt$id
      data[grep("^category", names(data))] <- subt$category
    }

    ord_data <- data[NULL,]
    for(i in c("CN", "NC", "CC", "ND", "NN")){
      cat_idx <- grep(paste0("^Category : ", i), data$category)
      if(length(cat_idx) > 0){
        cat_idx <- data[cat_idx,]
        w <- cat_idx[order(cat_idx$category, decreasing = TRUE),]
        ord_data <- rbind(ord_data, w)
      }
    }

    if(nrow(ord_data) != 0)
      data <- ord_data
    }
  else if(length(grep("^score", names(data)))){
    subt <- data[, c(1, grep("^score", names(data)))]
    subt <- as.data.frame(subt)
    colnames(subt) <- c("id", "score")
    subt$category <- paste("Score :", round(subt$score,4))  # keep same name for simplicity
    subt$score <- NULL
    rownames(subt) <- subt$id
    data <- data[order(data[,grep("^score", names(data))], decreasing = TRUE),]
    data <- data[,-grep("^score", names(data))]
  }
  else{
    subt <- NULL
  }

  data$description <- NULL
  data1 <- tidyr::gather(data[, grep("^id$|^\\d{2}C_", colnames(data))],
                         condition, reading, -id)
  if(length(treatmentlevel) != length(get_treat_level(data))){
    data1 <- data1 %>% 
      filter(grepl(paste(treatmentlevel, collapse = "|"), condition))
  }
  a <- data1$condition[1]
  if (length(unlist(strsplit(a, "_"))) == 3) {
    data1 <- tidyr::separate(data1, condition, into = c("temperature",
                                                        "replicate", 
                                                        "treatment"), sep = "_")
    temperature <- sort(unique(data1$temperature))
    temp_idx <- grep("^[0-9]", temperature)
    if(length(temp_idx) != length(temperature)){
      temperature <- c(sort(temperature[-temp_idx]), sort(temperature[temp_idx]))
    }
    data1$id <- factor(data1$id, levels = unique(data1$id), ordered = TRUE) # preserve order
    cdata <- plyr::ddply(data1, c("id", "temperature", "treatment"),
                         summarise, N = length(na.omit(reading)), mean = mean(reading,na.rm = T),
                         sd = sd(reading, na.rm = T), se = sd/sqrt(N))
    cdata$id <- as.character(cdata$id)

    if (length(layout) == 0) {
      layout <- c(3, 1)
    }
  }
  else {
    stop("make sure the namings of the columns of the dasaset are correct.")
  }
  cdata <- cdata %>% dplyr::rowwise() %>% 
   dplyr::mutate(condition = paste(temperature, treatment, sep = "_"))
  
  cdata$id <- factor(cdata$id, levels = unique(cdata$id), ordered = TRUE)
  
  cdata$treatment <- factor(as.character(cdata$treatment),
                            levels = treatmentlevel)
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = apply(expand.grid(temperature, treatmentlevel),
                                           1, paste, collapse = "_"))
  # if data with different temperatures, prevent from creating non sense factors
  cdata$condition <- factor(as.character(cdata$condition),
                            levels = levels(cdata$condition)[levels(cdata$condition)
                                                             %in% as.character(cdata$condition) 
                                                             ]
                            )

  message("Generating fitted plot, pls wait.")

  plots <- plyr::dlply(cdata, plyr::.(id), .fun = barplotting)

  if(save_pdf){
    message("Start saving plot")

    params <- list(nrow = layout[1], ncol = layout[2])
    n <- with(params, nrow * ncol)
    pages <- length(plots)%/%n + as.logical(length(plots)%%n)
    groups <- split(seq_along(plots), gl(pages, n, length(plots)))
    n_p <- length(names(groups))
    pl <- lapply(names(groups), function(i) {
      message(paste("Saving page", i, "/", n_p))
      do.call(gridExtra::arrangeGrob,
              c(plots[groups[[i]]], params,
                top = toplabel, left = leftlabel,
                bottom = bottomlabel)
              )
      })

    class(pl) <- c("arrangelist", "ggplot", class(pl))
    pdfname <- paste0(pdfname, ".pdf")
    ggsave(file = paste0(format(Sys.time(), "%y%m%d_%H%M_"),
                         dataname, "_", pdfname), 
           pl, height = pdfheight,
           width = pdfwidth)

    message("IMPRINTS-CETSA bar plot file generated successfully.")
  }

  if(ret_plot){
    message("IMPRINTS-CETSA bar plot generated successfully.")
    return(plots)
  }
  else{
    g <- ggplot(data.frame(x = c(0,1), y = c(0,1)), aes(x,y, label = "s")) +
      geom_text(x=0.5, y=0.5, label = "All the barplots has been saved succesfully !
                                       \nGo check your files", size = 6) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank()
            )

    return(g)
  }
}



### PaletteWithoutGrey function ###
#generates a color list depending on the number of element of a character vector
PaletteWithoutGrey <- function(treatment){
  
  n = length(unique(treatment))
  x <- grDevices::colors(distinct = TRUE)                           #all the color from R
  mycol <- x[-grep("gr(e|a)y", x)]   #keep only colors that are not grey
  
  listcolor <- c()
  for (i in 0:(n-1)){
    listcolor <- append(listcolor, mycol[((i*20 + 9) %% length(mycol)) + 1])      #save a color from the list (the number 20 and 9 were chosen in order to have distincts colors, this is empirical, can be changed)
  }
  
  return(listcolor)
}

getGeneName <- function (x){
  gene = strsplit(strsplit(x, "GN=")[[1]][2], " ")[[1]][1]
  if (length(gene) == 0) {
    return(" ")
  }
  else {
    return(gene)
  }
}

getProteinName <- function (x, pfdatabase = FALSE)
{
  if (pfdatabase) {
    protein = gsub("gene_product=", "", strsplit(x, "\\|")[[1]][4])
  }
  else {
    protein = strsplit(x, " OS=")[[1]][1]
  }
  if (length(protein) == 0) {
    return(" ")
  }
  else {
    return(protein)
  }
}


