#' imprints_phoQP_hit_peptide
#'
#' Function to categorize peptide according to their Expression change
#'
#'
#' @param data The input data set from which categorization is performed on and hitlist is produced from;
#'    i.e. the peptide dataset after normalization
#' @param data_diff The peptide output from imprints_caldiff; can also be a path to the data.
#' @param ctrl The name of the control.
#' @param phospho Is your data of phosphoproteomics type. If not, should be quantiative proteomics.
#' @param FC_cutoff The fold change cutoff
#' @param fixed_score_cutoff Logical to tell if you want to use a fixed cutoff for the I-score.
#'   If TRUE, the value IS_cutoff will directly be used as the cutoff and for all treatments. If FALSE,
#'   the I-score cutoff will be calculated as the value selected for IS_cutoff plus the median of the
#'   I-scores of the proteins which have a p-value lower than the median of all p-values for a given treatment.
#'   Default is TRUE.
#' @param curvature The curvature used for the curve on the volcano plot
#' @param FDR The FDR used for the BH corrected combined p-value
#' @param folder_name The name of the folder in which you want to save the results.
#'
#' @return A dataframe which contains the hits.
#'
#' @export

imprints_phoQP_hit_peptide <- function(data, data_diff, ctrl, phospho = TRUE,
                                       FC_cutoff = 0.5, fixed_score_cutoff = TRUE,
                                       FDR = 0.01, curvature = 0.1,
                                       folder_name = ifelse(phospho, "PhosphoHits_analysis",
                                                            "QPpeptideHits_analysis")
                                       ){
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  name_cond <- grep("^\\d{2}", colnames(data), value = TRUE)
  cond <- strsplit(name_cond, "_")
  cond <- unique(sapply(cond, "[[", 3))
  if(!(ctrl %in% cond)){
    stop("'ctrl' is not in the conditions. Please, check spelling or column names.")
  }
  cond <- cond[-which(cond == ctrl)]

  if("character" %in% class(data_diff)){
    if(grepl("\\.tsv$|\\.csv$|\\.txt$", data_diff))
      data_diff <- readr::read_tsv(data_diff)
    else
      stop("Format isn't recognize")
  }
  else if(!("data.frame" %in% class(data_diff))){
    stop("data_diff is neither a file or a data frame !")
  }

  if(!("sumPSMs" %in% colnames(data)))
    data$sumPSMs <- NA
  if(!("sumPSMs" %in% colnames(data_diff)))
    data_diff$sumPSMs <- NA

  message("Computing mean values...")
  # get average value among bioreplicate for each protein
  diff_FC <- data_diff[,-grep(paste0("_", ctrl, "$"), colnames(data_diff))]
  diff_FC$Gene <- sapply(diff_FC$description, IMPRINTS.CETSA.app:::getGeneName, USE.NAMES = FALSE)

  diff_FC <- diff_FC %>%
    dplyr::rename(Sequence = Annotated.Sequence) %>%
    tidyr::gather("key", "reading",
                  -Master.Protein.Accessions, -description, -Gene,
                  -Positions.in.Master.Proteins, -Sequence,
                  -Modifications, -countNum, -sumPSMs) %>%
    tidyr::separate(key,into = c("temperature", "replicate", "treatment"), sep = "_") %>%
    dplyr::group_by(Master.Protein.Accessions, description, Gene,
                    Positions.in.Master.Proteins, Sequence, Modifications,
                    temperature, treatment) %>%
    dplyr::summarise(reading.mean = mean(reading, na.rm = T)) %>%
    tidyr::unite(treatment, temperature, treatment, sep = "_") %>%
    tidyr::spread(treatment, reading.mean)

  diff_FC$Sequence <- gsub(".*\\]\\.|\\.\\[.*", "", diff_FC$Sequence)
  diff_FC$Annotated.Sequence <- apply(diff_FC[,c("Sequence", "Positions.in.Master.Proteins", "Modifications")],
                                     1, function(x){
                                       sequence <- strsplit(x[["Sequence"]], "")[[1]]
                                       phospho <- grep("Phospho",
                                                       strsplit(x[["Modifications"]], "\\]; ")[[1]],
                                                       value = TRUE)
                                       if(length(phospho)){
                                         start_pos <- as.numeric(gsub(".* \\[|-.*", "",
                                                                      sub(";.*", "", x[["Positions.in.Master.Proteins"]]))
                                                                 )
                                         phospho <- gsub("\\[|\\]", "",
                                                         strsplit(sub(".* \\[", "", phospho), "; ")[[1]])
                                         phospho <- gsub("\\)", "%)", phospho)

                                         phospho_idx <- as.numeric(gsub("S|T|Y|\\(.*", "", phospho)) - start_pos + 1
                                         phospho <- paste0("-", phospho, "-")
                                         sequence[phospho_idx] <- phospho
                                       }
                                       sequence <- paste(sequence, collapse = "")
                                       sequence <- gsub("--", "-", sequence)
                                       gsub("^-|-$", "", sequence)
                                     })

  message("Getting p-values...")
  for(k in cond){
    message(k)
    M <- data[,grep(paste0("_", ctrl, "$|_", k, "$"), colnames(data))]
    X <- M[,grep("^\\d{2}C_", colnames(M))]
    grp <- unname(sapply(colnames(X), function(x) strsplit(x, "_")[[1]][3]))

    res <- MKmisc::mod.t.test(as.matrix(X),
                                   group = factor(grp),
                                   adjust.method = "BH")$p.value

    diff_FC[[paste0("pval_", k)]] <- res
  }

  message("Computing cutoff...")
  diff_FC_plot <- diff_FC[,c("Annotated.Sequence", "Positions.in.Master.Proteins",
                             "Modifications", "Gene",
                             grep("^\\d{2}|^pval_", colnames(diff_FC), value = TRUE))]
  colnames(diff_FC_plot) <- gsub("\\d{2}", "F", colnames(diff_FC_plot))

  diff_FC_plot <- diff_FC_plot %>%
    tidyr::gather("treatment", "reading",
                  -Annotated.Sequence, -Positions.in.Master.Proteins,
                  -Modifications, -Gene) %>%
    tidyr::separate(treatment, into = c("value", "treatment"), sep = "_") %>%
    tidyr::spread(value, reading)
  diff_FC_plot$treatment <- factor(diff_FC_plot$treatment)

  if(fixed_score_cutoff){
    cutoff <- diff_FC_plot %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(pval))/length(pval))*FDR) %>%
      dplyr::reframe(BH = find_cutoff(pval, BH),
                     FC_pos = FC_cutoff,
                     FC_neg = -FC_cutoff)
  }
  else{
    cutoff <- diff_FC_plot %>% dplyr::group_by(treatment) %>%
      dplyr::mutate(BH = (order(order(pval))/length(pval))*FDR) %>%
      dplyr::summarise(BH = find_cutoff(pval, BH),
                       FC_pos = FC_cutoff + median(FC[which(pval < quantile(pval, 0.5, na.rm = TRUE))], na.rm = TRUE),
                       FC_neg = -FC_cutoff - median(FC[which(pval < quantile(pval, 0.5, na.rm = TRUE))], na.rm = TRUE))
  }

  cutoff_file <- paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "cutoff.txt")
  readr::write_tsv(cutoff, cutoff_file)
  extra_info <- paste0("\nParameters: \nFC cutoff=", FC_cutoff, ", FDR=", FDR*100, "%, curvature=", curvature)
  write(extra_info, cutoff_file, sep = "\n", append = TRUE)

  diff_FC_plot <- diff_FC_plot %>%
    dplyr::group_by(Annotated.Sequence, Positions.in.Master.Proteins,
                    Modifications, Gene, treatment) %>%
    dplyr::mutate(criteria = pval <= cutoff$BH[which(cutoff$treatment == treatment)] &
                    (FC >= cutoff$FC_pos[which(cutoff$treatment == treatment)] | FC <= cutoff$FC_neg[which(cutoff$treatment == treatment)]),
                  curve = curve(FC, cutoff$FC_neg[which(cutoff$treatment == treatment)],
                                cutoff$FC_pos[which(cutoff$treatment == treatment)],
                                cutoff$BH[which(cutoff$treatment == treatment)],
                                curvature = curvature),
                  criteria_curve = -log10(pval) >= curve
                  )
  diff_FC_plot$criteria_curve <- tidyr::replace_na(diff_FC_plot$criteria_curve, FALSE)

  cond <- unique(diff_FC_plot$treatment)
  n_cond <- length(cond)
  df_curve <- data.frame(FC = rep(seq(min(diff_FC_plot$FC, na.rm = TRUE),
                                      max(diff_FC_plot$FC, na.rm = TRUE), 0.01), n_cond))
  df_curve$treatment <- rep(cond, each = nrow(df_curve)/n_cond)
  df_curve <- df_curve %>% dplyr::group_by(treatment, rownames(df_curve)) %>%
    dplyr::mutate(curve = curve(FC, cutoff$FC_neg[which(cutoff$treatment == treatment)],
                                cutoff$FC_pos[which(cutoff$treatment == treatment)], cutoff$BH[which(cutoff$treatment == treatment)],
                                curvature = curvature)
                  )

  message("Creating and saving plot")
  g_h <- ggplot(diff_FC_plot, aes(FC, -log10(pval), color = criteria_curve)) +
    geom_point() +
    geom_line(data = df_curve, aes(x = FC, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(diff_FC_plot$pval)))) +
    labs(title = paste("Volcano plot -",
                       ifelse(phospho, "phospho peptides", "QP peptides")),
         y = "-log10(p-value)",
         x = "Fold Change") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    facet_wrap(~treatment) +
    ggrepel::geom_label_repel(data = diff_FC_plot[diff_FC_plot$criteria_curve,],
                              aes(FC, -log10(pval),
                                  label = paste(Gene, Annotated.Sequence,
                                                sep = "\n")),
                              min.segment.length = 0, show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold", size = rel(1.6)))

  ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plot.png"),
         plot = g_h, device = "png", path = outdir,
         dpi = 300, width = 18, height = 12)


  g_I <- ggplot(diff_FC_plot, aes(FC, -log10(pval), color = criteria_curve)) +
    geom_point(aes(text = paste0("p-value: ", diff_FC_plot$pval, "\n",
                                 "FC: ", diff_FC_plot$FC, "\n",
                                 "Gene: ", diff_FC_plot$Gene, "\n",
                                 "Position in master protein: ", diff_FC_plot$Positions.in.Master.Proteins, "\n",
                                 "Peptide sequence: ",  diff_FC_plot$Annotated.Sequence
                                 )
                   )) +
    geom_line(data = df_curve, aes(x = FC, y = curve), linetype = "dashed", color = "black") +
    ylim(c(0, max(-log10(diff_FC_plot$pval)))) +
    labs(title = paste("Volcano plot -",
                       ifelse(phospho, "phospho peptides", "QP peptides")),
         y = "-log10(p-value)",
         x = "Fold Change") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"))  +
    facet_wrap(~treatment) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold", size = rel(1.6)))

  g_I <- plotly::ggplotly(g_I, tooltip = "text", width = 1080, height = 560)
  htmltools::save_html(g_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plotInt.html"))

  message("Saving datas...")
  diff_FC_plot <- diff_FC_plot[diff_FC_plot$criteria_curve,]
  diff_FC_plot$criteria <- NULL
  diff_FC_plot$curve <- NULL
  diff_FC_plot$criteria_curve <- NULL

  diff_FC <- left_join(diff_FC, data[,c("Master.Protein.Accessions", "description",
                                        "Positions.in.Master.Proteins", "Modifications",
                                        "sumPSMs", "countNum")],
                       by = c("Master.Protein.Accessions", "description",
                              "Positions.in.Master.Proteins", "Modifications"))
  diff_FC <- diff_FC[,c("Gene", "description", "Master.Protein.Accessions", "Positions.in.Master.Proteins",
                        "Sequence", "Annotated.Sequence", "Modifications", "sumPSMs", "countNum",
                        grep("^\\d{2}C_", colnames(diff_FC), value = TRUE),
                        grep("^pval_", colnames(diff_FC), value = TRUE)
                        )]

  if(nrow(diff_FC_plot))
    diff_FC <- left_join(diff_FC,
                         diff_FC_plot %>%
                           group_by(Gene, Positions.in.Master.Proteins,
                                    Annotated.Sequence, Modifications) %>%
                           summarise(hit = paste(treatment, collapse = " & ")),
                         by = c("Gene", "Positions.in.Master.Proteins",
                                "Annotated.Sequence", "Modifications"))
  else
    diff_FC$hit <- NA

  diff_FC_plot <- diff_FC_plot[,c("Gene", "Positions.in.Master.Proteins",
                                  "Annotated.Sequence", "Modifications",
                                  "treatment", "FC", "pval")]

  diff_FC <- diff_FC[order(diff_FC$Gene),]
  diff_FC_plot <- diff_FC_plot[order(diff_FC_plot$pval),]

  openxlsx::write.xlsx(diff_FC, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_analysis_tab.xlsx"))
  openxlsx::write.xlsx(diff_FC_plot, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_summary.xlsx"))

  if(nrow(diff_FC_plot) > 1 & length(cond) > 1){
    diff_FC_plot$treatment <- factor(diff_FC_plot$treatment)
    vennlist <- diff_FC_plot
    vennlist <- (vennlist %>% dplyr::ungroup() %>%
                   tidyr::unite(id, Gene, Positions.in.Master.Proteins,
                                Annotated.Sequence, Modifications,
                                sep = "_") %>%
                   dplyr::group_by(treatment) %>%
                   dplyr::summarize(vennlist = list(id)) %>%
                   dplyr::select(vennlist))[[1]]
    names(vennlist) <- levels(diff_FC_plot$treatment)
    VennDiagram::venn.diagram(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "VennDiagram.png"),
                              output = TRUE, disable.logging = TRUE,
                              imagetype = "png",
                              fill = RColorBrewer::brewer.pal(length(names(vennlist)), "Pastel2")[1:length(vennlist)],
                              fontface = "bold",
                              fontfamiliy = "sans",
                              cat.cex = 1.6,
                              cat.fontface = "bold",
                              cat.fontfamily = "sans")
    vennlist <- IMPRINTS.CETSA.app::com_protein_loop(vennlist)
    too_long <- which(sapply(names(vennlist), nchar) > 31)
    if (length(too_long)) {
      for (n in too_long) {
        name_toolong <- names(vennlist)[n]
        name_toolong <- gsub(" ", "", name_toolong)
        if (nchar(name_toolong) > 31) {
          if (grepl("&", name_toolong)) {
            in_common <- Reduce(intersect, strsplit(strsplit(name_toolong,
                                                             "&")[[1]], ""))
            if (length(in_common)) {
              name_toolong <- gsub(paste(gsub("\\.",
                                              "\\\\.", in_common), collapse = "|"),
                                   "", name_toolong)
            }
            if (nchar(name_toolong) > 31) {
              name_toolong <- paste0("&", name_toolong,
                                     "&")
              name_toolong <- stringr::str_remove_all(name_toolong,
                                                      "(?<=&[a-zA-Z]).+?(?=&)")
              name_toolong <- gsub("^&|&$", "", name_toolong)
            }
          }
          else {
            name_toolong <- strsplit(name_toolong, "")[[1]]
            to_rm <- sample(1:length(name_toolong),
                            length(name_toolong) - 31)
            name_toolong <- name_toolong[-to_rm]
            name_toolong <- paste(name_toolong, collapse = "")
          }
        }
        names(vennlist)[n] <- name_toolong
      }
    }
    for(i in 1:length(vennlist)){
      pep <- t(sapply(strsplit(vennlist[[i]], "_"), function(p) t(p)))
      pep <- as.data.frame(pep)
      colnames(pep) <- c("Gene", "Positions.in.Master.Proteins",
                         "Annotated.Sequence", "Modifications")

      vennlist[[i]] <- dplyr::left_join(pep, diff_FC, by = c("Gene", "Positions.in.Master.Proteins",
                                                                "Annotated.Sequence", "Modifications"))
    }
    openxlsx::write.xlsx(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "Venn tab.xlsx"))
  }

  message("Calculation done !")
  return(diff_FC_plot)
}

curve <- function(x, cut_neg, cut_pos, cut_p, curvature = curvature){
  y <- rep(NA, length(x))
  neg <- which(x <= cut_neg)
  pos <- which(x >= cut_pos)
  y[neg] <- curvature/abs(x[neg] - cut_neg) + -log10(cut_p)
  y[pos] <- curvature/abs(x[pos] - cut_pos) + -log10(cut_p)

  return(y)
}

find_cutoff <- function(x,y){
  id <- order(x)
  x <- x[id]
  y <- y[id]

  x <- x[which(x > y)]
  return(x[1])
}
