#' imprints_phoQP_hit_peptide
#'
#' Function to categorize peptide according to their Expression change
#'
#'
#' @param data The input data set from which categorization is performed on and hitlist is produced from;
#'    i.e. the peptide dataset after normalization
#' @param data_diff The peptide output from imprints_caldiff; can be NULL and if so will compute it; can also be a path to the data.
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
#' @param peptide_count_col The name of the column that contain the unique peptide count.
#'                          If it is sumUniPeps, don't bother with this parameter.
#' @param species The species on which you did the experiment (not necessary if already present in your datas).
#'                Default is 'Homo Sapiens'.
#'
#' @return A dataframe which contains the hits.
#'
#' @export

imprints_phoQP_hit_peptide <- function(data, data_diff = NULL, ctrl,
                              phospho = TRUE,
                              FC_cutoff = 0.5, fixed_score_cutoff = TRUE,
                              FDR = 0.01, curvature = 0.1,
                              folder_name = ifelse(phospho, "PhosphoHits_analysis", "QPpeptideHits_analysis"),
                              peptide_count_col = "peptides_counts_all",
                              species = "Homo Sapiens"){
  wd <- getwd()
  outdir <- paste0(wd, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", folder_name)
  if (dir.exists(outdir) == FALSE) {
    dir.create(outdir)
  }

  n <- nrow(data)
  name_cond <- stringr::str_subset(colnames(data), "^\\d{2}")
  if(ncol(data) - length(name_cond) != 5){
    stop("Your data must have 5 information columns,
          like 'id', 'description', 'genes', 'peptide_count' and 'protein_names' for example.")
  }
  cond <- stringr::str_split(name_cond, "_")
  cond <- unique(unlist(lapply(cond, function(x) x[3])))
  if(!(ctrl %in% cond)){
    stop("'ctrl' is not in the conditions. Please, check spelling or column names.")
  }
  cond <- cond[-which(cond == ctrl)]

  temp <- stringr::str_split(name_cond, "_")
  temp <- unique(unlist(lapply(temp, function(x) x[1])))

  if(is.null(data_diff)){
    data_diff <- data
    if(!("description" %in% colnames(data_diff))){
      to_describe <- ""
      while(length(to_describe) != 2 & all(!(to_describe %in% colnames(data)))){
        to_describe <- readline(prompt = "Type the column name of the description you want to keep in first and the gene column name in second; separated by 2 spaces"
        )
        to_describe <- stringr::str_split(to_describe, "  ")[[1]]

        if(length(to_describe) != 2){
          print("Separate the columns names by 2 spaces")
        }
        if(all(!(to_describe %in% colnames(data)))){
          print("Check columns names; it doesn't match")
        }
      }

      data_diff$description <- paste(data[[to_describe[1]]], "OS=Homo sapiens", paste0("GN=", data[[to_describe[2]]]))
      data_diff[[to_describe[1]]] <- NULL
    }
    if(!("sumUniPeps" %in% colnames(data_diff))){
      colnames(data_diff)[which(colnames(data_diff) == peptide_count_col)] <- "sumUniPeps"
    }
    message("Getting fold change...")
    colnames(data_diff)[-stringr::str_which(colnames(data_diff), "^\\d{1,}|description|sumUniPeps")] <- c("id", "sumPSMs", "countNum")
    data_diff <- data_diff[, c("id", "description", name_cond, "sumUniPeps", "sumPSMs", "countNum")]
    data_diff <- IMPRINTS.CETS::imprints_caldiff_f(data_diff, reftreatment = ctrl)
  }
  else if("character" %in% class(data_diff)){
    if(stringr::str_detect(data_diff, "\\.tsv$|\\.csv$|\\.txt$"))
      data_diff <- readr::read_tsv(data_diff)
    else
      stop("Format isn't recognize")
  }
  else if(!("data.frame" %in% class(data_diff))){
    stop("data_diff is neither a file, neither a data frame !")
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


  message("Computing mean values...")
  # get average value among bioreplicate for each protein
  diff_FC <- data_diff[,-stringr::str_which(colnames(data_diff), paste0("_", ctrl, "$"))]
  diff_FC <- tidyr::gather(diff_FC, treatment, reading, -id, -description, -sumUniPeps, -sumPSMs, -countNum)
  diff_FC <- tidyr::separate(diff_FC, treatment, into = c("temperature",
                                                          "replicate", "treatment"), sep = "_")
  diff_FC <- diff_FC %>% dplyr::group_by(id, description, temperature, treatment) %>%
    dplyr::summarise(reading.mean = mean(reading, na.rm = T))
  diff_FC <- tidyr::unite(diff_FC, treatment, temperature, treatment,
                          sep = "_")
  diff_FC <- tidyr::spread(diff_FC, treatment, reading.mean)


  message("Getting p-values...")
  for(k in cond){
    message(k)
    res <- list()
    M <- data[,stringr::str_which(colnames(data), paste0("_", ctrl, "$|_", k, "$"))]
    for(i in temp){
      message(i)
      X <- M[,stringr::str_which(colnames(M), paste0("^", i, "_"))]
      grp <- unname(sapply(colnames(X), function(x) stringr::str_split(x, "_")[[1]][3]))

      res[[i]] <- MKmisc::mod.t.test(as.matrix(X),
                                     group = factor(grp),
                                     adjust.method = "BH")$p.value
    }
    diff_FC[[paste0("pval_", k)]] <- res[[1]]
  }


  diff_FC <- diff_FC[,-c(1:2)]
  info <- data_diff[,c("id", "description", "sumUniPeps", "countNum")]
  colnames(info) <- c("id", "Position.in.protein", "Modifications", "description")
  diff_FC <- as.data.frame(cbind(info, diff_FC))

  diff_FC_plot <- diff_FC[,c("id", "Position.in.protein", "Modifications", "description", stringr::str_subset(colnames(diff_FC), "^\\d{2}|^pval_"))]
  colnames(diff_FC_plot) <- stringr::str_replace_all(colnames(diff_FC_plot), "\\d{2}", "F")
  diff_FC_plot <- tidyr::gather(diff_FC_plot, treatment, reading, -id, -Position.in.protein, -Modifications, -description)
  diff_FC_plot <- tidyr::separate(diff_FC_plot, treatment, into = c("Value", "treatment"))
  diff_FC_plot <- tidyr::spread(diff_FC_plot, Value, reading)
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


  diff_FC_plot <- diff_FC_plot %>% dplyr::group_by(id, Position.in.protein, Modifications, description, treatment) %>%
    dplyr::mutate(criteria = pval <= cutoff$BH[which(cutoff$treatment == treatment)] &
                    (FC >= cutoff$FC_pos[which(cutoff$treatment == treatment)] | FC <= cutoff$FC_neg[which(cutoff$treatment == treatment)]),
                  curve = curve(FC, cutoff$FC_neg[which(cutoff$treatment == treatment)],
                                cutoff$FC_pos[which(cutoff$treatment == treatment)],
                                cutoff$BH[which(cutoff$treatment == treatment)],
                                curvature = curvature
                  ),
                  criteria_curve = -log10(pval) >= curve
    )
  diff_FC_plot$criteria_curve <- tidyr::replace_na(diff_FC_plot$criteria_curve, FALSE)

  cond <- unique(diff_FC_plot$treatment)
  n_cond <- length(cond)
  df_curve <- data.frame(FC = rep(seq(min(diff_FC_plot$FC, na.rm = TRUE), max(diff_FC_plot$FC, na.rm = TRUE), 0.01), n_cond))
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
                                  label = paste(sub("_.*", "", id),
                                                Position.in.protein
                                                )
                                  ), show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold"))

  ggsave(paste0(format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plot.png"),
         plot = g_h,
         device = "png",
         path = outdir,
         width = 14,
         height = 8)


  g_I <- ggplot(diff_FC_plot, aes(FC, -log10(pval), color = criteria_curve)) +
    geom_point(aes(text = paste0("p-value: ", diff_FC_plot$pval, "\n",
                                 "FC: ", diff_FC_plot$FC, "\n",
                                 "Peptide sequence: ", sub("_.*", "", diff_FC_plot$id), "\n",
                                 ifelse(phospho,
                                        paste0("Phospho site: ",
                                               unlist(lapply(strsplit(diff_FC_plot$Modifications, "]; "),
                                                             function(y){
                                                               y <- y[grep("Phospho", y)];
                                                               y
                                                             }
                                                             )
                                                      ), "\n"),
                                        ""
                                        ),
                                 "Position in master protein: ", diff_FC_plot$Position.in.protein
                                 )
                   )
    ) +
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
          strip.text = element_text(face = "bold"))

  g_I <- plotly::ggplotly(g_I, tooltip = "text", width = 1080, height = 560)
  htmltools::save_html(g_I, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_plotInt.html"))

  diff_FC_plot <- diff_FC_plot[diff_FC_plot$criteria_curve,]
  diff_FC_plot$criteria <- NULL
  diff_FC_plot$curve <- NULL
  diff_FC_plot$criteria_curve <- NULL

  message("Saving datas...")
  openxlsx::write.xlsx(diff_FC, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_analysis_tab.xlsx"))
  openxlsx::write.xlsx(diff_FC_plot, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "hits_summary.xlsx"))

  if(nrow(diff_FC_plot) > 1 & length(cond) > 1){
    diff_FC_plot$treatment <- factor(diff_FC_plot$treatment)
    vennlist <- diff_FC_plot
    vennlist <- (vennlist %>% dplyr::ungroup() %>% dplyr::group_by(treatment) %>%
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
    too_long <- which(sapply(names(vennlist), stringr::str_length) > 31)
    if(length(too_long)){
      for(n in too_long){
        name_toolong <- names(vennlist)[n]
        in_common <- Reduce(intersect, stringr::str_split(stringr::str_split(name_toolong, " & ")[[1]], ""))
        name_toolong <- stringr::str_remove_all(name_toolong,
                                       paste(in_common, collapse = "|"))
        if(stringr::str_length(in_common) > 31){
          name_toolong <- stringr::str_remove_all(name_toolong, " ")
        }
        if(stringr::str_length(in_common) > 31){
          name_toolong <- paste0("&", name_toolong, "&")
          name_toolong <- stringr::str_remove_all(name_toolong, "(?<=&[a-zA-Z]).+?(?=&)")
          name_toolong <-  stringr::str_remove_all(name_toolong, "^&|&$")
        }
        names(vennlist)[n] <- name_toolong
      }
    }
    for(i in 1:length(vennlist)){
      prot <- data.frame("id" = vennlist[[i]],
                         "Position.in.protein" = info$Position.in.protein[which(!is.na(match(info$id,
                                                                                              vennlist[[i]])))],
                         "Modifications" = info$Modifications[which(!is.na(match(info$id,
                                                                             vennlist[[i]])))],
                         "description" = info$description[which(!is.na(match(info$id,
                                                                               vennlist[[i]])))]
                         )
      score_info <- diff_FC[,c(1, stringr::str_which(colnames(diff_FC), "^FC_|^pval_"))]

      vennlist[[i]] <- dplyr::left_join(prot, score_info, by = "id")
    }
    openxlsx::write.xlsx(vennlist, paste0(outdir, "/", format(Sys.time(), "%y%m%d_%H%M"), "_", "Venn tab.xlsx"))
  }

  message("Calculation done !")
  return(diff_FC_plot)
}

