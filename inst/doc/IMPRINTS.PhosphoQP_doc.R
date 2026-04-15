## ----eval=FALSE---------------------------------------------------------------
# if(!requireNamespace("devtools", quietly = TRUE)){ #check if you already have the devtools package
#  install.packages("devtools")  #if not, install it
# }
# devtools::install_github("nkdailingyun/IMPRINTS.CETSA")
# devtools::install_github("mgerault/IMPRINTS.CETSA.app")

## ----eval=FALSE---------------------------------------------------------------
# devtools::install_github("mgerault/IMPRINTS.PhosphoQP")

## ----message=FALSE, eval=FALSE------------------------------------------------
# library("IMPRINTS.PhosphoQP")

## ----eval=FALSE---------------------------------------------------------------
# # the treatment names corresponding to each TMT channel
# treatment <- c("B1_ctrl", "B2_ctrl", "B3_ctrl",
#                "B1_A", "B2_A", "B3_A",
#                "B1_B", "B2_B", "B3_B",
#                "B1_C", "B2_C", "B3_C",
#                "B1_C", "B2_C", "B3_C",
#                "Mix"
#                )
# # the treatment names needs to follow this type of format --> 'Rep_Condition'
# 
# data_phos <- imprints_phospho_rawread("Path/to/Phospho_data.txt", treatment)

## ----eval=FALSE---------------------------------------------------------------
# # the treatment names corresponding to each TMT channel
# treatment <- c("B1_ctrl", "B2_ctrl", "B3_ctrl",
#                "B1_A", "B2_A", "B3_A",
#                "B1_B", "B2_B", "B3_B",
#                "B1_C", "B2_C", "B3_C",
#                "B1_C", "B2_C", "B3_C",
#                "Mix"
#                )
# # again, the treatment names needs to follow this type of format --> 'Rep_Condition'
# 
# data_QPpep <- imprints_QPpep_rawread("Path/to/QPpeptide_data.txt", treatment)

## ----eval=FALSE---------------------------------------------------------------
# library(IMPRINTS.CETSA.app)
# 
# # normalization
# data_phos_s <- imprints_normalize_peptides(data_phos, dataset_name = "phospho")
# 
# # caldiff
# data_phos_s1 <- imprints_sequence_peptides(data_phos_s, control = "ctrl", dataset_name = "phospho")

## ----eval=FALSE---------------------------------------------------------------
# # normalization
# data_QPpep_s <- imprints_normalize_peptides(data_QPpep, dataset_name = "QPpep")
# 
# # caldiff
# data_QPpep_s1 <- imprints_sequence_peptides(data_QPpep_s, control = "ctrl", dataset_name = "QPpep")

## ----eval=FALSE---------------------------------------------------------------
# # hits
# data_phos_hits <- imprints_phoQP_hit_peptide(data_phos_s, data_phos_s1, "ctrl",
#                                              FC_cutoff = 0.5, FDR = 0.05)

## ----eval=FALSE---------------------------------------------------------------
# # hits
# data_QPpep_hits <- imprints_phoQP_hit_peptide(data_QPpep_s, data_QPpep_s1, "ctrl",
#                                               phospho = FALSE,
#                                               FC_cutoff = 0.5,  FDR = 0.05)

## ----eval=FALSE---------------------------------------------------------------
# # data_QP_s, the normalized QP data
# # data_QP_s1, the fold change from the normalized QP data
# 
# # it can also be an IMPRINTS.CETSA dataset containing the QP data
# # you'll need to precise the argument 'QP_name'
# 
# # hits
# data_QP_hits <- imprints_phoQP_hit_protein(data_QP_s, data_QP_s1, ctrl = "ctrl",
#                                            FC_cutoff = 0.15, FDR = 0.05,
#                                            QP_name = "36C")

## ----eval=FALSE---------------------------------------------------------------
# #  plot phospho hits
# phospho_hits_toplt <- data_phos_s1 %>%
#   right_join(data_phos_hits[,c("Positions.in.Master.Proteins", "Modifications")],
#              by = c("Positions.in.Master.Proteins", "Modifications"))
# imprints_barplotting_peptides(phospho_hits_toplt, ret_plot = FALSE, save_pdf = TRUE,
#                               pdfname = "Phospho hits")

## ----eval=FALSE---------------------------------------------------------------
# # plot QP pep hits
# QPpep_hits_toplt <- data_QPpep_s1 %>%
#   right_join(data_QPpep_hits[,c("Positions.in.Master.Proteins", "Modifications")],
#              by = c("Positions.in.Master.Proteins", "Modifications"))
# imprints_barplotting_peptides(QPpep_hits_toplt, ret_plot = FALSE, save_pdf = TRUE,
#                               pdfname = "QP peptide hits")

## ----eval=FALSE---------------------------------------------------------------
# # phospho + QPpep
# data_QPpep_phos <- imprints_phoQP_join(data_phos_s1, data_QPpep_s1)
# 
# # phospho + QP
# data_phos_QP <- imprints_phoQP_join(data_phos_s1, QP = data_QP_s1)
# 
# # QPpep + QP
# data_QPpep_QP <- imprints_phoQP_join(QPpep = data_QPpep_s1, QP = data_QP_s1)
# 
# # phos + QPpep + QP
# data_QPpep_phos_QP <- imprints_phoQP_join(data_phos_s1, data_QPpep_s1, data_QP_s1)

## ----eval=FALSE---------------------------------------------------------------
# # phos + QPpep + QP --> inner
# data_QPpep_phos_QP_inner <- imprints_phoQP_join(data_phos_s1,
#                                                 data_QPpep_s1,
#                                                 data_QP_s1,
#                                                 method = "inner")

## ----eval=FALSE---------------------------------------------------------------
# # QPpep + QP
# imprints_phoQP_barplotting(data_QPpep_QP["23283",],
#                            ret_plot = T, save_pdf = F
#                            )
# 
# # phospho + QP
# imprints_phoQP_barplotting(data_phos_QP["3000",],
#                            ret_plot = T, save_pdf = F
#                            )
# 
# # QPpep + phospho
# imprints_phoQP_barplotting(data_QPpep_phos["3000",],
#                            ret_plot = T, save_pdf = F
#                            )
# 
# # QPpep + phospho + QP
# imprints_phoQP_barplotting(data_QPpep_phos_QP["3000",],
#                            ret_plot = T, save_pdf = F
#                            )

## ----eval=FALSE---------------------------------------------------------------
# imprints_phoQP_barplotting(data_QPpep_phos_QP_inner,
#                            treatmentlevel = c("A", "B", "C", "D"),
#                            colorpanel = c("red", "green", "purple", "yellow"),
#                            ret_plot = FALSE, save_pdf = TRUE,
#                            pdfname = "common_joined"
#                            )

## ----eval=FALSE---------------------------------------------------------------
# ksea <- imprints_phoQP_ksea("path_to/analysis_tab.xlsx", TRUE)

