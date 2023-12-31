% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_phoQP_barplotting.R
\name{imprints_phoQP_barplotting}
\alias{imprints_phoQP_barplotting}
\title{imprints_phoQP_barplotting}
\usage{
imprints_phoQP_barplotting(
  data,
  treatmentlevel = get_treat_level(data),
  printBothName = TRUE,
  printGeneName = FALSE,
  pfdatabase = FALSE,
  witherrorbar = TRUE,
  layout = NULL,
  colorpanel = PaletteWithoutGrey(treatmentlevel),
  ratio = 0.6,
  ret_plot = TRUE,
  save_pdf = FALSE,
  toplabel = "IMPRINTS-CETSA bar plotting",
  leftlabel = "",
  bottomlabel = "",
  pdfname = "barplot",
  pdfheight = 12,
  pdfwidth = 12
)
}
\arguments{
\item{data}{Dataframe obtained after \code{imprints_phoQP_join} or any output from 
\code{\link[IMPRINTS.CETSA::imprints_caldiff]{imprints_caldiff}}}

\item{treatmentlevel}{A vector of treatment labels, such as c("DMSO","TNFa","AT26533")
the order determines the arrangement, so in this case DMSO
group would be the first group}

\item{printBothName}{A logical to tell if you want to print the both protein names on the plot}

\item{printGeneName}{A logical to tell if you want to print the gene names on the plot}

\item{pfdatabase}{A logical for using pdf database or not}

\item{witherrorbar}{A logical to print or not the error bar on the plot}

\item{layout}{A vector indicating the panel layout for multi-panel plots per page,
default value is c(2,3) for set containing data, otherwise c(4,3), use when save_pdf = TRUE}

\item{colorpanel}{A vector of customizable color scheme provided by default with the function PaletteWithoutGrey}

\item{ratio}{Aspect ratio of the plot, default set to 0.6}

\item{ret_plot}{Logical to tell if you want to return the last plot}

\item{save_pdf}{A logical to tell if you want to save plots in a pdf file}

\item{toplabel}{Textual label at the top part of the page}

\item{leftlabel}{Textual label at the left side of the page}

\item{bottomlabel}{Textual label at the bottom part of the page}

\item{pdfname}{Textual label of the pdf file}

\item{pdfheight}{A number indicate the height of pdf file, default value 12}

\item{pdfwidth}{A number indicate the width of pdf file, default value 12}
}
\value{
The imprints barplot
}
\description{
Function to print or save barplots from peptide data and protein data like 
phosphoproteomics and quantitative proteomics and CETSA.
}
