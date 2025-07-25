% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROI_detect.R
\name{ROI_detect}
\alias{ROI_detect}
\title{ROI detect}
\usage{
ROI_detect(seuratobj, df.genome, target, tfs, fromTSS)
}
\arguments{
\item{seuratobj}{Seurat object}

\item{target}{target of interest (character)}

\item{tfs}{tfs}

\item{fromTSS}{Which peak is the nth peak above or below from TSS?}
}
\value{
txt file with TF motifs matching
}
\description{
\code{ROI_detect} Detection of Region Of Interest
}
