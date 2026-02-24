

library(CARTA)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(tidyverse)
library(GenomicRanges)
library(future)
library(ggplot2)
library(patchwork)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)
set.seed(1234)
library(data.table)
library(progress)
library(JASPAR2022)
library(TFBSTools)
"%not.in%" <- Negate("%in%")


# Load data
gexatac_merge <- readRDS("gexatac_merge_region.rds")
df.genome <- data.frame(gexatac_merge[["ATAC"]]@annotation)

CARTA_Net_summarize <- read_tsv("CARTA_Net_conservation_score.txt")

TFi <- "Transcription fator of interst"
goi <- "target gene of interest"
target <- goi
roi <- CARTA_Net_summarize %>% dplyr::filter(target == goi, TF == TFi)


visualize_region.tmp <- visualize_region(seuratobj = gexatac_merge, df.genome = df.genome,  target = target,
                                         ranges.links.final = roi)
visualize_highlight.tmp <- visualize_highlight(seuratobj = gexatac_merge, df.genome = df.genome,  target = target,
                                               ranges.links.final = roi)
# There are too many links, so they will be only included those connected to the target TSS.
seuratobj.tmp <- gexatac_merge # temporal copy
target.links.df <- data.frame(seuratobj.tmp[["ATAC"]]@links)
# Extract only those items in target.links.df that correspond to links.df.
target.links.df <- target.links.df %>% dplyr::filter(gene == target)

# Convert from dataframe to granges
links.granges <-makeGRangesFromDataFrame(target.links.df,
                                         keep.extra.columns=T,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field=c("seqnames", "seqname",
                                                          "chromosome", "chrom",
                                                          "chr", "chromosome_name",
                                                          "seqid"),
                                         start.field="start",
                                         end.field=c("end", "stop"),
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE)

Links(seuratobj.tmp) <- links.granges
CoveragePlot(seuratobj.tmp, region = visualize_region.tmp,
             region.highlight = visualize_highlight.tmp,
             peaks = T, links = T,
             features = goi,
             expression.assay = "RNA",
)

ggsave(paste(TFi, "_", goi, ".pdf", sep = ""), width = 6, height = 4)

gc();gc()




