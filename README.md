
# CARTA

<!-- badges: start -->
<!-- badges: end -->

("https://github.com/user-attachments/assets/fd59fa55-7a8b-4d83-b940-aed2b9faae06")

CARTA (Connected Accessible Regions by input of the combination between Transcription factor and tArget) , similar to the traditional Japanese game of karuta, where the upper phrase is read and then matched with the lower phrase, can extract candidate enhancement regions by inputting combinations of TF and target.

CARTA is R package to construct cis-regulatory networks (CARTA-Net), extracting TF-binding motifs within enhancer-like regions that positively correlate with gene expression, are co-accessible with TSSs, and are highly conserved across mammals, using only TF-target combinations of interest as input. <br>
CARTA is designed for single-cell multiome data analysis under Seurat and Signac environment.




## Installation

You can install the development version of CARTA from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("iaki-dev/CARTA")
```

## Tutorial1: Create CARTA-Net

This is a basic example which shows you how to create CARTA-Net.

``` r
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

# Load your data including scRNA-seq and scATAC-seq on each slot.
gexatac_merge <- readRDS("seuratobj.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca")
DefaultAssay(gexatac_merge) <- "ATAC"
df.genome <- data.frame(gexatac_merge[["ATAC"]]@annotation)
genome <- BSgenome.Mmusculus.UCSC.mm10


# Motif matching with peaks of scATAC-seq
search.peaks <- gexatac_merge[["ATAC"]]@ranges

opts <- list()
opts[["tax_group"]] <-  "vertebrates"
opts[["all_versions"]] <- FALSE
PFMatrixList <- getMatrixSet(JASPAR2022, opts)
PFMatrixList

motif_pos <- matchMotifs(PFMatrixList, search.peaks, genome = genome,  out = c("positions"))
saveRDS(motif_pos, "motif_pos.rds")
# motif_pos <- readRDS("motif_pos.rds")

tfs <- maketfmotiftable(PFMatrixList = PFMatrixList)
tfs2 <- cleantfmotiftable(tfs =tfs)

gexatac_merge <- RegionStats(gexatac_merge, genome = BSgenome.Mmusculus.UCSC.mm10)

# It takeks about 4 hours!
gexatac_merge <- LinkPeaks(
  object = gexatac_merge,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = rownames(gexatac_merge@assays$RNA)
)


# Markers detection
DefaultAssay(gexatac_merge) <- "RNA"
marker_rna <- FindAllMarkers(gexatac_merge, only.pos = T)
write_tsv(marker_rna, "markers_rna.txt")

DefaultAssay(gexatac_merge) <- "ATAC"
marker_atac<- FindAllMarkers(gexatac_merge, only.pos = T)
write_tsv(marker_atac, "markers_atac.txt")

DefaultAssay(gexatac_merge) <- "ATAC"



# Calculate the TF to be used for the network from here.
using_genes <- paste(substr(tfs2$name, 1, 1), tolower(substr(tfs2$name, 2, nchar(tfs2$name))), sep = "")
using_genes <- using_genes %>% unique() %>% sort()


# Limit to genes that are expressed, i.e., using.genes is limited to only transcription factors that are expressed.
df_avg <- AverageExpression(object = gexatac_merge, group.by = "seurat_clusters")$RNA %>% as.data.frame()
exp.genes <- df_avg %>%
  filter(if_any(everything(), ~ .x >= 0.1))
using_genes <- intersect(rownames(exp.genes), using_genes) %>% sort()
length(using_genes)

# Filtering marker genes
tmp.markers <- markers_rna %>% dplyr::filter(avg_log2FC > 0.2,
                                         p_val < 1e-3)

df.genome_genes <- df.genome$gene_name
markers_sub <- intersect(tmp.markers$gene, df.genome_genes)



# RNA-RNA correlation
goi <- c(using_genes, markers_sub) %>% unique()
exp.data <- df_avg[goi,] %>% t()
cor_matrix <- corSparse(exp.data)
rownames(cor_matrix) <- colnames(exp.data)
colnames(cor_matrix) <- colnames(exp.data)
df_cor <- as.data.frame(cor_matrix) %>%
  rownames_to_column("Variable")
# Convert data frame to tidy form
tidy_cor_matrix <- df_cor %>%
  gather(key="Variable2", value="Correlation", -Variable)
colnames(tidy_cor_matrix) <- c("TF", "target", "Correlation")

write_tsv(tidy_cor_matrix, "tidy_cor_matrix.txt")
nrow(tidy_cor_matrix)

tidy_cor_matrix <- tidy_cor_matrix %>% dplyr::filter(Correlation > 0,
                                                     TF %in% using_genes,
                                                     target %in% markers_sub)
nrow(tidy_cor_matrix)
write_tsv(tidy_cor_matrix, "tidy_cor_matrix_filt.txt")
# tidy_cor_matrix <- read_tsv("tidy_cor_matrix_filt.txt")

# Summarizing gene information
df.genome_group <- df.genome %>%
  dplyr::group_by(gene_name, strand, seqnames) %>%
  summarise(
    chr.min = min(start, end),
    chr.max = max(start, end),
  ) %>%
  ungroup() %>% dplyr::filter(strand %in% c("+", "-"),
                              !gene_name == "")


# Filtering by RNA-ATAC correlation of links
links.df <- data.frame(gexatac_merge[["ATAC"]]@links)
links.df <- links.df %>% filter(score > 0.1)
# Detection of specific peaks
tmp.markers_atac <- markers_atac %>% dplyr::filter(avg_log2FC > 0.1,
                                                   p_val < 0.01)
# Filter links in advance as one of the DAPs
links.df <- links.df %>% dplyr::filter(peak %in% tmp.markers_atac$gene)


tidy_cor_matrix <- tidy_cor_matrix %>% mutate(TFFORID = toupper(TF))


# Added progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed",
  total = nrow(tidy_cor_matrix), clear = FALSE, width = 60
)


# Creating CARTA-Net needs almost 6h.
CARTA_Net <- tidy_cor_matrix %>%
  mutate(MotifID = map(.x = TFFORID,
                       .f = ~tfs2 %>%
                         dplyr::filter(name == .x) %>%
                         pull(ID))) %>%
  unnest(MotifID) %>%
  mutate(result = pmap(list(seuratobj = list(gexatac_merge),
                            df.genome = list(df.genome),
                            target = target,
                            tfs2 = list(tfs2),
                            id = MotifID),
                       function(seuratobj, df.genome, target, tfs2, id) {
                         if (!pb$finished) {
                           pb$tick()
                         }
                         CreateCARTANet(seuratobj, df.genome, target, tfs2, id)
                       })) %>% unnest(result)
write_tsv(CARTA_Net, "CARTA_Net.txt")


# conservation filtering
# test_bw <-  "~/workspace/mm10.60way.phastCons60wayEuarchontoGlire.bw"
# gr <- import(test_bw)
# gr <- gr[seqnames(gr) %in% standardChromosomes(gr)]
# gr_df <- as.data.frame(gr)
# saveRDS(gr_df, "mm10.60way.phastCons60wayEuarchontoGlire.bw_df.rds")
# mm10.60way.phastCons60wayEuarchontoGlire.bw_df.rds is more than 10 GB data.
gr_df <- readRDS("mm10.60way.phastCons60wayEuarchontoGlire.bw_df.rds")
gr_df$seqnames <- as.character(gr_df$seqnames)
gr_df$start %>% class
# [1] "integer"
gr_df$end %>% class
# [1] "integer"

nrow(CARTA_Net)


library(GenomicRanges)
library(dplyr)

gr1 <- GRanges(
  seqnames = gr_df$seqnames,
  ranges = IRanges(start = gr_df$start, end = gr_df$end),
  score = gr_df$score
)

gr2 <- GRanges(
  seqnames = CARTA_Net$seqnames_motif,
  ranges = IRanges(start = CARTA_Net$start_motif, end = CARTA_Net$end_motif)
)

hits <- findOverlaps(gr2, gr1)
avg_scores <- tapply(gr1$score[subjectHits(hits)],
                     INDEX = queryHits(hits),
                     FUN = mean, na.rm = TRUE)

CARTA_Net$avg_conservation_score <- NA
CARTA_Net$avg_conservation_score[as.integer(names(avg_scores))] <- avg_scores

write.table(CARTA_Net, "CARTA_Net_conservation_score.txt", quote = F, sep = "\t", row.names = F)



```


## Tutorial2: Visualization of CARTA-Net output

This is a basic example which shows you how to visualization of CARTA-Net results.

``` r



```


