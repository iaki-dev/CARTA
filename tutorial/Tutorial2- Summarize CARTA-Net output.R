

setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered")
NCC_merge <- readRDS("gexatac_merge_NCC_mes.rds")
DimPlot(gexatac_merge, reduction = "umap.rpca")

links.df <- data.frame(gexatac_merge[["ATAC"]]@links)
links.df <- links.df %>% dplyr::filter(score > 0)

Idents(gexatac_merge) <- gexatac_merge@meta.data$region2
markers_atac <- FindAllMarkers(gexatac_merge, assay = "ATAC")
write_tsv(markers_atac, "markers_atac.txt")


# Pharyngeal
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Pharyngeal")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Pharyngeal", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)
CARTA_Net_Pharyngeal <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                    TF %in% rownames(df_avg_rna),
                                                    peak %in% rownames(df_avg_atac))
CARTA_Net_Pharyngeal$name_ID <- paste(CARTA_Net_Pharyngeal$TF, CARTA_Net_Pharyngeal$MotifID, sep = "_")
CARTA_Net_Pharyngeal$region <- "Pharyngeal"


# Transitional
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 141712     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Transitional")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Transitional", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Transitional <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                      TF %in% rownames(df_avg_rna),
                                                      peak %in% rownames(df_avg_atac))
CARTA_Net_Transitional$name_ID <- paste(CARTA_Net_Transitional$TF, CARTA_Net_Transitional$MotifID, sep = "_")
CARTA_Net_Transitional$region <- "Transitional"


# Cardiac_Cushion
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Cardiac_Cushion")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Cardiac_Cushion", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Cardiac_Cushion <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                         TF %in% rownames(df_avg_rna),
                                                         peak %in% rownames(df_avg_atac))
CARTA_Net_Cardiac_Cushion$name_ID <- paste(CARTA_Net_Cardiac_Cushion$TF, CARTA_Net_Cardiac_Cushion$MotifID, sep = "_")
CARTA_Net_Cardiac_Cushion$region <- "Cardiac_Cushion"


# Subvalvular
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "Subvalvular")

df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "Subvalvular", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_Subvalvular <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                     TF %in% rownames(df_avg_rna),
                                                     peak %in% rownames(df_avg_atac))
CARTA_Net_Subvalvular$name_ID <- paste(CARTA_Net_Subvalvular$TF, CARTA_Net_Subvalvular$MotifID, sep = "_")
CARTA_Net_Subvalvular$region <- "Subvalvular"


# AP_septum
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "AP_septum")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "AP_septum", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_AP_septum <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                   TF %in% rownames(df_avg_rna),
                                                   peak %in% rownames(df_avg_atac))
CARTA_Net_AP_septum$name_ID <- paste(CARTA_Net_AP_septum$TF, CARTA_Net_AP_septum$MotifID, sep = "_")
CARTA_Net_AP_septum$region <- "AP_septum"


# SMC_GA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_GA")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_GA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_GA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                TF %in% rownames(df_avg_rna),
                                                peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_GA$name_ID <- paste(CARTA_Net_SMC_GA$TF, CARTA_Net_SMC_GA$MotifID, sep = "_")
CARTA_Net_SMC_GA$region <- "SMC_GA"


# SMC_DA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_DA")

df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_DA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_DA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                TF %in% rownames(df_avg_rna),
                                                peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_DA$name_ID <- paste(CARTA_Net_SMC_DA$TF, CARTA_Net_SMC_DA$MotifID, sep = "_")
CARTA_Net_SMC_DA$region <- "SMC_DA"


# SMC_CA
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net <- read_tsv("CARTA_Net_conservation_score.txt")
dim(CARTA_Net)
# [1] 950035     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(avg_conservation_score > 0.8)
dim(CARTA_Net)
# [1] 14612     26
CARTA_Net <- CARTA_Net %>% dplyr::filter(peak %in% links.df$peak)
setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/")
markers <- read_tsv("markers_region2.txt")
markers <- markers %>% dplyr::filter(avg_log2FC > 0.5,
                                     p_val < 1e-4)
markers_sub <- markers %>% dplyr::filter(cluster == "SMC_CA")
df_avg_atac <- AverageExpression(object = gexatac_merge,  group.by = "region2")$ATAC %>% as.data.frame()
df_avg_atac <- df_avg_atac %>% dplyr::filter(`Cardiac-Cushion` > 0.5)
markers_atac <- read_tsv("markers_atac.txt")
markers_atac <- markers_atac %>% dplyr::filter(cluster == "SMC_CA", avg_log2FC > 0)
df_avg_atac <- df_avg_atac[markers_atac$gene,]

df_avg_rna <- AverageExpression(object = gexatac_merge,  group.by = "region2")$RNA %>% as.data.frame()
df_avg_rna <- df_avg_rna %>% dplyr::filter(`Cardiac-Cushion` > 0.3)

CARTA_Net_SMC_CA <- CARTA_Net %>% dplyr::filter(target %in% markers_sub$gene,
                                                TF %in% rownames(df_avg_rna),
                                                peak %in% rownames(df_avg_atac))
CARTA_Net_SMC_CA$name_ID <- paste(CARTA_Net_SMC_CA$TF, CARTA_Net_SMC_CA$MotifID, sep = "_")
CARTA_Net_SMC_CA$region <- "SMC_CA"


setwd("~/workspace/data/Multiome_Wnt1/Merge_Wnt1/v14/LogNormAll/1st/filtered/peaks/CARTA_Net/v6/")
CARTA_Net_summarize <- bind_rows(CARTA_Net_Pharyngeal,CARTA_Net_Transitional,CARTA_Net_Cardiac_Cushion,
                                 CARTA_Net_Subvalvular,CARTA_Net_AP_septum,
                                 CARTA_Net_SMC_GA, CARTA_Net_SMC_DA, CARTA_Net_SMC_CA)
write_tsv(CARTA_Net_summarize, "CARTA_Net_summarize_filt_conservation.txt")


# The results of CARTA_Net indicated the combinatation between TF and target gene associated with the locus of enhancer-like reigons matcing with TF-binding motifs.
# The peak colomn indicated the peak of scATAC-seq.
# The colomns of seqnames_motif, start_motif, end_motif, and sequence_motif indicatedd the TF-binding motif site.
# The avg_conservation_score indicated the homologous region among the sepcies you analyzed.

CARTA_Net_summarize %>%
  head(3) %>%
  print(width = Inf)

# A tibble: 3 × 28
TF     target Correlation TFFORID MotifID  seqnames    start      end  width strand score gene
<chr>  <chr>        <dbl> <chr>   <chr>    <chr>       <dbl>    <dbl>  <dbl> <chr>  <dbl> <chr>
  1 Plagl1 Barx1       0.622  PLAGL1  MA1615.1 chr13    48662998 48717198  54201 *      0.152 Barx1
2 Ctcf   Bnc2        0.531  CTCF    MA1930.1 chr4     84356981 84675275 318295 *      0.148 Bnc2
3 Elf1   Bnc2        0.0816 ELF1    MA0473.3 chr4     84454756 84675275 220520 *      0.100 Bnc2
peak                    zscore      pvalue peak_seqnames peak_start peak_end seqnames_motif start_motif
<chr>                    <dbl>       <dbl> <chr>              <dbl>    <dbl> <chr>                <dbl>
  1 chr13-48716699-48717696   4.82 0.000000708 chr13           48716699 48717696 chr13             48717127
2 chr4-84356452-84357509    4.08 0.0000229   chr4            84356452 84357509 chr4              84357439
3 chr4-84454286-84455225    2.38 0.00874     chr4            84454286 84455225 chr4              84454794
end_motif width_motif strand_motif score_motif sequence_motif                      avg_conservation_score
<dbl>       <dbl> <chr>              <dbl> <chr>                                                <dbl>
  1  48717139          13 -                  10.2  CCTGACCCCAGCA                                        0.875
2  84357473          35 -                   7.52 GGAAGAAGCCTTCTGCCAAAGACTGCAAACTGCAG                  1.000
3  84454807          14 +                  11.4  AAGAAGGAAGTCGG                                       0.988
name_ID         region
<chr>           <chr>
  1 Plagl1_MA1615.1 Pharyngeal
2 Ctcf_MA1930.1   Pharyngeal
3 Elf1_MA0473.3   Pharyngeal




