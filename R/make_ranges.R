
#' @title Make ranges
#' @description \code{make_ranges} Make peak ranges related TSS and gene expression
#'
#' @importFrom dplyr filter
#' @param seuratobj Seurat object
#' @param df.genome output of data.frame(seuratobj[["ATAC"]]@annotation)
#' @param target Target of interest
#' @param tfs2 output of cleantfmotiftable
#' @param id output of idcheck or interest of TF ID
#' @return ranges.links.final is peak ranges related TSS and gene expression
#' @export
#' @examples
#' # ranges.links.final <- make_ranges(seuratobj = gexatac_merge, df.genome = df.genome,  target = target, tfs2 = tfs2, id = id)


make_ranges <- function(seuratobj,
                               df.genome,
                               target,
                               tfs2,
                               id){


    chr.min <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% min()
    chr.max <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% max()
    chr.tmp <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
    chr.tmp <- paste(as.character(chr.tmp[1,1]), "-", chr.min, "-", chr.max, sep = "")
    chr.tmp

    strand <-  df.genome %>% dplyr::filter(gene_name == target)
    strand <-as.character(strand$strand)[1]
    strand

    target.chr <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
    target.chr <- target.chr$seqnames[1] %>% as.character()
    target.chr


    links.df <- data.frame(seuratobj[["ATAC"]]@links)

    ranges.links <- links.df %>% dplyr::filter(gene == target)
    peaks.name <- ranges.links$peak %>% str_split(pattern = "-", simplify = T) %>% data.frame()
    colnames(peaks.name) <- c("peak_seqnames", "peak_start", "peak_end")
    ranges.links <- cbind(ranges.links, peaks.name)


    motif.tmp <- motif_pos[[id]] %>% as.data.frame() %>% dplyr::filter(seqnames == target.chr)




    ranges.links.min <- ranges.links %>% dplyr::select(start, end) %>% min
    ranges.links.max <- ranges.links %>% dplyr::select(start, end) %>% max
    motif.tmp.filtered <- motif.tmp %>% dplyr::filter(start > ranges.links.min) %>% dplyr::filter(end < ranges.links.max) #フィルタリング

    using.ranges <- motif.tmp.filtered %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
    using.ranges <- StringToGRanges(using.ranges)
    range.seqs <- getSeq(genome, using.ranges) %>% data.frame()
    colnames(range.seqs) <- "sequence"
    motif.tmp.filtered <- cbind(motif.tmp.filtered, range.seqs)
    motif.filtered <- motif.tmp.filtered


    write.table(motif.filtered, paste(TF, "_", target, "_links_motif.ranges.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)


    ranges.links.final <- c()
    for(m in 1:nrow(ranges.links)){
      for(n in 1:nrow(motif.tmp.filtered)){
        if((ranges.links$peak_start[m] < motif.tmp.filtered$start[n]) & (ranges.links$peak_end[m] > motif.tmp.filtered$end[n])){
          tmp.ranges.links.motifs <- cbind(ranges.links[m,], motif.tmp.filtered[n,])

          colnames(tmp.ranges.links.motifs) <- c("seqnames", "start","end", "width", "strand","score", "gene","peak", "zscore","pvalue",
                                                 "peak_seqnames", "peak_start", "peak_end",
                                                 "seqnames_motif","start_motif", "end_motif", "width_motif",
                                                       "strand_motif", "score_motif", "sequence_motif")
          ranges.links.final <- rbind(ranges.links.final, tmp.ranges.links.motifs)
        }
      }
    }

    if(nrow(ranges.links.final) == 0){
      message(paste("There are no ranges.links.final"))
    }else{

      write.table(ranges.links.final, paste(TF, "_", target, "_ranges.links.final_pre.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)

      return(ranges.links.final)


    }

}




