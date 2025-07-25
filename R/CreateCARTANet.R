
#' @title Create CARTA Net
#' @description \code{CreateCARTANet} Create CARTA Network
#'
#' @importFrom tidyverse str_to_upper
#' @importFrom tidyverse str_detect
#' @importFrom dplyr filter
#' @param seuratobj Seurat object
#' @param target target of interest (character)
#' @param df.genome output of data.frame(seuratobj[["ATAC"]]@annotation)
#' @param tfs2 tfs2
#' @param id Motif ID
#' @param tfs Data frame of JASPAR motif database
#'
#' @return tfs2 which is cleaed TF table
#' @export
#' @examples
#' # tfs2 <- cleantfmotiftable(tfs = tfs)
#'
#'
#'
CreateCARTANet <- function(seuratobj,
                               df.genome,
                               target,
                               tfs2,
                               id){


    ranges.links <- links.df %>% dplyr::filter(gene == target)

    if(nrow(ranges.links) >= 1) {

      peaks.name <- ranges.links$peak %>% str_split(pattern = "-", simplify = T) %>% data.frame()
      colnames(peaks.name) <- c("peak_seqnames", "peak_start", "peak_end")
      ranges.links <- cbind(ranges.links, peaks.name)

      motif.tmp <- motif_pos[[id]] %>% as.data.frame() %>% dplyr::filter(seqnames == as.character(ranges.links$seqnames[1]))


      motif.tmp.filtered <- motif.tmp %>% dplyr::filter(start > min(ranges.links %>% dplyr::select(start, end)),
                                                        end < max(ranges.links %>% dplyr::select(start, end)))
      if(nrow(motif.tmp.filtered) > 0){

        using.ranges <- motif.tmp.filtered %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
        using.ranges <- StringToGRanges(using.ranges)
        range.seqs <- getSeq(genome, using.ranges) %>% data.frame()
        colnames(range.seqs) <- "sequence"
        motif.tmp.filtered <- cbind(motif.tmp.filtered, range.seqs)
        motif.filtered <- motif.tmp.filtered



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
        return(ranges.links.final)
      }

    }



    }




