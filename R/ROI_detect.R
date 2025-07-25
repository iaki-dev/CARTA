

#' @title ROI detect
#' @description \code{ROI_detect} Detection of Region Of Interest
#'
#' @importFrom dplyr filter
#' @param seuratobj Seurat object
#' @param target target of interest (character)
#' @param tfs tfs
#' @param fromTSS Which peak is the nth peak above or below from TSS?
#' @return txt file with TF motifs matching
#' @export
#' @examples
#'

ROI_detect <- function(seuratobj,
                        df.genome,
                        target,
                        tfs,
                        fromTSS){

  DefaultAssay(seuratobj) <- "ATAC"

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
  links.df <- links.df %>% filter(seqnames == target.chr)

  if(strand == "+"){

    ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
    ranges.df$number <- 1:nrow(ranges.df)
    ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(start <= chr.min) %>% dplyr::arrange(-start)
    ranges.TSS <- ranges.tmp[1,]


    links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end)
    links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end)

    links.df <- rbind(links.df3, links.df5)
  }else{

    ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
    ranges.df$number <- 1:nrow(ranges.df)
    ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(end >= chr.max) %>% dplyr::arrange(start)
    ranges.TSS <- ranges.tmp[1,]

    links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end)
    links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end)
    links.df <- rbind(links.df3, links.df5)


  }

  if(fromTSS == 0){
    ranges.interest <- ranges.TSS
  }else{
    if(strand == "+"){

      if(fromTSS < 0){
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start < chr.min) %>% dplyr::arrange(-start)
        ranges.interest <- ranges.interest[(-1*(fromTSS-1)),]
      } else{
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start > chr.min) %>% dplyr::arrange(start)
        ranges.interest <- ranges.interest[(fromTSS),]
      }
    } else{

      if(fromTSS < 0){
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start < chr.max) %>% dplyr::arrange(-start)
        ranges.interest <- ranges.interest[(-1*(fromTSS-1)),]
      } else{
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start > chr.max) %>% dplyr::arrange(start)
        ranges.interest <- ranges.interest[(fromTSS),]
      }
    }
  }


  ranges.interest


  gr.ranges <- ranges.interest
  gr.ranges <- seuratobj[["ATAC"]]@ranges[gr.ranges$number,]
  gr.ranges


  sub.chr.tmp <- chr.tmp %>% str_split(pattern = "-", simplify = T) %>% data.frame()
  colnames(sub.chr.tmp) <- c("seqnames", "start", "end")
  sub.chr.tmp <- rbind(sub.chr.tmp, data.frame(gr.ranges)[1:3], ranges.TSS[1:3])
  sub.chr.tmp$start <- as.numeric(sub.chr.tmp$start)
  sub.chr.tmp$end <- as.numeric(sub.chr.tmp$end)
  sub.chr.tmp <- paste(sub.chr.tmp$seqnames[1],
                       sub.chr.tmp %>% dplyr::select(start, end) %>% min(), sub.chr.tmp %>% dplyr::select(start, end) %>% max(), sep = "-")

  ranges.show <- gr.ranges
  ranges.show$color <- c("orange")


  matchtfs <- c()
  for(i in 1:nrow(tfs)){
    id <- tfs$ID[i]
    motif.tmp <- motif_pos[[id]] %>% as.data.frame() %>% dplyr::filter(seqnames == target.chr)
    if(motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end) %>% nrow() > 0){
      tmp.df <- motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end)
      tmp.df <- cbind(tfs$name[i], tmp.df)
      colnames(tmp.df)[1] <- "TF"
      using.ranges <- motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end) %>%
        dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
      using.ranges <- StringToGRanges(using.ranges)
      range.seqs <- getSeq(genome, using.ranges) %>% data.frame()
      colnames(range.seqs) <- "sequence"
      tmp.df <- cbind(tmp.df, range.seqs)


      matchtfs <- rbind(matchtfs, tmp.df)
    }
  }
  write.table(matchtfs, paste("matchtfs_", target, "_fromTSS_", fromTSS, ".txt", sep = ""), quote = F, sep = "\t", col.names = NA)

}
