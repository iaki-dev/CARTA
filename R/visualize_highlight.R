
#' @title Visualize highlight of CARTA output
#' @description \code{visualize_highlight} Visualize highlight of CARTA output
#'
#' @importFrom dplyr filter
#' @param seuratobj Seurat object
#' @param df.genome output of data.frame(seuratobj[["ATAC"]]@annotation)
#' @param target target of interest (character)
#' @param ranges.links.final outputs of make_ranges and/or filtered
#' @param scorethreshold Correlation value between RNA count and ATAC peak
#' @param assay_exp Assays RNA or SCT
#' @param goi Genes to show in violin plot
#' @return ranges.show
#' @export
#' @examples
#' # visualize_highlight.tmp <- visualize_highlight(seuratobj = gexatac_merge, df.genome = df.genome,  target = target, ranges.links.final = roi)

visualize_highlight <- function(seuratobj,
                               df.genome,
                               target,
                               # tfs2,
                               # id,
                               ranges.links.final,
                               scorethreshold,
                               #gr.ranges,
                               exp.assay,
                               goi){


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


    if(nrow(ranges.links.final) == 0){
      message(paste("There are no ranges.links.final after filtering"))
    }else{

      gr.ranges <- StringToGRanges(region = ranges.links.final$peak)
      ranges.show <- gr.ranges
      ranges.show$color <- "green"


      return(ranges.show)




    }
  }





