
#' @title Filter ranges
#' @description \code{filter_ranges} Filter ranges
#'
#' @importFrom dplyr filter
#' @param seuratobj Seurat object
#' @param target target of interest (character)
#' @param ranges.links.final outputs of make_ranges
#' @param corrtheshold Correlation value between RNA count and ATAC peak
#' @param assay_exp Assays RNA or SCT
#' @param meantatachreshold Minimum value within the maximum value for each cluster of ATAC counts)
#' @return ranges.links.final
#' @export
#' @examples
#'


filter_ranges <- function(seuratobj,
                        target,
                        ranges.links.final,
                        corrtheshold,
                        assay_exp,
                        meantatachreshold){

  if(length(unique(ranges.links.final$peak)) == 1){
    peak.data <- seuratobj[["ATAC"]]@counts[unique(ranges.links.final$peak),]  %>%  as.data.frame()
    colnames(peak.data) <- unique(ranges.links.final$peak)
  }else{
    peak.data <- seuratobj[["ATAC"]]@counts[unique(ranges.links.final$peak),] %>% t() %>%  as.data.frame()
  }

  if(assay_exp == "RNA"){
    exp.data <- seuratobj[["RNA"]]$data[target, ] %>% as.data.frame()
    colnames(exp.data) <- target
  }else if(assay_exp == "SCT"){
    exp.data <- seuratobj[["SCT"]]$data[target, ] %>% as.data.frame()
    colnames(exp.data) <- target
  }


      cluster <- Idents(seuratobj) %>% as.data.frame()
      colnames(cluster) <- "cluster"
      peak.data$cluster <- cluster$cluster
      summarize_peak <- peak.data %>% dplyr::group_by(cluster) %>% summarize_all(mean) %>%
        dplyr::select(-cluster) %>% apply(2, max) %>% as.data.frame()
      colnames(summarize_peak) <- "max"

      summarize_peak <- summarize_peak %>% dplyr::filter(max > meantatachreshold)


      cor.data <- cbind(exp.data, peak.data[, rownames(summarize_peak)])
      colnames(cor.data) <- c(target, rownames(summarize_peak))
      res_cor <- cor(cor.data, method = "pearson") %>% as.data.frame()

      res_cor <- res_cor[which(res_cor[,1] > corrtheshold), ]


      int_peak <- intersect(rownames(res_cor), rownames(summarize_peak))
      ranges.links.final <- ranges.links.final %>% dplyr::filter(peak %in% int_peak)

      write.table(ranges.links.final, paste(TF, "_", target, "_ranges.links.final_filtered.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)

      return(ranges.links.final)
    }
