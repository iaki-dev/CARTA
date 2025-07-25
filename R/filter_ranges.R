


# @ranges.final = ranges.final
# @corrtheshold = 0.03 (ATAC-SCTの相関)
# @meantatachreshold = 0.2 (ATACカウントのクラスタごとのmaxの最低)
# @assay_exp = RNA or SCT

filter_ranges <- function(seuratobj, 
                        target, 
                        ranges.links.final,
                        corrtheshold,
                        assay_exp,
                        meantatachreshold){
      # ranges.links.finalのうち、RNAとの相関を調べる
      # peak.data <- seuratobj[["ATAC"]]@counts %>% as.data.frame()
  if(length(unique(ranges.links.final$peak)) == 1){ # debug(行列の入れ替わりが起こるため)
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
  } # assay_expの条件分岐
  
      
      
      # クラスタごとの平均発現量を調べる
      cluster <- Idents(seuratobj) %>% as.data.frame()
      colnames(cluster) <- "cluster"
      peak.data$cluster <- cluster$cluster
      summarize_peak <- peak.data %>% dplyr::group_by(cluster) %>% summarize_all(mean) %>% 
        dplyr::select(-cluster) %>% apply(2, max) %>% as.data.frame()
      colnames(summarize_peak) <- "max"
      # ranges.links.finalのフィルタリング（ATACのカウント値のフィルタリング）
      # クラスタのどれかで0.1より発現していればOK
      summarize_peak <- summarize_peak %>% dplyr::filter(max > meantatachreshold)

      # correlationする前にフィルタリングすることで計算コストを削減する
      cor.data <- cbind(exp.data, peak.data[, rownames(summarize_peak)])
      colnames(cor.data) <- c(target, rownames(summarize_peak))
      res_cor <- cor(cor.data, method = "pearson") %>% as.data.frame()
      # ranges.links.finalのフィルタリング（ATACとSCTの相関）
      res_cor <- res_cor[which(res_cor[,1] > corrtheshold), ]
      # res_cor <- res_cor[-1,]  # 遺伝子名除去
      
      
      
      # 合わせてranges.links.finalをフィルタリング
      int_peak <- intersect(rownames(res_cor), rownames(summarize_peak))
      ranges.links.final <- ranges.links.final %>% dplyr::filter(peak %in% int_peak)
      
      # 保存
      write.table(ranges.links.final, paste(TF, "_", target, "_ranges.links.final_filtered.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)
    
      return(ranges.links.final)
    }
