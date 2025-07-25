
# @seuratobj = EYFP_merge
# @df.genome = df.genome
# @target = target
# @tfs2 = tfs2
# @id = NULL
# @corrtheshold = 0.03 (ATAC-SCTの相関)
# @meantatachreshold = 0.2 (ATACカウントのクラスタごとのmaxの最低)
# @assay_exp = RNA or SCT
CreateCARTANet <- function(seuratobj, 
                               df.genome, 
                               target, 
                               tfs2, 
                               id){

    ####
    # targetを対象として、TSSとlinksでつながる領域を特定し、rangesとmacs2の解析すべき対象と決める
    ###
    # gene_nameで任意の遺伝子にアクセス
    # annotations[annotations$gene_name == "Col2a1"]
    # chr.min <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% min()
    # chr.max <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% max()
    # chr.tmp <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
    # chr.tmp <- paste(as.character(chr.tmp[1,1]), "-", chr.min, "-", chr.max, sep = "")
    # chr.tmp
    # # chr.tmpはtargetの遺伝子座
    # 
    # strand <-  df.genome %>% dplyr::filter(gene_name == target)
    # strand <-as.character(strand$strand)[1]
    # strand
    # # strandの向き
    # 
    # target.chr <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
    # target.chr <- target.chr$seqnames[1] %>% as.character()
    # target.chr
    # target.chrがtargetの染色体
    
    
    # TSS付近はchr.minに相当するので、chr.minとLinksでつながる領域をモチーフの検索対象とする
    # links.df <- data.frame(seuratobj[["ATAC"]]@links)
    # links.df <- links.df %>% dplyr::filter(seqnames == target.chr
  

    ranges.links <- links.df %>% dplyr::filter(gene == target)    
    # 250712
    # 下記の条件分岐の代わりに事前のLinkPeaksの結果で事前にフィルタリングを使うことにする
    # if(length(unique(ranges.links$peak)) == 1){ # debug(行列の入れ替わりが起こるため)
    #   peak.data <- seuratobj[["ATAC"]]@counts[unique(ranges.links$peak),]  %>%  as.data.frame()
    #   colnames(peak.data) <- unique(ranges.links$peak)
    # }else{
    #   peak.data <- seuratobj[["ATAC"]]@counts[unique(ranges.links$peak),] %>% t() %>%  as.data.frame()
    # }
    # 
    # if(assay_exp == "RNA"){
    #   exp.data <- seuratobj[["RNA"]]$data[target, ] %>% as.data.frame()
    #   colnames(exp.data) <- target
    # }else if(assay_exp == "SCT"){
    #   exp.data <- seuratobj[["SCT"]]$data[target, ] %>% as.data.frame()
    #   colnames(exp.data) <- target
    # } # assay_expの条件分岐
    # 
    # # クラスタごとの平均発現量を調べる
    # cluster <- Idents(seuratobj) %>% as.data.frame()
    # colnames(cluster) <- "cluster"
    # peak.data$cluster <- cluster$cluster
    # summarize_peak <- peak.data %>% dplyr::group_by(cluster) %>% summarize_all(mean) %>% 
    #   dplyr::select(-cluster) %>% apply(2, max) %>% as.data.frame()
    # colnames(summarize_peak) <- "max"
    # # ranges.links.finalのフィルタリング（ATACのカウント値のフィルタリング）
    # # クラスタのどれかで0.1より発現していればOK
    # summarize_peak <- summarize_peak %>% dplyr::filter(max > meantatachreshold)
    # 
    # # correlationする前にフィルタリングすることで計算コストを削減する
    # cor.data <- cbind(exp.data, peak.data[, rownames(summarize_peak)])
    # colnames(cor.data) <- c(target, rownames(summarize_peak))
    # res_cor <- cor(cor.data, method = "pearson") %>% as.data.frame()
    # # ranges.links.finalのフィルタリング（ATACとSCTの相関）
    # res_cor <- res_cor[which(res_cor[,1] > corrtheshold), ]
    # # res_cor <- res_cor[-1,]  # 遺伝子名除去

    # # 合わせてranges.links.finalをフィルタリング
    # int_peak <- intersect(rownames(res_cor), rownames(summarize_peak))
    # ranges.links <- ranges.links %>% dplyr::filter(peak %in% int_peak)
    # 
    # 保存
    # write.table(ranges.links.final, paste(TF, "_", target, "_ranges.links.final_filtered.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)
    
    # return(ranges.links)
    
    
    
   
    
    if(nrow(ranges.links) >= 1) {
      
      # links.dfのpeakにある値の中心がstart/endのいずれかに対応する。peakのstart/endを抽出したい
      peaks.name <- ranges.links$peak %>% str_split(pattern = "-", simplify = T) %>% data.frame()
      colnames(peaks.name) <- c("peak_seqnames", "peak_start", "peak_end")
      ranges.links <- cbind(ranges.links, peaks.name)
      
      # # chr.minが遺伝子領域を超えて5'側にまたがったピークになることがあるので、chr.minが含まれるピークを抽出。strandの向きでTSSの定義が変わるので
      # # if文で処理する
      # if(strand == "+"){
      #   # + strandの場合
      #   ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
      #   ranges.df$number <- 1:nrow(ranges.df) # あとでgranges objectから抽出するときに使う
      #   ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(start <= chr.min) %>% arrange(-start)
      #   ranges.TSS <- ranges.tmp[1,]
      #   # ranges.tmpがATAC上のTSSに相当する
      #   # links.df3はTSSの3'とリンクする領域
      #   links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end) # links.dfのstartがTSS内部にある
      #   # links.df5はTSSの5'とリンクする領域
      #   links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end) # links.dfのendがTSS内部にある
      #   
      #   links.df <- rbind(links.df3, links.df5)
      # }else{
      #   # - strandの場合
      #   ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
      #   ranges.df$number <- 1:nrow(ranges.df) # あとでgranges objectから抽出するときに使う
      #   ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(end >= chr.max) %>% arrange(start)
      #   ranges.TSS <- ranges.tmp[1,]
      #   # ranges.tmpがATAC上のTSSに相当する
      #   # links.df3はTSSの3'とリンクする領域
      #   links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end) # links.dfのstartがTSS内部にある
      #   # links.df5はTSSの5'とリンクする領域
      #   links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end) # links.dfのendがTSS内部にある
      #   
      #   links.df <- rbind(links.df3, links.df5)
      # }
      # write.table(links.df, 
      #             paste("links.df_", target, ".txt", sep = ""), quote = F, sep = "\t", row.names = F)
      # # ranges.tmpに挟まれたlinkと相関するリンクが目的のTSSと相関するゲノム領域がlinks.dfであ
      # 
      # # ranges.TSS　# TSSのピーク
      # ranges.TSS$index <- "TSS"
      
      # links.dfのピーク内にmotif_posで指定された部分があるかどうか調べる
      # links.dfはリンク情報なので遠い関係を示しているので、links.dfの属する箇所をgrangesで特定する。
      # その範囲内にmotif_posがあるものを抽出する。それをranges.showにして可視化する（motif.posだけだと10bpほどなので可視化しても見えない）
      # つまり、ATACのピークでTSSとlinksでつながるもののうち、motif.posもある領域を特定する
      
      
      motif.tmp <- motif_pos[[id]] %>% as.data.frame() %>% dplyr::filter(seqnames == as.character(ranges.links$seqnames[1]))
      
      # links.df$end、つまりTSSのリンク先がranges.dfのどの領域内にあるかを抽出するために
      # links.df$end[l]より小さいstartのうち、startが最大のものが１番links.df$end[l]に近いゲノム領域
      # これをrangesで実施する
      
      # ranges.links <- c()
      # 実は、ここもTSSの5'か3'かで場合分けが必要
      # links.df5 <- links.df %>% dplyr::filter(end == chr.min)
      # links.df3 <- links.df %>% dplyr::filter(start == chr.min)
      # 
      # for(l in 1:nrow(links.df3)){
      #   # TSSより3'のデータに対して
      #   ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(start < links.df3$end[l]) %>% arrange(-start) 
      #   ranges.tmp <- ranges.tmp[1,]
      #   ranges.links <- rbind(ranges.links, ranges.tmp)
      # } #ここでlが終わる
      # for(l in 1:nrow(links.df5)){
      #   # TSSより5'のデータに対して
      #   ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(start < links.df5$start[l]) %>% arrange(-start) 
      #   ranges.tmp <- ranges.tmp[1,]
      #   ranges.links <- rbind(ranges.links, ranges.tmp)
      # } #ここでlが終わる
      # 
      # # TSSを追加
      # ranges.links$index <- "Peaks"
      # ranges.links <- rbind(ranges.TSS, ranges.links)
      
      # motif.tmpがある領域を取得する
      # ranges.links.min <- ranges.links %>% dplyr::select(start, end) %>% min
      # ranges.links.max <- ranges.links %>% dplyr::select(start, end) %>% max
      motif.tmp.filtered <- motif.tmp %>% dplyr::filter(start > min(ranges.links %>% dplyr::select(start, end)),
                                                        end < max(ranges.links %>% dplyr::select(start, end))) #フィルタリング
      
      if(nrow(motif.tmp.filtered) > 0){
        # motif.tmp.dplyr::filteredの10 bpほどの配列データを残しておく
        using.ranges <- motif.tmp.filtered %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
        using.ranges <- StringToGRanges(using.ranges)
        range.seqs <- getSeq(genome, using.ranges) %>% data.frame()
        colnames(range.seqs) <- "sequence"
        motif.tmp.filtered <- cbind(motif.tmp.filtered, range.seqs)
        motif.filtered <- motif.tmp.filtered
        
        # 保存
        # write.table(motif.filtered, paste(TF, "_", target, "_links_motif.ranges.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)
        
        
        ranges.links.final <- c()
        for(m in 1:nrow(ranges.links)){
          for(n in 1:nrow(motif.tmp.filtered)){
            if((ranges.links$peak_start[m] < motif.tmp.filtered$start[n]) & (ranges.links$peak_end[m] > motif.tmp.filtered$end[n])){
              tmp.ranges.links.motifs <- cbind(ranges.links[m,], motif.tmp.filtered[n,])
              # ここでstart, end, width, strandの名前が2つあるのがよくないのでmotif側を変更する
              colnames(tmp.ranges.links.motifs) <- c("seqnames", "start","end", "width", "strand","score", "gene","peak", "zscore","pvalue",
                                                     "peak_seqnames", "peak_start", "peak_end",
                                                     "seqnames_motif","start_motif", "end_motif", "width_motif", 
                                                     "strand_motif", "score_motif", "sequence_motif")
              ranges.links.final <- rbind(ranges.links.final, tmp.ranges.links.motifs)
              
            }
          }#ここでnが終わる
        } #ここでmが終わる
        return(ranges.links.final)
      } # if(nrow(motif.tmp.filtered) > 0)
      
    } # if(nrow(ranges.links) >= 1)
    
    
    
    
    # if(nrow(ranges.links.final) == 0){
    #   message(paste("There are no ranges.links.final"))
    # }else{
    # #   ranges.links.final <- ranges.links.final %>% dplyr::distinct(start, end, .keep_all = T)
    # #   ranges.links.final$feature <- ranges.links.final %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
    # #   ranges.links.final
    # # 保存
    #   write.table(ranges.links.final, paste(TF, "_", target, "_ranges.links.final_pre.txt", sep = ""),  quote = F, sep = "\t", col.names =NA)
    #   
    #   return(ranges.links.final)
      
      
    } # rangs.links.finalが存在するときの処理の条件分岐

# }    




