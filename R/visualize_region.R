

# @seuratobj = EYFP_merge
# @df.genome = df.genome
# @target = target
# @tfs2 = tfs2
# @id = NULL
# @scorethreshold = 0.1　(ATACのTSSとの相関)
# @corrtheshold = 0.03 (ATAC-SCTの相関)
# @meantatachreshold = 0.3 (ATACカウントのクラスタごとのmaxの最低)
# @exp.assayはSCTかRNA
# goiはviolin plotで示す遺伝子
visualize_region <- function(seuratobj, 
                               df.genome, 
                               target, 
                               # tfs2,
                               # id,
                               ranges.links.final,
                               scorethreshold,
                               gr.ranges,
                               exp.assay,
                               goi){

    ####
    # targetを対象として、TSSとlinksでつながる領域を特定し、rangesとmacs2の解析すべき対象と決める
    ###
    # gene_nameで任意の遺伝子にアクセス
    # annotations[annotations$gene_name == "Col2a1"]
  chr.min <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% min()
  chr.max <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(start, end) %>% max()
  chr.tmp <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
  chr.tmp <- paste(as.character(chr.tmp[1,1]), "-", chr.min, "-", chr.max, sep = "")
  chr.tmp
  # chr.tmpはtargetの遺伝子座
  
  strand <-  df.genome %>% dplyr::filter(gene_name == target)
  strand <-as.character(strand$strand)[1]
  strand
  # strandの向き
  
  target.chr <- df.genome %>% dplyr::filter(gene_name == target) %>% dplyr::select(seqnames)
  target.chr <- target.chr$seqnames[1] %>% as.character()
  target.chr
  # target.chrがtargetの染色体
  
  
  # TSS付近はchr.minに相当するので、chr.minとLinksでつながる領域をモチーフの検索対象とする
  links.df <- data.frame(seuratobj[["ATAC"]]@links)
  # links.df <- links.df %>% dplyr::filter(seqnames == target.chr)
  ranges.links <- links.df %>% dplyr::filter(gene == target)
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
    
   
    
    # フィルタリング後にrangs.links.finalが存在するかどうかの条件分岐
    if(nrow(ranges.links.final) == 0){
      message(paste("There are no ranges.links.final after filtering"))
    }else{
      # 全体を可視化
      # 10X 
      # gr.ranges <- ranges.links.final
      # gr.ranges <- read.table("")
      # gr.ranges <- seuratobj[["ATAC"]]@ranges[gr.ranges$number,]
      # gr.ranges <- seuratobj[["ATAC"]]@ranges
      # gr.ranges
      gr.ranges <- StringToGRanges(region = ranges.links.final$peak)
      ranges.show <- gr.ranges
      ranges.show$color <- "green"
      
      
      # Genesの可視化領域を増やす()
      sub.chr.tmp <- chr.tmp %>% str_split(pattern = "-", simplify = T) %>% data.frame()
      colnames(sub.chr.tmp) <- c("seqnames", "start", "end")
      sub.chr.tmp <- rbind(sub.chr.tmp, data.frame(gr.ranges)[1:3] )#, ranges.TSS[1:3])
      sub.chr.tmp$start <- as.numeric(sub.chr.tmp$start)
      sub.chr.tmp$end <- as.numeric(sub.chr.tmp$end)
      sub.chr.tmp <- paste(sub.chr.tmp$seqnames[1], 
                           sub.chr.tmp %>% select(start, end) %>% min(), sub.chr.tmp %>% select(start, end) %>% max(), sep = "-")
      
      # 
      # 
      # # goi <- c(target,
      # #          paste(TF %>% str_sub(start = 1, end = 1), TF %>% str_sub(start = 2) %>% str_to_lower(), sep = ""))
      # 
      # # linksが多すぎて見づらいのでtargetのTSSと繋がる箇所だけにする
      # seuratobj.tmp <- seuratobj # コピー
      # target.links.df <- data.frame(seuratobj.tmp[["ATAC"]]@links)
      # # target.links.dfのうち、links.dfに該当数するものだけ抽出する
      # target.links.df <- target.links.df %>% dplyr::filter(gene == target)
      # # バグの原因として、applyのcollapseでスペースがなぜか入ってしまう現象がある。おそらく数が多いことが原因のようなので
      # # ある程度フィルタリングをかけた後にapplyするようにする
      # # links.df$feature <- links.df %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
      # # 
      # # 
      # # target.links.df3 <- target.links.df %>% dplyr::filter(start %in% links.df3$start[1])
      # # target.links.df5 <- target.links.df %>% dplyr::filter(end %in% links.df5$end[1])
      # # target.links.df <- rbind(target.links.df3, target.links.df5)
      # # target.links.df$feature <- target.links.df %>% dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
      # # 
      # # target.links.df <- target.links.df %>% dplyr::filter(feature %in% links.df$feature)
      # # target.links.df <- target.links.df %>% dplyr::filter(score > scorethreshold) # score > 0.1に変更した
      # # target.links.df %>% head
      # 
      # 
      # # dataframeからgrangesへconvert
      # links.granges <-makeGRangesFromDataFrame(target.links.df,
      #                                          keep.extra.columns=T,
      #                                          ignore.strand=FALSE,
      #                                          seqinfo=NULL,
      #                                          seqnames.field=c("seqnames", "seqname",
      #                                                           "chromosome", "chrom",
      #                                                           "chr", "chromosome_name",
      #                                                           "seqid"),
      #                                          start.field="start",
      #                                          end.field=c("end", "stop"),
      #                                          strand.field="strand",
      #                                          starts.in.df.are.0based=FALSE)
      # 
      # Links(seuratobj.tmp) <- links.granges
      # 
      # g1 <- CoveragePlot(seuratobj.tmp, region = sub.chr.tmp,　
      #              # group.by = "seurat_clusters",
      #              peaks = T, links = T,
      #              features = goi,
      #              expression.assay = exp.assay,
      #              # extend.upstream = 3000,
      #              # extend.downstream = 3000,
      #              region.highlight = ranges.show
      # )
      # 
      return(sub.chr.tmp)
    
      # dir.create(paste(goi[2], "_", goi[1],  sep = ""))
      # ggsave(paste(goi[2], "_", goi[1], ".pdf", sep = ""), width = 16, height = 8)
      
      
      
    } # ranges.links.finalの条件分岐
  } # nrow(tfs2の条件分岐)





