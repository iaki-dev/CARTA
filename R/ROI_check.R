

# @seuratobj = EYFP_merge
# @df.genome = df.genome
# @target = target
# @tfs2 = tfs2
# @id = NULL

ROI_check <- function(seuratobj, 
                        df.genome, 
                        target, 
                        tfs,
                        fromTSS){
  
  ##############
  # TF-target-idの指定
  # TF <- "RUNX1"
  # target <- "Csf1r"
  
  # id <- tfs2 %>% filter(name == TF)
  # id
  # id <- "MA0774.1"
  ##############
  
  
  DefaultAssay(seuratobj) <- "ATAC"
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
  links.df <- links.df %>% filter(seqnames == target.chr)
  
  # chr.minが遺伝子領域を超えて5'側にまたがったピークになることがあるので、chr.minが含まれるピークを抽出。strandの向きでTSSの定義が変わるので
  # if文で処理する
  if(strand == "+"){
    # + strandの場合
    ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
    ranges.df$number <- 1:nrow(ranges.df) # あとでgranges objectから抽出するときに使う
    ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(start <= chr.min) %>% dplyr::arrange(-start)
    ranges.TSS <- ranges.tmp[1,]
    # ranges.tmpがATAC上のTSSに相当する
    # links.df3はTSSの3'とリンクする領域
    links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end) # links.dfのstartがTSS内部にある
    # links.df5はTSSの5'とリンクする領域
    links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end) # links.dfのendがTSS内部にある
    
    links.df <- rbind(links.df3, links.df5)
  }else{
    # - strandの場合
    ranges.df <- data.frame(seuratobj[["ATAC"]]@ranges)
    ranges.df$number <- 1:nrow(ranges.df) # あとでgranges objectから抽出するときに使う
    ranges.tmp <- ranges.df %>% dplyr::filter(seqnames == target.chr) %>% dplyr::filter(end >= chr.max) %>% dplyr::arrange(start)
    ranges.TSS <- ranges.tmp[1,]
    # ranges.tmpがATAC上のTSSに相当する
    # links.df3はTSSの3'とリンクする領域
    links.df3 <- links.df %>% dplyr::filter(start > ranges.TSS$start,  start < ranges.TSS$end) # links.dfのstartがTSS内部にある
    # links.df5はTSSの5'とリンクする領域
    links.df5 <- links.df %>% dplyr::filter(end > ranges.TSS$start,  end < ranges.TSS$end) # links.dfのendがTSS内部にある
    
    links.df <- rbind(links.df3, links.df5)
    
    
  }
  
  
  # write.table(links.df, 
  #             paste("links.df_", target, ".txt", sep = ""), quote = F, sep = "\t", row.names = F)
  # ranges.tmpに挟まれたlinkと相関するリンクが目的のTSSと相関するゲノム領域がlinks.dfである
  
  
  
  #####################
  # ranges.TSS　# TSSのピーク
  # ranges.TSS$index <- "TSS"
  # TSSから何番目のゲノムを調べるか指定する(strandの向きによって、また上流か下流かで書き方が変わる)
  # ranges.TSSのnumberはあくまで固有値で、順番とpositionが対応しないことに注意
  # fromTSSはTSSを1番目と数えて
  # TSSから-か+かで条件分岐
  # その前に+/- strandで条件分岐
  # fromTSS <- -3
  dir.create(paste(target, "_fromTSS_", fromTSS,  sep = ""))
  # change directory
  setwd(paste(target, "_fromTSS_", fromTSS,  sep = ""))
  
  if(fromTSS == 0){
    ranges.interest <- ranges.TSS
  }else{
    if(strand == "+"){
      #  + strandの場合
      if(fromTSS < 0){
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start < chr.min) %>% dplyr::arrange(-start) 
        ranges.interest <- ranges.interest[(-1*(fromTSS-1)),]
      } else{
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start > chr.min) %>% dplyr::arrange(start) 
        ranges.interest <- ranges.interest[(fromTSS),]
      }
    } else{
      # - strandの場合
      if(fromTSS < 0){
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start < chr.max) %>% dplyr::arrange(-start) 
        ranges.interest <- ranges.interest[(-1*(fromTSS-1)),]
      } else{
        ranges.interest <- ranges.df %>% dplyr::filter(seqnames == target.chr, start > chr.max) %>% dplyr::arrange(start) 
        ranges.interest <- ranges.interest[(fromTSS),]
      }
    }
  }
  
  #####################
  
  # 該当箇所を可視化
  ranges.interest
  
  # 10x
  gr.ranges <- ranges.interest
  gr.ranges <- seuratobj[["ATAC"]]@ranges[gr.ranges$number,]
  gr.ranges
  
  # Genesの可視化領域を増やす()
  sub.chr.tmp <- chr.tmp %>% str_split(pattern = "-", simplify = T) %>% data.frame()
  colnames(sub.chr.tmp) <- c("seqnames", "start", "end")
  sub.chr.tmp <- rbind(sub.chr.tmp, data.frame(gr.ranges)[1:3], ranges.TSS[1:3])
  sub.chr.tmp$start <- as.numeric(sub.chr.tmp$start)
  sub.chr.tmp$end <- as.numeric(sub.chr.tmp$end)
  sub.chr.tmp <- paste(sub.chr.tmp$seqnames[1], 
                       sub.chr.tmp %>% dplyr::select(start, end) %>% min(), sub.chr.tmp %>% dplyr::select(start, end) %>% max(), sep = "-")
  
  ranges.show <- gr.ranges
  ranges.show$color <- c("orange")
  
  plot1 <- CoveragePlot(seuratobj, region = sub.chr.tmp,　
               # group.by = "seurat_clusters",
               # peaks = T, links = T,
               features = target,
               # expression.assay = "SCT",f
               # extend.upstream = 3000,
               # extend.downstream = 3000,
               region.highlight = ranges.show
  )
  return(plot1)
  # ggsave(plot1, paste(target, "_fromTSS_", fromTSS, ".pdf", sep = ""), width = 8, height = 8)
  
  
  # # ranges.interest内にmotif_posのどれが入っているか調べる
  # # tic()
  # matchtfs <- c()
  # for(i in 1:nrow(tfs)){
  #   id <- tfs$ID[i]
  #   motif.tmp <- motif_pos[[id]] %>% as.data.frame() %>% dplyr::filter(seqnames == target.chr)
  #   if(motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end) %>% nrow() > 0){
  #     tmp.df <- motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end) 
  #     tmp.df <- cbind(tfs$name[i], tmp.df)
  #     colnames(tmp.df)[1] <- "TF"
  #     using.ranges <- motif.tmp %>% dplyr::filter(start > ranges.interest$start) %>% dplyr::filter(end < ranges.interest$end) %>% 
  #       dplyr::select(seqnames, start, end) %>% apply(1, paste, collapse = "-")
  #     using.ranges <- StringToGRanges(using.ranges)
  #     range.seqs <- getSeq(genome, using.ranges) %>% data.frame()
  #     colnames(range.seqs) <- "sequence"
  #     tmp.df <- cbind(tmp.df, range.seqs)
  #     
  #     
  #     matchtfs <- rbind(matchtfs, tmp.df)
  #   } 
  # } # ここでiが終わる
  # # toc()
  # # 14.61 sec elapsed
  # write.table(matchtfs, paste("matchtfs_", target, "_fromTSS_", fromTSS, ".txt", sep = ""), quote = F, sep = "\t", col.names = NA)
  
} 
