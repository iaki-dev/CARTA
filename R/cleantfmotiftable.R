



# @PFMatrixList = PFMatrixList

cleantfmotiftable <- function(tfs){
  # tfsが遺伝子名が大文字,小文字、()などばらばらなので調整する
  # すべて大文字にしてしまう
  tfs$name <- tfs$name %>% str_to_upper()
  
  # TF名に::や(var.2)などがいるので綺麗にする
  # まずnormal
  normal <-  tfs$name %>% str_detect(pattern = c("VAR."), negate = T) # 否定
  normaltfs <- tfs[normal,]
  normal <- normaltfs$name %>% str_detect(pattern = c("::"), negate = T) # 否定
  normaltfs <- normaltfs[normal,]
  normaltfs
  # irregularその4 EWSR1-FLI1
  irregular <-  tfs$name %>% str_detect(pattern = "EWSR1-FLI1")
  irregular4 <- tfs[irregular,]
  if(nrow(irregular4) > 0){
    tmp <- irregular4$name %>% str_split(pattern ="-", simplify = T) %>% as.data.frame()
    irregular4_1 <- cbind(tmp$V1, irregular4$ID)
    colnames(irregular4_1) <- c("name", "ID")
    irregular4_2 <- cbind(tmp$V2,  irregular4$ID)
    colnames(irregular4_2) <- c("name", "ID")
    irregular4 <- rbind(irregular4_1, irregular4_2) %>% as.data.frame()
    normaltfs <- normaltfs %>% dplyr::filter(name != "EWSR1-FLI1")
  }
  
  # irregularその2
  irregular <-  tfs$name %>% str_detect(pattern = "VAR.")
  irregular2 <- tfs[irregular,]
  # ::がつくものはirregular1で検出されるので、irregular2では除く
  irregular3 <- irregular2$name %>% str_detect(pattern = "::", negate = T)
  irregular2 <- irregular2[irregular3,]
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.2", replacement = "")
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.3", replacement = "")
  irregular2$name　<- irregular2$name %>% str_remove(pattern = "\\(\\)") # \\で特殊文字を表現する
  irregular2
  # irregularその1
  irregular <-  tfs$name %>% str_detect(pattern = "::")
  irregular1 <- tfs[irregular,]
  tmp <- irregular1$name %>% str_split(pattern ="::", simplify = T) %>% as.data.frame()
  irregular1_1 <- cbind(tmp$V1, irregular1$ID)
  irregular1_2 <- cbind(tmp$V2, irregular1$ID)
  irregular1 <- rbind(irregular1_1, irregular1_2) %>% as.data.frame()
  # irregular1のなかでもvarがつくものがあるので修正
  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.2", replacement = "")
  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.3", replacement = "")
  irregular1$V1　<- irregular1$V1 %>% str_remove(pattern = "\\(\\)") # \\で特殊文字を表現する
  colnames(irregular1) <- c("name", "ID")
  
  tfs2 <- rbind(normaltfs, irregular2, irregular1, irregular4)
  # これが::やvarを除去したtfsリスト。しかし重複があるので、以降で条件分岐で検証するようにする
  return(tfs2)
}
  
