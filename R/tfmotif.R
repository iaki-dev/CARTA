
#' @title Integration of maketfmotiftable and cleantfmotiftable
#' @description \code{tfmotif}  Make TF motif table
#'
#' @importFrom dplyr filter
#' @param PFMatrixList PFMatrixList from getMatrixSet output
#' @return tfs
#' @export
#' @examples
#' # tfs <- maketfmotiftable(PFMatrixList = PFMatrixList)




# @PFMatrixList = PFMatrixList
maketfmotiftable <- function(PFMatrixList){
  tfs <- c()
  for(i in 1:length(PFMatrixList@listData)){
    tmp <- PFMatrixList@listData[[i]]@name %>% as.data.frame()
    tmp2 <- PFMatrixList@listData[[i]]@ID %>% as.data.frame()
    tmp3 <- cbind(tmp, tmp2)
    tfs <- rbind(tfs, tmp3)
  }
  colnames(tfs) <- c("name", "ID")
  return(tfs)
}



# @tfs = tfs
cleantfmotiftable <- function(tfs){

  tfs$name <- tfs$name %>% str_to_upper()


  normal <-  tfs$name %>% str_detect(pattern = c("VAR."), negate = T) # 否定
  normaltfs <- tfs[normal,]
  normal <- normaltfs$name %>% str_detect(pattern = c("::"), negate = T) # 否定
  normaltfs <- normaltfs[normal,]

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


  irregular <-  tfs$name %>% str_detect(pattern = "VAR.")
  irregular2 <- tfs[irregular,]

  irregular3 <- irregular2$name %>% str_detect(pattern = "::", negate = T)
  irregular2 <- irregular2[irregular3,]
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.2", replacement = "")
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.3", replacement = "")
  irregular2$name　<- irregular2$name %>% str_remove(pattern = "\\(\\)")

  irregular <-  tfs$name %>% str_detect(pattern = "::")
  irregular1 <- tfs[irregular,]
  tmp <- irregular1$name %>% str_split(pattern ="::", simplify = T) %>% as.data.frame()
  irregular1_1 <- cbind(tmp$V1, irregular1$ID)
  irregular1_2 <- cbind(tmp$V2, irregular1$ID)
  irregular1 <- rbind(irregular1_1, irregular1_2) %>% as.data.frame()

  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.2", replacement = "")
  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.3", replacement = "")
  irregular1$V1　<- irregular1$V1 %>% str_remove(pattern = "\\(\\)")
  colnames(irregular1) <- c("name", "ID")

  tfs2 <- rbind(normaltfs, irregular2, irregular1, irregular4)


  return(tfs2)
}





