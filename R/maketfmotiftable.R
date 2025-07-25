
#' @title Make TF motif table
#' @description \code{maketfmotiftable}  Make TF motif table
#'
#' @importFrom
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

