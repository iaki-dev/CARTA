
#' @title ID check
#' @description \code{idcheck} Check ID of TF motifs. If there are some TF motifs, Select single ID by own.
#'
#' @importFrom dplyr  filter
#' @param tfs2 output of cleantfmotiftable
#' @param TF Interest TF
#' @return output of combination TF and ID
#' @export
#' @examples
#' # id <- idcheck(tfs2 = tfs2, TF = TF)


# @tfs2 = tfs2
# @TF = TF
idcheck <- function(tfs2, TF){
  id <- tfs2 %>% dplyr::filter(name == TF)
  if(nrow(id) == 0){
    warning(paste("There is no ID."))
    id <- NULL
  }else if(nrow(id) == 1){
    id <- id$ID
    # return(id)
  }else if(nrow(id) >= 2){
    print(id)
    warning(paste("There are multiple IDs. Please choose ID from tfs2."))
    id <- NULL
  }
}

