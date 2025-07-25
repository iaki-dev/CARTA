
#' @title Clean TF motif table
#' @description \code{cleantfmotiftable} set the TF table
#'
#' @importFrom tidyverse str_to_upper
#' @importFrom tidyverse str_detect
#' @param tfs Data frame of JASPAR motif database
#' @return tfs2 which is cleaed TF table
#' @export
#' @examples
#' # tfs2 <- cleantfmotiftable(tfs = tfs)



cleantfmotiftable <- function(tfs){
  # Adjust TFS because gene names are mixed case, with parentheses, etc.
  # Convert all to uppercase
  tfs$name <- tfs$name %>% str_to_upper()

  # Clean up TF names that contain :: or (var.2), etc.
  # Normal
  normal <-  tfs$name %>% str_detect(pattern = c("VAR."), negate = T)
  normaltfs <- tfs[normal,]
  normal <- normaltfs$name %>% str_detect(pattern = c("::"), negate = T)
  normaltfs <- normaltfs[normal,]
  normaltfs
  # irregular EWSR1-FLI1
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

  # irregular
  irregular <-  tfs$name %>% str_detect(pattern = "VAR.")
  irregular2 <- tfs[irregular,]
  # ::Items marked with “::” are detected as irregular1, so exclude them in irregular2.
  irregular3 <- irregular2$name %>% str_detect(pattern = "::", negate = T)
  irregular2 <- irregular2[irregular3,]
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.2", replacement = "")
  irregular2$name <- irregular2$name %>% str_replace(pattern =　"VAR.3", replacement = "")
  irregular2$name　<- irregular2$name %>% str_remove(pattern = "\\(\\)") # Use \\ to represent special characters
  irregular2
  # irregular1
  irregular <-  tfs$name %>% str_detect(pattern = "::")
  irregular1 <- tfs[irregular,]
  tmp <- irregular1$name %>% str_split(pattern ="::", simplify = T) %>% as.data.frame()
  irregular1_1 <- cbind(tmp$V1, irregular1$ID)
  irregular1_2 <- cbind(tmp$V2, irregular1$ID)
  irregular1 <- rbind(irregular1_1, irregular1_2) %>% as.data.frame()
  # Some of the irregular1 items have var, so please correct them.
  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.2", replacement = "")
  irregular1$V1 <- irregular1$V1 %>%  str_replace(pattern =　"VAR.3", replacement = "")
  irregular1$V1　<- irregular1$V1 %>% str_remove(pattern = "\\(\\)") # Use \\ to represent special characters
  colnames(irregular1) <- c("name", "ID")

  tfs2 <- rbind(normaltfs, irregular2, irregular1, irregular4)
  #
  return(tfs2)
}

