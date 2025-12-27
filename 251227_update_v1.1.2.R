
# update v1.1
# add tfs, tfs2
# update v1.1.1
# add mm10.60way.phastCons60wayEuarchontoGlire.bw_df

setwd("~/Documents/GitHub/CARTA")

usethis::proj_activate(getwd())
file.exists("DESCRIPTION")

tfs <- readRDS("data/tfs.rds")
usethis::use_data(tfs, overwrite = TRUE)

tfs2 <- readRDS("data/tfs2.rds")
usethis::use_data(tfs2, overwrite = TRUE)

list.files("data")


devtools::load_all()
tfs
tfs2
?tfs
?tfs2



setwd("~/Documents/GitHub/CARTA/data")
gr_df <- readRDS("mm10.60way.phastCons60wayEuarchontoGlire.bw_df.rds")
# rdaに変換
usethis::use_data(gr_df, overwrite = TRUE)

devtools::document()
