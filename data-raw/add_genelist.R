#' add_genelist.R
#'
#' Add a pre-processed list of genes and locations to the package, to use for
#' annotation. The list is the taken directly from FUSION, a snapshop from
#' the humane genome build 19 containing only chromosome, start- and end
#' position and name: clearly, more could be done.
#'
#' Main prupose for now is to make this available as a default annotation data
#' set in line with original FUSION, but with an eye on making it possible
#' to replace it with something more recent / comprehensive.
#'
#' Alexander.Ploner@ki.se  2022-08-10

#' Run from package base dir as usually
genelist <- read.table("data-raw/glist-hg19", as.is = TRUE)
colnames(genelist) <- c("CHR","P0","P1","ID")

usethis::use_data(genelist)
