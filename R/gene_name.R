

#' @title Add transcript version
#' @description Adds 1 to a Gene.name string
#'
#' @param input character, gene model name
#'
#' @export
#'
#' @description Add .1 at the end of the string
# Well, there are no other versions other than 1
# I checked:
# Zm <- corncyc_pathway$Gene.name[grep("Zm",corncyc_pathway$Gene.name)]
# Zm.ver <- gsub(".*\\.","",Zm) %>% as.integer()
# table(Zm.ver)
# Zm.ver
# 1
# 9120
#' @param input Biocyc Gene Name ending with .1
#'
#' @export
add_transcript_version<-function(x){

  paste(x,1, sep = '.' )

}

#' @title Remove transcript name suffix
#' @description Trim gene name so that you get only the gene model name Zm or GRM
#'
#' @param input Biocyc Gene Name ending with .1
#'
#' @export

drop_transcript_suffix <- function(x) {
  gsub("_[TP]\\d\\d.*", "", x)
}



