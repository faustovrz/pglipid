# The configuration is stored at "sysdata.rda"

# internal_data <- file.path("R","sysdata.rda")
# if(file.exists(internal_data)){
#   load(internal_data)
# }
#
# internal_data <- file.path("..","R","sysdata.rda")
# if(file.exists(internal_data)){
#   load(internal_data)
# }


get_chr_GR<- function (version = "AGPv4", config){

  chr_bed  <- config$ref[[version]]$bed$file


  chr <- regioneR::toGRanges(chr_bed)

  GenomeInfoDb::seqlengths(chr) <- IRanges::width(chr)
  GenomeInfoDb::genome(chr) <- version
  GenomeInfoDb::isCircular(chr) <- rep(FALSE,10)
  return(chr)
}



get_gene_annot<- function( version = "AGPv4", organellar = FALSE, config){
  gff_file <- config$ref[[version]]$gff$file
  chr <- get_chr_GR(version, config)
  genes <- rtracklayer::import(gff_file)

  if(organellar == FALSE){

    # I could not manage to make the subset function of the [
    # operator to work inside this function
    # genes <- subset(genes, GenomeInfoDb::seqnames(genes) %in% 1:10)
    # So I used the prune and the subsetting must be done in an outside script
    # that loads (attaches) GenomicRanges
    #
    GenomeInfoDb::seqlevels(genes, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(chr)
    GenomeInfoDb::seqinfo(genes) <- GenomeInfoDb::seqinfo(chr)
    genes <- GenomeInfoDb::sortSeqlevels(genes)
  }
  names(genes) <- genes$gene_id

  return(genes)
}


