config.yaml <- file.path('config.yaml')

file.exists(config.yaml)
configr::eval.config.sections(config.yaml)

config <- configr::read.config(file = config.yaml)

ref <-  configr::eval.config(
  config = "ref",
  file = config.yaml)

input <- configr::eval.config(
  config = "input",
  file = config.yaml)

output <- configr::eval.config(
  config = "output",
  file = config.yaml)



get_chr_GR<- function (version = "AGPv4"){

  chr_bed  <- ref[[version]]$bed$file


  chr <- regioneR::toGRanges(chr_bed)

  GenomeInfoDb::seqlengths(chr) <- IRanges::width(chr)
  GenomeInfoDb::genome(chr) <- version
  GenomeInfoDb::isCircular(chr) <- rep(FALSE,10)
  return(chr)
}



get_gene_annot<- function( version = "AGPv4", organellar = FALSE){
  gff_file <- ref[[version]]$gff$file
  chr <- get_chr_GR(version)
  genes <- rtracklayer::import(gff_file)

  if(organellar == FALSE){
    # genes <- subset(genes, GenomeInfoDb::seqnames(genes) %in% 1:10)
    GenomeInfoDb::seqlevels(genes, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(chr)
    GenomeInfoDb::seqinfo(genes) <- GenomeInfoDb::seqinfo(chr)
    genes <- GenomeInfoDb::sortSeqlevels(genes)
  }
  names(genes) <- genes$gene_id
  return(genes)
}


