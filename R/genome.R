config.yaml <- file.path('extdata','config.yaml')

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
  gff_file <- ref[[version]]$gff$file

  chr <- regioneR::toGRanges(chr_bed)

  GenomeInfoDb::seqlengths(chr) <- width(chr)
  GenomeInfoDb::genome(chr) <- version
  GenomeInfoDb::isCircular(chr) <- rep(FALSE,10)
  chr
}

get_chr_GR("AGPv4")

get_gene_annot<- function( version = "AGPv4", organellar = FALSE){

  chr <- get_chr_GR("AGPv4")

  genes <- rtracklayer::import(gff_file)

  if(organellar == FALSE){
    genes <- subset(genes, seqnames %in% 1:10)
    GenomeInfoDb::seqlevels(genes) <- GenomeInfoDb::seqlevels(chr)
    GenomeInfoDb::seqinfo(genes) <- GenomeInfoDb::seqinfo(chr)
    genes <- GenomeInfoDb::sortSeqlevels(genes)
  }

  return(genes)
}


