library(dplyr)
library(GenomicRanges)
library(regioneR)
library(liftOver)


# AGPv3 genome annotation ######################################################
#
# B73_RefGen_v3_Chr.bed
# file with each chromosome start and end coordinates
# this is useful for checking the genome version of the markers
#
# Modified from
# ftp://ftp.ensemblgenomes.org/pub/plants/release-31/fasta/zea_mays/dna_index/Zea_mays.AGPv3.dna.toplevel.fa.gz.fai

# Zea_mays.AGPv3.84_chr.gff3 gene annotation
#
# Modified from
# ftp://ftp.ensemblgenomes.org/pub/release-31/plants/gff3/zea_mays/Zea_mays.AGPv3.84.gff3.gz
#
# annot <- rtracklayer::import('/ref/zea/Zea_mays.AGPv3.84.gff3')
#
# genes <- subset(annot,
#                 type == "gene" &
#                   biotype != "transposable_element" &
#                   biotype != "pseudogene" &
#                   seqnames %in% 1:10)
#
# Add chr to the filname
#
# rtracklayer::export(genes,
#                     '/ref/zea/Zea_mays.AGPv3.84_chr.gff3',
#                     format = "gff3")

AGPv3 <- regioneR::toGRanges('/ref/zea/B73_RefGen_v3_Chr.bed')
seqlengths(AGPv3) <- width(AGPv3)
genome(AGPv3 ) <- "AGPv3"




# AGPv4 genome annotation #######################################################

# ZmB73_RefGen_v4_Chr.bed
# file with each chromosome start and end coordinates
# this is useful for checking the genome version of the markers
#
# Modified from
# ftp://ftp.ensemblgenomes.org/pub/plants/release-31/fasta/zea_mays/dna_index/Zea_mays.AGPv3.dna.toplevel.fa.gz.fai

AGPv4 <- regioneR::toGRanges("/ref/zea/ZmB73_RefGen_v4_Chr.bed")
seqlengths(AGPv4) <- width(AGPv4)
genome(AGPv4 ) <- "AGPv4"


# Zea_mays.B73_RefGen_v4.41.chr.gff3 gene annotation
#
# Downloaded from
# ftp://ftp.ensemblgenomes.org/pub/release-41/plants/gff3/zea_mays/Zea_mays.B73_RefGen_v4.41.chr.gff3.gz
#
#
annot <- rtracklayer::import('/ref/zea/Zea_mays.B73_RefGen_v4.41.chr.gff3')



genes <- subset(
  annot,
  type == "gene" &
    biotype != "transposable_element" &
    biotype != "pseudogene" &
    seqnames %in% 1:10
  ) + 10000 # + 10Kb upstream and downstream

names(genes) <- genes$gene_id


# 1. Get genes corresponding to the 5% outlier SNPs  by population

pbe_files <- c("SW_US.allPBE.txt",  "MexHigh.allPBE.txt",
               "GuaHigh.allPBE.txt", "Andes.allPBE.txt")
pops <- c("US","MH","GH","AN")

names(pbe_files) <- pops

chain_file <- "/ref/zea/chain_files/AGPv3_to_B73_RefGen_v4.chain"
ch <-import.chain(chain_file)

pbe_outlier <- data.frame( gene = genes$gene_id)
gene_hits <- data.frame()
outlier_hits <- data.frame()

for (pop in pops) {
  pbe_file <- file.path("PBE",pbe_files[pop])
  SNP <- GenomicRanges::makeGRangesFromDataFrame(
      read.delim(pbe_file, header = TRUE),
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = seqinfo(AGPv3),
      seqnames.field = "CHROM",
      start.field = "POS",
      end.field = "POS"
      ) %>%
    liftOver(ch) %>%
    unlist()

  seqinfo(AGPv3)

  seqlevels(SNP) <- as.character(1:10)
  seqinfo(SNP) <-  seqinfo(AGPv4)

  # Find genes overlapping with top 5% SNPs (outliers) in this population
  top <- 0.05
  thresh <- 1 - top

  outlier_SNP <- SNP[ SNP$PBE0 > quantile(SNP$PBE0, thresh)]
  outlier_olap <- findOverlaps(outlier_SNP, genes)
  with_outlier_SNP <- subjectHits(outlier_olap) %>% unique()

  write.csv( data.frame(
            gene = genes$gene_id[with_outlier_SNP]),
            file = paste0(pop,"_PBE_outlier_genes.csv"),
            row.names = FALSE,
            quote = FALSE)

  pbe_outlier[pop] <- 0
  pbe_outlier[with_outlier_SNP,pop] <- 1

  genic_olap <- findOverlaps(SNP,genes)
  with_PBE_SNP <- subjectHits(genic_olap) %>% unique()
  bg_genes <- genes$gene_id[with_PBE_SNP]

  write.csv( data.frame(gene = bg_genes),
             file = paste0(pop,"_PBE_bg_genes.csv"),
             row.names = FALSE,
             quote = FALSE)

  gene_hits <- rbind(
     gene_hits,
     data.frame(
       pop  = pop,
       chr  = seqnames(SNP)[queryHits(genic_olap)],
       pos  = start(SNP)[queryHits(genic_olap)],
       gene = names(genes)[subjectHits(genic_olap)]
       ) %>%
       dplyr::inner_join(
         as.data.frame(SNP),
         by=c(chr= "seqnames",pos = "start")) %>%
       arrange(chr,pos)
    )

  outlier_hits <- rbind(
    outlier_hits,
    data.frame( pop     = pop,
                gene    = names(genes)[with_outlier_SNP],
                outlier = 1)
    )

}


# Get PBE statistics for genes in these populations----------------------------
# this is really ugly but it  gets the results

pbe_outlier <- data.frame(gene = genes$gene_id) %>%
  dplyr::left_join(outlier_hits) %>%
  tidyr::pivot_wider(names_from = pop,
                     values_from = outlier,
                     values_fill = 0) %>%
  dplyr::select(-`NA`) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sum = rowSums(.[2:5])) %>%
  dplyr::left_join(
    gene_hits %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(max_PBE = max(PBE0),
                       max_PBE_pop = pop[which.max(PBE0)])) %>%
  dplyr::arrange(-sum, -max_PBE) #%>%
  # dplyr::left_join(
  #   as.data.frame(genes) %>%
  #    dplyr::select(gene =gene_id, description)
  #   )

write.csv(pbe_outlier,file ="pbe_outlier.csv", quote = FALSE, row.names = FALSE)


pbe_outlier%>%
  dplyr::left_join(outlier_hits) %>%
  dplyr::group_by(pop,gene) %>%
  dplyr::summarise(mean_PBE = mean(PBE), max_PBE = max(PBE),sum = max(sum)) %>%
  ggplot2::ggplot(aes(x=factor(sum), y=max_PBE, fill=factor(pop))) +
  ggplot2::geom_text(aes(color = factor(pop), label = shapes[factor(pop)]),
                     position = position_jitterdodge()) +
  ggplot2::xlab("# Population Selected") +
  ggplot2::ylab("Max PBE") +
  ggpubr::theme_classic2()


#-------------------------------------------------------------------------------
# Upset plot from the outliers

set_m <- pbe_outlier[,c("gene",pops,"sum")] %>%
  dplyr::filter( sum >0) %>%
  dplyr::select(gene:AN) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()


# n <- length(names(genes))
n <- length(testgeneID)
Result <- SuperExactTest::supertest(SET_input(set_m) ,n=n, degree = 2:4)

pvalues <- Result$P.value

pdf(file = "Figure2_A_SET_style.pdf")
plot(Result, Layout="landscape", sort.by="size", keep=FALSE)
dev.off()


# Blue plot

pdf(file = "Figure2_A.pdf")
plot_blue(set_m, Result$P.value)
dev.off()




