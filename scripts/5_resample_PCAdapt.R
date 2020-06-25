source("resample.R")
library(liftOver)
library(dplyr)
library(GenomicRanges)
library(regioneR)
library(ggplot2)


# AGPv4 genome annotation -----------------------------------------------------#

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

# AGPv2 genome annotation -----------------------------------------------------#
#
# ZmB73_5b.60_FGS_chr.gff modified from
# https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v2/ZmB73_5b.60_FGS.gff3.gz
#
#  B73v2 <- rtracklayer::import('/ref/zea/ZmB73_5b.60_FGS.gff3')
#
# genes <- subset(B73v2,
#                 type == "gene" &
#                   biotype != "transposable_element" &
#                   biotype != "pseudogene" &
#                   seqnames %in% paste0("chr",1:10))
#
# not_chr <- seqlevels(genes)[!grepl("\\d",seqlevels(genes), perl =TRUE)]
# dropSeqlevels(genes, not_chr)
# renameSeqlevels(genes, gsub("chr","",seqlevels(genes)))
#
# Add chr to the file name
# rtracklayer::export(genes,
#                     '/ref/zea/ZmB73_5b.60_FGS_chr.gff3',
#                     format = "gff3")
#
# annot <- rtracklayer::import('/ref/zea/ZmB73_5b.60_FGS_chr.gff')


genes <- subset(annot,
                type == "gene" &
                  biotype != "transposable_element" &
                  biotype != "pseudogene" &
                  seqnames %in% 1:10) + 10000 # + 10Kb upstream and downstream
genes$gene_id
names(genes) <- genes$gene_id


# Loading the test gene set and SNP statistic files ############################
pbe_outlier <- read.csv("pbe_outlier.csv")
# geneset <- read.csv("glycerolipid_pathways.csv")

 #geneset <- read.csv("glycerolipid_pathways_PBE_overlap_outlier.csv")
 geneset <- read.csv("pglipid_candidates.csv")
 geneset <-  geneset %>%
   dplyr::select(gene) %>%
   dplyr::left_join(pbe_outlier)

# Loading "colm" dataframe with PCAdapt pvalues

load("pcadapt_corrected.Rimage", verbose = TRUE)
PCAdapt <- colm
colm <- NULL


# Using AGPv2 because PCAdapt data loads into a GRanges object with warnings
# about out of range SNPS if AGPv3 is used
#
# B73_RefGen_v2_Chr.bed
# file with each chromosome start and end coordinates
# this is useful for checking the genome version of the markers
# I think we must uplift all coordinates to AGPv4
#
# Indexed from
# https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v2/B73_RefGen_v2.fa.gz

AGPv2 <- toGRanges("/ref/zea/B73_RefGen_v2_Chr.bed")

seqlengths(AGPv2) <- width(AGPv2)
genome(AGPv2 ) <- "AGPv2"

chain_file <- "/ref/zea/chain_files/AGPv2_to_B73_RefGen_v4.chain"
ch <-import.chain(chain_file)
SNP <- makeGRangesFromDataFrame(PCAdapt,
                                keep.extra.columns=TRUE,
                                ignore.strand=TRUE,
                                seqinfo = seqinfo(AGPv2),
                                seqnames.field="CHR",
                                start.field="BP",
                                end.field="BP")%>%
  liftOver(ch) %>% unlist()

SNP$negLogP <- -log10(SNP$P)
nrow(PCAdapt)
nrow(mcols(SNP))
testgeneID <- as.character(geneset$gene)

genic_olap <- findOverlaps(SNP,genes)

PBE_PCAdapt <- geneset %>%
  dplyr::left_join(
    as.data.frame(genes - 10000) %>%
    dplyr::mutate(gene = gene_id, chr = seqnames) %>%
    dplyr::select(gene,chr,start,end, description)
  ) %>%
  dplyr::left_join(
    data.frame( chr  =  seqnames(SNP)[queryHits(genic_olap)],
                pos  =   start(SNP)[queryHits(genic_olap)],
                gene = names(genes)[subjectHits(genic_olap)]
    ) %>%
    dplyr::inner_join(
      as.data.frame(SNP) %>%
        dplyr::select(chr=seqnames, pos = start, negLogP)
    )
  )  %>%
  dplyr::group_by(gene) %>%
  dplyr::slice(which.max(negLogP)) %>%
  dplyr::rename(PCAdapt_peak = pos) %>%
  dplyr::arrange(-negLogP) %>% print(n = 500)

write.csv(PBE_PCAdapt,
          file ="glycerolipid_pathways_PBE_PCAdapt.csv",
          row.names = FALSE)


test_genes <- geneset$gene
top <- 0.05
thresh <- 1 - top
length(test_genes)
length(subjectHits(genic_olap) %>% unique())
nrow(mcols(genes))

PCAdapt_resample <- resample(SNP,
                          gene_annotation = genes,
                          stat = "negLogP",
                          sampling = "nr",
                          bg = "genic",
                          test =  testgeneID)

pdf(file="glycerolipid_pcadapt_negLogP_nr_genic.pdf")
plot_resample(PCAdapt_resample) +
  ggplot2::ggtitle(paste("Glycerolipid Genes",
                         hist_title(PCAdapt_resample))) +
  ggplot2::xlab( negLogP() )
dev.off()

# Non redundant sampling of  PBE selected glycerolipid genes____________________

pops <- pops <- c("US","MH","GH","AN")

selected_rs <-  list()

for (pop in pops) {

  testgeneID <- geneset$gene[geneset[,pop] == 1] %>% as.character()


  selected_rs[[pop]] <- resample(SNP,
                           gene_annotation = genes,
                           stat = "negLogP",
                           sampling = "nr",
                           bg = "genic",
                           test =  testgeneID)

}

# Plotting

selected_rs_plots<- list()

for(pop in pops ){
selected_rs_plots[[pop]] <- plot_resample(selected_rs[[pop]]) +
  ggplot2::ggtitle(
    paste(pop, selected_rs[[pop]] %>% hist_title(verbose = TRUE))) +
  ggplot2::xlab(negLogP() ) +
  ggplot2::theme(plot.title = element_text(size=12))
}


pdf(file="glycerolipid_pcadapt_negLogP_nr_genic_PBE_selected.pdf")

ggpubr::ggarrange(plotlist =selected_rs_plots,
                  # labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2 # %>%
  # ggpubr::annotate_figure(
  #   top = ggpubr::text_grob(
  #     "PCAdapt for PBE outliers",
  #     face = "bold", size = 20
  # )
  )
dev.off()

# Plot summarizing info in PBE_PCAdapt -----------------------------------------
xref <- read.delim("/ref/zea/gene_model_xref_v4.txt", header = TRUE) %>%
  dplyr::inner_join(
    geneset %>%
      dplyr::rename(v4_gene_model = gene)
  ) %>%
  dplyr::select(gene = v4_gene_model, v4_locus, v3_gene_model, Entrez) %>%
  unique()

toplot <- geneset %>% dplyr::left_join(
  data.frame( pop = pop,
              chr     =  seqnames(SNP)[queryHits(genic_olap)],
              pos      =  start(SNP)[queryHits(genic_olap)],
              gene    = names(genes)[subjectHits(genic_olap)]
    )  %>%
  dplyr::inner_join(
      as.data.frame(SNP) %>%
        dplyr::mutate( chr = seqnames, pos = start)
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
    chr = chr[1],pos =pos[which.max(negLogP)],
      max_PCAdapt = max(negLogP))) %>%
    dplyr::mutate(NumPop = factor(sum),
                  max_PBE_pop = factor(max_PBE_pop, levels = pops)) %>%
    dplyr::arrange(-max_PCAdapt) %>%
   dplyr::left_join(xref)

toplot$label <- toplot$v4_locus
toplot <- within(toplot, label[max_PBE< 1] <- "")
toplot <- within(toplot, label[max_PCAdapt < 50] <- "")

topoint <- c("GRMZM2G050641","GRMZM2G154366",
                  "GRMZM2G159890" , "GRMZM2G353444","GRMZM2G481755")
names(topoint) <- c("dgatii2","pah","gpat15", "pla1.2", "lpcat1")

toplot[match(topoint,toplot$v3_gene_model), "label"] <- names(topoint)

toplot$direction <- "Undetermined"
metabolic_direction <- c("Synthesis","Breakdown", "Synthesis","Breakdown", "Synthesis")
toplot[match(topoint,toplot$v3_gene_model), "direction"] <- metabolic_direction
toplot$direction <- factor(toplot$direction, levels = c("Undetermined", "Synthesis", "Breakdown"))
dir_pal <- c("black", "darkgreen","darkred" )
names(dir_pal) <- c("Undetermined", "Synthesis", "Breakdown")
shapes <- c(24,22,21,25,NA)

pbe_pal <- viridis::magma(5)
names(pbe_pal) <- 0:4


shp <- c(21,24,25)
names(shp) <- names(dir_pal)

p <- toplot %>%
  ggplot2::ggplot(
    aes(x=max_PCAdapt,
        y=max_PBE,
        shape = max_PBE_pop,
        label = label,
        fill  =  NumPop,
        color = direction)
   ) +
   ggplot2::geom_point(key_glyph = draw_key_rect) +
   ggplot2::geom_point(size = 5) +
   ggplot2::scale_fill_manual(values = pbe_pal) +
   ggplot2::scale_shape_manual(values = shapes ) +
   ggplot2::scale_color_manual(
     name = "Direction",
     labels = paste("<span style='color:",
                    dir_pal,
                    "'>",
                    names(dir_pal),
                    "</span>",
                    sep = ""),
     values = dir_pal) +
   ggplot2::guides(
     fill = guide_legend(
       reverse = TRUE,
       override.aes =  list( shape = "")
       ),
     shape = guide_legend(
        override.aes =  list( shape = shapes)
     ),
     color = guide_legend(
       override.aes =  list(shape = 22, alpha =0)
     )
    ) +
   ggrepel::geom_text_repel(point.padding =0.3, force = 10) +
   ggplot2::labs(fill="# PBE pops", shape="Max PBE pop")  +
   ggplot2::xlab(expression( "PCAdapt" ~~ -"Log"[10] ( italic(P) ))) +
   ggplot2::ylab("Max PBE") +
   ggpubr::theme_classic2() +
   ggplot2::theme(legend.text = ggtext::element_markdown())





pdf(file ="glycerolipid_pathways_PBE_PCAdapt.pdf")
print(p)
dev.off()


# Multivariate oulier detection

toplot %>%
  dplyr::arrange(-sum, -max_PBE) %>% tibble() %>% print(n =500)


